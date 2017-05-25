/** YIELDING
 *
 * Runs challenging nonlinear Stokes problems that involve yielding.
 *
 ******************************************************************************
 * Author:             Johann Rudi <johann@ices.utexas.edu>
 *****************************************************************************/

#include <rhea.h>
#include <example_share_mesh.h>
#include <example_share_stokes.h>
#include <example_share_vtk.h>

#ifdef USE_CATALYST
#include <rhea_vis_adaptor.h>
#endif

/******************************************************************************
 * Main Program
 *****************************************************************************/

/**
 * Runs Stokes solver.
 */
static void
yielding_run_solver (ymir_vec_t *sol_vel_press,
                     const int iter_max, const double rel_tol,
                     rhea_stokes_problem_t *stokes_problem)
{
  const char         *this_fn_name = "yielding_run_solver";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* run solver */
  rhea_stokes_problem_solve (sol_vel_press, iter_max, rel_tol, stokes_problem);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Runs the program.
 */
int
main (int argc, char **argv)
{
  const char         *this_fn_name = "yielding:main";
  /* parallel environment */
  MPI_Comm            mpicomm = MPI_COMM_WORLD;
  int                 mpisize, mpirank, ompsize;
  /* options */
  ymir_options_t               *opt;
  rhea_domain_options_t         domain_options;
  rhea_temperature_options_t    temp_options;
  rhea_viscosity_options_t      visc_options;
  rhea_discretization_options_t discr_options;
  rhea_newton_options_t         newton_options;
  /* options local to this program */
  int                 solver_iter_max;
  double              solver_rel_tol;
  char               *vtk_write_input_path;
  char               *vtk_write_newton_itn_path;
  char               *vtk_write_solution_path;
#ifdef USE_CATALYST
  char               *vis_catalyst_script;
#endif
  /* mesh */
  p4est_t            *p4est;
  ymir_mesh_t        *ymir_mesh;
  ymir_pressure_elem_t  *press_elem;
  /* Stokes */
  rhea_stokes_problem_t *stokes_problem;
  ymir_vec_t         *sol_vel_press;

  /*
   * Initialize Program
   */

  /* begin program initialization */
  rhea_init_begin (&mpisize, &mpirank, &ompsize, argc, argv, mpicomm);

  /* create options */
  opt = ymir_options_global_new (argv[0] /* program path */);
  rhea_add_options_base (opt);

  /* add options of this program */
  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  /* solver options */
  YMIR_OPTIONS_I, "solver-iter-max", '\0',
    &solver_iter_max, 100,
    "Maximum number of iterations for Stokes solver",
  YMIR_OPTIONS_D, "solver-rel-tol", '\0',
    &solver_rel_tol, 1.0e-6,
    "Relative tolerance for Stokes solver",

  /* vtk output options */
  YMIR_OPTIONS_S, "vtk-write-input-path", '\0',
    &(vtk_write_input_path), NULL,
    "File path for vtk files for the input of the Stokes problem",
  YMIR_OPTIONS_S, "vtk-write-newton-itn-path", '\0',
    &(vtk_write_newton_itn_path), NULL,
    "File path for vtk files for iterations of Newton's method",
  YMIR_OPTIONS_S, "vtk-write-solution-path", '\0',
    &(vtk_write_solution_path), NULL,
    "File path for vtk files for the solution of the Stokes problem",

  /* visualization */
#ifdef USE_CATALYST
  YMIR_OPTIONS_S, "vis-catalyst-script", '\0',
    &(vis_catalyst_script), NULL,
    "Python script for ParaView-Catalyst",
#endif

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add sub-options */
  rhea_add_options_all (opt);
  ymir_options_add_suboptions_solver_stokes (opt);

  /* end program initialization */
  rhea_init_end (opt);

  /* initialize ParaView-Catalyst */
#ifdef USE_CATALYST
  CatalystInitialize (1, &vis_catalyst_script);
#endif

  /*
   * Print Environment and Options
   */

  RHEA_GLOBAL_PRODUCTIONF (
      "Into %s (production %i)\n", this_fn_name, rhea_get_production_run ());
  RHEA_GLOBAL_PRODUCTIONF (
      "Parallel environment: MPI size %i, OpenMP size %i\n", mpisize, ompsize);

  /* print & process options */
  ymir_options_print_summary (SC_LP_INFO, opt);
  rhea_process_options_all (&domain_options, &temp_options,
                            &visc_options, &discr_options,
                            &newton_options);

  /*
   * Setup Mesh
   */

  example_share_mesh_new (&p4est, &ymir_mesh, &press_elem, mpicomm,
                          &domain_options, &discr_options);

  /*
   * Setup Stokes Problem
   */

  example_share_stokes_new (&stokes_problem, ymir_mesh, press_elem,
                            &domain_options, &temp_options, &visc_options,
                            &newton_options, vtk_write_newton_itn_path);

  /* write vtk of input data */
  example_share_vtk_write_input_data (vtk_write_input_path, stokes_problem,
                                      &temp_options, &visc_options);

  /*
   * Solve Stokes Problem
   */

  /* initialize solution vector */
  sol_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);

  /* setup solver */
  rhea_stokes_problem_setup_solver (stokes_problem);

  /* run solver */
  yielding_run_solver (sol_vel_press, solver_iter_max, solver_rel_tol,
                       stokes_problem);

  /* write vtk of solution */
  example_share_vtk_write_solution (vtk_write_solution_path, sol_vel_press,
                                    stokes_problem);
#ifdef USE_CATALYST
  {
    const ymir_locidx_t n_elements = ymir_mesh_get_num_elems_loc (ymir_mesh);
    const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (
                                                                    ymir_mesh);
    const int           order = ymir_mesh_get_order (ymir_mesh);
    ymir_locidx_t       elid;
    int                 nodeid;

    unsigned int        n_coordinates = n_elements * n_nodes_per_el;
    double             *coordinates;
    unsigned int       *element_data;
    double             *velocity_data;
    double             *pressure_data;

    RHEA_ASSERT (1 <= order && order <= 2);

    coordinates = RHEA_ALLOC (double, n_coordinates * 3);
    element_data = RHEA_ALLOC (unsigned int, n_coordinates);
    velocity_data = RHEA_ALLOC (double, n_coordinates * 3);
    pressure_data = RHEA_ALLOC (double, n_elements);

    for (elid = 0; elid < n_elements; elid++) {
      const double *x = ymir_mesh_get_elem_coord_x (elid, ymir_mesh);
      const double *y = ymir_mesh_get_elem_coord_y (elid, ymir_mesh);
      const double *z = ymir_mesh_get_elem_coord_z (elid, ymir_mesh);

      for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
        const ymir_locidx_t idx = n_nodes_per_el*elid + nodeid;

        /* set coordinates data */
        coordinates[3*idx    ] = x[nodeid];
        coordinates[3*idx + 1] = y[nodeid];
        coordinates[3*idx + 2] = z[nodeid];
        element_data[idx] = idx;

        /* set velocity and pressure */
        velocity_data[3*idx    ] = sin (M_PI * x[nodeid]);
        velocity_data[3*idx + 1] = 0.0;
        velocity_data[3*idx + 2] = 0.0;
        pressure_data[elid] = 1.0;
      }
    }

    CatalystCoProcess (n_coordinates, coordinates,
                       (unsigned int) n_elements, element_data,
                       order,
                       velocity_data, pressure_data);

    RHEA_FREE (coordinates);
    RHEA_FREE (element_data);
    RHEA_FREE (velocity_data);
    RHEA_FREE (pressure_data);
  }
#endif

  /* destroy */
  rhea_velocity_pressure_destroy (sol_vel_press);

  /*
   * Finalize
   */

  /* destroy Stokes problem */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);
  example_share_stokes_destroy (stokes_problem);

  /* destroy mesh */
  example_share_mesh_destroy (ymir_mesh, press_elem, p4est, &discr_options);

  /* destroy options */
  ymir_options_global_destroy ();

  /* finalize ParaView-Catalyst */
#ifdef USE_CATALYST
  CatalystFinalize ();
#endif

  /* print that this function is ending */
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

  /* finalize rhea */
  rhea_finalize ();

  return 0;
}
