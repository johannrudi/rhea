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
  rhea_weakzone_options_t       weak_options;
  rhea_viscosity_options_t      visc_options;
  rhea_discretization_options_t discr_options;
  rhea_newton_options_t         newton_options;
  /* options local to this program */
  int                 solver_iter_max;
  double              solver_rel_tol;
  char               *vtk_write_input_path;
  char               *vtk_write_newton_itn_path;
  char               *vtk_write_solution_path;
#ifdef RHEA_USE_CATALYST
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
#ifdef RHEA_USE_CATALYST
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

  /*
   * Print Environment and Options
   */

  RHEA_GLOBAL_PRODUCTIONF (
      "Into %s (production %i)\n", this_fn_name, rhea_get_production_run ());
  RHEA_GLOBAL_PRODUCTIONF (
      "Parallel environment: MPI size %i, OpenMP size %i\n", mpisize, ompsize);

  /* print & process options */
  ymir_options_print_summary (SC_LP_INFO, opt);
  rhea_process_options_all (&domain_options, &temp_options, &weak_options,
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
                            &temp_options, &weak_options, &visc_options,
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

#ifdef RHEA_USE_CATALYST
  if (vis_catalyst_script != NULL) {
    /* initialize visualization */
    rhea_vis_initialize ((const char **) &vis_catalyst_script, 1);

    /* process visualization */
    {
      ymir_vec_t         *velocity, *pressure;

      velocity = rhea_velocity_new (ymir_mesh);
      pressure = rhea_pressure_new (ymir_mesh, press_elem);
      rhea_velocity_pressure_copy_components (velocity, pressure, sol_vel_press,
                                              press_elem);
      rhea_vis_process_primary (velocity, pressure);
      rhea_velocity_destroy (velocity);
      rhea_pressure_destroy (pressure);
    }

    /* finalize visualization */
    rhea_vis_finalize ();
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
  example_share_stokes_destroy (stokes_problem, &temp_options, &weak_options,
                                &visc_options);

  /* destroy mesh */
  example_share_mesh_destroy (ymir_mesh, press_elem, p4est, &discr_options);

  /* destroy options */
  ymir_options_global_destroy ();

  /* print that this function is ending */
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

  /* finalize rhea */
  rhea_finalize ();

  return 0;
}
