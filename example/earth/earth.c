/** EARTH
 *
 * Runs global mantle convection.
 *
 ******************************************************************************
 * Author:             Johann Rudi <johann@ices.utexas.edu>
 *****************************************************************************/

#include <rhea.h>
#include <rhea_domain_subset.h>
#include <rhea_stokes_problem_amr.h>
#include <example_share_mesh.h>
#include <example_share_stokes.h>
#include <example_share_io.h>
#include <example_share_vtk.h>
#include <ymir_vtk.h>

/******************************************************************************
 * Monitoring
 *****************************************************************************/

/* perfomance monitor tags and names */
typedef enum
{
  RHEA_MAIN_PERFMON_SETUP_MESH,
  RHEA_MAIN_PERFMON_SETUP_STOKES,
  RHEA_MAIN_PERFMON_SETUP_SOLVER,
  RHEA_MAIN_PERFMON_SOLVE,
  RHEA_MAIN_PERFMON_TOTAL,
  RHEA_MAIN_PERFMON_N
}
rhea_main_performance_monitor_idx_t;

static const char  *rhea_main_performance_monitor_name[RHEA_MAIN_PERFMON_N] =
{
  "Setup Mesh",
  "Setup Stokes",
  "Setup Solver",
  "Solve",
  "Total"
};

/******************************************************************************
 * Post-Processing
 *****************************************************************************/

/** Writes the viscous stress tensor as VTK and TXT output. */
static void
earth_write_viscstress_tensor (ymir_vec_t *sol_vel_press,
                               rhea_stokes_problem_t *stokes_problem,
                               const char *vtk_path,
                               const char *txt_aleutian_path)
{
  rhea_domain_options_t      *domain_options;
  rhea_temperature_options_t *temp_options;
  rhea_viscosity_options_t   *visc_options;
  ymir_mesh_t          *ymir_mesh;
  ymir_pressure_elem_t *press_elem;
  ymir_vec_t         *velocity, *viscosity;
  ymir_vec_t         *strainrate, *viscstress;

  /* exit if nothing to do */
  if (vtk_path == NULL && txt_aleutian_path == NULL) {
    return;
  }

  /* get options */
  domain_options = rhea_stokes_problem_get_domain_options (stokes_problem);
  temp_options = rhea_stokes_problem_get_temperature_options (stokes_problem);
  visc_options = rhea_stokes_problem_get_viscosity_options (stokes_problem);

  /* get mesh */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);

  /* get fields */
  velocity = rhea_velocity_new (ymir_mesh);
  viscosity = rhea_viscosity_new (ymir_mesh);
  rhea_velocity_pressure_copy_components (velocity, NULL, sol_vel_press,
                                          press_elem);
  rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

  /* compute strain rate tensor */
  strainrate = rhea_strainrate_new (ymir_mesh);
  rhea_strainrate_compute (strainrate, velocity);
  rhea_velocity_destroy (velocity);

  /* compute viscous stress tensor */
  viscstress = rhea_stress_new (ymir_mesh);
  rhea_stress_compute_viscstress (viscstress, strainrate, viscosity);
  rhea_viscosity_destroy (viscosity);
  rhea_strainrate_destroy (strainrate);

  /* convert to physical dimensions */
  rhea_stress_convert_to_dimensional_Pa (viscstress, domain_options,
                                         temp_options, visc_options);


  /* write vtk file */
  if (vtk_path != NULL) {
    ymir_vec_t         *viscstress_diag, *viscstress_offdiag;

    viscstress_diag = ymir_dvec_new (ymir_mesh, 3, viscstress->node_type);
    viscstress_offdiag = ymir_dvec_new (ymir_mesh, 3, viscstress->node_type);
    rhea_stress_separate_diag_offdiag (viscstress_diag, viscstress_offdiag,
                                       viscstress);

    ymir_vtk_write (ymir_mesh, vtk_path,
                    viscstress_diag, "viscstress_diag",
                    viscstress_offdiag, "viscstress_offdiag", NULL);
    ymir_vec_destroy (viscstress_diag);
    ymir_vec_destroy (viscstress_offdiag);
  }

  /* write txt file */
  if (txt_aleutian_path != NULL) {
    ymir_gloidx_t       n_points_in_subset_global;

    rhea_domain_subset_write_txt_aleutian (txt_aleutian_path, viscstress,
                                           &n_points_in_subset_global,
                                           domain_options);
  }

  /* destroy */
  rhea_stress_destroy (viscstress);
}

/******************************************************************************
 * Main Program
 *****************************************************************************/

/**
 * Runs the program.
 */
int
main (int argc, char **argv)
{
  static const char   func_name[] = "earth:main";
  /* parallel environment */
  MPI_Comm            mpicomm = sc_MPI_COMM_WORLD;
  int                 mpisize, mpirank, ompsize;
  /* options */
  ymir_options_t               *opt;
  rhea_domain_options_t         domain_options;
  rhea_temperature_options_t    temp_options;
  rhea_plate_options_t          plate_options;
  rhea_weakzone_options_t       weak_options;
  rhea_topography_options_t     topo_options;
  rhea_viscosity_options_t      visc_options;
  rhea_discretization_options_t discr_options;
  /* options local to this program */
  int                 solver_iter_max;
  double              solver_rel_tol;
  char               *velocity_file_path_bin;
  char               *pressure_file_path_bin;
  int                *write_solution_amr;
  char               *bin_solver_path;
  char               *vtk_input_path;
  char               *vtk_solution_path;
  char               *vtk_solver_path;
  char               *txt_solution_surf_path;
  char               *txt_aleutian_path;
  /* mesh */
  p4est_t              *p4est;
  ymir_mesh_t          *ymir_mesh;
  ymir_pressure_elem_t *press_elem;
  /* Stokes */
  rhea_stokes_problem_t  *stokes_problem;
  ymir_vec_t         *sol_vel_press;
  int                 nonzero_inital_guess;

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
    &(solver_iter_max), 100,
    "Maximum number of iterations for Stokes solver",
  YMIR_OPTIONS_D, "solver-rel-tol", '\0',
    &(solver_rel_tol), 1.0e-6,
    "Relative tolerance for Stokes solver",

  /* velocity & pressure (initial guess) intput */
  YMIR_OPTIONS_S, "velocity-file-path", '\0',
    &(velocity_file_path_bin), NULL,
    "Path to a binary file that contains a velocity field",
  YMIR_OPTIONS_S, "pressure-file-path", '\0',
    &(pressure_file_path_bin), NULL,
    "Path to a binary file that contains a pressure field",

  /* file output */
  YMIR_OPTIONS_B, "write-solution-amr", '\0',
    &(write_solution_amr), 0,
    "Perform AMR before writing output of the solution",

  /* binary file output */
  YMIR_OPTIONS_S, "bin-write-solver-path", '\0',
    &(bin_solver_path), NULL,
    "Bin file path for solver internals (e.g., iterations of Newton's method)",

  /* vtk file output */
  YMIR_OPTIONS_S, "vtk-write-input-path", '\0',
    &(vtk_input_path), NULL,
    "VTK file path for the input of the Stokes problem",
  YMIR_OPTIONS_S, "vtk-write-solution-path", '\0',
    &(vtk_solution_path), NULL,
    "VTK file path for the solution of the Stokes problem",
  YMIR_OPTIONS_S, "vtk-write-solver-path", '\0',
    &(vtk_solver_path), NULL,
    "VTK file path for solver internals (e.g., iterations of Newton's method)",

  /* text file output */
  YMIR_OPTIONS_S, "txt-write-solution-surf-path", '\0',
    &(txt_solution_surf_path), NULL,
    "TXT file path for velocity and surface stress at surface",
  YMIR_OPTIONS_S, "txt-write-aleutian-path", '\0',
    &(txt_aleutian_path), NULL,
    "TXT file path for the viscous stress tensor of the Aleutian subduction",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add sub-options */
  rhea_add_options_all (opt);
  ymir_options_add_suboptions_solver_stokes (opt);

  /* end program initialization */
  rhea_init_end (opt);

  /* initialize performance monitors */
  rhea_performance_monitor_init (rhea_main_performance_monitor_name,
                                 RHEA_MAIN_PERFMON_N);

  /* start performance monitors */
  rhea_performance_monitor_start_barrier (RHEA_MAIN_PERFMON_TOTAL);

  /*
   * Print Environment and Options
   */

  RHEA_GLOBAL_PRODUCTIONF (
      "Into %s (production %i)\n", func_name, rhea_production_run_get ());
  RHEA_GLOBAL_PRODUCTIONF (
      "Parallel environment: MPI size %i, OpenMP size %i\n", mpisize, ompsize);

  /* print & process options */
  ymir_options_print_summary (SC_LP_INFO, opt);
  rhea_process_options_all (&domain_options, &temp_options, &plate_options,
                            &weak_options, &topo_options, &visc_options,
                            &discr_options);

  /*
   * Setup Mesh
   */

  example_share_mesh_new (&p4est, &ymir_mesh, &press_elem, mpicomm,
                          &domain_options, &topo_options, &discr_options,
                          RHEA_MAIN_PERFMON_SETUP_MESH);

  /*
   * Setup Stokes Problem
   */

  example_share_stokes_new (&stokes_problem, &ymir_mesh, &press_elem,
                            &temp_options, &plate_options,
                            &weak_options, &visc_options,
                            p4est, &discr_options,
                            RHEA_MAIN_PERFMON_SETUP_MESH,
                            RHEA_MAIN_PERFMON_SETUP_STOKES,
                            bin_solver_path, vtk_solver_path);

  /* write vtk of input data */
  example_share_vtk_write_input_data (vtk_input_path, stokes_problem,
                                      &plate_options);

  /*
   * Solve Stokes Problem
   */

  /* setup solver */
  rhea_performance_monitor_start_barrier (RHEA_MAIN_PERFMON_SETUP_SOLVER);
  rhea_stokes_problem_setup_solver (stokes_problem);
  rhea_performance_monitor_stop_add (RHEA_MAIN_PERFMON_SETUP_SOLVER);

  /* initialize solution vector */
  sol_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);
  if (velocity_file_path_bin != NULL && pressure_file_path_bin != NULL) {
    int                 read_success;

    read_success = rhea_velocity_pressure_read (
        sol_vel_press, velocity_file_path_bin, pressure_file_path_bin,
        press_elem, mpicomm);
    if (!read_success) {
      RHEA_LERRORF (
        "%s: Failed reading velocity & pressure from files: \"%s\", \"%s\"\n",
        func_name, velocity_file_path_bin, pressure_file_path_bin);
    }
    nonzero_inital_guess = (read_success != 0);
  }
  else {
    nonzero_inital_guess = 0;
  }

  /* run solver */
  rhea_performance_monitor_start_barrier (RHEA_MAIN_PERFMON_SOLVE);
  rhea_stokes_problem_solve (&sol_vel_press, nonzero_inital_guess,
                             solver_iter_max, solver_rel_tol, stokes_problem);
  rhea_performance_monitor_stop_add (RHEA_MAIN_PERFMON_SOLVE);

  /* write vtk of solution */
  example_share_vtk_write_solution (vtk_solution_path, sol_vel_press,
                                    stokes_problem);

  /* write various output, possibly on a coarse mesh */
  {
    int                 amr_iter = 0;

    /* run AMR before writing solution */
    if (write_solution_amr) {
      rhea_stokes_problem_set_velocity_pressure (stokes_problem,
                                                 sol_vel_press);
      amr_iter = rhea_stokes_problem_solution_amr (stokes_problem, p4est,
                                                   &discr_options);
      sol_vel_press = rhea_stokes_problem_get_velocity_pressure (
          stokes_problem);
    }

    /* write vtk of coarse solution */
    if (vtk_solution_path != NULL && 0 < amr_iter) {
      char                path[BUFSIZ];

      snprintf (path, BUFSIZ, "%s_coarse", vtk_solution_path);
      example_share_vtk_write_solution (path, sol_vel_press, stokes_problem);
    }

    /* write vtk & txt of (coarse) viscous stress tensor */
    if (vtk_solution_path != NULL || txt_aleutian_path != NULL) {
      char                vtk_path[BUFSIZ];
      char                aleu_path[BUFSIZ];

      snprintf (vtk_path, BUFSIZ, "%s_viscstress", vtk_solution_path);
      snprintf (aleu_path, BUFSIZ, "%s_viscstress", txt_aleutian_path);
      earth_write_viscstress_tensor (
          sol_vel_press, stokes_problem,
          (vtk_solution_path != NULL ? vtk_path : NULL),
          (txt_aleutian_path != NULL ? aleu_path : NULL));
    }

    /* write txt of (coarse) solution at surface */
    example_share_io_write_solution_surf_txt (
        txt_solution_surf_path, RHEA_DOMAIN_COORDINATE_SPHERICAL_GEO,
        sol_vel_press, stokes_problem);
  }

  /*
   * Clear Stokes Problem & Mesh
   */

  /* destroy solution */
  rhea_velocity_pressure_destroy (sol_vel_press);

  /* destroy Stokes problem */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);
  example_share_stokes_destroy (stokes_problem, &temp_options, &plate_options,
                                &weak_options, &visc_options);

  /* destroy mesh */
  example_share_mesh_destroy (ymir_mesh, press_elem, p4est, &topo_options,
                              &discr_options);

  /*
   * Finalize
   */

  /* stop performance monitors */
  rhea_performance_monitor_stop_add_barrier (RHEA_MAIN_PERFMON_TOTAL);

  /* print performance statistics */
  rhea_performance_monitor_print (func_name,
                                  1 /* print wtime */,
                                  0 /* print #calls */,
                                  0 /* print flops */,
                                  1 /* print ymir */);
  rhea_performance_monitor_finalize ();

  /* destroy options */
  ymir_options_global_destroy ();

  /* print that this function is ending */
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", func_name);

  /* finalize rhea */
  rhea_finalize ();

  return 0;
}
