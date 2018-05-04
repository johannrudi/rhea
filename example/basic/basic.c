/** BASIC
 *
 * Runs rhea's elemental solver components for incompressible Stokes systems
 * with linear and nonlinear rheologies.
 *
 ******************************************************************************
 * Author:             Johann Rudi <johann@ices.utexas.edu>
 *****************************************************************************/

#include <rhea.h>
#include <example_share_mesh.h>
#include <example_share_stokes.h>
#include <example_share_vtk.h>

/******************************************************************************
 * Monitoring
 *****************************************************************************/

/* perfomance monitor tags and names */
typedef enum
{
  RHEA_MAIN_PERFMON_SETUP_MESH,
  RHEA_MAIN_PERFMON_SETUP_STOKES,
  RHEA_MAIN_PERFMON_TOTAL,
  RHEA_MAIN_PERFMON_N
}
rhea_main_performance_monitor_idx_t;

static const char  *rhea_main_performance_monitor_name[RHEA_MAIN_PERFMON_N] =
{
  "Setup Mesh",
  "Setup Stokes",
  "Total"
};

/******************************************************************************
 * Main Program
 *****************************************************************************/

/**
 * Runs the program.
 */
int
main (int argc, char **argv)
{
  static const char   func_name[] = "basic:main";
  /* parallel environment */
  MPI_Comm            mpicomm = sc_MPI_COMM_WORLD;
  int                 mpisize, mpirank, ompsize;
  /* options */
  ymir_options_t               *opt;
  rhea_domain_options_t         domain_options;
  rhea_temperature_options_t    temp_options;
  rhea_weakzone_options_t       weak_options;
  rhea_viscosity_options_t      visc_options;
  rhea_topography_options_t     topo_options;
  rhea_discretization_options_t discr_options;
  rhea_newton_options_t         newton_options;
  /* options local to this program */
  int                 solver_iter_max;
  double              solver_rel_tol;
  char               *bin_solver_path;
  char               *vtk_input_path;
  char               *vtk_solution_path;
  char               *vtk_solver_path;
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
  rhea_process_options_all (&domain_options, &temp_options, &weak_options,
                            &visc_options, &topo_options, &discr_options,
                            &newton_options);

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
                            &temp_options, &weak_options, &visc_options,
                            &newton_options, p4est, &discr_options,
                            RHEA_MAIN_PERFMON_SETUP_MESH,
                            RHEA_MAIN_PERFMON_SETUP_STOKES,
                            bin_solver_path, vtk_solver_path);

  /* write vtk of input data */
  example_share_vtk_write_input_data (vtk_input_path, stokes_problem,
                                      &temp_options, &visc_options);

  /*
   * Solve Stokes Problem
   */

  /* setup solver */
  rhea_stokes_problem_setup_solver (stokes_problem);

  /* initialize solution vector */
  sol_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);

  /* run solver */
  rhea_stokes_problem_solve (&sol_vel_press, solver_iter_max, solver_rel_tol,
                             stokes_problem);

  /* write vtk of solution */
  example_share_vtk_write_solution (vtk_solution_path, sol_vel_press,
                                    stokes_problem);

  /*
   * Clear Stokes Problem & Mesh
   */

  /* destroy solution */
  rhea_velocity_pressure_destroy (sol_vel_press);

  /* destroy Stokes problem */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);
  example_share_stokes_destroy (stokes_problem, &temp_options, &weak_options,
                                &visc_options);

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
                                  0 /* print ymir */);
  rhea_performance_monitor_finalize ();

  /* destroy options */
  ymir_options_global_destroy ();

  /* print that this function is ending */
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", func_name);

  /* finalize rhea */
  rhea_finalize ();

  return 0;
}
