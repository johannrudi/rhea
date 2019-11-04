/** PACIFIC_XSECTION_INVERSION
 *
 * Runs an inverse problem solver for a mantle convection model that is
 * supported on a cross section of Earth across the Pacific.
 */

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
  RHEA_MAIN_PERFMON_SOLVE_INVERSION,
  RHEA_MAIN_PERFMON_TOTAL,
  RHEA_MAIN_PERFMON_N
}
rhea_main_performance_monitor_idx_t;

static const char  *rhea_main_performance_monitor_name[RHEA_MAIN_PERFMON_N] =
{
  "Setup Mesh",
  "Setup Stokes",
  "Solve Inversion",
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
  static const char   func_name[] = "pacific_xsection_inversion:main";
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
  char               *bin_solver_path;
  char               *txt_inv_solver_path;
  char               *vtk_input_path;
  char               *vtk_solution_path;
  char               *vtk_solver_path;
  char               *vtk_inv_solver_vol_path;
  char               *vtk_inv_solver_surf_path;
  char               *vtk_inv_param_derivative_path;
  /* mesh */
  p4est_t                *p4est;
  ymir_mesh_t            *ymir_mesh;
  ymir_pressure_elem_t   *press_elem;
  /* Stokes */
  rhea_stokes_problem_t  *stokes_problem;
  /* Inversion */
  rhea_inversion_problem_t *inv_problem;

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

  /* binary file output */
  YMIR_OPTIONS_S, "bin-write-solver-path", '\0',
    &(bin_solver_path), NULL,
    "Bin file path for solver internals (e.g., iterations of Newton's method)",

  /* text file output */
  YMIR_OPTIONS_S, "txt-write-inversion-solver-path", '\0',
    &(txt_inv_solver_path), NULL,
    "Text file path for solver internals of the inversion",

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
  YMIR_OPTIONS_S, "vtk-write-inversion-solver-volume-path", '\0',
    &(vtk_inv_solver_vol_path), NULL,
    "VTK file path for solver internals of the inversion (volume fields).",
  YMIR_OPTIONS_S, "vtk-write-inversion-solver-surface-path", '\0',
    &(vtk_inv_solver_surf_path), NULL,
    "VTK file path for solver internals of the inversion (surface fields).",
  YMIR_OPTIONS_S, "vtk-write-inversion-param-derivative-path", '\0',
    &(vtk_inv_param_derivative_path), NULL,
    "VTK file path for testing derivatives of the viscosity w.r.t. parameters.",

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

  /* setup Stokes solver */
  rhea_performance_monitor_start_barrier (RHEA_MAIN_PERFMON_SETUP_STOKES);
  rhea_stokes_problem_setup_solver_ext (stokes_problem,
                                        RHEA_INVERSION_KRYLOV_SOLVER_N);
  rhea_performance_monitor_stop_add (RHEA_MAIN_PERFMON_SETUP_STOKES);

  /*
   * Solve Inverse Problem
   */

  /* setup inversion solver */
  inv_problem = rhea_inversion_new (stokes_problem);
  rhea_inversion_set_txt_output (inv_problem, txt_inv_solver_path);
  rhea_inversion_set_vtk_output (inv_problem, vtk_inv_solver_vol_path,
                                 vtk_inv_solver_surf_path);

  /* run inversion solver */
  rhea_performance_monitor_start_barrier (RHEA_MAIN_PERFMON_SOLVE_INVERSION);
  rhea_inversion_solve (inv_problem,
                        0 /* use model param's as initial guess */,
                        NULL /* no parameter vector */);
  rhea_performance_monitor_stop_add_barrier (RHEA_MAIN_PERFMON_SOLVE_INVERSION);

  /* destroy */
  rhea_inversion_destroy (inv_problem);

  /*
   * Clear Stokes Problem & Mesh
   */

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
                                  RHEA_PERFMON_PRINT_WTIME_ALL,
                                  RHEA_PERFMON_PRINT_NCALLS_ESSENTIAL,
                                  RHEA_PERFMON_PRINT_FLOPS_NONE,
                                  RHEA_PERFMON_PRINT_YMIR_NONE);
  rhea_performance_monitor_finalize ();

  /* destroy options */
  ymir_options_global_destroy ();

  /* print that this function is ending */
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", func_name);

  /* finalize rhea */
  rhea_finalize ();

  return 0;
}
