/** EARTH
 *
 * Runs global mantle convection.
 *
 ******************************************************************************
 * Author:             Johann Rudi <johann@ices.utexas.edu>
 *****************************************************************************/

#include <rhea.h>
#include <example_share_mesh.h>
#include <example_share_stokes.h>
#include <example_share_vtk.h>
#include <ymir_perf_counter.h>
#include <ymir_monitor.h>

/******************************************************************************
 * Monitoring
 *****************************************************************************/

/* perfomance counters */
typedef enum
{
  EARTH_PERF_COUNTER_SETUP_MESH,
  EARTH_PERF_COUNTER_SETUP_STOKES,
  EARTH_PERF_COUNTER_SETUP_SOLVER,
  EARTH_PERF_COUNTER_SOLVE,
//EARTH_PERF_COUNTER_MATVECS,
  EARTH_PERF_COUNTER_TOTAL,
  EARTH_PERF_COUNTER_N
}
earth_perf_counter_idx_t;
ymir_perf_counter_t earth_perf_counter[EARTH_PERF_COUNTER_N];
const char         *earth_perf_counter_name[EARTH_PERF_COUNTER_N] =
{
  "Setup Mesh",
  "Setup Stokes",
  "Setup Solver",
  "Solve",
//"All matvecs (setup + solve)",
  "Total"
};
sc_statinfo_t       earth_perf_stats[
                      EARTH_PERF_COUNTER_N * YMIR_PERF_COUNTER_N_STATS];
char                earth_perf_stats_name[
                      EARTH_PERF_COUNTER_N * YMIR_PERF_COUNTER_N_STATS][
                      YMIR_PERF_COUNTER_NAME_SIZE];
int                 earth_perf_n_stats;

/**
 * Initializes performance counters.
 */
static void
earth_perf_counter_init (const int active)
{
  ymir_perf_counter_init_all (earth_perf_counter, earth_perf_counter_name,
                              EARTH_PERF_COUNTER_N, active);
}

/**
 * Gathers statistics of performance counters.
 */
static void
earth_perf_counter_gather (MPI_Comm mpicomm,
                           const int print_wtime,
                           const int print_n_calls,
                           const int print_flops)
{
  earth_perf_n_stats = ymir_perf_counter_gather_stats (
      earth_perf_counter, EARTH_PERF_COUNTER_N,
      earth_perf_stats, earth_perf_stats_name,
      mpicomm, print_wtime, print_n_calls, print_flops);
}

/**
 * Prints statistics of performance counters.
 */
static void
earth_perf_counter_print ()
{
  ymir_perf_counter_print_stats (earth_perf_stats, earth_perf_n_stats,
                                 "Earth Simulation");
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
  MPI_Comm            mpicomm = MPI_COMM_WORLD;
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

  /* vtk output options */
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

  /* initialize performance counters */
  earth_perf_counter_init (rhea_get_monitor_performance ());

  /* start performance counters */
  ymir_perf_counter_start_barrier (
      &earth_perf_counter[EARTH_PERF_COUNTER_TOTAL], mpicomm);

  /*
   * Print Environment and Options
   */

  RHEA_GLOBAL_PRODUCTIONF (
      "Into %s (production %i)\n", func_name, rhea_get_production_run ());
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

  ymir_perf_counter_start_barrier (
      &earth_perf_counter[EARTH_PERF_COUNTER_SETUP_MESH], mpicomm);
  example_share_mesh_new (&p4est, &ymir_mesh, &press_elem, mpicomm,
                          &domain_options, &topo_options, &discr_options);
  ymir_perf_counter_stop_add (
      &earth_perf_counter[EARTH_PERF_COUNTER_SETUP_MESH]);

  /*
   * Setup Stokes Problem
   */

  ymir_perf_counter_start_barrier (
      &earth_perf_counter[EARTH_PERF_COUNTER_SETUP_STOKES], mpicomm);
  example_share_stokes_new (&stokes_problem, &ymir_mesh, &press_elem,
                            &temp_options, &weak_options, &visc_options,
                            &newton_options, p4est, &discr_options,
                            vtk_solver_path);
  ymir_perf_counter_stop_add (
      &earth_perf_counter[EARTH_PERF_COUNTER_SETUP_STOKES]);

  /* write vtk of input data */
  example_share_vtk_write_input_data (vtk_input_path, stokes_problem,
                                      &temp_options, &visc_options);

  /*
   * Solve Stokes Problem
   */

  /* setup solver */
  ymir_perf_counter_start_barrier (
      &earth_perf_counter[EARTH_PERF_COUNTER_SETUP_SOLVER], mpicomm);
  rhea_stokes_problem_setup_solver (stokes_problem);
  ymir_perf_counter_stop_add (
      &earth_perf_counter[EARTH_PERF_COUNTER_SETUP_SOLVER]);

  /* initialize solution vector */
  sol_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);

  /* run solver */
  ymir_perf_counter_start_barrier (
      &earth_perf_counter[EARTH_PERF_COUNTER_SOLVE], mpicomm);
  rhea_stokes_problem_solve (&sol_vel_press, solver_iter_max, solver_rel_tol,
                             stokes_problem);
  ymir_perf_counter_stop_add (
      &earth_perf_counter[EARTH_PERF_COUNTER_SOLVE]);

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

  /* stop performance counters */
  ymir_perf_counter_stop_add_barrier (
      &earth_perf_counter[EARTH_PERF_COUNTER_TOTAL], mpicomm);

  /* print performance statistics */
  if (rhea_get_monitor_performance ()) {
    /* print ymir performance statistics */
  //ymir_gmg_hierarchy_mesh_perf_counter_print ();   /* GMG mesh */
  //ymir_stress_op_perf_counter_print ();            /* Stress Op */
  //ymir_stress_pc_perf_counter_print ();            /* Stress PC */
  //ymir_gmg_hierarchy_stress_perf_counter_print (); /* GMG stress */
  //ymir_stiff_op_perf_counter_print ();             /* Stiffness Op */
  //ymir_stiff_pc_perf_counter_print ();             /* Stiffness PC */
  //ymir_gmg_hierarchy_stiff_perf_counter_print ();  /* GMG stiffness */
  //ymir_pressure_vec_perf_counter_print ();         /* B^T or B */
  //ymir_bbt_perf_counter_print ();                  /* BB^T */
  //ymir_bfbt_perf_counter_print ();                 /* BFBT */
  //ymir_stokes_op_perf_counter_print ();            /* Stokes Op */
  //ymir_stokes_pc_perf_counter_print ();            /* Stokes PC */

    /* gather & print main performance statistics */
    earth_perf_counter_gather (mpicomm, 1 /* wtime */, 0 /* #calls */,
                               0 /* flops */);
    earth_perf_counter_print ();
  }

  /* destroy options */
  ymir_options_global_destroy ();

  /* print that this function is ending */
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", func_name);

  /* finalize rhea */
  rhea_finalize ();

  return 0;
}
