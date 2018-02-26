/*
 */

#include <rhea.h>
#include <rhea_io_mpi.h>
#include <rhea_amr.h>
#include <ymir.h>
#include <ymir_perf_counter.h>
#include <ymir_stress_pc.h>
#include <ymir_stokes_pc.h>
#include <ymir_bbt.h>
#include <ymir_bfbt.h>

/* default options */
#define RHEA_OPTDEF_PRODUCTION_RUN (0)
#define RHEA_OPTDEF_MONITOR_PERFORMANCE (0)

/* initialization parameters */
int                 rhea_argc = 0;
char              **rhea_argv = NULL;
sc_MPI_Comm         rhea_mpicomm = sc_MPI_COMM_NULL;
int                 rhea_opt_production_run = RHEA_OPTDEF_PRODUCTION_RUN;
int                 rhea_opt_monitor_performance =
                      RHEA_OPTDEF_MONITOR_PERFORMANCE;

void
rhea_init_begin (int *mpisize, int *mpirank, int *ompsize,
                 int argc, char **argv, sc_MPI_Comm mpicomm)
{
  int                 mpiret;

  /* set parameters */
  rhea_argc = argc;
  rhea_argv = argv;
  rhea_mpicomm = mpicomm;

  /* initialize rhea and sub-packages */
  rhea_initialize (argc, argv, mpicomm);

  /* get MPI environment */
  if (mpisize != NULL) {
    mpiret = sc_MPI_Comm_size (mpicomm, mpisize); SC_CHECK_MPI (mpiret);
  }
  if (mpirank != NULL) {
    mpiret = sc_MPI_Comm_rank (mpicomm, mpirank); SC_CHECK_MPI (mpiret);
  }

  /* get OpenMP environment */
  if (ompsize != NULL) {
#ifdef RHEA_ENABLE_OPENMP
    *ompsize = omp_get_max_threads ();
#else
    *ompsize = 1;
#endif
  }
}

void
rhea_init_end (ymir_options_t *opt)
{
  const int           err_priority = SC_LP_INFO;
  int                 optret;

  /* check */
  RHEA_ASSERT (0 < rhea_argc);
  RHEA_ASSERT (rhea_argv != NULL);
  RHEA_ASSERT (rhea_mpicomm != sc_MPI_COMM_NULL);

  /* parse options */
  optret = ymir_options_parse (err_priority, opt, rhea_argc, rhea_argv);
  if (optret < 0) { /* if parsing was not successful */
    ymir_options_print_usage (err_priority, opt, NULL /* args usage */);
    RHEA_GLOBAL_PRODUCTION ("Option parsing failed\n");
    exit (0);
  }

  /* set up ymir */
  ymir_set_up (rhea_argc, rhea_argv, rhea_mpicomm, rhea_opt_production_run);
}

int
rhea_production_run_get ()
{
  RHEA_ASSERT (rhea_opt_production_run == ymir_get_production_run ());
  return rhea_opt_production_run;
}

void
rhea_production_run_set (const int is_production_run)
{
  if (is_production_run) {
    rhea_opt_production_run = 1;
  }
  else {
    rhea_opt_production_run = 0;
  }
  ymir_set_production_run (rhea_opt_production_run);
}

/******************************************************************************
 * Options
 *****************************************************************************/

void
rhea_add_options_base (ymir_options_t *opt)
{
  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  /* help message */
  YMIR_OPTIONS_CALLBACK, "help", 'h', 0 /* no callback fn args */,
    ymir_options_print_usage_and_exit_fn, NULL /* no arg usage */,
    "Print usage and exit",

  /* options file */
  YMIR_OPTIONS_INIFILE, "options-file", 'f',
    ".ini file with option values",

  /* performance */
  YMIR_OPTIONS_B, "production-run", 'p',
    &(rhea_opt_production_run), RHEA_OPTDEF_PRODUCTION_RUN,
    "Execute as a production run (to reduce some overhead and checks)",

  /* monitoring */
  YMIR_OPTIONS_B, "monitor-performance", 'm',
    &(rhea_opt_monitor_performance), RHEA_OPTDEF_MONITOR_PERFORMANCE,
    "Measure and print performance statistics (e.g., runtime or flops)",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */
}

void
rhea_add_options_all (ymir_options_t *options)
{
  rhea_domain_add_options (options);
  rhea_temperature_add_options (options);
  rhea_temperature_add_options_sinker (options);
  rhea_temperature_add_options_plume (options);
  rhea_weakzone_add_options (options);
  rhea_viscosity_add_options (options);
  rhea_topography_add_options (options);
  rhea_discretization_add_options (options);
  rhea_stokes_problem_add_options (options);
  rhea_newton_add_options (options);
  rhea_io_mpi_add_options (options);
}

void
rhea_add_options_newton (ymir_options_t *options)
{
  rhea_domain_add_options (options);
  rhea_discretization_add_options (options);
  rhea_newton_add_options (options);
}

void
rhea_process_options_all (rhea_domain_options_t *domain_options,
                          rhea_temperature_options_t *temperature_options,
                          rhea_weakzone_options_t *weakzone_options,
                          rhea_viscosity_options_t *viscosity_options,
                          rhea_topography_options_t *topography_options,
                          rhea_discretization_options_t *discr_options,
                          rhea_newton_options_t *newton_options)
{
  rhea_domain_process_options (domain_options);
  rhea_temperature_process_options (temperature_options, domain_options);
  rhea_weakzone_process_options (weakzone_options, domain_options);
  rhea_viscosity_process_options (viscosity_options, domain_options);
  rhea_topography_process_options (topography_options, domain_options);
  rhea_discretization_process_options (discr_options, domain_options,
                                       topography_options);
  rhea_newton_process_options (newton_options);
}

void
rhea_process_options_newton (rhea_domain_options_t *domain_options,
                             rhea_discretization_options_t *discr_options,
                             rhea_newton_options_t *newton_options)
{
  rhea_domain_process_options (domain_options);
  rhea_discretization_process_options (discr_options, domain_options, NULL);
  rhea_newton_process_options (newton_options);
}

/******************************************************************************
 * Monitoring
 *****************************************************************************/

int
rhea_performance_monitor_active ()
{
  return rhea_opt_monitor_performance;
}

/* perfomance monitors */
ymir_perf_counter_t  *rhea_performance_monitor = NULL;
int                   rhea_performance_monitor_n = 0;

void
rhea_performance_monitor_init (const char **monitor_name,
                               const int n_monitors)
{
  const int           active = rhea_performance_monitor_active ();

  RHEA_ASSERT (rhea_performance_monitor == NULL);
  RHEA_ASSERT (rhea_performance_monitor_n == 0);

  /* create and initialize main performance monitors */
  rhea_performance_monitor = RHEA_ALLOC (ymir_perf_counter_t, n_monitors);
  rhea_performance_monitor_n = n_monitors;
  ymir_perf_counter_init_all (rhea_performance_monitor, monitor_name,
                              rhea_performance_monitor_n, active);

  /* initialize rhea performance monitors */
  rhea_io_mpi_perfmon_init (active, 0);
  rhea_amr_perfmon_init (active, 0);
}

void
rhea_performance_monitor_finalize ()
{
  RHEA_ASSERT (rhea_performance_monitor != NULL);
  RHEA_ASSERT (0 < rhea_performance_monitor_n);

  RHEA_FREE (rhea_performance_monitor);
  rhea_performance_monitor = NULL;
  rhea_performance_monitor_n = 0;
}

void
rhea_performance_monitor_print (const char *title,
                                const int print_wtime,
                                const int print_n_calls,
                                const int print_flops,
                                const int print_ymir)
{
  const int           active = rhea_performance_monitor_active ();
  int                 n_stats = rhea_performance_monitor_n *
                                YMIR_PERF_COUNTER_N_STATS;
  sc_statinfo_t       stats[n_stats];
  char                stats_name[n_stats][YMIR_PERF_COUNTER_NAME_SIZE];

  /* exit if nothing to do */
  if (!active) {
    return;
  }

  RHEA_ASSERT (rhea_performance_monitor != NULL);
  RHEA_ASSERT (0 < rhea_performance_monitor_n);

  /* print ymir performance statistics */
  if (print_ymir) {
    ymir_gmg_hierarchy_mesh_perf_counter_print ();   /* GMG mesh */

    ymir_stiff_op_perf_counter_print ();             /* Stiffness Op */
    ymir_stiff_pc_perf_counter_print ();             /* Stiffness PC */
    ymir_gmg_hierarchy_stiff_perf_counter_print ();  /* GMG stiffness */

    ymir_stress_op_perf_counter_print ();            /* Stress Op */
    ymir_stress_pc_perf_counter_print ();            /* Stress PC */
    ymir_gmg_hierarchy_stress_perf_counter_print (); /* GMG stress */

    ymir_pressure_vec_perf_counter_print ();         /* B^T and B */
    ymir_bbt_perf_counter_print ();                  /* BB^T */
    ymir_bfbt_perf_counter_print ();                 /* BFBT */

    ymir_stokes_op_perf_counter_print ();            /* Stokes Op */
    ymir_stokes_pc_perf_counter_print ();            /* Stokes PC */
  }

  /* print rhea performance statistics */
  rhea_io_mpi_perfmon_print (rhea_mpicomm, print_wtime, print_n_calls,
                             print_flops);
  rhea_amr_perfmon_print (rhea_mpicomm, print_wtime, print_n_calls,
                          print_flops);

  /* gather main performance statistics */
  n_stats = ymir_perf_counter_gather_stats (
      rhea_performance_monitor, rhea_performance_monitor_n,
      stats, stats_name, rhea_mpicomm,
      print_wtime, print_n_calls, print_flops);

  /* print main performance statistics */
  ymir_perf_counter_print_stats (stats, n_stats, title);
}

void
rhea_performance_monitor_start (const int monitor_index)
{
  RHEA_ASSERT (rhea_performance_monitor != NULL);
  RHEA_ASSERT (monitor_index < rhea_performance_monitor_n);

  ymir_perf_counter_start (&rhea_performance_monitor[monitor_index]);
}

void
rhea_performance_monitor_stop_add (const int monitor_index)
{
  RHEA_ASSERT (rhea_performance_monitor != NULL);
  RHEA_ASSERT (monitor_index < rhea_performance_monitor_n);

  ymir_perf_counter_stop_add (&rhea_performance_monitor[monitor_index]);
}

void
rhea_performance_monitor_start_barrier (const int monitor_index)
{
  RHEA_ASSERT (rhea_performance_monitor != NULL);
  RHEA_ASSERT (monitor_index < rhea_performance_monitor_n);

  ymir_perf_counter_start_barrier (&rhea_performance_monitor[monitor_index],
                                   rhea_mpicomm);
}

void
rhea_performance_monitor_stop_add_barrier (const int monitor_index)
{
  RHEA_ASSERT (rhea_performance_monitor != NULL);
  RHEA_ASSERT (monitor_index < rhea_performance_monitor_n);

  ymir_perf_counter_stop_add_barrier (&rhea_performance_monitor[monitor_index],
                                      rhea_mpicomm);
}
