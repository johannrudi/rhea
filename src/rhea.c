#include <rhea.h>
#include <rhea_amr.h>
#include <rhea_io_mpi.h>
#include <ymir.h>
#include <ymir_perf_counter.h>
#include <ymir_stress_pc.h>
#include <ymir_stokes_pc.h>
#include <ymir_bbt.h>
#include <ymir_bfbt.h>

/* default options */
#define RHEA_OPTDEF_PRODUCTION_RUN (0)
#define RHEA_OPTDEF_MONITOR_PERFORMANCE (0)
#define RHEA_OPTDEF_MONITOR_PERFORMANCE_MATVEC (0)

int                 rhea_opt_production_run = RHEA_OPTDEF_PRODUCTION_RUN;
int                 rhea_opt_monitor_performance =
                      RHEA_OPTDEF_MONITOR_PERFORMANCE;
int                 rhea_opt_monitor_performance_matvec =
                      RHEA_OPTDEF_MONITOR_PERFORMANCE_MATVEC;

/* initialization parameters */
int                 rhea_argc = 0;
char              **rhea_argv = NULL;
sc_MPI_Comm         rhea_mpicomm = sc_MPI_COMM_NULL;

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
rhea_init_end (ymir_options_t *options)
{
  const int           err_priority = SC_LP_INFO;
  int                 optret;

  /* check */
  RHEA_ASSERT (0 < rhea_argc);
  RHEA_ASSERT (rhea_argv != NULL);
  RHEA_ASSERT (rhea_mpicomm != sc_MPI_COMM_NULL);

  /* parse options */
  optret = ymir_options_parse (err_priority, options, rhea_argc, rhea_argv);
  if (optret < 0) { /* if parsing was not successful */
    ymir_options_print_usage (err_priority, options, NULL /* args usage */);
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
rhea_add_options_base (ymir_options_t *options)
{
  /* *INDENT-OFF* */
  ymir_options_addv (options,

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
  YMIR_OPTIONS_B, "monitor-performance-matvec", '\0',
    &(rhea_opt_monitor_performance_matvec),
    RHEA_OPTDEF_MONITOR_PERFORMANCE_MATVEC,
    "Measure and print performance statistics of all Matvecs",

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
  rhea_composition_add_options (options);
  rhea_plate_add_options (options);
  rhea_weakzone_add_options (options);
  rhea_topography_add_options (options);
  rhea_viscosity_add_options (options);
  rhea_discretization_add_options (options);
  rhea_stokes_problem_add_options (options);
  rhea_inversion_add_options (options);
  rhea_io_mpi_add_options (options);
}

void
rhea_add_options_newton (ymir_options_t *options)
{
  rhea_domain_add_options (options);
  rhea_discretization_add_options (options);
  rhea_newton_add_options (NULL, options);
}

void
rhea_process_options_all (rhea_all_options_t *all_options)
{
  rhea_domain_process_options (all_options->domain_options);
  rhea_temperature_process_options (
      all_options->temperature_options,
      all_options->domain_options);
  rhea_composition_process_options (
      all_options->composition_options,
      all_options->domain_options);
  rhea_plate_process_options (
      all_options->plate_options,
      all_options->domain_options,
      all_options->temperature_options);
  rhea_weakzone_process_options (
      all_options->weakzone_options,
      all_options->domain_options);
  rhea_topography_process_options (
      all_options->topography_options,
      all_options->domain_options);
  rhea_viscosity_process_options (
      all_options->viscosity_options,
      all_options->domain_options,
      all_options->temperature_options);
  rhea_discretization_process_options (
      all_options->discr_options,
      all_options->domain_options,
      all_options->topography_options);
  rhea_stokes_problem_process_options ();

  /* print */
  rhea_print_physics_const_options (
      all_options->domain_options,
      all_options->temperature_options,
      all_options->viscosity_options);
}

void
rhea_process_options_newton (rhea_domain_options_t *domain_options,
                             rhea_discretization_options_t *discr_options,
                             rhea_newton_options_t *newton_options)
{
  rhea_domain_process_options (domain_options);
  rhea_discretization_process_options (discr_options, domain_options, NULL);
  rhea_newton_get_global_options (newton_options);

  /* print */
  rhea_domain_print_const_options (domain_options);
}

void
rhea_print_physics_const_options (
                              rhea_domain_options_t *domain_options,
                              rhea_temperature_options_t *temperature_options,
                              rhea_viscosity_options_t *viscosity_options)
{
  const double        rho   = domain_options->density_kg_m3;
  const double        g     = domain_options->gravity_m_s2;
  const double        b     = domain_options->radius_max_m;
  const double        alpha = temperature_options->thermal_expansivity_1_k;
  const double        kappa = temperature_options->thermal_diffusivity_m2_s;
  const double        dT    = temperature_options->temperature_difference_K;
  const double        R     = temperature_options->gas_constant_J_Kmol;
  const double        mu    = viscosity_options->representative_Pas;
  const double        rayleigh = (alpha*dT*rho*g*b*b*b) / (mu*kappa);
  const double        rhs_scal = temperature_options->rhs_scaling;

  /* print domain options */
  rhea_domain_print_const_options (domain_options);

  /* print options regarding nondimensionalization */
  RHEA_GLOBAL_INFO ("========================================\n");
  RHEA_GLOBAL_INFOF ("%s\n", __func__);
  RHEA_GLOBAL_INFO ("----------------------------------------\n");
  RHEA_GLOBAL_INFOF ("  thermal expansivity        [K^-1]    : %g\n", alpha);
  RHEA_GLOBAL_INFOF ("  density                    [kg/m^3]  : %g\n", rho);
  RHEA_GLOBAL_INFOF ("  temperature difference     [K]       : %g\n", dT);
  RHEA_GLOBAL_INFOF ("  thermal diffusivity        [m^2/s]   : %g\n", kappa);
  RHEA_GLOBAL_INFOF ("  gas constant               [J/K/mol] : %g\n", R);
  RHEA_GLOBAL_INFOF ("  representative viscosity   [Pa*s]    : %g\n", mu);
  RHEA_GLOBAL_INFOF ("  gravitational acceleration [m/s^2]   : %g\n", g);
  RHEA_GLOBAL_INFOF ("  depth of convecting region [m]       : %g\n", b);
  RHEA_GLOBAL_INFO ("----------------------------------------\n");
  RHEA_GLOBAL_INFOF ("  Rayleigh number                      : %g\n", rayleigh);
  RHEA_GLOBAL_INFOF ("  RHS scaling from temperature         : %g\n", rhs_scal);
  RHEA_GLOBAL_INFO ("========================================\n");
}

/******************************************************************************
 * Monitoring
 *****************************************************************************/

/* main perfomance monitors */
ymir_perf_counter_t  *rhea_perfmon_main = NULL;
int                   rhea_perfmon_main_n = 0;

/* matvec perfomance monitor tags and names */
typedef enum
{
  RHEA_PERFMON_MATVEC_POISSON,
  RHEA_PERFMON_MATVEC_STRESS,
  RHEA_PERFMON_MATVEC_STOKES_GRAD_DIV,
  RHEA_PERFMON_MATVEC_POISSON_PRESS,
  RHEA_PERFMON_MATVEC_STOKES,
  RHEA_PERFMON_MATVEC_TOTAL,
  RHEA_PERFMON_MATVEC_N
}
rhea_perfmon_matvec_idx_t;

static const char  *rhea_perfmon_matvec_name[RHEA_PERFMON_MATVEC_N] =
{
  "Matvec: Poisson",
  "Matvec: Viscous Stress",
  "Matvec: B and B^T",
  "Matvec: BB^T",
  "Matvec: Stokes",
  "Matvec (total)"
};
ymir_perf_counter_t rhea_perfmon_matvec[RHEA_PERFMON_MATVEC_N];

static void
rhea_perfmon_matvec_init (const int activate, const int skip_if_active)
{
  const int           active = activate && rhea_opt_monitor_performance_matvec;

  ymir_perf_counter_init_all_ext (rhea_perfmon_matvec, rhea_perfmon_matvec_name,
                                  RHEA_PERFMON_MATVEC_N, active,
                                  skip_if_active);
}

static void
rhea_perfmon_matvec_print (sc_MPI_Comm mpicomm,
                           const int print_wtime,
                           const int print_n_calls,
                           const int print_flops)
{
  const int           active = rhea_opt_monitor_performance_matvec;
  const int           print = (print_wtime || print_n_calls || print_flops);
  int                 n_stats = RHEA_PERFMON_MATVEC_N *
                                YMIR_PERF_COUNTER_N_STATS;
  sc_statinfo_t       stats[n_stats];
  char                stats_name[n_stats][YMIR_PERF_COUNTER_NAME_SIZE];

  /* exit if nothing to do */
  if (!active || !print) {
    return;
  }

  /* get each matvec */
  ymir_stiff_op_perf_counter_add_matvecs (
      &rhea_perfmon_matvec[RHEA_PERFMON_MATVEC_POISSON]);
  ymir_stress_op_perf_counter_add_matvecs (
      &rhea_perfmon_matvec[RHEA_PERFMON_MATVEC_STRESS]);
  ymir_pressure_vec_perf_counter_add_matvecs (
      &rhea_perfmon_matvec[RHEA_PERFMON_MATVEC_STOKES_GRAD_DIV]);
  ymir_bbt_perf_counter_add_matvecs (
      &rhea_perfmon_matvec[RHEA_PERFMON_MATVEC_POISSON_PRESS]);
  ymir_stokes_op_perf_counter_add_matvecs (
      &rhea_perfmon_matvec[RHEA_PERFMON_MATVEC_STOKES]);

  /* accumulate total from all matvecs */
  ymir_stiff_op_perf_counter_add_matvecs (
      &rhea_perfmon_matvec[RHEA_PERFMON_MATVEC_TOTAL]);
  ymir_stress_op_perf_counter_add_matvecs (
      &rhea_perfmon_matvec[RHEA_PERFMON_MATVEC_TOTAL]);
  ymir_pressure_vec_perf_counter_add_matvecs (
      &rhea_perfmon_matvec[RHEA_PERFMON_MATVEC_TOTAL]);
  ymir_bbt_perf_counter_add_matvecs (
      &rhea_perfmon_matvec[RHEA_PERFMON_MATVEC_TOTAL]);
  ymir_stokes_op_perf_counter_add_matvecs (
      &rhea_perfmon_matvec[RHEA_PERFMON_MATVEC_TOTAL]);

  /* gather performance statistics */
  n_stats = ymir_perf_counter_gather_stats (
      rhea_perfmon_matvec, RHEA_PERFMON_MATVEC_N,
      stats, stats_name, rhea_mpicomm,
      print_wtime, print_n_calls, print_flops);

  /* print performance statistics */
  ymir_perf_counter_print_stats (stats, n_stats, "Matvecs");
}

int
rhea_performance_monitor_active ()
{
  return rhea_opt_monitor_performance;
}

void
rhea_performance_monitor_init (const char **monitor_name,
                               const int n_monitors)
{
  const int           active = rhea_performance_monitor_active ();

  RHEA_ASSERT (rhea_perfmon_main == NULL);
  RHEA_ASSERT (rhea_perfmon_main_n == 0);

  /* create and initialize main performance monitors */
  rhea_perfmon_main = RHEA_ALLOC (ymir_perf_counter_t, n_monitors);
  rhea_perfmon_main_n = n_monitors;
  ymir_perf_counter_init_all (rhea_perfmon_main, monitor_name,
                              rhea_perfmon_main_n, active);

  /* initialize rhea's internal performance monitors */
  rhea_io_mpi_perfmon_init (active, 0);
  rhea_amr_perfmon_init (active, 0);
  rhea_plate_perfmon_init (active, 0);
  rhea_weakzone_perfmon_init (active, 0);
  rhea_stokes_problem_perfmon_init (active, 0);
  rhea_inversion_perfmon_init (active, 0);
  rhea_perfmon_matvec_init (active, 0);
}

void
rhea_performance_monitor_finalize ()
{
  RHEA_ASSERT (rhea_perfmon_main != NULL);
  RHEA_ASSERT (0 < rhea_perfmon_main_n);

  RHEA_FREE (rhea_perfmon_main);
  rhea_perfmon_main = NULL;
  rhea_perfmon_main_n = 0;
}

void
rhea_performance_monitor_print (
                        const char *title,
                        rhea_performance_monitor_print_wtime_t print_wtime,
                        rhea_performance_monitor_print_ncalls_t print_n_calls,
                        rhea_performance_monitor_print_flops_t print_flops,
                        rhea_performance_monitor_print_ymir_t print_ymir)
{
  const int           active = rhea_performance_monitor_active ();
  const int           print = (print_wtime || print_n_calls || print_flops);
  int                 n_stats = rhea_perfmon_main_n *
                                YMIR_PERF_COUNTER_N_STATS;
  sc_statinfo_t       stats[n_stats];
  char                stats_name[n_stats][YMIR_PERF_COUNTER_NAME_SIZE];

  /* exit if nothing to do */
  if (!active || !(print || print_ymir)) {
    return;
  }

  RHEA_ASSERT (rhea_perfmon_main != NULL);
  RHEA_ASSERT (0 < rhea_perfmon_main_n);

  /* print rhea's setup performance statistics */
  rhea_io_mpi_perfmon_print (
      rhea_mpicomm,
      RHEA_PERFMON_PRINT_WTIME_ESSENTIAL <= print_wtime,
      RHEA_PERFMON_PRINT_NCALLS_ESSENTIAL <= print_n_calls,
      RHEA_PERFMON_PRINT_FLOPS_ALL <= print_flops);
  rhea_amr_perfmon_print (
      rhea_mpicomm,
      RHEA_PERFMON_PRINT_WTIME_ESSENTIAL <= print_wtime,
      RHEA_PERFMON_PRINT_NCALLS_ESSENTIAL <= print_n_calls,
      RHEA_PERFMON_PRINT_FLOPS_ALL <= print_flops);
  rhea_plate_perfmon_print (
      rhea_mpicomm,
      RHEA_PERFMON_PRINT_WTIME_ESSENTIAL <= print_wtime,
      RHEA_PERFMON_PRINT_NCALLS_ESSENTIAL <= print_n_calls,
      RHEA_PERFMON_PRINT_FLOPS_ALL <= print_flops);
  rhea_weakzone_perfmon_print (
      rhea_mpicomm,
      RHEA_PERFMON_PRINT_WTIME_ESSENTIAL <= print_wtime,
      RHEA_PERFMON_PRINT_NCALLS_ESSENTIAL <= print_n_calls,
      RHEA_PERFMON_PRINT_FLOPS_ALL <= print_flops);
  rhea_stokes_problem_perfmon_print (
      rhea_mpicomm,
      RHEA_PERFMON_PRINT_WTIME_ESSENTIAL <= print_wtime,
      RHEA_PERFMON_PRINT_NCALLS_ESSENTIAL <= print_n_calls,
      RHEA_PERFMON_PRINT_FLOPS_ALL <= print_flops);
  rhea_inversion_perfmon_print (
      rhea_mpicomm,
      RHEA_PERFMON_PRINT_WTIME_ESSENTIAL <= print_wtime,
      RHEA_PERFMON_PRINT_NCALLS_ESSENTIAL <= print_n_calls,
      RHEA_PERFMON_PRINT_FLOPS_ALL <= print_flops);

  /* print ymir performance statistics */
  if (RHEA_PERFMON_PRINT_YMIR_ALL == print_ymir) {
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

  /* print matvec performance statistics */
  rhea_perfmon_matvec_print (
      rhea_mpicomm,
      RHEA_PERFMON_PRINT_WTIME_ESSENTIAL <= print_wtime,
      RHEA_PERFMON_PRINT_NCALLS_ESSENTIAL <= print_n_calls,
      RHEA_PERFMON_PRINT_FLOPS_ALL <= print_flops);

  /* gather main performance statistics */
  n_stats = ymir_perf_counter_gather_stats (
      rhea_perfmon_main, rhea_perfmon_main_n, stats, stats_name, rhea_mpicomm,
      RHEA_PERFMON_PRINT_WTIME_ESSENTIAL <= print_wtime,
      RHEA_PERFMON_PRINT_NCALLS_ALL <= print_n_calls,
      RHEA_PERFMON_PRINT_FLOPS_ALL <= print_flops);

  /* print main performance statistics */
  ymir_perf_counter_print_stats (stats, n_stats, title);
}

void
rhea_performance_monitor_start (const int monitor_index)
{
  RHEA_ASSERT (rhea_perfmon_main != NULL);
  RHEA_ASSERT (monitor_index < rhea_perfmon_main_n);

  ymir_perf_counter_start (&rhea_perfmon_main[monitor_index]);
}

void
rhea_performance_monitor_stop_add (const int monitor_index)
{
  RHEA_ASSERT (rhea_perfmon_main != NULL);
  RHEA_ASSERT (monitor_index < rhea_perfmon_main_n);

  ymir_perf_counter_stop_add (&rhea_perfmon_main[monitor_index]);
}

void
rhea_performance_monitor_start_barrier (const int monitor_index)
{
  RHEA_ASSERT (rhea_perfmon_main != NULL);
  RHEA_ASSERT (monitor_index < rhea_perfmon_main_n);

  ymir_perf_counter_start_barrier (&rhea_perfmon_main[monitor_index],
                                   rhea_mpicomm);
}

void
rhea_performance_monitor_stop_add_barrier (const int monitor_index)
{
  RHEA_ASSERT (rhea_perfmon_main != NULL);
  RHEA_ASSERT (monitor_index < rhea_perfmon_main_n);

  ymir_perf_counter_stop_add_barrier (&rhea_perfmon_main[monitor_index],
                                      rhea_mpicomm);
}
