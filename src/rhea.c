/*
 */

#include <rhea.h>
#include <ymir.h>

/* default options */
#define RHEA_DEFAULT_PRODUCTION_RUN (0)

/* initialization parameters */
int                 rhea_argc = 0;
char              **rhea_argv = NULL;
sc_MPI_Comm         rhea_mpicomm = sc_MPI_COMM_NULL;
int                 rhea_production_run = RHEA_DEFAULT_PRODUCTION_RUN;

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
  ymir_set_up (rhea_argc, rhea_argv, rhea_mpicomm, rhea_production_run);
}

int
rhea_get_production_run ()
{
  RHEA_ASSERT (rhea_production_run == ymir_get_production_run ());
  return rhea_production_run;
}

void
rhea_set_production_run (const int is_production_run)
{
  rhea_production_run = is_production_run;
  ymir_set_production_run (is_production_run);;
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
    &(rhea_production_run), 0,
    "Execute as a production run (to reduce some overhead and checks)",

  /* monitoring */
  //TODO

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
  rhea_discretization_add_options (options);
  rhea_amr_add_options (options);
  rhea_stokes_problem_add_options (options);
  rhea_newton_add_options (options);
}

void
rhea_add_options_newton (ymir_options_t *options)
{
  rhea_domain_add_options (options);
  rhea_discretization_add_options (options);
  rhea_amr_add_options (options);
  rhea_newton_add_options (options);
}

void
rhea_process_options_all (rhea_domain_options_t *domain_options,
                          rhea_temperature_options_t *temperature_options,
                          rhea_weakzone_options_t *weakzone_options,
                          rhea_viscosity_options_t *viscosity_options,
                          rhea_discretization_options_t *discr_options,
                          rhea_newton_options_t *newton_options)
{
  rhea_domain_process_options (domain_options);
  rhea_temperature_process_options (temperature_options, domain_options);
  rhea_weakzone_process_options (weakzone_options);
  rhea_viscosity_process_options (viscosity_options, domain_options);
  rhea_discretization_process_options (discr_options, domain_options);
  rhea_newton_process_options (newton_options);
}

void
rhea_process_options_newton (rhea_domain_options_t *domain_options,
                             rhea_discretization_options_t *discr_options,
                             rhea_newton_options_t *newton_options)
{
  rhea_domain_process_options (domain_options);
  rhea_discretization_process_options (discr_options, domain_options);
  rhea_newton_process_options (newton_options);
}
