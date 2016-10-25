/*
 */

#include <rhea.h>

/**
 * Sets up the mesh.
 */
static void
basic_setup_mesh (p4est_t **p4est,
                  ymir_mesh_t **ymir_mesh,
                  ymir_pressure_elem_t **press_elem,
                  MPI_Comm mpicomm,
                  rhea_domain_options_t *domain_options,
                  rhea_discretization_options_t *discr_options)
{
  const char         *this_fn_name = "basic_setup_mesh";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* create p4est */
  *p4est = rhea_discretization_p4est_new (mpicomm, discr_options,
                                          domain_options);

  /* set up boundary, store in `discr_options` */
  rhea_discretization_options_set_boundary (discr_options, *p4est,
                                            domain_options);

  /* create ymir mesh and pressure element */
  rhea_discretization_ymir_mesh_new_from_p4est (ymir_mesh, press_elem, *p4est,
                                                discr_options);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Sets up a linear Stokes problem.
 */
static void
basic_setup_stokes (rhea_stokes_linear_problem_t **lin_stokes,
                    ymir_mesh_t *ymir_mesh,
                    ymir_pressure_elem_t *press_elem,
                    rhea_domain_options_t *domain_options,
                    rhea_viscosity_options_t *viscosity_options)
{
  const char         *this_fn_name = "basic_setup_stokes";
  ymir_vec_t         *temperature, *viscosity;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* create temperature */
  temperature = rhea_temperature_new (ymir_mesh);
  //TODO set temp
  ymir_vec_set_value (temperature, 0.5);

  /* create Stokes coefficient (= 2 * viscosity) */
  viscosity = rhea_viscosity_new (ymir_mesh);
  //TODO set visc
  ymir_vec_set_value (viscosity, 1.0);

  /* create Stokes problem */
  *lin_stokes = rhea_stokes_linear_problem_new (
      temperature, viscosity, 1 /* pass viscosity vector as owned */,
      ymir_mesh, press_elem, domain_options);

  /* destroy */
  rhea_temperature_destroy (temperature);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Cleans up Stokes problem and mesh.
 */
static void
basic_setup_clear_all (rhea_stokes_linear_problem_t *lin_stokes,
                       p4est_t *p4est,
                       ymir_mesh_t *ymir_mesh,
                       ymir_pressure_elem_t *press_elem,
                       rhea_discretization_options_t *discr_options)
{
  const char         *this_fn_name = "basic_setup_clear_all";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* destroy Stokes problem */
  rhea_stokes_linear_problem_destroy (lin_stokes);

  /* destroy mesh */
  rhea_discretization_ymir_mesh_destroy (ymir_mesh, press_elem);
  rhea_discretization_p4est_destroy (p4est);

  /* destroy (some) options */
  rhea_discretization_options_clear (discr_options);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Runs Stokes solver.
 */
static void
basic_run_solver ()
{
  const char         *this_fn_name = "basic_run_solver";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* initialize velocity and pressure */
//slabs_stokes_state_init_vel_press (state, mesh, press_elem);

  /* run solver */
//slabs_solve_stokes (lin_stokes, &nl_stokes, p8est, &mesh, &press_elem,
//                    state, &physics_options, &discr_options,
//                    &solver_options, workload_out_path,
//                    bin_nl_filepath, vtk_nl_filepath);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Runs the program.
 */
int
main (int argc, char **argv)
{
  const char         *this_fn_name = "basic:main";
  /* MPI */
  MPI_Comm            mpicomm = MPI_COMM_WORLD;
  int                 mpisize, mpirank, ompsize;
  int                 mpiret;
  /* options */
  ymir_options_t     *opt;
  rhea_domain_options_t         domain_options;
  rhea_discretization_options_t discr_options;
  rhea_viscosity_options_t      viscosity_options;
  /* options local to this function */
  int                 production_run;
  /* mesh */
  p4est_t            *p4est;
  ymir_mesh_t        *ymir_mesh;
  ymir_pressure_elem_t  *press_elem;
  /* Stokes */
  rhea_stokes_linear_problem_t *lin_stokes;

  /*
   * Initialize Libraries
   */

  /* initialize rhea and sub-packages */
  rhea_initialize (argc, argv, mpicomm);

  /* get parallel environment */
  mpiret = MPI_Comm_size (mpicomm, &mpisize); YMIR_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &mpirank); YMIR_CHECK_MPI (mpiret);

#ifdef RHEA_ENABLE_OPENMP
  #pragma omp parallel
  {
    #pragma omp single
    {
      ompsize = omp_get_num_threads ();
    }
  }
#else
  ompsize = 1;
#endif

  /*
   * Define & Parse Options
   */

  opt = ymir_options_global_new (argv[0] /* program path */);

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

    /* basic options */
  YMIR_OPTIONS_CALLBACK, "help", 'h', 0 /* no callback fn args */,
    ymir_options_print_usage_and_exit_fn, NULL /* no arg usage */,
    "Print usage and exit",
  YMIR_OPTIONS_INIFILE, "options-file", 'f',
    ".ini file with option values",

  /* performance & monitoring options */
  YMIR_OPTIONS_B, "production-run", '\0',
    &production_run, 0,
    "Execute as a production run (to reduce some overhead and checks)",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add sub-options */
  rhea_add_options_all (opt);
  ymir_options_add_suboptions_solver_stokes (opt);

  /* parse options */
  {
    int                 optret;

    optret = ymir_options_parse (SC_LP_INFO, opt, argc, argv);
    if (optret < 0) { /* if parsing was not successful */
      ymir_options_print_usage (SC_LP_INFO, opt, NULL /* args usage */);
      RHEA_GLOBAL_INFO ("Option parsing failed\n");
      exit (0);
    }
  }

  /*
   * Initialize Main Program
   */

  RHEA_GLOBAL_PRODUCTIONF (
      "Into %s (production %i)\n", this_fn_name, production_run);
  RHEA_GLOBAL_PRODUCTIONF (
      "Parallel environment: MPI size %i, OpenMP size %i\n", mpisize, ompsize);
  ymir_set_up (argc, argv, mpicomm, production_run);

  /* print & process options */
  ymir_options_print_summary (SC_LP_INFO, opt);
  rhea_process_options_all (&domain_options, &discr_options,
                            &viscosity_options);

  /*
   * Setup Mesh
   */

  basic_setup_mesh (&p4est, &ymir_mesh, &press_elem, mpicomm,
                    &domain_options, &discr_options);

  /*
   * Setup Stokes Problem
   */

  basic_setup_stokes (&lin_stokes, ymir_mesh, press_elem,
                      &domain_options, &viscosity_options);

  /*
   * Solve Stokes Problem
   */

//basic_run_solver ()

  /*
   * Finalize
   */

  /* destroy Stokes problem and mesh */
  basic_setup_clear_all (lin_stokes, p4est, ymir_mesh, press_elem,
                         &discr_options);

  /* destroy options */
  ymir_options_global_destroy ();

  /* print that this function is ending */
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

  /* finalize rhea */
  rhea_finalize ();

  return 0;
}
