/* Collide Example:
 *
 * Cartesian domain.  A constant viscosity is used with discontinuous jumps.
 * The viscosity is defined with the coordinate.  For example, in the lower
 * mantle, visc=100, and in the upper mantle, the viscosity is 1 except a weak
 * zone in the center where the vsic is 0.1.  The temperature should also be
 * defined with the location.  Same thing for the boundary condition.  The
 * normal velocity boundary condition is prefixed with the location.  For
 * example, on the left side, there can be inflow at the upper half and outflow
 * in the lower half.
 */

#include <rhea.h>

/* enumerator for domain shapes */
typedef enum
{
  COLLIDE_VEL_DIR_BC_WALLSLIDE,
  COLLIDE_VEL_DIR_BC_INOUTFLOW
}
collide_vel_dir_bc_t;

/* options of collide example */
typedef struct collide_options
{
  collide_vel_dir_bc_t vel_dir_bc;
  double              wallslide_vel;
}
collide_options_t;

/**
 * Sets up the mesh.
 */
static void
collide_setup_mesh (p4est_t **p4est,
                    ymir_mesh_t **ymir_mesh,
                    ymir_pressure_elem_t **press_elem,
                    MPI_Comm mpicomm,
                    rhea_domain_options_t *domain_options,
                    rhea_discretization_options_t *discr_options)
{
  const char         *this_fn_name = "collide_setup_mesh";

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
 * Sets tangential velocity on one side of the domain.
 */
static ymir_dir_code_t
collide_set_vel_dir_wallslide (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
  if (face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE3 && 0.5 < Z) {
    return YMIR_VEL_DIRICHLET_ALL;
  }
  else {
    return YMIR_VEL_DIRICHLET_NORM;
  }
}

void
collide_set_rhs_vel_nonzero_dir_wallslide (
    double *vel, double x, double y, double z,
    ymir_locidx_t nid, void *data)
{
  collide_options_t  *collide_options = data;
  const double        wallslide_vel = collide_options->wallslide_vel;

  if (fabs (y) < SC_1000_EPS && 0.5 < z) {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = wallslide_vel * sin (2.0 * M_PI * z);
  }
  else {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }
}

//TODO
//  const double        inflow_length = 1.0/3.0;
//  double tx;
//  vel[0] = vel[2] = 0.0;
//  tx = 20.0*(z-0.75); 
//  vel[1] = 2.0 * inflow_length * lid_vel * 
//           (exp (tx) - exp (-tx)) / 
//           (exp (tx) + exp (-tx)) + 
//           inflow_length * lid_vel;
//  vel[1] = - sin (2.0 * M_PI * z);

/**
 * Sets up a linear Stokes problem.
 */
static void
collide_setup_stokes (rhea_stokes_problem_t **stokes_problem,
                      ymir_mesh_t *ymir_mesh,
                      ymir_pressure_elem_t *press_elem,
                      rhea_domain_options_t *domain_options,
                      rhea_temperature_options_t *temp_options,
                      rhea_viscosity_options_t *visc_options,
                      collide_options_t *collide_options,
                      const char *vtk_write_input_path)
{
  const char         *this_fn_name = "collide_setup_stokes";
  ymir_vec_t         *temperature, *weakzone;
  ymir_vec_t         *rhs_vel_nonzero_dirichlet = NULL;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* compute temperature */
  temperature = rhea_temperature_new (ymir_mesh);
  //TODO xi, is the following computation of temperature what you want?
  rhea_temperature_compute (temperature, temp_options);

  /* compute weak zone */
  weakzone = rhea_weakzone_new (ymir_mesh);
  //TODO xi, here you need to provide the weak zone factor
  ymir_vec_set_value (weakzone, 1.0);

  /* set velocity boundary conditions & nonzero Dirichlet values */
  switch (collide_options->vel_dir_bc) {
  case COLLIDE_VEL_DIR_BC_WALLSLIDE:
    rhea_domain_set_user_velocity_dirichlet_bc (
        collide_set_vel_dir_wallslide, NULL /* no data necessary */,
        0 /* TODO don't need this flag */);
    rhs_vel_nonzero_dirichlet = rhea_velocity_new (ymir_mesh);
    ymir_cvec_set_function (rhs_vel_nonzero_dirichlet, 
                            collide_set_rhs_vel_nonzero_dir_wallslide, 
                            collide_options);
    break;
  default: /* unknown boundary condition */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* write vtk of input data */ //TODO better move this into main fnc
  if (vtk_write_input_path != NULL) {
    const rhea_viscosity_t  viscosity_type = visc_options->type;
    ymir_vec_t         *background_temp = rhea_temperature_new (ymir_mesh);
    ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *rhs_vel = rhea_velocity_new (ymir_mesh);

    rhea_temperature_background_compute (background_temp, temp_options);
    if (viscosity_type == RHEA_VISCOSITY_NONLINEAR) {
      visc_options->type = RHEA_VISCOSITY_LINEAR;
    }
    rhea_viscosity_compute (viscosity,
                            NULL /* nl. Stokes output */,
                            NULL /* nl. Stokes output */,
                            NULL /* nl. Stokes output */,
                            temperature, weakzone,
                            NULL /* nl. Stokes input */,
                            visc_options);
    if (viscosity_type == RHEA_VISCOSITY_NONLINEAR) {
      visc_options->type = viscosity_type;
    }
    rhea_temperature_compute_rhs_vel (rhs_vel, temperature, temp_options);

    rhea_vtk_write_input_data (vtk_write_input_path, temperature,
                               background_temp, weakzone, viscosity, rhs_vel);

    rhea_temperature_destroy (background_temp);
    rhea_viscosity_destroy (viscosity);
    rhea_velocity_destroy (rhs_vel);
  }

  /* create Stokes problem */
  *stokes_problem = rhea_stokes_problem_new (
      temperature, weakzone, rhs_vel_nonzero_dirichlet, 
      ymir_mesh, press_elem,
      domain_options, temp_options, visc_options);
  rhea_stokes_problem_setup_solver (*stokes_problem);

  /* destroy */
  rhea_temperature_destroy (temperature);
  rhea_temperature_destroy (weakzone);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Cleans up Stokes problem and mesh.
 */
static void
collide_setup_clear_all (rhea_stokes_problem_t *stokes_problem,
                         p4est_t *p4est,
                         ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         rhea_discretization_options_t *discr_options)
{
  const char         *this_fn_name = "collide_setup_clear_all";
  ymir_vec_t         *rhs_vel_nonzero_dirichlet = NULL;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* destroy Stokes problem */
  rhs_vel_nonzero_dirichlet = 
    rhea_stokes_problem_get_rhs_vel_nonzero_dirichlet (stokes_problem);
  if (rhs_vel_nonzero_dirichlet != NULL) {
    rhea_velocity_destroy (rhs_vel_nonzero_dirichlet);
  }
  rhea_stokes_problem_destroy (stokes_problem);

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
collide_run_solver (ymir_vec_t *sol_vel_press,
                    ymir_mesh_t *ymir_mesh,
                    ymir_pressure_elem_t *press_elem,
                    rhea_stokes_problem_t *stokes_problem,
                    const int iter_max, const double rel_tol)
{
  const char         *this_fn_name = "collide_run_solver";
  ymir_vec_t         *rhs_vel_nonzero_dirichlet;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* run solver */
  rhea_stokes_problem_solve (sol_vel_press, iter_max, rel_tol, stokes_problem);

  /* add nonzero Dirichlet values of the velocity to the solution */
  rhs_vel_nonzero_dirichlet = 
    rhea_stokes_problem_get_rhs_vel_nonzero_dirichlet (stokes_problem);
  if (rhs_vel_nonzero_dirichlet != NULL) {
    ymir_vec_t         *sol_vel = rhea_velocity_new (ymir_mesh);

    ymir_stokes_vec_get_velocity (sol_vel_press, sol_vel, press_elem);
    ymir_vec_add (1.0, rhs_vel_nonzero_dirichlet, sol_vel);
    ymir_stokes_vec_set_velocity (sol_vel, sol_vel_press, press_elem);
    rhea_velocity_destroy (sol_vel);
  }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Runs the program.
 */
int
main (int argc, char **argv)
{
  const char         *this_fn_name = "collide:main";
  /* MPI */
  MPI_Comm            mpicomm = MPI_COMM_WORLD;
  int                 mpisize, mpirank, ompsize;
  int                 mpiret;
  /* options */
  ymir_options_t     *opt;
  rhea_domain_options_t         domain_options;
  rhea_temperature_options_t    temp_options;
  rhea_viscosity_options_t      visc_options;
  rhea_discretization_options_t discr_options;
  /* collide options */
  int                 vel_dir_bc;
  double              wallslide_vel;
  collide_options_t   collide_options;
  /* options local to this function */
  int                 production_run;
  int                 solver_iter_max;
  double              solver_rel_tol;
  char               *vtk_write_input_path;
  char               *vtk_write_solution_path;
  /* mesh */
  p4est_t            *p4est;
  ymir_mesh_t        *ymir_mesh;
  ymir_pressure_elem_t  *press_elem;
  /* Stokes */
  rhea_stokes_problem_t *stokes_problem;
  ymir_vec_t         *sol_vel_press;

  /*
   * Initialize Libraries
   */

  /* initialize rhea and sub-packages */
  rhea_initialize (argc, argv, mpicomm);

  /* get parallel environment */
  mpiret = MPI_Comm_size (mpicomm, &mpisize); YMIR_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &mpirank); YMIR_CHECK_MPI (mpiret);

#ifdef RHEA_ENABLE_OPENMP
  ompsize = omp_get_max_threads ();
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

  /* velocity Dirichlet BC's */
  YMIR_OPTIONS_I, "velocity-dirichlet-bc", '\0',
    &vel_dir_bc, COLLIDE_VEL_DIR_BC_WALLSLIDE,
    "Velocity Dirichlet boundary condition",
  YMIR_OPTIONS_D, "wallslide-velocity", '\0',
    &wallslide_vel, 1.0,
    "Tangential velocity",

  /* solver options */
  YMIR_OPTIONS_I, "solver-iter-max", '\0',
    &solver_iter_max, 100,
    "Maximum number of iterations for Stokes solver",
  YMIR_OPTIONS_D, "solver-rel-tol", '\0',
    &solver_rel_tol, 1.0e-6,
    "Relative tolerance for Stokes solver",

  /* performance & monitoring options */
  YMIR_OPTIONS_B, "production-run", '\0',
    &(production_run), 0,
    "Execute as a production run (to reduce some overhead and checks)",

  /* vtk output options */
  YMIR_OPTIONS_S, "vtk-write-input-path", '\0',
    &(vtk_write_input_path), NULL,
    "File path for vtk files for the input of the Stokes problem",
  YMIR_OPTIONS_S, "vtk-write-solution-path", '\0',
    &(vtk_write_solution_path), NULL,
    "File path for vtk files for the solution of the Stokes problem",

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
   * Process Collide Options
   */

  collide_options.vel_dir_bc = (collide_vel_dir_bc_t) vel_dir_bc;
  collide_options.wallslide_vel = wallslide_vel;

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
  rhea_process_options_all (&domain_options, &temp_options,
                            &visc_options, &discr_options);

  /*
   * Setup Mesh
   */

  collide_setup_mesh (&p4est, &ymir_mesh, &press_elem, mpicomm,
                      &domain_options, &discr_options);

  /*
   * Setup Stokes Problem
   */

  collide_setup_stokes (&stokes_problem, ymir_mesh, press_elem,
                        &domain_options, &temp_options, &visc_options,
                        &collide_options, vtk_write_input_path);

  /*
   * Solve Stokes Problem
   */

  /* initialize solution vector */
  sol_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);

  /* run solver */
  collide_run_solver (sol_vel_press, ymir_mesh, press_elem, stokes_problem,
                      solver_iter_max, solver_rel_tol);

  /* write vtk of solution */
  if (vtk_write_solution_path != NULL) {
    ymir_vec_t         *velocity = rhea_velocity_new (ymir_mesh);
    ymir_vec_t         *pressure = rhea_pressure_new (ymir_mesh, press_elem);
    ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);

    ymir_stokes_vec_get_components (sol_vel_press, velocity, pressure,
                                    press_elem);
    rhea_stokes_problem_get_viscosity (viscosity, stokes_problem);

    rhea_vtk_write_solution (vtk_write_solution_path, velocity, pressure,
                             viscosity);
  }

  /* destroy */
  rhea_velocity_pressure_destroy (sol_vel_press);

  /*
   * Finalize
   */

  /* destroy Stokes problem and mesh */
  collide_setup_clear_all (stokes_problem, p4est, ymir_mesh, press_elem,
                           &discr_options);

  /* destroy options */
  ymir_options_global_destroy ();

  /* print that this function is ending */
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

  /* finalize rhea */
  rhea_finalize ();

  return 0;
}
