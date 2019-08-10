/** SPHERE_INVERSION
 *
 * Runs an inverse problem solver for incompressible Stokes models with linear
 * and nonlinear rheologies in a spherical shell domain with two plates and one
 * subduction zone.
 */

#include <rhea.h>
#include <example_share_mesh.h>
#include <example_share_stokes.h>
#include <example_share_vtk.h>
#include <ymir_comm.h>
#include <ymir_velocity_vec.h>
#include <sc_functions.h>

/* include headers for testing purposes */
#include <rhea_viscosity_param_derivative.h>
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
  RHEA_MAIN_PERFMON_SOLVE_STOKES,
  RHEA_MAIN_PERFMON_SOLVE_INVERSION,
  RHEA_MAIN_PERFMON_TOTAL,
  RHEA_MAIN_PERFMON_N
}
rhea_main_performance_monitor_idx_t;

static const char  *rhea_main_performance_monitor_name[RHEA_MAIN_PERFMON_N] =
{
  "Setup Mesh",
  "Setup Stokes",
  "Setup Stokes Solver",
  "Solve Stokes",
  "Solve Inversion",
  "Total"
};

/******************************************************************************
 * Inverse Problem
 *****************************************************************************/

static void
sphere_inversion_test_euler_pole_computation (
                                        const int pid,
                                        rhea_stokes_problem_t *stokes_problem)
{
  rhea_plate_options_t *plate_options =
    rhea_stokes_problem_get_plate_options (stokes_problem);
  ymir_mesh_t        *ymir_mesh =
    rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  const double        rot_center[3] = {0.0, 0.0, 0.0};
  double              rot_axis_ref[3];
  double              rot_axis_eval[3];
  ymir_vec_t         *rot_surf;

  /* create work variables */
  rot_surf = rhea_velocity_surface_new (ymir_mesh);

  /* test rotation about unit vector (1,0,0) */
  rot_axis_ref[0] = 1.0;
  rot_axis_ref[1] = 0.0;
  rot_axis_ref[2] = 0.0;
  ymir_velocity_vec_generate_rotation (rot_surf, rot_center, rot_axis_ref);
  rhea_plate_velocity_evaluate_rotation (rot_axis_eval, rot_surf, pid,
                                         0 /* !project_out_mean_rot */,
                                         plate_options);
  RHEA_GLOBAL_INFOF_FN_TAG (
      __func__,
      "plate_idx=%i, rot_axis_ref=(%g, %g, %g), rot_axis_eval=(%g, %g, %g)",
      pid, rot_axis_ref[0], rot_axis_ref[1], rot_axis_ref[2],
      rot_axis_eval[0], rot_axis_eval[1], rot_axis_eval[2]);

  /* test rotation about unit vector (0,1,0) */
  rot_axis_ref[0] = 0.0;
  rot_axis_ref[1] = 1.0;
  rot_axis_ref[2] = 0.0;
  ymir_velocity_vec_generate_rotation (rot_surf, rot_center, rot_axis_ref);
  rhea_plate_velocity_evaluate_rotation (rot_axis_eval, rot_surf, pid,
                                         0 /* !project_out_mean_rot */,
                                         plate_options);
  RHEA_GLOBAL_INFOF_FN_TAG (
      __func__,
      "plate_idx=%i, rot_axis_ref=(%g, %g, %g), rot_axis_eval=(%g, %g, %g)",
      pid, rot_axis_ref[0], rot_axis_ref[1], rot_axis_ref[2],
      rot_axis_eval[0], rot_axis_eval[1], rot_axis_eval[2]);

  /* test rotation about unit vector (0,0,1) */
  rot_axis_ref[0] = 0.0;
  rot_axis_ref[1] = 0.0;
  rot_axis_ref[2] = 1.0;
  ymir_velocity_vec_generate_rotation (rot_surf, rot_center, rot_axis_ref);
  rhea_plate_velocity_evaluate_rotation (rot_axis_eval, rot_surf, pid,
                                         0 /* !project_out_mean_rot */,
                                         plate_options);
  RHEA_GLOBAL_INFOF_FN_TAG (
      __func__,
      "plate_idx=%i, rot_axis_ref=(%g, %g, %g), rot_axis_eval=(%g, %g, %g)",
      pid, rot_axis_ref[0], rot_axis_ref[1], rot_axis_ref[2],
      rot_axis_eval[0], rot_axis_eval[1], rot_axis_eval[2]);

  /* test rotation about random vector */
  rot_axis_ref[0] = sc_rand_uniform ();
  rot_axis_ref[1] = sc_rand_uniform ();
  rot_axis_ref[2] = sc_rand_uniform ();
  ymir_velocity_vec_generate_rotation (rot_surf, rot_center, rot_axis_ref);
  rhea_plate_velocity_evaluate_rotation (rot_axis_eval, rot_surf, pid,
                                         0 /* !project_out_mean_rot */,
                                         plate_options);
  RHEA_GLOBAL_INFOF_FN_TAG (
      __func__,
      "plate_idx=%i, rot_axis_ref=(%g, %g, %g), rot_axis_eval=(%g, %g, %g)",
      pid, rot_axis_ref[0], rot_axis_ref[1], rot_axis_ref[2],
      rot_axis_eval[0], rot_axis_eval[1], rot_axis_eval[2]);

  /* destroy */
  ymir_vec_destroy (rot_surf);
}

static void
sphere_inversion_solve_with_vel_obs (rhea_inversion_problem_t *inv_problem,
                                     ymir_vec_t *sol_vel_press,
                                     const int vel_obs_euler_pole,
                                     const double vel_obs_add_noise_stddev,
                                     rhea_stokes_problem_t *stokes_problem)
{
  rhea_plate_options_t *plate_options =
    rhea_stokes_problem_get_plate_options (stokes_problem);
  const int             n_plates = rhea_plate_get_n_plates (plate_options);
  ymir_mesh_t          *ymir_mesh;
  ymir_pressure_elem_t *press_elem;
  ymir_vec_t         *vel_sol;
  ymir_vec_t         *vel_obs_surf;

  /* get mesh data */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);

  /* generate synthetic velocity observations from previous Stokes solve */
  vel_sol = rhea_velocity_new (ymir_mesh);
  rhea_velocity_pressure_copy_components (vel_sol, NULL, sol_vel_press,
                                          press_elem);
  ymir_vec_share_owned (vel_sol);
  rhea_stokes_problem_velocity_boundary_set_zero (vel_sol, stokes_problem);
  RHEA_ASSERT (rhea_velocity_is_valid (vel_sol));

  /* project velocity from volume to surface */
  vel_obs_surf = rhea_velocity_surface_new (ymir_mesh);
  rhea_velocity_surface_interpolate (vel_obs_surf, vel_sol);

  /* run inversion */
  if (0 < n_plates && vel_obs_euler_pole) { /* set up Euler poles */
    double              rot_axis[3];
    int                 pid;

    /* calculate rotational axis */
    RHEA_ASSERT (plate_options != NULL);
    for (pid = 0; pid < n_plates; pid++) { /* loop over all plates */
#ifdef RHEA_ENABLE_DEBUG
      sphere_inversion_test_euler_pole_computation (pid, stokes_problem);
#endif
      rhea_plate_velocity_evaluate_rotation (rot_axis, vel_obs_surf, pid,
                                             0 /* !project_out_mean_rot */,
                                             plate_options);
      plate_options->angular_velocity[3*pid    ] = rot_axis[0];
      plate_options->angular_velocity[3*pid + 1] = rot_axis[1];
      plate_options->angular_velocity[3*pid + 2] = rot_axis[2];
      RHEA_GLOBAL_INFOF_FN_TAG (
          __func__, "plate_idx=%i, vel_obs_rot_axis=(%g, %g, %g)",
          pid, rot_axis[0], rot_axis[1], rot_axis[2]);
    }

    /* run solver */
    rhea_inversion_solve (inv_problem,
                          1 /* use model param's as initial guess */,
                          NULL /* no parameter vector */);
  }
  else { /* use spacially varying surface velocity */
    ymir_vec_t         *vel_obs_weight_surf;

    /* generate weight function */
    vel_obs_weight_surf = ymir_face_cvec_new (
        ymir_mesh, RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    ymir_vec_set_value (vel_obs_weight_surf, 1.0);
    if (0 < n_plates) { /* if filter subducting plate */
      rhea_plate_apply_filter_vec (vel_obs_weight_surf, 0 /* plate_label */,
                                   plate_options);
    }

    /* run solver */
    rhea_inversion_solve_with_vel_obs (
        inv_problem, 1 /* use model param's as initial guess */,
        NULL /* no parameter vector */, vel_obs_surf, vel_obs_weight_surf,
        vel_obs_add_noise_stddev);

    /* destroy */
    ymir_vec_destroy (vel_obs_weight_surf);
  }

  /* destroy */
  rhea_velocity_destroy (vel_sol);
  rhea_velocity_surface_destroy (vel_obs_surf);
}

static void
sphere_inversion_vtk_write_param_derivatives (
                                        const char *vtk_path,
                                        ymir_vec_t *sol_vel_press,
                                        rhea_stokes_problem_t *stokes_problem)
{
  rhea_weakzone_options_t  *weak_options =
    rhea_stokes_problem_get_weakzone_options (stokes_problem);
  rhea_viscosity_options_t *visc_options =
    rhea_stokes_problem_get_viscosity_options (stokes_problem);
  ymir_mesh_t          *ymir_mesh =
                          rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  ymir_pressure_elem_t *press_elem =
                          rhea_stokes_problem_get_press_elem (stokes_problem);
  ymir_vec_t         *temperature =
                        rhea_stokes_problem_get_temperature (stokes_problem);
  ymir_vec_t         *weakzone =
                        rhea_stokes_problem_get_weakzone (stokes_problem);
  ymir_vec_t         *vel_sol;
  ymir_vec_t         *viscosity, *bounds_marker, *yielding_marker;
  char                path[BUFSIZ];

  /* retrieve velocity */
  vel_sol = rhea_velocity_new (ymir_mesh);
  rhea_velocity_pressure_copy_components (vel_sol, NULL, sol_vel_press,
                                          press_elem);
  RHEA_ASSERT (rhea_velocity_is_valid (vel_sol));

  /* compute viscosity and related fields */
  viscosity = rhea_viscosity_new (ymir_mesh);
  bounds_marker = rhea_viscosity_new (ymir_mesh);
  yielding_marker = rhea_viscosity_new (ymir_mesh);
  rhea_viscosity_compute (
      /* out: */ viscosity, NULL, bounds_marker, yielding_marker,
      /* in:  */ temperature, weakzone, vel_sol, visc_options);

  /*
   * Viscosity Parameters
   */
  {
    ymir_vec_t         *deriv_min = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *deriv_max = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *deriv_um_scal = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *deriv_um_Ea = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *deriv_lm_scal = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *deriv_lm_Ea = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *deriv_stress_exp = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *deriv_yielding = rhea_viscosity_new (ymir_mesh);

    rhea_viscosity_param_derivative (
        deriv_min,        RHEA_VISCOSITY_PARAM_DERIVATIVE_MIN,
        viscosity, bounds_marker, yielding_marker, temperature, vel_sol,
        visc_options);
    rhea_viscosity_param_derivative (
        deriv_max,        RHEA_VISCOSITY_PARAM_DERIVATIVE_MAX,
        viscosity, bounds_marker, yielding_marker, temperature, vel_sol,
        visc_options);
    rhea_viscosity_param_derivative (
        deriv_um_scal,    RHEA_VISCOSITY_PARAM_DERIVATIVE_UPPER_MANTLE_SCALING,
        viscosity, bounds_marker, yielding_marker, temperature, vel_sol,
        visc_options);
    rhea_viscosity_param_derivative (
        deriv_um_Ea,
        RHEA_VISCOSITY_PARAM_DERIVATIVE_UPPER_MANTLE_ACTIVATION_ENERGY,
        viscosity, bounds_marker, yielding_marker, temperature, vel_sol,
        visc_options);
    rhea_viscosity_param_derivative (
        deriv_lm_scal,    RHEA_VISCOSITY_PARAM_DERIVATIVE_LOWER_MANTLE_SCALING,
        viscosity, bounds_marker, yielding_marker, temperature, vel_sol,
        visc_options);
    rhea_viscosity_param_derivative (
        deriv_lm_Ea,
        RHEA_VISCOSITY_PARAM_DERIVATIVE_LOWER_MANTLE_ACTIVATION_ENERGY,
        viscosity, bounds_marker, yielding_marker, temperature, vel_sol,
        visc_options);
    rhea_viscosity_param_derivative (
        deriv_stress_exp, RHEA_VISCOSITY_PARAM_DERIVATIVE_STRESS_EXPONENT,
        viscosity, bounds_marker, yielding_marker, temperature, vel_sol,
        visc_options);
    rhea_viscosity_param_derivative (
        deriv_yielding,   RHEA_VISCOSITY_PARAM_DERIVATIVE_YIELD_STRENGTH,
        viscosity, bounds_marker, yielding_marker, temperature, vel_sol,
        visc_options);

    snprintf (path, BUFSIZ, "%s_visc_params", vtk_path);
    ymir_vtk_write (ymir_mesh, path,
                    deriv_min,        "derivative_min",
                    deriv_max,        "derivative_max",
                    deriv_um_scal,    "derivative_um_scal",
                    deriv_um_Ea,      "derivative_um_Ea",
                    deriv_lm_scal,    "derivative_lm_scal",
                    deriv_lm_Ea,      "derivative_lm_Ea",
                    deriv_stress_exp, "derivative_stress_exp",
                    deriv_yielding,   "derivative_yielding",
                    temperature, "temperature", viscosity, "viscosity", NULL);

    rhea_viscosity_destroy (deriv_min);
    rhea_viscosity_destroy (deriv_max);
    rhea_viscosity_destroy (deriv_um_scal);
    rhea_viscosity_destroy (deriv_um_Ea);
    rhea_viscosity_destroy (deriv_lm_scal);
    rhea_viscosity_destroy (deriv_lm_Ea);
    rhea_viscosity_destroy (deriv_stress_exp);
    rhea_viscosity_destroy (deriv_yielding);
  }

  /*
   * Weak Zone Parameters
   */
  {
    const rhea_viscosity_param_derivative_t deriv_type =
      RHEA_VISCOSITY_PARAM_DERIVATIVE_WEAK_FACTOR_INTERIOR;
  //ymir_vec_t         *deriv_weak_thick = rhea_viscosity_new (ymir_mesh);
  //ymir_vec_t         *deriv_weak_thick_const = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *deriv_weak_factor_sl = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *deriv_weak_factor_ri = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *deriv_weak_factor_fz = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *deriv_weak_factor_1001 = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *deriv_weak_factor_2001 = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *deriv_weak_factor_3001 = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *deriv_weak_factor_3002 = rhea_viscosity_new (ymir_mesh);

    rhea_viscosity_param_derivative_weakzone (
        deriv_weak_factor_sl, deriv_type, RHEA_WEAKZONE_LABEL_CLASS_SLAB,
        viscosity, bounds_marker, yielding_marker, weakzone,
        weak_options, visc_options);
    rhea_viscosity_param_derivative_weakzone (
        deriv_weak_factor_ri, deriv_type, RHEA_WEAKZONE_LABEL_CLASS_RIDGE,
        viscosity, bounds_marker, yielding_marker, weakzone,
        weak_options, visc_options);
    rhea_viscosity_param_derivative_weakzone (
        deriv_weak_factor_fz, deriv_type, RHEA_WEAKZONE_LABEL_CLASS_FRACTURE,
        viscosity, bounds_marker, yielding_marker, weakzone,
        weak_options, visc_options);
    rhea_viscosity_param_derivative_weakzone (
        deriv_weak_factor_1001, deriv_type, (rhea_weakzone_label_t) 1001,
        viscosity, bounds_marker, yielding_marker, weakzone,
        weak_options, visc_options);
    rhea_viscosity_param_derivative_weakzone (
        deriv_weak_factor_2001, deriv_type, (rhea_weakzone_label_t) 2001,
        viscosity, bounds_marker, yielding_marker, weakzone,
        weak_options, visc_options);
    rhea_viscosity_param_derivative_weakzone (
        deriv_weak_factor_3001, deriv_type, (rhea_weakzone_label_t) 3001,
        viscosity, bounds_marker, yielding_marker, weakzone,
        weak_options, visc_options);
    rhea_viscosity_param_derivative_weakzone (
        deriv_weak_factor_3002, deriv_type, (rhea_weakzone_label_t) 3002,
        viscosity, bounds_marker, yielding_marker, weakzone,
        weak_options, visc_options);

    snprintf (path, BUFSIZ, "%s_weak_params", vtk_path);
    ymir_vtk_write (ymir_mesh, path,
                    deriv_weak_factor_sl,   "derivative_weak_factor_sl",
                    deriv_weak_factor_ri,   "derivative_weak_factor_ri",
                    deriv_weak_factor_fz,   "derivative_weak_factor_fz",
                    deriv_weak_factor_1001, "derivative_weak_factor_1001",
                    deriv_weak_factor_2001, "derivative_weak_factor_2001",
                    deriv_weak_factor_3001, "derivative_weak_factor_3001",
                    deriv_weak_factor_3002, "derivative_weak_factor_3002",
                    weakzone, "weakzone", viscosity, "viscosity", NULL);

    rhea_viscosity_destroy (deriv_weak_factor_sl);
    rhea_viscosity_destroy (deriv_weak_factor_ri);
    rhea_viscosity_destroy (deriv_weak_factor_fz);
    rhea_viscosity_destroy (deriv_weak_factor_1001);
    rhea_viscosity_destroy (deriv_weak_factor_2001);
    rhea_viscosity_destroy (deriv_weak_factor_3001);
    rhea_viscosity_destroy (deriv_weak_factor_3002);
  }

  /* destroy */
  rhea_velocity_destroy (vel_sol);
  rhea_viscosity_destroy (viscosity);
  rhea_viscosity_destroy (bounds_marker);
  rhea_viscosity_destroy (yielding_marker);
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
  static const char   func_name[] = "sphere_inversion:main";
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
  int                 vel_obs_euler_pole;
  double              vel_obs_add_noise_stddev;
  char               *bin_solver_path;
  char               *txt_inv_solver_path;
  char               *vtk_input_path;
  char               *vtk_solution_path;
  char               *vtk_solver_path;
  char               *vtk_inv_solver_path;
  char               *vtk_inv_param_derivative_path;
  /* mesh */
  p4est_t                *p4est;
  ymir_mesh_t            *ymir_mesh;
  ymir_pressure_elem_t   *press_elem;
  /* Stokes */
  rhea_stokes_problem_t  *stokes_problem;
  ymir_vec_t         *sol_vel_press;
  int                 nonzero_inital_guess;
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

  /* inversion options */
  YMIR_OPTIONS_B, "velocity-observations-euler-pole", '\0',
    &(vel_obs_euler_pole), 0,
    "Velocity observations are generated from rotation of manufactured solution",
  YMIR_OPTIONS_D, "velocity-observations-add-noise-stddev", '\0',
    &(vel_obs_add_noise_stddev), NAN,
    "Standard deviation of noise added to manufactured velocity observations",

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
  YMIR_OPTIONS_S, "vtk-write-inversion-solver-path", '\0',
    &(vtk_inv_solver_path), NULL,
    "VTK file path for solver internals of the inversion.",
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

  /*
   * Solve Stokes Problem
   */

  /* setup Stokes solver */
  rhea_performance_monitor_start_barrier (RHEA_MAIN_PERFMON_SETUP_SOLVER);
  rhea_stokes_problem_setup_solver_ext (stokes_problem,
                                        RHEA_INVERSION_KRYLOV_SOLVER_N);
  rhea_performance_monitor_stop_add (RHEA_MAIN_PERFMON_SETUP_SOLVER);

  /* initialize solution vector */
  sol_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);
  nonzero_inital_guess = 0;

  /* run Stokes solver */
  rhea_performance_monitor_start_barrier (RHEA_MAIN_PERFMON_SOLVE_STOKES);
  rhea_stokes_problem_solve (&sol_vel_press, nonzero_inital_guess,
                             solver_iter_max, solver_rel_tol, stokes_problem);
  rhea_performance_monitor_stop_add (RHEA_MAIN_PERFMON_SOLVE_STOKES);

  /* write vtk of solution */
  example_share_vtk_write_solution (vtk_solution_path, sol_vel_press,
                                    stokes_problem);

  /*
   * Test Derivatives of the Viscosity w.r.t. the Inversion Parameters
   */

  if (NULL != vtk_inv_param_derivative_path) {
    sphere_inversion_vtk_write_param_derivatives (
        vtk_inv_param_derivative_path, sol_vel_press, stokes_problem);
  }

  /*
   * Solve Inverse Problem
   */

  /* setup inversion solver */
  inv_problem = rhea_inversion_new (stokes_problem);
  rhea_inversion_set_txt_output (inv_problem, txt_inv_solver_path);
  rhea_inversion_set_vtk_output (inv_problem, vtk_inv_solver_path);

  /* run inversion solver */
  rhea_performance_monitor_start_barrier (RHEA_MAIN_PERFMON_SOLVE_INVERSION);
  sphere_inversion_solve_with_vel_obs (inv_problem, sol_vel_press,
                                       vel_obs_euler_pole,
                                       vel_obs_add_noise_stddev,
                                       stokes_problem);
  rhea_performance_monitor_stop_add_barrier (RHEA_MAIN_PERFMON_SOLVE_INVERSION);

  /* destroy */
  rhea_inversion_destroy (inv_problem);

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
