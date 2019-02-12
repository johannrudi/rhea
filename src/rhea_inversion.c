#include <rhea_inversion.h>
#include <rhea_base.h>
#include <rhea_inversion_param.h>
#include <rhea_inversion_obs_velocity.h>
#include <rhea_newton.h>
#include <rhea_velocity.h>
#include <rhea_velocity_pressure.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_INVERSION_DEFAULT_VEL_OBS_TYPE (RHEA_INVERSION_OBS_VELOCITY_NONE)
#define RHEA_INVERSION_DEFAULT_FORWARD_SOLVER_ITER_MAX (100)
#define RHEA_INVERSION_DEFAULT_FORWARD_SOLVER_RTOL (1.0e-6)
#define RHEA_INVERSION_DEFAULT_ADJOINT_SOLVER_ITER_MAX (100)
#define RHEA_INVERSION_DEFAULT_ADJOINT_SOLVER_RTOL (1.0e-6)
#define RHEA_INVERSION_DEFAULT_CHECK_GRADIENT (0)
#define RHEA_INVERSION_DEFAULT_CHECK_HESSIAN (0)

/* initialize options */
int                 rhea_inversion_vel_obs_type =
                      RHEA_INVERSION_DEFAULT_VEL_OBS_TYPE;
int                 rhea_inversion_forward_solver_iter_max =
                      RHEA_INVERSION_DEFAULT_FORWARD_SOLVER_ITER_MAX;
double              rhea_inversion_forward_solver_rtol =
                      RHEA_INVERSION_DEFAULT_FORWARD_SOLVER_RTOL;
int                 rhea_inversion_adjoint_solver_iter_max =
                      RHEA_INVERSION_DEFAULT_ADJOINT_SOLVER_ITER_MAX;
double              rhea_inversion_adjoint_solver_rtol =
                      RHEA_INVERSION_DEFAULT_ADJOINT_SOLVER_RTOL;
int                 rhea_inversion_check_gradient =
                      RHEA_INVERSION_DEFAULT_CHECK_GRADIENT;
int                 rhea_inversion_check_hessian =
                      RHEA_INVERSION_DEFAULT_CHECK_HESSIAN;

rhea_inversion_param_options_t  rhea_inversion_param_options;
rhea_newton_options_t           rhea_inversion_newton_options;

void
rhea_inversion_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Inversion";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  /* observational data options */
  YMIR_OPTIONS_I, "velocity-observations-type", '\0',
    &(rhea_inversion_vel_obs_type), RHEA_INVERSION_DEFAULT_VEL_OBS_TYPE,
    "Type of velocity observations at the surface",

  /* inner solver (i.e., forward/adjoint) options */
  YMIR_OPTIONS_I, "forward-solver-iter-max", '\0',
    &(rhea_inversion_forward_solver_iter_max),
    RHEA_INVERSION_DEFAULT_FORWARD_SOLVER_ITER_MAX,
    "Forward solver: Maximum number of iterations for Stokes solver",
  YMIR_OPTIONS_D, "forward-solver-rtol", '\0',
    &(rhea_inversion_forward_solver_rtol),
    RHEA_INVERSION_DEFAULT_FORWARD_SOLVER_RTOL,
    "Forward solver: Relative tolerance for Stokes solver",
  YMIR_OPTIONS_I, "adjoint-solver-iter-max", '\0',
    &(rhea_inversion_adjoint_solver_iter_max),
    RHEA_INVERSION_DEFAULT_ADJOINT_SOLVER_ITER_MAX,
    "Adjoint solver: Maximum number of iterations for Stokes solver",
  YMIR_OPTIONS_D, "adjoint-solver-rtol", '\0',
    &(rhea_inversion_adjoint_solver_rtol),
    RHEA_INVERSION_DEFAULT_ADJOINT_SOLVER_RTOL,
    "Adjoint solver: Relative tolerance for Stokes solver",

  /* outer solver (inversion) options */
  YMIR_OPTIONS_B, "check-gradient", '\0',
    &(rhea_inversion_check_gradient), RHEA_INVERSION_DEFAULT_CHECK_GRADIENT,
    "Check the gradient during Newton iterations",
  YMIR_OPTIONS_B, "check-hessian", '\0',
    &(rhea_inversion_check_hessian), RHEA_INVERSION_DEFAULT_CHECK_HESSIAN,
    "Check the Hessian during Newton iterations",


  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add sub-options */
  rhea_inversion_param_add_options (&rhea_inversion_param_options, opt);
  rhea_newton_add_options (&rhea_inversion_newton_options, opt);

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);
}

/******************************************************************************
 * Inverse Problem Structure
 *****************************************************************************/

/* inverse problem */
struct rhea_inversion_problem
{
  /* Stokes problem (not owned) */
  rhea_stokes_problem_t  *stokes_problem;

  /* inversion parameters */
  rhea_inversion_param_options_t *inv_param_options;
  rhea_inversion_param_t         *inv_param;

  /* observations: velocity at the surface */
  rhea_inversion_obs_velocity_t vel_obs_type;
  ymir_vec_t         *vel_obs_surf;
  ymir_vec_t         *vel_obs_weight_surf;

  /* forward and adjoint states */
  ymir_vec_t         *forward_vel_press;
  ymir_vec_t         *forward_vel;
  ymir_vec_t         *adjoint_vel_press;

  /* Newton problem */
  rhea_newton_options_t  *newton_options;
  rhea_newton_problem_t  *newton_problem;
  ymir_vec_t             *newton_neg_gradient_vec;
  ymir_vec_t             *newton_step_vec;
};

static int
rhea_inversion_solver_data_exists (rhea_inversion_problem_t *inv_problem)
{
  int                 fwd_adj_states_exist;
  int                 vel_obs_exist;

  fwd_adj_states_exist = (inv_problem->forward_vel_press != NULL &&
                          inv_problem->adjoint_vel_press != NULL &&
                          inv_problem->forward_vel != NULL);
  vel_obs_exist = (inv_problem->vel_obs_surf != NULL);

  return (fwd_adj_states_exist && vel_obs_exist);
}

/******************************************************************************
 * Callback Functions for Newton Solver
 *****************************************************************************/

static void
rhea_inversion_newton_create_solver_data_fn (ymir_vec_t *solution, void *data)
{
  rhea_inversion_problem_t *inv_problem = data;
  rhea_stokes_problem_t    *stokes_problem = inv_problem->stokes_problem;
  ymir_mesh_t              *ymir_mesh;
  ymir_pressure_elem_t     *press_elem;
  ymir_vec_t         *parameter_vec;
  const int           stokes_nonzero_inital_guess = 0;
  const int           stokes_iter_max = rhea_inversion_forward_solver_iter_max;
  const double        stokes_rtol = rhea_inversion_forward_solver_rtol;

  /* exit if nothing to do */
  if (rhea_inversion_solver_data_exists (inv_problem)) {
    return;
  }

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (inv_problem->forward_vel_press == NULL);
  RHEA_ASSERT (inv_problem->adjoint_vel_press == NULL);
  RHEA_ASSERT (inv_problem->forward_vel == NULL);
  RHEA_ASSERT (inv_problem->vel_obs_surf == NULL);
  RHEA_ASSERT (inv_problem->vel_obs_weight_surf == NULL);

  /* get mesh data */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);

  /* create observational data */
  inv_problem->vel_obs_surf = rhea_velocity_surface_new (ymir_mesh);
  inv_problem->vel_obs_weight_surf = ymir_face_cvec_new (
      ymir_mesh, RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);

  //TODO fill vel obs type, obs vec, and weight
  ymir_vec_set_value (inv_problem->vel_obs_surf, NAN);
  ymir_vec_set_value (inv_problem->vel_obs_weight_surf, 1.0);

  /* create forward and adjoint states */
  inv_problem->forward_vel_press = rhea_velocity_pressure_new (ymir_mesh,
                                                               press_elem);
  inv_problem->adjoint_vel_press = rhea_velocity_pressure_new (ymir_mesh,
                                                               press_elem);
  inv_problem->forward_vel = rhea_velocity_new (ymir_mesh);

  /* update inversion parameters */
  parameter_vec = rhea_inversion_param_get_vector (inv_problem->inv_param);
  if (solution != NULL && solution != parameter_vec) { /* if copy params */
    ymir_vec_copy (solution, parameter_vec);
  }
  rhea_inversion_param_push_to_model (inv_problem->inv_param);

  /* set up Stokes solver for forward problem */
  rhea_stokes_problem_setup_solver (stokes_problem);

  /* solve for forward state */
  rhea_stokes_problem_solve (
      &(inv_problem->forward_vel_press), stokes_nonzero_inital_guess,
      stokes_iter_max, stokes_rtol, inv_problem->stokes_problem);

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static void
rhea_inversion_newton_clear_solver_data_fn (void *data)
{
  rhea_inversion_problem_t *inv_problem = data;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (inv_problem->forward_vel_press != NULL);
  RHEA_ASSERT (inv_problem->adjoint_vel_press != NULL);
  RHEA_ASSERT (inv_problem->forward_vel != NULL);

  /* destroy forward and adjoint states */
  rhea_velocity_pressure_destroy (inv_problem->forward_vel_press);
  rhea_velocity_pressure_destroy (inv_problem->adjoint_vel_press);
  rhea_velocity_destroy (inv_problem->forward_vel);
  inv_problem->forward_vel_press = NULL;
  inv_problem->adjoint_vel_press = NULL;
  inv_problem->forward_vel = NULL;

  /* destroy observational data */
  if (inv_problem->vel_obs_surf != NULL) {
    rhea_velocity_surface_destroy (inv_problem->vel_obs_surf);
    inv_problem->vel_obs_surf = NULL;
  }
  if (inv_problem->vel_obs_weight_surf != NULL) {
    ymir_vec_destroy (inv_problem->vel_obs_weight_surf);
    inv_problem->vel_obs_weight_surf = NULL;
  }

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static void
rhea_inversion_newton_update_operator_fn (ymir_vec_t *solution, void *data)
{
  rhea_inversion_problem_t *inv_problem = data;
  ymir_vec_t         *parameter_vec;
  const int           stokes_nonzero_inital_guess = 1;
  const int           stokes_iter_max = rhea_inversion_forward_solver_iter_max;
  const double        stokes_rtol = rhea_inversion_forward_solver_rtol;

  //TODO add performance monitoring

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* update inversion parameters */
  parameter_vec = rhea_inversion_param_get_vector (inv_problem->inv_param);
  if (solution != parameter_vec) { /* if copy params */
    ymir_vec_copy (solution, parameter_vec);
  }
  rhea_inversion_param_push_to_model (inv_problem->inv_param);

  /* update Stokes solver for forward problem */
  rhea_stokes_problem_update_solver (
      inv_problem->stokes_problem, inv_problem->forward_vel_press,
      0 /* !override_rhs */, NULL, NULL, NULL, NULL);

  /* solve for forward state */
  rhea_stokes_problem_solve (
      &(inv_problem->forward_vel_press), stokes_nonzero_inital_guess,
      stokes_iter_max, stokes_rtol, inv_problem->stokes_problem);

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static void
rhea_inversion_newton_update_hessian_fn (ymir_vec_t *solution,
                                         ymir_vec_t *step_vec,
                                         const double step_length,
                                         void *data)
{
  //rhea_inversion_problem_t *inv_problem = data;

  //TODO add performance monitoring

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  //TODO

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static double
rhea_inversion_newton_evaluate_objective_fn (ymir_vec_t *solution, void *data)
{
  rhea_inversion_problem_t *inv_problem = data;
  rhea_stokes_problem_t    *stokes_problem = inv_problem->stokes_problem;
  ymir_pressure_elem_t     *press_elem;
  rhea_domain_options_t    *domain_options;
  double              obj_val;

  //TODO add performance monitoring

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* get mesh data */
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);
  domain_options = rhea_stokes_problem_get_domain_options (stokes_problem);

  /* retrieve velocity of the forward state */
  rhea_velocity_pressure_copy_components (
      inv_problem->forward_vel, NULL, inv_problem->forward_vel_press,
      press_elem);
  rhea_stokes_problem_velocity_boundary_set_zero (
      inv_problem->forward_vel, inv_problem->stokes_problem);
  RHEA_ASSERT (rhea_velocity_is_valid (inv_problem->forward_vel));

  /* compute misfit term, e.g., (squared norm of the) observation data misfit */
  obj_val = rhea_inversion_obs_velocity_misfit (
      inv_problem->forward_vel, inv_problem->vel_obs_surf,
      inv_problem->vel_obs_weight_surf, inv_problem->vel_obs_type,
      domain_options);
  obj_val *= 0.5;

  /* compute & add prior term */
  //TODO

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);

  /* return value of objective functional */
  return obj_val;
}

static void
rhea_inversion_newton_compute_negative_gradient_fn (
                                            ymir_vec_t *neg_gradient,
                                            ymir_vec_t *solution, void *data)
{
  rhea_inversion_problem_t *inv_problem = data;
  rhea_stokes_problem_t    *stokes_problem = inv_problem->stokes_problem;
  ymir_mesh_t              *ymir_mesh;
  ymir_pressure_elem_t     *press_elem;
  ymir_vec_t         *adjoint_vel_press = inv_problem->adjoint_vel_press;
  ymir_vec_t         *rhs_vel_mass;
  ymir_vec_t         *rhs_vel_press;
  const int           stokes_nonzero_inital_guess = 0;
  const int           stokes_iter_max = rhea_inversion_adjoint_solver_iter_max;
  const double        stokes_rtol = rhea_inversion_adjoint_solver_rtol;

  //TODO add performance monitoring

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* Note: Assume that the forward state `inv_problem->forward_vel_press` was
   *       computed before, during the operator update, and that the velocity
   *       `inv_problem->forward_vel` was extracted from it during the
   *       evaluation of the objective functional. */

  /* get mesh data */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);

  /* set Stokes right-hand side for adjoint problem */
  rhs_vel_mass = rhea_velocity_new (ymir_mesh);
  rhs_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);
  rhea_inversion_obs_velocity_adjoint_rhs (
      rhs_vel_mass, inv_problem->forward_vel, inv_problem->vel_obs_surf,
      inv_problem->vel_obs_weight_surf, inv_problem->vel_obs_type,
      rhea_stokes_problem_get_domain_options (stokes_problem));
  ymir_vec_set_zero (rhs_vel_press);
  rhea_velocity_pressure_set_components (rhs_vel_press, rhs_vel_mass, NULL,
                                         press_elem);
  rhea_velocity_destroy (rhs_vel_mass);

  /* update Stokes solver for adjoint problem */
  rhea_stokes_problem_update_solver (
      stokes_problem, inv_problem->forward_vel_press,
      1 /* override_rhs */, rhs_vel_press, NULL, NULL, NULL);
  rhea_velocity_pressure_destroy (rhs_vel_press);

  /* solve for adjoint state */
  rhea_stokes_problem_solve_ext (
      &adjoint_vel_press, stokes_nonzero_inital_guess, stokes_iter_max,
      stokes_rtol, stokes_problem, 1 /* force linear solve */);
  inv_problem->adjoint_vel_press = adjoint_vel_press;

  /* compute negative gradient */
  rhea_inversion_param_compute_gradient (
      inv_problem->newton_neg_gradient_vec,
      inv_problem->forward_vel_press,
      inv_problem->adjoint_vel_press,
      inv_problem->inv_param);
  ymir_vec_scale (-1.0, inv_problem->newton_neg_gradient_vec);

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static double
rhea_inversion_newton_compute_gradient_norm_fn (ymir_vec_t *neg_gradient,
                                                void *data, double *norm_comp)
{
  rhea_inversion_problem_t *inv_problem = data;

  return rhea_inversion_param_compute_gradient_norm (
      inv_problem->newton_neg_gradient_vec, inv_problem->inv_param);
}

static void
rhea_inversion_newton_apply_hessian_fn (ymir_vec_t *out, ymir_vec_t *in,
                                        void *data)
{
  //rhea_inversion_problem_t *inv_problem = data;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  //TODO

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static int
rhea_inversion_newton_solve_hessian_system_fn (
                                            ymir_vec_t *step,
                                            ymir_vec_t *neg_gradient,
                                            const int lin_iter_max,
                                            const double lin_res_norm_rtol,
                                            const int nonzero_initial_guess,
                                            void *data, int *lin_iter_count)
{
  //rhea_inversion_problem_t *inv_problem = data;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (
      __func__, "lin_iter_max=%i, lin_rtol=%.1e, nonzero_init_guess=%i",
      lin_iter_max, lin_res_norm_rtol, nonzero_initial_guess);

  //TODO

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);

  return 0;
}

static void
rhea_inversion_newton_output_prestep_fn (ymir_vec_t *solution, const int iter,
                                         void *data)
{
  //rhea_inversion_problem_t *inv_problem = data;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  //TODO

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

/******************************************************************************
 * Inverse Problem
 *****************************************************************************/

rhea_inversion_problem_t *
rhea_inversion_new (rhea_stokes_problem_t *stokes_problem)
{
  rhea_inversion_problem_t *inv_problem;
  ymir_vec_t         *parameter_vec;

  RHEA_GLOBAL_PRODUCTION_FN_BEGIN (__func__);

  /* check input */
  //TODO

  /* initialize inverse problem */
  inv_problem = RHEA_ALLOC (rhea_inversion_problem_t, 1);
  inv_problem->stokes_problem = stokes_problem;
  inv_problem->vel_obs_type = rhea_inversion_vel_obs_type;
  inv_problem->vel_obs_surf = NULL;
  inv_problem->vel_obs_weight_surf = NULL;
  inv_problem->forward_vel_press = NULL;
  inv_problem->forward_vel = NULL;
  inv_problem->adjoint_vel_press = NULL;

  /* create parameters */
  inv_problem->inv_param_options = &rhea_inversion_param_options;
  inv_problem->inv_param = rhea_inversion_param_new (
      stokes_problem, inv_problem->inv_param_options);
  parameter_vec = rhea_inversion_param_get_vector (inv_problem->inv_param);

  /* initialize Newton problem */
  inv_problem->newton_options = &rhea_inversion_newton_options;
  inv_problem->newton_neg_gradient_vec = ymir_vec_template (parameter_vec);
  inv_problem->newton_step_vec = ymir_vec_template (parameter_vec);

  /* create Newton problem */
  {
    const int           grad_norm_n_components = 2;//TODO
    rhea_newton_problem_t *newton_problem;

    newton_problem = rhea_newton_problem_new (
        rhea_inversion_newton_compute_negative_gradient_fn,
        rhea_inversion_newton_compute_gradient_norm_fn,
        grad_norm_n_components,
        rhea_inversion_newton_solve_hessian_system_fn);
    inv_problem->newton_problem = newton_problem;

    rhea_newton_problem_set_vectors (
        inv_problem->newton_neg_gradient_vec,
        inv_problem->newton_step_vec, newton_problem);

    rhea_newton_problem_set_data_fn (
        inv_problem,
        rhea_inversion_newton_create_solver_data_fn,
        rhea_inversion_newton_clear_solver_data_fn, newton_problem);
    rhea_newton_problem_set_evaluate_objective_fn (
        rhea_inversion_newton_evaluate_objective_fn, newton_problem);
    rhea_newton_problem_set_apply_hessian_fn (
        rhea_inversion_newton_apply_hessian_fn, newton_problem);
    rhea_newton_problem_set_update_fn (
        rhea_inversion_newton_update_operator_fn,
        rhea_inversion_newton_update_hessian_fn,
        NULL /*rhea_inversion_newton_modify_hessian_system_fn*/, newton_problem);
  //rhea_newton_problem_set_setup_poststep_fn (
  //    rhea_inversion_newton_amr_fn, newton_problem);
    rhea_newton_problem_set_output_fn (
        rhea_inversion_newton_output_prestep_fn, newton_problem);
    rhea_newton_problem_set_checks (
        rhea_inversion_check_gradient, rhea_inversion_check_hessian,
        newton_problem);
  }

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);

  return inv_problem;
}

void
rhea_inversion_destroy (rhea_inversion_problem_t *inv_problem)
{
  RHEA_GLOBAL_PRODUCTION_FN_BEGIN (__func__);

  /* destroy solver data */
  if (rhea_inversion_solver_data_exists (inv_problem)) {
    rhea_inversion_newton_clear_solver_data_fn (inv_problem);
  }

  /* destroy Newton problem */
  rhea_newton_problem_destroy (inv_problem->newton_problem);

  /* destroy vectors */
  ymir_vec_destroy (inv_problem->newton_neg_gradient_vec);
  ymir_vec_destroy (inv_problem->newton_step_vec);

  /* destroy parameters */
  rhea_inversion_param_destroy (inv_problem->inv_param);

  /* destroy structure */
  RHEA_FREE (inv_problem);

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);
}

/******************************************************************************
 * Inverse Problem Solver
 *****************************************************************************/

void
rhea_inversion_solve (rhea_inversion_problem_t *inv_problem)
{
  RHEA_GLOBAL_PRODUCTION_FN_BEGIN (__func__);

  //TODO

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);
}

void
rhea_inversion_solve_with_vel_obs (rhea_inversion_problem_t *inv_problem,
                                   ymir_vec_t *vel_obs_surf,
                                   ymir_vec_t *vel_obs_weight_surf)
{
  ymir_vec_t         *parameter_vec =
    rhea_inversion_param_get_vector (inv_problem->inv_param);
#ifdef RHEA_ENABLE_DEBUG
  ymir_mesh_t        *ymir_mesh;
#endif

  RHEA_GLOBAL_PRODUCTION_FN_BEGIN (__func__);

  /* create solver data */
  rhea_inversion_newton_create_solver_data_fn (parameter_vec, inv_problem);

  /* check input */
#ifdef RHEA_ENABLE_DEBUG
  ymir_mesh = ymir_vec_get_mesh (inv_problem->vel_obs_surf);
  RHEA_ASSERT (vel_obs_surf != NULL);
  RHEA_ASSERT (rhea_velocity_surface_check_vec_type (vel_obs_surf));
  RHEA_ASSERT (ymir_vec_get_mesh (vel_obs_surf) == ymir_mesh);
  RHEA_ASSERT (vel_obs_weight_surf == NULL ||
               ymir_vec_get_mesh (vel_obs_weight_surf) == ymir_mesh);
#endif

  /* override observational data */
  ymir_vec_copy (vel_obs_surf, inv_problem->vel_obs_surf);
  if (vel_obs_weight_surf != NULL) {
    ymir_vec_copy (vel_obs_weight_surf, inv_problem->vel_obs_weight_surf);
  }

  /* run solver */
  rhea_inversion_solve (inv_problem);

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);
}
