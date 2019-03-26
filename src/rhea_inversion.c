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
#define RHEA_INVERSION_DEFAULT_FIRST_ORDER_HESSIAN_APPROX (1)
#define RHEA_INVERSION_DEFAULT_ASSEMBLE_HESSIAN_MATRIX (1)
#define RHEA_INVERSION_DEFAULT_FORWARD_SOLVER_ITER_MAX (100)
#define RHEA_INVERSION_DEFAULT_FORWARD_SOLVER_RTOL (1.0e-6)
#define RHEA_INVERSION_DEFAULT_ADJOINT_SOLVER_ITER_MAX (100)
#define RHEA_INVERSION_DEFAULT_ADJOINT_SOLVER_RTOL (1.0e-6)
#define RHEA_INVERSION_DEFAULT_INCREMENTAL_FORWARD_SOLVER_ITER_MAX (100)
#define RHEA_INVERSION_DEFAULT_INCREMENTAL_FORWARD_SOLVER_RTOL (1.0e-6)
#define RHEA_INVERSION_DEFAULT_INCREMENTAL_ADJOINT_SOLVER_ITER_MAX (100)
#define RHEA_INVERSION_DEFAULT_INCREMENTAL_ADJOINT_SOLVER_RTOL (1.0e-6)
#define RHEA_INVERSION_DEFAULT_CHECK_GRADIENT (0)
#define RHEA_INVERSION_DEFAULT_CHECK_HESSIAN (0)

/* initialize options */
int                 rhea_inversion_vel_obs_type =
                      RHEA_INVERSION_DEFAULT_VEL_OBS_TYPE;
int                 rhea_inversion_first_order_hessian_approx =
                      RHEA_INVERSION_DEFAULT_FIRST_ORDER_HESSIAN_APPROX;
int                 rhea_inversion_assemble_hessian_matrix =
                      RHEA_INVERSION_DEFAULT_ASSEMBLE_HESSIAN_MATRIX;
int                 rhea_inversion_forward_solver_iter_max =
                      RHEA_INVERSION_DEFAULT_FORWARD_SOLVER_ITER_MAX;
double              rhea_inversion_forward_solver_rtol =
                      RHEA_INVERSION_DEFAULT_FORWARD_SOLVER_RTOL;
int                 rhea_inversion_adjoint_solver_iter_max =
                      RHEA_INVERSION_DEFAULT_ADJOINT_SOLVER_ITER_MAX;
double              rhea_inversion_adjoint_solver_rtol =
                      RHEA_INVERSION_DEFAULT_ADJOINT_SOLVER_RTOL;
int                 rhea_inversion_incremental_forward_solver_iter_max =
                      RHEA_INVERSION_DEFAULT_INCREMENTAL_FORWARD_SOLVER_ITER_MAX;
double              rhea_inversion_incremental_forward_solver_rtol =
                      RHEA_INVERSION_DEFAULT_INCREMENTAL_FORWARD_SOLVER_RTOL;
int                 rhea_inversion_incremental_adjoint_solver_iter_max =
                      RHEA_INVERSION_DEFAULT_INCREMENTAL_ADJOINT_SOLVER_ITER_MAX;
double              rhea_inversion_incremental_adjoint_solver_rtol =
                      RHEA_INVERSION_DEFAULT_INCREMENTAL_ADJOINT_SOLVER_RTOL;
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

  /* outer solver options (inversion) */
  YMIR_OPTIONS_B, "first-order-hessian-approximation", '\0',
    &(rhea_inversion_first_order_hessian_approx),
    RHEA_INVERSION_DEFAULT_FIRST_ORDER_HESSIAN_APPROX,
    "1st-oder (or Gauss-Newton) approx. of the Hessian drops higher derivatives",
  YMIR_OPTIONS_B, "assemble-hessian-matrix", '\0',
    &(rhea_inversion_assemble_hessian_matrix),
    RHEA_INVERSION_DEFAULT_ASSEMBLE_HESSIAN_MATRIX,
    "Toggle assembly of the Hessian matrix",

  /* inner solver options (i.e., forward/adjoint) */
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
  YMIR_OPTIONS_I, "incremental-forward-solver-iter-max", '\0',
    &(rhea_inversion_incremental_forward_solver_iter_max),
    RHEA_INVERSION_DEFAULT_INCREMENTAL_FORWARD_SOLVER_ITER_MAX,
    "Incremental forward solver: Maximum number of iterations for Stokes solver",
  YMIR_OPTIONS_D, "incremental-forward-solver-rtol", '\0',
    &(rhea_inversion_incremental_forward_solver_rtol),
    RHEA_INVERSION_DEFAULT_INCREMENTAL_FORWARD_SOLVER_RTOL,
    "Incremental forward solver: Relative tolerance for Stokes solver",
  YMIR_OPTIONS_I, "incremental-adjoint-solver-iter-max", '\0',
    &(rhea_inversion_incremental_adjoint_solver_iter_max),
    RHEA_INVERSION_DEFAULT_INCREMENTAL_ADJOINT_SOLVER_ITER_MAX,
    "Incremental adjoint solver: Maximum number of iterations for Stokes solver",
  YMIR_OPTIONS_D, "incremental-adjoint-solver-rtol", '\0',
    &(rhea_inversion_incremental_adjoint_solver_rtol),
    RHEA_INVERSION_DEFAULT_INCREMENTAL_ADJOINT_SOLVER_RTOL,
    "Incremental adjoint solver: Relative tolerance for Stokes solver",

  /* derivative checks */
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
  ymir_vec_t         *adjoint_vel_press;
  int                 forward_is_outdated;
  int                 adjoint_is_outdated;

  /* incremental forward and adjoint states */
  ymir_vec_t         *incremental_forward_vel_press;
  ymir_vec_t         *incremental_adjoint_vel_press;
  int                 incremental_forward_is_outdated;
  int                 incremental_adjoint_is_outdated;

  /* assembled Hessian matrix */
  sc_dmatrix_t       *hessian_matrix;

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
                          inv_problem->incremental_forward_vel_press != NULL &&
                          inv_problem->incremental_adjoint_vel_press != NULL);
  vel_obs_exist = (inv_problem->vel_obs_surf != NULL);

  return (fwd_adj_states_exist && vel_obs_exist);
}

/******************************************************************************
 * State Solvers
 *****************************************************************************/

/**
 * Runs solver for the forward state.
 */
static void
rhea_inversion_inner_solve_forward (rhea_inversion_problem_t *inv_problem)
{
  rhea_stokes_problem_t  *stokes_problem = inv_problem->stokes_problem;
  const int           nonzero_init = 0; //TODO 1?
  const int           iter_max = rhea_inversion_forward_solver_iter_max;
  const double        rtol = rhea_inversion_forward_solver_rtol;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* update Stokes right-hand-side for forward problem */
  if (!inv_problem->adjoint_is_outdated ||
      !inv_problem->incremental_forward_is_outdated ||
      !inv_problem->incremental_adjoint_is_outdated) {
    rhea_stokes_problem_update_solver (
        stokes_problem, 0 /* !update coeff */, NULL /* since no coeff update */,
        0 /* !override_rhs */, NULL, NULL, NULL, NULL);
  }

  /* run Stokes solver */
  rhea_stokes_problem_solve (
      &(inv_problem->forward_vel_press), nonzero_init, iter_max, rtol,
      stokes_problem);
  inv_problem->forward_is_outdated = 0;

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

/**
 * Runs solver for the adjoint state.
 */
static void
rhea_inversion_inner_solve_adjoint (rhea_inversion_problem_t *inv_problem)
{
  rhea_stokes_problem_t  *stokes_problem = inv_problem->stokes_problem;
  const int           nonzero_init = 0; //TODO 1?
  const int           iter_max = rhea_inversion_adjoint_solver_iter_max;
  const double        rtol = rhea_inversion_adjoint_solver_rtol;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* run Stokes solver */
  rhea_stokes_problem_solve_ext (
      &(inv_problem->adjoint_vel_press), nonzero_init, iter_max, rtol,
      stokes_problem, 1 /* force linear solve */);
  inv_problem->adjoint_is_outdated = 0;

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

/**
 * Runs solver for the incremental forward state.
 */
static void
rhea_inversion_inner_solve_incremental_forward (
                                        rhea_inversion_problem_t *inv_problem)
{
  rhea_stokes_problem_t  *stokes_problem = inv_problem->stokes_problem;
  const int           nonzero_init = 0; //TODO 1?
  const int           iter_max =
                        rhea_inversion_incremental_forward_solver_iter_max;
  const double        rtol =
                        rhea_inversion_incremental_forward_solver_rtol;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* run Stokes solver */
  rhea_stokes_problem_solve_ext (
      &(inv_problem->incremental_forward_vel_press), nonzero_init,
      iter_max, rtol, stokes_problem, 1 /* force linear solve */);
  inv_problem->incremental_forward_is_outdated = 0;

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

/**
 * Runs solver for the incremental adjoint state.
 */
static void
rhea_inversion_inner_solve_incremental_adjoint (
                                        rhea_inversion_problem_t *inv_problem)
{
  rhea_stokes_problem_t  *stokes_problem = inv_problem->stokes_problem;
  const int           nonzero_init = 0; //TODO 1?
  const int           iter_max =
                        rhea_inversion_incremental_adjoint_solver_iter_max;
  const double        rtol =
                        rhea_inversion_incremental_adjoint_solver_rtol;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* run Stokes solver */
  rhea_stokes_problem_solve_ext (
      &(inv_problem->incremental_adjoint_vel_press), nonzero_init,
      iter_max, rtol, stokes_problem, 1 /* force linear solve */);
  inv_problem->incremental_adjoint_is_outdated = 0;

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

/******************************************************************************
 * Callback Functions for Newton Solver
 *****************************************************************************/

/**
 * Updates inversion parameters in the Stokes model.
 */
static void
rhea_inversion_update_param (ymir_vec_t *solution,
                             rhea_inversion_param_t *inv_param)
{
  const int           solution_exists = (solution != NULL);
  ymir_vec_t         *parameter_vec;
#ifdef RHEA_ENABLE_DEBUG
  ymir_vec_t         *parameter_vec_prev =
                        rhea_inversion_param_vec_new (inv_param);

  /* get parameters before update */
  rhea_inversion_param_pull_from_model (parameter_vec_prev, inv_param);
#endif

  /* get parameter vector */
  if (solution_exists) { /* if parameters are given */
    parameter_vec = solution;
  }
  else { /* otherwise use zeros as parameters */
    parameter_vec = rhea_inversion_param_vec_new (inv_param);
    ymir_vec_set_zero (parameter_vec);
  }

  /* update parameters in Stokes model */
  rhea_inversion_param_push_to_model (parameter_vec, inv_param);

  /* print previous and new parameters */
#ifdef RHEA_ENABLE_DEBUG
  {
    const int          *active = rhea_inversion_param_get_active (inv_param);
    const double       *param = parameter_vec->meshfree->e[0];
    const double       *param_prev = parameter_vec_prev->meshfree->e[0];
    int                 i;

    RHEA_GLOBAL_VERBOSE ("========================================\n");
    RHEA_GLOBAL_VERBOSEF ("%s\n", __func__);
    RHEA_GLOBAL_VERBOSE ("----------------------------------------\n");
    for (i = 0; i < rhea_inversion_param_get_n_parameters (inv_param); i++) {
      if (active[i]) {
        RHEA_GLOBAL_VERBOSEF ("param# %3i: %.15e -> %.15e\n", i,
                              param_prev[i], param[i]);
      }
    }
    RHEA_GLOBAL_VERBOSE ("========================================\n");
  }
  rhea_inversion_param_vec_destroy (parameter_vec_prev);
#endif

  /* destroy */
  if (!solution_exists) {
    rhea_inversion_param_vec_destroy (parameter_vec);
  }
}

static void
rhea_inversion_newton_create_solver_data_fn (ymir_vec_t *solution, void *data)
{
  rhea_inversion_problem_t *inv_problem = data;
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;
  rhea_stokes_problem_t    *stokes_problem = inv_problem->stokes_problem;
  ymir_mesh_t              *ymir_mesh;
  ymir_pressure_elem_t     *press_elem;

  /* exit if nothing to do */
  if (rhea_inversion_solver_data_exists (inv_problem)) {
    return;
  }

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (solution == NULL ||
               rhea_inversion_param_vec_check_type (solution, inv_param));
  RHEA_ASSERT (solution == NULL ||
               rhea_inversion_param_vec_is_valid (solution, inv_param));

  /* get mesh data */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);

  /* create observational data */
  RHEA_ASSERT (inv_problem->vel_obs_surf == NULL);
  RHEA_ASSERT (inv_problem->vel_obs_weight_surf == NULL);
  inv_problem->vel_obs_surf = rhea_velocity_surface_new (ymir_mesh);
  inv_problem->vel_obs_weight_surf = ymir_face_cvec_new (
      ymir_mesh, RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);

  //TODO fill vel obs type, obs vec, and weight
  ymir_vec_set_value (inv_problem->vel_obs_surf, NAN);
  ymir_vec_set_value (inv_problem->vel_obs_weight_surf, 1.0);

  /* create (forward/adjoint) state variables */
  RHEA_ASSERT (inv_problem->forward_vel_press == NULL);
  RHEA_ASSERT (inv_problem->adjoint_vel_press == NULL);
  RHEA_ASSERT (inv_problem->incremental_forward_vel_press == NULL);
  RHEA_ASSERT (inv_problem->incremental_adjoint_vel_press == NULL);
  inv_problem->forward_vel_press =
    rhea_velocity_pressure_new (ymir_mesh, press_elem);
  inv_problem->adjoint_vel_press =
    rhea_velocity_pressure_new (ymir_mesh, press_elem);
  inv_problem->incremental_forward_vel_press =
    rhea_velocity_pressure_new (ymir_mesh, press_elem);
  inv_problem->incremental_adjoint_vel_press =
    rhea_velocity_pressure_new (ymir_mesh, press_elem);

  /* create Hessian matrix */
  if (rhea_inversion_assemble_hessian_matrix) {
    const int           n = rhea_inversion_param_get_n_active (inv_param);

    RHEA_ASSERT (inv_problem->hessian_matrix == NULL);
    inv_problem->hessian_matrix = sc_dmatrix_new (n, n);
  }

  /* update inversion parameters */
  rhea_inversion_update_param (solution, inv_param);

  /* set up Stokes solver for forward problem */
  rhea_stokes_problem_setup_solver (stokes_problem);

  /* reset flags for all states */
  inv_problem->forward_is_outdated = 1;
  inv_problem->adjoint_is_outdated = 1;
  inv_problem->incremental_forward_is_outdated = 1;
  inv_problem->incremental_adjoint_is_outdated = 1;

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static void
rhea_inversion_newton_clear_solver_data_fn (void *data)
{
  rhea_inversion_problem_t *inv_problem = data;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (rhea_inversion_solver_data_exists (inv_problem));

  /* destroy (forward/adjoint) state variables */
  rhea_velocity_pressure_destroy (inv_problem->forward_vel_press);
  rhea_velocity_pressure_destroy (inv_problem->adjoint_vel_press);
  rhea_velocity_pressure_destroy (inv_problem->incremental_forward_vel_press);
  rhea_velocity_pressure_destroy (inv_problem->incremental_adjoint_vel_press);
  inv_problem->forward_vel_press = NULL;
  inv_problem->adjoint_vel_press = NULL;
  inv_problem->incremental_forward_vel_press = NULL;
  inv_problem->incremental_adjoint_vel_press = NULL;

  /* destroy Hessian matrix */
  if (inv_problem->hessian_matrix != NULL) {
    sc_dmatrix_destroy (inv_problem->hessian_matrix);
    inv_problem->hessian_matrix = NULL;
  }

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
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;
  rhea_stokes_problem_t    *stokes_problem = inv_problem->stokes_problem;

  //TODO add performance monitoring

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (solution == NULL ||
               rhea_inversion_param_vec_check_type (solution, inv_param));
  RHEA_ASSERT (solution == NULL ||
               rhea_inversion_param_vec_is_valid (solution, inv_param));

  /* solve for forward state */
  if (inv_problem->forward_is_outdated) {
    rhea_inversion_inner_solve_forward (inv_problem);
  }

  /* update inversion parameters */
  rhea_inversion_update_param (solution, inv_param);

  /* update Stokes solver for forward problem
   * Note: The Stokes coefficient needs only to be updated for linear Stokes,
   *       since the nonlinear solver will re-compute the coefficient later. */
  RHEA_ASSERT (!inv_problem->forward_is_outdated);
  rhea_stokes_problem_update_solver (
      stokes_problem,
      !rhea_stokes_problem_is_nonlinear (stokes_problem) /* update if lin. */,
      inv_problem->forward_vel_press /* velocity-pressure for coeff update */,
      0 /* !override_rhs */, NULL, NULL, NULL, NULL);

  /* reset flags for all states */
  inv_problem->forward_is_outdated = 1;
  inv_problem->adjoint_is_outdated = 1;
  inv_problem->incremental_forward_is_outdated = 1;
  inv_problem->incremental_adjoint_is_outdated = 1;

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static void         rhea_inversion_newton_apply_hessian (
                                        ymir_vec_t *param_vec_out,
                                        ymir_vec_t *param_vec_in,
                                        rhea_inversion_problem_t *inv_problem,
                                        const int first_order_approx);

static void
rhea_inversion_newton_update_hessian_fn (ymir_vec_t *solution,
                                         ymir_vec_t *step_vec,
                                         const double step_length,
                                         void *data)
{
  rhea_inversion_problem_t *inv_problem = data;
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;
  ymir_vec_t         *param_vec_in;
  ymir_vec_t         *param_vec_out;
  sc_dmatrix_t       *param_vec_reduced;
  sc_dmatrix_t       *hessian_mat = inv_problem->hessian_matrix;
  const int           n_parameters =
                        rhea_inversion_param_get_n_parameters (inv_param);
  const int          *active = rhea_inversion_param_get_active (inv_param);
  int                 i, row, col;

  //TODO add performance monitoring

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (
      __func__, "first_order_approx=%i, assemble_matrix=%i",
      rhea_inversion_first_order_hessian_approx,
      rhea_inversion_assemble_hessian_matrix);

  /* check input */
  RHEA_ASSERT (!rhea_inversion_assemble_hessian_matrix ||
               inv_problem->hessian_matrix != NULL);

  /* exit if nothing to do */
  if (!rhea_inversion_assemble_hessian_matrix) {
    RHEA_GLOBAL_VERBOSE_FN_END (__func__);
    return;
  }

  /* check input */
  RHEA_ASSERT (hessian_mat->m == rhea_inversion_param_get_n_active (inv_param));
  RHEA_ASSERT (hessian_mat->n == rhea_inversion_param_get_n_active (inv_param));

  /* create work vectors */
  param_vec_in = rhea_inversion_param_vec_new (inv_param);
  param_vec_out = rhea_inversion_param_vec_new (inv_param);

  /* assemble Hessian matrix */
  row = 0;
  col = 0;
  for (i = 0; i < n_parameters; i++) { /* loop over all parameters */
    if (active[i]) {
      /* set unit vector */
      ymir_vec_set_zero (param_vec_in);
      param_vec_in->meshfree->e[0][i] = 1.0;

      /* apply unit vector to Hessian */
      rhea_inversion_newton_apply_hessian (
          param_vec_out, param_vec_in, inv_problem,
          rhea_inversion_first_order_hessian_approx);
      rhea_inversion_param_vec_is_valid (param_vec_out, inv_param);

      /* fill one column of the Hessian matrix */
      param_vec_reduced = rhea_inversion_param_vec_reduced_new (param_vec_out,
                                                                inv_param);
      RHEA_ASSERT (param_vec_reduced->m == hessian_mat->m);
      RHEA_ASSERT (param_vec_reduced->n == 1);
      RHEA_ASSERT (col < hessian_mat->n);
      for (row = 0; row < param_vec_reduced->m; row++) {
        hessian_mat->e[row][col] = param_vec_reduced->e[row][0];
      }
      rhea_inversion_param_vec_reduced_destroy (param_vec_reduced);

      /* increment column index of matrix */
      col++;
    }
  }

  /* print Hessian */
#ifdef RHEA_ENABLE_DEBUG
  RHEA_GLOBAL_VERBOSE ("========================================\n");
  RHEA_GLOBAL_VERBOSEF ("%s\n", __func__);
  RHEA_GLOBAL_VERBOSE ("----------------------------------------\n");
  for (row = 0; row < hessian_mat->m; row++) {
    char                line[BUFSIZ] = "";
    char               *pos = line;
    const char         *end = line + BUFSIZ;

    for (col = 0; col < hessian_mat->n; col++) {
      pos += snprintf (pos, end - pos, "  %+.6e", hessian_mat->e[row][col]);
    }
    RHEA_GLOBAL_VERBOSEF ("%s\n", line);
  }
  RHEA_GLOBAL_VERBOSE ("========================================\n");
#endif

  /* check output */
  sc_dmatrix_is_valid (hessian_mat);

  /* destroy */
  rhea_inversion_param_vec_destroy (param_vec_in);
  rhea_inversion_param_vec_destroy (param_vec_out);

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static double
rhea_inversion_newton_evaluate_objective_fn (ymir_vec_t *solution, void *data)
{
  rhea_inversion_problem_t *inv_problem = data;
  rhea_stokes_problem_t    *stokes_problem = inv_problem->stokes_problem;
  ymir_mesh_t              *ymir_mesh;
  ymir_pressure_elem_t     *press_elem;
  ymir_vec_t         *vel;
  double              obj_val;

  //TODO add performance monitoring

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* solve for forward state */
  if (inv_problem->forward_is_outdated) {
    rhea_inversion_inner_solve_forward (inv_problem);
  }

  /* get mesh data */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);

  /* retrieve velocity of the forward state */
  vel = rhea_velocity_new (ymir_mesh);
  rhea_velocity_pressure_copy_components (
      vel, NULL, inv_problem->forward_vel_press, press_elem);
  rhea_stokes_problem_velocity_boundary_set_zero (vel, stokes_problem);
  RHEA_ASSERT (rhea_velocity_is_valid (vel));

  /* compute misfit term, e.g., (squared norm of the) observation data misfit */
  obj_val = rhea_inversion_obs_velocity_misfit (
      vel, inv_problem->vel_obs_surf, inv_problem->vel_obs_weight_surf,
      inv_problem->vel_obs_type,
      rhea_stokes_problem_get_domain_options (stokes_problem));
  obj_val *= 0.5;
  rhea_velocity_destroy (vel);

  /* compute & add prior term */
  //TODO

  RHEA_GLOBAL_VERBOSEF_FN_END (__func__, "result=%g", obj_val);

  /* return value of objective functional */
  return obj_val;
}

static void
rhea_inversion_newton_compute_negative_gradient_fn (
                                            ymir_vec_t *neg_gradient,
                                            ymir_vec_t *solution, void *data)
{
  rhea_inversion_problem_t *inv_problem = data;
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;
  rhea_stokes_problem_t    *stokes_problem = inv_problem->stokes_problem;
  ymir_mesh_t              *ymir_mesh;
  ymir_pressure_elem_t     *press_elem;
  ymir_vec_t         *vel, *rhs_vel_mass, *rhs_vel_press;

  //TODO add performance monitoring

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (neg_gradient, inv_param));

  /* solve for forward state */
  if (inv_problem->forward_is_outdated) {
    rhea_inversion_inner_solve_forward (inv_problem);
  }

  /* get mesh data */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);

  /* create work variables */
  vel = rhea_velocity_new (ymir_mesh);
  rhs_vel_mass = rhea_velocity_new (ymir_mesh);
  rhs_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);

  /* retrieve velocity of the forward state */
  rhea_velocity_pressure_copy_components (
      vel, NULL, inv_problem->forward_vel_press, press_elem);
  rhea_stokes_problem_velocity_boundary_set_zero (vel, stokes_problem);
  RHEA_ASSERT (rhea_velocity_is_valid (vel));

  /* set the right-hand side of the momentum eq. for the adjoint problem */
  rhea_inversion_obs_velocity_adjoint_rhs (
      rhs_vel_mass, vel, inv_problem->vel_obs_surf,
      inv_problem->vel_obs_weight_surf, inv_problem->vel_obs_type,
      rhea_stokes_problem_get_domain_options (stokes_problem));

  /* enforce Dirichlet BC's on right-hand side */
  rhea_stokes_problem_velocity_boundary_set_zero (rhs_vel_mass,
                                                  stokes_problem);

  /* set Stokes right-hand side for the adjoint problem */
  ymir_vec_set_zero (rhs_vel_press);
  rhea_velocity_pressure_set_components (rhs_vel_press, rhs_vel_mass, NULL,
                                         press_elem);

  /* update Stokes solver for adjoint problem
   * Note: The Stokes coefficient is updated only for nonlinear Stokes, since
   *       the coefficient depends on the forward velocity. */
  rhea_stokes_problem_update_solver (
      stokes_problem,
      rhea_stokes_problem_is_nonlinear (stokes_problem) /* update if nl. */,
      inv_problem->forward_vel_press /* velocity-pressure for coeff update */,
      1 /* override_rhs */, rhs_vel_press, NULL, NULL, NULL);

  /* destroy */
  rhea_velocity_destroy (vel);
  rhea_velocity_destroy (rhs_vel_mass);
  rhea_velocity_pressure_destroy (rhs_vel_press);

  /* solve for adjoint state */
  rhea_inversion_inner_solve_adjoint (inv_problem);

  /* compute the (positive) gradient */
  rhea_inversion_param_compute_gradient (
      neg_gradient, inv_problem->forward_vel_press,
      inv_problem->adjoint_vel_press, inv_param);

  /* print gradient */
#ifdef RHEA_ENABLE_DEBUG
  RHEA_GLOBAL_VERBOSE ("========================================\n");
  RHEA_GLOBAL_VERBOSE ("Inversion gradient\n");
  RHEA_GLOBAL_VERBOSE ("----------------------------------------\n");
  rhea_inversion_param_print (neg_gradient, inv_param);
  RHEA_GLOBAL_VERBOSE ("========================================\n");
#endif

  /* flip the sign of the (positive) gradient */
  ymir_vec_scale (-1.0, neg_gradient);

  /* check output */
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (neg_gradient, inv_param));

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static double
rhea_inversion_newton_compute_gradient_norm_fn (ymir_vec_t *neg_gradient,
                                                void *data, double *norm_comp)
{
  rhea_inversion_problem_t *inv_problem = data;
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;
  double              grad_norm;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (neg_gradient, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (neg_gradient, inv_param));

  /* compute norm */
  grad_norm = rhea_inversion_param_compute_gradient_norm (neg_gradient,
                                                          inv_param);
  RHEA_GLOBAL_VERBOSEF_FN_TAG (__func__, "gradient_norm=%g", grad_norm);

  return grad_norm;
}

static void
rhea_inversion_newton_apply_hessian (ymir_vec_t *param_vec_out,
                                     ymir_vec_t *param_vec_in,
                                     rhea_inversion_problem_t *inv_problem,
                                     const int first_order_approx)
{
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;
  rhea_stokes_problem_t    *stokes_problem = inv_problem->stokes_problem;
  ymir_mesh_t              *ymir_mesh;
  ymir_pressure_elem_t     *press_elem;
  ymir_vec_t         *vel, *rhs_vel_mass, *rhs_vel_press;

  /* solve for forward state */
  if (inv_problem->forward_is_outdated) {
    rhea_inversion_inner_solve_forward (inv_problem);
  }

  /* get mesh data */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);

  /* create work variables */
  vel = rhea_velocity_new (ymir_mesh);
  rhs_vel_mass = rhea_velocity_new (ymir_mesh);
  rhs_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);

  { /* BEGIN: incremental forward */
    /* set right-hand side of the momentum eq. for the incr. forward problem */
    RHEA_ASSERT (!inv_problem->forward_is_outdated);
    rhea_inversion_param_incremental_forward_rhs (
        rhs_vel_mass, param_vec_in, inv_problem->forward_vel_press, inv_param);

    /* enforce Dirichlet BC's on right-hand side */
    rhea_stokes_problem_velocity_boundary_set_zero (rhs_vel_mass,
                                                    stokes_problem);

    /* set Stokes right-hand side for the incremental forward problem */
    ymir_vec_set_zero (rhs_vel_press);
    rhea_velocity_pressure_set_components (rhs_vel_press, rhs_vel_mass, NULL,
                                           press_elem);

    /* update Stokes right-hand side for incremental forward problem */
    rhea_stokes_problem_update_solver (
        stokes_problem, 0 /* !update coeff */, NULL /* since no coeff update */,
        1 /* override_rhs */, rhs_vel_press, NULL, NULL, NULL);

    /* solve for incremental forward state */
    rhea_inversion_inner_solve_incremental_forward (inv_problem);
  } /* END: incremental forward */

  { /* BEGIN: incremental adjoint */
    /* retrieve velocity of the incremental forward state */
    rhea_velocity_pressure_copy_components (
        vel, NULL, inv_problem->incremental_forward_vel_press, press_elem);
    rhea_stokes_problem_velocity_boundary_set_zero (vel, stokes_problem);
    RHEA_ASSERT (rhea_velocity_is_valid (vel));

    /* set right-hand side of the momentum eq. for the incr. adjoint problem */
    rhea_inversion_obs_velocity_incremental_adjoint_rhs (
        rhs_vel_mass, vel, inv_problem->vel_obs_type,
        rhea_stokes_problem_get_domain_options (stokes_problem));

    /* add 2nd-order terms to right-hand side */
    if (!first_order_approx) { /* if full Hessian */
      /* retrieve velocity of the adjont state */
      RHEA_ASSERT (!inv_problem->adjoint_is_outdated);
      rhea_velocity_pressure_copy_components (
          vel, NULL, inv_problem->adjoint_vel_press, press_elem);
      rhea_stokes_problem_velocity_boundary_set_zero (vel, stokes_problem);
      RHEA_ASSERT (rhea_velocity_is_valid (vel));

      RHEA_ABORT_NOT_REACHED (); //TODO
    }

    /* enforce Dirichlet BC's on right-hand side */
    rhea_stokes_problem_velocity_boundary_set_zero (rhs_vel_mass,
                                                    stokes_problem);

    /* set Stokes right-hand side for incremental adjoint problem */
    ymir_vec_set_zero (rhs_vel_press);
    rhea_velocity_pressure_set_components (rhs_vel_press, rhs_vel_mass, NULL,
                                           press_elem);

    /* update Stokes right-hand side for incremental adjoint problem */
    rhea_stokes_problem_update_solver (
        stokes_problem, 0 /* !update coeff */, NULL /* since no coeff update */,
        1 /* override_rhs */, rhs_vel_press, NULL, NULL, NULL);

    /* solve for incremental adjoint state */
    rhea_inversion_inner_solve_incremental_adjoint (inv_problem);
  } /* END: incremental adjoint */

  /* destroy */
  rhea_velocity_destroy (vel);
  rhea_velocity_destroy (rhs_vel_mass);
  rhea_velocity_pressure_destroy (rhs_vel_press);

  /* compute the Hessian application */
  if (!first_order_approx) { /* if full Hessian */
    rhea_inversion_param_apply_hessian (
        param_vec_out, param_vec_in,
        inv_problem->forward_vel_press,
        inv_problem->adjoint_vel_press,
        inv_problem->incremental_forward_vel_press,
        inv_problem->incremental_adjoint_vel_press, inv_param);
  }
  else {
    rhea_inversion_param_apply_hessian (
        param_vec_out, param_vec_in,
        inv_problem->forward_vel_press, NULL /* adjoint */,
        NULL /* incr. forward */, inv_problem->incremental_adjoint_vel_press,
        inv_param);
  }
}

static void
rhea_inversion_newton_apply_hessian_fn (ymir_vec_t *out, ymir_vec_t *in,
                                        void *data)
{
  rhea_inversion_problem_t *inv_problem = data;
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (
      __func__, "first_order_approx=%i, assemble_matrix=%i",
      rhea_inversion_first_order_hessian_approx,
      rhea_inversion_assemble_hessian_matrix);

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (out, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (in, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (in, inv_param));
  RHEA_ASSERT (!rhea_inversion_assemble_hessian_matrix ||
               inv_problem->hessian_matrix != NULL);

  /* apply the Hessian */
  if (!rhea_inversion_assemble_hessian_matrix) { /* if call apply function */
    rhea_inversion_newton_apply_hessian (
        out, in, inv_problem, rhea_inversion_first_order_hessian_approx);
  }
  else { /* otherwise use assembled Hessian matrix */
    const sc_dmatrix_t *hessian_mat = inv_problem->hessian_matrix;
    sc_dmatrix_t       *in_reduced, *out_reduced;

    /* create reduced parameter vectors */
    in_reduced = rhea_inversion_param_vec_reduced_new (in, inv_param);
    out_reduced = rhea_inversion_param_vec_reduced_new (out, inv_param);

    /* run matrix-vector multiply; copy entries to output */
    sc_dmatrix_vector (SC_NO_TRANS, SC_NO_TRANS, SC_NO_TRANS,
                       1.0 /* factor for in vec */, hessian_mat, in_reduced,
                       0.0 /* factor for out vec */, out_reduced);
    ymir_vec_set_zero (out);
    rhea_inversion_param_vec_reduced_copy (out, out_reduced, inv_param);

    /* destroy */
    rhea_inversion_param_vec_reduced_destroy (in_reduced);
    rhea_inversion_param_vec_reduced_destroy (out_reduced);
  }
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (out, inv_param));

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
  rhea_inversion_problem_t *inv_problem = data;
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;
  int                 itn, stop_reason;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (
      __func__, "lin_iter_max=%i, lin_rtol=%.1e, nonzero_init_guess=%i, "
      "assemble_matrix=%i", lin_iter_max, lin_res_norm_rtol,
      nonzero_initial_guess, rhea_inversion_assemble_hessian_matrix);

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (step, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (neg_gradient, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (neg_gradient, inv_param));
  RHEA_ASSERT (!rhea_inversion_assemble_hessian_matrix ||
               inv_problem->hessian_matrix != NULL);

  /* (approximately) invert the Hessian */
  if (!rhea_inversion_assemble_hessian_matrix) { /* if matrix not assembled */
    RHEA_ABORT_NOT_REACHED (); //TODO implement non-assembled case
    stop_reason = 1;
    itn = 0;
  }
  else { /* if Hessian matrix is assembled */
    const sc_dmatrix_t *hessian_mat = inv_problem->hessian_matrix;
    sc_dmatrix_t       *neg_grad_reduced, *step_reduced;

    /* create reduced parameter vectors */
    neg_grad_reduced = rhea_inversion_param_vec_reduced_new (neg_gradient,
                                                             inv_param);
    step_reduced = rhea_inversion_param_vec_reduced_new (step, inv_param);

    /* run direct solver; copy entries to output */
    sc_dmatrix_ldivide (SC_NO_TRANS, hessian_mat, neg_grad_reduced,
                        step_reduced);
    ymir_vec_set_zero (step);
    rhea_inversion_param_vec_reduced_copy (step, step_reduced, inv_param);
    stop_reason = 1;
    itn = 1;

    /* destroy */
    rhea_inversion_param_vec_reduced_destroy (neg_grad_reduced);
    rhea_inversion_param_vec_reduced_destroy (step_reduced);
  }
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (step, inv_param));

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);

  /* return iteraton count and "stopping" reason */
  if (lin_iter_count != NULL) {
    *lin_iter_count = itn;
  }
  return stop_reason;
}

static void
rhea_inversion_newton_output_prestep_fn (ymir_vec_t *solution, const int iter,
                                         void *data)
{
  rhea_inversion_problem_t *inv_problem = data;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (__func__, "newton_iter=%i", iter);

  /* print active inversion parameters */
  RHEA_GLOBAL_INFO ("========================================\n");
  RHEA_GLOBAL_INFO ("Inversion parameters\n");
  RHEA_GLOBAL_INFO ("----------------------------------------\n");
  rhea_inversion_param_print (solution, inv_problem->inv_param);
  RHEA_GLOBAL_INFO ("========================================\n");

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
  inv_problem->adjoint_vel_press = NULL;
  inv_problem->incremental_forward_vel_press = NULL;
  inv_problem->incremental_adjoint_vel_press = NULL;
  inv_problem->forward_is_outdated = 1;
  inv_problem->adjoint_is_outdated = 1;
  inv_problem->incremental_forward_is_outdated = 1;
  inv_problem->incremental_adjoint_is_outdated = 1;
  inv_problem->hessian_matrix = NULL;

  /* create parameters */
  inv_problem->inv_param_options = &rhea_inversion_param_options;
  inv_problem->inv_param = rhea_inversion_param_new (
      stokes_problem, inv_problem->inv_param_options);

  /* initialize Newton problem */
  inv_problem->newton_options = &rhea_inversion_newton_options;
  inv_problem->newton_neg_gradient_vec =
    rhea_inversion_param_vec_new (inv_problem->inv_param);
  inv_problem->newton_step_vec =
    rhea_inversion_param_vec_new (inv_problem->inv_param);

  /* create Newton problem */
  {
    const int           grad_norm_n_components = 0;//TODO
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
  rhea_inversion_param_vec_destroy (inv_problem->newton_neg_gradient_vec);
  rhea_inversion_param_vec_destroy (inv_problem->newton_step_vec);

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
rhea_inversion_solve (rhea_inversion_problem_t *inv_problem,
                      const int nonzero_initial_guess)
{
  rhea_inversion_param_t *inv_param = inv_problem->inv_param;
  rhea_newton_options_t *newton_options = inv_problem->newton_options;
  rhea_newton_problem_t *newton_problem = inv_problem->newton_problem;
  ymir_vec_t         *sol_parameter_vec;

  RHEA_GLOBAL_PRODUCTIONF_FN_BEGIN (__func__, "nonzero_initial_guess=%i",
                                    nonzero_initial_guess);

  /* create vector for inversion parameters */
  sol_parameter_vec = rhea_inversion_param_vec_new (inv_param);
  if (nonzero_initial_guess) {
    rhea_inversion_param_pull_from_model (sol_parameter_vec, inv_param);
  }

  /* run Newton solver */
  newton_options->nonzero_initial_guess = nonzero_initial_guess;
  rhea_newton_solve (&sol_parameter_vec, newton_problem, newton_options);

  /* destroy */
  rhea_inversion_param_vec_destroy (sol_parameter_vec);

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);
}

void
rhea_inversion_solve_with_vel_obs (rhea_inversion_problem_t *inv_problem,
                                   const int nonzero_initial_guess,
                                   ymir_vec_t *vel_obs_surf,
                                   ymir_vec_t *vel_obs_weight_surf)
{
  RHEA_GLOBAL_PRODUCTIONF_FN_BEGIN (__func__, "nonzero_initial_guess=%i",
                                    nonzero_initial_guess);

  /* create solver data */
  if (nonzero_initial_guess) {
    rhea_inversion_param_t *inv_param = inv_problem->inv_param;
    ymir_vec_t         *param_vec = rhea_inversion_param_vec_new (inv_param);

    rhea_inversion_param_pull_from_model (param_vec, inv_param);
    rhea_inversion_newton_create_solver_data_fn (param_vec, inv_problem);
    rhea_inversion_param_vec_destroy (param_vec);
  }
  else {
    rhea_inversion_newton_create_solver_data_fn (NULL, inv_problem);
  }

  /* check input */
  RHEA_ASSERT (vel_obs_surf != NULL);
  RHEA_ASSERT (rhea_velocity_surface_check_vec_type (vel_obs_surf));
  RHEA_ASSERT (ymir_vec_get_mesh (vel_obs_surf) ==
               ymir_vec_get_mesh (inv_problem->vel_obs_surf));
  RHEA_ASSERT (vel_obs_weight_surf == NULL ||
               ymir_vec_get_mesh (vel_obs_weight_surf) ==
               ymir_vec_get_mesh (inv_problem->vel_obs_surf));

  /* override observational data */
  ymir_vec_copy (vel_obs_surf, inv_problem->vel_obs_surf);
  if (vel_obs_weight_surf != NULL) {
    RHEA_ASSERT (inv_problem->vel_obs_weight_surf != NULL);
    ymir_vec_copy (vel_obs_weight_surf, inv_problem->vel_obs_weight_surf);
  }

  /* run solver */
  rhea_inversion_solve (inv_problem, nonzero_initial_guess);

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);
}
