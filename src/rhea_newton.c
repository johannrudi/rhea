/*
 */

#include <rhea_newton.h>
#include <rhea_newton_check.h>
#include <rhea_base.h>

/* Nonlinear problem */
struct rhea_newton_problem
{
  /* vectors (not owned) */
  ymir_vec_t         *neg_gradient_vec;
  ymir_vec_t         *step_vec;

  /* callback functions for Newton algorithm */
  rhea_newton_conv_criterion_t                conv_criterion;
  rhea_newton_evaluate_objective_fn_t         evaluate_objective;
  rhea_newton_compute_negative_gradient_fn_t  compute_neg_gradient;
  rhea_newton_compute_norm_of_gradient_fn_t   compute_gradient_norm;
  int                 grad_norm_multi_components;

  rhea_newton_apply_hessian_fn_t         apply_hessian;
  rhea_newton_solve_hessian_system_fn_t  solve_hessian_sys;

  rhea_newton_update_operator_fn_t  update_operator;
  rhea_newton_update_hessian_fn_t   update_hessian;

  /* data and related callback functions */
  void                        *data;
  rhea_newton_data_init_fn_t   data_init;
  rhea_newton_data_clear_fn_t  data_clear;

  /* options */
  int                 check_gradient;
  int                 check_hessian;
};

/* Newton step */
typedef struct rhea_newton_step
{
  ymir_vec_t         *vec;
  double              length;

  int                 iter;

  double              lin_res_norm_rtol;
  double              lin_res_norm_reduction;
  int                 lin_iter_count;
  double              lin_convergence;
}
rhea_newton_step_t;

/* Newton status */
typedef struct rhea_newton_status
{
  rhea_newton_conv_criterion_t  conv_criterion;

  /* objective value at (initial, previous, and current) iteration */
  double              obj_init;
  double              obj_prev;
  double              obj_curr;

  /* reduction of objective in just one iteration (previous and current) */
  double              obj_reduction_prev;
  double              obj_reduction_curr;

  /* reduction of objective over all iterations */
  double              obj_reduction;

  /* gradient norm at (initial, previous, and current) iteration */
  double              grad_norm_init;
  double              grad_norm_prev;
  double              grad_norm_curr;

  /* components of the gradient norm */
  int                 grad_norm_multi_components;
  double             *grad_norm_init_comp;
  double             *grad_norm_prev_comp;
  double             *grad_norm_curr_comp;

  /* reduction of gradient norm in just one iteration (previous and current) */
  double              grad_norm_reduction_prev;
  double              grad_norm_reduction_curr;

  /* reduction of gradient norm over all iterations */
  double              grad_norm_reduction;
}
rhea_newton_status_t;

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_NEWTON_DEFAULT_ITER_START (0)
#define RHEA_NEWTON_DEFAULT_ITER_MAX (10)
#define RHEA_NEWTON_DEFAULT_RTOL (1.0e-6)
#define RHEA_NEWTON_DEFAULT_LIN_ITER_MAX (100)
#define RHEA_NEWTON_DEFAULT_LIN_RTOL_INIT_N_ITER (1)
#define RHEA_NEWTON_DEFAULT_LIN_RTOL_INIT (1.0e-2)
#define RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_EXPONENT (1.618)
#define RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_MAX (1.0e-3)
#define RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_MIN_ACTIVE (1)
#define RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_MIN_THRESHOLD (0.1)
#define RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_PROGRESSIVE_N_ITER (0)
#define RHEA_NEWTON_DEFAULT_STEP_SEARCH_ITER_MAX (12)
#define RHEA_NEWTON_DEFAULT_STEP_LENGTH_MIN (1.0e-5)
#define RHEA_NEWTON_DEFAULT_STEP_LENGTH_MAX (1.0)
#define RHEA_NEWTON_DEFAULT_STEP_REDUCTION (0.5)
#define RHEA_NEWTON_DEFAULT_STEP_DESCEND_CONDITION_RELAXATION (1.0e-4)
#define RHEA_NEWTON_DEFAULT_PRINT_SUMMARY (0)

#define RHEA_NEWTON_DEFAULT_NONZERO_INITIAL_GUESS (0)
#define RHEA_NEWTON_DEFAULT_ABORT_FAILED_STEP_SEARCH (0)
#define RHEA_NEWTON_DEFAULT_STATUS_VERBOSITY (0)

/* initialize options */
int                 rhea_newton_iter_start = RHEA_NEWTON_DEFAULT_ITER_START;
int                 rhea_newton_iter_max = RHEA_NEWTON_DEFAULT_ITER_MAX;
double              rhea_newton_rtol = RHEA_NEWTON_DEFAULT_RTOL;
int                 rhea_newton_lin_iter_max = RHEA_NEWTON_DEFAULT_LIN_ITER_MAX;
int                 rhea_newton_lin_rtol_init_n_iter =
                      RHEA_NEWTON_DEFAULT_LIN_RTOL_INIT_N_ITER;
double              rhea_newton_lin_rtol_init =
                      RHEA_NEWTON_DEFAULT_LIN_RTOL_INIT;
double              rhea_newton_lin_rtol_adaptive_exponent =
                      RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_EXPONENT;
double              rhea_newton_lin_rtol_adaptive_max =
                      RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_MAX;
int                 rhea_newton_lin_rtol_adaptive_min_active =
                      RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_MIN_ACTIVE;
double              rhea_newton_lin_rtol_adaptive_min_threshold =
                      RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_MIN_THRESHOLD;
int                 rhea_newton_lin_rtol_adaptive_progressive_n_iter =
                      RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_PROGRESSIVE_N_ITER;
int                 rhea_newton_step_search_iter_max =
                      RHEA_NEWTON_DEFAULT_STEP_SEARCH_ITER_MAX;
double              rhea_newton_step_length_min =
                      RHEA_NEWTON_DEFAULT_STEP_LENGTH_MIN;
double              rhea_newton_step_length_max =
                      RHEA_NEWTON_DEFAULT_STEP_LENGTH_MAX;
double              rhea_newton_step_reduction =
                      RHEA_NEWTON_DEFAULT_STEP_REDUCTION;
double              rhea_newton_step_descend_condition_relaxation =
                      RHEA_NEWTON_DEFAULT_STEP_DESCEND_CONDITION_RELAXATION;
int                 rhea_newton_print_summary =
                      RHEA_NEWTON_DEFAULT_PRINT_SUMMARY;

void
rhea_newton_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Newton";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_I, "iter-start", '\0',
    &(rhea_newton_iter_start), RHEA_NEWTON_DEFAULT_ITER_START,
    "Start at this iteration number",
  YMIR_OPTIONS_I, "iter-max", '\0',
    &(rhea_newton_iter_max), RHEA_NEWTON_DEFAULT_ITER_MAX,
    "Maximum number of iterations",
  YMIR_OPTIONS_D, "rtol", '\0',
    &(rhea_newton_rtol), RHEA_NEWTON_DEFAULT_RTOL,
    "Relative tolerance",

  YMIR_OPTIONS_I, "lin-iter-max", '\0',
    &(rhea_newton_lin_iter_max), RHEA_NEWTON_DEFAULT_LIN_ITER_MAX,
    "Linear sub-solver: Maximum number of iterations",
  YMIR_OPTIONS_I, "lin-rtol-init-n-iter", '\0',
    &(rhea_newton_lin_rtol_init_n_iter),
    RHEA_NEWTON_DEFAULT_LIN_RTOL_INIT_N_ITER,
    "Linear sub-solver: Number of iterations using initial relative tolerance",
  YMIR_OPTIONS_D, "lin-rtol-init", '\0',
    &(rhea_newton_lin_rtol_init), RHEA_NEWTON_DEFAULT_LIN_RTOL_INIT,
    "Linear sub-solver: Initial relative tolerance",
  YMIR_OPTIONS_D, "lin-rtol-adaptive-exponent", '\0',
    &(rhea_newton_lin_rtol_adaptive_exponent),
    RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_EXPONENT,
    "Linear sub-solver: Exponent for adaptive relative tolerance",
  YMIR_OPTIONS_D, "lin-rtol-adaptive-max", '\0',
    &(rhea_newton_lin_rtol_adaptive_max),
    RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_MAX,
    "Linear sub-solver: Max relative tolerance",
  YMIR_OPTIONS_I, "lin-rtol-adaptive-min-active", '\0',
    &(rhea_newton_lin_rtol_adaptive_min_active),
    RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_MIN_ACTIVE,
    "Linear sub-solver: Use threshold for min relative tolerance",
  YMIR_OPTIONS_D, "lin-rtol-adaptive-min-threshold", '\0',
    &(rhea_newton_lin_rtol_adaptive_min_threshold),
    RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_MIN_THRESHOLD,
    "Linear sub-solver: Min relative tolerance",
  YMIR_OPTIONS_I, "lin-rtol-adaptive-progressive-n-iter", '\0',
    &(rhea_newton_lin_rtol_adaptive_progressive_n_iter),
    RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_PROGRESSIVE_N_ITER,
    "Linear sub-solver: #iter (~ variance) for progressive tightening of rtol",

  YMIR_OPTIONS_I, "step-search-iter-max", '\0',
    &(rhea_newton_step_search_iter_max),
    RHEA_NEWTON_DEFAULT_STEP_SEARCH_ITER_MAX,
    "Line search: Maximum number of step search iterations",
  YMIR_OPTIONS_D, "step-length-min", '\0',
    &(rhea_newton_step_length_min), RHEA_NEWTON_DEFAULT_STEP_LENGTH_MIN,
    "Line search: Minimum allowed step length",
  YMIR_OPTIONS_D, "step-length-max", '\0',
    &(rhea_newton_step_length_max), RHEA_NEWTON_DEFAULT_STEP_LENGTH_MAX,
    "Line search: Maximum step length at start of search (usually 1)",
  YMIR_OPTIONS_D, "step-reduction", '\0',
    &(rhea_newton_step_reduction), RHEA_NEWTON_DEFAULT_STEP_REDUCTION,
    "Line search: Reduction factor for step search",
  YMIR_OPTIONS_D, "step-descend-condition-relaxation", '\0',
    &(rhea_newton_step_descend_condition_relaxation),
    RHEA_NEWTON_DEFAULT_STEP_DESCEND_CONDITION_RELAXATION,
    "Line search: Relaxation factor for the descend condition",

  YMIR_OPTIONS_I, "print-summary", '\0',
    &(rhea_newton_print_summary), RHEA_NEWTON_DEFAULT_PRINT_SUMMARY,
    "Print summary of Newton iterations",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);
}

void
rhea_newton_process_options (rhea_newton_options_t *opt)
{
  rhea_newton_options_set_defaults (opt);

  opt->iter_start = rhea_newton_iter_start;
  opt->iter_max = rhea_newton_iter_max;
  opt->rtol = rhea_newton_rtol;

  opt->lin_iter_max = rhea_newton_lin_iter_max;
  opt->lin_rtol_init_n_iter = rhea_newton_lin_rtol_init_n_iter;
  opt->lin_rtol_init = rhea_newton_lin_rtol_init;
  opt->lin_rtol_adaptive_exponent = rhea_newton_lin_rtol_adaptive_exponent;
  opt->lin_rtol_adaptive_max = rhea_newton_lin_rtol_adaptive_max;
  opt->lin_rtol_adaptive_min_active = rhea_newton_lin_rtol_adaptive_min_active;
  opt->lin_rtol_adaptive_min_threshold =
    rhea_newton_lin_rtol_adaptive_min_threshold;
  opt->lin_rtol_adaptive_progressive_n_iter =
    rhea_newton_lin_rtol_adaptive_progressive_n_iter;

  opt->step_search_iter_max = rhea_newton_step_search_iter_max;
  opt->step_length_min = rhea_newton_step_length_min;
  opt->step_length_max = rhea_newton_step_length_max;
  opt->step_reduction = rhea_newton_step_reduction;
  opt->step_descend_condition_relaxation =
    rhea_newton_step_descend_condition_relaxation;

  opt->print_summary = rhea_newton_print_summary;
}

void
rhea_newton_options_set_defaults (rhea_newton_options_t *opt)
{
  opt->nonzero_initial_guess = RHEA_NEWTON_DEFAULT_NONZERO_INITIAL_GUESS;
  opt->abort_failed_step_search = RHEA_NEWTON_DEFAULT_ABORT_FAILED_STEP_SEARCH;

  opt->iter_start = RHEA_NEWTON_DEFAULT_ITER_START;
  opt->iter_max = RHEA_NEWTON_DEFAULT_ITER_MAX;
  opt->rtol = RHEA_NEWTON_DEFAULT_RTOL;

  opt->lin_iter_max = RHEA_NEWTON_DEFAULT_LIN_ITER_MAX;
  opt->lin_rtol_init_n_iter = RHEA_NEWTON_DEFAULT_LIN_RTOL_INIT_N_ITER;
  opt->lin_rtol_init = RHEA_NEWTON_DEFAULT_LIN_RTOL_INIT;
  opt->lin_rtol_adaptive_exponent =
    RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_EXPONENT;
  opt->lin_rtol_adaptive_max = RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_MAX;
  opt->lin_rtol_adaptive_min_active =
    RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_MIN_ACTIVE;
  opt->lin_rtol_adaptive_min_threshold =
    RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_MIN_THRESHOLD;
  opt->lin_rtol_adaptive_progressive_n_iter =
    RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_PROGRESSIVE_N_ITER;

  opt->step_search_iter_max = RHEA_NEWTON_DEFAULT_STEP_SEARCH_ITER_MAX;
  opt->step_length_min = RHEA_NEWTON_DEFAULT_STEP_LENGTH_MIN;
  opt->step_length_max = RHEA_NEWTON_DEFAULT_STEP_LENGTH_MAX;
  opt->step_reduction = RHEA_NEWTON_DEFAULT_STEP_REDUCTION;
  opt->step_descend_condition_relaxation =
    RHEA_NEWTON_DEFAULT_STEP_DESCEND_CONDITION_RELAXATION;

  opt->status_verbosity = RHEA_NEWTON_DEFAULT_STATUS_VERBOSITY;
  opt->print_summary = RHEA_NEWTON_DEFAULT_PRINT_SUMMARY;
}

/******************************************************************************
 * Nonlinear Problem
 *****************************************************************************/

rhea_newton_problem_t *
rhea_newton_problem_new (
              ymir_vec_t *neg_gradient_vec,
              ymir_vec_t *step_vec,
              rhea_newton_compute_negative_gradient_fn_t compute_neg_gradient,
              rhea_newton_solve_hessian_system_fn_t solve_hessian_sys)
{
  rhea_newton_problem_t  *nl_problem = RHEA_ALLOC (rhea_newton_problem_t, 1);

  /* check input */
  RHEA_ASSERT (neg_gradient_vec != NULL);
  RHEA_ASSERT (step_vec != NULL);
  RHEA_ASSERT (compute_neg_gradient != NULL);
  RHEA_ASSERT (solve_hessian_sys != NULL);

  nl_problem->neg_gradient_vec = neg_gradient_vec;
  nl_problem->step_vec = step_vec;

  nl_problem->data = NULL;
  nl_problem->data_init = NULL;
  nl_problem->data_clear = NULL;

  nl_problem->conv_criterion = RHEA_NEWTON_CONV_CRITERION_NONE;
  nl_problem->evaluate_objective = NULL;
  nl_problem->compute_neg_gradient = compute_neg_gradient;
  nl_problem->compute_gradient_norm = NULL;
  nl_problem->grad_norm_multi_components = 0;

  nl_problem->apply_hessian = NULL;
  nl_problem->solve_hessian_sys = solve_hessian_sys;

  nl_problem->update_operator = NULL;
  nl_problem->update_hessian = NULL;

  nl_problem->check_gradient = 0;
  nl_problem->check_hessian = 0;

  return nl_problem;
}

void
rhea_newton_problem_destroy (rhea_newton_problem_t *nl_problem)
{
  RHEA_FREE (nl_problem);
}

void
rhea_newton_problem_set_data_fn (
              void *data,
              rhea_newton_data_init_fn_t data_init,
              rhea_newton_data_clear_fn_t data_clear,
              rhea_newton_problem_t *nl_problem)
{
  nl_problem->data = data;
  nl_problem->data_init = data_init;
  nl_problem->data_clear = data_clear;
}

void
rhea_newton_problem_set_conv_criterion_fn (
              rhea_newton_conv_criterion_t conv_criterion,
              rhea_newton_evaluate_objective_fn_t evaluate_objective,
              rhea_newton_compute_norm_of_gradient_fn_t compute_gradient_norm,
              const int grad_norm_multi_components,
              rhea_newton_problem_t *nl_problem)
{
#if RHEA_ENABLE_DEBUG
  switch (conv_criterion) {
  case RHEA_NEWTON_CONV_CRITERION_OBJECTIVE:
    RHEA_ASSERT (evaluate_objective != NULL);
    break;
  case RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM:
  case RHEA_NEWTON_CONV_CRITERION_RESIDUAL_NORM:
    RHEA_ASSERT (compute_gradient_norm != NULL);
    break;
  default: /* unknown criterion */
    RHEA_ABORT_NOT_REACHED ();
  }
#endif

  nl_problem->conv_criterion = conv_criterion;
  nl_problem->evaluate_objective = evaluate_objective;
  nl_problem->compute_gradient_norm = compute_gradient_norm;
  nl_problem->grad_norm_multi_components = grad_norm_multi_components;
}

void
rhea_newton_problem_set_apply_hessian_fn (
              rhea_newton_apply_hessian_fn_t apply_hessian,
              rhea_newton_problem_t *nl_problem)
{
  nl_problem->apply_hessian = apply_hessian;
}

void
rhea_newton_problem_set_update_fn (
              rhea_newton_update_operator_fn_t update_operator,
              rhea_newton_update_hessian_fn_t update_hessian,
              rhea_newton_problem_t *nl_problem)
{
  nl_problem->update_operator = update_operator;
  nl_problem->update_hessian = update_hessian;
}

void
rhea_newton_problem_set_checks (const int check_gradient,
                                const int check_hessian,
                                rhea_newton_problem_t *nl_problem)
{
  nl_problem->check_gradient = (0 < check_gradient);
  nl_problem->check_hessian = (0 < check_hessian);
}

int
rhea_newton_problem_get_check_gradient (rhea_newton_problem_t *nl_problem)
{
  return nl_problem->check_gradient;
}

void
rhea_newton_problem_set_check_gradient (const int check_gradient,
                                        rhea_newton_problem_t *nl_problem)
{
  nl_problem->check_gradient = (0 < check_gradient);
}

int
rhea_newton_problem_get_check_hessian (rhea_newton_problem_t *nl_problem)
{
  return nl_problem->check_hessian;
}

void
rhea_newton_problem_set_check_hessian (const int check_hessian,
                                       rhea_newton_problem_t *nl_problem)
{
  nl_problem->check_hessian = (0 < check_hessian);
}

/******************************************************************************
 * Access to Data and Callback Functions of Nonlinear Problem
 *****************************************************************************/

ymir_vec_t *
rhea_newton_problem_get_neg_gradient_vec (rhea_newton_problem_t *nl_problem)
{
  return nl_problem->neg_gradient_vec;
}

ymir_vec_t *
rhea_newton_problem_get_step_vec (rhea_newton_problem_t *nl_problem)
{
  return nl_problem->step_vec;
}

void *
rhea_newton_problem_get_data (rhea_newton_problem_t *nl_problem)
{
  return nl_problem->data;
}

int
rhea_newton_problem_evaluate_objective_exists (
                                            rhea_newton_problem_t *nl_problem)
{
  return (NULL != nl_problem->evaluate_objective);
}

double
rhea_newton_problem_evaluate_objective (ymir_vec_t *solution,
                                        rhea_newton_problem_t *nl_problem)
{
  RHEA_ASSERT (rhea_newton_problem_evaluate_objective_exists (nl_problem));
  return nl_problem->evaluate_objective (solution, nl_problem->data);
}

int
rhea_newton_problem_compute_neg_gradient_exists (
                                            rhea_newton_problem_t *nl_problem)
{
  return (NULL != nl_problem->compute_neg_gradient);
}

void
rhea_newton_problem_compute_neg_gradient (ymir_vec_t *neg_gradient,
                                          ymir_vec_t *solution,
                                          rhea_newton_problem_t *nl_problem)
{
  RHEA_ASSERT (rhea_newton_problem_compute_neg_gradient_exists (nl_problem));
  nl_problem->compute_neg_gradient (neg_gradient, solution, nl_problem->data);
}

int
rhea_newton_problem_apply_hessian_exists (rhea_newton_problem_t *nl_problem)
{
  return (NULL != nl_problem->apply_hessian);
}

void
rhea_newton_problem_apply_hessian (ymir_vec_t *out, ymir_vec_t *in,
                                   rhea_newton_problem_t *nl_problem)
{
  RHEA_ASSERT (rhea_newton_problem_apply_hessian_exists (nl_problem));
  nl_problem->apply_hessian (out, in, nl_problem->data);
}

int
rhea_newton_problem_update_operator_exists (rhea_newton_problem_t *nl_problem)
{
  return (NULL != nl_problem->update_operator);
}

void
rhea_newton_problem_update_operator (ymir_vec_t *solution,
                                     rhea_newton_problem_t *nl_problem)
{
  if (NULL != nl_problem->update_operator) {
    nl_problem->update_operator (solution, nl_problem->data);
  }
}

/******************************************************************************
 * Newton Step
 *****************************************************************************/

/**
 * Initializes a step object.
 */
static void
rhea_newton_step_init (rhea_newton_step_t *step,
                       ymir_vec_t *step_vec)
{
  step->vec = step_vec;
  step->length = -1.0;

  step->iter = -1;

  step->lin_res_norm_rtol = -1.0;
  step->lin_res_norm_reduction = -1.0;
  step->lin_iter_count = -1;
  step->lin_convergence = -1.0;
}

/******************************************************************************
 * Newton Status
 *****************************************************************************/

/**
 * Initializes a status object.
 */
static void
rhea_newton_status_init (rhea_newton_status_t *status,
                         const rhea_newton_conv_criterion_t conv_criterion,
                         const int grad_norm_multi_components)
{
  /* check input */
  RHEA_ASSERT (conv_criterion != RHEA_NEWTON_CONV_CRITERION_NONE);

  /* init convergence criterion */
  status->conv_criterion = conv_criterion;

  /* init valus pertaining to objective values */
  status->obj_init = NAN;
  status->obj_prev = NAN;
  status->obj_curr = NAN;

  status->obj_reduction_prev = 1.0;
  status->obj_reduction_curr = 1.0;
  status->obj_reduction = 1.0;

  /* init valus pertaining to gradient norm(s) */
  status->grad_norm_init = NAN;
  status->grad_norm_prev = NAN;
  status->grad_norm_curr = NAN;

  status->grad_norm_multi_components = grad_norm_multi_components;
  if (0 < status->grad_norm_multi_components) {
    const int           n_components = status->grad_norm_multi_components;
    int                 compid;

    status->grad_norm_init_comp = RHEA_ALLOC (double, n_components);
    status->grad_norm_prev_comp = RHEA_ALLOC (double, n_components);
    status->grad_norm_curr_comp = RHEA_ALLOC (double, n_components);
    for (compid = 0; compid < n_components; compid++) {
      status->grad_norm_init_comp[compid] = NAN;
      status->grad_norm_prev_comp[compid] = NAN;
      status->grad_norm_curr_comp[compid] = NAN;
    }
  }
  else {
    status->grad_norm_init_comp = NULL;
    status->grad_norm_prev_comp = NULL;
    status->grad_norm_curr_comp = NULL;
  }

  status->grad_norm_reduction_prev = 1.0;
  status->grad_norm_reduction_curr = 1.0;
  status->grad_norm_reduction = 1.0;
}

/**
 * Destroys content of a status object.
 */
static void
rhea_newton_status_clear (rhea_newton_status_t *status)
{
  if (status->grad_norm_init_comp != NULL) {
    RHEA_FREE (status->grad_norm_init_comp);
  }
  if (status->grad_norm_prev_comp != NULL) {
    RHEA_FREE (status->grad_norm_prev_comp);
  }
  if (status->grad_norm_curr_comp != NULL) {
    RHEA_FREE (status->grad_norm_curr_comp);
  }
}

/**
 * Sets current status.
 */
static void
rhea_newton_status_set_curr (rhea_newton_status_t *status,
                             const double obj,
                             const double grad_norm,
                             const double *grad_norm_comp)
{
  /* set values pertaining to objective */
  status->obj_curr = obj;

  if (isfinite (status->obj_prev)) {
    status->obj_reduction_curr = status->obj_curr / status->obj_prev;
  }
  if (isfinite (status->obj_init)) {
    status->obj_reduction = status->obj_curr / status->obj_init;
  }

  /* set valus pertaining to gradient norm(s) */
  status->grad_norm_curr = grad_norm;
  if (0 < status->grad_norm_multi_components && grad_norm_comp != NULL) {
    int                 compid;

    for (compid = 0; compid < status->grad_norm_multi_components; compid++) {
      status->grad_norm_curr_comp[compid] = grad_norm_comp[compid];
    }
  }

  if (isfinite (status->grad_norm_prev)) {
    status->grad_norm_reduction_curr = status->grad_norm_curr /
                                       status->grad_norm_prev;
  }
  if (isfinite (status->grad_norm_init)) {
    status->grad_norm_reduction = status->grad_norm_curr /
                                  status->grad_norm_init;
  }
}

/**
 * Copies status of the current iteration to the initial iteration.
 */
static void
rhea_newton_status_copy_curr_to_init (rhea_newton_status_t *status)
{
  /* copy values pertaining to objective */
  status->obj_init = status->obj_curr;

  /* copy valus pertaining to gradient norm(s) */
  status->grad_norm_init = status->grad_norm_curr;
  if (0 < status->grad_norm_multi_components) {
    int                 compid;

    for (compid = 0; compid < status->grad_norm_multi_components; compid++) {
      status->grad_norm_init_comp[compid] = status->grad_norm_curr_comp[compid];
    }
  }
}

/**
 * Copies status of the current iteration to the previous iteration.
 */
static void
rhea_newton_status_copy_curr_to_prev (rhea_newton_status_t *status)
{
  /* copy values pertaining to objective */
  status->obj_prev = status->obj_curr;
  status->obj_reduction_prev = status->obj_reduction_curr;

  /* copy valus pertaining to gradient norm(s) */
  status->grad_norm_prev = status->grad_norm_curr;
  status->grad_norm_reduction_prev = status->grad_norm_reduction_curr;
  if (0 < status->grad_norm_multi_components) {
    int                 compid;

    for (compid = 0; compid < status->grad_norm_multi_components; compid++) {
      status->grad_norm_prev_comp[compid] = status->grad_norm_curr_comp[compid];
    }
  }
}

/**
 * Computes and sets the current status.
 */
static void
rhea_newton_status_compute_curr (rhea_newton_status_t *status,
                                 int *neg_gradient_updated,
                                 ymir_vec_t *solution,
                                 rhea_newton_problem_t *nl_problem,
                                 const int compute_all)
{
  double              obj_val = NAN;
  double              grad_norm = NAN;
  double             *grad_norm_comp = NULL;
  ymir_vec_t         *neg_gradient = nl_problem->neg_gradient_vec;

  /* check input */
  RHEA_ASSERT (status->conv_criterion != RHEA_NEWTON_CONV_CRITERION_NONE);

  /* evaluate objective functional */
  if (RHEA_NEWTON_CONV_CRITERION_OBJECTIVE == status->conv_criterion ||
      (compute_all && nl_problem->evaluate_objective != NULL)) {
    RHEA_ASSERT (nl_problem->evaluate_objective != NULL);
    obj_val = nl_problem->evaluate_objective (solution, nl_problem->data);
  }

  /* compute gradient norm */
  if (RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM == status->conv_criterion ||
      RHEA_NEWTON_CONV_CRITERION_RESIDUAL_NORM == status->conv_criterion ||
      (compute_all && nl_problem->compute_gradient_norm != NULL)) {
    const int           n_components = status->grad_norm_multi_components;

    RHEA_ASSERT (neg_gradient != NULL);
    RHEA_ASSERT (nl_problem->compute_neg_gradient != NULL);
    RHEA_ASSERT (nl_problem->compute_gradient_norm != NULL);

    /* compute (negative) gradient */
    nl_problem->compute_neg_gradient (neg_gradient, solution,
                                      nl_problem->data);
    *neg_gradient_updated = 1;
    if (nl_problem->check_gradient) {
      rhea_newton_check_gradient (solution, nl_problem);
    }

    /* compute norm of gradient */
    if (0 < n_components) {
      grad_norm_comp = RHEA_ALLOC (double, n_components);
    }
    grad_norm = nl_problem->compute_gradient_norm (neg_gradient,
                                                   nl_problem->data,
                                                   grad_norm_comp);
    if (0 < n_components) {
      RHEA_FREE (grad_norm_comp);
    }
  }
  else {
    *neg_gradient_updated = 0;
  }

  /* set current status */
  RHEA_ASSERT (isfinite (obj_val) || isfinite (grad_norm));
  rhea_newton_status_set_curr (status, obj_val, grad_norm, grad_norm_comp);
}

/**
 * Returns initial value that determines convergence.
 */
static double
rhea_newton_status_get_conv_value_init (rhea_newton_status_t *status)
{
  double              val;

  switch (status->conv_criterion) {
  case RHEA_NEWTON_CONV_CRITERION_OBJECTIVE:
    val = status->obj_init;
    break;
  case RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM:
  case RHEA_NEWTON_CONV_CRITERION_RESIDUAL_NORM:
    val = status->grad_norm_init;
    break;
  default: /* unknown criterion */
    RHEA_ABORT_NOT_REACHED ();
    val = -1.0;
  }

  return val;
}

/**
 * Returns previous value that determines convergence.
 */
static double
rhea_newton_status_get_conv_value_prev (rhea_newton_status_t *status)
{
  double              val;

  switch (status->conv_criterion) {
  case RHEA_NEWTON_CONV_CRITERION_OBJECTIVE:
    val = status->obj_prev;
    break;
  case RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM:
  case RHEA_NEWTON_CONV_CRITERION_RESIDUAL_NORM:
    val = status->grad_norm_prev;
    break;
  default: /* unknown criterion */
    RHEA_ABORT_NOT_REACHED ();
    val = -1.0;
  }

  return val;
}

/**
 * Returns current value that determines convergence.
 */
static double
rhea_newton_status_get_conv_value_curr (rhea_newton_status_t *status)
{
  double              val;

  switch (status->conv_criterion) {
  case RHEA_NEWTON_CONV_CRITERION_OBJECTIVE:
    val = status->obj_curr;
    break;
  case RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM:
  case RHEA_NEWTON_CONV_CRITERION_RESIDUAL_NORM:
    val = status->grad_norm_curr;
    break;
  default: /* unknown criterion */
    RHEA_ABORT_NOT_REACHED ();
    val = -1.0;
  }

  return val;
}

/**
 * Returns overall reduction of the value that determines convergence.
 */
static double
rhea_newton_status_get_reduction (rhea_newton_status_t *status)
{
  double              reduction;

  switch (status->conv_criterion) {
  case RHEA_NEWTON_CONV_CRITERION_OBJECTIVE:
    reduction = status->obj_reduction;
    break;
  case RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM:
  case RHEA_NEWTON_CONV_CRITERION_RESIDUAL_NORM:
    reduction = status->grad_norm_reduction;
    break;
  default: /* unknown criterion */
    RHEA_ABORT_NOT_REACHED ();
    reduction = -1.0;
  }

  return reduction;
}

/**
 * Returns previous reduction of the value that determines convergence.
 */
static double
rhea_newton_status_get_reduction_prev (rhea_newton_status_t *status)
{
  double              reduction;

  switch (status->conv_criterion) {
  case RHEA_NEWTON_CONV_CRITERION_OBJECTIVE:
    reduction = status->obj_reduction_prev;
    break;
  case RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM:
  case RHEA_NEWTON_CONV_CRITERION_RESIDUAL_NORM:
    reduction = status->grad_norm_reduction_prev;
    break;
  default: /* unknown criterion */
    RHEA_ABORT_NOT_REACHED ();
    reduction = -1.0;
  }

  return reduction;
}

/**
 * Returns current reduction of the value that determines convergence.
 */
static double
rhea_newton_status_get_reduction_curr (rhea_newton_status_t *status)
{
  double              reduction;

  switch (status->conv_criterion) {
  case RHEA_NEWTON_CONV_CRITERION_OBJECTIVE:
    reduction = status->obj_reduction_curr;
    break;
  case RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM:
  case RHEA_NEWTON_CONV_CRITERION_RESIDUAL_NORM:
    reduction = status->grad_norm_reduction_curr;
    break;
  default: /* unknown criterion */
    RHEA_ABORT_NOT_REACHED ();
    reduction = -1.0;
  }

  return reduction;
}

/**
 * Prints current Newton status.
 */
static void
rhea_newton_status_print_curr (rhea_newton_status_t *status,
                               const int print_all,
                               const int iter,
                               const char *name)
{
  switch (status->conv_criterion) {
  case RHEA_NEWTON_CONV_CRITERION_OBJECTIVE:
    if (!print_all) {
      RHEA_GLOBAL_INFOF (
          "Newton iter %i -- %s: Objective functional %.6e, reduction %.3e\n",
          iter, name, status->obj_curr, status->obj_reduction);
    }
    else {
      RHEA_GLOBAL_INFOF (
          "Newton iter %i -- %s: Objective functional %.6e, reduction %.3e "
          "(Gradient norm %.6e, reduction %.3e)\n",
          iter, name, status->obj_curr, status->obj_reduction,
          status->grad_norm_curr, status->grad_norm_reduction);
    }
    break;
  case RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM:
    if (!print_all) {
      RHEA_GLOBAL_INFOF (
          "Newton iter %i -- %s: Gradient norm %.6e, reduction %.3e\n",
          iter, name, status->grad_norm_curr, status->grad_norm_reduction);
    }
    else {
      RHEA_GLOBAL_INFOF (
          "Newton iter %i -- %s: Gradient norm %.6e, reduction %.3e "
          "(Objective functional %.6e, reduction %.3e)\n",
          iter, name, status->grad_norm_curr, status->grad_norm_reduction,
          status->obj_curr, status->obj_reduction);
    }
    break;
  case RHEA_NEWTON_CONV_CRITERION_RESIDUAL_NORM:
    RHEA_ASSERT (!print_all);
    RHEA_GLOBAL_INFOF (
        "Newton iter %i -- %s: Residual norm %.6e, reduction %.3e\n",
        iter, name, status->grad_norm_curr, status->grad_norm_reduction);
    break;
  default: /* unknown criterion */
    RHEA_ABORT_NOT_REACHED ();
  }
}

/******************************************************************************
 * Inexact Newton--Krylov Method
 *****************************************************************************/

/**
 * Calculates the accuracy for solving the linearized system, i.e., the Krylov
 * relative tolerance (or "forcing") for the linear solver of inexact Newton.
 */
static void
rhea_newton_set_accuracy (rhea_newton_step_t *step,
                          rhea_newton_status_t *status,
                          rhea_newton_options_t *opt)
{
  const char         *this_fn_name = "rhea_newton_set_accuracy";
  const int           iter = step->iter;

  /* check input */
  RHEA_ASSERT (0 <= step->iter);

  /*
   * Initial Relative Tolerance
   */
  if (iter < opt->lin_rtol_init_n_iter) {
    /* check input */
    RHEA_ASSERT (0.0 < opt->lin_rtol_init);

    /* set initial rtol and exit function */
    step->lin_res_norm_rtol = opt->lin_rtol_init;
    RHEA_GLOBAL_INFOF ("Newton iter %i -- %s: Use initial rtol %.3e\n",
                       iter, this_fn_name, step->lin_res_norm_rtol);
    return;
  }

  /*
   * Adaptive Relative Tolerance
   */
  {
    const int           prog_reduction_active =
                          (0 < opt->lin_rtol_adaptive_progressive_n_iter);
    double              prog_reduction;

    const int           lin_rtol_min_active =
                          opt->lin_rtol_adaptive_min_active;
    const double        lin_rtol_min_thresh =
                          opt->lin_rtol_adaptive_min_threshold;

    const double        exponent = opt->lin_rtol_adaptive_exponent;
    const double        reduction_prev =
      SC_MIN (rhea_newton_status_get_reduction_prev (status), 1.0);
    const double        reduction_curr =
      SC_MIN (rhea_newton_status_get_reduction_curr (status), 1.0);

    double              lin_rtol_max;
    double              lin_rtol_min, lin_rtol_min_effective;
    double              lin_rtol;

    /* check input */
    RHEA_ASSERT (0.0 < opt->lin_rtol_adaptive_max);
    RHEA_ASSERT (opt->lin_rtol_adaptive_max < 1.0);
    RHEA_ASSERT (0.0 < rhea_newton_status_get_reduction_prev (status));
    RHEA_ASSERT (0.0 < rhea_newton_status_get_reduction_curr (status));

    /* set progressive reduction factor:
     *
     *   (exp(-(#iter / prog #iter)^2) + eps/rtol_max) /
     *                              (1 + eps/rtol_max)
     */
    if (prog_reduction_active) {
      const double        it2 = (double) iter * iter;
      const double        itn2 = (double)
                            opt->lin_rtol_adaptive_progressive_n_iter *
                            opt->lin_rtol_adaptive_progressive_n_iter;
      const double        eps = SC_1000_EPS;

      prog_reduction = (exp (-it2/itn2) + eps/lin_rtol_max) /
                       (1.0 + eps/lin_rtol_max);
      RHEA_ASSERT (0.0 < prog_reduction && prog_reduction <= 1.0);
    }
    else {
      prog_reduction = 1.0;
    }

    /* set max rtol */
    lin_rtol_max = prog_reduction * opt->lin_rtol_adaptive_max;

    /* set min rtol, which serves as a saveguard to avoid oversolving
     *
     *   rtol_max * ( (prev reduction + current reduction) / 2 )^exponent
     */
    if (lin_rtol_min_active) {
      double              reduction_avg;

      reduction_avg = 0.5 * (reduction_prev + reduction_curr);
      lin_rtol_min = lin_rtol_max * pow (reduction_avg, exponent);

      /* activate min rtol above a given threshold */
      if (lin_rtol_min_thresh < lin_rtol_min) {
        lin_rtol_min_effective = lin_rtol_min;
      }
      else {
        lin_rtol_min_effective = 0.0;
      }
    }
    else {
      lin_rtol_min = -1.0;
      lin_rtol_min_effective = 0.0;
    }

    /* set the candidate for the relative tolerance according to "Choice 2" in
     * [Eisenstat, Walker, 1996]:
     *
     *   rtol_max * (current residual / prev residual)^exponent
     */
    lin_rtol = lin_rtol_max * pow (reduction_curr, exponent);

    /* set the effective adaptive rtol */
    step->lin_res_norm_rtol = SC_MAX (lin_rtol_min_effective, lin_rtol);

    RHEA_GLOBAL_INFOF (
        "Newton iter %i -- %s: Set adaptive rtol %.3e, candidate rtol %.3e, "
        "max %.3e (w/ progressive reduction %.3e), "
        "min %.3e (w/ threshold %.3e)\n",
        iter, this_fn_name, step->lin_res_norm_rtol, lin_rtol,
        lin_rtol_max, prog_reduction,
        lin_rtol_min, lin_rtol_min_thresh);
  }
}

/**
 * Solves the linearized system for an inexact Newton step.
 */
static void
rhea_newton_compute_step (rhea_newton_step_t *step,
                          rhea_newton_problem_t *nl_problem,
                          rhea_newton_options_t *opt)
{
  const char         *this_fn_name = "rhea_newton_compute_step";
  const int           iter = step->iter;
  ymir_vec_t         *rhs = nl_problem->neg_gradient_vec;
  const int           nonzero_initial_guess = 0;
  const int           lin_iter_max = opt->lin_iter_max;
  int                 lin_iter_count = -1;
  const double        lin_res_norm_rtol = step->lin_res_norm_rtol;
  double              lin_res_norm_reduction;
  double              lin_conv;

  RHEA_GLOBAL_INFOF ("Newton iter %i -- Into %s (max iter %i, rtol %.1e)\n",
                     iter, this_fn_name, lin_iter_max, lin_res_norm_rtol);

  /* check input */
  RHEA_ASSERT (0 <= step->iter);
  RHEA_ASSERT (0.0 < step->lin_res_norm_rtol && step->lin_res_norm_rtol < 1.0);

  /* run solver for the linearized system */
  nl_problem->solve_hessian_sys (
      step->vec, rhs, lin_iter_max, lin_res_norm_rtol,
      nonzero_initial_guess, nl_problem->data, &lin_iter_count);
  RHEA_ASSERT (0 <= lin_iter_count);

  /* calculate the residual reduction of the linearized solve */
  if (nl_problem->apply_hessian != NULL) { /* if function provided */
    ymir_vec_t         *lin_residual_vec;
    double              lin_res_norm_init;
    double              lin_res_norm_curr;

    /* compute the inital l^2-norm of the residual of the linearized system */
    lin_res_norm_init = ymir_vec_norm (rhs);

    /* compute the new l^2-norm of the residual of the linearized system */
    lin_residual_vec = ymir_vec_template (rhs);
    nl_problem->apply_hessian (lin_residual_vec, step->vec, nl_problem->data);
    ymir_vec_add (-1.0, rhs, lin_residual_vec);
    lin_res_norm_curr = ymir_vec_norm (lin_residual_vec);
    ymir_vec_destroy (lin_residual_vec);

    /* calculate residual reduction */
    lin_res_norm_reduction = lin_res_norm_curr / lin_res_norm_init;

    /* calculate convergence of the linear solver `rtol^(1/#iter)` */
    lin_conv = exp ( log (lin_res_norm_reduction) / ((double) lin_iter_count) );
  }
  else {
    lin_res_norm_reduction = -1.0;
    lin_conv = -1.0;
  }

  /* set the #iterations used by the linear solver and the residual reduction */
  step->lin_iter_count = lin_iter_count;
  step->lin_res_norm_reduction = lin_res_norm_reduction;
  step->lin_convergence = lin_conv;

  RHEA_GLOBAL_INFOF ("Newton iter %i -- Done %s (num iter %i, "
                     "residual reduction %3.e, prescribed rtol %.3e)\n",
                     iter, this_fn_name, lin_iter_count,
                     lin_res_norm_reduction, lin_res_norm_rtol);
}

static int
rhea_newton_search_step_length_check_descend (
                                          rhea_newton_status_t *status,
                                          rhea_newton_options_t *opt,
                                          const double lin_res_norm_reduction,
                                          const int iter,
                                          const int print_condition)
{
  const char         *this_fn_name = "rhea_newton_search_step_length";
  const double        relax = opt->step_descend_condition_relaxation;
  double              reduction, descend_reduction;

  /* check input */
  RHEA_ASSERT (0.0 < lin_res_norm_reduction);
  RHEA_ASSERT (lin_res_norm_reduction <= 1.0);
  RHEA_ASSERT (0.0 < rhea_newton_status_get_reduction_curr (status));

  /* compute actual and required reduction */
  reduction = rhea_newton_status_get_reduction_curr (status);
  descend_reduction = 1.0 - relax * (1.0 - lin_res_norm_reduction);

  /* print */
  if (print_condition) {
    RHEA_GLOBAL_INFOF (
        "Newton iter %i -- %s: Descend condition: reduction <= %.3e\n",
        iter, this_fn_name, descend_reduction);
  }

  /* check for descend */
  return (reduction <= descend_reduction);
}

/**
 * Performs line search for the step length (updates the nonlinear operator).
 */
static int
rhea_newton_search_step_length (ymir_vec_t *solution,
                                rhea_newton_step_t *step,
                                rhea_newton_status_t *status,
                                int *neg_gradient_updated,
                                rhea_newton_problem_t *nl_problem,
                                rhea_newton_options_t *opt)
{
  const char         *this_fn_name = "rhea_newton_search_step_length";
  const int           print_all_conv_criteria = (2 <= opt->status_verbosity);
  const int           iter = step->iter;
  ymir_vec_t         *solution_prev = ymir_vec_template (solution);
  const int           search_iter_max = opt->step_search_iter_max;
  int                 search_id;
  int                 search_success = 0;
  const double        step_length_min = opt->step_length_min;
  const double        step_reduction = opt->step_reduction;

  RHEA_GLOBAL_INFOF ("Newton iter %i -- Into %s\n", iter, this_fn_name);

  /* check input */
  RHEA_ASSERT (0 <= step->iter);

  /* copy previous solution and status */
  ymir_vec_copy (solution, solution_prev);
  rhea_newton_status_copy_curr_to_prev (status);

  /* initialize step length */
  step->length = opt->step_length_max;

  /* search for step length */
  for (search_id = 1; search_id <= search_iter_max; search_id++) {
    /* move solution into step direction */
    ymir_vec_add (step->length, step->vec, solution);

    /* update the nonlinear operator at new solution */
    if (nl_problem->update_operator != NULL) {
      nl_problem->update_operator (solution, nl_problem->data);
    }

    /* update status */
    rhea_newton_status_compute_curr (status, neg_gradient_updated,
                                     solution, nl_problem,
                                     print_all_conv_criteria);

    /* check descend condition */
    search_success = rhea_newton_search_step_length_check_descend (
        status, opt, step->lin_res_norm_reduction, iter,
        (0 == search_id) /* print only once */);

    RHEA_GLOBAL_INFOF ("Newton iter %i -- %s: #%2i, step length %.3e, "
                       "convergence criterion %i, "
                       "current %.3e, previous %.3e\n",
                       iter, this_fn_name, search_id, step->length,
                       status->conv_criterion,
                       rhea_newton_status_get_conv_value_curr (status),
                       rhea_newton_status_get_conv_value_prev (status));

    /* continue or terminate loop */
    if (search_success) { /* terminate if descend condition is satisfied */
      break;
    }
    else { /* otherwise continue searching */
      /* reduce step length */
      step->length *= step_reduction;

      /* stop search if step length became too low */
      if (step->length < step_length_min) {
        step->length = 0.0;
        break;
      }

      /* put back previous solution to continue loop */
      if (search_id < search_iter_max) {
        ymir_vec_copy (solution_prev, solution);
      }
    }
  } /* END: step length search */

  /* post-processing depending on whether a step length was found */
  if (search_success) { /* if search was successful */
    RHEA_GLOBAL_INFOF (
        "Newton iter %i -- %s: Line search successful, step length %.3e\n",
        iter, this_fn_name, step->length);
  }
  else { /* if search failed */
    /* check and print reason */
    if (0.0 < step->length) {
      step->length = -1.0;
      RHEA_GLOBAL_INFOF (
          "Newton iter %i -- %s: Line search failed, "
          "max number of step reductions reached\n", iter, this_fn_name);
    }
    else {
      step->length = 0.0;
      RHEA_GLOBAL_INFOF (
          "Newton iter %i -- %s: Line search failed, "
          "min step length reached\n", iter, this_fn_name);
    }

    /* reverse updates of solution and nonlinear operator */
    if (opt->abort_failed_step_search) {
      ymir_vec_copy (solution_prev, solution);
      if (nl_problem->update_operator != NULL) {
        nl_problem->update_operator (solution, nl_problem->data);
      }
    }
  }

  /* destroy */
  ymir_vec_destroy (solution_prev);

  RHEA_GLOBAL_INFOF ("Newton iter %i -- Done %s\n", iter, this_fn_name);

  /* return success of step length search */
  return search_success;
}

void
rhea_newton_solve (ymir_vec_t *solution,
                   rhea_newton_problem_t *nl_problem,
                   rhea_newton_options_t *opt)
{
  const char         *this_fn_name = "rhea_newton_solve";
  const int           print_all_conv_criteria = (2 <= opt->status_verbosity);
  const int           print_summary = opt->print_summary;
  const int           iter_start = opt->iter_start;
  const int           iter_max = opt->iter_max;
  int                 iter;
  const double        rtol = opt->rtol;
  int                 neg_gradient_updated;
  int                 step_search_success;
  rhea_newton_step_t  step;
  rhea_newton_status_t  status;

  if (iter_start == 0) {
    RHEA_GLOBAL_PRODUCTIONF (
        "Into %s (max iter %i, rtol %.1e, nonzero guess %i)\n",
        this_fn_name, iter_max, rtol, opt->nonzero_initial_guess);
  }
  else {
    RHEA_GLOBAL_PRODUCTIONF (
        "Into %s (iter start %i, iter max %i, rtol %.1e, nonzero guess %i)\n",
        this_fn_name, iter_start, iter_max, rtol, opt->nonzero_initial_guess);
  }

  /*
   * Initialize
   */
  {
    /* initialize data */
    if (nl_problem->data_init != NULL) {
      nl_problem->data_init ((opt->nonzero_initial_guess ? solution : NULL),
                             nl_problem->data);
    }

    /* create step */
    rhea_newton_step_init (&step, nl_problem->step_vec);

    /* create status */
    rhea_newton_status_init (&status, nl_problem->conv_criterion,
                             nl_problem->grad_norm_multi_components);

    /* initialize solution vector and status */
    if (opt->nonzero_initial_guess) { /* if nonzero initial guess */
      rhea_newton_status_compute_curr (&status, &neg_gradient_updated,
                                       solution, nl_problem,
                                       print_all_conv_criteria);
    }
    else { /* if zero initial guess */
      ymir_vec_set_zero (solution);
      rhea_newton_status_compute_curr (&status, &neg_gradient_updated,
                                       NULL, nl_problem,
                                       print_all_conv_criteria);
    }
    rhea_newton_status_copy_curr_to_init (&status);
  }

  /*
   * Iterations Loop
   */

  for (iter = iter_start; iter <= iter_max; iter++) { /* BEGIN: Newton iter */

    /*
     * Pre-Step Output
     */
    {
      rhea_newton_status_print_curr (&status, print_all_conv_criteria, iter,
                                     this_fn_name);

      //TODO add more output
    }

    /*
     * Check Stopping Criteria
     */
    if (rhea_newton_status_get_reduction (&status) < rtol) { /* if converged */
      RHEA_GLOBAL_PRODUCTIONF (
          "%s: Nonlinear problem converged to rtol (%.3e)\n",
          this_fn_name, rtol);
      break;
    }
    if (iter_max == iter) { /* if max #iterations reached */
      RHEA_GLOBAL_PRODUCTIONF (
          "%s: Maximum number of nonlinear iterations reached (%i)\n",
          this_fn_name, iter_max);
      break;
    }

    /*
     * Newton Step
     */
    {
      step.iter = iter;

      /* update Hessian operator */
      if (iter != iter_start) { /* if this is not the first Newton step */
        if (nl_problem->update_hessian != NULL) {
          nl_problem->update_hessian (solution, nl_problem->data);
        }
        if (nl_problem->check_hessian) {
          rhea_newton_check_hessian (solution, nl_problem);
        }
      }

      /* calculate the accuracy for computing the step, i.e., for solving the
       * linearized system */
      rhea_newton_set_accuracy (
          /* out: */ &step,
          /* in:  */ &status, opt);

      /* compute right-hand side for the step solve */
      if (!neg_gradient_updated) { /* if (neg.) gradient was not updated */
        nl_problem->compute_neg_gradient (nl_problem->neg_gradient_vec,
                                          solution, nl_problem->data);
        if (nl_problem->check_gradient) {
          rhea_newton_check_gradient (solution, nl_problem);
        }
      }

      /* solve the linearized system to get an inexact Newton step */
      rhea_newton_compute_step (
          /* out: */ &step,
          /* in:  */ nl_problem, opt);

      /* perform line search to get the step length (updates the solution and
       * the nonlinear operator) */
      step_search_success = rhea_newton_search_step_length (
          /* out: */ solution, &step, &status, &neg_gradient_updated,
          /* in:  */ nl_problem, opt);
    }

    /*
     * Post-Step Output
     */
    {
      RHEA_GLOBAL_INFOF (
          "Newton iter %i -- %s: Linear solver: num iter %i, "
          "residual reduction %.3e, prescribed rtol %.3e, "
          "convergence %.3e (= rtol^(1/#iter))\n",
          iter, this_fn_name, step.lin_iter_count,
          step.lin_res_norm_reduction, step.lin_res_norm_rtol,
          step.lin_convergence);

      RHEA_GLOBAL_INFOF (
          "Newton iter %i -- %s: Step length %g\n",
          iter, this_fn_name, step.length);

      //TODO move summary to the end
      if (print_summary) {
        RHEA_GLOBAL_PRODUCTIONF (
            "%s: Newton step summary: [%2d, %.6f, %4d, %.3e, %.3e, %.3e]\n",
            this_fn_name, iter, step.length,
            step.lin_iter_count, step.lin_res_norm_reduction,
            step.lin_res_norm_rtol, step.lin_convergence);
      }
    }

    /*
     * Abort Check
     */
    if (!step_search_success && opt->abort_failed_step_search) {
      RHEA_GLOBAL_PRODUCTIONF ("%s: Abort due to failed step length search\n",
                               this_fn_name);
      break;
    }

    /*
     * Setup Next Step
     */
    {
      //TODO
    }

  } /* END: Newton iter */

  /*
   * Finalize
   */
  {
    /* destroy */
    rhea_newton_status_clear (&status);

    /* clear data */
    if (nl_problem->data_clear != NULL) {
      nl_problem->data_clear (nl_problem->data);
    }

  }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}
