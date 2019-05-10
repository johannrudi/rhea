#include <rhea_newton.h>
#include <rhea_newton_check.h>
#include <rhea_base.h>

/* Newton status */
typedef struct rhea_newton_status rhea_newton_status_t;
typedef struct rhea_newton_status_summary rhea_newton_status_summary_t;

/* Newton step */
typedef struct rhea_newton_step rhea_newton_step_t;

/* Nonlinear problem */
struct rhea_newton_problem
{
  /* vectors (not owned) */
  ymir_vec_t         *neg_gradient_vec;
  ymir_vec_t         *step_vec;

  /* callback functions for Newton algorithm */
  rhea_newton_evaluate_objective_fn_t         evaluate_objective;
  rhea_newton_compute_negative_gradient_fn_t  compute_neg_gradient;
  rhea_newton_compute_norm_of_gradient_fn_t   compute_gradient_norm;
  int                 obj_multi_components;
  int                 grad_norm_multi_components;

  rhea_newton_apply_hessian_fn_t        apply_hessian;
  rhea_newton_solve_hessian_system_fn_t solve_hessian_sys;

  rhea_newton_update_operator_fn_t        update_operator;
  rhea_newton_update_hessian_fn_t         update_hessian;
  rhea_newton_modify_hessian_system_fn_t  modify_hessian_system;

  /* data and related callback functions */
  void                        *data;
  rhea_newton_data_init_fn_t   data_init;
  rhea_newton_data_clear_fn_t  data_clear;

  /* setup post Newton step */
  rhea_newton_setup_poststep_fn_t setup_poststep;

  /* user output */
  rhea_newton_output_prestep_fn_t output_prestep;

  /* status (not owned) */
  rhea_newton_status_t *status;
  int                 num_iterations;
  double              residual_reduction;

  /* options */
  int                 check_gradient;
  int                 check_hessian;
  sc_MPI_Comm         mpicomm;
};

/**
 * Calculates the reduction of an end value relative to a start value.  Avoids
 * ill-defined operations like division by zero.  Use, e.g., for checking the
 * reduction of a residual or objective functional.
 */
static double
rhea_newton_calculate_reduction (const double start_value,
                                 const double end_value)
{
  /* check input */
  RHEA_ASSERT (isfinite (start_value));
  RHEA_ASSERT (0.0 <= start_value);
  RHEA_ASSERT (isfinite (end_value));
  RHEA_ASSERT (0.0 <= end_value);

  if (end_value < DBL_MIN) { /* if end value is small */
    /* assume max reduction and set value below machine precision */
    return pow (SC_EPS, 2.0);
  }
  else if (start_value < SC_EPS * end_value) { /* if start value is small */
    /* assume no reduction */
    return 1.0;
  }
  else { /* if division is well-defined */
    /* compute reduction */
    return end_value / start_value;
  }
}

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
#define RHEA_NEWTON_DEFAULT_LIN_MONITOR_REDUCTION (0)
#define RHEA_NEWTON_DEFAULT_STEP_SEARCH_ITER_MAX (12)
#define RHEA_NEWTON_DEFAULT_STEP_LENGTH_MIN (1.0e-5)
#define RHEA_NEWTON_DEFAULT_STEP_LENGTH_MAX (1.0)
#define RHEA_NEWTON_DEFAULT_STEP_REDUCTION (0.5)
#define RHEA_NEWTON_DEFAULT_STEP_DESCEND_CONDITION_RELAXATION (1.0e-4)
#define RHEA_NEWTON_DEFAULT_PRINT_SUMMARY (0)

#define RHEA_NEWTON_DEFAULT_NONZERO_INITIAL_GUESS (0)
#define RHEA_NEWTON_DEFAULT_ABORT_FAILED_STEP_SEARCH (0)
#define RHEA_NEWTON_DEFAULT_STATUS_VERBOSITY (0)

/* global options */
rhea_newton_options_t rhea_newton_options;

void
rhea_newton_add_options (rhea_newton_options_t *newton_options,
                         ymir_options_t *opt_sup)
{
  const char         *opt_prefix = "Newton";
  ymir_options_t     *opt = ymir_options_new ();
  rhea_newton_options_t  *newton_opt;

  /* set options storage */
  if (newton_options != NULL) {
    /* choose provided options */
    newton_opt = newton_options;
  }
  else {
    /* choose global options */
    newton_opt = &rhea_newton_options;
  }

  /* initialize options */
  rhea_newton_options_set_defaults (newton_opt);

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_I, "iter-start", '\0',
    &(newton_opt->iter_start), RHEA_NEWTON_DEFAULT_ITER_START,
    "Start at this iteration number",
  YMIR_OPTIONS_I, "iter-max", '\0',
    &(newton_opt->iter_max), RHEA_NEWTON_DEFAULT_ITER_MAX,
    "Maximum number of iterations",
  YMIR_OPTIONS_D, "rtol", '\0',
    &(newton_opt->rtol), RHEA_NEWTON_DEFAULT_RTOL,
    "Relative tolerance",

  YMIR_OPTIONS_I, "lin-iter-max", '\0',
    &(newton_opt->lin_iter_max), RHEA_NEWTON_DEFAULT_LIN_ITER_MAX,
    "Linear sub-solver: Maximum number of iterations",
  YMIR_OPTIONS_I, "lin-rtol-init-n-iter", '\0',
    &(newton_opt->lin_rtol_init_n_iter),
    RHEA_NEWTON_DEFAULT_LIN_RTOL_INIT_N_ITER,
    "Linear sub-solver: Number of iterations using initial relative tolerance",
  YMIR_OPTIONS_D, "lin-rtol-init", '\0',
    &(newton_opt->lin_rtol_init), RHEA_NEWTON_DEFAULT_LIN_RTOL_INIT,
    "Linear sub-solver: Initial relative tolerance",
  YMIR_OPTIONS_D, "lin-rtol-adaptive-exponent", '\0',
    &(newton_opt->lin_rtol_adaptive_exponent),
    RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_EXPONENT,
    "Linear sub-solver: Exponent for adaptive relative tolerance",
  YMIR_OPTIONS_D, "lin-rtol-adaptive-max", '\0',
    &(newton_opt->lin_rtol_adaptive_max),
    RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_MAX,
    "Linear sub-solver: Max relative tolerance",
  YMIR_OPTIONS_I, "lin-rtol-adaptive-min-active", '\0',
    &(newton_opt->lin_rtol_adaptive_min_active),
    RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_MIN_ACTIVE,
    "Linear sub-solver: Use threshold for min relative tolerance",
  YMIR_OPTIONS_D, "lin-rtol-adaptive-min-threshold", '\0',
    &(newton_opt->lin_rtol_adaptive_min_threshold),
    RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_MIN_THRESHOLD,
    "Linear sub-solver: Min relative tolerance",
  YMIR_OPTIONS_I, "lin-rtol-adaptive-progressive-n-iter", '\0',
    &(newton_opt->lin_rtol_adaptive_progressive_n_iter),
    RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_PROGRESSIVE_N_ITER,
    "Linear sub-solver: #iter (~ variance) for progressive tightening of rtol",
  YMIR_OPTIONS_B, "lin-monitor-reduction", '\0',
    &(newton_opt->lin_monitor_reduction),
    RHEA_NEWTON_DEFAULT_LIN_MONITOR_REDUCTION,
    "Linear sub-solver: Monitor residual reduction of each linearized solve",

  YMIR_OPTIONS_I, "step-search-iter-max", '\0',
    &(newton_opt->step_search_iter_max),
    RHEA_NEWTON_DEFAULT_STEP_SEARCH_ITER_MAX,
    "Line search: Maximum number of step search iterations",
  YMIR_OPTIONS_D, "step-length-min", '\0',
    &(newton_opt->step_length_min), RHEA_NEWTON_DEFAULT_STEP_LENGTH_MIN,
    "Line search: Minimum allowed step length",
  YMIR_OPTIONS_D, "step-length-max", '\0',
    &(newton_opt->step_length_max), RHEA_NEWTON_DEFAULT_STEP_LENGTH_MAX,
    "Line search: Maximum step length at start of search (usually 1)",
  YMIR_OPTIONS_D, "step-reduction", '\0',
    &(newton_opt->step_reduction), RHEA_NEWTON_DEFAULT_STEP_REDUCTION,
    "Line search: Reduction factor for step search",
  YMIR_OPTIONS_D, "step-descend-condition-relaxation", '\0',
    &(newton_opt->step_descend_condition_relaxation),
    RHEA_NEWTON_DEFAULT_STEP_DESCEND_CONDITION_RELAXATION,
    "Line search: Relaxation factor for the descend condition",

  YMIR_OPTIONS_B, "print-summary", '\0',
    &(newton_opt->print_summary), RHEA_NEWTON_DEFAULT_PRINT_SUMMARY,
    "Print summary of Newton iterations",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);
}

void
rhea_newton_get_options (rhea_newton_options_t *opt)
{
  rhea_newton_options_set_defaults (opt);

  opt->iter_start = rhea_newton_options.iter_start;
  opt->iter_max   = rhea_newton_options.iter_max;
  opt->rtol       = rhea_newton_options.rtol;

  opt->lin_iter_max                 = rhea_newton_options.lin_iter_max;
  opt->lin_rtol_init_n_iter         = rhea_newton_options.lin_rtol_init_n_iter;
  opt->lin_rtol_init                = rhea_newton_options.lin_rtol_init;
  opt->lin_rtol_adaptive_exponent =
    rhea_newton_options.lin_rtol_adaptive_exponent;
  opt->lin_rtol_adaptive_max =
    rhea_newton_options.lin_rtol_adaptive_max;
  opt->lin_rtol_adaptive_min_active =
    rhea_newton_options.lin_rtol_adaptive_min_active;
  opt->lin_rtol_adaptive_min_threshold =
    rhea_newton_options.lin_rtol_adaptive_min_threshold;
  opt->lin_rtol_adaptive_progressive_n_iter =
    rhea_newton_options.lin_rtol_adaptive_progressive_n_iter;

  opt->step_search_iter_max = rhea_newton_options.step_search_iter_max;
  opt->step_length_min      = rhea_newton_options.step_length_min;
  opt->step_length_max      = rhea_newton_options.step_length_max;
  opt->step_reduction       = rhea_newton_options.step_reduction;
  opt->step_descend_condition_relaxation =
    rhea_newton_options.step_descend_condition_relaxation;

  opt->print_summary = rhea_newton_options.print_summary;
}

void
rhea_newton_options_set_defaults (rhea_newton_options_t *opt)
{
  opt->nonzero_initial_guess    = RHEA_NEWTON_DEFAULT_NONZERO_INITIAL_GUESS;
  opt->abort_failed_step_search = RHEA_NEWTON_DEFAULT_ABORT_FAILED_STEP_SEARCH;

  opt->iter_start = RHEA_NEWTON_DEFAULT_ITER_START;
  opt->iter_max   = RHEA_NEWTON_DEFAULT_ITER_MAX;
  opt->rtol       = RHEA_NEWTON_DEFAULT_RTOL;

  opt->lin_iter_max           = RHEA_NEWTON_DEFAULT_LIN_ITER_MAX;
  opt->lin_rtol_init_n_iter   = RHEA_NEWTON_DEFAULT_LIN_RTOL_INIT_N_ITER;
  opt->lin_rtol_init          = RHEA_NEWTON_DEFAULT_LIN_RTOL_INIT;
  opt->lin_rtol_adaptive_exponent =
    RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_EXPONENT;
  opt->lin_rtol_adaptive_max  = RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_MAX;
  opt->lin_rtol_adaptive_min_active =
    RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_MIN_ACTIVE;
  opt->lin_rtol_adaptive_min_threshold =
    RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_MIN_THRESHOLD;
  opt->lin_rtol_adaptive_progressive_n_iter =
    RHEA_NEWTON_DEFAULT_LIN_RTOL_ADAPTIVE_PROGRESSIVE_N_ITER;

  opt->step_search_iter_max = RHEA_NEWTON_DEFAULT_STEP_SEARCH_ITER_MAX;
  opt->step_length_min      = RHEA_NEWTON_DEFAULT_STEP_LENGTH_MIN;
  opt->step_length_max      = RHEA_NEWTON_DEFAULT_STEP_LENGTH_MAX;
  opt->step_reduction       = RHEA_NEWTON_DEFAULT_STEP_REDUCTION;
  opt->step_descend_condition_relaxation =
    RHEA_NEWTON_DEFAULT_STEP_DESCEND_CONDITION_RELAXATION;

  opt->status_verbosity = RHEA_NEWTON_DEFAULT_STATUS_VERBOSITY;
  opt->print_summary    = RHEA_NEWTON_DEFAULT_PRINT_SUMMARY;
}

/******************************************************************************
 * Nonlinear Problem
 *****************************************************************************/

rhea_newton_problem_t *
rhea_newton_problem_new (
              rhea_newton_compute_negative_gradient_fn_t compute_neg_gradient,
              rhea_newton_compute_norm_of_gradient_fn_t compute_gradient_norm,
              const int grad_norm_multi_components,
              rhea_newton_solve_hessian_system_fn_t solve_hessian_sys)
{
  rhea_newton_problem_t  *nl_problem = RHEA_ALLOC (rhea_newton_problem_t, 1);

  /* check input */
  RHEA_ASSERT (compute_neg_gradient != NULL);
  RHEA_ASSERT (compute_gradient_norm != NULL);
  RHEA_ASSERT (solve_hessian_sys != NULL);

  nl_problem->neg_gradient_vec = NULL;
  nl_problem->step_vec = NULL;
  nl_problem->compute_neg_gradient = compute_neg_gradient;
  nl_problem->compute_gradient_norm = compute_gradient_norm;
  nl_problem->grad_norm_multi_components = grad_norm_multi_components;
  nl_problem->solve_hessian_sys = solve_hessian_sys;

  rhea_newton_problem_set_data_fn (NULL, NULL, NULL, nl_problem);
  rhea_newton_problem_set_evaluate_objective_fn (NULL, 0, nl_problem);
  rhea_newton_problem_set_apply_hessian_fn (NULL, nl_problem);
  rhea_newton_problem_set_update_fn (NULL, NULL, NULL, nl_problem);
  rhea_newton_problem_set_setup_poststep_fn (NULL, nl_problem);
  rhea_newton_problem_set_output_fn (NULL, nl_problem);
  rhea_newton_problem_set_checks (0, 0, nl_problem);
  rhea_newton_problem_set_mpicomm (MPI_COMM_NULL, nl_problem);

  nl_problem->status = NULL;
  nl_problem->num_iterations = -1;
  nl_problem->residual_reduction = NAN;

  return nl_problem;
}

void
rhea_newton_problem_destroy (rhea_newton_problem_t *nl_problem)
{
  RHEA_FREE (nl_problem);
}

static int
rhea_newton_problem_is_valid (rhea_newton_problem_t *nl_problem)
{
  return (
    nl_problem->neg_gradient_vec != NULL &&
    nl_problem->step_vec != NULL &&
    nl_problem->compute_neg_gradient != NULL &&
    nl_problem->solve_hessian_sys != NULL
  );
}

void
rhea_newton_problem_set_vectors (
              ymir_vec_t *neg_gradient_vec,
              ymir_vec_t *step_vec,
              rhea_newton_problem_t *nl_problem)
{
  nl_problem->neg_gradient_vec = neg_gradient_vec;
  nl_problem->step_vec = step_vec;
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
rhea_newton_problem_set_evaluate_objective_fn (
              rhea_newton_evaluate_objective_fn_t evaluate_objective,
              const int obj_multi_components,
              rhea_newton_problem_t *nl_problem)
{
  nl_problem->evaluate_objective = evaluate_objective;
  nl_problem->obj_multi_components = obj_multi_components;
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
              rhea_newton_modify_hessian_system_fn_t modify_hessian_system,
              rhea_newton_problem_t *nl_problem)
{
  nl_problem->update_operator = update_operator;
  nl_problem->update_hessian = update_hessian;
  nl_problem->modify_hessian_system = modify_hessian_system;
}

void
rhea_newton_problem_set_setup_poststep_fn (
              rhea_newton_setup_poststep_fn_t setup_poststep,
              rhea_newton_problem_t *nl_problem)
{
  nl_problem->setup_poststep = setup_poststep;
}

void
rhea_newton_problem_set_output_fn (
              rhea_newton_output_prestep_fn_t output_prestep,
              rhea_newton_problem_t *nl_problem)
{
  nl_problem->output_prestep = output_prestep;
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

sc_MPI_Comm
rhea_newton_problem_get_mpicomm (rhea_newton_problem_t *nl_problem)
{
  return nl_problem->mpicomm;
}

void
rhea_newton_problem_set_mpicomm (sc_MPI_Comm mpicomm,
                                 rhea_newton_problem_t *nl_problem)
{
  nl_problem->mpicomm = mpicomm;
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

static int
rhea_newton_problem_data_init_exists (rhea_newton_problem_t *nl_problem)
{
  return (NULL != nl_problem->data_init);
}

static void
rhea_newton_problem_data_init (ymir_vec_t *solution,
                               rhea_newton_problem_t *nl_problem)
{
  if (rhea_newton_problem_data_init_exists (nl_problem)) {
    nl_problem->data_init (solution, nl_problem->data);
  }
}

static int
rhea_newton_problem_data_clear_exists (rhea_newton_problem_t *nl_problem)
{
  return (NULL != nl_problem->data_clear);
}

static void
rhea_newton_problem_data_clear (rhea_newton_problem_t *nl_problem)
{
  if (rhea_newton_problem_data_clear_exists (nl_problem)) {
    nl_problem->data_clear (nl_problem->data);
  }
}

int
rhea_newton_problem_evaluate_objective_exists (
                                            rhea_newton_problem_t *nl_problem)
{
  return (NULL != nl_problem->evaluate_objective);
}

double
rhea_newton_problem_evaluate_objective (ymir_vec_t *solution,
                                        rhea_newton_problem_t *nl_problem,
                                        double *obj_comp)
{
  RHEA_ASSERT (rhea_newton_problem_evaluate_objective_exists (nl_problem));
  return nl_problem->evaluate_objective (solution, nl_problem->data, obj_comp);
}

int
rhea_newton_problem_compute_gradient_norm_exists (
                                            rhea_newton_problem_t *nl_problem)
{
  return (NULL != nl_problem->compute_gradient_norm);
}

double
rhea_newton_problem_compute_gradient_norm (ymir_vec_t *neg_gradient,
                                           rhea_newton_problem_t *nl_problem,
                                           double *grad_norm_comp)
{
  RHEA_ASSERT (rhea_newton_problem_compute_gradient_norm_exists (nl_problem));
  return nl_problem->compute_gradient_norm (neg_gradient, nl_problem->data,
                                            grad_norm_comp);
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
  if (rhea_newton_problem_update_operator_exists (nl_problem)) {
    nl_problem->update_operator (solution, nl_problem->data);
  }
}

int
rhea_newton_problem_update_hessian_exists (rhea_newton_problem_t *nl_problem)
{
  return (NULL != nl_problem->update_hessian);
}

void
rhea_newton_problem_update_hessian (ymir_vec_t *solution,
                                    ymir_vec_t *step_vec,
                                    const double step_length,
                                    rhea_newton_problem_t *nl_problem)
{
  if (rhea_newton_problem_update_hessian_exists (nl_problem)) {
    nl_problem->update_hessian (solution, step_vec, step_length,
                                nl_problem->data);
  }
}

int
rhea_newton_problem_modify_hessian_system_exists (
                                            rhea_newton_problem_t *nl_problem)
{
  return (NULL != nl_problem->modify_hessian_system);
}

void
rhea_newton_problem_modify_hessian_system (ymir_vec_t *neg_gradient,
                                           ymir_vec_t *solution,
                                           rhea_newton_problem_t *nl_problem)
{
  if (rhea_newton_problem_modify_hessian_system_exists (nl_problem)) {
    nl_problem->modify_hessian_system (neg_gradient, solution,
                                       nl_problem->data);
  }
}

int
rhea_newton_problem_setup_poststep_exists (rhea_newton_problem_t *nl_problem)
{
  return (NULL != nl_problem->setup_poststep);
}

int
rhea_newton_problem_setup_poststep (ymir_vec_t **solution,
                                    const int iter,
                                    rhea_newton_problem_t *nl_problem)
{
  if (rhea_newton_problem_setup_poststep_exists (nl_problem)) {
    RHEA_ASSERT (solution != NULL);
    return nl_problem->setup_poststep (solution, iter, nl_problem->data);
  }
  else {
    return 0;
  }
}

int
rhea_newton_problem_output_prestep_exists (rhea_newton_problem_t *nl_problem)
{
  return (NULL != nl_problem->output_prestep);
}

void
rhea_newton_problem_output_prestep (ymir_vec_t *solution,
                                    const int iter,
                                    rhea_newton_problem_t *nl_problem)
{
  if (rhea_newton_problem_output_prestep_exists (nl_problem)) {
    nl_problem->output_prestep (solution, iter, nl_problem->data);
  }
}

/******************************************************************************
 * Newton Step
 *****************************************************************************/

struct rhea_newton_step
{
  double              length;

  int                 search_success;
  int                 search_iter_count;

  int                 iter;

  double              lin_res_norm_rtol;
  double              lin_res_norm_reduction;
  int                 lin_iter_count;
  double              lin_convergence;
};

/**
 * Initializes a step object.
 */
static void
rhea_newton_step_init (rhea_newton_step_t *step)
{
  step->length = 0.0;

  step->search_success = 1;
  step->search_iter_count = 0;

  step->iter = -1;

  step->lin_res_norm_rtol = NAN;
  step->lin_res_norm_reduction = NAN;
  step->lin_iter_count = 0;
  step->lin_convergence = NAN;
}

/******************************************************************************
 * Newton Status
 *****************************************************************************/

struct rhea_newton_status
{
  /* objective value at (initial, previous, and current) iteration */
  double              obj_init;
  double              obj_prev;
  double              obj_curr;

  /* components of the objective */
  int                 obj_multi_components;
  double             *obj_curr_comp;

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
  double             *grad_norm_curr_comp;

  /* reduction of gradient norm in just one iteration (previous and current) */
  double              grad_norm_reduction_prev;
  double              grad_norm_reduction_curr;

  /* reduction of gradient norm over all iterations */
  double              grad_norm_reduction;
};

/**
 * Initializes a status object.
 */
static void
rhea_newton_status_init (rhea_newton_status_t *status,
                         const int obj_multi_components,
                         const int grad_norm_multi_components)
{
  int                 compid;

  /* init valus pertaining to objective values */
  status->obj_init = NAN;
  status->obj_prev = NAN;
  status->obj_curr = NAN;

  status->obj_multi_components = obj_multi_components;
  if (0 < status->obj_multi_components) {
    status->obj_curr_comp = RHEA_ALLOC (double, status->obj_multi_components);
    for (compid = 0; compid < status->obj_multi_components; compid++) {
      status->obj_curr_comp[compid] = NAN;
    }
  }
  else {
    status->obj_curr_comp = NULL;
  }

  status->obj_reduction_prev = 1.0;
  status->obj_reduction_curr = 1.0;
  status->obj_reduction = 1.0;

  /* init valus pertaining to gradient norm(s) */
  status->grad_norm_init = NAN;
  status->grad_norm_prev = NAN;
  status->grad_norm_curr = NAN;

  status->grad_norm_multi_components = grad_norm_multi_components;
  if (0 < status->grad_norm_multi_components) {
    status->grad_norm_curr_comp =
      RHEA_ALLOC (double, status->grad_norm_multi_components);
    for (compid = 0; compid < status->grad_norm_multi_components; compid++) {
      status->grad_norm_curr_comp[compid] = NAN;
    }
  }
  else {
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
  if (status->obj_curr_comp != NULL) {
    RHEA_FREE (status->obj_curr_comp);
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
                             const double *obj_comp,
                             const double grad_norm,
                             const double *grad_norm_comp)
{
  int                 compid;

  /* set values pertaining to objective */
  if (isfinite (obj)) {
    status->obj_curr = obj;
    if (0 < status->obj_multi_components && obj_comp != NULL) {
      for (compid = 0; compid < status->obj_multi_components; compid++) {
        status->obj_curr_comp[compid] = obj_comp[compid];
      }
    }
    if (isfinite (status->obj_prev)) {
      status->obj_reduction_curr = rhea_newton_calculate_reduction (
          status->obj_prev, status->obj_curr);
    }
    if (isfinite (status->obj_init)) {
      status->obj_reduction = rhea_newton_calculate_reduction (
          status->obj_init, status->obj_curr);
    }
  }

  /* set valus pertaining to gradient norm(s) */
  if (isfinite (grad_norm)) {
    status->grad_norm_curr = grad_norm;
    if (0 < status->grad_norm_multi_components && grad_norm_comp != NULL) {
      for (compid = 0; compid < status->grad_norm_multi_components; compid++) {
        status->grad_norm_curr_comp[compid] = grad_norm_comp[compid];
      }
    }
    if (isfinite (status->grad_norm_prev)) {
      status->grad_norm_reduction_curr = rhea_newton_calculate_reduction (
          status->grad_norm_prev, status->grad_norm_curr);
    }
    if (isfinite (status->grad_norm_init)) {
      status->grad_norm_reduction = rhea_newton_calculate_reduction (
          status->grad_norm_init, status->grad_norm_curr);
    }
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
}

/* enumerator for convergence critera */
typedef enum
{
  RHEA_NEWTON_CONV_CRITERION_OBJECTIVE,
  RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM,
  RHEA_NEWTON_CONV_CRITERION_ALL
}
rhea_newton_conv_criterion_t;

/**
 * Computes and sets the current status.
 */
static void
rhea_newton_status_compute_curr (
                            rhea_newton_status_t *status,
                            int *neg_gradient_updated,
                            const rhea_newton_conv_criterion_t conv_criterion,
                            ymir_vec_t *solution,
                            const int iter,
                            rhea_newton_problem_t *nl_problem)
{
  int                 compute_obj = 0;
  int                 compute_grad = 0;
  double              obj_val = NAN;
  double             *obj_comp = NULL;
  double              grad_norm = NAN;
  double             *grad_norm_comp = NULL;

  /* set flags, which values should be computed */
  switch (conv_criterion) {
  case RHEA_NEWTON_CONV_CRITERION_OBJECTIVE:
    compute_obj = rhea_newton_problem_evaluate_objective_exists (nl_problem);
    compute_grad = !compute_obj;
    break;
  case RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM:
    compute_obj = 0;
    compute_grad = 1;
    break;
  case RHEA_NEWTON_CONV_CRITERION_ALL:
    compute_obj = rhea_newton_problem_evaluate_objective_exists (nl_problem);
    compute_grad = 1;
    break;
  default: /* unknown criterion */
    RHEA_ABORT_NOT_REACHED ();
  }
  RHEA_ASSERT (compute_obj || compute_grad);

  /* evaluate objective functional */
  if (compute_obj) {
    if (0 < status->obj_multi_components) {
      obj_comp = RHEA_ALLOC (double, status->obj_multi_components);
    }
    obj_val = rhea_newton_problem_evaluate_objective (solution, nl_problem,
                                                      obj_comp);
    RHEA_ASSERT (isfinite (obj_val));
  }

  /* compute gradient norm */
  if (compute_grad) {
    ymir_vec_t         *neg_gradient = nl_problem->neg_gradient_vec;

    /* compute (negative) gradient */
    RHEA_ASSERT (nl_problem->neg_gradient_vec != NULL);
    rhea_newton_problem_compute_neg_gradient (neg_gradient, solution,
                                              nl_problem);
    if (neg_gradient_updated != NULL) {
      *neg_gradient_updated = 1;
    }
    if (nl_problem->check_gradient) {
      rhea_newton_check_gradient (solution, neg_gradient, iter, nl_problem);
    }

    /* compute norm of gradient */
    if (0 < status->grad_norm_multi_components) {
      grad_norm_comp = RHEA_ALLOC (double, status->grad_norm_multi_components);
    }
    grad_norm = rhea_newton_problem_compute_gradient_norm (neg_gradient,
                                                           nl_problem,
                                                           grad_norm_comp);
    RHEA_ASSERT (isfinite (grad_norm));
    RHEA_ASSERT (0.0 <= grad_norm);

    /* set objective if its evaluation function does not exist */
    if (!rhea_newton_problem_evaluate_objective_exists (nl_problem)) {
      obj_val = grad_norm;
    }
  }
  else if (neg_gradient_updated != NULL) {
    *neg_gradient_updated = 0;
  }

  /* set current status */
  rhea_newton_status_set_curr (status, obj_val, obj_comp, grad_norm,
                               grad_norm_comp);

  /* destroy */
  if (obj_comp != NULL) {
    RHEA_FREE (obj_comp);
  }
  if (grad_norm_comp != NULL) {
    RHEA_FREE (grad_norm_comp);
  }
}

/**
 * Returns previous value that determines convergence.
 */
static double
rhea_newton_status_get_conv_value_prev (
                            rhea_newton_status_t *status,
                            const rhea_newton_conv_criterion_t conv_criterion)
{
  double              val;

  switch (conv_criterion) {
  case RHEA_NEWTON_CONV_CRITERION_OBJECTIVE:
    val = status->obj_prev;
    break;
  case RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM:
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
rhea_newton_status_get_conv_value_curr (
                            rhea_newton_status_t *status,
                            const rhea_newton_conv_criterion_t conv_criterion)
{
  double              val;

  switch (conv_criterion) {
  case RHEA_NEWTON_CONV_CRITERION_OBJECTIVE:
    val = status->obj_curr;
    break;
  case RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM:
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
rhea_newton_status_get_reduction (
                            rhea_newton_status_t *status,
                            const rhea_newton_conv_criterion_t conv_criterion)
{
  double              reduction;

  switch (conv_criterion) {
  case RHEA_NEWTON_CONV_CRITERION_OBJECTIVE:
    reduction = status->obj_reduction;
    break;
  case RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM:
    reduction = status->grad_norm_reduction;
    break;
  default: /* unknown criterion */
    RHEA_ABORT_NOT_REACHED ();
    reduction = -1.0;
  }

  return reduction;
}

static double
rhea_newton_status_get_grad_reduction (rhea_newton_status_t *status)
{
  return rhea_newton_status_get_reduction (
      status, RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM);
}

/**
 * Returns previous reduction of the value that determines convergence.
 */
static double
rhea_newton_status_get_reduction_prev (
                            rhea_newton_status_t *status,
                            const rhea_newton_conv_criterion_t conv_criterion)
{
  double              reduction;

  switch (conv_criterion) {
  case RHEA_NEWTON_CONV_CRITERION_OBJECTIVE:
    reduction = status->obj_reduction_prev;
    break;
  case RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM:
    reduction = status->grad_norm_reduction_prev;
    break;
  default: /* unknown criterion */
    RHEA_ABORT_NOT_REACHED ();
    reduction = -1.0;
  }

  return reduction;
}

static double
rhea_newton_status_get_grad_reduction_prev (rhea_newton_status_t *status)
{
  return rhea_newton_status_get_reduction_prev (
      status, RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM);
}

/**
 * Returns current reduction of the value that determines convergence.
 */
static double
rhea_newton_status_get_reduction_curr (
                            rhea_newton_status_t *status,
                            const rhea_newton_conv_criterion_t conv_criterion)
{
  double              reduction;

  switch (conv_criterion) {
  case RHEA_NEWTON_CONV_CRITERION_OBJECTIVE:
    reduction = status->obj_reduction_curr;
    break;
  case RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM:
    reduction = status->grad_norm_reduction_curr;
    break;
  default: /* unknown criterion */
    RHEA_ABORT_NOT_REACHED ();
    reduction = -1.0;
  }

  return reduction;
}

static double
rhea_newton_status_get_grad_reduction_curr (rhea_newton_status_t *status)
{
  return rhea_newton_status_get_reduction_curr (
      status, RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM);
}

/**
 * Writes an array of double values into a string.
 */
static void
rhea_newton_status_vals_to_str (char string[BUFSIZ],
                                const int n_values,
                                const double *value)
{
  char               *pos = string;
  const char         *end = string + BUFSIZ;
  int                 k;

  /* check input */
  RHEA_ASSERT (0 < n_values);
  RHEA_ASSERT (value != NULL);

  /* initialize empty string */
  snprintf (string, BUFSIZ, "");

  /* write array of values to string */
  pos += snprintf (pos, end - pos, "%.6e", value[0]);
  for (k = 1; k < n_values; k++) {
    pos += snprintf (pos, end - pos, ", %.6e", value[k]);
  }
}

/**
 * Helper functions for `rhea_newton_status_print_curr (...)` below.
 */
static void
rhea_newton_status_print_curr_obj (const double value, const double reduction,
                                   const int iter, const char *func_name)
{
  char                func_tag[BUFSIZ];

  snprintf (func_tag, BUFSIZ, "%s_status", func_name);
  RHEA_GLOBAL_INFOF_FN_TAG (
      func_tag, "newton_iter=%i, objective_functional=%.15e, reduction=%.3e",
      iter, value, reduction);
}

static void
rhea_newton_status_print_curr_grad (const double value, const double reduction,
                                    const int iter, const char *func_name)
{
  char                func_tag[BUFSIZ];

  snprintf (func_tag, BUFSIZ, "%s_status", func_name);
  RHEA_GLOBAL_INFOF_FN_TAG (
      func_tag, "newton_iter=%i, gradient_norm=%.15e, reduction=%.3e",
      iter, value, reduction);
}

static void
rhea_newton_status_print_curr_res (const double value, const double reduction,
                                   const int iter, const char *func_name)
{
  char                func_tag[BUFSIZ];

  snprintf (func_tag, BUFSIZ, "%s_status", func_name);
  RHEA_GLOBAL_INFOF_FN_TAG (
      func_tag, "newton_iter=%i, residual_norm=%.15e, reduction=%.3e",
      iter, value, reduction);
}

static void
rhea_newton_status_print_curr_comp (const int n_components,
                                    const double *comp_value,
                                    const int iter, const char *func_name,
                                    const char *component_label)
{
  char                func_tag[BUFSIZ];
  char                comp_vals[BUFSIZ];

  snprintf (func_tag, BUFSIZ, "%s_status", func_name);
  rhea_newton_status_vals_to_str (comp_vals, n_components, comp_value);
  RHEA_GLOBAL_INFOF_FN_TAG (
      func_tag, "newton_iter=%i, %s=[%s]", iter, component_label, comp_vals);
}

/**
 * Prints current Newton status.
 */
static void
rhea_newton_status_print_curr (rhea_newton_status_t *status,
                               const int evaluate_objective_exists,
                               const int iter,
                               const char *func_name)
{
  if (evaluate_objective_exists) {
    rhea_newton_status_print_curr_obj (status->obj_curr, status->obj_reduction,
                                       iter, func_name);
    if (0 < status->obj_multi_components) { /* if print components */
      rhea_newton_status_print_curr_comp (status->obj_multi_components,
                                          status->obj_curr_comp,
                                          iter, func_name,
                                          "objective_components");
    }
    rhea_newton_status_print_curr_grad (status->grad_norm_curr,
                                        status->grad_norm_reduction,
                                        iter, func_name);
    if (0 < status->grad_norm_multi_components) { /* if print components */
      rhea_newton_status_print_curr_comp (status->grad_norm_multi_components,
                                          status->grad_norm_curr_comp,
                                          iter, func_name,
                                          "gradient_norm_components");
    }
  }
  else {
    rhea_newton_status_print_curr_res (status->grad_norm_curr,
                                       status->grad_norm_reduction,
                                       iter, func_name);
    if (0 < status->grad_norm_multi_components) { /* if print components */
      rhea_newton_status_print_curr_comp (status->grad_norm_multi_components,
                                          status->grad_norm_curr_comp,
                                          iter, func_name,
                                          "residual_norm_components");
    }
  }
}

/* status summary */
struct rhea_newton_status_summary
{
  char              **iteration_stats;
  int                 iteration_length;
  int                 iteration_current_idx;
  double              total_obj_reduction;
  double              total_grad_norm_reduction;
  int                 total_step_search_success;
  int                 total_step_search_iter_count;
  int                 total_step_below_unit_count;
  double              total_step_length_avg;
  int                 total_lin_iter_count;
};

/**
 * Creates an object to store statistics of a Newton solve.
 */
static rhea_newton_status_summary_t *
rhea_newton_status_summary_new (const int max_num_iterations)
{
  rhea_newton_status_summary_t *summary;
  int                 k;

  /* create summary */
  summary = RHEA_ALLOC (rhea_newton_status_summary_t, 1);

  /* create storage for iterations */
  summary->iteration_stats = RHEA_ALLOC (char *, max_num_iterations);
  for (k = 0; k < max_num_iterations; k++) {
    summary->iteration_stats[k] = RHEA_ALLOC (char, BUFSIZ);
  }
  summary->iteration_length = max_num_iterations;
  summary->iteration_current_idx = 0;

  /* initialize totals */
  summary->total_obj_reduction = NAN;
  summary->total_grad_norm_reduction = NAN;
  summary->total_step_search_success = 1;
  summary->total_step_search_iter_count = 0;
  summary->total_step_below_unit_count = 0;
  summary->total_step_length_avg = 0.0;
  summary->total_lin_iter_count = 0;

  /* return summary */
  return summary;
}

/**
 * Destroys an object to store statistics of a Newton solve.
 */
static void
rhea_newton_status_summary_destroy (rhea_newton_status_summary_t *summary)
{
  int                 k;

  /* destroy storage for iterations */
  for (k = 0; k < summary->iteration_length; k++) {
    RHEA_FREE (summary->iteration_stats[k]);
  }
  RHEA_FREE (summary->iteration_stats);

  /* destroy summary */
  RHEA_FREE (summary);
}

/**
 * Adds current Newton status to the statistics summary.
 */
static void
rhea_newton_status_summary_write (rhea_newton_status_summary_t *summary,
                                  rhea_newton_status_t *status,
                                  const int evaluate_objective_exists,
                                  rhea_newton_step_t *step)
{
  double              obj, obj_reduction;
  double              grad_norm, grad_norm_reduction;
  char                obj_comp[BUFSIZ] = "";
  char                grad_norm_comp[BUFSIZ] = "";

  /* retrieve objective functional and gradient norm(s) */
  if (evaluate_objective_exists) {
    obj = status->obj_curr;
    obj_reduction = status->obj_reduction;
    if (0 < status->obj_multi_components) { /* if print components */
      rhea_newton_status_vals_to_str (
          obj_comp, status->obj_multi_components,
          status->obj_curr_comp);
    }
  }
  else {
    obj = NAN;
    obj_reduction = NAN;
  }
  grad_norm = status->grad_norm_curr;
  grad_norm_reduction = status->grad_norm_reduction;
  if (0 < status->grad_norm_multi_components) { /* if print components */
    rhea_newton_status_vals_to_str (
        grad_norm_comp, status->grad_norm_multi_components,
        status->grad_norm_curr_comp);
  }

  /* write status to string */
  snprintf (summary->iteration_stats[summary->iteration_current_idx], BUFSIZ,
            "%3d ; %.15e [%s], %.15e [%s] ; %i, %2i, %.15f ; %4d, %.3e (%.3e)",
            step->iter, obj, obj_comp, grad_norm, grad_norm_comp,
            step->search_success, step->search_iter_count, step->length,
            step->lin_iter_count,
            step->lin_res_norm_reduction, step->lin_res_norm_rtol);
  summary->iteration_current_idx++;

  /* update summary values */
  summary->total_obj_reduction = obj_reduction;
  summary->total_grad_norm_reduction = grad_norm_reduction;
  summary->total_step_search_success = summary->total_step_search_success &&
                                       step->search_success;
  summary->total_step_search_iter_count += SC_MAX (0, step->search_iter_count);
  summary->total_step_below_unit_count += (0 < step->search_iter_count);
  summary->total_step_length_avg += SC_MAX (0.0, step->length);
  summary->total_lin_iter_count += SC_MAX (0, step->lin_iter_count);
}

/**
 * Print statistics summary of a Newton solve.
 */
static void
rhea_newton_status_summary_print (rhea_newton_status_summary_t *summary)
{
  const int           n_iter = summary->iteration_current_idx - 1;
  int                 k;

  /* finalize summary values */
  summary->total_step_length_avg /= (double) n_iter;

  /* print summary */
  RHEA_GLOBAL_INFO ("========================================\n");
  RHEA_GLOBAL_INFO ("Newton solve summary\n");
  RHEA_GLOBAL_INFO ("----------------------------------------\n");
  for (k = 0; k < summary->iteration_current_idx; k++) {
    RHEA_GLOBAL_INFOF ("%s\n", summary->iteration_stats[k]);
  }
  RHEA_GLOBAL_INFO ("----------------------------------------\n");
  RHEA_GLOBAL_INFOF (
      "%3i ; %.15e, %.15e ; %i, %i, %i, %.15f ; %i\n",
      n_iter,
      summary->total_obj_reduction,
      summary->total_grad_norm_reduction,
      summary->total_step_search_success,
      summary->total_step_search_iter_count,
      summary->total_step_below_unit_count,
      summary->total_step_length_avg,
      summary->total_lin_iter_count);
  RHEA_GLOBAL_INFO ("========================================\n");
}

/******************************************************************************
 * Access to Newton Status
 *****************************************************************************/

double
rhea_newton_problem_get_reduction_curr (rhea_newton_problem_t *nl_problem)
{
  RHEA_ASSERT (nl_problem->status != NULL);
  return rhea_newton_status_get_reduction_curr (
      nl_problem->status, RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM);
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
  const int           iter = step->iter;

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "newton_iter=%i", iter);

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
    RHEA_GLOBAL_INFOF_FN_TAG (__func__, "type=\"initial\", lin_rtol=%.3e",
                              step->lin_res_norm_rtol);
    RHEA_GLOBAL_INFOF_FN_END (__func__, "newton_iter=%i", iter);
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
      SC_MIN (rhea_newton_status_get_grad_reduction_prev (status), 1.0);
    const double        reduction_curr =
      SC_MIN (rhea_newton_status_get_grad_reduction_curr (status), 1.0);

    double              lin_rtol_max;
    double              lin_rtol_min, lin_rtol_min_effective;
    double              lin_rtol;

    /* check input */
    RHEA_ASSERT (0.0 < opt->lin_rtol_adaptive_max);
    RHEA_ASSERT (opt->lin_rtol_adaptive_max < 1.0);
    RHEA_ASSERT (0.0 < rhea_newton_status_get_grad_reduction_prev (status));
    RHEA_ASSERT (0.0 < rhea_newton_status_get_grad_reduction_curr (status));

    /* set progressive reduction factor:
     *
     *   (exp(-(#iter / prog #iter)^2) + eps/rtol_max) /
     *                              (1 + eps/rtol_max)
     */
    if (prog_reduction_active) {
      const double        it2 = (double) iter * iter;
      const double        itp2 = (double)
                            opt->lin_rtol_adaptive_progressive_n_iter *
                            opt->lin_rtol_adaptive_progressive_n_iter;
      const double        eps = SC_1000_EPS;

      lin_rtol_max = opt->lin_rtol_adaptive_max;
      prog_reduction = (exp (-it2/itp2) + eps/lin_rtol_max) /
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

    RHEA_GLOBAL_INFOF_FN_TAG (
        __func__,
        "type=\"adaptive\" lin_rtol=%.3e, candidate_rtol=%.3e, "
        "max=%.3e (w/ progressive reduction=%.3e), "
        "min=%.3e (w/ threshold=%.3e)",
        step->lin_res_norm_rtol, lin_rtol,
        lin_rtol_max, prog_reduction,
        lin_rtol_min, lin_rtol_min_thresh);
    RHEA_GLOBAL_INFOF_FN_END (__func__, "newton_iter=%i", iter);
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
  const int           iter = step->iter;
  ymir_vec_t         *rhs = nl_problem->neg_gradient_vec;
  ymir_vec_t         *step_vec = nl_problem->step_vec;
  const int           nonzero_initial_guess = 0;
  const int           lin_iter_max = opt->lin_iter_max;
  int                 stop_reason, lin_iter_count = -1;
  const double        lin_res_norm_rtol = step->lin_res_norm_rtol;
  double              lin_res_norm_reduction;
  double              lin_conv;
  char                func_tag_stop[BUFSIZ];

  RHEA_GLOBAL_INFOF_FN_BEGIN (
      __func__, "newton_iter=%i, lin_iter_max=%i, lin_rtol=%.1e",
      iter, lin_iter_max, lin_res_norm_rtol);

  /* check input */
  RHEA_ASSERT (nl_problem->neg_gradient_vec != NULL);
  RHEA_ASSERT (nl_problem->step_vec != NULL);
  RHEA_ASSERT (0 <= step->iter);
  RHEA_ASSERT (0.0 < step->lin_res_norm_rtol && step->lin_res_norm_rtol < 1.0);

  /* run solver for the linearized system */
  stop_reason = nl_problem->solve_hessian_sys (
      step_vec, rhs, lin_iter_max, lin_res_norm_rtol,
      nonzero_initial_guess, nl_problem->data, &lin_iter_count);
  RHEA_ASSERT (0 <= lin_iter_count);

  /* calculate the residual reduction of the linearized solve */
  if (opt->lin_monitor_reduction &&
      rhea_newton_problem_apply_hessian_exists (nl_problem)) {
    ymir_vec_t         *lin_residual_vec;
    double              lin_res_norm_init;
    double              lin_res_norm_curr;

    /* compute the inital l^2-norm of the residual of the linearized system */
    lin_res_norm_init = ymir_vec_norm (rhs);

    /* compute the new l^2-norm of the residual of the linearized system */
    lin_residual_vec = ymir_vec_template (rhs);
    rhea_newton_problem_apply_hessian (lin_residual_vec, step_vec, nl_problem);
    ymir_vec_add (-1.0, rhs, lin_residual_vec);
    lin_res_norm_curr = ymir_vec_norm (lin_residual_vec);
    ymir_vec_destroy (lin_residual_vec);

    /* calculate linear residual reduction */
    lin_res_norm_reduction = rhea_newton_calculate_reduction (
        lin_res_norm_init, lin_res_norm_curr);

    /* calculate convergence of the linear solver `rtol^(1/#iter)` */
    lin_conv = exp ( log (lin_res_norm_reduction) / ((double) lin_iter_count) );
  }
  else { /* if actual linear reduction cannot be calculated */
    lin_res_norm_reduction = NAN;
    lin_conv = NAN;
  }

  /* set the #iterations used by the linear solver and the residual reduction */
  step->lin_iter_count = lin_iter_count;
  step->lin_res_norm_reduction = lin_res_norm_reduction;
  step->lin_convergence = lin_conv;

  snprintf (func_tag_stop, BUFSIZ, "%s_stop", __func__);
  RHEA_GLOBAL_INFOF_FN_TAG (
      func_tag_stop, "reason=%i, iterations=%i, "
      "residual_reduction=%.1e, lin_rtol=%.1e (prescribed)",
      stop_reason, lin_iter_count, lin_res_norm_reduction, lin_res_norm_rtol);
  RHEA_GLOBAL_INFOF_FN_END (__func__, "newton_iter=%i", iter);
}

static int
rhea_newton_search_step_length_check_descend (
                                          rhea_newton_status_t *status,
                                          rhea_newton_options_t *opt,
                                          const double lin_res_norm_reduction,
                                          const int print_condition,
                                          const char *func_name)
{
  const double        relax = opt->step_descend_condition_relaxation;
  double              reduction, descend_reduction;

  /* compute actual reduction */
  reduction = rhea_newton_status_get_reduction_curr (
      status, RHEA_NEWTON_CONV_CRITERION_OBJECTIVE);
  RHEA_ASSERT (0.0 <= reduction);

  /* compute required reduction */
  if (isfinite (lin_res_norm_reduction) &&
      0.0 < lin_res_norm_reduction && lin_res_norm_reduction <= 1.0) {
    descend_reduction = 1.0 - relax * (1.0 - lin_res_norm_reduction);
  }
  else { /* if linear reduction is out of bounds */
    descend_reduction = 1.0 - relax;
  }

  /* print */
  if (print_condition) {
    RHEA_GLOBAL_INFOF_FN_TAG (
        func_name, "descend_condition=\"reduction <= %.3e\"",
        descend_reduction);
  }

  /* check for descend */
  return (reduction <= descend_reduction);
}

/**
 * Performs line search for the step length (updates the nonlinear operator).
 */
static void
rhea_newton_search_step_length (ymir_vec_t *solution,
                                rhea_newton_step_t *step,
                                rhea_newton_status_t *status,
                                int *neg_gradient_updated,
                                rhea_newton_problem_t *nl_problem,
                                rhea_newton_options_t *opt)
{
  const int           iter = step->iter;
  ymir_vec_t         *step_vec = nl_problem->step_vec;
  ymir_vec_t         *solution_prev = ymir_vec_template (solution);
  const int           search_iter_start = 1;
  const int           search_iter_max = opt->step_search_iter_max;
  int                 k;
  int                 search_success = 0;
  double              conv_val_prev, conv_val_curr;
  const double        step_length_min = opt->step_length_min;
  const double        step_reduction = opt->step_reduction;
  char                func_tag_status[BUFSIZ];
  char                func_tag_stop[BUFSIZ];

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "newton_iter=%i", iter);

  /* check input */
  RHEA_ASSERT (nl_problem->step_vec != NULL);
  RHEA_ASSERT (0 <= step->iter);

  /* set function tags */
  snprintf (func_tag_status, BUFSIZ, "%s_status", __func__);
  snprintf (func_tag_stop, BUFSIZ, "%s_stop", __func__);

  /* copy previous solution and status */
  ymir_vec_copy (solution, solution_prev);
  rhea_newton_status_copy_curr_to_prev (status);

  /* initialize step length */
  step->length = opt->step_length_max;

  /* search for step length */
  for (k = search_iter_start; k <= search_iter_max; k++) {
    /* move solution in step direction */
    ymir_vec_add (step->length, step_vec, solution);

    /* update the nonlinear operator at new solution */
    rhea_newton_problem_update_operator (solution, nl_problem);

    /* update status */
    rhea_newton_status_compute_curr (status, neg_gradient_updated,
                                     RHEA_NEWTON_CONV_CRITERION_OBJECTIVE,
                                     solution, iter, nl_problem);

    /* check descend condition */
    search_success = rhea_newton_search_step_length_check_descend (
        status, opt, step->lin_res_norm_reduction,
        (search_iter_start == k) /* print only once */, __func__);

    conv_val_prev = rhea_newton_status_get_conv_value_prev (
        status, RHEA_NEWTON_CONV_CRITERION_OBJECTIVE);
    conv_val_curr = rhea_newton_status_get_conv_value_curr (
        status, RHEA_NEWTON_CONV_CRITERION_OBJECTIVE);
    RHEA_GLOBAL_INFOF_FN_TAG (
        func_tag_status,
        "trial#=%2i, step_length=%.3e, current=%.3e, previous=%.3e",
        k, step->length, conv_val_curr, conv_val_prev);

    /* terminate search if descend condition is satisfied */
    if (search_success) {
      RHEA_GLOBAL_INFOF_FN_TAG (
          func_tag_stop, "success=%i, step_length=%.3e",
          search_success, step->length);
      break;
    }
    /* terminate search if max number of search iterations reached */
    if (k == search_iter_max) {
      RHEA_GLOBAL_INFOF_FN_TAG (
          func_tag_stop, "success=%i, step_length=%.3e, "
          "reason=\"max number of step reductions reached\"",
          search_success, step->length);
      break;
    }
    /* terminate search if step length becomes too low */
    if ((step->length * step_reduction) < step_length_min) {
      RHEA_GLOBAL_INFOF_FN_TAG (
          func_tag_stop, "success=%i, step_length=%.3e, "
          "reason=\"min step length reached\"",
          search_success, step->length);
      break;
    }

    /* otherwise continue search loop with reduced step length and after
     * restoring the previous solution */
    step->length *= step_reduction;
    ymir_vec_copy (solution_prev, solution);
  }

  /* reverse updates of solution and nonlinear operator */
  if (!search_success && opt->abort_failed_step_search) {
    ymir_vec_copy (solution_prev, solution);
    rhea_newton_problem_update_operator (solution, nl_problem);
  }

  /* destroy */
  ymir_vec_destroy (solution_prev);

  RHEA_GLOBAL_INFOF_FN_END (__func__, "newton_iter=%i, success=%i",
                            iter, search_success);

  /* set output values */
  step->search_success = search_success;
  step->search_iter_count = k - search_iter_start;
}

int
rhea_newton_solve (ymir_vec_t **solution,
                   rhea_newton_problem_t *nl_problem,
                   rhea_newton_options_t *opt)
{
  const int           iter_start = opt->iter_start;
  const int           iter_max = opt->iter_max;
  int                 iter;
  const double        rtol = opt->rtol;
  int                 neg_gradient_updated;
  int                 solution_post_update = 0;
  int                 stop_reason = 0;
  const int           print_summary = opt->print_summary;
  rhea_newton_step_t            step;
  rhea_newton_status_t          status;
  rhea_newton_status_summary_t *summary = NULL;
  char                func_tag_lin_status[BUFSIZ];
  char                func_tag_step_status[BUFSIZ];
  char                func_tag_stop[BUFSIZ];

  RHEA_GLOBAL_PRODUCTIONF_FN_BEGIN (
      __func__, "iter_start=%i, iter_max=%i, rtol=%.1e, nonzero_init_guess=%i",
      iter_start, iter_max, rtol, opt->nonzero_initial_guess);

  /* check input */
  RHEA_ASSERT (rhea_newton_problem_is_valid (nl_problem));

  /*
   * Initialize
   */
  {
    /* initialize data */
    rhea_newton_problem_data_init (
        (opt->nonzero_initial_guess ? *solution : NULL), nl_problem);

    /* init step & status */
    rhea_newton_step_init (&step);
    rhea_newton_status_init (&status, nl_problem->obj_multi_components,
                             nl_problem->grad_norm_multi_components);
    nl_problem->status = &status;

    /* initialize solution vector and status */
    if (opt->nonzero_initial_guess) { /* if nonzero initial guess */
      rhea_newton_status_compute_curr (&status, &neg_gradient_updated,
                                       RHEA_NEWTON_CONV_CRITERION_ALL,
                                       *solution, iter_start, nl_problem);
    }
    else { /* if zero initial guess */
      ymir_vec_set_zero (*solution);
      rhea_newton_status_compute_curr (&status, &neg_gradient_updated,
                                       RHEA_NEWTON_CONV_CRITERION_ALL,
                                       NULL, iter_start, nl_problem);
    }
    rhea_newton_status_copy_curr_to_init (&status);

    /* create storage for summary */
    if (print_summary) {
      summary = rhea_newton_status_summary_new (iter_max - iter_start + 1);
    }
  }

  /* set function tags */
  snprintf (func_tag_lin_status, BUFSIZ, "%s_linear_solve_status", __func__);
  snprintf (func_tag_step_status, BUFSIZ, "%s_step_search_status", __func__);
  snprintf (func_tag_stop, BUFSIZ, "%s_stop", __func__);

  /*
   * Iterations Loop
   */
  for (iter = iter_start; iter <= iter_max; iter++) { /* BEGIN: Newton iter */
    step.iter = iter;

    /* check environment */
    RHEA_ASSERT (rhea_newton_problem_is_valid (nl_problem));

    /*
     * Pre-Step Output
     */
    {
      rhea_newton_status_print_curr (
          &status, rhea_newton_problem_evaluate_objective_exists (nl_problem),
          iter, __func__);
      if (print_summary) {
        rhea_newton_status_summary_write (
            summary, &status,
            rhea_newton_problem_evaluate_objective_exists (nl_problem), &step);
      }

      /* call user output function */
      rhea_newton_problem_output_prestep (*solution, iter, nl_problem);
    }

    /*
     * Check Stopping Criteria
     */

    /* stop if converged */
    if (rhea_newton_status_get_grad_reduction (&status) < rtol) {
      RHEA_GLOBAL_PRODUCTIONF_FN_TAG (
          func_tag_stop, "newton_iter=%i, reason=\"converged to rtol=%.3e\"",
          iter, rtol);
      stop_reason = 1;
      break;
    }

    /* stop if max #iterations reached */
    if (iter_max == iter) {
      RHEA_GLOBAL_PRODUCTIONF_FN_TAG (
          func_tag_stop,
          "newton_iter=%i, reason=\"reached maximum number of iterations=%i\"",
          iter, iter_max);
      stop_reason = 2;
      break;
    }

    /*
     * Newton Step
     */
    {
      /* compute right-hand side for the step solve */
      if (!neg_gradient_updated) { /* if (neg.) gradient not yet updated */
        rhea_newton_problem_compute_neg_gradient (nl_problem->neg_gradient_vec,
                                                  *solution, nl_problem);
        if (nl_problem->check_gradient) {
          rhea_newton_check_gradient (*solution, nl_problem->neg_gradient_vec,
                                      iter, nl_problem);
        }
      }

      /* update Hessian operator */
      if (iter == iter_start || solution_post_update) { /* if no step exists */
        rhea_newton_problem_update_hessian (*solution, NULL, NAN, nl_problem);
      }
      else { /* if step does exist */
        rhea_newton_problem_update_hessian (*solution, nl_problem->step_vec,
                                            step.length, nl_problem);
        if (nl_problem->check_hessian) {
          rhea_newton_check_hessian (*solution, nl_problem->neg_gradient_vec,
                                     iter, nl_problem);
        }
      }

      /* modify the right-hand side of the linear system, which originally is
       * the previously computed (negative) gradient  */
      if (iter == iter_start && !opt->nonzero_initial_guess) { /* if no sol. */
        rhea_newton_problem_modify_hessian_system (nl_problem->neg_gradient_vec,
                                                   NULL, nl_problem);
      }
      else { /* if solution vector exists and it is not zero */
        rhea_newton_problem_modify_hessian_system (nl_problem->neg_gradient_vec,
                                                   *solution, nl_problem);
      }

      /* calculate the accuracy for computing the step, i.e., for solving the
       * linearized system */
      rhea_newton_set_accuracy (
          /* out: */ &step,
          /* in:  */ &status, opt);

      /* solve the linearized system to get an inexact Newton step */
      rhea_newton_compute_step (
          /* out: */ &step,
          /* in:  */ nl_problem, opt);

      /* perform line search to get the step length (updates the solution and
       * the nonlinear operator) */
      rhea_newton_search_step_length (
          /* out: */ *solution, &step, &status, &neg_gradient_updated,
          /* in:  */ nl_problem, opt);
    }

    /*
     * Post-Step Output
     */
    {
      RHEA_GLOBAL_INFOF_FN_TAG (
          func_tag_lin_status,
          "newton_iter=%i, lin_iterations=%i, residual_reduction=%.3e, "
          "lin_rtol=%.3e (prescribed), convergence=%.3e (=rtol^(1/#iter))",
          iter, step.lin_iter_count, step.lin_res_norm_reduction,
          step.lin_res_norm_rtol, step.lin_convergence);

      RHEA_GLOBAL_INFOF_FN_TAG (
          func_tag_step_status, "newton_iter=%i, success=%i, step_length=%g",
          iter, step.search_success, step.length);
    }

    /*
     * Abort Check
     */
    if (!step.search_success && opt->abort_failed_step_search) {
      RHEA_GLOBAL_PRODUCTIONF_FN_TAG (
          func_tag_stop, "newton_iter=%i, reason=\"failed step length search\"",
          iter);
      stop_reason = 3;
      break;
    }

    /*
     * Setup Next Step
     */
    {
      /* run post-processing after completing one step */
      solution_post_update = rhea_newton_problem_setup_poststep (
          solution, iter, nl_problem);
      if (solution_post_update) {
        rhea_newton_status_compute_curr (&status, &neg_gradient_updated,
                                         RHEA_NEWTON_CONV_CRITERION_ALL,
                                         *solution, iter+1, nl_problem);
      }

      /* compute the gradient norm, if not already performed */
      if (!neg_gradient_updated) {
        rhea_newton_status_compute_curr (
                                      &status, &neg_gradient_updated,
                                      RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM,
                                      *solution, iter+1, nl_problem);
      }
    }

  } /* END: Newton iter */

  /*
   * Finalize
   */
  {
    /* store iterations count and residual reduction */
    nl_problem->num_iterations = iter;
    nl_problem->residual_reduction = status.grad_norm_reduction;

    /* destroy */
    rhea_newton_status_clear (&status);
    nl_problem->status = NULL;

    /* clear data */
    rhea_newton_problem_data_clear (nl_problem);

    /* print summary */
    if (print_summary) {
      rhea_newton_status_summary_print (summary);
      rhea_newton_status_summary_destroy (summary);
    }
  }

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);

  /* return stopping reason */
  return stop_reason;
}

int
rhea_newton_solve_get_num_iterations (rhea_newton_problem_t *nl_problem){
  return nl_problem->num_iterations;
}

double
rhea_newton_solve_get_residual_reduction (rhea_newton_problem_t *nl_problem)
{
  return nl_problem->residual_reduction;
}
