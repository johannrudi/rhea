/*
 */

#include <rhea_newton.h>
#include <rhea_base.h>

/* Nonlinear problem */
struct rhea_newton_problem
{
  ymir_vec_t         *neg_gradient_vec;
  ymir_vec_t         *step_vec;

  rhea_newton_problem_evaluate_objective_fn_t         evaluate_objective;
  rhea_newton_problem_compute_negative_gradient_fn_t  compute_neg_gradient;
  rhea_newton_problem_compute_norm_of_gradient_fn_t   compute_gradient_norm;
  int                 gradient_norm_multi_components;

  rhea_newton_problem_apply_hessian_fn_t         apply_hessian;
  rhea_newton_problem_solve_hessian_system_fn_t  solve_hessian_sys;

  rhea_newton_problem_update_operator_fn_t  update_operator;
  rhea_newton_problem_update_hessian_fn_t   update_hessian;

  void               *data;
  /*
  ymir_vec_t         *step_vec;
  ymir_vec_t         *residual_vec;

  rhea_newton_problem_apply_operator_fn_t   apply_operator;
  rhea_newton_problem_apply_jacobian_fn_t   apply_jacobian;
  rhea_newton_problem_update_operator_fn_t  update_operator;
  rhea_newton_problem_update_jacobian_fn_t  update_jacobian;

  rhea_newton_problem_compute_residual_fn_t       compute_residual;
  rhea_newton_problem_compute_residual_norm_fn_t  compute_residual_norm;
  int                 residual_norm_n_components;

  rhea_newton_problem_solve_linearized_fn_t  solve_linearized;

  void               *data;
  */
};

/* Newton options */
struct rhea_newton_options
{
  /* options for the nonlinear solver */
  int                 iter_start;
  int                 iter_max;
  double              rtol;

  /* options for the solver of the linearized system */
  int                 lin_iter_max;

  int                 lin_rtol_init_n_iter;
  double              lin_rtol_init;

  double              lin_rtol_adaptive_exponent;
  double              lin_rtol_adaptive_max;
  int                 lin_rtol_adaptive_min_active;
  double              lin_rtol_adaptive_min_threshold;
  int                 lin_rtol_adaptive_progressive_n_iter;

  /* options for line search */
  int                 step_search_iter_max;
  double              step_length_min;
  double              step_length_max;
  double              step_reduction;
  double              step_descend_condition_relaxation;

  int                 print_summary;
};

/* Newton residual */
typedef struct rhea_newton_residual
{
  ymir_vec_t         *vec;
  int                 norm_n_components;

  /* residual norm at an iteration (initial, previous, and current) */
  double              norm_init;
  double             *norm_init_comp;
  double              norm_prev;
  double             *norm_prev_comp;
  double              norm_curr;
  double             *norm_curr_comp;

  /* reduction of residual norm in just one iteration (previous and current) */
  double              norm_reduction_prev;
  double              norm_reduction_curr;

  /* reduction of residual norm over all iterations */
  double              norm_reduction;
}
rhea_newton_residual_t;

/* Newton step */
typedef struct rhea_newton_step
{
  ymir_vec_t         *rhs_vec;
  ymir_vec_t         *vec;
  double              length;

  int                 iter;

  double              lin_res_norm_rtol;
  double              lin_res_norm_reduction;
  int                 lin_iter_count;
  double              lin_convergence;
}
rhea_newton_step_t;

/* enumerator for convergence critera */
typedef enum
{
  RHEA_NEWTON_CONV_CRITERION_OBJECTIVE,
  RHEA_NEWTON_CONV_CRITERION_GRAD_NORM
}
rhea_newton_conv_criterion_t;

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

  /* number of components for a gradient with multiple compoents */
  int                 grad_norm_multi_components;

  /* gradient norm at (initial, previous, and current) iteration */
  double              grad_norm_init;
  double             *grad_norm_init_comp;
  double              grad_norm_prev;
  double             *grad_norm_prev_comp;
  double              grad_norm_curr;
  double             *grad_norm_curr_comp;

  /* reduction of gradient norm in just one iteration (previous and current) */
  double              grad_norm_reduction_prev;
  double              grad_norm_reduction_curr;

  /* reduction of gradient norm over all iterations */
  double              grad_norm_reduction;
}
rhea_newton_status_t;

/******************************************************************************
 * Nonlinear Problem
 *****************************************************************************/

rhea_newton_problem_t *
rhea_newton_problem_new (
    ymir_vec_t *neg_gradient_vec,
    ymir_vec_t *step_vec,
    rhea_newton_problem_evaluate_objective_fn_t evaluate_objective,
    rhea_newton_problem_compute_negative_gradient_fn_t compute_neg_gradient,
    rhea_newton_problem_compute_norm_of_gradient_fn_t compute_gradient_norm,
    rhea_newton_problem_apply_hessian_fn_t apply_hessian,
    rhea_newton_problem_solve_hessian_system_fn_t solve_hessian_sys,
    rhea_newton_problem_update_operator_fn_t update_operator,
    rhea_newton_problem_update_hessian_fn_t update_hessian,
    const int gradient_norm_multi_components,
    void *data)
{
  rhea_newton_problem_t  *nl_problem = RHEA_ALLOC (rhea_newton_problem_t, 1);

  nl_problem->neg_gradient_vec = neg_gradient_vec;
  nl_problem->step_vec = step_vec;

  nl_problem->evaluate_objective = evaluate_objective;
  nl_problem->compute_neg_gradient = compute_neg_gradient;
  nl_problem->compute_gradient_norm = compute_gradient_norm;
  nl_problem->apply_hessian = apply_hessian;
  nl_problem->solve_hessian_sys = solve_hessian_sys;
  nl_problem->update_operator = update_operator;
  nl_problem->update_hessian = update_hessian;

  nl_problem->gradient_norm_multi_components = gradient_norm_multi_components;

  nl_problem->data = data;

  return nl_problem;
}

void
rhea_newton_problem_destroy (rhea_newton_problem_t *nl_problem)
{
  RHEA_FREE (nl_problem);
}

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

/******************************************************************************
 * Newton Step
 *****************************************************************************/

/**
 * Initializes a step object.
 */
static void
rhea_newton_step_init (rhea_newton_step_t *step,
                       ymir_vec_t *rhs_vec,
                       ymir_vec_t *step_vec)
{
  step->rhs_vec = rhs_vec;
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
                         const int grad_norm_multi_components)
{
  /* init convergence criterion */
  status->conv_criterion = RHEA_NEWTON_CONV_CRITERION_OBJECTIVE; //TODO

  /* init valus pertaining to objective values */
  status->obj_init = -1.0;
  status->obj_prev = -1.0;
  status->obj_curr = -1.0;

  status->obj_reduction_prev = 1.0;
  status->obj_reduction_curr = 1.0;
  status->obj_reduction = 1.0;

  /* init valus pertaining to gradient norm(s) */
  status->grad_norm_multi_components = grad_norm_multi_components;

  status->grad_norm_init = -1.0;
  status->grad_norm_prev = -1.0;
  status->grad_norm_curr = -1.0;

  if (0 < status->grad_norm_multi_components) {
    const int           n_components = status->grad_norm_multi_components;
    int                 compid;

    status->grad_norm_init_comp = RHEA_ALLOC (double, n_components);
    status->grad_norm_prev_comp = RHEA_ALLOC (double, n_components);
    status->grad_norm_curr_comp = RHEA_ALLOC (double, n_components);
    for (compid = 0; compid < n_components; compid++) {
      status->grad_norm_init_comp[compid] = -1.0;
      status->grad_norm_prev_comp[compid] = -1.0;
      status->grad_norm_curr_comp[compid] = -1.0;
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
  if (0 < status->grad_norm_multi_components) {
    RHEA_FREE (status->grad_norm_init_comp);
    RHEA_FREE (status->grad_norm_prev_comp);
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

  if (0.0 < status->obj_prev) {
    status->obj_reduction_curr = status->obj_curr / status->obj_prev;
  }
  if (0.0 < status->obj_init) {
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

  if (0.0 < status->grad_norm_prev) {
    status->grad_norm_reduction_curr = status->grad_norm_curr /
                                       status->grad_norm_prev;
  }
  if (0.0 < status->grad_norm_init) {
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
 * Return overall reduction of the value that determines convergence.
 */
static double
rhea_newton_status_get_reduction (rhea_newton_status_t *status)
{
  double              reduction;

  switch (status->conv_criterion) {
  case RHEA_NEWTON_CONV_CRITERION_OBJECTIVE:
    reduction = status->obj_reduction;
  case RHEA_NEWTON_CONV_CRITERION_GRAD_NORM:
    reduction = status->grad_norm_reduction;
  default: /* unknown criterion */
    RHEA_ABORT_NOT_REACHED ();
    reduction = -1.0;
  }

  return reduction;
}

/**
 * Return previous reduction of the value that determines convergence.
 */
static double
rhea_newton_status_get_reduction_prev (rhea_newton_status_t *status)
{
  double              reduction;

  switch (status->conv_criterion) {
  case RHEA_NEWTON_CONV_CRITERION_OBJECTIVE:
    reduction = status->obj_reduction_prev;
  case RHEA_NEWTON_CONV_CRITERION_GRAD_NORM:
    reduction = status->grad_norm_reduction_prev;
  default: /* unknown criterion */
    RHEA_ABORT_NOT_REACHED ();
    reduction = -1.0;
  }

  return reduction;
}

/**
 * Return current reduction of the value that determines convergence.
 */
static double
rhea_newton_status_get_reduction_curr (rhea_newton_status_t *status)
{
  double              reduction;

  switch (status->conv_criterion) {
  case RHEA_NEWTON_CONV_CRITERION_OBJECTIVE:
    reduction = status->obj_reduction_curr;
  case RHEA_NEWTON_CONV_CRITERION_GRAD_NORM:
    reduction = status->grad_norm_reduction_curr;
  default: /* unknown criterion */
    RHEA_ABORT_NOT_REACHED ();
    reduction = -1.0;
  }

  return reduction;
}

/**
 *
 */
static void
rhea_newton_status_print_curr (rhea_newton_status_t *status, const int iter,
                               const char *name)
{
  switch (status->conv_criterion) {
  case RHEA_NEWTON_CONV_CRITERION_OBJECTIVE:
    RHEA_GLOBAL_INFOF (
        "%s: Newton iter %i, objective functional %.6e, reduction %.3e\n",
        name, iter, status->obj_curr, status->obj_reduction);
    break;
  case RHEA_NEWTON_CONV_CRITERION_GRAD_NORM:
    RHEA_GLOBAL_INFOF (
        "%s: Newton iter %i, gradient norm %.6e, reduction %.3e\n",
        name, iter, status->grad_norm_curr, status->grad_norm_reduction);
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
                          rhea_newton_residual_t *residual,
                          rhea_newton_problem_t *nl_problem,
                          rhea_newton_options_t *opt)
{
  const char         *this_fn_name = "rhea_newton_compute_step";
  const int           iter = step->iter;
  const int           nonzero_initial_guess = 0;
  const int           lin_iter_max = opt->lin_iter_max;
  int                 lin_iter_count = -1;
  const double        lin_res_norm_rtol = step->lin_res_norm_rtol;
  double              lin_res_norm_reduction;
  double              lin_conv;

  RHEA_GLOBAL_INFOF ("Newton iter %i -- Into %s\n", iter, this_fn_name);

  /* check input */
  RHEA_ASSERT (0 <= step->iter);
  RHEA_ASSERT (0.0 < step->lin_res_norm_rtol && step->lin_res_norm_rtol < 1.0);

  /* run solver for the linearized system */ //TODO
  nl_problem->solve_hessian_sys (
      step->vec, residual->vec, lin_iter_max, lin_res_norm_rtol,
      nonzero_initial_guess, nl_problem->data, &lin_iter_count);
  RHEA_ASSERT (0 <= lin_iter_count);

  /* calculate the residual reduction of the linearized solve */
  if (nl_problem->apply_hessian != NULL) { /* if function provided */
    ymir_vec_t         *lin_residual_vec;
    double              lin_res_norm_init;
    double              lin_res_norm_curr;

    /* compute the inital l^2-norm of the residual of the linearized system */ //TODO
    lin_res_norm_init = ymir_vec_norm (residual->vec);

    /* compute the new l^2-norm of the residual of the linearized system */ //TODO
    lin_residual_vec = ymir_vec_template (residual->vec);
    nl_problem->apply_hessian (lin_residual_vec, step->vec, nl_problem->data);
    ymir_vec_add (-1.0, residual->vec, lin_residual_vec);
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

//TODO
static int
rhea_newton_search_step_length_check_descend (rhea_newton_step_t *step,
                                              rhea_newton_residual_t *residual,
                                              rhea_newton_options_t *opt,
                                              const int print_condition)
{
  const char         *this_fn_name = "rhea_newton_search_step_length";
  const int           iter = step->iter;
  const double        lin_reduction = step->lin_res_norm_reduction;
  const double        relax = opt->step_descend_condition_relaxation;
  const double        reduction = residual->norm_reduction_curr;
  double              descend_reduction;

  /* check input */
  RHEA_ASSERT (0.0 < step->lin_res_norm_reduction);
  RHEA_ASSERT (step->lin_res_norm_reduction <= 1.0);
  RHEA_ASSERT (0.0 < residual->norm_reduction_curr);

  descend_reduction = 1.0 - relax * (1.0 - lin_reduction);

  if (print_condition) {
    RHEA_GLOBAL_INFOF (
        "Newton iter %i -- %s: Descend condition: residual norm <= %.3e\n",
        iter, this_fn_name, descend_reduction);
  }

  return (reduction <= descend_reduction);
}

/**
 * Performs line search for the step length (updates the Jacobian).
 */
static void
rhea_newton_search_step_length (ymir_vec_t *solution,
                                rhea_newton_step_t *step,
                                rhea_newton_residual_t *residual,
                                rhea_newton_problem_t *nl_problem,
                                rhea_newton_options_t *opt)
{
  const char         *this_fn_name = "rhea_newton_search_step_length";
  const int           iter = step->iter;
  ymir_vec_t         *solution_prev = ymir_vec_template (solution);
  const int           search_iter_max = opt->step_search_iter_max;
  int                 search_id;
  int                 search_success = 0;
  const double        step_length_min = opt->step_length_min;
  const double        step_reduction = opt->step_reduction;
  double              step_length = opt->step_length_max;
  double              res_norm, *res_norm_comp;

  RHEA_GLOBAL_INFOF ("Newton iter %i -- Into %s\n", iter, this_fn_name);

  /* check input */
  RHEA_ASSERT (0 <= step->iter);

  /* create work variables */
  if (0 < residual->norm_n_components) {
    res_norm_comp = RHEA_ALLOC (double, residual->norm_n_components);
  }
  else {
    res_norm_comp = NULL;
  }

  /* copy previous solution and residual norms */
  ymir_vec_copy (solution, solution_prev);
  rhea_newton_residual_norm_copy_curr_to_prev (residual);

  /* search for step length */
  for (search_id = 1; search_id <= search_iter_max; search_id++) {
    /* move solution into step direction */
    ymir_vec_add (step_length, step->vec, solution);

    /* update the nonlinear operator at new solution */
    nl_problem->update_operator (solution, nl_problem->data);

    /* compute new Newton right-hand side */ //TODO
    nl_problem->compute_neg_gradient (residual->vec, solution, nl_problem->data);
    res_norm = nl_problem->compute_gradient_norm (residual->vec,
                                                  nl_problem->data,
                                                  res_norm_comp);
    rhea_newton_residual_norm_set_curr (residual, res_norm, res_norm_comp);

    /* check descend condition */
    search_success = rhea_newton_search_step_length_check_descend (
        step, residual, opt, (0 == search_id) /* print only once */);

    RHEA_GLOBAL_INFOF ("Newton iter %i -- %s: check # %i, step length %.3e, "
                       "residual norm %.3e, previous %.3e\n",
                       iter, this_fn_name, search_id, step_length,
                       residual->norm_curr, residual->norm_prev);

    /* continue or terminate loop */
    if (search_success) {
      /* stop search successfully if descend condition is satisfied */
      break;
    }
    else {
      /* reduce step length */
      step_length *= step_reduction;

      /* stop search if step length is too low */
      if (step_length < step_length_min) {
        step_length = 0.0;
        break;
      }
    }
  }

  /* post-processing depending on whether a step length was found */
  if (search_success) {
    RHEA_GLOBAL_INFOF (
        "Newton iter %i -- %s: Line search successful, step length %.3e\n",
        iter, this_fn_name, step_length);

    /* update the linearized operator */
    nl_problem->update_hessian (solution, nl_problem->data);
  }
  else {
    /* check and print why search did not succeed */
    if (0.0 < step_length) {
      step_length = -1.0;
      RHEA_GLOBAL_INFOF (
          "Newton iter %i -- %s: Line search failed, "
          "max number of step reductions reached\n", iter, this_fn_name);
    }
    else {
      step_length = 0.0;
      RHEA_GLOBAL_INFOF (
          "Newton iter %i -- %s: Line search failed, "
          "min step length reached\n", iter, this_fn_name);
    }

    /* reverse updates of nonlinear operator */
    nl_problem->update_operator (solution_prev, nl_problem->data);
  }

  /* destroy */
  ymir_vec_destroy (solution_prev);
  if (0 < residual->norm_n_components) {
    RHEA_FREE (res_norm_comp);
  }

  RHEA_GLOBAL_INFOF ("Newton iter %i -- Done %s\n", iter, this_fn_name);
}

void
rhea_newton_solve (ymir_vec_t *solution,
                   rhea_newton_problem_t *nl_problem,
                   rhea_newton_options_t *opt)
{
  const char         *this_fn_name = "rhea_newton_solve";
  const int           print_summary = opt->print_summary;
  const int           iter_start = opt->iter_start;
  const int           iter_max = opt->iter_max;
  int                 iter;
  const double        rtol = opt->rtol;
  rhea_newton_step_t  step;
  rhea_newton_status_t  status;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /*
   * Initialize
   */
  {
    const int           n_components =
                          nl_problem->gradient_norm_multi_components;
    double              obj_val, grad_norm, *grad_norm_comp;

    /* set initial guess */
    ymir_vec_set_zero (solution);

    /* create step */
    rhea_newton_step_init (&step, nl_problem->neg_gradient_vec,
                           nl_problem->step_vec);

    /* compute initial right-hand side */
    nl_problem->compute_neg_gradient (step.rhs_vec, solution, nl_problem->data);

    /* create status */
    if (0 < n_components) {
      grad_norm_comp = RHEA_ALLOC (double, n_components);
    }
    else {
      grad_norm_comp = NULL;
    }
    rhea_newton_status_init (&status, n_components);

    /* set initial status */
    obj_val = nl_problem->evaluate_objective (solution, nl_problem->data);
    grad_norm = nl_problem->compute_gradient_norm (
        step.rhs_vec, nl_problem->data, grad_norm_comp);
    rhea_newton_status_set_curr (&status, obj_val, grad_norm, grad_norm_comp);
    rhea_newton_status_copy_curr_to_init (&status);

    /* destroy */
    if (0 < n_components) {
      RHEA_FREE (grad_norm_comp);
    }
  }

  /*
   * Iterations Loop
   */

  for (iter = iter_start; iter < iter_max; iter++) { /* BEGIN: Newton iter */

    /*
     * Pre-Step Output
     */
    {
      rhea_newton_status_print_curr (&status, iter, this_fn_name);

      //TODO
    }

    /*
     * Newton Step
     */
    {
      /* check for convergence */
      if (rhea_newton_status_get_reduction (&status) < rtol) {
        RHEA_GLOBAL_PRODUCTIONF (
            "%s: Nonlinear problem converged to rtol (%.3e)\n",
            this_fn_name, rtol);
        break;
      }

      /* calculate the accuracy for computing the step, i.e., solving the
       * linearized system */
      rhea_newton_set_accuracy (
          /* out: */ &step,
          /* in:  */ &status, opt);

      /* solve the linearized system for an inexact Newton step */
//TODO
//    rhea_newton_compute_step (
//        /* out: */ &step,
//        /* in:  */ &residual, nl_problem, opt);

      /* perform line search for the step length (updates the solution and the
       * nonlinear operator) */
//TODO
//    rhea_newton_search_step_length (
//        /* out: */ solution, &step, &residual,
//        /* in:  */ nl_problem, opt);
    }

    /*
     * Post-Step Output
     */
    {
      RHEA_GLOBAL_INFOF (
          "%s: Newton iter %i -- Step length %g\n",
          this_fn_name, iter, step.length);

      RHEA_GLOBAL_INFOF (
          "%s: Newton iter %i -- Linear solver: num iter %i, "
          "residual reduction %.3e, prescribed rtol %.3e, "
          "convergence %.3e (= rtol^(1/#iter))\n",
          this_fn_name, iter, step.lin_iter_count,
          step.lin_res_norm_reduction, step.lin_res_norm_rtol,
          step.lin_convergence);

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
     * Setup Next Step
     */
    {
      //TODO
    }

  } /* END: Newton iter */

  /* check if max #iterations was reached without convergence */
  if (iter_max == iter && rtol <= rhea_newton_status_get_reduction (&status)) {
    RHEA_GLOBAL_PRODUCTIONF (
        "%s: Maximum number of nonlinear iterations reached (%i)\n",
        this_fn_name, iter_max);
  }

  /*
   * Finalize
   */
  {
    /* destroy */
    rhea_newton_status_clear (&status);

    //TODO
  }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}
