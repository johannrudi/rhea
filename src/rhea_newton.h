/*
 */

#ifndef RHEA_NEWTON_H
#define RHEA_NEWTON_H

#include <ymir_vec_ops.h>

/* Nonlinear problem (opaque) */
typedef struct rhea_newton_problem rhea_newton_problem_t;

/* callback functions for Newton's method */
//TODO add comments explaining in/output
typedef double    (*rhea_newton_evaluate_objective_fn_t) (
                                            ymir_vec_t *solution, void *data);

typedef void      (*rhea_newton_compute_negative_gradient_fn_t) (
                                            ymir_vec_t *neg_gradient,
                                            ymir_vec_t *solution, void *data);

typedef double    (*rhea_newton_compute_norm_of_gradient_fn_t) (
                                            ymir_vec_t *neg_gradient,
                                            void *data, double *res_norm_comp);

typedef void      (*rhea_newton_apply_hessian_fn_t) (
                                            ymir_vec_t *out, ymir_vec_t *in,
                                            void *data);

typedef int       (*rhea_newton_solve_hessian_system_fn_t) (
                                            ymir_vec_t *step,
                                            ymir_vec_t *neg_gradient,
                                            const int lin_iter_max,
                                            const double lin_res_norm_rtol,
                                            const int nonzero_initial_guess,
                                            void *data,
                                            int *lin_iter_count);

typedef void      (*rhea_newton_update_operator_fn_t) (
                                            ymir_vec_t *solution, void *data);

typedef void      (*rhea_newton_update_hessian_fn_t) (
                                            ymir_vec_t *solution, void *data);

/* Newton options */
typedef struct rhea_newton_options
{
  int                 nonzero_initial_guess;
  int                 abort_failed_step_search;

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

  /* output options */
  int                 status_verbosity;
  int                 print_summary;
}
rhea_newton_options_t;

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_newton_add_options (ymir_options_t * opt_sup);

/**
 * Processes options and stores them.
 */
void                rhea_newton_process_options (rhea_newton_options_t *opt);

/**
 * Sets options to default values.
 */
void                rhea_newton_options_set_defaults (
                                                  rhea_newton_options_t *opt);
/**
 * Creates a new nonlinear problem.
 */
rhea_newton_problem_t *rhea_newton_problem_new (
              ymir_vec_t *neg_gradient_vec,
              ymir_vec_t *step_vec,
              rhea_newton_evaluate_objective_fn_t evaluate_objective,
              rhea_newton_compute_negative_gradient_fn_t compute_neg_gradient,
              rhea_newton_compute_norm_of_gradient_fn_t compute_gradient_norm,
              rhea_newton_apply_hessian_fn_t apply_hessian,
              rhea_newton_solve_hessian_system_fn_t solve_hessian_sys,
              rhea_newton_update_operator_fn_t update_operator,
              rhea_newton_update_hessian_fn_t update_hessian,
              const int grad_norm_multi_components,
              void *data);

/**
 * Destroys a nonlinear problem.
 */
void                rhea_newton_problem_destroy (
                                            rhea_newton_problem_t *nl_problem);

/**
 * Returns a pointer to the (negative) gradient vector.
 */
ymir_vec_t         *rhea_newton_problem_get_neg_gradient_vec (
                                            rhea_newton_problem_t *nl_problem);

/**
 * Returns a pointer to the step vector.
 */
ymir_vec_t         *rhea_newton_problem_get_step_vec (
                                            rhea_newton_problem_t *nl_problem);

/**
 * Returns a pointer to the user data.
 */
void               *rhea_newton_problem_get_data (
                                            rhea_newton_problem_t *nl_problem);

/**
 * Activates/deactivates gradient and Hessian checks.
 */
void                rhea_newton_problem_set_checks (
                                            const int check_gradient,
                                            const int check_hessian,
                                            rhea_newton_problem_t *nl_problem);

/**
 * Solves a nonlinear problem with inexact Newton--Krylov.
 */
void                rhea_newton_solve (ymir_vec_t *solution,
                                       rhea_newton_problem_t *nl_problem,
                                       rhea_newton_options_t *opt);

#endif /* RHEA_NEWTON_H */
