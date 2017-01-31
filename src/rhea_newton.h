/*
 */

#ifndef RHEA_NEWTON_H
#define RHEA_NEWTON_H

#include <ymir_vec_ops.h>

/**
 * Callback function for Newton's method.
 * Initializes user data of the nonlinear problem to be able to run the Newton
 * solver.
 *
 * \param [in] solution Solution vector
 * \param [in] data     User data
 */
typedef void      (*rhea_newton_data_init_fn_t) (ymir_vec_t *solution,
                                                 void *data);

/**
 * Callback function for Newton's method.
 * Clears user data of the nonlinear problem after the Newton solve.
 *
 * \param [in] data     User data
 */
typedef void      (*rhea_newton_data_clear_fn_t) (void *data);

/**
 * Callback function for Newton's method.
 * Evaluates the objective functional at the current solution vector.
 *
 * \return              Value of objective functional
 * \param [in] solution Current solution vector
 * \param [in] data     User data
 */
typedef double    (*rhea_newton_evaluate_objective_fn_t) (
                                            ymir_vec_t *solution, void *data);

/**
 * Callback function for Newton's method.
 * Computes the negative gradient at the current solution vector.
 *
 * \param [out] neg_gradient  Negative gradient vector
 * \param [in] solution       Current solution vector
 * \param [in] data           User data
 */
typedef void      (*rhea_newton_compute_negative_gradient_fn_t) (
                                            ymir_vec_t *neg_gradient,
                                            ymir_vec_t *solution, void *data);

/**
 * Callback function for Newton's method.
 * Computes norm of the (negative) gradient vector.
 *
 * \return                  Norm of the gradient vector
 * \param [in] neg_gradient Negative gradient vector
 * \param [in] data         User data
 * \param [out] norm_comp   Components of the gradient norm
 */
typedef double    (*rhea_newton_compute_norm_of_gradient_fn_t) (
                                            ymir_vec_t *neg_gradient,
                                            void *data, double *norm_comp);

/**
 * Callback function for Newton's method.
 * Applies the Hessian operator to a vector.
 *
 * \param [out] out Result of Hessian application
 * \param [in] in   Input vector
 * \param [in] data User data
 */
typedef void      (*rhea_newton_apply_hessian_fn_t) (
                                            ymir_vec_t *out, ymir_vec_t *in,
                                            void *data);

/**
 * Callback function for Newton's method.
 * Inverts the Hessian operator approximatively, i.e., up to a given tolerance.
 *
 * \return                            Stopping reason
 * \param [out] step                  Result of Hessian application
 * \param [in] neg_gradient           Negative gradient vector (serves as RHS)
 * \param [in] lin_iter_max           Max #iterations for iterative solver
 * \param [in] lin_res_norm_rtol      Relative tolerance for iterative solver
 * \param [in] nonzero_initial_guess  Flag if initial guess is nonzero
 * \param [in] data                   User data
 * \param [out] lin_iter_count        Number of iterations taken
 */
typedef int       (*rhea_newton_solve_hessian_system_fn_t) (
                                            ymir_vec_t *step,
                                            ymir_vec_t *neg_gradient,
                                            const int lin_iter_max,
                                            const double lin_res_norm_rtol,
                                            const int nonzero_initial_guess,
                                            void *data,
                                            int *lin_iter_count);

/**
 * Callback function for Newton's method.
 * Updates the current solution of the (nonlinear) operator.
 *
 * \param [in] solution Current solution vector
 * \param [in] data     User data
 */
typedef void      (*rhea_newton_update_operator_fn_t) (
                                            ymir_vec_t *solution, void *data);

/**
 * Callback function for Newton's method.
 * Updates the current solution of the Hessian operator.
 *
 * \param [in] solution Current solution vector
 * \param [in] data     User data
 */
typedef void      (*rhea_newton_update_hessian_fn_t) (
                                            ymir_vec_t *solution, void *data);

/* enumerator for convergence critera */
typedef enum
{
  RHEA_NEWTON_CONV_CRITERION_NONE = -1,
  RHEA_NEWTON_CONV_CRITERION_OBJECTIVE,
  RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM,
  RHEA_NEWTON_CONV_CRITERION_RESIDUAL_NORM
}
rhea_newton_conv_criterion_t;

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

/* Nonlinear problem (opaque) */
typedef struct rhea_newton_problem rhea_newton_problem_t;

/**
 * Creates a new nonlinear problem.
 */
rhea_newton_problem_t *rhea_newton_problem_new (
              ymir_vec_t *neg_gradient_vec,
              ymir_vec_t *step_vec,
              rhea_newton_compute_negative_gradient_fn_t compute_neg_gradient,
              rhea_newton_solve_hessian_system_fn_t solve_hessian_sys);

/**
 * Destroys a nonlinear problem.
 */
void                rhea_newton_problem_destroy (
              rhea_newton_problem_t *nl_problem);

/**
 * Sets callback function for initializing the Newton solver.
 */
void                rhea_newton_problem_set_data_fn (
              void *data,
              rhea_newton_data_init_fn_t data_init,
              rhea_newton_data_clear_fn_t data_clear,
              rhea_newton_problem_t *nl_problem);

/**
 * Sets callback functions related to the convergence criterion.
 */
void                rhea_newton_problem_set_conv_criterion_fn (
              rhea_newton_conv_criterion_t conv_criterion,
              rhea_newton_evaluate_objective_fn_t evaluate_objective,
              rhea_newton_compute_norm_of_gradient_fn_t compute_gradient_norm,
              const int grad_norm_multi_components,
              rhea_newton_problem_t *nl_problem);

/**
 * Sets callback function for applying the Hessian.
 */
void                rhea_newton_problem_set_apply_hessian_fn (
              rhea_newton_apply_hessian_fn_t apply_hessian,
              rhea_newton_problem_t *nl_problem);

/**
 * Sets callback function for updating the nonlinear and Hessian operators.
 */
void                rhea_newton_problem_set_update_fn (
              rhea_newton_update_operator_fn_t update_operator,
              rhea_newton_update_hessian_fn_t update_hessian,
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
