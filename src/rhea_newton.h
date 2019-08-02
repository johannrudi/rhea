/** RHEA_NEWTON
 *
 * Nonlinear solver using Newton's method.
 */

#ifndef RHEA_NEWTON_H
#define RHEA_NEWTON_H

#include <ymir_vec_ops.h>

/******************************************************************************
 * Callback Function for Newton's Method
 *****************************************************************************/

/**
 * Initializes user data of the nonlinear problem to be able to run the Newton
 * solver.
 *
 * \param [in] solution Solution vector (may be NULL)
 * \param [in] data     User data
 */
typedef void      (*rhea_newton_data_init_fn_t) (ymir_vec_t *solution,
                                                 void *data);

/**
 * Clears user data of the nonlinear problem after the Newton solve.
 *
 * \param [in] data     User data
 */
typedef void      (*rhea_newton_data_clear_fn_t) (void *data);

/**
 * Evaluates the objective functional.
 *
 * \return              Value of objective functional
 * \param [in] solution Current solution vector (may be NULL)
 * \param [in] data     User data
 */
typedef double    (*rhea_newton_evaluate_objective_fn_t) (
                                            ymir_vec_t *solution, void *data,
                                            double *obj_comp);

/**
 * Computes the negative gradient of the objective functional.
 *
 * \param [out] neg_gradient  Negative gradient vector
 * \param [in] solution       Current solution vector (may be NULL)
 * \param [in] data           User data
 */
typedef void      (*rhea_newton_compute_negative_gradient_fn_t) (
                                            ymir_vec_t *neg_gradient,
                                            ymir_vec_t *solution, void *data);

/**
 * Computes the norm of the gradient.
 *
 * \return                  Norm of the gradient vector
 * \param [in] neg_gradient Negative gradient vector
 * \param [in] data         User data
 * \param [out] norm_comp   Components of the gradient norm
 */
typedef double    (*rhea_newton_compute_norm_of_gradient_fn_t) (
                                            ymir_vec_t *neg_gradient,
                                            void *data, double *grad_norm_comp);

/**
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
 * Applies inexact inverse of the Hessian operator.
 *
 * \return                            Stopping reason
 * \param [out] step                  Newton step
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
 * Modfies the Newton step after it has been.computed.
 *
 * \param [in/out] step Newton step (previously computed)
 * \param [in] solution Current solution vector
 * \param [in] data     User data
 */
typedef void      (*rhea_newton_modify_step_fn_t) (ymir_vec_t *step,
                                                   ymir_vec_t *solution,
                                                   void *data);

/**
 * Updates the nonlinear operator at the current solution vector.
 *
 * \param [in] solution Current solution vector (may be NULL)
 * \param [in] data     User data
 */
typedef void      (*rhea_newton_update_operator_fn_t) (
                                            ymir_vec_t *solution, void *data);

/**
 * Updates the Hessian operator at the current solution vector.
 *
 * \param [in] solution     Current solution vector (may be NULL)
 * \param [in] step_vec     Recent step vector, which (together with step
 *                          length) resulted in current solution (may be NULL)
 * \param [in] step_length  Step length corresponding to step vector
 * \param [in] data         User data
 */
typedef void      (*rhea_newton_update_hessian_fn_t) (
                                            ymir_vec_t *solution,
                                            ymir_vec_t *step_vec,
                                            const double step_length,
                                            void *data);

/**
 * Modifies the Hessian system before launching the solver for this system.
 *
 * \param [in/out] neg_gradient Negative gradient vector (serves as RHS)
 * \param [in] solution         Current solution vector (may be NULL)
 * \param [in] data             User data
 */
typedef void      (*rhea_newton_modify_hessian_system_fn_t) (
                                            ymir_vec_t *neg_gradient,
                                            ymir_vec_t *solution,
                                            void *data);

/**
 * Performs additional setup for next Newton step, e.g., AMR.
 *
 * \return              Flag whether to recompute Newton status
 * \param [in] solution Pointer to solution vector
 * \param [in] iter     Iteration number of Newton's method
 * \param [in] data     User data
 */
typedef int       (*rhea_newton_setup_poststep_fn_t) (ymir_vec_t **solution,
                                                      const int iter,
                                                      void *data);

/**
 * Writes or prints output at the beginning of a Newton step.
 *
 * \param [in] solution Solution vector
 * \param [in] iter     Iteration number of Newton's method
 * \param [in] data     User data
 */
typedef void      (*rhea_newton_output_prestep_fn_t) (ymir_vec_t *solution,
                                                      const int iter,
                                                      void *data);

/******************************************************************************
 * Options
 *****************************************************************************/

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

  int                 lin_monitor_reduction;

  /* options for line search */
  int                 step_search_iter_max;
  double              step_length_min;
  double              step_length_max;
  double              step_reduction;
  double              step_descend_condition_relaxation;

  /* options for resuming a nonlinear solve */
  int                 resume;
  double              resume_obj_init;
  double              resume_obj_prev;
  double              resume_obj_reduction_prev;
  double              resume_grad_norm_init;
  double              resume_grad_norm_prev;
  double              resume_grad_norm_reduction_prev;

  /* output options */
  int                 status_verbosity;
  int                 print_summary;
}
rhea_newton_options_t;

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_newton_add_options (
                                        rhea_newton_options_t *newton_options,
                                        ymir_options_t *opt_sup);

/**
 * Gets global options values.
 */
void                rhea_newton_get_options (rhea_newton_options_t *opt);

/**
 * Sets options to default values.
 */
void                rhea_newton_options_set_defaults (
                                                  rhea_newton_options_t *opt);

/******************************************************************************
 * Newton Problem
 *****************************************************************************/

/* Nonlinear problem (opaque) */
typedef struct rhea_newton_problem rhea_newton_problem_t;

/**
 * Creates a new nonlinear problem.
 */
rhea_newton_problem_t *rhea_newton_problem_new (
              rhea_newton_compute_negative_gradient_fn_t compute_neg_gradient,
              rhea_newton_compute_norm_of_gradient_fn_t compute_gradient_norm,
              const int grad_norm_multi_components,
              rhea_newton_solve_hessian_system_fn_t solve_hessian_sys);

/**
 * Destroys a nonlinear problem.
 */
void                rhea_newton_problem_destroy (
              rhea_newton_problem_t *nl_problem);

/**
 * Sets vectors required by the Newton solver.
 */
void                rhea_newton_problem_set_vectors (
              ymir_vec_t *neg_gradient_vec,
              ymir_vec_t *step_vec,
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
 * Sets callback function for evaluating the objective functional.
 */
void                rhea_newton_problem_set_evaluate_objective_fn (
              rhea_newton_evaluate_objective_fn_t evaluate_objective,
              const int obj_multi_components,
              rhea_newton_problem_t *nl_problem);

/**
 * Sets callback function for applying the Hessian.
 */
void                rhea_newton_problem_set_apply_hessian_fn (
              rhea_newton_apply_hessian_fn_t apply_hessian,
              rhea_newton_problem_t *nl_problem);

/**
 * Sets callback function for modfiying the Newton step after it has
 * been.computed.
 */
void                rhea_newton_problem_set_modify_step_fn (
              rhea_newton_modify_step_fn_t modify_step_vec,
              rhea_newton_problem_t *nl_problem);

/**
 * Sets callback functions for updating the nonlinear and Hessian operators.
 */
void                rhea_newton_problem_set_update_fn (
              rhea_newton_update_operator_fn_t update_operator,
              rhea_newton_update_hessian_fn_t update_hessian,
              rhea_newton_modify_hessian_system_fn_t modify_hessian_system,
              rhea_newton_problem_t *nl_problem);

/**
 * Sets callback functions for performing additional setup for next Newton step.
 */
void                rhea_newton_problem_set_setup_poststep_fn (
              rhea_newton_setup_poststep_fn_t setup_poststep,
              rhea_newton_problem_t *nl_problem);

/**
 * Sets callback function for output.
 */
void                rhea_newton_problem_set_output_fn (
              rhea_newton_output_prestep_fn_t output_prestep,
              rhea_newton_problem_t *nl_problem);

/**
 * Activates/deactivates gradient and Hessian checks.
 */
void                rhea_newton_problem_set_checks (
                                            const int check_gradient,
                                            const int check_hessian,
                                            rhea_newton_problem_t *nl_problem);

int                 rhea_newton_problem_get_check_gradient (
                                            rhea_newton_problem_t *nl_problem);

void                rhea_newton_problem_set_check_gradient (
                                            const int check_gradient,
                                            rhea_newton_problem_t *nl_problem);

int                 rhea_newton_problem_get_check_hessian (
                                            rhea_newton_problem_t *nl_problem);

void                rhea_newton_problem_set_check_hessian (
                                            const int check_hessian,
                                            rhea_newton_problem_t *nl_problem);

/**
 * Get/set MPI communicator.
 */
sc_MPI_Comm         rhea_newton_problem_get_mpicomm (
                                            rhea_newton_problem_t *nl_problem);

void                rhea_newton_problem_set_mpicomm (
                                            sc_MPI_Comm mpicomm,
                                            rhea_newton_problem_t *nl_problem);

/**
 * Solves a nonlinear problem with inexact Newton--Krylov.
 */
int                 rhea_newton_solve (ymir_vec_t **solution,
                                       rhea_newton_problem_t *nl_problem,
                                       rhea_newton_options_t *opt);

int                 rhea_newton_solve_has_converged (const int stop_reason);
int                 rhea_newton_solve_get_num_iterations (
                                            rhea_newton_problem_t *nl_problem);
double              rhea_newton_solve_get_residual_reduction (
                                            rhea_newton_problem_t *nl_problem);

/**
 * Accesses data and callback functions of a nonlinear problem.
 */
ymir_vec_t         *rhea_newton_problem_get_neg_gradient_vec (
                                            rhea_newton_problem_t *nl_problem);

ymir_vec_t         *rhea_newton_problem_get_step_vec (
                                            rhea_newton_problem_t *nl_problem);

void               *rhea_newton_problem_get_data (
                                            rhea_newton_problem_t *nl_problem);

int                 rhea_newton_problem_evaluate_objective_exists (
                                            rhea_newton_problem_t *nl_problem);
double              rhea_newton_problem_evaluate_objective (
                                            ymir_vec_t *solution,
                                            rhea_newton_problem_t *nl_problem,
                                            double *obj_comp);

int                 rhea_newton_problem_compute_gradient_norm_exists (
                                            rhea_newton_problem_t *nl_problem);
double              rhea_newton_problem_compute_gradient_norm (
                                            ymir_vec_t *neg_gradient,
                                            rhea_newton_problem_t *nl_problem,
                                            double *grad_norm_comp);

int                 rhea_newton_problem_compute_neg_gradient_exists (
                                            rhea_newton_problem_t *nl_problem);
void                rhea_newton_problem_compute_neg_gradient (
                                            ymir_vec_t *neg_gradient,
                                            ymir_vec_t *solution,
                                            rhea_newton_problem_t *nl_problem);

int                 rhea_newton_problem_apply_hessian_exists (
                                            rhea_newton_problem_t *nl_problem);
void                rhea_newton_problem_apply_hessian (
                                            ymir_vec_t *out, ymir_vec_t *in,
                                            rhea_newton_problem_t *nl_problem);

int                 rhea_newton_problem_update_operator_exists (
                                            rhea_newton_problem_t *nl_problem);
void                rhea_newton_problem_update_operator (
                                            ymir_vec_t *solution,
                                            rhea_newton_problem_t *nl_problem);

int                 rhea_newton_problem_update_hessian_exists (
                                            rhea_newton_problem_t *nl_problem);
void                rhea_newton_problem_update_hessian (
                                            ymir_vec_t *solution,
                                            ymir_vec_t *step_vec,
                                            const double step_length,
                                            rhea_newton_problem_t *nl_problem);

int                 rhea_newton_problem_modify_hessian_system_exists (
                                            rhea_newton_problem_t *nl_problem);
void                rhea_newton_problem_modify_hessian_system (
                                            ymir_vec_t *neg_gradient,
                                            ymir_vec_t *solution,
                                            rhea_newton_problem_t *nl_problem);

int                 rhea_newton_problem_setup_poststep_exists (
                                            rhea_newton_problem_t *nl_problem);
int                 rhea_newton_problem_setup_poststep (
                                            ymir_vec_t **solution,
                                            const int iter,
                                            rhea_newton_problem_t *nl_problem);

int                 rhea_newton_problem_output_prestep_exists (
                                            rhea_newton_problem_t *nl_problem);
void                rhea_newton_problem_output_prestep (
                                            ymir_vec_t *solution,
                                            const int iter,
                                            rhea_newton_problem_t *nl_problem);

/**
 * Accesses to status of a nonlinear problem.
 */
double              rhea_newton_problem_get_reduction_curr (
                                            rhea_newton_problem_t *nl_problem);

#endif /* RHEA_NEWTON_H */
