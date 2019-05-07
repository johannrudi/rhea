/*
 */

#ifndef NEWTON_POLYNOMIAL_BASE_H
#define NEWTON_POLYNOMIAL_BASE_H

#include <ymir_vec_ops.h>

/**
 * Computes interpolating quadratic polynomial via Hermite interpolation.
 */
double             *newton_polynomial_new_p2_coeff_hermite (
                                                      const double start_node,
                                                      const double start_val,
                                                      const double start_deriv,
                                                      const double end_node,
                                                      const double end_val);

/**
 * Evaluates quadratic polynomial.
 */
double              newton_polynomial_eval_p2 (const double x,
                                               const double coeff[3]);

/**
 * Evaluates 1st derivative of quadratic polynomial.
 */
double              newton_polynomial_eval_p2d (const double x,
                                                const double coeff[3]);

/**
 * Evaluates 2nd derivative of quadratic polynomial.
 */
double              newton_polynomial_eval_p2dd (const double x,
                                                 const double coeff[3]);

/**
 * Evaluates quadratic polynomial element-wise over a vector as input values.
 */
void                newton_polynomial_evaluate (ymir_vec_t *f_vec,
                                                ymir_vec_t *x_vec,
                                                const double coeff[3],
                                                const int derivative);

/* Object containing information about the problem */
typedef struct newton_polynomial_problem
{
  double              coeff[3];
  ymir_vec_t         *data_a;
  ymir_vec_t         *data_b;
  ymir_vec_t         *sol_x;
}
newton_polynomial_problem_t;

/**
 * Computes the misfit between the data and a trial solution vector.
 */
void                newton_polynomial_compute_misfit (
                                    ymir_vec_t *misfit_x,
                                    ymir_vec_t *misfit_y,
                                    ymir_vec_t *sol_x,
                                    newton_polynomial_problem_t *poly_problem);


/**
 * Computes the distance between the data and the polynomial surface.
 */
void                newton_polynomial_compute_distance (
                                    ymir_vec_t *dist, ymir_vec_t *x_vec,
                                    newton_polynomial_problem_t *poly_problem);

/******************************************************************************
 * Callback Functions for Newton Solver
 *****************************************************************************/

/**
 * Initializes user data of the nonlinear problem to be able to run the Newton
 * solver.
 * (Callback function for Newton's method.)
 */
void                newton_polynomial_data_init (ymir_vec_t *solution,
                                                 void *data);

/**
 * Evaluates the objective functional.
 * (Callback function for Newton's method.)
 */
double              newton_polynomial_evaluate_objective (ymir_vec_t *solution,
                                                          void *data,
                                                          double *obj_comp);

/**
 * Computes the negative gradient of the objective functional.
 * (Callback function for Newton's method.)
 */
void                newton_polynomial_compute_negative_gradient (
                                                      ymir_vec_t *neg_gradient,
                                                      ymir_vec_t *solution,
                                                      void *data);

/**
 * Computes the norm of the gradient.
 * (Callback function for Newton's method.)
 */
double              newton_polynomial_compute_gradient_norm (
                                                      ymir_vec_t *neg_gradient,
                                                      void *data,
                                                      double *norm_comp);

/**
 * Applies Hessian operator.
 * (Callback function for Newton's method.)
 */
void                newton_polynomial_apply_hessian (ymir_vec_t *out,
                                                     ymir_vec_t *in,
                                                     void *data);

/**
 * Applies inverse of the Hessian operator.
 * (Callback function for Newton's method.)
 */
int                 newton_polynomial_solve_hessian_system (
                                              ymir_vec_t *step,
                                              ymir_vec_t *neg_gradient,
                                              const int lin_iter_max,
                                              const double lin_res_norm_rtol,
                                              const int nonzero_initial_guess,
                                              void *data,
                                              int *lin_iter_count);

/**
 * Updates the Hessian operator, given a trial solution vector.
 * (Callback function for Newton's method.)
 */
void                newton_polynomial_update_hessian (ymir_vec_t *solution,
                                                      ymir_vec_t *step_vec,
                                                      const double step_length,
                                                      void *data);

#endif /* NEWTON_POLYNOMIAL_BASE_H */
