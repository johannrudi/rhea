/*
 */

#ifndef RHEA_NEWTON_H
#define RHEA_NEWTON_H

#include <ymir_vec_ops.h>

/* Nonlinear problem (opaque) */
typedef struct rhea_newton_problem rhea_newton_problem_t;

/* Newton options */
typedef struct rhea_newton_options rhea_newton_options_t;

/* callback functions for Newton's method */
typedef void      (*rhea_newton_problem_apply_operator_fn_t) (
                                            ymir_vec_t *out, ymir_vec_t *in,
                                            void *data);

typedef void      (*rhea_newton_problem_apply_jacobian_fn_t) (
                                            ymir_vec_t *out, ymir_vec_t *in,
                                            void *data);

typedef void      (*rhea_newton_problem_update_operator_fn_t) (
                                            ymir_vec_t *solution, void *data);

typedef void      (*rhea_newton_problem_update_jacobian_fn_t) (
                                            ymir_vec_t *solution, void *data);

typedef void      (*rhea_newton_problem_compute_residual_fn_t) (
                                            ymir_vec_t *residual,
                                            ymir_vec_t *solution, void *data);

typedef double    (*rhea_newton_problem_compute_residual_norm_fn_t) (
                                            ymir_vec_t *residual, void *data,
                                            double *res_norm_comp);

typedef int       (*rhea_newton_problem_solve_linearized_fn_t) (
                                            ymir_vec_t *step,
                                            ymir_vec_t *residual,
                                            const int lin_iter_max,
                                            const double lin_res_norm_rtol,
                                            const int nonzero_initial_guess,
                                            void *data,
                                            int *lin_iter_count);

/**
 * Creates a new nonlinear problem.
 */
rhea_newton_problem_t *rhea_newton_problem_new (
    ymir_vec_t *step_vec,
    ymir_vec_t *residual_vec,
    rhea_newton_problem_apply_operator_fn_t apply_operator,
    rhea_newton_problem_apply_jacobian_fn_t apply_jacobian,
    rhea_newton_problem_update_operator_fn_t update_operator,
    rhea_newton_problem_update_jacobian_fn_t update_jacobian,
    rhea_newton_problem_compute_residual_fn_t compute_residual,
    rhea_newton_problem_compute_residual_norm_fn_t compute_residual_norm,
    const int residual_norm_n_components,
    rhea_newton_problem_solve_linearized_fn_t solve_linearized,
    void *data);

/**
 * Destroys a nonlinear problem.
 */
void                rhea_newton_problem_destroy (
                                            rhea_newton_problem_t *nl_problem);

/**
 * Solves a nonlinear problem with inexact Newton--Krylov.
 */
void                rhea_newton_solve (ymir_vec_t *solution,
                                       rhea_newton_problem_t *nl_problem,
                                       rhea_newton_options_t *opt);

#endif /* RHEA_NEWTON_H */
