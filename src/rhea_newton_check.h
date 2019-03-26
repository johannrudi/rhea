/*
 */

#ifndef RHEA_NEWTON_CHECK_H
#define RHEA_NEWTON_CHECK_H

#include <rhea_newton.h>

/**
 * Checks gradient of an objective functional by comparing to the finite
 * difference gradient of the functional.
 */
void                rhea_newton_check_gradient (
                                            ymir_vec_t *solution,
                                            ymir_vec_t *neg_gradient,
                                            const int iter,
                                            rhea_newton_problem_t *nl_problem);

/**
 * Checks Hessian of an objective functional by comparing to the finite
 * difference derivative of the functianal's gradient.
 */
void                rhea_newton_check_hessian (
                                            ymir_vec_t *solution,
                                            ymir_vec_t *neg_gradient,
                                            const int iter,
                                            rhea_newton_problem_t *nl_problem);

#endif /* RHEA_NEWTON_CHECK_H */
