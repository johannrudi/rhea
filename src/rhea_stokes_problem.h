/*
 */

#ifndef RHEA_STOKES_PROBLEM_H
#define RHEA_STOKES_PROBLEM_H

#include <rhea_domain.h>
#include <rhea_temperature.h>
#include <rhea_viscosity.h>
#include <ymir_pressure_elem.h>

/* Stokes problem (opaque) */
typedef struct rhea_stokes_problem rhea_stokes_problem_t;

/**
 * Creates a new Stokes problem.
 */
rhea_stokes_problem_t *rhea_stokes_problem_new (
                                    ymir_vec_t *temperature,
                                    ymir_vec_t *weakzone,
                                    ymir_mesh_t *ymir_mesh,
                                    ymir_pressure_elem_t *press_elem,
                                    rhea_domain_options_t *domain_options,
                                    rhea_temperature_options_t *temp_options,
                                    rhea_viscosity_options_t *visc_options);

/**
 * Destroys a Stokes problem.
 */
void                rhea_stokes_problem_destroy (
                                    rhea_stokes_problem_t *stokes_problem);

/**
 * Sets up the solver for a Stokes problem.
 */
void                rhea_stokes_problem_setup_solver (
                                    rhea_stokes_problem_t *stokes_problem);

/**
 * Solves a Stokes problem.
 */
void                rhea_stokes_problem_solve (
                                    ymir_vec_t *sol_vel_press,
                                    const double rel_tol,
                                    const int maxiter,
                                    rhea_stokes_problem_t *stokes_problem);

/**
 * Gets viscosity vector from Stokes problem.
 */
void                rhea_stokes_problem_get_viscosity (
                                    ymir_vec_t *viscosity,
                                    rhea_stokes_problem_t *stokes_problem);

/**
 * Sets velocity components on the boundary, which are constrained by Dirichlet
 * boundary conditions, to zero.
 */
void                rhea_stokes_problem_velocity_boundary_set_zero (
                                    ymir_vec_t *velocity,
                                    rhea_stokes_problem_t *stokes_problem);

/* Analogous function declarations for a linear Stokes problem */

rhea_stokes_problem_t *rhea_stokes_problem_linear_new (
                                    ymir_vec_t *temperature,
                                    ymir_vec_t *weakzone,
                                    ymir_mesh_t *ymir_mesh,
                                    ymir_pressure_elem_t *press_elem,
                                    rhea_domain_options_t *domain_options,
                                    rhea_temperature_options_t *temp_options,
                                    rhea_viscosity_options_t *visc_options);

void                rhea_stokes_problem_linear_destroy (
                                    rhea_stokes_problem_t *stokes_problem_lin);

void                rhea_stokes_problem_linear_setup_solver (
                                    rhea_stokes_problem_t *stokes_problem_lin);

void                rhea_stokes_problem_linear_solve (
                                    ymir_vec_t *sol_vel_press,
                                    const double rel_tol,
                                    const int maxiter,
                                    rhea_stokes_problem_t *stokes_problem_lin);

/* Analogous function declarations for a nonlinear Stokes problem */

rhea_stokes_problem_t *rhea_stokes_problem_nonlinear_new (
                                    ymir_vec_t *temperature,
                                    ymir_vec_t *weakzone,
                                    ymir_mesh_t *ymir_mesh,
                                    ymir_pressure_elem_t *press_elem,
                                    rhea_domain_options_t *domain_options,
                                    rhea_temperature_options_t *temp_options,
                                    rhea_viscosity_options_t *visc_options);

void                rhea_stokes_problem_nonlinear_destroy (
                                    rhea_stokes_problem_t *stokes_problem_nl);

void                rhea_stokes_problem_nonlinear_setup_solver (
                                    rhea_stokes_problem_t *stokes_problem_nl);

void                rhea_stokes_problem_nonlinear_solve (
                                    ymir_vec_t *sol_vel_press,
                                    const double rel_tol,
                                    const int maxiter,
                                    rhea_stokes_problem_t *stokes_problem_nl);

void                rhea_stokes_problem_nonlinear_setup_solver (
                                    rhea_stokes_problem_t *stokes_problem_nl);

void                rhea_stokes_problem_nonlinear_solve (
                                    ymir_vec_t *sol_vel_press,
                                    const double rel_tol,
                                    const int maxiter,
                                    rhea_stokes_problem_t *stokes_problem_nl);

#endif /* RHEA_STOKES_PROBLEM_H */
