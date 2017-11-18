/*
 */

#ifndef RHEA_STOKES_PROBLEM_H
#define RHEA_STOKES_PROBLEM_H

#include <rhea_domain.h>
#include <rhea_viscosity.h>
#include <rhea_velocity.h>
#include <rhea_pressure.h>
#include <rhea_velocity_pressure.h>
#include <ymir_stokes_op.h>

/* types of linearization */
typedef enum
{
  RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_PICARD,
  RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_REGULAR,
  RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL,
  RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL_SYMM,
  RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_DEV1, //TODO change name
  RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_DEV2  //TODO change name
}
rhea_stokes_problem_nonlinear_linearization_t;

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_stokes_problem_add_options (ymir_options_t * opt_sup);

/* Stokes problem (opaque) */
typedef struct rhea_stokes_problem rhea_stokes_problem_t;

/**
 * Creates/destroys a Stokes problem.
 */
rhea_stokes_problem_t *rhea_stokes_problem_new (
                                    ymir_vec_t *temperature,
                                    ymir_vec_t *weakzone,
                                    ymir_vec_t *rhs_vel,
                                    ymir_vec_t *rhs_vel_nonzero_dirichlet,
                                    ymir_mesh_t *ymir_mesh,
                                    ymir_pressure_elem_t *press_elem,
                                    rhea_domain_options_t *domain_options,
                                    rhea_viscosity_options_t *visc_options,
                                    void *solver_options);

void                rhea_stokes_problem_destroy (
                                    rhea_stokes_problem_t *stokes_problem);

/**
 * Sets up/clears objects that have dependencies on the mesh.
 */
void                rhea_stokes_problem_create_mesh_dependencies (
                                    rhea_stokes_problem_t *stokes_problem,
                                    ymir_mesh_t *ymir_mesh,
                                    ymir_pressure_elem_t *press_elem);

void                rhea_stokes_problem_clear_mesh_dependencies (
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
                                    const int iter_max,
                                    const double rel_tol,
                                    rhea_stokes_problem_t *stokes_problem);

/**
 * Sets output path for vtk output the iterations of a nonlinear solve.
 */
void                rhea_stokes_problem_nonlinear_set_output (
                              char *vtk_write_newton_iteration_path,
                              rhea_stokes_problem_t *stokes_problem_nl);

/**
 * Accesses data of a Stokes problem.
 */
ymir_mesh_t        *rhea_stokes_problem_get_ymir_mesh (
                                    rhea_stokes_problem_t *stokes_problem);

ymir_pressure_elem_t *rhea_stokes_problem_get_press_elem (
                                    rhea_stokes_problem_t *stokes_problem);

ymir_vec_t         *rhea_stokes_problem_get_temperature (
                                    rhea_stokes_problem_t *stokes_problem);

ymir_vec_t         *rhea_stokes_problem_get_weakzone (
                                    rhea_stokes_problem_t *stokes_problem);

ymir_vec_t         *rhea_stokes_problem_get_rhs_vel (
                                    rhea_stokes_problem_t *stokes_problem);

ymir_vec_t         *rhea_stokes_problem_get_rhs_vel_nonzero_dirichlet (
                                    rhea_stokes_problem_t *stokes_problem);

rhea_domain_options_t    *rhea_stokes_problem_get_domain_options (
                                    rhea_stokes_problem_t *stokes_problem);

rhea_viscosity_options_t *rhea_stokes_problem_get_viscosity_options (
                                    rhea_stokes_problem_t *stokes_problem);

void                rhea_stokes_problem_copy_viscosity (
                                    ymir_vec_t *viscosity,
                                    rhea_stokes_problem_t *stokes_problem);

ymir_stokes_op_t   *rhea_stokes_problem_get_stokes_op (
                                    rhea_stokes_problem_t *stokes_problem);

/**
 * Sets velocity components on the boundary, which are constrained by Dirichlet
 * boundary conditions, to zero.
 */
void                rhea_stokes_problem_velocity_boundary_set_zero (
                                    ymir_vec_t *velocity,
                                    rhea_stokes_problem_t *stokes_problem);

#endif /* RHEA_STOKES_PROBLEM_H */
