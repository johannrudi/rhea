/*
  This file is part of the ymir Library.
  ymir is a C library for modeling ice sheets

  Copyright (C) 2010, 2011 Carsten Burstedde, Toby Isaac, Georg Stadler,
                           Lucas Wilcox.

  The ymir Library is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The ymir Library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the ymir Library.  If not, see <http://www.gnu.org/licenses/>.

  ---

  This is the slabs example for global instantaneous mantle flow with plates.

*/

#ifndef SLABS_NONLINEAR_STOKES_PROBLEM_H
#define SLABS_NONLINEAR_STOKES_PROBLEM_H

#include <ymir_stokes_op.h>
#include <ymir_stokes_pc.h>
#include <slabs_stokes_state.h>
#include <slabs_physics.h>

/* structure for nonlinear Stokes problem */
typedef struct slabs_nl_stokes_problem
{
  ymir_mesh_t        *mesh;
  ymir_pressure_elem_t  *press_elem;
  ymir_vel_dir_t     *vel_dir;

  ymir_cvec_t        *rhs_u_point;
  ymir_dvec_t        *viscosity;
  ymir_dvec_t        *bounds_marker;
  ymir_dvec_t        *yielding_marker;

  ymir_stokes_op_t   *stokes_op;
  ymir_stokes_pc_t   *stokes_pc;
}
slabs_nl_stokes_problem_t;

/* enumerator for type of diagonal Schur complement approximation */
typedef enum
{
  SL_NL_STOKES_PROB_SCHUR_DIAG_INV_VISC_PMASS,
  SL_NL_STOKES_PROB_SCHUR_DIAG_INV_VISC_TENSOR_PMASS
}
slabs_nl_stokes_problem_schur_diag_t;

/* enumerator for type of diagonal scaling TODO deprecated */
typedef enum
{
  SL_NL_STOKES_PROB_SCALING_NONE,
  SL_NL_STOKES_PROB_SCALING_ADIAG_VISC_PMASS,
  SL_NL_STOKES_PROB_SCALING_ADIAG_VISC_TENSOR_PMASS,
  SL_NL_STOKES_PROB_SCALING_ADIAG_PMASS,
  SL_NL_STOKES_PROB_SCALING_ADIAG_BDIAG, /* for BFBT Schur PC only! */
  SL_NL_STOKES_PROB_SCALING_UMASS_PMASS,
  SL_NL_STOKES_PROB_SCALING_VISC_UMASS_PMASS,
  SL_NL_STOKES_PROB_SCALING_VISC_UMASS_VISC_PMASS,
  SL_NL_STOKES_PROB_SCALING_GRAD_VISC_PMASS,
  SL_NL_STOKES_PROB_SCALING_UUMASS_PPMASS,
  SL_NL_STOKES_PROB_SCALING_AU_BP
}
slabs_nl_stokes_problem_scaling_t;

/**
 * Creates a new nonlinear Stokes problem.
 */
slabs_nl_stokes_problem_t *
slabs_nonlinear_stokes_problem_new (slabs_stokes_state_t *state,
                                    ymir_mesh_t *mesh,
                                    ymir_pressure_elem_t *press_elem,
                                    slabs_physics_options_t *physics_options);

/**
 * Destroys nonlinear Stokes problem.
 */
void
slabs_nonlinear_stokes_problem_destroy (slabs_nl_stokes_problem_t *nl_stokes);

/**
 * Creates Stokes operator for nonlinear Stokes problem.
 */
void
slabs_nonlinear_stokes_op_init (
            slabs_nl_stokes_problem_t *nl_stokes,
            slabs_stokes_state_t *state,
            slabs_physics_options_t *physics_options,
            const int init,
            const slabs_nl_solver_type_t nl_solver_type,
            const slabs_nl_solver_primaldual_type_t primaldual_type,
            const slabs_nl_solver_primaldual_scal_type_t primaldual_scal_type,
            const int check_derivative);

/**
 * Updates Stokes operator for nonlinear Stokes problem.
 */
void
slabs_nonlinear_stokes_op_update (slabs_nl_stokes_problem_t *nl_stokes,
                                  slabs_stokes_state_t *state,
                                  slabs_physics_options_t *physics_options,
                                  const slabs_nl_solver_type_t nl_solver_type);

/**
 * Updates 4th order tensor of linearized Stokes coefficient.
 */
void
slabs_nonlinear_stokes_op_update_tensor (
            slabs_nl_stokes_problem_t *nl_stokes,
            slabs_stokes_state_t *state,
            slabs_physics_options_t *physics_options,
            const slabs_nl_solver_primaldual_type_t primaldual_type,
            const slabs_nl_solver_primaldual_scal_type_t primaldual_scal_type,
            const int check_derivative);

/**
 * Modifies nonlinear Stokes problem for Picard method.
 */
void
slabs_nonlinear_stokes_op_switch_picard (slabs_nl_stokes_problem_t *nl_stokes);

/**
 * Creates Stokes preconditioner for nonlinear Stokes problem.
 */
void
slabs_nonlinear_stokes_pc_init (
                    slabs_nl_stokes_problem_t *nl_stokes,
                    slabs_stokes_state_t *state,
                    const slabs_nl_stokes_problem_schur_diag_t schur_diag_type,
                    slabs_nl_stokes_problem_scaling_t scaling_type,
                    const int init);

/**
 * Updates Stokes preconditioner for nonlinear Stokes problem.
 */
ymir_stokes_pc_t *
slabs_nonlinear_stokes_pc_update (
                    slabs_nl_stokes_problem_t *nl_stokes,
                    slabs_stokes_state_t *state,
                    const slabs_nl_stokes_problem_schur_diag_t schur_diag_type,
                    const slabs_nl_stokes_problem_scaling_t scaling_type);

/*****************************************************************************
 * scaling stuff TODO deprecated
 ****************************************************************************/

/**
 * Computes diagonal pressure Schur complement approximation with inverse
 * viscosity tensor weighted lumped mass matrix.
 */
void
slabs_nl_stokes_problem_schur_diag_visc_tensor_mass (ymir_vec_t *schur_diag,
                                                     void *data);

/**
 * Computes velocity scaling with inverse lumped mass matrix.
 */
void
slabs_nl_stokes_problem_uscale_mass (ymir_vec_t *uscale, void *data);

/**
 * Computes velocity scaling with inverse viscosity weighted lumped mass matrix.
 */
void
slabs_nl_stokes_problem_uscale_visc_mass (ymir_vec_t *uscale, void *data);

/**
 * Computes velocity scaling with velocity weighted lumped mass matrix.
 */
void
slabs_nl_stokes_problem_uscale_vel_mass (ymir_vec_t *uscale, void *data);

/**
 * Computes velocity scaling by applying stress operator to velocity.
 */
void
slabs_nl_stokes_problem_uscale_Au (ymir_vec_t *uscale, void *data);

/**
 * Computes pressure scaling with inverse lumped mass matrix.
 */
void
slabs_nl_stokes_problem_pscale_mass (ymir_vec_t *pscale, void *data);

/**
 * Computes pressure scaling with inverse viscosity weighted lumped mass matrix.
 */
void
slabs_nl_stokes_problem_pscale_visc_mass (ymir_vec_t *pscale, void *data);

/**
 * Computes pressure scaling with inverse viscosity tensor weighted lumped mass
 * matrix.
 */
void
slabs_nl_stokes_problem_pscale_visc_tensor_mass (ymir_vec_t *pscale,
                                                 void *data);

/**
 * Computes pressure scaling with pressure weighted lumped mass matrix.
 */
void
slabs_nl_stokes_problem_pscale_press_mass (ymir_vec_t *pscale, void *data);

/**
 * Computes pressure scaling by applying gradient operator to pressure.
 */
void
slabs_nl_stokes_problem_pscale_Bp (ymir_vec_t *pscale, void *data);

/******************************************************************************
 * Tests
 *****************************************************************************/

/**
 * Tests Jacobian of Stokes operator by checking the equality:
 *
 *   [D F](x) * y = [d F](x;y)
 *
 * where
 *   Frechet derivative: [D F](x) = Jacobian of F(x)
 *   Gateaux derivative: [d F](x;y) = lim_{h -> 0} (F(x + h*y) - F(x)) / h
 *   `x` is current Stokes state
 *   `y` is random vector
 */
void
slabs_nonlinear_stokes_problem_test_current_deriv (
    ymir_stokes_op_t *stokes_op,
    slabs_stokes_state_t *state,
    ymir_mesh_t *mesh,
    ymir_pressure_elem_t *press_elem,
    ymir_vel_dir_t *vel_dir,
    ymir_dvec_t *viscosity,
    slabs_physics_options_t *physics_options);

#endif /* SLABS_NONLINEAR_STOKES_PROBLEM_H */
