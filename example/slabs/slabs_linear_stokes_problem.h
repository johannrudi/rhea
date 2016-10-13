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

#ifndef SLABS_LINEAR_STOKES_PROBLEM_H
#define SLABS_LINEAR_STOKES_PROBLEM_H

#include <ymir_stokes_op.h>
#include <ymir_stokes_pc.h>
#include <slabs_stokes_state.h>
#include <slabs_physics.h>

#define SL_DIRSCAL_MIN 1.0e-4

/* structure for linear Stokes problem */
typedef struct slabs_lin_stokes_problem
{
  ymir_mesh_t        *mesh;
  ymir_pressure_elem_t  *press_elem;
  ymir_vel_dir_t     *vel_dir;

  ymir_vec_t         *rhs_u_point;
  ymir_vec_t         *rhs;
  ymir_vec_t         *viscosity;

  ymir_stokes_op_t   *stokes_op;
  ymir_stokes_pc_t   *stokes_pc;
}
slabs_lin_stokes_problem_t;

/**
 * Sets Dirichlet boundary conditions.
 */
ymir_vel_dir_t *
slabs_set_dirichlet_bc (ymir_mesh_t *mesh, ymir_vec_t *dirscal, void *data);

/**
 * Creates a new Stokes operator.
 */
ymir_stokes_op_t *
slabs_stokes_op_new (ymir_vec_t *viscosity,
                     ymir_vel_dir_t *vel_dir,
                     ymir_vel_rob_t *vel_rob,
                     ymir_vec_t *usol,
                     ymir_vec_t *dvdIIe,
                     ymir_pressure_elem_t *press_elem,
                     slabs_physics_options_t *physics_options);

/**
 * Destroys a Stokes operator.
 */
void
slabs_stokes_op_destroy (ymir_stokes_op_t *stokes_op);

/**
 * Creates a new linear Stokes problem.
 */
slabs_lin_stokes_problem_t *
slabs_linear_stokes_problem_new (slabs_stokes_state_t *state,
                                 ymir_mesh_t *mesh,
                                 ymir_pressure_elem_t *press_elem,
                                 slabs_physics_options_t *physics_options);

/**
 * Destroys a linear Stokes problem.
 */
void
slabs_linear_stokes_problem_destroy (slabs_lin_stokes_problem_t *lin_stokes);

/**
 * Projects out null spaces of velocity and pressure.
 */
void
slabs_linear_stokes_problem_project_out_nullspace (ymir_vec_t *up,
                                                   ymir_stokes_op_t *stokes_op,
                                                   slabs_physics_options_t
                                                     *physics_options,
                                                   const int residual_space);

#endif /* SLABS_LINEAR_STOKES_PROBLEM_H */
