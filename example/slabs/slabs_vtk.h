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

#ifndef SLABS_VTK_H
#define SLABS_VTK_H

#include <slabs_physics.h>
#include <slabs_discretization.h>
#include <slabs_linear_stokes_problem.h>
#include <slabs_nonlinear_stokes_problem.h>

/**
 * Writes vtk output of the p8est mesh.
 */
void
slabs_vtk_write_p8est (const char *filepath, p8est_t *p4est,
                       const slabs_domain_shape_t domain_shape);

/**
 * Writes vtk output of mangll scalar fields.
 */
int
slabs_vtk_write_mangll (const char *filepath, mangll_t *ma, double *field,
                        char *name);

/**
 * Write VTK files for Stokes state fields and viscosity.
 */
void
slabs_vtk_write_state (const char *filepath,
                       slabs_stokes_state_t *state,
                       ymir_vec_t *viscosity,
                       const int compute_viscosity,
                       const int compute_stress_norm,
                       ymir_mesh_t *mesh,
                       ymir_pressure_elem_t *press_elem,
                       slabs_physics_options_t *physics_options);

/**
 * Writes vtk output of the input data for a linear Stokes problem.
 */
void
slabs_vtk_write_input_lin_stokes (const char *filepath,
                                  slabs_stokes_state_t *state,
                                  slabs_lin_stokes_problem_t *lin_stokes,
                                  slabs_physics_options_t *physics_options,
                                  slabs_discr_options_t *discr_options);

/**
 * Writes vtk output of the input data for a nonlinear Stokes problem.
 */
void
slabs_vtk_write_input_nl_stokes (const char *filepath,
                                 slabs_stokes_state_t *state,
                                 slabs_nl_stokes_problem_t *nl_stokes,
                                 slabs_physics_options_t *physics_options,
                                 slabs_discr_options_t *discr_options);

/**
 * Writes vtk output of the solution of a linear Stokes problem.
 */
void
slabs_vtk_write_solution_lin_stokes (const char *filepath,
                                     slabs_stokes_state_t *state,
                                     slabs_lin_stokes_problem_t *lin_stokes,
                                     const slabs_krylov_type_t krylov_type,
                                     slabs_physics_options_t *physics_options,
                                     slabs_discr_options_t *discr_options);

/**
 * Writes vtk output of the solution of a nonlinear Stokes problem.
 */
void
slabs_vtk_write_solution_nl_stokes (const char *filepath,
                                    slabs_stokes_state_t *state,
                                    slabs_nl_stokes_problem_t *nl_stokes,
                                    slabs_physics_options_t *physics_options,
                                    slabs_discr_options_t *discr_options);

/**
 * Writes vtk output of an iteration of the nonlinear solver.
 */
void
slabs_vtk_write_iteration_nl_stokes (const char *filepath,
                                     slabs_stokes_state_t *state,
                                     ymir_vec_t *step_up,
                                     ymir_vec_t *step_strain,
                                     ymir_vec_t *residual_lump_up,
                                     ymir_vec_t *residual_lump_strain,
                                     slabs_nl_stokes_problem_t *nl_stokes,
                                     slabs_physics_options_t *physics_options,
                                     slabs_discr_options_t *discr_options);

#endif /* SLABS_VTK_H */
