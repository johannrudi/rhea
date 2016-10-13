/*
  This file is part of the ymir Library.
  ymir is a C library for modeling ice sheets

  Copyright (C) 2015 Carsten Burstedde, Toby Isaac, Johann Rudi, Georg Stadler,
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

  This example runs perfomace tests.

*/

#ifndef SLABS_SETUP_H
#define SLABS_SETUP_H

#include <slabs_physics.h>
#include <slabs_discretization.h>
#include <slabs_linear_stokes_problem.h>
#include <slabs_nonlinear_stokes_problem.h>
#include <slabs_norm.h>
#include <slabs_io.h>

#include <ymir_base.h>
#include <ymir_mesh.h>
#include <ymir_pressure_elem.h>
#include <ymir_perf_counter.h>

//TODO this has to go into a nonlinear solver header (that doesn't exist yet)
/* enumerator for initial guess type */
typedef enum
{
  SL_NL_SOLVER_INITIAL_GUESS_ZERO,
  SL_NL_SOLVER_INITIAL_GUESS_NONZERO

//SL_NL_SOLVER_INITIAL_GUESS_URHS_PZERO,
//SL_NL_SOLVER_INITIAL_GUESS_URHS_PNONZERO
}
slabs_nl_solver_initial_guess_t;

/* enumerator for step length reduction type */
typedef enum
{
  SL_NL_SOLVER_STEP_REDUCTION_CONST,
  SL_NL_SOLVER_STEP_REDUCTION_ADAPTIVE,
  SL_NL_SOLVER_STEP_REDUCTION_QUAD_INTERP
}
slabs_nl_solver_step_reduction_type_t;

/* enumerator for nullspace projection type */
typedef enum
{
  SL_NL_SOLVER_PROJECT_NULLSPACE_NONE,
  SL_NL_SOLVER_PROJECT_NULLSPACE,
  SL_NL_SOLVER_PROJECT_NULLSPACE_RES,
  SL_NL_SOLVER_PROJECT_NULLSPACE_SYMM
}
slabs_nl_solver_project_nullspace_t;

/* options for nonlinear solver */
typedef struct slabs_nl_solver_options
{
  slabs_nl_solver_type_t  nl_solver_type;
  slabs_nl_solver_primaldual_type_t  nl_solver_primaldual_type;
  slabs_nl_solver_primaldual_scal_type_t  nl_solver_primaldual_scal_type;
  int                 nl_maxiter;
  double              nl_rtol;
  slabs_nl_solver_initial_guess_t  initial_guess_type;
  int                 nl_resume_at_iter;
  double              nl_resume_at_time;
  double              nl_resume_prev_res;
  double              nl_resume_init_res;

  double              nl_forcing_exponent;
  double              nl_forcing_max;
  int                 nl_forcing_max_progressive_iter;
  double              nl_forcing_total_min;
  int                 nl_forcing_saveguard;
  double              nl_forcing_saveguard_threshold;

  double              nl_step_length_descend_cond_relax;
  slabs_nl_solver_step_reduction_type_t  nl_step_length_reduction_type;
  double              nl_step_length_reduction_min;
  double              nl_step_length_reduction_max;
  double              nl_step_length_reduction_reg;
  double              nl_step_length_min;

  double              nl_switch_picard_step_length_min;
  int                 nl_switch_picard_after_amr;
  int                 nl_switch_picard_init;
  int                 nl_switch_picard_maxiter;
  double              nl_switch_picard_rtol;

  slabs_norm_type_t   norm_type;
  double              norm_Hminus1_mass_scaling;
  ymir_Hminus1_norm_op_t  *norm_op;

  slabs_nl_stokes_problem_schur_diag_t  schur_diag_type;
  slabs_nl_stokes_problem_scaling_t     scaling_type;
  slabs_nl_solver_project_nullspace_t   project_out_nullspace;
  int                 enforce_unscaled_reduction;

  double              grid_continuation_init_amr_threshold;
  double              grid_continuation_init_forcing;
  double              grid_continuation_init_forcing_exp;
  int                 grid_continuation_init_steps;
  int                 grid_continuation_skipsteps;
  int                 grid_continuation_maxsteps;
  double              viscosity_bounds_continuation_min;
  double              viscosity_bounds_continuation_max;
  int                 viscosity_bounds_continuation_steps;

  int                 nl_check_derivative;
  int                 log_physics_stats;

  slabs_krylov_type_t  krylov_type;
  int                 krylov_maxiter;
  double              krylov_atol;
  double              krylov_rtol;
  double              krylov_rtol_init;
  int                 krylov_gmres_num_vecs;
  double              krylov_residual_reduction; /* actual reduction */

  int                 lin_solve_stress_block_only;
  int                 lin_solve_press_bbt_only;
  int                 lin_solve_cnode_bbt_only;
}
slabs_nl_solver_options_t;

/**
 * Defines options and adds them as sub-options.
 */
void
slabs_setup_add_suboptions (ymir_options_t * opt_sup);

/**
 * Processes options.
 */
void
slabs_setup_process_options (slabs_physics_options_t *physics_options,
                             slabs_discr_options_t *discr_options,
                             slabs_nl_solver_options_t *solver_options);

/**
 * Sets up the mesh.
 */
void
slabs_setup_mesh (p8est_t **p8est,
                  ymir_mesh_t **mesh,
                  ymir_pressure_elem_t **press_elem,
                  slabs_stokes_state_t **state,
                  slabs_discr_enforce_refinement_data_t
                    **enforce_refinement_data,
                  slabs_physics_coarsen_stokes_coeff_data_t
                    **coarsen_coeff_data,
                  MPI_Comm mpicomm,
                  slabs_physics_options_t *physics_options,
                  slabs_discr_options_t *discr_options,
                  slabs_nl_solver_options_t *solver_options,
                  ymir_perf_counter_t * perf_counter,
                  const char *workload_filepath);

/**
 * Sets up a linear or nonlinear Stokes problem.
 */
void
slabs_setup_stokes (slabs_lin_stokes_problem_t **lin_stokes,
                    slabs_nl_stokes_problem_t **nl_stokes,
                    p8est_t *p8est,
                    ymir_mesh_t *mesh,
                    ymir_pressure_elem_t *press_elem,
                    slabs_stokes_state_t *state,
                    slabs_physics_options_t *physics_options,
                    slabs_nl_solver_options_t *solver_options,
                    ymir_perf_counter_t * perf_counter,
                    const char *workload_filepath);

/**
 * Cleans Up.
 */
void
slabs_clear (slabs_lin_stokes_problem_t *lin_stokes,
             slabs_nl_stokes_problem_t *nl_stokes,
             p8est_t *p8est,
             ymir_mesh_t *mesh,
             ymir_pressure_elem_t *press_elem,
             slabs_stokes_state_t *state,
             slabs_discr_enforce_refinement_data_t *enforce_refinement_data,
             slabs_physics_coarsen_stokes_coeff_data_t *coarsen_coeff_data,
             slabs_physics_options_t *physics_options,
             slabs_discr_options_t *discr_options);

#endif /* SLABS_SETUP_H */
