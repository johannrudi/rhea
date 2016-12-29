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

#include <slabs_nonlinear_stokes_problem.h>
#include <ymir_comm.h>
#include <ymir_mass_vec.h>
#include <ymir_velocity_vec.h>
#include <ymir_pressure_vec.h>
#include <ymir_interp_vec.h>
#include <slabs_linear_stokes_problem.h>
#include <slabs_norm.h>

/* declare local functions */
static void
slabs_nonlinear_stokes_problem_op_clear (slabs_nl_stokes_problem_t *nl_stokes);

static void
slabs_nonlinear_stokes_problem_pc_clear (slabs_nl_stokes_problem_t *nl_stokes);

static void
slabs_nl_stokes_problem_set_schur_diag (
                   slabs_nl_stokes_problem_t *nl_stokes,
                   slabs_stokes_state_t *state,
                   const slabs_nl_stokes_problem_schur_diag_t schur_diag_type);

static void
slabs_nl_stokes_problem_clear_schur_diag (slabs_nl_stokes_problem_t *nl_stokes);


/* scaling stuff TODO deprecated */
static void
slabs_nl_stokes_problem_set_scaling (slabs_nl_stokes_problem_t *nl_stokes,
                                     slabs_stokes_state_t *state,
                                     const slabs_nl_stokes_problem_scaling_t
                                       scaling_type);
static void
slabs_nl_stokes_problem_clear_scaling ();

/**
 * Creates a new nonlinear Stokes problem.
 */
slabs_nl_stokes_problem_t *
slabs_nonlinear_stokes_problem_new (slabs_stokes_state_t *state,
                                    ymir_mesh_t *mesh,
                                    ymir_pressure_elem_t *press_elem,
                                    slabs_physics_options_t *physics_options)
{
  const char         *this_fn_name = "slabs_nonlinear_stokes_problem_new";
  ymir_vec_t         *dirscal;
  ymir_vel_dir_t     *vel_dir;
  ymir_vec_t         *rhs_u_point;
  slabs_nl_stokes_problem_t  *nl_stokes;

  YMIR_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* initialize linear Stokes problem */
  nl_stokes = YMIR_ALLOC (slabs_nl_stokes_problem_t, 1);

  /* setup Dirichlet boundary conditions */
  if (physics_options->bc_default_dirichlet_scale) {
    dirscal = ymir_cvec_new (mesh, 1);
    ymir_vec_set_value (dirscal, 1.0);
  }
  else {
    dirscal = NULL;
  }
  vel_dir = slabs_set_dirichlet_bc (mesh, dirscal, physics_options);

  /* create right-hand side forcing */
  rhs_u_point = ymir_cvec_new (mesh, 3);
  slabs_physics_compute_rhs_u_point (rhs_u_point, state, physics_options);
  if (physics_options->rhs_random) { /* if set random rhs (for testing) */
#ifdef YMIR_PETSC
    ymir_petsc_vec_set_random (rhs_u_point, YMIR_MESH_PETSCLAYOUT_NONE);
#else
    YMIR_ABORT_NOT_REACHED ();
#endif
    if (physics_options->domain_shape == SL_DOMAIN_SHELL) { /* if shell */
      ymir_velocity_vec_project_out_mean_rotation (
          rhs_u_point, physics_options->domain_center,
          physics_options->domain_moment_of_inertia, 0);
    }
  }

  /* assign pointers to nonlinear Stokes problem structure */
  nl_stokes->mesh = mesh;
  nl_stokes->press_elem = press_elem;
  nl_stokes->vel_dir = vel_dir;
  nl_stokes->rhs_u_point = rhs_u_point;
  nl_stokes->viscosity = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
  nl_stokes->bounds_marker = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
  nl_stokes->yielding_marker = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);

  /* set NULL pointers */
  nl_stokes->stokes_op = NULL;
  nl_stokes->stokes_pc = NULL;

  YMIR_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

  return nl_stokes;
}

/**
 * Destroys nonlinear Stokes problem.
 */
void
slabs_nonlinear_stokes_problem_destroy (slabs_nl_stokes_problem_t *nl_stokes)
{
  const char         *this_fn_name = "slabs_nonlinear_stokes_problem_destroy";

  YMIR_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* destroy Stokes operator, preconditioner, and related variables */
  slabs_nonlinear_stokes_problem_pc_clear (nl_stokes);
  slabs_nonlinear_stokes_problem_op_clear (nl_stokes);

  /* destroy vectors */
  ymir_vec_destroy (nl_stokes->rhs_u_point);
  ymir_vec_destroy (nl_stokes->viscosity);
  ymir_vec_destroy (nl_stokes->bounds_marker);
  ymir_vec_destroy (nl_stokes->yielding_marker);

  /* destroy Dirichlet BC's */
  if (nl_stokes->vel_dir != NULL) {
    if (nl_stokes->vel_dir->scale != NULL) {
      ymir_vec_destroy (nl_stokes->vel_dir->scale);
    }
    ymir_vel_dir_destroy (nl_stokes->vel_dir);
  }

  YMIR_FREE (nl_stokes);

  YMIR_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 *
 */
static void
slabs_nonlinear_stokes_primaldual_nltensor (
                                  slabs_nl_stokes_problem_t *nl_stokes,
                                  slabs_stokes_state_t *state,
                                  const slabs_nl_solver_primaldual_type_t
                                    nl_solver_primaldual_type,
                                  const slabs_nl_solver_primaldual_scal_type_t
                                    nl_solver_primaldual_scal_type)
{
#ifdef YMIR_DEBUG
  const char         *this_fn_name =
                        "slabs_nonlinear_stokes_primaldual_nltensor";
  ymir_gloidx_t       n_scaled_nodes;
#endif
  ymir_stress_op_t   *stress_op = nl_stokes->stokes_op->stress_op;
  ymir_vec_t         *dual = state->dual_vec;
  ymir_vec_t         *nltensor;

  /* return if nothing to do */
  if (nl_solver_primaldual_type == SL_NL_SOLVER_PRIMALDUAL_NONE) {
    return;
  }

  /* initialize nonlinear tensor */
  nltensor = ymir_vec_template (dual);

  /* set nonlinear tensor for primal-dual method */
  switch (nl_solver_primaldual_type) {
  case SL_NL_SOLVER_PRIMALDUAL_NORMSTRAIN:
    {
      ymir_mesh_t        *mesh = nl_stokes->mesh;
      ymir_vec_t         *vel_vec;
      ymir_vec_t         *strain_frob = ymir_dvec_new (mesh, 1,
                                                       YMIR_GAUSS_NODE);
      ymir_velocity_elem_t  *vel_elem;

      /* set nonlinear tensor to scaled strain rate */
      ymir_dvec_copy (dual, nltensor);

      /* scale tensor */
      switch (nl_solver_primaldual_scal_type) {
      case SL_NL_SOLVER_PRIMALDUAL_SCAL_NONE:
#ifdef YMIR_DEBUG
        n_scaled_nodes = 0;
#endif
        break;

      case SL_NL_SOLVER_PRIMALDUAL_SCAL_NORMALIZE:
#ifdef YMIR_DEBUG
        n_scaled_nodes =
#endif
        slabs_norm_symtensor_normalize_frobenius (nltensor);
        break;

      case SL_NL_SOLVER_PRIMALDUAL_SCAL_THRESHOLD:
#ifdef YMIR_DEBUG
        n_scaled_nodes =
#endif
        slabs_norm_symtensor_normalize_threshold_frobenius (nltensor, 1.0);
        break;

      default: /* unknown primal-dual scaling type */
        YMIR_ABORT_NOT_REACHED ();
      }

      /* compute Frobenius-norm of the strain rate tensor */
      slabs_stokes_vec_get_components_view (&vel_vec, NULL,
                                            state->vel_press_vec);
      vel_elem = ymir_velocity_elem_new (mesh->ma->N, mesh->ma->ompsize);
      ymir_second_invariant_vec (vel_vec, strain_frob, vel_elem);
      ymir_dvec_scale (2.0, strain_frob);
      ymir_dvec_sqrt (strain_frob, strain_frob);

      /* multiply in Frobenius-norm */
      ymir_dvec_multiply_in1 (strain_frob, nltensor);

      /* destroy */
      ymir_vec_destroy (vel_vec);
      ymir_vec_destroy (strain_frob);
      ymir_velocity_elem_destroy (vel_elem);
    }
    break;

  case SL_NL_SOLVER_PRIMALDUAL_STRESS:
    {
      ymir_vec_t         *viscosity =
                            nl_stokes->stokes_op->stress_op->viscosity;
      ymir_mesh_t        *mesh = nl_stokes->mesh;
      ymir_vec_t         *vel_vec;
      ymir_vec_t         *strain_frob = ymir_dvec_new (mesh, 1,
                                                       YMIR_GAUSS_NODE);
      ymir_velocity_elem_t  *vel_elem;

      /* set nonlinear tensor to stress tensor */
      ymir_dvec_copy (dual, nltensor);

      /* divide in viscosity */
      ymir_dvec_divide_in1 (viscosity, nltensor);

      /* compute Frobenius-norm of the strain rate tensor */
      slabs_stokes_vec_get_components_view (&vel_vec, NULL,
                                            state->vel_press_vec);
      vel_elem = ymir_velocity_elem_new (mesh->cnodes->N, mesh->ma->ompsize);
      ymir_second_invariant_vec (vel_vec, strain_frob, vel_elem);
      ymir_dvec_scale (2.0, strain_frob);
      ymir_dvec_sqrt (strain_frob, strain_frob);

      /* scale tensor */
      switch (nl_solver_primaldual_scal_type) {
      case SL_NL_SOLVER_PRIMALDUAL_SCAL_NONE:
#ifdef YMIR_DEBUG
        n_scaled_nodes = 0;
#endif
        break;

      case SL_NL_SOLVER_PRIMALDUAL_SCAL_NORMALIZE:
#ifdef YMIR_DEBUG
        n_scaled_nodes =
#endif
        slabs_norm_symtensor_scale_frobenius (nltensor,
                                                               strain_frob);
        break;

      case SL_NL_SOLVER_PRIMALDUAL_SCAL_THRESHOLD:
#ifdef YMIR_DEBUG
        n_scaled_nodes =
#endif
        slabs_norm_symtensor_scale_threshold_frobenius (nltensor, strain_frob,
                                                        1.0);
        break;

      default: /* unknown primal-dual scaling type */
        YMIR_ABORT_NOT_REACHED ();
      }

      /* destroy */
      ymir_vec_destroy (vel_vec);
      ymir_vec_destroy (strain_frob);
      ymir_velocity_elem_destroy (vel_elem);
    }
    break;

  default: /* unknown primal-dual type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* print number of scaled nodes */
#ifdef YMIR_DEBUG
  {
    p4est_gloidx_t      n_elements = state->p8est->global_num_quadrants;
    int                 n_nodes_per_el = nl_stokes->mesh->ma->Np;
    double              rel_n_scaled_nodes;

    rel_n_scaled_nodes =   ((double) n_scaled_nodes)
                         / ((double) n_elements)
                         / ((double) n_nodes_per_el);

    YMIR_GLOBAL_VERBOSEF ("%s: %lli nodes (%.2f %%) scaled\n",
                          this_fn_name, (long long int) n_scaled_nodes,
                          100.0 * rel_n_scaled_nodes);
  }
#endif

  /* replace tensor in viscous stress operator */
  ymir_stress_op_set_nltens_rank2 (stress_op, nltensor);

  /* destroy */
  ymir_vec_destroy (nltensor);
}

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
            const int check_derivative)
{
  const char         *this_fn_name = "slabs_nonlinear_stokes_op_init";

  ymir_mesh_t        *mesh = nl_stokes->mesh;
  ymir_pressure_elem_t  *press_elem = nl_stokes->press_elem;
  ymir_vel_dir_t     *vel_dir = nl_stokes->vel_dir;

  ymir_vec_t         *viscosity = nl_stokes->viscosity;
  ymir_vec_t         *dvisc_dIIe, *rank1_tensor_scal;
  ymir_stokes_op_t   *stokes_op;

  YMIR_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* check input parameters */
  YMIR_ASSERT (nl_stokes->stokes_op == NULL);
  YMIR_ASSERT (nl_stokes->stokes_pc == NULL);

  /* create variables for the derivative of the nonlinear Stokes operator */
  dvisc_dIIe = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
  rank1_tensor_scal = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);

  /* compute physical viscosity and it's derivative
   * (updates `state->vel_bc_vec`) */
  if (init) {
    slabs_physics_compute_init_nl_stokes_coeff (
        viscosity, dvisc_dIIe, rank1_tensor_scal, nl_stokes->bounds_marker,
        nl_stokes->yielding_marker, state, press_elem, vel_dir,
        physics_options);
  }
  else {
    slabs_physics_compute_stokes_coeff (
        viscosity, dvisc_dIIe, rank1_tensor_scal, nl_stokes->bounds_marker,
        nl_stokes->yielding_marker, state, press_elem, vel_dir,
        physics_options);
  }

  /* set derivative of viscosity to zero if Picard method is used */
  if (nl_solver_type == SL_NL_SOLVER_PICARD) {
    ymir_vec_set_zero (dvisc_dIIe);
    ymir_vec_set_zero (rank1_tensor_scal);
  }

  /* update scaling of Dirichlet BC's */
  if (physics_options->bc_default_dirichlet_scale) {
    YMIR_ASSERT (vel_dir->scale != NULL);
    ymir_velocity_bc_default_dirscal (vel_dir->scale, viscosity);
    slabs_cvec_set_min_value (vel_dir->scale, SL_DIRSCAL_MIN);
    YMIR_ASSERT (sc_dmatrix_is_valid (vel_dir->scale->data));
  }

  /* create Stokes operator */
  stokes_op = slabs_stokes_op_new (viscosity, vel_dir, NULL, NULL, NULL,
                                   press_elem, physics_options);
  nl_stokes->stokes_op = stokes_op;

  /* update 4th order tensor */
  //TODO setting usol is deprecated
  ymir_stress_op_set_usol (stokes_op->stress_op, state->vel_bc_vec, dvisc_dIIe);
  ymir_stress_op_coeff_compute_rank1_tensor (stokes_op->stress_op,
                                             rank1_tensor_scal,
                                             state->vel_bc_vec);

  /* replace tensor for primal-dual method */
  slabs_nonlinear_stokes_primaldual_nltensor (nl_stokes, state,
                                              primaldual_type,
                                              primaldual_scal_type);

  /* test derivative of Stokes operator */
  if (check_derivative && nl_solver_type != SL_NL_SOLVER_PICARD) {
    slabs_nonlinear_stokes_problem_test_current_deriv (
        stokes_op, state, mesh, press_elem, vel_dir, viscosity,
        physics_options);
  }

  YMIR_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Destroys Stokes operator and preconditioner of nonlinear Stokes problem.
 */
static void
slabs_nonlinear_stokes_problem_op_clear (slabs_nl_stokes_problem_t *nl_stokes)
{
  const char         *this_fn_name = "slabs_nonlinear_stokes_problem_op_clear";

  /* destroy Stokes operator */
  if (nl_stokes->stokes_op != NULL) {
    ymir_stress_op_t   *stress_op = nl_stokes->stokes_op->stress_op;
    ymir_vec_t         *dvisc_dIIe = stress_op->dvdIIe;
    ymir_vec_t         *rank1_tensor_scal = stress_op->coeff_rank1_tensor_scal;

    YMIR_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

    /* destroy (assume viscosity in Stokes operator is not destroyed) */
    slabs_stokes_op_destroy (nl_stokes->stokes_op);
    nl_stokes->stokes_op = NULL;
    ymir_vec_destroy (dvisc_dIIe);
    ymir_vec_destroy (rank1_tensor_scal);

    YMIR_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
  }
}

/**
 * Updates Stokes operator for nonlinear Stokes problem.
 */
void
slabs_nonlinear_stokes_op_update (slabs_nl_stokes_problem_t *nl_stokes,
                                  slabs_stokes_state_t *state,
                                  slabs_physics_options_t *physics_options,
                                  const slabs_nl_solver_type_t
                                    nl_solver_type)
{
  const char         *this_fn_name = "slabs_nonlinear_stokes_op_update";
  ymir_stress_op_t   *stress_op = nl_stokes->stokes_op->stress_op;
  ymir_vel_dir_t     *vel_dir = nl_stokes->vel_dir;
  ymir_vec_t         *viscosity = stress_op->viscosity;
  ymir_vec_t         *dvisc_dIIe = stress_op->dvdIIe;
  ymir_vec_t         *rank1_tensor_scal = stress_op->coeff_rank1_tensor_scal;

  YMIR_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  /* check input */
  YMIR_ASSERT (nl_stokes->stokes_op != NULL);

  /* compute physical viscosity and its derivative
   * (updates `state->vel_bc_vec`) */
  slabs_physics_compute_stokes_coeff (viscosity, dvisc_dIIe, rank1_tensor_scal,
                                      nl_stokes->bounds_marker,
                                      nl_stokes->yielding_marker,
                                      state, nl_stokes->press_elem, vel_dir,
                                      physics_options);
  ymir_stress_op_setup_geo_coeff (stress_op);

  /* set derivative of viscosity to zero if Picard method is used */
  if (nl_solver_type == SL_NL_SOLVER_PICARD) {
    slabs_nonlinear_stokes_op_switch_picard (nl_stokes);
  }

  /* update scaling of Dirichlet BC's */
  if (physics_options->bc_default_dirichlet_scale) {
    YMIR_ASSERT (vel_dir->scale != NULL);
    ymir_velocity_bc_default_dirscal (vel_dir->scale, viscosity);
    slabs_cvec_set_min_value (vel_dir->scale, SL_DIRSCAL_MIN);
    YMIR_ASSERT (sc_dmatrix_is_valid (vel_dir->scale->cvec));
  }

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

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
            const int check_derivative)
{
  const char         *this_fn_name = "slabs_nonlinear_stokes_op_update_tensor";
  ymir_stress_op_t   *stress_op = nl_stokes->stokes_op->stress_op;

  YMIR_GLOBAL_INFOF ("%s: Update 4th order tensor\n", this_fn_name);

  /* check input */
  YMIR_ASSERT_HAS_CVEC (state->vel_bc_vec);

  /* update 4th order tensor */
  //TODO setting usol is deprecated
  ymir_stress_op_set_usol (stress_op, state->vel_bc_vec, stress_op->dvdIIe);
  ymir_stress_op_coeff_compute_rank1_tensor (stress_op,
                                             stress_op->coeff_rank1_tensor_scal,
                                             state->vel_bc_vec);

  /* replace tensor for primal-dual method */
  slabs_nonlinear_stokes_primaldual_nltensor (nl_stokes, state,
                                              primaldual_type,
                                              primaldual_scal_type);

  /* test derivative of Stokes operator */
  if (check_derivative) {
    slabs_nonlinear_stokes_problem_test_current_deriv (
        nl_stokes->stokes_op, state, nl_stokes->mesh, nl_stokes->press_elem,
        nl_stokes->vel_dir, nl_stokes->viscosity, physics_options);
  }
}

/**
 * Modifies nonlinear Stokes problem for Picard method.
 */
void
slabs_nonlinear_stokes_op_switch_picard (slabs_nl_stokes_problem_t *nl_stokes)
{
  ymir_vec_set_zero (nl_stokes->stokes_op->stress_op->dvdIIe);
  ymir_vec_set_zero (nl_stokes->stokes_op->stress_op->coeff_rank1_tensor_scal);
}

/**
 * Creates Stokes preconditioner for nonlinear Stokes problem.
 */
void
slabs_nonlinear_stokes_pc_init (
                    slabs_nl_stokes_problem_t *nl_stokes,
                    slabs_stokes_state_t *state,
                    const slabs_nl_stokes_problem_schur_diag_t schur_diag_type,
                    slabs_nl_stokes_problem_scaling_t scaling_type,
                    const int init)
{
  const char         *this_fn_name = "slabs_nonlinear_stokes_pc_init";
  ymir_stokes_op_t   *stokes_op = nl_stokes->stokes_op;

  YMIR_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* check input parameters */
  YMIR_ASSERT (nl_stokes->stokes_op != NULL);
  YMIR_ASSERT (nl_stokes->stokes_pc == NULL);

  /* set diagonal pressure Schur complement approximation */
  slabs_nl_stokes_problem_set_schur_diag (nl_stokes, state, schur_diag_type);

  /* switch off scaling for initial nonlinear step */
  if ( init &&
       (scaling_type == SL_NL_STOKES_PROB_SCALING_UUMASS_PPMASS ||
        scaling_type == SL_NL_STOKES_PROB_SCALING_AU_BP) ) {
    scaling_type = SL_NL_STOKES_PROB_SCALING_NONE;
  }

  /* set diagonal scaling */
  if (scaling_type != SL_NL_STOKES_PROB_SCALING_NONE) {
    slabs_nl_stokes_problem_set_scaling (nl_stokes, state, scaling_type);
  }

  /* create Stokes preconditioner */
  nl_stokes->stokes_pc = ymir_nlstokes_pc_new (stokes_op);

  YMIR_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Destroys Stokes preconditioner of nonlinear Stokes problem.
 */
static void
slabs_nonlinear_stokes_problem_pc_clear (slabs_nl_stokes_problem_t *nl_stokes)
{
  const char         *this_fn_name = "slabs_nonlinear_stokes_problem_pc_clear";

  YMIR_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* clear diagonal Schur complement variables */
  slabs_nl_stokes_problem_clear_schur_diag (nl_stokes);

  /* clear scaling variables */
  slabs_nl_stokes_problem_clear_scaling ();

  /* destroy Stokes preconditioner */
  if (nl_stokes->stokes_pc != NULL) {
    ymir_stokes_pc_destroy (nl_stokes->stokes_pc);
    nl_stokes->stokes_pc = NULL;
  }

  YMIR_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Updates Stokes preconditioner for nonlinear Stokes problem.
 */
ymir_stokes_pc_t *
slabs_nonlinear_stokes_pc_update (
                    slabs_nl_stokes_problem_t *nl_stokes,
                    slabs_stokes_state_t *state,
                    const slabs_nl_stokes_problem_schur_diag_t schur_diag_type,
                    const slabs_nl_stokes_problem_scaling_t scaling_type)
{
  const char         *this_fn_name = "slabs_nonlinear_stokes_pc_update";

  YMIR_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  /* check input parameters */
  YMIR_ASSERT (nl_stokes->stokes_pc != NULL);

#ifdef YMIR_PETSC
  /* clean up diagonal Schur complement variables */
  slabs_nl_stokes_problem_clear_schur_diag (nl_stokes);

  /* set diagonal pressure Schur complement approximation */
  slabs_nl_stokes_problem_set_schur_diag (nl_stokes, state, schur_diag_type);

  /* clean up scaling variables */
  slabs_nl_stokes_problem_clear_scaling ();

  /* set diagonal scaling */
  if (scaling_type != SL_NL_STOKES_PROB_SCALING_NONE) {
    slabs_nl_stokes_problem_set_scaling (nl_stokes, state, scaling_type);
  }

  /* update Stokes preconditioner */
  ymir_nlstokes_pc_recompute (nl_stokes->stokes_pc);
#else
  YMIR_ABORT_NOT_REACHED ();
#endif

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);

  /* return new Stokes preconditioner */
  return nl_stokes->stokes_pc;
}

/******************************************************************************
 * Tests
 *****************************************************************************/

/**
 * Tests Jacobian of nonlinear Stokes operator by checking the equality:
 *
 *   [D F](x) * y = [d F](x;y)
 *
 * where
 *
 *   `F` is Stokes operator
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
    ymir_vec_t  *viscosity,
    slabs_physics_options_t *physics_options)
{
  const char         *this_fn_name =
                        "slabs_nonlinear_stokes_problem_test_current_deriv";

  const ymir_stress_op_rot_project_t  stokes_op_project_rot =
    stokes_op->project_out_rot;
  const ymir_stokes_op_press_project_t  stokes_op_project_mean =
    stokes_op->project_out_mean;
  ymir_vec_t         *stokes_op_upscale = stokes_op->upscale;
  slabs_stokes_state_t  *test_state;
  ymir_vec_t         *test_visc = ymir_vec_template (viscosity);
  ymir_vec_t         *test_dvisc_dIIe = ymir_vec_template (viscosity);
  ymir_stokes_op_t   *test_stokes_op;
  ymir_stress_op_rot_project_t  test_stokes_op_project_rot;
  ymir_stokes_op_press_project_t  test_stokes_op_project_mean;
  ymir_vec_t         *test_stokes_op_upscale;
  ymir_vec_t         *x = ymir_vec_template (state->vel_press_vec);
  ymir_vec_t         *y = ymir_vec_template (state->vel_press_vec);
  ymir_vec_t         *y_vel;
  ymir_vec_t         *y_press;
  ymir_vec_t         *x_hy = ymir_vec_template (state->vel_press_vec);
  ymir_vec_t         *Fx = ymir_vec_template (state->vel_press_vec);
  ymir_vec_t         *gateaux_vec = ymir_vec_template (state->vel_press_vec);
  ymir_vec_t         *frechet_vec = ymir_vec_template (state->vel_press_vec);
  const int           n_tests = 6;
  const double        h_init = 1.0e-0;
  const double        h_step = 1.0e-2;
  double              h;
  int                 testid;
  double              abs_error, frechet_norm;

  YMIR_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  /* copy Stokes state */
  test_state = slabs_stokes_state_new (state->p8est);
  slabs_stokes_state_init_temp (test_state, mesh->cnodes);
  slabs_stokes_state_init_temp_vec (test_state, mesh);
  slabs_stokes_state_init_weakzone (test_state, mesh,
                                    slabs_physics_compute_weakzone,
                                    physics_options);
  ymir_vec_copy (state->weak_vec, test_state->weak_vec);
  slabs_stokes_state_init_vel_press (test_state, mesh, press_elem);
  sc_dmatrix_copy (state->temperature, test_state->temperature);
  ymir_vec_copy (state->vel_press_vec, test_state->vel_press_vec);
  ymir_vec_copy (state->vel_bc_vec, test_state->vel_bc_vec);

  /* create Stokes operator with nonlinear viscosity */
  ymir_vec_copy (viscosity, test_visc);
  test_stokes_op = slabs_stokes_op_new (test_visc, vel_dir, NULL, NULL, NULL,
                                        press_elem, physics_options);

  /* deactivate projection and scaling of Stokes operators */
  test_stokes_op_project_rot = test_stokes_op->project_out_rot;
  test_stokes_op_project_mean = test_stokes_op->project_out_mean;
  test_stokes_op_upscale = test_stokes_op->upscale;
  test_stokes_op->project_out_rot = YMIR_STRESS_OP_ROT_PROJECT_NONE;
  test_stokes_op->project_out_mean = YMIR_STOKES_OP_PRESS_PROJECT_NONE;
  test_stokes_op->upscale = NULL;
  stokes_op->project_out_rot = YMIR_STRESS_OP_ROT_PROJECT_NONE;
  stokes_op->project_out_mean = YMIR_STOKES_OP_PRESS_PROJECT_NONE;
  stokes_op->upscale = NULL;

  /* create vector views */
  slabs_stokes_vec_get_components_view (&y_vel, &y_press, y);

  /*
   * Test derivative of Stokes operator
   */

  /* compute random vector `y` */
#ifdef YMIR_PETSC
  ymir_petsc_vec_set_random (y, YMIR_MESH_PETSCLAYOUT_NONE);
#else
  YMIR_ABORT_NOT_REACHED ();
#endif

  /* apply random vector `y` to Frechet derivative */
  ymir_nlstokes_op_apply (y, frechet_vec, stokes_op);
  frechet_norm = ymir_vec_norm (frechet_vec);

  /* apply vector `x`, containing current state, to get `F(x)` */
  ymir_vec_copy (state->vel_press_vec, x);
  ymir_stokes_op_apply (x, Fx, test_stokes_op);

  h = h_init;
  for (testid = 0; testid < n_tests; testid++) {
    /* compute new position `x + h*y` */
    ymir_vec_copy (x, x_hy);
    ymir_vec_add (h, y, x_hy);

    /* compute viscosity */
    ymir_vec_copy (x_hy, test_state->vel_press_vec);
    slabs_physics_compute_stokes_coeff (test_visc, test_dvisc_dIIe,
                                        NULL, NULL, NULL,
                                        test_state, press_elem, vel_dir,
                                        physics_options);

    /* apply Stokes operator to get `F(x + h*y)` */
    ymir_stokes_op_apply (x_hy, gateaux_vec, test_stokes_op);

    /* compute Gateaux derivative */
    ymir_vec_add (-1.0, Fx, gateaux_vec);
    ymir_vec_scale (1.0 / h, gateaux_vec);

    /* check error */
    ymir_vec_add (-1.0, frechet_vec, gateaux_vec);
    abs_error = ymir_vec_norm (gateaux_vec);

    /* print error */
    if (0.0 < frechet_norm) {
      YMIR_GLOBAL_VERBOSEF ("%s: h %1.1e, abs error %1.3e, rel error %1.3e\n",
                            this_fn_name, h, abs_error,
                            abs_error / frechet_norm);
    }
    else {
      YMIR_GLOBAL_VERBOSEF ("%s: h %1.1e, abs error %1.3e\n",
                            this_fn_name, h, abs_error);
    }

    /* reduce Gateaux step length */
    h *= h_step;
  }

  /*
   * Clean up
   */

  /* restore projection and scaling of Stokes operators */
  test_stokes_op->project_out_rot = test_stokes_op_project_rot;
  test_stokes_op->project_out_mean = test_stokes_op_project_mean;
  test_stokes_op->upscale = test_stokes_op_upscale;
  stokes_op->project_out_rot = stokes_op_project_rot;
  stokes_op->project_out_mean = stokes_op_project_mean;
  stokes_op->upscale = stokes_op_upscale;

  /* destroy */
  slabs_stokes_op_destroy (test_stokes_op);
  slabs_stokes_state_destroy (test_state);
  ymir_vec_destroy (test_visc);
  ymir_vec_destroy (test_dvisc_dIIe);
  ymir_vec_destroy (x);
  ymir_vec_destroy (y);
  ymir_vec_destroy (y_vel);
  ymir_vec_destroy (y_press);
  ymir_vec_destroy (x_hy);
  ymir_vec_destroy (Fx);
  ymir_vec_destroy (gateaux_vec);
  ymir_vec_destroy (frechet_vec);

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

/* structure for computations with 4th order viscosity tensor & pressure mass */
typedef struct slabs_nl_stokes_problem_visc_tensor_mass_data
{
  ymir_vec_t         *velocity;
  ymir_vec_t         *viscosity;
  ymir_vec_t         *dvisc_dIIe;
  ymir_vec_t         *IIe;
  ymir_mesh_t        *mesh;
  ymir_pressure_elem_t  *press_elem;
}
slabs_nl_stokes_problem_visc_tensor_mass_data_t;

/**
 * Computes the Frobenius norm of the 4th order viscosity tensor at each node.
 *
 * Note: In order for the Frobenius norm of the identity to be one, we scale
 * the norm by 1/3.
 */
static void
slabs_nl_stokes_problem_visc_tensor_norm (double *norm, double x, double y,
                                          double z, ymir_locidx_t nid,
                                          void *data)
{
  slabs_nl_stokes_problem_visc_tensor_mass_data_t  *d =
    (slabs_nl_stokes_problem_visc_tensor_mass_data_t *) data;
  double             *viscosity_data = d->viscosity->dvec->e[0];
  double             *dvisc_dIIe_data = d->dvisc_dIIe->dvec->e[0];
  double             *IIe_data;

  YMIR_ASSERT (d->IIe != NULL);
  IIe_data = d->IIe->dvec->e[0];

  *norm = sqrt (
        9.0 * SC_SQR (viscosity_data[nid])
      + 4.0 * viscosity_data[nid] * dvisc_dIIe_data[nid] * IIe_data[nid]
      + 4.0 * SC_SQR (dvisc_dIIe_data[nid]) * SC_SQR (IIe_data[nid])
  ) / 3.0;
}

/**
 * Computes diagonal pressure Schur complement approximation with inverse
 * viscosity tensor weighted lumped mass matrix.
 */
void
slabs_nl_stokes_problem_schur_diag_visc_tensor_mass (ymir_vec_t *schur_diag,
                                                     void *data)
{
  slabs_nl_stokes_problem_visc_tensor_mass_data_t  *d =
    (slabs_nl_stokes_problem_visc_tensor_mass_data_t *) data;
  ymir_vec_t         *velocity = d->velocity;
  ymir_mesh_t        *mesh = d->mesh;
  ymir_pressure_elem_t  *press_elem = d->press_elem;
  ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (mesh->ma->N,
                                                            mesh->ma->ompsize);
  ymir_vec_t         *IIe = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
  ymir_vec_t         *visc_tensor_norm = ymir_dvec_new (mesh, 1,
                                                        YMIR_GAUSS_NODE);

  d->IIe = IIe;
  ymir_second_invariant_vec (velocity, IIe, vel_elem);
  ymir_dvec_set_function (visc_tensor_norm,
                          slabs_nl_stokes_problem_visc_tensor_norm, data);

  ymir_pressure_vec_lump_weighted_mass (schur_diag, visc_tensor_norm, -1.0,
                                        press_elem);
  ymir_evec_fabs (schur_diag, schur_diag);
  ymir_evec_reciprocal (schur_diag);

  ymir_velocity_elem_destroy (vel_elem);
  ymir_vec_destroy (visc_tensor_norm);
  ymir_vec_destroy (IIe);
  d->IIe = NULL;
}

/**
 * Sets the function used for computing the diagonal Schur complement approx.
 * and the corresponding function data.
 */
static void
slabs_nl_stokes_problem_set_schur_diag (slabs_nl_stokes_problem_t *nl_stokes,
    slabs_stokes_state_t *state,
    const slabs_nl_stokes_problem_schur_diag_t schur_diag_type)
{
  switch (schur_diag_type) {
  case SL_NL_STOKES_PROB_SCHUR_DIAG_INV_VISC_PMASS:
    /* will be set by ymir_stress_pc_new/recompute */
    break;

  case SL_NL_STOKES_PROB_SCHUR_DIAG_INV_VISC_TENSOR_PMASS:
    {
      ymir_vec_t         *vel_view;
      slabs_nl_stokes_problem_visc_tensor_mass_data_t  *data;

      /* get view of velocity and pressure components */
      slabs_stokes_vec_get_components_view (&vel_view, NULL,
                                            state->vel_press_vec);

      /* set diagonal Schur complement approx */
      data = YMIR_ALLOC (slabs_nl_stokes_problem_visc_tensor_mass_data_t, 1);
      data->velocity = vel_view;
      data->viscosity = nl_stokes->stokes_op->stress_op->viscosity;
      data->dvisc_dIIe = nl_stokes->stokes_op->stress_op->dvdIIe;
      data->IIe = NULL;
      data->mesh = nl_stokes->mesh;
      data->press_elem = nl_stokes->press_elem;

      ymir_stokes_pc_set_schur_diag_fn =
        slabs_nl_stokes_problem_schur_diag_visc_tensor_mass;
      ymir_stokes_pc_set_schur_diag_fn_data = data;
    }
    break;

  default: /* unknown diagonal Schur type */
    YMIR_ABORT_NOT_REACHED ();
  }
}

/**
 * Destroys data that was used for diagonal Schur complement approximation.
 */
static void
slabs_nl_stokes_problem_clear_schur_diag (slabs_nl_stokes_problem_t *nl_stokes)
{
  if (   ymir_stokes_pc_set_schur_diag_fn
      == slabs_nl_stokes_problem_schur_diag_visc_tensor_mass ) {
    slabs_nl_stokes_problem_visc_tensor_mass_data_t  *d =
      (slabs_nl_stokes_problem_visc_tensor_mass_data_t *)
      ymir_stokes_pc_set_schur_diag_fn_data;

    ymir_vec_destroy (d->velocity);
    YMIR_FREE (ymir_stokes_pc_set_schur_diag_fn_data);
  }
}

/*****************************************************************************
 * scaling stuff TODO deprecated
 ****************************************************************************/

/**
 * Computes velocity scaling with inverse lumped mass matrix.
 */
void
slabs_nl_stokes_problem_uscale_mass (ymir_vec_t *uscale, void *data)
{
  ymir_mass_lump (uscale);
  ymir_cvec_pow (-1.0 / 2.0, uscale);
}

/**
 * Computes velocity scaling with inverse viscosity weighted lumped mass matrix.
 */
void
slabs_nl_stokes_problem_uscale_visc_mass (ymir_vec_t *uscale, void *data)
{
  ymir_vec_t         *viscosity = (ymir_vec_t *) data;

  /* compute velocity scaling */
  ymir_mass_lump_weighted (uscale, viscosity, 1.0);
  ymir_cvec_fabs (uscale, uscale);
  ymir_cvec_pow (-1.0 / 2.0, uscale);
}

/**
 * Computes velocity scaling with quotient grad viscosity/viscosity.
 */
void
slabs_nl_stokes_problem_uscale_grad_visc (ymir_vec_t *uscale, void *data)
{
  ymir_vec_t         *viscosity = (ymir_vec_t *) data;
  ymir_vec_t         *grad_visc = ymir_dvec_new (viscosity->mesh, 3,
                                                 YMIR_GAUSS_NODE);
  const sc_bint_t     m = grad_visc->dvec->m;
  const sc_bint_t     n = grad_visc->dvec->n;
  sc_dmatrix_t       *visc_mat = sc_dmatrix_clone (viscosity->dvec);

  /* compute gradient of viscosity (let's accept the error that's produced
   * because the viscosity is defined on Gauss nodes and not GLL nodes TODO) */
  sc_dmatrix_reshape (grad_visc->dvec, viscosity->K * viscosity->Np, 3);
  sc_dmatrix_reshape (visc_mat, viscosity->K * viscosity->Np, 1);
  slabs_gradient_gll_to_gauss (visc_mat, grad_visc->dvec, viscosity->mesh->ma);
  sc_dmatrix_reshape (grad_visc->dvec, m, n);

  /* compute velocity scaling */
  ymir_dvec_fabs (grad_visc, grad_visc);
  slabs_matrix_bound_values (grad_visc->dvec, 0.1, DBL_MAX); // TODO
  ymir_dvec_divide_in1 (viscosity, grad_visc);
  ymir_interp_vec (grad_visc, uscale);
  ymir_cvec_fabs (uscale, uscale);
  ymir_cvec_sqrt (uscale, uscale);

  /* destroy */
  sc_dmatrix_destroy (visc_mat);
  ymir_vec_destroy (grad_visc);
}

/**
 * Computes velocity scaling with velocity weighted lumped mass matrix.
 */
typedef struct slabs_nl_stokes_problem_uscale_vel_mass_data
{
  ymir_vec_t         *velocity;
  ymir_mesh_t        *mesh;
}
slabs_nl_stokes_problem_uscale_vel_mass_data_t;

void
slabs_nl_stokes_problem_uscale_vel_mass (ymir_vec_t *uscale, void *data)
{
  slabs_nl_stokes_problem_uscale_vel_mass_data_t  *d =
    (slabs_nl_stokes_problem_uscale_vel_mass_data_t *) data;
  ymir_vec_t         *velocity = d->velocity;
  ymir_mesh_t        *mesh = d->mesh;
  ymir_vec_t         *vel_abs = ymir_cvec_new (mesh, 3);

  ymir_mass_lump (uscale);
  ymir_cvec_fabs (velocity, vel_abs);
  slabs_matrix_bound_values (vel_abs->cvec, 0.1, DBL_MAX); // TODO
  ymir_cvec_multiply_in (vel_abs, uscale);
  ymir_cvec_pow (-1.0 / 2.0, uscale);

  ymir_vec_destroy (vel_abs);
}

/**
 * Computes velocity scaling by applying stress operator to velocity.
 */
typedef struct slabs_nl_stokes_problem_uscale_Au_data
{
  ymir_vec_t         *velocity;
  ymir_stress_op_t   *stress_op;
}
slabs_nl_stokes_problem_uscale_Au_data_t;

void
slabs_nl_stokes_problem_uscale_Au (ymir_vec_t *uscale, void *data)
{
  slabs_nl_stokes_problem_uscale_Au_data_t  *d =
    (slabs_nl_stokes_problem_uscale_Au_data_t *) data;
  ymir_vec_t         *velocity = d->velocity;
  ymir_stress_op_t   *stress_op = d->stress_op;

  ymir_nlstress_op_apply (velocity, uscale, stress_op);
  ymir_cvec_fabs (uscale, uscale);
  slabs_matrix_bound_values (uscale->cvec, 0.1, DBL_MAX); // TODO
  ymir_cvec_pow (-1.0 / 2.0, uscale);
}

/**
 * Computes pressure scaling with inverse lumped mass matrix.
 */
void
slabs_nl_stokes_problem_pscale_mass (ymir_vec_t *pscale, void *data)
{
  ymir_pressure_elem_t  *press_elem = (ymir_pressure_elem_t *) data;

  ymir_pressure_vec_lump_mass (pscale, press_elem);
  ymir_evec_pow (-1.0 / 2.0, pscale);
}

/**
 * Computes pressure scaling with inverse viscosity weighted lumped mass matrix.
 */
typedef struct slabs_nl_stokes_problem_pscale_visc_mass_data
{
  ymir_vec_t         *viscosity;
  ymir_pressure_elem_t  *press_elem;
}
slabs_nl_stokes_problem_pscale_visc_mass_data_t;

void
slabs_nl_stokes_problem_pscale_visc_mass (ymir_vec_t *pscale, void *data)
{
  slabs_nl_stokes_problem_pscale_visc_mass_data_t  *d =
    (slabs_nl_stokes_problem_pscale_visc_mass_data_t *) data;
  ymir_vec_t         *viscosity = d->viscosity;
  ymir_pressure_elem_t  *press_elem = d->press_elem;

  ymir_pressure_vec_lump_weighted_mass (pscale, viscosity, 1.0, press_elem);
  ymir_evec_fabs (pscale, pscale);
  ymir_evec_pow (-1.0 / 2.0, pscale);
}

/**
 * Computes pressure scaling with inverse viscosity tensor weighted lumped mass
 * matrix.
 */
void
slabs_nl_stokes_problem_pscale_visc_tensor_mass (ymir_vec_t *pscale,
                                                 void *data)
{
  slabs_nl_stokes_problem_visc_tensor_mass_data_t  *d =
    (slabs_nl_stokes_problem_visc_tensor_mass_data_t *) data;
  ymir_vec_t         *velocity = d->velocity;
  ymir_mesh_t        *mesh = d->mesh;
  ymir_pressure_elem_t  *press_elem = d->press_elem;
  ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (mesh->ma->N,
                                                            mesh->ma->ompsize);
  ymir_vec_t         *IIe = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
  ymir_vec_t         *visc_tensor_norm = ymir_dvec_new (mesh, 1,
                                                        YMIR_GAUSS_NODE);

  d->IIe = IIe;
  ymir_second_invariant_vec (velocity, IIe, vel_elem);
  ymir_dvec_set_function (visc_tensor_norm,
                          slabs_nl_stokes_problem_visc_tensor_norm, data);

  ymir_pressure_vec_lump_weighted_mass (pscale, visc_tensor_norm, 1.0,
                                        press_elem);
  ymir_evec_fabs (pscale, pscale);
  ymir_evec_pow (-1.0 / 2.0, pscale);

  ymir_velocity_elem_destroy (vel_elem);
  ymir_vec_destroy (visc_tensor_norm);
  ymir_vec_destroy (IIe);
  d->IIe = NULL;
}

/**
 * Computes pressure scaling with pressure weighted lumped mass matrix.
 */
typedef struct slabs_nl_stokes_problem_pscale_press_mass_data
{
  ymir_evec_t        *pressure;
  ymir_mesh_t        *mesh;
  ymir_pressure_elem_t  *press_elem;
}
slabs_nl_stokes_problem_pscale_press_mass_data_t;

void
slabs_nl_stokes_problem_pscale_press_mass (ymir_vec_t *pscale, void *data)
{
  slabs_nl_stokes_problem_pscale_press_mass_data_t  *d =
    (slabs_nl_stokes_problem_pscale_press_mass_data_t *) data;
  ymir_evec_t        *pressure = d->pressure;
  ymir_mesh_t        *mesh = d->mesh;
  ymir_pressure_elem_t  *press_elem = d->press_elem;
  ymir_vec_t         *press_abs = ymir_pressure_vec_new (mesh, press_elem);

  ymir_pressure_vec_lump_mass (pscale, press_elem);
  ymir_evec_fabs (pressure, press_abs);
  slabs_matrix_bound_values (press_abs->evec, 0.1, DBL_MAX); // TODO
  ymir_evec_multiply_in (press_abs, pscale);
  ymir_evec_pow (-1.0 / 2.0, pscale);

  ymir_vec_destroy (press_abs);
}

/**
 * Computes pressure scaling by applying gradient operator to pressure.
 */
typedef struct slabs_nl_stokes_problem_pscale_Bp_data
{
  ymir_evec_t        *pressure;
  ymir_stokes_op_t   *stokes_op;
  ymir_mesh_t        *mesh;
}
slabs_nl_stokes_problem_pscale_Bp_data_t;

void
slabs_nl_stokes_problem_pscale_Bp (ymir_vec_t *pscale, void *data)
{
  slabs_nl_stokes_problem_pscale_Bp_data_t  *d =
    (slabs_nl_stokes_problem_pscale_Bp_data_t *) data;
  ymir_evec_t        *pressure = d->pressure;
  ymir_stokes_op_t   *stokes_op = d->stokes_op;
  ymir_mesh_t        *mesh = d->mesh;
  ymir_vec_t         *Bp = ymir_cvec_new (mesh, 3);
  ymir_vec_t         *Bp_norm = ymir_cvec_new (mesh, 1);

  ymir_stokes_grad_block_apply (pressure, Bp, stokes_op);
  slabs_cvec_compute_magnitude (Bp, Bp_norm);
  ymir_interp_vec (Bp_norm, pscale);
  ymir_evec_fabs (pscale, pscale);
  slabs_matrix_bound_values (pscale->evec, 0.1, DBL_MAX); // TODO
  ymir_evec_pow (-1.0 / 2.0, pscale);

  ymir_vec_destroy (Bp);
  ymir_vec_destroy (Bp_norm);
}

/**
 * Sets the scaling function and the corresponding function data.
 */
static void
slabs_nl_stokes_problem_set_scaling (slabs_nl_stokes_problem_t *nl_stokes,
                                     slabs_stokes_state_t *state,
                                     const slabs_nl_stokes_problem_scaling_t
                                       scaling_type)
{
  YMIR_ASSERT (ymir_stokes_pc_default_scale == 1);

  switch (scaling_type) {
  case SL_NL_STOKES_PROB_SCALING_ADIAG_VISC_PMASS:
    {
      slabs_nl_stokes_problem_pscale_visc_mass_data_t  *p_data;

      /* velocity scaling will be set by ymir_stress_pc_new/recompute */

      /* set pressure scaling */
      p_data = YMIR_ALLOC (slabs_nl_stokes_problem_pscale_visc_mass_data_t, 1);
      p_data->viscosity = nl_stokes->stokes_op->stress_op->viscosity;
      p_data->press_elem = nl_stokes->press_elem;

      ymir_stokes_pc_set_pscale_fn = slabs_nl_stokes_problem_pscale_visc_mass;
      ymir_stokes_pc_set_pscale_fn_data = p_data;
    }
    break;

  case SL_NL_STOKES_PROB_SCALING_ADIAG_VISC_TENSOR_PMASS:
    {
      ymir_vec_t         *vel_view;
      slabs_nl_stokes_problem_visc_tensor_mass_data_t  *data;

      /* get view of velocity and pressure components */
      slabs_stokes_vec_get_components_view (&vel_view, NULL,
                                            state->vel_press_vec);

      /* velocity scaling will be set by ymir_stress_pc_new/recompute */

      /* set pressure scaling */
      data = YMIR_ALLOC (slabs_nl_stokes_problem_visc_tensor_mass_data_t, 1);
      data->velocity = vel_view;
      data->viscosity = nl_stokes->stokes_op->stress_op->viscosity;
      data->dvisc_dIIe = nl_stokes->stokes_op->stress_op->dvdIIe;
      data->IIe = NULL;
      data->mesh = nl_stokes->mesh;
      data->press_elem = nl_stokes->press_elem;

      ymir_stokes_pc_set_pscale_fn =
        slabs_nl_stokes_problem_pscale_visc_tensor_mass;
      ymir_stokes_pc_set_pscale_fn_data = data;
    }
    break;

  case SL_NL_STOKES_PROB_SCALING_ADIAG_PMASS:
    /* velocity scaling will be set by ymir_stress_pc_new/recompute */

    /* set pressure scaling */
    ymir_stokes_pc_set_pscale_fn = slabs_nl_stokes_problem_pscale_mass;
    ymir_stokes_pc_set_pscale_fn_data = nl_stokes->press_elem;
    break;

  case SL_NL_STOKES_PROB_SCALING_ADIAG_BDIAG: /* BFBT Schur PC only! */
    /* velocity scaling will be set by ymir_stress_pc_new/recompute */

    YMIR_ASSERT (ymir_stokes_pc_schur_type == YMIR_STOKES_PC_SCHUR_BFBT ||
                 ymir_stokes_pc_schur_type == YMIR_STOKES_PC_SCHUR_LAPL);
    /* pressure scaling will be set by ymir_stokes_pc_new/recompute */
    break;

  case SL_NL_STOKES_PROB_SCALING_UMASS_PMASS:
    /* set velocity scaling */
    ymir_stokes_pc_set_uscale_fn = slabs_nl_stokes_problem_uscale_mass;

    /* set pressure scaling */
    ymir_stokes_pc_set_pscale_fn = slabs_nl_stokes_problem_pscale_mass;
    ymir_stokes_pc_set_pscale_fn_data = nl_stokes->press_elem;
    break;

  case SL_NL_STOKES_PROB_SCALING_VISC_UMASS_PMASS:
    /* set velocity scaling */
    ymir_stokes_pc_set_uscale_fn = slabs_nl_stokes_problem_uscale_visc_mass;
    ymir_stokes_pc_set_uscale_fn_data =
      nl_stokes->stokes_op->stress_op->viscosity;

    /* set pressure scaling */
    ymir_stokes_pc_set_pscale_fn = slabs_nl_stokes_problem_pscale_mass;
    ymir_stokes_pc_set_pscale_fn_data = nl_stokes->press_elem;
    break;

  case SL_NL_STOKES_PROB_SCALING_VISC_UMASS_VISC_PMASS:
    {
      slabs_nl_stokes_problem_pscale_visc_mass_data_t  *p_data;

      /* set velocity scaling */
      ymir_stokes_pc_set_uscale_fn = slabs_nl_stokes_problem_uscale_visc_mass;
      ymir_stokes_pc_set_uscale_fn_data =
        nl_stokes->stokes_op->stress_op->viscosity;

      /* set pressure scaling */
      p_data = YMIR_ALLOC (slabs_nl_stokes_problem_pscale_visc_mass_data_t, 1);
      p_data->viscosity = nl_stokes->stokes_op->stress_op->viscosity;
      p_data->press_elem = nl_stokes->press_elem;

      ymir_stokes_pc_set_pscale_fn = slabs_nl_stokes_problem_pscale_visc_mass;
      ymir_stokes_pc_set_pscale_fn_data = p_data;
    }
    break;

  case SL_NL_STOKES_PROB_SCALING_GRAD_VISC_PMASS:
    /* set velocity scaling */
    ymir_stokes_pc_set_uscale_fn = slabs_nl_stokes_problem_uscale_grad_visc;
    ymir_stokes_pc_set_uscale_fn_data =
      nl_stokes->stokes_op->stress_op->viscosity;

    /* set pressure scaling */
    ymir_stokes_pc_set_pscale_fn = slabs_nl_stokes_problem_pscale_mass;
    ymir_stokes_pc_set_pscale_fn_data = nl_stokes->press_elem;
    break;

  case SL_NL_STOKES_PROB_SCALING_UUMASS_PPMASS:
    {
      ymir_vec_t         *vel_view;
      ymir_evec_t        *press_view;
      slabs_nl_stokes_problem_uscale_vel_mass_data_t  *u_data;
      slabs_nl_stokes_problem_pscale_press_mass_data_t  *p_data;

      /* get view of velocity and pressure components */
      slabs_stokes_vec_get_components_view (&vel_view, &press_view,
                                            state->vel_press_vec);

      /* set velocity scaling */
      u_data = YMIR_ALLOC (slabs_nl_stokes_problem_uscale_vel_mass_data_t, 1);
      u_data->velocity = vel_view;
      u_data->mesh = nl_stokes->mesh;

      ymir_stokes_pc_set_uscale_fn = slabs_nl_stokes_problem_uscale_vel_mass;
      ymir_stokes_pc_set_uscale_fn_data = u_data;

      /* set pressure scaling */
      p_data = YMIR_ALLOC (slabs_nl_stokes_problem_pscale_press_mass_data_t, 1);
      p_data->pressure = press_view;
      p_data->mesh = nl_stokes->mesh;
      p_data->press_elem = nl_stokes->press_elem;

      ymir_stokes_pc_set_pscale_fn = slabs_nl_stokes_problem_pscale_press_mass;
      ymir_stokes_pc_set_pscale_fn_data = p_data;
    }
    break;

  case SL_NL_STOKES_PROB_SCALING_AU_BP:
    {
      ymir_vec_t         *vel_view;
      ymir_evec_t        *press_view;
      slabs_nl_stokes_problem_uscale_Au_data_t  *u_data;
      slabs_nl_stokes_problem_pscale_Bp_data_t  *p_data;

      /* get view of velocity and pressure components */
      slabs_stokes_vec_get_components_view (&vel_view, &press_view,
                                            state->vel_press_vec);

      /* set velocity scaling */
      u_data = YMIR_ALLOC (slabs_nl_stokes_problem_uscale_Au_data_t, 1);
      u_data->velocity = vel_view;
      u_data->stress_op = nl_stokes->stokes_op->stress_op;

      ymir_stokes_pc_set_uscale_fn = slabs_nl_stokes_problem_uscale_Au;
      ymir_stokes_pc_set_uscale_fn_data = u_data;

      /* set pressure scaling */
      p_data = YMIR_ALLOC (slabs_nl_stokes_problem_pscale_Bp_data_t, 1);
      p_data->pressure = press_view;
      p_data->stokes_op = nl_stokes->stokes_op;
      p_data->mesh = nl_stokes->mesh;

      ymir_stokes_pc_set_pscale_fn = slabs_nl_stokes_problem_pscale_Bp;
      ymir_stokes_pc_set_pscale_fn_data = p_data;
    }
    break;

//TODO deprecated
#if 0
  case SL_NL_STOKES_PROB_SCALING_STRESS_PRESS:
    YMIR_ASSERT (state != NULL);
    YMIR_ASSERT (state->vel_press_vec != NULL);
    {
      ymir_vec_t         *viscosity =
                            stokes_pc->stokes_op->stress_op->viscosity;
      ymir_vec_t         *vel_view;
      ymir_vec_t         *press_view;
      ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (mesh->ma->N,
                                                            mesh->ma->ompsize);
      ymir_vec_t         *IIe = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
      ymir_vec_t         *stress_cont = ymir_cvec_new (mesh, 1);
      ymir_vec_t         *press_abs = ymir_pressure_vec_new (mesh, press_elem);

      /* get view of velocity and pressure components */
      slabs_stokes_vec_get_components_view (&vel_view, &press_view,
                                            state->vel_press_vec);

      /* compute velocity scaling^2 */
      ymir_second_invariant_vec (vel_view, IIe, vel_elem);
      ymir_dvec_sqrt (IIe, IIe);
      ymir_dvec_multiply_in (viscosity, IIe);
      ymir_interp_vec (IIe, stress_cont);
      ymir_cvec_fabs (stress_cont, stress_cont);

      ymir_mass_lump (uscale_sq);
      ymir_cvec_multiply_in1 (stress_cont, uscale_sq);
      ymir_cvec_reciprocal (uscale_sq);

      /* compute pressure scaling */
      ymir_pressure_vec_lump_mass (pscale, press_elem);
      ymir_evec_fabs (press_view, press_abs);
      ymir_evec_multiply_in (press_abs, pscale);
      ymir_evec_pow (-1.0 / 2.0, pscale);

      /* destroy */
      ymir_vec_destroy (vel_view);
      ymir_vec_destroy (press_view);
      ymir_velocity_elem_destroy (vel_elem);
      ymir_vec_destroy (IIe);
      ymir_vec_destroy (stress_cont);
      ymir_vec_destroy (press_abs);
    }
    break;
#endif

  default: /* unknown scaling type */
    YMIR_ABORT_NOT_REACHED ();
  }
}

/**
 * Destroys data that was used for scaling functions.
 */
static void
slabs_nl_stokes_problem_clear_scaling ()
{
  /* clear velocity scaling data */
  if (   ymir_stokes_pc_set_uscale_fn
      == slabs_nl_stokes_problem_uscale_vel_mass ) {
    slabs_nl_stokes_problem_uscale_vel_mass_data_t  *d =
      (slabs_nl_stokes_problem_uscale_vel_mass_data_t *)
      ymir_stokes_pc_set_uscale_fn_data;

    ymir_vec_destroy (d->velocity);
    YMIR_FREE (ymir_stokes_pc_set_uscale_fn_data);
  }
  else if (   ymir_stokes_pc_set_uscale_fn
           == slabs_nl_stokes_problem_uscale_Au ) {
    slabs_nl_stokes_problem_uscale_Au_data_t  *d =
      (slabs_nl_stokes_problem_uscale_Au_data_t *)
      ymir_stokes_pc_set_uscale_fn_data;

    ymir_vec_destroy (d->velocity);
    YMIR_FREE (ymir_stokes_pc_set_uscale_fn_data);
  }

  /* clear pressure scaling data */
  if (   ymir_stokes_pc_set_pscale_fn
      == slabs_nl_stokes_problem_pscale_visc_mass ) {
    YMIR_FREE (ymir_stokes_pc_set_pscale_fn_data);
  }
  else if (   ymir_stokes_pc_set_pscale_fn
           == slabs_nl_stokes_problem_pscale_visc_tensor_mass ) {
    slabs_nl_stokes_problem_visc_tensor_mass_data_t  *d =
      (slabs_nl_stokes_problem_visc_tensor_mass_data_t *)
      ymir_stokes_pc_set_pscale_fn_data;

    ymir_vec_destroy (d->velocity);
    YMIR_FREE (ymir_stokes_pc_set_pscale_fn_data);
  }
  else if (   ymir_stokes_pc_set_pscale_fn
           == slabs_nl_stokes_problem_pscale_press_mass ) {
    slabs_nl_stokes_problem_pscale_press_mass_data_t  *d =
      (slabs_nl_stokes_problem_pscale_press_mass_data_t *)
      ymir_stokes_pc_set_pscale_fn_data;

    ymir_vec_destroy (d->pressure);
    YMIR_FREE (ymir_stokes_pc_set_pscale_fn_data);
  }
  else if (   ymir_stokes_pc_set_pscale_fn
           == slabs_nl_stokes_problem_pscale_Bp ) {
    slabs_nl_stokes_problem_pscale_Bp_data_t  *d =
      (slabs_nl_stokes_problem_pscale_Bp_data_t *)
      ymir_stokes_pc_set_pscale_fn_data;

    ymir_vec_destroy (d->pressure);
    YMIR_FREE (ymir_stokes_pc_set_pscale_fn_data);
  }
}

