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

#include <slabs_setup.h>
#include <rhea.h>
#include <ymir.h>
#include <ymir_comm.h>
#include <ymir_stokes_vec.h>
#include <ymir_perf_counter.h>
#include <ymir_monitor.h>

/* for performance statistics */
#include <ymir_gmg_hierarchy_mesh.h>
#include <ymir_gmg_hierarchy_stress.h>
#include <ymir_gmg_hierarchy_stiff.h>
#include <ymir_gmg_hierarchy_bbt.h>

// for nonlinear solver
#include <slabs_norm.h>
#include <slabs_vtk.h>
#include <ymir_velocity_vec.h>

#if defined(__bgq__)
#include <ymir_bgq.h>
#endif

//###DEV###
#include <ymir_interp_vec.h>
#include <ymir_mass_vec.h>
#include <ymir_pressure_vec.h>
#include <ymir_vtk.h>

/* perfomance counters */
typedef enum
{
  SLABS_PERF_COUNTER_SETUP_MESH,
  SLABS_PERF_COUNTER_SETUP_STOKES,
  SLABS_PERF_COUNTER_SOLVE_STOKES,
  SLABS_PERF_COUNTER_MATVECS,
  SLABS_PERF_COUNTER_TOTAL,
  SLABS_PERF_COUNTER_N
}
slabs_perf_counter_idx_t;
ymir_perf_counter_t slabs_perf_counter[SLABS_PERF_COUNTER_N];
const char         *slabs_perf_counter_name[SLABS_PERF_COUNTER_N] =
{
  "Setup mesh",
  "Setup Stokes",
  "Solve Stokes",
  "All matvecs (setup + solve)",
  "Total"
};
sc_statinfo_t       slabs_perf_stats[
                      SLABS_PERF_COUNTER_N * YMIR_PERF_COUNTER_N_STATS];
char                slabs_perf_stats_name[
                      SLABS_PERF_COUNTER_N * YMIR_PERF_COUNTER_N_STATS][
                      YMIR_PERF_COUNTER_NAME_SIZE];
int                 slabs_perf_n_stats;

/* declare local functions */
static void
slabs_test_stress (MPI_Comm mpicomm,
                   const slabs_domain_shape_t domain_shape,
                   const int minlevel, const int maxlevel,
                   char *refine, const int order);

static void
slabs_test_stiffness (MPI_Comm mpicomm,
                      const slabs_domain_shape_t domain_shape,
                      const int minlevel, const int maxlevel,
                      char *refine, const int order);

static void
slabs_test_mesh_independence_stokes_norm (MPI_Comm mpicomm,
                                          const slabs_domain_shape_t
                                            domain_shape,
                                          const int minlevel,
                                          const int maxlevel,
                                          char *refine,
                                          const int order,
                                          double mass_scaling,
                                          const int max_refine_steps);

static void
slabs_test_interpolation (MPI_Comm mpicomm,
                          const slabs_domain_shape_t domain_shape,
                          const int minlevel, const int maxlevel,
                          char *refine, const int order);

/**
 * Computes vector residual, its norm, and the residual "function."
 */
static double
slabs_nl_solver_compute_residual (ymir_vec_t *residual_up,
                                  ymir_vec_t *residual_lump_up,
                                  ymir_vec_t *residual_dual,
                                  double *norm_vel, double *norm_press,
                                  double *norm_dual,
                                  ymir_vec_t *up,
                                  ymir_dvec_t *dual,
                                  slabs_nl_stokes_problem_t *nl_stokes,
                                  slabs_nl_solver_options_t *solver_options)
{
  const slabs_nl_solver_primaldual_type_t  nl_solver_primaldual_type =
    solver_options->nl_solver_primaldual_type;
  const slabs_krylov_type_t  k_type = solver_options->krylov_type;
  const slabs_norm_type_t    norm_type = solver_options->norm_type;
  ymir_Hminus1_norm_op_t    *norm_op = solver_options->norm_op;
  ymir_stokes_op_t   *stokes_op = nl_stokes->stokes_op;
  ymir_stokes_pc_t   *stokes_pc = nl_stokes->stokes_pc;
  ymir_vec_t         *rhs_u_point = nl_stokes->rhs_u_point;
  double              norm_res;

  /* check input */
  YMIR_ASSERT (k_type == SL_KRYLOV_GMRES);

  /* compute residual for momentum and mass equations */
  norm_res = slabs_norm_compute_residual (residual_up, norm_vel, norm_press,
                                          up, rhs_u_point,
                                          stokes_op, stokes_pc, k_type,
                                          norm_type, norm_op);

  /* compute additional residual for primal-dual method */
  if (nl_solver_primaldual_type != SL_NL_SOLVER_PRIMALDUAL_NONE) {
    /* compute residual of dual equation */
    *norm_dual = slabs_norm_compute_residual_dual (
        residual_dual, up, dual, stokes_op, nl_solver_primaldual_type, k_type);

    /* update combined norm */
    norm_res = sqrt (SC_SQR (norm_res) + SC_SQR (*norm_dual));
  }

  /* weight residual with inverse mass matrix */
  if (residual_lump_up != NULL) {
    ymir_vec_copy (residual_up, residual_lump_up);
    slabs_norm_weight_residual (residual_lump_up, nl_stokes->press_elem);
  }

  /* return norm of residual */
  return norm_res;
}

/**
 * Calculate the forcing for the nonlinear solver, i.e., the Krylov relative
 * tolerance.
 */
static void
slabs_nl_solver_compute_forcing (const int nl_iter,
                                 const double res_new,
                                 const double res_prev,
                                 const double res_reduction_prev,
                                 slabs_nl_solver_options_t *solver_options)
{
  const int           grid_cont_init_steps =
                        solver_options->grid_continuation_init_steps;
  const double        grid_cont_init_forcing =
                        solver_options->grid_continuation_init_forcing;
  double              k_rtol;

  if (0 == nl_iter) {
    /* set Krylov rtol for first nl. iteration */
    k_rtol = solver_options->krylov_rtol_init;

    YMIR_GLOBAL_INFOF ("NL iter %d, set initial Krylov rtol %1.3e\n",
                       nl_iter, k_rtol);
  }
  else if (   0 < grid_cont_init_steps && nl_iter <= (grid_cont_init_steps + 1)
           && 0.0 < grid_cont_init_forcing ) {
    /* set Krylov rtol for grid continuation */
    k_rtol = grid_cont_init_forcing;

    YMIR_GLOBAL_INFOF ("NL iter %d, set grid continuation Krylov rtol "
                       "%1.3e\n", nl_iter, k_rtol);
  }
  else {
    double              forcing_exp = solver_options->nl_forcing_exponent;
    const double        forcing_max = solver_options->nl_forcing_max;
    const double        forcing_max_prog_iter = (double)
                          solver_options->nl_forcing_max_progressive_iter;
    const int           saveguard = solver_options->nl_forcing_saveguard;
    const double        saveguard_threshold =
                          solver_options->nl_forcing_saveguard_threshold;
    const double        grid_cont_init_forcing_exp =
                          solver_options->grid_continuation_init_forcing_exp;
    double              progressive_reduction;
    double              forcing_min;

    /* use different forcing exponent during initial grid continuation */
    if (   0 < grid_cont_init_steps && nl_iter <= (grid_cont_init_steps + 1)
        && 1.0 < grid_cont_init_forcing_exp ) {
      forcing_exp = grid_cont_init_forcing_exp;

      YMIR_GLOBAL_INFOF ("NL iter %d, set grid continuation forcing "
                         "exponent %1.3e\n", nl_iter, forcing_exp);
    }

    /* set Krylov rtol for inexact Newton method according to
     * "Choice 2" in [Eisenstat, Walker, 1996] */
    k_rtol = forcing_max * pow (res_new / res_prev, forcing_exp);

    /* set progressive reduction factor:
     *
     *   ( exp(-(#iter / prog_iter)^2) + eps/forcing ) / ( 1 + eps/forcing )
     */
    if (0.0 < forcing_max_prog_iter) {
      progressive_reduction =
        (  exp (- SC_SQR (((double) nl_iter) / forcing_max_prog_iter))
         + SC_1000_EPS / SC_MIN (k_rtol, forcing_max) )
        / (1.0 + SC_1000_EPS / SC_MIN (k_rtol, forcing_max));
    }
    else {
      progressive_reduction = 1.0;
    }

    /* set min Krylov rtol, a saveguard to avoid oversolving */
    if (saveguard) {
      forcing_min = forcing_max
        * pow (0.5 * (SC_MIN (res_reduction_prev, 1.0) + res_new / res_prev),
               forcing_exp);
    }
    else {
      forcing_min = 0.0;
    }

    YMIR_GLOBAL_INFOF ("NL iter %d, candidate for Krylov rtol %1.3e, "
                       "max %1.3e, progressive reduction %1.3e, "
                       "saveguard %1.3e (threshold %1.3e)\n",
                       nl_iter, k_rtol, forcing_max, progressive_reduction,
                       forcing_min, saveguard_threshold);

    /* deactivate saveguard below a given threshold */
    if (forcing_min <= saveguard_threshold) {
      forcing_min = 0.0;
    }

    /* set Krylov rtol */
    k_rtol = SC_MIN (k_rtol, forcing_max);
    k_rtol *= progressive_reduction;
    k_rtol = SC_MAX (k_rtol, forcing_min);

    YMIR_GLOBAL_INFOF ("NL iter %d, set adaptive Krylov rtol %1.3e\n",
                       nl_iter, k_rtol);
  }

  /* store the Krylov relative tolerance */
  solver_options->krylov_rtol = k_rtol;
}

#ifdef YMIR_PETSC
static PetscErrorCode
slabs_nl_solver_newton_step_stress_residual (KSP ksp, PetscInt n,
                                             PetscReal rnorm, void * data)
{
  PetscErrorCode      ierr;

  ierr = PetscPrintf (PETSC_COMM_WORLD,
                      "%3D KSP Viscous stress residual norm %1.12e\n",
                      n, rnorm);
  YMIR_CHECK_PETSC (ierr);

  return 0;
}
#endif

/**
 * Solves one Newton step for velocity and pressure.
 */
static int
slabs_nl_solver_newton_step (ymir_vec_t *step_up,
                             double *k_res_reduction,
                             ymir_vec_t *rhs_up,
                             slabs_nl_stokes_problem_t *nl_stokes,
                             slabs_nl_solver_options_t *solver_options)
{
  const int           enforce_unscaled_reduction =
                        solver_options->enforce_unscaled_reduction;
  const int           k_maxiter_total = solver_options->krylov_maxiter;
  const double        k_rtol_total = solver_options->krylov_rtol;
  const double        k_atol = solver_options->krylov_atol;
  const int           k_gmres_num_vecs = solver_options->krylov_gmres_num_vecs;
  int                 k_maxiter = k_maxiter_total;
  double              k_rtol = k_rtol_total;
  double              k_res_init, k_res_new;
  int                 k_num_iter = 0;
  int                 k_num_iter_total = 0;
  int                 k = 0;
  ymir_vec_t         *residual_up = ymir_stokes_vec_new (
                          nl_stokes->mesh, nl_stokes->press_elem);

  /* solve viscous stress block only */
#ifdef YMIR_PETSC
  if (solver_options->lin_solve_stress_block_only) {
    ymir_vec_t         *u, *rhs_u;
    KSP                 ksp = nl_stokes->stokes_pc->stress_pc->ksp;
    PetscInt            ksp_gmres_restart;
    PetscErrorCode      ierr;

    /* get GMRES restart from Stokes KSP */
    ierr = KSPGMRESGetRestart (nl_stokes->stokes_pc->ksp, &ksp_gmres_restart);
    YMIR_CHECK_PETSC (ierr);
    YMIR_ASSERT (0 < ksp_gmres_restart);

    /* modify viscous stress KSP parameters */
    ierr = KSPSetFromOptions (ksp); YMIR_CHECK_PETSC (ierr);
    ierr = KSPSetType (ksp, KSPGMRES); YMIR_CHECK_PETSC (ierr);
    ierr = KSPSetPCSide (ksp, PC_RIGHT); YMIR_CHECK_PETSC (ierr);
    ierr = KSPSetNormType (ksp, KSP_NORM_UNPRECONDITIONED);
    YMIR_CHECK_PETSC (ierr);
    ierr = KSPGMRESSetRestart (ksp, ksp_gmres_restart); YMIR_CHECK_PETSC (ierr);
    ierr = KSPSetTolerances (ksp, k_rtol_total, k_atol, 1.e6, k_maxiter);
    YMIR_CHECK_PETSC (ierr);
    ierr = KSPMonitorSet (ksp, slabs_nl_solver_newton_step_stress_residual,
                          NULL, NULL); YMIR_CHECK_PETSC (ierr);

    /* run solver */
    slabs_stokes_vec_get_components_view (&u, NULL, step_up);
    slabs_stokes_vec_get_components_view (&rhs_u, NULL, rhs_up);
    ymir_petsc_ksp_solve (u, rhs_u, ksp, 0 /* zero initial guess */,
                          YMIR_STOKES_PC_VEL_PETSCLAYOUT_TYPE, NULL);
    ymir_vec_destroy (u);
    ymir_vec_destroy (rhs_u);

    /* restore viscous stress KSP parameters */
    ierr = KSPSetPCSide (ksp, PC_SIDE_DEFAULT); YMIR_CHECK_PETSC (ierr);
    ierr = KSPSetNormType (ksp, KSP_NORM_NONE); YMIR_CHECK_PETSC (ierr);
    ierr = KSPMonitorCancel (ksp);
    ierr = KSPSetFromOptions (ksp); YMIR_CHECK_PETSC (ierr);
  }
#endif

  /* compute l^2 norm of unscaled residual before running Krylov solver */
  k_res_init = ymir_vec_norm (rhs_up);

  /* solve algebraic system */
  *k_res_reduction = 1.0;
  while ( 0 == k_num_iter_total ||
          (   enforce_unscaled_reduction
           && 0 < k_num_iter
           && ((double) k_num_iter_total) < (0.9 * (double) k_maxiter_total)
           && (k_rtol_total / *k_res_reduction) < 0.9 ) ) {

    YMIR_GLOBAL_INFOF ("Krylov step %i: l^2 norm of initial Krylov residual "
                       "%1.3e\n", k, k_res_init);

    /* run linear solver (use zero initial guess) */
    k_num_iter = -1;
    ymir_stokes_pc_solve (rhs_up, step_up, nl_stokes->stokes_pc,
                          (0 < k_num_iter_total) /* nonzero initial guess */,
                          k_rtol, k_atol, k_maxiter, k_gmres_num_vecs,
                          &k_num_iter);
    YMIR_ASSERT (-1 < k_num_iter);
    k_num_iter_total += k_num_iter;

    /* compute l^2 norm of actual (i.e., unscaled) residual */
    ymir_nlstokes_op_apply (step_up, residual_up, nl_stokes->stokes_op);
    ymir_vec_add (-1.0, rhs_up, residual_up);
    k_res_new = ymir_vec_norm (residual_up);

    YMIR_GLOBAL_INFOF ("Krylov step %i: l^2 norm of Krylov residual %1.3e\n",
                       k, k_res_new);

    /* update residual reduction */
    *k_res_reduction = k_res_new / k_res_init;

    YMIR_GLOBAL_INFOF ("Krylov step %i: %i Krylov iter for rtol %1.3e "
                       "(prescribed rtol %1.3e)\n",
                       k, k_num_iter_total, *k_res_reduction, k_rtol_total);

    /* update Krylov rtol and maxiter */
    k_rtol = k_rtol * k_rtol_total / *k_res_reduction;
    k_maxiter = k_maxiter_total - k_num_iter_total;
    k++;
  }

  /* destroy */
  ymir_vec_destroy (residual_up);

  /* set residual reduction and return the number of Krylov iterations */
  solver_options->krylov_residual_reduction = *k_res_reduction;
  return k_num_iter_total;
}

/**
 * Computes dual step for primal-dual Newton.
 */
static void
slabs_nl_solver_newton_dual_step (ymir_dvec_t *step_dual,
                                  ymir_dvec_t *dual,
                                  const double step_length,
                                  ymir_vec_t *up,
                                  ymir_vec_t *step_up,
                                  ymir_dvec_t *residual_dual,
                                  slabs_nl_stokes_problem_t *nl_stokes,
                                  slabs_physics_options_t *physics_options,
                                  slabs_nl_solver_options_t *solver_options)
{
  const char         *this_fn_name = "slabs_nl_solver_newton_dual_step";
  const slabs_nl_solver_primaldual_type_t  nl_solver_primaldual_type =
                        solver_options->nl_solver_primaldual_type;
  const slabs_nl_solver_primaldual_scal_type_t  nl_solver_primaldual_scal_type =
                        solver_options->nl_solver_primaldual_scal_type;

  YMIR_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  switch (nl_solver_primaldual_type) {
  case SL_NL_SOLVER_PRIMALDUAL_NONE:
    /* nothing to do */
    break;

  case SL_NL_SOLVER_PRIMALDUAL_NORMSTRAIN:
    {
      const double        visc_IIe_reg =
                            physics_options->viscosity_IIe_regularization;
      mangll_t           *mangll = nl_stokes->mesh->ma;
      const mangll_locidx_t  n_elements = mangll->mesh->K;
      const int           N = ymir_n (mangll->N);
      const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);
      ymir_cvec_t        *vel_vec, *step_vel_vec;
      sc_dmatrix_t       *vel_el_mat, *step_vel_el_mat;
      sc_dmatrix_t       *dual_el_mat, *step_dual_el_mat;
      sc_dmatrix_t       *grad_vel_el_mat, *tmp_dvel, *tmp_3d;
      sc_dmatrix_t       *strain_frob_el_mat;
      sc_dmatrix_t       *strain_vel_el_mat, *strain_step_vel_el_mat;
      sc_dmatrix_t       *tmp_1d, *tmp_6d;
      mangll_locidx_t     elid;

      /* create views */
      slabs_stokes_vec_get_components_view (&vel_vec, NULL, up);
      slabs_stokes_vec_get_components_view (&step_vel_vec, NULL, step_up);

      /* create work variables */
      vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
      step_vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
      dual_el_mat = sc_dmatrix_new (n_nodes_per_el, 6);
      step_dual_el_mat = sc_dmatrix_new (n_nodes_per_el, 6);
      grad_vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 9);
      tmp_dvel = sc_dmatrix_new (n_nodes_per_el, 3);
      tmp_3d = sc_dmatrix_new (n_nodes_per_el, 3);
      strain_frob_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
      strain_vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 6);
      strain_step_vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 6);
      tmp_1d = sc_dmatrix_new (n_nodes_per_el, 1);
      tmp_6d = sc_dmatrix_new (n_nodes_per_el, 6);

      for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
        /* get velocity of this element at GLL nodes */
        ymir_cvec_get_elem_interp (vel_vec, vel_el_mat, YMIR_STRIDE_NODE, elid,
                                   YMIR_GLL_NODE, YMIR_READ);

        /* get velocity step of this element at GLL nodes */
        ymir_cvec_get_elem_interp (step_vel_vec, step_vel_el_mat,
                                   YMIR_STRIDE_NODE, elid, YMIR_GLL_NODE,
                                   YMIR_READ);

        /* get dual tensor for this element and normalize it */
        ymir_dvec_get_elem (dual, dual_el_mat, YMIR_STRIDE_NODE, elid,
                            YMIR_COPY);

        /* get dual residual for this element
         * (store in `step_dual_el_mat`) */
        ymir_dvec_get_elem (residual_dual, step_dual_el_mat, YMIR_STRIDE_NODE,
                            elid, YMIR_COPY);

        /* compute strain rate of the velocity iterate */
        slabs_second_invariant_elem (vel_el_mat, strain_frob_el_mat, mangll,
                                     elid, grad_vel_el_mat, tmp_dvel, tmp_3d,
                                     SL_GAUSS_NODE);
        slabs_strain_rate_tensor_elem (strain_vel_el_mat, mangll,
                                       grad_vel_el_mat);

        /* compute Frobenius-norm of the strain rate tensor;
         * regularize Frobenius-norm to avoid division by zero */
        sc_dmatrix_scale (2.0, strain_frob_el_mat);
        sc_dmatrix_shift (visc_IIe_reg, strain_frob_el_mat);
        sc_dmatrix_sqrt (strain_frob_el_mat, strain_frob_el_mat);

        /* compute strain rate of the velocity step */
        slabs_second_invariant_elem (step_vel_el_mat, tmp_1d, mangll, elid,
                                     grad_vel_el_mat, tmp_dvel, tmp_3d,
                                     SL_GAUSS_NODE);
        slabs_strain_rate_tensor_elem (strain_step_vel_el_mat, mangll,
                                       grad_vel_el_mat);

        /* scale dual tensor */
        switch (nl_solver_primaldual_scal_type) {
        case SL_NL_SOLVER_PRIMALDUAL_SCAL_NONE:
          /* do nothing */
          break;

        case SL_NL_SOLVER_PRIMALDUAL_SCAL_NORMALIZE:
          slabs_norm_symtensor_normalize_frobenius_elem (dual_el_mat);
          break;

        case SL_NL_SOLVER_PRIMALDUAL_SCAL_THRESHOLD:
          slabs_norm_symtensor_normalize_threshold_frobenius_elem (dual_el_mat,
                                                                   1.0);
          break;

        default: /* unknown primal-dual scaling type */
          YMIR_ABORT_NOT_REACHED ();
        }

        /* compute dual step:
         *
         *   1 / ||grad_s (u)||_{F}
         *   * (
         *        residual_dual + grad_s (step_u)
         *      - 1 / (2 * ||grad_s (u)||_{F})
         *      * (  (dual : grad_s (step_u)) * grad_s (u)
         *         + (grad_s (u) : grad_s (step_u)) * dual )
         *     )
         */
        sc_dmatrix_add (1.0, strain_step_vel_el_mat, step_dual_el_mat);

        slabs_norm_symtensor_innerprod_frobenius_elem (
            dual_el_mat, strain_step_vel_el_mat, tmp_1d);
        sc_dmatrix_scale (0.5, tmp_1d);
        sc_dmatrix_dotdivide (strain_frob_el_mat, tmp_1d);
        sc_dmatrix_copy (strain_vel_el_mat, tmp_6d);
        slabs_matrix_multiply_in_1d (tmp_1d, tmp_6d);
        sc_dmatrix_add (-1.0, tmp_6d, step_dual_el_mat);

        slabs_norm_symtensor_innerprod_frobenius_elem (
            strain_vel_el_mat, strain_step_vel_el_mat, tmp_1d);
        sc_dmatrix_scale (0.5, tmp_1d);
        sc_dmatrix_dotdivide (strain_frob_el_mat, tmp_1d);
        sc_dmatrix_copy (dual_el_mat, tmp_6d);
        slabs_matrix_multiply_in_1d (tmp_1d, tmp_6d);
        sc_dmatrix_add (-1.0, tmp_6d, step_dual_el_mat);

        slabs_matrix_divide_in_1d (strain_frob_el_mat, step_dual_el_mat);

        /* set values of dual step vector */
        ymir_dvec_set_elem (step_dual, step_dual_el_mat, YMIR_STRIDE_NODE, elid,
                            YMIR_SET);
      }

      YMIR_ASSERT (sc_dmatrix_is_valid (step_dual->dvec));

      /* update dual iterate */
      ymir_vec_add (step_length, step_dual, dual);

      /* destroy work variables */
      sc_dmatrix_destroy (vel_el_mat);
      sc_dmatrix_destroy (step_vel_el_mat);
      sc_dmatrix_destroy (dual_el_mat);
      sc_dmatrix_destroy (step_dual_el_mat);
      sc_dmatrix_destroy (grad_vel_el_mat);
      sc_dmatrix_destroy (tmp_dvel);
      sc_dmatrix_destroy (tmp_3d);
      sc_dmatrix_destroy (strain_frob_el_mat);
      sc_dmatrix_destroy (strain_vel_el_mat);
      sc_dmatrix_destroy (strain_step_vel_el_mat);
      sc_dmatrix_destroy (tmp_1d);
      sc_dmatrix_destroy (tmp_6d);
      ymir_vec_destroy (vel_vec);
      ymir_vec_destroy (step_vel_vec);
    }
    break;

  case SL_NL_SOLVER_PRIMALDUAL_STRESS:
    {
      ymir_stress_op_t   *stress_op = nl_stokes->stokes_op->stress_op;
      ymir_dvec_t        *viscosity = stress_op->viscosity;
      ymir_dvec_t        *dvisc_dIIe = stress_op->dvdIIe;
      mangll_t           *mangll = nl_stokes->mesh->ma;
      const mangll_locidx_t  n_elements = mangll->mesh->K;
      const int           N = ymir_n (mangll->N);
      const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);
      ymir_cvec_t        *vel_vec, *step_vel_vec;
      sc_dmatrix_t       *vel_el_mat, *step_vel_el_mat;
      sc_dmatrix_t       *viscosity_el_mat, *dvisc_dIIe_el_mat;
      sc_dmatrix_t       *dual_el_mat, *step_dual_el_mat;
      sc_dmatrix_t       *grad_vel_el_mat, *tmp_dvel, *tmp_3d;
      sc_dmatrix_t       *strain_frob_el_mat;
      sc_dmatrix_t       *strain_vel_el_mat, *strain_step_vel_el_mat;
      sc_dmatrix_t       *tmp_1d, *tmp_6d;
      mangll_locidx_t     elid;

      /* create views */
      slabs_stokes_vec_get_components_view (&vel_vec, NULL, up);
      slabs_stokes_vec_get_components_view (&step_vel_vec, NULL, step_up);

      /* create work variables */
      vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
      step_vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
      viscosity_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
      dvisc_dIIe_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
      dual_el_mat = sc_dmatrix_new (n_nodes_per_el, 6);
      step_dual_el_mat = sc_dmatrix_new (n_nodes_per_el, 6);
      grad_vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 9);
      tmp_dvel = sc_dmatrix_new (n_nodes_per_el, 3);
      tmp_3d = sc_dmatrix_new (n_nodes_per_el, 3);
      strain_frob_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
      strain_vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 6);
      strain_step_vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 6);
      tmp_1d = sc_dmatrix_new (n_nodes_per_el, 1);
      tmp_6d = sc_dmatrix_new (n_nodes_per_el, 6);

      for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
        /* get velocity of this element at GLL nodes */
        ymir_cvec_get_elem_interp (vel_vec, vel_el_mat, YMIR_STRIDE_NODE, elid,
                                   YMIR_GLL_NODE, YMIR_READ);

        /* get velocity step of this element at GLL nodes */
        ymir_cvec_get_elem_interp (step_vel_vec, step_vel_el_mat,
                                   YMIR_STRIDE_NODE, elid, YMIR_GLL_NODE,
                                   YMIR_READ);

        /* get viscosity and viscosity derivative for this element */
        ymir_dvec_get_elem (viscosity, viscosity_el_mat, YMIR_STRIDE_NODE, elid,
                            YMIR_READ);
        ymir_dvec_get_elem (dvisc_dIIe, dvisc_dIIe_el_mat, YMIR_STRIDE_NODE,
                            elid, YMIR_READ);

        /* get dual tensor for this element */
        ymir_dvec_get_elem (dual, dual_el_mat, YMIR_STRIDE_NODE, elid,
                            YMIR_COPY);

        /* get dual residual for this element (store in `step_dual_el_mat`) */
        ymir_dvec_get_elem (residual_dual, step_dual_el_mat, YMIR_STRIDE_NODE,
                            elid, YMIR_COPY);

        /* compute strain rate of the velocity iterate */
        slabs_second_invariant_elem (vel_el_mat, tmp_1d, mangll, elid,
                                     grad_vel_el_mat, tmp_dvel, tmp_3d,
                                     SL_GAUSS_NODE);
        slabs_strain_rate_tensor_elem (strain_vel_el_mat, mangll,
                                       grad_vel_el_mat);

        /* compute Frobenius-norm of the strain rate tensor */
        sc_dmatrix_scale (2.0, tmp_1d);
        sc_dmatrix_sqrt (tmp_1d, strain_frob_el_mat);

        /* compute strain rate of the velocity step */
        slabs_second_invariant_elem (step_vel_el_mat, tmp_1d, mangll, elid,
                                     grad_vel_el_mat, tmp_dvel, tmp_3d,
                                     SL_GAUSS_NODE);
        slabs_strain_rate_tensor_elem (strain_step_vel_el_mat, mangll,
                                       grad_vel_el_mat);

        /* scale dual tensor */
        switch (nl_solver_primaldual_scal_type) {
        case SL_NL_SOLVER_PRIMALDUAL_SCAL_NONE:
          /* do nothing */
          break;

        case SL_NL_SOLVER_PRIMALDUAL_SCAL_NORMALIZE:
          sc_dmatrix_copy (strain_frob_el_mat, tmp_1d);
          sc_dmatrix_dotmultiply (viscosity_el_mat, tmp_1d);
          slabs_norm_symtensor_scale_frobenius_elem (dual_el_mat, tmp_1d);
          break;

        case SL_NL_SOLVER_PRIMALDUAL_SCAL_THRESHOLD:
          sc_dmatrix_copy (strain_frob_el_mat, tmp_1d);
          sc_dmatrix_dotmultiply (viscosity_el_mat, tmp_1d);
          slabs_norm_symtensor_scale_threshold_frobenius_elem (dual_el_mat,
                                                               tmp_1d, 1.0);
          break;

        default: /* unknown primal-dual scaling type */
          YMIR_ABORT_NOT_REACHED ();
        }

        /* compute dual step:
         *
         *   2 * viscosity * ( residual_dual + grad_s (step_u) )
         *   +
         *   (2 * dvisc / dIIe) / (2 * viscosity)
         *   * 0.5 * (  (dual : grad_s (step_u)) * grad_s (u)
         *            + (grad_s (u) : grad_s (step_u)) * dual )
         */
        sc_dmatrix_add (1.0, strain_step_vel_el_mat, step_dual_el_mat);
        slabs_matrix_multiply_in_1d (viscosity_el_mat, step_dual_el_mat);

        slabs_norm_symtensor_innerprod_frobenius_elem (
            dual_el_mat, strain_step_vel_el_mat, tmp_1d);
        sc_dmatrix_scale (0.5, tmp_1d);
        sc_dmatrix_dotmultiply (dvisc_dIIe_el_mat, tmp_1d);
        sc_dmatrix_dotdivide (viscosity_el_mat, tmp_1d);
        sc_dmatrix_copy (strain_vel_el_mat, tmp_6d);
        slabs_matrix_multiply_in_1d (tmp_1d, tmp_6d);
        sc_dmatrix_add (1.0, tmp_6d, step_dual_el_mat);

        slabs_norm_symtensor_innerprod_frobenius_elem (
            strain_vel_el_mat, strain_step_vel_el_mat, tmp_1d);
        sc_dmatrix_scale (0.5, tmp_1d);
        sc_dmatrix_dotmultiply (dvisc_dIIe_el_mat, tmp_1d);
        sc_dmatrix_dotdivide (viscosity_el_mat, tmp_1d);
        sc_dmatrix_copy (dual_el_mat, tmp_6d);
        slabs_matrix_multiply_in_1d (tmp_1d, tmp_6d);
        sc_dmatrix_add (1.0, tmp_6d, step_dual_el_mat);

        /* set values of dual step vector */
        ymir_dvec_set_elem (step_dual, step_dual_el_mat, YMIR_STRIDE_NODE, elid,
                            YMIR_SET);
      }

      YMIR_ASSERT (sc_dmatrix_is_valid (step_dual->dvec));

      /* update dual iterate */
      ymir_vec_add (step_length, step_dual, dual);

      /* destroy work variables */
      sc_dmatrix_destroy (vel_el_mat);
      sc_dmatrix_destroy (step_vel_el_mat);
      sc_dmatrix_destroy (viscosity_el_mat);
      sc_dmatrix_destroy (dvisc_dIIe_el_mat);
      sc_dmatrix_destroy (dual_el_mat);
      sc_dmatrix_destroy (step_dual_el_mat);
      sc_dmatrix_destroy (grad_vel_el_mat);
      sc_dmatrix_destroy (tmp_dvel);
      sc_dmatrix_destroy (tmp_3d);
      sc_dmatrix_destroy (strain_frob_el_mat);
      sc_dmatrix_destroy (strain_vel_el_mat);
      sc_dmatrix_destroy (strain_step_vel_el_mat);
      sc_dmatrix_destroy (tmp_1d);
      sc_dmatrix_destroy (tmp_6d);
      ymir_vec_destroy (vel_vec);
      ymir_vec_destroy (step_vel_vec);
    }
    break;

  default: /* unknown primal-dual type */
    YMIR_ABORT_NOT_REACHED ();
  }
}

/**
 * Solves one Picard step for velocity and pressure.
 */
static int
slabs_nl_solver_picard_step (ymir_vec_t *step_up,
                             double *k_res_reduction,
                             ymir_vec_t *rhs_up,
                             slabs_nl_stokes_problem_t *nl_stokes,
                             slabs_nl_solver_options_t *solver_options)
{
  /* solve for velocity and pressure step */
  return slabs_nl_solver_newton_step (step_up, k_res_reduction, rhs_up,
                                      nl_stokes, solver_options);
}

/**
 * Computes dual step for primal-dual Picard.
 */
static void
slabs_nl_solver_picard_dual_step (ymir_dvec_t *step_dual,
                                  ymir_dvec_t *dual,
                                  const double step_length,
                                  ymir_vec_t *up,
                                  ymir_vec_t *step_up,
                                  ymir_dvec_t *residual_dual,
                                  slabs_nl_stokes_problem_t *nl_stokes,
                                  slabs_physics_options_t *physics_options,
                                  slabs_nl_solver_options_t *solver_options)
{
  const char         *this_fn_name = "slabs_nl_solver_picard_dual_step";
  const slabs_nl_solver_primaldual_type_t  nl_solver_primaldual_type =
                        solver_options->nl_solver_primaldual_type;

  YMIR_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  switch (nl_solver_primaldual_type) {
  case SL_NL_SOLVER_PRIMALDUAL_NONE:
    /* nothing to do */
    break;

  case SL_NL_SOLVER_PRIMALDUAL_NORMSTRAIN:
    /* compute dual step:
     *
     *   1 / ||grad_s (u)||_{F} * ( residual_dual + grad_s (step_u) )
     */
    {
      const double        visc_IIe_reg =
                            physics_options->viscosity_IIe_regularization;
      ymir_mesh_t        *mesh = nl_stokes->mesh;
      ymir_cvec_t        *u, *step_u;
      ymir_dvec_t        *strain_frob = ymir_dvec_new (mesh, 1,
                                                       YMIR_GAUSS_NODE);
      ymir_velocity_elem_t  *vel_elem;

      /* compute Frobenius-norm of the strain rate tensor */
      slabs_stokes_vec_get_components_view (&u, NULL, up);
      vel_elem = ymir_velocity_elem_new (mesh->ma->N, mesh->ma->ompsize);
      ymir_second_invariant_vec (u, strain_frob, vel_elem);
      ymir_dvec_scale (2.0, strain_frob);
      ymir_dvec_shift (visc_IIe_reg, strain_frob);
      ymir_dvec_sqrt (strain_frob, strain_frob);

      /* compute strain rate of the velocity step */
      slabs_stokes_vec_get_components_view (&step_u, NULL, step_up);
      ymir_velocity_strain_rate (step_u, step_dual, 0);

      /* add dual residual */
      ymir_dvec_add (1.0, residual_dual, step_dual);

      /* scale by inverse Frobenius-norm of the strain rate tensor */
      ymir_dvec_divide_in1 (strain_frob, step_dual);

      /* update dual iterate */
      ymir_vec_add (step_length, step_dual, dual);

      /* destroy */
      ymir_vec_destroy (u);
      ymir_vec_destroy (step_u);
      ymir_vec_destroy (strain_frob);
      ymir_velocity_elem_destroy (vel_elem);
    }
    break;

  case SL_NL_SOLVER_PRIMALDUAL_STRESS:
    /* compute dual step:
     *
     *   2 * viscosity * ( residual_dual + grad_s (step_u) )
     */
    {
      ymir_stress_op_t   *stress_op = nl_stokes->stokes_op->stress_op;
      ymir_dvec_t        *viscosity = stress_op->viscosity;
      ymir_cvec_t        *step_u;

      /* compute strain rate of the velocity step */
      slabs_stokes_vec_get_components_view (&step_u, NULL, step_up);
      ymir_velocity_strain_rate (step_u, step_dual, 0);

      /* add dual residual */
      ymir_dvec_add (1.0, residual_dual, step_dual);

      /* scale by viscosity */
      ymir_dvec_multiply_in1 (viscosity, step_dual);

      /* update dual iterate */
      ymir_vec_add (step_length, step_dual, dual);

      /* destroy */
      ymir_vec_destroy (step_u);
    }
    break;

  default: /* unknown primal-dual type */
    YMIR_ABORT_NOT_REACHED ();
  }
}

/**
 * Computes slope of linear model in step direction:
 *
 *     d/d alpha [ phi (x + alpha * step) ]
 *   = d/d alpha [ ||res (x + alpha * step)|| ]
 *   = ( res' (x + alpha * step) * step , res (x + alpha * step) )
 *     / ||res (x + alpha * step)||
 *
 * where
 *
 *   x is previous iterate
 *   alpha is step length
 *   res (x) = OP (x) - RHS, for some operator OP and right-hand side RHS
 *   phi (x) = ||res (x)||, for some norm induced by inner product (. , .)
 */
static double
slabs_nl_solver_compute_lin_model_slope (ymir_stokes_op_t *stokes_op,
                                         ymir_vec_t *step_up,
                                         ymir_vec_t *residual_up,
                                         const double res_prev,
                                         slabs_nl_solver_options_t
                                           *solver_options)
{
  const slabs_norm_type_t  norm_type = solver_options->norm_type;
  ymir_Hminus1_norm_op_t  *norm_op = solver_options->norm_op;

  ymir_pressure_elem_t  *press_elem = stokes_op->press_elem;
  ymir_vec_t         *dir_derivative = ymir_vec_template (step_up);
  double              innerprod_vel, innerprod_press;

  /* compute directional derivative in step direction */
  ymir_nlstokes_op_apply (step_up, dir_derivative, stokes_op);

  /* compute inner product with residual */
  slabs_norm_innerprod_comp (&innerprod_vel, &innerprod_press,
                             dir_derivative, residual_up,
                             press_elem, norm_type, norm_op);

  /* destoy */
  ymir_vec_destroy (dir_derivative);

  /* return slope of linear model in step direction */
  return - (SC_SQR (innerprod_vel) + SC_SQR (innerprod_press)) / res_prev;
}

/**
 * Sets adaptive step length  element-wise, for velocity & pressure
 * individually.
 */
static void
slabs_nl_solver_steplength_set_adaptive_length (
                                            sc_dmatrix_t *adaptive_step_length,
                                            ymir_vec_t *rel_residual_up,
                                            ymir_pressure_elem_t *press_elem,
                                            const double step_ln_min,
                                            const double step_ln_max)
{
  const char         *this_fn_name =
                        "slabs_nl_solver_steplength_set_adaptive_length";

  ymir_mesh_t        *mesh = rel_residual_up->mesh;
  const ymir_locidx_t n_elements = mesh->cnodes->K;
  const int           Np = ymir_np (mesh->ma->N);
  ymir_locidx_t       elid;
  int                 nodeid;

  MPI_Comm            mpicomm = mesh->ma->mpicomm;
  int                 mpiret;

  sc_dmatrix_t       *stepln_u_mat = sc_dmatrix_new_view_offset (
                          0, 1, n_elements, adaptive_step_length);
  sc_dmatrix_t       *stepln_p_mat = sc_dmatrix_new_view_offset (
                          1, 1, n_elements, adaptive_step_length);
  double             *_sc_restrict stepln_u_data = stepln_u_mat->e[0];
  double             *_sc_restrict stepln_p_data = stepln_p_mat->e[0];
  double              scal_max_glo[2], scal_max_loc[2] = {0.0, 0.0};
  double              s_vel, s_press;

  ymir_vec_t         *rel_res_u, *rel_res_p;
  sc_dmatrix_t       *elemvel = sc_dmatrix_new (Np, 3);
  sc_dmatrix_t       *elempress = sc_dmatrix_new (1, Np);
  double             *_sc_restrict elemvel_data = elemvel->e[0];
  double             *_sc_restrict elempress_data = elempress->e[0];

  /* calculate mean values of relative residuals in each element */
  slabs_stokes_vec_get_components (&rel_res_u, &rel_res_p,
                                   rel_residual_up, press_elem);
  for (elid = 0; elid < n_elements; elid++) {
    ymir_vec_get_elem_interp (rel_res_u, elemvel, YMIR_STRIDE_NODE,
                              elid, YMIR_GLL_NODE, YMIR_COPY);
    ymir_vec_get_elem_interp (rel_res_p, elempress, YMIR_STRIDE_COMP,
                              elid, YMIR_GAUSS_NODE, YMIR_COPY);

    s_vel = 0.0;
    s_press = 0.0;
    for (nodeid = 0; nodeid < Np; nodeid++) {
      s_vel += sqrt ( elemvel_data[3*nodeid    ] * elemvel_data[3*nodeid    ] +
                      elemvel_data[3*nodeid + 1] * elemvel_data[3*nodeid + 1] +
                      elemvel_data[3*nodeid + 2] * elemvel_data[3*nodeid + 2] );
      s_press += elempress_data[nodeid];
    }
    stepln_u_data[elid] = s_vel / (double) Np;
    stepln_p_data[elid] = s_press / (double) Np;

    scal_max_loc[0] = SC_MAX (scal_max_loc[0], stepln_u_data[elid]);
    scal_max_loc[1] = SC_MAX (scal_max_loc[1], stepln_p_data[elid]);
  }

  /* communicate to get global max */
  mpiret = MPI_Allreduce (scal_max_loc, scal_max_glo, 2, MPI_DOUBLE, MPI_MAX,
                          mpicomm); YMIR_CHECK_MPI (mpiret);
  YMIR_GLOBAL_INFOF ("%s: Max relative residuals: vel %.3e, press %.3e\n",
                     this_fn_name, scal_max_glo[0], scal_max_glo[1]);
  YMIR_ASSERT (0.0 < scal_max_glo[0]);
  YMIR_ASSERT (0.0 < scal_max_glo[1]);

  /* scale relative residuals to get adaptive step length */
  sc_dmatrix_scale_shift ((step_ln_min - step_ln_max) / scal_max_glo[0],
                          step_ln_max, stepln_u_mat);
  sc_dmatrix_scale_shift ((step_ln_min - step_ln_max) / scal_max_glo[1],
                          step_ln_max, stepln_p_mat);

  /* destroy */
  sc_dmatrix_destroy (stepln_u_mat);
  sc_dmatrix_destroy (stepln_p_mat);
  sc_dmatrix_destroy (elemvel);
  sc_dmatrix_destroy (elempress);
  ymir_vec_destroy (rel_res_u);
  ymir_vec_destroy (rel_res_p);
}

/**
 *
 */
static void
slabs_nl_solver_steplength_apply_adaptive_length (
                                           ymir_vec_t *step_up,
                                           sc_dmatrix_t *adaptive_step_length,
                                           ymir_pressure_elem_t *press_elem)
{
  ymir_mesh_t        *mesh = step_up->mesh;
  const ymir_locidx_t n_elements = mesh->cnodes->K;
  const int           Np = ymir_np (mesh->ma->N);
  ymir_locidx_t       elid;

  sc_dmatrix_t       *stepln_u_mat = sc_dmatrix_new_view_offset (
                          0, 1, n_elements, adaptive_step_length);
  sc_dmatrix_t       *stepln_p_mat = sc_dmatrix_new_view_offset (
                          1, 1, n_elements, adaptive_step_length);
  double             *_sc_restrict stepln_u_data = stepln_u_mat->e[0];
  double             *_sc_restrict stepln_p_data = stepln_p_mat->e[0];

  ymir_vec_t         *step_u, *step_p;
  sc_dmatrix_t       *elemvel, *elempress;

  /* get velocity and pressure components of step */
  slabs_stokes_vec_get_components (&step_u, &step_p, step_up, press_elem);

  /* reduce step length of velocity step */
  elemvel = sc_dmatrix_new (Np, 3);
  for (elid = 0; elid < n_elements; elid++) {
    ymir_vec_get_elem_interp (step_u, elemvel, YMIR_STRIDE_NODE,
                              elid, YMIR_GLL_NODE, YMIR_COPY);
    sc_dmatrix_scale (stepln_u_data[elid], elemvel);
    ymir_vec_set_elem_interp (step_u, elemvel, YMIR_STRIDE_NODE,
                              elid, YMIR_GLL_NODE, YMIR_SET);
  }
  ymir_vec_share_owned (step_u);

  /* reduce step length of pressure step */
  if (press_elem->space == YMIR_PRESSURE_SPACE_STAB) {
    elempress = sc_dmatrix_new (1, Np);
    for (elid = 0; elid < n_elements; elid++) {
      ymir_vec_get_elem_interp (step_p, elempress, YMIR_STRIDE_COMP,
                                elid, YMIR_GLL_NODE, YMIR_COPY);
      sc_dmatrix_scale (stepln_p_data[elid], elempress);
      ymir_vec_set_elem_interp (step_p, elempress, YMIR_STRIDE_COMP,
                                elid, YMIR_GLL_NODE, YMIR_SET);
    }
    ymir_vec_share_owned (step_p);
  }
  else {
    elempress = sc_dmatrix_new (0, 0);
    for (elid = 0; elid < n_elements; elid++) {
      ymir_evec_get_elem (step_p, elempress, elid, YMIR_RW);
      sc_dmatrix_scale (stepln_p_data[elid], elempress);
      ymir_evec_set_elem (step_p, elempress, elid, YMIR_SET);
    }
  }

  /* set velocity and pressure components of step */
  slabs_stokes_vec_set_components (step_up, step_u, step_p, press_elem);

  /* destroy */
  sc_dmatrix_destroy (stepln_u_mat);
  sc_dmatrix_destroy (stepln_p_mat);
  sc_dmatrix_destroy (elemvel);
  sc_dmatrix_destroy (elempress);
  ymir_vec_destroy (step_u);
  ymir_vec_destroy (step_p);
}

/**
 * Moves iterate into the step direction with given step length.
 */
static void
slabs_nl_solver_steplength_move (ymir_vec_t *up,
                                 ymir_vec_t *dual,
                                 ymir_vec_t *step_up,
                                 ymir_vec_t *step_dual,
                                 ymir_vec_t *residual_dual,
                                 ymir_vec_t *up_prev,
                                 const double step_length,
                                 sc_dmatrix_t *adaptive_step_length,
                                 const int step_reduction_iter,
                                 slabs_nl_stokes_problem_t *nl_stokes,
                                 slabs_physics_options_t *physics_options,
                                 slabs_nl_solver_options_t *solver_options)
{
  const slabs_nl_solver_primaldual_type_t  primaldual_type =
                        solver_options->nl_solver_primaldual_type;
  const slabs_nl_solver_step_reduction_type_t  step_reduction_type =
                        solver_options->nl_step_length_reduction_type;

  /* check input */
  YMIR_ASSERT (up != NULL);
  YMIR_ASSERT (step_up != NULL);

  /* move (velocity,pressure) iterate */
  switch (step_reduction_type) {
  case SL_NL_SOLVER_STEP_REDUCTION_CONST:
  case SL_NL_SOLVER_STEP_REDUCTION_QUAD_INTERP:
    YMIR_ASSERT (0.0 < step_length);
    ymir_vec_add (step_length, step_up, up);
    break;
  case SL_NL_SOLVER_STEP_REDUCTION_ADAPTIVE:
    if (0 == step_reduction_iter) {
      ymir_vec_add (1.0, step_up, up);
    }
    else {
      ymir_vec_t         *reduced_step_up = ymir_vec_clone (step_up);

      YMIR_ASSERT (adaptive_step_length != NULL);

      slabs_nl_solver_steplength_apply_adaptive_length (
          reduced_step_up, adaptive_step_length, nl_stokes->press_elem);
      ymir_vec_add (1.0, reduced_step_up, up);

      ymir_vec_destroy (reduced_step_up);
    }
    break;
  default: /* unknown step length reduction type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* compute dual step and move dual iterate */
  if (primaldual_type == SL_NL_SOLVER_PRIMALDUAL_STRESS) {
    YMIR_ASSERT (dual != NULL);
    YMIR_ASSERT (step_dual != NULL);
    YMIR_ASSERT (residual_dual != NULL);
    YMIR_ASSERT (up_prev != NULL);

    switch (step_reduction_type) {
    case SL_NL_SOLVER_STEP_REDUCTION_CONST:
    case SL_NL_SOLVER_STEP_REDUCTION_QUAD_INTERP:
      slabs_nl_solver_newton_dual_step (step_dual, dual,
                                        step_length, up_prev, step_up,
                                        residual_dual, nl_stokes,
                                        physics_options, solver_options);
      break;
    default: /* unknown step length reduction type */
      YMIR_ABORT_NOT_REACHED ();
    }
  }
}

/**
 * Updates the nonlinear Stokes operator.
 */
static void
slabs_nl_solver_steplength_update_operator (
                                    slabs_stokes_state_t *state,
                                    ymir_vec_t *dual,
                                    ymir_vec_t *step_up,
                                    ymir_vec_t *step_dual,
                                    ymir_vec_t *residual_dual,
                                    const double step_length,
                                    slabs_nl_stokes_problem_t *nl_stokes,
                                    slabs_physics_options_t *physics_options,
                                    slabs_nl_solver_options_t *solver_options)
{
  const slabs_nl_solver_type_t  nl_solver_type = solver_options->nl_solver_type;
  const slabs_nl_solver_primaldual_type_t  primaldual_type =
                        solver_options->nl_solver_primaldual_type;
  const slabs_nl_solver_step_reduction_type_t  step_reduction_type =
                        solver_options->nl_step_length_reduction_type;

  /* update Stokes operator */
  slabs_nonlinear_stokes_op_update (nl_stokes, state, physics_options,
                                    nl_solver_type);

  /* compute dual step and update dual iterate */
  if (primaldual_type == SL_NL_SOLVER_PRIMALDUAL_NORMSTRAIN) {
    ymir_vec_t         *up = state->vel_press_vec;

    YMIR_ASSERT (dual != NULL);
    YMIR_ASSERT (step_dual != NULL);
    YMIR_ASSERT (residual_dual != NULL);

    switch (step_reduction_type) {
    case SL_NL_SOLVER_STEP_REDUCTION_CONST:
    case SL_NL_SOLVER_STEP_REDUCTION_QUAD_INTERP:
      slabs_nl_solver_newton_dual_step (step_dual, dual,
                                        step_length, up, step_up,
                                        residual_dual, nl_stokes,
                                        physics_options, solver_options);
      break;
    default: /* unknown step length reduction type */
      YMIR_ABORT_NOT_REACHED ();
    }
  }

}

/**
 * Computes residual norm at new iterate position.
 */
static double
slabs_nl_solver_steplength_compute_residual (
                                     double *res_norm_vel,
                                     double *res_norm_press,
                                     double *res_norm_dual,
                                     ymir_vec_t *residual_up,
                                     ymir_vec_t *residual_dual,
                                     ymir_vec_t *up,
                                     ymir_vec_t *dual,
                                     slabs_nl_stokes_problem_t *nl_stokes,
                                     slabs_physics_options_t *physics_options,
                                     slabs_nl_solver_options_t *solver_options)
{
  const slabs_nl_solver_primaldual_type_t  primaldual_type =
                        solver_options->nl_solver_primaldual_type;
  const slabs_krylov_type_t  k_type = solver_options->krylov_type;
  ymir_stokes_op_t   *stokes_op = nl_stokes->stokes_op;
  ymir_stokes_pc_t   *stokes_pc = nl_stokes->stokes_pc;
  ymir_cvec_t        *rhs_u_point = nl_stokes->rhs_u_point;
  const slabs_norm_type_t  norm_type = solver_options->norm_type;
  ymir_Hminus1_norm_op_t  *norm_op = solver_options->norm_op;
  double              res_norm;

  /* compute residual norm at new position */
  res_norm = slabs_norm_compute_residual (residual_up,
                                          res_norm_vel, res_norm_press,
                                          up, rhs_u_point,
                                          stokes_op, stokes_pc, k_type,
                                          norm_type, norm_op);

  /* update residual norm in case of primal-dual method */
  if (primaldual_type != SL_NL_SOLVER_PRIMALDUAL_NONE) {
    /* compute residual of dual equation */
    *res_norm_dual = slabs_norm_compute_residual_dual (
        residual_dual, up, dual, stokes_op, primaldual_type, k_type);

    /* update combined norm */
    res_norm = sqrt (SC_SQR (res_norm) + SC_SQR (*res_norm_dual));
  }

  /* return norm of residual */
  return res_norm;
}

/**
 * Reduces the step length.
 */
static int
slabs_nl_solver_steplength_reduce (double *step_length,
                                   sc_dmatrix_t *adaptive_step_length,
                                   double *search_lin_model_slope,
                                   ymir_vec_t *residual_up,
                                   ymir_vec_t *residual_up_prev,
                                   const double res_norm_new,
                                   const double res_norm_prev,
                                   const double res_norm_prev_vel,
                                   const double res_norm_prev_press,
                                   slabs_stokes_state_t *state,
                                   ymir_vec_t *step_up,
                                   ymir_vec_t *up_prev,
                                   const int step_reduction_iter,
                                   slabs_nl_stokes_problem_t *nl_stokes,
                                   slabs_physics_options_t *physics_options,
                                   slabs_nl_solver_options_t *solver_options)
{
  const char         *this_fn_name = "slabs_nl_solver_steplength_reduce";
  const slabs_nl_solver_type_t  nl_solver_type = solver_options->nl_solver_type;
  const slabs_nl_solver_step_reduction_type_t  step_reduction_type =
                        solver_options->nl_step_length_reduction_type;
  const double        step_reduction_min =
                        solver_options->nl_step_length_reduction_min;
  const double        step_reduction_max =
                        solver_options->nl_step_length_reduction_max;
  const double        step_length_min = solver_options->nl_step_length_min;
  const double        step_length_prev = *step_length;

  /* check input */
  YMIR_ASSERT (step_length != NULL);

  /* set new step length */
  switch (step_reduction_type) {
  case SL_NL_SOLVER_STEP_REDUCTION_CONST:
    /* reduce step length by a constant factor */
    YMIR_ASSERT (0.0 < *step_length && *step_length <= 1.0);
    *step_length *= 0.5;
    break;

  case SL_NL_SOLVER_STEP_REDUCTION_ADAPTIVE:
    /* compute adaptive step length based on residual reduction:
     *
     *   d_i  = |e_i^T * residual_up| /
     *          (|e_i^T * residual_up_prev| + eps * ||residual_up_prev||)
     *
     *   D_ii = (alpha_max - alpha_min) * (1 - d_i / max(d_j)) + alpha_min
     *        = (alpha_min - alpha_max) / max(d_j) * d_i + alpha_max
     *
     * then
     *
     *   adaptive_step_length = D = diag (D_ii)
     */
    {
      const double        step_reduction_reg =
                            solver_options->nl_step_length_reduction_reg;
      const double        iter_reduc = pow (0.5, (double) step_reduction_iter);
      const double        step_ln_min = iter_reduc * step_reduction_min;
      const double        step_ln_max = iter_reduc * step_reduction_max;
      ymir_pressure_elem_t  *press_elem = nl_stokes->press_elem;
      ymir_vec_t         *rel_residual_up = ymir_vec_template (residual_up);
      ymir_vec_t         *rel_res_u, *rel_res_p;

      YMIR_ASSERT (0.0 < step_reduction_min && step_reduction_min < 1.0);
      YMIR_ASSERT (step_reduction_min <= step_reduction_max &&
                   step_reduction_max < 1.0);
      YMIR_ASSERT (0.0 < res_norm_prev_vel);
      YMIR_ASSERT (0.0 < res_norm_prev_press);
      YMIR_ASSERT (residual_up_prev != NULL);
      YMIR_ASSERT (adaptive_step_length != NULL);

      /* get absolute value of previous residual */
      ymir_vec_copy (residual_up_prev, rel_residual_up);
      ymir_vec_fabs (rel_residual_up, rel_residual_up);

      /* add regularization (separately for velocity and pressure) */
      slabs_stokes_vec_get_components (&rel_res_u, &rel_res_p,
                                       rel_residual_up, press_elem);
      ymir_vec_shift (step_reduction_reg * res_norm_prev_vel, rel_res_u);
      ymir_vec_shift (step_reduction_reg * res_norm_prev_press, rel_res_p);
      slabs_stokes_vec_set_components (rel_residual_up, rel_res_u,
                                       rel_res_p, press_elem);

      /* invert regularized previous residual */
      ymir_vec_reciprocal (rel_residual_up);
      YMIR_ASSERT (sc_dmatrix_is_valid (rel_residual_up->dataown));
      YMIR_ASSERT (sc_dmatrix_is_valid (rel_residual_up->coff));

      /* multiply element-wise with the new residual */
      ymir_vec_multiply_in (residual_up, rel_residual_up);
      ymir_vec_fabs (rel_residual_up, rel_residual_up);

      /* transform to adaptive step length */
      slabs_nl_solver_steplength_set_adaptive_length (
          adaptive_step_length, rel_residual_up, press_elem,
          step_ln_min, step_ln_max);
      YMIR_ASSERT (sc_dmatrix_is_valid (adaptive_step_length));

      /* print statistics */
#ifdef YMIR_DEBUG
      {
        ymir_mesh_t        *mesh = nl_stokes->mesh;
      //ymir_vec_t         *mean = ymir_vec_template (rel_res_u);
      //ymir_vec_t         *unit = ymir_vec_template (rel_res_u);
      //const double        vol = ymir_mass_compute_volume (mesh);

        ymir_vec_t         *res_u, *res_p;
        ymir_vec_t         *res_prev_u, *res_prev_p;
        ymir_vec_t         *adap_step_ln_u, *adap_step_ln_p;
        sc_dmatrix_t       *stepln_u_mat = sc_dmatrix_new_view_offset (
                                0, 1, mesh->cnodes->K, adaptive_step_length);
        sc_dmatrix_t       *stepln_p_mat = sc_dmatrix_new_view_offset (
                                1, 1, mesh->cnodes->K,
                                adaptive_step_length);

        adap_step_ln_u = slabs_dvec_new_from_element_data (stepln_u_mat, mesh,
                                                           YMIR_GLL_NODE);
        adap_step_ln_p = slabs_dvec_new_from_element_data (stepln_p_mat, mesh,
                                                           YMIR_GLL_NODE);
        slabs_stokes_vec_get_components (&res_u, &res_p, residual_up,
                                         press_elem);
        slabs_stokes_vec_get_components (&res_prev_u, &res_prev_p,
                                         residual_up_prev, press_elem);

      //ymir_vec_set_value (unit, 1.0);
      //ymir_mass_apply (unit, mean);
      //ymir_vec_scale (1.0 / vol, mean);
      //*step_length = ymir_vec_innerprod (adap_step_ln_u, mean) / 3.0;

        YMIR_GLOBAL_INFOF (
            "%s: adaptive step length (vel):   min %.3e max %.3e mean %.3e\n",
            this_fn_name,
            ymir_vec_min_global (adap_step_ln_u),
            ymir_vec_max_global (adap_step_ln_u),
            -1.0);

        YMIR_GLOBAL_INFOF (
            "%s: adaptive step length (press): min %.3e max %.3e mean %.3e\n",
            this_fn_name,
            ymir_vec_min_global (adap_step_ln_p),
            ymir_vec_max_global (adap_step_ln_p),
            -1.0);

        ymir_vtk_write (mesh, this_fn_name,
                        adap_step_ln_u, "adap_step_ln_u",
                        adap_step_ln_p, "adap_step_ln_p",
                        res_u, "res_u",
                        res_p, "res_p",
                        res_prev_u, "res_prev_u",
                        res_prev_p, "res_prev_p",
                        NULL);

      //ymir_vec_destroy (mean);
      //ymir_vec_destroy (unit);
        ymir_vec_destroy (adap_step_ln_u);
        ymir_vec_destroy (adap_step_ln_p);
        ymir_vec_destroy (res_u);
        ymir_vec_destroy (res_p);
        ymir_vec_destroy (res_prev_u);
        ymir_vec_destroy (res_prev_p);
      }
#endif

      /* destroy */
      ymir_vec_destroy (rel_res_u);
      ymir_vec_destroy (rel_res_p);
    }
    break;

  case SL_NL_SOLVER_STEP_REDUCTION_QUAD_INTERP:
    /* compute slope of linear model in step direction at step length zero,
     * i.e., alpha = 0
     *
     *     d/d alpha [ phi (x + alpha * step) ] (alpha = 0)
     *   = d/d alpha [ ||res (x + alpha * step)|| ] (alpha = 0)
     *   = ( res' (x) * step , res (x) ) / ||res (x)||
     *
     * where
     *
     *   x is previous iterate
     *   alpha is step length
     *   res (x) = OP (x) - RHS, for some operator OP and right-hand side RHS
     *   phi (x) = ||res (x)||
     *   || . || is some norm induced by inner product ( . , . )
     */
    {
      ymir_vec_t         *up = state->vel_press_vec;
      ymir_stokes_op_t   *stokes_op = nl_stokes->stokes_op;
      double              step_reduction;

      YMIR_ASSERT (0.0 < step_reduction_min && step_reduction_min < 1.0);
      YMIR_ASSERT (step_reduction_min <= step_reduction_max &&
                   step_reduction_max < 1.0);
      YMIR_ASSERT (res_norm_new < res_norm_prev);
      YMIR_ASSERT (search_lin_model_slope != NULL);
      YMIR_ASSERT (residual_up_prev != NULL);
      YMIR_ASSERT (step_up != NULL);
      YMIR_ASSERT (up_prev != NULL);

      if (*search_lin_model_slope < DBL_MIN) { /* if initial step search */
        /* reverse update of Stokes operator */
        ymir_vec_copy (up_prev, up);
        slabs_nonlinear_stokes_op_update (nl_stokes, state, physics_options,
                                          nl_solver_type);

        /* compute slope of linear model in step direction */
        *search_lin_model_slope = slabs_nl_solver_compute_lin_model_slope (
            stokes_op, step_up, residual_up_prev, res_norm_prev,
            solver_options);

        YMIR_GLOBAL_INFOF ("%s: linear model slope for step search %.3e\n",
                           this_fn_name, *search_lin_model_slope);
      }

      /* compute step reduction via quadratic interpolation s.t.
       *
       *   alpha_new = - ( phi' (0) * alpha^2 )
       *               / ( 2.0 * ( phi (alpha) - phi (0) - alpha * phi' (0) ) )
       */
      step_reduction = - *search_lin_model_slope * *step_length / 2.0
                       / ( res_norm_new - res_norm_prev -
                           *step_length * *search_lin_model_slope );

      YMIR_GLOBAL_INFOF ("%s: Candidate for step reduction %.3e, "
                         "min %.3e, max %.3e\n", this_fn_name, step_reduction,
                         step_reduction_min, step_reduction_max);

      /* bound step reduction */
      step_reduction = SC_MIN (SC_MAX (step_reduction, step_reduction_min),
                               step_reduction_max);

      /* set new step length */
      *step_length *= step_reduction;
    }
    break;

  default: /* unknown step length reduction type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* return if the step length is below the prescribed minimum */
  if (0.0 < step_length_min && *step_length < step_length_min) {
    YMIR_GLOBAL_INFOF ("%s: Minimum step length reached (%.3e)\n",
                       this_fn_name, step_length_min);
    *step_length = step_length_prev;
    return 0;
  }
  else {
    return 1;
  }
}

/**
 * Performs line search to reduce the energy associated with the solution.
 */
static double
slabs_nl_solver_steplength_search (slabs_stokes_state_t *state,
                                   ymir_vec_t *step_up,
                                   ymir_vec_t *residual_up,
                                   ymir_dvec_t *step_dual,
                                   ymir_dvec_t *residual_dual,
                                   double *res_new,
                                   double *res_new_vel,
                                   double *res_new_press,
                                   double *res_new_dual,
                                   double *lin_model_slope,
                                   const int nl_iter,
                                   const double res_prev,
                                   const double res_prev_vel,
                                   const double res_prev_press,
                                   slabs_nl_stokes_problem_t *nl_stokes,
                                   slabs_physics_options_t *physics_options,
                                   slabs_nl_solver_options_t *solver_options)
{
  const slabs_nl_solver_type_t  nl_solver_type = solver_options->nl_solver_type;
  const slabs_nl_solver_primaldual_type_t  nl_solver_primaldual_type =
                        solver_options->nl_solver_primaldual_type;
  const slabs_nl_solver_step_reduction_type_t  step_reduction_type =
                        solver_options->nl_step_length_reduction_type;
  const double        k_res_reduction =
                        solver_options->krylov_residual_reduction;
  const double        descend_cond_relax =
                        solver_options->nl_step_length_descend_cond_relax;
  const int           step_search_maxiter = 12;

  ymir_vec_t         *up = state->vel_press_vec;
  ymir_vec_t         *dual = state->dual_vec;
  ymir_vec_t         *up_prev;
  ymir_vec_t         *dual_prev;
  ymir_vec_t         *residual_up_prev;
  sc_dmatrix_t       *adaptive_step_length;
  double              search_lin_model_slope = 0.0;
  double              step_length = 1.0;
  int                 step_reduction_success;
  int                 k;
#ifdef YMIR_DEBUG
  ymir_stokes_op_t   *stokes_op = nl_stokes->stokes_op;
  double              up_norm, step_up_norm, up_normratio;
  double              u_norm, step_u_norm, u_normratio;
  double              p_norm, step_p_norm, p_normratio;
  double              dual_norm, step_dual_norm, dual_normratio;
#endif

#ifdef YMIR_DEBUG
  /* compute norm ratio */
  up_norm = ymir_vec_norm (up);
  step_up_norm = ymir_vec_norm (step_up);
  up_normratio = (up_norm + step_up_norm) / step_up_norm;

  YMIR_GLOBAL_INFOF ("NL iter %d, l^2 norm of step before line search %1.3e\n",
                     nl_iter, step_up_norm);
  YMIR_GLOBAL_INFOF ("NL iter %d, l^2 norm of iterate %1.3e\n",
                     nl_iter, up_norm);
  YMIR_GLOBAL_INFOF ("NL iter %d, l^2 norm ratio: "
                     "(iterate + step)/step = %1.3e\n",
                     nl_iter, up_normratio);
  YMIR_GLOBAL_INFOF ("NL iter %d, l^2 norm ratio: step/iterate = %1.3e\n",
                     nl_iter, step_up_norm / up_norm);

  u_norm = ymir_cvec_norm (up);
  step_u_norm = ymir_cvec_norm (step_up);
  u_normratio = (u_norm + step_u_norm) / step_u_norm;

  YMIR_GLOBAL_INFOF ("NL iter %d, l^2 norm of velocity step before "
                     "line search %1.3e\n",
                     nl_iter, step_u_norm);
  YMIR_GLOBAL_INFOF ("NL iter %d, l^2 norm of velocity %1.3e\n",
                     nl_iter, u_norm);
  YMIR_GLOBAL_INFOF ("NL iter %d, l^2 norm ratio of velocity: "
                     "(iterate + step)/step = %1.3e\n",
                     nl_iter, u_normratio);
  YMIR_GLOBAL_INFOF ("NL iter %d, l^2 norm ratio of velocity: "
                     "step/iterate = %1.3e\n",
                     nl_iter, step_u_norm / u_norm);

  p_norm = ymir_evec_norm (up);
  step_p_norm = ymir_evec_norm (step_up);
  p_normratio = (p_norm + step_p_norm) / step_p_norm;

  YMIR_GLOBAL_INFOF ("NL iter %d, l^2 norm of pressure step before "
                     "line search %1.3e\n",
                     nl_iter, step_p_norm);
  YMIR_GLOBAL_INFOF ("NL iter %d, l^2 norm of pressure %1.3e\n",
                     nl_iter, p_norm);
  YMIR_GLOBAL_INFOF ("NL iter %d, l^2 norm ratio of pressure: "
                     "(iterate + step)/step = %1.3e\n",
                     nl_iter, p_normratio);
  YMIR_GLOBAL_INFOF ("NL iter %d, l^2 norm ratio of pressure: "
                     "step/iterate = %1.3e\n",
                     nl_iter, step_p_norm / p_norm);
#endif

  YMIR_GLOBAL_INFOF ("NL iter %d, descend condition factor %.6e\n", nl_iter,
                     1.0 - descend_cond_relax * (1.0 - k_res_reduction));

  /* store current iterate */
  up_prev = ymir_vec_clone (up);
  if (nl_solver_primaldual_type != SL_NL_SOLVER_PRIMALDUAL_NONE) {
    dual_prev = ymir_vec_clone (dual);
  }
  else {
    dual_prev = NULL;
  }

  /* store vectors for step length reduction */
  switch (step_reduction_type) {
  case SL_NL_SOLVER_STEP_REDUCTION_CONST:
    residual_up_prev = NULL;
    adaptive_step_length = NULL;
    break;
  case SL_NL_SOLVER_STEP_REDUCTION_ADAPTIVE:
    residual_up_prev = ymir_vec_clone (residual_up);
    adaptive_step_length = sc_dmatrix_new (2, nl_stokes->mesh->cnodes->K);
    break;
  case SL_NL_SOLVER_STEP_REDUCTION_QUAD_INTERP:
    residual_up_prev = ymir_vec_clone (residual_up);
    adaptive_step_length = NULL;
    break;
  default: /* unknown step length reduction type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* search for step length */
  for (k = 0; k < step_search_maxiter; k++) {
    /* move iterate into step direction */
    slabs_nl_solver_steplength_move (
        up, dual, step_up, step_dual, residual_dual, up_prev,
        step_length, adaptive_step_length, k,
        nl_stokes, physics_options, solver_options);

    /* update nonlinear operator at new iterate position */
    slabs_nl_solver_steplength_update_operator (
        state, dual, step_up, step_dual, residual_dual,
        step_length, nl_stokes, physics_options, solver_options);

    /* compute residual norm at new iterate position */
    *res_new = slabs_nl_solver_steplength_compute_residual (
        res_new_vel, res_new_press, res_new_dual, residual_up, residual_dual,
        up, dual, nl_stokes, physics_options, solver_options);

    switch (step_reduction_type) {
    case SL_NL_SOLVER_STEP_REDUCTION_CONST:
    case SL_NL_SOLVER_STEP_REDUCTION_QUAD_INTERP:
      YMIR_GLOBAL_INFOF (
          "NL iter %d, step length iter %d, constant length %.3e\n",
          nl_iter, k, step_length);
      break;
    case SL_NL_SOLVER_STEP_REDUCTION_ADAPTIVE:
      if (0 < k) {
        /*TODO
        YMIR_GLOBAL_INFOF (
            "NL iter %d, step length iter %d, adaptive length "
            "min %.3e max %.3e\n", nl_iter, k,
            ymir_vec_min_global (adaptive_step_length),
            ymir_vec_max_global (adaptive_step_length));
        */
      }
      else {
        YMIR_GLOBAL_INFOF (
            "NL iter %d, step length iter %d, adaptive length "
            "min %.3e max %.3e\n", nl_iter, k, 1.0, 1.0);
      }
      break;
    default: /* unknown step length reduction type */
      YMIR_ABORT_NOT_REACHED ();
    }
    YMIR_GLOBAL_INFOF (
        "NL iter %d, step length iter %d, res_new %.3e, res_prev %.3e\n",
        nl_iter, k, *res_new, res_prev);

    /* stop search if descend condition is satisfied */
    if ( (*res_new / res_prev) <=
         (1.0 - descend_cond_relax * (1.0 - k_res_reduction)) ) {
      break;
    }

    /* reduce step length */
    step_reduction_success = slabs_nl_solver_steplength_reduce (
        &step_length, adaptive_step_length, &search_lin_model_slope,
        residual_up, residual_up_prev, *res_new,
        res_prev, res_prev_vel, res_prev_press,
        state, step_up, up_prev, k, nl_stokes, physics_options, solver_options);

    /* stop search if step length reduction was not successful */
    if (!step_reduction_success) {
      break;
    }

    /* restore previous iterate for next iteration of the step length search */
    ymir_vec_copy (up_prev, up);
    if (nl_solver_primaldual_type != SL_NL_SOLVER_PRIMALDUAL_NONE) {
      ymir_vec_copy (dual_prev, dual);
    }
  }

  /* check if descend step was found */
  if (k == step_search_maxiter) { /* if max iterations reached */
    /* reverse updates of Stokes operator (assume iterate is where it
     * started) */
    slabs_nonlinear_stokes_op_update (nl_stokes, state, physics_options,
                                      nl_solver_type);

    /* set step length to invalid value and step to zero */
    step_length = -1.0;
    ymir_vec_set_zero (step_up);
  }
  else {
    /* update 4th order tensor in linearized Stokes operator */
    slabs_nonlinear_stokes_op_update_tensor (
        nl_stokes, state, physics_options, nl_solver_primaldual_type,
        solver_options->nl_solver_primaldual_scal_type,
        solver_options->nl_check_derivative);
  }

#ifdef YMIR_DEBUG
  /* compute norm of the actual step taken */
  YMIR_GLOBAL_INFOF ("NL iter %d, l^2 norm of step %1.3e\n",
                     nl_iter, step_length * step_up_norm);
  YMIR_GLOBAL_INFOF ("NL iter %d, l^2 norm of velocity step %1.3e\n",
                     nl_iter, step_length * step_u_norm);
  YMIR_GLOBAL_INFOF ("NL iter %d, l^2 norm of pressure step %1.3e\n",
                     nl_iter, step_length * step_p_norm);

  /* compute norm ratio of dual tensor */
  if (nl_solver_primaldual_type != SL_NL_SOLVER_PRIMALDUAL_NONE) {
    dual_norm = ymir_vec_norm (dual_prev);
    step_dual_norm = ymir_vec_norm (step_dual);
    dual_normratio = (dual_norm + step_dual_norm) / step_dual_norm;

    YMIR_GLOBAL_INFOF ("NL iter %d, l^2 norm of dual step %1.3e\n",
                       nl_iter, step_dual_norm);
    YMIR_GLOBAL_INFOF ("NL iter %d, l^2 norm of dual iterate %1.3e\n",
                       nl_iter, dual_norm);
    YMIR_GLOBAL_INFOF ("NL iter %d, l^2 norm ratio of dual: "
                       "(iterate + step)/step = %1.3e\n",
                       nl_iter, dual_normratio);
    YMIR_GLOBAL_INFOF ("NL iter %d, l^2 norm ratio of dual: "
                       "step/iterate = %1.3e\n",
                       nl_iter, step_dual_norm / dual_norm);
  }

  /* compute slope of linear model in step direction with updated Stokes op. */
  *lin_model_slope =
    slabs_nl_solver_compute_lin_model_slope (stokes_op, step_up, residual_up,
                                             *res_new, solver_options);

  YMIR_GLOBAL_INFOF ("NL iter %d, linear model slope at new "
                     "iterate %1.3e\n", nl_iter, *lin_model_slope);
#else
  *lin_model_slope = 0.0;
#endif

  /* destroy */
  ymir_vec_destroy (up_prev);
  if (dual_prev != NULL) {
    ymir_vec_destroy (dual_prev);
  }
  if (residual_up_prev != NULL) {
    ymir_vec_destroy (residual_up_prev);
  }
  if (adaptive_step_length != NULL) {
    sc_dmatrix_destroy (adaptive_step_length);
  }

  /* return step length */
  return step_length;
}

#if 0
/**
 *
 */
static void
slabs_nl_solver_bbt_apply (ymir_vec_t *in, ymir_vec_t *out, void *data)
{
  ymir_stokes_pc_t   *stokes_pc = (ymir_stokes_pc_t *) data;
  ymir_vel_dir_t     *vel_dir = stokes_pc->stokes_op->stress_op->vel_dir;
  ymir_pressure_elem_t  *press_elem = stokes_pc->stokes_op->press_elem;
  ymir_vec_t         *tmp = ymir_cvec_new_zero (in->mesh, 3);
  ymir_vec_t         *bfbt_uscale = stokes_pc->bfbt_inner_uscale;

  ymir_pressure_vec_weak_grad (in, tmp, press_elem, vel_dir, 0);

  if (bfbt_uscale != NULL) {
    ymir_vec_multiply_in (bfbt_uscale, tmp);
  }

  ymir_pressure_vec_weak_grad (out, tmp, press_elem, vel_dir, 1);

  ymir_vec_destroy (tmp);
}

/**
 * Computes initial guess for nonlinear solver:
 *
 *   u = f / visc
 *   p = (B * B^T)^(-1) * B * (f - A * u)
 *
 * where
 *
 *   f --- right-hand side
 *   A --- elliptic operator
 *   B --- gradient operator
 */
static void
slabs_nl_solver_set_initial_guess_from_rhs (
    ymir_vec_t *init_up,
    slabs_nl_stokes_problem_t *nl_stokes,
    ymir_vec_t *tmp_up,
    const slabs_nl_solver_initial_guess_t initial_guess_type)
{
  ymir_mesh_t        *mesh = nl_stokes->mesh;
  ymir_pressure_elem_t  *press_elem = nl_stokes->press_elem;
  ymir_stokes_op_t   *stokes_op = nl_stokes->stokes_op;
  ymir_stokes_pc_t   *stokes_pc = nl_stokes->stokes_pc;
  ymir_vel_dir_t     *vel_dir = nl_stokes->vel_dir;
  ymir_cvec_t        *rhs_u_point = nl_stokes->rhs_u_point;
  ymir_dvec_t        *viscosity = stokes_op->stress_op->viscosity;
  MPI_Comm            mpicomm = mesh->ma->mpicomm;

  ymir_cvec_t        *viscosity_cont = ymir_cvec_new (mesh, 1);
  ymir_cvec_t        *init_u, *tmp_u;
  ymir_evec_t        *init_p, *tmp_p;

  ymir_matrix_t      *bbt = NULL;
#ifdef YMIR_PETSC
  MatNullSpace        nsp;
  Mat                 bbt_shell;
  Mat                 bbt_mat = PETSC_NULL;
  KSP                 bbt_ksp;
  PetscErrorCode      ierr;
#endif

  /* create vector views */
  slabs_stokes_vec_get_components_view (&init_u, &init_p, init_up);
  slabs_stokes_vec_get_components_view (&tmp_u, &tmp_p, tmp_up);

  /* set initial velocity: u = f/visc
   * Assume: velocity and pressure was set to zero! */
  ymir_cvec_copy (rhs_u_point, init_u);
  ymir_interp_vec (viscosity, viscosity_cont);
  ymir_cvec_fabs (viscosity_cont, viscosity_cont);
  ymir_cvec_divide_in1 (viscosity_cont, init_u);
  ymir_vec_destroy (viscosity_cont);

  /* done for zero pressure initial guess */
  if (initial_guess_type == SL_NL_SOLVER_INITIAL_GUESS_URHS_PZERO) {
    ymir_vec_destroy (init_u);
    ymir_vec_destroy (init_p);
    ymir_vec_destroy (tmp_u);
    ymir_vec_destroy (tmp_p);
    return;
  }

  /* compute residual (f - A * u) */
  {
    ymir_vec_t         *rhs_up_point = ymir_stokes_vec_new (mesh, press_elem);
    ymir_vec_t         *rhs;

    ymir_vec_set_zero (rhs_up_point);
    ymir_stokes_vec_set_velocity (rhs_u_point, rhs_up_point, press_elem);
    rhs = ymir_stokes_vec_new (mesh, press_elem);
    ymir_stokes_op_construct_rhs (rhs_up_point, NULL, NULL, rhs, stokes_op);

    ymir_stokes_op_apply (init_up, tmp_up, stokes_op);
    ymir_vec_add (-1.0, rhs, tmp_up);
    ymir_vec_scale (-1.0, tmp_up);
    ymir_evec_set_zero (tmp_p);

    ymir_vec_destroy (rhs);
    ymir_vec_destroy (rhs_up_point);
  }

  /* apply divergence operator to residual: B * (f - A * u) */
  ymir_pressure_vec_weak_grad_dirty (tmp_p, tmp_u, press_elem, vel_dir, 1);

#ifdef YMIR_PETSC
  /* set pressure Laplacian apply operator */
  ierr = MatShellCreate_ymir (slabs_nl_solver_bbt_apply,
                              (void *) stokes_pc, init_p, init_p, &bbt_shell);
  YMIR_CHECK_PETSC (ierr);

  /* assemble pressure Laplacian matrix */
  ymir_pressure_laplacian_matrix_c2e (vel_dir, NULL, press_elem, &bbt, 0);

  bbt_mat = bbt->p;
  ierr = PetscObjectReference ((PetscObject) bbt_mat);
  YMIR_CHECK_PETSC (ierr);

  /* set null space of pressure gradient operator */
  /*
  if (stokes_op->press_elem->space == YMIR_PRESSURE_SPACE_POLY) {
    ierr = MatNullSpaceCreate (mpicomm, 0, 1, &(stokes_pc->bbt->nvec), &nsp);
    YMIR_CHECK_PETSC (ierr);
  }
  else {
    ierr = MatNullSpaceCreate (mpicomm, 1, 0, PETSC_NULL, &nsp);
    YMIR_CHECK_PETSC (ierr);
  }
  ierr = MatSetNearNullSpace (bbt_mat, nsp);
  YMIR_CHECK_PETSC (ierr);
  */

  /* create PETSc KSP operator for pressure Laplacian (B * B^T) */
  ierr = KSPCreate (mpicomm, &bbt_ksp);
  YMIR_CHECK_PETSC (ierr);
  ierr = KSPSetOptionsPrefix (bbt_ksp, "initial_guess_bbt_");
  YMIR_CHECK_PETSC (ierr);
  ierr = KSPSetOperators (bbt_ksp, bbt_shell, bbt_mat);
  YMIR_CHECK_PETSC (ierr);
  ierr = KSPSetFromOptions (bbt_ksp);
  YMIR_CHECK_PETSC (ierr);
  ierr = KSPSetUp (bbt_ksp);
  YMIR_CHECK_PETSC (ierr);

  /* destroy */
  ierr = MatDestroy (&bbt_shell);    YMIR_CHECK_PETSC (ierr);
  ierr = MatDestroy (&bbt_mat);      YMIR_CHECK_PETSC (ierr);
  ierr = MatNullSpaceDestroy (&nsp); YMIR_CHECK_PETSC (ierr);

  /* solve for initial pressure: p = (B * B^T)^(-1) * B * (f - A * u) */
  ymir_petsc_ksp_solve (init_p, tmp_p, bbt_ksp, 0, YMIR_MESH_PETSCLAYOUT_NONE,
                        NULL);

  /* destroy KSP operator for pressure Laplacian */
  ierr = KSPDestroy (&bbt_ksp); YMIR_CHECK_PETSC (ierr);
  ymir_matrix_destroy (bbt);
#endif

  /* destroy vector views */
  ymir_vec_destroy (init_u);
  ymir_vec_destroy (init_p);
  ymir_vec_destroy (tmp_u);
  ymir_vec_destroy (tmp_p);
}
#endif

/**
 * Solves nonlinear Stokes problem with Picard and/or Newton's method.
 */
static void
slabs_nonlinear_solver_picard_newton (slabs_stokes_state_t *state,
                                      slabs_nl_stokes_problem_t **nl_stokes,
                                      slabs_physics_options_t *physics_options,
                                      slabs_discr_options_t *discr_options,
                                      slabs_nl_solver_options_t *solver_options,
                                      const char *bin_filepath,
                                      const char *vtk_filepath)
{
  const char         *this_fn_name = "slabs_nonlinear_solver_picard_newton";
  int                 mpiret;
  /* nonlinear solver variables */
  const slabs_nl_solver_type_t  nl_solver_type = solver_options->nl_solver_type;
  const slabs_nl_solver_primaldual_type_t  nl_solver_primaldual_type =
                        solver_options->nl_solver_primaldual_type;
  const slabs_nl_solver_primaldual_scal_type_t  nl_solver_primaldual_scal_type =
                        solver_options->nl_solver_primaldual_scal_type;
  const int           nl_maxiter = solver_options->nl_maxiter;
  const double        nl_rtol = solver_options->nl_rtol;
  const slabs_nl_solver_initial_guess_t  initial_guess_type =
                        solver_options->initial_guess_type;
  const int           nl_resume_at_iter =
                        SC_MAX (solver_options->nl_resume_at_iter, 0);
  int                 nl_iter;
  int                 newton_iter = 0;
  int                 picard_iter = 0;
  int                 use_picard;
  int                 picard_maxiter = 0;
  const double        picard_rtol = solver_options->nl_switch_picard_rtol;
  const int           switch_picard_maxiter =
                        solver_options->nl_switch_picard_maxiter;
  const int           switch_picard_init =
                        solver_options->nl_switch_picard_init;
  const double        switch_picard_step_length_min =
                        solver_options->nl_switch_picard_step_length_min;
  const int           switch_picard_after_amr =
                        solver_options->nl_switch_picard_after_amr;
  const slabs_nl_stokes_problem_schur_diag_t  schur_diag_type =
                        solver_options->schur_diag_type;
  const slabs_nl_stokes_problem_scaling_t  scaling_type =
                        solver_options->scaling_type;
  const double        grid_cont_init_threshold =
                        solver_options->grid_continuation_init_amr_threshold;
  const int           grid_cont_init_steps =
                        solver_options->grid_continuation_init_steps;
  const int           grid_cont_skipsteps =
                        solver_options->grid_continuation_skipsteps;
  const int           grid_cont_maxsteps =
                        solver_options->grid_continuation_maxsteps;
  int                 grid_cont_steps_skipped = 0;
  const double        visc_min = physics_options->viscosity_min;
  const double        visc_max = physics_options->viscosity_max;
  const double        visc_bounds_cont_min =
                        solver_options->viscosity_bounds_continuation_min;
  const double        visc_bounds_cont_max =
                        solver_options->viscosity_bounds_continuation_max;
  const int           visc_bounds_cont_steps =
                        solver_options->viscosity_bounds_continuation_steps;
  double              lin_model_slope = 0.0;
  const int           check_derivative = solver_options->nl_check_derivative;
  const int           log_physics_stats = solver_options->log_physics_stats;
  const int           project_out_nullspace =
                        solver_options->project_out_nullspace;
  /* linear solver variables */
  double              lin_res_reduction;
  int                 lin_num_iter;
  /* discretization variables */
  ymir_mesh_t        *mesh = (*nl_stokes)->mesh;
  ymir_pressure_elem_t  *press_elem = (*nl_stokes)->press_elem;
  slabs_discr_amr_indicator_params_t  *indicator_params;
  const double        amr_threshold = discr_options->amr_rel_threshold;
  int                 run_amr = 0;
  int                 mesh_modified = 0;
  /* vectors and related variables */
  ymir_vec_t         *up;
  ymir_dvec_t        *dual;
  ymir_vec_t         *residual_up, *residual_lump_up;
  ymir_dvec_t        *residual_dual;
  double              res_new;
  double              res_new_vel = 0.0;
  double              res_new_press = 0.0;
  double              res_new_dual = 0.0;
  double              res_prev, res_prev_vel, res_prev_press;
  double              res_init;
  double              res_reduction_prev;
  double              res_new_last_mesh, res_new_last_mesh_vel,
                      res_new_last_mesh_press, res_new_last_mesh_dual;
  double              picard_res_init;
  ymir_vec_t         *step_up;
  ymir_dvec_t        *step_dual;
  double              step_length;
  int                 destroy_dual;
  /* miscellaneous variables */
  double              time_start, time;

  /*
   * Initialization
   */

  YMIR_GLOBAL_PRODUCTIONF ("Into %s: rtol %1.3e, maxiter %d\n",
                           this_fn_name, nl_rtol, nl_maxiter);

  /* record start time */
  mpiret = MPI_Barrier (mesh->ma->mpicomm); YMIR_CHECK_MPI (mpiret);
  time_start = MPI_Wtime () - SC_MAX (solver_options->nl_resume_at_time, 0.0);

  /* initialize AMR parameters from options */
  if (discr_options->amr_visc_indicator_type != SL_AMR_INDICATOR_NONE ||
      discr_options->amr_visc_dr_indicator_type != SL_AMR_INDICATOR_NONE ||
      discr_options->amr_strain_rate_indicator_type != SL_AMR_INDICATOR_NONE) {
    int                 n_indicators = 0;
    int                 indicator_idx = 0;

    /* allocate memory for AMR indicator parameters */
    if (discr_options->amr_visc_indicator_type != SL_AMR_INDICATOR_NONE) {
      n_indicators++;
    }
    if (discr_options->amr_visc_dr_indicator_type != SL_AMR_INDICATOR_NONE) {
      n_indicators++;
    }
    if (   discr_options->amr_strain_rate_indicator_type
        != SL_AMR_INDICATOR_NONE ) {
      n_indicators++;
    }
    indicator_params = slabs_discr_amr_indicator_params_new (n_indicators);

    /* set viscosity AMR parameters */
    if (discr_options->amr_visc_indicator_type != SL_AMR_INDICATOR_NONE) {
      indicator_params->type[indicator_idx] =
        discr_options->amr_visc_indicator_type;
      indicator_params->tol_min[indicator_idx] =
        discr_options->amr_visc_tol_min;
      indicator_params->tol_max[indicator_idx] =
        discr_options->amr_visc_tol_max;
      indicator_params->level_min[indicator_idx] = 0;
      indicator_params->level_max[indicator_idx] = 0;
      indicator_idx++;
    }

    /* set viscosity dynamic range AMR parameters */
    if (discr_options->amr_visc_dr_indicator_type != SL_AMR_INDICATOR_NONE) {
      indicator_params->type[indicator_idx] =
        discr_options->amr_visc_dr_indicator_type;
      indicator_params->tol_min[indicator_idx] =
        discr_options->amr_visc_dr_tol_min;
      indicator_params->tol_max[indicator_idx] =
        discr_options->amr_visc_dr_tol_max;
      indicator_params->level_min[indicator_idx] = 0;
      indicator_params->level_max[indicator_idx] = 0;
      indicator_idx++;
    }

    /* set strain rate AMR parameters */
    if (   discr_options->amr_strain_rate_indicator_type
        != SL_AMR_INDICATOR_NONE ) {
      indicator_params->type[indicator_idx] =
        discr_options->amr_strain_rate_indicator_type;
      indicator_params->tol_min[indicator_idx] =
        discr_options->amr_strain_rate_tol_min;
      indicator_params->tol_max[indicator_idx] =
        discr_options->amr_strain_rate_tol_max;
      indicator_params->level_min[indicator_idx] = 0;
      indicator_params->level_max[indicator_idx] = 0;
    }
  }
  else {
    indicator_params = NULL;
  }

  /* viscosity bounds continuation: set initial bounds */
  if (nl_resume_at_iter < visc_bounds_cont_steps) {
    if (0.0 < visc_bounds_cont_min && visc_min < visc_bounds_cont_min) {
      physics_options->viscosity_min = visc_bounds_cont_min;

      YMIR_GLOBAL_INFOF ("%s: Set viscosity bound min %1.3e\n",
                         this_fn_name, physics_options->viscosity_min);
    }
    if (0.0 < visc_bounds_cont_max && visc_bounds_cont_max < visc_max) {
      physics_options->viscosity_max = visc_bounds_cont_max;

      YMIR_GLOBAL_INFOF ("%s: Set viscosity bound max %1.3e\n",
                         this_fn_name, physics_options->viscosity_max);
    }
  }

  /* check Stokes state */
  YMIR_ASSERT (state->temp_vec != NULL);
  YMIR_ASSERT (state->vel_press_vec != NULL);
  if (nl_solver_primaldual_type != SL_NL_SOLVER_PRIMALDUAL_NONE) {
    if (state->dual_vec == NULL) {
      slabs_stokes_state_init_dual_tensor (state, mesh);
      destroy_dual = 1;
    }
    else {
      destroy_dual = 0;
    }
  }
  else {
    destroy_dual = 0;
  }

  /* create vectors */
  up = state->vel_press_vec;
  dual = state->dual_vec;
  residual_up = ymir_stokes_vec_new (mesh, press_elem);
  if (vtk_filepath != NULL) {
    residual_lump_up = ymir_stokes_vec_new (mesh, press_elem);
  }
  else {
    residual_lump_up = NULL;
  }
  step_up = ymir_stokes_vec_new (mesh, press_elem);
  if (nl_solver_primaldual_type != SL_NL_SOLVER_PRIMALDUAL_NONE) {
    residual_dual = ymir_vec_template (dual);
    step_dual = ymir_vec_template (dual);
  }
  else {
    residual_dual = NULL;
    step_dual = NULL;
  }

  /* set initial guess */
  switch (initial_guess_type) {
  case SL_NL_SOLVER_INITIAL_GUESS_ZERO:
    YMIR_GLOBAL_INFOF ("%s: Set initial guess to zero\n", this_fn_name);
    ymir_vec_set_zero (up);
    if (nl_solver_primaldual_type != SL_NL_SOLVER_PRIMALDUAL_NONE) {
      ymir_vec_set_zero (dual);
    }
    break;
  case SL_NL_SOLVER_INITIAL_GUESS_NONZERO:
    YMIR_GLOBAL_INFOF ("%s: Use nonzero initial guess\n", this_fn_name);
    if (nl_solver_primaldual_type != SL_NL_SOLVER_PRIMALDUAL_NONE) {
      ymir_vec_set_zero (dual); //TODO don't know what else to do here
    }
    break;
  default: /* unknown type of initial guess */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* create a new Stokes operator */
  slabs_nonlinear_stokes_op_init (
      *nl_stokes, state, physics_options, 1, SL_NL_SOLVER_NEWTON,
      nl_solver_primaldual_type, nl_solver_primaldual_scal_type,
      check_derivative);

  /* create a new Stokes preconditioner */
  slabs_nonlinear_stokes_pc_init (*nl_stokes, state, schur_diag_type,
                                  scaling_type, (0 == nl_resume_at_iter));

  /* create H^-1 norm operator */
  if (solver_options->norm_type == SL_NORM_FNC_HMINUS1_L2) {
    solver_options->norm_op = ymir_Hminus1_norm_op_new (
        mesh, solver_options->norm_Hminus1_mass_scaling);
  }
  else {
    solver_options->norm_op = NULL;
  }

  /* compute initial residual */
  res_new = slabs_nl_solver_compute_residual (
      residual_up, residual_lump_up, residual_dual,
      &res_new_vel, &res_new_press, &res_new_dual,
      up, dual, *nl_stokes, solver_options);
  res_prev = res_new;
  res_prev_vel = res_new_vel;
  res_prev_press = res_new_press;
  res_init = res_new;
  res_reduction_prev = 1.0;
  res_new_last_mesh = res_new;
  res_new_last_mesh_vel = res_new_vel;
  res_new_last_mesh_press = res_new_press;
  res_new_last_mesh_dual = res_new_dual;
  if (0 < nl_resume_at_iter) {
    if (0.0 < solver_options->nl_resume_prev_res) {
      res_prev = solver_options->nl_resume_prev_res;
    //res_prev_vel = TODO
    //res_prev_press = TODO
    }
    if (0.0 < solver_options->nl_resume_init_res) {
      res_init = solver_options->nl_resume_init_res;
    }
  }

  /* set type of nonlinear solver that we start with */
  if (switch_picard_init) {
    /* activate Picard for first step */
    use_picard = 1;
    picard_maxiter += switch_picard_maxiter;
    picard_res_init = res_new;

    /* modify Stokes operator for Picard */
    slabs_nonlinear_stokes_op_switch_picard (*nl_stokes);
  }
  else {
    /* deactivate Picard for first step */
    use_picard = 0;
    picard_res_init = 0.0;
  }

  /* set initial relative AMR threshold */
  if (0 < grid_cont_init_steps) {
    discr_options->amr_rel_threshold = grid_cont_init_threshold;
  }

  /*
   * Run Nonlinear Solver
   */

  for (nl_iter = nl_resume_at_iter; nl_iter <= nl_maxiter; nl_iter++) {

    /* output of elapsed time */
    mpiret = MPI_Barrier (mesh->ma->mpicomm); YMIR_CHECK_MPI (mpiret);
    time = MPI_Wtime ();
    YMIR_GLOBAL_INFOF ("NL iter %d, elapsed time %.16e sec\n", nl_iter,
                       time - time_start);

    /* output of current residual norm */
    YMIR_GLOBAL_INFOF ("NL iter %d, residual %.6e\n", nl_iter, res_new);
    YMIR_GLOBAL_INFOF ("NL iter %d, residual mom %.3e, mass %.3e, dual %.3e\n",
                       nl_iter, res_new_vel, res_new_press, res_new_dual);

    /* output of residual norm before AMR */
    if (mesh_modified && res_new_last_mesh != res_prev) {
      YMIR_GLOBAL_VERBOSEF ("NL iter %d, residual last mesh %.3e, "
                            "mom %.3e, mass %.3e, dual %.3e\n",
                            nl_iter, res_new_last_mesh, res_new_last_mesh_vel,
                            res_new_last_mesh_press, res_new_last_mesh_dual);
    }

    /* output of residual norm at previous and initial nonlinear iteration */
    YMIR_GLOBAL_VERBOSEF ("NL iter %d, residual %.3e, prev %.3e, init %.3e\n",
                          nl_iter, res_new, res_prev, res_init);

    /* output of l^2 residual norm */
    {
      double              nrm, nrm_u, nrm_p;

      nrm = slabs_norm_of_residual (&nrm_u, &nrm_p, residual_up,
                                    (*nl_stokes)->stokes_op,
                                    (*nl_stokes)->stokes_pc,
                                    solver_options->krylov_type,
                                    SL_NORM_VEC_L2, NULL);
      YMIR_GLOBAL_INFOF ("NL iter %d, residual l^2 norm %.3e, "
                         "mom %.3e, mass %.3e\n",
                         nl_iter, nrm, nrm_u, nrm_p);

      /* set absolute Krylov tolerance based on initial residual
       * to avoid oversolving */
      if ( 0.0 < solver_options->nl_forcing_total_min &&
           solver_options->krylov_atol <= 0.0 ) {
        solver_options->krylov_atol = nrm *
                                      solver_options->nl_forcing_total_min;
        YMIR_GLOBAL_INFOF ("%s: Set absolute Krylov tolerance to %.3e "
                           "(init residual l^2 norm %.3e)\n",
                           this_fn_name, solver_options->krylov_atol, nrm);
      }
    }

    /* output of other residual norms */
#if 0
    {
      double              nrm, nrm_u, nrm_p;

      nrm = slabs_norm_of_residual (&nrm_u, &nrm_p, residual_up,
                                    (*nl_stokes)->stokes_op,
                                    (*nl_stokes)->stokes_pc,
                                    solver_options->krylov_type,
                                    SL_NORM_FNC_L2, NULL);
      YMIR_GLOBAL_INFOF ("NL iter %d, residual L^2 norm %1.3e, "
                         "vel %1.3e, press %1.3e\n",
                         nl_iter, nrm, nrm_u, nrm_p);

      if (solver_options->norm_op != NULL) {
        nrm = slabs_norm_of_residual_comp (&nrm_u, &nrm_p, residual_up,
                                           (*nl_stokes)->stokes_op,
                                           (*nl_stokes)->stokes_pc,
                                           solver_options->krylov_type,
                                           SL_NORM_FNC_HMINUS1_L2,
                                           solver_options->norm_op);
        YMIR_GLOBAL_INFOF ("NL iter %d, residual (H^-1,L^2) norm %1.3e, "
                           "vel %1.3e, press %1.3e\n",
                           nl_iter, nrm, nrm_u, nrm_p);
      }
    }
#endif

    /* output of physics statistics */
    if (log_physics_stats) {
      /* print all major physics statistics */
      {
        ymir_vec_t         *viscosity_nondim =
          ymir_vec_clone ((*nl_stokes)->stokes_op->stress_op->viscosity);
        ymir_vec_t         *yielding_marker = (*nl_stokes)->yielding_marker;

        ymir_vec_scale (0.5, viscosity_nondim);
        slabs_physics_stats_print_all (up, viscosity_nondim, yielding_marker,
                                       (*nl_stokes)->rhs_u_point,
                                       (*nl_stokes)->stokes_op,
                                       physics_options);
        ymir_vec_destroy (viscosity_nondim);
      }

      /* print dual tensor statistics */
      if (nl_solver_primaldual_type != SL_NL_SOLVER_PRIMALDUAL_NONE) {
        double              min, max;

        min = slabs_norm_symtensor_frobenius_min (dual);
        max = slabs_norm_symtensor_frobenius_max (dual);
        YMIR_GLOBAL_STATISTICSF (
            "%s: Range of dual:            min %.3e, max %.3e, max/min %.3e\n",
            this_fn_name, min, max, max/min);
      }
    }

    /* write binary file */
    if (bin_filepath != NULL) {
      char            path[BUFSIZ];

      snprintf (path, BUFSIZ, "%s_nl_itn%02d", bin_filepath, nl_iter);
      slabs_stokes_state_save (state, path);

      /* create options to resume nonlinear iterations */
      if (mesh->ma->mpirank == 0) {
        FILE               *fp;
        int                 ret;

        snprintf (path, BUFSIZ, "%s_nl_resume.ini", bin_filepath);
        fp = fopen (path, "w");
        YMIR_CHECK_ABORT (fp != NULL, "file open");

        fprintf (fp, "[Mantle]\n");
        fprintf (fp, "  p4est-import = %s_nl_itn%02d_mesh\n",
                 bin_filepath, nl_iter);
        fprintf (fp, "  temperature = %s_nl_itn%02d_temperature.bin\n",
                 bin_filepath, nl_iter);
        fprintf (fp, "  velocity-import = %s_nl_itn%02d_velocity.bin\n",
                 bin_filepath, nl_iter);
        fprintf (fp, "  pressure-import = %s_nl_itn%02d_pressure.bin\n",
                 bin_filepath, nl_iter);
        fprintf (fp, "  viscosity-for-init-nl-stokes = 0\n");
        fprintf (fp, "  refine-surface-max-dist = 0\n");
        fprintf (fp, "  refine-lower-upper-interface-max-dist = 0\n");
        fprintf (fp, "  init-amr-visc-indicator = 0\n");
        fprintf (fp, "  init-amr-weak-import-indicator = 0\n");
        fprintf (fp, "  init-amr-rhs-indicator = 0\n");
        fprintf (fp, "  nonlinear-initial-guess-type = %i\n",
                 SL_NL_SOLVER_INITIAL_GUESS_NONZERO);
        fprintf (fp, "  nonlinear-resume-at-iter = %i\n", nl_iter);
        fprintf (fp, "  nonlinear-resume-at-time = %.16e\n", time-time_start);
        fprintf (fp, "  nonlinear-resume-prev-res = %.16e\n", res_prev);
        fprintf (fp, "  nonlinear-resume-init-res = %.16e\n", res_init);
        fprintf (fp, "  krylov-atol = %.16e\n", solver_options->krylov_atol);

        ret = fclose (fp);
        SC_CHECK_ABORT (ret == 0, "file close");
      }
    }

    /* write vtk file */
    if (vtk_filepath != NULL) {
      char            path[BUFSIZ];

      snprintf (path, BUFSIZ, "%s_nl_itn%02d", vtk_filepath, nl_iter);
      slabs_vtk_write_iteration_nl_stokes (
          path, state, step_up, step_dual, residual_lump_up, residual_dual,
          *nl_stokes, physics_options, discr_options);
    }

    /* check for convergence and exit conditions */
    if (0.0 < nl_rtol && res_new < (nl_rtol * res_init)) {
      YMIR_GLOBAL_PRODUCTIONF (
          "%s: Nonlinear problem converged to rtol\n", this_fn_name);
      break;
    }
    if (nl_maxiter <= nl_iter) {
      YMIR_GLOBAL_PRODUCTIONF (
          "%s: Maximum number of nonlinear iterations reached\n", this_fn_name);
      break;
    }

    /* calculate the forcing, i.e., the Krylov relative tolerance */
    slabs_nl_solver_compute_forcing (nl_iter, res_new, res_prev,
                                     res_reduction_prev, solver_options);

    /* update previous residual norms to current value */
    res_reduction_prev = res_new / res_prev;
    res_prev = res_new;
    res_prev_vel = res_new_vel;
    res_prev_press = res_new_press;

    /* setup the Stokes preconditioner */
    if (nl_resume_at_iter < nl_iter) { /* if not initial nonlinear iteration */
      if (0 == run_amr) { /* if no AMR in previous step */
        /* update old preconditioner */
        slabs_nonlinear_stokes_pc_update (*nl_stokes, state, schur_diag_type,
                                          scaling_type);
      }
      else { /* if AMR was performed in previous step */
        /* create a new preconditioner */
        slabs_nonlinear_stokes_pc_init (*nl_stokes, state, schur_diag_type,
                                        scaling_type, 0);
      }
    }

    /* print memory usage */
    ymir_monitor_print_global_mem_usage (mesh->ma->mpicomm);

    if (use_picard) {
      /*
       * Picard Step
       */

      picard_iter++;

      /* solve for Picard step `step_up` */
      lin_num_iter = slabs_nl_solver_picard_step (step_up, &lin_res_reduction,
                                                  residual_up, *nl_stokes,
                                                  solver_options);

      /* project out nullspaces of step */
      if (project_out_nullspace == SL_NL_SOLVER_PROJECT_NULLSPACE ||
          project_out_nullspace == SL_NL_SOLVER_PROJECT_NULLSPACE_SYMM) {
        slabs_linear_stokes_problem_project_out_nullspace (
            step_up, (*nl_stokes)->stokes_op, physics_options, 0);
      }

      /* set full step (no backtracking) */
      step_length = 1.0;
      res_new_last_mesh = res_new;

      /* compute dual step and update dual iterate */
      if (nl_solver_primaldual_type != SL_NL_SOLVER_PRIMALDUAL_NONE) {
        slabs_nl_solver_picard_dual_step (step_dual, dual,
                                          step_length, up, step_up,
                                          residual_dual, *nl_stokes,
                                          physics_options, solver_options);
      }

      /* update (velocity,pressure) iterate */
      ymir_vec_add (step_length, step_up, up);

      YMIR_GLOBAL_INFOF ("NL iter %d, Picard step, step length %g, "
                         "%d Krylov iter for reduction %1.3e "
                         "(prescribed rtol %1.3e)\n",
                         nl_iter, step_length, lin_num_iter,
                         lin_res_reduction, solver_options->krylov_rtol);

      /* decide to switch off Picard for next step */
      if (picard_iter == picard_maxiter) { /* if max num Picard iter reached */
        use_picard = 0;
      }
      else if (res_prev < picard_rtol * picard_res_init) { /* if residual below
                                                            * Picard rtol */
        use_picard = 0;
        picard_maxiter = picard_iter;
      }

      /* update Stokes operator */
      if (use_picard) { /* if next step will be Picard */
        slabs_nonlinear_stokes_op_update (*nl_stokes, state, physics_options,
                                          SL_NL_SOLVER_PICARD);
      }
      else { /* if next step will be Newton */
        YMIR_GLOBAL_INFOF ("%s: Switch to Newton\n", this_fn_name);

        slabs_nonlinear_stokes_op_update (*nl_stokes, state, physics_options,
                                          SL_NL_SOLVER_NEWTON);
        slabs_nonlinear_stokes_op_update_tensor (
            *nl_stokes, state, physics_options, nl_solver_primaldual_type,
            nl_solver_primaldual_scal_type, check_derivative);
      }
    }
    else {
      /*
       * Newton Step
       */

      newton_iter++;

      /* solve for Newton step `step_up` */
      lin_num_iter = slabs_nl_solver_newton_step (step_up, &lin_res_reduction,
                                                  residual_up, *nl_stokes,
                                                  solver_options);

      /* project out nullspaces of step */
      if (project_out_nullspace == SL_NL_SOLVER_PROJECT_NULLSPACE ||
          project_out_nullspace == SL_NL_SOLVER_PROJECT_NULLSPACE_SYMM) {
        slabs_linear_stokes_problem_project_out_nullspace (
            step_up, (*nl_stokes)->stokes_op, physics_options, 0);
      }

      /* update (velocity,pressure) iterate; perform line search to reduce
       * energy associated with the solution
       * (also updates the Stokes operator) */
      step_length = slabs_nl_solver_steplength_search (
          state, step_up, residual_up, step_dual, residual_dual,
          &res_new_last_mesh, &res_new_last_mesh_vel, &res_new_last_mesh_press,
          &res_new_last_mesh_dual, &lin_model_slope, nl_iter,
          res_prev, res_prev_vel, res_prev_press,
          *nl_stokes, physics_options, solver_options);

      YMIR_GLOBAL_INFOF (
          "NL iter %d, Newton step, step length %g, "
          "%d Krylov iter for reduction %.3e (prescribed rtol %.3e), "
          "Krylov convergence rate %.3e\n",
          nl_iter, step_length, lin_num_iter, lin_res_reduction,
          solver_options->krylov_rtol,
          exp (log (lin_res_reduction) / ((double) lin_num_iter)) );

      /* decide to switch on Picard for next step */
      if ( nl_solver_type == SL_NL_SOLVER_PICARD_NEWTON &&
           (step_length < 0.0 || step_length < switch_picard_step_length_min)
         ) {
        YMIR_GLOBAL_INFOF ("%s: Switch to Picard\n", this_fn_name);

        /* set Picard for next step */
        use_picard = 1;
        picard_maxiter += switch_picard_maxiter;
        picard_res_init = 0.0;

        /* modify Stokes operator for Picard */
        slabs_nonlinear_stokes_op_switch_picard (*nl_stokes);
      }
      else if ( nl_solver_type != SL_NL_SOLVER_PICARD_NEWTON &&
                step_length < 0.0 ) { /* if no progress made */
        YMIR_GLOBAL_PRODUCTIONF ("%s: Descend direction not found\n",
                                 this_fn_name);
        break;
      }
    }

    /*
     * Set Up Stokes Problem for Next Nonlinear Step
     */

    /* decide whether to run AMR or not */
    if ( indicator_params != NULL &&
         (grid_cont_maxsteps < 0 || nl_iter < grid_cont_maxsteps) &&
         (grid_cont_skipsteps <= grid_cont_steps_skipped || 0 == nl_iter) ) {
      run_amr = 1;
    }
    else {
      run_amr = 0;
    }

    /* restore relative AMR threshold */
    if (0 < grid_cont_init_steps && grid_cont_init_steps == nl_iter) {
      discr_options->amr_rel_threshold = amr_threshold;
    }

    /* viscosity bounds continuation: set bounds */
    if (0 < visc_bounds_cont_steps && nl_iter < visc_bounds_cont_steps) {
      if (0.0 < visc_bounds_cont_min && visc_min < visc_bounds_cont_min) {
        double              visc_cont_step;

        visc_cont_step = pow (visc_bounds_cont_min / visc_min,
                              1.0 / ((double) visc_bounds_cont_steps));
        physics_options->viscosity_min =
          visc_bounds_cont_min * pow (visc_cont_step, (double) (-1 - nl_iter));

        YMIR_GLOBAL_INFOF ("%s: Set viscosity bound min %1.3e\n",
                           this_fn_name, physics_options->viscosity_min);
      }
      if (0.0 < visc_bounds_cont_max && visc_bounds_cont_max < visc_max) {
        double              visc_cont_step;

        visc_cont_step = pow (visc_max / visc_bounds_cont_max,
                              1.0 / ((double) visc_bounds_cont_steps));
        physics_options->viscosity_max =
          visc_bounds_cont_max * pow (visc_cont_step, (double) (1 + nl_iter));

        YMIR_GLOBAL_INFOF ("%s: Set viscosity bound max %1.3e\n",
                           this_fn_name, physics_options->viscosity_max);
      }
    }
    else if (0 < visc_bounds_cont_steps && visc_bounds_cont_steps == nl_iter) {
      physics_options->viscosity_min = visc_min;
      physics_options->viscosity_max = visc_max;

      YMIR_GLOBAL_INFOF ("%s: Set original viscosity bounds: min %1.3e, "
                         "max %1.3e\n", this_fn_name,
                         physics_options->viscosity_min,
                         physics_options->viscosity_max);
    }

    /* perform AMR */
    if (run_amr) { /* if run AMR */
      /* destroy Stokes operator, PC, and nonlinear Stokes problem
       * (needs to be done before AMR) */
      slabs_nonlinear_stokes_problem_destroy (*nl_stokes);

      /* destroy H^-1 norm operator */
      ymir_Hminus1_norm_op_destroy (solver_options->norm_op);

      /* perform AMR */
      mesh_modified = slabs_discr_amr (state, &mesh, &press_elem, state->p8est,
                                       indicator_params, physics_options,
                                       discr_options, 0);

      /* create nonlinear Stokes problem */
      *nl_stokes = slabs_nonlinear_stokes_problem_new (state, mesh, press_elem,
                                                       physics_options);

      /* create a new Stokes operator */
      slabs_nonlinear_stokes_op_init (
          *nl_stokes, state, physics_options, 0,
          (use_picard ? SL_NL_SOLVER_PICARD : SL_NL_SOLVER_NEWTON),
          nl_solver_primaldual_type, nl_solver_primaldual_scal_type,
          check_derivative);

      /* create H^-1 norm operator */
      if (solver_options->norm_type == SL_NORM_FNC_HMINUS1_L2) {
        solver_options->norm_op = ymir_Hminus1_norm_op_new (
            mesh, solver_options->norm_Hminus1_mass_scaling);
      }

      /* update all mesh dependent variables */
      if (mesh_modified) {
        grid_cont_steps_skipped = 0;

        /* destroy vectors */
        ymir_vec_destroy (residual_up);
        if (vtk_filepath != NULL) {
          ymir_vec_destroy (residual_lump_up);
        }
        ymir_vec_destroy (step_up);
        if (nl_solver_primaldual_type != SL_NL_SOLVER_PRIMALDUAL_NONE) {
          ymir_vec_destroy (residual_dual);
          ymir_vec_destroy (step_dual);
        }

        /* create vectors */
        up = state->vel_press_vec;
        dual = state->dual_vec;
        residual_up = ymir_stokes_vec_new (mesh, press_elem);
        if (vtk_filepath != NULL) {
          residual_lump_up = ymir_stokes_vec_new (mesh, press_elem);
        }
        step_up = ymir_stokes_vec_new (mesh, press_elem);
        if (nl_solver_primaldual_type != SL_NL_SOLVER_PRIMALDUAL_NONE) {
          residual_dual = ymir_vec_template (dual);
          step_dual = ymir_vec_template (dual);
        }

        /* enforce Dirichlet boundary conditions for interpolated velocity */
        {
          ymir_vec_t         *u = ymir_vec_new_blank ();

          ymir_upvec_get_u (up, u, YMIR_RW);
          ymir_vel_dir_separate (u, NULL, NULL, NULL, (*nl_stokes)->vel_dir);
          ymir_upvec_set_u (up, u, YMIR_SET);
          ymir_vec_destroy (u);
        }

        /* project out nullspaces of interpolated velocity and pressure */
        slabs_linear_stokes_problem_project_out_nullspace (
            up, (*nl_stokes)->stokes_op, physics_options, 0);

        /* decide to switch on Picard for next nonlinear step */
        if (   nl_solver_type == SL_NL_SOLVER_PICARD_NEWTON
            && switch_picard_after_amr ) {
          /* set Picard parameters for next step */
          if (!use_picard) {
            YMIR_GLOBAL_INFOF ("%s: Switch to Picard\n", this_fn_name);

            use_picard = 1;
            picard_res_init = 0.0;
            picard_maxiter += switch_picard_maxiter;
          }
          else {
            picard_res_init = 0.0;
            picard_maxiter = picard_iter + switch_picard_maxiter;
          }

          /* modify Stokes operator for Picard */
          slabs_nonlinear_stokes_op_switch_picard (*nl_stokes);
        }
      }
      else {
        grid_cont_steps_skipped++;
      }
    }
    else { /* if no AMR */
      grid_cont_steps_skipped++;
    }

    /* update residual */
    res_new = slabs_nl_solver_compute_residual (
        residual_up, residual_lump_up, residual_dual,
        &res_new_vel, &res_new_press, &res_new_dual,
        up, dual, *nl_stokes, solver_options);

    /* project out null spaces of new residual */
    if (project_out_nullspace == SL_NL_SOLVER_PROJECT_NULLSPACE_RES ||
        project_out_nullspace == SL_NL_SOLVER_PROJECT_NULLSPACE_SYMM) {
      slabs_linear_stokes_problem_project_out_nullspace (
          residual_up, (*nl_stokes)->stokes_op, physics_options, 1);
    }

    /* store initial residual for set of consecutive Picard steps */
    if (use_picard && picard_res_init <= 0.0) { /* if starting Picard next */
      picard_res_init = res_new;
    }
  } /* END: loop of nonlinear iterations */

  /* enforce Dirichlet boundary conditions for solution velocity */
  {
    ymir_vec_t         *u = ymir_vec_new_blank ();

    ymir_upvec_get_u (up, u, YMIR_RW);
    ymir_vel_dir_separate (u, NULL, NULL, NULL, (*nl_stokes)->vel_dir);
    ymir_upvec_set_u (up, u, YMIR_SET);
    ymir_vec_destroy (u);
  }

  /* project out null spaces of nonlinear solution */
  slabs_linear_stokes_problem_project_out_nullspace (
      up, (*nl_stokes)->stokes_op, physics_options, 0);

  /*
   * Clean Up
   */

  /* hack to avoid SEGV when destroying stress PC (TODO fix) */
  /*TODO delete if ok
  if (indicator_params != NULL && !grid_cont_steps_skipped) { // if ran AMR
    ymir_nlstokes_pc_solve (residual_up, step_up, (*nl_stokes)->stokes_pc, 0,
                            0.1, 0.0, 1, 1, NULL);
  }
  */

  /* destroy vectors */
  if (destroy_dual == 1) {
    slabs_stokes_state_clear_dual_tensor (state);
  }
  ymir_vec_destroy (residual_up);
  if (vtk_filepath != NULL) {
    ymir_vec_destroy (residual_lump_up);
  }
  ymir_vec_destroy (step_up);
  if (nl_solver_primaldual_type != SL_NL_SOLVER_PRIMALDUAL_NONE) {
    ymir_vec_destroy (residual_dual);
    ymir_vec_destroy (step_dual);
  }

  /* destroy H^-1 norm operator */
  ymir_Hminus1_norm_op_destroy (solver_options->norm_op);

  /* destroy AMR parameters */
  slabs_discr_amr_indicator_params_destroy (indicator_params);

  /* restore relative AMR threshold */
  if (0 < grid_cont_init_steps) {
    discr_options->amr_rel_threshold = amr_threshold;
  }
}

/**
 * Solves nonlinear Stokes problem.
 */
static void
slabs_nonlinear_solver_solve_state (slabs_stokes_state_t *state,
                                    slabs_nl_stokes_problem_t **nl_stokes,
                                    slabs_physics_options_t *physics_options,
                                    slabs_discr_options_t *discr_options,
                                    slabs_nl_solver_options_t *solver_options,
                                    const char *bin_filepath,
                                    const char *vtk_filepath)
{
  const double        switch_picard_step_length_min =
                        solver_options->nl_switch_picard_step_length_min;
  const int           switch_picard_after_amr =
                        solver_options->nl_switch_picard_after_amr;
  const int           switch_picard_init =
                        solver_options->nl_switch_picard_init;
  const int           switch_picard_maxiter =
                        solver_options->nl_switch_picard_maxiter;
  const double        switch_picard_rtol =
                        solver_options->nl_switch_picard_rtol;

  /* change solver options according to type of nonlinear solver */
  switch (solver_options->nl_solver_type) {
  case SL_NL_SOLVER_PICARD:
    solver_options->nl_switch_picard_init = 1;
    solver_options->nl_switch_picard_maxiter = solver_options->nl_maxiter;
    solver_options->nl_switch_picard_rtol = solver_options->nl_rtol;
    break;

  case SL_NL_SOLVER_NEWTON:
    solver_options->nl_switch_picard_step_length_min = -1.0;
    solver_options->nl_switch_picard_after_amr = 0;
    solver_options->nl_switch_picard_init = 0;
    break;

  case SL_NL_SOLVER_PICARD_NEWTON:
    /* keep same parameters */
    break;

  default: /* unknown nonlinear solver */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* run nonlinear solver */
  slabs_nonlinear_solver_picard_newton (state, nl_stokes, physics_options,
                                        discr_options, solver_options,
                                        bin_filepath, vtk_filepath);

  /* restore solver options */
  switch (solver_options->nl_solver_type) {
  case SL_NL_SOLVER_PICARD:
    solver_options->nl_switch_picard_init = switch_picard_init;
    solver_options->nl_switch_picard_maxiter = switch_picard_maxiter;
    solver_options->nl_switch_picard_rtol = switch_picard_rtol;
    break;

  case SL_NL_SOLVER_NEWTON:
    solver_options->nl_switch_picard_step_length_min =
      switch_picard_step_length_min;
    solver_options->nl_switch_picard_after_amr = switch_picard_after_amr;
    solver_options->nl_switch_picard_init = switch_picard_init;
    break;

  case SL_NL_SOLVER_PICARD_NEWTON:
    /* nothing to restore */
    break;

  default: /* unknown nonlinear solver */
    YMIR_ABORT_NOT_REACHED ();
  }
}

/******************************************************************************
 * Main program
 *****************************************************************************/

/**
 * Solves a Stokes problem.
 */
static void
slabs_solve_stokes (slabs_lin_stokes_problem_t *lin_stokes,
                    slabs_nl_stokes_problem_t **nl_stokes,
                    p8est_t *p8est,
                    ymir_mesh_t **mesh,
                    ymir_pressure_elem_t **press_elem,
                    slabs_stokes_state_t *state,
                    slabs_physics_options_t *physics_options,
                    slabs_discr_options_t *discr_options,
                    slabs_nl_solver_options_t *solver_options,
                    const char *workload_out_path,
                    const char *bin_nl_filepath,
                    const char *vtk_nl_filepath)
{
  const char         *this_fn_name = "slabs_solve_stokes";
  const slabs_nl_solver_type_t  nl_solver_type = solver_options->nl_solver_type;
  MPI_Comm            mpicomm = p8est->mpicomm;
  char                path[BUFSIZ];

  YMIR_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

#if defined(__bgq__)

#if defined(__HAVE_MONEQ)
  // MonEQ start
  MonEQ_StartPowerTag("YMIR_LOOP_Block");
  if(p8est->mpirank == 0)
    printf("MONEQ LOOP START\n");
#endif

#if defined(__HAVE_HPM)
  HPM_Start("YMIR_LOOP_Block");
  if(p8est->mpirank == 0)
    printf("HPM LOOP START\n");
#endif

//  trace_start();
  summary_start();
#endif

  /* start performance counters */
  ymir_perf_counter_start_barrier (
      &slabs_perf_counter[SLABS_PERF_COUNTER_SOLVE_STOKES], mpicomm);

  /* run solver */
  if (nl_solver_type == SL_NL_SOLVER_NONE) { /* if linear solve */
    const int           krylov_maxiter = solver_options->krylov_maxiter;
    const double        krylov_atol = solver_options->krylov_atol;
    const double        krylov_rtol = solver_options->krylov_rtol;
    const int           krylov_gmres_num_vecs =
                          solver_options->krylov_gmres_num_vecs;
    const int           lin_solve_stress_block_only =
                          solver_options->lin_solve_stress_block_only;
    const int           lin_solve_press_bbt_only =
                          solver_options->lin_solve_press_bbt_only;
    const int           lin_solve_cnode_bbt_only =
                          solver_options->lin_solve_cnode_bbt_only;
#if !defined(_MINIMIZE_MEMORY)
    ymir_vec_t         *residual_up = ymir_stokes_vec_new (*mesh, *press_elem);
#endif
    int                 conv_reason;
    int                 n_iter;
    double              norm_res = 1.0;
    double              norm_res_init = -1.0;
    double              conv_rate;

    /*
     * Linear Solver
     */

    YMIR_ASSERT (lin_stokes != NULL);

    /* invoke linear solver */
    if (!lin_solve_stress_block_only && !lin_solve_press_bbt_only &&
        !lin_solve_cnode_bbt_only) {
      /*
       * Solve Stokes System
       */

      /* compute residual with zero inital guess */
#if !defined(_MINIMIZE_MEMORY)
      ymir_vec_set_zero (state->vel_press_vec);
      norm_res_init = slabs_norm_compute_residual (
          NULL, NULL, NULL, NULL, lin_stokes->rhs_u_point,
          lin_stokes->stokes_op, lin_stokes->stokes_pc,
          solver_options->krylov_type, SL_NORM_VEC_L2, NULL);
#endif

      /* run solver for Stokes */
      conv_reason = ymir_stokes_pc_solve (lin_stokes->rhs,
                                          state->vel_press_vec,
                                          lin_stokes->stokes_pc,
                                          0 /* zero inital guess */,
                                          krylov_rtol, krylov_atol,
                                          krylov_maxiter,
                                          krylov_gmres_num_vecs /* unused */,
                                          &n_iter);

      /* compute residual at solution */
#if !defined(_MINIMIZE_MEMORY)
      norm_res = slabs_norm_compute_residual (
          residual_up, NULL, NULL, state->vel_press_vec,
          lin_stokes->rhs_u_point,
          lin_stokes->stokes_op, lin_stokes->stokes_pc,
          solver_options->krylov_type, SL_NORM_VEC_L2, NULL);
#endif
    }
    else if (lin_solve_stress_block_only) {
      const int           N = ymir_n ((*mesh)->ma->N);
      ymir_vec_t         *u, *rhs_u;

      /*
       * Solve Stress Block
       */

      /* compute residual with zero inital guess */
#if !defined(_MINIMIZE_MEMORY)
      ymir_vec_set_zero (state->vel_press_vec);
      slabs_norm_compute_residual (
          NULL, &norm_res_init, NULL, NULL, lin_stokes->rhs_u_point,
          lin_stokes->stokes_op, lin_stokes->stokes_pc,
          solver_options->krylov_type, SL_NORM_VEC_L2, NULL);
#endif

      /* run solver for viscous stress block */
      if (N == 1) {
        u = ymir_cvec_new (*mesh, 3);
        rhs_u = ymir_cvec_new (*mesh, 3);
        ymir_stokes_vec_get_velocity (state->vel_press_vec, u, *press_elem);
        ymir_stokes_vec_get_velocity (lin_stokes->rhs, rhs_u, *press_elem);
      }
      else {
        slabs_stokes_vec_get_components_view (&u, NULL, state->vel_press_vec);
        slabs_stokes_vec_get_components_view (&rhs_u, NULL, lin_stokes->rhs);
      }
      conv_reason = ymir_stress_pc_solve_ext (rhs_u, u,
                                              lin_stokes->stokes_pc->stress_pc,
                                              0 /* zero initial guess */,
                                              krylov_rtol, krylov_atol,
                                              krylov_maxiter, &n_iter);
      if (N == 1) {
        ymir_stokes_vec_set_velocity (u, state->vel_press_vec, *press_elem);
      }
      ymir_vec_destroy (u);
      ymir_vec_destroy (rhs_u);

      /* compute residual at solution */
#if !defined(_MINIMIZE_MEMORY)
      slabs_norm_compute_residual (
          residual_up, &norm_res, NULL, state->vel_press_vec,
          lin_stokes->rhs_u_point,
          lin_stokes->stokes_op, lin_stokes->stokes_pc,
          solver_options->krylov_type, SL_NORM_VEC_L2, NULL);
#endif
    }
    else if (lin_solve_press_bbt_only) {
      const int           N = ymir_n ((*mesh)->ma->N);
      ymir_vec_t         *p;
      ymir_vec_t         *rhs_p = ymir_pressure_vec_new (*mesh, *press_elem);
      ymir_vec_t         *res_p = ymir_pressure_vec_new (*mesh, *press_elem);
      ymir_bbt_t         *bbt;

      /*
       * Solve pressure BB^T Block
       */

      YMIR_ASSERT (lin_stokes->stokes_pc->bfbt != NULL);
      YMIR_ASSERT (lin_stokes->stokes_pc->bfbt->right_bbt != NULL);
      bbt = lin_stokes->stokes_pc->bfbt->right_bbt;

      /* apply divergence operator to velocity right-hand side */
      ymir_pressure_vec_weak_grad (rhs_p, lin_stokes->rhs_u_point, *press_elem,
                                   lin_stokes->vel_dir, 1);
      ymir_bbt_project_out_nullspace (rhs_p, bbt,
                                      YMIR_STOKES_OP_PRESS_PROJECT_SYMM, 1);

      /* compute residual with zero inital guess */
      ymir_vec_set_zero (state->vel_press_vec);
      norm_res_init = ymir_vec_norm (rhs_p);

      /* solve pressure BB^T */
      if (N == 1) {
        p = ymir_pressure_vec_new (*mesh, *press_elem);
        ymir_stokes_vec_get_velocity (state->vel_press_vec, p, *press_elem);
      }
      else {
        slabs_stokes_vec_get_components_view (NULL, &p, state->vel_press_vec);
      }
#ifdef YMIR_PETSC
      {
        PetscErrorCode      ierr;

        ierr = KSPSetInitialGuessNonzero (bbt->ksp, PETSC_FALSE);
        YMIR_CHECK_PETSC (ierr);
      }
#endif
      conv_reason = ymir_bbt_solve (p, rhs_p, bbt, 0, &n_iter);

      /* hack to avoid SEGV when destroying stress PC (TODO fix) */
      {
        ymir_vec_t         *tmp = ymir_stokes_vec_new (*mesh, *press_elem);

        ymir_stokes_pc_solve (lin_stokes->rhs, tmp, lin_stokes->stokes_pc,
                              0, 0.1, 0.0, 1, 1, NULL);
        ymir_vec_destroy (tmp);
      }

      /* compute residual at solution */
      ymir_bbt_apply (p, res_p, bbt);
      ymir_vec_add (-1.0, rhs_p, res_p);
      ymir_vec_scale (-1.0, res_p);
      norm_res = ymir_vec_norm (res_p);

      /* set output vectors */
#if !defined(_MINIMIZE_MEMORY)
      ymir_vec_set_zero (residual_up);
      ymir_stokes_vec_set_pressure (residual_up, res_p, *press_elem);
      ymir_stokes_vec_set_pressure (state->vel_press_vec, p, *press_elem);
#endif

      /* destroy */
      ymir_vec_destroy (p);
      ymir_vec_destroy (rhs_p);
      ymir_vec_destroy (res_p);
    }
    else if (lin_solve_cnode_bbt_only) {
      ymir_vec_t         *sol = ymir_cvec_new (*mesh, 1);
      ymir_vec_t         *rhs = ymir_cvec_new (*mesh, 1);
      ymir_vec_t         *res = ymir_cvec_new (*mesh, 1);
      ymir_vec_t         *rhs_p = ymir_pressure_vec_new (*mesh, *press_elem);
      ymir_vec_t         *lump_mass = ymir_pressure_vec_new (*mesh,
                                                             *press_elem);
      ymir_bbt_t         *bbt;
      ymir_stiff_pc_t    *stiff_pc;

      /*
       * Solve Continuous Approximation of BB^T Block
       */

      YMIR_ASSERT (lin_stokes->stokes_pc->bfbt != NULL);
      YMIR_ASSERT (lin_stokes->stokes_pc->bfbt->right_bbt != NULL);
      YMIR_ASSERT (lin_stokes->stokes_pc->bfbt->right_bbt->stiff_pc != NULL);
      bbt = lin_stokes->stokes_pc->bfbt->right_bbt;
      stiff_pc = lin_stokes->stokes_pc->bfbt->right_bbt->stiff_pc;

      /* apply divergence operator to velocity right-hand side */
      ymir_pressure_vec_weak_grad (rhs_p, lin_stokes->rhs_u_point, *press_elem,
                                   lin_stokes->vel_dir, 1);
      ymir_bbt_project_out_nullspace (rhs_p, bbt,
                                      YMIR_STOKES_OP_PRESS_PROJECT_SYMM, 1);
      ymir_pressure_vec_lump_mass (lump_mass, *press_elem);
      ymir_vec_divide_in (lump_mass, rhs_p);
      ymir_mass_apply (rhs_p, rhs);
      ymir_stiff_op_project_out_nullspace (rhs, stiff_pc->stiff_op,
                                           YMIR_STIFF_OP_MEAN_PROJECT_SYMM, 1);

      /* compute residual with zero inital guess */
      ymir_vec_set_zero (sol);
      norm_res_init = ymir_vec_norm (rhs);

      /* solve pressure BB^T */
      conv_reason = ymir_stiff_pc_solve (sol, rhs, stiff_pc,
                                         0 /* zero initial guess */, &n_iter);

      /* hack to avoid SEGV when destroying stress PC (TODO fix) */
      {
        ymir_vec_t         *tmp = ymir_stokes_vec_new (*mesh, *press_elem);

        ymir_stokes_pc_solve (lin_stokes->rhs, tmp, lin_stokes->stokes_pc,
                              0, 0.1, 0.0, 1, 1, NULL);
        ymir_vec_destroy (tmp);
      }

      /* compute residual at solution */
      ymir_stiff_op_apply (sol, res, stiff_pc->stiff_op);
      ymir_vec_add (-1.0, rhs, res);
      ymir_vec_scale (-1.0, res);
      norm_res = ymir_vec_norm (res);

      /* set output vectors */
#if !defined(_MINIMIZE_MEMORY)
      ymir_vec_set_zero (residual_up);
      ymir_vec_set_zero (state->vel_press_vec);
      ymir_cvec_set_comp (residual_up, res->cvec, 0, YMIR_SET);
      ymir_cvec_set_comp (state->vel_press_vec, sol->cvec, 0, YMIR_SET);
#endif

      /* destroy */
      ymir_vec_destroy (sol);
      ymir_vec_destroy (rhs);
      ymir_vec_destroy (res);
      ymir_vec_destroy (rhs_p);
      ymir_vec_destroy (lump_mass);
    }
    else {
      YMIR_ABORT_NOT_REACHED ();
    }

    /* calculate rate of convergence */
    conv_rate = exp (log (norm_res / norm_res_init) / ((double) n_iter));

    YMIR_GLOBAL_PRODUCTIONF (
        "%s: Linear solver converged reason %d, "
        "%d iterations for reduction %.3e, convergence rate %.3e\n",
        this_fn_name, conv_reason, n_iter, norm_res/norm_res_init, conv_rate);

    /* project out nullspaces of linear solution */
    slabs_linear_stokes_problem_project_out_nullspace (
        state->vel_press_vec, lin_stokes->stokes_op, physics_options, 0);

    /* print all major physics statistics */
    if (solver_options->log_physics_stats) {
      ymir_vec_t         *viscosity_nondim =
        ymir_vec_clone (lin_stokes->stokes_op->stress_op->viscosity);

      ymir_vec_scale (0.5, viscosity_nondim);
      slabs_physics_stats_print_all (state->vel_press_vec, viscosity_nondim,
                                     NULL, lin_stokes->rhs_u_point,
                                     lin_stokes->stokes_op, physics_options);
      ymir_vec_destroy (viscosity_nondim);
    }

    /* destroy */
#if !defined(_MINIMIZE_MEMORY)
    ymir_vec_destroy (residual_up);
#endif
  }
  else {
    /*
     * Nonlinear Solver
     */

    YMIR_ASSERT (*nl_stokes != NULL);

    /* invoke nonlinear solver */
    slabs_nonlinear_solver_solve_state (state, nl_stokes, physics_options,
                                        discr_options, solver_options,
                                        bin_nl_filepath, vtk_nl_filepath);

    /* update ymir mesh and pressure element */
    *mesh = (*nl_stokes)->mesh;
    *press_elem = (*nl_stokes)->press_elem;
  }

  /* stop performance counters */
  ymir_perf_counter_stop_add (
      &slabs_perf_counter[SLABS_PERF_COUNTER_SOLVE_STOKES]);

#if defined(__bgq__)
  summary_stop();
//  trace_stop();

#if defined(__HAVE_HPM)
  HPM_Stop("YMIR_LOOP_Block");
  if(p8est->mpirank == 0)
    printf("HPM LOOP STOP\n");
#endif

#if defined(__HAVE_MONEQ)
  // Stop recording power
  MonEQ_EndPowerTag("YMIR_LOOP_Block");
  if(p8est->mpirank == 0)
    printf("MONEQ LOOP STOP\n");

  // MonEQ Finalize
  int moneq_status = MonEQ_Finalize();
  if (moneq_status) printf("Error finalizing MonEQ\n");
#endif

#endif

  /* print mesh statistics */
  ymir_monitor_print_global_mesh_stats (*mesh, *press_elem);

  /* print memory usage */
  ymir_monitor_print_global_mem_usage (mpicomm);

  /* write global workload to a file */
  if (workload_out_path != NULL) {
    snprintf (path, BUFSIZ, "%s_solve_stokes", workload_out_path);
    ymir_monitor_dump_global_mesh_mem_stats (*mesh, *press_elem, path);
    ymir_monitor_dump_global_mem_usage_petsc (mpicomm, path);
  }

  YMIR_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Initializes performance counters.
 */
static void
slabs_perf_counter_init (const int active)
{
  slabs_io_perf_counter_init (active);

  ymir_perf_counter_init_all (slabs_perf_counter, slabs_perf_counter_name,
                              SLABS_PERF_COUNTER_N, active);
}

/**
 * Gathers statistics of performance counters.
 */
static void
slabs_perf_counter_gather (MPI_Comm mpicomm,
                           const int print_wtime,
                           const int print_n_calls,
                           const int print_flops)
{
  slabs_perf_n_stats = ymir_perf_counter_gather_stats (
      slabs_perf_counter, SLABS_PERF_COUNTER_N,
      slabs_perf_stats, slabs_perf_stats_name,
      mpicomm, print_wtime, print_n_calls, print_flops);
}

/**
 * Prints statistics of performance counters.
 */
static void
slabs_perf_counter_print ()
{
  ymir_perf_counter_print_stats (slabs_perf_stats, slabs_perf_n_stats,
                                 "Slabs example main");
}

/**
 * Runs the program.
 */
int
main (int argc, char **argv)
{
  const char         *this_fn_name = "slabs:main";
  /* MPI */
  MPI_Comm            mpicomm = MPI_COMM_WORLD;
  int                 mpisize, mpirank, ompsize;
  int                 mpiret;
  /* options */
  ymir_options_t     *opt;
  rhea_domain_options_t         domain_options;
  rhea_temperature_options_t    temp_options;
  rhea_viscosity_options_t      viscosity_options;
  rhea_discretization_options_t discr_options;
  rhea_newton_options_t         newton_options;
  slabs_physics_options_t   slabs_physics_options;
  slabs_discr_options_t     slabs_discr_options;
  slabs_nl_solver_options_t slabs_solver_options;
  /* options local to this function */
  double              visc_override_const;
  double              weak_override_factor;
  int                 production_run;
  int                 monitor_performance;
  char               *workload_out_path;
  char               *bin_out_filepath;
  int                 bin_out_input, bin_out_nl_iter, bin_out_solution;
  char               *vtk_out_path;
  int                 vtk_input, vtk_p8est, vtk_nl_iter;
  int                 vtk_input_state;
  int                 exit_after_mesh_setup;
  int                 test_strain_rate, test_stress, test_stiffness;
  int                 test_norms, test_interpolation;
  /* mesh */
  p8est_t            *p8est;
  ymir_mesh_t        *mesh;
  ymir_pressure_elem_t  *press_elem;
  /* Stokes problem */
  slabs_stokes_state_t  *state;
  slabs_discr_enforce_refinement_data_t     *enforce_refinement_data;
  slabs_physics_coarsen_stokes_coeff_data_t *coarsen_coeff_data;
  slabs_lin_stokes_problem_t  *lin_stokes = NULL;
  slabs_nl_stokes_problem_t   *nl_stokes = NULL;
  /* other */
  char                path[BUFSIZ];

  /*
   * Initialize Libraries
   */

  /* initialize rhea and sub-packages */
  rhea_initialize (argc, argv, mpicomm);

  /* get parallel environment */
  mpiret = MPI_Comm_size (mpicomm, &mpisize); YMIR_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &mpirank); YMIR_CHECK_MPI (mpiret);

#ifdef YMIR_ENABLE_OPENMP
  #pragma omp parallel
  {
    #pragma omp single
    {
      ompsize = omp_get_num_threads ();
    }
  }
#else
  ompsize = 1;
#endif

  /*
   * Define & Parse Options
   */

  opt = ymir_options_global_new (argv[0] /* program path */);

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  /* basic options */
  YMIR_OPTIONS_CALLBACK, "help", 'h', 0 /* no callback fn args */,
    ymir_options_print_usage_and_exit_fn, NULL /* no arg usage */,
    "Print usage and exit",
  YMIR_OPTIONS_INIFILE, "options-file", 'f',
    ".ini file with option values",

  /* override physics options */
  YMIR_OPTIONS_D, "viscosity-override-const", '\0',
    &visc_override_const, -1.0,
    "Override viscosity to be constant after mesh was generated",
  YMIR_OPTIONS_D, "weakzone-override-factor", '\0',
    &weak_override_factor, -1.0,
    "Override weak zone factor after mesh was generated",

  /* performance & monitoring options */
  YMIR_OPTIONS_B, "production-run", '\0',
                  &production_run, 0,
    "Execute as a production run (to reduce some overhead and checks)",

  YMIR_OPTIONS_B, "monitor-performance", '\0',
                  &monitor_performance, 0,
    "Measure and print performance statistics, e.g., runtime or flops",

  YMIR_OPTIONS_S, "workload-out-path", '\0',
                  &workload_out_path, NULL,
    "File path for output of global workload for each processor",

  /* bin output options */
  YMIR_OPTIONS_S, "bin-out", '\0', &bin_out_filepath, NULL,
    "File path for binary output, activates output of solution",
  YMIR_OPTIONS_I, "bin-out-input", '\0', &bin_out_input, 0,
    "Write binary output of input data",
  YMIR_OPTIONS_I, "bin-out-nonlinear-iter", '\0', &bin_out_nl_iter, 0,
    "Write binary output of each nonlinear iteration",
  YMIR_OPTIONS_I, "bin-out-solution", '\0', &bin_out_solution, 0,
    "Write binary output of data at solution",

  /* vtk output options */
  YMIR_OPTIONS_S, "vtk-out-path", '\0', &vtk_out_path, NULL,
    "File path for vtk files, activates output of solution",
  YMIR_OPTIONS_I, "vtk-out-input-state", '\0', &vtk_input_state, 0,
    "Write vtk output of input Stokes state",
  YMIR_OPTIONS_I, "vtk-out-input", '\0', &vtk_input, 0,
    "Write vtk output of input data",
  YMIR_OPTIONS_I, "vtk-out-nonlinear-iter", '\0', &vtk_nl_iter, 0,
    "Write vtk output of each nonlinear iteration",
  YMIR_OPTIONS_I, "vtk-out-p8est", '\0', &vtk_p8est, 0,
    "Write vtk output of p8est mesh",

  YMIR_OPTIONS_I, "exit-after-mesh-setup", '\0', &exit_after_mesh_setup, 0,
    "Exit program after mesh setup and vtk output",

  /* testing options */
  YMIR_OPTIONS_I, "test-strain-rate", '\0', &test_strain_rate, 0,
    "Run test: strain rate tensor and 2nd invariant computation",
  YMIR_OPTIONS_I, "test-stress", '\0', &test_stress, 0,
    "Run test: viscous stress operator",
  YMIR_OPTIONS_I, "test-stiffness", '\0', &test_stiffness, 0,
    "Run test: stiffness operator with and without coefficient",
  YMIR_OPTIONS_I, "test-norms", '\0', &test_norms, 0,
    "Run test: mesh independence of norms",
  YMIR_OPTIONS_I, "test-interpolation", '\0', &test_interpolation, 0,
    "Run test: interpolation of fields from coarse to fine mesh",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add sub-options */
  slabs_setup_add_suboptions (opt);
  rhea_add_options_all (opt);
  ymir_options_add_suboptions_solver_stokes (opt);

  /* parse options */
  {
    int                 optret;

    optret = ymir_options_parse (SC_LP_INFO, opt, argc, argv);
    if (optret < 0) { /* if parsing was not successful */
      ymir_options_print_usage (SC_LP_INFO, opt, NULL /* args usage */);
      YMIR_GLOBAL_INFO ("Option parsing failed\n");
      exit (0);
    }
  }

  /*
   * Initialize Main Program
   */

  YMIR_GLOBAL_PRODUCTIONF ("Into %s (production %i)\n",
      this_fn_name, production_run);
  YMIR_GLOBAL_PRODUCTIONF (
      "Parallel environment: MPI size %i, OpenMP size %i\n",
      mpisize, ompsize);
  ymir_set_up (argc, argv, mpicomm, production_run);

  /* print memory usage */
  ymir_monitor_print_global_mem_usage (mpicomm);

  /* write memory usage to file */
  if (workload_out_path != NULL) {
    snprintf (path, BUFSIZ, "%s_main_begin", workload_out_path);
    ymir_monitor_dump_global_mem_usage_petsc (mpicomm, path);
  }

  /* initialize performance counters */
  slabs_perf_counter_init (monitor_performance);

  /* start performance counters */
  ymir_perf_counter_start_barrier (
      &slabs_perf_counter[SLABS_PERF_COUNTER_TOTAL], mpicomm);

  /* print & process options */
  ymir_options_print_summary (SC_LP_INFO, opt);
  rhea_process_options_all (&domain_options, &temp_options,
                            &viscosity_options, &discr_options,
                            &newton_options);
  slabs_setup_process_options (&slabs_physics_options, &slabs_discr_options,
                               &slabs_solver_options);
  slabs_physics_options.domain_options = &domain_options;
  slabs_physics_options.temp_options = &temp_options;
  slabs_physics_options.viscosity_options = &viscosity_options;
  slabs_discr_options.inspect_p4est = monitor_performance;

  /*
   * Setup Mesh
   */

  /* create mesh and Stokes state */
  slabs_setup_mesh (&p8est, &mesh, &press_elem, &state,
                    &enforce_refinement_data, &coarsen_coeff_data,
                    mpicomm, &slabs_physics_options, &slabs_discr_options,
                    &slabs_solver_options,
                    &slabs_perf_counter[SLABS_PERF_COUNTER_SETUP_MESH],
                    workload_out_path);

  /* binary output of Stokes state */
  if (bin_out_filepath != NULL && bin_out_input) {
    snprintf (path, BUFSIZ, "%s_input", bin_out_filepath);
    slabs_stokes_state_save (state, path);
  }

  /* vtk output of Stokes state and viscosity */
  if (vtk_out_path != NULL && vtk_input_state) {
    snprintf (path, BUFSIZ, "%s_input", vtk_out_path);
    slabs_vtk_write_state (path, state, NULL, 1, 1, mesh, press_elem,
                           &slabs_physics_options);
  }

  /* vtk output of p4est mesh */
  if (vtk_out_path != NULL && vtk_p8est) {
    snprintf (path, BUFSIZ, "%s_p4est", vtk_out_path);
    slabs_vtk_write_p8est (path, p8est, slabs_physics_options.domain_shape);
  }

  /* exit after mesh generation */
  if (exit_after_mesh_setup) {
    /* destroy Stokes problem, Stokes state, and mesh */
    slabs_clear (lin_stokes, nl_stokes, p8est, mesh, press_elem, state,
                 enforce_refinement_data, coarsen_coeff_data,
                 &slabs_physics_options, &slabs_discr_options);

    /* destroy options */
    ymir_options_global_destroy ();

    /* print that this function is ending */
    YMIR_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

    /* finalize rhea */
    rhea_finalize ();

    return 0;
  }

  /*
   * Run Tests
   */

  {
    const slabs_domain_shape_t  domain_shape =
                                  slabs_physics_options.domain_shape;
    const int8_t        minlevel = slabs_discr_options.minlevel;
    const int8_t        maxlevel = slabs_discr_options.maxlevel;
    char               *refine = slabs_discr_options.refine;
    const int           N = slabs_discr_options.order;
    const double        nl_norm_Hminus1_mass_scaling =
                          slabs_solver_options.norm_Hminus1_mass_scaling;

    /* test strain rate tensor and 2nd invariant computation */
    if (test_strain_rate) {
      if (domain_shape == SL_DOMAIN_CUBE) {
        /* test computation of 2nd invariant of the strain rate */
        slabs_base_test_second_invariant (mesh);

        /* test computation of the strain rate tensor */
        slabs_base_test_strain_rate_tensor (mesh);
      }
    }

    /* test viscous stress operator */
    if (test_stress) {
      if (domain_shape == SL_DOMAIN_CUBE) {
        slabs_test_stress (mpicomm, domain_shape, minlevel, maxlevel,
                           refine, N);
      }
    }

    /* test stiffness operator with and without coefficient */
    if (test_stiffness) {
      if (domain_shape == SL_DOMAIN_CUBE) {
        slabs_test_stiffness (mpicomm, domain_shape, minlevel, maxlevel,
                              refine, N);
      }
    }

    /* test norms */
    if (test_norms) {
      if (domain_shape == SL_DOMAIN_CUBE) {
        ymir_Hminus1_norm_op_t  *norm_op;

        /* test Frobenius-norm of the strain rate tensor */
        slabs_norm_test_frobenius_of_strain_rate (mesh);

        /* test H^-1 norm operator */
        norm_op = ymir_Hminus1_norm_op_new (mesh, nl_norm_Hminus1_mass_scaling);
        slabs_norm_Hminus1_stiff_test (norm_op, mesh, domain_shape);
        ymir_Hminus1_norm_op_destroy (norm_op);

        /* test norms for mesh-independence */
        slabs_norm_test_mesh_independence (mpicomm, domain_shape,
                                           minlevel, maxlevel, refine, N,
                                           nl_norm_Hminus1_mass_scaling, 3);

        /* test Stokes norms for mesh-independence */
        slabs_test_mesh_independence_stokes_norm (mpicomm, domain_shape,
                                                  minlevel, maxlevel, refine, N,
                                                  nl_norm_Hminus1_mass_scaling,
                                                  3);
      }
      else if (domain_shape == SL_DOMAIN_BRICK) {
        /* test Stokes norms for mesh-independence */
        slabs_test_mesh_independence_stokes_norm (mpicomm, domain_shape,
                                                  minlevel, maxlevel, refine, N,
                                                  nl_norm_Hminus1_mass_scaling,
                                                  3);
      }
    }

    /* test interpolation */
    if (test_norms) {
      if (domain_shape == SL_DOMAIN_CUBE) {
        slabs_test_interpolation (mpicomm, domain_shape, minlevel, maxlevel,
                                  refine, N);
      }
    }
  }

  /*
   * Setup Stokes Problem
   */

  /* override viscosity to be constant = 1 (for testing purposes) */
  if (0.0 < visc_override_const) {
    YMIR_GLOBAL_INFOF ("%s: Warning: Override viscosity to be const %g!\n",
                       this_fn_name, visc_override_const);
    slabs_physics_options.viscosity_type = SL_VISCOSITY_CONST;
    slabs_physics_options.viscosity_scaling = visc_override_const;
  }

  /* override weak zone factor (for testing purposes) */
  if (0.0 < weak_override_factor) {
    YMIR_GLOBAL_INFOF ("%s: Warning: Override weak zone factor to be %1.3e!\n",
                       this_fn_name, weak_override_factor);
    slabs_physics_options.weakzone_import_weak_factor = weak_override_factor;
    slabs_physics_options.weakzone_2plates_subdu_weak_factor =
      weak_override_factor;
    slabs_physics_options.weakzone_2plates_ridge_weak_factor =
      weak_override_factor;
    slabs_physics_compute_weakzone (state->weak_vec, &slabs_physics_options);
  }

  /* setup Stokes problem */
  slabs_setup_stokes (&lin_stokes, &nl_stokes, p8est, mesh, press_elem, state,
                      &slabs_physics_options, &slabs_solver_options,
                      &slabs_perf_counter[SLABS_PERF_COUNTER_SETUP_STOKES],
                      workload_out_path);

  /* vtk output of input data */
  if (vtk_out_path != NULL && vtk_input) {
    snprintf (path, BUFSIZ, "%s_input", vtk_out_path);

    /* write vtk file */
    if (slabs_solver_options.nl_solver_type == SL_NL_SOLVER_NONE) {
      /* if linear solve */
      slabs_vtk_write_input_lin_stokes (path, state, lin_stokes,
                                        &slabs_physics_options,
                                        &slabs_discr_options);
    }
    else { /* if nonlinear solve */
      slabs_vtk_write_input_nl_stokes (path, state, nl_stokes,
                                       &slabs_physics_options,
                                       &slabs_discr_options);
    }
  }

  //###DEV###
  //mpiret = MPI_Barrier (mpicomm); YMIR_CHECK_MPI (mpiret);
  //YMIR_GLOBAL_LERROR ("###DEV### stop :-[\n");
  //return 1;

  /*
   * Solve Stokes Problem
   */
  {
    char               *bin_nl_filepath = NULL;
    char               *vtk_nl_filepath = NULL;

    /* set output filenames for nonlinear iterations */
    if (slabs_solver_options.nl_solver_type != SL_NL_SOLVER_NONE) {
      if (bin_out_filepath != NULL && bin_out_nl_iter) {
        bin_nl_filepath = bin_out_filepath;
      }
      if (vtk_out_path != NULL && vtk_nl_iter) {
        vtk_nl_filepath = vtk_out_path;
      }
    }

    /* run solver */
    slabs_solve_stokes (lin_stokes, &nl_stokes, p8est, &mesh, &press_elem,
                        state, &slabs_physics_options, &slabs_discr_options,
                        &slabs_solver_options, workload_out_path,
                        bin_nl_filepath, vtk_nl_filepath);
  }

  /* stop performance counters */
  ymir_perf_counter_stop_add_barrier (
      &slabs_perf_counter[SLABS_PERF_COUNTER_TOTAL], mpicomm);

  /*
   * Output of Results
   */

  /* binary output of Stokes state at solution */
  if (bin_out_filepath != NULL && bin_out_solution) {
    snprintf (path, BUFSIZ, "%s_solution", bin_out_filepath);
    slabs_stokes_state_save (state, path);
  }

  /* vtk output of solution */
  if (vtk_out_path != NULL) {
    snprintf (path, BUFSIZ, "%s_solution", vtk_out_path);

    if (slabs_solver_options.nl_solver_type == SL_NL_SOLVER_NONE) {
      /* if linear solve */
      slabs_vtk_write_solution_lin_stokes (path, state, lin_stokes,
                                           slabs_solver_options.krylov_type,
                                           &slabs_physics_options,
                                           &slabs_discr_options);
    }
    else {
      /* if nonlinear solve */
      slabs_vtk_write_solution_nl_stokes (path, state, nl_stokes,
                                          &slabs_physics_options,
                                          &slabs_discr_options);
    }
  }

  //###DEV### papi test
#if 0
  if (slabs_solver_options.nl_solver_type == SL_NL_SOLVER_NONE) {
    ymir_vec_t         *a = ymir_dvec_new (mesh, 1, YMIR_GLL_NODE);
    ymir_vec_t         *b = ymir_dvec_new (mesh, 1, YMIR_GLL_NODE);
    ymir_vec_t         *x = ymir_cvec_new (mesh, 3);
    ymir_vec_t         *y = ymir_cvec_new (mesh, 3);
    ymir_perf_counter_t counter_ip, counter_stress;
    sc_statinfo_t       stats_ip[2], stats_stress[2];

    ymir_vec_set_value (a, 1.0);
    ymir_vec_set_value (b, 2.0);
    ymir_vec_set_value (x, 3.0);

    mpiret = MPI_Barrier (mpicomm); YMIR_CHECK_MPI (mpiret);
    ymir_perf_counter_init (&counter_ip, NULL, 1);
    ymir_perf_counter_init (&counter_stress, NULL, 1);

    mpiret = MPI_Barrier (mpicomm); YMIR_CHECK_MPI (mpiret);
    ymir_perf_counter_start (&counter_ip);
    ymir_dvec_innerprod_local (a, b);
    ymir_perf_counter_stop_add (&counter_ip);

    mpiret = MPI_Barrier (mpicomm); YMIR_CHECK_MPI (mpiret);
    ymir_perf_counter_start (&counter_stress);
    ymir_stress_op_apply_dirty (x, y, lin_stokes->stokes_op->stress_op);
    ymir_perf_counter_stop_add (&counter_stress);

    ymir_vec_destroy (a);
    ymir_vec_destroy (b);
    ymir_vec_destroy (x);
    ymir_vec_destroy (y);

    sc_stats_set1 (&stats_ip[0], counter_ip.wtime_cumul,
                   "###DEV### inner product Wtime");
    sc_stats_set1 (&stats_ip[1], counter_ip.flops_cumul,
                   "###DEV### inner product Flops");

    sc_stats_set1 (&stats_stress[0], counter_stress.wtime_cumul,
                   "###DEV### stress apply Wtime");
    sc_stats_set1 (&stats_stress[1], counter_stress.flops_cumul,
                   "###DEV### stress apply Flops");

    sc_stats_compute (mpicomm, 2, stats_ip);
    sc_stats_print (-1, SC_LP_STATISTICS, 2, stats_ip, 1, 1);

    sc_stats_compute (mpicomm, 2, stats_stress);
    sc_stats_print (-1, SC_LP_STATISTICS, 2, stats_stress, 1, 1);
  }
#endif
  //###DEV### end

  /*
   * Finalize
   */

  /* destroy Stokes problem, Stokes state, and mesh */
  slabs_clear (lin_stokes, nl_stokes, p8est, mesh, press_elem, state,
               enforce_refinement_data, coarsen_coeff_data,
               &slabs_physics_options, &slabs_discr_options);

  /* print memory usage */
  ymir_monitor_print_global_mem_usage (mpicomm);

  /* write memory usage to file */
  if (workload_out_path != NULL) {
    snprintf (path, BUFSIZ, "%s_main_end", workload_out_path);
    ymir_monitor_dump_global_mem_usage_petsc (mpicomm, path);
  }

  /* print performance statistics */
  ymir_gmg_hierarchy_mesh_perf_counter_print ();   /* GMG mesh */
  ymir_stiff_op_perf_counter_print ();             /* Stiffness Op */
  ymir_stiff_pc_perf_counter_print ();             /* Stiffness PC */
  ymir_gmg_hierarchy_stiff_perf_counter_print ();  /* GMG stiffness */
  ymir_stress_op_perf_counter_print ();            /* Stress Op */
  ymir_stress_pc_perf_counter_print ();            /* Stress PC */
  ymir_gmg_hierarchy_stress_perf_counter_print (); /* GMG stress */
  ymir_pressure_vec_perf_counter_print ();         /* B^T or B */
  ymir_bbt_perf_counter_print ();                  /* BB^T */
//ymir_gmg_hierarchy_bbt_perf_counter_print ();    /* GMG BB^T */
  ymir_bfbt_perf_counter_print ();                 /* BFBT */
  ymir_stokes_op_perf_counter_print ();            /* Stokes Op */
  ymir_stokes_pc_perf_counter_print ();            /* Stokes PC */

  /* print slabs performance statistics */
  if (monitor_performance) {
    const int           print_wtime = 1;
    const int           print_n_calls = 0;
    const int           print_flops = 1;

    /* add counters from all matvecs */
    ymir_stress_op_perf_counter_add_matvecs (
        &slabs_perf_counter[SLABS_PERF_COUNTER_MATVECS]);
    ymir_stiff_op_perf_counter_add_matvecs (
        &slabs_perf_counter[SLABS_PERF_COUNTER_MATVECS]);
    ymir_pressure_vec_perf_counter_add_matvecs (
        &slabs_perf_counter[SLABS_PERF_COUNTER_MATVECS]);
    ymir_bbt_perf_counter_add_matvecs (
        &slabs_perf_counter[SLABS_PERF_COUNTER_MATVECS]);
    ymir_stokes_op_perf_counter_add_matvecs (
        &slabs_perf_counter[SLABS_PERF_COUNTER_MATVECS]);

    /* gather and print slabs_io statistics */
    slabs_io_perf_counter_gather (mpicomm, print_wtime,
                                  0 /* #calls */, 0 /* flops */);
    slabs_io_perf_counter_print ();

    /* gather and print slabs statistics */
    slabs_perf_counter_gather (mpicomm, print_wtime, print_n_calls,
                               print_flops);
    slabs_perf_counter_print ();
  }

  /* destroy options */
  ymir_options_global_destroy ();

  /* print that this function is ending */
  YMIR_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

  /* finalize rhea */
  rhea_finalize ();

  return 0;
}

/******************************************************************************
 * Tests
 *****************************************************************************/

#define SL_TEST_BRICK_DX (1)
#define SL_TEST_BRICK_DY (16)
#define SL_TEST_BRICK_DZ (16)
#define SL_TEST_BRICK_INNER_PLUME_DECAY (0.2)
#define SL_TEST_BRICK_INNER_PLUME_SCALING (-0.5)

/**
 * Sets input field:
 *
 *   u (x,y,z) = sin(2*pi/lx * x) * sin(2*pi*ly * y) * sin(2*pi*lz * z)
 *
 * where lx, ly, and lz are domain side lengths.
 */
static void
slabs_test_stiffness_set_input_fn (double *val, double x, double y, double z,
                                   ymir_locidx_t nid, void *data)
{
  slabs_physics_options_t *physics_options = (slabs_physics_options_t *) data;
  const double        lx = (double) physics_options->domain_brick_dx;
  const double        ly = (double) physics_options->domain_brick_dy;
  const double        lz = (double) physics_options->domain_brick_dz;
  const double        sinx = sin (2.0*M_PI/lx * x);
  const double        siny = sin (2.0*M_PI/ly * y);
  const double        sinz = sin (2.0*M_PI/lz * z);

  *val = sinx*siny*sinz;
  YMIR_ASSERT (isfinite (*val));
}

/**
 * Sets output field for constant coefficient:
 *
 *   - div (grad u (x,y,z))
 *   =
 *   4*pi^2 * (1/lx^2 + 1/ly^2 + 1/lz^2) *
 *   sin(2*pi/lx * x) * sin(2*pi*ly * y) * sin(2*pi*lz * z)
 */
static void
slabs_test_stiffness_set_output_coeff_const_fn (double *val, double x,
                                                double y, double z,
                                                ymir_locidx_t nid, void *data)
{
  slabs_physics_options_t *physics_options = (slabs_physics_options_t *) data;
  const double        lx = (double) physics_options->domain_brick_dx;
  const double        ly = (double) physics_options->domain_brick_dy;
  const double        lz = (double) physics_options->domain_brick_dz;
  const double        sinx = sin (2.0*M_PI/lx * x);
  const double        siny = sin (2.0*M_PI/ly * y);
  const double        sinz = sin (2.0*M_PI/lz * z);

  *val = 4.0*M_PI*M_PI * (1.0/lx/lx + 1.0/ly/ly + 1.0/lz/lz) * sinx*siny*sinz;
  YMIR_ASSERT (isfinite (*val));
}

/**
 * Sets variable scalar coefficient:
 *
 *   k(x,y,z) = exp (E * sqrt((x/lx)^2 + (y/ly)^2 + (z/lz)^2))
 */
static void
slabs_test_stiffness_set_coeff_scal_fn (double *val, double x, double y,
                                        double z, ymir_locidx_t nid, void *data)
{
  slabs_physics_options_t *physics_options = (slabs_physics_options_t *) data;
  const double        lx = (double) physics_options->domain_brick_dx;
  const double        ly = (double) physics_options->domain_brick_dy;
  const double        lz = (double) physics_options->domain_brick_dz;
  const double        norm = sqrt (x*x/(lx*lx) + y*y/(ly*ly) + z*z/(lz*lz));
  const double        E = 5.3;

  *val = exp (E*norm);
  YMIR_ASSERT (isfinite (*val));
}

/**
 * Sets output field for variable scalar coefficient:
 *
 *   - div (k(x,y,z) grad u(x,y,z))
 *   =
 *   4*pi^2 * (1/lx^2 + 1/ly^2 + 1/lz^2) * k(x,y,z) *
 *   sin(2*pi/lx * x) * sin(2*pi*ly * y) * sin(2*pi*lz * z)
 *   -
 *   2*pi *
 *   (
 *     (d/dx k)/lx * cos(2*pi/lx * x) * sin(2*pi/ly * y) * sin(2*pi/lz * z) +
 *     (d/dx k)/ly * sin(2*pi/lx * x) * cos(2*pi/ly * y) * sin(2*pi/lz * z) +
 *     (d/dx k)/lz * sin(2*pi/lx * x) * sin(2*pi/ly * y) * cos(2*pi/lz * z)
 *   )
 */
static void
slabs_test_stiffness_set_output_coeff_scal_fn (double *val, double x,
                                               double y, double z,
                                               ymir_locidx_t nid, void *data)
{
  slabs_physics_options_t *physics_options = (slabs_physics_options_t *) data;
  const double        lx = (double) physics_options->domain_brick_dx;
  const double        ly = (double) physics_options->domain_brick_dy;
  const double        lz = (double) physics_options->domain_brick_dz;
  const double        sinx = sin (2.0*M_PI/lx * x);
  const double        siny = sin (2.0*M_PI/ly * y);
  const double        sinz = sin (2.0*M_PI/lz * z);
  const double        cosx = cos (2.0*M_PI/lx * x);
  const double        cosy = cos (2.0*M_PI/ly * y);
  const double        cosz = cos (2.0*M_PI/lz * z);
  const double        norm = sqrt (x*x/(lx*lx) + y*y/(ly*ly) + z*z/(lz*lz));
  const double        E = 5.3;

  /* return zero if (x,y,z) = (0,0,0) */
  if (norm <= SC_EPS) {
    *val = 0.0;
    return;
  }

  *val = 4.0*M_PI*M_PI * (1.0/lx/lx + 1.0/ly/ly + 1.0/lz/lz) *
         exp (E*norm) * sinx*siny*sinz
         -
         2.0*M_PI * (
           E*x/(lx*lx)/norm * exp (E*norm) / lx * cosx*siny*sinz +
           E*y/(ly*ly)/norm * exp (E*norm) / ly * sinx*cosy*sinz +
           E*y/(ly*ly)/norm * exp (E*norm) / ly * sinx*siny*cosz
         );
  YMIR_ASSERT (isfinite (*val));
}

/**
 * Sets variable diagonal tensor coefficient:
 *
 *   K(x,y,z) = diag(k1, k2, k3),
 *
 *   k1(x,y,z) = exp (E * sqrt(((x-lx)/lx)^2 + (y/ly)^2 + (z/lz)^2))
 *   k2(x,y,z) = exp (E * sqrt((x/lx)^2 + ((y-ly)/ly)^2 + (z/lz)^2))
 *   k3(x,y,z) = exp (E * sqrt((x/lx)^2 + (y/ly)^2 + ((z-lz)/lz)^2))
 */
static void
slabs_test_stiffness_set_coeff_diag_fn (double *val, double x, double y,
                                        double z, ymir_locidx_t nid, void *data)
{
  slabs_physics_options_t *physics_options = (slabs_physics_options_t *) data;
  const double        lx = (double) physics_options->domain_brick_dx;
  const double        ly = (double) physics_options->domain_brick_dy;
  const double        lz = (double) physics_options->domain_brick_dz;
  const double        sx = x - lx;
  const double        sy = y - ly;
  const double        sz = z - lz;
  const double        norm_sx = sqrt (sx*sx/(lx*lx) + y*y/(ly*ly) +
                                      z*z/(lz*lz));
  const double        norm_sy = sqrt (x*x/(lx*lx) + sy*sy/(ly*ly) +
                                      z*z/(lz*lz));
  const double        norm_sz = sqrt (x*x/(lx*lx) + y*y/(ly*ly) +
                                      sz*sz/(lz*lz));
  const double        E = 2.6;

  val[0] = exp (E*norm_sx);
  val[1] = exp (E*norm_sy);
  val[2] = exp (E*norm_sz);
  YMIR_ASSERT (isfinite (val[0]) && isfinite (val[1]) && isfinite (val[2]));
}

/**
 * Sets output field for variable diagonal tensor coefficient:
 *
 *   - div (diag(k1, k2, k3) grad u(x,y,z))
 *   =
 *   4*pi^2 * (k1/lx^2 + k2/ly^2 + k3/lz^2) *
 *   sin(2*pi/lx * x) * sin(2*pi*ly * y) * sin(2*pi*lz * z)
 *   -
 *   2*pi *
 *   (
 *     (d/dx k1)/lx * cos(2*pi/lx * x) * sin(2*pi/ly * y) * sin(2*pi/lz * z) +
 *     (d/dx k2)/ly * sin(2*pi/lx * x) * cos(2*pi/ly * y) * sin(2*pi/lz * z) +
 *     (d/dx k3)/lz * sin(2*pi/lx * x) * sin(2*pi/ly * y) * cos(2*pi/lz * z)
 *   )
 */
static void
slabs_test_stiffness_set_output_coeff_diag_fn (double *val, double x,
                                               double y, double z,
                                               ymir_locidx_t nid, void *data)
{
  slabs_physics_options_t *physics_options = (slabs_physics_options_t *) data;
  const double        lx = (double) physics_options->domain_brick_dx;
  const double        ly = (double) physics_options->domain_brick_dy;
  const double        lz = (double) physics_options->domain_brick_dz;
  const double        sinx = sin (2.0*M_PI/lx * x);
  const double        siny = sin (2.0*M_PI/ly * y);
  const double        sinz = sin (2.0*M_PI/lz * z);
  const double        cosx = cos (2.0*M_PI/lx * x);
  const double        cosy = cos (2.0*M_PI/ly * y);
  const double        cosz = cos (2.0*M_PI/lz * z);
  double              sx = x - lx;
  double              sy = y - ly;
  double              sz = z - lz;
  double              norm_sx = sqrt (sx*sx/(lx*lx) + y*y/(ly*ly) +
                                      z*z/(lz*lz));
  double              norm_sy = sqrt (x*x/(lx*lx) + sy*sy/(ly*ly) +
                                      z*z/(lz*lz));
  double              norm_sz = sqrt (x*x/(lx*lx) + y*y/(ly*ly) +
                                      sz*sz/(lz*lz));
  const double        E = 2.6;

  /* compute first part of value */
  *val = 4.0*M_PI*M_PI * (
           exp (E*norm_sx) / (lx*lx) +
           exp (E*norm_sy) / (ly*ly) +
           exp (E*norm_sz) / (lz*lz)
         ) * sinx*siny*sinz;
  YMIR_ASSERT (isfinite (*val));

  /* avoid division by zero */
  if (norm_sx <= SC_EPS) {
    sx = 0.0;
    norm_sx = 1.0;
  }
  if (norm_sy <= SC_EPS) {
    sy = 0.0;
    norm_sy = 1.0;
  }
  if (norm_sz <= SC_EPS) {
    sz = 0.0;
    norm_sz = 1.0;
  }

  /* compute second part of value */
  *val -= 2.0*M_PI * (
            E*sx/(lx*lx)/norm_sx * exp (E*norm_sx) / lx * cosx*siny*sinz +
            E*sy/(ly*ly)/norm_sy * exp (E*norm_sy) / ly * sinx*cosy*sinz +
            E*sy/(ly*ly)/norm_sz * exp (E*norm_sz) / ly * sinx*siny*cosz
          );
  YMIR_ASSERT (isfinite (*val));
}

/**
 *
 */
static void
slabs_test_stiffness (MPI_Comm mpicomm,
                      const slabs_domain_shape_t domain_shape,
                      const int minlevel,
                      const int maxlevel,
                      char *refine,
                      const int order)
{
  const char         *this_fn_name = "slabs_test_stiffness";

  slabs_physics_options_t  physics_options;
  slabs_discr_options_t  discr_options;
  p8est_t            *p8est;
  mangll_t           *mangll;
  mangll_cnodes_t    *cnodes;
  ymir_mesh_t        *mesh;
  ymir_pressure_elem_t  *press_elem;
  ymir_vec_t         *in, *out, *out_exact, *out_exact_mass, *error;
  int                 dirfaces[6];
  int                 faceid;
  ymir_stiff_op_t    *stiff_op;
  ymir_vec_t         *coeff;
#ifdef YMIR_PETSC
  ymir_matrix_t      *ymat;
#endif
  double              abs_error, rel_error;
  char                filename[BUFSIZ];

  /* test is only implemented for the cube and brick domain */
  if (domain_shape != SL_DOMAIN_CUBE && domain_shape != SL_DOMAIN_BRICK) {
    return;
  }

  YMIR_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  /* set required options */
  physics_options.domain_shape = domain_shape;
  if (domain_shape == SL_DOMAIN_CUBE) {
    physics_options.domain_brick_dx = 1;
    physics_options.domain_brick_dy = 1;
    physics_options.domain_brick_dz = 1;
  }
  else if (domain_shape == SL_DOMAIN_BRICK) {
    physics_options.domain_brick_dx = SL_TEST_BRICK_DX;
    physics_options.domain_brick_dy = SL_TEST_BRICK_DY;
    physics_options.domain_brick_dz = SL_TEST_BRICK_DZ;
  }
  //physics_options.bc_type = SL_VEL_BC_DIRICHLET_ALL;
  physics_options.bc_type = SL_VEL_BC_DIRICHLET_NORM;
  physics_options.temperature_type = SL_TEMP_NONE;
  physics_options.viscosity_type = SL_VISCOSITY_CONST;
  physics_options.viscosity_scaling = 1.0;
  physics_options.rhs_scaling = 1.0;
  physics_options.plume_type = SL_PLUME_INNER;
  physics_options.plume_center_x = 0.5 * physics_options.domain_brick_dx;
  physics_options.plume_center_y = 0.5 * physics_options.domain_brick_dy;
  physics_options.plume_center_z = 0.5 * physics_options.domain_brick_dz;
  if (domain_shape == SL_DOMAIN_CUBE) {
    physics_options.plume_decay = 10.0;
    physics_options.plume_scaling = -1.0;
  }
  else if (domain_shape == SL_DOMAIN_BRICK) {
    physics_options.plume_decay = SL_TEST_BRICK_INNER_PLUME_DECAY;
    physics_options.plume_scaling = SL_TEST_BRICK_INNER_PLUME_SCALING;
  }
  discr_options.minlevel = (int8_t) minlevel;
  discr_options.maxlevel = (int8_t) maxlevel;
  discr_options.refine = refine;
  discr_options.X_fn = slabs_discr_identity_X;
  discr_options.order = order;

  /* create p4est and boundary variables */
  p8est = slabs_discr_p8est_new (mpicomm, &physics_options, &discr_options);
  slabs_discr_options_set_boundary (&discr_options, p8est, &physics_options);

  /* create mesh */
  slabs_discr_mangll_and_cnodes_new (&mangll, &cnodes, p8est, &discr_options);
  slabs_discr_ymir_new (&mesh, &press_elem, mangll, cnodes, &discr_options);

  /* create work vectors */
  in = ymir_cvec_new (mesh, 1);
  out = ymir_cvec_new (mesh, 1);
  out_exact = ymir_cvec_new (mesh, 1);
  out_exact_mass = ymir_cvec_new (mesh, 1);
  error = ymir_cvec_new (mesh, 1);

  /* set input field */
  ymir_cvec_set_function (in, slabs_test_stiffness_set_input_fn,
                          &physics_options);
  YMIR_ASSERT (sc_dmatrix_is_valid (in->data));

  /* set Dirichlet boundaries */
  for (faceid = 0; faceid < 6; faceid++) {
    dirfaces[faceid] = 1;
  }

  /*
   * Constant Coefficient
   */

  /* create stiffness operator */
  stiff_op = ymir_stiff_op_new (mesh, dirfaces);

  /* apply stiffness operator */
  ymir_stiff_op_apply (in, out, stiff_op);
  YMIR_ASSERT (sc_dmatrix_is_valid (out->data));

  /* compute exact result of stiffness apply */
  ymir_cvec_set_function (
      out_exact, slabs_test_stiffness_set_output_coeff_const_fn,
      &physics_options);
  ymir_mass_apply (out_exact, out_exact_mass);
  YMIR_ASSERT (sc_dmatrix_is_valid (out_exact->data));
  YMIR_ASSERT (sc_dmatrix_is_valid (out_exact_mass->data));

  /* compute & print matrix-free error */
  ymir_vec_copy (out, error)
  ymir_vec_add (-1.0, out_exact_mass, error);
  abs_error = ymir_vec_norm (error);
  rel_error = abs_error / ymir_vec_norm (out_exact_mass);
  YMIR_GLOBAL_INFOF ("%s: coeff const: abs error %1.3e, rel error %1.3e\n",
                     this_fn_name, abs_error, rel_error);

  /* write vtk output */
  snprintf (filename, BUFSIZ, "%s_coeff_const", this_fn_name);
  ymir_vtk_write (mesh, filename, in, "in", out, "out", error, "error",
                  out_exact, "out_exact", out_exact_mass, "out_exact_mass",
                  NULL);

#ifdef YMIR_PETSC
  /* assemble stiffness matrix */
  ymat = NULL;
  if (ymir_submesh_use) {
    ymir_stiffness_submesh_matrix (mesh, mesh->submesh,
                                   YMIR_MESH_PETSCLAYOUT_NONE,
                                   stiff_op->dirscale, stiff_op->coeff,
                                   stiff_op->mass_scaling,
                                   stiff_op->mass_coeff, &ymat,
                                   YMIR_STIFF_OP_MATRIX_FULL);
  }
  else {
    ymir_stiffness_matrix (mesh, YMIR_MESH_PETSCLAYOUT_NONE,
                           stiff_op->dirscale, stiff_op->coeff,
                           stiff_op->mass_scaling, stiff_op->mass_coeff,
                           &ymat, YMIR_STIFF_OP_MATRIX_FULL, 1, 0);
  }

  /* apply matrix */
  ymir_matrix_apply (in, out, ymat);
  YMIR_ASSERT (sc_dmatrix_is_valid (out->dataown));
  YMIR_ASSERT (sc_dmatrix_is_valid (out->coff));

  /* compute & print matrix error */
  ymir_vec_copy (out, error)
  ymir_vec_add (-1.0, out_exact_mass, error);
  abs_error = ymir_vec_norm (error);
  rel_error = abs_error / ymir_vec_norm (out_exact_mass);
  YMIR_GLOBAL_INFOF ("%s: coeff const, mat: abs error %1.3e, rel error %1.3e\n",
                     this_fn_name, abs_error, rel_error);
#endif

  /* destroy */
  ymir_stiff_op_destroy (stiff_op);
#ifdef YMIR_PETSC
  ymir_matrix_destroy (ymat);
#endif

  /*
   * Variable Scalar Coefficient
   */

  /* create stiffness operator */
  stiff_op = ymir_stiff_op_new (mesh, dirfaces);

  /* set coefficient */
  coeff = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
  ymir_dvec_set_function (
      coeff, slabs_test_stiffness_set_coeff_scal_fn, &physics_options);
  YMIR_ASSERT (sc_dmatrix_is_valid (coeff->data));
  ymir_stiff_op_set_coeff (stiff_op, coeff);

  /* apply stiffness operator */
  ymir_stiff_op_apply (in, out, stiff_op);
  YMIR_ASSERT (sc_dmatrix_is_valid (out->data));

  /* compute exact result of stiffness apply */
  ymir_cvec_set_function (
      out_exact, slabs_test_stiffness_set_output_coeff_scal_fn,
      &physics_options);
  ymir_mass_apply (out_exact, out_exact_mass);
  YMIR_ASSERT (sc_dmatrix_is_valid (out_exact->data));
  YMIR_ASSERT (sc_dmatrix_is_valid (out_exact_mass->data));

  /* compute & print matrix-free error */
  ymir_vec_copy (out, error)
  ymir_vec_add (-1.0, out_exact_mass, error);
  abs_error = ymir_vec_norm (error);
  rel_error = abs_error / ymir_vec_norm (out_exact_mass);
  YMIR_GLOBAL_INFOF ("%s: coeff scal: abs error %1.3e, rel error %1.3e\n",
                     this_fn_name, abs_error, rel_error);

  /* write vtk output */
  snprintf (filename, BUFSIZ, "%s_coeff_scal", this_fn_name);
  ymir_vtk_write (mesh, filename, coeff, "coeff", in, "in", out, "out",
                  out_exact, "out_exact", out_exact_mass, "out_exact_mass",
                  error, "error", NULL);

#ifdef YMIR_PETSC
  /* assemble stiffness matrix */
  ymat = NULL;
  if (ymir_submesh_use) {
    ymir_stiffness_submesh_matrix (mesh, mesh->submesh,
                                   YMIR_MESH_PETSCLAYOUT_NONE,
                                   stiff_op->dirscale, stiff_op->coeff,
                                   stiff_op->mass_scaling,
                                   stiff_op->mass_coeff, &ymat,
                                   YMIR_STIFF_OP_MATRIX_FULL);
  }
  else {
    ymir_stiffness_matrix (mesh, YMIR_MESH_PETSCLAYOUT_NONE,
                           stiff_op->dirscale, stiff_op->coeff,
                           stiff_op->mass_scaling, stiff_op->mass_coeff,
                           &ymat, YMIR_STIFF_OP_MATRIX_FULL, 1, 0);
  }

  /* apply matrix */
  ymir_matrix_apply (in, out, ymat);
  YMIR_ASSERT (sc_dmatrix_is_valid (out->dataown));
  YMIR_ASSERT (sc_dmatrix_is_valid (out->coff));

  /* compute & print matrix error */
  ymir_vec_copy (out, error)
  ymir_vec_add (-1.0, out_exact_mass, error);
  abs_error = ymir_vec_norm (error);
  rel_error = abs_error / ymir_vec_norm (out_exact_mass);
  YMIR_GLOBAL_INFOF ("%s: coeff scal, mat: abs error %1.3e, rel error %1.3e\n",
                     this_fn_name, abs_error, rel_error);
#endif

  /* destroy */
  ymir_vec_destroy (coeff);
  ymir_stiff_op_destroy (stiff_op);
#ifdef YMIR_PETSC
  ymir_matrix_destroy (ymat);
#endif

  /*
   * Variable, anisotropic coefficient (diagonal tensor)
   */

  /* create stiffness operator */
  stiff_op = ymir_stiff_op_new (mesh, dirfaces);

  /* set coefficient */
  coeff = ymir_dvec_new (mesh, 3, YMIR_GAUSS_NODE);
  ymir_dvec_set_function (
      coeff, slabs_test_stiffness_set_coeff_diag_fn, &physics_options);
  YMIR_ASSERT (sc_dmatrix_is_valid (coeff->data));
  ymir_stiff_op_set_coeff (stiff_op, coeff);

  /* apply stiffness operator */
  ymir_stiff_op_apply (in, out, stiff_op);
  YMIR_ASSERT (sc_dmatrix_is_valid (out->data));

  /* compute exact result of stiffness apply */
  ymir_cvec_set_function (
      out_exact, slabs_test_stiffness_set_output_coeff_diag_fn,
      &physics_options);
  ymir_mass_apply (out_exact, out_exact_mass);
  YMIR_ASSERT (sc_dmatrix_is_valid (out_exact->data));
  YMIR_ASSERT (sc_dmatrix_is_valid (out_exact_mass->data));

  /* compute & print matrix-free error */
  ymir_vec_copy (out, error)
  ymir_vec_add (-1.0, out_exact_mass, error);
  abs_error = ymir_vec_norm (error);
  rel_error = abs_error / ymir_vec_norm (out_exact_mass);
  YMIR_GLOBAL_INFOF ("%s: coeff diag: abs error %1.3e, rel error %1.3e\n",
                     this_fn_name, abs_error, rel_error);

  /* write vtk output */
  snprintf (filename, BUFSIZ, "%s_coeff_diag", this_fn_name);
  ymir_vtk_write (mesh, filename, coeff, "coeff", in, "in", out, "out",
                  out_exact, "out_exact", out_exact_mass, "out_exact_mass",
                  error, "error", NULL);

#ifdef YMIR_PETSC
  /* assemble stiffness matrix */
  ymat = NULL;
  if (ymir_submesh_use) {
    ymir_stiffness_submesh_matrix (mesh, mesh->submesh,
                                   YMIR_MESH_PETSCLAYOUT_NONE,
                                   stiff_op->dirscale, stiff_op->coeff,
                                   stiff_op->mass_scaling,
                                   stiff_op->mass_coeff, &ymat,
                                   YMIR_STIFF_OP_MATRIX_FULL);
  }
  else {
    ymir_stiffness_matrix (mesh, YMIR_MESH_PETSCLAYOUT_NONE,
                           stiff_op->dirscale, stiff_op->coeff,
                           stiff_op->mass_scaling, stiff_op->mass_coeff,
                           &ymat, YMIR_STIFF_OP_MATRIX_FULL, 1, 0);
  }

  /* apply matrix */
  ymir_matrix_apply (in, out, ymat);
  YMIR_ASSERT (sc_dmatrix_is_valid (out->dataown));
  YMIR_ASSERT (sc_dmatrix_is_valid (out->coff));

  /* compute & print matrix error */
  ymir_vec_copy (out, error)
  ymir_vec_add (-1.0, out_exact_mass, error);
  abs_error = ymir_vec_norm (error);
  rel_error = abs_error / ymir_vec_norm (out_exact_mass);
  YMIR_GLOBAL_INFOF ("%s: coeff diag, mat: abs error %1.3e, rel error %1.3e\n",
                     this_fn_name, abs_error, rel_error);
#endif

  /* destroy */
  ymir_vec_destroy (coeff);
  ymir_stiff_op_destroy (stiff_op);
#ifdef YMIR_PETSC
  ymir_matrix_destroy (ymat);
#endif

  /*
   * Variable, anisotropic coefficient (full symmetric tensor)
   */

  //TODO

  /* destroy */
  ymir_vec_destroy (in);
  ymir_vec_destroy (out);
  ymir_vec_destroy (out_exact);
  ymir_vec_destroy (out_exact_mass);
  ymir_vec_destroy (error);
  slabs_discr_ymir_mangll_destroy (mesh, press_elem);
  slabs_discr_p8est_destroy (p8est);
  slabs_discr_options_clear_boundary (&discr_options);

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

/**
 *
 */
static void
slabs_test_stress (MPI_Comm mpicomm,
                   const slabs_domain_shape_t domain_shape,
                   const int minlevel,
                   const int maxlevel,
                   char *refine,
                   const int order)
{
  const char         *this_fn_name = "slabs_test_stress";

  slabs_physics_options_t  physics_options;
  slabs_discr_options_t  discr_options;
  p8est_t            *p8est;
  mangll_t           *mangll;
  mangll_cnodes_t    *cnodes;
  ymir_mesh_t        *mesh;
  ymir_pressure_elem_t  *press_elem;

  ymir_vec_t         *coeff;
  ymir_vel_dir_t     *vel_dir;
  ymir_stress_op_t   *stress_op;
  ymir_vec_t         *in, *out;
  double              min, max;

  /* test is only implemented for the cube and brick domain */
  if (domain_shape != SL_DOMAIN_CUBE && domain_shape != SL_DOMAIN_BRICK) {
    return;
  }

  YMIR_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  /* set required options */
  physics_options.domain_shape = domain_shape;
  if (domain_shape == SL_DOMAIN_CUBE) {
    physics_options.domain_brick_dx = 1;
    physics_options.domain_brick_dy = 1;
    physics_options.domain_brick_dz = 1;
  }
  else if (domain_shape == SL_DOMAIN_BRICK) {
    physics_options.domain_brick_dx = SL_TEST_BRICK_DX;
    physics_options.domain_brick_dy = SL_TEST_BRICK_DY;
    physics_options.domain_brick_dz = SL_TEST_BRICK_DZ;
  }
  physics_options.bc_type = SL_VEL_BC_DIRICHLET_NORM;
  physics_options.temperature_type = SL_TEMP_NONE;
  physics_options.viscosity_type = SL_VISCOSITY_CONST;
  physics_options.viscosity_scaling = 1.0;
  physics_options.rhs_scaling = 1.0;
  physics_options.plume_type = SL_PLUME_INNER;
  physics_options.plume_center_x = 0.5 * physics_options.domain_brick_dx;
  physics_options.plume_center_y = 0.5 * physics_options.domain_brick_dy;
  physics_options.plume_center_z = 0.5 * physics_options.domain_brick_dz;
  if (domain_shape == SL_DOMAIN_CUBE) {
    physics_options.plume_decay = 10.0;
    physics_options.plume_scaling = -1.0;
  }
  else if (domain_shape == SL_DOMAIN_BRICK) {
    physics_options.plume_decay = SL_TEST_BRICK_INNER_PLUME_DECAY;
    physics_options.plume_scaling = SL_TEST_BRICK_INNER_PLUME_SCALING;
  }
  discr_options.minlevel = (int8_t) minlevel;
  discr_options.maxlevel = (int8_t) maxlevel;
  discr_options.refine = refine;
  discr_options.X_fn = slabs_discr_identity_X;
  discr_options.order = order;

  /* create p4est and boundary variables */
  p8est = slabs_discr_p8est_new (mpicomm, &physics_options, &discr_options);
  slabs_discr_options_set_boundary (&discr_options, p8est, &physics_options);

  /* create mesh */
  slabs_discr_mangll_and_cnodes_new (&mangll, &cnodes, p8est, &discr_options);
  slabs_discr_ymir_new (&mesh, &press_elem, mangll, cnodes, &discr_options);

  /* create coefficient */
  coeff = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
  ymir_vec_set_value (coeff, 1.0);

  /* set Dirichlet boundaries */
  vel_dir = slabs_set_dirichlet_bc (mesh, NULL, &physics_options);

  /* create viscous stress operator */
  stress_op = ymir_stress_op_new (coeff, vel_dir, NULL, NULL, NULL);

  /* create work vectors */
  in = ymir_cvec_new (mesh, 3);
  out = ymir_cvec_new (mesh, 3);

  /*
   * Constant input, no Dirichlet boundary
   */

  /* set constant input vector */
  ymir_vec_set_value (in, 1.0);

  /* apply stress operator */
  stress_op->skip_dir = 1; //TODO deprecated
  ymir_stress_op_apply (in, out, stress_op);
  stress_op->skip_dir = 0; //TODO deprecated
  YMIR_ASSERT (sc_dmatrix_is_valid (out->data));

  /* compute min and max absolute value of output vector */
  ymir_vec_fabs (out, out);
  min = ymir_vec_min_global (out);
  max = ymir_vec_max_global (out);
  YMIR_GLOBAL_INFOF ("%s: rowsum, no Dir BC: min %1.3e, max %1.3e\n",
                     this_fn_name, min, max);

  /*
   * Constant input, with Dirichlet boundary in normal direction
   */

  /* set constant input vector */
  ymir_vec_set_value (in, 1.0);

  /* apply stress operator */
  ymir_stress_op_apply (in, out, stress_op);
  YMIR_ASSERT (sc_dmatrix_is_valid (out->data));

  /* zero out boundary nodes */
  {
    const ymir_locidx_t offset = 0;
    const ymir_locidx_t count = mesh->cnodes->owncount;
    ymir_locidx_t       cnid;
    sc_dmatrix_t       *node_mat = sc_dmatrix_new (1, out->ncfields);

    for (cnid = offset; cnid < (offset + count); cnid++) {
      if (vel_dir->rank[cnid]) {
        ymir_cvec_get_node (out, node_mat, cnid, YMIR_WRITE);
        sc_dmatrix_set_value (node_mat, 0.0);
        ymir_cvec_set_node (out, node_mat, cnid, YMIR_SET);
      }
    }

    sc_dmatrix_destroy (node_mat);
  }

  /* compute min and max absolute value of output vector */
  ymir_vec_fabs (out, out);
  min = ymir_vec_min_global (out);
  max = ymir_vec_max_global (out);
  YMIR_GLOBAL_INFOF ("%s: rowsum, norm Dir BC: min %1.3e, max %1.3e "
                     "(non-boundary nodes)\n", this_fn_name, min, max);

  /* destroy */
  ymir_vec_destroy (in);
  ymir_vec_destroy (out);
  ymir_stress_op_destroy (stress_op);
  ymir_vel_dir_destroy (vel_dir);
  ymir_vec_destroy (coeff);
  slabs_discr_ymir_mangll_destroy (mesh, press_elem);
  slabs_discr_p8est_destroy (p8est);
  slabs_discr_options_clear_boundary (&discr_options);

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

/**
 *
 */
static int
slabs_test_refine_all_fn (p8est_t *p8est, p4est_topidx_t tree,
                          p8est_quadrant_t *quadrant)
{
  return 1;
}

/* parameter list for test of norms */
typedef struct slabs_test_set_data
{
  int                 n_fields;
  int                 brick_dx, brick_dy, brick_dz;
  double              cx, cy, cz;
}
slabs_test_set_data_t;

/**
 *
 */
static void
slabs_test_norm_set_fn (double *val, double x, double y, double z,
                        ymir_locidx_t nid, void *data)
{
  slabs_test_set_data_t  *d = (slabs_test_set_data_t *) data;
  const double        cx = d->cx;
  const double        cy = d->cy;
  const double        cz = d->cz;
  const double        brick_dx = (double) d->brick_dx;
  const double        brick_dy = (double) d->brick_dy;
  const double        brick_dz = (double) d->brick_dz;

  val[0] = 0.0;
  val[1] = 0.0;
  val[2] = exp (M_PI * x) * sin (cx / brick_dx * M_PI * x)
                          * sin (cy / brick_dy * M_PI * y)
                          * sin (cz / brick_dz * M_PI * z);
}

/**
 *
 */
static void
slabs_test_mesh_independence_stokes_norm (MPI_Comm mpicomm,
                                          const slabs_domain_shape_t
                                            domain_shape,
                                          const int minlevel,
                                          const int maxlevel,
                                          char *refine,
                                          const int order,
                                          double mass_scaling,
                                          const int max_refine_steps)
{
  const char         *this_fn_name = "slabs_test_mesh_independence_stokes_norm";

  slabs_physics_options_t  physics_options;
  slabs_discr_options_t  discr_options;
  slabs_test_set_data_t  data;
  p8est_t            *p8est;
  mangll_t           *mangll;
  mangll_cnodes_t    *cnodes;
  ymir_mesh_t        *mesh;
  ymir_pressure_elem_t  *press_elem;
  slabs_stokes_state_t  *state;
  slabs_lin_stokes_problem_t  *lin_stokes;
  ymir_Hminus1_norm_op_t  *norm_op;
  ymir_vec_t         *guess_up, *guess_u, *guess_p;
  ymir_vec_t         *residual_up, *residual_u, *residual_p;
  int                 refine_step;
  sc_dmatrix_t       *results;
  char                filename[BUFSIZ];

  /* test is only implemented for the cube and brick domain */
  if (domain_shape != SL_DOMAIN_CUBE && domain_shape != SL_DOMAIN_BRICK) {
    return;
  }

  YMIR_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  /* set required options */
  physics_options.domain_shape = domain_shape;
  if (domain_shape == SL_DOMAIN_CUBE) {
    physics_options.domain_brick_dx = 1;
    physics_options.domain_brick_dy = 1;
    physics_options.domain_brick_dz = 1;
  }
  else if (domain_shape == SL_DOMAIN_BRICK) {
    physics_options.domain_brick_dx = SL_TEST_BRICK_DX;
    physics_options.domain_brick_dy = SL_TEST_BRICK_DY;
    physics_options.domain_brick_dz = SL_TEST_BRICK_DZ;
  }
  //physics_options.bc_type = SL_VEL_BC_DIRICHLET_ALL;
  physics_options.bc_type = SL_VEL_BC_DIRICHLET_NORM;
  physics_options.temperature_type = SL_TEMP_NONE;
  physics_options.viscosity_type = SL_VISCOSITY_CONST;
  physics_options.viscosity_scaling = 1.0;
  physics_options.rhs_scaling = 1.0;
  physics_options.plume_type = SL_PLUME_INNER;
  physics_options.plume_center_x = 0.5 * physics_options.domain_brick_dx;
  physics_options.plume_center_y = 0.5 * physics_options.domain_brick_dy;
  physics_options.plume_center_z = 0.5 * physics_options.domain_brick_dz;
  if (domain_shape == SL_DOMAIN_CUBE) {
    physics_options.plume_decay = 10.0;
    physics_options.plume_scaling = -1.0;
  }
  else if (domain_shape == SL_DOMAIN_BRICK) {
    physics_options.plume_decay = SL_TEST_BRICK_INNER_PLUME_DECAY;
    physics_options.plume_scaling = SL_TEST_BRICK_INNER_PLUME_SCALING;
  }
  discr_options.minlevel = (int8_t) minlevel;
  discr_options.maxlevel = (int8_t) maxlevel;
  discr_options.refine = refine;
  discr_options.X_fn = slabs_discr_identity_X;
  discr_options.order = order;

  /* set data for test function */
  data.n_fields = 3;
  data.cx = 2.0;
  data.cy = 4.0;
  data.cz = 1.0;
  data.brick_dx = physics_options.domain_brick_dx;
  data.brick_dy = physics_options.domain_brick_dy;
  data.brick_dz = physics_options.domain_brick_dz;

  /* create p4est and boundary variables */
  p8est = slabs_discr_p8est_new (mpicomm, &physics_options, &discr_options);
  slabs_discr_options_set_boundary (&discr_options, p8est, &physics_options);

  /* create vector to store results */
  results = sc_dmatrix_new (max_refine_steps, 5);

  for (refine_step = 0; refine_step < max_refine_steps; refine_step++) {
    /* refine p4est */
    p8est_refine (p8est, 0, slabs_test_refine_all_fn, NULL);

    /* create mesh */
    slabs_discr_mangll_and_cnodes_new (&mangll, &cnodes, p8est, &discr_options);
    slabs_discr_ymir_new (&mesh, &press_elem, mangll, cnodes, &discr_options);

    /* create Stokes state */
    state = slabs_stokes_state_new (p8est);
    slabs_stokes_state_init_temp (state, cnodes);
    slabs_stokes_state_init_temp_vec (state, mesh);
    ymir_cvec_set_function (state->temp_vec, slabs_physics_temperature_set_fn,
                            &physics_options);
    slabs_stokes_state_init_vel_press (state, mesh, press_elem);

    /* create linear Stokes problem */
    lin_stokes = slabs_linear_stokes_problem_new (state, mesh, press_elem,
                                                  &physics_options);

    /* create H^-1 norm operator */
    norm_op = ymir_Hminus1_norm_op_new (mesh, mass_scaling);

    /* create vectors */
    guess_up = ymir_stokes_vec_new (mesh, press_elem);
    residual_up = ymir_stokes_vec_new (mesh, press_elem);

    /* set guess */
    ymir_vec_set_zero (guess_up);
    ymir_cvec_set_function (guess_up, slabs_test_norm_set_fn, &data);

    /* compute residual */
    // version 1:
    //ymir_mass_apply (guess_up, residual_up);
    // version 2:
    //ymir_stokes_op_apply (guess_up, residual_up, lin_stokes->stokes_op);
    // version 3:
    ymir_stokes_op_apply (guess_up, residual_up, lin_stokes->stokes_op);
    ymir_vec_add (-1.0, lin_stokes->rhs, residual_up);
    ymir_vec_scale (-1.0, residual_up);

    /* [0] compute l^2 norm of guess */
    results->e[refine_step][0] = slabs_norm (NULL, NULL, guess_up, NULL,
                                             SL_NORM_VEC_L2, NULL);

    /* [1] compute l^2 norm */
    results->e[refine_step][1] = slabs_norm (NULL, NULL, residual_up, NULL,
                                             SL_NORM_VEC_L2, NULL);

    /* [2] compute L^2 norm */
    results->e[refine_step][2] = slabs_norm (NULL, NULL, residual_up,
                                             press_elem, SL_NORM_FNC_L2, NULL);

    /* [3] compute H^-1 norm component-wise */
    results->e[refine_step][3] = slabs_norm (NULL, NULL, residual_up, NULL,
                                             SL_NORM_FNC_HMINUS1, norm_op);

    /* [4] compute H^-1 norm by "inverting" the Stokes operator */
    ymir_stokes_op_apply (guess_up, residual_up, lin_stokes->stokes_op);
    results->e[refine_step][4] = sqrt (ymir_cvec_innerprod (residual_up,
                                                            guess_up));

    /* output */
    guess_u = ymir_cvec_new (mesh, 3);
    guess_p = ymir_pressure_vec_new (mesh, press_elem);
    ymir_stokes_vec_get_components (guess_up, guess_u, guess_p, press_elem);
    residual_u = ymir_cvec_new (mesh, 3);
    residual_p = ymir_pressure_vec_new (mesh, press_elem);
    ymir_stokes_vec_get_components (residual_up, residual_u, residual_p,
                                    press_elem);

    snprintf (filename, BUFSIZ, "%s_l%i", this_fn_name,
              minlevel + refine_step + 1);
    ymir_vtk_write (mesh, filename,
                    guess_u, "guess_u", guess_p, "guess_p",
                    residual_u, "residual_u", residual_p, "residual_p",
                    lin_stokes->viscosity, "viscosity",
                    lin_stokes->rhs_u_point, "rhs_u_point",
                    NULL);

    snprintf (filename, BUFSIZ, "%s_l%i_guess", this_fn_name,
              minlevel + refine_step + 1);
    slabs_io_write_cvec_to_textfile (filename, guess_u,
                                     SL_CARTESIAN_COORDINATE);

    snprintf (filename, BUFSIZ, "%s_l%i_residual", this_fn_name,
              minlevel + refine_step + 1);
    slabs_io_write_cvec_to_textfile (filename, residual_u,
                                     SL_CARTESIAN_COORDINATE);

    ymir_vec_destroy (guess_u);
    ymir_vec_destroy (guess_p);
    ymir_vec_destroy (residual_u);
    ymir_vec_destroy (residual_p);

    /* hack to avoid SEGV when destroying stress PC (TODO fix) */
    ymir_vec_set_zero (guess_up);
    ymir_cvec_set_function (guess_up, slabs_test_norm_set_fn, &data);
    ymir_stokes_pc_solve (lin_stokes->rhs, guess_up, lin_stokes->stokes_pc,
                          1, 0.1, 0.1, 1, 1, NULL);

    /* destroy */
    ymir_vec_destroy (guess_up);
    ymir_vec_destroy (residual_up);
    ymir_Hminus1_norm_op_destroy (norm_op);
    slabs_linear_stokes_problem_destroy (lin_stokes);
    slabs_stokes_state_destroy (state);
    slabs_discr_ymir_mangll_destroy (mesh, press_elem);
  }

  /* destroy p4est and options */
  slabs_discr_p8est_destroy (p8est);
  slabs_discr_options_clear_boundary (&discr_options);

  /* output of results */
  for (refine_step = 0; refine_step < max_refine_steps; refine_step++) {
    YMIR_GLOBAL_INFOF ("%s: level %i, l^2 guess %1.3e, "
                       "l^2 %1.3e, L^2 %1.3e, H^-1 comp %1.3e, "
                       "(H^-1 full %1.3e)\n",
                       this_fn_name, minlevel + refine_step + 1,
                       results->e[refine_step][0],
                       results->e[refine_step][1],
                       results->e[refine_step][2],
                       results->e[refine_step][3],
                       results->e[refine_step][4]);
  }

  /* destroy results */
  sc_dmatrix_destroy (results);

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

/**
 *
 */
static void
slabs_test_interp_set_vel_fn (double *val, double x, double y, double z,
                              ymir_locidx_t nid, void *data)
{
  slabs_test_set_data_t  *d = (slabs_test_set_data_t *) data;
  const double        cx = d->cx;
  const double        cy = d->cy;
  const double        cz = d->cz;
  const double        brick_dx = (double) d->brick_dx;
  const double        brick_dy = (double) d->brick_dy;
  const double        brick_dz = (double) d->brick_dz;

  val[0] = exp (M_PI * x);
  val[1] =   cos (cx / brick_dx * M_PI * x)
           * cos (cy / brick_dy * M_PI * y)
           * cos (cz / brick_dz * M_PI * z);
  val[2] =   sin (cx / brick_dx * M_PI * x)
           * sin (cy / brick_dy * M_PI * y)
           * sin (cz / brick_dz * M_PI * z);
}

/**
 *
 */
static void
slabs_test_interp_set_press_fn (double *val, double *x, double *y, double *z,
                                ymir_locidx_t elem_id, void *data)
{
  slabs_test_set_data_t  *d = (slabs_test_set_data_t *) data;
  const int           n_fields = d->n_fields;
  const double        cx = d->cx;
  const double        cy = d->cy;
  const double        cz = d->cz;
  const double        brick_dx = (double) d->brick_dx;
  const double        brick_dy = (double) d->brick_dy;
  const double        brick_dz = (double) d->brick_dz;

  int                 fieldid;

  for (fieldid = 0; fieldid < n_fields; fieldid++) {
    val[fieldid] =   cos (cx / brick_dx * M_PI * x[fieldid])
                   * sin (cy / brick_dy * M_PI * y[fieldid])
                   * cos (cz / brick_dz * M_PI * z[fieldid]);
  }
}

/**
 *
 */
static void
slabs_test_interp_set_dual_fn (double *val, double x, double y, double z,
                               ymir_locidx_t nid, void *data)
{
  slabs_test_set_data_t  *d = (slabs_test_set_data_t *) data;
  const int           n_fields = d->n_fields;
  const double        cx = d->cx;
  const double        cy = d->cy;
  const double        cz = d->cz;
  const double        brick_dx = (double) d->brick_dx;
  const double        brick_dy = (double) d->brick_dy;
  const double        brick_dz = (double) d->brick_dz;

  int                 fieldid;

  /* set value for first field */
  val[0] = exp (M_PI * x);

  /* set values for remaining fields */
  for (fieldid = 1; fieldid < n_fields; fieldid++) {
    if (fieldid % 2 == 0) { /* if even field-id */
      val[fieldid] =   cos (((double) fieldid) * cx / brick_dx * M_PI * x)
                     * cos (((double) fieldid) * cy / brick_dy * M_PI * y)
                     * cos (((double) fieldid) * cz / brick_dz * M_PI * z);
    }
    else { /* if odd field-id */
      val[fieldid] =   sin (((double) fieldid) * cx / brick_dx * M_PI * x)
                     * sin (((double) fieldid) * cy / brick_dy * M_PI * y)
                     * sin (((double) fieldid) * cz / brick_dz * M_PI * z);
    }
  }
}

/**
 *
 */
static void
slabs_test_interpolation (MPI_Comm mpicomm,
                          const slabs_domain_shape_t domain_shape,
                          const int minlevel, const int maxlevel,
                          char *refine, const int order)
{
  const char         *this_fn_name = "slabs_test_interpolation";

  slabs_physics_options_t  physics_options;
  slabs_discr_options_t  discr_options;
  p8est_t            *p8est;
  mangll_t           *mangll;
  mangll_cnodes_t    *cnodes;
  ymir_mesh_t        *mesh;
  ymir_pressure_elem_t  *press_elem;
  slabs_stokes_state_t  *state;
  slabs_test_set_data_t  data;
  slabs_discr_amr_indicator_params_t  *indicator_params;
  ymir_vec_t         *fine_vel_press;
  ymir_dvec_t        *fine_norm_dual_rate;
  double              fine_norm;
  double              fine_vel_norm, fine_press_norm, fine_dual_norm;
  double              abs_error;
  double              abs_error_vel, abs_error_press, abs_error_dual;

  /* test is only implemented for the cube and brick domain */
  if (domain_shape != SL_DOMAIN_CUBE && domain_shape != SL_DOMAIN_BRICK) {
    return;
  }

  YMIR_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  /* set required options */
  physics_options.domain_shape = domain_shape;
  if (domain_shape == SL_DOMAIN_CUBE) {
    physics_options.domain_brick_dx = 1;
    physics_options.domain_brick_dy = 1;
    physics_options.domain_brick_dz = 1;
  }
  else if (domain_shape == SL_DOMAIN_BRICK) {
    physics_options.domain_brick_dx = SL_TEST_BRICK_DX;
    physics_options.domain_brick_dy = SL_TEST_BRICK_DY;
    physics_options.domain_brick_dz = SL_TEST_BRICK_DZ;
  }
  //physics_options.bc_type = SL_VEL_BC_DIRICHLET_ALL;
  physics_options.bc_type = SL_VEL_BC_DIRICHLET_NORM;
  physics_options.temperature_type = SL_TEMP_NONE;
  physics_options.weakzone_type = SL_WEAKZONE_NONE;
  physics_options.viscosity_type = SL_VISCOSITY_CONST;
  physics_options.viscosity_scaling = 1.0;
  physics_options.rhs_scaling = 1.0;
  physics_options.plume_type = SL_PLUME_INNER;
  physics_options.plume_center_x = 0.5 * physics_options.domain_brick_dx;
  physics_options.plume_center_y = 0.5 * physics_options.domain_brick_dy;
  physics_options.plume_center_z = 0.5 * physics_options.domain_brick_dz;
  if (domain_shape == SL_DOMAIN_CUBE) {
    physics_options.plume_decay = 10.0;
    physics_options.plume_scaling = -1.0;
  }
  else if (domain_shape == SL_DOMAIN_BRICK) {
    physics_options.plume_decay = SL_TEST_BRICK_INNER_PLUME_DECAY;
    physics_options.plume_scaling = SL_TEST_BRICK_INNER_PLUME_SCALING;
  }
  discr_options.minlevel = (int8_t) minlevel;
  discr_options.maxlevel = (int8_t) maxlevel;
  discr_options.X_fn = slabs_discr_identity_X;
  discr_options.order = order;
  discr_options.refine = refine;
  discr_options.amr_max_steps = 1;
  discr_options.amr_rel_threshold = 0.0;
  discr_options.mesh_partitioning_type = SL_MESH_PARTITIONING_ELEM;

  /* create p4est and boundary variables */
  p8est = slabs_discr_p8est_new (mpicomm, &physics_options, &discr_options);
  slabs_discr_options_set_boundary (&discr_options, p8est, &physics_options);

  /* create mesh */
  slabs_discr_mangll_and_cnodes_new (&mangll, &cnodes, p8est, &discr_options);
  slabs_discr_ymir_new (&mesh, &press_elem, mangll, cnodes, &discr_options);

  /* create Stokes state */
  state = slabs_stokes_state_new (p8est);
  slabs_stokes_state_init_temp (state, cnodes);
  slabs_stokes_state_init_temp_vec (state, mesh);
  ymir_cvec_set_function (state->temp_vec, slabs_physics_temperature_set_fn,
                          &physics_options);
  slabs_stokes_state_init_vel_press (state, mesh, press_elem);
  slabs_stokes_state_init_dual_tensor (state, mesh);

  /* set function data which is used to set the exact values */
  data.brick_dx = physics_options.domain_brick_dx;
  data.brick_dy = physics_options.domain_brick_dy;
  data.brick_dz = physics_options.domain_brick_dz;

  /* set exact values on coarse mesh */
  data.cx = 2.0;
  data.cy = 4.0;
  data.cz = 1.0;
  data.n_fields = state->vel_press_vec->ncfields;
  ymir_cvec_set_function (state->vel_press_vec, slabs_test_interp_set_vel_fn,
                          &data);
  data.n_fields = state->vel_press_vec->nefields;
  ymir_evec_set_function (state->vel_press_vec, slabs_test_interp_set_press_fn,
                          &data);
  data.cx = 1.0;
  data.cy = 1.0;
  data.cz = 1.0;
  data.n_fields = state->dual_vec->ndfields;
  ymir_dvec_set_function (state->dual_vec,
                          slabs_test_interp_set_dual_fn, &data);

  /* refine mesh uniformly and interpolate Stokes state */
  indicator_params = slabs_discr_amr_indicator_params_new (1);
  indicator_params->type[0] = SL_AMR_INDICATOR_REFINE_ALL;
  indicator_params->tol_min[0] = 0.0;
  indicator_params->tol_max[0] = 0.5;
  indicator_params->level_min[0] = 0;
  indicator_params->level_max[0] = 0;
  slabs_discr_amr (state, &mesh, &press_elem, p8est, indicator_params,
                   &physics_options, &discr_options, 0);
  slabs_discr_amr_indicator_params_destroy (indicator_params);

  /* set exact values on fine mesh */
  data.cx = 2.0;
  data.cy = 4.0;
  data.cz = 1.0;
  fine_vel_press = ymir_stokes_vec_new (mesh, press_elem);
  data.n_fields = fine_vel_press->ncfields;
  ymir_cvec_set_function (fine_vel_press, slabs_test_interp_set_vel_fn, &data);
  data.n_fields = fine_vel_press->nefields;
  ymir_evec_set_function (fine_vel_press, slabs_test_interp_set_press_fn,
                          &data);
  fine_norm_dual_rate = ymir_dvec_new (mesh, 6, YMIR_GAUSS_NODE);
  data.cx = 1.0;
  data.cy = 1.0;
  data.cz = 1.0;
  data.n_fields = fine_norm_dual_rate->ndfields;
  ymir_dvec_set_function (fine_norm_dual_rate,
                          slabs_test_interp_set_dual_fn, &data);

  /* compute norm of exact vector */
  fine_vel_norm = ymir_cvec_norm (fine_vel_press);
  fine_press_norm = ymir_evec_norm (fine_vel_press);
  fine_dual_norm = ymir_dvec_norm (fine_norm_dual_rate);
  fine_norm = sqrt (  SC_SQR (fine_vel_norm) + SC_SQR (fine_press_norm)
                    + SC_SQR (fine_dual_norm) );

  /* check error */
  ymir_vec_add (-1.0, state->vel_press_vec, fine_vel_press);
  abs_error_vel = ymir_cvec_norm (fine_vel_press);
  abs_error_press = ymir_evec_norm (fine_vel_press);
  ymir_vec_add (-1.0, state->dual_vec, fine_norm_dual_rate);
  abs_error_dual = ymir_dvec_norm (fine_norm_dual_rate);
  abs_error = sqrt (  SC_SQR (abs_error_vel) + SC_SQR (abs_error_press)
                    + SC_SQR (abs_error_dual) );

  /* print errors */
  if (0.0 < fine_vel_norm) {
    YMIR_GLOBAL_INFOF ("%s: velocity abs error %1.3e, rel error %1.3e\n",
                       this_fn_name, abs_error_vel,
                       abs_error_vel / fine_vel_norm);
  }
  else {
    YMIR_GLOBAL_INFOF ("%s: velocity abs error %1.3e\n", this_fn_name,
                       abs_error_vel);
  }
  if (0.0 < fine_press_norm) {
    YMIR_GLOBAL_INFOF ("%s: pressure abs error %1.3e, rel error %1.3e\n",
                       this_fn_name, abs_error_press,
                       abs_error_press / fine_press_norm);
  }
  else {
    YMIR_GLOBAL_INFOF ("%s: pressure abs error %1.3e\n", this_fn_name,
                       abs_error_press);
  }
  if (0.0 < fine_dual_norm) {
    YMIR_GLOBAL_INFOF ("%s: dual abs error %1.3e, rel error %1.3e\n",
                       this_fn_name, abs_error_dual,
                       abs_error_dual / fine_dual_norm);
  }
  else {
    YMIR_GLOBAL_INFOF ("%s: dual abs error %1.3e\n",
                       this_fn_name, abs_error_dual);
  }
  if (0.0 < fine_norm) {
    YMIR_GLOBAL_INFOF ("%s: total abs error %1.3e, rel error %1.3e\n",
                       this_fn_name, abs_error, abs_error / fine_norm);
  }
  else {
    YMIR_GLOBAL_INFOF ("%s: total abs error %1.3e\n", this_fn_name,
                       abs_error);
  }

  /* destroy vectors */
  ymir_vec_destroy (fine_vel_press);
  ymir_vec_destroy (fine_norm_dual_rate);
  slabs_stokes_state_destroy (state);

  /* destroy mesh, p4est, and options */
  slabs_discr_ymir_mangll_destroy (mesh, press_elem);
  slabs_discr_p8est_destroy (p8est);
  slabs_discr_options_clear_boundary (&discr_options);

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

