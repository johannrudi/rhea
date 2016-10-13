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

#ifndef SLABS_NORM_H
#define SLABS_NORM_H

#include <ymir_Hminus1_norm.h>
#include <ymir_stokes_op.h>
#include <ymir_stokes_pc.h>
#include <slabs_base.h>

/* enumerator for norms */
typedef enum
{
  SL_NORM_VEC_L2,
  SL_NORM_VEC_L2_MASS_WEIGHTED,
  SL_NORM_FNC_L2,
  SL_NORM_FNC_HMINUS1,   /* H^-1 norm for velocity */
  SL_NORM_FNC_HMINUS1_L2 /* mixed norm: H^-1 for velocity & L^2 for pressure */
}
slabs_norm_type_t;

/**
 *
 */
void
slabs_norm_innerprod_comp (double *innerprod_vel, double *innerprod_press,
                           ymir_vec_t *up_l, ymir_vec_t *up_r,
                           ymir_pressure_elem_t *press_elem,
                           slabs_norm_type_t norm_type,
                           ymir_Hminus1_norm_op_t *norm_op);

/**
 *
 */
double
slabs_norm (double *norm_vel, double *norm_press,
            ymir_vec_t *up, ymir_pressure_elem_t *press_elem,
            slabs_norm_type_t norm_type, ymir_Hminus1_norm_op_t *norm_op);

/**
 *
 */
void
slabs_norm_symtensor_innerprod_frobenius_elem (const sc_dmatrix_t *tensor_l,
                                               const sc_dmatrix_t *tensor_r,
                                               sc_dmatrix_t *innerprod);

/**
 *
 */
void
slabs_norm_symtensor_frobenius_elem (const sc_dmatrix_t *tensor,
                                     sc_dmatrix_t *norm);

/**
 *
 */
double
slabs_norm_symtensor_frobenius_min (const ymir_dvec_t *tensor);

/**
 *
 */
double
slabs_norm_symtensor_frobenius_max (const ymir_dvec_t *tensor);

/**
 *
 */
void
slabs_norm_symtensor_normalize_threshold_frobenius_elem (sc_dmatrix_t *tensor,
                                                         const double
                                                         threshold);
/**
 *
 */
void
slabs_norm_symtensor_normalize_frobenius_elem (sc_dmatrix_t *tensor);

/**
 *
 */
ymir_gloidx_t
slabs_norm_symtensor_normalize_threshold_frobenius (ymir_dvec_t *tensor,
                                                    const double threshold);

/**
 *
 */
ymir_gloidx_t
slabs_norm_symtensor_normalize_frobenius (ymir_dvec_t *tensor);

/**
 *
 */
void
slabs_norm_symtensor_scale_threshold_frobenius_elem (sc_dmatrix_t *tensor,
                                                     sc_dmatrix_t *scaling,
                                                     const double threshold);

/**
 *
 */
void
slabs_norm_symtensor_scale_frobenius_elem (sc_dmatrix_t *tensor,
                                           sc_dmatrix_t *scaling);

/**
 *
 */
ymir_gloidx_t
slabs_norm_symtensor_scale_threshold_frobenius (ymir_dvec_t *tensor,
                                                const ymir_dvec_t *scaling,
                                                const double threshold);

/**
 *
 */
ymir_gloidx_t
slabs_norm_symtensor_scale_frobenius (ymir_dvec_t *tensor,
                                      const ymir_dvec_t *scaling);

/**
 *
 */
double
slabs_norm_of_residual (double *norm_vel, double *norm_press,
                        ymir_vec_t *residual_up,
                        ymir_stokes_op_t *stokes_op,
                        ymir_stokes_pc_t *stokes_pc,
                        slabs_krylov_type_t krylov_type,
                        slabs_norm_type_t norm_type,
                        ymir_Hminus1_norm_op_t *norm_op);

/**
 *
 */
double
slabs_norm_of_residual_l2 (ymir_vec_t *residual_up,
                           ymir_stokes_op_t *stokes_op,
                           ymir_stokes_pc_t *stokes_pc,
                           slabs_krylov_type_t krylov_type);

/* given the current solution up and the current operator stokes_op, compute
 * the residual for the state stokes system, and potentially
 * return the norms of the residuals in the velocity and pressure
 * equations.  Note, this actually returns the negative residual, approrpriate
 * as a right hand side in a Newton step */
double
slabs_norm_compute_residual (ymir_vec_t *residual_up,
                             double *norm_vel, double *norm_press,
                             ymir_vec_t *up,
                             ymir_cvec_t *rhs_u_point,
                             ymir_stokes_op_t *stokes_op,
                             ymir_stokes_pc_t *stokes_pc,
                             slabs_krylov_type_t krylov_type,
                             slabs_norm_type_t norm_type,
                             ymir_Hminus1_norm_op_t *norm_op);

/**
 *
 */
double
slabs_norm_compute_residual_l2 (ymir_vec_t *residual_up,
                                ymir_vec_t *up,
                                ymir_cvec_t *rhs_u_point,
                                ymir_stokes_op_t *stokes_op,
                                ymir_stokes_pc_t *stokes_pc,
                                slabs_krylov_type_t krylov_type);

/**
 *
 */
double
slabs_norm_compute_residual_dual (ymir_dvec_t *residual_dual,
                                  ymir_vec_t *up,
                                  ymir_dvec_t *dual,
                                  ymir_stokes_op_t *stokes_op,
                                  const slabs_nl_solver_primaldual_type_t
                                    nl_solver_primaldual_type,
                                  const slabs_krylov_type_t krylov_type);

/**
 *
 */
void
slabs_norm_weight_residual (ymir_vec_t *residual_up,
                            ymir_pressure_elem_t *press_elem);

/******************************************************************************
 * Tests
 *****************************************************************************/

/**
 * Tests computation of Frobenius-norm of the strain rate tensor.
 */
void
slabs_norm_test_frobenius_of_strain_rate (ymir_mesh_t *mesh);

/**
 * Tests the Operator of the H^-1 norm by solving a Poisson problem and
 * comparing the error to the exact solution.
 */
void
slabs_norm_Hminus1_stiff_test (ymir_Hminus1_norm_op_t *norm_op,
                               ymir_mesh_t *mesh,
                               slabs_domain_shape_t domain_shape);

/**
 * Tests norms for mesh-independence.
 */
void
slabs_norm_test_mesh_independence (MPI_Comm mpicomm,
                                   const slabs_domain_shape_t domain_shape,
                                   const int minlevel, const int maxlevel,
                                   char *refine, const int order,
                                   double mass_scaling,
                                   const int max_refine_steps);

#endif /* SLABS_NORM_H */
