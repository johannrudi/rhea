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

#ifndef SLABS_BASE_H
#define SLABS_BASE_H

#include <p8est.h>
#include <p4est_to_p8est.h>
#include <ymir_mesh.h>
#include <ymir_pressure_elem.h>
#include <ymir_velocity_dirichlet.h>

/* enumerator for types of nodes */
typedef enum
{
  SL_GLL_CONTINUOUS_NODE,
  SL_GLL_DISCONTINUOUS_NODE,
  SL_GAUSS_NODE
}
slabs_node_type_t;

/* enumerator for domain shapes */
typedef enum
{
  SL_DOMAIN_CUBE,
  SL_DOMAIN_BRICK,
  SL_DOMAIN_SHELL,
  SL_DOMAIN_SHELL_CHUNK,
  SL_DOMAIN_SHELL_SLICE
}
slabs_domain_shape_t;

/* enumerator for velocity boundray conditions */
typedef enum
{
  /* Dirichlet all */
  SL_VEL_BC_DIRICHLET_ALL = 0,

  /* Dirchlet in normal direction */
  SL_VEL_BC_DIRICHLET_NORM = 1,

  /* Dirchlet norm (if shell: fix two points on inner sphere for rotation
   * invariance) */
  SL_VEL_BC_DIRICHLET_NORM_FIXDOF = 2,

  /* Dirichlet all on inner sphere and Dirichlet norm on outer sphere of a shell
   * (for shell domain only) */
  SL_VEL_BC_DIRICHLET_ALL_INNER = 3,

  /* Dirichlet all on side faces and Neumann on top and bottom faces
   * (for rectangular domains only) */
  SL_VEL_BC_DIR_SIDES_NEU_TB = 4,

  /* Dirichlet norm on side faces and Dirichlet all on top and bottom faces
   * (for rectangular domains only) */
  SL_VEL_BC_DIR_NORM_SIDES_DIR_ALL_TB = 5,

  /* Neumann */
  SL_VEL_BC_NEUMANN_ALL = 6
}
slabs_vel_bc_t;

/* enumerator for boundary faces */
typedef enum
{
  SL_NONE = -1,
  SL_BASE = 0,
  SL_TOP,
  SL_SIDE1,
  SL_SIDE2,
  SL_SIDE3,
  SL_SIDE4
}
slabs_boundary_face_t;

/* enumerator for type of nonlinear solver */
typedef enum
{
  SL_NL_SOLVER_NONE = -1,
  SL_NL_SOLVER_PICARD = 0,
  SL_NL_SOLVER_NEWTON,
  SL_NL_SOLVER_PICARD_NEWTON
}
slabs_nl_solver_type_t;

/* enumerator for type of primal-dual method */
typedef enum
{
  SL_NL_SOLVER_PRIMALDUAL_NONE,

  /* dual = grad_s (u) / ||grad_s (u)||_{F} */
  SL_NL_SOLVER_PRIMALDUAL_NORMSTRAIN,

  /* dual = 2 * viscosity * grad_s (u) */
  SL_NL_SOLVER_PRIMALDUAL_STRESS
}
slabs_nl_solver_primaldual_type_t;

/* enumerator for type of scaling for primal-dual method */
typedef enum
{
  SL_NL_SOLVER_PRIMALDUAL_SCAL_NONE,
  SL_NL_SOLVER_PRIMALDUAL_SCAL_NORMALIZE, /* normalize by Frobenius-norm */
  SL_NL_SOLVER_PRIMALDUAL_SCAL_THRESHOLD  /* normalize by Frob.-norm if > 1 */
}
slabs_nl_solver_primaldual_scal_type_t;

/* enumerator for type of Krylov solver */
typedef enum
{
  SL_KRYLOV_NONE = -1,
  SL_KRYLOV_MINRES = 0,
  SL_KRYLOV_GMRES = 1
}
slabs_krylov_type_t;

/**
 * Assigns boundary faces of enumerator `slabs_boundary_face_t` to tree faces.
 */
ymir_topidx_t *
slabs_cube_tree_to_bf (p8est_connectivity_t *conn);

ymir_topidx_t *
slabs_shell_tree_to_bf (p8est_connectivity_t *conn);

/**
 * Sets the type of Dirichlet boundary condition for boundary faces.
 */
ymir_dir_code_t
slabs_cube_vel_dir_fn (double X, double Y, double Z,
                       double nx, double ny, double nz,
                       ymir_topidx_t face, ymir_locidx_t node_id,
                       void *dir_data);

ymir_dir_code_t
slabs_shell_vel_dir_fn (double X, double Y, double Z,
                        double nx, double ny, double nz,
                        ymir_topidx_t face, ymir_locidx_t node_id,
                        void *dir_data);

/**
 * Restricts degrees of freedom by setting further Dirichlet boundary
 * conditions.
 */
void
slabs_shell_vel_dir_restrict_dof (ymir_vel_dir_t * vel_dir);

/**
 * Gets maximum absolute value of all entries of a matrix.
 */
double
slabs_matrix_compute_abs_max (const sc_dmatrix_t *mat);

/**
 * Computes minimum and maximum absolute value of all elements of a matrix.
 */
void
slabs_matrix_compute_abs_min_max (double *min, double *max,
                                  const sc_dmatrix_t *mat);

/**
 * Bounds all entries of a matrix to a minimum and a maximum.
 */
void
slabs_matrix_bound_values (sc_dmatrix_t *mat, const double min,
                           const double max);

/**
 * Rounds all entries of a matrix.
 */
void
slabs_matrix_round (sc_dmatrix_t *mat);

/**
 * Applies floor to all entries of a matrix.
 */
void
slabs_matrix_floor (sc_dmatrix_t *mat);

/**
 * Applies ceil to all entries of a matrix.
 */
void
slabs_matrix_ceil (sc_dmatrix_t *mat);

/**
 * Perform element-wise multiplication, Y[i,:] := Y[i,:] .* X[i,0], for all i.
 */
void
slabs_matrix_multiply_in_1d (const sc_dmatrix_t *X, sc_dmatrix_t *Y);

/**
 * Perform element-wise division, Y[i,:] := Y[i,:] ./ X[i,0], for all i.
 */
void
slabs_matrix_divide_in_1d (const sc_dmatrix_t *X, sc_dmatrix_t *Y);

/**
 *
 */
double
slabs_elem_volume (mangll_t *mangll, mangll_locidx_t elid);

/**
 *
 */
void
slabs_elem_get_gauss_coordinates (double *x, double *y, double *z,
                                  const mangll_locidx_t elid,
                                  mangll_t *mangll, double *tmp_el);

/**
 *
 */
void
slabs_interp_gll_to_gauss (sc_dmatrix_t *mat_gll, sc_dmatrix_t *mat_gauss,
                           mangll_t *mangll);

/**
 *
 */
void
slabs_interp_gauss_to_gll (sc_dmatrix_t *mat_gauss, sc_dmatrix_t *mat_gll,
                           mangll_t *mangll);

/**
 *
 */
void
slabs_gradient_gll_to_gll_elem (sc_dmatrix_t *u, sc_dmatrix_t *grad_u,
                                mangll_t *mangll, const mangll_locidx_t elid,
                                sc_dmatrix_t *tmp_du);

/**
 *
 */
void
slabs_gradient_gll_to_gauss_elem (sc_dmatrix_t *u, sc_dmatrix_t *grad_u,
                                  mangll_t *mangll, mangll_locidx_t elid,
                                  sc_dmatrix_t *tmp_du, sc_dmatrix_t *tmp);

/**
 *
 */
void
slabs_gradient_gll_to_gauss (sc_dmatrix_t *u, sc_dmatrix_t *grad_u,
                             mangll_t *mangll);

/**
 *
 */
void
slabs_second_invariant_elem (sc_dmatrix_t *u, sc_dmatrix_t *IIe,
                             mangll_t *mangll, mangll_locidx_t elid,
                             sc_dmatrix_t *tmp_grad_u,
                             sc_dmatrix_t *tmp_du,
                             sc_dmatrix_t *tmp,
                             slabs_node_type_t node_type);

/**
 *
 */
void
slabs_second_invariant (sc_dmatrix_t *u, sc_dmatrix_t *IIe,
                        mangll_t *mangll, slabs_node_type_t node_type);

/**
 *
 */
void
slabs_strain_rate_tensor_elem (sc_dmatrix_t *strain_rate_tensor,
                               mangll_t *mangll,
                               sc_dmatrix_t *grad_u);

/**
 *
 */
int
slabs_strain_rate_tensor_idx (const int row, const int col);

/**
 *
 */
int
slabs_vec_is_finite (ymir_vec_t *vec);

/**
 *
 */
void
slabs_cvec_set_min_value (ymir_cvec_t *cvec, const double min);

/**
 *
 */
void
slabs_cvec_bound_values (ymir_vec_t *cvec, const double min,
                         const double max);

void
slabs_dvec_bound_values (ymir_vec_t *dvec, const double min,
                         const double max);

/**
 *
 */
void
slabs_cvec_round (ymir_vec_t *cvec);

void
slabs_dvec_round (ymir_vec_t *dvec);

/**
 *
 */
void
slabs_cvec_floor (ymir_vec_t *cvec);

void
slabs_dvec_floor (ymir_vec_t *dvec);

/**
 *
 */
void
slabs_cvec_ceil (ymir_vec_t *cvec);

void
slabs_dvec_ceil (ymir_vec_t *dvec);

/**
 *
 */
void
slabs_cvec_compute_magnitude (const ymir_vec_t *cvec, ymir_vec_t *magn);

void
slabs_dvec_compute_magnitude (const ymir_vec_t *dvec, ymir_vec_t *magn);

/**
 *
 */
ymir_dvec_t *
slabs_dvec_new_from_element_data (sc_dmatrix_t *data,
                                  ymir_mesh_t *mesh,
                                  ymir_node_type_t node_type);

/**
 *
 */
void
slabs_stokes_vec_get_components_view (ymir_cvec_t **u, ymir_evec_t **p,
                                      ymir_vec_t *up);

/**
 *
 */
void
slabs_stokes_vec_get_components (ymir_vec_t **u, ymir_vec_t **p,
                                 ymir_vec_t *up,
                                 ymir_pressure_elem_t *press_elem);

/**
 *
 */
void
slabs_stokes_vec_get_components_update (ymir_vec_t *u, ymir_vec_t *p,
                                        ymir_vec_t *up,
                                        ymir_pressure_elem_t *press_elem);

/**
 *
 */
void
slabs_stokes_vec_set_components (ymir_vec_t *up,
                                 ymir_vec_t *u, ymir_vec_t *p,
                                 ymir_pressure_elem_t *press_elem);

/******************************************************************************
 * Tests
 *****************************************************************************/

/**
 *
 */
void
slabs_base_test_second_invariant (ymir_mesh_t *mesh);

/**
 *
 */
void
slabs_base_test_strain_rate_tensor (ymir_mesh_t *mesh);

#endif /* SLABS_BASE_H */
