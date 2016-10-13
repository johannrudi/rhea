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

#ifndef SLABS_PHYSICS_EXTENDED_H
#define SLABS_PHYSICS_EXTENDED_H

#include <slabs_physics.h>

/* global variables */
extern WeakzoneKDTree *slabs_physics_weakzone_kdtree;

/**
 * Computes the radius of a shell domain or the corresponding value for a
 * rectangular domain.
 */
double
slabs_compute_radius (const double x, const double y, const double z,
                      slabs_physics_options_t *physics_options);

/**
 * Computes the longitude of a shell domain or the corresponding value for a
 * rectangular domain.
 */
double
slabs_compute_longitude (const double x, const double y, const double z,
                         slabs_physics_options_t *physics_options);

/**
 * Returns whether the element with coordinates (x,y,z) is located in the
 * upper mantle.
 */
int
slabs_physics_elem_in_upper_mantle (const double *x, const double *y,
                                    const double *z, const int *Vmask,
                                    slabs_physics_options_t *physics_options);

/**
 * Returns whether the element with coordinates (x,y,z) is located close enough
 * to weak zones.
 */
int
slabs_physics_elem_has_weakzone (const double *x, const double *y,
                                 const double *z, const int *Vmask,
                                 slabs_physics_options_t *physics_options);

/**
 * Computes the background temperature (e.g. for right-hand side forcing)
 * according to physics options as `data`.
 * (Callback function for `ymir_vec_set_function`)
 *
 * \param [out] back_temp    Background temperature
 * \param [in] x             x-coordinate
 * \param [in] y             y-coordinate
 * \param [in] z             z-coordinate
 * \param [in] nid           Local node index
 * \param [in] data          physics_options
 */
void
slabs_physics_background_temp_set_fn (double *back_temp, double x, double y,
                                      double z, ymir_locidx_t nid, void *data);

/**
 * Computes distance to weak zone between plates.
 */
double
slabs_2plates_weakzone_poly2_subdu_dist (double r, double lon,
                                         double start_node,
                                         double start_val,
                                         double start_deriv,
                                         double end_node,
                                         double end_val);

/**
 * Computes distance to weak zone at ridge.
 */
double
slabs_2plates_weakzone_poly2_ridge_dist (double r, double lon,
                                         double end_node,
                                         double end_val);

/**
 *
 */
void
slabs_physics_weak_dist_set_fn (double *dist, double x, double y, double z,
                                ymir_locidx_t nid, void *data);

/**
 *
 */
void
slabs_physics_weak_dist_set_face_fn (double *dist,
                                     double x, double y, double z,
                                     double nx, double ny, double nz,
                                     ymir_topidx_t face, ymir_locidx_t node_id,
                                     void *data);

/**
 * Processes temperature before it can be used for computations.
 */
void
slabs_temp_postprocess_elem (sc_dmatrix_t *temp_el_mat,
                             const double *x, const double *y, const double *z,
                             slabs_physics_options_t *physics_options);

/**
 * Computes weak zone factor of an element.
 */
void
slabs_weak_elem (sc_dmatrix_t *weak_el_mat,
                 const double *x, const double *y, const double *z,
                 const int *Vmask, slabs_physics_options_t *physics_options);

/**
 * Calculates viscosity for temperature dependent (only) viscosity law.
 *
 *   visc (T) = exp (E * (0.5 - T))
 */
double
slabs_visc_temp_fn (const double visc_temp_decay, const double temp);

/**
 * Computes temperature dependent viscosity in an element.
 */
void
slabs_visc_temp_elem (sc_dmatrix_t *visc_el_mat,
                      const double *x, const double *y, const double *z,
                      const int *Vmask,
                      const sc_dmatrix_t *temp_el_mat,
                      const sc_dmatrix_t *weak_el_mat,
                      slabs_physics_options_t *physics_options,
                      const int restrict_to_bounds);

/**
 * Computes nonlinear viscosity in an element.
 */
void
slabs_visc_nl_elem (sc_dmatrix_t *visc_el_mat,
                    sc_dmatrix_t *dvisc_dIIe_el_mat,
                    sc_dmatrix_t *rank1_scal_el_mat,
                    sc_dmatrix_t *bounds_el_mat,
                    sc_dmatrix_t *yielding_el_mat,
                    const double *x, const double *y, const double *z,
                    const int *Vmask,
                    const sc_dmatrix_t *temp_el_mat,
                    const sc_dmatrix_t *weak_el_mat,
                    const sc_dmatrix_t *IIe_el_mat,
                    slabs_physics_options_t *physics_options,
                    const int prescribe_bounds_yielding);

/**
 * Computes viscosity for linear or nonlinear rheology.
 *
 * \param [out] viscosity         Viscosity dnodes vector
 * \param [out] dvisc_dIIe        Derviative of viscosity w.r.t. 2nd invariant
 * \param [out] rank1_tensor_scal Scaling for rank 1 component of coefficient
 * \param [out] bounds_marker     Viscosity bounds marker
 * \param [out] yielding_marker   Yielding viscosity marker
 * \param [in] state              Stokes state
 * \param [in] press_elem         Pressure element information
 * \param [in] vel_dir            Velocity Dirichlet BC information
 * \param [in] physics_options    Parameter list for mantle flow physics
 */
void
slabs_physics_compute_viscosity (ymir_vec_t *viscosity,
                                 ymir_vec_t *dvisc_dIIe,
                                 ymir_vec_t *rank1_tensor_scal,
                                 ymir_vec_t *bounds_marker,
                                 ymir_vec_t *yielding_marker,
                                 slabs_stokes_state_t *state,
                                 ymir_pressure_elem_t *press_elem,
                                 ymir_vel_dir_t *vel_dir,
                                 slabs_physics_options_t *physics_options,
                                 const int prescribe_bounds_yielding);

/**
 * Computes viscosity for linear Stokes operator.
 *
 * \param [out] viscosity       Viscosity dnodes vector
 * \param [in] state            Stokes state
 * \param [in] physics_options  Parameter list for mantle flow physics
 */
void
slabs_physics_compute_linear_viscosity (ymir_dvec_t *viscosity,
                                        slabs_stokes_state_t *state,
                                        slabs_physics_options_t
                                          *physics_options);

/**
 * Computes linear viscosity for nonlinear Stokes operator at initial iteration.
 *
 * \param [out] viscosity         Viscosity dnodes vector
 * \param [out] dvisc_dIIe        Derviative of viscosity w.r.t. 2nd invariant
 * \param [out] rank1_tensor_scal Scaling for rank 1 component of coefficient
 * \param [out] bounds_marker     Viscosity bounds marker
 * \param [out] yielding_marker   Yielding viscosity marker
 * \param [in] state              Stokes state
 * \param [in] press_elem         Pressure element information
 * \param [in] vel_dir            Velocity Dirichlet BC information
 * \param [in] physics_options    Parameter list for mantle flow physics
 */
void
slabs_physics_compute_init_nl_viscosity (ymir_vec_t *viscosity,
                                         ymir_vec_t *dvisc_dIIe,
                                         ymir_vec_t *rank1_tensor_scal,
                                         ymir_vec_t *bounds_marker,
                                         ymir_vec_t *yielding_marker,
                                         slabs_stokes_state_t *state,
                                         ymir_pressure_elem_t *press_elem,
                                         ymir_vel_dir_t *vel_dir,
                                         slabs_physics_options_t
                                           *physics_options);

/**
 * Computes right-hand side in an element.
 */
void
slabs_rhs_elem (sc_dmatrix_t *rhs_el_mat,
                const double *x, const double *y, const double *z,
                const sc_dmatrix_t *temp_el_mat,
                slabs_physics_options_t *physics_options);

#endif /* SLABS_PHYSICS_EXTENDED_H */
