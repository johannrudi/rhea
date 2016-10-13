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

#ifndef SLABS_PHYSICS_H
#define SLABS_PHYSICS_H

#include <slabs_base.h>
#include <slabs_stokes_state.h>
#include <weakzone_wrapper.h>
#include <ymir_stokes_op.h>

/* basic constants */
#define SL_SEC_PER_YEAR (31557600.0)    /* seconds in a year (365.25*24*3600) */
#define SL_EARTH_RADIUS (6371.0e3)      /* mean radius of the Earth [m]       */
#define SL_UPPER_MANTLE_DEPTH (660.0e3) /* approx. depth of upper mantle [m]  */
#define SL_THERM_DIFFUS (1.0e-6)        /* thermal diffusivity [m^2 / s]      */
#define SL_TEMP_DIFF (1400.0)           /* temperature difference [K]         */
#define SL_VISC_REP (1.0e20)            /* representative viscosity [Pa s]    */

/* shell parameters */
#define SL_SHELL_RADIUS_BOTTOM 0.55
#define SL_SHELL_RADIUS_TOP 1.0

/* temperature parameters */
#define SL_DEFAULT_CONST_TEMP 0.5

/* set factor to estimate viscosity in the upper mantle from the lower mantle
 * for computing the viscosity for the initial nonlinear iteration;
 * LM is supposed to have a viscosity 10-100 times higher than UM */
#define SL_VISCOSITY_UM_REL_TO_LM 1.0 /* 0.1 */

/* enumerator for types of temperature */
typedef enum
{
  SL_TEMP_NONE,
  SL_TEMP_HOT_CORE,
  SL_TEMP_COLD_PLATE,
  SL_TEMP_HOT_CORE_COLD_PLATE,
  SL_TEMP_2PLATES_POLY2,
  SL_TEMP_2PLATES_RHEA1,
  SL_TEMP_IMPORT_FILE
}
slabs_temp_t;

/* enumerator for types of weak fault zones */
typedef enum
{
  SL_WEAKZONE_NONE,
  SL_WEAKZONE_2PLATES_POLY2,
  SL_WEAKZONE_2PLATES_RHEA1,
  SL_WEAKZONE_IMPORT_FILE
}
slabs_weakzone_t;

/* enumerator for types of weak zone coarsening */
typedef enum
{
  SL_WEAKZONE_COARSEN_EVAL,
  SL_WEAKZONE_COARSEN_EVAL_PIECEWISE_CONST,
  SL_WEAKZONE_COARSEN_LIN_INTERP
}
slabs_weakzone_coarsen_t;

/* enumerator for types of smoothing of weak fault zones */
typedef enum
{
  SL_WEAKZONE_SMOOTHING_CONTINUOUS,
  SL_WEAKZONE_SMOOTHING_DISCONTINUOUS
}
slabs_weakzone_smoothing_t;

/* enumerator for types of viscosities */
typedef enum
{
  SL_VISCOSITY_CONST,
  SL_VISCOSITY_LINEAR,
  SL_VISCOSITY_NONLINEAR
}
slabs_viscosity_t;

/* enumerator for types of initial viscosities for nonlinear Stokes problem*/
typedef enum
{
  SL_VISCOSITY_INIT_NL_STOKES_DEFAULT,
  SL_VISCOSITY_INIT_NL_STOKES_CONST_CENTERED_TO_BOUNDS,
  SL_VISCOSITY_INIT_NL_STOKES_CONST_MIN_BOUND,
  SL_VISCOSITY_INIT_NL_STOKES_TEMP,
  SL_VISCOSITY_INIT_NL_STOKES_TEMP_UM_REL_TO_LM
}
slabs_viscosity_init_nl_stokes_t;

/* enumerator for types of viscosity models */
typedef enum
{
  /* (1) weak zone, (2) yielding, (3) upper bound, (4) lower bound
   * viscosity bounds via cut-off */
  SL_VISCOSITY_MODEL_WYUL,

  /* (1) upper bound for temp. dep. viscosity, (2) weak zone, (3) yielding,
   * (4) upper bound for full viscosity, (5) lower bound
   * viscosity bounds via cut-off */
  SL_VISCOSITY_MODEL_UWYUL,

  /* (1) upper bound, (2) weak zone, (3) yielding, (4) lower bound
   * viscosity bounds via cut-off */
  SL_VISCOSITY_MODEL_UWYL,

  /* (1) upper bound, (2) weak zone, (3) yielding, (4) lower bound
   * upper viscosity bound via cut-off, lower viscosity bound via addition */
  SL_VISCOSITY_MODEL_UWYL_LREG,

  /* (1) upper bound + shift, (2) weak zone, (3) yielding, (4) lower bound
   * upper viscosity bound via shift, lower viscosity bound via addition */
  SL_VISCOSITY_MODEL_UWYL_SHIFT_LREG,

  /* (1) upper bound, (2) yielding, (3) weak zone, (4) lower bound
   * viscosity bounds via cut-off */
  SL_VISCOSITY_MODEL_UYWL,

  /* (1) upper bound, (2) yielding, (3) weak zone, (4) lower bound
   * viscosity bounds via right-shift */
  SL_VISCOSITY_MODEL_UYWL_SHIFT,

  /* (1) upper bound, (2) weak zone, (3) lower bound
   * adds regularization to 2nd invariant */
  SL_VISCOSITY_MODEL_UWL_IIE_REG
}
slabs_viscosity_model_t;

/* enumerator for types of viscosity coarsening */
typedef enum
{
  SL_VISCOSITY_COARSEN_EVAL,
  SL_VISCOSITY_COARSEN_EVAL_PIECEWISE_DISCONT
}
slabs_viscosity_coarsen_t;

/* enumerator for types of plumes */
typedef enum
{
  SL_PLUME_NONE,
  SL_PLUME_INNER,
  SL_PLUME_RISING
}
slabs_plume_t;

/* parameter list for mantle flow physics */
typedef struct slabs_physics_options
{
  /* domain parameters */
  slabs_domain_shape_t  domain_shape;
  int                 domain_brick_dx;
  int                 domain_brick_dy;
  int                 domain_brick_dz;

  char               *p4est_import_filename;

  /* domain properties */
  double              domain_x_min;
  double              domain_x_max;
  double              domain_y_min;
  double              domain_y_max;
  double              domain_z_min;
  double              domain_z_max;
  double              domain_lon_min;
  double              domain_lon_max;
  double              domain_radius_min;
  double              domain_radius_max;
  double              domain_volume;
  double              domain_center[3];
  double              domain_moment_of_inertia[3];

  /* boundary conditions */
  slabs_vel_bc_t      bc_type;
  int                 bc_default_dirichlet_scale;

  /* temperature parameters */
  slabs_temp_t        temperature_type;

  double              temp_background_plate_age;

  double              temp_import_plate_age_min;
  char               *temp_import_filename_txt;
  char               *temp_import_filename_bin;
  char               *temp_import_verification_out;

  double              temp_2plates_trench_longitude;
  double              temp_2plates_dip_angle;
  double              temp_2plates_subd_depth;
  double              temp_2plates_subd_width;
  double              temp_2plates_subd_edge_width;
  double              temp_2plates_subd_edge_smoothwidth;
  double              temp_2plates_subd_plate_velocity;
  double              temp_2plates_subd_plate_initial_age;
  double              temp_2plates_over_plate_age;

  /* weak zone parameters */
  slabs_weakzone_t    weakzone_type;

  char               *weakzone_import_filename_txt;
  char               *weakzone_import_filename_bin;
  char               *weakzone_import_verification_out;
  int                 weakzone_import_n_points;
  double              weakzone_import_thickness;
  double              weakzone_import_thickness_const;
  double              weakzone_import_weak_factor;

  double              weakzone_2plates_subdu_longitude;  /* 2plates weak zone */
  double              weakzone_2plates_subdu_dip_angle;  /* at subduction */
  double              weakzone_2plates_subdu_depth;
  double              weakzone_2plates_subdu_width;
  double              weakzone_2plates_subdu_thickness;
  double              weakzone_2plates_subdu_thickness_const;
  double              weakzone_2plates_subdu_weak_factor;
  double              weakzone_2plates_ridge_depth;      /* 2plates weak zone */
  double              weakzone_2plates_ridge_width;      /* at ridge */
  double              weakzone_2plates_ridge_smoothwidth;
  double              weakzone_2plates_ridge_weak_factor;

  /* velocity & pressure parameters */
  char               *velocity_import_filename;
  char               *pressure_import_filename;

  /* viscosity parameters */
  slabs_viscosity_t                 viscosity_type;
  slabs_viscosity_init_nl_stokes_t  viscosity_type_for_init_nl_stokes;
  slabs_viscosity_model_t           viscosity_model_type;
  double              viscosity_IIe_regularization;
  double              viscosity_min;
  double              viscosity_max;
  double              viscosity_temp_max;
  double              viscosity_scaling;
  double              viscosity_temp_decay;
  double              viscosity_stress_exponent;
  double              viscosity_stress_yield;
  double              viscosity_yielding_reg;
  double              viscosity_upper_mantle_radius;
  double              viscosity_lower_mantle_scaling;
  double              viscosity_lower_mantle_temp_decay;
  double              viscosity_lower_upper_transition_zone;

  /* viscosity coarsening */
  int                 viscosity_coarsen_eval;
  int                 viscosity_p_coarsen_eval;
  slabs_weakzone_coarsen_t  weakzone_coarsen_type;
  slabs_viscosity_coarsen_t viscosity_coarsen_type;
  double              viscosity_lower_upper_transition_zone_incr;
  ymir_vec_t         *current_fine_bounds;
  ymir_vec_t         *current_fine_yielding;

  /* right-hand side parameters */
  double              rhs_scaling;
  int                 rhs_random;
  int                 rhs_multiply_in_weak_zone;

  slabs_plume_t       plume_type;
  double              plume_center_x;
  double              plume_center_y;
  double              plume_center_z;
  double              plume_decay;
  double              plume_scaling;
}
slabs_physics_options_t;

/**
 * Computes and stores boundaries of a domain.
 *
 * \param [in/out] physics_options  Parameter list for mantle flow physics
 */
void
slabs_physics_compute_domain_bounds (slabs_physics_options_t *physics_options);

/**
 * Computes and stores the volume of a domain.
 *
 * \param [in/out] physics_options  Parameter list for mantle flow physics
 */
void
slabs_physics_compute_domain_volume (slabs_physics_options_t *physics_options);

/**
 * Computes and stores center of domain.
 *
 * \param [in/out] physics_options  Parameter list for mantle flow physics
 */
void
slabs_physics_compute_domain_center (slabs_physics_options_t *physics_options);

/**
 * Computes and stores the moment of inertia of a domain.
 * It is assumed that the density is uniformly one, hence the mass of the
 * domain equals the volume.
 *
 * \param [in/out] physics_options  Parameter list for mantle flow physics
 */
void
slabs_physics_compute_domain_moment_of_inertia (slabs_physics_options_t
                                                 *physics_options);

/**
 * Sets up domain variables.
 *
 * \param [in/out] physics_options  Parameter list for mantle flow physics
 */
void
slabs_physics_compute_domain_setup (slabs_physics_options_t *physics_options);

/**
 * Converts Cartesian coordinates (x, y, z) to spherical coordinates
 * (r, theta, phi), where the mathematical convention is used for the
 * spherical coordinates (see below).
 *
 * \param [in] x                     x-coordinate
 * \param [in] y                     y-coordinate
 * \param [in] z                     z-coordinate
 * \param [out] r = ||[x,y,z]||_2    radius
 * \param [out] theta = arctan(y/x)  azimutal angle
 * \param [out] phi = arccos(z/r)    polar angle
 */
void
slabs_convert_cartesian_to_spherical_math_conv (const double x,
                                                const double y,
                                                const double z,
                                                double *r,
                                                double *theta,
                                                double *phi);

/**
 * Converts Cartesian coordinates (x, y, z) to spherical coordinates
 * (r, phi, theta), where the geophysical convention is used for the
 * spherical coordinates (see below).
 *
 * \param [in] x                     x-coordinate
 * \param [in] y                     y-coordinate
 * \param [in] z                     z-coordinate
 * \param [out] r = ||[x,y,z]||_2    radius
 * \param [out] phi = arctan(y/x)    azimutal angle
 * \param [out] theta = arccos(z/r)  colatitude
 */
void
slabs_convert_cartesian_to_spherical_geo_conv (const double x,
                                               const double y,
                                               const double z,
                                               double *r,
                                               double *phi,
                                               double *theta);

/**
 *
 */
void
slabs_physics_postprocess_temperature (slabs_stokes_state_t *state,
                                       slabs_physics_options_t
                                         *physics_options);

/**
 *
 */
void
slabs_physics_verify_temperature (ymir_vec_t *temp_vec,
                                  slabs_physics_options_t *physics_options);

/**
 * Computes the temperature according to physics options as `data`.
 * (Callback function for `ymir_vec_set_function`)
 *
 * \param [out] temp         Temperature
 * \param [in] x             x-coordinate
 * \param [in] y             y-coordinate
 * \param [in] z             z-coordinate
 * \param [in] nid           Local node index
 * \param [in] data          physics_options
 */
void
slabs_physics_temperature_set_fn (double *temp, double x, double y, double z,
                                  ymir_locidx_t nid, void *data);

/**
 * Computes the weak zone parameters longitude, dip angle, and weak zone width
 * from given temperature and weak zone data for the `2plates_poly2` model
 * problem.
 *
 * \param [in/out] physics_options  Parameter list for mantle flow physics
 */
void
slabs_2plates_poly2_set_weakzone_params_from_temp (slabs_physics_options_t
                                                     *physics_options);

/**
 *
 */
void
slabs_physics_init_weakzone (MPI_Comm mpicomm,
                             slabs_physics_options_t *physics_options);

/**
 *
 */
void
slabs_physics_clear_weakzone (slabs_physics_options_t *physics_options);

/**
 * Computes the weak zone factor according to physics options as `data`.
 * (Callback function of type slabs_stokes_state_set_weakzone_t)
 *
 * \param [out] weakzone        Weak zone factor
 * \param [in] data             physics_options_t, parameter list for mantle
 *                              flow physics
 */
void
slabs_physics_compute_weakzone (ymir_dvec_t *weakzone, void *data);

/**
 * Computes coefficient for linear or nonlinear Stokes operator.
 *
 * \param [out] coeff             Coefficient vector at Gauss nodes
 * \param [out] coeff_deriv       Derviative of coefficient vector
 * \param [out] rank1_tensor_scal Scaling for rank 1 component of coefficient
 * \param [out] bounds_marker     Viscosity bounds marker
 * \param [out] yielding_marker   Yielding viscosity marker
 * \param [in] state              Stokes state
 * \param [in] press_elem         Pressure element information
 * \param [in] vel_dir            Velocity Dirichlet BC information
 * \param [in] physics_options    Parameter list for mantle flow physics
 */
void
slabs_physics_compute_stokes_coeff (ymir_vec_t *coeff,
                                    ymir_vec_t *coeff_deriv,
                                    ymir_vec_t *rank1_tensor_scal,
                                    ymir_vec_t *bounds_marker,
                                    ymir_vec_t *yielding_marker,
                                    slabs_stokes_state_t *state,
                                    ymir_pressure_elem_t *press_elem,
                                    ymir_vel_dir_t *vel_dir,
                                    slabs_physics_options_t *physics_options);

/**
 * Computes coefficient for linear Stokes operator.
 *
 * \param [out] coeff           Coefficient vector at Gauss nodes
 * \param [in] state            Stokes state
 * \param [in] physics_options  Parameter list for mantle flow physics
 */
void
slabs_physics_compute_lin_stokes_coeff (ymir_vec_t *coeff,
                                        slabs_stokes_state_t *state,
                                        slabs_physics_options_t
                                          *physics_options);

/**
 * Computes coefficient for nonlinear Stokes operator at initial nonlinear step.
 *
 * \param [out] coeff             Coefficient vector at Gauss nodes
 * \param [out] coeff_deriv       Derviative of coefficient vector
 * \param [out] rank1_tensor_scal Scaling for rank 1 component of coefficient
 * \param [out] bounds_marker     Viscosity bounds marker
 * \param [out] yielding_marker   Yielding viscosity marker
 * \param [in] state              Stokes state
 * \param [in] press_elem         Pressure element information
 * \param [in] vel_dir            Velocity Dirichlet BC information
 * \param [in] physics_options    Parameter list for mantle flow physics
 */
void
slabs_physics_compute_init_nl_stokes_coeff (ymir_vec_t *coeff,
                                            ymir_vec_t *coeff_deriv,
                                            ymir_vec_t *rank1_tensor_scal,
                                            ymir_vec_t *bounds_marker,
                                            ymir_vec_t *yielding_marker,
                                            slabs_stokes_state_t *state,
                                            ymir_pressure_elem_t *press_elem,
                                            ymir_vel_dir_t *vel_dir,
                                            slabs_physics_options_t
                                              *physics_options);

/* data for coarsening the Stokes coefficient */
typedef struct slabs_physics_coarsen_stokes_coeff_data
{
  slabs_stokes_state_t    *state;
  slabs_physics_options_t *physics_options;
  int                 init_nl_stokes_coeff;
  int                 coarsen_count;

  ymir_vec_t         *buffer_fine_weak;
  ymir_vec_t         *buffer_fine_temp;
  ymir_vec_t         *buffer_fine_vel;

  ymir_vec_t         *buffer_fine_bounds;
  ymir_vec_t         *buffer_fine_yielding;
}
slabs_physics_coarsen_stokes_coeff_data_t;

/**
 *
 */
slabs_physics_coarsen_stokes_coeff_data_t *
slabs_physics_coarsen_stokes_coeff_data_new (slabs_stokes_state_t *state,
                                             slabs_physics_options_t
                                               *physics_options);

/**
 *
 */
void
slabs_physics_coarsen_stokes_coeff_data_reset (void *data);

/**
 *
 */
void
slabs_physics_coarsen_stokes_coeff_data_destroy (
                               slabs_physics_coarsen_stokes_coeff_data_t *data);

/**
 * Computes coefficient for linear or nonlinear Stokes operator.
 */
void
slabs_physics_coarsen_stokes_coeff (ymir_vec_t *coarse_coeff,
                                    ymir_vec_t *coarse_coeff_deriv,
                                    ymir_vec_t *coarse_vel,
                                    mangll_t * upstream_interp_mangll,
                                    mangll_t * upstream_partition_mangll,
                                    ymir_mesh_t *coarse_ymir_mesh,
                                    ymir_pressure_elem_t *coarse_press_elem,
                                    ymir_vel_dir_t *coarse_vel_dir,
                                    p8est_t *coarse_p8est,
                                    const int is_fine_level,
                                    void *data);

/**
 * Computes right-hand side in (primal) function space.
 *
 * \param [out] rhs_u_point     Right-hand side vector
 * \param [in] state            Stokes state
 * \param [in] physics_options  Parameter list for mantle flow physics
 */
void
slabs_physics_compute_rhs_u_point (ymir_vec_t *rhs_u_point,
                                   slabs_stokes_state_t *state,
                                   slabs_physics_options_t *physics_options);

/**
 *
 */
void
slabs_physics_compute_normal_boundary_stress (ymir_vec_t *stress_bndr_norm,
                                              ymir_vec_t *up,
                                              ymir_vec_t *rhs_u_point,
                                              ymir_stokes_op_t *stokes_op);

/**
 *
 */
void
slabs_physics_set_filter_lower_mantle (
                                     ymir_vec_t *filter,
                                     slabs_physics_options_t *physics_options);

/**
 *
 */
void
slabs_physics_set_filter_upper_mantle (
                                     ymir_vec_t *filter,
                                     slabs_physics_options_t *physics_options);

/**
 *
 */
void
slabs_physics_set_filter_lithosphere_layer (
                                     ymir_vec_t *filter,
                                     slabs_physics_options_t *physics_options);

/**
 *
 */
void
slabs_physics_set_filter_plates (ymir_vec_t *filter,
                                 ymir_vec_t *viscosity,
                                 slabs_physics_options_t *physics_options);

/**
 *
 */
void
slabs_physics_set_filter_asthenosphere_away_from_subd (
                                     ymir_vec_t *filter,
                                     ymir_vec_t *viscosity,
                                     slabs_physics_options_t *physics_options);

/**
 *
 */
void
slabs_physics_set_filter_away_from_plate_bndr (
                                     ymir_vec_t *filter,
                                     slabs_physics_options_t *physics_options);

/**
 *
 */
void
slabs_physics_stats_quantity_in_lower_mantle (
                                     double *min, double *max, double *mean,
                                     ymir_vec_t *quantity,
                                     slabs_physics_options_t *physics_options);

/**
 *
 */
void
slabs_physics_stats_quantity_in_plates (
                                     double *min, double *max, double *mean,
                                     ymir_vec_t *quantity,
                                     ymir_vec_t *viscosity_nondim,
                                     slabs_physics_options_t *physics_options);

/**
 *
 */
void
slabs_physics_stats_quantity_in_asthenosphere (
                                     double *min, double *max, double *mean,
                                     ymir_vec_t *quantity,
                                     ymir_vec_t *viscosity_nondim,
                                     slabs_physics_options_t *physics_options);

/**
 *
 */
void
slabs_physics_stats_velocity_away_from_plate_bndr (
                                     double *min, double *max, double *mean,
                                     ymir_vec_t *velocity,
                                     slabs_physics_options_t *physics_options);

/**
 *
 */
double
slabs_physics_stats_plates_vol (ymir_vec_t *viscosity,
                                slabs_physics_options_t *physics_options);

/**
 *
 */
double
slabs_physics_stats_yielding_vol (ymir_vec_t *yielding_marker);

/**
 *
 */
double
slabs_physics_stats_plates_with_yielding_vol (
                                     ymir_vec_t *viscosity,
                                     ymir_vec_t *yielding_marker,
                                     slabs_physics_options_t *physics_options);

/**
 *
 */
void
slabs_physics_stats_print_all (ymir_vec_t *vel_press,
                               ymir_vec_t *viscosity,
                               ymir_vec_t *yielding_marker,
                               ymir_vec_t *rhs_u_point,
                               ymir_stokes_op_t *stokes_op,
                               slabs_physics_options_t *physics_options);

/**
 * TODO del
 */
void
slabs_physics_get_plate_velocity (double *vel_magn_subdu,
                                  double *vel_magn_overr,
                                  slabs_stokes_state_t *state,
                                  slabs_physics_options_t *physics_options);

/**
 *
 */
void
slabs_physics_transform_to_dimensional_temperature (ymir_vec_t * temp);

/**
 *
 */
void
slabs_physics_transform_to_dimensional_viscosity (ymir_vec_t * visc);

/**
 *
 */
void
slabs_physics_transform_to_dimensional_velocity (ymir_vec_t * vel);

/**
 *
 */
void
slabs_physics_transform_to_dimensional_strain_rate (ymir_vec_t * strain_rate);

/**
 *
 */
void
slabs_physics_transform_to_dimensional_stress (ymir_vec_t * stress);

/**
 *
 */
void
slabs_physics_transform_to_dimensional_pressure (ymir_vec_t * press);

#endif /* SLABS_PHYSICS_H */
