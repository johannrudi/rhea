/*
 */

#ifndef RHEA_VELOCITY_H
#define RHEA_VELOCITY_H

#include <rhea_domain.h>
#include <rhea_temperature.h>

/******************************************************************************
 * Velocity Vector
 *****************************************************************************/

/**
 * Creates a new velocity vector.
 */
ymir_vec_t         *rhea_velocity_new (ymir_mesh_t *ymir_mesh);

/**
 * Destroys a velocity vector.
 */
void                rhea_velocity_destroy (ymir_vec_t *velocity);

/**
 * Converts entries of a nondimensional velocity vector into dimensional
 * values: [m/s]
 */
void                rhea_velocity_convert_to_dimensional (
                                    ymir_vec_t * velocity,
                                    rhea_domain_options_t *domain_options,
                                    rhea_temperature_options_t *temp_options);

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_velocity_check_vec_type (ymir_vec_t *vec);

/**
 * Checks entries of a vector.
 */
int                 rhea_velocity_is_valid (ymir_vec_t *vec);

/**
 * Creates a new velocity vector at surface.
 */
ymir_vec_t         *rhea_velocity_surface_new (ymir_mesh_t *ymir_mesh);

/**
 * Destroys a velocity vector at surface.
 */
void                rhea_velocity_surface_destroy (ymir_vec_t *vel_surf);

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_velocity_surface_check_vec_type (ymir_vec_t *vec);

/**
 * Checks entries of a vector.
 */
int                 rhea_velocity_surface_is_valid (ymir_vec_t *vec);

/**
 * Creates a new velocity vector at surface with values interpolated from a
 * velocity volume vector.
 */
ymir_vec_t         *rhea_velocity_surface_new_from_vol (ymir_vec_t *vel_vol);

/**
 * Interpolates velocity from volume to surface.
 */
void                rhea_velocity_surface_interpolate (ymir_vec_t *vel_surf,
                                                       ymir_vec_t *vel_vol);

/******************************************************************************
 * Get & Set Values
 *****************************************************************************/

/**
 * Gets the velocity of one element at GLL nodes.
 */
void                rhea_velocity_get_elem_gll (sc_dmatrix_t *vel_el_mat,
                                                ymir_vec_t *vel_vec,
                                                const ymir_locidx_t elid);

/******************************************************************************
 * Right-Hand Side Computation
 *****************************************************************************/

typedef void      (*rhea_velocity_rhs_compute_fn_t) (ymir_vec_t *rhs_vel,
                                                     ymir_vec_t *temperature,
                                                     void *data);

typedef void      (*rhea_velocity_rhs_nz_dir_compute_fn_t) (
                                        ymir_vec_t *rhs_vel_nonzero_dirichlet,
                                        void *data);

void                rhea_velocity_rhs_compute (ymir_vec_t *rhs_vel,
                                              ymir_vec_t *temperature,
                                              void *data);

/******************************************************************************
 * Statistics
 *****************************************************************************/

void                rhea_velocity_stats_get_global (
                                double *magn_cm_yr_min, double *magn_cm_yr_max,
                                double *magn_cm_yr_mean, ymir_vec_t *velocity,
                                rhea_domain_options_t *domain_options,
                                rhea_temperature_options_t *temp_options);

void                rhea_velocity_stats_get_global_surface (
                                double *magn_cm_yr_min, double *magn_cm_yr_max,
                                double *magn_cm_yr_mean, ymir_vec_t *vel_vol,
                                rhea_domain_options_t *domain_options,
                                rhea_temperature_options_t *temp_options);

#endif /* RHEA_VELOCITY_H */
