/*
 */

#ifndef RHEA_PRESSURE_H
#define RHEA_PRESSURE_H

#include <rhea_domain.h>
#include <rhea_temperature.h>
#include <rhea_viscosity.h>
#include <ymir_pressure_elem.h>

/******************************************************************************
 * Pressure Vector
 *****************************************************************************/

/**
 * Creates a new pressure vector.
 */
ymir_vec_t         *rhea_pressure_new (ymir_mesh_t *ymir_mesh,
                                       ymir_pressure_elem_t *press_elem);

/**
 * Destroys a pressure vector.
 */
void                rhea_pressure_destroy (ymir_vec_t *pressure);

/**
 * Converts entries of a nondimensional pressure vector into dimensional
 * values:
 *
 *   [Pa]
 */
void                rhea_pressure_convert_to_dimensional_Pa (
                                      ymir_vec_t * pressure,
                                      rhea_domain_options_t *domain_options,
                                      rhea_temperature_options_t *temp_options,
                                      rhea_viscosity_options_t *visc_options);

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_pressure_check_vec_type (
                                            ymir_vec_t *vec,
                                            ymir_pressure_elem_t *press_elem);

/**
 * Checks entries of a vector.
 */
int                 rhea_pressure_is_valid (ymir_vec_t *vec);

/**
 * Gets rank-global offsets or rank-local sizes of a distributed vector for
 * each MPI-rank.
 */
MPI_Offset         *rhea_pressure_segment_offset_create (ymir_vec_t *vec);
MPI_Offset          rhea_pressure_segment_offset_get (ymir_vec_t *vec);
int                 rhea_pressure_segment_size_get (ymir_vec_t *vec);

/**
 * Computes global mean value.
 */
double              rhea_pressure_compute_mean (
                                            ymir_vec_t *pressure,
                                            ymir_pressure_elem_t *press_elem);

/**
 * Projects out global mean value.
 */
int                 rhea_pressure_project_out_mean (
                                            ymir_vec_t *pressure,
                                            ymir_pressure_elem_t *press_elem);

/**
 * Inverts (lumped) mass matrix.
 */
void                rhea_pressure_remove_mass (
                                            ymir_vec_t *press,
                                            ymir_pressure_elem_t *press_elem);

/******************************************************************************
 * Statistics
 *****************************************************************************/

/**
 * Computes global pressure statistics.
 */
void                rhea_pressure_stats_get_global (
                                double *abs_min_Pa, double *abs_max_Pa,
                                double *mean_Pa, ymir_vec_t *pressure,
                                ymir_pressure_elem_t *press_elem,
                                rhea_domain_options_t *domain_options,
                                rhea_temperature_options_t *temp_options,
                                rhea_viscosity_options_t *visc_options);

#endif /* RHEA_PRESSURE_H */

