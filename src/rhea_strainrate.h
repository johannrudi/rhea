/* RHEA_STRAINRATE  Strain rate tensor and its second invariant. */

#ifndef RHEA_STRAINRATE_H
#define RHEA_STRAINRATE_H

#include <rhea_domain.h>
#include <rhea_temperature.h>

/* constant: neutral/default value for the sqrt of the 2nd invariant of the
 * strain rate */
#define RHEA_STRAINRATE_2INV_NEUTRAL_VALUE (1.0)

/******************************************************************************
 * Options
 *****************************************************************************/

/**
 * Gets the scaling factor to convert nondimensional values to corresponding
 * dimensional quantities.
 *   Unit: [1/s]
 */
double              rhea_strainrate_get_dim_1_s (
                                    rhea_domain_options_t *domain_options,
                                    rhea_temperature_options_t *temp_options);

/******************************************************************************
 * Strain Rate Vector
 *****************************************************************************/

/**
 * Creates a new (second invariant of) strain rate tensor.
 */
ymir_vec_t         *rhea_strainrate_new (ymir_mesh_t *ymir_mesh);

ymir_vec_t         *rhea_strainrate_2inv_new (ymir_mesh_t *ymir_mesh);

/**
 * Destroys a (second invariant of) strain rate tensor.
 */
void                rhea_strainrate_destroy (ymir_vec_t *strainrate);

void                rhea_strainrate_2inv_destroy (ymir_vec_t *strainrate_2inv);

/**
 * Converts entries of a nondimensional strain rate vector into dimensional
 * values:
 *   Unit: [1/s]
 */
void                rhea_strainrate_convert_to_dimensional_1_s (
                                    ymir_vec_t * strainrate,
                                    rhea_domain_options_t *domain_options,
                                    rhea_temperature_options_t *temp_options);

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_strainrate_check_vec_type (ymir_vec_t *vec);

int                 rhea_strainrate_2inv_check_vec_type (ymir_vec_t *vec);

/**
 * Checks entries of a vector.
 */
int                 rhea_strainrate_is_valid (ymir_vec_t *vec);

int                 rhea_strainrate_2inv_is_valid (ymir_vec_t *vec);

/**
 * Computes the strain rate tensor.
 */
void                rhea_strainrate_compute (ymir_vec_t *strainrate,
                                             ymir_vec_t *velocity);

/**
 * Computes the squre root of the second invariant of the strain rate.
 */
void                rhea_strainrate_compute_sqrt_of_2inv (
                                              ymir_vec_t *strainrate_sqrt_2inv,
                                              ymir_vec_t *velocity);

/******************************************************************************
 * Get & Set Values
 *****************************************************************************/

/**
 * Computes the square root of the second invariant of the strain rate at Gauss
 * nodes, given the velocity at GLL nodes.
 */
void                rhea_strainrate_compute_sqrt_of_2inv_elem (
                                        sc_dmatrix_t *strainrate_2inv_el_mat,
                                        sc_dmatrix_t *vel_el_mat,
                                        ymir_vec_t *vel_vec,
                                        const ymir_locidx_t elid,
                                        sc_dmatrix_t *tmp_grad_vel,
                                        sc_dmatrix_t *tmp_dvel,
                                        sc_dmatrix_t *tmp_vel);

/******************************************************************************
 * Statistics
 *****************************************************************************/

/**
 * Computes global strain rate statistics.
 */
void                rhea_strainrate_stats_get_global (
                                    double *min_1_s, double *max_1_s,
                                    double *mean_1_s, ymir_vec_t *velocity,
                                    rhea_domain_options_t *domain_options,
                                    rhea_temperature_options_t *temp_options);

#endif /* RHEA_STRAINRATE_H */
