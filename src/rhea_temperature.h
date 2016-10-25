/*
 */

#ifndef RHEA_TEMPERATURE_H
#define RHEA_TEMPERATURE_H

#include <rhea_domain.h>
#include <ymir_vec_ops.h>

/* enumerator for types of temperature */
typedef enum
{
  RHEA_TEMPERATURE_NONE,
  RHEA_TEMPERATURE_COLD_PLATE,
  RHEA_TEMPERATURE_2PLATES_POLY2,
  RHEA_TEMPERATURE_IMPORT
}
rhea_temperature_t;

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_temperature_add_options (ymir_options_t * opt_sup);

/**
 * Processes options and stores them.
 */
void                rhea_temperature_process_options ();

/**
 * Creates a new temperature vector.
 */
ymir_vec_t         *rhea_temperature_new (ymir_mesh_t *ymir_mesh);

/**
 * Destroys a temperature vector.
 */
void                rhea_temperature_destroy (ymir_vec_t *temperature);

/**
 * Gets the temperature of one element at Gauss nodes.
 */
double             *rhea_temperature_get_elem_gauss (sc_dmatrix_t *temp_el_mat,
                                                     ymir_vec_t *temp_vec,
                                                     const ymir_locidx_t elid);

/**
 * Computes velocity right-hand side in (primal) function space, given a
 * temperature vector.
 */
void                rhea_temperature_compute_rhs_vel (
                                        ymir_vec_t *rhs_vel,
                                        ymir_vec_t *temperature,
                                        rhea_domain_options_t *domain_options);

#endif /* RHEA_TEMPERATURE_H */
