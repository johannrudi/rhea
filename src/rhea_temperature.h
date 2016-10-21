/*
 */

#ifndef RHEA_TEMPERATURE_H
#define RHEA_TEMPERATURE_H

#include <ymir_vec_ops.h>

/**
 * Gets the temperature of one element at Gauss nodes.
 */
double             *rhea_temperature_get_elem_gauss (sc_dmatrix_t *temp_el_mat,
                                                     ymir_vec_t *temp_vec,
                                                     const ymir_locidx_t elid);

#endif /* RHEA_TEMPERATURE_H */
