/*
 */

#ifndef RHEA_WEAKZONE_H
#define RHEA_WEAKZONE_H

#include <ymir_vec_ops.h>

/**
 * Gets the weak zone of one element at Gauss nodes.
 */
void                rhea_weakzone_get_elem_gauss (sc_dmatrix_t *weak_el_mat,
                                                  ymir_vec_t *weak_vec,
                                                  const ymir_locidx_t elid);

/**
 * Sets the weak zone of one element at Gauss nodes.
 */
void                rhea_weakzone_set_elem_gauss (ymir_vec_t *weak_vec,
                                                  sc_dmatrix_t *weak_el_mat,
                                                  const ymir_locidx_t elid);

#endif /* RHEA_WEAKZONE_H */
