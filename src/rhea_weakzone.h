/*
 */

#ifndef RHEA_WEAKZONE_H
#define RHEA_WEAKZONE_H

#include <ymir_vec_ops.h>

/* constant: default value for weak zone (i.e., no weakening) */
#define RHEA_WEAKZONE_DEFAULT_VALUE (1.0)

/**
 * Creates a new weak zone vector.
 */
ymir_vec_t         *rhea_weakzone_new (ymir_mesh_t *ymir_mesh);

/**
 * Destroys a weak zone vector.
 */
void                rhea_weakzone_destroy (ymir_vec_t *weakzone);

/**
 * Gets the weak zone of one element at Gauss nodes.
 */
double             *rhea_weakzone_get_elem_gauss (sc_dmatrix_t *weak_el_mat,
                                                  ymir_vec_t *weak_vec,
                                                  const ymir_locidx_t elid);

/**
 * Sets the weak zone of one element at Gauss nodes.
 */
void                rhea_weakzone_set_elem_gauss (ymir_vec_t *weak_vec,
                                                  sc_dmatrix_t *weak_el_mat,
                                                  const ymir_locidx_t elid);

#endif /* RHEA_WEAKZONE_H */
