/*
 */

#ifndef RHEA_STRAINRATE_H
#define RHEA_STRAINRATE_H

#include <ymir_vec_ops.h>

/**
 * Creates a new second invariant of the strain rate vector.
 */
ymir_vec_t         *rhea_strainrate_2inv_new (ymir_mesh_t *ymir_mesh);

/**
 * Destroys a second invariant of the strain rate vector.
 */
void                rhea_strainrate_2inv_destroy (ymir_vec_t *strainrate_2inv);

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_strainrate_2inv_check_vec_type (ymir_vec_t *vec);

/**
 * Checks entries of a vector.
 */
int                 rhea_strainrate_2inv_is_valid (ymir_vec_t *vec);

/**
 * Computes the squre root of the second invariant of the strain rate.
 */
void                rhea_strainrate_compute_sqrt_of_2inv (
                                              ymir_vec_t *strainrate_sqrt_2inv,
                                              ymir_vec_t *velocity);

#endif /* RHEA_STRAINRATE_H */
