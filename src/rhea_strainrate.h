/*
 */

#ifndef RHEA_STRAINRATE_H
#define RHEA_STRAINRATE_H

#include <ymir_vec_ops.h>

/* constant: neutral/default value for the sqrt of the 2nd invariant of the
 * strain rate */
#define RHEA_STRAINRATE_2INV_NEUTRAL_VALUE (1.0)

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
 * Computes the squre root of the second invariant of the strain rate.
 */
void                rhea_strainrate_compute_sqrt_of_2inv (
                                              ymir_vec_t *strainrate_sqrt_2inv,
                                              ymir_vec_t *velocity);

/**
 * Computes the square root of the second invariant of the strain rate at Gauss
 * nodes, given the velocity at GLL nodes.
 */
void                rhea_strainrate_compute_sqrt_of_2inv_elem_gauss (
                                        sc_dmatrix_t *strainrate_2inv_el_mat,
                                        sc_dmatrix_t *vel_el_mat,
                                        ymir_vec_t *vel_vec,
                                        const ymir_locidx_t elid,
                                        sc_dmatrix_t *tmp_grad_vel,
                                        sc_dmatrix_t *tmp_dvel,
                                        sc_dmatrix_t *tmp_vel);

#endif /* RHEA_STRAINRATE_H */
