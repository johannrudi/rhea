/*
 */

#ifndef RHEA_VELOCITY_H
#define RHEA_VELOCITY_H

#include <ymir_vec_ops.h>

/**
 * Creates a new velocity vector.
 */
ymir_vec_t         *rhea_velocity_new (ymir_mesh_t *ymir_mesh);

/**
 * Destroys a velocity vector.
 */
void                rhea_velocity_destroy (ymir_vec_t *velocity);

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_velocity_check_vec_type (ymir_vec_t *vec);

/**
 * Checks entries of a vector.
 */
int                 rhea_velocity_is_valid (ymir_vec_t *vec);

/**
 * Gets the velocity of one element at GLL nodes.
 */
void                rhea_velocity_get_elem_gll (sc_dmatrix_t *vel_el_mat,
                                                ymir_vec_t *vel_vec,
                                                const ymir_locidx_t elid);

/**
 * Gets the velocity of one element at GLL nodes, then computes the square root
 * of 2nd invariant of the strain rate at Gauss nodes.
 */
void                rhea_velocity_get_elem_compute_strain_rate_2inv_gauss (
                                        sc_dmatrix_t *strain_rate_2inv_el_mat,
                                        sc_dmatrix_t *vel_el_mat,
                                        ymir_vec_t *vel_vec,
                                        const ymir_locidx_t elid,
                                        sc_dmatrix_t *tmp_grad_vel,
                                        sc_dmatrix_t *tmp_dvel,
                                        sc_dmatrix_t *tmp_vel);

#endif /* RHEA_VELOCITY_H */
