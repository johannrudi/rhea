/*
 */

#ifndef RHEA_VELOCITY_H
#define RHEA_VELOCITY_H

#include <ymir_vec_ops.h>

/**
 * Gets the velocity of one element at GLL nodes
 */
void                rhea_velocity_get_elem_gll (sc_dmatrix_t *vel_el_mat,
                                                ymir_vec_t *vel_vec,
                                                const ymir_locidx_t elid);

/**
 * Computes the square root of 2nd invariant of the strain rate,
 *
 *   2nd inv. = sqrt ( 1/2 * grad_sym (velocity) : grad_sym (velocity) ),
 *
 * of one element at Gauss nodes.
 */
void                rhea_velocity_compute_strain_rate_2inv_elem_gauss (
                                        sc_dmatrix_t *strain_rate_2inv_el_mat,
                                        sc_dmatrix_t *vel_el_mat,
                                        ymir_vec_t *vel_vec,
                                        const ymir_locidx_t elid,
                                        sc_dmatrix_t *tmp_grad_vel,
                                        sc_dmatrix_t *tmp_dvel,
                                        sc_dmatrix_t *tmp_vel);

#endif /* RHEA_VELOCITY_H */
