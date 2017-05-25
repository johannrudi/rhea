/*
 */

#ifndef RHEA_STRESS_H
#define RHEA_STRESS_H

#include <ymir_vec_ops.h>

/**
 * Creates a new (second invariant of) stress tensor.
 */
ymir_vec_t         *rhea_stress_new (ymir_mesh_t *ymir_mesh);

ymir_vec_t         *rhea_stress_2inv_new (ymir_mesh_t *ymir_mesh);

/**
 * Destroys a (second invariant of) stress tensor.
 */
void                rhea_stress_destroy (ymir_vec_t *stress);

void                rhea_stress_2inv_destroy (ymir_vec_t *stress_2inv);

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_stress_check_vec_type (ymir_vec_t *vec);

int                 rhea_stress_2inv_check_vec_type (ymir_vec_t *vec);

/**
 * Checks entries of a vector.
 */
int                 rhea_stress_is_valid (ymir_vec_t *vec);

int                 rhea_stress_2inv_is_valid (ymir_vec_t *vec);

/**
 * Computes the square root of the second invariant of the viscous stress
 * tensor.
 */
void                rhea_stress_compute_viscstress_sqrt_of_2inv (
                                             ymir_vec_t *viscstress_sqrt_2inv,
                                             ymir_vec_t *strainrate_sqrt_2inv,
                                             ymir_vec_t *viscosity);

/**
 * Computes the norm of a stress tensor.
 */
double              rhea_stress_compute_norm (ymir_vec_t *stress);

#endif /* RHEA_STRESS_H */
