/*
 */

#ifndef RHEA_VELOCITY_PRESSURE_H
#define RHEA_VELOCITY_PRESSURE_H

#include <ymir_vec_ops.h>
#include <ymir_pressure_elem.h>

/**
 * Creates a new (velocity,pressure) vector.
 */
ymir_vec_t         *rhea_velocity_pressure_new (
                                            ymir_mesh_t *ymir_mesh,
                                            ymir_pressure_elem_t *press_elem);

/**
 * Destroys a (velocity,pressure) vector.
 */
void                rhea_velocity_pressure_destroy (
                                            ymir_vec_t *velocity_pressure);

#endif /* RHEA_VELOCITY_PRESSURE_H */
