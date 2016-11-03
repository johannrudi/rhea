/*
 */

#ifndef RHEA_PRESSURE_H
#define RHEA_PRESSURE_H

#include <ymir_vec_ops.h>
#include <ymir_pressure_elem.h>

/**
 * Creates a new pressure vector.
 */
ymir_vec_t         *rhea_pressure_new (ymir_mesh_t *ymir_mesh,
                                       ymir_pressure_elem_t *press_elem);

/**
 * Destroys a pressure vector.
 */
void                rhea_pressure_destroy (ymir_vec_t *pressure);

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_pressure_check_vec_type (
                                            ymir_vec_t *vec,
                                            ymir_pressure_elem_t *press_elem);

/**
 * Checks entries of a vector.
 */
int                 rhea_pressure_is_valid (ymir_vec_t *vec);

#endif /* RHEA_PRESSURE_H */

