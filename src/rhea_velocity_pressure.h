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

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_velocity_pressure_check_vec_type (ymir_vec_t *vec);

/**
 * Checks entries of a vector.
 */
int                 rhea_velocity_pressure_is_valid (ymir_vec_t *vec);

/**
 * Creates velocity and pressure vectors and gets data of each component from a
 * combined velocity-pressure vector.  Tries to obtain a view onto the data if
 * possible.
 */
int                 rhea_velocity_pressure_create_components (
                                            ymir_vec_t **vel,
                                            ymir_vec_t **press,
                                            ymir_vec_t *vel_press,
                                            ymir_pressure_elem_t *press_elem);

/**
 * Copies velocity and pressure from a combined velocity-pressure vector.
 */
void                rhea_velocity_pressure_copy_components (
                                            ymir_vec_t *vel,
                                            ymir_vec_t *press,
                                            ymir_vec_t *vel_press,
                                            ymir_pressure_elem_t *press_elem);

/**
 * Copies velocity and pressure to a combined velocity-pressure vector.
 */
void                rhea_velocity_pressure_set_components (
                                            ymir_vec_t *vel_press,
                                            ymir_vec_t *vel, ymir_vec_t *press,
                                            ymir_pressure_elem_t *press_elem);

/**
 * Reads velocity and pressure from binary files.
 */
int                 rhea_velocity_pressure_read (
                                            ymir_vec_t *vel_press,
                                            char *vel_file_path_bin,
                                            char *press_file_path_bin,
                                            ymir_pressure_elem_t *press_elem,
                                            sc_MPI_Comm mpicomm);

/**
 * Writes velocity and pressure to binary files.
 */
int                 rhea_velocity_pressure_write (
                                            char *vel_file_path_bin,
                                            char *press_file_path_bin,
                                            ymir_vec_t *vel_press,
                                            ymir_pressure_elem_t *press_elem,
                                            sc_MPI_Comm mpicomm);

#endif /* RHEA_VELOCITY_PRESSURE_H */
