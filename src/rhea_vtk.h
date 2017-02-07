/*
 */

#ifndef RHEA_VTK_H
#define RHEA_VTK_H

#include <ymir_vec_ops.h>
#include <ymir_pressure_elem.h>

void                rhea_vtk_write_input_data (const char *filepath,
                                               ymir_vec_t *temperature,
                                               ymir_vec_t *background_temp,
                                               ymir_vec_t *weakzone,
                                               ymir_vec_t *viscosity,
                                               ymir_vec_t *rhs_vel);

void                rhea_vtk_write_solution (const char *filepath,
                                             ymir_vec_t *velocity,
                                             ymir_vec_t *pressure,
                                             ymir_vec_t *viscosity);

void                rhea_vtk_write_nonlinear_stokes_iteration (
                                            const char *filepath,
                                            ymir_vec_t *velocity,
                                            ymir_vec_t *pressure,
                                            ymir_vec_t *viscosity,
                                            ymir_vec_t *bounds_marker,
                                            ymir_vec_t *yielding_marker);

#endif /* RHEA_VTK_H */
