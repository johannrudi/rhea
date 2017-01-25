/*
 */

#ifndef RHEA_VTK_H
#define RHEA_VTK_H

#include <ymir_vec_ops.h>

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
int                rhea_output_pressure (const char *filepath,
                                             ymir_vec_t *pressure);

#endif /* RHEA_VTK_H */
