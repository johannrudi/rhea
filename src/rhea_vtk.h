/*
 */

#ifndef RHEA_VTK_H
#define RHEA_VTK_H

#include <ymir_vec_ops.h>

void                rhea_vtk_write_simple (const char *filepath,
                                           ymir_vec_t *temperature,
                                           ymir_vec_t *background_temp,
                                           ymir_vec_t *weakzone,
                                           ymir_vec_t *viscosity,
                                           ymir_vec_t *rhs_vel);

#endif /* RHEA_VTK_H */
