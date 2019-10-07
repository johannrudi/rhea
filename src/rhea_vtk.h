/* RHEA_VTK  Writes VTK files to disk. */

#ifndef RHEA_VTK_H
#define RHEA_VTK_H

#include <ymir_vec_ops.h>

int                 rhea_vtk_write_input_data (const char *filepath,
                                               ymir_vec_t *temperature,
                                               ymir_vec_t *background_temp,
                                               ymir_vec_t *weakzone,
                                               ymir_vec_t *viscosity,
                                               ymir_vec_t *marker,
                                               ymir_vec_t *rhs_vel);

int                 rhea_vtk_write_observation_data (const char *filepath,
                                                     ymir_vec_t *plate_label,
                                                     ymir_vec_t *plate_vel);

int                 rhea_vtk_write_solution (const char *filepath,
                                             ymir_vec_t *velocity,
                                             ymir_vec_t *pressure,
                                             ymir_vec_t *viscosity,
                                             ymir_vec_t *marker);

int                 rhea_vtk_write_solution_surf (const char *filepath,
                                                  ymir_vec_t *velocity_surf,
                                                  ymir_vec_t *stress_norm_surf,
                                                  ymir_vec_t *viscosity_surf,
                                                  ymir_vec_t *marker_surf);

int                 rhea_vtk_write_nonlinear_stokes_iteration (
                                                  const char *filepath,
                                                  ymir_vec_t *velocity,
                                                  ymir_vec_t *pressure,
                                                  ymir_vec_t *viscosity,
                                                  ymir_vec_t *marker);

int                 rhea_vtk_write_nonlinear_stokes_iteration_surf (
                                                  const char *filepath,
                                                  ymir_vec_t *velocity_surf,
                                                  ymir_vec_t *stress_norm_surf,
                                                  ymir_vec_t *viscosity_surf);

int                 rhea_vtk_write_inversion_iteration_surf (
                                                  const char *filepath,
                                                  ymir_vec_t *velocity_surf,
                                                  ymir_vec_t *velocity_obs_surf,
                                                  ymir_vec_t *misfit_surf);

#endif /* RHEA_VTK_H */
