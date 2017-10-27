#ifndef RHEA_STOKES_PROBLEM_AMR_H
#define RHEA_STOKES_PROBLEM_AMR_H

#include <rhea_discretization.h>
#include <rhea_viscosity.h>

int                 rhea_stokes_problem_amr (
                                  p4est_t *p4est,
                                  ymir_mesh_t **ymir_mesh,
                                  ymir_pressure_elem_t **press_elem,
                                  ymir_vec_t **temperature,
                                  ymir_vec_t **velocity_pressure,
                                  rhea_discretization_options_t *discr_options,
                                  rhea_viscosity_options_t *visc_options);

#endif /* RHEA_STOKES_PROBLEM_AMR_H */
