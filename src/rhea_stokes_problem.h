/*
 */

#ifndef RHEA_STOKES_H
#define RHEA_STOKES_H

#include <rhea_domain.h>
#include <rhea_temperature.h>
#include <rhea_viscosity.h>
#include <ymir_pressure_elem.h>

/* linear Stokes problem (opaque) */
typedef struct rhea_stokes_problem rhea_stokes_problem_t;

/**
 * Creates a new linear Stokes problem.
 */
rhea_stokes_problem_t  *rhea_stokes_problem_linear_new (
                                    ymir_vec_t *temperature,
                                    ymir_vec_t *weakzone,
                                    ymir_mesh_t *ymir_mesh,
                                    ymir_pressure_elem_t *press_elem,
                                    rhea_domain_options_t *domain_options,
                                    rhea_temperature_options_t *temp_options,
                                    rhea_viscosity_options_t *visc_options);

/**
 * Destroys a linear Stokes problem.
 */
void                rhea_stokes_problem_linear_destroy (
                                    rhea_stokes_problem_t *lin_stokes);

#endif /* RHEA_STOKES_H */
