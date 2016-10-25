/*
 */

#ifndef RHEA_STOKES_H
#define RHEA_STOKES_H

#include <rhea_domain.h>
#include <ymir_pressure_elem.h>

/* linear Stokes problem (opaque) */
typedef struct rhea_stokes_linear_problem rhea_stokes_linear_problem_t;

/**
 * Creates a new linear Stokes problem.
 */
rhea_stokes_linear_problem_t *rhea_stokes_linear_problem_new (
                                    ymir_vec_t *temperature,
                                    ymir_vec_t *viscosity,
                                    const int viscosity_is_owned,
                                    ymir_mesh_t *ymir_mesh,
                                    ymir_pressure_elem_t *press_elem,
                                    rhea_domain_options_t *domain_options);

/**
 * Destroys a linear Stokes problem.
 */
void                rhea_stokes_linear_problem_destroy (
                                    rhea_stokes_linear_problem_t *lin_stokes);

#endif /* RHEA_STOKES_H */
