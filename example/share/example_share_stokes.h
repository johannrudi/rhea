/**
 * Shared functions for rhea examples.
 *
 ******************************************************************************
 * Author:             Johann Rudi <johann@ices.utexas.edu>
 *****************************************************************************/

#ifndef EXAMPLE_SHARE_STOKES_H
#define EXAMPLE_SHARE_STOKES_H

#include <rhea_domain.h>
#include <rhea_newton.h>
#include <rhea_stokes_problem.h>
#include <rhea_temperature.h>
#include <rhea_viscosity.h>

/**
 * Creates a new Stokes problem.
 */
void                example_share_stokes_new (
                                      rhea_stokes_problem_t **stokes_problem,
                                      ymir_mesh_t *ymir_mesh,
                                      ymir_pressure_elem_t *press_elem,
                                      rhea_domain_options_t *domain_options,
                                      rhea_temperature_options_t *temp_options,
                                      rhea_viscosity_options_t *visc_options,
                                      rhea_newton_options_t *newton_options,
                                      char *vtk_write_newton_itn_path);


/**
 * Destroys a Stokes problem.
 */
void                example_share_stokes_destroy (
                                      rhea_stokes_problem_t *stokes_problem);

#endif /* EXAMPLE_SHARE_STOKES_H */
