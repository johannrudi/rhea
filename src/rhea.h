/*
 */

#ifndef RHEA_H
#define RHEA_H

#include <rhea_base.h>
#include <rhea_domain.h>
#include <rhea_discretization.h>
#include <rhea_temperature.h>
#include <rhea_viscosity.h>
#include <rhea_velocity.h>
#include <rhea_pressure.h>
#include <rhea_velocity_pressure.h>
#include <rhea_stokes_problem.h>
#include <rhea_vtk.h>

/**
 * Defines all rhea options and adds them as sub-options.
 */
void                rhea_add_options_all (ymir_options_t * options);

/**
 * Processes all rhea options and stores them.
 */
void                rhea_process_options_all (
                              rhea_domain_options_t *domain_options,
                              rhea_temperature_options_t *temperature_options,
                              rhea_viscosity_options_t *viscosity_options,
                              rhea_discretization_options_t *discr_options);

/**
 * Get whether current program execution is flagged as a production run.
 */
int                 rhea_get_production_run ();

/**
 * Set current program execution as a production run.
 */
void                rhea_set_production_run (const int is_production_run);

#endif /* RHEA_H */
