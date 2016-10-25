/*
 */

#ifndef RHEA_H
#define RHEA_H

#include <rhea_base.h>
#include <rhea_domain.h>
#include <rhea_discretization.h>
#include <rhea_temperature.h>
#include <rhea_viscosity.h>
#include <rhea_stokes.h>

/**
 * Defines all rhea options and adds them as sub-options.
 */
void                rhea_add_options_all (ymir_options_t * options);

/**
 * Processes all rhea options and stores them.
 */
void                rhea_process_options_all (
                                  rhea_domain_options_t *domain_options,
                                  rhea_discretization_options_t *discr_options,
                                  rhea_viscosity_options_t *viscosity_options);

#endif /* RHEA_H */
