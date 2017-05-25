/*
 */

#ifndef RHEA_H
#define RHEA_H

#include <rhea_base.h>
#include <rhea_domain.h>
#include <rhea_discretization.h>
#include <rhea_temperature.h>
#include <rhea_weakzone.h>
#include <rhea_viscosity.h>
#include <rhea_velocity.h>
#include <rhea_pressure.h>
#include <rhea_velocity_pressure.h>
#include <rhea_stokes_problem.h>
#include <rhea_newton.h>
#include <rhea_vtk.h>

/**
 * Begin the initialization of a program powered by rhea.  Should be followed
 * by calling `rhea_init_end (...)`.
 *
 * Initializes the rhea library and dependent libraries.  Retrieves parameters
 * of the parallel envirionment.
 */
void                rhea_init_begin (int *mpisize, int *mpirank, int *ompsize,
                                     int argc, char **argv, MPI_Comm mpicomm);

/**
 * Ends the initialization of a program powered by rhea.  Should follow after
 * calling `rhea_init_begin (...)`.
 *
 * Parses options and sets up ymir library.
 */
void                rhea_init_end (ymir_options_t *opt);

/**
 * Get whether current program execution is flagged as a production run.
 */
int                 rhea_get_production_run ();

/**
 * Set current program execution as a production run.
 */
void                rhea_set_production_run (const int is_production_run);

/******************************************************************************
 * Options
 *****************************************************************************/

/**
 * Defines rhea options and adds them as sub-options.
 */
void                rhea_add_options_base (ymir_options_t *opt);

void                rhea_add_options_all (ymir_options_t * options);

void                rhea_add_options_newton (ymir_options_t *options);

/**
 * Processes all rhea options and stores them.
 */
void                rhea_process_options_all (
                              rhea_domain_options_t *domain_options,
                              rhea_temperature_options_t *temperature_options,
                              rhea_viscosity_options_t *viscosity_options,
                              rhea_discretization_options_t *discr_options,
                              rhea_newton_options_t *newton_options);

/**
 * Processes a subset of options and stores them.
 */
void                rhea_process_options_newton (
                              rhea_domain_options_t *domain_options,
                              rhea_discretization_options_t *discr_options,
                              rhea_newton_options_t *newton_options);

#endif /* RHEA_H */
