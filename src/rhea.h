/* RHEA  Main header file that should be included in applications/examples. */

#ifndef RHEA_H
#define RHEA_H

#include <rhea_base.h>
#include <rhea_error_stats.h>
#include <rhea_domain.h>
#include <rhea_discretization.h>
#include <rhea_temperature.h>
#include <rhea_composition.h>
#include <rhea_plate.h>
#include <rhea_weakzone.h>
#include <rhea_topography.h>
#include <rhea_viscosity.h>
#include <rhea_velocity.h>
#include <rhea_pressure.h>
#include <rhea_velocity_pressure.h>
#include <rhea_strainrate.h>
#include <rhea_stress.h>
#include <rhea_stokes_problem.h>
#include <rhea_inversion.h>
#include <rhea_newton.h>
#include <rhea_vis.h>
#include <rhea_vtk.h>

/**
 * Begin the initialization of a program powered by rhea.  Should be followed
 * by calling `rhea_init_end (...)`.
 *
 * Initializes the rhea library and dependent libraries.  Retrieves parameters
 * of the parallel envirionment.
 */
void                rhea_init_begin (int *mpisize, int *mpirank, int *ompsize,
                                     int argc, char **argv,
                                     sc_MPI_Comm mpicomm);

/**
 * Ends the initialization of a program powered by rhea.  Should follow after
 * calling `rhea_init_begin (...)`.
 *
 * Parses options and sets up ymir library.
 */
void                rhea_init_end (ymir_options_t *options);

/**
 * Get whether the program execution is flagged as a production run.
 */
int                 rhea_production_run_get ();

/**
 * Set the program execution as a production run.
 */
void                rhea_production_run_set (const int is_production_run);

/******************************************************************************
 * Options
 *****************************************************************************/

/* collection of Rhea's options */
typedef struct rhea_all_options
{
  rhea_domain_options_t          *domain_options;
  rhea_temperature_options_t     *temperature_options;
  rhea_composition_options_t     *composition_options;
  rhea_plate_options_t           *plate_options;
  rhea_weakzone_options_t        *weakzone_options;
  rhea_topography_options_t      *topography_options;
  rhea_viscosity_options_t       *viscosity_options;
  rhea_discretization_options_t  *discr_options;
}
rhea_all_options_t;

/**
 * Defines rhea options and adds them as sub-options.
 */
void                rhea_add_options_base (ymir_options_t *options);

void                rhea_add_options_all (ymir_options_t * options);

void                rhea_add_options_newton (ymir_options_t *options);

/**
 * Processes all rhea options and stores them.
 */
void                rhea_process_options_all (rhea_all_options_t *all_options);

/**
 * Processes a subset of options and stores them.
 */
void                rhea_process_options_newton (
                              rhea_domain_options_t *domain_options,
                              rhea_discretization_options_t *discr_options,
                              rhea_newton_options_t *newton_options);

/**
 * Print options that describe flow physics.
 */
void                rhea_print_physics_const_options (
                              rhea_domain_options_t *domain_options,
                              rhea_temperature_options_t *temperature_options,
                              rhea_viscosity_options_t *viscosity_options);

/******************************************************************************
 * Monitoring
 *****************************************************************************/

/**
 * Get whether performance monitoring is active.
 */
int                 rhea_performance_monitor_active ();

/**
 * Initializes performance monitors.
 */
void                rhea_performance_monitor_init (const char **monitor_name,
                                                   const int n_monitors);

/**
 * Finalizes performance monitors.
 */
void                rhea_performance_monitor_finalize ();

/* output verbosity for performance monitors */
typedef enum
{
  RHEA_PERFMON_PRINT_WTIME_NONE = 0,
  RHEA_PERFMON_PRINT_WTIME_ESSENTIAL,
  RHEA_PERFMON_PRINT_WTIME_ALL
}
rhea_performance_monitor_print_wtime_t;

typedef enum
{
  RHEA_PERFMON_PRINT_NCALLS_NONE = 0,
  RHEA_PERFMON_PRINT_NCALLS_ESSENTIAL,
  RHEA_PERFMON_PRINT_NCALLS_ALL
}
rhea_performance_monitor_print_ncalls_t;

typedef enum
{
  RHEA_PERFMON_PRINT_FLOPS_NONE = 0,
  RHEA_PERFMON_PRINT_FLOPS_ALL
}
rhea_performance_monitor_print_flops_t;

typedef enum
{
  RHEA_PERFMON_PRINT_YMIR_NONE = 0,
  RHEA_PERFMON_PRINT_YMIR_ALL
}
rhea_performance_monitor_print_ymir_t;

/**
 * Prints statistics collected by performance monitors.
 */
void                rhea_performance_monitor_print (
                        const char *title,
                        rhea_performance_monitor_print_wtime_t print_wtime,
                        rhea_performance_monitor_print_ncalls_t print_n_calls,
                        rhea_performance_monitor_print_flops_t print_flops,
                        rhea_performance_monitor_print_ymir_t print_ymir);

/**
 * Starts/stops a single performance monitor.
 */
void                rhea_performance_monitor_start (const int monitor_index);
void                rhea_performance_monitor_stop_add (const int monitor_index);

/**
 * Starts/stops a single performance monitor with MPI barrier.
 */
void                rhea_performance_monitor_start_barrier (
                                                      const int monitor_index);
void                rhea_performance_monitor_stop_add_barrier (
                                                      const int monitor_index);

#endif /* RHEA_H */
