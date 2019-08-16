/** RHEA_STOKES_PROBLEM_AMR
 *
 * Provides the adaptive mesh refinement (AMR) for Stokes problems.
 */

#ifndef RHEA_STOKES_PROBLEM_AMR_H
#define RHEA_STOKES_PROBLEM_AMR_H

#include <rhea_stokes_problem.h>
#include <rhea_discretization.h>

/* set verbosity:
 *   0: all off
 *   1: only essentials
 *   2: concise information
 *   3: detailed information
 *   4: all on
 */
#define RHEA_STOKES_PROBLEM_AMR_VERBOSE (1)
#define RHEA_STOKES_PROBLEM_AMR_VERBOSE_VTK (0)

/******************************************************************************
 * Options
 *****************************************************************************/

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_stokes_problem_amr_add_options (
                                                    ymir_options_t * opt_sup);

/**
 * Processes options and stores them.
 */
void                rhea_stokes_problem_amr_process_options ();

/**
 * Sets AMR types.
 */
void                rhea_stokes_problem_amr_set_init_type_name (
                                                      const char *type_name);

void                rhea_stokes_problem_amr_set_nonlinear_type_name (
                                                      const char *type_name);

void                rhea_stokes_problem_amr_set_solution_type_name (
                                                      const char *type_name);

/******************************************************************************
 * Perform AMR
 *****************************************************************************/

/**
 * Runs AMR for the initial mesh.
 */
int                 rhea_stokes_problem_init_amr (
                                rhea_stokes_problem_t *stokes_problem,
                                p4est_t *p4est,
                                rhea_discretization_options_t *discr_options);

/**
 * Runs AMR at nonlinear iterations.
 */
int                 rhea_stokes_problem_nonlinear_amr (
                                rhea_stokes_problem_t *stokes_problem,
                                p4est_t *p4est,
                                rhea_discretization_options_t *discr_options,
                                const int nonlinear_iter);

/**
 * Runs AMR for the solution.
 */
int                 rhea_stokes_problem_solution_amr (
                                rhea_stokes_problem_t *stokes_problem,
                                p4est_t *p4est,
                                rhea_discretization_options_t *discr_options);

/**
 *
 */
int                 rhea_stokes_problem_amr (
                                rhea_stokes_problem_t *stokes_problem,
                                p4est_t *p4est,
                                rhea_discretization_options_t *discr_options);

#endif /* RHEA_STOKES_PROBLEM_AMR_H */
