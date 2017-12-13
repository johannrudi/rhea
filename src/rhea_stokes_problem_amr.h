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
#define RHEA_STOKES_PROBLEM_AMR_VERBOSE (3)
#define RHEA_STOKES_PROBLEM_AMR_VERBOSE_VTK (3)

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_stokes_problem_amr_add_options (
                                                    ymir_options_t * opt_sup);

/**
 *
 */
int                 rhea_stokes_problem_init_amr (
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
