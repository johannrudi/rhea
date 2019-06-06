/** RHEA_INVERSION
 *
 * Inversion for Stokes parameters.
 */

#ifndef RHEA_INVERSION_H
#define RHEA_INVERSION_H

#include <rhea_stokes_problem.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_inversion_add_options (ymir_options_t * opt_sup);

/******************************************************************************
 * Monitoring
 *****************************************************************************/

/**
 * Initializes performance counters.
 */
void                rhea_inversion_perfmon_init (const int activate,
                                                 const int skip_if_active);

/**
 * Prints statistics collected by performance monitors.
 */
void                rhea_inversion_perfmon_print (sc_MPI_Comm mpicomm,
                                                  const int print_wtime,
                                                  const int print_n_calls,
                                                  const int print_flops);

/******************************************************************************
 * Inverse Problem
 *****************************************************************************/

/* inverse problem (opaque) */
typedef struct rhea_inversion_problem rhea_inversion_problem_t;

/**
 * Creates/destroys an inverse problem.
 */
rhea_inversion_problem_t *rhea_inversion_new (
                                        rhea_stokes_problem_t *stokes_problem);

void                      rhea_inversion_destroy (
                                        rhea_inversion_problem_t *inv_problem);

/**
 * Sets paths for writing output files.
 */
void                      rhea_inversion_set_txt_output (
                                        rhea_inversion_problem_t *inv_problem,
                                        char *txt_path);

void                      rhea_inversion_set_vtk_output (
                                        rhea_inversion_problem_t *inv_problem,
                                        char *vtk_path);

/******************************************************************************
 * Inverse Problem Solver
 *****************************************************************************/

/**
 * Solves an inverse problem.
 */
void                rhea_inversion_solve (
                                        rhea_inversion_problem_t *inv_problem,
                                        const int use_initial_guess,
                                        ymir_vec_t *inv_parameter_vec);

/**
 * Solves an inverse problem with given velocity observations at the surface.
 */
void                rhea_inversion_solve_with_vel_obs (
                                        rhea_inversion_problem_t *inv_problem,
                                        const int use_initial_guess,
                                        ymir_vec_t *inv_parameter_vec,
                                        ymir_vec_t *vel_obs_surf,
                                        ymir_vec_t *vel_obs_weight_surf,
                                        const double vel_obs_add_noise_stddev);

#endif /* RHEA_INVERSION_H */
