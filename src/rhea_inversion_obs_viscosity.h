/*
 */

#ifndef RHEA_INVERSION_OBS_VISCOSITY_H
#define RHEA_INVERSION_OBS_VISCOSITY_H

#include <rhea_domain_subset.h>
#include <rhea_stokes_problem.h>

/* enumerator for viscosity observational data */
typedef enum
{
  RHEA_INVERSION_OBS_VISCOSITY_NONE = -1,
  RHEA_INVERSION_OBS_VISCOSITY_AVERAGE_REGION = 0,
  RHEA_INVERSION_OBS_VISCOSITY_AVERAGE_UNDER_PLATES
}
rhea_inversion_obs_viscosity_t;

/**
 * Creates observational data of avarage viscosities.
 */
rhea_domain_subset_column_t **rhea_inversion_obs_viscosity_new (
                                  int *n_columns,
                                  double **value,
                                  double **weight,
                                  const rhea_inversion_obs_viscosity_t obs_type,
                                  char *value_Pas_list,
                                  char *stddev_rel_list,
                                  rhea_plate_options_t *plate_options,
                                  rhea_viscosity_options_t *visc_options);

/**
 * Destroys observational data of avarage viscosities.
 */
void                rhea_inversion_obs_viscosity_destroy (
                                  rhea_domain_subset_column_t **column,
                                  const int n_columns,
                                  double *value,
                                  double *weight);

/**
 * Computes the misfit term of average viscosities.
 */
double              rhea_inversion_obs_viscosity_misfit (
                                ymir_vec_t *forward_vel,
                                const rhea_inversion_obs_viscosity_t obs_type,
                                const int n_columns,
                                rhea_domain_subset_column_t **column,
                                const double *obs_val,
                                const double *weight,
                                rhea_stokes_problem_t *stokes_problem);

/**
 * Adds contribution to the right-hand side for adjoint equations.
 */
void                rhea_inversion_obs_viscosity_add_adjoint_rhs (
                                ymir_vec_t *rhs_vel_mass,
                                ymir_vec_t *forward_vel,
                                const rhea_inversion_obs_viscosity_t obs_type,
                                const int n_columns,
                                rhea_domain_subset_column_t **column,
                                const double *obs_val,
                                const double *weight,
                                rhea_stokes_problem_t *stokes_problem);

#endif /* RHEA_INVERSION_OBS_VISCOSITY_H */
