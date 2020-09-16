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

double              rhea_inversion_obs_viscosity_misfit (
                                ymir_vec_t *forward_vel,
                                const rhea_inversion_obs_viscosity_t obs_type,
                                const int n_columns,
                                rhea_domain_subset_column_t **column,
                                const double *obs_val,
                                const double *weight,
                                rhea_stokes_problem_t *stokes_problem);

#endif /* RHEA_INVERSION_OBS_VISCOSITY_H */
