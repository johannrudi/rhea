/*
 */

#ifndef RHEA_INVERSION_OBS_STRESS_H
#define RHEA_INVERSION_OBS_STRESS_H

#include <rhea_stokes_problem.h>

/* enumerator for stress observational data */
typedef enum
{
  RHEA_INVERSION_OBS_STRESS_NONE = -1,
  RHEA_INVERSION_OBS_STRESS_VOLUME = 0,
  RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_NORMAL,
  RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_0,
  RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_1,
  RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_2
}
rhea_inversion_obs_stress_t;

/**
 * Compute the difference of the stresses:
 *   weight * (ObsOp(stress(vel,press)) - stress_obs)
 * where
 *   ObsOp(stress(vel,press)) = 2 * visc * strainrate(vel) - press * I
 */
void                rhea_inversion_obs_stress_diff (
                                  ymir_vec_t *misfit_stress,
                                  ymir_vec_t *forward_vel,
                                  ymir_vec_t *obs_stress,
                                  ymir_vec_t *weight,
                                  rhea_stokes_problem_t *stokes_problem);

/**
 * Computes misfit term (squared norm of the difference) of the of
 * the difference in stresses.
 */
double              rhea_inversion_obs_stress_misfit (
                                  ymir_vec_t *forward_vel,
                                  ymir_vec_t *obs_stress,
                                  ymir_vec_t *weight,
                                  const rhea_inversion_obs_stress_t obs_type,
                                  rhea_stokes_problem_t *stokes_problem);

/**
 * Adds contribution to the right-hand side for adjoint equations.
 */
void                rhea_inversion_obs_stress_add_adjoint_rhs (
                                  ymir_vec_t *rhs_vel_mass,
                                  ymir_vec_t *forward_vel,
                                  ymir_vec_t *obs_stress,
                                  ymir_vec_t *weight,
                                  const rhea_inversion_obs_stress_t obs_type,
                                  rhea_stokes_problem_t *stokes_problem);

#endif /* RHEA_INVERSION_OBS_STRESS_H */
