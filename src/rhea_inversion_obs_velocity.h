/*
 */

#ifndef RHEA_INVERSION_OBS_VELOCITY_H
#define RHEA_INVERSION_OBS_VELOCITY_H

#include <rhea_plate.h>
#include <ymir_vec_ops.h>

/* enumerator for velocity observational data */
typedef enum
{
  RHEA_INVERSION_OBS_VELOCITY_NONE = -1,
  RHEA_INVERSION_OBS_VELOCITY_NORMAL = 0,
  RHEA_INVERSION_OBS_VELOCITY_TANGENTIAL,
  RHEA_INVERSION_OBS_VELOCITY_TANGENTIAL_ROTFREE,
  RHEA_INVERSION_OBS_VELOCITY_ALL,
  RHEA_INVERSION_OBS_VELOCITY_ALL_ROTFREE
}
rhea_inversion_obs_velocity_t;

/* enumerator for weights of observational data */
typedef enum
{
  RHEA_INVERSION_OBS_VELOCITY_WEIGHT_UNIFORM = 0,
  RHEA_INVERSION_OBS_VELOCITY_WEIGHT_INV_AREA_SQRT,
  RHEA_INVERSION_OBS_VELOCITY_WEIGHT_INV_AREA_LOG,
  RHEA_INVERSION_OBS_VELOCITY_WEIGHT_INV_AREA_LIN
}
rhea_inversion_obs_velocity_weight_t;

/**
 * Generates the observational data.
 */
void                rhea_inversion_obs_velocity_generate (
                        ymir_vec_t *vel_obs_surf,
                        ymir_vec_t *weight_surf,
                        const rhea_inversion_obs_velocity_t obs_type,
                        const rhea_inversion_obs_velocity_weight_t weight_type,
                        rhea_plate_options_t *plate_options,
                        double *calculated_vel_stddev);

/**
 * Compute the misfit of the velocity at the surface:
 *   ObsOp(vel) - vel_obs
 */
void                rhea_inversion_obs_velocity_misfit_vec (
                                  ymir_vec_t *misfit_surf,
                                  ymir_vec_t *vel_fwd_vol,
                                  ymir_vec_t *vel_obs_surf,
                                  ymir_vec_t *weight_surf,
                                  const rhea_inversion_obs_velocity_t obs_type,
                                  rhea_domain_options_t *domain_options);

/**
 * Computes the (squared norm of the) misfit of velocities at the surface.
 */
double              rhea_inversion_obs_velocity_misfit (
                                  ymir_vec_t *vel_fwd_vol,
                                  ymir_vec_t *vel_obs_surf,
                                  ymir_vec_t *weight_surf,
                                  const rhea_inversion_obs_velocity_t obs_type,
                                  rhea_domain_options_t *domain_options);

/**
 * Computes the right-hand side for adjoint equations.
 */
void                rhea_inversion_obs_velocity_adjoint_rhs (
                                  ymir_vec_t *rhs_vel_mass,
                                  ymir_vec_t *vel_fwd_vol,
                                  ymir_vec_t *vel_obs_surf,
                                  ymir_vec_t *weight_surf,
                                  const rhea_inversion_obs_velocity_t obs_type,
                                  rhea_domain_options_t *domain_options);

/**
 * Computes the right-hand side (component) for incremental adjoint equations.
 */
void                rhea_inversion_obs_velocity_incremental_adjoint_rhs (
                                  ymir_vec_t *rhs_vel_mass,
                                  ymir_vec_t *vel_incrfwd_vol,
                                  const rhea_inversion_obs_velocity_t obs_type,
                                  rhea_domain_options_t *domain_options);

/**
 * Removes normal/tangential components of a surface velocity field.
 */
void                rhea_inversion_obs_velocity_remove_normal (
                                  ymir_vec_t * vel_surf);
void                rhea_inversion_obs_velocity_remove_tangential (
                                  ymir_vec_t * vel_surf);

#endif /* RHEA_INVERSION_OBS_VELOCITY_H */
