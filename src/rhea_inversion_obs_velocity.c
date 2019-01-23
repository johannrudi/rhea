#include <rhea_inversion_obs_velocity.h>
#include <rhea_base.h>
#include <rhea_velocity.h>
#include <rhea_velocity_pressure.h>
#include <ymir_velocity_vec.h>

static int
rhea_inversion_obs_velocity_weight_check_vec_type (ymir_vec_t *vec)
{
  if (1 == vec->ncfields) {
    return (
        ymir_vec_is_cvec (vec) &&
        vec->node_type == YMIR_GLL_NODE &&
        ymir_vec_is_face_vec (vec) &&
        vec->meshnum == RHEA_DOMAIN_BOUNDARY_FACE_TOP
    );
  }
  else {
    return rhea_velocity_surface_check_vec_type (vec);
  }
}

static int
rhea_inversion_obs_velocity_weight_is_valid (ymir_vec_t *vec)
{
  if (1 == vec->ncfields) {
    return (sc_dmatrix_is_valid (vec->dataown) &&
            sc_dmatrix_is_valid (vec->coff));
  }
  else {
    return rhea_velocity_surface_is_valid (vec);
  }
}

double
rhea_inversion_obs_velocity_misfit (ymir_vec_t *vel_fwd_vol,
                                    ymir_vec_t *vel_obs_surf,
                                    ymir_vec_t *weight_surf,
                                    rhea_inversion_obs_velocity_t obs_type,
                                    rhea_domain_options_t *domain_options)
{
  ymir_vec_t         *vel_fwd_surf =
                        rhea_velocity_surface_new_from_vol (vel_fwd_vol);
  double              misfit = NAN;

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (vel_fwd_vol));
  RHEA_ASSERT (rhea_velocity_is_valid (vel_fwd_vol));
  RHEA_ASSERT (rhea_velocity_surface_check_vec_type (vel_obs_surf));
  RHEA_ASSERT (rhea_velocity_surface_is_valid (vel_obs_surf));
  RHEA_ASSERT (weight_surf == NULL ||
               rhea_inversion_obs_velocity_weight_check_vec_type (weight_surf));
  RHEA_ASSERT (weight_surf == NULL ||
               rhea_inversion_obs_velocity_weight_is_valid (weight_surf));

  /* project out mean rotation */
  switch (obs_type) {
  case RHEA_INVERSION_OBS_VELOCITY_NORMAL:
  case RHEA_INVERSION_OBS_VELOCITY_TANGENTIAL:
  case RHEA_INVERSION_OBS_VELOCITY_ALL:
    break;
  case RHEA_INVERSION_OBS_VELOCITY_TANGENTIAL_ROTFREE:
  case RHEA_INVERSION_OBS_VELOCITY_ALL_ROTFREE:
    RHEA_ASSERT (domain_options != NULL);
    ymir_velocity_vec_project_out_mean_rotation (
        vel_fwd_surf, domain_options->center,
        domain_options->moment_of_inertia_surface, 0 /* !residual_space */);
    break;
  default: /* unknown observation type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* compute difference between forward and observed velocities */
  switch (obs_type) {
  case RHEA_INVERSION_OBS_VELOCITY_NORMAL:
    RHEA_ABORT_NOT_REACHED (); //TODO
    break;
  case RHEA_INVERSION_OBS_VELOCITY_TANGENTIAL:
  case RHEA_INVERSION_OBS_VELOCITY_TANGENTIAL_ROTFREE:
    //TODO
    break;
  case RHEA_INVERSION_OBS_VELOCITY_ALL:
  case RHEA_INVERSION_OBS_VELOCITY_ALL_ROTFREE:
    RHEA_ABORT_NOT_REACHED (); //TODO
    break;
  default: /* unknown observation type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* apply weight to difference */
  if (weight_surf != NULL) {
    //TODO
  }

  /* compute L2-norm */
  //TODO

  /* destroy */
  rhea_velocity_surface_destroy (vel_fwd_surf);

  /* return misfit value */
  return misfit;
}

void
rhea_inversion_obs_velocity_adjoint_rhs (
                                      ymir_vec_t *rhs_vel_press,
                                      ymir_vec_t *vel_fwd_vol,
                                      ymir_vec_t *vel_obs_surf,
                                      ymir_vec_t *weight_surf,
                                      rhea_inversion_obs_velocity_t obs_type,
                                      rhea_domain_options_t *domain_options)
{
  /* check input */
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (rhs_vel_press));
  RHEA_ASSERT (rhea_velocity_check_vec_type (vel_fwd_vol));
  RHEA_ASSERT (rhea_velocity_is_valid (vel_fwd_vol));
  RHEA_ASSERT (rhea_velocity_surface_check_vec_type (vel_obs_surf));
  RHEA_ASSERT (rhea_velocity_surface_is_valid (vel_obs_surf));
  RHEA_ASSERT (weight_surf == NULL ||
               rhea_inversion_obs_velocity_weight_check_vec_type (weight_surf));
  RHEA_ASSERT (weight_surf == NULL ||
               rhea_inversion_obs_velocity_weight_is_valid (weight_surf));

  //TODO
}
