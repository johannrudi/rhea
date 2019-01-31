#include <rhea_inversion_obs_velocity.h>
#include <rhea_base.h>
#include <rhea_velocity.h>
#include <ymir_velocity_vec.h>
#include <ymir_mass_vec.h>

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

static void
rhea_inversion_obs_velocity_remove_normal_fn (double *vec,
                                       double X, double Y, double Z,
                                       double nx, double ny, double nz,
                                       ymir_topidx_t face,
                                       ymir_locidx_t node_id, void *data)
{
  double              vr;

  /* remove the normal component of the face vector */
  vr = nx * vec[0] + ny * vec[1] + nz * vec[2];
  vec[0] -= vr * nx;
  vec[1] -= vr * ny;
  vec[2] -= vr * nz;
}

static void
rhea_inversion_obs_velocity_remove_normal (ymir_vec_t * vel_surf)
{
  RHEA_ASSERT (rhea_inversion_obs_velocity_weight_check_vec_type (vel_surf));

  ymir_face_cvec_set_function (
      vel_surf, rhea_inversion_obs_velocity_remove_normal_fn, NULL);
}

/**
 * Compute the misfit of the velocity at the surface:
 *   ObsOp(vel) - vel_obs
 */
static void
rhea_inversion_obs_velocity_misfit_vec (
                                  ymir_vec_t *misfit_surf,
                                  ymir_vec_t *vel_fwd_vol,
                                  ymir_vec_t *vel_obs_surf,
                                  ymir_vec_t *weight_surf,
                                  const rhea_inversion_obs_velocity_t obs_type,
                                  rhea_domain_options_t *domain_options)
{
  /* check input */
  RHEA_ASSERT (rhea_velocity_surface_check_vec_type (misfit_surf));
  RHEA_ASSERT (rhea_velocity_check_vec_type (vel_fwd_vol));
  RHEA_ASSERT (rhea_velocity_is_valid (vel_fwd_vol));
  RHEA_ASSERT (rhea_velocity_surface_check_vec_type (vel_obs_surf));
  RHEA_ASSERT (rhea_velocity_surface_is_valid (vel_obs_surf));
  RHEA_ASSERT (weight_surf == NULL ||
               rhea_inversion_obs_velocity_weight_check_vec_type (weight_surf));
  RHEA_ASSERT (weight_surf == NULL ||
               rhea_inversion_obs_velocity_weight_is_valid (weight_surf));

  /* project velocity from volume to surface */
  rhea_velocity_surface_interpolate (misfit_surf, vel_fwd_vol);

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
        misfit_surf, domain_options->center,
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
    ymir_vec_add (-1.0, vel_obs_surf, misfit_surf);
    rhea_inversion_obs_velocity_remove_normal (misfit_surf);
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
    if (1 == weight_surf->ncfields) {
      ymir_vec_multiply_in1 (weight_surf, misfit_surf);
    }
    else {
      RHEA_ASSERT (3 == weight_surf->ncfields);
      ymir_vec_multiply_in (weight_surf, misfit_surf);
    }
  }
}

double
rhea_inversion_obs_velocity_misfit (
                                  ymir_vec_t *vel_fwd_vol,
                                  ymir_vec_t *vel_obs_surf,
                                  ymir_vec_t *weight_surf,
                                  const rhea_inversion_obs_velocity_t obs_type,
                                  rhea_domain_options_t *domain_options)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vel_fwd_vol);
  ymir_vec_t         *misfit_surf = rhea_velocity_surface_new (ymir_mesh);
  ymir_vec_t         *misfit_surf_mass = rhea_velocity_surface_new (ymir_mesh);
  double              misfit_norm_sq;

  /* compute misfit vector */
  rhea_inversion_obs_velocity_misfit_vec (
      misfit_surf, vel_fwd_vol, vel_obs_surf, weight_surf, obs_type,
      domain_options);

  /* compute squared L2-norm */
  ymir_mass_apply (misfit_surf, misfit_surf_mass);
  misfit_norm_sq = ymir_vec_innerprod (misfit_surf, misfit_surf_mass);

  /* destroy */
  rhea_velocity_surface_destroy (misfit_surf);
  rhea_velocity_surface_destroy (misfit_surf_mass);

  /* return squared L2-norm of misfit */
  return misfit_norm_sq;
}

void
rhea_inversion_obs_velocity_adjoint_rhs (
                                  ymir_vec_t *rhs_vel_mass,
                                  ymir_vec_t *vel_fwd_vol,
                                  ymir_vec_t *vel_obs_surf,
                                  ymir_vec_t *weight_surf,
                                  const rhea_inversion_obs_velocity_t obs_type,
                                  rhea_domain_options_t *domain_options)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vel_fwd_vol);
  ymir_vec_t         *misfit_surf = rhea_velocity_surface_new (ymir_mesh);

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (vel_fwd_vol));
  RHEA_ASSERT (rhea_velocity_is_valid (vel_fwd_vol));
  RHEA_ASSERT (rhea_velocity_surface_check_vec_type (vel_obs_surf));
  RHEA_ASSERT (rhea_velocity_surface_is_valid (vel_obs_surf));
  RHEA_ASSERT (weight_surf == NULL ||
               rhea_inversion_obs_velocity_weight_check_vec_type (weight_surf));
  RHEA_ASSERT (weight_surf == NULL ||
               rhea_inversion_obs_velocity_weight_is_valid (weight_surf));

  /* compute misfit vector */
  rhea_inversion_obs_velocity_misfit_vec (
      misfit_surf, vel_fwd_vol, vel_obs_surf, weight_surf, obs_type,
      domain_options);

  /* project misfit into volume */
  rhea_velocity_interpolate_from_surface (rhs_vel_mass, misfit_surf,
                                          1 /* mass_weighted */);

  /* change sign to obtain the right-hand side */
  ymir_vec_scale (-1.0, rhs_vel_mass);

  /* destroy */
  rhea_velocity_surface_destroy (misfit_surf);
}
