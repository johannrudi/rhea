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
                                              ymir_locidx_t node_id,
                                              void *data)
{
  double              vn;

  /* remove the normal component of the face vector */
  vn = nx * vec[0] + ny * vec[1] + nz * vec[2];
  vec[0] -= vn * nx;
  vec[1] -= vn * ny;
  vec[2] -= vn * nz;
}

static void
rhea_inversion_obs_velocity_remove_tangential_fn (
                                              double *vec,
                                              double X, double Y, double Z,
                                              double nx, double ny, double nz,
                                              ymir_topidx_t face,
                                              ymir_locidx_t node_id,
                                              void *data)
{
  double              vn;

  /* keep only the normal component of the face vector */
  vn = nx * vec[0] + ny * vec[1] + nz * vec[2];
  vec[0] = vn * nx;
  vec[1] = vn * ny;
  vec[2] = vn * nz;
}

void
rhea_inversion_obs_velocity_remove_normal (ymir_vec_t * vel_surf)
{
  RHEA_ASSERT (rhea_inversion_obs_velocity_weight_check_vec_type (vel_surf));

  ymir_face_cvec_set_function (
      vel_surf, rhea_inversion_obs_velocity_remove_normal_fn, NULL);
}

void
rhea_inversion_obs_velocity_remove_tangential (ymir_vec_t * vel_surf)
{
  RHEA_ASSERT (rhea_inversion_obs_velocity_weight_check_vec_type (vel_surf));

  ymir_face_cvec_set_function (
      vel_surf, rhea_inversion_obs_velocity_remove_tangential_fn, NULL);
}

static double
_area_to_weight_sqrt (double plate_area, double total_area)
{
  RHEA_ASSERT (0 < plate_area && plate_area <= total_area);
  return sqrt (total_area/plate_area);
}

static double
_area_to_weight_log (double plate_area, double total_area)
{
  RHEA_ASSERT (0 < plate_area && plate_area <= total_area);
  return 1.0 + log (total_area/plate_area);
}

static double
_area_to_weight_lin (double plate_area, double total_area)
{
  RHEA_ASSERT (0 < plate_area && plate_area <= total_area);
  return total_area/plate_area;
}

void
rhea_inversion_obs_velocity_generate (
                        ymir_vec_t *vel_obs_surf,
                        ymir_vec_t *weight_surf,
                        const rhea_inversion_obs_velocity_t obs_type,
                        const rhea_inversion_obs_velocity_weight_t weight_type,
                        rhea_plate_options_t *plate_options,
                        double *calculated_vel_rel_stddev)
{
  const int           n_plates = rhea_plate_get_n_plates (plate_options);

  /* generate velocity data */
  switch (obs_type) {
  case RHEA_INVERSION_OBS_VELOCITY_NORMAL:
    RHEA_ABORT_NOT_REACHED (); //TODO
    break;
  case RHEA_INVERSION_OBS_VELOCITY_TANGENTIAL:
  case RHEA_INVERSION_OBS_VELOCITY_TANGENTIAL_ROTFREE:
    ymir_vec_set_value (vel_obs_surf, 0.0);
    if (0 < n_plates) { /* if plates exist */
      rhea_plate_velocity_generate_all (vel_obs_surf, plate_options);
    }
    break;
  case RHEA_INVERSION_OBS_VELOCITY_ALL:
  case RHEA_INVERSION_OBS_VELOCITY_ALL_ROTFREE:
    RHEA_ABORT_NOT_REACHED (); //TODO
    break;
  default: /* unknown observation type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* generate weights */
  switch (weight_type) {
  case RHEA_INVERSION_OBS_VELOCITY_WEIGHT_UNIFORM:
    ymir_vec_set_value (weight_surf, 1.0);
    if (0 < n_plates) { /* if plates exist */
      rhea_plate_apply_filter_all_vec (weight_surf, plate_options);
    }
    break;
  case RHEA_INVERSION_OBS_VELOCITY_WEIGHT_INV_AREA_SQRT:
    if (0 < n_plates) { /* if plates exist */
      rhea_plate_set_weight_vec (weight_surf, _area_to_weight_sqrt,
                                 plate_options, calculated_vel_rel_stddev);
    }
    else {
      ymir_vec_set_value (weight_surf, 1.0);
    }
    break;
  case RHEA_INVERSION_OBS_VELOCITY_WEIGHT_INV_AREA_LOG:
    if (0 < n_plates) { /* if plates exist */
      rhea_plate_set_weight_vec (weight_surf, _area_to_weight_log,
                                 plate_options, calculated_vel_rel_stddev);
    }
    else {
      ymir_vec_set_value (weight_surf, 1.0);
    }
    break;
  case RHEA_INVERSION_OBS_VELOCITY_WEIGHT_INV_AREA_LIN:
    if (0 < n_plates) { /* if plates exist */
      rhea_plate_set_weight_vec (weight_surf, _area_to_weight_lin,
                                 plate_options, calculated_vel_rel_stddev);
    }
    else {
      ymir_vec_set_value (weight_surf, 1.0);
    }
    break;
  default: /* unknown weight type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* transform plate weights to standard deviations */
  if (NULL != calculated_vel_rel_stddev) {
    int                 pid;

    if (weight_type ==  RHEA_INVERSION_OBS_VELOCITY_WEIGHT_UNIFORM) {
      /* set same std. dev. for all plates = 1.0 */
      for (pid = 0; pid < n_plates; pid++) {
        calculated_vel_rel_stddev[pid] = 1.0;
      }
    }
    else {
      /* set std. dev. from calculated weights */
      for (pid = 0; pid < n_plates; pid++) {
        RHEA_ASSERT (isfinite (calculated_vel_rel_stddev[pid]) &&
                     0.0 < calculated_vel_rel_stddev[pid]);
        calculated_vel_rel_stddev[pid] =
          1.0 / sqrt(calculated_vel_rel_stddev[pid]);
      }
    }
  }
}

void
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
  RHEA_ASSERT (vel_fwd_vol == NULL ||
               rhea_velocity_check_vec_type (vel_fwd_vol));
  RHEA_ASSERT (vel_fwd_vol == NULL ||
               rhea_velocity_is_valid (vel_fwd_vol));
  RHEA_ASSERT (vel_obs_surf == NULL ||
               rhea_velocity_surface_check_vec_type (vel_obs_surf));
  RHEA_ASSERT (vel_obs_surf == NULL ||
               rhea_velocity_surface_is_valid (vel_obs_surf));
  RHEA_ASSERT (weight_surf == NULL ||
               rhea_inversion_obs_velocity_weight_check_vec_type (weight_surf));
  RHEA_ASSERT (weight_surf == NULL ||
               rhea_inversion_obs_velocity_weight_is_valid (weight_surf));

  /* project velocity from volume to surface */
  if (vel_fwd_vol != NULL) {
    rhea_velocity_surface_interpolate (misfit_surf, vel_fwd_vol);
  }
  else {
    ymir_vec_set_zero (misfit_surf);
  }

  /* compute difference between forward and observed velocities */
  if (vel_obs_surf != NULL) {
    ymir_vec_add (-1.0, vel_obs_surf, misfit_surf);
  }

  /* project out mean rotation (globally) */
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

  /* remove normal/tangential components (pointwise) */
  switch (obs_type) {
  case RHEA_INVERSION_OBS_VELOCITY_NORMAL:
    rhea_inversion_obs_velocity_remove_tangential (misfit_surf);
    break;
  case RHEA_INVERSION_OBS_VELOCITY_TANGENTIAL:
  case RHEA_INVERSION_OBS_VELOCITY_TANGENTIAL_ROTFREE:
    rhea_inversion_obs_velocity_remove_normal (misfit_surf);
    break;
  case RHEA_INVERSION_OBS_VELOCITY_ALL:
  case RHEA_INVERSION_OBS_VELOCITY_ALL_ROTFREE:
    break;
  default: /* unknown observation type */
    RHEA_ABORT_NOT_REACHED ();
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
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vel_obs_surf);
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
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (rhs_vel_mass);
  ymir_vec_t         *misfit_surf = rhea_velocity_surface_new (ymir_mesh);

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (rhs_vel_mass));
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

void
rhea_inversion_obs_velocity_incremental_adjoint_rhs (
                                  ymir_vec_t *rhs_vel_mass,
                                  ymir_vec_t *vel_incrfwd_vol,
                                  const rhea_inversion_obs_velocity_t obs_type,
                                  rhea_domain_options_t *domain_options)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (rhs_vel_mass);
  ymir_vec_t         *vel_incrfwd_surf = rhea_velocity_surface_new (ymir_mesh);

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (rhs_vel_mass));
  RHEA_ASSERT (rhea_velocity_check_vec_type (vel_incrfwd_vol));
  RHEA_ASSERT (rhea_velocity_is_valid (vel_incrfwd_vol));

  /* compute misfit vector */
  rhea_inversion_obs_velocity_misfit_vec (
      vel_incrfwd_surf, vel_incrfwd_vol,
      NULL /* vel_obs_surf */, NULL /* weight_surf */,
      obs_type, domain_options);

  /* project misfit into volume */
  rhea_velocity_interpolate_from_surface (rhs_vel_mass, vel_incrfwd_surf,
                                          1 /* mass_weighted */);

  /* change sign to obtain the right-hand side */
  ymir_vec_scale (-1.0, rhs_vel_mass);

  /* destroy */
  rhea_velocity_surface_destroy (vel_incrfwd_surf);
}
