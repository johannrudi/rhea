/*
 */

#include <rhea_velocity.h>
#include <rhea_base.h>
#include <ymir_vec_getset.h>
#include <ymir_velocity_elem.h>
#include <ymir_interp_vec.h>
#include <ymir_mass_vec.h>

/******************************************************************************
 * Velocity Vector
 *****************************************************************************/

ymir_vec_t *
rhea_velocity_new (ymir_mesh_t *ymir_mesh)
{
  return ymir_cvec_new (ymir_mesh, 3);
}

void
rhea_velocity_destroy (ymir_vec_t *velocity)
{
  ymir_vec_destroy (velocity);
}

void
rhea_velocity_convert_to_dimensional (ymir_vec_t * velocity,
                                      rhea_domain_options_t *domain_options,
                                      rhea_temperature_options_t *temp_options)
{
  const double        dim_scal = temp_options->thermal_diffusivity_m2_s /
                                 domain_options->radius_max_m;

  ymir_vec_scale (dim_scal, velocity);
}

int
rhea_velocity_check_vec_type (ymir_vec_t *vec)
{
  return (
      ymir_vec_is_cvec (vec) &&
      vec->ncfields == 3 &&
      vec->node_type == YMIR_GLL_NODE
  );
}

int
rhea_velocity_is_valid (ymir_vec_t *vec)
{
  return sc_dmatrix_is_valid (vec->dataown) && sc_dmatrix_is_valid (vec->coff);
}

ymir_vec_t *
rhea_velocity_surface_new (ymir_mesh_t *ymir_mesh)
{
  return ymir_face_cvec_new (ymir_mesh, RHEA_DOMAIN_BOUNDARY_FACE_TOP, 3);
}

void
rhea_velocity_surface_destroy (ymir_vec_t *velocity)
{
  ymir_vec_destroy (velocity);
}

int
rhea_velocity_surface_check_vec_type (ymir_vec_t *vec)
{
  return (
      rhea_velocity_check_vec_type (vec) &&
      vec->meshnum == RHEA_DOMAIN_BOUNDARY_FACE_TOP
  );
}

int
rhea_velocity_surface_is_valid (ymir_vec_t *vec)
{
  return rhea_velocity_is_valid (vec);
}

ymir_vec_t *
rhea_velocity_surface_new_from_vol (ymir_vec_t *vel_vol)
{
  ymir_mesh_t        *ymir_mesh;
  ymir_vec_t         *vel_surf;

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (vel_vol));

  /* create surface vector */
  ymir_mesh = ymir_vec_get_mesh (vel_vol);
  vel_surf = rhea_velocity_surface_new (ymir_mesh);

  /* interpolate velocity from volume to surface */
  rhea_velocity_surface_interpolate (vel_surf, vel_vol);
  return vel_surf;
}

void
rhea_velocity_surface_interpolate (ymir_vec_t *vel_surf, ymir_vec_t *vel_vol)
{
  /* check input */
  RHEA_ASSERT (rhea_velocity_surface_check_vec_type (vel_surf));
  RHEA_ASSERT (rhea_velocity_check_vec_type (vel_vol));
  RHEA_ASSERT (rhea_velocity_is_valid (vel_vol));

  /* interpolate */
  ymir_interp_vec (vel_vol, vel_surf);
  RHEA_ASSERT (rhea_velocity_is_valid (vel_surf));
}

/******************************************************************************
 * Get & Set Values
 *****************************************************************************/

void
rhea_velocity_get_elem_gll (sc_dmatrix_t *vel_el_mat,
                            ymir_vec_t *vel_vec,
                            const ymir_locidx_t elid)
{
#ifdef RHEA_ENABLE_DEBUG
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (vel_vec);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (vel_vec));
  RHEA_ASSERT (vel_el_mat->m == n_nodes_per_el);
  RHEA_ASSERT (vel_el_mat->n == 3);
#endif

  /* interpolate from continuous nodes to discontinuous GLL nodes */
  ymir_cvec_get_elem_interp (vel_vec, vel_el_mat, YMIR_STRIDE_NODE,
                             elid, YMIR_GLL_NODE, YMIR_READ);
}

/******************************************************************************
 * Right-Hand Side Computation
 *****************************************************************************/

void
rhea_velocity_rhs_compute (ymir_vec_t *rhs_vel, ymir_vec_t *temperature,
                           void *data)
{
  rhea_temperature_options_t *temp_options = data;

  rhea_temperature_compute_rhs_vel (rhs_vel, temperature, temp_options);
}

/******************************************************************************
 * Statistics
 *****************************************************************************/

static void
rhea_velocity_stats_compute_magnitude (ymir_vec_t *magnitude,
                                       ymir_vec_t *velocity)
{
  ymir_cvec_innerprod_pointwise (magnitude, velocity, velocity);
  ymir_cvec_sqrt (magnitude, magnitude);
}

static double
rhea_velocity_stats_compute_mean (ymir_vec_t *velocity,
                                  rhea_domain_options_t *domain_options)
{
  ymir_vec_t         *velocity_mass = ymir_vec_template (velocity);
  ymir_vec_t         *unit = ymir_vec_template (velocity);
  double              mean;

  ymir_cvec_set_value (unit, 1.0);
  ymir_mass_apply (velocity, velocity_mass);
  mean = ymir_cvec_innerprod (velocity_mass, unit) / domain_options->volume;

  ymir_vec_destroy (velocity_mass);
  ymir_vec_destroy (unit);

  return mean;
}

void
rhea_velocity_stats_get_global (double *magn_cm_yr_min, double *magn_cm_yr_max,
                                double *magn_cm_yr_mean, ymir_vec_t *velocity,
                                rhea_domain_options_t *domain_options,
                                rhea_temperature_options_t *temp_options)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (velocity);
  ymir_vec_t         *magnitude_cm_yr = ymir_cvec_new (ymir_mesh, 1);

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (velocity));
  RHEA_ASSERT (rhea_velocity_is_valid (velocity));

  /* compute pointwise magnitude */
  rhea_velocity_stats_compute_magnitude (magnitude_cm_yr, velocity);

  /* convert to dimensional values [cm/yr] */
  rhea_velocity_convert_to_dimensional (magnitude_cm_yr, domain_options,
                                        temp_options);
  ymir_cvec_scale (100.0, magnitude_cm_yr);

  /* find global values */
  if (magn_cm_yr_min != NULL) {
    *magn_cm_yr_min = ymir_vec_min_global (magnitude_cm_yr);
  }
  if (magn_cm_yr_max != NULL) {
    *magn_cm_yr_max = ymir_vec_max_global (magnitude_cm_yr);
  }
  if (magn_cm_yr_mean != NULL) {
    *magn_cm_yr_mean = rhea_velocity_stats_compute_mean (magnitude_cm_yr,
                                                         domain_options);
  }

  /* destroy */
  ymir_vec_destroy (magnitude_cm_yr);
}

void
rhea_velocity_stats_get_global_surface (
                                double *magn_cm_yr_min, double *magn_cm_yr_max,
                                double *magn_cm_yr_mean, ymir_vec_t *vel_vol,
                                rhea_domain_options_t *domain_options,
                                rhea_temperature_options_t *temp_options)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vel_vol);
  ymir_vec_t         *vel_surf, *magn_surf_cm_yr;

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (vel_vol));
  RHEA_ASSERT (rhea_velocity_is_valid (vel_vol));

  /* compute pointwise magnitude at surface */
  RHEA_INFO ("###DEV### into rhea_velocity_surface_new_from_vol\n");
  vel_surf = rhea_velocity_surface_new_from_vol (vel_vol);
  magn_surf_cm_yr = ymir_face_cvec_new (ymir_mesh,
                                        RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
  RHEA_INFO ("###DEV### into rhea_velocity_stats_compute_magnitude\n");
  rhea_velocity_stats_compute_magnitude (magn_surf_cm_yr, vel_surf);
  ymir_vec_destroy (vel_surf);

  /* convert to dimensional values [cm/yr] */
  RHEA_INFO ("###DEV### into rhea_velocity_convert_to_dimensional\n");
  rhea_velocity_convert_to_dimensional (magn_surf_cm_yr, domain_options,
                                        temp_options);
  ymir_cvec_scale (100.0, magn_surf_cm_yr);

  /* find global values */
  RHEA_INFO ("###DEV### into find global\n");
  if (magn_cm_yr_min != NULL) {
    *magn_cm_yr_min = ymir_vec_min_global (magn_surf_cm_yr);
  }
  if (magn_cm_yr_max != NULL) {
    *magn_cm_yr_max = ymir_vec_max_global (magn_surf_cm_yr);
  }
  if (magn_cm_yr_mean != NULL) {
    *magn_cm_yr_mean = rhea_velocity_stats_compute_mean (magn_surf_cm_yr,
                                                         domain_options);
  }

  /* destroy */
  ymir_vec_destroy (magn_surf_cm_yr);
}
