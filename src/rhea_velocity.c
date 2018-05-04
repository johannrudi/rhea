/*
 */

#include <rhea_velocity.h>
#include <rhea_base.h>
#include <rhea_viscosity.h>
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

static double
rhea_velocity_get_dim_scal (rhea_domain_options_t *domain_options,
                            rhea_temperature_options_t *temp_options)
{
  return temp_options->thermal_diffusivity_m2_s /
         domain_options->radius_max_m;
}

void
rhea_velocity_convert_to_dimensional_m_s (
                                      ymir_vec_t * velocity,
                                      rhea_domain_options_t *domain_options,
                                      rhea_temperature_options_t *temp_options)
{
  ymir_vec_scale (rhea_velocity_get_dim_scal (domain_options, temp_options),
                  velocity);
}

void
rhea_velocity_convert_to_dimensional_cm_yr (
                                      ymir_vec_t * velocity,
                                      rhea_domain_options_t *domain_options,
                                      rhea_temperature_options_t *temp_options)
{
  ymir_vec_scale (rhea_velocity_get_dim_scal (domain_options, temp_options) *
                  100.0 * RHEA_TEMPERATURE_SECONDS_PER_YEAR, velocity);
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

MPI_Offset *
rhea_velocity_segment_offset_create (ymir_vec_t *vec)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vec);
  const mangll_locidx_t *n_nodes = ymir_mesh->cnodes->Ngo;
  const int           n_fields = vec->ncfields;
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  int                 mpisize, mpiret;
  MPI_Offset         *segment_offset;
  int                 r;

  /* get parallel environment */
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize); SC_CHECK_MPI (mpiret);

  /* create segment offsets */
  segment_offset = RHEA_ALLOC (MPI_Offset, mpisize + 1);
  segment_offset[0] = 0;
  for (r = 0; r < mpisize; r++) {
    segment_offset[r+1] = segment_offset[r] +
                          (MPI_Offset) (n_fields * n_nodes[r]);
  }

  return segment_offset;
}

MPI_Offset
rhea_velocity_segment_offset_get (ymir_vec_t *vec)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vec);
  const mangll_locidx_t *n_nodes = ymir_mesh->cnodes->Ngo;
  const int           n_fields = vec->ncfields;
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  int                 mpirank, mpiret;
  MPI_Offset          segment_offset;
  int                 r;

  /* get parallel environment */
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);

  /* add up segment offsets */
  segment_offset = 0;
  for (r = 0; r < mpirank; r++) {
    segment_offset += (MPI_Offset) (n_fields * n_nodes[r]);
  }

  return segment_offset;
}

int
rhea_velocity_segment_size_get (ymir_vec_t *vec)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vec);
  const mangll_locidx_t *n_nodes = ymir_mesh->cnodes->Ngo;
  const int           n_fields = vec->ncfields;
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  int                 mpirank, mpiret;

  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);
  return (int) (n_fields * n_nodes[mpirank]);
}

int
rhea_velocity_is_valid (ymir_vec_t *vec)
{
  return sc_dmatrix_is_valid (vec->dataown) && sc_dmatrix_is_valid (vec->coff);
}

/******************************************************************************
 * Velocity Surface Vector
 *****************************************************************************/

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
rhea_velocity_stats_compute_mean (ymir_vec_t *vec, ymir_vec_t *filter,
                                  rhea_domain_options_t *domain_options)
{
  ymir_vec_t         *vec_mass = ymir_vec_template (vec);
  ymir_vec_t         *unit = ymir_vec_template (vec);
  double              volume, mean;

  /* check input */
  RHEA_ASSERT (ymir_vec_is_cvec (vec));
  RHEA_ASSERT (filter == NULL || ymir_vec_is_cvec (filter));
  RHEA_ASSERT (filter == NULL || filter->ncfields == 1);

  /* set unit vector; set/calculate volume */
  ymir_cvec_set_value (unit, 1.0);
  if (filter != NULL) {
    ymir_mass_apply (unit, vec_mass);
    ymir_cvec_multiply_in1 (filter, vec_mass);
    volume = ymir_cvec_innerprod (vec_mass, unit);
  }
  else if (domain_options == NULL) {
    ymir_mass_apply (unit, vec_mass);
    volume = ymir_cvec_innerprod (vec_mass, unit);
  }
  else {
    volume = domain_options->volume;
  }

  /* compute mean */
  ymir_mass_apply (vec, vec_mass);
  if (filter != NULL) {
    ymir_cvec_multiply_in1 (filter, vec_mass);
  }
  mean = ymir_cvec_innerprod (vec_mass, unit) / volume;

  /* destroy */
  ymir_vec_destroy (vec_mass);
  ymir_vec_destroy (unit);

  return mean;
}

void
rhea_velocity_stats_get_global (double *magn_min_cm_yr, double *magn_max_cm_yr,
                                double *magn_mean_cm_yr, ymir_vec_t *velocity,
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
  rhea_velocity_convert_to_dimensional_cm_yr (magnitude_cm_yr, domain_options,
                                              temp_options);

  /* find global values */
  if (magn_min_cm_yr != NULL) {
    *magn_min_cm_yr = ymir_vec_min_global (magnitude_cm_yr);
  }
  if (magn_max_cm_yr != NULL) {
    *magn_max_cm_yr = ymir_vec_max_global (magnitude_cm_yr);
  }
  if (magn_mean_cm_yr != NULL) {
    *magn_mean_cm_yr = rhea_velocity_stats_compute_mean (magnitude_cm_yr, NULL,
                                                         domain_options);
  }

  /* destroy */
  ymir_vec_destroy (magnitude_cm_yr);
}

void
rhea_velocity_stats_get_global_lithosphere (
                                    double *magn_max_cm_yr,
                                    double *magn_mean_cm_yr,
                                    ymir_vec_t *velocity,
                                    ymir_vec_t *viscosity,
                                    rhea_domain_options_t *domain_options,
                                    rhea_temperature_options_t *temp_options)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (velocity);
  ymir_vec_t         *filter = ymir_cvec_new (ymir_mesh, 1);
  ymir_vec_t         *magnitude_cm_yr = ymir_cvec_new (ymir_mesh, 1);

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (velocity));
  RHEA_ASSERT (rhea_velocity_is_valid (velocity));

  /* compute pointwise magnitude */
  rhea_velocity_stats_compute_magnitude (magnitude_cm_yr, velocity);

  /* filter */
  rhea_viscosity_stats_filter_lithosphere (filter, viscosity, NULL, NAN);
  ymir_vec_multiply_in1 (filter, magnitude_cm_yr);

  /* convert to dimensional values [cm/yr] */
  rhea_velocity_convert_to_dimensional_cm_yr (magnitude_cm_yr, domain_options,
                                              temp_options);

  /* find global values */
  if (magn_max_cm_yr != NULL) {
    *magn_max_cm_yr = ymir_vec_max_global (magnitude_cm_yr);
  }
  if (magn_mean_cm_yr != NULL) {
    *magn_mean_cm_yr = rhea_velocity_stats_compute_mean (magnitude_cm_yr,
                                                         filter,
                                                         domain_options);
  }

  /* destroy */
  ymir_vec_destroy (filter);
  ymir_vec_destroy (magnitude_cm_yr);
}

void
rhea_velocity_stats_get_global_surface (
                                double *magn_min_cm_yr, double *magn_max_cm_yr,
                                double *magn_mean_cm_yr, ymir_vec_t *vel_vol,
                                rhea_domain_options_t *domain_options,
                                rhea_temperature_options_t *temp_options)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vel_vol);
  ymir_vec_t         *magn_vol_cm_yr = ymir_cvec_new (ymir_mesh, 1);
  ymir_vec_t         *magn_surf_cm_yr = ymir_face_cvec_new (
                          ymir_mesh, RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (vel_vol));
  RHEA_ASSERT (rhea_velocity_is_valid (vel_vol));

  /* compute pointwise magnitude */
  rhea_velocity_stats_compute_magnitude (magn_vol_cm_yr, vel_vol);

  /* map magnitude to surface */
  ymir_interp_vec (magn_vol_cm_yr, magn_surf_cm_yr);
  RHEA_ASSERT (rhea_velocity_is_valid (magn_surf_cm_yr));

  /* convert to dimensional values [cm/yr] */
  rhea_velocity_convert_to_dimensional_cm_yr (magn_surf_cm_yr, domain_options,
                                              temp_options);

  /* find global values */
  if (magn_min_cm_yr != NULL) {
    *magn_min_cm_yr = ymir_vec_min_global (magn_surf_cm_yr);
  }
  if (magn_max_cm_yr != NULL) {
    *magn_max_cm_yr = ymir_vec_max_global (magn_surf_cm_yr);
  }
  if (magn_mean_cm_yr != NULL) {
    *magn_mean_cm_yr = rhea_velocity_stats_compute_mean (magn_surf_cm_yr,
                                                         NULL, NULL);
  }

  /* destroy */
  ymir_vec_destroy (magn_vol_cm_yr);
  ymir_vec_destroy (magn_surf_cm_yr);
}

void
rhea_velocity_stats_get_global_surface_lithosphere (
                                    double *magn_max_cm_yr,
                                    double *magn_mean_cm_yr,
                                    ymir_vec_t *vel_vol,
                                    ymir_vec_t *visc_vol,
                                    rhea_domain_options_t *domain_options,
                                    rhea_temperature_options_t *temp_options)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vel_vol);
  ymir_vec_t         *filter_vol = ymir_cvec_new (ymir_mesh, 1);
  ymir_vec_t         *filter_surf = ymir_face_cvec_new (
                          ymir_mesh, RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
  ymir_vec_t         *magn_vol_cm_yr = ymir_cvec_new (ymir_mesh, 1);
  ymir_vec_t         *magn_surf_cm_yr = ymir_face_cvec_new (
                          ymir_mesh, RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (vel_vol));
  RHEA_ASSERT (rhea_velocity_is_valid (vel_vol));

  /* compute pointwise magnitude */
  rhea_velocity_stats_compute_magnitude (magn_vol_cm_yr, vel_vol);

  /* filter in volume */
  rhea_viscosity_stats_filter_lithosphere (filter_vol, visc_vol, NULL, NAN);
  ymir_vec_multiply_in1 (filter_vol, magn_vol_cm_yr);

  /* map magnitude to surface */
  ymir_interp_vec (magn_vol_cm_yr, magn_surf_cm_yr);
  RHEA_ASSERT (rhea_velocity_is_valid (magn_surf_cm_yr));

  /* filter at surface */
  rhea_viscosity_stats_filter_lithosphere_surf (filter_surf, visc_vol, NULL,
                                                NAN);
  ymir_vec_multiply_in1 (filter_surf, magn_surf_cm_yr);

  /* convert to dimensional values [cm/yr] */
  rhea_velocity_convert_to_dimensional_cm_yr (magn_surf_cm_yr, domain_options,
                                              temp_options);

  /* find global values */
  if (magn_max_cm_yr != NULL) {
    *magn_max_cm_yr = ymir_vec_max_global (magn_surf_cm_yr);
  }
  if (magn_mean_cm_yr != NULL) {
    *magn_mean_cm_yr = rhea_velocity_stats_compute_mean (magn_surf_cm_yr,
                                                         filter_surf, NULL);
  }

  /* destroy */
  ymir_vec_destroy (filter_vol);
  ymir_vec_destroy (filter_surf);
  ymir_vec_destroy (magn_vol_cm_yr);
  ymir_vec_destroy (magn_surf_cm_yr);
}
