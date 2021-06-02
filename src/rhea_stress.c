#include <rhea_stress.h>
#include <rhea_base.h>
#include <rhea_strainrate.h>
#include <rhea_velocity.h>
#include <rhea_pressure.h>
#include <ymir_interp_vec.h>
#include <ymir_mass_vec.h>

/******************************************************************************
 * Options
 *****************************************************************************/

double
rhea_stress_get_dim_Pa (rhea_domain_options_t *domain_options,
                        rhea_temperature_options_t *temp_options,
                        rhea_viscosity_options_t *visc_options)
{
  return visc_options->representative_Pas *
         temp_options->thermal_diffusivity_m2_s /
         (domain_options->radius_max_m * domain_options->radius_max_m);
}

/******************************************************************************
 * Stress Vector
 *****************************************************************************/

ymir_vec_t *
rhea_stress_new (ymir_mesh_t *ymir_mesh)
{
  return ymir_dvec_new (ymir_mesh, 6, YMIR_GAUSS_NODE);
}

void
rhea_stress_destroy (ymir_vec_t *stress)
{
  ymir_vec_destroy (stress);
}

void
rhea_stress_convert_to_dimensional_Pa (ymir_vec_t * stress,
                                       rhea_domain_options_t *domain_options,
                                       rhea_temperature_options_t *temp_options,
                                       rhea_viscosity_options_t *visc_options)
{
  ymir_vec_scale (rhea_stress_get_dim_Pa (domain_options, temp_options,
                                          visc_options),
                  stress);
}

int
rhea_stress_check_vec_type (ymir_vec_t *vec)
{
  return (
      ymir_vec_is_dvec (vec) &&
      vec->ndfields == 6 &&
      vec->node_type == YMIR_GAUSS_NODE
  );
}

int
rhea_stress_is_valid (ymir_vec_t *vec)
{
  return sc_dmatrix_is_valid (vec->dataown);
}

void
rhea_stress_compute_viscstress (ymir_vec_t *viscstress,
                                ymir_vec_t *strainrate,
                                ymir_vec_t *viscosity)
{
  /* check input */
  RHEA_ASSERT (rhea_stress_check_vec_type (viscstress));
  RHEA_ASSERT (rhea_strainrate_check_vec_type (strainrate));
  RHEA_ASSERT (rhea_strainrate_is_valid (strainrate));
  RHEA_ASSERT (rhea_viscosity_check_vec_type (viscosity));
  RHEA_ASSERT (rhea_viscosity_is_valid (viscosity));

  /* compute viscous stress tensor */
  ymir_dvec_copy (strainrate, viscstress);
  ymir_dvec_multiply_in1 (viscosity, viscstress);
  ymir_dvec_scale (2.0, viscstress);
  RHEA_ASSERT (rhea_stress_is_valid (viscstress));
}

void
rhea_stress_combine_stresses (ymir_vec_t *stress,
                              ymir_vec_t *pressure,
                              ymir_pressure_elem_t *press_elem)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (pressure);
  ymir_vec_t         *press_gauss = ymir_dvec_new (ymir_mesh, 1,
                                                   YMIR_GAUSS_NODE);
  const ymir_locidx_t n_elements = stress->K;
  const int           n_nodes_per_el = stress->Np;
  ymir_locidx_t       elid;
  int                 nodeid;

  /* check input */
  RHEA_ASSERT (rhea_stress_check_vec_type (stress));
  RHEA_ASSERT (rhea_stress_is_valid (stress));
  RHEA_ASSERT (rhea_pressure_check_vec_type (pressure, press_elem));
  RHEA_ASSERT (rhea_pressure_is_valid (pressure));

  /* interpolate pressure onto Gauss nodes */
  ymir_interp_vec (pressure, press_gauss);

  /* combine viscous stress and isotropic stress from negative pressure */
  for (elid = 0; elid < n_elements; elid++) {
    for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
      const double       *P = ymir_dvec_index (press_gauss, elid, nodeid, 0);
      double             *S = ymir_dvec_index (stress, elid, nodeid, 0);

      S[0] -= P[0];
      S[3] -= P[0];
      S[5] -= P[0];
    }
  }
  RHEA_ASSERT (rhea_stress_is_valid (stress));

  /* destroy */
  ymir_vec_destroy (press_gauss);
}

double
rhea_stress_compute_norm (ymir_vec_t *stress)
{
  ymir_vec_t         *stress_mass = ymir_vec_clone (stress);
  double              ip;

  /* check input */
  RHEA_ASSERT (rhea_stress_check_vec_type (stress));
  RHEA_ASSERT (rhea_stress_is_valid (stress));

  /* compute inner product */
  ymir_mass_apply_gauss (stress_mass);
  ip = ymir_dvec_innerprod (stress_mass, stress);
  ymir_vec_destroy (stress_mass);

  /* return norm */
  return sqrt (ip);
}

void
rhea_stress_separate_diag_offdiag (ymir_vec_t *stress_diag,
                                   ymir_vec_t *stress_offdiag,
                                   ymir_vec_t *stress)
{
  const ymir_locidx_t n_elements = stress->K;
  const int           n_nodes_per_el = stress->Np;
  ymir_locidx_t       elid;
  int                 nodeid;

  /* check input */
  RHEA_ASSERT (ymir_vec_is_dvec (stress_diag) &&
               stress_diag->ndfields == 3 &&
               stress_diag->node_type == stress->node_type);
  RHEA_ASSERT (ymir_vec_is_dvec (stress_offdiag) &&
               stress_offdiag->ndfields == 3 &&
               stress_offdiag->node_type == stress->node_type);
  RHEA_ASSERT (rhea_stress_check_vec_type (stress));
  RHEA_ASSERT (rhea_stress_is_valid (stress));

  /* separate at each node; assume: stress tensor is upper-triangular */
  for (elid = 0; elid < n_elements; elid++) {
    for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
      const double       *S = ymir_dvec_index (stress, elid, nodeid, 0);
      double             *D = ymir_dvec_index (stress_diag, elid, nodeid, 0);
      double             *E = ymir_dvec_index (stress_offdiag, elid, nodeid, 0);

      D[0] = S[0];
      D[1] = S[3];
      D[2] = S[5];

      E[0] = S[1];
      E[1] = S[2];
      E[2] = S[4];
    }
  }
}

ymir_vec_t *
rhea_stress_nonsymmetric_new (ymir_mesh_t *ymir_mesh)
{
  return ymir_dvec_new (ymir_mesh, 9, YMIR_GAUSS_NODE);
}

void
rhea_stress_nonsymmetric_destroy (ymir_vec_t *stress)
{
  ymir_vec_destroy (stress);
}

ymir_vec_t *
rhea_stress_normal_new (ymir_mesh_t *ymir_mesh)
{
  return ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
}

void
rhea_stress_normal_destroy (ymir_vec_t *stress_norm)
{
  ymir_vec_destroy (stress_norm);
}

ymir_vec_t *
rhea_stress_tangential_new (ymir_mesh_t *ymir_mesh)
{
  return ymir_dvec_new (ymir_mesh, 3, YMIR_GAUSS_NODE);
}

void
rhea_stress_tangential_destroy (ymir_vec_t *stress_tang)
{
  ymir_vec_destroy (stress_tang);
}

int
rhea_stress_normal_check_vec_type (ymir_vec_t *vec)
{
  return (
      ymir_vec_is_dvec (vec) &&
      vec->ndfields == 1 &&
      vec->node_type == YMIR_GAUSS_NODE
  );
}

int
rhea_stress_normal_is_valid (ymir_vec_t *vec)
{
  return sc_dmatrix_is_valid (vec->dataown);
}

int
rhea_stress_tangential_check_vec_type (ymir_vec_t *vec)
{
  return (
      ymir_vec_is_dvec (vec) &&
      vec->ndfields == 3 &&
      vec->node_type == YMIR_GAUSS_NODE
  );
}

int
rhea_stress_tangential_is_valid (ymir_vec_t *vec)
{
  return sc_dmatrix_is_valid (vec->dataown);
}

void
rhea_stress_normal_compute_normal (ymir_vec_t *stress_normal_normal,
                                   ymir_vec_t *stress,
                                   ymir_vec_t *normal)
{
  const ymir_locidx_t n_elements     = stress->K;
  const int           n_nodes_per_el = stress->Np;
  ymir_locidx_t       elid;
  int                 nodeid;
  double              Sn[3];

  /* check input */
  RHEA_ASSERT (rhea_stress_normal_check_vec_type (stress_normal_normal));
  RHEA_ASSERT (rhea_stress_check_vec_type (stress));
  RHEA_ASSERT (rhea_stress_is_valid (stress));
  RHEA_ASSERT (ymir_vec_is_dvec (normal) &&
               normal->ndfields == 3 &&
               normal->node_type == stress->node_type);
  RHEA_ASSERT (sc_dmatrix_is_valid (normal->dataown));

  /* apply normal and then get its normal component;
     assume: stress tensor is upper-triangular */
  for (elid = 0; elid < n_elements; elid++) {
    for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
      const double *S = ymir_dvec_index (stress, elid, nodeid, 0);
      const double *N = ymir_dvec_index (normal, elid, nodeid, 0);
      double *nSn = ymir_dvec_index (stress_normal_normal, elid, nodeid, 0);

      /* compute Sn = S * N */
      Sn[0] = S[0]*N[0] + S[1]*N[1] + S[2]*N[2];
      Sn[1] = S[1]*N[0] + S[3]*N[1] + S[4]*N[2];
      Sn[2] = S[2]*N[0] + S[4]*N[1] + S[5]*N[2];

      /* compute nSn = N^T * Sn = N^T * S * N */
      nSn[0] = Sn[0]*N[0] + Sn[1]*N[1] + Sn[2]*N[2];
    }
  }
  RHEA_ASSERT (rhea_stress_normal_is_valid (stress_normal_normal));
}

void
rhea_stress_normal_compute_tangential (ymir_vec_t *stress_normal_tangential,
                                       ymir_vec_t *stress,
                                       ymir_vec_t *normal)
{
  const ymir_locidx_t n_elements     = stress->K;
  const int           n_nodes_per_el = stress->Np;
  ymir_locidx_t       elid;
  int                 nodeid;
  double              Sn[3], nSn;

  /* check input */
  RHEA_ASSERT (
      rhea_stress_tangential_check_vec_type (stress_normal_tangential));
  RHEA_ASSERT (rhea_stress_check_vec_type (stress));
  RHEA_ASSERT (rhea_stress_is_valid (stress));
  RHEA_ASSERT (ymir_vec_is_dvec (normal) &&
               normal->ndfields == 3 &&
               normal->node_type == stress->node_type);
  RHEA_ASSERT (sc_dmatrix_is_valid (normal->dataown));

  /* apply normal and then get its tangential component;
     assume: stress tensor is upper-triangular */
  for (elid = 0; elid < n_elements; elid++) {
    for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
      const double *S = ymir_dvec_index (stress, elid, nodeid, 0);
      const double *N = ymir_dvec_index (normal, elid, nodeid, 0);
      double *tSn = ymir_dvec_index (stress_normal_tangential, elid, nodeid, 0);

      /* compute Sn = S * N */
      Sn[0] = S[0]*N[0] + S[1]*N[1] + S[2]*N[2];
      Sn[1] = S[1]*N[0] + S[3]*N[1] + S[4]*N[2];
      Sn[2] = S[2]*N[0] + S[4]*N[1] + S[5]*N[2];

      /* compute nSn = N^T * Sn = N^T * S * N */
      nSn = Sn[0]*N[0] + Sn[1]*N[1] + Sn[2]*N[2];

      /* compute tSn = Sn - nSn*N */
      tSn[0] = Sn[0] - nSn*N[0];
      tSn[1] = Sn[1] - nSn*N[1];
      tSn[2] = Sn[2] - nSn*N[2];
    }
  }
  RHEA_ASSERT (rhea_stress_tangential_is_valid (stress_normal_tangential));
}

ymir_vec_t *
rhea_stress_2inv_new (ymir_mesh_t *ymir_mesh)
{
  return ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
}

void
rhea_stress_2inv_destroy (ymir_vec_t *stress_2inv)
{
  ymir_vec_destroy (stress_2inv);
}

int
rhea_stress_2inv_check_vec_type (ymir_vec_t *vec)
{
  return (
      ymir_vec_is_dvec (vec) &&
      vec->ndfields == 1 &&
      vec->node_type == YMIR_GAUSS_NODE
  );
}

int
rhea_stress_2inv_is_valid (ymir_vec_t *vec)
{
  return (
      sc_dmatrix_is_valid (vec->dataown) &&
      0.0 <= ymir_dvec_min_global (vec)
  );
}

void
rhea_stress_compute_viscstress_sqrt_of_2inv (ymir_vec_t *viscstress_sqrt_2inv,
                                             ymir_vec_t *strainrate_sqrt_2inv,
                                             ymir_vec_t *viscosity)
{
  /* check input */
  RHEA_ASSERT (rhea_stress_2inv_check_vec_type (viscstress_sqrt_2inv));
  RHEA_ASSERT (rhea_strainrate_2inv_check_vec_type (strainrate_sqrt_2inv));
  RHEA_ASSERT (rhea_strainrate_2inv_is_valid (strainrate_sqrt_2inv));
  RHEA_ASSERT (rhea_viscosity_check_vec_type (viscosity));
  RHEA_ASSERT (rhea_viscosity_is_valid (viscosity));

  /* compute sqrt of the 2nd invariant of the viscous stress */
  ymir_dvec_copy (strainrate_sqrt_2inv, viscstress_sqrt_2inv);
  ymir_dvec_multiply_in (viscosity, viscstress_sqrt_2inv);
  ymir_dvec_scale (2.0, viscstress_sqrt_2inv);
  RHEA_ASSERT (rhea_stress_2inv_is_valid (viscstress_sqrt_2inv));
}

/******************************************************************************
 * Stress Surface Vector
 *****************************************************************************/

ymir_vec_t *
rhea_stress_surface_new (ymir_mesh_t *ymir_mesh)
{
  return ymir_face_cvec_new (ymir_mesh, RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
}

void
rhea_stress_surface_destroy (ymir_vec_t *stress)
{
  ymir_vec_destroy (stress);
}

int
rhea_stress_surface_check_vec_type (ymir_vec_t *vec)
{
  return (
      ymir_vec_is_cvec (vec) &&
      vec->ncfields == 1 &&
      vec->node_type == YMIR_GLL_NODE &&
      vec->meshnum == RHEA_DOMAIN_BOUNDARY_FACE_TOP
  );
}

int
rhea_stress_surface_is_valid (ymir_vec_t *vec)
{
  return sc_dmatrix_is_valid (vec->dataown) && sc_dmatrix_is_valid (vec->coff);
}

static void
rhea_stress_surface_extract_normal_fn (double *stress_norm,
                                       double X, double Y, double Z,
                                       double nx, double ny, double nz,
                                       ymir_topidx_t face,
                                       ymir_locidx_t node_id, void *data)
{
  ymir_vec_t         *vec_surf = (ymir_vec_t *) data;
  double             *v = ymir_cvec_index (vec_surf, node_id, 0);

  YMIR_ASSERT (vec_surf->ncfields == 3);

  /* compute inner product with the boundary's outer normal vector */
  *stress_norm = nx * v[0] + ny * v[1] + nz * v[2];
}

void
rhea_stress_surface_extract_from_residual (ymir_vec_t *stress_norm_surf,
                                           ymir_vec_t *residual_mom)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (stress_norm_surf);
  ymir_vec_t         *residual_mom_surf;
  ymir_vec_t         *mass_lump_surf;

  /* check input */
  RHEA_ASSERT (rhea_stress_surface_check_vec_type (stress_norm_surf));
  RHEA_ASSERT (rhea_velocity_check_vec_type (residual_mom));
  RHEA_ASSERT (rhea_velocity_is_valid (residual_mom));

  /* interpolate (velocity component of the) residual onto surface */
  residual_mom_surf = rhea_velocity_surface_new_from_vol (residual_mom);
  RHEA_ASSERT (rhea_velocity_surface_is_valid (residual_mom_surf));

  /* extract the normal component of the residual at the surface */
  RHEA_ASSERT (residual_mom_surf->meshnum == stress_norm_surf->meshnum);
  ymir_face_cvec_set_function (
      stress_norm_surf, rhea_stress_surface_extract_normal_fn,
      residual_mom_surf);
  rhea_velocity_surface_destroy (residual_mom_surf);

  /* remove the surface mass from residual (i.e., invert mass matrix) */
  mass_lump_surf = rhea_stress_surface_new (ymir_mesh);
  ymir_mass_lump (mass_lump_surf);
  ymir_vec_divide_in (mass_lump_surf, stress_norm_surf);
  rhea_stress_surface_destroy (mass_lump_surf);
}

/******************************************************************************
 * Statistics
 *****************************************************************************/

static double
rhea_stress_stats_compute_mean (ymir_vec_t *vec, ymir_vec_t *filter,
                                rhea_domain_options_t *domain_options)
{
  ymir_vec_t         *vec_mass = ymir_vec_template (vec);
  ymir_vec_t         *unit = ymir_vec_template (vec);
  double              volume, mean;

  /* check input */
  RHEA_ASSERT (ymir_vec_is_dvec (vec));
  RHEA_ASSERT (filter == NULL || ymir_vec_is_dvec (filter));
  RHEA_ASSERT (filter == NULL || filter->ndfields == 1);
  RHEA_ASSERT (filter == NULL || filter->node_type == vec->node_type);

  /* set unit vector; set/calculate volume */
  ymir_dvec_set_value (unit, 1.0);
  if (filter != NULL) {
    ymir_mass_apply (unit, vec_mass);
    ymir_dvec_multiply_in1 (filter, vec_mass);
    volume = ymir_dvec_innerprod (vec_mass, unit);
  }
  else if (domain_options == NULL) {
    ymir_mass_apply (unit, vec_mass);
    volume = ymir_dvec_innerprod (vec_mass, unit);
  }
  else {
    volume = domain_options->volume;
  }

  /* compute mean */
  ymir_mass_apply (vec, vec_mass);
  if (filter != NULL) {
    ymir_dvec_multiply_in1 (filter, vec_mass);
  }
  mean = ymir_dvec_innerprod (vec_mass, unit) / volume;

  /* destroy */
  ymir_vec_destroy (vec_mass);
  ymir_vec_destroy (unit);

  return mean;
}

void
rhea_stress_stats_get_global (double *min_Pa, double *max_Pa, double *mean_Pa,
                              ymir_vec_t *velocity, ymir_vec_t *viscosity,
                              rhea_domain_options_t *domain_options,
                              rhea_temperature_options_t *temp_options,
                              rhea_viscosity_options_t *visc_options)
{
  const double        dim_scal = rhea_stress_get_dim_Pa (domain_options,
                                                         temp_options,
                                                         visc_options);
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (velocity);
  ymir_vec_t         *sr_sqrt_2inv = rhea_strainrate_2inv_new (ymir_mesh);
  ymir_vec_t         *vs_sqrt_2inv = rhea_stress_2inv_new (ymir_mesh);

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (velocity));
  RHEA_ASSERT (rhea_velocity_is_valid (velocity));
  RHEA_ASSERT (rhea_viscosity_check_vec_type (viscosity));
  RHEA_ASSERT (rhea_viscosity_is_valid (viscosity));

  /* compute sqrt of the 2nd invariant */
  rhea_strainrate_compute_sqrt_of_2inv (sr_sqrt_2inv, velocity);
  RHEA_ASSERT (rhea_strainrate_2inv_is_valid (sr_sqrt_2inv));
  rhea_stress_compute_viscstress_sqrt_of_2inv (vs_sqrt_2inv, sr_sqrt_2inv,
                                               viscosity);
  RHEA_ASSERT (rhea_stress_2inv_is_valid (vs_sqrt_2inv));

  /* find global values */
  if (min_Pa != NULL) {
    *min_Pa = dim_scal * ymir_vec_min_global (vs_sqrt_2inv);
  }
  if (max_Pa != NULL) {
    *max_Pa = dim_scal * ymir_vec_max_global (vs_sqrt_2inv);
  }
  if (mean_Pa != NULL) {
    *mean_Pa = dim_scal * rhea_stress_stats_compute_mean (vs_sqrt_2inv, NULL,
                                                          domain_options);
  }

  /* destroy */
  rhea_strainrate_2inv_destroy (sr_sqrt_2inv);
  rhea_stress_2inv_destroy (vs_sqrt_2inv);
}

static double
rhea_stress_surface_stats_compute_mean (ymir_vec_t *vec, ymir_vec_t *filter)
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
  else {
    ymir_mass_apply (unit, vec_mass);
    volume = ymir_cvec_innerprod (vec_mass, unit);
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
rhea_stress_surface_stats_get_global (double *min_Pa, double *max_Pa,
                                      double *mean_Pa,
                                      ymir_vec_t *stress_norm_surf,
                                      rhea_domain_options_t *domain_options,
                                      rhea_temperature_options_t *temp_options,
                                      rhea_viscosity_options_t *visc_options)
{
  const double        dim_scal = rhea_stress_get_dim_Pa (domain_options,
                                                         temp_options,
                                                         visc_options);

  /* check input */
  RHEA_ASSERT (rhea_stress_surface_check_vec_type (stress_norm_surf));
  RHEA_ASSERT (rhea_stress_surface_is_valid (stress_norm_surf));

  /* find global values */
  if (min_Pa != NULL) {
    *min_Pa = dim_scal * ymir_vec_min_global (stress_norm_surf);
  }
  if (max_Pa != NULL) {
    *max_Pa = dim_scal * ymir_vec_max_global (stress_norm_surf);
  }
  if (mean_Pa != NULL) {
    *mean_Pa = dim_scal * rhea_stress_surface_stats_compute_mean (
        stress_norm_surf, NULL);
  }
}
