/*
 */

#include <rhea_stress.h>
#include <rhea_base.h>
#include <rhea_strainrate.h>
#include <rhea_velocity.h>
#include <ymir_mass_vec.h>

/******************************************************************************
 * Viscous Stress Vector
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

static double
rhea_stress_get_dim_scal (rhea_domain_options_t *domain_options,
                          rhea_temperature_options_t *temp_options,
                          rhea_viscosity_options_t *visc_options)
{
  return visc_options->representative_Pas *
         temp_options->thermal_diffusivity_m2_s /
         (domain_options->radius_max_m * domain_options->radius_max_m);
}

void
rhea_stress_convert_to_dimensional_Pa (ymir_vec_t * stress,
                                       rhea_domain_options_t *domain_options,
                                       rhea_temperature_options_t *temp_options,
                                       rhea_viscosity_options_t *visc_options)
{
  ymir_vec_scale (rhea_stress_get_dim_scal (domain_options, temp_options,
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
  RHEA_ASSERT (rhea_stress_2inv_check_vec_type (viscstress_sqrt_2inv));
  RHEA_ASSERT (rhea_strainrate_2inv_check_vec_type (viscstress_sqrt_2inv));
  RHEA_ASSERT (rhea_viscosity_check_vec_type (viscstress_sqrt_2inv));

  ymir_dvec_copy (strainrate_sqrt_2inv, viscstress_sqrt_2inv);
  ymir_dvec_multiply_in (viscosity, viscstress_sqrt_2inv);
  ymir_dvec_scale (2.0, viscstress_sqrt_2inv);
}

double
rhea_stress_compute_norm (ymir_vec_t *stress)
{
  ymir_vec_t         *stress_mass = ymir_vec_clone (stress);
  double              ip;

  RHEA_ASSERT (rhea_stress_check_vec_type (stress));
  RHEA_ASSERT (rhea_stress_is_valid (stress));

  ymir_mass_apply_gauss (stress_mass);
  ip = ymir_dvec_innerprod (stress_mass, stress);
  ymir_vec_destroy (stress_mass);

  return sqrt (ip);
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
  const double        dim_scal = rhea_stress_get_dim_scal (domain_options,
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
