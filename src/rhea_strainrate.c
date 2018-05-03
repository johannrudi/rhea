/*
 */

#include <rhea_strainrate.h>
#include <rhea_base.h>
#include <rhea_velocity.h>
#include <ymir_velocity_elem.h>
#include <ymir_velocity_vec.h>
#include <ymir_mass_vec.h>

/******************************************************************************
 * Strain Rate Vector
 *****************************************************************************/

ymir_vec_t *
rhea_strainrate_new (ymir_mesh_t *ymir_mesh)
{
  return ymir_dvec_new (ymir_mesh, 6, YMIR_GAUSS_NODE);
}

void
rhea_strainrate_destroy (ymir_vec_t *strainrate)
{
  ymir_vec_destroy (strainrate);
}

static double
rhea_strainrate_get_dim_scal (rhea_domain_options_t *domain_options,
                              rhea_temperature_options_t *temp_options)
{
  return temp_options->thermal_diffusivity_m2_s /
         (domain_options->radius_max_m * domain_options->radius_max_m);
}

void
rhea_strainrate_convert_to_dimensional_1_s (
                                    ymir_vec_t * strainrate,
                                    rhea_domain_options_t *domain_options,
                                    rhea_temperature_options_t *temp_options)
{
  ymir_vec_scale (rhea_strainrate_get_dim_scal (domain_options, temp_options),
                  strainrate);
}

int
rhea_strainrate_check_vec_type (ymir_vec_t *vec)
{
  return (
      ymir_vec_is_dvec (vec) &&
      vec->ndfields == 6 &&
      vec->node_type == YMIR_GAUSS_NODE
  );
}

int
rhea_strainrate_is_valid (ymir_vec_t *vec)
{
  return sc_dmatrix_is_valid (vec->dataown);
}

ymir_vec_t *
rhea_strainrate_2inv_new (ymir_mesh_t *ymir_mesh)
{
  return ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
}

void
rhea_strainrate_2inv_destroy (ymir_vec_t *strainrate_2inv)
{
  ymir_vec_destroy (strainrate_2inv);
}

int
rhea_strainrate_2inv_check_vec_type (ymir_vec_t *vec)
{
  return (
      ymir_vec_is_dvec (vec) &&
      vec->ndfields == 1 &&
      vec->node_type == YMIR_GAUSS_NODE
  );
}

int
rhea_strainrate_2inv_is_valid (ymir_vec_t *vec)
{
  return (
      sc_dmatrix_is_valid (vec->dataown) &&
      0.0 <= ymir_dvec_min_global (vec)
  );
}

void
rhea_strainrate_compute_sqrt_of_2inv (ymir_vec_t *strainrate_sqrt_2inv,
                                      ymir_vec_t *velocity)
{
  ymir_mesh_t           *ymir_mesh = ymir_vec_get_mesh (velocity);
  ymir_velocity_elem_t  *vel_elem;

  /* check input */
  RHEA_ASSERT (rhea_strainrate_2inv_check_vec_type (strainrate_sqrt_2inv));
  RHEA_ASSERT (rhea_velocity_check_vec_type (velocity));

  /* create work variables */
  vel_elem = ymir_velocity_elem_new (ymir_mesh->ma->N, ymir_mesh->ma->ompsize);

  /* compute the square root of the second invariant of the strain rate */
  ymir_second_invariant_vec (velocity, strainrate_sqrt_2inv, vel_elem);
  ymir_dvec_sqrt (strainrate_sqrt_2inv, strainrate_sqrt_2inv);

  /* destroy */
  ymir_velocity_elem_destroy (vel_elem);
}

/******************************************************************************
 * Get & Set Values
 *****************************************************************************/

void
rhea_strainrate_compute_sqrt_of_2inv_elem (
                                        sc_dmatrix_t *strainrate_2inv_el_mat,
                                        sc_dmatrix_t *vel_el_mat,
                                        ymir_vec_t *vel_vec,
                                        const ymir_locidx_t elid,
                                        sc_dmatrix_t *tmp_grad_vel,
                                        sc_dmatrix_t *tmp_dvel,
                                        sc_dmatrix_t *tmp_vel)
{
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (vel_vec);

  /* get velocity at GLL nodes */
  rhea_velocity_get_elem_gll (vel_el_mat, vel_vec, elid);

  /* compute the 2nd invariant of the strain rate at Gauss nodes */
  ymir_velocity_elem_compute_strain_rate_2inv (
      strainrate_2inv_el_mat, vel_el_mat, mesh, elid,
      tmp_grad_vel, tmp_dvel, tmp_vel);
}

/******************************************************************************
 * Statistics
 *****************************************************************************/

static double
rhea_strainrate_stats_compute_mean (ymir_vec_t *vec, ymir_vec_t *filter,
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
rhea_strainrate_stats_get_global (double *min_1_s, double *max_1_s,
                                  double *mean_1_s, ymir_vec_t *velocity,
                                  rhea_domain_options_t *domain_options,
                                  rhea_temperature_options_t *temp_options)
{
  const double        dim_scal = rhea_strainrate_get_dim_scal (domain_options,
                                                               temp_options);
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (velocity);
  ymir_vec_t         *sr_sqrt_2inv = rhea_strainrate_2inv_new (ymir_mesh);

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (velocity));
  RHEA_ASSERT (rhea_velocity_is_valid (velocity));

  /* compute sqrt of the 2nd invariant */
  rhea_strainrate_compute_sqrt_of_2inv (sr_sqrt_2inv, velocity);

  /* find global values */
  if (min_1_s != NULL) {
    *min_1_s = dim_scal * ymir_vec_min_global (sr_sqrt_2inv);
  }
  if (max_1_s != NULL) {
    *max_1_s = dim_scal * ymir_vec_max_global (sr_sqrt_2inv);
  }
  if (mean_1_s != NULL) {
    *mean_1_s = dim_scal * rhea_strainrate_stats_compute_mean (
        sr_sqrt_2inv, NULL, domain_options);
  }

  /* destroy */
  rhea_strainrate_2inv_destroy (sr_sqrt_2inv);
}
