/*
 */

#include <rhea_pressure.h>
#include <rhea_base.h>
#include <ymir_pressure_vec.h>

ymir_vec_t *
rhea_pressure_new (ymir_mesh_t *ymir_mesh, ymir_pressure_elem_t *press_elem)
{
  return ymir_pressure_vec_new (ymir_mesh, press_elem);
}

void
rhea_pressure_destroy (ymir_vec_t *pressure)
{
  ymir_vec_destroy (pressure);
}

void
rhea_pressure_convert_to_dimensional (ymir_vec_t * pressure,
                                      rhea_domain_options_t *domain_options,
                                      rhea_temperature_options_t *temp_options,
                                      rhea_viscosity_options_t *visc_options)
{
  const double        dim_scal = visc_options->representative_Pas *
                                 temp_options->thermal_diffusivity_m2_s /
                                 domain_options->radius_max_m /
                                 domain_options->radius_max_m;

  ymir_vec_scale (dim_scal, pressure);
}

int
rhea_pressure_check_vec_type (ymir_vec_t *vec, ymir_pressure_elem_t *press_elem)
{
  if (press_elem != NULL) {
    return ymir_pressure_vec_is_press_vec (vec, press_elem);
  }
  else {
    return ymir_pressure_vec_is_press_vec_simple (vec);
  }
}

int
rhea_pressure_is_valid (ymir_vec_t *vec)
{
  return sc_dmatrix_is_valid (vec->dataown);
}
