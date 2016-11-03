/*
 */

#include <rhea_velocity_pressure.h>
#include <rhea_base.h>
#include <ymir_stokes_vec.h>

ymir_vec_t *
rhea_velocity_pressure_new (ymir_mesh_t *ymir_mesh,
                            ymir_pressure_elem_t *press_elem)
{
  return ymir_stokes_vec_new (ymir_mesh, press_elem);
}

void
rhea_velocity_pressure_destroy (ymir_vec_t *velocity_pressure)
{
  ymir_vec_destroy (velocity_pressure);
}

int
rhea_velocity_pressure_check_vec_type (ymir_vec_t *vec)
{
  return ymir_stokes_vec_is_stokes_vec (vec);
}

int
rhea_velocity_pressure_is_valid (ymir_vec_t *vec)
{
  return sc_dmatrix_is_valid (vec->dataown) && sc_dmatrix_is_valid (vec->coff);
}
