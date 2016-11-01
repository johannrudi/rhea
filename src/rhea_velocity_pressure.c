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


