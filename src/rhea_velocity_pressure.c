/*
 */

#include <rhea_velocity_pressure.h>
#include <rhea_base.h>
#include <ymir_pressure_vec.h>
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

void
rhea_velocity_pressure_get_components (ymir_vec_t **vel, ymir_vec_t **press,
                                       ymir_vec_t *vel_press,
                                       ymir_pressure_elem_t *press_elem)
{
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (vel_press);

  /* check input */
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (vel_press));

  if (press_elem->space == YMIR_PRESSURE_SPACE_STAB) { /* if copy data */
    if (vel != NULL && press != NULL) {
      *vel = ymir_cvec_new (mesh, 3);
      *press = ymir_pressure_vec_new (mesh, press_elem);
      ymir_stokes_vec_get_components (vel_press, *vel, *press, press_elem);
    }
    else if (vel != NULL) {
      *vel = ymir_cvec_new (mesh, 3);
      ymir_stokes_vec_get_velocity (vel_press, *vel, press_elem);
    }
    else if (press != NULL) {
      *press = ymir_pressure_vec_new (mesh, press_elem);
      ymir_stokes_vec_get_pressure (vel_press, *press, press_elem);
    }
  }
  else { /* otherwise create view onto data */
    RHEA_ASSERT (press_elem->space == YMIR_PRESSURE_SPACE_POLY ||
                 press_elem->space == YMIR_PRESSURE_SPACE_TENS);
    if (vel != NULL) {
      *vel = ymir_cvec_new_data (mesh, vel_press->ncfields, vel_press->cvec);
    }
    if (press != NULL) {
      *press = ymir_evec_new_data (mesh, vel_press->nefields,
                                   vel_press->e_to_d_fn,
                                   vel_press->e_to_d_data,
                                   vel_press->evec);
    }
  }
}
