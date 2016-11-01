/*
 */

#include <rhea_vtk.h>
#include <rhea_base.h>
#include <rhea_temperature.h>
#include <rhea_weakzone.h>
#include <rhea_viscosity.h>
#include <ymir_vtk.h>

void
rhea_vtk_write_simple (const char *filepath,
                       ymir_vec_t *temperature,
                       ymir_vec_t *background_temp,
                       ymir_vec_t *weakzone,
                       ymir_vec_t *viscosity,
                       ymir_vec_t *rhs_vel)
{
  const char         *this_fn_name = "rhea_vtk_write_simple";
  const int           in_temp = (temperature != NULL ? 1 : 0);
  const int           in_back = (background_temp != NULL ? 1 : 0);
  const int           in_weak = (weakzone != NULL ? 1 : 0);
  const int           in_visc = (viscosity != NULL ? 1 : 0);
  ymir_mesh_t        *mesh;

  RHEA_GLOBAL_INFOF ("Into %s: Write VTK file \"%s\"\n",
                     this_fn_name, filepath);

  /* check input */
  RHEA_ASSERT (in_temp || in_back || in_weak || in_visc || rhs_vel != NULL);

  /* get ymir mesh */
  if (in_temp) {
    mesh = ymir_vec_get_mesh (temperature);
  }
  else if (in_back) {
    mesh = ymir_vec_get_mesh (background_temp);
  }
  else if (in_weak) {
    mesh = ymir_vec_get_mesh (weakzone);
  }
  else if (in_visc) {
    mesh = ymir_vec_get_mesh (viscosity);
  }
  else {
    mesh = ymir_vec_get_mesh (rhs_vel);
  }

  /* create vectors that were not given as input */
  if (!in_temp) {
    temperature = rhea_temperature_new (mesh);
    ymir_vec_set_value (temperature, -1.0);
  }
  if (!in_back) {
    background_temp = rhea_temperature_new (mesh);
    ymir_vec_set_value (background_temp, -1.0);
  }
  if (!in_weak) {
    weakzone = rhea_weakzone_new (mesh);
    ymir_vec_set_value (weakzone, -1.0);
  }
  if (!in_visc) {
    viscosity = rhea_viscosity_new (mesh);
    ymir_vec_set_value (viscosity, -1.0);
  }

  /* write vtk file */
  ymir_vtk_write (
      mesh, filepath,
      temperature, "temperature", background_temp, "background_temp",
      weakzone, "weakzone", viscosity, "viscosity",
      rhs_vel, "rhs_vel", NULL);

  /* destroy */
  if (!in_temp) {
    rhea_temperature_destroy (temperature);
  }
  if (!in_back) {
    rhea_temperature_destroy (background_temp);
  }
  if (!in_weak) {
    rhea_weakzone_destroy (weakzone);
  }
  if (!in_visc) {
    rhea_viscosity_destroy (viscosity);
  }

  RHEA_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}
