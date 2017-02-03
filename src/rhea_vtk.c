/*
 */

#include <rhea_vtk.h>
#include <rhea_base.h>
#include <rhea_temperature.h>
#include <rhea_weakzone.h>
#include <rhea_viscosity.h>
#include <rhea_velocity.h>
#include <rhea_pressure.h>
#include <rhea_strainrate.h>
#include <ymir_vtk.h>

void
rhea_vtk_write_input_data (const char *filepath,
                           ymir_vec_t *temperature,
                           ymir_vec_t *background_temp,
                           ymir_vec_t *weakzone,
                           ymir_vec_t *viscosity,
                           ymir_vec_t *rhs_vel)
{
  const char         *this_fn_name = "rhea_vtk_write_input_data";
  const int           in_temp = (temperature != NULL ? 1 : 0);
  const int           in_back = (background_temp != NULL ? 1 : 0);
  const int           in_weak = (weakzone != NULL ? 1 : 0);
  const int           in_visc = (viscosity != NULL ? 1 : 0);
  ymir_mesh_t        *ymir_mesh;

  RHEA_GLOBAL_INFOF ("Into %s: Write VTK file \"%s\"\n",
                     this_fn_name, filepath);

  /* check input */
  RHEA_ASSERT (in_temp || in_back || in_weak || in_visc || rhs_vel != NULL);

  /* get ymir mesh */
  if (in_temp) {
    ymir_mesh = ymir_vec_get_mesh (temperature);
  }
  else if (in_back) {
    ymir_mesh = ymir_vec_get_mesh (background_temp);
  }
  else if (in_weak) {
    ymir_mesh = ymir_vec_get_mesh (weakzone);
  }
  else if (in_visc) {
    ymir_mesh = ymir_vec_get_mesh (viscosity);
  }
  else {
    ymir_mesh = ymir_vec_get_mesh (rhs_vel);
  }

  /* create vectors that were not given as input */
  if (!in_temp) {
    temperature = rhea_temperature_new (ymir_mesh);
    ymir_vec_set_value (temperature, -1.0);
  }
  if (!in_back) {
    background_temp = rhea_temperature_new (ymir_mesh);
    ymir_vec_set_value (background_temp, -1.0);
  }
  if (!in_weak) {
    weakzone = rhea_weakzone_new (ymir_mesh);
    ymir_vec_set_value (weakzone, -1.0);
  }
  if (!in_visc) {
    viscosity = rhea_viscosity_new (ymir_mesh);
    ymir_vec_set_value (viscosity, -1.0);
  }

  /* write vtk file */
  ymir_vtk_write (ymir_mesh, filepath,
                  temperature, "temperature",
                  background_temp, "background_temp",
                  weakzone, "weakzone",
                  viscosity, "viscosity",
                  rhs_vel, "rhs_vel",
                  NULL);

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

static void
rhea_vtk_write_primary (const char *filepath,
                        ymir_vec_t *velocity,
                        ymir_vec_t *pressure)
{
  const char         *this_fn_name = "rhea_vtk_write_primary";
  const int           in_vel = (velocity != NULL ? 1 : 0);
  const int           in_press = (pressure != NULL ? 1 : 0);
  ymir_mesh_t        *ymir_mesh;

  RHEA_GLOBAL_INFOF ("Into %s (file path \"%s\")\n",
                     this_fn_name, filepath);

  /* check input */
  RHEA_ASSERT (in_vel || in_press);

  /* get ymir mesh */
  if (in_vel) {
    ymir_mesh = ymir_vec_get_mesh (velocity);
  }
  else if (in_press) {
    ymir_mesh = ymir_vec_get_mesh (pressure);
  }
  else {
    RHEA_ABORT_NOT_REACHED ();
  }

  /* create vectors that were not given as input */
  if (!in_vel) {
    velocity = rhea_velocity_new (ymir_mesh);
    ymir_vec_set_value (velocity, 0.0);
  }
  if (!in_press) {
    pressure = ymir_dvec_new (ymir_mesh, 1, YMIR_GLL_NODE); /* since press_elem
                                                             * not available */
    ymir_vec_set_value (pressure, 0.0);
  }

  /* write vtk file */
  ymir_vtk_write (ymir_mesh, filepath,
                  velocity, "velocity",
                  pressure, "pressure", NULL);

  /* destroy */
  if (!in_vel) {
    rhea_velocity_destroy (velocity);
  }
  if (!in_press) {
    ymir_vec_destroy (pressure);
  }

  RHEA_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

static void
rhea_vtk_write_secondary (const char *filepath,
                          ymir_vec_t *viscosity,
                          ymir_vec_t *velocity)
{
  const char         *this_fn_name = "rhea_vtk_write_secondary";
  const int           in_vel = (velocity != NULL ? 1 : 0);
  ymir_mesh_t        *ymir_mesh;
  ymir_vec_t         *strainrate_sqrt_2inv;

  RHEA_GLOBAL_INFOF ("Into %s (file path \"%s\")\n",
                     this_fn_name, filepath);

  /* check input */
  RHEA_ASSERT (viscosity != NULL);

  /* get ymir mesh */
  ymir_mesh = ymir_vec_get_mesh (viscosity);

  /* compute the square root of the 2nd inv. of the strain rate */
  if (in_vel) {
    strainrate_sqrt_2inv = rhea_strainrate_2inv_new (ymir_mesh);
    rhea_strainrate_compute_sqrt_of_2inv (strainrate_sqrt_2inv, velocity);
  }

  /* write vtk file */
  if (in_vel) {
    ymir_vtk_write (ymir_mesh, filepath,
                    strainrate_sqrt_2inv, "strain_rate_sqrt_2inv",
                    viscosity, "viscosity", NULL);
  }
  else {
    ymir_vtk_write (ymir_mesh, filepath,
                    viscosity, "viscosity", NULL);
  }

  /* destroy */
  if (in_vel) {
    rhea_strainrate_2inv_destroy (strainrate_sqrt_2inv);
  }

  RHEA_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

void
rhea_vtk_write_solution (const char *filepath,
                         ymir_vec_t *velocity,
                         ymir_vec_t *pressure,
                         ymir_vec_t *viscosity)
{
  const char         *this_fn_name = "rhea_vtk_write_solution";
  char                path[BUFSIZ];

  RHEA_GLOBAL_INFOF ("Into %s (file path \"%s\")\n",
                     this_fn_name, filepath);

  snprintf (path, BUFSIZ, "%s_primary", filepath);
  rhea_vtk_write_primary (path, velocity, pressure);

  snprintf (path, BUFSIZ, "%s_secondary", filepath);
  rhea_vtk_write_secondary (path, viscosity, velocity);

  RHEA_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

void
rhea_vtk_write_nonlinear_stokes_iteration (const char *filepath,
                                           ymir_vec_t *velocity_pressure,
                                           ymir_vec_t *viscosity,
                                           ymir_vec_t *bounds_marker,
                                           ymir_vec_t *yielding_marker,
                                           ymir_pressure_elem_t *press_elem)
{
  const char         *this_fn_name =
                        "rhea_vtk_write_nonlinear_stokes_iteration";
  ymir_mesh_t        *ymir_mesh;
  ymir_vec_t         *velocity, *pressure;
  ymir_vec_t         *strainrate_sqrt_2inv;

  RHEA_GLOBAL_INFOF ("Into %s (file path \"%s\")\n",
                     this_fn_name, filepath);

  /* check input */
  RHEA_ASSERT (velocity_pressure != NULL);
  RHEA_ASSERT (viscosity != NULL);
  RHEA_ASSERT (bounds_marker != NULL);
  RHEA_ASSERT (yielding_marker != NULL);

  /* create work variables */
  ymir_mesh = ymir_vec_get_mesh (velocity_pressure);
  velocity = rhea_velocity_new (ymir_mesh);
  pressure = rhea_pressure_new (ymir_mesh, press_elem);
  strainrate_sqrt_2inv = rhea_strainrate_2inv_new (ymir_mesh);

  /* get velocity and pressure components */
  rhea_velocity_pressure_copy_components (velocity, pressure, velocity_pressure,
                                          press_elem);

  /* compute sqrt of the 2nd inv. of the strain rate */
  rhea_strainrate_compute_sqrt_of_2inv (strainrate_sqrt_2inv, velocity);

  /* write vtk file */
  ymir_vtk_write (ymir_mesh, filepath,
                  velocity, "velocity",
                  pressure, "pressure",
                  strainrate_sqrt_2inv, "strain_rate_sqrt_2inv",
                  viscosity, "viscosity",
                  bounds_marker, "bounds_marker",
                  yielding_marker, "yielding_marker", NULL);

  /* destroy */
  rhea_velocity_destroy (velocity);
  rhea_pressure_destroy (pressure);
  rhea_strainrate_2inv_destroy (strainrate_sqrt_2inv);

  RHEA_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}
