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

/* constants: field names for vtk files */
#define RHEA_VTK_NAME_TEMPERATURE "temperature"
#define RHEA_VTK_NAME_BACKGROUND_TEMPERATURE "background_temp"
#define RHEA_VTK_NAME_WEAKZONE "weakzone"
#define RHEA_VTK_NAME_VISCOSITY "viscosity"
#define RHEA_VTK_NAME_BOUNDS_MARKER "bounds_marker"
#define RHEA_VTK_NAME_YIELDING_MARKER "yielding_marker"
#define RHEA_VTK_NAME_VELOCITY "velocity"
#define RHEA_VTK_NAME_PRESSURE "pressure"
#define RHEA_VTK_NAME_VELOCITY_RHS "rhs_vel"
#define RHEA_VTK_NAME_STRAINRATE_SQRT_2INV "strain_rate_sqrt_2inv"

void
rhea_vtk_write_input_data (const char *filepath,
                           ymir_vec_t *temperature,
                           ymir_vec_t *background_temp,
                           ymir_vec_t *weakzone,
                           ymir_vec_t *viscosity,
                           ymir_vec_t *bounds_marker,
                           ymir_vec_t *rhs_vel)
{
  const char         *this_fn_name = "rhea_vtk_write_input_data";
  const int           in_temp = (temperature != NULL);
  const int           in_back = (background_temp != NULL);
  const int           in_weak = (weakzone != NULL);
  const int           in_bounds = (bounds_marker != NULL);
  ymir_mesh_t        *ymir_mesh;

  RHEA_GLOBAL_INFOF ("Into %s: Write VTK file \"%s\"\n",
                     this_fn_name, filepath);

  /* check input */
  RHEA_ASSERT (viscosity != NULL);
  RHEA_ASSERT (rhs_vel != NULL);

  /* get ymir mesh */
  ymir_mesh = ymir_vec_get_mesh (viscosity);

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
  if (!in_bounds) {
    bounds_marker = rhea_viscosity_new (ymir_mesh);
    ymir_vec_set_value (bounds_marker, 0.0);
  }

  /* write vtk file (reduce output for common use cases) */
  if (!in_temp && !in_back && !in_weak && !in_bounds) {
    ymir_vtk_write (ymir_mesh, filepath,
                    viscosity, RHEA_VTK_NAME_VISCOSITY,
                    rhs_vel, RHEA_VTK_NAME_VELOCITY_RHS, NULL);
  }
  else if (!in_temp && !in_back) {
    ymir_vtk_write (ymir_mesh, filepath,
                    viscosity, RHEA_VTK_NAME_VISCOSITY,
                    rhs_vel, RHEA_VTK_NAME_VELOCITY_RHS,
                    weakzone, RHEA_VTK_NAME_WEAKZONE,
                    bounds_marker, RHEA_VTK_NAME_BOUNDS_MARKER, NULL);
  }
  else if (!in_weak && !in_bounds) {
    ymir_vtk_write (ymir_mesh, filepath,
                    viscosity, RHEA_VTK_NAME_VISCOSITY,
                    rhs_vel, RHEA_VTK_NAME_VELOCITY_RHS,
                    temperature, RHEA_VTK_NAME_TEMPERATURE,
                    background_temp, RHEA_VTK_NAME_BACKGROUND_TEMPERATURE,
                    NULL);
  }
  else if (!in_bounds) {
    ymir_vtk_write (ymir_mesh, filepath,
                    viscosity, RHEA_VTK_NAME_VISCOSITY,
                    rhs_vel, RHEA_VTK_NAME_VELOCITY_RHS,
                    temperature, RHEA_VTK_NAME_TEMPERATURE,
                    background_temp, RHEA_VTK_NAME_BACKGROUND_TEMPERATURE,
                    weakzone, RHEA_VTK_NAME_WEAKZONE, NULL);
  }
  else { /* otherwise write all fields */
    ymir_vtk_write (ymir_mesh, filepath,
                    viscosity, RHEA_VTK_NAME_VISCOSITY,
                    rhs_vel, RHEA_VTK_NAME_VELOCITY_RHS,
                    temperature, RHEA_VTK_NAME_TEMPERATURE,
                    background_temp, RHEA_VTK_NAME_BACKGROUND_TEMPERATURE,
                    weakzone, RHEA_VTK_NAME_WEAKZONE,
                    bounds_marker, RHEA_VTK_NAME_BOUNDS_MARKER, NULL);
  }

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
  if (!in_bounds) {
    rhea_viscosity_destroy (bounds_marker);
  }

  RHEA_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

static void
rhea_vtk_write_primary (const char *filepath,
                        ymir_vec_t *velocity,
                        ymir_vec_t *pressure)
{
  const char         *this_fn_name = "rhea_vtk_write_primary";
  const int           in_vel = (velocity != NULL);
  const int           in_press = (pressure != NULL);
  ymir_mesh_t        *ymir_mesh;

  RHEA_GLOBAL_INFOF ("Into %s (file path \"%s\")\n",
                     this_fn_name, filepath);

  /* check input */
  RHEA_ASSERT (in_vel || in_press);

  /* write vtk file */
  if (in_vel && in_press) {
    ymir_mesh = ymir_vec_get_mesh (velocity);
    ymir_vtk_write (ymir_mesh, filepath,
                    velocity, RHEA_VTK_NAME_VELOCITY,
                    pressure, RHEA_VTK_NAME_PRESSURE, NULL);
  }
  else if (in_vel) {
    ymir_mesh = ymir_vec_get_mesh (velocity);
    ymir_vtk_write (ymir_mesh, filepath,
                    velocity, RHEA_VTK_NAME_VELOCITY, NULL);
  }
  else if (in_press) {
    ymir_mesh = ymir_vec_get_mesh (pressure);
    ymir_vtk_write (ymir_mesh, filepath,
                    pressure, RHEA_VTK_NAME_PRESSURE, NULL);
  }
  else {
    RHEA_ABORT_NOT_REACHED ();
  }

  RHEA_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

static void
rhea_vtk_write_secondary (const char *filepath,
                          ymir_vec_t *viscosity,
                          ymir_vec_t *velocity)
{
  const char         *this_fn_name = "rhea_vtk_write_secondary";
  const int           in_vel = (velocity != NULL);
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
                    strainrate_sqrt_2inv, RHEA_VTK_NAME_STRAINRATE_SQRT_2INV,
                    viscosity, RHEA_VTK_NAME_VISCOSITY, NULL);
  }
  else {
    ymir_vtk_write (ymir_mesh, filepath,
                    viscosity, RHEA_VTK_NAME_VISCOSITY, NULL);
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
                                           ymir_vec_t *velocity,
                                           ymir_vec_t *pressure,
                                           ymir_vec_t *viscosity,
                                           ymir_vec_t *bounds_marker,
                                           ymir_vec_t *yielding_marker)
{
  const char         *this_fn_name =
                        "rhea_vtk_write_nonlinear_stokes_iteration";
  ymir_mesh_t        *ymir_mesh;
  ymir_vec_t         *strainrate_sqrt_2inv;

  RHEA_GLOBAL_INFOF ("Into %s (file path \"%s\")\n",
                     this_fn_name, filepath);

  /* check input */
  RHEA_ASSERT (velocity != NULL);
  RHEA_ASSERT (pressure != NULL);
  RHEA_ASSERT (viscosity != NULL);
  RHEA_ASSERT (bounds_marker != NULL);
  RHEA_ASSERT (yielding_marker != NULL);

  /* create work variables */
  ymir_mesh = ymir_vec_get_mesh (viscosity);
  strainrate_sqrt_2inv = rhea_strainrate_2inv_new (ymir_mesh);

  /* compute sqrt of the 2nd inv. of the strain rate */
  rhea_strainrate_compute_sqrt_of_2inv (strainrate_sqrt_2inv, velocity);

  /* write vtk file */
  ymir_vtk_write (ymir_mesh, filepath,
                  velocity, RHEA_VTK_NAME_VELOCITY,
                  pressure, RHEA_VTK_NAME_PRESSURE,
                  strainrate_sqrt_2inv, RHEA_VTK_NAME_STRAINRATE_SQRT_2INV,
                  viscosity, RHEA_VTK_NAME_VISCOSITY,
                  bounds_marker, RHEA_VTK_NAME_BOUNDS_MARKER,
                  yielding_marker, RHEA_VTK_NAME_YIELDING_MARKER, NULL);

  /* destroy */
  rhea_strainrate_2inv_destroy (strainrate_sqrt_2inv);

  RHEA_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}
