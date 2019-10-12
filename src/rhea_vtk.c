#include <rhea_vtk.h>
#include <rhea_base.h>
#include <rhea_temperature.h>
#include <rhea_plate.h>
#include <rhea_weakzone.h>
#include <rhea_viscosity.h>
#include <rhea_velocity.h>
#include <rhea_pressure.h>
#include <rhea_strainrate.h>
#include <ymir_vtk.h>

/* constants: field names for vtk files */
#define RHEA_VTK_NAME_TEMPERATURE "temperature"
#define RHEA_VTK_NAME_VELOCITY "velocity"
#define RHEA_VTK_NAME_PRESSURE "pressure"

#define RHEA_VTK_NAME_BACKGROUND_TEMPERATURE "background_temp"
#define RHEA_VTK_NAME_VELOCITY_RHS "rhs_vel"
#define RHEA_VTK_NAME_STRAINRATE_SQRT_2INV "strainrate_sqrt_2inv"
#define RHEA_VTK_NAME_VISCSTRESS_SQRT_2INV "viscstress_sqrt_2inv"
#define RHEA_VTK_NAME_STRESS_NORM "stress_norm"

#define RHEA_VTK_NAME_WEAKZONE "weakzone"
#define RHEA_VTK_NAME_VISCOSITY "viscosity"
#define RHEA_VTK_NAME_MARKER "marker"

#define RHEA_VTK_NAME_PLATE_LABEL "plate_label"
#define RHEA_VTK_NAME_PLATE_VEL "plate_velocity"

#define RHEA_VTK_NAME_VELOCITY_FWD "velocity_fwd"
#define RHEA_VTK_NAME_PRESSURE_FWD "pressure_fwd"
#define RHEA_VTK_NAME_VELOCITY_ADJ "velocity_adj"
#define RHEA_VTK_NAME_PRESSURE_ADJ "pressure_adj"
#define RHEA_VTK_NAME_VELOCITY_OBS "velocity_obs"
#define RHEA_VTK_NAME_VELOCITY_MISFIT "velocity_misfit"

int
rhea_vtk_write_input_data (const char *filepath,
                           ymir_vec_t *temperature,
                           ymir_vec_t *background_temp,
                           ymir_vec_t *weakzone,
                           ymir_vec_t *viscosity,
                           ymir_vec_t *marker,
                           ymir_vec_t *rhs_vel)
{
  const int           in_temp = (temperature != NULL);
  const int           in_back = (background_temp != NULL);
  const int           in_weak = (weakzone != NULL);
  const int           in_marker = (marker != NULL);
  ymir_mesh_t        *ymir_mesh;

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "path=\"%s\"", filepath);

  /* check input */
  RHEA_ASSERT (viscosity != NULL);
  RHEA_ASSERT (rhs_vel != NULL);

  /* get ymir mesh */
  if (viscosity != NULL) {
    ymir_mesh = ymir_vec_get_mesh (viscosity);
  }
  else { /* otherwise return failure */
    return 0;
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
  if (!in_marker) {
    marker = rhea_viscosity_new (ymir_mesh);
    ymir_vec_set_value (marker, 0.0);
  }

  /* write vtk file fields (reduce output for common use cases) */
  if (!in_temp && !in_back && !in_weak && !in_marker) {
    ymir_vtk_write (ymir_mesh, filepath,
                    viscosity, RHEA_VTK_NAME_VISCOSITY,
                    rhs_vel, RHEA_VTK_NAME_VELOCITY_RHS, NULL);
  }
  else if (!in_temp && !in_back) {
    ymir_vtk_write (ymir_mesh, filepath,
                    viscosity, RHEA_VTK_NAME_VISCOSITY,
                    rhs_vel, RHEA_VTK_NAME_VELOCITY_RHS,
                    weakzone, RHEA_VTK_NAME_WEAKZONE,
                    marker, RHEA_VTK_NAME_MARKER, NULL);
  }
  else if (!in_weak && !in_marker) {
    ymir_vtk_write (ymir_mesh, filepath,
                    viscosity, RHEA_VTK_NAME_VISCOSITY,
                    rhs_vel, RHEA_VTK_NAME_VELOCITY_RHS,
                    temperature, RHEA_VTK_NAME_TEMPERATURE,
                    background_temp, RHEA_VTK_NAME_BACKGROUND_TEMPERATURE,
                    NULL);
  }
  else if (!in_marker) {
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
                    marker, RHEA_VTK_NAME_MARKER, NULL);
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
  if (!in_marker) {
    rhea_viscosity_destroy (marker);
  }

  RHEA_GLOBAL_INFO_FN_END (__func__);

  /* return success */
  return 1;
}

int
rhea_vtk_write_observation_data (const char *filepath,
                                 ymir_vec_t *plate_label,
                                 ymir_vec_t *plate_vel)
{
  const int           in_pl_label = (plate_label != NULL);
  const int           in_pl_vel = (plate_vel != NULL);
  ymir_mesh_t        *ymir_mesh;

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "path=\"%s\"", filepath);

  /* check input */
  RHEA_ASSERT (in_pl_label || in_pl_vel);

  /* get ymir mesh */
  if (in_pl_label) {
    ymir_mesh = ymir_vec_get_mesh (plate_label);
  }
  else if (in_pl_vel) {
    ymir_mesh = ymir_vec_get_mesh (plate_vel);
  }
  else { /* otherwise return failure */
    return 0;
  }

  /* create vectors that were not given as input */
  if (!in_pl_label) {
    plate_label = rhea_viscosity_surface_new (ymir_mesh);
    ymir_vec_set_value (plate_label, RHEA_PLATE_NONE);
  }
  if (!in_pl_vel) {
    plate_vel = rhea_velocity_surface_new (ymir_mesh);
    ymir_vec_set_value (plate_vel, NAN);
  }

  /* write vtk file fields (reduce output for common use cases) */
  ymir_vtk_write (ymir_mesh, filepath,
                  plate_label, RHEA_VTK_NAME_PLATE_LABEL,
                  plate_vel, RHEA_VTK_NAME_PLATE_VEL, NULL);

  /* destroy */
  if (!in_pl_label) {
    rhea_viscosity_surface_destroy (plate_label);
  }
  if (!in_pl_vel) {
    rhea_velocity_surface_destroy (plate_vel);
  }

  RHEA_GLOBAL_INFO_FN_END (__func__);

  /* return success */
  return 1;
}

static int
rhea_vtk_write_primary (const char *filepath,
                        ymir_vec_t *velocity,
                        ymir_vec_t *pressure)
{
  const int           in_vel = (velocity != NULL);
  const int           in_press = (pressure != NULL);
  ymir_mesh_t        *ymir_mesh;

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "path=\"%s\"", filepath);

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
  else { /* otherwise return failure */
    return 0;
  }

  RHEA_GLOBAL_INFO_FN_END (__func__);

  /* return success */
  return 1;
}

static int
rhea_vtk_write_secondary (const char *filepath,
                          ymir_vec_t *viscosity,
                          ymir_vec_t *marker,
                          ymir_vec_t *velocity)
{
  const int           in_marker = (marker != NULL);
  const int           in_vel = (velocity != NULL);
  ymir_mesh_t        *ymir_mesh;
  ymir_vec_t         *strainrate_sqrt_2inv;

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "path=\"%s\"", filepath);

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
  if (in_marker && in_vel) {
    ymir_vtk_write (ymir_mesh, filepath,
                    strainrate_sqrt_2inv, RHEA_VTK_NAME_STRAINRATE_SQRT_2INV,
                    viscosity, RHEA_VTK_NAME_VISCOSITY,
                    marker, RHEA_VTK_NAME_MARKER, NULL);
  }
  else if (in_marker) {
    ymir_vtk_write (ymir_mesh, filepath,
                    viscosity, RHEA_VTK_NAME_VISCOSITY,
                    marker, RHEA_VTK_NAME_MARKER, NULL);
  }
  else {
    ymir_vtk_write (ymir_mesh, filepath,
                    viscosity, RHEA_VTK_NAME_VISCOSITY, NULL);
  }

  /* destroy */
  if (in_vel) {
    rhea_strainrate_2inv_destroy (strainrate_sqrt_2inv);
  }

  RHEA_GLOBAL_INFO_FN_END (__func__);

  /* return success */
  return 1;
}

int
rhea_vtk_write_solution (const char *filepath,
                         ymir_vec_t *velocity,
                         ymir_vec_t *pressure,
                         ymir_vec_t *viscosity,
                         ymir_vec_t *marker)
{
  char                path[BUFSIZ];
  int                 success = 0;

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "path=\"%s\"", filepath);

  snprintf (path, BUFSIZ, "%s_primary", filepath);
  success += rhea_vtk_write_primary (path, velocity, pressure);

  snprintf (path, BUFSIZ, "%s_secondary", filepath);
  success += rhea_vtk_write_secondary (path, viscosity, marker, velocity);

  RHEA_GLOBAL_INFO_FN_END (__func__);

  return success;
}

int
rhea_vtk_write_solution_surf (const char *filepath,
                              ymir_vec_t *velocity_surf,
                              ymir_vec_t *stress_norm_surf,
                              ymir_vec_t *viscosity_surf,
                              ymir_vec_t *marker_surf)
{
  const int           in_marker = (marker_surf != NULL);
  ymir_mesh_t        *ymir_mesh;

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "path=\"%s\"", filepath);

  /* check input */
  RHEA_ASSERT (velocity_surf != NULL);
  RHEA_ASSERT (stress_norm_surf != NULL);
  RHEA_ASSERT (viscosity_surf != NULL);

  /* get ymir mesh */
  ymir_mesh = ymir_vec_get_mesh (viscosity_surf);

  /* write vtk file */
  if (in_marker) {
    ymir_vtk_write (ymir_mesh, filepath,
                    velocity_surf, RHEA_VTK_NAME_VELOCITY,
                    stress_norm_surf, RHEA_VTK_NAME_STRESS_NORM,
                    viscosity_surf, RHEA_VTK_NAME_VISCOSITY,
                    marker_surf, RHEA_VTK_NAME_MARKER, NULL);
  }
  else {
    ymir_vtk_write (ymir_mesh, filepath,
                    velocity_surf, RHEA_VTK_NAME_VELOCITY,
                    stress_norm_surf, RHEA_VTK_NAME_STRESS_NORM,
                    viscosity_surf, RHEA_VTK_NAME_VISCOSITY, NULL);
  }

  RHEA_GLOBAL_INFO_FN_END (__func__);

  return 1;
}

int
rhea_vtk_write_nonlinear_stokes_iteration (const char *filepath,
                                           ymir_vec_t *velocity,
                                           ymir_vec_t *pressure,
                                           ymir_vec_t *viscosity,
                                           ymir_vec_t *marker)
{
  ymir_mesh_t        *ymir_mesh;
  ymir_vec_t         *strainrate_sqrt_2inv;

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "path=\"%s\"", filepath);

  /* check input */
  RHEA_ASSERT (velocity != NULL);
  RHEA_ASSERT (pressure != NULL);
  RHEA_ASSERT (viscosity != NULL);
  RHEA_ASSERT (marker != NULL);

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
                  marker, RHEA_VTK_NAME_MARKER, NULL);

  /* destroy */
  rhea_strainrate_2inv_destroy (strainrate_sqrt_2inv);

  RHEA_GLOBAL_INFO_FN_END (__func__);

  return 1;
}

int
rhea_vtk_write_nonlinear_stokes_iteration_surf (const char *filepath,
                                                ymir_vec_t *velocity_surf,
                                                ymir_vec_t *stress_norm_surf,
                                                ymir_vec_t *viscosity_surf)
{
  ymir_mesh_t        *ymir_mesh;

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "path=\"%s\"", filepath);

  /* check input */
  RHEA_ASSERT (velocity_surf != NULL);
  RHEA_ASSERT (stress_norm_surf != NULL);
  RHEA_ASSERT (viscosity_surf != NULL);

  /* create work variables */
  ymir_mesh = ymir_vec_get_mesh (viscosity_surf);

  /* write vtk file */
  ymir_vtk_write (ymir_mesh, filepath,
                  velocity_surf, RHEA_VTK_NAME_VELOCITY,
                  stress_norm_surf, RHEA_VTK_NAME_STRESS_NORM,
                  viscosity_surf, RHEA_VTK_NAME_VISCOSITY, NULL);

  RHEA_GLOBAL_INFO_FN_END (__func__);

  return 1;
}

int
rhea_vtk_write_inversion_iteration (const char *filepath,
                                    ymir_vec_t *velocity_fwd,
                                    ymir_vec_t *pressure_fwd,
                                    ymir_vec_t *velocity_adj,
                                    ymir_vec_t *pressure_adj)
{
  ymir_mesh_t        *ymir_mesh;

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "path=\"%s\"", filepath);

  /* check input */
  RHEA_ASSERT (NULL != velocity_fwd);
  RHEA_ASSERT (NULL != pressure_fwd);
  RHEA_ASSERT (NULL != velocity_adj);
  RHEA_ASSERT (NULL != pressure_adj);

  /* create work variables */
  ymir_mesh = ymir_vec_get_mesh (velocity_fwd);

  /* write vtk file */
  ymir_vtk_write (ymir_mesh, filepath,
                  velocity_fwd, RHEA_VTK_NAME_VELOCITY_FWD,
                  pressure_fwd, RHEA_VTK_NAME_PRESSURE_FWD,
                  velocity_adj, RHEA_VTK_NAME_VELOCITY_ADJ,
                  pressure_adj, RHEA_VTK_NAME_PRESSURE_ADJ, NULL);

  RHEA_GLOBAL_INFO_FN_END (__func__);

  return 1;
}

int
rhea_vtk_write_inversion_iteration_surf (const char *filepath,
                                         ymir_vec_t *velocity_fwd_surf,
                                         ymir_vec_t *velocity_adj_surf,
                                         ymir_vec_t *velocity_obs_surf,
                                         ymir_vec_t *misfit_surf)
{
  ymir_mesh_t        *ymir_mesh;

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "path=\"%s\"", filepath);

  /* check input */
  RHEA_ASSERT (NULL != velocity_fwd_surf);
  RHEA_ASSERT (NULL != velocity_adj_surf);
  RHEA_ASSERT (NULL != velocity_obs_surf);
  RHEA_ASSERT (NULL != misfit_surf);

  /* create work variables */
  ymir_mesh = ymir_vec_get_mesh (velocity_fwd_surf);

  /* write vtk file */
  ymir_vtk_write (ymir_mesh, filepath,
                  velocity_fwd_surf, RHEA_VTK_NAME_VELOCITY_FWD,
                  velocity_adj_surf, RHEA_VTK_NAME_VELOCITY_ADJ,
                  velocity_obs_surf, RHEA_VTK_NAME_VELOCITY_OBS,
                  misfit_surf, RHEA_VTK_NAME_VELOCITY_MISFIT, NULL);

  RHEA_GLOBAL_INFO_FN_END (__func__);

  return 1;
}
