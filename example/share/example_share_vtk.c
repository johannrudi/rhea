#include <example_share_vtk.h>
#include <rhea_base.h>
#include <rhea_temperature.h>
#include <rhea_composition.h>
#include <rhea_viscosity.h>
#include <rhea_velocity.h>
#include <rhea_pressure.h>
#include <rhea_velocity_pressure.h>
#include <rhea_strainrate.h>
#include <rhea_stress.h>
#include <rhea_vtk.h>
#include <ymir_interp_vec.h>

void
example_share_vtk_write_input_data (const char *vtk_write_input_path,
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_plate_options_t *plate_options)
{
  rhea_domain_options_t      *domain_options;
  rhea_temperature_options_t *temp_options;
  rhea_composition_options_t *comp_options;
  rhea_viscosity_options_t   *visc_options;
  ymir_mesh_t        *ymir_mesh;
  ymir_vec_t         *temperature, *background_temp;
  ymir_vec_t         *composition;
  ymir_vec_t         *weakzone, *viscosity, *marker;
  ymir_vec_t         *rhs_vel;
//ymir_vec_t         *rhs_vel_nonzero_dirichlet;
  const int           plates_exist =
                        (0 < rhea_plate_get_n_plates (plate_options));
  ymir_vec_t         *plate_label, *plate_weight, *plate_vel;
  char                path[BUFSIZ];

  /* exit if nothing to do */
  if (vtk_write_input_path == NULL) {
    return;
  }

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (__func__, "path=\"%s\"",
                                 vtk_write_input_path);

  /* get options */
  domain_options = rhea_stokes_problem_get_domain_options (stokes_problem);
  temp_options = rhea_stokes_problem_get_temperature_options (stokes_problem);
  comp_options = rhea_stokes_problem_get_composition_options (stokes_problem);
  visc_options = rhea_stokes_problem_get_viscosity_options (stokes_problem);

  /* get mesh */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);

  /* get vectors */
  temperature = ymir_vec_clone (
      rhea_stokes_problem_get_temperature (stokes_problem));
  rhs_vel = ymir_vec_clone (
      rhea_stokes_problem_get_rhs_vel (stokes_problem));
//rhs_vel_nonzero_dirichlet =
//  rhea_stokes_problem_get_rhs_vel_nonzero_dirichlet (stokes_problem);
  if (rhea_composition_exists (comp_options)) {
    composition = ymir_vec_clone (
        rhea_stokes_problem_get_composition_density (stokes_problem));
  }
  else {
    composition = NULL;
  }

  /* compute background temperature */
  background_temp = rhea_temperature_new (ymir_mesh);
  rhea_temperature_background_compute (background_temp, temp_options);


  /* get plate labels */
  if (plates_exist) {
    plate_label  = rhea_viscosity_surface_new (ymir_mesh);
    plate_weight = rhea_viscosity_surface_new (ymir_mesh);
    plate_vel    = rhea_velocity_surface_new (ymir_mesh);
    rhea_plate_set_label_vec (plate_label, plate_options);
    rhea_plate_set_weight_vec (plate_weight, NULL, NULL, plate_options, NULL);
    rhea_plate_velocity_generate_all (plate_vel, plate_options);
  }
  else {
    plate_label = NULL;
    plate_vel = NULL;
  }

  /* compute weak zone */
  weakzone = rhea_weakzone_new (ymir_mesh);
  rhea_stokes_problem_weakzone_compute (weakzone, stokes_problem);

  /* compute viscosity */
  viscosity = rhea_viscosity_new (ymir_mesh);
  marker = rhea_viscosity_new (ymir_mesh);
  switch (visc_options->type) {
  case RHEA_VISCOSITY_LINEAR:
    rhea_stokes_problem_viscosity_compute (
        /* out: */ viscosity, NULL, marker,
        /* in:  */ temperature, weakzone, NULL, stokes_problem);
    break;
  case RHEA_VISCOSITY_NONLINEAR:
    rhea_viscosity_compute_nonlinear_init (
        /* out: */ viscosity, NULL, marker,
        /* in:  */ temperature, weakzone, visc_options);
    break;
  default: /* unknown viscosity type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* convert to physical dimensions */
  rhea_temperature_convert_to_dimensional_K (temperature, temp_options);
  rhea_temperature_convert_to_dimensional_K (background_temp, temp_options);
  rhea_viscosity_convert_to_dimensional_Pas (viscosity, visc_options);
  rhea_velocity_convert_to_dimensional_mm_yr (rhs_vel, domain_options,
                                              temp_options);
  if (plate_vel != NULL) {
    rhea_velocity_convert_to_dimensional_mm_yr (plate_vel, domain_options,
                                                temp_options);
  }

  /* write vtk */
  rhea_vtk_write_input_data (vtk_write_input_path, temperature,
                             background_temp, composition,
                             weakzone, viscosity, marker, rhs_vel);
  if (plates_exist) {
    snprintf (path, BUFSIZ, "%s_obs", vtk_write_input_path);
    rhea_vtk_write_observation_data (path, plate_label, plate_weight,
                                     plate_vel);
  }

  /* destroy */
  ymir_vec_destroy (temperature);
  rhea_temperature_destroy (background_temp);
  rhea_weakzone_destroy (weakzone);
  rhea_viscosity_destroy (viscosity);
  rhea_viscosity_destroy (marker);
  ymir_vec_destroy (rhs_vel);
  if (rhea_composition_exists (comp_options)) {
    ymir_vec_destroy (composition);
  }
  if (plates_exist) {
    rhea_viscosity_surface_destroy (plate_label);
    rhea_viscosity_surface_destroy (plate_weight);
    rhea_velocity_surface_destroy (plate_vel);
  }

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

void
example_share_vtk_write_solution (const char *vtk_write_solution_path,
                                  ymir_vec_t *sol_vel_press,
                                  rhea_stokes_problem_t *stokes_problem)
{
  rhea_domain_options_t      *domain_options;
  rhea_temperature_options_t *temp_options;
  rhea_viscosity_options_t   *visc_options;
  rhea_weakzone_options_t    *weak_options;
  ymir_mesh_t          *ymir_mesh;
  ymir_pressure_elem_t *press_elem;
  ymir_vec_t         *velocity, *pressure, *viscosity, *marker, *stress;
  ymir_vec_t         *weak_normal;
  ymir_vec_t         *velocity_surf, *stress_norm_surf, *viscosity_surf,
                     *marker_surf;
  double              strainrate_dim_1_s;

  /* exit if nothing to do */
  if (vtk_write_solution_path == NULL) {
    return;
  }

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (__func__, "path=\"%s\"",
                                 vtk_write_solution_path);

  /* get options */
  domain_options = rhea_stokes_problem_get_domain_options (stokes_problem);
  temp_options = rhea_stokes_problem_get_temperature_options (stokes_problem);
  visc_options = rhea_stokes_problem_get_viscosity_options (stokes_problem);
  weak_options = rhea_stokes_problem_get_weakzone_options (stokes_problem);

  /* get mesh */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);

  /* create volume variables */
  velocity = rhea_velocity_new (ymir_mesh);
  pressure = rhea_pressure_new (ymir_mesh, press_elem);
  viscosity = rhea_viscosity_new (ymir_mesh);
  marker = rhea_viscosity_new (ymir_mesh);
  stress = rhea_stress_new (ymir_mesh);

  /* get volume fields */
  rhea_velocity_pressure_copy_components (velocity, pressure, sol_vel_press,
                                          press_elem);
  rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);
  rhea_stokes_problem_copy_marker (marker, stokes_problem);
  rhea_stokes_problem_stress_compute (
      stress, sol_vel_press, stokes_problem,
      NULL /* !override_stress_op */, 0 /* !linearized */,
      0 /* !skip_pressure */);

  /* create normals to weak zones */
  if (rhea_weakzone_exists (weak_options)) {
    weak_normal = rhea_weakzone_normal_new (ymir_mesh);
    rhea_weakzone_compute_normal (weak_normal, weak_options);
  }
  else {
    weak_normal = NULL;
  }

  /* create surface variables */
  velocity_surf = rhea_velocity_surface_new (ymir_mesh);
  stress_norm_surf = rhea_stress_surface_new (ymir_mesh);
  viscosity_surf = rhea_viscosity_surface_new (ymir_mesh);
  marker_surf = rhea_viscosity_surface_new (ymir_mesh);

  /* get surface fields */
  rhea_velocity_surface_interpolate (velocity_surf, velocity);
  rhea_stokes_problem_stress_compute_normal_at_surface (stress_norm_surf,
                                                        sol_vel_press,
                                                        stokes_problem);
  rhea_viscosity_surface_interpolate (viscosity_surf, viscosity,
                                      visc_options->min, visc_options->max);
  ymir_interp_vec (marker, marker_surf);

  /* convert to physical dimensions */
  rhea_velocity_convert_to_dimensional_mm_yr (velocity, domain_options,
                                              temp_options);
  rhea_velocity_convert_to_dimensional_mm_yr (velocity_surf, domain_options,
                                              temp_options);
  rhea_pressure_convert_to_dimensional_Pa (pressure, domain_options,
                                           temp_options, visc_options);
  rhea_stress_convert_to_dimensional_Pa (stress_norm_surf, domain_options,
                                         temp_options, visc_options);
  rhea_viscosity_convert_to_dimensional_Pas (viscosity, visc_options);
  rhea_viscosity_convert_to_dimensional_Pas (viscosity_surf, visc_options);
  rhea_stress_convert_to_dimensional_Pa (stress, domain_options, temp_options,
                                         visc_options);
  strainrate_dim_1_s =
    rhea_strainrate_get_dim_1_s (domain_options, temp_options) /
    rhea_velocity_get_dim_mm_yr (domain_options, temp_options);

  /* write vtk */
  rhea_vtk_write_solution (vtk_write_solution_path, velocity, pressure,
                           viscosity, marker, stress, weak_normal,
                           strainrate_dim_1_s);
  rhea_vtk_write_solution_surf (vtk_write_solution_path, velocity_surf,
                                stress_norm_surf, viscosity_surf, marker_surf);

  /* destroy */
  rhea_velocity_destroy (velocity);
  rhea_pressure_destroy (pressure);
  rhea_viscosity_destroy (viscosity);
  rhea_viscosity_destroy (marker);
  rhea_stress_destroy (stress);
  rhea_velocity_surface_destroy (velocity_surf);
  rhea_stress_surface_destroy (stress_norm_surf);
  rhea_viscosity_surface_destroy (viscosity_surf);
  rhea_viscosity_surface_destroy (marker_surf);
  if (rhea_weakzone_exists (weak_options)) {
    rhea_weakzone_normal_destroy (weak_normal);
  }

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}
