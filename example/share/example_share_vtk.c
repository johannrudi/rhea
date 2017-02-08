/**
 * Shared functions for rhea examples.
 *
 ******************************************************************************
 * Author:             Johann Rudi <johann@ices.utexas.edu>
 *****************************************************************************/

#include <example_share_vtk.h>
#include <rhea_base.h>
#include <rhea_vtk.h>
#include <rhea_temperature.h>
#include <rhea_weakzone.h>
#include <rhea_viscosity.h>
#include <rhea_velocity.h>
#include <rhea_pressure.h>
#include <rhea_velocity_pressure.h>

void
example_share_vtk_write_input_data (const char *vtk_write_input_path,
                                    ymir_mesh_t *ymir_mesh,
                                    rhea_temperature_options_t *temp_options,
                                    rhea_viscosity_options_t *visc_options)
{
  ymir_vec_t         *temperature, *background_temp;
  ymir_vec_t         *weakzone, *viscosity;
  ymir_vec_t         *rhs_vel;

  /* exit if nothing to do */
  if (vtk_write_input_path == NULL) {
    return;
  }

  /* create work variables */
  temperature = rhea_temperature_new (ymir_mesh);
  background_temp = rhea_temperature_new (ymir_mesh);
  weakzone = NULL; //TODO process weakzone
  viscosity = rhea_viscosity_new (ymir_mesh);
  rhs_vel = rhea_velocity_new (ymir_mesh);

  /* compute fields */
  rhea_temperature_compute (temperature, temp_options);
  rhea_temperature_background_compute (background_temp, temp_options);
  if (visc_options->type == RHEA_VISCOSITY_NONLINEAR) {
    rhea_viscosity_compute_init_nonlinear (viscosity,
                                           NULL /* nl. Stokes output */,
                                           NULL /* nl. Stokes output */,
                                           NULL /* nl. Stokes output */,
                                           temperature, weakzone,
                                           visc_options);
  }
  else {
    rhea_viscosity_compute (viscosity,
                            NULL /* nl. Stokes output */,
                            NULL /* nl. Stokes output */,
                            NULL /* nl. Stokes output */,
                            temperature, weakzone,
                            NULL /* nl. Stokes input */,
                            visc_options);
  }
  rhea_temperature_compute_rhs_vel (rhs_vel, temperature, temp_options);

  /* write vtk */
  rhea_vtk_write_input_data (vtk_write_input_path, temperature,
                             background_temp, weakzone, viscosity, rhs_vel);

  /* destroy */
  rhea_temperature_destroy (temperature);
  rhea_temperature_destroy (background_temp);
  rhea_viscosity_destroy (viscosity);
  rhea_velocity_destroy (rhs_vel);
}

void
example_share_vtk_write_solution (const char *vtk_write_solution_path,
                                  ymir_vec_t *sol_vel_press,
                                  rhea_stokes_problem_t *stokes_problem)
{
  ymir_mesh_t        *ymir_mesh;
  ymir_pressure_elem_t *press_elem;
  ymir_vec_t         *velocity, *pressure, *viscosity;

  /* exit if nothing to do */
  if (vtk_write_solution_path != NULL) {
    return;
  }

  /* create work variables */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);
  velocity = rhea_velocity_new (ymir_mesh);
  pressure = rhea_pressure_new (ymir_mesh, press_elem);
  viscosity = rhea_viscosity_new (ymir_mesh);

  /* get fields */
  rhea_velocity_pressure_copy_components (velocity, pressure, sol_vel_press,
                                          press_elem);
  rhea_stokes_problem_get_viscosity (viscosity, stokes_problem);

  /* write vtk */
  rhea_vtk_write_solution (vtk_write_solution_path, velocity, pressure,
                           viscosity);

  /* destroy */
  rhea_velocity_destroy (velocity);
  rhea_pressure_destroy (pressure);
  rhea_viscosity_destroy (viscosity);
}
