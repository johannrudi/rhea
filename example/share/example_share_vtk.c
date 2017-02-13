/**
 * Shared functions for rhea examples.
 *
 ******************************************************************************
 * Author:             Johann Rudi <johann@ices.utexas.edu>
 *****************************************************************************/

#include <example_share_vtk.h>
#include <rhea_base.h>
#include <rhea_velocity.h>
#include <rhea_pressure.h>
#include <rhea_velocity_pressure.h>
#include <rhea_vtk.h>

void
example_share_vtk_write_input_data (const char *vtk_write_input_path,
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_temperature_options_t *temp_options,
                                    rhea_viscosity_options_t *visc_options)
{
  ymir_mesh_t        *ymir_mesh;
  ymir_vec_t         *temperature, *weakzone;
  ymir_vec_t         *rhs_vel, *rhs_vel_nonzero_dirichlet;
  ymir_vec_t         *background_temp, *viscosity, *bounds_marker;

  /* exit if nothing to do */
  if (vtk_write_input_path == NULL) {
    return;
  }

  /* get mesh */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);

  /* get vectors */
  temperature = rhea_stokes_problem_get_temperature (stokes_problem);
  weakzone = rhea_stokes_problem_get_weakzone (stokes_problem);
  rhs_vel = rhea_stokes_problem_get_rhs_vel (stokes_problem);
  rhs_vel_nonzero_dirichlet =
    rhea_stokes_problem_get_rhs_vel_nonzero_dirichlet (stokes_problem);

  /* compute background temperature */
  background_temp = rhea_temperature_new (ymir_mesh);
  rhea_temperature_background_compute (background_temp, temp_options);

  /* compute viscosity */
  viscosity = rhea_viscosity_new (ymir_mesh);
  bounds_marker = rhea_viscosity_new (ymir_mesh);
  switch (visc_options->type) {
  case RHEA_VISCOSITY_LINEAR:
    rhea_viscosity_compute (viscosity,
                            NULL /* nl. Stokes output */,
                            bounds_marker,
                            NULL /* nl. Stokes output */,
                            temperature, weakzone,
                            NULL /* nl. Stokes input */,
                            visc_options);
    break;
  case RHEA_VISCOSITY_NONLINEAR:
    rhea_viscosity_compute_nonlinear_init (viscosity,
                                           NULL /* nl. Stokes output */,
                                           bounds_marker,
                                           NULL /* nl. Stokes output */,
                                           temperature, weakzone,
                                           visc_options);
    break;
  default: /* unknown viscosity type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* write vtk */
  rhea_vtk_write_input_data (vtk_write_input_path, temperature,
                             background_temp, weakzone, viscosity,
                             bounds_marker, rhs_vel);

  /* destroy */
  rhea_temperature_destroy (background_temp);
  rhea_viscosity_destroy (viscosity);
  rhea_viscosity_destroy (bounds_marker);
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
  if (vtk_write_solution_path == NULL) {
    return;
  }

  /* get mesh */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);

  /* create work variables */
  velocity = rhea_velocity_new (ymir_mesh);
  pressure = rhea_pressure_new (ymir_mesh, press_elem);
  viscosity = rhea_viscosity_new (ymir_mesh);

  /* get fields */
  rhea_velocity_pressure_copy_components (velocity, pressure, sol_vel_press,
                                          press_elem);
  rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

  /* write vtk */
  rhea_vtk_write_solution (vtk_write_solution_path, velocity, pressure,
                           viscosity);

  /* destroy */
  rhea_velocity_destroy (velocity);
  rhea_pressure_destroy (pressure);
  rhea_viscosity_destroy (viscosity);
}
