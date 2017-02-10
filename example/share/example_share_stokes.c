/**
 * Shared functions for rhea examples.
 *
 ******************************************************************************
 * Author:             Johann Rudi <johann@ices.utexas.edu>
 *****************************************************************************/

#include <example_share_stokes.h>
#include <rhea_base.h>
#include <rhea_velocity.h>

void
example_share_stokes_new (rhea_stokes_problem_t **stokes_problem,
                          ymir_mesh_t *ymir_mesh,
                          ymir_pressure_elem_t *press_elem,
                          rhea_domain_options_t *domain_options,
                          rhea_temperature_options_t *temp_options,
                          rhea_viscosity_options_t *visc_options,
                          rhea_newton_options_t *newton_options,
                          char *vtk_write_newton_itn_path)
{
  const char         *this_fn_name = "example_share_stokes_new";
  ymir_vec_t         *temperature, *weakzone;
  ymir_vec_t         *rhs_vel, *rhs_vel_nonzero_dirichlet;
  void               *solver_options;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* compute temperature */
  temperature = rhea_temperature_new (ymir_mesh);
  rhea_temperature_compute (temperature, temp_options);

  /* compute weak zone */
  weakzone = NULL;  /* TODO create weak zones */

  /* create velocity right-hand side volume forcing */
  rhs_vel = rhea_velocity_new (ymir_mesh);
  rhea_temperature_compute_rhs_vel (rhs_vel, temperature, temp_options);

  /* create the nonzero velocity values at the Dirichlet boundary */
  rhs_vel_nonzero_dirichlet = NULL; //TODO?

  /* set solver-specific variables */
  switch (visc_options->type) {
  case RHEA_VISCOSITY_LINEAR:
    solver_options = NULL;
    break;
  case RHEA_VISCOSITY_NONLINEAR:
    solver_options = newton_options;
    break;
  default: /* unknown viscosity type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* create Stokes problem */
  *stokes_problem = rhea_stokes_problem_new (
      temperature, weakzone, rhs_vel, rhs_vel_nonzero_dirichlet,
      ymir_mesh, press_elem, domain_options, visc_options, solver_options);
  if (RHEA_VISCOSITY_NONLINEAR == visc_options->type) {
    rhea_stokes_problem_nonlinear_set_output (vtk_write_newton_itn_path,
                                              *stokes_problem);
  }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

void
example_share_stokes_destroy (rhea_stokes_problem_t *stokes_problem)
{
  const char         *this_fn_name = "example_share_stokes_destroy";
  ymir_vec_t         *temperature, *weakzone;
  ymir_vec_t         *rhs_vel, *rhs_vel_nonzero_dirichlet;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* get vectors */
  temperature = rhea_stokes_problem_get_temperature (stokes_problem);
  weakzone = rhea_stokes_problem_get_weakzone (stokes_problem);
  rhs_vel = rhea_stokes_problem_get_rhs_vel (stokes_problem);
  rhs_vel_nonzero_dirichlet =
    rhea_stokes_problem_get_rhs_vel_nonzero_dirichlet (stokes_problem);

  /* destroy Stokes problem */
  rhea_stokes_problem_destroy (stokes_problem);

  /* destroy vectors */
  if (temperature != NULL) {
    rhea_temperature_destroy (temperature);
  }
  if (weakzone != NULL) {
    rhea_weakzone_destroy (weakzone);
  }
  if (rhs_vel != NULL) {
    rhea_velocity_destroy (rhs_vel);
  }
  if (rhs_vel_nonzero_dirichlet != NULL) {
    rhea_velocity_destroy (rhs_vel_nonzero_dirichlet);
  }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}
