/**
 * Shared functions for rhea examples.
 *
 ******************************************************************************
 * Author:             Johann Rudi <johann@ices.utexas.edu>
 *****************************************************************************/

#include <example_share_stokes.h>
#include <rhea_base.h>

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
  void               *solver_options;

  RHEA_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  /* compute temperature */
  temperature = rhea_temperature_new (ymir_mesh);
  rhea_temperature_compute (temperature, temp_options);

  /* compute weak zone */
  weakzone = NULL;  /* TODO create weak zones */

  /* set solver-specific variables */
  if (RHEA_VISCOSITY_NONLINEAR == visc_options->type) {
    solver_options = newton_options;
  }
  else {
    solver_options = NULL;
  }

  /* create Stokes problem */
  *stokes_problem = rhea_stokes_problem_new (
      temperature, weakzone, NULL /* nonzero vel. Dirichlet values */,
      ymir_mesh, press_elem,
      domain_options, temp_options, visc_options, solver_options);
  if (RHEA_VISCOSITY_NONLINEAR == visc_options->type) {
    rhea_stokes_problem_nonlinear_set_output (vtk_write_newton_itn_path,
                                              *stokes_problem);
  }

  /* destroy */
  rhea_temperature_destroy (temperature);

  RHEA_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

void
example_share_stokes_destroy (rhea_stokes_problem_t *stokes_problem)
{
  const char         *this_fn_name = "example_share_stokes_destroy";

  RHEA_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  rhea_stokes_problem_destroy (stokes_problem);

  RHEA_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}
