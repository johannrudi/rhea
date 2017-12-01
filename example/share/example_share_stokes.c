/**
 * Shared functions for rhea examples.
 *
 ******************************************************************************
 * Author:             Johann Rudi <johann@ices.utexas.edu>
 *****************************************************************************/

#include <example_share_stokes.h>
#include <rhea_base.h>
#include <rhea_velocity.h>
#include <rhea_weakzone.h>

void
example_share_stokes_new (rhea_stokes_problem_t **stokes_problem,
                          ymir_mesh_t *ymir_mesh,
                          ymir_pressure_elem_t *press_elem,
                          rhea_temperature_options_t *temp_options,
                          rhea_weakzone_options_t *weak_options,
                          rhea_viscosity_options_t *visc_options,
                          rhea_newton_options_t *newton_options,
                          char *vtk_write_newton_itn_path)
{
  const char         *this_fn_name = "example_share_stokes_new";
  rhea_domain_options_t *domain_options = visc_options->domain_options;
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  ymir_vec_t         *temperature;
  void               *solver_options;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* set up data */
  rhea_weakzone_data_create (weak_options, mpicomm);

  /* compute temperature */
  temperature = rhea_temperature_new (ymir_mesh);
  rhea_temperature_compute (temperature, temp_options);

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
      ymir_mesh, press_elem, temperature, domain_options, temp_options,
      weak_options, visc_options, solver_options);
  if (RHEA_VISCOSITY_NONLINEAR == visc_options->type) {
    rhea_stokes_problem_nonlinear_set_output (vtk_write_newton_itn_path,
                                              *stokes_problem);
  }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

void
example_share_stokes_destroy (rhea_stokes_problem_t *stokes_problem,
                              rhea_temperature_options_t *temp_options,
                              rhea_weakzone_options_t *weak_options,
                              rhea_viscosity_options_t *visc_options)
{
  const char         *this_fn_name = "example_share_stokes_destroy";
  ymir_vec_t         *temperature;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* get temperature */
  temperature = rhea_stokes_problem_get_temperature (stokes_problem);

  /* destroy Stokes problem */
  rhea_stokes_problem_destroy (stokes_problem);

  /* destroy vectors */
  if (temperature != NULL) {
    rhea_temperature_destroy (temperature);
  }

  /* destroy data */
  rhea_weakzone_data_clear (weak_options);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}
