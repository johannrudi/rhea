#include <example_share_stokes.h>
#include <rhea.h>
#include <rhea_stokes_problem_amr.h>

void
example_share_stokes_new (rhea_stokes_problem_t **stokes_problem,
                          ymir_mesh_t **ymir_mesh,
                          ymir_pressure_elem_t **press_elem,
                          rhea_temperature_options_t *temp_options,
                          rhea_plate_options_t *plate_options,
                          rhea_weakzone_options_t *weak_options,
                          rhea_viscosity_options_t *visc_options,
						  rhea_composition_options_t *comp_options,
                          p4est_t *p4est,
                          rhea_discretization_options_t *discr_options,
                          const int performance_monitor_index_mesh,
                          const int performance_monitor_index_stokes,
                          char *solver_bin_path,
                          char *solver_vtk_path)
{
  rhea_domain_options_t *domain_options = visc_options->domain_options;
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (*ymir_mesh);
  ymir_vec_t         *temperature;
  ymir_vec_t		 *composition;

  RHEA_GLOBAL_PRODUCTION_FN_BEGIN (__func__);

  /* set up data */
  rhea_plate_data_create (plate_options, mpicomm);
  rhea_weakzone_data_create (weak_options, mpicomm);

  /* compute temperature */
  temperature = rhea_temperature_new (*ymir_mesh);
  rhea_temperature_compute (temperature, temp_options);

  /* read in composition */
  composition = rhea_composition_new (*ymir_mesh);
  rhea_composition_read (composition, comp_options);

  /* create Stokes problem */
  rhea_performance_monitor_start_barrier (performance_monitor_index_stokes);
  *stokes_problem = rhea_stokes_problem_new (
      *ymir_mesh, *press_elem, temperature, composition, domain_options, temp_options,
      weak_options, visc_options, comp_options);
  rhea_stokes_problem_set_plate_options (*stokes_problem, plate_options);
  rhea_stokes_problem_set_solver_amr (*stokes_problem, p4est, discr_options);
  rhea_stokes_problem_set_solver_bin_output (*stokes_problem, solver_bin_path);
  rhea_stokes_problem_set_solver_vtk_output (*stokes_problem, solver_vtk_path);
  rhea_performance_monitor_stop_add (performance_monitor_index_stokes);

  /* perform initial AMR */
  if (p4est != NULL && discr_options != NULL) {
    rhea_performance_monitor_start_barrier (performance_monitor_index_mesh);
    rhea_stokes_problem_init_amr (*stokes_problem, p4est, discr_options);
    rhea_performance_monitor_stop_add (performance_monitor_index_mesh);

    /* retrieve adapted mesh */
    *ymir_mesh = rhea_stokes_problem_get_ymir_mesh (*stokes_problem);
    *press_elem = rhea_stokes_problem_get_press_elem (*stokes_problem);
  }

  /* destroy vector composition */
  rhea_composition_destroy (composition);

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);
}

void
example_share_stokes_destroy (rhea_stokes_problem_t *stokes_problem,
                              rhea_temperature_options_t *temp_options,
                              rhea_plate_options_t *plate_options,
                              rhea_weakzone_options_t *weak_options,
                              rhea_viscosity_options_t *visc_options)
{
  ymir_vec_t         *temperature;
  ymir_vec_t		 *compositional_density, *compositional_viscosity;

  RHEA_GLOBAL_PRODUCTION_FN_BEGIN (__func__);

  /* get temperature */
  temperature = rhea_stokes_problem_get_temperature (stokes_problem);
  /* get compositional density */
  compositional_density = rhea_stokes_problem_get_compositional_density (stokes_problem);
  /* get compositional viscosity */
  compositional_viscosity = rhea_stokes_problem_get_compositional_viscosity (stokes_problem);

  /* destroy Stokes problem */
  rhea_stokes_problem_destroy (stokes_problem);

  /* destroy vectors */
  if (temperature != NULL) {
    rhea_temperature_destroy (temperature);
  }
  if (compositional_density != NULL) {
    rhea_composition_destroy (compositional_density);
  }
  if (compositional_viscosity != NULL) {
    rhea_composition_destroy (compositional_viscosity);
  }

  /* destroy data */
  rhea_plate_data_clear (plate_options);
  rhea_weakzone_data_clear (weak_options);

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);
}
