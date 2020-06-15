#include <subduction_essential.h>


/************************************************************************************************
 * Sets up the mesh.
 ************************************************************************************************/

void
subd_setup_mesh (p4est_t **p4est,
                    ymir_mesh_t **ymir_mesh,
                    ymir_pressure_elem_t **press_elem,
                    MPI_Comm mpicomm,
                    rhea_domain_options_t *domain_options,
                    rhea_topography_options_t *topo_options,
                    rhea_discretization_options_t *discr_options,
                    subd_options_t *subd_options)

{
  const char         *this_fn_name = "subd_setup_mesh";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* set up data */
  rhea_topography_data_create (topo_options, mpicomm);

  /* create p4est */
  *p4est = rhea_discretization_p4est_new (mpicomm, discr_options,
                                          domain_options);

  /* set up boundary, store in `discr_options` */
  rhea_discretization_boundary_create (discr_options, *p4est,
                                       domain_options);

  /* create ymir mesh and pressure element */
  rhea_discretization_ymir_mesh_new_from_p4est (ymir_mesh, press_elem, *p4est,
                                                discr_options);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}


/*********************************************************************************************
 * Sets up a linear Stokes problem.
 *********************************************************************************************/
/*Basic Stokes problem */
void
subd_setup_stokes (rhea_stokes_problem_t **stokes_problem,
                    p4est_t     *p4est,
                    ymir_mesh_t **ymir_mesh,
                    ymir_pressure_elem_t **press_elem,
                    rhea_discretization_options_t *discr_options,
                    rhea_domain_options_t *domain_options,
                    rhea_temperature_options_t *temp_options,
                    rhea_weakzone_options_t *weak_options,
                    rhea_viscosity_options_t *visc_options,
					rhea_composition_options_t *comp_options,
                    subd_options_t *subd_options,
                    const char *vtk_write_input_path,
                    const char *solver_bin_path,
                    const char *solver_vtk_path)
{
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (*ymir_mesh);
  ymir_vec_t         *temperature;
  ymir_vec_t		 *composition;

  RHEA_GLOBAL_PRODUCTION_FN_BEGIN (__func__);

  /* set up data */
  rhea_weakzone_data_create (weak_options, mpicomm);

  /* compute temperature */
  temperature = rhea_temperature_new (*ymir_mesh);
  subd_compute_temperature (temperature, temp_options, subd_options);

  /* read composition */
  composition = rhea_composition_new (*ymir_mesh);
  rhea_composition_read (composition, comp_options);

  subd_set_velocity_dirichlet_bc (domain_options, subd_options);

  /* create Stokes problem */
  *stokes_problem = rhea_stokes_problem_new (
      *ymir_mesh, *press_elem, temperature, composition, domain_options,
      temp_options, comp_options, weak_options, visc_options);

  /* destroy vector composition */
  rhea_composition_destroy (composition);

  /*call back to recompute rhs_vel, rhs_velbc_nonzero_dir, rhs_velbc_nonzero_neumann, weakzone*/
  subd_compute_rhs_vel (*stokes_problem, subd_options);
  subd_compute_rhs_velbc_dirichlet (*stokes_problem, subd_options);
  subd_compute_rhs_velbc_neumann (*stokes_problem, subd_options);
  subd_compute_weakzone (*stokes_problem, subd_options);

  /* reset function of viscosity */
  subd_viscosity_set_function (*stokes_problem, subd_options);

  rhea_stokes_problem_set_solver_bin_output (*stokes_problem, solver_bin_path);
  rhea_stokes_problem_set_solver_vtk_output (*stokes_problem, solver_vtk_path);

  /* about amr */
  {
    rhea_stokes_problem_set_solver_amr (*stokes_problem, p4est, discr_options);
    /* perform initial AMR */
    if (p4est != NULL && discr_options != NULL) {
      rhea_stokes_problem_init_amr (*stokes_problem, p4est, discr_options);

      /* retrieve adapted mesh */
      *ymir_mesh = rhea_stokes_problem_get_ymir_mesh (*stokes_problem);
      *press_elem = rhea_stokes_problem_get_press_elem (*stokes_problem);
    }
  }

  /* set up Stokes solver */
  rhea_stokes_problem_setup_solver (*stokes_problem);

    /* write vtk of problem input */
  if (vtk_write_input_path != NULL) {
    subd_write_input_basic (*stokes_problem, temp_options, visc_options,
                            vtk_write_input_path);
  }

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);
}

/********************************************************************************************
 * Cleans up Stokes problem and mesh.
 ********************************************************************************************/

void
subd_setup_clear_all (rhea_stokes_problem_t *stokes_problem,
                         p4est_t *p4est,
                         ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         rhea_temperature_options_t *temp_options,
                         rhea_viscosity_options_t *visc_options,
                         rhea_weakzone_options_t *weak_options,
                         rhea_topography_options_t *topo_options,
                         rhea_plate_options_t *plate_options,
                         rhea_discretization_options_t *discr_options)
{
  const char         *this_fn_name = "subd_setup_clear_all";
  ymir_vec_t         *visc_TI_svisc;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* destroy anisotropic viscosity */
  {
    ymir_stokes_op_t    *stokes_op;
    ymir_vec_t          *visc_TI_svisc;

    stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
    visc_TI_svisc = stokes_op->stress_op->coeff_TI_svisc;

    if (visc_TI_svisc != NULL) {
      rhea_viscosity_destroy (visc_TI_svisc);
    }
  }

  example_share_stokes_destroy (stokes_problem, temp_options, NULL /* comp_options */,
                                plate_options, weak_options, visc_options);

  /* destroy mesh */
  example_share_mesh_destroy (ymir_mesh, press_elem, p4est, topo_options,
                              discr_options);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**************************************************************************************************
 * runs stokes solver.
 *************************************************************************************************/

void
subd_run_solver (ymir_vec_t *sol_vel_press,
                    ymir_mesh_t *ymir_mesh,
                    ymir_pressure_elem_t *press_elem,
                    rhea_stokes_problem_t *stokes_problem,
                    const int nonzero_initial_guess,
                    const int iter_max, const double rel_tol)
{
  const char         *this_fn_name = "subd_run_solver";
  ymir_vec_t         *rhs_vel_nonzero_dirichlet;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* run solver */
  rhea_stokes_problem_solve (&sol_vel_press, nonzero_initial_guess,
                iter_max, rel_tol, stokes_problem);

  /* add nonzero dirichlet values of the velocity to the solution */
  rhs_vel_nonzero_dirichlet =
    rhea_stokes_problem_get_rhs_vel_nonzero_dirichlet (stokes_problem);
  if (rhs_vel_nonzero_dirichlet != NULL) {
    ymir_vec_t         *sol_vel = rhea_velocity_new (ymir_mesh);

    ymir_stokes_vec_get_velocity (sol_vel_press, sol_vel, press_elem);
    ymir_vec_add (1.0, rhs_vel_nonzero_dirichlet, sol_vel);
    ymir_stokes_vec_set_velocity (sol_vel, sol_vel_press, press_elem);
    rhea_velocity_destroy (sol_vel);
  }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}


