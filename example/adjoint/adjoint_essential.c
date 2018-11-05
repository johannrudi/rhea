#include <adjoint_essential.h>

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

void
adjoint_setup_clear_all (rhea_stokes_problem_t *stokes_problem,
                   rhea_newton_problem_t *newton_problem,
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
  const char         *this_fn_name = "adjoint_setup_clear_all";
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

  adjoint_problem_destroy (newton_problem);
  example_share_stokes_destroy (stokes_problem, temp_options, plate_options, weak_options,
                                visc_options);

  /* destroy mesh */
  example_share_mesh_destroy (ymir_mesh, press_elem, p4est, topo_options,
                              discr_options);
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}


