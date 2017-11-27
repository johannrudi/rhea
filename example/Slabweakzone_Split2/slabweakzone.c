/* Slabs_Weakzone Example:
 *
 * Cartesian domain.
*/

#include <rhea.h>
#include <ymir_velocity_vec.h>
#include <ymir_derivative_elem.h>
#include <ymir_vec_ops.h>
#include <ymir_stokes_pc.h>
#include <ymir_stokes_vec.h>
#include <ymir_stress_op.h>
#include <ymir_pressure_vec.h>

/**
 * Sets up the mesh.
 **/
static void
slabs_setup_mesh (p4est_t **p4est,
                    ymir_mesh_t **ymir_mesh,
                    ymir_pressure_elem_t **press_elem,
                    MPI_Comm mpicomm,
                    rhea_domain_options_t *domain_options,
                    rhea_discretization_options_t *discr_options,
                    slabs_options_t *slabs_options)

{
  const char         *this_fn_name = "slabs_setup_mesh";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

//  slabs_surface_location (slabs_options, discr_options);

  /* create p4est */
  *p4est = rhea_discretization_p4est_new (mpicomm, discr_options,
                                          domain_options);

  /* set up boundary, store in `discr_options` */
  rhea_discretization_options_set_boundary (discr_options, *p4est,
                                            domain_options);

  /* create ymir mesh and pressure element */
  rhea_discretization_ymir_mesh_new_from_p4est (ymir_mesh, press_elem, *p4est,
                                                discr_options);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Sets up a linear Stokes problem.
 */
static void
slabs_setup_stokes (rhea_stokes_problem_t **stokes_problem,
                    ymir_mesh_t *ymir_mesh,
                    ymir_pressure_elem_t *press_elem,
                    rhea_domain_options_t *domain_options,
                    rhea_temperature_options_t *temp_options,
                    rhea_viscosity_options_t *visc_options,
                    slabs_options_t *slabs_options,
                    const char *vtk_write_input_path)
{
  const char         *this_fn_name = "slabs_setup_stokes";
  ymir_vec_t         *temperature, *weakzone;
  ymir_vec_t         *coeff_TI_svisc=NULL, *TI_rotate=NULL;
  ymir_vec_t         *rhs_vel, *rhs_vel_nonzero_dirichlet=NULL;
  void               *solver_options = NULL;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* compute temperature */
  temperature = rhea_temperature_new (ymir_mesh);
  /* compute weak zone */
  weakzone = rhea_viscosity_new (ymir_mesh);
  ymir_vec_set_value (weakzone, 1.0);

  switch (slabs_options->buoyancy_type)  {
  /* if non-specified, use: rhea_temperature_compute (temperature, temp_options); */
    case SINKER:
      rhea_temperature_compute (temperature, temp_options);
    break;

    case SLAB:
      slabs_temperature_compute (temperature, slabs_options);
      if (slabs_options->slabs_visc_options->viscosity_anisotropy
         == SLABS_VISC_ISOTROPY) {
        slabs_weakzone_compute (weakzone, slabs_options);
      }
    break;

    case COLLIDE:
      rhea_temperature_compute (temperature, temp_options);
      if (slabs_options->slabs_visc_options->viscosity_anisotropy
        == SLABS_VISC_ISOTROPY) {
        slabs_weakzone_compute (weakzone, slabs_options);
      }
      /* set custom function to compute viscosity */
      rhea_viscosity_set_viscosity_compute_fn (slabs_viscosity_compute,
                                               slabs_options);
    break;

    case DRAG:
      ymir_cvec_set_function (temperature, drag_temperature_set_fn, slabs_options);
      if (slabs_options->slabs_visc_options->viscosity_anisotropy
        == SLABS_VISC_ISOTROPY) {
        slabs_weakzone_compute (weakzone, slabs_options);
      }
      /* set custom function to compute viscosity */
      rhea_viscosity_set_viscosity_compute_fn (slabs_viscosity_compute,
                                               slabs_options);
    break;

    case TEST_MANUFACTURED:
      RHEA_GLOBAL_PRODUCTIONF ("buoyancy from %i\n", slabs_options->buoyancy_type);
      rhea_temperature_compute (temperature, temp_options);  // neutral value: T=0.5
      rhea_viscosity_set_viscosity_compute_fn (slabs_viscosity_compute,
                                               slabs_options);
    break;

    case TESTTOPO:
      ymir_cvec_set_function (temperature, testtopo_temperature_set_fn, slabs_options);
    break;

    case TESTTOPO2:
      ymir_cvec_set_function (temperature, testtopo2_temperature_set_fn, slabs_options->slabs_surf_options->topo_profile);
    break;

    case TESTNONE:
      ymir_vec_set_value (temperature, 0.5);
    break;

    default:
      RHEA_ABORT_NOT_REACHED ();
  }

  rhs_vel = rhea_velocity_new (ymir_mesh);
  /* for the test using manufactured solution,
   * overwrite rhs_vel with estimated forcing term from given velocity and pressure field.
   * don't apply mass now */
  if (slabs_options->buoyancy_type == TEST_MANUFACTURED ||
      slabs_options->slabs_test_options->test_manufactured != SLABS_TEST_MANUFACTURED_NONE) {
    slabs_test_manufactured_rhs_compute (rhs_vel, slabs_options);
  }
  else {
    /* compute velocity right-hand side volume forcing */
    rhea_temperature_compute_rhs_vel (rhs_vel, temperature, temp_options);
  }

  /* set velocity boundary conditions & nonzero Dirichlet values */
  if (domain_options->velocity_bc_type == RHEA_DOMAIN_VELOCITY_BC_USER) {
    rhs_vel_nonzero_dirichlet = rhea_velocity_new (ymir_mesh);
    if (slabs_options->buoyancy_type == TEST_MANUFACTURED ||
        slabs_options->slabs_test_options->test_manufactured != SLABS_TEST_MANUFACTURED_NONE) {
      slabs_test_manufactured_velbc_compute (rhs_vel_nonzero_dirichlet,
                                             slabs_options);
    }
    else {
      slabs_vel_nonzero_dirichlet_compute (rhs_vel_nonzero_dirichlet,
                                           slabs_options);
    }
  }

  /* create Stokes problem */
  *stokes_problem = rhea_stokes_problem_new (
      temperature, weakzone, rhs_vel, rhs_vel_nonzero_dirichlet,
      ymir_mesh, press_elem, domain_options, visc_options, solver_options);

  /* add the anisotropic viscosity to the viscous stress operator */
  if (slabs_options->slabs_visc_options->viscosity_anisotropy
      == SLABS_VISC_TRANSVERSELY_ISOTROPY) {
    coeff_TI_svisc = rhea_viscosity_new (ymir_mesh);
    TI_rotate = rhea_viscosity_new (ymir_mesh);
    if (slabs_options->buoyancy_type == TEST_MANUFACTURED ||
        slabs_options->slabs_test_options->test_manufactured != SLABS_TEST_MANUFACTURED_NONE) {
      slabs_stokes_problem_setup_TI_manufactured (ymir_mesh, *stokes_problem, slabs_options,
                                     coeff_TI_svisc, TI_rotate);
    }
    else {
      slabs_stokes_problem_setup_TI (ymir_mesh, *stokes_problem, slabs_options,
                                     coeff_TI_svisc, TI_rotate);
    }
  }

  /* write vtk of problem input */
  if (vtk_write_input_path != NULL) {
    slabs_write_input (ymir_mesh, *stokes_problem, temp_options,
                       temperature, weakzone, coeff_TI_svisc, TI_rotate,
                       vtk_write_input_path);
  }

  /* set up Stokes solver */
  rhea_stokes_problem_setup_solver (*stokes_problem);

  /* destroy */
  if (slabs_options->slabs_visc_options->viscosity_anisotropy
      == SLABS_VISC_TRANSVERSELY_ISOTROPY)  {
      rhea_viscosity_destroy (TI_rotate);
  }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Cleans up Stokes problem and mesh.
 */
static void
slabs_setup_clear_all (rhea_stokes_problem_t *stokes_problem,
                         p4est_t *p4est,
                         ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         rhea_discretization_options_t *discr_options)
{
  const char         *this_fn_name = "slabs_setup_clear_all";
  ymir_vec_t         *temperature, *weakzone;
  ymir_vec_t         *visc_TI_svisc;
  ymir_vec_t         *rhs_vel, *rhs_vel_nonzero_dirichlet;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* get vectors */
  temperature = rhea_stokes_problem_get_temperature (stokes_problem);
  weakzone = rhea_stokes_problem_get_weakzone (stokes_problem);
  rhs_vel = rhea_stokes_problem_get_rhs_vel (stokes_problem);
  rhs_vel_nonzero_dirichlet =
    rhea_stokes_problem_get_rhs_vel_nonzero_dirichlet (stokes_problem);

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

  /* destroy mesh */
  rhea_discretization_ymir_mesh_destroy (ymir_mesh, press_elem);
  rhea_discretization_p4est_destroy (p4est);

  /* destroy (some) options */
  rhea_discretization_options_clear (discr_options);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * runs stokes solver.
 */
static void
slabs_run_solver (ymir_vec_t *sol_vel_press,
                    ymir_mesh_t *ymir_mesh,
                    ymir_pressure_elem_t *press_elem,
                    rhea_stokes_problem_t *stokes_problem,
                    const int iter_max, const double rel_tol)
{
  const char         *this_fn_name = "slabs_run_solver";
  ymir_vec_t         *rhs_vel_nonzero_dirichlet;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* run solver */
  rhea_stokes_problem_solve (sol_vel_press, iter_max, rel_tol, stokes_problem);

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

/**
 * Runs the program.
 */
int
main (int argc, char **argv)
{
  const char         *this_fn_name = "slabweakzone:main";
  /* MPI */
  MPI_Comm            mpicomm = MPI_COMM_WORLD;
  int                 mpisize, mpirank, ompsize;
  int                 mpiret;
  /* options */
  ymir_options_t     *opt;
  rhea_domain_options_t         domain_options;
  rhea_temperature_options_t    temp_options;
  rhea_viscosity_options_t      visc_options;
  rhea_discretization_options_t discr_options;
  rhea_newton_options_t         newton_options;

  int                     buoyancy_type;
  int                     viscosity_anisotropy;
  int                     x_func;
  int                     vel_dir_bc;
  int                     test_manufactured;
  int                     test_stress_op;
  int                     test_stress_comp;

  /* slabs options */
  slabs_domain_options_t     slabs_domain_options;
  slabs_temp_options_t     slabs_temp_options;
  slabs_visc_options_t     slabs_visc_options;
  slabs_weak_options_t     slabs_weak_options;
  slabs_surf_options_t     slabs_surf_options;
  slabs_velbc_options_t    slabs_velbc_options;
  slabs_test_options_t    slabs_test_options;
  slabs_options_t          slabs_options;

  /* options local to this function */
  int                 production_run;
  int                 solver_iter_max;
  double              solver_rel_tol;
  char               *vtk_write_input_path;
  char               *vtk_write_solution_path;
  char               *vtk_write_stress_path;
  char               *vtk_write_postp_path;
  char               *vtk_write_test_path;
  char               *vtk_write_freesurface_path;
  char               *vtk_write_input2_path;
  char               *vtk_write_solution2_path;
  char               *vtk_write_input3_path;
  char               *vtk_write_solution3_path;

  /* mesh */
  p4est_t            *p4est;
  ymir_mesh_t        *ymir_mesh;
  ymir_pressure_elem_t  *press_elem;
  /* Stokes */
  rhea_stokes_problem_t *stokes_problem;
  ymir_vec_t         *sol_vel_press;

  /*
   * Initialize Libraries
   */

  /* initialize rhea and sub-packages */
  rhea_initialize (argc, argv, mpicomm);

  /* get parallel environment */
  mpiret = MPI_Comm_size (mpicomm, &mpisize); YMIR_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &mpirank); YMIR_CHECK_MPI (mpiret);

#ifdef RHEA_ENABLE_OPENMP
  ompsize = omp_get_max_threads ();
#else
  ompsize = 1;
#endif

  /*
   * Define & Parse Options
   */

  opt = ymir_options_global_new (argv[0] /* program path */);

  /*add slabweakzone options*/
  slabweakzone_options (opt);

  /* parse options */
  {
    int                 optret;

    optret = ymir_options_parse (SC_LP_INFO, opt, argc, argv);
    if (optret < 0) { /* if parsing was not successful */
      ymir_options_print_usage (SC_LP_INFO, opt, NULL /* args usage */);
      RHEA_GLOBAL_INFO ("Option parsing failed\n");
      exit (0);
    }

  /*
   * Process Slabs Options
   */
  slabweakzone_options_process (&slabs_options,
                                &slabs_temp_options, &slabs_visc_options,
                                &slabs_weak_options, &slabs_surf_options,
                                &slabs_velbc_options, &slabs_test_options);

  /*
   * Initialize Main Program
   */

  RHEA_GLOBAL_PRODUCTIONF (
      "Into %s (production %i)\n", this_fn_name, production_run);
  RHEA_GLOBAL_PRODUCTIONF (
      "Parallel environment: MPI size %i, OpenMP size %i\n", mpisize, ompsize);
  ymir_set_up (argc, argv, mpicomm, production_run);

  /* print & process options */
  ymir_options_print_summary (SC_LP_INFO, opt);
  rhea_process_options_all (&domain_options, &temp_options,
                            &visc_options, &discr_options,
                            &newton_options);

  /* copy rhea domain options into local example domain options */
  slabweakzone_options_domain_opt (&slabs_options,
                                   &slabs_domain_options, &domain_options);


  /*
   * Setup Mesh
   */

  rhea_discretization_set_user_X_fn (&discr_options,
                                     slabs_X_fn_identity, NULL);
  slabs_setup_mesh (&p4est, &ymir_mesh, &press_elem, mpicomm,
                      &domain_options, &discr_options, &slabs_options);

  /*
   * Setup Stokes Problem
   */

  slabs_setup_stokes (&stokes_problem, ymir_mesh, press_elem,
                        &domain_options, &temp_options, &visc_options,
                        &slabs_options, vtk_write_input_path);

  /*
   * Solve Stokes Problem
   */

  /* initialize solution vector */
  sol_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);

  /* run solver */
  slabs_run_solver (sol_vel_press, ymir_mesh, press_elem, stokes_problem,
                      solver_iter_max, solver_rel_tol);

  /* write vtk of solution */
  if (vtk_write_solution_path != NULL) {
    ymir_vec_t         *velocity = rhea_velocity_new (ymir_mesh);
    ymir_vec_t         *pressure = rhea_pressure_new (ymir_mesh, press_elem);
    ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);

    ymir_stokes_vec_get_components (sol_vel_press, velocity, pressure,
                                    press_elem);
    rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

    rhea_vtk_write_solution (vtk_write_solution_path, velocity, pressure,
                             viscosity);

    rhea_pressure_destroy (pressure);
    rhea_viscosity_destroy (viscosity);
    rhea_velocity_destroy (velocity);
  }

  if (vtk_write_test_path != NULL)  {
    if(slabs_test_options.test_manufactured != SLABS_TEST_MANUFACTURED_NONE) {
      ymir_vec_t         *vel_ref = rhea_velocity_new (ymir_mesh);
      ymir_vec_t         *vel_chk = rhea_velocity_new (ymir_mesh);
      ymir_vec_t         *pres_ref = ymir_cvec_new (ymir_mesh, 1);
      ymir_vec_t         *pres_chk = rhea_pressure_new (ymir_mesh, press_elem);
      double              vel_abs_err, vel_rel_err;
      ymir_vec_t         *vel_err = rhea_velocity_new (ymir_mesh);
      ymir_vec_t         *pres_err = ymir_cvec_new (ymir_mesh, 1);
      const               slabs_test_manufactured_t
                          test_type = slabs_test_options.test_manufactured;
      ymir_stokes_op_t   *stokes_op;
      ymir_stress_op_t   *stress_op;

      stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
      stress_op = stokes_op->stress_op;
      ymir_vec_set_value (pres_ref, .0);

      /* compute velocity fields */
      switch (test_type) {
        case SLABS_TEST_MANUFACTURED_SINCOS1_ISO:
        case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT90:
        case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT45:
        case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT60:
        case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT60_VISCEXP60:
          /* compute reference velocity field (output) */
          ymir_cvec_set_function (vel_ref, slabs_test_sincos1_vel_in_fn,
                                  NULL);
          break;

        case SLABS_TEST_MANUFACTURED_POLY1_TIROT90:
        case SLABS_TEST_MANUFACTURED_POLY1_TIROT90_VISCEXP:
          /* compute reference velocity field (output) */
          ymir_cvec_set_function (vel_ref, slabs_test_poly1_vel_in_fn,
                                  NULL);

          break;

       default:
          RHEA_ABORT_NOT_REACHED ();
      }

      ymir_stokes_vec_get_components (sol_vel_press, vel_chk, pres_chk,
                                        press_elem);
      slabs_test_manufactured_compute_vel_err (&vel_abs_err, &vel_rel_err,
                                                 vel_err, vel_ref, vel_chk,
                                                 stress_op);

      RHEA_GLOBAL_INFOF ("manufactured solution test (test type %i): abs error %1.3e rel error %1.3e\n",
                          test_type, vel_abs_err, vel_rel_err);


      char      path[BUFSIZ];
      snprintf (path, BUFSIZ, "%s_manufactured", vtk_write_test_path);
      ymir_vtk_write (ymir_mesh, path,
                      vel_ref, "vel_reference",
                      vel_chk, "vel_check",
                      vel_err, "vel_error",
                      pres_ref, "pressure_reference",
                      pres_chk, "pressure_check",
                      NULL);
      rhea_velocity_destroy (vel_ref);
      rhea_velocity_destroy (vel_chk);
      rhea_velocity_destroy (vel_err);
      rhea_pressure_destroy (pres_chk);
      ymir_vec_destroy (pres_ref);
      ymir_vec_destroy (pres_err);
    }

    if (slabs_test_options.test_stress_comp != SLABS_TEST_MANUFACTURED_NONE)
    {
      const               slabs_test_manufactured_t
                          test_type = slabs_test_options.test_stress_comp;

      ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (ymir_mesh->ma->N,
                                                                ymir_mesh->ma->ompsize);
      ymir_vec_t         *velocity = rhea_velocity_new (ymir_mesh);
      ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);
      ymir_vec_t         *shear_visc = rhea_viscosity_new (ymir_mesh);
      ymir_vec_t         *TI_tensor = ymir_dvec_new (ymir_mesh, 9, YMIR_GAUSS_NODE);
      ymir_vec_t         *trac_n = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
      ymir_vec_t         *trac_s = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
      ymir_vec_t         *normal_force = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
      ymir_vec_t         *shear_force = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
      ymir_vec_t         *traction_ref = rhea_velocity_new (ymir_mesh);
      ymir_vec_t         *stress = ymir_dvec_new (ymir_mesh, 3, YMIR_GAUSS_NODE);
      ymir_vec_t         *stress_ref = rhea_velocity_new (ymir_mesh);
      RHEA_GLOBAL_PRODUCTIONF ("In %s: Start test stress_component\n", this_fn_name);

      ymir_stokes_vec_get_velocity (sol_vel_press, velocity,
                                      press_elem);
      rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

      if (slabs_visc_options.viscosity_anisotropy == SLABS_VISC_TRANSVERSELY_ISOTROPY)  {
        ymir_stokes_op_t      *stokes_op;
        ymir_stress_op_t      *stress_op;

        /* get the viscous stress operator */
        stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
        stress_op = stokes_op->stress_op;

        /* copy shear viscosity */
        slabs_stress_op_copy_shear_visc (shear_visc, stress_op);
        slabs_stress_op_copy_TI_tensor (TI_tensor, stress_op);

        slabs_manufactured_strainvec_coupling_compute (velocity, trac_n, trac_s,
                                               shear_visc, viscosity);

        slabs_manufactured_stressvec_coupling_compute (velocity, viscosity,
                                              shear_visc, TI_tensor,
                                              normal_force, shear_force,
                                              vel_elem, &slabs_options);

        slabs_stressvec_TI (velocity, stress, viscosity, shear_visc, TI_tensor, vel_elem);
        switch (test_type) {
          case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT90:
            ymir_cvec_set_function (traction_ref, slabs_test_sincos1_TIrot90_traction_fn,
                                      NULL);
            ymir_cvec_set_function (stress_ref, slabs_test_sincos1_TIrot90_stress_fn,
                                      NULL);
          case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT45:
            ymir_cvec_set_function (traction_ref, slabs_test_sincos1_TIrot45_traction_fn,
                                      NULL);
            ymir_cvec_set_function (stress_ref, slabs_test_sincos1_TIrot45_stress_fn,
                                      NULL);
          case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT60:
            ymir_cvec_set_function (traction_ref, slabs_test_sincos1_TIrot60_traction_fn,
                                      NULL);
            ymir_cvec_set_function (stress_ref, slabs_test_sincos1_TIrot60_stress_fn,
                                      NULL);
            break;

          default:
            RHEA_ABORT_NOT_REACHED ();
        }
      }
      char      path[BUFSIZ];
      snprintf (path, BUFSIZ, "%s_stress_component", vtk_write_test_path);
      ymir_vtk_write (ymir_mesh, path,
                      stress_ref, "stress_reference (yy,zz,yz)",
                      stress, "stress_check",
                      traction_ref, "traction_reference (0,fn,fs)",
                      normal_force, "normal_force (from stress)",
                      shear_force, "shear_force (from stress)",
                      trac_n, "normal_force (from strain)",
                      trac_s, "shear_force (from strain)",
                      NULL);

      ymir_velocity_elem_destroy (vel_elem);
      rhea_velocity_destroy (velocity);
      rhea_viscosity_destroy (viscosity);
      rhea_viscosity_destroy (shear_visc);
      ymir_vec_destroy (TI_tensor);
      ymir_vec_destroy (trac_n);
      ymir_vec_destroy (trac_s);
      ymir_vec_destroy (normal_force);
      ymir_vec_destroy (shear_force);
      ymir_vec_destroy (traction_ref);
      ymir_vec_destroy (stress);
      ymir_vec_destroy (stress_ref);
    } /* end of stress component*/
    /*stress operator test, TODO*/
  }

  /* compute and output second invariant strain_rate, stress, and surface normal stress  */
  if (vtk_write_stress_path != NULL)  {
    ymir_vec_t         *velocity = rhea_velocity_new (ymir_mesh);
    ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);
    ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (ymir_mesh->ma->N,
                                                              ymir_mesh->ma->ompsize);
    ymir_vec_t            *edotII = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
    ymir_vec_t            *tauII = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
    ymir_vec_t            *surf_normal_stress = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);

    RHEA_GLOBAL_PRODUCTIONF ("In %s: Start vtk_write_stress\n", this_fn_name);

    ymir_stokes_vec_get_velocity (sol_vel_press, velocity,
                                    press_elem);
    rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

    /* compute 2nd invariant of the strain rate */
    ymir_second_invariant_vec (velocity, edotII, vel_elem);
    ymir_vec_sqrt (edotII, edotII);

    /* compute 2nd invariant of deviatoric stress tau = 2* (2nd invariant of strain_rate * viscosity )
      and its projection on the surface */
    if (slabs_visc_options.viscosity_anisotropy == SLABS_VISC_TRANSVERSELY_ISOTROPY)  {
      ymir_stokes_op_t      *stokes_op;
      ymir_stress_op_t      *stress_op;
      ymir_vec_t            *shear_visc = rhea_viscosity_new (ymir_mesh);
      ymir_vec_t            *TI_tensor = ymir_dvec_new (ymir_mesh, 9,
                                                        YMIR_GAUSS_NODE);

      /* get the viscous stress operator */
      stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
      stress_op = stokes_op->stress_op;

      /* copy shear viscosity */
      slabs_stress_op_copy_shear_visc (shear_visc, stress_op);
      slabs_stress_op_copy_TI_tensor (TI_tensor, stress_op);

      slabs_2inv_stress_TI (velocity, tauII,
                              viscosity, shear_visc, TI_tensor, vel_elem);

      rhea_viscosity_destroy (shear_visc);
      ymir_vec_destroy (TI_tensor);
    }
    else  {
      ymir_vec_copy (edotII, tauII)
      ymir_vec_multiply_in (viscosity, tauII);
      ymir_vec_scale (2.0, tauII);
    }

    /* compute surface normal stress sigma */
    slabs_physics_compute_normal_boundary_stress (
                   surf_normal_stress, sol_vel_press,
                   rhea_stokes_problem_get_rhs_vel (stokes_problem),
                   rhea_stokes_problem_get_stokes_op (stokes_problem));

    {
      char            path[BUFSIZ];

      snprintf (path, BUFSIZ, "%s", vtk_write_stress_path);
      ymir_vtk_write (ymir_mesh, path,
                      edotII, "edotII",
                      tauII, "tauII",
                      surf_normal_stress, "surf_normal_stress",
                      NULL);
    }

    /* destroy */
    rhea_viscosity_destroy (viscosity);
    rhea_velocity_destroy (velocity);
    ymir_vec_destroy (edotII);
    ymir_vec_destroy (tauII);
    ymir_vec_destroy (surf_normal_stress);
    ymir_velocity_elem_destroy (vel_elem);

    RHEA_GLOBAL_PRODUCTIONF ("In %s: Done vtk_write_stress\n", this_fn_name);
  }

    /* compute and output analysis of stress */
  if (vtk_write_postp_path != NULL)  {
    ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (ymir_mesh->ma->N,
                                                              ymir_mesh->ma->ompsize);
    ymir_vec_t         *velocity = rhea_velocity_new (ymir_mesh);
    ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *shear_visc = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *TI_tensor = ymir_dvec_new (ymir_mesh, 9, YMIR_GAUSS_NODE);
    ymir_vec_t         *traction = ymir_dvec_new (ymir_mesh, 3, YMIR_GAUSS_NODE);
    ymir_vec_t         *normal_force = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
    ymir_vec_t         *shear_force = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
    RHEA_GLOBAL_PRODUCTIONF ("In %s: Start vtk_write_postp\n", this_fn_name);

    ymir_stokes_vec_get_velocity (sol_vel_press, velocity,
                                    press_elem);
    rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

    if (slabs_visc_options.viscosity_anisotropy == SLABS_VISC_TRANSVERSELY_ISOTROPY)  {
      ymir_stokes_op_t      *stokes_op;
      ymir_stress_op_t      *stress_op;

      /* get the viscous stress operator */
      stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
      stress_op = stokes_op->stress_op;

      /* copy shear viscosity */
      slabs_stress_op_copy_shear_visc (shear_visc, stress_op);
      slabs_stress_op_copy_TI_tensor (TI_tensor, stress_op);

      slabs_postp_weakzone_traction_compute (velocity, traction,
                                             shear_visc, viscosity,
                                             &slabs_options);

      slabs_postp_weakzone_coupling_compute (velocity, viscosity,
                                            shear_visc, TI_tensor,
                                            normal_force, shear_force,
                                            vel_elem, &slabs_options);
    }
    else  {
      slabs_postp_weakzone_traction_compute (velocity, traction,
                                             viscosity, viscosity,
                                             &slabs_options);

      slabs_postp_weakzone_coupling_compute (velocity, viscosity,
                                             shear_visc, TI_tensor,
                                             normal_force, shear_force,
                                             vel_elem, &slabs_options);
   }
   rhea_viscosity_destroy (shear_visc);
   ymir_vec_destroy (TI_tensor);

   {
     char            path[BUFSIZ];

     snprintf (path, BUFSIZ, "%s", vtk_write_postp_path);
     ymir_vtk_write (ymir_mesh, path,
                     normal_force, "normal_force from sigma",
                     shear_force, "shear_force from sigma",
                     traction, "trac_n from edot",
                     NULL);
   }

    /* destroy */
    rhea_viscosity_destroy (viscosity);
    rhea_velocity_destroy (velocity);
    ymir_vec_destroy (normal_force);
    ymir_vec_destroy (shear_force);
    ymir_vec_destroy (traction);
    ymir_velocity_elem_destroy (vel_elem);

    RHEA_GLOBAL_PRODUCTIONF ("In %s: Done vtk_write_postp\n", this_fn_name);
  }

  /* compute and output second invariant strain_rate, stress, and surface normal stress  */
  if (vtk_write_freesurface_path != NULL)  {
    ymir_vec_t            *surf_normal_stress = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    ymir_vec_t            *surf_normal_stress2;
    p4est_t               *p4est2;
    ymir_mesh_t           *ymir_mesh2;
    ymir_pressure_elem_t  *press_elem2;
    ymir_vec_t            *sol_vel_press2;
    rhea_stokes_problem_t *stokes_problem2;
    rhea_discretization_options_t discr_options2;

    ymir_vec_t            *surf_normal_stress3;
    ymir_vec_t            *sol_vel_press3;

    ymir_topidx_t         fm;
    ymir_face_mesh_t      *fmesh;
    mangll_cnodes_t       *cnodes = ymir_mesh->cnodes;
    const int             N = cnodes->N;
    const int             Nrp = N + 1;
    int                   N3;
    int                   Np;
    int                   Ntotal;
    int                   K;
    int                   i, il, ik;
    sc_dmatrix_t          *elem;
    double                *elemd;
    double                *tX, *tY, *tZ;
    double                *Xd, *Yd, *Zd;
    double                avg_stress, topo_nondim;
    slabs_topo_profile_t  topo;

    RHEA_GLOBAL_PRODUCTIONF ("In %s: Start vtk_write_freesurface\n", this_fn_name);

    /* compute surface normal stress sigma */
    slabs_physics_compute_normal_boundary_stress (
                   surf_normal_stress, sol_vel_press,
                   rhea_stokes_problem_get_rhs_vel (stokes_problem),
                   rhea_stokes_problem_get_stokes_op (stokes_problem));

    fm = surf_normal_stress->meshnum;
    fmesh = &(ymir_mesh->fmeshes[fm]);
    N3 = N * N;
    Np = (N + 1) * (N + 1);
    K = fmesh->K;
    Ntotal = Np * K;
    Xd = fmesh->X->e[0];
    Yd = fmesh->Y->e[0];
    Zd = fmesh->Z->e[0];

    tX = (double *) malloc(Ntotal * sizeof(double));
    tY = (double *) malloc(Ntotal * sizeof(double));
    tZ = (double *) malloc(Ntotal * sizeof(double));

    for (il = 0; il < Ntotal; ++il)  {
      tX[il] = Xd[il];
      tY[il] = Yd[il];
    }

    elem = sc_dmatrix_new (Np, 1);
    avg_stress = 0.0;
    for (ik = 0; ik < K; ik++)  {
      ymir_vec_get_elem_interp (surf_normal_stress, elem, YMIR_STRIDE_NODE, ik, YMIR_GLL_NODE, YMIR_COPY);
      elemd = elem->e[0];
      for (i = 0; i < Np; i++)  {
        avg_stress += elemd[i];
      }
    }
    avg_stress /= Ntotal;
    avg_stress = 0.0;

    for (ik = 0; ik < K; ik++)  {
      ymir_vec_get_elem_interp (surf_normal_stress, elem, YMIR_STRIDE_NODE, ik, YMIR_GLL_NODE, YMIR_COPY);
      elemd = elem->e[0];
      for (i = 0; i < Np; i++)  {
        topo_nondim = (elemd[i] - avg_stress) * (-0.084/7.0e7);
        tZ[Np * ik + i] =  topo_nondim + 1.0;

        RHEA_GLOBAL_INFOF ("element %d, nodeid %d, topo[%d]=%lf\n",ik, i, Np*ik+i, tZ[Np * ik + i]);
      }
    }


    topo.tX = tX;
    topo.tY = tY;
    topo.tZ = tZ;
    topo.nsurf = Ntotal;
    slabs_surf_options.topo_profile = &topo;
    rhea_discretization_process_options (&discr_options2, &domain_options);
    rhea_discretization_set_user_X_fn (&discr_options2,
                                       slabs_X_fn_profile, &topo);

    slabs_setup_mesh (&p4est2, &ymir_mesh2, &press_elem2, mpicomm,
                      &domain_options, &discr_options2, &slabs_options);

    buoyancy_type = 7;
    slabs_options.buoyancy_type = (slabs_buoyancy_type_t) buoyancy_type;
    slabs_setup_stokes (&stokes_problem2, ymir_mesh2, press_elem2,
                        &domain_options, &temp_options, &visc_options,
                        &slabs_options, vtk_write_input2_path);

    /* initialize solution vector */
    sol_vel_press2 = rhea_velocity_pressure_new (ymir_mesh2, press_elem2);

   /* run solver */
    slabs_run_solver (sol_vel_press2, ymir_mesh2, press_elem2, stokes_problem2,
                      solver_iter_max, solver_rel_tol);

    {
      ymir_vec_t         *velocity2 = rhea_velocity_new (ymir_mesh2);
      ymir_vec_t         *pressure2 = rhea_pressure_new (ymir_mesh2, press_elem2);
      ymir_vec_t         *viscosity2 = rhea_viscosity_new (ymir_mesh2);

      ymir_stokes_vec_get_components (sol_vel_press2, velocity2, pressure2,
                                      press_elem2);
      rhea_stokes_problem_copy_viscosity (viscosity2, stokes_problem2);

      rhea_vtk_write_solution (vtk_write_solution2_path, velocity2, pressure2,
                               viscosity2);

      rhea_pressure_destroy (pressure2);
      rhea_viscosity_destroy (viscosity2);
      rhea_velocity_destroy (velocity2);
    }

  /* compute surface normal stress sigma */
    surf_normal_stress2 = ymir_face_cvec_new (ymir_mesh2,
                                         RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);

    slabs_physics_compute_normal_boundary_stress (
                   surf_normal_stress2, sol_vel_press2,
                   rhea_stokes_problem_get_rhs_vel (stokes_problem2),
                   rhea_stokes_problem_get_stokes_op (stokes_problem2));

#if 0
/* test a weird case */
    {

    buoyancy_type = 6;
    slabs_options.buoyancy_type = (slabs_buoyancy_type_t) buoyancy_type;
    slabs_setup_stokes (&stokes_problem, ymir_mesh, press_elem,
                        &domain_options, &temp_options, &visc_options,
                        &slabs_options, vtk_write_input3_path);

    /* initialize solution vector */
    sol_vel_press3 = rhea_velocity_pressure_new (ymir_mesh, press_elem);

   /* run solver */
    slabs_run_solver (sol_vel_press3, ymir_mesh, press_elem, stokes_problem,
                      solver_iter_max, solver_rel_tol);

    {
      ymir_vec_t         *velocity3 = rhea_velocity_new (ymir_mesh);
      ymir_vec_t         *pressure3 = rhea_pressure_new (ymir_mesh, press_elem);
      ymir_vec_t         *viscosity3 = rhea_viscosity_new (ymir_mesh);

      ymir_stokes_vec_get_components (sol_vel_press3, velocity3, pressure3,
                                      press_elem);
      rhea_stokes_problem_copy_viscosity (viscosity3, stokes_problem);

      rhea_vtk_write_solution (vtk_write_solution3_path, velocity3, pressure3,
                               viscosity3);

      rhea_pressure_destroy (pressure3);
      rhea_viscosity_destroy (viscosity3);
      rhea_velocity_destroy (velocity3);
    }

  /* compute surface normal stress sigma */
    surf_normal_stress3 = ymir_face_cvec_new (ymir_mesh,
                                         RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);

    slabs_physics_compute_normal_boundary_stress (
                   surf_normal_stress3, sol_vel_press3,
                   rhea_stokes_problem_get_rhs_vel (stokes_problem),
                   rhea_stokes_problem_get_stokes_op (stokes_problem));
    }
#endif

    {
      char            path[BUFSIZ];

      snprintf (path, BUFSIZ, "%s", vtk_write_freesurface_path);
      ymir_vtk_write (ymir_mesh, path,
                      surf_normal_stress, "surf_normal_stress",
                      surf_normal_stress2, "surf_normal_stress with f-surface",
//                      surf_normal_stress3, "surf_normal_stress with air",
                      NULL);
    }

    /* destroy */
    sc_dmatrix_destroy (elem);
    free (tX);
    free (tY);
    free (tZ);
    ymir_vec_destroy (surf_normal_stress);
    ymir_vec_destroy (surf_normal_stress2);
//    ymir_vec_destroy (surf_normal_stress3);
    rhea_velocity_pressure_destroy (sol_vel_press2);
//    rhea_velocity_pressure_destroy (sol_vel_press3);
    slabs_surf_options.topo_profile = NULL;

    slabs_setup_clear_all (stokes_problem2, p4est2, ymir_mesh2, press_elem2,
                           &discr_options2);

    RHEA_GLOBAL_PRODUCTIONF ("In %s: Done vtk_write_freesurface\n", this_fn_name);
  }

  /* destroy */
  rhea_velocity_pressure_destroy (sol_vel_press);

   /*
    * Finalize
    */

  /* destroy Stokes problem and mesh */
  slabs_setup_clear_all (stokes_problem, p4est, ymir_mesh, press_elem,
                           &discr_options);

  /* destroy options */
  ymir_options_global_destroy ();

  /* print that this function is ending */
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

  /* finalize rhea */
  rhea_finalize ();

  return 0;
}
