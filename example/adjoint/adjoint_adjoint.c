#include <adjoint_adjoint.h>

/*********************************************************************************************
 * Sets up a linear Stokes problem.
 *********************************************************************************************/
/*Basic Stokes problem */
void
adjoint_stokes_new (rhea_stokes_problem_t **stokes_problem,
                    ymir_mesh_t **ymir_mesh,
                    ymir_pressure_elem_t **press_elem,
                    rhea_domain_options_t *domain_options,
                    rhea_temperature_options_t *temp_options,
                    rhea_weakzone_options_t *weak_options,
                    rhea_viscosity_options_t *visc_options,
                    subd_options_t *subd_options)
{
  const char         *this_fn_name = "adjoint_stokes_new";
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (*ymir_mesh);
  ymir_vec_t         *temperature;
  void               *solver_options = NULL;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* set up data */
  rhea_weakzone_data_create (weak_options, mpicomm);

  /* compute temperature */
  temperature = rhea_temperature_new (*ymir_mesh);
  subd_compute_temperature (temperature, temp_options, subd_options);

  subd_set_velocity_dirichlet_bc (domain_options, subd_options);

  /* create Stokes problem */
  *stokes_problem = rhea_stokes_problem_new (
      *ymir_mesh, *press_elem, temperature, domain_options, temp_options,
      weak_options, visc_options, solver_options);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

void
adjoint_stokes_update_init (rhea_stokes_problem_t **stokes_problem,
                          p4est_t     *p4est,
                          ymir_mesh_t **ymir_mesh,
                          ymir_pressure_elem_t **press_elem,
                          rhea_discretization_options_t *discr_options,
                          rhea_temperature_options_t *temp_options,
                          subd_options_t *subd_options,
                          const char *vtk_write_input_path)

{
  const char *this_fn_name = "adjoint_stokes_update";
  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /*call back to recompute rhs_vel,
   * rhs_velbc_nonzero_dir, rhs_velbc_nonzero_neumann, weakzone*/

  subd_options->rhs_type = SUBD_RHS_DENSITY;
  subd_compute_rhs_vel (*stokes_problem, subd_options);

//  subd_options->velbc_options->velbc_nonzero_neu = SUBD_VELBC_NEU_ZERO;
//  subd_compute_rhs_velbc_neumann (*stokes_problem, subd_options);

  subd_compute_weakzone (*stokes_problem, subd_options);
  /* reset function of viscosity */
  subd_viscosity_set_function (*stokes_problem, subd_options);

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
    subd_write_input_basic (*stokes_problem, temp_options,
                            vtk_write_input_path);
  }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

void
adjoint_run_solver (ymir_vec_t *usol,
                    ymir_vec_t *sol_vel_press,
                    ymir_mesh_t *ymir_mesh,
                    ymir_pressure_elem_t *press_elem,
                    rhea_stokes_problem_t *stokes_problem,
                    const int nonzero_initial_guess,
                    const int iter_max, const double rel_tol)
{
  ymir_vec_set_zero (sol_vel_press);
  subd_run_solver (sol_vel_press, ymir_mesh, press_elem, stokes_problem,
                   nonzero_initial_guess, iter_max, rel_tol);

  ymir_upvec_get_u (sol_vel_press, usol, YMIR_COPY);
}

/********************************************************************************************
 * Newton callback functions
 * ******************************************************************************************/

/* initialize solution m0 and solve for usol(m0) */
void
adjoint_data_init (ymir_vec_t *solution,
                         void *data)
{
  const char    *this_fn_name = "adjoint_data_init";

  adjoint_problem_t *adjoint = data;

  rhea_stokes_problem_t *stokes_problem = adjoint->stokes_problem;
  p4est_t   *p4est = adjoint->p4est;
  ymir_mesh_t   *ymir_mesh = adjoint->ymir_mesh;
  ymir_pressure_elem_t *press_elem = adjoint->press_elem;
  rhea_discretization_options_t *discr_options = adjoint->discr_options;
  rhea_temperature_options_t  *temp_options = adjoint->temp_options;
  subd_options_t  *subd_options = adjoint->subd_options;
  const char *vtk_write_input_path = adjoint->vtk_write_input_path;
  ymir_vec_t *usol = adjoint->usol;

  int solver_nonzero_initial_guess = 0;
  ymir_vec_t *sol_vel_press = adjoint->sol_vel_press;
  int solver_iter_max = adjoint->solver_iter_max;
  double solver_rel_tol = adjoint->solver_rel_tol;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  solution->meshfree->e[0][0] = subd_options->visc_options->visc_lith; /*TODO a function: set visc parameters*/
  ymir_vec_copy (solution, adjoint->msol);
RHEA_GLOBAL_PRODUCTIONF ("msol=%f\n", adjoint->msol->meshfree->e[0][0]);

  adjoint_stokes_update_init (&stokes_problem, p4est, &ymir_mesh, &press_elem,
                              discr_options, temp_options, subd_options,
                              vtk_write_input_path);

  /* solve Stokes problem*/
  adjoint_run_solver (usol, sol_vel_press, ymir_mesh, press_elem, stokes_problem,
                      solver_nonzero_initial_guess, solver_iter_max, solver_rel_tol);

  {
    char  path[BUFSIZ];
    char  *vtk_write_init_path = "/scratch/02600/xiliu/rhea/Adjoint/Test_newton/check_gradient/init";
      snprintf (path, BUFSIZ, "%s", vtk_write_init_path);
      ymir_vtk_write (ymir_mesh, path,
                      usol, "usol",
                      NULL);
  }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

double
adjoint_evaluate_objective (ymir_vec_t *solution, void *data)
{
  const char    *this_fn_name = "adjoint_evaluate_objective";
  adjoint_problem_t *adjoint_problem = data;
  double        obj_val;
  ymir_vec_t    *usol = adjoint_problem->usol;
  ymir_mesh_t   *ymir_mesh = adjoint_problem->ymir_mesh;
  ymir_vec_t    *usol_surf = ymir_face_cvec_new (ymir_mesh,
                          RHEA_DOMAIN_BOUNDARY_FACE_TOP, 3);
  ymir_vec_t    *usol_surfZ = ymir_face_cvec_new (ymir_mesh,
                          RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
  sc_dmatrix_t  *usol_surfZ_mat = usol_surfZ->cvec;
  ymir_cvec_t  *mass = ymir_vec_template (usol_surfZ);

  ymir_interp_vec (usol, usol_surf);
  ymir_cvec_get_comp (usol_surf, usol_surfZ_mat, 2, YMIR_COPY);

  ymir_mass_apply (usol_surfZ, mass);
  obj_val = ymir_vec_innerprod (usol_surfZ, mass);
  obj_val *= 0.5;

  {
    char  path[BUFSIZ];
    char  *vtk_write_object_path = "/scratch/02600/xiliu/rhea/Adjoint/Test_newton/check_gradient/object";
      snprintf (path, BUFSIZ, "%s", vtk_write_object_path);
      ymir_vtk_write (ymir_mesh, path,
                      usol, "usol",
                      NULL);
      RHEA_GLOBAL_PRODUCTIONF ("check_gradient: compute_objective %d\n", 0);

  }


  ymir_vec_destroy (usol_surf);
  ymir_vec_destroy (usol_surfZ);
  ymir_vec_destroy (mass);

  return (obj_val);
}

static void
adjoint_strainrate_dotprod (ymir_vec_t *dotprod,
                    ymir_vec_t *usol, ymir_vec_t *vsol)
{
  ymir_mesh_t *ymir_mesh = usol->mesh;
  ymir_vec_t *edotu = ymir_dvec_new (ymir_mesh, 6, YMIR_GAUSS_NODE);
  ymir_vec_t *edotv = ymir_dvec_new (ymir_mesh, 6, YMIR_GAUSS_NODE);


  ymir_velocity_strain_rate (usol, edotu, 0);
  ymir_velocity_strain_rate (vsol, edotv, 0);

  ymir_velocity_symtens_dotprod (edotu, edotv, dotprod);

  ymir_vec_destroy (edotu);
  ymir_vec_destroy (edotv);
}

void
adjoint_neg_gradient_from_velo (ymir_vec_t *neg_gradient,
                              ymir_vec_t *usol, ymir_vec_t *vsol,
                              subd_options_t *subd_opt)
{
  ymir_mesh_t *ymir_mesh = usol->mesh;
  ymir_vec_t *dotprod = ymir_dvec_new(ymir_mesh, 1, YMIR_GAUSS_NODE);
  ymir_dvec_t *vec_e = ymir_vec_template (dotprod);
  ymir_dvec_t *vec_mass = ymir_vec_template (dotprod);
  ymir_dvec_t *dotprod_mass = ymir_vec_template (dotprod);
  ymir_dvec_t *visc_grad = ymir_vec_template (dotprod);
  ymir_dvec_t *stencil_visc = ymir_vec_template (dotprod);

  /*TODO multiply with d_eta/d_m for different m */
  ymir_vec_set_value (visc_grad, 1.0);
  adjoint_stencil_visc (stencil_visc, subd_opt);
//  ymir_vec_set_value (stencil_visc, 1.0);
  ymir_vec_multiply_in (stencil_visc, visc_grad);
  ymir_vec_scale (2.0, visc_grad);

  adjoint_strainrate_dotprod (dotprod, usol, vsol);
  ymir_mass_apply (dotprod, dotprod_mass);
  neg_gradient->meshfree->e[0][0] = ymir_vec_innerprod (visc_grad, dotprod_mass);

//  ymir_vec_set_value (vec_e, 1.0);
//  ymir_mass_apply (vec_e, vec_mass);
//  ymir_vec_multiply_in (stencil_visc, dotprod);
//  neg_gradient->meshfree->e[0][0] = ymir_vec_innerprod (dotprod, vec_mass);

  ymir_vec_scale (-1.0, neg_gradient);

  ymir_vec_destroy (dotprod);
  ymir_vec_destroy (dotprod_mass);

  ymir_vec_destroy (vec_e);
  ymir_vec_destroy (vec_mass);
  ymir_vec_destroy (stencil_visc);
  ymir_vec_destroy (visc_grad);
}

void
adjoint_compute_negative_gradient (ymir_vec_t *neg_gradient,
                                   ymir_vec_t *solution, void *data)
{
  const char      *this_fn_name = "adjoint_compute_negative_gradient";
  adjoint_problem_t *adjoint_problem = data;
  ymir_vec_t  *usol = adjoint_problem->usol;
  ymir_vec_t  *vsol = adjoint_problem->vsol;
  ymir_vec_t  *sol_vel_press = adjoint_problem->sol_vel_press;

  ymir_pressure_elem_t *press_elem = adjoint_problem->press_elem;
  ymir_mesh_t   *ymir_mesh = adjoint_problem->ymir_mesh;
  subd_options_t *subd_options = adjoint_problem->subd_options;
  rhea_stokes_problem_t *stokes_problem = adjoint_problem->stokes_problem;
  ymir_stokes_op_t    *stokes_op;
  ymir_vec_t **rhs_vel_neumann;
  ymir_vec_t *rhs_vel_press;

  int solver_iter_max = adjoint_problem->solver_iter_max;
  double solver_rel_tol = adjoint_problem->solver_rel_tol;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);
  int solver_nonzero_initial_guess = 0;
  ymir_vec_set_zero (sol_vel_press);

//  subd_options->rhs_type = SUBD_RHS_ADJOINT_NONE;

  subd_options->velbc_options->velbc_nonzero_neu = SUBD_VELBC_NEU_VEC;
  subd_options->data = usol;
  subd_compute_rhs_velbc_neumann (stokes_problem, subd_options);

  rhs_vel_press = rhea_stokes_problem_get_rhs_vel_press (stokes_problem);
  rhs_vel_neumann = rhea_stokes_problem_get_rhs_vel_neumann (stokes_problem);
  stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);

  ymir_stokes_pc_construct_rhs (
      rhs_vel_press, //out and set to stokes
      NULL, //rhs_vel: volume forcing is zero
      rhs_vel_neumann, // from forward Stokes
      NULL, //rhs_vel_nonzero_dirichlet is zero
      1, //incompressible
      stokes_op,
      0);

  /* solve Stokes problem*/
  adjoint_run_solver (vsol, sol_vel_press, ymir_mesh, press_elem, stokes_problem,
                   solver_nonzero_initial_guess, solver_iter_max, solver_rel_tol);

  {
    char  path[BUFSIZ];
    char  *vtk_write_gradient_path = "/scratch/02600/xiliu/rhea/Adjoint/Test_newton/check_gradient/gradient";
      snprintf (path, BUFSIZ, "%s", vtk_write_gradient_path);
      ymir_vtk_write (ymir_mesh, path,
                      usol, "usol",
                      vsol, "vsol",
                      NULL);
RHEA_GLOBAL_PRODUCTIONF ("check_gradient: compute_gradient %d\n", 0);
  }

  adjoint_neg_gradient_from_velo (neg_gradient, usol, vsol, subd_options);

  subd_options->data = NULL;

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

double
adjoint_compute_gradient_norm (ymir_vec_t *neg_gradient,
                               void *data, double *norm_comp)
{
  double            norm;
  adjoint_problem_t *adjoint_problem = data;

  norm = ymir_vec_norm (neg_gradient);

  return (norm);
}

/*modify from rhea_stokes_problem_nonlinear_update_operator_fn */
void
adjoint_update_operator_fn (ymir_vec_t *solution, void *data)
{
  const char    *this_fn_name = "adjoint_update_operator";
  adjoint_problem_t *adjoint_problem = data;
  rhea_stokes_problem_t *stokes_problem = adjoint_problem->stokes_problem;
  p4est_t   *p4est = adjoint_problem->p4est;
  ymir_mesh_t   *ymir_mesh = adjoint_problem->ymir_mesh;
  ymir_pressure_elem_t *press_elem = adjoint_problem->press_elem;
  subd_options_t  *subd_options = adjoint_problem->subd_options;
  ymir_stress_op_t *stress_op;
  ymir_stokes_op_t *stokes_op;
  ymir_vec_t      *coeff;
  ymir_vec_t      *msol = adjoint_problem->msol;
  ymir_vec_t      *usol = adjoint_problem->usol;
  ymir_vec_t *rhs_vel;
  ymir_vec_t *rhs_vel_press;

  ymir_vec_t *sol_vel_press = adjoint_problem->sol_vel_press;
  int solver_nonzero_initial_guess = 0;
  int solver_iter_max = adjoint_problem->solver_iter_max;
  double solver_rel_tol = adjoint_problem->solver_rel_tol;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

//  if (msol->meshfree->e[0][0] == solution->meshfree->e[0][0]) {
//    return;
//  }
//  ymir_vec_copy (solution, msol);

  subd_options->visc_options->visc_lith = solution->meshfree->e[0][0]; /*TODO a function: get visc para*/
  RHEA_GLOBAL_PRODUCTIONF ("visc_lith=%12f\n", subd_options->visc_options->visc_lith);
  subd_viscosity_set_function (stokes_problem, subd_options);

  rhea_stokes_problem_compute_coefficient (stokes_problem, 0 /* !init */);
  coeff = rhea_stokes_problem_get_coeff (stokes_problem);
  stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
  stress_op = stokes_op->stress_op;
  ymir_stress_op_set_coeff_scal (stress_op, coeff);

  subd_options->rhs_type = SUBD_RHS_DENSITY;
  subd_compute_rhs_vel (stokes_problem, subd_options);

//  subd_compute_rhs_velbc_dirichlet (stokes_problem, subd_options);
//  subd_options->velbc_options->velbc_nonzero_neu = SUBD_VELBC_NEU_ZERO;
//  subd_compute_rhs_velbc_neumann (stokes_problem, subd_options);

  rhs_vel_press = rhea_stokes_problem_get_rhs_vel_press (stokes_problem);
  rhs_vel = rhea_stokes_problem_get_rhs_vel (stokes_problem);
  ymir_stokes_pc_construct_rhs (
      rhs_vel_press, //out and set to stokes
      rhs_vel, //rhs_vel: volume forcing is zero
      NULL, // from forward Stokes
      NULL, //rhs_vel_nonzero_dirichlet is zero
      1, //incompressible
      stokes_op,
      0);

  /* solve Stokes problem*/
  adjoint_run_solver (usol, sol_vel_press, ymir_mesh, press_elem, stokes_problem,
                     solver_nonzero_initial_guess, solver_iter_max, solver_rel_tol);
}

int
adjoint_solve_hessian_system_fn (ymir_vec_t *step, ymir_vec_t *neg_gradient,
                                const int lin_iter_max, const double lin_res_norm_rtol,
                                const int nonzero_initial_guess, void *data,
                                int *lin_iter_count)
{
  adjoint_problem_t *adjoint_problem = data;
  rhea_stokes_problem_t *stokes_problem = adjoint_problem->stokes_problem;
  p4est_t   *p4est = adjoint_problem->p4est;
  ymir_mesh_t   *ymir_mesh = adjoint_problem->ymir_mesh;
  ymir_pressure_elem_t *press_elem = adjoint_problem->press_elem;
  subd_options_t  *subd_options = adjoint_problem->subd_options;
  ymir_vec_t    *viscosity = rhea_viscosity_new (ymir_mesh);

  ymir_vec_t *sol_vel_press = adjoint_problem->sol_vel_press;
  ymir_vec_t *Husol = adjoint_problem->Husol;
  ymir_vec_t *Hvsol = adjoint_problem->Hvsol;
  ymir_vec_t *usol = adjoint_problem->usol;

  int solver_nonzero_initial_guess = 0;
  int solver_iter_max = adjoint_problem->solver_iter_max;
  double solver_rel_tol = adjoint_problem->solver_rel_tol;

  /*solve for Hessian forward*/
  {
    adjoint_set_rhs_hessian_forward (adjoint_problem);

//    subd_options->rhs_type = SUBD_RHS_ADJOINT_NONE;
//    subd_options->velbc_options->velbc_nonzero_neu = SUBD_VELBC_NEU_ZERO;

    /* solve Stokes problem*/
    adjoint_run_solver (Husol, sol_vel_press, ymir_mesh, press_elem, stokes_problem,
                     solver_nonzero_initial_guess, solver_iter_max, solver_rel_tol);
  }

  {
    ymir_vec_t    *rhs_vel_press;
    ymir_vec_t    **rhs_vel_neumann;
    ymir_stokes_op_t  *stokes_op;

//    subd_options->rhs_type = SUBD_RHS_ADJOINT_NONE;
//    subd_compute_rhs_vel (stokes_problem, subd_options);

    subd_options->velbc_options->velbc_nonzero_neu = SUBD_VELBC_NEU_VEC;
    subd_options->data = Husol;
    subd_compute_rhs_velbc_neumann (stokes_problem, subd_options);

    rhs_vel_press = rhea_stokes_problem_get_rhs_vel_press (stokes_problem);
    rhs_vel_neumann = rhea_stokes_problem_get_rhs_vel_neumann (stokes_problem);
    stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);

    ymir_stokes_pc_construct_rhs (
      rhs_vel_press, //out and set to stokes
      NULL, //rhs_vel: volume forcing is zero
      rhs_vel_neumann, // from stokes
      NULL, //rhs_vel_nonzero_dirichlet is zero
      1, //incompressible
      stokes_op,
      0);

    /* solve Stokes problem*/
    adjoint_run_solver (Hvsol, sol_vel_press, ymir_mesh, press_elem, stokes_problem,
                       solver_nonzero_initial_guess, solver_iter_max, solver_rel_tol);
  }

  rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);
  adjoint_hessian_from_velo (step, usol, Hvsol, viscosity, subd_options);
  ymir_vec_reciprocal (step);
  ymir_vec_multiply_in (neg_gradient, step);

  rhea_viscosity_destroy (viscosity);

  /* return iteraton count and "stopping" reason */
  if (lin_iter_count != NULL) {
    *lin_iter_count = 1;
  }

  return 1;
}

void
adjoint_update_hessian_fn (ymir_vec_t *solution, ymir_vec_t *step,
                          const double step_length, void *data)
{
  adjoint_problem_t   *adjoint_problem = (adjoint_problem_t *) data;

  ymir_vec_copy (solution, adjoint_problem->msol);
}


/********************************************************************************************
 * create and destroy adjoint_problem
 ********************************************************************************************/
adjoint_problem_t *
adjoint_problem_new (rhea_stokes_problem_t *stokes_problem,
                       p4est_t *p4est, ymir_mesh_t *ymir_mesh,
                       ymir_pressure_elem_t *press_elem,
                       rhea_discretization_options_t *discr_options,
                       rhea_temperature_options_t *temp_options,
                       subd_options_t *subd_options,
                       const char *vtk_write_input_path,
                       int       solver_iter_max,
                       double    solver_rel_tol)
{
  adjoint_problem_t *adjoint = RHEA_ALLOC (adjoint_problem_t, 1);

  adjoint->stokes_problem = stokes_problem;
  adjoint->p4est = p4est;
  adjoint->ymir_mesh = ymir_mesh;
  adjoint->press_elem = press_elem;
  adjoint->discr_options = discr_options;
  adjoint->temp_options = temp_options;
  adjoint->subd_options = subd_options;
  adjoint->vtk_write_input_path = vtk_write_input_path;

  adjoint->solver_iter_max = solver_iter_max;
  adjoint->solver_rel_tol = solver_rel_tol;

  adjoint->msol = ymir_vec_new_meshfree (1);
  adjoint->usol = rhea_velocity_new (ymir_mesh);
  adjoint->vsol = rhea_velocity_new (ymir_mesh);
  adjoint->Husol = rhea_velocity_new (ymir_mesh);
  adjoint->Hvsol = rhea_velocity_new (ymir_mesh);
  adjoint->sol_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);

  return (adjoint);
}

void
adjoint_problem_destroy (rhea_newton_problem_t *newton_problem)
{
  adjoint_problem_t *adjoint_problem = rhea_newton_problem_get_data (newton_problem);

  RHEA_ASSERT (adjoint_problem->msol != NULL);
  RHEA_ASSERT (adjoint_problem->usol != NULL);
  RHEA_ASSERT (adjoint_problem->vsol != NULL);
  RHEA_ASSERT (adjoint_problem->Husol != NULL);
  RHEA_ASSERT (adjoint_problem->Hvsol != NULL);
  RHEA_ASSERT (adjoint_problem->sol_vel_press != NULL);

  ymir_vec_destroy (adjoint_problem->msol);
  rhea_velocity_destroy (adjoint_problem->usol);
  rhea_velocity_destroy (adjoint_problem->vsol);
  rhea_velocity_destroy (adjoint_problem->Husol);
  rhea_velocity_destroy (adjoint_problem->Hvsol);
  rhea_velocity_pressure_destroy (adjoint_problem->sol_vel_press);
  RHEA_FREE (adjoint_problem);

  /* destroy Newton problem */
  ymir_vec_destroy (rhea_newton_problem_get_neg_gradient_vec (newton_problem));
  ymir_vec_destroy (rhea_newton_problem_get_step_vec (newton_problem));
  rhea_newton_problem_destroy (newton_problem);
}


void
adjoint_setup_newton (rhea_newton_problem_t **newton_problem,
                      adjoint_problem_t *adjoint_problem)
{
  char      *this_fn_name = "adjoint_setup_newton";

 RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);


  /*create newton problem */
  {
    ymir_vec_t    *neg_gradient_vec = ymir_vec_new_meshfree (1);
//    ymir_vec_t    *step_vec = ymir_vec_template (neg_gradient_vec);
    ymir_vec_t    *step_vec = ymir_vec_new_meshfree (1);

    ymir_vec_set_zero (neg_gradient_vec);
    ymir_vec_set_zero (step_vec);

    *newton_problem = rhea_newton_problem_new (
        adjoint_compute_negative_gradient,
        adjoint_solve_hessian_system_fn);

    rhea_newton_problem_set_vectors (
        neg_gradient_vec, step_vec, *newton_problem);

    rhea_newton_problem_set_data_fn (
        adjoint_problem, adjoint_data_init, NULL /* no clear fnc. */,
        *newton_problem);

    rhea_newton_problem_set_conv_criterion_fn (
//        RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM,
        RHEA_NEWTON_CONV_CRITERION_OBJECTIVE,
        adjoint_evaluate_objective,
        adjoint_compute_gradient_norm,
        0 /* no multi-component norms */, *newton_problem);

//    rhea_newton_problem_set_apply_hessian_fn (
//        adjoint_apply_hessian_fn, *newton_problem);

    rhea_newton_problem_set_update_fn (
        adjoint_update_operator_fn,
        adjoint_update_hessian_fn,
        NULL, /* adjoint_modify_hessian_system_fn, is not necessary */
        *newton_problem);
  }

#ifdef RHEA_ENABLE_DEBUG
  rhea_newton_problem_set_checks (1 /* grad */, 0 /* Hessian */, *newton_problem);
#endif

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**************************************************************
 * Adjoint related
***************************************************************/
static void
subd_adjoint_stencil_visc_custom_layers_um_elem (double *_sc_restrict visc_elem,
                                            const double *_sc_restrict x,
                                            const double *_sc_restrict y,
                                            const double *_sc_restrict z,
                                            const int n_nodes,
                                            subd_options_t *subd_options)
{
  int                 nodeid;
  double              z_mid = subd_options->visc_options->z_lith;
  double              stencil_um = subd_options->adjoint_options->stencil_options->value;
  double              stencil_lm = .0;
  int                 m;

  /* compute viscosity in this element */
    for (nodeid = 0; nodeid < n_nodes; nodeid++) {
      if (z[nodeid] > z_mid)  {
        visc_elem[nodeid] = stencil_um;
      }
      else {
        visc_elem[nodeid] = stencil_lm;
      }
    }
    /* check viscosity for `nan`, `inf`, and positivity */
    RHEA_ASSERT (isfinite (visc_elem[nodeid]));
}

void
adjoint_stencil_visc  (ymir_vec_t *stencil_visc,
                         subd_options_t *subd_options)
{
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (stencil_visc);
  const ymir_locidx_t  n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);
  mangll_t           *mangll = mesh->ma;
  const int           N = ymir_n (mangll->N);
  char                *this_fn_name = "subd_adjoint_visc_stencil";

  sc_dmatrix_t       *stencil_el_mat;
  double             *x, *y, *z, *tmp_el,*stencil_el_data;
  ymir_locidx_t       elid;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);
  /* create work variables */
  stencil_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  stencil_el_data = stencil_el_mat->e[0];
  x = RHEA_ALLOC (double, n_nodes_per_el);
  y = RHEA_ALLOC (double, n_nodes_per_el);
  z = RHEA_ALLOC (double, n_nodes_per_el);
  tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) {
    /* get coordinates of this element at Gauss nodes */
    ymir_mesh_get_elem_coord_gauss (x, y, z, elid, mesh, tmp_el);

    /* compute stencilosity stencil*/ //TODO: later on add visc_rhea, and other viscosity
    if (subd_options->adjoint_options->stencil_options->type == SUBD_ADJOINT_STENCIL_VISC_CUSTOM) {
      switch (subd_options->adjoint_options->stencil_options->custom_type) {
        case SUBD_ADJOINT_STENCIL_VISC_CUSTOM_LAYERS_UM:
          subd_adjoint_stencil_visc_custom_layers_um_elem (stencil_el_data, x, y, z, n_nodes_per_el,
                                           subd_options);
        break;

        default:
        RHEA_ABORT_NOT_REACHED ();
      }
    }
    /* set viscosity of this element */
    rhea_viscosity_set_elem_gauss (stencil_visc, stencil_el_mat, elid);
  }

  /* destroy */
  sc_dmatrix_destroy (stencil_el_mat);
  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);
}

void
adjoint_set_rhs_hessian_forward (adjoint_problem_t *adjoint_problem)
{
  subd_options_t        *subd_opt = adjoint_problem->subd_options;
  rhea_domain_options_t    *domain_opt = subd_opt->domain_options;
  ymir_mesh_t           *ymir_mesh = adjoint_problem->ymir_mesh;
  ymir_pressure_elem_t  *press_elem = adjoint_problem->press_elem;
  rhea_stokes_problem_t *stokes_problem = adjoint_problem->stokes_problem;
  ymir_vec_t            *usol = adjoint_problem->usol;

  ymir_vec_t            *rhs_vel_press =
                              rhea_velocity_pressure_new (ymir_mesh, press_elem);
  ymir_vec_t            *viscosity = rhea_viscosity_new (ymir_mesh);
  ymir_vec_t            *stencil_visc = rhea_viscosity_new (ymir_mesh);
  ymir_stokes_op_t      *stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
  ymir_stress_op_t      *stress_op;
  ymir_vel_dir_t        *vel_dir = rhea_stokes_problem_get_vel_dir (stokes_problem);

  ymir_vec_t            *rhs_vel = rhea_velocity_new (ymir_mesh);

  const char            *this_fn_name = "subd_adjoint_rhs_hessian_forward";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);
  adjoint_stencil_visc (stencil_visc, subd_opt);
  ymir_vec_multiply_in (stencil_visc, viscosity);
  ymir_vec_scale (2.0, viscosity); // coeff

//  stress_op = stokes_op->stress_op;
  stress_op = ymir_stress_op_new (viscosity, vel_dir, NULL, NULL, NULL);
  ymir_stress_pc_apply_stress_op (usol, rhs_vel, stress_op, 0, 0);
  ymir_vec_scale (-1.0, rhs_vel);

  ymir_vec_set_zero (rhs_vel_press);
  ymir_upvec_set_u (rhs_vel_press, rhs_vel, YMIR_SET);
  rhea_stokes_problem_set_rhs_vel_press (stokes_problem, rhs_vel_press);

  ymir_stress_op_destroy (stress_op);
  rhea_viscosity_destroy (stencil_visc);
  rhea_viscosity_destroy (viscosity);
  rhea_velocity_destroy (rhs_vel);
  rhea_velocity_pressure_destroy (rhs_vel_press);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

void
adjoint_hessian_from_velo (ymir_vec_t *step,
                        ymir_vec_t *usol, ymir_vec_t *Hvsol,
                        ymir_vec_t *viscosity,
                        subd_options_t *subd_opt)
{
  ymir_mesh_t *ymir_mesh = usol->mesh;
  ymir_vec_t *dotprod = ymir_dvec_new(ymir_mesh, 1, YMIR_GAUSS_NODE);
  ymir_dvec_t *dotprod_mass = ymir_vec_template (dotprod);
  ymir_dvec_t *vec_e = ymir_vec_template (dotprod);
  ymir_dvec_t *vec_mass = ymir_vec_template (dotprod);
  ymir_vec_t *stencil_visc = ymir_vec_template (viscosity);

  adjoint_stencil_visc (stencil_visc, subd_opt);
  ymir_vec_multiply_in (stencil_visc, viscosity);
  ymir_vec_scale (2.0, viscosity);

  adjoint_strainrate_dotprod (dotprod, usol, Hvsol);
  ymir_mass_apply (dotprod, dotprod_mass);

  step->meshfree->e[0][0] = ymir_vec_innerprod (viscosity, dotprod_mass);

  ymir_vec_destroy (dotprod);
  ymir_vec_destroy (dotprod_mass);

  ymir_vec_destroy (vec_e);
  ymir_vec_destroy (vec_mass);
  ymir_vec_destroy (stencil_visc);
}

double
subd_adjoint_gradient (ymir_vec_t *edot0, ymir_vec_t *edot1,
                            ymir_vec_t *temp, ymir_vec_t *visc)
{
  ymir_mesh_t *ymir_mesh = ymir_vec_get_mesh (edot0);
  ymir_dvec_t *dotprod = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
  ymir_dvec_t *integrant = ymir_vec_template (dotprod);
  ymir_dvec_t *vec_e = ymir_vec_template (dotprod);
  ymir_dvec_t *masse = ymir_vec_template (dotprod);

  const int           N  = ymir_n (ymir_mesh->cnodes->N);
  const int           Np = ymir_np (ymir_mesh->cnodes->N);
  const int           K  = ymir_mesh->cnodes->K;
  sc_dmatrix_t       *elemtemp = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemvisc = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemdotprod = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemintegrant = sc_dmatrix_new (1, Np);

  int     elid, gp;
  double integral = .0;

  ymir_velocity_symtens_dotprod (edot0, edot1, dotprod);

  for (elid = 0; elid < K; elid++) {
    ymir_cvec_get_elem_interp (temp, elemtemp, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_COPY);
    ymir_dvec_get_elem (visc, elemvisc, YMIR_STRIDE_NODE, elid,
                               YMIR_COPY);
    ymir_dvec_get_elem (dotprod, elemdotprod, YMIR_STRIDE_NODE, elid,
                               YMIR_COPY);
    ymir_dvec_get_elem (integrant, elemintegrant, YMIR_STRIDE_NODE, elid,
                               YMIR_WRITE);
    double     *_sc_restrict tempd = elemtemp->e[0];
    double     *_sc_restrict viscd = elemvisc->e[0];
    double     *_sc_restrict dotprodd = elemdotprod->e[0];
    double     *_sc_restrict integrantd = elemintegrant->e[0];


    for (gp = 0; gp < Np; gp++) {
      integrantd[gp] = 2.0 * viscd[gp] * dotprodd[gp];

    }
    ymir_dvec_set_elem (integrant, elemintegrant, YMIR_STRIDE_NODE, elid,
                        YMIR_SET);
  }
  ymir_vec_set_value (vec_e, 1.0);
  ymir_mass_apply (vec_e, masse);
  integral = ymir_vec_innerprod (integrant, masse);
//  integral = 2.0 * ymir_vec_innerprod (dotprod, masse);

  sc_dmatrix_destroy (elemtemp);
  sc_dmatrix_destroy (elemvisc);
  sc_dmatrix_destroy (elemdotprod);
  sc_dmatrix_destroy (elemintegrant);
  ymir_vec_destroy (dotprod);
  ymir_vec_destroy (integrant);
  ymir_vec_destroy (vec_e);
  ymir_vec_destroy (masse);

  return (integral);
}
