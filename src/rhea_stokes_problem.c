/*
 */

#include <rhea_stokes_problem.h>
#include <rhea_base.h>
#include <rhea_velocity.h>
#include <rhea_velocity_pressure.h>
#include <ymir_stokes_op.h>
#include <ymir_stokes_pc.h>

/******************************************************************************
 * General Stokes Problem
 *****************************************************************************/

/* enumerator for types of Stokes problems */
typedef enum
{
  RHEA_STOKES_PROBLEM_LINEAR,
  RHEA_STOKES_PROBLEM_NONLINEAR
}
rhea_stokes_problem_type_t;

/* Stokes problem */
struct rhea_stokes_problem
{
  rhea_stokes_problem_type_t  type;

  /* mesh (not owned) */
  ymir_mesh_t           *ymir_mesh;
  ymir_pressure_elem_t  *press_elem;

  /* data of the Stokes problem */
  ymir_vec_t         *coeff;
  ymir_vec_t         *rank1_tensor_scal;
  ymir_vec_t         *bounds_marker;
  ymir_vec_t         *yielding_marker;
  ymir_vel_dir_t     *vel_dir;
  ymir_vec_t         *rhs_vel;
  ymir_vec_t         *rhs_vel_press;

  /* Stokes operator and preconditioner */
  ymir_stokes_op_t   *stokes_op;
  ymir_stokes_pc_t   *stokes_pc;
};

/**
 * Creates and initializes a new structure of a Stokes problem.
 */
static rhea_stokes_problem_t *
rhea_stokes_problem_struct_new (const rhea_stokes_problem_type_t type,
                                ymir_mesh_t *ymir_mesh,
                                ymir_pressure_elem_t *press_elem)
{
  rhea_stokes_problem_t *stokes_problem = RHEA_ALLOC (rhea_stokes_problem_t, 1);

  stokes_problem->type = type;
  stokes_problem->ymir_mesh = ymir_mesh;
  stokes_problem->press_elem = press_elem;

  stokes_problem->coeff = NULL;
  stokes_problem->rank1_tensor_scal = NULL;
  stokes_problem->bounds_marker = NULL;
  stokes_problem->yielding_marker = NULL;
  stokes_problem->vel_dir = NULL;
  stokes_problem->rhs_vel = NULL;
  stokes_problem->rhs_vel_press = NULL;

  stokes_problem->stokes_op = NULL;
  stokes_problem->stokes_pc = NULL;

  return stokes_problem;
}

/**
 * Destroys a structure of a Stokes problem.
 */
static void
rhea_stokes_problem_struct_destroy (rhea_stokes_problem_t *stokes_problem)
{
  RHEA_FREE (stokes_problem);
}

rhea_stokes_problem_t *
rhea_stokes_problem_new (ymir_vec_t *temperature,
                         ymir_vec_t *weakzone,
                         ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         rhea_domain_options_t *domain_options,
                         rhea_temperature_options_t *temp_options,
                         rhea_viscosity_options_t *visc_options)
{
  if (RHEA_VISCOSITY_NONLINEAR == visc_options->type) {
    return rhea_stokes_problem_nonlinear_new (
        temperature, weakzone, ymir_mesh, press_elem,
        domain_options, temp_options, visc_options);
  }
  else {
    return rhea_stokes_problem_linear_new (
        temperature, weakzone, ymir_mesh, press_elem,
        domain_options, temp_options, visc_options);
  }
}

void
rhea_stokes_problem_destroy (rhea_stokes_problem_t *stokes_problem)
{
  if (RHEA_STOKES_PROBLEM_NONLINEAR == stokes_problem->type) {
    rhea_stokes_problem_nonlinear_destroy (stokes_problem);
  }
  else {
    rhea_stokes_problem_linear_destroy (stokes_problem);
  }
}

void
rhea_stokes_problem_setup_solver (rhea_stokes_problem_t *stokes_problem)
{
  if (RHEA_STOKES_PROBLEM_NONLINEAR == stokes_problem->type) {
    rhea_stokes_problem_nonlinear_setup_solver (stokes_problem);
  }
  else {
    rhea_stokes_problem_linear_setup_solver (stokes_problem);
  }
}

void
rhea_stokes_problem_solve (ymir_vec_t *sol_vel_press,
                           const int iter_max,
                           const double rel_tol,
                           rhea_stokes_problem_t *stokes_problem)
{
  if (RHEA_STOKES_PROBLEM_NONLINEAR == stokes_problem->type) {
    rhea_stokes_problem_nonlinear_solve (sol_vel_press, rel_tol, iter_max,
                                         stokes_problem);
  }
  else {
    rhea_stokes_problem_linear_solve (sol_vel_press, rel_tol, iter_max,
                                      stokes_problem);
  }
}

void
rhea_stokes_problem_get_viscosity (ymir_vec_t *viscosity,
                                   rhea_stokes_problem_t *stokes_problem)
{
  /* check input */
  RHEA_ASSERT (stokes_problem->coeff != NULL);

  /* copy Stokes coefficient and divide by 2 */
  ymir_vec_copy (stokes_problem->coeff, viscosity);
  ymir_vec_scale (0.5, viscosity);
}

/**
 * Creates new velocity boundary conditions.
 */
static void
rhea_stokes_problem_velocity_boundary_create (
                                        rhea_stokes_problem_t *stokes_problem,
                                        ymir_mesh_t *ymir_mesh,
                                        rhea_domain_options_t *domain_options)
{
  ymir_vec_t         *dirscal = NULL; //TODO seems to be best choice, isn't it?
  ymir_vel_dir_t     *vel_dir;

  stokes_problem->vel_dir = rhea_domain_create_velocity_dirichlet_bc (
      ymir_mesh, dirscal, domain_options);
}

/**
 * Destroys velocity boundary conditions (if they exist).
 */
static void
rhea_stokes_problem_velocity_boundary_clear (
                                        rhea_stokes_problem_t *stokes_problem)
{
  if (stokes_problem->vel_dir != NULL) {
    if (stokes_problem->vel_dir->scale != NULL) {
      ymir_vec_destroy (stokes_problem->vel_dir->scale);
    }
    ymir_vel_dir_destroy (stokes_problem->vel_dir);
  }
}

void
rhea_stokes_problem_velocity_boundary_set_zero (
                                        ymir_vec_t *velocity,
                                        rhea_stokes_problem_t *stokes_problem)
{
  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (velocity));

  /* set Dirichlet components to zero if Dirichlet BC's exist */
  if (stokes_problem->vel_dir != NULL) {
    ymir_vel_dir_separate (velocity, NULL, NULL, NULL, stokes_problem->vel_dir);
  }
}

/******************************************************************************
 * Linear Stokes Problem
 *****************************************************************************/

rhea_stokes_problem_t *
rhea_stokes_problem_linear_new (ymir_vec_t *temperature,
                                ymir_vec_t *weakzone,
                                ymir_mesh_t *ymir_mesh,
                                ymir_pressure_elem_t *press_elem,
                                rhea_domain_options_t *domain_options,
                                rhea_temperature_options_t *temp_options,
                                rhea_viscosity_options_t *visc_options)
{
  const char         *this_fn_name = "rhea_stokes_problem_linear_new";
  rhea_stokes_problem_t *stokes_problem_lin;
  ymir_vec_t         *coeff;
  ymir_vec_t         *rhs_vel;
  ymir_vec_t         *rhs_vel_press;
  ymir_stokes_op_t   *stokes_op;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (visc_options->type == RHEA_VISCOSITY_LINEAR);

  /* initialize Stokes problem structure */
  stokes_problem_lin = rhea_stokes_problem_struct_new (
      RHEA_STOKES_PROBLEM_LINEAR, ymir_mesh, press_elem);

  /* create velocity boundary conditions */
  rhea_stokes_problem_velocity_boundary_create (stokes_problem_lin,
                                                ymir_mesh, domain_options);

  /* create coefficient */
  {
    ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);

    /* compute viscosity */
    rhea_viscosity_compute (viscosity,
                            NULL /* nl. Stokes output */,
                            NULL /* nl. Stokes output */,
                            NULL /* nl. Stokes output */,
                            temperature, weakzone,
                            NULL /* nl. Stokes input */,
                            visc_options);

    /* set Stokes coefficient */
    coeff = viscosity;
    ymir_vec_scale (2.0, coeff);
  }

  /* create Stokes operator */
  stokes_op = ymir_stokes_op_new_ext (coeff, stokes_problem_lin->vel_dir,
                                      NULL /* Robin BC's */,
                                      NULL /* nl. Stokes input */,
                                      NULL /* nl. Stokes input */,
                                      press_elem, domain_options->center,
                                      domain_options->moment_of_inertia);

  /* create velocity right-hand side forcing */
  rhs_vel = rhea_velocity_new (ymir_mesh);
  rhea_temperature_compute_rhs_vel (rhs_vel, temperature, temp_options);

  /* construct right-hand side for incompressible Stokes system */
  rhs_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);
  ymir_stokes_op_construct_rhs_ext (rhs_vel /* Dirichlet forcing */,
                                    NULL /* Neumann forcing */,
                                    NULL /* Dirichlet lift */,
                                    rhs_vel_press,
                                    1 /* incompressible */,
                                    stokes_op);

  /* fill, and return the structure of the linear Stokes problem */
  stokes_problem_lin->coeff = coeff;
  stokes_problem_lin->rhs_vel = rhs_vel;
  stokes_problem_lin->rhs_vel_press = rhs_vel_press;
  stokes_problem_lin->stokes_op = stokes_op;

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

  return stokes_problem_lin;
}

void
rhea_stokes_problem_linear_destroy (rhea_stokes_problem_t *stokes_problem_lin)
{
  const char         *this_fn_name = "rhea_stokes_problem_linear_destroy";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);
  RHEA_ASSERT (stokes_problem_lin->stokes_op != NULL);
  RHEA_ASSERT (stokes_problem_lin->coeff != NULL);
  RHEA_ASSERT (stokes_problem_lin->rank1_tensor_scal == NULL);
  RHEA_ASSERT (stokes_problem_lin->bounds_marker == NULL);
  RHEA_ASSERT (stokes_problem_lin->yielding_marker == NULL);
  RHEA_ASSERT (stokes_problem_lin->rhs_vel != NULL);
  RHEA_ASSERT (stokes_problem_lin->rhs_vel_press != NULL);

  /* destroy Stokes operator and preconditioner */
  if (stokes_problem_lin->stokes_pc != NULL) {
    ymir_stokes_pc_destroy (stokes_problem_lin->stokes_pc);
  }
  ymir_stokes_op_destroy (stokes_problem_lin->stokes_op);

  /* destroy vectors */
  rhea_viscosity_destroy (stokes_problem_lin->coeff);
  rhea_velocity_destroy (stokes_problem_lin->rhs_vel);
  rhea_velocity_pressure_destroy (stokes_problem_lin->rhs_vel_press);

  /* destroy velocity boundary conditions */
  rhea_stokes_problem_velocity_boundary_clear (stokes_problem_lin);

  /* destroy problem structure */
  rhea_stokes_problem_struct_destroy (stokes_problem_lin);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

void
rhea_stokes_problem_linear_setup_solver (
                                    rhea_stokes_problem_t *stokes_problem_lin)
{
  const char         *this_fn_name = "rhea_stokes_problem_linear_setup_solver";
  ymir_stokes_op_t   *stokes_op = stokes_problem_lin->stokes_op;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);
  RHEA_ASSERT (stokes_problem_lin->stokes_op != NULL);

  /* build Stokes preconditioner */
  stokes_problem_lin->stokes_pc = ymir_stokes_pc_new (stokes_op);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

void
rhea_stokes_problem_linear_solve (ymir_vec_t *sol_vel_press,
                                  const int iter_max,
                                  const double rel_tol,
                                  rhea_stokes_problem_t *stokes_problem_lin)
{
  const char         *this_fn_name = "rhea_stokes_problem_linear_solve";
  const int           nonzero_initial_guess = 0;
  const double        abs_tol = 1.0e-100;
  const int           krylov_gmres_n_vecs = 100;
  const int           out_residual = !rhea_get_production_run ();

  double              norm_res_init, norm_res;
  int                 itn, stop_reason;

  RHEA_GLOBAL_INFOF ("Into %s (rel tol %.3e, max iter %i)\n",
                     this_fn_name, rel_tol, iter_max);

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);
  RHEA_ASSERT (stokes_problem_lin->stokes_pc != NULL);
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (sol_vel_press));

  /* compute residual with zero inital guess */
  if (out_residual) {
    norm_res_init = -1.0;
  //TODO
  //ymir_vec_set_zero (state->vel_press_vec);
  //norm_res_init = slabs_norm_compute_residual (
  //    NULL, NULL, NULL, NULL, lin_stokes->rhs_u_point,
  //    lin_stokes->stokes_op, lin_stokes->stokes_pc,
  //    solver_options->krylov_type, SL_NORM_VEC_L2, NULL);
  }

  /* run solver for Stokes */
  stop_reason = ymir_stokes_pc_solve (stokes_problem_lin->rhs_vel_press,
                                      sol_vel_press,
                                      stokes_problem_lin->stokes_pc,
                                      nonzero_initial_guess,
                                      rel_tol, abs_tol, iter_max,
                                      krylov_gmres_n_vecs /* unused */,
                                      &itn);

  /* compute residual at solution */
  if (out_residual) {
    double              res_reduction, conv_factor;

    norm_res = 1.0;
  //TODO
  //norm_res = slabs_norm_compute_residual (
  //    residual_up, NULL, NULL, state->vel_press_vec,
  //    lin_stokes->rhs_u_point,
  //    lin_stokes->stokes_op, lin_stokes->stokes_pc,
  //    solver_options->krylov_type, SL_NORM_VEC_L2, NULL);

    /* calculate convergence */
    res_reduction = norm_res/norm_res_init;
    conv_factor = exp (log (res_reduction) / ((double) itn));

    /* project out nullspaces of solution */
    //TODO
//  slabs_linear_stokes_problem_project_out_nullspace (
//      state->vel_press_vec, lin_stokes->stokes_op, physics_options, 0);

    RHEA_GLOBAL_INFOF ("Done %s (solver stopping reason %i, iterations %i, "
                       "residual reduction %.3e, convergence %.3e)\n",
                       this_fn_name, stop_reason, itn,
                       res_reduction, conv_factor);
  }
  else {
    /* project out nullspaces of solution */
    //TODO
//  slabs_linear_stokes_problem_project_out_nullspace (
//      state->vel_press_vec, lin_stokes->stokes_op, physics_options, 0);

    RHEA_GLOBAL_INFOF ("Done %s (solver stopping reason %i, num iter %i)\n",
                       this_fn_name, stop_reason, itn);
  }
}

/******************************************************************************
 * Nonlinear Stokes Problem
 *****************************************************************************/

rhea_stokes_problem_t *
rhea_stokes_problem_nonlinear_new (ymir_vec_t *temperature,
                                   ymir_vec_t *weakzone,
                                   ymir_mesh_t *ymir_mesh,
                                   ymir_pressure_elem_t *press_elem,
                                   rhea_domain_options_t *domain_options,
                                   rhea_temperature_options_t *temp_options,
                                   rhea_viscosity_options_t *visc_options)
{
  const char         *this_fn_name = "rhea_stokes_problem_nonlinear_new";
  rhea_stokes_problem_t *stokes_problem_nl;
  ymir_vec_t         *coeff;
  ymir_vec_t         *rank1_tensor_scal;
  ymir_vec_t         *bounds_marker;
  ymir_vec_t         *yielding_marker;
  ymir_vec_t         *rhs_vel;
  ymir_vec_t         *rhs_vel_press;
//ymir_stokes_op_t   *stokes_op; //TODO

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (visc_options->type == RHEA_VISCOSITY_NONLINEAR);

  /* initialize Stokes problem structure */
  stokes_problem_nl = rhea_stokes_problem_struct_new (
      RHEA_STOKES_PROBLEM_NONLINEAR, ymir_mesh, press_elem);

  /* create velocity boundary conditions */
  rhea_stokes_problem_velocity_boundary_create (stokes_problem_nl,
                                                ymir_mesh, domain_options);

  /* create coefficient vectors */
  coeff = rhea_viscosity_new (ymir_mesh);
  rank1_tensor_scal = rhea_viscosity_new (ymir_mesh);
  bounds_marker = rhea_viscosity_new (ymir_mesh);
  yielding_marker = rhea_viscosity_new (ymir_mesh);

  /* create Stokes operator */
//TODO
//stokes_op = ymir_stokes_op_new_ext (coeff, stokes_problem_nl->vel_dir,
//                                    NULL /* Robin BC's */,
//                                    NULL /* nl. Stokes input */,
//                                    NULL /* nl. Stokes input */,
//                                    press_elem, domain_options->center,
//                                    domain_options->moment_of_inertia);

  /* create velocity right-hand side forcing */
  rhs_vel = rhea_velocity_new (ymir_mesh);
  rhea_temperature_compute_rhs_vel (rhs_vel, temperature, temp_options);

  /* construct right-hand side for incompressible Stokes system */
  rhs_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);
//TODO
//ymir_stokes_op_construct_rhs_ext (rhs_vel /* Dirichlet forcing */,
//                                  NULL /* Neumann forcing */,
//                                  NULL /* Dirichlet lift */,
//                                  rhs_vel_press,
//                                  1 /* incompressible */,
//                                  stokes_op);

  /* fill, and return the structure of the nonlinear Stokes problem */
  stokes_problem_nl->coeff = coeff;
  stokes_problem_nl->rank1_tensor_scal = rank1_tensor_scal;
  stokes_problem_nl->bounds_marker = bounds_marker;
  stokes_problem_nl->yielding_marker = yielding_marker;
  stokes_problem_nl->rhs_vel = rhs_vel;
  stokes_problem_nl->rhs_vel_press = rhs_vel_press;
//stokes_problem_nl->stokes_op = stokes_op; //TODO

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

  return stokes_problem_nl;
}

void
rhea_stokes_problem_nonlinear_destroy (
                                    rhea_stokes_problem_t *stokes_problem_nl)
{
  const char         *this_fn_name = "rhea_stokes_problem_nonlinear_destroy";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
//RHEA_ASSERT (stokes_problem_nl->stokes_op != NULL); //TODO
  RHEA_ASSERT (stokes_problem_nl->coeff != NULL);
  RHEA_ASSERT (stokes_problem_nl->rank1_tensor_scal != NULL);
  RHEA_ASSERT (stokes_problem_nl->bounds_marker != NULL);
  RHEA_ASSERT (stokes_problem_nl->yielding_marker != NULL);
  RHEA_ASSERT (stokes_problem_nl->rhs_vel != NULL);
  RHEA_ASSERT (stokes_problem_nl->rhs_vel_press != NULL);

  /* destroy Stokes operator and preconditioner */
  if (stokes_problem_nl->stokes_pc != NULL) {
    ymir_stokes_pc_destroy (stokes_problem_nl->stokes_pc);
  }
//ymir_stokes_op_destroy (stokes_problem_nl->stokes_op); //TODO

  /* destroy vectors */
  rhea_viscosity_destroy (stokes_problem_nl->coeff);
  rhea_viscosity_destroy (stokes_problem_nl->rank1_tensor_scal);
  rhea_viscosity_destroy (stokes_problem_nl->bounds_marker);
  rhea_viscosity_destroy (stokes_problem_nl->yielding_marker);
  rhea_velocity_destroy (stokes_problem_nl->rhs_vel);
  rhea_velocity_pressure_destroy (stokes_problem_nl->rhs_vel_press);

  /* destroy velocity boundary conditions */
  rhea_stokes_problem_velocity_boundary_clear (stokes_problem_nl);

  /* destroy problem structure */
  rhea_stokes_problem_struct_destroy (stokes_problem_nl);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

void
rhea_stokes_problem_nonlinear_setup_solver (
                                    rhea_stokes_problem_t *stokes_problem_nl)
{
  const char         *this_fn_name =
                        "rhea_stokes_problem_nonlinear_setup_solver";
  ymir_stokes_op_t   *stokes_op = stokes_problem_nl->stokes_op;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);

  //TODO

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

static void
rhea_stokes_problem_nonlinear_apply_operator (ymir_vec_t *out, ymir_vec_t *in,
                                              void *data)
{
  rhea_stokes_problem_t *stokes_problem_nl = data;
  ymir_stokes_op_t   *stokes_op = stokes_problem_nl->stokes_op;
  const int           nl = 0;
  const int           dirty = 1;

  RHEA_ASSERT (stokes_problem_nl->stokes_op != NULL);
  ymir_stokes_pc_apply_stokes_op (in, out, stokes_op, nl, dirty);
}

static void
rhea_stokes_problem_nonlinear_apply_jacobian (ymir_vec_t *out, ymir_vec_t *in,
                                              void *data)
{
  rhea_stokes_problem_t *stokes_problem_nl = data;
  ymir_stokes_op_t   *stokes_op = stokes_problem_nl->stokes_op;
  const int           nl = 1;
  const int           dirty = 1;

  RHEA_ASSERT (stokes_problem_nl->stokes_op != NULL);
  ymir_stokes_pc_apply_stokes_op (in, out, stokes_op, nl, dirty);
}

static void
rhea_stokes_problem_nonlinear_update_operator (ymir_vec_t *sol_vel_press,
                                               void *data)
{
  const char         *this_fn_name =
                        "rhea_stokes_problem_nonlinear_update_operator";
  rhea_stokes_problem_t *stokes_problem_nl = data;
  ymir_stress_op_t   *stress_op = stokes_problem_nl->stokes_op->stress_op;

  RHEA_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->stokes_op != NULL);
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (sol_vel_press));

  /* compute physical viscosity and its derivative
   * (updates `state->vel_bc_vec`) */
//TODO
//slabs_physics_compute_stokes_coeff (viscosity, dvisc_dIIe, rank1_tensor_scal,
//                                    nl_stokes->bounds_marker,
//                                    nl_stokes->yielding_marker,
//                                    state, nl_stokes->press_elem, vel_dir,
//                                    physics_options);
//ymir_stress_op_setup_geo_coeff (stress_op);

  /* set derivative of viscosity to zero if Picard method is used */
//if (nl_solver_type == SL_NL_SOLVER_PICARD) {
//  slabs_nonlinear_stokes_op_switch_picard (nl_stokes);
//}

  RHEA_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

static void
rhea_stokes_problem_nonlinear_update_jacobian (ymir_vec_t *sol_vel_press,
                                               void *data)
{
  const char         *this_fn_name =
                        "rhea_stokes_problem_nonlinear_update_jacobian";
  rhea_stokes_problem_t *stokes_problem_nl = data;
  ymir_stress_op_t   *stress_op = stokes_problem_nl->stokes_op->stress_op;

  RHEA_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->stokes_op != NULL);
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (sol_vel_press));

  /* update 4th order tensor */
  //TODO setting usol is deprecated
//ymir_stress_op_set_usol (stress_op, state->vel_bc_vec, stress_op->dvdIIe);
//ymir_stress_op_coeff_compute_rank1_tensor (stress_op,
//                                           stress_op->coeff_rank1_tensor_scal,
//                                           state->vel_bc_vec);

  /* test derivative of Stokes operator */
//TODO
//if (check_derivative) {
//  slabs_nonlinear_stokes_problem_test_current_deriv (
//      nl_stokes->stokes_op, state, nl_stokes->mesh, nl_stokes->press_elem,
//      nl_stokes->vel_dir, nl_stokes->viscosity, physics_options);
//}

  RHEA_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

static void
rhea_stokes_problem_nonlinear_compute_residual (ymir_vec_t *residual,
                                                ymir_vec_t *solution,
                                                void *data)
{
  //TODO
}

static double
rhea_stokes_problem_nonlinear_compute_residual_norm (ymir_vec_t *residual,
                                                     void *data,
                                                     double *res_norm_comp)
{
  return 1.; //TODO
}

static int
rhea_stokes_problem_nonlinear_solve_linearized (
                                      ymir_vec_t *sol_vel_press,
                                      ymir_vec_t *rhs_vel_press,
                                      const int iter_max,
                                      const double rel_tol,
                                      const int nonzero_initial_guess,
                                      void *data,
                                      int *n_iter)
{
  const char         *this_fn_name =
                        "rhea_stokes_problem_nonlinear_solve_linearized";
  const double        abs_tol = 1.0e-100;
  const int           krylov_gmres_n_vecs = 100;
  rhea_stokes_problem_t *stokes_problem_nl = data;

  int                 itn, stop_reason;

  RHEA_GLOBAL_INFOF ("Into %s (rel tol %.3e, max iter %i)\n",
                     this_fn_name, rel_tol, iter_max);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->stokes_pc != NULL);
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (sol_vel_press));
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (rhs_vel_press));

  /* run solver for the linearized Stokes system */
  stop_reason = ymir_stokes_pc_solve (rhs_vel_press, sol_vel_press,
                                      stokes_problem_nl->stokes_pc,
                                      nonzero_initial_guess,
                                      rel_tol, abs_tol, iter_max,
                                      krylov_gmres_n_vecs /* unused */, &itn);

  /* set number of iterations as a return value */
  if (n_iter != NULL) {
    *n_iter = itn;
  }

  RHEA_GLOBAL_INFOF ("Done %s (solver stopping reason %i, num iter %i)\n",
                     this_fn_name, stop_reason, itn);

  /* return the termination reason of the solver */
  return stop_reason;
}

void
rhea_stokes_problem_nonlinear_solve (ymir_vec_t *sol_vel_press,
                                     const int iter_max,
                                     const double rel_tol,
                                     rhea_stokes_problem_t *stokes_problem_nl)
{
  const char         *this_fn_name = "rhea_stokes_problem_nonlinear_solve";
  const int           nonzero_initial_guess = 0;
  const double        abs_tol = 1.0e-100;
  const int           krylov_gmres_n_vecs = 100;
  const int           out_residual = !rhea_get_production_run ();

  double              norm_res_init, norm_res;
  int                 conv_reason;
  int                 n_iter;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s (rel tol %.3e, max iter %i)\n",
                           this_fn_name, rel_tol, iter_max);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (sol_vel_press));

  //TODO

  /* project out nullspaces of solution */
  //TODO
//slabs_linear_stokes_problem_project_out_nullspace (
//    state->vel_press_vec, lin_stokes->stokes_op, physics_options, 0);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}
