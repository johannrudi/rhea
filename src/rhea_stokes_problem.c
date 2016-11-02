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
    return NULL;
    //TODO
    //return rhea_stokes_problem_nonlinear_new (...);
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
    //TODO
    //rhea_stokes_problem_nonlinear_destroy (stokes_problem);
  }
  else {
    rhea_stokes_problem_linear_destroy (stokes_problem);
  }
}

void
rhea_stokes_problem_solve (ymir_vec_t *sol_vel_press,
                           const double rel_tol,
                           const int maxiter,
                           rhea_stokes_problem_t *stokes_problem)
{
  if (RHEA_STOKES_PROBLEM_NONLINEAR == stokes_problem->type) {
    //TODO
    //rhea_stokes_problem_nonlinear_solve (stokes_problem);
  }
  else {
    rhea_stokes_problem_linear_solve (sol_vel_press, rel_tol, maxiter,
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
  ymir_vec_t         *viscosity, *coeff;
  ymir_vec_t         *dirscal = NULL; //TODO seems to be best choice, isn't it?
  ymir_vel_dir_t     *vel_dir;
  ymir_vec_t         *rhs_vel;
  ymir_vec_t         *rhs_vel_press;
  ymir_stokes_op_t   *stokes_op;
  ymir_stokes_pc_t   *stokes_pc;
  rhea_stokes_problem_t *stokes_problem_lin;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (visc_options->type == RHEA_VISCOSITY_CONST ||
               visc_options->type == RHEA_VISCOSITY_LINEAR);

  /* create viscosity */
  viscosity = rhea_viscosity_new (ymir_mesh);
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

  /* create Dirichlet boundary conditions */
  vel_dir = rhea_domain_create_velocity_dirichlet_bc (ymir_mesh, dirscal,
                                                      domain_options);

  /* create Stokes operator */
  stokes_op = ymir_stokes_op_new_ext (coeff, vel_dir,
                                      NULL /* Robin BC's */,
                                      NULL /* nl. Stokes input */,
                                      NULL /* nl. Stokes input */,
                                      press_elem, domain_options->center,
                                      domain_options->moment_of_inertia);

  /* build Stokes preconditioner */
  stokes_pc = ymir_stokes_pc_new (stokes_op);

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

  /* create, fill, and return the structure of the linear Stokes problem */
  stokes_problem_lin = rhea_stokes_problem_struct_new (
      RHEA_STOKES_PROBLEM_LINEAR, ymir_mesh, press_elem);
  stokes_problem_lin->coeff = coeff;
  stokes_problem_lin->vel_dir = vel_dir;
  stokes_problem_lin->rhs_vel = rhs_vel;
  stokes_problem_lin->rhs_vel_press = rhs_vel_press;
  stokes_problem_lin->stokes_op = stokes_op;
  stokes_problem_lin->stokes_pc = stokes_pc;

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

  /* destroy Stokes operator and preconditioner */
  ymir_stokes_pc_destroy (stokes_problem_lin->stokes_pc);
  ymir_stokes_op_destroy (stokes_problem_lin->stokes_op);

  /* destroy vectors */
  rhea_viscosity_destroy (stokes_problem_lin->coeff);
  rhea_velocity_destroy (stokes_problem_lin->rhs_vel);
  rhea_velocity_pressure_destroy (stokes_problem_lin->rhs_vel_press);

  /* destroy Dirichlet BC's */
  if (stokes_problem_lin->vel_dir != NULL) {
    if (stokes_problem_lin->vel_dir->scale != NULL) {
      ymir_vec_destroy (stokes_problem_lin->vel_dir->scale);
    }
    ymir_vel_dir_destroy (stokes_problem_lin->vel_dir);
  }

  /* destroy linear Stokes problem */
  rhea_stokes_problem_struct_destroy (stokes_problem_lin);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

void
rhea_stokes_problem_linear_solve (ymir_vec_t *sol_vel_press,
                                  const double rel_tol,
                                  const int maxiter,
                                  rhea_stokes_problem_t *stokes_problem_lin)
{
  const char         *this_fn_name = "rhea_stokes_problem_linear_solve";
  const int           nonzero_initial_guess = 0;
  const double        abs_tol = 1.0e-16;
  const int           krylov_gmres_n_vecs = 100;

  double              norm_res_init, norm_res;
  int                 conv_reason;
  int                 n_iter;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s (rel tol %.3e, max iter %i)\n",
                           this_fn_name, rel_tol, maxiter);

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);
  RHEA_ASSERT (stokes_problem_lin->stokes_pc != NULL);
  //TODO check if `sol_vel_press` is of correct type

  /* compute residual with zero inital guess */
  if (rhea_get_production_run ()) {
    norm_res_init = -1.0;
  //TODO
  //ymir_vec_set_zero (state->vel_press_vec);
  //norm_res_init = slabs_norm_compute_residual (
  //    NULL, NULL, NULL, NULL, lin_stokes->rhs_u_point,
  //    lin_stokes->stokes_op, lin_stokes->stokes_pc,
  //    solver_options->krylov_type, SL_NORM_VEC_L2, NULL);
  }

  /* run solver for Stokes */
  conv_reason = ymir_stokes_pc_solve (stokes_problem_lin->rhs_vel_press,
                                      sol_vel_press,
                                      stokes_problem_lin->stokes_pc,
                                      nonzero_initial_guess,
                                      rel_tol, abs_tol, maxiter,
                                      krylov_gmres_n_vecs /* unused */,
                                      &n_iter);

  /* compute residual at solution */
  if (rhea_get_production_run ()) {
    double              res_reduction, conv_factor;

    norm_res = 1.0;
  //TODO
  //norm_res = slabs_norm_compute_residual (
  //    residual_up, NULL, NULL, state->vel_press_vec,
  //    lin_stokes->rhs_u_point,
  //    lin_stokes->stokes_op, lin_stokes->stokes_pc,
  //    solver_options->krylov_type, SL_NORM_VEC_L2, NULL);

    /* calculate rate of convergence */
    res_reduction = norm_res/norm_res_init;
    conv_factor = exp (log (res_reduction) / ((double) n_iter));

    RHEA_GLOBAL_PRODUCTIONF (
        "%s: Linear solver converged reason %d, %d iterations, "
        "residual reduction %.3e, convergence %.3e\n",
        this_fn_name, conv_reason, n_iter, res_reduction, conv_factor);
  }
  else {
    RHEA_GLOBAL_PRODUCTIONF (
        "%s: Linear solver converged reason %d, %d iterations\n",
        this_fn_name, conv_reason, n_iter);
  }

  /* project out nullspaces of solution */
  //TODO
//slabs_linear_stokes_problem_project_out_nullspace (
//    state->vel_press_vec, lin_stokes->stokes_op, physics_options, 0);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/******************************************************************************
 * Nonlinear Stokes Problem
 *****************************************************************************/

