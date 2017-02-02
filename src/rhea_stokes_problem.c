/*
 */

#include <rhea_stokes_problem.h>
#include <rhea_base.h>
#include <rhea_newton.h>
#include <rhea_stokes_norm.h>
#include <rhea_velocity.h>
#include <rhea_velocity_pressure.h>
#include <ymir_stokes_pc.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_NORM_TYPE \
  (RHEA_STOKES_NORM_HMINUS1_L2)
#define RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_MASS_SCALING (0.0)
#define RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_CHECK_JACOBIAN (0)

/* initialize options */
int                 rhea_stokes_problem_nonlinear_norm_type =
                      RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_NORM_TYPE;
double              rhea_stokes_problem_nonlinear_norm_mass_scaling =
                      RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_MASS_SCALING;
int                 rhea_stokes_problem_nonlinear_check_jacobian =
                      RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_CHECK_JACOBIAN;

void
rhea_stokes_problem_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "StokesProblem";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_I, "nonlinear-norm-type", '\0',
    &(rhea_stokes_problem_nonlinear_norm_type),
    RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_NORM_TYPE,
    "Norm for convergence checks within nonlinear solver",
  YMIR_OPTIONS_D, "nonlinear-norm-Hminus1-mass-scaling", '\0',
    &(rhea_stokes_problem_nonlinear_norm_mass_scaling),
    RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_MASS_SCALING,
    "Scalar weight factor multiplying mass matrix of H^-1 norm operator",
  YMIR_OPTIONS_I, "nonlinear-check-jacobian", '\0',
    &(rhea_stokes_problem_nonlinear_check_jacobian),
    RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_CHECK_JACOBIAN,
    "Check Jacobians during the solve of a nonlinear problem with Newton",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);
}

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
  ymir_vec_t         *sol_vel;            /* nonlinear Stokes only */
  ymir_vec_t         *rank1_tensor_scal;  /* nonlinear Stokes only */
  ymir_vec_t         *bounds_marker;      /* nonlinear Stokes only */
  ymir_vec_t         *yielding_marker;    /* nonlinear Stokes only */
  ymir_vel_dir_t     *vel_dir;
  ymir_vec_t         *rhs_vel;
  ymir_vec_t         *rhs_vel_press;

  /* additional data of the Stokes problem (not owned) */
  ymir_vec_t         *temperature;
  ymir_vec_t         *weakzone;
  ymir_vec_t         *rhs_vel_nonzero_dirichlet;
  rhea_domain_options_t *domain_options;
  rhea_viscosity_options_t *visc_options;

  /* Stokes operator and preconditioner */
  ymir_stokes_op_t   *stokes_op;
  ymir_stokes_pc_t   *stokes_pc;

  /* Newton problem (nonlinear Stokes only) */
  rhea_newton_options_t  *newton_options;
  rhea_newton_problem_t  *newton_problem;
  rhea_stokes_norm_type_t norm_type;
  ymir_Hminus1_norm_op_t *norm_op;
  double                  norm_op_mass_scaling;
};

/**
 * Creates and initializes a new structure of a Stokes problem.
 */
static rhea_stokes_problem_t *
rhea_stokes_problem_struct_new (const rhea_stokes_problem_type_t type,
                                ymir_mesh_t *ymir_mesh,
                                ymir_pressure_elem_t *press_elem,
                                ymir_vec_t *rhs_vel_nonzero_dirichlet,
                                ymir_vec_t *temperature,
                                ymir_vec_t *weakzone,
                                rhea_domain_options_t *domain_options,
                                rhea_viscosity_options_t *visc_options)
{
  rhea_stokes_problem_t *stokes_problem = RHEA_ALLOC (rhea_stokes_problem_t, 1);

  stokes_problem->type = type;
  stokes_problem->ymir_mesh = ymir_mesh;
  stokes_problem->press_elem = press_elem;

  stokes_problem->coeff = NULL;
  stokes_problem->sol_vel = NULL;
  stokes_problem->rank1_tensor_scal = NULL;
  stokes_problem->bounds_marker = NULL;
  stokes_problem->yielding_marker = NULL;
  stokes_problem->vel_dir = NULL;
  stokes_problem->rhs_vel = NULL;
  stokes_problem->rhs_vel_press = NULL;

  stokes_problem->temperature = temperature;
  stokes_problem->weakzone = weakzone;
  stokes_problem->rhs_vel_nonzero_dirichlet = rhs_vel_nonzero_dirichlet;
  stokes_problem->domain_options = domain_options;
  stokes_problem->visc_options = visc_options;

  stokes_problem->stokes_op = NULL;
  stokes_problem->stokes_pc = NULL;

  stokes_problem->newton_options = NULL;
  stokes_problem->newton_problem = NULL;
  stokes_problem->norm_type = RHEA_STOKES_NORM_NONE;
  stokes_problem->norm_op = NULL;
  stokes_problem->norm_op_mass_scaling = NAN;

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
                         ymir_vec_t *rhs_vel_nonzero_dirichlet,
                         ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         rhea_domain_options_t *domain_options,
                         rhea_temperature_options_t *temp_options,
                         rhea_viscosity_options_t *visc_options,
                         void *solver_options)
{
  if (RHEA_VISCOSITY_NONLINEAR == visc_options->type) {
    return rhea_stokes_problem_nonlinear_new (
        temperature, weakzone, rhs_vel_nonzero_dirichlet,
        ymir_mesh, press_elem,
        domain_options, temp_options, visc_options, solver_options);
  }
  else {
    return rhea_stokes_problem_linear_new (
        temperature, weakzone, rhs_vel_nonzero_dirichlet,
        ymir_mesh, press_elem,
        domain_options, temp_options, visc_options, solver_options);
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
    rhea_stokes_problem_nonlinear_solve (sol_vel_press, iter_max, rel_tol,
                                         stokes_problem);
  }
  else {
    rhea_stokes_problem_linear_solve (sol_vel_press, iter_max, rel_tol,
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

ymir_vec_t *
rhea_stokes_problem_get_rhs_vel_nonzero_dirichlet (
                                        rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->rhs_vel_nonzero_dirichlet;
}

ymir_vec_t *
rhea_stokes_problem_get_rhs_vel (rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->rhs_vel;
}

ymir_stokes_op_t *
rhea_stokes_problem_get_stokes_op (rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->stokes_op;
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
                                ymir_vec_t *rhs_vel_nonzero_dirichlet,
                                ymir_mesh_t *ymir_mesh,
                                ymir_pressure_elem_t *press_elem,
                                rhea_domain_options_t *domain_options,
                                rhea_temperature_options_t *temp_options,
                                rhea_viscosity_options_t *visc_options,
                                void *solver_options)
{
  const char         *this_fn_name = "rhea_stokes_problem_linear_new";
  rhea_stokes_problem_t *stokes_problem_lin;
  ymir_vec_t         *coeff;
  ymir_vec_t         *rhs_vel;
  ymir_vec_t         *rhs_vel_press;
  ymir_stokes_op_t   *stokes_op;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (visc_options->type != RHEA_VISCOSITY_NONLINEAR);

  /* initialize Stokes problem structure */
  stokes_problem_lin = rhea_stokes_problem_struct_new (
      RHEA_STOKES_PROBLEM_LINEAR, ymir_mesh, press_elem,
      temperature, weakzone, rhs_vel_nonzero_dirichlet,
      domain_options, visc_options);

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

    /* transform viscosity to Stokes coefficient */
    coeff = viscosity;
    ymir_vec_scale (2.0, coeff);
  }

  /* create Stokes operator */
  stokes_op = ymir_stokes_op_new_ext (coeff, stokes_problem_lin->vel_dir,
                                      NULL /* Robin BC's */,
                                      NULL /* deprecated */,
                                      NULL /* deprecated */,
                                      press_elem, domain_options->center,
                                      domain_options->moment_of_inertia);

  /* create velocity right-hand side forcing */
  rhs_vel = rhea_velocity_new (ymir_mesh);
  rhea_temperature_compute_rhs_vel (rhs_vel, temperature, temp_options);

  /* construct right-hand side for incompressible Stokes system */
  rhs_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);
  ymir_stokes_op_construct_rhs_ext (
      rhs_vel /* Dirichlet forcing */,
      NULL /* Neumann forcing */,
      rhs_vel_nonzero_dirichlet /* nonzero Dirichlet bndr. */,
      rhs_vel_press,
      1 /* incompressible */,
      stokes_op);
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (rhs_vel_press));
  RHEA_ASSERT (ymir_vec_is_not_dirty (rhs_vel_press));

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
  RHEA_ASSERT (stokes_problem_lin->sol_vel == NULL);
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
  ymir_vec_t         *rhs_vel_press = stokes_problem_lin->rhs_vel_press;

  int                 stop_reason, itn;
  double              norm_res_init = NAN, norm_res = NAN;
  double              norm_res_init_vel = NAN, norm_res_vel = NAN;
  double              norm_res_init_press = NAN, norm_res_press = NAN;
  double              res_reduction = NAN, conv_factor = NAN;
  double              res_reduction_vel = NAN, res_reduction_press = NAN;

  RHEA_GLOBAL_INFOF ("Into %s (max iter %i, rel tol %.3e)\n",
                     this_fn_name, iter_max, rel_tol);

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);
  RHEA_ASSERT (stokes_problem_lin->stokes_pc != NULL);
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (sol_vel_press));

  /* compute residual with zero inital guess */
  if (out_residual) {
    norm_res_init = rhea_stokes_norm_compute (
        &norm_res_init_vel, &norm_res_init_press, rhs_vel_press,
        RHEA_STOKES_NORM_L2_VEC_SP, NULL, stokes_problem_lin->press_elem);
  }

  /* run solver for Stokes */
  stop_reason = ymir_stokes_pc_solve (rhs_vel_press, sol_vel_press,
                                      stokes_problem_lin->stokes_pc,
                                      nonzero_initial_guess,
                                      rel_tol, abs_tol, iter_max,
                                      krylov_gmres_n_vecs /* unused */,
                                      &itn);

  /* compute residual at solution */
  if (out_residual) {
    ymir_vec_t         *residual_vel_press = ymir_vec_template (rhs_vel_press);

    ymir_stokes_pc_apply_stokes_op (sol_vel_press, residual_vel_press,
                                    stokes_problem_lin->stokes_op, 0, 0);
    ymir_vec_add (-1.0, rhs_vel_press, residual_vel_press);
    ymir_vec_scale (-1.0, residual_vel_press);
    norm_res = rhea_stokes_norm_compute (
        &norm_res_vel, &norm_res_press, residual_vel_press,
        RHEA_STOKES_NORM_L2_VEC_SP, NULL, stokes_problem_lin->press_elem);
    ymir_vec_destroy (residual_vel_press);

    /* calculate convergence */
    res_reduction_vel = norm_res_vel/norm_res_init_vel;
    res_reduction_press = norm_res_press/norm_res_init_press;
    res_reduction = norm_res/norm_res_init;
    conv_factor = exp (log (res_reduction) / ((double) itn));
  }

  /* project out nullspaces of solution */
  //TODO
//slabs_linear_stokes_problem_project_out_nullspace (
//    state->vel_press_vec, lin_stokes->stokes_op, physics_options, 0);

  /* print output */
  if (out_residual) {
    RHEA_GLOBAL_INFOF ("Done %s (solver stopping reason %i, iterations %i, "
                       "residual reduction %.3e, convergence %.3e, "
                       "reduction vel %.3e, reduction press %.3e)\n",
                       this_fn_name, stop_reason, itn,
                       res_reduction, conv_factor,
                       res_reduction_vel, res_reduction_press);
  }
  else {
    RHEA_GLOBAL_INFOF ("Done %s (solver stopping reason %i, num iter %i)\n",
                       this_fn_name, stop_reason, itn);
  }
}

/******************************************************************************
 * Callback Functions for Newton Solver Used by Nonlinear Stokes Problem
 *****************************************************************************/

/**
 * Initializes user data of the nonlinear problem to be able to run the Newton
 * solver.
 * (Callback function for Newton's method.)
 */
static void
rhea_stokes_problem_nonlinear_data_init (ymir_vec_t *solution, void *data)
{
  const char         *this_fn_name = "rhea_stokes_problem_nonlinear_data_init";
  rhea_stokes_problem_t *stokes_problem_nl = data;
  const int           nonzero_initial_guess = (solution != NULL);
  ymir_vec_t         *coeff = stokes_problem_nl->coeff;
  ymir_vec_t         *sol_vel = stokes_problem_nl->sol_vel;
  ymir_vec_t         *rank1_tensor_scal = stokes_problem_nl->rank1_tensor_scal;
  ymir_vec_t         *bounds_marker = stokes_problem_nl->bounds_marker;
  ymir_vec_t         *yielding_marker = stokes_problem_nl->yielding_marker;
  ymir_vec_t         *temperature = stokes_problem_nl->temperature;
  ymir_vec_t         *weakzone = stokes_problem_nl->weakzone;

  RHEA_GLOBAL_VERBOSEF ("Into %s (nonzero initial guess %i)\n",
                        this_fn_name, nonzero_initial_guess);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->stokes_op == NULL);
  RHEA_ASSERT (stokes_problem_nl->stokes_pc == NULL);
  RHEA_ASSERT (stokes_problem_nl->norm_op == NULL);

  /* compute viscosity */
  if (nonzero_initial_guess) { /* if a velocity field is given */
    rhea_velocity_pressure_copy_components (sol_vel, NULL, solution,
                                            stokes_problem_nl->press_elem);
    rhea_viscosity_compute (
        coeff, rank1_tensor_scal, bounds_marker, yielding_marker,
        temperature, weakzone, sol_vel, stokes_problem_nl->visc_options);
  }
  else { /* otherwise assume zero velocity */
    ymir_vec_set_zero (sol_vel);
    rhea_viscosity_compute_init_nonlinear (
        coeff, rank1_tensor_scal, bounds_marker, yielding_marker,
        temperature, weakzone, stokes_problem_nl->visc_options);
  }

  /* transform viscosity to Stokes coefficient */
  ymir_vec_scale (2.0, coeff);

  /* create Stokes operator */
  stokes_problem_nl->stokes_op = ymir_stokes_op_new_ext (
      coeff, stokes_problem_nl->vel_dir,
      NULL /* Robin BC's */,
      NULL /* deprecated */,
      NULL /* deprecated */,
      stokes_problem_nl->press_elem,
      stokes_problem_nl->domain_options->center,
      stokes_problem_nl->domain_options->moment_of_inertia);

  /* set 4th-order tensor part of the coefficient */
  ymir_stress_op_coeff_compute_rank1_tensor (
      stokes_problem_nl->stokes_op->stress_op, rank1_tensor_scal, sol_vel);
  RHEA_ASSERT (ymir_stress_op_is_nl (stokes_problem_nl->stokes_op->stress_op));

  /* create Stokes preconditioner */
  stokes_problem_nl->stokes_pc = ymir_nlstokes_pc_new (
      stokes_problem_nl->stokes_op);

  /* create H^-1 norm operator */
  if (stokes_problem_nl->norm_type == RHEA_STOKES_NORM_HMINUS1_L2) {
    stokes_problem_nl->norm_op = ymir_Hminus1_norm_op_new (
        stokes_problem_nl->ymir_mesh,
        stokes_problem_nl->norm_op_mass_scaling);
  }

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
}

/**
 * Clears the data of the nonlinear problem after Newton solve.
 * (Callback function for Newton's method.)
 */
static void
rhea_stokes_problem_nonlinear_data_clear (void *data)
{
  const char         *this_fn_name = "rhea_stokes_problem_nonlinear_data_clear";
  rhea_stokes_problem_t *stokes_problem_nl = data;

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", this_fn_name);

  /* destroy Stokes preconditioner */
  if (stokes_problem_nl->stokes_pc != NULL) {
    ymir_stokes_pc_destroy (stokes_problem_nl->stokes_pc);
    stokes_problem_nl->stokes_pc = NULL;
  }

  /* destroy Stokes operator */
  if (stokes_problem_nl->stokes_op != NULL) {
    ymir_stokes_op_destroy (stokes_problem_nl->stokes_op);
    stokes_problem_nl->stokes_op = NULL;
  }

  /* destroy H^-1 norm operator */
  if (stokes_problem_nl->norm_type == RHEA_STOKES_NORM_HMINUS1_L2) {
    RHEA_ASSERT (stokes_problem_nl->norm_op != NULL);
    ymir_Hminus1_norm_op_destroy (stokes_problem_nl->norm_op);
  }

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
}

/**
 * Computes the negative gradient of the objective functional.
 * (Callback function for Newton's method.)
 */
static void
rhea_stokes_problem_nonlinear_compute_negative_gradient (
                                            ymir_vec_t *neg_gradient,
                                            ymir_vec_t *solution, void *data)
{
  const char         *this_fn_name =
    "rhea_stokes_problem_nonlinear_compute_negative_gradient";
  rhea_stokes_problem_t *stokes_problem_nl = data;
  ymir_stokes_op_t   *stokes_op = stokes_problem_nl->stokes_op;
  ymir_vec_t         *rhs_vel = stokes_problem_nl->rhs_vel;
  ymir_vec_t         *rhs_vel_nonzero_dirichlet =
                        stokes_problem_nl->rhs_vel_nonzero_dirichlet;
  ymir_vec_t         *rhs_vel_press = stokes_problem_nl->rhs_vel_press;

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", this_fn_name);

  /* check input */
#if RHEA_ENABLE_DEBUG
  RHEA_ASSERT (neg_gradient != NULL);
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (neg_gradient));
  if (solution != NULL) {
    RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (solution));
    RHEA_ASSERT (rhea_velocity_pressure_is_valid (solution));
    RHEA_ASSERT (ymir_vec_is_not_dirty (solution));
  }
  RHEA_ASSERT (stokes_op != NULL);
  RHEA_ASSERT (rhs_vel_press != NULL);
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (rhs_vel_press));
#endif

  /* construct the right-hand side */
  ymir_stokes_op_construct_rhs_ext (
      rhs_vel /* Dirichlet forcing */,
      NULL /* Neumann forcing */,
      rhs_vel_nonzero_dirichlet /* nonzero Dirichlet bndr. */,
      rhs_vel_press,
      1 /* incompressible */,
      stokes_op);
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (rhs_vel_press));
  RHEA_ASSERT (ymir_vec_is_not_dirty (rhs_vel_press));

  /* compute the residual (which assumes the role of the negative gradient) */
  if (solution != NULL) {
    /* compute Stokes residual
     *   r = b - F * x
     * where
     *   r --- residual
     *   b --- right-hand side
     *   F --- Stokes operator
     *   x --- velocity-pressure vector */
    ymir_stokes_pc_apply_stokes_op (solution, neg_gradient, stokes_op, 0, 0);
    ymir_vec_add (-1.0, rhs_vel_press, neg_gradient);
    ymir_vec_scale (-1.0, neg_gradient);
  }
  else {
    /* set residual to be the right-hand side */
    ymir_vec_copy (rhs_vel_press, neg_gradient);
  }
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (neg_gradient));
  RHEA_ASSERT (ymir_vec_is_not_dirty (neg_gradient));

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
}

/**
 * Computes the norm of the gradient.
 * (Callback function for Newton's method.)
 */
static double
rhea_stokes_problem_nonlinear_compute_gradient_norm (ymir_vec_t *neg_gradient,
                                                     void *data,
                                                     double *norm_comp)
{
  const char         *this_fn_name =
                        "rhea_stokes_problem_nonlinear_compute_gradient_norm";
  rhea_stokes_problem_t *stokes_problem_nl = data;
  rhea_stokes_norm_type_t norm_type = stokes_problem_nl->norm_type;
  ymir_Hminus1_norm_op_t *norm_op = stokes_problem_nl->norm_op;
  ymir_vec_t         *innerprod_arg_left, *innerprod_arg_right;
  double              ip_vel, ip_press;

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (neg_gradient != NULL);
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (neg_gradient));
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (neg_gradient));
  RHEA_ASSERT (stokes_problem_nl->norm_type != RHEA_STOKES_NORM_NONE);
  RHEA_ASSERT (stokes_problem_nl->press_elem != NULL);

  /* compute inner products */
  if (ymir_stokes_pc_block_type == YMIR_STOKES_PC_BLOCK_L) { /* if left PC */
    ymir_stokes_pc_t   *stokes_pc = stokes_problem_nl->stokes_pc;

    RHEA_ASSERT (stokes_problem_nl->stokes_pc != NULL);
    RHEA_ASSERT (norm_type != RHEA_STOKES_NORM_HMINUS1_L2); // not implemented

    innerprod_arg_right = ymir_vec_template (neg_gradient);
    if (ymir_stress_op_is_nl (stokes_pc->stokes_op->stress_op)) {
      ymir_nlstokes_pc_apply_lower (neg_gradient, innerprod_arg_right,
                                    stokes_pc);
    }
    else {
      ymir_stokes_pc_apply_lower (neg_gradient, innerprod_arg_right, stokes_pc);
    }
    innerprod_arg_left = innerprod_arg_right;
  }
  else { /* otherwise right or symmetric preconditioning */
    innerprod_arg_left = neg_gradient;
    innerprod_arg_right = neg_gradient;
  }
  rhea_stokes_norm_innerprod (
      &ip_vel, &ip_press, innerprod_arg_left, innerprod_arg_right,
      norm_type, norm_op, stokes_problem_nl->press_elem);

  /* destroy */
  if (ymir_stokes_pc_block_type == YMIR_STOKES_PC_BLOCK_L) {
    ymir_vec_destroy (innerprod_arg_right);
  }

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);

  /* calculate norms */
  if (norm_comp != NULL) {
    norm_comp[0] = sqrt (ip_vel);
    norm_comp[1] = sqrt (ip_press);
  }
  return sqrt (ip_vel + ip_press);
}

/**
 * Applies Hessian operator.
 * (Callback function for Newton's method.)
 */
static void
rhea_stokes_problem_nonlinear_apply_hessian (ymir_vec_t *out, ymir_vec_t *in,
                                             void *data)
{
  rhea_stokes_problem_t *stokes_problem_nl = data;
  ymir_stokes_op_t   *stokes_op = stokes_problem_nl->stokes_op;
  const int           nl = 1;
  const int           dirty = 0;

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->stokes_op != NULL);
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (out));
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (in));
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (in));

  /* apply */
  ymir_stokes_pc_apply_stokes_op (in, out, stokes_op, nl, dirty);
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (out));
}

/**
 * Applies inexact inverse of the Hessian operator.
 * (Callback function for Newton's method.)
 */
static int
rhea_stokes_problem_nonlinear_solve_hessian_system (
                                            ymir_vec_t *step,
                                            ymir_vec_t *neg_gradient,
                                            const int lin_iter_max,
                                            const double lin_res_norm_rtol,
                                            const int nonzero_initial_guess,
                                            void *data,
                                            int *lin_iter_count)
{
  const char         *this_fn_name =
                        "rhea_stokes_problem_nonlinear_solve_hessian_system";
  const double        lin_res_abs_tol = 1.0e-100;
  const int           krylov_gmres_restart = 100;
  rhea_stokes_problem_t *stokes_problem_nl = data;
  int                 itn, stop_reason;

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->stokes_pc != NULL);
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (step));
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (neg_gradient));
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (neg_gradient));

  /* run solver for the linearized Stokes system */
  stop_reason = ymir_stokes_pc_solve (neg_gradient, step,
                                      stokes_problem_nl->stokes_pc,
                                      nonzero_initial_guess,
                                      lin_res_norm_rtol, lin_res_abs_tol,
                                      lin_iter_max,
                                      krylov_gmres_restart /* unused */, &itn);
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (step));

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);

  /* return iteraton count and "stopping" reason */
  if (lin_iter_count != NULL) {
    *lin_iter_count = itn;
  }
  return stop_reason;
}

/**
 * Updates the nonlinear operator at the current solution vector.
 * (Callback function for Newton's method.)
 */
static void
rhea_stokes_problem_nonlinear_update_operator (ymir_vec_t *solution, void *data)
{
  const char         *this_fn_name =
                        "rhea_stokes_problem_nonlinear_update_operator";
  rhea_stokes_problem_t *stokes_problem_nl = data;
  ymir_vec_t         *coeff = stokes_problem_nl->coeff;
  ymir_vec_t         *sol_vel = stokes_problem_nl->sol_vel;
  ymir_vec_t         *rank1_tensor_scal = stokes_problem_nl->rank1_tensor_scal;
  ymir_vec_t         *bounds_marker = stokes_problem_nl->bounds_marker;
  ymir_vec_t         *yielding_marker = stokes_problem_nl->yielding_marker;
  ymir_vec_t         *temperature = stokes_problem_nl->temperature;
  ymir_vec_t         *weakzone = stokes_problem_nl->weakzone;

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->stokes_op != NULL);
  RHEA_ASSERT (ymir_stress_op_is_nl (stokes_problem_nl->stokes_op->stress_op));
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (solution));

  /* retrieve velocity */
  rhea_velocity_pressure_copy_components (sol_vel, NULL, solution,
                                          stokes_problem_nl->press_elem);

  /* compute viscosity; transform viscosity to Stokes coefficient */
  rhea_viscosity_compute (
      coeff, rank1_tensor_scal, bounds_marker, yielding_marker,
      temperature, weakzone, sol_vel, stokes_problem_nl->visc_options);
  ymir_vec_scale (2.0, coeff);

  /* re-setup coefficient data of viscous stress operator */
  ymir_stress_op_setup_geo_coeff (stokes_problem_nl->stokes_op->stress_op);

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
}

/**
 * Updates the Hessian operator at the current solution vector.
 * (Callback function for Newton's method.)
 */
static void
rhea_stokes_problem_nonlinear_update_hessian (ymir_vec_t *solution, void *data)
{
  const char         *this_fn_name =
                        "rhea_stokes_problem_nonlinear_update_hessian";
  rhea_stokes_problem_t *stokes_problem_nl = data;
  ymir_stress_op_t   *stress_op = stokes_problem_nl->stokes_op->stress_op;
  ymir_vec_t         *sol_vel = stokes_problem_nl->sol_vel;
  ymir_vec_t         *rank1_tensor_scal = stokes_problem_nl->rank1_tensor_scal;

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->stokes_op != NULL);
  RHEA_ASSERT (stokes_problem_nl->stokes_pc != NULL);
  RHEA_ASSERT (ymir_stress_op_is_nl (stokes_problem_nl->stokes_op->stress_op));
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (solution));

  /* retrieve velocity */
  rhea_velocity_pressure_copy_components (sol_vel, NULL, solution,
                                          stokes_problem_nl->press_elem);

  /* set 4th-order tensor part of the coefficient */
  ymir_stress_op_coeff_compute_rank1_tensor (
      stokes_problem_nl->stokes_op->stress_op, rank1_tensor_scal, sol_vel);

  /* update Stokes preconditioner */
#if YMIR_WITH_PETSC
  ymir_nlstokes_pc_recompute (stokes_problem_nl->stokes_pc);
#else
  RHEA_ABORT_NOT_REACHED ();
#endif

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
}

/******************************************************************************
 * Nonlinear Stokes Problem
 *****************************************************************************/

rhea_stokes_problem_t *
rhea_stokes_problem_nonlinear_new (ymir_vec_t *temperature,
                                   ymir_vec_t *weakzone,
                                   ymir_vec_t *rhs_vel_nonzero_dirichlet,
                                   ymir_mesh_t *ymir_mesh,
                                   ymir_pressure_elem_t *press_elem,
                                   rhea_domain_options_t *domain_options,
                                   rhea_temperature_options_t *temp_options,
                                   rhea_viscosity_options_t *visc_options,
                                   void *solver_options)
{
  const char         *this_fn_name = "rhea_stokes_problem_nonlinear_new";
  rhea_stokes_problem_t *stokes_problem_nl;
  ymir_vec_t         *coeff;
  ymir_vec_t         *sol_vel;
  ymir_vec_t         *rank1_tensor_scal;
  ymir_vec_t         *bounds_marker;
  ymir_vec_t         *yielding_marker;
  ymir_vec_t         *rhs_vel;
  ymir_vec_t         *rhs_vel_press;
  rhea_newton_options_t *newton_options = solver_options;
  rhea_newton_problem_t *newton_problem;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (visc_options->type == RHEA_VISCOSITY_NONLINEAR);

  /* initialize Stokes problem structure */
  stokes_problem_nl = rhea_stokes_problem_struct_new (
      RHEA_STOKES_PROBLEM_NONLINEAR, ymir_mesh, press_elem,
      temperature, weakzone, rhs_vel_nonzero_dirichlet,
      domain_options, visc_options);

  /* create velocity boundary conditions */
  rhea_stokes_problem_velocity_boundary_create (stokes_problem_nl,
                                                ymir_mesh, domain_options);

  /* create coefficient vectors */
  coeff = rhea_viscosity_new (ymir_mesh);
  sol_vel = rhea_velocity_new (ymir_mesh);
  rank1_tensor_scal = rhea_viscosity_new (ymir_mesh);
  bounds_marker = rhea_viscosity_new (ymir_mesh);
  yielding_marker = rhea_viscosity_new (ymir_mesh);

  /* create velocity right-hand side forcing */
  rhs_vel = rhea_velocity_new (ymir_mesh);
  rhea_temperature_compute_rhs_vel (rhs_vel, temperature, temp_options);

  /* create right-hand side vector */
  rhs_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);

  /* create Newton problem */
  {
    ymir_vec_t         *step_vec = rhea_velocity_pressure_new (
                            ymir_mesh, press_elem);
    ymir_vec_t         *neg_gradient_vec = rhea_velocity_pressure_new (
                            ymir_mesh, press_elem);

    newton_problem = rhea_newton_problem_new (
        neg_gradient_vec, step_vec,
        rhea_stokes_problem_nonlinear_compute_negative_gradient,
        rhea_stokes_problem_nonlinear_solve_hessian_system);

    rhea_newton_problem_set_data_fn (
        stokes_problem_nl,
        rhea_stokes_problem_nonlinear_data_init,
        rhea_stokes_problem_nonlinear_data_clear, newton_problem);
    rhea_newton_problem_set_conv_criterion_fn (
        RHEA_NEWTON_CONV_CRITERION_RESIDUAL_NORM,
        NULL /* objective functional is not provided */,
        rhea_stokes_problem_nonlinear_compute_gradient_norm,
        2 /* gradient norms have 2 components */, newton_problem);
    rhea_newton_problem_set_apply_hessian_fn (
        rhea_stokes_problem_nonlinear_apply_hessian, newton_problem);
    rhea_newton_problem_set_update_fn (
        rhea_stokes_problem_nonlinear_update_operator,
        rhea_stokes_problem_nonlinear_update_hessian, newton_problem);

    if (rhea_stokes_problem_nonlinear_check_jacobian) {
      rhea_newton_problem_set_checks (
          0 /* grad */, 1 /* Hessian */, newton_problem);
    }
  }

  /* fill and return the structure of the nonlinear Stokes problem */
  stokes_problem_nl->coeff = coeff;
  stokes_problem_nl->sol_vel = sol_vel;
  stokes_problem_nl->rank1_tensor_scal = rank1_tensor_scal;
  stokes_problem_nl->bounds_marker = bounds_marker;
  stokes_problem_nl->yielding_marker = yielding_marker;
  stokes_problem_nl->rhs_vel = rhs_vel;
  stokes_problem_nl->rhs_vel_nonzero_dirichlet = rhs_vel_nonzero_dirichlet;
  stokes_problem_nl->rhs_vel_press = rhs_vel_press;
  stokes_problem_nl->newton_options = newton_options;
  stokes_problem_nl->newton_problem = newton_problem;
  stokes_problem_nl->norm_type =
    (rhea_stokes_norm_type_t) rhea_stokes_problem_nonlinear_norm_type;
  stokes_problem_nl->norm_op_mass_scaling =
    rhea_stokes_problem_nonlinear_norm_mass_scaling;

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
  RHEA_ASSERT (stokes_problem_nl->newton_problem != NULL);
  RHEA_ASSERT (stokes_problem_nl->coeff != NULL);
  RHEA_ASSERT (stokes_problem_nl->sol_vel != NULL);
  RHEA_ASSERT (stokes_problem_nl->rank1_tensor_scal != NULL);
  RHEA_ASSERT (stokes_problem_nl->bounds_marker != NULL);
  RHEA_ASSERT (stokes_problem_nl->yielding_marker != NULL);
  RHEA_ASSERT (stokes_problem_nl->rhs_vel != NULL);
  RHEA_ASSERT (stokes_problem_nl->rhs_vel_press != NULL);

  /* destroy Newton problem */
  {
    rhea_newton_problem_t *newton_problem = stokes_problem_nl->newton_problem;

    ymir_vec_destroy (
        rhea_newton_problem_get_neg_gradient_vec (newton_problem));
    ymir_vec_destroy (
        rhea_newton_problem_get_step_vec (newton_problem));
    rhea_newton_problem_destroy (newton_problem);
  }

  /* destroy vectors */
  rhea_viscosity_destroy (stokes_problem_nl->coeff);
  rhea_viscosity_destroy (stokes_problem_nl->sol_vel);
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
  return;
}

void
rhea_stokes_problem_nonlinear_solve (ymir_vec_t *sol_vel_press,
                                     const int iter_max,
                                     const double rel_tol,
                                     rhea_stokes_problem_t *stokes_problem_nl)
{
  const char         *this_fn_name = "rhea_stokes_problem_nonlinear_solve";
  const int           nonzero_initial_guess = 0;
  const int           status_verbosity = 1;
  rhea_newton_options_t *newton_options = stokes_problem_nl->newton_options;
  rhea_newton_problem_t *newton_problem = stokes_problem_nl->newton_problem;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (sol_vel_press));

  /* run Newton solver */
  newton_options->nonzero_initial_guess = nonzero_initial_guess;
  newton_options->iter_max = newton_options->iter_start + iter_max;
  newton_options->rtol = rel_tol;
  newton_options->status_verbosity = status_verbosity;
  rhea_newton_solve (sol_vel_press, newton_problem, newton_options);

  /* project out nullspaces of solution */
  //TODO
//slabs_linear_stokes_problem_project_out_nullspace (
//    state->vel_press_vec, lin_stokes->stokes_op, physics_options, 0);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}
