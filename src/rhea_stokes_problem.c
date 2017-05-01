/*
 */

#include <rhea_stokes_problem.h>
#include <rhea_base.h>
#include <rhea_newton.h>
#include <rhea_stokes_norm.h>
#include <rhea_velocity.h>
#include <rhea_velocity_pressure.h>
#include <rhea_strainrate.h>
#include <rhea_stress.h>
#include <rhea_vtk.h>
#include <ymir_stokes_pc.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_LINEARIZATION_TYPE \
  (RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON)
#define RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_NORM_TYPE \
  (RHEA_STOKES_NORM_HMINUS1_L2)
#define RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_MASS_SCALING (0.0)
#define RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_CHECK_JACOBIAN (0)

/* initialize options */
int                 rhea_stokes_problem_nonlinear_linearization_type =
                      RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_LINEARIZATION_TYPE;
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

  YMIR_OPTIONS_I, "nonlinear-linearization-type", '\0',
    &(rhea_stokes_problem_nonlinear_linearization_type),
    RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_LINEARIZATION_TYPE,
    "Type of linearization for nonlinear solver",
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
 * General Stokes Problem Structure
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
  /* Stokes problem type */
  rhea_stokes_problem_type_t  type;
  int                         incompressible;

  /* mesh (not owned) */
  ymir_mesh_t           *ymir_mesh;
  ymir_pressure_elem_t  *press_elem;

  /* input data of the Stokes problem (not owned) */
  ymir_vec_t         *temperature;
  ymir_vec_t         *weakzone;
  ymir_vec_t         *rhs_vel;
  ymir_vec_t         *rhs_vel_nonzero_dirichlet;
  rhea_domain_options_t *domain_options;
  rhea_viscosity_options_t *visc_options;

  /* data of the Stokes problem */
  ymir_vec_t         *coeff;
  ymir_vec_t         *bounds_marker;      /* nonlinear Stokes only */
  ymir_vec_t         *yielding_marker;    /* nonlinear Stokes only */
  ymir_vec_t         *sol_vel;            /* nonlinear Stokes only */
  ymir_vec_t         *rank1_tensor_scal;  /* nonlinear Stokes only */
  ymir_vec_t         *rank1_tensor;       /* nonlinear Stokes only */
  ymir_vel_dir_t     *vel_dir;
  ymir_vec_t         *rhs_vel_press;

  /* Stokes operator and preconditioner */
  ymir_stokes_op_t   *stokes_op;
  ymir_stokes_pc_t   *stokes_pc;

  /* Newton problem (nonlinear Stokes only) */
  rhea_stokes_problem_nonlinear_linearization_t linearization_type;
  rhea_newton_options_t  *newton_options;
  rhea_newton_problem_t  *newton_problem;
  rhea_stokes_norm_type_t norm_type;
  ymir_Hminus1_norm_op_t *norm_op;
  double                  norm_op_mass_scaling;
  char                   *vtk_write_newton_itn_path;
};

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

/**
 * Creates and initializes a new structure of a Stokes problem.
 */
static rhea_stokes_problem_t *
rhea_stokes_problem_struct_new (const rhea_stokes_problem_type_t type,
                                ymir_mesh_t *ymir_mesh,
                                ymir_pressure_elem_t *press_elem,
                                ymir_vec_t *temperature,
                                ymir_vec_t *weakzone,
                                ymir_vec_t *rhs_vel,
                                ymir_vec_t *rhs_vel_nonzero_dirichlet,
                                rhea_domain_options_t *domain_options,
                                rhea_viscosity_options_t *visc_options)
{
  rhea_stokes_problem_t *stokes_problem = RHEA_ALLOC (rhea_stokes_problem_t, 1);

  /* initialize */
  stokes_problem->type = type;
  stokes_problem->incompressible = 1;
  stokes_problem->ymir_mesh = ymir_mesh;
  stokes_problem->press_elem = press_elem;

  stokes_problem->temperature = temperature;
  stokes_problem->weakzone = weakzone;
  stokes_problem->rhs_vel = rhs_vel;
  stokes_problem->rhs_vel_nonzero_dirichlet = rhs_vel_nonzero_dirichlet;
  stokes_problem->domain_options = domain_options;
  stokes_problem->visc_options = visc_options;

  stokes_problem->coeff = NULL;
  stokes_problem->bounds_marker = NULL;
  stokes_problem->yielding_marker = NULL;
  stokes_problem->sol_vel = NULL;
  stokes_problem->rank1_tensor_scal = NULL;
  stokes_problem->rank1_tensor = NULL;
  stokes_problem->vel_dir = NULL;
  stokes_problem->rhs_vel_press = NULL;

  stokes_problem->stokes_op = NULL;
  stokes_problem->stokes_pc = NULL;

  stokes_problem->linearization_type =
    RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_LINEARIZATION_TYPE;
  stokes_problem->newton_options = NULL;
  stokes_problem->newton_problem = NULL;
  stokes_problem->norm_type = RHEA_STOKES_NORM_NONE;
  stokes_problem->norm_op = NULL;
  stokes_problem->norm_op_mass_scaling = NAN;
  stokes_problem->vtk_write_newton_itn_path = NULL;

  /* create velocity boundary conditions */
  rhea_stokes_problem_velocity_boundary_create (stokes_problem, ymir_mesh,
                                                domain_options);

  return stokes_problem;
}

/**
 * Destroys a structure of a Stokes problem.
 */
static void
rhea_stokes_problem_struct_destroy (rhea_stokes_problem_t *stokes_problem)
{
  /* destroy velocity boundary conditions */
  rhea_stokes_problem_velocity_boundary_clear (stokes_problem);

  /* destroy structure */
  RHEA_FREE (stokes_problem);
}

/******************************************************************************
 * Linear Stokes Problem
 *****************************************************************************/

static rhea_stokes_problem_t *
rhea_stokes_problem_linear_new (ymir_vec_t *temperature,
                                ymir_vec_t *weakzone,
                                ymir_vec_t *rhs_vel,
                                ymir_vec_t *rhs_vel_nonzero_dirichlet,
                                ymir_mesh_t *ymir_mesh,
                                ymir_pressure_elem_t *press_elem,
                                rhea_domain_options_t *domain_options,
                                rhea_viscosity_options_t *visc_options,
                                void *solver_options)
{
  const char         *this_fn_name = "rhea_stokes_problem_linear_new";
  rhea_stokes_problem_t *stokes_problem_lin;
  ymir_vec_t         *coeff;
  ymir_vec_t         *rhs_vel_press;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (visc_options->type == RHEA_VISCOSITY_LINEAR);

  /* initialize Stokes problem structure */
  stokes_problem_lin = rhea_stokes_problem_struct_new (
      RHEA_STOKES_PROBLEM_LINEAR, ymir_mesh, press_elem,
      temperature, weakzone, rhs_vel, rhs_vel_nonzero_dirichlet,
      domain_options, visc_options);

  /* create coefficient and right-hand side vectors */
  coeff = rhea_viscosity_new (ymir_mesh);
  rhs_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);

  /* fill the structure of the linear Stokes problem */
  stokes_problem_lin->coeff = coeff;
  stokes_problem_lin->rhs_vel_press = rhs_vel_press;

  /* compute viscosity */
  rhea_viscosity_compute (coeff,
                          NULL /* nl. Stokes output */,
                          NULL /* nl. Stokes output */,
                          NULL /* nl. Stokes output */,
                          temperature, weakzone,
                          NULL /* nl. Stokes input */,
                          visc_options);

  /* transform viscosity to viscous stress coefficient */
  ymir_vec_scale (2.0, coeff);

  /* create Stokes operator */
  stokes_problem_lin->stokes_op = ymir_stokes_op_new_ext (
      coeff, stokes_problem_lin->vel_dir,
      NULL /* Robin BC's */,
      NULL /* deprecated */,
      NULL /* deprecated */,
      press_elem, domain_options->center,
      domain_options->moment_of_inertia);

  /* set viscous stress coefficient */
  {
    ymir_stress_op_t   *stress_op = stokes_problem_lin->stokes_op->stress_op;

    ymir_stress_op_set_coeff_type (
        stress_op, YMIR_STRESS_OP_COEFF_ISOTROPIC_SCAL);
    ymir_stress_op_set_coeff_shift (
        stress_op, rhea_viscosity_get_visc_shift (visc_options));
  }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

  return stokes_problem_lin;
}

static void
rhea_stokes_problem_linear_destroy (rhea_stokes_problem_t *stokes_problem_lin)
{
  const char         *this_fn_name = "rhea_stokes_problem_linear_destroy";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);
  RHEA_ASSERT (stokes_problem_lin->stokes_op != NULL);
  RHEA_ASSERT (stokes_problem_lin->coeff != NULL);
  RHEA_ASSERT (stokes_problem_lin->bounds_marker == NULL);
  RHEA_ASSERT (stokes_problem_lin->yielding_marker == NULL);
  RHEA_ASSERT (stokes_problem_lin->sol_vel == NULL);
  RHEA_ASSERT (stokes_problem_lin->rank1_tensor_scal == NULL);
  RHEA_ASSERT (stokes_problem_lin->rank1_tensor == NULL);
  RHEA_ASSERT (stokes_problem_lin->rhs_vel_press != NULL);

  /* destroy Stokes preconditioner and operator */
  if (stokes_problem_lin->stokes_pc != NULL) {
    ymir_stokes_pc_destroy (stokes_problem_lin->stokes_pc);
  }
  ymir_stokes_op_destroy (stokes_problem_lin->stokes_op);

  /* destroy vectors */
  rhea_viscosity_destroy (stokes_problem_lin->coeff);
  rhea_velocity_pressure_destroy (stokes_problem_lin->rhs_vel_press);

  /* destroy problem structure */
  rhea_stokes_problem_struct_destroy (stokes_problem_lin);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

static void
rhea_stokes_problem_linear_setup_solver (
                                    rhea_stokes_problem_t *stokes_problem_lin)
{
  const char         *this_fn_name = "rhea_stokes_problem_linear_setup_solver";
  ymir_stokes_op_t   *stokes_op = stokes_problem_lin->stokes_op;
  ymir_vec_t         *rhs_vel_press = stokes_problem_lin->rhs_vel_press;
  ymir_vec_t         *rhs_vel = stokes_problem_lin->rhs_vel;
  ymir_vec_t         *rhs_vel_nonzero_dirichlet =
                        stokes_problem_lin->rhs_vel_nonzero_dirichlet;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);
  RHEA_ASSERT (stokes_problem_lin->stokes_op != NULL);
  RHEA_ASSERT (stokes_problem_lin->rhs_vel_press != NULL);

  /* construct right-hand side for Stokes system */
  ymir_stokes_pc_construct_rhs (
      rhs_vel_press             /* output: right-hand side */,
      rhs_vel                   /* input: volume forcing */,
      NULL                      /* input: Neumann forcing */,
      rhs_vel_nonzero_dirichlet /* input: nonzero Dirichlet boundary */,
      stokes_problem_lin->incompressible, stokes_op, 0 /* !linearized */);
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (rhs_vel_press));

  /* build Stokes preconditioner */
  stokes_problem_lin->stokes_pc = ymir_stokes_pc_new (stokes_op);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

static void
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

    /* calculate convergence (note: division by zero is possible) */
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
    if (stokes_problem_lin->incompressible) {
      RHEA_GLOBAL_INFOF ("Done %s (solver stopping reason %i, iterations %i, "
                         "residual reduction %.3e, convergence %.3e, "
                         "reduction momentum %.3e, norm mass %.3e)\n",
                         this_fn_name, stop_reason, itn,
                         res_reduction, conv_factor,
                         res_reduction_vel, norm_res_press);
    }
    else {
      RHEA_GLOBAL_INFOF ("Done %s (solver stopping reason %i, iterations %i, "
                         "residual reduction %.3e, convergence %.3e, "
                         "reduction momentum %.3e, reduction mass %.3e)\n",
                         this_fn_name, stop_reason, itn,
                         res_reduction, conv_factor,
                         res_reduction_vel, res_reduction_press);
    }
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
 * Sets the linearized part of the viscous stress coefficient.
 */
static void
rhea_stokes_problem_nonlinear_set_coeff_linearized (
                                      rhea_stokes_problem_t *stokes_problem_nl)
{
  rhea_viscosity_options_t *visc_options = stokes_problem_nl->visc_options;
  ymir_vec_t         *coeff = stokes_problem_nl->coeff;
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (coeff);
  ymir_stress_op_t   *stress_op;

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->stokes_op != NULL);

  /* get viscous stress operator */
  stress_op = stokes_problem_nl->stokes_op->stress_op;

  /* set linearized viscous stress coefficient */
  switch (stokes_problem_nl->linearization_type) {
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON:
    {
      ymir_vec_t         *sol_vel = stokes_problem_nl->sol_vel;
      ymir_vec_t         *tensor_scal = stokes_problem_nl->rank1_tensor_scal;
      ymir_vec_t         *rank1_tensor = stokes_problem_nl->rank1_tensor;

      /* check input */
      RHEA_ASSERT (stokes_problem_nl->sol_vel != NULL);
      RHEA_ASSERT (stokes_problem_nl->rank1_tensor_scal != NULL);
      RHEA_ASSERT (stokes_problem_nl->rank1_tensor != NULL);
      RHEA_ASSERT (ymir_stress_op_get_coeff_type_linearized (stress_op) ==
                   YMIR_STRESS_OP_COEFF_ANISOTROPIC_PPR1);

      /* compute and normalize strain rate tensor */
      ymir_stress_op_optimized_compute_strain_rate (
          rank1_tensor, sol_vel, stress_op);
      ymir_stress_op_tensor_normalize (rank1_tensor);

      /* pass coefficient data to viscous stress operator */
      ymir_stress_op_set_coeff_ppr1 (
          stress_op, coeff, tensor_scal, rank1_tensor);
    }
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL:
    {
      const double        yield_strength =
                            rhea_viscosity_get_yield_strength (visc_options);
      ymir_vec_t         *sol_vel = stokes_problem_nl->sol_vel;
      ymir_vec_t         *tensor_scal = stokes_problem_nl->rank1_tensor_scal;
      ymir_vec_t         *strain_tens = rhea_strainrate_new (ymir_mesh);
      ymir_vec_t         *stress_tens = rhea_stress_new (ymir_mesh);

      /* check input */
      RHEA_ASSERT (stokes_problem_nl->sol_vel != NULL);
      RHEA_ASSERT (stokes_problem_nl->rank1_tensor_scal != NULL);
      RHEA_ASSERT (ymir_stress_op_get_coeff_type_linearized (stress_op) ==
                   YMIR_STRESS_OP_COEFF_ANISOTROPIC_PPR2);

      /* compute and normalize strain rate tensor */
      ymir_stress_op_optimized_compute_strain_rate (
          strain_tens, sol_vel, stress_op);
      ymir_stress_op_tensor_normalize (strain_tens);

      /* compute and bound viscous stress tensor */
      ymir_stress_op_optimized_compute_visc_stress (
          stress_tens, sol_vel, stress_op, 1 /* linearized */);
      ymir_stress_op_tensor_bound_norm (stress_tens, yield_strength);
      ymir_vec_scale (1.0/yield_strength, stress_tens);

      /* pass coefficient data to viscous stress operator */
      ymir_stress_op_set_coeff_ppr2 (
          stress_op, coeff, tensor_scal, strain_tens, stress_tens);

      rhea_strainrate_destroy (strain_tens);
      rhea_stress_destroy (stress_tens);
    }
    break;
  default: /* unknown linearization type */
    RHEA_ABORT_NOT_REACHED ();
  }
}

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
  rhea_viscosity_options_t *visc_options = stokes_problem_nl->visc_options;
  const int           nonzero_initial_guess = (solution != NULL);
  ymir_vec_t         *coeff = stokes_problem_nl->coeff;
  ymir_vec_t         *bounds_marker = stokes_problem_nl->bounds_marker;
  ymir_vec_t         *yielding_marker = stokes_problem_nl->yielding_marker;
  ymir_vec_t         *sol_vel = stokes_problem_nl->sol_vel;
  ymir_vec_t         *tensor_scal = stokes_problem_nl->rank1_tensor_scal;
  ymir_vec_t         *temperature = stokes_problem_nl->temperature;
  ymir_vec_t         *weakzone = stokes_problem_nl->weakzone;
  ymir_stress_op_coeff_t  coeff_type, coeff_type_linearized;

  RHEA_GLOBAL_VERBOSEF ("Into %s (nonzero initial guess %i)\n",
                        this_fn_name, nonzero_initial_guess);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->stokes_op == NULL);
  RHEA_ASSERT (stokes_problem_nl->stokes_pc == NULL);
  RHEA_ASSERT (stokes_problem_nl->norm_op == NULL);

  /* set coefficient types */
  coeff_type = YMIR_STRESS_OP_COEFF_ISOTROPIC_SCAL;
  switch (stokes_problem_nl->linearization_type) {
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON:
    coeff_type_linearized = YMIR_STRESS_OP_COEFF_ANISOTROPIC_PPR1;
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL:
    coeff_type_linearized = YMIR_STRESS_OP_COEFF_ANISOTROPIC_PPR2;
    break;
  default: /* unknown linearization type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* compute viscosity */
  if (nonzero_initial_guess) { /* if a velocity field is given */
    rhea_velocity_pressure_copy_components (sol_vel, NULL, solution,
                                            stokes_problem_nl->press_elem);
    rhea_viscosity_compute (
        coeff, tensor_scal, bounds_marker, yielding_marker,
        temperature, weakzone, sol_vel, visc_options);
  }
  else { /* otherwise assume zero velocity */
    ymir_vec_set_zero (sol_vel);
    rhea_viscosity_compute_nonlinear_init (
        coeff, tensor_scal, bounds_marker, yielding_marker,
        temperature, weakzone, visc_options);
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

  /* set coefficient */
  {
    ymir_stress_op_t   *stress_op = stokes_problem_nl->stokes_op->stress_op;

    /* set viscous stress coefficient */
    ymir_stress_op_set_coeff_type (stress_op, coeff_type);
    ymir_stress_op_set_coeff_shift (
        stress_op, rhea_viscosity_get_visc_shift (visc_options));

    /* set the linearized part of the viscous stress coefficient */
    ymir_stress_op_set_coeff_type_linearized (
        stress_op, coeff_type_linearized);
    ymir_stress_op_set_coeff_shift_proj (
        stress_op, rhea_viscosity_get_visc_shift_proj (visc_options));
    rhea_stokes_problem_nonlinear_set_coeff_linearized (stokes_problem_nl);
  }

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

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);

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
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
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
  ymir_stokes_pc_construct_rhs (
      rhs_vel_press             /* output: right-hand side */,
      rhs_vel                   /* input: volume forcing */,
      NULL                      /* input: Neumann forcing */,
      rhs_vel_nonzero_dirichlet /* input: nonzero Dirichlet boundary */,
      stokes_problem_nl->incompressible, stokes_op, 0 /* !linearized */);
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (rhs_vel_press));

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
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
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
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->stokes_op != NULL);
  RHEA_ASSERT (ymir_stress_op_is_nl (stokes_problem_nl->stokes_op->stress_op));
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
  rhea_viscosity_options_t *visc_options = stokes_problem_nl->visc_options;
  ymir_vec_t         *coeff = stokes_problem_nl->coeff;
  ymir_vec_t         *bounds_marker = stokes_problem_nl->bounds_marker;
  ymir_vec_t         *yielding_marker = stokes_problem_nl->yielding_marker;
  ymir_vec_t         *sol_vel = stokes_problem_nl->sol_vel;
  ymir_vec_t         *tensor_scal = stokes_problem_nl->rank1_tensor_scal;
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
      coeff, tensor_scal, bounds_marker, yielding_marker,
      temperature, weakzone, sol_vel, visc_options);

  /* set viscous stress coefficient */
  {
    ymir_stress_op_t   *stress_op = stokes_problem_nl->stokes_op->stress_op;

    RHEA_ASSERT (ymir_stress_op_get_coeff_type (stress_op) ==
                 YMIR_STRESS_OP_COEFF_ISOTROPIC_SCAL);
    RHEA_ASSERT (ymir_stress_op_get_coeff_shift (stress_op) ==
                 rhea_viscosity_get_visc_shift (visc_options));

    ymir_vec_scale (2.0, coeff);
    ymir_stress_op_set_coeff_scal (stress_op, coeff);
  }

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
  rhea_viscosity_options_t *visc_options = stokes_problem_nl->visc_options;
  ymir_vec_t         *sol_vel = stokes_problem_nl->sol_vel;
  ymir_stress_op_t   *stress_op;

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->press_elem != NULL);
  RHEA_ASSERT (stokes_problem_nl->stokes_op != NULL);
  RHEA_ASSERT (stokes_problem_nl->stokes_pc != NULL);
  RHEA_ASSERT (ymir_stress_op_is_nl (stokes_problem_nl->stokes_op->stress_op));
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (solution));

  /* retrieve velocity */
  rhea_velocity_pressure_copy_components (sol_vel, NULL, solution,
                                          stokes_problem_nl->press_elem);

  /* set the linearized part of the viscous stress coefficient */
  stress_op = stokes_problem_nl->stokes_op->stress_op;
  RHEA_ASSERT (ymir_stress_op_get_coeff_shift_proj (stress_op) ==
               rhea_viscosity_get_visc_shift_proj (visc_options));
  rhea_stokes_problem_nonlinear_set_coeff_linearized (stokes_problem_nl);

  /* update Stokes preconditioner */
#if YMIR_WITH_PETSC
  ymir_nlstokes_pc_recompute (stokes_problem_nl->stokes_pc);
#else
  RHEA_ABORT_NOT_REACHED ();
#endif

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
}

/**
 * Writes or prints output at the beginning of a Newton step.
 * (Callback function for Newton's method.)
 */
static void
rhea_stokes_problem_nonlinear_output_prestep (ymir_vec_t *solution,
                                              const int iter, void *data)
{
  const char         *this_fn_name =
                        "rhea_stokes_problem_nonlinear_output_prestep";
  rhea_stokes_problem_t *stokes_problem_nl = data;
  const char         *filepath = stokes_problem_nl->vtk_write_newton_itn_path;
  ymir_mesh_t        *ymir_mesh = stokes_problem_nl->ymir_mesh;
  ymir_pressure_elem_t *press_elem = stokes_problem_nl->press_elem;
  ymir_vec_t         *velocity, *pressure, *viscosity;
  const double        domain_vol = stokes_problem_nl->domain_options->volume;
  double              bounds_vol_min, bounds_vol_max;
  double              yielding_vol;

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->bounds_marker != NULL);
  RHEA_ASSERT (stokes_problem_nl->yielding_marker != NULL);

  /* get fields */
  rhea_velocity_pressure_create_components (&velocity, &pressure, solution,
                                            press_elem);
  viscosity = rhea_viscosity_new (ymir_mesh);
  rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem_nl);

  /* compute statistics */
  rhea_viscosity_stats_get_bounds_volume (&bounds_vol_min, &bounds_vol_max,
                                          stokes_problem_nl->bounds_marker);
  yielding_vol = rhea_viscosity_stats_get_yielding_volume (
                                          stokes_problem_nl->yielding_marker);

  /* print statistics */
  RHEA_GLOBAL_INFOF (
      "%s: Bounds volume abs: min %.8e, max %.8e\n",
      this_fn_name, bounds_vol_min, bounds_vol_max);
  RHEA_GLOBAL_INFOF (
      "%s: Bounds volume rel: min %.3f, max %.3f\n",
      this_fn_name, bounds_vol_min/domain_vol, bounds_vol_max/domain_vol);
  RHEA_GLOBAL_INFOF (
      "%s: Yielding volume abs: %.8e\n", this_fn_name, yielding_vol);
  RHEA_GLOBAL_INFOF (
      "%s: Yielding volume rel: %.3f\n", this_fn_name, yielding_vol/domain_vol);

  /* write vtk */
  if (filepath != NULL) {
    char                path[BUFSIZ];

    snprintf (path, BUFSIZ, "%s_itn%02i", filepath, iter);
    rhea_vtk_write_nonlinear_stokes_iteration (
        path, velocity, pressure, viscosity,
        stokes_problem_nl->bounds_marker, stokes_problem_nl->yielding_marker);
  }

  /* destroy */
  rhea_velocity_destroy (velocity);
  rhea_pressure_destroy (pressure);
  rhea_viscosity_destroy (viscosity);

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
}

/******************************************************************************
 * Nonlinear Stokes Problem
 *****************************************************************************/

static rhea_stokes_problem_t *
rhea_stokes_problem_nonlinear_new (ymir_vec_t *temperature,
                                   ymir_vec_t *weakzone,
                                   ymir_vec_t *rhs_vel,
                                   ymir_vec_t *rhs_vel_nonzero_dirichlet,
                                   ymir_mesh_t *ymir_mesh,
                                   ymir_pressure_elem_t *press_elem,
                                   rhea_domain_options_t *domain_options,
                                   rhea_viscosity_options_t *visc_options,
                                   void *solver_options)
{
  const char         *this_fn_name = "rhea_stokes_problem_nonlinear_new";
  rhea_stokes_problem_t *stokes_problem_nl;
  ymir_vec_t         *coeff;
  ymir_vec_t         *bounds_marker;
  ymir_vec_t         *yielding_marker;
  ymir_vec_t         *sol_vel;
  ymir_vec_t         *rank1_tensor_scal = NULL;
  ymir_vec_t         *rank1_tensor = NULL;
  ymir_vec_t         *rhs_vel_press;
  rhea_newton_options_t *newton_options = solver_options;
  rhea_newton_problem_t *newton_problem;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (visc_options->type == RHEA_VISCOSITY_NONLINEAR);

  /* initialize Stokes problem structure */
  stokes_problem_nl = rhea_stokes_problem_struct_new (
      RHEA_STOKES_PROBLEM_NONLINEAR, ymir_mesh, press_elem,
      temperature, weakzone, rhs_vel, rhs_vel_nonzero_dirichlet,
      domain_options, visc_options);

  /* create vectors related to the coefficient */
  coeff = rhea_viscosity_new (ymir_mesh);
  bounds_marker = rhea_viscosity_new (ymir_mesh);
  yielding_marker = rhea_viscosity_new (ymir_mesh);
  sol_vel = rhea_velocity_new (ymir_mesh);

  /* create linearization data */
  switch (stokes_problem_nl->linearization_type) {
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON:
    rank1_tensor_scal = rhea_viscosity_new (ymir_mesh);
    rank1_tensor = rhea_strainrate_new (ymir_mesh);
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL:
    break;
  default: /* unknown linearization type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* create the right-hand side vector */
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
    rhea_newton_problem_set_output_fn (
        rhea_stokes_problem_nonlinear_output_prestep, newton_problem);

    if (rhea_stokes_problem_nonlinear_check_jacobian) {
      rhea_newton_problem_set_checks (
          0 /* grad */, 1 /* Hessian */, newton_problem);
    }
  }

  /* fill the structure of the nonlinear Stokes problem */
  stokes_problem_nl->coeff = coeff;
  stokes_problem_nl->bounds_marker = bounds_marker;
  stokes_problem_nl->yielding_marker = yielding_marker;
  stokes_problem_nl->sol_vel = sol_vel;
  stokes_problem_nl->rank1_tensor_scal = rank1_tensor_scal;
  stokes_problem_nl->rank1_tensor = rank1_tensor;
  stokes_problem_nl->rhs_vel_press = rhs_vel_press;
  stokes_problem_nl->linearization_type =
    (rhea_stokes_problem_nonlinear_linearization_t)
    rhea_stokes_problem_nonlinear_linearization_type;
  stokes_problem_nl->newton_options = newton_options;
  stokes_problem_nl->newton_problem = newton_problem;
  stokes_problem_nl->norm_type =
    (rhea_stokes_norm_type_t) rhea_stokes_problem_nonlinear_norm_type;
  stokes_problem_nl->norm_op_mass_scaling =
    rhea_stokes_problem_nonlinear_norm_mass_scaling;

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

  return stokes_problem_nl;
}

static void
rhea_stokes_problem_nonlinear_destroy (
                                    rhea_stokes_problem_t *stokes_problem_nl)
{
  const char         *this_fn_name = "rhea_stokes_problem_nonlinear_destroy";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->newton_problem != NULL);
  RHEA_ASSERT (stokes_problem_nl->coeff != NULL);
  RHEA_ASSERT (stokes_problem_nl->bounds_marker != NULL);
  RHEA_ASSERT (stokes_problem_nl->yielding_marker != NULL);
  RHEA_ASSERT (stokes_problem_nl->sol_vel != NULL);
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

  /* destroy linearization data */
  switch (stokes_problem_nl->linearization_type) {
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON:
    RHEA_ASSERT (stokes_problem_nl->rank1_tensor_scal != NULL);
    RHEA_ASSERT (stokes_problem_nl->rank1_tensor != NULL);
    rhea_viscosity_destroy (stokes_problem_nl->rank1_tensor_scal);
    rhea_strainrate_destroy (stokes_problem_nl->rank1_tensor);
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL:
    break;
  default: /* unknown linearization type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* destroy vectors */
  rhea_viscosity_destroy (stokes_problem_nl->coeff);
  rhea_viscosity_destroy (stokes_problem_nl->bounds_marker);
  rhea_viscosity_destroy (stokes_problem_nl->yielding_marker);
  rhea_velocity_destroy (stokes_problem_nl->sol_vel);
  rhea_velocity_pressure_destroy (stokes_problem_nl->rhs_vel_press);

  /* destroy problem structure */
  rhea_stokes_problem_struct_destroy (stokes_problem_nl);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

static void
rhea_stokes_problem_nonlinear_setup_solver (
                                    rhea_stokes_problem_t *stokes_problem_nl)
{
  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);

  /* setup happens during nonlinear solve */
  return;
}

static void
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

void
rhea_stokes_problem_nonlinear_set_output (
                              char *vtk_write_newton_iteration_path,
                              rhea_stokes_problem_t *stokes_problem_nl)
{
  stokes_problem_nl->vtk_write_newton_itn_path =
    vtk_write_newton_iteration_path;
}

/******************************************************************************
 * General Stokes Problem
 *****************************************************************************/

rhea_stokes_problem_t *
rhea_stokes_problem_new (ymir_vec_t *temperature,
                         ymir_vec_t *weakzone,
                         ymir_vec_t *rhs_vel,
                         ymir_vec_t *rhs_vel_nonzero_dirichlet,
                         ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         rhea_domain_options_t *domain_options,
                         rhea_viscosity_options_t *visc_options,
                         void *solver_options)
{
  switch (visc_options->type) {
  case RHEA_VISCOSITY_LINEAR:
    return rhea_stokes_problem_linear_new (
        temperature, weakzone, rhs_vel, rhs_vel_nonzero_dirichlet,
        ymir_mesh, press_elem, domain_options, visc_options, solver_options);
  case RHEA_VISCOSITY_NONLINEAR:
    return rhea_stokes_problem_nonlinear_new (
        temperature, weakzone, rhs_vel, rhs_vel_nonzero_dirichlet,
        ymir_mesh, press_elem, domain_options, visc_options, solver_options);
  default: /* not clear which Stokes type to choose */
    RHEA_ABORT_NOT_REACHED ();
  }
}

void
rhea_stokes_problem_destroy (rhea_stokes_problem_t *stokes_problem)
{
  switch (stokes_problem->type) {
  case RHEA_STOKES_PROBLEM_LINEAR:
    rhea_stokes_problem_linear_destroy (stokes_problem);
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR:
    rhea_stokes_problem_nonlinear_destroy (stokes_problem);
    break;
  default: /* unknown Stokes type */
    RHEA_ABORT_NOT_REACHED ();
  }
}

void
rhea_stokes_problem_setup_solver (rhea_stokes_problem_t *stokes_problem)
{
  switch (stokes_problem->type) {
  case RHEA_STOKES_PROBLEM_LINEAR:
    rhea_stokes_problem_linear_setup_solver (stokes_problem);
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR:
    rhea_stokes_problem_nonlinear_setup_solver (stokes_problem);
    break;
  default: /* unknown Stokes type */
    RHEA_ABORT_NOT_REACHED ();
  }
}

void
rhea_stokes_problem_solve (ymir_vec_t *sol_vel_press,
                           const int iter_max,
                           const double rel_tol,
                           rhea_stokes_problem_t *stokes_problem)
{
  switch (stokes_problem->type) {
  case RHEA_STOKES_PROBLEM_LINEAR:
    rhea_stokes_problem_linear_solve (sol_vel_press, iter_max, rel_tol,
                                      stokes_problem);
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR:
    rhea_stokes_problem_nonlinear_solve (sol_vel_press, iter_max, rel_tol,
                                         stokes_problem);
    break;
  default: /* unknown Stokes type */
    RHEA_ABORT_NOT_REACHED ();
  }
}

/******************************************************************************
 * Access Data of a Stokes Problem
 *****************************************************************************/

ymir_mesh_t *
rhea_stokes_problem_get_ymir_mesh (rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->ymir_mesh;
}

ymir_pressure_elem_t *
rhea_stokes_problem_get_press_elem (rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->press_elem;
}

ymir_vec_t *
rhea_stokes_problem_get_temperature (rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->temperature;
}

ymir_vec_t *
rhea_stokes_problem_get_weakzone (rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->weakzone;
}

ymir_vec_t *
rhea_stokes_problem_get_rhs_vel (rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->rhs_vel;
}

ymir_vec_t *
rhea_stokes_problem_get_rhs_vel_nonzero_dirichlet (
                                        rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->rhs_vel_nonzero_dirichlet;
}

void
rhea_stokes_problem_copy_viscosity (ymir_vec_t *viscosity,
                                    rhea_stokes_problem_t *stokes_problem)
{
  /* check input */
  RHEA_ASSERT (stokes_problem->coeff != NULL);

  /* copy Stokes coefficient and divide by 2 */
  ymir_vec_copy (stokes_problem->coeff, viscosity);
  ymir_vec_scale (0.5, viscosity);
}

ymir_stokes_op_t *
rhea_stokes_problem_get_stokes_op (rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->stokes_op;
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
