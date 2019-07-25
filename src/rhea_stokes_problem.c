#include <rhea_stokes_problem.h>
#include <rhea_base.h>
#include <rhea_strainrate.h>
#include <rhea_stress.h>
#include <rhea_stokes_norm.h>
#include <rhea_newton.h>
#include <rhea_stokes_problem_amr.h>
#include <rhea.h>
#include <ymir_mass_vec.h>
#include <ymir_stress_op_optimized.h>
#include <ymir_stokes_pc.h>
#include <ymir_vtk.h>
#include <p8est_extended.h>
#ifdef RHEA_ENABLE_DEBUG
# include <ymir_comm.h>
#endif

/* Stokes solver */
#define RHEA_STOKES_PROBLEM_NONLINEAR_GRAD_NORM_N_COMPONENTS 2

/* I/O labels that are attached at the end of file names */
#define RHEA_STOKES_PROBLEM_IO_LABEL_P4EST ".p4est"
#define RHEA_STOKES_PROBLEM_IO_LABEL_TEMPERATURE "_temp"
#define RHEA_STOKES_PROBLEM_IO_LABEL_VELOCITY "_vel"
#define RHEA_STOKES_PROBLEM_IO_LABEL_PRESSURE "_press"
#define RHEA_STOKES_PROBLEM_IO_LABEL_NL_ITER "_itn"

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_LINEARIZATION_TYPE \
  (RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_REGULAR)
#define RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_NORM_TYPE \
  (RHEA_STOKES_NORM_HMINUS1_L2)
#define RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_MASS_SCALING (0.0)
#define RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_CHECK_JACOBIAN (0)
#define RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_LIN_ANISO_CHECK (0)
#define RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_LIN_ANISO_PC (0)

/* initialize options */
int                 rhea_stokes_problem_nonlinear_linearization_type =
                      RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_LINEARIZATION_TYPE;
int                 rhea_stokes_problem_nonlinear_norm_type =
                      RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_NORM_TYPE;
double              rhea_stokes_problem_nonlinear_norm_mass_scaling =
                      RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_MASS_SCALING;
int                 rhea_stokes_problem_nonlinear_check_jacobian =
                      RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_CHECK_JACOBIAN;
int                 rhea_stokes_problem_nonlinear_lin_aniso_check =
                      RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_LIN_ANISO_CHECK;
int                 rhea_stokes_problem_nonlinear_lin_aniso_pc =
                      RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_LIN_ANISO_PC;

rhea_newton_options_t rhea_stokes_problem_newton_options;

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
  YMIR_OPTIONS_B, "nonlinear-check-jacobian", '\0',
    &(rhea_stokes_problem_nonlinear_check_jacobian),
    RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_CHECK_JACOBIAN,
    "Check Jacobians during the solve of a nonlinear problem with Newton",
  YMIR_OPTIONS_B, "nonlinear-linearization-anisotropy-check", '\0',
    &(rhea_stokes_problem_nonlinear_lin_aniso_check),
    RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_LIN_ANISO_CHECK,
    "Check linearization-induced anisotropy in gradient (=RHS) and Newton step",
  YMIR_OPTIONS_B, "nonlinear-linearization-anisotropy-pc", '\0',
    &(rhea_stokes_problem_nonlinear_lin_aniso_pc),
    RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_LIN_ANISO_PC,
    "Precondition linearization-induced anisotropy",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add sub-options */
  rhea_stokes_problem_amr_add_options (opt);
  rhea_newton_add_options (&rhea_stokes_problem_newton_options, opt);

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
  ymir_mesh_t          *ymir_mesh;
  ymir_pressure_elem_t *press_elem;
  int                   mesh_modified_by_solver;

  /* state variables of the Stokes problem (not owned) */
  ymir_vec_t         *temperature;
  ymir_vec_t         *velocity_pressure;

  /* options (not owned) */
  rhea_domain_options_t      *domain_options;
  rhea_temperature_options_t *temp_options;
  rhea_plate_options_t       *plate_options;
  rhea_weakzone_options_t    *weak_options;
  rhea_viscosity_options_t   *visc_options;

  /* callback functions for vector computations */
  rhea_weakzone_compute_fn_t  weakzone_compute_fn;
  void                       *weakzone_compute_fn_data;
  rhea_viscosity_compute_fn_t viscosity_compute_fn;
  void                       *viscosity_compute_fn_data;
  rhea_velocity_rhs_compute_fn_t         rhs_vel_compute_fn;
  void                                  *rhs_vel_compute_fn_data;
  rhea_velocity_rhs_nz_dir_compute_fn_t  rhs_vel_nonzero_dir_compute_fn;
  void                                  *rhs_vel_nonzero_dir_compute_fn_data;
  rhea_velocity_rhs_nz_neu_compute_fn_t  rhs_vel_nonzero_neu_compute_fn;
  void                                  *rhs_vel_nonzero_neu_compute_fn_data;

  /* boundary conditions (mesh data) */
  ymir_vel_dir_t     *vel_dir;

  /* coefficient (mesh data) */
  ymir_vec_t         *coeff;            /* out: linear & nonlinear Stokes */
  ymir_vec_t         *proj_scal;        /* out: nonlinear Stokes only */
  ymir_vec_t         *bounds_marker;    /* out: nonlinear Stokes only */
  ymir_vec_t         *yielding_marker;  /* out: nonlinear Stokes only */
  ymir_vec_t         *weakzone;         /* in: linear & nonlinear Stokes */
  ymir_vec_t         *sol_vel;          /* in: nonlinear Stokes only */
  ymir_vec_t         *proj_tens;        /* in: deprecated TODO */

  /* right-hand side (mesh data) */
  ymir_vec_t         *rhs_vel;
  ymir_vec_t         *rhs_vel_nonzero_dirichlet;
  ymir_vec_t        **rhs_vel_face_nonzero_neumann;
  ymir_vec_t         *rhs_vel_press;

  /* Stokes operator and preconditioner (solver data) */
  ymir_stokes_op_t   *stokes_op;
  ymir_stokes_pc_t   *stokes_pc;
  int                 num_krylov_solvers;

  /* Newton problem (nonlinear Stokes only) */
  rhea_stokes_problem_nonlinear_linearization_t linearization_type;
  rhea_newton_options_t  *newton_options;
  rhea_newton_problem_t  *newton_problem;
  ymir_vec_t             *newton_neg_gradient_vec;
  ymir_vec_t             *newton_step_vec;
  rhea_stokes_norm_type_t norm_type;
  ymir_Hminus1_norm_op_t *norm_op;
  double                  norm_op_mass_scaling;

  /* flag for clearing/creating mesh dependencies */
  int                 recreate_solver_data;

  /* data for solver AMR/grid continuation (not owned) */
  p4est_t            *p4est;
  rhea_discretization_options_t *discr_options;

  /* solver output paths */
  char               *solver_bin_path;
  char               *solver_vtk_path;
};

/**
 * Creates and initializes a new structure of a Stokes problem.
 */
static rhea_stokes_problem_t *
rhea_stokes_problem_struct_new (const rhea_stokes_problem_type_t type,
                                rhea_domain_options_t *domain_options,
                                rhea_temperature_options_t *temp_options,
                                rhea_plate_options_t *plate_options,
                                rhea_weakzone_options_t *weak_options,
                                rhea_viscosity_options_t *visc_options)
{
  rhea_stokes_problem_t *stokes_problem = RHEA_ALLOC (rhea_stokes_problem_t, 1);

  /* initialize */
  stokes_problem->type = type;
  stokes_problem->incompressible = 1;
  stokes_problem->ymir_mesh = NULL;
  stokes_problem->press_elem = NULL;
  stokes_problem->mesh_modified_by_solver = 0;

  stokes_problem->temperature = NULL;
  stokes_problem->velocity_pressure = NULL;

  stokes_problem->domain_options = domain_options;
  stokes_problem->temp_options = temp_options;
  stokes_problem->plate_options = plate_options;
  stokes_problem->weak_options = weak_options;
  stokes_problem->visc_options = visc_options;

  stokes_problem->vel_dir = NULL;

  stokes_problem->coeff = NULL;
  stokes_problem->proj_scal = NULL;
  stokes_problem->bounds_marker = NULL;
  stokes_problem->yielding_marker = NULL;
  stokes_problem->weakzone = NULL;
  stokes_problem->sol_vel = NULL;
  stokes_problem->proj_tens = NULL;

  stokes_problem->rhs_vel = NULL;
  stokes_problem->rhs_vel_nonzero_dirichlet = NULL;
  stokes_problem->rhs_vel_face_nonzero_neumann = NULL;
  stokes_problem->rhs_vel_press = NULL;

  stokes_problem->stokes_op = NULL;
  stokes_problem->stokes_pc = NULL;
  stokes_problem->num_krylov_solvers = 1;

  stokes_problem->linearization_type =
    RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_LINEARIZATION_TYPE;
  stokes_problem->newton_options = NULL;
  stokes_problem->newton_problem = NULL;
  stokes_problem->newton_neg_gradient_vec = NULL;
  stokes_problem->newton_step_vec = NULL;
  stokes_problem->norm_type = RHEA_STOKES_NORM_NONE;
  stokes_problem->norm_op = NULL;
  stokes_problem->norm_op_mass_scaling = NAN;

  stokes_problem->p4est = NULL;
  stokes_problem->discr_options = NULL;

  stokes_problem->recreate_solver_data = 0;

  stokes_problem->solver_bin_path = NULL;
  stokes_problem->solver_vtk_path = NULL;

  /* set default callback functions */
  if (rhea_weakzone_exists (weak_options)) {
    stokes_problem->weakzone_compute_fn = rhea_weakzone_compute;
    stokes_problem->weakzone_compute_fn_data = weak_options;
  }
  else {
    stokes_problem->weakzone_compute_fn = NULL;
    stokes_problem->weakzone_compute_fn_data = NULL;
  }
  rhea_stokes_problem_set_viscosity_compute_fn (
      stokes_problem, rhea_viscosity_compute, visc_options);
  stokes_problem->rhs_vel_compute_fn = rhea_velocity_rhs_compute;
  stokes_problem->rhs_vel_compute_fn_data = temp_options;
  stokes_problem->rhs_vel_nonzero_dir_compute_fn = NULL;
  stokes_problem->rhs_vel_nonzero_dir_compute_fn_data = NULL;
  stokes_problem->rhs_vel_nonzero_neu_compute_fn = NULL;
  stokes_problem->rhs_vel_nonzero_neu_compute_fn_data = NULL;

  return stokes_problem;
}

/**
 * Destroys a structure of a Stokes problem.
 */
static void
rhea_stokes_problem_struct_destroy (rhea_stokes_problem_t *stokes_problem)
{
  /* destroy structure */
  RHEA_FREE (stokes_problem);
}

/**
 * Creates new velocity boundary conditions.
 */
static void
rhea_stokes_problem_velocity_boundary_create (
                                        rhea_stokes_problem_t *stokes_problem)
{
  ymir_mesh_t        *ymir_mesh = stokes_problem->ymir_mesh;
  rhea_domain_options_t *domain_options = stokes_problem->domain_options;
  ymir_vec_t         *dirscal = NULL; //TODO seems to be best choice, isn't it?

  /* check input */
  RHEA_ASSERT (stokes_problem->ymir_mesh != NULL);
  RHEA_ASSERT (stokes_problem->domain_options != NULL);
  RHEA_ASSERT (stokes_problem->vel_dir == NULL);

  /* create Dirichlet boundary conditions */
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
  ymir_vel_dir_t     *vel_dir = stokes_problem->vel_dir;

  /* destroy Dirichlet boundary conditions */
  if (vel_dir != NULL) {
    if (vel_dir->scale != NULL) {
      ymir_vec_destroy (vel_dir->scale);
    }
    ymir_vel_dir_destroy (vel_dir);
    stokes_problem->vel_dir = NULL;
  }
}

/**
 * Retrieves velocity.
 */
static ymir_vec_t *
rhea_stokes_problem_retrieve_velocity (ymir_vec_t *velocity_pressure,
                                       rhea_stokes_problem_t *stokes_problem)
{
  ymir_pressure_elem_t *press_elem = stokes_problem->press_elem;
  ymir_vec_t         *vel = stokes_problem->sol_vel;

  /* check input */
  RHEA_ASSERT (stokes_problem->press_elem != NULL);
  RHEA_ASSERT (stokes_problem->sol_vel != NULL);
  RHEA_ASSERT (rhea_velocity_check_vec_type (stokes_problem->sol_vel));
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (velocity_pressure));

  /* retrieve velocity */
  if (velocity_pressure != NULL) { /* if velocity is given */
    RHEA_ASSERT (rhea_velocity_pressure_is_valid (velocity_pressure));
    rhea_velocity_pressure_copy_components (vel, NULL, velocity_pressure,
                                            press_elem);
    RHEA_ASSERT (rhea_velocity_is_valid (vel));
    rhea_stokes_problem_velocity_boundary_set_zero (vel, stokes_problem);
    RHEA_ASSERT (rhea_velocity_is_valid (vel));
  }
  else { /* otherwise assume zero velocity */
    ymir_vec_set_zero (vel);
  }

  /* return velocity */
  return vel;
}

void
rhea_stokes_problem_compute_coefficient (rhea_stokes_problem_t *stokes_problem,
                                         ymir_vec_t *velocity_pressure,
                                         const int nonlinear_init)
{
  rhea_viscosity_options_t *visc_options = stokes_problem->visc_options;
  /* input vectors */
  ymir_vec_t         *temperature     = stokes_problem->temperature;
  ymir_vec_t         *weakzone        = stokes_problem->weakzone;
  ymir_vec_t         *vel             = NULL;
  /* output vectors */
  ymir_vec_t         *coeff           = stokes_problem->coeff;
  ymir_vec_t         *proj_scal       = stokes_problem->proj_scal;
  ymir_vec_t         *bounds_marker   = stokes_problem->bounds_marker;
  ymir_vec_t         *yielding_marker = stokes_problem->yielding_marker;

  /* check input */
  RHEA_ASSERT (stokes_problem->coeff != NULL);
  RHEA_ASSERT (stokes_problem->visc_options != NULL);

  /* compute viscosity */
  switch (stokes_problem->type) {
  case RHEA_STOKES_PROBLEM_LINEAR:
    RHEA_ASSERT (visc_options->type == RHEA_VISCOSITY_LINEAR);
    rhea_stokes_problem_viscosity_compute (
        /* out: */ coeff, NULL, NULL, NULL,
        /* in:  */ temperature, weakzone, NULL, stokes_problem);
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR:
    RHEA_ASSERT (visc_options->type == RHEA_VISCOSITY_NONLINEAR);
    RHEA_ASSERT (stokes_problem->bounds_marker != NULL);
    RHEA_ASSERT (stokes_problem->yielding_marker != NULL);
    if (nonlinear_init) {
      rhea_viscosity_compute_nonlinear_init (
          /* out: */ coeff, proj_scal, bounds_marker, yielding_marker,
          /* in:  */ temperature, weakzone, visc_options);
    }
    else {
      RHEA_ASSERT (velocity_pressure != NULL);
      RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (velocity_pressure));
      RHEA_ASSERT (rhea_velocity_pressure_is_valid (velocity_pressure));
      rhea_stokes_problem_set_velocity_pressure (
          stokes_problem, velocity_pressure);
      vel = rhea_stokes_problem_retrieve_velocity (
          stokes_problem->velocity_pressure, stokes_problem);
      rhea_stokes_problem_viscosity_compute (
          /* out: */ coeff, proj_scal, bounds_marker, yielding_marker,
          /* in:  */ temperature, weakzone, vel, stokes_problem);
    }
    break;
  default: /* unknown Stokes type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* transform viscosity to viscous stress coefficient */
  ymir_vec_scale (2.0, coeff);
}

void
rhea_stokes_problem_compute_and_update_coefficient (
                                        rhea_stokes_problem_t *stokes_problem,
                                        ymir_vec_t *velocity_pressure,
                                        const int nonlinear_init)
{
  /* check input */
  RHEA_ASSERT (stokes_problem->stokes_op != NULL);
  RHEA_ASSERT (stokes_problem->visc_options != NULL);
  RHEA_ASSERT (stokes_problem->coeff != NULL);

  /* compute coefficient */
  rhea_stokes_problem_compute_coefficient (stokes_problem, velocity_pressure,
                                           nonlinear_init);

  /* set viscous stress coefficient */
  {
    ymir_stress_op_t   *stress_op = stokes_problem->stokes_op->stress_op;
    rhea_viscosity_options_t *visc_options = stokes_problem->visc_options;

    RHEA_ASSERT (stokes_problem->type == RHEA_STOKES_PROBLEM_LINEAR ||
                 ymir_stress_op_has_linearized (stress_op));
    RHEA_ASSERT (ymir_stress_op_get_coeff_type (stress_op) ==
                 YMIR_STRESS_OP_COEFF_ISOTROPIC_SCAL);
    RHEA_ASSERT (fabs (ymir_stress_op_get_coeff_shift (stress_op) -
                       rhea_viscosity_get_visc_shift (visc_options)) < SC_EPS);

    ymir_stress_op_set_coeff_scal (stress_op, stokes_problem->coeff);
  }
}

/**
 * Calculates the reduction of an end value relative to a start value.  Avoids
 * ill-defined operations like division by zero.
 */
static double
rhea_stokes_problem_calculate_reduction (const double start_value,
                                         const double end_value)
{
  /* check input */
  RHEA_ASSERT (isfinite (start_value));
  RHEA_ASSERT (0.0 <= start_value);
  RHEA_ASSERT (isfinite (end_value));
  RHEA_ASSERT (0.0 <= end_value);

  if (end_value < DBL_MIN) { /* if end value is small */
    /* assume max reduction and set value below machine precision */
    return pow (SC_EPS, 2.0);
  }
  else if (start_value < SC_EPS * end_value) { /* if start value is small */
    /* assume no reduction */
    return 1.0;
  }
  else { /* if division is well-defined */
    /* compute reduction */
    return end_value / start_value;
  }
}

/**
 * Runs the Krylov solver to solve a Stokes system.
 */
static int
rhea_stokes_problem_run_krylov_solver (ymir_vec_t *sol_vel_press,
                                       ymir_vec_t *rhs_vel_press,
                                       const int nonzero_initial_guess,
                                       const int iter_max,
                                       const double rtol,
                                       rhea_stokes_problem_t *stokes_problem,
                                       const int krylov_solver_idx,
                                       int *num_iterations,
                                       double residual_reduction[3])
{
  const double        atol = pow (SC_EPS, 2);
  const int           check_residual = (NULL != residual_reduction);
  ymir_vec_t         *residual_vel_press = NULL;
  int                 stop_reason, itn;
  double              norm_res_init = NAN, norm_res = NAN;
  double              norm_res_init_vel = NAN, norm_res_vel = NAN;
  double              norm_res_init_press = NAN, norm_res_press = NAN;
  double              res_reduction;
  double              res_reduction_vel, res_reduction_press;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (
      __func__, "iter_max=%i, rtol=%.1e, nonzero_init_guess=%i, krylov_idx=%i",
      iter_max, rtol, nonzero_initial_guess, krylov_solver_idx);

  /* check input */
  RHEA_ASSERT (stokes_problem->stokes_op != NULL);
  RHEA_ASSERT (stokes_problem->stokes_pc != NULL);
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (sol_vel_press));

  /* compute initial norm of (negative) residual (A*u - f) */
  if (check_residual) {
    residual_vel_press = ymir_vec_template (rhs_vel_press);

    if (nonzero_initial_guess) {
      RHEA_ASSERT (residual_vel_press != NULL);
      ymir_stokes_pc_apply_stokes_op (sol_vel_press, residual_vel_press,
                                      stokes_problem->stokes_op, 0, 1);
      ymir_vec_add (-1.0, rhs_vel_press, residual_vel_press);
      norm_res_init = rhea_stokes_norm_compute (
          &norm_res_init_vel, &norm_res_init_press, residual_vel_press,
          RHEA_STOKES_NORM_L2_VEC_SP, NULL, stokes_problem->press_elem);
    }
    else { /* if zero inital guess */
      norm_res_init = rhea_stokes_norm_compute (
          &norm_res_init_vel, &norm_res_init_press, rhs_vel_press,
          RHEA_STOKES_NORM_L2_VEC_SP, NULL, stokes_problem->press_elem);
    }
  }

  /* run solver for Stokes */
  stop_reason = ymir_stokes_pc_solve (rhs_vel_press, sol_vel_press,
                                      stokes_problem->stokes_pc,
                                      nonzero_initial_guess, rtol, atol,
                                      iter_max, krylov_solver_idx, &itn);
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (sol_vel_press));

  /* compute norm of (negative) residual (A*u - f) at solution */
  if (check_residual) {
    RHEA_ASSERT (residual_vel_press != NULL);
    ymir_stokes_pc_apply_stokes_op (sol_vel_press, residual_vel_press,
                                    stokes_problem->stokes_op, 0, 1);
    ymir_vec_add (-1.0, rhs_vel_press, residual_vel_press);
    norm_res = rhea_stokes_norm_compute (
        &norm_res_vel, &norm_res_press, residual_vel_press,
        RHEA_STOKES_NORM_L2_VEC_SP, NULL, stokes_problem->press_elem);

    /* calculate convergence (note: division by zero is possible) */
    res_reduction =
      rhea_stokes_problem_calculate_reduction (norm_res_init, norm_res);
    res_reduction_vel =
      rhea_stokes_problem_calculate_reduction (norm_res_init_vel, norm_res_vel);
    if (stokes_problem->incompressible) {
      res_reduction_press = norm_res_press;
    }
    else {
      res_reduction_press =
        rhea_stokes_problem_calculate_reduction (norm_res_init_press,
                                                 norm_res_press);
    }

    ymir_vec_destroy (residual_vel_press);
  }

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);

  /* return iterations count, residual reduction, and stopping reason */
  if (NULL != num_iterations) {
    *num_iterations = itn;
  }
  if (NULL != residual_reduction) {
    residual_reduction[0] = res_reduction;
    residual_reduction[1] = res_reduction_vel;
    residual_reduction[2] = res_reduction_press;
  }
  return stop_reason;
}

/******************************************************************************
 * Linear Stokes Problem
 *****************************************************************************/

static int
rhea_stokes_problem_linear_mesh_data_exists (
                                    rhea_stokes_problem_t *stokes_problem_lin)
{
  int                 mesh_exists, vectors_exist;
  int                 weak_exists;
  int                 rhs_nonzero_dir_exists, rhs_nonzero_neu_exists;

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);

  /* check if mesh and corresponding data exist */
  mesh_exists = (
      stokes_problem_lin->ymir_mesh != NULL &&
      stokes_problem_lin->press_elem != NULL);
  vectors_exist = (
      stokes_problem_lin->coeff != NULL &&
      stokes_problem_lin->rhs_vel != NULL &&
      stokes_problem_lin->rhs_vel_press != NULL);

  /* check if not required vectors exist */
  weak_exists = (
      stokes_problem_lin->weakzone_compute_fn == NULL ||
      stokes_problem_lin->weakzone != NULL);
  rhs_nonzero_dir_exists = (
      stokes_problem_lin->rhs_vel_nonzero_dir_compute_fn == NULL ||
      stokes_problem_lin->rhs_vel_nonzero_dirichlet != NULL);
  rhs_nonzero_neu_exists = (
      stokes_problem_lin->rhs_vel_nonzero_neu_compute_fn == NULL ||
      stokes_problem_lin->rhs_vel_face_nonzero_neumann != NULL);

  return (mesh_exists && vectors_exist && weak_exists &&
          rhs_nonzero_dir_exists && rhs_nonzero_neu_exists);
}

static void
rhea_stokes_problem_linear_create_mesh_data (
                                    rhea_stokes_problem_t *stokes_problem_lin,
                                    ymir_mesh_t *ymir_mesh,
                                    ymir_pressure_elem_t *press_elem)
{
  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);
  RHEA_ASSERT (
      !rhea_stokes_problem_linear_mesh_data_exists (stokes_problem_lin));
  RHEA_ASSERT (stokes_problem_lin->domain_options != NULL);

  /* assign mesh objects */
  stokes_problem_lin->ymir_mesh = ymir_mesh;
  stokes_problem_lin->press_elem = press_elem;

  /* create velocity boundary conditions */
  rhea_stokes_problem_velocity_boundary_create (stokes_problem_lin);

  /* create coefficient and vectors related to it */
  if (stokes_problem_lin->weakzone_compute_fn != NULL) {
    stokes_problem_lin->weakzone = rhea_weakzone_new (ymir_mesh);
  }
  stokes_problem_lin->coeff = rhea_viscosity_new (ymir_mesh);

  /* create right-hand side vectors */
  stokes_problem_lin->rhs_vel = rhea_velocity_new (ymir_mesh);
  stokes_problem_lin->rhs_vel_compute_fn (
      stokes_problem_lin->rhs_vel,
      stokes_problem_lin->temperature,
      stokes_problem_lin->rhs_vel_compute_fn_data);
  if (stokes_problem_lin->rhs_vel_nonzero_dir_compute_fn != NULL) {
    stokes_problem_lin->rhs_vel_nonzero_dirichlet = rhea_velocity_new (
        ymir_mesh);
    stokes_problem_lin->rhs_vel_nonzero_dir_compute_fn (
        stokes_problem_lin->rhs_vel_nonzero_dirichlet,
        stokes_problem_lin->rhs_vel_nonzero_dir_compute_fn_data);
  }
  if (stokes_problem_lin->rhs_vel_nonzero_neu_compute_fn != NULL) {
    ymir_topidx_t       fm;

    stokes_problem_lin->rhs_vel_face_nonzero_neumann =
      RHEA_ALLOC (ymir_vec_t *, ymir_mesh->num_face_meshes);
    for (fm = 0; fm < ymir_mesh->num_face_meshes; fm++) { /* for all faces */
      stokes_problem_lin->rhs_vel_face_nonzero_neumann[fm] =
        ymir_face_cvec_new (ymir_mesh, fm, 3);
    }
    stokes_problem_lin->rhs_vel_nonzero_neu_compute_fn (
        stokes_problem_lin->rhs_vel_face_nonzero_neumann,
        stokes_problem_lin->rhs_vel_nonzero_neu_compute_fn_data);
  }
  stokes_problem_lin->rhs_vel_press = rhea_velocity_pressure_new (
      ymir_mesh, press_elem);

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static void
rhea_stokes_problem_linear_clear_mesh_data (
                                    rhea_stokes_problem_t *stokes_problem_lin)
{
  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);
  RHEA_ASSERT (
      rhea_stokes_problem_linear_mesh_data_exists (stokes_problem_lin));

  /* destroy vectors */
  rhea_viscosity_destroy (stokes_problem_lin->coeff);
  rhea_velocity_destroy (stokes_problem_lin->rhs_vel);
  rhea_velocity_pressure_destroy (stokes_problem_lin->rhs_vel_press);
  stokes_problem_lin->coeff = NULL;
  stokes_problem_lin->rhs_vel = NULL;
  stokes_problem_lin->rhs_vel_press = NULL;
  if (stokes_problem_lin->weakzone != NULL) {
    rhea_weakzone_destroy (stokes_problem_lin->weakzone);
    stokes_problem_lin->weakzone = NULL;
  }
  if (stokes_problem_lin->rhs_vel_nonzero_dirichlet != NULL) {
    rhea_velocity_destroy (stokes_problem_lin->rhs_vel_nonzero_dirichlet);
    stokes_problem_lin->rhs_vel_nonzero_dirichlet = NULL;
  }
  if (stokes_problem_lin->rhs_vel_face_nonzero_neumann != NULL) {
    ymir_mesh_t        *ymir_mesh = stokes_problem_lin->ymir_mesh;
    ymir_topidx_t       fm;

    for (fm = 0; fm < ymir_mesh->num_face_meshes; fm++) { /* for all faces */
      ymir_vec_destroy (stokes_problem_lin->rhs_vel_face_nonzero_neumann[fm]);
    }
    RHEA_FREE (stokes_problem_lin->rhs_vel_face_nonzero_neumann);
    stokes_problem_lin->rhs_vel_face_nonzero_neumann = NULL;
  }

  /* destroy velocity boundary conditions */
  rhea_stokes_problem_velocity_boundary_clear (stokes_problem_lin);

  /* remove vector references */
  rhea_stokes_problem_remove_temperature (stokes_problem_lin);
  rhea_stokes_problem_remove_velocity_pressure (stokes_problem_lin);

  /* remove mesh references */
  stokes_problem_lin->ymir_mesh = NULL;
  stokes_problem_lin->press_elem = NULL;

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static int
rhea_stokes_problem_linear_solver_data_exists (
                                    rhea_stokes_problem_t *stokes_problem_lin)
{
  int                 mesh_data_exists, stokes_exists;

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);

  /* check if mesh data exists */
  mesh_data_exists =
    rhea_stokes_problem_linear_mesh_data_exists (stokes_problem_lin);

  /* check if Stokes operator exists */
  stokes_exists = (stokes_problem_lin->stokes_op != NULL);

  return (mesh_data_exists && stokes_exists);
}

static void
rhea_stokes_problem_linear_create_solver_data (
                                    rhea_stokes_problem_t *stokes_problem_lin)
{
  const char         *vtk_path = stokes_problem_lin->solver_vtk_path;
  rhea_domain_options_t *domain_options = stokes_problem_lin->domain_options;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);
  RHEA_ASSERT (
      rhea_stokes_problem_linear_mesh_data_exists (stokes_problem_lin));
  RHEA_ASSERT (
      !rhea_stokes_problem_linear_solver_data_exists (stokes_problem_lin));
  RHEA_ASSERT (stokes_problem_lin->domain_options != NULL);

  /* set VTK path for debuggin */
  if (vtk_path != NULL) {
    ymir_vtk_set_debug_path (vtk_path);
  }

  /* compute weak zone */
  rhea_stokes_problem_weakzone_compute (stokes_problem_lin->weakzone,
                                        stokes_problem_lin);

  /* compute viscosity */
  rhea_stokes_problem_compute_coefficient (stokes_problem_lin,
                                           NULL /* unused */, 0 /* unused */);

  /* create Stokes operator */
  stokes_problem_lin->stokes_op = ymir_stokes_op_new_ext (
      stokes_problem_lin->coeff, stokes_problem_lin->vel_dir,
      NULL /* Robin BC's */,
      NULL /* deprecated */,
      NULL /* deprecated */,
      stokes_problem_lin->press_elem,
      domain_options->center,
      domain_options->moment_of_inertia);

  /* set viscous stress coefficient */
  {
    ymir_stress_op_t   *stress_op = stokes_problem_lin->stokes_op->stress_op;
    rhea_viscosity_options_t *visc_options = stokes_problem_lin->visc_options;

    ymir_stress_op_set_coeff_type (
        stress_op, YMIR_STRESS_OP_COEFF_ISOTROPIC_SCAL);
    ymir_stress_op_set_coeff_shift (
        stress_op, rhea_viscosity_get_visc_shift (visc_options));
  }

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static void
rhea_stokes_problem_linear_clear_solver_data (
                                    rhea_stokes_problem_t *stokes_problem_lin)
{
  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);

  /* destroy Stokes preconditioner */
  if (stokes_problem_lin->stokes_pc != NULL) {
    ymir_stokes_pc_destroy (stokes_problem_lin->stokes_pc);
    stokes_problem_lin->stokes_pc = NULL;
  }

  /* destroy Stokes operator */
  if (stokes_problem_lin->stokes_op != NULL) {
    ymir_stokes_op_destroy (stokes_problem_lin->stokes_op);
    stokes_problem_lin->stokes_op = NULL;
  }

  /* clear VTK path for debugging */
  ymir_vtk_set_debug_path ("");

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static rhea_stokes_problem_t *
rhea_stokes_problem_linear_new (ymir_mesh_t *ymir_mesh,
                                ymir_pressure_elem_t *press_elem,
                                ymir_vec_t *temperature,
                                rhea_domain_options_t *domain_options,
                                rhea_temperature_options_t *temp_options,
                                rhea_plate_options_t *plate_options,
                                rhea_weakzone_options_t *weak_options,
                                rhea_viscosity_options_t *visc_options)
{
  rhea_stokes_problem_t *stokes_problem_lin;

  RHEA_GLOBAL_PRODUCTION_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (visc_options->type == RHEA_VISCOSITY_LINEAR);

  /* initialize general Stokes problem structure */
  stokes_problem_lin = rhea_stokes_problem_struct_new (
      RHEA_STOKES_PROBLEM_LINEAR,
      domain_options, temp_options, plate_options, weak_options, visc_options);

  /* set temperature */
  rhea_stokes_problem_set_temperature (stokes_problem_lin, temperature);

  /* create mesh data */
  rhea_stokes_problem_linear_create_mesh_data (stokes_problem_lin, ymir_mesh,
                                               press_elem);

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);

  return stokes_problem_lin;
}

static void
rhea_stokes_problem_linear_destroy (rhea_stokes_problem_t *stokes_problem_lin)
{
  RHEA_GLOBAL_PRODUCTION_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);

  /* destroy solver data */
  rhea_stokes_problem_linear_clear_solver_data (stokes_problem_lin);

  /* destroy mesh data */
  rhea_stokes_problem_linear_clear_mesh_data (stokes_problem_lin);

  /* destroy problem structure */
  rhea_stokes_problem_struct_destroy (stokes_problem_lin);

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);
}

static void
rhea_stokes_problem_linear_setup_solver (
                                    rhea_stokes_problem_t *stokes_problem_lin)
{
  ymir_vec_t         *rhs_vel;
  ymir_vec_t         *rhs_vel_nonzero_dirichlet;
  ymir_vec_t        **rhs_vel_face_nonzero_neumann;
  ymir_vec_t         *rhs_vel_press;
  ymir_stokes_op_t   *stokes_op;

  RHEA_GLOBAL_PRODUCTION_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);
  RHEA_ASSERT (
      rhea_stokes_problem_linear_mesh_data_exists (stokes_problem_lin));
  RHEA_ASSERT (stokes_problem_lin->stokes_pc == NULL);

  /* create solver data */
  rhea_stokes_problem_linear_create_solver_data (stokes_problem_lin);
  RHEA_ASSERT (stokes_problem_lin->stokes_op != NULL);

  /* construct right-hand side for Stokes system */
  rhs_vel = stokes_problem_lin->rhs_vel;
  rhs_vel_nonzero_dirichlet = stokes_problem_lin->rhs_vel_nonzero_dirichlet;
  rhs_vel_face_nonzero_neumann =
    stokes_problem_lin->rhs_vel_face_nonzero_neumann;
  rhs_vel_press = stokes_problem_lin->rhs_vel_press;
  stokes_op = stokes_problem_lin->stokes_op;
  ymir_stokes_pc_construct_rhs (
      rhs_vel_press,                /* output: right-hand side */
      rhs_vel,                      /* input: volume forcing */
      rhs_vel_face_nonzero_neumann, /* input: nonzero Neumann forcing */
      rhs_vel_nonzero_dirichlet,    /* input: nonzero Dirichlet boundary */
      stokes_problem_lin->incompressible, stokes_op, 0 /* !linearized */);
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (rhs_vel_press));

  /* build Stokes preconditioner */
  stokes_problem_lin->stokes_pc = ymir_stokes_pc_new_multi_solver (
      stokes_op, 0 /* !linearized */, stokes_problem_lin->num_krylov_solvers);

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);
}

static void
rhea_stokes_problem_linear_update_solver (
                                    rhea_stokes_problem_t *stokes_problem_lin,
                                    const int update_coeff,
                                    const int override_rhs,
                                    ymir_vec_t *rhs_vel_press,
                                    ymir_vec_t *rhs_vel,
                                    ymir_vec_t *rhs_vel_nonzero_dirichlet,
                                    ymir_vec_t **rhs_vel_face_nonzero_neumann)
{
  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (__func__, "update_coeff=%i, override_rhs=%i",
                                 update_coeff, override_rhs);

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);
  RHEA_ASSERT (
      rhea_stokes_problem_linear_mesh_data_exists (stokes_problem_lin));
  RHEA_ASSERT (
      rhea_stokes_problem_linear_solver_data_exists (stokes_problem_lin));
  RHEA_ASSERT (stokes_problem_lin->stokes_pc != NULL);
  RHEA_ASSERT (rhs_vel_press == NULL ||
               rhea_velocity_pressure_check_vec_type (rhs_vel_press));
  RHEA_ASSERT (rhs_vel == NULL || rhea_velocity_check_vec_type (rhs_vel));
  RHEA_ASSERT (rhs_vel_nonzero_dirichlet == NULL ||
               rhea_velocity_check_vec_type (rhs_vel_nonzero_dirichlet));
#ifdef RHEA_ENABLE_DEBUG
  if (rhs_vel_face_nonzero_neumann != NULL) {
    ymir_mesh_t        *ymir_mesh = stokes_problem_lin->ymir_mesh;
    ymir_topidx_t       fm;

    for (fm = 0; fm < ymir_mesh->num_face_meshes; fm++) {
      RHEA_ASSERT (
          rhea_velocity_check_vec_type (rhs_vel_face_nonzero_neumann[fm]) &&
          ymir_vec_is_face_vec (rhs_vel_face_nonzero_neumann[fm]) &&
          rhs_vel_face_nonzero_neumann[fm]->meshnum == fm );
    }
  }
#endif

  /* update the Stokes coefficient and dependent objects */
  if (update_coeff) {
    /* compute weak zone */
    rhea_stokes_problem_weakzone_compute (stokes_problem_lin->weakzone,
                                          stokes_problem_lin);

    /* update coefficient of Stokes operator */
    rhea_stokes_problem_compute_and_update_coefficient (
        stokes_problem_lin, NULL /* vel_press */, 0 /* !init */);

    /* update Stokes preconditioner */
    ymir_stokes_pc_recompute (stokes_problem_lin->stokes_pc);
  }

  /* re-construct right-hand side for Stokes system */
  if (override_rhs) {
    RHEA_GLOBAL_VERBOSEF_FN_TAG (
        __func__,
        "rhs_vel_press=%i, rhs_vel=%i, rhs_vel_nz_dir=%i, rhs_vel_nz_neu=%i",
        (rhs_vel_press != NULL), (rhs_vel != NULL),
        (rhs_vel_nonzero_dirichlet != NULL),
        (rhs_vel_face_nonzero_neumann != NULL));
    if (rhs_vel_press != NULL) {
      ymir_vec_copy (rhs_vel_press, stokes_problem_lin->rhs_vel_press);
    }
    else {
      ymir_stokes_pc_construct_rhs (
          stokes_problem_lin->rhs_vel_press, /* output: right-hand side */
          rhs_vel,                      /* input: volume forcing */
          rhs_vel_face_nonzero_neumann, /* input: nonzero Neumann forcing */
          rhs_vel_nonzero_dirichlet,    /* input: nonzero Dirichlet boundary */
          stokes_problem_lin->incompressible,
          stokes_problem_lin->stokes_op, 0 /* !linearized */);
    }
  }
  else {
    rhs_vel = stokes_problem_lin->rhs_vel;
    rhs_vel_nonzero_dirichlet =
      stokes_problem_lin->rhs_vel_nonzero_dirichlet;
    rhs_vel_face_nonzero_neumann =
      stokes_problem_lin->rhs_vel_face_nonzero_neumann;
    ymir_stokes_pc_construct_rhs (
        stokes_problem_lin->rhs_vel_press, /* output: right-hand side */
        rhs_vel,                      /* input: volume forcing */
        rhs_vel_face_nonzero_neumann, /* input: nonzero Neumann forcing */
        rhs_vel_nonzero_dirichlet,    /* input: nonzero Dirichlet boundary */
        stokes_problem_lin->incompressible,
        stokes_problem_lin->stokes_op, 0 /* !linearized */);
  }
  RHEA_ASSERT (
      rhea_velocity_pressure_is_valid (stokes_problem_lin->rhs_vel_press));

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static int
rhea_stokes_problem_linear_solve (ymir_vec_t **sol_vel_press,
                                  const int nonzero_initial_guess,
                                  const int iter_max,
                                  const double rtol,
                                  rhea_stokes_problem_t *stokes_problem_lin,
                                  const int krylov_solver_idx,
                                  int *num_iterations,
                                  double *residual_reduction)
{
  const int           check_residual = (NULL != residual_reduction) ||
                                       !rhea_production_run_get ();
  ymir_vec_t         *rhs_vel_press = stokes_problem_lin->rhs_vel_press;
  int                 stop_reason, itn;
  double              res_reduction[3], conv_factor;
  char                func_name_stop[BUFSIZ];

  /* setup solver if it has not been done yet */
  if (stokes_problem_lin->stokes_pc == NULL) {
    rhea_stokes_problem_linear_setup_solver (stokes_problem_lin);
  }

  RHEA_GLOBAL_PRODUCTIONF_FN_BEGIN (
      __func__, "iter_max=%i, rtol=%.1e, nonzero_init_guess=%i",
      iter_max, rtol, nonzero_initial_guess);

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);
  RHEA_ASSERT (stokes_problem_lin->rhs_vel_press != NULL);
  RHEA_ASSERT (stokes_problem_lin->stokes_op != NULL);
  RHEA_ASSERT (stokes_problem_lin->stokes_pc != NULL);
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (*sol_vel_press));

  /* run Krylov solver */
  stop_reason = rhea_stokes_problem_run_krylov_solver (
      *sol_vel_press, rhs_vel_press, nonzero_initial_guess, iter_max, rtol,
      stokes_problem_lin, krylov_solver_idx,
      &itn, (check_residual ?  res_reduction : NULL));

  /* project out nullspaces of solution */
//TODO add option
//rhea_stokes_problem_project_out_nullspace (*sol_vel_press, stokes_problem_nl);

  /* print status */
  snprintf (func_name_stop, BUFSIZ, "%s_status", __func__);
  if (check_residual) {
    char                mass_label[16];

    conv_factor = exp (log (res_reduction[0]) / ((double) itn));
    if (stokes_problem_lin->incompressible) {
      snprintf (mass_label, 16, "%s", "norm_mass");
    }
    else {
      snprintf (mass_label, 16, "%s", "reduction_mass");
    }

    RHEA_GLOBAL_PRODUCTIONF_FN_TAG (
        func_name_stop, "reason=%i, iterations=%i, "
        "residual_reduction=%.3e, convergence=%.3e, "
        "reduction_momentum=%.3e, %s=%.3e",
        stop_reason, itn, res_reduction[0], conv_factor,
        res_reduction[1], mass_label, res_reduction[2]);
  }
  else {
    RHEA_GLOBAL_PRODUCTIONF_FN_TAG (
        func_name_stop, "<reason=%i, iterations=%i", stop_reason, itn);
  }
  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);

  /* return iterations count, residual reduction, and stopping reason */
  if (NULL != num_iterations) {
    *num_iterations = itn;
  }
  if (NULL != residual_reduction) {
    *residual_reduction = res_reduction[0];
  }
  return stop_reason;
}

/******************************************************************************
 * Callback Functions for Newton Solver Used by Nonlinear Stokes Problem
 *****************************************************************************/

/**
 * Updates the nonlinear operator at the current solution vector.
 * (Callback function for Newton's method)
 */
static void
rhea_stokes_problem_nonlinear_update_operator_fn (ymir_vec_t *solution,
                                                  void *data)
{
  const int           sol_exists = (solution != NULL);
  rhea_stokes_problem_t *stokes_problem_nl = data;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->stokes_op != NULL);
  RHEA_ASSERT (!sol_exists || rhea_velocity_pressure_check_vec_type (solution));
  RHEA_ASSERT (!sol_exists || rhea_velocity_pressure_is_valid (solution));

  /* compute coefficient */
  rhea_stokes_problem_compute_and_update_coefficient (
      stokes_problem_nl, solution, !sol_exists /* is init? */);

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

/**
 * Updates the Hessian operator at the current solution vector.
 * (Callback function for Newton's method)
 */
static void
rhea_stokes_problem_nonlinear_update_hessian_fn (ymir_vec_t *solution,
                                                 ymir_vec_t *step_vec,
                                                 const double step_length,
                                                 void *data)
{
  rhea_stokes_problem_t *stokes_problem_nl = data;
  rhea_stokes_problem_nonlinear_linearization_t linearization_type =
    stokes_problem_nl->linearization_type;
  rhea_viscosity_options_t *visc_options = stokes_problem_nl->visc_options;
  ymir_mesh_t        *ymir_mesh = stokes_problem_nl->ymir_mesh;
  const int           solution_exists = (solution != NULL);
  const int           step_exists =
    (step_vec != NULL && isfinite (step_length) && 0.0 < step_length);
  ymir_stokes_op_t   *stokes_op;
  ymir_stress_op_t   *stress_op;
  ymir_vec_t         *coeff;
  ymir_vec_t         *sol_vel;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (
      __func__, "linearization=%i, sol=%i, step=%i",
      (int) linearization_type, solution_exists, step_exists);

  /* check input */
#ifdef RHEA_ENABLE_DEBUG
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->stokes_op != NULL);
  if (solution_exists) {
    RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (solution));
    RHEA_ASSERT (rhea_velocity_pressure_is_valid (solution));
  }
  if (step_exists) {
    RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (step_vec));
    RHEA_ASSERT (rhea_velocity_pressure_is_valid (step_vec));
  }
#endif

  /* get data from nonlinear Stokes problem object */
  stokes_op = stokes_problem_nl->stokes_op;
  coeff = stokes_problem_nl->coeff;

  /* get viscous stress operator */
  stress_op = stokes_op->stress_op;
  RHEA_ASSERT (ymir_stress_op_has_linearized (stress_op));
  RHEA_ASSERT (
      fabs (ymir_stress_op_get_coeff_shift_proj (stress_op) -
            rhea_viscosity_get_visc_shift_proj (visc_options)) < SC_EPS );

  /*
   * Set Linearized Part of the Viscous Stress Coefficient
   */

  switch (linearization_type) {
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_PICARD:
    RHEA_ASSERT (ymir_stress_op_get_coeff_type_linearized (stress_op) ==
                 YMIR_STRESS_OP_COEFF_ISOTROPIC_SCAL);
    break; /* RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_PICARD */

  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_REGULAR:
    {
      ymir_vec_t         *proj_scal = stokes_problem_nl->proj_scal;
      ymir_vec_t         *proj_tens = stokes_problem_nl->proj_tens;

      /* check input */
      RHEA_ASSERT (ymir_stress_op_get_coeff_type_linearized (stress_op) ==
                   YMIR_STRESS_OP_COEFF_ANISOTROPIC_PPR1);
      RHEA_ASSERT (stokes_problem_nl->proj_scal != NULL);
      RHEA_ASSERT (stokes_problem_nl->proj_tens != NULL);
      RHEA_ASSERT (rhea_viscosity_check_vec_type (proj_scal));
      RHEA_ASSERT (rhea_strainrate_check_vec_type (proj_tens));

      /* set coefficient data */
      if (solution_exists) { /* if solution is provided */
        /* retrieve velocity */
        sol_vel = rhea_stokes_problem_retrieve_velocity (
            solution, stokes_problem_nl);

        /* compute and normalize strain rate tensor */
        ymir_stress_op_optimized_compute_strain_rate (
            proj_tens, sol_vel, stress_op);
        ymir_stress_op_tensor_normalize (proj_tens);
      }
      else { /* if solution is not provided */
        ymir_vec_set_zero (proj_tens);
      }

      /* pass coefficient data to viscous stress operator */
      ymir_stress_op_set_coeff_ppr1 (
          stress_op, coeff, proj_scal, proj_tens);
    }
    break; /* RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_REGULAR */

  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL:
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL_SYMM:
    {
      ymir_vec_t         *proj_scal = stokes_problem_nl->proj_scal;
      ymir_vec_t         *proj_tens_stress = rhea_stress_new (ymir_mesh);
      ymir_vec_t         *proj_tens_strain = rhea_strainrate_new (ymir_mesh);
#ifdef RHEA_ENABLE_DEBUG
      ymir_stress_op_coeff_t  coeff_type_linearized;

      /* check coefficient type */
      if (linearization_type ==
          RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL) {
        coeff_type_linearized = YMIR_STRESS_OP_COEFF_ANISOTROPIC_PPR2;
      }
      else if (linearization_type ==
          RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL_SYMM) {
        coeff_type_linearized = YMIR_STRESS_OP_COEFF_ANISOTROPIC_PPR2S;
      }
      else {
        RHEA_ABORT_NOT_REACHED ();
      }
      RHEA_ASSERT (ymir_stress_op_get_coeff_type_linearized (stress_op) ==
                   coeff_type_linearized);
#endif

      /* check input */
      RHEA_ASSERT (stokes_problem_nl->proj_scal != NULL);
      RHEA_ASSERT (rhea_viscosity_check_vec_type (proj_scal));

      /* set coefficient data */
      if (solution_exists) { /* if solution is provided */
        ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);

        /* retrieve velocity of current solution */
        sol_vel = rhea_stokes_problem_retrieve_velocity (
            solution, stokes_problem_nl);

        /* in order to check alignment of primal-dual Newton with regular
         * Newton, store current viscosity */
        ymir_stress_op_get_viscosity (viscosity, stress_op);

        /*
         * Primal-Dual Stress Tensor Update:
         *
         *   stress = stress + stress step
         *          = (linearized stress of velocity step) +
         *            (stress of previous velocity)
         */

        if (step_exists) { /* if step is provided */
          ymir_vec_t         *step_vel = rhea_velocity_new (ymir_mesh);
          ymir_vec_t         *prev_vel = rhea_velocity_new (ymir_mesh);
          ymir_vec_t         *prev_coeff = rhea_viscosity_new (ymir_mesh);
          ymir_vec_t         *temperature = stokes_problem_nl->temperature;
          ymir_vec_t         *weakzone = stokes_problem_nl->weakzone;
          ymir_vec_t         *prev_stress = proj_tens_strain; /* reuse alloc */

          /* retrieve velocity of step */
          rhea_velocity_pressure_copy_components (
              step_vel, NULL, step_vec, stokes_problem_nl->press_elem);
          rhea_stokes_problem_velocity_boundary_set_zero (
              step_vel, stokes_problem_nl);
          RHEA_ASSERT (rhea_velocity_is_valid (step_vel));

          /* compute previous solution velocity (before step) */
          ymir_vec_copy (sol_vel, prev_vel);
          ymir_vec_scale (step_length, step_vel);
          ymir_vec_add (-1.0, step_vel, prev_vel);
          RHEA_ASSERT (rhea_velocity_is_valid (prev_vel));

          /* compute viscosity corresponding to prev solution */
          rhea_stokes_problem_viscosity_compute (
              prev_coeff, NULL, NULL, NULL,
              temperature, weakzone, prev_vel, stokes_problem_nl);

          /* replace coefficient in viscous stress operator (reverted below) */
          ymir_vec_scale (2.0, prev_coeff);
          ymir_stress_op_set_coeff_scal (stress_op, prev_coeff);

          /* compute stress tensor at previous solution velocity */
          ymir_stress_op_optimized_compute_visc_stress (
              prev_stress, prev_vel, stress_op, 0 /* !linearized */);
          RHEA_ASSERT (rhea_stress_is_valid (prev_stress));

          /* compute stress tensor of step velocity (scaled with step length) */
          ymir_stress_op_optimized_compute_visc_stress (
              proj_tens_stress, step_vel, stress_op, 1 /* linearized */);
          RHEA_ASSERT (rhea_stress_is_valid (proj_tens_stress));

          /* combine to get the updated stress tensor */
          ymir_vec_add (1.0, prev_stress, proj_tens_stress);
          RHEA_ASSERT (rhea_stress_is_valid (proj_tens_stress));

          /* destroy */
          rhea_velocity_destroy (step_vel);
          rhea_velocity_destroy (prev_vel);
          rhea_viscosity_destroy (prev_coeff);
        }
        else { /* if step is not provided */
          /* compute stress tensor at current solution velocity */
          ymir_stress_op_optimized_compute_visc_stress (
              proj_tens_stress, sol_vel, stress_op, 0 /* !linearized */);
          RHEA_ASSERT (rhea_stress_is_valid (proj_tens_stress));
        }

        /* calculate how well the primal-dual Newton is aligned with regular
         * Newton */
        {
          double              alignment;

          ymir_stress_op_optimized_compute_strain_rate (
              proj_tens_strain, sol_vel, stress_op);
          ymir_vec_multiply_in1 (viscosity, proj_tens_strain);
          ymir_vec_scale (2.0, proj_tens_strain);
          ymir_mass_apply_gauss (proj_tens_strain);
          alignment = ymir_vec_innerprod (proj_tens_stress, proj_tens_strain);

          RHEA_GLOBAL_INFOF (
              "%s: Alignment <primal-dual stress, regular stress> = %.2e\n",
              __func__, alignment);
        }
        rhea_viscosity_destroy (viscosity);

        /*
         * Primal-Dual Projection Tensors
         */

        /* bound stress tensor for projection within viscous stress operator */
        if (rhea_viscosity_has_yielding (visc_options)) { /* if yielding */
          const double        yield_strength =
            rhea_viscosity_get_yield_strength (visc_options);

          RHEA_ASSERT (0.0 < yield_strength);
          ymir_vec_scale (1.0/yield_strength, proj_tens_stress);
        }
        else { /* otherwise */
#if 0
          ymir_vec_t         *strainrate_sqrt_2inv;

          strainrate_sqrt_2inv = rhea_strainrate_2inv_new (ymir_mesh);
          rhea_strainrate_compute_sqrt_of_2inv (strainrate_sqrt_2inv, sol_vel);
          ymir_vec_divide_in1 (strainrate_sqrt_2inv, proj_tens_stress);
          rhea_strainrate_2inv_destroy (strainrate_sqrt_2inv);
#else
          //TODO code above might be incorrect
          RHEA_ABORT_NOT_REACHED ();
#endif
        }
        ymir_stress_op_tensor_bound_norm (proj_tens_stress, 1.0);

        /* compute and normalize strain rate tensor */
        ymir_stress_op_optimized_compute_strain_rate (
            proj_tens_strain, sol_vel, stress_op);
        ymir_stress_op_tensor_normalize (proj_tens_strain);
      }
      else { /* if solution is not provided */
        ymir_vec_set_zero (proj_tens_stress);
        ymir_vec_set_zero (proj_tens_strain);
      }

      /* pass coefficient data to viscous stress operator */
      ymir_stress_op_set_coeff_ppr2 (
          stress_op, coeff, proj_scal, proj_tens_strain, proj_tens_stress);

      /* destroy */
      rhea_stress_destroy (proj_tens_stress);
      rhea_strainrate_destroy (proj_tens_strain);
    }
    break;
    /* RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL      *
     * RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL_SYMM */

  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_DEV1:
    {
      ymir_vec_t         *proj_scal = stokes_problem_nl->proj_scal;
      ymir_vec_t         *proj_tens_model = rhea_strainrate_new (ymir_mesh);
      ymir_vec_t         *proj_tens_strain = rhea_strainrate_new (ymir_mesh);
#if 0 //#ifdef RHEA_ENABLE_DEBUG
      ymir_stress_op_coeff_t  coeff_type_linearized;

      /* check coefficient type */
      if (linearization_type ==
          RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL) {
        coeff_type_linearized = YMIR_STRESS_OP_COEFF_ANISOTROPIC_PPR2;
      }
      else if (linearization_type ==
          RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL_SYMM) {
        coeff_type_linearized = YMIR_STRESS_OP_COEFF_ANISOTROPIC_PPR2S;
      }
      else {
        RHEA_ABORT_NOT_REACHED ();
      }
      RHEA_ASSERT (ymir_stress_op_get_coeff_type_linearized (stress_op) ==
                   coeff_type_linearized);
#endif

      /* check input */
      RHEA_ASSERT (stokes_problem_nl->proj_scal != NULL);
      RHEA_ASSERT (rhea_viscosity_check_vec_type (proj_scal));

      /* set coefficient data */
      if (solution_exists) { /* if solution is provided */
        /* retrieve velocity of current solution */
        sol_vel = rhea_stokes_problem_retrieve_velocity (
            solution, stokes_problem_nl);

        /*
         * Primal-Dual Stress Tensor Update:
         *
         *   stress = stress + stress step
         *          = (linearized stress of velocity step) +
         *            (stress of previous velocity)
         */

        if (step_exists) { /* if step is provided */
          ymir_vec_t         *step_vel = rhea_velocity_new (ymir_mesh);
          ymir_vec_t         *prev_vel = rhea_velocity_new (ymir_mesh);
          ymir_vec_t         *prev_strain = proj_tens_strain; /* reuse alloc */
          ymir_vec_t         *visc_invisible;

          /* retrieve velocity of step */
          rhea_velocity_pressure_copy_components (
              step_vel, NULL, step_vec, stokes_problem_nl->press_elem);
          rhea_stokes_problem_velocity_boundary_set_zero (
              step_vel, stokes_problem_nl);
          RHEA_ASSERT (rhea_velocity_is_valid (step_vel));

          /* compute previous solution velocity (before step) */
          ymir_vec_copy (sol_vel, prev_vel);
          ymir_vec_scale (step_length, step_vel);
          ymir_vec_add (-1.0, step_vel, prev_vel);
          RHEA_ASSERT (rhea_velocity_is_valid (prev_vel));

          /* compute strain rate at previous solution velocity */
          ymir_stress_op_optimized_compute_strain_rate (
              prev_strain, prev_vel, stress_op);
          RHEA_ASSERT (rhea_stress_is_valid (prev_strain));

          /* compute "strain rate of step velocity" (scaled with step length) */
          visc_invisible = rhea_viscosity_new (ymir_mesh);
          ymir_vec_set_value (visc_invisible, 2.0/2.0); //TODO rem. scaling by 2
          ymir_stress_op_optimized_compute_visc_stress_override_viscosity (
              proj_tens_model, step_vel, stress_op, 1 /* linearized */,
              visc_invisible);
          rhea_viscosity_destroy (visc_invisible);
          RHEA_ASSERT (rhea_stress_is_valid (proj_tens_model));

          /* combine to get the updated stress tensor */
          ymir_vec_add (1.0, prev_strain, proj_tens_model);
          RHEA_ASSERT (rhea_stress_is_valid (proj_tens_model));

          /* destroy */
          rhea_velocity_destroy (step_vel);
          rhea_velocity_destroy (prev_vel);
        }
        else { /* if step is not provided */
          /* compute strain rate at current solution velocity */
          ymir_stress_op_optimized_compute_strain_rate (
              proj_tens_model, sol_vel, stress_op);
          RHEA_ASSERT (rhea_stress_is_valid (proj_tens_model));
        }

        /*
         * Primal-Dual Projection Tensors
         */

        /* compute strain rate tensor */
        ymir_stress_op_optimized_compute_strain_rate (
            proj_tens_strain, sol_vel, stress_op);

        /* normalize tensors */
        ymir_stress_op_tensor_normalize (proj_tens_model);
        ymir_stress_op_tensor_normalize (proj_tens_strain);

        /* calculate how well the model prediction is aligned with the true
         * strain rate tensor */
        {
          ymir_vec_t         *alignment_ptw, *alignment_ptw_mass;
          double              ip, mean, min, max;

          alignment_ptw = rhea_viscosity_new (ymir_mesh);
          ymir_dvec_innerprod_pointwise (alignment_ptw, proj_tens_model,
                                         proj_tens_strain);
          ymir_dvec_fabs (alignment_ptw, alignment_ptw);

          alignment_ptw_mass = rhea_viscosity_new (ymir_mesh);
          ymir_mass_apply (alignment_ptw, alignment_ptw_mass);
          ip = ymir_dvec_innerprod (alignment_ptw_mass, alignment_ptw);
          rhea_viscosity_destroy (alignment_ptw_mass);

          mean = ip / stokes_problem_nl->domain_options->volume;
          min = ymir_dvec_min_global (alignment_ptw);
          max = ymir_dvec_max_global (alignment_ptw);
          rhea_viscosity_destroy (alignment_ptw);

          RHEA_GLOBAL_INFOF (
              "%s: Alignment of strain rate tensors |<model(x),true(x)>_F|, "
              "mean %.2e, min %.2e, max %.2e, max/min %.2e\n",
              __func__, mean, min, max, max/min);
        }
      }
      else { /* if solution is not provided */
        ymir_vec_set_zero (proj_tens_model);
        ymir_vec_set_zero (proj_tens_strain);
      }

      /* pass coefficient data to viscous stress operator */
      ymir_stress_op_set_coeff_ppr2 (
          stress_op, coeff, proj_scal, proj_tens_strain, proj_tens_model);

      /* destroy */
      rhea_strainrate_destroy (proj_tens_model);
      rhea_strainrate_destroy (proj_tens_strain);
    }
    break; /* RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_DEV1 */

  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_DEV2:
    {
      ymir_vec_t         *proj_scal = stokes_problem_nl->proj_scal;
      ymir_vec_t         *proj_tens = stokes_problem_nl->proj_tens;
      ymir_vec_t         *proj_tens_model = rhea_strainrate_new (ymir_mesh);
      ymir_vec_t         *proj_tens_strain = rhea_strainrate_new (ymir_mesh);
#if 0 //#ifdef RHEA_ENABLE_DEBUG
      ymir_stress_op_coeff_t  coeff_type_linearized;

      /* check coefficient type */
      if (linearization_type ==
          RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL) {
        coeff_type_linearized = YMIR_STRESS_OP_COEFF_ANISOTROPIC_PPR2;
      }
      else if (linearization_type ==
          RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL_SYMM) {
        coeff_type_linearized = YMIR_STRESS_OP_COEFF_ANISOTROPIC_PPR2S;
      }
      else {
        RHEA_ABORT_NOT_REACHED ();
      }
      RHEA_ASSERT (ymir_stress_op_get_coeff_type_linearized (stress_op) ==
                   coeff_type_linearized);
#endif

      /* check input */
      RHEA_ASSERT (stokes_problem_nl->proj_scal != NULL);
      RHEA_ASSERT (rhea_viscosity_check_vec_type (proj_scal));

      /* set coefficient data */
      if (solution_exists) { /* if solution is provided */
        ymir_vec_t         *alignment_ptws;

        /* retrieve velocity of current solution */
        sol_vel = rhea_stokes_problem_retrieve_velocity (
            solution, stokes_problem_nl);

        /*
         * Primal-Dual Stress Tensor Update:
         *
         *   stress = stress + stress step
         *          = (linearized stress of velocity step) +
         *            (stress of previous velocity)
         */

        if (step_exists) { /* if step is provided */
          ymir_vec_t         *step_vel = rhea_velocity_new (ymir_mesh);
          ymir_vec_t         *prev_vel = rhea_velocity_new (ymir_mesh);
          ymir_vec_t         *prev_strain = proj_tens_strain; /* reuse alloc */
          ymir_vec_t         *visc_invisible;

          /* retrieve velocity of step */
          rhea_velocity_pressure_copy_components (
              step_vel, NULL, step_vec, stokes_problem_nl->press_elem);
          rhea_stokes_problem_velocity_boundary_set_zero (
              step_vel, stokes_problem_nl);
          RHEA_ASSERT (rhea_velocity_is_valid (step_vel));

          /* compute previous solution velocity (before step) */
          ymir_vec_copy (sol_vel, prev_vel);
          ymir_vec_scale (step_length, step_vel);
          ymir_vec_add (-1.0, step_vel, prev_vel);
          RHEA_ASSERT (rhea_velocity_is_valid (prev_vel));

          /* compute strain rate at previous solution velocity */
          ymir_stress_op_optimized_compute_strain_rate (
              prev_strain, prev_vel, stress_op);
          RHEA_ASSERT (rhea_stress_is_valid (prev_strain));

          /* compute "strain rate of step velocity" (scaled with step length) */
          visc_invisible = rhea_viscosity_new (ymir_mesh);
          ymir_vec_set_value (visc_invisible, 2.0/2.0); //TODO rem. scaling by 2
          ymir_stress_op_optimized_compute_visc_stress_override_viscosity (
              proj_tens_model, step_vel, stress_op, 1 /* linearized */,
              visc_invisible);
          rhea_viscosity_destroy (visc_invisible);
          RHEA_ASSERT (rhea_stress_is_valid (proj_tens_model));

          /* combine to get the updated stress tensor */
          ymir_vec_add (1.0, prev_strain, proj_tens_model);
          RHEA_ASSERT (rhea_stress_is_valid (proj_tens_model));

          /* destroy */
          rhea_velocity_destroy (step_vel);
          rhea_velocity_destroy (prev_vel);
        }
        else { /* if step is not provided */
          /* compute strain rate at current solution velocity */
          ymir_stress_op_optimized_compute_strain_rate (
              proj_tens_model, sol_vel, stress_op);
          RHEA_ASSERT (rhea_stress_is_valid (proj_tens_model));
        }

        /*
         * Primal-Dual Projection Tensors
         */

        /* compute strain rate tensor */
        ymir_stress_op_optimized_compute_strain_rate (
            proj_tens_strain, sol_vel, stress_op);

        /* normalize tensors */
        ymir_stress_op_tensor_normalize (proj_tens_model);
        ymir_stress_op_tensor_normalize (proj_tens_strain);

        /* multiply in alignment of model prediction with true normalized
         * strain rate tensor */
        alignment_ptws = rhea_viscosity_new (ymir_mesh);
        ymir_dvec_innerprod_pointwise (alignment_ptws, proj_tens_model,
                                       proj_tens_strain);
        ymir_dvec_multiply_in (alignment_ptws, proj_scal);
        rhea_viscosity_destroy (alignment_ptws);

        /* print how well the model prediction is aligned with the truth */
        {
          double              alignment_avg;

          alignment_avg = ymir_vec_innerprod (proj_tens_model,
                                              proj_tens_strain);
          RHEA_GLOBAL_INFOF (
              "%s: Alignment of strain rate tensors <model, true> = %.2e\n",
              __func__, alignment_avg);
        }
      }
      else { /* if solution is not provided */
        ymir_vec_set_zero (proj_tens_model);
        ymir_vec_set_zero (proj_tens_strain);
      }

      /* pass coefficient data to viscous stress operator */
      ymir_vec_copy (proj_tens_strain, proj_tens); //TODO change interface
      ymir_stress_op_set_coeff_ppr1 (
          stress_op, coeff, proj_scal, proj_tens);

      /* destroy */
      rhea_strainrate_destroy (proj_tens_model);
      rhea_strainrate_destroy (proj_tens_strain);
    }
    break; /* RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_DEV2 */

  default: /* unknown linearization type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* set up linearization-induced anisotropy */
  {
    ymir_vec_t         *vel_aniso;

    switch (linearization_type) {
    case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_PICARD:
      break;
    case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_REGULAR:
    case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_DEV1:
      /* destroy previous velocity of linearization-induced anisotropy */
      vel_aniso = ymir_stress_op_nsp_lin_aniso_get_velocity (NULL, stress_op);
      if (vel_aniso != NULL) {
        ymir_vec_destroy (vel_aniso);
        ymir_stress_op_nsp_lin_aniso_set_velocity (NULL, stress_op);
      }

      /* set velocity of linearization-induced anisotropy */
      if (solution_exists) { /* if solution is provided */
        rhea_velocity_pressure_create_components (
            &vel_aniso, NULL, solution, stokes_problem_nl->press_elem);
        rhea_stokes_problem_velocity_boundary_set_zero (
            vel_aniso, stokes_problem_nl);
        RHEA_ASSERT (rhea_velocity_is_valid (vel_aniso));
        ymir_stress_op_nsp_lin_aniso_set_velocity (vel_aniso, stress_op);
      }
      break;
  //TODO
  //case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL:
  //case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL_SYMM:
  //case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_DEV2:
    default: /* unknown linearization type */
      RHEA_ABORT_NOT_REACHED ();
    }
  }

  /*
   * Set Up Stokes Preconditioner
   */

#if YMIR_WITH_PETSC
  if (stokes_problem_nl->stokes_pc == NULL) {
    /* create Stokes preconditioner */
    stokes_problem_nl->stokes_pc = ymir_stokes_pc_new_multi_solver (
        stokes_op, 1 /* linearized */, stokes_problem_nl->num_krylov_solvers);
  }
  else {
    /* update Stokes preconditioner */
    ymir_nlstokes_pc_recompute (stokes_problem_nl->stokes_pc);
  }
#else
  RHEA_ABORT_NOT_REACHED ();
#endif

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static int          rhea_stokes_problem_nonlinear_mesh_data_exists (
                                    rhea_stokes_problem_t *stokes_problem_nl);

static int
rhea_stokes_problem_nonlinear_solver_data_exists (
                                    rhea_stokes_problem_t *stokes_problem_nl)
{
  int                 mesh_data_exists, stokes_exists, norm_exists;

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);

  /* check if mesh data exists */
  mesh_data_exists =
    rhea_stokes_problem_nonlinear_mesh_data_exists (stokes_problem_nl);

  /* check if Stokes operator exists */
  stokes_exists = (stokes_problem_nl->stokes_op != NULL);

  /* check if norm operator exists */
  if (stokes_problem_nl->norm_type == RHEA_STOKES_NORM_HMINUS1_L2) {
    norm_exists = (stokes_problem_nl->norm_op != NULL);
  }
  else {
    norm_exists = 1;
  }

  return (mesh_data_exists && stokes_exists && norm_exists);
}

/**
 * Initializes user data of the nonlinear problem to be able to run the Newton
 * solver.
 * (Callback function for Newton's method)
 */
static void
rhea_stokes_problem_nonlinear_create_solver_data_fn (ymir_vec_t *solution,
                                                     void *data)
{
  rhea_stokes_problem_t *stokes_problem_nl = data;
  const int           solver_data_exists =
    rhea_stokes_problem_nonlinear_solver_data_exists (stokes_problem_nl);
  const char         *vtk_path = stokes_problem_nl->solver_vtk_path;
  const int           nonzero_init_guess = (solution != NULL);

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (__func__, "nonzero init guess=%i",
                                 nonzero_init_guess);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (
      rhea_stokes_problem_nonlinear_mesh_data_exists (stokes_problem_nl));

  /* set VTK path for debuggin */
  if (vtk_path != NULL) {
    char                path[BUFSIZ];

    snprintf (path, BUFSIZ, "%s%s%02i", vtk_path,
              RHEA_STOKES_PROBLEM_IO_LABEL_NL_ITER, 0 /* iteration */);
    ymir_vtk_set_debug_path (path);
  }

  /* compute weak zone */
  rhea_stokes_problem_weakzone_compute (stokes_problem_nl->weakzone,
                                        stokes_problem_nl);

  /* set up Stokes operator */
  if (solver_data_exists) {
    /* update Stokes operator */
    rhea_stokes_problem_nonlinear_update_operator_fn (solution,
                                                      stokes_problem_nl);
  }
  else {
    /* compute viscosity */
    rhea_stokes_problem_compute_coefficient (stokes_problem_nl, solution,
                                             !nonzero_init_guess);

    /* create Stokes operator */
    stokes_problem_nl->stokes_op = ymir_stokes_op_new_ext (
        stokes_problem_nl->coeff, stokes_problem_nl->vel_dir,
        NULL /* Robin BC's */,
        NULL /* deprecated */,
        NULL /* deprecated */,
        stokes_problem_nl->press_elem,
        stokes_problem_nl->domain_options->center,
        stokes_problem_nl->domain_options->moment_of_inertia);
  }

  /* set coefficient types and shifts */
  {
    ymir_stress_op_t   *stress_op = stokes_problem_nl->stokes_op->stress_op;
    ymir_stress_op_coeff_t  coeff_type_linearized;
    rhea_viscosity_options_t *visc_options = stokes_problem_nl->visc_options;

    /* set viscous stress coefficient */
    ymir_stress_op_set_coeff_type (
        stress_op, YMIR_STRESS_OP_COEFF_ISOTROPIC_SCAL);
    ymir_stress_op_set_coeff_shift (
        stress_op, rhea_viscosity_get_visc_shift (visc_options));

    /* set the linearized part of the viscous stress coefficient */
    switch (stokes_problem_nl->linearization_type) {
    case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_PICARD:
      coeff_type_linearized = YMIR_STRESS_OP_COEFF_ISOTROPIC_SCAL;
      break;
    case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_REGULAR:
      coeff_type_linearized = YMIR_STRESS_OP_COEFF_ANISOTROPIC_PPR1;
      break;
    case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL:
      coeff_type_linearized = YMIR_STRESS_OP_COEFF_ANISOTROPIC_PPR2;
      break;
    case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL_SYMM:
      coeff_type_linearized = YMIR_STRESS_OP_COEFF_ANISOTROPIC_PPR2S;
      break;
    case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_DEV1:
      coeff_type_linearized = YMIR_STRESS_OP_COEFF_ANISOTROPIC_PPR2S;
      break;
    case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_DEV2:
      coeff_type_linearized = YMIR_STRESS_OP_COEFF_ANISOTROPIC_PPR1;
      break;
    default: /* unknown linearization type */
      RHEA_ABORT_NOT_REACHED ();
    }
    ymir_stress_op_set_coeff_type_linearized (
        stress_op, coeff_type_linearized);
    ymir_stress_op_set_coeff_shift_proj (
        stress_op, rhea_viscosity_get_visc_shift_proj (visc_options));
  }

  /* create H^-1 norm operator */
  if (!solver_data_exists &&
      stokes_problem_nl->norm_type == RHEA_STOKES_NORM_HMINUS1_L2) {
    stokes_problem_nl->norm_op = ymir_Hminus1_norm_op_new (
        stokes_problem_nl->ymir_mesh,
        stokes_problem_nl->norm_op_mass_scaling);
  }

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

/**
 * Clears the data of the nonlinear problem after Newton solve.
 * (Callback function for Newton's method)
 */
static void
rhea_stokes_problem_nonlinear_clear_solver_data_fn (void *data)
{
  rhea_stokes_problem_t *stokes_problem_nl = data;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);

  /* destroy H^-1 norm operator */
  if (stokes_problem_nl->norm_type == RHEA_STOKES_NORM_HMINUS1_L2) {
    RHEA_ASSERT (stokes_problem_nl->norm_op != NULL);
    ymir_Hminus1_norm_op_destroy (stokes_problem_nl->norm_op);
    stokes_problem_nl->norm_op = NULL;
  }

  /* destroy Stokes preconditioner */
  if (stokes_problem_nl->stokes_pc != NULL) {
    ymir_stokes_pc_destroy (stokes_problem_nl->stokes_pc);
    stokes_problem_nl->stokes_pc = NULL;
  }

  /* destroy Stokes operator */
  if (stokes_problem_nl->stokes_op != NULL) {
    ymir_stress_op_t   *stress_op = stokes_problem_nl->stokes_op->stress_op;
    ymir_vec_t         *vel_aniso;

    vel_aniso = ymir_stress_op_nsp_lin_aniso_get_velocity (NULL, stress_op);
    if (vel_aniso != NULL) {
      ymir_vec_destroy (vel_aniso);
      ymir_stress_op_nsp_lin_aniso_set_velocity (NULL, stress_op);
    }

    ymir_stokes_op_destroy (stokes_problem_nl->stokes_op);
    stokes_problem_nl->stokes_op = NULL;
  }

  /* clear VTK path for debugging */
  ymir_vtk_set_debug_path ("");

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

/**
 * Computes the negative gradient of the objective functional.
 * (Callback function for Newton's method)
 */
static void
rhea_stokes_problem_nonlinear_compute_negative_gradient_fn (
                                            ymir_vec_t *neg_gradient,
                                            ymir_vec_t *solution, void *data)
{
  rhea_stokes_problem_t *stokes_problem_nl = data;
  ymir_vec_t         *rhs_vel = stokes_problem_nl->rhs_vel;
  ymir_vec_t         *rhs_vel_nonzero_dirichlet =
                        stokes_problem_nl->rhs_vel_nonzero_dirichlet;
  ymir_vec_t        **rhs_vel_face_nonzero_neumann =
                        stokes_problem_nl->rhs_vel_face_nonzero_neumann;
  ymir_vec_t         *rhs_vel_press = stokes_problem_nl->rhs_vel_press;
  ymir_stokes_op_t   *stokes_op = stokes_problem_nl->stokes_op;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* check input */
#ifdef RHEA_ENABLE_DEBUG
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
      rhs_vel_press,                /* output: right-hand side */
      rhs_vel,                      /* input: volume forcing */
      rhs_vel_face_nonzero_neumann, /* input: nonzero Neumann forcing */
      rhs_vel_nonzero_dirichlet,    /* input: nonzero Dirichlet boundary */
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
    ymir_stokes_pc_apply_stokes_op (solution, neg_gradient, stokes_op, 0, 1);
    ymir_vec_add (-1.0, rhs_vel_press, neg_gradient);
    ymir_vec_scale (-1.0, neg_gradient);
  }
  else {
    /* set residual to be the right-hand side */
    ymir_vec_copy (rhs_vel_press, neg_gradient);
  }
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (neg_gradient));

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

/**
 * Computes the norm of the gradient.
 * (Callback function for Newton's method)
 */
static double
rhea_stokes_problem_nonlinear_compute_gradient_norm_fn (
                                                      ymir_vec_t *neg_gradient,
                                                      void *data,
                                                      double *norm_comp)
{
  rhea_stokes_problem_t *stokes_problem_nl = data;
  rhea_stokes_norm_type_t norm_type = stokes_problem_nl->norm_type;
  ymir_Hminus1_norm_op_t *norm_op = stokes_problem_nl->norm_op;
  ymir_vec_t         *innerprod_arg_left, *innerprod_arg_right;
  double              ip_vel, ip_press;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (neg_gradient != NULL);
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (neg_gradient));
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (neg_gradient));
  RHEA_ASSERT (stokes_problem_nl->norm_type != RHEA_STOKES_NORM_NONE);
  RHEA_ASSERT (stokes_problem_nl->press_elem != NULL);

  /* compute inner products */
  if (ymir_stokes_pc_block_type == YMIR_STOKES_PC_BLOCK_L) { /* if left PC */
#if 0
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
#else
    //TODO code above might not work, needs implementation
    RHEA_ABORT_NOT_REACHED ();
#endif
  }
  else { /* otherwise right or symmetric preconditioning */
    innerprod_arg_left = neg_gradient;
    innerprod_arg_right = neg_gradient;
  }
  rhea_stokes_norm_innerprod (
      &ip_vel, &ip_press, innerprod_arg_left, innerprod_arg_right,
      norm_type, norm_op, stokes_problem_nl->press_elem);
  RHEA_ASSERT (isfinite (ip_vel));
  RHEA_ASSERT (isfinite (ip_press));
  RHEA_ASSERT (0.0 <= ip_vel);
  RHEA_ASSERT (0.0 <= ip_press);

  /* destroy */
  if (ymir_stokes_pc_block_type == YMIR_STOKES_PC_BLOCK_L) {
    ymir_vec_destroy (innerprod_arg_right);
  }

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);

  /* calculate norms */
  if (norm_comp != NULL) {
    norm_comp[0] = sqrt (ip_vel);
    norm_comp[1] = sqrt (ip_press);
  }
  return sqrt (ip_vel + ip_press);
}

/**
 * Applies Hessian operator.
 * (Callback function for Newton's method)
 */
static void
rhea_stokes_problem_nonlinear_apply_hessian_fn (ymir_vec_t *out, ymir_vec_t *in,
                                                void *data)
{
  rhea_stokes_problem_t *stokes_problem_nl = data;
  ymir_stokes_op_t   *stokes_op = stokes_problem_nl->stokes_op;
  const int           linearized = 1;
  const int           dirty = 1;

  /* check input */
#ifdef RHEA_ENABLE_DEBUG
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->stokes_op != NULL);
  {
    ymir_stress_op_t   *stress_op = stokes_problem_nl->stokes_op->stress_op;

    RHEA_ASSERT (ymir_stress_op_has_linearized (stress_op));
  }
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (out));
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (in));
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (in));
#endif

  /* apply */
  ymir_stokes_pc_apply_stokes_op (in, out, stokes_op, linearized, dirty);
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (out));
}

/**
 * Modifies the Hessian system before launching the solver for this system.
 * (Callback function for Newton's method)
 */
static void
rhea_stokes_problem_nonlinear_modify_hessian_system_fn (
                                                      ymir_vec_t *neg_gradient,
                                                      ymir_vec_t *solution,
                                                      void *data)
{
  const int           lin_aniso_pc =
    rhea_stokes_problem_nonlinear_lin_aniso_pc;
  const int           check_lin_aniso =
    rhea_stokes_problem_nonlinear_lin_aniso_check;
  rhea_stokes_problem_t *stokes_problem_nl = data;
  const rhea_stokes_problem_nonlinear_linearization_t linearization_type =
    stokes_problem_nl->linearization_type;
  ymir_pressure_elem_t *press_elem = stokes_problem_nl->press_elem;
  ymir_stress_op_t   *stress_op;

  /* get viscous stress operator */
  RHEA_ASSERT (stokes_problem_nl->stokes_op != NULL);
  stress_op = stokes_problem_nl->stokes_op->stress_op;

  /* exit if nothing to do */
  if (solution == NULL) {
    return;
  }

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (neg_gradient));
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (solution));
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (neg_gradient));
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (solution));

  /* precondition linearization-induced anisotropy */
  switch (linearization_type) {
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_PICARD:
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_REGULAR:
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_DEV1:
#ifdef RHEA_ENABLE_DEBUG
    if (lin_aniso_pc || check_lin_aniso) {
      ymir_vec_t         *vel_aniso;

      vel_aniso = ymir_stress_op_nsp_lin_aniso_get_velocity (NULL, stress_op);
      RHEA_ASSERT (vel_aniso != NULL);
    }
#endif

    /* check for linearization-induced anisotropy in (neg.) gradient */
    if (check_lin_aniso) {
      ymir_vec_t         *neg_grad_vel;
      double              neg_grad_aniso_part, normL2;

      rhea_velocity_pressure_create_components (
          &neg_grad_vel, NULL, neg_gradient, press_elem);
      neg_grad_aniso_part = ymir_stress_op_nsp_lin_aniso_compute_aniso_part (
          neg_grad_vel, stress_op, 1 /* dual/residual space */);
      ymir_vec_destroy (neg_grad_vel);

      if (isfinite (neg_grad_aniso_part)) {
        rhea_stokes_norm_compute (
            &normL2, NULL, neg_gradient, RHEA_STOKES_NORM_L2_DUAL, NULL,
            press_elem);
        RHEA_GLOBAL_INFOF (
            "%s: Linearization-induced anisotropy part in gradient: "
            "abs %+.2e, rel %+.2e\n",
            __func__, neg_grad_aniso_part, neg_grad_aniso_part/normL2);
      }
    }

    /* apply linearization-induced anisotropy correction */
    if (lin_aniso_pc) {
      ymir_vec_t         *neg_grad_vel;
      double              aniso_corr;

      /* apply correction to velocity of gradient */
      rhea_velocity_pressure_create_components (
          &neg_grad_vel, NULL, neg_gradient, press_elem);
      aniso_corr = ymir_stress_op_nsp_lin_aniso_set_correction (
          neg_grad_vel, stress_op, 1 /* dual/residual space */);
      ymir_stress_op_nsp_lin_aniso_apply_correction (
          neg_grad_vel, stress_op, 1 /* dual/residual space */);
      RHEA_GLOBAL_INFOF (
          "%s: Linearization-induced anisotropy correction: %+.3e\n",
          __func__, aniso_corr);

      /* set new velocity of gradient */
      rhea_stokes_problem_velocity_boundary_set_zero (
          neg_grad_vel, stokes_problem_nl);
      rhea_velocity_pressure_set_components (
          neg_gradient, neg_grad_vel, NULL, press_elem);

      /* activate anisotropy correction in Stokes operator */
      ymir_stokes_op_nsp_lin_aniso_set_proj (
          stokes_problem_nl->stokes_op,
          YMIR_STRESS_OP_NSP_LIN_ANISO_PROJECT_RES);

      /* destroy */
      ymir_vec_destroy (neg_grad_vel);
    }
    break;
//TODO
//case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL:
//case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL_SYMM:
//case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_DEV2:
  default: /* unknown linearization type */
    RHEA_ABORT_NOT_REACHED ();
  }

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

/**
 * Applies inexact inverse of the Hessian operator.
 * (Callback function for Newton's method)
 */
static int
rhea_stokes_problem_nonlinear_solve_hessian_system_fn (
                                            ymir_vec_t *step,
                                            ymir_vec_t *neg_gradient,
                                            const int lin_iter_max,
                                            const double lin_res_norm_rtol,
                                            const int nonzero_initial_guess,
                                            void *data,
                                            int *lin_iter_count)
{
  const int           krylov_solver_idx = 0;
  const int           check_lin_aniso =
    rhea_stokes_problem_nonlinear_lin_aniso_check;
  rhea_stokes_problem_t *stokes_problem_nl = data;
  int                 stop_reason;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (
      __func__, "lin_iter_max=%i, lin_rtol=%.1e, nonzero_init_guess=%i",
      lin_iter_max, lin_res_norm_rtol, nonzero_initial_guess);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->stokes_pc != NULL);
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (step));
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (neg_gradient));
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (neg_gradient));

  /* run Krylov solver for the linearized system */
  stop_reason = rhea_stokes_problem_run_krylov_solver (
      step, neg_gradient, nonzero_initial_guess, lin_iter_max,
      lin_res_norm_rtol, stokes_problem_nl, krylov_solver_idx,
      lin_iter_count, NULL /* residual_reduction */);
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (step));

  /* check for linearization-induced anisotropy in (neg.) gradient & step */
  if (check_lin_aniso) {
    ymir_pressure_elem_t *press_elem = stokes_problem_nl->press_elem;
    ymir_stress_op_t   *stress_op = stokes_problem_nl->stokes_op->stress_op;
    ymir_vec_t         *neg_grad_vel, *step_vel;
    double              neg_grad_aniso_part, step_aniso_part, normL2;

    rhea_velocity_pressure_create_components (
        &neg_grad_vel, NULL, neg_gradient, press_elem);
    neg_grad_aniso_part = ymir_stress_op_nsp_lin_aniso_compute_aniso_part (
        neg_grad_vel, stress_op, 1 /* dual/residual space */);
    ymir_vec_destroy (neg_grad_vel);

    rhea_velocity_pressure_create_components (
        &step_vel, NULL, step, press_elem);
    step_aniso_part = ymir_stress_op_nsp_lin_aniso_compute_aniso_part (
        step_vel, stress_op, 0 /* primal space */);
    ymir_vec_destroy (step_vel);

    if (isfinite (neg_grad_aniso_part)) {
      rhea_stokes_norm_compute (
          &normL2, NULL, neg_gradient, RHEA_STOKES_NORM_L2_DUAL, NULL,
          press_elem);
      RHEA_GLOBAL_INFOF (
          "%s: Linearization-induced anisotropy part in gradient: "
          "abs %+.2e, rel %+.2e\n",
          __func__, neg_grad_aniso_part, neg_grad_aniso_part/normL2);
    }

    if (isfinite (step_aniso_part)) {
      rhea_stokes_norm_compute (
          &normL2, NULL, step, RHEA_STOKES_NORM_L2_PRIMAL, NULL, press_elem);
      RHEA_GLOBAL_INFOF (
          "%s: Linearization-induced anisotropy part in step:     "
          "abs %+.2e, rel %+.2e\n",
          __func__, step_aniso_part, step_aniso_part/normL2);
    }
  }

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);

  /* return stopping reason */
  return stop_reason;
}

/**
 * Performs AMR/grid continuation in between Newton steps.
 * (Callback function for Newton's method)
 */
static int
rhea_stokes_problem_nonlinear_amr_fn (ymir_vec_t **solution, const int iter,
                                      void *data)
{
  rhea_stokes_problem_t *stokes_problem_nl = data;
  int                 amr_iter;

  /* exit if nothing to do */
  if (stokes_problem_nl->p4est == NULL ||
      stokes_problem_nl->discr_options == NULL) {
    return 0;
  }

  /* run AMR */
  rhea_stokes_problem_set_velocity_pressure (stokes_problem_nl, *solution);
  amr_iter = rhea_stokes_problem_nonlinear_amr (
      stokes_problem_nl, stokes_problem_nl->p4est,
      stokes_problem_nl->discr_options, iter);

  /* (possibly) recover new vector for solution */
  if (0 < amr_iter) { /* if mesh was adapted */
    stokes_problem_nl->mesh_modified_by_solver = 1;
    *solution = rhea_stokes_problem_get_velocity_pressure (stokes_problem_nl);

    /* project out null spaces in solution */
    rhea_stokes_problem_project_out_nullspace (*solution, stokes_problem_nl);

    return 1;
  }
  else { /* if mesh stayed the same */
    return 0;
  }
}

/**
 * Writes or prints output at the beginning of a Newton step.
 * (Callback function for Newton's method)
 */
static void
rhea_stokes_problem_nonlinear_output_prestep_fn (ymir_vec_t *solution,
                                                 const int iter, void *data)
{
  rhea_stokes_problem_t *stokes_problem_nl = data;
  const char         *bin_path = stokes_problem_nl->solver_bin_path;
  const char         *vtk_path = stokes_problem_nl->solver_vtk_path;
  const double        domain_vol = stokes_problem_nl->domain_options->volume;
  ymir_mesh_t          *ymir_mesh = stokes_problem_nl->ymir_mesh;
  ymir_pressure_elem_t *press_elem = stokes_problem_nl->press_elem;
  ymir_vec_t         *velocity, *pressure, *viscosity;
  ymir_vec_t         *velocity_surf, *stress_norm_surf, *viscosity_surf;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (__func__, "newton_iter=%i", iter);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->bounds_marker != NULL);
  RHEA_ASSERT (stokes_problem_nl->yielding_marker != NULL);

  /* get volume fields */
  rhea_velocity_pressure_create_components (&velocity, &pressure, solution,
                                            press_elem);
  viscosity = rhea_viscosity_new (ymir_mesh);
  rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem_nl);

  /* get surface fields */
  velocity_surf = rhea_velocity_surface_new (ymir_mesh);
  rhea_velocity_surface_interpolate (velocity_surf, velocity);
  stress_norm_surf = rhea_stress_surface_new (ymir_mesh);
  rhea_stokes_problem_stress_compute_normal_at_surface (stress_norm_surf,
                                                        solution,
                                                        stokes_problem_nl);
  viscosity_surf = rhea_viscosity_surface_new (ymir_mesh);
  rhea_viscosity_surface_interpolate (viscosity_surf, viscosity,
                                      stokes_problem_nl->visc_options->min,
                                      stokes_problem_nl->visc_options->max);

  RHEA_GLOBAL_STATISTICS ("========================================\n");
  RHEA_GLOBAL_STATISTICSF ("Stokes flow statistics, newton_iter=%i\n", iter);
  RHEA_GLOBAL_STATISTICS ("----------------------------------------\n");

  /* print velocity statistics */
  {
    double              magn_min_mm_yr, magn_max_mm_yr, magn_mean_mm_yr;
    double              lith_magn_max_mm_yr, lith_magn_mean_mm_yr;
    double              surf_magn_min_mm_yr, surf_magn_max_mm_yr,
                        surf_magn_mean_mm_yr;
    double              surf_lith_magn_max_mm_yr, surf_lith_magn_mean_mm_yr;
    double              mean_rot_axis[3];
    int                 has_mean_rot = 0;

    rhea_velocity_stats_get_global (
        &magn_min_mm_yr, &magn_max_mm_yr, &magn_mean_mm_yr, velocity,
        stokes_problem_nl->domain_options, stokes_problem_nl->temp_options);
    rhea_velocity_stats_get_global_lithosphere (
        &lith_magn_max_mm_yr, &lith_magn_mean_mm_yr, velocity, viscosity,
        stokes_problem_nl->domain_options, stokes_problem_nl->temp_options);
    rhea_velocity_stats_get_global_surface (
        &surf_magn_min_mm_yr, &surf_magn_max_mm_yr, &surf_magn_mean_mm_yr,
        velocity, stokes_problem_nl->domain_options,
        stokes_problem_nl->temp_options);
    rhea_velocity_stats_get_global_surface_lithosphere (
        &surf_lith_magn_max_mm_yr, &surf_lith_magn_mean_mm_yr,
        velocity, viscosity, stokes_problem_nl->domain_options,
        stokes_problem_nl->temp_options);
    has_mean_rot = rhea_stokes_problem_velocity_compute_mean_rotation (
        mean_rot_axis, velocity, stokes_problem_nl);

    RHEA_GLOBAL_STATISTICSF (
        "Velocity magn [mm/yr] vol:  global min %.3e, max %.3e, mean %.3e\n",
        magn_min_mm_yr, magn_max_mm_yr, magn_mean_mm_yr);
    RHEA_GLOBAL_STATISTICSF (
        "Velocity magn [mm/yr] vol:  lithosphere %*s max %.3e, mean %.3e\n",
        /*space*/ 9, "", lith_magn_max_mm_yr, lith_magn_mean_mm_yr);
    RHEA_GLOBAL_STATISTICSF (
        "Velocity magn [mm/yr] surf: global min %.3e, max %.3e, mean %.3e\n",
        surf_magn_min_mm_yr, surf_magn_max_mm_yr,
        surf_magn_mean_mm_yr);
    RHEA_GLOBAL_STATISTICSF (
        "Velocity magn [mm/yr] surf: lithosphere %*s max %.3e, mean %.3e\n",
        /*space*/ 9, "", surf_lith_magn_max_mm_yr,
        surf_lith_magn_mean_mm_yr);
    if (has_mean_rot) {
      RHEA_GLOBAL_STATISTICSF (
        "Velocity mean rot axis: x,y,z = %+.3e , %+.3e , %+.3e\n",
        mean_rot_axis[0], mean_rot_axis[1], mean_rot_axis[2]);
    }
  }

  /* print pressure statistics */
  {
    double              abs_min_Pa, abs_max_Pa, mean_Pa;

    rhea_pressure_stats_get_global (
        &abs_min_Pa, &abs_max_Pa, &mean_Pa, pressure,
        press_elem, stokes_problem_nl->domain_options,
        stokes_problem_nl->temp_options, stokes_problem_nl->visc_options);

    RHEA_GLOBAL_STATISTICSF (
        "Pressure [Pa]: global abs min %.3e, abs max %.3e, mean %.3e\n",
        abs_min_Pa, abs_max_Pa, mean_Pa);
  }

  /* print strain rate statistics */
  {
    double              min_1_s, max_1_s, mean_1_s;

    rhea_strainrate_stats_get_global (
        &min_1_s, &max_1_s, &mean_1_s, velocity,
        stokes_problem_nl->domain_options, stokes_problem_nl->temp_options);
    RHEA_GLOBAL_STATISTICSF (
        "Strain rate (sqrt 2nd inv) [1/s]: global min %.3e, max %.3e, "
        "max/min %.3e, mean %.3e\n",
        min_1_s, max_1_s, max_1_s/min_1_s, mean_1_s);
  }

  /* print stress statistics */
  {
    double              min_Pa, max_Pa, mean_Pa;
    double              surf_min_Pa, surf_max_Pa, surf_mean_Pa;

    rhea_stress_stats_get_global (
        &min_Pa, &max_Pa, &mean_Pa, velocity, viscosity,
        stokes_problem_nl->domain_options, stokes_problem_nl->temp_options,
        stokes_problem_nl->visc_options);
    rhea_stress_surface_stats_get_global (
        &surf_min_Pa, &surf_max_Pa, &surf_mean_Pa, stress_norm_surf,
        stokes_problem_nl->domain_options, stokes_problem_nl->temp_options,
        stokes_problem_nl->visc_options);

    RHEA_GLOBAL_STATISTICSF (
        "Visc stress (sqrt 2nd inv) [Pa]:  global min %.3e, max %.3e, "
        "max/min %.3e, mean %.3e\n",
        min_Pa, max_Pa, max_Pa/min_Pa, mean_Pa);
    RHEA_GLOBAL_STATISTICSF (
        "Stress normal at surface [Pa]:    global min %+.3e, max %+.3e, "
        "mean %+.3e\n",
        surf_min_Pa, surf_max_Pa, surf_mean_Pa);
  }

  /* print viscosity statistics */
  {
    double              min_Pas, max_Pas, mean_Pas;
    double              upper_mantle_mean_Pas, lower_mantle_mean_Pas;
    double              lith_mean_Pas, asth_mean_Pas;
    double              bounds_vol_min, bounds_vol_max;
    double              yielding_vol;
    double              lith_vol, asth_vol;

    rhea_viscosity_stats_get_global (&min_Pas, &max_Pas, &mean_Pas, viscosity,
                                     stokes_problem_nl->visc_options);
    rhea_viscosity_stats_get_regional (&upper_mantle_mean_Pas,
                                       &lower_mantle_mean_Pas,
                                       &lith_mean_Pas,
                                       &asth_mean_Pas, viscosity,
                                       stokes_problem_nl->visc_options);

    rhea_viscosity_stats_get_bounds_volume (
        &bounds_vol_min, &bounds_vol_max, stokes_problem_nl->bounds_marker);
    yielding_vol = rhea_viscosity_stats_get_yielding_volume (
        stokes_problem_nl->yielding_marker);
    lith_vol = rhea_viscosity_stats_get_lithosphere_volume (
        viscosity, stokes_problem_nl->visc_options);
    asth_vol = rhea_viscosity_stats_get_asthenosphere_volume (
        viscosity, stokes_problem_nl->visc_options);

    RHEA_GLOBAL_STATISTICSF (
        "Viscosity [Pa*s]: global min %.3e, max %.3e, max/min %.3e, "
        "mean %.3e\n",
        min_Pas, max_Pas, max_Pas/min_Pas, mean_Pas);
    RHEA_GLOBAL_STATISTICSF (
        "Viscosity [Pa*s]: upp-m mean %.3e, low-m mean %.3e\n",
        upper_mantle_mean_Pas, lower_mantle_mean_Pas);
    RHEA_GLOBAL_STATISTICSF (
        "Viscosity [Pa*s]: asth. mean %.3e, lith. mean %.3e\n",
        asth_mean_Pas, lith_mean_Pas);

    RHEA_GLOBAL_STATISTICSF (
        "Visc min bounds volume: abs %.8e, rel %.6f\n",
        bounds_vol_min, bounds_vol_min/domain_vol);
    RHEA_GLOBAL_STATISTICSF (
        "Visc max bounds volume: abs %.8e, rel %.6f\n",
        bounds_vol_max, bounds_vol_max/domain_vol);
    RHEA_GLOBAL_STATISTICSF (
        "Yielding volume:        abs %.8e, rel %.6f\n",
        yielding_vol, yielding_vol/domain_vol);
    RHEA_GLOBAL_STATISTICSF (
        "Lithoshpere volume:     abs %.8e, rel %.6f\n",
        lith_vol, lith_vol/domain_vol);
    RHEA_GLOBAL_STATISTICSF (
        "Asthenoshpere volume:   abs %.8e, rel %.6f\n",
        asth_vol, asth_vol/domain_vol);
  }

  /* print plate velocities */
  if (stokes_problem_nl->plate_options != NULL &&
      0 < rhea_plate_get_n_plates (stokes_problem_nl->plate_options)) {
    rhea_plate_options_t *plate_options = stokes_problem_nl->plate_options;
    const int           n_plates = rhea_plate_get_n_plates (plate_options);
    int                 n_plates_print;
    int                 pid;
    ymir_vec_t         *velocity_surf_mm_yr = ymir_vec_clone (velocity_surf);
    double             *mean_vel_magn_mm_yr = RHEA_ALLOC (double, n_plates);

    rhea_velocity_convert_to_dimensional_mm_yr (
        velocity_surf_mm_yr, stokes_problem_nl->domain_options,
        stokes_problem_nl->temp_options);
    rhea_plate_velocity_get_mean_magnitude_all (
        mean_vel_magn_mm_yr, velocity_surf_mm_yr, NULL /* plate_label */,
        1 /* project_out_mean_rot */, plate_options);
    ymir_vec_destroy (velocity_surf_mm_yr);

    /* set #plates to print */
    switch (stokes_problem_nl->domain_options->shape) {
    case RHEA_DOMAIN_SHELL:
      n_plates_print = RHEA_PLATE_EARTH_YZ+1; /* MORVEL(25) only */
      break;
    default:
      n_plates_print = n_plates;
    }
    n_plates_print = SC_MIN (n_plates_print, n_plates);

    /* print */
    for (pid = 0; pid < n_plates_print; pid++) {
      RHEA_GLOBAL_STATISTICSF (
          "Plate velocity: idx %i, mean velocity magn \"%g mm/yr\"\n",
          pid, mean_vel_magn_mm_yr[pid]);
    }
    RHEA_FREE (mean_vel_magn_mm_yr);
  }

  /* write binary files */
  if (bin_path != NULL) {
    char                path[BUFSIZ];

    snprintf (path, BUFSIZ, "%s%s%02i", bin_path,
              RHEA_STOKES_PROBLEM_IO_LABEL_NL_ITER, iter);
    rhea_stokes_problem_write (path, stokes_problem_nl);
  }

  /* create visualization */
  if (vtk_path != NULL) {
    char                path[BUFSIZ];

    /* set path */
    snprintf (path, BUFSIZ, "%s%s%02i", vtk_path,
              RHEA_STOKES_PROBLEM_IO_LABEL_NL_ITER, iter);

    /* write VTK */
    rhea_vtk_write_nonlinear_stokes_iteration (
        path, velocity, pressure, viscosity,
        stokes_problem_nl->bounds_marker, stokes_problem_nl->yielding_marker);
    rhea_vtk_write_nonlinear_stokes_iteration_surf (
        path, velocity_surf, stress_norm_surf, viscosity_surf);

    /* set VTK path for debugging within ymir */
    ymir_vtk_set_debug_path (path);
  }

  RHEA_GLOBAL_STATISTICS ("========================================\n");

  /* destroy */
  rhea_velocity_destroy (velocity);
  rhea_pressure_destroy (pressure);
  rhea_viscosity_destroy (viscosity);
  rhea_velocity_surface_destroy (velocity_surf);
  rhea_stress_surface_destroy (stress_norm_surf);
  rhea_viscosity_surface_destroy (viscosity_surf);

  RHEA_GLOBAL_VERBOSEF_FN_END (__func__, "newton_iter=%i", iter);
}

/******************************************************************************
 * Nonlinear Stokes Problem
 *****************************************************************************/

static int
rhea_stokes_problem_nonlinear_mesh_data_exists (
                                    rhea_stokes_problem_t *stokes_problem_nl)
{
  int                 mesh_exists, vectors_exist;
  int                 weak_exists;
  int                 rhs_nonzero_dir_exists, rhs_nonzero_neu_exists;
  int                 nl_vec_exist;

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);

  /* check if mesh and corresponding data exist */
  mesh_exists = (
      stokes_problem_nl->ymir_mesh != NULL &&
      stokes_problem_nl->press_elem != NULL);
  vectors_exist = (
      stokes_problem_nl->coeff != NULL &&
      stokes_problem_nl->bounds_marker != NULL &&
      stokes_problem_nl->yielding_marker != NULL &&
      stokes_problem_nl->rhs_vel != NULL &&
      stokes_problem_nl->rhs_vel_press != NULL);

  /* check if not required vectors exist */
  weak_exists = (
      stokes_problem_nl->weakzone_compute_fn == NULL ||
      stokes_problem_nl->weakzone != NULL);
  rhs_nonzero_dir_exists = (
      stokes_problem_nl->rhs_vel_nonzero_dir_compute_fn == NULL ||
      stokes_problem_nl->rhs_vel_nonzero_dirichlet != NULL);
  rhs_nonzero_neu_exists = (
      stokes_problem_nl->rhs_vel_nonzero_neu_compute_fn == NULL ||
      stokes_problem_nl->rhs_vel_face_nonzero_neumann != NULL);

  /* check if nonlinear solver data exists */
  nl_vec_exist = (
      stokes_problem_nl->sol_vel != NULL &&
      stokes_problem_nl->newton_neg_gradient_vec != NULL &&
      stokes_problem_nl->newton_step_vec != NULL);
  switch (stokes_problem_nl->linearization_type) {
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_PICARD:
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_REGULAR:
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_DEV2:
    nl_vec_exist = nl_vec_exist &&
                   stokes_problem_nl->proj_scal != NULL &&
                   stokes_problem_nl->proj_tens != NULL;
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL:
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL_SYMM:
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_DEV1:
    nl_vec_exist = nl_vec_exist &&
                   stokes_problem_nl->proj_scal != NULL;
    break;
  default: /* unknown linearization type */
    RHEA_ABORT_NOT_REACHED ();
  }

  return (mesh_exists && vectors_exist && weak_exists &&
          rhs_nonzero_dir_exists && rhs_nonzero_neu_exists && nl_vec_exist);
}

static void
rhea_stokes_problem_nonlinear_create_mesh_data (
                                    rhea_stokes_problem_t *stokes_problem_nl,
                                    ymir_mesh_t *ymir_mesh,
                                    ymir_pressure_elem_t *press_elem)
{
  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (
      !rhea_stokes_problem_nonlinear_mesh_data_exists (stokes_problem_nl));
  RHEA_ASSERT (stokes_problem_nl->domain_options != NULL);

  /* assign mesh objects */
  stokes_problem_nl->ymir_mesh = ymir_mesh;
  stokes_problem_nl->press_elem = press_elem;

  /* create velocity boundary conditions */
  rhea_stokes_problem_velocity_boundary_create (stokes_problem_nl);

  /* create coefficient and vectors related to it */
  if (stokes_problem_nl->weakzone_compute_fn != NULL) {
    stokes_problem_nl->weakzone = rhea_weakzone_new (ymir_mesh);
  }
  stokes_problem_nl->coeff = rhea_viscosity_new (ymir_mesh);
  stokes_problem_nl->bounds_marker = rhea_viscosity_new (ymir_mesh);
  stokes_problem_nl->yielding_marker = rhea_viscosity_new (ymir_mesh);

  /* create right-hand side vectors */
  stokes_problem_nl->rhs_vel = rhea_velocity_new (ymir_mesh);
  stokes_problem_nl->rhs_vel_compute_fn (
      stokes_problem_nl->rhs_vel,
      stokes_problem_nl->temperature,
      stokes_problem_nl->rhs_vel_compute_fn_data);
  if (stokes_problem_nl->rhs_vel_nonzero_dir_compute_fn != NULL) {
    stokes_problem_nl->rhs_vel_nonzero_dirichlet = rhea_velocity_new (
        ymir_mesh);
    stokes_problem_nl->rhs_vel_nonzero_dir_compute_fn (
        stokes_problem_nl->rhs_vel_nonzero_dirichlet,
        stokes_problem_nl->rhs_vel_nonzero_dir_compute_fn_data);
  }
  if (stokes_problem_nl->rhs_vel_nonzero_neu_compute_fn != NULL) {
    ymir_topidx_t       fm;

    stokes_problem_nl->rhs_vel_face_nonzero_neumann =
      RHEA_ALLOC (ymir_vec_t *, ymir_mesh->num_face_meshes);
    for (fm = 0; fm < ymir_mesh->num_face_meshes; fm++) { /* for all faces */
      stokes_problem_nl->rhs_vel_face_nonzero_neumann[fm] =
        ymir_face_cvec_new (ymir_mesh, fm, 3);
    }
    stokes_problem_nl->rhs_vel_nonzero_neu_compute_fn (
        stokes_problem_nl->rhs_vel_face_nonzero_neumann,
        stokes_problem_nl->rhs_vel_nonzero_neu_compute_fn_data);
  }
  stokes_problem_nl->rhs_vel_press = rhea_velocity_pressure_new (
      ymir_mesh, press_elem);

  /* create linearization data */
  stokes_problem_nl->sol_vel = rhea_velocity_new (ymir_mesh);
  switch (stokes_problem_nl->linearization_type) {
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_PICARD:
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_REGULAR:
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_DEV2:
    stokes_problem_nl->proj_scal = rhea_viscosity_new (ymir_mesh);
    stokes_problem_nl->proj_tens = rhea_strainrate_new (ymir_mesh);
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL:
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL_SYMM:
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_DEV1:
    stokes_problem_nl->proj_scal = rhea_viscosity_new (ymir_mesh);
    break;
  default: /* unknown linearization type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* create vectors for Newton problem */
  stokes_problem_nl->newton_neg_gradient_vec = rhea_velocity_pressure_new (
      ymir_mesh, press_elem);
  stokes_problem_nl->newton_step_vec = rhea_velocity_pressure_new (
      ymir_mesh, press_elem);

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static void
rhea_stokes_problem_nonlinear_clear_mesh_data (
                                    rhea_stokes_problem_t *stokes_problem_nl)
{
  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (
      rhea_stokes_problem_nonlinear_mesh_data_exists (stokes_problem_nl));
  RHEA_ASSERT (stokes_problem_nl->newton_problem != NULL);

  /* destroy vectors of Newton problem */
  rhea_velocity_pressure_destroy (
      stokes_problem_nl->newton_neg_gradient_vec);
  rhea_velocity_pressure_destroy (
      stokes_problem_nl->newton_step_vec);
  stokes_problem_nl->newton_neg_gradient_vec = NULL;
  stokes_problem_nl->newton_step_vec = NULL;
  rhea_newton_problem_set_vectors (
      NULL, NULL, stokes_problem_nl->newton_problem);

  /* destroy linearization data */
  rhea_velocity_destroy (stokes_problem_nl->sol_vel);
  switch (stokes_problem_nl->linearization_type) {
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_PICARD:
    RHEA_ASSERT (stokes_problem_nl->proj_scal == NULL);
    RHEA_ASSERT (stokes_problem_nl->proj_tens == NULL);
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_REGULAR:
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_DEV2:
    RHEA_ASSERT (stokes_problem_nl->proj_scal != NULL);
    RHEA_ASSERT (stokes_problem_nl->proj_tens != NULL);
    rhea_viscosity_destroy (stokes_problem_nl->proj_scal);
    rhea_strainrate_destroy (stokes_problem_nl->proj_tens);
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL:
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL_SYMM:
  case RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_DEV1:
    RHEA_ASSERT (stokes_problem_nl->proj_scal != NULL);
    RHEA_ASSERT (stokes_problem_nl->proj_tens == NULL);
    rhea_viscosity_destroy (stokes_problem_nl->proj_scal);
    break;
  default: /* unknown linearization type */
    RHEA_ABORT_NOT_REACHED ();
  }
  stokes_problem_nl->sol_vel = NULL;
  stokes_problem_nl->proj_scal = NULL;
  stokes_problem_nl->proj_tens = NULL;

  /* destroy vectors */
  rhea_viscosity_destroy (stokes_problem_nl->coeff);
  rhea_viscosity_destroy (stokes_problem_nl->bounds_marker);
  rhea_viscosity_destroy (stokes_problem_nl->yielding_marker);
  rhea_velocity_destroy (stokes_problem_nl->rhs_vel);
  rhea_velocity_pressure_destroy (stokes_problem_nl->rhs_vel_press);
  stokes_problem_nl->coeff = NULL;
  stokes_problem_nl->bounds_marker = NULL;
  stokes_problem_nl->yielding_marker = NULL;
  stokes_problem_nl->rhs_vel = NULL;
  stokes_problem_nl->rhs_vel_press = NULL;
  if (stokes_problem_nl->weakzone != NULL) {
    rhea_weakzone_destroy (stokes_problem_nl->weakzone);
    stokes_problem_nl->weakzone = NULL;
  }
  if (stokes_problem_nl->rhs_vel_nonzero_dirichlet != NULL) {
    rhea_velocity_destroy (stokes_problem_nl->rhs_vel_nonzero_dirichlet);
    stokes_problem_nl->rhs_vel_nonzero_dirichlet = NULL;
  }
  if (stokes_problem_nl->rhs_vel_face_nonzero_neumann != NULL) {
    ymir_mesh_t        *ymir_mesh = stokes_problem_nl->ymir_mesh;
    ymir_topidx_t       fm;

    for (fm = 0; fm < ymir_mesh->num_face_meshes; fm++) { /* for all faces */
      ymir_vec_destroy (stokes_problem_nl->rhs_vel_face_nonzero_neumann[fm]);
    }
    RHEA_FREE (stokes_problem_nl->rhs_vel_face_nonzero_neumann);
    stokes_problem_nl->rhs_vel_face_nonzero_neumann = NULL;
  }

  /* destroy velocity boundary conditions */
  rhea_stokes_problem_velocity_boundary_clear (stokes_problem_nl);

  /* remove vector references */
  rhea_stokes_problem_remove_temperature (stokes_problem_nl);
  rhea_stokes_problem_remove_velocity_pressure (stokes_problem_nl);

  /* remove mesh references */
  stokes_problem_nl->ymir_mesh = NULL;
  stokes_problem_nl->press_elem = NULL;

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static rhea_stokes_problem_t *
rhea_stokes_problem_nonlinear_new (ymir_mesh_t *ymir_mesh,
                                   ymir_pressure_elem_t *press_elem,
                                   ymir_vec_t *temperature,
                                   rhea_domain_options_t *domain_options,
                                   rhea_temperature_options_t *temp_options,
                                   rhea_plate_options_t *plate_options,
                                   rhea_weakzone_options_t *weak_options,
                                   rhea_viscosity_options_t *visc_options)
{
  rhea_stokes_problem_t *stokes_problem_nl;

  RHEA_GLOBAL_PRODUCTIONF_FN_BEGIN (
      __func__, "linearization=%i",
      rhea_stokes_problem_nonlinear_linearization_type);

  /* check input */
  RHEA_ASSERT (visc_options->type == RHEA_VISCOSITY_NONLINEAR);

  /* initialize general Stokes problem structure */
  stokes_problem_nl = rhea_stokes_problem_struct_new (
      RHEA_STOKES_PROBLEM_NONLINEAR,
      domain_options, temp_options, plate_options, weak_options, visc_options);

  /* set temperature */
  rhea_stokes_problem_set_temperature (stokes_problem_nl, temperature);

  /* initialize nonlinear structure */
  stokes_problem_nl->linearization_type =
    (rhea_stokes_problem_nonlinear_linearization_t)
    rhea_stokes_problem_nonlinear_linearization_type;
  stokes_problem_nl->newton_options = &rhea_stokes_problem_newton_options;
  stokes_problem_nl->norm_type =
    (rhea_stokes_norm_type_t) rhea_stokes_problem_nonlinear_norm_type;
  stokes_problem_nl->norm_op_mass_scaling =
    rhea_stokes_problem_nonlinear_norm_mass_scaling;

  /* create mesh data */
  rhea_stokes_problem_nonlinear_create_mesh_data (stokes_problem_nl, ymir_mesh,
                                                  press_elem);

  /* create Newton problem */
  {
    rhea_newton_problem_t *newton_problem;

    newton_problem = rhea_newton_problem_new (
        rhea_stokes_problem_nonlinear_compute_negative_gradient_fn,
        rhea_stokes_problem_nonlinear_compute_gradient_norm_fn,
        RHEA_STOKES_PROBLEM_NONLINEAR_GRAD_NORM_N_COMPONENTS,
        rhea_stokes_problem_nonlinear_solve_hessian_system_fn);
    stokes_problem_nl->newton_problem = newton_problem;

    RHEA_ASSERT (stokes_problem_nl->newton_neg_gradient_vec != NULL);
    RHEA_ASSERT (stokes_problem_nl->newton_step_vec != NULL);
    rhea_newton_problem_set_vectors (
        stokes_problem_nl->newton_neg_gradient_vec,
        stokes_problem_nl->newton_step_vec,
        newton_problem);

    rhea_newton_problem_set_data_fn (
        stokes_problem_nl,
        rhea_stokes_problem_nonlinear_create_solver_data_fn,
        NULL /* keep solver data for postprocessing */, newton_problem);
    rhea_newton_problem_set_apply_hessian_fn (
        rhea_stokes_problem_nonlinear_apply_hessian_fn, newton_problem);
    rhea_newton_problem_set_update_fn (
        rhea_stokes_problem_nonlinear_update_operator_fn,
        rhea_stokes_problem_nonlinear_update_hessian_fn,
        rhea_stokes_problem_nonlinear_modify_hessian_system_fn, newton_problem);
    rhea_newton_problem_set_setup_poststep_fn (
        rhea_stokes_problem_nonlinear_amr_fn, newton_problem);
    rhea_newton_problem_set_output_fn (
        rhea_stokes_problem_nonlinear_output_prestep_fn, newton_problem);

    if (rhea_stokes_problem_nonlinear_check_jacobian) {
      rhea_newton_problem_set_checks (
          0 /* grad */, 1 /* Hessian */, newton_problem);
    }
  }

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);

  return stokes_problem_nl;
}

static void
rhea_stokes_problem_nonlinear_destroy (
                                    rhea_stokes_problem_t *stokes_problem_nl)
{
  RHEA_GLOBAL_PRODUCTION_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->newton_problem != NULL);

  /* destroy solver data */
  if (rhea_stokes_problem_nonlinear_solver_data_exists (stokes_problem_nl)) {
    rhea_stokes_problem_nonlinear_clear_solver_data_fn (stokes_problem_nl);
  }

  /* destroy mesh data */
  rhea_stokes_problem_nonlinear_clear_mesh_data (stokes_problem_nl);

  /* destroy Newton problem */
  rhea_newton_problem_destroy (stokes_problem_nl->newton_problem);

  /* destroy problem structure */
  rhea_stokes_problem_struct_destroy (stokes_problem_nl);

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);
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
rhea_stokes_problem_nonlinear_update_solver (
                                    rhea_stokes_problem_t *stokes_problem_nl,
                                    const int update_coeff,
                                    ymir_vec_t *vel_press,
                                    const int override_rhs,
                                    ymir_vec_t *rhs_vel_press,
                                    ymir_vec_t *rhs_vel,
                                    ymir_vec_t *rhs_vel_nonzero_dirichlet,
                                    ymir_vec_t **rhs_vel_face_nonzero_neumann)
{
  const int           vel_press_exists = (vel_press != NULL);

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (
      __func__, "update_coeff=%i, vel_press=%i, override_rhs=%i",
      update_coeff, vel_press_exists, override_rhs);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (
      rhea_stokes_problem_nonlinear_mesh_data_exists (stokes_problem_nl));
  RHEA_ASSERT (
      rhea_stokes_problem_nonlinear_solver_data_exists (stokes_problem_nl));
  RHEA_ASSERT (stokes_problem_nl->stokes_pc != NULL);
  RHEA_ASSERT (rhs_vel_press == NULL ||
               rhea_velocity_pressure_check_vec_type (rhs_vel_press));
  RHEA_ASSERT (rhs_vel == NULL || rhea_velocity_check_vec_type (rhs_vel));
  RHEA_ASSERT (rhs_vel_nonzero_dirichlet == NULL ||
               rhea_velocity_check_vec_type (rhs_vel_nonzero_dirichlet));
#ifdef RHEA_ENABLE_DEBUG
  if (rhs_vel_face_nonzero_neumann != NULL) {
    ymir_mesh_t        *ymir_mesh = stokes_problem_nl->ymir_mesh;
    ymir_topidx_t       fm;

    for (fm = 0; fm < ymir_mesh->num_face_meshes; fm++) {
      RHEA_ASSERT (
          rhea_velocity_check_vec_type (rhs_vel_face_nonzero_neumann[fm]) &&
          ymir_vec_is_face_vec (rhs_vel_face_nonzero_neumann[fm]) &&
          rhs_vel_face_nonzero_neumann[fm]->meshnum == fm );
    }
  }
#endif

  /* update the Stokes coefficient and dependent objects */
  if (update_coeff) {
    /* compute weak zone */
    rhea_stokes_problem_weakzone_compute (stokes_problem_nl->weakzone,
                                          stokes_problem_nl);

    /* update Stokes operator and preconditioner */
    rhea_stokes_problem_nonlinear_update_operator_fn (vel_press,
                                                      stokes_problem_nl);
    rhea_stokes_problem_nonlinear_update_hessian_fn (vel_press, NULL, NAN,
                                                     stokes_problem_nl);
  }

  /* re-construct right-hand side for Stokes system */
  if (override_rhs) {
    RHEA_GLOBAL_VERBOSEF_FN_TAG (
        __func__,
        "rhs_vel_press=%i, rhs_vel=%i, rhs_vel_nz_dir=%i, rhs_vel_nz_neu=%i",
        (rhs_vel_press != NULL), (rhs_vel != NULL),
        (rhs_vel_nonzero_dirichlet != NULL),
        (rhs_vel_face_nonzero_neumann != NULL));
    if (rhs_vel_press != NULL) {
      ymir_vec_copy (rhs_vel_press, stokes_problem_nl->rhs_vel_press);
    }
    else {
      ymir_stokes_pc_construct_rhs (
          stokes_problem_nl->rhs_vel_press, /* output: right-hand side */
          rhs_vel,                      /* input: volume forcing */
          rhs_vel_face_nonzero_neumann, /* input: nonzero Neumann forcing */
          rhs_vel_nonzero_dirichlet,    /* input: nonzero Dirichlet boundary */
          stokes_problem_nl->incompressible,
          stokes_problem_nl->stokes_op, 1 /* linearized */);
    }
  }
  else {
    rhs_vel = stokes_problem_nl->rhs_vel;
    rhs_vel_nonzero_dirichlet =
      stokes_problem_nl->rhs_vel_nonzero_dirichlet;
    rhs_vel_face_nonzero_neumann =
      stokes_problem_nl->rhs_vel_face_nonzero_neumann;
    ymir_stokes_pc_construct_rhs (
        stokes_problem_nl->rhs_vel_press, /* output: right-hand side */
        rhs_vel,                      /* input: volume forcing */
        rhs_vel_face_nonzero_neumann, /* input: nonzero Neumann forcing */
        rhs_vel_nonzero_dirichlet,    /* input: nonzero Dirichlet boundary */
        stokes_problem_nl->incompressible,
        stokes_problem_nl->stokes_op, 1 /* linearized */);
  }
  RHEA_ASSERT (
      rhea_velocity_pressure_is_valid (stokes_problem_nl->rhs_vel_press));

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static int
rhea_stokes_problem_nonlinear_solve (ymir_vec_t **sol_vel_press,
                                     const int nonzero_initial_guess,
                                     const int iter_max,
                                     const double rtol,
                                     rhea_stokes_problem_t *stokes_problem_nl,
                                     int *num_iterations,
                                     double *residual_reduction)
{
  const int           status_verbosity = 1;
  rhea_newton_options_t *newton_options = stokes_problem_nl->newton_options;
  rhea_newton_problem_t *newton_problem = stokes_problem_nl->newton_problem;
  int                 stop_reason;

  RHEA_GLOBAL_PRODUCTION_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (*sol_vel_press));

  /* initialize flag that tracks whether AMR was performed during solve */
  stokes_problem_nl->mesh_modified_by_solver = 0;

  /* run Newton solver */
  newton_options->nonzero_initial_guess = nonzero_initial_guess;
  newton_options->iter_max = newton_options->iter_start + iter_max;
  newton_options->rtol = rtol;
  newton_options->status_verbosity = status_verbosity;
  stop_reason = rhea_newton_solve (sol_vel_press, newton_problem,
                                   newton_options);

  /* project out null spaces in solution */
  rhea_stokes_problem_project_out_nullspace (*sol_vel_press, stokes_problem_nl);

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);

  /* return iterations count, residual reduction, and stopping reason */
  if (NULL != num_iterations) {
    *num_iterations = rhea_newton_solve_get_num_iterations (newton_problem);
  }
  if (NULL != residual_reduction) {
    *residual_reduction =
      rhea_newton_solve_get_residual_reduction (newton_problem);
  }
  return stop_reason;
}

/******************************************************************************
 * Stokes Problem
 *****************************************************************************/

rhea_stokes_problem_t *
rhea_stokes_problem_new (ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         ymir_vec_t *temperature,
                         rhea_domain_options_t *domain_options,
                         rhea_temperature_options_t *temp_options,
                         rhea_weakzone_options_t *weak_options,
                         rhea_viscosity_options_t *visc_options)
{
  switch (visc_options->type) {
  case RHEA_VISCOSITY_LINEAR:
    return rhea_stokes_problem_linear_new (
        ymir_mesh, press_elem, temperature, domain_options, temp_options,
        NULL /* plate_options */, weak_options, visc_options);
  case RHEA_VISCOSITY_NONLINEAR:
    return rhea_stokes_problem_nonlinear_new (
        ymir_mesh, press_elem, temperature, domain_options, temp_options,
        NULL /* plate_options */, weak_options, visc_options);
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
rhea_stokes_problem_create_mesh_dependencies (
                                        rhea_stokes_problem_t *stokes_problem,
                                        ymir_mesh_t *ymir_mesh,
                                        ymir_pressure_elem_t *press_elem)
{
  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);

  switch (stokes_problem->type) {
  case RHEA_STOKES_PROBLEM_LINEAR:
    rhea_stokes_problem_linear_create_mesh_data (
        stokes_problem, ymir_mesh, press_elem);
    if (stokes_problem->recreate_solver_data) {
      rhea_stokes_problem_linear_create_solver_data (stokes_problem);
    }
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR:
    rhea_stokes_problem_nonlinear_create_mesh_data (
        stokes_problem, ymir_mesh, press_elem);
    RHEA_ASSERT (stokes_problem->newton_neg_gradient_vec != NULL);
    RHEA_ASSERT (stokes_problem->newton_step_vec != NULL);
    RHEA_ASSERT (stokes_problem->newton_problem != NULL);
    rhea_newton_problem_set_vectors (
        stokes_problem->newton_neg_gradient_vec,
        stokes_problem->newton_step_vec,
        stokes_problem->newton_problem);
    if (stokes_problem->recreate_solver_data) {
      rhea_stokes_problem_nonlinear_create_solver_data_fn (
          stokes_problem->velocity_pressure, stokes_problem);
    }
    break;
  default: /* unknown Stokes type */
    RHEA_ABORT_NOT_REACHED ();
  }

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

void
rhea_stokes_problem_clear_mesh_dependencies (
                                        rhea_stokes_problem_t *stokes_problem)
{
  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);

  switch (stokes_problem->type) {
  case RHEA_STOKES_PROBLEM_LINEAR:
    if (rhea_stokes_problem_linear_solver_data_exists (stokes_problem)) {
      rhea_stokes_problem_linear_clear_solver_data (stokes_problem);
      stokes_problem->recreate_solver_data = 1;
    }
    else {
      stokes_problem->recreate_solver_data = 0;
    }
    rhea_stokes_problem_linear_clear_mesh_data (stokes_problem);
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR:
    if (rhea_stokes_problem_nonlinear_solver_data_exists (stokes_problem)) {
      rhea_stokes_problem_nonlinear_clear_solver_data_fn (stokes_problem);
      stokes_problem->recreate_solver_data = 1;
    }
    else {
      stokes_problem->recreate_solver_data = 0;
    }
    rhea_stokes_problem_nonlinear_clear_mesh_data (stokes_problem);
    break;
  default: /* unknown Stokes type */
    RHEA_ABORT_NOT_REACHED ();
  }

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

void
rhea_stokes_problem_setup_solver_ext (rhea_stokes_problem_t *stokes_problem,
                                      const int num_krylov_solvers)
{
  if (1 < num_krylov_solvers) {
    stokes_problem->num_krylov_solvers = num_krylov_solvers;
  }

  switch (stokes_problem->type) {
  case RHEA_STOKES_PROBLEM_LINEAR:
    if (rhea_stokes_problem_linear_solver_data_exists (stokes_problem) &&
        stokes_problem->stokes_pc != NULL) {
      rhea_stokes_problem_linear_update_solver (
          stokes_problem, 1 /* update_coeff */,
          0 /* !override_rhs */, NULL, NULL, NULL, NULL);
    }
    else {
      rhea_stokes_problem_linear_setup_solver (stokes_problem);
    }
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR:
    rhea_stokes_problem_nonlinear_setup_solver (stokes_problem);
    break;
  default: /* unknown Stokes type */
    RHEA_ABORT_NOT_REACHED ();
  }
}

void
rhea_stokes_problem_setup_solver (rhea_stokes_problem_t *stokes_problem)
{
  rhea_stokes_problem_setup_solver_ext (stokes_problem, -1 /* [ignored] */);
}

void
rhea_stokes_problem_update_solver (rhea_stokes_problem_t *stokes_problem,
                                   const int update_coeff,
                                   ymir_vec_t *vel_press,
                                   const int override_rhs,
                                   ymir_vec_t *rhs_vel_press,
                                   ymir_vec_t *rhs_vel,
                                   ymir_vec_t *rhs_vel_nonzero_dirichlet,
                                   ymir_vec_t **rhs_vel_face_nonzero_neumann)
{
  switch (stokes_problem->type) {
  case RHEA_STOKES_PROBLEM_LINEAR:
    rhea_stokes_problem_linear_update_solver (
        stokes_problem, update_coeff, override_rhs, rhs_vel_press, rhs_vel,
        rhs_vel_nonzero_dirichlet, rhs_vel_face_nonzero_neumann);
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR:
    rhea_stokes_problem_nonlinear_update_solver (
        stokes_problem, update_coeff, vel_press, override_rhs, rhs_vel_press,
        rhs_vel, rhs_vel_nonzero_dirichlet, rhs_vel_face_nonzero_neumann);
    break;
  default: /* unknown Stokes type */
    RHEA_ABORT_NOT_REACHED ();
  }
}

int
rhea_stokes_problem_solve_ext (ymir_vec_t **sol_vel_press,
                               const int nonzero_initial_guess,
                               const int iter_max,
                               const double rtol,
                               rhea_stokes_problem_t *stokes_problem,
                               const int force_linear_solve,
                               const int krylov_solver_idx,
                               int *num_iterations,
                               double *residual_reduction,
                               int *mesh_modified_by_solver)
{
  int                 stop_reason;

  switch (stokes_problem->type) {
  case RHEA_STOKES_PROBLEM_LINEAR:
    stop_reason = rhea_stokes_problem_linear_solve (
        sol_vel_press, nonzero_initial_guess, iter_max, rtol, stokes_problem,
        krylov_solver_idx, num_iterations, residual_reduction);
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR:
    if (force_linear_solve) { /* if run only a single linear solve */
      const int           check_residual = (NULL != residual_reduction);
      double              res_reduction[3];

      stop_reason = rhea_stokes_problem_run_krylov_solver (
          *sol_vel_press, stokes_problem->rhs_vel_press,
          nonzero_initial_guess, iter_max, rtol,
          stokes_problem, krylov_solver_idx,
          num_iterations, (check_residual ? res_reduction : NULL));
      if (check_residual) {
        *residual_reduction = res_reduction[0];
      }
    }
    else { /* otherwise run full nonlinear solve */
      stop_reason = rhea_stokes_problem_nonlinear_solve (
          sol_vel_press, nonzero_initial_guess, iter_max, rtol,
          stokes_problem, num_iterations, residual_reduction);
    }
    break;
  default: /* unknown Stokes type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* return whether AMR was performed during solve */
  if (mesh_modified_by_solver != NULL) {
    *mesh_modified_by_solver = stokes_problem->mesh_modified_by_solver;
  }

  /* return stopping reason */
  return stop_reason;
}

int
rhea_stokes_problem_has_converged_ext (const int stop_reason,
                                       rhea_stokes_problem_t *stokes_problem,
                                       const int force_linear_solve)
{
  switch (stokes_problem->type) {
  case RHEA_STOKES_PROBLEM_LINEAR:
    return (2 == stop_reason || 3 == stop_reason);
  case RHEA_STOKES_PROBLEM_NONLINEAR:
    if (force_linear_solve) {
      return (2 == stop_reason || 3 == stop_reason);
    }
    else {
      return rhea_newton_solve_has_converged (stop_reason);
    }
  default: /* unknown Stokes type */
    RHEA_ABORT_NOT_REACHED ();
  }
}

int
rhea_stokes_problem_solve (ymir_vec_t **sol_vel_press,
                           const int nonzero_initial_guess,
                           const int iter_max,
                           const double rtol,
                           rhea_stokes_problem_t *stokes_problem)
{
  return rhea_stokes_problem_solve_ext (
      sol_vel_press, nonzero_initial_guess, iter_max, rtol, stokes_problem,
      0 /* !force_linear_solve */, 0 /* krylov_solver_idx */,
      NULL /* num_iter */, NULL /* res_reduc */,
      NULL /* mesh_modified_by_solver */);
}

int
rhea_stokes_problem_has_converged (const int stop_reason,
                                   rhea_stokes_problem_t *stokes_problem)
{
  return rhea_stokes_problem_has_converged_ext (stop_reason, stokes_problem,
                                                0 /* !force_linear_solve */);
}

void
rhea_stokes_problem_set_solver_amr (
                                rhea_stokes_problem_t *stokes_problem,
                                p4est_t *p4est,
                                rhea_discretization_options_t *discr_options)
{
  stokes_problem->p4est = p4est;
  stokes_problem->discr_options = discr_options;
}

void
rhea_stokes_problem_set_solver_bin_output (
                                    rhea_stokes_problem_t *stokes_problem,
                                    char *bin_path)
{
  stokes_problem->solver_bin_path = bin_path;
}

void
rhea_stokes_problem_set_solver_vtk_output (
                                    rhea_stokes_problem_t *stokes_problem,
                                    char *vtk_path)
{
  stokes_problem->solver_vtk_path = vtk_path;
}

/******************************************************************************
 * Data Access
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

void
rhea_stokes_problem_set_temperature (rhea_stokes_problem_t *stokes_problem,
                                     ymir_vec_t *temperature)
{
  stokes_problem->temperature = temperature;
}

ymir_vec_t *
rhea_stokes_problem_get_temperature (rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->temperature;
}

void
rhea_stokes_problem_remove_temperature (rhea_stokes_problem_t *stokes_problem)
{
  rhea_stokes_problem_set_temperature (stokes_problem, NULL);
}

void
rhea_stokes_problem_set_velocity_pressure (rhea_stokes_problem_t *stokes_problem,
                                           ymir_vec_t *velocity_pressure)
{
  stokes_problem->velocity_pressure = velocity_pressure;
}

ymir_vec_t *
rhea_stokes_problem_get_velocity_pressure (rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->velocity_pressure;
}

void
rhea_stokes_problem_remove_velocity_pressure (
                                        rhea_stokes_problem_t *stokes_problem)
{
  rhea_stokes_problem_set_velocity_pressure (stokes_problem, NULL);
}

int
rhea_stokes_problem_is_nonlinear (rhea_stokes_problem_t *stokes_problem)
{
  switch (stokes_problem->type) {
  case RHEA_STOKES_PROBLEM_LINEAR:    return 0;
  case RHEA_STOKES_PROBLEM_NONLINEAR: return 1;
  default: /* unknown Stokes type */
    RHEA_ABORT_NOT_REACHED ();
    return -1;
  }
}

rhea_domain_options_t *
rhea_stokes_problem_get_domain_options (rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->domain_options;
}

void
rhea_stokes_problem_set_domain_options (rhea_stokes_problem_t *stokes_problem,
                                        rhea_domain_options_t *domain_options)
{
  stokes_problem->domain_options = domain_options;
}

rhea_temperature_options_t *
rhea_stokes_problem_get_temperature_options (
                                        rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->temp_options;
}

void
rhea_stokes_problem_set_temperature_options (
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_temperature_options_t *temp_options)
{
  stokes_problem->temp_options = temp_options;
}

rhea_plate_options_t *
rhea_stokes_problem_get_plate_options (rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->plate_options;
}

void
rhea_stokes_problem_set_plate_options (rhea_stokes_problem_t *stokes_problem,
                                       rhea_plate_options_t *plate_options)
{
  stokes_problem->plate_options = plate_options;
}

rhea_weakzone_options_t *
rhea_stokes_problem_get_weakzone_options (
                                        rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->weak_options;
}

void
rhea_stokes_problem_set_weakzone_options (
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_weakzone_options_t *weak_options)
{
  stokes_problem->weak_options = weak_options;
}

rhea_viscosity_options_t *
rhea_stokes_problem_get_viscosity_options (
                                        rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->visc_options;
}

void
rhea_stokes_problem_set_viscosity_options (
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_viscosity_options_t *visc_options)
{
  stokes_problem->visc_options = visc_options;
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

ymir_vec_t *
rhea_stokes_problem_get_weakzone (rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->weakzone;
}

//TODO check wheter this is necessary
ymir_vec_t *
rhea_stokes_problem_get_rhs_vel_press (rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->rhs_vel_press;
}//XI

//TODO check wheter this is necessary
void
rhea_stokes_problem_set_rhs_vel_press (
                                  rhea_stokes_problem_t *stokes_problem,
                                  ymir_vec_t *rhs_vel_press)
{
  /* check input */
  RHEA_ASSERT (rhs_vel_press != NULL);

  ymir_vec_copy (rhs_vel_press, stokes_problem->rhs_vel_press);
}//XI for hessian forward rhs.

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

ymir_vec_t **
rhea_stokes_problem_get_rhs_vel_nonzero_neumann (
                                        rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->rhs_vel_face_nonzero_neumann;
}

ymir_vec_t *
rhea_stokes_problem_get_rhs_vel_nonzero_neumann_surface (
                                        rhea_stokes_problem_t *stokes_problem)
{
  const ymir_topidx_t face_surf = RHEA_DOMAIN_BOUNDARY_FACE_TOP;

  if (stokes_problem->rhs_vel_face_nonzero_neumann != NULL) {
    return stokes_problem->rhs_vel_face_nonzero_neumann[face_surf];
  }
  else {
    return NULL;
  }
}

void
rhea_stokes_problem_set_weakzone_compute_fn (
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_weakzone_compute_fn_t fn,
                                    void *data)
{
  /* set function and corresponding data */
  stokes_problem->weakzone_compute_fn = fn;
  stokes_problem->weakzone_compute_fn_data = data;

  /* create or destroy weak zone vector */
  if (NULL != stokes_problem->weakzone_compute_fn &&
      NULL == stokes_problem->weakzone) { /* if vector does not exist */
      RHEA_ASSERT (stokes_problem->ymir_mesh != NULL);
      stokes_problem->weakzone = rhea_weakzone_new (stokes_problem->ymir_mesh);
  }
  else if (NULL == stokes_problem->weakzone_compute_fn &&
           NULL != stokes_problem->weakzone) { /* if vector is not required */
      rhea_weakzone_destroy (stokes_problem->weakzone);
      stokes_problem->weakzone = NULL;
  }
}

void
rhea_stokes_problem_set_viscosity_compute_fn (
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_viscosity_compute_fn_t fn,
                                    void *data)
{
  /* check input */
  RHEA_ASSERT (fn != NULL);

  /* set function and corresponding data */
  stokes_problem->viscosity_compute_fn = fn;
  stokes_problem->viscosity_compute_fn_data = data;
}

void
rhea_stokes_problem_set_rhs_vel_compute_fn (
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_velocity_rhs_compute_fn_t fn,
                                    void *data)
{
  /* check input */
  RHEA_ASSERT (stokes_problem->rhs_vel != NULL);
  RHEA_ASSERT (fn != NULL);

  /* set function and corresponding data */
  stokes_problem->rhs_vel_compute_fn = fn;
  stokes_problem->rhs_vel_compute_fn_data = data;

  /* recompute */
  stokes_problem->rhs_vel_compute_fn (
      stokes_problem->rhs_vel,
      stokes_problem->temperature,
      stokes_problem->rhs_vel_compute_fn_data);
}

void
rhea_stokes_problem_set_rhs_vel_nonzero_dir_compute_fn (
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_velocity_rhs_nz_dir_compute_fn_t fn,
                                    void *data)
{
  /* check input */
  RHEA_ASSERT (stokes_problem->ymir_mesh != NULL);

  /* set function and corresponding data */
  stokes_problem->rhs_vel_nonzero_dir_compute_fn = fn;
  stokes_problem->rhs_vel_nonzero_dir_compute_fn_data = data;

  /* recompute */
  if (stokes_problem->rhs_vel_nonzero_dir_compute_fn != NULL) { /* if fn. */
    if (stokes_problem->rhs_vel_nonzero_dirichlet == NULL) {
      stokes_problem->rhs_vel_nonzero_dirichlet =
        rhea_velocity_new (stokes_problem->ymir_mesh);
    }
    stokes_problem->rhs_vel_nonzero_dir_compute_fn (
        stokes_problem->rhs_vel_nonzero_dirichlet,
        stokes_problem->rhs_vel_nonzero_dir_compute_fn_data);
  }
  else { /* if no function provided */
    if (stokes_problem->rhs_vel_nonzero_dirichlet != NULL) {
      rhea_velocity_destroy (stokes_problem->rhs_vel_nonzero_dirichlet);
      stokes_problem->rhs_vel_nonzero_dirichlet = NULL;
    }
  }
}

void
rhea_stokes_problem_set_rhs_vel_nonzero_neu_compute_fn (
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_velocity_rhs_nz_neu_compute_fn_t fn,
                                    void *data)
{
  ymir_mesh_t        *ymir_mesh = stokes_problem->ymir_mesh;
  ymir_topidx_t       fm;

  /* check input */
  RHEA_ASSERT (stokes_problem->ymir_mesh != NULL);

  /* set function and corresponding data */
  stokes_problem->rhs_vel_nonzero_neu_compute_fn = fn;
  stokes_problem->rhs_vel_nonzero_neu_compute_fn_data = data;

  /* recompute */
  if (stokes_problem->rhs_vel_nonzero_neu_compute_fn != NULL) { /* if fn. */
    if (stokes_problem->rhs_vel_face_nonzero_neumann == NULL) {
      stokes_problem->rhs_vel_face_nonzero_neumann =
        RHEA_ALLOC (ymir_vec_t *, ymir_mesh->num_face_meshes);
      for (fm = 0; fm < ymir_mesh->num_face_meshes; fm++) { /* for all faces */
        stokes_problem->rhs_vel_face_nonzero_neumann[fm] =
          ymir_face_cvec_new (ymir_mesh, fm, 3);
      }
    }
    stokes_problem->rhs_vel_nonzero_neu_compute_fn (
        stokes_problem->rhs_vel_face_nonzero_neumann,
        stokes_problem->rhs_vel_nonzero_neu_compute_fn_data);
  }
  else { /* if no function provided */
    if (stokes_problem->rhs_vel_face_nonzero_neumann != NULL) {
      for (fm = 0; fm < ymir_mesh->num_face_meshes; fm++) { /* for all faces */
        ymir_vec_destroy (stokes_problem->rhs_vel_face_nonzero_neumann[fm]);
      }
      RHEA_FREE (stokes_problem->rhs_vel_face_nonzero_neumann);
      stokes_problem->rhs_vel_face_nonzero_neumann = NULL;
    }
  }
}

void
rhea_stokes_problem_weakzone_compute (ymir_vec_t *weakzone,
                                      rhea_stokes_problem_t *stokes_problem)
{
  if (NULL != stokes_problem->weakzone_compute_fn) { /* if fn. exists */
    if (NULL != weakzone) { /* if use vector from fn. input */
      stokes_problem->weakzone_compute_fn (
          weakzone, stokes_problem->weakzone_compute_fn_data);
    }
    else if (NULL != stokes_problem->weakzone) { /* if  use internal vector */
      stokes_problem->weakzone_compute_fn (
          stokes_problem->weakzone, stokes_problem->weakzone_compute_fn_data);
    }
    else { /* unkown output vector */
      RHEA_ABORT_NOT_REACHED ();
    }
  }
  else if (weakzone != NULL) { /* if no fn., then set neutral value */
    ymir_vec_set_value (weakzone, RHEA_WEAKZONE_NEUTRAL_VALUE);
  }
}

void
rhea_stokes_problem_viscosity_compute (ymir_vec_t *viscosity,
                                       ymir_vec_t *proj_scal,
                                       ymir_vec_t *bounds_marker,
                                       ymir_vec_t *yielding_marker,
                                       ymir_vec_t *temperature,
                                       ymir_vec_t *weakzone,
                                       ymir_vec_t *velocity,
                                       rhea_stokes_problem_t *stokes_problem)
{
  RHEA_ASSERT (stokes_problem->viscosity_compute_fn != NULL);
  stokes_problem->viscosity_compute_fn (
      /* out: */ viscosity, proj_scal, bounds_marker, yielding_marker,
      /* in:  */ temperature, weakzone, velocity,
      stokes_problem->viscosity_compute_fn_data);
}

ymir_stokes_op_t *
rhea_stokes_problem_get_stokes_op (rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->stokes_op;
}

int
rhea_stokes_problem_get_num_krylov_solvers (
                                        rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->num_krylov_solvers;
}

/******************************************************************************
 * I/O
 *****************************************************************************/

int
rhea_stokes_problem_write (char *base_path_bin,
                           rhea_stokes_problem_t *stokes_problem)
{
  ymir_mesh_t          *ymir_mesh = rhea_stokes_problem_get_ymir_mesh (
                            stokes_problem);
  ymir_pressure_elem_t *press_elem = rhea_stokes_problem_get_press_elem (
                            stokes_problem);
  ymir_vec_t         *temperature = rhea_stokes_problem_get_temperature (
                          stokes_problem);
  ymir_vec_t         *vel_press = rhea_stokes_problem_get_velocity_pressure (
                          stokes_problem);
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  int                 success = 0;

  /* write p4est */
  if (stokes_problem->p4est != NULL) {
    char                path[BUFSIZ];

    snprintf (path, BUFSIZ, "%s%s", base_path_bin,
              RHEA_STOKES_PROBLEM_IO_LABEL_P4EST);
    p4est_save_ext (path, stokes_problem->p4est, 0 /* do not save data */,
                    1 /* save partition dependent */);
    success++;
  }

  /* write temperature */
  if (temperature != NULL) {
    char                path[BUFSIZ];

    snprintf (path, BUFSIZ, "%s%s.bin", base_path_bin,
              RHEA_STOKES_PROBLEM_IO_LABEL_TEMPERATURE);
    success += rhea_temperature_write (path, temperature, mpicomm);
  }

  /* write velocity & pressure */
  if (vel_press != NULL) {
    char                vel_path[BUFSIZ], press_path[BUFSIZ];

    snprintf (vel_path, BUFSIZ, "%s%s.bin", base_path_bin,
              RHEA_STOKES_PROBLEM_IO_LABEL_VELOCITY);
    snprintf (press_path, BUFSIZ, "%s%s.bin", base_path_bin,
              RHEA_STOKES_PROBLEM_IO_LABEL_PRESSURE);
    success += rhea_velocity_pressure_write (vel_path, press_path, vel_press,
                                             press_elem, mpicomm);
  }

  return success;
}

/******************************************************************************
 * Vector Operations
 *****************************************************************************/

int
rhea_stokes_problem_velocity_boundary_set_zero (
                                        ymir_vec_t *velocity,
                                        rhea_stokes_problem_t *stokes_problem)
{
  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (velocity));

  /* set Dirichlet components to zero if Dirichlet BC's exist */
  if (stokes_problem->vel_dir != NULL) {
    ymir_vel_dir_separate (velocity, NULL, NULL, NULL, stokes_problem->vel_dir);
    return 1;
  }
  else {
    return 0;
  }
}

int
rhea_stokes_problem_velocity_compute_mean_rotation (
                                        double mean_rot_axis[3],
                                        ymir_vec_t *velocity,
                                        rhea_stokes_problem_t *stokes_problem)
{
  ymir_stress_op_t   *stress_op;

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (velocity));

  /* compute mean rotation */
  switch (stokes_problem->domain_options->shape) {
  case RHEA_DOMAIN_SHELL:
    RHEA_ASSERT (stokes_problem->stokes_op != NULL);
    stress_op = stokes_problem->stokes_op->stress_op;
    ymir_stress_op_compute_mean_rotation (mean_rot_axis, velocity,
                                          stress_op, 0);
    return 1;
  default: /* no mean rotation */
    mean_rot_axis[0] = NAN;
    mean_rot_axis[1] = NAN;
    mean_rot_axis[2] = NAN;
    return 0;
  }
}

int
rhea_stokes_problem_velocity_project_out_mean_rotation (
                                        ymir_vec_t *velocity,
                                        const int residual_space,
                                        rhea_stokes_problem_t *stokes_problem)
{
  ymir_stress_op_t   *stress_op;

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (velocity));

  /* compute mean rotation */
  switch (stokes_problem->domain_options->shape) {
  case RHEA_DOMAIN_SHELL:
    RHEA_ASSERT (stokes_problem->stokes_op != NULL);
    stress_op = stokes_problem->stokes_op->stress_op;
    ymir_stress_op_project_out_mean_rotation (velocity, stress_op,
                                              residual_space);
    return 1;
  default: /* no mean rotation */
    return 0;
  }
}

int
rhea_stokes_problem_project_out_nullspace (
                                        ymir_vec_t *vel_press,
                                        rhea_stokes_problem_t *stokes_problem)
{
  ymir_pressure_elem_t *press_elem = stokes_problem->press_elem;
  ymir_vec_t         *velocity, *pressure;
  int                 is_view;
  int                 proj_vel_rot, proj_press_mean;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (vel_press));

  /* get fields */
  is_view = rhea_velocity_pressure_create_components (&velocity, &pressure,
                                                      vel_press, press_elem);

  /* project out null spaces */
  proj_vel_rot = rhea_stokes_problem_velocity_project_out_mean_rotation (
      velocity, 0 /* !residual_space */, stokes_problem);
  proj_press_mean = rhea_pressure_project_out_mean (pressure, press_elem);

  /* set fields */
  if (!is_view) { /* if vector was copied */
    rhea_velocity_pressure_set_components (vel_press, velocity, pressure,
                                           press_elem);
  }

  /* destroy */
  ymir_vec_destroy (velocity);
  ymir_vec_destroy (pressure);

  RHEA_GLOBAL_VERBOSEF_FN_END (__func__, "proj_vel_rot=%i, proj_press_mean=%i",
                               proj_vel_rot, proj_press_mean);

  return (proj_vel_rot + proj_press_mean);
}

int
rhea_stokes_problem_stress_compute_normal_at_surface (
                                    ymir_vec_t *stress_norm_surf,
                                    ymir_vec_t *vel_press,
                                    rhea_stokes_problem_t *stokes_problem)
{
  ymir_mesh_t          *ymir_mesh = rhea_stokes_problem_get_ymir_mesh (
                            stokes_problem);
  ymir_pressure_elem_t *press_elem = rhea_stokes_problem_get_press_elem (
                            stokes_problem);

  int                 data_exists_lin, data_exists_nl;
  ymir_stokes_op_t   *stokes_op;
  ymir_vec_t         *rhs_vel = stokes_problem->rhs_vel;
  ymir_vec_t         *rhs_vel_nonzero_dirichlet =
                        stokes_problem->rhs_vel_nonzero_dirichlet;
  ymir_vec_t         *rhs_vel_press;
  ymir_vec_t         *vel_press_nspfree;
  ymir_vec_t         *residual_vel_press, *residual_vel;

  /* check input */
  RHEA_ASSERT (rhea_stress_surface_check_vec_type (stress_norm_surf));
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (vel_press));
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (vel_press));

  /* get Stokes operator */
  data_exists_lin = (
      stokes_problem->type == RHEA_STOKES_PROBLEM_LINEAR &&
      rhea_stokes_problem_linear_solver_data_exists (stokes_problem));
  data_exists_nl = (
      stokes_problem->type == RHEA_STOKES_PROBLEM_NONLINEAR &&
      rhea_stokes_problem_nonlinear_solver_data_exists (stokes_problem));
  if (data_exists_lin || data_exists_nl) {
    stokes_op = stokes_problem->stokes_op;
  }
  else {
    return 0;
  }

  /* construct right-hand side for Stokes system */
  rhs_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);
  ymir_stokes_pc_construct_rhs (
      rhs_vel_press,              /* output: right-hand side */
      rhs_vel,                    /* input: volume forcing */
      NULL,                       /* input: Neumann forcing */
      rhs_vel_nonzero_dirichlet,  /* input: nonzero Dirichlet boundary */
      stokes_problem->incompressible, stokes_op, 0 /* !linearized */);
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (rhs_vel_press));

  /* project out null spaces */
  vel_press_nspfree = rhea_velocity_pressure_new (ymir_mesh, press_elem);
  ymir_vec_copy (vel_press, vel_press_nspfree);
  rhea_stokes_problem_project_out_nullspace (vel_press_nspfree,
                                             stokes_problem);

  /* compute residual without constraining Dirichlet BC's
   *   r_mom  = A * u + B^T * p - f
   *   r_mass = B * u
   */
  residual_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);
  ymir_stokes_pc_apply_skip_dir_stokes_op (vel_press_nspfree,
                                           residual_vel_press, stokes_op,
                                           0 /* !linearized */, 1 /* dirty */);
  ymir_vec_add (-1.0, rhs_vel_press, residual_vel_press);
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (residual_vel_press));
  rhea_velocity_pressure_destroy (rhs_vel_press);
  rhea_velocity_pressure_destroy (vel_press_nspfree);

  /* communicate shared node values */
  ymir_vec_share_owned (residual_vel_press);

  /* get the velocity component of the residual */
  residual_vel = rhea_velocity_new (ymir_mesh);
  rhea_velocity_pressure_copy_components (residual_vel, NULL,
                                          residual_vel_press, press_elem);
  RHEA_ASSERT (rhea_velocity_is_valid (residual_vel));
  rhea_velocity_pressure_destroy (residual_vel_press);

  /* extract the normal component of the residual at the surface */
  rhea_stress_surface_extract_from_residual (stress_norm_surf, residual_vel);
  RHEA_ASSERT (rhea_velocity_surface_is_valid (stress_norm_surf));
  rhea_velocity_destroy (residual_vel);

  return 1;
}
