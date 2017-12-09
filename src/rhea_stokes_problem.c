/*
 */

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
#ifdef RHEA_ENABLE_DEBUG
# include <ymir_comm.h>
#endif

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

  /* state variables of the Stokes problem (not owned) */
  ymir_vec_t         *temperature;
  ymir_vec_t         *velocity_pressure;

  /* options (not owned) */
  rhea_domain_options_t      *domain_options;
  rhea_temperature_options_t *temp_options;
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
  ymir_vec_t         *rhs_vel_press;

  /* Stokes operator and preconditioner (solver data) */
  ymir_stokes_op_t   *stokes_op;
  ymir_stokes_pc_t   *stokes_pc;

  /* Newton problem (nonlinear Stokes only) */
  rhea_stokes_problem_nonlinear_linearization_t linearization_type;
  rhea_newton_options_t  *newton_options;
  rhea_newton_problem_t  *newton_problem;
  ymir_vec_t             *newton_neg_gradient_vec;
  ymir_vec_t             *newton_step_vec;
  rhea_stokes_norm_type_t norm_type;
  ymir_Hminus1_norm_op_t *norm_op;
  double                  norm_op_mass_scaling;

  /* solver VTK output path */
  char                   *solver_vtk_path;
};

/**
 * Creates and initializes a new structure of a Stokes problem.
 */
static rhea_stokes_problem_t *
rhea_stokes_problem_struct_new (const rhea_stokes_problem_type_t type,
                                rhea_domain_options_t *domain_options,
                                rhea_temperature_options_t *temp_options,
                                rhea_weakzone_options_t *weak_options,
                                rhea_viscosity_options_t *visc_options)
{
  rhea_stokes_problem_t *stokes_problem = RHEA_ALLOC (rhea_stokes_problem_t, 1);

  /* initialize */
  stokes_problem->type = type;
  stokes_problem->incompressible = 1;
  stokes_problem->ymir_mesh = NULL;
  stokes_problem->press_elem = NULL;

  stokes_problem->temperature = NULL;
  stokes_problem->velocity_pressure = NULL;

  stokes_problem->domain_options = domain_options;
  stokes_problem->temp_options = temp_options;
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
  stokes_problem->rhs_vel_press = NULL;

  stokes_problem->stokes_op = NULL;
  stokes_problem->stokes_pc = NULL;

  stokes_problem->linearization_type =
    RHEA_STOKES_PROBLEM_NONLINEAR_DEFAULT_LINEARIZATION_TYPE;
  stokes_problem->newton_options = NULL;
  stokes_problem->newton_problem = NULL;
  stokes_problem->newton_neg_gradient_vec = NULL;
  stokes_problem->newton_step_vec = NULL;
  stokes_problem->norm_type = RHEA_STOKES_NORM_NONE;
  stokes_problem->norm_op = NULL;
  stokes_problem->norm_op_mass_scaling = NAN;

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
  rhea_stokes_prbolem_set_viscosity_compute_fn (
      stokes_problem, rhea_viscosity_compute, visc_options);
  stokes_problem->rhs_vel_compute_fn = rhea_velocity_rhs_compute;
  stokes_problem->rhs_vel_compute_fn_data = temp_options;
  stokes_problem->rhs_vel_nonzero_dir_compute_fn = NULL;
  stokes_problem->rhs_vel_nonzero_dir_compute_fn_data = NULL;

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

/**
 * Computes viscosity.
 */
static void
rhea_stokes_problem_compute_coefficient (rhea_stokes_problem_t *stokes_problem,
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
    stokes_problem->viscosity_compute_fn (
        /* out: */ coeff, NULL, NULL, NULL,
        /* in:  */ temperature, weakzone, NULL,
        stokes_problem->viscosity_compute_fn_data);
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR:
    RHEA_ASSERT (visc_options->type == RHEA_VISCOSITY_NONLINEAR);
    RHEA_ASSERT (stokes_problem->bounds_marker != NULL);
    RHEA_ASSERT (stokes_problem->yielding_marker != NULL);
    if (nonlinear_init) {
      rhea_viscosity_compute_nonlinear_init (
          /* out: */ coeff, proj_scal, bounds_marker, yielding_marker,
          /* in:  */ temperature, weakzone,
          visc_options);
    }
    else {
      vel = rhea_stokes_problem_retrieve_velocity (
          stokes_problem->velocity_pressure, stokes_problem);
      stokes_problem->viscosity_compute_fn (
          /* out: */ coeff, proj_scal, bounds_marker, yielding_marker,
          /* in:  */ temperature, weakzone, vel,
          stokes_problem->viscosity_compute_fn_data);
    }
    break;
  default: /* unknown Stokes type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* transform viscosity to viscous stress coefficient */
  ymir_vec_scale (2.0, coeff);
}

/******************************************************************************
 * Linear Stokes Problem
 *****************************************************************************/

static void
rhea_stokes_problem_linear_create_mesh_data (
                                    rhea_stokes_problem_t *stokes_problem_lin,
                                    ymir_mesh_t *ymir_mesh,
                                    ymir_pressure_elem_t *press_elem)
{
  const char         *this_fn_name =
                        "rhea_stokes_problem_linear_create_mesh_data";

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);
  RHEA_ASSERT (stokes_problem_lin->ymir_mesh  == NULL);
  RHEA_ASSERT (stokes_problem_lin->press_elem == NULL);
  RHEA_ASSERT (stokes_problem_lin->coeff         == NULL);
  RHEA_ASSERT (stokes_problem_lin->weakzone      == NULL);
  RHEA_ASSERT (stokes_problem_lin->rhs_vel       == NULL);
  RHEA_ASSERT (stokes_problem_lin->rhs_vel_press == NULL);
  RHEA_ASSERT (stokes_problem_lin->domain_options != NULL);

  /* assign mesh objects */
  stokes_problem_lin->ymir_mesh = ymir_mesh;
  stokes_problem_lin->press_elem = press_elem;

  /* create velocity boundary conditions */
  rhea_stokes_problem_velocity_boundary_create (stokes_problem_lin);

  /* create coefficient and vectors related to it */
  if (stokes_problem_lin->weakzone_compute_fn != NULL) {
    stokes_problem_lin->weakzone = rhea_weakzone_new (ymir_mesh);
    stokes_problem_lin->weakzone_compute_fn (
        stokes_problem_lin->weakzone,
        stokes_problem_lin->weakzone_compute_fn_data);
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
  stokes_problem_lin->rhs_vel_press = rhea_velocity_pressure_new (
      ymir_mesh, press_elem);

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
}

static void
rhea_stokes_problem_linear_clear_mesh_data (
                                    rhea_stokes_problem_t *stokes_problem_lin)
{
  const char         *this_fn_name =
                        "rhea_stokes_problem_linear_clear_mesh_data";

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);
  RHEA_ASSERT (stokes_problem_lin->ymir_mesh  != NULL);
  RHEA_ASSERT (stokes_problem_lin->press_elem != NULL);
  RHEA_ASSERT (stokes_problem_lin->coeff         != NULL);
  RHEA_ASSERT (stokes_problem_lin->rhs_vel       != NULL);
  RHEA_ASSERT (stokes_problem_lin->rhs_vel_press != NULL);

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

  /* destroy velocity boundary conditions */
  rhea_stokes_problem_velocity_boundary_clear (stokes_problem_lin);

  /* remove vector references */
  rhea_stokes_problem_remove_temperature (stokes_problem_lin);
  rhea_stokes_problem_remove_velocity_pressure (stokes_problem_lin);

  /* remove mesh references */
  stokes_problem_lin->ymir_mesh = NULL;
  stokes_problem_lin->press_elem = NULL;

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
}

static void
rhea_stokes_problem_linear_create_solver_data (
                                    rhea_stokes_problem_t *stokes_problem_lin)
{
  const char         *this_fn_name =
                        "rhea_stokes_problem_linear_create_solver_data";
  const char         *vtk_path = stokes_problem_lin->solver_vtk_path;
  rhea_domain_options_t *domain_options = stokes_problem_lin->domain_options;

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->ymir_mesh  != NULL);
  RHEA_ASSERT (stokes_problem_lin->press_elem != NULL);
  RHEA_ASSERT (stokes_problem_lin->coeff   != NULL);
  RHEA_ASSERT (stokes_problem_lin->rhs_vel != NULL);
  RHEA_ASSERT (stokes_problem_lin->domain_options != NULL);
  RHEA_ASSERT (stokes_problem_lin->stokes_op == NULL);

  /* set VTK path for debuggin */
  if (vtk_path != NULL) {
    ymir_vtk_set_debug_path (vtk_path);
  }

  /* compute viscosity */
  rhea_stokes_problem_compute_coefficient (stokes_problem_lin, 0 /* unused */);

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

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
}

static void
rhea_stokes_problem_linear_clear_solver_data (
                                    rhea_stokes_problem_t *stokes_problem_lin)
{
  const char         *this_fn_name =
                        "rhea_stokes_problem_linear_clear_solver_data";

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", this_fn_name);

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

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
}

static rhea_stokes_problem_t *
rhea_stokes_problem_linear_new (ymir_mesh_t *ymir_mesh,
                                ymir_pressure_elem_t *press_elem,
                                ymir_vec_t *temperature,
                                rhea_domain_options_t *domain_options,
                                rhea_temperature_options_t *temp_options,
                                rhea_weakzone_options_t *weak_options,
                                rhea_viscosity_options_t *visc_options,
                                void *solver_options)
{
  const char         *this_fn_name = "rhea_stokes_problem_linear_new";
  rhea_stokes_problem_t *stokes_problem_lin;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (visc_options->type == RHEA_VISCOSITY_LINEAR);

  /* initialize general Stokes problem structure */
  stokes_problem_lin = rhea_stokes_problem_struct_new (
      RHEA_STOKES_PROBLEM_LINEAR,
      domain_options, temp_options, weak_options, visc_options);

  /* set temperature */
  rhea_stokes_problem_set_temperature (stokes_problem_lin, temperature);

  /* create mesh data */
  rhea_stokes_problem_linear_create_mesh_data (stokes_problem_lin, ymir_mesh,
                                               press_elem);

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

  /* destroy solver data */
  rhea_stokes_problem_linear_clear_solver_data (stokes_problem_lin);

  /* destroy mesh data */
  rhea_stokes_problem_linear_clear_mesh_data (stokes_problem_lin);

  /* destroy problem structure */
  rhea_stokes_problem_struct_destroy (stokes_problem_lin);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

static void
rhea_stokes_problem_linear_setup_solver (
                                    rhea_stokes_problem_t *stokes_problem_lin)
{
  const char         *this_fn_name = "rhea_stokes_problem_linear_setup_solver";
  ymir_stokes_op_t   *stokes_op;
  ymir_vec_t         *rhs_vel, *rhs_vel_nonzero_dirichlet;
  ymir_vec_t         *rhs_vel_press;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);
  RHEA_ASSERT (stokes_problem_lin->rhs_vel_press != NULL);
  RHEA_ASSERT (stokes_problem_lin->stokes_op == NULL);
  RHEA_ASSERT (stokes_problem_lin->stokes_pc == NULL);

  /* create solver data */
  rhea_stokes_problem_linear_create_solver_data (stokes_problem_lin);
  RHEA_ASSERT (stokes_problem_lin->stokes_op != NULL);

  /* construct right-hand side for Stokes system */
  stokes_op = stokes_problem_lin->stokes_op;
  rhs_vel = stokes_problem_lin->rhs_vel;
  rhs_vel_nonzero_dirichlet = stokes_problem_lin->rhs_vel_nonzero_dirichlet;
  rhs_vel_press = stokes_problem_lin->rhs_vel_press;
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

  /* setup solver if it has not been done yet */
  if (stokes_problem_lin->stokes_pc == NULL) {
    rhea_stokes_problem_linear_setup_solver (stokes_problem_lin);
  }

  RHEA_GLOBAL_PRODUCTIONF ("Into %s (max iter %i, rel tol %.3e)\n",
                           this_fn_name, iter_max, rel_tol);

  /* check input */
  RHEA_ASSERT (stokes_problem_lin->type == RHEA_STOKES_PROBLEM_LINEAR);
  RHEA_ASSERT (stokes_problem_lin->rhs_vel_press != NULL);
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
                                    stokes_problem_lin->stokes_op, 0, 1);
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
      RHEA_GLOBAL_PRODUCTIONF (
          "Done %s (solver stopping reason %i, iterations %i, "
          "residual reduction %.3e, convergence %.3e, "
          "reduction momentum %.3e, norm mass %.3e)\n",
          this_fn_name, stop_reason, itn, res_reduction, conv_factor,
          res_reduction_vel, norm_res_press);
    }
    else {
      RHEA_GLOBAL_PRODUCTIONF (
          "Done %s (solver stopping reason %i, iterations %i, "
          "residual reduction %.3e, convergence %.3e, "
          "reduction momentum %.3e, reduction mass %.3e)\n",
          this_fn_name, stop_reason, itn, res_reduction, conv_factor,
          res_reduction_vel, res_reduction_press);
    }
  }
  else {
    RHEA_GLOBAL_PRODUCTIONF (
        "Done %s (solver stopping reason %i, num iter %i)\n",
        this_fn_name, stop_reason, itn);
  }
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
  const char         *this_fn_name =
                        "rhea_stokes_problem_nonlinear_update_operator_fn";
  rhea_stokes_problem_t *stokes_problem_nl = data;

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->stokes_op != NULL);
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (solution));

  /* compute coefficient */
  rhea_stokes_problem_set_velocity_pressure (stokes_problem_nl, solution);
  rhea_stokes_problem_compute_coefficient (stokes_problem_nl, 0 /* !init */);

  /* set viscous stress coefficient */
  {
    ymir_stress_op_t   *stress_op = stokes_problem_nl->stokes_op->stress_op;
    rhea_viscosity_options_t *visc_options = stokes_problem_nl->visc_options;

    RHEA_ASSERT (ymir_stress_op_has_linearized (stress_op));
    RHEA_ASSERT (ymir_stress_op_get_coeff_type (stress_op) ==
                 YMIR_STRESS_OP_COEFF_ISOTROPIC_SCAL);
    RHEA_ASSERT (fabs (ymir_stress_op_get_coeff_shift (stress_op) -
                       rhea_viscosity_get_visc_shift (visc_options)) < SC_EPS);
    RHEA_ASSERT (stokes_problem_nl->coeff != NULL);

    ymir_stress_op_set_coeff_scal (stress_op, stokes_problem_nl->coeff);
  }

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
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
  const char         *this_fn_name =
                        "rhea_stokes_problem_nonlinear_update_hessian_fn";
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

  RHEA_GLOBAL_VERBOSEF (
      "Into %s (linearization %i, sol %i, step %i)\n",
      this_fn_name, (int) linearization_type, solution_exists, step_exists);

  /* check input */
#ifdef RHEA_ENABLE_DEBUG
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->stokes_op != NULL);
  RHEA_ASSERT (stokes_problem_nl->sol_vel != NULL);
  if (solution != NULL) {
    RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (solution));
  }
  if (step_vec != NULL) {
    RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (step_vec));
  }
#endif

  /* get data from nonlinear Stokes problem object */
  stokes_op = stokes_problem_nl->stokes_op;
  coeff = stokes_problem_nl->coeff;
  sol_vel = stokes_problem_nl->sol_vel;

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
        rhea_velocity_pressure_copy_components (
            sol_vel, NULL, solution, stokes_problem_nl->press_elem);
        rhea_stokes_problem_velocity_boundary_set_zero (
            sol_vel, stokes_problem_nl);
        RHEA_ASSERT (rhea_velocity_is_valid (sol_vel));

        /* compute and normalize strain rate tensor */
        ymir_stress_op_optimized_compute_strain_rate (
            proj_tens, sol_vel, stress_op);
        ymir_stress_op_tensor_normalize (proj_tens);
      }
      else { /* if solution is not provided */
        ymir_vec_set_zero (sol_vel);
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
        rhea_velocity_pressure_copy_components (
            sol_vel, NULL, solution, stokes_problem_nl->press_elem);
        rhea_stokes_problem_velocity_boundary_set_zero (
            sol_vel, stokes_problem_nl);
        RHEA_ASSERT (rhea_velocity_is_valid (sol_vel));

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
          stokes_problem_nl->viscosity_compute_fn (
              prev_coeff, NULL, NULL, NULL,
              temperature, weakzone, prev_vel,
              stokes_problem_nl->viscosity_compute_fn_data);

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
              this_fn_name, alignment);
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
        rhea_velocity_pressure_copy_components (
            sol_vel, NULL, solution, stokes_problem_nl->press_elem);
        rhea_stokes_problem_velocity_boundary_set_zero (
            sol_vel, stokes_problem_nl);
        RHEA_ASSERT (rhea_velocity_is_valid (sol_vel));

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
              this_fn_name, mean, min, max, max/min);
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
        rhea_velocity_pressure_copy_components (
            sol_vel, NULL, solution, stokes_problem_nl->press_elem);
        rhea_stokes_problem_velocity_boundary_set_zero (
            sol_vel, stokes_problem_nl);
        RHEA_ASSERT (rhea_velocity_is_valid (sol_vel));

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
              this_fn_name, alignment_avg);
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
    stokes_problem_nl->stokes_pc = ymir_nlstokes_pc_new (stokes_op);
  }
  else {
    /* update Stokes preconditioner */
    ymir_nlstokes_pc_recompute (stokes_problem_nl->stokes_pc);
  }
#else
  RHEA_ABORT_NOT_REACHED ();
#endif

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
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
  const char         *this_fn_name =
                        "rhea_stokes_problem_nonlinear_create_solver_data_fn";
  rhea_stokes_problem_t *stokes_problem_nl = data;
  const char         *vtk_path = stokes_problem_nl->solver_vtk_path;
  const int           nonzero_init_guess = (solution != NULL);

  RHEA_GLOBAL_VERBOSEF ("Into %s (nonzero init guess %i)\n",
                        this_fn_name, nonzero_init_guess);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->ymir_mesh  != NULL);
  RHEA_ASSERT (stokes_problem_nl->press_elem != NULL);
  RHEA_ASSERT (stokes_problem_nl->coeff   != NULL);
  RHEA_ASSERT (stokes_problem_nl->rhs_vel != NULL);
  RHEA_ASSERT (stokes_problem_nl->stokes_op == NULL);
  RHEA_ASSERT (stokes_problem_nl->stokes_pc == NULL);
  RHEA_ASSERT (stokes_problem_nl->norm_op   == NULL);

  /* set VTK path for debuggin */
  if (vtk_path != NULL) {
    char                path[BUFSIZ];

    snprintf (path, BUFSIZ, "%s_itn%02i", vtk_path, 0);
    ymir_vtk_set_debug_path (path);
  }

  /* compute viscosity */
  if (nonzero_init_guess) { /* if a velocity is given as an initial guess */
    RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (solution));
    RHEA_ASSERT (rhea_velocity_pressure_is_valid (solution));
    rhea_stokes_problem_set_velocity_pressure (stokes_problem_nl, solution);
    rhea_stokes_problem_compute_coefficient (stokes_problem_nl, 0 /* !init */);
  }
  else { /* otherwise assume zero velocity */
    rhea_stokes_problem_set_velocity_pressure (stokes_problem_nl, NULL);
    rhea_stokes_problem_compute_coefficient (stokes_problem_nl, 1 /* init */);
  }

  /* create Stokes operator */
  stokes_problem_nl->stokes_op = ymir_stokes_op_new_ext (
      stokes_problem_nl->coeff, stokes_problem_nl->vel_dir,
      NULL /* Robin BC's */,
      NULL /* deprecated */,
      NULL /* deprecated */,
      stokes_problem_nl->press_elem,
      stokes_problem_nl->domain_options->center,
      stokes_problem_nl->domain_options->moment_of_inertia);

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
  if (stokes_problem_nl->norm_type == RHEA_STOKES_NORM_HMINUS1_L2) {
    stokes_problem_nl->norm_op = ymir_Hminus1_norm_op_new (
        stokes_problem_nl->ymir_mesh,
        stokes_problem_nl->norm_op_mass_scaling);
  }

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
}

/**
 * Clears the data of the nonlinear problem after Newton solve.
 * (Callback function for Newton's method)
 */
static void
rhea_stokes_problem_nonlinear_clear_solver_data_fn (void *data)
{
  const char         *this_fn_name =
                        "rhea_stokes_problem_nonlinear_clear_solver_data_fn";
  rhea_stokes_problem_t *stokes_problem_nl = data;

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", this_fn_name);

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

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
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
  const char         *this_fn_name =
    "rhea_stokes_problem_nonlinear_compute_negative_gradient_fn";
  rhea_stokes_problem_t *stokes_problem_nl = data;
  ymir_stokes_op_t   *stokes_op = stokes_problem_nl->stokes_op;
  ymir_vec_t         *rhs_vel = stokes_problem_nl->rhs_vel;
  ymir_vec_t         *rhs_vel_nonzero_dirichlet =
                        stokes_problem_nl->rhs_vel_nonzero_dirichlet;
  ymir_vec_t         *rhs_vel_press = stokes_problem_nl->rhs_vel_press;

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", this_fn_name);

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
    ymir_stokes_pc_apply_stokes_op (solution, neg_gradient, stokes_op, 0, 1);
    ymir_vec_add (-1.0, rhs_vel_press, neg_gradient);
    ymir_vec_scale (-1.0, neg_gradient);
  }
  else {
    /* set residual to be the right-hand side */
    ymir_vec_copy (rhs_vel_press, neg_gradient);
  }
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (neg_gradient));

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
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
  const char         *this_fn_name =
    "rhea_stokes_problem_nonlinear_compute_gradient_norm_fn";
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
  const char         *this_fn_name =
    "rhea_stokes_problem_nonlinear_modify_hessian_system_fn";
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

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", this_fn_name);

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
            this_fn_name, neg_grad_aniso_part, neg_grad_aniso_part/normL2);
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
          this_fn_name, aniso_corr);

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

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
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
  const char         *this_fn_name =
                        "rhea_stokes_problem_nonlinear_solve_hessian_system_fn";
  const int           check_lin_aniso =
    rhea_stokes_problem_nonlinear_lin_aniso_check;
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
          this_fn_name, neg_grad_aniso_part, neg_grad_aniso_part/normL2);
    }

    if (isfinite (step_aniso_part)) {
      rhea_stokes_norm_compute (
          &normL2, NULL, step, RHEA_STOKES_NORM_L2_PRIMAL, NULL, press_elem);
      RHEA_GLOBAL_INFOF (
          "%s: Linearization-induced anisotropy part in step:     "
          "abs %+.2e, rel %+.2e\n",
          this_fn_name, step_aniso_part, step_aniso_part/normL2);
    }
  }

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);

  /* return iteraton count and "stopping" reason */
  if (lin_iter_count != NULL) {
    *lin_iter_count = itn;
  }
  return stop_reason;
}

/**
 * Writes or prints output at the beginning of a Newton step.
 * (Callback function for Newton's method)
 */
static void
rhea_stokes_problem_nonlinear_output_prestep_fn (ymir_vec_t *solution,
                                                 const int iter, void *data)
{
  const char         *this_fn_name =
                        "rhea_stokes_problem_nonlinear_output_prestep_fn";
  rhea_stokes_problem_t *stokes_problem_nl = data;
  const char         *vtk_path = stokes_problem_nl->solver_vtk_path;
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

  /* create visualization */
  if (vtk_path != NULL) {
    char                path[BUFSIZ];

    /* set path */
    snprintf (path, BUFSIZ, "%s_itn%02i", vtk_path, iter);

    /* write VTK */
    rhea_vtk_write_nonlinear_stokes_iteration (
        path, velocity, pressure, viscosity,
        stokes_problem_nl->bounds_marker, stokes_problem_nl->yielding_marker);

    /* set VTK path for debugging within ymir */
    ymir_vtk_set_debug_path (path);
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

static void
rhea_stokes_problem_nonlinear_create_mesh_data (
                                    rhea_stokes_problem_t *stokes_problem_nl,
                                    ymir_mesh_t *ymir_mesh,
                                    ymir_pressure_elem_t *press_elem)
{
  const char         *this_fn_name =
                        "rhea_stokes_problem_nonlinear_create_mesh_data";

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->ymir_mesh  == NULL);
  RHEA_ASSERT (stokes_problem_nl->press_elem == NULL);
  RHEA_ASSERT (stokes_problem_nl->coeff           == NULL);
  RHEA_ASSERT (stokes_problem_nl->proj_scal       == NULL);
  RHEA_ASSERT (stokes_problem_nl->bounds_marker   == NULL);
  RHEA_ASSERT (stokes_problem_nl->yielding_marker == NULL);
  RHEA_ASSERT (stokes_problem_nl->weakzone  == NULL);
  RHEA_ASSERT (stokes_problem_nl->sol_vel   == NULL);
  RHEA_ASSERT (stokes_problem_nl->proj_tens == NULL);
  RHEA_ASSERT (stokes_problem_nl->rhs_vel       == NULL);
  RHEA_ASSERT (stokes_problem_nl->rhs_vel_press == NULL);
  RHEA_ASSERT (stokes_problem_nl->newton_neg_gradient_vec == NULL);
  RHEA_ASSERT (stokes_problem_nl->newton_step_vec         == NULL);
  RHEA_ASSERT (stokes_problem_nl->domain_options != NULL);

  /* assign mesh objects */
  stokes_problem_nl->ymir_mesh = ymir_mesh;
  stokes_problem_nl->press_elem = press_elem;

  /* create velocity boundary conditions */
  rhea_stokes_problem_velocity_boundary_create (stokes_problem_nl);

  /* create coefficient and vectors related to it */
  if (stokes_problem_nl->weakzone_compute_fn != NULL) {
    stokes_problem_nl->weakzone = rhea_weakzone_new (ymir_mesh);
    stokes_problem_nl->weakzone_compute_fn (
        stokes_problem_nl->weakzone,
        stokes_problem_nl->weakzone_compute_fn_data);
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

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
}

static void
rhea_stokes_problem_nonlinear_clear_mesh_data (
                                    rhea_stokes_problem_t *stokes_problem_nl)
{
  const char         *this_fn_name =
                        "rhea_stokes_problem_nonlinear_clear_mesh_data";

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (stokes_problem_nl->type == RHEA_STOKES_PROBLEM_NONLINEAR);
  RHEA_ASSERT (stokes_problem_nl->coeff           != NULL);
  RHEA_ASSERT (stokes_problem_nl->bounds_marker   != NULL);
  RHEA_ASSERT (stokes_problem_nl->yielding_marker != NULL);
  RHEA_ASSERT (stokes_problem_nl->sol_vel         != NULL);
  RHEA_ASSERT (stokes_problem_nl->rhs_vel         != NULL);
  RHEA_ASSERT (stokes_problem_nl->rhs_vel_press   != NULL);
  RHEA_ASSERT (stokes_problem_nl->newton_problem          != NULL);
  RHEA_ASSERT (stokes_problem_nl->newton_neg_gradient_vec != NULL);
  RHEA_ASSERT (stokes_problem_nl->newton_step_vec         != NULL);

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

  /* destroy velocity boundary conditions */
  rhea_stokes_problem_velocity_boundary_clear (stokes_problem_nl);

  /* remove vector references */
  rhea_stokes_problem_remove_temperature (stokes_problem_nl);
  rhea_stokes_problem_remove_velocity_pressure (stokes_problem_nl);

  /* remove mesh references */
  stokes_problem_nl->ymir_mesh = NULL;
  stokes_problem_nl->press_elem = NULL;

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
}

static rhea_stokes_problem_t *
rhea_stokes_problem_nonlinear_new (ymir_mesh_t *ymir_mesh,
                                   ymir_pressure_elem_t *press_elem,
                                   ymir_vec_t *temperature,
                                   rhea_domain_options_t *domain_options,
                                   rhea_temperature_options_t *temp_options,
                                   rhea_weakzone_options_t *weak_options,
                                   rhea_viscosity_options_t *visc_options,
                                   void *solver_options)
{
  const char         *this_fn_name = "rhea_stokes_problem_nonlinear_new";
  rhea_stokes_problem_t *stokes_problem_nl;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s (linearization %i)\n", this_fn_name,
                           rhea_stokes_problem_nonlinear_linearization_type);

  /* check input */
  RHEA_ASSERT (visc_options->type == RHEA_VISCOSITY_NONLINEAR);

  /* initialize general Stokes problem structure */
  stokes_problem_nl = rhea_stokes_problem_struct_new (
      RHEA_STOKES_PROBLEM_NONLINEAR,
      domain_options, temp_options, weak_options, visc_options);

  /* set temperature */
  rhea_stokes_problem_set_temperature (stokes_problem_nl, temperature);

  /* initialize nonlinear structure */
  stokes_problem_nl->linearization_type =
    (rhea_stokes_problem_nonlinear_linearization_t)
    rhea_stokes_problem_nonlinear_linearization_type;
  stokes_problem_nl->newton_options = solver_options;
  stokes_problem_nl->norm_type =
    (rhea_stokes_norm_type_t) rhea_stokes_problem_nonlinear_norm_type;
  stokes_problem_nl->norm_op_mass_scaling =
    rhea_stokes_problem_nonlinear_norm_mass_scaling;

  /* create mesh data */
  rhea_stokes_problem_nonlinear_create_mesh_data (stokes_problem_nl, ymir_mesh,
                                                  press_elem);

  /* create Newton problem */
  {
    const int           grad_norm_n_components = 2;
    rhea_newton_problem_t *newton_problem;

    newton_problem = rhea_newton_problem_new (
        rhea_stokes_problem_nonlinear_compute_negative_gradient_fn,
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
        rhea_stokes_problem_nonlinear_clear_solver_data_fn, newton_problem);
    rhea_newton_problem_set_conv_criterion_fn (
        RHEA_NEWTON_CONV_CRITERION_RESIDUAL_NORM,
        NULL /* objective functional is not provided */,
        rhea_stokes_problem_nonlinear_compute_gradient_norm_fn,
        grad_norm_n_components, newton_problem);
    rhea_newton_problem_set_apply_hessian_fn (
        rhea_stokes_problem_nonlinear_apply_hessian_fn, newton_problem);
    rhea_newton_problem_set_update_fn (
        rhea_stokes_problem_nonlinear_update_operator_fn,
        rhea_stokes_problem_nonlinear_update_hessian_fn,
        rhea_stokes_problem_nonlinear_modify_hessian_system_fn, newton_problem);
    rhea_newton_problem_set_output_fn (
        rhea_stokes_problem_nonlinear_output_prestep_fn, newton_problem);

    if (rhea_stokes_problem_nonlinear_check_jacobian) {
      rhea_newton_problem_set_checks (
          0 /* grad */, 1 /* Hessian */, newton_problem);
    }
  }

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

  /* destroy Newton problem */
  rhea_newton_problem_destroy (stokes_problem_nl->newton_problem);

  /* destroy mesh data */
  rhea_stokes_problem_nonlinear_clear_mesh_data (stokes_problem_nl);

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

/******************************************************************************
 * General Stokes Problem
 *****************************************************************************/

rhea_stokes_problem_t *
rhea_stokes_problem_new (ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         ymir_vec_t *temperature,
                         rhea_domain_options_t *domain_options,
                         rhea_temperature_options_t *temp_options,
                         rhea_weakzone_options_t *weak_options,
                         rhea_viscosity_options_t *visc_options,
                         void *solver_options)
{
  switch (visc_options->type) {
  case RHEA_VISCOSITY_LINEAR:
    return rhea_stokes_problem_linear_new (
        ymir_mesh, press_elem, temperature, domain_options, temp_options,
        weak_options, visc_options, solver_options);
  case RHEA_VISCOSITY_NONLINEAR:
    return rhea_stokes_problem_nonlinear_new (
        ymir_mesh, press_elem, temperature, domain_options, temp_options,
        weak_options, visc_options, solver_options);
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
  const char         *this_fn_name =
                        "rhea_stokes_problem_create_mesh_dependencies";

  RHEA_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  switch (stokes_problem->type) {
  case RHEA_STOKES_PROBLEM_LINEAR:
    rhea_stokes_problem_linear_create_mesh_data (
        stokes_problem, ymir_mesh, press_elem);
    rhea_stokes_problem_linear_create_solver_data (stokes_problem);
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
    rhea_stokes_problem_nonlinear_create_solver_data_fn (
        stokes_problem->velocity_pressure, stokes_problem);
    break;
  default: /* unknown Stokes type */
    RHEA_ABORT_NOT_REACHED ();
  }

  RHEA_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

void
rhea_stokes_problem_clear_mesh_dependencies (
                                        rhea_stokes_problem_t *stokes_problem)
{
  const char         *this_fn_name =
                        "rhea_stokes_problem_clear_mesh_dependencies";

  RHEA_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  switch (stokes_problem->type) {
  case RHEA_STOKES_PROBLEM_LINEAR:
    rhea_stokes_problem_linear_clear_solver_data (stokes_problem);
    rhea_stokes_problem_linear_clear_mesh_data (stokes_problem);
    break;
  case RHEA_STOKES_PROBLEM_NONLINEAR:
    rhea_stokes_problem_nonlinear_clear_solver_data_fn (stokes_problem);
    rhea_stokes_problem_nonlinear_clear_mesh_data (stokes_problem);
    break;
  default: /* unknown Stokes type */
    RHEA_ABORT_NOT_REACHED ();
  }

  RHEA_GLOBAL_INFOF ("Done %s\n", this_fn_name);
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

void
rhea_stokes_problem_set_solver_vtk_output (
                                    rhea_stokes_problem_t *stokes_problem,
                                    char *vtk_path)
{
  stokes_problem->solver_vtk_path = vtk_path;
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

rhea_domain_options_t *
rhea_stokes_problem_get_domain_options (rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->domain_options;
}

rhea_temperature_options_t *
rhea_stokes_problem_get_temperature_options (
                                        rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->temp_options;
}

rhea_weakzone_options_t *
rhea_stokes_problem_get_weakzone_options (
                                        rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->weak_options;
}

rhea_viscosity_options_t *
rhea_stokes_problem_get_viscosity_options (
                                        rhea_stokes_problem_t *stokes_problem)
{
  return stokes_problem->visc_options;
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
rhea_stokes_prbolem_set_weakzone_compute_fn (
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_weakzone_compute_fn_t fn,
                                    void *data)
{
  /* check input */
  RHEA_ASSERT (stokes_problem->ymir_mesh != NULL);

  /* set function and corresponding data */
  stokes_problem->weakzone_compute_fn = fn;
  stokes_problem->weakzone_compute_fn_data = data;

  /* recompute */
  if (stokes_problem->weakzone_compute_fn != NULL) { /* if fn. provided */
    if (stokes_problem->weakzone == NULL) {
      stokes_problem->weakzone = rhea_weakzone_new (stokes_problem->ymir_mesh);
    }
    stokes_problem->weakzone_compute_fn (
        stokes_problem->weakzone,
        stokes_problem->weakzone_compute_fn_data);
  }
  else { /* if no function provided */
    if (stokes_problem->weakzone != NULL) {
      rhea_weakzone_destroy (stokes_problem->weakzone);
      stokes_problem->weakzone = NULL;
    }
  }
}

void
rhea_stokes_prbolem_set_viscosity_compute_fn (
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
rhea_stokes_prbolem_set_rhs_vel_compute_fn (
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
rhea_stokes_prbolem_set_rhs_vel_nonzero_dir_compute_fn (
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
