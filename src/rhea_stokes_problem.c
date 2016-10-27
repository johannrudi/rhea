/*
 */

#include <rhea_stokes_problem.h>
#include <rhea_base.h>
#include <rhea_velocity.h>
#include <ymir_stokes_op.h>
#include <ymir_stokes_pc.h>
#include <ymir_stokes_vec.h>

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
  rhs_vel_press = ymir_stokes_vec_new (ymir_mesh, press_elem);
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
  ymir_vec_destroy (stokes_problem_lin->coeff);
  ymir_vec_destroy (stokes_problem_lin->rhs_vel);
  ymir_vec_destroy (stokes_problem_lin->rhs_vel_press);

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

/******************************************************************************
 * Nonlinear Stokes Problem
 *****************************************************************************/

