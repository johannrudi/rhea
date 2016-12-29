/*
  This file is part of the ymir Library.
  ymir is a C library for modeling ice sheets

  Copyright (C) 2010, 2011 Carsten Burstedde, Toby Isaac, Georg Stadler,
                           Lucas Wilcox.

  The ymir Library is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The ymir Library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the ymir Library.  If not, see <http://www.gnu.org/licenses/>.

  ---

  This is the slabs example for global instantaneous mantle flow with plates.

*/

#include <slabs_linear_stokes_problem.h>
#include <ymir_comm.h>
#include <ymir_mass_vec.h>
#include <ymir_velocity_vec.h>
#include <ymir_pressure_vec.h>
#include <ymir_stokes_vec.h>
#include <ymir_interp_vec.h>

/**
 * Sets Dirichlet boundary conditions.
 */
ymir_vel_dir_t *
slabs_set_dirichlet_bc (ymir_mesh_t *mesh, ymir_vec_t *dirscal, void *data)
{
  slabs_physics_options_t *physics_options = (slabs_physics_options_t *) data;
  ymir_vel_dir_t     *vel_dir;

  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_BRICK:
  case SL_DOMAIN_SHELL_CHUNK:
  case SL_DOMAIN_SHELL_SLICE:
    vel_dir = ymir_vel_dir_new_from_faces (mesh, NULL, slabs_cube_vel_dir_fn,
                                           &(physics_options->bc_type),
                                           dirscal);
    break;

  case SL_DOMAIN_SHELL:
    vel_dir = ymir_vel_dir_new_from_faces (mesh, NULL, slabs_shell_vel_dir_fn,
                                           &(physics_options->bc_type),
                                           dirscal);
    if (physics_options->bc_type == SL_VEL_BC_DIRICHLET_NORM_FIXDOF) {
      slabs_shell_vel_dir_restrict_dof (vel_dir);
    }
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* set zero Dirichlet boundary */
  vel_dir->nonzero = 0;

  return vel_dir;
}

/**
 * Creates a new Stokes operator.
 */
ymir_stokes_op_t *
slabs_stokes_op_new (ymir_vec_t *viscosity,
                     ymir_vel_dir_t *vel_dir,
                     ymir_vel_rob_t *vel_rob,
                     ymir_vec_t *usol,
                     ymir_vec_t *dvdIIe,
                     ymir_pressure_elem_t *press_elem,
                     slabs_physics_options_t *physics_options)
{
  /* check input */
#ifdef YMIR_DEBUG
  if (viscosity != NULL) {
    YMIR_ASSERT_IS_DVEC (viscosity);
  }
  if (usol != NULL) {
    YMIR_ASSERT_IS_CVEC (usol);
  }
  if (dvdIIe != NULL) {
    YMIR_ASSERT_IS_DVEC (dvdIIe);
  }
#endif

  /* create Stokes operator */
  return ymir_stokes_op_new_ext (
      viscosity, vel_dir, vel_rob, usol, dvdIIe, press_elem,
      physics_options->domain_center,
      physics_options->domain_moment_of_inertia);
}

/**
 * Destroys a Stokes operator.
 */
void
slabs_stokes_op_destroy (ymir_stokes_op_t *stokes_op)
{
  ymir_stokes_op_destroy (stokes_op);
}

/**
 * Creates a new linear Stokes problem.
 */
slabs_lin_stokes_problem_t *
slabs_linear_stokes_problem_new (slabs_stokes_state_t *state,
                                 ymir_mesh_t *mesh,
                                 ymir_pressure_elem_t *press_elem,
                                 slabs_physics_options_t *physics_options)
{
  const char         *this_fn_name = "slabs_linear_stokes_problem_new";
  ymir_vec_t         *viscosity;
  ymir_vel_dir_t     *vel_dir;
  ymir_vec_t         *dirscal;
  ymir_stokes_op_t   *stokes_op;
  ymir_stokes_pc_t   *stokes_pc;
  ymir_vec_t         *rhs_u_point;
  ymir_vec_t         *rhs;
  slabs_lin_stokes_problem_t  *lin_stokes;

  YMIR_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /*
   * Create Stokes Operator
   */

  /* set physical viscosity */
  viscosity = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
  slabs_physics_compute_lin_stokes_coeff (viscosity, state, physics_options);

  /* setup Dirichlet boundary conditions */
  if (physics_options->bc_default_dirichlet_scale) {
    dirscal = ymir_cvec_new (mesh, 1);
    ymir_velocity_bc_default_dirscal (dirscal, viscosity);
    slabs_cvec_set_min_value (dirscal, SL_DIRSCAL_MIN);
    YMIR_ASSERT (sc_dmatrix_is_valid (dirscal->cvec));
  }
  else {
    dirscal = NULL;
  }
  vel_dir = slabs_set_dirichlet_bc (mesh, dirscal, physics_options);

  /* create stokes operator */
  stokes_op = slabs_stokes_op_new (viscosity, vel_dir, NULL, NULL, NULL,
                                   press_elem, physics_options);

  /*
   * Create Stokes Preconditioner
   */

  /* build Stokes preconditioner */
  stokes_pc = ymir_stokes_pc_new (stokes_op);

  /*
   * Create Right-Hand Side
   */

  /* create right-hand side forcing */
  rhs_u_point = ymir_cvec_new (mesh, 3);
  slabs_physics_compute_rhs_u_point (rhs_u_point, state, physics_options);
  if (physics_options->rhs_random) { /* if set random rhs (for testing) */
#ifdef YMIR_PETSC
    ymir_petsc_vec_set_random (rhs_u_point, YMIR_MESH_PETSCLAYOUT_NONE);
#else
    YMIR_ABORT_NOT_REACHED ();
#endif
    if (physics_options->domain_shape == SL_DOMAIN_SHELL) { /* if shell */
      ymir_velocity_vec_project_out_mean_rotation (
          rhs_u_point, physics_options->domain_center,
          physics_options->domain_moment_of_inertia, 0);
    }
  }

  /* construct right-hand side for incompressible Stokes system */
  rhs = ymir_stokes_vec_new (mesh, press_elem);
  ymir_stokes_op_construct_rhs_ext (rhs_u_point, NULL, NULL, rhs,
                                    1 /* incompressible */, stokes_op);

  /*
   * Create Linear Stokes Problem
   */

  /* initialize linear Stokes problem */
  lin_stokes = YMIR_ALLOC (slabs_lin_stokes_problem_t, 1);

  /* assign pointers to linear Stokes problem struct */
  lin_stokes->mesh = mesh;
  lin_stokes->press_elem = press_elem;
  lin_stokes->vel_dir = vel_dir;
  lin_stokes->rhs_u_point = rhs_u_point;
  lin_stokes->rhs = rhs;
  lin_stokes->viscosity = viscosity;
  lin_stokes->stokes_op = stokes_op;
  lin_stokes->stokes_pc = stokes_pc;

  YMIR_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

  return lin_stokes;
}

/**
 * Destroys a linear Stokes problem.
 */
void
slabs_linear_stokes_problem_destroy (slabs_lin_stokes_problem_t *lin_stokes)
{
  const char         *this_fn_name = "slabs_linear_stokes_problem_destroy";

  YMIR_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* destroy Stokes operator and preconditioner */
  ymir_stokes_pc_destroy (lin_stokes->stokes_pc);
  slabs_stokes_op_destroy (lin_stokes->stokes_op);

  /* destroy vectors */
  ymir_vec_destroy (lin_stokes->viscosity);
  ymir_vec_destroy (lin_stokes->rhs_u_point);
  ymir_vec_destroy (lin_stokes->rhs);

  /* destroy Dirichlet BC's */
  if (lin_stokes->vel_dir != NULL) {
    if (lin_stokes->vel_dir->scale != NULL) {
      ymir_vec_destroy (lin_stokes->vel_dir->scale);
    }
    ymir_vel_dir_destroy (lin_stokes->vel_dir);
  }

  /* destroy linear Stokes problem */
  YMIR_FREE (lin_stokes);

  YMIR_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Projects out null spaces of velocity and pressure.
 */
void
slabs_linear_stokes_problem_project_out_nullspace (ymir_vec_t *up,
                                                   ymir_stokes_op_t *stokes_op,
                                                   slabs_physics_options_t
                                                     *physics_options,
                                                   const int residual_space)
{
  const char         *this_fn_name =
                        "slabs_linear_stokes_problem_project_out_nullspace";

  YMIR_GLOBAL_INFOF ("Into %s (residual space %i)\n", this_fn_name,
                     residual_space);

  /*
   * Project Out Mean Velocity Rotation
   */

  if (physics_options->domain_shape == SL_DOMAIN_SHELL) { /* if shell domain */
    ymir_vec_t         *u;
#ifdef YMIR_DEBUG
    ymir_pressure_elem_t *press_elem = stokes_op->press_elem;
    double              rot_axis[3], rot_axis_proj[3];

    /* compute mean rotation */
    u = ymir_cvec_new (up->mesh, 3);
    ymir_stokes_vec_get_velocity (up, u, press_elem);
    ymir_stress_op_compute_mean_rotation (rot_axis, u,
                                          stokes_op->stress_op,
                                          residual_space);
    ymir_vec_destroy (u);
#endif

    /* remove rotation */
    u = ymir_vec_new_blank ();
    ymir_upvec_get_u (up, u, YMIR_RW);
    ymir_stress_op_project_out_nullspace (u, stokes_op->stress_op,
                                          YMIR_STRESS_OP_TRANSL_PROJECT_NONE,
                                          YMIR_STRESS_OP_ROT_PROJECT_SYMM,
                                          residual_space);
    ymir_upvec_set_u (up, u, YMIR_SET);
    ymir_vec_destroy (u);

#ifdef YMIR_DEBUG
    /* compute and print mean rotation */
    u = ymir_cvec_new (up->mesh, 3);
    ymir_stokes_vec_get_velocity (up, u, press_elem);
    ymir_stress_op_compute_mean_rotation (rot_axis_proj, u,
                                          stokes_op->stress_op,
                                          residual_space);
    ymir_vec_destroy (u);

    YMIR_GLOBAL_INFOF (
        "%s: Mean rotatation %1.3e, %1.3e, %1.3e, "
        "after projection %1.3e, %1.3e, %1.3e\n",
        this_fn_name, rot_axis[0], rot_axis[1], rot_axis[2],
        rot_axis_proj[0], rot_axis_proj[1], rot_axis_proj[2]);
#endif
  }

  /*
   * Project Out Mean Pressure
   */

  {
    ymir_pressure_elem_t *press_elem = stokes_op->press_elem;
    ymir_vec_t         *p;
#ifdef YMIR_DEBUG
    double              mean_press, mean_press_proj;

    /* compute mean pressure */
    p = ymir_pressure_vec_new (up->mesh, press_elem);
    ymir_stokes_vec_get_pressure (up, p, press_elem);
    mean_press = ymir_pressure_vec_mean (p, press_elem, residual_space);
    ymir_vec_destroy (p);
#endif

    /* remove mean pressure */
    p = ymir_vec_new_blank ();
    ymir_upvec_get_p (up, p, YMIR_RW);
    ymir_pressure_vec_project_out_mean (p, press_elem, residual_space);
    ymir_upvec_set_p (up, p, YMIR_SET);
    ymir_vec_destroy (p);

#ifdef YMIR_DEBUG
    /* compute and print mean pressure */
    p = ymir_pressure_vec_new (up->mesh, press_elem);
    ymir_stokes_vec_get_pressure (up, p, press_elem);
    mean_press_proj = ymir_pressure_vec_mean (p, press_elem, residual_space);
    ymir_vec_destroy (p);

    YMIR_GLOBAL_INFOF ("%s: Mean pressure %1.3e, after projection %1.3e\n",
                       this_fn_name, mean_press, mean_press_proj);
#endif
  }

  YMIR_GLOBAL_INFOF ("Done %s (residual space %i)\n", this_fn_name,
                     residual_space);
}

