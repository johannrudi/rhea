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

#include <slabs_vtk.h>
#include <p8est_vtk.h>
#include <mangll_vtk.h>
#include <ymir_vtk.h>
#include <ymir_mass_vec.h>
#include <ymir_velocity_vec.h>
#include <ymir_pressure_vec.h>
#include <ymir_stokes_vec.h>
#include <ymir_interp_vec.h>
#include <slabs_norm.h>
#include <slabs_physics_extended.h>
#include <slabs_discretization_extended.h>

void
slabs_vtk_write_p8est (const char *filepath, p8est_t *p4est,
                       const slabs_domain_shape_t domain_shape)
{
  const char         *this_fn_name = "slabs_vtk_write_p8est";
  p8est_connectivity_t *conn;
  p8est_geometry_t   *geom;

  YMIR_GLOBAL_INFOF ("Into %s: Write VTK file \"%s\"\n",
                     this_fn_name, filepath);

  switch (domain_shape) {
  case SL_DOMAIN_SHELL:
  case SL_DOMAIN_SHELL_CHUNK:
    conn = p8est_connectivity_new_shell ();
    geom = p8est_geometry_new_shell (conn, SL_SHELL_RADIUS_TOP,
                                     SL_SHELL_RADIUS_BOTTOM);
    break;

  default: /* all other domain types */
    conn = p8est_connectivity_new_unitcube ();
    geom = p8est_geometry_new_connectivity (conn);
  }

  /* write vtk file */
  p8est_vtk_write_file (p4est, geom, filepath);

  /* destroy */
  p8est_geometry_destroy (geom);
  p8est_connectivity_destroy (conn);

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

int
slabs_vtk_write_mangll (const char *filepath, mangll_t *ma, double *field,
                        char *name)
{
  const char         *this_fn_name = "slabs_vtk_write_mangll";
  int                 error;

  YMIR_GLOBAL_INFOF ("Into %s: Write VTK file \"%s\"\n",
                     this_fn_name, filepath);

  /* write header */
  error = mangll_vtk_write_header (ma, filepath, name, NULL, 1);
  if (error) {
    YMIR_GLOBAL_LERRORF ("%s: Error writing mangll header in vtk file \"%s\". "
                         "Error no: %i", this_fn_name, filepath, error);
    return error;
  }

  /* write scalar field */
  error = mangll_vtk_write_point_scalar (ma, filepath, field, 1, name);
  if (error) {
    YMIR_GLOBAL_LERRORF ("%s: Error writing mangll scalar field %s "
                         "in vtk file \"%s\". Error no: %i",
                         this_fn_name, name, filepath, error);
    return error;
  }

  /* write footer */
  error = mangll_vtk_write_footer (ma, filepath);
  if (error) {
    YMIR_GLOBAL_LERRORF ("%s: Error writing mangll footer in vtk file \"%s\". "
                         "Error no: %i", this_fn_name, filepath, error);
    return error;
  }

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);

  return 0;
}

/**
 *
 */
static ymir_vec_t *
slabs_vtk_create_vec_from_indicator (
                                   const slabs_discr_amr_indicator_type_t type,
                                   slabs_stokes_state_t *state,
                                   ymir_mesh_t *mesh,
                                   slabs_physics_options_t *physics_options,
                                   slabs_discr_options_t *discr_options,
                                   const int init_amr,
                                   const int import_data)
{
  mangll_t           *mangll = mesh->ma;
  mangll_cnodes_t    *cnodes = mesh->cnodes;
  const mangll_locidx_t  n_quadrants_loc = cnodes->K;
  const double        tol_min = 0.0;
  const double        tol_max = 1.0;
  const int8_t        level_min = 1;
  const int8_t        level_max = 2;
  slabs_discr_amr_indicator_t *ind;
  sc_dmatrix_t       *ind_val;
  ymir_vec_t         *ind_vec;

  /* create indicator */
  if (init_amr && !import_data) {
    ind_val = slabs_discr_init_amr_indicator_new (mangll, physics_options,
                                                  discr_options, type);
  }
  else {
    ind = slabs_discr_amr_indicator_new (type, tol_min, tol_max,
                                         level_min, level_max,
                                         n_quadrants_loc);
    slabs_discr_amr_indicator_set (ind, state, mangll, cnodes,
                                   physics_options, discr_options, init_amr);
    ind_val = ind->val;
  }

  /* create vector from indicator */
  ind_vec = slabs_dvec_new_from_element_data (ind_val, mesh, YMIR_GAUSS_NODE);

  /* destroy */
  if (init_amr && !import_data) {
    sc_dmatrix_destroy (ind_val);
  }
  else {
    slabs_discr_amr_indicator_destroy (ind);
  }

  /* return indicator vector */
  return ind_vec;
}

/**
 * Writes VTK files of input data.
 */
static void
slabs_vtk_write_input (const char *filepath,
                       slabs_stokes_state_t *state,
                       ymir_vec_t *rhs_u_point,
                       ymir_vel_dir_t *vel_dir,
                       slabs_physics_options_t *physics_options,
                       slabs_discr_options_t *discr_options)
{
  const char         *this_fn_name = "slabs_vtk_write_input";
  ymir_mesh_t        *mesh = rhs_u_point->mesh;
  char                path[BUFSIZ];

  YMIR_GLOBAL_INFOF ("Into %s: Write VTK file \"%s\"\n",
                     this_fn_name, filepath);

  /* check input */
  YMIR_ASSERT (state->temp_vec != NULL);
  YMIR_ASSERT (state->weak_vec != NULL);

  /*
   * Write main fields
   */

  {
    const slabs_viscosity_t  viscosity_type = physics_options->viscosity_type;
    ymir_vec_t         *background_temp = ymir_cvec_new (mesh, 1);
    ymir_vec_t         *visc_temp = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);

    /* compute background temperature */
    ymir_cvec_set_function (background_temp,
                            slabs_physics_background_temp_set_fn,
                            physics_options);

    /* compute temperature (only) dependent viscosity */
    if (viscosity_type == SL_VISCOSITY_NONLINEAR) {
      physics_options->viscosity_type = SL_VISCOSITY_LINEAR;
    }
    slabs_physics_compute_linear_viscosity (visc_temp, state, physics_options);
    physics_options->viscosity_type = viscosity_type;

    /* set vtk path */
    snprintf (path, BUFSIZ, "%s_main", filepath);

    /* write vtk file */
    ymir_vtk_write (
        mesh, path,
        state->temp_vec, "temperature", background_temp, "background_temp",
        state->weak_vec, "weakzone", visc_temp, "visc_temp",
        rhs_u_point, "rhs_u", NULL);

    /* destroy */
    ymir_vec_destroy (background_temp);
    ymir_vec_destroy (visc_temp);
  }

  /*
   * Write boundary conditions
   */

  if (vel_dir != NULL) {
    /* set vtk path */
    snprintf (path, BUFSIZ, "%s_bc", filepath);

    if (physics_options->domain_shape == SL_DOMAIN_SHELL) { /* if shell shape */
      ymir_vec_t         *vel_dir_nvec_surf, *vel_dir_nvec_inner;
      ymir_vec_t         *vel_dir_scale_surf, *vel_dir_scale_inner;

      /* create vectors */
      vel_dir_nvec_surf = ymir_face_cvec_new (mesh, SL_TOP, 3);
      vel_dir_nvec_inner = ymir_face_cvec_new (mesh, SL_BASE, 3);
      vel_dir_scale_surf = ymir_face_cvec_new (mesh, SL_TOP, 1);
      vel_dir_scale_inner = ymir_face_cvec_new (mesh, SL_BASE, 1);

      /* interpolate normal vectors onto faces */
      ymir_interp_vec (vel_dir->nvec, vel_dir_nvec_surf);
      ymir_interp_vec (vel_dir->nvec, vel_dir_nvec_inner);

      /* set scaling of Dirichlet boundary nodes */
      if (vel_dir->scale != NULL) {
        ymir_interp_vec (vel_dir->scale, vel_dir_scale_surf);
        ymir_interp_vec (vel_dir->scale, vel_dir_scale_inner);
      }
      else {
        ymir_vec_set_value (vel_dir_scale_surf, -1.0);
        ymir_vec_set_value (vel_dir_scale_inner, -1.0);
      }

      /* write vtk file */
      ymir_vtk_write (
          mesh, path,
          vel_dir_nvec_surf, "dir_nvec_surf",
          vel_dir_nvec_inner, "dir_nvec_inner",
          vel_dir_scale_surf, "dir_scale_surf",
          vel_dir_scale_inner, "dir_scale_inner", NULL);

      /* destroy */
      ymir_vec_destroy (vel_dir_nvec_surf);
      ymir_vec_destroy (vel_dir_nvec_inner);
      ymir_vec_destroy (vel_dir_scale_surf);
      ymir_vec_destroy (vel_dir_scale_inner);
    }
    else { /* if domain has any other shape than shell */
      ymir_vec_t         *vel_dir_nvec = vel_dir->nvec;
      ymir_vec_t         *vel_dir_scale = ymir_cvec_new (mesh, 1);;

      /* set scaling of Dirichlet boundary nodes */
      if (vel_dir->scale != NULL) {
        ymir_cvec_copy (vel_dir->scale, vel_dir_scale);
      }
      else {
        ymir_cvec_set_value (vel_dir_scale, -1.0);
      }

      /* write vtk file */
      ymir_vtk_write (
          mesh, path,
          vel_dir_nvec, "dir_nvec", vel_dir_scale, "dir_scale", NULL);

      /* destroy */
      ymir_vec_destroy (vel_dir_scale);
    }
  }

  /*
   * Write init AMR indicators
   */

  {
    const int           import_data =
      ( physics_options->temperature_type == SL_TEMP_IMPORT_FILE ||
        physics_options->weakzone_type == SL_WEAKZONE_IMPORT_FILE );
    const slabs_discr_amr_indicator_type_t  visc_indicator_type =
      discr_options->init_amr_visc_indicator_type;
    const slabs_discr_amr_indicator_type_t  weak_subdu_indicator_type =
      discr_options->init_amr_weak_subdu_indicator_type;
    const slabs_discr_amr_indicator_type_t  weak_ridge_indicator_type =
      discr_options->init_amr_weak_ridge_indicator_type;
    const slabs_discr_amr_indicator_type_t  weak_import_indicator_type =
      discr_options->init_amr_weak_import_indicator_type;
    const slabs_discr_amr_indicator_type_t  rhs_indicator_type =
      discr_options->init_amr_rhs_indicator_type;
    ymir_vec_t         *visc_amr_indicator_vec;
    ymir_vec_t         *weak_amr_indicator_vec;
    ymir_vec_t         *rhs_amr_indicator_vec;

    /* compute viscosity AMR indicator */
    if (visc_indicator_type != SL_AMR_INDICATOR_NONE) {
      visc_amr_indicator_vec = slabs_vtk_create_vec_from_indicator (
          visc_indicator_type, state, mesh, physics_options, discr_options,
          1 /* init_amr */, import_data);
    }
    else {
      visc_amr_indicator_vec = ymir_dvec_new_zero (mesh, 1, YMIR_GAUSS_NODE);
    }

    /* compute weak zone AMR indicator */
    if (import_data && weak_import_indicator_type != SL_AMR_INDICATOR_NONE) {
      weak_amr_indicator_vec = slabs_vtk_create_vec_from_indicator (
          weak_import_indicator_type, state, mesh,
          physics_options, discr_options, 1 /* init_amr */, import_data);
    }
    else if (!import_data) {
      mangll_t           *mangll = mesh->ma;
      ymir_vec_t         *weak_subdu_amr_indicator_vec;
      ymir_vec_t         *weak_ridge_amr_indicator_vec;

      if (weak_subdu_indicator_type != SL_AMR_INDICATOR_NONE) {
        sc_dmatrix_t       *weak_subdu_amr_indicator_mat;

        weak_subdu_amr_indicator_mat = slabs_discr_init_amr_indicator_new (
            mangll, physics_options, discr_options, weak_subdu_indicator_type);
        weak_subdu_amr_indicator_vec = slabs_dvec_new_from_element_data (
            weak_subdu_amr_indicator_mat, mesh, YMIR_GAUSS_NODE);
        sc_dmatrix_destroy (weak_subdu_amr_indicator_mat);
      }
      else {
        weak_subdu_amr_indicator_vec = ymir_dvec_new_zero (mesh, 1,
                                                           YMIR_GAUSS_NODE);
      }

      if (weak_ridge_indicator_type != SL_AMR_INDICATOR_NONE) {
        sc_dmatrix_t       *weak_ridge_amr_indicator_mat;

        weak_ridge_amr_indicator_mat = slabs_discr_init_amr_indicator_new (
            mangll, physics_options, discr_options, weak_ridge_indicator_type);
        weak_ridge_amr_indicator_vec = slabs_dvec_new_from_element_data (
            weak_ridge_amr_indicator_mat, mesh, YMIR_GAUSS_NODE);
        sc_dmatrix_destroy (weak_ridge_amr_indicator_mat);
      }
      else {
        weak_ridge_amr_indicator_vec = ymir_dvec_new_zero (mesh, 1,
                                                           YMIR_GAUSS_NODE);
      }

      weak_amr_indicator_vec = ymir_dvec_new_zero (mesh, 1, YMIR_GAUSS_NODE);
      ymir_dvec_add (1.0, weak_subdu_amr_indicator_vec, weak_amr_indicator_vec);
      ymir_dvec_add (1.0, weak_ridge_amr_indicator_vec, weak_amr_indicator_vec);

      ymir_vec_destroy (weak_subdu_amr_indicator_vec);
      ymir_vec_destroy (weak_ridge_amr_indicator_vec);
    }
    else {
      weak_amr_indicator_vec = ymir_dvec_new_zero (mesh, 1, YMIR_GAUSS_NODE);
    }

    /* compute right-hand side AMR indicator */
    if (rhs_indicator_type != SL_AMR_INDICATOR_NONE) {
      rhs_amr_indicator_vec = slabs_vtk_create_vec_from_indicator (
          rhs_indicator_type, state, mesh, physics_options, discr_options,
          1 /* init_amr */, import_data);
    }
    else {
      rhs_amr_indicator_vec = ymir_dvec_new_zero (mesh, 1, YMIR_GAUSS_NODE);
    }

    /* set vtk path */
    snprintf (path, BUFSIZ, "%s_init_amr", filepath);

    /* write vtk file */
    ymir_vtk_write (
        mesh, path,
        visc_amr_indicator_vec, "indicator_visc",
        weak_amr_indicator_vec, "indicator_weak",
        rhs_amr_indicator_vec, "indicator_rhs", NULL);

    /* destroy */
    ymir_vec_destroy (visc_amr_indicator_vec);
    ymir_vec_destroy (weak_amr_indicator_vec);
    ymir_vec_destroy (rhs_amr_indicator_vec);
  }

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

void
slabs_vtk_write_input_lin_stokes (const char *filepath,
                                  slabs_stokes_state_t *state,
                                  slabs_lin_stokes_problem_t *lin_stokes,
                                  slabs_physics_options_t *physics_options,
                                  slabs_discr_options_t *discr_options)
{
  ymir_vec_t         *rhs_u_point = lin_stokes->rhs_u_point;
  ymir_vel_dir_t     *vel_dir = lin_stokes->vel_dir;

  slabs_vtk_write_input (filepath, state, rhs_u_point, vel_dir,
                         physics_options, discr_options);
}

void
slabs_vtk_write_input_nl_stokes (const char *filepath,
                                 slabs_stokes_state_t *state,
                                 slabs_nl_stokes_problem_t *nl_stokes,
                                 slabs_physics_options_t *physics_options,
                                 slabs_discr_options_t *discr_options)
{
  ymir_vec_t         *rhs_u_point = nl_stokes->rhs_u_point;
  ymir_vel_dir_t     *vel_dir = nl_stokes->vel_dir;

  slabs_vtk_write_input (filepath, state, rhs_u_point, vel_dir,
                         physics_options, discr_options);
}

/**
 * Write VTK files for Stokes state fields and viscosity.
 */
void
slabs_vtk_write_state (const char *filepath,
                       slabs_stokes_state_t *state,
                       ymir_vec_t *viscosity_nondim,
                       const int compute_viscosity,
                       const int compute_stress_norm,
                       ymir_mesh_t *mesh,
                       ymir_pressure_elem_t *press_elem,
                       slabs_physics_options_t *physics_options)
{
  const slabs_viscosity_t  visc_type = physics_options->viscosity_type;
  const int           out_vel_press = (state->vel_press_vec != NULL);
  const int           out_visc = (viscosity_nondim != NULL ||
                                  compute_viscosity || compute_stress_norm);

  ymir_vec_t         *temp_nondim = state->temp_vec;
  ymir_vec_t         *u_surf = NULL;
  ymir_vec_t         *strain_rate_dim = NULL;
  ymir_vec_t         *stress_norm_dim = NULL;
  ymir_vec_t         *viscosity_dim = NULL;
  ymir_vec_t         *coeff = NULL;

  char                path[BUFSIZ];

  /* check input */
  YMIR_ASSERT (state->temp_vec != NULL);

  /*
   * Write VTK Files for State Variables
   */

  /* set vtk path */
  snprintf (path, BUFSIZ, "%s_state", filepath);

  if (!out_vel_press) { /* in output of temperature only */
    /* write vtk file */
    ymir_vtk_write (mesh, path, temp_nondim, "temperature", NULL);
  }
  else { /* if output of temperature, velocity, pressure */
    ymir_vec_t         *u_dim = ymir_cvec_new (mesh, 3);
    ymir_vec_t         *p_dim = ymir_pressure_vec_new (mesh, press_elem);
    ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (
                                          mesh->ma->N, mesh->ma->ompsize);

    /* get velocity and pressure from solution */
    ymir_stokes_vec_get_components (state->vel_press_vec, u_dim, p_dim,
                                    press_elem);

    /* compute 2nd invariant of the strain rate (used later) */
    strain_rate_dim = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
    ymir_second_invariant_vec (u_dim, strain_rate_dim, vel_elem);
    ymir_dvec_sqrt (strain_rate_dim, strain_rate_dim);

    /* interpolate fields to surface (used later) */
    u_surf = ymir_face_cvec_new (mesh, SL_TOP, 3);
    ymir_interp_vec (u_dim, u_surf);

    /* transform variables to physical dimensions */
    slabs_physics_transform_to_dimensional_velocity (u_dim);
    slabs_physics_transform_to_dimensional_velocity (u_surf);
    slabs_physics_transform_to_dimensional_pressure (p_dim);
    slabs_physics_transform_to_dimensional_strain_rate (strain_rate_dim);

    /* write vtk file */
    ymir_vtk_write (mesh, path, temp_nondim, "temperature", u_dim, "velocity",
                    p_dim, "pressure", NULL);

    /* destroy */
    ymir_vec_destroy (u_dim);
    ymir_vec_destroy (p_dim);
    ymir_velocity_elem_destroy (vel_elem);
  }

  /*
   * Write VTK Files for Derived Variables
   */

  /* set vtk path */
  snprintf (path, BUFSIZ, "%s_derived", filepath);

  /* compute normal stress at the surface */
  if (out_vel_press && compute_stress_norm) {
    ymir_vec_t         *coeff_deriv, *rank1_tensor_scal;
    ymir_vel_dir_t     *vel_dir;
    ymir_stokes_op_t   *stokes_op;
    ymir_vec_t         *rhs_u_point;

    /* set boudnary conditions */
    vel_dir = slabs_set_dirichlet_bc (mesh, NULL, physics_options);

    /* compute Stokes coefficient */
    coeff             = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
    coeff_deriv       = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
    rank1_tensor_scal = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
    slabs_physics_compute_stokes_coeff (
        coeff, coeff_deriv, rank1_tensor_scal, NULL, NULL,
        state, press_elem, vel_dir, physics_options);

    /* create Stokes operator */
    stokes_op = slabs_stokes_op_new (coeff, vel_dir, NULL, NULL, NULL,
                                     press_elem, physics_options);
    if (visc_type == SL_VISCOSITY_NONLINEAR) {
      ymir_stress_op_coeff_compute_rank1_tensor (stokes_op->stress_op,
                                                 rank1_tensor_scal,
                                                 state->vel_bc_vec);
    }

    /* create right-hand side forcing */
    rhs_u_point = ymir_cvec_new (mesh, 3);
    slabs_physics_compute_rhs_u_point (rhs_u_point, state, physics_options);

    /* compute surface stress */
    stress_norm_dim = ymir_face_cvec_new (mesh, SL_TOP, 1);
    slabs_physics_compute_normal_boundary_stress (
        stress_norm_dim, state->vel_press_vec, rhs_u_point, stokes_op);

    /* destroy */
    slabs_stokes_op_destroy (stokes_op);
    ymir_vec_destroy (rhs_u_point);
    ymir_vec_destroy (coeff_deriv);
    ymir_vec_destroy (rank1_tensor_scal);
    ymir_vel_dir_destroy (vel_dir);
  }

  /* compute viscosity if necessary */
  if (viscosity_nondim != NULL) {
    viscosity_dim = ymir_vec_clone (viscosity_nondim);
  }
  else if (viscosity_nondim == NULL && compute_viscosity) {
    viscosity_dim = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);

    if (coeff != NULL) {
      ymir_vec_copy (coeff, viscosity_dim);
      ymir_vec_scale (0.5, viscosity_dim);
    }
    else {
      if (visc_type != SL_VISCOSITY_NONLINEAR || state->vel_press_vec == NULL) {
        if (visc_type == SL_VISCOSITY_NONLINEAR) {
          physics_options->viscosity_type = SL_VISCOSITY_LINEAR;
        }

        /* compute linear (only temperature dependent) viscosity */
        slabs_physics_compute_linear_viscosity (viscosity_dim, state,
                                                physics_options);

        physics_options->viscosity_type = visc_type;
      }
      else {
        ymir_vec_t         *dvisc_dIIe;
        ymir_vel_dir_t     *vel_dir;

        dvisc_dIIe = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
        vel_dir = slabs_set_dirichlet_bc (mesh, NULL, physics_options);

        /* compute nonlinear viscosity */
        slabs_physics_compute_viscosity (
            viscosity_dim, dvisc_dIIe, NULL, NULL, NULL,
            state, press_elem, vel_dir, physics_options, 0);

        ymir_vec_destroy (dvisc_dIIe);
        ymir_vel_dir_destroy (vel_dir);
      }
    }
  }

  if (out_visc) { /* if output of viscosity */
    ymir_vec_t         *visc_surf = ymir_face_dvec_new (mesh, SL_TOP, 1,
                                                        YMIR_GAUSS_NODE);

    /* interpolate fields to surface */
    ymir_interp_vec (viscosity_dim, visc_surf);

    /* restrict surface viscosity to valid range */
    if (0.0 < physics_options->viscosity_min) {
      ymir_vec_bound_min (visc_surf, physics_options->viscosity_min);
    }
    else {
      ymir_vec_bound_min (visc_surf, ymir_vec_min_global (viscosity_dim));
    }
    if (0.0 < physics_options->viscosity_max) {
      ymir_vec_bound_max (visc_surf, physics_options->viscosity_max);
    }
    else {
      ymir_vec_bound_max (visc_surf, ymir_vec_max_global (viscosity_dim));
    }

    /* transform variables to physical dimensions */
    slabs_physics_transform_to_dimensional_viscosity (viscosity_dim);
    slabs_physics_transform_to_dimensional_viscosity (visc_surf);
    slabs_physics_transform_to_dimensional_stress (stress_norm_dim);

    /* write vtk file */
    if (out_vel_press) {
      ymir_vtk_write (mesh, path, strain_rate_dim, "strain_rate",
                      viscosity_dim, "viscosity",
                      u_surf, "u_surf", visc_surf, "visc_surf",
                      stress_norm_dim, "stress_normal",
                      NULL);
    }
    else {
      ymir_vtk_write (mesh, path, viscosity_dim, "viscosity",
                      visc_surf, "visc_surf", NULL);
    }

    /* destroy */
    ymir_vec_destroy (visc_surf);
  }

  /* destroy */
  if (u_surf != NULL) {
    ymir_vec_destroy (u_surf);
  }
  if (strain_rate_dim != NULL) {
    ymir_vec_destroy (strain_rate_dim);
  }
  if (stress_norm_dim != NULL) {
    ymir_vec_destroy (stress_norm_dim);
  }
  if (coeff != NULL) {
    ymir_vec_destroy (coeff);
  }
  if (viscosity_dim != NULL) {
    ymir_vec_destroy (viscosity_dim);
  }
}

static void
slabs_vtk_write_solution (const char *filepath,
                          slabs_stokes_state_t *state,
                          ymir_vec_t *viscosity,
                          ymir_pressure_elem_t *press_elem,
                          ymir_vec_t *rhs_u_point,
                          ymir_stokes_pc_t *stokes_pc,
                          const slabs_krylov_type_t krylov_type,
                          slabs_physics_options_t *physics_options,
                          slabs_discr_options_t *discr_options)
{
  const char         *this_fn_name = "slabs_vtk_write_solution";
  ymir_mesh_t        *mesh;
  char                path[BUFSIZ];

  YMIR_GLOBAL_INFOF ("Into %s: Write VTK file \"%s\"\n",
                     this_fn_name, filepath);

  /* check input */
  YMIR_ASSERT (viscosity != NULL);
  YMIR_ASSERT (press_elem != NULL);

  /* set mesh */
  mesh = viscosity->mesh;

  /*
   * Write main fields
   */

  {
    ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (
                                          mesh->ma->N, mesh->ma->ompsize);
    ymir_vec_t         *u_dim = ymir_cvec_new (mesh, 3);
    ymir_vec_t         *p_dim = ymir_pressure_vec_new (mesh, press_elem);
    ymir_vec_t         *viscosity_dim = ymir_vec_clone (viscosity);
    ymir_vec_t         *strain_rate_dim = ymir_dvec_new (mesh, 1,
                                                         YMIR_GAUSS_NODE);
    ymir_vec_t         *u_surf = ymir_face_cvec_new (mesh, SL_TOP, 3);
    ymir_vec_t         *strain_surf = ymir_face_dvec_new (mesh, SL_TOP, 1,
                                                          YMIR_GAUSS_NODE);
    ymir_vec_t         *visc_surf = ymir_face_dvec_new (mesh, SL_TOP, 1,
                                                        YMIR_GAUSS_NODE);

    /* get velocity and pressure from solution */
    ymir_stokes_vec_get_components (state->vel_press_vec, u_dim, p_dim,
                                    press_elem);

    /* compute 2nd invariant of the strain rate */
    ymir_second_invariant_vec (u_dim, strain_rate_dim, vel_elem);
    ymir_dvec_sqrt (strain_rate_dim, strain_rate_dim);

    /* interpolate fields to surface */
    ymir_interp_vec (u_dim, u_surf);
    ymir_interp_vec (strain_rate_dim, strain_surf);
    ymir_interp_vec (viscosity, visc_surf);

    /* restrict surface 2nd invariant to range of volume 2nd invariant */
    ymir_vec_bound_min (strain_surf, ymir_vec_min_global (strain_rate_dim));
    ymir_vec_bound_max (strain_surf, ymir_vec_max_global (strain_rate_dim));

    /* restrict surface viscosity to valid range */
    if (0.0 < physics_options->viscosity_min) {
      ymir_vec_bound_min (visc_surf, physics_options->viscosity_min);
    }
    else {
      ymir_vec_bound_min (visc_surf, ymir_vec_min_global (viscosity));
    }
    if (0.0 < physics_options->viscosity_max) {
      ymir_vec_bound_max (visc_surf, physics_options->viscosity_max);
    }
    else {
      ymir_vec_bound_max (visc_surf, ymir_vec_max_global (viscosity));
    }

    /* transform variables to physical dimensions */
    slabs_physics_transform_to_dimensional_velocity (u_dim);
    slabs_physics_transform_to_dimensional_velocity (u_surf);
    slabs_physics_transform_to_dimensional_pressure (p_dim);
    slabs_physics_transform_to_dimensional_viscosity (viscosity_dim);
    slabs_physics_transform_to_dimensional_viscosity (visc_surf);
    slabs_physics_transform_to_dimensional_strain_rate (strain_rate_dim);
    slabs_physics_transform_to_dimensional_strain_rate (strain_surf);

    /* set vtk path */
    snprintf (path, BUFSIZ, "%s_main", filepath);

    if (rhs_u_point == NULL) {
      /* write vtk file */
      ymir_vtk_write (
          mesh, path,
          u_dim, "velocity", p_dim, "pressure",
          strain_rate_dim, "strain_rate", viscosity_dim, "viscosity",
          u_surf, "u_surf", strain_surf, "strain_rate_surf",
          visc_surf, "visc_surf", NULL);
    }
    else {
      ymir_vec_t         *stress_norm = ymir_face_cvec_new (mesh, SL_TOP, 1);

      /* compute surface stress */
      YMIR_ASSERT (stokes_pc != NULL);
      slabs_physics_compute_normal_boundary_stress (
          stress_norm, state->vel_press_vec, rhs_u_point, stokes_pc->stokes_op);

      /* transform variables to physical dimensions */
      slabs_physics_transform_to_dimensional_stress (stress_norm);

      /* write vtk file */
      ymir_vtk_write (
          mesh, path,
          u_dim, "velocity", p_dim, "pressure",
          strain_rate_dim, "strain_rate", viscosity_dim, "viscosity",
          u_surf, "u_surf", strain_surf, "strain_rate_surf",
          visc_surf, "visc_surf", stress_norm, "stress_normal", NULL);

      ymir_vec_destroy (stress_norm);
    }

    /* destroy */
    ymir_vec_destroy (u_dim);
    ymir_vec_destroy (p_dim);
    ymir_vec_destroy (viscosity_dim);
    ymir_vec_destroy (strain_rate_dim);
    ymir_vec_destroy (u_surf);
    ymir_vec_destroy (strain_surf);
    ymir_vec_destroy (visc_surf);
    ymir_velocity_elem_destroy (vel_elem);
  }

  /*
   * Write residuals
   */

  if (rhs_u_point != NULL && krylov_type != SL_KRYLOV_NONE) {
    ymir_stokes_op_t   *stokes_op;
    ymir_vec_t         *residual_up = ymir_stokes_vec_new (mesh, press_elem);
    ymir_vec_t         *residual_lump_up =
                          ymir_stokes_vec_new (mesh, press_elem);
    ymir_vec_t         *residual_lump_u = ymir_cvec_new (mesh, 3);
    ymir_vec_t         *residual_lump_p = ymir_pressure_vec_new (mesh,
                                                                 press_elem);

    /* set Stokes operator */
    YMIR_ASSERT (stokes_pc != NULL);
    stokes_op = stokes_pc->stokes_op;

    /* compute residual vector */
    slabs_norm_compute_residual_l2 (residual_up, state->vel_press_vec,
                                    rhs_u_point, stokes_op, stokes_pc,
                                    krylov_type);

    /* compute weighted residual */
    ymir_vec_copy (residual_up, residual_lump_up);
    slabs_norm_weight_residual (residual_lump_up, press_elem);
    ymir_stokes_vec_get_components (residual_lump_up, residual_lump_u,
                                    residual_lump_p, press_elem);

    /* set vtk path */
    snprintf (path, BUFSIZ, "%s_residual", filepath);

    /* write vtk file */
    ymir_vtk_write (
        mesh, path,
        residual_lump_u, "residual_mom",
        residual_lump_p, "residual_mass", NULL);

    /* destroy vectors */
    ymir_vec_destroy (residual_up);
    ymir_vec_destroy (residual_lump_up);
    ymir_vec_destroy (residual_lump_u);
    ymir_vec_destroy (residual_lump_p);
  }

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

void
slabs_vtk_write_solution_lin_stokes (const char *filepath,
                                     slabs_stokes_state_t *state,
                                     slabs_lin_stokes_problem_t *lin_stokes,
                                     const slabs_krylov_type_t krylov_type,
                                     slabs_physics_options_t *physics_options,
                                     slabs_discr_options_t *discr_options)
{
  ymir_vec_t         *viscosity = ymir_vec_clone (lin_stokes->viscosity);
  ymir_pressure_elem_t  *press_elem = lin_stokes->press_elem;
  ymir_vec_t         *rhs_u_point = lin_stokes->rhs_u_point;
  ymir_stokes_pc_t   *stokes_pc = lin_stokes->stokes_pc;

  /* scale viscosity to get the actual viscosity */
  ymir_dvec_scale (0.5, viscosity);

  /* write vtk */
  slabs_vtk_write_solution (filepath, state, viscosity, press_elem,
                            rhs_u_point, stokes_pc, krylov_type,
                            physics_options, discr_options);

  /* destroy */
  ymir_vec_destroy (viscosity);
}

void
slabs_vtk_write_solution_nl_stokes (const char *filepath,
                                    slabs_stokes_state_t *state,
                                    slabs_nl_stokes_problem_t *nl_stokes,
                                    slabs_physics_options_t *physics_options,
                                    slabs_discr_options_t *discr_options)
{
  ymir_vec_t         *viscosity = ymir_vec_clone (nl_stokes->viscosity);
  ymir_pressure_elem_t  *press_elem = nl_stokes->press_elem;
  ymir_vec_t         *rhs_u_point = nl_stokes->rhs_u_point;
  ymir_stokes_pc_t   *stokes_pc = nl_stokes->stokes_pc;

  /* scale viscosity to get the actual viscosity */
  ymir_dvec_scale (0.5, viscosity);

  /* write vtk */
  slabs_vtk_write_solution (filepath, state, viscosity, press_elem,
                            rhs_u_point, stokes_pc, SL_KRYLOV_NONE,
                            physics_options, discr_options);

  /* destroy */
  ymir_vec_destroy (viscosity);
}

void
slabs_vtk_write_iteration_nl_stokes (const char *filepath,
                                     slabs_stokes_state_t *state,
                                     ymir_vec_t *step_up,
                                     ymir_vec_t *step_dual,
                                     ymir_vec_t *residual_lump_up,
                                     ymir_vec_t *residual_dual,
                                     slabs_nl_stokes_problem_t *nl_stokes,
                                     slabs_physics_options_t *physics_options,
                                     slabs_discr_options_t *discr_options)
{
  const char         *this_fn_name = "slabs_vtk_write_iteration_nl_stokes";
  ymir_vec_t         *up = state->vel_press_vec;
  ymir_vec_t         *dual = state->dual_vec;
  ymir_mesh_t        *mesh = nl_stokes->mesh;
  ymir_pressure_elem_t  *press_elem = nl_stokes->press_elem;
  ymir_stokes_op_t   *stokes_op = nl_stokes->stokes_op;
  ymir_dvec_t        *viscosity =
                        ymir_vec_clone (stokes_op->stress_op->viscosity);
//ymir_dvec_t        *dvisc_dIIe = stokes_op->stress_op->dvdIIe;
  ymir_dvec_t        *bounds_marker = nl_stokes->bounds_marker;
  ymir_dvec_t        *yielding_marker = nl_stokes->yielding_marker;
  ymir_cvec_t        *u, *step_u, *residual_lump_u;
  ymir_evec_t        *p, *step_p, *residual_lump_p;
  ymir_dvec_t        *IIe = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
  ymir_dvec_t        *grad_visc = ymir_dvec_new (mesh, 3, YMIR_GAUSS_NODE);
  ymir_dvec_t        *peclet_no = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
  ymir_dvec_t        *visc_dr;
//ymir_dvec_t        *dvisc_dIIe_dr;
  ymir_dvec_t        *bounds_diff;
  ymir_cvec_t        *uscale = ymir_cvec_new (mesh, 3);
  ymir_vec_t         *pscale = ymir_pressure_vec_new (mesh, press_elem);
  ymir_dvec_t        *amr_visc_indicator_vec;
  ymir_dvec_t        *amr_strain_rate_indicator_vec;

  YMIR_GLOBAL_INFOF ("Into %s: Write VTK file \"%s\"\n",
                     this_fn_name, filepath);

  /* scale viscosity to get the physical viscosity */
  ymir_dvec_scale (0.5, viscosity);

  /* get view of velocity and pressure components */
  slabs_stokes_vec_get_components_view (&u, &p, up);
  slabs_stokes_vec_get_components_view (&step_u, &step_p, step_up);
  slabs_stokes_vec_get_components_view (&residual_lump_u, &residual_lump_p,
                                        residual_lump_up);

  /* compute 2nd invariant of the strain rate */
  {
    ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (
                                          mesh->ma->N, mesh->ma->ompsize);

    ymir_second_invariant_vec (u, IIe, vel_elem);
    ymir_velocity_elem_destroy (vel_elem);
  }

  /* compute approx. gradient of viscosity and Peclet number */
  {
    const sc_bint_t     m = grad_visc->dvec->m;
    const sc_bint_t     n = grad_visc->dvec->n;
    sc_dmatrix_t       *visc_mat = sc_dmatrix_clone (viscosity->dvec);
    ymir_dvec_t        *lump = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);

    /* compute gradient of viscosity (let's accept the error that's produced
     * because the viscosity is defined on Gauss nodes and not GLL nodes) */
    sc_dmatrix_reshape (visc_mat, viscosity->K * viscosity->Np, 1);
    sc_dmatrix_reshape (grad_visc->dvec, viscosity->K * viscosity->Np, 3);
    slabs_gradient_gll_to_gauss (visc_mat, grad_visc->dvec, mesh->ma);
    sc_dmatrix_reshape (grad_visc->dvec, m, n);

    /* compute Peclet number */
    slabs_dvec_compute_magnitude (grad_visc, peclet_no);
    ymir_dvec_divide_in (viscosity, peclet_no);
    ymir_mass_lump (lump);
    ymir_dvec_multiply_in (lump, peclet_no);
    ymir_dvec_scale (discr_options->domain_size_normalization, peclet_no);

    /* destroy */
    sc_dmatrix_destroy (visc_mat);
    ymir_vec_destroy (lump);
  }

  /* compute dynamic range of viscosity per element */
  {
    sc_dmatrix_t       *visc_dr_mat;

    visc_dr_mat = slabs_compute_dynamic_range_per_element (viscosity);
    visc_dr = slabs_dvec_new_from_element_data (visc_dr_mat, mesh,
                                                YMIR_GAUSS_NODE);
    sc_dmatrix_destroy (visc_dr_mat);
  }
//{
//  sc_dmatrix_t       *visc_dr_mat, *dvisc_dIIe_dr_mat;

//  visc_dr_mat = slabs_compute_dynamic_range_per_element (viscosity);
//  visc_dr = slabs_dvec_new_from_element_data (visc_dr_mat, mesh,
//                                              YMIR_GAUSS_NODE);

//  dvisc_dIIe_dr_mat = slabs_compute_dynamic_range_per_element (dvisc_dIIe);
//  dvisc_dIIe_dr = slabs_dvec_new_from_element_data (dvisc_dIIe_dr_mat, mesh,
//                                                    YMIR_GAUSS_NODE);

//  /* destroy */
//  sc_dmatrix_destroy (visc_dr_mat);
//  sc_dmatrix_destroy (dvisc_dIIe_dr_mat);
//}

  /* compute max difference of bounds per element */
  {
    sc_dmatrix_t       *bounds_diff_mat;

    bounds_diff_mat = slabs_compute_max_difference_per_element (bounds_marker);
    bounds_diff = slabs_dvec_new_from_element_data (bounds_diff_mat, mesh,
                                                    YMIR_GAUSS_NODE);

    sc_dmatrix_destroy (bounds_diff_mat);
  }

  /* compute velocity scaling */
  if (ymir_stokes_pc_set_uscale_fn != NULL) {
    ymir_stokes_pc_set_uscale_fn (uscale, ymir_stokes_pc_set_uscale_fn_data);
  }
  else {
    ymir_vec_set_zero (uscale);
  }

  /* compute pressure scaling */
  if (ymir_stokes_pc_set_pscale_fn != NULL) {
    ymir_stokes_pc_set_pscale_fn (pscale, ymir_stokes_pc_set_pscale_fn_data);
  }
  else {
    ymir_vec_set_zero (pscale);
  }

  {
    const slabs_discr_amr_indicator_type_t  visc_indicator_type =
      discr_options->amr_visc_indicator_type;
    const slabs_discr_amr_indicator_type_t  strain_rate_indicator_type =
      discr_options->amr_strain_rate_indicator_type;

    /* compute AMR indictor based on viscosity */
    if (visc_indicator_type != SL_AMR_INDICATOR_NONE) {
      amr_visc_indicator_vec = slabs_vtk_create_vec_from_indicator (
          visc_indicator_type, state, mesh, physics_options, discr_options,
          0 /* !init_amr */, 0 /* import_data? */);
    }
    else {
      amr_visc_indicator_vec = ymir_dvec_new_zero (mesh, 1, YMIR_GAUSS_NODE);
    }

    /* compute AMR indictor based on strain rate */
    if (strain_rate_indicator_type != SL_AMR_INDICATOR_NONE) {
      amr_strain_rate_indicator_vec = slabs_vtk_create_vec_from_indicator (
          strain_rate_indicator_type, state, mesh,
          physics_options, discr_options,
          0 /* !init_amr */, 0 /* import_data? */);
    }
    else {
      amr_strain_rate_indicator_vec = ymir_dvec_new_zero (mesh, 1,
                                                          YMIR_GAUSS_NODE);
    }
  }

  /* write vtk file */
  if (   dual != NULL
      && residual_dual != NULL
      && step_dual != NULL ) {
    ymir_dvec_t        *dual_magn = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
    ymir_dvec_t        *step_dual_magn = ymir_dvec_new (mesh, 1,
                                                        YMIR_GAUSS_NODE);
    ymir_dvec_t        *residual_dual_magn = ymir_dvec_new (mesh, 1,
                                                            YMIR_GAUSS_NODE);

    slabs_dvec_compute_magnitude (dual, dual_magn);
    slabs_dvec_compute_magnitude (step_dual, step_dual_magn);
    slabs_dvec_compute_magnitude (residual_dual, residual_dual_magn);

    ymir_vtk_write (
        mesh, filepath,
        u, "u", p, "p", step_u, "step_u", step_p, "step_p", IIe, "IIe",
        residual_lump_u, "residual_lump_mom",
        residual_lump_p, "residual_lump_mass",
        dual_magn, "dual_magn", step_dual_magn, "step_dual_magn",
        residual_dual_magn, "residual_dual_magn",
        viscosity, "viscosity", visc_dr, "visc_dr",
        grad_visc, "grad_visc", peclet_no, "peclet_no_approx",
      //dvisc_dIIe, "dvisc_dIIe", dvisc_dIIe_dr, "dvisc_dIIe_dr",
        bounds_marker, "bounds_marker", bounds_diff, "bounds_diff",
        yielding_marker, "yielding_marker",
        uscale, "uscale", pscale, "pscale",
        amr_visc_indicator_vec, "amr_visc_indicator",
        amr_strain_rate_indicator_vec, "amr_strain_rate_indicator",
        NULL);

    ymir_vec_destroy (dual_magn);
    ymir_vec_destroy (step_dual_magn);
    ymir_vec_destroy (residual_dual_magn);
  }
  else {
    ymir_vtk_write (
        mesh, filepath,
        u, "u", p, "p", step_u, "step_u", step_p, "step_p", IIe, "IIe",
        residual_lump_u, "residual_lump_mom",
        residual_lump_p, "residual_lump_mass",
        viscosity, "viscosity", visc_dr, "visc_dr",
        grad_visc, "grad_visc", peclet_no, "peclet_no_approx",
      //dvisc_dIIe, "dvisc_dIIe", dvisc_dIIe_dr, "dvisc_dIIe_dr",
        bounds_marker, "bounds_marker", bounds_diff, "bounds_diff",
        yielding_marker, "yielding_marker",
        uscale, "uscale", pscale, "pscale",
        amr_visc_indicator_vec, "amr_visc_indicator",
        amr_strain_rate_indicator_vec, "amr_strain_rate_indicator",
        NULL);
  }

  /* destroy vectors */
  ymir_vec_destroy (viscosity);
  ymir_vec_destroy (u);
  ymir_vec_destroy (p);
  ymir_vec_destroy (step_u);
  ymir_vec_destroy (step_p);
  ymir_vec_destroy (residual_lump_u);
  ymir_vec_destroy (residual_lump_p);
  ymir_vec_destroy (IIe);
  ymir_vec_destroy (grad_visc);
  ymir_vec_destroy (peclet_no);
  ymir_vec_destroy (visc_dr);
//ymir_vec_destroy (dvisc_dIIe_dr);
  ymir_vec_destroy (bounds_diff);
  ymir_vec_destroy (uscale);
  ymir_vec_destroy (pscale);
  ymir_vec_destroy (amr_visc_indicator_vec);
  ymir_vec_destroy (amr_strain_rate_indicator_vec);

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

