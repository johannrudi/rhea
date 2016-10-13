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

#ifndef SLABS_DISCRETIZATION_H
#define SLABS_DISCRETIZATION_H

#include <p8est.h>
#include <p4est_to_p8est.h>
#include <slabs_stokes_state.h>
#include <slabs_physics.h>

/* enumerator for AMR indicator types, refinement w.r.t. viscosity */
typedef enum
{
  /* [0] general */
  SL_AMR_INDICATOR_NONE        = 0,
  SL_AMR_INDICATOR_COARSEN_ALL = 1,
  SL_AMR_INDICATOR_REFINE_ALL  = 2,

  /* [3] viscosity */
  SL_AMR_INDICATOR_VISC_DR            = 3,
  SL_AMR_INDICATOR_VISC_GRAD          = 4,
  SL_AMR_INDICATOR_VISC_PECLET        = 5,
  SL_AMR_INDICATOR_VISC_TENSOR_PECLET = 6,
  SL_AMR_INDICATOR_RHS_OVER_VISC      = 7,
  SL_AMR_INDICATOR_RHS_OVER_VISC_DR   = 8,
  SL_AMR_INDICATOR_RHS_OVER_VISC_GRAD = 9,

  /* [10] strain rate */
  SL_AMR_INDICATOR_STRAIN_RATE_PECLET = 10,

  /* [11] weak zone */
  SL_AMR_INDICATOR_WEAK_DR          = 11,
  SL_AMR_INDICATOR_WEAK_GRAD        = 12,
  SL_AMR_INDICATOR_WEAK_SUBDU_DIST  = 13,
  SL_AMR_INDICATOR_WEAK_RIDGE_DIST  = 14,
  SL_AMR_INDICATOR_WEAK_IMPORT_DIST = 15,

  /* [16] right-hand side */
  SL_AMR_INDICATOR_RHS_MAGNITUDE  = 16,
  SL_AMR_INDICATOR_RHS_PECLET     = 17,

  /* [18] misc */
  SL_AMR_INDICATOR_REFINE_SURFACE = 18,
  SL_AMR_INDICATOR_REFINE_LAYER   = 19
}
slabs_discr_amr_indicator_type_t;

/* parameters of a set of AMR indicators */
typedef struct slabs_discr_amr_indicator_params
{
  int                 n_indicators;
  slabs_discr_amr_indicator_type_t  *type;
  double             *tol_min;    /* coarsen below this value */
  double             *tol_max;    /* refine above this value */
  int8_t             *level_min;  /* min level for coarsening */
  int8_t             *level_max;  /* max level for refinement */
}
slabs_discr_amr_indicator_params_t;

/* enumerator for partitioning types */
typedef enum
{
  SL_MESH_PARTITIONING_ELEM,
  SL_MESH_PARTITIONING_DOF_VEL,
  SL_MESH_PARTITIONING_DOF_VEL_PRESS
}
slabs_mesh_partitioning_type_t;

/* options for the discretization */
typedef struct slabs_discr_options
{
  /* boundary information */
  ymir_mesh_e_to_fm_t e_to_fm_fn;
  ymir_topidx_t      *tree_to_bf; /* data for `e_to_fm_fn` */

  /* transformation function from element space to physical space */
  mangll_X_t          X_fn;

  /* basis fuctions properties */
  int                 order;
  int                 n_vel_dnodes_per_el;
  int                 n_vel_dnodes_per_el_interior;
  int                 n_vel_dnodes_per_face_interior;
  int                 n_vel_dnodes_per_edge_interior;
  int                 n_vel_dnodes_per_corner;

  int                 order_press;
  int                 n_press_dnodes_per_el;

  /* refinement levels */
  int8_t              minlevel;
  int8_t              maxlevel;

  /* refinement type for initial refinement */
  char               *refine;

  /* radii for which the domain is increasingly refined */
  double             *refine_radius;
  int                 refine_n_radii;

  /* mesh parameters for data import */
  int                 import_mesh_order;
  int                 import_mesh_minlevel;

  /* refinement at the surface */
  double              refine_surface_maxdist;
  int8_t              refine_surface_maxlevel;
  int                 refine_surface_n_trees;

  /* refinement about a layer at a fixed radius */
  double              refine_layer_radius;
  double              refine_layer_maxdist;
  int8_t              refine_layer_maxlevel;
  int                 refine_layer_n_trees;
  int                 enforce_refinement_at_layer;

  /* domain AMR normalization factor */
  double              domain_size_normalization;

  /* initial AMR */
  int                 init_amr_max_steps;       /* max num of AMR cycles     */
  double              init_amr_rel_threshold;   /* min #quads marked for AMR */
  double              init_amr_n_elements_max;  /* max #quadrants in mesh    */
  int                 init_amr_override_order;  /* override poly discr order */
  int                 init_amr_lower_mantle;    /* init AMR in lower mantle  */

  slabs_discr_amr_indicator_type_t  init_amr_visc_indicator_type;
  double              init_amr_visc_tol_min;  /* coarsen below this value */
  double              init_amr_visc_tol_max;  /* refine above this value  */
  double              init_amr_visc_inside_plates_tol_min;
  double              init_amr_visc_inside_plates_tol_max;

  slabs_discr_amr_indicator_type_t  init_amr_weak_subdu_indicator_type;
  double              init_amr_weak_subdu_tol_min;
  double              init_amr_weak_subdu_tol_max;
  double              init_amr_weak_subdu_maxdist;
  int8_t              init_amr_weak_subdu_maxlevel;

  slabs_discr_amr_indicator_type_t  init_amr_weak_ridge_indicator_type;
  double              init_amr_weak_ridge_tol_min;
  double              init_amr_weak_ridge_tol_max;
  double              init_amr_weak_ridge_maxdist;
  int8_t              init_amr_weak_ridge_maxlevel;

  slabs_discr_amr_indicator_type_t  init_amr_weak_import_indicator_type;
  double              init_amr_weak_import_tol_min;
  double              init_amr_weak_import_tol_max;
  double              init_amr_weak_import_maxdist;
  int8_t              init_amr_weak_import_maxlevel;

  slabs_discr_amr_indicator_type_t  init_amr_rhs_indicator_type;
  double              init_amr_rhs_tol_min;
  double              init_amr_rhs_tol_max;
  double              init_amr_rhs_norm_shift;

  int                 init_amr_post_uniform_n_steps;

  /* AMR */
  int                 amr_max_steps;       /* max num of AMR cycles         */
  double              amr_rel_threshold;   /* min #quadrants marked for AMR */
  double              amr_n_elements_max;  /* max #quadrants in mesh        */
  int                 amr_lower_mantle;    /* AMR in lower mantle           */

  slabs_discr_amr_indicator_type_t  amr_visc_indicator_type;
  double              amr_visc_tol_min;  /* coarsen below this value */
  double              amr_visc_tol_max;  /* refine above this value */

  slabs_discr_amr_indicator_type_t  amr_visc_dr_indicator_type;
  double              amr_visc_dr_tol_min;
  double              amr_visc_dr_tol_max;

  slabs_discr_amr_indicator_type_t  amr_strain_rate_indicator_type;
  double              amr_strain_rate_tol_min;
  double              amr_strain_rate_tol_max;

  int                 amr_log_maxlevel;

  /* type of partitioning of the mesh among processsors */
  slabs_mesh_partitioning_type_t  mesh_partitioning_type;

  int                 inspect_p4est;
}
slabs_discr_options_t;

/**
 *
 */
void
slabs_discr_set_order (slabs_discr_options_t *discr_options, const int N);

/**
 *
 */
slabs_discr_amr_indicator_params_t *
slabs_discr_amr_indicator_params_new (const int n_indicators);

/**
 *
 */
void
slabs_discr_amr_indicator_params_destroy (slabs_discr_amr_indicator_params_t
                                           *indicator_params);

/**
 * Computes the level required to reach a desired resolution.
 *
 * \param [in] element_resolution  Desired side length of single elem at surface
 * \return                         level
 */
double
slabs_discr_resolution_to_level (double element_density,
                                 slabs_physics_options_t *physics_options);

/**
 * Computes width and depth of weak zone in left corner (at mid ocean ridge)
 * such that the boundaries of the weak zone are aligned with element
 * boundaries.
 */
void
slabs_discr_align_weak_ridge_2plates (slabs_physics_options_t *physics_options,
                                      slabs_discr_options_t *discr_options);

/**
 * Computes radius of interface between lower and upper mantle such that the
 * interface is aligned with element boundaries.
 */
void
slabs_discr_set_upper_mantle_radius (slabs_physics_options_t *physics_options,
                                     slabs_discr_options_t *discr_options);

/**
 * Computes a normalization factor that ought to make mesh elements with the
 * same refinement level comparable in size.
 */
void
slabs_discr_compute_domain_size_normalization (slabs_discr_options_t
                                                 *discr_options,
                                               slabs_physics_options_t
                                                 *physics_options);

/**
 * Identity geometry transformation.
 */
void
slabs_discr_identity_X (mangll_tag_t tag, mangll_locidx_t np,
                        const double *_sc_restrict EX,
                        const double *_sc_restrict EY,
                        const double *_sc_restrict EZ,
                        double *_sc_restrict X,
                        double *_sc_restrict Y,
                        double *_sc_restrict Z, void *data);

/**
 * Geometry transformation for 24-octree hollow sphere.
 */
void
slabs_discr_shell_X (mangll_tag_t tag, mangll_locidx_t np,
                     const double *_sc_restrict EX,
                     const double *_sc_restrict EY,
                     const double *_sc_restrict EZ,
                     double *_sc_restrict X,
                     double *_sc_restrict Y,
                     double *_sc_restrict Z, void *data);

/**
 * Geometry transformation for chunk and slice of hollow sphere.
 */
void
slabs_discr_shell_chunk_and_slice_X (mangll_tag_t tag, mangll_locidx_t np,
                                     const double *_sc_restrict EX,
                                     const double *_sc_restrict EY,
                                     const double *_sc_restrict EZ,
                                     double *_sc_restrict X,
                                     double *_sc_restrict Y,
                                     double *_sc_restrict Z, void *data);

/**
 * Creates new p4est object.
 */
p8est_t *
slabs_discr_p8est_new (MPI_Comm mpicomm,
                       slabs_physics_options_t *physics_options,
                       slabs_discr_options_t *discr_options);

/**
 * Destroys p4est object.
 */
void
slabs_discr_p8est_destroy (p8est_t *p8est);

/**
 * Initializes user data of a p4est object.
 */
void
slabs_discr_p8est_init_data (p8est_t *p8est, const int inspect_p4est);

/**
 * Clears user data from a p4est object.
 */
void
slabs_discr_p8est_clear_data (p8est_t *p8est);

/**
 *
 */
void
slabs_discr_options_set_boundary (slabs_discr_options_t *discr_options,
                                  p8est_t *p8est,
                                  slabs_physics_options_t *physics_options);

/**
 *
 */
void
slabs_discr_options_clear_boundary (slabs_discr_options_t *discr_options);

/**
 * Creates new mangll mesh structure.
 */
mangll_mesh_t *
slabs_discr_mangll_mesh_new (p8est_t *p8est,
                             slabs_discr_options_t *discr_options);

/**
 * Creates new mangll and cnodes structures.
 */
void
slabs_discr_mangll_and_cnodes_new (mangll_t **mangll,
                                   mangll_cnodes_t **cnodes,
                                   p8est_t *p8est,
                                   slabs_discr_options_t *discr_options);

/**
 * Creates new ymir mesh and pressure element structures.
 */
void
slabs_discr_ymir_new (ymir_mesh_t **mesh, ymir_pressure_elem_t **press_elem,
                      mangll_t *mangll, mangll_cnodes_t *cnodes,
                      slabs_discr_options_t *discr_options);

/**
 * Destroys ymir mesh and corresponding mangll structures.
 */
void
slabs_discr_ymir_mangll_destroy (ymir_mesh_t *mesh,
                                 ymir_pressure_elem_t *press_elem);

/* data for coarsening the Stokes coefficient */
typedef struct slabs_discr_enforce_refinement_data_t
{
  slabs_physics_options_t *physics_options;
  slabs_discr_options_t   *discr_options;
}
slabs_discr_enforce_refinement_data_t;

/**
 *
 */
slabs_discr_enforce_refinement_data_t *
slabs_discr_enforce_refinement_data_new (
                                      slabs_physics_options_t *physics_options,
                                      slabs_discr_options_t *discr_options);

/**
 *
 */
void
slabs_discr_enforce_refinement_data_destroy (
                                  slabs_discr_enforce_refinement_data_t *data);

/**
 * Enforces (keeps) mesh refinement about a layer that is located at a certain
 * depth.
 */
void
slabs_discr_enforce_refinement_at_layer (p8est_t *p8est, void *data);

/**
 * Initial AMR and partitioning of p4est mesh.
 */
void
slabs_discr_initial_amr_no_interp (p8est_t *p8est,
                                   slabs_physics_options_t *physics_options,
                                   slabs_discr_options_t *discr_options);

/**
 *
 */
int
slabs_discr_amr (slabs_stokes_state_t *state,
                 ymir_mesh_t **mesh, ymir_pressure_elem_t **press_elem,
                 p8est_t *p8est,
                 slabs_discr_amr_indicator_params_t *indicator_params,
                 slabs_physics_options_t *physics_options,
                 slabs_discr_options_t *discr_options,
                 const int init_amr);

#endif /* SLABS_DISCRETIZATION_H */
