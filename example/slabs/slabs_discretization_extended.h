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

#ifndef SLABS_DISCRETIZATION_EXTENDED_H
#define SLABS_DISCRETIZATION_EXTENDED_H

#include <slabs_discretization.h>

/* AMR indicator */
typedef struct slabs_discr_amr_indicator
{
  sc_dmatrix_t       *val;        /* values of indicator (#cols = #quads) */
  slabs_discr_amr_indicator_type_t  type;  /* type of indicator */
  double              tol_min;    /* coarsen below this value */
  double              tol_max;    /* refine above this value */
  int8_t              level_min;  /* min level for coarsening */
  int8_t              level_max;  /* max level for refinement */
}
slabs_discr_amr_indicator_t;

/* enumerator for AMR markers flags */
typedef enum
{
  SLABS_DISCR_AMR_MARKER_COARSEN = -1,
  SLABS_DISCR_AMR_MARKER_NONE    =  0,
  SLABS_DISCR_AMR_MARKER_REFINE  =  1
}
slabs_discr_amr_marker_flag_t;

/* AMR marker */
typedef struct slabs_discr_amr_marker
{
  slabs_discr_amr_marker_flag_t *flag;  /* refine/coarsen flags */
  p4est_locidx_t      n_coarsen_quadrants_loc;
  p4est_locidx_t      n_refine_quadrants_loc;

  p4est_t            *p4est;  /* associated p4est (not owned) */
}
slabs_discr_amr_marker_t;

/* user data of a p4est quadrant for AMR */
typedef struct
{
  int                 flag;
}
slabs_discr_p4est_qdata_t;

/*
 * Refines mesh radially (for spherical geometry) or in z-direction (for cubes).
 */
void
slabs_discr_refine_wrt_depth (p8est_t *p8est,
                              slabs_physics_options_t *physics_options,
                              slabs_discr_options_t *discr_options);

/**
 * Refines mesh about a layer that is located at a certain depth.
 */
void
slabs_discr_refine_at_layer (p8est_t *p8est,
                             slabs_physics_options_t *physics_options,
                             slabs_discr_options_t *discr_options);

/**
 * Refines mesh at the surface.
 */
void
slabs_discr_refine_at_surface (p8est_t *p8est,
                               slabs_physics_options_t *physics_options,
                               slabs_discr_options_t *discr_options);

/**
 * Creates a new AMR indicator.
 */
slabs_discr_amr_indicator_t *
slabs_discr_amr_indicator_new (slabs_discr_amr_indicator_type_t type,
                               const double tol_min, const double tol_max,
                               const int8_t level_min, const int8_t level_max,
                               const p4est_locidx_t n_quadrants_loc);

/**
 * Destroys an AMR indicator.
 */
void
slabs_discr_amr_indicator_destroy (slabs_discr_amr_indicator_t *indicator);

/**
 * Creates new initial AMR indicator values. TODO use _indicator_t structure
 */
sc_dmatrix_t *
slabs_discr_init_amr_indicator_new (mangll_t *mangll,
                                    slabs_physics_options_t *physics_options,
                                    slabs_discr_options_t *discr_options,
                                    int indicator_type);

/**
 * Sets values of AMR indicator.
 */
void
slabs_discr_amr_indicator_set (slabs_discr_amr_indicator_t *indicator,
                               slabs_stokes_state_t *state,
                               mangll_t *mangll, mangll_cnodes_t *cnodes,
                               slabs_physics_options_t *physics_options,
                               slabs_discr_options_t *discr_options,
                               const int init_amr);

/**
 * Computes max difference per element.
 */
sc_dmatrix_t *
slabs_compute_max_difference_per_element (ymir_dvec_t *vec);

/**
 * Computes dynamic range per element.
 */
sc_dmatrix_t *
slabs_compute_dynamic_range_per_element (ymir_dvec_t *vec);

#endif /* SLABS_DISCRETIZATION_EXTENDED_H */
