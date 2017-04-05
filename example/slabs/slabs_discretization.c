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

#include <slabs_discretization.h>
#include <slabs_discretization_extended.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <mangll_fields.h>
#include <mangll_p8est.h>
#include <mangll_tensor.h>
#include <ymir_stress_pc.h>
#include <ymir_monitor.h>
#include <slabs_physics_extended.h>

//#define SLABS_DISCR_AMR_VTK

#ifdef SLABS_DISCR_AMR_VTK
#include <slabs_vtk.h>
#endif

#define SLABS_DISCR_AMR_REL_THRESH_MIN 1.0e-4

/**
 *
 */
void
slabs_discr_set_order (slabs_discr_options_t *discr_options, const int N)
{
  /* set polynomial order of velocity basis functions */
  discr_options->order = ymir_n (N);

  /* set number of velocity dnodes per element */
  discr_options->n_vel_dnodes_per_el = (N + 1) * (N + 1) * (N + 1);

  /* set number of velocity dnodes per element interior (volume dnodes) */
  discr_options->n_vel_dnodes_per_el_interior = N * N * N;

  /* set number of velocity dnodes per face interior (face dnodes) */
  discr_options->n_vel_dnodes_per_face_interior = N * N;

  /* set number of velocity dnodes per edge interior (edge dnodes) */
  discr_options->n_vel_dnodes_per_edge_interior = N;

  /* set number of velocity dnodes per corner (corner dnodes) */
  discr_options->n_vel_dnodes_per_corner = 1;

  /* set pressure order and number of dnodes per element */
  switch (ymir_pressure_default_space) {
  case YMIR_PRESSURE_SPACE_POLY:
    YMIR_ASSERT (1 < N);
    discr_options->order_press = N - 1;
    discr_options->n_press_dnodes_per_el = N * (N + 1) * (N + 2) / 6;
    break;

  case YMIR_PRESSURE_SPACE_TENS:
    YMIR_ASSERT (1 < N);
    discr_options->order_press = N - 0;
    discr_options->n_press_dnodes_per_el = (N - 1) * (N - 1) * (N - 1);
    break;

  case YMIR_PRESSURE_SPACE_STAB:
    YMIR_ASSERT (1 == N);
    discr_options->order_press = 1;
    discr_options->n_press_dnodes_per_el = (N + 1) * (N + 1) * (N + 1);
    break;

  default:
    YMIR_ABORT_NOT_REACHED ();
  }
}

/**
 *
 */
slabs_discr_amr_indicator_params_t *
slabs_discr_amr_indicator_params_new (const int n_indicators)
{
  slabs_discr_amr_indicator_params_t  *indicator_params;

  if (n_indicators <= 0) {
    return NULL;
  }

  indicator_params = YMIR_ALLOC (slabs_discr_amr_indicator_params_t, 1);

  indicator_params->n_indicators = n_indicators;
  indicator_params->type =
    YMIR_ALLOC (slabs_discr_amr_indicator_type_t, n_indicators);
  indicator_params->tol_min = YMIR_ALLOC (double, n_indicators);
  indicator_params->tol_max = YMIR_ALLOC (double, n_indicators);
  indicator_params->level_min = YMIR_ALLOC (int8_t, n_indicators);
  indicator_params->level_max = YMIR_ALLOC (int8_t, n_indicators);

  return indicator_params;
}

/**
 *
 */
void
slabs_discr_amr_indicator_params_destroy (slabs_discr_amr_indicator_params_t
                                           *indicator_params)
{
  if (indicator_params == NULL) {
    return;
  }

  YMIR_FREE (indicator_params->type);
  YMIR_FREE (indicator_params->tol_min);
  YMIR_FREE (indicator_params->tol_max);
  YMIR_FREE (indicator_params->level_min);
  YMIR_FREE (indicator_params->level_max);

  YMIR_FREE (indicator_params);
}

/**
 * Computes the level required to reach a desired resolution.
 */
 double
slabs_discr_resolution_to_level (double element_resolution,
                                 slabs_physics_options_t *physics_options)
{
  const double        lon_min = physics_options->domain_lon_min;
  const double        lon_max = physics_options->domain_lon_max;
  double              surf_length_per_tree;

  /* compute the length that is covert by one tree at the domain surface */
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_SHELL_CHUNK:
    /* has width pi/4 */
    surf_length_per_tree = (lon_max - lon_min) * SL_EARTH_RADIUS;
    break;

  case SL_DOMAIN_BRICK:
    /* treat like shell slice */
    surf_length_per_tree = (lon_max - lon_min) * SL_EARTH_RADIUS
                           / ((double) physics_options->domain_brick_dy);
    break;

  case SL_DOMAIN_SHELL:
    /* has 8 trees for one circle */
    surf_length_per_tree = M_PI * SL_EARTH_RADIUS / 4.0;
    break;

  case SL_DOMAIN_SHELL_SLICE:
    /* has width pi/4 */
    surf_length_per_tree = (lon_max - lon_min) * SL_EARTH_RADIUS
                           / ((double) physics_options->domain_brick_dy);
    break;

  default: /* invalid domain type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* return number of levels */
  if (1.0 <= (surf_length_per_tree / element_resolution)) {
    /* return level required to reach desired resolution */
    return log2 (surf_length_per_tree / element_resolution);
  }
  else {
    /* given resolution is coarser than resolution of one tree */
    return -1.0;
  }
}

/**
 *
 */
static void
slabs_align_cube_z_direction_with_mesh (double *z,
                                        double *z_ref_domain,
                                        double max_z, int max_level,
                                        slabs_physics_options_t
                                          *physics_options,
                                        slabs_discr_options_t *discr_options)
{
  int                 level, n;
  double              h;

  /* check input parameters */
  YMIR_ASSERT (SL_SHELL_RADIUS_BOTTOM <= max_z);
  YMIR_ASSERT (max_z <= SL_SHELL_RADIUS_TOP);

  /* set number of subdivisions given by the number of trees */
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
    n = 1;
    break;

  case SL_DOMAIN_BRICK:
    n = physics_options->domain_brick_dz;
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }
  YMIR_ASSERT (n >= 1);

  /* set level of refinement at `max_z` */
  level = SC_MAX (discr_options->minlevel, max_level);

  /* calculate element size `h` */
  h = (SL_SHELL_RADIUS_TOP - SL_SHELL_RADIUS_BOTTOM) / ((double) (n << level));

  /* calculate closest aligned z-coordinate below `max_z` */
  *z = SL_SHELL_RADIUS_BOTTOM
            + h * floor ((max_z - SL_SHELL_RADIUS_BOTTOM) / h);
  if (z_ref_domain != NULL) {
    *z_ref_domain = (*z - SL_SHELL_RADIUS_BOTTOM)
                    / (SL_SHELL_RADIUS_TOP - SL_SHELL_RADIUS_BOTTOM);
  }
}

/**
 *
 */
static void
slabs_align_shell_radius_with_mesh (double *radius, double *radius_ref_domain,
                                    double max_radius, int max_level,
                                    slabs_physics_options_t *physics_options,
                                    slabs_discr_options_t *discr_options)
{
  const double        Rt_Rb = SL_SHELL_RADIUS_TOP / SL_SHELL_RADIUS_BOTTOM;
  const double        Rbsq_Rt = SL_SHELL_RADIUS_BOTTOM * SL_SHELL_RADIUS_BOTTOM
                                / SL_SHELL_RADIUS_TOP;
  int                 level, n;
  double              h;
  double              max_radius_ref, radius_ref;

  /* check input parameters */
  YMIR_ASSERT (SL_SHELL_RADIUS_BOTTOM <= max_radius);
  YMIR_ASSERT (max_radius <= SL_SHELL_RADIUS_TOP);

  /* map max radius from physical space to reference space */
  max_radius_ref = log (max_radius / Rbsq_Rt) / log (Rt_Rb);

  /* set number of subdivisions given by the number of trees */
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_SHELL:
  case SL_DOMAIN_SHELL_CHUNK:
    n = 1;
    break;

  case SL_DOMAIN_SHELL_SLICE:
    n = physics_options->domain_brick_dz;
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }
  YMIR_ASSERT (n >= 1);

  /* set level of refinement at `max_radius_ref` */
  level = SC_MAX (discr_options->minlevel, max_level);

  /* calculate element size `h` */
  h = 1.0 / ((double) (n << level));

  /* calculate closest aligned radius below `max_radius_ref` */
  radius_ref = 1.0 + h * floor ((max_radius_ref - 1.0) / h);
  if (radius_ref_domain != NULL) {
    *radius_ref_domain = radius_ref;
  }

  /* map radius from reference space to physical space */
  *radius = Rbsq_Rt * pow (Rt_Rb, radius_ref);
}

/**
 * Computes width and depth of weak zone in left corner (at mid ocean ridge)
 * such that the boundaries of the weak zone are aligned with element
 * boundaries.
 */
void
slabs_discr_align_weak_ridge_2plates (slabs_physics_options_t *physics_options,
                                      slabs_discr_options_t *discr_options)
{
  const double        lon_min = physics_options->domain_lon_min;
  const double        lon_max = physics_options->domain_lon_max;
  int                 level = discr_options->minlevel
                              + discr_options->refine_n_radii;
  double              weakzone_max_radius_bottom;
  double              weakzone_radius_bottom;
  double              weakzone_min_longitude;
  double              weakzone_longitude;
  double              hy;

  weakzone_max_radius_bottom = SL_SHELL_RADIUS_TOP
    - physics_options->weakzone_2plates_ridge_depth / SL_EARTH_RADIUS;
  weakzone_min_longitude = lon_min
    + physics_options->weakzone_2plates_ridge_width / SL_EARTH_RADIUS
    * SL_SHELL_RADIUS_TOP;

  /* set initial element size `h` */
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_BRICK:
    hy = (lon_max - lon_min)
         / ((double) physics_options->domain_brick_dy);
    slabs_align_cube_z_direction_with_mesh (&weakzone_radius_bottom, NULL,
        weakzone_max_radius_bottom, level, physics_options, discr_options);
    break;

  case SL_DOMAIN_SHELL_SLICE:
    hy = (lon_max - lon_min)
         / ((double) physics_options->domain_brick_dy);
    slabs_align_shell_radius_with_mesh (&weakzone_radius_bottom, NULL,
        weakzone_max_radius_bottom, level, physics_options, discr_options);
    break;

  default: /* invalid domain type */
    YMIR_ABORT_NOT_REACHED ();
  }
  hy /= (double) (1 << level);

  /* compute weak zone longitude */
  weakzone_longitude = lon_min
    + hy * ceil ((weakzone_min_longitude - lon_min) / hy);

  /* set aligned weak zone depth */
  physics_options->weakzone_2plates_ridge_depth = SL_EARTH_RADIUS
    * (SL_SHELL_RADIUS_TOP - weakzone_radius_bottom);

  /* set aligned weak zone width */
  physics_options->weakzone_2plates_ridge_width = SL_EARTH_RADIUS
    * (weakzone_longitude - lon_min) / SL_SHELL_RADIUS_TOP;
}

/**
 * Computes radius of interface between lower and upper mantle such that the
 * interface is aligned with element boundaries.
 */
void
slabs_discr_set_upper_mantle_radius (slabs_physics_options_t *physics_options,
                                     slabs_discr_options_t *discr_options)
{
  const double        max_radius = SL_SHELL_RADIUS_TOP
                        - SL_UPPER_MANTLE_DEPTH / SL_EARTH_RADIUS;
  const double        refine_maxdist = discr_options->refine_layer_maxdist;
  int                 refine_maxlevel = discr_options->refine_layer_maxlevel;
  double              radius, radius_ref_domain;

  /* set to no additional refinement at lower/upper mantle interface */
  if (refine_maxdist <= 0.0) {
    refine_maxlevel = 0;
  }

  /* set element size `h` given the number of trees */
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_BRICK:
    slabs_align_cube_z_direction_with_mesh (&radius, &radius_ref_domain,
        max_radius, refine_maxlevel, physics_options, discr_options);
    break;

  case SL_DOMAIN_SHELL:
  case SL_DOMAIN_SHELL_CHUNK:
  case SL_DOMAIN_SHELL_SLICE:
    slabs_align_shell_radius_with_mesh (&radius, &radius_ref_domain,
        max_radius, refine_maxlevel, physics_options, discr_options);
    radius_ref_domain -= 1.0;
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* set aligned radius of upper mantle */
  physics_options->viscosity_upper_mantle_radius = radius;
  discr_options->refine_layer_radius = radius_ref_domain;
}

/**
 * Computes a normalization factor that ought to make mesh elements with the
 * same refinement level comparable in size.
 */
void
slabs_discr_compute_domain_size_normalization (slabs_discr_options_t
                                                 *discr_options,
                                               slabs_physics_options_t
                                                 *physics_options)
{
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_BRICK:
    discr_options->domain_size_normalization =
      1.0 / physics_options->domain_z_max;
    break;

  case SL_DOMAIN_SHELL:
  case SL_DOMAIN_SHELL_CHUNK:
  case SL_DOMAIN_SHELL_SLICE:
    discr_options->domain_size_normalization =
      1.0 / (  physics_options->domain_radius_max
             - physics_options->domain_radius_min );
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }
}

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
                        double *_sc_restrict Z, void *data)
{
  mangll_locidx_t     il;

  for (il = 0; il < np; ++il) {
    X[il] = EX[il];
    Y[il] = EY[il];
    Z[il] = EZ[il];
  }
}

#define SLABS_DISCR_SHELL_X_TAG_TOP 3
#define SLABS_DISCR_SHELL_X_TAG_LEFT 2
#define SLABS_DISCR_SHELL_X_TAG_RIGHT 0
#define SLABS_DISCR_SHELL_X_TAG_BOTTOM 1
#define SLABS_DISCR_SHELL_X_TAG_BACK 4
#define SLABS_DISCR_SHELL_X_TAG_FRONT 5

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
                     double *_sc_restrict Z, void *data)
{
  mangll_locidx_t     il;
  double              x, y, R, q;

  const double        R1 = SL_SHELL_RADIUS_BOTTOM;
  const double        R2 = SL_SHELL_RADIUS_TOP;
  const double        R2byR1 = R2 / R1;
  const double        R1sqrbyR2 = R1 * R1 / R2;

  MANGLL_ASSERT (0 <= tag && tag < 24);

#ifdef MANGLL_DEBUG
  for (il = 0; il < np; ++il) {
    MANGLL_ASSERT (EX[il] < 1.0 + SC_1000_EPS && EX[il] > -1.0 - SC_1000_EPS);
    MANGLL_ASSERT (EY[il] < 1.0 + SC_1000_EPS && EY[il] > -1.0 - SC_1000_EPS);
    MANGLL_ASSERT (EZ[il] < 2.0 + SC_1000_EPS && EZ[il] > 1.0 - SC_1000_EPS);
  }
#endif

  switch (tag / 4) {
  case SLABS_DISCR_SHELL_X_TAG_TOP:  /* top (+z) */
    for (il = 0; il < np; ++il) {
      /* transform abc[0] and y in-place for nicer grading */
      x = tan (EX[il] * M_PI_4);
      y = tan (EY[il] * M_PI_4);

      /* compute transformation ingredients */
      R = R1sqrbyR2 * pow (R2byR1, EZ[il]);
      q = R / sqrt (x * x + y * y + 1.);

      /* assign correct coordinates based on patch id */
      X[il] = +q * y;
      Y[il] = -q * x;
      Z[il] = +q;
    }
    break;
  case SLABS_DISCR_SHELL_X_TAG_LEFT:  /* left (-y) */
    for (il = 0; il < np; ++il) {
      x = tan (EX[il] * M_PI_4);
      y = tan (EY[il] * M_PI_4);
      R = R1sqrbyR2 * pow (R2byR1, EZ[il]);
      q = R / sqrt (x * x + y * y + 1.);
      X[il] = -q;
      Y[il] = -q * x;
      Z[il] = +q * y;
    }
    break;
  case SLABS_DISCR_SHELL_X_TAG_BOTTOM:  /* bottom (-z) */
    for (il = 0; il < np; ++il) {
      x = tan (EX[il] * M_PI_4);
      y = tan (EY[il] * M_PI_4);
      R = R1sqrbyR2 * pow (R2byR1, EZ[il]);
      q = R / sqrt (x * x + y * y + 1.);
      X[il] = -q * y;
      Y[il] = -q * x;
      Z[il] = -q;
    }
    break;
  case SLABS_DISCR_SHELL_X_TAG_RIGHT:  /* right (+y) */
    for (il = 0; il < np; ++il) {
      x = tan (EX[il] * M_PI_4);
      y = tan (EY[il] * M_PI_4);
      R = R1sqrbyR2 * pow (R2byR1, EZ[il]);
      q = R / sqrt (x * x + y * y + 1.);
      X[il] = +q;
      Y[il] = -q * x;
      Z[il] = -q * y;
    }
    break;
  case SLABS_DISCR_SHELL_X_TAG_BACK:  /* back (-x) */
    for (il = 0; il < np; ++il) {
      x = tan (EX[il] * M_PI_4);
      y = tan (EY[il] * M_PI_4);
      R = R1sqrbyR2 * pow (R2byR1, EZ[il]);
      q = R / sqrt (x * x + y * y + 1.);
      X[il] = -q * x;
      Y[il] = +q;
      Z[il] = +q * y;
    }
    break;
  case SLABS_DISCR_SHELL_X_TAG_FRONT:  /* front (+x) */
    for (il = 0; il < np; ++il) {
      x = tan (EX[il] * M_PI_4);
      y = tan (EY[il] * M_PI_4);
      R = R1sqrbyR2 * pow (R2byR1, EZ[il]);
      q = R / sqrt (x * x + y * y + 1.);
      X[il] = +q * x;
      Y[il] = -q;
      Z[il] = +q * y;
    }
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

/**
 * Gets tag corresponding to a part of the spherical shell domain from element
 * coordinates:
 *   If (EX,EY) \in [-1,+1) x [-1,+1), then return `top` tag.
 *   If (EX,EY) \in [-3,-1) x [-1,+1), then return `back` tag.
 *   If (EX,EY) \in [+1,+3) x [-1,+1), then return `front` tag.
 *   If (EX,EY) \in [-1,+1) x [-3,-1), then return `left` tag.
 *   If (EX,EY) \in [-1,+1) x [+1,+3), then return `right` tag.
 */
static mangll_tag_t
slabs_discr_shell_chunk_and_slice_get_tag (const double EX, const double EY)
{
  const mangll_tag_t  n_tags_per_part = 4;

  if ( (-1.0 <= EX && EX < 1.0) && (-1.0 <= EY && EY < 1.0) ) {
    /* return `top` if (EX,EY) \in [-1,+1] x [-1,+1] */
    return SLABS_DISCR_SHELL_X_TAG_TOP * n_tags_per_part;
  }
  else if ( (-3.0 <= EX && EX < -1.0) && (-1.0 <= EY && EY < 1.0) ) {
    /* return `back` if (EX,EY) \in [-3,-1) x [-1,+1] */
    return SLABS_DISCR_SHELL_X_TAG_BACK * n_tags_per_part;
  }
  else if ( (1.0 <= EX && EX < 3.0) && (-1.0 <= EY && EY < 1.0) ) {
    /* return `front` if (EX,EY) \in (+1,+3] x [-1,+1] */
    return SLABS_DISCR_SHELL_X_TAG_FRONT * n_tags_per_part;
  }
  else if ( (-1.0 <= EX && EX < 1.0) && (-3.0 <= EY && EY < -1.0) ) {
    /* return `left` if (EX,EY) \in [-1,+1] x [-3,-1) */
    return SLABS_DISCR_SHELL_X_TAG_LEFT * n_tags_per_part;
  }
  else if ( (-1.0 <= EX && EX < 1.0) && (1.0 <= EY && EY < 3.0) ) {
    /* return `right` if (EX,EY) \in [-1,+1] x (+1,+3] */
    return SLABS_DISCR_SHELL_X_TAG_RIGHT * n_tags_per_part;
  }
  else { /* unknown part of spherical shell domain */
    YMIR_ABORT_NOT_REACHED ();
  }
}

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
                                     double *_sc_restrict Z, void *data)
{
  const mangll_tag_t  n_tags_per_part = 4;
  mangll_tag_t        shell_tag;
  double             *E_replace;
  mangll_locidx_t     il;

  /* set tag for shell X-function */
  shell_tag = slabs_discr_shell_chunk_and_slice_get_tag (EX[0], EY[0]);

  /* init replacement for a component of element coordinates */
  if (shell_tag/n_tags_per_part != SLABS_DISCR_SHELL_X_TAG_TOP) {
    E_replace = YMIR_ALLOC (double, np);
  }

  /* call shell X-function */
  switch (shell_tag/n_tags_per_part) {
  case SLABS_DISCR_SHELL_X_TAG_TOP:  /* top */
    slabs_discr_shell_X (shell_tag, np, EX, EY, EZ, X, Y, Z, data);
    break;
  case SLABS_DISCR_SHELL_X_TAG_BACK:  /* back */
    for (il = 0; il < np; ++il) { /* shift EX from [-3,-1] to [-1,+1] */
      E_replace[il] = EX[il] + 2.0;
    }
    slabs_discr_shell_X (shell_tag, np, E_replace, EY, EZ, X, Y, Z, data);
    break;
  case SLABS_DISCR_SHELL_X_TAG_FRONT:  /* front */
    for (il = 0; il < np; ++il) { /* shift EX from [+1,+3] to [-1,+1] */
      E_replace[il] = EX[il] - 2.0;
    }
    slabs_discr_shell_X (shell_tag, np, E_replace, EY, EZ, X, Y, Z, data);
    break;
  case SLABS_DISCR_SHELL_X_TAG_LEFT:  /* left */
    for (il = 0; il < np; ++il) { /* shift EY from [-3,-1] to [-1,+1] */
      E_replace[il] = EY[il] + 2.0;
    }
    slabs_discr_shell_X (shell_tag, np, EX, E_replace, EZ, X, Y, Z, data);
    break;
  case SLABS_DISCR_SHELL_X_TAG_RIGHT:  /* right */
    for (il = 0; il < np; ++il) { /* shift EY from [+1,+3] to [-1,+1] */
      E_replace[il] = EY[il] - 2.0;
    }
    slabs_discr_shell_X (shell_tag, np, EX, E_replace, EZ, X, Y, Z, data);
    break;
  default: /* unknown part of spherical shell domain */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* destroy */
  if (shell_tag/n_tags_per_part != SLABS_DISCR_SHELL_X_TAG_TOP) {
    YMIR_FREE (E_replace);
  }
}

/**
 * Creates new p4est object.
 */
p8est_t *
slabs_discr_p8est_new (MPI_Comm mpicomm,
                       slabs_physics_options_t *physics_options,
                       slabs_discr_options_t *discr_options)
{
  const int           minlevel = (int) discr_options->minlevel;
  const char         *refine = discr_options->refine;
  p8est_refine_t      refine_fn = NULL;
  int                 level_ext;  /* passed to p8est_new_ext */
  int                 level_ref;  /* used inside first refinement */
  int                 level_max;  /* maximum level for p8est_refine */

  p8est_connectivity_t *conn;
  p8est_t            *p8est;

  /*
   * create connectivity
   */

  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
    conn = p8est_connectivity_new_unitcube ();
    break;

  case SL_DOMAIN_BRICK:
    conn = p8est_connectivity_new_brick (physics_options->domain_brick_dx,
                                         physics_options->domain_brick_dy,
                                         physics_options->domain_brick_dz,
                                         0, 0, 0);
    /* scale coordinates of vertices, s.t. the domain has unit height:
     *   x component: [0, dx] -> 1/dz * [0, dx]
     *   y component: [0, dy] -> 1/dz * [0, dy]
     *   z component: [0, dz] -> 1/dz * [0, dz]
     */
    {
      p4est_topidx_t  vi;
      double         *vertices = conn->vertices;
      double          scale = 1.0 / ((double) physics_options->domain_brick_dz);

      for (vi = 0; vi < conn->num_vertices; vi++) { /* loop over all vertices */
        vertices[3 * vi    ] *= scale;
        vertices[3 * vi + 1] *= scale;
        vertices[3 * vi + 2] *= scale;
      }
    }
    break;

  case SL_DOMAIN_SHELL:
    conn = p8est_connectivity_new_shell ();
    break;

  case SL_DOMAIN_SHELL_CHUNK:
    conn = p8est_connectivity_new_unitcube ();
    /* shift coordinates of vertices, so (modified) shell_X fnc can be used:
     *   x component: [0, 1] -> [-0.5, 0.5]
     *   y component: [0, 1] -> [-0.5, 0.5]
     *   z component: [0, 1] -> [ 1  , 2  ]
     */
    {
      p4est_topidx_t  vi;
      double         *vertices = conn->vertices;

      for (vi = 0; vi < conn->num_vertices; vi++) { /* loop over all vertices */
        vertices[3 * vi    ] -= 0.5;
        vertices[3 * vi + 1] -= 0.5;
        vertices[3 * vi + 2] += 1.0;
      }
    }
    break;

  case SL_DOMAIN_SHELL_SLICE:
    conn = p8est_connectivity_new_brick (physics_options->domain_brick_dx,
                                         physics_options->domain_brick_dy,
                                         physics_options->domain_brick_dz,
                                         0, 0, 0);
    /* shift coordinates of vertices, so (modified) shell_X fnc can be used:
     *   x component: [0, dx] -> 1/dz * [-dx/2, dx/2]
     *   y component: [0, dy] -> 1/dz * [-dy/2, dy/2]
     *   z component: [0, dz] ->        [ 1   , 2   ]
     */
    {
      p4est_topidx_t  vi;
      double         *vertices = conn->vertices;
      double          scale = 1.0 / ((double) physics_options->domain_brick_dz);

      for (vi = 0; vi < conn->num_vertices; vi++) { /* loop over all vertices */
        vertices[3 * vi] =
          scale * (  vertices[3 * vi]
                   - 0.5 * ((double) physics_options->domain_brick_dx) );
        vertices[3 * vi + 1] =
          scale * (  vertices[3 * vi + 1]
                   - 0.5 * ((double) physics_options->domain_brick_dy) );
        vertices[3 * vi + 2] = scale * vertices[3 * vi + 2] + 1.0;
      }
    }
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /*
   * set up initial refinement
   */

  level_ext = level_ref = level_max = minlevel;
  if (refine != NULL) {
    if (!strcmp (refine, "uniform")) {
      /* uniform refinement is the default action */
    }
    else if (!strcmp (refine, "evenodd")) {
      refine_fn = mangll_p8est_refine_evenodd;
      level_ext = minlevel - 1;
    }
    else if (!strcmp (refine, "oddeven")) {
      refine_fn = mangll_p8est_refine_oddeven;
      level_ext = minlevel - 1;
    }
    else if (!strcmp (refine, "origin")) {
      refine_fn = mangll_p8est_refine_origin;
      level_ext = 0;
    }
    else if (!strcmp (refine, "fractal")) {
      refine_fn = mangll_p8est_refine_fractal;
      level_ext = minlevel - 4;
    }
    else if (!strcmp (refine, "halfhalf")) {
      refine_fn = mangll_p8est_refine_halfhalf;
      level_ext = minlevel - 1;
    }
    else if (!strcmp (refine, "finercenter2")) {
      refine_fn = mangll_p8est_refine_finercenter2;
      level_ext = level_ref = minlevel - 1;
    }
    else if (!strcmp (refine, "finercenter4")) {
      refine_fn = mangll_p8est_refine_finercenter4;
      level_ext = level_ref = minlevel - 1;
    }
    else {
      YMIR_GLOBAL_LERRORF ("Invalid refinement type `%s`\n", refine);
      return NULL;
    }
  }
  if (discr_options->maxlevel >= 0) {
    level_ext = SC_MIN (discr_options->maxlevel, level_ext);
    level_max = SC_MIN (discr_options->maxlevel, level_max);
  }

  /*
   * create p8est mesh
   */

  /* create p8est structure */
  p8est = p8est_new_ext (mpicomm, conn, 0, level_ext, 1, 0, NULL, &level_ref);

  /* refine p8est */
  if (refine != NULL && refine_fn != NULL) {
    p8est_refine_ext (p8est, 1, level_max, refine_fn, NULL, NULL);
    p8est_partition_ext (p8est, 1, NULL);
    p8est_balance (p8est, P8EST_CONNECT_FULL, NULL);
    p8est_partition_ext (p8est, 1, NULL);
  }

  /* initialize user data */
  slabs_discr_p8est_init_data (p8est, discr_options->inspect_p4est);

  /* return p8est mesh */
  return p8est;
}

/**
 * Destroys p4est object.
 */
void
slabs_discr_p8est_destroy (p8est_t *p8est)
{
  p8est_connectivity_t *conn = p8est->connectivity;

  slabs_discr_p8est_clear_data (p8est);
  p8est_destroy (p8est);
  p8est_connectivity_destroy (conn);
}

/**
 * Initializes user data of a p4est object.
 */
void
slabs_discr_p8est_init_data (p8est_t *p8est, const int inspect_p4est)
{
  /* set p4est user pointer */
  p8est->user_pointer = NULL;

  /* set p4est user data for each quadrant */
  p8est_reset_data (p8est, sizeof (slabs_discr_p4est_qdata_t), NULL,
                    p8est->user_pointer);

  /* collect p4est statistics */
  //TODO delete deprecated code below; remove fnc param inspect_p4est
//if (inspect_p4est) {
//  p8est->inspect = YMIR_ALLOC_ZERO (p4est_inspect_t, 1);
//  p8est->inspect->owned = 0;

//  p8est->inspect->stats = sc_statistics_new (p8est->mpicomm);
//  sc_statistics_add_empty (p8est->inspect->stats, "time_balance");
//  sc_statistics_add_empty (p8est->inspect->stats, "time_partition_ext");
//  sc_statistics_add_empty (p8est->inspect->stats, "time_partition_given");
//  sc_statistics_add_empty (p8est->inspect->stats,
//                           "time_partition_for_coarsening");
//}
}

/**
 * Clears user data from a p4est object.
 */
void
slabs_discr_p8est_clear_data (p8est_t *p8est)
{
  //TODO delete deprecated code below
#if 0
  if (p8est->inspect != NULL && p8est->inspect->stats != NULL) {
    /* gather & print p4est statistics */
    sc_statistics_compute (p8est->inspect->stats);
    YMIR_GLOBAL_INFO ("===================================================\n");
    sc_statistics_print (p8est->inspect->stats, -1, SC_LP_STATISTICS, 1, 1);
    YMIR_GLOBAL_INFO ("===================================================\n");

    /* destroy p4est statistics */
    sc_statistics_destroy (p8est->inspect->stats);
    YMIR_FREE (p8est->inspect);
    p8est->inspect = NULL;
  }
#ifdef YMIR_DEBUG
  else {
    YMIR_ASSERT (p8est->inspect == NULL);
  }
#endif
#endif
}

/**
 *
 */
void
slabs_discr_options_set_boundary (slabs_discr_options_t *discr_options,
                                  p8est_t *p8est,
                                  slabs_physics_options_t *physics_options)
{
  /* set element to boundary function */
  discr_options->e_to_fm_fn = ymir_mesh_e_to_fm_tag;

  /* set boundary faces */
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_BRICK:
  case SL_DOMAIN_SHELL_CHUNK:
  case SL_DOMAIN_SHELL_SLICE:
    discr_options->tree_to_bf = slabs_cube_tree_to_bf (p8est->connectivity);
    break;

  case SL_DOMAIN_SHELL:
    discr_options->tree_to_bf = slabs_shell_tree_to_bf (p8est->connectivity);
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }
}

/**
 *
 */
void
slabs_discr_options_clear_boundary (slabs_discr_options_t *discr_options)
{
  /* destroy boundary faces */
  YMIR_FREE (discr_options->tree_to_bf);
}

/**
 * Creates new mangll mesh structure.
 */
mangll_mesh_t *
slabs_discr_mangll_mesh_new (p8est_t *p8est,
                             slabs_discr_options_t *discr_options)
{
  p8est_ghost_t      *ghost;
  mangll_mesh_t      *mangll_mesh;

  /* create p8est ghost */
  ghost = p8est_ghost_new (p8est, P8EST_CONNECT_FULL);

  /* create mangll mesh structure */
  mangll_mesh = mangll_p8est_mesh_new_full (p8est, ghost);

  /* assign mapping to physical space */
  mangll_mesh->X_fn = discr_options->X_fn;
  mangll_mesh->X_data = NULL;

  /* destroy ghost */
  p8est_ghost_destroy (ghost);

  /* return mangll mesh */
  return mangll_mesh;
}

/**
 * Creates new mangll and cnodes objects.
 */
void
slabs_discr_mangll_and_cnodes_new (mangll_t **mangll,
                                   mangll_cnodes_t **cnodes,
                                   p8est_t *p8est,
                                   slabs_discr_options_t *discr_options)
{
  MPI_Comm            mpicomm = p8est->mpicomm;
  const int           order = discr_options->order;
  const mangll_refel_quadrature_type_t  quad_type = MANGLL_REFEL_QUAD_GAUSS;

  p8est_ghost_t      *ghost;
  mangll_mesh_t      *mangll_mesh;

  /* create p8est ghost */
  if (cnodes != NULL) {
    ghost = p8est_ghost_new (p8est, P8EST_CONNECT_FULL);
  }
  else {
    ghost = p8est_ghost_new (p8est, P8EST_CONNECT_FACE);
  }

  /* create mangll mesh structure */
  mangll_mesh = mangll_p8est_mesh_new_full (p8est, ghost);

  /* assign mapping to physical space */
  mangll_mesh->X_fn = discr_options->X_fn;
  mangll_mesh->X_data = NULL;

  YMIR_ASSERT (mangll != NULL);
  if (cnodes != NULL) { /* if mesh with continuous nodes should be created */
    /* create continuous node structure */
    *cnodes = mangll_p8est_cnodes_new (p8est, ghost, discr_options->order);

    /* create mangll structure with continuous geometry */
    *mangll = mangll_new_ext (mpicomm, order, quad_type, mangll_mesh, *cnodes);
  }
  else {
    /* create mangll structure with discontinuous geometry */
    *mangll = mangll_new_ext (mpicomm, order, quad_type, mangll_mesh, NULL);
  }

  /* destroy ghost */
  p8est_ghost_destroy (ghost);
}

/**
 * Creates new mangll object for interpolation and partitioning during AMR.
 */
static void
slabs_discr_mangll_for_interp_new (mangll_t **mangll,
                                   p8est_t *p8est,
                                   slabs_discr_options_t *discr_options)
{
  MPI_Comm            mpicomm = p8est->mpicomm;
  const int           order = discr_options->order;
  const mangll_refel_quadrature_type_t  quad_type = MANGLL_REFEL_QUAD_GAUSS;

  mangll_mesh_t      *mangll_mesh;

  /* create mangll mesh structure */
  mangll_mesh = mangll_p8est_mesh_new_partial (p8est);

  /* create mangll structure with discontinuous geometry */
  *mangll = mangll_for_interpolation_new (mpicomm, order, quad_type,
                                          mangll_mesh);
}

/**
 * Destroys mangll object for interpolation and partitioning during AMR.
 */
static void
slabs_discr_mangll_for_interp_destroy (mangll_t *mangll)
{
  mangll_for_interpolation_destroy (mangll);
}

/**
 * Creates new ymir mesh and pressure element structures.
 */
void
slabs_discr_ymir_new (ymir_mesh_t **mesh, ymir_pressure_elem_t **press_elem,
                      mangll_t *mangll, mangll_cnodes_t *cnodes,
                      slabs_discr_options_t *discr_options)
{
  /* create ymir mesh */
  *mesh = ymir_mesh_new_ext (mangll, cnodes, discr_options->e_to_fm_fn,
                             discr_options->tree_to_bf,
                             1 /* skip face prealloc */,
                             ymir_stress_pc_gmg /* skip diag prealloc */);

  /* create pressure element */
  if (press_elem != NULL) {
    *press_elem = ymir_pressure_elem_new (mangll->refel, mangll->ompsize);
  }
}

/**
 * Destroys ymir mesh and corresponding mangll structures.
 */
void
slabs_discr_ymir_mangll_destroy (ymir_mesh_t *mesh,
                                 ymir_pressure_elem_t *press_elem)
{
  /* destroy pressure element */
  if (press_elem != NULL) {
    ymir_pressure_elem_destroy (press_elem);
  }

  /* destroy mangll, cnodes, and ymir_mesh */
  if (mesh != NULL) {
    mangll_destroy (mesh->ma);
    mangll_p8est_cnodes_destroy (mesh->cnodes);
    ymir_mesh_destroy (mesh);
  }
}

/**
 * Callback function for p8est refinement. Increases the refinement radially.
 */
static int
slabs_cube_refine_z_direction (p8est_t *p8est, p4est_topidx_t tree,
                               p8est_quadrant_t *quadrant)
{
  slabs_discr_options_t  *discr_options =
                            (slabs_discr_options_t *) p8est->user_pointer;
  const int           minlevel = discr_options->minlevel;
  const int           maxlevel = discr_options->maxlevel;
  double             *radius = discr_options->refine_radius;
  const int           n_radii = discr_options->refine_n_radii;
  const int           l = quadrant->level;
  double              coord[P8EST_DIM];
  double              r;
  int                 j;

  /* refine if level is below minimum level */
  if (l < minlevel) {
    return 1;
  }

  /* get z-coordinate of quadrant */
  p8est_qcoord_to_vertex (p8est->connectivity, tree, quadrant->x,
                          quadrant->y, quadrant->z, coord);
  r = coord[P8EST_DIM - 1];

  /* refine more the further away we are from the shell's center */
  for (j = 1; j <= n_radii; j++) {
    if (   (r > radius[j-1]) && (l < minlevel + j)
        && (maxlevel == 0 || (maxlevel > 0 && l < maxlevel)) ) {
      /* refine quadrant */
      return 1;
    }
  }

  /* do not refine quadrant */
  return 0;
}

/**
 * Callback function for p8est refinement. Increases the refinement radially.
 */
static int
slabs_shell_refine_radially (p8est_t *p8est, p4est_topidx_t tree,
                             p8est_quadrant_t *quadrant)
{
  slabs_discr_options_t  *discr_options =
                            (slabs_discr_options_t *) p8est->user_pointer;
  const int           minlevel = discr_options->minlevel;
  const int           maxlevel = discr_options->maxlevel;
  double             *radius = discr_options->refine_radius;
  const int           n_radii = discr_options->refine_n_radii;
  const int           l = quadrant->level;
  double              coord[P8EST_DIM];
  double              r;
  int                 j;

  /* refine if level is below minimum level */
  if (l < minlevel) {
    return 1;
  }

  /* get z-coordinate of quadrant; shift it from interval [1,2] to [0,1] */
  p8est_qcoord_to_vertex (p8est->connectivity, tree, quadrant->x,
                          quadrant->y, quadrant->z, coord);
  r = coord[P8EST_DIM - 1] - 1.0;

  /* refine more the further away we are from the shell's center */
  for (j = 1; j <= n_radii; j++) {
    if (   (r > radius[j-1]) && (l < minlevel + j)
        && (maxlevel == 0 || (maxlevel > 0 && l < maxlevel)) ) {
      /* refine quadrant */
      return 1;
    }
  }

  /* do not refine quadrant */
  return 0;
}

/**
 * Refines mesh radially (for spherical geometry) or in z-direction (for cubes).
 */
void
slabs_discr_refine_wrt_depth (p8est_t *p8est,
                              slabs_physics_options_t *physics_options,
                              slabs_discr_options_t *discr_options)
{
  void               *forest_user_pointer;

  /* check if there is anything to refine */
  if (discr_options->refine_n_radii <= 0) {
    return;
  }

  /* swap user pointer */
  forest_user_pointer = p8est->user_pointer;
  p8est->user_pointer = discr_options;

  /* call refine function */
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_BRICK:
    p8est_refine (p8est, 1, slabs_cube_refine_z_direction, NULL);
    break;

  case SL_DOMAIN_SHELL:
  case SL_DOMAIN_SHELL_CHUNK:
  case SL_DOMAIN_SHELL_SLICE:
    p8est_refine (p8est, 1, slabs_shell_refine_radially, NULL);
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* partition p8est and restore old user pointer */
  p8est_partition_ext (p8est, 1, NULL);
  p8est->user_pointer = forest_user_pointer;
}

/**
 * Refines within a certain distance about a layer, which is a surface at a
 * fixed radius.
 */
static inline int
slabs_discr_refine_at_layer_elem (p8est_quadrant_t *quadrant,
                                  const double radius,
                                  const double layer_radius,
                                  const double layer_maxdist,
                                  const int8_t layer_maxlevel,
                                  const int n_trees,
                                  const int8_t maxlevel)
{
  const int8_t        l = quadrant->level;
  const double        dist_above = radius - layer_radius;
  const double        dist_below = -dist_above;
  const double        maxdist_above = layer_maxdist;
  const double        maxdist_below = layer_maxdist +
                        1.0 / ((double) (1 << l)) / ((double) n_trees);

  /* refine if within a certain distance and if max level not reached */
  if ( ( (0 <= dist_above && dist_above <= maxdist_above) ||
         (0 <= dist_below && dist_below <= maxdist_below) )
       &&
       l < SC_MIN (maxlevel, layer_maxlevel) ) {
    /* refine quadrant */
    return 1;
  }
  else {
    /* do not refine quadrant */
    return 0;
  }
}

/**
 * Callback function for p8est refinement. Increases the refinement about
 * a layer that is located at a certain radius.
 */
static int
slabs_discr_refine_at_layer_cube (p8est_t *p8est, p4est_topidx_t tree,
                                  p8est_quadrant_t *quadrant)
{
  slabs_discr_options_t  *discr_options =
                            (slabs_discr_options_t *) p8est->user_pointer;
  const double        layer_radius = discr_options->refine_layer_radius;
  const double        layer_maxdist = discr_options->refine_layer_maxdist;
  const int8_t        layer_maxlevel = discr_options->refine_layer_maxlevel;
  const int           n_trees = discr_options->refine_layer_n_trees;
  const int8_t        maxlevel = discr_options->maxlevel;
  double              coord[P8EST_DIM];
  double              radius;

  /* get z-coordinate of quadrant */
  p8est_qcoord_to_vertex (p8est->connectivity, tree, quadrant->x,
                          quadrant->y, quadrant->z, coord);
  radius = coord[P8EST_DIM - 1];

  /* refine if within a certain distance */
  return slabs_discr_refine_at_layer_elem (quadrant, radius, layer_radius,
                                           layer_maxdist, layer_maxlevel,
                                           n_trees, maxlevel);
}

/**
 * Callback function for p8est refinement. Increases the refinement about
 * a layer that is located at a certain radius.
 */
static int
slabs_discr_refine_at_layer_shell (p8est_t *p8est, p4est_topidx_t tree,
                                   p8est_quadrant_t *quadrant)
{
  slabs_discr_options_t  *discr_options =
                            (slabs_discr_options_t *) p8est->user_pointer;
  const double        layer_radius = discr_options->refine_layer_radius;
  const double        layer_maxdist = discr_options->refine_layer_maxdist;
  const int8_t        layer_maxlevel = discr_options->refine_layer_maxlevel;
  const int           n_trees = discr_options->refine_layer_n_trees;
  const int8_t        maxlevel = discr_options->maxlevel;
  double              coord[P8EST_DIM];
  double              radius;

  /* get z-coordinate of quadrant; shift it from interval [1,2] to [0,1] */
  p8est_qcoord_to_vertex (p8est->connectivity, tree, quadrant->x,
                          quadrant->y, quadrant->z, coord);
  radius = coord[P8EST_DIM - 1] - 1.0;

  /* refine if within a certain distance */
  return slabs_discr_refine_at_layer_elem (quadrant, radius, layer_radius,
                                           layer_maxdist, layer_maxlevel,
                                           n_trees, maxlevel);
}

/**
 * Get number of trees in z- or radial direction;
 */
static int
slabs_discr_refine_at_layer_get_n_trees (slabs_physics_options_t
                                           *physics_options)
{
  int                 n_trees;

  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_SHELL:
  case SL_DOMAIN_SHELL_CHUNK:
    n_trees = 1;
    break;
  case SL_DOMAIN_BRICK:
  case SL_DOMAIN_SHELL_SLICE:
    n_trees = physics_options->domain_brick_dz;
    break;
  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }

  return n_trees;
}

/**
 * Refines mesh about a layer that is located at a certain depth.
 */
void
slabs_discr_refine_at_layer (p8est_t *p8est,
                             slabs_physics_options_t *physics_options,
                             slabs_discr_options_t *discr_options)
{
  void               *forest_user_pointer;

  /* quit if there is nothing to refine */
  if ( discr_options->refine_layer_radius <= 0.0 ||
       discr_options->refine_layer_maxdist <= 0.0 ) {
    return;
  }

  /* swap user pointer */
  forest_user_pointer = p8est->user_pointer;
  p8est->user_pointer = discr_options;

  /* set number of trees */
  discr_options->refine_layer_n_trees =
    slabs_discr_refine_at_layer_get_n_trees (physics_options);

  /* call refine function */
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_BRICK:
    p8est_refine (p8est, 1, slabs_discr_refine_at_layer_cube, NULL);
    break;
  case SL_DOMAIN_SHELL:
  case SL_DOMAIN_SHELL_CHUNK:
  case SL_DOMAIN_SHELL_SLICE:
    p8est_refine (p8est, 1, slabs_discr_refine_at_layer_shell, NULL);
    break;
  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* balance p8est */
  p8est_balance (p8est, P8EST_CONNECT_FULL, NULL);

  /* partition p8est and restore old user pointer */
  p8est_partition_ext (p8est, 1, NULL);
  p8est->user_pointer = forest_user_pointer;
}

/**
 *
 */
slabs_discr_enforce_refinement_data_t *
slabs_discr_enforce_refinement_data_new (
                                      slabs_physics_options_t *physics_options,
                                      slabs_discr_options_t *discr_options)
{
  slabs_discr_enforce_refinement_data_t *d;

  d = YMIR_ALLOC (slabs_discr_enforce_refinement_data_t, 1);
  d->physics_options = physics_options;
  d->discr_options = discr_options;

  return d;
}

/**
 *
 */
void
slabs_discr_enforce_refinement_data_destroy (
                                  slabs_discr_enforce_refinement_data_t *data)
{
  if (data != NULL) {
    YMIR_FREE (data);
  }
}

/**
 *
 */
void
slabs_discr_enforce_refinement_at_layer (p8est_t *p8est, void *data)
{
  slabs_discr_enforce_refinement_data_t *d =
    (slabs_discr_enforce_refinement_data_t *) data;
  slabs_physics_options_t *physics_options = d->physics_options;
  slabs_discr_options_t   *discr_options = d->discr_options;
  void               *forest_user_pointer;

  /* quit if there is nothing to refine */
  if ( !discr_options->enforce_refinement_at_layer ||
       discr_options->refine_layer_radius <= 0.0 ||
       discr_options->refine_layer_maxdist <= 0.0 ) {
    return;
  }

  /* swap user pointer */
  forest_user_pointer = p8est->user_pointer;
  p8est->user_pointer = discr_options;

  /* set number of trees */
  discr_options->refine_layer_n_trees =
    slabs_discr_refine_at_layer_get_n_trees (physics_options);

  /* call refine function */
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_BRICK:
    p8est_refine (p8est, 1, slabs_discr_refine_at_layer_cube, NULL);
    break;
  case SL_DOMAIN_SHELL:
  case SL_DOMAIN_SHELL_CHUNK:
  case SL_DOMAIN_SHELL_SLICE:
    p8est_refine (p8est, 1, slabs_discr_refine_at_layer_shell, NULL);
    break;
  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* restore old user pointer */
  p8est->user_pointer = forest_user_pointer;
}

/**
 * Callback function for p8est refinement. Increases the refinement at the
 * surface.
 */
static int
slabs_discr_refine_at_surface_cube (p8est_t *p8est, p4est_topidx_t tree,
                                    p8est_quadrant_t *quadrant)
{
  slabs_discr_options_t  *discr_options =
                            (slabs_discr_options_t *) p8est->user_pointer;
  const double        layer_radius = 1.0;
  const double        layer_maxdist = discr_options->refine_surface_maxdist;
  const int8_t        layer_maxlevel = discr_options->refine_surface_maxlevel;
  const int           n_trees = discr_options->refine_surface_n_trees;
  const int8_t        maxlevel = discr_options->maxlevel;
  double              coord[P8EST_DIM];
  double              radius;

  /* get z-coordinate of quadrant */
  p8est_qcoord_to_vertex (p8est->connectivity, tree, quadrant->x,
                          quadrant->y, quadrant->z, coord);
  radius = coord[P8EST_DIM - 1];

  /* refine if within a certain distance */
  return slabs_discr_refine_at_layer_elem (quadrant, radius, layer_radius,
                                           layer_maxdist, layer_maxlevel,
                                           n_trees, maxlevel);
}

/**
 * Callback function for p8est refinement. Increases the refinement at the
 * surface.
 */
static int
slabs_discr_refine_at_surface_shell (p8est_t *p8est, p4est_topidx_t tree,
                                     p8est_quadrant_t *quadrant)
{
  slabs_discr_options_t  *discr_options =
                            (slabs_discr_options_t *) p8est->user_pointer;
  const double        layer_radius = 1.0;
  const double        layer_maxdist = discr_options->refine_surface_maxdist;
  const int8_t        layer_maxlevel = discr_options->refine_surface_maxlevel;
  const int           n_trees = discr_options->refine_surface_n_trees;
  const int8_t        maxlevel = discr_options->maxlevel;
  double              coord[P8EST_DIM];
  double              radius;

  /* get z-coordinate of quadrant; shift it from interval [1,2] to [0,1] */
  p8est_qcoord_to_vertex (p8est->connectivity, tree, quadrant->x,
                          quadrant->y, quadrant->z, coord);
  radius = coord[P8EST_DIM - 1] - 1.0;

  /* refine if within a certain distance */
  return slabs_discr_refine_at_layer_elem (quadrant, radius, layer_radius,
                                           layer_maxdist, layer_maxlevel,
                                           n_trees, maxlevel);
}

/**
 * Refines mesh at the surface.
 */
void
slabs_discr_refine_at_surface (p8est_t *p8est,
                               slabs_physics_options_t *physics_options,
                               slabs_discr_options_t *discr_options)
{
  void               *forest_user_pointer;

  /* quit if there is nothing to refine */
  if (discr_options->refine_surface_maxdist <= 0.0) {
    return;
  }

  /* swap user pointer */
  forest_user_pointer = p8est->user_pointer;
  p8est->user_pointer = discr_options;

  /* set number of trees */
  discr_options->refine_surface_n_trees =
    slabs_discr_refine_at_layer_get_n_trees (physics_options);

  /* call refine function */
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_BRICK:
    p8est_refine (p8est, 1, slabs_discr_refine_at_surface_cube, NULL);
    break;
  case SL_DOMAIN_SHELL:
  case SL_DOMAIN_SHELL_CHUNK:
  case SL_DOMAIN_SHELL_SLICE:
    p8est_refine (p8est, 1, slabs_discr_refine_at_surface_shell, NULL);
    break;
  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* balance p8est */
  p8est_balance (p8est, P8EST_CONNECT_FULL, NULL);

  /* partition p8est and restore old user pointer */
  p8est_partition_ext (p8est, 1, NULL);
  p8est->user_pointer = forest_user_pointer;
}

/**
 * Creates a new AMR indicator.
 */
slabs_discr_amr_indicator_t *
slabs_discr_amr_indicator_new (slabs_discr_amr_indicator_type_t type,
                               const double tol_min, const double tol_max,
                               const int8_t level_min, const int8_t level_max,
                               const p4est_locidx_t n_quadrants_loc)
{
  slabs_discr_amr_indicator_t *indicator;

  /* create new indicator */
  indicator = YMIR_ALLOC (slabs_discr_amr_indicator_t, 1);

  /* initialize indicator */
  indicator->val = sc_dmatrix_new (1, n_quadrants_loc);;
  indicator->type = type;
  indicator->tol_min = tol_min;
  indicator->tol_max = tol_max;
  indicator->level_min = level_min;
  indicator->level_max = level_max;

  return indicator;
}

/**
 * Destroys an AMR indicator.
 */
void
slabs_discr_amr_indicator_destroy (slabs_discr_amr_indicator_t *indicator)
{
  /* destroy values of indicator */
  if (indicator->val != NULL) {
    sc_dmatrix_destroy (indicator->val);
  }

  /* destroy indicator */
  YMIR_FREE (indicator);
}

/**
 * Callback function for sorting indicator values.
 */
int
slabs_discr_amr_indicator_compare_fn (const void *val1, const void *val2)
{
  const double       *v1 = (const double *) val1;
  const double       *v2 = (const double *) val2;

  if (*v1 < *v2) {
    return -1;
  }
  else if (*v1 == *v2) {
    return  0;
  }
  else {
    return  1;
  }
}

/**
 * Reduce max tolerance of indicator in order to limit refinement.
 */
void
slabs_discr_amr_indicator_increase_tol_max (
                                   slabs_discr_amr_indicator_t *indicator,
                                   const p4est_locidx_t n_refine_quadrants_loc,
                                   const p4est_gloidx_t n_refine_quadrants_glo,
                                   const p4est_gloidx_t n_refine_quadrants_max)
{
  const double        n_refine_loc = (double) n_refine_quadrants_loc;
  const double        n_refine_glo = (double) n_refine_quadrants_glo;
  const double        n_refine_max = (double) n_refine_quadrants_max;
  double              n_refine_next;
  const size_t        indicator_size = indicator->val->n;
  double             *indicator_sort_data;
  sc_array_t          indicator_sort_array;
  size_t              idx;
#ifdef YMIR_DEBUG
  const double        tol_max = indicator->tol_max;
#endif

  /* check input */
  YMIR_ASSERT (sc_dmatrix_is_valid (indicator->val));
  YMIR_ASSERT (n_refine_quadrants_loc <= n_refine_quadrants_glo);
  YMIR_ASSERT (n_refine_quadrants_max <= n_refine_quadrants_glo);

  /* copy indicator values */
  indicator_sort_data = YMIR_ALLOC (double, indicator_size);
  memcpy (indicator_sort_data, indicator->val->e[0],
          (size_t) indicator_size * sizeof (double));

  /* sort indicator values in ascending order */
  sc_array_init_data (&indicator_sort_array, indicator_sort_data,
                      (size_t) sizeof (double), (size_t) indicator_size);
  sc_array_sort (&indicator_sort_array, slabs_discr_amr_indicator_compare_fn);
  sc_array_reset (&indicator_sort_array);

  /* set new #quadrants that are allowed to be marked for refinement */
  n_refine_next = floor (n_refine_loc / n_refine_glo * n_refine_max);
  YMIR_ASSERT (((size_t) n_refine_next) < indicator_size);

  /* increase max tol s.t. only allowd #quadrants are marked for refinement */
  idx = indicator_size - ((size_t) n_refine_next);
  if (0 < idx && idx < indicator_size) {
    const double        ind_val_low = indicator_sort_data[idx-1];
    const double        ind_val_high = indicator_sort_data[idx];

    YMIR_ASSERT (isfinite (indicator_sort_data[idx-1]));
    YMIR_ASSERT (isfinite (indicator_sort_data[idx]));

    if ((ind_val_high - ind_val_low) < SC_EPS) { /* if no separation */
      indicator->tol_max = 2.0 * ind_val_high;
    }
    else { /* if indicator values are well separated */
      indicator->tol_max = 0.5 * (ind_val_low + ind_val_high);
    }
  }
  else if (indicator_size < idx) { /* if no refinement allowed */
    YMIR_ASSERT (isfinite (indicator_sort_data[indicator_size-1]));
    indicator->tol_max = 2.0 * indicator_sort_data[indicator_size-1];
  }
  YMIR_ASSERT (tol_max <= indicator->tol_max);

  /* destroy sorted indicator values */
  YMIR_FREE (indicator_sort_data);
}

/**
 * Creates a new AMR marker and initializes to zero marked elements.
 */
slabs_discr_amr_marker_t *
slabs_discr_amr_marker_new (p4est_t *p4est)
{
  const p4est_locidx_t  n_flags = p4est->local_num_quadrants;
  slabs_discr_amr_marker_t *marker;

  /* create new marker */
  marker = YMIR_ALLOC (slabs_discr_amr_marker_t, 1);

  /* initialize marker */
  marker->flag = YMIR_ALLOC_ZERO (slabs_discr_amr_marker_flag_t, n_flags);
  marker->n_coarsen_quadrants_loc = 0;
  marker->n_refine_quadrants_loc = 0;
  marker->p4est = p4est;

  return marker;
}

/**
 * Creates a copy of an existing AMR marker.
 */
slabs_discr_amr_marker_t *
slabs_discr_amr_marker_duplicate (slabs_discr_amr_marker_t *marker)
{
  const p4est_locidx_t  n_flags = marker->p4est->local_num_quadrants;
  slabs_discr_amr_marker_t *marker_clone;

  /* create new marker */
  marker_clone = YMIR_ALLOC (slabs_discr_amr_marker_t, 1);

  /* initialize marker */
  marker_clone->flag = YMIR_ALLOC (slabs_discr_amr_marker_flag_t, n_flags);
  marker_clone->n_coarsen_quadrants_loc = marker->n_coarsen_quadrants_loc;
  marker_clone->n_refine_quadrants_loc = marker->n_refine_quadrants_loc;
  marker_clone->p4est = marker->p4est;

  /* copy marker flags */
  memcpy (marker_clone->flag, marker->flag,
          (size_t) n_flags * sizeof (p4est_locidx_t));

  return marker_clone;
}

/**
 * Destroys an AMR marker.
 */
void
slabs_discr_amr_marker_destroy (slabs_discr_amr_marker_t *marker)
{
  /* destroy flags of marker */
  if (marker->flag != NULL) {
    YMIR_FREE (marker->flag);
  }

  /* destroy marker */
  YMIR_FREE (marker);
}

/**
 * Sets or merges an AMR indicator into an AMR marker.
 */
static void
slabs_discr_amr_marker_set_merge_indicator (
                                        slabs_discr_amr_marker_t *marker,
                                        slabs_discr_amr_indicator_t *indicator,
                                        const int set)
{
  slabs_discr_amr_marker_flag_t *marker_flag = marker->flag;
  p4est_t            *p4est = marker->p4est;
  const double        tol_min = indicator->tol_min;
  const double        tol_max = indicator->tol_max;
  const int8_t        level_min = indicator->level_min;
  const int8_t        level_max = indicator->level_max;
  const double       *indicator_data = indicator->val->e[0];

  p4est_topidx_t      ti;
  size_t              tqi;

  /* check input */
  YMIR_ASSERT (0 <= level_min && level_min <= level_max);
  YMIR_ASSERT (0 < indicator->val->m);
  YMIR_ASSERT (((p4est_locidx_t) indicator->val->n) ==
               p4est->local_num_quadrants);

  /* reset quadrant counters if setting mode */
  if (set) {
    marker->n_coarsen_quadrants_loc = 0;
    marker->n_refine_quadrants_loc = 0;
  }

  for (ti = p4est->first_local_tree; ti <= p4est->last_local_tree; ++ti) {
    p4est_tree_t       *t = p4est_tree_array_index (p4est->trees, ti);
    sc_array_t         *tquadrants = &t->quadrants;
    const size_t        tqoffset = t->quadrants_offset;

    for (tqi = 0; tqi < tquadrants->elem_count; ++tqi) {
      p4est_quadrant_t   *q = p4est_quadrant_array_index (tquadrants, tqi);
      const p4est_locidx_t  qi = tqoffset + tqi;
      const double        val = indicator_data[qi];
      slabs_discr_amr_marker_flag_t  m;

      /*
       * suggest coarsening or refinement from indicator value and level
       */
      if (0 < level_max && level_max < q->level) {
        /* mark for coarsening if above max level */
        m = SLABS_DISCR_AMR_MARKER_COARSEN;
      }
      else if (q->level < level_min) {
        /* mark for refinement if below min level */
        m = SLABS_DISCR_AMR_MARKER_REFINE;
      }
      else if (0.0 <= val) {
        if (val < tol_min && level_min < q->level) {
          /* mark for coarsening if value is below min tolerance and
           * min level is not reached */
          m = SLABS_DISCR_AMR_MARKER_COARSEN;
          YMIR_ASSERT (0.0 < tol_min);
        }
        else if ( 0.0 < tol_max && tol_max < val &&
                  (level_max <= 0 || q->level < level_max) ) {
          /* mark for refinement if value is above max tolerance and
           * max level is not reached */
          m = SLABS_DISCR_AMR_MARKER_REFINE;
        }
        else {
          /* mark for nothing (no coarsening or refinement) */
          m = SLABS_DISCR_AMR_MARKER_NONE;
        }
      }
      else {
        /* mark for nothing (no coarsening or refinement) */
        m = SLABS_DISCR_AMR_MARKER_NONE;
      }

      /*
       * set marker for coarsening or refinement
       */
      if (set) { /* if setting mode */
        marker_flag[qi] = m;
        if (SLABS_DISCR_AMR_MARKER_COARSEN == m) {
          marker->n_coarsen_quadrants_loc++;
        }
        else if (SLABS_DISCR_AMR_MARKER_REFINE == m) {
          marker->n_refine_quadrants_loc++;
        }
#ifdef YMIR_DEBUG
        else { YMIR_ASSERT (SLABS_DISCR_AMR_MARKER_NONE == m); }
#endif
      }
      else { /* if merging mode */
        if (SLABS_DISCR_AMR_MARKER_REFINE == m) {
          if (SLABS_DISCR_AMR_MARKER_COARSEN == marker_flag[qi]) {
            marker->n_coarsen_quadrants_loc--;
            marker->n_refine_quadrants_loc++;
          }
          else if (SLABS_DISCR_AMR_MARKER_NONE == marker_flag[qi]) {
            marker->n_refine_quadrants_loc++;
          }
          marker_flag[qi] = SLABS_DISCR_AMR_MARKER_REFINE;
        }
        else if (SLABS_DISCR_AMR_MARKER_NONE == m) {
          if (SLABS_DISCR_AMR_MARKER_COARSEN == marker_flag[qi]) {
            marker->n_coarsen_quadrants_loc--;
            marker_flag[qi] = SLABS_DISCR_AMR_MARKER_NONE;
          }
        }
#ifdef YMIR_DEBUG
        else { YMIR_ASSERT (SLABS_DISCR_AMR_MARKER_COARSEN == m); }
#endif
      }
    }
  }
}

/**
 * Sets an AMR indicator into an AMR marker.
 */
void
slabs_discr_amr_marker_set_indicator (slabs_discr_amr_marker_t *marker,
                                      slabs_discr_amr_indicator_t *indicator)
{
  const int           set = 1;

  slabs_discr_amr_marker_set_merge_indicator (marker, indicator, set);
}

/**
 * Merges an AMR indicator into an AMR marker.
 */
void
slabs_discr_amr_marker_merge_indicator (slabs_discr_amr_marker_t *marker,
                                        slabs_discr_amr_indicator_t *indicator)
{
  const int           set = 0;

  slabs_discr_amr_marker_set_merge_indicator (marker, indicator, set);
}

/**
 * Computes the global number of quadrants marked for coarsening or
 * refinement; estimate the new number of quadrants with current AMR marker.
 */
void
slabs_discr_amr_marker_get_quadrant_counts (
                                       p4est_gloidx_t *n_coarsen_quadrants_glo,
                                       p4est_gloidx_t *n_refine_quadrants_glo,
                                       p4est_gloidx_t *n_quadrants_next,
                                       slabs_discr_amr_marker_t *marker)
{
  MPI_Comm            mpicomm = marker->p4est->mpicomm;
  int                 mpiret;
  int64_t             n_quads_loc[2];
  int64_t             n_quads_glo[2];

  /* set local counters */
  n_quads_loc[0] = (int64_t) marker->n_coarsen_quadrants_loc;
  n_quads_loc[1] = (int64_t) marker->n_refine_quadrants_loc;

  /* get processor-global counters */
  mpiret = MPI_Allreduce (&n_quads_loc, &n_quads_glo, 2, MPI_INT64_T, MPI_SUM,
                          mpicomm); YMIR_CHECK_MPI (mpiret);

  /* set output */
  if (n_coarsen_quadrants_glo != NULL) {
    *n_coarsen_quadrants_glo = (p4est_gloidx_t) n_quads_glo[0];
  }
  if (n_refine_quadrants_glo != NULL) {
    *n_refine_quadrants_glo = (p4est_gloidx_t) n_quads_glo[1];
  }
  if (n_quadrants_next != NULL) {
    const double        coarsen_factor = 1.0 / ((double) P4EST_CHILDREN) - 1.0;
    const double        n_coarsen = (double) n_quads_glo[0];
    const p4est_gloidx_t  refine_factor = P4EST_CHILDREN - 1;
    const p4est_gloidx_t  n_refine = (p4est_gloidx_t) n_quads_glo[1];

    *n_quadrants_next = marker->p4est->global_num_quadrants;
    *n_quadrants_next += (p4est_gloidx_t) round (coarsen_factor * n_coarsen);
    *n_quadrants_next += refine_factor * n_refine;
  }
}

/**
 * Sets the user data of each p4est quadrant to the corresponding AMR marker.
 */
void
slabs_discr_amr_marker_flag_p4est (slabs_discr_amr_marker_t *marker)
{
  slabs_discr_amr_marker_flag_t *marker_flag = marker->flag;
  p4est_t            *p4est = marker->p4est;

  p4est_topidx_t      ti;
  size_t              tqi;

  for (ti = p4est->first_local_tree; ti <= p4est->last_local_tree; ++ti) {
    p4est_tree_t       *t = p4est_tree_array_index (p4est->trees, ti);
    sc_array_t         *tquadrants = &t->quadrants;
    const size_t        tqoffset = t->quadrants_offset;

    for (tqi = 0; tqi < tquadrants->elem_count; ++tqi) {
      p4est_quadrant_t   *q = p4est_quadrant_array_index (tquadrants, tqi);
      slabs_discr_p4est_qdata_t *qdata =
        (slabs_discr_p4est_qdata_t *) q->p.user_data;

      qdata->flag = marker_flag[tqoffset + tqi];
    }
  }
}

/**
 * Callback function for initializing new p8est quadrants.
 */
static void
slabs_discr_amr_p4est_init (p4est_t *p4est, p4est_topidx_t tree,
                            p4est_quadrant_t *quadrant)
{
  slabs_discr_p4est_qdata_t *qdata =
    (slabs_discr_p4est_qdata_t *) quadrant->p.user_data;

  /* set to no coarsening or refinement */
  qdata->flag = SLABS_DISCR_AMR_MARKER_NONE;
}

/**
 * Callback function for p4est coarsening.
 */
static int
slabs_discr_amr_p4est_coarsen (p4est_t *p8est, p4est_topidx_t tree,
                               p4est_quadrant_t *quadrants[])
{
  int                 k;

  SC_CHECK_ABORT (p8est_quadrant_is_familypv (quadrants), "Coarsen invocation");

  for (k = 0; k < P4EST_CHILDREN; k++) {
    slabs_discr_p4est_qdata_t *qdata =
      (slabs_discr_p4est_qdata_t *) quadrants[k]->p.user_data;

    if (qdata->flag != SLABS_DISCR_AMR_MARKER_COARSEN) {
      /* if at least one child is not marked for coarsening */
      return 0;
    }
  }

  /* if all of the children are marked for coarsening */
  return 1;
}

/**
 * Callback function for p4est refinement.
 */
static int
slabs_discr_amr_p4est_refine (p4est_t *p8est, p4est_topidx_t tree,
                              p4est_quadrant_t *quadrant)
{
  slabs_discr_p4est_qdata_t *qdata =
    (slabs_discr_p4est_qdata_t *) quadrant->p.user_data;

  return (qdata->flag == SLABS_DISCR_AMR_MARKER_REFINE);
}

/**
 * TODO deprecated
 */
static p4est_gloidx_t
slabs_discr_amr_mark_wrt_indicator_overwrite (p4est_t *p4est,
                                              sc_dmatrix_t *indicator,
                                              double tol_min, double tol_max,
                                              int minlevel, int maxlevel,
                                              int overwrite)
{
  MPI_Comm            mpicomm = p4est->mpicomm;
  int                 mpiret;

  p4est_topidx_t      ti;
  p4est_tree_t       *tree;
  sc_array_t         *tquadrants;
  size_t              qoffset;
  size_t              qi;
  p4est_quadrant_t   *q;
  slabs_discr_p4est_qdata_t *qdata;
  double              value;
  p4est_gloidx_t      n_marked_quadrants_local = 0;
  p4est_gloidx_t      n_marked_quadrants_global;

  /* check input */
  YMIR_ASSERT (0 < indicator->m);
  YMIR_ASSERT (((p4est_locidx_t) indicator->n) == p4est->local_num_quadrants);
  YMIR_ASSERT (0 <= minlevel && minlevel <= maxlevel);

  for (ti = p4est->first_local_tree; ti <= p4est->last_local_tree; ++ti) {
    /* get tree */
    tree = p4est_tree_array_index (p4est->trees, ti);

    /* set number of quadrants of previous trees */
    qoffset = tree->quadrants_offset;

    /* get quadrants array for this tree */
    tquadrants = &tree->quadrants;

    for (qi = 0; qi < tquadrants->elem_count; ++qi) {
      /* get quadrant and its user data */
      q = p4est_quadrant_array_index (tquadrants, qi);
      qdata = (slabs_discr_p4est_qdata_t *) q->p.user_data;

      /* check if flag is valid */
      YMIR_ASSERT (overwrite ||
                   ( SLABS_DISCR_AMR_MARKER_COARSEN <= qdata->flag &&
                     qdata->flag <= SLABS_DISCR_AMR_MARKER_REFINE ));

      /* set value from `indicator` for this quadrant */
      value = indicator->e[0][qoffset + qi];

      if (0.0 <= value) {
        /* mark quadrant for coarsening or refinement */
        if ( value < tol_min && 0.0 < tol_min && minlevel < q->level ) {
          /* mark for coarsening if value is below tolerance and
           * min level is not reached */
          if (overwrite) {
            qdata->flag = SLABS_DISCR_AMR_MARKER_COARSEN;
            n_marked_quadrants_local++;
          }
        }
        else if ( tol_max < value && 0.0 < tol_max &&
                  (maxlevel <= 0 || q->level < maxlevel) ) {
          /* mark for refinement if value is above tolerance and
           * max level is not reached */
          if (overwrite || qdata->flag == SLABS_DISCR_AMR_MARKER_NONE) {
            n_marked_quadrants_local++;
          }
          qdata->flag = SLABS_DISCR_AMR_MARKER_REFINE;
        }
        else {
          if (overwrite) {
            qdata->flag = SLABS_DISCR_AMR_MARKER_NONE;
          }
          else if (!overwrite &&
                   qdata->flag == SLABS_DISCR_AMR_MARKER_COARSEN) {
            qdata->flag = SLABS_DISCR_AMR_MARKER_NONE;
            n_marked_quadrants_local--;
          }
        }
      }
    }
  }

  /* get processor-global number of marked quadrants */
  mpiret = MPI_Allreduce (&n_marked_quadrants_local, &n_marked_quadrants_global,
                          1, MPI_INT64_T, MPI_SUM, mpicomm);
  YMIR_CHECK_MPI (mpiret);

  /* return number of marked quadrants */
  return n_marked_quadrants_global;
}

/**
 *
 */
static p4est_gloidx_t
slabs_discr_amr_mark_wrt_indicator (p4est_t *p4est, sc_dmatrix_t *indicator,
                                    double tol_min, double tol_max,
                                    int minlevel, int maxlevel)
{
  return slabs_discr_amr_mark_wrt_indicator_overwrite (p4est, indicator,
                                                       tol_min, tol_max,
                                                       minlevel, maxlevel, 1);
}

/**
 *
 */
static inline void
slabs_discr_init_amr_visc_temp_elem (sc_dmatrix_t *visc_el_mat,
                                     const mangll_locidx_t elid,
                                     mangll_t *mangll,
                                     sc_dmatrix_t *temp_el_mat,
                                     sc_dmatrix_t *weak_el_mat,
                                     slabs_physics_options_t *physics_options,
                                     const int restrict_to_bounds)
{
  const int           visc_type = physics_options->viscosity_type;
  const int           visc_type_init =
                        physics_options->viscosity_type_for_init_nl_stokes;
  const double        visc_scaling = physics_options->viscosity_scaling;
  const double        upper_mantle_radius =
                        physics_options->viscosity_upper_mantle_radius;
  const double        visc_lower_mantle_scaling =
                        physics_options->viscosity_lower_mantle_scaling;
  const int           n_nodes_per_el = visc_el_mat->m;
  const double       *_sc_restrict x = mangll->X->e[elid];
  const double       *_sc_restrict y = mangll->Y->e[elid];
  const double       *_sc_restrict z = mangll->Z->e[elid];
  double             *temp_el_data = temp_el_mat->e[0];
  int                 nodeid;

  /* compute temperature field of this element at discontinuous GLL nodes */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    slabs_physics_temperature_set_fn (&temp_el_data[nodeid],
                                      x[nodeid], y[nodeid], z[nodeid], 0,
                                      physics_options);
  }

  /* get weak zone of this element */
  if (weak_el_mat != NULL) {
    slabs_weak_elem (weak_el_mat, x, y, z, mangll->refel->Vmask,
                     physics_options);
  }

  /* change scaling of viscosity in upper mantle
   * (snippet copied from slabs_physics.c) */
  if (   visc_type == SL_VISCOSITY_NONLINEAR
      && visc_type_init == SL_VISCOSITY_INIT_NL_STOKES_TEMP_UM_REL_TO_LM
      && 0.0 < upper_mantle_radius) {
    physics_options->viscosity_scaling = SL_VISCOSITY_UM_REL_TO_LM
                                         * visc_lower_mantle_scaling;
  }

  /* compute temperature dependent visosity for this element
   * (assume here: coordinates are used by `slabs_visc_temp_elem` to determine
   * position in upper/lower mantle only) */
  slabs_visc_temp_elem (visc_el_mat, x, y, z, mangll->refel->Vmask,
                        temp_el_mat, weak_el_mat, physics_options,
                        restrict_to_bounds);

  /* restore viscosity scaling (snippet copied from slabs_physics.c) */
  if (   visc_type == SL_VISCOSITY_NONLINEAR
      && visc_type_init == SL_VISCOSITY_INIT_NL_STOKES_TEMP_UM_REL_TO_LM
      && 0.0 < upper_mantle_radius) {
    physics_options->viscosity_scaling = visc_scaling;
  }
}

/**
 *
 */
static inline void
slabs_discr_init_amr_rhs_elem (sc_dmatrix_t *rhs_el_mat,
                               const double *x, const double *y,
                               const double *z,
                               sc_dmatrix_t *temp_el_mat,
                               slabs_physics_options_t *physics_options)
{
  double             *temp_el_data = temp_el_mat->e[0];
  const int           n_nodes_per_el = rhs_el_mat->m;
  int                 nodeid;

  /* compute temperature field of this element at discontinuous GLL nodes */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    slabs_physics_temperature_set_fn (&temp_el_data[nodeid],
                                      x[nodeid], y[nodeid], z[nodeid], 0,
                                      physics_options);
  }

  /* compute right-hand side for this element */
  slabs_rhs_elem (rhs_el_mat, x, y, z, temp_el_mat, physics_options);
}

/**
 *
 */
sc_dmatrix_t *
slabs_discr_init_amr_indicator_new (mangll_t *mangll,
                                    slabs_physics_options_t *physics_options,
                                    slabs_discr_options_t *discr_options,
                                    int indicator_type)
{
  const double        domain_norm = discr_options->domain_size_normalization;
  const int           amr_lower_mantle = discr_options->init_amr_lower_mantle;
  const double        visc_max = SC_MAX (0.0, physics_options->viscosity_max);
  const double        rhs_norm_shift = discr_options->init_amr_rhs_norm_shift;
  const mangll_locidx_t  n_elements = mangll->mesh->K;
  const int           N = ymir_n (mangll->N);
  const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);
  const int          *Vmask = mangll->refel->Vmask;

  sc_dmatrix_t       *indicator, *indicator_plates;
  double             *indicator_data;
  mangll_locidx_t     elid;
  int                 nodeid, fieldid;
  double             *x, *y, *z;
  double              min, max;

  sc_dmatrix_t       *temp_el_mat = NULL;
  sc_dmatrix_t       *weak_el_mat = NULL;
  sc_dmatrix_t       *visc_el_mat = NULL;
  sc_dmatrix_t       *rhs_el_mat = NULL;
  sc_dmatrix_t       *rhs_norm_el_mat = NULL;
  sc_dmatrix_t       *grad_el_mat = NULL;
  sc_dmatrix_t       *grad_tens_el_mat = NULL;
  sc_dmatrix_t       *tmp_vec = NULL;

  /* check input */
  YMIR_ASSERT (physics_options->temperature_type != SL_TEMP_IMPORT_FILE);
  YMIR_ASSERT (physics_options->weakzone_type != SL_WEAKZONE_IMPORT_FILE);

  /* create new refinement indicator vector:
   *   [0] indicator value, [1] inside plate
   */
  indicator = sc_dmatrix_new (2, n_elements);
  indicator_data = indicator->e[0];
  indicator_plates = sc_dmatrix_new_view_offset (1, 1, n_elements, indicator);
  sc_dmatrix_set_value (indicator_plates, -1.0);

  /* return zero indicator if nothing to do */
  if (   indicator_type == SL_AMR_INDICATOR_NONE
      || indicator_type == SL_AMR_INDICATOR_COARSEN_ALL ) {
    sc_dmatrix_set_zero (indicator);
    sc_dmatrix_destroy (indicator_plates);
    return indicator;
  }

  /* return one indicator if refine all */
  if (indicator_type == SL_AMR_INDICATOR_REFINE_ALL) {
    sc_dmatrix_set_value (indicator, 1.0);
    sc_dmatrix_destroy (indicator_plates);
    return indicator;
  }

  /* create work variables */
  if (   indicator_type == SL_AMR_INDICATOR_VISC_DR
      || indicator_type == SL_AMR_INDICATOR_VISC_GRAD
      || indicator_type == SL_AMR_INDICATOR_VISC_PECLET
      || indicator_type == SL_AMR_INDICATOR_RHS_OVER_VISC
      || indicator_type == SL_AMR_INDICATOR_RHS_OVER_VISC_DR
      || indicator_type == SL_AMR_INDICATOR_RHS_OVER_VISC_GRAD
      || indicator_type == SL_AMR_INDICATOR_WEAK_GRAD ) {
    temp_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
    visc_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  }
  else if (   indicator_type == SL_AMR_INDICATOR_RHS_MAGNITUDE
           || indicator_type == SL_AMR_INDICATOR_RHS_PECLET ) {
    temp_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  }
  if (   indicator_type == SL_AMR_INDICATOR_VISC_DR
      || indicator_type == SL_AMR_INDICATOR_VISC_GRAD
      || indicator_type == SL_AMR_INDICATOR_VISC_PECLET
      || indicator_type == SL_AMR_INDICATOR_RHS_OVER_VISC
      || indicator_type == SL_AMR_INDICATOR_RHS_OVER_VISC_DR
      || indicator_type == SL_AMR_INDICATOR_RHS_OVER_VISC_GRAD
      || indicator_type == SL_AMR_INDICATOR_WEAK_DR
      || indicator_type == SL_AMR_INDICATOR_WEAK_GRAD ) {
    weak_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  }
  if (   indicator_type == SL_AMR_INDICATOR_RHS_OVER_VISC
      || indicator_type == SL_AMR_INDICATOR_RHS_OVER_VISC_DR
      || indicator_type == SL_AMR_INDICATOR_RHS_OVER_VISC_GRAD ) {
    rhs_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
    rhs_norm_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  }
  else if (   indicator_type == SL_AMR_INDICATOR_RHS_MAGNITUDE
           || indicator_type == SL_AMR_INDICATOR_RHS_PECLET ) {
    rhs_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
  }
  if (   indicator_type == SL_AMR_INDICATOR_VISC_GRAD
      || indicator_type == SL_AMR_INDICATOR_VISC_PECLET
      || indicator_type == SL_AMR_INDICATOR_RHS_OVER_VISC_GRAD
      || indicator_type == SL_AMR_INDICATOR_WEAK_GRAD ) {
    grad_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
  }
  if (indicator_type == SL_AMR_INDICATOR_RHS_PECLET) {
    grad_tens_el_mat = sc_dmatrix_new (n_nodes_per_el, 9);
    tmp_vec = sc_dmatrix_new (n_nodes_per_el, 3);
  }

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    /* set element coordinates of physical space at GLL nodes */
    x = mangll->X->e[elid];
    y = mangll->Y->e[elid];
    z = mangll->Z->e[elid];

    /* if no AMR in lower mantle */
    if (!amr_lower_mantle &&
        !slabs_physics_elem_in_upper_mantle (x, y, z, Vmask, physics_options)) {
      indicator_data[elid] = 0.0;
      continue;
    }

    /* compute indicator */
    switch (indicator_type) {
    case SL_AMR_INDICATOR_VISC_DR: /* dynamic range of viscosity */
      /* compute temperature dependent visosity for this element */
      slabs_discr_init_amr_visc_temp_elem (visc_el_mat, elid, mangll,
                                           temp_el_mat, weak_el_mat,
                                           physics_options, 1);

      /* find min and max viscosity */
      slabs_matrix_compute_abs_min_max (&min, &max, visc_el_mat);

      /* set dynamic range for this element */
      indicator_data[elid] = max / min;
      break;

    case SL_AMR_INDICATOR_VISC_GRAD: /* max norm of viscosity gradient */
      /* compute temperature dependent visosity for this element */
      slabs_discr_init_amr_visc_temp_elem (visc_el_mat, elid, mangll,
                                           temp_el_mat, weak_el_mat,
                                           physics_options, 1);

      /* compute gradient (overwrite `temp_el_mat`) */
      slabs_gradient_gll_to_gll_elem (visc_el_mat, grad_el_mat, mangll,
                                      elid, temp_el_mat);

      /* calculate max norm of gradient */
      max = slabs_matrix_compute_abs_max (grad_el_mat);

      /* set max norm weighted with element size */
      indicator_data[elid] = max * slabs_elem_volume (mangll, elid)
                                 * domain_norm;
      break;

    case SL_AMR_INDICATOR_VISC_PECLET: /* max norm of visc grad over visc */
      /* compute temperature dependent visosity for this element */
      slabs_discr_init_amr_visc_temp_elem (visc_el_mat, elid, mangll,
                                           temp_el_mat, weak_el_mat,
                                           physics_options, 1);

      /* compute gradient (overwrite `temp_el_mat`) */
      slabs_gradient_gll_to_gll_elem (visc_el_mat, grad_el_mat, mangll,
                                      elid, temp_el_mat);

      /* divide in viscosity (node-wise) */
      for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
        for (fieldid = 0; fieldid < 3; fieldid++) {
          grad_el_mat->e[nodeid][fieldid] /= visc_el_mat->e[nodeid][0];
        }
      }

      /* calculate max norm of gradient */
      max = slabs_matrix_compute_abs_max (grad_el_mat);

      /* set max norm weighted with element size */
      indicator_data[elid] = max * slabs_elem_volume (mangll, elid)
                                 * domain_norm;
      break;

    case SL_AMR_INDICATOR_RHS_OVER_VISC: /* quotient (rhs + 1)/visc */
      /* compute temperature dependent visosity for this element */
      slabs_discr_init_amr_visc_temp_elem (visc_el_mat, elid, mangll,
                                           temp_el_mat, NULL,
                                           physics_options, 0);

      /* compute right-hand side for this element */
      slabs_discr_init_amr_rhs_elem (rhs_el_mat, x, y, z, temp_el_mat,
                                     physics_options);

      /* compute max-norm of right-hand side */
      for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
        const double          *rhs = rhs_el_mat->e[nodeid];

        rhs_norm_el_mat->e[nodeid][0] =
          SC_MAX ( fabs (rhs[0]), SC_MAX (fabs (rhs[1]), fabs (rhs[2])) );

        if (0.0 < rhs_norm_shift) {
          rhs_norm_el_mat->e[nodeid][0] += rhs_norm_shift;
        }
      }

      /* divide in viscosity */
      sc_dmatrix_dotdivide (rhs_norm_el_mat, visc_el_mat);

      /* calculate max norm of quotient */
      indicator_data[elid] = slabs_matrix_compute_abs_max (rhs_norm_el_mat);
      break;

    case SL_AMR_INDICATOR_RHS_OVER_VISC_DR: /* DR of (rhs + 1)/visc */
      /* compute temperature dependent visosity for this element */
      slabs_discr_init_amr_visc_temp_elem (visc_el_mat, elid, mangll,
                                           temp_el_mat, NULL,
                                           physics_options, 0);

      /* compute right-hand side for this element */
      slabs_discr_init_amr_rhs_elem (rhs_el_mat, x, y, z, temp_el_mat,
                                     physics_options);

      /* compute max-norm of right-hand side */
      for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
        const double          *rhs = rhs_el_mat->e[nodeid];

        rhs_norm_el_mat->e[nodeid][0] =
          SC_MAX ( fabs (rhs[0]), SC_MAX (fabs (rhs[1]), fabs (rhs[2])) );

        if (0.0 < rhs_norm_shift) {
          rhs_norm_el_mat->e[nodeid][0] += rhs_norm_shift;
        }
      }

      /* divide in viscosity */
      sc_dmatrix_dotdivide (rhs_norm_el_mat, visc_el_mat);

      /* find min and max */
      slabs_matrix_compute_abs_min_max (&min, &max, rhs_norm_el_mat);

      /* set dynamic range for this element */
      indicator_data[elid] = max / min;
      break;

    case SL_AMR_INDICATOR_RHS_OVER_VISC_GRAD: /* grad of (rhs + 1)/visc */
      /* compute temperature dependent visosity for this element */
      slabs_discr_init_amr_visc_temp_elem (visc_el_mat, elid, mangll,
                                           temp_el_mat, NULL,
                                           physics_options, 0);

      /* compute right-hand side for this element */
      slabs_discr_init_amr_rhs_elem (rhs_el_mat, x, y, z, temp_el_mat,
                                     physics_options);

      /* compute max-norm of right-hand side */
      for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
        const double          *rhs = rhs_el_mat->e[nodeid];

        rhs_norm_el_mat->e[nodeid][0] =
          SC_MAX ( fabs (rhs[0]), SC_MAX (fabs (rhs[1]), fabs (rhs[2])) );

        if (0.0 < rhs_norm_shift) {
          rhs_norm_el_mat->e[nodeid][0] += rhs_norm_shift;
        }
      }

      /* divide in viscosity */
      sc_dmatrix_dotdivide (rhs_norm_el_mat, visc_el_mat);

      /* compute gradient (overwrite `temp_el_mat`) */
      slabs_gradient_gll_to_gll_elem (rhs_norm_el_mat, grad_el_mat, mangll,
                                      elid, temp_el_mat);

      /* calculate max norm of gradient */
      max = slabs_matrix_compute_abs_max (grad_el_mat);

      /* set max norm weighted with element size */
      indicator_data[elid] = max * slabs_elem_volume (mangll, elid)
                                 * domain_norm;

      /* set indicator if inside plates */
      if (   0.0 < visc_max
          && visc_max < slabs_matrix_compute_abs_max (visc_el_mat)) {
        indicator_plates->e[0][elid] = indicator_data[elid];
      }
      break;

    case SL_AMR_INDICATOR_WEAK_DR: /* dynamic range of weak zone */
      /* compute weak zone for this element */
      slabs_weak_elem (weak_el_mat, x, y, z, Vmask, physics_options);

      /* find min and max weak zone factor */
      slabs_matrix_compute_abs_min_max (&min, &max, weak_el_mat);

      /* set dynamic range for this element */
      indicator_data[elid] = max / min;
      break;

    case SL_AMR_INDICATOR_WEAK_GRAD: /* max norm of weak zone gradient */
      /* compute weak zone for this element */
      slabs_weak_elem (weak_el_mat, x, y, z, Vmask, physics_options);

      /* compute gradient (overwrite `temp_el_mat`) */
      slabs_gradient_gll_to_gll_elem (weak_el_mat, grad_el_mat, mangll,
                                      elid, temp_el_mat);

      /* calculate max norm of gradient */
      max = slabs_matrix_compute_abs_max (grad_el_mat);

      /* set max norm weighted with element size */
      indicator_data[elid] = max * slabs_elem_volume (mangll, elid)
                                 * domain_norm;
      break;

    case SL_AMR_INDICATOR_WEAK_SUBDU_DIST: /* refine within distance */
      if (physics_options->weakzone_type == SL_WEAKZONE_2PLATES_POLY2) {
        const double        weak_subdu_maxdist =
                              discr_options->init_amr_weak_subdu_maxdist;
        double              subdu_lon;
        double              subdu_dip_angle;
        double              subdu_depth, subdu_width;
        double              subdu_thickness, subdu_thickness_const;
        double              courtesy_width;
        double              total_thickness;
        double              start_node, start_val, start_deriv;
        double              end_node, end_val;
        double              r, lon, dist;

        /* set parameters according to physics options */
        subdu_lon = physics_options->weakzone_2plates_subdu_longitude;
        subdu_dip_angle = physics_options->weakzone_2plates_subdu_dip_angle;
        subdu_depth = physics_options->weakzone_2plates_subdu_depth;
        subdu_width = physics_options->weakzone_2plates_subdu_width;
        subdu_thickness = physics_options->weakzone_2plates_subdu_thickness;
        subdu_thickness_const =
          physics_options->weakzone_2plates_subdu_thickness_const;
        courtesy_width = subdu_thickness / SL_EARTH_RADIUS;
        total_thickness = (2.0 * subdu_thickness - subdu_thickness_const)
                          / SL_EARTH_RADIUS;

        /* set points for polynomial interpolation */
        start_node = subdu_lon;
        start_val = SL_SHELL_RADIUS_TOP;
        start_deriv = tan (-subdu_dip_angle / 180.0 * M_PI);
        end_node = start_node + subdu_width / SL_EARTH_RADIUS
                              * SL_SHELL_RADIUS_TOP;
        end_val = start_val - subdu_depth / SL_EARTH_RADIUS;

        /* find minimum distance from a node of this element to curve */
        min = DBL_MAX;
        for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
          /* compute radius and longitude */
          r = slabs_compute_radius (x[nodeid], y[nodeid], z[nodeid],
                                    physics_options);
          lon = slabs_compute_longitude (x[nodeid], y[nodeid], z[nodeid],
                                         physics_options);

          /* only consider point in a rectangle containing the weak zone */
          if (   (  start_node - 0.5 * total_thickness
                  / sin (subdu_dip_angle / 180.0 * M_PI) - courtesy_width )
                 <= lon
              && lon <= (end_node + 0.5 * total_thickness + courtesy_width)
              && (end_val - total_thickness - courtesy_width) <= r ) {
            /* compute distance to weak zone */
            dist = slabs_2plates_weakzone_poly2_subdu_dist (
                r, lon, start_node, start_val, start_deriv, end_node,
                end_val);
          }
          else {
            dist = weak_subdu_maxdist;
          }

          /* update min distance */
          min = SC_MIN (min, dist);
        }

        /* set minimum distance to boundary of refinement layer */
        indicator_data[elid] = SC_MAX (0.0, weak_subdu_maxdist - min);
      }
      else {
        /* this AMR type is not available for other weak zones */
        YMIR_ABORT_NOT_REACHED ();
      }
      break;

    case SL_AMR_INDICATOR_WEAK_RIDGE_DIST: /* refine within distance */
      if (physics_options->weakzone_type == SL_WEAKZONE_2PLATES_POLY2) {
        const double        weak_ridge_maxdist =
                              discr_options->init_amr_weak_ridge_maxdist;
        double              ridge_depth, ridge_width;
        double              ridge_smoothwidth;
        double              courtesy_width;
        double              lon_min = physics_options->domain_lon_min;
        double              end_node, end_val;
        double              r, lon, dist;

        /* set parameters according to physics options */
        ridge_depth = physics_options->weakzone_2plates_ridge_depth;
        ridge_width = physics_options->weakzone_2plates_ridge_width;
        ridge_smoothwidth =
          physics_options->weakzone_2plates_ridge_smoothwidth;
        courtesy_width = 2.0 * ridge_smoothwidth / SL_EARTH_RADIUS;

        /* set bottom left corner of weak zone */
        end_node = lon_min + ridge_width / SL_EARTH_RADIUS
                           * SL_SHELL_RADIUS_TOP;
        end_val = SL_SHELL_RADIUS_TOP - ridge_depth / SL_EARTH_RADIUS;

        /* find minimum distance from a node of this element to curve */
        min = DBL_MAX;
        for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
          /* compute radius and longitude */
          r = slabs_compute_radius (x[nodeid], y[nodeid], z[nodeid],
                                    physics_options);
          lon = slabs_compute_longitude (x[nodeid], y[nodeid], z[nodeid],
                                         physics_options);

          /* only consider points close to weak zone */
          if (   lon <= (end_node + courtesy_width)
              && (end_val - courtesy_width) <= r) {
            /* compute distance to weak zone */
            dist = slabs_2plates_weakzone_poly2_ridge_dist (
                r, lon, end_node, end_val);

            if (dist < DBL_MIN) { /* if inside weak zone */
              dist = -1.0;
            }
          }
          else {
            dist = weak_ridge_maxdist;
          }

          /* update min distance */
          min = SC_MIN (min, dist);
        }

        /* set minimum distance to boundary of refinement layer */
        indicator_data[elid] = SC_MAX (0.0, weak_ridge_maxdist - min);
      }
      else {
        /* this AMR type is not available for other weak zones */
        YMIR_ABORT_NOT_REACHED ();
      }
      break;

    case SL_AMR_INDICATOR_RHS_MAGNITUDE: /* magnitude of rhs */
      {
        const double        rhs_scaling = fabs (physics_options->rhs_scaling);

        /* compute right-hand side for this element */
        slabs_discr_init_amr_rhs_elem (rhs_el_mat, x, y, z, temp_el_mat,
                                       physics_options);

        /* compute max indicator for this element */
        max = 0.0;
        for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
          const double       *rhs = rhs_el_mat->e[nodeid];
          double              rhs_norm;

          /* compute Euclidian norm of right-hand side at this node */
          rhs_norm = sqrt (rhs[0]*rhs[0] + rhs[1]*rhs[1] + rhs[2]*rhs[2]);

          /* normalize right-hand side by scaling value; now assume
           * 0 <= rhs_norm <= 1 */
          rhs_norm /= rhs_scaling;

          /* update max */
          max = SC_MAX (max, rhs_norm);
        }

        /* set max norm weighted with element size */
        indicator_data[elid] = max * slabs_elem_volume (mangll, elid)
                                   * domain_norm;
      }
      break;

    case SL_AMR_INDICATOR_RHS_PECLET: /* max norm of rhs grad over rhs */
      /* compute right-hand side for this element */
      slabs_discr_init_amr_rhs_elem (rhs_el_mat, x, y, z, temp_el_mat,
                                     physics_options);

      /* compute gradient */
      slabs_gradient_gll_to_gll_elem (rhs_el_mat, grad_tens_el_mat, mangll,
                                      elid, tmp_vec);

      /* compute max indicator for this element */
      max = 0.0;
      for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
        const double       *rhs = rhs_el_mat->e[nodeid];
        double              rhs_norm;

        /* compute Euclidian norm of right-hand side at this node */
        rhs_norm = sqrt (rhs[0]*rhs[0] + rhs[1]*rhs[1] + rhs[2]*rhs[2]);

        if (SC_1000_EPS < rhs_norm) {
          const double       *grad = grad_tens_el_mat->e[nodeid];
          double              grad_norm;

          /* shift norm of right-hand side */
          if (0.0 < rhs_norm_shift) {
            rhs_norm += rhs_norm_shift;
          }

          /* compute Frobenius norm of gradient matrix at this node */
          grad_norm = 0.0;
          for (fieldid = 0; fieldid < 9; fieldid++) {
            grad_norm += grad[fieldid] * grad[fieldid];
          }
          grad_norm = sqrt (grad_norm);

          /* update max */
          max = SC_MAX (max, grad_norm / rhs_norm);
        }
      }

      /* set max norm weighted with element size */
      indicator_data[elid] = max * slabs_elem_volume (mangll, elid)
                                 * domain_norm;
      break;

    default: /* unknown AMR indicator type */
      YMIR_ABORT_NOT_REACHED ();
    }
  }

  /* destroy work variables */
  if (temp_el_mat != NULL) {
    sc_dmatrix_destroy (temp_el_mat);
  }
  if (weak_el_mat != NULL) {
    sc_dmatrix_destroy (weak_el_mat);
  }
  if (visc_el_mat != NULL) {
    sc_dmatrix_destroy (visc_el_mat);
  }
  if (rhs_el_mat != NULL) {
    sc_dmatrix_destroy (rhs_el_mat);
  }
  if (rhs_norm_el_mat != NULL) {
    sc_dmatrix_destroy (rhs_norm_el_mat);
  }
  if (grad_el_mat != NULL) {
    sc_dmatrix_destroy (grad_el_mat);
  }
  if (grad_tens_el_mat != NULL) {
    sc_dmatrix_destroy (grad_tens_el_mat);
  }
  if (tmp_vec != NULL) {
    sc_dmatrix_destroy (tmp_vec);
  }

  /* destroy view */
  sc_dmatrix_destroy (indicator_plates);

  /* return refinement indicator vector */
  return indicator;
}

/**
 * AMR and partitioning with respect to a refinement indicator vector.
 */
static void
slabs_discr_init_amr_indicator (p8est_t *p8est,
                                slabs_physics_options_t *physics_options,
                                slabs_discr_options_t *discr_options,
                                const int indicator_type,
                                const double tol_min, const double tol_max,
                                const double inside_plates_tol_min,
                                const double inside_plates_tol_max,
                                int max_steps,
                                const int minlevel, const int maxlevel)
{
  const char         *this_fn_name = "slabs_discr_init_amr_indicator";

  const int           uniform =
    (indicator_type == SL_AMR_INDICATOR_REFINE_ALL);
  mangll_t           *mangll;
  int                 amr_step;
  p4est_locidx_t      n_quadrants_loc;
  p4est_gloidx_t      n_marked_quadrants;
  sc_dmatrix_t       *indicator, *indicator_global, *indicator_plates;

  YMIR_GLOBAL_INFOF ("Into %s (indicator %i)\n", this_fn_name, indicator_type);

  /* check number of max AMR steps */
  if (max_steps <= 0) {
    max_steps = P4EST_MAXLEVEL;
  }

  for (amr_step = 0; amr_step < max_steps; amr_step++) {
    n_quadrants_loc = p8est->local_num_quadrants;

    /* print mesh statistics */
    ymir_monitor_print_global_element_stats (p8est);

    /* print memory usage */
    ymir_monitor_print_global_mem_usage (p8est->mpicomm);

    /*
     * Create AMR indicator and mark elements for refinement or coarsening
     */

    /* create mangll */
    slabs_discr_mangll_and_cnodes_new (&mangll, NULL, p8est, discr_options);

    /* compute refinement indicators */
    indicator = slabs_discr_init_amr_indicator_new (mangll, physics_options,
                                                    discr_options,
                                                    indicator_type);
    indicator_global = sc_dmatrix_new_view_offset (0, 1, n_quadrants_loc,
                                                   indicator);
    indicator_plates = sc_dmatrix_new_view_offset (1, 1, n_quadrants_loc,
                                                   indicator);

    /* destroy mangll */
    mangll_destroy (mangll);

    /* mark quadrants for coarsening or refinement */
    if ( (0.0 <= tol_min || 0.0 <= tol_max) && tol_min < tol_max ) {
      YMIR_GLOBAL_INFOF ("%s: Step %i, mark globally w.r.t. indicator %i, "
                         "tol [%g,%g]\n", this_fn_name, amr_step,
                         indicator_type, tol_min, tol_max);

      n_marked_quadrants = slabs_discr_amr_mark_wrt_indicator (
          p8est, indicator_global, tol_min, tol_max, minlevel, maxlevel);
    }
    else {
      n_marked_quadrants = 0;
    }

    /* mark quadrants inside plates for coarsening or refinement */
    if (   (0.0 <= inside_plates_tol_min || 0.0 <= inside_plates_tol_max)
        && inside_plates_tol_min < inside_plates_tol_max ) {
      YMIR_GLOBAL_INFOF ("%s: Step %i, mark in plates w.r.t. indicator %i, "
                         "tol [%g,%g]\n",
                         this_fn_name, amr_step, indicator_type,
                         inside_plates_tol_min, inside_plates_tol_max);

      n_marked_quadrants += slabs_discr_amr_mark_wrt_indicator (
          p8est, indicator_plates, inside_plates_tol_min, inside_plates_tol_max,
          minlevel, maxlevel);
    }

    YMIR_GLOBAL_INFOF ("%s: Step %i, %lli quadrants marked for AMR\n",
                       this_fn_name, amr_step,
                       (long long int) n_marked_quadrants);

    /* exit AMR loop if no elements were marked for refinement or coarsening */
    if (!n_marked_quadrants) {
      break;
    }
    else {
      /* destroy indicators if this was not the last AMR step */
      if (amr_step < (max_steps - 1)) {
        sc_dmatrix_destroy (indicator_global);
        sc_dmatrix_destroy (indicator_plates);
        sc_dmatrix_destroy (indicator);
      }
    }

    /*
     * Adapt & partition p8est mesh
     */

    /* coarsen and refine p8est */
    p8est_coarsen (p8est, 0, slabs_discr_amr_p4est_coarsen,
                   slabs_discr_amr_p4est_init);
    p8est_refine (p8est, 0, slabs_discr_amr_p4est_refine, NULL);

    if (!uniform) {
      /* balance p8est */
      p8est_balance (p8est, P8EST_CONNECT_FULL, NULL);

      /* partition p8est */
      switch (discr_options->mesh_partitioning_type) {
      case SL_MESH_PARTITIONING_ELEM:
        p8est_partition_ext (p8est, 1, NULL);
        break;

      case SL_MESH_PARTITIONING_DOF_VEL:
        YMIR_ABORT_NOT_REACHED ();
        /*TODO delete deprecated code below and option mesh_partitioning_type
        p8est_partition_lnodes_ext (
            p8est, NULL,
            3 * discr_options->n_vel_dnodes_per_el_interior,
            3 * discr_options->n_vel_dnodes_per_face_interior,
            3 * discr_options->n_vel_dnodes_per_edge_interior,
            3 * discr_options->n_vel_dnodes_per_corner,
            1);
        */
        break;

      case SL_MESH_PARTITIONING_DOF_VEL_PRESS:
        YMIR_ABORT_NOT_REACHED ();
        /*TODO delete deprecated code below and option mesh_partitioning_type
        p8est_partition_lnodes_ext (
            p8est, NULL,
            3 * discr_options->n_vel_dnodes_per_el_interior +
            discr_options->n_press_dnodes_per_el,
            3 * discr_options->n_vel_dnodes_per_face_interior,
            3 * discr_options->n_vel_dnodes_per_edge_interior,
            3 * discr_options->n_vel_dnodes_per_corner,
            1);
        */
        break;

      default: /* unknown partitioning type */
        YMIR_ABORT_NOT_REACHED ();
      }
    }
  }

  /*
   * Status output
   */

#ifdef YMIR_DEBUG
  {
    MPI_Comm            mpicomm = p8est->mpicomm;
    int                 mpiret;

    p4est_topidx_t      ti;
    p4est_tree_t       *tree;
    mangll_locidx_t     elid;
    int8_t              maxlevel_local;
    int8_t              maxlevel_global = 0;
    double              max_local;
    double              max_global = 0;
    double             *indicator_data;

    /* find processor local max level of all forests */
    maxlevel_local = 0;
    for (ti = p8est->first_local_tree; ti <= p8est->last_local_tree; ++ti) {
      tree = p4est_tree_array_index (p8est->trees, ti);
      maxlevel_local = SC_MAX (maxlevel_local, tree->maxlevel);
    }

    /* get processor global max level */
    mpiret = MPI_Allreduce (&maxlevel_local, &maxlevel_global, 1, MPI_INT8_T,
                            MPI_MAX, mpicomm);
    YMIR_CHECK_MPI (mpiret);

    /* output */
    YMIR_GLOBAL_INFOF ("%s: %i AMR steps, indicator %i, max level %i\n",
                       this_fn_name, amr_step, indicator_type,
                       (int) maxlevel_global);

    /* get processor local max indicator value */
    indicator_data = indicator_global->e[0];
    max_local = 0.0;
    for (elid = 0; elid < indicator_global->n; elid++) {
      max_local = SC_MAX (max_local, indicator_data[elid]);
    }

    /* get processor global max indicator value */
    mpiret = MPI_Allreduce (&max_local, &max_global, 1, MPI_DOUBLE, MPI_MAX,
                            mpicomm);
    YMIR_CHECK_MPI (mpiret);

    /* output */
    YMIR_GLOBAL_INFOF ("%s: %i AMR steps, indicator %i, "
                       "max global indicator %1.3e, tol [%g,%g]\n",
                       this_fn_name, amr_step, indicator_type, max_global,
                       tol_min, tol_max);

    if (   (0.0 < inside_plates_tol_min || 0.0 < inside_plates_tol_max)
        && inside_plates_tol_min < inside_plates_tol_max ) {
      /* get processor local max indicator value */
      indicator_data = indicator_plates->e[0];
      max_local = 0.0;
      for (elid = 0; elid < indicator_plates->n; elid++) {
        max_local = SC_MAX (max_local, indicator_data[elid]);
      }

      /* get processor global max indicator value */
      mpiret = MPI_Allreduce (&max_local, &max_global, 1, MPI_DOUBLE, MPI_MAX,
                              mpicomm);
      YMIR_CHECK_MPI (mpiret);

      /* output */
      YMIR_GLOBAL_INFOF ("%s: %i AMR steps, indicator %i, "
                         "max plates indicator %1.3e, tol [%g,%g]\n",
                         this_fn_name, amr_step, indicator_type, max_global,
                         inside_plates_tol_min, inside_plates_tol_max);
    }
  }
#endif

  YMIR_GLOBAL_INFOF ("Done %s (indicator %i) after %i AMR steps\n",
                     this_fn_name, indicator_type, amr_step);

  /* destroy */
  sc_dmatrix_destroy (indicator_global);
  sc_dmatrix_destroy (indicator_plates);
  sc_dmatrix_destroy (indicator);
}

/**
 * Inital AMR and partitioning of p4est mesh.
 */
void
slabs_discr_initial_amr_no_interp (p8est_t *p8est,
                                   slabs_physics_options_t *physics_options,
                                   slabs_discr_options_t *discr_options)
{
  const char         *this_fn_name = "slabs_discr_initial_amr_no_interp";
  const int           order = discr_options->order;
  const int           minlevel = discr_options->minlevel;
  const int           maxlevel = discr_options->maxlevel;

  const int           amr_override_order =
                        discr_options->init_amr_override_order;
  const int           amr_max_steps = discr_options->init_amr_max_steps;

  const int           visc_indicator =
                        discr_options->init_amr_visc_indicator_type;
  const double        visc_tol_min = discr_options->init_amr_visc_tol_min;
  const double        visc_tol_max = discr_options->init_amr_visc_tol_max;
  const double        visc_in_plates_tol_min =
                        discr_options->init_amr_visc_inside_plates_tol_min;
  const double        visc_in_plates_tol_max =
                        discr_options->init_amr_visc_inside_plates_tol_max;

  const int           weak_subdu_indicator =
                        discr_options->init_amr_weak_subdu_indicator_type;
  double              weak_subdu_tol_min =
                        discr_options->init_amr_weak_subdu_tol_min;
  double              weak_subdu_tol_max =
                        discr_options->init_amr_weak_subdu_tol_max;

  const int           weak_ridge_indicator =
                        discr_options->init_amr_weak_ridge_indicator_type;
  double              weak_ridge_tol_min =
                        discr_options->init_amr_weak_ridge_tol_min;
  double              weak_ridge_tol_max =
                        discr_options->init_amr_weak_ridge_tol_max;

  const int           rhs_indicator =
                        discr_options->init_amr_rhs_indicator_type;
  const double        rhs_tol_min = discr_options->init_amr_rhs_tol_min;
  const double        rhs_tol_max = discr_options->init_amr_rhs_tol_max;

  const int           post_uniform_n_steps =
                        discr_options->init_amr_post_uniform_n_steps;

  /* change polynomial order of discretization */
  if (0 < amr_override_order && order != amr_override_order) {
    discr_options->order = amr_override_order;
    YMIR_GLOBAL_INFOF ("%s: Override polynomial discr order for init AMR "
                       "(%i instead of %i)\n",
                       this_fn_name, amr_override_order, order);
  }

  /* refine mesh radially or in z-direction */
  slabs_discr_refine_wrt_depth (p8est, physics_options, discr_options);

  /* refine mesh at lower/upper mantle interface */
  slabs_discr_refine_at_layer (p8est, physics_options, discr_options);

  /* refine mesh at surface */
  slabs_discr_refine_at_surface (p8est, physics_options, discr_options);

  /* refine mesh w.r.t. subduction weak zone */
  if (   weak_subdu_indicator != SL_AMR_INDICATOR_NONE
      && SL_AMR_INDICATOR_WEAK_DR <= weak_subdu_indicator
      && 0.0 <= weak_subdu_tol_max
      && weak_subdu_tol_min <= weak_subdu_tol_max ) {
    const int           maxl =
      SC_MIN (maxlevel, discr_options->init_amr_weak_subdu_maxlevel);

    /* modify AMR parameters for refinement within distance */
    if (weak_subdu_indicator == SL_AMR_INDICATOR_WEAK_SUBDU_DIST) {
      weak_subdu_tol_min = -1.0;
      weak_subdu_tol_max = 0.1;
    }

    /* refine mesh */
    slabs_discr_init_amr_indicator (p8est, physics_options, discr_options,
                                    weak_subdu_indicator,
                                    weak_subdu_tol_min, weak_subdu_tol_max,
                                    0.0, 0.0,
                                    amr_max_steps, minlevel, maxl);
  }

  /* refine mesh w.r.t. ridge weak zone */
  if (   weak_subdu_indicator == SL_AMR_INDICATOR_WEAK_SUBDU_DIST
      && weak_ridge_indicator == SL_AMR_INDICATOR_WEAK_RIDGE_DIST
      && 0.0 <= weak_ridge_tol_max
      && weak_ridge_tol_min <= weak_ridge_tol_max ) {
    const int           maxl =
      SC_MIN (maxlevel, discr_options->init_amr_weak_subdu_maxlevel);

    /* modify AMR parameters for refinement within distance */
    weak_ridge_tol_min = -1.0;
    weak_ridge_tol_max = 0.1;

    /* refine mesh */
    slabs_discr_init_amr_indicator (p8est, physics_options, discr_options,
                                    weak_ridge_indicator,
                                    weak_ridge_tol_min, weak_ridge_tol_max,
                                    0.0, 0.0,
                                    amr_max_steps, minlevel, maxl);
  }

  /* refine mesh w.r.t. viscosity */
  if (   visc_indicator != SL_AMR_INDICATOR_NONE
      && visc_indicator <= SL_AMR_INDICATOR_RHS_OVER_VISC_GRAD
      && (
            (   (0.0 < visc_tol_min || 0.0 < visc_tol_max)
             && visc_tol_min < visc_tol_max )
          ||
            (   (0.0 < visc_in_plates_tol_min || 0.0 < visc_in_plates_tol_max)
             && visc_in_plates_tol_min < visc_in_plates_tol_max )
         )
     ) {
    slabs_discr_init_amr_indicator (p8est, physics_options, discr_options,
                                    visc_indicator,
                                    visc_tol_min, visc_tol_max,
                                    visc_in_plates_tol_min,
                                    visc_in_plates_tol_max,
                                    amr_max_steps, minlevel, maxlevel);
  }

  /* refine mesh w.r.t. right-hand side */
  if (   rhs_indicator != SL_AMR_INDICATOR_NONE
      && (   (0.0 < rhs_tol_min || 0.0 < rhs_tol_max)
          && rhs_tol_min < rhs_tol_max )
     ) {
    slabs_discr_init_amr_indicator (p8est, physics_options, discr_options,
                                    rhs_indicator,
                                    rhs_tol_min, rhs_tol_max,
                                    0.0, 0.0,
                                    amr_max_steps, minlevel, maxlevel);
  }

  /* uniform refinement after AMR */
  if (0 < post_uniform_n_steps) {
    slabs_discr_init_amr_indicator (p8est, physics_options, discr_options,
                                    SL_AMR_INDICATOR_REFINE_ALL,
                                    0.0, 0.5, 0.0, 0.0,
                                    post_uniform_n_steps, minlevel, maxlevel);
  }

  /* restore polynomial order of discretization */
  if (0 < amr_override_order && order != amr_override_order) {
    discr_options->order = order;
  }
}

/**
 * Gets temperature and velocity field of this element at discont. GLL nodes.
 */
static void
slabs_discr_amr_get_state_fields (sc_dmatrix_t *temp_el_mat,
                                  sc_dmatrix_t *weak_el_mat,
                                  sc_dmatrix_t *vel_el_mat,
                                  const mangll_locidx_t elid,
                                  slabs_stokes_state_t *state,
                                  mangll_t *mangll, mangll_cnodes_t *cnodes,
                                  slabs_physics_options_t *physics_options)
{
  const int           N = ymir_n (mangll->N);
  const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);
  const double       *x = mangll->X->e[elid];
  const double       *y = mangll->Y->e[elid];
  const double       *z = mangll->Z->e[elid];

  /* get temperature and velocity of this element */
  if (cnodes != NULL) {
    mangll_refel_t     *refel = mangll->refel;
    mangll_locidx_t    *el2cnode_map;
    int                 nodeid;
    int                 fieldid;

    /* get nodes indices that are relevant for transform: cont. -> discont. */
    el2cnode_map = cnodes->EtoCn + elid * n_nodes_per_el;

    /* get continuous temperature field */
    if (temp_el_mat != NULL) {
      sc_dmatrix_t       *temperature = state->temperature;

      YMIR_ASSERT (temperature != NULL);

      /* copy values from continuous nodes */
      for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
        temp_el_mat->e[nodeid][0] = temperature->e[el2cnode_map[nodeid]][0];
      }

      /* transform: continuous -> discontinuous */
      mangll_cdinterp_single (temp_el_mat->e[0], 1, cnodes->F[elid], refel, 0);
    }

    /* get continuous velocity field */
    if (vel_el_mat != NULL) {
      sc_dmatrix_t       *velocity = state->velocity;

      YMIR_ASSERT (velocity != NULL);

      /* copy values from continuous nodes */
      for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
        for (fieldid = 0; fieldid < 3; fieldid++) {
          vel_el_mat->e[nodeid][fieldid] =
            velocity->e[el2cnode_map[nodeid]][fieldid];
        }
      }

      /* transform: continuous -> discontinuous */
      mangll_cdinterp_single (vel_el_mat->e[0], 3, cnodes->F[elid], refel, 0);
    }
  }
  else {
    /* get discontinuous temperature field */
    if (temp_el_mat != NULL) {
      YMIR_ASSERT (state->amr_buffer_coarse_temp != NULL);

      memmove (temp_el_mat->e[0], state->amr_buffer_coarse_temp->e[elid],
               n_nodes_per_el * sizeof (double));
    }

    /* get discontinuous velocity field */
    if (vel_el_mat != NULL) {
      YMIR_ASSERT (state->amr_buffer_coarse_vel != NULL);

      memmove (vel_el_mat->e[0], state->amr_buffer_coarse_vel->e[elid],
               3 * n_nodes_per_el * sizeof (double));
    }
  }

  /* post-process temperature */
  slabs_temp_postprocess_elem (temp_el_mat, x, y, z, physics_options);

  /* get weak zone of this element */
  if (weak_el_mat != NULL) {
    YMIR_ASSERT (weak_el_mat->m == n_nodes_per_el && weak_el_mat->n == 1);

    if (state->weakzone != NULL) { /* if weak zone stored in state */
      double             *Brinv_data = mangll->refel->Brinv->e[0];
      sc_dmatrix_t       *weak_gauss_el_mat =
                            sc_dmatrix_new (n_nodes_per_el, 1);
      double             *weak_el_data = weak_el_mat->e[0];
      double             *weak_gauss_el_data = weak_gauss_el_mat->e[0];

      /* get weak zone on Gauss nodes from Stokes state */
      memmove (weak_gauss_el_data, state->weakzone->e[elid],
               n_nodes_per_el * sizeof (double));

      /* interpolate from Gauss to GLL nodes */
      MANGLL_TENSOR_IIAX_APPLY_ELEM (N + 1, Brinv_data,
                                     weak_gauss_el_data, weak_el_data);
      MANGLL_TENSOR_IAIX_APPLY_ELEM (N + 1, Brinv_data,
                                     weak_el_data, weak_gauss_el_data);
      MANGLL_TENSOR_AIIX_APPLY_ELEM (N + 1, Brinv_data,
                                     weak_gauss_el_data, weak_el_data);

      /* bound weak zone factor to valid range */
      slabs_matrix_bound_values (weak_el_mat, SC_EPS, 1.0);

      /* destroy */
      sc_dmatrix_destroy (weak_gauss_el_mat);
    }
    else { /* if weak zone has to be computed */
      /* compute weak zone of this element */
      slabs_weak_elem (weak_el_mat, x, y, z, mangll->refel->Vmask,
                       physics_options);
    }
  }
}

/**
 * Computes linear viscosity for an element.
 */
static void
slabs_discr_amr_visc_lin_elem (sc_dmatrix_t *temp_el_mat,
                               sc_dmatrix_t *weak_el_mat,
                               sc_dmatrix_t *visc_el_mat,
                               const mangll_locidx_t elid,
                               mangll_t *mangll, mangll_cnodes_t *cnodes,
                               slabs_stokes_state_t *state,
                               slabs_physics_options_t *physics_options)
{
  const int           visc_type = physics_options->viscosity_type;
  const int           visc_type_init =
                        physics_options->viscosity_type_for_init_nl_stokes;
  const double        visc_scaling = physics_options->viscosity_scaling;
  const double        upper_mantle_radius =
                        physics_options->viscosity_upper_mantle_radius;
  const double        visc_lower_mantle_scaling =
                        physics_options->viscosity_lower_mantle_scaling;
  const double       *_sc_restrict x = mangll->X->e[elid];
  const double       *_sc_restrict y = mangll->Y->e[elid];
  const double       *_sc_restrict z = mangll->Z->e[elid];

  /* check input */
  YMIR_ASSERT (mangll != NULL);
  YMIR_ASSERT (temp_el_mat != NULL);
  YMIR_ASSERT (weak_el_mat != NULL);
  YMIR_ASSERT (visc_el_mat != NULL);

  /* get temperature and weak zone fields of this element at
   * discontinuous GLL nodes */
  slabs_discr_amr_get_state_fields (temp_el_mat, weak_el_mat, NULL, elid,
                                    state, mangll, cnodes, physics_options);

  /* change scaling of viscosity in upper mantle
   * (snippet copied from slabs_physics.c) */
  if (   visc_type == SL_VISCOSITY_NONLINEAR
      && visc_type_init == SL_VISCOSITY_INIT_NL_STOKES_TEMP_UM_REL_TO_LM
      && 0.0 < upper_mantle_radius) {
    physics_options->viscosity_scaling = SL_VISCOSITY_UM_REL_TO_LM
                                         * visc_lower_mantle_scaling;
  }

  /* compute temperature dependent visosity for this element
   * (assume here: coordinates are used by `slabs_visc_temp_elem` to determine
   * position in upper/lower mantle only) */
  slabs_visc_temp_elem (visc_el_mat, x, y, z, mangll->refel->Vmask,
                        temp_el_mat, weak_el_mat, physics_options, 1);

  /* restore viscosity scaling (snippet copied from slabs_physics.c) */
  if (   visc_type == SL_VISCOSITY_NONLINEAR
      && visc_type_init == SL_VISCOSITY_INIT_NL_STOKES_TEMP_UM_REL_TO_LM
      && 0.0 < upper_mantle_radius) {
    physics_options->viscosity_scaling = visc_scaling;
  }
}

/**
 *
 */
static void
slabs_discr_amr_visc_nl_elem (sc_dmatrix_t *temp_el_mat,
                              sc_dmatrix_t *weak_el_mat,
                              sc_dmatrix_t *vel_el_mat,
                              sc_dmatrix_t *grad_vel_el_mat,
                              sc_dmatrix_t *IIe_el_mat,
                              sc_dmatrix_t *strain_rate_tensor_el_mat,
                              sc_dmatrix_t *visc_el_mat,
                              sc_dmatrix_t *dvisc_dIIe_el_mat,
                              sc_dmatrix_t *rank1_scal_el_mat,
                              sc_dmatrix_t *bounds_el_mat,
                              sc_dmatrix_t *yielding_el_mat,
                              const mangll_locidx_t elid,
                              mangll_t *mangll, mangll_cnodes_t *cnodes,
                              slabs_stokes_state_t *state,
                              slabs_physics_options_t *physics_options,
                              sc_dmatrix_t *tmp_dvel,
                              sc_dmatrix_t *tmp)
{
  /* check input */
  YMIR_ASSERT (mangll != NULL);

  /* get temperature, weak zone, and velocity fields of this element at
   * discontinuous GLL nodes */
  slabs_discr_amr_get_state_fields (temp_el_mat, weak_el_mat, vel_el_mat, elid,
                                    state, mangll, cnodes, physics_options);

  /* compute 2nd invariant of the strain rate at GLL nodes */
  if (IIe_el_mat != NULL) {
    YMIR_ASSERT (vel_el_mat != NULL);
    YMIR_ASSERT (grad_vel_el_mat != NULL);
    YMIR_ASSERT (tmp_dvel != NULL);
    YMIR_ASSERT (tmp != NULL);

    slabs_second_invariant_elem (vel_el_mat, IIe_el_mat, mangll, elid,
                                 grad_vel_el_mat, tmp_dvel, tmp,
                                 SL_GLL_DISCONTINUOUS_NODE);
  }

  /* compute strain rate tensor for this element */
  if (strain_rate_tensor_el_mat != NULL) {
    YMIR_ASSERT (grad_vel_el_mat != NULL);
    YMIR_ASSERT (IIe_el_mat != NULL);

    slabs_strain_rate_tensor_elem (strain_rate_tensor_el_mat, mangll,
                                   grad_vel_el_mat);
  }

  /* compute nonlinear visosity for this element */
  if (visc_el_mat != NULL) {
    const double       *_sc_restrict x = mangll->X->e[elid];
    const double       *_sc_restrict y = mangll->Y->e[elid];
    const double       *_sc_restrict z = mangll->Z->e[elid];

    YMIR_ASSERT (temp_el_mat != NULL);
    YMIR_ASSERT (IIe_el_mat != NULL);
    YMIR_ASSERT (weak_el_mat != NULL);
    YMIR_ASSERT (dvisc_dIIe_el_mat != NULL);

    /* compute viscosity (assume here: coordinates are used by
     * `slabs_visc_nl_elem` to determine position in upper/lower mantle only) */
    slabs_visc_nl_elem (visc_el_mat, dvisc_dIIe_el_mat, rank1_scal_el_mat,
                        bounds_el_mat, yielding_el_mat,
                        x, y, z, mangll->refel->Vmask,
                        temp_el_mat, weak_el_mat, IIe_el_mat,
                        physics_options, 0);
  }
}

/**
 * Computes derivative of 4th order viscosity tensor, (using index notation)
 *
 *   d_j Lambda_{i,j,k,l}
 *   =
 *   d_j visc * delta_{i,k} * delta_{j,l}
 *   +
 *   (d_j dvisc/dIIe) * e_{i,j} * e_{k,l}
 *   +
 *   dvisc/dIIe * ( (d_j e_{i,j}) * e_{k,l} + e_{i,j} * (d_j e_{k,l}) )
 *
 * where
 *
 *   delta_{i,k} --- Kronecker delta
 *   e_{i,j}     --- strain rate tensor
 *   IIe         --- 2nd invariant of the strain rate tensor
 */
static double
slabs_discr_amr_visc_tensor_deriv (const int i, const int k, const int l,
                                   const int nodeid,
                                   sc_dmatrix_t *grad_visc_el_mat,
                                   sc_dmatrix_t *dvisc_dIIe_el_mat,
                                   sc_dmatrix_t *grad_dvisc_dIIe_el_mat,
                                   sc_dmatrix_t *strain_rate_tensor_el_mat,
                                   sc_dmatrix_t *grad_strain_rate_el_mat)
{
  double             *_sc_restrict grad_visc_data =
                        grad_visc_el_mat->e[0] + 3 * nodeid;
  double             *_sc_restrict dvisc_dIIe_data =
                        dvisc_dIIe_el_mat->e[0] + nodeid;
  double             *_sc_restrict grad_dvisc_dIIe_data =
                        grad_dvisc_dIIe_el_mat->e[0] + 3 * nodeid;
  double             *_sc_restrict strain_rate_data =
                        strain_rate_tensor_el_mat->e[0] + 6 * nodeid;
  double             *_sc_restrict grad_strain_rate_data =
                        grad_strain_rate_el_mat->e[0] + 18 * nodeid;
  int                 j;
  double              sum, deriv;

  YMIR_ASSERT (0 <= i && i < 3);
  YMIR_ASSERT (0 <= k && k < 3);
  YMIR_ASSERT (0 <= l && l < 3);

  /* set 1st summand: `d_l visc * delta_{i,k}` */
  if (i == k) {
    deriv = grad_visc_data[l];
  }
  else {
    deriv = 0.0;
  }

  /* add 2nd summand: `(d_j dvisc/dIIe) * e_{i,j} * e_{k,l}` */
  sum = 0.0;
  for (j = 0; j < 3; j++) {
    sum += grad_dvisc_dIIe_data[j]
           * strain_rate_data[slabs_strain_rate_tensor_idx (i, j)];
  }
  deriv += sum * strain_rate_data[slabs_strain_rate_tensor_idx (k, l)];

  /* add 3rd summand:
   *   `dvisc/dIIe * ( (d_j e_{i,j}) * e_{k,l} + e_{i,j} * (d_j e_{k,l}) )` */
  sum = 0.0;
  for (j = 0; j < 3; j++) {
    sum +=   grad_strain_rate_data[3 * slabs_strain_rate_tensor_idx (i, j) + j]
           * strain_rate_data[slabs_strain_rate_tensor_idx (k, l)]
           + strain_rate_data[slabs_strain_rate_tensor_idx (i, j)]
           * grad_strain_rate_data[3 * slabs_strain_rate_tensor_idx (k, l) + j];
  }
  deriv += dvisc_dIIe_data[0] * sum;

  /* return result */
  return deriv;
}

/**
 * Computes right-hand side for an element.
 */
static void
slabs_discr_amr_rhs_elem (sc_dmatrix_t *temp_el_mat,
                          sc_dmatrix_t *rhs_el_mat,
                          const mangll_locidx_t elid,
                          mangll_t *mangll, mangll_cnodes_t *cnodes,
                          slabs_stokes_state_t *state,
                          slabs_physics_options_t *physics_options)
{
  const double       *x = mangll->X->e[elid];
  const double       *y = mangll->Y->e[elid];
  const double       *z = mangll->Z->e[elid];

  /* check input */
  YMIR_ASSERT (mangll != NULL);
  YMIR_ASSERT (temp_el_mat != NULL);
  YMIR_ASSERT (rhs_el_mat != NULL);

  /* get temperature field of this element at discontinuous GLL nodes */
  slabs_discr_amr_get_state_fields (temp_el_mat, NULL, NULL, elid,
                                    state, mangll, cnodes, physics_options);

  /* compute right-hand side for this element */
  slabs_rhs_elem (rhs_el_mat, x, y, z, temp_el_mat, physics_options);
}

/**
 * Sets values of AMR indicator.
 */
void
slabs_discr_amr_indicator_set (slabs_discr_amr_indicator_t *indicator,
                               slabs_stokes_state_t *state,
                               mangll_t *mangll, mangll_cnodes_t *cnodes,
                               slabs_physics_options_t *physics_options,
                               slabs_discr_options_t *discr_options,
                               const int init_amr)
{
#ifdef YMIR_DEBUG
  const char         *this_fn_name = "slabs_discr_amr_indicator_set";
#endif

  const slabs_discr_amr_indicator_type_t  indicator_type = indicator->type;
  double             *indicator_data = indicator->val->e[0];
  const double        domain_norm = discr_options->domain_size_normalization;
  int                 amr_lower_mantle;
  const int          *Vmask = mangll->refel->Vmask;
  const mangll_locidx_t  n_elements = mangll->mesh->K;
  const int           N = ymir_n (mangll->N);
  const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);

  mangll_locidx_t     elid;
  int                 nodeid, fieldid;
  double              min, max;

  sc_dmatrix_t       *temp_el_mat = NULL, *vel_el_mat = NULL;
  sc_dmatrix_t       *grad_vel_el_mat = NULL, *tmp_dvel = NULL, *tmp = NULL;
  sc_dmatrix_t       *IIe_el_mat = NULL, *strain_rate_tensor_el_mat = NULL;
  sc_dmatrix_t       *weak_el_mat = NULL;
  sc_dmatrix_t       *visc_el_mat = NULL, *dvisc_dIIe_el_mat = NULL;
  sc_dmatrix_t       *yielding_el_mat = NULL;

  sc_dmatrix_t       *grad_visc_el_mat = NULL;

  sc_dmatrix_t       *grad_dvisc_dIIe_el_mat = NULL;
  sc_dmatrix_t       *grad_strain_rate_tens_el_mat = NULL;
  sc_dmatrix_t       *tmp_strain_rate_tens = NULL;

  sc_dmatrix_t       *strain_rate_el_mat = NULL;
  sc_dmatrix_t       *grad_strain_rate_el_mat = NULL;

  sc_dmatrix_t       *rhs_el_mat = NULL;
  sc_dmatrix_t       *grad_rhs_el_mat = NULL;
  sc_dmatrix_t       *tmp_rhs = NULL;
#ifdef YMIR_DEBUG
  double              indicator_loc_min = DBL_MAX;
  double              indicator_loc_max = 0.0;
  double              indicator_loc_sum = 0.0;
#endif

  /* set indicator paramters */
  if (init_amr) {
    amr_lower_mantle = discr_options->init_amr_lower_mantle;
  }
  else {
    amr_lower_mantle = discr_options->amr_lower_mantle;
  }

  /* return zero indicator if nothing to do */
  if ( indicator_type == SL_AMR_INDICATOR_NONE ||
       indicator_type == SL_AMR_INDICATOR_COARSEN_ALL ) {
    sc_dmatrix_set_zero (indicator->val);
    return;
  }

  /* return one indicator if refine all */
  if (indicator_type == SL_AMR_INDICATOR_REFINE_ALL) {
    sc_dmatrix_set_value (indicator->val, 1.0);
    return;
  }

  /* create work variables */
  if ( indicator_type == SL_AMR_INDICATOR_VISC_DR ||
       indicator_type == SL_AMR_INDICATOR_VISC_GRAD ||
       indicator_type == SL_AMR_INDICATOR_VISC_PECLET ||
       indicator_type == SL_AMR_INDICATOR_VISC_TENSOR_PECLET ) {
    temp_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
    weak_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);

    if (!init_amr) {
      vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
      grad_vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 9);
      IIe_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
    }

    visc_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
    if (!init_amr) {
      dvisc_dIIe_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
    }

    if (!init_amr) {
      tmp_dvel = sc_dmatrix_new (n_nodes_per_el, 3);
    }
    tmp = sc_dmatrix_new (n_nodes_per_el, 1);
  }
  if ( indicator_type == SL_AMR_INDICATOR_VISC_GRAD ||
       indicator_type == SL_AMR_INDICATOR_VISC_PECLET ) {
    grad_visc_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
  }
  if (indicator_type == SL_AMR_INDICATOR_VISC_TENSOR_PECLET) {
    YMIR_ASSERT (!init_amr);

    yielding_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
    grad_visc_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
    grad_dvisc_dIIe_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
    strain_rate_tensor_el_mat = sc_dmatrix_new (n_nodes_per_el, 6);
    grad_strain_rate_tens_el_mat = sc_dmatrix_new (n_nodes_per_el, 18);
    tmp_strain_rate_tens = sc_dmatrix_new (n_nodes_per_el, 6);
  }

  if (indicator_type == SL_AMR_INDICATOR_STRAIN_RATE_PECLET) {
    YMIR_ASSERT (!init_amr);

    temp_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
    vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
    grad_vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 9);
    strain_rate_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);

    tmp_dvel = sc_dmatrix_new (n_nodes_per_el, 3);
    tmp = sc_dmatrix_new (n_nodes_per_el, 1);

    grad_strain_rate_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
  }

  if ( indicator_type == SL_AMR_INDICATOR_WEAK_IMPORT_DIST &&
       physics_options->weakzone_type != SL_WEAKZONE_IMPORT_FILE ) {
    weak_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  }

  if (indicator_type == SL_AMR_INDICATOR_RHS_MAGNITUDE) {
    temp_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
    rhs_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
  }

  if (indicator_type == SL_AMR_INDICATOR_RHS_PECLET) {
    temp_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
    rhs_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
    grad_rhs_el_mat = sc_dmatrix_new (n_nodes_per_el, 9);
    tmp_rhs = sc_dmatrix_new (n_nodes_per_el, 3);
  }

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    const double       *_sc_restrict x = mangll->X->e[elid];
    const double       *_sc_restrict y = mangll->Y->e[elid];
    const double       *_sc_restrict z = mangll->Z->e[elid];

    /* if no AMR in lower mantle */
    if (!amr_lower_mantle && indicator_type != SL_AMR_INDICATOR_REFINE_LAYER &&
        !slabs_physics_elem_in_upper_mantle (x, y, z, Vmask, physics_options)) {
      indicator_data[elid] = 0.0;
      continue;
    }

    /* compute indicator */
    switch (indicator_type) {
    case SL_AMR_INDICATOR_VISC_DR: /* dynamic range of viscosity */
      /* compute visosity for this element */
      if (init_amr) {
        slabs_discr_amr_visc_lin_elem (temp_el_mat, weak_el_mat, visc_el_mat,
                                       elid, mangll, cnodes, state,
                                       physics_options);
      }
      else {
        slabs_discr_amr_visc_nl_elem (temp_el_mat, weak_el_mat, vel_el_mat,
            grad_vel_el_mat, IIe_el_mat, NULL,
            visc_el_mat, dvisc_dIIe_el_mat, NULL, NULL, NULL,
            elid, mangll, cnodes, state, physics_options, tmp_dvel, tmp);
      }

      /* find min and max viscosity */
      slabs_matrix_compute_abs_min_max (&min, &max, visc_el_mat);

      /* set indicator to dynamic range for this element */
      indicator_data[elid] = max / min;
      break;

    case SL_AMR_INDICATOR_VISC_GRAD: /* max norm of viscosity gradient */
      /* compute visosity for this element */
      if (init_amr) {
        slabs_discr_amr_visc_lin_elem (temp_el_mat, weak_el_mat, visc_el_mat,
                                       elid, mangll, cnodes, state,
                                       physics_options);
      }
      else {
        slabs_discr_amr_visc_nl_elem (temp_el_mat, weak_el_mat, vel_el_mat,
            grad_vel_el_mat, IIe_el_mat, NULL,
            visc_el_mat, dvisc_dIIe_el_mat, NULL, NULL, NULL,
            elid, mangll, cnodes, state, physics_options, tmp_dvel, tmp);
      }

      /* compute gradient */
      slabs_gradient_gll_to_gll_elem (visc_el_mat, grad_visc_el_mat, mangll,
                                      elid, tmp);

      /* calculate max norm of gradient */
      max = slabs_matrix_compute_abs_max (grad_visc_el_mat);

      /* set indicator to max weighted with element size */
      indicator_data[elid] = max * slabs_elem_volume (mangll, elid)
                                 * domain_norm;
      break;

    case SL_AMR_INDICATOR_VISC_PECLET: /* grad viscosity / viscosity */
      /* compute visosity for this element */
      if (init_amr) {
        slabs_discr_amr_visc_lin_elem (temp_el_mat, weak_el_mat, visc_el_mat,
                                       elid, mangll, cnodes, state,
                                       physics_options);
      }
      else {
        slabs_discr_amr_visc_nl_elem (temp_el_mat, weak_el_mat, vel_el_mat,
            grad_vel_el_mat, IIe_el_mat, NULL,
            visc_el_mat, dvisc_dIIe_el_mat, NULL, NULL, NULL,
            elid, mangll, cnodes, state, physics_options, tmp_dvel, tmp);
      }

      /* compute gradient */
      slabs_gradient_gll_to_gll_elem (visc_el_mat, grad_visc_el_mat, mangll,
                                      elid, tmp);

      /* compute max indicator for this element */
      max = 0.0;
      for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
        const double        visc = visc_el_mat->e[nodeid][0];
        const double       *grad = grad_visc_el_mat->e[nodeid];
        double              grad_norm;

        /* compute Frobenius norm of gradient matrix at this node */
        grad_norm = 0.0;
        for (fieldid = 0; fieldid < 3; fieldid++) {
          grad_norm += grad[fieldid] * grad[fieldid];
        }
        grad_norm = sqrt (grad_norm);

        /* update max */
        max = SC_MAX (max, grad_norm / visc);
      }

      /* set indicator to max weighted with element size */
      indicator_data[elid] = max * slabs_elem_volume (mangll, elid)
                                 * domain_norm;
      break;

    case SL_AMR_INDICATOR_VISC_TENSOR_PECLET: /* grad visc tens / visc tens */
      {
        const double        stress_exp =
                              physics_options->viscosity_stress_exponent;
        int                 i, k, l;
        double              entry, val;

        /* compute visosity, its derivative w.r.t. 2nd invariant,
         * and strain rate tensor */
        slabs_discr_amr_visc_nl_elem (temp_el_mat, weak_el_mat, vel_el_mat,
            grad_vel_el_mat, IIe_el_mat, strain_rate_tensor_el_mat,
            visc_el_mat, dvisc_dIIe_el_mat, NULL, NULL, yielding_el_mat,
            elid, mangll, cnodes, state, physics_options, tmp_dvel, tmp);

        /* compute gradient of viscosity */
        slabs_gradient_gll_to_gll_elem (visc_el_mat, grad_visc_el_mat, mangll,
                                        elid, tmp);

        /* compute gradient of viscosity derivative */
        slabs_gradient_gll_to_gll_elem (dvisc_dIIe_el_mat,
                                        grad_dvisc_dIIe_el_mat, mangll, elid,
                                        tmp);

        /* compute gradient of strain rate tensor */
        slabs_gradient_gll_to_gll_elem (strain_rate_tensor_el_mat,
                                        grad_strain_rate_tens_el_mat, mangll,
                                        elid, tmp_strain_rate_tens);

        max = 0.0;
        for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
          /* compute (scaled) Frobenius norm of gradient of viscosity tensor */
          val = 0.0;
          for (k = 0; k < 3; k++) {
            for (i = 0; i < 3; i++) {
              entry = slabs_discr_amr_visc_tensor_deriv (i, k, k, nodeid,
                  grad_visc_el_mat, dvisc_dIIe_el_mat, grad_dvisc_dIIe_el_mat,
                  strain_rate_tensor_el_mat, grad_strain_rate_tens_el_mat);
              val += entry * entry;
            }
          }
          for (k = 0; k < 3; k++) {
            for (l = k + 1; l < 3; l++) {
              for (i = 0; i < 3; i++) {
                entry = slabs_discr_amr_visc_tensor_deriv (i, k, l, nodeid,
                    grad_visc_el_mat, dvisc_dIIe_el_mat, grad_dvisc_dIIe_el_mat,
                    strain_rate_tensor_el_mat, grad_strain_rate_tens_el_mat);
                val += 2.0 * entry * entry;
              }
            }
          }
          val = sqrt (val) / 3.0;

          /* divide by the (scaled) Frobenius norm of the viscosity tensor */
          if (yielding_el_mat->e[nodeid][0]) {
          val /= visc_el_mat->e[nodeid][0]
                 * sqrt ( 8.0 + SC_SQR (1.0 + 4.0 *
                                        (1.0 - stress_exp)/(2.0*stress_exp)) )
                 / 3.0;
          }
          else {
            val /= visc_el_mat->e[nodeid][0]; /* `* 3 / 3` cancels */
          }

          /* update max */
          max = SC_MAX (max, val);
        }

        /* set indicator to max over all nodes weighted with element size */
        indicator_data[elid] = max * slabs_elem_volume (mangll, elid)
                                   * domain_norm;
      }
      break;

    case SL_AMR_INDICATOR_STRAIN_RATE_PECLET: /* grad strain rate / strain r. */
      /* compute strain rate (i.e., sqrt of 2nd invariant) for this element */
      slabs_discr_amr_visc_nl_elem (temp_el_mat, NULL, vel_el_mat,
          grad_vel_el_mat, strain_rate_el_mat, NULL,
          NULL, NULL, NULL, NULL, NULL,
          elid, mangll, cnodes, state, physics_options, tmp_dvel, tmp);
      sc_dmatrix_sqrt (strain_rate_el_mat, strain_rate_el_mat);

      /* compute gradient (overwrite `temp_el_mat`) */
      slabs_gradient_gll_to_gll_elem (strain_rate_el_mat,
                                      grad_strain_rate_el_mat, mangll, elid,
                                      tmp);

      /* compute Euclidean norm of gradient and divide by strain rate;
       * and find max of this element */
      max = 0.0;
      for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
        const double        strain = strain_rate_el_mat->e[nodeid][0];

        if (SC_1000_EPS < strain) {
          const double       *grad = grad_strain_rate_el_mat->e[nodeid];
          double              grad_norm;

          /* compute Frobenius norm of gradient matrix at this node */
          grad_norm = 0.0;
          for (fieldid = 0; fieldid < 3; fieldid++) {
            grad_norm += grad[fieldid] * grad[fieldid];
          }
          grad_norm = sqrt (grad_norm);

          /* update max */
          max = SC_MAX (max, grad_norm / strain);
        }
      }

      /* set indicator to max weighted with element size */
      indicator_data[elid] = max * slabs_elem_volume (mangll, elid)
                                 * domain_norm;
      break;

    case SL_AMR_INDICATOR_WEAK_IMPORT_DIST: /* distance to weak zone */
      if (physics_options->weakzone_type == SL_WEAKZONE_IMPORT_FILE) {
        const double        weak_maxdist =
                              discr_options->init_amr_weak_import_maxdist;
        WeakzoneKDTree     *tree = slabs_physics_weakzone_kdtree;
        double              pt[3];
        double              dist;
        int                 i;

        YMIR_ASSERT (tree != NULL);

        /* if element is not under the influence of weak zones */
        if ( !slabs_physics_elem_has_weakzone (x, y, z, Vmask,
                                               physics_options) ) {
          indicator_data[elid] = 0.0;
          continue;
        }

        /* find min distance from a vertex point to point cloud */
        indicator_data[elid] = 0.0;
        for (i = 0; i < 8; i++) { /* loop over all vertices */
          const int           nodeid = Vmask[i];

          /* set target point */
          pt[0] = x[nodeid];
          pt[1] = y[nodeid];
          pt[2] = z[nodeid];

          /* compute distance to weak zone via kd-tree nearest neighbor
           * search on point cloud */
          dist = WeakzoneKDTree_find_shortest_distance_single (tree, pt);

          /* set indicator if distance below tolerance,
           * stop distance computation */
          if (dist < weak_maxdist) {
            indicator_data[elid] = 1.0;
            break;
          }
        }
      }
      else { /* if kd-tree does not exist */
        /* get weak zone of this element at discontinuous GLL nodes */
        slabs_discr_amr_get_state_fields (NULL, weak_el_mat, NULL, elid,
                                          state, mangll, cnodes,
                                          physics_options);

        /* check weak zone factor */
        indicator_data[elid] = 0.0;
        for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
          /* set indicator if weak zone factor below tolerance, stop loop */
          if (weak_el_mat->e[nodeid][0] < 0.5) {
            indicator_data[elid] = 1.0;
            break;
          }
        }
      }
      break;

    case SL_AMR_INDICATOR_RHS_MAGNITUDE: /* magnitude of rhs */
      {
        const double        rhs_scaling = fabs (physics_options->rhs_scaling);

        /* compute right-hand side for this element */
        slabs_discr_amr_rhs_elem (temp_el_mat, rhs_el_mat, elid, mangll, cnodes,
                                  state, physics_options);

        /* compute max indicator for this element */
        max = 0.0;
        for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
          const double       *rhs = rhs_el_mat->e[nodeid];
          double              rhs_norm;

          /* compute Euclidian norm of right-hand side at this node */
          rhs_norm = sqrt (rhs[0]*rhs[0] + rhs[1]*rhs[1] + rhs[2]*rhs[2]);

          /* normalize right-hand side by scaling value; now assume
           * 0 <= rhs_norm <= 1 */
          rhs_norm /= rhs_scaling;

          /* update max */
          max = SC_MAX (max, rhs_norm);
        }

        /* set max norm weighted with element size */
        indicator_data[elid] = max * slabs_elem_volume (mangll, elid)
                                   * domain_norm;
      }
      break;

    case SL_AMR_INDICATOR_RHS_PECLET: /* grad rhs / rhs */
      {
        const double        rhs_norm_shift =
                              discr_options->init_amr_rhs_norm_shift;

        /* compute right-hand side for this element */
        slabs_discr_amr_rhs_elem (temp_el_mat, rhs_el_mat, elid, mangll, cnodes,
                                  state, physics_options);

        /* compute gradient of right-hand side */
        slabs_gradient_gll_to_gll_elem (rhs_el_mat, grad_rhs_el_mat, mangll,
                                        elid, tmp_rhs);

        /* compute max indicator for this element */
        max = 0.0;
        for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
          const double       *rhs = rhs_el_mat->e[nodeid];
          double              rhs_norm;

          /* compute Euclidian norm of right-hand side at this node */
          rhs_norm = sqrt (rhs[0]*rhs[0] + rhs[1]*rhs[1] + rhs[2]*rhs[2]);

          if (SC_1000_EPS < rhs_norm) {
            const double       *grad = grad_rhs_el_mat->e[nodeid];
            double              grad_norm;

            /* shift norm of right-hand side */
            if (0.0 < rhs_norm_shift) {
              rhs_norm += rhs_norm_shift;
            }

            /* compute Frobenius norm of gradient matrix at this node */
            grad_norm = 0.0;
            for (fieldid = 0; fieldid < 9; fieldid++) {
              grad_norm += grad[fieldid] * grad[fieldid];
            }
            grad_norm = sqrt (grad_norm);

            /* update max */
            max = SC_MAX (max, grad_norm / rhs_norm);
          }
        }

        /* set indicator to max weighted with element size */
        indicator_data[elid] = max * slabs_elem_volume (mangll, elid)
                                   * domain_norm;
      }
      break;

    case SL_AMR_INDICATOR_REFINE_SURFACE: /* distance from surface */
      {
        const double        surface_radius = 1.0;
        const double        maxdist =
                              discr_options->refine_surface_maxdist *
                              (SL_SHELL_RADIUS_TOP - SL_SHELL_RADIUS_BOTTOM);
        double              r_first, r_last;

        YMIR_ASSERT (0.0 < maxdist);

        /* compute radii */
        r_first = slabs_compute_radius (x[0], y[0], z[0], physics_options);
        r_last = slabs_compute_radius (x[n_nodes_per_el-1],
                                       y[n_nodes_per_el-1],
                                       z[n_nodes_per_el-1], physics_options);

        /* if either node at opposite element corners are within max dist */
        if (fabs (surface_radius - r_first) <= maxdist ||
            fabs (surface_radius - r_last) <= maxdist) {
          indicator_data[elid] = 1.0;
        }
        else {
          indicator_data[elid] = 0.0;
        }
      }
      break;

    case SL_AMR_INDICATOR_REFINE_LAYER: /* distance layer w/ constant radius */
      {
        const double        upper_mantle_radius =
                              physics_options->viscosity_upper_mantle_radius;
        const double        maxdist =
                              discr_options->refine_layer_maxdist *
                              (SL_SHELL_RADIUS_TOP - SL_SHELL_RADIUS_BOTTOM);
        double              r_first, r_last;

        YMIR_ASSERT (0.0 < upper_mantle_radius);
        YMIR_ASSERT (0.0 < maxdist);

        /* compute radii */
        r_first = slabs_compute_radius (x[0], y[0], z[0], physics_options);
        r_last = slabs_compute_radius (x[n_nodes_per_el-1],
                                       y[n_nodes_per_el-1],
                                       z[n_nodes_per_el-1], physics_options);

        /* if either node at opposite element corners are within max dist */
        if (fabs (upper_mantle_radius - r_first) <= maxdist ||
            fabs (upper_mantle_radius - r_last) <= maxdist) {
          indicator_data[elid] = 1.0;
        }
        else {
          indicator_data[elid] = 0.0;
        }
      }
      break;

    default: /* unknown AMR indicator type */
      YMIR_ABORT_NOT_REACHED ();
    }

    /* update minimum and maximum indicator value */
#ifdef YMIR_DEBUG
    indicator_loc_min = SC_MIN (indicator_loc_min, indicator_data[elid]);
    indicator_loc_max = SC_MAX (indicator_loc_max, indicator_data[elid]);
    indicator_loc_sum += indicator_data[elid];
#endif
  }

  /* destroy work variables */
  if ( indicator_type == SL_AMR_INDICATOR_VISC_DR ||
       indicator_type == SL_AMR_INDICATOR_VISC_GRAD ||
       indicator_type == SL_AMR_INDICATOR_VISC_PECLET ||
       indicator_type == SL_AMR_INDICATOR_VISC_TENSOR_PECLET ) {
    sc_dmatrix_destroy (temp_el_mat);
    sc_dmatrix_destroy (weak_el_mat);

    if (!init_amr) {
      sc_dmatrix_destroy (vel_el_mat);
      sc_dmatrix_destroy (grad_vel_el_mat);
      sc_dmatrix_destroy (IIe_el_mat);
    }

    sc_dmatrix_destroy (visc_el_mat);
    if (!init_amr) {
      sc_dmatrix_destroy (dvisc_dIIe_el_mat);
    }

    if (!init_amr) {
      sc_dmatrix_destroy (tmp_dvel);
    }
    sc_dmatrix_destroy (tmp);
  }
  if ( indicator_type == SL_AMR_INDICATOR_VISC_GRAD ||
       indicator_type == SL_AMR_INDICATOR_VISC_PECLET ) {
    sc_dmatrix_destroy (grad_visc_el_mat);
  }
  if (indicator_type == SL_AMR_INDICATOR_VISC_TENSOR_PECLET ) {
    sc_dmatrix_destroy (yielding_el_mat);
    sc_dmatrix_destroy (grad_visc_el_mat);
    sc_dmatrix_destroy (grad_dvisc_dIIe_el_mat);
    sc_dmatrix_destroy (strain_rate_tensor_el_mat);
    sc_dmatrix_destroy (grad_strain_rate_tens_el_mat);
    sc_dmatrix_destroy (tmp_strain_rate_tens);
  }

  if (indicator_type == SL_AMR_INDICATOR_STRAIN_RATE_PECLET) {
    sc_dmatrix_destroy (temp_el_mat);
    sc_dmatrix_destroy (vel_el_mat);
    sc_dmatrix_destroy (grad_vel_el_mat);
    sc_dmatrix_destroy (strain_rate_el_mat);

    sc_dmatrix_destroy (tmp_dvel);
    sc_dmatrix_destroy (tmp);

    sc_dmatrix_destroy (grad_strain_rate_el_mat);
  }

  if ( indicator_type == SL_AMR_INDICATOR_WEAK_IMPORT_DIST &&
       physics_options->weakzone_type != SL_WEAKZONE_IMPORT_FILE ) {
    sc_dmatrix_destroy (weak_el_mat);
  }

  if (indicator_type == SL_AMR_INDICATOR_RHS_MAGNITUDE) {
    sc_dmatrix_destroy (temp_el_mat);
    sc_dmatrix_destroy (rhs_el_mat);
  }

  if (indicator_type == SL_AMR_INDICATOR_RHS_PECLET) {
    sc_dmatrix_destroy (temp_el_mat);

    sc_dmatrix_destroy (rhs_el_mat);
    sc_dmatrix_destroy (grad_rhs_el_mat);
    sc_dmatrix_destroy (tmp_rhs);
  }

  /* communicate global minimum and maximum indicator value */
#ifdef YMIR_DEBUG
  {
    MPI_Comm            mpicomm = state->p8est->mpicomm;
    int                 mpiret;
    const double        n_quadrants =
                          (double) state->p8est->global_num_quadrants;
    double              indicator_glo_min = 0.0;
    double              indicator_glo_max = 0.0;
    double              indicator_glo_sum = 0.0;

    mpiret = MPI_Allreduce (&indicator_loc_min, &indicator_glo_min, 1,
                            MPI_DOUBLE, MPI_MIN, mpicomm);
    YMIR_CHECK_MPI (mpiret);

    mpiret = MPI_Allreduce (&indicator_loc_max, &indicator_glo_max, 1,
                            MPI_DOUBLE, MPI_MAX, mpicomm);
    YMIR_CHECK_MPI (mpiret);

    mpiret = MPI_Allreduce (&indicator_loc_sum, &indicator_glo_sum, 1,
                            MPI_DOUBLE, MPI_SUM, mpicomm);
    YMIR_CHECK_MPI (mpiret);

    YMIR_GLOBAL_INFOF ("%s: Indicator %i, min %1.3e, max %1.3e, mean %1.3e\n",
                       this_fn_name, indicator_type, indicator_glo_min,
                       indicator_glo_max, indicator_glo_sum / n_quadrants);
  }
#endif
}

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
                 const int init_amr)
{
  const char         *this_fn_name = "slabs_discr_amr";

  const double        n_children = (double) P4EST_CHILDREN;
  const double        coarsen_factor = 1.0 / n_children - 1.0;
  const double        refine_factor = n_children - 1.0;
  const int           n_indicators = indicator_params->n_indicators;
  const int           uniform =
    ( n_indicators == 1 &&
      indicator_params->type[0] == SL_AMR_INDICATOR_REFINE_ALL );
  int                 amr_max_steps;
  double              amr_rel_threshold;
  p4est_gloidx_t      amr_n_quadrants_max;

  mangll_t           *mangll_coarse = (*mesh)->ma;
  mangll_t           *mangll_fine = NULL;
  mangll_t           *mangll_partitioned = NULL;
  mangll_cnodes_t    *cnodes = (*mesh)->cnodes;
  int                 amr_step;
#ifdef SLABS_DISCR_AMR_VTK
  char                filename[BUFSIZ];
#endif

  YMIR_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  /* set AMR parameters */
  if (init_amr) {
    amr_max_steps = discr_options->init_amr_max_steps;
    amr_rel_threshold = SC_MAX (0.0, discr_options->init_amr_rel_threshold);
    amr_n_quadrants_max =
      (p4est_gloidx_t) discr_options->init_amr_n_elements_max;
  }
  else {
    amr_max_steps = discr_options->amr_max_steps;
    amr_rel_threshold = SC_MAX (0.0, discr_options->amr_rel_threshold);
    amr_n_quadrants_max = (p4est_gloidx_t) discr_options->amr_n_elements_max;
  }

  /* check number of max AMR steps */
  if (amr_max_steps <= 0) {
    amr_max_steps = P4EST_MAXLEVEL;
  }

  for (amr_step = 0; amr_step < amr_max_steps; amr_step++) {
    const p4est_gloidx_t  n_quadrants = p8est->global_num_quadrants;
    p4est_gloidx_t      n_coarsen_quadrants_glo = -1;
    p4est_gloidx_t      n_refine_quadrants_glo = -1;
    slabs_discr_amr_marker_t *marker;
    int                 indicator_skipped = 0;
    int                 indicator_id;

    YMIR_GLOBAL_INFOF ("%s: Into AMR step %i\n", this_fn_name, amr_step);

    /* check variables */
    YMIR_ASSERT (mangll_coarse != NULL);

    /* print mesh statistics */
    ymir_monitor_print_global_element_stats (p8est);

    /* print memory usage */
    //ymir_monitor_print_global_mem_usage (p8est->mpicomm);

    /*
     * Compute AMR Indicator and Mark Elements for Refinement or Coarsening
     */

    /* create new AMR marker */
    marker = slabs_discr_amr_marker_new (p8est);

    /* mark elements */
    for (indicator_id = 0; indicator_id < n_indicators; indicator_id++) {
      const slabs_discr_amr_indicator_type_t  indicator_type =
        indicator_params->type[indicator_id];
      const double        tol_min = indicator_params->tol_min[indicator_id];
      const double        tol_max = indicator_params->tol_max[indicator_id];
      int8_t              level_min, level_max;
      slabs_discr_amr_indicator_t *indicator;
      slabs_discr_amr_marker_t *marker_next;
      p4est_gloidx_t      n_quadrants_next;

      /* skip if indicator type is not ok or tolerances are invalid */
      //if ( indicator_type == SL_AMR_INDICATOR_NONE ||
      //     tol_min < 0.0 || tol_max <= 0.0 )
      // TODO interpolation for coarsening is not implemented, instead:
      if ( indicator_type == SL_AMR_INDICATOR_NONE ||
           tol_min > 0.0 || tol_max <= 0.0 ) {
        YMIR_GLOBAL_INFOF (
            "%s: Step %i, skipped indicator %i, tol [%g,%g]\n",
            this_fn_name, amr_step, (int) indicator_type, tol_min, tol_max);

        indicator_skipped++;
        continue;
      }

      /* set min level */
      level_min = indicator_params->level_min[indicator_id];
      if (level_min <= 0 && 0 < discr_options->minlevel) {
        level_min = discr_options->minlevel;
      }
      else if (0 < level_min && 0 < discr_options->minlevel) {
        level_min = SC_MAX (level_min, discr_options->minlevel);
      }

      /* set max level */
      level_max = indicator_params->level_max[indicator_id];
      if (level_max <= 0 && 0 < discr_options->maxlevel) {
        level_max = discr_options->maxlevel;
      }
      else if (0 < level_max && 0 < discr_options->maxlevel) {
        level_max = SC_MIN (level_max, discr_options->maxlevel);
      }

      YMIR_GLOBAL_INFOF (
          "%s: Step %i, indicator %i, tol [%g,%g], level bounds [%i,%i]\n",
          this_fn_name, amr_step, (int) indicator_type,
          tol_min, tol_max, level_min, level_max);

      /* create new AMR indicator */
      indicator = slabs_discr_amr_indicator_new (
          indicator_type, tol_min, tol_max, level_min, level_max,
          p8est->local_num_quadrants);

      /* set the current AMR indicator */
      slabs_discr_amr_indicator_set (indicator, state, mangll_coarse, cnodes,
                                     physics_options, discr_options, init_amr);

      /* update a copy of AMR marker with the current AMR indicator */
      marker_next = slabs_discr_amr_marker_duplicate (marker);
      if (indicator_id <= indicator_skipped) { /* if first valid indicator */
        slabs_discr_amr_marker_set_indicator (marker_next, indicator);
      }
      else { /* if a subsequent valid indicator */
        slabs_discr_amr_marker_merge_indicator (marker_next, indicator);
      }

      if (!uniform) {
        /* estimate the new number of quadrants with current AMR marker */
        slabs_discr_amr_marker_get_quadrant_counts (
            &n_coarsen_quadrants_glo, &n_refine_quadrants_glo,
            &n_quadrants_next, marker_next);

        YMIR_GLOBAL_INFOF (
            "%s: Step %i, indicator %i: #quadrants marked for "
            "coarsening %lli, refinement %lli\n",
            this_fn_name, amr_step, (int) indicator_type,
            (long long int) n_coarsen_quadrants_glo,
            (long long int) n_refine_quadrants_glo);
        YMIR_GLOBAL_INFOF (
            "%s: Step %i, indicator %i: #quadrants next %lli (approx)\n",
            this_fn_name, amr_step, (int) indicator_type,
            (long long int) n_quadrants_next);
      }

      /* check if the next number of quadrants (appox) is above max limit */
      if ( !uniform && 0 < amr_n_quadrants_max &&
           amr_n_quadrants_max < n_quadrants_next ) {
        const double        diff = (double) (amr_n_quadrants_max - n_quadrants);
        const double        n_coarsen = (double) n_coarsen_quadrants_glo;
        p4est_gloidx_t      n_refine_quadrants_max;

        /* set max #quadrants that are allowed to be refined */
        n_refine_quadrants_max = (p4est_gloidx_t) round (diff / refine_factor);
        n_refine_quadrants_max -= (p4est_gloidx_t) coarsen_factor * n_coarsen;

        YMIR_GLOBAL_INFOF (
            "%s: Step %i, indicator %i: limit max #quadrants marked for "
            "refinement to %lli\n",
            this_fn_name, amr_step, (int) indicator_type,
            (long long int) n_refine_quadrants_max);

        /* increase max tol of current indicator to meet limit for #quads */
        slabs_discr_amr_indicator_increase_tol_max (
            indicator, marker_next->n_refine_quadrants_loc,
            n_refine_quadrants_glo, n_refine_quadrants_max);

        //YMIR_LDEBUGF (
        //    "%s: Step %i, indicator %i: increase max tol to %g\n",
        //    this_fn_name, amr_step, (int) indicator_type, indicator->tol_max);

        /* update AMR marker with the modified AMR indicator */
        if (indicator_id <= indicator_skipped) { /* if first valid indicator */
          slabs_discr_amr_marker_set_indicator (marker, indicator);
        }
        else { /* if a subsequent valid indicator */
          slabs_discr_amr_marker_merge_indicator (marker, indicator);
        }

        /* update the new number of quadrants with reduced AMR marker */
        slabs_discr_amr_marker_get_quadrant_counts (
            &n_coarsen_quadrants_glo, &n_refine_quadrants_glo,
            &n_quadrants_next, marker);

        YMIR_GLOBAL_INFOF (
            "%s: Step %i, indicator %i: new #quadrants marked for "
            "coarsening %lli, refinement %lli\n",
            this_fn_name, amr_step, (int) indicator_type,
            (long long int) n_coarsen_quadrants_glo,
            (long long int) n_refine_quadrants_glo);

        /* no more new AMR steps */
        amr_max_steps = amr_step;

        /* exit indicator loop */
        slabs_discr_amr_marker_destroy (marker_next);
        slabs_discr_amr_indicator_destroy (indicator);
        break;
      }
      else { /* if max element limit is not reached */
        slabs_discr_amr_marker_destroy (marker);
        marker = marker_next;
      }

      /* destroy AMR indicator */
      slabs_discr_amr_indicator_destroy (indicator);
    } /* end indicator loop */

    /* mark p4est quadrants for coarsening or refinement; destroy marker */
    slabs_discr_amr_marker_flag_p4est (marker);
    slabs_discr_amr_marker_destroy (marker);

    if (!uniform) {
      p4est_gloidx_t      n_marked_quadrants;
      double              rel_n_marked_quadrants;

      /* calculate absolute and relative number of marked quadrants */
      n_marked_quadrants = n_coarsen_quadrants_glo + n_refine_quadrants_glo;
      rel_n_marked_quadrants = ((double) n_marked_quadrants) /
                               ((double) n_quadrants);

      YMIR_GLOBAL_INFOF (
          "%s: Step %i, %lli quadrants (%.2f %%) marked for AMR "
          "(threshold %.2f %%)\n",
          this_fn_name, amr_step, (long long int) n_marked_quadrants,
          100.0 * rel_n_marked_quadrants, 100.0 * amr_rel_threshold);

      /* exit AMR loop if not enough elements were marked */
      if (rel_n_marked_quadrants <= amr_rel_threshold) {
        YMIR_GLOBAL_INFOF (
            "%s: Done AMR step %i since too few elements were marked\n",
            this_fn_name, amr_step);
        break;
      }
    }

    /*
     * On First Loop Iteration, Prepare Stokes State Fields for AMR
     */

    if (amr_step == 0) {
      /* setup state fields for AMR */
      slabs_stokes_state_amr_prepare (state, *mesh);

      /* destroy ymir mesh and pressure element */
      ymir_mesh_destroy (*mesh);
      ymir_pressure_elem_destroy (*press_elem);
      *mesh = NULL;
      *press_elem = NULL;

      /* destroy continuous mangll geometry */
      mangll_p8est_cnodes_destroy (cnodes);
      cnodes = NULL;

      /* set new AMR relative threshold */
      amr_rel_threshold = SLABS_DISCR_AMR_REL_THRESH_MIN;
    }

    /*
     * Refine & Coarsen p4est Mesh and Interpolate Stokes State Fields
     */

#ifdef SLABS_DISCR_AMR_VTK
    /* vtk output of p4est mesh */
    snprintf (filename, BUFSIZ, "slabs_discr_amr_step%02d_coarse", amr_step);
    slabs_vtk_write_p8est (filename, p8est, physics_options->domain_shape);
#endif

    /* coarsen and refine p8est */
    //TODO implement coarsening
    //p8est_coarsen (p8est, 0, slabs_discr_amr_p4est_coarsen,
    //               slabs_discr_amr_p4est_init);
    p8est_refine (p8est, 0, slabs_discr_amr_p4est_refine,
                  slabs_discr_amr_p4est_init);

#ifdef SLABS_DISCR_AMR_VTK
    /* vtk output of p4est mesh */
    snprintf (filename, BUFSIZ, "slabs_discr_amr_step%02d_fine", amr_step);
    slabs_vtk_write_p8est (filename, p8est, physics_options->domain_shape);
#endif

    /* balance p8est */
    if (!uniform) {
      p8est_balance (p8est, P8EST_CONNECT_FULL, NULL);
    }

#ifdef SLABS_DISCR_AMR_VTK
    /* vtk output of p4est mesh */
    snprintf (filename, BUFSIZ, "slabs_discr_amr_step%02d_bal", amr_step);
    slabs_vtk_write_p8est (filename, p8est, physics_options->domain_shape);
#endif

    /* create refined mangll structures */
    YMIR_ASSERT (mangll_fine == NULL);
    slabs_discr_mangll_for_interp_new (&mangll_fine, p8est, discr_options);

    /* interpolate fields */
    slabs_stokes_state_project (state, mangll_coarse, mangll_fine);

    /* destroy coarse mangll */
    mangll_destroy (mangll_coarse);
    mangll_coarse = NULL;

    /*
     * Partition p4est Mesh and Partition Stokes State Fields
     */

    /* partition p8est */
    if (!uniform) {
      switch (discr_options->mesh_partitioning_type) {
      case SL_MESH_PARTITIONING_ELEM:
        p8est_partition_ext (p8est, 1, NULL);
        break;

      case SL_MESH_PARTITIONING_DOF_VEL:
        YMIR_ABORT_NOT_REACHED ();
        /*TODO delete deprecated code below and option mesh_partitioning_type
        p8est_partition_lnodes_ext (
            p8est, NULL,
            3 * discr_options->n_vel_dnodes_per_el_interior,
            3 * discr_options->n_vel_dnodes_per_face_interior,
            3 * discr_options->n_vel_dnodes_per_edge_interior,
            3 * discr_options->n_vel_dnodes_per_corner,
            1);
        */
        break;

      case SL_MESH_PARTITIONING_DOF_VEL_PRESS:
        YMIR_ABORT_NOT_REACHED ();
        /*TODO delete deprecated code below and option mesh_partitioning_type
        p8est_partition_lnodes_ext (
            p8est, NULL,
            3 * discr_options->n_vel_dnodes_per_el_interior +
            discr_options->n_press_dnodes_per_el,
            3 * discr_options->n_vel_dnodes_per_face_interior,
            3 * discr_options->n_vel_dnodes_per_edge_interior,
            3 * discr_options->n_vel_dnodes_per_corner,
            1);
        */
        break;

      default: /* unknown partitioning type */
        YMIR_ABORT_NOT_REACHED ();
      }
    }

#ifdef SLABS_DISCR_AMR_VTK
    /* vtk output of p4est mesh */
    snprintf (filename, BUFSIZ, "slabs_discr_amr_step%02d_part", amr_step);
    slabs_vtk_write_p8est (filename, p8est, physics_options->domain_shape);

    /* save complete p4est data to disk */
    snprintf (filename, BUFSIZ, "slabs_discr_amr_step%02d_part_data", amr_step);
    p8est_save_ext (filename, p8est, 0, 0);
#endif

    /* create partitioned mangll structures */
    YMIR_ASSERT (mangll_partitioned == NULL);
    slabs_discr_mangll_and_cnodes_new (&mangll_partitioned, NULL, p8est,
                                       discr_options);

    /* partition fields */
    slabs_stokes_state_partition (state, mangll_fine, mangll_partitioned);

    /* destroy unpartitioned mangll */
    slabs_discr_mangll_for_interp_destroy (mangll_fine);
    mangll_fine = NULL;

    /* reassign partitioned mangll to be the new coarse mangll in next loop
     * iteration */
    mangll_coarse = mangll_partitioned;
    mangll_partitioned = NULL;

    YMIR_GLOBAL_INFOF ("%s: Done AMR step %i\n", this_fn_name, amr_step);
  } /* end AMR step loop */

  /*
   * Finalize AMR
   */

  /* create refined mesh */
  if (0 < amr_step) { /* if ymir mesh was destroyed before */
    YMIR_ASSERT (mangll_coarse != NULL);
    YMIR_ASSERT (cnodes == NULL);
    YMIR_ASSERT (*mesh == NULL);
    YMIR_ASSERT (*press_elem == NULL);

    /* destroy discontinuous mangll */
    mangll_destroy (mangll_coarse);

    /* create continuous mangll structures */
    slabs_discr_mangll_and_cnodes_new (&mangll_coarse, &cnodes, p8est,
                                       discr_options);

    /* create new ymir mesh and pressure element */
    slabs_discr_ymir_new (mesh, press_elem, mangll_coarse, cnodes,
                          discr_options);

    /* finalize state fields after AMR */
    slabs_stokes_state_amr_finalize (state, *mesh, *press_elem);

    /* postprocess temperature in Stokes state*/
    slabs_physics_postprocess_temperature (state, physics_options);
  }

  /* log new max level */
  if (discr_options->amr_log_maxlevel) {
    MPI_Comm            mpicomm = p8est->mpicomm;
    int                 mpiret;
    p4est_topidx_t      ti;
    int8_t              level_max_loc;
    int8_t              level_max_glo = 0;

    /* find processor local max level of all forests */
    level_max_loc = 0;
    for (ti = p8est->first_local_tree; ti <= p8est->last_local_tree; ++ti) {
      p4est_tree_t       *tree = p4est_tree_array_index (p8est->trees, ti);

      level_max_loc = SC_MAX (level_max_loc, tree->maxlevel);
    }

    /* get processor-global max level */
    mpiret = MPI_Allreduce (&level_max_loc, &level_max_glo, 1, MPI_INT8_T,
                            MPI_MAX, mpicomm); YMIR_CHECK_MPI (mpiret);

    /* output */
    YMIR_GLOBAL_INFOF ("%s: Max level %i\n", this_fn_name, (int) level_max_glo);
  }

  YMIR_GLOBAL_INFOF ("Done %s (%i AMR steps)\n", this_fn_name, amr_step);

  /* print mesh statistics */
  if (0 < amr_step) {
    ymir_monitor_print_global_mesh_stats (*mesh, *press_elem);
  }

  /* return the number of performed AMR steps */
  return amr_step;
}

/**
 * Computes max difference per element.
 */
sc_dmatrix_t *
slabs_compute_max_difference_per_element (ymir_dvec_t *vec)
{
  sc_dmatrix_t       *vec_e = sc_dmatrix_new (vec->ndfields, vec->Np);
  sc_dmatrix_t       *vec_diff = sc_dmatrix_new (3, vec->K);
  double              min, max;
  ymir_locidx_t       elid;
  int                 fieldid, nodeid;

  for (elid = 0; elid < vec->K; elid++) { /* loop over all elements */
    /* find min and max for this element */
    ymir_dvec_get_elem (vec, vec_e, YMIR_STRIDE_COMP, elid, YMIR_READ);
    min = vec_e->e[0][0];
    max = vec_e->e[0][0];
    for (fieldid = 0; fieldid < vec->ndfields; fieldid++) { /* loop over all
                                                             * fields */
      for (nodeid = 0; nodeid < vec->Np; nodeid++) { /* loop over all nodes */
        min = SC_MIN (min, vec_e->e[fieldid][nodeid]);
        max = SC_MAX (max, vec_e->e[fieldid][nodeid]);
      }
    }

    /* set difference, min, and max */
    vec_diff->e[0][elid] = max - min;
    vec_diff->e[1][elid] = min;
    vec_diff->e[2][elid] = max;
  }

  /* destroy */
  sc_dmatrix_destroy (vec_e);

  /* return vector with max difference for each element */
  return vec_diff;
}

/**
 * Computes dynamic range per element.
 */
sc_dmatrix_t *
slabs_compute_dynamic_range_per_element (ymir_dvec_t *vec)
{
  sc_dmatrix_t       *vec_e = sc_dmatrix_new (vec->ndfields, vec->Np);
  sc_dmatrix_t       *vec_dr = sc_dmatrix_new (3, vec->K);
  double              min, max;
  ymir_locidx_t       elid;
  int                 fieldid, nodeid;

  for (elid = 0; elid < vec->K; elid++) { /* loop over all elements */
    /* find min and max for this element */
    ymir_dvec_get_elem (vec, vec_e, YMIR_STRIDE_COMP, elid, YMIR_READ);
    min = fabs (vec_e->e[0][0]);
    max = fabs (vec_e->e[0][0]);
    for (fieldid = 0; fieldid < vec->ndfields; fieldid++) { /* loop over all
                                                             * fields */
      for (nodeid = 0; nodeid < vec->Np; nodeid++) { /* loop over all nodes */
        min = SC_MIN (min, fabs (vec_e->e[fieldid][nodeid]));
        max = SC_MAX (max, fabs (vec_e->e[fieldid][nodeid]));
      }
    }

    /* set dynamic range, min, and max */
    if (min > 0.0) {
      vec_dr->e[0][elid] = max / min;
    }
    else {
      vec_dr->e[0][elid] = -1.0;
    }
    vec_dr->e[1][elid] = min;
    vec_dr->e[2][elid] = max;
  }

  /* destroy */
  sc_dmatrix_destroy (vec_e);

  /* return vector with dynamic range values for each element */
  return vec_dr;
}

