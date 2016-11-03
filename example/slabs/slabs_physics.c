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

#include <slabs_physics.h>
#include <slabs_physics_extended.h>
#include <ymir_comm.h>
#include <ymir_velocity_vec.h>
#include <ymir_stokes_vec.h>
#include <ymir_interp_vec.h>
#include <ymir_mass_vec.h>
#include <ymir_stokes_pc.h>
#include <ymir_gmg_intergrid.h>
#include <slabs_io.h>
#ifdef YMIR_DEBUG
# include <ymir_stress_op.h>
#endif

#ifdef YMIR_DEBUG
#define SL_COARSEN_STOKES_COEFF_VTK
#endif

#ifdef SL_COARSEN_STOKES_COEFF_VTK
#include <ymir_vtk.h>
#endif

#define SL_CLOSEST_PT_NEWTON_RTOL (1.0e8 * SC_EPS)
#define SL_CLOSEST_PT_NEWTON_MAXITER 40

/* pointcloud data for weakzones */
WeakzonePointCloud *slabs_physics_weakzone_pointcloud = NULL;
WeakzoneKDTree     *slabs_physics_weakzone_kdtree = NULL;

/* enumerator for orientations of a point w.r.t. a curve (for 2plates) */
typedef enum
{
  SL_CURVE_ORIENT_TOP,
  SL_CURVE_ORIENT_BOTTOM,
  SL_CURVE_ORIENT_TOP_RIGHT,
  SL_CURVE_ORIENT_BOTTOM_LEFT
}
slabs_subd_edge_orient_enum_t;

/* data for setting nonlinear viscosity */
typedef struct slabs_visc_nl_set_fn_data
{
  ymir_dvec_t        *dvisc_dIIe;
  ymir_dvec_t        *bounds_marker;
  ymir_dvec_t        *yielding_marker;
  ymir_dvec_t        *IIe;
  ymir_dvec_t        *weakzone;
  slabs_physics_options_t  *physics_options;
}
slabs_visc_nl_set_fn_data_t;

/**
 * Computes and stores boundaries of a domain.
 */
void
slabs_physics_compute_domain_bounds (slabs_physics_options_t *physics_options)
{
  const double        brick_dx_dbl = (double) physics_options->domain_brick_dx;
  const double        brick_dy_dbl = (double) physics_options->domain_brick_dy;
  const double        brick_dz_dbl = (double) physics_options->domain_brick_dz;

  /* set radius bounds for all domain types */
  physics_options->domain_radius_min = SL_SHELL_RADIUS_BOTTOM;
  physics_options->domain_radius_max = SL_SHELL_RADIUS_TOP;

  /* set domain-specific bounds */
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
    physics_options->domain_x_min = 0.0;
    physics_options->domain_x_max = 1.0;
    physics_options->domain_y_min = 0.0;
    physics_options->domain_y_max = 1.0;
    physics_options->domain_z_min = 0.0;
    physics_options->domain_z_max = 1.0;
    physics_options->domain_lon_min = -M_PI/8.0;
    physics_options->domain_lon_max = +M_PI/8.0;
    break;

  case SL_DOMAIN_BRICK:
    physics_options->domain_x_min = 0.0;
    physics_options->domain_x_max = brick_dx_dbl / brick_dz_dbl;
    physics_options->domain_y_min = 0.0;
    physics_options->domain_y_max = brick_dy_dbl / brick_dz_dbl;
    physics_options->domain_z_min = 0.0;
    physics_options->domain_z_max = 1.0;
    physics_options->domain_lon_min = -M_PI/8.0 * brick_dy_dbl / brick_dz_dbl;
    physics_options->domain_lon_max = +M_PI/8.0 * brick_dy_dbl / brick_dz_dbl;
    break;

  case SL_DOMAIN_SHELL:
    break;

  case SL_DOMAIN_SHELL_CHUNK:
    physics_options->domain_lon_min = -M_PI/8.0;
    physics_options->domain_lon_max = +M_PI/8.0;
    break;

  case SL_DOMAIN_SHELL_SLICE:
    physics_options->domain_lon_min = -M_PI/8.0 * brick_dy_dbl / brick_dz_dbl;
    physics_options->domain_lon_max = +M_PI/8.0 * brick_dy_dbl / brick_dz_dbl;
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }
}

/**
 * Computes and stores the volume of a domain.
 */
void
slabs_physics_compute_domain_volume (slabs_physics_options_t *physics_options)
{
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_BRICK:
    {
      const double        dx = physics_options->domain_x_max -
                               physics_options->domain_x_min;
      const double        dy = physics_options->domain_y_max -
                               physics_options->domain_y_min;
      const double        dz = physics_options->domain_z_max -
                               physics_options->domain_z_min;

      YMIR_ASSERT (isfinite (physics_options->domain_x_min));
      YMIR_ASSERT (isfinite (physics_options->domain_x_max));
      YMIR_ASSERT (isfinite (physics_options->domain_y_min));
      YMIR_ASSERT (isfinite (physics_options->domain_y_max));
      YMIR_ASSERT (isfinite (physics_options->domain_z_min));
      YMIR_ASSERT (isfinite (physics_options->domain_z_max));
      physics_options->domain_volume = dx * dy * dz;
    }
    break;

  case SL_DOMAIN_SHELL: /* (4/3 * pi * (radius_top^3 - radius_bottom^3)) */
    {
      const double        rmin = physics_options->domain_radius_min;
      const double        rmax = physics_options->domain_radius_max;

      YMIR_ASSERT (isfinite (physics_options->domain_radius_min));
      YMIR_ASSERT (isfinite (physics_options->domain_radius_max));

      physics_options->domain_volume = 4.0 / 3.0 * M_PI *
                                       (rmax*rmax*rmax - rmin*rmin*rmin);
    }
    break;

  case SL_DOMAIN_SHELL_CHUNK: /* (volume shell / 24) */
    {
      const double        rmin = physics_options->domain_radius_min;
      const double        rmax = physics_options->domain_radius_max;

      YMIR_ASSERT (isfinite (physics_options->domain_radius_min));
      YMIR_ASSERT (isfinite (physics_options->domain_radius_max));

      physics_options->domain_volume = 1.0 / 18.0 * M_PI *
                                       (rmax*rmax*rmax - rmin*rmin*rmin);
    }
    break;

  case SL_DOMAIN_SHELL_SLICE: /* (volume shell / 24 * (dx*dy)/dz^2 */
    {
      const double        rmin = physics_options->domain_radius_min;
      const double        rmax = physics_options->domain_radius_max;
      const double        brick_dx_dbl =
                            (double) physics_options->domain_brick_dx;
      const double        brick_dy_dbl =
                            (double) physics_options->domain_brick_dy;
      const double        brick_dz_dbl =
                            (double) physics_options->domain_brick_dz;

      YMIR_ASSERT (isfinite (physics_options->domain_radius_min));
      YMIR_ASSERT (isfinite (physics_options->domain_radius_max));
      YMIR_ASSERT (isfinite (physics_options->domain_brick_dx));
      YMIR_ASSERT (isfinite (physics_options->domain_brick_dy));
      YMIR_ASSERT (isfinite (physics_options->domain_brick_dz));

      physics_options->domain_volume = 1.0 / 18.0 * M_PI *
                                       (rmax*rmax*rmax - rmin*rmin*rmin) *
                                       (brick_dx_dbl * brick_dy_dbl) /
                                       (brick_dz_dbl * brick_dz_dbl);
    }
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }

  YMIR_ASSERT (isfinite (physics_options->domain_volume));
}

/**
 * Computes and stores center of mass.
 */
void
slabs_physics_compute_domain_center (slabs_physics_options_t *physics_options)
{
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_BRICK:
    YMIR_ASSERT (isfinite (physics_options->domain_x_min));
    YMIR_ASSERT (isfinite (physics_options->domain_x_max));
    YMIR_ASSERT (isfinite (physics_options->domain_y_min));
    YMIR_ASSERT (isfinite (physics_options->domain_y_max));
    YMIR_ASSERT (isfinite (physics_options->domain_z_min));
    YMIR_ASSERT (isfinite (physics_options->domain_z_max));
    physics_options->domain_center[0] = 0.5 * (physics_options->domain_x_min +
                                               physics_options->domain_x_max);
    physics_options->domain_center[1] = 0.5 * (physics_options->domain_y_min +
                                               physics_options->domain_y_max);
    physics_options->domain_center[2] = 0.5 * (physics_options->domain_z_min +
                                               physics_options->domain_z_max);
    break;

  case SL_DOMAIN_SHELL:
    physics_options->domain_center[0] = 0.0;
    physics_options->domain_center[1] = 0.0;
    physics_options->domain_center[2] = 0.0;
    break;

  case SL_DOMAIN_SHELL_CHUNK:
  case SL_DOMAIN_SHELL_SLICE:
    YMIR_ASSERT (isfinite (physics_options->domain_radius_min));
    YMIR_ASSERT (isfinite (physics_options->domain_radius_max));
    physics_options->domain_center[0] = 0.0;
    physics_options->domain_center[1] = 0.0;
    physics_options->domain_center[2] = 0.5 *
      (physics_options->domain_radius_min + physics_options->domain_radius_max);
    //TODO taking the middle of the radii is just a crude approximation
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }

  YMIR_ASSERT (isfinite (physics_options->domain_center[0]));
  YMIR_ASSERT (isfinite (physics_options->domain_center[1]));
  YMIR_ASSERT (isfinite (physics_options->domain_center[2]));
}

/**
 * Computes and stores the moment of inertia of a domain.
 */
void
slabs_physics_compute_domain_moment_of_inertia (slabs_physics_options_t
                                                 *physics_options)
{
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE: /* (mass * side_length^2 / 6) */
    {
      double              moment_of_inertia;

      moment_of_inertia = 1.0 / 6.0;
      physics_options->domain_moment_of_inertia[0] = moment_of_inertia;
      physics_options->domain_moment_of_inertia[1] = moment_of_inertia;
      physics_options->domain_moment_of_inertia[2] = moment_of_inertia;
    }
    break;

  case SL_DOMAIN_BRICK: /* x: (mass/12 * (dy^2 + dz^2))
                         * y: (mass/12 * (dx^2 + dz^2))
                         * z: (mass/12 * (dx^2 + dy^2)) */
    {
      const double        vol = physics_options->domain_volume;
      const double        dx = physics_options->domain_x_max -
                               physics_options->domain_x_min;
      const double        dy = physics_options->domain_y_max -
                               physics_options->domain_y_min;
      const double        dz = physics_options->domain_z_max -
                               physics_options->domain_z_min;

      YMIR_ASSERT (isfinite (physics_options->domain_volume));
      YMIR_ASSERT (isfinite (physics_options->domain_x_min));
      YMIR_ASSERT (isfinite (physics_options->domain_x_max));
      YMIR_ASSERT (isfinite (physics_options->domain_y_min));
      YMIR_ASSERT (isfinite (physics_options->domain_y_max));
      YMIR_ASSERT (isfinite (physics_options->domain_z_min));
      YMIR_ASSERT (isfinite (physics_options->domain_z_max));

      physics_options->domain_moment_of_inertia[0] = vol/12.0 * (dy*dy + dz*dz);
      physics_options->domain_moment_of_inertia[1] = vol/12.0 * (dx*dx + dz*dz);
      physics_options->domain_moment_of_inertia[2] = vol/12.0 * (dx*dx + dy*dy);
    }
    break;

  case SL_DOMAIN_SHELL: /* (2/5 * mass * (radius_top^5 - radius_bottom^5)
                         *             / (radius_top^3 - radius_bottom^3)) */
    {
      const double        rmin = physics_options->domain_radius_min;
      const double        rmax = physics_options->domain_radius_max;
      double              moment_of_inertia;

      YMIR_ASSERT (isfinite (physics_options->domain_radius_min));
      YMIR_ASSERT (isfinite (physics_options->domain_radius_max));

      moment_of_inertia = (8.0 / 15.0) * M_PI * (pow(rmax, 5) - pow(rmin, 5));
      physics_options->domain_moment_of_inertia[0] = moment_of_inertia;
      physics_options->domain_moment_of_inertia[1] = moment_of_inertia;
      physics_options->domain_moment_of_inertia[2] = moment_of_inertia;
    }
    break;

  case SL_DOMAIN_SHELL_CHUNK:
    //TODO
    physics_options->domain_moment_of_inertia[0] = 0.0;
    physics_options->domain_moment_of_inertia[1] = 0.0;
    physics_options->domain_moment_of_inertia[2] = 0.0;
    break;

  case SL_DOMAIN_SHELL_SLICE:
    //TODO
    physics_options->domain_moment_of_inertia[0] = 0.0;
    physics_options->domain_moment_of_inertia[1] = 0.0;
    physics_options->domain_moment_of_inertia[2] = 0.0;
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }

  YMIR_ASSERT (isfinite (physics_options->domain_moment_of_inertia[0]));
  YMIR_ASSERT (isfinite (physics_options->domain_moment_of_inertia[1]));
  YMIR_ASSERT (isfinite (physics_options->domain_moment_of_inertia[2]));
}

/**
 * Set up domain variables.
 */
void
slabs_physics_compute_domain_setup (slabs_physics_options_t *physics_options)
{
  slabs_physics_compute_domain_bounds (physics_options);
  slabs_physics_compute_domain_volume (physics_options);
  slabs_physics_compute_domain_center (physics_options);
  slabs_physics_compute_domain_moment_of_inertia (physics_options);
}

/**
 * Computes the radius of a shell domain or the corresponding value for a
 * rectangular domain.
 */
double
slabs_compute_radius (const double x, const double y, const double z,
                      slabs_physics_options_t *physics_options)
{
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_BRICK:
    /* transform interval of `z` to the interval of the shell slice's radius */
    {
      const double    z_max = physics_options->domain_z_max;
      const double    radius_min = physics_options->domain_radius_min;
      const double    radius_max = physics_options->domain_radius_max;

      return z / z_max * (radius_max - radius_min) + radius_min;
    }
    break;

  case SL_DOMAIN_SHELL:
  case SL_DOMAIN_SHELL_CHUNK:
  case SL_DOMAIN_SHELL_SLICE:
    return sqrt (x * x + y * y + z * z);
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }
}

/**
 * Computes the longitude of a shell domain or the corresponding value for a
 * rectangular domain.
 */
 double
slabs_compute_longitude (const double x, const double y, const double z,
                         slabs_physics_options_t *physics_options)
{
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_BRICK:
    /* transform interval of `y` to the interval of the slice's longitude */
    {
      const double        lon_min = physics_options->domain_lon_min;
      const double        lon_max = physics_options->domain_lon_max;

      return y / physics_options->domain_y_max * (lon_max - lon_min) + lon_min;
    }
    break;

  case SL_DOMAIN_SHELL:
  case SL_DOMAIN_SHELL_CHUNK:
  case SL_DOMAIN_SHELL_SLICE:
    return 0.5 * M_PI - atan2 (z, x);
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }
}

/**
 * Computes the radius at the center of an element in a shell domain or the
 * corresponding value for a rectangular domain.
 */
static double
slabs_compute_radius_at_elem_center (const double *x,
                                     const double *y,
                                     const double *z,
                                     const int *_sc_restrict Vmask,
                                     slabs_physics_options_t *physics_options)
{
  double              radii_sum = 0;
  int                 i;

  for (i = 0; i < 8; i++) { /* loop over all vertices */
    const int           nodeid = Vmask[i];

    radii_sum += slabs_compute_radius (x[nodeid], y[nodeid], z[nodeid],
                                       physics_options);
  }

  /* return mean value of radii at vertices */
  return radii_sum / 8.0;
}

/**
 * Converts Cartesian coordinates (x, y, z) to spherical coordinates
 * (r, theta, phi), where the mathematical convention is used.
 */
void
slabs_convert_cartesian_to_spherical_math_conv (const double x,
                                                const double y,
                                                const double z,
                                                double *r,
                                                double *theta,
                                                double *phi)
{
  /* compute radius `r = sqrt(x^2 + y^2 + z^2)` */
  *r = sqrt (x*x + y*y + z*z);
  YMIR_ASSERT (0.0 <= *r);

  /* compute azimuthal angle `theta = arctan(y/x)` with `-pi < theta <= pi` */
  *theta = atan2 (y, x);
  YMIR_ASSERT (isfinite (*theta));
  YMIR_ASSERT (-M_PI - SC_1000_EPS < *theta && *theta <= M_PI + SC_1000_EPS);

  /* compute polar angle `phi = arccos(z/r)` with `0 <= phi <= pi` */
  if (0.0 < *r) {
    YMIR_ASSERT (fabs (z / *r) <= 1.0);
    *phi = acos (z / *r);
  }
  else {
    *phi = 0.0;
  }
  YMIR_ASSERT (isfinite (*phi));
  YMIR_ASSERT (0.0 <= *phi && *phi <= M_PI + SC_1000_EPS);
}

/**
 * Converts Cartesian coordinates (x, y, z) to spherical coordinates
 * (r, phi, theta), where the geophysical convention is used.
 */
void
slabs_convert_cartesian_to_spherical_geo_conv (const double x,
                                               const double y,
                                               const double z,
                                               double *r,
                                               double *phi,
                                               double *theta)
{
  slabs_convert_cartesian_to_spherical_math_conv (x, y, z, r, phi, theta);
}

/**
 * Computes Cartesian coordianates (x, y, z) of the normalized spherical
 * position vector `e_r` given by (r, theta, phi).
 */
static inline void
slabs_compute_spherical_unit_pos_vector (double r, double theta, double phi,
                                         double *pos_vec)
{
  pos_vec[0] = cos (theta) * sin (phi);
  pos_vec[1] = sin (theta) * sin (phi);
  pos_vec[2] = cos (phi);
}

/**
 * Returns whether the element with coordinates (x,y,z) is located in the
 * upper mantle.
 */
int
slabs_physics_elem_in_upper_mantle (const double *x, const double *y,
                                    const double *z, const int *Vmask,
                                    slabs_physics_options_t *physics_options)
{
  const double        radius_um =
    physics_options->viscosity_upper_mantle_radius;
  const double        radius_el_center =
    slabs_compute_radius_at_elem_center (x, y, z, Vmask, physics_options);

  /* return if radius of element center is in upper mantle */
  return (radius_um <= 0.0 || radius_um <= radius_el_center);
}

/**
 * Returns whether the element with coordinates (x,y,z) is located close enough
 * to weak zones.
 */
int
slabs_physics_elem_has_weakzone (const double *x, const double *y,
                                 const double *z, const int *Vmask,
                                 slabs_physics_options_t *physics_options)
{
  const double        radius_um =
    physics_options->viscosity_upper_mantle_radius;
  const double        radius_el_center =
    slabs_compute_radius_at_elem_center (x, y, z, Vmask, physics_options);

  /* return if radius of element center is below 1/2 of upper mantle depth */
  return ( radius_um <= 0.0 ||
           (0.5 * (SL_SHELL_RADIUS_TOP + radius_um)) <= radius_el_center );
}

/**
 * Evaluates linear polynomial.
 */
#define slabs_poly1(x,poly_coeff) ( (poly_coeff)[0] + (poly_coeff)[1]*(x) )

/**
 * Evaluates quadratic polynomial.
 */
#define slabs_poly2(x,poly_coeff) \
  ( (poly_coeff)[0] + (poly_coeff)[1]*(x) + (poly_coeff)[2]*(x)*(x) )

/**
 * Evaluates derivative of quadratic polynomial.
 */
#define slabs_poly2_deriv(x,poly_coeff) \
  ( (poly_coeff)[1] + 2*(poly_coeff)[2]*(x) )

/**
 * Computes interpolating polynomial in 2D via Hermite interpolation.
 */
static double *
slabs_compute_poly2_interpolation (double start_node, double start_val,
                                   double start_deriv,
                                   double end_node, double end_val)
{
  double             *poly_coeff;

  /* allocate coefficients */
  poly_coeff = YMIR_ALLOC (double, 3);

  /* calculate coefficients via Hermite interpolation */
  poly_coeff[2] = ( (end_val - start_val) / (end_node - start_node)
                    - start_deriv ) / (end_node - start_node);
  poly_coeff[1] = start_deriv - 2 * poly_coeff[2] * start_node;
  poly_coeff[0] = start_val - (start_deriv * start_node)
                  + (poly_coeff[2] * start_node * start_node);

  /* return coefficients */
  return poly_coeff;
}

/**
 * Runs Newton's method to find the closest point on the curve of a
 * piecewise linear/quadratic polynomial to a given point.
 */
static int
slabs_closest_pt_newton (double *x, double point_x, double point_y,
                         double rtol, double maxiter,
                         double *poly1_left_coeff, double stitch_node_left,
                         double *poly2_coeff, double stitch_node_right,
                         double *poly1_right_coeff) {
  double              x_prev;
  int                 k;

  /* run Newton iterations */
  for (k = 0; k < maxiter; k++) {
    /* store previous step */
    x_prev = *x;

    /* compute one Newton step */
    if (*x <= stitch_node_left) { /* if at left lin. poly. */
      *x = *x - ( (*x - point_x) + poly1_left_coeff[1]
                  * (slabs_poly1 (*x, poly1_left_coeff) - point_y) )
              / ( 1.0 + poly1_left_coeff[1] * poly1_left_coeff[1] );
    }
    else if (stitch_node_right <= *x) { /* if at right lin. poly. */
      *x = *x - ( (*x - point_x) + poly1_right_coeff[1]
                  * (slabs_poly1 (*x, poly1_right_coeff) - point_y) )
              / ( 1.0 + poly1_right_coeff[1] * poly1_right_coeff[1] );
    }
    else { /* if in middle at quadratic poly. */
      *x = *x - ( (*x - point_x) + (slabs_poly2 (*x, poly2_coeff) - point_y)
                                 * slabs_poly2_deriv (*x, poly2_coeff) )
              / ( 1.0 + slabs_poly2_deriv (*x, poly2_coeff)
                      * slabs_poly2_deriv (*x, poly2_coeff)
                      + (slabs_poly2 (*x, poly2_coeff) - point_y)
                      * 2.0 * poly2_coeff[2] );
    }

    /* check for convergence */
    YMIR_ASSERT (*x != 0.0);
    if (fabs ((x_prev - *x) / *x) < rtol) {
      break;
    }
  }

  /* return number of iterations */
  return k;
}

/**
 * Computes distance and orientation of a point to a quadratic polynomial in 2D.
 */
static double *
slabs_compute_closest_pt_on_poly2 (double point_x, double point_y,
                                   double *poly2_coeff, double start_node,
                                   double start_val, double start_deriv,
                                   double end_node, double end_val,
                                   int *orientation_wrt_curve)
{
  double              poly1_left_coeff[2], poly1_right_coeff[2];
  double              tangent_start[2], tangent_end[2];
  double              orth_start[2], orth_end[2];
  double              line_orientation_start, line_orientation_end;
  double              proj;
  double              dist_sq_mid, dist_sq_start, dist_sq_end;
  //int                 newton_num_iter;
  double             *closest_pt;

  /* check input parameters */
  YMIR_ASSERT (orientation_wrt_curve != NULL);

  /* compute coefficients of linear polynomial extending the quadratic
   * polynomial to the left and right */
  poly1_left_coeff[1] = slabs_poly2_deriv (start_node, poly2_coeff);
  poly1_left_coeff[0] = slabs_poly2 (start_node, poly2_coeff)
                        - poly1_left_coeff[1] * start_node;
  poly1_right_coeff[1] = slabs_poly2_deriv (end_node, poly2_coeff);
  poly1_right_coeff[0] = slabs_poly2 (end_node, poly2_coeff)
                         - poly1_right_coeff[1] * end_node;

  /* compute normalized tangent vectors of curve at start and end points */
  tangent_start[1] = slabs_poly2_deriv (start_node, poly2_coeff);
  tangent_start[0] = 1.0 / sqrt (1.0 + tangent_start[1] * tangent_start[1]);
  tangent_start[1] *= tangent_start[0];
  tangent_end[1] = slabs_poly2_deriv (end_node, poly2_coeff);
  tangent_end[0] = 1.0 / sqrt (1.0 + tangent_end[1] * tangent_end[1]);
  tangent_end[1] *= tangent_end[0];

  /* compute normal vectors to slopes at start and end points */
  orth_start[0] = slabs_poly2_deriv (start_node, poly2_coeff);
  orth_start[1] = -1.0;
  orth_end[0] = slabs_poly2_deriv (end_node, poly2_coeff);
  orth_end[1] = -1.0;

  /* flip signs if x-value is negative*/
  if (orth_start[0] < 0.0) {
    orth_start[0] = -orth_start[0];
    orth_start[1] = -orth_start[1];
  }
  if (orth_end[0] < 0.0) {
    orth_end[0] = -orth_end[0];
    orth_end[1] = -orth_end[1];
  }

  /* compute orientation of target point w.r.t. lines through start and
   * end points along the resp. normals (use cross product) */
  line_orientation_start = orth_start[0] * (point_y - start_val) -
                           orth_start[1] * (point_x - start_node);
  line_orientation_end = orth_end[0] * (point_y - end_val) -
                         orth_end[1] * (point_x - end_node);

  /* compute closest point on curve and set orientation of point w.r.t. curve */
  closest_pt = YMIR_ALLOC (double, 2);
  if (line_orientation_start > 0) { /* if pt. is left of line through start pt.
                                     * with direction normal to start deriv. */
    proj = - (start_node - point_x) * tangent_start[0]
           - (start_val - point_y) * tangent_start[1];
    closest_pt[0] = start_node + proj * tangent_start[0];
    closest_pt[1] = start_val + proj * tangent_start[1];

    *orientation_wrt_curve = SL_CURVE_ORIENT_BOTTOM;
  }
  else if (line_orientation_end < 0) { /* if pt. is right of line through end pt
                                        * with dir. normal to end deriv. */
    proj = - (end_node - point_x) * tangent_end[0]
           - (end_val - point_y) * tangent_end[1];
    closest_pt[0] = end_node + proj * tangent_end[0];
    closest_pt[1] = end_val + proj * tangent_end[1];

    if (point_y <= slabs_poly1 (point_x, poly1_right_coeff)) { /* if pt below */
      *orientation_wrt_curve = SL_CURVE_ORIENT_BOTTOM_LEFT;
    }
    else { /* if point above linear poly. */
      *orientation_wrt_curve = SL_CURVE_ORIENT_TOP_RIGHT;
    }
  }
  else { /* if pt. is between the lines */
    /* set initial guess for Newton's method */
    closest_pt[0] = (start_node + end_node) / 2.0;
    closest_pt[1] = slabs_poly2 (closest_pt[0], poly2_coeff);
    dist_sq_mid = (point_x - closest_pt[0]) * (point_x - closest_pt[0])
                  + (point_y - closest_pt[1]) * (point_y - closest_pt[1]);
    dist_sq_start = (point_x - start_node) * (point_x - start_node)
                    + (point_y - start_val) * (point_y - start_val);
    dist_sq_end = (point_x - end_node) * (point_x - end_node)
                  + (point_y - end_val) * (point_y - end_val);
    if (dist_sq_mid > dist_sq_start || dist_sq_mid > dist_sq_end) {
      if (dist_sq_start < dist_sq_end) {
        closest_pt[0] = start_node;
      }
      else {
        closest_pt[0] = end_node;
      }
    }

    /* find x-value of closest point via Newton's method in 1D */
    /*newton_num_iter = */
    slabs_closest_pt_newton (closest_pt, point_x, point_y,
                             SL_CLOSEST_PT_NEWTON_RTOL,
                             SL_CLOSEST_PT_NEWTON_MAXITER,
                             poly1_left_coeff, start_node,
                             poly2_coeff, end_node,
                             poly1_right_coeff);

    /* compute corresponding y-value on curve */
    if (closest_pt[0] <= start_node) { /* if at left lin. poly. */
      closest_pt[1] = slabs_poly1 (closest_pt[0], poly1_left_coeff);
    }
    else if (end_node <= closest_pt[0]) { /* if at right lin. poly. */
      closest_pt[1] = slabs_poly1 (closest_pt[0], poly1_right_coeff);
    }
    else { /* if in middle at quadratic poly. */
      closest_pt[1] = slabs_poly2 (closest_pt[0], poly2_coeff);
    }

    /*
    if (newton_num_iter == SL_CLOSEST_PT_NEWTON_MAXITER) {
      YMIR_INFOF ("Warning: slabs_compute_closest_pt_on_poly2: "
                  "Max number of Newton iterations reached for "
                  "target point (%g;%g). "
                  "Closest point (%g;%g) possibly wrong.\n",
                  point_x, point_y, closest_pt[0], closest_pt[1]);
    }
    */

    /* set orientation */
    if (point_y <= slabs_poly2 (point_x, poly2_coeff)) { /* if pt. below curve*/
      *orientation_wrt_curve = SL_CURVE_ORIENT_BOTTOM;
    }
    else { /* if pt. above curve */
      *orientation_wrt_curve = SL_CURVE_ORIENT_TOP;
    }
  }

  /* return distance */
  return closest_pt;
}

/**
 * Computes the temperature for brick with two plates, one subducting and
 * one overriding plate by using polynomial interpolation.
 */
static double
slabs_2plates_temperature_poly2 (double r, double lon,
                                 slabs_physics_options_t *physics_options)
{
  double              trench_lon;
  double              dip_angle;
  double              subd_depth, subd_width;
  double              subd_edge_width, subd_edge_smoothwidth;
  double              subd_plate_vel, subd_plate_init_age;
  double              over_plate_age;

  double              lon_min = physics_options->domain_lon_min;
  double              start_node, start_val, start_deriv;
  double              end_node, end_val;
  double             *poly2_coeff;
  int                 orientation_wrt_curve;
  double             *closest_pt;

  double              depth;
  double              ridge_dist;
  double              dist;
  double              vel;
  double              plate_time;
  double              subd_length;
  double              temp_subd = 1.0;
  double              temp_tip;
  double              temp;

  /* set parameters according to physics options */
  trench_lon = physics_options->temp_2plates_trench_longitude;
  dip_angle = physics_options->temp_2plates_dip_angle;
  subd_depth = physics_options->temp_2plates_subd_depth;
  subd_width = physics_options->temp_2plates_subd_width;
  subd_edge_width = physics_options->temp_2plates_subd_edge_width;
  subd_edge_smoothwidth = physics_options->temp_2plates_subd_edge_smoothwidth;
  subd_plate_vel = physics_options->temp_2plates_subd_plate_velocity;
  subd_plate_init_age = physics_options->temp_2plates_subd_plate_initial_age;
  over_plate_age = physics_options->temp_2plates_over_plate_age;

  /* check parameters */
  YMIR_ASSERT (0.0 < trench_lon);
  YMIR_ASSERT (0.0 < dip_angle && dip_angle < 90.0);
  YMIR_ASSERT (0.0 < subd_depth && 0.0 < subd_width);
  YMIR_ASSERT (0.0 <= subd_edge_width && 0.0 <= subd_edge_smoothwidth);
  YMIR_ASSERT (0.0 < subd_plate_vel);
  YMIR_ASSERT (0.0 <= subd_plate_init_age);
  YMIR_ASSERT (0.0 < over_plate_age);

  /* set points for polynomial interpolation */
  start_node = trench_lon;
  start_val = SL_SHELL_RADIUS_TOP;
  start_deriv = tan (-dip_angle / 180.0 * M_PI);
  end_node = start_node + subd_width / SL_EARTH_RADIUS;
  end_val = start_val - subd_depth / SL_EARTH_RADIUS;
  // TODO missing here and below: shell radius top for longitudinal scaling

  /* compute interpolating quadratic polynomial */
  poly2_coeff = slabs_compute_poly2_interpolation (start_node, start_val,
                                                   start_deriv,
                                                   end_node, end_val);

  /* compute closest point on curve and orientation w.r.t. curve */
  closest_pt = slabs_compute_closest_pt_on_poly2 (lon, r,
                                                  poly2_coeff, start_node,
                                                  start_val, start_deriv,
                                                  end_node, end_val,
                                                  &orientation_wrt_curve);

  /* compute "fake" (TODO) distance */
  dist = SL_EARTH_RADIUS * sqrt ( (lon - closest_pt[0]) * (lon - closest_pt[0])
                                  + (r - closest_pt[1]) * (r - closest_pt[1]) );

  /* For computing the temperature, use temperature profile of the plates
   * from halfspace cooling model, erf(z/(2*sqrt(t*kappa))), as base
   * temperature field. */

  /* compute temperature of plates at surface */
  if (orientation_wrt_curve == SL_CURVE_ORIENT_BOTTOM) { /* if subd plate */
    /* calculate age of subducting plate */
    ridge_dist = fabs (lon - lon_min) * SL_EARTH_RADIUS;
    vel = subd_plate_vel / SL_SEC_PER_YEAR;
    plate_time = ridge_dist / vel + subd_plate_init_age * SL_SEC_PER_YEAR;

    /* avoid division by zero */
    plate_time = SC_MAX (1.0e-3, plate_time);
  }
  else  { /* if overriding plate */
    /* calculate age of overriding plate */
    plate_time = over_plate_age * SL_SEC_PER_YEAR;
  }
  depth = (SL_SHELL_RADIUS_TOP - r) * SL_EARTH_RADIUS;
  temp = erf ( depth / (2.0 * sqrt (plate_time * SL_THERM_DIFFUS)) );

  /* compute temperature of the plate subducted inside of mantle */
  if (   orientation_wrt_curve == SL_CURVE_ORIENT_BOTTOM
      || orientation_wrt_curve == SL_CURVE_ORIENT_BOTTOM_LEFT ) {
    /* calculate total lenght of subducted plate */
    subd_length = SL_EARTH_RADIUS * sqrt (
                    (start_node - end_node) * (start_node - end_node)
                    + (start_val - end_val) * (start_val - end_val) );

    /* calculate age of plate */
    if (lon <= trench_lon) {
      ridge_dist = fabs (lon - lon_min) * SL_EARTH_RADIUS;
    }
    else {
      ridge_dist = fabs (trench_lon - lon_min) * SL_EARTH_RADIUS
                   + subd_length * depth / subd_depth;
    }
    vel = subd_plate_vel / SL_SEC_PER_YEAR;
    plate_time = ridge_dist / vel + subd_plate_init_age * SL_SEC_PER_YEAR;

    /* avoid division by zero */
    plate_time = SC_MAX (1.0e-3, plate_time);

    /* compute temperature of subducting plate */
    temp_subd = erf ( dist / (2.0 * sqrt (plate_time * SL_THERM_DIFFUS)) );
  }

  /* smooth the subducting plate's top edge */
  if (   orientation_wrt_curve == SL_CURVE_ORIENT_TOP
      || orientation_wrt_curve == SL_CURVE_ORIENT_TOP_RIGHT ) {
    if (dist < subd_edge_width) {
      /* set constant temperature inside edge width */
      temp_subd = 0.0;
    }
    else {
      /* smooth edge with Gaussian */
      temp_subd = 1.0 - exp (
        - (dist - subd_edge_width) * (dist - subd_edge_width)
        / (0.5 * subd_edge_smoothwidth * subd_edge_smoothwidth) );
    }
  }

  /* smooth tip of subducting plate */
  if (   orientation_wrt_curve == SL_CURVE_ORIENT_TOP_RIGHT
      || orientation_wrt_curve == SL_CURVE_ORIENT_BOTTOM_LEFT ) {
    /* compute distance between tip and closest point on curve */
    dist = SL_EARTH_RADIUS * sqrt (
      (end_node - closest_pt[0]) * (end_node - closest_pt[0])
      + (end_val - closest_pt[1]) * (end_val - closest_pt[1]) );

    /* compute temperature */
    temp_tip = 1.0 -
      exp (- dist*dist / (0.5 * subd_edge_smoothwidth*subd_edge_smoothwidth) );
    temp_subd = SC_MAX (temp_subd, temp_tip);
  }

  /* update global temperature */
  temp = SC_MIN (temp, temp_subd);

  /* destroy */
  YMIR_FREE (closest_pt);
  YMIR_FREE (poly2_coeff);

  /* return temperature */
  YMIR_ASSERT (isfinite (temp));
  return temp;
}

/**
 * Computes the temperature for brick with two plates, one subducting and
 * one overriding plate by using polynomial interpolation.
 */
static double
slabs_brick_2plates_temperature_poly2 (double x, double y, double z,
                                     slabs_physics_options_t *physics_options)
{
  double              lon;
  double              r;

  /* compute radius and longitude */
  r = slabs_compute_radius (x, y, z, physics_options);
  lon = slabs_compute_longitude (x, y, z, physics_options);

  /* compute temperature */
  return slabs_2plates_temperature_poly2 (r, lon, physics_options);
}

/**
 * Computes the background temperature for the brick with two plates,
 * one subducting and one overriding plate.
 *
 * Corresponds to temperature function:
 *   slabs_brick_2plates_temperature_poly2 (...)
 */
static double
slabs_brick_2plates_background_temp_poly2 (double x, double y, double z,
                                     slabs_physics_options_t *physics_options)
{
  /* retrieve original temperature from brick center */
  return slabs_brick_2plates_temperature_poly2 (
      0.5 * physics_options->domain_x_max,
      0.5 * physics_options->domain_y_max,
      z, physics_options);
}

/**
 * Computes the temperature for brick with two plates, one subducting and
 * one overriding plate.
 */
static double
slabs_brick_2plates_temperature_rhea1 (double x, double y, double z,
                                     slabs_physics_options_t *physics_options)
{
  double              y_max = physics_options->domain_y_max;
  double              r;
  double              ridge_dist;
  double              subd_vel;
  double              over_age;
  double              plate_time;
  double              depth;
  double              lon;
  double              lon_min, lon_max, lon_trench;

  double              slab_center, slab_dip;
  double              max_slab_depth;
  double              slab_dist;
  double              Tplate, Tslab;
  double              Terf, Tgauss, weight;
  double              gauss_width, Tcenter;
  double              depth_corr;
  double              value = 1.0;

  double              subd_plate_vel;
  double              over_plate_age;

  /* set parameters according to physics options */
  lon_trench = physics_options->temp_2plates_trench_longitude;
  max_slab_depth = physics_options->temp_2plates_subd_depth;
  subd_plate_vel = physics_options->temp_2plates_subd_plate_velocity;
  over_plate_age = physics_options->temp_2plates_over_plate_age;

  /* check parameters */
  YMIR_ASSERT (0.0 < lon_trench);
  YMIR_ASSERT (0.0 < max_slab_depth);

  /* transform interval of `z` to the interval of the slice's radius */
  r = z * (SL_SHELL_RADIUS_TOP - SL_SHELL_RADIUS_BOTTOM)
        + SL_SHELL_RADIUS_BOTTOM;

  /* transform interval of `z` to the interval of the slice's lon */
  lon = y / y_max * (M_PI / 4.0) - (M_PI / 8.0);

  depth = (SL_SHELL_RADIUS_TOP - r) * SL_EARTH_RADIUS;
  lon_min = - M_PI / 8.0;
  lon_max = M_PI / 8.0;
  slab_dip = (M_PI / 6.0) + (M_PI / 6.0) * (depth / 500.0e3);
  slab_center = lon_trench + (1.0 - r) / tan (slab_dip);
  slab_dist = sin (slab_dip) * fabs (lon - slab_center) * SL_EARTH_RADIUS;

  /* Use temperature profile of the plates from
   * halfspace cooling model: erf(z/(2*sqrt(t*kappa)))
   * as base temperature field */
  if (lon < slab_center) { /* if subducting plate */
    ridge_dist = fabs (lon - lon_min) * SL_EARTH_RADIUS;
    subd_vel = subd_plate_vel / SL_SEC_PER_YEAR;
    plate_time = ridge_dist / subd_vel;
    plate_time = SC_MAX (1.0e-3, plate_time); // avoid division by zero
    Tplate = erf (depth / (2.0 * sqrt (plate_time * SL_THERM_DIFFUS)));
  }
  else { /* if overriding plate */
    ridge_dist = (lon_max - lon) * SL_EARTH_RADIUS;
    over_age = over_plate_age * SL_SEC_PER_YEAR;
    plate_time = over_age;
    Tplate = erf (depth / (2.0 * sqrt (plate_time * SL_THERM_DIFFUS)));
  }

  /* Superimpose slab on base temperature field. */
  /* Compute error function of plate before trench */
  if (lon <= slab_center) {
    /* calculate age of slab */
    ridge_dist = fabs (lon_trench - lon_min) * SL_EARTH_RADIUS + depth;
    subd_vel = subd_plate_vel / SL_SEC_PER_YEAR;
    plate_time = ridge_dist / subd_vel;
    //plate_time = 80.0e6 * SL_SEC_PER_YEAR;
    Terf = erf (slab_dist / (2.0 * sqrt (plate_time * SL_THERM_DIFFUS)));
  }
  else {
    plate_time = 1.0e6 * SL_SEC_PER_YEAR;
    Terf = erf (slab_dist / (2.0 * sqrt (plate_time * SL_THERM_DIFFUS)));
  }

  /* Compute Gaussian function of slab in mantle; conserve buoyancy! */
  gauss_width = 50.0e3; //80e3; //100e3;
  Tcenter = 0.91;       //0.57; //0.45;
  Tgauss = 1.0 - Tcenter * exp (-slab_dist * slab_dist /
                                (2.0 * gauss_width * gauss_width));

  /* Weighted average of Gaussian and error function */
  //weight = exp (-depth / 250e3);
  //weight = 1.0 - depth / max_slab_depth;
  weight = ((1.0 - depth / (2.0 * max_slab_depth)) <= 1.0) ?
    (1.0 - depth / (2.0 * max_slab_depth)) : 1.0;
  weight = (weight > 0.0) ? (weight) : 0.0;
  Tslab = weight * Terf + (1.0 - weight) * Tgauss;
  value = (Tslab < Tplate) ? (Tslab) : (Tplate);

  /* Apply smoothing at the tip of the slab with 2D Gaussian */
  if (lon > slab_center) {
    depth_corr = depth + (slab_dist * cos (slab_dip));
  }
  else {
    depth_corr = depth - (slab_dist * cos (slab_dip));
  }
  if ((depth_corr > max_slab_depth) && (depth > (max_slab_depth - 50.0e3))) {
    Tslab = 1.0 - exp (-((depth_corr - max_slab_depth) *
                         (depth_corr - max_slab_depth) /
                         (2.0 * gauss_width * gauss_width)) -
                       ((slab_dist * slab_dist) /
                        (2.0 * gauss_width * gauss_width)));
    value = (value < Tslab) ? (Tslab) : (value);
  }

  YMIR_ASSERT (isfinite (value));

  return value;
}

/**
 * Computes the background temperature for the brick with two plates,
 * one subducting and one overriding plate.
 *
 * Corresponds to temperature function:
 *   slabs_brick_2plates_temperature_rhea1 (...)
 */
static double
slabs_brick_2plates_background_temp_rhea1 (double x, double y, double z,
                                     slabs_physics_options_t *physics_options)
{
  /* retrieve original temperature from brick center */
  return slabs_brick_2plates_temperature_rhea1 (
      0.5 * physics_options->domain_x_max,
      0.5 * physics_options->domain_y_max,
      z, physics_options);
}

/**
 * Computes the nondimensional scaling factor that is necessary to calculate
 * the temperature distribution with the half-space cooling model (HSCM).
 *
 *   scaling = R / (2 * sqrt(kappa * t))
 *
 * where
 *   R     ... earth radius [m]
 *   kappa ... thermal diffusivity [m^2 / s]
 *   t     ... plate age [yr]
 */
static double
slabs_shell_hscm_plate_age_to_scaling (const double plate_age)
{
  const double        plate_time_sec = plate_age * SL_SEC_PER_YEAR;
  const double        thermal_diffus = SL_THERM_DIFFUS;
  const double        earth_radius = SL_EARTH_RADIUS;

  return earth_radius / (2.0 * sqrt (plate_time_sec * thermal_diffus));
}

/**
 * Computes radially dependent temperature distribution according to the
 * half-space cooling model.
 *
 *   T(r) = erf( c * (r_t - r)),  0 <= r <= r_t
 *
 * Range(T) = [0, 1]
 */
static double
slabs_shell_cold_plate_temperature (const double radius, const double plate_age)
{
  const double        c = slabs_shell_hscm_plate_age_to_scaling (plate_age);

  return erf (c * (SL_SHELL_RADIUS_TOP - radius));
}

/**
 * Computes the radius dependent temperature for the shell with hot core:
 *
 *   T(r) = 1 - 1/2 * erf(c * (r - r_b)),  0 <= r <= r_t
 *
 * Range(T) = [0.5, 1]
 */
static double
slabs_shell_hot_core_temperature (const double radius, const double plate_age)
{
  const double        c = slabs_shell_hscm_plate_age_to_scaling (plate_age);

  return 1.0 - 0.5 * erf (c * (radius - SL_SHELL_RADIUS_BOTTOM));
}

/**
 * Computes the radius dependent temperature for the shell with hot core and
 * cold plate:
 *
 *   T(r) = 1/2 + 1/2 * ( erf(c * (r - r_b)) - erf(c * (r_t - r)) )
 *
 * Range(T) = [0, 1]
 */
static double
slabs_shell_hot_core_cold_plate_temperature (const double radius,
                                             const double plate_age)
{
  const double        c = slabs_shell_hscm_plate_age_to_scaling (plate_age);

  return 0.5 + 0.5 * ( erf (c * (SL_SHELL_RADIUS_TOP - radius)) -
                       erf (c * (radius - SL_SHELL_RADIUS_BOTTOM)) );
}

/**
 * Computes the temperature for shell slice with two plates, one subducting and
 * one overriding plate by using polynomial interpolation.
 */
static double
slabs_shell_slice_2plates_temperature_poly2 (double x, double y, double z,
                                      slabs_physics_options_t *physics_options)
{
  double              lon;
  double              r;

  /* compute radius and longitude */
  r = slabs_compute_radius (x, y, z, physics_options);
  lon = slabs_compute_longitude (x, y, z, physics_options);

  /* compute temperature */
  return slabs_2plates_temperature_poly2 (r, lon, physics_options);
}

/**
 * Computes the background temperature for the shell slice with two plates,
 * one subducting and one overriding plate.
 *
 * Corresponds to temperature function:
 *   slabs_shell_slice_2plates_temperature_poly2 (...)
 */
static double
slabs_shell_slice_2plates_background_temp_poly2 (double x, double y, double z,
                                      slabs_physics_options_t *physics_options)
{
  /* retrieve original temperature from slice center */
  return slabs_shell_slice_2plates_temperature_poly2 (0.0, 0.0,
                                                      sqrt (x*x + y*y + z*z),
                                                      physics_options);
}

/**
 * Computes the temperature for shell slice with two plates, one subducting and
 * one overriding plate.
 *
 * (Copied from rhea, file `example/mantleS/trilinear_mc.c`,
 * function `tilinear_mc_slab_init_temp_function`)
 */
static double
slabs_shell_slice_2plates_temperature_rhea1 (double x, double y, double z,
                                     slabs_physics_options_t * physics_options)
{
  double              r;
  double              ridge_dist;
  double              subd_vel;
  double              over_age;
  double              plate_time;
  double              depth;
  double              y_new, z_new, lon;
  double              lon_min, lon_max, lon_trench;

  double              slab_center, slab_dip;
  double              max_slab_depth;
  double              slab_dist;
  double              Tplate, Tslab;
  double              Terf, Tgauss, weight;
  double              gauss_width, Tcenter;
  double              depth_corr;
  double              value = 1.0;

  double              subd_plate_vel, subd_plate_init_age;
  double              over_plate_age;

  /* set parameters according to physics options */
  lon_trench = physics_options->temp_2plates_trench_longitude;
  max_slab_depth = physics_options->temp_2plates_subd_depth;
  subd_plate_vel = physics_options->temp_2plates_subd_plate_velocity;
  subd_plate_init_age = physics_options->temp_2plates_subd_plate_initial_age;
  over_plate_age = physics_options->temp_2plates_over_plate_age;

  /* check parameters */
  YMIR_ASSERT (0.0 < lon_trench);
  YMIR_ASSERT (0.0 < max_slab_depth);

  /* coordinate transformation: from cartesian to spherical,
     and from chunk at pole to chunk at equator */
  y_new = z;
  z_new = y;

  r = sqrt (x * x + y_new * y_new + z_new * z_new);
  depth = (1.0 - r) * SL_EARTH_RADIUS;
  lon = 0.5 * M_PI - atan2 (y_new, x);
  lon_min = - M_PI / 8.0;
  lon_max = M_PI / 8.0;
  slab_dip = (M_PI / 6.0) + (M_PI / 6.0) * (depth / 500.0e3);
  slab_center = lon_trench + (1.0 - r) / tan (slab_dip);
  slab_dist = sin (slab_dip) * fabs (lon - slab_center) * SL_EARTH_RADIUS;

  /* Use temperature profile of the plates from
   * halfspace cooling model: erf(z/(2*sqrt(t*kappa)))
   * as base temperature field */
  if (lon < slab_center) { /* if subducting plate */
    ridge_dist = fabs (lon - lon_min) * SL_EARTH_RADIUS;
    subd_vel = subd_plate_vel / SL_SEC_PER_YEAR;
    plate_time = ridge_dist / subd_vel + subd_plate_init_age * SL_SEC_PER_YEAR;
    plate_time = SC_MAX (1.0e-3, plate_time); // avoid division by zero
    Tplate = erf (depth / (2.0 * sqrt (plate_time * SL_THERM_DIFFUS)));
  }
  else { /* if overriding plate */
    ridge_dist = (lon_max - lon) * SL_EARTH_RADIUS;
    over_age = over_plate_age * SL_SEC_PER_YEAR;
    plate_time = over_age;
    Tplate = erf (depth / (2.0 * sqrt (plate_time * SL_THERM_DIFFUS)));
  }

  /* Superimpose slab on base temperature field. */
  /* Compute error function of plate before trench */
  if (lon <= slab_center) {
    /* calculate age of slab */
    ridge_dist = fabs (lon_trench - lon_min) * SL_EARTH_RADIUS + depth;
    subd_vel = subd_plate_vel / SL_SEC_PER_YEAR;
    plate_time = ridge_dist / subd_vel;
    //plate_time = 80.0e6 * SL_SEC_PER_YEAR;
    Terf = erf (slab_dist / (2.0 * sqrt (plate_time * SL_THERM_DIFFUS)));
  }
  else {
    plate_time = 1.0e6 * SL_SEC_PER_YEAR;
    Terf = erf (slab_dist / (2.0 * sqrt (plate_time * SL_THERM_DIFFUS)));
  }

  /* Compute Gaussian function of slab in mantle; conserve buoyancy! */
  gauss_width = 50.0e3; //80e3; //100e3;
  Tcenter = 0.91;       //0.57; //0.45;
  Tgauss = 1.0 - Tcenter * exp (-slab_dist * slab_dist /
                                (2.0 * gauss_width * gauss_width));

  /* Weighted average of Gaussian and error function */
  //weight = exp (-depth / 250e3);
  //weight = 1.0 - depth / max_slab_depth;
  weight = ((1.0 - depth / (2.0 * max_slab_depth)) <= 1.0) ?
    (1.0 - depth / (2.0 * max_slab_depth)) : 1.0;
  weight = (weight > 0.0) ? (weight) : 0.0;
  Tslab = weight * Terf + (1.0 - weight) * Tgauss;
  value = (Tslab < Tplate) ? (Tslab) : (Tplate);

  /* Apply smoothing at the tip of the slab with 2D Gaussian */
  if (lon > slab_center) {
    depth_corr = depth + (slab_dist * cos (slab_dip));
  }
  else {
    depth_corr = depth - (slab_dist * cos (slab_dip));
  }
  if ((depth_corr > max_slab_depth) && (depth > (max_slab_depth - 50.0e3))) {
    Tslab = 1.0 - exp (-((depth_corr - max_slab_depth) *
                         (depth_corr - max_slab_depth) /
                         (2.0 * gauss_width * gauss_width)) -
                       ((slab_dist * slab_dist) /
                        (2.0 * gauss_width * gauss_width)));
    value = (value < Tslab) ? (Tslab) : (value);
  }

  YMIR_ASSERT (isfinite (value));

  return value;
}

/**
 * Computes the background temperature for the shell slice with two plates,
 * one subducting and one overriding plate.
 *
 * Corresponds to temperature function:
 *   slabs_shell_slice_2plates_temperature (...)
 */
static double
slabs_shell_slice_2plates_background_temp_rhea1 (double x, double y, double z,
                                     slabs_physics_options_t * physics_options)
{
  /* retrieve original temperature from slice center */
  return slabs_shell_slice_2plates_temperature_rhea1 (0.0, 0.0,
                                                      sqrt (x*x + y*y + z*z),
                                                      physics_options);
}

/**
 * Computes the temperature of a plume (e.g., for right-hand side forcing):
 *
 *   T(x) = exp(-decay * ||center - x||^2)
 */
static double
slabs_plume_temp (double x, double y, double z,
                  double center_x, double center_y, double center_z,
                  double decay, double scaling)
{
  double              dist;

  /* compute square distance from plume center */
  dist = (center_x - x) * (center_x - x) +
         (center_y - y) * (center_y - y) +
         (center_z - z) * (center_z - z);

  /* compute temperature */
  return scaling * exp (-decay * dist);
}

/**
 * Temperature
 */
void
slabs_physics_temperature_set_fn (double *temp, double x, double y, double z,
                                  ymir_locidx_t nid, void *data)
{
  slabs_physics_options_t  *physics_options = (slabs_physics_options_t *) data;
  const double        plate_age = physics_options->temp_background_plate_age;
  double              radius;

  /* compute temperature */
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
    switch (physics_options->temperature_type) {
    case SL_TEMP_NONE:
      *temp = SL_DEFAULT_CONST_TEMP;
      break;

    case SL_TEMP_COLD_PLATE:
      radius = slabs_compute_radius (x, y, z, physics_options);
      *temp = slabs_shell_cold_plate_temperature (radius, plate_age);
      break;

    default: /* unknown temperature type */
      YMIR_ABORT_NOT_REACHED ();
    }
    break;

  case SL_DOMAIN_BRICK:
    switch (physics_options->temperature_type) {
    case SL_TEMP_NONE:
      *temp = SL_DEFAULT_CONST_TEMP;
      break;

    case SL_TEMP_COLD_PLATE:
      radius = slabs_compute_radius (x, y, z, physics_options);
      *temp = slabs_shell_cold_plate_temperature (radius, plate_age);
      break;

    case SL_TEMP_2PLATES_POLY2:
      *temp = slabs_brick_2plates_temperature_poly2 (x, y, z, physics_options);
      break;

    case SL_TEMP_2PLATES_RHEA1:
      *temp = slabs_brick_2plates_temperature_rhea1 (x, y, z, physics_options);
      break;

    default: /* unknown temperature type */
      YMIR_ABORT_NOT_REACHED ();
    }
    break;

  case SL_DOMAIN_SHELL:
    switch (physics_options->temperature_type) {
    case SL_TEMP_NONE:
      *temp = SL_DEFAULT_CONST_TEMP;
      break;

    case SL_TEMP_COLD_PLATE:
      radius = slabs_compute_radius (x, y, z, physics_options);
      *temp = slabs_shell_cold_plate_temperature (radius, plate_age);
      break;

    case SL_TEMP_HOT_CORE:
      radius = slabs_compute_radius (x, y, z, physics_options);
      *temp = slabs_shell_hot_core_temperature (radius, plate_age);
      break;

    case SL_TEMP_HOT_CORE_COLD_PLATE:
      radius = slabs_compute_radius (x, y, z, physics_options);
      *temp = slabs_shell_hot_core_cold_plate_temperature (radius, plate_age);
      break;

    default: /* unknown temperature type */
      YMIR_ABORT_NOT_REACHED ();
    }
    break;

  case SL_DOMAIN_SHELL_CHUNK:
  case SL_DOMAIN_SHELL_SLICE:
    switch (physics_options->temperature_type) {
    case SL_TEMP_NONE:
      *temp = SL_DEFAULT_CONST_TEMP;
      break;

    case SL_TEMP_COLD_PLATE:
      radius = slabs_compute_radius (x, y, z, physics_options);
      *temp = slabs_shell_cold_plate_temperature (radius, plate_age);
      break;

    case SL_TEMP_HOT_CORE:
      radius = slabs_compute_radius (x, y, z, physics_options);
      *temp = slabs_shell_hot_core_temperature (radius, plate_age);
      break;

    case SL_TEMP_HOT_CORE_COLD_PLATE:
      radius = slabs_compute_radius (x, y, z, physics_options);
      *temp = slabs_shell_hot_core_cold_plate_temperature (radius, plate_age);
      break;

    case SL_TEMP_2PLATES_POLY2:
      *temp = slabs_shell_slice_2plates_temperature_poly2 (x, y, z,
                                                           physics_options);
      break;

    case SL_TEMP_2PLATES_RHEA1:
      *temp = slabs_shell_slice_2plates_temperature_rhea1 (x, y, z,
                                                           physics_options);
      break;

    default: /* unknown temperature type */
      YMIR_ABORT_NOT_REACHED ();
    }
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* add plume temperature */
  if (physics_options->plume_type != SL_PLUME_NONE) {
    const double        plume_center_x = physics_options->plume_center_x;
    const double        plume_center_y = physics_options->plume_center_y;
    const double        plume_center_z = physics_options->plume_center_z;
    const double        plume_decay = physics_options->plume_decay;
    const double        plume_scaling = physics_options->plume_scaling;

    *temp += slabs_plume_temp (x, y, z, plume_center_x, plume_center_y,
                               plume_center_z, plume_decay, plume_scaling);
  }

  /* check temperature for `nan` and `inf` */
  YMIR_ASSERT (isfinite (*temp));

  /* bound temperature to valid interval */
  *temp = SC_MIN (1.0, *temp);
  *temp = SC_MAX (0.0, *temp);
}

/**
 * Background temperature
 */
void
slabs_physics_background_temp_set_fn (double *back_temp, double x, double y,
                                      double z, ymir_locidx_t nid, void *data)
{
  slabs_physics_options_t  *physics_options = (slabs_physics_options_t *) data;
  const double        plate_age = physics_options->temp_background_plate_age;
  double              radius;

  /* compute background temperature */
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
    switch (physics_options->temperature_type) {
    case SL_TEMP_NONE:
      *back_temp = SL_DEFAULT_CONST_TEMP;
      break;

    case SL_TEMP_COLD_PLATE:
    case SL_TEMP_IMPORT_FILE:
      radius = slabs_compute_radius (x, y, z, physics_options);
      *back_temp = slabs_shell_cold_plate_temperature (radius, plate_age);
      break;

    default: /* unknown temperature type */
      YMIR_ABORT_NOT_REACHED ();
    }
    break;

  case SL_DOMAIN_BRICK:
    switch (physics_options->temperature_type) {
    case SL_TEMP_NONE:
      *back_temp = SL_DEFAULT_CONST_TEMP;
      break;

    case SL_TEMP_COLD_PLATE:
    case SL_TEMP_IMPORT_FILE:
      radius = slabs_compute_radius (x, y, z, physics_options);
      *back_temp = slabs_shell_cold_plate_temperature (radius, plate_age);
      break;

    case SL_TEMP_2PLATES_POLY2:
      *back_temp = slabs_brick_2plates_background_temp_poly2 (x, y, z,
                                                              physics_options);
      break;

    case SL_TEMP_2PLATES_RHEA1:
      *back_temp = slabs_brick_2plates_background_temp_rhea1 (x, y, z,
                                                              physics_options);
      break;

    default: /* unknown temperature type */
      YMIR_ABORT_NOT_REACHED ();
    }
    break;

  case SL_DOMAIN_SHELL:
    switch (physics_options->temperature_type) {
    case SL_TEMP_NONE:
      *back_temp = SL_DEFAULT_CONST_TEMP;
      break;

    case SL_TEMP_COLD_PLATE:
    case SL_TEMP_IMPORT_FILE:
      radius = slabs_compute_radius (x, y, z, physics_options);
      *back_temp = slabs_shell_cold_plate_temperature (radius, plate_age);
      break;

    case SL_TEMP_HOT_CORE:
      radius = slabs_compute_radius (x, y, z, physics_options);
      *back_temp = slabs_shell_hot_core_temperature (radius, plate_age);
      break;

    case SL_TEMP_HOT_CORE_COLD_PLATE:
      radius = slabs_compute_radius (x, y, z, physics_options);
      *back_temp = slabs_shell_hot_core_cold_plate_temperature (radius,
                                                                plate_age);
      break;

    default: /* unknown temperature type */
      YMIR_ABORT_NOT_REACHED ();
    }
    break;

  case SL_DOMAIN_SHELL_CHUNK:
  case SL_DOMAIN_SHELL_SLICE:
    switch (physics_options->temperature_type) {
    case SL_TEMP_NONE:
      *back_temp = SL_DEFAULT_CONST_TEMP;
      break;

    case SL_TEMP_COLD_PLATE:
    case SL_TEMP_IMPORT_FILE:
      radius = slabs_compute_radius (x, y, z, physics_options);
      *back_temp = slabs_shell_cold_plate_temperature (radius, plate_age);
      break;

    case SL_TEMP_HOT_CORE:
      radius = slabs_compute_radius (x, y, z, physics_options);
      *back_temp = slabs_shell_hot_core_temperature (radius, plate_age);
      break;

    case SL_TEMP_HOT_CORE_COLD_PLATE:
      radius = slabs_compute_radius (x, y, z, physics_options);
      *back_temp = slabs_shell_hot_core_cold_plate_temperature (radius,
                                                                plate_age);
      break;

    case SL_TEMP_2PLATES_POLY2:
      *back_temp = slabs_shell_slice_2plates_background_temp_poly2 (
          x, y, z, physics_options);
      break;

    case SL_TEMP_2PLATES_RHEA1:
      *back_temp = slabs_shell_slice_2plates_background_temp_rhea1 (
          x, y, z, physics_options);
      break;

    default: /* unknown temperature type */
      YMIR_ABORT_NOT_REACHED ();
    }
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* check background temperature for `nan` and `inf` */
  YMIR_ASSERT (isfinite (*back_temp));

  /* bound background temperature to valid interval */
  *back_temp = SC_MIN (1.0, *back_temp);
  *back_temp = SC_MAX (0.0, *back_temp);
}

/**
 *
 */
static void
slabs_temp_enforce_min_plate_age_fn (double *temp, double x, double y, double z,
                                     ymir_locidx_t nid, void *data)
{
  slabs_physics_options_t *physics_options = (slabs_physics_options_t *) data;
  const double        plate_age_min =
                        physics_options->temp_import_plate_age_min;
  double              radius, temp_min;

  /* compute radius */
  radius = slabs_compute_radius (x, y, z, physics_options);

  /* compute min temperature */
  temp_min = slabs_shell_cold_plate_temperature (radius, plate_age_min);
  YMIR_ASSERT (isfinite (temp_min));

  /* update temperature */
  *temp = SC_MIN (*temp, temp_min);
}

/**
 *
 */
void
slabs_temp_postprocess_elem (sc_dmatrix_t *temp_el_mat,
                             const double *x, const double *y, const double *z,
                             slabs_physics_options_t *physics_options)
{
  const int           n_nodes_per_el = temp_el_mat->m;
  double             *temp_el_data = temp_el_mat->e[0];
  int                 nodeid;

  /* check input */
  YMIR_ASSERT (temp_el_mat->n == 1);

  /* compute weak zone factor for each node */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) { /* loop over all
                                                         * nodes */
    /* enforce min plate age */
    if (physics_options->temperature_type == SL_TEMP_IMPORT_FILE &&
        0.0 < physics_options->temp_import_plate_age_min) {
      slabs_temp_enforce_min_plate_age_fn (
          &temp_el_data[nodeid], x[nodeid], y[nodeid], z[nodeid], 0,
          physics_options);
    }

    /* bound temperature to valid range */
    temp_el_data[nodeid] = SC_MAX (0.0, temp_el_data[nodeid]);
    temp_el_data[nodeid] = SC_MIN (1.0, temp_el_data[nodeid]);
  }
}

/**
 *
 */
void
slabs_physics_postprocess_temperature (slabs_stokes_state_t *state,
                                       slabs_physics_options_t
                                         *physics_options)
{
  /* check input */
  YMIR_ASSERT (state->temp_vec != NULL);

  /* enforce min plate age */
  if (physics_options->temperature_type == SL_TEMP_IMPORT_FILE &&
      0.0 < physics_options->temp_import_plate_age_min) {
    ymir_cvec_set_function (state->temp_vec,
                            slabs_temp_enforce_min_plate_age_fn,
                            physics_options);
  }

  /* bound temperature to valid range */
  slabs_cvec_bound_values (state->temp_vec, 0.0, 1.0);
}

/**
 *
 */
void
slabs_physics_verify_temperature (ymir_vec_t *temp_vec,
                                  slabs_physics_options_t *physics_options)
{
  if (physics_options->temperature_type == SL_TEMP_IMPORT_FILE &&
      physics_options->temp_import_verification_out != NULL) {
    /* write temperature to file */
    slabs_io_write_cvec_to_textfile (
        physics_options->temp_import_verification_out,
        temp_vec, SL_CARTESIAN_COORDINATE);
  }
}

/**
 * Computes the weak zone factor depending on the distance to a weak zone.
 * The edges of the weak zone are smoothed by a Gaussian.
 *
 *   1 - (1 - weak_factor) * exp ( - dist^2 / (2 * (0.5*thickness)^2) )
 */
static inline double
slabs_weakzone_factor_fn (const double distance,
                          const double thickness,
                          const double thickness_const,
                          const double weak_factor)
{
  const double        d = distance - 0.5 * thickness_const;
  const double        std_dev = 0.5 * (thickness - thickness_const);

  YMIR_ASSERT (thickness_const <= thickness);

  if (d <= 0.0) {
    /* return value inside zone with constant weak factor */
    return weak_factor;
  }
  else {
    /* return smoothed weak zone */
    return 1.0 - (1.0 - weak_factor) * exp (-d*d / (2.0 * std_dev*std_dev));
  }
}

/**
 *
 */
static int
slabs_weakzone_params_newton (double *x, double y, double rtol, double maxiter,
                              double *poly2_coeff, double weak_dist) {
  double              x_prev;
  int                 k;

  /* run Newton iterations */
  for (k = 0; k < maxiter; k++) {
    /* store previous step */
    x_prev = *x;

    /* compute one Newton step */
    *x = *x - ( slabs_poly2 (*x, poly2_coeff) + weak_dist
                / sqrt (1.0 + slabs_poly2_deriv (*x, poly2_coeff)
                            * slabs_poly2_deriv (*x, poly2_coeff)) - y )
            / ( slabs_poly2_deriv (*x, poly2_coeff) - weak_dist
                * slabs_poly2_deriv (*x, poly2_coeff) * 2.0*poly2_coeff[2]
                / (1.0 + slabs_poly2_deriv (*x, poly2_coeff)
                       * slabs_poly2_deriv (*x, poly2_coeff))
                / sqrt (1.0 + slabs_poly2_deriv (*x, poly2_coeff)
                            * slabs_poly2_deriv (*x, poly2_coeff)) );

    /* check for convergence */
    YMIR_ASSERT (*x != 0.0);
    if (fabs ((x_prev - *x) / *x) < rtol) {
      break;
    }
  }

  /* return number of iterations */
  return k;
}

/**
 * Computes the weak zone parameters longitude, dip angle, and weak zone width
 * from given temperature and weak zone data for the `2plates_poly2` model
 * problem.
 */
void
slabs_2plates_poly2_set_weakzone_params_from_temp (slabs_physics_options_t
                                                     *physics_options)
{
  double              start_node, start_val, start_deriv;
  double              end_node, end_val;
  double             *poly2_coeff;
  double              weak_dist;
  double              radius;
  double              start_lon, end_lon;
#ifdef YMIR_DEBUG
  int                 newton_num_iter;
#endif

  double              weak_start_lon, weak_end_lon;
  double              weak_dip;
  double              weakzone_width;

  /* check input parameters */
  YMIR_ASSERT (physics_options->weakzone_type == SL_WEAKZONE_2PLATES_POLY2);
  YMIR_ASSERT (physics_options->temperature_type == SL_TEMP_2PLATES_POLY2 ||
               physics_options->temperature_type == SL_TEMP_IMPORT_FILE);

  /* set points for polynomial interpolation */
  start_node = physics_options->temp_2plates_trench_longitude;
  start_val = SL_SHELL_RADIUS_TOP;
  start_deriv = tan (-physics_options->temp_2plates_dip_angle / 180.0 * M_PI);
  end_node = start_node
             + physics_options->temp_2plates_subd_width / SL_EARTH_RADIUS;
  end_val = start_val
            - physics_options->temp_2plates_subd_depth / SL_EARTH_RADIUS;

  /* compute interpolating quadratic polynomial */
  poly2_coeff = slabs_compute_poly2_interpolation (start_node, start_val,
                                                   start_deriv,
                                                   end_node, end_val);

  /*
   * compute weak zone parameters at the start
   */

  /* calculate distance between subducting plate and weak zone */
  weak_dist = (  physics_options->temp_2plates_subd_edge_width
             //+ physics_options->temp_2plates_subd_edge_smoothwidth
               + 0.5 * physics_options->weakzone_2plates_subdu_thickness
              ) / SL_EARTH_RADIUS;

  /* run Newton to find start lon. of plate corresponding to weak zone */
  radius = SL_SHELL_RADIUS_TOP;
  start_lon = start_node;
#ifdef YMIR_DEBUG
  newton_num_iter =
#endif
  slabs_weakzone_params_newton (&start_lon, radius, SC_1000_EPS,
                                SL_CLOSEST_PT_NEWTON_MAXITER,
                                poly2_coeff, weak_dist);
  YMIR_ASSERT (newton_num_iter < SL_CLOSEST_PT_NEWTON_MAXITER);

  /* compute longitude of start of weak zone */
  weak_start_lon = start_lon
                   - weak_dist * slabs_poly2_deriv (start_lon, poly2_coeff)
                   / sqrt ( 1.0 + slabs_poly2_deriv (start_lon, poly2_coeff)
                                * slabs_poly2_deriv (start_lon, poly2_coeff) );

  /* compute dip angle of weak zone */
  weak_dip = acos (1.0 / sqrt (1.0 + slabs_poly2_deriv (start_lon, poly2_coeff)
                                   * slabs_poly2_deriv (start_lon, poly2_coeff))
             ) / M_PI * 180.0;

  /*
   * compute weak zone parameters at the end
   */

  /* calculate distance between subducting plate and weak zone */
  weak_dist = (  physics_options->temp_2plates_subd_edge_width
               + physics_options->temp_2plates_subd_edge_smoothwidth
             //- 0.5 * physics_options->weakzone_2plates_subdu_thickness
              ) / SL_EARTH_RADIUS;

  /* run Newton's method to find end lon. of plate corresponding to weak zone */
  radius = SL_SHELL_RADIUS_TOP
           - physics_options->weakzone_2plates_subdu_depth / SL_EARTH_RADIUS;
  end_lon = start_lon;
#ifdef YMIR_DEBUG
  newton_num_iter =
#endif
  slabs_weakzone_params_newton (&end_lon, radius, SC_1000_EPS,
                                SL_CLOSEST_PT_NEWTON_MAXITER,
                                poly2_coeff, weak_dist);
  YMIR_ASSERT (newton_num_iter < SL_CLOSEST_PT_NEWTON_MAXITER);

  /* compute longitude of end of weak zone */
  weak_end_lon = end_lon - weak_dist * slabs_poly2_deriv (end_lon, poly2_coeff)
                 / sqrt ( 1.0 + slabs_poly2_deriv (end_lon, poly2_coeff)
                              * slabs_poly2_deriv (end_lon, poly2_coeff) );

  /* compute width of weak zone */
  weakzone_width = (weak_end_lon - weak_start_lon) * SL_EARTH_RADIUS;

  /* update weak zone parameters */
  physics_options->weakzone_2plates_subdu_longitude = weak_start_lon;
  physics_options->weakzone_2plates_subdu_dip_angle = weak_dip;
  physics_options->weakzone_2plates_subdu_width = weakzone_width;

  /* destroy */
  YMIR_FREE (poly2_coeff);
}

/**
 * Computes distance to weak zone between plates.
 */
 double
slabs_2plates_weakzone_poly2_subdu_dist (double r, double lon,
                                         double start_node,
                                         double start_val,
                                         double start_deriv,
                                         double end_node,
                                         double end_val)
{
  double             *poly2_coeff;
  int                 orientation_wrt_curve;
  double             *closest_pt;
  double              dist;

  /* compute interpolating quadratic polynomial */
  poly2_coeff = slabs_compute_poly2_interpolation (start_node, start_val,
                                                   start_deriv,
                                                   end_node, end_val);

  /* compute closest point on curve and orientation w.r.t. curve */
  closest_pt = slabs_compute_closest_pt_on_poly2 (lon, r,
                                                  poly2_coeff, start_node,
                                                  start_val, start_deriv,
                                                  end_node, end_val,
                                                  &orientation_wrt_curve);

  /* compute distance */
  if (   orientation_wrt_curve == SL_CURVE_ORIENT_TOP
      || orientation_wrt_curve == SL_CURVE_ORIENT_BOTTOM ) {
    /* compute distance to closest point on curve */
    dist = SL_EARTH_RADIUS * sqrt (  SC_SQR (lon - closest_pt[0])
                                   + SC_SQR (r - closest_pt[1]) );
  }
  else if (   orientation_wrt_curve == SL_CURVE_ORIENT_TOP_RIGHT
           || orientation_wrt_curve == SL_CURVE_ORIENT_BOTTOM_LEFT ) {
    /* compute distance to tip of curve */
    dist = SL_EARTH_RADIUS * sqrt (  SC_SQR (end_node - lon)
                                   + SC_SQR (end_val - r) );
  }
  else {
    double              zero = 0.0;  /* no const to avoid warning */
    dist = 1.0 / zero;
    YMIR_ABORT_NOT_REACHED ();
  }

  /* destroy */
  YMIR_FREE (closest_pt);
  YMIR_FREE (poly2_coeff);

  /* return distance to weak zone curve */
  return dist;
}

/**
 * Computes distance to weak zone at ridge.
 */
 double
slabs_2plates_weakzone_poly2_ridge_dist (double r, double lon,
                                         double end_node,
                                         double end_val)
{
  double              dist;

  /* compute distance to weak zone */
  if (lon <= end_node && end_val <= r) { /* if inside weak zone */
    dist = 0.0;
  }
  else { /* if outside weak zone */
    if ( end_node < lon && end_val <= r ) {
      /* compute distance to right edge */
      dist = SL_EARTH_RADIUS * (lon - end_node);
    }
    else if (lon <= end_node && r < end_val) {
      /* compute distance to bottom edge */
      dist = SL_EARTH_RADIUS * (end_val - r);
    }
    else {
      /* compute distance to bottom right corner */
      dist = SL_EARTH_RADIUS * sqrt (  SC_SQR (lon - end_node)
                                     + SC_SQR (end_val - r) );
    }
  }

  /* return distance to ridge weak zone */
  return dist;
}

/**
 * Computes value of a smoothing function for the weak zone.
 * The smoothing is such that the dynamic range of the viscosity per element
 * is as small as possible.
 *
 * This smoothing does not improve convergence of linear solver, we rather use
 * the Gaussian smoothing above.
 */
static inline double
slabs_2plates_weakzone_poly2_smoothing_const_dr (double dist,
                                                 double thickness,
                                                 double smoothwidth,
                                                 double weak_factor)
{
  if ( dist <= (0.5 * thickness) ) {
    /* return value inside of weak zone */
    return weak_factor;
  }
  else {
    double              atan_ratio = 0.1; /* ratio that sets transition dist */
    double              atan_dist;
    double              atan_val, atan_deriv;
    double              exp_factor;
    double              weak;

    /* smooth edge with exp/atan composite function */

    /* calculate distance for transition between `exp` and `atan` */
    atan_dist = 0.5 * thickness + (1.0 - atan_ratio) * smoothwidth;

    /* calculate exponent factor */
    exp_factor = log (1.0 / weak_factor) / smoothwidth;

    if (dist <= atan_dist) { /* if `exp` component */
      return weak_factor * exp ( exp_factor * (dist - 0.5 * thickness) );
    }
    else {
      /* set value and derivative at transition */
      atan_val = weak_factor
                 * exp ( exp_factor * (atan_dist - 0.5 * thickness) );
      atan_deriv = exp_factor * atan_val;

      /* compute `atan` component */
      weak = 1.0 - atan_val - (1.0 - atan_val) / M_PI_2
                            * atan (  atan_deriv * M_PI_2 / (1.0 - atan_val)
                                    * (dist - atan_dist) );

      /* smooth `atan` component with Gaussian */
      weak *= exp (- SC_SQR (dist - atan_dist)
                   / (2.0 * SC_SQR (atan_ratio * smoothwidth)) );

      /* return weak zone value */
      return 1.0 - weak;
    }
  }
}

/**
 * Computes the weak fault zone for the shell slice with two plates.
 */
static double
slabs_2plates_weakzone_poly2 (double r, double lon,
                              slabs_physics_options_t * physics_options)
{
  double              subdu_lon;
  double              subdu_dip_angle;
  double              subdu_depth, subdu_width;
  double              subdu_thickness, subdu_thickness_const;
  double              subdu_weak_factor;
  double              ridge_depth, ridge_width;
  double              ridge_smoothwidth;
  double              ridge_weak_factor;

  double              courtesy_width;
  double              total_thickness;
  double              lon_min = physics_options->domain_lon_min;
  double              start_node, start_val, start_deriv;
  double              end_node, end_val;
  double              dist;
  double              weak = 1.0;

  /* set parameters according to physics options */
  subdu_lon = physics_options->weakzone_2plates_subdu_longitude;
  subdu_dip_angle = physics_options->weakzone_2plates_subdu_dip_angle;
  subdu_depth = physics_options->weakzone_2plates_subdu_depth;
  subdu_width = physics_options->weakzone_2plates_subdu_width;
  subdu_thickness = physics_options->weakzone_2plates_subdu_thickness;
  subdu_thickness_const =
    physics_options->weakzone_2plates_subdu_thickness_const;
  subdu_weak_factor = physics_options->weakzone_2plates_subdu_weak_factor;
  ridge_depth = physics_options->weakzone_2plates_ridge_depth;
  ridge_width = physics_options->weakzone_2plates_ridge_width;
  ridge_smoothwidth = physics_options->weakzone_2plates_ridge_smoothwidth;
  ridge_weak_factor = physics_options->weakzone_2plates_ridge_weak_factor;

  /* check parameters */
  YMIR_ASSERT (0.0 < subdu_lon);
  YMIR_ASSERT (0.0 < subdu_dip_angle);
  YMIR_ASSERT (0.0 < subdu_depth && 0.0 < subdu_width);
  YMIR_ASSERT (0.0 < subdu_thickness);
  YMIR_ASSERT (subdu_thickness_const <= subdu_thickness);
  YMIR_ASSERT (0.0 < subdu_weak_factor && subdu_weak_factor <= 1.0);
  YMIR_ASSERT (0.0 < ridge_depth && 0.0 < ridge_width);
  YMIR_ASSERT (0.0 <= ridge_smoothwidth);
  YMIR_ASSERT (0.0 < ridge_weak_factor && ridge_weak_factor <= 1.0);

  /*
   * set subduction weak zone between plates
   */

  /* set points for polynomial interpolation */
  start_node = subdu_lon;
  start_val = SL_SHELL_RADIUS_TOP;
  start_deriv = tan (-subdu_dip_angle / 180.0 * M_PI);
  end_node = start_node + subdu_width / SL_EARTH_RADIUS * SL_SHELL_RADIUS_TOP;
  end_val = start_val - subdu_depth / SL_EARTH_RADIUS;

  /* only consider point in a rectangle containing the weak zone */
  courtesy_width = subdu_thickness / SL_EARTH_RADIUS;
  total_thickness = (2.0 * subdu_thickness - subdu_thickness_const)
                    / SL_EARTH_RADIUS;
  if (   (  start_node - 0.5 * total_thickness
          / sin (subdu_dip_angle / 180.0 * M_PI) - courtesy_width ) <= lon
      && lon <= (end_node + 0.5 * total_thickness + courtesy_width)
      && (end_val - total_thickness - courtesy_width) <= r ) {
    /* compute distance to subduction weak zone */
    dist = slabs_2plates_weakzone_poly2_subdu_dist (
        r, lon, start_node, start_val, start_deriv, end_node, end_val);

    /* set weak zone factor */
    weak = slabs_weakzone_factor_fn (dist, subdu_thickness,
                                     subdu_thickness_const, subdu_weak_factor);
  }

  /*
   * set ridge weak zone in the left corner of domain, so that subducting
   * plate can "move"
   */

  /* set bottom left corner of weak zone */
  end_node = lon_min + ridge_width / SL_EARTH_RADIUS * SL_SHELL_RADIUS_TOP;
  end_val = SL_SHELL_RADIUS_TOP - ridge_depth / SL_EARTH_RADIUS;

  /* only consider points close to weak zone */
  courtesy_width = 2.0 * ridge_smoothwidth / SL_EARTH_RADIUS;
  if (lon <= (end_node + courtesy_width) && (end_val - courtesy_width) <= r) {
    /* compute distance to ridge weak zone */
    dist = slabs_2plates_weakzone_poly2_ridge_dist (r, lon, end_node, end_val);

    /* set weak zone factor */
    weak = slabs_weakzone_factor_fn (dist, ridge_smoothwidth, 0.0,
                                     ridge_weak_factor);
  }

  /*
   * return weak zone factor
   */

  return weak;
}

/**
 * Computes the weak fault zone for the shell slice with two plates.
 */
static double
slabs_brick_2plates_weakzone_poly2 (const double x, const double y,
                                    const double z,
                                    slabs_physics_options_t *physics_options)
{
  double              lon;
  double              r;

  /* compute radius and longitude */
  r = slabs_compute_radius (x, y, z, physics_options);
  lon = slabs_compute_longitude (x, y, z, physics_options);

  /* compute weak zone factor */
  return slabs_2plates_weakzone_poly2 (r, lon, physics_options);
}

/**
 * Computes the weak fault zone for the shell slice with two plates.
 *
 * (Copied from rhea, file `example/mantleS/trilinear_mc.c`,
 * function `tilinear_mc_weak_zone`)
 */
static double
slabs_brick_2plates_weakzone_rhea1 (const double x, const double y,
                                    const double z,
                                    slabs_physics_options_t *physics_options)
{
  double              y_max = physics_options->domain_y_max;
  double              dummy1, dummy2;
  double              r;
  double              slab_dip;
  double              center_width, smooth_width, weak_depth, weak_val;
  double              subd_start_weak_width, subd_start_weak_depth;
  double              dist_trench, dist_weak;
  double              depth;
  double              lon;
  double              lon_min, lon_trench;
  double              func1, func2;
  double              value = 1.0;

  double              slab_center, slab_dist;

  lon_trench = physics_options->weakzone_2plates_subdu_longitude;
  weak_depth = physics_options->weakzone_2plates_subdu_depth;
  center_width = physics_options->weakzone_2plates_subdu_thickness_const;
  smooth_width = physics_options->weakzone_2plates_subdu_thickness -
                 physics_options->weakzone_2plates_subdu_thickness_const;
  weak_val = physics_options->weakzone_2plates_subdu_weak_factor;

  subd_start_weak_width = physics_options->weakzone_2plates_ridge_width;
  subd_start_weak_depth = physics_options->weakzone_2plates_ridge_depth;

  /* check parameters */
  YMIR_ASSERT (0.0 < lon_trench);
  YMIR_ASSERT (0.0 < weak_depth);

  slab_dip = M_PI / 6.0 * 1.65; /* slab dip angle, clockwise from surface */

  /* transform interval of `z` to the interval of the shell slice's radius */
  r = z * (SL_SHELL_RADIUS_TOP - SL_SHELL_RADIUS_BOTTOM)
        + SL_SHELL_RADIUS_BOTTOM;

  /* transform interval of `z` to the interval of the shell slice's lon */
  lon = y / y_max * (M_PI / 4.0) - (M_PI / 8.0);

  lon_min = - M_PI / 8.0;

  depth = (SL_SHELL_RADIUS_TOP - r) * SL_EARTH_RADIUS;
  slab_dip = slab_dip + slab_dip * (depth / 500.0e3);
  slab_center = lon_trench + (1.0 - r) / tan (slab_dip);
  slab_dist = sin (slab_dip) * fabs (lon - slab_center) * SL_EARTH_RADIUS;
  dist_trench = slab_dist;

  if (lon >= lon_trench) { /* weak zone at tip of subducting plate */
    if ((depth <= dist_trench * tan (slab_dip)) &&
        (depth >= (dist_trench - smooth_width) * tan (slab_dip))) {
      /* smoothing zone, subducting plate side */
      dist_weak =
        fabs ((depth / tan (slab_dip)) + smooth_width - dist_trench);

      func1 = cos (M_PI * dist_weak / smooth_width);
      func2 = cos (M_PI * (depth - weak_depth) / smooth_width);

      if (depth <= weak_depth) {
        value = weak_val + (1.0 - weak_val) * (0.5 - 0.5 * func1);
      }
      else if ((depth > weak_depth) && (depth <= (weak_depth + smooth_width))) {
        dummy1 = weak_val + (1.0 - weak_val) * (0.5 - 0.5 * func1);
        dummy2 = weak_val + (1.0 - weak_val) * (0.5 - 0.5 * func2);
        value = (dummy1 > dummy2) ? (dummy1) : (dummy2);
      }

    }
    else
      if ((depth <=
           (dist_trench - center_width + smooth_width) * tan (slab_dip))
          && (depth >= (dist_trench - center_width) * tan (slab_dip))) {
      /* smoothing zone, overriding plate side */
      dist_weak =
        fabs (smooth_width -
              ((depth / tan (slab_dip) + center_width - dist_trench)));

      func1 = cos (M_PI * dist_weak / smooth_width);
      func2 = cos (M_PI * (depth - weak_depth) / smooth_width);

      if (depth <= weak_depth) {
        value = weak_val + (1.0 - weak_val) * (0.5 - 0.5 * func1);
      }
      else if ((depth > weak_depth) && (depth <= (weak_depth + smooth_width))) {
        dummy1 = weak_val + (1.0 - weak_val) * (0.5 - 0.5 * func1);
        dummy2 = weak_val + (1.0 - weak_val) * (0.5 - 0.5 * func2);
        value = (dummy1 > dummy2) ? (dummy1) : (dummy2);
      }

    }
    else if ((depth <= (dist_trench - smooth_width) * tan (slab_dip)) &&
             (depth >=
              (dist_trench - center_width + smooth_width) * tan (slab_dip))) {
      /* center of weak zone */
      func2 = cos (M_PI * (depth - weak_depth) / smooth_width);

      if (depth <= weak_depth) {
        value = weak_val;
      }
      else if ((depth > weak_depth) && (depth <= (weak_depth + smooth_width))) {
        value = weak_val + (1.0 - weak_val) * (0.5 - 0.5 * func2);
      }
    }
  }
  else if (lon <= lon_min + (subd_start_weak_width + smooth_width)
                            / SL_EARTH_RADIUS) {
    /* weak zone at far end of subducting plate */

    //smooth_width = 5.0e3;
    //weak_val = 1.0e-5;

    if ( depth <= subd_start_weak_depth &&
         lon <= lon_min + subd_start_weak_width / SL_EARTH_RADIUS ) {
      /* center of weak zone */
      dist_weak = 0.0;
    }
    else if ( depth <= (subd_start_weak_depth + smooth_width) &&
              lon <= lon_min + subd_start_weak_width / SL_EARTH_RADIUS ) {
      /* smoothing at bottom of weak zone */
      dist_weak = depth - subd_start_weak_depth;
    }
    else if ( depth < subd_start_weak_depth &&
              lon <= lon_min + (subd_start_weak_width + smooth_width)
                               / SL_EARTH_RADIUS ) {
      /* smoothing at side of weak zone */
      dist_weak = fabs (lon - lon_min) * SL_EARTH_RADIUS
                  - subd_start_weak_width ;
    }
    else if ( depth <= (subd_start_weak_depth + smooth_width) &&
              lon <= lon_min + (subd_start_weak_width + smooth_width)
                               / SL_EARTH_RADIUS ) {
      /* bottom corner of weak zone */
      dummy1 = depth - subd_start_weak_depth;
      dummy2 = fabs (lon - lon_min) * SL_EARTH_RADIUS - subd_start_weak_width;
      dist_weak = (dummy1 > dummy2) ? (dummy1) : (dummy2);
    }
    else {
      /* no weak zone below (subd_start_weak_depth + smooth_width) */
      dist_weak = smooth_width;
    }
    func1 = cos (M_PI * dist_weak / smooth_width);

    value = weak_val + (1.0 - weak_val) * (0.5 - 0.5 * func1);
  }
  else {
    value = 1.0;
  }
  return value;
}

/**
 * Computes the weak fault zone for the shell slice with two plates.
 */
static double
slabs_shell_slice_2plates_weakzone_poly2 (
    const double x, const double y, const double z,
    slabs_physics_options_t *physics_options)
{
  double              lon;
  double              r;

  /* compute radius and longitude */
  r = slabs_compute_radius (x, y, z, physics_options);
  lon = slabs_compute_longitude (x, y, z, physics_options);

  /* compute weak zone factor */
  return slabs_2plates_weakzone_poly2 (r, lon, physics_options);
}

/**
 * Computes the weak fault zone for the shell slice with two plates.
 *
 * (Copied from rhea, file `example/mantleS/trilinear_mc.c`,
 * function `tilinear_mc_weak_zone`)
 */
static double
slabs_shell_slice_2plates_weakzone_rhea1 (
    const double x, const double y, const double z,
    slabs_physics_options_t *physics_options)
{
  double              dummy1, dummy2;
  double              r;
  double              slab_dip;
  double              center_width, smooth_width, weak_depth, weak_val;
  double              dist_trench, dist_weak;
  double              depth;
  double              y_new, z_new, lon;
  double              lon_min, lon_trench;
  double              func1, func2;
  double              value = 1.0;

  double              slab_center, slab_dist;

  lon_trench = physics_options->weakzone_2plates_subdu_longitude;
  weak_depth = physics_options->weakzone_2plates_subdu_depth;
  center_width = physics_options->weakzone_2plates_subdu_thickness_const;
  smooth_width = physics_options->weakzone_2plates_subdu_thickness -
                 physics_options->weakzone_2plates_subdu_thickness_const;
  weak_val = physics_options->weakzone_2plates_subdu_weak_factor;
  //lon_trench = 0.128; //0.12;              /* x position of trench, nondim */
  //center_width = 10.0e3 * 4.0;  /* width of maximum value in weak zone, m */
  //smooth_width = 5.0e3 * 4.0; /* width of smoothing on edges of max zone, m */
  //weak_depth = 50.0e3;          /* max depth of max zone, m */
  //weak_val = 1.0e-3; //1.0e-5;  /* strongest value of weak zone */

  /* check parameters */
  YMIR_ASSERT (0.0 < lon_trench);
  YMIR_ASSERT (0.0 < weak_depth);

  slab_dip = M_PI / 6.0 * 1.65; /* slab dip angle, clockwise from surface */

  /* coordinate transformation: from cartesian to spherical,
     and from chunk at pole to chunk at equator */
  y_new = z;
  z_new = y;

  r = sqrt (x * x + y_new * y_new + z_new * z_new);
  lon = 0.5 * M_PI - atan2 (y_new, x);
  lon_min = - M_PI / 8.0;

  depth = (1.0 - r) * SL_EARTH_RADIUS;
  slab_dip = slab_dip + slab_dip * (depth / 500.0e3);
  slab_center = lon_trench + (1.0 - r) / tan (slab_dip);
  slab_dist = sin (slab_dip) * fabs (lon - slab_center) * SL_EARTH_RADIUS;
  dist_trench = slab_dist;
  //dist_trench = fabs (lon - lon_trench) * SL_EARTH_RADIUS;

  if (lon >= lon_trench) {      /* weak zone at tip of subducting plate */
    if ((depth <= dist_trench * tan (slab_dip)) &&
        (depth >= (dist_trench - smooth_width) * tan (slab_dip))) {
      /* smoothing zone, subducting plate side */
      dist_weak =
        fabs ((depth / tan (slab_dip)) + smooth_width - dist_trench);

      func1 = cos (M_PI * dist_weak / smooth_width);
      //func1 = exp (- dist_weak * dist_weak / (2.0 * smooth_width * smooth_width));
      func2 = cos (M_PI * (depth - weak_depth) / smooth_width);

      if (depth <= weak_depth) {
        value = weak_val + (1.0 - weak_val) * (0.5 - 0.5 * func1);
      }
      else if ((depth > weak_depth) && (depth <= (weak_depth + smooth_width))) {
        dummy1 = weak_val + (1.0 - weak_val) * (0.5 - 0.5 * func1);
        dummy2 = weak_val + (1.0 - weak_val) * (0.5 - 0.5 * func2);
        value = (dummy1 > dummy2) ? (dummy1) : (dummy2);
      }

    }
    else
      if ((depth <=
           (dist_trench - center_width + smooth_width) * tan (slab_dip))
          && (depth >= (dist_trench - center_width) * tan (slab_dip))) {
      /* smoothing zone, overriding plate side */
      dist_weak =
        fabs (smooth_width -
              ((depth / tan (slab_dip) + center_width - dist_trench)));

      func1 = cos (M_PI * dist_weak / smooth_width);
      //func1 = exp (- dist_weak * dist_weak / (2.0 * smooth_width * smooth_width));
      func2 = cos (M_PI * (depth - weak_depth) / smooth_width);

      if (depth <= weak_depth) {
        value = weak_val + (1.0 - weak_val) * (0.5 - 0.5 * func1);
      }
      else if ((depth > weak_depth) && (depth <= (weak_depth + smooth_width))) {
        dummy1 = weak_val + (1.0 - weak_val) * (0.5 - 0.5 * func1);
        dummy2 = weak_val + (1.0 - weak_val) * (0.5 - 0.5 * func2);
        value = (dummy1 > dummy2) ? (dummy1) : (dummy2);
      }

    }
    else if ((depth <= (dist_trench - smooth_width) * tan (slab_dip)) &&
             (depth >=
              (dist_trench - center_width + smooth_width) * tan (slab_dip))) {
      /* center of weak zone */
      func2 = cos (M_PI * (depth - weak_depth) / smooth_width);

      if (depth <= weak_depth) {
        value = weak_val;
      }
      else if ((depth > weak_depth) && (depth <= (weak_depth + smooth_width))) {
        value = weak_val + (1.0 - weak_val) * (0.5 - 0.5 * func2);
      }
    }
  }
  else if (lon <= (lon_min + 40e3 / SL_EARTH_RADIUS)) {        /* weak zone at far end of subducting plate */
    center_width = 30.0e3; //20.0e3;
    smooth_width = 5.0e3;
    weak_depth = 30.0e3;
    //weak_val = 1.0e-5;

    if ((depth <= weak_depth)
        && (lon <= center_width / SL_EARTH_RADIUS)) {
      /* center of weak zone */
      dist_weak = 0.0;
    }
    else if ((depth <= (weak_depth + smooth_width)) &&
             (lon <= center_width / SL_EARTH_RADIUS)) {
      /* smoothing at bottom of weak zone */
      dist_weak = depth - weak_depth;
    }
    else if ((depth < weak_depth) && (lon <= (center_width +
                                              smooth_width) /
                                      SL_EARTH_RADIUS)) {
      /* smoothing at side of weak zone */
      dist_weak = (fabs (lon - lon_min) * SL_EARTH_RADIUS) - center_width;
    }
    else if ((depth <= (weak_depth + smooth_width)) &&
             (lon <= (center_width + smooth_width) / SL_EARTH_RADIUS)) {
      /* bottom corner of weak zone */
      dummy1 = depth - weak_depth;
      dummy2 = (fabs (lon - lon_min) * SL_EARTH_RADIUS) - center_width;
      dist_weak = (dummy1 > dummy2) ? (dummy1) : (dummy2);
    }
    else {
      dist_weak = smooth_width; /* no weak zone below (weak_depth+smooth_width) */
    }
    func1 = cos (M_PI * dist_weak / smooth_width);
    //func1 = exp (- dist_weak * dist_weak / (2.0 * smooth_width * smooth_width));

    value = weak_val + (1.0 - weak_val) * (0.5 - 0.5 * func1);
  }
  else {
    value = 1.0;
  }
  return value;
}

/**
 *
 */
void
slabs_physics_init_weakzone (MPI_Comm mpicomm,
                             slabs_physics_options_t *physics_options)
{
  if (physics_options->weakzone_type == SL_WEAKZONE_IMPORT_FILE) {
#ifndef SLABS_IO_USE_MPI_IO
    const char         *filename_txt =
                          physics_options->weakzone_import_filename_txt;
    const char         *filename_bin =
                          physics_options->weakzone_import_filename_bin;
#else
    char               *filename_txt =
                          physics_options->weakzone_import_filename_txt;
    char               *filename_bin =
                          physics_options->weakzone_import_filename_bin;
#endif
    const int           n_points_bin =
                          physics_options->weakzone_import_n_points;
    int                 n_points_txt, n_points;
    double             *point;
    WeakzonePointCloud *cloud;
    WeakzoneKDTree     *tree;

    /* check state */
    YMIR_ASSERT (slabs_physics_weakzone_pointcloud == NULL);
    YMIR_ASSERT (slabs_physics_weakzone_kdtree == NULL);

    /* read weak zone point coordinates from file */
    point = slabs_io_read_weakzone (mpicomm, filename_txt, filename_bin,
                                    &n_points_txt, n_points_bin);

    /* do not read txt file any more, since content was copied to bin file */
    physics_options->weakzone_import_filename_txt = NULL;

    /* set number of points */
    if (filename_txt != NULL) {
      n_points = n_points_txt;
      physics_options->weakzone_import_n_points = n_points_txt;
    }
    else {
      n_points = n_points_bin;
    }
    YMIR_ASSERT (0 < n_points);

    /* create weak zone point cloud with kd-tree interface */
    cloud = WeakzonePointCloud_new ();
    WeakzonePointCloud_set_points (cloud, point, n_points);
    slabs_physics_weakzone_pointcloud = cloud;

    /* create weak zone kd-tree */
    tree = WeakzoneKDTree_new (cloud);
    slabs_physics_weakzone_kdtree = tree;

    /* write weak zone points to textfile */
    if (physics_options->weakzone_import_verification_out != NULL) {
      int                 mpirank, mpiret;

      mpiret = MPI_Comm_rank (mpicomm, &mpirank); YMIR_CHECK_MPI (mpiret);

      if (mpirank == 0) {
        slabs_io_write_double_vec_to_textfile (
            physics_options->weakzone_import_verification_out,
            point, n_points, 3, 0);
      }
    }

    /* destroy point coordinates */
    YMIR_FREE (point);
  }
}

/**
 *
 */
void
slabs_physics_clear_weakzone (slabs_physics_options_t *physics_options)
{
  if (physics_options->weakzone_type == SL_WEAKZONE_IMPORT_FILE) {
    /* check state */
    YMIR_ASSERT (slabs_physics_weakzone_pointcloud != NULL);
    YMIR_ASSERT (slabs_physics_weakzone_kdtree != NULL);

    /* destroy weak zone kd-tree */
    WeakzoneKDTree_destroy (slabs_physics_weakzone_kdtree);
    slabs_physics_weakzone_kdtree = NULL;

    /* destroy weak zone point cloud */
    WeakzonePointCloud_destroy (slabs_physics_weakzone_pointcloud);
    slabs_physics_weakzone_pointcloud = NULL;
  }
}

/**
 * Computes weak zone factor at a node.
 */
static inline double
slabs_weak_node (const double x, const double y, const double z,
                 slabs_physics_options_t *physics_options)
{
  /* compute weak zone factor based on point cloud data */
  if (physics_options->weakzone_type == SL_WEAKZONE_IMPORT_FILE) {
    const double        thickness =
      physics_options->weakzone_import_thickness / SL_EARTH_RADIUS;
    const double        thickness_const =
      physics_options->weakzone_import_thickness_const / SL_EARTH_RADIUS;
    const double        weak_factor =
      physics_options->weakzone_import_weak_factor;
    WeakzoneKDTree     *tree = slabs_physics_weakzone_kdtree;
    const double        pt[3] = {x, y, z};
    double              dist;

    YMIR_ASSERT (tree != NULL);

    /* compute distance to weak zone via kd-tree nearest neighbor search on
     * point cloud */
    dist = WeakzoneKDTree_find_shortest_distance_single (tree, pt);

    /* compute and return weak zone factor at this node */
    return slabs_weakzone_factor_fn (dist, thickness, thickness_const,
                                     weak_factor);
  }

  /* compute weak zone factor if the weak zones are given by functions */
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_SHELL:
    switch (physics_options->weakzone_type) {
    case SL_WEAKZONE_NONE:
      return 1.0;

    default: /* unknown weak zone type */
      YMIR_ABORT_NOT_REACHED ();
    }
    break;

  case SL_DOMAIN_BRICK:
    switch (physics_options->weakzone_type) {
    case SL_WEAKZONE_NONE:
      return 1.0;

    case SL_WEAKZONE_2PLATES_POLY2:
      return slabs_brick_2plates_weakzone_poly2 (x, y, z, physics_options);

    case SL_WEAKZONE_2PLATES_RHEA1:
      return slabs_brick_2plates_weakzone_rhea1 (x, y, z, physics_options);

    default: /* unknown weak zone type */
      YMIR_ABORT_NOT_REACHED ();
    }
    break;

  case SL_DOMAIN_SHELL_CHUNK:
  case SL_DOMAIN_SHELL_SLICE:
    switch (physics_options->weakzone_type) {
    case SL_WEAKZONE_NONE:
      return 1.0;

    case SL_WEAKZONE_2PLATES_POLY2:
      return slabs_shell_slice_2plates_weakzone_poly2 (x, y, z,
                                                       physics_options);

    case SL_WEAKZONE_2PLATES_RHEA1:
      return slabs_shell_slice_2plates_weakzone_rhea1 (x, y, z,
                                                       physics_options);

    default: /* unknown weak zone type */
      YMIR_ABORT_NOT_REACHED ();
    }
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }
}

/**
 * Computes weak zone factor of an element.
 */
void
slabs_weak_elem (sc_dmatrix_t *weak_el_mat,
                 const double *x, const double *y, const double *z,
                 const int *Vmask, slabs_physics_options_t *physics_options)
{
  const int           n_nodes_per_el = weak_el_mat->m;
  double             *weak_el_data = weak_el_mat->e[0];
  int                 nodeid;

  /* check input */
  YMIR_ASSERT (weak_el_mat->n == 1);

  /* set weak factor to 1 if element is not under the influence of weak zones */
  if (!slabs_physics_elem_has_weakzone (x, y, z, Vmask, physics_options)) {
    sc_dmatrix_set_value (weak_el_mat, 1.0);
    return;
  }

  /* compute weak zone factor for each node */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) { /* loop over all
                                                         * nodes */
    weak_el_data[nodeid] = slabs_weak_node (x[nodeid], y[nodeid], z[nodeid],
                                            physics_options);
  }
}

/**
 * Computes weak zone factor.
 */
void
slabs_physics_compute_weakzone (ymir_dvec_t *weakzone, void *data)
{
  slabs_physics_options_t *physics_options = (slabs_physics_options_t *) data;
  ymir_mesh_t        *mesh = weakzone->mesh;
  mangll_t           *mangll = mesh->ma;
  const mangll_locidx_t  n_elements = mesh->cnodes->K;
  const unsigned int  N = ymir_n (mangll->N);
  const unsigned int  n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);

  sc_dmatrix_t       *weak_el_mat;
  double             *x, *y, *z, *tmp_el;
  mangll_locidx_t     elid;

  /* check input */
  YMIR_ASSERT (weakzone->ndfields == 1);

  /* create work variables */
  weak_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  x = YMIR_ALLOC (double, n_nodes_per_el);
  y = YMIR_ALLOC (double, n_nodes_per_el);
  z = YMIR_ALLOC (double, n_nodes_per_el);
  tmp_el = YMIR_ALLOC (double, n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    /* get coordinates of this element at Gauss nodes */
    slabs_elem_get_gauss_coordinates (x, y, z, elid, mangll, tmp_el);

    /* compute weak zone factor */
    slabs_weak_elem (weak_el_mat, x, y, z, mangll->refel->Vmask,
                     physics_options);

    /* set weak zone of this element */
    ymir_dvec_set_elem (weakzone, weak_el_mat, YMIR_STRIDE_NODE, elid,
                        YMIR_SET);
  }

  /* destroy */
  sc_dmatrix_destroy (weak_el_mat);
  YMIR_FREE (x);
  YMIR_FREE (y);
  YMIR_FREE (z);
  YMIR_FREE (tmp_el);
}

/**
 * Computes min distance to a weak zone at a node.
 */
static inline double
slabs_physics_weak_dist_node (const double x, const double y, const double z,
                              slabs_physics_options_t *physics_options)
{
  /* compute weak zone distance based on point cloud data */
  if (physics_options->weakzone_type == SL_WEAKZONE_IMPORT_FILE) {
    WeakzoneKDTree     *tree = slabs_physics_weakzone_kdtree;
    const double        pt[3] = {x, y, z};

    YMIR_ASSERT (tree != NULL);

    /* compute distance to weak zone via kd-tree nearest neighbor search on
     * point cloud */
    return WeakzoneKDTree_find_shortest_distance_single (tree, pt);
  }

  /* compute weak zone distance if the weak zones are given by functions */
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_SHELL:
    switch (physics_options->weakzone_type) {
    case SL_WEAKZONE_NONE:
      return 1.0;

    default: /* unknown weak zone type */
      YMIR_ABORT_NOT_REACHED ();
    }
    break;

  case SL_DOMAIN_BRICK:
  case SL_DOMAIN_SHELL_CHUNK:
  case SL_DOMAIN_SHELL_SLICE:
    switch (physics_options->weakzone_type) {
    case SL_WEAKZONE_NONE:
      return 1.0;

    case SL_WEAKZONE_2PLATES_POLY2:
      {
        double              subdu_lon;
        double              subdu_dip_angle;
        double              subdu_depth, subdu_width;
        double              subdu_thickness, subdu_thickness_const;
        double              total_thickness;

        double              ridge_depth, ridge_width;
        double              ridge_smoothwidth;
        double              lon_min = physics_options->domain_lon_min;

        double              courtesy_width;
        double              start_node, start_val, start_deriv;
        double              end_node, end_val;
        double              r, lon;
        double              d, dist = DBL_MAX;

        /* compute radius and longitude */
        r = slabs_compute_radius (x, y, z, physics_options);
        lon = slabs_compute_longitude (x, y, z, physics_options);

        /* set subduction parameters according to physics options */
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

        /* split by points in a rectangle containing the weak zone */
        if (   (  start_node - 0.5 * total_thickness
                / sin (subdu_dip_angle / 180.0 * M_PI) - courtesy_width )
               <= lon
            && lon <= (end_node + 0.5 * total_thickness + courtesy_width)
            && (end_val - total_thickness - courtesy_width) <= r ) {
          /* compute distance to subduction weak zone */
          d = slabs_2plates_weakzone_poly2_subdu_dist (
              r, lon, start_node, start_val, start_deriv, end_node, end_val);
          dist = SC_MIN (d / SL_EARTH_RADIUS, dist);
        }
        else {
          /* compute approx. distance to subduction weak zone */
          d = sqrt ( (start_node - lon) * (start_node - lon) +
                     (start_val - r) * (start_val - r) );
          dist = SC_MIN (d, dist);
          d = sqrt ( (end_node - lon) * (end_node - lon) +
                     (end_val - r) * (end_val - r) );
          dist = SC_MIN (d, dist);
        }

        /* set ridge parameters according to physics options */
        ridge_depth = physics_options->weakzone_2plates_ridge_depth;
        ridge_width = physics_options->weakzone_2plates_ridge_width;
        ridge_smoothwidth =
          physics_options->weakzone_2plates_ridge_smoothwidth;
        courtesy_width = 2.0 * ridge_smoothwidth / SL_EARTH_RADIUS;

        /* set bottom left corner of weak zone */
        end_node = lon_min + ridge_width / SL_EARTH_RADIUS
                           * SL_SHELL_RADIUS_TOP;
        end_val = SL_SHELL_RADIUS_TOP - ridge_depth / SL_EARTH_RADIUS;

        /* compute distance to ridge weak zone */
        d = slabs_2plates_weakzone_poly2_ridge_dist (
            r, lon, end_node, end_val);
        dist = SC_MIN (d / SL_EARTH_RADIUS, dist);

        /* return min distance to any weak zone */
        return dist;
      }

    default: /* unknown weak zone type */
      YMIR_ABORT_NOT_REACHED ();
    }
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }
}

/**
 *
 */
typedef struct slabs_physics_weak_dist_set_fn_data
{
  int                 n_fields;
  slabs_physics_options_t  *physics_options;
}
slabs_physics_weak_dist_set_fn_data_t;

/**
 *
 */
void
slabs_physics_weak_dist_set_fn (double *dist, double x, double y, double z,
                                ymir_locidx_t nid, void *data)
{
  slabs_physics_weak_dist_set_fn_data_t  *d =
    (slabs_physics_weak_dist_set_fn_data_t *) data;
  int                 fieldid;
  double              distance;

  distance = slabs_physics_weak_dist_node (x, y, z, d->physics_options);
  for (fieldid = 0; fieldid < d->n_fields; fieldid++) {
    dist[fieldid] = distance;
  }
}

/**
 *
 */
void
slabs_physics_weak_dist_set_face_fn (double *dist,
                                     double x, double y, double z,
                                     double nx, double ny, double nz,
                                     ymir_topidx_t face, ymir_locidx_t node_id,
                                     void *data)
{
  slabs_physics_weak_dist_set_fn_data_t  *d =
    (slabs_physics_weak_dist_set_fn_data_t *) data;
  int                 fieldid;
  double              distance;

  distance = slabs_physics_weak_dist_node (x, y, z, d->physics_options);

  for (fieldid = 0; fieldid < d->n_fields; fieldid++) {
    dist[fieldid] = distance;
  }
}

/**
 * Calculates viscosity for temperature dependent (only) viscosity law.
 *
 *   visc (T) = exp (E * (0.5 - T))
 */
 double
slabs_visc_temp_fn (const double visc_temp_decay, const double temp)
{
  return exp (visc_temp_decay * (0.5 - temp));
}

/**
 * Computes temperature dependent viscosity at a node.
 */
static inline double
slabs_visc_temp_node (const double temp, const double weak,
                      const double scaling, const double visc_temp_decay,
                      slabs_physics_options_t *physics_options,
                      const int restrict_to_bounds)
{
  const slabs_viscosity_model_t  viscosity_model =
                        physics_options->viscosity_model_type;
  const double        visc_min = physics_options->viscosity_min;
  const double        visc_max = physics_options->viscosity_max;
  const double        visc_temp_max = physics_options->viscosity_temp_max;

  double              visc_node;

  /* compute viscosity */
  visc_node = scaling * slabs_visc_temp_fn (visc_temp_decay, temp);

  /* apply weak zone and restrict to bounds */
  switch (viscosity_model) {
  case SL_VISCOSITY_MODEL_WYUL:
  case SL_VISCOSITY_MODEL_UWYUL:
    /* restrict viscosity to upper bound for temp. dep. viscosity */
    if (   viscosity_model == SL_VISCOSITY_MODEL_UWYUL && restrict_to_bounds
        && 0.0 < visc_temp_max && visc_temp_max < visc_node ) {
      visc_node = visc_temp_max;
    }

    /* multiply by weak zone */
    visc_node *= weak;

    /* restrict viscosity to upper bound */
    if ( restrict_to_bounds && 0.0 < visc_max && visc_max < visc_node ) {
      visc_node = visc_max;
    }

    /* restrict viscosity to lower bound */
    if ( restrict_to_bounds && 0.0 < visc_min && visc_node < visc_min ) {
      visc_node = visc_min;
    }
    break;

  case SL_VISCOSITY_MODEL_UWYL:
  case SL_VISCOSITY_MODEL_UYWL:
  case SL_VISCOSITY_MODEL_UYWL_SHIFT:
  case SL_VISCOSITY_MODEL_UWL_IIE_REG:
    /* restrict viscosity to upper bound */
    if ( restrict_to_bounds && 0.0 < visc_max && visc_max < visc_node ) {
      visc_node = visc_max;
    }

    /* multiply by weak zone */
    visc_node *= weak;

    /* restrict viscosity to lower bound */
    if ( restrict_to_bounds && 0.0 < visc_min && visc_node < visc_min ) {
      visc_node = visc_min;
    }
    break;

  case SL_VISCOSITY_MODEL_UWYL_LREG:
  case SL_VISCOSITY_MODEL_UWYL_SHIFT_LREG:
    /* restrict viscosity to upper bound */
    if ( restrict_to_bounds && 0.0 < visc_max && visc_max < visc_node ) {
      visc_node = visc_max;
    }

    /* multiply by weak zone */
    visc_node *= weak;

    /* restrict viscosity to lower bound */
    if ( restrict_to_bounds && 0.0 < visc_min ) {
      visc_node += visc_min;
    }
    break;

  default: /* unknown viscosity model type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* return temperature dependent viscosity */
  return visc_node;
}

/**
 * Computes temperature dependent viscosity in an element.
 */
void
slabs_visc_temp_elem (sc_dmatrix_t *visc_el_mat,
                      const double *x, const double *y, const double *z,
                      const int *Vmask,
                      const sc_dmatrix_t *temp_el_mat,
                      const sc_dmatrix_t *weak_el_mat,
                      slabs_physics_options_t *physics_options,
                      const int restrict_to_bounds)
{
  const double        upper_mantle_radius =
                        physics_options->viscosity_upper_mantle_radius;
  const double        transition_zone =
                        physics_options->viscosity_lower_upper_transition_zone;
  const int           n_nodes_per_el = visc_el_mat->m;
  const double       *temp_el_data = temp_el_mat->e[0];
  double             *visc_el_data = visc_el_mat->e[0];
  double              scaling, visc_temp_decay;
  int                 nodeid;

  /* check input parameters */
  YMIR_ASSERT (visc_el_mat->n == 1);
  YMIR_ASSERT (temp_el_mat != NULL);
  YMIR_ASSERT (temp_el_mat->m == n_nodes_per_el);
  YMIR_ASSERT (temp_el_mat->n == visc_el_mat->n);
  YMIR_ASSERT (weak_el_mat == NULL || weak_el_mat->m == n_nodes_per_el);
  YMIR_ASSERT (weak_el_mat == NULL || weak_el_mat->n == visc_el_mat->n);
  YMIR_ASSERT (sc_dmatrix_is_valid (temp_el_mat));

  /* set parameters depending on location in lower or upper mantle */
  if (slabs_physics_elem_in_upper_mantle (x, y, z, Vmask, physics_options)) {
    /* set upper mantle parameters */
    scaling = physics_options->viscosity_scaling;
    visc_temp_decay = physics_options->viscosity_temp_decay;
  }
  else {
    /* set lower mantle parameters */
    scaling = physics_options->viscosity_lower_mantle_scaling;
    visc_temp_decay = physics_options->viscosity_lower_mantle_temp_decay;
  }

  /* compute viscosity in this element */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) { /* loop over all
                                                         * nodes */
    const double        r =
      slabs_compute_radius (x[nodeid], y[nodeid], z[nodeid], physics_options);
    const double        temp = temp_el_data[nodeid];
    double              weak;

    /* check temperature for valid range [0,1] */
    YMIR_ASSERT (0.0 <= temp && temp <= 1.0);

    /* set weak zone */
    if (weak_el_mat != NULL) {
      weak = weak_el_mat->e[nodeid][0];
      YMIR_ASSERT (isfinite (weak));
      YMIR_ASSERT (0.0 < weak && weak <= 1.0);
    }
    else {
      weak = 1.0;
    }

    /* compute viscosity */
    if (upper_mantle_radius <= 0.0 || transition_zone <= 0.0 ||
        transition_zone < fabs (r - upper_mantle_radius)) {
      /* compute viscosity disregarding LM/UM transition zone */
      visc_el_data[nodeid] = slabs_visc_temp_node (temp, weak, scaling,
                                                   visc_temp_decay,
                                                   physics_options,
                                                   restrict_to_bounds);
    }
    else {
      const double        factor =
        (transition_zone + (r - upper_mantle_radius)) / (2.0*transition_zone);
      double              visc_lm, visc_um;

      /* compute viscosity with LM/UM transition zone (convex combination) */
      visc_lm = slabs_visc_temp_node (
          temp, weak, physics_options->viscosity_lower_mantle_scaling,
          physics_options->viscosity_lower_mantle_temp_decay,
          physics_options, restrict_to_bounds);
      visc_um = slabs_visc_temp_node (
          temp, weak, physics_options->viscosity_scaling,
          physics_options->viscosity_temp_decay,
          physics_options, restrict_to_bounds);
      visc_el_data[nodeid] = (1.0 - factor) * visc_lm + factor * visc_um;
    }

    /* check viscosity for `nan`, `inf`, and positivity */
    YMIR_ASSERT (isfinite (visc_el_data[nodeid]));
    YMIR_ASSERT (0.0 < visc_el_data[nodeid]);
  }
}

/**
 * Computes the temperature dependent viscosity vector.
 */
static void
slabs_viscosity_linear (ymir_dvec_t *viscosity,
                        slabs_stokes_state_t *state,
                        slabs_physics_options_t *physics_options)
{
  ymir_cvec_t        *temp_vec = state->temp_vec;
  ymir_dvec_t        *weak_vec = state->weak_vec;
  rhea_domain_options_t     domain_opt;
  rhea_viscosity_options_t  visc_opt;

  //TODO remove this copying of options
  domain_opt.shape = physics_options->domain_shape;
  domain_opt.lm_um_interface_radius =
    physics_options->viscosity_upper_mantle_radius;
  domain_opt.lm_um_interface_smoothing_width =
    physics_options->viscosity_lower_upper_transition_zone;
  domain_opt.x_min = physics_options->domain_x_min;
  domain_opt.x_max = physics_options->domain_x_max;
  domain_opt.y_min = physics_options->domain_y_min;
  domain_opt.y_max = physics_options->domain_y_max;
  domain_opt.z_min = physics_options->domain_z_min;
  domain_opt.z_max = physics_options->domain_z_max;
  domain_opt.lon_min = physics_options->domain_lon_min;
  domain_opt.lon_max = physics_options->domain_lon_max;
  domain_opt.radius_min = physics_options->domain_radius_min;
  domain_opt.radius_max = physics_options->domain_radius_max;
//domain_opt.volume = ;
//domain_opt.center = ;
//domain_opt.moment_of_inertia = ;

  visc_opt.type = physics_options->viscosity_type;
  visc_opt.type_init_nonlinear =
    physics_options->viscosity_type_for_init_nl_stokes;
  visc_opt.model = RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT;
  visc_opt.min = physics_options->viscosity_min;
  visc_opt.max = physics_options->viscosity_max;
  visc_opt.upper_mantle_scaling =
    physics_options->viscosity_scaling;
  visc_opt.upper_mantle_activation_energy =
    physics_options->viscosity_temp_decay;
  visc_opt.lower_mantle_scaling =
    physics_options->viscosity_lower_mantle_scaling;
  visc_opt.lower_mantle_activation_energy =
    physics_options->viscosity_lower_mantle_temp_decay;
  visc_opt.stress_exponent = physics_options->viscosity_stress_exponent;
  visc_opt.yield_stress = physics_options->viscosity_stress_yield;
  visc_opt.domain_options = &domain_opt;

  rhea_viscosity_compute (viscosity, NULL, NULL, NULL, temp_vec, weak_vec,
                          NULL, &visc_opt);
}

/**
 * Computes the nonlinear viscosity with left shift of the 2nd invariant.
 * Computes the derivative of the viscosity w.r.t. the 2nd invariant.
 */
static inline void
slabs_visc_nl_fn (double *viscosity,
                  double *dvisc_dIIe,
                  double *rank1_scal,
                  const double viscosity_temp,
                  const double IIe,
                  const double IIe_reg,
                  const double stress_exp)
{
  YMIR_ASSERT (0.0 < viscosity_temp);
  YMIR_ASSERT (0.0 <= IIe);
  YMIR_ASSERT (0.0 <= IIe_reg);
  YMIR_ASSERT (1.0 <= stress_exp);

  /* compute viscosity
   *
   *   viscosity = visc_temp * (IIe + IIe_reg)^((1 - n) / (2*n))
   */
  *viscosity = viscosity_temp * pow (IIe + IIe_reg,
                                     (1.0 - stress_exp) / (2.0 * stress_exp));

  /* compute derivative w.r.t. 2nd invariant
   *
   *   dvisc/dIIe = visc_temp * (1 - n) / (2*n)
   *                          * (IIe + IIe_reg)^((1 - 3*n) / (2*n))
   */
  *dvisc_dIIe = *viscosity * (1.0 - stress_exp) / (2.0 * stress_exp)
                           / (IIe + IIe_reg);

  /* compute rank-1 tensor scaling
   *
   *   rank1_scal = (1 - n) / n
   */
  *rank1_scal = (1.0 - stress_exp) / stress_exp;
}

/**
 * Computes the nonlinear viscosity with shift of the strain rate.
 * Computes the derivative of the viscosity w.r.t. the 2nd invariant.
 */
static inline void
slabs_visc_nl_shift_fn (double *viscosity,
                        double *dvisc_dIIe,
                        double *rank1_scal,
                        const double viscosity_temp,
                        const double strain_rate,
                        const double shift,
                        const double stress_exp)
{
  YMIR_ASSERT (0.0 < viscosity_temp);
  YMIR_ASSERT (0.0 <= strain_rate);
  YMIR_ASSERT (0.0 <= shift);
  YMIR_ASSERT (1.0 <= stress_exp);

  /* compute viscosity
   *
   *   viscosity = visc_temp * (strain_rate - shift)^(1/n) / strain_rate
   */
  *viscosity = viscosity_temp
               * pow (strain_rate - shift, 1.0 / stress_exp)
               / strain_rate;

  /* compute derivative w.r.t. 2nd invariant
   *
   *   dvisc/dIIe = visc_temp * (strain_rate - shift)^((1 - n) / n)
   *                          * (1/n - 1 + shift/strain_rate)
   *                          / (2 * strain_rate^2)
   */
  *dvisc_dIIe = *viscosity / (strain_rate - shift) / (2.0 * strain_rate)
                           * (1.0 / stress_exp - 1.0 + shift / strain_rate);

  /* compute rank-1 tensor scaling
   *
   *   rank1_scal = (strain_rate/n - strain_rate + shift)
   *                / (strain_rate - shift)
   *              = (1 - n) / n + shift / (n * (strain_rate - shift))
   */
  *rank1_scal = (1.0 - stress_exp) / stress_exp +
                shift / (stress_exp * (strain_rate - shift));
  *rank1_scal = SC_MIN ( SC_MAX (-1.0, *rank1_scal), 0.0);
}

/**
 * Computes nonlinear viscosity at a node.
 */
static inline void
slabs_visc_nl_node (double *viscosity, double *dvisc_dIIe, double *rank1_scal,
                    double *bounds_active, double *yielding_active,
                    const double temp, const double weak, const double IIe,
                    const double scaling, const double visc_temp_decay,
                    slabs_physics_options_t *physics_options)
{
  const slabs_viscosity_model_t  viscosity_model =
                        physics_options->viscosity_model_type;
  const double        stress_exp = physics_options->viscosity_stress_exponent;
  const double        visc_min = SC_MAX (0.0, physics_options->viscosity_min);
  const double        visc_max = SC_MAX (0.0, physics_options->viscosity_max);
  const double        stress_yield = physics_options->viscosity_stress_yield;
  const double        yield_reg =
    SC_MIN (SC_MAX (0.0, physics_options->viscosity_yielding_reg), 1.0);
  const double        strain_rate = sqrt (IIe);

  double              viscosity_temp =
                        scaling * slabs_visc_temp_fn (visc_temp_decay, temp);

  /* initialize bounds and yielding markers */
  *bounds_active = 0.0;
  *yielding_active = 0.0;

  switch (viscosity_model) {
  case SL_VISCOSITY_MODEL_WYUL: //TODO outdated
  case SL_VISCOSITY_MODEL_UWYUL: //TODO outdated
    /* Viscosity and its derivative w.r.t. the 2nd invariant of the strain rate
     * are computed, using the model
     *
     *   (1) upper bound for temperature dependent viscosity (UWYUL only),
     *   (2) weak zone,
     *   (3) yielding,
     *   (4) upper bound for full viscosity,
     *   (5) lower bound,
     *
     * as follows:
     *
     *   a_u = a_U (T) = min(visc_max, a (T))
     *
     *   visc (e) = w * a_U * e^((1 - n)/n)     , if e < e_yield
     *   visc (e) = stress_yield / (2 * e)
     *              + visc_yield * (1 - e_yield / e), if e_yield <= e
     *
     *   dvisc/dIIe (e) = w*a_U * (1-n)/(2*n) * e^((1 - 3*n)/n), if e < e_yield
     *   dvisc/dIIe (e) = (  visc_yield * e_yield / 2
     *                     - stress_yield / 4 ) / e^3          , if e_yield <= e
     *
     *   IF visc_max < visc (e)
     *     visc (e) = visc_max, dvisc/dIIe (e) = 0
     *
     *   IF visc (e) < visc_min
     *     visc (e) = visc_min, dvisc/dIIe (e) = 0
     *
     * where
     *
     *   e     --- strain rate,
     *   a (T) --- temperature dependent viscosity (prefactor),
     *   w     --- weak zone factor,
     *   n     --- stress exponent.
     */
    {
      const double        visc_temp_max = physics_options->viscosity_temp_max;
      double              strain_rate_yield;
      double              visc_yield;

      /* apply upper bound to temperature dependent viscosity */
      if (   viscosity_model == SL_VISCOSITY_MODEL_UWYUL
          && 0.0 < visc_temp_max && visc_temp_max < viscosity_temp ) {
        viscosity_temp = visc_temp_max;
        *bounds_active = 0.5;
      }

      /* set yielding parameters */
      if (0.0 < stress_yield) {
        /* compute strain rate at which yielding starts:
         *
         *   e_yield = (stress_yield / (2 * w * a))^n
         */
        strain_rate_yield = pow (stress_yield / (2.0 * weak * viscosity_temp),
                                 stress_exp);

        /* set viscosity for yielding */
        visc_yield = visc_min;
      }
      else { /* if no yielding */
        strain_rate_yield = DBL_MAX;
        visc_yield = 0.0;
      }

      /* compute viscosity and derivative w.r.t. 2nd invariant */
      if ( !(0.0 < stress_yield && strain_rate_yield <= strain_rate) ) {
        /* compute nonlinear viscosity */
        slabs_visc_nl_fn (viscosity, dvisc_dIIe, rank1_scal,
                          weak * viscosity_temp, IIe, 0.0, stress_exp);
      }
      else {
        /* compute viscosity when yielding */
        *viscosity =
          (0.5 * stress_yield / strain_rate)
          + visc_yield * SC_MAX (0.0, 1.0 - strain_rate_yield / strain_rate);
        *dvisc_dIIe =
          SC_MIN (0.0, 0.5*visc_yield*strain_rate_yield - 0.25*stress_yield)
          / (strain_rate * strain_rate * strain_rate);
        *yielding_active = 1.0;
      }

      /* apply upper bound to viscosity */
      YMIR_ASSERT (0.0 < visc_max);
      if (visc_max < *viscosity) {
        *viscosity = visc_max;
        *dvisc_dIIe = 0.0;
        *bounds_active = 1.0;
        *yielding_active = 0.0;
      }

      /* apply lower bound to viscosity */
      if (0.0 < visc_min && *viscosity < visc_min) {
        *viscosity = visc_min;
        *dvisc_dIIe = 0.0;
        *bounds_active = -1.0;
      }
    }
    break;

  case SL_VISCOSITY_MODEL_UWYL:
    /* Viscosity and its derivative w.r.t. the 2nd invariant of the strain rate
     * are computed, using the model
     *
     *   (1) upper bound, (2) weak zone, (3) yielding, (4) lower bound,
     *
     * as follows:
     *
     *   visc (e) = max(w*visc_max, visc_min), if e < min(e_min, e_yield)
     *   visc (e) = w*a * e^((1 - n)/n)    , if e_min <= e < min(e_yield, e_max)
     *   visc (e) = stress_yield / (2 * e)
     *              + visc_yield
     *              * (1 - e_yield / e)    , if e_yield <= e < e_max
     *   visc (e) = visc_min               , if e_max <= e
     *
     *   dvisc/dIIe (e) = 0                , if e < min(e_min, e_yield)
     *   dvisc/dIIe (e) = w*a * (1-n)/(2*n)
     *                    * e^((1 - 3*n)/n), if e_min <= e < min(e_yield, e_max)
     *   dvisc/dIIe (e) = (visc_yield*e_yield/2
     *                     - stress_yield/4)
     *                    / e^3            , if e_yield <= e < e_max
     *   dvisc/dIIe (e) = 0                , if e_max <= e
     *
     * where
     *
     *   e --- strain rate,
     *   a --- prefactor,
     *   w --- weak zone factor,
     *   n --- stress exponent.
     */
    {
      /* compute nonlinear viscosity */
      slabs_visc_nl_fn (viscosity, dvisc_dIIe, rank1_scal,
                        viscosity_temp, IIe, 0.0, stress_exp);

      /* apply upper bound to viscosity */
      YMIR_ASSERT (0.0 < visc_max);
      if (visc_max < *viscosity) {
        *viscosity = visc_max;
        *dvisc_dIIe = 0.0;
        *rank1_scal = 0.0;
        *bounds_active = 1.0;
      }

      /* multiply in weak factor */
      *viscosity *= weak;
      *dvisc_dIIe *= weak;

      /* apply yielding */
      if (0.0 < stress_yield) {
        double              strain_rate_yield;
        double              visc_yield;

        /* compute strain rate at which yielding starts:
         *
         *   e_yield = (stress_yield / (2 * w * a))^n
         */
        strain_rate_yield = pow (stress_yield / (2.0 * weak * viscosity_temp),
                                 stress_exp);

        /* set viscosity for yielding */
        visc_yield = visc_min;

        /* compute viscosity when yielding */
        if (strain_rate_yield <= strain_rate) {
          double              v, vy, dv, dvy, r1, r1y;

          v = *viscosity;
          vy = (0.5 * stress_yield / strain_rate) +
               visc_yield * SC_MAX (0.0, 1.0 - strain_rate_yield / strain_rate);
          dv = *dvisc_dIIe;
          dvy = SC_MIN ( 0.5*visc_yield*strain_rate_yield - 0.25*stress_yield,
                         0.0 ) / (strain_rate * strain_rate * strain_rate);
          r1 = *rank1_scal;
          r1y = -1.0;

          if (0.0 < yield_reg) {
            *viscosity = yield_reg * v + (1.0 - yield_reg) * vy;
            *dvisc_dIIe = yield_reg * dv + (1.0 - yield_reg) * dvy;
            *rank1_scal = (yield_reg * v * r1 + (1.0 - yield_reg) * vy * r1y) /
                          *viscosity;
          }
          else {
            *viscosity = vy;
            *dvisc_dIIe = dvy;
            *rank1_scal = r1y;
          }
          *yielding_active = 1.0;
        }
      }

      /* apply upper bound to viscosity (to enforce bound when yielding) */
      if ((weak * visc_max) < *viscosity) {
        *viscosity = weak * visc_max;
        *dvisc_dIIe = 0.0;
        *rank1_scal = 0.0;
        *bounds_active = 1.0;
        *yielding_active = 0.0;
      }

      /* change upper bound marker if weak zone is present */
      if (fabs (*bounds_active - 1.0) < SC_EPS && weak < 0.5) {
        *bounds_active = 0.5;
      }

      /* apply lower bound to viscosity */
      if (0.0 < visc_min && *viscosity < visc_min) {
        *viscosity = visc_min;
        *dvisc_dIIe = 0.0;
        *rank1_scal = 0.0;
        *bounds_active = -1.0;
      }
    }
    break;

  case SL_VISCOSITY_MODEL_UWYL_LREG:
    /* Viscosity and its derivative w.r.t. the 2nd invariant of the strain rate
     * are computed, using the model
     *
     *   (1) upper bound, (2) weak zone, (3) yielding, (4) lower bound,
     *
     * as follows:
     *
     *   visc (e) = w * a * e^((1 - n) / n)           , if e < e_yield
     *   visc (e) = (stress_yield / 2) / e            , if e_yield <= e
     *
     *   dvisc/dIIe (e) = w * a * ((1 - n) / n)
     *                    * e^((1 - n) / n)
     *                    / (2 * e^2)                 , if e < e_yield
     *   dvisc/dIIe (e) = - (stress_yield / 2) / e
     *                    / (2 * e^2)                 , if e_yield <= e
     *
     * where
     *
     *   e --- strain rate,
     *   a --- prefactor,
     *   w --- weak zone factor,
     *   n --- stress exponent.
     */
    {
      /* compute nonlinear viscosity */
      slabs_visc_nl_fn (viscosity, dvisc_dIIe, rank1_scal,
                        viscosity_temp, IIe, 0.0, stress_exp);

      /* apply upper bound to viscosity */
      YMIR_ASSERT (0.0 < visc_max);
      if (visc_max < *viscosity) {
        *viscosity = visc_max;
        *dvisc_dIIe = 0.0;
        *rank1_scal = 0.0;
        *bounds_active = 1.0;
      }

      /* multiply in weak factor */
      *viscosity *= weak;
      *dvisc_dIIe *= weak;

      /* change upper bound marker if weak zone is present */
      if (fabs (*bounds_active - 1.0) < SC_EPS && weak < 0.5) {
        *bounds_active = 0.5;
      }

      /* apply yielding */
      if ( 0.0 < stress_yield &&
           stress_yield < (2.0 * (*viscosity) * strain_rate) ) {
        double              v, vy, dv, dvy, r1, r1y;

        v = *viscosity;
        vy = (stress_yield / 2.0) / strain_rate;
        dv = *dvisc_dIIe;
        dvy = - (stress_yield / 2.0) / strain_rate /
                (2.0 * strain_rate * strain_rate);
        r1 = *rank1_scal;
        r1y = -1.0;

        if (0.0 < yield_reg) {
          *viscosity = yield_reg * v + (1.0 - yield_reg) * vy;
          *dvisc_dIIe = yield_reg * dv + (1.0 - yield_reg) * dvy;
          *rank1_scal = (yield_reg * v * r1 + (1.0 - yield_reg) * vy * r1y) /
                        *viscosity;
        }
        else {
          *viscosity = vy;
          *dvisc_dIIe = dvy;
          *rank1_scal = r1y;
        }
        *yielding_active = 1.0;
      }

      /* apply lower bound to viscosity */
      if (0.0 < visc_min) {
        *viscosity += visc_min;
        YMIR_ASSERT (
            fabs (ymir_nlstress_op_coeff_tensor_add + 2.0 * visc_min) < SC_EPS);
      }
    }
    break;

  case SL_VISCOSITY_MODEL_UWYL_SHIFT_LREG:
    /* Viscosity and its derivative w.r.t. the 2nd invariant of the strain rate
     * are computed, using the model
     *
     *   (1) upper bound + shift, (2) weak zone, (3) yielding, (4) lower bound,
     *
     * as follows:
     *
     *   visc (e) = w * a * (e - d)^(1/n) / e         , if e < e_yield
     *   visc (e) = (stress_yield / 2) / e            , if e_yield <= e
     *
     *   dvisc/dIIe (e) = w * a * (e - d)^((1 - n)/n)
     *                    * (1 / n - 1 + d / e) / e
     *                    / (2 * e^2)                 , if e < e_yield
     *   dvisc/dIIe (e) = - (stress_yield / 2) / e
     *                    / (2 * e^2)                 , if e_yield <= e
     *
     * where
     *
     *   e --- strain rate,
     *   a --- prefactor,
     *   w --- weak zone factor,
     *   n --- stress exponent,
     *   d --- shift (greater zero).
     */
    {
      double              strain_rate_min;
      double              shift;

      /* compute minimal strain rate for nonlinear viscosity
       *
       *   e_min = a / visc_max * (visc_max * n / a)^(1 / (1 - n))
       */
      YMIR_ASSERT (0.0 < visc_max);
      strain_rate_min = viscosity_temp / visc_max
                        * pow ( visc_max * stress_exp / viscosity_temp,
                                1.0 / (1.0 - stress_exp) );

      /* compute shift of strain rate
       *
       *   d = e_min - (visc_max * n / a)^(n / (1 - n))
       */
      shift = strain_rate_min
              - pow ( visc_max * stress_exp / viscosity_temp,
                      stress_exp / (1.0 - stress_exp) );
      shift = SC_MAX (shift, 0.0);

      /* compute nonlinear viscosity */
      slabs_visc_nl_shift_fn (viscosity, dvisc_dIIe, rank1_scal,
                              viscosity_temp, strain_rate, shift, stress_exp);

      /* apply upper bound */
      if (strain_rate < strain_rate_min) {
        *viscosity = visc_max;
        *dvisc_dIIe = 0.0;
        *rank1_scal = 0.0;
        *bounds_active = 1.0;
      }

      /* multiply in weak factor */
      *viscosity *= weak;
      *dvisc_dIIe *= weak;

      /* change upper bound marker if weak zone is present */
      if (fabs (*bounds_active - 1.0) < SC_EPS && weak < 0.5) {
        *bounds_active = 0.5;
      }

      /* apply yielding */
      if ( 0.0 < stress_yield &&
           stress_yield < (2.0 * (*viscosity) * strain_rate) ) {
        double              v, vy, dv, dvy, r1, r1y;

        v = *viscosity;
        vy = (stress_yield / 2.0) / strain_rate;
        dv = *dvisc_dIIe;
        dvy = - (stress_yield / 2.0) / strain_rate /
                (2.0 * strain_rate * strain_rate);
        r1 = *rank1_scal;
        r1y = -1.0;

        if (0.0 < yield_reg) {
          *viscosity = yield_reg * v + (1.0 - yield_reg) * vy;
          *dvisc_dIIe = yield_reg * dv + (1.0 - yield_reg) * dvy;
          *rank1_scal = (yield_reg * v * r1 + (1.0 - yield_reg) * vy * r1y) /
                        *viscosity;
        }
        else {
          *viscosity = vy;
          *dvisc_dIIe = dvy;
          *rank1_scal = r1y;
        }
        *yielding_active = 1.0;
      }

      /* apply lower bound to viscosity */
      if (0.0 < visc_min) {
        *viscosity += visc_min;
        YMIR_ASSERT (
            fabs (ymir_nlstress_op_coeff_tensor_add + 2.0 * visc_min) < SC_EPS);
      }
    }
    break;

  case SL_VISCOSITY_MODEL_UYWL: //TODO outdated
    /* Viscosity and its derivative w.r.t. the 2nd invariant of the strain rate
     * are computed, using the model
     *
     *   (1) upper bound, (2) yielding, (3) weak zone, (4) lower bound,
     *
     * as follows:
     *
     *   visc (e) = max(w * visc_max,
     *                  visc_min)          , if e < min(e_min, e_yield, e_max)
     *   visc (e) = w*a * e^((1 - n)/n)    , if e_min <= e < min(e_yield, e_max)
     *   visc (e) = w*(visc_yield
     *                 + b_y/(2*e))        , if e_yield <= e < e_max
     *   visc (e) = visc_min               , if e_max <= e
     *
     *   dvisc/dIIe (e) = 0                , if e < min(e_min, e_yield, e_max)
     *   dvisc/dIIe (e) = w*a * (1 - n)/(2*n)
     *                    * e^((1 - 3*n)/n), if e_min <= e < min(e_yield, e_max)
     *   dvisc/dIIe (e) = -w*b_y / 4 / e^3 , if e_yield <= e < e_max
     *   dvisc/dIIe (e) = 0                , if e_max <= e
     *
     * where
     *
     *   e   --- strain rate,
     *   a   --- prefactor,
     *   w   --- weak zone factor,
     *   n   --- stress exponent,
     *   b   --- lift of linear viscosity for lower bound,
     *   b_y --- lift of linear viscosity for yielding.
     */
    {
      double              strain_rate_min, strain_rate_yield, strain_rate_max;
      double              visc_yield;
      double              lift_yield;

      /*
       * set parameters
       */

      /* set upper bound parameters */
      YMIR_ASSERT (0.0 < visc_max);
      YMIR_ASSERT (visc_min < visc_max);

      /* compute minimal strain rate, where rheology is non-Newtoninan:
       *
       *   e_min = (visc_max / a)^(n / (1 - n))
       */
      strain_rate_min = pow (visc_max / viscosity_temp,
                             stress_exp / (1.0 - stress_exp));

      /* set yielding parameters */
      if (0.0 < stress_yield) {
        /* compute strain rate at which yielding starts:
         *
         *   e_yield = (stress_yield / (2*a))^n
         */
        strain_rate_yield = pow (stress_yield / (2.0 * viscosity_temp),
                                 stress_exp);

        /* adjust if yielding within upper bound
         *
         *   e_yield = stress_yield / (2*visc_max), if e_yield < e_min
         */
        if (strain_rate_yield < strain_rate_min) {
          strain_rate_yield = stress_yield / (2.0 * visc_max);
        }

        /* set viscosity for yielding */
        visc_yield = visc_min / weak;

        /* compute lift for linear rheology at strain rates where yielding
         * occurs (require to be non-negative) */
        lift_yield = SC_MAX (0.0, stress_yield
                                  - 2.0 * visc_yield * strain_rate_yield);
      }
      else { /* if no yielding */
        strain_rate_yield = DBL_MAX;
        visc_yield = 0.0;
        lift_yield = 0.0;
      }

      /* set lower bound parameters */
      if (0.0 < visc_min) {
        /* compute maximal strain rate, where rheology is non-Newtoninan:
         *
         *   e_max = (visc_min / w*a)^(n / (1 - n))
         */
        strain_rate_max = pow (visc_min / (weak * viscosity_temp),
                               stress_exp / (1.0 - stress_exp));

        /* modify if lower bound is after yielding */
        if (strain_rate_yield < strain_rate_max) {
          strain_rate_max = DBL_MAX;
        }

        /* modify min strain rate if lower bound is before upper bound */
        if (strain_rate_max <= strain_rate_min) {
          strain_rate_min = DBL_MAX;
        }
      }
      else { /* if no lower bound for viscosity */
        strain_rate_max = DBL_MAX;
      }

      /*
       * compute viscosity and its derivative
       */

      if (  strain_rate
          < SC_MIN (strain_rate_min,
                    SC_MIN (strain_rate_yield, strain_rate_max)) ) {
        /* if linear rheology for low strain rates, corresponding to
         * upper bound */
        *viscosity = weak * visc_max;
        *dvisc_dIIe = 0.0;
        *bounds_active = 1.0;

        /* enforce min viscosity bound if enabled */
        if (0.0 < visc_min && *viscosity < visc_min) {
          *viscosity = visc_min;
          *bounds_active = -1.0;
        }
      }
      else if (   strain_rate_min <= strain_rate
               && strain_rate < SC_MIN (strain_rate_yield, strain_rate_max) ) {
        /* if nonlinear rheology */
        slabs_visc_nl_fn (viscosity, dvisc_dIIe, rank1_scal,
                          weak * viscosity_temp, IIe, 0.0, stress_exp);
      }
      else if (   strain_rate_yield <= strain_rate
               && strain_rate < strain_rate_max ) {
        /* if rheology caused by yielding */
        *viscosity = weak * (visc_yield + lift_yield / (2.0 * strain_rate));
        *dvisc_dIIe = - (weak * lift_yield / 4.0)
                      / (strain_rate * strain_rate * strain_rate);
        *yielding_active = 1.0;
      }
      else if (strain_rate_max <= strain_rate) {
        /* if linear rheology for high strain rates corresponding to
         * lower bound */
        *viscosity = visc_min;
        *dvisc_dIIe = 0.0;
        *bounds_active = -1.0;
      }
      else { /* if strain rate is not valid */
        YMIR_ABORT_NOT_REACHED ();
      }
    }
    break;

  case SL_VISCOSITY_MODEL_UYWL_SHIFT: //TODO outdated
    /* Viscosity and its derivative w.r.t. the 2nd invariant of the strain rate
     * are computed, using the model
     *
     *   (1) upper bound, (2) yielding, (3) weak zone, (4) lower bound,
     *
     * as follows:
     *
     *   visc (e) = max(w * visc_max,
     *                  visc_min)          , if e < min(e_min, e_yield, e_max)
     *   visc (e) = w*a / e * (e - d)^(1/n), if e_min <= e < min(e_yield, e_max)
     *   visc (e) = w*(visc_yield
     *                 + b_y/(2*e))        , if e_yield <= e < e_max
     *   visc (e) = visc_min + b / e       , if e_max <= e
     *
     *   dvisc/dIIe (e) = 0                , if e < min(e_min, e_yield, e_max)
     *   dvisc/dIIe (e) = w*a / 2 / e^2
     *                    * (e - d)^((1 - n)/n)
     *                    * (1/n - 1 + d/e), if e_min <= e < min(e_yield, e_max)
     *   dvisc/dIIe (e) = -w*b_y / 4 / e^3 , if e_yield <= e < e_max
     *   dvisc/dIIe (e) = -b / 4 / e^3     , if e_max <= e
     *
     * where
     *
     *   e   --- strain rate,
     *   a   --- prefactor,
     *   w   --- weak zone factor,
     *   n   --- stress exponent,
     *   d   --- right-shift of strain rate,
     *   b   --- lift of linear viscosity for lower bound,
     *   b_y --- lift of linear viscosity for yielding.
     */
    {
      double              rshift;
      double              strain_rate_min, strain_rate_yield, strain_rate_max;
      double              visc_yield;
      double              lift_yield, lift_max;

      /*
       * set parameters
       */

      /* set upper bound parameters */
      YMIR_ASSERT (0.0 < visc_max);
      YMIR_ASSERT (visc_min < visc_max);

      /* compute minimal strain rate, where rheology is non-Newtonian:
       *
       *   e_min = (2*a / visc_max) * (n * visc_max / (2*a))^(1 / (1 - n))
       */
      strain_rate_min = (2.0 * viscosity_temp / visc_max)
                        * pow (stress_exp * visc_max / (2.0 * viscosity_temp),
                               1.0 / (1.0 - stress_exp));

      /* compute right shift:
       *
       *   d = e_min - (n * visc_max / (2*a))^(n / (1 - n)
       */
      rshift = strain_rate_min
               - pow (stress_exp * visc_max / (2.0 * viscosity_temp),
                      stress_exp / (1.0 - stress_exp));

      /* set yielding parameters */
      if (0.0 < stress_yield) {
        /* compute strain rate at which yielding starts:
         *
         *   e_yield = (stress_yield / (2*a))^n + d
         */
        strain_rate_yield = pow (stress_yield / (2.0 * viscosity_temp),
                                 stress_exp)
                            + rshift;

        /* adjust if yielding within upper bound
         *
         *   e_yield = stress_yield / (2*visc_max), if e_yield < e_min
         */
        if (strain_rate_yield < strain_rate_min) {
          strain_rate_yield = stress_yield / (2.0 * visc_max);
        }

        /* set viscosity for yielding */
        visc_yield = visc_min / weak;

        /* compute lift for linear rheology at strain rates where yielding
         * occurs (require to be non-negative) */
        lift_yield = SC_MAX (0.0,
                             stress_yield - 2.0*visc_yield*strain_rate_yield);
      }
      else { /* if no yielding */
        strain_rate_yield = DBL_MAX;
        visc_yield = 0.0;
        lift_yield = 0.0;
      }

      /* set lower bound parameters */
      if (0.0 < visc_min) {
        /* compute maximal strain rate, where rheology is non-Newtonian:
         *
         *   e_max = (n * visc_min / (2 * w * a))^(n / (1 - n)) + d
         */
        strain_rate_max = pow (stress_exp*visc_min / (2.0*weak*viscosity_temp),
                               stress_exp / (1.0 - stress_exp)) + rshift;

        /* set lift (require to be non-negative) */
        lift_max = SC_MAX ( 0.0,
                            2.0 * weak * viscosity_temp
                            * pow (strain_rate_max - rshift, 1.0 / stress_exp)
                            - 2.0 * visc_min * strain_rate_max );

        /* modify if lower bound is after yielding */
        if (strain_rate_yield < strain_rate_max) {
          strain_rate_max = DBL_MAX;
        }

        /* modify min strain rate if lower bound is before upper bound */
        if (strain_rate_max <= strain_rate_min) {
          strain_rate_min = DBL_MAX;
        }
      }
      else { /* if no lower bound for viscosity */
        strain_rate_max = DBL_MAX;
        lift_max = 0.0;
      }

      /*
       * compute viscosity and its derivative
       */

      if (  strain_rate
          < SC_MIN (strain_rate_min,
                    SC_MIN (strain_rate_yield, strain_rate_max)) ) {
        /* if linear rheology for low strain rates, corresponding to
         * upper bound */
        *viscosity = weak * visc_max;
        *dvisc_dIIe = 0.0;
        *bounds_active = 1.0;

        /* enforce min viscosity bound if enabled */
        if (0.0 < visc_min && *viscosity < visc_min) {
          *viscosity = visc_min;
          *bounds_active = -1.0;
        }
      }
      else if (   strain_rate_min <= strain_rate
               && strain_rate < SC_MIN (strain_rate_yield, strain_rate_max) ) {
        /* if nonlinear rheology */
        slabs_visc_nl_shift_fn (viscosity, dvisc_dIIe, rank1_scal,
                                weak * viscosity_temp, strain_rate, rshift,
                                stress_exp);
      }
      else if (   strain_rate_yield <= strain_rate
               && strain_rate < strain_rate_max ) {
        /* if rheology caused by yielding */
        *viscosity = weak * (visc_yield + lift_yield / (2.0 * strain_rate));
        *dvisc_dIIe = - (weak * lift_yield / 4.0)
                      / (strain_rate * strain_rate * strain_rate);
        *yielding_active = 1.0;
      }
      else if (strain_rate_max <= strain_rate) {
        /* if linear rheology for high strain rates corresponding to
         * lower bound */
        *viscosity = visc_min + lift_max / (2.0 * strain_rate);
        *dvisc_dIIe = - lift_max / 4.0
                      / (strain_rate * strain_rate * strain_rate);
        *bounds_active = -1.0;
      }
      else { /* if strain rate is not valid */
        YMIR_ABORT_NOT_REACHED ();
      }
    }
    break;

  case SL_VISCOSITY_MODEL_UWL_IIE_REG: //TODO outdated
    /* Viscosity and its derivative w.r.t. the 2nd invariant of the strain rate
     * are computed, using the model
     *
     *   (1) upper bound, (2) weak zone, (3) lower bound,
     *
     * as follows:
     *
     *   visc (IIe) = w*a * (IIe + d)^((1 - n)/(2 * n))
     *
     *   dvisc/dIIe (IIe) = w*a * (1 - n)/(2*n) * (IIe + d)^((1 - 3*n)/(2*n))
     *
     * where
     *
     *   IIe --- 2nd invariant of the strain rate,
     *   a   --- prefactor,
     *   w   --- weak zone factor,
     *   n   --- stress exponent,
     *   d   --- (left) shift of the 2nd invariant.
     */
    {
      const double        IIe_reg =
                            physics_options->viscosity_IIe_regularization;

      YMIR_ASSERT (stress_yield <= 0.0); /* because it's not implemented */

      /* compute the nonlinear viscosity and the derivative of the viscosity
       * w.r.t. the second invariant */
      slabs_visc_nl_fn (viscosity, dvisc_dIIe, rank1_scal,
                        viscosity_temp, IIe, IIe_reg, stress_exp);

      /* bound viscosity from above (and set derivative w.r.t. IIe to zero) */
      if (0.0 < visc_max && visc_max < *viscosity) {
        *viscosity = visc_max;
        *dvisc_dIIe = 0.0;
        *bounds_active = 1.0;
      }

      /* multiply viscosity and derivative by weak zone */
      *viscosity *= weak;
      *dvisc_dIIe *= weak;

      /* bound viscosity from below (and set derivative w.r.t. IIe to zero) */
      if (0.0 < visc_min && *viscosity < visc_min) {
        *viscosity = visc_min;
        *dvisc_dIIe = 0.0;
        *bounds_active = -1.0;
      }
    }
    break;

  default: /* unknown viscosity model type */
    YMIR_ABORT_NOT_REACHED ();
  }
}

static inline void
slabs_visc_nl_node_prescribe_bounds_yielding (double *viscosity,
                                              double *dvisc_dIIe,
                                              const double *bounds_active,
                                              const double *yielding_active,
                                              const double temp,
                                              const double weak,
                                              const double IIe,
                                              const double scaling,
                                              const double visc_temp_decay,
                                              slabs_physics_options_t
                                                *physics_options)
{
  const slabs_viscosity_model_t  viscosity_model =
                        physics_options->viscosity_model_type;
  const double        stress_exp = physics_options->viscosity_stress_exponent;
  const double        visc_min = SC_MAX (0.0, physics_options->viscosity_min);
  const double        visc_max = SC_MAX (0.0, physics_options->viscosity_max);
  const double        stress_yield = physics_options->viscosity_stress_yield;
  const double        strain_rate = sqrt (IIe);

  double              viscosity_temp =
                        scaling * slabs_visc_temp_fn (visc_temp_decay, temp);

  double              dummy; //TODO update s.t. dummy is not necessary

  switch (viscosity_model) {
  //case SL_VISCOSITY_MODEL_WYUL: TODO
  //case SL_VISCOSITY_MODEL_UWYUL: TODO

  case SL_VISCOSITY_MODEL_UWYL: //TODO outdated
    {
      /* compute nonlinear viscosity */
      slabs_visc_nl_fn (viscosity, dvisc_dIIe, &dummy, viscosity_temp,
                        IIe, 0.0, stress_exp);

      /* apply upper bound to viscosity */
      YMIR_ASSERT (0.0 < visc_max);
      *viscosity = SC_MIN (visc_max, *viscosity);
      if (0.0 < *bounds_active) {
        YMIR_ASSERT (*bounds_active <= 1.0);
        *dvisc_dIIe *= (1.0 - *bounds_active);
      }

      /* multiply in weak factor */
      *viscosity *= weak;
      *dvisc_dIIe *= weak;

      /* apply yielding */
      if (0.0 < *yielding_active) {
        double              strain_rate_yield;
        double              visc_yield;

        YMIR_ASSERT (*yielding_active <= 1.0);

        /* compute strain rate at which yielding starts:
         *
         *   e_yield = (stress_yield / (2 * w * a))^n
         */
        strain_rate_yield = pow (stress_yield / (2.0 * weak * viscosity_temp),
                                 stress_exp);

        /* set viscosity for yielding */
        visc_yield = visc_min;

        /* compute viscosity when yielding */
        *viscosity =
          *viscosity * (1.0 - *yielding_active) +
          *yielding_active * (0.5 * stress_yield / strain_rate) +
          visc_yield * SC_MAX (0.0, 1.0 - strain_rate_yield / strain_rate);
        *dvisc_dIIe =
          *dvisc_dIIe * (1.0 - *yielding_active) -
          *yielding_active *
          SC_MIN (0.0, 0.5*visc_yield*strain_rate_yield - 0.25*stress_yield)
          / (strain_rate * strain_rate * strain_rate);
      }

      /* apply lower bound to viscosity */
      if (0.0 < visc_min) {
        *viscosity = SC_MAX (visc_min, *viscosity);
        if (*bounds_active < 0.0) {
          YMIR_ASSERT (-1.0 <= *bounds_active);
          *dvisc_dIIe *= (1.0 + *bounds_active);
        }
      }
    }
    break;

  case SL_VISCOSITY_MODEL_UWYL_LREG: //TODO outdated
  case SL_VISCOSITY_MODEL_UWYL_SHIFT_LREG: //TODO outdated
    {
      /* compute nonlinear viscosity */
      slabs_visc_nl_fn (viscosity, dvisc_dIIe, &dummy, viscosity_temp,
                        IIe, 0.0, stress_exp);

      /* apply upper bound to viscosity */
      YMIR_ASSERT (0.0 < visc_max);
      *viscosity = SC_MIN (visc_max, *viscosity);
      if (0.0 < *bounds_active) {
        YMIR_ASSERT (*bounds_active <= 1.0);
        *dvisc_dIIe *= (1.0 - *bounds_active);
      }

      /* multiply in weak factor */
      *viscosity *= weak;
      *dvisc_dIIe *= weak;

      /* apply yielding */
      if (0.0 < *yielding_active) {
        YMIR_ASSERT (*yielding_active <= 1.0);
        *viscosity = *viscosity * (1.0 - *yielding_active) +
                     *yielding_active * (stress_yield / 2.0) / strain_rate;
        *dvisc_dIIe = *dvisc_dIIe * (1.0 - *yielding_active) -
                      *yielding_active * (stress_yield / 2.0) / strain_rate
                      / (2.0 * strain_rate * strain_rate);
      }

      /* apply lower bound to viscosity */
      if (0.0 < visc_min) {
        *viscosity += visc_min;
      }
    }
    break;

  //case SL_VISCOSITY_MODEL_UYWL: TODO
  //case SL_VISCOSITY_MODEL_UYWL_SHIFT: TODO
  //case SL_VISCOSITY_MODEL_UWL_IIE_REG: TODO

  default: /* unknown viscosity model type */
    YMIR_ABORT_NOT_REACHED ();
  }
}

/**
 * Computes nonlinear viscosity in an element.
 */
void
slabs_visc_nl_elem (sc_dmatrix_t *visc_el_mat,
                    sc_dmatrix_t *dvisc_dIIe_el_mat,
                    sc_dmatrix_t *rank1_scal_el_mat,
                    sc_dmatrix_t *bounds_el_mat,
                    sc_dmatrix_t *yielding_el_mat,
                    const double *x, const double *y, const double *z,
                    const int *Vmask,
                    const sc_dmatrix_t *temp_el_mat,
                    const sc_dmatrix_t *weak_el_mat,
                    const sc_dmatrix_t *IIe_el_mat,
                    slabs_physics_options_t *physics_options,
                    const int prescribe_bounds_yielding)
{
  const int           n_nodes_per_el = visc_el_mat->m;
  const double       *temp_el_data = temp_el_mat->e[0];
  const double       *weak_el_data = weak_el_mat->e[0];
  const double       *IIe_el_data = IIe_el_mat->e[0];
  double             *visc_el_data = visc_el_mat->e[0];
  double             *dvisc_dIIe_el_data = dvisc_dIIe_el_mat->e[0];
  double             *rank1_scal_el_data;
  double             *bounds_el_data;
  double             *yielding_el_data;
  int                 nodeid;

  /* check input */
  YMIR_ASSERT (visc_el_mat->n == 1);
  YMIR_ASSERT (dvisc_dIIe_el_mat->m == n_nodes_per_el);
  YMIR_ASSERT (dvisc_dIIe_el_mat->n == visc_el_mat->n);
  YMIR_ASSERT (rank1_scal_el_mat == NULL ||
               rank1_scal_el_mat->m == n_nodes_per_el);
  YMIR_ASSERT (rank1_scal_el_mat == NULL ||
               rank1_scal_el_mat->n == visc_el_mat->n);
  YMIR_ASSERT (bounds_el_mat == NULL || bounds_el_mat->m == n_nodes_per_el);
  YMIR_ASSERT (bounds_el_mat == NULL || bounds_el_mat->n == visc_el_mat->n);
  YMIR_ASSERT (yielding_el_mat == NULL || yielding_el_mat->m == n_nodes_per_el);
  YMIR_ASSERT (yielding_el_mat == NULL || yielding_el_mat->n == visc_el_mat->n);
  YMIR_ASSERT (temp_el_mat->m == n_nodes_per_el);
  YMIR_ASSERT (temp_el_mat->n == visc_el_mat->n);
  YMIR_ASSERT (weak_el_mat->m == n_nodes_per_el);
  YMIR_ASSERT (weak_el_mat->n == visc_el_mat->n);
  YMIR_ASSERT (IIe_el_mat->m == n_nodes_per_el);
  YMIR_ASSERT (IIe_el_mat->n == visc_el_mat->n);
  YMIR_ASSERT (sc_dmatrix_is_valid (temp_el_mat));
  YMIR_ASSERT (sc_dmatrix_is_valid (weak_el_mat));
  YMIR_ASSERT (sc_dmatrix_is_valid (IIe_el_mat));

  /* get pointers to entries of rank-1 tensor scaling */
  if (rank1_scal_el_mat != NULL) {
    rank1_scal_el_data = rank1_scal_el_mat->e[0];
  }
  else {
    rank1_scal_el_data = NULL;
  }

  /* get pointers to entries of bounds and yielding markers */
  if (bounds_el_mat != NULL) {
    bounds_el_data = bounds_el_mat->e[0];
  }
  else {
    bounds_el_data = NULL;
  }
  if (yielding_el_mat != NULL) {
    yielding_el_data = yielding_el_mat->e[0];
  }
  else {
    yielding_el_data = NULL;
  }

  /*
   * Compute Viscosity
   */

  /* compute viscosity depending on location in lower or upper mantle */
  if (slabs_physics_elem_in_upper_mantle (x, y, z, Vmask, physics_options)) {
    const double        scaling = physics_options->viscosity_scaling;
    const double        visc_temp_decay =
                          physics_options->viscosity_temp_decay;
    double             *rank1_scal, *bounds, *yielding, dummy;

    for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) { /* loop over all
                                                           * nodes */
      const double        temp = temp_el_data[nodeid];
      const double        weak = weak_el_data[nodeid];
      const double        IIe = IIe_el_data[nodeid];

      /* check temperature for valid range [0,1] */
      YMIR_ASSERT (0.0 <= temp && temp <= 1.0);

      /* check weak zone for valid range (0,1] */
      YMIR_ASSERT (0.0 < weak && weak <= 1.0);

      /* check second invariant for non-negativity */
      YMIR_ASSERT (0.0 <= IIe);

      /* set pointer to rank-1 tensor scaling */
      if (rank1_scal_el_mat != NULL) {
        rank1_scal = &rank1_scal_el_data[nodeid];
      }
      else {
        rank1_scal = &dummy;
      }

      /* set pointers to bounds and yielding markers */
      if (bounds_el_mat != NULL) {
        bounds = &bounds_el_data[nodeid];
      }
      else {
        bounds = &dummy;
      }
      if (yielding_el_mat != NULL) {
        yielding = &yielding_el_data[nodeid];
      }
      else {
        yielding = &dummy;
      }

      /* compute nonlinear viscosity in upper mantle */
      if (!prescribe_bounds_yielding) {
        slabs_visc_nl_node (
            &visc_el_data[nodeid], &dvisc_dIIe_el_data[nodeid],
            rank1_scal, bounds, yielding,
            temp, weak, IIe, scaling, visc_temp_decay, physics_options);
      }
      else {
        slabs_visc_nl_node_prescribe_bounds_yielding (
            &visc_el_data[nodeid], &dvisc_dIIe_el_data[nodeid], bounds,
            yielding, temp, weak, IIe, scaling, visc_temp_decay,
            physics_options);
      }
    }
  }
  else { /* if element is located in lower mantle */
    const double        scaling =
                          physics_options->viscosity_lower_mantle_scaling;
    const double        visc_temp_decay =
                          physics_options->viscosity_lower_mantle_temp_decay;

    for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) { /* loop over all
                                                           * nodes */
      const double        temp = temp_el_data[nodeid];
      const double        weak = weak_el_data[nodeid];

      /* check temperature for valid range [0,1] */
      YMIR_ASSERT (0.0 <= temp && temp <= 1.0);

      /* check weak zone for valid range (0,1] */
      YMIR_ASSERT (0.0 < weak && weak <= 1.0);

      /* compute linear viscosity in lower mantle */
      visc_el_data[nodeid] = slabs_visc_temp_node (
          temp, weak, scaling, visc_temp_decay, physics_options, 1);

      /* set derivative to zero */
      dvisc_dIIe_el_data[nodeid] = 0.0;

      /* set rank-1 tensor scaling to zero */
      if (rank1_scal_el_mat != NULL) {
        rank1_scal_el_data[nodeid] = 0.0;
      }

      /* set bounds and yielding marker to zero */
      if (bounds_el_mat != NULL && !prescribe_bounds_yielding) {
        bounds_el_data[nodeid] = 0.0;
      }
      if (yielding_el_mat != NULL && !prescribe_bounds_yielding) {
        yielding_el_data[nodeid] = 0.0;
      }
    }
  }

//TODO del
#if 0
  /* compute viscosity in this element */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) { /* loop over all
                                                         * nodes */
    const double        temp = temp_el_data[nodeid];
    const double        weak = weak_el_data[nodeid];

    /* check temperature for valid range [0,1] */
    YMIR_ASSERT (0.0 <= temp && temp <= 1.0);

    /* check weak zone for valid range (0,1] */
    YMIR_ASSERT (0.0 < weak && weak <= 1.0);

    /*
     * Compute Viscosity
     */

    /* compute viscosity depending on location in lower or upper mantle */
    if (0.0 < upper_mantle_radius
        &&
        slabs_compute_radius (x[nodeid], y[nodeid], z[nodeid], physics_options)
        < upper_mantle_radius) {
      const double scaling = physics_options->viscosity_lower_mantle_scaling;
      const double visc_temp_decay =
        physics_options->viscosity_lower_mantle_temp_decay;

      /* compute linear viscosity in lower mantle */
      visc_el_data[nodeid] = slabs_visc_temp_node (temp, weak, scaling,
                                                   visc_temp_decay,
                                                   physics_options, 1);

      /* set derivative to zero */
      dvisc_dIIe_el_data[nodeid] = 0.0;

      /* set rank-1 tensor scaling to zero */
      if (rank1_scal_el_mat != NULL) {
        rank1_scal_el_data[nodeid] = 0.0;
      }

      /* set bounds and yielding marker to zero */
      if (bounds_el_mat != NULL && !prescribe_bounds_yielding) {
        bounds_el_data[nodeid] = 0.0;
      }
      if (yielding_el_mat != NULL && !prescribe_bounds_yielding) {
        yielding_el_data[nodeid] = 0.0;
      }
    }
    else {
      const double        scaling = physics_options->viscosity_scaling;
      const double        visc_temp_decay =
                            physics_options->viscosity_temp_decay;
      const double        IIe = IIe_el_data[nodeid];
      double             *rank1_scal, *bounds, *yielding, dummy;

      /* check second invariant for non-negativity */
      YMIR_ASSERT (0.0 <= IIe);

      /* set pointers to rank-1 tensor scaling */
      if (rank1_scal_el_mat != NULL) {
        rank1_scal = &rank1_scal_el_data[nodeid];
      }
      else {
        rank1_scal = &dummy;
      }

      /* set pointers to bounds and yielding markers */
      if (bounds_el_mat != NULL) {
        bounds = &bounds_el_data[nodeid];
      }
      else {
        bounds = &dummy;
      }
      if (yielding_el_mat != NULL) {
        yielding = &yielding_el_data[nodeid];
      }
      else {
        yielding = &dummy;
      }

      /* compute nonlinear viscosity in upper mantle */
      if (!prescribe_bounds_yielding) {
        slabs_visc_nl_node (
            &visc_el_data[nodeid], &dvisc_dIIe_el_data[nodeid],
            rank1_scal, bounds, yielding,
            temp, weak, IIe, scaling, visc_temp_decay, physics_options);
      }
      else {
        slabs_visc_nl_node_prescribe_bounds_yielding (
            &visc_el_data[nodeid], &dvisc_dIIe_el_data[nodeid], bounds,
            yielding, temp, weak, IIe, scaling, visc_temp_decay,
            physics_options);
      }
    }
  }
#endif

  /*
   * Check Results
   */

#ifdef YMIR_DEBUG
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) { /* loop over all
                                                         * nodes */
    /* check viscosity for `nan`, `inf`, and positivity */
    YMIR_ASSERT (isfinite (visc_el_data[nodeid]));
    YMIR_ASSERT (0.0 < visc_el_data[nodeid]);

    /* check derivative of viscosity for `nan`, `inf`, and non-positivity */
    YMIR_ASSERT (isfinite (dvisc_dIIe_el_data[nodeid]));
    YMIR_ASSERT (dvisc_dIIe_el_data[nodeid] <= 0.0);

    /* check rank-1 tensor scaling for `nan`, `inf`, and valid range [-1,0] */
    if (rank1_scal_el_mat != NULL) {
      YMIR_ASSERT (isfinite (rank1_scal_el_data[nodeid]));
      YMIR_ASSERT (-1.0 <= rank1_scal_el_data[nodeid]);
      YMIR_ASSERT (rank1_scal_el_data[nodeid] <= 0.0);
    }

    /* check bounds marker for `nan`, `inf`, and valid range [-1,1] */
    if (bounds_el_mat != NULL) {
      YMIR_ASSERT (isfinite (bounds_el_data[nodeid]));
      YMIR_ASSERT (-1.0 <= bounds_el_data[nodeid]);
      YMIR_ASSERT (bounds_el_data[nodeid] <= 1.0);
    }

    /* check yielding marker for `nan`, `inf` and valid range [0,1] */
    if (yielding_el_mat != NULL) {
      YMIR_ASSERT (isfinite (yielding_el_data[nodeid]));
      YMIR_ASSERT (0.0 <= yielding_el_data[nodeid]);
      YMIR_ASSERT (yielding_el_data[nodeid] <= 1.0);
    }
  }
#endif
}

/**
 * Computes the temperature and velocity dependent viscosity vector.
 */
static void
slabs_viscosity_nonlinear (ymir_vec_t *viscosity,
                           ymir_vec_t *dvisc_dIIe,
                           ymir_vec_t *rank1_tensor_scal,
                           ymir_vec_t *bounds_marker,
                           ymir_vec_t *yielding_marker,
                           slabs_stokes_state_t *state,
                           ymir_pressure_elem_t *press_elem,
                           ymir_vel_dir_t *vel_dir,
                           slabs_physics_options_t *physics_options,
                           const int prescribe_bounds_yielding)
{
  ymir_vec_t         *temp_vec = state->temp_vec;
  ymir_vec_t         *weak_vec = state->weak_vec;
  ymir_vec_t         *vel_press_vec = state->vel_press_vec;
  ymir_vec_t         *vel_bc_vec = state->vel_bc_vec;
  ymir_mesh_t        *mesh = viscosity->mesh;
  mangll_t           *mangll = mesh->ma;
  const mangll_locidx_t  n_elements = mesh->cnodes->K;
  const int           N = ymir_n (mangll->N);
  const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);

  double             *x, *y, *z, *tmp_el;
  sc_dmatrix_t       *temp_el_mat, *weak_el_mat;
  sc_dmatrix_t       *vel_el_mat, *IIe_el_mat;
  sc_dmatrix_t       *tmp_grad_vel, *tmp_dvel, *tmp_vel;
  sc_dmatrix_t       *visc_el_mat, *dvisc_dIIe_el_mat, *rank1_scal_el_mat;
  sc_dmatrix_t       *bounds_el_mat, *yielding_el_mat;
  mangll_locidx_t     elid;

  /* check input */
  YMIR_ASSERT (dvisc_dIIe != NULL);
  YMIR_ASSERT (press_elem != NULL);
  YMIR_ASSERT (vel_dir != NULL);
  YMIR_ASSERT (state->temperature != NULL);
  YMIR_ASSERT (state->temp_vec != NULL);
  YMIR_ASSERT (state->weakzone != NULL);
  YMIR_ASSERT (state->weak_vec != NULL);
  YMIR_ASSERT (state->velocity != NULL);
  YMIR_ASSERT (state->vel_press_vec != NULL);
#ifdef YMIR_DEBUG
  switch (physics_options->viscosity_model_type) {
  case SL_VISCOSITY_MODEL_WYUL:
  case SL_VISCOSITY_MODEL_UWYUL:
  case SL_VISCOSITY_MODEL_UWYL:
  case SL_VISCOSITY_MODEL_UWYL_LREG:
  case SL_VISCOSITY_MODEL_UWYL_SHIFT_LREG:
  case SL_VISCOSITY_MODEL_UYWL:
  case SL_VISCOSITY_MODEL_UYWL_SHIFT:
    YMIR_ASSERT (0.0 < physics_options->viscosity_max);
    break;
  case SL_VISCOSITY_MODEL_UWL_IIE_REG:
    YMIR_ASSERT (0.0 < physics_options->viscosity_IIe_regularization);
    break;
  default: /* unknown viscosity model */
    YMIR_ABORT_NOT_REACHED ();
  }
#endif

  /* check if it is necessary to compute the strain rate dependent viscosity */
  if (fabs (physics_options->viscosity_stress_exponent - 1.0) < SC_EPS) {
    /* compute just temperature dependent viscosity */
    slabs_viscosity_linear (viscosity, state, physics_options);

    /* set derivative to zero */
    ymir_dvec_set_zero (dvisc_dIIe);

    /* set rank-1 tensor scaling to zero */
    if (rank1_tensor_scal != NULL) {
      ymir_dvec_set_zero (rank1_tensor_scal);
    }

    /* set bounds and yielding marker to zero */
    if (bounds_marker != NULL) {
      ymir_dvec_set_zero (bounds_marker);
    }
    if (yielding_marker != NULL) {
      ymir_dvec_set_zero (yielding_marker);
    }

    /* end execution */
    return;
  }

  /* get velocity */
  ymir_stokes_vec_get_velocity (vel_press_vec, vel_bc_vec, press_elem);
  ymir_vel_dir_separate (vel_bc_vec, NULL, NULL, NULL, vel_dir);

  /* create work variables */
  x = YMIR_ALLOC (double, n_nodes_per_el);
  y = YMIR_ALLOC (double, n_nodes_per_el);
  z = YMIR_ALLOC (double, n_nodes_per_el);
  tmp_el = YMIR_ALLOC (double, n_nodes_per_el);
  temp_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  weak_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
  IIe_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  tmp_grad_vel = sc_dmatrix_new (n_nodes_per_el, 9);
  tmp_dvel = sc_dmatrix_new (n_nodes_per_el, 3);
  tmp_vel = sc_dmatrix_new (n_nodes_per_el, 3);
  visc_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  dvisc_dIIe_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  if (rank1_tensor_scal != NULL) {
    rank1_scal_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  }
  else {
    rank1_scal_el_mat = NULL;
  }
  if (bounds_marker != NULL) {
    bounds_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  }
  else {
    bounds_el_mat = NULL;
  }
  if (yielding_marker != NULL) {
    yielding_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  }
  else {
    yielding_el_mat = NULL;
  }

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    /* get coordinates of this element at Gauss nodes */
    slabs_elem_get_gauss_coordinates (x, y, z, elid, mangll, tmp_el);

    /* get temperature field of this element from state at Gauss nodes */
    ymir_cvec_get_elem_interp (temp_vec, temp_el_mat, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_READ);
    slabs_matrix_bound_values (temp_el_mat, 0.0, 1.0);

    /* get weak zone of this element */
    ymir_dvec_get_elem (weak_vec, weak_el_mat, YMIR_STRIDE_NODE, elid,
                        YMIR_READ);

    /* get velocity field of this element from state at GLL nodes */
    ymir_cvec_get_elem_interp (vel_bc_vec, vel_el_mat, YMIR_STRIDE_NODE,
                               elid, YMIR_GLL_NODE, YMIR_READ);

    /* compute 2nd invariant of the strain rate at Gauss nodes */
    slabs_second_invariant_elem (vel_el_mat, IIe_el_mat, mangll, elid,
                                 tmp_grad_vel, tmp_dvel, tmp_vel,
                                 SL_GAUSS_NODE);

    /* get bounds marker for this element */
    if (bounds_marker != NULL) {
      ymir_dvec_get_elem (bounds_marker, bounds_el_mat, YMIR_STRIDE_NODE,
                          elid, YMIR_RW);
    }

    /* get yielding marker for this element */
    if (yielding_marker != NULL) {
      ymir_dvec_get_elem (yielding_marker, yielding_el_mat, YMIR_STRIDE_NODE,
                          elid, YMIR_RW);
    }

    /* compute viscosity */
    slabs_visc_nl_elem (visc_el_mat, dvisc_dIIe_el_mat, rank1_scal_el_mat,
                        bounds_el_mat, yielding_el_mat,
                        x, y, z, mangll->refel->Vmask,
                        temp_el_mat, weak_el_mat, IIe_el_mat,
                        physics_options, prescribe_bounds_yielding);

    /* set viscosity for this element */
    ymir_dvec_set_elem (viscosity, visc_el_mat,
                        YMIR_STRIDE_NODE, elid, YMIR_SET);

    /* set derivative for viscosity of this element */
    ymir_dvec_set_elem (dvisc_dIIe, dvisc_dIIe_el_mat,
                        YMIR_STRIDE_NODE, elid, YMIR_SET);

    /* set rank-1 tensor scaling for this element */
    if (rank1_tensor_scal != NULL) {
      ymir_dvec_set_elem (rank1_tensor_scal, rank1_scal_el_mat,
                          YMIR_STRIDE_NODE, elid, YMIR_SET);
    }

    /* set bounds marker for this element */
    if (bounds_marker != NULL) {
      ymir_dvec_set_elem (bounds_marker, bounds_el_mat,
                          YMIR_STRIDE_NODE, elid, YMIR_SET);
    }

    /* set yielding marker for this element */
    if (yielding_marker != NULL) {
      ymir_dvec_set_elem (yielding_marker, yielding_el_mat,
                          YMIR_STRIDE_NODE, elid, YMIR_SET);
    }
  }

  /* destroy */
  sc_dmatrix_destroy (temp_el_mat);
  sc_dmatrix_destroy (weak_el_mat);
  sc_dmatrix_destroy (vel_el_mat);
  sc_dmatrix_destroy (IIe_el_mat);
  sc_dmatrix_destroy (tmp_grad_vel);
  sc_dmatrix_destroy (tmp_dvel);
  sc_dmatrix_destroy (tmp_vel);
  sc_dmatrix_destroy (visc_el_mat);
  sc_dmatrix_destroy (dvisc_dIIe_el_mat);
  if (rank1_tensor_scal != NULL) {
    sc_dmatrix_destroy (rank1_scal_el_mat);
  }
  if (bounds_marker != NULL) {
    sc_dmatrix_destroy (bounds_el_mat);
  }
  if (yielding_marker != NULL) {
    sc_dmatrix_destroy (yielding_el_mat);
  }
  YMIR_FREE (x);
  YMIR_FREE (y);
  YMIR_FREE (z);
  YMIR_FREE (tmp_el);
}

/**
 * Computes viscosity.
 */
void
slabs_physics_compute_viscosity (ymir_vec_t *viscosity,
                                 ymir_vec_t *dvisc_dIIe,
                                 ymir_vec_t *rank1_tensor_scal,
                                 ymir_vec_t *bounds_marker,
                                 ymir_vec_t *yielding_marker,
                                 slabs_stokes_state_t *state,
                                 ymir_pressure_elem_t *press_elem,
                                 ymir_vel_dir_t *vel_dir,
                                 slabs_physics_options_t *physics_options,
                                 const int prescribe_bounds_yielding)
{
  /* check input */
  YMIR_ASSERT (viscosity != NULL);
  YMIR_ASSERT_IS_DVEC (viscosity);
  YMIR_ASSERT (viscosity->ndfields == 1);
  YMIR_ASSERT (viscosity->node_type == YMIR_GAUSS_NODE);
  if (dvisc_dIIe != NULL) {
    YMIR_ASSERT_IS_DVEC (dvisc_dIIe);
    YMIR_ASSERT (dvisc_dIIe->ndfields == 1);
    YMIR_ASSERT (dvisc_dIIe->node_type == YMIR_GAUSS_NODE);
  }
  if (rank1_tensor_scal != NULL) {
    YMIR_ASSERT_IS_DVEC (rank1_tensor_scal);
    YMIR_ASSERT (rank1_tensor_scal->ndfields == 1);
    YMIR_ASSERT (rank1_tensor_scal->node_type == YMIR_GAUSS_NODE);
  }
  if (bounds_marker != NULL) {
    YMIR_ASSERT_IS_DVEC (bounds_marker);
    YMIR_ASSERT (bounds_marker->ndfields == 1);
    YMIR_ASSERT (bounds_marker->node_type == YMIR_GAUSS_NODE);
  }
  if (yielding_marker != NULL) {
    YMIR_ASSERT_IS_DVEC (yielding_marker);
    YMIR_ASSERT (yielding_marker->ndfields == 1);
    YMIR_ASSERT (yielding_marker->node_type == YMIR_GAUSS_NODE);
  }
  YMIR_ASSERT (!prescribe_bounds_yielding || bounds_marker != NULL);
  YMIR_ASSERT (!prescribe_bounds_yielding || yielding_marker != NULL);

  /* set viscosity and its derivative */
  switch (physics_options->viscosity_type) {
  case SL_VISCOSITY_CONST: /* constant viscosity */
    ymir_vec_set_value (viscosity, physics_options->viscosity_scaling);
    break;

  case SL_VISCOSITY_LINEAR: /* linear, temperature dependent viscosity */
    slabs_viscosity_linear (viscosity, state, physics_options);
    break;

  case SL_VISCOSITY_NONLINEAR: /* nonlinear viscosity */
    slabs_viscosity_nonlinear (viscosity, dvisc_dIIe, rank1_tensor_scal,
                               bounds_marker, yielding_marker,
                               state, press_elem, vel_dir, physics_options,
                               prescribe_bounds_yielding);
    break;

  default: /* unknown viscosity type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* check viscosity for `nan` and `inf` */
  YMIR_ASSERT (sc_dmatrix_is_valid (viscosity->dvec));
  if (physics_options->viscosity_type == SL_VISCOSITY_NONLINEAR) {
    if (dvisc_dIIe != NULL) {
      YMIR_ASSERT (sc_dmatrix_is_valid (dvisc_dIIe->dvec));
    }
    if (rank1_tensor_scal != NULL) {
      YMIR_ASSERT (sc_dmatrix_is_valid (rank1_tensor_scal->dvec));
    }
  }
}

/**
 * Computes viscosity for linear Stokes operator.
 */
void
slabs_physics_compute_linear_viscosity (ymir_dvec_t *viscosity,
                                        slabs_stokes_state_t *state,
                                        slabs_physics_options_t
                                          *physics_options)
{
  slabs_physics_compute_viscosity (viscosity, NULL, NULL, NULL, NULL, state,
                                   NULL, NULL, physics_options, 0);
}

/**
 * Sets constant viscosity in upper mantle.
 */
static void
slabs_viscosity_set_const_upper_mantle_fn (double *viscosity,
                                           double x, double y, double z,
                                           ymir_locidx_t nid, void *data)
{
  slabs_physics_options_t *physics_options = (slabs_physics_options_t *) data;
  const double        upper_mantle_radius =
                        physics_options->viscosity_upper_mantle_radius;
  const double        visc_min = SC_MAX (0.0, physics_options->viscosity_min);
  const double        visc_max = SC_MAX (0.0, physics_options->viscosity_max);
  double              upper_mantle_visc;

  /* pick viscosity for upper mantle */
  switch (physics_options->viscosity_type_for_init_nl_stokes) {
  case SL_VISCOSITY_INIT_NL_STOKES_CONST_CENTERED_TO_BOUNDS:
    upper_mantle_visc = sqrt (visc_min * visc_max);
    break;

  case SL_VISCOSITY_INIT_NL_STOKES_CONST_MIN_BOUND:
    upper_mantle_visc = visc_min;
    break;

  default: /* unknown initial viscosity type for nonlinear Stokes */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* set constant viscosity in upper mantle */
  if (0.0 < upper_mantle_radius &&
      upper_mantle_radius <= slabs_compute_radius (x, y, z, physics_options)) {
    *viscosity = upper_mantle_visc;
  }
}

/**
 *
 */
void
slabs_physics_compute_init_nl_viscosity (ymir_vec_t *viscosity,
                                         ymir_vec_t *dvisc_dIIe,
                                         ymir_vec_t *rank1_tensor_scal,
                                         ymir_vec_t *bounds_marker,
                                         ymir_vec_t *yielding_marker,
                                         slabs_stokes_state_t *state,
                                         ymir_pressure_elem_t *press_elem,
                                         ymir_vel_dir_t *vel_dir,
                                         slabs_physics_options_t
                                           *physics_options)
{
  ymir_vec_t         *vel_press_vec = state->vel_press_vec;
  ymir_cvec_t        *vel_bc_vec = state->vel_bc_vec;
  const slabs_viscosity_t  visc_type = physics_options->viscosity_type;
  const slabs_viscosity_init_nl_stokes_t  visc_type_init =
                        physics_options->viscosity_type_for_init_nl_stokes;
  const double        visc_scaling = physics_options->viscosity_scaling;
  const double        upper_mantle_radius =
                        physics_options->viscosity_upper_mantle_radius;
  const double        visc_lower_mantle_scaling =
                        physics_options->viscosity_lower_mantle_scaling;
#ifdef YMIR_DEBUG
  const double        visc_min = physics_options->viscosity_min;
  const double        visc_max = physics_options->viscosity_max;
#endif

  /* check input */
  YMIR_ASSERT (viscosity != NULL);
  YMIR_ASSERT_IS_DVEC (viscosity);
  YMIR_ASSERT (viscosity->ndfields == 1);
  YMIR_ASSERT (viscosity->node_type == YMIR_GAUSS_NODE);
  YMIR_ASSERT (dvisc_dIIe != NULL);
  YMIR_ASSERT_IS_DVEC (dvisc_dIIe);
  YMIR_ASSERT (dvisc_dIIe->ndfields == 1);
  YMIR_ASSERT (dvisc_dIIe->node_type == YMIR_GAUSS_NODE);
  if (rank1_tensor_scal != NULL) {
    YMIR_ASSERT_IS_DVEC (rank1_tensor_scal);
    YMIR_ASSERT (rank1_tensor_scal->ndfields == 1);
    YMIR_ASSERT (rank1_tensor_scal->node_type == YMIR_GAUSS_NODE);
  }
  if (bounds_marker != NULL) {
    YMIR_ASSERT_IS_DVEC (bounds_marker);
    YMIR_ASSERT (bounds_marker->ndfields == 1);
    YMIR_ASSERT (bounds_marker->node_type == YMIR_GAUSS_NODE);
  }
  if (yielding_marker != NULL) {
    YMIR_ASSERT_IS_DVEC (yielding_marker);
    YMIR_ASSERT (yielding_marker->ndfields == 1);
    YMIR_ASSERT (yielding_marker->node_type == YMIR_GAUSS_NODE);
  }

  switch (visc_type_init) {
  case SL_VISCOSITY_INIT_NL_STOKES_DEFAULT:
    YMIR_ASSERT (press_elem != NULL);
    YMIR_ASSERT (vel_dir != NULL);

    /* set viscosity and it's derivative */
    slabs_physics_compute_viscosity (viscosity, dvisc_dIIe, rank1_tensor_scal,
                                     bounds_marker, yielding_marker,
                                     state, press_elem, vel_dir,
                                     physics_options, 0);
    break;

  case SL_VISCOSITY_INIT_NL_STOKES_CONST_CENTERED_TO_BOUNDS:
  case SL_VISCOSITY_INIT_NL_STOKES_CONST_MIN_BOUND:
    /* check input parameters */
    YMIR_ASSERT (0.0 < visc_min);
    YMIR_ASSERT (0.0 < visc_max);

    /* get velocity (for nonlinear Stokes operator) */
    ymir_stokes_vec_get_velocity (vel_press_vec, vel_bc_vec, press_elem);
    ymir_vel_dir_separate (vel_bc_vec, NULL, NULL, NULL, vel_dir);

    /* set viscosity (everywhere) */
    slabs_physics_compute_linear_viscosity (viscosity, state, physics_options);

    /* modify viscosity in upper mantle */
    ymir_dvec_set_function (
        viscosity, slabs_viscosity_set_const_upper_mantle_fn, physics_options);

    /* set derivative to zero */
    ymir_dvec_set_zero (dvisc_dIIe);

    /* set rank-1 tensor scaling to zero */
    if (rank1_tensor_scal != NULL) {
      ymir_dvec_set_zero (rank1_tensor_scal);
    }

    /* set bounds and yielding marker to zero */
    if (bounds_marker != NULL) {
      ymir_dvec_set_zero (bounds_marker);
    }
    if (yielding_marker != NULL) {
      ymir_dvec_set_zero (yielding_marker);
    }
    break;

  case SL_VISCOSITY_INIT_NL_STOKES_TEMP:
  case SL_VISCOSITY_INIT_NL_STOKES_TEMP_UM_REL_TO_LM:
    /* get velocity (for nonlinear Stokes operator) */
    ymir_stokes_vec_get_velocity (vel_press_vec, vel_bc_vec, press_elem);
    ymir_vel_dir_separate (vel_bc_vec, NULL, NULL, NULL, vel_dir);

    /* change scaling of viscosity in upper mantle */
    if (   visc_type_init == SL_VISCOSITY_INIT_NL_STOKES_TEMP_UM_REL_TO_LM
        && 0.0 < upper_mantle_radius) {
      physics_options->viscosity_scaling = SL_VISCOSITY_UM_REL_TO_LM
                                           * visc_lower_mantle_scaling;
    }

    /* compute temperature dependent viscosity */
    physics_options->viscosity_type = SL_VISCOSITY_LINEAR;
    slabs_physics_compute_linear_viscosity (viscosity, state, physics_options);
    physics_options->viscosity_type = visc_type;

    /* restore viscosity scaling */
    if (   visc_type_init == SL_VISCOSITY_INIT_NL_STOKES_TEMP_UM_REL_TO_LM
        && 0.0 < upper_mantle_radius) {
      physics_options->viscosity_scaling = visc_scaling;
    }

    /* set derivative of viscosity */
    ymir_dvec_set_zero (dvisc_dIIe);

    /* set rank-1 tensor scaling to zero */
    if (rank1_tensor_scal != NULL) {
      ymir_dvec_set_zero (rank1_tensor_scal);
    }

    /* set bounds and yielding marker to zero */
    if (bounds_marker != NULL) {
      ymir_dvec_set_zero (bounds_marker);
    }
    if (yielding_marker != NULL) {
      ymir_dvec_set_zero (yielding_marker);
    }
    break;

  default: /* unknown initial viscosity type for nonlinear Stokes */
    YMIR_ABORT_NOT_REACHED ();
  }
}

/**
 * Computes coefficient for linear or nonlinear Stokes operator.
 */
static void
slabs_physics_compute_stokes_coeff_int (ymir_vec_t *coeff,
                                        ymir_vec_t *coeff_deriv,
                                        ymir_vec_t *rank1_tensor_scal,
                                        ymir_vec_t *bounds_marker,
                                        ymir_vec_t *yielding_marker,
                                        slabs_stokes_state_t *state,
                                        ymir_pressure_elem_t *press_elem,
                                        ymir_vel_dir_t *vel_dir,
                                        slabs_physics_options_t
                                          *physics_options,
                                        const int prescribe_bounds_yielding)
{
  /* compute viscosity */
  slabs_physics_compute_viscosity (coeff, coeff_deriv, rank1_tensor_scal,
                                   bounds_marker, yielding_marker,
                                   state, press_elem, vel_dir,
                                   physics_options, prescribe_bounds_yielding);

  /* scale viscosity to get the coefficient */
  ymir_vec_scale (2.0, coeff);
  ymir_vec_scale (2.0, coeff_deriv);

  /* store marker for coarsening function; assume: this function is called
   * right before building the GMG hierarchy */
  if (physics_options->viscosity_coarsen_eval && !prescribe_bounds_yielding &&
      bounds_marker != NULL) {
    physics_options->current_fine_bounds = bounds_marker;

  }
  if (physics_options->viscosity_coarsen_eval && !prescribe_bounds_yielding &&
      yielding_marker != NULL) {
    physics_options->current_fine_yielding = yielding_marker;
  }
}

void
slabs_physics_compute_stokes_coeff (ymir_vec_t *coeff,
                                    ymir_vec_t *coeff_deriv,
                                    ymir_vec_t *rank1_tensor_scal,
                                    ymir_vec_t *bounds_marker,
                                    ymir_vec_t *yielding_marker,
                                    slabs_stokes_state_t *state,
                                    ymir_pressure_elem_t *press_elem,
                                    ymir_vel_dir_t *vel_dir,
                                    slabs_physics_options_t *physics_options)
{
  slabs_physics_compute_stokes_coeff_int (
      coeff, coeff_deriv, rank1_tensor_scal, bounds_marker, yielding_marker,
      state, press_elem, vel_dir, physics_options, 0);
}

/**
 * Computes coefficient for linear Stokes operator.
 */
void
slabs_physics_compute_lin_stokes_coeff (ymir_vec_t *coeff,
                                        slabs_stokes_state_t *state,
                                        slabs_physics_options_t
                                          *physics_options)
{
  /* compute viscosity */
  slabs_physics_compute_viscosity (coeff, NULL, NULL, NULL, NULL,
                                   state, NULL, NULL, physics_options, 0);

  /* scale viscosity to get the coefficient */
  ymir_vec_scale (2.0, coeff);

  /* remove marker for coarsening function */
  if (physics_options->viscosity_coarsen_eval) {
    physics_options->current_fine_bounds = NULL;
    physics_options->current_fine_yielding = NULL;
  }
}

/**
 * Computes coefficient for nonlinear Stokes operator at initial nonlinear step.
 */
void
slabs_physics_compute_init_nl_stokes_coeff (ymir_vec_t *coeff,
                                            ymir_vec_t *coeff_deriv,
                                            ymir_vec_t *rank1_tensor_scal,
                                            ymir_vec_t *bounds_marker,
                                            ymir_vec_t *yielding_marker,
                                            slabs_stokes_state_t *state,
                                            ymir_pressure_elem_t *press_elem,
                                            ymir_vel_dir_t *vel_dir,
                                            slabs_physics_options_t
                                              *physics_options)
{
  /* compute viscosity */
  slabs_physics_compute_init_nl_viscosity (coeff, coeff_deriv,
                                           rank1_tensor_scal, bounds_marker,
                                           yielding_marker, state, press_elem,
                                           vel_dir, physics_options);

  /* scale viscosity to get the coefficient */
  ymir_vec_scale (2.0, coeff);
  ymir_vec_scale (2.0, coeff_deriv);

  /* remove marker for coarsening function */
  if (physics_options->viscosity_coarsen_eval) {
    physics_options->current_fine_bounds = NULL;
    physics_options->current_fine_yielding = NULL;
  }
}

/**
 *
 */
slabs_physics_coarsen_stokes_coeff_data_t *
slabs_physics_coarsen_stokes_coeff_data_new (slabs_stokes_state_t *state,
                                             slabs_physics_options_t
                                               *physics_options)
{
  slabs_physics_coarsen_stokes_coeff_data_t *d;

  d = YMIR_ALLOC (slabs_physics_coarsen_stokes_coeff_data_t, 1);
  d->state = state;
  d->physics_options = physics_options;
  d->init_nl_stokes_coeff = -1;
  d->coarsen_count = 0;
  d->buffer_fine_weak = NULL;
  d->buffer_fine_temp = NULL;
  d->buffer_fine_vel = NULL;
  d->buffer_fine_bounds = NULL;
  d->buffer_fine_yielding = NULL;

  return d;
}

/**
 *
 */
void
slabs_physics_coarsen_stokes_coeff_data_reset (void *data)
{
  slabs_physics_coarsen_stokes_coeff_data_t *d =
    (slabs_physics_coarsen_stokes_coeff_data_t *) data;

  /* reset counter */
  d->coarsen_count = 0;

  /* destroy buffers */
  if (d->buffer_fine_weak != NULL) {
    ymir_vec_destroy (d->buffer_fine_weak);
    d->buffer_fine_weak = NULL;
  }
  if (d->buffer_fine_temp != NULL) {
    ymir_vec_destroy (d->buffer_fine_temp);
    d->buffer_fine_temp = NULL;
  }
  if (d->buffer_fine_vel != NULL) {
    ymir_vec_destroy (d->buffer_fine_vel);
    d->buffer_fine_vel = NULL;
  }
  if (d->buffer_fine_bounds != NULL) {
    ymir_vec_destroy (d->buffer_fine_bounds);
    d->buffer_fine_bounds = NULL;
  }
  if (d->buffer_fine_yielding != NULL) {
    ymir_vec_destroy (d->buffer_fine_yielding);
    d->buffer_fine_yielding = NULL;
  }
}

/**
 *
 */
void
slabs_physics_coarsen_stokes_coeff_data_destroy (
                               slabs_physics_coarsen_stokes_coeff_data_t *data)
{
  if (data == NULL) {
    return;
  }

  slabs_physics_coarsen_stokes_coeff_data_reset (data);
  YMIR_FREE (data);
}

/**
 *
 */
static void
slabs_physics_coarsen_bounds_preprocess_fine (ymir_vec_t *bounds_marker)
{
  sc_dmatrix_t       *bounds_mat = bounds_marker->dvec;
  const int           totalsize = bounds_mat->m * bounds_mat->n;
  double             *bounds_data = bounds_mat->e[0];
  int                 i;

  /* check input */
  YMIR_ASSERT_HAS_DVEC (bounds_marker);

  for (i = 0; i < totalsize; i++) {
    if (0.0 < bounds_data[i]) { /* if upper bound is active */
      bounds_data[i] = 1.0;
    }
    else if (bounds_data[i] < 0.0) { /* if lower bound is active */
      bounds_data[i] = -1.0;
    }
  }
}

/**
 *
 */
static void
slabs_physics_coarsen_bounds_yielding_preprocess (ymir_vec_t *bounds_marker,
                                                  ymir_vec_t *yielding_marker)
{
  const ymir_locidx_t n_elements = bounds_marker->mesh->cnodes->K;
  const int           n_nodes_per_el = bounds_marker->mesh->cnodes->Np;
  const int           n_fields = bounds_marker->ndfields;

  sc_dmatrix_t       *bounds_mat = sc_dmatrix_new (n_fields, n_nodes_per_el);
  sc_dmatrix_t       *yielding_mat = sc_dmatrix_new (n_fields, n_nodes_per_el);
  double             *bounds_data, *yielding_data;
  int                 upper_bound_active, lower_bound_active/*, yielding_active*/;
  ymir_locidx_t       elid;
  int                 nodeid;

  /* check input */
  YMIR_ASSERT_HAS_DVEC (bounds_marker);
  YMIR_ASSERT_HAS_DVEC (yielding_marker);
  YMIR_ASSERT (bounds_marker->ndfields == 1);
  YMIR_ASSERT (yielding_marker->ndfields == 1);
  YMIR_ASSERT (yielding_marker->mesh->cnodes->K == n_elements);
  YMIR_ASSERT (yielding_marker->mesh->cnodes->Np == n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) {
    /* get values at this element */
    ymir_dvec_get_elem (bounds_marker, bounds_mat, YMIR_STRIDE_COMP, elid,
                        YMIR_RW);
    ymir_dvec_get_elem (yielding_marker, yielding_mat, YMIR_STRIDE_COMP, elid,
                        YMIR_RW);
    bounds_data = bounds_mat->e[0];
    yielding_data = yielding_mat->e[0];

    /* check if bounds or yielding is active on any of the nodes */
    upper_bound_active = 0;
    lower_bound_active = 0;
    //yielding_active = 0;
    for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
      if (0.0 < bounds_data[nodeid]) { /* if upper bound is active */
        upper_bound_active = 1;
      }
      else if (bounds_data[nodeid] < 0.0) { /* if lower bound is active */
        lower_bound_active = 1;
      }

      if (0.0 < yielding_data[nodeid]) { /* if yielding is active */
        //yielding_active = 1;
      }
    }

    /* activate upper/lower bound uniformly at all nodes */
    if (upper_bound_active) {
      for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
        bounds_data[nodeid] = 1.0;
      }
    }
    else if (lower_bound_active) {
      for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
        bounds_data[nodeid] = -1.0;
      }
    }

    /* activate yielding uniformly at all nodes */
    /*
    if (yielding_active) {
      for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
        yielding_data[nodeid] = 1.0;
      }
    }
    */

    /* set values at this element */
    ymir_dvec_set_elem (bounds_marker, bounds_mat, YMIR_STRIDE_COMP, elid,
                        YMIR_SET);
    ymir_dvec_set_elem (yielding_marker, yielding_mat, YMIR_STRIDE_COMP, elid,
                        YMIR_SET);
  }

  /* destroy */
  sc_dmatrix_destroy (bounds_mat);
  sc_dmatrix_destroy (yielding_mat);
}

/**
 * Corrects bounds and yielding markers after they were coarsened.
 * It is important that no over-/undershooting was introduced by coarsening,
 * which should be the case for linear interpolation.
 */
static void
slabs_physics_coarsen_bounds_yielding_postprocess (ymir_vec_t *bounds_marker,
                                                   ymir_vec_t *yielding_marker)
{
  sc_dmatrix_t       *bounds_mat = bounds_marker->dvec;
  sc_dmatrix_t       *yielding_mat = yielding_marker->dvec;
  const int           totalsize = bounds_mat->m * bounds_mat->n;
  double             *bounds_data = bounds_mat->e[0];
  double             *yielding_data = yielding_mat->e[0];
  int                 i;

  /* check input */
  YMIR_ASSERT_HAS_DVEC (bounds_marker);
  YMIR_ASSERT_HAS_DVEC (yielding_marker);
  YMIR_ASSERT ((yielding_mat->m * yielding_mat->n) == totalsize);

  /* enforce valid values for markers */
  for (i = 0; i < totalsize; i++) {
    bounds_data[i] = SC_MIN (1.0, bounds_data[i]);
    bounds_data[i] = SC_MAX (-1.0, bounds_data[i]);

    yielding_data[i] = SC_MIN (1.0, yielding_data[i]);
    yielding_data[i] = SC_MAX (0.0, yielding_data[i]);
  }
}

/**
 * Computes coefficient for linear or nonlinear Stokes operator.
 */
void
slabs_physics_coarsen_stokes_coeff (ymir_vec_t *coarse_coeff,
                                    ymir_vec_t *coarse_coeff_deriv,
                                    ymir_vec_t *coarse_vel,
                                    mangll_t *upstream_interp_mangll,
                                    mangll_t *upstream_partition_mangll,
                                    ymir_mesh_t *coarse_ymir_mesh,
                                    ymir_pressure_elem_t *coarse_press_elem,
                                    ymir_vel_dir_t *coarse_vel_dir,
                                    p8est_t *coarse_p8est,
                                    const int is_fine_level,
                                    void *data)
{
  slabs_physics_coarsen_stokes_coeff_data_t *d =
    (slabs_physics_coarsen_stokes_coeff_data_t *) data;
  slabs_stokes_state_t *fine_state = d->state;
  slabs_physics_options_t  *physics_options = d->physics_options;
  slabs_viscosity_t         visc_type = physics_options->viscosity_type;
  slabs_weakzone_coarsen_t  weakzone_coarsen_type =
                              physics_options->weakzone_coarsen_type;
  slabs_viscosity_coarsen_t viscosity_coarsen_type =
                              physics_options->viscosity_coarsen_type;
  const double        transition_zone_incr =
    physics_options->viscosity_lower_upper_transition_zone_incr;
  const int           p_projection =
    (upstream_interp_mangll == NULL && upstream_partition_mangll == NULL);

  slabs_stokes_state_t *coarse_state;
  ymir_vec_t         *fine_weak, *fine_temp, *fine_vel;
#ifdef SL_COARSEN_STOKES_COEFF_VTK
  char                filename[BUFSIZ];

  /* set filename for vtk output */
  snprintf (filename, BUFSIZ, "slabs_physics_coarsen_stokes_coeff_%02i",
            d->coarsen_count + 1);
#endif

  /* set fine fields */
  if (is_fine_level) { /* if at finest level in hierarchy */
    fine_weak = fine_state->weak_vec;
    fine_temp = fine_state->temp_vec;
    fine_vel = fine_state->vel_bc_vec;

    /* set whether this is the first call (after new) or not */
    if (d->init_nl_stokes_coeff == -1) {
      d->init_nl_stokes_coeff = 1;
    }
    else {
      d->init_nl_stokes_coeff = 0;
    }
  }
  else { /* if at intermediate level in hierarchy */
    fine_weak = d->buffer_fine_weak;
    fine_temp = d->buffer_fine_temp;
    fine_vel = d->buffer_fine_vel;
  }

  /* create Stokes state on coarse mesh */
  coarse_state = slabs_stokes_state_new (coarse_p8est);

  /*
   * Compute Weak Zone for Coarse Stokes State
   */

  slabs_stokes_state_init_weakzone (coarse_state, coarse_ymir_mesh,
                                    slabs_physics_compute_weakzone,
                                    physics_options);

  switch (weakzone_coarsen_type) {
  case SL_WEAKZONE_COARSEN_EVAL:
    /* evaluate exact weak zone field */
    slabs_physics_compute_weakzone (coarse_state->weak_vec, physics_options);
    break;

  case SL_WEAKZONE_COARSEN_EVAL_PIECEWISE_CONST:
    /* evaluate exact weak zone field and set const inside of elements */
    slabs_physics_compute_weakzone (coarse_state->weak_vec, physics_options);
    ymir_dvec_set_const_inside_elems (coarse_state->weak_vec);
    break;

  case SL_WEAKZONE_COARSEN_LIN_INTERP:
    {
      const double        visc_min = physics_options->viscosity_min;
      const double        visc_max = physics_options->viscosity_max;

      /* restrict weak zone field from fine to coarse mesh (lin interp) */
      YMIR_ASSERT (fine_weak != NULL);
      ymir_gmg_intergrid_restrict_coeff (
          fine_weak, upstream_interp_mangll,
          upstream_partition_mangll, coarse_state->weak_vec,
          1, 2.0 * visc_min, 2.0 * visc_max, 1, 0);

      /* store coarse weak zone field in buffer (for next coarsening) */
      if (d->buffer_fine_weak != NULL) {
        ymir_vec_destroy (d->buffer_fine_weak);
      }
      d->buffer_fine_weak = ymir_vec_clone (coarse_state->weak_vec);
    }
    break;

  default: /* unknown weak zone coarsening type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /*
   * Restrict Temperature Field from Fine to Coarse Mesh
   */

  slabs_stokes_state_init_temp (coarse_state, coarse_ymir_mesh->cnodes);
  slabs_stokes_state_init_temp_vec (coarse_state, coarse_ymir_mesh);
  YMIR_ASSERT (fine_temp != NULL);
  if (p_projection) {
    ymir_gmg_intergrid_p_restrict_cnode_accurate (
        fine_temp, coarse_state->temp_vec);
  }
  else {
    ymir_gmg_intergrid_restrict_partition_cnode_accurate (
        fine_temp, upstream_interp_mangll,
        upstream_partition_mangll, coarse_state->temp_vec);
  }

  /* bound temperature to valid interval */
  slabs_matrix_bound_values (coarse_state->temperature, 0.0, 1.0);

  /* store coarse temperature field in buffer (for next coarsening) */
  if (d->buffer_fine_temp != NULL) {
    ymir_vec_destroy (d->buffer_fine_temp);
  }
  d->buffer_fine_temp = ymir_vec_clone (coarse_state->temp_vec);

  /*
   * Compute Coefficient of Stokes Operator on Coarse Mesh
   */

  if (visc_type != SL_VISCOSITY_NONLINEAR) { /* if coarsen linear coefficient */
    physics_options->viscosity_lower_upper_transition_zone =
      ((double) (d->coarsen_count + 1)) * transition_zone_incr;
    slabs_physics_compute_lin_stokes_coeff (coarse_coeff, coarse_state,
                                            physics_options);
    physics_options->viscosity_lower_upper_transition_zone = 0.0;

#ifdef SL_COARSEN_STOKES_COEFF_VTK
    /* vtk output */
    ymir_vtk_write (coarse_ymir_mesh, filename, coarse_coeff, "coeff", NULL);
#endif
  }
  else { /* if coarsen nonlinear coefficient */
    ymir_vec_t         *coarse_coeff_deriv_int, *coarse_vel_int;

    YMIR_ASSERT (fine_vel != NULL);
    YMIR_ASSERT (coarse_press_elem != NULL);
    YMIR_ASSERT (coarse_vel_dir != NULL);

    /* set internal vectors for the case of Picard PC for nonlinear viscosity */
    if (coarse_coeff_deriv != NULL) {
      coarse_coeff_deriv_int = coarse_coeff_deriv;
    }
    else {
      coarse_coeff_deriv_int = ymir_dvec_new (coarse_ymir_mesh, 1,
                                              YMIR_GAUSS_NODE);
    }
    if (coarse_vel != NULL) {
      coarse_vel_int = coarse_vel;
    }
    else {
      coarse_vel_int = ymir_cvec_new (coarse_ymir_mesh, 3);
    }

    /* restrict velocity field from fine to coarse mesh */
    slabs_stokes_state_init_vel_press (coarse_state, coarse_ymir_mesh,
                                       coarse_press_elem);
    if (p_projection) {
      ymir_gmg_intergrid_p_restrict_cnode_accurate (
          fine_vel, coarse_state->vel_bc_vec);
    }
    else {
      ymir_gmg_intergrid_restrict_partition_cnode_accurate (
          fine_vel, upstream_interp_mangll,
          upstream_partition_mangll, coarse_state->vel_bc_vec);
    }
    ymir_vel_dir_separate (coarse_state->vel_bc_vec, NULL, NULL, NULL,
                           coarse_vel_dir);
    ymir_stokes_vec_set_velocity (coarse_state->vel_bc_vec,
                                  coarse_state->vel_press_vec,
                                  coarse_press_elem);
    ymir_vec_copy (coarse_state->vel_bc_vec, coarse_vel_int);

    /* store coarse velocity field in buffer (for next coarsening) */
    if (d->buffer_fine_vel != NULL) {
      ymir_vec_destroy (d->buffer_fine_vel);
    }
    d->buffer_fine_vel = ymir_vec_clone (coarse_state->vel_bc_vec);

    /* compute coefficient of Stokes operator on coarse mesh */
    if (d->init_nl_stokes_coeff) {
      slabs_physics_compute_init_nl_stokes_coeff (
          coarse_coeff, coarse_coeff_deriv_int, NULL, NULL, NULL,
          coarse_state, coarse_press_elem, coarse_vel_dir, physics_options);

#ifdef SL_COARSEN_STOKES_COEFF_VTK
      /* vtk output */
      ymir_vtk_write (coarse_ymir_mesh, filename, coarse_coeff, "coeff",
                      coarse_coeff_deriv_int, "coeff_deriv",
                      coarse_vel_int, "coarse_vel", NULL);
#endif
    }
    else {
      ymir_vec_t         *fine_bounds, *fine_yielding;
      ymir_vec_t         *coarse_bounds = ymir_dvec_new (coarse_ymir_mesh, 1,
                                                         YMIR_GAUSS_NODE);
      ymir_vec_t         *coarse_yielding = ymir_dvec_new (coarse_ymir_mesh, 1,
                                                           YMIR_GAUSS_NODE);

      /* set fine bounds and yielding markers */
      if (is_fine_level) { /* if at finest level in hierarchy */
        YMIR_ASSERT (physics_options->current_fine_bounds != NULL);
        YMIR_ASSERT (physics_options->current_fine_yielding != NULL);

        fine_bounds = ymir_vec_clone (physics_options->current_fine_bounds);
        fine_yielding = ymir_vec_clone (physics_options->current_fine_yielding);
        slabs_physics_coarsen_bounds_preprocess_fine (fine_bounds);
      }
      else { /* if at intermediate level in hierarchy */
        YMIR_ASSERT (d->buffer_fine_bounds != NULL);
        YMIR_ASSERT (d->buffer_fine_yielding != NULL);

        fine_bounds = d->buffer_fine_bounds;
        fine_yielding = d->buffer_fine_yielding;
      }

      /* pre-process fine markers */
      slabs_physics_coarsen_bounds_yielding_preprocess (fine_bounds,
                                                        fine_yielding);

      /* coarsen bounds and yielding marker */
      if (p_projection) {
        YMIR_ABORT_NOT_REACHED ();
        //TODO implement this
      }
      else {
        ymir_gmg_intergrid_restrict_partition_gauss_lin_simple (
            fine_bounds, upstream_interp_mangll,
            upstream_partition_mangll, coarse_bounds);
        ymir_gmg_intergrid_restrict_partition_gauss_lin_simple (
            fine_yielding, upstream_interp_mangll,
            upstream_partition_mangll, coarse_yielding);
      }

      /* post-process coarse markers */
      slabs_physics_coarsen_bounds_yielding_postprocess (coarse_bounds,
                                                         coarse_yielding);

      /* destroy fine markers */
      if (is_fine_level) {
        ymir_vec_destroy (fine_bounds);
        ymir_vec_destroy (fine_yielding);
      }

      /* store coarse markers in buffer (for next coarsening) */
      if (d->buffer_fine_bounds != NULL) {
        ymir_vec_destroy (d->buffer_fine_bounds);
      }
      d->buffer_fine_bounds = coarse_bounds;
      if (d->buffer_fine_yielding != NULL) {
        ymir_vec_destroy (d->buffer_fine_yielding);
      }
      d->buffer_fine_yielding = coarse_yielding;

      /* compute nonlinear coefficient */
      slabs_physics_compute_stokes_coeff_int (
          coarse_coeff, coarse_coeff_deriv_int, NULL,
          coarse_bounds, coarse_yielding,
          coarse_state, coarse_press_elem, coarse_vel_dir, physics_options, 1);

#ifdef SL_COARSEN_STOKES_COEFF_VTK
      /* vtk output */
      ymir_vtk_write (coarse_ymir_mesh, filename, coarse_coeff, "coeff",
                      coarse_coeff_deriv_int, "coeff_deriv",
                      coarse_vel_int, "coarse_vel",
                      coarse_bounds, "bounds_marker",
                      coarse_yielding, "yielding_marker", NULL);
#endif
    }

    /* destroy */
    if (coarse_coeff_deriv == NULL) {
      ymir_vec_destroy (coarse_coeff_deriv_int);
    }
    if (coarse_vel == NULL) {
      ymir_vec_destroy (coarse_vel_int);
    }
  }

  /*
   * Post-process Coefficient
   */

  switch (viscosity_coarsen_type) {
  case SL_VISCOSITY_COARSEN_EVAL:
    break;

  case SL_VISCOSITY_COARSEN_EVAL_PIECEWISE_DISCONT:
    ymir_dvec_set_const_inside_elems (coarse_coeff);
    if (coarse_coeff_deriv != NULL) {
      ymir_dvec_set_const_inside_elems (coarse_coeff_deriv);
    }
    break;

  default: /* unknown viscosity coarsening type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* destroy coarse Stokes state */
  slabs_stokes_state_destroy (coarse_state);

  /* increase counter */
  d->coarsen_count++;
}

/**
 * Computes right-hand side from temperature and background temperature.
 */
void
slabs_rhs_fn (double *rhs, const double x, const double y, const double z,
              const double temp, const double back_temp,
              slabs_physics_options_t *physics_options)
{
  const double        scaling = physics_options->rhs_scaling;

  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_BRICK:
    /**
     * Computes right-hand side from temperature and background temperature
     * on the cube domain:
     *
     *   f(x) = e_z * (T - T_0)
     *
     * where `e_z` is unit vector in z-direction, `T` is temperature,
     * `T_0` is background temperature.
     */
    rhs[0] = 0.0;
    rhs[1] = 0.0;
    rhs[2] = scaling * (temp - back_temp);
    break;

  case SL_DOMAIN_SHELL:
  case SL_DOMAIN_SHELL_CHUNK:
  case SL_DOMAIN_SHELL_SLICE:
    /**
     * Computes right-hand side from temperature and background temperature
     * on the shell domain:
     *
     *   f(x) = e_r * (T - T_0)
     *
     * where `e_r` is normalized spherical position vector, `T` is temperature,
     * `T_0` is background temperature.
     */
    {
      const double        radius = sqrt (x*x + y*y + z*z);

      rhs[0] = scaling * x / radius * (temp - back_temp);
      rhs[1] = scaling * y / radius * (temp - back_temp);
      rhs[2] = scaling * z / radius * (temp - back_temp);
    }
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }
}

/**
 * Computes right-hand side in an element.
 */
 void
slabs_rhs_elem (sc_dmatrix_t *rhs_el_mat,
                const double *x, const double *y, const double *z,
                const sc_dmatrix_t *temp_el_mat,
                slabs_physics_options_t *physics_options)
{
  const double       *temp_el_data = temp_el_mat->e[0];
  const int           n_nodes_per_el = temp_el_mat->m;
  int                 nodeid;

  /* check input */
  YMIR_ASSERT (temp_el_mat->m == rhs_el_mat->m);
  YMIR_ASSERT (temp_el_mat->n == 1);
  YMIR_ASSERT (rhs_el_mat->n == 3);

  /* compute right-hand side for this element */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    const double        temp = temp_el_data[nodeid];
    double              back_temp;
    double             *rhs = rhs_el_mat->e[nodeid];

    /* compute background temperature at this node */
    slabs_physics_background_temp_set_fn (&back_temp, x[nodeid], y[nodeid],
                                          z[nodeid], 0, physics_options);

    /* compute right-hand side at this node */
    slabs_rhs_fn (rhs, x[nodeid], y[nodeid], z[nodeid], temp, back_temp,
                  physics_options);
  }
}

typedef struct slabs_rhs_set_fn_data
{
  ymir_vec_t         *temp_vec;
  slabs_physics_options_t  *physics_options;
}
slabs_rhs_set_fn_data_t;

/**
 * Computes right-hand side at a node.
 */
static void
slabs_rhs_set_fn (double *rhs, double x, double y, double z,
                  ymir_locidx_t nid, void *data)
{
  slabs_rhs_set_fn_data_t  *d = (slabs_rhs_set_fn_data_t *) data;
  slabs_physics_options_t  *physics_options = d->physics_options;
  const double        temp = d->temp_vec->cvec->e[0][nid];
  double              back_temp;

  /* compute background temperature */
  slabs_physics_background_temp_set_fn (&back_temp, x, y, z, nid,
                                        physics_options);

  /* compute right-hand side */
  slabs_rhs_fn (rhs, x, y, z, temp, back_temp, physics_options);
}

/**
 * Computes right-hand side in (primal) function space.
 */
void
slabs_physics_compute_rhs_u_point (ymir_vec_t *rhs_u_point,
                                   slabs_stokes_state_t *state,
                                   slabs_physics_options_t *physics_options)
{
  const char         *this_fn_name = "slabs_physics_compute_rhs_u_point";
  const int           multiply_in_weak_zone =
                        physics_options->rhs_multiply_in_weak_zone;
  slabs_rhs_set_fn_data_t  data;

  /* check input parameters */
  YMIR_ASSERT (state->temperature != NULL);
  YMIR_ASSERT (state->temp_vec != NULL);

  /* set right-hand side */
  data.temp_vec = state->temp_vec;
  data.physics_options = physics_options;
  ymir_cvec_set_function (rhs_u_point, slabs_rhs_set_fn, &data);

  /* multiply in weak zone */
  if (multiply_in_weak_zone && state->weak_vec != NULL) {
    ymir_mesh_t        *mesh = rhs_u_point->mesh;
    ymir_vec_t         *weak_cvec = ymir_cvec_new (mesh, 1);
    ymir_vec_t         *inv_mass_lump = ymir_cvec_new (mesh, 1);

    YMIR_GLOBAL_INFOF ("%s: Multiply in weak zone to right-hand side\n",
                       this_fn_name);

    /* create inverse of lumped mass matrix */
    ymir_mass_lump (inv_mass_lump);
    ymir_vec_reciprocal (inv_mass_lump);
    ymir_vec_fabs (inv_mass_lump, inv_mass_lump);

    /* interpolate weak zone: discont. Gauss nodes -> cont. GLL nodes */
    ymir_mass_apply (state->weak_vec, weak_cvec);
    ymir_vec_multiply_in (inv_mass_lump, weak_cvec);

    /* modify right-hand side */
    ymir_cvec_multiply_in1 (weak_cvec, rhs_u_point);

    /* destroy */
    ymir_vec_destroy (weak_cvec);
    ymir_vec_destroy (inv_mass_lump);
  }
  else if (multiply_in_weak_zone && state->weak_vec == NULL) {
    YMIR_GLOBAL_LERRORF ("%s: ERROR: Multipling in weak zone failed!\n",
                         this_fn_name);
  }
}

/**
 *
 */
static void
slabs_physics_normal_boundary_stress_fn (double *stress_norm,
                                         double X, double Y, double Z,
                                         double nx, double ny, double nz,
                                         ymir_topidx_t face,
                                         ymir_locidx_t node_id,
                                         void *data)
{
  ymir_vec_t         *vec_bndr = (ymir_vec_t *) data;
  double             *v = ymir_cvec_index (vec_bndr, node_id, 0);

  YMIR_ASSERT (vec_bndr->ncfields == 3);

  /* compute inner product with boundary outer normal vector */
  *stress_norm = nx * v[0] + ny * v[1] + nz * v[2];
}

void
slabs_physics_compute_normal_boundary_stress (ymir_vec_t *stress_bndr_norm,
                                              ymir_vec_t *up,
                                              ymir_vec_t *rhs_u_point,
                                              ymir_stokes_op_t *stokes_op)
{
  ymir_mesh_t        *mesh = up->mesh;
  ymir_pressure_elem_t  *press_elem = stokes_op->press_elem;
  ymir_stress_op_t   *stress_op = stokes_op->stress_op;
  const int           skip_dir = stress_op->skip_dir;
  const ymir_topidx_t face_id = stress_bndr_norm->meshnum;

  ymir_vec_t         *rhs = ymir_stokes_vec_new (mesh, press_elem);
  ymir_vec_t         *residual_up = ymir_stokes_vec_new (mesh, press_elem);
  ymir_vec_t         *residual_u = ymir_cvec_new (mesh, 3);
  ymir_vec_t         *residual_bndr = ymir_face_cvec_new (mesh, face_id, 3);
  ymir_vec_t         *mass_lump_boundary;

  /* check input */
  YMIR_ASSERT_IS_CVEC (stress_bndr_norm);
  YMIR_ASSERT (stress_bndr_norm->ncfields == 1);
  YMIR_ASSERT (ymir_stokes_vec_is_stokes_vec (up));
  YMIR_ASSERT_IS_CVEC (rhs_u_point);
  YMIR_ASSERT (ymir_vec_is_not_dirty (up));
  YMIR_ASSERT (ymir_vec_is_not_dirty (rhs_u_point));

  /* construct the right-hand side */
  ymir_stokes_op_construct_rhs_ext (rhs_u_point, NULL, NULL, rhs,
                                    1 /* incompressible */, stokes_op);
  YMIR_ASSERT (sc_dmatrix_is_valid (rhs->dataown));
  YMIR_ASSERT (sc_dmatrix_is_valid (rhs->coff));
  YMIR_ASSERT (ymir_vec_is_not_dirty (rhs));

  /* turn off boundary constraints */
  stress_op->skip_dir = 1;

  /* compute (unconstrained) residual
   *   r_mom  = A * u + B^T * p - f
   *   r_mass = B * u
   */
  ymir_stokes_pc_apply_stokes_op (up, residual_up, stokes_op, 0, 0);
  ymir_vec_add (-1.0, rhs, residual_up);
  YMIR_ASSERT (sc_dmatrix_is_valid (residual_up->dataown));
  YMIR_ASSERT (sc_dmatrix_is_valid (residual_up->coff));
  ymir_vec_destroy (rhs);

  /* restore boundary constraints */
  stress_op->skip_dir = skip_dir;

  /* get the velocity component of the residual */
  ymir_stokes_vec_get_velocity (residual_up, residual_u, stokes_op->press_elem);

  /* interpolate residual onto boundary */
  ymir_interp_vec (residual_u, residual_bndr);
  YMIR_ASSERT (sc_dmatrix_is_valid (residual_bndr->dataown));
  YMIR_ASSERT (sc_dmatrix_is_valid (residual_bndr->coff));
  ymir_vec_destroy (residual_up);
  ymir_vec_destroy (residual_u);

  /* get the normal part of the residual */
  ymir_face_cvec_set_function (stress_bndr_norm,
                               slabs_physics_normal_boundary_stress_fn,
                               residual_bndr);
  ymir_vec_destroy (residual_bndr);

  /* invert mass matrix on boundary */
  mass_lump_boundary = ymir_face_cvec_new (mesh, face_id, 1);
  ymir_mass_lump (mass_lump_boundary);
  ymir_vec_divide_in (mass_lump_boundary, stress_bndr_norm);
  ymir_vec_destroy (mass_lump_boundary);
}

/**
 *
 */
static void
slabs_physics_separate_filter_values (ymir_vec_t *filter,
                                      const double threshold)
{
  /* separate filter values in case of cvec */
  if (ymir_vec_has_cvec (filter)) {
    ymir_face_mesh_t   *fmesh = &filter->mesh->fmeshes[filter->meshnum];
    const ymir_locidx_t n_cnodes = fmesh->Ncn;
    const int           n_fields = filter->ncfields;
    ymir_locidx_t       cnid;
    int                 fieldid;

    for (cnid = 0; cnid < n_cnodes; cnid++) {
      for (fieldid = 0; fieldid < n_fields; fieldid++) {
        double *val = ymir_cvec_index (filter, cnid, fieldid);

        if (threshold < *val) { /* if filter is including this node */
          *val = 1.0;
        }
        else { /* if filter is excluding this node */
          *val = 0.0;
        }
      }
    }
  }

  /* separate filter values in case of dvec */
  if (ymir_vec_has_dvec (filter)) {
    const ymir_locidx_t n_elements = filter->K;
    const int           n_nodes_per_el = filter->Np;
    const int           n_fields = filter->ndfields;
    ymir_locidx_t       elid;
    int                 nodeid, fieldid;

    for (elid = 0; elid < n_elements; elid++) {
      for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
        for (fieldid = 0; fieldid < n_fields; fieldid++) {
          double *val = ymir_dvec_index (filter, elid, nodeid, fieldid);

          if (threshold < *val) { /* if filter is including this node */
            *val = 1.0;
          }
          else { /* if filter is excluding this node */
            *val = 0.0;
          }
        }
      }
    }
  }

  /* separate filter values in case of evec */
  if (ymir_vec_has_evec (filter)) {
    YMIR_ABORT_NOT_REACHED (); //TODO implement this
  }
}

/**
 *
 */
static void
slabs_physics_add_filter (ymir_vec_t *filter_sum, ymir_vec_t *filter_add)
{
  /* compute pointwise: filter_sum += filter_add */
  ymir_vec_add (1.0, filter_add, filter_sum);
  slabs_physics_separate_filter_values (filter_sum, 0.0);
}

/**
 *
 */
static void
slabs_physics_multiply_filter (ymir_vec_t *filter_prod, ymir_vec_t *filter_mult)
{
  /* compute pointwise: filter_prod *= filter_mult */
  ymir_vec_multiply_in (filter_mult, filter_prod);
}

/**
 *
 */
static void
slabs_physics_invert_filter (ymir_vec_t *filter)
{
  /* compute pointwise: filter = - (filter - 1) */
  ymir_vec_shift (-1.0, filter);
  ymir_vec_scale (-1.0, filter);
}

/**
 *
 */
static double
slabs_physics_get_filter_volume (ymir_vec_t *filter)
{
  ymir_vec_t         *unit = ymir_vec_template (filter);
  ymir_vec_t         *mass_out = ymir_vec_template (filter);
  double              vol;

  /* set unit vector */
  ymir_vec_set_value (unit, 1.0);

  /* integrate to get volume */
  ymir_mass_apply (filter, mass_out);
  vol = ymir_vec_innerprod (unit, mass_out);

  /* destroy */
  ymir_vec_destroy (unit);
  ymir_vec_destroy (mass_out);

  /* return volume of filter */
  return vol;
}

/**
 *
 */
static double
slabs_physics_get_filtered_quantity_dmatrix_min (sc_dmatrix_t *quantity,
                                                 sc_dmatrix_t *filter)
{
  const sc_bint_t     totalsize = quantity->m * quantity->n;
  const double       *quantity_data = quantity->e[0];
  const double       *filter_data = filter->e[0];
  sc_bint_t           i;
  double              min = DBL_MAX;

  YMIR_ASSERT ((filter->m * filter->n) == totalsize );

  for (i = 0; i < totalsize; i++) {
    if (0.0 < filter_data[i]) {
      min = SC_MIN (min, quantity_data[i]);
    }
  }

  return min;
}

static double
slabs_physics_get_filtered_quantity_dmatrix_max (sc_dmatrix_t *quantity,
                                                 sc_dmatrix_t *filter)
{
  const sc_bint_t     totalsize = quantity->m * quantity->n;
  const double       *quantity_data = quantity->e[0];
  const double       *filter_data = filter->e[0];
  sc_bint_t           i;
  double              max = -DBL_MAX;

  YMIR_ASSERT ((filter->m * filter->n) == totalsize );

  for (i = 0; i < totalsize; i++) {
    if (0.0 < filter_data[i]) {
      max = SC_MAX (max, quantity_data[i]);
    }
  }

  return max;
}

static double
slabs_physics_get_filtered_quantity_min (ymir_vec_t *quantity,
                                         ymir_vec_t *filter)
{
  int                 mpiret;
  double              min_local, min_global;

  /* get local min */
  min_local = slabs_physics_get_filtered_quantity_dmatrix_min (
      quantity->dataown, filter->dataown);
  if (quantity->coff) {
    min_local = slabs_physics_get_filtered_quantity_dmatrix_min (
        quantity->dataown, filter->dataown);
  }

  /* get global min */
  mpiret = MPI_Allreduce (&min_local, &min_global, 1, MPI_DOUBLE, MPI_MIN,
                          quantity->mesh->ma->mpicomm);
  YMIR_CHECK_MPI (mpiret);

  return min_global;
}

static double
slabs_physics_get_filtered_quantity_max (ymir_vec_t *quantity,
                                         ymir_vec_t *filter)
{
  int                 mpiret;
  double              max_local, max_global;

  /* get local max */
  max_local = slabs_physics_get_filtered_quantity_dmatrix_max (
      quantity->dataown, filter->dataown);
  if (quantity->coff) {
    max_local = slabs_physics_get_filtered_quantity_dmatrix_max (
        quantity->dataown, filter->dataown);
  }

  /* get global max */
  mpiret = MPI_Allreduce (&max_local, &max_global, 1, MPI_DOUBLE, MPI_MAX,
                          quantity->mesh->ma->mpicomm);
  YMIR_CHECK_MPI (mpiret);

  return max_global;
}

/**
 *
 */
static void
slabs_physics_get_filtered_quantity_stats (double *min,
                                           double *max,
                                           double *mean,
                                           ymir_vec_t *quantity,
                                           ymir_vec_t *filter)
{
  ymir_vec_t         *unit = ymir_vec_template (quantity);
  ymir_vec_t         *mass_out = ymir_vec_template (quantity);

  /* compute min, max */
  *min = slabs_physics_get_filtered_quantity_min (quantity, filter);
  *max = slabs_physics_get_filtered_quantity_max (quantity, filter);

  /* set unit vector */
  ymir_vec_set_value (unit, 1.0);

  /* integrate filtered quantity */
  ymir_mass_apply (quantity, mass_out);
  *mean = ymir_vec_innerprod (unit, mass_out);

  /* divide by filtered volume */
  ymir_mass_apply (filter, mass_out);
  *mean *= 1.0 / ymir_vec_innerprod (unit, mass_out);

  /* destroy */
  ymir_vec_destroy (unit);
  ymir_vec_destroy (mass_out);
}

/**
 *
 */
void
slabs_physics_set_filter_lower_mantle (
                                     ymir_vec_t *filter,
                                     slabs_physics_options_t *physics_options)
{
  slabs_physics_set_filter_upper_mantle (filter, physics_options);
  slabs_physics_invert_filter (filter);
}

/**
 *
 */
void
slabs_physics_set_filter_upper_mantle (
                                     ymir_vec_t *filter,
                                     slabs_physics_options_t *physics_options)
{
  ymir_mesh_t        *mesh = filter->mesh;
  mangll_t           *mangll = mesh->ma;
  const int          *Vmask = mangll->refel->Vmask;
  const mangll_locidx_t  n_elements = mangll->mesh->K;
  const int           n_nodes_per_el = ymir_np (mangll->N);
  mangll_locidx_t     elid;
  sc_dmatrix_t       *filter_el_mat = sc_dmatrix_new (n_nodes_per_el,
                                                      filter->ndfields);
  double             *x, *y, *z;

  /* check input */
  YMIR_ASSERT_IS_DVEC (filter);

  /* initialize filter */
  ymir_dvec_set_zero (filter);
  sc_dmatrix_set_value (filter_el_mat, 1.0);

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    /* set element coordinates of physical space at GLL nodes */
    x = mangll->X->e[elid];
    y = mangll->Y->e[elid];
    z = mangll->Z->e[elid];

    /* set filter if element is in lower mantle */
    if (slabs_physics_elem_in_upper_mantle (x, y, z, Vmask, physics_options)) {
      ymir_dvec_set_elem (filter, filter_el_mat, YMIR_STRIDE_NODE, elid,
                          YMIR_SET);
    }
  }

  /* destroy */
  sc_dmatrix_destroy (filter_el_mat);
}

/**
 *
 */
void
slabs_physics_set_filter_lithosphere_layer (
                                     ymir_vec_t *filter,
                                     slabs_physics_options_t *physics_options)
{
  const double        layer_radius = SL_SHELL_RADIUS_TOP -
                                     200.0e3 / SL_EARTH_RADIUS; /* 200 km */
  ymir_mesh_t        *mesh = filter->mesh;
  mangll_t           *mangll = mesh->ma;
  const int          *Vmask = mangll->refel->Vmask;
  const mangll_locidx_t  n_elements = mangll->mesh->K;
  const int           n_nodes_per_el = ymir_np (mangll->N);
  mangll_locidx_t     elid;
  sc_dmatrix_t       *filter_el_mat = sc_dmatrix_new (n_nodes_per_el,
                                                      filter->ndfields);
  double             *x, *y, *z;

  /* check input */
  YMIR_ASSERT_IS_DVEC (filter);

  /* initialize filter */
  ymir_dvec_set_zero (filter);
  sc_dmatrix_set_value (filter_el_mat, 1.0);

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    /* set element coordinates of physical space at GLL nodes */
    x = mangll->X->e[elid];
    y = mangll->Y->e[elid];
    z = mangll->Z->e[elid];

    /* set filter if element is in lower mantle */
    if (layer_radius <=
        slabs_compute_radius_at_elem_center (x, y, z, Vmask, physics_options)) {
      ymir_dvec_set_elem (filter, filter_el_mat, YMIR_STRIDE_NODE, elid,
                          YMIR_SET);
    }
  }

  /* destroy */
  sc_dmatrix_destroy (filter_el_mat);
}

/**
 *
 */
void
slabs_physics_set_filter_plates (ymir_vec_t *filter,
                                 ymir_vec_t *viscosity_nondim,
                                 slabs_physics_options_t *physics_options)
{
  const double        visc_max_tol = 0.1;
  double              visc_max;

  /* check input */
  YMIR_ASSERT_IS_DVEC (filter);
  YMIR_ASSERT_IS_DVEC (viscosity_nondim);
  YMIR_ASSERT (filter->ndfields == viscosity_nondim->ndfields);
  YMIR_ASSERT (filter->node_type == viscosity_nondim->node_type);

  /* set max viscosity */
  if (0.0 < physics_options->viscosity_max) {
    visc_max = physics_options->viscosity_max;
  }
  else {
    visc_max = ymir_dvec_max_global (viscosity_nondim);
  }

  /* set filter for plates */
  ymir_vec_copy (viscosity_nondim, filter);
  slabs_physics_separate_filter_values (filter, visc_max * visc_max_tol);
}

/**
 *
 */
void
slabs_physics_set_filter_asthenosphere_away_from_subd (
                                     ymir_vec_t *filter,
                                     ymir_vec_t *viscosity_nondim,
                                     slabs_physics_options_t *physics_options)
{
  const double        subd_radius = 800.0e3 / SL_EARTH_RADIUS; /* 800 km */
  ymir_vec_t         *filter_asth = ymir_vec_template (filter);
  ymir_vec_t         *filter_no_plates = ymir_vec_template (filter);

  /* check input */
  YMIR_ASSERT_IS_DVEC (filter);

  /* set filter for asthenosphere layer */
  {
    ymir_vec_t         *filter_tmp = filter_no_plates;

    slabs_physics_set_filter_lithosphere_layer (filter_asth, physics_options);
    slabs_physics_invert_filter (filter_asth);

    slabs_physics_set_filter_upper_mantle (filter_tmp, physics_options);
    slabs_physics_multiply_filter (filter_asth, filter_tmp);
  }

  /* set filter for no plates */
  slabs_physics_set_filter_plates (filter_no_plates, viscosity_nondim,
                                   physics_options);
  slabs_physics_invert_filter (filter_no_plates);

  /* find distances to weak zones */
  {
    slabs_physics_weak_dist_set_fn_data_t  data;

    data.n_fields = filter->ndfields;
    data.physics_options = physics_options;
    ymir_dvec_set_function (filter, slabs_physics_weak_dist_set_fn, &data);
  }

  /* set filter for asthenosphere sufficiently far away from subduction zones */
  slabs_physics_separate_filter_values (filter, subd_radius);
  slabs_physics_multiply_filter (filter, filter_asth);
  slabs_physics_multiply_filter (filter, filter_no_plates);

  /* destroy */
  ymir_vec_destroy (filter_asth);
  ymir_vec_destroy (filter_no_plates);
}

/**
 *
 */
void
slabs_physics_set_filter_away_from_plate_bndr (
                                     ymir_vec_t *filter,
                                     slabs_physics_options_t *physics_options)
{
  double              weak_radius;
  slabs_physics_weak_dist_set_fn_data_t  data_cvec;
  slabs_physics_weak_dist_set_fn_data_t  data_dvec;

  /* set weak zone radius, which will be filtered out */
  if (physics_options->weakzone_type == SL_WEAKZONE_IMPORT_FILE) {
    weak_radius = physics_options->weakzone_import_thickness;
  }
  else {
    weak_radius = physics_options->weakzone_2plates_subdu_thickness;
  }
  weak_radius *= 5.0 / SL_EARTH_RADIUS;

  /* set data for callback function */
  data_cvec.n_fields = filter->ncfields;
  data_cvec.physics_options = physics_options;
  data_dvec.n_fields = filter->ndfields;
  data_dvec.physics_options = physics_options;

  /* find distances to weak zones */
  if (filter->meshnum == YMIR_VOL_MESH) { /* if volume vector */
    ymir_vec_set_function (
        filter,
        slabs_physics_weak_dist_set_fn, &data_cvec,
        slabs_physics_weak_dist_set_fn, &data_dvec,
        NULL, NULL);
  }
  else { /* if face vector */
    ymir_face_vec_set_function (
        filter,
        slabs_physics_weak_dist_set_face_fn, &data_cvec,
        slabs_physics_weak_dist_set_face_fn, &data_dvec,
        NULL, NULL);
  }

  /* set filter */
  slabs_physics_separate_filter_values (filter, weak_radius);
}

/**
 *
 */
void
slabs_physics_stats_quantity_in_lower_mantle (
                                     double *min, double *max, double *mean,
                                     ymir_vec_t *quantity,
                                     slabs_physics_options_t *physics_options)
{
  ymir_vec_t         *quantity_filtered = ymir_vec_clone (quantity);
  ymir_vec_t         *filter = ymir_vec_template (quantity);

  /* return if nothing to do */
  if (physics_options->viscosity_upper_mantle_radius <= 0.0) {
    *min = 0.0;
    *max = 0.0;
    *mean = 0.0;

    ymir_vec_destroy (quantity_filtered);
    ymir_vec_destroy (filter);
    return;
  }

  /* filter viscosity */
  slabs_physics_set_filter_lower_mantle (filter, physics_options);
  ymir_vec_multiply_in (filter, quantity_filtered);

  /* compute min, max, mean of filtered viscosity */
  slabs_physics_get_filtered_quantity_stats (min, max, mean,
                                             quantity_filtered, filter);

  /* destroy */
  ymir_vec_destroy (quantity_filtered);
  ymir_vec_destroy (filter);
}

/**
 *
 */
void
slabs_physics_stats_quantity_in_plates (
                                     double *min, double *max, double *mean,
                                     ymir_vec_t *quantity,
                                     ymir_vec_t *viscosity_nondim,
                                     slabs_physics_options_t *physics_options)
{
  ymir_vec_t         *quantity_filtered = ymir_vec_clone (quantity);
  ymir_vec_t         *filter = ymir_vec_template (quantity);

  /* filter viscosity */
  slabs_physics_set_filter_plates (filter, viscosity_nondim, physics_options);
  ymir_vec_multiply_in (filter, quantity_filtered);

  /* compute min, max, mean of filtered viscosity */
  slabs_physics_get_filtered_quantity_stats (min, max, mean,
                                             quantity_filtered, filter);

  /* destroy */
  ymir_vec_destroy (quantity_filtered);
  ymir_vec_destroy (filter);
}

/**
 *
 */
void
slabs_physics_stats_quantity_in_asthenosphere (
                                     double *min, double *max, double *mean,
                                     ymir_vec_t *quantity,
                                     ymir_vec_t *viscosity_nondim,
                                     slabs_physics_options_t *physics_options)
{
  ymir_vec_t         *quantity_filtered = ymir_vec_clone (quantity);
  ymir_vec_t         *filter = ymir_vec_template (quantity);

  /* filter viscosity */
  slabs_physics_set_filter_asthenosphere_away_from_subd (
      filter, viscosity_nondim, physics_options);
  ymir_vec_multiply_in (filter, quantity_filtered);

  /* compute min, max, mean of filtered viscosity */
  slabs_physics_get_filtered_quantity_stats (min, max, mean,
                                             quantity_filtered, filter);

  /* destroy */
  ymir_vec_destroy (quantity_filtered);
  ymir_vec_destroy (filter);
}

/**
 *
 */
void
slabs_physics_stats_velocity_away_from_plate_bndr (
                                     double *min, double *max, double *mean,
                                     ymir_vec_t *velocity,
                                     slabs_physics_options_t *physics_options)
{
  ymir_vec_t         *velocity_magn, *filter;

  /* check input */
  YMIR_ASSERT_IS_CVEC (velocity);

  /* compute mangitude of velocity */
  if (velocity->meshnum == YMIR_VOL_MESH) {
    velocity_magn = ymir_cvec_new (velocity->mesh, 1);
  }
  else {
    velocity_magn = ymir_face_cvec_new (velocity->mesh, velocity->meshnum, 1);
  }
  slabs_cvec_compute_magnitude (velocity, velocity_magn);

  /* filter out velocity at relevant points */
  filter = ymir_vec_template (velocity_magn);
  slabs_physics_set_filter_away_from_plate_bndr (filter, physics_options);
  ymir_vec_multiply_in (filter, velocity_magn);

  /* compute min, max, mean of filtered velocity */
  slabs_physics_get_filtered_quantity_stats (min, max, mean,
                                             velocity_magn, filter);

  /* destroy */
  ymir_vec_destroy (velocity_magn);
  ymir_vec_destroy (filter);
}

/**
 *
 */
double
slabs_physics_stats_plates_vol (ymir_vec_t *viscosity_nondim,
                                slabs_physics_options_t *physics_options)
{
  ymir_vec_t         *filter = ymir_vec_template (viscosity_nondim);
  double              vol;

  /* filter plates */
  slabs_physics_set_filter_plates (filter, viscosity_nondim, physics_options);

  /* compute volume of filter */
  vol = slabs_physics_get_filter_volume (filter);

  /* destroy */
  ymir_vec_destroy (filter);

  /* return volume of plates */
  return vol;
}

/**
 *
 */
double
slabs_physics_stats_yielding_vol (ymir_vec_t *yielding_marker)
{
  /* return volume of yielding regions */
  return slabs_physics_get_filter_volume (yielding_marker);
}

/**
 *
 */
double
slabs_physics_stats_plates_with_yielding_vol (
                                     ymir_vec_t *viscosity_nondim,
                                     ymir_vec_t *yielding_marker,
                                     slabs_physics_options_t *physics_options)
{
  ymir_vec_t         *filter = ymir_vec_template (viscosity_nondim);
  double              vol;

  /* filter plates and add yielding marker */
  slabs_physics_set_filter_plates (filter, viscosity_nondim, physics_options);
  slabs_physics_add_filter (filter, yielding_marker);

  /* compute volume of filter */
  vol = slabs_physics_get_filter_volume (filter);

  /* destroy */
  ymir_vec_destroy (filter);

  /* return volume of plates with yielding */
  return vol;
}

/**
 *
 */
void
slabs_physics_stats_print_all (ymir_vec_t *vel_press_nondim,
                               ymir_vec_t *viscosity_nondim,
                               ymir_vec_t *yielding_marker,
                               ymir_vec_t *rhs_u_point,
                               ymir_stokes_op_t *stokes_op,
                               slabs_physics_options_t *physics_options)
{
  const char         *this_fn_name = "slabs_physics_stats_print_all";
  const double        domain_vol = physics_options->domain_volume;
  ymir_mesh_t        *mesh = vel_press_nondim->mesh;
  ymir_pressure_elem_t  *press_elem = stokes_op->press_elem;
  ymir_vec_t         *u_dim = ymir_cvec_new (mesh, 3);
  ymir_vec_t         *p_dim = ymir_pressure_vec_new (mesh, press_elem);
  ymir_vec_t         *visc_dim = ymir_vec_clone (viscosity_nondim);
  ymir_vec_t         *strain_rate_dim = ymir_dvec_new (mesh, 1,
                                                       YMIR_GAUSS_NODE);
  ymir_vec_t         *stress_dim = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
  double              min, max, mean, vol;

  /* get velocity and pressure components */
  ymir_stokes_vec_get_components (vel_press_nondim, u_dim, p_dim, press_elem);

  /* compute 2nd invariant of the strain rate tensor (def. with sqrt) */
  {
    ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (
                                          mesh->ma->N, mesh->ma->ompsize);

    ymir_second_invariant_vec (u_dim, strain_rate_dim, vel_elem);
    ymir_dvec_sqrt (strain_rate_dim, strain_rate_dim);

    ymir_velocity_elem_destroy (vel_elem);
  }

  /* compute 2nd invariant of the stress tensor (def. with sqrt) */
  ymir_dvec_copy (strain_rate_dim, stress_dim);
  ymir_dvec_multiply_in (visc_dim, stress_dim);
  ymir_dvec_scale (2.0, stress_dim);

  /* tranform to physical dimensions */
  slabs_physics_transform_to_dimensional_velocity (u_dim);
  slabs_physics_transform_to_dimensional_pressure (p_dim);
  slabs_physics_transform_to_dimensional_viscosity (visc_dim);
  slabs_physics_transform_to_dimensional_strain_rate (strain_rate_dim);
  slabs_physics_transform_to_dimensional_stress (stress_dim);

  /* scale to geophysically practical values */
  ymir_vec_scale (100.0 * SL_SEC_PER_YEAR, u_dim);

  /*
   * Global Statistics
   */

  /* print velocity statistics */
  {
    ymir_vec_t         *u_magn = ymir_cvec_new (mesh, 1);

    slabs_cvec_compute_magnitude (u_dim, u_magn);
    min = ymir_vec_min_global (u_magn);
    max = ymir_vec_max_global (u_magn);
    YMIR_GLOBAL_STATISTICSF (
        "%s: Global velocity magn [cm/yr]:  min %.3e, max %.3e, max/min %.3e\n",
        this_fn_name, min, max, max/min);

    ymir_vec_destroy (u_magn);
  }

  /* print mean rotation of velocity */
  if (physics_options->domain_shape == SL_DOMAIN_SHELL) {
    double              rot_axis[3];

    ymir_stress_op_compute_mean_rotation (rot_axis, u_dim,
                                          stokes_op->stress_op, 0);
    YMIR_GLOBAL_STATISTICSF (
        "%s: Global velocity mean rot axis: x,y,z = %+.3e , %+.3e , %+.3e\n",
        this_fn_name, rot_axis[0], rot_axis[1], rot_axis[2]);
  }

  /* print pressure statistics */
  min = ymir_vec_min_global (p_dim);
  max = ymir_vec_max_global (p_dim);
  mean = ymir_pressure_vec_mean (p_dim, press_elem, 0);
  YMIR_GLOBAL_STATISTICSF (
      "%s: Global pressure [kg/m^3]:      (min %+.1e, max %+.1e), mean %.3e\n",
      this_fn_name, min, max, mean);

  /* print viscosity statistics */
  min = ymir_vec_min_global (visc_dim);
  max = ymir_vec_max_global (visc_dim);
  YMIR_GLOBAL_STATISTICSF (
      "%s: Global viscosity [Pa s]:       min %.3e, max %.3e, max/min %.3e\n",
      this_fn_name, min, max, max/min);

  /* print strain rate statistics */
  min = ymir_vec_min_global (strain_rate_dim);
  max = ymir_vec_max_global (strain_rate_dim);
  YMIR_GLOBAL_STATISTICSF (
      "%s: Global strain rate (II) [1/s]: min %.3e, max %.3e, max/min %.3e\n",
      this_fn_name, min, max, max/min);

  /* print stress statistics */
  min = ymir_vec_min_global (stress_dim);
  max = ymir_vec_max_global (stress_dim);
  YMIR_GLOBAL_STATISTICSF (
      "%s: Global stress (II) [Pa]:       min %.3e, max %.3e, max/min %.3e\n",
      this_fn_name, min, max, max/min);

  /* print statistics of normal stress at surface */
  if (rhs_u_point != NULL) {
    ymir_vec_t         *stress_norm_dim = ymir_face_cvec_new (mesh, SL_TOP, 1);
    ymir_vec_t         *unit = ymir_vec_template (stress_norm_dim);
    ymir_vec_t         *mass_out = ymir_vec_template (stress_norm_dim);

    slabs_physics_compute_normal_boundary_stress (
        stress_norm_dim, vel_press_nondim, rhs_u_point, stokes_op);
    slabs_physics_transform_to_dimensional_stress (stress_norm_dim);

    min = ymir_vec_min_global (stress_norm_dim);
    max = ymir_vec_max_global (stress_norm_dim);

    ymir_vec_set_value (unit, 1.0);
    ymir_mass_apply (stress_norm_dim, mass_out);
    mean = ymir_vec_innerprod (unit, mass_out);

    YMIR_GLOBAL_STATISTICSF (
        "%s: Global stress norm surf [Pa]:  min %+.2e, max %+.2e, mean %+.2e\n",
        this_fn_name, min, max, mean);

    ymir_vec_destroy (stress_norm_dim);
    ymir_vec_destroy (unit);
    ymir_vec_destroy (mass_out);
  }

  /*
   * Regional Statistics
   */

  /* print statistics of viscosity in lower mantle */
  slabs_physics_stats_quantity_in_lower_mantle (&min, &max, &mean,
                                                visc_dim, physics_options);
  YMIR_GLOBAL_STATISTICSF (
      "%s: Viscosity in LM [Pa s]:        min %.3e, max %.3e, mean %.3e\n",
      this_fn_name, min, max, mean);

  /* print statistics of viscosity in asthenosphere */
  slabs_physics_stats_quantity_in_asthenosphere (&min, &max, &mean,
                                                 visc_dim, viscosity_nondim,
                                                 physics_options);
  YMIR_GLOBAL_STATISTICSF (
      "%s: Viscosity in asthenosp [Pa s]: min %.3e, max %.3e, mean %.3e\n",
      this_fn_name, min, max, mean);

  /* print statistics of strain rate in plates */
  slabs_physics_stats_quantity_in_plates (&min, &max, &mean,
                                          strain_rate_dim, viscosity_nondim,
                                          physics_options);
  YMIR_GLOBAL_STATISTICSF (
      "%s: Strain rate in plates [1/s]:   min %.3e, max %.3e, mean %.3e\n",
      this_fn_name, min, max, mean);

  /* print statistics of stress in plates */
  slabs_physics_stats_quantity_in_plates (&min, &max, &mean,
                                          stress_dim, viscosity_nondim,
                                          physics_options);
  YMIR_GLOBAL_STATISTICSF (
      "%s: Stress in plates [Pa]:         min %.3e, max %.3e, mean %.3e\n",
      this_fn_name, min, max, mean);

  /* print statistics of surface velocity away from plate boundaries */
  {
    ymir_vec_t         *u_surf = ymir_face_cvec_new (mesh, SL_TOP, 3);

    ymir_interp_vec (u_dim, u_surf);
    slabs_physics_stats_velocity_away_from_plate_bndr (&min, &max, &mean,
                                                       u_surf, physics_options);
    YMIR_GLOBAL_STATISTICSF (
        "%s: Pl. velocity at surf [cm/yr]:  min %.3e, max %.3e, mean %.3e\n",
        this_fn_name, min, max, mean);

    ymir_vec_destroy (u_surf);
  }

  /*
   * Volume Statistics
   */

  /* print statistics of relative plates volume */
  vol = slabs_physics_stats_plates_vol (viscosity_nondim, physics_options);
  YMIR_GLOBAL_STATISTICSF (
      "%s: Volume of all plates:          rel %.7f\n",
      this_fn_name, vol/domain_vol);

  /* print statistics of relative yielding volume */
  if (yielding_marker != NULL) {
    vol = slabs_physics_stats_yielding_vol (yielding_marker);
    YMIR_GLOBAL_STATISTICSF (
        "%s: Volume of yielding regions:    rel %.7f\n",
        this_fn_name, vol/domain_vol);
  }

  /* print statistics of relative plates+yielding volume */
  if (yielding_marker != NULL) {
    vol = slabs_physics_stats_plates_with_yielding_vol (
        viscosity_nondim, yielding_marker, physics_options);
    YMIR_GLOBAL_STATISTICSF (
        "%s: Volume of plates+yielding:     rel %.7f\n",
        this_fn_name, vol/domain_vol);
  }

//TODO del
#if 0
  double              vel_subdu[3];
  double              vel_overr[3];

  /* compute range of plate velocities */
  slabs_physics_get_plate_velocity (vel_subdu, vel_overr, state,
                                    physics_options);

  if (!isnan (vel_subdu[0]) || !isnan (vel_subdu[1])) {
    YMIR_GLOBAL_STATISTICSF (
        "%s: Range of subdu plate vel: min %1.3e, max %1.3e, mean %1.3e\n",
        this_fn_name, vel_subdu[0], vel_subdu[1], vel_subdu[2]);
  }
  if (!isnan (vel_overr[0]) || !isnan (vel_overr[1])) {
    YMIR_GLOBAL_STATISTICSF (
        "%s: Range of overr plate vel: min %1.3e, max %1.3e, mean %1.3e\n",
        this_fn_name, vel_overr[0], vel_overr[1], vel_overr[2]);
  }
#endif

  /* destroy */
  ymir_vec_destroy (visc_dim);
  ymir_vec_destroy (u_dim);
  ymir_vec_destroy (p_dim);
  ymir_vec_destroy (strain_rate_dim);
  ymir_vec_destroy (stress_dim);
}

/**
 * TODO del
 */
typedef struct slabs_2plates_poly2_plate_velocity_set_fn_data
{
  double             *vel_magn_subdu;
  double             *vel_magn_overr;
  int                 n_nodes_subdu;
  int                 n_nodes_overr;
  ymir_face_mesh_t   *fmeshes;
  slabs_physics_options_t  *physics_options;
}
slabs_2plates_poly2_plate_velocity_set_fn_data_t;

static void
slabs_2plates_poly2_plate_velocity_set_fn (double *out, double x, double y,
                                           double z, double nx, double ny,
                                           double nz, ymir_topidx_t face,
                                           ymir_locidx_t node_id, void *data)
{
  slabs_2plates_poly2_plate_velocity_set_fn_data_t  *d =
    (slabs_2plates_poly2_plate_velocity_set_fn_data_t *) data;
  ymir_face_mesh_t   *fmesh = &(d->fmeshes[face]);
  slabs_physics_options_t  *physics_options = d->physics_options;
  const double        weak_subdu_longitude =
                        physics_options->weakzone_2plates_subdu_longitude;
  const double        weak_subdu_thickness =
                        physics_options->weakzone_2plates_subdu_thickness;
  const double        weak_ridge_width =
                        physics_options->weakzone_2plates_ridge_width;
  const double        weak_ridge_smoothwidth =
                        physics_options->weakzone_2plates_ridge_smoothwidth;
  double             *vel_magn_subdu = d->vel_magn_subdu;
  int                *n_nodes_subdu = &(d->n_nodes_subdu);
  double             *vel_magn_overr = d->vel_magn_overr;
  int                *n_nodes_overr = &(d->n_nodes_overr);
  double              lon, magn;

  /* compute longitude from cartesian coordinates */
  lon = slabs_compute_longitude (x, y, z, physics_options);

  /* compute magnitude of velocity */
  magn = sqrt (out[0]*out[0] + out[1]*out[1] + out[2]*out[2]);

  if (node_id < fmesh->owncount) { /* if node is owned by MPI process */
    /* set min, max, and add magnitude to average */
    if ( ( physics_options->domain_lon_min +
           (weak_ridge_width + weak_ridge_smoothwidth)/SL_EARTH_RADIUS ) < lon
         &&
         lon < (weak_subdu_longitude - weak_subdu_thickness/SL_EARTH_RADIUS)
       ) {
      /* if at subducting plate */
      vel_magn_subdu[0] = SC_MIN (vel_magn_subdu[0], magn);
      vel_magn_subdu[1] = SC_MAX (vel_magn_subdu[1], magn);
      vel_magn_subdu[2] += magn;
      (*n_nodes_subdu)++;
    }
    else if (  (weak_subdu_longitude + weak_subdu_thickness/SL_EARTH_RADIUS)
             < lon ) {
      /* if at overriding plate */
      vel_magn_overr[0] = SC_MIN (vel_magn_overr[0], magn);
      vel_magn_overr[1] = SC_MAX (vel_magn_overr[1], magn);
      vel_magn_overr[2] += magn;
      (*n_nodes_overr)++;
    }
  }
}

/**
 * TODO del
 */
static void
slabs_2plates_poly2_plate_velocity (double *vel_magn_subdu,
                                    double *vel_magn_overr,
                                    slabs_stokes_state_t *state,
                                    slabs_physics_options_t *physics_options)
{
  MPI_Comm            mpicomm = state->p8est->mpicomm;
  int                 mpisize = state->p8est->mpisize;
  int                 mpiret;

  ymir_mesh_t        *mesh = state->vel_press_vec->mesh;
  ymir_cvec_t        *u;
  ymir_evec_t        *p;
  ymir_cvec_t        *u_surface = ymir_face_cvec_new (mesh, SL_TOP, 3);
  slabs_2plates_poly2_plate_velocity_set_fn_data_t  data;
  double              vel_local[8];
  double             *vel_global;
  double              n_nodes_subdu;
  double              n_nodes_overr;
  int                 pid;

  /* get velocity and pressure from solution */
  slabs_stokes_vec_get_components_view (&u, &p, state->vel_press_vec);

  /* extract surface velocity */
  ymir_interp_vec (u, u_surface);

  /* initialize function data */
  vel_magn_subdu[0] = DBL_MAX;
  vel_magn_subdu[1] = 0.0;
  vel_magn_subdu[2] = 0.0;
  vel_magn_overr[0] = DBL_MAX;
  vel_magn_overr[1] = 0.0;
  vel_magn_overr[2] = 0.0;
  data.vel_magn_subdu = vel_magn_subdu;
  data.n_nodes_subdu = 0;
  data.vel_magn_overr = vel_magn_overr;
  data.n_nodes_overr = 0;
  data.fmeshes = mesh->fmeshes;
  data.physics_options = physics_options;

  /* compute min, max, and average velocity magnitude of subducting and
   * overriding plate */
  ymir_face_cvec_set_function (
      u_surface, slabs_2plates_poly2_plate_velocity_set_fn, &data);

  /* communicate local min, max, and mean velocities */
  vel_local[0] = vel_magn_subdu[0];
  vel_local[1] = vel_magn_subdu[1];
  vel_local[2] = vel_magn_subdu[2];
  vel_local[3] = (double) data.n_nodes_subdu;
  vel_local[4] = vel_magn_overr[0];
  vel_local[5] = vel_magn_overr[1];
  vel_local[6] = vel_magn_overr[2];
  vel_local[7] = (double) data.n_nodes_overr;
  vel_global = YMIR_ALLOC (double, 8 * mpisize);
  mpiret = MPI_Allgather (vel_local, 8, MPI_DOUBLE, vel_global, 8, MPI_DOUBLE,
                          mpicomm);
  YMIR_CHECK_MPI (mpiret);

  /* compute global min, max, and average velocity */
  vel_magn_subdu[0] = DBL_MAX;
  vel_magn_subdu[1] = 0.0;
  vel_magn_subdu[2] = 0.0;
  vel_magn_overr[0] = DBL_MAX;
  vel_magn_overr[1] = 0.0;
  vel_magn_overr[2] = 0.0;
  n_nodes_subdu = 0.0;
  n_nodes_overr = 0.0;
  for (pid = 0; pid < mpisize; pid++) {
    if (0.0 < vel_global[8 * pid + 3]) { /* if has nodes of subdu plate */
      vel_magn_subdu[0] = SC_MIN (vel_magn_subdu[0], vel_global[8 * pid]);
      vel_magn_subdu[1] = SC_MAX (vel_magn_subdu[1], vel_global[8 * pid + 1]);
      vel_magn_subdu[2] += vel_global[8 * pid + 2];
      n_nodes_subdu += vel_global[8 * pid + 3];
    }
    if (0.0 < vel_global[8 * pid + 7]) { /* if has nodes of overr plate */
      vel_magn_overr[0] = SC_MIN (vel_magn_overr[0], vel_global[8 * pid + 4]);
      vel_magn_overr[1] = SC_MAX (vel_magn_overr[1], vel_global[8 * pid + 5]);
      vel_magn_overr[2] += vel_global[8 * pid + 6];
      n_nodes_overr += vel_global[8 * pid + 7];
    }
  }
  vel_magn_subdu[2] /= n_nodes_subdu;
  vel_magn_overr[2] /= n_nodes_overr;

  /* destroy */
  ymir_vec_destroy (u);
  ymir_vec_destroy (p);
  ymir_vec_destroy (u_surface);
  YMIR_FREE (vel_global);
}

/**
 * TODO del
 */
void
slabs_physics_get_plate_velocity (double *vel_magn_subdu,
                                  double *vel_magn_overr,
                                  slabs_stokes_state_t *state,
                                  slabs_physics_options_t *physics_options)
{
  double              zero = 0.0;  /* no const to avoid warning */

  switch (physics_options->domain_shape) {
  case SL_DOMAIN_BRICK:
  case SL_DOMAIN_SHELL_SLICE:
    if (physics_options->weakzone_type == SL_WEAKZONE_2PLATES_POLY2) {
      slabs_2plates_poly2_plate_velocity (vel_magn_subdu, vel_magn_overr,
                                          state, physics_options);
    }
    else {
      /* otherwise set to NAN */
      vel_magn_subdu[0] = 0.0 / zero;
      vel_magn_subdu[1] = 0.0 / zero;
      vel_magn_subdu[2] = 0.0 / zero;
      vel_magn_overr[0] = 0.0 / zero;
      vel_magn_overr[1] = 0.0 / zero;
      vel_magn_overr[2] = 0.0 / zero;
    }
    break;

  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_SHELL:
  case SL_DOMAIN_SHELL_CHUNK:
    /* set to NAN */
    vel_magn_subdu[0] = 0.0 / zero;
    vel_magn_subdu[1] = 0.0 / zero;
    vel_magn_subdu[2] = 0.0 / zero;
    vel_magn_overr[0] = 0.0 / zero;
    vel_magn_overr[1] = 0.0 / zero;
    vel_magn_overr[2] = 0.0 / zero;
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }
}

/**
 *
 */
void
slabs_physics_transform_to_dimensional_temperature (ymir_vec_t * temp)
{
  YMIR_ABORT_NOT_REACHED (); //TODO implement this
}

/**
 *
 */
void
slabs_physics_transform_to_dimensional_viscosity (ymir_vec_t * visc)
{
  ymir_vec_scale (SL_VISC_REP, visc);
}

/**
 *
 */
void
slabs_physics_transform_to_dimensional_velocity (ymir_vec_t * vel)
{
  const double        dim_scal = SL_THERM_DIFFUS / SL_EARTH_RADIUS;

  ymir_vec_scale (dim_scal, vel);
}

/**
 *
 */
void
slabs_physics_transform_to_dimensional_strain_rate (ymir_vec_t * strain_rate)
{
  const double        dim_scal = SL_THERM_DIFFUS /
                                 (SL_EARTH_RADIUS * SL_EARTH_RADIUS);

  ymir_vec_scale (dim_scal, strain_rate);
}

/**
 *
 */
void
slabs_physics_transform_to_dimensional_stress (ymir_vec_t * stress)
{
  const double        dim_scal = SL_VISC_REP * SL_THERM_DIFFUS /
                                 (SL_EARTH_RADIUS * SL_EARTH_RADIUS);

  ymir_vec_scale (dim_scal, stress);
}

/**
 *
 */
void
slabs_physics_transform_to_dimensional_pressure (ymir_vec_t * press)
{
  const double        dim_scal = SL_VISC_REP * SL_THERM_DIFFUS /
                                 (SL_EARTH_RADIUS * SL_EARTH_RADIUS);

  ymir_vec_scale (dim_scal, press);
}

