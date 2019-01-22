#ifndef ADJOINT_GEOMETRY_H
#define ADJOINT_GEOMETRY_H

#include <rhea_base.h>
#include <adjoint_options.h>
#include <adjoint_physics.h>


/* 2plates_poly2 parameters */
#define SUBD_SLAB_CLOSEST_PT_NEWTON_RTOL (1.0e8 * SC_EPS)
#define SUBD_SLAB_CLOSEST_PT_NEWTON_MAXITER 40

/* Evaluates linear polynomial. */
#define subd_poly1(x,poly_coeff) ( (poly_coeff)[0] + (poly_coeff)[1]*(x) )

/* Evaluates quadratic polynomial.*/
#define subd_poly2(x,poly_coeff) \
  ( (poly_coeff)[0] + (poly_coeff)[1]*(x) + (poly_coeff)[2]*(x)*(x) )

/* Evaluates derivative of quadratic polynomial.*/
#define subd_poly2_deriv(x,poly_coeff) \
  ( (poly_coeff)[1] + 2*(poly_coeff)[2]*(x) )

/* enumerator for orientations of a point w.r.t. a curve (for 2plates) */
typedef enum
{
  SUBD_SLAB_CURVE_ORIENT_TOP,
  SUBD_SLAB_CURVE_ORIENT_BOTTOM,
  SUBD_SLAB_CURVE_ORIENT_BOTTOM_BACK,
  SUBD_SLAB_CURVE_ORIENT_TOP_RIGHT,
  SUBD_SLAB_CURVE_ORIENT_BOTTOM_LEFT
}
subd_subd_edge_orient_enum_t;

/* Computes interpolating polynomial in 2D via Hermite interpolation.*/
double *
subd_compute_poly2_interpolation (double start_node, double start_val,
                                   double start_deriv,
                                   double end_node, double end_val);

/* Runs Newton's method to find the closest point on the curve of a
 * piecewise linear/quadratic polynomial to a given point.*/
int
subd_closest_pt_newton (double *x, double point_x, double point_y,
                         double rtol, double maxiter,
                         double *poly1_left_coeff, double stitch_node_left,
                         double *poly2_coeff, double stitch_node_right,
                         double *poly1_right_coeff);

/* Computes distance and orientation of a point to a quadratic polynomial in 2D.*/
double *
subd_compute_closest_pt_on_poly2 (double point_x, double point_y,
                                   double *poly2_coeff, double start_node,
                                   double start_val, double start_deriv,
                                   double end_node, double end_val,
                                   int *orientation_wrt_curve);

/* Computes distance and orientation of a point to a quadratic polynomial in 2D.*/
double
subd_compute_rot_on_poly2 (double *poly2_coeff, double *closest_pt);

double
subd_compute_radius (const double x, const double y, const double z,
                      subd_options_t *subd_options);

double
subd_compute_longitude (const double x, const double y, const double z,
                         subd_options_t *subd_options);

void
subd_surface_location (subd_options_t *subd_options,
                       rhea_discretization_options_t *discr_options);

void
subd_coordx_set_fn (double *coordx,
                     double x, double y, double z,
                     double nx, double ny, double nz,
                     ymir_topidx_t face, ymir_locidx_t nid,
                     void *data);

void
subd_coordx_set_fn (double *coordx,
                     double x, double y, double z,
                     double nx, double ny, double nz,
                     ymir_topidx_t face, ymir_locidx_t nid,
                     void *data);


#endif /*SUBDUCTION_GEOMETRY_H*/
