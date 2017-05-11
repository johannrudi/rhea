/* Slabs_Weakzone Example:
 *
 * Cartesian domain.
*/

#include <rhea.h>
#include <ymir_velocity_vec.h>
#include <ymir_derivative_elem.h>
#include <ymir_vec_ops.h>
#include <ymir_stokes_pc.h>
#include <ymir_stokes_vec.h>
#include <ymir_stress_op.h>
#include <ymir_pressure_vec.h>

/* basic constants */
#define SLABS_SEC_PER_YEAR (31557600.0)    /* seconds in a year (365.25*24*3600) */
#define SLABS_EARTH_RADIUS (6371.0e3)      /* mean radius of the Earth [m]       */
#define SLABS_MANTLE_DEPTH (2871.0e3)      /* mean radius of the Earth [m]       */
#define SLABS_UPPER_MANTLE_DEPTH (660.0e3) /* approx. depth of upper mantle [m]  */
#define SLABS_THERM_DIFFUS (1.0e-6)        /* thermal diffusivity [m^2 / s]      */
#define SLABS_TEMP_DIFF (1400.0)           /* temperature difference [K]         */
#define SLABS_VISC_REP (1.0e20)            /* representative viscosity [Pa s]    */

/* shell parameters */
#define SLABS_SHELL_RADIUS_BOTTOM 0.55
#define SLABS_SHELL_RADIUS_TOP 1.0

/* temperature parameters */
#define SLABS_DEFAULT_CONST_TEMP 0.5

/* 2plates_poly2 parameters */
#define SLABS_CLOSEST_PT_NEWTON_RTOL (1.0e8 * SC_EPS)
#define SLABS_CLOSEST_PT_NEWTON_MAXITER 40

/* temperature 2plates_poly2 */
#define SLABS_DEFAULT_TEMP_TYPE_NAME "NONE"
#define SLABS_DEFAULT_TEMP_BACKGROUND_PLATE_AGE (50.0e6)    /* rhea1: 60 Myr */
#define SLABS_DEFAULT_TEMP_2PL_TRENCH_LONGITUDE (0.13)      /* rhea1: 0.13   */
#define SLABS_DEFAULT_TEMP_2PL_DIP_ANGLE (5.0)              /* rhea1: --     */
#define SLABS_DEFAULT_TEMP_2PL_SUBD_DEPTH (400.0e3)         /* rhea1: 400 km */
#define SLABS_DEFAULT_TEMP_2PL_SUBD_WIDTH (300.0e3)         /* rhea1: --     */
#define SLABS_DEFAULT_TEMP_2PL_SUBD_EDGE_WIDTH (1.0e3)      /* rhea1: --     */
#define SLABS_DEFAULT_TEMP_2PL_SUBD_EDGE_SMOOTHWIDTH (40.0e3)/* rhea1: --     */
#define SLABS_DEFAULT_TEMP_2PL_SUBD_PLATE_VELOCITY (4.0e-2)  /* rhea1: 4 cm/y */
#define SLABS_DEFAULT_TEMP_2PL_SUBD_PLATE_INITIAL_AGE (1.0e6)/* rhea1: --     */
#define SLABS_DEFAULT_TEMP_2PL_OVER_PLATE_AGE (40.0e6)       /* rhea1: 40 Myr */

/* weak zone 2plates_poly2 */
#define SLABS_WEAKZONE_2PLATES_SUBDU_LONGITUDE (-100.0)      /* rhea1:  0.13 */
#define SLABS_WEAKZONE_2PLATES_SUBDU_DIP_ANGLE (-1.0)        /* rhea1:   N/A */
#define SLABS_WEAKZONE_2PLATES_SUBDU_DEPTH (80.0e3)          /* rhea1: 50 km */
#define SLABS_WEAKZONE_2PLATES_SUBDU_WIDTH (-1.0)            /* rhea1:   N/A */
#define SLABS_WEAKZONE_2PLATES_SUBDU_THICKNESS (20.0e3)      /* rhea1: 20 km */
#define SLABS_WEAKZONE_2PLATES_SUBDU_THICKNESS_CONST (5.0e3) /* rhea1: 10 km */
#define SLABS_WEAKZONE_2PLATES_SUBDU_WEAK_FACTOR (1.0e-5)    /* rhea1:  1e-5 */
#define SLABS_WEAKZONE_2PLATES_RIDGE_DEPTH (30.0e3)          /* rhea1: 30 km */
#define SLABS_WEAKZONE_2PLATES_RIDGE_WIDTH (30.0e3)          /* rhea1: 30 km */
#define SLABS_WEAKZONE_2PLATES_RIDGE_SMOOTHWIDTH (10.0e3)    /* rhea1:  5 km */
#define SLABS_WEAKZONE_2PLATES_RIDGE_WEAK_FACTOR (1.0e-5)    /* rhea1:  1e-5 */

/* initialize slabs options: temperature */
double              temp_back_plate_age =
  SLABS_DEFAULT_TEMP_BACKGROUND_PLATE_AGE;
double              temp_2pl_trench_lon =
  SLABS_DEFAULT_TEMP_2PL_TRENCH_LONGITUDE;
double              temp_2pl_dip_angle = SLABS_DEFAULT_TEMP_2PL_DIP_ANGLE;
double              temp_2pl_subd_depth = SLABS_DEFAULT_TEMP_2PL_SUBD_DEPTH;
double              temp_2pl_subd_width = SLABS_DEFAULT_TEMP_2PL_SUBD_WIDTH;
double              temp_2pl_subd_edge_width =
  SLABS_DEFAULT_TEMP_2PL_SUBD_EDGE_WIDTH;
double              temp_2pl_subd_edge_smoothwidth =
  SLABS_DEFAULT_TEMP_2PL_SUBD_EDGE_SMOOTHWIDTH;
double              temp_2pl_subd_plate_vel =
  SLABS_DEFAULT_TEMP_2PL_SUBD_PLATE_VELOCITY;
double              temp_2pl_subd_plate_init_age =
  SLABS_DEFAULT_TEMP_2PL_SUBD_PLATE_INITIAL_AGE;
double              temp_2pl_over_plate_age =
  SLABS_DEFAULT_TEMP_2PL_OVER_PLATE_AGE;

/* initialize slabs options: weak zone */
double              weakzone_2pl_subdu_lon =
  SLABS_WEAKZONE_2PLATES_SUBDU_LONGITUDE;
double              weakzone_2pl_subdu_dip_angle =
  SLABS_WEAKZONE_2PLATES_SUBDU_DIP_ANGLE;
double              weakzone_2pl_subdu_depth = SLABS_WEAKZONE_2PLATES_SUBDU_DEPTH;
double              weakzone_2pl_subdu_width = SLABS_WEAKZONE_2PLATES_SUBDU_WIDTH;
double              weakzone_2pl_subdu_thickness =
  SLABS_WEAKZONE_2PLATES_SUBDU_THICKNESS;
double              weakzone_2pl_subdu_thickness_const =
  SLABS_WEAKZONE_2PLATES_SUBDU_THICKNESS_CONST;
double              weakzone_2pl_subdu_weak_factor =
  SLABS_WEAKZONE_2PLATES_SUBDU_WEAK_FACTOR;
double              weakzone_2pl_ridge_depth = SLABS_WEAKZONE_2PLATES_RIDGE_DEPTH;
double              weakzone_2pl_ridge_width = SLABS_WEAKZONE_2PLATES_RIDGE_WIDTH;
double              weakzone_2pl_ridge_smoothwidth =
  SLABS_WEAKZONE_2PLATES_RIDGE_SMOOTHWIDTH;
double              weakzone_2pl_ridge_weak_factor =
  SLABS_WEAKZONE_2PLATES_RIDGE_WEAK_FACTOR;

typedef struct slabs_domain_options
{
  double        x_min;
  double        x_max;
  double        y_min;
  double        y_max;
  double        z_min;
  double        z_max;
  double        lon_min;
  double        lon_max;
  double        radius_min;
  double        radius_max;
}
slabs_domain_options_t;


typedef struct slabs_2plates_poly2_geo_coeff
{
  double        start_node;
  double        start_val;
  double        end_node;
  double        end_val;
  double        start_deriv;
  double        *poly2_coeff;
}
slabs_2plates_poly2_geo_coeff_t;

/* enumerator for orientations of a point w.r.t. a curve (for 2plates) */
typedef enum
{
  SLABS_CURVE_ORIENT_TOP,
  SLABS_CURVE_ORIENT_BOTTOM,
  SLABS_CURVE_ORIENT_TOP_RIGHT,
  SLABS_CURVE_ORIENT_BOTTOM_LEFT
}
slabs_subd_edge_orient_enum_t;

/* struct for temperature options in slabs_options_t */
typedef struct slabs_temp_options
{
  double              temp_background_plate_age;
  double              temp_2plates_trench_longitude;
  double              temp_2plates_dip_angle;
  double              temp_2plates_subd_depth;
  double              temp_2plates_subd_width;
  double              temp_2plates_subd_edge_width;
  double              temp_2plates_subd_edge_smoothwidth;
  double              temp_2plates_subd_plate_velocity;
  double              temp_2plates_subd_plate_initial_age;
  double              temp_2plates_over_plate_age;
  slabs_2plates_poly2_geo_coeff_t *temp_2plates_geo_coeff;
}
slabs_temp_options_t;

/* enumerator for viscosity types */
typedef enum
{
  SLABS_VISC_ISOTROPY,
  SLABS_VISC_TRANSVERSELY_ISOTROPY
}
slabs_viscosity_anisotropy_t;

/* struct for viscosity options in slabs_options_t */
typedef struct slabs_visc_options
{
  slabs_viscosity_anisotropy_t  viscosity_anisotropy;
  double              viscosity_width;
  double              viscosity_factor;
  double              viscosity_TI_shear;
  double              viscosity_lith;
  double              viscosity_mantle;
}
slabs_visc_options_t;

/* struct for weak zone options in slabs_options_t */
typedef struct slabs_weak_options
{
  double              weakzone_2plates_subdu_longitude;  /* 2plates weak zone */
  double              weakzone_2plates_subdu_dip_angle;  /* at subduction */
  double              weakzone_2plates_subdu_depth;
  double              weakzone_2plates_subdu_width;
  double              weakzone_2plates_subdu_thickness;
  double              weakzone_2plates_subdu_thickness_const;
  double              weakzone_2plates_subdu_weak_factor;
  double              weakzone_2plates_ridge_depth;      /* 2plates weak zone */
  double              weakzone_2plates_ridge_width;      /* at ridge */
  double              weakzone_2plates_ridge_smoothwidth;
  double              weakzone_2plates_ridge_weak_factor;
  slabs_2plates_poly2_geo_coeff_t *weak_2plates_geo_coeff;
}
slabs_weak_options_t;

/* enumerator for boundary conditions */
typedef enum
{
  SLABS_VEL_DIR_BC_INOUTFLOW_SIN,
  SLABS_VEL_DIR_BC_INOUTFLOW_TANH,
}
slabs_vel_dir_bc_t;

/* struct for velocity boundary condition options in slabs_options_t */
typedef struct slabs_velbc_options
{
  slabs_vel_dir_bc_t  vel_dir_bc;
  double              flow_scale;

  double              vel_dir_bc_upper;
  double              vel_dir_bc_lower;
}
slabs_velbc_options_t;

/* options of slabs example */
typedef struct slabs_options
{
  slabs_domain_options_t   * slabs_domain_options;
  slabs_temp_options_t   * slabs_temp_options;
  slabs_visc_options_t   * slabs_visc_options;
  slabs_weak_options_t   * slabs_weak_options;
  slabs_velbc_options_t  * slabs_velbc_options;
}
slabs_options_t;

/****************************************************
 * Curvature Slabs 2plate_poly2
*****************************************************/

/* Computes the radius of a shell domain or the corresponding value for a
 * rectangular domain. */
double
slabs_compute_radius (const double x, const double y, const double z,
                      slabs_options_t *slabs_options)
{
      slabs_domain_options_t    * dom_opt = slabs_options->slabs_domain_options;
      const double    z_max = dom_opt->z_max;
      const double    radius_min = dom_opt->radius_min;
      const double    radius_max = dom_opt->radius_max;

      return z / z_max * (radius_max - radius_min) + radius_min;
}

/* Computes the longitude of a shell domain or the corresponding value for a
 * rectangular domain. */
double
slabs_compute_longitude (const double x, const double y, const double z,
                         slabs_options_t *slabs_options)
{
      slabs_domain_options_t    * dom_opt = slabs_options->slabs_domain_options;
      const double        lon_min = dom_opt->lon_min;
      const double        lon_max = dom_opt->lon_min;
      const double        y_max = dom_opt->y_max;

      return y / y_max * (lon_max - lon_min) + lon_min;
}


/* Evaluates linear polynomial. */
#define slabs_poly1(x,poly_coeff) ( (poly_coeff)[0] + (poly_coeff)[1]*(x) )

/* Evaluates quadratic polynomial.*/
#define slabs_poly2(x,poly_coeff) \
  ( (poly_coeff)[0] + (poly_coeff)[1]*(x) + (poly_coeff)[2]*(x)*(x) )

/* Evaluates derivative of quadratic polynomial.*/
#define slabs_poly2_deriv(x,poly_coeff) \
  ( (poly_coeff)[1] + 2*(poly_coeff)[2]*(x) )

/* Computes interpolating polynomial in 2D via Hermite interpolation.*/
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

/* Runs Newton's method to find the closest point on the curve of a
 * piecewise linear/quadratic polynomial to a given point.*/
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

/* Computes distance and orientation of a point to a quadratic polynomial in 2D.*/
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

    *orientation_wrt_curve = SLABS_CURVE_ORIENT_BOTTOM;
  }
  else if (line_orientation_end < 0) { /* if pt. is right of line through end pt
                                        * with dir. normal to end deriv. */
    proj = - (end_node - point_x) * tangent_end[0]
           - (end_val - point_y) * tangent_end[1];
    closest_pt[0] = end_node + proj * tangent_end[0];
    closest_pt[1] = end_val + proj * tangent_end[1];

    if (point_y <= slabs_poly1 (point_x, poly1_right_coeff)) { /* if pt below */
      *orientation_wrt_curve = SLABS_CURVE_ORIENT_BOTTOM_LEFT;
    }
    else { /* if point above linear poly. */
      *orientation_wrt_curve = SLABS_CURVE_ORIENT_TOP_RIGHT;
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
                             SLABS_CLOSEST_PT_NEWTON_RTOL,
                             SLABS_CLOSEST_PT_NEWTON_MAXITER,
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

    /* set orientation */
    if (point_y <= slabs_poly2 (point_x, poly2_coeff)) { /* if pt. below curve*/
      *orientation_wrt_curve = SLABS_CURVE_ORIENT_BOTTOM;
    }
    else { /* if pt. above curve */
      *orientation_wrt_curve = SLABS_CURVE_ORIENT_TOP;
    }
  }

  /* return distance */
  return closest_pt;
}

/**************************************
 * Temperature Computation
 *************************************/
double
slabs_temperature_brick_2plates_poly2 (double r, double lon,
                                       slabs_options_t * slabs_options)
{
  double              trench_lon;
  double              subd_depth, subd_width;
  double              subd_edge_width, subd_edge_smoothwidth;
  double              subd_plate_vel, subd_plate_init_age;
  double              over_plate_age;

  double              y_min = slabs_options->slabs_domain_options->y_min;
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

  slabs_temp_options_t * temp_options = slabs_options->slabs_temp_options;
  slabs_2plates_poly2_geo_coeff_t *geo = temp_options->temp_2plates_geo_coeff;

  /* set parameters according to temperature options */
  trench_lon = temp_options->temp_2plates_trench_longitude;
  subd_depth = temp_options->temp_2plates_subd_depth;
  subd_width = temp_options->temp_2plates_subd_width;
  subd_edge_width = temp_options->temp_2plates_subd_edge_width;
  subd_edge_smoothwidth = temp_options->temp_2plates_subd_edge_smoothwidth;
  subd_plate_vel = temp_options->temp_2plates_subd_plate_velocity;
  subd_plate_init_age = temp_options->temp_2plates_subd_plate_initial_age;
  over_plate_age = temp_options->temp_2plates_over_plate_age;

  /* check parameters */
  RHEA_ASSERT (0.0 < trench_lon);
  RHEA_ASSERT (0.0 < subd_depth && 0.0 < subd_width);
  RHEA_ASSERT (0.0 <= subd_edge_width && 0.0 <= subd_edge_smoothwidth);
  RHEA_ASSERT (0.0 < subd_plate_vel);
  RHEA_ASSERT (0.0 <= subd_plate_init_age);
  RHEA_ASSERT (0.0 < over_plate_age);

  start_node = geo->start_node;
  start_val = geo->start_val;
  start_deriv = geo->start_deriv;
  end_node = geo->end_node;
  end_val = geo->end_val;
  poly2_coeff = geo->poly2_coeff;

  /* compute closest point on curve and orientation w.r.t. curve */
  closest_pt = slabs_compute_closest_pt_on_poly2 (lon, r,
                                                  poly2_coeff, start_node,
                                                  start_val, start_deriv,
                                                  end_node, end_val,
                                                  &orientation_wrt_curve);


  /* compute "fake" (TODO) distance */
  dist = SLABS_MANTLE_DEPTH * sqrt ( (lon - closest_pt[0]) * (lon - closest_pt[0])
                                  + (r - closest_pt[1]) * (r - closest_pt[1]) );

  /* For computing the temperature, use temperature profile of the plates
   * from halfspace cooling model, erf(z/(2*sqrt(t*kappa))), as base
   * temperature field. */

  /* compute temperature of plates at surface */
  if (orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM) { /* if subd plate */
    /* calculate age of subducting plate */
    ridge_dist = fabs (lon - y_min) * SLABS_MANTLE_DEPTH;
    vel = subd_plate_vel / SLABS_SEC_PER_YEAR;
    plate_time = ridge_dist / vel + subd_plate_init_age * SLABS_SEC_PER_YEAR;

    /* avoid division by zero */
    plate_time = SC_MAX (1.0e-3, plate_time);
  }
  else  { /* if overriding plate */
    /* calculate age of overriding plate */
    plate_time = over_plate_age * SLABS_SEC_PER_YEAR;
  }
  depth = (SLABS_SHELL_RADIUS_TOP - r) * SLABS_MANTLE_DEPTH;
  temp = erf ( depth / (2.0 * sqrt (plate_time * SLABS_THERM_DIFFUS)) );

  /* compute temperature of the plate subducted inside of mantle */
  if (   orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM
      || orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_LEFT ) {
    /* calculate total lenght of subducted plate */
    subd_length = SLABS_MANTLE_DEPTH * sqrt (
                    (start_node - end_node) * (start_node - end_node)
                    + (start_val - end_val) * (start_val - end_val) );

    /* calculate age of plate */
    if (lon <= trench_lon) {
      ridge_dist = fabs (lon - y_min) * SLABS_MANTLE_DEPTH;
    }
    else {
      ridge_dist = fabs (trench_lon - y_min) * SLABS_MANTLE_DEPTH
                   + subd_length * depth / subd_depth;
    }
    vel = subd_plate_vel / SLABS_SEC_PER_YEAR;
    plate_time = ridge_dist / vel + subd_plate_init_age * SLABS_SEC_PER_YEAR;

    /* avoid division by zero */
    plate_time = SC_MAX (1.0e-3, plate_time);

    /* compute temperature of subducting plate */
    temp_subd = erf ( dist / (2.0 * sqrt (plate_time * SLABS_THERM_DIFFUS)) );
  }

  /* smooth the subducting plate's top edge */
  if (   orientation_wrt_curve == SLABS_CURVE_ORIENT_TOP
      || orientation_wrt_curve == SLABS_CURVE_ORIENT_TOP_RIGHT ) {
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
  if (   orientation_wrt_curve == SLABS_CURVE_ORIENT_TOP_RIGHT
      || orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_LEFT ) {
    /* compute distance between tip and closest point on curve */
    dist = SLABS_MANTLE_DEPTH * sqrt (
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

  /* return temperature */
  YMIR_ASSERT (isfinite (temp));
  return temp;
}

void
slabs_temperature_set_fn (double *temp, double x, double y, double z,
                          ymir_locidx_t nid, void *data)
{
  slabs_options_t    *slabs_options = (slabs_options_t *) data;
  double              lon, r;

  /* compute radius and longitude */
  lon = y;
  r = z;

  *temp = slabs_temperature_brick_2plates_poly2 (r, lon, slabs_options);

  /* check temperature for `nan` and `inf` */
  YMIR_ASSERT (isfinite (*temp));

  /* bound temperature to valid interval */
  *temp = SC_MIN (1.0, *temp);
  *temp = SC_MAX (0.0, *temp);
}

static void
slabs_temperature_compute (ymir_vec_t *temperature,
                           slabs_options_t *slabs_options)
{
  const char         *this_fn_name = "slabs_temperature_compute";
  double              trench_lon;
  double              dip_angle;
  double              subd_depth, subd_width;
  double              start_node, start_val, start_deriv, end_node, end_val;
  double              *poly2_coeff;
  slabs_2plates_poly2_geo_coeff_t geo;
  slabs_temp_options_t * temp_options = slabs_options->slabs_temp_options;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* set parameters according to physics options */
  trench_lon = temp_options->temp_2plates_trench_longitude;
  dip_angle = temp_options->temp_2plates_dip_angle;
  subd_depth = temp_options->temp_2plates_subd_depth;
  subd_width = temp_options->temp_2plates_subd_width;

  /* check parameters */
  RHEA_ASSERT (0.0 < trench_lon);
  RHEA_ASSERT (0.0 < dip_angle && dip_angle < 90.0);
  RHEA_ASSERT (0.0 < subd_depth && 0.0 < subd_width);

  /* set points for polynomial interpolation */
  start_node = trench_lon;
  start_val = SLABS_SHELL_RADIUS_TOP;
  start_deriv = tan (-dip_angle / 180.0 * M_PI);
  end_node = start_node + subd_width / SLABS_MANTLE_DEPTH;
  end_val = start_val - subd_depth / SLABS_MANTLE_DEPTH;
  /* compute interpolating quadratic polynomial */
  poly2_coeff = slabs_compute_poly2_interpolation (start_node, start_val,
                                                   start_deriv,
                                                   end_node, end_val);

  RHEA_GLOBAL_PRODUCTIONF("temperature: brick_2plates_poly2 downgoing slabs geometry\nstart_node: (%f %f), start_deriv: %f, end_node: (%f %f)\npoly2_coeff: (%f %f %f)\n", start_node, start_val, start_deriv, end_node, end_val, poly2_coeff[2], poly2_coeff[1], poly2_coeff[0]);

  geo.start_node = start_node;
  geo.start_val = start_val;
  geo.start_deriv = start_deriv;
  geo.end_node = end_node;
  geo.end_val = end_val;
  geo.poly2_coeff = poly2_coeff;
  temp_options->temp_2plates_geo_coeff = &geo;

  ymir_cvec_set_function (temperature,
                          slabs_temperature_set_fn,
                          slabs_options);
  YMIR_FREE (poly2_coeff);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**************************************
 * Weakzone Computation
 *************************************/

/* Computes the weak zone factor depending on the distance to a weak zone.
 * The edges of the weak zone are smoothed by a Gaussian.
 * 1 - (1 - weak_factor) * exp ( - dist^2 / (2 * (0.5*thickness)^2) ) */
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

/* Computes distance to weak zone between plates. */
 double
slabs_weakzone_subduct_dist_2plates_poly2 (double r, double lon,
                                         double *poly2_coeff,
                                         double start_node,
                                         double start_val,
                                         double start_deriv,
                                         double end_node,
                                         double end_val)
{
  int                 orientation_wrt_curve;
  double             *closest_pt;
  double              dist;

  /* compute closest point on curve and orientation w.r.t. curve */
  closest_pt = slabs_compute_closest_pt_on_poly2 (lon, r,
                                                  poly2_coeff, start_node,
                                                  start_val, start_deriv,
                                                  end_node, end_val,
                                                  &orientation_wrt_curve);

  /* compute distance */
  if (   orientation_wrt_curve == SLABS_CURVE_ORIENT_TOP
      || orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM ) {
    /* compute distance to closest point on curve */
    dist = SLABS_MANTLE_DEPTH * sqrt (  SC_SQR (lon - closest_pt[0])
                                   + SC_SQR (r - closest_pt[1]) );
  }
  else if (   orientation_wrt_curve == SLABS_CURVE_ORIENT_TOP_RIGHT
           || orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_LEFT ) {
    /* compute distance to tip of curve */
    dist = SLABS_MANTLE_DEPTH * sqrt (  SC_SQR (end_node - lon)
                                   + SC_SQR (end_val - r) );
  }
  else {
    double              zero = 0.0;  /* no const to avoid warning */
    dist = 1.0 / zero;
    YMIR_ABORT_NOT_REACHED ();
  }

  /* destroy */
  YMIR_FREE (closest_pt);

  /* return distance to weak zone curve */
  return dist;
}

/* Computes distance to weak zone at ridge. */
 double
slabs_weakzone_ridge_dist_2plates_poly2 (double r, double lon,
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
      dist = SLABS_MANTLE_DEPTH * (lon - end_node);
    }
    else if (lon <= end_node && r < end_val) {
      /* compute distance to bottom edge */
      dist = SLABS_MANTLE_DEPTH * (end_val - r);
    }
    else {
      /* compute distance to bottom right corner */
      dist = SLABS_MANTLE_DEPTH * sqrt (  SC_SQR (lon - end_node)
                                     + SC_SQR (end_val - r) );
    }
  }

  /* return distance to ridge weak zone */
  return dist;
}

/* Computes the weak fault zone for the shell slice with two plates. */
static double
slabs_weakzone_brick_2plates_poly2 (double r, double lon,
                              slabs_options_t * slabs_options)
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
  double              y_min = slabs_options->slabs_domain_options->y_min;
  double              start_node, start_val, start_deriv, end_node, end_val;
  double              *poly2_coeff;
  double              dist;
  double              weak = 1.0;

  slabs_weak_options_t *weak_options = slabs_options->slabs_weak_options;
  slabs_2plates_poly2_geo_coeff_t *geo = weak_options->weak_2plates_geo_coeff;

  /* set parameters according to weakzone options */
  subdu_lon = weak_options->weakzone_2plates_subdu_longitude;
  subdu_dip_angle = weak_options->weakzone_2plates_subdu_dip_angle;
  subdu_depth = weak_options->weakzone_2plates_subdu_depth;
  subdu_width = weak_options->weakzone_2plates_subdu_width;
  subdu_thickness = weak_options->weakzone_2plates_subdu_thickness;
  subdu_thickness_const =
    weak_options->weakzone_2plates_subdu_thickness_const;
  subdu_weak_factor = weak_options->weakzone_2plates_subdu_weak_factor;
  ridge_depth = weak_options->weakzone_2plates_ridge_depth;
  ridge_width = weak_options->weakzone_2plates_ridge_width;
  ridge_smoothwidth = weak_options->weakzone_2plates_ridge_smoothwidth;
  ridge_weak_factor = weak_options->weakzone_2plates_ridge_weak_factor;

  /* check parameters */
  YMIR_ASSERT (0.0 < subdu_lon);
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

  start_node = geo->start_node;
  start_val = geo->start_val;
  start_deriv = geo->start_deriv;
  end_node = geo->end_node;
  end_val = geo->end_val;
  poly2_coeff = geo->poly2_coeff;

  /* only consider point in a rectangle containing the weak zone */
  courtesy_width = subdu_thickness / SLABS_MANTLE_DEPTH;
  total_thickness = (2.0 * subdu_thickness - subdu_thickness_const)
                    / SLABS_MANTLE_DEPTH;
  if (   (  start_node - 0.5 * total_thickness
          / sin (subdu_dip_angle / 180.0 * M_PI) - courtesy_width ) <= lon
      && lon <= (end_node + 0.5 * total_thickness + courtesy_width)
      && (end_val - total_thickness - courtesy_width) <= r ) {
    /* compute distance to subduction weak zone */
    dist = slabs_weakzone_subduct_dist_2plates_poly2 (r, lon, poly2_coeff,
                                                      start_node, start_val, start_deriv,
                                                      end_node, end_val);

    /* set weak zone factor */
    weak = slabs_weakzone_factor_fn (dist, subdu_thickness,
                                     subdu_thickness_const, subdu_weak_factor);
  }

  /*
   * set ridge weak zone in the left corner of domain, so that subducting
   * plate can "move"
   */

  /* set bottom left corner of weak zone */
  end_node = y_min + ridge_width / SLABS_MANTLE_DEPTH;
  end_val = SLABS_SHELL_RADIUS_TOP - ridge_depth / SLABS_MANTLE_DEPTH;

  /* only consider points close to weak zone */
  courtesy_width = 2.0 * ridge_smoothwidth / SLABS_MANTLE_DEPTH;
  if (lon <= (end_node + courtesy_width) && (end_val - courtesy_width) <= r) {
    /* compute distance to ridge weak zone */
    dist = slabs_weakzone_ridge_dist_2plates_poly2 (r, lon, end_node, end_val);

    /* set weak zone factor */
    weak = slabs_weakzone_factor_fn (dist, ridge_smoothwidth, 0.0,
                                     ridge_weak_factor);
  }

  /*
   * return weak zone factor
   */
  return weak;
}

static double
slabs_weakzone_node (const double x, const double y, const double z,
                     slabs_options_t *slabs_options)
{
  double              lon;
  double              r;

  /* compute radius and longitude */
  lon = y;
  r = z;

  /* compute weak zone factor */
  return slabs_weakzone_brick_2plates_poly2 (r, lon, slabs_options);
}

/* compute weak zone factor of an element */
void
slabs_weakzone_elem (double *_sc_restrict weak_elem,
                     const double *x,
                     const double *y,
                     const double *z,
                     const int n_nodes_per_el,
                     slabs_options_t *slabs_options)
{
  int            nodeid;

  /* compute weak zone factor for each node */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) { /* loop over all
                                                         * nodes */
    weak_elem[nodeid] = slabs_weakzone_node (x[nodeid], y[nodeid], z[nodeid],
                                             slabs_options);
  }
}

void
slabs_weakzone_compute (ymir_dvec_t *weakzone, slabs_options_t *slabs_options)
{
  const char         *this_fn_name = "slabs_weakzone_compute";
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (weakzone);
  const ymir_locidx_t n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  sc_dmatrix_t       *weak_el_mat;
  double             *weak_el_data;
  double             *x, *y, *z, *tmp_el;
  ymir_locidx_t      elid;

  double              start_node, start_val, start_deriv, end_node, end_val;
  double              subdu_lon;
  double              subdu_dip_angle;
  double              subdu_depth, subdu_width;
  double              *poly2_coeff;
  slabs_2plates_poly2_geo_coeff_t geo;
  slabs_weak_options_t *weak_options = slabs_options->slabs_weak_options;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* set parameters according to weakzone options */
  subdu_lon = weak_options->weakzone_2plates_subdu_longitude;
  subdu_dip_angle = weak_options->weakzone_2plates_subdu_dip_angle;
  subdu_depth = weak_options->weakzone_2plates_subdu_depth;
  subdu_width = weak_options->weakzone_2plates_subdu_width;

  /* check parameters */
  RHEA_ASSERT (0.0 < subdu_lon);
  RHEA_ASSERT (0.0 < subdu_dip_angle);
  RHEA_ASSERT (0.0 < subdu_depth && 0.0 < subdu_width);

  /* set points for polynomial interpolation */
  start_node = subdu_lon;
  start_val = SLABS_SHELL_RADIUS_TOP;
  start_deriv = tan (-subdu_dip_angle / 180.0 * M_PI);
  end_node = start_node + subdu_width / SLABS_MANTLE_DEPTH;
  end_val = start_val - subdu_depth / SLABS_MANTLE_DEPTH;
  /* compute interpolating quadratic polynomial */
  poly2_coeff = slabs_compute_poly2_interpolation (start_node, start_val,
                                                   start_deriv,
                                                   end_node, end_val);

  RHEA_GLOBAL_PRODUCTIONF("weakzone: brick_2plates_poly2 downgoing slabs geometry\nstart_node: (%f %f), start_deriv: %f, end_node: (%f %f)\npoly2_coeff: (%f %f %f)\n", start_node, start_val, start_deriv, end_node, end_val, poly2_coeff[2], poly2_coeff[1], poly2_coeff[0]);

  geo.start_node = start_node;
  geo.start_val = start_val;
  geo.start_deriv = start_deriv;
  geo.end_node = end_node;
  geo.end_val = end_val;
  geo.poly2_coeff = poly2_coeff;
  weak_options->weak_2plates_geo_coeff = &geo;

  /* check input */
  RHEA_ASSERT (weakzone->ndfields == 1);

  /* create work variables */
  weak_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  weak_el_data = weak_el_mat->e[0];
  x = RHEA_ALLOC (double, n_nodes_per_el);
  y = RHEA_ALLOC (double, n_nodes_per_el);
  z = RHEA_ALLOC (double, n_nodes_per_el);
  tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    /* get coordinates of this element at Gauss nodes */
    ymir_mesh_get_elem_coord_gauss (x, y, z, elid, mesh, tmp_el);

    /* compute weak zone factor */
    slabs_weakzone_elem (weak_el_data, x, y, z,
                         n_nodes_per_el, slabs_options);

    /* set weak zone of this element */
    rhea_viscosity_set_elem_gauss (weakzone, weak_el_mat, elid);
  }

  /* destroy */
  sc_dmatrix_destroy (weak_el_mat);
  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);
  YMIR_FREE (poly2_coeff);  //TO FIND OUT: if I use RHEA_FREE, there is memory leak stderr

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**************************************
 * Viscosity Computation
 *************************************/


/******************************************************
 *  Transversely Isotropic Viscosity
*******************************************************/

/* get shear viscosity for transversely isotropy model from stress operator*/
void
slabs_stress_op_copy_shear_visc (ymir_vec_t *shear_visc,
                                   ymir_stress_op_t *stress_op)
{
  /* check input */
  RHEA_ASSERT (stress_op->coeff_TI_svisc != NULL);

  /* copy Stokes coefficient and divide by 2 */
  ymir_vec_copy (stress_op->coeff_TI_svisc, shear_visc);
  ymir_vec_scale (0.5, shear_visc);
}

/* get TI tensor for transversely isotropy model from stress operator*/
void
slabs_stress_op_copy_TI_tensor (ymir_vec_t *TI_tensor,
                                   ymir_stress_op_t *stress_op)
{
  /* check input */
  RHEA_ASSERT (stress_op->coeff_TI_tensor != NULL);

  ymir_vec_copy (stress_op->coeff_TI_tensor, TI_tensor);
}

/* Computes distance and orientation of a point to a quadratic polynomial in 2D.*/
static double
slabs_compute_rot_on_poly2 (double *poly2_coeff, double *closet_pt)
{
  double              tangent[2];
  double              horiz[2] = {1.0, .0};
  double              rot;

  /* compute vectors of curve at start and end points */
  tangent[1] = slabs_poly2_deriv (closet_pt[0], poly2_coeff);
  tangent[0] = 1.0 / sqrt (1.0 + tangent[1] * tangent[1]);
  tangent[1] *= tangent[0];

  rot = acos(horiz[0]*tangent[0] + horiz[1]*tangent[1]);

  /* return distance */
  return rot;
}

/* Computes the rotation from x axis of weak fault zone with two plates. */
static double
slabs_TI_rotation_brick_2plates_poly2 (double r, double lon,
                                       slabs_options_t * slabs_options)
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
  double              y_min = slabs_options->slabs_domain_options->y_min;
  double              start_node, start_val, start_deriv, end_node, end_val;
  double              *poly2_coeff;
  double              rot = 0.0;
  int                 orientation_wrt_curve;
  double             *closest_pt, dist=1.0;

  slabs_weak_options_t *weak_options = slabs_options->slabs_weak_options;
  slabs_2plates_poly2_geo_coeff_t *geo = weak_options->weak_2plates_geo_coeff;

  /* set parameters according to weakzone options */
  subdu_lon = weak_options->weakzone_2plates_subdu_longitude;
  subdu_dip_angle = weak_options->weakzone_2plates_subdu_dip_angle;
  subdu_depth = weak_options->weakzone_2plates_subdu_depth;
  subdu_width = weak_options->weakzone_2plates_subdu_width;
  subdu_thickness = weak_options->weakzone_2plates_subdu_thickness;
  subdu_thickness_const =
    weak_options->weakzone_2plates_subdu_thickness_const;
  ridge_depth = weak_options->weakzone_2plates_ridge_depth;
  ridge_width = weak_options->weakzone_2plates_ridge_width;
  ridge_smoothwidth = weak_options->weakzone_2plates_ridge_smoothwidth;

  /* check parameters */
  YMIR_ASSERT (0.0 < subdu_lon);
  YMIR_ASSERT (0.0 < subdu_dip_angle);
  YMIR_ASSERT (0.0 < subdu_depth && 0.0 < subdu_width);
  YMIR_ASSERT (0.0 < subdu_thickness);
  YMIR_ASSERT (subdu_thickness_const <= subdu_thickness);
  YMIR_ASSERT (0.0 < ridge_depth && 0.0 < ridge_width);
  YMIR_ASSERT (0.0 <= ridge_smoothwidth);

  /*
   * set subduction weak zone between plates
   */
  start_node = geo->start_node;
  start_val = geo->start_val;
  start_deriv = geo->start_deriv;
  end_node = geo->end_node;
  end_val = geo->end_val;
  poly2_coeff = geo->poly2_coeff;

  /* only consider point in a rectangle containing the weak zone */
  courtesy_width = subdu_thickness / SLABS_MANTLE_DEPTH;
  total_thickness = (2.0 * subdu_thickness - subdu_thickness_const)
                    / SLABS_MANTLE_DEPTH;
  if (   (  start_node - 0.5 * total_thickness
          / sin (subdu_dip_angle / 180.0 * M_PI) - courtesy_width ) <= lon
      && lon <= (end_node + 0.5 * total_thickness + courtesy_width)
      && (end_val - total_thickness - courtesy_width) <= r ) {
    /*
     * compute rotation of subduction weak zone from y=0 axis
     * */
    /* compute closest point on curve and orientation w.r.t. curve */
    closest_pt = slabs_compute_closest_pt_on_poly2 (lon, r,
                                                  poly2_coeff, start_node,
                                                  start_val, start_deriv,
                                                  end_node, end_val,
                                                  &orientation_wrt_curve);

    dist = slabs_weakzone_subduct_dist_2plates_poly2 (r, lon, poly2_coeff,
                                                      start_node, start_val, start_deriv,
                                                      end_node, end_val);
    if (dist < 0.5 * total_thickness *  SLABS_MANTLE_DEPTH)  {
      rot = slabs_compute_rot_on_poly2 (poly2_coeff, closest_pt);
    }
    else {
      rot = 0.0;
    }

    YMIR_FREE (closest_pt);
  }

  /*
   * set ridge weak zone in the left corner of domain, so that subducting
   * plate can "move"
   */

  /* set bottom left corner of weak zone */
  end_node = y_min + ridge_width / SLABS_MANTLE_DEPTH;
  end_val = SLABS_SHELL_RADIUS_TOP - ridge_depth / SLABS_MANTLE_DEPTH;

  /* only consider points close to weak zone */
  courtesy_width = 2.0 * ridge_smoothwidth / SLABS_MANTLE_DEPTH;
  if (lon <= (end_node + courtesy_width) && (end_val - courtesy_width) <= r) {
    /* compute rotation of ridge weak zone */
//TODO:    rot = slabs_ridge_rot_2plates_poly2 (r, lon, end_node, end_val);
    rot = .0;  // now assume the rotation is 0, horizontal
  }

  /*
   * return weak zone factor
   */
  return rot;
}

static double
slabs_TI_rotation_node (const double x, const double y, const double z,
                        slabs_options_t *slabs_options)
{
  double              lon;
  double              r;

  /* compute radius and longitude */
  lon = y;
  r = z;

  /* compute weak zone factor */
  return slabs_TI_rotation_brick_2plates_poly2 (r, lon, slabs_options);
}

/* compute weak zone factor of an element */
void
slabs_TI_rotation_elem (double *_sc_restrict rot_elem,
                        double *_sc_restrict weak_elem,
                         const double *x,
                         const double *y,
                         const double *z,
                         const int n_nodes_per_el,
                         slabs_options_t *slabs_options)
{
  int            nodeid;

  /* compute weak zone factor for each node */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    if (fabs(1.0 - weak_elem[nodeid]) < SC_EPS)  {
      rot_elem[nodeid] = .0;
    }
    else {
      rot_elem[nodeid] = slabs_TI_rotation_node (x[nodeid], y[nodeid], z[nodeid],
                                                 slabs_options);
    }
  }
}

static void
slabs_TI_rotation_compute (ymir_vec_t *rotate, ymir_vec_t *weak,
                           slabs_options_t *slabs_options)
{
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (rotate);
  const ymir_locidx_t n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  sc_dmatrix_t       *rot_el_mat, *weak_el_mat;
  double             *rot_el_data, *weak_el_data;
  double             *x, *y, *z, *tmp_el;
  ymir_locidx_t      elid;

  double              start_node, start_val, start_deriv;
  double              end_node, end_val;
  double              subdu_lon;
  double              subdu_dip_angle;
  double              subdu_depth, subdu_width;
  double              *poly2_coeff;
  slabs_2plates_poly2_geo_coeff_t geo;
  slabs_weak_options_t *weak_options = slabs_options->slabs_weak_options;
  const char         *this_fn_name = "slabs_TI_rotation_compute";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* set parameters according to weakzone options */
  subdu_lon = weak_options->weakzone_2plates_subdu_longitude;
  subdu_dip_angle = weak_options->weakzone_2plates_subdu_dip_angle;
  subdu_depth = weak_options->weakzone_2plates_subdu_depth;
  subdu_width = weak_options->weakzone_2plates_subdu_width;

  /* check parameters */
  RHEA_ASSERT (0.0 < subdu_lon);
  RHEA_ASSERT (0.0 < subdu_dip_angle);
  RHEA_ASSERT (0.0 < subdu_depth && 0.0 < subdu_width);

  /*
   * set subduction weak zone between plates
   */

  /* set points for polynomial interpolation */
  start_node = subdu_lon;
  start_val = SLABS_SHELL_RADIUS_TOP;
  start_deriv = tan (-subdu_dip_angle / 180.0 * M_PI);
  end_node = start_node + subdu_width / SLABS_MANTLE_DEPTH;
  end_val = start_val - subdu_depth / SLABS_MANTLE_DEPTH;
  /* compute interpolating quadratic polynomial */
  poly2_coeff = slabs_compute_poly2_interpolation (start_node, start_val,
                                                   start_deriv,
                                                   end_node, end_val);

  geo.start_node = start_node;
  geo.start_val = start_val;
  geo.start_deriv = start_deriv;
  geo.end_node = end_node;
  geo.end_val = end_val;
  geo.poly2_coeff = poly2_coeff;
  weak_options->weak_2plates_geo_coeff = &geo;

  /* check input */
  RHEA_ASSERT (rotate->ndfields == 1);

  /* create work variables */
  rot_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  rot_el_data = rot_el_mat->e[0];
  weak_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  x = RHEA_ALLOC (double, n_nodes_per_el);
  y = RHEA_ALLOC (double, n_nodes_per_el);
  z = RHEA_ALLOC (double, n_nodes_per_el);
  tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    /* get coordinates of this element at Gauss nodes */
    ymir_mesh_get_elem_coord_gauss (x, y, z, elid, mesh, tmp_el);

    weak_el_data = rhea_viscosity_get_elem_gauss (weak_el_mat, weak, elid);

    /* compute weak zone factor */
    slabs_TI_rotation_elem (rot_el_data, weak_el_data, x, y, z,
                            n_nodes_per_el, slabs_options);

    /* set weak zone of this element */
    rhea_viscosity_set_elem_gauss (rotate, rot_el_mat, elid);
  }

  /* destroy */
  sc_dmatrix_destroy (rot_el_mat);
  sc_dmatrix_destroy (weak_el_mat);
  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);
  YMIR_FREE (poly2_coeff);
}

/* Computes shear viscosity and rotation angle.*/
static void
slabs_TI_viscosity_compute ( ymir_mesh_t *ymir_mesh,  ymir_vec_t *TI_svisc,
                            ymir_vec_t *viscosity,
                            ymir_vec_t *weakzone,
                            slabs_options_t *slabs_options)
{
  ymir_vec_multiply (viscosity, weakzone, TI_svisc);
}

/* setup the TI shear viscosity and tensor in stress operator */
void
slabs_stokes_problem_setup_TI (ymir_mesh_t *ymir_mesh,
                               rhea_stokes_problem_t *stokes_problem,
                               slabs_options_t *slabs_options,
                               ymir_vec_t *coeff_TI_svisc,
                               ymir_vec_t *TI_rotate)
{
  const char         *this_fn_name = "slabs_stokes_problem_setup_TI";
  ymir_stokes_op_t   *stokes_op;
  ymir_stress_op_t   *stress_op;
  ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);
  ymir_vec_t         *weakzone = rhea_viscosity_new (ymir_mesh);

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* copy viscosity */
  rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

  /* compute the shear viscosity and rotation angles */
  slabs_weakzone_compute (weakzone, slabs_options);
  slabs_TI_viscosity_compute (ymir_mesh, coeff_TI_svisc, viscosity, weakzone, slabs_options);

  ymir_vec_scale (2.0, coeff_TI_svisc);

  slabs_TI_rotation_compute (TI_rotate, weakzone, slabs_options);

  /* get the viscous stress operator */
  stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
  stress_op = stokes_op->stress_op;

 /* update viscous stress operator providing the anisotropic viscosity */
  ymir_stress_op_coeff_compute_TI_tensor (stress_op, coeff_TI_svisc,
                                          TI_rotate);
  /* destroy */
  rhea_viscosity_destroy (viscosity);
  rhea_viscosity_destroy (weakzone);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**************************************
 * Output Functions
 *************************************/

/* Write vtk of input data. */
slabs_write_input (ymir_mesh_t *ymir_mesh,
                   rhea_stokes_problem_t *stokes_problem,
                   rhea_temperature_options_t *temp_options,
                   ymir_vec_t *temperature,
                   ymir_vec_t *weakzone,
                   ymir_vec_t *visc_TI_svisc,
                   ymir_vec_t *visc_TI_rotate,
                   const char *vtk_write_input_path)
{
  const char         *this_fn_name = "slabs_write_input";
  ymir_vec_t         *background_temp = rhea_temperature_new (ymir_mesh);
  ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);
  ymir_vec_t         *rhs_vel;
  char                path[BUFSIZ];

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);
  rhs_vel = rhea_stokes_problem_get_rhs_vel (stokes_problem);

  rhea_temperature_background_compute (background_temp, temp_options);


  rhea_vtk_write_input_data (vtk_write_input_path, temperature,
                             background_temp, weakzone, viscosity, NULL,
                             rhs_vel);

  rhea_temperature_destroy (background_temp);
  rhea_viscosity_destroy (viscosity);

  /* output transversely isotropy parameters*/
  {
    if (visc_TI_svisc != NULL && visc_TI_rotate != NULL) {
      ymir_vec_t         *shear_visc = ymir_vec_clone (visc_TI_svisc);
      ymir_vec_scale (0.5, shear_visc);
      ymir_vec_t *angle = ymir_vec_clone (visc_TI_rotate);
      ymir_vec_scale (180.0/M_PI, angle);
      snprintf (path, BUFSIZ, "%s_anisotropic_viscosity", vtk_write_input_path);
      ymir_vtk_write (ymir_mesh, path,
                      shear_visc, "shear_viscosity",
                      angle, "rotation",
                      NULL);
      rhea_viscosity_destroy (shear_visc);
      rhea_viscosity_destroy (angle);
    }
  }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/***********************************************************
 * Post-processing for 2nd invariant stress, traction, .etc.
 ***********************************************************/

/* compute traction as well as normal/shear stress at each element*/
void
slabs_2inv_stress_TI_elem (sc_dmatrix_t * in, sc_dmatrix_t * out,
                            sc_dmatrix_t * visc, sc_dmatrix_t * svisc,
                            sc_dmatrix_t * TItens, ymir_velocity_elem_t * vel_elem,
                            double *_sc_restrict rxd, double *_sc_restrict sxd,
                            double *_sc_restrict txd, double *_sc_restrict ryd,
                            double *_sc_restrict syd, double *_sc_restrict tyd,
                            double *_sc_restrict rzd, double *_sc_restrict szd,
                            double *_sc_restrict tzd, sc_dmatrix_t * drst, sc_dmatrix_t * brst)
{
  const int           N = ymir_n (vel_elem->N);
  const int           Np = ymir_n (vel_elem->Np);
  int                 gp, i, j, k, l;
  double              S[9], E[9];
  double              cs, ct;
  double             *_sc_restrict viscd  = visc->e[0];
  double             *_sc_restrict sviscd = svisc->e[0];
  double             *_sc_restrict outd   = out->e[0];

  sc_dmatrix_t       *tempvec1 = vel_elem->tempvec1;
  sc_dmatrix_t       *tempvec2 = vel_elem->tempvec2;
  sc_dmatrix_t       *tempvec3 = vel_elem->tempvec3;
  sc_dmatrix_t       *temptens = vel_elem->temptens1;
  ymir_derivative_elem_grad (N, 3, drst->e[0], brst->e[0],
                             rxd, sxd, txd, ryd, syd, tyd,
                             rzd, szd, tzd, in, temptens, tempvec1, tempvec2,
                             tempvec3, 0);

  /* create stress from duvw/dxyz * viscosity */
  for (gp = 0; gp < Np; gp++) {
    double               outv = .0;
    double               *_sc_restrict temptensd = temptens->e[0] + 9 * gp;
    double               *_sc_restrict T   = TItens->e[0] + 9 * gp;

    S[0] = temptensd[0];
    S[1] = (temptensd[1] + temptensd[3]) * (1. / 2.);
    S[2] = (temptensd[2] + temptensd[6]) * (1. / 2.);
    S[3] = S[1];
    S[4] = temptensd[4];
    S[5] = (temptensd[5] + temptensd[7]) * (1. / 2.);
    S[6] = S[2];
    S[7] = S[5];
    S[8] = temptensd[8];

    for (k = 0; k < 9; k++) {
      E[k] = 0.0;
    }
    for (k=0; k<3; k++) {
      for (l=0; l<3; l++) {
        E[0] += T[    k] * S[3*k + l] * T[3*l    ];
        E[1] += T[    k] * S[3*k + l] * T[3*l + 1];
        E[2] += T[    k] * S[3*k + l] * T[3*l + 2];
        E[4] += T[3 + k] * S[3*k + l] * T[3*l + 1];
        E[5] += T[3 + k] * S[3*k + l] * T[3*l + 2];
        E[8] += T[6 + k] * S[3*k + l] * T[3*l + 2];
      }
    }
    E[3] = E[1];
    E[6] = E[2];
    E[7] = E[5];

    cs = viscd[gp] + sviscd[gp];
    ct = viscd[gp] - sviscd[gp];

    /* compute linear combination of viscous stress and 3x3 tensor */
    S[0] = cs * S[0] + ct * E[0];
    S[1] = cs * S[1] + ct * E[1];
    S[2] = cs * S[2] + ct * E[2];
    S[3] = cs * S[3] + ct * E[3];
    S[4] = cs * S[4] + ct * E[4];
    S[5] = cs * S[5] + ct * E[5];
    S[6] = cs * S[6] + ct * E[6];
    S[7] = cs * S[7] + ct * E[7];
    S[8] = cs * S[8] + ct * E[8];

    for (k = 0; k < 9; k++) {
      outv += SC_SQR (S[k]);
    }
    outd[gp] = sqrt (0.5 * outv);
  }
}

/* compute traction as well as normal and shear stress*/
void
slabs_2inv_stress_TI (ymir_cvec_t * vel, ymir_dvec_t * tauII,
                       ymir_dvec_t * visc, ymir_dvec_t *svisc,
                       ymir_dvec_t * TItens, ymir_velocity_elem_t * vel_elem)
{
  ymir_mesh_t         *mesh = vel->mesh;
  ymir_locidx_t       elid;
  const int           N  = ymir_n (mesh->cnodes->N);
  const int           Np = ymir_np (mesh->cnodes->N);
  const int           K  = mesh->cnodes->K;
  sc_dmatrix_t       *elemin = sc_dmatrix_new (1, 3 * Np);
  sc_dmatrix_t       *elemout = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemvisc = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemsvisc = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemTItens = sc_dmatrix_new (1, 9 * Np);

  ymir_dvec_set_zero (tauII);

  for (elid = 0; elid < K; elid++)  {
    double             *_sc_restrict rxd = mesh->drdx->e[elid];
    double             *_sc_restrict sxd = mesh->dsdx->e[elid];
    double             *_sc_restrict txd = mesh->dtdx->e[elid];
    double             *_sc_restrict ryd = mesh->drdy->e[elid];
    double             *_sc_restrict syd = mesh->dsdy->e[elid];
    double             *_sc_restrict tyd = mesh->dtdy->e[elid];
    double             *_sc_restrict rzd = mesh->drdz->e[elid];
    double             *_sc_restrict szd = mesh->dsdz->e[elid];
    double             *_sc_restrict tzd = mesh->dtdz->e[elid];
    double             *_sc_restrict Jdetd = mesh->Jdet->e[elid];

    ymir_cvec_get_elem_interp (vel, elemin, YMIR_STRIDE_NODE, elid,
                               YMIR_GLL_NODE, YMIR_COPY);
    ymir_dvec_get_elem_interp (visc, elemvisc, YMIR_STRIDE_COMP, elid,
                              YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (svisc, elemsvisc, YMIR_STRIDE_COMP, elid,
                              YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (TItens, elemTItens, YMIR_STRIDE_COMP, elid,
                              YMIR_GAUSS_NODE, YMIR_READ);


    slabs_2inv_stress_TI_elem (elemin, elemout,
                                elemvisc, elemsvisc, elemTItens,
                                vel_elem, rxd, sxd, txd, ryd,
                                syd, tyd, rzd, szd, tzd, mesh->drst, mesh->brst);

    ymir_dvec_set_elem_interp (tauII, elemout, YMIR_STRIDE_COMP, elid,
                             YMIR_GAUSS_NODE, YMIR_SET);

    ymir_read_view_release (elemvisc);
    ymir_read_view_release (elemsvisc);
    ymir_read_view_release (elemTItens);
  }

  sc_dmatrix_destroy (elemin);
  sc_dmatrix_destroy (elemout);
  sc_dmatrix_destroy (elemvisc);
  sc_dmatrix_destroy (elemsvisc);
  sc_dmatrix_destroy (elemTItens);
}

/* compute traction. It is an alternative approach that takes advantage of an existing
   subroutine ymir_velocity_strain_rate and directly compute traction on each node*/
slabs_traction (ymir_cvec_t * vel, ymir_dvec_t *traction,
                  double *n_dir, ymir_dvec_t *visc)
{
  ymir_mesh_t        *mesh = vel->mesh;
  const int           N = ymir_n (mesh->cnodes->N);
  const int           Np = (N + 1)*(N+1)*(N+1);
  ymir_locidx_t       K = mesh->cnodes->K;
  ymir_dvec_t         *tau_tensor = ymir_dvec_new (mesh, 6,
                                                     YMIR_GAUSS_NODE);
  ymir_locidx_t       elid;
  sc_dmatrix_t       *an = sc_dmatrix_new (0, 0);
  sc_dmatrix_t       *bn = sc_dmatrix_new (0, 0);
  int                 i, j;

  ymir_velocity_strain_rate (vel, tau_tensor, 0);
  ymir_dvec_multiply_in1 (visc, tau_tensor);
  for (elid = 0; elid < K; elid++) {
    for (j = 0; j < Np; j++) {
      double              val = 0.;
      ymir_dvec_get_node (tau_tensor, an, elid, j, YMIR_READ);
      ymir_dvec_get_node (traction,   bn, elid, j, YMIR_WRITE);
      bn->e[0][0] = an->e[0][0] * n_dir[0] + an->e[0][1] * n_dir[1] + an->e[0][2] * n_dir[2];
      bn->e[0][1] = an->e[0][1] * n_dir[0] + an->e[0][3] * n_dir[1] + an->e[0][4] * n_dir[2];
      bn->e[0][2] = an->e[0][2] * n_dir[0] + an->e[0][4] * n_dir[1] + an->e[0][5] * n_dir[2];
      ymir_read_view_release (an);
      ymir_dvec_set_node (traction, bn, elid, j, YMIR_SET);
    }
  }
  sc_dmatrix_destroy (an);
  sc_dmatrix_destroy (bn);
}

/* compute traction as well as normal/shear stress at each element*/
void
slabs_normal_stress_elem (sc_dmatrix_t * in, sc_dmatrix_t * out1, sc_dmatrix_t * out2,
                            sc_dmatrix_t * out3, double * n_dir,
                  sc_dmatrix_t * visc_mat, ymir_velocity_elem_t * vel_elem,
                  double *_sc_restrict rxd, double *_sc_restrict sxd,
                  double *_sc_restrict txd, double *_sc_restrict ryd,
                  double *_sc_restrict syd, double *_sc_restrict tyd,
                  double *_sc_restrict rzd, double *_sc_restrict szd,
                  double *_sc_restrict tzd, sc_dmatrix_t * drst, sc_dmatrix_t * brst)
{
  const int           N = ymir_n (vel_elem->N);
  int                 gp, i, j, k, l;
  double              tempmatd[9];
  double              tempvec[3];
  double              temp;
  double             *_sc_restrict viscd = visc_mat->e[0];
  sc_dmatrix_t       *tempvec1 = vel_elem->tempvec1;
  sc_dmatrix_t       *tempvec2 = vel_elem->tempvec2;
  sc_dmatrix_t       *tempvec3 = vel_elem->tempvec3;
  sc_dmatrix_t       *temptens = vel_elem->temptens1;

  ymir_derivative_elem_grad (N, 3, drst->e[0], brst->e[0],
                             rxd, sxd, txd, ryd, syd, tyd,
                             rzd, szd, tzd, in, temptens, tempvec1, tempvec2,
                             tempvec3, 0);

  /* create stress from duvw/dxyz * viscosity */
  for (gp = 0, k = 0; k < N + 1; k++) {
    for (j = 0; j < N + 1; j++) {
      for (i = 0; i < N + 1; i++, gp++) {
        double               *_sc_restrict temptensd = temptens->e[0] + 9 * gp;
        double               *_sc_restrict normal = out1->e[0] + gp;
        double               *_sc_restrict shear  = out2->e[0] + gp;
        double               *_sc_restrict trac   = out3->e[0] + 3 * gp;

        tempmatd[0] = temptensd[0];
        tempmatd[1] = (temptensd[1] + temptensd[3]) * (1. / 2.);
        tempmatd[2] = (temptensd[2] + temptensd[6]) * (1. / 2.);
        tempmatd[3] = tempmatd[1];
        tempmatd[4] = temptensd[4];
        tempmatd[5] = (temptensd[5] + temptensd[7]) * (1. / 2.);
        tempmatd[6] = tempmatd[2];
        tempmatd[7] = tempmatd[5];
        tempmatd[8] = temptensd[8];
        *(trac++) = tempvec[0] = ( tempmatd[0] * n_dir[0]
                                 + tempmatd[1] * n_dir[1]
                                 + tempmatd[2] * n_dir[2]) * 2.0 * viscd[gp];
        *(trac++) = tempvec[1] = ( tempmatd[3] * n_dir[0]
                                 + tempmatd[4] * n_dir[1]
                                 + tempmatd[5] * n_dir[2]) * 2.0 * viscd[gp];
        *(trac++) = tempvec[2] = ( tempmatd[6] * n_dir[0]
                                 + tempmatd[7] * n_dir[1]
                                 + tempmatd[8] * n_dir[2]) * 2.0 * viscd[gp];
        *normal = temp = tempvec[0] * n_dir[0] + tempvec[1] * n_dir[1] + tempvec[2] * n_dir[2];
        *shear  = sqrt( tempvec[0] * tempvec[0] + tempvec[1] * tempvec[1] + tempvec[2] * tempvec[2]
                        - temp * temp);
      }
    }
  }
}

/* compute traction as well as normal and shear stress*/
void
slabs_normal_stress (ymir_cvec_t * vel, ymir_dvec_t * n_tau, ymir_dvec_t * s_tau,
                        ymir_dvec_t * traction, double * n_dir,
                       ymir_dvec_t * visc, ymir_velocity_elem_t * vel_elem)
{
  ymir_mesh_t         *mesh = vel->mesh;
  ymir_locidx_t       elid;
  const int           N  = ymir_n (mesh->cnodes->N);
  const int           Np = ymir_np (mesh->cnodes->N);
  const int           K  = mesh->cnodes->K;
  sc_dmatrix_t       *elemin = sc_dmatrix_new (1, 3 * Np);
  sc_dmatrix_t       *elemout1 = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemout2 = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemout3 = sc_dmatrix_new (1, 3 * Np);
  sc_dmatrix_t       *elemvisc = sc_dmatrix_new (1, Np);

  ymir_dvec_set_zero (traction);
  ymir_dvec_set_zero (n_tau);
  ymir_dvec_set_zero (s_tau);

  for (elid = 0; elid < K; elid++)  {
    double             *_sc_restrict rxd = mesh->drdx->e[elid];
    double             *_sc_restrict sxd = mesh->dsdx->e[elid];
    double             *_sc_restrict txd = mesh->dtdx->e[elid];
    double             *_sc_restrict ryd = mesh->drdy->e[elid];
    double             *_sc_restrict syd = mesh->dsdy->e[elid];
    double             *_sc_restrict tyd = mesh->dtdy->e[elid];
    double             *_sc_restrict rzd = mesh->drdz->e[elid];
    double             *_sc_restrict szd = mesh->dsdz->e[elid];
    double             *_sc_restrict tzd = mesh->dtdz->e[elid];
    double             *_sc_restrict Jdetd = mesh->Jdet->e[elid];

    ymir_cvec_get_elem_interp (vel, elemin, YMIR_STRIDE_NODE, elid,
                               YMIR_GLL_NODE, YMIR_COPY);
    ymir_dvec_get_elem_interp (visc, elemvisc, YMIR_STRIDE_COMP, elid,
                              YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (n_tau, elemout1, YMIR_STRIDE_COMP, elid,
                             YMIR_GAUSS_NODE, YMIR_WRITE);
    ymir_dvec_get_elem_interp (s_tau, elemout2, YMIR_STRIDE_COMP, elid,
                             YMIR_GAUSS_NODE, YMIR_WRITE);

    slabs_normal_stress_elem (elemin, elemout1, elemout2, elemout3, n_dir, elemvisc,
                            vel_elem, rxd, sxd, txd, ryd,
                            syd, tyd, rzd, szd, tzd, mesh->drst, mesh->brst);

    ymir_dvec_set_elem_interp (n_tau, elemout1, YMIR_STRIDE_COMP, elid,
                             YMIR_GAUSS_NODE, YMIR_SET);
    ymir_dvec_set_elem_interp (s_tau, elemout2, YMIR_STRIDE_COMP, elid,
                             YMIR_GAUSS_NODE, YMIR_SET);
    ymir_dvec_set_elem_interp (traction, elemout3, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);

    ymir_read_view_release (elemvisc);
  }

  sc_dmatrix_destroy (elemin);
  sc_dmatrix_destroy (elemout1);
  sc_dmatrix_destroy (elemout2);
  sc_dmatrix_destroy (elemout3);
  sc_dmatrix_destroy (elemvisc);
}

static void
slabs_physics_normal_boundary_stress_fn (double *stress_norm,
                                         double x, double y, double z,
                                         double nx, double ny, double nz,
                                         ymir_topidx_t face,
                                         ymir_locidx_t node_id,
                                         void *data)
{
  ymir_vec_t         *vec_bndr = (ymir_vec_t *) data;
  double             *v = ymir_cvec_index (vec_bndr, node_id, 0);

  RHEA_ASSERT (vec_bndr->ncfields == 3);

  /* compute inner product with boundary outer normal vector */
  *stress_norm = nx * v[0] + ny * v[1] + nz * v[2];
}

/* normal stress on the surface */
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
  RHEA_ASSERT (stress_bndr_norm->ncfields == 1);
  RHEA_ASSERT (ymir_stokes_vec_is_stokes_vec (up));
  YMIR_ASSERT_IS_CVEC (rhs_u_point);
  RHEA_ASSERT (ymir_vec_is_not_dirty (up));
  RHEA_ASSERT (ymir_vec_is_not_dirty (rhs_u_point));

  /* construct the right-hand side */
  ymir_stokes_pc_construct_rhs (rhs, rhs_u_point, NULL, NULL,
                                1 /* incompressible */, stokes_op, 0);
  RHEA_ASSERT (sc_dmatrix_is_valid (rhs->dataown));
  RHEA_ASSERT (sc_dmatrix_is_valid (rhs->coff));
  RHEA_ASSERT (ymir_vec_is_not_dirty (rhs));

  /* turn off boundary constraints */
  stress_op->skip_dir = 1;

  /* compute (unconstrained) residual
   *   r_mom  = a * u + b^t * p - f
   *   r_mass = b * u
   */
  ymir_stokes_pc_apply_stokes_op (up, residual_up, stokes_op, 0, 0);
  ymir_vec_add (-1.0, rhs, residual_up);
  RHEA_ASSERT (sc_dmatrix_is_valid (residual_up->dataown));
  RHEA_ASSERT (sc_dmatrix_is_valid (residual_up->coff));
  ymir_vec_destroy (rhs);

  /* restore boundary constraints */
  stress_op->skip_dir = skip_dir;

  /* get the velocity component of the residual */
  ymir_stokes_vec_get_velocity (residual_up, residual_u, stokes_op->press_elem);

  /* interpolate residual onto boundary */
  ymir_interp_vec (residual_u, residual_bndr);
  RHEA_ASSERT (sc_dmatrix_is_valid (residual_bndr->dataown));
  RHEA_ASSERT (sc_dmatrix_is_valid (residual_bndr->coff));
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

/**************************************************************
 * Non-zero Dirichlet boundary conditions
***************************************************************/

/* Dirichlet all and tangential on one side: A test case 'wallslide'*/
static ymir_dir_code_t
slabs_set_vel_dir_freeslip (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
    return YMIR_VEL_DIRICHLET_NORM;
}

/*
 * Dirichlet all on one side of the domain: SIDE3 (y=0).
 */
static ymir_dir_code_t
slabs_set_vel_dir_inoutflow (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
  if (face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE3) {
    return YMIR_VEL_DIRICHLET_ALL;
  }
  else if (face == RHEA_DOMAIN_BOUNDARY_FACE_BASE) {
     return YMIR_VEL_DIRICHLET_NORM;
  }
  else {
    return YMIR_VEL_DIRICHLET_NORM;
  }
}

/**
 * In-out flow sine velocity on one side of the domain.
 */
void
slabs_set_rhs_vel_nonzero_dir_inoutflow_sin (
    double *vel, double x, double y, double z,
    ymir_locidx_t nid, void *data)
{
  slabs_options_t  *slabs_options = data;
  const double     flow_scale = slabs_options->slabs_velbc_options->flow_scale;

  if (fabs (y) < SC_1000_EPS) {
    vel[0] = 0.0;
    vel[1] = -flow_scale * sin (2.0 * M_PI * z);
    vel[2] = 0.0;
  }
  else {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }
}

/**
 * In-out flow tanh velocity on one side of the domain.
 */
void
slabs_set_rhs_vel_nonzero_dir_inoutflow_tanh (
    double *vel, double x, double y, double z,
    ymir_locidx_t nid, void *data)
{
  slabs_options_t  *slabs_options = data;
  slabs_velbc_options_t *velbc_options = slabs_options->slabs_velbc_options;
  const double        flow_scale = velbc_options->flow_scale;
  const double        zU = velbc_options->vel_dir_bc_upper;
  const double        zL = velbc_options->vel_dir_bc_lower;
  const double        a = zU-zL, b = 1.0-a, c = 0.5*(zU+zL);
  const double        shape = 2.0 * M_PI, scaling = 0.5*(b+a)*flow_scale, shift = 0.5*(b-a)*flow_scale;
  double txL = shape*(z-zL),txU = shape*(z-zU);

  if (fabs (y) < SC_1000_EPS) {
    vel[0] = 0.0;
    vel[2] = 0.0;
    if (y < SC_1000_EPS)  {
      if (z<=c)
        vel[1] = shift + scaling *
               ( (exp (txL) - exp (-txL)) /
               (exp (txL) + exp (-txL)) );
      else
        vel[1] = shift - scaling *
               ( (exp (txU) - exp (-txU)) /
               (exp (txU) + exp (-txU)) );
    }
  }
  else if ((2.0 - y) < SC_1000_EPS) {
    vel[0] = 0.0;
    vel[2] = 0.0;
    if (z<=c)
      vel[1] = -shift - scaling *
             ( (exp (txL) - exp (-txL)) /
             (exp (txL) + exp (-txL)) );
    else
      vel[1] = -shift + scaling *
             ( (exp (txU) - exp (-txU)) /
             (exp (txU) + exp (-txU)) );
  }
  else {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }
}

void
slabs_vel_nonzero_dirichlet_compute ( ymir_vec_t * rhs_vel_nonzero_dirichlet,
                                      slabs_options_t * slabs_options,
                                      ymir_mesh_t * ymir_mesh)
{
  switch (slabs_options->slabs_velbc_options->vel_dir_bc) {
    case SLABS_VEL_DIR_BC_INOUTFLOW_SIN:
      rhea_domain_set_user_velocity_dirichlet_bc (
          slabs_set_vel_dir_inoutflow, NULL /* no data necessary */,
          0 /* TODO don't need this flag */);
      rhs_vel_nonzero_dirichlet = rhea_velocity_new (ymir_mesh);
      ymir_cvec_set_function (rhs_vel_nonzero_dirichlet,
                              slabs_set_rhs_vel_nonzero_dir_inoutflow_sin,
                              slabs_options);
      break;

    case SLABS_VEL_DIR_BC_INOUTFLOW_TANH:
      rhea_domain_set_user_velocity_dirichlet_bc (
          slabs_set_vel_dir_inoutflow, NULL /* no data necessary */,
          0 /* TODO don't need this flag */);
      rhs_vel_nonzero_dirichlet = rhea_velocity_new (ymir_mesh);
      ymir_cvec_set_function (rhs_vel_nonzero_dirichlet,
                              slabs_set_rhs_vel_nonzero_dir_inoutflow_tanh,
                              slabs_options);
      break;

   default: /* BC not set */
      RHEA_ABORT_NOT_REACHED ();
  }
}

/**
 * Sets up the mesh.
 **/
static void
slabs_setup_mesh (p4est_t **p4est,
                    ymir_mesh_t **ymir_mesh,
                    ymir_pressure_elem_t **press_elem,
                    MPI_Comm mpicomm,
                    rhea_domain_options_t *domain_options,
                    rhea_discretization_options_t *discr_options,
                    slabs_options_t *slabs_options)

{
  const char         *this_fn_name = "slabs_setup_mesh";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* create p4est */
  *p4est = rhea_discretization_p4est_new (mpicomm, discr_options,
                                          domain_options);

  /* set up boundary, store in `discr_options` */
  rhea_discretization_options_set_boundary (discr_options, *p4est,
                                            domain_options);

  /* create ymir mesh and pressure element */
  rhea_discretization_ymir_mesh_new_from_p4est (ymir_mesh, press_elem, *p4est,
                                                discr_options);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Sets up a linear Stokes problem.
 */
static void
slabs_setup_stokes (rhea_stokes_problem_t **stokes_problem,
                      ymir_mesh_t *ymir_mesh,
                      ymir_pressure_elem_t *press_elem,
                      rhea_domain_options_t *domain_options,
                      rhea_temperature_options_t *temp_options,
                      rhea_viscosity_options_t *visc_options,
                      slabs_options_t *slabs_options,
                      const char *vtk_write_input_path)
{
  const char         *this_fn_name = "slabs_setup_stokes";
  ymir_vec_t         *temperature, *weakzone;
  ymir_vec_t         *coeff_TI_svisc=NULL, *TI_rotate=NULL;
  ymir_vec_t         *rhs_vel, *rhs_vel_nonzero_dirichlet=NULL;
  void               *solver_options = NULL;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* compute temperature */
  temperature = rhea_temperature_new (ymir_mesh);
  slabs_temperature_compute (temperature, slabs_options); /* if non-specified, use:
                                                             rhea_temperature_compute
                                                             (temperature, temp_options); */
  /* compute weak zone */
  weakzone = rhea_viscosity_new (ymir_mesh);
  if (slabs_options->slabs_visc_options->viscosity_anisotropy
    == SLABS_VISC_TRANSVERSELY_ISOTROPY) {
    RHEA_GLOBAL_INFO ("Into TI\n");
    ymir_vec_set_value (weakzone, 1.0);
  }
  else {
    slabs_weakzone_compute (weakzone, slabs_options);
  }

  /* compute velocity right-hand side volume forcing */
  rhs_vel = rhea_velocity_new (ymir_mesh);
  rhea_temperature_compute_rhs_vel (rhs_vel, temperature, temp_options);

  /* set velocity boundary conditions & nonzero Dirichlet values */
  if (domain_options->velocity_bc_type == RHEA_DOMAIN_VELOCITY_BC_USER) {
    rhs_vel_nonzero_dirichlet = rhea_velocity_new (ymir_mesh);
    slabs_vel_nonzero_dirichlet_compute (rhs_vel_nonzero_dirichlet,
                                         slabs_options, ymir_mesh);
  }

  /* create Stokes problem */
  *stokes_problem = rhea_stokes_problem_new (
      temperature, weakzone, rhs_vel, rhs_vel_nonzero_dirichlet,
      ymir_mesh, press_elem, domain_options, visc_options, solver_options);

  /* add the anisotropic viscosity to the viscous stress operator */
  if (slabs_options->slabs_visc_options->viscosity_anisotropy
      == SLABS_VISC_TRANSVERSELY_ISOTROPY) {
    coeff_TI_svisc = rhea_viscosity_new (ymir_mesh);
    TI_rotate = rhea_viscosity_new (ymir_mesh);
    slabs_stokes_problem_setup_TI (ymir_mesh, *stokes_problem, slabs_options,
                                   coeff_TI_svisc, TI_rotate);
  }

  /* write vtk of problem input */
  if (vtk_write_input_path != NULL) {
    slabs_write_input (ymir_mesh, *stokes_problem, temp_options,
                       temperature, weakzone, coeff_TI_svisc, TI_rotate,
                       vtk_write_input_path);
  }

  /* set up Stokes solver */
  rhea_stokes_problem_setup_solver (*stokes_problem);

  /* destroy */
  if ( slabs_options->slabs_visc_options->viscosity_anisotropy
      == SLABS_VISC_TRANSVERSELY_ISOTROPY) {
    rhea_viscosity_destroy (TI_rotate);
  }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Cleans up Stokes problem and mesh.
 */
static void
slabs_setup_clear_all (rhea_stokes_problem_t *stokes_problem,
                         p4est_t *p4est,
                         ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         rhea_discretization_options_t *discr_options)
{
  const char         *this_fn_name = "slabs_setup_clear_all";
  ymir_vec_t         *temperature, *weakzone;
  ymir_vec_t         *visc_TI_svisc;
  ymir_vec_t         *rhs_vel, *rhs_vel_nonzero_dirichlet;

  RHEA_GLOBAL_PRODUCTIONF ("into %s\n", this_fn_name);

  /* get vectors */
  temperature = rhea_stokes_problem_get_temperature (stokes_problem);
  weakzone = rhea_stokes_problem_get_weakzone (stokes_problem);
  rhs_vel = rhea_stokes_problem_get_rhs_vel (stokes_problem);
  rhs_vel_nonzero_dirichlet =
    rhea_stokes_problem_get_rhs_vel_nonzero_dirichlet (stokes_problem);

  /* destroy anisotropic viscosity */
  {
    ymir_stokes_op_t    *stokes_op;
    ymir_vec_t          *visc_TI_svisc;

    stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
    visc_TI_svisc = stokes_op->stress_op->coeff_TI_svisc;

    if (visc_TI_svisc != NULL) {
      rhea_viscosity_destroy (visc_TI_svisc);
    }
  }

  /* destroy Stokes problem */
  rhea_stokes_problem_destroy (stokes_problem);

  /* destroy vectors */
  if (temperature != NULL) {
    rhea_temperature_destroy (temperature);
  }
  if (weakzone != NULL) {
    rhea_weakzone_destroy (weakzone);
  }
  if (rhs_vel != NULL) {
    rhea_velocity_destroy (rhs_vel);
  }
  if (rhs_vel_nonzero_dirichlet != NULL) {
    rhea_velocity_destroy (rhs_vel_nonzero_dirichlet);
  }

  /* destroy mesh */
  rhea_discretization_ymir_mesh_destroy (ymir_mesh, press_elem);
  rhea_discretization_p4est_destroy (p4est);

  /* destroy (some) options */
  rhea_discretization_options_clear (discr_options);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * runs stokes solver.
 */
static void
slabs_run_solver (ymir_vec_t *sol_vel_press,
                    ymir_mesh_t *ymir_mesh,
                    ymir_pressure_elem_t *press_elem,
                    rhea_stokes_problem_t *stokes_problem,
                    const int iter_max, const double rel_tol)
{
  const char         *this_fn_name = "slabs_run_solver";
  ymir_vec_t         *rhs_vel_nonzero_dirichlet;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* run solver */
  rhea_stokes_problem_solve (sol_vel_press, iter_max, rel_tol, stokes_problem);

  /* add nonzero dirichlet values of the velocity to the solution */
  rhs_vel_nonzero_dirichlet =
    rhea_stokes_problem_get_rhs_vel_nonzero_dirichlet (stokes_problem);
  if (rhs_vel_nonzero_dirichlet != NULL) {
    ymir_vec_t         *sol_vel = rhea_velocity_new (ymir_mesh);

    ymir_stokes_vec_get_velocity (sol_vel_press, sol_vel, press_elem);
    ymir_vec_add (1.0, rhs_vel_nonzero_dirichlet, sol_vel);
    ymir_stokes_vec_set_velocity (sol_vel, sol_vel_press, press_elem);
    rhea_velocity_destroy (sol_vel);
  }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Runs the program.
 */
int
main (int argc, char **argv)
{
  const char         *this_fn_name = "slabweakzone:main";
  /* MPI */
  MPI_Comm            mpicomm = MPI_COMM_WORLD;
  int                 mpisize, mpirank, ompsize;
  int                 mpiret;
  /* options */
  ymir_options_t     *opt;
  rhea_domain_options_t         domain_options;
  rhea_temperature_options_t    temp_options;
  rhea_viscosity_options_t      visc_options;
  rhea_discretization_options_t discr_options;
  rhea_newton_options_t         newton_options;

  /* slabs options */
  slabs_domain_options_t     slabs_domain_options;
  slabs_temp_options_t     slabs_temp_options;
  slabs_visc_options_t     slabs_visc_options;
  slabs_weak_options_t     slabs_weak_options;
  slabs_velbc_options_t    slabs_velbc_options;
  slabs_options_t          slabs_options;

  /* temperature */
  int                 vel_dir_bc;
  double              flow_scale;
  double              velocity_bc_upper;
  double              velocity_bc_lower;
  int                 viscosity_anisotropy;

  /* options local to this function */
  int                 production_run;
  int                 solver_iter_max;
  double              solver_rel_tol;
  char               *vtk_write_input_path;
  char               *vtk_write_solution_path;
  char               *vtk_write_stress_path;
  char               *vtk_write_traction_path;
  /* mesh */
  p4est_t            *p4est;
  ymir_mesh_t        *ymir_mesh;
  ymir_pressure_elem_t  *press_elem;
  /* Stokes */
  rhea_stokes_problem_t *stokes_problem;
  ymir_vec_t         *sol_vel_press;

  /*
   * Initialize Libraries
   */

  /* initialize rhea and sub-packages */
  rhea_initialize (argc, argv, mpicomm);

  /* get parallel environment */
  mpiret = MPI_Comm_size (mpicomm, &mpisize); YMIR_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &mpirank); YMIR_CHECK_MPI (mpiret);

#ifdef RHEA_ENABLE_OPENMP
  ompsize = omp_get_max_threads ();
#else
  ompsize = 1;
#endif

  /*
   * Define & Parse Options
   */

  opt = ymir_options_global_new (argv[0] /* program path */);

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  /* basic options */
  YMIR_OPTIONS_CALLBACK, "help", 'h', 0 /* no callback fn args */,
    ymir_options_print_usage_and_exit_fn, NULL /* no arg usage */,
    "Print usage and exit",
  YMIR_OPTIONS_INIFILE, "options-file", 'f',
    ".ini file with option values",

  /* velocity Dirichlet BC's */
  YMIR_OPTIONS_I, "velocity-dirichlet-bc", '\0',
    &vel_dir_bc, SLABS_VEL_DIR_BC_INOUTFLOW_SIN,
    "Velocity Dirichlet boundary condition",
  YMIR_OPTIONS_D, "flow-scaling", '\0',
    &flow_scale, 1.0,
    "scaling of velocity BC.",
  YMIR_OPTIONS_D, "velocity-bc-upper", '\0',
    &velocity_bc_upper, 1.0,
    "location of velocity BC: upper bound",
  YMIR_OPTIONS_D, "velocity-bc-lower", '\0',
    &velocity_bc_lower, 0.0,
    "location of velocity BC: lower bound",

  /* temperature */
  YMIR_OPTIONS_D, "temp-background-plate-age", '\0',
    &temp_back_plate_age, SLABS_DEFAULT_TEMP_BACKGROUND_PLATE_AGE,
    "Bachground temperature descibed by plate age [yr]",
  YMIR_OPTIONS_D, "temp-2plates-trench-longitude", '\0',
    &temp_2pl_trench_lon, SLABS_DEFAULT_TEMP_2PL_TRENCH_LONGITUDE,
    "2plates temp: Longitude of trench in interval (-pi/8, pi/8)",
  YMIR_OPTIONS_D, "temp-2plates-dip-angle", '\0',
    &temp_2pl_dip_angle, SLABS_DEFAULT_TEMP_2PL_DIP_ANGLE,
    "2plates temp: Dip angle of subducting plate (in degrees < 0)",
  YMIR_OPTIONS_D, "temp-2plates-subd-depth", '\0',
    &temp_2pl_subd_depth, SLABS_DEFAULT_TEMP_2PL_SUBD_DEPTH,
    "2plates temp: Maximal depth of subducting plate inside of mantle [m]",
  YMIR_OPTIONS_D, "temp-2plates-subd-width", '\0',
    &temp_2pl_subd_width, SLABS_DEFAULT_TEMP_2PL_SUBD_WIDTH,
    "2plates temp: Maximal width of subducting zone inside of mantle [m]",
  YMIR_OPTIONS_D, "temp-2plates-subd-edge-width", '\0',
    &temp_2pl_subd_edge_width, SLABS_DEFAULT_TEMP_2PL_SUBD_EDGE_WIDTH,
    "2plates temp: Width for subducting plate's top edge [m]",
  YMIR_OPTIONS_D, "temp-2plates-subd-edge-smoothwidth", '\0',
    &temp_2pl_subd_edge_smoothwidth,
    SLABS_DEFAULT_TEMP_2PL_SUBD_EDGE_SMOOTHWIDTH,
    "2plates weak zone: Width of smoothing of subd. plate's top edge [m]",
  YMIR_OPTIONS_D, "temp-2plates-subd-plate-velocity", '\0',
    &temp_2pl_subd_plate_vel, SLABS_DEFAULT_TEMP_2PL_SUBD_PLATE_VELOCITY,
    "2plates temp: Velocity of subducting plate [m/y]",
  YMIR_OPTIONS_D, "temp-2plates-subd-plate-initial-age", '\0',
    &temp_2pl_subd_plate_init_age,
    SLABS_DEFAULT_TEMP_2PL_SUBD_PLATE_INITIAL_AGE,
    "2plates temp: Age of subducting plate at left boundary [y]",
  YMIR_OPTIONS_D, "temp-2plates-over-plate-age", '\0',
    &temp_2pl_over_plate_age, SLABS_DEFAULT_TEMP_2PL_OVER_PLATE_AGE,
    "2plates temp: Age of overriding plate [y]",

  /* viscosity */
  YMIR_OPTIONS_I, "viscosity-anisotropy",'\0',
    &(viscosity_anisotropy), SLABS_VISC_ISOTROPY,
    "0: isotropy, 1: transversely isotropy",

  /* weakzone */
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-longitude", '\0',
    &weakzone_2pl_subdu_lon, SLABS_WEAKZONE_2PLATES_SUBDU_LONGITUDE,
    "2plates weak zone: Longitude in interval (-pi/8, pi/8), "
    "where weak zone begins",
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-dip-angle", '\0',
    &weakzone_2pl_subdu_dip_angle, SLABS_WEAKZONE_2PLATES_SUBDU_DIP_ANGLE,
    "2plates weak zone: Dip angle of weak zone (in degrees > 0)",
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-depth", '\0',
    &weakzone_2pl_subdu_depth, SLABS_WEAKZONE_2PLATES_SUBDU_DEPTH,
    "2plates weak zone: Depth of weak zone [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-width", '\0',
    &weakzone_2pl_subdu_width, SLABS_WEAKZONE_2PLATES_SUBDU_WIDTH,
    "2plates weak zone: Width of weak zone [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-thickness", '\0',
    &weakzone_2pl_subdu_thickness, SLABS_WEAKZONE_2PLATES_SUBDU_THICKNESS,
    "2plates weak zone: Width at center of weak zone [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-thickness-const", '\0',
    &weakzone_2pl_subdu_thickness_const,
    SLABS_WEAKZONE_2PLATES_SUBDU_THICKNESS_CONST,
    "2plates weak zone: Width of smoothing of edges of weak zone [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-weak-factor", '\0',
    &weakzone_2pl_subdu_weak_factor, SLABS_WEAKZONE_2PLATES_SUBDU_WEAK_FACTOR,
    "2plates weak zone: Value of weak zone factor",
  YMIR_OPTIONS_D, "weakzone-2plates-ridge-depth", '\0',
    &weakzone_2pl_ridge_depth, SLABS_WEAKZONE_2PLATES_RIDGE_DEPTH,
    "2plates weak zone: Depth of weak zone in left corner of domain [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-ridge-width", '\0',
    &weakzone_2pl_ridge_width, SLABS_WEAKZONE_2PLATES_RIDGE_WIDTH,
    "2plates weak zone: Width of weak zone in left corner of domain [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-ridge-smoothwidth", '\0',
    &weakzone_2pl_ridge_smoothwidth, SLABS_WEAKZONE_2PLATES_RIDGE_SMOOTHWIDTH,
    "2plates weak zone: Smoothing width of edges of weak zone in corner [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-ridge-weak-factor", '\0',
    &weakzone_2pl_ridge_weak_factor, SLABS_WEAKZONE_2PLATES_RIDGE_WEAK_FACTOR,
    "2plates weak zone: Value of weak zone factor for weak zone in corner",

  /* solver options */
  YMIR_OPTIONS_I, "solver-iter-max", '\0',
    &solver_iter_max, 100,
    "Maximum number of iterations for Stokes solver",
  YMIR_OPTIONS_D, "solver-rel-tol", '\0',
    &solver_rel_tol, 1.0e-6,
    "Relative tolerance for Stokes solver",

  /* performance & monitoring options */
  YMIR_OPTIONS_B, "production-run", '\0',
    &(production_run), 0,
    "Execute as a production run (to reduce some overhead and checks)",

  /* vtk output options */
  YMIR_OPTIONS_S, "vtk-write-input-path", '\0',
    &(vtk_write_input_path), NULL,
    "File path for vtk files for the input of the Stokes problem",
  YMIR_OPTIONS_S, "vtk-write-solution-path", '\0',
    &(vtk_write_solution_path), NULL,
    "File path for vtk files for the solution of the Stokes problem",
 YMIR_OPTIONS_S, "vtk-write-stress-path", '\0',
    &(vtk_write_stress_path), NULL,
    "File path for vtk files for stress results",
  YMIR_OPTIONS_S, "vtk-write-traction-path", '\0',
    &(vtk_write_traction_path), NULL,
    "File path for vtk files for traction results",


  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add sub-options */
  rhea_add_options_all (opt);
  ymir_options_add_suboptions_solver_stokes (opt);

  /* parse options */
  {
    int                 optret;

    optret = ymir_options_parse (SC_LP_INFO, opt, argc, argv);
    if (optret < 0) { /* if parsing was not successful */
      ymir_options_print_usage (SC_LP_INFO, opt, NULL /* args usage */);
      RHEA_GLOBAL_INFO ("Option parsing failed\n");
      exit (0);
    }
  }

  /*
   * Process Slabs Options
   */
  /* temperature */
  slabs_temp_options.temp_background_plate_age = temp_back_plate_age;
  slabs_temp_options.temp_2plates_trench_longitude = temp_2pl_trench_lon;
  slabs_temp_options.temp_2plates_dip_angle = temp_2pl_dip_angle;
  slabs_temp_options.temp_2plates_subd_depth = temp_2pl_subd_depth;
  slabs_temp_options.temp_2plates_subd_width = temp_2pl_subd_width;
  slabs_temp_options.temp_2plates_subd_edge_width = temp_2pl_subd_edge_width;
  slabs_temp_options.temp_2plates_subd_edge_smoothwidth =
    temp_2pl_subd_edge_smoothwidth;
  slabs_temp_options.temp_2plates_subd_plate_velocity = temp_2pl_subd_plate_vel;
  slabs_temp_options.temp_2plates_subd_plate_initial_age =
    temp_2pl_subd_plate_init_age;
  slabs_temp_options.temp_2plates_over_plate_age = temp_2pl_over_plate_age;


  /* viscosity */
  slabs_visc_options.viscosity_anisotropy = (slabs_viscosity_anisotropy_t) viscosity_anisotropy;

  /* weak zone */
  slabs_weak_options.weakzone_2plates_subdu_longitude =
    weakzone_2pl_subdu_lon;
  slabs_weak_options.weakzone_2plates_subdu_dip_angle =
    weakzone_2pl_subdu_dip_angle;
  slabs_weak_options.weakzone_2plates_subdu_depth =
    weakzone_2pl_subdu_depth;
  slabs_weak_options.weakzone_2plates_subdu_width =
    weakzone_2pl_subdu_width;
  slabs_weak_options.weakzone_2plates_subdu_thickness =
    weakzone_2pl_subdu_thickness;
  slabs_weak_options.weakzone_2plates_subdu_thickness_const =
    weakzone_2pl_subdu_thickness_const;
  slabs_weak_options.weakzone_2plates_subdu_weak_factor =
    weakzone_2pl_subdu_weak_factor;
  slabs_weak_options.weakzone_2plates_ridge_depth =
    weakzone_2pl_ridge_depth;
  slabs_weak_options.weakzone_2plates_ridge_width =
    weakzone_2pl_ridge_width;
  slabs_weak_options.weakzone_2plates_ridge_smoothwidth =
    weakzone_2pl_ridge_smoothwidth;
  slabs_weak_options.weakzone_2plates_ridge_weak_factor =
    weakzone_2pl_ridge_weak_factor;

  /* velocity B.C. condition */
  slabs_velbc_options.vel_dir_bc = (slabs_vel_dir_bc_t) vel_dir_bc;
  slabs_velbc_options.flow_scale = flow_scale;
  slabs_velbc_options.vel_dir_bc_upper = velocity_bc_upper;
  slabs_velbc_options.vel_dir_bc_lower = velocity_bc_lower;

  /* assign slabs_options */
  slabs_options.slabs_temp_options = &slabs_temp_options;
  slabs_options.slabs_visc_options = &slabs_visc_options;
  slabs_options.slabs_weak_options = &slabs_weak_options;
  slabs_options.slabs_velbc_options = &slabs_velbc_options;

  /*
   * Initialize Main Program
   */

  RHEA_GLOBAL_PRODUCTIONF (
      "Into %s (production %i)\n", this_fn_name, production_run);
  RHEA_GLOBAL_PRODUCTIONF (
      "Parallel environment: MPI size %i, OpenMP size %i\n", mpisize, ompsize);
  ymir_set_up (argc, argv, mpicomm, production_run);

  /* print & process options */
  ymir_options_print_summary (SC_LP_INFO, opt);
  rhea_process_options_all (&domain_options, &temp_options,
                            &visc_options, &discr_options,
                            &newton_options);

  /* copy rhea domain options into local example domain options */
  slabs_domain_options.x_min = domain_options.x_min;
  slabs_domain_options.x_max = domain_options.x_max;
  slabs_domain_options.y_min = domain_options.y_min;
  slabs_domain_options.y_max = domain_options.y_max;
  slabs_domain_options.z_min = domain_options.z_min;
  slabs_domain_options.z_max = domain_options.z_max;
  slabs_domain_options.lon_min = domain_options.lon_min;
  slabs_domain_options.lon_max = domain_options.lon_max;
  slabs_domain_options.radius_min = domain_options.radius_min;
  slabs_domain_options.radius_max = domain_options.radius_max;
  slabs_options.slabs_domain_options = &slabs_domain_options;

  /*
   * Setup Mesh
   */

  slabs_setup_mesh (&p4est, &ymir_mesh, &press_elem, mpicomm,
                      &domain_options, &discr_options, &slabs_options);

  /*
   * Setup Stokes Problem
   */

  slabs_setup_stokes (&stokes_problem, ymir_mesh, press_elem,
                        &domain_options, &temp_options, &visc_options,
                        &slabs_options, vtk_write_input_path);

  /*
   * Solve Stokes Problem
   */

  /* initialize solution vector */
  sol_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);

  /* run solver */
  slabs_run_solver (sol_vel_press, ymir_mesh, press_elem, stokes_problem,
                      solver_iter_max, solver_rel_tol);

  /* write vtk of solution */
  if (vtk_write_solution_path != NULL) {
    ymir_vec_t         *velocity = rhea_velocity_new (ymir_mesh);
    ymir_vec_t         *pressure = rhea_pressure_new (ymir_mesh, press_elem);
    ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);

    ymir_stokes_vec_get_components (sol_vel_press, velocity, pressure,
                                    press_elem);
    rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

    rhea_vtk_write_solution (vtk_write_solution_path, velocity, pressure,
                             viscosity);

    rhea_pressure_destroy (pressure);
    rhea_viscosity_destroy (viscosity);
    rhea_velocity_destroy (velocity);
  }

  /* compute and output second invariant strain_rate, stress, and surface normal stress  */
  if (vtk_write_stress_path != NULL)  {
    ymir_vec_t         *velocity = rhea_velocity_new (ymir_mesh);
    ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);
    ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (
                                                    ymir_mesh->ma->N, ymir_mesh->ma->ompsize);
    ymir_vec_t            *edotII = ymir_dvec_new (ymir_mesh, 1,
                                                        YMIR_GAUSS_NODE);
    ymir_vec_t            *tauII = ymir_dvec_new (ymir_mesh, 1,
                                                        YMIR_GAUSS_NODE);

    ymir_vec_t            *surf_normal_stress = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);

    ymir_stokes_vec_get_velocity (sol_vel_press, velocity,
                                    press_elem);
    rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

    /* compute 2nd invariant of the strain rate */
    ymir_second_invariant_vec (velocity, edotII, vel_elem);
    ymir_vec_sqrt (edotII, edotII);

    /* compute 2nd invariant of deviatoric stress tau = 2* (2nd invariant of strain_rate * viscosity )
      and its projection on the surface */
    if (slabs_visc_options.viscosity_anisotropy == SLABS_VISC_TRANSVERSELY_ISOTROPY)  {
      ymir_stokes_op_t      *stokes_op;
      ymir_stress_op_t      *stress_op;
      ymir_vec_t            *shear_visc = rhea_viscosity_new (ymir_mesh);
      ymir_vec_t            *TI_tensor = ymir_dvec_new (ymir_mesh, 9,
                                                        YMIR_GAUSS_NODE);

      /* get the viscous stress operator */
      stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
      stress_op = stokes_op->stress_op;

      /* copy shear viscosity */
      slabs_stress_op_copy_shear_visc (shear_visc, stress_op);
      slabs_stress_op_copy_TI_tensor (TI_tensor, stress_op);

      slabs_2inv_stress_TI (velocity, tauII,
                              viscosity, shear_visc, TI_tensor, vel_elem);

      rhea_viscosity_destroy (shear_visc);
      ymir_vec_destroy (TI_tensor);
    }
    else  {
      ymir_vec_copy (edotII, tauII)
      ymir_vec_multiply_in (viscosity, tauII);
      ymir_vec_scale (2.0, tauII);
    }

    /* compute surface normal stress sigma */
    slabs_physics_compute_normal_boundary_stress (
                   surf_normal_stress, sol_vel_press,
                   rhea_stokes_problem_get_rhs_vel (stokes_problem),
                   rhea_stokes_problem_get_stokes_op (stokes_problem));

    {
      char            path[BUFSIZ];

      snprintf (path, BUFSIZ, "%s", vtk_write_stress_path);
      ymir_vtk_write (ymir_mesh, path,
                      edotII, "edotII",
                      tauII, "tauII",
                      surf_normal_stress, "surf_normal_stress",
                      NULL);
    }

    /* destroy */
    rhea_viscosity_destroy (viscosity);
    rhea_velocity_destroy (velocity);
    ymir_vec_destroy (edotII);
    ymir_vec_destroy (tauII);
    ymir_vec_destroy (surf_normal_stress);
    ymir_velocity_elem_destroy (vel_elem);
  }

    /* compute and output analysis of stress */
  if (vtk_write_traction_path != NULL)  {
    ymir_vec_t         *velocity = rhea_velocity_new (ymir_mesh);
    ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);
    ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (
                                                    ymir_mesh->ma->N, ymir_mesh->ma->ompsize);
    ymir_vec_t         *vert_n_tau = ymir_dvec_new (ymir_mesh, 1,
                                                     YMIR_GAUSS_NODE);
    ymir_vec_t         *vert_s_tau = ymir_dvec_new (ymir_mesh, 1,
                                                     YMIR_GAUSS_NODE);
    ymir_vec_t         *hori_n_tau = ymir_dvec_new (ymir_mesh, 1,
                                                     YMIR_GAUSS_NODE);
    ymir_vec_t         *hori_s_tau = ymir_dvec_new (ymir_mesh, 1,
                                                     YMIR_GAUSS_NODE);
    ymir_vec_t         *vert_traction = ymir_dvec_new (ymir_mesh, 3,
                                                     YMIR_GAUSS_NODE);
    ymir_vec_t         *hori_traction = ymir_dvec_new (ymir_mesh, 3,
                                                     YMIR_GAUSS_NODE);

    ymir_stokes_vec_get_velocity (sol_vel_press, velocity,
                                    press_elem);
    rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

    double vert_n_dir[3] = {0.0, 0.0, 1.0};
    slabs_normal_stress (velocity, vert_n_tau, vert_s_tau, vert_traction, vert_n_dir,
                            viscosity, vel_elem);
    double hori_n_dir[3] = {0.0, 1.0, 0.0};
    slabs_normal_stress (velocity, hori_n_tau, hori_s_tau, hori_traction, hori_n_dir,
                            viscosity, vel_elem);
    {
      char            path[BUFSIZ];

      snprintf (path, BUFSIZ, "%s", vtk_write_traction_path);
      ymir_vtk_write (ymir_mesh, path,
                      vert_traction, "vertical tau_vec",
                      vert_n_tau, "vertical normal tau",
                      vert_s_tau, "vertical shear tau",
                      hori_traction, "horizontal tau_vec",
                      hori_n_tau, "horizontal normal tau",
                      hori_s_tau, "horizontal shear tau",
                      NULL);
    }

    /* destroy */
    rhea_viscosity_destroy (viscosity);
    rhea_velocity_destroy (velocity);
    ymir_vec_destroy (vert_traction);
    ymir_vec_destroy (vert_n_tau);
    ymir_vec_destroy (vert_s_tau);
    ymir_vec_destroy (hori_traction);
    ymir_vec_destroy (hori_n_tau);
    ymir_vec_destroy (hori_s_tau);
    ymir_velocity_elem_destroy (vel_elem);
  }

  /* destroy */
  rhea_velocity_pressure_destroy (sol_vel_press);

  /*
   * Finalize
   */

  /* destroy Stokes problem and mesh */
  slabs_setup_clear_all (stokes_problem, p4est, ymir_mesh, press_elem,
                           &discr_options);

  /* destroy options */
  ymir_options_global_destroy ();

  /* print that this function is ending */
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

  /* finalize rhea */
  rhea_finalize ();

  return 0;
}
