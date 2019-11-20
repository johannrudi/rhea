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
#include <ymir_interp_vec.h>
#include <ymir_mass_vec.h>
#include <ymir_vtk.h>
#include <ymir.h>

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

/* viscosity */
#define VISCOSITY_LITHOSPHERE (100.0)
#define VISCOSITY_ASTHENOSPHERE (1.0)
#define VISCOSITY_LITHOSPHERE_RADIUS_LOCATION (0.9)
#define VISCOSITY_ASTHENOSPHERE_RADIUS_LOCATION (0.2)
#define VISCOSITY_SLAB_RADIUS_LOCATION (0.7)
#define VISCOSITY_SLAB_WIDTH (0.2)

/* Dirichlet velocity B.C. */
#define COLLIDE_FLOW_SCALE (10.0)
#define COLLIDE_ZERO_POINT_LOCATION_MIDDLE (0.5)
#define COLLIDE_ZERO_POINT_LOCATION_UPPER (0.5)
#define COLLIDE_ZERO_POINT_LOCATION_LOWER (0.5)

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

double              center_x;
double              center_y;
double              center_z;
double              diameter;
double              sinker_scaling;
double              sinker_decay;

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

/* initialize slabs options: viscosity */
double              visc_lith =
  VISCOSITY_LITHOSPHERE;
double              visc_asthen =
  VISCOSITY_ASTHENOSPHERE;
double              visc_z_lith =
  VISCOSITY_LITHOSPHERE_RADIUS_LOCATION;
double              visc_z_asthen =
  VISCOSITY_ASTHENOSPHERE_RADIUS_LOCATION;
double              visc_z_slab =
  VISCOSITY_SLAB_RADIUS_LOCATION;
double              visc_slab_width =
  VISCOSITY_SLAB_WIDTH;

/* initialize slabs options: velocity boundary condition */
double              flow_scale =
  COLLIDE_FLOW_SCALE;
double              velocity_bc_middle =
  COLLIDE_ZERO_POINT_LOCATION_MIDDLE;
double              velocity_bc_upper =
  COLLIDE_ZERO_POINT_LOCATION_UPPER;
double              velocity_bc_lower =
  COLLIDE_ZERO_POINT_LOCATION_LOWER;

typedef enum
{
  SINKER,
  SLAB,
  COLLIDE,
  DRAG,
  TEST_MANUFACTURED,
  CUSTOM,
  TESTNONE
}
slabs_buoyancy_type_t;


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
  SLABS_CURVE_ORIENT_BOTTOM_BACK,
  SLABS_CURVE_ORIENT_TOP_RIGHT,
  SLABS_CURVE_ORIENT_BOTTOM_LEFT
}
slabs_subd_edge_orient_enum_t;

typedef enum
{
  CUSTOM_THINBOX,
  CUSTOM_SINKER
}
slabs_temp_custom_t;

typedef struct slabs_custom_sinker
{
  double center_x;
  double center_y;
  double center_z;
  double diameter;
  double scaling;
  double decay;
}
slabs_custom_sinker_t;

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

  slabs_temp_custom_t custom_type;
  slabs_custom_sinker_t *sinker_options;
}
slabs_temp_options_t;

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
  SLABS_VEL_DIR_BC_INOUTFLOW_TANH_TWOLAYER,
  SLABS_VEL_DIR_BC_INOUTFLOW_TANH_THREELAYER,
  SLABS_VEL_DIR_BC_INOUTFLOW_DOUBLE_TANH_THREELAYER,
  SLABS_VEL_DIR_BC_MANUFACTURED_SINCOS_ISO,
  SLABS_VEL_DIR_BC_MANUFACTURED_POLY_ANISO,
  SLABS_FREESURFACE,
  SLABS_FREESLIP
}
slabs_vel_dir_bc_t;

/* enumerator for viscosity types */
typedef enum
{
  SLABS_VISC_ISOTROPY,
  SLABS_VISC_TRANSVERSELY_ISOTROPY
}
slabs_viscosity_anisotropy_t;

typedef enum
{
  SLABS_VISC_LAYERS,
  SLABS_VISC_LAYERS_COUPLING
}
slabs_viscosity_geometry_t;


/* struct for viscosity options in slabs_options_t */
typedef struct slabs_visc_options
{
  slabs_viscosity_anisotropy_t  viscosity_anisotropy;
  slabs_viscosity_geometry_t        viscosity_geometry;
  double              z_lith;
  double              z_asthen;
  double              z_slab;
  double              slab_width;
  double              z_mantle;
  double              visc_lith;
  double              visc_asthen;
  double              visc_mantle;
}
slabs_visc_options_t;

/* enumerator for domain shapes */
typedef enum
{
  SLABS_X_FUNCTION_IDENTITY,
  SLABS_X_FUNCTION_SINE,
  SLABS_X_FUNCTION_PROFILE
}
slabs_x_func_t;

typedef struct slabs_topo_profile
{
  int               nsurf;
  double            *tX;
  double            *tY;
  double            *tZ;
}
slabs_topo_profile_t;

typedef struct slabs_surf_options
{
  slabs_x_func_t    x_func;
  slabs_topo_profile_t *topo_profile;
}
slabs_surf_options_t;

/* struct for velocity boundary condition options in slabs_options_t */
typedef struct slabs_velbc_options
{
  slabs_vel_dir_bc_t  vel_dir_bc;
  double              flow_scale;

  double              vel_dir_bc_middle;
  double              vel_dir_bc_upper;
  double              vel_dir_bc_lower;
}
slabs_velbc_options_t;

/* enumerator for tests */
typedef enum
{
  SLABS_TEST_STRESS_OP_NONE = 0,
  SLABS_TEST_STRESS_OP_SINCOS_SAME_OUTPUT,
  SLABS_TEST_STRESS_OP_SINCOS_ISO,
  SLABS_TEST_STRESS_OP_SINCOS_ANISO
}
slabs_test_stress_op_t;

typedef enum
{
  SLABS_TEST_MANUFACTURED_NONE = 0,
  SLABS_TEST_MANUFACTURED_SINCOS1_ISO,
  SLABS_TEST_MANUFACTURED_SINCOS1_TIROT90,
  SLABS_TEST_MANUFACTURED_SINCOS1_TIROT45,
  SLABS_TEST_MANUFACTURED_SINCOS1_TIROT60,
  SLABS_TEST_MANUFACTURED_SINCOS1_TIROT60_VISCEXP60,
  SLABS_TEST_MANUFACTURED_POLY1_TIROT90,
  SLABS_TEST_MANUFACTURED_POLY1_TIROT90_VISCEXP
}
slabs_test_manufactured_t;

typedef struct slabs_test_options
{
  slabs_test_stress_op_t    test_stress_op;
  slabs_test_manufactured_t test_manufactured;
  slabs_test_manufactured_t test_stress_comp;
}
slabs_test_options_t;

/* options of slabs example */
typedef struct slabs_options
{
  slabs_buoyancy_type_t    buoyancy_type;
  slabs_domain_options_t   * slabs_domain_options;
  slabs_temp_options_t   * slabs_temp_options;
  slabs_visc_options_t   * slabs_visc_options;
  slabs_weak_options_t   * slabs_weak_options;
  slabs_surf_options_t   * slabs_surf_options;
  slabs_velbc_options_t  * slabs_velbc_options;
  slabs_test_options_t   * slabs_test_options;
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
  poly_coeff = RHEA_ALLOC (double, 3);

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
    RHEA_ASSERT (*x != 0.0);
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
  RHEA_ASSERT (orientation_wrt_curve != NULL);

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
  closest_pt = RHEA_ALLOC (double, 2);
  if (line_orientation_start > 0) { /* if pt. is left of line through start pt.
                                     * with direction normal to start deriv. */
    proj = - (start_node - point_x) * tangent_start[0]
           - (start_val - point_y) * tangent_start[1];
    closest_pt[0] = start_node + proj * tangent_start[0];
    closest_pt[1] = start_val + proj * tangent_start[1];

    *orientation_wrt_curve = SLABS_CURVE_ORIENT_BOTTOM_BACK;
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

/* Computes distance and orientation of a point to a quadratic polynomial in 2D.*/
static double
slabs_compute_rot_on_poly2 (double *poly2_coeff, double *closest_pt)
{
  double              tangent[2];
  double              horiz[2] = {1.0, .0};
  double              rot;

  /* compute vectors of curve at start and end points */
  tangent[1] = slabs_poly2_deriv (closest_pt[0], poly2_coeff);
  tangent[0] = 1.0 / sqrt (1.0 + tangent[1] * tangent[1]);
  tangent[1] *= tangent[0];

  rot = acos(horiz[0]*tangent[0] + horiz[1]*tangent[1]);
  /* return distance */
  return rot;
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
  if (orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM ||
      orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_BACK) { /* if subd plate */
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
      || orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_BACK
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
  RHEA_FREE (closest_pt);

  /* return temperature */
  RHEA_ASSERT (isfinite (temp));
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
  RHEA_ASSERT (isfinite (*temp));

  /* bound temperature to valid interval */
  *temp = SC_MIN (1.0, *temp);
  *temp = SC_MAX (0.0, *temp);
}

static void
slabs_poly2_temperature_compute (ymir_vec_t *temperature,
                           slabs_options_t *slabs_options)
{
  const char         *this_fn_name = "slabs_poly2_temperature_compute";
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
  RHEA_FREE (poly2_coeff);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

void
drag_temperature_set_fn (double *temp, double x, double y, double z,
                          ymir_locidx_t nid, void *data)
{
  slabs_options_t    *slabs_options = (slabs_options_t *) data;
  double              subdu_lon, subdu_dip_angle, subdu_depth, subdu_width;
  double              start_node, start_val, start_deriv, end_node, end_val;
  slabs_weak_options_t *weak_options = slabs_options->slabs_weak_options;

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
  end_node = start_node + subdu_width / SLABS_MANTLE_DEPTH;
  end_val = start_val - subdu_depth / SLABS_MANTLE_DEPTH;

  if ((y<end_node) && (y>end_node-0.4) && (z<end_val) && (z>end_val-0.2) )
    *temp = 1.0;
  else
    *temp = 0.5;

  /* check temperature for `nan` and `inf` */
  RHEA_ASSERT (isfinite (*temp));

  /* bound temperature to valid interval */
  *temp = SC_MIN (1.0, *temp);
  *temp = SC_MAX (0.0, *temp);
}

void
thinbox_temperature_set_fn (double *temp, double x, double y, double z,
                          ymir_locidx_t nid, void *data)
{
  if ( (fabs(y-0.5) < 0.2)  && (z > 0.95 && z < 0.97))
    *temp = 0.0;
  else
    *temp = 0.5;

  /* check temperature for `nan` and `inf` */
  RHEA_ASSERT (isfinite (*temp));

  /* bound temperature to valid interval */
  *temp = SC_MIN (1.0, *temp);
  *temp = SC_MAX (0.0, *temp);
}

void
sinker_temperature_set_fn (double *temp, double x, double y, double z,
                          ymir_locidx_t nid, void *data)
{
  slabs_options_t *opt = data;
  slabs_custom_sinker_t *sinker = opt->slabs_temp_options->sinker_options;
  double center_x = sinker->center_x;
  double center_y = sinker->center_y;
  double center_z = sinker->center_z;
  double radius = 0.5 * sinker->diameter;
  double scaling = sinker->scaling;
  double decay = sinker->decay;
  double val;

  val = sqrt ( (center_x - x) * (center_x - x) +
               (center_y - y) * (center_y - y) +
               (center_z - z) * (center_z - z) );

  val = SC_MAX (0.0, val - radius);

  val = exp (-decay * val * val);
  RHEA_ASSERT (0.0 <= val && val <= 1.0);

  val = scaling * (1.0 - val);
  val = SC_MAX (0.0, SC_MIN (val, 1.0));
  *temp = val;
  *temp += (0.5 - scaling);

  /* check temperature for `nan` and `inf` */
  RHEA_ASSERT (isfinite (*temp));

  /* bound temperature to valid interval */
  *temp = SC_MIN (1.0, *temp);
  *temp = SC_MAX (0.0, *temp);
}


void
slabs_custom_temperature_compute (ymir_vec_t *temp,
                                    slabs_options_t *slabs_options)
{
  slabs_temp_custom_t type = slabs_options->slabs_temp_options->custom_type;

  switch (type)  {
    case CUSTOM_THINBOX:
      ymir_cvec_set_function (temp, thinbox_temperature_set_fn, slabs_options);
    break;

    case CUSTOM_SINKER:
      ymir_cvec_set_function (temp, sinker_temperature_set_fn, slabs_options);
    break;

    default:
      RHEA_ABORT_NOT_REACHED ();
  }

}


/**************************************
 * Weakzone Computation
 *************************************/

/* Computes the weak zone factor depending on the distance to a weak zone.
 * The edges of the weak zone are smoothed by a Gaussian.
 * 1 - (1 - weak_factor) * exp ( - dist^2 / (2 * (0.5*thickness)^2) )
 * I added a scale of 20.0 so that the smoothing is sharper, at designited width, weak is aoubt exp(-10)~1-0.001~0.999*/
static inline double
slabs_weakzone_factor_fn (const double distance,
                          const double thickness,
                          const double thickness_const,
                          const double weak_factor)
{
  const double        d = distance - 0.5 * thickness_const;
  const double        std_dev = 0.5 * (thickness - thickness_const);

  RHEA_ASSERT (thickness_const <= thickness);

  if (d <= 0.0) {
    /* return value inside zone with constant weak factor */
    return weak_factor;
  }
  else {
    /* return smoothed weak zone */
    return 1.0 - (1.0 - weak_factor) * exp (-20.0 * d*d / (2.0 * std_dev*std_dev));
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
      || orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM
      || orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_BACK ) {
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
    RHEA_ABORT_NOT_REACHED ();
  }

  /* destroy */
  RHEA_FREE (closest_pt);

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

  double              consider_thickness, courtesy_width;
  double              y_min = slabs_options->slabs_domain_options->y_min;
  double              start_node, start_val, start_deriv, end_node, end_val;
  double              *poly2_coeff;
  double              dist;
  double              weak = 1.0;
  double              s, a, b;

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
  RHEA_ASSERT (0.0 < subdu_lon);
  RHEA_ASSERT (0.0 < subdu_depth && 0.0 < subdu_width);
  RHEA_ASSERT (0.0 < subdu_thickness);
  RHEA_ASSERT (subdu_thickness_const <= subdu_thickness);
  RHEA_ASSERT (0.0 < subdu_weak_factor && subdu_weak_factor <= 1.0);
  RHEA_ASSERT (0.0 < ridge_depth && 0.0 < ridge_width);
  RHEA_ASSERT (0.0 <= ridge_smoothwidth);
  RHEA_ASSERT (0.0 < ridge_weak_factor && ridge_weak_factor <= 1.0);

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
  consider_thickness = (2.0 * subdu_thickness)
                    / SLABS_MANTLE_DEPTH;
  if (   (  start_node -  0.5 * consider_thickness
          / sin (subdu_dip_angle / 180.0 * M_PI) ) <= lon
      && lon <= (end_node + 0.5 * consider_thickness )
      && (end_val - 0.5 * consider_thickness ) <= r ) {

    /* compute distance to subduction weak zone */
    dist = slabs_weakzone_subduct_dist_2plates_poly2 (r, lon, poly2_coeff,
                                                      start_node, start_val, start_deriv,
                                                      end_node, end_val);
    if (r > 0.95) {
      /*from r=0.9 to r=1, shrink the width of weakzone by a factor of 2*/
      s = - 10.0 * r + 10.5;
      a = subdu_thickness * s;
      b = subdu_thickness_const * s;
      weak = slabs_weakzone_factor_fn (dist, a,
                                       b, subdu_weak_factor);
    }
    else {
    /* set weak zone factor */
    weak = slabs_weakzone_factor_fn (dist, subdu_thickness,
                                     subdu_thickness_const, subdu_weak_factor);
    }
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
  double         z_mid = slabs_options->slabs_visc_options->z_lith;

  /* compute weak zone factor for each node */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) { /* loop over all
                                                         * nodes */
    weak_elem[nodeid] = slabs_weakzone_node (x[nodeid], y[nodeid], z[nodeid],
                                             slabs_options);

    /*temporary, to make sure that the viscosity along the weak zone does not have a sharp change
     * at the lith-athen boundary*/
    if (z[nodeid] < z_mid)  {
      weak_elem[nodeid] = 1.0;
    }
  }
}

void
slabs_weakzone_compute (ymir_dvec_t *weakzone, void *data)
{
  const char         *this_fn_name = "slabs_weakzone_compute";
  slabs_options_t    *slabs_options = data;
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
  RHEA_FREE (poly2_coeff);  //TO FIND OUT: if I use RHEA_FREE, there is memory leak stderr

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**************************************
 * Viscosity Computation
 *************************************/
static void
slabs_layers_viscosity_old_elem (double *_sc_restrict visc_elem,
                            const double *_sc_restrict x,
                            const double *_sc_restrict y,
                            const double *_sc_restrict z,
                            const double *_sc_restrict weak_elem,
                            const int n_nodes_per_el,
                            slabs_options_t *slabs_options)
{
  int                 nodeid;
  double              z_mid = slabs_options->slabs_visc_options->z_lith;
  double              visc_lith = slabs_options->slabs_visc_options->visc_lith;
  double              visc_asthen = slabs_options->slabs_visc_options->visc_asthen;
  double              visc_smooth;
  double              factor = 1.0;
  slabs_topo_profile_t *topo = slabs_options->slabs_surf_options->topo_profile;
  double *tX = topo->tX;
  double *tY = topo->tY;
  double *tZ = topo->tZ;
  int     m, nsurf = topo->nsurf;

  /* compute viscosity in this element */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    if (topo->tZ != NULL) {
      for (m = 0; m < nsurf; m++)  {
        if (fabs(x[nodeid] - tX[m]) < SC_1000_EPS &&
            fabs(y[nodeid] - tY[m]) < SC_1000_EPS)  {
          factor = tZ[m];
        }
      }
    }
    visc_smooth = 10.0 * (z[nodeid]/factor - 0.45) * (visc_lith - visc_asthen) + visc_asthen;
    if (z[nodeid]/factor >= z_mid)  {
      visc_elem[nodeid] = (z[nodeid]/factor-z_mid)>0.05? visc_lith: visc_smooth;
    }
    else {
      visc_elem[nodeid] = (z[nodeid]/factor-z_mid)<-0.05? visc_asthen: visc_smooth;
    }
    visc_elem[nodeid] *= weak_elem[nodeid];

    /* check viscosity for `nan`, `inf`, and positivity */
    RHEA_ASSERT (isfinite (visc_elem[nodeid]));
    RHEA_ASSERT (0.0 < visc_elem[nodeid]);
  }
}

static void
slabs_layers_viscosity_elem (double *_sc_restrict visc_elem,
                            const double *_sc_restrict x,
                            const double *_sc_restrict y,
                            const double *_sc_restrict z,
                            const double *_sc_restrict weak_elem,
                            const int n_nodes_per_el,
                            slabs_options_t *slabs_options)
{
  int                 nodeid;
  double              z_mid = slabs_options->slabs_visc_options->z_lith;
  double              visc_lith = slabs_options->slabs_visc_options->visc_lith;
  double              visc_asthen = slabs_options->slabs_visc_options->visc_asthen;
  int                 m;

  /* compute viscosity in this element */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    if (z[nodeid] >= z_mid)  {
      visc_elem[nodeid] = visc_lith;
    }
    else {
      visc_elem[nodeid] = visc_asthen;
    }
    visc_elem[nodeid] *= weak_elem[nodeid];

    /* check viscosity for `nan`, `inf`, and positivity */
    RHEA_ASSERT (isfinite (visc_elem[nodeid]));
    RHEA_ASSERT (0.0 < visc_elem[nodeid]);
  }
}

static void
slabs_layers_coupling_viscosity_elem (double *_sc_restrict visc_elem,
                                      const double *_sc_restrict x,
                                      const double *_sc_restrict y,
                                      const double *_sc_restrict z,
                                      const double *_sc_restrict weak_elem,
                                      const int n_nodes_per_el,
                                      slabs_options_t *slabs_options)
{
  int                 nodeid;
  double              z_mid = slabs_options->slabs_visc_options->z_lith;
  double              z_slab = slabs_options->slabs_visc_options->z_slab;
  double              slab_width = slabs_options->slabs_visc_options->slab_width;
  double              y_center = 0.5;

  /* compute viscosity in this element */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    if (z[nodeid] >= z_mid)  {
      visc_elem[nodeid] = slabs_options->slabs_visc_options->visc_lith;
    }
    else if (z[nodeid] >= z_slab && fabs(y[nodeid] - y_center) <= 0.5 * slab_width) {
      visc_elem[nodeid] = slabs_options->slabs_visc_options->visc_lith;
    }
    else {
      visc_elem[nodeid] = slabs_options->slabs_visc_options->visc_asthen;
    }
    visc_elem[nodeid] *= weak_elem[nodeid];

    /* check viscosity for `nan`, `inf`, and positivity */
    RHEA_ASSERT (isfinite (visc_elem[nodeid]));
    RHEA_ASSERT (0.0 < visc_elem[nodeid]);
  }
}

static void
slabs_viscosity_compute (ymir_vec_t *viscosity,
                         ymir_vec_t *rank1_tensor_scal,
                         ymir_vec_t *bounds_marker,
                         ymir_vec_t *yielding_marker,
                         ymir_vec_t *temperature,
                         ymir_vec_t *weakzone,
                         ymir_vec_t *velocity,
                         void *data)
{
  slabs_options_t  *slabs_options = data;
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (viscosity);
  const ymir_locidx_t  n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);
  mangll_t           *mangll = mesh->ma;
  const int           N = ymir_n (mangll->N);

  sc_dmatrix_t       *visc_el_mat, *weak_el_mat;
  double             *x, *y, *z, *tmp_el,*visc_el_data, *weak_el_data;
  ymir_locidx_t       elid;

  /* create work variables */
  weak_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  visc_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  visc_el_data = visc_el_mat->e[0];
  x = RHEA_ALLOC (double, n_nodes_per_el);
  y = RHEA_ALLOC (double, n_nodes_per_el);
  z = RHEA_ALLOC (double, n_nodes_per_el);
  tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) {
    /* get coordinates of this element at Gauss nodes */
    ymir_mesh_get_elem_coord_gauss (x, y, z, elid, mesh, tmp_el);

    weak_el_data = rhea_viscosity_get_elem_gauss (weak_el_mat, weakzone, elid);
    /* compute user defined weak zone viscosity*/
    switch (slabs_options->slabs_visc_options->viscosity_geometry) {
      case SLABS_VISC_LAYERS:
        slabs_layers_viscosity_elem (visc_el_data, x, y, z, weak_el_data, n_nodes_per_el,
                                     slabs_options);
      break;

      case SLABS_VISC_LAYERS_COUPLING:
        slabs_layers_coupling_viscosity_elem (visc_el_data, x, y, z, weak_el_data, n_nodes_per_el,
                                             slabs_options);
      break;

      default:
      RHEA_ABORT_NOT_REACHED ();
    }

    /* set viscosity of this element */
    rhea_viscosity_set_elem_gauss (viscosity, visc_el_mat, elid);
  }

  /* destroy */
  sc_dmatrix_destroy (weak_el_mat);
  sc_dmatrix_destroy (visc_el_mat);
  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);
}

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

  double              consider_thickness, courtesy_width;
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
  RHEA_ASSERT (0.0 < subdu_lon);
  RHEA_ASSERT (0.0 < subdu_dip_angle);
  RHEA_ASSERT (0.0 < subdu_depth && 0.0 < subdu_width);
  RHEA_ASSERT (0.0 < subdu_thickness);
  RHEA_ASSERT (subdu_thickness_const <= subdu_thickness);
  RHEA_ASSERT (0.0 < ridge_depth && 0.0 < ridge_width);
  RHEA_ASSERT (0.0 <= ridge_smoothwidth);

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
  consider_thickness = (2.0 * subdu_thickness)
                    / SLABS_MANTLE_DEPTH;
  if (   (  start_node - 0.5 * consider_thickness
          / sin (subdu_dip_angle / 180.0 * M_PI)) <= lon
      && lon <= (end_node + 0.5 * consider_thickness)
      && (end_val - consider_thickness) <= r ) {
    /*
     * compute rotation of subduction weak zone from y=0 axis
     * */

    /* compute closest point on curve and orientation w.r.t. curve */
    closest_pt = slabs_compute_closest_pt_on_poly2 (lon, r,
                                                  poly2_coeff, start_node,
                                                  start_val, start_deriv,
                                                  end_node, end_val,
                                                  &orientation_wrt_curve);

    if (orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_BACK)
      closest_pt[0] = start_node;
    else if (orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_LEFT ||
             orientation_wrt_curve == SLABS_CURVE_ORIENT_TOP_RIGHT)
      closest_pt[0] = end_node;

    rot = slabs_compute_rot_on_poly2 (poly2_coeff, closest_pt);

    RHEA_FREE (closest_pt);
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
    if (fabs(1.0 - weak_elem[nodeid]) < 1e-3)  {
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
  RHEA_FREE (poly2_coeff);
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
  ymir_vec_t         *TI_weakzone = rhea_viscosity_new (ymir_mesh);

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* copy viscosity */
  rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

  /* compute the shear viscosity and rotation angles */
  slabs_weakzone_compute (TI_weakzone, slabs_options);
  slabs_TI_viscosity_compute (ymir_mesh, coeff_TI_svisc, viscosity, TI_weakzone, slabs_options);

  ymir_vec_scale (2.0, coeff_TI_svisc);

  slabs_TI_rotation_compute (TI_rotate, TI_weakzone, slabs_options);

  /* get the viscous stress operator */
  stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
  stress_op = stokes_op->stress_op;

 /* update viscous stress operator providing the anisotropic viscosity */
  ymir_stress_op_coeff_compute_TI_tensor (stress_op, coeff_TI_svisc,
                                          TI_rotate);
  /* destroy */
  rhea_viscosity_destroy (viscosity);
  rhea_viscosity_destroy (TI_weakzone);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/* setup the TI shear viscosity and tensor in stress operator */
void
slabs_stokes_problem_setup_TI_manufactured (ymir_mesh_t *ymir_mesh,
                                           rhea_stokes_problem_t *stokes_problem,
                                           slabs_options_t *slabs_options,
                                           ymir_vec_t *coeff_TI_svisc,
                                           ymir_vec_t *TI_rotate)
{
  const char         *this_fn_name = "slabs_stokes_problem_setup_TI_manufactured";
  ymir_stokes_op_t   *stokes_op;
  ymir_stress_op_t   *stress_op;
  ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);
  ymir_vec_t         *TI_weakzone = rhea_viscosity_new (ymir_mesh);
  const               slabs_test_manufactured_t
                      test_type = slabs_options->slabs_test_options->test_manufactured;
  double              rot, s_n_ratio = 0.2; /*eta_n=5, eta_s=1 */

  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (ymir_mesh);
  const ymir_locidx_t n_elements = ymir_mesh_get_num_elems_loc (ymir_mesh);
  sc_dmatrix_t       *weak_el_mat;
  double             *weak_el_data;
  double             *x, *y, *z, *tmp_el;
  ymir_locidx_t       elid;
  int                 nodeid;


  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* copy viscosity */
  rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

  /* compute the shear viscosity and rotation angles */
  switch (test_type) {
    case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT90:
    case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT45:
    case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT60:
    case SLABS_TEST_MANUFACTURED_POLY1_TIROT90:
      ymir_vec_set_value (TI_weakzone, s_n_ratio);
      break;

     /* eta_n=5, eta_s=s_n_ratio * eta_n * 0.5*(exp(y)+exp(z)) */
    case SLABS_TEST_MANUFACTURED_POLY1_TIROT90_VISCEXP:
      weak_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
      weak_el_data = weak_el_mat->e[0];
      x = RHEA_ALLOC (double, n_nodes_per_el);
      y = RHEA_ALLOC (double, n_nodes_per_el);
      z = RHEA_ALLOC (double, n_nodes_per_el);
      tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

      for (elid = 0; elid < n_elements; elid++) {
        ymir_mesh_get_elem_coord_gauss (x, y, z, elid, ymir_mesh, tmp_el);
        for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
          weak_el_data[nodeid] = s_n_ratio * 0.5 * (exp(y[nodeid]) + exp(z[nodeid]));
        }
        rhea_viscosity_set_elem_gauss (TI_weakzone, weak_el_mat, elid);
      }
      sc_dmatrix_destroy (weak_el_mat);
      RHEA_FREE (x);
      RHEA_FREE (y);
      RHEA_FREE (z);
      RHEA_FREE (tmp_el);
      break;

    case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT60_VISCEXP60:
      weak_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
      weak_el_data = weak_el_mat->e[0];
      x = RHEA_ALLOC (double, n_nodes_per_el);
      y = RHEA_ALLOC (double, n_nodes_per_el);
      z = RHEA_ALLOC (double, n_nodes_per_el);
      tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

      for (elid = 0; elid < n_elements; elid++) {
        ymir_mesh_get_elem_coord_gauss (x, y, z, elid, ymir_mesh, tmp_el);
        for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
          weak_el_data[nodeid] = s_n_ratio * exp(0.5 * (sqrt(3.0) * y[nodeid] + z[nodeid]));
        }
        rhea_viscosity_set_elem_gauss (TI_weakzone, weak_el_mat, elid);
      }
      sc_dmatrix_destroy (weak_el_mat);
      RHEA_FREE (x);
      RHEA_FREE (y);
      RHEA_FREE (z);
      RHEA_FREE (tmp_el);
      break;


    default: /* BC not set */
      RHEA_ABORT_NOT_REACHED ();
  }
  slabs_TI_viscosity_compute (ymir_mesh, coeff_TI_svisc, viscosity, TI_weakzone, slabs_options);
  ymir_vec_scale (2.0, coeff_TI_svisc);

  /* rotation angle */
  switch (test_type) {
    case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT90:
    case SLABS_TEST_MANUFACTURED_POLY1_TIROT90:
    case SLABS_TEST_MANUFACTURED_POLY1_TIROT90_VISCEXP:
      rot = 0.5 * M_PI;
      break;

    case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT45:
      rot = 0.25 * M_PI;
      break;

    case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT60:
    case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT60_VISCEXP60:
      rot = 1.0/3.0 * M_PI;
      break;

    default: /* BC not set */
      RHEA_ABORT_NOT_REACHED ();
  }
  ymir_vec_set_value (TI_rotate, rot);

  /* get the viscous stress operator */
  stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
  stress_op = stokes_op->stress_op;

 /* update viscous stress operator providing the anisotropic viscosity */
  ymir_stress_op_coeff_compute_TI_tensor (stress_op, coeff_TI_svisc,
                                          TI_rotate);
  /* destroy */
  rhea_viscosity_destroy (viscosity);
  rhea_viscosity_destroy (TI_weakzone);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**************************************
 * Geometry Transformations
 *************************************/

static void
slabs_X_fn_identity (mangll_tag_t tag, mangll_locidx_t np,
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

static void
slabs_X_fn_sine (mangll_tag_t tag, mangll_locidx_t np,
                       const double *_sc_restrict EX,
                       const double *_sc_restrict EY,
                       const double *_sc_restrict EZ,
                       double *_sc_restrict X,
                       double *_sc_restrict Y,
                       double *_sc_restrict Z, void *data)
{
  mangll_locidx_t     il;
  const double        slope = 0.5;

  for (il = 0; il < np; ++il) {
    X[il] = EX[il];
    Y[il] = EY[il];
    Z[il] = EZ[il] * (1 + slope * sin(M_PI * EY[il]));
  }
}

static void
slabs_X_fn_function (mangll_tag_t tag, mangll_locidx_t np,
                       const double *_sc_restrict EX,
                       const double *_sc_restrict EY,
                       const double *_sc_restrict EZ,
                       double *_sc_restrict X,
                       double *_sc_restrict Y,
                       double *_sc_restrict Z, void *data)
{
  mangll_locidx_t     il;
  int           k;
  double        factor;
//  double        coeff[11] = {1.0146, -0.0030374, -0.11728, -1.788, 14.385,
//                            -64.742, 177.44, -283.49, 258.25, -124.93, 24.987};
  double        coeff[11] = {1.0511, 0.97075, -14.253, 78.626, -214.55,
                              324.36, -289.11, 154.51, -47.999, 7.7874, -0.48331};


  for (il = 0; il < np; ++il) {
    X[il] = EX[il];
    Y[il] = EY[il];
    factor = .0;
    for (k = 0; k <= 10; k++) {
      factor += coeff[k] * pow(EY[il],(double) k);
    }
    Z[il] = EZ[il] * factor;
  }
}

/*distort mesh following the imposed topography on the surface
 * the topography information is given in void *data */
static void
slabs_X_fn_profile (mangll_tag_t tag, mangll_locidx_t np,
                   const double *_sc_restrict EX,
                   const double *_sc_restrict EY,
                   const double *_sc_restrict EZ,
                   double *_sc_restrict X,
                   double *_sc_restrict Y,
                   double *_sc_restrict Z, void *data)
{
  mangll_locidx_t     il;
  double factor;
  slabs_topo_profile_t *topo = (slabs_topo_profile_t *) data;
  double *tX = topo->tX;
  double *tY = topo->tY;
  double *tZ = topo->tZ;
  int     m, nsurf = topo->nsurf;

  /*loop for nodes in each element*/
  for (il = 0; il < np; ++il) {
    X[il] = EX[il];
    Y[il] = EY[il];
    /*loop over all the topography information
     * to find the corresponding distortion factor */
    for (m = 0; m < nsurf; m++)  {
      if (fabs(EX[il] - tX[m]) < SC_1000_EPS &&
          fabs(EY[il] - tY[m]) < SC_1000_EPS)  {
        factor = tZ[m];
      }
    }
    Z[il] = EZ[il] * factor;
  }
}

#if 0
static void
slabs_X_fn_profile_interpolation (mangll_tag_t tag, mangll_locidx_t np,
                   const double *_sc_restrict EX,
                   const double *_sc_restrict EY,
                   const double *_sc_restrict EZ,
                   double *_sc_restrict X,
                   double *_sc_restrict Y,
                   double *_sc_restrict Z, void *data)
{
  mangll_locidx_t     il;
  double factor;
  slabs_topo_profile_t *topo = (slabs_topo_profile_t *) data;
  double *tX = topo->tX;
  double *tY = topo->tY;
  double *tZ = topo->tZ;
  int     m, nsurf = topo->nsurf;

  /*loop for nodes in each element*/
  for (il = 0; il < np; ++il) {
    X[il] = EX[il];
    Y[il] = EY[il];
    /*loop over all the topography information
     * to find the corresponding distortion factor */
    for (m = 0; m < nsurf; m++)  {
      if ((tX[m] - EX[il]) >= SC_1000_EPS)
        continue;
      else if (fabs(EX[il] - tX[m]) < SC_1000_EPS) {
        if ((tY[m] - EY[il]) >= SC_1000_EPS)
          continue;
        else if (fabs(EY[il] - tY[m]) < SC_1000_EPS)
          factor = tZ[m];
          break;
        }
        else {
          interpolation;
          factor = **;
          break;
        }
      }
      else {

        interpolation of tX
      }

    }
    Z[il] = EZ[il] * factor;
  }
}

#endif

void
slabs_surface_location (slabs_options_t *slabs_options,
                       rhea_discretization_options_t *discr_options)
{
 slabs_surf_options_t  *slabs_surf_options = slabs_options->slabs_surf_options;
 slabs_topo_profile_t *topo = slabs_surf_options->topo_profile;

    /* set custom X-function */
  switch (slabs_surf_options->x_func) {
    case SLABS_X_FUNCTION_IDENTITY:
      rhea_discretization_set_X_fn (discr_options, slabs_X_fn_identity, NULL);
      break;

    case SLABS_X_FUNCTION_SINE:
      rhea_discretization_set_X_fn (discr_options, slabs_X_fn_sine, NULL);
      break;

    case SLABS_X_FUNCTION_PROFILE:
      rhea_discretization_set_X_fn (discr_options, slabs_X_fn_profile, topo);
      break;

    default:
      RHEA_ABORT_NOT_REACHED ();
    break;
  }
}

void
slabs_coordx_set_fn (double *coordx,
                     double x, double y, double z,
                     double nx, double ny, double nz,
                     ymir_topidx_t face, ymir_locidx_t nid,
                     void *data)
{
  *coordx = x;
}

void
slabs_coordy_set_fn (double *coordy,
                     double x, double y, double z,
                     double nx, double ny, double nz,
                     ymir_topidx_t face, ymir_locidx_t nid,
                     void *data)
{
  *coordy = y;
}

/**************************************************************
 * Non-zero Dirichlet boundary conditions
***************************************************************/

/* Dirichlet all */
static ymir_dir_code_t
slabs_set_vel_dir_all (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
     return YMIR_VEL_DIRICHLET_ALL;
}

static ymir_dir_code_t
slabs_set_vel_dir_all_2D (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
  if (face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE3 ||
      face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE4 ||
      face == RHEA_DOMAIN_BOUNDARY_FACE_BASE  ||
      face == RHEA_DOMAIN_BOUNDARY_FACE_TOP) {
      return YMIR_VEL_DIRICHLET_ALL;
  }
  else
    return YMIR_VEL_DIRICHLET_NORM;
}

/* free-slip */
static ymir_dir_code_t
slabs_set_vel_dir_freeslip (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
    return YMIR_VEL_DIRICHLET_NORM;
}

/* free-surface */
static ymir_dir_code_t
slabs_set_vel_dir_freesurface (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
  if (face == RHEA_DOMAIN_BOUNDARY_FACE_TOP)
    return YMIR_VEL_DIRICHLET_NONE;
  else
    return YMIR_VEL_DIRICHLET_NORM;
}



/* Dirichlet all on one side of the domain: SIDE3 (y=0).*/
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
  else {
    return YMIR_VEL_DIRICHLET_NORM;
  }
}

/* Dirichlet all on both left and right sides of the domain: SIDE3 (y=0), and SIDE4 (y=ymax).*/
static ymir_dir_code_t
slabs_set_vel_dir_inoutflow_double (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
  if (face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE3 ||
      face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE4) {
    return YMIR_VEL_DIRICHLET_ALL;
  }
  else {
    return YMIR_VEL_DIRICHLET_NORM;
  }
}

/* Dirichlet all on one side of the domain: SIDE3 (y=0), and Neumann at the base*/
static ymir_dir_code_t
slabs_set_vel_dir_inoutflow_basefree (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
  if (face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE3) {
    return YMIR_VEL_DIRICHLET_ALL;
  }
  else if (face == RHEA_DOMAIN_BOUNDARY_FACE_BASE) {
     return YMIR_VEL_DIRICHLET_NONE;
  }
  else {
    return YMIR_VEL_DIRICHLET_NORM;
  }
}

/* In-out flow sine velocity on one side of the domain. */
void
slabs_set_rhs_vel_nonzero_dir_inoutflow_sin (
    double *vel, double x, double y, double z,
    ymir_locidx_t nid, void *data)
{
  slabs_options_t  *slabs_options = data;
  const double     flow_scale = slabs_options->slabs_velbc_options->flow_scale;

  if (y < SC_1000_EPS) {
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

/* In-out flow tanh velocity on one side of the domain. 2 layers */
void
slabs_set_rhs_vel_nonzero_dir_inoutflow_tanh_2layer (double *vel, double x, double y, double z,
                                              ymir_locidx_t nid, void *data)
{
  slabs_options_t  *slabs_options = data;
  slabs_velbc_options_t *velbc_options = slabs_options->slabs_velbc_options;
  const double        z_max = slabs_options->slabs_domain_options->z_max;
  const double        flow_scale = velbc_options->flow_scale;
  const double        zM = velbc_options->vel_dir_bc_middle;
  const double        a = z_max-zM, b = z_max-a;
  const double        shape = 2.0 * M_PI, scaling = 0.5*flow_scale, shift = 0.5*(b-a)*flow_scale;
  double txM = shape * (z-zM);

  if (y < SC_1000_EPS) {
    vel[0] = 0.0;
    vel[2] = 0.0;
    vel[1] = shift + scaling *
             ( (exp (txM) - exp (-txM)) /
             (exp (txM) + exp (-txM)) );
  }
  else {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }
}

/* In-out flow tanh velocity on one side of the domain.*/
void
slabs_set_rhs_vel_nonzero_dir_inoutflow_tanh_3layer (double *vel, double x, double y, double z,
                                              ymir_locidx_t nid, void *data)
{
  slabs_options_t  *slabs_options = data;
  slabs_velbc_options_t *velbc_options = slabs_options->slabs_velbc_options;
  const double        z_max = slabs_options->slabs_domain_options->z_max;
  const double        flow_scale = velbc_options->flow_scale;
  const double        zU = velbc_options->vel_dir_bc_upper;
  const double        zL = velbc_options->vel_dir_bc_lower;
  const double        a = zU-zL, b = z_max-a, c = 0.5*(zU+zL);
  const double        shape = 2.0 * M_PI, scaling = 0.5*(b+a)*flow_scale, shift = 0.5*(b-a)*flow_scale;
  double txL = shape*(z-zL),txU = shape*(z-zU);

  if (y < SC_1000_EPS) {
    vel[0] = 0.0;
    vel[2] = 0.0;
    if (z<=c)
      vel[1] = shift + scaling *
             ( (exp (txL) - exp (-txL)) /
             (exp (txL) + exp (-txL)) );
    else
      vel[1] = shift - scaling *
             ( (exp (txU) - exp (-txU)) /
             (exp (txU) + exp (-txU)) );
  }
  else {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }
}

/* In-out flow tanh velocity on both left and right sides of the domain.*/
void
slabs_set_rhs_vel_nonzero_dir_inoutflow_double_tanh_3layer (
    double *vel, double x, double y, double z,
    ymir_locidx_t nid, void *data)
{
  slabs_options_t  *slabs_options = data;
  slabs_velbc_options_t *velbc_options = slabs_options->slabs_velbc_options;
  const double        z_max = slabs_options->slabs_domain_options->z_max;
  const double        y_max = slabs_options->slabs_domain_options->y_max;
  const double        flow_scale = velbc_options->flow_scale;
  const double        zU = velbc_options->vel_dir_bc_upper;
  const double        zL = velbc_options->vel_dir_bc_lower;
  const double        a = zU-zL, b = z_max-a, c = 0.5*(zU+zL);
  const double        shape = 2.0 * M_PI, scaling = 0.5*(b+a)*flow_scale, shift = 0.5*(b-a)*flow_scale;
  double txL = shape*(z-zL),txU = shape*(z-zU);

  if (y < SC_1000_EPS) {
    vel[0] = 0.0;
    vel[2] = 0.0;
    if (z<=c)
      vel[1] = shift + scaling *
               ( (exp (txL) - exp (-txL)) /
               (exp (txL) + exp (-txL)) );
    else
      vel[1] = shift - scaling *
               ( (exp (txU) - exp (-txU)) /
               (exp (txU) + exp (-txU)) );
  }
  else if ((y_max - y) < SC_1000_EPS) {
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
slabs_vel_nonzero_dirichlet_compute (ymir_vec_t * rhs_vel_nonzero_dirichlet,
                                     void * data)
{
  slabs_options_t    *slabs_options = data;

  switch (slabs_options->slabs_velbc_options->vel_dir_bc) {
    case SLABS_VEL_DIR_BC_INOUTFLOW_SIN:
      rhea_domain_set_user_velocity_dirichlet_bc (
          slabs_set_vel_dir_inoutflow, NULL /* no data necessary */,
          0 /* TODO don't need this flag */);
      ymir_cvec_set_function (rhs_vel_nonzero_dirichlet,
                              slabs_set_rhs_vel_nonzero_dir_inoutflow_sin,
                              slabs_options);
      break;

    case SLABS_VEL_DIR_BC_INOUTFLOW_TANH_TWOLAYER:
      rhea_domain_set_user_velocity_dirichlet_bc (
          slabs_set_vel_dir_inoutflow, NULL /* no data necessary */,
          0 /* TODO don't need this flag */);
      ymir_cvec_set_function (rhs_vel_nonzero_dirichlet,
                              slabs_set_rhs_vel_nonzero_dir_inoutflow_tanh_2layer,
                              slabs_options);
      break;

    case SLABS_VEL_DIR_BC_INOUTFLOW_TANH_THREELAYER:
      rhea_domain_set_user_velocity_dirichlet_bc (
          slabs_set_vel_dir_inoutflow, NULL /* no data necessary */,
          0 /* TODO don't need this flag */);
      ymir_cvec_set_function (rhs_vel_nonzero_dirichlet,
                              slabs_set_rhs_vel_nonzero_dir_inoutflow_tanh_3layer,
                              slabs_options);
      break;


   case SLABS_FREESURFACE:
      rhea_domain_set_user_velocity_dirichlet_bc (
          slabs_set_vel_dir_freesurface, NULL /* no data necessary */,
          0 );
   break;

   case SLABS_FREESLIP:
      rhea_domain_set_user_velocity_dirichlet_bc (
          slabs_set_vel_dir_freeslip, NULL /* no data necessary */,
          0 );
   break;

   default: /* BC not set */
      RHEA_ABORT_NOT_REACHED ();
  }
}


/***********************************************************
 * Test using manufactured solution
 ***********************************************************/
/*This flow field is divergence free
 * use in both stress op test and manufactured solution test
 * for TI case with 90 degree rotation the -grad(u) is the same with that for ISO case*/
static void
slabs_test_sincos1_vel_in_fn (double * vel, double x, double y,
                                          double z, ymir_locidx_t nodeid,
                                          void *data)
{
  vel[0] = 0.0;
  vel[1] = + sin (M_PI * y) * cos (M_PI * z);
  vel[2] = - cos (M_PI * y) * sin (M_PI * z);
}

static void
slabs_test_sincos1_ISO_vel_out_fn (double * vel, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data)
{
  vel[0] = 0.0;
  vel[1] = +10.0 * M_PI * M_PI * sin (M_PI * y) * cos (M_PI * z);
  vel[2] = -10.0 * M_PI * M_PI * cos (M_PI * y) * sin (M_PI * z);
}

static void
slabs_test_sincos1_TIrot90_vel_out_fn (double * vel, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data)
{
  vel[0] = 0.0;
  vel[1] = +10.0 * M_PI * M_PI * sin (M_PI * y) * cos (M_PI * z);
  vel[2] = -10.0 * M_PI * M_PI * cos (M_PI * y) * sin (M_PI * z);
}


/* 2eta_s*edots, 2eta_n*edotn, along the 'fault plane'*/
static void
slabs_test_sincos1_TIrot90_traction_fn (double * trac, double x, double y,
                                          double z, ymir_locidx_t nodeid,
                                          void *data)
{
  double tt = 90.0 * M_PI / 180.0;
  double zz = - cos(tt) * cos(tt) + sin(tt) * sin(tt);
  double yz = 2 * sin(tt) * cos(tt);
  double nvisc = 5.0, svisc = 1.0;

  /*fn=trac1=10, fs=trac2=0 */
  trac[0] = .0;
  trac[1] = 2.0 * nvisc * zz;
  trac[2] = 2.0 * svisc * yz;
  trac[1] *= ( M_PI * cos(M_PI * y) * cos(M_PI * z) );
  trac[2] *= ( M_PI * cos(M_PI * y) * cos(M_PI * z) );
}

static void
slabs_test_sincos1_TIrot90_stress_fn (double * stress, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data)
{
  double tt = 90.0 * M_PI / 180.0;
  double a = - cos(tt) * cos(tt) + sin(tt) * sin(tt);
  double b = 2 * sin(tt) * cos(tt);
  double nvisc = 5.0, svisc = 1.0;

  stress[0] = nvisc * a * a + svisc * b * b;
  stress[1] = -stress[0];
  stress[2] = (nvisc - svisc) * a * b;
  stress[0] *= ( 2.0 * M_PI * cos(M_PI * y) * cos(M_PI * z) );
  stress[1] *= ( 2.0 * M_PI * cos(M_PI * y) * cos(M_PI * z) );
  stress[2] *= ( 2.0 * M_PI * cos(M_PI * y) * cos(M_PI * z) );
}

static void
slabs_test_sincos1_TIrot45_vel_out_fn (double * vel, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data)
{
  vel[0] = 0.0;
  vel[1] = sin (M_PI * y) * cos (M_PI * z);
  vel[2] = - cos (M_PI * y) * sin (M_PI * z);
  vel[1] *= 2 * M_PI * M_PI;
  vel[2] *= 2 * M_PI * M_PI;
}

static void
slabs_test_sincos1_TIrot45_stress_fn (double * stress, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data)
{
  double tt = 45.0 * M_PI / 180.0;
  double a = - cos(tt) * cos(tt) + sin(tt) * sin(tt);
  double b = 2 * sin(tt) * cos(tt);
  double nvisc = 5.0, svisc = 1.0;

  stress[0] = nvisc * a * a + svisc * b * b;
  stress[1] = -stress[0];
  stress[2] = (nvisc - svisc) * a * b;
  stress[0] *= 2.0 * M_PI * cos(M_PI * y) * cos(M_PI * z);
  stress[1] *= 2.0 * M_PI * cos(M_PI * y) * cos(M_PI * z);
  stress[2] *= 2.0 * M_PI * cos(M_PI * y) * cos(M_PI * z);
}

/* 2eta_s*edots, 2eta_n*edotn, along the 'fault plane'*/
static void
slabs_test_sincos1_TIrot45_traction_fn (double * trac, double x, double y,
                                          double z, ymir_locidx_t nodeid,
                                          void *data)
{
  double tt = 45.0 * M_PI / 180.0;
  double zz = - cos(tt) * cos(tt) + sin(tt) * sin(tt);
  double yz = 2 * sin(tt) * cos(tt);
  double nvisc = 5.0, svisc = 1.0;

  /*fn=trac1=0, fs=trac2=2 */
  trac[0] = .0;
  trac[1] = 2.0 * nvisc * zz;
  trac[2] = 2.0 * svisc * yz;
  trac[1] *= ( M_PI * cos(M_PI * y) * cos(M_PI * z) );
  trac[2] *= ( M_PI * cos(M_PI * y) * cos(M_PI * z) );
}

static void
slabs_test_sincos1_TIrot60_vel_out_fn (double * vel, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data)
{
  vel[0] = 0.0;
  vel[1] = 4.0 * sin (M_PI * y) * cos (M_PI * z) - 2.0 * sqrt(3.0) * cos (M_PI * y) * sin (M_PI * z);
  vel[2] = -2.0 * sqrt(3.0) * sin (M_PI * y) * cos (M_PI * z) - 4.0 * cos (M_PI * y) * sin (M_PI * z);
  vel[1] *= (M_PI * M_PI);
  vel[2] *= (M_PI * M_PI);
}

static void
slabs_test_sincos1_TIrot60_stress_fn (double * stress, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data)
{
  double tt = 60.0 * M_PI / 180.0;
  double a = - cos(tt) * cos(tt) + sin(tt) * sin(tt);
  double b = 2 * sin(tt) * cos(tt);
  double nvisc = 5.0, svisc = 1.0;

  stress[0] = nvisc * a * a + svisc * b * b;
  stress[1] = -stress[0];
  stress[2] = (nvisc - svisc) * a * b;
  stress[0] *= 2.0 * M_PI * cos(M_PI * y) * cos(M_PI * z);
  stress[1] *= 2.0 * M_PI * cos(M_PI * y) * cos(M_PI * z);
  stress[2] *= 2.0 * M_PI * cos(M_PI * y) * cos(M_PI * z);
}

/* 2eta_s*edots, 2eta_n*edotn, along the 'fault plane'*/
static void
slabs_test_sincos1_TIrot60_traction_fn (double * trac, double x, double y,
                                          double z, ymir_locidx_t nodeid,
                                          void *data)
{
  double tt = 60.0 * M_PI / 180.0;
  double zz = - cos(tt) * cos(tt) + sin(tt) * sin(tt);
  double yz = 2 * sin(tt) * cos(tt);
  double nvisc = 5.0, svisc = 1.0;

  /*fn=trac1=5, fs=trac2=sqrt(3) */
  trac[0] = .0;
  trac[1] = 2.0 * nvisc * zz;
  trac[2] = 2.0 * svisc * yz;
  trac[1] *= ( M_PI * cos(M_PI * y) * cos(M_PI * z) );
  trac[2] *= ( M_PI * cos(M_PI * y) * cos(M_PI * z) );
}

static void
slabs_test_sincos1_manufactured_set_velbc (double * vel, double x, double y,
                                              double z, ymir_locidx_t nodeid,
                                              void *data)
{
  slabs_options_t  *slabs_options = data;
  const double  z_max = slabs_options->slabs_domain_options->z_max;
  const double  y_max = slabs_options->slabs_domain_options->y_max;
  const double  x_max = slabs_options->slabs_domain_options->x_max;

  if (y < SC_1000_EPS || (y_max - y) < SC_1000_EPS ||
      z < SC_1000_EPS || (z_max - z) < SC_1000_EPS ||
      x < SC_1000_EPS || (x_max - x) < SC_1000_EPS)  {
    vel[0] = 0.0;
    vel[1] = + sin (M_PI * y) * cos (M_PI * z);
    vel[2] = - cos (M_PI * y) * sin (M_PI * z);
  }
  else {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }
}

/*This velocity field is not divergence-free,
 * currently only used in stress op test, not manufactured solution test*/
static void
slabs_test_sincos2_vel_in_fn (double * vel, double x, double y,
                                        double z, ymir_locidx_t nodeid,
                                        void *data)
{
  vel[0] = 0.0;
  vel[1] = sin (M_PI * y) * cos (M_PI * z);
  vel[2] = cos (M_PI * y) * sin (M_PI * z);
}

static void
slabs_test_sincos2_TIrot90_vel_out_fn (double * vel, double x, double y,
                                         double z, ymir_locidx_t nodeid,
                                         void *data)
{
  vel[0] = 0.0;
  vel[1] = 12.0 * M_PI * M_PI * sin (M_PI * y) * cos (M_PI * z);
  vel[2] = 12.0 * M_PI * M_PI * cos (M_PI * y) * sin (M_PI * z);
}

/*This flow field is from Worthen et al., 2014, PEPI, divergence-free*/
static void
slabs_test_poly1_vel_in_fn (double * vel, double x, double y,
                                        double z, ymir_locidx_t nodeid,
                                        void *data)
{
  vel[0] = 0.0;
  vel[1] =  y + y*y - 2.0*y*z + y*y*y - 3.0*y*z*z + y*y*z;
  vel[2] = -z - 2*y*z + z*z - 3.0*y*y*z + z*z*z - y*z*z;
}

static void
slabs_test_poly1_TIrot90_vel_out_fn (double * vel, double x, double y,
                                  double z, ymir_locidx_t nodeid,
                                  void *data)
{
  vel[0] = 0.0;
  vel[1] = - (18.0 + 48.0 * y + 18.0 * z);
  vel[2] = - (18.0 - 18.0 * y + 48.0 * z);
}

static void
slabs_test_poly1_TIrot90_viscexp_vel_out_fn (double * vel, double x, double y,
                                  double z, ymir_locidx_t nodeid,
                                  void *data)
{
  vel[0] = 0.0;
  vel[1] = - (20.0 + 60.0 * y + 20.0 * z + (exp(y) + exp(z)) * (-z - 6 * y - 1.0));
  vel[2] = - (20.0 - 20.0 * y + 60.0 * z + (exp(y) + exp(z)) * ( y - 6 * z - 1.0));
}

static void
slabs_test_sincos1_TIrot60_viscexp60_vel_out_fn (double * vel, double x, double y,
                                                double z, ymir_locidx_t nodeid,
                                                void *data)
{
  double visc = exp(0.5 * (sqrt(3.0) * y + z));
  double a = (2.5 + 1.5 * visc) * M_PI * M_PI;
  double b = 0.5 * sqrt(3.0) * (visc - 5.0) * M_PI * M_PI;

  vel[0] = 0.0;
  vel[1] = a * sin(M_PI * y) * cos(M_PI * z) + b * cos(M_PI * y) * sin(M_PI * z);
  vel[2] = b * sin(M_PI * y) * cos(M_PI * z) - a * cos(M_PI * y) * sin(M_PI * z);
}

static void
slabs_test_poly1_manufactured_set_velbc (double * vel, double x, double y,
                                              double z, ymir_locidx_t nodeid,
                                              void *data)
{
  slabs_options_t  *slabs_options = data;
  const double  z_max = slabs_options->slabs_domain_options->z_max;
  const double  y_max = slabs_options->slabs_domain_options->y_max;
  const double  x_max = slabs_options->slabs_domain_options->x_max;

  if (y < SC_1000_EPS || (y_max - y) < SC_1000_EPS ||
      z < SC_1000_EPS || (z_max - z) < SC_1000_EPS ||
      x < SC_1000_EPS || (x_max - x) < SC_1000_EPS) {
    vel[0] = 0.0;
    vel[1] =  y + y*y - 2.0*y*z + y*y*y - 3.0*y*z*z + y*y*z;
    vel[2] = -z - 2*y*z + z*z - 3.0*y*y*z + z*z*z - y*z*z;
  }
  else {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }
}

/* manufactured solution */
static void
slabs_test_manufactured_rhs_compute (ymir_vec_t *rhs_vel,
                                     ymir_vec_t *temperature /* unused */,
                                     void *data)
{
  const char          *this_fn_name = "slabs_test_manufactured_rhs_compute";
  slabs_options_t     *slabs_options = data;
  const               slabs_test_manufactured_t
                      test_type = slabs_options->slabs_test_options->test_manufactured;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);
  /* check input */
  RHEA_ASSERT (test_type != SLABS_TEST_MANUFACTURED_NONE);

  /* compute velocity fields */
  switch (test_type) {
    case SLABS_TEST_MANUFACTURED_SINCOS1_ISO:
      /* compute reference velocity field (output) */
      ymir_cvec_set_function (rhs_vel, slabs_test_sincos1_ISO_vel_out_fn,
                              NULL);
      break;

    case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT90:
      /* compute reference velocity field (output) */
      ymir_cvec_set_function (rhs_vel, slabs_test_sincos1_TIrot90_vel_out_fn,
                              NULL);
      break;

    case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT45:
      /* compute reference velocity field (output) */
      ymir_cvec_set_function (rhs_vel, slabs_test_sincos1_TIrot45_vel_out_fn,
                              NULL);
      break;

    case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT60:
      /* compute reference velocity field (output) */
      ymir_cvec_set_function (rhs_vel, slabs_test_sincos1_TIrot60_vel_out_fn,
                              NULL);
      break;

    case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT60_VISCEXP60:
      /* compute reference velocity field (output) */
      ymir_cvec_set_function (rhs_vel, slabs_test_sincos1_TIrot60_viscexp60_vel_out_fn,
                              NULL);
      break;


    case SLABS_TEST_MANUFACTURED_POLY1_TIROT90:
      /* compute reference velocity field (output) */
      ymir_cvec_set_function (rhs_vel, slabs_test_poly1_TIrot90_vel_out_fn,
                              NULL);
      break;

    case SLABS_TEST_MANUFACTURED_POLY1_TIROT90_VISCEXP:
      /* compute reference velocity field (output) */
      ymir_cvec_set_function (rhs_vel, slabs_test_poly1_TIrot90_viscexp_vel_out_fn,
                              NULL);
      break;

  default:
      RHEA_ABORT_NOT_REACHED ();
  }
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

void
slabs_test_manufactured_velbc_compute (ymir_vec_t * rhs_vel_nonzero_dirichlet,
                                       void * data)
{
  const char          *this_fn_name = "slabs_test_manufactured_velbc_compute";
  slabs_options_t     *slabs_options = data;
  const               slabs_test_manufactured_t
                      test_type = slabs_options->slabs_test_options->test_manufactured;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);
  switch (test_type) {
    case SLABS_TEST_MANUFACTURED_SINCOS1_ISO:
    case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT90:
    case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT45:
    case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT60:
    case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT60_VISCEXP60:
      rhea_domain_set_user_velocity_dirichlet_bc (slabs_set_vel_dir_all, NULL, 0);
      ymir_cvec_set_function (rhs_vel_nonzero_dirichlet,
                              slabs_test_sincos1_manufactured_set_velbc,
                              slabs_options);
      break;

    case SLABS_TEST_MANUFACTURED_POLY1_TIROT90:
    case SLABS_TEST_MANUFACTURED_POLY1_TIROT90_VISCEXP:
      rhea_domain_set_user_velocity_dirichlet_bc (slabs_set_vel_dir_all, NULL, 0);
      ymir_cvec_set_function (rhs_vel_nonzero_dirichlet,
                              slabs_test_poly1_manufactured_set_velbc,
                              slabs_options);
      break;

   default: /* BC not set */
      RHEA_ABORT_NOT_REACHED ();
  }
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

static void
slabs_test_manufactured_compute_vel_err (double * abs_error, double * rel_error,
                                ymir_vec_t *vel_error, ymir_vec_t *vel_ref,
                                ymir_vec_t *vel_chk,
                                ymir_stress_op_t * stress_op)
{
  ymir_vec_t         *vel_ref_zero_bndr = ymir_vec_clone (vel_ref);

  /* calculate error vector */
  ymir_vec_copy (vel_chk, vel_error);
  ymir_vec_add (-1.0, vel_ref, vel_error);

  /* set boundary values to zero */
  if (ymir_stress_op_has_dirichlet (stress_op)) {
    ymir_vel_dir_separate (vel_error, NULL, NULL, NULL, stress_op->vel_dir);
    ymir_vel_dir_separate (vel_ref_zero_bndr, NULL, NULL, NULL,
                           stress_op->vel_dir);
  }

  /* take norm of error vector to get a value for absolute and relative errors */
  *abs_error = ymir_vec_norm (vel_error);
  *rel_error = *abs_error / ymir_vec_norm (vel_ref_zero_bndr);

  /* destroy */
  ymir_vec_destroy (vel_ref_zero_bndr);
}


/***********************************************************
 * Post-processing for 2nd invariant stress, traction, .etc.
 ***********************************************************/

/* compute traction as well as normal/shear stress at each element*/
void
slabs_stress_elem (sc_dmatrix_t * in, sc_dmatrix_t * out,
                   sc_dmatrix_t * visc, ymir_velocity_elem_t * vel_elem,
                   double *_sc_restrict rxd, double *_sc_restrict sxd,
                   double *_sc_restrict txd, double *_sc_restrict ryd,
                   double *_sc_restrict syd, double *_sc_restrict tyd,
                   double *_sc_restrict rzd, double *_sc_restrict szd,
                   double *_sc_restrict tzd, sc_dmatrix_t * drst, sc_dmatrix_t * brst)
{
  const int           N = ymir_n (vel_elem->N);
  const int           Np = ymir_n (vel_elem->Np);
  int                 gp, i, j, k, l;
  double              S[9];
  double             *_sc_restrict viscd  = visc->e[0];

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
    double               *_sc_restrict temptensd = temptens->e[0] + 9 * gp;
    double               *_sc_restrict outd  = out->e[0] + 9 * gp;

    S[0] = temptensd[0];
    S[1] = (temptensd[1] + temptensd[3]) * (1. / 2.);
    S[2] = (temptensd[2] + temptensd[6]) * (1. / 2.);
    S[3] = S[1];
    S[4] = temptensd[4];
    S[5] = (temptensd[5] + temptensd[7]) * (1. / 2.);
    S[6] = S[2];
    S[7] = S[5];
    S[8] = temptensd[8];

    /* compute linear combination of viscous stress and 3x3 tensor */
    for (i = 0; i < 9; i++)  {
      outd[i] = 2.0 * viscd[gp] * S[i];
    }
  }
}

/* compute traction as well as normal and shear stress*/
void
slabs_stress (ymir_cvec_t * vel, ymir_dvec_t * tau,
              ymir_dvec_t * visc, ymir_velocity_elem_t * vel_elem)
{
  ymir_mesh_t         *mesh = vel->mesh;
  ymir_locidx_t       elid;
  const int           N  = ymir_n (mesh->cnodes->N);
  const int           Np = ymir_np (mesh->cnodes->N);
  const int           K  = mesh->cnodes->K;
  sc_dmatrix_t       *elemin = sc_dmatrix_new (1, 3 * Np);
  sc_dmatrix_t       *elemout = sc_dmatrix_new (1,9 * Np);
  sc_dmatrix_t       *elemvisc = sc_dmatrix_new (1, Np);
  const char          *this_fn_name = "slabs_stress";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  ymir_dvec_set_zero (tau);

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
    ymir_dvec_get_elem_interp (tau, elemout, YMIR_STRIDE_NODE, elid,
                             YMIR_GAUSS_NODE, YMIR_WRITE);

    slabs_stress_elem (elemin, elemout,
                       elemvisc,
                       vel_elem, rxd, sxd, txd, ryd,
                       syd, tyd, rzd, szd, tzd, mesh->drst, mesh->brst);

    ymir_dvec_set_elem_interp (tau, elemout, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);

    ymir_read_view_release (elemvisc);
  }

  sc_dmatrix_destroy (elemin);
  sc_dmatrix_destroy (elemout);
  sc_dmatrix_destroy (elemvisc);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/* compute traction as well as normal/shear stress at each element*/
void
slabs_stress_TI_elem (sc_dmatrix_t * in, sc_dmatrix_t * out,
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
    double               *_sc_restrict temptensd = temptens->e[0] + 9 * gp;
    double               *_sc_restrict T   = TItens->e[0] + 9 * gp;
    double               *_sc_restrict outd   = out->e[0] + 9 * gp;

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
    for (i = 0; i < 9; i++)  {
      outd[i] = cs * S[i] + ct * E[i];
    }
  }
}

/* compute traction as well as normal and shear stress*/
void
slabs_stress_TI (ymir_cvec_t * vel, ymir_dvec_t * tau,
              ymir_dvec_t * visc, ymir_dvec_t *svisc,
              ymir_dvec_t * TItens, ymir_velocity_elem_t * vel_elem)
{
  ymir_mesh_t         *mesh = vel->mesh;
  ymir_locidx_t       elid;
  const int           N  = ymir_n (mesh->cnodes->N);
  const int           Np = ymir_np (mesh->cnodes->N);
  const int           K  = mesh->cnodes->K;
  sc_dmatrix_t       *elemin = sc_dmatrix_new (1, 3 * Np);
  sc_dmatrix_t       *elemout = sc_dmatrix_new (1,9 * Np);
  sc_dmatrix_t       *elemvisc = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemsvisc = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemTItens = sc_dmatrix_new (1, 9 * Np);
  const char          *this_fn_name = "slabs_stress_TI";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  ymir_dvec_set_zero (tau);

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
    ymir_dvec_get_elem_interp (TItens, elemTItens, YMIR_STRIDE_NODE, elid,
                              YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (tau, elemout, YMIR_STRIDE_NODE, elid,
                             YMIR_GAUSS_NODE, YMIR_WRITE);

    slabs_stress_TI_elem (elemin, elemout,
                       elemvisc, elemsvisc, elemTItens,
                       vel_elem, rxd, sxd, txd, ryd,
                       syd, tyd, rzd, szd, tzd, mesh->drst, mesh->brst);

    ymir_dvec_set_elem_interp (tau, elemout, YMIR_STRIDE_NODE, elid,
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

/* compute traction as well as normal/shear stress at each element*/
void
slabs_stressvec_TI_elem (sc_dmatrix_t * in, sc_dmatrix_t * out,
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
    double               *_sc_restrict temptensd = temptens->e[0] + 9 * gp;
    double               *_sc_restrict T   = TItens->e[0] + 9 * gp;
    double               *_sc_restrict outd   = out->e[0] + 3 * gp;

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
    outd[0] = cs * S[4] + ct * E[4];
    outd[1] = cs * S[8] + ct * E[8];
    outd[2] = cs * S[5] + ct * E[5];
  }
}

/* compute traction as well as normal and shear stress*/
void
slabs_stressvec_TI (ymir_cvec_t * vel, ymir_dvec_t * tauvec,
              ymir_dvec_t * visc, ymir_dvec_t *svisc,
              ymir_dvec_t * TItens, ymir_velocity_elem_t * vel_elem)
{
  ymir_mesh_t         *mesh = vel->mesh;
  ymir_locidx_t       elid;
  const int           N  = ymir_n (mesh->cnodes->N);
  const int           Np = ymir_np (mesh->cnodes->N);
  const int           K  = mesh->cnodes->K;
  sc_dmatrix_t       *elemin = sc_dmatrix_new (1, 3 * Np);
  sc_dmatrix_t       *elemout = sc_dmatrix_new (1,3 * Np);
  sc_dmatrix_t       *elemvisc = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemsvisc = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemTItens = sc_dmatrix_new (1, 9 * Np);
  const char          *this_fn_name = "slabs_stress_TI";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  ymir_dvec_set_zero (tauvec);

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
    ymir_dvec_get_elem_interp (TItens, elemTItens, YMIR_STRIDE_NODE, elid,
                              YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (tauvec, elemout, YMIR_STRIDE_NODE, elid,
                             YMIR_GAUSS_NODE, YMIR_WRITE);

    slabs_stressvec_TI_elem (elemin, elemout,
                       elemvisc, elemsvisc, elemTItens,
                       vel_elem, rxd, sxd, txd, ryd,
                       syd, tyd, rzd, szd, tzd, mesh->drst, mesh->brst);

    ymir_dvec_set_elem_interp (tauvec, elemout, YMIR_STRIDE_NODE, elid,
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
    for (i = 0; i < 9; i++)  {
      outv += SC_SQR (cs * S[i] + ct * E[i]) ;
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
  const char          *this_fn_name = "slabs_2inv_stress_TI";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

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
    ymir_dvec_get_elem_interp (TItens, elemTItens, YMIR_STRIDE_NODE, elid,
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

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/* Computes the shear and normal traction along 2plates_poly2 weakzone in 2D Cartesian domain*/
void
slabs_postp_weakzone_coupling_brick_2plates_poly2 (double r, double lon,
                                                   double *shear, double *normal,
                                                   double *tau,
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

  double              consider_thickness;
  double              y_min = slabs_options->slabs_domain_options->y_min;
  double              start_node, start_val, start_deriv, end_node, end_val;
  double              *poly2_coeff;
  double              orth[3], tangent[3], tempvec[3];
  int                 orientation_wrt_curve;
  double             *closest_pt, dist=1.0, s;


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
  RHEA_ASSERT (0.0 < subdu_lon);
  RHEA_ASSERT (0.0 < subdu_dip_angle);
  RHEA_ASSERT (0.0 < subdu_depth && 0.0 < subdu_width);
  RHEA_ASSERT (0.0 < subdu_thickness);
  RHEA_ASSERT (subdu_thickness_const <= subdu_thickness);
  RHEA_ASSERT (0.0 < ridge_depth && 0.0 < ridge_width);
  RHEA_ASSERT (0.0 <= ridge_smoothwidth);

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
  consider_thickness = (2.0 * subdu_thickness) / SLABS_MANTLE_DEPTH;
  * shear = * normal = 0.0;
  orth[0] = tangent[0] = 0.0;

  if (   (  start_node - 0.5 * consider_thickness
          / sin (subdu_dip_angle / 180.0 * M_PI) ) <= lon
      && lon <= (end_node + 0.5 * consider_thickness )
      && (end_val - consider_thickness ) <= r ) {
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

    if (orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_BACK)
      closest_pt[0] = start_node;
    else if (orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_LEFT ||
             orientation_wrt_curve == SLABS_CURVE_ORIENT_TOP_RIGHT)
      closest_pt[0] = end_node;

    if (r > 0.95) {
      /*from r=0.9 to r=1, shrink the width of weakzone by a factor of 2*/
      s = - 10.0 * r + 10.5;
    }
    else
      s = 1.0;

    if (dist < 0.5 * s * subdu_thickness)  {
      orth[1] = - slabs_poly2_deriv (closest_pt[0], poly2_coeff);
      orth[2] = 1.0 / sqrt(1.0 + orth[1] * orth[1]);
      orth[1] *= orth[2];
      tangent[2] = slabs_poly2_deriv (closest_pt[0], poly2_coeff);
      tangent[1] = 1.0 / sqrt (1.0 + tangent[2] * tangent[2]);
      tangent[2] *= tangent[1];

      tempvec[0] = ( tau[0] * orth[0] + tau[1] * orth[1] + tau[2] * orth[2]);
      tempvec[1] = ( tau[3] * orth[0] + tau[4] * orth[1] + tau[5] * orth[2]);
      tempvec[2] = ( tau[6] * orth[0] + tau[7] * orth[1] + tau[8] * orth[2]);

      * shear  = (tempvec[0] * tangent[0] + tempvec[1] * tangent[1] + tempvec[2] * tangent[2]);
      * normal = (tempvec[0] * orth[0] + tempvec[1] * orth[1] + tempvec[2] * orth[2]);
    }
    RHEA_FREE (closest_pt);
  }
}

void
slabs_postp_weakzone_coupling_node (const double x, const double y, const double z,
                                    double *shear, double *normal,
                                    double *tau,
                                    slabs_options_t *slabs_options)
{
  double              lon;
  double              r;

  /* compute radius and longitude */
  lon = y;
  r = z;

  /* compute weak zone factor */
  slabs_postp_weakzone_coupling_brick_2plates_poly2 (r, lon, shear, normal,
                                                     tau, slabs_options);
}

/* compute weak zone factor of an element */
void
slabs_postp_weakzone_coupling_elem (sc_dmatrix_t * tau,
                                    sc_dmatrix_t * shear,
                                    sc_dmatrix_t * normal,
                                    const double *x,
                                    const double *y,
                                    const double *z,
                                    const int n_nodes_per_el,
                                    slabs_options_t *slabs_options)
{
  int            nodeid;

  /* compute weak zone factor for each node */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    double               *_sc_restrict tau_data = tau->e[0] + 9 * nodeid;
    double               *_sc_restrict shear_data = shear->e[0] + nodeid;
    double               *_sc_restrict normal_data = normal->e[0] + nodeid;
    slabs_postp_weakzone_coupling_node (x[nodeid], y[nodeid], z[nodeid],
                                        shear_data, normal_data, tau_data,
                                        slabs_options);
  }
}

static void
slabs_postp_weakzone_coupling_compute (ymir_cvec_t *vel, ymir_dvec_t *visc,
                                       ymir_dvec_t *svisc, ymir_dvec_t *TItens,
                                       ymir_dvec_t *normal, ymir_dvec_t *shear,
                                       ymir_velocity_elem_t *vel_elem,
                                       slabs_options_t *slabs_options)
{
  ymir_mesh_t         *mesh = ymir_vec_get_mesh (vel);
  const ymir_locidx_t n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  ymir_dvec_t         *tau = ymir_dvec_new (mesh, 9, YMIR_GAUSS_NODE);
  sc_dmatrix_t        *tau_el_mat, *shear_el_mat, *normal_el_mat;
  sc_dmatrix_t        *elemtau, *elemout_s, *elemout_n;
  double              *tau_el_data, *shear_el_data, *normal_el_data;
  double              *x, *y, *z, *tmp_el;
  ymir_locidx_t       elid;

  double              start_node, start_val, start_deriv;
  double              end_node, end_val;
  double              subdu_lon;
  double              subdu_dip_angle;
  double              subdu_depth, subdu_width;
  double              *poly2_coeff;
  slabs_2plates_poly2_geo_coeff_t geo;
  slabs_weak_options_t *weak_options = slabs_options->slabs_weak_options;
  const char          *this_fn_name = "slabs_postp_weakzone_coupling_compute";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  if (slabs_options->slabs_visc_options->viscosity_anisotropy
      == SLABS_VISC_TRANSVERSELY_ISOTROPY)  {
    slabs_stress_TI (vel, tau, visc, svisc, TItens, vel_elem);
  }
  else  {
    slabs_stress (vel, tau, visc, vel_elem);
  }

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

  /* create work variables */
  elemtau = sc_dmatrix_new (1, 9 * n_nodes_per_el);
  elemout_s = sc_dmatrix_new (1, n_nodes_per_el);
  elemout_n = sc_dmatrix_new (1, n_nodes_per_el);
  x = RHEA_ALLOC (double, n_nodes_per_el);
  y = RHEA_ALLOC (double, n_nodes_per_el);
  z = RHEA_ALLOC (double, n_nodes_per_el);
  tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    /* get coordinates of this element at Gauss nodes */
    ymir_mesh_get_elem_coord_gauss (x, y, z, elid, mesh, tmp_el);

    ymir_dvec_get_elem_interp (tau, elemtau, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (shear, elemout_s, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_WRITE);
    ymir_dvec_get_elem_interp (normal, elemout_n, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_WRITE);

    /* compute traction on each element */
    slabs_postp_weakzone_coupling_elem (elemtau, elemout_s, elemout_n,
                                        x, y, z, n_nodes_per_el, slabs_options);

    /* set traction of this element */
    ymir_dvec_set_elem_interp (shear, elemout_s, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);
    ymir_dvec_set_elem_interp (normal, elemout_n, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);
  }

  /* destroy */
  sc_dmatrix_destroy (elemtau);
  sc_dmatrix_destroy (elemout_s);
  sc_dmatrix_destroy (elemout_n);
  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);
  RHEA_FREE (poly2_coeff);

  ymir_vec_destroy (tau);
}


/* Computes the shear and normal traction along 2plates_poly2 weakzone in 2D Cartesian domain*/
void
slabs_postp_weakzone_traction_brick_2plates_poly2 (double r, double lon,
                                                   double *trac, double *gradv,
                                                   double *svisc, double *nvisc,
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

  double              consider_thickness;
  double              y_min = slabs_options->slabs_domain_options->y_min;
  double              start_node, start_val, start_deriv, end_node, end_val;
  double              *poly2_coeff;
  double              orth[3], tangent[3], tempvec[3];
  int                 orientation_wrt_curve;
  double             *closest_pt, dist=1.0, s;
  double              tempt, tt;


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
  RHEA_ASSERT (0.0 < subdu_lon);
  RHEA_ASSERT (0.0 < subdu_dip_angle);
  RHEA_ASSERT (0.0 < subdu_depth && 0.0 < subdu_width);
  RHEA_ASSERT (0.0 < subdu_thickness);
  RHEA_ASSERT (subdu_thickness_const <= subdu_thickness);
  RHEA_ASSERT (0.0 < ridge_depth && 0.0 < ridge_width);
  RHEA_ASSERT (0.0 <= ridge_smoothwidth);

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
  consider_thickness = (2.0 * subdu_thickness) / SLABS_MANTLE_DEPTH;

  trac[0] = 0.0;
  trac[1] = 0.0;
  trac[2] = 0.0;
  orth[0] = tangent[0] = 0.0;

  if (   (  start_node - 0.5 * consider_thickness
          / sin (subdu_dip_angle / 180.0 * M_PI) ) <= lon
      && lon <= (end_node + 0.5 * consider_thickness )
      && (end_val - consider_thickness ) <= r ) {
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
    if (orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_BACK)
      closest_pt[0] = start_node;
    else if (orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_LEFT ||
             orientation_wrt_curve == SLABS_CURVE_ORIENT_TOP_RIGHT)
      closest_pt[0] = end_node;

    if (r > 0.95) {
      /*from r=0.9 to r=1, shrink the width of weakzone by a factor of 2*/
      s = - 10.0 * r + 10.5;
    }
    else
      s = 1.0;

    if (dist < 0.5 * s * subdu_thickness)  {
//      orth[1] = - slabs_poly2_deriv (closest_pt[0], poly2_coeff);
//      orth[2] = 1.0 / sqrt(1.0 + orth[1] * orth[1]);
//      orth[1] *= orth[2];
//      tangent[2] = slabs_poly2_deriv (closest_pt[0], poly2_coeff);
//      tangent[1] = 1.0 / sqrt (1.0 + tangent[2] * tangent[2]);
//      tangent[2] *= tangent[1];

      tempt = - slabs_poly2_deriv (closest_pt[0], poly2_coeff);
      tt = atan(tempt);
      orth[1] = sin(tt);
      orth[2] = cos(tt);
      tangent[1] = cos(tt);
      tangent[2] = -sin(tt);

      tempvec[0] = ( gradv[0] * orth[0] + gradv[1] * orth[1] + gradv[2] * orth[2]);
      tempvec[1] = ( gradv[1] * orth[0] + gradv[3] * orth[1] + gradv[4] * orth[2]);
      tempvec[2] = ( gradv[2] * orth[0] + gradv[4] * orth[1] + gradv[5] * orth[2]);

      trac[1] = 2.0 * svisc[0] *
                (tempvec[0] * tangent[0] + tempvec[1] * tangent[1] + tempvec[2] * tangent[2]);
      trac[2] = 2.0 * nvisc[0] *
                (tempvec[0] * orth[0] + tempvec[1] * orth[1] + tempvec[2] * orth[2]);

    }
    RHEA_FREE (closest_pt);
  }
}

void
slabs_postp_weakzone_traction_node (const double x, const double y, const double z,
                                    double *trac, double *gradv,
                                    double *svisc, double *nvisc,
                                    slabs_options_t *slabs_options)
{
  double              lon;
  double              r;

  /* compute radius and longitude */
  lon = y;
  r = z;

  /* compute weak zone factor */
  slabs_postp_weakzone_traction_brick_2plates_poly2 (r, lon, trac, gradv,
                                                     svisc, nvisc, slabs_options);
}

/* compute weak zone factor of an element */
void
slabs_postp_weakzone_traction_elem (sc_dmatrix_t *trac,
                                    sc_dmatrix_t *gradv,
                                    sc_dmatrix_t *svisc,
                                    sc_dmatrix_t *nvisc,
                                    const double *x,
                                    const double *y,
                                    const double *z,
                                    const int n_nodes_per_el,
                                    slabs_options_t *slabs_options)
{
  int            nodeid;

  /* compute weak zone factor for each node */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    double               *_sc_restrict trac_data = trac->e[0] + 3 * nodeid;
    double               *_sc_restrict gradv_data = gradv->e[0] + 6 * nodeid;
    double               *_sc_restrict svisc_data = svisc->e[0] + nodeid;
    double               *_sc_restrict nvisc_data = nvisc->e[0] + nodeid;
    slabs_postp_weakzone_traction_node (x[nodeid], y[nodeid], z[nodeid],
                                        trac_data, gradv_data,
                                        svisc_data, nvisc_data,
                                        slabs_options);
  }
}

static void
slabs_postp_weakzone_traction_compute (ymir_vec_t *vel, ymir_vec_t *traction,
                                       ymir_vec_t *svisc, ymir_vec_t *nvisc,
                                       slabs_options_t *slabs_options)
{
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (vel);
  const ymir_locidx_t n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  ymir_dvec_t        *gradvel = ymir_dvec_new (mesh, 6, YMIR_GAUSS_NODE);
  sc_dmatrix_t       *elemtrac, *elemgradv, *elemsvisc, *elemnvisc;
  double             *x, *y, *z, *tmp_el;
  ymir_locidx_t       elid;

  double              start_node, start_val, start_deriv;
  double              end_node, end_val;
  double              subdu_lon;
  double              subdu_dip_angle;
  double              subdu_depth, subdu_width;
  double              *poly2_coeff;
  slabs_2plates_poly2_geo_coeff_t geo;
  slabs_weak_options_t *weak_options = slabs_options->slabs_weak_options;
  const char         *this_fn_name = "slabs_postp_weakzone_traction_compute";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  ymir_velocity_strain_rate (vel, gradvel, 0);

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

  /* create work variables */
  elemtrac = sc_dmatrix_new (1, 3 * n_nodes_per_el);
  elemgradv = sc_dmatrix_new (1, 6 * n_nodes_per_el);
  elemsvisc = sc_dmatrix_new (1, n_nodes_per_el);
  elemnvisc = sc_dmatrix_new (1, n_nodes_per_el);
  x = RHEA_ALLOC (double, n_nodes_per_el);
  y = RHEA_ALLOC (double, n_nodes_per_el);
  z = RHEA_ALLOC (double, n_nodes_per_el);
  tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    /* get coordinates of this element at Gauss nodes */
    ymir_mesh_get_elem_coord_gauss (x, y, z, elid, mesh, tmp_el);

    ymir_dvec_get_elem_interp (gradvel, elemgradv, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (svisc, elemsvisc, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (nvisc, elemnvisc, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_READ);

    ymir_dvec_get_elem_interp (traction, elemtrac, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_WRITE);
    /* compute traction on each element */
    slabs_postp_weakzone_traction_elem (elemtrac, elemgradv, elemsvisc, elemnvisc,
                                        x, y, z, n_nodes_per_el, slabs_options);

    /* set traction of this element */
    ymir_dvec_set_elem_interp (traction, elemtrac, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);

    ymir_read_view_release (elemgradv);
    ymir_read_view_release (elemsvisc);
    ymir_read_view_release (elemnvisc);
  }

  /* destroy */
  sc_dmatrix_destroy (elemtrac);
  sc_dmatrix_destroy (elemgradv);
  sc_dmatrix_destroy (elemsvisc);
  sc_dmatrix_destroy (elemnvisc);
  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);
  RHEA_FREE (poly2_coeff);

  ymir_vec_destroy (gradvel);
}

void
slabs_manufactured_stressvec_coupling_node (double *shear, double *normal,
                                            double *tau)
{
  double              orth[3], tangent[3], tempvec[3];
  double              tempt, tt;

  orth[0] = tangent[0] = 0.0;

   tt = 60.0 * M_PI / 180.0;
   orth[1] = sin(tt);
   orth[2] = cos(tt);
   tangent[1] = cos(tt);
   tangent[2] = -sin(tt);

   tempvec[0] = ( tau[0] * orth[0] + tau[1] * orth[1] + tau[2] * orth[2]);
   tempvec[1] = ( tau[3] * orth[0] + tau[4] * orth[1] + tau[5] * orth[2]);
   tempvec[2] = ( tau[6] * orth[0] + tau[7] * orth[1] + tau[8] * orth[2]);

   * shear  = (tempvec[0] * tangent[0] + tempvec[1] * tangent[1] + tempvec[2] * tangent[2]);
   * normal = (tempvec[0] * orth[0] + tempvec[1] * orth[1] + tempvec[2] * orth[2]);
}

static void
slabs_manufactured_stressvec_coupling_compute (ymir_cvec_t *vel, ymir_dvec_t *visc,
                                       ymir_dvec_t *svisc, ymir_dvec_t *TItens,
                                       ymir_dvec_t *normal, ymir_dvec_t *shear,
                                       ymir_velocity_elem_t *vel_elem,
                                       slabs_options_t *slabs_options)
{
  ymir_mesh_t         *mesh = ymir_vec_get_mesh (vel);
  const ymir_locidx_t n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  ymir_dvec_t         *tau = ymir_dvec_new (mesh, 9, YMIR_GAUSS_NODE);
  sc_dmatrix_t        *tau_el_mat, *shear_el_mat, *normal_el_mat;
  sc_dmatrix_t        *elemtau, *elemout_s, *elemout_n;
  double              *tau_el_data, *shear_el_data, *normal_el_data;
  double              *x, *y, *z, *tmp_el;
  ymir_locidx_t       elid;
  int                 nodeid;

  const char          *this_fn_name = "slabs_manufactured_stressvec_coupling_compute";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  if (slabs_options->slabs_visc_options->viscosity_anisotropy
      == SLABS_VISC_TRANSVERSELY_ISOTROPY)  {
    slabs_stress_TI (vel, tau, visc, svisc, TItens, vel_elem);
  }
  else  {
    slabs_stress (vel, tau, visc, vel_elem);
  }

 /* create work variables */
  elemtau = sc_dmatrix_new (1, 9 * n_nodes_per_el);
  elemout_s = sc_dmatrix_new (1, n_nodes_per_el);
  elemout_n = sc_dmatrix_new (1, n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    ymir_dvec_get_elem_interp (tau, elemtau, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (shear, elemout_s, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_WRITE);
    ymir_dvec_get_elem_interp (normal, elemout_n, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_WRITE);

    /* compute weak zone factor for each node */
    for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
      double               *_sc_restrict tau_data = elemtau->e[0] + 9 * nodeid;
      double               *_sc_restrict shear_data = elemout_s->e[0] + nodeid;
      double               *_sc_restrict normal_data = elemout_n->e[0] + nodeid;
      slabs_manufactured_stressvec_coupling_node (shear_data, normal_data, tau_data);
    }

    /* set traction of this element */
    ymir_dvec_set_elem_interp (shear, elemout_s, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);
    ymir_dvec_set_elem_interp (normal, elemout_n, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);
  }

  /* destroy */
  sc_dmatrix_destroy (elemtau);
  sc_dmatrix_destroy (elemout_s);
  sc_dmatrix_destroy (elemout_n);
  ymir_vec_destroy (tau);
}

/* Computes the shear and normal traction along 2plates_poly2 weakzone in 2D Cartesian domain*/
void
slabs_manufactured_strainvec_coupling_node (double *tracn, double *tracs, double *gradv,
                                           double *svisc, double *nvisc)
{
  double              orth[3], tangent[3], tempvec[3];
  double              tempt, tt;

  orth[0] = tangent[0] = 0.0;

  tt = 60.0 * M_PI / 180.0;
  orth[1] = sin(tt);
  orth[2] = cos(tt);
  tangent[1] = cos(tt);
  tangent[2] = -sin(tt);

  tempvec[0] = ( gradv[0] * orth[0] + gradv[1] * orth[1] + gradv[2] * orth[2]);
  tempvec[1] = ( gradv[1] * orth[0] + gradv[3] * orth[1] + gradv[4] * orth[2]);
  tempvec[2] = ( gradv[2] * orth[0] + gradv[4] * orth[1] + gradv[5] * orth[2]);

  tracn[0] = 2.0 * nvisc[0] *
            (tempvec[0] * orth[0] + tempvec[1] * orth[1] + tempvec[2] * orth[2]);
  tracs[0] = 2.0 * svisc[0] *
            (tempvec[0] * tangent[0] + tempvec[1] * tangent[1] + tempvec[2] * tangent[2]);

}

static void
slabs_manufactured_strainvec_coupling_compute (ymir_vec_t *vel, ymir_vec_t *tracn, ymir_vec_t *tracs,
                                       ymir_vec_t *svisc, ymir_vec_t *nvisc)
{
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (vel);
  const ymir_locidx_t n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  ymir_dvec_t        *gradvel = ymir_dvec_new (mesh, 6, YMIR_GAUSS_NODE);
  sc_dmatrix_t       *elemtracn, *elemtracs, *elemgradv, *elemsvisc, *elemnvisc;
  ymir_locidx_t       elid;
  int                 nodeid;

 const char         *this_fn_name = "slabs_manufactured_strainvec_coupling_compute";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  ymir_velocity_strain_rate (vel, gradvel, 0);

  /* create work variables */
  elemtracn = sc_dmatrix_new (1, n_nodes_per_el);
  elemtracs = sc_dmatrix_new (1, n_nodes_per_el);
  elemgradv = sc_dmatrix_new (1, 6 * n_nodes_per_el);
  elemsvisc = sc_dmatrix_new (1, n_nodes_per_el);
  elemnvisc = sc_dmatrix_new (1, n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    ymir_dvec_get_elem_interp (gradvel, elemgradv, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (svisc, elemsvisc, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (nvisc, elemnvisc, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_READ);

    ymir_dvec_get_elem_interp (tracn, elemtracn, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_WRITE);
    ymir_dvec_get_elem_interp (tracs, elemtracs, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_WRITE);

  /* compute weak zone factor for each node */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    double               *_sc_restrict tracn_data = elemtracn->e[0] + nodeid;
    double               *_sc_restrict tracs_data = elemtracs->e[0] + nodeid;
    double               *_sc_restrict gradv_data = elemgradv->e[0] + 6 * nodeid;
    double               *_sc_restrict svisc_data = elemsvisc->e[0] + nodeid;
    double               *_sc_restrict nvisc_data = elemnvisc->e[0] + nodeid;
    slabs_manufactured_strainvec_coupling_node (tracn_data, tracs_data, gradv_data,
                                                svisc_data, nvisc_data);
  }

    /* set traction of this element */
    ymir_dvec_set_elem_interp (tracn, elemtracn, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);
    ymir_dvec_set_elem_interp (tracs, elemtracs, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);

    ymir_read_view_release (elemgradv);
    ymir_read_view_release (elemsvisc);
    ymir_read_view_release (elemnvisc);
  }

  /* destroy */
  sc_dmatrix_destroy (elemtracn);
  sc_dmatrix_destroy (elemtracs);
  sc_dmatrix_destroy (elemgradv);
  sc_dmatrix_destroy (elemsvisc);
  sc_dmatrix_destroy (elemnvisc);
  ymir_vec_destroy (gradvel);
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
  ymir_vec_t         *rhs_vel, *rhs_vel_nonzero_dirichlet;
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

 // {
 //   rhs_vel_nonzero_dirichlet = rhea_stokes_problem_get_rhs_vel_nonzero_dirichlet (stokes_problem);
 //   if (rhs_vel_nonzero_dirichlet != NULL)
 //     snprintf (path, BUFSIZ, "%s_nonzero_dirichlet", vtk_write_input_path);
 //     ymir_vtk_write (ymir_mesh, path,
 //                     rhs_vel_nonzero_dirichlet, "velocity B.C.",
 //                     NULL);
 // }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

static int
slabs_output_pressure (const char *filepath, ymir_vec_t *pressure)
{
  const char *this_fn_name = "slabs_output_pressure";

  ymir_mesh_t        *mesh = ymir_vec_get_mesh (pressure);
  mangll_cnodes_t    *cnodes = mesh->cnodes;
  const int           N = cnodes->N;
  int                 Np;
  ymir_locidx_t       n_elements;
  ymir_locidx_t       elid;
  ymir_face_mesh_t   *fmesh;
  ymir_locidx_t       Ntotal;
  int                 i, j, k;
  double             *x, *y, *z, *tmp_el;
  ymir_topidx_t       fm = pressure->meshnum;
  int                 mpirank = mesh->ma->mpirank;
  sc_dmatrix_t       *elem;
  double             *elemd;
  FILE               *prsfile;
  char                prsfilename[BUFSIZ];


  if (fm == YMIR_VOL_MESH) {
    Np = (N + 1) * (N + 1) * (N + 1);
    n_elements = cnodes->K;
  }
  else {
    fmesh = &(mesh->fmeshes[fm]);
    Np = (N + 1) * (N + 1);
    n_elements = fmesh->K;
  }
  Ntotal = Np * n_elements;

  /* Have each proc write to its own file */
  if (fm == YMIR_VOL_MESH) {
    snprintf (prsfilename, BUFSIZ, "%s_pressure_%04d", filepath, mpirank);
  }
  else {
    snprintf (prsfilename, BUFSIZ, "%s_pressure_%04d.face%d", filepath, mpirank,
              (int) fm);
  }
  prsfile = fopen (prsfilename, "w");
  if (prsfile == NULL) {
    YMIR_LERRORF ("Could not open %s for output!\n", prsfilename);
    return -1;
  }
  x = RHEA_ALLOC (double, Np);
  y = RHEA_ALLOC (double, Np);
  z = RHEA_ALLOC (double, Np);
  tmp_el = RHEA_ALLOC (double, Np);

  elem = sc_dmatrix_new (Np, 1);
  for (elid = 0; elid < n_elements; elid++) {
    if (pressure->nefields)
      ymir_vec_get_elem_interp (pressure, elem, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_COPY);
    else if (pressure->ndfields)
      ymir_dvec_get_elem (pressure, elem, YMIR_STRIDE_NODE, elid,
                          YMIR_COPY);

    ymir_mesh_get_elem_coord_gauss (x, y, z, elid, mesh, tmp_el);

    elemd = elem->e[0];
    for (i = 0; i < Np; ++i)  {
      fprintf (prsfile, "%24.16e    %24.16e    %24.16e    %24.16e\n", x[i],
             y[i], z[i], elemd[i]);
    }
  }
  if (fclose (prsfile)) {
    YMIR_LERROR ("ymir_vtk: Error closing footer\n");
    return -1;
  }
  prsfile = NULL;

  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);
  sc_dmatrix_destroy (elem);

  return 0;
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

//  slabs_surface_location (slabs_options, discr_options);

  /* create p4est */
  *p4est = rhea_discretization_p4est_new (mpicomm, discr_options,
                                          domain_options);

  /* set up boundary, store in `discr_options` */
  rhea_discretization_boundary_create (discr_options, *p4est,
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
                    rhea_weakzone_options_t *weak_options,
                    rhea_viscosity_options_t *visc_options,
                    slabs_options_t *slabs_options,
                    const char *vtk_write_input_path)
{
  const char         *this_fn_name = "slabs_setup_stokes";
  ymir_vec_t         *temperature;
  ymir_vec_t         *coeff_TI_svisc = NULL, *TI_rotate = NULL;
  int                 mpirank = ymir_mesh->ma->mpirank;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* compute temperature */
  temperature = rhea_temperature_new (ymir_mesh);
  switch (slabs_options->buoyancy_type)  {
  /* if non-specified, use: rhea_temperature_compute (temperature, temp_options); */
    case SINKER:
      rhea_temperature_compute (temperature, temp_options);
      break;
    case SLAB:
      slabs_poly2_temperature_compute (temperature, slabs_options);
      break;
    case COLLIDE:
      rhea_temperature_compute (temperature, temp_options);
      break;
    case DRAG:
      ymir_cvec_set_function (temperature, drag_temperature_set_fn, slabs_options);
      break;
    case TEST_MANUFACTURED:
      RHEA_GLOBAL_PRODUCTIONF ("buoyancy from %i\n", slabs_options->buoyancy_type);
      rhea_temperature_compute (temperature, temp_options);  // neutral value: T=0.5
      break;
    case CUSTOM:
      slabs_custom_temperature_compute (temperature, slabs_options);
      break;
    case TESTNONE:
      ymir_vec_set_value (temperature, 0.5);
      break;
    default:
      RHEA_ABORT_NOT_REACHED ();
  }

  /* create Stokes problem */
  *stokes_problem = rhea_stokes_problem_new (
      ymir_mesh, press_elem, temperature, domain_options, temp_options,
      weak_options, visc_options);

  /* provide own function to compute weak zones and viscosity */
  switch (slabs_options->buoyancy_type) {
    case SINKER:
    case TESTNONE:
      break;
    case SLAB:
      rhea_stokes_problem_set_weakzone_compute_fn (
          *stokes_problem, slabs_weakzone_compute, slabs_options);
      break;
    case DRAG:
    case COLLIDE:
      if (slabs_options->slabs_visc_options->viscosity_anisotropy
        == SLABS_VISC_ISOTROPY) {
        rhea_stokes_problem_set_weakzone_compute_fn (
            *stokes_problem, slabs_weakzone_compute, slabs_options);
      }
      rhea_stokes_problem_set_viscosity_compute_fn (
          *stokes_problem, slabs_viscosity_compute, slabs_options);
      break;
    case TEST_MANUFACTURED:
      rhea_stokes_problem_set_viscosity_compute_fn (
          *stokes_problem, slabs_viscosity_compute, slabs_options);
      break;
    case CUSTOM:
      rhea_stokes_problem_set_viscosity_compute_fn (
          *stokes_problem, slabs_viscosity_compute, slabs_options);
      break;
    default:
      RHEA_ABORT_NOT_REACHED ();
  }

  /* for the test using manufactured solution,
   * overwrite rhs_vel with estimated forcing term from given velocity and pressure field.
   * don't apply mass now */
  if (slabs_options->buoyancy_type == TEST_MANUFACTURED ||
      slabs_options->slabs_test_options->test_manufactured != SLABS_TEST_MANUFACTURED_NONE) {
    rhea_stokes_problem_set_rhs_vel_compute_fn (
        *stokes_problem, slabs_test_manufactured_rhs_compute, slabs_options);
  }

  /* set velocity boundary conditions & nonzero Dirichlet values */
  if (domain_options->velocity_bc_type == RHEA_DOMAIN_VELOCITY_BC_USER) {
    if (slabs_options->buoyancy_type == TEST_MANUFACTURED ||
        slabs_options->slabs_test_options->test_manufactured != SLABS_TEST_MANUFACTURED_NONE) {
      rhea_stokes_problem_set_rhs_vel_nonzero_dir_compute_fn (
          *stokes_problem, slabs_test_manufactured_velbc_compute, slabs_options);
    }
    else {
      rhea_stokes_problem_set_rhs_vel_nonzero_dir_compute_fn (
          *stokes_problem, slabs_vel_nonzero_dirichlet_compute, slabs_options);
    }
  }

  /* add the anisotropic viscosity to the viscous stress operator */
  if (slabs_options->slabs_visc_options->viscosity_anisotropy
      == SLABS_VISC_TRANSVERSELY_ISOTROPY) {
    coeff_TI_svisc = rhea_viscosity_new (ymir_mesh);
    TI_rotate = rhea_viscosity_new (ymir_mesh);
    if (slabs_options->buoyancy_type == TEST_MANUFACTURED ||
        slabs_options->slabs_test_options->test_manufactured != SLABS_TEST_MANUFACTURED_NONE) {
      slabs_stokes_problem_setup_TI_manufactured (ymir_mesh, *stokes_problem, slabs_options,
                                     coeff_TI_svisc, TI_rotate);
    }
    else {
      slabs_stokes_problem_setup_TI (ymir_mesh, *stokes_problem, slabs_options,
                                     coeff_TI_svisc, TI_rotate);
    }
  }

  /* write vtk of problem input */
  if (vtk_write_input_path != NULL) {
    slabs_write_input (ymir_mesh, *stokes_problem, temp_options,
                       temperature, NULL, coeff_TI_svisc, TI_rotate,
                       vtk_write_input_path);
  }

  /* set up Stokes solver */
  rhea_stokes_problem_setup_solver (*stokes_problem);

if (mpirank == 0) {
printf("Got here!");
fflush (stdout);
}
  /* destroy */
  if (slabs_options->slabs_visc_options->viscosity_anisotropy
      == SLABS_VISC_TRANSVERSELY_ISOTROPY)  {
      rhea_viscosity_destroy (TI_rotate);
  }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

static void
slabs_setup_stokes_surfdist (rhea_stokes_problem_t **stokes_problem,
                    ymir_mesh_t *ymir_mesh,
                    ymir_pressure_elem_t *press_elem,
                    rhea_domain_options_t *domain_options,
                    rhea_temperature_options_t *temp_options,
                    rhea_viscosity_options_t *visc_options,
                    slabs_options_t *slabs_options,
                    const char *vtk_write_input_path,
                    ymir_vec_t *temperature, ymir_vec_t *weakzone)
{
  const char         *this_fn_name = "slabs_setup_stokes";
//  ymir_vec_t         *temperature, *weakzone;
  ymir_vec_t         *coeff_TI_svisc = NULL, *TI_rotate = NULL;
  ymir_vec_t         *rhs_vel, *rhs_vel_nonzero_dirichlet = NULL;
  int                 mpirank = ymir_mesh->ma->mpirank;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  rhs_vel = rhea_velocity_new (ymir_mesh);
  /* for the test using manufactured solution,
   * overwrite rhs_vel with estimated forcing term from given velocity and pressure field.
   * don't apply mass now */
  if (slabs_options->buoyancy_type == TEST_MANUFACTURED ||
      slabs_options->slabs_test_options->test_manufactured != SLABS_TEST_MANUFACTURED_NONE) {
    slabs_test_manufactured_rhs_compute (rhs_vel, NULL, slabs_options);
  }
  else {
    /* compute velocity right-hand side volume forcing */
    rhea_temperature_compute_rhs_vel (rhs_vel, temperature, temp_options);
  }

  /* set velocity boundary conditions & nonzero Dirichlet values */
  if (domain_options->velocity_bc_type == RHEA_DOMAIN_VELOCITY_BC_USER) {
    rhs_vel_nonzero_dirichlet = rhea_velocity_new (ymir_mesh);
    if (slabs_options->buoyancy_type == TEST_MANUFACTURED ||
        slabs_options->slabs_test_options->test_manufactured != SLABS_TEST_MANUFACTURED_NONE) {
      slabs_test_manufactured_velbc_compute (rhs_vel_nonzero_dirichlet,
                                             slabs_options);
    }
    else {
      slabs_vel_nonzero_dirichlet_compute (rhs_vel_nonzero_dirichlet,
                                           slabs_options);
    }
  }

  /* create Stokes problem */
  //TODO change function arguments as above
//*stokes_problem = rhea_stokes_problem_new (
//    temperature, weakzone, rhs_vel, rhs_vel_nonzero_dirichlet,
//    ymir_mesh, press_elem, domain_options, visc_options);

  /* add the anisotropic viscosity to the viscous stress operator */
  if (slabs_options->slabs_visc_options->viscosity_anisotropy
      == SLABS_VISC_TRANSVERSELY_ISOTROPY) {
    coeff_TI_svisc = rhea_viscosity_new (ymir_mesh);
    TI_rotate = rhea_viscosity_new (ymir_mesh);
    if (slabs_options->buoyancy_type == TEST_MANUFACTURED ||
        slabs_options->slabs_test_options->test_manufactured != SLABS_TEST_MANUFACTURED_NONE) {
      slabs_stokes_problem_setup_TI_manufactured (ymir_mesh, *stokes_problem, slabs_options,
                                     coeff_TI_svisc, TI_rotate);
    }
    else {
      slabs_stokes_problem_setup_TI (ymir_mesh, *stokes_problem, slabs_options,
                                     coeff_TI_svisc, TI_rotate);
    }
  }

  /* write vtk of problem input */
  if (vtk_write_input_path != NULL) {
    slabs_write_input (ymir_mesh, *stokes_problem, temp_options,
                       temperature, weakzone, coeff_TI_svisc, TI_rotate,
                       vtk_write_input_path);
  }

  /* set up Stokes solver */
  rhea_stokes_problem_setup_solver (*stokes_problem);

if (mpirank == 0) {
printf("Got here!");
fflush (stdout);
}
  /* destroy */
  if (slabs_options->slabs_visc_options->viscosity_anisotropy
      == SLABS_VISC_TRANSVERSELY_ISOTROPY)  {
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
  ymir_vec_t         *temperature;
  ymir_vec_t         *visc_TI_svisc;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* get vectors */
  temperature = rhea_stokes_problem_get_temperature (stokes_problem);

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

  /* destroy mesh */
  rhea_discretization_ymir_mesh_destroy (ymir_mesh, press_elem);
  rhea_discretization_p4est_destroy (p4est);

  /* destroy (some) options */
  rhea_discretization_boundary_clear (discr_options);

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
  rhea_stokes_problem_solve (&sol_vel_press, 0, iter_max, rel_tol,
                             stokes_problem);

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
  rhea_weakzone_options_t       weak_options;
  rhea_viscosity_options_t      visc_options;
  rhea_topography_options_t     topo_options;
  rhea_discretization_options_t discr_options;

  double                  *tX, *tY, *tZ;
  slabs_topo_profile_t    topo = {.tZ=NULL};
  double                  surf_dist;
  double                  visc_trial;
  int                     buoyancy_type;
  int                     temp_custom;
  int                     viscosity_anisotropy;
  int                     viscosity_geometry;
  int                     x_func;
  int                     vel_dir_bc;
  int                     test_manufactured;
  int                     test_stress_op;
  int                     test_stress_comp;

  /* slabs options */
  slabs_domain_options_t   slabs_domain_options;
  slabs_temp_options_t     slabs_temp_options;
  slabs_custom_sinker_t    slabs_sinker_options;
  slabs_visc_options_t     slabs_visc_options;
  slabs_weak_options_t     slabs_weak_options;
  slabs_surf_options_t     slabs_surf_options;
  slabs_velbc_options_t    slabs_velbc_options;
  slabs_test_options_t     slabs_test_options;
  slabs_options_t          slabs_options;

  /* options local to this function */
  int                 production_run;
  int                 solver_iter_max;
  double              solver_rel_tol;
  char               *vtk_write_input_path;
  char               *vtk_write_solution_path;
  char               *vtk_write_stress_path;
  char               *vtk_write_postp_path;
  char               *vtk_write_test_path;
  char               *vtk_write_freesurface_path;
  char               *vtk_write_input2_path;
  char               *vtk_write_solution2_path;
  char               *vtk_write_io_path;
  char               *vtk_write_mpiio_path;
  char               *vtk_read_io_path;
  char               *vtk_read_mpiio_path;
  char               *vtk_write_ioface_path;
  char               *vtk_read_ioface_path;
  char               *ascii_read_topo_path;
  char               *vtk_write_testdist_path;

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

  /* buoyancy */
  YMIR_OPTIONS_I, "buoyancy-type",'\0',
    &(buoyancy_type),SINKER,
    "0: sinker, 1: 2plates_poly2 slabs, 2: collide",

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

  YMIR_OPTIONS_I, "temp-custom-type",'\0',
    &(temp_custom),CUSTOM_SINKER,
    "0: thin_box, 1: sinker",
  YMIR_OPTIONS_D, "temp-sinker-center-x", '\0',
    &center_x, 0.0,
    "center of the sinker: x",
  YMIR_OPTIONS_D, "temp-sinker-center-y", '\0',
    &center_y, 0.5,
    "center of the sinker: y",
  YMIR_OPTIONS_D, "temp-sinker-center-z", '\0',
    &center_z, 0.5,
    "center of the sinker: z",
  YMIR_OPTIONS_D, "temp-sinker-diameter", '\0',
    &diameter, 0.1,
    "diameter of the sinker",
  YMIR_OPTIONS_D, "temp-sinker-scaling", '\0',
    &sinker_scaling, 0.5,
    "scaling of the sinker temperature",
  YMIR_OPTIONS_D, "temp-sinker-decay", '\0',
    &sinker_decay, 100.0,
    "decay (exponent) of the sinker temperature",


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

  /* viscosity */
  YMIR_OPTIONS_I, "viscosity-anisotropy",'\0',
    &(viscosity_anisotropy),SLABS_VISC_ISOTROPY,
    "0: isotropy, 1: transversely isotropy",
  YMIR_OPTIONS_I, "viscosity-geometry",'\0',
    &(viscosity_geometry),SLABS_VISC_LAYERS,
    "0: layers, 1: layers-coupling",
  YMIR_OPTIONS_D, "viscosity-lithosphere", '\0',
    &visc_lith, VISCOSITY_LITHOSPHERE,
    "Viscosity in the lithosphere",
  YMIR_OPTIONS_D, "viscosity-lith-trial", '\0',
    &visc_trial, 1.0,
    "Viscosity in the lithosphere",
  YMIR_OPTIONS_D, "viscosity-asthenosphere", '\0',
    &visc_asthen, VISCOSITY_ASTHENOSPHERE,
    "Viscosity in the asthenosphere",
  YMIR_OPTIONS_D, "visc-zlocation-lithosphere", '\0',
    &visc_z_lith, VISCOSITY_LITHOSPHERE_RADIUS_LOCATION,
    "Viscosity lithosphere location",
  YMIR_OPTIONS_D, "visc-zlocation-asthenosphere", '\0',
    &visc_z_asthen, VISCOSITY_ASTHENOSPHERE_RADIUS_LOCATION,
    "Viscosity asthenosphere location",
  YMIR_OPTIONS_D, "visc-zlocation-slab", '\0',
    &visc_z_slab, VISCOSITY_SLAB_RADIUS_LOCATION,
    "Viscosity slab location",
  YMIR_OPTIONS_D, "visc-slab-width", '\0',
    &visc_slab_width, VISCOSITY_SLAB_WIDTH,
    "Viscosity slab width",

  /* surface location  */
  YMIR_OPTIONS_I, "bound-x-function", '\0',
    &x_func, SLABS_X_FUNCTION_IDENTITY,
    "boundary location: surface topography",
  YMIR_OPTIONS_D, "surface-distortion-factor", '\0',
    &surf_dist, 1.0,
    "surface distortion factor.",

  /* velocity Dirichlet BC's */
  YMIR_OPTIONS_I, "velocity-dirichlet-bc", '\0',
    &vel_dir_bc, SLABS_VEL_DIR_BC_INOUTFLOW_SIN,
    "Velocity Dirichlet boundary condition",
  YMIR_OPTIONS_D, "flow-scaling", '\0',
    &flow_scale, COLLIDE_FLOW_SCALE,
    "scaling of velocity BC.",
  YMIR_OPTIONS_D, "velocity-bc-middle", '\0',
    &velocity_bc_middle, COLLIDE_ZERO_POINT_LOCATION_MIDDLE,
    "location of velocity BC: middle bound",
  YMIR_OPTIONS_D, "velocity-bc-upper", '\0',
    &velocity_bc_upper, COLLIDE_ZERO_POINT_LOCATION_UPPER,
    "location of velocity BC: upper bound",
  YMIR_OPTIONS_D, "velocity-bc-lower", '\0',
    &velocity_bc_lower, COLLIDE_ZERO_POINT_LOCATION_LOWER,
    "location of velocity BC: lower bound",

  /* test options: stress operator test and manufactured solution test */
  YMIR_OPTIONS_I, "test-stress-operator", '\0',
    &test_stress_op, SLABS_TEST_STRESS_OP_NONE,
    "the input velocity for stress operator test",
  YMIR_OPTIONS_I, "test-manufactured-solution", '\0',
    &test_manufactured, SLABS_TEST_MANUFACTURED_NONE,
    "the input for velocity and pressure field for manufactured solution test",
  YMIR_OPTIONS_I, "test-stress-component", '\0',
    &test_stress_comp, SLABS_TEST_MANUFACTURED_NONE,
    "the input for velocity and pressure field for manufactured solution test",

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
  YMIR_OPTIONS_S, "vtk-write-postp-path", '\0',
    &(vtk_write_postp_path), NULL,
    "File path for vtk files of post processing data",
  YMIR_OPTIONS_S, "vtk-write-test-path", '\0',
    &(vtk_write_test_path), NULL,
    "File path for vtk files of test: stress_op or manufactured",
  YMIR_OPTIONS_S, "vtk-write-freesurface-path", '\0',
    &(vtk_write_freesurface_path), NULL,
    "File path for vtk files of free-surface data",
  YMIR_OPTIONS_S, "vtk-write-input2-path", '\0',
    &(vtk_write_input2_path), NULL,
    "File path for vtk files for the input of the Stokes problem",
  YMIR_OPTIONS_S, "vtk-write-solution2-path", '\0',
    &(vtk_write_solution2_path), NULL,
    "File path for vtk files for the solution of the Stokes problem",
  YMIR_OPTIONS_S, "vtk-write-io-path", '\0',
    &(vtk_write_io_path), NULL,
    "File path for the test of io",
  YMIR_OPTIONS_S, "vtk-write-mpiio-path", '\0',
    &(vtk_write_mpiio_path), NULL,
    "File path for the test of mpiio",
  YMIR_OPTIONS_S, "vtk-read-io-path", '\0',
    &(vtk_read_io_path), NULL,
    "File path for the test of io",
  YMIR_OPTIONS_S, "vtk-read-mpiio-path", '\0',
    &(vtk_read_mpiio_path), NULL,
    "File path for the test of mpiio",
  YMIR_OPTIONS_S, "vtk-write-ioface-path", '\0',
    &(vtk_write_ioface_path), NULL,
    "File path for the test of face_mesh io",
  YMIR_OPTIONS_S, "vtk-read-ioface-path", '\0',
    &(vtk_read_ioface_path), NULL,
    "File path for the test of face_mesh io",
  YMIR_OPTIONS_S, "ascii-read-topo-path", '\0',
    &(ascii_read_topo_path), NULL,
    "File path for reading topography path",
  YMIR_OPTIONS_S, "vtk-write-testdist-path", '\0',
    &(vtk_write_testdist_path), NULL,
    "File path for the test of distorted surface",

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
  /* buoyancy type */
  if (test_manufactured)
    buoyancy_type = 4;
  slabs_options.buoyancy_type = (slabs_buoyancy_type_t) buoyancy_type;

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

  slabs_temp_options.custom_type = (slabs_temp_custom_t) temp_custom;
  slabs_sinker_options.center_x = center_x;
  slabs_sinker_options.center_y = center_y;
  slabs_sinker_options.center_z = center_z;
  slabs_sinker_options.diameter = diameter;
  slabs_sinker_options.scaling = sinker_scaling;
  slabs_sinker_options.decay = sinker_decay;
  slabs_temp_options.sinker_options = &slabs_sinker_options;

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

  /* viscosity */
  slabs_visc_options.viscosity_anisotropy = (slabs_viscosity_anisotropy_t) viscosity_anisotropy;
  slabs_visc_options.viscosity_geometry = (slabs_viscosity_geometry_t) viscosity_geometry;
  slabs_visc_options.visc_lith = visc_lith;
  slabs_visc_options.visc_asthen = visc_asthen;
  slabs_visc_options.z_lith = visc_z_lith;
  slabs_visc_options.z_asthen = visc_z_asthen;
  slabs_visc_options.z_slab = visc_z_slab;
  slabs_visc_options.slab_width = visc_slab_width;

  /*geometry transformation: surface location*/
  slabs_surf_options.x_func = (slabs_x_func_t) x_func;

  /* velocity B.C. condition */
  slabs_velbc_options.vel_dir_bc = (slabs_vel_dir_bc_t) vel_dir_bc;
  slabs_velbc_options.flow_scale = flow_scale;
  slabs_velbc_options.vel_dir_bc_middle = velocity_bc_middle;
  slabs_velbc_options.vel_dir_bc_upper = velocity_bc_upper;
  slabs_velbc_options.vel_dir_bc_lower = velocity_bc_lower;

  /* test */
  slabs_test_options.test_stress_op = (slabs_test_stress_op_t) test_stress_op;
  slabs_test_options.test_manufactured = (slabs_test_manufactured_t) test_manufactured;
  slabs_test_options.test_stress_comp = (slabs_test_manufactured_t) test_stress_comp;

  /* assign slabs_options */
  slabs_options.slabs_temp_options = &slabs_temp_options;
  slabs_options.slabs_visc_options = &slabs_visc_options;
  slabs_options.slabs_weak_options = &slabs_weak_options;
  slabs_options.slabs_surf_options = &slabs_surf_options;
  slabs_options.slabs_velbc_options = &slabs_velbc_options;
  slabs_options.slabs_test_options = &slabs_test_options;

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
  rhea_process_options_all (&domain_options, &temp_options, NULL,
                            &weak_options, &topo_options, &visc_options,
                            &discr_options);

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
  if (ascii_read_topo_path != NULL) {
    int                   Ncn;
    int                   i;
    FILE                  *infile;
    char                  infilename[BUFSIZ];

    RHEA_GLOBAL_PRODUCTIONF ("Read topography from %s\n", ascii_read_topo_path);

    tX = RHEA_ALLOC (double, Ncn);
    tY = RHEA_ALLOC (double, Ncn);
    tZ = RHEA_ALLOC (double, Ncn);

    snprintf (infilename, BUFSIZ, "%s_visc", ascii_read_topo_path);
    infile = fopen (infilename, "r");
    if (infile == NULL) {
      YMIR_LERRORF ("Could not open %s for reading!\n", infilename);
      return -1;
    }
    fscanf (infile, "%d\n", &Ncn);
    for (i = 0; i < Ncn; i++) {
      fscanf (infile, "%lf %lf %lf\n", &tX[i], &tY[i], &tZ[i]);
    }
    if (fclose (infile)) {
      YMIR_LERROR ("main: Error closing footer\n");
      return -1;
    }

      /*store topography information in slabs_surf_options*/
    topo.tX = tX;
    topo.tY = tY;
    topo.tZ = tZ;
    topo.nsurf = Ncn;
    slabs_surf_options.topo_profile = &topo;

  /* discr influences domain ? */
    rhea_discretization_process_options (&discr_options, &domain_options,
                                         &topo_options);

    rhea_discretization_set_X_fn (&discr_options, slabs_X_fn_profile, &topo);
    slabs_setup_mesh (&p4est, &ymir_mesh, &press_elem, mpicomm,
                          &domain_options, &discr_options, &slabs_options);

    RHEA_GLOBAL_PRODUCTIONF ("Done read topography from %s\n", ascii_read_topo_path);
  }
  else if (vtk_write_testdist_path != NULL) {
    p4est_t            *p4est_I;
    ymir_mesh_t        *ymir_mesh_I;
    ymir_pressure_elem_t  *press_elem_I;
    rhea_discretization_options_t discr_options_I;
    rhea_stokes_problem_t *stokes_problem_I;
    ymir_vec_t         *temperature_I, *weakzone_I;
    ymir_vec_t         *temperature, *weakzone;

    /*have to use a new discr_option, otherwise the code crashes with memory balance error*/
    rhea_discretization_process_options (&discr_options_I, &domain_options,
                                         &topo_options);
    slabs_setup_mesh (&p4est_I, &ymir_mesh_I, &press_elem_I, mpicomm,
                      &domain_options, &discr_options_I, &slabs_options);

//    /* compute temperature */
    temperature_I = rhea_temperature_new (ymir_mesh_I);
    weakzone_I = rhea_viscosity_new (ymir_mesh_I);

    slabs_poly2_temperature_compute (temperature_I, &slabs_options);
    slabs_weakzone_compute (weakzone_I, &slabs_options);

    rhea_discretization_set_X_fn (&discr_options, slabs_X_fn_function, NULL);

    slabs_setup_mesh (&p4est, &ymir_mesh, &press_elem, mpicomm,
                      &domain_options, &discr_options, &slabs_options);

    temperature = rhea_temperature_new (ymir_mesh);
    weakzone = rhea_viscosity_new (ymir_mesh);
    temperature->cd = temperature_I->cd;
    temperature->cdo = temperature_I->cdo;
    weakzone->dd = weakzone_I->dd;
    slabs_setup_stokes_surfdist (&stokes_problem, ymir_mesh, press_elem,
                        &domain_options, &temp_options, &visc_options,
                        &slabs_options, vtk_write_input_path, temperature, weakzone);
//
//      char            path[BUFSIZ];
//      snprintf (path, BUFSIZ, "%s", vtk_write_testdist_path);
//      ymir_vtk_write (ymir_mesh, path,
//                    temperature_I, "temperature",
//                    temperature_D, "temperature_D",
//                    NULL);

    rhea_temperature_destroy (temperature_I);
    rhea_temperature_destroy (weakzone_I);
//    rhea_temperature_destroy (temperature_D);
//    rhea_temperature_destroy (weakzone_D);
    rhea_discretization_ymir_mesh_destroy (ymir_mesh_I, press_elem_I);
    rhea_discretization_p4est_destroy (p4est_I);
    rhea_discretization_boundary_clear (&discr_options_I);
  }
  else {
    rhea_discretization_set_X_fn (&discr_options, slabs_X_fn_identity, NULL);
    slabs_setup_mesh (&p4est, &ymir_mesh, &press_elem, mpicomm,
                      &domain_options, &discr_options, &slabs_options);
 /*
   * Setup Stokes Problem
   */

  /*try different viscosity on the top layer*/
//  slabs_visc_options.visc_lith *= visc_trial;

  slabs_setup_stokes (&stokes_problem, ymir_mesh, press_elem,
                      &domain_options, &temp_options, &weak_options,
                      &visc_options, &slabs_options, vtk_write_input_path);

  }
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
                             viscosity, NULL, NAN);

    rhea_pressure_destroy (pressure);
    rhea_viscosity_destroy (viscosity);
    rhea_velocity_destroy (velocity);
  }

  if (vtk_write_test_path != NULL)  {
    if(slabs_test_options.test_manufactured != SLABS_TEST_MANUFACTURED_NONE) {
      ymir_vec_t         *vel_ref = rhea_velocity_new (ymir_mesh);
      ymir_vec_t         *vel_chk = rhea_velocity_new (ymir_mesh);
      ymir_vec_t         *pres_ref = ymir_cvec_new (ymir_mesh, 1);
      ymir_vec_t         *pres_chk = rhea_pressure_new (ymir_mesh, press_elem);
      double              vel_abs_err, vel_rel_err;
      ymir_vec_t         *vel_err = rhea_velocity_new (ymir_mesh);
      ymir_vec_t         *pres_err = ymir_cvec_new (ymir_mesh, 1);
      const               slabs_test_manufactured_t
                          test_type = slabs_test_options.test_manufactured;
      ymir_stokes_op_t   *stokes_op;
      ymir_stress_op_t   *stress_op;

      stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
      stress_op = stokes_op->stress_op;
      ymir_vec_set_value (pres_ref, .0);

      /* compute velocity fields */
      switch (test_type) {
        case SLABS_TEST_MANUFACTURED_SINCOS1_ISO:
        case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT90:
        case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT45:
        case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT60:
        case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT60_VISCEXP60:
          /* compute reference velocity field (output) */
          ymir_cvec_set_function (vel_ref, slabs_test_sincos1_vel_in_fn,
                                  NULL);
          break;

        case SLABS_TEST_MANUFACTURED_POLY1_TIROT90:
        case SLABS_TEST_MANUFACTURED_POLY1_TIROT90_VISCEXP:
          /* compute reference velocity field (output) */
          ymir_cvec_set_function (vel_ref, slabs_test_poly1_vel_in_fn,
                                  NULL);

          break;

       default:
          RHEA_ABORT_NOT_REACHED ();
      }

      ymir_stokes_vec_get_components (sol_vel_press, vel_chk, pres_chk,
                                        press_elem);
      slabs_test_manufactured_compute_vel_err (&vel_abs_err, &vel_rel_err,
                                                 vel_err, vel_ref, vel_chk,
                                                 stress_op);

      RHEA_GLOBAL_INFOF ("manufactured solution test (test type %i): abs error %1.3e rel error %1.3e\n",
                          test_type, vel_abs_err, vel_rel_err);


      char      path[BUFSIZ];
      snprintf (path, BUFSIZ, "%s_manufactured", vtk_write_test_path);
      ymir_vtk_write (ymir_mesh, path,
                      vel_ref, "vel_reference",
                      vel_chk, "vel_check",
                      vel_err, "vel_error",
                      pres_ref, "pressure_reference",
                      pres_chk, "pressure_check",
                      NULL);
      rhea_velocity_destroy (vel_ref);
      rhea_velocity_destroy (vel_chk);
      rhea_velocity_destroy (vel_err);
      rhea_pressure_destroy (pres_chk);
      ymir_vec_destroy (pres_ref);
      ymir_vec_destroy (pres_err);
    }

    if (slabs_test_options.test_stress_comp != SLABS_TEST_MANUFACTURED_NONE)
    {
      const               slabs_test_manufactured_t
                          test_type = slabs_test_options.test_stress_comp;

      ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (ymir_mesh->ma->N,
                                                                ymir_mesh->ma->ompsize);
      ymir_vec_t         *velocity = rhea_velocity_new (ymir_mesh);
      ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);
      ymir_vec_t         *shear_visc = rhea_viscosity_new (ymir_mesh);
      ymir_vec_t         *TI_tensor = ymir_dvec_new (ymir_mesh, 9, YMIR_GAUSS_NODE);
      ymir_vec_t         *trac_n = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
      ymir_vec_t         *trac_s = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
      ymir_vec_t         *normal_force = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
      ymir_vec_t         *shear_force = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
      ymir_vec_t         *traction_ref = rhea_velocity_new (ymir_mesh);
      ymir_vec_t         *stress = ymir_dvec_new (ymir_mesh, 3, YMIR_GAUSS_NODE);
      ymir_vec_t         *stress_ref = rhea_velocity_new (ymir_mesh);
      RHEA_GLOBAL_PRODUCTIONF ("In %s: Start test stress_component\n", this_fn_name);

      ymir_stokes_vec_get_velocity (sol_vel_press, velocity,
                                      press_elem);
      rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

      if (slabs_visc_options.viscosity_anisotropy == SLABS_VISC_TRANSVERSELY_ISOTROPY)  {
        ymir_stokes_op_t      *stokes_op;
        ymir_stress_op_t      *stress_op;

        /* get the viscous stress operator */
        stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
        stress_op = stokes_op->stress_op;

        /* copy shear viscosity */
        slabs_stress_op_copy_shear_visc (shear_visc, stress_op);
        slabs_stress_op_copy_TI_tensor (TI_tensor, stress_op);

        slabs_manufactured_strainvec_coupling_compute (velocity, trac_n, trac_s,
                                               shear_visc, viscosity);

        slabs_manufactured_stressvec_coupling_compute (velocity, viscosity,
                                              shear_visc, TI_tensor,
                                              normal_force, shear_force,
                                              vel_elem, &slabs_options);

        slabs_stressvec_TI (velocity, stress, viscosity, shear_visc, TI_tensor, vel_elem);
        switch (test_type) {
          case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT90:
            ymir_cvec_set_function (traction_ref, slabs_test_sincos1_TIrot90_traction_fn,
                                      NULL);
            ymir_cvec_set_function (stress_ref, slabs_test_sincos1_TIrot90_stress_fn,
                                      NULL);
          case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT45:
            ymir_cvec_set_function (traction_ref, slabs_test_sincos1_TIrot45_traction_fn,
                                      NULL);
            ymir_cvec_set_function (stress_ref, slabs_test_sincos1_TIrot45_stress_fn,
                                      NULL);
          case SLABS_TEST_MANUFACTURED_SINCOS1_TIROT60:
            ymir_cvec_set_function (traction_ref, slabs_test_sincos1_TIrot60_traction_fn,
                                      NULL);
            ymir_cvec_set_function (stress_ref, slabs_test_sincos1_TIrot60_stress_fn,
                                      NULL);
            break;

          default:
            RHEA_ABORT_NOT_REACHED ();
        }
      }
      char      path[BUFSIZ];
      snprintf (path, BUFSIZ, "%s_stress_component", vtk_write_test_path);
      ymir_vtk_write (ymir_mesh, path,
                      stress_ref, "stress_reference (yy,zz,yz)",
                      stress, "stress_check",
                      traction_ref, "traction_reference (0,fn,fs)",
                      normal_force, "normal_force (from stress)",
                      shear_force, "shear_force (from stress)",
                      trac_n, "normal_force (from strain)",
                      trac_s, "shear_force (from strain)",
                      NULL);

      ymir_velocity_elem_destroy (vel_elem);
      rhea_velocity_destroy (velocity);
      rhea_viscosity_destroy (viscosity);
      rhea_viscosity_destroy (shear_visc);
      ymir_vec_destroy (TI_tensor);
      ymir_vec_destroy (trac_n);
      ymir_vec_destroy (trac_s);
      ymir_vec_destroy (normal_force);
      ymir_vec_destroy (shear_force);
      ymir_vec_destroy (traction_ref);
      ymir_vec_destroy (stress);
      ymir_vec_destroy (stress_ref);
    } /* end of stress component*/
    /*stress operator test, TODO*/
  }

  /* compute and output second invariant strain_rate, stress, and surface normal stress  */
  if (vtk_write_stress_path != NULL)  {
    ymir_vec_t         *velocity = rhea_velocity_new (ymir_mesh);
    ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);
    ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (ymir_mesh->ma->N,
                                                              ymir_mesh->ma->ompsize);
    ymir_vec_t            *edotII = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
    ymir_vec_t            *tauII = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
    ymir_vec_t            *surf_normal_stress = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    ymir_vec_t            *vec_e = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    ymir_vec_t            *masse = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    ymir_vec_t            *massu = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    ymir_vec_t            *vec_topo = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    double                avg_stress, tempe, tempu;

    RHEA_GLOBAL_PRODUCTIONF ("In %s: Start vtk_write_stress\n", this_fn_name);

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

    /*compute average surface normal stress
     * vec_e is 1.0
     * average = (e, M*u)/(e, M*e)*/
    ymir_vec_set_value(vec_e, 1.0)
    ymir_mass_apply (surf_normal_stress, massu);
    ymir_mass_apply (vec_e, masse);
    tempe = ymir_vec_innerprod (vec_e, masse);
    tempu = ymir_vec_innerprod (vec_e, massu);
    avg_stress = tempu/tempe;

    /*scale normal stress to topography*/
    ymir_vec_copy (surf_normal_stress, vec_topo);
    ymir_vec_shift (-avg_stress, vec_topo);
    /*vec_topo = (-2.0/Ra) * (nstress-avg_stress) + 1.0, here Ra=1 */
    ymir_vec_scale_shift (-2.0, 1.0, vec_topo);

    {
      char            path[BUFSIZ];

      snprintf (path, BUFSIZ, "%s", vtk_write_stress_path);
      ymir_vtk_write (ymir_mesh, path,
                      edotII, "edotII",
                      tauII, "tauII",
                      surf_normal_stress, "surf_normal_stress",
                      vec_topo, "topography",
                      NULL);
    }

    /* destroy */
    rhea_viscosity_destroy (viscosity);
    rhea_velocity_destroy (velocity);
    ymir_vec_destroy (edotII);
    ymir_vec_destroy (tauII);
    ymir_vec_destroy (surf_normal_stress);
    ymir_vec_destroy (vec_e);
    ymir_vec_destroy (masse);
    ymir_vec_destroy (massu);
    ymir_vec_destroy (vec_topo);
    ymir_velocity_elem_destroy (vel_elem);

    RHEA_GLOBAL_PRODUCTIONF ("In %s: Done vtk_write_stress\n", this_fn_name);
  }

    /* compute and output analysis of stress */
  if (vtk_write_postp_path != NULL)  {
    ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (ymir_mesh->ma->N,
                                                              ymir_mesh->ma->ompsize);
    ymir_vec_t         *velocity = rhea_velocity_new (ymir_mesh);
    ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *shear_visc = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *TI_tensor = ymir_dvec_new (ymir_mesh, 9, YMIR_GAUSS_NODE);
    ymir_vec_t         *traction = ymir_dvec_new (ymir_mesh, 3, YMIR_GAUSS_NODE);
    ymir_vec_t         *normal_force = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
    ymir_vec_t         *shear_force = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
    RHEA_GLOBAL_PRODUCTIONF ("In %s: Start vtk_write_postp\n", this_fn_name);

    ymir_stokes_vec_get_velocity (sol_vel_press, velocity,
                                    press_elem);
    rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

    if (slabs_visc_options.viscosity_anisotropy == SLABS_VISC_TRANSVERSELY_ISOTROPY)  {
      ymir_stokes_op_t      *stokes_op;
      ymir_stress_op_t      *stress_op;

      /* get the viscous stress operator */
      stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
      stress_op = stokes_op->stress_op;

      /* copy shear viscosity */
      slabs_stress_op_copy_shear_visc (shear_visc, stress_op);
      slabs_stress_op_copy_TI_tensor (TI_tensor, stress_op);

      slabs_postp_weakzone_traction_compute (velocity, traction,
                                             shear_visc, viscosity,
                                             &slabs_options);

      slabs_postp_weakzone_coupling_compute (velocity, viscosity,
                                            shear_visc, TI_tensor,
                                            normal_force, shear_force,
                                            vel_elem, &slabs_options);
    }
    else  {
      slabs_postp_weakzone_traction_compute (velocity, traction,
                                             viscosity, viscosity,
                                             &slabs_options);

      slabs_postp_weakzone_coupling_compute (velocity, viscosity,
                                             shear_visc, TI_tensor,
                                             normal_force, shear_force,
                                             vel_elem, &slabs_options);
   }
   rhea_viscosity_destroy (shear_visc);
   ymir_vec_destroy (TI_tensor);

   {
     char            path[BUFSIZ];

     snprintf (path, BUFSIZ, "%s", vtk_write_postp_path);
     ymir_vtk_write (ymir_mesh, path,
                     normal_force, "normal_force from sigma",
                     shear_force, "shear_force from sigma",
                     traction, "trac_n from edot",
                     NULL);
   }

    /* destroy */
    rhea_viscosity_destroy (viscosity);
    rhea_velocity_destroy (velocity);
    ymir_vec_destroy (normal_force);
    ymir_vec_destroy (shear_force);
    ymir_vec_destroy (traction);
    ymir_velocity_elem_destroy (vel_elem);

    RHEA_GLOBAL_PRODUCTIONF ("In %s: Done vtk_write_postp\n", this_fn_name);
  }

  /* compute and output second invariant strain_rate, stress, and surface normal stress  */
  if (vtk_write_freesurface_path != NULL)  {
    ymir_vec_t            *surf_normal_stress = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    ymir_vec_t            *surf_normal_stress2;
    ymir_vec_t            *vec_e = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    ymir_vec_t            *masse = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    ymir_vec_t            *massu = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    ymir_vec_t            *vec_topo = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    ymir_vec_t            *vec_x = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    ymir_vec_t            *vec_y = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    double                tempe, tempu;
    p4est_t               *p4est2;
    ymir_mesh_t           *ymir_mesh2;
    ymir_pressure_elem_t  *press_elem2;
    ymir_vec_t            *sol_vel_press2;
    rhea_stokes_problem_t *stokes_problem2;
    rhea_discretization_options_t discr_options2;

    ymir_topidx_t         fm;
    ymir_face_mesh_t      *fmesh;
    int                   Ncn;
    int                   i;
    double                *tX, *tY, *tZ;
    double                avg_stress, topo_nondim;
    slabs_topo_profile_t  topo;
    FILE                  *outfile;
    char                  outfilename[BUFSIZ];

    RHEA_GLOBAL_PRODUCTIONF ("In %s: Start vtk_write_freesurface\n", this_fn_name);

    /* compute surface normal stress sigma */
    slabs_physics_compute_normal_boundary_stress (
                   surf_normal_stress, sol_vel_press,
                   rhea_stokes_problem_get_rhs_vel (stokes_problem),
                   rhea_stokes_problem_get_stokes_op (stokes_problem));

    /*surface mesh information*/
    fm = surf_normal_stress->meshnum;
    fmesh = &(ymir_mesh->fmeshes[fm]);
    Ncn = fmesh->Ncn;

    /*compute average surface normal stress
     * vec_e is 1.0
     * average = (e, M*u)/(e, M*e)*/
    ymir_vec_set_value(vec_e, 1.0)
    ymir_mass_apply (surf_normal_stress, massu);
    ymir_mass_apply (vec_e, masse);
    tempe = ymir_vec_innerprod (vec_e, masse);
    tempu = ymir_vec_innerprod (vec_e, massu);
    avg_stress = tempu/tempe;

    /*scale normal stress to topography*/
    ymir_vec_copy (surf_normal_stress, vec_topo);
    ymir_vec_shift (-avg_stress, vec_topo);
    /*vec_topo = (-2.0/Ra) * (nstress-avg_stress) + 1.0, here Ra=1 */
    ymir_vec_scale_shift (-2.0, 1.0, vec_topo);

    /*obtain x and y coordinates*/
    ymir_face_cvec_set_function (vec_x, slabs_coordx_set_fn, NULL);
    ymir_face_cvec_set_function (vec_y, slabs_coordy_set_fn, NULL);
    tX = vec_x->cvec->e[0];
    tY = vec_y->cvec->e[0];
    tZ = vec_topo->cvec->e[0];

    snprintf (outfilename, BUFSIZ, "%s_visc_%04d", vtk_write_freesurface_path, mpirank);
    outfile = fopen (outfilename, "w");
    if (outfile == NULL) {
      YMIR_LERRORF ("Could not open %s for output!\n", outfilename);
      return -1;
    }
    for (i = 0; i < Ncn; i++) {
      fprintf (outfile, "%10.6e   %10.6e   %10.6e\n", tX[i], tY[i], tZ[i]);
    }
    if (fclose (outfile)) {
      YMIR_LERROR ("main: Error closing footer\n");
      return -1;
    }

    /*store topograhy information in slabs_surf_options*/
    topo.tX = tX;
    topo.tY = tY;
    topo.tZ = tZ;
    topo.nsurf = Ncn;
    slabs_surf_options.topo_profile = &topo;
    /*creat new discr and domain options with the distorted surface*/
    rhea_discretization_process_options (&discr_options2, &domain_options,
                                         &topo_options);
    rhea_discretization_set_X_fn (&discr_options2, slabs_X_fn_profile, &topo);

    slabs_setup_mesh (&p4est2, &ymir_mesh2, &press_elem2, mpicomm,
                      &domain_options, &discr_options2, &slabs_options);

    /*use free-surface B.C. on the surface*/
    vel_dir_bc = 6;
    slabs_velbc_options.vel_dir_bc = (slabs_vel_dir_bc_t) vel_dir_bc;

    /*try different viscosity on the top layer*/
    slabs_visc_options.visc_lith *= visc_trial;

    slabs_setup_stokes (&stokes_problem2, ymir_mesh2, press_elem2,
                        &domain_options, &temp_options, &weak_options,
                        &visc_options, &slabs_options, vtk_write_input2_path);

    /* initialize solution vector */
    sol_vel_press2 = rhea_velocity_pressure_new (ymir_mesh2, press_elem2);

    /* run solver */
    slabs_run_solver (sol_vel_press2, ymir_mesh2, press_elem2, stokes_problem2,
                      solver_iter_max, solver_rel_tol);

    {
      ymir_vec_t         *velocity2 = rhea_velocity_new (ymir_mesh2);
      ymir_vec_t         *pressure2 = rhea_pressure_new (ymir_mesh2, press_elem2);
      ymir_vec_t         *viscosity2 = rhea_viscosity_new (ymir_mesh2);

      ymir_stokes_vec_get_components (sol_vel_press2, velocity2, pressure2,
                                      press_elem2);
      rhea_stokes_problem_copy_viscosity (viscosity2, stokes_problem2);

      rhea_vtk_write_solution (vtk_write_solution2_path, velocity2, pressure2,
                               viscosity2, NULL, NAN);

      rhea_pressure_destroy (pressure2);
      rhea_viscosity_destroy (viscosity2);
      rhea_velocity_destroy (velocity2);
    }
#if 0
  /* compute surface normal stress sigma */
    surf_normal_stress2 = ymir_face_cvec_new (ymir_mesh2,
                                         RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);

    slabs_physics_compute_normal_boundary_stress (
                   surf_normal_stress2, sol_vel_press2,
                   rhea_stokes_problem_get_rhs_vel (stokes_problem2),
                   rhea_stokes_problem_get_stokes_op (stokes_problem2));

    {
      char            path[BUFSIZ];

      snprintf (path, BUFSIZ, "%s", vtk_write_freesurface_path);
      ymir_vtk_write (ymir_mesh, path,
                      surf_normal_stress, "surf_normal_stress",
                      surf_normal_stress2, "surf_normal_stress with f-surface",
                      vec_topo, "topography vector",
                      NULL);
    }
#endif
    /* destroy */
    ymir_vec_destroy (surf_normal_stress);
//    ymir_vec_destroy (surf_normal_stress2);
    ymir_vec_destroy (vec_e);
    ymir_vec_destroy (masse);
    ymir_vec_destroy (massu);
    ymir_vec_destroy (vec_topo);
    ymir_vec_destroy (vec_x);
    ymir_vec_destroy (vec_y);
    rhea_velocity_pressure_destroy (sol_vel_press2);
    slabs_surf_options.topo_profile = NULL;

    slabs_setup_clear_all (stokes_problem2, p4est2, ymir_mesh2, press_elem2,
                           &discr_options2);

    RHEA_GLOBAL_PRODUCTIONF ("In %s: Done vtk_write_freesurface\n", this_fn_name);
  }


  if (vtk_write_io_path != NULL) {
    mangll_cnodes_t       *cnodes = ymir_mesh->cnodes;
    const int             N = cnodes->N;
    int                   Np, mpirank;
    ymir_locidx_t         n_elements, elid, Ntotal, fm;
    int                   i, j, k;
    double                *x, *y, *z, *tmp_el, *elemd;
    sc_dmatrix_t          *elem;
    FILE                  *outfile;
    char                  outfilename[BUFSIZ];

    ymir_vec_t            *viscosity = rhea_viscosity_new (ymir_mesh);

    RHEA_GLOBAL_PRODUCTIONF ("In %s: Start vtk_write_io\n", this_fn_name);

    rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

    mpirank = ymir_mesh->ma->mpirank;

    Np = (N + 1) * (N + 1) * (N + 1);
    n_elements = cnodes->K;
    Ntotal = Np * n_elements;

    snprintf (outfilename, BUFSIZ, "%s_visc_%04d", vtk_write_io_path, mpirank);

    outfile = fopen (outfilename, "w");
    if (outfile == NULL) {
      YMIR_LERRORF ("Could not open %s for output!\n", outfilename);
      return -1;
    }

    x = RHEA_ALLOC (double, Np);
    y = RHEA_ALLOC (double, Np);
    z = RHEA_ALLOC (double, Np);
    tmp_el = RHEA_ALLOC (double, Np);

    elem = sc_dmatrix_new (Np, 1);
    elemd = elem->e[0];
    for (elid = 0; elid < n_elements; elid++) {
      elemd =  rhea_viscosity_get_elem_gauss (elem, viscosity, elid);

      ymir_mesh_get_elem_coord_gauss (x, y, z, elid, ymir_mesh, tmp_el);

      for (i = 0; i < Np; ++i)  {
        fprintf (outfile, "%24.16e    %24.16e    %24.16e    %24.16e\n", x[i],
               y[i], z[i], elemd[i]);
      }
    }
    if (fclose (outfile)) {
      YMIR_LERROR ("main: Error closing footer\n");
      return -1;
    }

    RHEA_FREE (x);
    RHEA_FREE (y);
    RHEA_FREE (z);
    RHEA_FREE (tmp_el);
    sc_dmatrix_destroy (elem);
    ymir_vec_destroy (viscosity);

  }

  if (vtk_write_mpiio_path != NULL) {
    mangll_cnodes_t       *cnodes = ymir_mesh->cnodes;
    const int             N = cnodes->N;
    int                   Np, mpirank, mpisize;
    ymir_locidx_t         n_elements, elid, Ntotal, fm;
    int                   i, j, k;
    double                *x, *y, *z, *tmp_el, *elemd;
    double                *Ax, *Ay, *Az, *Avisc;
    double                *Px, *Py, *Pz, *Pvisc;
    sc_dmatrix_t          *elem;
    FILE                  *outfile;
    char                  outfilename[BUFSIZ];

    ymir_vec_t            *viscosity = rhea_viscosity_new (ymir_mesh);

    RHEA_GLOBAL_PRODUCTIONF ("In %s: Start vtk_write_mpiio\n", this_fn_name);

    rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

    mpirank = ymir_mesh->ma->mpirank;
    mpisize = ymir_mesh->ma->mpisize;

    Np = (N + 1) * (N + 1) * (N + 1);
    n_elements = cnodes->K;
    Ntotal = Np * n_elements;

    snprintf (outfilename, BUFSIZ, "%s_visc_merged", vtk_write_mpiio_path);
    if (mpirank == 0) {
      outfile = fopen (outfilename, "w");
      if (outfile == NULL) {
        YMIR_LERRORF ("Could not open %s for output!\n", outfilename);
        return -1;
      }
      Ax = RHEA_ALLOC (double, mpisize * n_elements * Np);
      Ay = RHEA_ALLOC (double, mpisize * n_elements * Np);
      Az = RHEA_ALLOC (double, mpisize * n_elements * Np);
      Avisc = RHEA_ALLOC (double, mpisize * n_elements * Np);
    }

    Px = RHEA_ALLOC (double, n_elements * Np);
    Py = RHEA_ALLOC (double, n_elements * Np);
    Pz = RHEA_ALLOC (double, n_elements * Np);
    Pvisc = RHEA_ALLOC (double, n_elements * Np);

    x = RHEA_ALLOC (double, Np);
    y = RHEA_ALLOC (double, Np);
    z = RHEA_ALLOC (double, Np);
    tmp_el = RHEA_ALLOC (double, Np);

    elem = sc_dmatrix_new (Np, 1);
    elemd = elem->e[0];
    for (elid = 0; elid < n_elements; elid++) {
      elemd =  rhea_viscosity_get_elem_gauss (elem, viscosity, elid);

      ymir_mesh_get_elem_coord_gauss (x, y, z, elid, ymir_mesh, tmp_el);

      for (i = 0; i < Np; ++i)  {
        Px[elid * Np + i] = x[i];
        Py[elid * Np + i] = y[i];
        Pz[elid * Np + i] = z[i];
        Pvisc[elid * Np + i] = elemd[i];
      }
    }
    MPI_Gather (Px, n_elements*Np, MPI_DOUBLE,
                Ax, n_elements*Np, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Gather (Py, n_elements*Np, MPI_DOUBLE,
                Ay, n_elements*Np, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Gather (Pz, n_elements*Np, MPI_DOUBLE,
                Az, n_elements*Np, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Gather (Pvisc, n_elements*Np, MPI_DOUBLE,
                Avisc, n_elements*Np, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    if (mpirank == 0) {
      for (i = 0; i < mpisize*n_elements*Np; i++) {
        fprintf (outfile, "%24.16e    %24.16e    %24.16e    %24.16e\n", Ax[i],
                 Ay[i], Az[i], Avisc[i]);
      }
      if (fclose (outfile)) {
        YMIR_LERROR ("main: Error closing footer\n");
        return -1;
      }
      RHEA_FREE (Ax);
      RHEA_FREE (Ay);
      RHEA_FREE (Az);
      RHEA_FREE (Avisc);
    }

    char            path[BUFSIZ];

    snprintf (path, BUFSIZ, "%s_original", vtk_write_mpiio_path);
    ymir_vtk_write (ymir_mesh, path,
                    viscosity, "original_visc",
                    NULL);

    RHEA_FREE (x);
    RHEA_FREE (y);
    RHEA_FREE (z);
    RHEA_FREE (tmp_el);
    RHEA_FREE (Px);
    RHEA_FREE (Py);
    RHEA_FREE (Pz);
    RHEA_FREE (Pvisc);
    sc_dmatrix_destroy (elem);
    ymir_vec_destroy (viscosity);
  }


  if (vtk_read_io_path != NULL)  {
    mangll_cnodes_t       *cnodes = ymir_mesh->cnodes;
    const int             N = cnodes->N;
    int                   Np, mpirank;
    ymir_locidx_t         n_elements, elid, Ntotal, fm;
    int                   i, j, k;
    double                *elemd;
    double                tmp;
    sc_dmatrix_t          *elem;
    FILE                  *infile;
    char                  infilename[BUFSIZ];
    ymir_vec_t            *readin_visc = rhea_viscosity_new (ymir_mesh);

    RHEA_GLOBAL_PRODUCTIONF ("In %s: Start vtk_read_io\n", this_fn_name);

    mpirank = ymir_mesh->ma->mpirank;

    Np = (N + 1) * (N + 1) * (N + 1);
    n_elements = cnodes->K;
    Ntotal = Np * n_elements;

    snprintf (infilename, BUFSIZ, "%s_visc_%04d", vtk_read_io_path, mpirank);

    infile = fopen (infilename, "r");
    if (infile == NULL) {
      YMIR_LERRORF ("Could not open %s for reading!\n", infilename);
      return -1;
    }

    elem = sc_dmatrix_new (Np, 1);
    elemd = elem->e[0];
    for (elid = 0; elid < n_elements; elid++) {
      for (i = 0; i < Np; ++i)  {
        fscanf (infile, "%lf %lf %lf %lf\n", &tmp,
               &tmp, &tmp, &elemd[i]);
      }
      rhea_viscosity_set_elem_gauss (readin_visc, elem, elid);

    }
    if (fclose (infile)) {
      YMIR_LERROR ("main: Error closing footer\n");
      return -1;
    }


    sc_dmatrix_destroy (elem);
    ymir_vec_destroy (readin_visc);
  }

  if (vtk_read_mpiio_path != NULL)  {
    mangll_cnodes_t       *cnodes = ymir_mesh->cnodes;
    const int             N = cnodes->N;
    int                   Np, mpirank;
    ymir_locidx_t         n_elements, elid, Ntotal, fm;
    int                   i, j, k;
    double                *elemd;
    double                *Avisc, *Pvisc;
    double                tmp, norm_diff;
    sc_dmatrix_t          *elem;
    FILE                  *infile;
    char                  infilename[BUFSIZ];
    ymir_vec_t            *readin_visc = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t            *viscosity = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t            *diff_visc = rhea_viscosity_new (ymir_mesh);

    RHEA_GLOBAL_PRODUCTIONF ("In %s: Start vtk_read_mpiio\n", this_fn_name);

    rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

    mpirank = ymir_mesh->ma->mpirank;

    Np = (N + 1) * (N + 1) * (N + 1);
    n_elements = cnodes->K;
    Ntotal = Np * n_elements;

    snprintf (infilename, BUFSIZ, "%s_visc_merged", vtk_read_mpiio_path);
    if (mpirank == 0) {
      infile = fopen (infilename, "r");
      if (infile == NULL) {
        YMIR_LERRORF ("Could not open %s for reading!\n", infilename);
        return -1;
      }
      Avisc = RHEA_ALLOC (double, mpisize * n_elements * Np);
      for (i = 0; i < mpisize*n_elements*Np; ++i)  {
        fscanf (infile, "%lf %lf %lf %lf\n", &tmp,
               &tmp, &tmp, &Avisc[i]);
      }
    }

    Pvisc = RHEA_ALLOC (double, n_elements * Np);
    MPI_Scatter (Avisc, n_elements*Np, MPI_DOUBLE,
                 Pvisc, n_elements*Np, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    elem = sc_dmatrix_new (Np, 1);
    elemd = elem->e[0];
    for (elid = 0; elid < n_elements; elid++) {
      for (i = 0; i < Np; ++i)  {
        elemd[i] = Pvisc[elid * Np + i];
      }
      rhea_viscosity_set_elem_gauss (readin_visc, elem, elid);
    }

    if (mpirank == 0) {
      if (fclose (infile)) {
        YMIR_LERROR ("main: Error closing footer\n");
        return -1;
      }
      RHEA_FREE (Avisc);
    }

    ymir_vec_copy (readin_visc, diff_visc);
    ymir_vec_add (-1.0, viscosity, diff_visc);
    norm_diff = ymir_vec_norm (diff_visc);
    RHEA_GLOBAL_INFOF ("mpiio_read: norm_diff=%f\n", norm_diff);

    char            path[BUFSIZ];

    snprintf (path, BUFSIZ, "%s_readin", vtk_read_mpiio_path);
    ymir_vtk_write (ymir_mesh, path,
                    readin_visc, "readin_visc",
                    diff_visc, "diff_visc",
                    NULL);


    RHEA_FREE (Pvisc);
    sc_dmatrix_destroy (elem);
    ymir_vec_destroy (readin_visc);
    ymir_vec_destroy (viscosity);
    ymir_vec_destroy (diff_visc);
  }

  if (vtk_write_ioface_path != NULL) {
    ymir_vec_t            *surf_normal_stress = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    mangll_cnodes_t       *cnodes = ymir_mesh->cnodes;
    const int             N = cnodes->N;
    int                   Np;
    int                   mpirank = ymir_mesh->ma->mpirank;
    ymir_locidx_t         n_elements, elid, Ntotal;
    ymir_locidx_t         fm = surf_normal_stress->meshnum;
    ymir_face_mesh_t      *fmesh = &(ymir_mesh->fmeshes[fm]);
    int                   i, j, k;
    double                *Xd, *Yd, *Zd;
    double                *elemd;
    sc_dmatrix_t          *elem;
    FILE                  *outfile;
    char                  outfilename[BUFSIZ];


    RHEA_GLOBAL_PRODUCTIONF ("In %s: Start vtk_write_ioface\n", this_fn_name);

    /* compute surface normal stress sigma */
    slabs_physics_compute_normal_boundary_stress (
                   surf_normal_stress, sol_vel_press,
                   rhea_stokes_problem_get_rhs_vel (stokes_problem),
                   rhea_stokes_problem_get_stokes_op (stokes_problem));

    Np = (N + 1) * (N + 1);
    n_elements = fmesh->K;
    Ntotal = Np * n_elements;
    Xd = fmesh->X->e[0];
    Yd = fmesh->Y->e[0];
    Zd = fmesh->Z->e[0];

    snprintf (outfilename, BUFSIZ, "%s_nstress_%04d.face%d", vtk_write_ioface_path, mpirank,
              (int) fm);

    outfile = fopen (outfilename, "w");
    if (outfile == NULL) {
      YMIR_LERRORF ("Could not open %s for output!\n", outfilename);
      return -1;
    }

    elem = sc_dmatrix_new (Np, 1);
    elemd = elem->e[0];
    for (elid = 0; elid < n_elements; elid++) {
      ymir_cvec_get_elem_interp (surf_normal_stress, elem, YMIR_STRIDE_NODE, elid,
                                 YMIR_GAUSS_NODE, YMIR_COPY);

      for (i = 0; i < Np; ++i)  {
        j = elid * Np + i;
        fprintf (outfile, "%24.16e    %24.16e    %24.16e    %24.16e\n", Xd[j],
               Yd[j], Zd[j], elemd[i]);
      }
    }
    if (fclose (outfile)) {
      YMIR_LERROR ("main: Error closing footer\n");
      return -1;
    }

    char            path[BUFSIZ];

    snprintf (path, BUFSIZ, "%s_original", vtk_write_ioface_path);
    ymir_vtk_write (ymir_mesh, path,
                    surf_normal_stress, "surf_normal_stress",
                    NULL);

    sc_dmatrix_destroy (elem);
    ymir_vec_destroy (surf_normal_stress);
  }

#if 0
  if (vtk_write_mpiioface_path != NULL) {
    ymir_vec_t            *surf_normal_stress = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    mangll_cnodes_t       *cnodes = ymir_mesh->cnodes;
    const int             N = cnodes->N;
    int                   Np;
    int                   mpirank = ymir_mesh->ma->mpirank;
    int                   mpirank = ymir_mesh->ma->mpisize;
    ymir_locidx_t         n_elements, elid, Ntotal;
    ymir_locidx_t         fm = surf_normal_stress->meshnum;
    ymir_face_mesh_t      *fmesh = &(ymir_mesh->fmeshes[fm]);
    int                   i, j, k;
    double                *Xd, *Yd, *Zd;
    double                *Ax, *Ay, *Az, *Asurf;
    double                *Px, *Py, *Pz, *Psurf;
    double                *elemd;
    sc_dmatrix_t          *elem;
    FILE                  *outfile;
    char                  outfilename[BUFSIZ];


    RHEA_GLOBAL_PRODUCTIONF ("In %s: Start vtk_write_ioface\n", this_fn_name);

    /* compute surface normal stress sigma */
    slabs_physics_compute_normal_boundary_stress (
                   surf_normal_stress, sol_vel_press,
                   rhea_stokes_problem_get_rhs_vel (stokes_problem),
                   rhea_stokes_problem_get_stokes_op (stokes_problem));

    Np = (N + 1) * (N + 1);
    n_elements = fmesh->K;
    Ntotal = Np * n_elements;
    Xd = fmesh->X->e[0];
    Yd = fmesh->Y->e[0];
    Zd = fmesh->Z->e[0];

    snprintf (outfilename, BUFSIZ, "%s_nstress_%04d.face%d", vtk_write_ioface_path, mpirank,
              (int) fm);

    outfile = fopen (outfilename, "w");
    if (outfile == NULL) {
      YMIR_LERRORF ("Could not open %s for output!\n", outfilename);
      return -1;
    }

    elem = sc_dmatrix_new (Np, 1);
    elemd = elem->e[0];
    for (elid = 0; elid < n_elements; elid++) {
      ymir_cvec_get_elem_interp (surf_normal_stress, elem, YMIR_STRIDE_NODE, elid,
                                 YMIR_GAUSS_NODE, YMIR_COPY);

      for (i = 0; i < Np; ++i)  {
        j = elid * Np + i;
        fprintf (outfile, "%24.16e    %24.16e    %24.16e    %24.16e\n", Xd[j],
               Yd[j], Zd[j], elemd[i]);
      }
    }
    if (fclose (outfile)) {
      YMIR_LERROR ("main: Error closing footer\n");
      return -1;
    }

    char            path[BUFSIZ];

    snprintf (path, BUFSIZ, "%s_original", vtk_write_ioface_path);
    ymir_vtk_write (ymir_mesh, path,
                    surf_normal_stress, "surf_normal_stress",
                    NULL);

    sc_dmatrix_destroy (elem);
    ymir_vec_destroy (surf_normal_stress);
  }
#endif

  if (vtk_read_ioface_path != NULL) {
    ymir_vec_t            *surf_normal_stress = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    mangll_cnodes_t       *cnodes = ymir_mesh->cnodes;
    const int             N = cnodes->N;
    int                   Np;
    int                   mpirank = ymir_mesh->ma->mpirank;
    ymir_locidx_t         n_elements, elid, Ntotal;
    ymir_locidx_t         fm = surf_normal_stress->meshnum;
    ymir_face_mesh_t      *fmesh = &(ymir_mesh->fmeshes[fm]);
    int                   i, j, k;
    double                *elemd;
    double                tmp;
    sc_dmatrix_t          *elem;
    FILE                  *infile;
    char                  infilename[BUFSIZ];


    RHEA_GLOBAL_PRODUCTIONF ("In %s: Start vtk_read_ioface\n", this_fn_name);

    Np = (N + 1) * (N + 1);
    n_elements = fmesh->K;
    Ntotal = Np * n_elements;

    snprintf (infilename, BUFSIZ, "%s_nstress_%04d.face%d", vtk_read_ioface_path, mpirank,
              (int) fm);

    infile = fopen (infilename, "r");
    if (infile == NULL) {
      YMIR_LERRORF ("Could not open %s for output!\n", infilename);
      return -1;
    }

    elem = sc_dmatrix_new (Np, 1);
    elemd = elem->e[0];
    for (elid = 0; elid < n_elements; elid++) {
      for (i = 0; i < Np; ++i) {
        fscanf(infile, "%lf %lf %lf %lf\n", &tmp, &tmp, &tmp, &elemd[i]);
      }
      ymir_cvec_set_elem_interp (surf_normal_stress, elem, YMIR_STRIDE_NODE, elid,
                                 YMIR_GAUSS_NODE, YMIR_SET);
    }
    if (fclose (infile)) {
      YMIR_LERROR ("main: Error closing footer\n");
      return -1;
    }

    char        path[BUFSIZ];
    snprintf (path, BUFSIZ, "%s_readin", vtk_read_ioface_path);
    ymir_vtk_write (ymir_mesh, path,
                    surf_normal_stress, "readin_surf_normal_stress",
                    NULL);

    sc_dmatrix_destroy (elem);
    ymir_vec_destroy (surf_normal_stress);
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

  if (ascii_read_topo_path != NULL) {
    RHEA_FREE(tX);
    RHEA_FREE(tY);
    RHEA_FREE(tZ);
    slabs_surf_options.topo_profile = NULL;
  }


  /* print that this function is ending */
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

  /* finalize rhea */
  rhea_finalize ();

  return 0;
}
