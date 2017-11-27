
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
  TEST_MANUFACTURED
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
  SLABS_VEL_DIR_BC_MANUFACTURED_POLY_ANISO
}
slabs_vel_dir_bc_t;

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
  double              z_lith;
  double              z_asthen;
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

typedef struct slabs_surf_options
{
  slabs_x_func_t    x_func;
}
slabs_surf_options_t;

typedef struct slabs_topo_profile
{
  int               nsurf;
  double            *tX;
  double            *tY;
  double            *tZ;
}
slabs_topo_profile_t;

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


