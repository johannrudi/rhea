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
#define SLABS_VISCOSITY_LITHOSPHERE (100.0)
#define SLABS_VISCOSITY_ASTHENOSPHERE (1.0)
#define SLABS_VISCOSITY_LITHOSPHERE_RADIUS_LOCATION (0.9)
#define SLABS_VISCOSITY_ASTHENOSPHERE_RADIUS_LOCATION (0.2)

/* Dirichlet velocity B.C. */
#define SLABS_COLLIDE_FLOW_SCALE (10.0)
#define SLABS_COLLIDE_ZERO_POINT_LOCATION_MIDDLE (0.5)
#define SLABS_COLLIDE_ZERO_POINT_LOCATION_UPPER (0.5)
#define SLABS_COLLIDE_ZERO_POINT_LOCATION_LOWER (0.5)

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
  SLABS_VISCOSITY_LITHOSPHERE;
double              visc_asthen =
  SLABS_VISCOSITY_ASTHENOSPHERE;
double              visc_z_lith =
  SLABS_VISCOSITY_LITHOSPHERE_RADIUS_LOCATION;
double              visc_z_asthen =
  SLABS_VISCOSITY_ASTHENOSPHERE_RADIUS_LOCATION;

/* initialize slabs options: velocity boundary condition */
double              flow_scale =
  SLABS_COLLIDE_FLOW_SCALE;
double              velocity_bc_middle =
  SLABS_COLLIDE_ZERO_POINT_LOCATION_MIDDLE;
double              velocity_bc_upper =
  SLABS_COLLIDE_ZERO_POINT_LOCATION_UPPER;
double              velocity_bc_lower =
  SLABS_COLLIDE_ZERO_POINT_LOCATION_LOWER;

typedef enum
{
  SINKER,
  SLAB,
  COLLIDE,
  DRAG,
  TEST_MANUFACTURED,
  TESTTOPO,
  TESTTOPO2,
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

/**
 *  * Defines options and adds them as sub-options.
 *   */
void
slabweakzone_options (ymir_options_t * opt)
{

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
  YMIR_OPTIONS_D, "viscosity-lithosphere", '\0',
    &visc_lith, SLABS_VISCOSITY_LITHOSPHERE,
    "Viscosity in the lithosphere",
  YMIR_OPTIONS_D, "viscosity-asthenosphere", '\0',
    &visc_asthen, SLABS_VISCOSITY_ASTHENOSPHERE,
    "Viscosity in the asthenosphere",
  YMIR_OPTIONS_D, "visc-zlocation-lithosphere", '\0',
    &visc_z_lith, SLABS_VISCOSITY_LITHOSPHERE_RADIUS_LOCATION,
    "Viscosity in the lithosphere",
  YMIR_OPTIONS_D, "visc-zlocation-asthenosphere", '\0',
    &visc_z_asthen, SLABS_VISCOSITY_ASTHENOSPHERE_RADIUS_LOCATION,
    "Viscosity in the asthenosphere",

  /* surface location  */
  YMIR_OPTIONS_I, "bound-x-function", '\0',
    &x_func, SLABS_X_FUNCTION_IDENTITY,
    "boundary location: surface topography",

  /* velocity Dirichlet BC's */
  YMIR_OPTIONS_I, "velocity-dirichlet-bc", '\0',
    &vel_dir_bc, SLABS_VEL_DIR_BC_INOUTFLOW_SIN,
    "Velocity Dirichlet boundary condition",
  YMIR_OPTIONS_D, "flow-scaling", '\0',
    &flow_scale, SLABS_COLLIDE_FLOW_SCALE,
    "scaling of velocity BC.",
  YMIR_OPTIONS_D, "velocity-bc-middle", '\0',
    &velocity_bc_middle, SLABS_COLLIDE_ZERO_POINT_LOCATION_MIDDLE,
    "location of velocity BC: middle bound",
  YMIR_OPTIONS_D, "velocity-bc-upper", '\0',
    &velocity_bc_upper, SLABS_COLLIDE_ZERO_POINT_LOCATION_UPPER,
    "location of velocity BC: upper bound",
  YMIR_OPTIONS_D, "velocity-bc-lower", '\0',
    &velocity_bc_lower, SLABS_COLLIDE_ZERO_POINT_LOCATION_LOWER,
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
  YMIR_OPTIONS_S, "vtk-write-input3-path", '\0',
    &(vtk_write_input3_path), NULL,
    "File path for vtk files for the input of the Stokes problem",
  YMIR_OPTIONS_S, "vtk-write-solution3-path", '\0',
    &(vtk_write_solution3_path), NULL,
    "File path for vtk files for the solution of the Stokes problem",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add sub-options */
  rhea_add_options_all (opt);
  ymir_options_add_suboptions_solver_stokes (opt);
}
  /*
   * Process Slabs Options
   */
slabweakzone_options_process (slabs_options_t *slabs_options,
                              slabs_temp_options_t *temp_options,
                              slabs_visc_options_t *visc_options,
                              slabs_weak_options_t *weak_options,
                              slabs_surf_options_t *surf_options,
                              slabs_velbc_options_t *velbc_options,
                              slabs_test_options_t *test_options)
{
  /* temperature */
  temp_options.temp_background_plate_age = temp_back_plate_age;
  temp_options->temp_2plates_trench_longitude = temp_2pl_trench_lon;
  temp_options->temp_2plates_dip_angle = temp_2pl_dip_angle;
  temp_options->temp_2plates_subd_depth = temp_2pl_subd_depth;
  temp_options->temp_2plates_subd_width = temp_2pl_subd_width;
  temp_options->temp_2plates_subd_edge_width = temp_2pl_subd_edge_width;
  temp_options->temp_2plates_subd_edge_smoothwidth =
    temp_2pl_subd_edge_smoothwidth;
  temp_options->temp_2plates_subd_plate_velocity = temp_2pl_subd_plate_vel;
  temp_options->temp_2plates_subd_plate_initial_age =
    temp_2pl_subd_plate_init_age;
  temp_options->temp_2plates_over_plate_age = temp_2pl_over_plate_age;

  /* weak zone */
  weak_options->weakzone_2plates_subdu_longitude =
    weakzone_2pl_subdu_lon;
  weak_options->weakzone_2plates_subdu_dip_angle =
    weakzone_2pl_subdu_dip_angle;
  weak_options->weakzone_2plates_subdu_depth =
    weakzone_2pl_subdu_depth;
  weak_options->weakzone_2plates_subdu_width =
    weakzone_2pl_subdu_width;
  weak_options->weakzone_2plates_subdu_thickness =
    weakzone_2pl_subdu_thickness;
  weak_options->weakzone_2plates_subdu_thickness_const =
    weakzone_2pl_subdu_thickness_const;
  weak_options->weakzone_2plates_subdu_weak_factor =
    weakzone_2pl_subdu_weak_factor;
  weak_options->weakzone_2plates_ridge_depth =
    weakzone_2pl_ridge_depth;
  weak_options->weakzone_2plates_ridge_width =
    weakzone_2pl_ridge_width;
  weak_options->weakzone_2plates_ridge_smoothwidth =
    weakzone_2pl_ridge_smoothwidth;
  weak_options->weakzone_2plates_ridge_weak_factor =
    weakzone_2pl_ridge_weak_factor;

  /* viscosity */
  visc_options->viscosity_anisotropy = (slabs_viscosity_anisotropy_t) viscosity_anisotropy;
  visc_options->visc_lith = visc_lith;
  visc_options->visc_asthen = visc_asthen;
  visc_options->z_lith = visc_z_lith;
  visc_options->z_asthen = visc_z_asthen;

  /*geometry transformation: surface location*/
  surf_options->x_func = (slabs_x_func_t) x_func;

  /* velocity B.C. condition */
  velbc_options->vel_dir_bc = (slabs_vel_dir_bc_t) vel_dir_bc;
  velbc_options->flow_scale = flow_scale;
  velbc_options->vel_dir_bc_middle = velocity_bc_middle;
  velbc_options->vel_dir_bc_upper = velocity_bc_upper;
  velbc_options->vel_dir_bc_lower = velocity_bc_lower;

  /* test */
  test_options->test_stress_op = (slabs_test_stress_op_t) test_stress_op;
  test_options->test_manufactured = (slabs_test_manufactured_t) test_manufactured;
  test_options->test_stress_comp = (slabs_test_manufactured_t) test_stress_comp;

  /* assign slabs_options */
  slabs_options->slabs_temp_options = temp_options;
  slabs_options->slabs_visc_options = visc_options;
  slabs_options->slabs_weak_options = weak_options;
  slabs_options->slabs_surf_options = surf_options;
  slabs_options->slabs_velbc_options = velbc_options;
  slabs_options->slabs_test_options = test_options;

  /* buoyancy type */
  if (test_manufactured)
    buoyancy_type = 4;
  slabs_options->buoyancy_type = (slabs_buoyancy_type_t) buoyancy_type;

}

/* copy rhea domain options into local example domain options */
slabweakzone_options_domain_opt (slabs_options_t *slabs_options,
                                 slabs_domain_options_t *slabs_domain_options,
                                 rhea_domain_options_t  *domain_options)
{
  slabs_domain_options->x_min = domain_options->x_min;
  slabs_domain_options->x_max = domain_options->x_max;
  slabs_domain_options->y_min = domain_options->y_min;
  slabs_domain_options->y_max = domain_options->y_max;
  slabs_domain_options->z_min = domain_options->z_min;
  slabs_domain_options->z_max = domain_options->z_max;
  slabs_domain_options->lon_min = domain_options->lon_min;
  slabs_domain_options->lon_max = domain_options->lon_max;
  slabs_domain_options->radius_min = domain_options->radius_min;
  slabs_domain_options->radius_max = domain_options->radius_max;
  slabs_options->slabs_domain_options = &slabs_domain_options;

}
