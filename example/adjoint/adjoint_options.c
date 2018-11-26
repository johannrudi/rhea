#include <adjoint_options.h>

double              rayleigh = SUBD_RAYLEIGH;
double              ref_density = SUBD_REF_DENSITY;
double              ref_gravity = SUBD_REF_GRAVITY;
double              ref_therm_expa = SUBD_REF_THERM_EXPA;
double              ref_temp_diff = SUBD_REF_TEMP_DIFF;
double              ref_spatial = SUBD_REF_SPATIAL;
double              ref_visc = SUBD_REF_VISC;
double              ref_therm_diffus = SUBD_REF_THERM_DIFFUS;

int                 rhs_type;

  /* declare or initialize subd_temp options*/
int                 temp_type;

double              temp_back_plate_age =
  SUBD_TEMP_BACKGROUND_PLATE_AGE;
double              temp_2pl_trench_lon =
  SUBD_TEMP_2PL_TRENCH_LONGITUDE;
double              temp_2pl_dip_angle = SUBD_TEMP_2PL_DIP_ANGLE;
double              temp_2pl_subd_depth = SUBD_TEMP_2PL_SUBD_DEPTH;
double              temp_2pl_subd_width = SUBD_TEMP_2PL_SUBD_WIDTH;
double              temp_2pl_subd_edge_width =
  SUBD_TEMP_2PL_SUBD_EDGE_WIDTH;
double              temp_2pl_subd_edge_smoothwidth =
  SUBD_TEMP_2PL_SUBD_EDGE_SMOOTHWIDTH;
double              temp_2pl_subd_plate_vel =
  SUBD_TEMP_2PL_SUBD_PLATE_VELOCITY;
double              temp_2pl_subd_plate_init_age =
  SUBD_TEMP_2PL_SUBD_PLATE_INITIAL_AGE;
double              temp_2pl_over_plate_age =
  SUBD_TEMP_2PL_OVER_PLATE_AGE;

int                 temp_custom;
subd_temp_2plates_slab_t   slab_options;
subd_temp_custom_sinker_t  sinker_options;
subd_temp_custom_lithblock_t  lithblock_options;
subd_temp_custom_hscm2block_t  hscm2block_options;

/*subd_temp: sinker_options*/
double              center_x;
double              center_y;
double              center_z;
double              diameter;
double              sinker_scaling;
double              sinker_decay;

/*subd_temp: lithblock_options*/
double              dist_r;
double              dist_deg;

/*subd_temp: hscm2block_options*/
double              plate1_age_yr;
double              plate2_age_yr;
double              plates12_bound;

/* initialize subd options: weak zone */
int                     weakzone_type;

double              weakzone_2pl_subdu_lon =
  SUBD_WEAK_2PLATES_SUBDU_LONGITUDE;
double              weakzone_2pl_subdu_dip_angle =
  SUBD_WEAK_2PLATES_SUBDU_DIP_ANGLE;
double              weakzone_2pl_subdu_depth = SUBD_WEAK_2PLATES_SUBDU_DEPTH;
double              weakzone_2pl_subdu_width = SUBD_WEAK_2PLATES_SUBDU_WIDTH;
double              weakzone_2pl_subdu_thickness =
  SUBD_WEAK_2PLATES_SUBDU_THICKNESS;
double              weakzone_2pl_subdu_thickness_const =
  SUBD_WEAK_2PLATES_SUBDU_THICKNESS_CONST;
double              weakzone_2pl_subdu_weak_factor =
  SUBD_WEAK_2PLATES_SUBDU_WEAK_FACTOR;
double              weakzone_2pl_ridge_depth = SUBD_WEAK_2PLATES_RIDGE_DEPTH;
double              weakzone_2pl_ridge_width = SUBD_WEAK_2PLATES_RIDGE_WIDTH;
double              weakzone_2pl_ridge_smoothwidth =
  SUBD_WEAK_2PLATES_RIDGE_SMOOTHWIDTH;
double              weakzone_2pl_ridge_weak_factor =
  SUBD_WEAK_2PLATES_RIDGE_WEAK_FACTOR;

/* initialize subd options: viscosity */
int                     visc_anisotropy;
int                     visc_type;
int                     visc_custom;
double              visc_lith =
  SUBD_VISC_LITH;
double              visc_asthen =
  SUBD_VISC_ASTHEN;
double              visc_z_lith =
  SUBD_VISC_LITH_RADIUS_LOCATION;
double              visc_z_asthen =
  SUBD_VISC_ASTHEN_RADIUS_LOCATION;
double              visc_z_slab =
  SUBD_VISC_SLAB_RADIUS_LOCATION;
double              visc_slab_width =
  SUBD_VISC_SLAB_WIDTH;

/* initialize subd options: velocity boundary condition */
double              flow_scale =
  SUBD_VELBC_COLLIDE_FLOW_SCALE;
double              velocity_bc_middle =
  SUBD_VELBC_COLLIDE_ZERO_POINT_LOCATION_MIDDLE;
double              velocity_bc_upper =
  SUBD_VELBC_COLLIDE_ZERO_POINT_LOCATION_UPPER;
double              velocity_bc_lower =
  SUBD_VELBC_COLLIDE_ZERO_POINT_LOCATION_LOWER;

double                  surf_dist;
double                  visc_trial;
int                     x_func;
int                     velbc;
int                     velbc_nonzero_dir;
int                     velbc_nonzero_neu;
char                    *txt_read_surfvelo_path;

int                     test_manufactured;
int                     test_stress_op;
int                     test_stress_comp;

int                     adjoint_visc_n_components;
int                     adjoint_visc_type;
int                     adjoint_stencil_visc_type;
int                     adjoint_stencil_visc_custom_type;
double                  adjoint_stencil_visc_value;
subd_adjoint_stencil_options_t stencil_options;

void
subduction_add_options (ymir_options_t * opt)
{
  ymir_options_addv (opt,

//  /* physics parameters */
  YMIR_OPTIONS_D, "rayleigh", '\0',
    &rayleigh, SUBD_RAYLEIGH,
    "rayleigh number",
  YMIR_OPTIONS_D, "reference-density", '\0',
    &ref_density, SUBD_REF_DENSITY,
    "reference density",
  YMIR_OPTIONS_D, "reference-gravity", '\0',
    &ref_gravity, SUBD_REF_GRAVITY,
    "reference gravitational accelaration",
  YMIR_OPTIONS_D, "reference-thermal-expansion", '\0',
    &ref_therm_expa, SUBD_REF_THERM_EXPA,
    "reference:thermal expansion",
  YMIR_OPTIONS_D, "temperature-difference", '\0',
    &ref_temp_diff, SUBD_REF_TEMP_DIFF,
    "reference:temperature difference across the mantle",
  YMIR_OPTIONS_D, "reference-spatial", '\0',
    &ref_spatial, SUBD_REF_SPATIAL,
    "spatial reference, radius_max or radius_diff",
  YMIR_OPTIONS_D, "reference-viscosity", '\0',
    &ref_visc, SUBD_REF_VISC,
    "reference viscosity",
  YMIR_OPTIONS_D, "reference-thermal-diffusion", '\0',
    &ref_therm_diffus, SUBD_REF_THERM_DIFFUS,
    "reference thermal diffusion",

  /* right-hand-side */
  YMIR_OPTIONS_I, "rhs-type",'\0',
    &(rhs_type),SUBD_RHS_BUOY,
    "0: density anomaly,  2: full density, 3: manufactured test",

  /* temperature */
  YMIR_OPTIONS_I, "temp-type",'\0',
    &(temp_type),SUBD_TEMP_RHEA,
    "0: sinker, 1: 2plates_poly2 subd, 2: collide",

  YMIR_OPTIONS_D, "temp-background-plate-age", '\0',
    &temp_back_plate_age, SUBD_TEMP_BACKGROUND_PLATE_AGE,
    "Bachground temperature descibed by plate age [yr]",
  YMIR_OPTIONS_D, "temp-2plates-trench-longitude", '\0',
    &temp_2pl_trench_lon, SUBD_TEMP_2PL_TRENCH_LONGITUDE,
    "2plates temp: Longitude of trench in interval (-pi/8, pi/8)",
  YMIR_OPTIONS_D, "temp-2plates-dip-angle", '\0',
    &temp_2pl_dip_angle, SUBD_TEMP_2PL_DIP_ANGLE,
    "2plates temp: Dip angle of subducting plate (in degrees < 0)",
  YMIR_OPTIONS_D, "temp-2plates-subd-depth", '\0',
    &temp_2pl_subd_depth, SUBD_TEMP_2PL_SUBD_DEPTH,
    "2plates temp: Maximal depth of subducting plate inside of mantle [m]",
  YMIR_OPTIONS_D, "temp-2plates-subd-width", '\0',
    &temp_2pl_subd_width, SUBD_TEMP_2PL_SUBD_WIDTH,
    "2plates temp: Maximal width of subducting zone inside of mantle [m]",
  YMIR_OPTIONS_D, "temp-2plates-subd-edge-width", '\0',
    &temp_2pl_subd_edge_width, SUBD_TEMP_2PL_SUBD_EDGE_WIDTH,
    "2plates temp: Width for subducting plate's top edge [m]",
  YMIR_OPTIONS_D, "temp-2plates-subd-edge-smoothwidth", '\0',
    &temp_2pl_subd_edge_smoothwidth,
    SUBD_TEMP_2PL_SUBD_EDGE_SMOOTHWIDTH,
    "2plates weak zone: Width of smoothing of subd. plate's top edge [m]",
  YMIR_OPTIONS_D, "temp-2plates-subd-plate-velocity", '\0',
    &temp_2pl_subd_plate_vel, SUBD_TEMP_2PL_SUBD_PLATE_VELOCITY,
    "2plates temp: Velocity of subducting plate [m/y]",
  YMIR_OPTIONS_D, "temp-2plates-subd-plate-initial-age", '\0',
    &temp_2pl_subd_plate_init_age,
    SUBD_TEMP_2PL_SUBD_PLATE_INITIAL_AGE,
    "2plates temp: Age of subducting plate at left boundary [y]",
  YMIR_OPTIONS_D, "temp-2plates-over-plate-age", '\0',
    &temp_2pl_over_plate_age, SUBD_TEMP_2PL_OVER_PLATE_AGE,
    "2plates temp: Age of overriding plate [y]",

  YMIR_OPTIONS_I, "temp-custom-type",'\0',
    &(temp_custom),SUBD_TEMP_CUSTOM_SINKER,
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

  YMIR_OPTIONS_D, "temp-lithblock-dist-radius", '\0',
    &dist_r, 0.001,
    "cold lithosphere depth in 0~0.45",
  YMIR_OPTIONS_D, "temp-lithblock-dist-degree", '\0',
    &dist_deg, 10.0,
    "cold lithosphere distance from right in degree",

  YMIR_OPTIONS_D, "plate1-age-yr", '\0',
    &plate1_age_yr, 5.0e7,
    "the age in year of the lithosphere on the left in half-space-cooling-model",
  YMIR_OPTIONS_D, "plate2-age-yr", '\0',
    &plate2_age_yr, 1.0e8,
    "the age in year of the lithosphere on the right in half-space-cooling-model",
  YMIR_OPTIONS_D, "plate1-plate2-boundary-x", '\0',
    &plates12_bound, 1.0,
    "the boundary (nondimensional) of plates 1 and 2 at x direction",

  /* weakzone */
  YMIR_OPTIONS_I, "weakzone-type",'\0',
    &(weakzone_type), SUBD_WEAK_NONE,
    "0: none, 1: poly2 slab",
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-longitude", '\0',
    &weakzone_2pl_subdu_lon, SUBD_WEAK_2PLATES_SUBDU_LONGITUDE,
    "2plates weak zone: Longitude in interval (-pi/8, pi/8), "
    "where weak zone begins",
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-dip-angle", '\0',
    &weakzone_2pl_subdu_dip_angle, SUBD_WEAK_2PLATES_SUBDU_DIP_ANGLE,
    "2plates weak zone: Dip angle of weak zone (in degrees > 0)",
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-depth", '\0',
    &weakzone_2pl_subdu_depth, SUBD_WEAK_2PLATES_SUBDU_DEPTH,
    "2plates weak zone: Depth of weak zone [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-width", '\0',
    &weakzone_2pl_subdu_width, SUBD_WEAK_2PLATES_SUBDU_WIDTH,
    "2plates weak zone: Width of weak zone [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-thickness", '\0',
    &weakzone_2pl_subdu_thickness, SUBD_WEAK_2PLATES_SUBDU_THICKNESS,
    "2plates weak zone: Width at center of weak zone [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-thickness-const", '\0',
    &weakzone_2pl_subdu_thickness_const,
    SUBD_WEAK_2PLATES_SUBDU_THICKNESS_CONST,
    "2plates weak zone: Width of smoothing of edges of weak zone [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-weak-factor", '\0',
    &weakzone_2pl_subdu_weak_factor, SUBD_WEAK_2PLATES_SUBDU_WEAK_FACTOR,
    "2plates weak zone: Value of weak zone factor",
  YMIR_OPTIONS_D, "weakzone-2plates-ridge-depth", '\0',
    &weakzone_2pl_ridge_depth, SUBD_WEAK_2PLATES_RIDGE_DEPTH,
    "2plates weak zone: Depth of weak zone in left corner of domain [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-ridge-width", '\0',
    &weakzone_2pl_ridge_width, SUBD_WEAK_2PLATES_RIDGE_WIDTH,
    "2plates weak zone: Width of weak zone in left corner of domain [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-ridge-smoothwidth", '\0',
    &weakzone_2pl_ridge_smoothwidth, SUBD_WEAK_2PLATES_RIDGE_SMOOTHWIDTH,
    "2plates weak zone: Smoothing width of edges of weak zone in corner [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-ridge-weak-factor", '\0',
    &weakzone_2pl_ridge_weak_factor, SUBD_WEAK_2PLATES_RIDGE_WEAK_FACTOR,
    "2plates weak zone: Value of weak zone factor for weak zone in corner",

  /* viscosity */
  YMIR_OPTIONS_I, "viscosity-anisotropy",'\0',
    &(visc_anisotropy),SUBD_VISC_ISOTROPY,
    "0: isotropy, 1: transversely isotropy",
  YMIR_OPTIONS_I, "viscosity-type",'\0',
    &(visc_type),SUBD_VISC_RHEA,
    "0: rhea, 1: custom",
  YMIR_OPTIONS_I, "custom-viscosity",'\0',
    &(visc_custom),SUBD_VISC_CUSTOM_LAYERS,
    "1: layers, 2: layers-coupling",
  YMIR_OPTIONS_D, "viscosity-lithosphere", '\0',
    &visc_lith, SUBD_VISC_LITH,
    "Viscosity in the lithosphere",
  YMIR_OPTIONS_D, "viscosity-lith-trial", '\0',
    &visc_trial, 1.0,
    "Viscosity in the lithosphere",
  YMIR_OPTIONS_D, "viscosity-asthenosphere", '\0',
    &visc_asthen, SUBD_VISC_ASTHEN,
    "Viscosity in the asthenosphere",
  YMIR_OPTIONS_D, "visc-zlocation-lithosphere", '\0',
    &visc_z_lith, SUBD_VISC_LITH_RADIUS_LOCATION,
    "Viscosity lithosphere location",
  YMIR_OPTIONS_D, "visc-zlocation-asthenosphere", '\0',
    &visc_z_asthen, SUBD_VISC_ASTHEN_RADIUS_LOCATION,
    "Viscosity asthenosphere location",
  YMIR_OPTIONS_D, "visc-zlocation-slab", '\0',
    &visc_z_slab, SUBD_VISC_SLAB_RADIUS_LOCATION,
    "Viscosity slab location",
  YMIR_OPTIONS_D, "visc-slab-width", '\0',
    &visc_slab_width, SUBD_VISC_SLAB_WIDTH,
    "Viscosity slab width",

  /* surface location  */
  YMIR_OPTIONS_I, "bound-x-function", '\0',
    &x_func, SUBD_X_FUNCTION_IDENTITY,
    "boundary location: surface topography",
  YMIR_OPTIONS_D, "surface-distortion-factor", '\0',
    &surf_dist, 1.0,
    "surface distortion factor.",

  /* velocity Dirichlet BC's */
  YMIR_OPTIONS_I, "velocity-bc", '\0',
    &velbc, SUBD_FREESURFACE,
    "Velocity boundary condition",

  YMIR_OPTIONS_I, "velocity-nonzero-dirichlet", '\0',
    &velbc_nonzero_dir, SUBD_VELBC_DIR_ZERO,
    "Nonezero Dirichlet boundary condition",
  YMIR_OPTIONS_D, "flow-scaling", '\0',
    &flow_scale, SUBD_VELBC_COLLIDE_FLOW_SCALE,
    "scaling of velocity BC.",
  YMIR_OPTIONS_D, "velocity-bc-middle", '\0',
    &velocity_bc_middle, SUBD_VELBC_COLLIDE_ZERO_POINT_LOCATION_MIDDLE,
    "location of velocity BC: middle bound",
  YMIR_OPTIONS_D, "velocity-bc-upper", '\0',
    &velocity_bc_upper, SUBD_VELBC_COLLIDE_ZERO_POINT_LOCATION_UPPER,
    "location of velocity BC: upper bound",
  YMIR_OPTIONS_D, "velocity-bc-lower", '\0',
    &velocity_bc_lower, SUBD_VELBC_COLLIDE_ZERO_POINT_LOCATION_LOWER,
    "location of velocity BC: lower bound",

  YMIR_OPTIONS_I, "velocity-nonzero-neumann", '\0',
    &velbc_nonzero_neu, SUBD_VELBC_NEU_ZERO,
    "Nonezero Neumann boundary condition",
  YMIR_OPTIONS_S, "txt-read-surfvelo-path", '\0',
      &(txt_read_surfvelo_path), NULL,
          "File path for writing txt files of the surface velocity",

  /* test options: stress operator test and manufactured solution test */
  YMIR_OPTIONS_I, "test-stress-operator", '\0',
    &test_stress_op, SUBD_TEST_STRESS_OP_NONE,
    "the input velocity for stress operator test",
  YMIR_OPTIONS_I, "test-manufactured-solution", '\0',
    &test_manufactured, SUBD_TEST_MANUFACTURED_NONE,
    "the input for velocity and pressure field for manufactured solution test",
  YMIR_OPTIONS_I, "test-stress-component", '\0',
    &test_stress_comp, SUBD_TEST_MANUFACTURED_NONE,
    "the input for velocity and pressure field for manufactured solution test",

  YMIR_OPTIONS_I, "adjoint-visc-n-components", '\0',
    &adjoint_visc_n_components, 1,
    "the number of unknown parameters in the inversion",

  YMIR_OPTIONS_I, "adjoint-visc-type", '\0',
    &adjoint_visc_type, SUBD_ADJOINT_VISC_PRE,
    "which viscosity parameter in adjoint problem, 0: prefactor; 1: activation energy",

  YMIR_OPTIONS_I, "adjoint-stencil-visc-type", '\0',
    &adjoint_stencil_visc_type, SUBD_ADJOINT_STENCIL_VISC_CUSTOM,
    "stencil of viscosity in adjoint problem, 0: rhea; 1: custom",

  YMIR_OPTIONS_I, "adjoint-stencil-visc-custom-type", '\0',
    &adjoint_stencil_visc_custom_type, SUBD_ADJOINT_STENCIL_VISC_CUSTOM_LAYERS_UM,
    "stencil of customer defined viscosity in adjoint problem, 1: layers; 2: layers-coupling",

  YMIR_OPTIONS_D, "adjoint-stencil-visc-value", '\0',
    &adjoint_stencil_visc_value, 1.0,
    "value of viscosity stencil  in adjoint problem",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

}

void
subduction_process_options (subd_options_t *subd_options,
                            rhea_domain_options_t  *domain_options,
                            rhea_viscosity_options_t  *rhea_visc_options,
                            subd_para_options_t  *para_options,
                            subd_temp_options_t  *temp_options,
                            subd_visc_options_t  *visc_options,
                            subd_weak_options_t  *weak_options,
                            subd_surf_options_t  *surf_options,
                            subd_velbc_options_t  *velbc_options,
                            subd_test_options_t  *test_options,
                            subd_adjoint_options_t *adjoint_options)
{
  para_options->rayleigh = rayleigh;
  para_options->ref_density = ref_density;
  para_options->ref_gravity = ref_gravity;
  para_options->ref_therm_expa = ref_therm_expa;
  para_options->ref_temp_diff = ref_temp_diff;
  para_options->ref_spatial = ref_spatial;
  para_options->ref_visc = ref_visc;
  para_options->ref_therm_diffus = ref_therm_diffus;

  /* temperature */
  if (test_manufactured)
    temp_type = 4;
  temp_options->type = (subd_temp_type_t) temp_type;

  slab_options.temp_background_plate_age = temp_back_plate_age;
  slab_options.temp_2plates_trench_longitude = temp_2pl_trench_lon;
  slab_options.temp_2plates_dip_angle = temp_2pl_dip_angle;
  slab_options.temp_2plates_subd_depth = temp_2pl_subd_depth;
  slab_options.temp_2plates_subd_width = temp_2pl_subd_width;
  slab_options.temp_2plates_subd_edge_width = temp_2pl_subd_edge_width;
  slab_options.temp_2plates_subd_edge_smoothwidth =
    temp_2pl_subd_edge_smoothwidth;
  slab_options.temp_2plates_subd_plate_velocity = temp_2pl_subd_plate_vel;
  slab_options.temp_2plates_subd_plate_initial_age =
    temp_2pl_subd_plate_init_age;
  slab_options.temp_2plates_over_plate_age = temp_2pl_over_plate_age;
  temp_options->slab_options = &slab_options;

  temp_options->custom_type = (subd_temp_custom_t) temp_custom;

  sinker_options.center_x = center_x;
  sinker_options.center_y = center_y;
  sinker_options.center_z = center_z;
  sinker_options.diameter = diameter;
  sinker_options.scaling = sinker_scaling;
  sinker_options.decay = sinker_decay;
  temp_options->sinker_options = &sinker_options;

  lithblock_options.dist_r = dist_r;
  lithblock_options.dist_deg = dist_deg;
  temp_options->lithblock_options = &lithblock_options;

  hscm2block_options.plate1_age_yr = plate1_age_yr;
  hscm2block_options.plate2_age_yr = plate2_age_yr;
  hscm2block_options.plates12_bound = plates12_bound;
  temp_options->hscm2block_options = &hscm2block_options;

  /* weak zone */
  weak_options->type = (subd_weakzone_type_t) weakzone_type;
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
  visc_options->anisotropy_type = (subd_visc_anisotropy_t) visc_anisotropy;
  visc_options->type = (subd_visc_type_t) visc_type;
  visc_options->custom_type = (subd_visc_custom_t) visc_custom;
  visc_options->visc_lith = visc_lith;
  visc_options->visc_asthen = visc_asthen;
  visc_options->z_lith = visc_z_lith;
  visc_options->z_asthen = visc_z_asthen;
  visc_options->z_slab = visc_z_slab;
  visc_options->slab_width = visc_slab_width;

  /*geometry transformation: surface location*/
  surf_options->x_func = (subd_x_func_t) x_func;

  /* velocity B.C.-> condition */
  velbc_options->velbc = (subd_velbc_t) velbc;
  velbc_options->velbc_nonzero_dir = (subd_velbc_nonzero_dir_t) velbc_nonzero_dir;
  velbc_options->flow_scale = flow_scale;
  velbc_options->vel_dir_bc_middle = velocity_bc_middle;
  velbc_options->vel_dir_bc_upper = velocity_bc_upper;
  velbc_options->vel_dir_bc_lower = velocity_bc_lower;

  velbc_options->velbc_nonzero_neu = (subd_velbc_nonzero_neu_t) velbc_nonzero_neu;
  velbc_options->txt_read_path_nonzero_neu = txt_read_surfvelo_path;

  /* test */
  test_options->test_stress_op = (subd_test_stress_op_t) test_stress_op;
  test_options->test_manufactured = (subd_test_manufactured_t) test_manufactured;
  test_options->test_stress_comp = (subd_test_manufactured_t) test_stress_comp;

  adjoint_options->n_components = adjoint_visc_n_components;
  adjoint_options->visc_type = (subd_adjoint_visc_type_t) adjoint_visc_type;
  stencil_options.type = (subd_adjoint_stencil_visc_type_t) adjoint_stencil_visc_type;
  stencil_options.custom_type = (subd_adjoint_stencil_visc_custom_t) adjoint_stencil_visc_custom_type;
  stencil_options.value = adjoint_stencil_visc_value;
  adjoint_options->stencil_options = &stencil_options;

  /* assign subd_options */
  subd_options->rhs_type = (subd_rhs_type_t) rhs_type;
  subd_options->para_options = para_options;
  subd_options->temp_options = temp_options;
  subd_options->visc_options = visc_options;
  subd_options->weak_options = weak_options;
  subd_options->surf_options = surf_options;
  subd_options->velbc_options = velbc_options;
  subd_options->test_options = test_options;
  subd_options->adjoint_options = adjoint_options;

  subd_options->domain_options = domain_options;
  subd_options->rhea_visc_options = rhea_visc_options;
}


