#ifndef ADJOINT_OPTIONS_H
#define ADJOINT_OPTIONS_H

#include <rhea.h>

/* basic constants */
#define SUBD_RAYLEIGH (1.0)
#define SUBD_REF_DENSITY (3300.0)  /* (kg/m^3) */
#define SUBD_REF_GRAVITY (9.8)  /* (m/s^2) */
#define SUBD_REF_THERM_EXPA (3.0e-5)
#define SUBD_REF_TEMP_DIFF (3000.0) /* K */
#define SUBD_REF_SPATIAL (3500.0) /* r_max-r_min (km) */
#define SUBD_REF_VISC (1.0e20) /* Pa s */
#define SUBD_REF_THERM_DIFFUS (1.0e-6)        /* thermal diffusivity [m^2 / s]      */
#define SUBD_REF_VISC (1.0e20)            /* representative viscosity [Pa s]    */

#define SUBD_SEC_PER_YEAR (31557600.0)    /* seconds in a year (365.25*24*3600) */
#define SUBD_EARTH_RADIUS (6371.0e3)      /* mean radius of the Earth [m]       */
#define SUBD_MANTLE_DEPTH (2871.0e3)      /* mean radius of the Earth [m]       */
#define SUBD_UPPER_MANTLE_DEPTH (660.0e3) /* approx. depth of upper mantle [m]  */

/* temperature parameters */
#define SUBD_CONST_TEMP 0.5

/* shell parameters */
#define SUBD_SHELL_RADIUS_BOTTOM 0.55
#define SUBD_SHELL_RADIUS_TOP 1.0

/* temperature 2plates_poly2 */
#define SUBD_TEMP_BACKGROUND_PLATE_AGE (50.0e6)    /* rhea1: 60 Myr */
#define SUBD_TEMP_2PL_TRENCH_LONGITUDE (0.13)      /* rhea1: 0.13   */
#define SUBD_TEMP_2PL_DIP_ANGLE (5.0)              /* rhea1: --     */
#define SUBD_TEMP_2PL_SUBD_DEPTH (400.0e3)         /* rhea1: 400 km */
#define SUBD_TEMP_2PL_SUBD_WIDTH (300.0e3)         /* rhea1: --     */
#define SUBD_TEMP_2PL_SUBD_EDGE_WIDTH (1.0e3)      /* rhea1: --     */
#define SUBD_TEMP_2PL_SUBD_EDGE_SMOOTHWIDTH (40.0e3)/* rhea1: --     */
#define SUBD_TEMP_2PL_SUBD_PLATE_VELOCITY (4.0e-2)  /* rhea1: 4 cm/y */
#define SUBD_TEMP_2PL_SUBD_PLATE_INITIAL_AGE (1.0e6)/* rhea1: --     */
#define SUBD_TEMP_2PL_OVER_PLATE_AGE (40.0e6)       /* rhea1: 40 Myr */

/* weak zone 2plates_poly2 */
#define SUBD_WEAK_2PLATES_SUBDU_LONGITUDE (-100.0)      /* rhea1:  0.13 */
#define SUBD_WEAK_2PLATES_SUBDU_DIP_ANGLE (-1.0)        /* rhea1:   N/A */
#define SUBD_WEAK_2PLATES_SUBDU_DEPTH (80.0e3)          /* rhea1: 50 km */
#define SUBD_WEAK_2PLATES_SUBDU_WIDTH (-1.0)            /* rhea1:   N/A */
#define SUBD_WEAK_2PLATES_SUBDU_THICKNESS (20.0e3)      /* rhea1: 20 km */
#define SUBD_WEAK_2PLATES_SUBDU_THICKNESS_CONST (5.0e3) /* rhea1: 10 km */
#define SUBD_WEAK_2PLATES_SUBDU_WEAK_FACTOR (1.0e-5)    /* rhea1:  1e-5 */
#define SUBD_WEAK_2PLATES_RIDGE_DEPTH (30.0e3)          /* rhea1: 30 km */
#define SUBD_WEAK_2PLATES_RIDGE_WIDTH (30.0e3)          /* rhea1: 30 km */
#define SUBD_WEAK_2PLATES_RIDGE_SMOOTHWIDTH (10.0e3)    /* rhea1:  5 km */
#define SUBD_WEAK_2PLATES_RIDGE_WEAK_FACTOR (1.0e-5)    /* rhea1:  1e-5 */

/* viscosity */
#define SUBD_VISC_LITH (100.0)
#define SUBD_VISC_ASTHEN (1.0)
#define SUBD_VISC_LITH_RADIUS_LOCATION (0.9)
#define SUBD_VISC_ASTHEN_RADIUS_LOCATION (0.2)
#define SUBD_VISC_SLAB_RADIUS_LOCATION (0.7)
#define SUBD_VISC_SLAB_WIDTH (0.2)

/* Dirichlet velocity B.C. */
#define SUBD_VELBC_COLLIDE_FLOW_SCALE (10.0)
#define SUBD_VELBC_COLLIDE_ZERO_POINT_LOCATION_MIDDLE (0.5)
#define SUBD_VELBC_COLLIDE_ZERO_POINT_LOCATION_UPPER (0.5)
#define SUBD_VELBC_COLLIDE_ZERO_POINT_LOCATION_LOWER (0.5)


typedef struct subd_2plates_poly2_geo_coeff
{
  double        start_node;
  double        start_val;
  double        end_node;
  double        end_val;
  double        start_deriv;
  double        *poly2_coeff;
}
subd_2plates_poly2_geo_coeff_t;

typedef struct subd_para_options
{
  double    rayleigh;
  double    ref_density;
  double    ref_gravity;
  double    ref_therm_expa;
  double    ref_temp_diff;
  double    ref_spatial;
  double    ref_visc;
  double    ref_therm_diffus;
}
subd_para_options_t;

/*
 * right-hand-side options
 */
typedef enum
{
  SUBD_RHS_BUOY,
  SUBD_RHS_DENSITY,
  SUBD_RHS_ADJOINT_NONE,
  SUBD_RHS_TEST_MANUFACTURED
}
subd_rhs_type_t;

/*temperature options*/
typedef enum
{
  SUBD_TEMP_RHEA,
  SUBD_TEMP_SLAB,
  SUBD_TEMP_CUSTOM,
  SUBD_TEMP_NONE
}
subd_temp_type_t;

typedef struct subd_temp_2plates_slab
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
  subd_2plates_poly2_geo_coeff_t *temp_2plates_geo_coeff;
}
subd_temp_2plates_slab_t;

typedef enum
{
  SUBD_TEMP_CUSTOM_NONE,
  SUBD_TEMP_CUSTOM_THINBOX = 1,
  SUBD_TEMP_CUSTOM_SINKER = 2,
  SUBD_TEMP_CUSTOM_DRAG = 3,
  SUBD_TEMP_CUSTOM_LITHBLOCK = 4,
  SUBD_TEMP_CUSTOM_HSCM2BLOCK = 5,
  SUBD_TEMP_CUSTOM_HSCM2BLOCKTANH = 6,
  SUBD_TEMP_CUSTOM_HSCM1PLATE = 7
}
subd_temp_custom_t;

typedef struct subd_temp_custom_sinker
{
  double center_x;
  double center_y;
  double center_z;
  double diameter;
  double scaling;
  double decay;
}
subd_temp_custom_sinker_t;

typedef struct subd_temp_custom_lithblock
{
  double dist_r;
  double dist_deg;
}
subd_temp_custom_lithblock_t;

typedef struct subd_temp_custom_hscm2block
{
  double plate1_age_yr;
  double plate2_age_yr;
  double plates12_bound;
}
subd_temp_custom_hscm2block_t;

/* struct for temperature options in subd_options_t */
typedef struct subd_temp_options
{
  subd_temp_type_t type;
  subd_temp_2plates_slab_t *slab_options;

  subd_temp_custom_t custom_type;
  subd_temp_custom_sinker_t *sinker_options;
  subd_temp_custom_lithblock_t *lithblock_options;
  subd_temp_custom_hscm2block_t *hscm2block_options;
}
subd_temp_options_t;

/*
 * weakzone options
 */
typedef enum
{
  SUBD_WEAK_NONE,
  SUBD_WEAK_SLAB
}
subd_weakzone_type_t;

/* struct for weak zone options in subd_options_t */
typedef struct subd_weak_options
{
  subd_weakzone_type_t                 type;
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
  subd_2plates_poly2_geo_coeff_t *weak_2plates_geo_coeff;
}
subd_weak_options_t;

/* enumerator for boundary conditions */
typedef enum
{
  SUBD_FREESURFACE,
  SUBD_MIXSURFACE,
  SUBD_SIDEWALL_YMIN_DIRALL,
  SUBD_SIDEWALLS_Y_DIRALL,
  SUBD_SIDEWALL_YMIN_DIRALL_FREEBASE,
}
subd_velbc_t;

/* enumerator for boundary conditions */
typedef enum
{
  SUBD_VELBC_DIR_ZERO,
  SUBD_VELBC_DIR_INOUTFLOW_SIN,
  SUBD_VELBC_DIR_INOUTFLOW_TANH_TWOLAYER,
  SUBD_VELBC_DIR_INOUTFLOW_TANH_THREELAYER,
  SUBD_VELBC_DIR_INOUTFLOW_DOUBLE_TANH_THREELAYER,
  SUBD_VELBC_DIR_MANUFACTURED_SINCOS_ISO,
  SUBD_VELBC_DIR_MANUFACTURED_POLY_ANISO,
}
subd_velbc_nonzero_dir_t;

typedef enum
{
  SUBD_VELBC_NEU_ZERO,
  SUBD_VELBC_NEU_VEC,
  SUBD_VELBC_NEU_SINE,
  SUBD_VELBC_NEU_READ,
}
subd_velbc_nonzero_neu_t;

/* struct for velocity boundary condition options in subd_options_t */
typedef struct subd_velbc_options
{
  subd_velbc_t  velbc;
  subd_velbc_nonzero_dir_t velbc_nonzero_dir;
  subd_velbc_nonzero_neu_t velbc_nonzero_neu;
  char                *txt_read_path_nonzero_neu;

  double              flow_scale;
  double              vel_dir_bc_middle;
  double              vel_dir_bc_upper;
  double              vel_dir_bc_lower;
}
subd_velbc_options_t;

/* enumerator for viscosity types */
typedef enum
{
  SUBD_VISC_ISOTROPY,
  SUBD_VISC_TRANSVERSELY_ISOTROPY
}
subd_visc_anisotropy_t;

typedef enum
{
 SUBD_VISC_RHEA,
 SUBD_VISC_CUSTOM
}
subd_visc_type_t;

typedef enum
{
  SUBD_VISC_CUSTOM_NONE,
  SUBD_VISC_CUSTOM_LAYERS,
  SUBD_VISC_CUSTOM_LAYERS_COUPLING
}
subd_visc_custom_t;


/* struct for viscosity options in subd_options_t */
typedef struct subd_visc_options
{
  subd_visc_anisotropy_t    anisotropy_type;
  subd_visc_type_t          type;
  subd_visc_custom_t        custom_type;
  double              z_lith;
  double              z_asthen;
  double              z_slab;
  double              slab_width;
  double              z_mantle;
  double              visc_lith;
  double              visc_asthen;
  double              visc_mantle;
}
subd_visc_options_t;

/* enumerator for domain shapes */
typedef enum
{
  SUBD_X_FUNCTION_IDENTITY,
  SUBD_X_FUNCTION_SINE,
  SUBD_X_FUNCTION_PROFILE
}
subd_x_func_t;

typedef struct subd_topo_profile
{
  int               nsurf;
  double            *tX;
  double            *tY;
  double            *tZ;
}
subd_topo_profile_t;

typedef struct subd_surf_options
{
  subd_x_func_t    x_func;
  subd_topo_profile_t *topo_profile;
}
subd_surf_options_t;


/* enumerator for tests */
typedef enum
{
  SUBD_TEST_STRESS_OP_NONE = 0,
  SUBD_TEST_STRESS_OP_SINCOS_SAME_OUTPUT,
  SUBD_TEST_STRESS_OP_SINCOS_ISO,
  SUBD_TEST_STRESS_OP_SINCOS_ANISO
}
subd_test_stress_op_t;

typedef enum
{
  SUBD_TEST_MANUFACTURED_NONE = 0,
  SUBD_TEST_MANUFACTURED_SINCOS1_ISO,
  SUBD_TEST_MANUFACTURED_SINCOS1_TIROT90,
  SUBD_TEST_MANUFACTURED_SINCOS1_TIROT45,
  SUBD_TEST_MANUFACTURED_SINCOS1_TIROT60,
  SUBD_TEST_MANUFACTURED_SINCOS1_TIROT60_VISCEXP60,
  SUBD_TEST_MANUFACTURED_POLY1_TIROT90,
  SUBD_TEST_MANUFACTURED_POLY1_TIROT90_VISCEXP
}
subd_test_manufactured_t;

typedef struct subd_test_options
{
  subd_test_stress_op_t    test_stress_op;
  subd_test_manufactured_t test_manufactured;
  subd_test_manufactured_t test_stress_comp;
}
subd_test_options_t;

typedef enum
{
  SUBD_ADJOINT_STENCIL_VISC_RHEA,
  SUBD_ADJOINT_STENCIL_VISC_CUSTOM
}
subd_adjoint_stencil_visc_type_t;

typedef enum
{
  SUBD_ADJOINT_STENCIL_VISC_CUSTOM_LAYERS_UM,
  SUBD_ADJOINT_STENCIL_VISC_CUSTOM_LAYERS_LM
}
subd_adjoint_stencil_visc_custom_t;

typedef struct subd_adjoint_stencil_options
{
  subd_adjoint_stencil_visc_type_t         type;
  subd_adjoint_stencil_visc_custom_t  custom_type;
  double                              value;
}
subd_adjoint_stencil_options_t;

typedef struct subd_adjoint_options
{
  subd_adjoint_stencil_options_t    *stencil_options;
}
subd_adjoint_options_t;

/* options of subd example */
typedef struct subd_options
{
  subd_rhs_type_t    rhs_type;
  rhea_domain_options_t   * domain_options;
  subd_para_options_t   *para_options;
  subd_temp_options_t   * temp_options;
  subd_visc_options_t   * visc_options;
  subd_weak_options_t   * weak_options;
  subd_surf_options_t   * surf_options;
  subd_velbc_options_t  * velbc_options;
  subd_test_options_t   * test_options;
  subd_adjoint_options_t * adjoint_options;
  void                  *data;
}
subd_options_t;

void
subduction_add_options (ymir_options_t * opt);

void
subduction_process_options (subd_options_t *subd_options,
                              rhea_domain_options_t *domain_options,
                              subd_para_options_t *para_options,
                              subd_temp_options_t *temp_options,
                              subd_visc_options_t *visc_options,
                              subd_weak_options_t *weak_options,
                              subd_surf_options_t *surf_options,
                              subd_velbc_options_t *velbc_options,
                              subd_test_options_t *test_options,
                              subd_adjoint_options_t *adjoint_options);

void
subd_options_destroy (subd_options_t *subd_options,
                       subd_temp_options_t  *temp_options,
                       subd_visc_options_t  *visc_options,
                       subd_weak_options_t  *weak_options,
                       subd_surf_options_t  *surf_options,
                       subd_velbc_options_t  *velbc_options,
                       subd_test_options_t  *test_options);

#endif /* SUBDUCTION_OPTIONS_H*/
