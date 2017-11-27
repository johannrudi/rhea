

#ifndef SLABS_OPTIONS_H
#define SLABS_OPTIONS_H

#include <rhea_domain.h>

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


void
slabweakzone_options (ymir_options_t * opt);

void
slabweakzone_options_process (slabs_options_t *slabs_options,
                              slabs_temp_options_t *temp_options,
                              slabs_visc_options_t *visc_options,
                              slabs_weak_options_t *weak_options,
                              slabs_surf_options_t *surf_options,
                              slabs_velbc_options_t *velbc_options,
                              slabs_test_options_t *test_options);

void
slabweakzone_options_domain_opt (slabs_options_t *slabs_options,
                                 slabs_domain_options_t *slabs_domain_options,
                                 rhea_domain_options_t  *domain_options);

#endif /* SLABS_OPTIONS_H */
