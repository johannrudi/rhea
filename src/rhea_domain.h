/*
 */

#ifndef RHEA_DOMAIN_H
#define RHEA_DOMAIN_H

/* enumerator for domain shapes */
typedef enum
{
  SL_DOMAIN_CUBE,
  SL_DOMAIN_BRICK,
  SL_DOMAIN_SHELL,
  SL_DOMAIN_SHELL_CHUNK,
  SL_DOMAIN_SHELL_SLICE
}
slabs_domain_shape_t;

/* parameter list for mantle flow physics */
typedef struct slabs_physics_options
{
  double              domain_x_min;
  double              domain_x_max;
  double              domain_y_min;
  double              domain_y_max;
  double              domain_z_min;
  double              domain_z_max;
  double              domain_lon_min;
  double              domain_lon_max;
  double              domain_radius_min;
  double              domain_radius_max;

  double              viscosity_upper_mantle_radius;
  double              viscosity_lower_upper_transition_zone;

  double              domain_volume;
  double              domain_center[3];
  double              domain_moment_of_inertia[3];
}
slabs_physics_options_t;

#endif /* RHEA_DOMAIN_H */
