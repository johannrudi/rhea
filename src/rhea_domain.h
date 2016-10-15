/*
 */

#ifndef RHEA_DOMAIN_H
#define RHEA_DOMAIN_H

/* enumerator for domain shapes */
typedef enum
{
  RHEA_DOMAIN_CUBE,
  RHEA_DOMAIN_BRICK,
  RHEA_DOMAIN_SHELL,
  RHEA_DOMAIN_SHELL_CHUNK,
  RHEA_DOMAIN_SHELL_SLICE
}
rhea_domain_shape_t;

/* options & properties of domains */
typedef struct rhea_domain_options
{
  /* shape of the domain */
  rhea_domain_shape_t shape;

  /* the domain knows the location of the lower-upper mantle interface,
   * which causes discontinuous material properties */
  double              lm_um_interface_radius;
  double              lm_um_interface_smooth_transition_width;

  /* properties of the domain */
  double              x_min;
  double              x_max;
  double              y_min;
  double              y_max;
  double              z_min;
  double              z_max;
  double              lon_min;
  double              lon_max;
  double              radius_min;
  double              radius_max;

  double              volume;
  double              center[3];
  double              moment_of_inertia[3];
}
rhea_domain_options_t;

/**
 * Computes the radius of a shell domain or the corresponding value for a
 * rectangular domain.
 */
double              rhea_domain_compute_radius (const double x, const double y,
                                                const double z,
                                                rhea_domain_options_t *opt);

#endif /* RHEA_DOMAIN_H */
