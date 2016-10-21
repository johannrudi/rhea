/*
 */

#ifndef RHEA_DOMAIN_H
#define RHEA_DOMAIN_H

#include <ymir_options.h>

/* enumerator for domain shapes */
typedef enum
{
  RHEA_DOMAIN_CUBE,
  RHEA_DOMAIN_BOX,
  RHEA_DOMAIN_SHELL,
  RHEA_DOMAIN_CUBE_SPHERICAL,
  RHEA_DOMAIN_BOX_SPHERICAL
}
rhea_domain_shape_t;

/* options & properties of a computational domain */
typedef struct rhea_domain_options
{
  /* shape of the domain */
  rhea_domain_shape_t shape;

  /* extension in each Cartesian direction for `box` domain */
  double              box_x_extension;
  double              box_y_extension;
  double              box_z_extension;

  /* the domain knows the location of the lower-upper mantle interface,
   * which causes discontinuous material properties */
  double              lm_um_interface_radius;
  double              lm_um_interface_smooth_transition_width;

  /* properties of the domain (reference coordinate system) */
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
 * Defines options and adds them as sub-options.
 */
void                rhea_domain_add_options (ymir_options_t * opt_sup);

/**
 * Processes options and stores them.
 */
void                rhea_domain_process_options (rhea_domain_options_t *opt);

/**
 * Computes the radius of a shell domain or the corresponding value for a
 * rectangular domain.
 */
double              rhea_domain_compute_radius (const double x, const double y,
                                                const double z,
                                                rhea_domain_options_t *opt);

/**
 * Returns whether the element's center is located in the upper mantle.
 */
int                 rhea_domain_elem_is_in_upper_mantle (
                                                  const double *x,
                                                  const double *y,
                                                  const double *z,
                                                  const int *Vmask,
                                                  rhea_domain_options_t *opt);

#endif /* RHEA_DOMAIN_H */
