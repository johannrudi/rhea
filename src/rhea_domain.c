/*
 */

#include <rhea_domain.h>
#include <rhea_base.h>

/**
 * Computes the radius of a shell domain or the corresponding value for a
 * rectangular domain.
 */
double
slabs_compute_radius (const double x, const double y, const double z,
                      slabs_physics_options_t *physics_options)
{
  switch (physics_options->domain_shape) {
  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_BRICK:
    /* transform interval of `z` to the interval of the shell slice's radius */
    {
      const double    z_max = physics_options->domain_z_max;
      const double    radius_min = physics_options->domain_radius_min;
      const double    radius_max = physics_options->domain_radius_max;

      return z / z_max * (radius_max - radius_min) + radius_min;
    }
    break;

  case SL_DOMAIN_SHELL:
  case SL_DOMAIN_SHELL_CHUNK:
  case SL_DOMAIN_SHELL_SLICE:
    return sqrt (x * x + y * y + z * z);
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }
}

/**
 * Computes the radius at the center of an element in a shell domain or the
 * corresponding value for a rectangular domain.
 */
static double
slabs_compute_radius_at_elem_center (const double *x,
                                     const double *y,
                                     const double *z,
                                     const int *_sc_restrict Vmask,
                                     slabs_physics_options_t *physics_options)
{
  double              radii_sum = 0;
  int                 i;

  for (i = 0; i < 8; i++) { /* loop over all vertices */
    const int           nodeid = Vmask[i];

    radii_sum += slabs_compute_radius (x[nodeid], y[nodeid], z[nodeid],
                                       physics_options);
  }

  /* return mean value of radii at vertices */
  return radii_sum / 8.0;
}
