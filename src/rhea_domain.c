/*
 */

#include <rhea_domain.h>
#include <rhea_base.h>

double
rhea_domain_compute_radius (const double x, const double y, const double z,
                            rhea_domain_options_t *opt)
{
  switch (opt->shape) {
  case RHEA_DOMAIN_CUBE:
  case RHEA_DOMAIN_BRICK:
    /* transform interval of `z` to the interval of the shell slice's radius */
    {
      const double    z_max = opt->z_max;
      const double    radius_min = opt->radius_min;
      const double    radius_max = opt->radius_max;

      return z / z_max * (radius_max - radius_min) + radius_min;
    }
    break;

  case RHEA_DOMAIN_SHELL:
  case RHEA_DOMAIN_SHELL_CHUNK:
  case RHEA_DOMAIN_SHELL_SLICE:
    return sqrt (x * x + y * y + z * z);
    break;

  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }
}

/**
 * Computes the radius at the center of a hexahedral element in a shell domain
 * (or the corresponding value in case of a rectangular domain).
 */
static double
rhea_domain_compute_radius_at_elem_center (const double *x,
                                           const double *y,
                                           const double *z,
                                           const int *_sc_restrict Vmask,
                                           rhea_domain_options_t *opt)
{
  double              radii_sum = 0;
  int                 i;

  for (i = 0; i < 8; i++) { /* loop over all corners/vertices */
    const int           nodeid = Vmask[i];

    radii_sum += rhea_domain_compute_radius (x[nodeid], y[nodeid], z[nodeid],
                                             opt);
  }

  /* return mean value of radii at vertices */
  return radii_sum / 8.0;
}
