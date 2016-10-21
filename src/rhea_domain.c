/*
 */

#include <rhea_domain.h>
#include <rhea_base.h>
#include <rhea_discretization.h>

/* constants */
#define RHEA_DOMAIN_REFERENCE_SHELL_RADIUS (1.0)

/* default options */
#define RHEA_DOMAIN_DEFAULT_SHAPE_NAME "cube"
#define RHEA_DOMAIN_DEFAULT_BOX_X_EXTENSION (1)
#define RHEA_DOMAIN_DEFAULT_BOX_Y_EXTENSION (16)
#define RHEA_DOMAIN_DEFAULT_BOX_Z_EXTENSION (16)
#define RHEA_DOMAIN_DEFAULT_EARTH_RADIUS (6371.0e3)
#define RHEA_DOMAIN_DEFAULT_MANTLE_DEPTH (2866.95e3)
#define RHEA_DOMAIN_DEFAULT_LM_UM_INTERFACE_DEPTH (660.0e3)
#define RHEA_DOMAIN_DEFAULT_LM_UM_INTERFACE_SMOOTH_TRANSITION_WIDTH (0.0)

/* initialize options */
char               *rhea_domain_shape_name = RHEA_DOMAIN_DEFAULT_SHAPE_NAME;
int                 rhea_domain_box_x_extension =
  RHEA_DOMAIN_DEFAULT_BOX_X_EXTENSION;
int                 rhea_domain_box_y_extension =
  RHEA_DOMAIN_DEFAULT_BOX_Y_EXTENSION;
int                 rhea_domain_box_z_extension =
  RHEA_DOMAIN_DEFAULT_BOX_Z_EXTENSION;
double              rhea_domain_earth_radius = RHEA_DOMAIN_DEFAULT_EARTH_RADIUS;
double              rhea_domain_mantle_depth = RHEA_DOMAIN_DEFAULT_MANTLE_DEPTH;
double              rhea_domain_lm_um_interface_depth =
  RHEA_DOMAIN_DEFAULT_LM_UM_INTERFACE_DEPTH;
double              rhea_domain_lm_um_interface_smooth_transition_width =
  RHEA_DOMAIN_DEFAULT_LM_UM_INTERFACE_SMOOTH_TRANSITION_WIDTH;

void
rhea_domain_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Domain";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_S, "shape", '\0',
    &rhea_domain_shape_name, RHEA_DOMAIN_DEFAULT_SHAPE_NAME,
    "Shape of domain: cube, box, shell, cube_spherical, box_spherical",

  YMIR_OPTIONS_I, "box-x-extension", '\0',
    &rhea_domain_box_x_extension, RHEA_DOMAIN_DEFAULT_BOX_X_EXTENSION,
    "For 'box' domain: Integer that determines extension in x-direction",
  YMIR_OPTIONS_I, "box-y-extension", '\0',
    &rhea_domain_box_y_extension, RHEA_DOMAIN_DEFAULT_BOX_Y_EXTENSION,
    "For 'box' domain: Integer that determines extension in y-direction",
  YMIR_OPTIONS_I, "box-z-extension", '\0',
    &rhea_domain_box_z_extension, RHEA_DOMAIN_DEFAULT_BOX_Z_EXTENSION,
    "For 'box' domain: Integer that determines extension in z-direction",

  YMIR_OPTIONS_D, "earth-radius", '\0',
    &rhea_domain_earth_radius, RHEA_DOMAIN_DEFAULT_EARTH_RADIUS,
    "Mean radius of the earth [m]",
  YMIR_OPTIONS_D, "mantle-depth", '\0',
    &rhea_domain_mantle_depth, RHEA_DOMAIN_DEFAULT_MANTLE_DEPTH,
    "Mean depth of the (whole) mantle [m]",
  YMIR_OPTIONS_D, "lower-upper-mantle-interface-depth", '\0',
    &rhea_domain_lm_um_interface_depth,
    RHEA_DOMAIN_DEFAULT_LM_UM_INTERFACE_DEPTH,
    "Depth of interface between lower and upper mantle [m]",
  YMIR_OPTIONS_D, "lower-upper-mantle-interface-smooth-transition-width", '\0',
    &rhea_domain_lm_um_interface_smooth_transition_width,
    RHEA_DOMAIN_DEFAULT_LM_UM_INTERFACE_SMOOTH_TRANSITION_WIDTH,
    "Width of smooth transition zone between lower and upper mantle [m]",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);
}

/**
 * Computes and stores boundaries of a domain.
 */
static void
rhea_domain_compute_bounds (rhea_domain_options_t *opt)
{
  const double        earth_radius = rhea_domain_earth_radius;
  const double        mantle_depth = rhea_domain_mantle_depth;
  const double        box_x_ext = (double) opt->box_x_extension;
  const double        box_y_ext = (double) opt->box_y_extension;
  const double        box_z_ext = (double) opt->box_z_extension;

  /* set radius bounds for all domain shapes */
  opt->radius_max = RHEA_DOMAIN_REFERENCE_SHELL_RADIUS;
  opt->radius_min = opt->radius_max *
                    SC_MAX (0.0, 1.0 - mantle_depth/earth_radius);

  /* initialize other values */
  opt->x_min = NAN;
  opt->x_max = NAN;
  opt->y_min = NAN;
  opt->y_max = NAN;
  opt->z_min = NAN;
  opt->z_max = NAN;
  opt->lon_min = NAN;
  opt->lon_max = NAN;

  /* set domain-specific bounds */
  switch (opt->shape) {
  case RHEA_DOMAIN_CUBE:
    opt->x_min = 0.0;
    opt->x_max = 1.0;
    opt->y_min = 0.0;
    opt->y_max = 1.0;
    opt->z_min = 0.0;
    opt->z_max = 1.0;
    opt->lon_min = -M_PI/8.0;
    opt->lon_max = +M_PI/8.0;
    break;
  case RHEA_DOMAIN_BOX:
    opt->x_min = 0.0;
    opt->x_max = box_x_ext / box_z_ext;
    opt->y_min = 0.0;
    opt->y_max = box_y_ext / box_z_ext;
    opt->z_min = 0.0;
    opt->z_max = 1.0;
    opt->lon_min = -M_PI/8.0 * box_y_ext / box_z_ext;
    opt->lon_max = +M_PI/8.0 * box_y_ext / box_z_ext;
    break;
  case RHEA_DOMAIN_SHELL:
    break;
  case RHEA_DOMAIN_CUBE_SPHERICAL:
    opt->lon_min = -M_PI/8.0;
    opt->lon_max = +M_PI/8.0;
    break;
  case RHEA_DOMAIN_BOX_SPHERICAL:
    opt->lon_min = -M_PI/8.0 * box_y_ext / box_z_ext;
    opt->lon_max = +M_PI/8.0 * box_y_ext / box_z_ext;
    break;
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }
}

/**
 * Computes and stores the volume of a domain.
 */
static void
rhea_domain_compute_volume (rhea_domain_options_t *opt)
{
  switch (opt->shape) {
  case RHEA_DOMAIN_CUBE:
  case RHEA_DOMAIN_BOX:
    RHEA_ASSERT (isfinite (opt->x_min));
    RHEA_ASSERT (isfinite (opt->x_max));
    RHEA_ASSERT (isfinite (opt->y_min));
    RHEA_ASSERT (isfinite (opt->y_max));
    RHEA_ASSERT (isfinite (opt->z_min));
    RHEA_ASSERT (isfinite (opt->z_max));
    {
      const double        dx = opt->x_max - opt->x_min;
      const double        dy = opt->y_max - opt->y_min;
      const double        dz = opt->z_max - opt->z_min;

      opt->volume = dx * dy * dz;
    }
    break;
  case RHEA_DOMAIN_SHELL: /* (4/3 * pi * (radius_top^3 - radius_bottom^3)) */
    RHEA_ASSERT (isfinite (opt->radius_min));
    RHEA_ASSERT (isfinite (opt->radius_max));
    {
      const double        rmin = opt->radius_min;
      const double        rmax = opt->radius_max;

      opt->volume = 4.0 / 3.0 * M_PI * (rmax*rmax*rmax - rmin*rmin*rmin);
    }
    break;
  case RHEA_DOMAIN_CUBE_SPHERICAL: /* (volume shell / 24) */
    RHEA_ASSERT (isfinite (opt->radius_min));
    RHEA_ASSERT (isfinite (opt->radius_max));
    {
      const double        rmin = opt->radius_min;
      const double        rmax = opt->radius_max;

      opt->volume = 1.0 / 18.0 * M_PI * (rmax*rmax*rmax - rmin*rmin*rmin);
    }
    break;
  case RHEA_DOMAIN_BOX_SPHERICAL: /* (volume shell / 24 * (dx*dy)/dz^2 */
    RHEA_ASSERT (isfinite (opt->radius_min));
    RHEA_ASSERT (isfinite (opt->radius_max));
    RHEA_ASSERT (isfinite (opt->box_x_extension));
    RHEA_ASSERT (isfinite (opt->box_y_extension));
    RHEA_ASSERT (isfinite (opt->box_z_extension));
    {
      const double        rmin = opt->radius_min;
      const double        rmax = opt->radius_max;
      const double        box_x_ext = (double) opt->box_x_extension;
      const double        box_y_ext = (double) opt->box_y_extension;
      const double        box_z_ext = (double) opt->box_z_extension;

      opt->volume = 1.0 / 18.0 * M_PI * (rmax*rmax*rmax - rmin*rmin*rmin) *
                    (box_x_ext * box_y_ext) / (box_z_ext * box_z_ext);
    }
    break;
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }
  RHEA_ASSERT (isfinite (opt->volume));
}

/**
 * Computes and stores the center of mass of a domain.
 */
static void
rhea_domain_compute_center_of_mass (rhea_domain_options_t *opt)
{
  switch (opt->shape) {
  case RHEA_DOMAIN_CUBE:
  case RHEA_DOMAIN_BOX:
    RHEA_ASSERT (isfinite (opt->x_min));
    RHEA_ASSERT (isfinite (opt->x_max));
    RHEA_ASSERT (isfinite (opt->y_min));
    RHEA_ASSERT (isfinite (opt->y_max));
    RHEA_ASSERT (isfinite (opt->z_min));
    RHEA_ASSERT (isfinite (opt->z_max));
    opt->center[0] = 0.5 * (opt->x_min + opt->x_max);
    opt->center[1] = 0.5 * (opt->y_min + opt->y_max);
    opt->center[2] = 0.5 * (opt->z_min + opt->z_max);
    break;
  case RHEA_DOMAIN_SHELL:
    opt->center[0] = 0.0;
    opt->center[1] = 0.0;
    opt->center[2] = 0.0;
    break;
  case RHEA_DOMAIN_CUBE_SPHERICAL:
  case RHEA_DOMAIN_BOX_SPHERICAL:
    RHEA_ASSERT (isfinite (opt->radius_min));
    RHEA_ASSERT (isfinite (opt->radius_max));
    opt->center[0] = 0.0;
    opt->center[1] = 0.0;
    opt->center[2] = 0.5 * (opt->radius_min + opt->radius_max);
    //TODO taking the middle of the radii is just a crude approximation
    break;
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }

  RHEA_ASSERT (isfinite (opt->center[0]));
  RHEA_ASSERT (isfinite (opt->center[1]));
  RHEA_ASSERT (isfinite (opt->center[2]));
}

/**
 * Computes and stores the moment of inertia of a domain.
 */
static void
rhea_domain_compute_moment_of_inertia (rhea_domain_options_t *opt)
{
  switch (opt->shape) {
  case RHEA_DOMAIN_CUBE: /* (mass * side_length^2 / 6) */
    {
      const double        moment_of_inertia = 1.0 / 6.0;

      opt->moment_of_inertia[0] = moment_of_inertia;
      opt->moment_of_inertia[1] = moment_of_inertia;
      opt->moment_of_inertia[2] = moment_of_inertia;
    }
    break;
  case RHEA_DOMAIN_BOX: /* x: (mass/12 * (dy^2 + dz^2))
                         * y: (mass/12 * (dx^2 + dz^2))
                         * z: (mass/12 * (dx^2 + dy^2)) */
    RHEA_ASSERT (isfinite (opt->volume));
    RHEA_ASSERT (isfinite (opt->x_min));
    RHEA_ASSERT (isfinite (opt->x_max));
    RHEA_ASSERT (isfinite (opt->y_min));
    RHEA_ASSERT (isfinite (opt->y_max));
    RHEA_ASSERT (isfinite (opt->z_min));
    RHEA_ASSERT (isfinite (opt->z_max));
    {
      const double        vol = opt->volume;
      const double        dx = opt->x_max - opt->x_min;
      const double        dy = opt->y_max - opt->y_min;
      const double        dz = opt->z_max - opt->z_min;

      opt->moment_of_inertia[0] = vol/12.0 * (dy*dy + dz*dz);
      opt->moment_of_inertia[1] = vol/12.0 * (dx*dx + dz*dz);
      opt->moment_of_inertia[2] = vol/12.0 * (dx*dx + dy*dy);
    }
    break;
  case RHEA_DOMAIN_SHELL: /* (2/5 * mass * (radius_top^5 - radius_bottom^5)
                           *             / (radius_top^3 - radius_bottom^3)) */
    RHEA_ASSERT (isfinite (opt->radius_min));
    RHEA_ASSERT (isfinite (opt->radius_max));
    {
      const double        rmin = opt->radius_min;
      const double        rmax = opt->radius_max;
      double              moment_of_inertia;

      moment_of_inertia = (8.0 / 15.0) * M_PI * (pow(rmax, 5) - pow(rmin, 5));
      opt->moment_of_inertia[0] = moment_of_inertia;
      opt->moment_of_inertia[1] = moment_of_inertia;
      opt->moment_of_inertia[2] = moment_of_inertia;
    }
    break;
  case RHEA_DOMAIN_CUBE_SPHERICAL:
  case RHEA_DOMAIN_BOX_SPHERICAL:
    //TODO
    opt->moment_of_inertia[0] = 0.0;
    opt->moment_of_inertia[1] = 0.0;
    opt->moment_of_inertia[2] = 0.0;
    break;

  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }
  RHEA_ASSERT (isfinite (opt->moment_of_inertia[0]));
  RHEA_ASSERT (isfinite (opt->moment_of_inertia[1]));
  RHEA_ASSERT (isfinite (opt->moment_of_inertia[2]));
}

/**
 *
 */
static double
rhea_domain_align_radius_to_mesh (double radius, const int level,
                                  rhea_domain_options_t *opt)
{
  const double        rmin = opt->radius_min;
  const double        rmax = opt->radius_max;
  const double        rmax_rmin = rmax / rmin;
  const double        rmin2_rmax = rmin*rmin / rmax;
  int                 n_domain_subdiv;
  double              h;

  /* check input */
  RHEA_ASSERT (isfinite (opt->radius_min));
  RHEA_ASSERT (isfinite (opt->radius_max));
  RHEA_ASSERT (0.0 <= opt->radius_min);
  RHEA_ASSERT (opt->radius_min < opt->radius_max);

  /* exit if nothing to do */
  if ( !(rmin < radius && radius < rmax) ) {
    return -1.0;
  }

  /* set number of subdivisions given by the extension in z-direction */
  switch (opt->shape) {
  case RHEA_DOMAIN_CUBE:
  case RHEA_DOMAIN_SHELL:
  case RHEA_DOMAIN_CUBE_SPHERICAL:
    n_domain_subdiv = 1;
    break;
  case RHEA_DOMAIN_BOX:
  case RHEA_DOMAIN_BOX_SPHERICAL:
    n_domain_subdiv = opt->box_z_extension;
    break;
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }
  RHEA_ASSERT (1 <= n_domain_subdiv);

  /* calculate aligned radius */
  switch (opt->shape) {
  case RHEA_DOMAIN_CUBE:
  case RHEA_DOMAIN_BOX:
    /* calculate element size: h = D / (n * 2^l) */
    h = (rmax - rmin) / ((double) (n_domain_subdiv << level));

    /* calculate closest aligned radius below `radius` */
    radius = rmin + h * floor ((radius - rmin) / h);
    break;

  case RHEA_DOMAIN_SHELL:
  case RHEA_DOMAIN_CUBE_SPHERICAL:
  case RHEA_DOMAIN_BOX_SPHERICAL:
    /* map radius from physical space to reference space */
    radius = log (radius / rmin2_rmax) / log (rmax_rmin);

    /* calculate element size: h = 1 / (n * 2^l) */
    h = 1.0 / ((double) (n_domain_subdiv << level));

    /* calculate closest aligned radius below `radius` */
    radius = 1.0 + h * floor ((radius - 1.0) / h);

    /* map radius from reference space to physical space */
    radius = rmin2_rmax * pow (rmax_rmin, radius);
    break;

  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* return aligned radius */
  return radius;
}

void
rhea_domain_process_options (rhea_domain_options_t *opt)
{
  const char         *this_fn_name = "rhea_domain_process_options";
  double              radius;

  /* set shape of domain */
  if (strcmp (rhea_domain_shape_name, "cube") == 0) {
    opt->shape = RHEA_DOMAIN_CUBE;
  }
  else if (strcmp (rhea_domain_shape_name, "box") == 0) {
    opt->shape = RHEA_DOMAIN_BOX;
  }
  else if (strcmp (rhea_domain_shape_name, "shell") == 0) {
    opt->shape = RHEA_DOMAIN_SHELL;
  }
  else if (strcmp (rhea_domain_shape_name, "cube_spherical") == 0) {
    opt->shape = RHEA_DOMAIN_CUBE_SPHERICAL;
  }
  else if (strcmp (rhea_domain_shape_name, "box_spherical") == 0) {
    opt->shape = RHEA_DOMAIN_BOX_SPHERICAL;
  }
  else { /* unknown shape name */
    RHEA_ABORT ("Unknown domain shape name");
  }

  /* set extension of box domain */
  switch (opt->shape) {
  case RHEA_DOMAIN_BOX:
  case RHEA_DOMAIN_BOX_SPHERICAL:
    /* for box-type domains, set from input  */
    opt->box_x_extension = rhea_domain_box_x_extension;
    opt->box_y_extension = rhea_domain_box_y_extension;
    opt->box_z_extension = rhea_domain_box_z_extension;
    break;
  case RHEA_DOMAIN_CUBE:
  case RHEA_DOMAIN_CUBE_SPHERICAL:
  case RHEA_DOMAIN_SHELL:
    /* for all other domains, set to a default value */
    opt->box_x_extension = 1;
    opt->box_y_extension = 1;
    opt->box_z_extension = 1;
    break;
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* compute domain properties */
  rhea_domain_compute_bounds (opt);
  rhea_domain_compute_volume (opt);
  rhea_domain_compute_center_of_mass (opt);
  rhea_domain_compute_moment_of_inertia (opt);

  /* align lower-upper mantle interface with mesh elements */
  RHEA_ASSERT (isfinite (opt->radius_max) && 0.0 < opt->radius_max);
  radius = opt->radius_max *
           (1.0 - rhea_domain_lm_um_interface_depth / rhea_domain_earth_radius);
  opt->lm_um_interface_radius = rhea_domain_align_radius_to_mesh (
      radius, rhea_discretization_level_min, opt);

  /* set lower-upper mantle interface */
  if (0.0 < rhea_domain_lm_um_interface_smooth_transition_width) {
    opt->lm_um_interface_smooth_transition_width =
      rhea_domain_lm_um_interface_smooth_transition_width /
      rhea_domain_earth_radius;
  }
  else {
    opt->lm_um_interface_smooth_transition_width = 0.0;
  }

  /* print derived domain properties */
  RHEA_GLOBAL_INFO ("===================================================\n");
  RHEA_GLOBAL_INFOF ("%s\n", this_fn_name);
  RHEA_GLOBAL_INFO ("---------------------------------------------------\n");
  RHEA_GLOBAL_INFOF ("  x min:      %g\n", opt->x_min);
  RHEA_GLOBAL_INFOF ("  x max:      %g\n", opt->x_max);
  RHEA_GLOBAL_INFOF ("  y min:      %g\n", opt->y_min);
  RHEA_GLOBAL_INFOF ("  y max:      %g\n", opt->y_max);
  RHEA_GLOBAL_INFOF ("  z min:      %g\n", opt->z_min);
  RHEA_GLOBAL_INFOF ("  z max:      %g\n", opt->z_max);
  RHEA_GLOBAL_INFOF ("  lon min:    %g\n", opt->lon_min);
  RHEA_GLOBAL_INFOF ("  lon max:    %g\n", opt->lon_max);
  RHEA_GLOBAL_INFOF ("  radius min: %g\n", opt->radius_min);
  RHEA_GLOBAL_INFOF ("  radius max: %g\n", opt->radius_max);
  RHEA_GLOBAL_INFOF ("  LM-UM interface radius: %g\n",
                     opt->lm_um_interface_radius);
  RHEA_GLOBAL_INFOF ("  LM-UM interface depth:  %g (%g km)\n",
                     opt->radius_max - opt->lm_um_interface_radius,
                     (opt->radius_max - opt->lm_um_interface_radius) *
                     rhea_domain_earth_radius / 1.0e3);
  RHEA_GLOBAL_INFOF ("  volume:            %g\n", opt->volume);
  RHEA_GLOBAL_INFOF ("  center of mass:    %g, %g, %g\n",
                     opt->center[0], opt->center[1], opt->center[2]);
  RHEA_GLOBAL_INFOF ("  moment of inertia: %g, %g, %g\n",
                     opt->moment_of_inertia[0],  opt->moment_of_inertia[1],
                     opt->moment_of_inertia[2]);
  RHEA_GLOBAL_INFO ("===================================================\n");
}

double
rhea_domain_compute_radius (const double x, const double y, const double z,
                            rhea_domain_options_t *opt)
{
  switch (opt->shape) {
  case RHEA_DOMAIN_CUBE:
  case RHEA_DOMAIN_BOX:
    /* transform interval of `z` to the interval of the shell slice's radius */
    {
      const double    z_max = opt->z_max;
      const double    radius_min = opt->radius_min;
      const double    radius_max = opt->radius_max;

      return z / z_max * (radius_max - radius_min) + radius_min;
    }
    break;

  case RHEA_DOMAIN_SHELL:
  case RHEA_DOMAIN_CUBE_SPHERICAL:
  case RHEA_DOMAIN_BOX_SPHERICAL:
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
rhea_domain_compute_radius_at_elem_center (const double *_sc_restrict x,
                                           const double *_sc_restrict y,
                                           const double *_sc_restrict z,
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

int
rhea_domain_elem_is_in_upper_mantle (const double *x, const double *y,
                                     const double *z, const int *Vmask,
                                     rhea_domain_options_t *opt)
{
  const double        lm_um_interface_radius = opt->lm_um_interface_radius;
  const double        r_center =
    rhea_domain_compute_radius_at_elem_center (x, y, z, Vmask, opt);

  /* return if radius of the element's center is in upper mantle */
  return (lm_um_interface_radius <= r_center);
}
