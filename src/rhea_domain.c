/*
 */

#include <rhea_domain.h>
#include <rhea_base.h>
#include <rhea_discretization.h>

/* constant: reference value for shell radius */
#define RHEA_DOMAIN_REFERENCE_SHELL_RADIUS (1.0)

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_DOMAIN_DEFAULT_SHAPE_NAME "cube"
#define RHEA_DOMAIN_DEFAULT_BOX_LENGTH_X (NAN)
#define RHEA_DOMAIN_DEFAULT_BOX_LENGTH_Y (NAN)
#define RHEA_DOMAIN_DEFAULT_BOX_LENGTH_Z (NAN)
#define RHEA_DOMAIN_DEFAULT_BOX_SUBDIVISION_X (16)
#define RHEA_DOMAIN_DEFAULT_BOX_SUBDIVISION_Y (1)
#define RHEA_DOMAIN_DEFAULT_BOX_SUBDIVISION_Z (16)
#define RHEA_DOMAIN_DEFAULT_BOX_SPHERICAL_DISTORTION_CORR (1)
#define RHEA_DOMAIN_DEFAULT_EARTH_RADIUS_M (6371.0e3)
#define RHEA_DOMAIN_DEFAULT_MANTLE_DEPTH_M (2866.95e3)
#define RHEA_DOMAIN_DEFAULT_LM_UM_INTERFACE_DEPTH_M (660.0e3)
#define RHEA_DOMAIN_DEFAULT_LM_UM_INTERFACE_SMOOTHING_WIDTH_M (0.0)
#define RHEA_DOMAIN_DEFAULT_VELOCITY_BC_TYPE \
  RHEA_DOMAIN_VELOCITY_BC_DIRICHLET_ALL

/* initialize options */
char               *rhea_domain_shape_name = RHEA_DOMAIN_DEFAULT_SHAPE_NAME;
double              rhea_domain_box_length_x =
  RHEA_DOMAIN_DEFAULT_BOX_LENGTH_X;
double              rhea_domain_box_length_y =
  RHEA_DOMAIN_DEFAULT_BOX_LENGTH_Y;
double              rhea_domain_box_length_z =
  RHEA_DOMAIN_DEFAULT_BOX_LENGTH_Z;
int                 rhea_domain_box_subdivision_x =
  RHEA_DOMAIN_DEFAULT_BOX_SUBDIVISION_X;
int                 rhea_domain_box_subdivision_y =
  RHEA_DOMAIN_DEFAULT_BOX_SUBDIVISION_Y;
int                 rhea_domain_box_subdivision_z =
  RHEA_DOMAIN_DEFAULT_BOX_SUBDIVISION_Z;
int                 rhea_domain_box_spherical_distortion_corr =
  RHEA_DOMAIN_DEFAULT_BOX_SPHERICAL_DISTORTION_CORR;
double              rhea_domain_earth_radius_m =
  RHEA_DOMAIN_DEFAULT_EARTH_RADIUS_M;
double              rhea_domain_mantle_depth_m =
  RHEA_DOMAIN_DEFAULT_MANTLE_DEPTH_M;
double              rhea_domain_lm_um_interface_depth_m =
  RHEA_DOMAIN_DEFAULT_LM_UM_INTERFACE_DEPTH_M;
double              rhea_domain_lm_um_interface_smoothing_width_m =
  RHEA_DOMAIN_DEFAULT_LM_UM_INTERFACE_SMOOTHING_WIDTH_M;
int                 rhea_domain_velocity_bc_type =
  RHEA_DOMAIN_DEFAULT_VELOCITY_BC_TYPE;

void
rhea_domain_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Domain";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_S, "shape", '\0',
    &(rhea_domain_shape_name), RHEA_DOMAIN_DEFAULT_SHAPE_NAME,
    "Shape of domain: cube, box, shell, cube_spherical, box_spherical",

  YMIR_OPTIONS_D, "box-length-x", '\0',
    &(rhea_domain_box_length_x), RHEA_DOMAIN_DEFAULT_BOX_LENGTH_X,
    "For 'box' domain: Length along x-axis; if not given, value is set "
    "automatically to avoid mesh anisotropies",
  YMIR_OPTIONS_D, "box-length-y", '\0',
    &(rhea_domain_box_length_y), RHEA_DOMAIN_DEFAULT_BOX_LENGTH_Y,
    "For 'box' domain: Length along y-axis; if not given, value is set "
    "automatically to avoid mesh anisotropies",
  YMIR_OPTIONS_D, "box-length-z", '\0',
    &(rhea_domain_box_length_z), RHEA_DOMAIN_DEFAULT_BOX_LENGTH_Z,
    "For 'box' domain: Length along z-axis; if not given, value is set "
    "automatically to avoid mesh anisotropies",

  YMIR_OPTIONS_I, "box-subdivision-x", '\0',
    &(rhea_domain_box_subdivision_x), RHEA_DOMAIN_DEFAULT_BOX_SUBDIVISION_X,
    "For 'box' domain: Integer determining subdivision into cubes along x-axis",
  YMIR_OPTIONS_I, "box-subdivision-y", '\0',
    &(rhea_domain_box_subdivision_y), RHEA_DOMAIN_DEFAULT_BOX_SUBDIVISION_Y,
    "For 'box' domain: Integer determining subdivision into cubes along x-axis",
  YMIR_OPTIONS_I, "box-subdivision-z", '\0',
    &(rhea_domain_box_subdivision_z), RHEA_DOMAIN_DEFAULT_BOX_SUBDIVISION_Z,
    "For 'box' domain: Integer determining subdivision into cubes along x-axis",

  YMIR_OPTIONS_B, "box-spherical-distortion-correction", '\0',
    &(rhea_domain_box_spherical_distortion_corr),
    RHEA_DOMAIN_DEFAULT_BOX_SPHERICAL_DISTORTION_CORR,
    "Correct a pillow-type distortion of domain shape `box_spherical`",

  YMIR_OPTIONS_D, "earth-radius", '\0',
    &(rhea_domain_earth_radius_m), RHEA_DOMAIN_DEFAULT_EARTH_RADIUS_M,
    "Mean radius of the earth [m]",
  YMIR_OPTIONS_D, "mantle-depth", '\0',
    &(rhea_domain_mantle_depth_m), RHEA_DOMAIN_DEFAULT_MANTLE_DEPTH_M,
    "Mean depth of the (whole) mantle [m]",
  YMIR_OPTIONS_D, "lower-upper-mantle-interface-depth", '\0',
    &(rhea_domain_lm_um_interface_depth_m),
    RHEA_DOMAIN_DEFAULT_LM_UM_INTERFACE_DEPTH_M,
    "Depth of interface between lower and upper mantle [m] "
    "(0: lower mantle everywhere, -1: upper mantle everywhere)",
  YMIR_OPTIONS_D, "lower-upper-mantle-interface-smooth-transition-width", '\0',
    &(rhea_domain_lm_um_interface_smoothing_width_m),
    RHEA_DOMAIN_DEFAULT_LM_UM_INTERFACE_SMOOTHING_WIDTH_M,
    "Width of smooth transition zone between lower and upper mantle [m]",

  YMIR_OPTIONS_I, "velocity-boundary-condition", '\0',
    &(rhea_domain_velocity_bc_type), RHEA_DOMAIN_DEFAULT_VELOCITY_BC_TYPE,
    "Velocity boundary condition: "
    "-1: Provided by user; "
    "0: Dirichlet all; "
    "1: Dirichlet in normal direction; "
    "2: Dirichlet norm (if shell: two inner points fixed for rotation inv); "
    "3: Dirichlet all on inner sphere & Dirichlet norm on outer sphere; "
    "4: Dirichlet all on side faces & Neumann on top, bottom faces; "
    "5: Dirichlet norm on side faces & Dirichlet all on top, bottom faces; "
    "6: Neumann all",

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
  const double        radius_min_m = opt->radius_min_m;
  const double        radius_max_m = opt->radius_max_m;

  /* set radius bounds for all domain shapes */
  opt->radius_max = RHEA_DOMAIN_REFERENCE_SHELL_RADIUS;
  opt->radius_min = opt->radius_max * radius_min_m/radius_max_m;

  /* set domain-specific bounds */
  switch (opt->shape) {
  case RHEA_DOMAIN_CUBE:
    opt->x_min = 0.0;
    opt->y_min = 0.0;
    opt->z_min = 0.0;
    opt->x_max = 1.0;
    opt->y_max = 1.0;
    opt->z_max = 1.0;
    opt->lon_min = -M_PI/8.0;
    opt->lon_max = +M_PI/8.0;
    break;
  case RHEA_DOMAIN_BOX:
    RHEA_ASSERT (isfinite (opt->box_length_x) && 0.0 < opt->box_length_x);
    RHEA_ASSERT (isfinite (opt->box_length_y) && 0.0 < opt->box_length_y);
    RHEA_ASSERT (isfinite (opt->box_length_z) && 0.0 < opt->box_length_z);
    opt->x_min = 0.0;
    opt->y_min = 0.0;
    opt->z_min = 0.0;
    opt->x_max = opt->box_length_x;
    opt->y_max = opt->box_length_y;
    opt->z_max = opt->box_length_z;
    opt->lon_min = -M_PI/8.0 * opt->box_length_y;
    opt->lon_max = +M_PI/8.0 * opt->box_length_y;
    break;
  case RHEA_DOMAIN_SHELL:
    opt->x_min = -1.0;
    opt->y_min = -1.0;
    opt->z_min = -1.0;
    opt->x_max = +1.0;
    opt->y_max = +1.0;
    opt->z_max = +1.0;
    opt->lon_min = NAN;
    opt->lon_max = NAN;
    break;
  case RHEA_DOMAIN_CUBE_SPHERICAL:
    opt->x_min = -1.0;
    opt->y_min = -1.0;
    opt->z_min = -1.0;
    opt->x_max = +1.0;
    opt->y_max = +1.0;
    opt->z_max = +1.0;
    opt->lon_min = -M_PI/8.0;
    opt->lon_max = +M_PI/8.0;
    break;
  case RHEA_DOMAIN_BOX_SPHERICAL:
    opt->x_min = -1.0;
    opt->y_min = -1.0;
    opt->z_min = -1.0;
    opt->x_max = +1.0;
    opt->y_max = +1.0;
    opt->z_max = +1.0;
    RHEA_ASSERT (0 < opt->box_subdivision_y);
    RHEA_ASSERT (0 < opt->box_subdivision_z);
    {
      const double        suby = (double) opt->box_subdivision_y;
      const double        subz = (double) opt->box_subdivision_z;

      opt->lon_min = -M_PI/8.0 * suby / subz;
      opt->lon_max = +M_PI/8.0 * suby / subz;
    }
    break;
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }
}

/**
 * Computes and stores the representative depth of the domain.
 */
static void
rhea_domain_compute_representative_depth (rhea_domain_options_t *opt)
{
  switch (opt->shape) {
  case RHEA_DOMAIN_CUBE:
  case RHEA_DOMAIN_BOX:
    RHEA_ASSERT (isfinite (opt->z_min));
    RHEA_ASSERT (isfinite (opt->z_max));
    opt->depth = opt->z_max - opt->z_min;
    break;
  case RHEA_DOMAIN_SHELL:
  case RHEA_DOMAIN_CUBE_SPHERICAL:
  case RHEA_DOMAIN_BOX_SPHERICAL:
    RHEA_ASSERT (isfinite (opt->radius_min));
    RHEA_ASSERT (isfinite (opt->radius_max));
    opt->depth = opt->radius_max - opt->radius_min;
    break;
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }
  RHEA_ASSERT (isfinite (opt->depth));
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
    RHEA_ASSERT (isfinite (opt->box_subdivision_x));
    RHEA_ASSERT (isfinite (opt->box_subdivision_y));
    RHEA_ASSERT (isfinite (opt->box_subdivision_z));
    {
      const double        rmin = opt->radius_min;
      const double        rmax = opt->radius_max;
      const double        subx = (double) opt->box_subdivision_x;
      const double        suby = (double) opt->box_subdivision_y;
      const double        subz = (double) opt->box_subdivision_z;

      opt->volume = 1.0 / 18.0 * M_PI * (rmax*rmax*rmax - rmin*rmin*rmin) *
                    (subx*suby) / (subz*subz);
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
rhea_domain_align_radius_with_mesh (double radius, const int level,
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
  if (radius <= rmin) {
    return rmin;
  }
  else if (rmax <= radius) {
    return rmax;
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
    n_domain_subdiv = opt->box_subdivision_z;
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
  double              scal_len_to_km, radius;

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

  /* set length and subdivision of the domain */
  switch (opt->shape) {
  case RHEA_DOMAIN_CUBE:
    /* set default values */
    opt->box_subdivision_x = 1;
    opt->box_subdivision_y = 1;
    opt->box_subdivision_z = 1;
    opt->box_spherical_distortion_corr = 0;
    opt->box_length_x = 1.0;
    opt->box_length_y = 1.0;
    opt->box_length_z = 1.0;
    break;
  case RHEA_DOMAIN_BOX:
    /* set all values from input */
    opt->box_subdivision_x = rhea_domain_box_subdivision_x;
    opt->box_subdivision_y = rhea_domain_box_subdivision_y;
    opt->box_subdivision_z = rhea_domain_box_subdivision_z;
    opt->box_spherical_distortion_corr = 0;
    if (isfinite (rhea_domain_box_length_x) && 0.0 < rhea_domain_box_length_x) {
      opt->box_length_x = rhea_domain_box_length_x;
    }
    else { /* if invalid, set default value */
      opt->box_length_x = (double) opt->box_subdivision_x /
                          (double) opt->box_subdivision_z;
    }
    if (isfinite (rhea_domain_box_length_y) && 0.0 < rhea_domain_box_length_y) {
      opt->box_length_y = rhea_domain_box_length_y;
    }
    else { /* if invalid, set default value */
      opt->box_length_y = (double) opt->box_subdivision_y /
                          (double) opt->box_subdivision_z;
    }
    if (isfinite (rhea_domain_box_length_z) && 0.0 < rhea_domain_box_length_z) {
      opt->box_length_z = rhea_domain_box_length_z;
    }
    else { /* if invalid, set default value */
      opt->box_length_z = 1.0;
    }
    break;
  case RHEA_DOMAIN_SHELL:
    /* set default values */
    opt->box_subdivision_x = 0;
    opt->box_subdivision_y = 0;
    opt->box_subdivision_z = 0;
    opt->box_spherical_distortion_corr = 0;
    opt->box_length_x = NAN;
    opt->box_length_y = NAN;
    opt->box_length_z = NAN;
    break;
  case RHEA_DOMAIN_CUBE_SPHERICAL:
    /* set default values */
    opt->box_subdivision_x = 1;
    opt->box_subdivision_y = 1;
    opt->box_subdivision_z = 1;
    opt->box_spherical_distortion_corr = 0;
    opt->box_length_x = NAN;
    opt->box_length_y = NAN;
    opt->box_length_z = NAN;
    break;
  case RHEA_DOMAIN_BOX_SPHERICAL:
    /* set partially from input */
    RHEA_ASSERT (rhea_domain_box_subdivision_y <=
                 rhea_domain_box_subdivision_x);
    RHEA_ASSERT (rhea_domain_box_subdivision_x <=
                 rhea_domain_box_subdivision_z * 6);
    opt->box_subdivision_x = rhea_domain_box_subdivision_x;
    opt->box_subdivision_y = rhea_domain_box_subdivision_y;
    opt->box_subdivision_z = rhea_domain_box_subdivision_z;
    opt->box_spherical_distortion_corr =
      rhea_domain_box_spherical_distortion_corr;
    opt->box_length_x = NAN;
    opt->box_length_y = NAN;
    opt->box_length_z = NAN;
    break;
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* set dimensional domain properties */
  opt->radius_min_m = SC_MAX (0.0, rhea_domain_earth_radius_m -
                                   rhea_domain_mantle_depth_m);
  opt->radius_max_m = rhea_domain_earth_radius_m;
  scal_len_to_km = opt->radius_max_m / 1.0e3;

  /* compute domain properties */
  rhea_domain_compute_bounds (opt);
  rhea_domain_compute_representative_depth (opt);
  rhea_domain_compute_volume (opt);
  rhea_domain_compute_center_of_mass (opt);
  rhea_domain_compute_moment_of_inertia (opt);
  RHEA_ASSERT (isfinite (opt->radius_min) && 0.0 <= opt->radius_min);
  RHEA_ASSERT (isfinite (opt->radius_max) && opt->radius_min < opt->radius_max);

  /* set lower-upper mantle interface */
  radius = 1.0 - rhea_domain_lm_um_interface_depth_m /
                 rhea_domain_earth_radius_m;
  radius = opt->radius_max * SC_MIN (radius, 1.0);
  if (0.0 <= rhea_domain_lm_um_interface_depth_m &&
      opt->radius_min < radius) { /* if lower mantle exists */
    /* align interface with mesh elements */
    opt->lm_um_interface_radius = rhea_domain_align_radius_with_mesh (
        radius, rhea_discretization_level_min, opt);
  }
  else { /* otherwise upper mantle everywhere */
    opt->lm_um_interface_radius = -1.0;
  }

  /* set smoothing of lower-upper mantle interface (if it exists) */
  if (0.0 < rhea_domain_lm_um_interface_smoothing_width_m &&
      opt->radius_min < opt->lm_um_interface_radius &&
      opt->lm_um_interface_radius < opt->radius_max) {
    opt->lm_um_interface_smoothing_width =
      rhea_domain_lm_um_interface_smoothing_width_m /
      rhea_domain_earth_radius_m;
  }
  else { /* otherwise deactivate smoothing */
    opt->lm_um_interface_smoothing_width = 0.0;
  }

  /* set velocity boundary conditions */
  opt->velocity_bc_type =
    (rhea_domain_velocity_bc_t) rhea_domain_velocity_bc_type;

  /* print domain properties */
  RHEA_GLOBAL_INFO ("========================================\n");
  RHEA_GLOBAL_INFOF ("%s\n", __func__);
  RHEA_GLOBAL_INFO ("----------------------------------------\n");
  RHEA_GLOBAL_INFOF ("  shape:                    %s (ID: %i)\n",
                     rhea_domain_shape_name, opt->shape);
  RHEA_GLOBAL_INFOF ("  box subdivision (x,y,z):  (%i %i %i)\n",
                     opt->box_subdivision_x, opt->box_subdivision_y,
                     opt->box_subdivision_z);
  RHEA_GLOBAL_INFOF ("  box length (x,y,z):       (%g %g %g)\n",
                     opt->box_length_x, opt->box_length_y, opt->box_length_z);
  RHEA_GLOBAL_INFOF ("  box_spherical distortion correction:  %i\n",
                     opt->box_spherical_distortion_corr);
  RHEA_GLOBAL_INFO ("----------------------------------------\n");
  RHEA_GLOBAL_INFOF ("  x min:      %g\n", opt->x_min);
  RHEA_GLOBAL_INFOF ("  x max:      %g\n", opt->x_max);
  RHEA_GLOBAL_INFOF ("  y min:      %g\n", opt->y_min);
  RHEA_GLOBAL_INFOF ("  y max:      %g\n", opt->y_max);
  RHEA_GLOBAL_INFOF ("  z min:      %g\n", opt->z_min);
  RHEA_GLOBAL_INFOF ("  z max:      %g\n", opt->z_max);
  RHEA_GLOBAL_INFOF ("  lon min:    %+g\n", opt->lon_min);
  RHEA_GLOBAL_INFOF ("  lon max:    %+g\n", opt->lon_max);
  RHEA_GLOBAL_INFOF ("  radius min: %g (%g km)\n",
                     opt->radius_min, opt->radius_min_m/1.0e3);
  RHEA_GLOBAL_INFOF ("  radius max: %g (%g km)\n",
                     opt->radius_max, opt->radius_max_m/1.0e3);
  RHEA_GLOBAL_INFOF ("  depth:      %g (%g km)\n",
                     opt->depth, opt->depth * scal_len_to_km);
  if (0.0 < opt->lm_um_interface_radius) {
    const double        r_lm_um = opt->lm_um_interface_radius;
    const double        r_max = opt->radius_max;

    RHEA_GLOBAL_INFOF ("  LM-UM interface radius:    %g (%g km)\n",
                       r_lm_um, r_lm_um * scal_len_to_km);
    RHEA_GLOBAL_INFOF ("  LM-UM interface depth:     %g (%g km)\n",
                       (r_max - r_lm_um), (r_max - r_lm_um) * scal_len_to_km);
  }
  else {
    RHEA_GLOBAL_INFOF ("  LM-UM interface radius:    %g\n",
                       opt->lm_um_interface_radius);
  }
  RHEA_GLOBAL_INFOF ("  LM-UM interface smoothing: %g (%g km)\n",
                     opt->lm_um_interface_smoothing_width,
                     opt->lm_um_interface_smoothing_width * scal_len_to_km);
  RHEA_GLOBAL_INFOF ("  volume:            %g\n", opt->volume);
  RHEA_GLOBAL_INFOF ("  center of mass:    %g, %g, %g\n",
                     opt->center[0], opt->center[1], opt->center[2]);
  RHEA_GLOBAL_INFOF ("  moment of inertia: %g, %g, %g\n",
                     opt->moment_of_inertia[0],  opt->moment_of_inertia[1],
                     opt->moment_of_inertia[2]);
  RHEA_GLOBAL_INFO ("========================================\n");
}

/******************************************************************************
 * Basic Calculations
 *****************************************************************************/

double
rhea_domain_depth_m_to_depth (const double depth_m,
                              rhea_domain_options_t *opt)
{
  return opt->radius_max * depth_m / opt->radius_max_m;
}

double
rhea_domain_depth_m_to_radius (const double depth_m,
                               rhea_domain_options_t *opt)
{
  return opt->radius_max * (1.0 - depth_m / opt->radius_max_m);
}

double
rhea_domain_radius_to_radius_m (const double radius,
                                rhea_domain_options_t *opt)
{
  return opt->radius_max_m * radius / opt->radius_max;
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
      const double        z_max = opt->z_max;
      const double        radius_min = opt->radius_min;
      const double        radius_max = opt->radius_max;

      return z / z_max * (radius_max - radius_min) + radius_min;
    }
  case RHEA_DOMAIN_SHELL:
  case RHEA_DOMAIN_CUBE_SPHERICAL:
  case RHEA_DOMAIN_BOX_SPHERICAL:
    return sqrt (x * x + y * y + z * z);
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
  //RHEA_ASSERT ((radii_sum / 8.0) < opt->radius_max); //TODO does not work for custom X-fnc's
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

void
rhea_domain_project_to_surface (double *x, double *y, double *z,
                                rhea_domain_options_t *opt)
{
  const double        radius_max = opt->radius_max;
  double              radius;

  /* compute radius */
  radius = rhea_domain_compute_radius (*x, *y, *z, opt);

  /* map coordinates onto surface parallel to radius */
  switch (opt->shape) {
  case RHEA_DOMAIN_CUBE:
  case RHEA_DOMAIN_BOX:
    *z *= radius_max / radius;
    break;
  case RHEA_DOMAIN_SHELL:
  case RHEA_DOMAIN_CUBE_SPHERICAL:
  case RHEA_DOMAIN_BOX_SPHERICAL:
    *x *= radius_max / radius;
    *y *= radius_max / radius;
    *z *= radius_max / radius;
    break;
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }
}

/******************************************************************************
 * Domain Boundary
 *****************************************************************************/

/**
 * Assigns boundary faces of to p4est-tree faces.
 */
static ymir_topidx_t *
rhea_domain_get_tree_to_bf (p4est_connectivity_t *conn,
                            rhea_domain_shape_t domain_shape)
{
  const ymir_topidx_t n_trees = conn->num_trees;
  ymir_topidx_t       treeid;
  int                 faceid;
  ymir_topidx_t      *tree_to_bf = RHEA_ALLOC (ymir_topidx_t, 6 * n_trees);

  switch (domain_shape) {
  case RHEA_DOMAIN_CUBE:
  case RHEA_DOMAIN_CUBE_SPHERICAL:
  case RHEA_DOMAIN_BOX:
  case RHEA_DOMAIN_BOX_SPHERICAL:
    for (treeid = 0; treeid < n_trees; treeid++) { /* loop over all trees */
      /* set sides */
      for (faceid = 0; faceid < 4; faceid++) { /* loop over all side faces */
        if (conn->tree_to_tree[6 * treeid + faceid] == treeid) {
          switch (faceid) {
          case 0:
            tree_to_bf[6 * treeid + faceid] = RHEA_DOMAIN_BOUNDARY_FACE_SIDE1;
            break;
          case 1:
            tree_to_bf[6 * treeid + faceid] = RHEA_DOMAIN_BOUNDARY_FACE_SIDE2;
            break;
          case 2:
            tree_to_bf[6 * treeid + faceid] = RHEA_DOMAIN_BOUNDARY_FACE_SIDE3;
            break;
          case 3:
            tree_to_bf[6 * treeid + faceid] = RHEA_DOMAIN_BOUNDARY_FACE_SIDE4;
            break;
          default: /* unknown face type */
            RHEA_ABORT_NOT_REACHED ();
          }
        }
        else {
          tree_to_bf[6 * treeid + faceid] = RHEA_DOMAIN_BOUNDARY_FACE_NONE;
        }
      }
      /* set top and bottom */
      tree_to_bf[6 * treeid + 4] = RHEA_DOMAIN_BOUNDARY_FACE_BASE;
      tree_to_bf[6 * treeid + 5] = RHEA_DOMAIN_BOUNDARY_FACE_TOP;
    }
    break;

  case RHEA_DOMAIN_SHELL:
    for (treeid = 0; treeid < n_trees; treeid++) { /* loop over all trees */
      /* set sides */
      for (faceid = 0; faceid < 4; faceid++) { /* loop over all side faces */
        tree_to_bf[6 * treeid + faceid] = RHEA_DOMAIN_BOUNDARY_FACE_NONE;
      }
      /* set top and bottom */
      tree_to_bf[6 * treeid + 4] = RHEA_DOMAIN_BOUNDARY_FACE_BASE;
      tree_to_bf[6 * treeid + 5] = RHEA_DOMAIN_BOUNDARY_FACE_TOP;
    }
    break;

  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }

  return tree_to_bf;
}

rhea_domain_boundary_t *
rhea_domain_boundary_new (p4est_t *p4est, rhea_domain_options_t *opt)
{
  rhea_domain_boundary_t  *boundary;

  /* create new object */
  boundary = RHEA_ALLOC (rhea_domain_boundary_t, 1);

  /* set element to boundary function */
  boundary->e_to_fm_fn = ymir_mesh_e_to_fm_tag;

  /* set boundary faces */
  boundary->tree_to_bf = rhea_domain_get_tree_to_bf (p4est->connectivity,
                                                     opt->shape);

  /* return boundary */
  return boundary;
}

void
rhea_domain_boundary_destroy (rhea_domain_boundary_t *boundary)
{
  RHEA_FREE (boundary->tree_to_bf);
  RHEA_FREE (boundary);
}

/******************************************************************************
 * Boundary Conditions
 *****************************************************************************/

/* user provided velocity boundary conditions */
ymir_vel_dir_fn_t   rhea_domain_velocity_bc_user_fn = NULL;
void               *rhea_domain_velocity_bc_user_data = NULL;
int                 rhea_domain_velocity_bc_user_nonzero = -1;

/**
 *
 */
static ymir_dir_code_t
rhea_domain_vel_dir_fn_box (double X, double Y, double Z,
                            double nx, double ny, double nz,
                            ymir_topidx_t face, ymir_locidx_t node_id,
                            void *data)
{
  int                 vel_bc_type = *((int *) data);

  switch (vel_bc_type) {
  case RHEA_DOMAIN_VELOCITY_BC_DIRICHLET_ALL:
    /* set Dirichlet in all directions for all points */
    return YMIR_VEL_DIRICHLET_ALL;
  case RHEA_DOMAIN_VELOCITY_BC_DIRICHLET_NORM:
    /* set Dirichlet in normal direction for all points */
    return YMIR_VEL_DIRICHLET_NORM;
  case RHEA_DOMAIN_VELOCITY_BC_DIR_SIDES_NEU_TB:
    if ( (face != RHEA_DOMAIN_BOUNDARY_FACE_BASE) &&
         (face != RHEA_DOMAIN_BOUNDARY_FACE_TOP) ) {
      /* set Dirichlet in all directions on side faces */
      return YMIR_VEL_DIRICHLET_ALL;
    }
    else {
      /* bottom and top faces will get Neumann BC's */
      return YMIR_VEL_DIRICHLET_NONE;
    }
  case RHEA_DOMAIN_VELOCITY_BC_DIR_NORM_SIDES_DIR_ALL_TB:
    if ( (face != RHEA_DOMAIN_BOUNDARY_FACE_BASE) &&
         (face != RHEA_DOMAIN_BOUNDARY_FACE_TOP) ) {
      /* set Dirichlet in normal direction on side faces */
      return YMIR_VEL_DIRICHLET_NORM;
    }
    else {
      /* set Dirichlet in all directions on bottom and top faces */
      return YMIR_VEL_DIRICHLET_ALL;
    }
  case RHEA_DOMAIN_VELOCITY_BC_NEUMANN_ALL:
    /* set Neumann for all points */
    return YMIR_VEL_DIRICHLET_NONE;
  default: /* unknown boundary condition */
    RHEA_ABORT_NOT_REACHED ();
  }
}

/**
 *
 */
static ymir_dir_code_t
rhea_domain_vel_dir_fn_shell (double X, double Y, double Z,
                              double nx, double ny, double nz,
                              ymir_topidx_t face, ymir_locidx_t node_id,
                              void *data)
{
  int                 vel_bc_type = *((int *) data);

  switch (vel_bc_type) {
  case RHEA_DOMAIN_VELOCITY_BC_DIRICHLET_ALL:
    /* set Dirichlet BC in all directions for all points */
    return YMIR_VEL_DIRICHLET_ALL;
  case RHEA_DOMAIN_VELOCITY_BC_DIRICHLET_NORM:
  case RHEA_DOMAIN_VELOCITY_BC_DIRICHLET_NORM_FIXDOF:
    /* set Dirichlet BC in normal direction for all points */
    return YMIR_VEL_DIRICHLET_NORM;
  case RHEA_DOMAIN_VELOCITY_BC_DIRICHLET_ALL_INNER:
    if (face == RHEA_DOMAIN_BOUNDARY_FACE_TOP) {
      /* set Dirichlet BC in normal direction on outer sphere of a shell */
      return YMIR_VEL_DIRICHLET_NORM;
    }
    else {
      /* set Dirichlet BC in all direction on inner sphere of a shell */
      return YMIR_VEL_DIRICHLET_ALL;
    }
  case RHEA_DOMAIN_VELOCITY_BC_NEUMANN_ALL:
    /* set Neumann for all points */
    return YMIR_VEL_DIRICHLET_NONE;
  default: /* unknown boundary condition */
    RHEA_ABORT_NOT_REACHED ();
  }
}

/**
 * Callback function for `slabs_shell_vel_dir_restrict_dof` below.
 */
static void
rhea_domain_vel_dir_shell_restrict_dof_fn (double *dummy,
                                           double X, double Y, double Z,
                                           double nx, double ny, double nz,
                                           ymir_topidx_t face,
                                           ymir_locidx_t node_id,
                                           void *data)
{
  ymir_vel_dir_t     *vel_dir = data;

  /* set one point to all Dirichlet */
  if (face == RHEA_DOMAIN_BOUNDARY_FACE_BASE && fabs (nz - 1.0) < SC_EPS) {
    ymir_mesh_t        *ymir_mesh = vel_dir->mesh;
    ymir_face_mesh_t   *fmesh = &ymir_mesh->fmeshes[face];
    ymir_locidx_t       vcnid = fmesh->CntoVCn[node_id];
    double             *vn = ymir_cvec_index (vel_dir->nvec, vcnid, 0);

    /* set Dirichlet in all directions */
    vel_dir->rank[vcnid] = YMIR_VEL_DIRICHLET_ALL;
    vn[0] = 0;
    vn[1] = 0;
    vn[2] = 0;

    RHEA_INFOF ("Set boundary node %d with coord. (%g %g %g) "
                "and normal (%g %g %g) to all Dirichlet\n",
                node_id, X, Y, Z, nx, ny, nz);
  }

  /* set another point to Dirichlet in normal and one tangential direction */
  if (face == RHEA_DOMAIN_BOUNDARY_FACE_BASE && fabs (ny - 1.0) < SC_EPS) {
    ymir_mesh_t        *ymir_mesh = vel_dir->mesh;
    ymir_face_mesh_t   *fmesh = &ymir_mesh->fmeshes[face];
    ymir_locidx_t       vcnid = fmesh->CntoVCn[node_id];
    double             *vn = ymir_cvec_index (vel_dir->nvec, vcnid, 0);

    /* set the way the vector `vn` is interpreted */
    vel_dir->rank[vcnid] = YMIR_VEL_DIRICHLET_TANG;

    /* set Dirichlet in y- and z-direction */
    vn[0] = 1;
    vn[1] = 0;
    vn[2] = 0;

    RHEA_INFOF ("Set boundary node %d with coord. (%g %g %g) "
                "and normal (%g %g %g) to Dirichlet in normal "
                "and one tangential direction\n",
                node_id, X, Y, Z, nx, ny, nz);
  }
}

/**
 *
 */
static void
rhea_domain_vel_dir_shell_restrict_dof (ymir_vel_dir_t * vel_dir)
{
  ymir_mesh_t        *ymir_mesh = vel_dir->mesh;
  const ymir_topidx_t n_faces = ymir_mesh->num_face_meshes;
  ymir_topidx_t       faceid;

  for (faceid = 0; faceid < n_faces; faceid++) {
    ymir_vec_t         *dummy_face_cvec;

    /* create dummy vector for current face */
    dummy_face_cvec = ymir_face_cvec_new (ymir_mesh, faceid, 1);

    /* set additional Dirichlet boundary conditions for current face */
    ymir_face_cvec_set_function (
        dummy_face_cvec, rhea_domain_vel_dir_shell_restrict_dof_fn, vel_dir);

    /* destroy vector */
    ymir_vec_destroy (dummy_face_cvec);
  }
}

void
rhea_domain_set_user_velocity_dirichlet_bc (ymir_vel_dir_fn_t vel_dir_bc_fn,
                                            void *vel_dir_bc_data,
                                            const int vel_dir_bc_nonzero)
{
  rhea_domain_velocity_bc_user_fn = vel_dir_bc_fn;
  rhea_domain_velocity_bc_user_data = vel_dir_bc_data;
  rhea_domain_velocity_bc_user_nonzero = vel_dir_bc_nonzero;
}

ymir_vel_dir_t *
rhea_domain_create_velocity_dirichlet_bc (ymir_mesh_t *ymir_mesh,
                                          ymir_vec_t *dirscal,
                                          void *data)
{
  rhea_domain_options_t  *opt = data;
  int                 vel_bc_type = (int) opt->velocity_bc_type;
  ymir_vel_dir_t     *vel_dir;

  /* call user provided function and return */
  if (RHEA_DOMAIN_VELOCITY_BC_USER == vel_bc_type) {
    RHEA_CHECK_ABORT (
        rhea_domain_velocity_bc_user_fn != NULL &&
        rhea_domain_velocity_bc_user_nonzero != -1,
        "Function for defining velocity boundary conditions was not provided");

    vel_dir = ymir_vel_dir_new_from_faces (
        ymir_mesh, NULL, rhea_domain_velocity_bc_user_fn,
        rhea_domain_velocity_bc_user_data, dirscal);
    vel_dir->nonzero = rhea_domain_velocity_bc_user_nonzero;
    return vel_dir;
  }

  /* create velocity Dirichlet BC's */
  switch (opt->shape) {
  case RHEA_DOMAIN_CUBE:
  case RHEA_DOMAIN_CUBE_SPHERICAL:
  case RHEA_DOMAIN_BOX:
  case RHEA_DOMAIN_BOX_SPHERICAL:
    vel_dir = ymir_vel_dir_new_from_faces (
        ymir_mesh, NULL, rhea_domain_vel_dir_fn_box, &vel_bc_type, dirscal);
    break;
  case RHEA_DOMAIN_SHELL:
    vel_dir = ymir_vel_dir_new_from_faces (
        ymir_mesh, NULL, rhea_domain_vel_dir_fn_shell, &vel_bc_type, dirscal);
    if (RHEA_DOMAIN_VELOCITY_BC_DIRICHLET_NORM_FIXDOF == vel_bc_type) {
      rhea_domain_vel_dir_shell_restrict_dof (vel_dir);
    }
    break;
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* set zero Dirichlet boundary */
  vel_dir->nonzero = 0;

  /* return BC's */
  return vel_dir;
}
