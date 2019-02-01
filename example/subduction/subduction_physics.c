#include <subduction_physics.h>

/**************************************
 * Temperature Computation
 *************************************/
double
subd_temperature_brick_2plates_poly2 (double r, double lon,
                                       subd_options_t * subd_options)
{
  double              trench_lon;
  double              subd_depth, subd_width;
  double              subd_edge_width, subd_edge_smoothwidth;
  double              subd_plate_vel, subd_plate_init_age;
  double              over_plate_age;

  double              y_min = subd_options->domain_options->y_min;
  double              start_node, start_val, start_deriv;
  double              end_node, end_val;
  double             *poly2_coeff;
  int                 orientation_wrt_curve;
  double             *closest_pt;

  double              depth;
  double              ridge_dist;
  double              dist;
  double              vel;
  double              plate_time;
  double              subd_length;
  double              temp_subd = 1.0;
  double              temp_tip;
  double              temp;

  subd_temp_2plates_slab_t * slab_options = subd_options->temp_options->slab_options;
  subd_2plates_poly2_geo_coeff_t *geo = slab_options->temp_2plates_geo_coeff;

  /* set parameters according to temperature options */
  trench_lon = slab_options->temp_2plates_trench_longitude;
  subd_depth = slab_options->temp_2plates_subd_depth;
  subd_width = slab_options->temp_2plates_subd_width;
  subd_edge_width = slab_options->temp_2plates_subd_edge_width;
  subd_edge_smoothwidth = slab_options->temp_2plates_subd_edge_smoothwidth;
  subd_plate_vel = slab_options->temp_2plates_subd_plate_velocity;
  subd_plate_init_age = slab_options->temp_2plates_subd_plate_initial_age;
  over_plate_age = slab_options->temp_2plates_over_plate_age;

  /* check parameters */
  RHEA_ASSERT (0.0 < trench_lon);
  RHEA_ASSERT (0.0 < subd_depth && 0.0 < subd_width);
  RHEA_ASSERT (0.0 <= subd_edge_width && 0.0 <= subd_edge_smoothwidth);
  RHEA_ASSERT (0.0 < subd_plate_vel);
  RHEA_ASSERT (0.0 <= subd_plate_init_age);
  RHEA_ASSERT (0.0 < over_plate_age);

  start_node = geo->start_node;
  start_val = geo->start_val;
  start_deriv = geo->start_deriv;
  end_node = geo->end_node;
  end_val = geo->end_val;
  poly2_coeff = geo->poly2_coeff;

  /* compute closest point on curve and orientation w.r.t. curve */
  closest_pt = subd_compute_closest_pt_on_poly2 (lon, r,
                                                  poly2_coeff, start_node,
                                                  start_val, start_deriv,
                                                  end_node, end_val,
                                                  &orientation_wrt_curve);


  /* compute "fake" (TODO) distance */
  dist = SUBD_MANTLE_DEPTH * sqrt ( (lon - closest_pt[0]) * (lon - closest_pt[0])
                                  + (r - closest_pt[1]) * (r - closest_pt[1]) );

  /* For computing the temperature, use temperature profile of the plates
   * from halfspace cooling model, erf(z/(2*sqrt(t*kappa))), as base
   * temperature field. */

  /* compute temperature of plates at surface */
  if (orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_BOTTOM ||
      orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_BOTTOM_BACK) { /* if subd plate */
    /* calculate age of subducting plate */
    ridge_dist = fabs (lon - y_min) * SUBD_MANTLE_DEPTH;
    vel = subd_plate_vel / SUBD_SEC_PER_YEAR;
    plate_time = ridge_dist / vel + subd_plate_init_age * SUBD_SEC_PER_YEAR;

    /* avoid division by zero */
    plate_time = SC_MAX (1.0e-3, plate_time);
  }
  else  { /* if overriding plate */
    /* calculate age of overriding plate */
    plate_time = over_plate_age * SUBD_SEC_PER_YEAR;
  }
  depth = (SUBD_SHELL_RADIUS_TOP - r) * SUBD_MANTLE_DEPTH;
  temp = erf ( depth / (2.0 * sqrt (plate_time * SUBD_REF_THERM_DIFFUS)) );

  /* compute temperature of the plate subducted inside of mantle */
  if (   orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_BOTTOM
      || orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_BOTTOM_BACK
      || orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_BOTTOM_LEFT ) {
    /* calculate total lenght of subducted plate */
    subd_length = SUBD_MANTLE_DEPTH * sqrt (
                    (start_node - end_node) * (start_node - end_node)
                    + (start_val - end_val) * (start_val - end_val) );

    /* calculate age of plate */
    if (lon <= trench_lon) {
      ridge_dist = fabs (lon - y_min) * SUBD_MANTLE_DEPTH;
    }
    else {
      ridge_dist = fabs (trench_lon - y_min) * SUBD_MANTLE_DEPTH
                   + subd_length * depth / subd_depth;
    }
    vel = subd_plate_vel / SUBD_SEC_PER_YEAR;
    plate_time = ridge_dist / vel + subd_plate_init_age * SUBD_SEC_PER_YEAR;

    /* avoid division by zero */
    plate_time = SC_MAX (1.0e-3, plate_time);

    /* compute temperature of subducting plate */
    temp_subd = erf ( dist / (2.0 * sqrt (plate_time * SUBD_REF_THERM_DIFFUS)) );
  }

  /* smooth the subducting plate's top edge */
  if (   orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_TOP
      || orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_TOP_RIGHT ) {
    if (dist < subd_edge_width) {
      /* set constant temperature inside edge width */
      temp_subd = 0.0;
    }
    else {
      /* smooth edge with Gaussian */
      temp_subd = 1.0 - exp (
        - (dist - subd_edge_width) * (dist - subd_edge_width)
        / (0.5 * subd_edge_smoothwidth * subd_edge_smoothwidth) );
    }
  }

  /* smooth tip of subducting plate */
  if (   orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_TOP_RIGHT
      || orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_BOTTOM_LEFT ) {
    /* compute distance between tip and closest point on curve */
    dist = SUBD_MANTLE_DEPTH * sqrt (
      (end_node - closest_pt[0]) * (end_node - closest_pt[0])
      + (end_val - closest_pt[1]) * (end_val - closest_pt[1]) );

    /* compute temperature */
    temp_tip = 1.0 -
      exp (- dist*dist / (0.5 * subd_edge_smoothwidth*subd_edge_smoothwidth) );
    temp_subd = SC_MAX (temp_subd, temp_tip);
  }

  /* update global temperature */
  temp = SC_MIN (temp, temp_subd);

  /* destroy */
  RHEA_FREE (closest_pt);

  /* return temperature */
  RHEA_ASSERT (isfinite (temp));
  return temp;
}

void
subd_poly2_temperature_set_fn (double *temp, double x, double y, double z,
                          ymir_locidx_t nid, void *data)
{
  subd_options_t    *subd_options = (subd_options_t *) data;
  double              lon, r;

  /* compute radius and longitude */
  lon = x;
  r = z;

  *temp = subd_temperature_brick_2plates_poly2 (r, lon, subd_options);

  /* check temperature for `nan` and `inf` */
  RHEA_ASSERT (isfinite (*temp));

  /* bound temperature to valid interval */
  *temp = SC_MIN (1.0, *temp);
  *temp = SC_MAX (0.0, *temp);
}

void
subd_poly2_temperature_compute (ymir_vec_t *temperature,
                           subd_options_t *subd_options)
{
  const char         *this_fn_name = "subd_poly2_temperature_compute";
  double              trench_lon;
  double              dip_angle;
  double              subd_depth, subd_width;
  double              start_node, start_val, start_deriv, end_node, end_val;
  double              *poly2_coeff;
  subd_2plates_poly2_geo_coeff_t geo;
  subd_temp_2plates_slab_t * slab_options = subd_options->temp_options->slab_options;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* set parameters according to physics options */
  trench_lon = slab_options->temp_2plates_trench_longitude;
  dip_angle = slab_options->temp_2plates_dip_angle;
  subd_depth = slab_options->temp_2plates_subd_depth;
  subd_width = slab_options->temp_2plates_subd_width;

  /* check parameters */
  RHEA_ASSERT (0.0 < trench_lon);
  RHEA_ASSERT (0.0 < dip_angle && dip_angle < 90.0);
  RHEA_ASSERT (0.0 < subd_depth && 0.0 < subd_width);

  /* set points for polynomial interpolation */
  start_node = trench_lon;
  start_val = SUBD_SHELL_RADIUS_TOP;
  start_deriv = tan (-dip_angle / 180.0 * M_PI);
  end_node = start_node + subd_width / SUBD_MANTLE_DEPTH;
  end_val = start_val - subd_depth / SUBD_MANTLE_DEPTH;
  /* compute interpolating quadratic polynomial */
  poly2_coeff = subd_compute_poly2_interpolation (start_node, start_val,
                                                   start_deriv,
                                                   end_node, end_val);

  RHEA_GLOBAL_PRODUCTIONF("temperature: brick_2plates_poly2 downgoing subd geometry\nstart_node: (%f %f), start_deriv: %f, end_node: (%f %f)\npoly2_coeff: (%f %f %f)\n", start_node, start_val, start_deriv, end_node, end_val, poly2_coeff[2], poly2_coeff[1], poly2_coeff[0]);

  geo.start_node = start_node;
  geo.start_val = start_val;
  geo.start_deriv = start_deriv;
  geo.end_node = end_node;
  geo.end_val = end_val;
  geo.poly2_coeff = poly2_coeff;
  slab_options->temp_2plates_geo_coeff = &geo;

  ymir_cvec_set_function (temperature,
                          subd_poly2_temperature_set_fn,
                          subd_options);
  RHEA_FREE (poly2_coeff);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}


/*half space cooling model*/
static double
subd_temperature_cold_plate_hscm_plate_age_to_scaling (
                                 const double radius_max_m,
                                 const double plate_age_yr,
                                 const double thermal_diffus_m2_s)
{
    const double        plate_age_s =
                         plate_age_yr * SUBD_SEC_PER_YEAR;

      return radius_max_m / (2.0 * sqrt (plate_age_s * thermal_diffus_m2_s));
}

double
subd_temperature_cold_plate_hscm (const double radius,
                                  const double radius_max,
                                  const double radius_max_m,
                                  const double plate_age_yr,
                                  const double thermal_diffus_m2_s)
{
    double              scaling_nondim;

    scaling_nondim = subd_temperature_cold_plate_hscm_plate_age_to_scaling (
                radius_max_m, plate_age_yr, thermal_diffus_m2_s);

    return erf (scaling_nondim * (radius_max - radius));
}


/* customer designed temperture field, mostly for case test purposes*/
static void
thinbox_temperature_set_fn (double *temp, double x, double y, double z,
                          ymir_locidx_t nid, void *data)
{
  if ( (fabs(y-0.5) < 0.2)  && (z > 0.95 && z < 0.97))
    *temp = 0.0;
  else
    *temp = 0.5;

  /* check temperature for `nan` and `inf` */
  RHEA_ASSERT (isfinite (*temp));

  /* bound temperature to valid interval */
  *temp = SC_MIN (1.0, *temp);
  *temp = SC_MAX (0.0, *temp);
}

static void
sinker_temperature_set_fn (double *temp, double x, double y, double z,
                          ymir_locidx_t nid, void *data)
{
  subd_options_t *opt = data;
  subd_temp_custom_sinker_t *sinker = opt->temp_options->sinker_options;
  double center_x = sinker->center_x;
  double center_y = sinker->center_y;
  double center_z = sinker->center_z;
  double radius = 0.5 * sinker->diameter;
  double scaling = sinker->scaling;
  double decay = sinker->decay;
  double val;

  val = sqrt ( (center_x - x) * (center_x - x) +
               (center_y - y) * (center_y - y) +
               (center_z - z) * (center_z - z) );

  val = SC_MAX (0.0, val - radius);

  val = exp (-decay * val * val);
  RHEA_ASSERT (0.0 <= val && val <= 1.0);

  val = scaling * (1.0 - val);
  val = SC_MAX (0.0, SC_MIN (val, 1.0));
  *temp = val;
  *temp += (0.5 - scaling);

  /* check temperature for `nan` and `inf` */
  RHEA_ASSERT (isfinite (*temp));

  /* bound temperature to valid interval */
  *temp = SC_MIN (1.0, *temp);
  *temp = SC_MAX (0.0, *temp);
}

static void
drag_temperature_set_fn (double *temp, double x, double y, double z,
                          ymir_locidx_t nid, void *data)
{
  subd_options_t    *subd_options = (subd_options_t *) data;
  double              subdu_lon, subdu_dip_angle, subdu_depth, subdu_width;
  double              start_node, start_val, start_deriv, end_node, end_val;
  subd_weak_options_t *weak_options = subd_options->weak_options;

  /* set parameters according to weakzone options */
  subdu_lon = weak_options->weakzone_2plates_subdu_longitude;
  subdu_dip_angle = weak_options->weakzone_2plates_subdu_dip_angle;
  subdu_depth = weak_options->weakzone_2plates_subdu_depth;
  subdu_width = weak_options->weakzone_2plates_subdu_width;

  /* check parameters */
  RHEA_ASSERT (0.0 < subdu_lon);
  RHEA_ASSERT (0.0 < subdu_dip_angle);
  RHEA_ASSERT (0.0 < subdu_depth && 0.0 < subdu_width);

  /* set points for polynomial interpolation */
  start_node = subdu_lon;
  start_val = SUBD_SHELL_RADIUS_TOP;
  end_node = start_node + subdu_width / SUBD_MANTLE_DEPTH;
  end_val = start_val - subdu_depth / SUBD_MANTLE_DEPTH;

  if ((y<end_node) && (y>end_node-0.4) && (z<end_val) && (z>end_val-0.2) )
    *temp = 1.0;
  else
    *temp = 0.5;

  /* check temperature for `nan` and `inf` */
  RHEA_ASSERT (isfinite (*temp));

  /* bound temperature to valid interval */
  *temp = SC_MIN (1.0, *temp);
  *temp = SC_MAX (0.0, *temp);
}

static void
lithblock_temperature_set_fn (double *temp, double x, double y, double z,
                          ymir_locidx_t nid, void *data)
{
  subd_options_t    *opt = (subd_options_t *) data;
  subd_temp_custom_lithblock_t *lithblock = opt->temp_options->lithblock_options;
  double r, deg;
  double dist_r, dist_deg;

  dist_r = 0.0141;
  dist_deg = 10.0;
  r = sqrt (x*x + z*z);
  deg = acos(x/r) * 180.0 / M_PI;

  if ( (1.0 - r < dist_r)  && (deg - 67.5 < dist_deg))  {
    *temp = .0;
  }
  else
    *temp = 1.0;

  /* check temperature for `nan` and `inf` */
  RHEA_ASSERT (isfinite (*temp));

  /* bound temperature to valid interval */
  *temp = SC_MIN (1.0, *temp);
  *temp = SC_MAX (0.0, *temp);
}

static void
hscm2block_temperature_set_fn (double *temp, double x, double y, double z,
                              ymir_locidx_t nid, void *data)
{
  subd_options_t    *opt = (subd_options_t *) data;
  subd_temp_custom_hscm2block_t *hscm2block = opt->temp_options->hscm2block_options;
  rhea_domain_options_t *domain_opt = opt->domain_options;
  const double        radius_max = domain_opt->radius_max;
  const double        radius_max_m = domain_opt->radius_max_m;
//  const double        plate1_age_yr = hscm2block->plate1_age_yr;
//  const double        plate2_age_yr = hscm2block->plate2_age_yr;
//  const double        p12bound = hscm2block->plates12_bound;
  double              radius;

  double plate1_age_yr = 1.0e8;
  double plate2_age_yr = 5.0e8;
  double p12bound = 1.0;
  double subsurf = .0;

  radius = rhea_domain_compute_radius (x, y, z, domain_opt);
  if (radius > radius_max - subsurf)  {
    *temp = .0;
  }
  else {
    radius += subsurf;
    if (x < p12bound)  {
      *temp = subd_temperature_cold_plate_hscm (
          radius, radius_max, radius_max_m, plate1_age_yr, SUBD_REF_THERM_DIFFUS);
    }
    else  {
      *temp = subd_temperature_cold_plate_hscm (
          radius, radius_max, radius_max_m, plate2_age_yr, SUBD_REF_THERM_DIFFUS);
    }
  }
/* check temperature for `nan` and `inf` */
RHEA_ASSERT (isfinite (*temp));

/* bound temperature to valid interval */
  *temp = SC_MIN (1.0, *temp);
  *temp = SC_MAX (0.0, *temp);

}

static void
hscm2blocktanh_temperature_set_fn (double *temp, double x, double y, double z,
                              ymir_locidx_t nid, void *data)
{
  subd_options_t    *opt = (subd_options_t *) data;
  subd_temp_custom_hscm2block_t *hscm2block = opt->temp_options->hscm2block_options;
  rhea_domain_options_t *domain_opt = opt->domain_options;
  const double        radius_max = domain_opt->radius_max;
  const double        radius_max_m = domain_opt->radius_max_m;
//  const double        plate1_age_yr = hscm2block->plate1_age_yr;
//  const double        plate2_age_yr = hscm2block->plate2_age_yr;
//  const double        p12bound = hscm2block->plates12_bound;
  double              radius, age;

  double plate1_age_yr = 1.0e8;
  double plate2_age_yr = 5.0e8;
  double p12bound = 1.0;
  double scaling = 100.0;

  age = (plate2_age_yr - plate1_age_yr) / 2.0 * tanh (scaling * (x - 1.0)) +
        (plate2_age_yr + plate1_age_yr) / 2.0;
  radius = rhea_domain_compute_radius (x, y, z, domain_opt);
  *temp = subd_temperature_cold_plate_hscm (
        radius, radius_max, radius_max_m, age, SUBD_REF_THERM_DIFFUS);

/* check temperature for `nan` and `inf` */
RHEA_ASSERT (isfinite (*temp));

/* bound temperature to valid interval */
*temp = SC_MIN (1.0, *temp);
  *temp = SC_MAX (0.0, *temp);

}

static void
hscm1plate_temperature_set_fn (double *temp, double x, double y, double z,
                              ymir_locidx_t nid, void *data)
{
  subd_options_t    *opt = (subd_options_t *) data;
  rhea_domain_options_t *domain_opt = opt->domain_options;
  const double        x_max = domain_opt->x_max;
  const double        radius_max = domain_opt->radius_max;
  const double        radius_max_m = domain_opt->radius_max_m;
  double              radius;
  double              age;
  double plate_age0_yr = 1.0e3;
  double plate_age1_yr = 5.0e8;

  age = x/x_max * plate_age1_yr + plate_age0_yr;
  radius = rhea_domain_compute_radius (x, y, z, domain_opt);

  *temp = subd_temperature_cold_plate_hscm (
          radius, radius_max, radius_max_m, age, SUBD_REF_THERM_DIFFUS);

  /* check temperature for `nan` and `inf` */
  RHEA_ASSERT (isfinite (*temp));

/* bound temperature to valid interval */
  *temp = SC_MIN (1.0, *temp);
  *temp = SC_MAX (0.0, *temp);

}


void
subd_custom_temperature_compute (ymir_vec_t *temp,
                                    subd_options_t *subd_options)
{
  subd_temp_custom_t type = subd_options->temp_options->custom_type;

  switch (type)  {
    case SUBD_TEMP_CUSTOM_THINBOX:
      ymir_cvec_set_function (temp, thinbox_temperature_set_fn, subd_options);
    break;

    case SUBD_TEMP_CUSTOM_SINKER:
      ymir_cvec_set_function (temp, sinker_temperature_set_fn, subd_options);
    break;

    case SUBD_TEMP_CUSTOM_DRAG:
      ymir_cvec_set_function (temp, drag_temperature_set_fn, subd_options);
    break;

    case SUBD_TEMP_CUSTOM_LITHBLOCK:
      ymir_cvec_set_function (temp, lithblock_temperature_set_fn, subd_options);
    break;

    case SUBD_TEMP_CUSTOM_HSCM2BLOCK:
      ymir_cvec_set_function (temp, hscm2block_temperature_set_fn, subd_options);
    break;

    case SUBD_TEMP_CUSTOM_HSCM2BLOCKTANH:
      ymir_cvec_set_function (temp, hscm2blocktanh_temperature_set_fn, subd_options);
    break;

    case SUBD_TEMP_CUSTOM_HSCM1PLATE:
      ymir_cvec_set_function (temp, hscm1plate_temperature_set_fn, subd_options);
    break;
    default:
      RHEA_ABORT_NOT_REACHED ();
  }

}

void
subd_compute_temperature (ymir_vec_t *temp,
                          rhea_temperature_options_t *temp_options,
                          subd_options_t *subd_options)
{
  switch (subd_options->temp_options->type)  {
  /* if non-specified, use: rhea_temperature_compute (temperature, temp_options); */
    case SUBD_TEMP_RHEA:
      rhea_temperature_compute (temp, temp_options);
      break;
    case SUBD_TEMP_SLAB:
      subd_poly2_temperature_compute (temp, subd_options);
      break;
    case SUBD_TEMP_CUSTOM:
      subd_custom_temperature_compute (temp, subd_options);
      break;
    case SUBD_TEMP_NONE:
      ymir_vec_set_value (temp, 0.5);
      break;
    default:
      RHEA_ABORT_NOT_REACHED ();
  }

}

/**************************************
 * Weakzone Computation
 *************************************/

/* Computes the weak zone factor depending on the distance to a weak zone.
 * The edges of the weak zone are smoothed by a Gaussian.
 * 1 - (1 - weak_factor) * exp ( - dist^2 / (2 * (0.5*thickness)^2) )
 * I added a scale of 20.0 so that the smoothing is sharper, at designited width, weak is aoubt exp(-10)~1-0.001~0.999*/
//static inline double
double
subd_weakzone_factor_fn (const double distance,
                          const double thickness,
                          const double thickness_const,
                          const double weak_factor)
{
  const double        d = distance - 0.5 * thickness_const;
  const double        std_dev = 0.5 * (thickness - thickness_const);

//  RHEA_ASSERT (thickness_const <= thickness);

  if (d <= 0.0) {
    /* return value inside zone with constant weak factor */
    return weak_factor;
  }
  else {
    /* return smoothed weak zone */
    return 1.0 - (1.0 - weak_factor) * exp (-20.0 * d*d / (2.0 * std_dev*std_dev));
  }
}

/* Computes distance to weak zone between plates. */
double
subd_weakzone_subduct_dist_2plates_poly2 (double r, double lon,
                                         double *poly2_coeff,
                                         double start_node,
                                         double start_val,
                                         double start_deriv,
                                         double end_node,
                                         double end_val)
{
  int                 orientation_wrt_curve;
  double             *closest_pt;
  double              dist;

  /* compute closest point on curve and orientation w.r.t. curve */
  closest_pt = subd_compute_closest_pt_on_poly2 (lon, r,
                                                  poly2_coeff, start_node,
                                                  start_val, start_deriv,
                                                  end_node, end_val,
                                                  &orientation_wrt_curve);

  /* compute distance */
  if (   orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_TOP
      || orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_BOTTOM
      || orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_BOTTOM_BACK) {
    /* compute distance to closest point on curve */
    dist = SUBD_MANTLE_DEPTH * sqrt (  SC_SQR (lon - closest_pt[0])
                                   + SC_SQR (r - closest_pt[1]) );
  }
  else if (   orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_TOP_RIGHT
           || orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_BOTTOM_LEFT ) {
    /* compute distance to tip of curve */
    dist = SUBD_MANTLE_DEPTH * sqrt (  SC_SQR (end_node - lon)
                                   + SC_SQR (end_val - r) );
  }
  else {
    double              zero = 0.0;  /* no const to avoid warning */
    dist = 1.0 / zero;
    RHEA_ABORT_NOT_REACHED ();
  }

  /* destroy */
  RHEA_FREE (closest_pt);

  /* return distance to weak zone curve */
  return dist;
}

/* Computes distance to weak zone at ridge. */
 double
subd_weakzone_ridge_dist_2plates_poly2 (double r, double lon,
                                         double end_node,
                                         double end_val)
{
  double              dist;

  /* compute distance to weak zone */
  if (lon <= end_node && end_val <= r) { /* if inside weak zone */
    dist = 0.0;
  }
  else { /* if outside weak zone */
    if ( end_node < lon && end_val <= r ) {
      /* compute distance to right edge */
      dist = SUBD_MANTLE_DEPTH * (lon - end_node);
    }
    else if (lon <= end_node && r < end_val) {
      /* compute distance to bottom edge */
      dist = SUBD_MANTLE_DEPTH * (end_val - r);
    }
    else {
      /* compute distance to bottom right corner */
      dist = SUBD_MANTLE_DEPTH * sqrt (  SC_SQR (lon - end_node)
                                     + SC_SQR (end_val - r) );
    }
  }

  /* return distance to ridge weak zone */
  return dist;
}

/* Computes the weak fault zone for the shell slice with two plates. */
//static double
double
subd_weakzone_brick_2plates_poly2 (double r, double lon,
                              subd_options_t * subd_options)
{
  double              subdu_lon;
  double              subdu_dip_angle;
  double              subdu_depth, subdu_width;
  double              subdu_thickness, subdu_thickness_const;
  double              subdu_weak_factor;
  double              ridge_depth, ridge_width;
  double              ridge_smoothwidth;
  double              ridge_weak_factor;

  double              weakzone_thickness, courtesy_width;
  double              y_min = subd_options->domain_options->y_min;
  double              start_node, start_val, start_deriv, end_node, end_val;
  double              *poly2_coeff;
  double              dist;
  double              weak = 1.0;
  double              s, a, b;

  subd_weak_options_t *weak_options = subd_options->weak_options;
  subd_2plates_poly2_geo_coeff_t *geo = weak_options->weak_2plates_geo_coeff;

  /* set parameters according to weakzone options */
  subdu_lon = weak_options->weakzone_2plates_subdu_longitude;
  subdu_dip_angle = weak_options->weakzone_2plates_subdu_dip_angle;
  subdu_depth = weak_options->weakzone_2plates_subdu_depth;
  subdu_width = weak_options->weakzone_2plates_subdu_width;
  subdu_thickness = weak_options->weakzone_2plates_subdu_thickness;
  subdu_thickness_const =
    weak_options->weakzone_2plates_subdu_thickness_const;
  subdu_weak_factor = weak_options->weakzone_2plates_subdu_weak_factor;
  ridge_depth = weak_options->weakzone_2plates_ridge_depth;
  ridge_width = weak_options->weakzone_2plates_ridge_width;
  ridge_smoothwidth = weak_options->weakzone_2plates_ridge_smoothwidth;
  ridge_weak_factor = weak_options->weakzone_2plates_ridge_weak_factor;

  /* check parameters */
  RHEA_ASSERT (0.0 < subdu_lon);
  RHEA_ASSERT (0.0 < subdu_depth && 0.0 < subdu_width);
  RHEA_ASSERT (0.0 < subdu_thickness);
  RHEA_ASSERT (subdu_thickness_const <= subdu_thickness);
  RHEA_ASSERT (0.0 < subdu_weak_factor && subdu_weak_factor <= 1.0);
  RHEA_ASSERT (0.0 < ridge_depth && 0.0 < ridge_width);
  RHEA_ASSERT (0.0 <= ridge_smoothwidth);
  RHEA_ASSERT (0.0 < ridge_weak_factor && ridge_weak_factor <= 1.0);

  /*
   * set subduction weak zone between plates
   */

  start_node = geo->start_node;
  start_val = geo->start_val;
  start_deriv = geo->start_deriv;
  end_node = geo->end_node;
  end_val = geo->end_val;
  poly2_coeff = geo->poly2_coeff;

  /* only consider point in a rectangle containing the weak zone */
  weakzone_thickness = subdu_thickness / SUBD_MANTLE_DEPTH;
  courtesy_width = 0.5 * weakzone_thickness;
  if (   (start_node - 0.5 * weakzone_thickness - courtesy_width) <= lon
      && lon <= (end_node + 0.5 * weakzone_thickness + courtesy_width)
      && (end_val - 0.5 * weakzone_thickness - courtesy_width) <= r ) {

    /* compute distance to subduction weak zone */
    dist = subd_weakzone_subduct_dist_2plates_poly2 (r, lon, poly2_coeff,
                                                      start_node, start_val, start_deriv,
                                                      end_node, end_val);
    if (r > 0.95) {
      /*from r=0.9 to r=1, shrink the width of weakzone by a factor of 2*/
      s = - 10.0 * r + 10.5;
      a = subdu_thickness * s;
      b = subdu_thickness_const * s;
      weak = subd_weakzone_factor_fn (dist, a,
                                       b, subdu_weak_factor);
    }
    else {
    /* set weak zone factor */
    weak = subd_weakzone_factor_fn (dist, subdu_thickness,
                                     subdu_thickness_const, subdu_weak_factor);
    }
  }

  /*
   * set ridge weak zone in the left corner of domain, so that subducting
   * plate can "move"
   */

  /* set bottom left corner of weak zone */
  end_node = y_min + ridge_width / SUBD_MANTLE_DEPTH;
  end_val = SUBD_SHELL_RADIUS_TOP - ridge_depth / SUBD_MANTLE_DEPTH;

  /* only consider points close to weak zone */
  courtesy_width = 2.0 * ridge_smoothwidth / SUBD_MANTLE_DEPTH;
  if (lon <= (end_node + courtesy_width) && (end_val - courtesy_width) <= r) {
    /* compute distance to ridge weak zone */
    dist = subd_weakzone_ridge_dist_2plates_poly2 (r, lon, end_node, end_val);

    /* set weak zone factor */
    weak = subd_weakzone_factor_fn (dist, ridge_smoothwidth, 0.0,
                                     ridge_weak_factor);
  }

  /*
   * return weak zone factor
   */
  return weak;
}

//static double
double
subd_weakzone_node (const double x, const double y, const double z,
                     subd_options_t *subd_options)
{
  double              lon;
  double              r;

  /* compute radius and longitude */
  lon = y;
  r = z;

  /* compute weak zone factor */
  return subd_weakzone_brick_2plates_poly2 (r, lon, subd_options);
}

/* compute weak zone factor of an element */
void
subd_weakzone_elem (double *_sc_restrict weak_elem,
                     const double *x,
                     const double *y,
                     const double *z,
                     const int n_nodes_per_el,
                     subd_options_t *subd_options)
{
  int            nodeid;
  double         z_mid = subd_options->visc_options->z_lith;

  /* compute weak zone factor for each node */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) { /* loop over all
                                                         * nodes */
    weak_elem[nodeid] = subd_weakzone_node (x[nodeid], y[nodeid], z[nodeid],
                                             subd_options);

    /*temporary, to make sure that the viscosity along the weak zone does not have a sharp change
     * at the lith-athen boundary*/
    if (z[nodeid] < z_mid)  {
      weak_elem[nodeid] = 1.0;
    }
  }
}

void
subd_poly2_weakzone_compute (ymir_dvec_t *weakzone, void *data)
{
  const char         *this_fn_name = "subd_poly2_weakzone_compute";
  subd_options_t    *subd_options = data;
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (weakzone);
  const ymir_locidx_t n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  sc_dmatrix_t       *weak_el_mat;
  double             *weak_el_data;
  double             *x, *y, *z, *tmp_el;
  ymir_locidx_t      elid;

  double              start_node, start_val, start_deriv, end_node, end_val;
  double              subdu_lon;
  double              subdu_dip_angle;
  double              subdu_depth, subdu_width;
  double              *poly2_coeff;
  subd_2plates_poly2_geo_coeff_t geo;
  subd_weak_options_t *weak_options = subd_options->weak_options;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* set parameters according to weakzone options */
  subdu_lon = weak_options->weakzone_2plates_subdu_longitude;
  subdu_dip_angle = weak_options->weakzone_2plates_subdu_dip_angle;
  subdu_depth = weak_options->weakzone_2plates_subdu_depth;
  subdu_width = weak_options->weakzone_2plates_subdu_width;

  /* check parameters */
  RHEA_ASSERT (0.0 < subdu_lon);
  RHEA_ASSERT (0.0 < subdu_dip_angle);
  RHEA_ASSERT (0.0 < subdu_depth && 0.0 < subdu_width);

  /* set points for polynomial interpolation */
  start_node = subdu_lon;
  start_val = SUBD_SHELL_RADIUS_TOP;
  start_deriv = tan (-subdu_dip_angle / 180.0 * M_PI);
  end_node = start_node + subdu_width / SUBD_MANTLE_DEPTH;
  end_val = start_val - subdu_depth / SUBD_MANTLE_DEPTH;
  /* compute interpolating quadratic polynomial */
  poly2_coeff = subd_compute_poly2_interpolation (start_node, start_val,
                                                   start_deriv,
                                                   end_node, end_val);

  RHEA_GLOBAL_PRODUCTIONF("weakzone: brick_2plates_poly2 downgoing subd geometry\nstart_node: (%f %f), start_deriv: %f, end_node: (%f %f)\npoly2_coeff: (%f %f %f)\n", start_node, start_val, start_deriv, end_node, end_val, poly2_coeff[2], poly2_coeff[1], poly2_coeff[0]);

  geo.start_node = start_node;
  geo.start_val = start_val;
  geo.start_deriv = start_deriv;
  geo.end_node = end_node;
  geo.end_val = end_val;
  geo.poly2_coeff = poly2_coeff;
  weak_options->weak_2plates_geo_coeff = &geo;

  /* check input */
  RHEA_ASSERT (weakzone->ndfields == 1);

  /* create work variables */
  weak_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  weak_el_data = weak_el_mat->e[0];
  x = RHEA_ALLOC (double, n_nodes_per_el);
  y = RHEA_ALLOC (double, n_nodes_per_el);
  z = RHEA_ALLOC (double, n_nodes_per_el);
  tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    /* get coordinates of this element at Gauss nodes */
    ymir_mesh_get_elem_coord_gauss (x, y, z, elid, mesh, tmp_el);

    /* compute weak zone factor */
    subd_weakzone_elem (weak_el_data, x, y, z,
                         n_nodes_per_el, subd_options);

    /* set weak zone of this element */
    rhea_viscosity_set_elem_gauss (weakzone, weak_el_mat, elid);
  }

  /* destroy */
  sc_dmatrix_destroy (weak_el_mat);
  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);
  RHEA_FREE (poly2_coeff);  //TO FIND OUT: if I use RHEA_FREE, there is memory leak stderr

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

void
subd_compute_weakzone (rhea_stokes_problem_t *stokes_problem,
                        subd_options_t *subd_options)
{
  const char      *this_fn_name = "subd_compute_weakzone";

  if (subd_options->weak_options->type != SUBD_WEAK_NONE) {
    RHEA_GLOBAL_PRODUCTIONF ("%s: compute weakzone\n", this_fn_name);
    rhea_stokes_problem_set_weakzone_compute_fn (
           stokes_problem, subd_poly2_weakzone_compute, subd_options);
  }
}


/**************************************
 * Viscosity Computation
 *************************************/
//TODO: check whether this subroutine works
void
subd_layers_viscosity_smooth_elem (double *_sc_restrict visc_elem,
                            const double *_sc_restrict x,
                            const double *_sc_restrict y,
                            const double *_sc_restrict z,
                            const int n_nodes_per_el,
                            subd_options_t *subd_options)
{
  int                 nodeid;
  double              z_mid = subd_options->visc_options->z_lith;
  double              visc_lith = subd_options->visc_options->visc_lith;
  double              visc_asthen = subd_options->visc_options->visc_asthen;
  double              visc_smooth;
  double              factor = 1.0;
  subd_topo_profile_t *topo = subd_options->surf_options->topo_profile;
  double *tX = topo->tX;
  double *tY = topo->tY;
  double *tZ = topo->tZ;
  int     m, nsurf = topo->nsurf;

  /* compute viscosity in this element */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    visc_smooth = 10.0 * (z[nodeid]/factor - 0.45) * (visc_lith - visc_asthen) + visc_asthen;
    if (z[nodeid]/factor >= z_mid)  {
      visc_elem[nodeid] = (z[nodeid]/factor-z_mid)>0.05? visc_lith: visc_smooth;
    }
    else {
      visc_elem[nodeid] = (z[nodeid]/factor-z_mid)<-0.05? visc_asthen: visc_smooth;
    }

    /* check viscosity for `nan`, `inf`, and positivity */
    RHEA_ASSERT (isfinite (visc_elem[nodeid]));
    RHEA_ASSERT (0.0 < visc_elem[nodeid]);
  }
}

static void
subd_layers_viscosity_elem (double *_sc_restrict visc_elem,
                            const double *_sc_restrict x,
                            const double *_sc_restrict y,
                            const double *_sc_restrict z,
                            const int n_nodes_per_el,
                            subd_options_t *subd_options)
{
  int                 nodeid;
  double              z_mid = subd_options->visc_options->z_lith;
  double              visc_lith = subd_options->visc_options->visc_lith;
  double              visc_asthen = subd_options->visc_options->visc_asthen;
  int                 m;

  /* compute viscosity in this element */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    if (z[nodeid] >= z_mid)  {
      visc_elem[nodeid] = visc_lith;
    }
    else {
      visc_elem[nodeid] = visc_asthen;
    }

    /* check viscosity for `nan`, `inf`, and positivity */
    RHEA_ASSERT (isfinite (visc_elem[nodeid]));
    RHEA_ASSERT (0.0 < visc_elem[nodeid]);
  }
}


static void
subd_layers_coupling_viscosity_elem (double *_sc_restrict visc_elem,
                                      const double *_sc_restrict x,
                                      const double *_sc_restrict y,
                                      const double *_sc_restrict z,
                                      const int n_nodes_per_el,
                                      subd_options_t *subd_options)
{
  int                 nodeid;
  double              z_mid = subd_options->visc_options->z_lith;
  double              z_slab = subd_options->visc_options->z_slab;
  double              slab_width = subd_options->visc_options->slab_width;
  double              y_center = 0.5;

  /* compute viscosity in this element */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    if (z[nodeid] >= z_mid)  {
      visc_elem[nodeid] = subd_options->visc_options->visc_lith;
    }
    else if (z[nodeid] >= z_slab && fabs(y[nodeid] - y_center) <= 0.5 * slab_width) {
      visc_elem[nodeid] = subd_options->visc_options->visc_lith;
    }
    else {
      visc_elem[nodeid] = subd_options->visc_options->visc_asthen;
    }

    /* check viscosity for `nan`, `inf`, and positivity */
    RHEA_ASSERT (isfinite (visc_elem[nodeid]));
    RHEA_ASSERT (0.0 < visc_elem[nodeid]);
  }
}

void
subd_custom_visc  (ymir_vec_t *viscosity,
                         ymir_vec_t *rank1_tensor_scal,
                         ymir_vec_t *bounds_marker,
                         ymir_vec_t *yielding_marker,
                         ymir_vec_t *temperature,
                         ymir_vec_t *weakzone,
                         ymir_vec_t *velocity,
                         void *data)
{
  subd_options_t  *subd_options = data;
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (viscosity);
  const ymir_locidx_t  n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);
  mangll_t           *mangll = mesh->ma;
  const int           N = ymir_n (mangll->N);

  sc_dmatrix_t       *visc_el_mat;
  double             *x, *y, *z, *tmp_el,*visc_el_data, *weak_el_data;
  ymir_locidx_t       elid;
  char                *this_fn_name = "subd_custom_visc";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* create work variables */
  visc_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  visc_el_data = visc_el_mat->e[0];
  x = RHEA_ALLOC (double, n_nodes_per_el);
  y = RHEA_ALLOC (double, n_nodes_per_el);
  z = RHEA_ALLOC (double, n_nodes_per_el);
  tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) {
    /* get coordinates of this element at Gauss nodes */
    ymir_mesh_get_elem_coord_gauss (x, y, z, elid, mesh, tmp_el);
    /* compute user defined weak zone viscosity*/
    switch (subd_options->visc_options->custom_type) {
      case SUBD_VISC_CUSTOM_LAYERS:
        subd_layers_viscosity_elem (visc_el_data, x, y, z, n_nodes_per_el,
                                     subd_options);
      break;

      case SUBD_VISC_CUSTOM_LAYERS_COUPLING:
        subd_layers_coupling_viscosity_elem (visc_el_data, x, y, z, n_nodes_per_el,
                                             subd_options);
      break;

      default:
      RHEA_ABORT_NOT_REACHED ();
    }

    /* set viscosity of this element */
    rhea_viscosity_set_elem_gauss (viscosity, visc_el_mat, elid);
  }

  RHEA_GLOBAL_PRODUCTIONF ("hello1 %s\n", this_fn_name);

  /* destroy */
  sc_dmatrix_destroy (visc_el_mat);
  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

void
subd_viscosity_set_function (rhea_stokes_problem_t *stokes_problem,
                            subd_options_t *subd_options)
{
  const char      *this_fn_name = "subd_viscosity_set_function";
  if (subd_options->visc_options->type != SUBD_VISC_RHEA) {
    RHEA_GLOBAL_PRODUCTIONF ("%s: reset viscosity function\n", this_fn_name);
    rhea_stokes_problem_set_viscosity_compute_fn (
                stokes_problem,
                subd_custom_visc,
                subd_options);
  }
}

/**************************************************************
 * Right-hand-side
***************************************************************/

/* data for callback function to compute the velocity right-hand side */
typedef struct subd_density_rhs_fn_data
{
    ymir_vec_t          *temperature;
    subd_options_t      *subd_options;
}
subd_density_rhs_fn_data_t;

static void
subd_rhs_full_density_fn (double *rhs, double x, double y, double z,
                        ymir_locidx_t nodeid, void *data)
{
  subd_density_rhs_fn_data_t *d = data;
  subd_options_t             *subd_opt = d->subd_options;
  rhea_domain_options_t *domain_options = subd_opt->domain_options;
  const double       scaling = subd_opt->para_options->rayleigh;
  const double       therm_expa = subd_opt->para_options->ref_therm_expa;
  const double       temp_diff = subd_opt->para_options->ref_temp_diff;
  double             shift= -1.0 / therm_expa / temp_diff;
  double             temp = *ymir_cvec_index (d->temperature, nodeid, 0);

  switch (domain_options->shape) {
  case RHEA_DOMAIN_CUBE:
  case RHEA_DOMAIN_BOX:
    /*   f(x) = Ra * (T - 1/(alpha*DeltaT) ) */
    rhs[0] = 0.0;
    rhs[1] = 0.0;
    rhs[2] = scaling * (temp + shift);
    break;
  case RHEA_DOMAIN_SHELL:
  case RHEA_DOMAIN_CUBE_SPHERICAL:
  case RHEA_DOMAIN_BOX_SPHERICAL:
    {
      const double        radius = rhea_domain_compute_radius (x, y, z,
                                                               domain_options);

      rhs[0] = scaling * x / radius * (temp + shift);
      rhs[1] = scaling * y / radius * (temp + shift);
      rhs[2] = scaling * z / radius * (temp + shift);
    }
    break;
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }
}

static void
subd_rhs_full_density (ymir_vec_t *rhs_vel,
                              ymir_vec_t *temp, void *data)
{
  subd_options_t      *subd_opt = data;
  subd_density_rhs_fn_data_t d;
  d.temperature = temp;
  d.subd_options = subd_opt;
  ymir_cvec_set_function (rhs_vel, subd_rhs_full_density_fn, &d);
}

static void
subd_rhs_adjoint_none (ymir_vec_t *rhs_vel,
                              ymir_vec_t *temp, void *data)
{
  ymir_cvec_set_value (rhs_vel, .0);
}

void
subd_compute_rhs_vel (rhea_stokes_problem_t *stokes_problem,
                      void *data)
{
  subd_options_t *subd_opt = (subd_options_t *) data;
  ymir_vec_t *temp = rhea_stokes_problem_get_temperature (stokes_problem);
  ymir_stokes_op_t *stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
  const char  *this_fn_name = "subd_compute_rhs_vel";

  if (subd_opt->rhs_type != SUBD_RHS_BUOY)  {
    RHEA_GLOBAL_PRODUCTIONF ("%s: recompute rhs_vel\n", this_fn_name);
    switch (subd_opt->rhs_type) {
    case SUBD_RHS_DENSITY:
      rhea_stokes_problem_set_rhs_vel_compute_fn (stokes_problem,
                  subd_rhs_full_density, subd_opt);
      break;

    case SUBD_RHS_ADJOINT_NONE:
      rhea_stokes_problem_set_rhs_vel_compute_fn (stokes_problem,
                  subd_rhs_adjoint_none, NULL);
      break;

    case SUBD_RHS_TEST_MANUFACTURED:
      rhea_stokes_problem_set_rhs_vel_compute_fn (stokes_problem,
                  subd_test_manufactured_rhs, subd_opt);
      break;

    default:
      RHEA_ABORT_NOT_REACHED ();
    }
  }
}


/**************************************************************
 * Boundary conditions
***************************************************************/

static ymir_dir_code_t
subd_set_vel_dir_all_2D (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
  if (face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE3 ||
      face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE4 ||
      face == RHEA_DOMAIN_BOUNDARY_FACE_BASE  ||
      face == RHEA_DOMAIN_BOUNDARY_FACE_TOP) {
      return YMIR_VEL_DIRICHLET_ALL;
  }
  else
    return YMIR_VEL_DIRICHLET_NORM;
}


/* free-surface */
static ymir_dir_code_t
subd_set_vel_dir_freesurface (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
  if (face == RHEA_DOMAIN_BOUNDARY_FACE_TOP)
    return YMIR_VEL_DIRICHLET_NONE;
  else
    return YMIR_VEL_DIRICHLET_NORM;
}

static ymir_dir_code_t
subd_set_vel_dir_mixsurface (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
  if (face == RHEA_DOMAIN_BOUNDARY_FACE_TOP) {
    if (Y >= 0.25 && Y <= 1.25)
      return YMIR_VEL_DIRICHLET_NONE;
    else
      return YMIR_VEL_DIRICHLET_NORM;
  }
  else
    return YMIR_VEL_DIRICHLET_NORM;
}

/* Dirichlet all on one side of the domain: SIDE3 (y=0).*/
static ymir_dir_code_t
subd_set_vel_dir_sidewall_ymin_dirall (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
  if (face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE3) {
    return YMIR_VEL_DIRICHLET_ALL;
  }
  else {
    return YMIR_VEL_DIRICHLET_NORM;
  }
}

/* Dirichlet all on both left and right sides of the domain: SIDE3 (y=0), and SIDE4 (y=ymax).*/
static ymir_dir_code_t
subd_set_vel_dir_sidewalls_y_dirall (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
  if (face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE3 ||
      face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE4) {
    return YMIR_VEL_DIRICHLET_ALL;
  }
  else {
    return YMIR_VEL_DIRICHLET_NORM;
  }
}

/* Dirichlet all on one side of the domain: SIDE3 (y=0), and Neumann at the base*/
static ymir_dir_code_t
subd_set_vel_dir_sidewall_ymin_dirall_freebase (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
  if (face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE3) {
    return YMIR_VEL_DIRICHLET_ALL;
  }
  else if (face == RHEA_DOMAIN_BOUNDARY_FACE_BASE) {
     return YMIR_VEL_DIRICHLET_NONE;
  }
  else {
    return YMIR_VEL_DIRICHLET_NORM;
  }
}


void
subd_set_velocity_dirichlet_bc (rhea_domain_options_t *domain_options,
                                subd_options_t *subd_options)
{
  /*set velocity_dirichlet_bc function*/
  if (domain_options->velocity_bc_type == RHEA_DOMAIN_VELOCITY_BC_USER) {
    int       nonzero_dir = (int) subd_options->velbc_options->velbc_nonzero_dir;

    switch (subd_options->velbc_options->velbc) {
      case SUBD_FREESURFACE:
         rhea_domain_set_user_velocity_dirichlet_bc (
             subd_set_vel_dir_freesurface, NULL /* no data necessary */,
             0 );
      break;

      case SUBD_MIXSURFACE:
         rhea_domain_set_user_velocity_dirichlet_bc (
             subd_set_vel_dir_mixsurface, NULL /* no data necessary */,
             0 );
      break;

      case SUBD_SIDEWALL_YMIN_DIRALL:
         rhea_domain_set_user_velocity_dirichlet_bc (
             subd_set_vel_dir_sidewall_ymin_dirall, NULL /* no data necessary */,
             0 );

      case SUBD_SIDEWALLS_Y_DIRALL:
         rhea_domain_set_user_velocity_dirichlet_bc (
             subd_set_vel_dir_sidewalls_y_dirall, NULL /* no data necessary */,
             0 );

      case SUBD_SIDEWALL_YMIN_DIRALL_FREEBASE:
         rhea_domain_set_user_velocity_dirichlet_bc (
             subd_set_vel_dir_sidewall_ymin_dirall_freebase, NULL /* no data necessary */,
             0 );

      break;

      default: /* BC not set */
         RHEA_ABORT_NOT_REACHED ();
    }
  }
}

/**************************************************************
* Non-zero Dirichlet
***************************************************************/

/* In-out flow sine velocity on one side of the domain. */
static void
subd_set_rhs_vel_nonzero_dir_inoutflow_sin (
    double *vel, double x, double y, double z,
    ymir_locidx_t nid, void *data)
{
  subd_options_t  *subd_options = data;
  const double     flow_scale = subd_options->velbc_options->flow_scale;

  if (y < SC_1000_EPS) {
    vel[0] = 0.0;
    vel[1] = -flow_scale * sin (2.0 * M_PI * z);
    vel[2] = 0.0;
  }
  else {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }
}

/* In-out flow tanh velocity on one side of the domain. 2 layers */
static void
subd_set_rhs_vel_nonzero_dir_inoutflow_tanh_2layer (double *vel, double x, double y, double z,
                                              ymir_locidx_t nid, void *data)
{
  subd_options_t  *subd_options = data;
  subd_velbc_options_t *velbc_options = subd_options->velbc_options;
  const double        z_max = subd_options->domain_options->z_max;
  const double        flow_scale = velbc_options->flow_scale;
  const double        zM = velbc_options->vel_dir_bc_middle;
  const double        a = z_max-zM, b = z_max-a;
  const double        shape = 2.0 * M_PI, scaling = 0.5*flow_scale, shift = 0.5*(b-a)*flow_scale;
  double txM = shape * (z-zM);

  if (y < SC_1000_EPS) {
    vel[0] = 0.0;
    vel[2] = 0.0;
    vel[1] = shift + scaling *
             ( (exp (txM) - exp (-txM)) /
             (exp (txM) + exp (-txM)) );
  }
  else {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }
}

/* In-out flow tanh velocity on one side of the domain.*/
static void
subd_set_rhs_vel_nonzero_dir_inoutflow_tanh_3layer (double *vel, double x, double y, double z,
                                              ymir_locidx_t nid, void *data)
{
  subd_options_t  *subd_options = data;
  subd_velbc_options_t *velbc_options = subd_options->velbc_options;
  const double        z_max = subd_options->domain_options->z_max;
  const double        flow_scale = velbc_options->flow_scale;
  const double        zU = velbc_options->vel_dir_bc_upper;
  const double        zL = velbc_options->vel_dir_bc_lower;
  const double        a = zU-zL, b = z_max-a, c = 0.5*(zU+zL);
  const double        shape = 2.0 * M_PI, scaling = 0.5*(b+a)*flow_scale, shift = 0.5*(b-a)*flow_scale;
  double txL = shape*(z-zL),txU = shape*(z-zU);

  if (y < SC_1000_EPS) {
    vel[0] = 0.0;
    vel[2] = 0.0;
    if (z<=c)
      vel[1] = shift + scaling *
             ( (exp (txL) - exp (-txL)) /
             (exp (txL) + exp (-txL)) );
    else
      vel[1] = shift - scaling *
             ( (exp (txU) - exp (-txU)) /
             (exp (txU) + exp (-txU)) );
  }
  else {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }
}

/* In-out flow tanh velocity on both left and right sides of the domain.*/
static void
subd_set_rhs_vel_nonzero_dir_inoutflow_double_tanh_3layer (
    double *vel, double x, double y, double z,
    ymir_locidx_t nid, void *data)
{
  subd_options_t  *subd_options = data;
  subd_velbc_options_t *velbc_options = subd_options->velbc_options;
  const double        z_max = subd_options->domain_options->z_max;
  const double        y_max = subd_options->domain_options->y_max;
  const double        flow_scale = velbc_options->flow_scale;
  const double        zU = velbc_options->vel_dir_bc_upper;
  const double        zL = velbc_options->vel_dir_bc_lower;
  const double        a = zU-zL, b = z_max-a, c = 0.5*(zU+zL);
  const double        shape = 2.0 * M_PI, scaling = 0.5*(b+a)*flow_scale, shift = 0.5*(b-a)*flow_scale;
  double txL = shape*(z-zL),txU = shape*(z-zU);

  if (y < SC_1000_EPS) {
    vel[0] = 0.0;
    vel[2] = 0.0;
    if (z<=c)
      vel[1] = shift + scaling *
               ( (exp (txL) - exp (-txL)) /
               (exp (txL) + exp (-txL)) );
    else
      vel[1] = shift - scaling *
               ( (exp (txU) - exp (-txU)) /
               (exp (txU) + exp (-txU)) );
  }
  else if ((y_max - y) < SC_1000_EPS) {
    vel[0] = 0.0;
    vel[2] = 0.0;
    if (z<=c)
      vel[1] = -shift - scaling *
             ( (exp (txL) - exp (-txL)) /
             (exp (txL) + exp (-txL)) );
    else
      vel[1] = -shift + scaling *
             ( (exp (txU) - exp (-txU)) /
             (exp (txU) + exp (-txU)) );
  }
  else {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }
}

static void
subd_velbc_nonzero_dirichlet (ymir_vec_t * rhs_vel_nonzero_dirichlet,
                                     void * data)
{
  subd_options_t    *subd_options = data;

  switch (subd_options->velbc_options->velbc_nonzero_dir) {
    case SUBD_VELBC_DIR_INOUTFLOW_SIN:
      ymir_cvec_set_function (rhs_vel_nonzero_dirichlet,
                              subd_set_rhs_vel_nonzero_dir_inoutflow_sin,
                              subd_options);
      break;

    case SUBD_VELBC_DIR_INOUTFLOW_TANH_TWOLAYER:
       ymir_cvec_set_function (rhs_vel_nonzero_dirichlet,
                              subd_set_rhs_vel_nonzero_dir_inoutflow_tanh_2layer,
                              subd_options);
      break;

    case SUBD_VELBC_DIR_INOUTFLOW_TANH_THREELAYER:
      ymir_cvec_set_function (rhs_vel_nonzero_dirichlet,
                              subd_set_rhs_vel_nonzero_dir_inoutflow_tanh_3layer,
                              subd_options);
      break;

    case SUBD_VELBC_DIR_INOUTFLOW_DOUBLE_TANH_THREELAYER:
      ymir_cvec_set_function (rhs_vel_nonzero_dirichlet,
                              subd_set_rhs_vel_nonzero_dir_inoutflow_double_tanh_3layer,
                              subd_options);
      break;

   default: /* BC not set */
      RHEA_ABORT_NOT_REACHED ();
  }
}

void
subd_compute_rhs_velbc_dirichlet (rhea_stokes_problem_t *stokes_problem,
                              subd_options_t *subd_options)
{
  const char    *this_fn_name = "subd_compute_rhs_velbc_dirichlet";

  if (subd_options->velbc_options->velbc_nonzero_dir != SUBD_VELBC_DIR_ZERO) {
    RHEA_GLOBAL_PRODUCTIONF ("%s: compute nonzero dirichelt vel_bc\n", this_fn_name);
    rhea_stokes_problem_set_rhs_vel_nonzero_dir_compute_fn (
          stokes_problem, subd_velbc_nonzero_dirichlet, subd_options);
  }
  if (subd_options->test_options->test_manufactured != SUBD_TEST_MANUFACTURED_NONE) {
    RHEA_GLOBAL_PRODUCTIONF ("%s: compute nonzero dirichelt vel_bc for manufactured test\n", this_fn_name);
      rhea_stokes_problem_set_rhs_vel_nonzero_dir_compute_fn (
          stokes_problem, subd_test_manufactured_velbc_dir, subd_options);
  }
}

/**************************************************************
* Neumann boundary condition
***************************************************************/

static void
subd_nz_neumann_sine_set_fn (double *neu_surf, double x, double y, double z,
                              double nx, double ny, double nz,
                              ymir_topidx_t face, ymir_locidx_t nid, void *data)
{
  neu_surf[0] = sin(M_PI * x);
  neu_surf[2] = cos(M_PI * x);
  neu_surf[1] = .0;
}

static void
subd_velbc_nonezero_neumann (ymir_vec_t **vel_neumann, void *data)
{
  const char         *this_fn_name = "subd_velbc_nonezero_neumann";
  subd_options_t    *subd_options = data;
  int                 fm;
  char              *file_path = subd_options->velbc_options->txt_read_path_nonzero_neu;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  for (fm = 0; fm < 6; fm++)  {
    ymir_vec_set_value (vel_neumann[fm], .0);
  }

  fm = (int) RHEA_DOMAIN_BOUNDARY_FACE_TOP;
  switch (subd_options->velbc_options->velbc_nonzero_neu) {
    case SUBD_VELBC_NEU_SINE:
      ymir_face_cvec_set_function (vel_neumann[fm],
                              subd_nz_neumann_sine_set_fn,
                              NULL);
    break;

    case SUBD_VELBC_NEU_READ:
      if (file_path != NULL) {
        /*read in surfvelo*/
        subd_facevec_read (vel_neumann[fm], file_path);
        /*In adjoint B.C., neumann + surfvelo = 0*/
        ymir_vec_scale (-1.0, vel_neumann[fm]);
      }
    break;

    default:
      RHEA_ABORT_NOT_REACHED ();
  }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

}

void
subd_compute_rhs_velbc_neumann (rhea_stokes_problem_t *stokes_problem,
                              subd_options_t *subd_options)
{
  const char      *this_fn_name = "subd_compute_rhs_velbc_neumann";

  if (subd_options->velbc_options->velbc_nonzero_neu != SUBD_VELBC_NEU_ZERO) {
    RHEA_GLOBAL_PRODUCTIONF ("%s: compute nonzero neumann vel_bc\n", this_fn_name);
      rhea_stokes_problem_set_rhs_vel_nonzero_neu_compute_fn (
          stokes_problem, subd_velbc_nonezero_neumann, subd_options);
  }
}

