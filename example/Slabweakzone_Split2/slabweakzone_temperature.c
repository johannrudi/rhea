/**************************************
 * Temperature Computation
 *************************************/
double
slabs_temperature_brick_2plates_poly2 (double r, double lon,
                                       slabs_options_t * slabs_options)
{
  double              trench_lon;
  double              subd_depth, subd_width;
  double              subd_edge_width, subd_edge_smoothwidth;
  double              subd_plate_vel, subd_plate_init_age;
  double              over_plate_age;

  double              y_min = slabs_options->slabs_domain_options->y_min;
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

  slabs_temp_options_t * temp_options = slabs_options->slabs_temp_options;
  slabs_2plates_poly2_geo_coeff_t *geo = temp_options->temp_2plates_geo_coeff;

  /* set parameters according to temperature options */
  trench_lon = temp_options->temp_2plates_trench_longitude;
  subd_depth = temp_options->temp_2plates_subd_depth;
  subd_width = temp_options->temp_2plates_subd_width;
  subd_edge_width = temp_options->temp_2plates_subd_edge_width;
  subd_edge_smoothwidth = temp_options->temp_2plates_subd_edge_smoothwidth;
  subd_plate_vel = temp_options->temp_2plates_subd_plate_velocity;
  subd_plate_init_age = temp_options->temp_2plates_subd_plate_initial_age;
  over_plate_age = temp_options->temp_2plates_over_plate_age;

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
  closest_pt = slabs_compute_closest_pt_on_poly2 (lon, r,
                                                  poly2_coeff, start_node,
                                                  start_val, start_deriv,
                                                  end_node, end_val,
                                                  &orientation_wrt_curve);


  /* compute "fake" (TODO) distance */
  dist = SLABS_MANTLE_DEPTH * sqrt ( (lon - closest_pt[0]) * (lon - closest_pt[0])
                                  + (r - closest_pt[1]) * (r - closest_pt[1]) );

  /* For computing the temperature, use temperature profile of the plates
   * from halfspace cooling model, erf(z/(2*sqrt(t*kappa))), as base
   * temperature field. */

  /* compute temperature of plates at surface */
  if (orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM ||
      orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_BACK) { /* if subd plate */
    /* calculate age of subducting plate */
    ridge_dist = fabs (lon - y_min) * SLABS_MANTLE_DEPTH;
    vel = subd_plate_vel / SLABS_SEC_PER_YEAR;
    plate_time = ridge_dist / vel + subd_plate_init_age * SLABS_SEC_PER_YEAR;

    /* avoid division by zero */
    plate_time = SC_MAX (1.0e-3, plate_time);
  }
  else  { /* if overriding plate */
    /* calculate age of overriding plate */
    plate_time = over_plate_age * SLABS_SEC_PER_YEAR;
  }
  depth = (SLABS_SHELL_RADIUS_TOP - r) * SLABS_MANTLE_DEPTH;
  temp = erf ( depth / (2.0 * sqrt (plate_time * SLABS_THERM_DIFFUS)) );

  /* compute temperature of the plate subducted inside of mantle */
  if (   orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM
      || orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_BACK
      || orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_LEFT ) {
    /* calculate total lenght of subducted plate */
    subd_length = SLABS_MANTLE_DEPTH * sqrt (
                    (start_node - end_node) * (start_node - end_node)
                    + (start_val - end_val) * (start_val - end_val) );

    /* calculate age of plate */
    if (lon <= trench_lon) {
      ridge_dist = fabs (lon - y_min) * SLABS_MANTLE_DEPTH;
    }
    else {
      ridge_dist = fabs (trench_lon - y_min) * SLABS_MANTLE_DEPTH
                   + subd_length * depth / subd_depth;
    }
    vel = subd_plate_vel / SLABS_SEC_PER_YEAR;
    plate_time = ridge_dist / vel + subd_plate_init_age * SLABS_SEC_PER_YEAR;

    /* avoid division by zero */
    plate_time = SC_MAX (1.0e-3, plate_time);

    /* compute temperature of subducting plate */
    temp_subd = erf ( dist / (2.0 * sqrt (plate_time * SLABS_THERM_DIFFUS)) );
  }

  /* smooth the subducting plate's top edge */
  if (   orientation_wrt_curve == SLABS_CURVE_ORIENT_TOP
      || orientation_wrt_curve == SLABS_CURVE_ORIENT_TOP_RIGHT ) {
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
  if (   orientation_wrt_curve == SLABS_CURVE_ORIENT_TOP_RIGHT
      || orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_LEFT ) {
    /* compute distance between tip and closest point on curve */
    dist = SLABS_MANTLE_DEPTH * sqrt (
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
slabs_temperature_set_fn (double *temp, double x, double y, double z,
                          ymir_locidx_t nid, void *data)
{
  slabs_options_t    *slabs_options = (slabs_options_t *) data;
  double              lon, r;

  /* compute radius and longitude */
  lon = y;
  r = z;

  *temp = slabs_temperature_brick_2plates_poly2 (r, lon, slabs_options);

  /* check temperature for `nan` and `inf` */
  RHEA_ASSERT (isfinite (*temp));

  /* bound temperature to valid interval */
  *temp = SC_MIN (1.0, *temp);
  *temp = SC_MAX (0.0, *temp);
}

static void
slabs_temperature_compute (ymir_vec_t *temperature,
                           slabs_options_t *slabs_options)
{
  const char         *this_fn_name = "slabs_temperature_compute";
  double              trench_lon;
  double              dip_angle;
  double              subd_depth, subd_width;
  double              start_node, start_val, start_deriv, end_node, end_val;
  double              *poly2_coeff;
  slabs_2plates_poly2_geo_coeff_t geo;
  slabs_temp_options_t * temp_options = slabs_options->slabs_temp_options;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* set parameters according to physics options */
  trench_lon = temp_options->temp_2plates_trench_longitude;
  dip_angle = temp_options->temp_2plates_dip_angle;
  subd_depth = temp_options->temp_2plates_subd_depth;
  subd_width = temp_options->temp_2plates_subd_width;

  /* check parameters */
  RHEA_ASSERT (0.0 < trench_lon);
  RHEA_ASSERT (0.0 < dip_angle && dip_angle < 90.0);
  RHEA_ASSERT (0.0 < subd_depth && 0.0 < subd_width);

  /* set points for polynomial interpolation */
  start_node = trench_lon;
  start_val = SLABS_SHELL_RADIUS_TOP;
  start_deriv = tan (-dip_angle / 180.0 * M_PI);
  end_node = start_node + subd_width / SLABS_MANTLE_DEPTH;
  end_val = start_val - subd_depth / SLABS_MANTLE_DEPTH;
  /* compute interpolating quadratic polynomial */
  poly2_coeff = slabs_compute_poly2_interpolation (start_node, start_val,
                                                   start_deriv,
                                                   end_node, end_val);

  RHEA_GLOBAL_PRODUCTIONF("temperature: brick_2plates_poly2 downgoing slabs geometry\nstart_node: (%f %f), start_deriv: %f, end_node: (%f %f)\npoly2_coeff: (%f %f %f)\n", start_node, start_val, start_deriv, end_node, end_val, poly2_coeff[2], poly2_coeff[1], poly2_coeff[0]);

  geo.start_node = start_node;
  geo.start_val = start_val;
  geo.start_deriv = start_deriv;
  geo.end_node = end_node;
  geo.end_val = end_val;
  geo.poly2_coeff = poly2_coeff;
  temp_options->temp_2plates_geo_coeff = &geo;

  ymir_cvec_set_function (temperature,
                          slabs_temperature_set_fn,
                          slabs_options);
  RHEA_FREE (poly2_coeff);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

void
drag_temperature_set_fn (double *temp, double x, double y, double z,
                          ymir_locidx_t nid, void *data)
{
  slabs_options_t    *slabs_options = (slabs_options_t *) data;
  double              subdu_lon, subdu_dip_angle, subdu_depth, subdu_width;
  double              start_node, start_val, start_deriv, end_node, end_val;
  slabs_weak_options_t *weak_options = slabs_options->slabs_weak_options;

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
  start_val = SLABS_SHELL_RADIUS_TOP;
  end_node = start_node + subdu_width / SLABS_MANTLE_DEPTH;
  end_val = start_val - subdu_depth / SLABS_MANTLE_DEPTH;

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

void
testtopo_temperature_set_fn (double *temp, double x, double y, double z,
                          ymir_locidx_t nid, void *data)
{
  if ( (fabs(y-0.5) < 0.2)  && (z > 0.8 && z < 0.9))
    *temp = 1.0;
  else
    *temp = 0.5;

  /* check temperature for `nan` and `inf` */
  RHEA_ASSERT (isfinite (*temp));

  /* bound temperature to valid interval */
  *temp = SC_MIN (1.0, *temp);
  *temp = SC_MAX (0.0, *temp);
}

void
testtopo2_temperature_set_fn (double *temp, double x, double y, double z,
                          ymir_locidx_t nid, void *data)
{
  slabs_topo_profile_t    *topo = (slabs_topo_profile_t *) data;
  int m, nsurf = topo->nsurf;
  double *tX = topo->tX;
  double *tY = topo->tY;
  double *tZ = topo->tZ;
  double factor;

  if (z > 0.95) {
    for (m = 0; m < nsurf; m++)  {
      if (fabs(x - tX[m]) < SC_1000_EPS &&
          fabs(y - tY[m]) < SC_1000_EPS)  {
        factor = tZ[m];
        if (z > factor)
          *temp = 100/8.4 + 0.5;
        else
          *temp = 0.5;
      }
    }
  }
  else if ( (fabs(y-0.5) < 0.2)  && (z > 0.8 && z < 0.9))
    *temp = .0;
  else
    *temp = 0.5;

  /* check temperature for `nan` and `inf` */
  RHEA_ASSERT (isfinite (*temp));

  /* bound temperature to valid interval */
  *temp = SC_MAX (0.0, *temp);
}


