/**************************************
 * Weakzone Computation
 *************************************/

/* Computes the weak zone factor depending on the distance to a weak zone.
 * The edges of the weak zone are smoothed by a Gaussian.
 * 1 - (1 - weak_factor) * exp ( - dist^2 / (2 * (0.5*thickness)^2) )
 * I added a scale of 20.0 so that the smoothing is sharper, at designited width, weak is aoubt exp(-10)~1-0.001~0.999*/
static inline double
slabs_weakzone_factor_fn (const double distance,
                          const double thickness,
                          const double thickness_const,
                          const double weak_factor)
{
  const double        d = distance - 0.5 * thickness_const;
  const double        std_dev = 0.5 * (thickness - thickness_const);

  RHEA_ASSERT (thickness_const <= thickness);

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
slabs_weakzone_subduct_dist_2plates_poly2 (double r, double lon,
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
  closest_pt = slabs_compute_closest_pt_on_poly2 (lon, r,
                                                  poly2_coeff, start_node,
                                                  start_val, start_deriv,
                                                  end_node, end_val,
                                                  &orientation_wrt_curve);

  /* compute distance */
  if (   orientation_wrt_curve == SLABS_CURVE_ORIENT_TOP
      || orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM
      || orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_BACK ) {
    /* compute distance to closest point on curve */
    dist = SLABS_MANTLE_DEPTH * sqrt (  SC_SQR (lon - closest_pt[0])
                                   + SC_SQR (r - closest_pt[1]) );
  }
  else if (   orientation_wrt_curve == SLABS_CURVE_ORIENT_TOP_RIGHT
           || orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_LEFT ) {
    /* compute distance to tip of curve */
    dist = SLABS_MANTLE_DEPTH * sqrt (  SC_SQR (end_node - lon)
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
slabs_weakzone_ridge_dist_2plates_poly2 (double r, double lon,
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
      dist = SLABS_MANTLE_DEPTH * (lon - end_node);
    }
    else if (lon <= end_node && r < end_val) {
      /* compute distance to bottom edge */
      dist = SLABS_MANTLE_DEPTH * (end_val - r);
    }
    else {
      /* compute distance to bottom right corner */
      dist = SLABS_MANTLE_DEPTH * sqrt (  SC_SQR (lon - end_node)
                                     + SC_SQR (end_val - r) );
    }
  }

  /* return distance to ridge weak zone */
  return dist;
}

/* Computes the weak fault zone for the shell slice with two plates. */
static double
slabs_weakzone_brick_2plates_poly2 (double r, double lon,
                              slabs_options_t * slabs_options)
{
  double              subdu_lon;
  double              subdu_dip_angle;
  double              subdu_depth, subdu_width;
  double              subdu_thickness, subdu_thickness_const;
  double              subdu_weak_factor;
  double              ridge_depth, ridge_width;
  double              ridge_smoothwidth;
  double              ridge_weak_factor;

  double              consider_thickness, courtesy_width;
  double              y_min = slabs_options->slabs_domain_options->y_min;
  double              start_node, start_val, start_deriv, end_node, end_val;
  double              *poly2_coeff;
  double              dist;
  double              weak = 1.0;
  double              s, a, b;

  slabs_weak_options_t *weak_options = slabs_options->slabs_weak_options;
  slabs_2plates_poly2_geo_coeff_t *geo = weak_options->weak_2plates_geo_coeff;

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
  consider_thickness = (2.0 * subdu_thickness)
                    / SLABS_MANTLE_DEPTH;
  if (   (  start_node -  0.5 * consider_thickness
          / sin (subdu_dip_angle / 180.0 * M_PI) ) <= lon
      && lon <= (end_node + 0.5 * consider_thickness )
      && (end_val - 0.5 * consider_thickness ) <= r ) {

    /* compute distance to subduction weak zone */
    dist = slabs_weakzone_subduct_dist_2plates_poly2 (r, lon, poly2_coeff,
                                                      start_node, start_val, start_deriv,
                                                      end_node, end_val);
    if (r > 0.95) {
      /*from r=0.9 to r=1, shrink the width of weakzone by a factor of 2*/
      s = - 10.0 * r + 10.5;
      a = subdu_thickness * s;
      b = subdu_thickness_const * s;
      weak = slabs_weakzone_factor_fn (dist, a,
                                       b, subdu_weak_factor);
    }
    else {
    /* set weak zone factor */
    weak = slabs_weakzone_factor_fn (dist, subdu_thickness,
                                     subdu_thickness_const, subdu_weak_factor);
    }
  }

  /*
   * set ridge weak zone in the left corner of domain, so that subducting
   * plate can "move"
   */

  /* set bottom left corner of weak zone */
  end_node = y_min + ridge_width / SLABS_MANTLE_DEPTH;
  end_val = SLABS_SHELL_RADIUS_TOP - ridge_depth / SLABS_MANTLE_DEPTH;

  /* only consider points close to weak zone */
  courtesy_width = 2.0 * ridge_smoothwidth / SLABS_MANTLE_DEPTH;
  if (lon <= (end_node + courtesy_width) && (end_val - courtesy_width) <= r) {
    /* compute distance to ridge weak zone */
    dist = slabs_weakzone_ridge_dist_2plates_poly2 (r, lon, end_node, end_val);

    /* set weak zone factor */
    weak = slabs_weakzone_factor_fn (dist, ridge_smoothwidth, 0.0,
                                     ridge_weak_factor);
  }

  /*
   * return weak zone factor
   */
  return weak;
}

static double
slabs_weakzone_node (const double x, const double y, const double z,
                     slabs_options_t *slabs_options)
{
  double              lon;
  double              r;

  /* compute radius and longitude */
  lon = y;
  r = z;

  /* compute weak zone factor */
  return slabs_weakzone_brick_2plates_poly2 (r, lon, slabs_options);
}

/* compute weak zone factor of an element */
void
slabs_weakzone_elem (double *_sc_restrict weak_elem,
                     const double *x,
                     const double *y,
                     const double *z,
                     const int n_nodes_per_el,
                     slabs_options_t *slabs_options)
{
  int            nodeid;
  double         z_mid = slabs_options->slabs_visc_options->z_lith;

  /* compute weak zone factor for each node */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) { /* loop over all
                                                         * nodes */
    weak_elem[nodeid] = slabs_weakzone_node (x[nodeid], y[nodeid], z[nodeid],
                                             slabs_options);

    /*temporary, to make sure that the viscosity along the weak zone does not have a sharp change
     * at the lith-athen boundary*/
    if (z[nodeid] < z_mid)  {
      weak_elem[nodeid] = 1.0;
    }
  }
}

void
slabs_weakzone_compute (ymir_dvec_t *weakzone, slabs_options_t *slabs_options)
{
  const char         *this_fn_name = "slabs_weakzone_compute";
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
  slabs_2plates_poly2_geo_coeff_t geo;
  slabs_weak_options_t *weak_options = slabs_options->slabs_weak_options;

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
  start_val = SLABS_SHELL_RADIUS_TOP;
  start_deriv = tan (-subdu_dip_angle / 180.0 * M_PI);
  end_node = start_node + subdu_width / SLABS_MANTLE_DEPTH;
  end_val = start_val - subdu_depth / SLABS_MANTLE_DEPTH;
  /* compute interpolating quadratic polynomial */
  poly2_coeff = slabs_compute_poly2_interpolation (start_node, start_val,
                                                   start_deriv,
                                                   end_node, end_val);

  RHEA_GLOBAL_PRODUCTIONF("weakzone: brick_2plates_poly2 downgoing slabs geometry\nstart_node: (%f %f), start_deriv: %f, end_node: (%f %f)\npoly2_coeff: (%f %f %f)\n", start_node, start_val, start_deriv, end_node, end_val, poly2_coeff[2], poly2_coeff[1], poly2_coeff[0]);

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
    slabs_weakzone_elem (weak_el_data, x, y, z,
                         n_nodes_per_el, slabs_options);

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
