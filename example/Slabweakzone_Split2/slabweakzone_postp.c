/* on each node, computes the shear and normal traction along
 * 2plates_poly2 weakzone in 2D Cartesian domain:
 *
 * compute the force vector acting on the fault plane,
 * then project the vector onto the normal and shear direction,
 * resulting in normal (plate coupling) and shear forces*/
void
slabs_postp_weakzone_coupling_brick_2plates_poly2 (double r, double lon,
                                                   double *shear, double *normal,
                                                   double *tau,
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

  double              consider_thickness;
  double              y_min = slabs_options->slabs_domain_options->y_min;
  double              start_node, start_val, start_deriv, end_node, end_val;
  double              *poly2_coeff;
  double              orth[3], tangent[3], tempvec[3];
  int                 orientation_wrt_curve;
  double             *closest_pt, dist=1.0, s;


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
  ridge_depth = weak_options->weakzone_2plates_ridge_depth;
  ridge_width = weak_options->weakzone_2plates_ridge_width;
  ridge_smoothwidth = weak_options->weakzone_2plates_ridge_smoothwidth;

  /* check parameters */
  RHEA_ASSERT (0.0 < subdu_lon);
  RHEA_ASSERT (0.0 < subdu_dip_angle);
  RHEA_ASSERT (0.0 < subdu_depth && 0.0 < subdu_width);
  RHEA_ASSERT (0.0 < subdu_thickness);
  RHEA_ASSERT (subdu_thickness_const <= subdu_thickness);
  RHEA_ASSERT (0.0 < ridge_depth && 0.0 < ridge_width);
  RHEA_ASSERT (0.0 <= ridge_smoothwidth);

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
  consider_thickness = (2.0 * subdu_thickness) / SLABS_MANTLE_DEPTH;
  * shear = * normal = 0.0;
  orth[0] = tangent[0] = 0.0;

  if (   (  start_node - 0.5 * consider_thickness
          / sin (subdu_dip_angle / 180.0 * M_PI) ) <= lon
      && lon <= (end_node + 0.5 * consider_thickness )
      && (end_val - consider_thickness ) <= r ) {
    /*
     * compute rotation of subduction weak zone from y=0 axis
     * */

    /* compute closest point on curve and orientation w.r.t. curve */
    closest_pt = slabs_compute_closest_pt_on_poly2 (lon, r,
                                                  poly2_coeff, start_node,
                                                  start_val, start_deriv,
                                                  end_node, end_val,
                                                  &orientation_wrt_curve);

    dist = slabs_weakzone_subduct_dist_2plates_poly2 (r, lon, poly2_coeff,
                                                        start_node, start_val, start_deriv,
                                                        end_node, end_val);

    if (orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_BACK)
      closest_pt[0] = start_node;
    else if (orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_LEFT ||
             orientation_wrt_curve == SLABS_CURVE_ORIENT_TOP_RIGHT)
      closest_pt[0] = end_node;

    if (r > 0.95) {
      /*from r=0.9 to r=1, shrink the width of weakzone by a factor of 2*/
      s = - 10.0 * r + 10.5;
    }
    else
      s = 1.0;

    if (dist < 0.5 * s * subdu_thickness)  {
      orth[1] = - slabs_poly2_deriv (closest_pt[0], poly2_coeff);
      orth[2] = 1.0 / sqrt(1.0 + orth[1] * orth[1]);
      orth[1] *= orth[2];
      tangent[2] = slabs_poly2_deriv (closest_pt[0], poly2_coeff);
      tangent[1] = 1.0 / sqrt (1.0 + tangent[2] * tangent[2]);
      tangent[2] *= tangent[1];

      tempvec[0] = ( tau[0] * orth[0] + tau[1] * orth[1] + tau[2] * orth[2]);
      tempvec[1] = ( tau[3] * orth[0] + tau[4] * orth[1] + tau[5] * orth[2]);
      tempvec[2] = ( tau[6] * orth[0] + tau[7] * orth[1] + tau[8] * orth[2]);

      * shear  = (tempvec[0] * tangent[0] + tempvec[1] * tangent[1] + tempvec[2] * tangent[2]);
      * normal = (tempvec[0] * orth[0] + tempvec[1] * orth[1] + tempvec[2] * orth[2]);
    }
    RHEA_FREE (closest_pt);
  }
}

/* compute weak zone factor for each node */
void
slabs_postp_weakzone_coupling_node (const double x, const double y, const double z,
                                    double *shear, double *normal,
                                    double *tau,
                                    slabs_options_t *slabs_options)
{
  double              lon;
  double              r;

  /* compute radius and longitude */
  lon = y;
  r = z;

  /* compute weak zone factor */
  slabs_postp_weakzone_coupling_brick_2plates_poly2 (r, lon, shear, normal,
                                                     tau, slabs_options);
}

/* compute weak zone factor of an element */
void
slabs_postp_weakzone_coupling_elem (sc_dmatrix_t * tau,
                                    sc_dmatrix_t * shear,
                                    sc_dmatrix_t * normal,
                                    const double *x,
                                    const double *y,
                                    const double *z,
                                    const int n_nodes_per_el,
                                    slabs_options_t *slabs_options)
{
  int            nodeid;

  /* compute weak zone factor for each node */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    double               *_sc_restrict tau_data = tau->e[0] + 9 * nodeid;
    double               *_sc_restrict shear_data = shear->e[0] + nodeid;
    double               *_sc_restrict normal_data = normal->e[0] + nodeid;
    slabs_postp_weakzone_coupling_node (x[nodeid], y[nodeid], z[nodeid],
                                        shear_data, normal_data, tau_data,
                                        slabs_options);
  }
}

/* compute the force vector acting on the fault plane,
 * then project the vector onto the normal and shear direction,
 * resulting in normal (plate coupling) and shear forces*/
static void
slabs_postp_weakzone_coupling_compute (ymir_cvec_t *vel, ymir_dvec_t *visc,
                                       ymir_dvec_t *svisc, ymir_dvec_t *TItens,
                                       ymir_dvec_t *normal, ymir_dvec_t *shear,
                                       ymir_velocity_elem_t *vel_elem,
                                       slabs_options_t *slabs_options)
{
  ymir_mesh_t         *mesh = ymir_vec_get_mesh (vel);
  const ymir_locidx_t n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  ymir_dvec_t         *tau = ymir_dvec_new (mesh, 9, YMIR_GAUSS_NODE);
  sc_dmatrix_t        *tau_el_mat, *shear_el_mat, *normal_el_mat;
  sc_dmatrix_t        *elemtau, *elemout_s, *elemout_n;
  double              *tau_el_data, *shear_el_data, *normal_el_data;
  double              *x, *y, *z, *tmp_el;
  ymir_locidx_t       elid;

  double              start_node, start_val, start_deriv;
  double              end_node, end_val;
  double              subdu_lon;
  double              subdu_dip_angle;
  double              subdu_depth, subdu_width;
  double              *poly2_coeff;
  slabs_2plates_poly2_geo_coeff_t geo;
  slabs_weak_options_t *weak_options = slabs_options->slabs_weak_options;
  const char          *this_fn_name = "slabs_postp_weakzone_coupling_compute";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  if (slabs_options->slabs_visc_options->viscosity_anisotropy
      == SLABS_VISC_TRANSVERSELY_ISOTROPY)  {
    slabs_stress_TI (vel, tau, visc, svisc, TItens, vel_elem);
  }
  else  {
    slabs_stress (vel, tau, visc, vel_elem);
  }

  /* set parameters according to weakzone options */
  subdu_lon = weak_options->weakzone_2plates_subdu_longitude;
  subdu_dip_angle = weak_options->weakzone_2plates_subdu_dip_angle;
  subdu_depth = weak_options->weakzone_2plates_subdu_depth;
  subdu_width = weak_options->weakzone_2plates_subdu_width;

  /* check parameters */
  RHEA_ASSERT (0.0 < subdu_lon);
  RHEA_ASSERT (0.0 < subdu_dip_angle);
  RHEA_ASSERT (0.0 < subdu_depth && 0.0 < subdu_width);

  /*
   * set subduction weak zone between plates
   */

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

  geo.start_node = start_node;
  geo.start_val = start_val;
  geo.start_deriv = start_deriv;
  geo.end_node = end_node;
  geo.end_val = end_val;
  geo.poly2_coeff = poly2_coeff;
  weak_options->weak_2plates_geo_coeff = &geo;

  /* create work variables */
  elemtau = sc_dmatrix_new (1, 9 * n_nodes_per_el);
  elemout_s = sc_dmatrix_new (1, n_nodes_per_el);
  elemout_n = sc_dmatrix_new (1, n_nodes_per_el);
  x = RHEA_ALLOC (double, n_nodes_per_el);
  y = RHEA_ALLOC (double, n_nodes_per_el);
  z = RHEA_ALLOC (double, n_nodes_per_el);
  tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    /* get coordinates of this element at Gauss nodes */
    ymir_mesh_get_elem_coord_gauss (x, y, z, elid, mesh, tmp_el);

    ymir_dvec_get_elem_interp (tau, elemtau, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (shear, elemout_s, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_WRITE);
    ymir_dvec_get_elem_interp (normal, elemout_n, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_WRITE);

    /* compute traction on each element */
    slabs_postp_weakzone_coupling_elem (elemtau, elemout_s, elemout_n,
                                        x, y, z, n_nodes_per_el, slabs_options);

    /* set traction of this element */
    ymir_dvec_set_elem_interp (shear, elemout_s, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);
    ymir_dvec_set_elem_interp (normal, elemout_n, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);
  }

  /* destroy */
  sc_dmatrix_destroy (elemtau);
  sc_dmatrix_destroy (elemout_s);
  sc_dmatrix_destroy (elemout_n);
  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);
  RHEA_FREE (poly2_coeff);

  ymir_vec_destroy (tau);
}


/* on each node, computes the shear and normal traction along
 * 2plates_poly2 weakzone in 2D Cartesian domain:
 *
 * compute the gradient of velocity (gradv) normal to and along the weakzone
 * by multiplying strain-rate (epsilon-dot)
 * with the unit director normal to and shear with the fault plane,
 * then multiplying the resulting gradvs with normal and shear viscosity, respectively
 * resulting in normal (plate coupling) and shear forces*/
void
slabs_postp_weakzone_traction_brick_2plates_poly2 (double r, double lon,
                                                   double *trac, double *gradv,
                                                   double *svisc, double *nvisc,
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

  double              consider_thickness;
  double              y_min = slabs_options->slabs_domain_options->y_min;
  double              start_node, start_val, start_deriv, end_node, end_val;
  double              *poly2_coeff;
  double              orth[3], tangent[3], tempvec[3];
  int                 orientation_wrt_curve;
  double             *closest_pt, dist=1.0, s;
  double              tempt, tt;


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
  ridge_depth = weak_options->weakzone_2plates_ridge_depth;
  ridge_width = weak_options->weakzone_2plates_ridge_width;
  ridge_smoothwidth = weak_options->weakzone_2plates_ridge_smoothwidth;

  /* check parameters */
  RHEA_ASSERT (0.0 < subdu_lon);
  RHEA_ASSERT (0.0 < subdu_dip_angle);
  RHEA_ASSERT (0.0 < subdu_depth && 0.0 < subdu_width);
  RHEA_ASSERT (0.0 < subdu_thickness);
  RHEA_ASSERT (subdu_thickness_const <= subdu_thickness);
  RHEA_ASSERT (0.0 < ridge_depth && 0.0 < ridge_width);
  RHEA_ASSERT (0.0 <= ridge_smoothwidth);

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
  consider_thickness = (2.0 * subdu_thickness) / SLABS_MANTLE_DEPTH;

  trac[0] = 0.0;
  trac[1] = 0.0;
  trac[2] = 0.0;
  orth[0] = tangent[0] = 0.0;

  if (   (  start_node - 0.5 * consider_thickness
          / sin (subdu_dip_angle / 180.0 * M_PI) ) <= lon
      && lon <= (end_node + 0.5 * consider_thickness )
      && (end_val - consider_thickness ) <= r ) {
    /*
     * compute rotation of subduction weak zone from y=0 axis
     * */

    /* compute closest point on curve and orientation w.r.t. curve */
    closest_pt = slabs_compute_closest_pt_on_poly2 (lon, r,
                                                    poly2_coeff, start_node,
                                                    start_val, start_deriv,
                                                    end_node, end_val,
                                                    &orientation_wrt_curve);
    dist = slabs_weakzone_subduct_dist_2plates_poly2 (r, lon, poly2_coeff,
                                                      start_node, start_val, start_deriv,
                                                      end_node, end_val);
    if (orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_BACK)
      closest_pt[0] = start_node;
    else if (orientation_wrt_curve == SLABS_CURVE_ORIENT_BOTTOM_LEFT ||
             orientation_wrt_curve == SLABS_CURVE_ORIENT_TOP_RIGHT)
      closest_pt[0] = end_node;

    if (r > 0.95) {
      /*from r=0.9 to r=1, shrink the width of weakzone by a factor of 2*/
      s = - 10.0 * r + 10.5;
    }
    else
      s = 1.0;

    if (dist < 0.5 * s * subdu_thickness)  {
//      orth[1] = - slabs_poly2_deriv (closest_pt[0], poly2_coeff);
//      orth[2] = 1.0 / sqrt(1.0 + orth[1] * orth[1]);
//      orth[1] *= orth[2];
//      tangent[2] = slabs_poly2_deriv (closest_pt[0], poly2_coeff);
//      tangent[1] = 1.0 / sqrt (1.0 + tangent[2] * tangent[2]);
//      tangent[2] *= tangent[1];

      tempt = - slabs_poly2_deriv (closest_pt[0], poly2_coeff);
      tt = atan(tempt);
      orth[1] = sin(tt);
      orth[2] = cos(tt);
      tangent[1] = cos(tt);
      tangent[2] = -sin(tt);

      tempvec[0] = ( gradv[0] * orth[0] + gradv[1] * orth[1] + gradv[2] * orth[2]);
      tempvec[1] = ( gradv[1] * orth[0] + gradv[3] * orth[1] + gradv[4] * orth[2]);
      tempvec[2] = ( gradv[2] * orth[0] + gradv[4] * orth[1] + gradv[5] * orth[2]);

      trac[1] = 2.0 * svisc[0] *
                (tempvec[0] * tangent[0] + tempvec[1] * tangent[1] + tempvec[2] * tangent[2]);
      trac[2] = 2.0 * nvisc[0] *
                (tempvec[0] * orth[0] + tempvec[1] * orth[1] + tempvec[2] * orth[2]);

    }
    RHEA_FREE (closest_pt);
  }
}

/* compute weak zone factor of an element */
void
slabs_postp_weakzone_traction_node (const double x, const double y, const double z,
                                    double *trac, double *gradv,
                                    double *svisc, double *nvisc,
                                    slabs_options_t *slabs_options)
{
  double              lon;
  double              r;

  /* compute radius and longitude */
  lon = y;
  r = z;

  /* compute weak zone factor */
  slabs_postp_weakzone_traction_brick_2plates_poly2 (r, lon, trac, gradv,
                                                     svisc, nvisc, slabs_options);
}

/* compute weak zone factor of an element */
void
slabs_postp_weakzone_traction_elem (sc_dmatrix_t *trac,
                                    sc_dmatrix_t *gradv,
                                    sc_dmatrix_t *svisc,
                                    sc_dmatrix_t *nvisc,
                                    const double *x,
                                    const double *y,
                                    const double *z,
                                    const int n_nodes_per_el,
                                    slabs_options_t *slabs_options)
{
  int            nodeid;

  /* compute weak zone factor for each node */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    double               *_sc_restrict trac_data = trac->e[0] + 3 * nodeid;
    double               *_sc_restrict gradv_data = gradv->e[0] + 6 * nodeid;
    double               *_sc_restrict svisc_data = svisc->e[0] + nodeid;
    double               *_sc_restrict nvisc_data = nvisc->e[0] + nodeid;
    slabs_postp_weakzone_traction_node (x[nodeid], y[nodeid], z[nodeid],
                                        trac_data, gradv_data,
                                        svisc_data, nvisc_data,
                                        slabs_options);
  }
}

/* compute the gradient of velocity (gradv) normal to and along the weakzone
 * by multiplying strain-rate (epsilon-dot)
 * with the unit director normal to and shear with the fault plane,
 * then multiplying the resulting gradvs with normal and shear viscosity, respectively
 * resulting in normal (plate coupling) and shear forces*/
static void
slabs_postp_weakzone_traction_compute (ymir_vec_t *vel, ymir_vec_t *traction,
                                       ymir_vec_t *svisc, ymir_vec_t *nvisc,
                                       slabs_options_t *slabs_options)
{
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (vel);
  const ymir_locidx_t n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  ymir_dvec_t        *gradvel = ymir_dvec_new (mesh, 6, YMIR_GAUSS_NODE);
  sc_dmatrix_t       *elemtrac, *elemgradv, *elemsvisc, *elemnvisc;
  double             *x, *y, *z, *tmp_el;
  ymir_locidx_t       elid;

  double              start_node, start_val, start_deriv;
  double              end_node, end_val;
  double              subdu_lon;
  double              subdu_dip_angle;
  double              subdu_depth, subdu_width;
  double              *poly2_coeff;
  slabs_2plates_poly2_geo_coeff_t geo;
  slabs_weak_options_t *weak_options = slabs_options->slabs_weak_options;
  const char         *this_fn_name = "slabs_postp_weakzone_traction_compute";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  ymir_velocity_strain_rate (vel, gradvel, 0);

  /* set parameters according to weakzone options */
  subdu_lon = weak_options->weakzone_2plates_subdu_longitude;
  subdu_dip_angle = weak_options->weakzone_2plates_subdu_dip_angle;
  subdu_depth = weak_options->weakzone_2plates_subdu_depth;
  subdu_width = weak_options->weakzone_2plates_subdu_width;

  /* check parameters */
  RHEA_ASSERT (0.0 < subdu_lon);
  RHEA_ASSERT (0.0 < subdu_dip_angle);
  RHEA_ASSERT (0.0 < subdu_depth && 0.0 < subdu_width);

  /*
   * set subduction weak zone between plates
   */

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

  geo.start_node = start_node;
  geo.start_val = start_val;
  geo.start_deriv = start_deriv;
  geo.end_node = end_node;
  geo.end_val = end_val;
  geo.poly2_coeff = poly2_coeff;
  weak_options->weak_2plates_geo_coeff = &geo;

  /* create work variables */
  elemtrac = sc_dmatrix_new (1, 3 * n_nodes_per_el);
  elemgradv = sc_dmatrix_new (1, 6 * n_nodes_per_el);
  elemsvisc = sc_dmatrix_new (1, n_nodes_per_el);
  elemnvisc = sc_dmatrix_new (1, n_nodes_per_el);
  x = RHEA_ALLOC (double, n_nodes_per_el);
  y = RHEA_ALLOC (double, n_nodes_per_el);
  z = RHEA_ALLOC (double, n_nodes_per_el);
  tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    /* get coordinates of this element at Gauss nodes */
    ymir_mesh_get_elem_coord_gauss (x, y, z, elid, mesh, tmp_el);

    ymir_dvec_get_elem_interp (gradvel, elemgradv, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (svisc, elemsvisc, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (nvisc, elemnvisc, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_READ);

    ymir_dvec_get_elem_interp (traction, elemtrac, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_WRITE);
    /* compute traction on each element */
    slabs_postp_weakzone_traction_elem (elemtrac, elemgradv, elemsvisc, elemnvisc,
                                        x, y, z, n_nodes_per_el, slabs_options);

    /* set traction of this element */
    ymir_dvec_set_elem_interp (traction, elemtrac, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);

    ymir_read_view_release (elemgradv);
    ymir_read_view_release (elemsvisc);
    ymir_read_view_release (elemnvisc);
  }

  /* destroy */
  sc_dmatrix_destroy (elemtrac);
  sc_dmatrix_destroy (elemgradv);
  sc_dmatrix_destroy (elemsvisc);
  sc_dmatrix_destroy (elemnvisc);
  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);
  RHEA_FREE (poly2_coeff);

  ymir_vec_destroy (gradvel);
}
