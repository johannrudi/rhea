#include <adjoint_postp.h>


/***********************************************************
 * Post-processing for 2nd invariant stress, traction, .etc.
 ***********************************************************/

/* compute traction as well as normal/shear stress at each element*/
void
subd_stress_elem (sc_dmatrix_t * in, sc_dmatrix_t * out,
                   sc_dmatrix_t * visc, ymir_velocity_elem_t * vel_elem,
                   double *_sc_restrict rxd, double *_sc_restrict sxd,
                   double *_sc_restrict txd, double *_sc_restrict ryd,
                   double *_sc_restrict syd, double *_sc_restrict tyd,
                   double *_sc_restrict rzd, double *_sc_restrict szd,
                   double *_sc_restrict tzd, sc_dmatrix_t * drst, sc_dmatrix_t * brst)
{
  const int           N = ymir_n (vel_elem->N);
  const int           Np = ymir_n (vel_elem->Np);
  int                 gp, i, j, k, l;
  double              S[9];
  double             *_sc_restrict viscd  = visc->e[0];

  sc_dmatrix_t       *tempvec1 = vel_elem->tempvec1;
  sc_dmatrix_t       *tempvec2 = vel_elem->tempvec2;
  sc_dmatrix_t       *tempvec3 = vel_elem->tempvec3;
  sc_dmatrix_t       *temptens = vel_elem->temptens1;
  ymir_derivative_elem_grad (N, 3, drst->e[0], brst->e[0],
                             rxd, sxd, txd, ryd, syd, tyd,
                             rzd, szd, tzd, in, temptens, tempvec1, tempvec2,
                             tempvec3, 0);

  /* create stress from duvw/dxyz * viscosity */
  for (gp = 0; gp < Np; gp++) {
    double               *_sc_restrict temptensd = temptens->e[0] + 9 * gp;
    double               *_sc_restrict outd  = out->e[0] + 9 * gp;

    S[0] = temptensd[0];
    S[1] = (temptensd[1] + temptensd[3]) * (1. / 2.);
    S[2] = (temptensd[2] + temptensd[6]) * (1. / 2.);
    S[3] = S[1];
    S[4] = temptensd[4];
    S[5] = (temptensd[5] + temptensd[7]) * (1. / 2.);
    S[6] = S[2];
    S[7] = S[5];
    S[8] = temptensd[8];

    /* compute linear combination of viscous stress and 3x3 tensor */
    for (i = 0; i < 9; i++)  {
      outd[i] = 2.0 * viscd[gp] * S[i];
    }
  }
}

/* compute traction as well as normal and shear stress*/
void
subd_stress (ymir_cvec_t * vel, ymir_dvec_t * tau,
              ymir_dvec_t * visc, ymir_velocity_elem_t * vel_elem)
{
  ymir_mesh_t         *mesh = vel->mesh;
  ymir_locidx_t       elid;
  const int           N  = ymir_n (mesh->cnodes->N);
  const int           Np = ymir_np (mesh->cnodes->N);
  const int           K  = mesh->cnodes->K;
  sc_dmatrix_t       *elemin = sc_dmatrix_new (1, 3 * Np);
  sc_dmatrix_t       *elemout = sc_dmatrix_new (1,9 * Np);
  sc_dmatrix_t       *elemvisc = sc_dmatrix_new (1, Np);
  const char          *this_fn_name = "subd_stress";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  ymir_dvec_set_zero (tau);

  for (elid = 0; elid < K; elid++)  {
    double             *_sc_restrict rxd = mesh->drdx->e[elid];
    double             *_sc_restrict sxd = mesh->dsdx->e[elid];
    double             *_sc_restrict txd = mesh->dtdx->e[elid];
    double             *_sc_restrict ryd = mesh->drdy->e[elid];
    double             *_sc_restrict syd = mesh->dsdy->e[elid];
    double             *_sc_restrict tyd = mesh->dtdy->e[elid];
    double             *_sc_restrict rzd = mesh->drdz->e[elid];
    double             *_sc_restrict szd = mesh->dsdz->e[elid];
    double             *_sc_restrict tzd = mesh->dtdz->e[elid];
    double             *_sc_restrict Jdetd = mesh->Jdet->e[elid];

    ymir_cvec_get_elem_interp (vel, elemin, YMIR_STRIDE_NODE, elid,
                               YMIR_GLL_NODE, YMIR_COPY);
    ymir_dvec_get_elem_interp (visc, elemvisc, YMIR_STRIDE_COMP, elid,
                              YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (tau, elemout, YMIR_STRIDE_NODE, elid,
                             YMIR_GAUSS_NODE, YMIR_WRITE);

    subd_stress_elem (elemin, elemout,
                       elemvisc,
                       vel_elem, rxd, sxd, txd, ryd,
                       syd, tyd, rzd, szd, tzd, mesh->drst, mesh->brst);

    ymir_dvec_set_elem_interp (tau, elemout, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);

    ymir_read_view_release (elemvisc);
  }

  sc_dmatrix_destroy (elemin);
  sc_dmatrix_destroy (elemout);
  sc_dmatrix_destroy (elemvisc);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/* Computes the shear and normal traction along 2plates_poly2 weakzone in 2D Cartesian domain*/
void
subd_postp_weakzone_coupling_brick_2plates_poly2 (double r, double lon,
                                                   double *shear, double *normal,
                                                   double *tau,
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

  double              consider_thickness;
  double              y_min = subd_options->domain_options->y_min;
  double              start_node, start_val, start_deriv, end_node, end_val;
  double              *poly2_coeff;
  double              orth[3], tangent[3], tempvec[3];
  int                 orientation_wrt_curve;
  double             *closest_pt, dist=1.0, s;


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
  consider_thickness = (2.0 * subdu_thickness) / SUBD_MANTLE_DEPTH;
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
    closest_pt = subd_compute_closest_pt_on_poly2 (lon, r,
                                                  poly2_coeff, start_node,
                                                  start_val, start_deriv,
                                                  end_node, end_val,
                                                  &orientation_wrt_curve);

    dist = subd_weakzone_subduct_dist_2plates_poly2 (r, lon, poly2_coeff,
                                                        start_node, start_val, start_deriv,
                                                        end_node, end_val);

    if (orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_BOTTOM_BACK)
      closest_pt[0] = start_node;
    else if (orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_BOTTOM_LEFT ||
             orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_TOP_RIGHT)
      closest_pt[0] = end_node;

    if (r > 0.95) {
      /*from r=0.9 to r=1, shrink the width of weakzone by a factor of 2*/
      s = - 10.0 * r + 10.5;
    }
    else
      s = 1.0;

    if (dist < 0.5 * s * subdu_thickness)  {
      orth[1] = - subd_poly2_deriv (closest_pt[0], poly2_coeff);
      orth[2] = 1.0 / sqrt(1.0 + orth[1] * orth[1]);
      orth[1] *= orth[2];
      tangent[2] = subd_poly2_deriv (closest_pt[0], poly2_coeff);
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

void
subd_postp_weakzone_coupling_node (const double x, const double y, const double z,
                                    double *shear, double *normal,
                                    double *tau,
                                    subd_options_t *subd_options)
{
  double              lon;
  double              r;

  /* compute radius and longitude */
  lon = y;
  r = z;

  /* compute weak zone factor */
  subd_postp_weakzone_coupling_brick_2plates_poly2 (r, lon, shear, normal,
                                                     tau, subd_options);
}

/* compute weak zone factor of an element */
void
subd_postp_weakzone_coupling_elem (sc_dmatrix_t * tau,
                                    sc_dmatrix_t * shear,
                                    sc_dmatrix_t * normal,
                                    const double *x,
                                    const double *y,
                                    const double *z,
                                    const int n_nodes_per_el,
                                    subd_options_t *subd_options)
{
  int            nodeid;

  /* compute weak zone factor for each node */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    double               *_sc_restrict tau_data = tau->e[0] + 9 * nodeid;
    double               *_sc_restrict shear_data = shear->e[0] + nodeid;
    double               *_sc_restrict normal_data = normal->e[0] + nodeid;
    subd_postp_weakzone_coupling_node (x[nodeid], y[nodeid], z[nodeid],
                                        shear_data, normal_data, tau_data,
                                        subd_options);
  }
}

void
subd_postp_weakzone_coupling_compute (ymir_cvec_t *vel, ymir_dvec_t *visc,
                                       ymir_dvec_t *svisc, ymir_dvec_t *TItens,
                                       ymir_dvec_t *normal, ymir_dvec_t *shear,
                                       ymir_velocity_elem_t *vel_elem,
                                       subd_options_t *subd_options)
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
  subd_2plates_poly2_geo_coeff_t geo;
  subd_weak_options_t *weak_options = subd_options->weak_options;
  const char          *this_fn_name = "subd_postp_weakzone_coupling_compute";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  if (subd_options->visc_options->anisotropy_type
      == SUBD_VISC_TRANSVERSELY_ISOTROPY)  {
    subd_stress_TI (vel, tau, visc, svisc, TItens, vel_elem);
  }
  else  {
    subd_stress (vel, tau, visc, vel_elem);
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
  start_val = SUBD_SHELL_RADIUS_TOP;
  start_deriv = tan (-subdu_dip_angle / 180.0 * M_PI);
  end_node = start_node + subdu_width / SUBD_MANTLE_DEPTH;
  end_val = start_val - subdu_depth / SUBD_MANTLE_DEPTH;
  /* compute interpolating quadratic polynomial */
  poly2_coeff = subd_compute_poly2_interpolation (start_node, start_val,
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
    subd_postp_weakzone_coupling_elem (elemtau, elemout_s, elemout_n,
                                        x, y, z, n_nodes_per_el, subd_options);

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


/* Computes the shear and normal traction along 2plates_poly2 weakzone in 2D Cartesian domain*/
void
subd_postp_weakzone_traction_brick_2plates_poly2 (double r, double lon,
                                                   double *trac, double *gradv,
                                                   double *svisc, double *nvisc,
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

  double              consider_thickness;
  double              y_min = subd_options->domain_options->y_min;
  double              start_node, start_val, start_deriv, end_node, end_val;
  double              *poly2_coeff;
  double              orth[3], tangent[3], tempvec[3];
  int                 orientation_wrt_curve;
  double             *closest_pt, dist=1.0, s;
  double              tempt, tt;


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
  consider_thickness = (2.0 * subdu_thickness) / SUBD_MANTLE_DEPTH;

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
    closest_pt = subd_compute_closest_pt_on_poly2 (lon, r,
                                                    poly2_coeff, start_node,
                                                    start_val, start_deriv,
                                                    end_node, end_val,
                                                    &orientation_wrt_curve);
    dist = subd_weakzone_subduct_dist_2plates_poly2 (r, lon, poly2_coeff,
                                                      start_node, start_val, start_deriv,
                                                      end_node, end_val);
    if (orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_BOTTOM_BACK)
      closest_pt[0] = start_node;
    else if (orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_BOTTOM_LEFT ||
             orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_TOP_RIGHT)
      closest_pt[0] = end_node;

    if (r > 0.95) {
      /*from r=0.9 to r=1, shrink the width of weakzone by a factor of 2*/
      s = - 10.0 * r + 10.5;
    }
    else
      s = 1.0;

    if (dist < 0.5 * s * subdu_thickness)  {
//      orth[1] = - subd_poly2_deriv (closest_pt[0], poly2_coeff);
//      orth[2] = 1.0 / sqrt(1.0 + orth[1] * orth[1]);
//      orth[1] *= orth[2];
//      tangent[2] = subd_poly2_deriv (closest_pt[0], poly2_coeff);
//      tangent[1] = 1.0 / sqrt (1.0 + tangent[2] * tangent[2]);
//      tangent[2] *= tangent[1];

      tempt = - subd_poly2_deriv (closest_pt[0], poly2_coeff);
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

void
subd_postp_weakzone_traction_node (const double x, const double y, const double z,
                                    double *trac, double *gradv,
                                    double *svisc, double *nvisc,
                                    subd_options_t *subd_options)
{
  double              lon;
  double              r;

  /* compute radius and longitude */
  lon = y;
  r = z;

  /* compute weak zone factor */
  subd_postp_weakzone_traction_brick_2plates_poly2 (r, lon, trac, gradv,
                                                     svisc, nvisc, subd_options);
}

/* compute weak zone factor of an element */
void
subd_postp_weakzone_traction_elem (sc_dmatrix_t *trac,
                                    sc_dmatrix_t *gradv,
                                    sc_dmatrix_t *svisc,
                                    sc_dmatrix_t *nvisc,
                                    const double *x,
                                    const double *y,
                                    const double *z,
                                    const int n_nodes_per_el,
                                    subd_options_t *subd_options)
{
  int            nodeid;

  /* compute weak zone factor for each node */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    double               *_sc_restrict trac_data = trac->e[0] + 3 * nodeid;
    double               *_sc_restrict gradv_data = gradv->e[0] + 6 * nodeid;
    double               *_sc_restrict svisc_data = svisc->e[0] + nodeid;
    double               *_sc_restrict nvisc_data = nvisc->e[0] + nodeid;
    subd_postp_weakzone_traction_node (x[nodeid], y[nodeid], z[nodeid],
                                        trac_data, gradv_data,
                                        svisc_data, nvisc_data,
                                        subd_options);
  }
}

void
subd_postp_weakzone_traction_compute (ymir_vec_t *vel, ymir_vec_t *traction,
                                       ymir_vec_t *svisc, ymir_vec_t *nvisc,
                                       subd_options_t *subd_options)
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
  subd_2plates_poly2_geo_coeff_t geo;
  subd_weak_options_t *weak_options = subd_options->weak_options;
  const char         *this_fn_name = "subd_postp_weakzone_traction_compute";

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
  start_val = SUBD_SHELL_RADIUS_TOP;
  start_deriv = tan (-subdu_dip_angle / 180.0 * M_PI);
  end_node = start_node + subdu_width / SUBD_MANTLE_DEPTH;
  end_val = start_val - subdu_depth / SUBD_MANTLE_DEPTH;
  /* compute interpolating quadratic polynomial */
  poly2_coeff = subd_compute_poly2_interpolation (start_node, start_val,
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
    subd_postp_weakzone_traction_elem (elemtrac, elemgradv, elemsvisc, elemnvisc,
                                        x, y, z, n_nodes_per_el, subd_options);

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

/* compute traction. It is an alternative approach that takes advantage of an existing
   subroutine ymir_velocity_strain_rate and directly compute traction on each node*/
void
subd_traction (ymir_cvec_t * vel, ymir_dvec_t *traction,
                  double *n_dir, ymir_dvec_t *visc)
{
  ymir_mesh_t        *mesh = vel->mesh;
  const int           N = ymir_n (mesh->cnodes->N);
  const int           Np = (N + 1)*(N+1)*(N+1);
  ymir_locidx_t       K = mesh->cnodes->K;
  ymir_dvec_t         *tau_tensor = ymir_dvec_new (mesh, 6,
                                                     YMIR_GAUSS_NODE);
  ymir_locidx_t       elid;
  sc_dmatrix_t       *an = sc_dmatrix_new (0, 0);
  sc_dmatrix_t       *bn = sc_dmatrix_new (0, 0);
  int                 i, j;

  ymir_velocity_strain_rate (vel, tau_tensor, 0);
  ymir_dvec_multiply_in1 (visc, tau_tensor);
  for (elid = 0; elid < K; elid++) {
    for (j = 0; j < Np; j++) {
      double              val = 0.;
      ymir_dvec_get_node (tau_tensor, an, elid, j, YMIR_READ);
      ymir_dvec_get_node (traction,   bn, elid, j, YMIR_WRITE);
      bn->e[0][0] = an->e[0][0] * n_dir[0] + an->e[0][1] * n_dir[1] + an->e[0][2] * n_dir[2];
      bn->e[0][1] = an->e[0][1] * n_dir[0] + an->e[0][3] * n_dir[1] + an->e[0][4] * n_dir[2];
      bn->e[0][2] = an->e[0][2] * n_dir[0] + an->e[0][4] * n_dir[1] + an->e[0][5] * n_dir[2];
      ymir_read_view_release (an);
      ymir_dvec_set_node (traction, bn, elid, j, YMIR_SET);
    }
  }
  sc_dmatrix_destroy (an);
  sc_dmatrix_destroy (bn);
}

/* compute traction as well as normal/shear stress at each element*/
void
subd_normal_stress_elem (sc_dmatrix_t * in, sc_dmatrix_t * out1, sc_dmatrix_t * out2,
                            sc_dmatrix_t * out3, double * n_dir,
                  sc_dmatrix_t * visc_mat, ymir_velocity_elem_t * vel_elem,
                  double *_sc_restrict rxd, double *_sc_restrict sxd,
                  double *_sc_restrict txd, double *_sc_restrict ryd,
                  double *_sc_restrict syd, double *_sc_restrict tyd,
                  double *_sc_restrict rzd, double *_sc_restrict szd,
                  double *_sc_restrict tzd, sc_dmatrix_t * drst, sc_dmatrix_t * brst)
{
  const int           N = ymir_n (vel_elem->N);
  int                 gp, i, j, k, l;
  double              tempmatd[9];
  double              tempvec[3];
  double              temp;
  double             *_sc_restrict viscd = visc_mat->e[0];
  sc_dmatrix_t       *tempvec1 = vel_elem->tempvec1;
  sc_dmatrix_t       *tempvec2 = vel_elem->tempvec2;
  sc_dmatrix_t       *tempvec3 = vel_elem->tempvec3;
  sc_dmatrix_t       *temptens = vel_elem->temptens1;

  ymir_derivative_elem_grad (N, 3, drst->e[0], brst->e[0],
                             rxd, sxd, txd, ryd, syd, tyd,
                             rzd, szd, tzd, in, temptens, tempvec1, tempvec2,
                             tempvec3, 0);

  /* create stress from duvw/dxyz * viscosity */
  for (gp = 0, k = 0; k < N + 1; k++) {
    for (j = 0; j < N + 1; j++) {
      for (i = 0; i < N + 1; i++, gp++) {
        double               *_sc_restrict temptensd = temptens->e[0] + 9 * gp;
        double               *_sc_restrict normal = out1->e[0] + gp;
        double               *_sc_restrict shear  = out2->e[0] + gp;
        double               *_sc_restrict trac   = out3->e[0] + 3 * gp;

        tempmatd[0] = temptensd[0];
        tempmatd[1] = (temptensd[1] + temptensd[3]) * (1. / 2.);
        tempmatd[2] = (temptensd[2] + temptensd[6]) * (1. / 2.);
        tempmatd[3] = tempmatd[1];
        tempmatd[4] = temptensd[4];
        tempmatd[5] = (temptensd[5] + temptensd[7]) * (1. / 2.);
        tempmatd[6] = tempmatd[2];
        tempmatd[7] = tempmatd[5];
        tempmatd[8] = temptensd[8];
        *(trac++) = tempvec[0] = ( tempmatd[0] * n_dir[0]
                                 + tempmatd[1] * n_dir[1]
                                 + tempmatd[2] * n_dir[2]) * 2.0 * viscd[gp];
        *(trac++) = tempvec[1] = ( tempmatd[3] * n_dir[0]
                                 + tempmatd[4] * n_dir[1]
                                 + tempmatd[5] * n_dir[2]) * 2.0 * viscd[gp];
        *(trac++) = tempvec[2] = ( tempmatd[6] * n_dir[0]
                                 + tempmatd[7] * n_dir[1]
                                 + tempmatd[8] * n_dir[2]) * 2.0 * viscd[gp];
        *normal = temp = tempvec[0] * n_dir[0] + tempvec[1] * n_dir[1] + tempvec[2] * n_dir[2];
        *shear  = sqrt( tempvec[0] * tempvec[0] + tempvec[1] * tempvec[1] + tempvec[2] * tempvec[2]
                        - temp * temp);
      }
    }
  }
}

/* compute traction as well as normal and shear stress*/
void
subd_normal_stress (ymir_cvec_t * vel, ymir_dvec_t * n_tau, ymir_dvec_t * s_tau,
                        ymir_dvec_t * traction, double * n_dir,
                       ymir_dvec_t * visc, ymir_velocity_elem_t * vel_elem)
{
  ymir_mesh_t         *mesh = vel->mesh;
  ymir_locidx_t       elid;
  const int           N  = ymir_n (mesh->cnodes->N);
  const int           Np = ymir_np (mesh->cnodes->N);
  const int           K  = mesh->cnodes->K;
  sc_dmatrix_t       *elemin = sc_dmatrix_new (1, 3 * Np);
  sc_dmatrix_t       *elemout1 = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemout2 = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemout3 = sc_dmatrix_new (1, 3 * Np);
  sc_dmatrix_t       *elemvisc = sc_dmatrix_new (1, Np);

  ymir_dvec_set_zero (traction);
  ymir_dvec_set_zero (n_tau);
  ymir_dvec_set_zero (s_tau);

  for (elid = 0; elid < K; elid++)  {
    double             *_sc_restrict rxd = mesh->drdx->e[elid];
    double             *_sc_restrict sxd = mesh->dsdx->e[elid];
    double             *_sc_restrict txd = mesh->dtdx->e[elid];
    double             *_sc_restrict ryd = mesh->drdy->e[elid];
    double             *_sc_restrict syd = mesh->dsdy->e[elid];
    double             *_sc_restrict tyd = mesh->dtdy->e[elid];
    double             *_sc_restrict rzd = mesh->drdz->e[elid];
    double             *_sc_restrict szd = mesh->dsdz->e[elid];
    double             *_sc_restrict tzd = mesh->dtdz->e[elid];
    double             *_sc_restrict Jdetd = mesh->Jdet->e[elid];

    ymir_cvec_get_elem_interp (vel, elemin, YMIR_STRIDE_NODE, elid,
                               YMIR_GLL_NODE, YMIR_COPY);
    ymir_dvec_get_elem_interp (visc, elemvisc, YMIR_STRIDE_COMP, elid,
                              YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (n_tau, elemout1, YMIR_STRIDE_COMP, elid,
                             YMIR_GAUSS_NODE, YMIR_WRITE);
    ymir_dvec_get_elem_interp (s_tau, elemout2, YMIR_STRIDE_COMP, elid,
                             YMIR_GAUSS_NODE, YMIR_WRITE);

    subd_normal_stress_elem (elemin, elemout1, elemout2, elemout3, n_dir, elemvisc,
                            vel_elem, rxd, sxd, txd, ryd,
                            syd, tyd, rzd, szd, tzd, mesh->drst, mesh->brst);

    ymir_dvec_set_elem_interp (n_tau, elemout1, YMIR_STRIDE_COMP, elid,
                             YMIR_GAUSS_NODE, YMIR_SET);
    ymir_dvec_set_elem_interp (s_tau, elemout2, YMIR_STRIDE_COMP, elid,
                             YMIR_GAUSS_NODE, YMIR_SET);
    ymir_dvec_set_elem_interp (traction, elemout3, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);

    ymir_read_view_release (elemvisc);
  }

  sc_dmatrix_destroy (elemin);
  sc_dmatrix_destroy (elemout1);
  sc_dmatrix_destroy (elemout2);
  sc_dmatrix_destroy (elemout3);
  sc_dmatrix_destroy (elemvisc);
}

static void
subd_physics_normal_boundary_stress_fn (double *stress_norm,
                                         double x, double y, double z,
                                         double nx, double ny, double nz,
                                         ymir_topidx_t face,
                                         ymir_locidx_t node_id,
                                         void *data)
{
  ymir_vec_t         *vec_bndr = (ymir_vec_t *) data;
  double             *v = ymir_cvec_index (vec_bndr, node_id, 0);

  RHEA_ASSERT (vec_bndr->ncfields == 3);

  /* compute inner product with boundary outer normal vector */
  *stress_norm = nx * v[0] + ny * v[1] + nz * v[2];
}

/* normal stress on the surface */
void
subd_physics_compute_normal_boundary_stress (ymir_vec_t *stress_bndr_norm,
                                              ymir_vec_t *up,
                                              ymir_vec_t *rhs_u_point,
                                              ymir_stokes_op_t *stokes_op)
{
  ymir_mesh_t        *mesh = up->mesh;
  ymir_pressure_elem_t  *press_elem = stokes_op->press_elem;
  ymir_stress_op_t   *stress_op = stokes_op->stress_op;
  const int           skip_dir = stress_op->skip_dir;
  const ymir_topidx_t face_id = stress_bndr_norm->meshnum;

  ymir_vec_t         *rhs = ymir_stokes_vec_new (mesh, press_elem);
  ymir_vec_t         *residual_up = ymir_stokes_vec_new (mesh, press_elem);
  ymir_vec_t         *residual_u = ymir_cvec_new (mesh, 3);
  ymir_vec_t         *residual_bndr = ymir_face_cvec_new (mesh, face_id, 3);
  ymir_vec_t         *mass_lump_boundary;

  /* check input */
  YMIR_ASSERT_IS_CVEC (stress_bndr_norm);
  RHEA_ASSERT (stress_bndr_norm->ncfields == 1);
  RHEA_ASSERT (ymir_stokes_vec_is_stokes_vec (up));
  YMIR_ASSERT_IS_CVEC (rhs_u_point);
  RHEA_ASSERT (ymir_vec_is_not_dirty (up));
//  RHEA_ASSERT (ymir_vec_is_not_dirty (rhs_u_point));

  /* construct the right-hand side */
  ymir_stokes_pc_construct_rhs (rhs, rhs_u_point, NULL, NULL,
                                1 /* incompressible */, stokes_op, 0);
  RHEA_ASSERT (sc_dmatrix_is_valid (rhs->dataown));
  RHEA_ASSERT (sc_dmatrix_is_valid (rhs->coff));
  RHEA_ASSERT (ymir_vec_is_not_dirty (rhs));

  /* turn off boundary constraints */
  stress_op->skip_dir = 1;

  /* compute (unconstrained) residual
   *   r_mom  = a * u + b^t * p - f
   *   r_mass = b * u
   */
  ymir_stokes_pc_apply_stokes_op (up, residual_up, stokes_op, 0, 0);
  ymir_vec_add (-1.0, rhs, residual_up);
  RHEA_ASSERT (sc_dmatrix_is_valid (residual_up->dataown));
  RHEA_ASSERT (sc_dmatrix_is_valid (residual_up->coff));
  ymir_vec_destroy (rhs);

  /* restore boundary constraints */
  stress_op->skip_dir = skip_dir;

  /* get the velocity component of the residual */
  ymir_stokes_vec_get_velocity (residual_up, residual_u, stokes_op->press_elem);

  /* interpolate residual onto boundary */
  ymir_interp_vec (residual_u, residual_bndr);
  RHEA_ASSERT (sc_dmatrix_is_valid (residual_bndr->dataown));
  RHEA_ASSERT (sc_dmatrix_is_valid (residual_bndr->coff));
  ymir_vec_destroy (residual_up);
  ymir_vec_destroy (residual_u);

 /* get the normal part of the residual */
  ymir_face_cvec_set_function (stress_bndr_norm,
                               subd_physics_normal_boundary_stress_fn,
                               residual_bndr);
  ymir_vec_destroy (residual_bndr);

  /* invert mass matrix on boundary */
  mass_lump_boundary = ymir_face_cvec_new (mesh, face_id, 1);
  ymir_mass_lump (mass_lump_boundary);
  ymir_vec_divide_in (mass_lump_boundary, stress_bndr_norm);
  ymir_vec_destroy (mass_lump_boundary);
}

void
subd_postp_topography (ymir_vec_t *topography,
                      rhea_stokes_problem_t *stokes_problem,
                      subd_options_t *subd_opt)
{
  ymir_mesh_t           *ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  ymir_vec_t            *surf_normal_stress = ymir_face_cvec_new (ymir_mesh,
                                                   RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
  ymir_vec_t            *vec_e = ymir_face_cvec_new (ymir_mesh,
                                                   RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
  ymir_vec_t            *masse = ymir_face_cvec_new (ymir_mesh,
                                                   RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
  ymir_vec_t            *massu = ymir_face_cvec_new (ymir_mesh,
                                                   RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
  double                avg_stress, tempe, tempu;
  double                therm_expa = subd_opt->para_options->ref_therm_expa;
  double                temp_diff = subd_opt->para_options->ref_temp_diff;
  double                rayleigh = subd_opt->para_options->rayleigh;
  double                radius_max = subd_opt->domain_options->radius_max;


      /* compute surface normal stress sigma */
  subd_physics_compute_normal_boundary_stress (
                 surf_normal_stress,
                 rhea_stokes_problem_get_velocity_pressure (stokes_problem),
                 rhea_stokes_problem_get_rhs_vel (stokes_problem),
                 rhea_stokes_problem_get_stokes_op (stokes_problem));

  /*compute average surface normal stress
   * vec_e is 1.0
   * average = (e, M*u)/(e, M*e)*/
  ymir_vec_set_value(vec_e, 1.0)
  ymir_mass_apply (surf_normal_stress, massu);
  ymir_mass_apply (vec_e, masse);
  tempu = ymir_vec_innerprod (vec_e, massu);
  tempe = ymir_vec_innerprod (vec_e, masse);
  avg_stress = tempu/tempe;

  /*scale normal stress to topography*/
  ymir_vec_copy (surf_normal_stress, topography);
  ymir_vec_shift (-avg_stress, topography);
  /*vec_topo = (-alpha*DeltaT/Ra) * (nstress-avg_stress) + 1.0, here Ra=1 */
  ymir_vec_scale_shift (-(therm_expa * temp_diff / rayleigh), radius_max, topography);

  ymir_vec_destroy (surf_normal_stress);
  ymir_vec_destroy (vec_e);
  ymir_vec_destroy (masse);
  ymir_vec_destroy (massu);

}


