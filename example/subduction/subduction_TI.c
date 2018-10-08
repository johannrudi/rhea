#include <subduction_TI.h>

/******************************************************
 *  Transversely Isotropic Viscosity
*******************************************************/

/* get shear viscosity for transversely isotropy model from stress operator*/
void
subd_stress_op_copy_shear_visc (ymir_vec_t *shear_visc,
                                   ymir_stress_op_t *stress_op)
{
  /* check input */
  RHEA_ASSERT (stress_op->coeff_TI_svisc != NULL);

  /* copy Stokes coefficient and divide by 2 */
  ymir_vec_copy (stress_op->coeff_TI_svisc, shear_visc);
  ymir_vec_scale (0.5, shear_visc);
}

/* get TI tensor for transversely isotropy model from stress operator*/
void
subd_stress_op_copy_TI_tensor (ymir_vec_t *TI_tensor,
                                   ymir_stress_op_t *stress_op)
{
  /* check input */
  RHEA_ASSERT (stress_op->coeff_TI_tensor != NULL);

  ymir_vec_copy (stress_op->coeff_TI_tensor, TI_tensor);
}

/* Computes the rotation from x axis of weak fault zone with two plates. */
static double
subd_TI_rotation_brick_2plates_poly2 (double r, double lon,
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

  double              consider_thickness, courtesy_width;
  double              y_min = subd_options->domain_options->y_min;
  double              start_node, start_val, start_deriv, end_node, end_val;
  double              *poly2_coeff;
  double              rot = 0.0;
  int                 orientation_wrt_curve;
  double             *closest_pt, dist=1.0;

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
  consider_thickness = (2.0 * subdu_thickness)
                    / SUBD_MANTLE_DEPTH;
  if (   (  start_node - 0.5 * consider_thickness
          / sin (subdu_dip_angle / 180.0 * M_PI)) <= lon
      && lon <= (end_node + 0.5 * consider_thickness)
      && (end_val - consider_thickness) <= r ) {
    /*
     * compute rotation of subduction weak zone from y=0 axis
     * */

    /* compute closest point on curve and orientation w.r.t. curve */
    closest_pt = subd_compute_closest_pt_on_poly2 (lon, r,
                                                  poly2_coeff, start_node,
                                                  start_val, start_deriv,
                                                  end_node, end_val,
                                                  &orientation_wrt_curve);

    if (orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_BOTTOM_BACK)
      closest_pt[0] = start_node;
    else if (orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_BOTTOM_LEFT ||
             orientation_wrt_curve == SUBD_SLAB_CURVE_ORIENT_TOP_RIGHT)
      closest_pt[0] = end_node;

    rot = subd_compute_rot_on_poly2 (poly2_coeff, closest_pt);

    RHEA_FREE (closest_pt);
  }

  /*
   * return weak zone factor
   */
  return rot;
}

static double
subd_TI_rotation_node (const double x, const double y, const double z,
                        subd_options_t *subd_options)
{
  double              lon;
  double              r;

  /* compute radius and longitude */
  lon = y;
  r = z;

  /* compute weak zone factor */
  return subd_TI_rotation_brick_2plates_poly2 (r, lon, subd_options);
}

/* compute weak zone factor of an element */
void
subd_TI_rotation_elem (double *_sc_restrict rot_elem,
                        double *_sc_restrict weak_elem,
                         const double *x,
                         const double *y,
                         const double *z,
                         const int n_nodes_per_el,
                         subd_options_t *subd_options)
{
  int            nodeid;

  /* compute weak zone factor for each node */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    if (fabs(1.0 - weak_elem[nodeid]) < 1e-3)  {
      rot_elem[nodeid] = .0;
    }
    else {
      rot_elem[nodeid] = subd_TI_rotation_node (x[nodeid], y[nodeid], z[nodeid],
                                                 subd_options);
    }
  }
}

static void
subd_TI_rotation_compute (ymir_vec_t *rotate, ymir_vec_t *weak,
                           subd_options_t *subd_options)
{
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (rotate);
  const ymir_locidx_t n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  sc_dmatrix_t       *rot_el_mat, *weak_el_mat;
  double             *rot_el_data, *weak_el_data;
  double             *x, *y, *z, *tmp_el;
  ymir_locidx_t      elid;

  double              start_node, start_val, start_deriv;
  double              end_node, end_val;
  double              subdu_lon;
  double              subdu_dip_angle;
  double              subdu_depth, subdu_width;
  double              *poly2_coeff;
  subd_2plates_poly2_geo_coeff_t geo;
  subd_weak_options_t *weak_options = subd_options->weak_options;
  const char         *this_fn_name = "subd_TI_rotation_compute";

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

  /* check input */
  RHEA_ASSERT (rotate->ndfields == 1);

  /* create work variables */
  rot_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  rot_el_data = rot_el_mat->e[0];
  weak_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  x = RHEA_ALLOC (double, n_nodes_per_el);
  y = RHEA_ALLOC (double, n_nodes_per_el);
  z = RHEA_ALLOC (double, n_nodes_per_el);
  tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    /* get coordinates of this element at Gauss nodes */
    ymir_mesh_get_elem_coord_gauss (x, y, z, elid, mesh, tmp_el);

    weak_el_data = rhea_viscosity_get_elem_gauss (weak_el_mat, weak, elid);

    /* compute weak zone factor */
    subd_TI_rotation_elem (rot_el_data, weak_el_data, x, y, z,
                            n_nodes_per_el, subd_options);

    /* set weak zone of this element */
    rhea_viscosity_set_elem_gauss (rotate, rot_el_mat, elid);
  }

  /* destroy */
  sc_dmatrix_destroy (rot_el_mat);
  sc_dmatrix_destroy (weak_el_mat);
  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);
  RHEA_FREE (poly2_coeff);
}

/* Computes shear viscosity and rotation angle.*/
static void
subd_TI_viscosity_compute ( ymir_mesh_t *ymir_mesh,  ymir_vec_t *TI_svisc,
                            ymir_vec_t *viscosity,
                            ymir_vec_t *weakzone,
                            subd_options_t *subd_options)
{
  ymir_vec_multiply (viscosity, weakzone, TI_svisc);
}

/* setup the TI shear viscosity and tensor in stress operator */
void
subd_stokes_problem_setup_TI (ymir_mesh_t *ymir_mesh,
                               rhea_stokes_problem_t *stokes_problem,
                               subd_options_t *subd_options,
                               ymir_vec_t *coeff_TI_svisc,
                               ymir_vec_t *TI_rotate)
{
  const char         *this_fn_name = "subd_stokes_problem_setup_TI";
  ymir_stokes_op_t   *stokes_op;
  ymir_stress_op_t   *stress_op;
  ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);
  ymir_vec_t         *TI_weakzone = rhea_viscosity_new (ymir_mesh);

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* copy viscosity */
  rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

  /* compute the shear viscosity and rotation angles */
  subd_poly2_weakzone_compute (TI_weakzone, subd_options);
  subd_TI_viscosity_compute (ymir_mesh, coeff_TI_svisc, viscosity, TI_weakzone, subd_options);

  ymir_vec_scale (2.0, coeff_TI_svisc);

  subd_TI_rotation_compute (TI_rotate, TI_weakzone, subd_options);

  /* get the viscous stress operator */
  stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
  stress_op = stokes_op->stress_op;

 /* update viscous stress operator providing the anisotropic viscosity */
  ymir_stress_op_coeff_compute_TI_tensor (stress_op, coeff_TI_svisc,
                                          TI_rotate);
  /* destroy */
  rhea_viscosity_destroy (viscosity);
  rhea_viscosity_destroy (TI_weakzone);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/* setup the TI shear viscosity and tensor in stress operator */
void
subd_stokes_problem_setup_TI_manufactured (ymir_mesh_t *ymir_mesh,
                                           rhea_stokes_problem_t *stokes_problem,
                                           subd_options_t *subd_options,
                                           ymir_vec_t *coeff_TI_svisc,
                                           ymir_vec_t *TI_rotate)
{
  const char         *this_fn_name = "subd_stokes_problem_setup_TI_manufactured";
  ymir_stokes_op_t   *stokes_op;
  ymir_stress_op_t   *stress_op;
  ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);
  ymir_vec_t         *TI_weakzone = rhea_viscosity_new (ymir_mesh);
  const               subd_test_manufactured_t
                      test_type = subd_options->test_options->test_manufactured;
  double              rot, s_n_ratio = 0.2; /*eta_n=5, eta_s=1 */

  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (ymir_mesh);
  const ymir_locidx_t n_elements = ymir_mesh_get_num_elems_loc (ymir_mesh);
  sc_dmatrix_t       *weak_el_mat;
  double             *weak_el_data;
  double             *x, *y, *z, *tmp_el;
  ymir_locidx_t       elid;
  int                 nodeid;


  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* copy viscosity */
  rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

  /* compute the shear viscosity and rotation angles */
  switch (test_type) {
    case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT90:
    case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT45:
    case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT60:
    case SUBD_TEST_MANUFACTURED_POLY1_TIROT90:
      ymir_vec_set_value (TI_weakzone, s_n_ratio);
      break;

     /* eta_n=5, eta_s=s_n_ratio * eta_n * 0.5*(exp(y)+exp(z)) */
    case SUBD_TEST_MANUFACTURED_POLY1_TIROT90_VISCEXP:
      weak_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
      weak_el_data = weak_el_mat->e[0];
      x = RHEA_ALLOC (double, n_nodes_per_el);
      y = RHEA_ALLOC (double, n_nodes_per_el);
      z = RHEA_ALLOC (double, n_nodes_per_el);
      tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

      for (elid = 0; elid < n_elements; elid++) {
        ymir_mesh_get_elem_coord_gauss (x, y, z, elid, ymir_mesh, tmp_el);
        for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
          weak_el_data[nodeid] = s_n_ratio * 0.5 * (exp(y[nodeid]) + exp(z[nodeid]));
        }
        rhea_viscosity_set_elem_gauss (TI_weakzone, weak_el_mat, elid);
      }
      sc_dmatrix_destroy (weak_el_mat);
      RHEA_FREE (x);
      RHEA_FREE (y);
      RHEA_FREE (z);
      RHEA_FREE (tmp_el);
      break;

    case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT60_VISCEXP60:
      weak_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
      weak_el_data = weak_el_mat->e[0];
      x = RHEA_ALLOC (double, n_nodes_per_el);
      y = RHEA_ALLOC (double, n_nodes_per_el);
      z = RHEA_ALLOC (double, n_nodes_per_el);
      tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

      for (elid = 0; elid < n_elements; elid++) {
        ymir_mesh_get_elem_coord_gauss (x, y, z, elid, ymir_mesh, tmp_el);
        for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
          weak_el_data[nodeid] = s_n_ratio * exp(0.5 * (sqrt(3.0) * y[nodeid] + z[nodeid]));
        }
        rhea_viscosity_set_elem_gauss (TI_weakzone, weak_el_mat, elid);
      }
      sc_dmatrix_destroy (weak_el_mat);
      RHEA_FREE (x);
      RHEA_FREE (y);
      RHEA_FREE (z);
      RHEA_FREE (tmp_el);
      break;


    default: /* BC not set */
      RHEA_ABORT_NOT_REACHED ();
  }
  subd_TI_viscosity_compute (ymir_mesh, coeff_TI_svisc, viscosity, TI_weakzone, subd_options);
  ymir_vec_scale (2.0, coeff_TI_svisc);

  /* rotation angle */
  switch (test_type) {
    case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT90:
    case SUBD_TEST_MANUFACTURED_POLY1_TIROT90:
    case SUBD_TEST_MANUFACTURED_POLY1_TIROT90_VISCEXP:
      rot = 0.5 * M_PI;
      break;

    case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT45:
      rot = 0.25 * M_PI;
      break;

    case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT60:
    case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT60_VISCEXP60:
      rot = 1.0/3.0 * M_PI;
      break;

    default: /* BC not set */
      RHEA_ABORT_NOT_REACHED ();
  }
  ymir_vec_set_value (TI_rotate, rot);

  /* get the viscous stress operator */
  stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
  stress_op = stokes_op->stress_op;

 /* update viscous stress operator providing the anisotropic viscosity */
  ymir_stress_op_coeff_compute_TI_tensor (stress_op, coeff_TI_svisc,
                                          TI_rotate);
  /* destroy */
  rhea_viscosity_destroy (viscosity);
  rhea_viscosity_destroy (TI_weakzone);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/* compute traction as well as normal/shear stress at each element*/
void
subd_stress_TI_elem (sc_dmatrix_t * in, sc_dmatrix_t * out,
                   sc_dmatrix_t * visc, sc_dmatrix_t * svisc,
                   sc_dmatrix_t * TItens, ymir_velocity_elem_t * vel_elem,
                   double *_sc_restrict rxd, double *_sc_restrict sxd,
                   double *_sc_restrict txd, double *_sc_restrict ryd,
                   double *_sc_restrict syd, double *_sc_restrict tyd,
                   double *_sc_restrict rzd, double *_sc_restrict szd,
                   double *_sc_restrict tzd, sc_dmatrix_t * drst, sc_dmatrix_t * brst)
{
  const int           N = ymir_n (vel_elem->N);
  const int           Np = ymir_n (vel_elem->Np);
  int                 gp, i, j, k, l;
  double              S[9], E[9];
  double              cs, ct;
  double             *_sc_restrict viscd  = visc->e[0];
  double             *_sc_restrict sviscd = svisc->e[0];

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
    double               *_sc_restrict T   = TItens->e[0] + 9 * gp;
    double               *_sc_restrict outd   = out->e[0] + 9 * gp;

    S[0] = temptensd[0];
    S[1] = (temptensd[1] + temptensd[3]) * (1. / 2.);
    S[2] = (temptensd[2] + temptensd[6]) * (1. / 2.);
    S[3] = S[1];
    S[4] = temptensd[4];
    S[5] = (temptensd[5] + temptensd[7]) * (1. / 2.);
    S[6] = S[2];
    S[7] = S[5];
    S[8] = temptensd[8];

    for (k = 0; k < 9; k++) {
      E[k] = 0.0;
    }
    for (k=0; k<3; k++) {
      for (l=0; l<3; l++) {
        E[0] += T[    k] * S[3*k + l] * T[3*l    ];
        E[1] += T[    k] * S[3*k + l] * T[3*l + 1];
        E[2] += T[    k] * S[3*k + l] * T[3*l + 2];
        E[4] += T[3 + k] * S[3*k + l] * T[3*l + 1];
        E[5] += T[3 + k] * S[3*k + l] * T[3*l + 2];
        E[8] += T[6 + k] * S[3*k + l] * T[3*l + 2];
      }
    }
    E[3] = E[1];
    E[6] = E[2];
    E[7] = E[5];

    cs = viscd[gp] + sviscd[gp];
    ct = viscd[gp] - sviscd[gp];

    /* compute linear combination of viscous stress and 3x3 tensor */
    for (i = 0; i < 9; i++)  {
      outd[i] = cs * S[i] + ct * E[i];
    }
  }
}

/* compute traction as well as normal and shear stress*/
void
subd_stress_TI (ymir_cvec_t * vel, ymir_dvec_t * tau,
              ymir_dvec_t * visc, ymir_dvec_t *svisc,
              ymir_dvec_t * TItens, ymir_velocity_elem_t * vel_elem)
{
  ymir_mesh_t         *mesh = vel->mesh;
  ymir_locidx_t       elid;
  const int           N  = ymir_n (mesh->cnodes->N);
  const int           Np = ymir_np (mesh->cnodes->N);
  const int           K  = mesh->cnodes->K;
  sc_dmatrix_t       *elemin = sc_dmatrix_new (1, 3 * Np);
  sc_dmatrix_t       *elemout = sc_dmatrix_new (1,9 * Np);
  sc_dmatrix_t       *elemvisc = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemsvisc = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemTItens = sc_dmatrix_new (1, 9 * Np);
  const char          *this_fn_name = "subd_stress_TI";

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
    ymir_dvec_get_elem_interp (svisc, elemsvisc, YMIR_STRIDE_COMP, elid,
                              YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (TItens, elemTItens, YMIR_STRIDE_NODE, elid,
                              YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (tau, elemout, YMIR_STRIDE_NODE, elid,
                             YMIR_GAUSS_NODE, YMIR_WRITE);

    subd_stress_TI_elem (elemin, elemout,
                       elemvisc, elemsvisc, elemTItens,
                       vel_elem, rxd, sxd, txd, ryd,
                       syd, tyd, rzd, szd, tzd, mesh->drst, mesh->brst);

    ymir_dvec_set_elem_interp (tau, elemout, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);

    ymir_read_view_release (elemvisc);
    ymir_read_view_release (elemsvisc);
    ymir_read_view_release (elemTItens);
  }

  sc_dmatrix_destroy (elemin);
  sc_dmatrix_destroy (elemout);
  sc_dmatrix_destroy (elemvisc);
  sc_dmatrix_destroy (elemsvisc);
  sc_dmatrix_destroy (elemTItens);
}

/* compute traction as well as normal/shear stress at each element*/
void
subd_stressvec_TI_elem (sc_dmatrix_t * in, sc_dmatrix_t * out,
                   sc_dmatrix_t * visc, sc_dmatrix_t * svisc,
                   sc_dmatrix_t * TItens, ymir_velocity_elem_t * vel_elem,
                   double *_sc_restrict rxd, double *_sc_restrict sxd,
                   double *_sc_restrict txd, double *_sc_restrict ryd,
                   double *_sc_restrict syd, double *_sc_restrict tyd,
                   double *_sc_restrict rzd, double *_sc_restrict szd,
                   double *_sc_restrict tzd, sc_dmatrix_t * drst, sc_dmatrix_t * brst)
{
  const int           N = ymir_n (vel_elem->N);
  const int           Np = ymir_n (vel_elem->Np);
  int                 gp, i, j, k, l;
  double              S[9], E[9];
  double              cs, ct;
  double             *_sc_restrict viscd  = visc->e[0];
  double             *_sc_restrict sviscd = svisc->e[0];

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
    double               *_sc_restrict T   = TItens->e[0] + 9 * gp;
    double               *_sc_restrict outd   = out->e[0] + 3 * gp;

    S[0] = temptensd[0];
    S[1] = (temptensd[1] + temptensd[3]) * (1. / 2.);
    S[2] = (temptensd[2] + temptensd[6]) * (1. / 2.);
    S[3] = S[1];
    S[4] = temptensd[4];
    S[5] = (temptensd[5] + temptensd[7]) * (1. / 2.);
    S[6] = S[2];
    S[7] = S[5];
    S[8] = temptensd[8];

    for (k = 0; k < 9; k++) {
      E[k] = 0.0;
    }
    for (k=0; k<3; k++) {
      for (l=0; l<3; l++) {
        E[0] += T[    k] * S[3*k + l] * T[3*l    ];
        E[1] += T[    k] * S[3*k + l] * T[3*l + 1];
        E[2] += T[    k] * S[3*k + l] * T[3*l + 2];
        E[4] += T[3 + k] * S[3*k + l] * T[3*l + 1];
        E[5] += T[3 + k] * S[3*k + l] * T[3*l + 2];
        E[8] += T[6 + k] * S[3*k + l] * T[3*l + 2];
      }
    }
    E[3] = E[1];
    E[6] = E[2];
    E[7] = E[5];

    cs = viscd[gp] + sviscd[gp];
    ct = viscd[gp] - sviscd[gp];

    /* compute linear combination of viscous stress and 3x3 tensor */
    outd[0] = cs * S[4] + ct * E[4];
    outd[1] = cs * S[8] + ct * E[8];
    outd[2] = cs * S[5] + ct * E[5];
  }
}

/* compute traction as well as normal and shear stress*/
void
subd_stressvec_TI (ymir_cvec_t * vel, ymir_dvec_t * tauvec,
              ymir_dvec_t * visc, ymir_dvec_t *svisc,
              ymir_dvec_t * TItens, ymir_velocity_elem_t * vel_elem)
{
  ymir_mesh_t         *mesh = vel->mesh;
  ymir_locidx_t       elid;
  const int           N  = ymir_n (mesh->cnodes->N);
  const int           Np = ymir_np (mesh->cnodes->N);
  const int           K  = mesh->cnodes->K;
  sc_dmatrix_t       *elemin = sc_dmatrix_new (1, 3 * Np);
  sc_dmatrix_t       *elemout = sc_dmatrix_new (1,3 * Np);
  sc_dmatrix_t       *elemvisc = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemsvisc = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemTItens = sc_dmatrix_new (1, 9 * Np);
  const char          *this_fn_name = "subd_stress_TI";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  ymir_dvec_set_zero (tauvec);

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
    ymir_dvec_get_elem_interp (svisc, elemsvisc, YMIR_STRIDE_COMP, elid,
                              YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (TItens, elemTItens, YMIR_STRIDE_NODE, elid,
                              YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (tauvec, elemout, YMIR_STRIDE_NODE, elid,
                             YMIR_GAUSS_NODE, YMIR_WRITE);

    subd_stressvec_TI_elem (elemin, elemout,
                       elemvisc, elemsvisc, elemTItens,
                       vel_elem, rxd, sxd, txd, ryd,
                       syd, tyd, rzd, szd, tzd, mesh->drst, mesh->brst);

    ymir_dvec_set_elem_interp (tauvec, elemout, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);

    ymir_read_view_release (elemvisc);
    ymir_read_view_release (elemsvisc);
    ymir_read_view_release (elemTItens);
  }

  sc_dmatrix_destroy (elemin);
  sc_dmatrix_destroy (elemout);
  sc_dmatrix_destroy (elemvisc);
  sc_dmatrix_destroy (elemsvisc);
  sc_dmatrix_destroy (elemTItens);
}
/* compute traction as well as normal/shear stress at each element*/
void
subd_2inv_stress_TI_elem (sc_dmatrix_t * in, sc_dmatrix_t * out,
                            sc_dmatrix_t * visc, sc_dmatrix_t * svisc,
                            sc_dmatrix_t * TItens, ymir_velocity_elem_t * vel_elem,
                            double *_sc_restrict rxd, double *_sc_restrict sxd,
                            double *_sc_restrict txd, double *_sc_restrict ryd,
                            double *_sc_restrict syd, double *_sc_restrict tyd,
                            double *_sc_restrict rzd, double *_sc_restrict szd,
                            double *_sc_restrict tzd, sc_dmatrix_t * drst, sc_dmatrix_t * brst)
{
  const int           N = ymir_n (vel_elem->N);
  const int           Np = ymir_n (vel_elem->Np);
  int                 gp, i, j, k, l;
  double              S[9], E[9];
  double              cs, ct;
  double             *_sc_restrict viscd  = visc->e[0];
  double             *_sc_restrict sviscd = svisc->e[0];
  double             *_sc_restrict outd   = out->e[0];

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
    double               outv = .0;
    double               *_sc_restrict temptensd = temptens->e[0] + 9 * gp;
    double               *_sc_restrict T   = TItens->e[0] + 9 * gp;

    S[0] = temptensd[0];
    S[1] = (temptensd[1] + temptensd[3]) * (1. / 2.);
    S[2] = (temptensd[2] + temptensd[6]) * (1. / 2.);
    S[3] = S[1];
    S[4] = temptensd[4];
    S[5] = (temptensd[5] + temptensd[7]) * (1. / 2.);
    S[6] = S[2];
    S[7] = S[5];
    S[8] = temptensd[8];

    for (k = 0; k < 9; k++) {
      E[k] = 0.0;
    }
    for (k=0; k<3; k++) {
      for (l=0; l<3; l++) {
        E[0] += T[    k] * S[3*k + l] * T[3*l    ];
        E[1] += T[    k] * S[3*k + l] * T[3*l + 1];
        E[2] += T[    k] * S[3*k + l] * T[3*l + 2];
        E[4] += T[3 + k] * S[3*k + l] * T[3*l + 1];
        E[5] += T[3 + k] * S[3*k + l] * T[3*l + 2];
        E[8] += T[6 + k] * S[3*k + l] * T[3*l + 2];
      }
    }
    E[3] = E[1];
    E[6] = E[2];
    E[7] = E[5];

    cs = viscd[gp] + sviscd[gp];
    ct = viscd[gp] - sviscd[gp];

    /* compute linear combination of viscous stress and 3x3 tensor */
    for (i = 0; i < 9; i++)  {
      outv += SC_SQR (cs * S[i] + ct * E[i]) ;
    }

    outd[gp] = sqrt (0.5 * outv);
  }
}

/* compute traction as well as normal and shear stress*/
void
subd_2inv_stress_TI (ymir_cvec_t * vel, ymir_dvec_t * tauII,
                       ymir_dvec_t * visc, ymir_dvec_t *svisc,
                       ymir_dvec_t * TItens, ymir_velocity_elem_t * vel_elem)
{
  ymir_mesh_t         *mesh = vel->mesh;
  ymir_locidx_t       elid;
  const int           N  = ymir_n (mesh->cnodes->N);
  const int           Np = ymir_np (mesh->cnodes->N);
  const int           K  = mesh->cnodes->K;
  sc_dmatrix_t       *elemin = sc_dmatrix_new (1, 3 * Np);
  sc_dmatrix_t       *elemout = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemvisc = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemsvisc = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemTItens = sc_dmatrix_new (1, 9 * Np);
  const char          *this_fn_name = "subd_2inv_stress_TI";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  ymir_dvec_set_zero (tauII);

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
    ymir_dvec_get_elem_interp (svisc, elemsvisc, YMIR_STRIDE_COMP, elid,
                              YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (TItens, elemTItens, YMIR_STRIDE_NODE, elid,
                              YMIR_GAUSS_NODE, YMIR_READ);

    subd_2inv_stress_TI_elem (elemin, elemout,
                                elemvisc, elemsvisc, elemTItens,
                                vel_elem, rxd, sxd, txd, ryd,
                                syd, tyd, rzd, szd, tzd, mesh->drst, mesh->brst);

    ymir_dvec_set_elem_interp (tauII, elemout, YMIR_STRIDE_COMP, elid,
                             YMIR_GAUSS_NODE, YMIR_SET);

    ymir_read_view_release (elemvisc);
    ymir_read_view_release (elemsvisc);
    ymir_read_view_release (elemTItens);
  }

  sc_dmatrix_destroy (elemin);
  sc_dmatrix_destroy (elemout);
  sc_dmatrix_destroy (elemvisc);
  sc_dmatrix_destroy (elemsvisc);
  sc_dmatrix_destroy (elemTItens);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}


