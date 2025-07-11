#include <subduction_test.h>

/***********************************************************
 * Test using manufactured solution
 ***********************************************************/
/*This flow field is divergence free
 * use in both stress op test and manufactured solution test
 * for TI case with 90 degree rotation the -grad(u) is the same with that for ISO case */

void
subd_test_sincos1_vel_in_fn (double * vel, double x, double y,
                                          double z, ymir_locidx_t nodeid,
                                          void *data)
{
  vel[0] = 0.0;
  vel[1] = + sin (M_PI * y) * cos (M_PI * z);
  vel[2] = - cos (M_PI * y) * sin (M_PI * z);
}

void
subd_test_sincos1_ISO_vel_out_fn (double * vel, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data)
{
  vel[0] = 0.0;
  vel[1] = +10.0 * M_PI * M_PI * sin (M_PI * y) * cos (M_PI * z);
  vel[2] = -10.0 * M_PI * M_PI * cos (M_PI * y) * sin (M_PI * z);
}

void
subd_test_sincos1_TIrot90_vel_out_fn (double * vel, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data)
{
  vel[0] = 0.0;
  vel[1] = +10.0 * M_PI * M_PI * sin (M_PI * y) * cos (M_PI * z);
  vel[2] = -10.0 * M_PI * M_PI * cos (M_PI * y) * sin (M_PI * z);
}


/* 2eta_s*edots, 2eta_n*edotn, along the 'fault plane'*/
void
subd_test_sincos1_TIrot90_traction_fn (double * trac, double x, double y,
                                          double z, ymir_locidx_t nodeid,
                                          void *data)
{
  double tt = 90.0 * M_PI / 180.0;
  double zz = - cos(tt) * cos(tt) + sin(tt) * sin(tt);
  double yz = 2 * sin(tt) * cos(tt);
  double nvisc = 5.0, svisc = 1.0;

  /*fn=trac1=10, fs=trac2=0 */
  trac[0] = .0;
  trac[1] = 2.0 * nvisc * zz;
  trac[2] = 2.0 * svisc * yz;
  trac[1] *= ( M_PI * cos(M_PI * y) * cos(M_PI * z) );
  trac[2] *= ( M_PI * cos(M_PI * y) * cos(M_PI * z) );
}

void
subd_test_sincos1_TIrot90_stress_fn (double * stress, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data)
{
  double tt = 90.0 * M_PI / 180.0;
  double a = - cos(tt) * cos(tt) + sin(tt) * sin(tt);
  double b = 2 * sin(tt) * cos(tt);
  double nvisc = 5.0, svisc = 1.0;

  stress[0] = nvisc * a * a + svisc * b * b;
  stress[1] = -stress[0];
  stress[2] = (nvisc - svisc) * a * b;
  stress[0] *= ( 2.0 * M_PI * cos(M_PI * y) * cos(M_PI * z) );
  stress[1] *= ( 2.0 * M_PI * cos(M_PI * y) * cos(M_PI * z) );
  stress[2] *= ( 2.0 * M_PI * cos(M_PI * y) * cos(M_PI * z) );
}

void
subd_test_sincos1_TIrot45_vel_out_fn (double * vel, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data)
{
  vel[0] = 0.0;
  vel[1] = sin (M_PI * y) * cos (M_PI * z);
  vel[2] = - cos (M_PI * y) * sin (M_PI * z);
  vel[1] *= 2 * M_PI * M_PI;
  vel[2] *= 2 * M_PI * M_PI;
}

void
subd_test_sincos1_TIrot45_stress_fn (double * stress, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data)
{
  double tt = 45.0 * M_PI / 180.0;
  double a = - cos(tt) * cos(tt) + sin(tt) * sin(tt);
  double b = 2 * sin(tt) * cos(tt);
  double nvisc = 5.0, svisc = 1.0;

  stress[0] = nvisc * a * a + svisc * b * b;
  stress[1] = -stress[0];
  stress[2] = (nvisc - svisc) * a * b;
  stress[0] *= 2.0 * M_PI * cos(M_PI * y) * cos(M_PI * z);
  stress[1] *= 2.0 * M_PI * cos(M_PI * y) * cos(M_PI * z);
  stress[2] *= 2.0 * M_PI * cos(M_PI * y) * cos(M_PI * z);
}

/* 2eta_s*edots, 2eta_n*edotn, along the 'fault plane'*/
void
subd_test_sincos1_TIrot45_traction_fn (double * trac, double x, double y,
                                          double z, ymir_locidx_t nodeid,
                                          void *data)
{
  double tt = 45.0 * M_PI / 180.0;
  double zz = - cos(tt) * cos(tt) + sin(tt) * sin(tt);
  double yz = 2 * sin(tt) * cos(tt);
  double nvisc = 5.0, svisc = 1.0;

  /*fn=trac1=0, fs=trac2=2 */
  trac[0] = .0;
  trac[1] = 2.0 * nvisc * zz;
  trac[2] = 2.0 * svisc * yz;
  trac[1] *= ( M_PI * cos(M_PI * y) * cos(M_PI * z) );
  trac[2] *= ( M_PI * cos(M_PI * y) * cos(M_PI * z) );
}

void
subd_test_sincos1_TIrot60_vel_out_fn (double * vel, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data)
{
  vel[0] = 0.0;
  vel[1] = 4.0 * sin (M_PI * y) * cos (M_PI * z) - 2.0 * sqrt(3.0) * cos (M_PI * y) * sin (M_PI * z);
  vel[2] = -2.0 * sqrt(3.0) * sin (M_PI * y) * cos (M_PI * z) - 4.0 * cos (M_PI * y) * sin (M_PI * z);
  vel[1] *= (M_PI * M_PI);
  vel[2] *= (M_PI * M_PI);
}

void
subd_test_sincos1_TIrot60_stress_fn (double * stress, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data)
{
  double tt = 60.0 * M_PI / 180.0;
  double a = - cos(tt) * cos(tt) + sin(tt) * sin(tt);
  double b = 2 * sin(tt) * cos(tt);
  double nvisc = 5.0, svisc = 1.0;

  stress[0] = nvisc * a * a + svisc * b * b;
  stress[1] = -stress[0];
  stress[2] = (nvisc - svisc) * a * b;
  stress[0] *= 2.0 * M_PI * cos(M_PI * y) * cos(M_PI * z);
  stress[1] *= 2.0 * M_PI * cos(M_PI * y) * cos(M_PI * z);
  stress[2] *= 2.0 * M_PI * cos(M_PI * y) * cos(M_PI * z);
}

/* 2eta_s*edots, 2eta_n*edotn, along the 'fault plane'*/
void
subd_test_sincos1_TIrot60_traction_fn (double * trac, double x, double y,
                                          double z, ymir_locidx_t nodeid,
                                          void *data)
{
  double tt = 60.0 * M_PI / 180.0;
  double zz = - cos(tt) * cos(tt) + sin(tt) * sin(tt);
  double yz = 2 * sin(tt) * cos(tt);
  double nvisc = 5.0, svisc = 1.0;

  /*fn=trac1=5, fs=trac2=sqrt(3) */
  trac[0] = .0;
  trac[1] = 2.0 * nvisc * zz;
  trac[2] = 2.0 * svisc * yz;
  trac[1] *= ( M_PI * cos(M_PI * y) * cos(M_PI * z) );
  trac[2] *= ( M_PI * cos(M_PI * y) * cos(M_PI * z) );
}

void
subd_test_sincos1_manufactured_set_velbc (double * vel, double x, double y,
                                              double z, ymir_locidx_t nodeid,
                                              void *data)
{
  subd_options_t  *subd_options = data;
  const double  z_max = subd_options->domain_options->z_max;
  const double  y_max = subd_options->domain_options->y_max;
  const double  x_max = subd_options->domain_options->x_max;

  if (y < SC_1000_EPS || (y_max - y) < SC_1000_EPS ||
      z < SC_1000_EPS || (z_max - z) < SC_1000_EPS ||
      x < SC_1000_EPS || (x_max - x) < SC_1000_EPS)  {
    vel[0] = 0.0;
    vel[1] = + sin (M_PI * y) * cos (M_PI * z);
    vel[2] = - cos (M_PI * y) * sin (M_PI * z);
  }
  else {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }
}

/*This velocity field is not divergence-free,
 * currently only used in stress op test, not manufactured solution test*/
void
subd_test_sincos2_vel_in_fn (double * vel, double x, double y,
                                        double z, ymir_locidx_t nodeid,
                                        void *data)
{
  vel[0] = 0.0;
  vel[1] = sin (M_PI * y) * cos (M_PI * z);
  vel[2] = cos (M_PI * y) * sin (M_PI * z);
}

void
subd_test_sincos2_TIrot90_vel_out_fn (double * vel, double x, double y,
                                         double z, ymir_locidx_t nodeid,
                                         void *data)
{
  vel[0] = 0.0;
  vel[1] = 12.0 * M_PI * M_PI * sin (M_PI * y) * cos (M_PI * z);
  vel[2] = 12.0 * M_PI * M_PI * cos (M_PI * y) * sin (M_PI * z);
}

/*This flow field is from Worthen et al., 2014, PEPI, divergence-free*/
void
subd_test_poly1_vel_in_fn (double * vel, double x, double y,
                                        double z, ymir_locidx_t nodeid,
                                        void *data)
{
  vel[0] = 0.0;
  vel[1] =  y + y*y - 2.0*y*z + y*y*y - 3.0*y*z*z + y*y*z;
  vel[2] = -z - 2*y*z + z*z - 3.0*y*y*z + z*z*z - y*z*z;
}

void
subd_test_poly1_TIrot90_vel_out_fn (double * vel, double x, double y,
                                  double z, ymir_locidx_t nodeid,
                                  void *data)
{
  vel[0] = 0.0;
  vel[1] = - (18.0 + 48.0 * y + 18.0 * z);
  vel[2] = - (18.0 - 18.0 * y + 48.0 * z);
}

void
subd_test_poly1_TIrot90_viscexp_vel_out_fn (double * vel, double x, double y,
                                  double z, ymir_locidx_t nodeid,
                                  void *data)
{
  vel[0] = 0.0;
  vel[1] = - (20.0 + 60.0 * y + 20.0 * z + (exp(y) + exp(z)) * (-z - 6 * y - 1.0));
  vel[2] = - (20.0 - 20.0 * y + 60.0 * z + (exp(y) + exp(z)) * ( y - 6 * z - 1.0));
}

void
subd_test_sincos1_TIrot60_viscexp60_vel_out_fn (double * vel, double x, double y,
                                                double z, ymir_locidx_t nodeid,
                                                void *data)
{
  double visc = exp(0.5 * (sqrt(3.0) * y + z));
  double a = (2.5 + 1.5 * visc) * M_PI * M_PI;
  double b = 0.5 * sqrt(3.0) * (visc - 5.0) * M_PI * M_PI;

  vel[0] = 0.0;
  vel[1] = a * sin(M_PI * y) * cos(M_PI * z) + b * cos(M_PI * y) * sin(M_PI * z);
  vel[2] = b * sin(M_PI * y) * cos(M_PI * z) - a * cos(M_PI * y) * sin(M_PI * z);
}

void
subd_test_poly1_manufactured_set_velbc (double * vel, double x, double y,
                                              double z, ymir_locidx_t nodeid,
                                              void *data)
{
  subd_options_t  *subd_options = data;
  const double  z_max = subd_options->domain_options->z_max;
  const double  y_max = subd_options->domain_options->y_max;
  const double  x_max = subd_options->domain_options->x_max;

  if (y < SC_1000_EPS || (y_max - y) < SC_1000_EPS ||
      z < SC_1000_EPS || (z_max - z) < SC_1000_EPS ||
      x < SC_1000_EPS || (x_max - x) < SC_1000_EPS) {
    vel[0] = 0.0;
    vel[1] =  y + y*y - 2.0*y*z + y*y*y - 3.0*y*z*z + y*y*z;
    vel[2] = -z - 2*y*z + z*z - 3.0*y*y*z + z*z*z - y*z*z;
  }
  else {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }
}

/* manufactured solution */
void
subd_test_manufactured_rhs (ymir_vec_t *rhs_vel,
                                     ymir_vec_t *temperature /* unused */,
                                     void *data)
{
  const char          *this_fn_name = "subd_test_manufactured_rhs_compute";
  subd_options_t     *subd_options = data;
  const               subd_test_manufactured_t
                      test_type = subd_options->test_options->test_manufactured;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);
  /* check input */
  RHEA_ASSERT (test_type != SUBD_TEST_MANUFACTURED_NONE);

  /* compute velocity fields */
  switch (test_type) {
    case SUBD_TEST_MANUFACTURED_SINCOS1_ISO:
      /* compute reference velocity field (output) */
      ymir_cvec_set_function (rhs_vel, subd_test_sincos1_ISO_vel_out_fn,
                              NULL);
      break;

    case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT90:
      /* compute reference velocity field (output) */
      ymir_cvec_set_function (rhs_vel, subd_test_sincos1_TIrot90_vel_out_fn,
                              NULL);
      break;

    case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT45:
      /* compute reference velocity field (output) */
      ymir_cvec_set_function (rhs_vel, subd_test_sincos1_TIrot45_vel_out_fn,
                              NULL);
      break;

    case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT60:
      /* compute reference velocity field (output) */
      ymir_cvec_set_function (rhs_vel, subd_test_sincos1_TIrot60_vel_out_fn,
                              NULL);
      break;

    case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT60_VISCEXP60:
      /* compute reference velocity field (output) */
      ymir_cvec_set_function (rhs_vel, subd_test_sincos1_TIrot60_viscexp60_vel_out_fn,
                              NULL);
      break;


    case SUBD_TEST_MANUFACTURED_POLY1_TIROT90:
      /* compute reference velocity field (output) */
      ymir_cvec_set_function (rhs_vel, subd_test_poly1_TIrot90_vel_out_fn,
                              NULL);
      break;

    case SUBD_TEST_MANUFACTURED_POLY1_TIROT90_VISCEXP:
      /* compute reference velocity field (output) */
      ymir_cvec_set_function (rhs_vel, subd_test_poly1_TIrot90_viscexp_vel_out_fn,
                              NULL);
      break;

  default:
      RHEA_ABORT_NOT_REACHED ();
  }
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

void
subd_test_manufactured_velbc_dir (ymir_vec_t * rhs_vel_nonzero_dirichlet,
                                       void * data)
{
  const char          *this_fn_name = "subd_test_manufactured_velbc_compute";
  subd_options_t     *subd_options = data;
  const               subd_test_manufactured_t
                      test_type = subd_options->test_options->test_manufactured;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);
  switch (test_type) {
    case SUBD_TEST_MANUFACTURED_SINCOS1_ISO:
    case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT90:
    case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT45:
    case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT60:
    case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT60_VISCEXP60:
      ymir_cvec_set_function (rhs_vel_nonzero_dirichlet,
                              subd_test_sincos1_manufactured_set_velbc,
                              subd_options);
      break;

    case SUBD_TEST_MANUFACTURED_POLY1_TIROT90:
    case SUBD_TEST_MANUFACTURED_POLY1_TIROT90_VISCEXP:
      ymir_cvec_set_function (rhs_vel_nonzero_dirichlet,
                              subd_test_poly1_manufactured_set_velbc,
                              subd_options);
      break;

   default: /* BC not set */
      RHEA_ABORT_NOT_REACHED ();
  }
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

void
subd_test_manufactured_compute_vel_err (double * abs_error, double * rel_error,
                                ymir_vec_t *vel_error, ymir_vec_t *vel_ref,
                                ymir_vec_t *vel_chk,
                                ymir_stress_op_t * stress_op)
{
  ymir_vec_t         *vel_ref_zero_bndr = ymir_vec_clone (vel_ref);

  /* calculate error vector */
  ymir_vec_copy (vel_chk, vel_error);
  ymir_vec_add (-1.0, vel_ref, vel_error);

  /* set boundary values to zero */
  if (ymir_stress_op_has_dirichlet (stress_op)) {
    ymir_vel_dir_separate (vel_error, NULL, NULL, NULL, stress_op->vel_dir);
    ymir_vel_dir_separate (vel_ref_zero_bndr, NULL, NULL, NULL,
                           stress_op->vel_dir);
  }

  /* take norm of error vector to get a value for absolute and relative errors */
  *abs_error = ymir_vec_norm (vel_error);
  *rel_error = *abs_error / ymir_vec_norm (vel_ref_zero_bndr);

  /* destroy */
  ymir_vec_destroy (vel_ref_zero_bndr);
}

void
subd_manufactured_stressvec_coupling_node (double *shear, double *normal,
                                            double *tau)
{
  double              orth[3], tangent[3], tempvec[3];
  double              tempt, tt;

  orth[0] = tangent[0] = 0.0;

   tt = 60.0 * M_PI / 180.0;
   orth[1] = sin(tt);
   orth[2] = cos(tt);
   tangent[1] = cos(tt);
   tangent[2] = -sin(tt);

   tempvec[0] = ( tau[0] * orth[0] + tau[1] * orth[1] + tau[2] * orth[2]);
   tempvec[1] = ( tau[3] * orth[0] + tau[4] * orth[1] + tau[5] * orth[2]);
   tempvec[2] = ( tau[6] * orth[0] + tau[7] * orth[1] + tau[8] * orth[2]);

   * shear  = (tempvec[0] * tangent[0] + tempvec[1] * tangent[1] + tempvec[2] * tangent[2]);
   * normal = (tempvec[0] * orth[0] + tempvec[1] * orth[1] + tempvec[2] * orth[2]);
}



/* Computes the shear and normal traction along 2plates_poly2 weakzone in 2D Cartesian domain*/
void
subd_manufactured_strainvec_coupling_node (double *tracn, double *tracs, double *gradv,
                                           double *svisc, double *nvisc)
{
  double              orth[3], tangent[3], tempvec[3];
  double              tempt, tt;

  orth[0] = tangent[0] = 0.0;

  tt = 60.0 * M_PI / 180.0;
  orth[1] = sin(tt);
  orth[2] = cos(tt);
  tangent[1] = cos(tt);
  tangent[2] = -sin(tt);

  tempvec[0] = ( gradv[0] * orth[0] + gradv[1] * orth[1] + gradv[2] * orth[2]);
  tempvec[1] = ( gradv[1] * orth[0] + gradv[3] * orth[1] + gradv[4] * orth[2]);
  tempvec[2] = ( gradv[2] * orth[0] + gradv[4] * orth[1] + gradv[5] * orth[2]);

  tracn[0] = 2.0 * nvisc[0] *
            (tempvec[0] * orth[0] + tempvec[1] * orth[1] + tempvec[2] * orth[2]);
  tracs[0] = 2.0 * svisc[0] *
            (tempvec[0] * tangent[0] + tempvec[1] * tangent[1] + tempvec[2] * tangent[2]);

}

void
subd_manufactured_strainvec_coupling_compute (ymir_vec_t *vel, ymir_vec_t *tracn, ymir_vec_t *tracs,
                                       ymir_vec_t *svisc, ymir_vec_t *nvisc)
{
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (vel);
  const ymir_locidx_t n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  ymir_dvec_t        *gradvel = ymir_dvec_new (mesh, 6, YMIR_GAUSS_NODE);
  sc_dmatrix_t       *elemtracn, *elemtracs, *elemgradv, *elemsvisc, *elemnvisc;
  ymir_locidx_t       elid;
  int                 nodeid;

 const char         *this_fn_name = "subd_manufactured_strainvec_coupling_compute";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  ymir_velocity_strain_rate (vel, gradvel, 0);

  /* create work variables */
  elemtracn = sc_dmatrix_new (1, n_nodes_per_el);
  elemtracs = sc_dmatrix_new (1, n_nodes_per_el);
  elemgradv = sc_dmatrix_new (1, 6 * n_nodes_per_el);
  elemsvisc = sc_dmatrix_new (1, n_nodes_per_el);
  elemnvisc = sc_dmatrix_new (1, n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    ymir_dvec_get_elem_interp (gradvel, elemgradv, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (svisc, elemsvisc, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (nvisc, elemnvisc, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_READ);

    ymir_dvec_get_elem_interp (tracn, elemtracn, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_WRITE);
    ymir_dvec_get_elem_interp (tracs, elemtracs, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_WRITE);

  /* compute weak zone factor for each node */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    double               *_sc_restrict tracn_data = elemtracn->e[0] + nodeid;
    double               *_sc_restrict tracs_data = elemtracs->e[0] + nodeid;
    double               *_sc_restrict gradv_data = elemgradv->e[0] + 6 * nodeid;
    double               *_sc_restrict svisc_data = elemsvisc->e[0] + nodeid;
    double               *_sc_restrict nvisc_data = elemnvisc->e[0] + nodeid;
    subd_manufactured_strainvec_coupling_node (tracn_data, tracs_data, gradv_data,
                                                svisc_data, nvisc_data);
  }

    /* set traction of this element */
    ymir_dvec_set_elem_interp (tracn, elemtracn, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);
    ymir_dvec_set_elem_interp (tracs, elemtracs, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);

    ymir_read_view_release (elemgradv);
    ymir_read_view_release (elemsvisc);
    ymir_read_view_release (elemnvisc);
  }

  /* destroy */
  sc_dmatrix_destroy (elemtracn);
  sc_dmatrix_destroy (elemtracs);
  sc_dmatrix_destroy (elemgradv);
  sc_dmatrix_destroy (elemsvisc);
  sc_dmatrix_destroy (elemnvisc);
  ymir_vec_destroy (gradvel);
}

void
subd_manufactured_stressvec_coupling_compute (ymir_cvec_t *vel, ymir_dvec_t *visc,
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
  int                 nodeid;

  const char          *this_fn_name = "subd_manufactured_stressvec_coupling_compute";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  if (subd_options->visc_options->anisotropy_type
      == SUBD_VISC_TRANSVERSELY_ISOTROPY)  {
    subd_stress_TI (vel, tau, visc, svisc, TItens, vel_elem);
  }
  else  {
    subd_stress (vel, tau, visc, vel_elem);
  }

 /* create work variables */
  elemtau = sc_dmatrix_new (1, 9 * n_nodes_per_el);
  elemout_s = sc_dmatrix_new (1, n_nodes_per_el);
  elemout_n = sc_dmatrix_new (1, n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    ymir_dvec_get_elem_interp (tau, elemtau, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (shear, elemout_s, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_WRITE);
    ymir_dvec_get_elem_interp (normal, elemout_n, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_WRITE);

    /* compute weak zone factor for each node */
    for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
      double               *_sc_restrict tau_data = elemtau->e[0] + 9 * nodeid;
      double               *_sc_restrict shear_data = elemout_s->e[0] + nodeid;
      double               *_sc_restrict normal_data = elemout_n->e[0] + nodeid;
      subd_manufactured_stressvec_coupling_node (shear_data, normal_data, tau_data);
    }

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
  ymir_vec_destroy (tau);
}


