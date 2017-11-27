#include <slabweakzone.h>
#include <rhea.h>
#include <ymir_velocity_vec.h>
#include <ymir_stress_op.h>

/**
 * Compute TI_tensor and rotation angle
 */
void
slabs_stokes_problem_setup_TI_manufactured (ymir_mesh_t *ymir_mesh,
                                           rhea_stokes_problem_t *stokes_problem,
                                           slabs_options_t *slabs_options,
                                           ymir_vec_t *coeff_TI_svisc,
                                           ymir_vec_t *TI_rotate);

/*This flow field is divergence free
 * used in both stress_op test and manufactured solution test */
static void
slabs_test_sincos1_vel_in_fn (double * vel, double x, double y,
                                          double z, ymir_locidx_t nodeid,
                                          void *data);
/*-grad(u) for isotropic viscous model*/
static void
slabs_test_sincos1_ISO_vel_out_fn (double * vel, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data);

/* for TI case with 90 degree rotation the -grad(u) is the same with that for ISO case*/
static void
slabs_test_sincos1_TIrot90_vel_out_fn (double * vel, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data);

/* 2eta_s*edots, 2eta_n*edotn, along the 'fault plane'*/
static void
slabs_test_sincos1_TIrot90_traction_fn (double * trac, double x, double y,
                                          double z, ymir_locidx_t nodeid,
                                          void *data);

/* on original coordiniate, sigma_yy, sigma_zz, sigma_yz */
static void
slabs_test_sincos1_TIrot90_stress_fn (double * stress, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data);

static void
slabs_test_sincos1_TIrot45_vel_out_fn (double * vel, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data);

static void
slabs_test_sincos1_TIrot45_stress_fn (double * stress, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data);

static void
slabs_test_sincos1_TIrot45_traction_fn (double * trac, double x, double y,
                                          double z, ymir_locidx_t nodeid,
                                          void *data);

static void
slabs_test_sincos1_TIrot60_vel_out_fn (double * vel, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data);

static void
slabs_test_sincos1_TIrot60_stress_fn (double * stress, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data);

static void
slabs_test_sincos1_TIrot60_traction_fn (double * trac, double x, double y,
                                          double z, ymir_locidx_t nodeid,
                                          void *data);

static void
slabs_test_sincos1_TIrot60_viscexp60_vel_out_fn (double * vel, double x, double y,
                                                double z, ymir_locidx_t nodeid,
                                                void *data);

static void
slabs_test_sincos1_manufactured_set_velbc (double * vel, double x, double y,
                                              double z, ymir_locidx_t nodeid,
                                              void *data);

/*This velocity field is not divergence-free,
 * currently only used in stress op test, not manufactured solution test*/
static void
slabs_test_sincos2_vel_in_fn (double * vel, double x, double y,
                                        double z, ymir_locidx_t nodeid,
                                        void *data);

static void
slabs_test_sincos2_TIrot90_vel_out_fn (double * vel, double x, double y,
                                         double z, ymir_locidx_t nodeid,
                                         void *data);

/*This flow field is from Worthen et al., 2014, PEPI, divergence-free*/
static void
slabs_test_poly1_vel_in_fn (double * vel, double x, double y,
                                        double z, ymir_locidx_t nodeid,
                                        void *data);

static void
slabs_test_poly1_TIrot90_vel_out_fn (double * vel, double x, double y,
                                  double z, ymir_locidx_t nodeid,
                                  void *data);

static void
slabs_test_poly1_TIrot90_viscexp_vel_out_fn (double * vel, double x, double y,
                                  double z, ymir_locidx_t nodeid,
                                  void *data);

/*impose velocity boundary condition*/
static void
slabs_test_poly1_manufactured_set_velbc (double * vel, double x, double y,
                                              double z, ymir_locidx_t nodeid,
                                              void *data);

/*compute rhs*/
static void
slabs_test_manufactured_rhs_compute (ymir_vec_t *rhs_vel,
                                     slabs_test_options_t *slabs_test_options);

static ymir_dir_code_t
slabs_test_manufactured_set_vel_dir_all (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data);

static ymir_dir_code_t
slabs_test_manufactured_set_vel_dir_all_2D (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data);

void
slabs_test_manufactured_velbc_compute (ymir_vec_t * rhs_vel_nonzero_dirichlet,
                                       slabs_options_t * slabs_options);

static void
slabs_test_manufactured_compute_vel_err (double * abs_error, double * rel_error,
                                ymir_vec_t *vel_error, ymir_vec_t *vel_ref,
                                ymir_vec_t *vel_chk,
                                ymir_stress_op_t * stress_op);

/*post-processing, compute sigma_yy, sigma_zz, sigma_yz  */
void
slabs_manufactured_stressvec_coupling_node (double *shear, double *normal,
                                            double *tau);

static void
slabs_manufactured_stressvec_coupling_compute (ymir_cvec_t *vel, ymir_dvec_t *visc,
                                       ymir_dvec_t *svisc, ymir_dvec_t *TItens,
                                       ymir_dvec_t *normal, ymir_dvec_t *shear,
                                       ymir_velocity_elem_t *vel_elem,
                                       slabs_options_t *slabs_options);

/*post-processing, compute shear and normal stress component on fault plane*/
void
slabs_manufactured_strainvec_coupling_node (double *tracn, double *tracs, double *gradv,
                                           double *svisc, double *nvisc);


static void
slabs_manufactured_strainvec_coupling_compute (ymir_vec_t *vel, ymir_vec_t *tracn, ymir_vec_t *tracs,
                                       ymir_vec_t *svisc, ymir_vec_t *nvisc);
