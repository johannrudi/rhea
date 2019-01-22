#ifndef ADJOINT_TEST_H
#define ADJOINT_TEST_H

#include <adjoint_options.h>
#include <adjoint_physics.h>
#include <adjoint_postp.h>
#include <adjoint_TI.h>
#include <rhea_domain.h>
#include <rhea_base.h>
#include <ymir_vec_ops.h>
#include <ymir_vec_getset.h>
#include <ymir_stress_op.h>
#include <ymir_velocity_vec.h>
#include <ymir_comm.h>

void
subd_test_sincos1_vel_in_fn (double * vel, double x, double y,
                                          double z, ymir_locidx_t nodeid,
                                          void *data);
void
subd_test_sincos1_ISO_vel_out_fn (double * vel, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data);

void
subd_test_sincos1_TIrot90_vel_out_fn (double * vel, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data);

void
subd_test_sincos1_TIrot90_traction_fn (double * trac, double x, double y,
                                          double z, ymir_locidx_t nodeid,
                                          void *data);

void
subd_test_sincos1_TIrot90_stress_fn (double * stress, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data);
void
subd_test_sincos1_TIrot45_vel_out_fn (double * vel, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data);
void
subd_test_sincos1_TIrot45_stress_fn (double * stress, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data);
void
subd_test_sincos1_TIrot45_traction_fn (double * trac, double x, double y,
                                          double z, ymir_locidx_t nodeid,
                                          void *data);
void
subd_test_sincos1_TIrot60_vel_out_fn (double * vel, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data);
void
subd_test_sincos1_TIrot60_stress_fn (double * stress, double x, double y,
                                      double z, ymir_locidx_t nodeid,
                                      void *data);
void
subd_test_sincos1_TIrot60_traction_fn (double * trac, double x, double y,
                                          double z, ymir_locidx_t nodeid,
                                          void *data);
void
subd_test_sincos1_manufactured_set_velbc (double * vel, double x, double y,
                                              double z, ymir_locidx_t nodeid,
                                              void *data);
void
subd_test_sincos2_vel_in_fn (double * vel, double x, double y,
                                        double z, ymir_locidx_t nodeid,
                                        void *data);
void
subd_test_sincos2_TIrot90_vel_out_fn (double * vel, double x, double y,
                                         double z, ymir_locidx_t nodeid,
                                         void *data);
void
subd_test_poly1_vel_in_fn (double * vel, double x, double y,
                                        double z, ymir_locidx_t nodeid,
                                        void *data);
void
subd_test_poly1_TIrot90_vel_out_fn (double * vel, double x, double y,
                                  double z, ymir_locidx_t nodeid,
                                  void *data);
void
subd_test_poly1_TIrot90_viscexp_vel_out_fn (double * vel, double x, double y,
                                  double z, ymir_locidx_t nodeid,
                                  void *data);
void
subd_test_sincos1_TIrot60_viscexp60_vel_out_fn (double * vel, double x, double y,
                                                double z, ymir_locidx_t nodeid,
                                                void *data);

void
subd_test_poly1_manufactured_set_velbc (double * vel, double x, double y,
                                              double z, ymir_locidx_t nodeid,
                                              void *data);

void
subd_test_manufactured_rhs (ymir_vec_t *rhs_vel,
                                     ymir_vec_t *temperature /* unused */,
                                     void *data);

void
subd_test_manufactured_velbc_dir (ymir_vec_t * rhs_vel_nonzero_dirichlet,
                                       void * data);

void
subd_test_manufactured_compute_vel_err (double * abs_error, double * rel_error,
                                ymir_vec_t *vel_error, ymir_vec_t *vel_ref,
                                ymir_vec_t *vel_chk,
                                ymir_stress_op_t * stress_op);

void
subd_manufactured_strainvec_coupling_compute (ymir_vec_t *vel, ymir_vec_t *tracn, ymir_vec_t *tracs,
                                       ymir_vec_t *svisc, ymir_vec_t *nvisc);


void
subd_manufactured_stressvec_coupling_node (double *shear, double *normal,
                                            double *tau);

void
subd_manufactured_stressvec_coupling_compute (ymir_cvec_t *vel, ymir_dvec_t *visc,
                                       ymir_dvec_t *svisc, ymir_dvec_t *TItens,
                                       ymir_dvec_t *normal, ymir_dvec_t *shear,
                                       ymir_velocity_elem_t *vel_elem,
                                       subd_options_t *subd_options);

#endif /*SUBDUCTION_TEST_H*/
