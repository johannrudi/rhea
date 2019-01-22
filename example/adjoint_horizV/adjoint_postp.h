#ifndef ADJOINT_POSTP_H
#define ADJOINT_POSTP_H

#include <adjoint_options.h>
#include <adjoint_physics.h>
#include <adjoint_TI.h>
#include <rhea_domain.h>
#include <rhea_base.h>
#include <ymir_vec_ops.h>
#include <ymir_vec_getset.h>
#include <ymir_stress_op.h>
#include <ymir_velocity_vec.h>
#include <ymir_comm.h>
#include <ymir_stokes_vec.h>
#include <ymir_derivative_elem.h>
#include <ymir_stokes_pc.h>
#include <ymir_interp_vec.h>
#include <ymir_mass_vec.h>
#include <ymir_vtk.h>

/* compute traction as well as normal/shear stress at each element*/
void
subd_stress_elem (sc_dmatrix_t * in, sc_dmatrix_t * out,
                   sc_dmatrix_t * visc, ymir_velocity_elem_t * vel_elem,
                   double *_sc_restrict rxd, double *_sc_restrict sxd,
                   double *_sc_restrict txd, double *_sc_restrict ryd,
                   double *_sc_restrict syd, double *_sc_restrict tyd,
                   double *_sc_restrict rzd, double *_sc_restrict szd,
                   double *_sc_restrict tzd, sc_dmatrix_t * drst, sc_dmatrix_t * brst);

/* compute traction as well as normal and shear stress*/
void
subd_stress (ymir_cvec_t * vel, ymir_dvec_t * tau,
              ymir_dvec_t * visc, ymir_velocity_elem_t * vel_elem);


/* Computes the shear and normal traction along 2plates_poly2 weakzone in 2D Cartesian domain*/
void
subd_postp_weakzone_coupling_brick_2plates_poly2 (double r, double lon,
                                                   double *shear, double *normal,
                                                   double *tau,
                                                   subd_options_t * subd_options);

void
subd_postp_weakzone_coupling_node (const double x, const double y, const double z,
                                    double *shear, double *normal,
                                    double *tau,
                                    subd_options_t *subd_options);

/* compute weak zone factor of an element */
void
subd_postp_weakzone_coupling_elem (sc_dmatrix_t * tau,
                                    sc_dmatrix_t * shear,
                                    sc_dmatrix_t * normal,
                                    const double *x,
                                    const double *y,
                                    const double *z,
                                    const int n_nodes_per_el,
                                    subd_options_t *subd_options);

void
subd_postp_weakzone_coupling_compute (ymir_cvec_t *vel, ymir_dvec_t *visc,
                                       ymir_dvec_t *svisc, ymir_dvec_t *TItens,
                                       ymir_dvec_t *normal, ymir_dvec_t *shear,
                                       ymir_velocity_elem_t *vel_elem,
                                       subd_options_t *subd_options);

/* Computes the shear and normal traction along 2plates_poly2 weakzone in 2D Cartesian domain*/
void
subd_postp_weakzone_traction_brick_2plates_poly2 (double r, double lon,
                                                   double *trac, double *gradv,
                                                   double *svisc, double *nvisc,
                                                   subd_options_t * subd_options);

void
subd_postp_weakzone_traction_node (const double x, const double y, const double z,
                                    double *trac, double *gradv,
                                    double *svisc, double *nvisc,
                                    subd_options_t *subd_options);

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
                                    subd_options_t *subd_options);

void
subd_postp_weakzone_traction_compute (ymir_vec_t *vel, ymir_vec_t *traction,
                                       ymir_vec_t *svisc, ymir_vec_t *nvisc,
                                       subd_options_t *subd_options);

/* compute traction. It is an alternative approach that takes advantage of an existing
   subroutine ymir_velocity_strain_rate and directly compute traction on each node*/
void
subd_traction (ymir_cvec_t * vel, ymir_dvec_t *traction,
                  double *n_dir, ymir_dvec_t *visc);

/* compute traction as well as normal/shear stress at each element*/
void
subd_normal_stress_elem (sc_dmatrix_t * in, sc_dmatrix_t * out1, sc_dmatrix_t * out2,
                            sc_dmatrix_t * out3, double * n_dir,
                  sc_dmatrix_t * visc_mat, ymir_velocity_elem_t * vel_elem,
                  double *_sc_restrict rxd, double *_sc_restrict sxd,
                  double *_sc_restrict txd, double *_sc_restrict ryd,
                  double *_sc_restrict syd, double *_sc_restrict tyd,
                  double *_sc_restrict rzd, double *_sc_restrict szd,
                  double *_sc_restrict tzd, sc_dmatrix_t * drst, sc_dmatrix_t * brst);

/* compute traction as well as normal and shear stress*/
void
subd_normal_stress (ymir_cvec_t * vel, ymir_dvec_t * n_tau, ymir_dvec_t * s_tau,
                        ymir_dvec_t * traction, double * n_dir,
                       ymir_dvec_t * visc, ymir_velocity_elem_t * vel_elem);

/* normal stress on the surface */
void
subd_physics_compute_normal_boundary_stress (ymir_vec_t *stress_bndr_norm,
                                              ymir_vec_t *up,
                                              ymir_vec_t *rhs_u_point,
                                              ymir_stokes_op_t *stokes_op);

void
subd_postp_topography (ymir_vec_t *topography,
                      rhea_stokes_problem_t *stokes_problem,
                      subd_options_t *subd_options);


#endif
