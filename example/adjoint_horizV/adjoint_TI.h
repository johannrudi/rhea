#ifndef ADJOINT_TI_H
#define ADJOINT_TI_H

#include <adjoint_options.h>
#include <adjoint_physics.h>
#include <rhea.h>
#include <ymir_vec_getset.h>
#include <ymir_vec_ops.h>
#include <ymir_comm.h>
#include <ymir_stress_op.h>
#include <ymir_velocity_vec.h>
#include <ymir_derivative_elem.h>

void
subd_stress_op_copy_shear_visc (ymir_vec_t *shear_visc,
                                   ymir_stress_op_t *stress_op);

void
subd_stress_op_copy_TI_tensor (ymir_vec_t *TI_tensor,
                                   ymir_stress_op_t *stress_op);


static double
subd_TI_rotation_brick_2plates_poly2 (double r, double lon,
                                       subd_options_t * subd_options);


static double
subd_TI_rotation_node (const double x, const double y, const double z,
                        subd_options_t *subd_options);


void
subd_TI_rotation_elem (double *_sc_restrict rot_elem,
                        double *_sc_restrict weak_elem,
                         const double *x,
                         const double *y,
                         const double *z,
                         const int n_nodes_per_el,
                         subd_options_t *subd_options);


static void
subd_TI_rotation_compute (ymir_vec_t *rotate, ymir_vec_t *weak,
                           subd_options_t *subd_options);


/* Computes shear viscosity and rotation angle.*/
static void
subd_TI_viscosity_compute ( ymir_mesh_t *ymir_mesh,  ymir_vec_t *TI_svisc,
                            ymir_vec_t *viscosity,
                            ymir_vec_t *weakzone,
                            subd_options_t *subd_options);

/* setup the TI shear viscosity and tensor in stress operator */
void
subd_stokes_problem_setup_TI (ymir_mesh_t *ymir_mesh,
                               rhea_stokes_problem_t *stokes_problem,
                               subd_options_t *subd_options,
                               ymir_vec_t *coeff_TI_svisc,
                               ymir_vec_t *TI_rotate);

/* setup the TI shear viscosity and tensor in stress operator */
void
subd_stokes_problem_setup_TI_manufactured (ymir_mesh_t *ymir_mesh,
                                           rhea_stokes_problem_t *stokes_problem,
                                           subd_options_t *subd_options,
                                           ymir_vec_t *coeff_TI_svisc,
                                           ymir_vec_t *TI_rotate);

/* compute traction as well as normal/shear stress at each element*/
void
subd_stress_TI_elem (sc_dmatrix_t * in, sc_dmatrix_t * out,
                   sc_dmatrix_t * visc, sc_dmatrix_t * svisc,
                   sc_dmatrix_t * TItens, ymir_velocity_elem_t * vel_elem,
                   double *_sc_restrict rxd, double *_sc_restrict sxd,
                   double *_sc_restrict txd, double *_sc_restrict ryd,
                   double *_sc_restrict syd, double *_sc_restrict tyd,
                   double *_sc_restrict rzd, double *_sc_restrict szd,
                   double *_sc_restrict tzd, sc_dmatrix_t * drst, sc_dmatrix_t * brst);

/* compute traction as well as normal and shear stress*/
void
subd_stress_TI (ymir_cvec_t * vel, ymir_dvec_t * tau,
              ymir_dvec_t * visc, ymir_dvec_t *svisc,
              ymir_dvec_t * TItens, ymir_velocity_elem_t * vel_elem);

/* compute traction as well as normal/shear stress at each element*/
void
subd_stressvec_TI_elem (sc_dmatrix_t * in, sc_dmatrix_t * out,
                   sc_dmatrix_t * visc, sc_dmatrix_t * svisc,
                   sc_dmatrix_t * TItens, ymir_velocity_elem_t * vel_elem,
                   double *_sc_restrict rxd, double *_sc_restrict sxd,
                   double *_sc_restrict txd, double *_sc_restrict ryd,
                   double *_sc_restrict syd, double *_sc_restrict tyd,
                   double *_sc_restrict rzd, double *_sc_restrict szd,
                   double *_sc_restrict tzd, sc_dmatrix_t * drst, sc_dmatrix_t * brst);

/* compute traction as well as normal and shear stress*/
void
subd_stressvec_TI (ymir_cvec_t * vel, ymir_dvec_t * tauvec,
              ymir_dvec_t * visc, ymir_dvec_t *svisc,
              ymir_dvec_t * TItens, ymir_velocity_elem_t * vel_elem);

/* compute traction as well as normal/shear stress at each element*/
void
subd_2inv_stress_TI_elem (sc_dmatrix_t * in, sc_dmatrix_t * out,
                            sc_dmatrix_t * visc, sc_dmatrix_t * svisc,
                            sc_dmatrix_t * TItens, ymir_velocity_elem_t * vel_elem,
                            double *_sc_restrict rxd, double *_sc_restrict sxd,
                            double *_sc_restrict txd, double *_sc_restrict ryd,
                            double *_sc_restrict syd, double *_sc_restrict tyd,
                            double *_sc_restrict rzd, double *_sc_restrict szd,
                            double *_sc_restrict tzd, sc_dmatrix_t * drst, sc_dmatrix_t * brst);

/* compute traction as well as normal and shear stress*/
void
subd_2inv_stress_TI (ymir_cvec_t * vel, ymir_dvec_t * tauII,
                       ymir_dvec_t * visc, ymir_dvec_t *svisc,
                       ymir_dvec_t * TItens, ymir_velocity_elem_t * vel_elem);

#endif /*SUBDUCTION_TI_H*/
