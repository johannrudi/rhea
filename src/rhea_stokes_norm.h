/*
 */

#ifndef RHEA_STOKES_NORM_H
#define RHEA_STOKES_NORM_H

#include <ymir_vec_ops.h>
#include <ymir_pressure_elem.h>
#include <ymir_Hminus1_norm.h>

/* enumerator for types of Stokes norms */
typedef enum
{
  RHEA_STOKES_NORM_NONE = -1,
  RHEA_STOKES_NORM_L2_VEC_SP,  /* l2-norm for vectors */
  RHEA_STOKES_NORM_L2_PRIMAL,  /* L2-norm for functions (primal space) */
  RHEA_STOKES_NORM_L2_DUAL,    /* L2-norm for functions (dual/residual space) */
  RHEA_STOKES_NORM_HMINUS1_L2  /* mixed norm (H^-1, L2) for vel-press pairing */
}
rhea_stokes_norm_type_t;

/**
 * Computes the norm of a Stokes (i.e., velocity-pressure) vector.
 */
double              rhea_stokes_norm_compute (
                                            ymir_vec_t *vec,
                                            rhea_stokes_norm_type_t norm_type,
                                            ymir_Hminus1_norm_op_t *norm_op,
                                            ymir_pressure_elem_t *press_elem,
                                            double *norm_vel,
                                            double *norm_press);

/**
 * Computes the inner product of the individual velocity and pressure
 * components of a Stokes vector.
 */
void                rhea_stokes_norm_innerprod (
                                            double *innerprod_vel,
                                            double *innerprod_press,
                                            ymir_vec_t *arg_left,
                                            ymir_vec_t *arg_right,
                                            rhea_stokes_norm_type_t norm_type,
                                            ymir_Hminus1_norm_op_t *norm_op,
                                            ymir_pressure_elem_t *press_elem);

#endif /* RHEA_STOKES_NORM_H */
