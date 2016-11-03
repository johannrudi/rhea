/*
 */

#ifndef RHEA_VISCOSITY_H
#define RHEA_VISCOSITY_H

#include <rhea_domain.h>
#include <ymir_vec_ops.h>

/* enumerator for types of viscosities */
typedef enum
{
  RHEA_VISCOSITY_CONST,
  RHEA_VISCOSITY_LINEAR,
  RHEA_VISCOSITY_NONLINEAR
}
rhea_viscosity_t;

/* enumerator for types of initial viscosities for nonlinear Stokes problems */
typedef enum
{
  RHEA_VISCOSITY_INIT_NONLINEAR_DEFAULT,
  RHEA_VISCOSITY_INIT_NONLINEAR_CONST_CENTERED_TO_BOUNDS,
  RHEA_VISCOSITY_INIT_NONLINEAR_CONST_MIN_BOUND,
  RHEA_VISCOSITY_INIT_NONLINEAR_TEMP,
  RHEA_VISCOSITY_INIT_NONLINEAR_TEMP_UM_REL_TO_LM
}
rhea_viscosity_init_nonlinear_t;

/* enumerator for types of viscosity models */
typedef enum
{
  /* (1) upper bound, (2) weak zone, (3) yielding, (4) lower bound
   * viscosity bounds via cut-off */
  RHEA_VISCOSITY_MODEL_UWYL,

  /* (1) upper bound, (2) weak zone, (3) yielding, (4) lower bound
   * upper viscosity bound via cut-off, lower viscosity bound via addition */
  RHEA_VISCOSITY_MODEL_UWYL_LADD_UCUT,

  /* (1) upper bound + shift, (2) weak zone, (3) yielding, (4) lower bound
   * upper viscosity bound via shift, lower viscosity bound via addition */
  RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT
}
rhea_viscosity_model_t;

/* options of the mantle's viscosity */
typedef struct rhea_viscosity_options
{
  /* type & model of the viscosity */
  rhea_viscosity_t                 type;
  rhea_viscosity_init_nonlinear_t  type_init_nonlinear;
  rhea_viscosity_model_t           model;

  /* lower and upper bounds for the viscosity */
  double              min;
  double              max;

  /* scaling factor and activation energy in Arrhenius relationship */
  double              upper_mantle_scaling;
  double              upper_mantle_activation_energy;
  double              lower_mantle_scaling;
  double              lower_mantle_activation_energy;

  /* stress exponent that governs strain rate weakening (aka. `n`) */
  double              stress_exponent;

  /* value of viscous stress above which plastic yielding occurs */
  double              yield_stress;
  double              yielding_regularization;

  /* options & properties of the computational domain */
  rhea_domain_options_t  *domain_options;
}
rhea_viscosity_options_t;

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_viscosity_add_options (ymir_options_t * opt_sup);

/**
 * Processes options and stores them.
 */
void                rhea_viscosity_process_options (
                                        rhea_viscosity_options_t *opt,
                                        rhea_domain_options_t *domain_options);
/**
 * Creates a new viscosity vector.
 */
ymir_vec_t         *rhea_viscosity_new (ymir_mesh_t *ymir_mesh);

/**
 * Destroys a viscosity vector.
 */
void                rhea_viscosity_destroy (ymir_vec_t *viscosity);

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_viscosity_check_vec_type (ymir_vec_t *vec);

/**
 * Checks entries of a vector.
 */
int                 rhea_viscosity_is_valid (ymir_vec_t *vec);

/**
 * Computes viscosity.
 */
void                rhea_viscosity_compute (ymir_vec_t *viscosity,
                                            ymir_vec_t *rank1_tensor_scal,
                                            ymir_vec_t *bounds_marker,
                                            ymir_vec_t *yielding_marker,
                                            ymir_vec_t *temperature,
                                            ymir_vec_t *weakzone,
                                            ymir_vec_t *velocity,
                                            rhea_viscosity_options_t *opt);

/**
 * Gets the viscosity of one element at Gauss nodes.
 */
double             *rhea_viscosity_get_elem_gauss (sc_dmatrix_t *visc_el_mat,
                                                   ymir_vec_t *visc_vec,
                                                   const ymir_locidx_t elid);

/**
 * Sets the viscosity of one element at Gauss nodes.
 */
void                rhea_viscosity_set_elem_gauss (ymir_vec_t *visc_vec,
                                                   sc_dmatrix_t *visc_el_mat,
                                                   const ymir_locidx_t elid);

/**
 * Gets the rank-1 4th-order tensor of one element at Gauss nodes.
 */
double             *rhea_viscosity_rank1_scal_get_elem_gauss (
                                              sc_dmatrix_t *rank1_scal_el_mat,
                                              ymir_vec_t *rank1_scal_vec,
                                              const ymir_locidx_t elid);

/**
 * Sets the rank-1 4th-order tensor of one element at Gauss nodes.
 */
void                rhea_viscosity_rank1_scal_set_elem_gauss (
                                              ymir_vec_t *rank1_scal_vec,
                                              sc_dmatrix_t *rank1_scal_el_mat,
                                              const ymir_locidx_t elid);

/**
 * Gets a marker of one element at Gauss nodes.
 */
double             *rhea_viscosity_marker_get_elem_gauss (
                                                  sc_dmatrix_t *marker_el_mat,
                                                  ymir_vec_t *marker_vec,
                                                  const ymir_locidx_t elid);

/**
 * Sets a marker of one element at Gauss nodes.
 */
void                rhea_viscosity_marker_set_elem_gauss (
                                                  ymir_vec_t *marker_vec,
                                                  sc_dmatrix_t *marker_el_mat,
                                                  const ymir_locidx_t elid);

#endif /* RHEA_VISCOSITY_H */
