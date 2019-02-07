/** RHEA_INVERSION_PARAM
 *
 * Parameters of inverse problems for mantle convection.
 */

#ifndef RHEA_INVERSION_PARAM_H
#define RHEA_INVERSION_PARAM_H

#include <rhea_stokes_problem.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* parameter options regarding the inverse problem */
typedef struct rhea_inversion_param_options
{
  /* flags activating inversion for viscosity parameters */
  int                 min_a;
  int                 max_a;
  int                 upper_mantle_scaling_a;
  int                 lower_mantle_scaling_a;
  int                 upper_mantle_arrhenius_activation_energy_a;
  int                 lower_mantle_arrhenius_activation_energy_a;
  int                 stress_exponent_a;
  int                 yield_strength_a;

  /* flags activating inversion for weak zone parameters */
  int                 thickness_a;
  int                 thickness_generic_slab_a;
  int                 thickness_generic_ridge_a;
  int                 thickness_generic_fracture_a;
  int                 thickness_const_a;
  int                 thickness_const_generic_slab_a;
  int                 thickness_const_generic_ridge_a;
  int                 thickness_const_generic_fracture_a;
  int                 weak_factor_interior_a;
  int                 weak_factor_interior_generic_slab_a;
  int                 weak_factor_interior_generic_ridge_a;
  int                 weak_factor_interior_generic_fracture_a;
  int                 weak_factor_interior_earth_slab_a;
  int                 weak_factor_interior_earth_ridge_a;
  int                 weak_factor_interior_earth_fracture_a;
}
rhea_inversion_param_options_t;

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_inversion_param_add_options (
                              rhea_inversion_param_options_t *inv_param_options,
                              ymir_options_t *opt_sup);

/******************************************************************************
 * Inversion Parameters
 *****************************************************************************/

/* inversion parameters (opaque) */
typedef struct rhea_inversion_param rhea_inversion_param_t;

/**
 * Creates a new set of inversion parameters.
 */
rhea_inversion_param_t *rhea_inversion_param_new (
                            rhea_stokes_problem_t *stokes_problem,
                            rhea_inversion_param_options_t *inv_param_options);

/**
 * Destroys a set of inversion parameters.
 */
void                rhea_inversion_param_destroy (
                            rhea_inversion_param_t *inv_param);

/**
 * Sets the inversion parameters to the values from a model.
 */
void                rhea_inversion_param_pull_from_model (
                            rhea_inversion_param_t *inv_param);

/**
 * Overwrites model parameters with values from the parameter inversion.
 */
void                rhea_inversion_param_push_to_model (
                            rhea_inversion_param_t *inv_param);

/**
 * Computes the gradient vector of the Stokes model w.r.t. model parameters.
 */
void                rhea_inversion_param_compute_gradient (
                            ymir_vec_t *gradient,
                            ymir_vec_t *forward_vel_press,
                            ymir_vec_t *adjoint_vel_press,
                            rhea_inversion_param_t *inv_param);

/******************************************************************************
 * Data Access
 *****************************************************************************/

/**
 * Gets the pointer to the parameter vector.
 */
ymir_vec_t         *rhea_inversion_param_get_vector (
                            rhea_inversion_param_t *inv_param);

#endif /* RHEA_INVERSION_PARAM_H */
