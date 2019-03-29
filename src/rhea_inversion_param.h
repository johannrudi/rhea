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
 * Parameters
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

/* verbosity for paramter output */
typedef enum
{
  RHEA_INVERSION_PARAM_VERBOSE_NONE,
  RHEA_INVERSION_PARAM_VERBOSE_REAL,
  RHEA_INVERSION_PARAM_VERBOSE_REAL_NONDIM,
  RHEA_INVERSION_PARAM_VERBOSE_REAL_NONDIM_DIM
}
rhea_inversion_param_verbosity_t;

/**
 * Prints the active parameters with specified verbosity.
 */
void                rhea_inversion_param_print (
                              ymir_vec_t *parameter_vec,
                              const rhea_inversion_param_verbosity_t verbosity,
                              rhea_inversion_param_t *inv_param);

/******************************************************************************
 * Mapping between Parameter Values and Model Values
 *****************************************************************************/

/**
 * Sets the inversion parameters to the values from a model.
 */
void                rhea_inversion_param_pull_from_model (
                                            ymir_vec_t *parameter_vec,
                                            rhea_inversion_param_t *inv_param);

/**
 * Overwrites model parameters with values from the parameter inversion.
 */
void                rhea_inversion_param_push_to_model (
                                            ymir_vec_t *parameter_vec,
                                            rhea_inversion_param_t *inv_param);

/**
 * Sets the inversion parameters neutral feasible values.
 */
void                rhea_inversion_param_set_neutral (
                                            ymir_vec_t *parameter_vec,
                                            rhea_inversion_param_t *inv_param);

/******************************************************************************
 * Parameter Related Computations
 *****************************************************************************/

/**
 * Computes the gradient vector of the Stokes model w.r.t. model parameters.
 */
void                rhea_inversion_param_compute_gradient (
                                            ymir_vec_t *gradient,
                                            ymir_vec_t *forward_vel_press,
                                            ymir_vec_t *adjoint_vel_press,
                                            rhea_inversion_param_t *inv_param);

/**
 * Computes the norm of the provided gradient vector.
 */
double              rhea_inversion_param_compute_gradient_norm (
                                            ymir_vec_t *gradient,
                                            rhea_inversion_param_t *inv_param);

/**
 * Computes the right-hand side for incremental forward equations.
 */
void                rhea_inversion_param_incremental_forward_rhs (
                                            ymir_vec_t *rhs_vel_mass,
                                            ymir_vec_t *gradient_direction,
                                            ymir_vec_t *forward_vel_press,
                                            rhea_inversion_param_t *inv_param);

/**
 * Applies the Hessian of the inverse problem to a parameter vector.
 */
void                rhea_inversion_param_apply_hessian (
                                            ymir_vec_t *param_vec_out,
                                            ymir_vec_t *param_vec_in,
                                            ymir_vec_t *forward_vel_press,
                                            ymir_vec_t *adjoint_vel_press,
                                            ymir_vec_t *incr_forward_vel_press,
                                            ymir_vec_t *incr_adjoint_vel_press,
                                            rhea_inversion_param_t *inv_param);

/******************************************************************************
 * Parameter Vector
 *****************************************************************************/

/**
 * Creates/destroys a parameter vector.
 */
ymir_vec_t         *rhea_inversion_param_vec_new (
                                            rhea_inversion_param_t *inv_param);

void                rhea_inversion_param_vec_destroy (ymir_vec_t *vec);

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_inversion_param_vec_check_type (
                                            ymir_vec_t *vec,
                                            rhea_inversion_param_t *inv_param);

/**
 * Checks entries of a vector.
 */
int                 rhea_inversion_param_vec_is_valid (
                                            ymir_vec_t *vec,
                                            rhea_inversion_param_t *inv_param);

/**
 * Prints the active values of a parameter vector.
 */
void                rhea_inversion_param_vec_print (
                                            ymir_vec_t *vec,
                                            rhea_inversion_param_t *inv_param);

/**
 * Creates/destroys a reduced parameter vector, where only the active entries
 * are copied from the full parameter vector in the input argument.
 */
sc_dmatrix_t       *rhea_inversion_param_vec_reduced_new (
                                            ymir_vec_t *vec,
                                            rhea_inversion_param_t *inv_param);

void                rhea_inversion_param_vec_reduced_destroy (
                                            sc_dmatrix_t *vec_reduced);

/**
 * Copy entries from a reduced to a full parameter vector.
 */
void                rhea_inversion_param_vec_reduced_copy (
                                            ymir_vec_t *vec,
                                            sc_dmatrix_t *vec_reduced,
                                            rhea_inversion_param_t *inv_param);

/******************************************************************************
 * Data Access
 *****************************************************************************/

/**
 * Returns the number of parameters.
 */
int                 rhea_inversion_param_get_n_parameters (
                                            rhea_inversion_param_t *inv_param);

/**
 * Returns the number of active parameters.
 */
int                 rhea_inversion_param_get_n_active (
                                            rhea_inversion_param_t *inv_param);

/**
 * Gets the pointer to the activation mask.
 */
int                *rhea_inversion_param_get_active (
                                            rhea_inversion_param_t *inv_param);

#endif /* RHEA_INVERSION_PARAM_H */
