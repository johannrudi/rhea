/** RHEA_INVERSION_PARAM
 *
 * Parameters of inverse problems for mantle convection.
 */

#ifndef RHEA_INVERSION_PARAM_H
#define RHEA_INVERSION_PARAM_H

#include <rhea_stokes_problem.h>
#include <rhea_weakzone_label.h>

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
  int                 upper_mantle_arrhenius_activation_energy_a;
  int                 lower_mantle_scaling_a;
  int                 lower_mantle_arrhenius_activation_energy_a;
  int                 stress_exponent_a;
  int                 yield_strength_a;

  /* flags activating inversion for weak zone parameters */
  int                 thickness_a;
  int                 thickness_class_slab_a;
  int                 thickness_class_ridge_a;
  int                 thickness_class_fracture_a;
  int                 thickness_const_a;
  int                 thickness_const_class_slab_a;
  int                 thickness_const_class_ridge_a;
  int                 thickness_const_class_fracture_a;
  int                 weak_factor_interior_a;
  int                 weak_factor_interior_class_slab_a;
  int                 weak_factor_interior_class_ridge_a;
  int                 weak_factor_interior_class_fracture_a;
  int                 weak_factor_interior_label_slab_a;
  int                 weak_factor_interior_label_ridge_a;
  int                 weak_factor_interior_label_fracture_a;
  char               *weak_factor_interior_label_file_path_txt_a;

  /* mean of Gaussian prior */
  double              prior_mean_perturb_stddev;

  /* standard deviation of Gaussian prior for viscosity parameters */
  double              prior_stddev_min;
  double              prior_stddev_max;
  double              prior_stddev_upper_mantle_scaling;
  double              prior_stddev_upper_mantle_arrhenius_activation_energy;
  double              prior_stddev_lower_mantle_scaling;
  double              prior_stddev_lower_mantle_arrhenius_activation_energy;
  double              prior_stddev_stress_exponent;
  double              prior_stddev_yield_strength;

  /* standard deviation of Gaussian prior for weak zone parameters */
  double              prior_stddev_thickness;
  double              prior_stddev_thickness_const;
  double              prior_stddev_weak_factor_interior;

  /* initial guess */
  char               *initial_guess_file_path_txt;
  double              initial_guess_perturb_stddev;
  double              initial_guess_shift_by_prior_stddev;
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
 * Transfer between Parameter Values and Model Values
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
 * Sets the inversion parameters to initial values.
 */
void                rhea_inversion_param_set_initial_guess (
                                          ymir_vec_t *parameter_vec,
                                          rhea_inversion_param_t *inv_param,
                                          rhea_inversion_param_options_t *opt);

/**
 * Computes the derivative of the parametrization function.
 *
 *  (d/dp *_convert_to_model_*) (*_convert_from_model_* (m))
 *
 * where
 *   *_min_max () Used for lower and upper bounds of viscosity.
 *   *_scal ()    Used for scaling factors (pre-factors) of viscosity.
 *   *_arrh ()    Used for Arrhenius activation energy.
 *   *_n ()       Used for stress exponent.
 *   *_yield ()   Used for yield strength.
 *   *_weak ()    Used for weak zone factors.
 */
double              rhea_inversion_param_derivative_min_max (
                                                  const double model_val);

double              rhea_inversion_param_derivative_scal (
                                                  const double model_val);

double              rhea_inversion_param_derivative_arrh (
                                                  const double model_val);

double              rhea_inversion_param_derivative_n (
                                                  const double model_val);

double              rhea_inversion_param_derivative_yield (
                                                  const double model_val);

double              rhea_inversion_param_derivative_weak (
                                                  const double model_val);

/******************************************************************************
 * Parameter Related Computations
 *****************************************************************************/

/**
 * Callback function that computes the derivative with respect to a single
 * parameter.
 *
 * \param [in] parameter_idx  Index of parameter.
 */
typedef double    (*rhea_inversion_param_derivative_fn_t) (
                                  ymir_vec_t *forward_vel_press,
                                  ymir_vec_t *adjoint_vel_press,
                                  const int parameter_idx,
                                  const int derivative_type,
                                  const rhea_weakzone_label_t weak_label,
                                  ymir_vec_t *coeff_param_derivative,
                                  ymir_stress_op_t *stress_op_param_derivative,
                                  void *data);

/**
 * Computes the (squared norm of the) prior term for the objective functional.
 */
double              rhea_inversion_param_prior (
                                        ymir_vec_t *parameter_vec,
                                        rhea_inversion_param_t *inv_param);
/**
 * Computes the gradient vector of the Stokes model w.r.t. model parameters.
 */
void                rhea_inversion_param_compute_gradient (
                            ymir_vec_t *gradient_vec,
                            ymir_vec_t *parameter_vec,
                            ymir_vec_t *forward_vel_press,
                            ymir_vec_t *adjoint_vel_press,
                            const double prior_weight,
                            rhea_inversion_param_t *inv_param,
                            ymir_vec_t *gradient_adjoint_comp_vec,
                            ymir_vec_t *gradient_prior_comp_vec,
                            ymir_vec_t *gradient_callback_comp_vec,
                            rhea_inversion_param_derivative_fn_t derivative_fn,
                            void *derivative_fn_data);

/**
 * Computes the norm of the provided gradient vector.
 */
double              rhea_inversion_param_compute_gradient_norm (
                                        ymir_vec_t *gradient_vec,
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
                                        const double prior_weight,
                                        rhea_inversion_param_t *inv_param);

/**
 * Computes the (diagonal) inverse covariance of the prior.
 */
void                rhea_inversion_param_prior_inv_cov (
                                        ymir_vec_t *inverse_covariance,
                                        rhea_inversion_param_t *inv_param);

/**
 * Restricts the parameters to their feasible values (elementwise).
 */
int                 rhea_inversion_param_restrict_to_feasible (
                                        ymir_vec_t *parameter_vec,
                                        const double restrict_to_prior_stddev,
                                        rhea_inversion_param_t *inv_param);

/**
 * Modifies the step such that the new parameter vector,
 *   parameter_vec + step_vec,
 * lies within the feasible set associated with the inversion parameters.
 */
int                 rhea_inversion_param_restrict_step_length_to_feasible (
                                        ymir_vec_t *step_vec,
                                        ymir_vec_t *parameter_vec,
                                        const double restrict_to_prior_stddev,
                                        rhea_inversion_param_t *inv_param,
                                        double *step_length_new);

/******************************************************************************
 * Parameter Vector
 *****************************************************************************/

/**
 * Creates/destroys a parameter vector.
 */
ymir_vec_t         *rhea_inversion_param_vec_new (
                                            rhea_inversion_param_t *inv_param);

ymir_vec_t         *rhea_inversion_param_vec_new_perturb (
                                            rhea_inversion_param_t *inv_param,
                                            const double perturb_stddev);

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

/**
 * Applies a reduced matrix to a full parameter vector.
 */
void                rhea_inversion_param_vec_reduced_apply_matrix (
                                            ymir_vec_t *out,
                                            sc_dmatrix_t *matrix_reduced,
                                            ymir_vec_t *in,
                                            rhea_inversion_param_t *inv_param);

/**
 * Solves a linear system with a reduced matrix and a right-hand side given by
 * a full parameter vector.
 */
int                 rhea_inversion_param_vec_reduced_solve_matrix (
                                            ymir_vec_t *sol,
                                            sc_dmatrix_t *matrix_reduced,
                                            ymir_vec_t *rhs,
                                            rhea_inversion_param_t *inv_param,
                                            int *n_iterations);

/**
 * Computes inner product of two vectors.
 */
double              rhea_inversion_param_vec_reduced_inner_product (
                                            ymir_vec_t *vecL,
                                            ymir_vec_t *vecR,
                                            ymir_vec_t *weight,
                                            rhea_inversion_param_t *inv_param,
                                            const int normalize_wrt_size);

double              rhea_inversion_param_vec_reduced_ip (
                                            sc_dmatrix_t *vecL_reduced,
                                            sc_dmatrix_t *vecR_reduced,
                                            sc_dmatrix_t *weight_reduced,
                                            const int normalize_wrt_size);

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
