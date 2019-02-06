/** RHEA_VISCOSITY_PARAM_DERIVATIVE
 *
 * Viscosity modeling mantle convection.
 */

#ifndef RHEA_VISCOSITY_PARAM_DERIVATIVE_H
#define RHEA_VISCOSITY_PARAM_DERIVATIVE_H

#include <rhea_viscosity.h>

/* list of parameters w.r.t. which viscosity derivatives exist */
typedef enum
{
  RHEA_VISCOSITY_PARAM_DERIVATIVE_MIN,
  RHEA_VISCOSITY_PARAM_DERIVATIVE_MAX,
  RHEA_VISCOSITY_PARAM_DERIVATIVE_UPPER_MANTLE_SCALING,
  RHEA_VISCOSITY_PARAM_DERIVATIVE_LOWER_MANTLE_SCALING,
  RHEA_VISCOSITY_PARAM_DERIVATIVE_UPPER_MANTLE_ACTIVATION_ENERGY,
  RHEA_VISCOSITY_PARAM_DERIVATIVE_LOWER_MANTLE_ACTIVATION_ENERGY,
  RHEA_VISCOSITY_PARAM_DERIVATIVE_STRESS_EXPONENT,
  RHEA_VISCOSITY_PARAM_DERIVATIVE_YIELD_STRENGTH,
  RHEA_VISCOSITY_PARAM_DERIVATIVE_WEAK_THICKNESS,
  RHEA_VISCOSITY_PARAM_DERIVATIVE_WEAK_THICKNESS_CONST,
  RHEA_VISCOSITY_PARAM_DERIVATIVE_WEAK_FACTOR_INTERIOR
}
rhea_viscosity_param_derivative_t;

/**
 * Computes a derivative of the viscosity w.r.t. the parameter that is specified
 * by `derivative_type`.
 */
void                rhea_viscosity_param_derivative (
                            ymir_vec_t *derivative,
                            rhea_viscosity_param_derivative_t derivative_type,
                            ymir_vec_t *viscosity,
                            ymir_vec_t *bounds_marker,
                            ymir_vec_t *yielding_marker,
                            ymir_vec_t *temperature,
                            ymir_vec_t *weakzone,
                            rhea_viscosity_options_t *visc_options);

#endif /* RHEA_VISCOSITY_PARAM_DERIVATIVE_H */
