/*
 */

#ifndef RHEA_VISCOSITY_H
#define RHEA_VISCOSITY_H

#include <rhea_domain.h>

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
  RHEA_VISCOSITY_MODEL_UWYL_LREG,

  /* (1) upper bound + shift, (2) weak zone, (3) yielding, (4) lower bound
   * upper viscosity bound via shift, lower viscosity bound via addition */
  RHEA_VISCOSITY_MODEL_UWYL_SHIFT_LREG
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

  /* strain rate weakening is governed by this exponent (aka. `n`) */
  double              stress_exponent;

  /* value of viscous stress above which plastic yielding occurs */
  double              yield_stress;

  /* options & properties of the computational domain */
  rhea_domain_options_t  *domain_options;
}
rhea_viscosity_options_t;

#endif /* RHEA_VISCOSITY_H */
