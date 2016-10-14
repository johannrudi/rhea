/*
 */

#ifndef RHEA_VISCOSITY_H
#define RHEA_VISCOSITY_H

#include <rhea_domain.h>

/* enumerator for types of viscosities */
typedef enum
{
  SL_VISCOSITY_CONST,
  SL_VISCOSITY_LINEAR,
  SL_VISCOSITY_NONLINEAR
}
slabs_viscosity_t;

/* enumerator for types of initial viscosities for nonlinear Stokes problem*/
typedef enum
{
  SL_VISCOSITY_INIT_NL_STOKES_DEFAULT,
  SL_VISCOSITY_INIT_NL_STOKES_CONST_CENTERED_TO_BOUNDS,
  SL_VISCOSITY_INIT_NL_STOKES_CONST_MIN_BOUND,
  SL_VISCOSITY_INIT_NL_STOKES_TEMP,
  SL_VISCOSITY_INIT_NL_STOKES_TEMP_UM_REL_TO_LM
}
slabs_viscosity_init_nl_stokes_t;

/* enumerator for types of viscosity models */
typedef enum
{
  /* (1) upper bound, (2) weak zone, (3) yielding, (4) lower bound
   * viscosity bounds via cut-off */
  SL_VISCOSITY_MODEL_UWYL,

  /* (1) upper bound, (2) weak zone, (3) yielding, (4) lower bound
   * upper viscosity bound via cut-off, lower viscosity bound via addition */
  SL_VISCOSITY_MODEL_UWYL_LREG,

  /* (1) upper bound + shift, (2) weak zone, (3) yielding, (4) lower bound
   * upper viscosity bound via shift, lower viscosity bound via addition */
  SL_VISCOSITY_MODEL_UWYL_SHIFT_LREG
}
slabs_viscosity_model_t;

/* parameter list for mantle flow physics */
typedef struct slabs_physics_options
{
  slabs_viscosity_t                 type;
  slabs_viscosity_init_nl_stokes_t  type_for_init_nonlinear;
  slabs_viscosity_model_t           model_type;

  double              min;
  double              max;

  double              upper_mantle_scaling;
  double              upper_mantle_activation_energy;
  double              lower_mantle_scaling;
  double              lower_mantle_activation_energy;

  double              stress_exponent;
  double              stress_yield;

  rhea_domain_options_t  *domain_options;
}
slabs_physics_options_t;

#endif /* RHEA_VISCOSITY_H */
