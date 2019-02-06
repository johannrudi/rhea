#include <rhea_viscosity_param_derivative.h>
#include <rhea_base.h>
#include <rhea_temperature.h>
#include <ymir_interp_vec.h>

/**
 * Computes derivative w.r.t. scaling factor of upper/lower mantle:
 *   deriv = visc
 */
static void
rhea_viscosity_param_derivative_scal (
                        ymir_vec_t *derivative,
                        ymir_vec_t *viscosity,
                        ymir_vec_t *bounds_marker,
                        ymir_vec_t *yielding_marker,
                        rhea_viscosity_options_t *visc_options,
                        const int filter_upper_mantle,
                        const int filter_lower_mantle)
{
  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (derivative));
  RHEA_ASSERT (!filter_upper_mantle || !filter_lower_mantle);

    /* initialize derivative */
  ymir_vec_copy (viscosity, derivative);

  /* filter upper/lower mantle */
  if (filter_upper_mantle || filter_lower_mantle) {
    ymir_vec_t         *filter = ymir_vec_template (derivative);

    /* create filter */
    if (filter_upper_mantle) {
      rhea_viscosity_stats_filter_upper_mantle (filter,
                                                visc_options->domain_options);
    }
    else if (filter_lower_mantle) {
      rhea_viscosity_stats_filter_lower_mantle (filter,
                                                visc_options->domain_options);
    }
    else {
      RHEA_ABORT_NOT_REACHED ();
    }

    /* apply filter */
    ymir_vec_multiply_in (filter, derivative);
    ymir_vec_destroy (filter);
  }

  /* set to zero where max bound reached */
  rhea_viscosity_filter_where_max (derivative, bounds_marker,
                                   1 /* invert_filter */);

  /* set to zero where yielding */
  rhea_viscosity_filter_where_yielding (derivative, yielding_marker,
                                        1 /* invert_filter */);

  /* check output */
  RHEA_ASSERT (rhea_viscosity_is_valid (derivative));
}

/**
 * Computes derivative w.r.t. activation energy in Arrhenius relationship:
 *   deriv = E_a * (0.5 - T) * visc
 */
static void
rhea_viscosity_param_derivative_arrh (
                        ymir_vec_t *derivative,
                        ymir_vec_t *viscosity,
                        ymir_vec_t *bounds_marker,
                        ymir_vec_t *yielding_marker,
                        ymir_vec_t *temperature,
                        rhea_viscosity_options_t *visc_options,
                        const int filter_upper_mantle,
                        const int filter_lower_mantle)
{
  ymir_vec_t         *scal_filter = ymir_vec_template (derivative);
  ymir_vec_t         *scal = ymir_vec_template (derivative);

  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (derivative));
  RHEA_ASSERT (!filter_upper_mantle || !filter_lower_mantle);

  /* initialize derivative */
  ymir_vec_copy (viscosity, derivative);

  /* multiply derivative by temperature difference term */
  ymir_interp_vec (temperature, scal);
  ymir_vec_scale_shift (-1.0, -RHEA_TEMPERATURE_NEUTRAL_VALUE, scal);
  ymir_vec_multiply_in (scal, derivative);

  /* initialize scaling factor for the derivative */
  ymir_vec_set_zero (scal);

  /* add contribution to the scaling factor from the upper mantle */
  if (!filter_lower_mantle) {
    rhea_viscosity_stats_filter_upper_mantle (scal_filter,
                                              visc_options->domain_options);
  }
  else { /* if discard upper mantle values */
    ymir_vec_set_zero (scal_filter);
  }
  ymir_vec_add (visc_options->upper_mantle_arrhenius_activation_energy,
                scal_filter, scal);

  /* add contribution to the scaling factor from the lower mantle */
  if (!filter_upper_mantle) {
    rhea_viscosity_stats_filter_lower_mantle (scal_filter,
                                              visc_options->domain_options);
  }
  else { /* if discard lower mantle values */
    ymir_vec_set_zero (scal_filter);
  }
  ymir_vec_add (visc_options->lower_mantle_arrhenius_activation_energy,
                scal_filter, scal);

  /* apply scaling factor to the derivative */
  ymir_vec_multiply_in (scal, derivative);
  ymir_vec_destroy (scal_filter);
  ymir_vec_destroy (scal);

  /* set to zero where max bound reached */
  rhea_viscosity_filter_where_max (derivative, bounds_marker,
                                   1 /* invert_filter */);

  /* set to zero where yielding */
  rhea_viscosity_filter_where_yielding (derivative, yielding_marker,
                                        1 /* invert_filter */);

  /* check output */
  RHEA_ASSERT (rhea_viscosity_is_valid (derivative));
}

void
rhea_viscosity_param_derivative (
                            ymir_vec_t *derivative,
                            rhea_viscosity_param_derivative_t derivative_type,
                            ymir_vec_t *viscosity,
                            ymir_vec_t *bounds_marker,
                            ymir_vec_t *yielding_marker,
                            ymir_vec_t *temperature,
                            ymir_vec_t *weakzone,
                            rhea_viscosity_options_t *visc_options)
{
  switch (derivative_type) {
//case RHEA_VISCOSITY_PARAM_DERIVATIVE_MIN:
//case RHEA_VISCOSITY_PARAM_DERIVATIVE_MAX:
  case RHEA_VISCOSITY_PARAM_DERIVATIVE_UPPER_MANTLE_SCALING:
    rhea_viscosity_param_derivative_scal (
        derivative, viscosity, bounds_marker, yielding_marker, visc_options,
        1 /* upper mantle */, 0 /* not lower mantle */);
    break;
  case RHEA_VISCOSITY_PARAM_DERIVATIVE_LOWER_MANTLE_SCALING:
    rhea_viscosity_param_derivative_scal (
        derivative, viscosity, bounds_marker, yielding_marker, visc_options,
        0 /* not upper mantle */, 1 /* lower mantle */);
    break;
  case RHEA_VISCOSITY_PARAM_DERIVATIVE_UPPER_MANTLE_ACTIVATION_ENERGY:
    rhea_viscosity_param_derivative_arrh (
        derivative, viscosity, bounds_marker, yielding_marker, temperature,
        visc_options, 1 /* upper mantle */, 0 /* not lower mantle */);
    break;
  case RHEA_VISCOSITY_PARAM_DERIVATIVE_LOWER_MANTLE_ACTIVATION_ENERGY:
    rhea_viscosity_param_derivative_arrh (
        derivative, viscosity, bounds_marker, yielding_marker, temperature,
        visc_options, 0 /* not upper mantle */, 1 /* lower mantle */);
    break;
//case RHEA_VISCOSITY_PARAM_DERIVATIVE_STRESS_EXPONENT:
//case RHEA_VISCOSITY_PARAM_DERIVATIVE_YIELD_STRENGTH:
//case RHEA_VISCOSITY_PARAM_DERIVATIVE_WEAK_THICKNESS:
//case RHEA_VISCOSITY_PARAM_DERIVATIVE_WEAK_THICKNESS_CONST:
//case RHEA_VISCOSITY_PARAM_DERIVATIVE_WEAK_FACTOR_INTERIOR:
  default: /* unknown derivative type */
    RHEA_ABORT_NOT_REACHED ();
  }
}
