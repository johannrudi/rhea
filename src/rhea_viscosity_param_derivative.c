#include <rhea_viscosity_param_derivative.h>
#include <rhea_base.h>
#include <rhea_temperature.h>
#include <rhea_strainrate.h>
#include <ymir_interp_vec.h>

static void
rhea_viscosity_param_derivative_init (ymir_vec_t *derivative,
                                      ymir_vec_t *viscosity,
                                      ymir_vec_t *bounds_marker,
                                      ymir_vec_t *yielding_marker,
                                      rhea_viscosity_options_t *visc_options,
                                      const int remove_bound_min,
                                      const int remove_bound_max,
                                      const int remove_yielding,
                                      const int remove_lower_mantle,
                                      const int remove_upper_mantle)
{
  /* check input */
  RHEA_ASSERT (viscosity != NULL);
  RHEA_ASSERT (visc_options != NULL);

  /* initialize derivative */
  ymir_vec_copy (viscosity, derivative);

  /* remove lower viscosity bound */
  if (remove_bound_min) {
    const double        visc_min = visc_options->min;

    switch (visc_options->model) {
    case RHEA_VISCOSITY_MODEL_UWYL:
      RHEA_ASSERT (bounds_marker != NULL);
      rhea_viscosity_filter_where_min (derivative, bounds_marker,
                                       1 /* invert_filter */);
      break;
    case RHEA_VISCOSITY_MODEL_UWYL_LADD_UCUT:
    case RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT:
      if (0.0 < visc_min) {
        ymir_vec_shift (-visc_min, derivative);
      }
      break;
    default: /* unknown viscosity model */
      RHEA_ABORT_NOT_REACHED ();
    }
  }

  /* remove upper viscosity bound */
  if (remove_bound_max) {
    RHEA_ASSERT (bounds_marker != NULL);
    rhea_viscosity_filter_where_max (derivative, bounds_marker,
                                     1 /* invert_filter */);
  }

  /* remove yielding */
  if (remove_yielding) {
    RHEA_ASSERT (yielding_marker != NULL);
    rhea_viscosity_filter_where_yielding (derivative, yielding_marker,
                                          1 /* invert_filter */);
  }

  /* remove lower mantle */
  if (remove_lower_mantle) {
    ymir_vec_t         *filter = ymir_vec_template (derivative);

    rhea_viscosity_stats_filter_upper_mantle (filter,
                                              visc_options->domain_options);
    ymir_vec_multiply_in (filter, derivative);
    ymir_vec_destroy (filter);
  }

  /* remove upper mantle */
  if (remove_upper_mantle) {
    ymir_vec_t         *filter = ymir_vec_template (derivative);

    rhea_viscosity_stats_filter_lower_mantle (filter,
                                              visc_options->domain_options);
    ymir_vec_multiply_in (filter, derivative);
    ymir_vec_destroy (filter);
  }
}

/**
 * Computes derivative w.r.t. the lower viscosity bound.
 */
static void
rhea_viscosity_param_derivative_min (ymir_vec_t *derivative,
                                     ymir_vec_t *viscosity,
                                     ymir_vec_t *bounds_marker,
                                     ymir_vec_t *yielding_marker,
                                     rhea_viscosity_options_t *visc_options)
{
  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (derivative));

  /* set derivative */
  switch (visc_options->model) {
  case RHEA_VISCOSITY_MODEL_UWYL:
    rhea_viscosity_param_derivative_init (
        derivative, viscosity, bounds_marker, yielding_marker, visc_options,
        0 /* !remove_bound_min */,
        0 /* !remove_bound_max */,
        1 /* remove_yielding */,
        0 /* !remove_lower_mantle */,
        0 /* !remove_upper_mantle */);
    rhea_viscosity_filter_where_min (derivative, bounds_marker,
                                     0 /* !invert_filter */);
    break;
  case RHEA_VISCOSITY_MODEL_UWYL_LADD_UCUT:
  case RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT:
    ymir_vec_set_value (derivative, visc_options->min);
    break;
  default: /* unknown viscosity model */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* check output */
  RHEA_ASSERT (sc_dmatrix_is_valid (derivative->dataown));
}

/**
 * Computes derivative w.r.t. the upper viscosity bound.
 */
static void
rhea_viscosity_param_derivative_max (ymir_vec_t *derivative,
                                     ymir_vec_t *viscosity,
                                     ymir_vec_t *bounds_marker,
                                     ymir_vec_t *yielding_marker,
                                     rhea_viscosity_options_t *visc_options)
{
  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (derivative));

  /* set derivative */
  rhea_viscosity_param_derivative_init (
      derivative, viscosity, bounds_marker, yielding_marker, visc_options,
      1 /* remove_bound_min */,
      0 /* !remove_bound_max */,
      1 /* remove_yielding */,
      0 /* !remove_lower_mantle */,
      0 /* !remove_upper_mantle */);
  rhea_viscosity_filter_where_max (derivative, bounds_marker,
                                   0 /* !invert_filter */);

  /* check output */
  RHEA_ASSERT (sc_dmatrix_is_valid (derivative->dataown));
}

/**
 * Computes derivative w.r.t. the scaling factor of upper/lower mantle:
 *   c = exp(param_c)
 *   deriv = visc
 */
static void
rhea_viscosity_param_derivative_scal (ymir_vec_t *derivative,
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

  /* set derivative */
  rhea_viscosity_param_derivative_init (
      derivative, viscosity, bounds_marker, yielding_marker, visc_options,
      1 /* remove_bound_min */,
      1 /* remove_bound_max */,
      1 /* remove_yielding */,
      filter_upper_mantle,
      filter_lower_mantle);

  /* check output */
  RHEA_ASSERT (sc_dmatrix_is_valid (derivative->dataown));
}

/**
 * Computes derivative w.r.t. the activation energy in Arrhenius relationship:
 *   E_a = exp(param_E_a)
 *   deriv = E_a * (0.5 - T) * visc
 */
static void
rhea_viscosity_param_derivative_arrh (ymir_vec_t *derivative,
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
  rhea_viscosity_param_derivative_init (
      derivative, viscosity, bounds_marker, yielding_marker, visc_options,
      1 /* remove_bound_min */,
      1 /* remove_bound_max */,
      1 /* remove_yielding */,
      filter_upper_mantle,
      filter_lower_mantle);

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

  /* check output */
  RHEA_ASSERT (sc_dmatrix_is_valid (derivative->dataown));
}

static void
rhea_viscosity_param_derivative_log_fn (double *val, double x, double y,
                                        double z, ymir_locidx_t nodeid,
                                        void *data)
{
  *val = log (*val);
}

/**
 * Computes derivative w.r.t. the stress exponent:
 *   n = 1 + exp(param_n)
 *   deriv = - log(strainrate - shift) * exp(param_n) / n^2 * visc
 */
static void
rhea_viscosity_param_derivative_stress_exp (
                                        ymir_vec_t *derivative,
                                        ymir_vec_t *viscosity,
                                        ymir_vec_t *bounds_marker,
                                        ymir_vec_t *yielding_marker,
                                        ymir_vec_t *velocity,
                                        rhea_viscosity_options_t *visc_options)
{
  const double        stress_exp = visc_options->stress_exponent;
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (derivative);
  ymir_vec_t         *strt_scal = rhea_strainrate_2inv_new (ymir_mesh);
  double              scal;

  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (derivative));

  /* initialize derivative */
  rhea_viscosity_param_derivative_init (
      derivative, viscosity, bounds_marker, yielding_marker, visc_options,
      1 /* remove_bound_min */,
      1 /* remove_bound_max */,
      1 /* remove_yielding */,
      1 /* remove_lower_mantle */,
      0 /* !remove_upper_mantle */);

  /* compute sqrt of the 2nd invariant of the strain rate */
  switch (visc_options->model) {
  case RHEA_VISCOSITY_MODEL_UWYL:
  case RHEA_VISCOSITY_MODEL_UWYL_LADD_UCUT:
    rhea_strainrate_compute_sqrt_of_2inv (strt_scal, velocity);
    break;
  case RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT:
    RHEA_ABORT_NOT_REACHED (); //TODO implement this; need a shift apply fnc.
    break;
  default: /* unknown viscosity model */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* multiply-in the strain rate dependent scaling */
  ymir_dvec_set_function (strt_scal, rhea_viscosity_param_derivative_log_fn,
                          NULL);
  ymir_vec_multiply_in (strt_scal, derivative);
  rhea_strainrate_2inv_destroy (strt_scal);

  /* multiply by constant scaling */
  scal = - SC_MAX (0.0, stress_exp - 1.0) / (stress_exp*stress_exp);
  ymir_vec_scale (scal, derivative);

  /* check output */
  RHEA_ASSERT (sc_dmatrix_is_valid (derivative->dataown));
}

/**
 * Computes derivative w.r.t. the yield strength:
 *   n = exp(param_yield)
 *   deriv = visc
 */
static void
rhea_viscosity_param_derivative_yield_strength (
                                        ymir_vec_t *derivative,
                                        ymir_vec_t *viscosity,
                                        ymir_vec_t *bounds_marker,
                                        ymir_vec_t *yielding_marker,
                                        rhea_viscosity_options_t *visc_options)
{
  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (derivative));

  /* set derivative */
  rhea_viscosity_param_derivative_init (
      derivative, viscosity, bounds_marker, yielding_marker, visc_options,
      1 /* remove_bound_min */,
      0 /* !remove_bound_max */,
      0 /* !remove_yielding */,
      1 /* remove_lower_mantle */,
      0 /* !remove_upper_mantle */);
  rhea_viscosity_filter_where_yielding (derivative, yielding_marker,
                                        0 /* !invert_filter */);

  /* check output */
  RHEA_ASSERT (sc_dmatrix_is_valid (derivative->dataown));
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
                            ymir_vec_t *velocity,
                            rhea_viscosity_options_t *visc_options)
{
  switch (derivative_type) {

  case RHEA_VISCOSITY_PARAM_DERIVATIVE_MIN:
    if (rhea_viscosity_restrict_min (visc_options)) {
      rhea_viscosity_param_derivative_min (
          derivative, viscosity, bounds_marker, yielding_marker, visc_options);
    }
    else {
      ymir_vec_set_zero (derivative);
    }
    break;

  case RHEA_VISCOSITY_PARAM_DERIVATIVE_MAX:
    if (rhea_viscosity_restrict_max (visc_options)) {
      rhea_viscosity_param_derivative_max (
          derivative, viscosity, bounds_marker, yielding_marker, visc_options);
    }
    else {
      ymir_vec_set_zero (derivative);
    }
    break;

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
    if (rhea_viscosity_has_arrhenius (visc_options)) {
      rhea_viscosity_param_derivative_arrh (
          derivative, viscosity, bounds_marker, yielding_marker, temperature,
          visc_options, 1 /* upper mantle */, 0 /* not lower mantle */);
    }
    else {
      ymir_vec_set_zero (derivative);
    }
    break;

  case RHEA_VISCOSITY_PARAM_DERIVATIVE_LOWER_MANTLE_ACTIVATION_ENERGY:
    if (rhea_viscosity_has_arrhenius (visc_options)) {
      rhea_viscosity_param_derivative_arrh (
          derivative, viscosity, bounds_marker, yielding_marker, temperature,
          visc_options, 0 /* not upper mantle */, 1 /* lower mantle */);
    }
    else {
      ymir_vec_set_zero (derivative);
    }
    break;

  case RHEA_VISCOSITY_PARAM_DERIVATIVE_STRESS_EXPONENT:
    if (rhea_viscosity_has_strain_rate_weakening (visc_options)) {
      rhea_viscosity_param_derivative_stress_exp (
          derivative, viscosity, bounds_marker, yielding_marker, velocity,
          visc_options);
    }
    else {
      ymir_vec_set_zero (derivative);
    }
    break;

  case RHEA_VISCOSITY_PARAM_DERIVATIVE_YIELD_STRENGTH:
    if (rhea_viscosity_has_yielding (visc_options)) {
      rhea_viscosity_param_derivative_yield_strength (
          derivative, viscosity, bounds_marker, yielding_marker, visc_options);
    }
    else {
      ymir_vec_set_zero (derivative);
    }
    break;

//case RHEA_VISCOSITY_PARAM_DERIVATIVE_WEAK_THICKNESS:
//case RHEA_VISCOSITY_PARAM_DERIVATIVE_WEAK_THICKNESS_CONST:
//case RHEA_VISCOSITY_PARAM_DERIVATIVE_WEAK_FACTOR_INTERIOR:

  default: /* unknown derivative type */
    RHEA_ABORT_NOT_REACHED ();
  }
}
