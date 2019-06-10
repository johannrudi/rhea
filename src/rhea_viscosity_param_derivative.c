#include <rhea_viscosity_param_derivative.h>
#include <rhea_base.h>
#include <rhea_inversion_param.h>
#include <rhea_temperature.h>
#include <rhea_strainrate.h>
#include <ymir_interp_vec.h>

#define RHEA_VISCOSITY_PARAM_DERIVATIVE_VTK (0)

#if (0 < RHEA_VISCOSITY_PARAM_DERIVATIVE_VTK)
# include <ymir_vtk.h>
#endif

static void
rhea_viscosity_param_derivative_init (ymir_vec_t *derivative,
                                      ymir_vec_t *viscosity,
                                      ymir_vec_t *bounds_marker,
                                      ymir_vec_t *yielding_marker,
                                      rhea_viscosity_options_t *visc_options,
                                      const int remove_bound_min,
                                      const int remove_bound_max,
                                      const int remove_yielding,
                                      const int remove_upper_mantle,
                                      const int remove_lower_mantle)
{
  rhea_domain_options_t  *domain_options = visc_options->domain_options;

  /* check input */
  RHEA_ASSERT (viscosity != NULL);
  RHEA_ASSERT (visc_options != NULL);

  /* initialize derivative */
  ymir_vec_copy (viscosity, derivative);

  /* remove lower viscosity bound */
  if (remove_bound_min) {
    switch (visc_options->model) {
    case RHEA_VISCOSITY_MODEL_UWYL:
      RHEA_ASSERT (bounds_marker != NULL);
      rhea_viscosity_filter_where_min (derivative, bounds_marker,
                                       1 /* invert_filter */);
      break;
    case RHEA_VISCOSITY_MODEL_UWYL_LADD_UCUT:
    case RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT:
      if (0.0 < visc_options->min) {
        ymir_vec_shift (-visc_options->min, derivative);
        ymir_vec_bound_min (derivative, 0.0);
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

    rhea_viscosity_stats_filter_upper_mantle (filter, domain_options);
    ymir_vec_multiply_in (filter, derivative);
    ymir_vec_destroy (filter);
  }

  /* remove upper mantle */
  if (remove_upper_mantle) {
    ymir_vec_t         *filter = ymir_vec_template (derivative);

    rhea_viscosity_stats_filter_lower_mantle (filter, domain_options);
    ymir_vec_multiply_in (filter, derivative);
    ymir_vec_destroy (filter);
  }
}

static void
rhea_viscosity_param_derivative_scale_mantle (
                                      ymir_vec_t *derivative,
                                      const double scaling_upper_mantle,
                                      const double scaling_lower_mantle,
                                      rhea_viscosity_options_t *visc_options)
{
  rhea_domain_options_t  *domain_options = visc_options->domain_options;
  const int           process_upper = isfinite (scaling_upper_mantle);
  const int           process_lower = isfinite (scaling_lower_mantle);


  if (process_upper && !process_lower) { /* if upper mantle */
    ymir_vec_scale (scaling_upper_mantle, derivative);
  }
  else if (!process_upper && process_lower) { /* if lower mantle */
    ymir_vec_scale (scaling_lower_mantle, derivative);
  }
  else if (process_upper && process_lower) { /* if whole mantle */
    ymir_vec_t         *scal_filter = ymir_vec_template (derivative);
    ymir_vec_t         *scal = ymir_vec_template (derivative);

    ymir_vec_set_zero (scal);

    rhea_viscosity_stats_filter_upper_mantle (scal_filter, domain_options);
    ymir_vec_add (scaling_upper_mantle, scal_filter, scal);
    rhea_viscosity_stats_filter_lower_mantle (scal_filter, domain_options);
    ymir_vec_add (scaling_lower_mantle, scal_filter, scal);

    ymir_vec_multiply_in (scal, derivative);

    ymir_vec_destroy (scal_filter);
    ymir_vec_destroy (scal);
  }
  else {
    RHEA_ABORT_NOT_REACHED ();
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
  const double        scaling_prev = visc_options->min;
  const double        scaling_deriv =
    rhea_inversion_param_convert_from_model_pos_deriv (scaling_prev);

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
        0 /* !remove_upper_mantle */,
        0 /* !remove_lower_mantle */);
    rhea_viscosity_filter_where_min (derivative, bounds_marker,
                                     0 /* !invert_filter */);
    ymir_vec_scale (scaling_deriv/scaling_prev, derivative);
    break;
  case RHEA_VISCOSITY_MODEL_UWYL_LADD_UCUT:
  case RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT:
    ymir_vec_set_value (derivative, scaling_deriv);
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
  const double        scaling_prev = visc_options->max;
  const double        scaling_deriv =
    rhea_inversion_param_convert_from_model_pos_deriv (scaling_prev);

  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (derivative));

  /* set derivative */
  rhea_viscosity_param_derivative_init (
      derivative, viscosity, bounds_marker, yielding_marker, visc_options,
      1 /* remove_bound_min */,
      0 /* !remove_bound_max */,
      1 /* remove_yielding */,
      0 /* !remove_upper_mantle */,
      0 /* !remove_lower_mantle */);
  rhea_viscosity_filter_where_max (derivative, bounds_marker,
                                   0 /* !invert_filter */);

  /* remove previous scaling and multiply by value corresponding to deriv. */
  ymir_vec_scale (scaling_deriv/scaling_prev, derivative);

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
                                      const int process_upper_mantle,
                                      const int process_lower_mantle)
{
  double              scaling_prev, scaling_deriv;
  double              scaling_upper_mantle, scaling_lower_mantle;

  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (derivative));
  RHEA_ASSERT (process_upper_mantle || process_lower_mantle);

  /* set derivative */
  rhea_viscosity_param_derivative_init (
      derivative, viscosity, bounds_marker, yielding_marker, visc_options,
      1 /* remove_bound_min */,
      1 /* remove_bound_max */,
      1 /* remove_yielding */,
      !process_upper_mantle,
      !process_lower_mantle);

  /* remove previous scaling and multiply by value corresponding to deriv. */
  if (process_upper_mantle) { /* if upper mantle */
    scaling_prev = visc_options->upper_mantle_scaling;
    scaling_deriv =
      rhea_inversion_param_convert_from_model_pos_deriv (scaling_prev);
    scaling_upper_mantle = scaling_deriv/scaling_prev;
  }
  else {
    scaling_upper_mantle = NAN;
  }
  if (process_lower_mantle) { /* if lower mantle */
    scaling_prev = visc_options->lower_mantle_scaling;
    scaling_deriv =
      rhea_inversion_param_convert_from_model_pos_deriv (scaling_prev);
    scaling_lower_mantle = scaling_deriv/scaling_prev;
  }
  else {
    scaling_lower_mantle = NAN;
  }
  rhea_viscosity_param_derivative_scale_mantle (
      derivative, scaling_upper_mantle, scaling_lower_mantle, visc_options);

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
                                      const int process_upper_mantle,
                                      const int process_lower_mantle)
{
  double              scaling_prev, scaling_deriv;
  double              scaling_upper_mantle, scaling_lower_mantle;
  ymir_vec_t         *scal = ymir_vec_template (derivative);

  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (derivative));
  RHEA_ASSERT (process_upper_mantle || process_lower_mantle);

  /* initialize derivative */
  rhea_viscosity_param_derivative_init (
      derivative, viscosity, bounds_marker, yielding_marker, visc_options,
      1 /* remove_bound_min */,
      1 /* remove_bound_max */,
      1 /* remove_yielding */,
      !process_upper_mantle,
      !process_lower_mantle);

  /* multiply derivative by temperature difference term */
  ymir_interp_vec (temperature, scal);
  rhea_temperature_bound (scal);
  ymir_vec_scale_shift (-1.0, RHEA_TEMPERATURE_NEUTRAL_VALUE, scal);
  ymir_vec_multiply_in (scal, derivative);
  ymir_vec_destroy (scal);

  /* multiply derivative by activation energy */
  if (process_upper_mantle) { /* if upper mantle */
    scaling_prev = visc_options->upper_mantle_arrhenius_activation_energy;
    scaling_deriv =
      rhea_inversion_param_convert_from_model_pos_deriv (scaling_prev);
    scaling_upper_mantle = scaling_deriv;
  }
  else {
    scaling_upper_mantle = NAN;
  }
  if (process_lower_mantle) { /* if lower mantle */
    scaling_prev = visc_options->lower_mantle_arrhenius_activation_energy;
    scaling_deriv =
      rhea_inversion_param_convert_from_model_pos_deriv (scaling_prev);
    scaling_lower_mantle = scaling_deriv;
  }
  else {
    scaling_lower_mantle = NAN;
  }
  rhea_viscosity_param_derivative_scale_mantle (
      derivative, scaling_upper_mantle, scaling_lower_mantle, visc_options);

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
  const double        stress_exp_deriv =
    rhea_inversion_param_convert_from_model_n_deriv (stress_exp);
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
      0 /* !remove_upper_mantle */,
      1 /* remove_lower_mantle */);

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
  scal = - SC_MAX (0.0, stress_exp_deriv) / (stress_exp*stress_exp);
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
  const double        scaling_prev = visc_options->yield_strength;
  const double        scaling_deriv =
    rhea_inversion_param_convert_from_model_pos_deriv (scaling_prev);

  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (derivative));

  /* set derivative */
  rhea_viscosity_param_derivative_init (
      derivative, viscosity, bounds_marker, yielding_marker, visc_options,
      1 /* remove_bound_min */,
      0 /* !remove_bound_max */,
      0 /* !remove_yielding */,
      0 /* !remove_upper_mantle */,
      1 /* remove_lower_mantle */);
  rhea_viscosity_filter_where_yielding (derivative, yielding_marker,
                                        0 /* !invert_filter */);

  /* remove previous scaling and multiply by value corresponding to deriv. */
  ymir_vec_scale (scaling_deriv/scaling_prev, derivative);

  /* check output */
  RHEA_ASSERT (sc_dmatrix_is_valid (derivative->dataown));
}

/**
 * Computes derivative w.r.t. the weak factor in the interior of a weak zone:
 *   weak_int = exp(-param_weak_int^2)
 *   weak = 1 - (1 - weak_int) * weak_indicator
 *   deriv = -2*param_weak_int * weak_indicator * visc/weak
 */
static void
rhea_viscosity_param_derivative_weak_factor_interior (
                                        ymir_vec_t *derivative,
                                        rhea_weakzone_label_t weak_label,
                                        ymir_vec_t *viscosity,
                                        ymir_vec_t *bounds_marker,
                                        ymir_vec_t *yielding_marker,
                                        ymir_vec_t *weakzone,
                                        rhea_weakzone_options_t *weak_options,
                                        rhea_viscosity_options_t *visc_options)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (weakzone);
  ymir_vec_t         *weak_indicator = rhea_weakzone_new (ymir_mesh);
  double              weak_factor, scaling_deriv;

  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (derivative));
  RHEA_ASSERT (rhea_weakzone_label_is_valid (weak_label));

  /* get weak zone factor */
  weak_factor = rhea_weakzone_lookup_factor_interior ((int) weak_label,
                                                      weak_options);

  /* set scaling factor for the derivative */
  scaling_deriv =
    rhea_inversion_param_convert_from_model_weak_deriv (weak_factor);

  /* init derivative */
  rhea_viscosity_param_derivative_init (
      derivative, viscosity, bounds_marker, yielding_marker, visc_options,
      1 /* remove_bound_min */,
      0 /* !remove_bound_max */,
      1 /* remove_yielding */,
      0 /* !remove_upper_mantle */,
      1 /* remove_lower_mantle */);

  /* get indicator for weak zones */
  rhea_weakzone_compute_indicator (weak_indicator, (int) weak_label,
                                   weak_options);

  /* compute derivative */
  ymir_vec_scale (scaling_deriv, weak_indicator);
  ymir_vec_multiply_in (weak_indicator, derivative);
  ymir_vec_divide_in (weakzone, derivative);

  /* destroy */
  rhea_weakzone_destroy (weak_indicator);

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
                            ymir_vec_t *velocity,
                            rhea_viscosity_options_t *visc_options)
{
  /* initialize to zero */
  ymir_vec_set_zero (derivative);

  /* compute derivative */
  switch (derivative_type) {
  case RHEA_VISCOSITY_PARAM_DERIVATIVE_MIN:
    if (rhea_viscosity_restrict_min (visc_options)) {
      rhea_viscosity_param_derivative_min (
          derivative, viscosity, bounds_marker, yielding_marker, visc_options);
    }
    break;

  case RHEA_VISCOSITY_PARAM_DERIVATIVE_MAX:
    if (rhea_viscosity_restrict_max (visc_options)) {
      rhea_viscosity_param_derivative_max (
          derivative, viscosity, bounds_marker, yielding_marker, visc_options);
    }
    break;

  case RHEA_VISCOSITY_PARAM_DERIVATIVE_UPPER_MANTLE_SCALING:
    rhea_viscosity_param_derivative_scal (
        derivative, viscosity, bounds_marker, yielding_marker, visc_options,
        1 /* upper mantle */, 0 /* not lower mantle */);
    break;

  case RHEA_VISCOSITY_PARAM_DERIVATIVE_UPPER_MANTLE_ACTIVATION_ENERGY:
    if (rhea_viscosity_has_arrhenius (visc_options)) {
      rhea_viscosity_param_derivative_arrh (
          derivative, viscosity, bounds_marker, yielding_marker, temperature,
          visc_options, 1 /* upper mantle */, 0 /* not lower mantle */);
    }
    break;

  case RHEA_VISCOSITY_PARAM_DERIVATIVE_LOWER_MANTLE_SCALING:
    rhea_viscosity_param_derivative_scal (
        derivative, viscosity, bounds_marker, yielding_marker, visc_options,
        0 /* not upper mantle */, 1 /* lower mantle */);
    break;

  case RHEA_VISCOSITY_PARAM_DERIVATIVE_LOWER_MANTLE_ACTIVATION_ENERGY:
    if (rhea_viscosity_has_arrhenius (visc_options)) {
      rhea_viscosity_param_derivative_arrh (
          derivative, viscosity, bounds_marker, yielding_marker, temperature,
          visc_options, 0 /* not upper mantle */, 1 /* lower mantle */);
    }
    break;

  case RHEA_VISCOSITY_PARAM_DERIVATIVE_STRESS_EXPONENT:
    if (rhea_viscosity_has_strain_rate_weakening (visc_options)) {
      rhea_viscosity_param_derivative_stress_exp (
          derivative, viscosity, bounds_marker, yielding_marker, velocity,
          visc_options);
    }
    break;

  case RHEA_VISCOSITY_PARAM_DERIVATIVE_YIELD_STRENGTH:
    if (rhea_viscosity_has_yielding (visc_options)) {
      rhea_viscosity_param_derivative_yield_strength (
          derivative, viscosity, bounds_marker, yielding_marker, visc_options);
    }
    break;

  default: /* unknown derivative type */
    RHEA_ABORT_NOT_REACHED ();
  }

#if (0 < RHEA_VISCOSITY_PARAM_DERIVATIVE_VTK)
  {
    char                path[BUFSIZ];

    snprintf (path, BUFSIZ, "%s_type%i", __func__, derivative_type);
    ymir_vtk_write (ymir_vec_get_mesh (derivative), path,
                    derivative, "derivative",
                    temperature, "temperature",
                    viscosity, "viscosity", NULL);
  }
#endif
}

void
rhea_viscosity_param_derivative_weakzone (
                            ymir_vec_t *derivative,
                            rhea_viscosity_param_derivative_t derivative_type,
                            rhea_weakzone_label_t weak_label,
                            ymir_vec_t *viscosity,
                            ymir_vec_t *bounds_marker,
                            ymir_vec_t *yielding_marker,
                            ymir_vec_t *weakzone,
                            rhea_weakzone_options_t *weak_options,
                            rhea_viscosity_options_t *visc_options)
{
  /* initialize to zero */
  ymir_vec_set_zero (derivative);

  /* exit if nothing to do */
  if (!rhea_weakzone_exists (weak_options)) {
    return;
  }

  /* compute derivative */
  switch (derivative_type) {
  case RHEA_VISCOSITY_PARAM_DERIVATIVE_WEAK_THICKNESS:
    RHEA_ABORT_NOT_REACHED (); //TODO
    break;
  case RHEA_VISCOSITY_PARAM_DERIVATIVE_WEAK_THICKNESS_CONST:
    RHEA_ABORT_NOT_REACHED (); //TODO
    break;
  case RHEA_VISCOSITY_PARAM_DERIVATIVE_WEAK_FACTOR_INTERIOR:
    rhea_viscosity_param_derivative_weak_factor_interior (
        derivative, weak_label, viscosity, bounds_marker, yielding_marker,
        weakzone, weak_options, visc_options);
    break;
  default: /* unknown derivative type */
    RHEA_ABORT_NOT_REACHED ();
  }

#if (0 < RHEA_VISCOSITY_PARAM_DERIVATIVE_VTK)
  {
    char                path[BUFSIZ];

    snprintf (path, BUFSIZ, "%s_type%i", __func__, derivative_type);
    ymir_vtk_write (ymir_vec_get_mesh (derivative), path,
                    derivative, "derivative",
                    viscosity, "viscosity", NULL);
  }
#endif
}
