#include <rhea_viscosity_param_derivative.h>
#include <rhea_base.h>
#include <rhea_math.h>
#include <rhea_inversion_param.h>
#include <rhea_temperature.h>
#include <rhea_strainrate.h>
#include <ymir_interp_vec.h>

#define RHEA_VISCOSITY_PARAM_DERIVATIVE_VTK (0)

#if (0 < RHEA_VISCOSITY_PARAM_DERIVATIVE_VTK)
# include <ymir_vtk.h>
#endif

/* types of markers for different viscosity regiemes */
typedef enum
{
  RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_MIN,
  RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_YLD,
  RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_WKZ,
  RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_MAX,
  RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_SRW,
  RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_LIN_UM,
  RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_LIN_LM
}
rhea_viscosity_param_derivative_filter_t;

/**
 * Filters upper mantle, i.e, removes lower mantle.
 */
static void
rhea_viscosity_param_derivative_filter_upper_mantle (
                                        ymir_vec_t *derivative,
                                        rhea_domain_options_t *domain_options)
{
  ymir_vec_t         *filter = ymir_vec_template (derivative);

  rhea_viscosity_stats_filter_upper_mantle (filter, domain_options);
  ymir_vec_multiply_in (filter, derivative);
  ymir_vec_destroy (filter);
}

/**
 * Filters lower mantle, i.e, removes upper mantle.
 */
static void
rhea_viscosity_param_derivative_filter_lower_mantle (
                                        ymir_vec_t *derivative,
                                        rhea_domain_options_t *domain_options)
{
  ymir_vec_t         *filter = ymir_vec_template (derivative);

  rhea_viscosity_stats_filter_lower_mantle (filter, domain_options);
  ymir_vec_multiply_in (filter, derivative);
  ymir_vec_destroy (filter);
}

/**
 * Computes the derivative of the smooth application of the upper viscosity
 * bound.
 */
static void
rhea_viscosity_param_derivative_visc_max_smooth (
                                        ymir_vec_t *derivative,
                                        ymir_vec_t *weakzone,
                                        rhea_viscosity_options_t *visc_options)
{
  const double        visc_max = visc_options->max;
  const double        smoothness_param = visc_options->max_smoothness_param;
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (derivative);
  const ymir_locidx_t n_elements = ymir_mesh->cnodes->K;
  const int           n_nodes = ymir_mesh->cnodes->N;
  ymir_locidx_t       elid;
  int                 nodeid;

  for (elid = 0; elid < n_elements; elid++) {
    for (nodeid = 0; nodeid < n_nodes; nodeid++) {
      const double        w = (weakzone == NULL ? 1.0 :
                               *ymir_dvec_index (weakzone, elid, nodeid, 0));
      double             *d = ymir_dvec_index (derivative, elid, nodeid, 0);

      *d *= 1.0/w;
      *d = rhea_math_smin_gpm_dx_impl_nondim (*d, visc_max, smoothness_param,
                                              visc_max);
      *d *= w;
    }
  }
}

/**
 * Initializes the derivative vector by reversing the operations of the
 * viscosity computation.
 *
 *           yield    visc_max  lower mantle lin visc
 *             ^            ^    ^
 *             |            |    |
 * visc_min -> * -> weak -> * -> * -> strain -> upper mantle lin visc
 */
static void
rhea_viscosity_param_derivative_init (
                        ymir_vec_t *derivative,
                        const rhea_viscosity_param_derivative_filter_t filter,
                        ymir_vec_t *viscosity,
                        ymir_vec_t *marker,
                        ymir_vec_t *weakzone,
                        rhea_viscosity_options_t *visc_options)
{
  /* check input */
  RHEA_ASSERT (NULL != viscosity);
  RHEA_ASSERT (NULL != marker);
  RHEA_ASSERT (NULL != visc_options);
  RHEA_ASSERT (visc_options->model == RHEA_VISCOSITY_MODEL_UWYL_LCUT_UCUT ||
               visc_options->model == RHEA_VISCOSITY_MODEL_UWYL_LADD_UCUT ||
               visc_options->model == RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT ||
               visc_options->model == RHEA_VISCOSITY_MODEL_UWYL_LADD_USMOOTH);

  /* initialize derivative */
  ymir_vec_copy (viscosity, derivative);

  /* lower viscosity bound */
  if (filter == RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_MIN) {
    if (rhea_viscosity_restrict_min (visc_options)) {
      rhea_viscosity_marker_filter_min (derivative, marker, 0 /* !invert */);
    }
    else { /* otherwise set to zero */
      ymir_vec_set_zero (derivative);
    }
    return;
  }
  else if (rhea_viscosity_restrict_min (visc_options)) { /* remove min bound */
    switch (visc_options->model) {
    case RHEA_VISCOSITY_MODEL_UWYL_LCUT_UCUT:
      rhea_viscosity_marker_filter_min (derivative, marker, 1 /* invert */);
      break;
    case RHEA_VISCOSITY_MODEL_UWYL_LADD_UCUT:
    case RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT:
    case RHEA_VISCOSITY_MODEL_UWYL_LADD_USMOOTH:
      ymir_vec_shift (-visc_options->min, derivative);
      ymir_vec_bound_min (derivative, 0.0);
      break;
    default: /* unknown viscosity model */
      RHEA_ABORT_NOT_REACHED ();
    }
  }

  /* yielding */
  if (filter == RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_YLD) {
    if (rhea_viscosity_has_yielding (visc_options)) {
      rhea_viscosity_marker_filter_yielding (derivative, marker,
                                             0 /* !invert */);
    }
    else { /* otherwise set to zero */
      ymir_vec_set_zero (derivative);
    }
    return;
  }
  else if (rhea_viscosity_has_yielding (visc_options)) { /* remove yielding */
    rhea_viscosity_marker_filter_yielding (derivative, marker,
                                           1 /* invert */);
  }

  /* weak zone */
  if (filter == RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_WKZ) {
    return;
  }

  /* upper viscosity bound */
  if (filter == RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_MAX) {
    RHEA_ASSERT (visc_options->model != RHEA_VISCOSITY_MODEL_UWYL_LADD_USMOOTH);
    //TODO implement for *_USMOOTH
    if (rhea_viscosity_restrict_max (visc_options)) {
      rhea_viscosity_marker_filter_max (derivative, marker, 0 /* !invert */);
    }
    else { /* otherwise set to zero */
      ymir_vec_set_zero (derivative);
    }
    return;
  }
  else if (rhea_viscosity_restrict_max (visc_options)) { /* remove max bound */
    switch (visc_options->model) {
    case RHEA_VISCOSITY_MODEL_UWYL_LCUT_UCUT:
    case RHEA_VISCOSITY_MODEL_UWYL_LADD_UCUT:
    case RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT:
      rhea_viscosity_marker_filter_max (derivative, marker, 1 /* invert */);
      break;
    case RHEA_VISCOSITY_MODEL_UWYL_LADD_USMOOTH:
      rhea_viscosity_marker_filter_max (derivative, marker, 1 /* invert */);
      rhea_viscosity_param_derivative_visc_max_smooth (
          derivative, weakzone, visc_options);
      break;
    default: /* unknown viscosity model */
      RHEA_ABORT_NOT_REACHED ();
    }
  }

  /* strain rate weakening */
  if (filter == RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_SRW) {
    if (rhea_viscosity_has_strain_rate_weakening (visc_options)) {
      rhea_viscosity_param_derivative_filter_upper_mantle (
          derivative, visc_options->domain_options);
    }
    else { /* otherwise set to zero */
      ymir_vec_set_zero (derivative);
    }
    return;
  }

  /* upper mantle */
  if (filter == RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_LIN_UM) {
    rhea_viscosity_param_derivative_filter_upper_mantle (
        derivative, visc_options->domain_options);
    return;
  }

  /* lower mantle */
  if (filter == RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_LIN_LM) {
    rhea_viscosity_param_derivative_filter_lower_mantle (
        derivative, visc_options->domain_options);
    return;
  }

  /* unknown filter */
  RHEA_ABORT_NOT_REACHED ();
}

/**
 * Computes derivative w.r.t. the lower viscosity bound.
 */
static void
rhea_viscosity_param_derivative_min (ymir_vec_t *derivative,
                                     ymir_vec_t *viscosity,
                                     ymir_vec_t *marker,
                                     rhea_viscosity_options_t *visc_options)
{
  const double        scaling_curr = visc_options->min;
  const double        scaling_deriv =
    rhea_inversion_param_derivative_min_max (scaling_curr);

  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (derivative));

  /* set derivative */
  switch (visc_options->model) {
  case RHEA_VISCOSITY_MODEL_UWYL_LCUT_UCUT:
    rhea_viscosity_param_derivative_init (
        derivative, RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_MIN,
        viscosity, marker, NULL /* weakzone */, visc_options);
    ymir_vec_scale (scaling_deriv/scaling_curr, derivative);
    break;
  case RHEA_VISCOSITY_MODEL_UWYL_LADD_UCUT:
  case RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT:
  case RHEA_VISCOSITY_MODEL_UWYL_LADD_USMOOTH:
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
                                     ymir_vec_t *marker,
                                     ymir_vec_t *weakzone,
                                     rhea_viscosity_options_t *visc_options)
{
  const double        scaling_curr = visc_options->max;
  const double        scaling_deriv =
    rhea_inversion_param_derivative_min_max (scaling_curr);

  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (derivative));

  /* init derivative */
  rhea_viscosity_param_derivative_init (
      derivative, RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_MAX,
      viscosity, marker, weakzone, visc_options);

  /* remove previous scaling and multiply by value corresponding to deriv. */
  ymir_vec_scale (scaling_deriv/scaling_curr, derivative);

  /* check output */
  RHEA_ASSERT (sc_dmatrix_is_valid (derivative->dataown));
}

/**
 * Computes derivative w.r.t. the scaling factor of upper/lower mantle:
 *   c = X(param_c)
 *   deriv = (d X(param_c)) * 1/(c*n) * visc
 */
static void
rhea_viscosity_param_derivative_upper_mantle_scaling (
                                        ymir_vec_t *derivative,
                                        ymir_vec_t *viscosity,
                                        ymir_vec_t *marker,
                                        ymir_vec_t *weakzone,
                                        rhea_viscosity_options_t *visc_options)
{
  double              scaling_curr, scaling_deriv, srw_exp;

  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (derivative));

  /* init derivative */
  rhea_viscosity_param_derivative_init (
      derivative, RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_LIN_UM,
      viscosity, marker, weakzone, visc_options);

  /* remove previous scaling and multiply by value corresponding to deriv. */
  scaling_curr = visc_options->upper_mantle_scaling;
  scaling_deriv = rhea_inversion_param_derivative_scal (scaling_curr);
  srw_exp = rhea_viscosity_get_strain_rate_weakening_exp (
      visc_options, 1 /* in upper mantle */);
  ymir_vec_scale (scaling_deriv/scaling_curr*srw_exp, derivative);

  /* check output */
  RHEA_ASSERT (sc_dmatrix_is_valid (derivative->dataown));
}

static void
rhea_viscosity_param_derivative_lower_mantle_scaling (
                                        ymir_vec_t *derivative,
                                        ymir_vec_t *viscosity,
                                        ymir_vec_t *marker,
                                        ymir_vec_t *weakzone,
                                        rhea_viscosity_options_t *visc_options)
{
  double              scaling_curr, scaling_deriv, srw_exp;

  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (derivative));

  /* init derivative */
  rhea_viscosity_param_derivative_init (
      derivative, RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_LIN_LM,
      viscosity, marker, weakzone, visc_options);

  /* remove previous scaling and multiply by value corresponding to deriv. */
  scaling_curr = visc_options->lower_mantle_scaling;
  scaling_deriv = rhea_inversion_param_derivative_scal (scaling_curr);
  srw_exp = rhea_viscosity_get_strain_rate_weakening_exp (
      visc_options, 0 /* in lower mantle */);
  ymir_vec_scale (scaling_deriv/scaling_curr*srw_exp, derivative);

  /* check output */
  RHEA_ASSERT (sc_dmatrix_is_valid (derivative->dataown));
}

/**
 * Computes derivative w.r.t. the activation energy in Arrhenius relationship:
 *   E_a = X(param_E_a)
 *   deriv = (d X(param_E_a)) * (0.5 - T)/n * visc
 */
static void
rhea_viscosity_param_derivative_scale_by_temp_diff (
                                        ymir_vec_t *derivative,
                                        ymir_vec_t *temperature,
                                        rhea_viscosity_options_t *visc_options)
{
  const double        temp_neutral = visc_options->temp_options->neutral;
  ymir_vec_t         *scaling_vec = ymir_vec_template (derivative);

  ymir_interp_vec (temperature, scaling_vec);
  rhea_temperature_bound (scaling_vec);
  ymir_vec_scale_shift (-1.0, temp_neutral, scaling_vec);
  ymir_vec_multiply_in (scaling_vec, derivative);

  ymir_vec_destroy (scaling_vec);
}

static void
rhea_viscosity_param_derivative_upper_mantle_activation_energy (
                                        ymir_vec_t *derivative,
                                        ymir_vec_t *viscosity,
                                        ymir_vec_t *marker,
                                        ymir_vec_t *temperature,
                                        ymir_vec_t *weakzone,
                                        rhea_viscosity_options_t *visc_options)
{
  double              scaling_curr, scaling_deriv, srw_exp;

  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (derivative));

  /* init derivative */
  rhea_viscosity_param_derivative_init (
      derivative, RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_LIN_UM,
      viscosity, marker, weakzone, visc_options);

  /* multiply by temperature difference term */
  rhea_viscosity_param_derivative_scale_by_temp_diff (derivative, temperature,
                                                      visc_options);

  /* multiply by the derivative of the activation energy */
  scaling_curr = visc_options->upper_mantle_arrhenius_activation_energy;
  scaling_deriv = rhea_inversion_param_derivative_arrh (scaling_curr);
  srw_exp = rhea_viscosity_get_strain_rate_weakening_exp (
      visc_options, 1 /* in upper mantle */);
  ymir_vec_scale (scaling_deriv*srw_exp, derivative);

  /* check output */
  RHEA_ASSERT (sc_dmatrix_is_valid (derivative->dataown));
}

static void
rhea_viscosity_param_derivative_lower_mantle_activation_energy (
                                        ymir_vec_t *derivative,
                                        ymir_vec_t *viscosity,
                                        ymir_vec_t *marker,
                                        ymir_vec_t *temperature,
                                        ymir_vec_t *weakzone,
                                        rhea_viscosity_options_t *visc_options)
{
  double              scaling_curr, scaling_deriv, srw_exp;

  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (derivative));

  /* init derivative */
  rhea_viscosity_param_derivative_init (
      derivative, RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_LIN_LM,
      viscosity, marker, weakzone, visc_options);

  /* multiply by temperature difference term */
  rhea_viscosity_param_derivative_scale_by_temp_diff (derivative, temperature,
                                                      visc_options);

  /* multiply by the derivative of the activation energy */
  scaling_curr = visc_options->lower_mantle_arrhenius_activation_energy;
  scaling_deriv = rhea_inversion_param_derivative_arrh (scaling_curr);
  srw_exp = rhea_viscosity_get_strain_rate_weakening_exp (
      visc_options, 0 /* in lower mantle */);
  ymir_vec_scale (scaling_deriv*srw_exp, derivative);

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
 *   n = X(param_n)
 *   deriv = (d X(param_n)) * (- log(c*a(T)*(strainrate-shift))/n^2 ) * visc
 */
static void
rhea_viscosity_param_derivative_stress_exp (
                                        ymir_vec_t *derivative,
                                        ymir_vec_t *viscosity,
                                        ymir_vec_t *marker,
                                        ymir_vec_t *temperature,
                                        ymir_vec_t *weakzone,
                                        ymir_vec_t *velocity,
                                        rhea_viscosity_options_t *visc_options)
{
  const double        stress_exp = visc_options->stress_exponent;
  const double        stress_exp_deriv =
    rhea_inversion_param_derivative_n (stress_exp);
  const rhea_viscosity_t  visc_type = visc_options->type;
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (derivative);
  ymir_vec_t         *scal_vec, *visc_lin;
  double              scal_val;

  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (derivative));

  /* init derivative */
  rhea_viscosity_param_derivative_init (
      derivative, RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_SRW,
      viscosity, marker, weakzone, visc_options);

  /* compute the spatially variing scaling */
  scal_vec = rhea_strainrate_2inv_new (ymir_mesh);
  visc_lin = rhea_viscosity_new (ymir_mesh);
  switch (visc_options->model) {
  case RHEA_VISCOSITY_MODEL_UWYL_LCUT_UCUT:
  case RHEA_VISCOSITY_MODEL_UWYL_LADD_UCUT:
  case RHEA_VISCOSITY_MODEL_UWYL_LADD_USMOOTH:
    rhea_strainrate_compute_sqrt_of_2inv (scal_vec, velocity);
    visc_options->type = RHEA_VISCOSITY_LINEAR;
    rhea_viscosity_compute (
        /* out: */ visc_lin, NULL, NULL,
        /* in:  */ temperature, NULL /* w/o weakzone */, NULL /* velocity */,
        visc_options);
    visc_options->type = visc_type;
    ymir_vec_multiply_in (visc_lin, scal_vec);
    break;
  case RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT:
    ymir_vec_set_value (scal_vec, NAN);
    RHEA_ABORT_NOT_REACHED (); //TODO implement this; need a shift apply fnc.
    break;
  default: /* unknown viscosity model */
    RHEA_ABORT_NOT_REACHED ();
  }
  ymir_dvec_set_function (scal_vec, rhea_viscosity_param_derivative_log_fn,
                          NULL);
  rhea_viscosity_destroy (visc_lin);

  /* multiply-in the spatially variing scaling */
  ymir_vec_multiply_in (scal_vec, derivative);
  rhea_strainrate_2inv_destroy (scal_vec);

  /* multiply by constant scaling */
  scal_val = SC_MAX (0.0, stress_exp_deriv) / (-stress_exp*stress_exp);
  ymir_vec_scale (scal_val, derivative);

  /* check output */
  RHEA_ASSERT (sc_dmatrix_is_valid (derivative->dataown));
}

/**
 * Computes derivative w.r.t. the yield strength:
 *   Y = X(param_Y)
 *   deriv = (d X(param_Y) * visc / Y
 */
static void
rhea_viscosity_param_derivative_yield_strength (
                                        ymir_vec_t *derivative,
                                        ymir_vec_t *viscosity,
                                        ymir_vec_t *marker,
                                        rhea_viscosity_options_t *visc_options)
{
  const double        scaling_curr = visc_options->yield_strength;
  const double        scaling_deriv =
    rhea_inversion_param_derivative_yield (scaling_curr);

  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (derivative));

  /* init derivative */
  rhea_viscosity_param_derivative_init (
      derivative, RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_YLD,
      viscosity, marker, NULL /* weakzone */, visc_options);

  /* remove previous scaling and multiply by value corresponding to deriv. */
  ymir_vec_scale (scaling_deriv/scaling_curr, derivative);

  /* check output */
  RHEA_ASSERT (sc_dmatrix_is_valid (derivative->dataown));
}

/**
 * Computes derivative w.r.t. the weak factor in the interior of a weak zone:
 *   weak_int = X(param_weak_int)
 *   weak = 1 - (1 - weak_int) * weak_indicator
 *   deriv = (d X(param_weak_int)) * weak_indicator * visc / weak
 */
static void
rhea_viscosity_param_derivative_weak_factor_interior (
                                        ymir_vec_t *derivative,
                                        rhea_weakzone_label_t weak_label,
                                        ymir_vec_t *viscosity,
                                        ymir_vec_t *marker,
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
  scaling_deriv = rhea_inversion_param_derivative_weak (weak_factor);

  /* init derivative */
  rhea_viscosity_param_derivative_init (
      derivative, RHEA_VISCOSITY_PARAM_DERIVATIVE_FILTER_WKZ,
      viscosity, marker, weakzone, visc_options);

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
                            ymir_vec_t *marker,
                            ymir_vec_t *temperature,
                            ymir_vec_t *weakzone,
                            ymir_vec_t *velocity,
                            rhea_viscosity_options_t *visc_options)
{
  /* initialize */
  ymir_vec_set_zero (derivative);

  /* compute derivative */
  switch (derivative_type) {
  case RHEA_VISCOSITY_PARAM_DERIVATIVE_MIN:
    if (rhea_viscosity_restrict_min (visc_options)) {
      rhea_viscosity_param_derivative_min (
          derivative, viscosity, marker, visc_options);
    }
    break;

  case RHEA_VISCOSITY_PARAM_DERIVATIVE_MAX:
    if (rhea_viscosity_restrict_max (visc_options)) {
      rhea_viscosity_param_derivative_max (
          derivative, viscosity, marker, weakzone, visc_options);
    }
    break;

  case RHEA_VISCOSITY_PARAM_DERIVATIVE_UPPER_MANTLE_SCALING:
    rhea_viscosity_param_derivative_upper_mantle_scaling (
        derivative, viscosity, marker, weakzone, visc_options);
    break;

  case RHEA_VISCOSITY_PARAM_DERIVATIVE_UPPER_MANTLE_ACTIVATION_ENERGY:
    if (rhea_viscosity_has_arrhenius (visc_options)) {
      rhea_viscosity_param_derivative_upper_mantle_activation_energy (
          derivative, viscosity, marker, temperature, weakzone, visc_options);
    }
    break;

  case RHEA_VISCOSITY_PARAM_DERIVATIVE_LOWER_MANTLE_SCALING:
    rhea_viscosity_param_derivative_lower_mantle_scaling (
        derivative, viscosity, marker, weakzone, visc_options);
    break;

  case RHEA_VISCOSITY_PARAM_DERIVATIVE_LOWER_MANTLE_ACTIVATION_ENERGY:
    if (rhea_viscosity_has_arrhenius (visc_options)) {
      rhea_viscosity_param_derivative_lower_mantle_activation_energy (
          derivative, viscosity, marker, temperature, weakzone, visc_options);
    }
    break;

  case RHEA_VISCOSITY_PARAM_DERIVATIVE_STRESS_EXPONENT:
    if (rhea_viscosity_has_strain_rate_weakening (visc_options)) {
      rhea_viscosity_param_derivative_stress_exp (
          derivative, viscosity, marker, temperature, weakzone, velocity,
          visc_options);
    }
    break;

  case RHEA_VISCOSITY_PARAM_DERIVATIVE_YIELD_STRENGTH:
    if (rhea_viscosity_has_yielding (visc_options)) {
      rhea_viscosity_param_derivative_yield_strength (
          derivative, viscosity, marker, visc_options);
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
                    viscosity, "viscosity",
                    marker, "marker",
                    temperature, "temperature", NULL);
  }
#endif
}

void
rhea_viscosity_param_derivative_weakzone (
                            ymir_vec_t *derivative,
                            rhea_viscosity_param_derivative_t derivative_type,
                            rhea_weakzone_label_t weak_label,
                            ymir_vec_t *viscosity,
                            ymir_vec_t *marker,
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
        derivative, weak_label, viscosity, marker, weakzone,
        weak_options, visc_options);
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
                    viscosity, "viscosity",
                    marker, "marker",
                    weakzone, "weakzone", NULL);
  }
#endif
}
