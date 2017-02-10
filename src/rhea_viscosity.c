/*
 */

#include <rhea_viscosity.h>
#include <rhea_base.h>
#include <rhea_temperature.h>
#include <rhea_weakzone.h>
#include <rhea_velocity.h>
#include <rhea_strainrate.h>
#include <ymir_vec_getset.h>
#include <ymir_stress_op.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_VISCOSITY_DEFAULT_TYPE (RHEA_VISCOSITY_LINEAR)
#define RHEA_VISCOSITY_DEFAULT_TYPE_INIT_NONLINEAR \
  (RHEA_VISCOSITY_INIT_NONLINEAR_DEFAULT)
#define RHEA_VISCOSITY_DEFAULT_MODEL_NAME "UWYL_LADD_USHIFT"
#define RHEA_VISCOSITY_DEFAULT_MIN (1.0e-2)
#define RHEA_VISCOSITY_DEFAULT_MAX (1.0e+4)
#define RHEA_VISCOSITY_DEFAULT_UPPER_MANTLE_SCALING (4.0e+3)
#define RHEA_VISCOSITY_DEFAULT_UPPER_MANTLE_ACTIVATION_ENERGY (17.5)
#define RHEA_VISCOSITY_DEFAULT_LOWER_MANTLE_SCALING (4.0e+5)
#define RHEA_VISCOSITY_DEFAULT_LOWER_MANTLE_ACTIVATION_ENERGY (17.5)
#define RHEA_VISCOSITY_DEFAULT_STRESS_EXPONENT (3.0)
#define RHEA_VISCOSITY_DEFAULT_YIELD_STRESS (-1.0)
#define RHEA_VISCOSITY_DEFAULT_YIELDING_REGULARIZATION (0.0)

/* initialize options */
int                 rhea_viscosity_type = RHEA_VISCOSITY_DEFAULT_TYPE;
int                 rhea_viscosity_type_init_nonlinear =
                      RHEA_VISCOSITY_DEFAULT_TYPE_INIT_NONLINEAR;
char               *rhea_viscosity_model_name =
                      RHEA_VISCOSITY_DEFAULT_MODEL_NAME;
double              rhea_viscosity_min = RHEA_VISCOSITY_DEFAULT_MIN;
double              rhea_viscosity_max = RHEA_VISCOSITY_DEFAULT_MAX;
double              rhea_viscosity_upper_mantle_scaling =
  RHEA_VISCOSITY_DEFAULT_UPPER_MANTLE_SCALING;
double              rhea_viscosity_upper_mantle_activation_energy =
  RHEA_VISCOSITY_DEFAULT_UPPER_MANTLE_ACTIVATION_ENERGY;
double              rhea_viscosity_lower_mantle_scaling =
  RHEA_VISCOSITY_DEFAULT_LOWER_MANTLE_SCALING;
double              rhea_viscosity_lower_mantle_activation_energy =
  RHEA_VISCOSITY_DEFAULT_LOWER_MANTLE_ACTIVATION_ENERGY;
double              rhea_viscosity_stress_exponent =
  RHEA_VISCOSITY_DEFAULT_STRESS_EXPONENT;
double              rhea_viscosity_yield_stress =
  RHEA_VISCOSITY_DEFAULT_YIELD_STRESS;
double              rhea_viscosity_yielding_regularization =
  RHEA_VISCOSITY_DEFAULT_YIELDING_REGULARIZATION;

void
rhea_viscosity_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Viscosity";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_I, "type", '\0',
    &(rhea_viscosity_type), RHEA_VISCOSITY_DEFAULT_TYPE,
    "Viscosity type: 0: constant; 1: linear; 2: nonlinear",
  YMIR_OPTIONS_I, "init-nonlinear-type", '\0',
    &(rhea_viscosity_type_init_nonlinear),
    RHEA_VISCOSITY_DEFAULT_TYPE_INIT_NONLINEAR,
    "Viscosity type for initial viscosity when solving a nonlinear problem",
  YMIR_OPTIONS_S, "model", '\0',
    &(rhea_viscosity_model_name), RHEA_VISCOSITY_DEFAULT_MODEL_NAME,
    "Viscosity model name, e.g., 'UWYL'",

  YMIR_OPTIONS_D, "min", '\0',
    &(rhea_viscosity_min), RHEA_VISCOSITY_DEFAULT_MIN,
    "Lower bound for viscosity",
  YMIR_OPTIONS_D, "max", '\0',
    &(rhea_viscosity_max), RHEA_VISCOSITY_DEFAULT_MAX,
    "Upper bound for viscosity",

  YMIR_OPTIONS_D, "upper-mantle-scaling", '\0',
    &(rhea_viscosity_upper_mantle_scaling),
    RHEA_VISCOSITY_DEFAULT_UPPER_MANTLE_SCALING,
    "UM-Scaling factor for viscosity",
  YMIR_OPTIONS_D, "upper-mantle-activation-energy", '\0',
    &(rhea_viscosity_upper_mantle_activation_energy),
    RHEA_VISCOSITY_DEFAULT_UPPER_MANTLE_ACTIVATION_ENERGY,
    "UM-Activation energy (or exp. decay) of temperature dependent viscosity",
  YMIR_OPTIONS_D, "lower-mantle-scaling", '\0',
    &(rhea_viscosity_lower_mantle_scaling),
    RHEA_VISCOSITY_DEFAULT_LOWER_MANTLE_SCALING,
    "LM-Scaling factor for viscosity",
  YMIR_OPTIONS_D, "lower-mantle-activation-energy", '\0',
    &(rhea_viscosity_lower_mantle_activation_energy),
    RHEA_VISCOSITY_DEFAULT_LOWER_MANTLE_ACTIVATION_ENERGY,
    "LM-Activation energy (or exp. decay) of temperature dependent viscosity",

  YMIR_OPTIONS_D, "stress-exponent", '\0',
    &(rhea_viscosity_stress_exponent), RHEA_VISCOSITY_DEFAULT_STRESS_EXPONENT,
    "Stress exponent that governs strain rate weakening (aka. 'n')",

  YMIR_OPTIONS_D, "yield-stress", '\0',
    &(rhea_viscosity_yield_stress), RHEA_VISCOSITY_DEFAULT_YIELD_STRESS,
    "Value of viscous stress above which plastic yielding occurs",
  YMIR_OPTIONS_D, "yielding-regularization", '\0',
    &(rhea_viscosity_yielding_regularization),
    RHEA_VISCOSITY_DEFAULT_YIELDING_REGULARIZATION,
    "Regularization for yielding in [0,1] (0: full yielding .. 1: no yielding)",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);
}

void
rhea_viscosity_process_options (rhea_viscosity_options_t *opt,
                                rhea_domain_options_t *domain_options)
{
  const char         *this_fn_name = "rhea_viscosity_process_options";

  /* set viscosity type */
  opt->type = (rhea_viscosity_t) rhea_viscosity_type;
  opt->type_init_nonlinear =
    (rhea_viscosity_init_nonlinear_t) rhea_viscosity_type_init_nonlinear;

  /* set viscosity model */
  if (strcmp (rhea_viscosity_model_name, "UWYL") == 0) {
    opt->model = RHEA_VISCOSITY_MODEL_UWYL;
  }
  else if (strcmp (rhea_viscosity_model_name, "UWYL_LADD_UCUT") == 0) {
    opt->model = RHEA_VISCOSITY_MODEL_UWYL_LADD_UCUT;
  }
  else if (strcmp (rhea_viscosity_model_name, "UWYL_LADD_USHIFT") == 0) {
    opt->model = RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT;
  }
  else { /* unknown model name */
    RHEA_ABORT ("Unknown viscosity model name");
  }

  /* set viscosity bounds */
  opt->min = rhea_viscosity_min;
  opt->max = rhea_viscosity_max;
  RHEA_CHECK_ABORT (opt->min <= 0.0 || opt->max <= 0.0 || opt->min < opt->max,
                    "Invalid viscosity lower/upper bounds");

  /* store linear viscosity options */
  opt->upper_mantle_scaling = rhea_viscosity_upper_mantle_scaling;
  opt->upper_mantle_activation_energy =
    rhea_viscosity_upper_mantle_activation_energy;
  if (0.0 < rhea_viscosity_lower_mantle_scaling) {
    opt->lower_mantle_scaling = rhea_viscosity_lower_mantle_scaling;
  }
  else {
    opt->lower_mantle_scaling = rhea_viscosity_upper_mantle_scaling;
  }
  if (0.0 < rhea_viscosity_lower_mantle_activation_energy) {
    opt->lower_mantle_activation_energy =
      rhea_viscosity_lower_mantle_activation_energy;
  }
  else {
    opt->lower_mantle_activation_energy =
      rhea_viscosity_upper_mantle_activation_energy;
  }

  /* store nonlinear viscosity options */
  opt->stress_exponent = rhea_viscosity_stress_exponent;
  opt->yield_stress = rhea_viscosity_yield_stress;
  opt->yielding_regularization = rhea_viscosity_yielding_regularization;

  /* set additive component of nonlinear viscous stress coefficient */
  if (opt->type == RHEA_VISCOSITY_NONLINEAR) {
    switch (opt->model) {
    case RHEA_VISCOSITY_MODEL_UWYL:
      break;
    case RHEA_VISCOSITY_MODEL_UWYL_LADD_UCUT:
    case RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT:
      if (0 < opt->min) {
        ymir_nlstress_op_coeff_tensor_add = -2.0 * opt->min;
        RHEA_GLOBAL_INFOF (
            "%s: Overriding option ymir_nlstress_op_coeff_tensor_add = %g\n",
            this_fn_name, ymir_nlstress_op_coeff_tensor_add);
      }
      break;
    default: /* unknown viscosity model */
      RHEA_ABORT_NOT_REACHED ();
    }
  }

  /* store domain options */
  opt->domain_options = domain_options;
}

/******************************************************************************
 * Viscosity Vector
 *****************************************************************************/

ymir_vec_t *
rhea_viscosity_new (ymir_mesh_t *ymir_mesh)
{
  return ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
}

void
rhea_viscosity_destroy (ymir_vec_t *viscosity)
{
  ymir_vec_destroy (viscosity);
}

int
rhea_viscosity_check_vec_type (ymir_vec_t *vec)
{
  return (
      ymir_vec_is_dvec (vec) &&
      vec->ndfields == 1 &&
      vec->node_type == YMIR_GAUSS_NODE
  );
}

int
rhea_viscosity_is_valid (ymir_vec_t *vec)
{
  return (
      sc_dmatrix_is_valid (vec->dataown) &&
      0.0 < ymir_dvec_min_global (vec)
  );
}

/******************************************************************************
 * Linear Viscosity
 *****************************************************************************/

/**
 * Calculates linear (i.e., temperature dependent) viscosity term from
 * Arrhenius relationship.
 *
 *   visc (T) = exp (E * (0.5 - T))
 *
 *   E --- activation energy > 0
 *   T --- temperature in [0,1]
 */
static double
rhea_viscosity_linear_arrhenius (const double activation_energy,
                                 const double temp)
{
  RHEA_ASSERT (0.0 <= activation_energy);
  RHEA_ASSERT (0.0 <= temp && temp <= 1.0);

  return exp (activation_energy * (0.5 - temp));
}

/**
 * Computes linear viscosity at one node.  It incorporates temperature
 * dependence and weak zones.
 */
static double
rhea_viscosity_linear_node (const double temp,
                            const double weak,
                            const double scaling,
                            const double activation_energy,
                            rhea_viscosity_options_t *opt,
                            const int restrict_to_bounds)
{
  const rhea_viscosity_model_t  model = opt->model;
  const double        visc_min = opt->min;
  const double        visc_max = opt->max;

  double              viscosity;

  /* compute viscosity from Arrhenius relationship */
  viscosity = scaling * rhea_viscosity_linear_arrhenius (activation_energy,
                                                         temp);

  /* apply weak zone and restrict to bounds */
  switch (model) {
  case RHEA_VISCOSITY_MODEL_UWYL:
    /* restrict viscosity to upper bound */
    if (restrict_to_bounds && 0.0 < visc_max && visc_max < viscosity) {
      viscosity = visc_max;
    }

    /* multiply by weak zone */
    viscosity *= weak;

    /* restrict viscosity to lower bound */
    if (restrict_to_bounds && 0.0 < visc_min && viscosity < visc_min) {
      viscosity = visc_min;
    }
    break;

  case RHEA_VISCOSITY_MODEL_UWYL_LADD_UCUT:
  case RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT:
    /* restrict viscosity to upper bound */
    if (restrict_to_bounds && 0.0 < visc_max && visc_max < viscosity) {
      viscosity = visc_max;
    }

    /* multiply by weak zone */
    viscosity *= weak;

    /* restrict viscosity to lower bound */
    if (restrict_to_bounds && 0.0 < visc_min) {
      viscosity += visc_min;
    }
    break;

  default: /* unknown viscosity model */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* return viscosity */
  return viscosity;
}

/**
 * Computes linear viscosity in an element.
 */
static void
rhea_viscosity_linear_elem (double *_sc_restrict visc_elem,
                            const double *_sc_restrict temp_elem,
                            const double *_sc_restrict weak_elem,
                            const double *_sc_restrict x,
                            const double *_sc_restrict y,
                            const double *_sc_restrict z,
                            const int n_nodes,
                            const int *_sc_restrict Vmask,
                            rhea_viscosity_options_t *opt,
                            const int restrict_to_bounds)
{
  const int           in_temp = (temp_elem != NULL ? 1 : 0);
  const int           in_weak = (weak_elem != NULL ? 1 : 0);
  const double        temp_default = RHEA_TEMPERATURE_DEFAULT_VALUE;
  const double        weak_default = RHEA_WEAKZONE_DEFAULT_VALUE;
  rhea_domain_options_t  *domain_options = opt->domain_options;
  const double        interface_radius =
                        domain_options->lm_um_interface_radius;
  const double        interface_smoothing_width =
                        domain_options->lm_um_interface_smoothing_width;
  const double        um_scaling = opt->upper_mantle_scaling;
  const double        um_activ_energy = opt->upper_mantle_activation_energy;
  const double        lm_scaling = opt->lower_mantle_scaling;
  const double        lm_activ_energy = opt->lower_mantle_activation_energy;

  double              scaling, activ_energy;
  int                 nodeid;

  /* check input */
  RHEA_ASSERT (visc_elem != NULL);
  RHEA_ASSERT (x != NULL && y != NULL && z != NULL);
  RHEA_ASSERT (Vmask != NULL);

  /* set parameters depending on location in lower or upper mantle */
  if (rhea_domain_elem_is_in_upper_mantle (x, y, z, Vmask, domain_options)) {
    scaling = um_scaling;
    activ_energy = um_activ_energy;
  }
  else {
    scaling = lm_scaling;
    activ_energy = lm_activ_energy;
  }

  /* compute viscosity at each node of this element */
  for (nodeid = 0; nodeid < n_nodes; nodeid++) {
    const double        r = rhea_domain_compute_radius (x[nodeid], y[nodeid],
                                                        z[nodeid],
                                                        domain_options);
    const int           interface_smoothing =
      ( 0.0 < interface_smoothing_width &&
        fabs (r - interface_radius) < interface_smoothing_width );
    const double        temp = (in_temp ? temp_elem[nodeid] : temp_default);
    const double        weak = (in_weak ? weak_elem[nodeid] : weak_default);

    /* check temperature for valid range [0,1] */
    RHEA_ASSERT (isfinite (temp));
    RHEA_ASSERT (0.0 <= temp && temp <= 1.0);
    /* check weak zone for valid range (0,1] */
    RHEA_ASSERT (isfinite (weak));
    RHEA_ASSERT (0.0 < weak && weak <= 1.0);

    /* compute viscosity */
    if (!interface_smoothing) {
      /* compute viscosity "sufficiently far" from LM-UM interface or with a
       * discontinuous LM-UM interface */
      visc_elem[nodeid] = rhea_viscosity_linear_node (
          temp, weak, scaling, activ_energy, opt, restrict_to_bounds);
    }
    else { /* if close to LM-UM interface and must apply smoothing */
      const double        c =
        (interface_smoothing_width + (r - interface_radius)) /
        (2.0 * interface_smoothing_width);
      double              visc_lm, visc_um;

      RHEA_ASSERT (0.0 <= c && c <= 1.0);

      /* compute viscosity with smooth transition around LM-UM interface (via
       * convex combination) */
      visc_lm = rhea_viscosity_linear_node (
          temp, weak, lm_scaling, lm_activ_energy, opt, restrict_to_bounds);
      visc_um = rhea_viscosity_linear_node (
          temp, weak, um_scaling, um_activ_energy, opt, restrict_to_bounds);
      visc_elem[nodeid] = (1.0 - c) * visc_lm + c * visc_um;
    }

    /* check viscosity for `nan`, `inf`, and positivity */
    RHEA_ASSERT (isfinite (visc_elem[nodeid]));
    RHEA_ASSERT (0.0 < visc_elem[nodeid]);
  }
}

/**
 * Computes the linear viscosity.
 */
static void
rhea_viscosity_linear_vec (ymir_vec_t *visc_vec,
                           ymir_vec_t *temp_vec,
                           ymir_vec_t *weak_vec,
                           rhea_viscosity_options_t *opt)
{
  const int           restrict_to_bounds = 1;
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (visc_vec);
  const ymir_locidx_t  n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);
  const int          *Vmask = ymir_mesh_get_vertex_indices (mesh);
  const int           in_temp = (temp_vec != NULL ? 1 : 0);
  const int           in_weak = (weak_vec != NULL ? 1 : 0);

  sc_dmatrix_t       *temp_el_mat, *weak_el_mat, *visc_el_mat;
  double             *temp_el_data, *weak_el_data, *visc_el_data;
  double             *x, *y, *z, *tmp_el;
  ymir_locidx_t       elid;

  /* check input */
  RHEA_ASSERT (visc_vec != NULL);

  /* create work variables */
  temp_el_mat = (in_temp ? sc_dmatrix_new (n_nodes_per_el, 1) : NULL);
  weak_el_mat = (in_weak ? sc_dmatrix_new (n_nodes_per_el, 1) : NULL);
  visc_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  x = RHEA_ALLOC (double, n_nodes_per_el);
  y = RHEA_ALLOC (double, n_nodes_per_el);
  z = RHEA_ALLOC (double, n_nodes_per_el);
  tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

  temp_el_data = NULL;
  weak_el_data = NULL;
  visc_el_data = visc_el_mat->e[0];

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    /* get coordinates at Gauss nodes */
    ymir_mesh_get_elem_coord_gauss (x, y, z, elid, mesh, tmp_el);

    /* get temperature and weak zone at Gauss nodes */
    if (in_temp) {
      temp_el_data = rhea_temperature_get_elem_gauss (temp_el_mat, temp_vec,
                                                      elid);
    }
    if (in_weak) {
      weak_el_data = rhea_weakzone_get_elem_gauss (weak_el_mat, weak_vec, elid);
    }

    /* compute linear viscosity */
    rhea_viscosity_linear_elem (visc_el_data, temp_el_data, weak_el_data,
                                x, y, z, n_nodes_per_el, Vmask, opt,
                                restrict_to_bounds);

    /* set viscosity */
    rhea_viscosity_set_elem_gauss (visc_vec, visc_el_mat, elid);
  }

  /* destroy */
  if (in_temp) {
    sc_dmatrix_destroy (temp_el_mat);
  }
  if (in_weak) {
    sc_dmatrix_destroy (weak_el_mat);
  }
  sc_dmatrix_destroy (visc_el_mat);
  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);
}

/******************************************************************************
 * Nonlinear Viscosity
 *****************************************************************************/

/**
 * Calculates the (nonlinear) strain rate weakening viscosity.
 */
static void
rhea_viscosity_nonlinear_strain_rate_weakening (double *viscosity,
                                                double *rank1_scal,
                                                const double visc_in,
                                                const double strain_rate_2inv,
                                                const double stress_exp)
{
  /* check input */
  RHEA_ASSERT (0.0 < visc_in);
  RHEA_ASSERT (0.0 <= strain_rate_2inv);
  RHEA_ASSERT (1.0 <= stress_exp);

  /* compute viscosity
   *
   *   visc = visc_in * strain_rate_2inv^(1/n) / strain_rate_2inv
   */
  *viscosity = visc_in * pow (strain_rate_2inv, 1.0/stress_exp) /
               strain_rate_2inv;

  /* compute scaling of the rank-1 4th-order tensor
   *
   *   rank1_scal = (1 - n) / n
   */
  *rank1_scal = (1.0 - stress_exp) / stress_exp;
}

/**
 * Calculates the (nonlinear) strain rate weakening viscosity with a shift of
 * the 2nd invariant of the strain rate.
 */
static void
rhea_viscosity_nonlinear_strain_rate_weakening_shift (
                                                double *viscosity,
                                                double *rank1_scal,
                                                const double visc_in,
                                                const double strain_rate_2inv,
                                                const double shift,
                                                const double stress_exp)
{
  /* check input */
  RHEA_ASSERT (0.0 < visc_in);
  RHEA_ASSERT (0.0 <= strain_rate_2inv);
  RHEA_ASSERT (0.0 <= shift);
  RHEA_ASSERT (1.0 <= stress_exp);

  /* compute viscosity
   *
   *   visc = visc_in * (strain_rate_2inv - shift)^(1/n) / strain_rate_2inv
   */
  *viscosity = visc_in * pow (strain_rate_2inv - shift, 1.0/stress_exp) /
               strain_rate_2inv;

  /* compute scaling of the rank-1 4th-order tensor
   *
   *   rank1_scal = (strain_rate_2inv/n - strain_rate_2inv + shift)
   *                / (strain_rate_2inv - shift)
   *              = (1 - n) / n + shift / (n * (strain_rate_2inv - shift))
   */
  *rank1_scal = (1.0 - stress_exp) / stress_exp +
                shift / (stress_exp * (strain_rate_2inv - shift));
  *rank1_scal = SC_MIN ( SC_MAX (-1.0, *rank1_scal), 0.0);
}

/**
 * Applies plastic yielding to viscosity.
 */
static void
rhea_viscosity_nonlinear_yielding (double *viscosity,
                                   double *rank1_scal,
                                   double *yielding_active,
                                   const double strain_rate_2inv,
                                   const double yield_stress,
                                   const double yield_reg)
{
  const double        visc_stress = 2.0 * (*viscosity) * strain_rate_2inv;

  /* check input */
  RHEA_ASSERT (0.0 < *viscosity);
  RHEA_ASSERT (-1.0 < *rank1_scal && *rank1_scal <= 0.0);
  RHEA_ASSERT (0.0 <= strain_rate_2inv);

  if (0.0 < yield_stress && yield_stress < visc_stress) {
    const double        vy = (yield_stress / 2.0) / strain_rate_2inv;
    const double        r1y = -1.0;

    if (0.0 < yield_reg) { /* if regularize yielding via convex combination */
      const double        v = *viscosity;
      const double        r1 = *rank1_scal;

      RHEA_ASSERT (yield_reg <= 1.0);

      *viscosity = yield_reg * v + (1.0 - yield_reg) * vy;
      *rank1_scal = (yield_reg * v * r1 + (1.0 - yield_reg) * vy * r1y) /
                    *viscosity;
    }
    else { /* if standard yielding without regularization is applied */
      *viscosity = vy;
      *rank1_scal = r1y;
    }

    /* set marker that yielding is active */
    *yielding_active = 1.0;
  }
  else {
    /* set marker that yielding is not active */
    *yielding_active = 0.0;
  }
}

/**
 * Computes nonlinear viscosity at one node.
 */
static void
rhea_viscosity_nonlinear_node (double *viscosity, double *rank1_scal,
                               double *bounds_active, double *yielding_active,
                               const double temp, const double weak,
                               const double strain_rate_2inv,
                               const double scaling, const double activ_energy,
                               rhea_viscosity_options_t *opt)
{
  const rhea_viscosity_model_t  model = opt->model;
  const double        visc_min = SC_MAX (0.0, opt->min);
  const double        visc_max = SC_MAX (0.0, opt->max);
  const double        stress_exp = opt->stress_exponent;
  const double        yield_stress = opt->yield_stress;
  const double        yield_reg =
    SC_MIN (SC_MAX (0.0, opt->yielding_regularization), 1.0);

  double              visc_lin =
    scaling * rhea_viscosity_linear_arrhenius (activ_energy, temp);

  /* initialize marker that viscosity bounds are reached */
  *bounds_active = 0.0;

  switch (model) {
  case RHEA_VISCOSITY_MODEL_UWYL:
    /* Viscosity and its derivative w.r.t. the 2nd invariant of the strain rate
     * are computed, using the model
     *
     *   (1) upper bound, (2) weak zone, (3) yielding, (4) lower bound,
     *
     * as follows: TODO update comment
     *
     *   visc (e) = max(w*visc_max, visc_min), if e < min(e_min, e_yield)
     *   visc (e) = w*a * e^((1 - n)/n)    , if e_min <= e < min(e_yield, e_max)
     *   visc (e) = yield_stress / (2 * e)
     *              + visc_yield
     *              * (1 - e_yield / e)    , if e_yield <= e < e_max
     *   visc (e) = visc_min               , if e_max <= e
     *
     *   dvisc/dIIe (e) = 0                , if e < min(e_min, e_yield)
     *   dvisc/dIIe (e) = w*a * (1-n)/(2*n)
     *                    * e^((1 - 3*n)/n), if e_min <= e < min(e_yield, e_max)
     *   dvisc/dIIe (e) = (visc_yield*e_yield/2
     *                     - yield_stress/4)
     *                    / e^3            , if e_yield <= e < e_max
     *   dvisc/dIIe (e) = 0                , if e_max <= e
     *
     * where
     *
     *   e --- strain rate,
     *   a --- prefactor,
     *   w --- weak zone factor,
     *   n --- stress exponent.
     */
    {
      /* compute strain rate weakening viscosity */
      rhea_viscosity_nonlinear_strain_rate_weakening (
          viscosity, rank1_scal, visc_lin, strain_rate_2inv, stress_exp);

      /* (1) apply upper bound to viscosity */
      RHEA_ASSERT (0.0 < visc_max);
      if (visc_max < *viscosity) {
        *viscosity = visc_max;
        *rank1_scal = 0.0;
        *bounds_active = 1.0;
      }

      /* (2) multiply in weak factor */
      *viscosity *= weak;

      /* (3) apply yielding */
      rhea_viscosity_nonlinear_yielding (
          viscosity, rank1_scal, yielding_active, strain_rate_2inv,
          yield_stress, yield_reg);

      /* change upper bound marker if weak zone is present */
      if (fabs (*bounds_active - 1.0) < SC_EPS && weak < 0.5) {
        *bounds_active = 0.5;
      }

      /* (4) apply lower bound to viscosity */
      if (0.0 < visc_min && *viscosity < visc_min) {
        *viscosity = visc_min;
        *rank1_scal = 0.0;
        *bounds_active = -1.0;
      }
    }
    break;

  case RHEA_VISCOSITY_MODEL_UWYL_LADD_UCUT:
    /* Viscosity and its derivative w.r.t. the 2nd invariant of the strain rate
     * are computed, using the model
     *
     *   (1) upper bound, (2) weak zone, (3) yielding, (4) lower bound,
     *
     * as follows: TODO update comment
     *
     *   visc (e) = w * a * e^((1 - n) / n)           , if e < e_yield
     *   visc (e) = (yield_stress / 2) / e            , if e_yield <= e
     *
     *   dvisc/dIIe (e) = w * a * ((1 - n) / n)
     *                    * e^((1 - n) / n)
     *                    / (2 * e^2)                 , if e < e_yield
     *   dvisc/dIIe (e) = - (yield_stress / 2) / e
     *                    / (2 * e^2)                 , if e_yield <= e
     *
     * where
     *
     *   e --- strain rate,
     *   a --- prefactor,
     *   w --- weak zone factor,
     *   n --- stress exponent.
     */
    {
      /* compute strain rate weakening viscosity */
      rhea_viscosity_nonlinear_strain_rate_weakening (
          viscosity, rank1_scal, visc_lin, strain_rate_2inv, stress_exp);

      /* (1) apply upper bound to viscosity */
      RHEA_ASSERT (0.0 < visc_max);
      if (visc_max < *viscosity) {
        *viscosity = visc_max;
        *rank1_scal = 0.0;
        *bounds_active = 1.0;
      }

      /* (2) multiply in weak factor */
      *viscosity *= weak;

      /* change upper bound marker if weak zone is present */
      if (fabs (*bounds_active - 1.0) < SC_EPS && weak < 0.5) {
        *bounds_active = 0.5;
      }

      /* (3) apply yielding */
      rhea_viscosity_nonlinear_yielding (
          viscosity, rank1_scal, yielding_active, strain_rate_2inv,
          yield_stress, yield_reg);

      /* (4) apply lower bound to viscosity */
      if (0.0 < visc_min) {
        *viscosity += visc_min;
        RHEA_ASSERT (
            fabs (ymir_nlstress_op_coeff_tensor_add + 2.0 * visc_min) < SC_EPS);
      }
    }
    break;

  case RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT:
    /* Viscosity and its derivative w.r.t. the 2nd invariant of the strain rate
     * are computed, using the model
     *
     *   (1) upper bound + shift, (2) weak zone, (3) yielding, (4) lower bound,
     *
     * as follows:
     *
     *   visc (e) = w * a * (e - d)^(1/n) / e         , if e < e_yield
     *   visc (e) = (yield_stress / 2) / e            , if e_yield <= e
     *
     *   dvisc/dIIe (e) = w * a * (e - d)^((1 - n)/n)
     *                    * (1 / n - 1 + d / e) / e
     *                    / (2 * e^2)                 , if e < e_yield
     *   dvisc/dIIe (e) = - (yield_stress / 2) / e
     *                    / (2 * e^2)                 , if e_yield <= e
     *
     * where
     *
     *   e --- strain rate,
     *   a --- prefactor,
     *   w --- weak zone factor,
     *   n --- stress exponent,
     *   d --- shift (greater zero).
     */
    {
      double              strain_rate_min;
      double              shift;

      /* compute minimal strain rate for nonlinear viscosity
       *
       *   e_min = a / visc_max * (visc_max * n / a)^(1 / (1 - n))
       */
      RHEA_ASSERT (0.0 < visc_max);
      strain_rate_min = visc_lin / visc_max
                        * pow ( visc_max * stress_exp / visc_lin,
                                1.0 / (1.0 - stress_exp) );

      /* compute shift of strain rate
       *
       *   d = e_min - (visc_max * n / a)^(n / (1 - n))
       */
      shift = strain_rate_min -
              pow ( visc_max * stress_exp / visc_lin,
                    stress_exp / (1.0 - stress_exp) );
      shift = SC_MAX (shift, 0.0);

      /* compute strain rate weakening viscosity */
      rhea_viscosity_nonlinear_strain_rate_weakening_shift (
          viscosity, rank1_scal, visc_lin, strain_rate_2inv, shift, stress_exp);

      /* (1) apply upper bound */
      if (strain_rate_2inv < strain_rate_min) {
        *viscosity = visc_max;
        *rank1_scal = 0.0;
        *bounds_active = 1.0;
      }

      /* (2) multiply in weak factor */
      *viscosity *= weak;

      /* change upper bound marker if weak zone is present */
      if (fabs (*bounds_active - 1.0) < SC_EPS && weak < 0.5) {
        *bounds_active = 0.5;
      }

      /* (3) apply yielding */
      rhea_viscosity_nonlinear_yielding (
          viscosity, rank1_scal, yielding_active, strain_rate_2inv,
          yield_stress, yield_reg);

      /* (4) apply lower bound to viscosity */
      if (0.0 < visc_min) {
        *viscosity += visc_min;
        RHEA_ASSERT (
            fabs (ymir_nlstress_op_coeff_tensor_add + 2.0 * visc_min) < SC_EPS);
      }
    }
    break;

  default: /* unknown viscosity model */
    RHEA_ABORT_NOT_REACHED ();
  }
}

/**
 * Computes nonlinear viscosity in an element.
 */
static void
rhea_viscosity_nonlinear_elem (double *_sc_restrict visc_elem,
                               double *_sc_restrict rank1_scal_elem,
                               double *_sc_restrict bounds_elem,
                               double *_sc_restrict yielding_elem,
                               const double *_sc_restrict temp_elem,
                               const double *_sc_restrict weak_elem,
                               const double *_sc_restrict strain_rate_2inv_elem,
                               const double *_sc_restrict x,
                               const double *_sc_restrict y,
                               const double *_sc_restrict z,
                               const int n_nodes,
                               const int *_sc_restrict Vmask,
                               rhea_viscosity_options_t *opt)
{
  const int           in_temp = (temp_elem != NULL ? 1 : 0);
  const int           in_weak = (weak_elem != NULL ? 1 : 0);
  const double        temp_default = RHEA_TEMPERATURE_DEFAULT_VALUE;
  const double        weak_default = RHEA_WEAKZONE_DEFAULT_VALUE;
  rhea_domain_options_t  *domain_options = opt->domain_options;
  const double        um_scaling = opt->upper_mantle_scaling;
  const double        um_activ_energy = opt->upper_mantle_activation_energy;
  const double        lm_scaling = opt->lower_mantle_scaling;
  const double        lm_activ_energy = opt->lower_mantle_activation_energy;
  int                 nodeid;

  /* check input */
  RHEA_ASSERT (visc_elem != NULL);
  RHEA_ASSERT (x != NULL && y != NULL && z != NULL);
  RHEA_ASSERT (Vmask != NULL);
  //TODO not implemented:
  RHEA_ASSERT (domain_options->lm_um_interface_smoothing_width <= 0.0);

  /* compute viscosity depending on location in lower or upper mantle */
  if (rhea_domain_elem_is_in_upper_mantle (x, y, z, Vmask, domain_options)) {
    double              r1, b, y;

    for (nodeid = 0; nodeid < n_nodes; nodeid++) {
      const double        temp = (in_temp ? temp_elem[nodeid] : temp_default);
      const double        weak = (in_weak ? weak_elem[nodeid] : weak_default);
      const double        sr2 = ( strain_rate_2inv_elem != NULL ?
                                  strain_rate_2inv_elem[nodeid] : 1.0 );

      /* check temperature for valid range [0,1] */
      RHEA_ASSERT (isfinite (temp));
      RHEA_ASSERT (0.0 <= temp && temp <= 1.0);
      /* check weak zone for valid range (0,1] */
      RHEA_ASSERT (isfinite (weak));
      RHEA_ASSERT (0.0 < weak && weak <= 1.0);
      /* check 2nd invariant of the strain rate for non-negativity */
      RHEA_ASSERT (isfinite (sr2));
      RHEA_ASSERT (0.0 <= sr2);

      /* compute nonlinear viscosity in upper mantle */
      rhea_viscosity_nonlinear_node (
          &visc_elem[nodeid], &r1, &b, &y, temp, weak, sr2,
          um_scaling, um_activ_energy, opt);

      /* store values in output arrays if they exist */
      if (rank1_scal_elem != NULL) {
        rank1_scal_elem[nodeid] = r1;
      }
      if (bounds_elem != NULL) {
        bounds_elem[nodeid] = b;
      }
      if (yielding_elem != NULL) {
        yielding_elem[nodeid] = y;
      }
    }
  }
  else { /* if element is located in lower mantle */
    const int           restrict_to_bounds = 1;
    const double        r1 = 0.0;
    const double        b = 0.0;
    const double        y = 0.0;

    for (nodeid = 0; nodeid < n_nodes; nodeid++) {
      const double        temp = (in_temp ? temp_elem[nodeid] : temp_default);
      const double        weak = (in_weak ? weak_elem[nodeid] : weak_default);

      /* check temperature for valid range [0,1] */
      RHEA_ASSERT (isfinite (temp));
      RHEA_ASSERT (0.0 <= temp && temp <= 1.0);
      /* check weak zone for valid range (0,1] */
      RHEA_ASSERT (isfinite (weak));
      RHEA_ASSERT (0.0 < weak && weak <= 1.0);

      /* compute linear viscosity in lower mantle */
      visc_elem[nodeid] = rhea_viscosity_linear_node (
          temp, weak, lm_scaling, lm_activ_energy, opt, restrict_to_bounds);

      /* store (default) values in output arrays if they exist */
      if (rank1_scal_elem != NULL) {
        rank1_scal_elem[nodeid] = r1;
      }
      if (bounds_elem != NULL) {
        bounds_elem[nodeid] = b;
      }
      if (yielding_elem != NULL) {
        yielding_elem[nodeid] = y;
      }
    }
  }

  /* check results */
#ifdef RHEA_ENABLE_DEBUG
  for (nodeid = 0; nodeid < n_nodes; nodeid++) {
    /* check viscosity for `nan`, `inf`, and positivity */
    RHEA_ASSERT (isfinite (visc_elem[nodeid]));
    RHEA_ASSERT (0.0 < visc_elem[nodeid]);

    /* check rank-1 tensor scaling for `nan`, `inf`, and valid range [-1,0] */
    if (rank1_scal_elem != NULL) {
      RHEA_ASSERT (isfinite (rank1_scal_elem[nodeid]));
      RHEA_ASSERT (-1.0 <= rank1_scal_elem[nodeid]);
      RHEA_ASSERT (rank1_scal_elem[nodeid] <= 0.0);
    }

    /* check bounds marker for `nan`, `inf`, and valid range [-1,1] */
    if (bounds_elem != NULL) {
      RHEA_ASSERT (isfinite (bounds_elem[nodeid]));
      RHEA_ASSERT (-1.0 <= bounds_elem[nodeid]);
      RHEA_ASSERT (bounds_elem[nodeid] <= 1.0);
    }

    /* check yielding marker for `nan`, `inf` and valid range [0,1] */
    if (yielding_elem != NULL) {
      RHEA_ASSERT (isfinite (yielding_elem[nodeid]));
      RHEA_ASSERT (0.0 <= yielding_elem[nodeid]);
      RHEA_ASSERT (yielding_elem[nodeid] <= 1.0);
    }
  }
#endif
}

/**
 * Computes the nonlinear viscosity.
 */
static void
rhea_viscosity_nonlinear_vec (ymir_vec_t *visc_vec,
                              ymir_vec_t *rank1_scal_vec,
                              ymir_vec_t *bounds_vec,
                              ymir_vec_t *yielding_vec,
                              ymir_vec_t *temp_vec,
                              ymir_vec_t *weak_vec,
                              ymir_vec_t *vel_vec,
                              rhea_viscosity_options_t *opt)
{
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (visc_vec);
  const ymir_locidx_t  n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);
  const int          *Vmask = ymir_mesh_get_vertex_indices (mesh);
  const int           in_temp = (temp_vec != NULL ? 1 : 0);
  const int           in_weak = (weak_vec != NULL ? 1 : 0);
  const int           out_rank1 = (rank1_scal_vec != NULL ? 1 : 0);
  const int           out_bounds = (bounds_vec != NULL ? 1 : 0);
  const int           out_yielding = (yielding_vec != NULL ? 1 : 0);

  sc_dmatrix_t       *temp_el_mat, *weak_el_mat, *vel_el_mat,
                     *strain_rate_2inv_el_mat;
  double             *temp_el_data, *weak_el_data, *strain_rate_2inv_el_data;
  sc_dmatrix_t       *visc_el_mat, *rank1_scal_el_mat,
                     *bounds_el_mat, *yielding_el_mat;
  double             *visc_el_data, *rank1_scal_el_data,
                     *bounds_el_data, *yielding_el_data;
  sc_dmatrix_t       *tmp_grad_vel, *tmp_dvel, *tmp_vel;
  double             *x, *y, *z, *tmp_el;
  ymir_locidx_t       elid;

  /* check input */
  RHEA_ASSERT (vel_vec != NULL);
  RHEA_ASSERT (visc_vec != NULL);
  RHEA_ASSERT (0.0 < opt->max);

  /* create work variables */
  /* *INDENT-OFF* */
  temp_el_mat = (in_temp ? sc_dmatrix_new (n_nodes_per_el, 1) : NULL);
  weak_el_mat = (in_weak ? sc_dmatrix_new (n_nodes_per_el, 1) : NULL);
  vel_el_mat  = sc_dmatrix_new (n_nodes_per_el, 3);
  strain_rate_2inv_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  visc_el_mat        = sc_dmatrix_new (n_nodes_per_el, 1);
  rank1_scal_el_mat  = (out_rank1 ? sc_dmatrix_new (n_nodes_per_el, 1) : NULL);
  bounds_el_mat      = (out_bounds ? sc_dmatrix_new (n_nodes_per_el, 1) : NULL);
  yielding_el_mat    = (out_yielding ? sc_dmatrix_new (n_nodes_per_el, 1) :
                                       NULL);
  tmp_grad_vel = sc_dmatrix_new (n_nodes_per_el, 9);
  tmp_dvel     = sc_dmatrix_new (n_nodes_per_el, 3);
  tmp_vel      = sc_dmatrix_new (n_nodes_per_el, 3);
  x = RHEA_ALLOC (double, n_nodes_per_el);
  y = RHEA_ALLOC (double, n_nodes_per_el);
  z = RHEA_ALLOC (double, n_nodes_per_el);
  tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

  temp_el_data             = NULL;
  weak_el_data             = NULL;
  strain_rate_2inv_el_data = strain_rate_2inv_el_mat->e[0];
  visc_el_data       = visc_el_mat->e[0];
  rank1_scal_el_data = (out_rank1 ? rank1_scal_el_mat->e[0] : NULL);
  bounds_el_data     = (out_bounds ? bounds_el_mat->e[0] : NULL);
  yielding_el_data   = (out_yielding ? yielding_el_mat->e[0] : NULL);
  /* *INDENT-ON* */

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    /* get coordinates of this element at Gauss nodes */
    ymir_mesh_get_elem_coord_gauss (x, y, z, elid, mesh, tmp_el);

    /* get temperature and weak zone at Gauss nodes */
    if (in_temp) {
      temp_el_data = rhea_temperature_get_elem_gauss (temp_el_mat, temp_vec,
                                                      elid);
    }
    if (in_weak) {
      weak_el_data = rhea_weakzone_get_elem_gauss (weak_el_mat, weak_vec, elid);
    }

    /* get velocity; compute 2nd invariant of the strain rate at Gauss nodes */
    rhea_strainrate_compute_sqrt_of_2inv_elem_gauss (
        strain_rate_2inv_el_mat, vel_el_mat, vel_vec, elid,
        tmp_grad_vel, tmp_dvel, tmp_vel);

    /* compute nonlinear viscosity */
    rhea_viscosity_nonlinear_elem (
        visc_el_data, rank1_scal_el_data, bounds_el_data, yielding_el_data,
        temp_el_data, weak_el_data, strain_rate_2inv_el_data, x, y, z,
        n_nodes_per_el, Vmask, opt);

    /* set viscosity and other output vectors */
    rhea_viscosity_set_elem_gauss (visc_vec, visc_el_mat, elid);
    if (out_rank1) {
      rhea_viscosity_rank1_scal_set_elem_gauss (rank1_scal_vec,
                                                rank1_scal_el_mat, elid);
    }
    if (out_bounds) {
      rhea_viscosity_marker_set_elem_gauss (bounds_vec, bounds_el_mat, elid);
    }
    if (out_yielding) {
      rhea_viscosity_marker_set_elem_gauss (yielding_vec, yielding_el_mat,
                                            elid);
    }
  }

  /* destroy */
  if (in_temp) {
    sc_dmatrix_destroy (temp_el_mat);
  }
  if (in_weak) {
    sc_dmatrix_destroy (weak_el_mat);
  }
  sc_dmatrix_destroy (vel_el_mat);
  sc_dmatrix_destroy (strain_rate_2inv_el_mat);
  sc_dmatrix_destroy (visc_el_mat);
  if (out_rank1) {
    sc_dmatrix_destroy (rank1_scal_el_mat);
  }
  if (out_bounds) {
    sc_dmatrix_destroy (bounds_el_mat);
  }
  if (out_yielding) {
    sc_dmatrix_destroy (yielding_el_mat);
  }
  sc_dmatrix_destroy (tmp_grad_vel);
  sc_dmatrix_destroy (tmp_dvel);
  sc_dmatrix_destroy (tmp_vel);
  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);
}

/******************************************************************************
 * Viscosity Computation
 *****************************************************************************/

void
rhea_viscosity_compute (ymir_vec_t *viscosity,
                        ymir_vec_t *rank1_tensor_scal,
                        ymir_vec_t *bounds_marker,
                        ymir_vec_t *yielding_marker,
                        ymir_vec_t *temperature,
                        ymir_vec_t *weakzone,
                        ymir_vec_t *velocity,
                        rhea_viscosity_options_t *opt)
{
  int                 is_nonlinear = 0;

  switch (opt->type) {
  case RHEA_VISCOSITY_LINEAR:
    rhea_viscosity_linear_vec (viscosity, temperature, weakzone, opt);
    break;

  case RHEA_VISCOSITY_NONLINEAR:
    if ( fabs (opt->stress_exponent - 1.0) < SC_EPS &&
         opt->yield_stress <= 0.0 ) { /* if viscosity would be linear */
      rhea_viscosity_linear_vec (viscosity, temperature, weakzone, opt);
    }
    else { /* otherwise viscosity is nonlinear */
      is_nonlinear = 1;
      rhea_viscosity_nonlinear_vec (viscosity, rank1_tensor_scal,
                                    bounds_marker, yielding_marker,
                                    temperature, weakzone, velocity, opt);
    }
    break;

  default: /* unknown viscosity type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* set (default) values for output pertaining to nonlinear viscosity */
  if (!is_nonlinear) {
    if (rank1_tensor_scal != NULL) {
      ymir_dvec_set_zero (rank1_tensor_scal);
    }
    if (bounds_marker != NULL) {
      ymir_dvec_set_zero (bounds_marker);
    }
    if (yielding_marker != NULL) {
      ymir_dvec_set_zero (yielding_marker);
    }
  }
}

/**
 * Computes viscosity to initialize a nonlinear solver with zero velocity.
 */
void
rhea_viscosity_compute_init_nonlinear (ymir_vec_t *viscosity,
                                       ymir_vec_t *rank1_tensor_scal,
                                       ymir_vec_t *bounds_marker,
                                       ymir_vec_t *yielding_marker,
                                       ymir_vec_t *temperature,
                                       ymir_vec_t *weakzone,
                                       rhea_viscosity_options_t *opt)
{
  const rhea_viscosity_init_nonlinear_t  type = opt->type_init_nonlinear;

  switch (type) {
  case RHEA_VISCOSITY_INIT_NONLINEAR_DEFAULT:
    {
      ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (viscosity);
      ymir_vec_t         *velocity = rhea_velocity_new (ymir_mesh);

      /* compute nonlinear viscosity with zero velocity */
      ymir_vec_set_zero (velocity);
      rhea_viscosity_nonlinear_vec (
          viscosity, rank1_tensor_scal, bounds_marker, yielding_marker,
          temperature, weakzone, velocity, opt);

      ymir_vec_destroy (velocity);
    }
    break;

  case RHEA_VISCOSITY_INIT_NONLINEAR_LIN:
  case RHEA_VISCOSITY_INIT_NONLINEAR_LIN_RESCALE_UM:
    {
      const double        upper_mantle_scaling = opt->upper_mantle_scaling;

      /* compute linear viscosity (possibly rescale upper mantle) */
      if (RHEA_VISCOSITY_INIT_NONLINEAR_LIN_RESCALE_UM == type) {
        opt->upper_mantle_scaling = opt->lower_mantle_scaling;
      }
      rhea_viscosity_linear_vec (viscosity, temperature, weakzone, opt);
      if (RHEA_VISCOSITY_INIT_NONLINEAR_LIN_RESCALE_UM == type) {
        opt->upper_mantle_scaling = upper_mantle_scaling;
      }

      /* set all other output vectors to zero */
      if (rank1_tensor_scal != NULL) {
        ymir_dvec_set_zero (rank1_tensor_scal);
      }
      if (bounds_marker != NULL) {
        ymir_dvec_set_zero (bounds_marker);
      }
      if (yielding_marker != NULL) {
        ymir_dvec_set_zero (yielding_marker);
      }
    }
    break;

  default: /* unknown initial nonlinear viscosity type */
    RHEA_ABORT_NOT_REACHED ();
  }
}

/******************************************************************************
 * Get & Set Functions
 *****************************************************************************/

double *
rhea_viscosity_get_elem_gauss (sc_dmatrix_t *visc_el_mat, ymir_vec_t *visc_vec,
                               const ymir_locidx_t elid)
{
#ifdef RHEA_ENABLE_DEBUG
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (visc_vec);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (visc_vec));
  RHEA_ASSERT (visc_el_mat->m == n_nodes_per_el);
  RHEA_ASSERT (visc_el_mat->n == 1);
#endif

  ymir_dvec_get_elem (visc_vec, visc_el_mat, YMIR_STRIDE_NODE, elid, YMIR_READ);
  return visc_el_mat->e[0];
}

void
rhea_viscosity_set_elem_gauss (ymir_vec_t *visc_vec, sc_dmatrix_t *visc_el_mat,
                               const ymir_locidx_t elid)
{
#ifdef RHEA_ENABLE_DEBUG
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (visc_vec);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  /* check input */
  RHEA_ASSERT (rhea_viscosity_check_vec_type (visc_vec));
  RHEA_ASSERT (visc_el_mat->m == n_nodes_per_el);
  RHEA_ASSERT (visc_el_mat->n == 1);
#endif

  ymir_dvec_set_elem (visc_vec, visc_el_mat, YMIR_STRIDE_NODE, elid, YMIR_SET);
}

double *
rhea_viscosity_rank1_scal_get_elem_gauss (sc_dmatrix_t *rank1_scal_el_mat,
                                          ymir_vec_t *rank1_scal_vec,
                                          const ymir_locidx_t elid)
{
#ifdef RHEA_ENABLE_DEBUG
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (rank1_scal_vec);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  /* check input */
  YMIR_ASSERT_IS_DVEC (rank1_scal_vec);
  RHEA_ASSERT (rank1_scal_vec->node_type == YMIR_GAUSS_NODE);
  RHEA_ASSERT (rank1_scal_el_mat->m == n_nodes_per_el);
  RHEA_ASSERT (rank1_scal_el_mat->n == 1);
#endif

  ymir_dvec_get_elem (rank1_scal_vec, rank1_scal_el_mat, YMIR_STRIDE_NODE,
                      elid, YMIR_READ);
  return rank1_scal_el_mat->e[0];
}

void
rhea_viscosity_rank1_scal_set_elem_gauss (ymir_vec_t *rank1_scal_vec,
                                          sc_dmatrix_t *rank1_scal_el_mat,
                                          const ymir_locidx_t elid)
{
#ifdef RHEA_ENABLE_DEBUG
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (rank1_scal_vec);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  /* check input */
  YMIR_ASSERT_IS_DVEC (rank1_scal_vec);
  RHEA_ASSERT (rank1_scal_vec->node_type == YMIR_GAUSS_NODE);
  RHEA_ASSERT (rank1_scal_el_mat->m == n_nodes_per_el);
  RHEA_ASSERT (rank1_scal_el_mat->n == 1);
#endif

  ymir_dvec_set_elem (rank1_scal_vec, rank1_scal_el_mat, YMIR_STRIDE_NODE,
                      elid, YMIR_SET);
}

double *
rhea_viscosity_marker_get_elem_gauss (sc_dmatrix_t *marker_el_mat,
                                      ymir_vec_t *marker_vec,
                                      const ymir_locidx_t elid)
{
#ifdef RHEA_ENABLE_DEBUG
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (marker_vec);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  /* check input */
  YMIR_ASSERT_IS_DVEC (marker_vec);
  RHEA_ASSERT (marker_vec->node_type == YMIR_GAUSS_NODE);
  RHEA_ASSERT (marker_el_mat->m == n_nodes_per_el);
  RHEA_ASSERT (marker_el_mat->n == 1);
#endif

  ymir_dvec_get_elem (marker_vec, marker_el_mat, YMIR_STRIDE_NODE, elid,
                      YMIR_READ);
  return marker_el_mat->e[0];
}

void
rhea_viscosity_marker_set_elem_gauss (ymir_vec_t *marker_vec,
                                      sc_dmatrix_t *marker_el_mat,
                                      const ymir_locidx_t elid)
{
#ifdef RHEA_ENABLE_DEBUG
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (marker_vec);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  /* check input */
  YMIR_ASSERT_IS_DVEC (marker_vec);
  RHEA_ASSERT (marker_vec->node_type == YMIR_GAUSS_NODE);
  RHEA_ASSERT (marker_el_mat->m == n_nodes_per_el);
  RHEA_ASSERT (marker_el_mat->n == 1);
#endif

  ymir_dvec_set_elem (marker_vec, marker_el_mat, YMIR_STRIDE_NODE, elid,
                      YMIR_SET);
}
