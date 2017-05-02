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

/* definition of viscosity bounds and yielding markers */
#define RHEA_VISCOSITY_BOUNDS_OFF (0.0)
#define RHEA_VISCOSITY_BOUNDS_MIN (-1.0)
#define RHEA_VISCOSITY_BOUNDS_MAX (+1.0)
#define RHEA_VISCOSITY_BOUNDS_MAX_WEAK (0.5)
#define RHEA_VISCOSITY_YIELDING_OFF (0.0)
#define RHEA_VISCOSITY_YIELDING_ACTIVE (1.0)

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_VISCOSITY_DEFAULT_TYPE (RHEA_VISCOSITY_LINEAR)
#define RHEA_VISCOSITY_DEFAULT_TYPE_LINEAR (RHEA_VISCOSITY_LINEAR_ARRHENIUS)
#define RHEA_VISCOSITY_DEFAULT_TYPE_NONLINEAR (RHEA_VISCOSITY_NONLINEAR_SRW_YLD)
#define RHEA_VISCOSITY_DEFAULT_TYPE_NONLINEAR_INIT \
  (RHEA_VISCOSITY_NONLINEAR_INIT_DEFAULT)
#define RHEA_VISCOSITY_DEFAULT_MODEL_NAME "UWYL_LADD_USHIFT"
#define RHEA_VISCOSITY_DEFAULT_MIN (1.0e-2)
#define RHEA_VISCOSITY_DEFAULT_MAX (1.0e+4)
#define RHEA_VISCOSITY_DEFAULT_SHIFT (0.0)
#define RHEA_VISCOSITY_DEFAULT_UPPER_MANTLE_SCALING (4.0e+3)
#define RHEA_VISCOSITY_DEFAULT_UPPER_MANTLE_ARRHENIUS_ACTIVATION_ENERGY (17.5)
#define RHEA_VISCOSITY_DEFAULT_LOWER_MANTLE_SCALING (4.0e+5)
#define RHEA_VISCOSITY_DEFAULT_LOWER_MANTLE_ARRHENIUS_ACTIVATION_ENERGY (17.5)
#define RHEA_VISCOSITY_DEFAULT_STRESS_EXPONENT (3.0)
#define RHEA_VISCOSITY_DEFAULT_YIELD_STRENGTH (NAN)
#define RHEA_VISCOSITY_DEFAULT_YIELDING_REGULARIZATION (0.0)

/* initialize options */
int                 rhea_viscosity_type = RHEA_VISCOSITY_DEFAULT_TYPE;
int                 rhea_viscosity_type_linear =
  RHEA_VISCOSITY_DEFAULT_TYPE_LINEAR;
int                 rhea_viscosity_type_nonlinear =
  RHEA_VISCOSITY_DEFAULT_TYPE_NONLINEAR;
int                 rhea_viscosity_type_nonlinear_init =
  RHEA_VISCOSITY_DEFAULT_TYPE_NONLINEAR_INIT;
char               *rhea_viscosity_model_name =
  RHEA_VISCOSITY_DEFAULT_MODEL_NAME;
double              rhea_viscosity_min = RHEA_VISCOSITY_DEFAULT_MIN;
double              rhea_viscosity_max = RHEA_VISCOSITY_DEFAULT_MAX;
double              rhea_viscosity_shift = RHEA_VISCOSITY_DEFAULT_SHIFT;
double              rhea_viscosity_upper_mantle_scaling =
  RHEA_VISCOSITY_DEFAULT_UPPER_MANTLE_SCALING;
double              rhea_viscosity_upper_mantle_arrhenius_activation_energy =
  RHEA_VISCOSITY_DEFAULT_UPPER_MANTLE_ARRHENIUS_ACTIVATION_ENERGY;
double              rhea_viscosity_lower_mantle_scaling =
  RHEA_VISCOSITY_DEFAULT_LOWER_MANTLE_SCALING;
double              rhea_viscosity_lower_mantle_arrhenius_activation_energy =
  RHEA_VISCOSITY_DEFAULT_LOWER_MANTLE_ARRHENIUS_ACTIVATION_ENERGY;
double              rhea_viscosity_stress_exponent =
  RHEA_VISCOSITY_DEFAULT_STRESS_EXPONENT;
double              rhea_viscosity_yield_strength =
  RHEA_VISCOSITY_DEFAULT_YIELD_STRENGTH;
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
    "Viscosity type: 0: linear; 1: nonlinear",
  YMIR_OPTIONS_I, "type-linear", '\0',
    &(rhea_viscosity_type_linear), RHEA_VISCOSITY_DEFAULT_TYPE_LINEAR,
    "Linear viscosity type: 0: constant; 1: reversed temperature; 2: Arrhenius",
  YMIR_OPTIONS_I, "type-nonlinear", '\0',
    &(rhea_viscosity_type_nonlinear), RHEA_VISCOSITY_DEFAULT_TYPE_NONLINEAR,
    "Nonlinear viscosity type: 0: strain rate weakening; 1: yielding; "
    "2: strain rate weakening & yielding",
  YMIR_OPTIONS_I, "type-nonlinear-init", '\0',
    &(rhea_viscosity_type_nonlinear_init),
    RHEA_VISCOSITY_DEFAULT_TYPE_NONLINEAR_INIT,
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

  YMIR_OPTIONS_D, "shift", '\0',
    &(rhea_viscosity_shift), RHEA_VISCOSITY_DEFAULT_SHIFT,
    "Shift of viscosity to introduce convexity to energy min. problem",

  YMIR_OPTIONS_D, "upper-mantle-scaling", '\0',
    &(rhea_viscosity_upper_mantle_scaling),
    RHEA_VISCOSITY_DEFAULT_UPPER_MANTLE_SCALING,
    "UM scaling factor for viscosity",
  YMIR_OPTIONS_D, "upper-mantle-arrhenius-activation-energy", '\0',
    &(rhea_viscosity_upper_mantle_arrhenius_activation_energy),
    RHEA_VISCOSITY_DEFAULT_UPPER_MANTLE_ARRHENIUS_ACTIVATION_ENERGY,
    "UM activation energy (or exp. decay) in Arrhenius relationship",
  YMIR_OPTIONS_D, "lower-mantle-scaling", '\0',
    &(rhea_viscosity_lower_mantle_scaling),
    RHEA_VISCOSITY_DEFAULT_LOWER_MANTLE_SCALING,
    "LM scaling factor for viscosity",
  YMIR_OPTIONS_D, "lower-mantle-arrhenius-activation-energy", '\0',
    &(rhea_viscosity_lower_mantle_arrhenius_activation_energy),
    RHEA_VISCOSITY_DEFAULT_LOWER_MANTLE_ARRHENIUS_ACTIVATION_ENERGY,
    "LM activation energy (or exp. decay) in Arrhenius relationship",

  YMIR_OPTIONS_D, "stress-exponent", '\0',
    &(rhea_viscosity_stress_exponent), RHEA_VISCOSITY_DEFAULT_STRESS_EXPONENT,
    "Stress exponent that governs strain rate weakening (aka. 'n')",

  YMIR_OPTIONS_D, "yield-strength", '\0',
    &(rhea_viscosity_yield_strength), RHEA_VISCOSITY_DEFAULT_YIELD_STRENGTH,
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
  opt->type_linear = (rhea_viscosity_linear_t) rhea_viscosity_type_linear;
  opt->type_nonlinear =
    (rhea_viscosity_nonlinear_t) rhea_viscosity_type_nonlinear;
  opt->type_nonlinear_init =
    (rhea_viscosity_nonlinear_init_t) rhea_viscosity_type_nonlinear_init;

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

  /* set viscosity shift */
  opt->shift = rhea_viscosity_shift;

  /* store linear viscosity options */
  opt->upper_mantle_scaling = rhea_viscosity_upper_mantle_scaling;
  opt->upper_mantle_arrhenius_activation_energy =
    rhea_viscosity_upper_mantle_arrhenius_activation_energy;
  if (0.0 < rhea_viscosity_lower_mantle_scaling) {
    opt->lower_mantle_scaling = rhea_viscosity_lower_mantle_scaling;
  }
  else {
    opt->lower_mantle_scaling = rhea_viscosity_upper_mantle_scaling;
  }
  if (0.0 < rhea_viscosity_lower_mantle_arrhenius_activation_energy) {
    opt->lower_mantle_arrhenius_activation_energy =
      rhea_viscosity_lower_mantle_arrhenius_activation_energy;
  }
  else {
    opt->lower_mantle_arrhenius_activation_energy =
      rhea_viscosity_upper_mantle_arrhenius_activation_energy;
  }

  /* store nonlinear viscosity options */
  opt->stress_exponent = rhea_viscosity_stress_exponent;
  opt->yield_strength = rhea_viscosity_yield_strength;
  opt->yielding_regularization = rhea_viscosity_yielding_regularization;

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
 * Linear Viscosity Components and Model
 *****************************************************************************/

/**
 * Calculates a viscosity that depends linearly on the temperature:
 *
 *   visc (T) = 1 - T
 *
 *   T --- temperature in [0,1]
 */
static double
rhea_viscosity_linear_tempreverse (const double temp)
{
  RHEA_ASSERT (0.0 <= temp && temp <= 1.0);

  return (1.0 - temp);
}

/**
 * Calculates a temperature dependent viscosity term from the Arrhenius
 * relationship:
 *
 *   visc (T) = exp (E * (T_neutral - T))
 *
 *   E         --- activation energy > 0
 *   T         --- temperature in [0,1]
 *   T_neutral --- neutral/default temperature value (e.g., 0,5)
 */
static double
rhea_viscosity_linear_arrhenius (const double activation_energy,
                                 const double temp)
{
  RHEA_ASSERT (0.0 <= activation_energy);
  RHEA_ASSERT (0.0 <= temp && temp <= 1.0);

  return exp (activation_energy * (RHEA_TEMPERATURE_NEUTRAL_VALUE - temp));
}

/**
 * Calculates the linear viscosity component.
 */
static double
rhea_viscosity_linear_comp (const double temp,
                            rhea_viscosity_options_t *opt,
                            const int is_in_upper_mantle)
{
  switch (opt->type_linear) {
  case RHEA_VISCOSITY_LINEAR_CONST:
    return RHEA_VISCOSITY_NEUTRAL_VALUE;
  case RHEA_VISCOSITY_LINEAR_TEMPREVERSE:
    return rhea_viscosity_linear_tempreverse (temp);
  case RHEA_VISCOSITY_LINEAR_ARRHENIUS:
    {
      double              activation_energy;

      if (is_in_upper_mantle) {
        activation_energy = opt->upper_mantle_arrhenius_activation_energy;
      }
      else {
        activation_energy = opt->lower_mantle_arrhenius_activation_energy;
      }
      return rhea_viscosity_linear_arrhenius (activation_energy, temp);
    }
  default: /* unknown linear viscosity type */
    RHEA_ABORT_NOT_REACHED ();
  }
}

/**
 * Calculates the linear viscosity according to a specified model.
 * Incorporates temperature, weak zones, and viscosity bounds.
 */
static void
rhea_viscosity_linear_model (double *viscosity, double *bounds_active,
                             const double temp, const double weak,
                             rhea_viscosity_options_t *opt,
                             const int is_in_upper_mantle,
                             const int restrict_to_bounds)
{
  const double        visc_min = opt->min;
  const double        visc_max = opt->max;
  double              scaling;

  /* set scaling */
  if (is_in_upper_mantle) {
    scaling = opt->upper_mantle_scaling;
  }
  else {
    scaling = opt->lower_mantle_scaling;
  }

  /* initialize marker that viscosity bounds are reached */
  *bounds_active = RHEA_VISCOSITY_BOUNDS_OFF;

  /* compute linear viscosity component */
  *viscosity = rhea_viscosity_linear_comp (temp, opt, is_in_upper_mantle);

  /* compose viscosity */
  switch (opt->model) {
  case RHEA_VISCOSITY_MODEL_UWYL:
    {
      /* multiply by scaling factor */
      *viscosity *= scaling;

      /* (U) restrict viscosity to upper bound */
      if (restrict_to_bounds && 0.0 < visc_max && visc_max < *viscosity) {
        *viscosity = visc_max;
        *bounds_active = RHEA_VISCOSITY_BOUNDS_MAX;
      }

      /* (W) multiply by weak zone */
      *viscosity *= weak;
      /* change upper bound marker if weak zone is present */
      if (fabs (*bounds_active - RHEA_VISCOSITY_BOUNDS_MAX) < SC_EPS &&
          weak < 0.5) {
        *bounds_active = RHEA_VISCOSITY_BOUNDS_MAX_WEAK;
      }

      /* (L) restrict viscosity to lower bound */
      if (restrict_to_bounds && 0.0 < visc_min && *viscosity < visc_min) {
        *viscosity = visc_min;
        *bounds_active = RHEA_VISCOSITY_BOUNDS_MIN;
      }
    }
    break;

  case RHEA_VISCOSITY_MODEL_UWYL_LADD_UCUT:
  case RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT:
    {
      /* multiply by scaling factor */
      if (restrict_to_bounds && 0.0 < visc_min && visc_min < scaling) {
        *viscosity *= (scaling - visc_min);
      }
      else {
        *viscosity *= scaling;
      }

      /* (U) restrict viscosity to upper bound */
      if (restrict_to_bounds && 0.0 < visc_max && visc_max < *viscosity) {
        *viscosity = visc_max;
        *bounds_active = RHEA_VISCOSITY_BOUNDS_MAX;
      }

      /* (W) multiply by weak zone */
      *viscosity *= weak;
      /* change upper bound marker if weak zone is present */
      if (fabs (*bounds_active - RHEA_VISCOSITY_BOUNDS_MAX) < SC_EPS &&
          weak < 0.5) {
        *bounds_active = RHEA_VISCOSITY_BOUNDS_MAX_WEAK;
      }

      /* (L) restrict viscosity to lower bound */
      if (restrict_to_bounds && 0.0 < visc_min) {
        if (*viscosity < visc_min) {
          *bounds_active = RHEA_VISCOSITY_BOUNDS_MIN;
        }
        *viscosity += visc_min;
      }
    }
    break;

  default: /* unknown viscosity model */
    RHEA_ABORT_NOT_REACHED ();
  }
}

/******************************************************************************
 * Linear Viscosity Vector Computation
 *****************************************************************************/

/**
 * Computes linear viscosity in an element.
 */
static void
rhea_viscosity_linear_elem (double *_sc_restrict visc_elem,
                            double *_sc_restrict bounds_elem,
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
  const int           in_temp = (temp_elem != NULL);
  const int           in_weak = (weak_elem != NULL);
  const double        temp_default = RHEA_TEMPERATURE_NEUTRAL_VALUE;
  const double        weak_default = RHEA_WEAKZONE_NEUTRAL_VALUE;
  rhea_domain_options_t  *domain_options = opt->domain_options;
  const double        interface_radius =
                        domain_options->lm_um_interface_radius;
  const double        interface_smoothing_width =
                        domain_options->lm_um_interface_smoothing_width;

  int                 is_in_upper_mantle;
  double              bd;
  int                 nodeid;

  /* check input */
  RHEA_ASSERT (visc_elem != NULL);
  RHEA_ASSERT (x != NULL && y != NULL && z != NULL);
  RHEA_ASSERT (Vmask != NULL);

  /* set flag if location is in lower or upper mantle */
  is_in_upper_mantle = rhea_domain_elem_is_in_upper_mantle (x, y, z, Vmask,
                                                            domain_options);

  /* compute viscosity at each node of this element */
  for (nodeid = 0; nodeid < n_nodes; nodeid++) {
    const double        temp = (in_temp ? temp_elem[nodeid] : temp_default);
    const double        weak = (in_weak ? weak_elem[nodeid] : weak_default);
    const double        r = rhea_domain_compute_radius (x[nodeid], y[nodeid],
                                                        z[nodeid],
                                                        domain_options);

    /* check temperature for valid range [0,1] */
    RHEA_ASSERT (isfinite (temp));
    RHEA_ASSERT (0.0 <= temp && temp <= 1.0);
    /* check weak zone for valid range (0,1] */
    RHEA_ASSERT (isfinite (weak));
    RHEA_ASSERT (0.0 < weak && weak <= 1.0);

    /* compute viscosity */
    rhea_viscosity_linear_model (&visc_elem[nodeid], &bd, temp, weak, opt,
                                 is_in_upper_mantle, restrict_to_bounds);

    /* update viscosity by applying a smooth transition close to the LM-UM
     * interface (via a convex combination) */
    if (0.0 < interface_smoothing_width &&
        fabs (r - interface_radius) < 0.5 * interface_smoothing_width) {
      const double        v1 = visc_elem[nodeid];
      double              v2, bd2, c;

      rhea_viscosity_linear_model (&v2, &bd2, temp, weak, opt,
                                   !is_in_upper_mantle, restrict_to_bounds);

      c = (fabs (r - interface_radius) - 0.5 * interface_smoothing_width) /
          interface_smoothing_width;
      RHEA_ASSERT (0.0 <= c && c <= 1.0);

      if (is_in_upper_mantle) {
        visc_elem[nodeid] = c * v2 + (1.0 - c) * v1;
      }
      else {
        visc_elem[nodeid] = c * v1 + (1.0 - c) * v2;
      }
    }

    /* store value in output array if it exists */
    if (bounds_elem != NULL) {
      bounds_elem[nodeid] = bd;
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
                           ymir_vec_t *bounds_vec,
                           ymir_vec_t *temp_vec,
                           ymir_vec_t *weak_vec,
                           rhea_viscosity_options_t *opt)
{
  const int           restrict_to_bounds = 1;
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (visc_vec);
  const ymir_locidx_t  n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);
  const int          *Vmask = ymir_mesh_get_vertex_indices (mesh);
  const int           in_temp = (temp_vec != NULL);
  const int           in_weak = (weak_vec != NULL);
  const int           out_bounds = (bounds_vec != NULL);

  sc_dmatrix_t       *temp_el_mat, *weak_el_mat;
  double             *temp_el_data, *weak_el_data;
  sc_dmatrix_t       *visc_el_mat, *bounds_el_mat;
  double             *visc_el_data, *bounds_el_data;
  double             *x, *y, *z, *tmp_el;
  ymir_locidx_t       elid;

  /* check input */
  RHEA_ASSERT (visc_vec != NULL);

  /* create work variables */
  temp_el_mat = (in_temp ? sc_dmatrix_new (n_nodes_per_el, 1) : NULL);
  weak_el_mat = (in_weak ? sc_dmatrix_new (n_nodes_per_el, 1) : NULL);
  visc_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  bounds_el_mat = (out_bounds ? sc_dmatrix_new (n_nodes_per_el, 1) : NULL);
  x = RHEA_ALLOC (double, n_nodes_per_el);
  y = RHEA_ALLOC (double, n_nodes_per_el);
  z = RHEA_ALLOC (double, n_nodes_per_el);
  tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

  temp_el_data = NULL;
  weak_el_data = NULL;
  visc_el_data = visc_el_mat->e[0];
  bounds_el_data = (out_bounds ? bounds_el_mat->e[0] : NULL);

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
    rhea_viscosity_linear_elem (visc_el_data, bounds_el_data,
                                temp_el_data, weak_el_data,
                                x, y, z, n_nodes_per_el, Vmask, opt,
                                restrict_to_bounds);

    /* set viscosity and other output vectors */
    rhea_viscosity_set_elem_gauss (visc_vec, visc_el_mat, elid);
    if (out_bounds) {
      rhea_viscosity_marker_set_elem_gauss (bounds_vec, bounds_el_mat, elid);
    }
  }

  /* destroy */
  if (in_temp) {
    sc_dmatrix_destroy (temp_el_mat);
  }
  if (in_weak) {
    sc_dmatrix_destroy (weak_el_mat);
  }
  sc_dmatrix_destroy (visc_el_mat);
  if (out_bounds) {
    sc_dmatrix_destroy (bounds_el_mat);
  }
  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);
}

/******************************************************************************
 * Nonlinear Viscosity Components and Model
 *****************************************************************************/

/**
 * Calculates the (nonlinear) strain rate weakening viscosity.
 */
static void
rhea_viscosity_nonlinear_strain_rate_weakening (
                                            double *viscosity,
                                            double *rank1_scal,
                                            const double visc_in,
                                            const double strainrate_sqrt_2inv,
                                            const double stress_exp)
{
  /* check input */
  RHEA_ASSERT (0.0 <= visc_in);
  RHEA_ASSERT (0.0 <= strainrate_sqrt_2inv);
  RHEA_ASSERT (1.0 <= stress_exp);

  /* compute viscosity
   *
   *   visc = visc_in * sqrt(strainrate_2inv)^(1/n) / sqrt(strainrate_2inv)
   */
  *viscosity = visc_in * pow (strainrate_sqrt_2inv, 1.0/stress_exp) /
               strainrate_sqrt_2inv;

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
                                            const double strainrate_sqrt_2inv,
                                            const double shift,
                                            const double stress_exp)
{
  /* check input */
  RHEA_ASSERT (0.0 <= visc_in);
  RHEA_ASSERT (0.0 <= strainrate_sqrt_2inv);
  RHEA_ASSERT (0.0 <= shift);
  RHEA_ASSERT (1.0 <= stress_exp);

  /* compute viscosity
   *
   *   visc = visc_in * (sqrt(strainrate_2inv) - shift)^(1/n)
   *                  / sqrt(strainrate_2inv)
   */
  *viscosity = visc_in * pow (strainrate_sqrt_2inv - shift, 1.0/stress_exp) /
               strainrate_sqrt_2inv;

  /* compute scaling of the rank-1 4th-order tensor
   *
   *   rank1_scal = (sqrt(strainrate_2inv)/n - sqrt(strainrate_2inv) + shift) /
   *                (sqrt(strainrate_2inv) - shift)
   *              = (1 - n) / n + shift / (n * (sqrt(strainrate_2inv) - shift))
   */
  *rank1_scal = (1.0 - stress_exp) / stress_exp +
                shift / (stress_exp * (strainrate_sqrt_2inv - shift));
  *rank1_scal = SC_MIN ( SC_MAX (-1.0, *rank1_scal), 0.0);
}

/**
 * Applies plastic yielding to viscosity.
 */
static void
rhea_viscosity_nonlinear_yielding (double *viscosity,
                                   double *rank1_scal,
                                   double *yielding_active,
                                   const double strainrate_sqrt_2inv,
                                   const double yield_strength,
                                   const double yield_reg)
{
  double              visc_stress;

  /* check input */
  RHEA_ASSERT (0.0 <= *viscosity);
  RHEA_ASSERT (-1.0 < *rank1_scal && *rank1_scal <= 0.0);
  RHEA_ASSERT (0.0 <= strainrate_sqrt_2inv);
  RHEA_ASSERT (yield_reg <= 0.0 || yield_reg <= 1.0);

  /* exit if nothing to do */
  if ( !(isfinite (yield_strength) && 0.0 < yield_strength) ) {
    *yielding_active = RHEA_VISCOSITY_YIELDING_OFF;
    return;
  }

  /* compute sqrt of the 2nd invariant of the viscous stress tensor */
  if (0.0 < *viscosity) {
    visc_stress = 2.0 * (*viscosity) * strainrate_sqrt_2inv;
  }
  else {
    visc_stress = 0.0;
  }

  /* apply yielding if yield strenght is exceeded */
  if (yield_strength < visc_stress) {
    const double        v = *viscosity;
    const double        r1 = *rank1_scal;
    const double        vy = (yield_strength / 2.0) / strainrate_sqrt_2inv;
    const double        r1y = -1.0;

    if (0.0 < yield_reg) { /* if regularize yielding via convex combination */
      *viscosity = yield_reg * v + (1.0 - yield_reg) * vy;
      RHEA_ASSERT (0.0 < *viscosity);
      *rank1_scal = (yield_reg * v * r1 + (1.0 - yield_reg) * vy * r1y) /
                    *viscosity;
    }
    else { /* if standard yielding without regularization is applied */
      *viscosity = vy;
      *rank1_scal = r1y;
    }

    /* set yielding marker to active */
    *yielding_active = RHEA_VISCOSITY_YIELDING_ACTIVE;
  }
}

/**
 * Calculates the nonlinear viscosity according to a specified model.
 * Incorporates temperature, weak zones, viscosity bounds, strain rate
 * weakening, and yielding.
 */
static void
rhea_viscosity_nonlinear_model (double *viscosity, double *rank1_scal,
                                double *bounds_active, double *yielding_active,
                                const double temp, const double weak,
                                const double strainrate_sqrt_2inv,
                                rhea_viscosity_options_t *opt)
{
  const double        visc_min = SC_MAX (0.0, opt->min);
  const double        visc_max = SC_MAX (0.0, opt->max);
  const double        scaling = opt->upper_mantle_scaling;
  const double        stress_exp = opt->stress_exponent;
  const double        yield_strength = opt->yield_strength;
  const double        yield_reg =
    SC_MIN (SC_MAX (0.0, opt->yielding_regularization), 1.0);
  const int           has_srw = rhea_viscosity_has_strain_rate_weakening (opt);
  const int           has_yld = rhea_viscosity_has_yielding (opt);
  double              visc_lin = rhea_viscosity_linear_comp (temp, opt, 1);

  RHEA_ASSERT (0.0 < visc_max);
  RHEA_ASSERT (0.0 <= visc_lin);
  RHEA_ASSERT (0.0 < visc_lin || 0.0 < visc_min);

  /* initialize marker that viscosity bounds are reached */
  *bounds_active = RHEA_VISCOSITY_BOUNDS_OFF;

  /* compose viscosity */
  switch (opt->model) {
  case RHEA_VISCOSITY_MODEL_UWYL:
    {
      /* multiply by scaling factor */
      visc_lin *= scaling;

      /* compute strain rate weakening viscosity */
      if (has_srw) {
        rhea_viscosity_nonlinear_strain_rate_weakening (
            viscosity, rank1_scal, visc_lin, strainrate_sqrt_2inv, stress_exp);
      }
      else {
        *viscosity = visc_lin;
        *rank1_scal = 0.0;
      }

      /* (U) apply upper bound to viscosity */
      if (visc_max < *viscosity) {
        *viscosity = visc_max;
        *rank1_scal = 0.0;
        *bounds_active = RHEA_VISCOSITY_BOUNDS_MAX;
      }

      /* (W) multiply in weak factor */
      *viscosity *= weak;
      /* change upper bound marker if weak zone is present */
      if (fabs (*bounds_active - RHEA_VISCOSITY_BOUNDS_MAX) < SC_EPS &&
          weak < 0.5) {
        *bounds_active = RHEA_VISCOSITY_BOUNDS_MAX_WEAK;
      }

      /* (Y) apply yielding */
      if (has_yld) {
        rhea_viscosity_nonlinear_yielding (
            viscosity, rank1_scal, yielding_active, strainrate_sqrt_2inv,
            yield_strength, yield_reg);
      }
      else {
        *yielding_active = RHEA_VISCOSITY_YIELDING_OFF;
      }

      /* (L) apply lower bound to viscosity */
      if (0.0 < visc_min && *viscosity < visc_min) {
        *viscosity = visc_min;
        *rank1_scal = 0.0;
        *bounds_active = RHEA_VISCOSITY_BOUNDS_MIN;
      }
    }
    break;

  case RHEA_VISCOSITY_MODEL_UWYL_LADD_UCUT:
    {
      /* multiply by scaling factor */
      if (0.0 < visc_min && visc_min < scaling) {
        visc_lin *= (scaling - visc_min);
      }
      else {
        visc_lin *= scaling;
      }

      /* compute strain rate weakening viscosity */
      if (has_srw) {
        rhea_viscosity_nonlinear_strain_rate_weakening (
            viscosity, rank1_scal, visc_lin, strainrate_sqrt_2inv, stress_exp);
      }
      else {
        *viscosity = visc_lin;
        *rank1_scal = 0.0;
      }

      /* (U) apply upper bound to viscosity */
      if (visc_max < *viscosity) {
        *viscosity = visc_max;
        *rank1_scal = 0.0;
        *bounds_active = RHEA_VISCOSITY_BOUNDS_MAX;
      }

      /* (W) multiply in weak factor */
      *viscosity *= weak;
      /* change upper bound marker if weak zone is present */
      if (fabs (*bounds_active - RHEA_VISCOSITY_BOUNDS_MAX) < SC_EPS &&
          weak < 0.5) {
        *bounds_active = RHEA_VISCOSITY_BOUNDS_MAX_WEAK;
      }

      /* (Y) apply yielding */
      if (has_yld) {
        rhea_viscosity_nonlinear_yielding (
            viscosity, rank1_scal, yielding_active, strainrate_sqrt_2inv,
            yield_strength, yield_reg);
      }
      else {
        *yielding_active = RHEA_VISCOSITY_YIELDING_OFF;
      }

      /* (L) apply lower bound to viscosity */
      if (0.0 < visc_min) {
        if (*viscosity < visc_min) {
          *bounds_active = RHEA_VISCOSITY_BOUNDS_MIN;
        }
        *viscosity += visc_min;
      }
    }
    break;

  case RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT:
    {
      double              strainrate_min;
      double              shift;

      /* multiply by scaling factor */
      if (0.0 < visc_min && visc_min < scaling) {
        visc_lin *= (scaling - visc_min);
      }
      else {
        visc_lin *= scaling;
      }

      /* compute minimal strain rate for nonlinear viscosity
       *
       *   e_min = a / visc_max * (visc_max * n / a)^(1 / (1 - n))
       */
      if (0.0 < visc_lin) {
        strainrate_min = visc_lin / visc_max *
                         pow ( visc_max * stress_exp / visc_lin,
                               1.0 / (1.0 - stress_exp) );
      }
      else {
        strainrate_min = 0.0;
      }

      /* compute shift of strain rate
       *
       *   d = e_min - (visc_max * n / a)^(n / (1 - n))
       */
      if (0.0 < visc_lin) {
        shift = strainrate_min -
                pow ( visc_max * stress_exp / visc_lin,
                      stress_exp / (1.0 - stress_exp) );
        shift = SC_MAX (shift, 0.0);
      }
      else {
        shift = 0.0;
      }

      /* compute strain rate weakening viscosity */
      if (has_srw) {
        rhea_viscosity_nonlinear_strain_rate_weakening_shift (
            viscosity, rank1_scal, visc_lin, strainrate_sqrt_2inv, shift,
            stress_exp);
      }
      else {
        *viscosity = visc_lin;
        *rank1_scal = 0.0;
      }

      /* (U) apply upper bound */
      if (strainrate_sqrt_2inv < strainrate_min) {
        *viscosity = visc_max;
        *rank1_scal = 0.0;
        *bounds_active = RHEA_VISCOSITY_BOUNDS_MAX;
      }

      /* (W) multiply in weak factor */
      *viscosity *= weak;
      /* change upper bound marker if weak zone is present */
      if (fabs (*bounds_active - RHEA_VISCOSITY_BOUNDS_MAX) < SC_EPS &&
          weak < 0.5) {
        *bounds_active = RHEA_VISCOSITY_BOUNDS_MAX_WEAK;
      }

      /* (Y) apply yielding */
      if (has_yld) {
        rhea_viscosity_nonlinear_yielding (
            viscosity, rank1_scal, yielding_active, strainrate_sqrt_2inv,
            yield_strength, yield_reg);
      }
      else {
        *yielding_active = RHEA_VISCOSITY_YIELDING_OFF;
      }

      /* (L) apply lower bound to viscosity */
      if (0.0 < visc_min) {
        if (*viscosity < visc_min) {
          *bounds_active = RHEA_VISCOSITY_BOUNDS_MIN;
        }
        *viscosity += visc_min;
      }
    }
    break;

  default: /* unknown viscosity model */
    RHEA_ABORT_NOT_REACHED ();
  }
}

/******************************************************************************
 * Nonlinear Viscosity Vector Computation
 *****************************************************************************/

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
                               const double *_sc_restrict strt_sqrt_2inv_elem,
                               const double *_sc_restrict x,
                               const double *_sc_restrict y,
                               const double *_sc_restrict z,
                               const int n_nodes,
                               const int *_sc_restrict Vmask,
                               rhea_viscosity_options_t *opt)
{
  const int           in_temp = (temp_elem != NULL);
  const int           in_weak = (weak_elem != NULL);
  const int           in_strt = (strt_sqrt_2inv_elem != NULL);
  const double        temp_default = RHEA_TEMPERATURE_NEUTRAL_VALUE;
  const double        weak_default = RHEA_WEAKZONE_NEUTRAL_VALUE;
  const double        strt_default = RHEA_STRAINRATE_2INV_NEUTRAL_VALUE;
  rhea_domain_options_t  *domain_options = opt->domain_options;
  const double        interface_radius =
                        domain_options->lm_um_interface_radius;
  const double        interface_smoothing_width =
                        domain_options->lm_um_interface_smoothing_width;

  double              r1 = 0.0;
  double              bd = RHEA_VISCOSITY_BOUNDS_OFF;
  double              yld = RHEA_VISCOSITY_YIELDING_OFF;
  int                 nodeid;

  /* check input */
  RHEA_ASSERT (visc_elem != NULL);
  RHEA_ASSERT (x != NULL && y != NULL && z != NULL);
  RHEA_ASSERT (Vmask != NULL);

  /* compute viscosity depending on location in lower or upper mantle */
  if (rhea_domain_elem_is_in_upper_mantle (x, y, z, Vmask, domain_options)) {
    for (nodeid = 0; nodeid < n_nodes; nodeid++) {
      const double        temp = (in_temp ? temp_elem[nodeid] : temp_default);
      const double        weak = (in_weak ? weak_elem[nodeid] : weak_default);
      const double        strt = (in_strt ? strt_sqrt_2inv_elem[nodeid] :
                                            strt_default);
      const double        r = rhea_domain_compute_radius (x[nodeid], y[nodeid],
                                                          z[nodeid],
                                                          domain_options);

      /* check temperature for valid range [0,1] */
      RHEA_ASSERT (isfinite (temp));
      RHEA_ASSERT (0.0 <= temp && temp <= 1.0);
      /* check weak zone for valid range (0,1] */
      RHEA_ASSERT (isfinite (weak));
      RHEA_ASSERT (0.0 < weak && weak <= 1.0);
      /* check 2nd invariant of the strain rate for non-negativity */
      RHEA_ASSERT (isfinite (strt));
      RHEA_ASSERT (0.0 <= strt);

      /* compute nonlinear viscosity in upper mantle */
      rhea_viscosity_nonlinear_model (&visc_elem[nodeid], &r1, &bd, &yld,
                                      temp, weak, strt, opt);

      /* update viscosity by applying a smooth transition close to the LM-UM
       * interface (via a convex combination) */
      if (0.0 < interface_smoothing_width &&
          fabs (r - interface_radius) < 0.5 * interface_smoothing_width) {
        const double        v1 = visc_elem[nodeid];
        double              v2, bd2, c;

        rhea_viscosity_linear_model (&v2, &bd2, temp, weak, opt,
                                     0 /* lower mantle */, 1 /* bounds on */);

        c = (fabs (r - interface_radius) - 0.5 * interface_smoothing_width) /
            interface_smoothing_width;
        RHEA_ASSERT (0.0 <= c && c <= 1.0);

        visc_elem[nodeid] = c * v2 + (1.0 - c) * v1;
      }

      /* store values in output arrays if they exist */
      if (rank1_scal_elem != NULL) {
        rank1_scal_elem[nodeid] = r1;
      }
      if (bounds_elem != NULL) {
        bounds_elem[nodeid] = bd;
      }
      if (yielding_elem != NULL) {
        yielding_elem[nodeid] = yld;
      }
    }
  }
  else { /* if element is located in lower mantle */
    const int           restrict_to_bounds = 1;

    for (nodeid = 0; nodeid < n_nodes; nodeid++) {
      const double        temp = (in_temp ? temp_elem[nodeid] : temp_default);
      const double        weak = (in_weak ? weak_elem[nodeid] : weak_default);
      const double        r = rhea_domain_compute_radius (x[nodeid], y[nodeid],
                                                          z[nodeid],
                                                          domain_options);

      /* check temperature for valid range [0,1] */
      RHEA_ASSERT (isfinite (temp));
      RHEA_ASSERT (0.0 <= temp && temp <= 1.0);
      /* check weak zone for valid range (0,1] */
      RHEA_ASSERT (isfinite (weak));
      RHEA_ASSERT (0.0 < weak && weak <= 1.0);

      /* compute linear viscosity in lower mantle */
      rhea_viscosity_linear_model (&visc_elem[nodeid], &bd, temp, weak, opt,
                                   0 /* lower mantle */, 1 /* bounds on */);

      /* update viscosity by applying a smooth transition close to the LM-UM
       * interface (via a convex combination) */
      if (0.0 < interface_smoothing_width &&
          fabs (r - interface_radius) < 0.5 * interface_smoothing_width) {
        const double        strt = (in_strt ? strt_sqrt_2inv_elem[nodeid] :
                                              strt_default);
        const double        v1 = visc_elem[nodeid];
        double              v2, r1_2, bd2, yld2, c;

        RHEA_ASSERT (isfinite (strt));
        RHEA_ASSERT (0.0 <= strt);
        rhea_viscosity_nonlinear_model (&v2, &r1_2, &bd2, &yld2,
                                        temp, weak, strt, opt);

        c = (fabs (r - interface_radius) - 0.5 * interface_smoothing_width) /
            interface_smoothing_width;
        RHEA_ASSERT (0.0 <= c && c <= 1.0);

        visc_elem[nodeid] = c * v1 + (1.0 - c) * v2;
      }

      /* store (default) values in output arrays if they exist */
      if (rank1_scal_elem != NULL) {
        rank1_scal_elem[nodeid] = r1;
      }
      if (bounds_elem != NULL) {
        bounds_elem[nodeid] = bd;
      }
      if (yielding_elem != NULL) {
        yielding_elem[nodeid] = yld;
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
  const int           in_temp = (temp_vec != NULL);
  const int           in_weak = (weak_vec != NULL);
  const int           out_rank1 = (rank1_scal_vec != NULL);
  const int           out_bounds = (bounds_vec != NULL);
  const int           out_yielding = (yielding_vec != NULL);

  sc_dmatrix_t       *temp_el_mat, *weak_el_mat, *vel_el_mat,
                     *strt_sqrt_2inv_el_mat;
  double             *temp_el_data, *weak_el_data, *strt_sqrt_2inv_el_data;
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
  strt_sqrt_2inv_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
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

  temp_el_data           = NULL;
  weak_el_data           = NULL;
  strt_sqrt_2inv_el_data = strt_sqrt_2inv_el_mat->e[0];
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
        strt_sqrt_2inv_el_mat, vel_el_mat, vel_vec, elid,
        tmp_grad_vel, tmp_dvel, tmp_vel);

    /* compute nonlinear viscosity */
    rhea_viscosity_nonlinear_elem (
        visc_el_data, rank1_scal_el_data, bounds_el_data, yielding_el_data,
        temp_el_data, weak_el_data, strt_sqrt_2inv_el_data, x, y, z,
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
  sc_dmatrix_destroy (strt_sqrt_2inv_el_mat);
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

/* callback function for computing the viscosity */
rhea_viscosity_compute_fn_t   rhea_viscosity_compute_fn = NULL;
void                         *rhea_viscosity_compute_fn_data = NULL;

void
rhea_viscosity_set_viscosity_compute_fn (rhea_viscosity_compute_fn_t fn,
                                         void *data)
{
  rhea_viscosity_compute_fn = fn;
  rhea_viscosity_compute_fn_data = data;
}

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
  /* compute custom viscosity and exit */
  if (rhea_viscosity_compute_fn != NULL) {
    rhea_viscosity_compute_fn (viscosity, rank1_tensor_scal, bounds_marker,
                               yielding_marker, temperature, weakzone,
                               velocity, opt, rhea_viscosity_compute_fn_data);
    RHEA_ASSERT (rhea_viscosity_is_valid (viscosity));
    return;
  }

  switch (opt->type) {
  case RHEA_VISCOSITY_LINEAR:
    /* compute linear viscosity and set bounds marker */
    rhea_viscosity_linear_vec (viscosity, bounds_marker,
                               temperature, weakzone, opt);

    /* set default values for vectors pertaining to nonlinear viscosity */
    if (rank1_tensor_scal != NULL) {
      ymir_dvec_set_zero (rank1_tensor_scal);
    }
    if (yielding_marker != NULL) {
      ymir_dvec_set_value (yielding_marker, RHEA_VISCOSITY_YIELDING_OFF);
    }
    break;

  case RHEA_VISCOSITY_NONLINEAR:
    /* compute nonlinear viscosity */
    rhea_viscosity_nonlinear_vec (viscosity, rank1_tensor_scal,
                                  bounds_marker, yielding_marker,
                                  temperature, weakzone, velocity, opt);
    break;

  default: /* unknown viscosity type */
    RHEA_ABORT_NOT_REACHED ();
  }
  RHEA_ASSERT (rhea_viscosity_is_valid (viscosity));
}

/**
 * Computes viscosity to initialize a nonlinear solver with zero velocity.
 */
void
rhea_viscosity_compute_nonlinear_init (ymir_vec_t *viscosity,
                                       ymir_vec_t *rank1_tensor_scal,
                                       ymir_vec_t *bounds_marker,
                                       ymir_vec_t *yielding_marker,
                                       ymir_vec_t *temperature,
                                       ymir_vec_t *weakzone,
                                       rhea_viscosity_options_t *opt)
{
  const rhea_viscosity_nonlinear_init_t  type = opt->type_nonlinear_init;

  switch (type) {
  case RHEA_VISCOSITY_NONLINEAR_INIT_DEFAULT:
    {
      ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (viscosity);
      ymir_vec_t         *velocity = rhea_velocity_new (ymir_mesh);

      /* compute nonlinear viscosity with zero velocity */
      ymir_vec_set_zero (velocity);
      rhea_viscosity_compute (
          viscosity, rank1_tensor_scal, bounds_marker, yielding_marker,
          temperature, weakzone, velocity, opt);

      ymir_vec_destroy (velocity);
    }
    break;

  case RHEA_VISCOSITY_NONLINEAR_INIT_LIN:
  case RHEA_VISCOSITY_NONLINEAR_INIT_LIN_RESCALE_UM:
    {
      const double        upper_mantle_scaling = opt->upper_mantle_scaling;

      /* compute linear viscosity (possibly rescale upper mantle) */
      if (RHEA_VISCOSITY_NONLINEAR_INIT_LIN_RESCALE_UM == type) {
        opt->upper_mantle_scaling = opt->lower_mantle_scaling;
      }
      rhea_viscosity_linear_vec (viscosity, bounds_marker,
                                 temperature, weakzone, opt);
      if (RHEA_VISCOSITY_NONLINEAR_INIT_LIN_RESCALE_UM == type) {
        opt->upper_mantle_scaling = upper_mantle_scaling;
      }

      /* set all other output vectors to zero */
      if (rank1_tensor_scal != NULL) {
        ymir_dvec_set_zero (rank1_tensor_scal);
      }
      if (yielding_marker != NULL) {
        ymir_dvec_set_value (yielding_marker, RHEA_VISCOSITY_YIELDING_OFF);
      }
    }
    break;

  default: /* unknown initial nonlinear viscosity type */
    RHEA_ABORT_NOT_REACHED ();
  }
}

double
rhea_viscosity_get_visc_shift (rhea_viscosity_options_t *opt)
{
  return opt->shift;
}

double
rhea_viscosity_get_visc_shift_proj (rhea_viscosity_options_t *opt)
{
  double              shift_proj;

  if (opt->type == RHEA_VISCOSITY_NONLINEAR) {
    switch (opt->model) {
    case RHEA_VISCOSITY_MODEL_UWYL:
      shift_proj = 0.0;
      break;
    case RHEA_VISCOSITY_MODEL_UWYL_LADD_UCUT:
    case RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT:
      if (0 < opt->min) {
        shift_proj = -2.0 * opt->min;
      }
      else {
        shift_proj = 0.0;
      }
      break;
    default: /* unknown viscosity model */
      RHEA_ABORT_NOT_REACHED ();
    }
  }

  return shift_proj;
}

int
rhea_viscosity_has_strain_rate_weakening (rhea_viscosity_options_t *opt)
{
  const int           is_valid = (isfinite (opt->stress_exponent) &&
                                  0.0 < opt->stress_exponent &&
                                  SC_EPS < fabs (1.0 - opt->stress_exponent));

  switch (opt->type_nonlinear) {
  case RHEA_VISCOSITY_NONLINEAR_SRW:      return is_valid;
  case RHEA_VISCOSITY_NONLINEAR_YLD:      return 0;
  case RHEA_VISCOSITY_NONLINEAR_SRW_YLD:  return is_valid;
  default: /* unknown nonlinear viscosity type */
    RHEA_ABORT_NOT_REACHED ();
  }
}

int
rhea_viscosity_has_yielding (rhea_viscosity_options_t *opt)
{
  const int           is_valid = (isfinite (opt->yield_strength) &&
                                  0.0 < opt->yield_strength);

  switch (opt->type_nonlinear) {
  case RHEA_VISCOSITY_NONLINEAR_SRW:      return 0;
  case RHEA_VISCOSITY_NONLINEAR_YLD:      return is_valid;
  case RHEA_VISCOSITY_NONLINEAR_SRW_YLD:  return is_valid;
  default: /* unknown nonlinear viscosity type */
    RHEA_ABORT_NOT_REACHED ();
  }
}

double
rhea_viscosity_get_yield_strength (rhea_viscosity_options_t *opt)
{
  if (rhea_viscosity_has_yielding (opt)) {
    return opt->yield_strength;
  }
  else {
    return NAN;
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

/******************************************************************************
 * Statistics
 *****************************************************************************/

/**
 * Computes the volume of a filter.  A filter is understood as a vector with
 * ones where the filter is active and zeros otherwise.
 */
static double
rhea_viscosity_stats_filter_compute_volume (ymir_vec_t *filter)
{
  ymir_vec_t         *unit = ymir_vec_template (filter);
  ymir_vec_t         *mass_out = ymir_vec_template (filter);
  double              vol;

  /* check input */
  RHEA_ASSERT (ymir_vec_is_dvec (filter) && filter->ndfields == 1);

  /* set unit vector */
  ymir_vec_set_value (unit, 1.0);

  /* integrate to get volume */
  ymir_mass_apply (filter, mass_out);
  vol = ymir_vec_innerprod (unit, mass_out);

  /* destroy */
  ymir_vec_destroy (unit);
  ymir_vec_destroy (mass_out);

  /* return volume of filter */
  return vol;
}

void
rhea_viscosity_stats_get_bounds_volume (double *vol_min, double *vol_max,
                                        ymir_vec_t *bounds_marker)
{
  ymir_vec_t         *bounds = ymir_vec_template (bounds_marker);

  /* compute volume of active min bounds */
  if (vol_min != NULL) {
    ymir_vec_copy (bounds_marker, bounds);
    ymir_vec_bound_max (bounds, RHEA_VISCOSITY_BOUNDS_OFF);
    ymir_vec_scale (1.0/RHEA_VISCOSITY_BOUNDS_MIN, bounds);
    *vol_min = rhea_viscosity_stats_filter_compute_volume (bounds);
  }

  /* compute volume of active max bounds */
  if (vol_max != NULL) {
    ymir_vec_copy (bounds_marker, bounds);
    ymir_vec_shift (-1.0/RHEA_VISCOSITY_BOUNDS_MAX_WEAK, bounds);
    ymir_vec_bound_min (bounds, RHEA_VISCOSITY_BOUNDS_OFF);
    ymir_vec_scale (
        1.0/(RHEA_VISCOSITY_BOUNDS_MAX - RHEA_VISCOSITY_BOUNDS_MAX_WEAK),
        bounds);
    *vol_max = rhea_viscosity_stats_filter_compute_volume (bounds);
  }

  /* destroy */
  ymir_vec_destroy (bounds);
}

double
rhea_viscosity_stats_get_yielding_volume (ymir_vec_t *yielding_marker)
{
  RHEA_ASSERT (fabs (1.0 - RHEA_VISCOSITY_YIELDING_ACTIVE) <= 0.0);
  return rhea_viscosity_stats_filter_compute_volume (yielding_marker);
}
