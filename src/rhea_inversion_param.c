#include <rhea_inversion_param.h>
#include <rhea_base.h>
#include <rhea_weakzone_label.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_INVERSION_PARAM_DEFAULT_ACTIVATE 1
#define RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE 0

/* global options */
rhea_inversion_param_options_t rhea_inversion_param_options;

void
rhea_inversion_param_add_options (
                              rhea_inversion_param_options_t *inv_param_options,
                              ymir_options_t *opt_sup)
{
  const char         *opt_prefix = "Parameters";
  ymir_options_t     *opt = ymir_options_new ();
  rhea_inversion_param_options_t *inv_param_opt;

  /* set options storage */
  if (inv_param_options != NULL) {
    /* choose provided options */
    inv_param_opt = inv_param_options;
  }
  else {
    /* choose global options */
    inv_param_opt = &rhea_inversion_param_options;
  }

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_B, "activate-min", '\0',
    &(inv_param_opt->min_a), RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate min viscosity bound",
  YMIR_OPTIONS_B, "activate-max", '\0',
    &(inv_param_opt->max_a), RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate max viscosity bound",

  YMIR_OPTIONS_B, "activate-upper-mantle-scaling", '\0',
    &(inv_param_opt->upper_mantle_scaling_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate scaling factor of upper mantle",
  YMIR_OPTIONS_B, "activate-lower-mantle-scaling", '\0',
    &(inv_param_opt->lower_mantle_scaling_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate scaling factor of lower mantle",

  YMIR_OPTIONS_B, "activate-upper-mantle-arrhenius-activation-energy", '\0',
    &(inv_param_opt->upper_mantle_arrhenius_activation_energy_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate the activation energy in Arrhenius relationship of upper mantle",
  YMIR_OPTIONS_B, "activate-lower-mantle-arrhenius-activation-energy", '\0',
    &(inv_param_opt->lower_mantle_arrhenius_activation_energy_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate the activation energy in Arrhenius relationship of lower mantle",

  YMIR_OPTIONS_B, "activate-stress-exponent", '\0',
    &(inv_param_opt->stress_exponent_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate the stress exponent that governs strain rate weakening",
  YMIR_OPTIONS_B, "activate-yield-strength", '\0',
    &(inv_param_opt->yield_strength_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate plastic yielding",

  YMIR_OPTIONS_B, "activate-weakzone-thickness", '\0',
    &(inv_param_opt->thickness_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone thickness",
  YMIR_OPTIONS_B, "activate-weakzone-thickness-generic-slab", '\0',
    &(inv_param_opt->thickness_generic_slab_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone thickness of slabs",
  YMIR_OPTIONS_B, "activate-weakzone-thickness-generic-ridge", '\0',
    &(inv_param_opt->thickness_generic_ridge_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone thickness of slabs",
  YMIR_OPTIONS_B, "activate-weakzone-thickness-generic-fracture", '\0',
    &(inv_param_opt->thickness_generic_fracture_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone thickness of slabs",

  YMIR_OPTIONS_B, "activate-weakzone-thickness-const", '\0',
    &(inv_param_opt->thickness_const_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone interior thickness (where min factor is reached)",
  YMIR_OPTIONS_B, "activate-weakzone-thickness-const-generic-slab", '\0',
    &(inv_param_opt->thickness_const_generic_slab_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate interior thickness for slabs (where min factor is reached)",
  YMIR_OPTIONS_B, "activate-weakzone-thickness-const-generic-ridge", '\0',
    &(inv_param_opt->thickness_const_generic_ridge_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate interior thickness for ridges (where min factor is reached)",
  YMIR_OPTIONS_B, "activate-weakzone-thickness-const-generic-fracture", '\0',
    &(inv_param_opt->thickness_const_generic_fracture_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate interior thickness for fractures (where min factor is reached)",

  YMIR_OPTIONS_B, "activate-weak-factor-interior", '\0',
    &(inv_param_opt->weak_factor_interior_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone factor",
  YMIR_OPTIONS_B, "activate-weak-factor-interior-generic-slab", '\0',
    &(inv_param_opt->weak_factor_interior_generic_slab_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone factor of generic slabs",
  YMIR_OPTIONS_B, "activate-weak-factor-interior-generic-ridge", '\0',
    &(inv_param_opt->weak_factor_interior_generic_ridge_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone factor of generic ridges",
  YMIR_OPTIONS_B, "activate-weak-factor-interior-generic-fracture", '\0',
    &(inv_param_opt->weak_factor_interior_generic_fracture_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone factor of generic fractures",

  YMIR_OPTIONS_B, "activate-weak-factor-interior-earth-slab", '\0',
    &(inv_param_opt->weak_factor_interior_earth_slab_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone factor of earth's slabs",
  YMIR_OPTIONS_B, "activate-weak-factor-interior-earth-ridge", '\0',
    &(inv_param_opt->weak_factor_interior_earth_ridge_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone factor of earth's ridges",
  YMIR_OPTIONS_B, "activate-weak-factor-interior-earth-fracture", '\0',
    &(inv_param_opt->weak_factor_interior_earth_fracture_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone factor of earth's fractures",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);
}

/******************************************************************************
 * Indices for Inversion Parameters
 *****************************************************************************/

/* indices of parameter vector */
typedef enum
{
  /*
   * Parameterization of Viscosity Options
   */

  /* parameterization of lower and upper bounds for the viscosity */
  RHEA_INVERSION_PARAM_VISC_MIN,
  RHEA_INVERSION_PARAM_VISC_MAX,

  /* parameterization of scaling factors */
  RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING,
  RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING,

  /* parameterization of activation energy in Arrhenius relationship */
  RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY,
  RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY,

  /* parameterization of stress exponent that governs strain rate weakening */
  RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT,

  /* parameterization of plastic yielding */
  RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH,

  /*
   * Parameterization of Weak Zone Options
   */

  /* parameterization of weak zone geometry */
  RHEA_INVERSION_PARAM_WEAK_THICKNESS,
  RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_SLAB,
  RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_RIDGE,
  RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_FRACTURE,
  RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST,
  RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_SLAB,
  RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_RIDGE,
  RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_FRACTURE,

  /* parameterization of generic max weakening in the interior of weak zones */
  RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR,
  RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_SLAB,
  RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_RIDGE,
  RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_FRACTURE,

  /* max number of global parameters (is greater than actual #parameters) */
  RHEA_INVERSION_PARAM_N_GLOBAL = 32,

  /* parameterization of earth's max weakening in the interior of weak zones */
  RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB,
  RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_RIDGE =
    RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB +
    RHEA_WEAKZONE_LABEL_EARTH_N_SL,
  RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_FRACTURE =
    RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_RIDGE +
    RHEA_WEAKZONE_LABEL_EARTH_N_RI,

  /* total number of global and local (i.e., weak zone) parameters */
  RHEA_INVERSION_PARAM_N =
    RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_FRACTURE +
    RHEA_WEAKZONE_LABEL_EARTH_N_FZ
}
rhea_inversion_param_idx_t;

/******************************************************************************
 * Parameter Activation
 *****************************************************************************/

static int *
rhea_inversion_param_activation_mask_new (
                              rhea_inversion_param_options_t *inv_param_options)
{
  int                *active = RHEA_ALLOC_ZERO (int, RHEA_INVERSION_PARAM_N);

  /* set activation mask of viscosity parameters */
  if (inv_param_options->min_a) {
    active[RHEA_INVERSION_PARAM_VISC_MIN] = 1;
  }
  if (inv_param_options->max_a) {
    active[RHEA_INVERSION_PARAM_VISC_MAX] = 1;
  }
  if (inv_param_options->upper_mantle_scaling_a) {
    active[RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING] = 1;
  }
  if (inv_param_options->lower_mantle_scaling_a) {
    active[RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING] = 1;
  }
  if (inv_param_options->upper_mantle_arrhenius_activation_energy_a) {
    active[RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY] = 1;
  }
  if (inv_param_options->lower_mantle_arrhenius_activation_energy_a) {
    active[RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY] = 1;
  }
  if (inv_param_options->stress_exponent_a) {
    active[RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT] = 1;
  }
  if (inv_param_options->yield_strength_a) {
    active[RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH] = 1;
  }

  /* set activation mask of weak zone thickness parameters */
  if (inv_param_options->thickness_a) {
    active[RHEA_INVERSION_PARAM_WEAK_THICKNESS] = 1;

    if (inv_param_options->thickness_generic_slab_a) {
      active[RHEA_INVERSION_PARAM_WEAK_THICKNESS] = 0;
      active[RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_SLAB] = 1;
    }
    if (inv_param_options->thickness_generic_ridge_a) {
      active[RHEA_INVERSION_PARAM_WEAK_THICKNESS] = 0;
      active[RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_RIDGE] = 1;
    }
    if (inv_param_options->thickness_generic_fracture_a) {
      active[RHEA_INVERSION_PARAM_WEAK_THICKNESS] = 0;
      active[RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_FRACTURE] = 1;
    }
  }
  if (inv_param_options->thickness_const_a) {
    active[RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST] = 1;

    if (inv_param_options->thickness_const_generic_slab_a) {
      active[RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST] = 0;
      active[RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_SLAB] = 1;
    }
    if (inv_param_options->thickness_const_generic_ridge_a) {
      active[RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST] = 0;
      active[RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_RIDGE] = 1;
    }
    if (inv_param_options->thickness_const_generic_fracture_a) {
      active[RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST] = 0;
      active[RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_FRACTURE] = 1;
    }
  }

  /* set activation mask of weak factors */
  if (inv_param_options->weak_factor_interior_a) {
    int                 i;

    active[RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR] = 1;

    if (inv_param_options->weak_factor_interior_earth_slab_a) {
      active[RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR] = 0;
      for (i = 0; i < RHEA_WEAKZONE_LABEL_EARTH_N_SL; i++) {
        active[RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB+i] = 1;
      }
    }
    else if (inv_param_options->weak_factor_interior_generic_slab_a) {
      active[RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR] = 0;
      active[RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_SLAB] = 1;
    }

    if (inv_param_options->weak_factor_interior_earth_ridge_a) {
      active[RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR] = 0;
      for (i = 0; i < RHEA_WEAKZONE_LABEL_EARTH_N_RI; i++) {
        active[RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_RIDGE+i] = 1;
      }
    }
    else if (inv_param_options->weak_factor_interior_generic_ridge_a) {
      active[RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR] = 0;
      active[RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_RIDGE] = 1;
    }

    if (inv_param_options->weak_factor_interior_earth_fracture_a) {
      active[RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR] = 0;
      for (i = 0; i < RHEA_WEAKZONE_LABEL_EARTH_N_FZ; i++) {
        active[RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_FRACTURE+i] = 1;
      }
    }
    else if (inv_param_options->weak_factor_interior_generic_fracture_a) {
      active[RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR] = 0;
      active[RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_FRACTURE] = 1;
    }
  }

  return active;
}

static void
rhea_inversion_param_activation_mask_destroy (int *active)
{
  RHEA_FREE (active);
}

static int
rhea_inversion_param_activation_mask_count (int *active)
{
  int                 n_active = 0;
  int                 i;

  for (i = 0; i < RHEA_INVERSION_PARAM_N; i++) {
    if (active[i]) {
      n_active++;
    }
  }

  return n_active;
}

/******************************************************************************
 * Parameters
 *****************************************************************************/

/* inversion parameters */
struct rhea_inversion_param
{
  /* vector containing all inversion parameters */
  ymir_vec_t         *parameter_vec;
  int                 n_parameters;

  /* activation mask for inversion parameters */
  int                *active;
  int                 n_active;

  //TODO add weights?

  /* options (not owned) */
  rhea_weakzone_options_t    *weak_options;
  rhea_viscosity_options_t   *visc_options;
};

rhea_inversion_param_t *
rhea_inversion_param_new (rhea_weakzone_options_t *weak_options,
                          rhea_viscosity_options_t *visc_options,
                          rhea_inversion_param_options_t *inv_param_options)
{
  rhea_inversion_param_options_t *inv_param_opt;
  rhea_inversion_param_t         *inv_param;

  /* set options storage */
  if (inv_param_options != NULL) {
    /* choose provided options */
    inv_param_opt = inv_param_options;
  }
  else {
    /* choose global options */
    inv_param_opt = &rhea_inversion_param_options;
  }

  /* initialize inversion parameters */
  inv_param = RHEA_ALLOC (rhea_inversion_param_t, 1);
  inv_param->n_parameters = RHEA_INVERSION_PARAM_N;
  inv_param->parameter_vec = ymir_vec_new_meshfree (inv_param->n_parameters);
  inv_param->weak_options = weak_options;
  inv_param->visc_options = visc_options;

  /* determine active parameters */
  inv_param->active =
    rhea_inversion_param_activation_mask_new (inv_param_opt);
  inv_param->n_active =
    rhea_inversion_param_activation_mask_count (inv_param->active);
  RHEA_ASSERT (0 < inv_param->n_active);

  /* set up parameter vector */
  {
    double             *param_data = inv_param->parameter_vec->meshfree->e[0];
    int                 i;

    for (i = 0; i < inv_param->n_parameters; i++) {
      param_data[i] = NAN;
    }
  }
  rhea_inversion_param_pull_from_model (inv_param);

  /* return inversion parameters */
  return inv_param;
}

void
rhea_inversion_param_destroy (rhea_inversion_param_t *inv_param)
{
  ymir_vec_destroy (inv_param->parameter_vec);
  rhea_inversion_param_activation_mask_destroy (inv_param->active);
  RHEA_FREE (inv_param);
}

static int
rhea_inversion_param_calculate_inversion_val_if_active_exp (
                                            const double model_val,
                                            const int param_idx,
                                            rhea_inversion_param_t *inv_param)
{
  double             *param_data = inv_param->parameter_vec->meshfree->e[0];

  if (inv_param->active[param_idx]) {
    param_data[param_idx] = log (model_val);
    return 1;
  }
  else {
    return 0;
  }
}

static int
rhea_inversion_param_calculate_inversion_val_if_active_exp2 (
                                            const double model_val,
                                            const int param_idx,
                                            rhea_inversion_param_t *inv_param)
{
  double             *param_data = inv_param->parameter_vec->meshfree->e[0];

  if (inv_param->active[param_idx]) {
    RHEA_ASSERT (0.0 < model_val && model_val <= 1.0);
    param_data[param_idx] = sqrt (-log (model_val));
    return 1;
  }
  else {
    return 0;
  }
}

void
rhea_inversion_param_pull_from_model (rhea_inversion_param_t *inv_param)
{
  rhea_weakzone_options_t  *weak_options = inv_param->weak_options;
  rhea_viscosity_options_t *visc_options = inv_param->visc_options;

  /*
   * Pull from Viscosity Options
   */

  /* pull lower and upper bounds of the viscosity:
   *   bound_p = log (bound)
   */
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      visc_options->min, RHEA_INVERSION_PARAM_VISC_MIN, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      visc_options->max, RHEA_INVERSION_PARAM_VISC_MAX, inv_param);

  /* pull scaling factors:
   *   scaling_p = log (scaling)
   */
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      visc_options->upper_mantle_scaling,
      RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      visc_options->lower_mantle_scaling,
      RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING, inv_param);

  /* pull activation energy in Arrhenius relationship:
   *   arrhenius_activation_energy_p = log (arrhenius_activation_energy)
   */
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      visc_options->upper_mantle_arrhenius_activation_energy,
      RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      visc_options->lower_mantle_arrhenius_activation_energy,
      RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY, inv_param);

  /* pull stress exponent that governs strain rate weakening:
   *   stress_exponent_p = log (stress_exponent - 1)
   */
  if (inv_param->active[RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT]) {
    double             *param_data = inv_param->parameter_vec->meshfree->e[0];

    RHEA_ASSERT (1.0 < visc_options->stress_exponent);
    param_data[RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT] =
      log (visc_options->stress_exponent - 1.0);
  }

  /* pull yield strength:
   *   yield_strength_p = log (yield_strength)
   */
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      visc_options->yield_strength,
      RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH, inv_param);

  /*
   * Pull from Weak Zone Options
   */

  /* pull weak zone geometry:
   *   thickness_p = log (thickness)
   */
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      weak_options->thickness,
      RHEA_INVERSION_PARAM_WEAK_THICKNESS, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      weak_options->thickness_generic_slab,
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_SLAB, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      weak_options->thickness_generic_ridge,
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_RIDGE, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      weak_options->thickness_generic_fracture,
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_FRACTURE, inv_param);

  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      weak_options->thickness_const,
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      weak_options->thickness_const_generic_slab,
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_SLAB, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      weak_options->thickness_const_generic_ridge,
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_RIDGE, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      weak_options->thickness_const_generic_fracture,
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_FRACTURE, inv_param);

  /* pull generic max weakening in the interior of weak zones:
   *   weak_factor_interior_p = sqrt (-log (weak_factor_interior))
   */
  rhea_inversion_param_calculate_inversion_val_if_active_exp2 (
      weak_options->weak_factor_interior,
      RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp2 (
      weak_options->weak_factor_interior_generic_slab,
      RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_SLAB, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp2 (
      weak_options->weak_factor_interior_generic_ridge,
      RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_RIDGE, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp2 (
      weak_options->weak_factor_interior_generic_fracture,
      RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_FRACTURE, inv_param);

  /* check earth's weak factors */
#ifdef RHEA_ENABLE_DEBUG
  {
    int                 offset, i;
    int                 n = 0;

    offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB;
    for (i = 0; i < RHEA_WEAKZONE_LABEL_EARTH_N; i++) {
      n += inv_param->active[offset+i];
    }

    RHEA_ASSERT (!n || weak_options->weak_factor_interior_earth != NULL);
  }
#endif

  /* pull earth's max weakening in the interior of weak zones:
   *   weak_factor_interior_p = sqrt (-log (weak_factor_interior))
   */
  {
    int                 offset, i;

    offset = 0;
    for (i = 0; i < RHEA_WEAKZONE_LABEL_EARTH_N_SL; i++) {
      rhea_inversion_param_calculate_inversion_val_if_active_exp2 (
          weak_options->weak_factor_interior_earth[offset+i],
          RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB+i, inv_param);
    }

    offset = RHEA_WEAKZONE_LABEL_EARTH_N_SL;
    for (i = 0; i < RHEA_WEAKZONE_LABEL_EARTH_N_RI; i++) {
      rhea_inversion_param_calculate_inversion_val_if_active_exp2 (
          weak_options->weak_factor_interior_earth[offset+i],
          RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_RIDGE+i, inv_param);
    }

    offset = RHEA_WEAKZONE_LABEL_EARTH_N_SL +
             RHEA_WEAKZONE_LABEL_EARTH_N_RI;
    for (i = 0; i < RHEA_WEAKZONE_LABEL_EARTH_N_FZ; i++) {
      rhea_inversion_param_calculate_inversion_val_if_active_exp2 (
          weak_options->weak_factor_interior_earth[offset+i],
          RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_FRACTURE+i,
          inv_param);
    }
  }
}

static int
rhea_inversion_param_calculate_model_val_if_active_exp (
                                            double *model_val,
                                            const int param_idx,
                                            rhea_inversion_param_t *inv_param)
{
  double             *param_data = inv_param->parameter_vec->meshfree->e[0];

  if (inv_param->active[param_idx]) {
    *model_val = exp (param_data[param_idx]);
    return 1;
  }
  else {
    return 0;
  }
}

static int
rhea_inversion_param_calculate_model_val_if_active_exp2 (
                                            double *model_val,
                                            const int param_idx,
                                            rhea_inversion_param_t *inv_param)
{
  double             *param_data = inv_param->parameter_vec->meshfree->e[0];

  if (inv_param->active[param_idx]) {
    *model_val = exp (-param_data[param_idx]*param_data[param_idx]);
    return 1;
  }
  else {
    return 0;
  }
}

void
rhea_inversion_param_push_to_model (rhea_inversion_param_t *inv_param)
{
  rhea_weakzone_options_t  *weak_options = inv_param->weak_options;
  rhea_viscosity_options_t *visc_options = inv_param->visc_options;

  /*
   * Push to Viscosity Options
   */

  /* push lower and upper bounds of the viscosity:
   *   bound = exp (bound_p)
   */
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(visc_options->min), RHEA_INVERSION_PARAM_VISC_MIN, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(visc_options->max), RHEA_INVERSION_PARAM_VISC_MAX, inv_param);

  /* push scaling factors:
   *   scaling = exp (scaling_p)
   */
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(visc_options->upper_mantle_scaling),
      RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(visc_options->lower_mantle_scaling),
      RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING, inv_param);

  /* push activation energy in Arrhenius relationship:
   *   arrhenius_activation_energy = exp (arrhenius_activation_energy_p)
   */
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(visc_options->upper_mantle_arrhenius_activation_energy),
      RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(visc_options->lower_mantle_arrhenius_activation_energy),
      RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY, inv_param);

  /* push stress exponent that governs strain rate weakening:
   *   stress_exponent = 1 + exp (stress_exponent_p)
   */
  if (inv_param->active[RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT]) {
    double             *param_data = inv_param->parameter_vec->meshfree->e[0];

    visc_options->stress_exponent =
      1.0 + exp (param_data[RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT]);
  }

  /* push yield strength:
   *   yield_strength = exp (yield_strength_p)
   */
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(visc_options->yield_strength),
      RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH, inv_param);

  /*
   * Push to Weak Zone Options
   */

  /* push weak zone geometry:
   *   thickness = exp (thickness_p)
   */
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(weak_options->thickness),
      RHEA_INVERSION_PARAM_WEAK_THICKNESS, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(weak_options->thickness_generic_slab),
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_SLAB, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(weak_options->thickness_generic_ridge),
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_RIDGE, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(weak_options->thickness_generic_fracture),
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_FRACTURE, inv_param);

  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(weak_options->thickness_const),
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(weak_options->thickness_const_generic_slab),
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_SLAB, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(weak_options->thickness_const_generic_ridge),
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_RIDGE, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(weak_options->thickness_const_generic_fracture),
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_FRACTURE, inv_param);

  /* push generic max weakening in the interior of weak zones:
   *   weak_factor_interior = exp (-weak_factor_interior_p^2)
   */
  rhea_inversion_param_calculate_model_val_if_active_exp2 (
      &(weak_options->weak_factor_interior),
      RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp2 (
      &(weak_options->weak_factor_interior_generic_slab),
      RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_SLAB, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp2 (
      &(weak_options->weak_factor_interior_generic_ridge),
      RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_RIDGE, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp2 (
      &(weak_options->weak_factor_interior_generic_fracture),
      RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_FRACTURE, inv_param);

  /* check earth's weak factors */
#ifdef RHEA_ENABLE_DEBUG
  {
    int                 offset, i;
    int                 n = 0;

    offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB;
    for (i = 0; i < RHEA_WEAKZONE_LABEL_EARTH_N; i++) {
      n += inv_param->active[offset+i];
    }

    RHEA_ASSERT (!n || weak_options->weak_factor_interior_earth != NULL);
  }
#endif

  /* push earth's max weakening in the interior of weak zones:
   *   weak_factor_interior = exp (-weak_factor_interior_p^2)
   */
  {
    int                 offset, i;

    offset = 0;
    for (i = 0; i < RHEA_WEAKZONE_LABEL_EARTH_N_SL; i++) {
      rhea_inversion_param_calculate_model_val_if_active_exp2 (
          &(weak_options->weak_factor_interior_earth[offset+i]),
          RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB+i, inv_param);
    }

    offset = RHEA_WEAKZONE_LABEL_EARTH_N_SL;
    for (i = 0; i < RHEA_WEAKZONE_LABEL_EARTH_N_RI; i++) {
      rhea_inversion_param_calculate_model_val_if_active_exp2 (
          &(weak_options->weak_factor_interior_earth[offset+i]),
          RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_RIDGE+i, inv_param);
    }

    offset = RHEA_WEAKZONE_LABEL_EARTH_N_SL +
             RHEA_WEAKZONE_LABEL_EARTH_N_RI;
    for (i = 0; i < RHEA_WEAKZONE_LABEL_EARTH_N_FZ; i++) {
      rhea_inversion_param_calculate_model_val_if_active_exp2 (
          &(weak_options->weak_factor_interior_earth[offset+i]),
          RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_FRACTURE+i,
          inv_param);
    }
  }
}

ymir_vec_t *
rhea_inversion_param_get_vector (rhea_inversion_param_t *inv_param)
{
  return inv_param->parameter_vec;
}
