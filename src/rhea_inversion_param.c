#include <rhea_inversion_param.h>
#include <rhea_base.h>
#include <rhea_viscosity_param_derivative.h>
#include <rhea_weakzone_label.h>
#include <ymir_stress_pc.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_INVERSION_PARAM_DEFAULT_ACTIVATE 1
#define RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE 0

/* global options */
rhea_inversion_param_options_t rhea_inversion_param_options_global;

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
    inv_param_opt = &rhea_inversion_param_options_global;
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
                              rhea_inversion_param_options_t *inv_param_options,
                              rhea_weakzone_options_t *weak_options,
                              rhea_viscosity_options_t *visc_options)
{
  const int           weak_exists = rhea_weakzone_exists (weak_options);
  int                *active = RHEA_ALLOC_ZERO (int, RHEA_INVERSION_PARAM_N);

  /* set activation mask of viscosity parameters */
  if (rhea_viscosity_restrict_min (visc_options) && inv_param_options->min_a) {
    active[RHEA_INVERSION_PARAM_VISC_MIN] = 1;
  }
  if (rhea_viscosity_restrict_max (visc_options) && inv_param_options->max_a) {
    active[RHEA_INVERSION_PARAM_VISC_MAX] = 1;
  }
  if (inv_param_options->upper_mantle_scaling_a) {
    active[RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING] = 1;
  }
  if (inv_param_options->lower_mantle_scaling_a) {
    active[RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING] = 1;
  }
  if (rhea_viscosity_has_arrhenius (visc_options) &&
      inv_param_options->upper_mantle_arrhenius_activation_energy_a) {
    active[RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY] = 1;
  }
  if (rhea_viscosity_has_arrhenius (visc_options) &&
      inv_param_options->lower_mantle_arrhenius_activation_energy_a) {
    active[RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY] = 1;
  }
  if (rhea_viscosity_has_strain_rate_weakening (visc_options) &&
      inv_param_options->stress_exponent_a) {
    active[RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT] = 1;
  }
  if (rhea_viscosity_has_yielding (visc_options) &&
      inv_param_options->yield_strength_a) {
    active[RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH] = 1;
  }

  /* set activation mask of weak zone thickness parameters */
  if (weak_exists && inv_param_options->thickness_a) {
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
  if (weak_exists && inv_param_options->thickness_const_a) {
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
  if (weak_exists && inv_param_options->weak_factor_interior_a) {
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
  /* activation mask for inversion parameters */
  int                *active;
  int                 n_active;
  int                 n_parameters;

  //TODO add weights?

  /* Stokes problem (not owned) */
  rhea_stokes_problem_t  *stokes_problem;

  /* options (not owned) */
  rhea_weakzone_options_t  *weak_options;
  rhea_viscosity_options_t *visc_options;
};

rhea_inversion_param_t *
rhea_inversion_param_new (rhea_stokes_problem_t *stokes_problem,
                          rhea_inversion_param_options_t *inv_param_options)
{
  rhea_inversion_param_options_t *inv_param_opt;
  rhea_inversion_param_t         *inv_param;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* set options storage */
  if (inv_param_options != NULL) {
    /* choose provided options */
    inv_param_opt = inv_param_options;
  }
  else {
    /* choose global options */
    inv_param_opt = &rhea_inversion_param_options_global;
  }

  /* initialize inversion parameters */
  inv_param = RHEA_ALLOC (rhea_inversion_param_t, 1);
  inv_param->n_parameters = RHEA_INVERSION_PARAM_N;
  inv_param->stokes_problem = stokes_problem;
  inv_param->weak_options =
    rhea_stokes_problem_get_weakzone_options (stokes_problem);
  inv_param->visc_options =
    rhea_stokes_problem_get_viscosity_options (stokes_problem);

  /* determine active parameters */
  inv_param->active = rhea_inversion_param_activation_mask_new (
      inv_param_opt, inv_param->weak_options, inv_param->visc_options);
  inv_param->n_active = rhea_inversion_param_activation_mask_count (
      inv_param->active);

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);

  /* return inversion parameters */
  return inv_param;
}

void
rhea_inversion_param_destroy (rhea_inversion_param_t *inv_param)
{
  rhea_inversion_param_activation_mask_destroy (inv_param->active);
  RHEA_FREE (inv_param);
}

static int
rhea_inversion_param_calculate_inversion_val_if_active_exp (
                                            const double model_val,
                                            const int param_idx,
                                            ymir_vec_t *parameter_vec,
                                            rhea_inversion_param_t *inv_param)
{
  double             *param = parameter_vec->meshfree->e[0];

  if (inv_param->active[param_idx]) {
    param[param_idx] = log (model_val);
    RHEA_ASSERT (isfinite (param[param_idx]));
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
                                            ymir_vec_t *parameter_vec,
                                            rhea_inversion_param_t *inv_param)
{
  double             *param = parameter_vec->meshfree->e[0];

  if (inv_param->active[param_idx]) {
    RHEA_ASSERT (0.0 < model_val && model_val <= 1.0);
    param[param_idx] = sqrt (-log (model_val));
    RHEA_ASSERT (isfinite (param[param_idx]));
    return 1;
  }
  else {
    return 0;
  }
}

void
rhea_inversion_param_pull_from_model (ymir_vec_t *parameter_vec,
                                      rhea_inversion_param_t *inv_param)
{
  rhea_weakzone_options_t  *weak_options = inv_param->weak_options;
  rhea_viscosity_options_t *visc_options = inv_param->visc_options;
  int                 pull_weak_earth;
  int                 offset, i;

  RHEA_GLOBAL_VERBOSEF_FN_TAG (__func__, "n_active=%i", inv_param->n_active);

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (parameter_vec, inv_param));

  /* initialize to zero */
  ymir_vec_set_zero (parameter_vec);

  /*
   * Pull from Viscosity Options
   */

  /* pull lower and upper bounds of the viscosity:
   *   bound_p = log (bound)
   */
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      visc_options->min, RHEA_INVERSION_PARAM_VISC_MIN,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      visc_options->max, RHEA_INVERSION_PARAM_VISC_MAX,
      parameter_vec, inv_param);

  /* pull scaling factors:
   *   scaling_p = log (scaling)
   */
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      visc_options->upper_mantle_scaling,
      RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      visc_options->lower_mantle_scaling,
      RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING,
      parameter_vec, inv_param);

  /* pull activation energy in Arrhenius relationship:
   *   arrhenius_activation_energy_p = log (arrhenius_activation_energy)
   */
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      visc_options->upper_mantle_arrhenius_activation_energy,
      RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      visc_options->lower_mantle_arrhenius_activation_energy,
      RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY,
      parameter_vec, inv_param);

  /* pull stress exponent that governs strain rate weakening:
   *   stress_exponent_p = log (stress_exponent - 1)
   */
  if (inv_param->active[RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT]) {
    double             *param = parameter_vec->meshfree->e[0];

    RHEA_ASSERT (1.0 < visc_options->stress_exponent);
    param[RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT] =
      log (visc_options->stress_exponent - 1.0);
    RHEA_ASSERT (isfinite (param[RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT]));
  }

  /* pull yield strength:
   *   yield_strength_p = log (yield_strength)
   */
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      visc_options->yield_strength,
      RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH, parameter_vec, inv_param);

  /*
   * Pull from Weak Zone Options
   */

  /* pull weak zone geometry:
   *   thickness_p = log (thickness)
   */
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      weak_options->thickness,
      RHEA_INVERSION_PARAM_WEAK_THICKNESS,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      weak_options->thickness_generic_slab,
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_SLAB,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      weak_options->thickness_generic_ridge,
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_RIDGE,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      weak_options->thickness_generic_fracture,
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_FRACTURE,
      parameter_vec, inv_param);

  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      weak_options->thickness_const,
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      weak_options->thickness_const_generic_slab,
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_SLAB,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      weak_options->thickness_const_generic_ridge,
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_RIDGE,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp (
      weak_options->thickness_const_generic_fracture,
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_FRACTURE,
      parameter_vec, inv_param);

  /* pull generic max weakening in the interior of weak zones:
   *   weak_factor_interior_p = sqrt (-log (weak_factor_interior))
   */
  rhea_inversion_param_calculate_inversion_val_if_active_exp2 (
      weak_options->weak_factor_interior,
      RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp2 (
      weak_options->weak_factor_interior_generic_slab,
      RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_SLAB,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp2 (
      weak_options->weak_factor_interior_generic_ridge,
      RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_RIDGE,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_inversion_val_if_active_exp2 (
      weak_options->weak_factor_interior_generic_fracture,
      RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_FRACTURE,
      parameter_vec, inv_param);

  /* check if earth's weak factors are active */
  {
    int                 n = 0;

    offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB;
    for (i = 0; i < RHEA_WEAKZONE_LABEL_EARTH_N; i++) {
      n += inv_param->active[offset+i];
    }
    pull_weak_earth = (0 < n);
  }

  /* pull earth's max weakening in the interior of weak zones:
   *   weak_factor_interior_p = sqrt (-log (weak_factor_interior))
   */
  if (pull_weak_earth) {
    const double       *weak_earth = weak_options->weak_factor_interior_earth;
    int                 offset_weak;

    RHEA_ASSERT (weak_options->weak_factor_interior_earth != NULL);

    offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB;
    offset_weak = 0;
    for (i = 0; i < RHEA_WEAKZONE_LABEL_EARTH_N_SL; i++) {
      rhea_inversion_param_calculate_inversion_val_if_active_exp2 (
          weak_earth[offset_weak+i], offset+i, parameter_vec, inv_param);
    }

    offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_RIDGE;
    offset_weak = RHEA_WEAKZONE_LABEL_EARTH_N_SL;
    for (i = 0; i < RHEA_WEAKZONE_LABEL_EARTH_N_RI; i++) {
      rhea_inversion_param_calculate_inversion_val_if_active_exp2 (
          weak_earth[offset_weak+i], offset+i, parameter_vec, inv_param);
    }

    offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_FRACTURE;
    offset_weak = RHEA_WEAKZONE_LABEL_EARTH_N_SL +
                  RHEA_WEAKZONE_LABEL_EARTH_N_RI;
    for (i = 0; i < RHEA_WEAKZONE_LABEL_EARTH_N_FZ; i++) {
      rhea_inversion_param_calculate_inversion_val_if_active_exp2 (
          weak_earth[offset_weak+i], offset+i, parameter_vec, inv_param);
    }
  }

  /* check output */
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (parameter_vec, inv_param));
}

static int
rhea_inversion_param_calculate_model_val_if_active_exp (
                                            double *model_val,
                                            const int param_idx,
                                            ymir_vec_t *parameter_vec,
                                            rhea_inversion_param_t *inv_param)
{
  double             *param = parameter_vec->meshfree->e[0];

  if (inv_param->active[param_idx]) {
    *model_val = exp (param[param_idx]);
    RHEA_ASSERT (isfinite (*model_val));
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
                                            ymir_vec_t *parameter_vec,
                                            rhea_inversion_param_t *inv_param)
{
  double             *param = parameter_vec->meshfree->e[0];

  if (inv_param->active[param_idx]) {
    *model_val = exp (-param[param_idx]*param[param_idx]);
    RHEA_ASSERT (isfinite (*model_val));
    return 1;
  }
  else {
    return 0;
  }
}

void
rhea_inversion_param_push_to_model (ymir_vec_t *parameter_vec,
                                    rhea_inversion_param_t *inv_param)
{
  rhea_weakzone_options_t  *weak_options = inv_param->weak_options;
  rhea_viscosity_options_t *visc_options = inv_param->visc_options;
  int                 pull_weak_earth;
  int                 offset, i;

  RHEA_GLOBAL_VERBOSEF_FN_TAG (__func__, "n_active=%i", inv_param->n_active);

  /*
   * Push to Viscosity Options
   */

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (parameter_vec, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (parameter_vec, inv_param));

  /* push lower and upper bounds of the viscosity:
   *   bound = exp (bound_p)
   */
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(visc_options->min), RHEA_INVERSION_PARAM_VISC_MIN,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(visc_options->max), RHEA_INVERSION_PARAM_VISC_MAX,
      parameter_vec, inv_param);

  /* push scaling factors:
   *   scaling = exp (scaling_p)
   */
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(visc_options->upper_mantle_scaling),
      RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(visc_options->lower_mantle_scaling),
      RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING,
      parameter_vec, inv_param);

  /* push activation energy in Arrhenius relationship:
   *   arrhenius_activation_energy = exp (arrhenius_activation_energy_p)
   */
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(visc_options->upper_mantle_arrhenius_activation_energy),
      RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(visc_options->lower_mantle_arrhenius_activation_energy),
      RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY,
      parameter_vec, inv_param);

  /* push stress exponent that governs strain rate weakening:
   *   stress_exponent = 1 + exp (stress_exponent_p)
   */
  if (inv_param->active[RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT]) {
    double             *param = parameter_vec->meshfree->e[0];

    visc_options->stress_exponent =
      1.0 + exp (param[RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT]);
    RHEA_ASSERT (isfinite (visc_options->stress_exponent));
  }

  /* push yield strength:
   *   yield_strength = exp (yield_strength_p)
   */
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(visc_options->yield_strength),
      RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH, parameter_vec, inv_param);

  /*
   * Push to Weak Zone Options
   */

  /* push weak zone geometry:
   *   thickness = exp (thickness_p)
   */
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(weak_options->thickness),
      RHEA_INVERSION_PARAM_WEAK_THICKNESS,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(weak_options->thickness_generic_slab),
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_SLAB,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(weak_options->thickness_generic_ridge),
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_RIDGE,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(weak_options->thickness_generic_fracture),
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_FRACTURE,
      parameter_vec, inv_param);

  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(weak_options->thickness_const),
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(weak_options->thickness_const_generic_slab),
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_SLAB,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(weak_options->thickness_const_generic_ridge),
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_RIDGE,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp (
      &(weak_options->thickness_const_generic_fracture),
      RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_FRACTURE,
      parameter_vec, inv_param);

  /* push generic max weakening in the interior of weak zones:
   *   weak_factor_interior = exp (-weak_factor_interior_p^2)
   */
  rhea_inversion_param_calculate_model_val_if_active_exp2 (
      &(weak_options->weak_factor_interior),
      RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp2 (
      &(weak_options->weak_factor_interior_generic_slab),
      RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_SLAB,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp2 (
      &(weak_options->weak_factor_interior_generic_ridge),
      RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_RIDGE,
      parameter_vec, inv_param);
  rhea_inversion_param_calculate_model_val_if_active_exp2 (
      &(weak_options->weak_factor_interior_generic_fracture),
      RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_FRACTURE,
      parameter_vec, inv_param);

  /* check if earth's weak factors are active */
  {
    int                 n = 0;

    offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB;
    for (i = 0; i < RHEA_WEAKZONE_LABEL_EARTH_N; i++) {
      n += inv_param->active[offset+i];
    }
    pull_weak_earth = (0 < n);
  }

  /* push earth's max weakening in the interior of weak zones:
   *   weak_factor_interior = exp (-weak_factor_interior_p^2)
   */
  if (pull_weak_earth) {
    int                 offset_weak;

    RHEA_ASSERT (weak_options->weak_factor_interior_earth != NULL);

    offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB;
    offset_weak = 0;
    for (i = 0; i < RHEA_WEAKZONE_LABEL_EARTH_N_SL; i++) {
      rhea_inversion_param_calculate_model_val_if_active_exp2 (
          &(weak_options->weak_factor_interior_earth[offset_weak+i]),
          offset+i, parameter_vec, inv_param);
    }

    offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_RIDGE;
    offset_weak = RHEA_WEAKZONE_LABEL_EARTH_N_SL;
    for (i = 0; i < RHEA_WEAKZONE_LABEL_EARTH_N_RI; i++) {
      rhea_inversion_param_calculate_model_val_if_active_exp2 (
          &(weak_options->weak_factor_interior_earth[offset_weak+i]),
          offset+i, parameter_vec, inv_param);
    }

    offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_FRACTURE;
    offset_weak = RHEA_WEAKZONE_LABEL_EARTH_N_SL +
                  RHEA_WEAKZONE_LABEL_EARTH_N_RI;
    for (i = 0; i < RHEA_WEAKZONE_LABEL_EARTH_N_FZ; i++) {
      rhea_inversion_param_calculate_model_val_if_active_exp2 (
          &(weak_options->weak_factor_interior_earth[offset_weak+i]),
          offset+i, parameter_vec, inv_param);
    }
  }
}

static rhea_viscosity_param_derivative_t
rhea_inversion_param_get_derivative_type (
                                    const rhea_inversion_param_idx_t param_idx)
{
  switch (param_idx) {
  case RHEA_INVERSION_PARAM_VISC_MIN:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_MIN;
  case RHEA_INVERSION_PARAM_VISC_MAX:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_MAX;
  case RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_UPPER_MANTLE_SCALING;
  case RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_LOWER_MANTLE_SCALING;
  case RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_UPPER_MANTLE_ACTIVATION_ENERGY;
  case RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_LOWER_MANTLE_ACTIVATION_ENERGY;
  case RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_STRESS_EXPONENT;
  case RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_YIELD_STRENGTH;

  case RHEA_INVERSION_PARAM_WEAK_THICKNESS:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_WEAK_THICKNESS;
//case RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_SLAB:
//case RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_RIDGE:
//case RHEA_INVERSION_PARAM_WEAK_THICKNESS_GENERIC_FRACTURE:

  case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_WEAK_THICKNESS_CONST;
//case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_SLAB:
//case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_RIDGE:
//case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_GENERIC_FRACTURE:

  case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_WEAK_FACTOR_INTERIOR;
//case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_SLAB:
//case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_RIDGE:
//case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_GENERIC_FRACTURE:

//case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB:
//case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_RIDGE:
//case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_FRACTURE:
  default: /* unknown derivative type */
    RHEA_ABORT_NOT_REACHED ();
  }
}

void
rhea_inversion_param_compute_gradient (ymir_vec_t *gradient,
                                       ymir_vec_t *forward_vel_press,
                                       ymir_vec_t *adjoint_vel_press,
                                       rhea_inversion_param_t *inv_param)
{
  rhea_stokes_problem_t    *stokes_problem = inv_param->stokes_problem;
  rhea_domain_options_t    *domain_options =
    rhea_stokes_problem_get_domain_options (stokes_problem);
  ymir_mesh_t          *ymir_mesh =
                          rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  ymir_pressure_elem_t *press_elem =
                          rhea_stokes_problem_get_press_elem (stokes_problem);
  ymir_vec_t         *temperature =
                        rhea_stokes_problem_get_temperature (stokes_problem);
  ymir_vec_t         *weakzone =
                        rhea_stokes_problem_get_weakzone (stokes_problem);
  ymir_vec_t         *forward_vel, *adjoint_vel, *op_out_vel;
  ymir_vec_t         *viscosity, *bounds_marker, *yielding_marker;
  ymir_vel_dir_t     *vel_dir;
  ymir_stress_op_t   *stress_op;

  ymir_vec_t         *derivative;
  const int           n_parameters = inv_param->n_parameters;
  const int          *active = inv_param->active;
  double             *grad = gradient->meshfree->e[0];
  int                 i;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (gradient, inv_param));
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (forward_vel_press));
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (adjoint_vel_press));

  /* retrieve forward and adjoint velocities */
  forward_vel = rhea_velocity_new (ymir_mesh);
  adjoint_vel = rhea_velocity_new (ymir_mesh);
  op_out_vel = rhea_velocity_new (ymir_mesh);
  rhea_velocity_pressure_copy_components (forward_vel, NULL, forward_vel_press,
                                          press_elem);
  rhea_velocity_pressure_copy_components (adjoint_vel, NULL, adjoint_vel_press,
                                          press_elem);

  /* compute viscosity and related fields */
  viscosity = rhea_viscosity_new (ymir_mesh);
  bounds_marker = rhea_viscosity_new (ymir_mesh);
  yielding_marker = rhea_viscosity_new (ymir_mesh);
  rhea_viscosity_compute (
      /* out: */ viscosity, NULL, bounds_marker, yielding_marker,
      /* in:  */ temperature, weakzone, forward_vel, inv_param->visc_options);

  /* init derivative to use as viscous stress coefficient */
  derivative = rhea_viscosity_new (ymir_mesh);
  ymir_vec_set_value (derivative, 1.0);

  /* create stress operator */
  vel_dir = rhea_domain_create_velocity_dirichlet_bc (
      ymir_mesh, NULL /* dirscal */, domain_options);
  stress_op = ymir_stress_op_new_ext (
      derivative, vel_dir,
      NULL /* Robin BC's */,
      NULL /* deprecated */,
      NULL /* deprecated */,
      domain_options->center, domain_options->moment_of_inertia);

  /* compute gradient entry for each parameter */
  ymir_vec_set_zero (gradient);
  for (i = 0; i < n_parameters; i++) { /* loop over all (possible) parameters */
    if (active[i]) {
      /* compute parameter derivative of viscosity */
      rhea_viscosity_param_derivative (
          derivative, rhea_inversion_param_get_derivative_type (i),
          viscosity, bounds_marker, yielding_marker, temperature, weakzone,
          forward_vel, inv_param->visc_options);

      /* transform to viscous stress coefficient */
      ymir_vec_scale (2.0, derivative);

      /* set derivative as the coefficient of the viscous stress operator */
      ymir_stress_op_set_coeff_scal (stress_op, derivative);

      /* apply viscous stress operator to forward velocity */
      ymir_stress_pc_apply_stress_op (forward_vel, op_out_vel, stress_op,
                                      0 /* !linearized */, 1 /* dirty */);

      /* compute inner product with adjoint velocity */
      grad[i] = ymir_vec_innerprod (op_out_vel, adjoint_vel);
    }
  }
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (gradient, inv_param));

  /* compute & add prior term */
  //TODO

  /* destroy */
  ymir_stress_op_destroy (stress_op);
  ymir_vel_dir_destroy (vel_dir);
  rhea_viscosity_destroy (derivative);
  rhea_viscosity_destroy (viscosity);
  rhea_viscosity_destroy (bounds_marker);
  rhea_viscosity_destroy (yielding_marker);
  rhea_velocity_destroy (forward_vel);
  rhea_velocity_destroy (adjoint_vel);
  rhea_velocity_destroy (op_out_vel);
}

double
rhea_inversion_param_compute_gradient_norm (ymir_vec_t *gradient,
                                            rhea_inversion_param_t *inv_param)
{
  const int           n_parameters = inv_param->n_parameters;
  const int          *active = inv_param->active;
  const double       *grad = gradient->meshfree->e[0];
  double              sum_of_squares = 0.0;
  int                 i;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (gradient, inv_param));

  /* sum up (squares of) entries of the gradient vector that are active */
  for (i = 0; i < n_parameters; i++) { /* loop over all (possible) parameters */
    if (active[i]) {
      sum_of_squares += grad[i] * grad[i];
    }
  }

  /* return l2-norm of active entries of the gradient vector */
  return sqrt (sum_of_squares);
}

void
rhea_inversion_param_incremental_forward_rhs (ymir_vec_t *rhs_vel_mass,
                                              ymir_vec_t *gradient_direction,
                                              ymir_vec_t *forward_vel_press,
                                              rhea_inversion_param_t *inv_param)
{
  rhea_stokes_problem_t    *stokes_problem = inv_param->stokes_problem;
  rhea_domain_options_t    *domain_options =
    rhea_stokes_problem_get_domain_options (stokes_problem);
  ymir_mesh_t          *ymir_mesh =
                          rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  ymir_pressure_elem_t *press_elem =
                          rhea_stokes_problem_get_press_elem (stokes_problem);
  ymir_vec_t         *temperature =
                        rhea_stokes_problem_get_temperature (stokes_problem);
  ymir_vec_t         *weakzone =
                        rhea_stokes_problem_get_weakzone (stokes_problem);
  ymir_vec_t         *forward_vel, *op_out_vel;
  ymir_vec_t         *viscosity, *bounds_marker, *yielding_marker;
  ymir_vel_dir_t     *vel_dir;
  ymir_stress_op_t   *stress_op;

  ymir_vec_t         *derivative;
  const int           n_parameters = inv_param->n_parameters;
  const int          *active = inv_param->active;
  double             *grad_dir = gradient_direction->meshfree->e[0];
  int                 i;

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (rhs_vel_mass));
  RHEA_ASSERT (
      rhea_inversion_param_vec_check_type (gradient_direction, inv_param));
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (forward_vel_press));

  /* retrieve forward velocity */
  forward_vel = rhea_velocity_new (ymir_mesh);
  op_out_vel = rhea_velocity_new (ymir_mesh);
  rhea_velocity_pressure_copy_components (forward_vel, NULL, forward_vel_press,
                                          press_elem);

  /* compute viscosity and related fields */
  viscosity = rhea_viscosity_new (ymir_mesh);
  bounds_marker = rhea_viscosity_new (ymir_mesh);
  yielding_marker = rhea_viscosity_new (ymir_mesh);
  rhea_viscosity_compute (
      /* out: */ viscosity, NULL, bounds_marker, yielding_marker,
      /* in:  */ temperature, weakzone, forward_vel, inv_param->visc_options);

  /* init derivative to use as viscous stress coefficient */
  derivative = rhea_viscosity_new (ymir_mesh);
  ymir_vec_set_value (derivative, 1.0);

  /* create stress operator */
  vel_dir = rhea_domain_create_velocity_dirichlet_bc (
      ymir_mesh, NULL /* dirscal */, domain_options);
  stress_op = ymir_stress_op_new_ext (
      derivative, vel_dir,
      NULL /* Robin BC's */,
      NULL /* deprecated */,
      NULL /* deprecated */,
      domain_options->center, domain_options->moment_of_inertia);

  /* add gradient w.r.t. each parameter */
  ymir_vec_set_zero (rhs_vel_mass);
  for (i = 0; i < n_parameters; i++) { /* loop over all (possible) parameters */
    if (active[i] && DBL_MIN < fabs (grad_dir[i])) {
      /* compute parameter derivative of viscosity */
      rhea_viscosity_param_derivative (
          derivative, rhea_inversion_param_get_derivative_type (i),
          viscosity, bounds_marker, yielding_marker, temperature, weakzone,
          forward_vel, inv_param->visc_options);

      /* transform to viscous stress coefficient */
      ymir_vec_scale (2.0, derivative);

      /* set derivative as the coefficient of the viscous stress operator */
      ymir_stress_op_set_coeff_scal (stress_op, derivative);

      /* apply viscous stress operator to forward velocity */
      ymir_stress_pc_apply_stress_op (forward_vel, op_out_vel, stress_op,
                                      0 /* !linearized */, 1 /* dirty */);

      /* scale and add gradient to right-hand side */
      ymir_vec_add (-grad_dir[i], op_out_vel, rhs_vel_mass);
    }
  }
  RHEA_ASSERT (rhea_velocity_is_valid (rhs_vel_mass));

  /* destroy */
  ymir_stress_op_destroy (stress_op);
  ymir_vel_dir_destroy (vel_dir);
  rhea_viscosity_destroy (derivative);
  rhea_viscosity_destroy (viscosity);
  rhea_viscosity_destroy (bounds_marker);
  rhea_viscosity_destroy (yielding_marker);
  rhea_velocity_destroy (forward_vel);
  rhea_velocity_destroy (op_out_vel);
}

void
rhea_inversion_param_apply_hessian (ymir_vec_t *param_vec_out,
                                    ymir_vec_t *param_vec_in,
                                    ymir_vec_t *forward_vel_press,
                                    ymir_vec_t *adjoint_vel_press,
                                    ymir_vec_t *incr_forward_vel_press,
                                    ymir_vec_t *incr_adjoint_vel_press,
                                    rhea_inversion_param_t *inv_param)
{
  const int           first_order_approx = (adjoint_vel_press == NULL);

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (param_vec_out, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (param_vec_in, inv_param));
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (forward_vel_press));
  RHEA_ASSERT (
      rhea_velocity_pressure_check_vec_type (incr_adjoint_vel_press));

  /*
   * First-Order Derivative Terms (aka. Gauss-Newton Hessian)
   */

  /* init output with parameter gradient function using incr. adj. vel. */
  rhea_inversion_param_compute_gradient (param_vec_out, forward_vel_press,
                                         incr_adjoint_vel_press, inv_param);

  /* compute & add prior term */
  //TODO

  /*
   * Second-Order Derivative Terms
   */

  if (!first_order_approx) { /* if full Hessian */
    RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (adjoint_vel_press));
    RHEA_ASSERT (
        rhea_velocity_pressure_check_vec_type (incr_forward_vel_press));

    //TODO
    RHEA_ABORT_NOT_REACHED ();
  }
}

void
rhea_inversion_param_print (ymir_vec_t *parameter_vec,
                            rhea_inversion_param_t *inv_param)
{
  const double       *param = parameter_vec->meshfree->e[0];
  const int          *active = inv_param->active;
  int                 i;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (parameter_vec, inv_param));

  /* exit if nothing to do */
  if (!inv_param->n_active) {
    return;
  }

  /* print each parameter */
  for (i = 0; i < inv_param->n_parameters; i++) {
    if (active[i]) {
      RHEA_GLOBAL_INFOF ("param# %3i: %g\n", i, param[i]);
    }
  }
}

/******************************************************************************
 * Parameter Vector
 *****************************************************************************/

ymir_vec_t *
rhea_inversion_param_vec_new (rhea_inversion_param_t *inv_param)
{
  return ymir_vec_new_meshfree (inv_param->n_parameters);
}

void
rhea_inversion_param_vec_destroy (ymir_vec_t *vec)
{
  ymir_vec_destroy (vec);
}

int
rhea_inversion_param_vec_check_type (ymir_vec_t *vec,
                                     rhea_inversion_param_t *inv_param)
{
  return (
      ymir_vec_is_meshfree (vec) &&
      vec->n_meshfree == inv_param->n_parameters &&
      RHEA_INVERSION_PARAM_N == inv_param->n_parameters
  );
}

int
rhea_inversion_param_vec_is_valid (ymir_vec_t *vec,
                                   rhea_inversion_param_t *inv_param)
{
  const int           n_parameters = inv_param->n_parameters;
  const int          *active = inv_param->active;
  const double       *v = vec->meshfree->e[0];
  int                 i;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (vec, inv_param));

  for (i = 0; i < n_parameters; i++) { /* loop over all (possible) parameters */
    if (active[i] && !isfinite (v[i])) {
      return 0;
    }
  }
  return 1;
}

sc_dmatrix_t *
rhea_inversion_param_vec_reduced_new (ymir_vec_t *vec,
                                      rhea_inversion_param_t *inv_param)
{
  const int           n_parameters = inv_param->n_parameters;
  const int           n_active = inv_param->n_active;
  const int          *active = inv_param->active;
  sc_dmatrix_t       *vec_reduced;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (vec, inv_param));

  /* create reduced parameter vector */
  vec_reduced = sc_dmatrix_new (n_active, 1);

  /* fill entries */
  {
    const double       *v = vec->meshfree->e[0];
    double             *r = vec_reduced->e[0];
    int                 i, row;

    row = 0;
    for (i = 0; i < n_parameters; i++) { /* loop over all parameters */
      if (active[i]) {
        r[row] = v[i];
        row++;
      }
    }
  }

  return vec_reduced;
}

void
rhea_inversion_param_vec_reduced_destroy (sc_dmatrix_t *vec_reduced)
{
  sc_dmatrix_destroy (vec_reduced);
}

void
rhea_inversion_param_vec_reduced_copy (ymir_vec_t *vec,
                                       sc_dmatrix_t *vec_reduced,
                                       rhea_inversion_param_t *inv_param)
{
  const int           n_parameters = inv_param->n_parameters;
  const int           n_active = inv_param->n_active;
  const int          *active = inv_param->active;
  const double       *r = vec_reduced->e[0];
  double             *v = vec->meshfree->e[0];
  int                 i, row;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (vec, inv_param));
  RHEA_ASSERT (vec_reduced->m == n_active);
  RHEA_ASSERT (vec_reduced->n == 1);

  /* fill entries */
  row = 0;
  for (i = 0; i < n_parameters; i++) { /* loop over all parameters */
    if (active[i]) {
      v[i] = r[row];
      row++;
    }
  }
}

/******************************************************************************
 * Data Access
 *****************************************************************************/

int
rhea_inversion_param_get_n_parameters (rhea_inversion_param_t *inv_param)
{
  return inv_param->n_parameters;
}

int
rhea_inversion_param_get_n_active (rhea_inversion_param_t *inv_param)
{
  return inv_param->n_active;
}

int *
rhea_inversion_param_get_active (rhea_inversion_param_t *inv_param)
{
  return inv_param->active;
}
