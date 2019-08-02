#include <rhea_inversion_param.h>
#include <rhea_base.h>
#include <rhea_viscosity_param_derivative.h>
#include <rhea_weakzone_label.h>
#include <ymir_stress_pc.h>
#include <fenv.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_INVERSION_PARAM_DEFAULT_ACTIVATE 1
#define RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE 0
#define RHEA_INVERSION_PARAM_DEFAULT_PRMN_PERTURB_STDDEV NAN
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_MIN 2.3
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_MAX 2.3
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_UM_SCALING 4.6
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_UM_ARRHENIUS_ACTIVATION_ENERGY 2.3
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_LM_SCALING 4.6
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_LM_ARRHENIUS_ACTIVATION_ENERGY 2.3
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_STRESS_EXPONENT 2.3
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_YIELD_STRENGTH 4.6
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_THICKNESS 2.3
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_THICKNESS_CONST 2.3
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_WEAK_FACTOR_INTERIOR 4.6
#define RHEA_INVERSION_PARAM_DEFAULT_GUESS_PERTURB_STDDEV 0.0
#define RHEA_INVERSION_PARAM_DEFAULT_GUESS_SHIFT_BY_PRIOR_STDDEV 0.0

#define RHEA_INVERSION_PARAM_DEFAULT_WEAK_FACTOR_MIN SC_EPS

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

  /****** Activation Flags ******/

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
  YMIR_OPTIONS_B, "activate-upper-mantle-arrhenius-activation-energy", '\0',
    &(inv_param_opt->upper_mantle_arrhenius_activation_energy_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate the activation energy in Arrhenius relationship of upper mantle",
  YMIR_OPTIONS_B, "activate-lower-mantle-scaling", '\0',
    &(inv_param_opt->lower_mantle_scaling_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate scaling factor of lower mantle",
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
  YMIR_OPTIONS_B, "activate-weakzone-thickness-class-slab", '\0',
    &(inv_param_opt->thickness_class_slab_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone thickness of all slabs",
  YMIR_OPTIONS_B, "activate-weakzone-thickness-class-ridge", '\0',
    &(inv_param_opt->thickness_class_ridge_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone thickness of all ridges",
  YMIR_OPTIONS_B, "activate-weakzone-thickness-class-fracture", '\0',
    &(inv_param_opt->thickness_class_fracture_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone thickness of all fractures",

  YMIR_OPTIONS_B, "activate-weakzone-thickness-const", '\0',
    &(inv_param_opt->thickness_const_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone interior thickness (where min factor is reached)",
  YMIR_OPTIONS_B, "activate-weakzone-thickness-const-class-slab", '\0',
    &(inv_param_opt->thickness_const_class_slab_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate interior thickness for all slabs",
  YMIR_OPTIONS_B, "activate-weakzone-thickness-const-class-ridge", '\0',
    &(inv_param_opt->thickness_const_class_ridge_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate interior thickness for all ridges",
  YMIR_OPTIONS_B, "activate-weakzone-thickness-const-class-fracture", '\0',
    &(inv_param_opt->thickness_const_class_fracture_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate interior thickness for all fractures",

  YMIR_OPTIONS_B, "activate-weak-factor-interior", '\0',
    &(inv_param_opt->weak_factor_interior_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone factor",
  YMIR_OPTIONS_B, "activate-weak-factor-interior-class-slab", '\0',
    &(inv_param_opt->weak_factor_interior_class_slab_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone factor of all slabs",
  YMIR_OPTIONS_B, "activate-weak-factor-interior-class-ridge", '\0',
    &(inv_param_opt->weak_factor_interior_class_ridge_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone factor of all ridges",
  YMIR_OPTIONS_B, "activate-weak-factor-interior-class-fracture", '\0',
    &(inv_param_opt->weak_factor_interior_class_fracture_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone factor of all fractures",

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

  /****** Prior ******/

  YMIR_OPTIONS_D, "prior-mean-perturbation-stddev", '\0',
    &(inv_param_opt->prior_mean_perturb_stddev),
    RHEA_INVERSION_PARAM_DEFAULT_PRMN_PERTURB_STDDEV,
    "Prior mean: standard deviation for perturbations of mean",

  YMIR_OPTIONS_D, "prior-stddev-min", '\0',
    &(inv_param_opt->prior_stddev_min),
    RHEA_INVERSION_PARAM_DEFAULT_PRSD_MIN,
    "Prior std dev: min viscosity bound",
  YMIR_OPTIONS_D, "prior-stddev-max", '\0',
    &(inv_param_opt->prior_stddev_max),
    RHEA_INVERSION_PARAM_DEFAULT_PRSD_MAX,
    "Prior std dev: max viscosity bound",

  YMIR_OPTIONS_D, "prior-stddev-upper-mantle-scaling", '\0',
    &(inv_param_opt->prior_stddev_upper_mantle_scaling),
    RHEA_INVERSION_PARAM_DEFAULT_PRSD_UM_SCALING,
    "Prior std dev: scaling factor of upper mantle",
  YMIR_OPTIONS_D, "prior-stddev-upper-mantle-arrhenius-activation-energy", '\0',
    &(inv_param_opt->prior_stddev_upper_mantle_arrhenius_activation_energy),
    RHEA_INVERSION_PARAM_DEFAULT_PRSD_UM_ARRHENIUS_ACTIVATION_ENERGY,
    "Prior std dev: activation energy in Arrhenius term of upper mantle",
  YMIR_OPTIONS_D, "prior-stddev-lower-mantle-scaling", '\0',
    &(inv_param_opt->prior_stddev_lower_mantle_scaling),
    RHEA_INVERSION_PARAM_DEFAULT_PRSD_LM_SCALING,
    "Prior std dev: scaling factor of lower mantle",
  YMIR_OPTIONS_D, "prior-stddev-lower-mantle-arrhenius-activation-energy", '\0',
    &(inv_param_opt->prior_stddev_lower_mantle_arrhenius_activation_energy),
    RHEA_INVERSION_PARAM_DEFAULT_PRSD_LM_ARRHENIUS_ACTIVATION_ENERGY,
    "Prior std dev: activation energy in Arrhenius term of lower mantle",

  YMIR_OPTIONS_D, "prior-stddev-stress-exponent", '\0',
    &(inv_param_opt->prior_stddev_stress_exponent),
    RHEA_INVERSION_PARAM_DEFAULT_PRSD_STRESS_EXPONENT,
    "Prior std dev: stress exponent that governs strain rate weakening",
  YMIR_OPTIONS_D, "prior-stddev-yield-strength", '\0',
    &(inv_param_opt->prior_stddev_yield_strength),
    RHEA_INVERSION_PARAM_DEFAULT_PRSD_YIELD_STRENGTH,
    "Prior std dev: plastic yielding",

  YMIR_OPTIONS_D, "prior-stddev-weakzone-thickness", '\0',
    &(inv_param_opt->prior_stddev_thickness),
    RHEA_INVERSION_PARAM_DEFAULT_PRSD_THICKNESS,
    "Prior std dev: weak zone thickness",
  YMIR_OPTIONS_D, "prior-stddev-weakzone-thickness-const", '\0',
    &(inv_param_opt->prior_stddev_thickness_const),
    RHEA_INVERSION_PARAM_DEFAULT_PRSD_THICKNESS_CONST,
    "Prior std dev: weak zone interior thickness (where min factor is reached)",
  YMIR_OPTIONS_D, "prior-stddev-weak-factor-interior", '\0',
    &(inv_param_opt->prior_stddev_weak_factor_interior),
    RHEA_INVERSION_PARAM_DEFAULT_PRSD_WEAK_FACTOR_INTERIOR,
    "Prior std dev: weak zone factor",

  /****** Initial Guess ******/

  YMIR_OPTIONS_D, "initial-guess-perturbation-stddev", '\0',
    &(inv_param_opt->initial_guess_perturb_stddev),
    RHEA_INVERSION_PARAM_DEFAULT_GUESS_PERTURB_STDDEV,
    "Initial guess: standard deviation for perturbations",

  YMIR_OPTIONS_D, "initial-guess-shift-by-prior-stddev", '\0',
    &(inv_param_opt->initial_guess_shift_by_prior_stddev),
    RHEA_INVERSION_PARAM_DEFAULT_GUESS_SHIFT_BY_PRIOR_STDDEV,
    "Initial guess: scaling for elementwise shift by std. dev. of prior",

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

  /* parameterization of scaling factors and of the activation energy in the
   * Arrhenius relationship */
  RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING,
  RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY,
  RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING,
  RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY,

  /* parameterization of stress exponent that governs strain rate weakening */
  RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT,

  /* parameterization of plastic yielding */
  RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH,

  /*
   * Parameterization of Weak Zone Options
   */

  /* parameterization of weak zone geometry */
  RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_NONE,
  RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_SLAB,
  RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_RIDGE,
  RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_FRACTURE,
  RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_NONE,
  RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_SLAB,
  RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_RIDGE,
  RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_FRACTURE,

  /* parameterization of max weakening factor (in the interior of weak zones) */
  RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_NONE,
  RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_SLAB,
  RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_RIDGE,
  RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE,

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
                                        rhea_inversion_param_options_t *opt,
                                        rhea_weakzone_options_t *weak_options,
                                        rhea_viscosity_options_t *visc_options)
{
  int                *active = RHEA_ALLOC_ZERO (int, RHEA_INVERSION_PARAM_N);
  int                 offset, idx;

  /* set activation mask of viscosity parameters */
  idx = RHEA_INVERSION_PARAM_VISC_MIN;
  active[idx] = (rhea_viscosity_restrict_min (visc_options) && opt->min_a);
  idx = RHEA_INVERSION_PARAM_VISC_MAX;
  active[idx] = (rhea_viscosity_restrict_max (visc_options) && opt->max_a);

  idx = RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING;
  active[idx] = opt->upper_mantle_scaling_a;
  idx = RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY;
  active[idx] = (rhea_viscosity_has_arrhenius (visc_options) &&
                 opt->upper_mantle_arrhenius_activation_energy_a);

  idx = RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING;
  active[idx] = opt->lower_mantle_scaling_a;
  idx = RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY;
  active[idx] = (rhea_viscosity_has_arrhenius (visc_options) &&
                 opt->lower_mantle_arrhenius_activation_energy_a);

  idx = RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT;
  active[idx] = (rhea_viscosity_has_strain_rate_weakening (visc_options) &&
                 opt->stress_exponent_a);
  idx = RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH;
  active[idx] = (rhea_viscosity_has_yielding (visc_options) &&
                 opt->yield_strength_a);

  /* set activation mask of weak zone parameters */
  if (rhea_weakzone_exists (weak_options)) {
    /* process weak zone thickness */
    if (!opt->thickness_class_slab_a &&
        !opt->thickness_class_ridge_a &&
        !opt->thickness_class_fracture_a) {
      idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_NONE;
      active[idx] = opt->thickness_a;
    }
    idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_SLAB;
    active[idx] = opt->thickness_class_slab_a;
    idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_RIDGE;
    active[idx] = opt->thickness_class_ridge_a;
    idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_FRACTURE;
    active[idx] = opt->thickness_class_fracture_a;

    /* process weak zone thickness of the inerior */
    if (!opt->thickness_const_class_slab_a &&
        !opt->thickness_const_class_ridge_a &&
        !opt->thickness_const_class_fracture_a) {
      idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_NONE;
      active[idx] = opt->thickness_const_a;
    }
    idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_SLAB;
    active[idx] = opt->thickness_const_class_slab_a;
    idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_RIDGE;
    active[idx] = opt->thickness_const_class_ridge_a;
    idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_FRACTURE;
    active[idx] = opt->thickness_const_class_fracture_a;

    /* process weak zone factors */
    if (!opt->weak_factor_interior_class_slab_a &&
        !opt->weak_factor_interior_class_ridge_a &&
        !opt->weak_factor_interior_class_fracture_a &&
        !opt->weak_factor_interior_earth_slab_a &&
        !opt->weak_factor_interior_earth_ridge_a &&
        !opt->weak_factor_interior_earth_fracture_a) {
      idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_NONE;
      active[idx] = opt->weak_factor_interior_a;
    }
    if (!opt->weak_factor_interior_earth_slab_a) {
      idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_SLAB;
      active[idx] = opt->weak_factor_interior_class_slab_a;
    }
    else {
      offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB;
      for (idx = offset; idx < offset + RHEA_WEAKZONE_LABEL_EARTH_N_SL; idx++) {
        active[idx] = 1;
      }
    }
    if (!opt->weak_factor_interior_earth_ridge_a) {
      idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_RIDGE;
      active[idx] = opt->weak_factor_interior_class_ridge_a;
    }
    else {
      offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_RIDGE;
      for (idx = offset; idx < offset + RHEA_WEAKZONE_LABEL_EARTH_N_RI; idx++) {
        active[idx] = 1;
      }
    }
    if (!opt->weak_factor_interior_earth_fracture_a) {
      idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE;
      active[idx] = opt->weak_factor_interior_class_fracture_a;
    }
    else {
      offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_FRACTURE;
      for (idx = offset; idx < offset + RHEA_WEAKZONE_LABEL_EARTH_N_FZ; idx++) {
        active[idx] = 1;
      }
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

/* required function declarations */
static void         rhea_inversion_param_prior_get_mean (
                                          ymir_vec_t *prior_mean_vec,
                                          rhea_inversion_param_t *inv_param,
                                          rhea_inversion_param_options_t *opt);
static void         rhea_inversion_param_prior_get_stddev (
                                          ymir_vec_t *prior_stddev_vec,
                                          rhea_inversion_param_t *inv_param,
                                          rhea_inversion_param_options_t *opt);

/* inversion parameters */
struct rhea_inversion_param
{
  /* activation mask for inversion parameters */
  int                *active;
  int                 n_active;
  int                 n_parameters;

  /* prior data */
  ymir_vec_t         *prior_mean;
  ymir_vec_t         *prior_stddev;

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

  /* create prior data */
  inv_param->prior_mean = rhea_inversion_param_vec_new (inv_param);
  inv_param->prior_stddev = rhea_inversion_param_vec_new (inv_param);

  /* set mean and standard deviation of the Gaussian prior */
  rhea_inversion_param_prior_get_mean (inv_param->prior_mean, inv_param,
                                       inv_param_opt);
  rhea_inversion_param_prior_get_stddev (inv_param->prior_stddev, inv_param,
                                         inv_param_opt);

  /* return inversion parameters */
  return inv_param;
}

void
rhea_inversion_param_destroy (rhea_inversion_param_t *inv_param)
{
  rhea_inversion_param_vec_destroy (inv_param->prior_mean);
  rhea_inversion_param_vec_destroy (inv_param->prior_stddev);
  rhea_inversion_param_activation_mask_destroy (inv_param->active);
  RHEA_FREE (inv_param);
}

/**
 * Checks if earth's weak factors are active parameters.
 */
static int
rhea_inversion_param_earth_weak_is_active (rhea_inversion_param_t *inv_param)
{
  int                 offset, idx;
  int                 n = 0;

  offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB;
  for (idx = offset; idx < offset + RHEA_WEAKZONE_LABEL_EARTH_N; idx++) {
    n += inv_param->active[idx];
  }

  return n;
}

/******************************************************************************
 * Mapping between Parameter Values and Model Values
 *****************************************************************************/

static void
rhea_inversion_param_get_model_vals (ymir_vec_t *parameter_vec,
                                     rhea_inversion_param_t *inv_param)
{
  const int          *active = inv_param->active;
  double             *p = parameter_vec->meshfree->e[0];
  rhea_weakzone_options_t  *weak_options = inv_param->weak_options;
  rhea_viscosity_options_t *visc_options = inv_param->visc_options;
  int                 idx;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (parameter_vec, inv_param));

  /* initialize to zero */
  ymir_vec_set_zero (parameter_vec);

  /*
   * Viscosity Options
   */

  /* get lower and upper bounds of the viscosity */
  idx = RHEA_INVERSION_PARAM_VISC_MIN;
  if (active[idx]) p[idx] = visc_options->min;
  idx = RHEA_INVERSION_PARAM_VISC_MAX;
  if (active[idx]) p[idx] = visc_options->max;

  /* get scaling factors and activation energy in Arrhenius relationship */
  idx = RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING;
  if (active[idx]) p[idx] = visc_options->upper_mantle_scaling;
  idx = RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY;
  if (active[idx]) {
    p[idx] = visc_options->upper_mantle_arrhenius_activation_energy;
  }
  idx = RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING;
  if (active[idx]) p[idx] = visc_options->lower_mantle_scaling;
  idx = RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY;
  if (active[idx]) {
    p[idx] = visc_options->lower_mantle_arrhenius_activation_energy;
  }

  /* get stress exponent that governs strain rate weakening */
  idx = RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT;
  if (active[idx]) p[idx] = visc_options->stress_exponent;

  /* get yield strength */
  idx = RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH;
  if (active[idx]) p[idx] = visc_options->yield_strength;

  /*
   * Weak Zone Options
   */

  /* get weak zone geometry */
  idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_NONE;
  if (active[idx]) p[idx] = weak_options->thickness;
  idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_SLAB;
  if (active[idx]) p[idx] = weak_options->thickness_class_slab;
  idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_RIDGE;
  if (active[idx]) p[idx] = weak_options->thickness_class_ridge;
  idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_FRACTURE;
  if (active[idx]) p[idx] = weak_options->thickness_class_fracture;

  idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_NONE;
  if (active[idx]) p[idx] = weak_options->thickness_const;
  idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_SLAB;
  if (active[idx]) p[idx] = weak_options->thickness_const_class_slab;
  idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_RIDGE;
  if (active[idx]) p[idx] = weak_options->thickness_const_class_ridge;
  idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_FRACTURE;
  if (active[idx]) p[idx] = weak_options->thickness_const_class_fracture;

  /* get max weakening factor */
  idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_NONE;
  if (active[idx]) p[idx] = weak_options->weak_factor_interior;
  idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_SLAB;
  if (active[idx]) p[idx] = weak_options->weak_factor_interior_class_slab;
  idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_RIDGE;
  if (active[idx]) p[idx] = weak_options->weak_factor_interior_class_ridge;
  idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE;
  if (active[idx]) p[idx] = weak_options->weak_factor_interior_class_fracture;

  /* get earth's max weakening in the interior of weak zones */
  if (rhea_inversion_param_earth_weak_is_active (inv_param)) {
    const double       *weak_earth = weak_options->weak_factor_interior_earth;
    int                 offset;

    RHEA_ASSERT (weak_options->weak_factor_interior_earth != NULL);

    offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB;
    for (idx = 0; idx < RHEA_WEAKZONE_LABEL_EARTH_N; idx++) {
      if (active[offset+idx]) p[offset+idx] = weak_earth[idx];
    }
  }

  /* check output */
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (parameter_vec, inv_param));
}

static void
rhea_inversion_param_set_model_vals (ymir_vec_t *parameter_vec,
                                     rhea_inversion_param_t *inv_param)
{
  const int          *active = inv_param->active;
  const double       *p = parameter_vec->meshfree->e[0];
  rhea_weakzone_options_t  *weak_options = inv_param->weak_options;
  rhea_viscosity_options_t *visc_options = inv_param->visc_options;
  int                 idx;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (parameter_vec, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (parameter_vec, inv_param));

  /*
   * Viscosity Options
   */

  /* set lower and upper bounds of the viscosity */
  idx = RHEA_INVERSION_PARAM_VISC_MIN;
  if (active[idx]) visc_options->min = p[idx];
  idx = RHEA_INVERSION_PARAM_VISC_MAX;
  if (active[idx]) visc_options->max = p[idx];

  /* set scaling factors and activation energy in Arrhenius relationship */
  idx = RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING;
  if (active[idx]) visc_options->upper_mantle_scaling = p[idx];
  idx = RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY;
  if (active[idx]) {
    visc_options->upper_mantle_arrhenius_activation_energy = p[idx];
  }
  idx = RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING;
  if (active[idx]) visc_options->lower_mantle_scaling = p[idx];
  idx = RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY;
  if (active[idx]) {
    visc_options->lower_mantle_arrhenius_activation_energy = p[idx];
  }

  /* set stress exponent that governs strain rate weakening */
  idx = RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT;
  if (active[idx]) visc_options->stress_exponent = p[idx];

  /* set yield strength */
  idx = RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH;
  if (active[idx]) visc_options->yield_strength = p[idx];

  /*
   * Weak Zone Options
   */

  /* set weak zone geometry */
  idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_NONE;
  if (active[idx]) weak_options->thickness = p[idx];
  idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_SLAB;
  if (active[idx]) weak_options->thickness_class_slab = p[idx];
  idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_RIDGE;
  if (active[idx]) weak_options->thickness_class_ridge = p[idx];
  idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_FRACTURE;
  if (active[idx]) weak_options->thickness_class_fracture = p[idx];

  idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_NONE;
  if (active[idx]) weak_options->thickness_const = p[idx];
  idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_SLAB;
  if (active[idx]) weak_options->thickness_const_class_slab = p[idx];
  idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_RIDGE;
  if (active[idx]) weak_options->thickness_const_class_ridge = p[idx];
  idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_FRACTURE;
  if (active[idx]) weak_options->thickness_const_class_fracture = p[idx];

  /* set max weakening factor */
  idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_NONE;
  if (active[idx]) weak_options->weak_factor_interior = p[idx];
  idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_SLAB;
  if (active[idx]) weak_options->weak_factor_interior_class_slab = p[idx];
  idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_RIDGE;
  if (active[idx]) weak_options->weak_factor_interior_class_ridge = p[idx];
  idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE;
  if (active[idx]) weak_options->weak_factor_interior_class_fracture = p[idx];

  /* set earth's max weakening in the interior of weak zones */
  if (rhea_inversion_param_earth_weak_is_active (inv_param)) {
    double             *weak_earth = weak_options->weak_factor_interior_earth;
    int                 offset;

    RHEA_ASSERT (weak_options->weak_factor_interior_earth != NULL);

    offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB;
    for (idx = 0; idx < RHEA_WEAKZONE_LABEL_EARTH_N; idx++) {
      if (active[offset+idx]) weak_earth[idx] = p[offset+idx];
    }
  }
}

#if 0 //TODO remove unused code
double
rhea_inversion_param_convert_to_model_pos (const double inv_param_val)
{
  double              model_val;

  feclearexcept (FE_ALL_EXCEPT);
  model_val = exp (inv_param_val);
  if (fetestexcept (FE_OVERFLOW)) { /* if argument of exp too large */
    model_val = DBL_MAX;
  }

  return model_val;
}

double
rhea_inversion_param_convert_to_model_pos_deriv (const double inv_param_val)
{
  return rhea_inversion_param_convert_to_model_pos (inv_param_val);
}

double
rhea_inversion_param_convert_from_model_pos (const double model_val)
{
  double              inv_param_val;

  RHEA_ASSERT (0.0 <= model_val);

  inv_param_val = log (model_val);
  if (isfinite (inv_param_val)) {
    return inv_param_val;
  }
  else {
    return -DBL_MAX;
  }
}

double
rhea_inversion_param_convert_from_model_pos_deriv (const double model_val)
{
  RHEA_ASSERT (0.0 <= model_val);
  return model_val;
}

double
rhea_inversion_param_convert_to_model_n (const double inv_param_val)
{
  return 1.0 + rhea_inversion_param_convert_to_model_pos (inv_param_val);
}

double
rhea_inversion_param_convert_to_model_n_deriv (const double inv_param_val)
{
  return rhea_inversion_param_convert_to_model_pos (inv_param_val);
}

double
rhea_inversion_param_convert_from_model_n (const double model_val)
{
  RHEA_ASSERT (1.0 <= model_val);
  return rhea_inversion_param_convert_from_model_pos (model_val - 1.0);
}

double
rhea_inversion_param_convert_from_model_n_deriv (const double model_val)
{
  RHEA_ASSERT (1.0 <= model_val);
  return rhea_inversion_param_convert_from_model_pos_deriv (model_val - 1.0);
}

double
rhea_inversion_param_convert_to_model_weak (const double inv_param_val)
{
  const double        weak_min = RHEA_INVERSION_PARAM_DEFAULT_WEAK_FACTOR_MIN;

  return SC_MAX (weak_min, exp (-inv_param_val*inv_param_val));
}

double
rhea_inversion_param_convert_to_model_weak_deriv (const double inv_param_val)
{
  const double        weak_min = RHEA_INVERSION_PARAM_DEFAULT_WEAK_FACTOR_MIN;

  if (weak_min < exp (-inv_param_val*inv_param_val)) {
    return -2.0 * inv_param_val * exp (-inv_param_val*inv_param_val);
  }
  else {
    return 0.0;
  }
}

double
rhea_inversion_param_convert_from_model_weak (const double model_val)
{
  const double        weak_min = RHEA_INVERSION_PARAM_DEFAULT_WEAK_FACTOR_MIN;

  RHEA_ASSERT (0.0 <= model_val && model_val <= 1.0);

  if (weak_min < model_val) {
    return -sqrt (-log (model_val));
  }
  else {
    return -sqrt (-log (weak_min));
  }
}

double
rhea_inversion_param_convert_from_model_weak_deriv (const double model_val)
{
  const double        weak_min = RHEA_INVERSION_PARAM_DEFAULT_WEAK_FACTOR_MIN;

  RHEA_ASSERT (0.0 <= model_val && model_val <= 1.0);

  if (weak_min < model_val) {
    return 2.0 * sqrt (-log (model_val)) * model_val;
  }
  else {
    return 0.0;
  }
}
#endif

double
rhea_inversion_param_convert_to_model_pos (const double inv_param_val)
{
  RHEA_ASSERT (0.0 < inv_param_val);
  return inv_param_val;
}

double
rhea_inversion_param_convert_to_model_pos_deriv (const double inv_param_val)
{
  return 1.0;
}

double
rhea_inversion_param_convert_from_model_pos (const double model_val)
{
  RHEA_ASSERT (0.0 < model_val);
  return model_val;
}

double
rhea_inversion_param_convert_from_model_pos_deriv (const double model_val)
{
  return 1.0;
}

static double
rhea_inversion_param_project_to_feasible_pos (const double inv_param_val)
{
  if (SC_EPS <= inv_param_val) {
    return 0.0;
  }
  else {
    return SC_EPS - inv_param_val;
  }
}

double
rhea_inversion_param_convert_to_model_n (const double inv_param_val)
{
  RHEA_ASSERT (1.0 <= inv_param_val);
  return inv_param_val;
}

double
rhea_inversion_param_convert_to_model_n_deriv (const double inv_param_val)
{
  return 1.0;
}

double
rhea_inversion_param_convert_from_model_n (const double model_val)
{
  RHEA_ASSERT (1.0 <= model_val);
  return model_val;
}

double
rhea_inversion_param_convert_from_model_n_deriv (const double model_val)
{
  return 1.0;
}

static double
rhea_inversion_param_project_to_feasible_n (const double inv_param_val)
{
  if (1.0 <= inv_param_val) {
    return 0.0;
  }
  else {
    return 1.0 - inv_param_val;
  }
}

double
rhea_inversion_param_convert_to_model_weak (const double inv_param_val)
{
  RHEA_ASSERT (0.0 < inv_param_val && inv_param_val <= 1.0);
  return inv_param_val;
}

double
rhea_inversion_param_convert_to_model_weak_deriv (const double inv_param_val)
{
  return 1.0;
}

double
rhea_inversion_param_convert_from_model_weak (const double model_val)
{
  RHEA_ASSERT (0.0 < model_val && model_val <= 1.0);
  return model_val;
}

double
rhea_inversion_param_convert_from_model_weak_deriv (const double model_val)
{
  return 1.0;
}

static double
rhea_inversion_param_project_to_feasible_weak (const double inv_param_val)
{
  if (SC_EPS <= inv_param_val && inv_param_val <= 1.0) {
    return 0.0;
  }
  else if (inv_param_val < SC_EPS) {
    return SC_EPS - inv_param_val;
  }
  else if (1.0 < inv_param_val) {
    return 1.0 - inv_param_val;
  }
  else {
    RHEA_ABORT_NOT_REACHED ();
  }
}

static void
rhea_inversion_param_convert_model_vals_to_params (
                                            ymir_vec_t *parameter_vec,
                                            rhea_inversion_param_t *inv_param)
{
  double             *p = parameter_vec->meshfree->e[0];
  const int          *active = inv_param->active;
  int                 offset, idx;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (parameter_vec, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (parameter_vec, inv_param));

  /* convert the values of active parameters */
  for (idx = 0; idx < inv_param->n_parameters; idx++) {
    int                 success = 0;

    if (active[idx]) {
      /* look up single index */
      switch (idx) {
      case RHEA_INVERSION_PARAM_VISC_MIN:
      case RHEA_INVERSION_PARAM_VISC_MAX:
      case RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING:
      case RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY:
      case RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING:
      case RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY:
      case RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_FRACTURE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_FRACTURE:
        p[idx] = rhea_inversion_param_convert_from_model_pos (p[idx]);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT:
        p[idx] = rhea_inversion_param_convert_from_model_n (p[idx]);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE:
        p[idx] = rhea_inversion_param_convert_from_model_weak (p[idx]);
        success = 1;
        break;
      default:
        break;
      }

      /* look up indices in a range */
      if (!success) {
        offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB;
        if (offset <= idx && idx < offset + RHEA_WEAKZONE_LABEL_EARTH_N) {
          p[idx] = rhea_inversion_param_convert_from_model_weak (p[idx]);
          success = 1;
        }
      }

      if (!success) { /* otherwise parameter index is unknown */
        RHEA_ABORT_NOT_REACHED ();
      }
    }
  }

  /* check output */
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (parameter_vec, inv_param));
}

static void
rhea_inversion_param_convert_params_to_model_vals (
                                            ymir_vec_t *parameter_vec,
                                            rhea_inversion_param_t *inv_param)
{
  double             *p = parameter_vec->meshfree->e[0];
  const int          *active = inv_param->active;
  int                 offset, idx;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (parameter_vec, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (parameter_vec, inv_param));

  /* convert the values of active parameters */
  for (idx = 0; idx < inv_param->n_parameters; idx++) {
    int                 success = 0;

    if (active[idx]) {
      /* look up single index */
      switch (idx) {
      case RHEA_INVERSION_PARAM_VISC_MIN:
      case RHEA_INVERSION_PARAM_VISC_MAX:
      case RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING:
      case RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY:
      case RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING:
      case RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY:
      case RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_FRACTURE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_FRACTURE:
        p[idx] = rhea_inversion_param_convert_to_model_pos (p[idx]);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT:
        p[idx] = rhea_inversion_param_convert_to_model_n (p[idx]);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE:
        p[idx] = rhea_inversion_param_convert_to_model_weak (p[idx]);
        success = 1;
        break;
      default:
        break;
      }

      /* look up indices in a range */
      if (!success) {
        offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB;
        if (offset <= idx && idx < offset + RHEA_WEAKZONE_LABEL_EARTH_N) {
          p[idx] = rhea_inversion_param_convert_to_model_weak (p[idx]);
          success = 1;
        }
      }

      if (!success) { /* otherwise parameter index is unknown */
        RHEA_ABORT_NOT_REACHED ();
      }
    }
  }

  /* check output */
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (parameter_vec, inv_param));
}

double
rhea_inversion_param_reduce_step_length_to_feasible (
                                        ymir_vec_t *step_vec,
                                        ymir_vec_t *parameter_vec,
                                        rhea_inversion_param_t *inv_param)
{
  const double       *param = parameter_vec->meshfree->e[0];
  const double       *step = step_vec->meshfree->e[0];
  double              step_length = 1.0;
  const int          *active = inv_param->active;
  int                 offset, idx;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (step_vec, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (parameter_vec, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (step_vec, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (parameter_vec, inv_param));

  for (idx = 0; idx < inv_param->n_parameters; idx++) {
    const double        p = param[idx];
    const double        s = step[idx];
    double              diff = 0.0;
    int                 success = 0;

    if (active[idx]) {
      /* look up single index */
      switch (idx) {
      case RHEA_INVERSION_PARAM_VISC_MIN:
      case RHEA_INVERSION_PARAM_VISC_MAX:
      case RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING:
      case RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY:
      case RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING:
      case RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY:
      case RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_FRACTURE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_FRACTURE:
        RHEA_ASSERT (
            fabs (rhea_inversion_param_project_to_feasible_pos (p)) <= 0.0);
        diff = rhea_inversion_param_project_to_feasible_pos (p + s);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT:
        RHEA_ASSERT (
            fabs (rhea_inversion_param_project_to_feasible_n (p)) <= 0.0);
        diff = rhea_inversion_param_project_to_feasible_n (p + s);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE:
        RHEA_ASSERT (
            fabs (rhea_inversion_param_project_to_feasible_weak (p)) <= 0.0);
        diff = rhea_inversion_param_project_to_feasible_weak (p + s);
        success = 1;
        break;
      default:
        break;
      }

      /* look up indices in a range */
      if (!success) {
        offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB;
        if (offset <= idx && idx < offset + RHEA_WEAKZONE_LABEL_EARTH_N) {
          RHEA_ASSERT (
              fabs (rhea_inversion_param_project_to_feasible_weak (p)) <= 0.0);
          diff = rhea_inversion_param_project_to_feasible_weak (p + s);
          success = 1;
        }
      }

      /* otherwise the parameter index is unknown */
      if (!success) {
        RHEA_ABORT_NOT_REACHED ();
      }

      /* reduce step length */
      if (0.0 < fabs (diff)) {
        step_length = SC_MIN ((s + diff)/s, step_length);
        RHEA_ASSERT (0.0 <= step_length);
      }
    }
  }

  /* reduce step */
  if (step_length < 1.0) {
    ymir_vec_scale (step_length, step_vec);
  }

  /* return new step length */
  return step_length;
}

void
rhea_inversion_param_pull_from_model (ymir_vec_t *parameter_vec,
                                      rhea_inversion_param_t *inv_param)
{
  RHEA_GLOBAL_VERBOSEF_FN_TAG (__func__, "n_active=%i", inv_param->n_active);

  rhea_inversion_param_get_model_vals (parameter_vec, inv_param);
  rhea_inversion_param_convert_model_vals_to_params (parameter_vec, inv_param);
}

void
rhea_inversion_param_push_to_model (ymir_vec_t *parameter_vec,
                                    rhea_inversion_param_t *inv_param)
{
  ymir_vec_t         *model_vals = rhea_inversion_param_vec_new (inv_param);

  RHEA_GLOBAL_VERBOSEF_FN_TAG (__func__, "n_active=%i", inv_param->n_active);

  ymir_vec_copy (parameter_vec, model_vals);
  rhea_inversion_param_convert_params_to_model_vals (model_vals, inv_param);
  rhea_inversion_param_set_model_vals (model_vals, inv_param);
  rhea_inversion_param_vec_destroy (model_vals);
}

void
rhea_inversion_param_set_initial_from_model (
                                          ymir_vec_t *parameter_vec,
                                          rhea_inversion_param_t *inv_param,
                                          rhea_inversion_param_options_t *opt)
{
  const double        perturb_stddev = opt->initial_guess_perturb_stddev;
  const double        shift_by_prior_stddev =
                        opt->initial_guess_shift_by_prior_stddev;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (
      __func__, "perturb_stddev=%g, shift_by_prior_stddev=%g",
      perturb_stddev, shift_by_prior_stddev);

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (parameter_vec, inv_param));

  /* get parameter values from model */
  rhea_inversion_param_pull_from_model (parameter_vec, inv_param);

  /* apply perturbation */
  if (isfinite (perturb_stddev) && 0.0 < fabs (perturb_stddev)) {
    ymir_vec_t         *perturb;

    perturb = rhea_inversion_param_vec_new_perturb (inv_param, perturb_stddev);
    ymir_vec_multiply_in (perturb, parameter_vec);
    rhea_inversion_param_vec_destroy (perturb);
  }

  /* apply shift */
  if (isfinite (shift_by_prior_stddev) && 0.0 < fabs (shift_by_prior_stddev)) {
    ymir_vec_t         *shift_vec = rhea_inversion_param_vec_new (inv_param);

    /* copy prior standard deviations; adjust the sign of some values */
    ymir_vec_copy (inv_param->prior_stddev, shift_vec);
    {
      double             *shift = shift_vec->meshfree->e[0];
      const int          *active = inv_param->active;
      int                 idx;

      /* adjust upper bounds of the viscosity */
      idx = RHEA_INVERSION_PARAM_VISC_MAX;
      if (active[idx]) shift[idx] *= -1.0;

      /* adjust activation energy in Arrhenius relationship */
      idx = RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY;
      if (active[idx]) shift[idx] *= -1.0;
      idx = RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY;
      if (active[idx]) shift[idx] *= -1.0;

      /* adjust stress exponent that governs strain rate weakening */
      idx = RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT;
      if (active[idx]) shift[idx] *= -1.0;
    }

    /* shift parameters */
    ymir_vec_add (shift_by_prior_stddev, shift_vec, parameter_vec);
    rhea_inversion_param_vec_destroy (shift_vec);
  }

  /* print parameters */
#ifdef RHEA_ENABLE_DEBUG
  RHEA_GLOBAL_VERBOSE ("========================================\n");
  RHEA_GLOBAL_VERBOSEF ("%s\n", __func__);
  RHEA_GLOBAL_VERBOSE ("----------------------------------------\n");
  rhea_inversion_param_vec_print (parameter_vec, inv_param);
  RHEA_GLOBAL_VERBOSE ("========================================\n");
#endif

  /* check output */
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (parameter_vec, inv_param));

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static rhea_viscosity_param_derivative_t
rhea_inversion_param_get_derivative_type (
                                    const rhea_inversion_param_idx_t param_idx,
                                    rhea_weakzone_label_t *weak_label)
{
  int                 offset;

  /* look up indices related to viscosity parameters */
  *weak_label = RHEA_WEAKZONE_LABEL_UNKNOWN;
  switch (param_idx) {
  case RHEA_INVERSION_PARAM_VISC_MIN:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_MIN;
  case RHEA_INVERSION_PARAM_VISC_MAX:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_MAX;

  case RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_UPPER_MANTLE_SCALING;
  case RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_UPPER_MANTLE_ACTIVATION_ENERGY;
  case RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_LOWER_MANTLE_SCALING;
  case RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_LOWER_MANTLE_ACTIVATION_ENERGY;

  case RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_STRESS_EXPONENT;
  case RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_YIELD_STRENGTH;

  default:
    break;
  }

  /* look up indices related to (general) weak zone parameters */
  switch (param_idx) {
  case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_NONE:
  case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_NONE:
  case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_NONE:
    *weak_label = RHEA_WEAKZONE_LABEL_CLASS_NONE;
    break;
  case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_SLAB:
  case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_SLAB:
  case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_SLAB:
    *weak_label = RHEA_WEAKZONE_LABEL_CLASS_SLAB;
    break;
  case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_RIDGE:
  case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_RIDGE:
  case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_RIDGE:
    *weak_label = RHEA_WEAKZONE_LABEL_CLASS_RIDGE;
    break;
  case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_FRACTURE:
  case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_FRACTURE:
  case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE:
    *weak_label = RHEA_WEAKZONE_LABEL_CLASS_FRACTURE;
    break;
  default:
    break;
  }
  switch (param_idx) {
  case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_NONE:
  case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_SLAB:
  case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_RIDGE:
  case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_FRACTURE:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_WEAK_THICKNESS;
  case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_NONE:
  case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_SLAB:
  case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_RIDGE:
  case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_FRACTURE:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_WEAK_THICKNESS_CONST;
  case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_NONE:
  case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_SLAB:
  case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_RIDGE:
  case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE:
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_WEAK_FACTOR_INTERIOR;
  default:
    break;
  }

  /* look up indices related to earth's weak zone parameters */
  offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB;
  if (offset <= param_idx && param_idx < offset + RHEA_WEAKZONE_LABEL_EARTH_N) {
    *weak_label = rhea_weakzone_label_earth_get_label (param_idx - offset);
    return RHEA_VISCOSITY_PARAM_DERIVATIVE_WEAK_FACTOR_INTERIOR;
  }

  /* unknown derivative type */
  RHEA_ABORT_NOT_REACHED ();
  return RHEA_VISCOSITY_PARAM_DERIVATIVE_NONE;
}

void
rhea_inversion_param_print (ymir_vec_t *parameter_vec,
                            const rhea_inversion_param_verbosity_t verbosity,
                            rhea_inversion_param_t *inv_param)
{
  const double       *param = parameter_vec->meshfree->e[0];
  const int          *active = inv_param->active;
  int                 i;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (parameter_vec, inv_param));

  /* exit if nothing to do */
  if (!inv_param->n_active || RHEA_INVERSION_PARAM_VERBOSE_NONE == verbosity) {
    return;
  }

  /* print parameters */
  switch (verbosity) {
  case RHEA_INVERSION_PARAM_VERBOSE_REAL:
    rhea_inversion_param_vec_print (parameter_vec, inv_param);
    break;
  case RHEA_INVERSION_PARAM_VERBOSE_REAL_NONDIM:
    {
      ymir_vec_t         *model_vals = rhea_inversion_param_vec_new (inv_param);
      double             *model = model_vals->meshfree->e[0];

      rhea_inversion_param_get_model_vals (model_vals, inv_param);
      for (i = 0; i < inv_param->n_parameters; i++) {
        if (active[i]) {
          RHEA_GLOBAL_INFOF ("param# %3i: %g, %g\n", i, param[i], model[i]);
        }
      }
      rhea_inversion_param_vec_destroy (model_vals);
    }
    break;
  case RHEA_INVERSION_PARAM_VERBOSE_REAL_NONDIM_DIM:
    RHEA_ABORT_NOT_REACHED (); //TODO
    break;
  default: /* unknown verbosity */
    RHEA_ABORT_NOT_REACHED ();
  }
}

/******************************************************************************
 * Prior Generation
 *****************************************************************************/

/**
 * Perturbs the mean of prior multiplicatively.
 */
static void
rhea_inversion_param_prior_perturb_mean (rhea_inversion_param_t *inv_param,
                                         const double perturb_stddev)
{
  ymir_vec_t         *perturb;

  /* create perturbation vector */
  perturb = rhea_inversion_param_vec_new_perturb (inv_param, perturb_stddev);

  /* print perturbation */
#ifdef RHEA_ENABLE_DEBUG
  RHEA_GLOBAL_VERBOSE ("========================================\n");
  RHEA_GLOBAL_VERBOSEF ("%s\n", __func__);
  RHEA_GLOBAL_VERBOSE ("----------------------------------------\n");
  rhea_inversion_param_vec_print (perturb, inv_param);
  RHEA_GLOBAL_VERBOSE ("========================================\n");
#endif

  /* multiply prior mean by perturbation */
  ymir_vec_multiply_in (perturb, inv_param->prior_mean);
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (inv_param->prior_mean,
                                                  inv_param));

  /* destroy */
  rhea_inversion_param_vec_destroy (perturb);
}

static void
rhea_inversion_param_prior_get_mean (ymir_vec_t *prior_mean_vec,
                                     rhea_inversion_param_t *inv_param,
                                     rhea_inversion_param_options_t *opt)
{
  const double        perturb_stddev = opt->prior_mean_perturb_stddev;

  /* set prior mean of parameters to model values */
  rhea_inversion_param_get_model_vals (prior_mean_vec, inv_param);
  rhea_inversion_param_convert_model_vals_to_params (prior_mean_vec,
                                                     inv_param);

  /* perturb mean */
  if (isfinite (perturb_stddev) && 0.0 < fabs (perturb_stddev)) {
    rhea_inversion_param_prior_perturb_mean (inv_param, perturb_stddev);
  }

  /* print mean */
#ifdef RHEA_ENABLE_DEBUG
  RHEA_GLOBAL_VERBOSE ("========================================\n");
  RHEA_GLOBAL_VERBOSEF ("%s\n", __func__);
  RHEA_GLOBAL_VERBOSE ("----------------------------------------\n");
  rhea_inversion_param_vec_print (prior_mean_vec, inv_param);
  RHEA_GLOBAL_VERBOSE ("========================================\n");
#endif
}

static void
rhea_inversion_param_prior_get_stddev (ymir_vec_t *prior_stddev_vec,
                                       rhea_inversion_param_t *inv_param,
                                       rhea_inversion_param_options_t *opt)
{
  double             *p = prior_stddev_vec->meshfree->e[0];
  int                 idx;

  /* initialize all parameters to zero */
  ymir_vec_set_zero (prior_stddev_vec);

  /*
   * Viscosity Options
   */

  /* set lower and upper bounds of the viscosity */
  idx = RHEA_INVERSION_PARAM_VISC_MIN;
  p[idx] = opt->prior_stddev_min;
  idx = RHEA_INVERSION_PARAM_VISC_MAX;
  p[idx] = opt->prior_stddev_max;

  /* set scaling factors and activation energy in Arrhenius relationship */
  idx = RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING;
  p[idx] = opt->prior_stddev_upper_mantle_scaling;
  idx = RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY;
  p[idx] = opt->prior_stddev_upper_mantle_arrhenius_activation_energy;
  idx = RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING;
  p[idx] = opt->prior_stddev_lower_mantle_scaling;
  idx = RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY;
  p[idx] = opt->prior_stddev_lower_mantle_arrhenius_activation_energy;

  /* set stress exponent that governs strain rate weakening */
  idx = RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT;
  p[idx] = opt->prior_stddev_stress_exponent;

  /* set yield strength */
  idx = RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH;
  p[idx] = opt->prior_stddev_yield_strength;

  /*
   * Weak Zone Options
   */
  {
    const double        stddev_thickness = opt->prior_stddev_thickness;
    const double        stddev_thickness_const =
                          opt->prior_stddev_thickness_const;
    const double        stddev_weak_factor_interior =
                          opt->prior_stddev_weak_factor_interior;
    int                 offset;

    /* set weak zone thickness */
    idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_NONE;
    p[idx] = stddev_thickness;
    idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_SLAB;
    p[idx] = stddev_thickness;
    idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_RIDGE;
    p[idx] = stddev_thickness;
    idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_FRACTURE;
    p[idx] = stddev_thickness;

    /* set weak zone interior thickness */
    idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_NONE;
    p[idx] = stddev_thickness_const;
    idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_SLAB;
    p[idx] = stddev_thickness_const;
    idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_RIDGE;
    p[idx] = stddev_thickness_const;
    idx = RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_FRACTURE;
    p[idx] = stddev_thickness_const;

    /* set weak zone factors */
    idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_NONE;
    p[idx] = stddev_weak_factor_interior;
    idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_SLAB;
    p[idx] = stddev_weak_factor_interior;
    idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_RIDGE;
    p[idx] = stddev_weak_factor_interior;
    idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE;
    p[idx] = stddev_weak_factor_interior;

    /* set earth's weak zone factors */
    offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_EARTH_SLAB;
    for (idx = offset; idx < offset + RHEA_WEAKZONE_LABEL_EARTH_N; idx++) {
      p[idx] = stddev_weak_factor_interior;
    }
  }

  /* print standard deviation */
#ifdef RHEA_ENABLE_DEBUG
  RHEA_GLOBAL_VERBOSE ("========================================\n");
  RHEA_GLOBAL_VERBOSEF ("%s\n", __func__);
  RHEA_GLOBAL_VERBOSE ("----------------------------------------\n");
  rhea_inversion_param_vec_print (prior_stddev_vec, inv_param);
  RHEA_GLOBAL_VERBOSE ("========================================\n");
#endif
}

/******************************************************************************
 * Parameter Related Computations
 *****************************************************************************/

/**
 * Calculates misfit w.r.t. prior mean:
 *   misfit = param_prior_mean - param
 */
static void
rhea_inversion_param_prior_misfit (ymir_vec_t *misfit_vec,
                                   ymir_vec_t *parameter_vec,
                                   rhea_inversion_param_t *inv_param)
{
  if (NULL != parameter_vec) {
    RHEA_ASSERT (rhea_inversion_param_vec_check_type (parameter_vec,
                                                      inv_param));
    ymir_vec_copy (parameter_vec, misfit_vec);
  }
  else {
    ymir_vec_set_zero (misfit_vec);
  }
  ymir_vec_add (-1.0, inv_param->prior_mean, misfit_vec);
}

double
rhea_inversion_param_prior (ymir_vec_t *parameter_vec,
                            rhea_inversion_param_t *inv_param)
{
  ymir_vec_t         *prior_misfit = rhea_inversion_param_vec_new (inv_param);
  ymir_vec_t         *prior_icov = rhea_inversion_param_vec_new (inv_param);
  const int           n_parameters = inv_param->n_parameters;
  const int          *active = inv_param->active;
  const double       *misfit = prior_misfit->meshfree->e[0];
  const double       *icov = prior_icov->meshfree->e[0];
  int                 i;
  double              obj_val = 0.0;

  /* calculate misfit w.r.t. prior mean */
  rhea_inversion_param_prior_misfit (prior_misfit, parameter_vec, inv_param);

  /* calculate inverse prior covariance */
  ymir_vec_copy (inv_param->prior_stddev, prior_icov);
  ymir_vec_multiply_in (inv_param->prior_stddev, prior_icov);
  ymir_vec_reciprocal (prior_icov);

  /* calculate inner product */
  for (i = 0; i < n_parameters; i++) {
    if (active[i]) {
      RHEA_ASSERT (isfinite (misfit[i]));
      RHEA_ASSERT (isfinite (icov[i]));
      obj_val += misfit[i] * icov[i] * misfit[i];
    }
  }

  /* destroy */
  rhea_inversion_param_vec_destroy (prior_misfit);
  rhea_inversion_param_vec_destroy (prior_icov);

  /* return squared norm of prior */
  return obj_val;
}

void
rhea_inversion_param_compute_gradient (ymir_vec_t *gradient_vec,
                                       ymir_vec_t *parameter_vec,
                                       ymir_vec_t *forward_vel_press,
                                       ymir_vec_t *adjoint_vel_press,
                                       const double prior_weight,
                                       rhea_inversion_param_t *inv_param,
                                       ymir_vec_t *gradient_adjoint_comp_vec,
                                       ymir_vec_t *gradient_prior_comp_vec)
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
  double             *grad = gradient_vec->meshfree->e[0];
  double             *grad_adj = NULL;
  double             *grad_prior = NULL;
  double              g;
  int                 i;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (gradient_vec, inv_param));
  RHEA_ASSERT (parameter_vec == NULL ||
               rhea_inversion_param_vec_check_type (parameter_vec, inv_param));
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
  ymir_vec_set_zero (gradient_vec);
  if (gradient_adjoint_comp_vec != NULL) {
    ymir_vec_set_zero (gradient_adjoint_comp_vec);
    grad_adj = gradient_adjoint_comp_vec->meshfree->e[0];
  }
  if (gradient_prior_comp_vec != NULL) {
    ymir_vec_set_zero (gradient_prior_comp_vec);
    grad_prior = gradient_prior_comp_vec->meshfree->e[0];
  }
  for (i = 0; i < n_parameters; i++) { /* loop over all (possible) parameters */
    if (active[i]) {
      const rhea_inversion_param_idx_t  idx = (rhea_inversion_param_idx_t) i;
      rhea_viscosity_param_derivative_t deriv_type;
      rhea_weakzone_label_t             weak_label;

      RHEA_GLOBAL_VERBOSEF_FN_TAG (__func__, "param_idx=%i", i);

      /* compute parameter derivative of viscosity */
      deriv_type = rhea_inversion_param_get_derivative_type (idx, &weak_label);
      if (RHEA_WEAKZONE_LABEL_UNKNOWN == weak_label) {
        rhea_viscosity_param_derivative (
            derivative, deriv_type, viscosity, bounds_marker, yielding_marker,
            temperature, forward_vel, inv_param->visc_options);
      }
      else {
        rhea_viscosity_param_derivative_weakzone (
            derivative, deriv_type, weak_label,
            viscosity, bounds_marker, yielding_marker, weakzone,
            inv_param->weak_options, inv_param->visc_options);
      }

      /* transform to viscous stress coefficient */
      ymir_vec_scale (2.0, derivative);

      /* set derivative as the coefficient of the viscous stress operator */
      ymir_stress_op_set_coeff_scal (stress_op, derivative);

      /* apply viscous stress operator to forward velocity */
      ymir_stress_pc_apply_stress_op (forward_vel, op_out_vel, stress_op,
                                      0 /* !linearized */, 1 /* dirty */);

      /* compute inner product with adjoint velocity */
      g = ymir_vec_innerprod (op_out_vel, adjoint_vel);
      grad[i] = g;
      if (grad_adj != NULL) {
        grad_adj[i] = g;
      }
    }
  }
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (gradient_vec, inv_param));

  /* compute & add prior term */
  if (isfinite (prior_weight) && 0.0 < prior_weight) {
    ymir_vec_t         *prior_misfit = rhea_inversion_param_vec_new (inv_param);
    ymir_vec_t         *prior_icov = rhea_inversion_param_vec_new (inv_param);
    const double       *misfit = prior_misfit->meshfree->e[0];
    const double       *icov = prior_icov->meshfree->e[0];

    /* calculate misfit w.r.t. prior mean */
    rhea_inversion_param_prior_misfit (prior_misfit, parameter_vec, inv_param);

    /* calculate inverse prior covariance */
    ymir_vec_copy (inv_param->prior_stddev, prior_icov);
    ymir_vec_multiply_in (inv_param->prior_stddev, prior_icov);
    ymir_vec_reciprocal (prior_icov);

    /* add prior term to gradient */
    for (i = 0; i < n_parameters; i++) {
      if (active[i]) {
        RHEA_ASSERT (isfinite (misfit[i]));
        g = prior_weight * icov[i] * misfit[i];
        grad[i] += g;
        if (grad_prior != NULL) {
          grad_prior[i] = g;
        }
      }
    }

    /* destroy */
    rhea_inversion_param_vec_destroy (prior_misfit);
    rhea_inversion_param_vec_destroy (prior_icov);
  }
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (gradient_vec, inv_param));

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
rhea_inversion_param_compute_gradient_norm (ymir_vec_t *gradient_vec,
                                            rhea_inversion_param_t *inv_param)
{
  const int           n_parameters = inv_param->n_parameters;
  const int          *active = inv_param->active;
  const double       *grad = gradient_vec->meshfree->e[0];
  double              sum_of_squares = 0.0;
  int                 i;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (gradient_vec, inv_param));

  /* sum up (squares of) entries of the gradient vector that are active */
  for (i = 0; i < n_parameters; i++) { /* loop over all (possible) parameters */
    if (active[i]) {
      sum_of_squares += grad[i] * grad[i];
    }
  }

  /* return l2-norm of active entries of the gradient vector */
  return sqrt (sum_of_squares / (double) inv_param->n_active);
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
    if (active[i] && 0.0 < fabs (grad_dir[i])) {
      const rhea_inversion_param_idx_t  idx = (rhea_inversion_param_idx_t) i;
      rhea_viscosity_param_derivative_t deriv_type;
      rhea_weakzone_label_t             weak_label;

      /* compute parameter derivative of viscosity */
      deriv_type = rhea_inversion_param_get_derivative_type (idx, &weak_label);
      if (RHEA_WEAKZONE_LABEL_UNKNOWN == weak_label) {
        rhea_viscosity_param_derivative (
            derivative, deriv_type, viscosity, bounds_marker, yielding_marker,
            temperature, forward_vel, inv_param->visc_options);
      }
      else {
        rhea_viscosity_param_derivative_weakzone (
            derivative, deriv_type, weak_label,
            viscosity, bounds_marker, yielding_marker, weakzone,
            inv_param->weak_options, inv_param->visc_options);
      }

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
                                    const double prior_weight,
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

  /* init output by "hijacking" the gradient function using incr. adj. vel. */
  rhea_inversion_param_compute_gradient (
      param_vec_out, param_vec_in, forward_vel_press, incr_adjoint_vel_press,
      NAN /* no prior term */, inv_param, NULL, NULL);

  /* compute & add prior term */
  if (isfinite (prior_weight) && 0.0 < prior_weight) {
    ymir_vec_t         *prior_term = rhea_inversion_param_vec_new (inv_param);
    ymir_vec_t         *prior_icov = rhea_inversion_param_vec_new (inv_param);
    const int           n_parameters = inv_param->n_parameters;
    const int          *active = inv_param->active;
    const double       *prterm = prior_term->meshfree->e[0];
    const double       *icov = prior_icov->meshfree->e[0];
    double             *param_out = param_vec_out->meshfree->e[0];
    int                 i;

    /* copy input */
    ymir_vec_copy (param_vec_in, prior_term);

    /* apply inverse prior covariance */
    ymir_vec_copy (inv_param->prior_stddev, prior_icov);
    ymir_vec_multiply_in (inv_param->prior_stddev, prior_icov);
    ymir_vec_reciprocal (prior_icov);

    /* add prior term to Hessian apply output */
    for (i = 0; i < n_parameters; i++) {
      if (active[i]) {
        RHEA_ASSERT (isfinite (prterm[i]));
        param_out[i] += prior_weight * icov[i] * prterm[i];
      }
    }

    /* destroy */
    rhea_inversion_param_vec_destroy (prior_term);
    rhea_inversion_param_vec_destroy (prior_icov);
  }

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

/******************************************************************************
 * Parameter Vector
 *****************************************************************************/

ymir_vec_t *
rhea_inversion_param_vec_new (rhea_inversion_param_t *inv_param)
{
  return ymir_vec_new_meshfree (inv_param->n_parameters);
}

ymir_vec_t *
rhea_inversion_param_vec_new_perturb (rhea_inversion_param_t *inv_param,
                                      const double perturb_stddev)
{
  ymir_mesh_t        *ymir_mesh =
    rhea_stokes_problem_get_ymir_mesh (inv_param->stokes_problem);
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  ymir_vec_t         *perturb = rhea_inversion_param_vec_new (inv_param);

  /* set random values; add one */
  ymir_vec_set_random_normal (perturb, fabs (perturb_stddev), 0.0 /* mean */);
  ymir_vec_shift (1.0, perturb);

  /* bound from below to a value greater zero */
  ymir_vec_bound_min (perturb, SC_1000_EPS);

  /* sync values across all processors */
  ymir_vec_meshfree_sync (perturb, 0 /* mpirank_master */, mpicomm);

  return perturb;
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
  double             *v = vec->meshfree->e[0];
  int                 i;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (vec, inv_param));

  /* inspect processor-local values */
  for (i = 0; i < n_parameters; i++) { /* loop over all (possible) parameters */
    if (active[i] && !isfinite (v[i])) {
      return 0;
    }
  }

  /* inspect values across processors (using a ring topology) */
  {
    ymir_mesh_t        *ymir_mesh =
      rhea_stokes_problem_get_ymir_mesh (inv_param->stokes_problem);
    sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
    const int           mpisize = ymir_mesh_get_MPI_Comm_size (ymir_mesh);
    const int           mpirank = ymir_mesh_get_MPI_Comm_rank (ymir_mesh);
    int                 mpiret;
    int                 mpirank_prev, mpirank_next;
    const int           THIS_MPI_TAG = 1;
    sc_MPI_Request      request_recv, request_send;
    sc_MPI_Status       status;
    double             *v_prev;

    /* return success if nothing to do */
    if (mpisize <= 1) {
      return 1;
    }

    /* set previous and next rank neighbors */
    mpirank_prev = mpirank - 1;
    mpirank_next = mpirank + 1;
    if (0 == mpirank) {
      mpirank_prev = mpisize - 1;
    }
    if (mpirank == mpisize - 1) {
      mpirank_next = 0;
    }

    /* receive values from previous neighbor */
    v_prev = RHEA_ALLOC (double, n_parameters);
    mpiret = sc_MPI_Irecv (v_prev, n_parameters, sc_MPI_DOUBLE, mpirank_prev,
                           THIS_MPI_TAG, mpicomm, &request_recv);
    SC_CHECK_MPI (mpiret);

    /* send values to next neighbor */
    mpiret = sc_MPI_Isend (v, n_parameters, sc_MPI_DOUBLE, mpirank_next,
                           THIS_MPI_TAG, mpicomm, &request_send);
    SC_CHECK_MPI (mpiret);

    /* wait until receive is complete */
    mpiret = sc_MPI_Waitall (1, &request_recv, &status);
    SC_CHECK_MPI (mpiret);

    /* compare values of previous neighbor */
    for (i = 0; i < n_parameters; i++) { /* loop over all (possible) params */
      const double        rel_err = fabs (v[i] - v_prev[i]) / fabs (v[i]);

      if (active[i] && SC_1000_EPS < rel_err) {
        return 0;
      }
    }
    RHEA_FREE (v_prev);

    /* wait until send is complete */
    mpiret = sc_MPI_Waitall (1, &request_send, &status);
    SC_CHECK_MPI (mpiret);
  }

  /* return success */
  return 1;
}

void
rhea_inversion_param_vec_print (ymir_vec_t *parameter_vec,
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
