#include <rhea_inversion_param.h>
#include <rhea_base.h>
#include <rhea_viscosity_param_derivative.h>
#include <rhea_weakzone_label.h>
#include <rhea_stress.h>
#include <rhea_io_mpi.h>
#include <ymir_stress_pc.h>
#include <fenv.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_INVERSION_PARAM_DEFAULT_ACTIVATE 1
#define RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE 0
#define RHEA_INVERSION_PARAM_DEFAULT_ACTIVATE_WEAK_LABEL_FILE_PATH_TXT NULL
#define RHEA_INVERSION_PARAM_DEFAULT_PRMN_PERTURB_STDDEV (NAN)
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_MIN (1.0e-1)
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_MAX (1.0e+1)
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_UM_SCALING (100.0)
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_UM_ARRHENIUS_ACTIVATION_ENERGY (4.0)
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_LM_SCALING (100.0)
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_LM_ARRHENIUS_ACTIVATION_ENERGY (4.0)
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_STRESS_EXPONENT (0.5)
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_YIELD_STRENGTH (4.6)
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_THICKNESS (10.0e3)
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_THICKNESS_CONST (5.0e3)
#define RHEA_INVERSION_PARAM_DEFAULT_PRSD_WEAK_FACTOR_INTERIOR (1.0e-1)
#define RHEA_INVERSION_PARAM_DEFAULT_GUESS_FILE_PATH_TXT NULL
#define RHEA_INVERSION_PARAM_DEFAULT_GUESS_PERTURB_STDDEV (0.0)
#define RHEA_INVERSION_PARAM_DEFAULT_GUESS_SHIFT_BY_PRIOR_STDDEV (0.0)

#define RHEA_INVERSION_PARAM_DEFAULT_DIMENSIONAL_SCALING (1)
#define RHEA_INVERSION_PARAM_DEFAULT_GRADIENT_NORM_WEIGHT_ICOV (0)

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

  YMIR_OPTIONS_B, "activate-weak-factor-interior-label-slab", '\0',
    &(inv_param_opt->weak_factor_interior_label_slab_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone factors of slabs",
  YMIR_OPTIONS_B, "activate-weak-factor-interior-label-ridge", '\0',
    &(inv_param_opt->weak_factor_interior_label_ridge_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone factors of ridges",
  YMIR_OPTIONS_B, "activate-weak-factor-interior-label-fracture", '\0',
    &(inv_param_opt->weak_factor_interior_label_fracture_a),
    RHEA_INVERSION_PARAM_DEFAULT_DEACTIVATE,
    "Activate weak zone factors of fractures",

  YMIR_OPTIONS_S, "activate-weak-factor-interior-label-file-path-txt", '\0',
    &(inv_param_opt->weak_factor_interior_label_file_path_txt_a),
    RHEA_INVERSION_PARAM_DEFAULT_ACTIVATE_WEAK_LABEL_FILE_PATH_TXT,
    "Path to a text file with weakzone activation flags "
    "(one entry for each label)",

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

  YMIR_OPTIONS_S, "initial-guess-file-path-txt", '\0',
    &(inv_param_opt->initial_guess_file_path_txt),
    RHEA_INVERSION_PARAM_DEFAULT_GUESS_FILE_PATH_TXT,
    "Path to a text file of inversion parameters for the initial guess",

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

  /* parameterization of max interior weakening factor based on classes */
  RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_NONE,
  RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_SLAB,
  RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_RIDGE,
  RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE,

  /* max number of global parameters (is greater than actual #parameters) */
  RHEA_INVERSION_PARAM_N_GLOBAL = 32,

  /* parameterization of max interior weakening factor based on labels */
  RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_NONE =
    1 + RHEA_INVERSION_PARAM_N_GLOBAL,
  RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_SLAB =
    1 + RHEA_INVERSION_PARAM_N_GLOBAL +
    RHEA_WEAKZONE_LABEL_MAX_N_NONE,
  RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_RIDGE =
    1 + RHEA_INVERSION_PARAM_N_GLOBAL +
    RHEA_WEAKZONE_LABEL_MAX_N_NONE +
    RHEA_WEAKZONE_LABEL_MAX_N_SLAB,
  RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_FRACTURE =
    1 + RHEA_INVERSION_PARAM_N_GLOBAL +
    RHEA_WEAKZONE_LABEL_MAX_N_NONE +
    RHEA_WEAKZONE_LABEL_MAX_N_SLAB +
    RHEA_WEAKZONE_LABEL_MAX_N_RIDGE,

  /* total number of global and local (i.e., weak zone) parameters */
  RHEA_INVERSION_PARAM_N =
    1 + RHEA_INVERSION_PARAM_N_GLOBAL +
    RHEA_WEAKZONE_LABEL_MAX_N_NONE +
    RHEA_WEAKZONE_LABEL_MAX_N_SLAB +
    RHEA_WEAKZONE_LABEL_MAX_N_RIDGE +
    RHEA_WEAKZONE_LABEL_MAX_N_FRACTURE
}
rhea_inversion_param_idx_t;

/* set offset and count of weak factors by label */
#define RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_OFFSET \
  RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_NONE
#define RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_COUNT \
  RHEA_WEAKZONE_LABEL_MAX_N_NONE + \
  RHEA_WEAKZONE_LABEL_MAX_N_SLAB + \
  RHEA_WEAKZONE_LABEL_MAX_N_RIDGE + \
  RHEA_WEAKZONE_LABEL_MAX_N_FRACTURE

static int
rhea_inversion_param_idx_is_weak_label (const int idx)
{
  const int offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_OFFSET;
  const int count  = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_COUNT;

  return (offset <= idx && idx < offset + count);
}

static int
rhea_inversion_param_idx_lookup_offset (const int idx)
{
  int                 range_begin, range_end;

  range_begin = 0;
  range_end   = RHEA_INVERSION_PARAM_N_GLOBAL;
  if (range_begin <= idx && idx < range_end) {
    return range_begin;
  }

  range_begin = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_NONE;
  range_end   = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_SLAB;
  if (range_begin <= idx && idx < range_end) {
    return range_begin;
  }

  range_begin = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_SLAB;
  range_end   = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_RIDGE;
  if (range_begin <= idx && idx < range_end) {
    return range_begin;
  }

  range_begin = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_RIDGE;
  range_end   = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_FRACTURE;
  if (range_begin <= idx && idx < range_end) {
    return range_begin;
  }

  range_begin = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_FRACTURE;
  range_end   = RHEA_INVERSION_PARAM_N;
  if (range_begin <= idx && idx < range_end) {
    return range_begin;
  }

  /* otherwise unknown parameter index */
  return -1;
}

static int
rhea_inversion_param_lookup_offset_from_weak_label (
                                            const rhea_weakzone_label_t label)
{
  const rhea_weakzone_label_t class_id = rhea_weakzone_label_get_class (label);

  switch (class_id) {
  case RHEA_WEAKZONE_LABEL_CLASS_NONE:
    return RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_NONE;
  case RHEA_WEAKZONE_LABEL_CLASS_SLAB:
    return RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_SLAB;
  case RHEA_WEAKZONE_LABEL_CLASS_RIDGE:
    return RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_RIDGE;
  case RHEA_WEAKZONE_LABEL_CLASS_FRACTURE:
    return RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_FRACTURE;
  default: /* unknown label */
    return -1;
  }
}

/******************************************************************************
 * Parameter Activation
 *****************************************************************************/

static int *
rhea_inversion_param_activation_mask_new (
                                        rhea_inversion_param_options_t *opt,
                                        rhea_weakzone_options_t *weak_options,
                                        rhea_viscosity_options_t *visc_options,
                                        sc_MPI_Comm mpicomm)
{
  int                *active = RHEA_ALLOC_ZERO (int, RHEA_INVERSION_PARAM_N);
  int                 offset, count, idx;

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

    /* process weak zone thickness of the interior */
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
    if (NULL != opt->weak_factor_interior_label_file_path_txt_a) {
      const char         *file_path_txt =
        opt->weak_factor_interior_label_file_path_txt_a;
      const int           n_entries =
        rhea_weakzone_get_total_n_labels (weak_options);
      int                 n_read, idx_weak;
      int                *activate_weakzone;

      /* read file */
      activate_weakzone = RHEA_ALLOC (int, n_entries);
      n_read = rhea_io_mpi_read_broadcast_int (
          activate_weakzone, n_entries, NULL /* path bin */, file_path_txt,
          mpicomm);
      RHEA_ASSERT (n_read == n_entries);

      /* set activation for each weakzone label individually */
      idx_weak = 0;

      offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_SLAB;
      count  = weak_options->n_labels[RHEA_WEAKZONE_LABEL_CLASS_SLAB];
      for (idx = offset; idx < offset + count; idx++) {
        RHEA_ASSERT (idx_weak < n_read);
        if (activate_weakzone[idx_weak]) {
          active[idx] = 1;
        }
        idx_weak++;
      }

      offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_RIDGE;
      count  = weak_options->n_labels[RHEA_WEAKZONE_LABEL_CLASS_RIDGE];
      for (idx = offset; idx < offset + count; idx++) {
        RHEA_ASSERT (idx_weak < n_read);
        if (activate_weakzone[idx_weak]) {
          active[idx] = 1;
        }
        idx_weak++;
      }

      offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_FRACTURE;
      count  = weak_options->n_labels[RHEA_WEAKZONE_LABEL_CLASS_FRACTURE];
      for (idx = offset; idx < offset + count; idx++) {
        RHEA_ASSERT (idx_weak < n_read);
        if (activate_weakzone[idx_weak]) {
          active[idx] = 1;
        }
        idx_weak++;
      }

      RHEA_FREE (activate_weakzone);
    }
    else {
      /* process weak zone factors as classes
       * (options for specific classes override general options) */
      if (!opt->weak_factor_interior_class_slab_a &&
          !opt->weak_factor_interior_class_ridge_a &&
          !opt->weak_factor_interior_class_fracture_a &&
          !opt->weak_factor_interior_label_slab_a &&
          !opt->weak_factor_interior_label_ridge_a &&
          !opt->weak_factor_interior_label_fracture_a) {
        idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_NONE;
        active[idx] = opt->weak_factor_interior_a;
      }
      if (!opt->weak_factor_interior_label_slab_a) {
        idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_SLAB;
        active[idx] = opt->weak_factor_interior_class_slab_a;
      }
      else {
        offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_SLAB;
        count  = weak_options->n_labels[RHEA_WEAKZONE_LABEL_CLASS_SLAB];
        for (idx = offset; idx < offset + count; idx++) {
          active[idx] = 1;
        }
      }
      if (!opt->weak_factor_interior_label_ridge_a) {
        idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_RIDGE;
        active[idx] = opt->weak_factor_interior_class_ridge_a;
      }
      else {
        offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_RIDGE;
        count  = weak_options->n_labels[RHEA_WEAKZONE_LABEL_CLASS_RIDGE];
        for (idx = offset; idx < offset + count; idx++) {
          active[idx] = 1;
        }
      }
      if (!opt->weak_factor_interior_label_fracture_a) {
        idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE;
        active[idx] = opt->weak_factor_interior_class_fracture_a;
      }
      else {
        offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_FRACTURE;
        count  = weak_options->n_labels[RHEA_WEAKZONE_LABEL_CLASS_FRACTURE];
        for (idx = offset; idx < offset + count; idx++) {
          active[idx] = 1;
        }
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
static void         rhea_inversion_param_set_dimensions (
                                          ymir_vec_t *dimensional_scaling_vec,
                                          rhea_inversion_param_t *inv_param);
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

  /* dimensional scaling factors */
  ymir_vec_t         *dimensional_scaling;

  /* prior data */
  ymir_vec_t         *prior_mean;
  ymir_vec_t         *prior_stddev;

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
  ymir_mesh_t        *ymir_mesh =
    rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  const int           use_dim =
    RHEA_INVERSION_PARAM_DEFAULT_DIMENSIONAL_SCALING;
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
      inv_param_opt, inv_param->weak_options, inv_param->visc_options,
      ymir_mesh_get_MPI_Comm (ymir_mesh));
  inv_param->n_active = rhea_inversion_param_activation_mask_count (
      inv_param->active);
  RHEA_ASSERT (0 < inv_param->n_active);

  /* create dimensional scaling factors from model parameters */
  if (use_dim) {
    inv_param->dimensional_scaling = rhea_inversion_param_vec_new (inv_param);
    rhea_inversion_param_set_dimensions (inv_param->dimensional_scaling,
                                         inv_param);
  }
  else {
    inv_param->dimensional_scaling = NULL;
  }

  /* create mean and standard deviation of the Gaussian prior */
  inv_param->prior_mean = rhea_inversion_param_vec_new (inv_param);
  inv_param->prior_stddev = rhea_inversion_param_vec_new (inv_param);
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
  if (inv_param->dimensional_scaling != NULL) {
    rhea_inversion_param_vec_destroy (inv_param->dimensional_scaling);
  }
  rhea_inversion_param_activation_mask_destroy (inv_param->active);
  RHEA_FREE (inv_param);
}

/**
 * Checks if weak factors that are distinguished by labels are active
 * parameters.
 */
static int
rhea_inversion_param_weak_label_is_active (rhea_inversion_param_t *inv_param)
{
  int                 offset, count, idx;
  int                 n = 0;

  offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_OFFSET;
  count  = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_COUNT;
  for (idx = offset; idx < offset + count; idx++) {
    n += inv_param->active[idx];
  }

  return n;
}

/******************************************************************************
 * Conversion Between Model Values and Their Corresponding Inversion Parameters.
 *
 * Denote a model value by `m` and an inversion parameter value by `p`, then
 *   m <- *_convert_to_model_* (p)
 *   p <- *_convert_from_model_* (m)
 *****************************************************************************/

static double
_diff_to_feasible (const double inv_param_val, const int idx,
                   double feasible_min, double feasible_max,
                   const double restrict_to_prior_stddev,
                   rhea_inversion_param_t *inv_param)
{
  const double       *mean = inv_param->prior_mean->meshfree->e[0];
  const double       *stddev = inv_param->prior_stddev->meshfree->e[0];
  double              prior_min, prior_max;

  if (isfinite (restrict_to_prior_stddev) && 0.0 < restrict_to_prior_stddev) {
    prior_min = mean[idx] - restrict_to_prior_stddev * stddev[idx];
    prior_max = mean[idx] + restrict_to_prior_stddev * stddev[idx];
    RHEA_ASSERT (isfinite (prior_min));
    RHEA_ASSERT (isfinite (prior_max));
    if (!isfinite (feasible_min) || feasible_min < prior_min) {
      feasible_min = prior_min;
    }
    if (!isfinite (feasible_max) || prior_max < feasible_max) {
      feasible_max = prior_max;
    }
  }

  if (isfinite (feasible_min) && inv_param_val < feasible_min) {
    return feasible_min - inv_param_val;
  }
  else if (isfinite (feasible_max) && feasible_max < inv_param_val) {
    return feasible_max - inv_param_val;
  }
  else {
    return 0.0;
  }
}

/****** General Conversion: positive number, linear scaling ******/

static double
_poslin_convert_to_model (const double inv_param_val)
{
  RHEA_ASSERT (isfinite (inv_param_val));
  RHEA_ASSERT (0.0 < inv_param_val);
  return inv_param_val;
}

static double
_poslin_convert_from_model (const double model_val)
{
  RHEA_ASSERT (isfinite (model_val));
  RHEA_ASSERT (0.0 < model_val);
  return model_val;
}

static double
_poslin_derivative (const double model_val)
{
  RHEA_ASSERT (isfinite (model_val));
  RHEA_ASSERT (0.0 < model_val);
  return 1.0;
}

/****** General Conversion: positive number, exponential scaling ******/

static double
_posexp_convert_to_model (const double inv_param_val)
{
  double              model_val;

  RHEA_ASSERT (isfinite (inv_param_val));

  feclearexcept (FE_ALL_EXCEPT);
  model_val = exp (inv_param_val);
  if (fetestexcept (FE_OVERFLOW)) { /* if argument of exp too large */
    model_val = DBL_MAX;
  }

  return model_val;
}

static double
_posexp_convert_from_model (const double model_val)
{
  double              inv_param_val;

  RHEA_ASSERT (0.0 <= model_val);

  inv_param_val = log (model_val);
  if (!isfinite (inv_param_val)) {
    inv_param_val = -DBL_MAX;
  }

  return inv_param_val;
}

static double
_posexp_derivative (const double model_val)
{
  RHEA_ASSERT (0.0 <= model_val);
  return model_val;
}

/****** Conversion: viscosity bounds ******/

double
rhea_inversion_param_derivative_min_max (const double model_val)
{
  return _poslin_derivative (model_val);
}

/****** Conversion: viscosity scaling ******/

static double
rhea_inversion_param_convert_to_model_scal (const double inv_param_val)
{
  return _posexp_convert_to_model (inv_param_val);
}

static double
rhea_inversion_param_convert_from_model_scal (const double model_val)
{
  return _posexp_convert_from_model (model_val);
}

double
rhea_inversion_param_derivative_scal (const double model_val)
{
  return _posexp_derivative (model_val);
}

static double
rhea_inversion_param_diff_to_feasible_scal (
                                        const double inv_param_val,
                                        const int idx,
                                        const double restrict_to_prior_stddev,
                                        rhea_inversion_param_t *inv_param)
{
  return _diff_to_feasible (inv_param_val, idx, NAN, NAN,
                            restrict_to_prior_stddev, inv_param);
}

/****** Conversion: Arrhenius activation energy ******/

static double
rhea_inversion_param_convert_to_model_arrh (const double inv_param_val)
{
  return _poslin_convert_to_model (inv_param_val);
}

static double
rhea_inversion_param_convert_from_model_arrh (const double model_val)
{
  return _poslin_convert_from_model (model_val);
}

double
rhea_inversion_param_derivative_arrh (const double model_val)
{
  return _poslin_derivative (model_val);
}

static double
rhea_inversion_param_diff_to_feasible_arrh (
                                        const double inv_param_val,
                                        const int idx,
                                        const double restrict_to_prior_stddev,
                                        rhea_inversion_param_t *inv_param)
{
  return _diff_to_feasible (inv_param_val, idx, SC_1000_EPS, NAN,
                            restrict_to_prior_stddev, inv_param);
}

/****** Conversion: stress exponent ******/

static double
rhea_inversion_param_convert_to_model_n (const double inv_param_val)
{
  return _poslin_convert_to_model (1.0/inv_param_val);
}

static double
rhea_inversion_param_convert_from_model_n (const double model_val)
{
  return _poslin_convert_from_model (1.0/model_val);
}

double
rhea_inversion_param_derivative_n (const double model_val)
{
  return _poslin_derivative (1.0/model_val);
}

static double
rhea_inversion_param_diff_to_feasible_n (const double inv_param_val,
                                         const int idx,
                                         const double restrict_to_prior_stddev,
                                         rhea_inversion_param_t *inv_param)
{
  return _diff_to_feasible (inv_param_val, idx, SC_1000_EPS, NAN,
                            restrict_to_prior_stddev, inv_param);
}

/****** Conversion: yield strenght ******/

static double
rhea_inversion_param_convert_to_model_yield (const double inv_param_val)
{
  return _poslin_convert_to_model (inv_param_val);
}

static double
rhea_inversion_param_convert_from_model_yield (const double model_val)
{
  return _poslin_convert_from_model (model_val);
}

double
rhea_inversion_param_derivative_yield (const double model_val)
{
  return _poslin_derivative (model_val);
}

static double
rhea_inversion_param_diff_to_feasible_yield (
                                        const double inv_param_val,
                                        const int idx,
                                        const double restrict_to_prior_stddev,
                                        rhea_inversion_param_t *inv_param)
{
  return _diff_to_feasible (inv_param_val, idx, SC_1000_EPS, NAN,
                            restrict_to_prior_stddev, inv_param);
}

/****** Conversion: weak zone factor ******/

static double
rhea_inversion_param_convert_to_model_weak (const double inv_param_val)
{
  RHEA_ASSERT (isfinite (inv_param_val));
  RHEA_ASSERT (inv_param_val <= 0.0);
  return _posexp_convert_to_model (inv_param_val);
}

static double
rhea_inversion_param_convert_from_model_weak (const double model_val)
{
  RHEA_ASSERT (0.0 < model_val && model_val <= 1.0);
  return _posexp_convert_from_model (model_val);
}

double
rhea_inversion_param_derivative_weak (const double model_val)
{
  RHEA_ASSERT (0.0 < model_val && model_val <= 1.0);
  return _posexp_derivative (model_val);
}

static double
rhea_inversion_param_diff_to_feasible_weak (
                                        const double inv_param_val,
                                        const int idx,
                                        const double restrict_to_prior_stddev,
                                        rhea_inversion_param_t *inv_param)
{
  return _diff_to_feasible (inv_param_val, idx, -DBL_MAX, -SC_1000_EPS,
                            restrict_to_prior_stddev, inv_param);
}

/******************************************************************************
 * Transfer between Parameter Values and Model Values
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

  /* get max weakening factor by class */
  idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_NONE;
  if (active[idx]) p[idx] = weak_options->weak_factor_interior;
  idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_SLAB;
  if (active[idx]) p[idx] = weak_options->weak_factor_interior_class_slab;
  idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_RIDGE;
  if (active[idx]) p[idx] = weak_options->weak_factor_interior_class_ridge;
  idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE;
  if (active[idx]) p[idx] = weak_options->weak_factor_interior_class_fracture;

  /* get max weakening factor by label */
  if (rhea_inversion_param_weak_label_is_active (inv_param)) {
    const double       *weak_label = weak_options->weak_factor_interior_label;
    int                 class_id, model_id;
    int                 offset, count;

    RHEA_ASSERT (weak_options->weak_factor_interior_label != NULL);

    model_id = 0;
    for (class_id = 0; class_id < RHEA_WEAKZONE_LABEL_CLASS_N; class_id++) {
      offset = rhea_inversion_param_lookup_offset_from_weak_label (
          (rhea_weakzone_label_t) class_id);
      count = weak_options->n_labels[class_id];
      for (idx = 0; idx < count; idx++) {
        RHEA_ASSERT (model_id <
                     rhea_weakzone_get_total_n_labels (weak_options));
        if (active[offset+idx]) p[offset+idx] = weak_label[model_id];
        model_id++;
      }
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

  /* set max weakening factor by class */
  idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_NONE;
  if (active[idx]) weak_options->weak_factor_interior = p[idx];
  idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_SLAB;
  if (active[idx]) weak_options->weak_factor_interior_class_slab = p[idx];
  idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_RIDGE;
  if (active[idx]) weak_options->weak_factor_interior_class_ridge = p[idx];
  idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE;
  if (active[idx]) weak_options->weak_factor_interior_class_fracture = p[idx];

  /* set max weakening factor by label */
  if (rhea_inversion_param_weak_label_is_active (inv_param)) {
    double             *weak_label = weak_options->weak_factor_interior_label;
    int                 class_id, model_id;
    int                 offset, count;

    RHEA_ASSERT (weak_options->weak_factor_interior_label != NULL);

    model_id = 0;
    for (class_id = 0; class_id < RHEA_WEAKZONE_LABEL_CLASS_N; class_id++) {
      offset = rhea_inversion_param_lookup_offset_from_weak_label (
          (rhea_weakzone_label_t) class_id);
      count = weak_options->n_labels[class_id];
      for (idx = 0; idx < count; idx++) {
        RHEA_ASSERT (model_id <
                     rhea_weakzone_get_total_n_labels (weak_options));
        if (active[offset+idx]) weak_label[model_id] = p[offset+idx];
        model_id++;
      }
    }
  }
}

static void
rhea_inversion_param_convert_model_vals_to_params (
                                            ymir_vec_t *parameter_vec,
                                            rhea_inversion_param_t *inv_param)
{
  double             *p = parameter_vec->meshfree->e[0];
  const int           dim_exists = (inv_param->dimensional_scaling != NULL);
  const double       *dim =
    (dim_exists ? inv_param->dimensional_scaling->meshfree->e[0] : NULL);
  const int          *active = inv_param->active;
  int                 idx;

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
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_FRACTURE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_FRACTURE:
        p[idx] = _poslin_convert_from_model (p[idx]);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING:
      case RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING:
        p[idx] = rhea_inversion_param_convert_from_model_scal (p[idx]);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY:
      case RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY:
        p[idx] = rhea_inversion_param_convert_from_model_arrh (p[idx]);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT:
        p[idx] = rhea_inversion_param_convert_from_model_n (p[idx]);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH:
        p[idx] = rhea_inversion_param_convert_from_model_yield (p[idx]);
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
      if (!success && rhea_inversion_param_idx_is_weak_label (idx)) {
        p[idx] = rhea_inversion_param_convert_from_model_weak (p[idx]);
        success = 1;
      }

      /* otherwise the parameter index is unknown */
      RHEA_CHECK_ABORT (success, "Parameter index is unknown.");

      /* remove dimensional scaling */
      if (dim_exists) {
        RHEA_ASSERT (0.0 < dim[idx]);
        p[idx] *= 1.0/dim[idx];
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
  const int           dim_exists = (inv_param->dimensional_scaling != NULL);
  const double       *dim =
    (dim_exists ? inv_param->dimensional_scaling->meshfree->e[0] : NULL);
  const int          *active = inv_param->active;
  int                 idx;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (parameter_vec, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (parameter_vec, inv_param));

  /* convert the values of active parameters */
  for (idx = 0; idx < inv_param->n_parameters; idx++) {
    int                 success = 0;

    if (active[idx]) {
      /* apply dimensional scaling */
      if (dim_exists) {
        RHEA_ASSERT (0.0 < dim[idx]);
        p[idx] *= dim[idx];
      }

      /* look up single index */
      switch (idx) {
      case RHEA_INVERSION_PARAM_VISC_MIN:
      case RHEA_INVERSION_PARAM_VISC_MAX:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_FRACTURE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_FRACTURE:
        p[idx] = _poslin_convert_to_model (p[idx]);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING:
      case RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING:
        p[idx] = rhea_inversion_param_convert_to_model_scal (p[idx]);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY:
      case RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY:
        p[idx] = rhea_inversion_param_convert_to_model_arrh (p[idx]);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT:
        p[idx] = rhea_inversion_param_convert_to_model_n (p[idx]);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH:
        p[idx] = rhea_inversion_param_convert_to_model_yield (p[idx]);
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
      if (!success && rhea_inversion_param_idx_is_weak_label (idx)) {
        p[idx] = rhea_inversion_param_convert_to_model_weak (p[idx]);
        success = 1;
      }

      /* otherwise the parameter index is unknown */
      RHEA_CHECK_ABORT (success, "Parameter index is unknown.");
    }
  }

  /* check output */
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (parameter_vec, inv_param));
}

static void
rhea_inversion_param_convert_model_vals_to_dim_model_vals (
                                            ymir_vec_t *parameter_vec,
                                            rhea_inversion_param_t *inv_param)
{
  rhea_viscosity_options_t   *visc_options = inv_param->visc_options;
  rhea_temperature_options_t *temp_options = visc_options->temp_options;
  rhea_domain_options_t      *domain_options = visc_options->domain_options;
  double             *p = parameter_vec->meshfree->e[0];
  const int          *active = inv_param->active;
  int                 idx;

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
      case RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_FRACTURE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_FRACTURE:
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE:
        /* keep same value */
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY:
      case RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY:
        p[idx] *= rhea_temperature_activation_energy_get_dim_J_mol (
            temp_options);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING:
        p[idx] *= rhea_viscosity_scaling_get_dim (
            visc_options->upper_mantle_arrhenius_activation_energy,
            visc_options);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING:
        p[idx] *= rhea_viscosity_scaling_get_dim (
            visc_options->lower_mantle_arrhenius_activation_energy,
            visc_options);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH:
        p[idx] *= rhea_stress_get_dim_Pa (
            domain_options, temp_options, visc_options);
        success = 1;
        break;
      default:
        break;
      }

      /* look up indices in a range */
      if (!success && rhea_inversion_param_idx_is_weak_label (idx)) {
        /* keep same value */
        success = 1;
      }

      /* otherwise the parameter index is unknown */
      RHEA_CHECK_ABORT (success, "Parameter index is unknown.");
    }
  }

  /* check output */
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (parameter_vec, inv_param));
}

static void
rhea_inversion_param_set_dimensions (ymir_vec_t *dimensional_scaling_vec,
                                     rhea_inversion_param_t *inv_param)
{
  double             *dim = dimensional_scaling_vec->meshfree->e[0];
  const int          *active = inv_param->active;
  int                 idx;
  ymir_vec_t         *restore_dimensional_scaling;

  /* initialize dimensional scaling to model values */
  restore_dimensional_scaling = inv_param->dimensional_scaling;
  inv_param->dimensional_scaling = NULL;
  rhea_inversion_param_get_model_vals (dimensional_scaling_vec, inv_param);
  rhea_inversion_param_convert_model_vals_to_params (dimensional_scaling_vec,
                                                     inv_param);
  inv_param->dimensional_scaling = restore_dimensional_scaling;

  /* set dimensional scaling */
  for (idx = 0; idx < inv_param->n_parameters; idx++) {
    if (active[idx] && 0.0 < fabs (dim[idx])) {
      dim[idx] = fabs (dim[idx]);
    }
    else {
      dim[idx] = 1.0;
    }
  }

  /* print dimensions */
  RHEA_GLOBAL_INFO ("========================================\n");
  RHEA_GLOBAL_INFOF ("%s\n", __func__);
  RHEA_GLOBAL_INFO ("----------------------------------------\n");
  rhea_inversion_param_vec_print (dimensional_scaling_vec, inv_param);
  RHEA_GLOBAL_INFO ("========================================\n");
}

int
rhea_inversion_param_restrict_to_feasible (
                                        ymir_vec_t *parameter_vec,
                                        const double restrict_to_prior_stddev,
                                        rhea_inversion_param_t *inv_param)
{
  double             *param = parameter_vec->meshfree->e[0];
  int                 param_modified = 0;
  const int          *active = inv_param->active;
  int                 idx;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (parameter_vec, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (parameter_vec, inv_param));

  for (idx = 0; idx < inv_param->n_parameters; idx++) {
    const double        p = param[idx];
    double              diff = 0.0;
    int                 success = 0;

    if (active[idx]) {
      /* look up single index */
      switch (idx) {
      case RHEA_INVERSION_PARAM_VISC_MIN:
      case RHEA_INVERSION_PARAM_VISC_MAX:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_FRACTURE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_FRACTURE:
        diff = _diff_to_feasible (
            p, idx, SC_1000_EPS, NAN, restrict_to_prior_stddev, inv_param);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING:
      case RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING:
        diff = rhea_inversion_param_diff_to_feasible_scal (
            p, idx, restrict_to_prior_stddev, inv_param);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY:
      case RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY:
        diff = rhea_inversion_param_diff_to_feasible_arrh (
            p, idx, restrict_to_prior_stddev, inv_param);
        success = 1;
        break;
        break;
      case RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT:
        diff = rhea_inversion_param_diff_to_feasible_n (
            p, idx, restrict_to_prior_stddev, inv_param);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH:
        diff = rhea_inversion_param_diff_to_feasible_yield (
            p, idx, restrict_to_prior_stddev, inv_param);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE:
        diff = rhea_inversion_param_diff_to_feasible_weak (
            p, idx, restrict_to_prior_stddev, inv_param);
        success = 1;
        break;
      default:
        break;
      }

      /* look up indices in a range */
      if (!success && rhea_inversion_param_idx_is_weak_label (idx)) {
        diff = rhea_inversion_param_diff_to_feasible_weak (
            p, idx, restrict_to_prior_stddev, inv_param);
        success = 1;
      }

      /* otherwise the parameter index is unknown */
      RHEA_CHECK_ABORT (success, "Parameter index is unknown.");

      /* project parameter value to feasible range */
      if (0.0 < fabs (diff)) {
        param[idx] = p + diff;
        param_modified = 1;
      }
    }
  }

  /* return whether parameters were modified */
  return param_modified;
}

int
rhea_inversion_param_restrict_step_length_to_feasible (
                                        ymir_vec_t *step_vec,
                                        ymir_vec_t *parameter_vec,
                                        const double restrict_to_prior_stddev,
                                        rhea_inversion_param_t *inv_param,
                                        double *step_length_new)
{
  const double       *param = parameter_vec->meshfree->e[0];
  double             *step = step_vec->meshfree->e[0];
  double              step_length = 1.0;
  int                 step_modified = 0;
  const int          *active = inv_param->active;
  int                 idx;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (step_vec, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (parameter_vec, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (step_vec, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (parameter_vec, inv_param));

  for (idx = 0; idx < inv_param->n_parameters; idx++) {
    const double        p = param[idx];
    const double        s = step[idx];
#ifdef RHEA_ENABLE_DEBUG
    double              p_diff;
#endif
    double              s_diff = 0.0;
    int                 success = 0;

    if (active[idx]) {
      /* look up single index */
      switch (idx) {
      case RHEA_INVERSION_PARAM_VISC_MIN:
      case RHEA_INVERSION_PARAM_VISC_MAX:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CLASS_FRACTURE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_THICKNESS_CONST_CLASS_FRACTURE:
#ifdef RHEA_ENABLE_DEBUG
        p_diff = _diff_to_feasible (
            p, idx, SC_1000_EPS, NAN, NAN /* !restrict to prior */, inv_param);
        RHEA_ASSERT (p_diff <= 0.0);
#endif
        s_diff = _diff_to_feasible (
            p + s, idx, SC_1000_EPS, NAN, restrict_to_prior_stddev, inv_param);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_SCALING:
      case RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_SCALING:
#ifdef RHEA_ENABLE_DEBUG
        p_diff = rhea_inversion_param_diff_to_feasible_scal (
            p, idx, NAN /* !restrict to prior */, inv_param);
        RHEA_ASSERT (p_diff <= 0.0);
#endif
        s_diff = rhea_inversion_param_diff_to_feasible_scal (
            p + s, idx, restrict_to_prior_stddev, inv_param);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_UPPER_MANTLE_ACTIVATION_ENERGY:
      case RHEA_INVERSION_PARAM_VISC_LOWER_MANTLE_ACTIVATION_ENERGY:
#ifdef RHEA_ENABLE_DEBUG
        p_diff = rhea_inversion_param_diff_to_feasible_arrh (
            p, idx, NAN /* !restrict to prior */, inv_param);
        RHEA_ASSERT (p_diff <= 0.0);
#endif
        s_diff = rhea_inversion_param_diff_to_feasible_arrh (
            p + s, idx, restrict_to_prior_stddev, inv_param);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_STRESS_EXPONENT:
#ifdef RHEA_ENABLE_DEBUG
        p_diff = rhea_inversion_param_diff_to_feasible_n (
            p, idx, NAN /* !restrict to prior */, inv_param);
        RHEA_ASSERT (p_diff <= 0.0);
#endif
        s_diff = rhea_inversion_param_diff_to_feasible_n (
            p + s, idx, restrict_to_prior_stddev, inv_param);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_VISC_YIELD_STRENGTH:
#ifdef RHEA_ENABLE_DEBUG
        p_diff = rhea_inversion_param_diff_to_feasible_yield (
            p, idx, NAN /* !restrict to prior */, inv_param);
        RHEA_ASSERT (p_diff <= 0.0);
#endif
        s_diff = rhea_inversion_param_diff_to_feasible_yield (
            p + s, idx, restrict_to_prior_stddev, inv_param);
        success = 1;
        break;
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_NONE:
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_SLAB:
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_RIDGE:
      case RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE:
#ifdef RHEA_ENABLE_DEBUG
        p_diff = rhea_inversion_param_diff_to_feasible_weak (
            p, idx, NAN /* !restrict to prior */, inv_param);
        RHEA_ASSERT (p_diff <= 0.0);
#endif
        s_diff = rhea_inversion_param_diff_to_feasible_weak (
            p + s, idx, restrict_to_prior_stddev, inv_param);
        success = 1;
        break;
      default:
        break;
      }

      /* look up indices in a range */
      if (!success && rhea_inversion_param_idx_is_weak_label (idx)) {
#ifdef RHEA_ENABLE_DEBUG
        p_diff = rhea_inversion_param_diff_to_feasible_weak (
            p, idx, NAN /* !restrict to prior */, inv_param);
        RHEA_ASSERT (p_diff <= 0.0);
#endif
        s_diff = rhea_inversion_param_diff_to_feasible_weak (
            p + s, idx, restrict_to_prior_stddev, inv_param);
        success = 1;
      }

      /* otherwise the parameter index is unknown */
      RHEA_CHECK_ABORT (success, "Parameter index is unknown.");

      /* reduce step length */
      if (0.0 < fabs (s_diff)) {
        const double        sl_new = (s + s_diff)/s;

        if (1.0e6*SC_EPS < sl_new) { /* if positive step length feasible */
          step_length = SC_MIN (sl_new, step_length);
        }
        else { /* otherwise positive step length is not possible */
          step[idx] = 0.0;
        }
        step_modified = 1;
      }
    }
  }

  /* reduce step */
  if (step_length < 1.0) {
    ymir_vec_scale (step_length, step_vec);
  }

  /* return new step length */
  if (step_length_new != NULL) {
    *step_length_new = step_length;
  }
  return step_modified;
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
rhea_inversion_param_set_initial_guess (ymir_vec_t *parameter_vec,
                                        rhea_inversion_param_t *inv_param,
                                        rhea_inversion_param_options_t *opt)
{
  const char         *file_path_txt = opt->initial_guess_file_path_txt;
  const double        perturb_stddev = opt->initial_guess_perturb_stddev;
  const double        shift_by_prior_stddev =
                        opt->initial_guess_shift_by_prior_stddev;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (
      __func__, "file_path=%s, perturb_stddev=%g, shift_by_prior_stddev=%g",
      file_path_txt, perturb_stddev, shift_by_prior_stddev);

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (parameter_vec, inv_param));

  /* get parameters */
  if (file_path_txt != NULL) { /* if get values from file */
    ymir_mesh_t        *ymir_mesh =
      rhea_stokes_problem_get_ymir_mesh (inv_param->stokes_problem);
    sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
    sc_dmatrix_t       *param_reduced;
    int                 success_n = 0;

    param_reduced = rhea_inversion_param_vec_reduced_new (parameter_vec,
                                                          inv_param);
    success_n = rhea_io_mpi_read_broadcast_double (
        param_reduced->e[0], param_reduced->m * param_reduced->n,
        NULL /* path bin */, file_path_txt, mpicomm);
    RHEA_ASSERT (success_n == inv_param->n_active);
    rhea_inversion_param_vec_reduced_copy (parameter_vec, param_reduced,
                                           inv_param);
    rhea_inversion_param_vec_reduced_destroy (param_reduced);
  }
  else { /* otherwise get values from model */
    rhea_inversion_param_pull_from_model (parameter_vec, inv_param);
  }

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

  /* restrict parameters to feasible values */
  rhea_inversion_param_restrict_to_feasible (
      parameter_vec, NAN /* !restrict to prior */, inv_param);

  /* print parameters */
  RHEA_GLOBAL_INFO ("========================================\n");
  RHEA_GLOBAL_INFOF ("%s\n", __func__);
  RHEA_GLOBAL_INFO ("----------------------------------------\n");
  rhea_inversion_param_vec_print (parameter_vec, inv_param);
  RHEA_GLOBAL_INFO ("========================================\n");

  /* check output */
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (parameter_vec, inv_param));

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static rhea_viscosity_param_derivative_t
rhea_inversion_param_get_derivative_type (
                                    const rhea_inversion_param_idx_t param_idx,
                                    rhea_weakzone_label_t *weak_label,
                                    rhea_weakzone_options_t *weak_options)
{
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

  /* look up indices related to parameters corresponding weak factor by label */
  if (rhea_inversion_param_idx_is_weak_label (param_idx)) {
    RHEA_ASSERT (0 <= rhea_inversion_param_idx_lookup_offset (param_idx));
    *weak_label = rhea_weakzone_lookup_label_from_index (
        param_idx - rhea_inversion_param_idx_lookup_offset (param_idx),
        weak_options);
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
  int                 idx;

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
      for (idx = 0; idx < inv_param->n_parameters; idx++) {
        if (active[idx]) {
          RHEA_GLOBAL_INFOF ("param# %3i: %g, %g\n",
                             idx, param[idx], model[idx]);
        }
      }
      rhea_inversion_param_vec_destroy (model_vals);
    }
    break;
  case RHEA_INVERSION_PARAM_VERBOSE_REAL_NONDIM_DIM:
    {
      ymir_vec_t     *model_vals = rhea_inversion_param_vec_new (inv_param);
      ymir_vec_t     *model_dim_vals = rhea_inversion_param_vec_new (inv_param);
      double         *model = model_vals->meshfree->e[0];
      double         *model_dim = model_dim_vals->meshfree->e[0];

      rhea_inversion_param_get_model_vals (model_vals, inv_param);
      ymir_vec_copy (model_vals, model_dim_vals);
      rhea_inversion_param_convert_model_vals_to_dim_model_vals (
          model_dim_vals, inv_param);
      for (idx = 0; idx < inv_param->n_parameters; idx++) {
        if (active[idx]) {
          RHEA_GLOBAL_INFOF ("param# %3i: %g, %g, %g\n",
                             idx, param[idx], model[idx], model_dim[idx]);
        }
      }
      rhea_inversion_param_vec_destroy (model_vals);
      rhea_inversion_param_vec_destroy (model_dim_vals);
    }
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
  RHEA_GLOBAL_INFO ("========================================\n");
  RHEA_GLOBAL_INFOF ("%s\n", __func__);
  RHEA_GLOBAL_INFO ("----------------------------------------\n");
  rhea_inversion_param_vec_print (prior_mean_vec, inv_param);
  RHEA_GLOBAL_INFO ("========================================\n");
}

static void
rhea_inversion_param_prior_get_stddev (ymir_vec_t *prior_stddev_vec,
                                       rhea_inversion_param_t *inv_param,
                                       rhea_inversion_param_options_t *opt)
{
  double             *p = prior_stddev_vec->meshfree->e[0];
  int                 idx;

  /* initialize all values to zero */
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
    int                 offset, count;

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

    /* set weak zone factors by class */
    idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_NONE;
    p[idx] = stddev_weak_factor_interior;
    idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_SLAB;
    p[idx] = stddev_weak_factor_interior;
    idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_RIDGE;
    p[idx] = stddev_weak_factor_interior;
    idx = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE;
    p[idx] = stddev_weak_factor_interior;

    /* set weak zone factors by label */
    offset = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_OFFSET;
    count  = RHEA_INVERSION_PARAM_WEAK_FACTOR_INTERIOR_LABEL_COUNT;
    for (idx = offset; idx < offset + count; idx++) {
      p[idx] = stddev_weak_factor_interior;
    }
  }

  /* convert std. deviations from units of model values to inversion param's */
  rhea_inversion_param_convert_model_vals_to_params (prior_stddev_vec,
                                                     inv_param);
  ymir_vec_fabs (prior_stddev_vec, prior_stddev_vec);

  /* print standard deviation */
  RHEA_GLOBAL_INFO ("========================================\n");
  RHEA_GLOBAL_INFOF ("%s\n", __func__);
  RHEA_GLOBAL_INFO ("----------------------------------------\n");
  rhea_inversion_param_vec_print (prior_stddev_vec, inv_param);
  RHEA_GLOBAL_INFO ("========================================\n");

  /* check output */
#ifdef RHEA_ENABLE_DEBUG
  for (idx = 0; idx < inv_param->n_parameters; idx++) {
    if (inv_param->active[idx]) {
      RHEA_ASSERT (isfinite (p[idx]));
      RHEA_ASSERT (0.0 < p[idx]);
    }
  }
#endif
}

/******************************************************************************
 * Parameter Related Computations
 *****************************************************************************/

/**
 * Calculates difference between parameters and prior mean:
 *   difference = param - param_prior_mean
 */
static void
rhea_inversion_param_prior_diff (ymir_vec_t *diff_vec,
                                 ymir_vec_t *parameter_vec,
                                 rhea_inversion_param_t *inv_param)
{
  if (NULL != parameter_vec) {
    RHEA_ASSERT (rhea_inversion_param_vec_check_type (parameter_vec,
                                                      inv_param));
    ymir_vec_copy (parameter_vec, diff_vec);
  }
  else {
    ymir_vec_set_zero (diff_vec);
  }
  ymir_vec_add (-1.0, inv_param->prior_mean, diff_vec);
}

void
rhea_inversion_param_prior_inv_cov (ymir_vec_t *inverse_covariance,
                                    rhea_inversion_param_t *inv_param)
{
#ifdef RHEA_ENABLE_DEBUG
  const double       *stddev = inv_param->prior_stddev->meshfree->e[0];
  const double       *icov = inverse_covariance->meshfree->e[0];
  int                 idx;

  /* check input */
  for (idx = 0; idx < inv_param->n_parameters; idx++) {
    if (inv_param->active[idx]) {
      RHEA_ASSERT (isfinite (stddev[idx]));
      RHEA_ASSERT (0.0 < stddev[idx]);
    }
  }
#endif

  /* compute inverse covariance of the prior */
  ymir_vec_copy (inv_param->prior_stddev, inverse_covariance);
  ymir_vec_multiply_in (inv_param->prior_stddev, inverse_covariance);
  ymir_vec_reciprocal (inverse_covariance);

  /* check output */
#ifdef RHEA_ENABLE_DEBUG
  for (idx = 0; idx < inv_param->n_parameters; idx++) {
    if (inv_param->active[idx]) {
      RHEA_ASSERT (isfinite (icov[idx]));
      RHEA_ASSERT (0.0 < icov[idx]);
    }
  }
#endif
}

double
rhea_inversion_param_prior (ymir_vec_t *parameter_vec,
                            rhea_inversion_param_t *inv_param)
{
  ymir_vec_t         *prior_diff = rhea_inversion_param_vec_new (inv_param);
  ymir_vec_t         *prior_icov = rhea_inversion_param_vec_new (inv_param);
  const int          *active = inv_param->active;
  const double       *diff = prior_diff->meshfree->e[0];
  const double       *icov = prior_icov->meshfree->e[0];
  int                 idx;
  double              prior_val = 0.0;

  /* calculate difference between parameters and prior mean */
  rhea_inversion_param_prior_diff (prior_diff, parameter_vec, inv_param);

  /* calculate inverse prior covariance */
  rhea_inversion_param_prior_inv_cov (prior_icov, inv_param);

  /* calculate inner product */
  for (idx = 0; idx < inv_param->n_parameters; idx++) {
    if (active[idx]) {
      RHEA_ASSERT (isfinite (diff[idx]));
      prior_val += diff[idx] * icov[idx] * diff[idx];
    }
  }

  /* destroy */
  rhea_inversion_param_vec_destroy (prior_diff);
  rhea_inversion_param_vec_destroy (prior_icov);

  /* return squared norm of prior */
  return 0.5*prior_val;
}

static double
_param_derivative_from_adjoint (ymir_vec_t *forward_vel,
                                ymir_vec_t *adjoint_vel,
                                ymir_stress_op_t *stress_op_param_derivative,
                                ymir_vec_t *tmp_vel)
{
  /* apply viscous stress operator to forward velocity */
  ymir_stress_pc_apply_stress_op (forward_vel, tmp_vel,
                                  stress_op_param_derivative,
                                  0 /* !linearized */, 1 /* dirty */);

  /* compute inner product with adjoint velocity */
  return ymir_vec_innerprod (tmp_vel, adjoint_vel);
}

void
rhea_inversion_param_compute_gradient (
                              ymir_vec_t *gradient_vec,
                              ymir_vec_t *parameter_vec,
                              ymir_vec_t *forward_vel_press,
                              ymir_vec_t *adjoint_vel_press,
                              const double prior_weight,
                              rhea_inversion_param_t *inv_param,
                              ymir_vec_t *gradient_adjoint_comp_vec,
                              ymir_vec_t *gradient_prior_comp_vec,
                              ymir_vec_t *gradient_callback_comp_vec,
                              rhea_inversion_param_derivative_fn_t derivative_fn,
                              void *derivative_fn_data)
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
  ymir_vec_t         *forward_vel, *adjoint_vel, *tmp_vel;
  ymir_vec_t         *viscosity, *marker;
  ymir_vel_dir_t     *vel_dir;
  ymir_stress_op_t   *stress_op;

  ymir_vec_t         *derivative;
  const int           dim_exists = (inv_param->dimensional_scaling != NULL);
  const double       *dim =
    (dim_exists ? inv_param->dimensional_scaling->meshfree->e[0] : NULL);
  const int          *active = inv_param->active;
  double             *grad = gradient_vec->meshfree->e[0];
  double             *grad_adj=NULL, *grad_prior=NULL, *grad_callback=NULL;
  int                 idx;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (__func__, "prior_weight=%g", prior_weight);

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (gradient_vec, inv_param));
  RHEA_ASSERT (parameter_vec == NULL ||
               rhea_inversion_param_vec_check_type (parameter_vec, inv_param));
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (forward_vel_press));
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (adjoint_vel_press));

  /* retrieve forward and adjoint velocities */
  forward_vel = rhea_velocity_new (ymir_mesh);
  adjoint_vel = rhea_velocity_new (ymir_mesh);
  tmp_vel     = rhea_velocity_new (ymir_mesh);
  rhea_velocity_pressure_copy_components (forward_vel, NULL, forward_vel_press,
                                          press_elem);
  rhea_velocity_pressure_copy_components (adjoint_vel, NULL, adjoint_vel_press,
                                          press_elem);

  /* compute viscosity and related fields */
  viscosity = rhea_viscosity_new (ymir_mesh);
  marker    = rhea_viscosity_new (ymir_mesh);
  rhea_viscosity_compute (
      /* out: */ viscosity, NULL, marker,
      /* in:  */ temperature, weakzone, forward_vel, inv_param->visc_options);

  /* init derivative to use as viscous stress coefficient */
  derivative = rhea_viscosity_new (ymir_mesh);
  ymir_vec_set_value (derivative, 1.0);

  /* create new viscous stress operator */
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
  if (gradient_callback_comp_vec != NULL) {
    ymir_vec_set_zero (gradient_callback_comp_vec);
    grad_callback = gradient_callback_comp_vec->meshfree->e[0];
  }
  for (idx = 0; idx < inv_param->n_parameters; idx++) {
    if (active[idx]) {
      rhea_viscosity_param_derivative_t deriv_type;
      rhea_weakzone_label_t             weak_label;
      double                            ga, gc;

      RHEA_GLOBAL_VERBOSEF_FN_TAG (__func__, "param_idx=%i", idx);

      /* compute parameter derivative of viscosity */
      deriv_type = rhea_inversion_param_get_derivative_type (
          (rhea_inversion_param_idx_t) idx, &weak_label,
          inv_param->weak_options);
      if (RHEA_WEAKZONE_LABEL_UNKNOWN == weak_label) {
        rhea_viscosity_param_derivative (
            derivative, deriv_type, viscosity, marker, temperature, weakzone,
            forward_vel, inv_param->visc_options);
      }
      else {
        rhea_viscosity_param_derivative_weakzone (
            derivative, deriv_type, weak_label, viscosity, marker, weakzone,
            inv_param->weak_options, inv_param->visc_options);
      }

      /* apply dimensional scaling */
      if (dim_exists) {
        RHEA_ASSERT (0.0 < dim[idx]);
        ymir_vec_scale (dim[idx], derivative);
      }

      /* transform to viscous stress coefficient */
      ymir_vec_scale (2.0, derivative);

      /* set derivative as the coefficient of the viscous stress operator */
      ymir_stress_op_set_coeff_scal (stress_op, derivative);

      /* compute derivative from forward and adjoint Stokes solutions */
      ga = _param_derivative_from_adjoint (forward_vel, adjoint_vel, stress_op,
                                           tmp_vel);

      /* compute gradient from callback */
      if (NULL != derivative_fn) {
        gc = derivative_fn (forward_vel_press, adjoint_vel_press,
                            idx, (int) deriv_type, (int) weak_label,
                            derivative, stress_op, derivative_fn_data);
      }
      else {
        gc = 0.0;
      }

      /* set gradient entry */
      grad[idx] = ga + gc;
      if (grad_adj != NULL) {
        grad_adj[idx] = ga;
      }
      if (grad_callback != NULL) {
        grad_callback[idx] = gc;
      }
    }
  }
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (gradient_vec, inv_param));

  /* compute & add prior term */
  if (isfinite (prior_weight) && 0.0 < prior_weight) {
    ymir_vec_t         *prior_diff = rhea_inversion_param_vec_new (inv_param);
    ymir_vec_t         *prior_icov = rhea_inversion_param_vec_new (inv_param);
    const double       *diff = prior_diff->meshfree->e[0];
    const double       *icov = prior_icov->meshfree->e[0];

    /* calculate difference between parameters and prior mean */
    rhea_inversion_param_prior_diff (prior_diff, parameter_vec, inv_param);

    /* calculate inverse prior covariance */
    rhea_inversion_param_prior_inv_cov (prior_icov, inv_param);

    /* add prior term to gradient */
    for (idx = 0; idx < inv_param->n_parameters; idx++) {
      if (active[idx]) {
        double              g;

        RHEA_ASSERT (isfinite (diff[idx]));
        g = prior_weight * icov[idx] * diff[idx];
        grad[idx] += g;
        if (grad_prior != NULL) {
          grad_prior[idx] = g;
        }
      }
    }

    /* destroy */
    rhea_inversion_param_vec_destroy (prior_diff);
    rhea_inversion_param_vec_destroy (prior_icov);
  }
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (gradient_vec, inv_param));

  /* destroy */
  ymir_stress_op_destroy (stress_op);
  ymir_vel_dir_destroy (vel_dir);
  rhea_viscosity_destroy (derivative);
  rhea_viscosity_destroy (viscosity);
  rhea_viscosity_destroy (marker);
  rhea_velocity_destroy (forward_vel);
  rhea_velocity_destroy (adjoint_vel);
  rhea_velocity_destroy (tmp_vel);

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

double
rhea_inversion_param_compute_gradient_norm (ymir_vec_t *gradient_vec,
                                            rhea_inversion_param_t *inv_param)
{
  const int           weight_icov =
                        RHEA_INVERSION_PARAM_DEFAULT_GRADIENT_NORM_WEIGHT_ICOV;
  ymir_vec_t         *prior_icov = NULL;
  double              norm;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (gradient_vec, inv_param));

  /* calculate inverse prior covariance */
  if (weight_icov) {
    prior_icov = rhea_inversion_param_vec_new (inv_param);
    rhea_inversion_param_prior_inv_cov (prior_icov, inv_param);
  }

  /* compute inner product */
  norm = sqrt (rhea_inversion_param_vec_reduced_inner_product (
      gradient_vec, gradient_vec, prior_icov, inv_param,
      1 /* normalize wrt. size */));

  /* destroy */
  if (weight_icov) {
    rhea_inversion_param_vec_destroy (prior_icov);
  }

  /* return l2-norm of active entries of the gradient vector */
  return norm;
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
  ymir_vec_t         *viscosity, *marker;
  ymir_vel_dir_t     *vel_dir;
  ymir_stress_op_t   *stress_op;

  ymir_vec_t         *derivative;
  const int           dim_exists = (inv_param->dimensional_scaling != NULL);
  const double       *dim =
    (dim_exists ? inv_param->dimensional_scaling->meshfree->e[0] : NULL);
  const int          *active = inv_param->active;
  double             *grad_dir = gradient_direction->meshfree->e[0];
  int                 idx;

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
  marker = rhea_viscosity_new (ymir_mesh);
  rhea_viscosity_compute (
      /* out: */ viscosity, NULL, marker,
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
  for (idx = 0; idx < inv_param->n_parameters; idx++) {
    if (active[idx] && 0.0 < fabs (grad_dir[idx])) {
      rhea_viscosity_param_derivative_t deriv_type;
      rhea_weakzone_label_t             weak_label;

      /* compute parameter derivative of viscosity */
      deriv_type = rhea_inversion_param_get_derivative_type (
          (rhea_inversion_param_idx_t) idx, &weak_label,
          inv_param->weak_options);
      if (RHEA_WEAKZONE_LABEL_UNKNOWN == weak_label) {
        rhea_viscosity_param_derivative (
            derivative, deriv_type, viscosity, marker, temperature, weakzone,
            forward_vel, inv_param->visc_options);
      }
      else {
        rhea_viscosity_param_derivative_weakzone (
            derivative, deriv_type, weak_label, viscosity, marker, weakzone,
            inv_param->weak_options, inv_param->visc_options);
      }

      /* apply dimensional scaling */
      if (dim_exists) {
        RHEA_ASSERT (0.0 < dim[idx]);
        ymir_vec_scale (dim[idx], derivative);
      }

      /* transform to viscous stress coefficient */
      ymir_vec_scale (2.0, derivative);

      /* set derivative as the coefficient of the viscous stress operator */
      ymir_stress_op_set_coeff_scal (stress_op, derivative);

      /* apply viscous stress operator to forward velocity */
      ymir_stress_pc_apply_stress_op (forward_vel, op_out_vel, stress_op,
                                      0 /* !linearized */, 1 /* dirty */);

      /* scale and add gradient to right-hand side */
      ymir_vec_add (-grad_dir[idx], op_out_vel, rhs_vel_mass);
    }
  }
  RHEA_ASSERT (rhea_velocity_is_valid (rhs_vel_mass));

  /* destroy */
  ymir_stress_op_destroy (stress_op);
  ymir_vel_dir_destroy (vel_dir);
  rhea_viscosity_destroy (derivative);
  rhea_viscosity_destroy (viscosity);
  rhea_viscosity_destroy (marker);
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
  const int           use_data_term = (forward_vel_press != NULL &&
                                       incr_adjoint_vel_press != NULL);
  const int           use_prior_term = (isfinite (prior_weight) &&
                                        0.0 < prior_weight);
  const int           use_first_order_approx = (adjoint_vel_press == NULL);

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (param_vec_out, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (param_vec_in, inv_param));

  /*
   * First-Order Derivative Terms (aka. Gauss-Newton Hessian)
   */

  /* compute data term */
  if (use_data_term) {
    RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (forward_vel_press));
    RHEA_ASSERT (
        rhea_velocity_pressure_check_vec_type (incr_adjoint_vel_press));

    /* "hijack" the gradient function using incr. adj. vel. */
    rhea_inversion_param_compute_gradient (
        param_vec_out, param_vec_in, forward_vel_press, incr_adjoint_vel_press,
        NAN /* no prior term */, inv_param, NULL, NULL, NULL,
        NULL, NULL); //TODO add derivatives from callback function
  }
  else {
    ymir_vec_set_zero (param_vec_out);
  }

  /* compute & add prior term */
  if (use_prior_term) {
    ymir_vec_t         *prior_term = rhea_inversion_param_vec_new (inv_param);
    ymir_vec_t         *prior_icov = rhea_inversion_param_vec_new (inv_param);
    const int          *active = inv_param->active;
    const double       *prterm = prior_term->meshfree->e[0];
    const double       *icov = prior_icov->meshfree->e[0];
    double             *param_out = param_vec_out->meshfree->e[0];
    int                 idx;

    /* copy input */
    ymir_vec_copy (param_vec_in, prior_term);

    /* calculate inverse prior covariance */
    rhea_inversion_param_prior_inv_cov (prior_icov, inv_param);

    /* add prior term to Hessian apply output */
    for (idx = 0; idx < inv_param->n_parameters; idx++) {
      if (active[idx]) {
        RHEA_ASSERT (isfinite (prterm[idx]));
        param_out[idx] += prior_weight * icov[idx] * prterm[idx];
      }
    }

    /* destroy */
    rhea_inversion_param_vec_destroy (prior_term);
    rhea_inversion_param_vec_destroy (prior_icov);
  }

  /*
   * Second-Order Derivative Terms
   */

  /* compute & add data term */
  if (use_data_term && !use_first_order_approx) { /* if full Hessian */
    RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (adjoint_vel_press));
    RHEA_ASSERT (
        rhea_velocity_pressure_check_vec_type (incr_forward_vel_press));

    RHEA_ABORT_NOT_REACHED (); //TODO
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

  /* set random values; add one to generate a (multiplicative) perturbation */
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
  int                 idx;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (vec, inv_param));

  /* inspect processor-local values */
  for (idx = 0; idx < n_parameters; idx++) {
    if (active[idx] && !isfinite (v[idx])) {
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
    for (idx = 0; idx < n_parameters; idx++) {
      const double        rel_err = fabs (v[idx] - v_prev[idx]) / fabs (v[idx]);

      if (active[idx] && SC_1000_EPS < rel_err) {
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
  int                 idx;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (parameter_vec, inv_param));

  /* exit if nothing to do */
  if (!inv_param->n_active) {
    return;
  }

  /* print each parameter */
  for (idx = 0; idx < inv_param->n_parameters; idx++) {
    if (active[idx]) {
      RHEA_GLOBAL_INFOF ("param# %3i: %g\n", idx, param[idx]);
    }
  }
}

sc_dmatrix_t *
rhea_inversion_param_vec_reduced_new (ymir_vec_t *vec,
                                      rhea_inversion_param_t *inv_param)
{
  const int          *active = inv_param->active;
  sc_dmatrix_t       *vec_reduced;

  /* check input */
  RHEA_ASSERT (NULL == vec ||
               rhea_inversion_param_vec_check_type (vec, inv_param));

  /* create reduced parameter vector */
  vec_reduced = sc_dmatrix_new (inv_param->n_active, 1);

  /* fill entries */
  if (NULL != vec) {
    const double       *v = vec->meshfree->e[0];
    double             *r = vec_reduced->e[0];
    int                 idx, row;

    row = 0;
    for (idx = 0; idx < inv_param->n_parameters; idx++) {
      if (active[idx]) {
        RHEA_ASSERT (row < vec_reduced->m * vec_reduced->n);
        r[row] = v[idx];
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
  const int          *active = inv_param->active;
  const double       *r = vec_reduced->e[0];
  double             *v = vec->meshfree->e[0];
  int                 idx, row;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (vec, inv_param));
  RHEA_ASSERT (vec_reduced->m == inv_param->n_active);
  RHEA_ASSERT (vec_reduced->n == 1);

  /* fill entries */
  row = 0;
  for (idx = 0; idx < inv_param->n_parameters; idx++) {
    if (active[idx]) {
      RHEA_ASSERT (row < vec_reduced->m * vec_reduced->n);
      v[idx] = r[row];
      row++;
    }
  }
}

void
rhea_inversion_param_vec_reduced_apply_matrix (
                                            ymir_vec_t *out,
                                            sc_dmatrix_t *matrix_reduced,
                                            ymir_vec_t *in,
                                            rhea_inversion_param_t *inv_param)
{
  sc_dmatrix_t       *in_reduced, *out_reduced;

  /* create reduced parameter vectors */
  in_reduced = rhea_inversion_param_vec_reduced_new (in, inv_param);
  out_reduced = rhea_inversion_param_vec_reduced_new (out, inv_param);

  /* run matrix-vector multiply; copy entries to output */
  sc_dmatrix_vector (SC_NO_TRANS, SC_NO_TRANS, SC_NO_TRANS,
                     1.0 /* factor for in vec */, matrix_reduced, in_reduced,
                     0.0 /* factor for out vec */, out_reduced);
  ymir_vec_set_zero (out);
  rhea_inversion_param_vec_reduced_copy (out, out_reduced, inv_param);

  /* destroy */
  rhea_inversion_param_vec_reduced_destroy (in_reduced);
  rhea_inversion_param_vec_reduced_destroy (out_reduced);
}

int
rhea_inversion_param_vec_reduced_solve_matrix (
                                            ymir_vec_t *sol,
                                            sc_dmatrix_t *matrix_reduced,
                                            ymir_vec_t *rhs,
                                            rhea_inversion_param_t *inv_param,
                                            int *iter_count)
{
  sc_dmatrix_t       *sol_reduced, *rhs_reduced;

  /* create work variables */
  sol_reduced = rhea_inversion_param_vec_reduced_new (sol, inv_param);
  rhs_reduced = rhea_inversion_param_vec_reduced_new (rhs, inv_param);

  /* run direct solver; copy entries to output */
  sc_dmatrix_ldivide (SC_NO_TRANS, matrix_reduced, rhs_reduced, sol_reduced);
  ymir_vec_set_zero (sol);
  rhea_inversion_param_vec_reduced_copy (sol, sol_reduced, inv_param);

  /* destroy */
  rhea_inversion_param_vec_reduced_destroy (sol_reduced);
  rhea_inversion_param_vec_reduced_destroy (rhs_reduced);

  /* return iteraton count and "stopping" reason */
  if (iter_count != NULL) {
    *iter_count = 1;
  }
  return 1;
}

double
rhea_inversion_param_vec_reduced_inner_product (
                                            ymir_vec_t *vecL,
                                            ymir_vec_t *vecR,
                                            ymir_vec_t *weight,
                                            rhea_inversion_param_t *inv_param,
                                            const int normalize_wrt_size)
{
  const double       *vecL_data = vecL->meshfree->e[0];
  const double       *vecR_data = vecR->meshfree->e[0];
  const double       *weight_data =
                        (weight != NULL ? weight->meshfree->e[0] : NULL);
  const int          *active = inv_param->active;
  int                 idx;
  double              sum_of_products = 0.0;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (vecL, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (vecR, inv_param));
  RHEA_ASSERT (NULL == weight ||
               rhea_inversion_param_vec_check_type (weight, inv_param));

  /* sum up elementwise products of vector entries that are active */
  if (NULL != weight) { /* if weight exists */
    for (idx = 0; idx < inv_param->n_parameters; idx++) {
      if (active[idx]) {
        RHEA_ASSERT (isfinite (vecL_data[idx]));
        RHEA_ASSERT (isfinite (vecR_data[idx]));
        RHEA_ASSERT (isfinite (weight_data[idx]));
        sum_of_products += vecL_data[idx] * weight_data[idx] * vecR_data[idx];
      }
    }
  }
  else {
    for (idx = 0; idx < inv_param->n_parameters; idx++) {
      if (active[idx]) {
        RHEA_ASSERT (isfinite (vecL_data[idx]));
        RHEA_ASSERT (isfinite (vecR_data[idx]));
        sum_of_products += vecL_data[idx] * vecR_data[idx];
      }
    }
  }

  /* return scaled inner product */
  if (normalize_wrt_size) {
    return sum_of_products / (double) inv_param->n_active;
  }
  else {
    return sum_of_products;
  }
}

double
rhea_inversion_param_vec_reduced_ip (sc_dmatrix_t *vecL_reduced,
                                     sc_dmatrix_t *vecR_reduced,
                                     sc_dmatrix_t *weight_reduced,
                                     const int normalize_wrt_size)
{
  const double       *vecL_data = vecL_reduced->e[0];
  const double       *vecR_data = vecR_reduced->e[0];
  const double       *weight_data =
                        (weight_reduced != NULL ? weight_reduced->e[0] : NULL);
  const int           size = vecL_reduced->m * vecL_reduced->n;
  int                 idx;
  double              sum_of_products = 0.0;

  /* check input */
  RHEA_ASSERT (vecL_reduced->m * vecL_reduced->n ==
               vecR_reduced->m * vecR_reduced->n);
  RHEA_ASSERT (NULL == weight_reduced ||
               size == weight_reduced->m * weight_reduced->n);

  /* sum up elementwise products of vector entries */
  if (NULL != weight_reduced) { /* if weight exists */
    for (idx = 0; idx < size; idx++) {
      RHEA_ASSERT (isfinite (vecL_data[idx]));
      RHEA_ASSERT (isfinite (vecR_data[idx]));
      RHEA_ASSERT (isfinite (weight_data[idx]));
      sum_of_products += vecL_data[idx] * weight_data[idx] * vecR_data[idx];
    }
  }
  else {
    for (idx = 0; idx < size; idx++) {
      RHEA_ASSERT (isfinite (vecL_data[idx]));
      RHEA_ASSERT (isfinite (vecR_data[idx]));
      sum_of_products += vecL_data[idx] * vecR_data[idx];
    }
  }

  /* return scaled inner product */
  if (normalize_wrt_size) {
    return sum_of_products / (double) size;
  }
  else {
    return sum_of_products;
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
