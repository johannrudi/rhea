/*
 */

#include <rhea_stokes_problem_amr.h>
#include <rhea_base.h>
#include <rhea_amr.h>
#include <rhea_temperature.h>
#include <rhea_weakzone.h>
#include <ymir_interp_vec.h>
#include <ymir_hmg_intergrid_h.h>
#include <ymir_mass_vec.h>
#include <ymir_pressure_vec.h>
#include <ymir_vec_optimized.h>
#include <ymir_monitor.h>
#include <mangll_fields.h>

#if (0 < RHEA_STOKES_PROBLEM_AMR_VERBOSE_VTK)
# include <ymir_vtk.h>
#endif

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_STOKES_PROBLEM_AMR_DEFAULT_N_CYCLES (P4EST_MAXLEVEL)
#define RHEA_STOKES_PROBLEM_AMR_DEFAULT_FLAGGED_ELEMS_THRESH_BEGIN (NAN)
#define RHEA_STOKES_PROBLEM_AMR_DEFAULT_FLAGGED_ELEMS_THRESH_CYCLE (NAN)
#define RHEA_STOKES_PROBLEM_AMR_DEFAULT_INIT_TYPE_NAME "NONE"
#define RHEA_STOKES_PROBLEM_AMR_DEFAULT_INIT_TOL_MIN (NAN)
#define RHEA_STOKES_PROBLEM_AMR_DEFAULT_INIT_TOL_MAX (NAN)
#define RHEA_STOKES_PROBLEM_AMR_DEFAULT_NONLINEAR_TYPE_NAME "NONE"
#define RHEA_STOKES_PROBLEM_AMR_DEFAULT_NONLINEAR_TOL_MIN (NAN)
#define RHEA_STOKES_PROBLEM_AMR_DEFAULT_NONLINEAR_TOL_MAX (NAN)
#define RHEA_STOKES_PROBLEM_AMR_DEFAULT_NONLINEAR_ITER_FIRST (0)
#define RHEA_STOKES_PROBLEM_AMR_DEFAULT_NONLINEAR_ITER_LAST (-1)

/* initialize options */
int                 rhea_stokes_problem_amr_n_cycles =
  RHEA_STOKES_PROBLEM_AMR_DEFAULT_N_CYCLES;
double              rhea_stokes_problem_amr_flagged_elements_thresh_begin =
  RHEA_STOKES_PROBLEM_AMR_DEFAULT_FLAGGED_ELEMS_THRESH_BEGIN;
double              rhea_stokes_problem_amr_flagged_elements_thresh_cycle =
  RHEA_STOKES_PROBLEM_AMR_DEFAULT_FLAGGED_ELEMS_THRESH_CYCLE;
char               *rhea_stokes_problem_amr_init_type_name =
  RHEA_STOKES_PROBLEM_AMR_DEFAULT_INIT_TYPE_NAME;
double              rhea_stokes_problem_amr_init_tol_min =
  RHEA_STOKES_PROBLEM_AMR_DEFAULT_INIT_TOL_MIN;
double              rhea_stokes_problem_amr_init_tol_max =
  RHEA_STOKES_PROBLEM_AMR_DEFAULT_INIT_TOL_MAX;
char               *rhea_stokes_problem_amr_nonlinear_type_name =
  RHEA_STOKES_PROBLEM_AMR_DEFAULT_NONLINEAR_TYPE_NAME;
double              rhea_stokes_problem_amr_nonlinear_tol_min =
  RHEA_STOKES_PROBLEM_AMR_DEFAULT_NONLINEAR_TOL_MIN;
double              rhea_stokes_problem_amr_nonlinear_tol_max =
  RHEA_STOKES_PROBLEM_AMR_DEFAULT_NONLINEAR_TOL_MAX;
int                 rhea_stokes_problem_amr_nonlinear_iter_first =
  RHEA_STOKES_PROBLEM_AMR_DEFAULT_NONLINEAR_ITER_FIRST;
int                 rhea_stokes_problem_amr_nonlinear_iter_last =
  RHEA_STOKES_PROBLEM_AMR_DEFAULT_NONLINEAR_ITER_LAST;

void
rhea_stokes_problem_amr_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "AMR";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_I, "num-cycles", '\0',
    &(rhea_stokes_problem_amr_n_cycles),
    RHEA_STOKES_PROBLEM_AMR_DEFAULT_N_CYCLES,
    "#Cycles to perform AMR (each cycle coarsens/refines by mostly one level)",
  YMIR_OPTIONS_D, "flagged-elements-threshold-begin", '\0',
    &(rhea_stokes_problem_amr_flagged_elements_thresh_begin),
    RHEA_STOKES_PROBLEM_AMR_DEFAULT_FLAGGED_ELEMS_THRESH_BEGIN,
    "Threshold for relative min #elements flagged to begin AMR cycle(s)",
  YMIR_OPTIONS_D, "flagged-elements-threshold-cycle", '\0',
    &(rhea_stokes_problem_amr_flagged_elements_thresh_cycle),
    RHEA_STOKES_PROBLEM_AMR_DEFAULT_FLAGGED_ELEMS_THRESH_CYCLE,
    "Threshold for relative min #elements flagged to continue AMR cycles "
    "(if not provided it is assumed to be same value as at the beginning)",

  YMIR_OPTIONS_S, "init-amr-name", '\0',
    &(rhea_stokes_problem_amr_init_type_name),
    RHEA_STOKES_PROBLEM_AMR_DEFAULT_INIT_TYPE_NAME,
    "AMR for initial mesh: name of coarsening/refinemen indicator",
  YMIR_OPTIONS_D, "init-amr-tol-min", '\0',
    &(rhea_stokes_problem_amr_init_tol_min),
    RHEA_STOKES_PROBLEM_AMR_DEFAULT_INIT_TOL_MIN,
    "AMR for initial mesh: tolerance below which the mesh is coarsened",
  YMIR_OPTIONS_D, "init-amr-tol-max", '\0',
    &(rhea_stokes_problem_amr_init_tol_max),
    RHEA_STOKES_PROBLEM_AMR_DEFAULT_INIT_TOL_MAX,
    "AMR for initial mesh: tolerance above which the mesh is refined",

  YMIR_OPTIONS_S, "nonlinear-amr-name", '\0',
    &(rhea_stokes_problem_amr_nonlinear_type_name),
    RHEA_STOKES_PROBLEM_AMR_DEFAULT_NONLINEAR_TYPE_NAME,
    "AMR for nl. solver/grid cont.: name of coarsening/refinemen indicator",
  YMIR_OPTIONS_D, "nonlinear-amr-tol-min", '\0',
    &(rhea_stokes_problem_amr_nonlinear_tol_min),
    RHEA_STOKES_PROBLEM_AMR_DEFAULT_NONLINEAR_TOL_MIN,
    "AMR for nl. solver/grid cont.: tol. below which the mesh is coarsened",
  YMIR_OPTIONS_D, "nonlinear-amr-tol-max", '\0',
    &(rhea_stokes_problem_amr_nonlinear_tol_max),
    RHEA_STOKES_PROBLEM_AMR_DEFAULT_NONLINEAR_TOL_MAX,
    "AMR for nl. solver/grid cont.: tol. above which the mesh is refined",
  YMIR_OPTIONS_I, "nonlinear-amr-iter-first", '\0',
    &(rhea_stokes_problem_amr_nonlinear_iter_first),
    RHEA_STOKES_PROBLEM_AMR_DEFAULT_NONLINEAR_ITER_FIRST,
    "AMR for nl. solver/grid cont.: start AMR at this nonlinear iteration",
  YMIR_OPTIONS_I, "nonlinear-amr-iter-last", '\0',
    &(rhea_stokes_problem_amr_nonlinear_iter_last),
    RHEA_STOKES_PROBLEM_AMR_DEFAULT_NONLINEAR_ITER_LAST,
    "AMR for nl. solver/grid cont.: max # of nonlinear iterations with AMR",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);
}

/******************************************************************************
 * AMR Data
 *****************************************************************************/

typedef struct rhea_stokes_problem_amr_data
{
  /* mesh */
  ymir_mesh_t            *ymir_mesh;
  ymir_pressure_elem_t   *press_elem;

  /* Stokes problem */
  rhea_stokes_problem_t  *stokes_problem;

  /* AMR iterations */
  int                 first_amr;
  mangll_t           *mangll_original;
  mangll_t           *mangll_adapted;
  mangll_t           *mangll_partitioned;

  /* buffers for projection & partitioning */
  sc_dmatrix_t       *temperature_original;
  sc_dmatrix_t       *temperature_adapted;

  sc_dmatrix_t       *velocity_original;
  sc_dmatrix_t       *velocity_adapted;

  sc_dmatrix_t       *pressure_original;
  sc_dmatrix_t       *pressure_adapted;

  /* options (not owned) */
  rhea_discretization_options_t  *discr_options;

  /* thresholds for coarsening/refinement */
  double              tol_min;
  double              tol_max;

  /* boolean to enable merging of coarsening/refinement flags */
  int                 merge_mode;
}
rhea_stokes_problem_amr_data_t;

static void
rhea_stokes_problem_amr_data_clear_fields (
                                      rhea_stokes_problem_amr_data_t *amr_data)
{
  rhea_stokes_problem_t *stokes_problem = amr_data->stokes_problem;
  ymir_vec_t         *temperature, *velocity_pressure;

  temperature = rhea_stokes_problem_get_temperature (stokes_problem);
  if (temperature != NULL) {
    rhea_temperature_destroy (temperature);
  }
  rhea_stokes_problem_remove_temperature (stokes_problem);

  velocity_pressure = rhea_stokes_problem_get_velocity_pressure (
      stokes_problem);
  if (velocity_pressure != NULL) {
    rhea_velocity_pressure_destroy (velocity_pressure);
  }
  rhea_stokes_problem_remove_velocity_pressure (stokes_problem);
}

static void
rhea_stokes_problem_amr_data_clear_buffer (
                                      rhea_stokes_problem_amr_data_t *amr_data)
{
  if (amr_data->temperature_original != NULL) {
    sc_dmatrix_destroy (amr_data->temperature_original);
  }
  if (amr_data->temperature_adapted != NULL) {
    sc_dmatrix_destroy (amr_data->temperature_adapted);
  }

  if (amr_data->velocity_original != NULL) {
    sc_dmatrix_destroy (amr_data->velocity_original);
  }
  if (amr_data->velocity_adapted != NULL) {
    sc_dmatrix_destroy (amr_data->velocity_adapted);
  }

  if (amr_data->pressure_original != NULL) {
    sc_dmatrix_destroy (amr_data->pressure_original);
  }
  if (amr_data->pressure_adapted != NULL) {
    sc_dmatrix_destroy (amr_data->pressure_adapted);
  }

  amr_data->temperature_original = NULL;
  amr_data->temperature_adapted = NULL;
  amr_data->velocity_original = NULL;
  amr_data->velocity_adapted = NULL;
  amr_data->pressure_original = NULL;
  amr_data->pressure_adapted = NULL;
}

rhea_stokes_problem_amr_data_t *
rhea_stokes_problem_amr_data_new (rhea_stokes_problem_t *stokes_problem,
                                  rhea_discretization_options_t *discr_options)
{
  rhea_stokes_problem_amr_data_t *amr_data;

  /* create */
  amr_data = RHEA_ALLOC (rhea_stokes_problem_amr_data_t, 1);

  /* init original mesh */
  amr_data->ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  amr_data->press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);

  /* init Stokes problem */
  amr_data->stokes_problem = stokes_problem;

  /* init AMR meshes */
  amr_data->first_amr = -1;
  amr_data->mangll_original = amr_data->ymir_mesh->ma;
  amr_data->mangll_adapted = NULL;
  amr_data->mangll_partitioned = NULL;

  /* init buffers */
  amr_data->temperature_original = NULL;
  amr_data->temperature_adapted = NULL;
  amr_data->velocity_original = NULL;
  amr_data->velocity_adapted = NULL;
  amr_data->pressure_original = NULL;
  amr_data->pressure_adapted = NULL;

  /* init options */
  amr_data->discr_options = discr_options;
  amr_data->tol_min = NAN;
  amr_data->tol_max = NAN;
  amr_data->merge_mode = 0;

  return amr_data;
}

void
rhea_stokes_problem_amr_data_destroy (rhea_stokes_problem_amr_data_t *amr_data)
{
  rhea_stokes_problem_amr_data_clear_buffer (amr_data);
  RHEA_FREE (amr_data);
}

/******************************************************************************
 * Flagging for Coarsening/Refinement
 *
 * To flag an element for coarsening or refinement we balance two terms:
 *
 *   indicator (element) * relative_size (element),
 *
 * where the relative element size is always
 *
 *   relative_size (element) = (element size / domain size).
 *
 * For Peclet-type indicators the indicator term takes on the form
 *
 *   indicator_peclet (element) = || |grad fn| / fn ||_element,
 *
 * where `fn` is a function with positive (scalar) values and
 * the norm is the L2-norm taken over an element.
 *****************************************************************************/

static double
rhea_stokes_problem_amr_norm_L2_elem (sc_dmatrix_t *ind_mat,
                                      sc_dmatrix_t *ind_mass_mat,
                                      sc_dmatrix_t *tmp_mat,
                                      const mangll_locidx_t elid,
                                      mangll_t *mangll)
{
  double             *_sc_restrict ind_data = ind_mat->e[0];
  double             *_sc_restrict ind_mass_data = ind_mass_mat->e[0];
  const size_t        total_size = ind_mat->m * ind_mat->n;
  size_t              i;
  double              sum;

  /* check input */
  RHEA_ASSERT (sc_dmatrix_is_valid (ind_mat));

  /* apply mass */
  sc_dmatrix_copy (ind_mat, ind_mass_mat);
  ymir_mass_apply_gll_elem (ind_mass_mat, tmp_mat, elid, mangll);
  RHEA_ASSERT (sc_dmatrix_is_valid (ind_mass_mat));

  /* sum elementwise products */
  sum = 0.0;
  for (i = 0; i < total_size; i++) {
    sum += ind_data[i] * ind_mass_data[i];
  }

  /* return L2-norm */
  return sqrt (sum);
}

static int
rhea_stokes_problem_amr_check_flag (rhea_p4est_quadrant_data_t *qd)
{
  int               check_success;

  if (rhea_amr_flag_is_valid (qd->amr_flag)) { /* if flagged previously */
    switch (qd->amr_flag) {
    case RHEA_AMR_FLAG_COARSEN: /* assume: 1:2-balancing cancelled coarsening */
      qd->amr_flag = RHEA_AMR_FLAG_NO_CHANGE;
      check_success = 1;
      break;
    case RHEA_AMR_FLAG_NO_CHANGE:
      check_success = 1;
      break;
    case RHEA_AMR_FLAG_REFINE: /* error: quadrant should have been refined */
      RHEA_ABORT_NOT_REACHED ();
      check_success = 0;
      break;
    default: /* unknown flag */
      RHEA_ABORT_NOT_REACHED ();
    }
  }
  else {
    /* return that the flag for this quadrant needs to be computed */
    check_success = 0;
  }

  return check_success;
}

static rhea_amr_flag_t
rhea_stokes_problem_amr_get_flag (const double ind,
                                  const double tol_min,
                                  const double tol_max,
                                  const int level,
                                  const int level_min,
                                  const int level_max)
{
  int                 level_min_reached = 0;
  int                 level_max_reached = 0;
  rhea_amr_flag_t     flag;

  /* check input */
  RHEA_ASSERT (0 <= level);

  /* no change at coarsest level */
  if (0 == level) {
    level_min_reached = 1;
  }

  /* coarsen if the level is above max */
  if (0 < level_max) {
    if (level_max < level) {
      return RHEA_AMR_FLAG_COARSEN;
    }
    else if (level_max == level) {
      level_max_reached = 1;
    }
  }

  /* refine if the level is below min */
  if (0 < level_min) {
    if (level < level_min) {
      return RHEA_AMR_FLAG_REFINE;
    }
    else if (level == level_min) {
      level_min_reached = 1;
    }
  }

  /* set flag based on indicator value */
  flag = RHEA_AMR_FLAG_NO_CHANGE;
  if (isfinite (ind) && 0.0 <= ind) { /* if indicator exists */
    const int           tol_min_exists = (isfinite (tol_min) && 0.0 < tol_min);
    const int           tol_max_exists = (isfinite (tol_max) && 0.0 < tol_max);

    RHEA_ASSERT (!(tol_min_exists && tol_max_exists) || tol_min < tol_max);

    if (tol_min_exists && ind < tol_min && !level_min_reached) {
      flag = RHEA_AMR_FLAG_COARSEN;
    }
    else if (tol_max_exists && tol_max < ind && !level_max_reached) {
      flag = RHEA_AMR_FLAG_REFINE;
    }
  }

  /* return flag */
  return flag;
}

static int
rhea_stokes_problem_amr_assign_flag (rhea_p4est_quadrant_data_t *qd,
                                     const rhea_amr_flag_t amr_flag,
                                     const int merge_mode,
                                     p4est_locidx_t *n_flagged_coarsen,
                                     p4est_locidx_t *n_flagged_refine)
{
  int                 flagged;

  if (!merge_mode) { /* if overwrite previous AMR flag */
    switch (amr_flag) {
    case RHEA_AMR_FLAG_COARSEN:
      (*n_flagged_coarsen)++;
      flagged = 1;
      break;
    case RHEA_AMR_FLAG_NO_CHANGE:
      flagged = 0;
      break;
    case RHEA_AMR_FLAG_REFINE:
      (*n_flagged_refine)++;
      flagged = 1;
      break;
    default: /* unknown flag */
      RHEA_ABORT_NOT_REACHED ();
    }
    qd->amr_flag = amr_flag;
  }
  else { /* if merge AMR flag with previous flag */
    switch (amr_flag) {
    case RHEA_AMR_FLAG_COARSEN: /* coarsen + (COAR,NO,REF) => (COAR,NO,REF) */
      flagged = 0;
      break;

    case RHEA_AMR_FLAG_NO_CHANGE:
      switch (qd->amr_flag) {
      case RHEA_AMR_FLAG_COARSEN: /* no change + coarsen => no change */
        qd->amr_flag = amr_flag;
        (*n_flagged_coarsen)--;
        flagged = -1;
        break;
      case RHEA_AMR_FLAG_NO_CHANGE: /* no change + (NO,REF) => (NO,REF) */
      case RHEA_AMR_FLAG_REFINE:
        flagged = 0;
        break;
      default: /* unknown flag */
        RHEA_ABORT_NOT_REACHED ();
      }
      break;

    case RHEA_AMR_FLAG_REFINE:
      switch (amr_flag) {
      case RHEA_AMR_FLAG_COARSEN: /* refine + coarsen => refine */
        qd->amr_flag = amr_flag;
        (*n_flagged_coarsen)--;
        (*n_flagged_refine)++;
        flagged = 0;
        break;
      case RHEA_AMR_FLAG_NO_CHANGE: /* refine + no change => refine */
        qd->amr_flag = amr_flag;
        (*n_flagged_refine)++;
        flagged = 1;
        break;
      case RHEA_AMR_FLAG_REFINE: /* refine + refine => refine */
        flagged = 0;
        break;
      default: /* unknown flag */
        RHEA_ABORT_NOT_REACHED ();
      }
      break;

    default: /* unknown flag */
      RHEA_ABORT_NOT_REACHED ();
    }
  }

  return flagged;
}

/**
 * Prints statistics of AMR indicators and number of flagged elements.
 */
#if (1 <= RHEA_STOKES_PROBLEM_AMR_VERBOSE)
static void
rhea_stokes_problem_amr_print_indicator_statistics (
                                        const double ind_loc_min,
                                        const double ind_loc_max,
                                        const double ind_loc_sum,
                                        const p4est_locidx_t n_flagged_coarsen,
                                        const p4est_locidx_t n_flagged_refine,
                                        const char *func_name,
                                        p4est_t *p4est)
{
  sc_MPI_Comm         mpicomm = p4est->mpicomm;
  int                 mpiret;
  double              ind_glo_min, ind_glo_max, ind_glo_sum;
  int64_t             n_loc[2], n_glo[2];

  mpiret = sc_MPI_Allreduce (&ind_loc_min, &ind_glo_min, 1, sc_MPI_DOUBLE,
                             sc_MPI_MIN, mpicomm); SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Allreduce (&ind_loc_max, &ind_glo_max, 1, sc_MPI_DOUBLE,
                             sc_MPI_MAX, mpicomm); SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Allreduce (&ind_loc_sum, &ind_glo_sum, 1, sc_MPI_DOUBLE,
                             sc_MPI_SUM, mpicomm); SC_CHECK_MPI (mpiret);
  RHEA_GLOBAL_INFOF ("%s: min %.3e, max %.3e, mean %.3e\n",
                     func_name, ind_glo_min, ind_glo_max,
                     ind_glo_sum / (double) p4est->global_num_quadrants);

  n_loc[0] = (int64_t) n_flagged_coarsen;
  n_loc[1] = (int64_t) n_flagged_refine;
  mpiret = sc_MPI_Allreduce (n_loc, &n_glo, 2, MPI_INT64_T, sc_MPI_SUM,
                             mpicomm); SC_CHECK_MPI (mpiret);
  //TODO `sc_MPI_INT64_T` does not exist
  RHEA_GLOBAL_INFOF ("%s: #flagged coarsen %li, refine %li, sum %li\n",
                     func_name, n_glo[0], n_glo[1], n_glo[0] + n_glo[1]);
}
#endif

static double
rhea_stokes_problem_amr_flag_weakzone_peclet_fn (p4est_t *p4est, void *data)
{
  rhea_stokes_problem_amr_data_t *d = data;
  rhea_stokes_problem_t *stokes_problem = d->stokes_problem;
  mangll_t           *mangll = d->mangll_original;
  int                 n_nodes_per_el, nodeid;

  rhea_amr_flag_t     amr_flag;
  p4est_locidx_t      n_flagged_coarsen, n_flagged_refine, n_flagged;
  p4est_topidx_t      ti;
  size_t              tqi;

  sc_dmatrix_t       *ind_mat, *ind_mass_mat, *tmp_mat;
  double             *ind_data;
  double              ind, elem_vol;
  const double        tol_min = d->tol_min;
  const double        tol_max = d->tol_max;
  const int           level_min = d->discr_options->level_min;
  const int           level_max = d->discr_options->level_max;

  rhea_domain_options_t *domain_options =
    rhea_stokes_problem_get_domain_options (stokes_problem);
  const double        domain_size = pow (domain_options->depth, 3.0);
  rhea_weakzone_options_t *weak_options =
    rhea_stokes_problem_get_weakzone_options (stokes_problem);
  const double        weak_thickness = weak_options->thickness;
  const double        weak_thickness_const = weak_options->thickness_const;
  double              weak_factor_interior = weak_options->weak_factor_interior;
#if (1 <= RHEA_STOKES_PROBLEM_AMR_VERBOSE)
  double              ind_loc_min = DBL_MAX;
  double              ind_loc_max = 0.0;
  double              ind_loc_sum = 0.0;
#endif
#if (3 <= RHEA_STOKES_PROBLEM_AMR_VERBOSE_VTK)
  int                 write_vtk = 0;
  sc_dmatrix_t       *indicator_vtk;
  char                debug_path[BUFSIZ], path[BUFSIZ];

  /* set up VTK output of AMR indicator */
  if (0 < ymir_vtk_get_debug_path (debug_path) && d->ymir_mesh != NULL) {
    write_vtk = 1;
    indicator_vtk = sc_dmatrix_new (d->ymir_mesh->cnodes->K, 1);
    snprintf (path, BUFSIZ, "%s_%s", debug_path, __func__);
  }
#endif

  /* check input */
  RHEA_ASSERT (d->mangll_original != NULL);

  /* create work variables */
  n_nodes_per_el = mangll->Np;
  tmp_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  ind_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  ind_data = ind_mat->e[0];
  ind_mass_mat = sc_dmatrix_new (n_nodes_per_el, 1);

  /* flag quadrants */
  n_flagged_coarsen = 0;
  n_flagged_refine = 0;
  for (ti = p4est->first_local_tree; ti <= p4est->last_local_tree; ++ti) {
    p4est_tree_t       *t = p4est_tree_array_index (p4est->trees, ti);
    sc_array_t         *tquadrants = &(t->quadrants);
    const size_t        tqoffset = t->quadrants_offset;

    for (tqi = 0; tqi < tquadrants->elem_count; ++tqi) {
      p4est_quadrant_t   *q = p4est_quadrant_array_index (tquadrants, tqi);
      rhea_p4est_quadrant_data_t *qd = q->p.user_data;
      const mangll_locidx_t elid = (mangll_locidx_t) (tqoffset + tqi);

      /* check if quadrant was flagged previously */
      if (rhea_stokes_problem_amr_check_flag (qd)) {
        continue;
      }

      /* compute element size */
      sc_dmatrix_set_value (ind_mat, 1.0);
      elem_vol = rhea_stokes_problem_amr_norm_L2_elem (ind_mat, ind_mass_mat,
                                                       tmp_mat, elid, mangll);

      /* compute indicator field for this element */
      {
        const double *_sc_restrict x = mangll_get_elem_coord_x (elid, mangll);
        const double *_sc_restrict y = mangll_get_elem_coord_y (elid, mangll);
        const double *_sc_restrict z = mangll_get_elem_coord_z (elid, mangll);

        /* compute quotient |grad fn|/fn directly */
        for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
          int                 label;
          const double        dist =
            rhea_weakzone_dist_node (
                &label, &weak_factor_interior,
                x[nodeid], y[nodeid], z[nodeid], weak_options);
          const double        w =
            rhea_weakzone_factor_node (
                dist, weak_thickness, weak_thickness_const,
                weak_factor_interior);
          const double        dw =
            rhea_weakzone_factor_deriv_node (
                dist, weak_thickness, weak_thickness_const,
                weak_factor_interior);

          ind_data[nodeid] = dw/w;
        }
      }

      /* set indicator to be the L2-norm of the quotient |grad fn|/fn */
      ind = rhea_stokes_problem_amr_norm_L2_elem (ind_mat, ind_mass_mat,
                                                  tmp_mat, elid, mangll);

      /* multiply indicator by relative element size */
      ind *= elem_vol / domain_size;
#if (1 <= RHEA_STOKES_PROBLEM_AMR_VERBOSE)
      ind_loc_min = SC_MIN (ind_loc_min, ind);
      ind_loc_max = SC_MAX (ind_loc_max, ind);
      ind_loc_sum += ind;
#endif
#if (3 <= RHEA_STOKES_PROBLEM_AMR_VERBOSE_VTK)
      if (write_vtk) {
        indicator_vtk->e[elid][0] = ind;
      }
#endif

      /* set flag for coarsening/refinement */
      amr_flag = rhea_stokes_problem_amr_get_flag (ind, tol_min, tol_max,
                                                   (int) q->level,
                                                   level_min, level_max);
      rhea_stokes_problem_amr_assign_flag (qd, amr_flag, d->merge_mode,
                                           &n_flagged_coarsen,
                                           &n_flagged_refine);
    }
  }
  n_flagged = n_flagged_coarsen + n_flagged_refine;

  /* destroy work variables */
  sc_dmatrix_destroy (tmp_mat);
  sc_dmatrix_destroy (ind_mat);
  sc_dmatrix_destroy (ind_mass_mat);

  /* VTK output of AMR indicator */
#if (3 <= RHEA_STOKES_PROBLEM_AMR_VERBOSE_VTK)
  if (write_vtk) {
    ymir_vec_t         *indicator_vec;

    indicator_vec = ymir_dvec_new_element_data (d->ymir_mesh, YMIR_GLL_NODE,
                                                indicator_vtk);
    ymir_vtk_write (d->ymir_mesh, path, indicator_vec, "indicator", NULL);
    ymir_vec_destroy (indicator_vec);
    sc_dmatrix_destroy (indicator_vtk);
  }
#endif

  /* print statistics */
#if (1 <= RHEA_STOKES_PROBLEM_AMR_VERBOSE)
  rhea_stokes_problem_amr_print_indicator_statistics (
      ind_loc_min, ind_loc_max, ind_loc_sum, n_flagged_coarsen,
      n_flagged_refine, __func__, p4est);
#endif

  /* return relative number of flagged quadrants */
  return rhea_amr_get_relative_global_num_flagged (n_flagged, p4est);
}

static double
rhea_stokes_problem_amr_flag_viscosity_peclet_fn (p4est_t *p4est, void *data)
{
  rhea_stokes_problem_amr_data_t *d = data;
  rhea_stokes_problem_t *stokes_problem = d->stokes_problem;
  mangll_t           *mangll = d->mangll_original;
  int                *Vmask;
  int                 n_nodes_per_el;

  rhea_amr_flag_t     amr_flag;
  p4est_locidx_t      n_flagged_coarsen, n_flagged_refine, n_flagged;
  p4est_topidx_t      ti;
  size_t              tqi;

  sc_dmatrix_t       *ind_mat, *ind_mass_mat, *tmp_mat;
  double              ind, elem_vol;
  const double        tol_min = d->tol_min;
  const double        tol_max = d->tol_max;
  const int           level_min = d->discr_options->level_min;
  const int           level_max = d->discr_options->level_max;

  rhea_domain_options_t *domain_options =
    rhea_stokes_problem_get_domain_options (stokes_problem);
  const double        domain_size = pow (domain_options->depth, 3.0);
  rhea_weakzone_options_t *weak_options =
    rhea_stokes_problem_get_weakzone_options (stokes_problem);
  rhea_viscosity_options_t *visc_options =
    rhea_stokes_problem_get_viscosity_options (stokes_problem);
  const int           has_visc = (d->ymir_mesh != NULL);
  const int           has_temp = (d->temperature_original != NULL);
  const int           has_weak = rhea_weakzone_exists (weak_options);
  const int           has_vel = (d->velocity_original != NULL);
  const int           nonlinear_init =
    (visc_options->type == RHEA_VISCOSITY_NONLINEAR && !has_vel);

  ymir_vec_t         *viscosity = NULL;
  sc_dmatrix_t       *weak_el_mat = NULL, *vel_el_mat = NULL,
                     *strt_sqrt_2inv_el_mat = NULL;
  double             *temp_el_data = NULL, *weak_el_data = NULL,
                     *strt_sqrt_2inv_el_data = NULL;
  sc_dmatrix_t       *strt_el_mat = NULL, *geo_el_mat = NULL;
  sc_dmatrix_t       *ur_el_mat, *us_el_mat, *ut_el_mat;
  sc_dmatrix_t       *visc_el_mat;
  double             *visc_el_data, *proj_scal_el_data,
                     *bounds_el_data, *yielding_el_data;
  double             *x = NULL, *y = NULL, *z = NULL;
  double             *tmp_data = NULL;
#if (1 <= RHEA_STOKES_PROBLEM_AMR_VERBOSE)
  double              ind_loc_min = DBL_MAX;
  double              ind_loc_max = 0.0;
  double              ind_loc_sum = 0.0;
#endif
#if (3 <= RHEA_STOKES_PROBLEM_AMR_VERBOSE_VTK)
  int                 write_vtk = 0;
  sc_dmatrix_t       *indicator_vtk;
  char                debug_path[BUFSIZ], path[BUFSIZ];

  /* set up VTK output of AMR indicator */
  if (0 < ymir_vtk_get_debug_path (debug_path) && d->ymir_mesh != NULL) {
    write_vtk = 1;
    indicator_vtk = sc_dmatrix_new (d->ymir_mesh->cnodes->K, 1);
    snprintf (path, BUFSIZ, "%s_%s", debug_path, __func__);
  }
#endif

  /* check input */
  RHEA_ASSERT (d->mangll_original != NULL);

  /* create work variables */
  Vmask = mangll->refel->Vmask;
  n_nodes_per_el = mangll->Np;
  tmp_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  ind_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  ind_mass_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  /* *INDENT-OFF* */
  if (has_visc) {
    viscosity = ymir_dvec_new (d->ymir_mesh, 1, YMIR_GAUSS_NODE);
    rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);
  }
  else {
    if (has_weak) {
      weak_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
      weak_el_data = weak_el_mat->e[0];
    }
    if (has_vel) {
      vel_el_mat            = sc_dmatrix_new (n_nodes_per_el, 3);
      strt_sqrt_2inv_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
      strt_el_mat           = sc_dmatrix_new (n_nodes_per_el, 6);
      geo_el_mat            = sc_dmatrix_new (n_nodes_per_el, 24);
      strt_sqrt_2inv_el_data = strt_sqrt_2inv_el_mat->e[0];
    }
    x = RHEA_ALLOC (double, n_nodes_per_el);
    y = RHEA_ALLOC (double, n_nodes_per_el);
    z = RHEA_ALLOC (double, n_nodes_per_el);
    tmp_data = tmp_mat->e[0];
  }
  ur_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
  us_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
  ut_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
  visc_el_mat       = sc_dmatrix_new (n_nodes_per_el, 1);
  visc_el_data      = visc_el_mat->e[0];
  proj_scal_el_data = NULL;
  bounds_el_data    = NULL;
  yielding_el_data  = NULL;
  /* *INDENT-ON* */

  /* flag quadrants */
  n_flagged_coarsen = 0;
  n_flagged_refine = 0;
  for (ti = p4est->first_local_tree; ti <= p4est->last_local_tree; ++ti) {
    p4est_tree_t       *t = p4est_tree_array_index (p4est->trees, ti);
    sc_array_t         *tquadrants = &(t->quadrants);
    const size_t        tqoffset = t->quadrants_offset;

    for (tqi = 0; tqi < tquadrants->elem_count; ++tqi) {
      p4est_quadrant_t   *q = p4est_quadrant_array_index (tquadrants, tqi);
      rhea_p4est_quadrant_data_t *qd = q->p.user_data;
      const mangll_locidx_t elid = (mangll_locidx_t) (tqoffset + tqi);

      /* check if quadrant was flagged previously */
      if (rhea_stokes_problem_amr_check_flag (qd)) {
        continue;
      }

      /* compute element size */
      sc_dmatrix_set_value (ind_mat, 1.0);
      elem_vol = rhea_stokes_problem_amr_norm_L2_elem (ind_mat, ind_mass_mat,
                                                       tmp_mat, elid, mangll);

      /* get/compute viscosity */
      if (has_visc) {
        rhea_viscosity_get_elem_gauss (visc_el_mat, viscosity,
                                       (ymir_locidx_t) elid);
      }
      else {
        /* get coordinates at Gauss nodes */
        mangll_get_elem_coord_gauss (x, y, z, elid, mangll, tmp_data);

        /* get temperature and weak zone at Gauss nodes */
        if (has_temp) {
          temp_el_data = d->temperature_original->e[elid];
        }
        if (has_weak) {
          rhea_weakzone_compute_elem (weak_el_data, x, y, z, n_nodes_per_el,
                                      weak_options);
        }

        /* get velocity; compute 2nd inv. of the strain rate at Gauss nodes */
        if (has_vel) {
          memcpy (vel_el_mat->e[0], d->velocity_original->e[elid],
                  3 * n_nodes_per_el * sizeof (double));
          ymir_elem_strain_rate_2inv_gauss_node (
              strt_sqrt_2inv_el_mat, vel_el_mat, elid, mangll,
              ur_el_mat, us_el_mat, ut_el_mat, geo_el_mat, strt_el_mat);
        }

        /* compute visosity */
        switch (nonlinear_init) {
        case 0:
          rhea_viscosity_compute_elem (
              visc_el_data, proj_scal_el_data, bounds_el_data, yielding_el_data,
              temp_el_data, weak_el_data, strt_sqrt_2inv_el_data, x, y, z,
              n_nodes_per_el, Vmask, visc_options);
          break;
        case 1:
          rhea_viscosity_compute_nonlinear_init_elem (
              visc_el_data, proj_scal_el_data, bounds_el_data, yielding_el_data,
              temp_el_data, weak_el_data, x, y, z,
              n_nodes_per_el, Vmask, visc_options);
          break;
        default:
          RHEA_ABORT_NOT_REACHED ();
        }
      }
      RHEA_ASSERT (sc_dmatrix_is_valid (visc_el_mat));

      /* compute the (pointwise-norm of the) gradient of the viscosity */
      ymir_elem_gradient_norm_gauss_node (ind_mat, visc_el_mat, elid, mangll,
                                          ur_el_mat, us_el_mat, ut_el_mat);
      RHEA_ASSERT (sc_dmatrix_is_valid (ind_mat));

      /* set indicator to be the L2-norm of the quotient |grad fn|/fn */
      sc_dmatrix_dotdivide (visc_el_mat, ind_mat);
      ind = rhea_stokes_problem_amr_norm_L2_elem (ind_mat, ind_mass_mat,
                                                  tmp_mat, elid, mangll);

      /* multiply indicator by relative element size */
      ind *= elem_vol / domain_size;
#if (1 <= RHEA_STOKES_PROBLEM_AMR_VERBOSE)
      ind_loc_min = SC_MIN (ind_loc_min, ind);
      ind_loc_max = SC_MAX (ind_loc_max, ind);
      ind_loc_sum += ind;
#endif
#if (3 <= RHEA_STOKES_PROBLEM_AMR_VERBOSE_VTK)
      if (write_vtk) {
        indicator_vtk->e[elid][0] = ind;
      }
#endif

      /* set flag for coarsening/refinement */
      amr_flag = rhea_stokes_problem_amr_get_flag (ind, tol_min, tol_max,
                                                   (int) q->level,
                                                   level_min, level_max);
      rhea_stokes_problem_amr_assign_flag (qd, amr_flag, d->merge_mode,
                                           &n_flagged_coarsen,
                                           &n_flagged_refine);
    }
  }
  n_flagged = n_flagged_coarsen + n_flagged_refine;

  /* destroy work variables */
  sc_dmatrix_destroy (tmp_mat);
  sc_dmatrix_destroy (ind_mat);
  sc_dmatrix_destroy (ind_mass_mat);
  if (has_visc) {
    ymir_vec_destroy (viscosity);
  }
  else {
    if (has_weak) {
      sc_dmatrix_destroy (weak_el_mat);
    }
    if (has_vel) {
      sc_dmatrix_destroy (vel_el_mat);
      sc_dmatrix_destroy (strt_sqrt_2inv_el_mat);
      sc_dmatrix_destroy (strt_el_mat);
      sc_dmatrix_destroy (geo_el_mat);
    }
    RHEA_FREE (x);
    RHEA_FREE (y);
    RHEA_FREE (z);
  }
  sc_dmatrix_destroy (ur_el_mat);
  sc_dmatrix_destroy (us_el_mat);
  sc_dmatrix_destroy (ut_el_mat);
  sc_dmatrix_destroy (visc_el_mat);

  /* VTK output of AMR indicator */
#if (3 <= RHEA_STOKES_PROBLEM_AMR_VERBOSE_VTK)
  if (write_vtk) {
    ymir_vec_t         *indicator_vec;

    indicator_vec = ymir_dvec_new_element_data (d->ymir_mesh, YMIR_GLL_NODE,
                                                indicator_vtk);
    ymir_vtk_write (d->ymir_mesh, path, indicator_vec, "indicator", NULL);
    ymir_vec_destroy (indicator_vec);
    sc_dmatrix_destroy (indicator_vtk);
  }
#endif

  /* print statistics */
#if (1 <= RHEA_STOKES_PROBLEM_AMR_VERBOSE)
  rhea_stokes_problem_amr_print_indicator_statistics (
      ind_loc_min, ind_loc_max, ind_loc_sum, n_flagged_coarsen,
      n_flagged_refine, __func__, p4est);
#endif

  /* return relative number of flagged quadrants */
  return rhea_amr_get_relative_global_num_flagged (n_flagged, p4est);
}

/******************************************************************************
 * AMR for a Single Field
 *****************************************************************************/

static sc_dmatrix_t *
rhea_stokes_problem_amr_field_to_buffer (ymir_vec_t *vector)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vector);
  const mangll_locidx_t n_elements = ymir_mesh->ma->mesh->K;
  int                 n_fields, n_entries_per_el;
  sc_dmatrix_t       *buffer_mat;
  ymir_vec_t         *buffer_view;

  /* find out #fields in vector */
  if (ymir_vec_is_cvec (vector)) {
    n_fields = vector->ncfields;
  }
  else if (ymir_vec_is_dvec (vector)) {
    n_fields = vector->ndfields;
  }
  else if (ymir_vec_is_evec (vector)) {
    n_fields = 1;
  }
  else { /* unsupported vector type */
    n_fields = 0;
    RHEA_ABORT_NOT_REACHED ();
  }
  n_entries_per_el = n_fields * ymir_mesh->ma->Np;

  /* create buffer and view onto that buffer */
  buffer_mat = sc_dmatrix_new (n_elements, n_entries_per_el);
  buffer_view = ymir_dvec_new_data (ymir_mesh, n_fields, YMIR_GAUSS_NODE,
                                    buffer_mat);

  /* project fields onto Gauss nodes, store in buffer */
  ymir_interp_vec (vector, buffer_view);

  /* destroy */
  ymir_vec_destroy (buffer_view);

  /* return buffer */
  return buffer_mat;
}

static void
rhea_stokes_problem_amr_buffer_to_field (ymir_vec_t *vector,
                                         sc_dmatrix_t *buffer_original,
                                         ymir_pressure_elem_t *press_elem)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vector);
  int                 n_fields;
  ymir_vec_t         *vector_mass = ymir_vec_template (vector);
  ymir_vec_t         *buffer_view;

  /* find out #fields in vector */
  if (ymir_vec_is_cvec (vector)) {
    n_fields = vector->ncfields;
  }
  else if (ymir_vec_is_dvec (vector)) {
    n_fields = vector->ndfields;
  }
  else if (ymir_vec_is_evec (vector)) {
    n_fields = 1;
  }
  else { /* unsupported vector type */
    n_fields = 0;
    RHEA_ABORT_NOT_REACHED ();
  }

  /* create view onto buffer */
  buffer_view = ymir_dvec_new_data (ymir_mesh, n_fields, YMIR_GAUSS_NODE,
                                    buffer_original);

  /* project buffer (Gauss nodes) onto fields */
  ymir_mass_apply_gauss (buffer_view);
  if (!ymir_vec_is_evec (vector)) { /* if apply accurate mass inversion */
    ymir_interp_vec (buffer_view, vector_mass);
    ymir_mass_invert (vector_mass, vector);
  }
  else { /* otherwise invert appoximately */
    RHEA_ASSERT (press_elem != NULL);
    ymir_interp_vec (buffer_view, vector);
    ymir_pressure_vec_lump_mass (vector_mass, press_elem);
    ymir_vec_divide_in (vector_mass, vector);
  }

  /* destroy */
  ymir_vec_destroy (vector_mass);
  ymir_vec_destroy (buffer_view);
}

static sc_dmatrix_t *
rhea_stokes_problem_amr_project_field (sc_dmatrix_t *buffer_original,
                                       const int n_fields,
                                       mangll_t *mangll_original,
                                       mangll_t *mangll_adapted)
{
  /* adapted mesh */
  const mangll_locidx_t  n_elements = mangll_adapted->mesh->K;
  const int           n_nodes_per_el = mangll_adapted->Np;
  sc_dmatrix_t       *buffer_adapted;

  /* create buffer for adapted data */
  buffer_adapted = sc_dmatrix_new (n_elements, n_fields * n_nodes_per_el);

  /* project original data */
  ymir_hmg_intergrid_h_project_gauss (
      buffer_adapted, mangll_adapted,
      buffer_original, mangll_original, n_fields);

  /* destroy buffer for original data */
  sc_dmatrix_destroy (buffer_original);

  /* return buffer for adapted data */
  return buffer_adapted;
}

static sc_dmatrix_t *
rhea_stokes_problem_amr_partition_field (sc_dmatrix_t *buffer_adapted,
                                         const int n_fields,
                                         mangll_t *mangll_adapted,
                                         mangll_t *mangll_partitioned)
{
  /* adapted mesh */
  mangll_gloidx_t    *RtoGEO_adapted = mangll_adapted->mesh->RtoGEO;
  sc_MPI_Comm         mpicomm = mangll_adapted->mpicomm;
  const int           mpisize = mangll_adapted->mpisize;
  const int           mpirank = mangll_adapted->mpirank;
  /* partitoned mesh */
  mangll_gloidx_t    *RtoGEO_partitioned = mangll_partitioned->mesh->RtoGEO;
  const mangll_locidx_t  n_elements = mangll_partitioned->mesh->K;
  const int           n_nodes_per_el = mangll_partitioned->Np;
  sc_dmatrix_t       *buffer_partitioned;

  /* create buffer for partitioned data */
  buffer_partitioned = sc_dmatrix_new (n_elements, n_fields * n_nodes_per_el);

  /* partition adapted data */
  mangll_field_partition (n_fields, n_nodes_per_el, mpirank, mpisize, mpicomm,
                          RtoGEO_adapted, buffer_adapted,
                          RtoGEO_partitioned, buffer_partitioned);

  /* destroy buffer for adapted data */
  sc_dmatrix_destroy (buffer_adapted);

  /* return buffer for partitioned data */
  return buffer_partitioned;
}

/******************************************************************************
 * AMR Callback Functions
 *****************************************************************************/

static void
rhea_stokes_problem_amr_data_initialize_fn (p4est_t *p4est, void *data)
{
  rhea_stokes_problem_amr_data_t *d = data;
  ymir_vec_t         *temperature =
    rhea_stokes_problem_get_temperature (d->stokes_problem);
  ymir_vec_t         *velocity_pressure =
    rhea_stokes_problem_get_velocity_pressure (d->stokes_problem);
  const int           has_temp = (temperature != NULL);
  const int           has_vel_press = (velocity_pressure != NULL);

  RHEA_GLOBAL_INFOF ("Into %s (temperature %i, velocity-pressure %i)\n",
                     __func__, has_temp, has_vel_press);

  /* check input */
  RHEA_ASSERT (d->mangll_original == d->ymir_mesh->ma);
  RHEA_ASSERT (d->mangll_adapted == NULL);
  RHEA_ASSERT (d->mangll_partitioned == NULL);
  RHEA_ASSERT (d->ymir_mesh != NULL);
  RHEA_ASSERT (d->press_elem != NULL);

  /*
   * Fields & Stokes Problem
   */

  /* store fields in buffers `*_original` */
  if (has_temp) {
    RHEA_ASSERT (rhea_temperature_check_vec_type (temperature));
    RHEA_ASSERT (d->temperature_original == NULL);

    d->temperature_original = rhea_stokes_problem_amr_field_to_buffer (
        temperature);
    RHEA_ASSERT (d->temperature_original->n == d->ymir_mesh->ma->Np);
  }
  if (has_vel_press) {
    ymir_vec_t         *velocity, *pressure;

    RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (velocity_pressure));
    RHEA_ASSERT (d->velocity_original == NULL);
    RHEA_ASSERT (d->pressure_original == NULL);

    rhea_velocity_pressure_create_components (
        &velocity, &pressure, velocity_pressure, d->press_elem);
    d->velocity_original = rhea_stokes_problem_amr_field_to_buffer (velocity);
    d->pressure_original = rhea_stokes_problem_amr_field_to_buffer (pressure);
    RHEA_ASSERT (d->velocity_original->n == 3 * d->ymir_mesh->ma->Np);
    RHEA_ASSERT (d->pressure_original->n == d->ymir_mesh->ma->Np);

    rhea_velocity_destroy (velocity);
    rhea_pressure_destroy (pressure);
  }

  /* destroy */
  rhea_stokes_problem_amr_data_clear_fields (d);
  rhea_stokes_problem_clear_mesh_dependencies (d->stokes_problem);

  /*
   * Mesh
   */

  /* set AMR parameters */
  d->first_amr = 1;

  /* destroy mesh partially (keep only the mangll object) */
  rhea_discretization_mangll_continuous_destroy (
      NULL, d->ymir_mesh->cnodes);
  ymir_mesh_destroy (d->ymir_mesh);
  ymir_pressure_elem_destroy (d->press_elem);
  d->ymir_mesh = NULL;
  d->press_elem = NULL;

  RHEA_GLOBAL_INFOF ("Done %s\n", __func__);
}

static void
rhea_stokes_problem_amr_data_finalize_fn (p4est_t *p4est, void *data)
{
  rhea_stokes_problem_amr_data_t *d = data;
  const int           has_temp = (d->temperature_original != NULL);
  const int           has_vel = (d->velocity_original != NULL);
  const int           has_press = (d->pressure_original != NULL);
  ymir_vec_t         *temperature, *velocity_pressure;

  RHEA_GLOBAL_INFOF ("Into %s (temperature %i, velocity-pressure %i-%i)\n",
                     __func__, has_temp, has_vel, has_press);
  /* check input */
  RHEA_ASSERT (d->mangll_original != NULL);
  RHEA_ASSERT (d->mangll_adapted == NULL);
  RHEA_ASSERT (d->mangll_partitioned == NULL);
  RHEA_ASSERT (d->ymir_mesh == NULL);
  RHEA_ASSERT (d->press_elem == NULL);

  /*
   * Mesh
   */

  /* destroy discontinuous mangll */
  RHEA_ASSERT (d->first_amr == 0);
  rhea_discretization_mangll_discontinuous_destroy (d->mangll_original);
  d->mangll_original = NULL;

  /* create ymir mesh */
  rhea_discretization_ymir_mesh_new_from_p4est (
      &(d->ymir_mesh), &(d->press_elem), p4est, d->discr_options);

  /*
   * Fields & Stokes Problem
   */

  /* retrieve fields from buffers `*_original` */
  if (has_temp) {
    temperature = rhea_temperature_new (d->ymir_mesh);
    rhea_stokes_problem_amr_buffer_to_field (
        temperature, d->temperature_original, d->press_elem);
  }
  if (has_vel || has_press) {
    ymir_vec_t         *velocity = rhea_velocity_new (d->ymir_mesh);
    ymir_vec_t         *pressure = rhea_pressure_new (d->ymir_mesh,
                                                      d->press_elem);

    if (has_vel) {
      rhea_stokes_problem_amr_buffer_to_field (
          velocity, d->velocity_original, d->press_elem);
    }
    else {
      ymir_vec_set_zero (velocity);
    }
    if (has_press) {
      rhea_stokes_problem_amr_buffer_to_field (
          pressure, d->pressure_original, d->press_elem);
    }
    else {
      ymir_vec_set_zero (pressure);
    }

    velocity_pressure = rhea_velocity_pressure_new (d->ymir_mesh,
                                                    d->press_elem);
    rhea_velocity_pressure_set_components (
        velocity_pressure, velocity, pressure, d->press_elem);

    rhea_velocity_destroy (velocity);
    rhea_pressure_destroy (pressure);
  }

  /* destroy buffers */
  rhea_stokes_problem_amr_data_clear_buffer (d);

  /* init Stokes problem */
  if (has_temp) {
    rhea_stokes_problem_set_temperature (d->stokes_problem, temperature);
  }
  if (has_vel || has_press) {
    rhea_stokes_problem_set_velocity_pressure (d->stokes_problem,
                                               velocity_pressure);
  }
  rhea_stokes_problem_create_mesh_dependencies (d->stokes_problem,
                                                d->ymir_mesh, d->press_elem);

  RHEA_GLOBAL_INFOF ("Done %s\n", __func__);
}

static void
rhea_stokes_problem_amr_data_project_fn (p4est_t *p4est, void *data)
{
  rhea_stokes_problem_amr_data_t *d = data;
  const int           has_temp = (d->temperature_original != NULL);
  const int           has_vel = (d->velocity_original != NULL);
  const int           has_press = (d->pressure_original != NULL);

  RHEA_GLOBAL_INFOF ("Into %s (temperature %i, velocity-pressure %i-%i)\n",
                     __func__, has_temp, has_vel, has_press);

  /* check input */
  RHEA_ASSERT (d->mangll_original != NULL);
  RHEA_ASSERT (d->mangll_adapted == NULL);
  RHEA_ASSERT (d->mangll_partitioned == NULL);
  RHEA_ASSERT (d->first_amr == 0 || d->first_amr == 1);

  /* create adapted mangll object */
  d->mangll_adapted = rhea_discretization_mangll_discontinuous_new (
      p4est, d->discr_options);
  RHEA_ASSERT (d->mangll_adapted != NULL);
  RHEA_ASSERT (d->mangll_adapted->N == d->mangll_original->N);

  /* project from original to adapted */
  if (has_temp) {
    RHEA_ASSERT (d->temperature_adapted == NULL);
    d->temperature_adapted = rhea_stokes_problem_amr_project_field (
        d->temperature_original, 1 /* #fields */,
        d->mangll_original, d->mangll_adapted);
    d->temperature_original = NULL;

    /* bound values */
    rhea_temperature_bound_mat (d->temperature_adapted);
  }
  if (has_vel) {
    RHEA_ASSERT (d->velocity_adapted == NULL);
    d->velocity_adapted = rhea_stokes_problem_amr_project_field (
        d->velocity_original, 3 /* #fields */,
        d->mangll_original, d->mangll_adapted);
    d->velocity_original = NULL;
  }
  if (has_press) {
    RHEA_ASSERT (d->pressure_adapted == NULL);
    d->pressure_adapted = rhea_stokes_problem_amr_project_field (
        d->pressure_original, 1 /* #fields */,
        d->mangll_original, d->mangll_adapted);
    d->pressure_original = NULL;
  }

  /* destroy original mangll */
  if (d->first_amr) {
    rhea_discretization_mangll_continuous_destroy (d->mangll_original, NULL);
    d->first_amr = 0;
  }
  else {
    rhea_discretization_mangll_discontinuous_destroy (d->mangll_original);
  }
  d->mangll_original = NULL;

  RHEA_GLOBAL_INFOF ("Done %s\n", __func__);
}

static void
rhea_stokes_problem_amr_data_partition_fn (p4est_t *p4est, void *data)
{
  rhea_stokes_problem_amr_data_t *d = data;
  const int           has_temp = (d->temperature_adapted != NULL);
  const int           has_vel = (d->velocity_adapted != NULL);
  const int           has_press = (d->pressure_adapted != NULL);

  RHEA_GLOBAL_INFOF ("Into %s (temperature %i, velocity-pressure %i-%i)\n",
                     __func__, has_temp, has_vel, has_press);

  /* check input */
  RHEA_ASSERT (d->mangll_original == NULL);
  RHEA_ASSERT (d->mangll_adapted != NULL);
  RHEA_ASSERT (d->mangll_partitioned == NULL);

  /* create partitioned mangll object */
  d->mangll_partitioned = rhea_discretization_mangll_discontinuous_new (
      p4est, d->discr_options);
  RHEA_ASSERT (d->mangll_partitioned != NULL);
  RHEA_ASSERT (d->mangll_partitioned->N == d->mangll_adapted->N);

  /* partition from adapted to original */
  if (has_temp) {
    RHEA_ASSERT (d->temperature_original == NULL);
    d->temperature_original = rhea_stokes_problem_amr_partition_field (
        d->temperature_adapted, 1 /* #fields*/,
        d->mangll_adapted, d->mangll_partitioned);
    d->temperature_adapted = NULL;
  }
  if (has_vel) {
    RHEA_ASSERT (d->velocity_original == NULL);
    d->velocity_original = rhea_stokes_problem_amr_partition_field (
        d->velocity_adapted, 3 /* #fields*/,
        d->mangll_adapted, d->mangll_partitioned);
    d->velocity_adapted = NULL;
  }
  if (has_press) {
    RHEA_ASSERT (d->pressure_original == NULL);
    d->pressure_original = rhea_stokes_problem_amr_partition_field (
        d->pressure_adapted, 1 /* #fields*/,
        d->mangll_adapted, d->mangll_partitioned);
    d->pressure_adapted = NULL;
  }

  /* destroy adapted but unpartitioned mangll */
  rhea_discretization_mangll_discontinuous_destroy (d->mangll_adapted);
  d->mangll_adapted = NULL;

  /* reassign original mangll for next projection */
  d->mangll_original = d->mangll_partitioned;
  d->mangll_partitioned = NULL;

  RHEA_GLOBAL_INFOF ("Done %s\n", __func__);
}

/******************************************************************************
 * AMR Main Functions
 *****************************************************************************/

int
rhea_stokes_problem_init_amr (rhea_stokes_problem_t *stokes_problem,
                              p4est_t *p4est,
                              rhea_discretization_options_t *discr_options)
{
  const char         *type_name = rhea_stokes_problem_amr_init_type_name;
  const double        tol_min = rhea_stokes_problem_amr_init_tol_min;
  const double        tol_max = rhea_stokes_problem_amr_init_tol_max;
  const int           n_cycles = rhea_stokes_problem_amr_n_cycles;
  const double        flagged_elements_thresh_begin =
    rhea_stokes_problem_amr_flagged_elements_thresh_begin;
  const double        flagged_elements_thresh_cycle =
    rhea_stokes_problem_amr_flagged_elements_thresh_cycle;

  rhea_stokes_problem_amr_data_t *amr_data;
  rhea_amr_flag_elements_fn_t     flag_fn;
  void                           *flag_fn_data;
  int                 amr_iter;

  /* exit if nothing to do */
  if (strcmp (type_name, "NONE") == 0) {
    return 0;
  }

  RHEA_GLOBAL_INFOF ("Into %s (%s)\n", __func__, type_name);

  /* check input */
  RHEA_ASSERT (p4est != NULL);
  RHEA_ASSERT (discr_options != NULL);

  /* create AMR data */
  amr_data = rhea_stokes_problem_amr_data_new (stokes_problem, discr_options);

  /* set up flagging of elements for coarsening/refinement */
  if (strcmp (type_name, "weakzone_peclet") == 0) {
    flag_fn = rhea_stokes_problem_amr_flag_weakzone_peclet_fn;
    flag_fn_data = amr_data;
  }
  else if (strcmp (type_name, "viscosity_peclet") == 0) {
    flag_fn = rhea_stokes_problem_amr_flag_viscosity_peclet_fn;
    flag_fn_data = amr_data;
    rhea_stokes_problem_compute_coefficient (stokes_problem, 1 /* nl. init */);
  }
  else { /* unknown type name */
    RHEA_ABORT ("Unknown initial AMR name");
    flag_fn = NULL;
    flag_fn_data = NULL;
  }
  amr_data->tol_min = tol_min;
  amr_data->tol_max = tol_max;

  /* perform AMR */
  amr_iter = rhea_amr (p4est, n_cycles, flagged_elements_thresh_begin,
                       flagged_elements_thresh_cycle, flag_fn, flag_fn_data,
                       rhea_stokes_problem_amr_data_initialize_fn,
                       rhea_stokes_problem_amr_data_finalize_fn,
                       rhea_stokes_problem_amr_data_project_fn,
                       rhea_stokes_problem_amr_data_partition_fn, amr_data);

  /* destroy */
  rhea_stokes_problem_amr_data_destroy (amr_data);

  /* print mesh statistics */
  if (0 < amr_iter) {
    ymir_monitor_print_global_mesh_stats (
        rhea_stokes_problem_get_ymir_mesh (stokes_problem),
        rhea_stokes_problem_get_press_elem (stokes_problem));
  }

  RHEA_GLOBAL_INFOF ("Done %s (%s)\n", __func__, type_name);

  /* return number of performed AMR iterations */
  return amr_iter;
}

int
rhea_stokes_problem_nonlinear_amr (rhea_stokes_problem_t *stokes_problem,
                                   p4est_t *p4est,
                                   rhea_discretization_options_t *discr_options,
                                   const int nonlinear_iter)
{
  const char         *type_name = rhea_stokes_problem_amr_nonlinear_type_name;
  const double        tol_min = rhea_stokes_problem_amr_nonlinear_tol_min;
  const double        tol_max = rhea_stokes_problem_amr_nonlinear_tol_max;
  const int           n_cycles = rhea_stokes_problem_amr_n_cycles;
  const double        flagged_elements_thresh_begin =
    rhea_stokes_problem_amr_flagged_elements_thresh_begin;
  const double        flagged_elements_thresh_cycle =
    rhea_stokes_problem_amr_flagged_elements_thresh_cycle;
  const int           iter_first = rhea_stokes_problem_amr_nonlinear_iter_first;
  const int           iter_last = rhea_stokes_problem_amr_nonlinear_iter_last;

  rhea_stokes_problem_amr_data_t *amr_data;
  rhea_amr_flag_elements_fn_t     flag_fn;
  void                           *flag_fn_data;
  int                 amr_iter;

  /* exit if nothing to do */
  if (strcmp (type_name, "NONE") == 0 ||
      nonlinear_iter < iter_first ||
      (0 <= iter_last && iter_last < nonlinear_iter)) {
    return 0;
  }

  RHEA_GLOBAL_INFOF ("Into %s (%s)\n", __func__, type_name);

  /* check input */
  RHEA_ASSERT (p4est != NULL);
  RHEA_ASSERT (discr_options != NULL);

  /* create AMR data */
  amr_data = rhea_stokes_problem_amr_data_new (stokes_problem, discr_options);

  /* set up flagging of elements for coarsening/refinement */
  if (strcmp (type_name, "viscosity_peclet") == 0) {
    flag_fn = rhea_stokes_problem_amr_flag_viscosity_peclet_fn;
    flag_fn_data = amr_data;
  }
  else { /* unknown type name */
    RHEA_ABORT ("Unknown nonlinear AMR name");
    flag_fn = NULL;
    flag_fn_data = NULL;
  }
  amr_data->tol_min = tol_min;
  amr_data->tol_max = tol_max;

  /* perform AMR */
  amr_iter = rhea_amr (p4est, n_cycles, flagged_elements_thresh_begin,
                       flagged_elements_thresh_cycle, flag_fn, flag_fn_data,
                       rhea_stokes_problem_amr_data_initialize_fn,
                       rhea_stokes_problem_amr_data_finalize_fn,
                       rhea_stokes_problem_amr_data_project_fn,
                       rhea_stokes_problem_amr_data_partition_fn, amr_data);

  /* destroy */
  rhea_stokes_problem_amr_data_destroy (amr_data);

  /* print mesh statistics */
  if (0 < amr_iter) {
    ymir_monitor_print_global_mesh_stats (
        rhea_stokes_problem_get_ymir_mesh (stokes_problem),
        rhea_stokes_problem_get_press_elem (stokes_problem));
  }

  RHEA_GLOBAL_INFOF ("Done %s (%s)\n", __func__, type_name);

  /* return number of performed AMR iterations */
  return amr_iter;
}

int
rhea_stokes_problem_amr (rhea_stokes_problem_t *stokes_problem,
                         p4est_t *p4est,
                         rhea_discretization_options_t *discr_options)
{
  //TODO
  const double        tol_min = rhea_stokes_problem_amr_init_tol_min;
  const double        tol_max = rhea_stokes_problem_amr_init_tol_max;
  const int           n_cycles = 2;
  const double        flagged_elements_thresh_begin = NAN;
  const double        flagged_elements_thresh_cycle = NAN;

  rhea_stokes_problem_amr_data_t *amr_data;
  rhea_amr_flag_elements_fn_t     flag_fn;
  void                           *flag_fn_data;
  int                 amr_iter;

  RHEA_GLOBAL_INFOF ("Into %s\n", __func__);

  /* create AMR data */
  amr_data = rhea_stokes_problem_amr_data_new (stokes_problem, discr_options);

  /* set function that flags elements for coarsening/refinement */
  //TODO
  flag_fn = rhea_stokes_problem_amr_flag_viscosity_peclet_fn;
  flag_fn_data = amr_data;
  amr_data->tol_min = tol_min;
  amr_data->tol_max = tol_max;

  /* perform AMR */
  amr_iter = rhea_amr (p4est, n_cycles, flagged_elements_thresh_begin,
                       flagged_elements_thresh_cycle, flag_fn, flag_fn_data,
                       rhea_stokes_problem_amr_data_initialize_fn,
                       rhea_stokes_problem_amr_data_finalize_fn,
                       rhea_stokes_problem_amr_data_project_fn,
                       rhea_stokes_problem_amr_data_partition_fn, amr_data);

  /* destroy */
  rhea_stokes_problem_amr_data_destroy (amr_data);

  /* print mesh statistics */
  if (0 < amr_iter) {
    ymir_monitor_print_global_mesh_stats (
        rhea_stokes_problem_get_ymir_mesh (stokes_problem),
        rhea_stokes_problem_get_press_elem (stokes_problem));
  }

  RHEA_GLOBAL_INFOF ("Done %s\n", __func__);

  /* return number of performed AMR iterations */
  return amr_iter;
}
