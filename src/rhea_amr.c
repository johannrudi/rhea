/*
 */

#include <rhea_amr.h>
#include <rhea_base.h>
#include <rhea_discretization.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <ymir_perf_counter.h>

#define RHEA_AMR_P4EST_REPLACE_FN NULL
#define RHEA_AMR_P4EST_PARTITION_WEIGHT_FN NULL

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_AMR_DEFAULT_N_ELEMENTS_MAX (NAN)
#define RHEA_AMR_DEFAULT_INIT_REFINE_NAME "uniform"
#define RHEA_AMR_DEFAULT_INIT_REFINE_DEPTH_M NULL
#define RHEA_AMR_DEFAULT_INIT_REFINE_LEVEL_MIN (-1)
#define RHEA_AMR_DEFAULT_INIT_REFINE_LEVEL_MAX (-1)
#define RHEA_AMR_DEFAULT_LOG_LEVEL (1)
#define RHEA_AMR_DEFAULT_LOG_LEVEL_MAX (1)
#define RHEA_AMR_DEFAULT_MONITOR_PERFORMANCE (0)

double              rhea_amr_n_elements_max = RHEA_AMR_DEFAULT_N_ELEMENTS_MAX;
char               *rhea_amr_init_refine_name =
  RHEA_AMR_DEFAULT_INIT_REFINE_NAME;
char               *rhea_amr_init_refine_depth_m =
  RHEA_AMR_DEFAULT_INIT_REFINE_DEPTH_M;
int                 rhea_amr_init_refine_level_min =
  RHEA_AMR_DEFAULT_INIT_REFINE_LEVEL_MIN;
int                 rhea_amr_init_refine_level_max =
  RHEA_AMR_DEFAULT_INIT_REFINE_LEVEL_MAX;
int                 rhea_amr_log_level = RHEA_AMR_DEFAULT_LOG_LEVEL;
int                 rhea_amr_log_level_max = RHEA_AMR_DEFAULT_LOG_LEVEL_MAX;
int                 rhea_amr_monitor_performance =
                      RHEA_AMR_DEFAULT_MONITOR_PERFORMANCE;

void
rhea_amr_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "AMR";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_D, "num-elements-max", '\0',
    &(rhea_amr_n_elements_max), RHEA_AMR_DEFAULT_N_ELEMENTS_MAX,
    "Max #elements above which the mesh is not refined further",

  YMIR_OPTIONS_S, "init-refine", '\0',
    &(rhea_amr_init_refine_name), RHEA_AMR_DEFAULT_INIT_REFINE_NAME,
    "Refinement of new/initial mesh: uniform, half, radial",
  YMIR_OPTIONS_S, "init-refine-depth", '\0',
    &(rhea_amr_init_refine_depth_m), RHEA_AMR_DEFAULT_INIT_REFINE_DEPTH_M,
    "Refinement 'depth': Sorted list of depths [m] (put shallowest at last)",

  YMIR_OPTIONS_I, "init-refine-level-min", '\0',
    &(rhea_amr_init_refine_level_min), RHEA_AMR_DEFAULT_INIT_REFINE_LEVEL_MIN,
    "Refinement of new/initial mesh: Minumum level of mesh refinement",
  YMIR_OPTIONS_I, "init-refine-level-max", '\0',
    &(rhea_amr_init_refine_level_max), RHEA_AMR_DEFAULT_INIT_REFINE_LEVEL_MAX,
    "Refinement of new/initial mesh: Maximum level of mesh refinement",

  YMIR_OPTIONS_B, "log-level", '\0',
    &(rhea_amr_log_level), RHEA_AMR_DEFAULT_LOG_LEVEL,
    "Print histogram of mesh refinement levels",
  YMIR_OPTIONS_B, "log-level-max", '\0',
    &(rhea_amr_log_level_max), RHEA_AMR_DEFAULT_LOG_LEVEL_MAX,
    "Print max mesh refinement level",

  YMIR_OPTIONS_B, "monitor-performance", '\0',
    &(rhea_amr_monitor_performance), RHEA_AMR_DEFAULT_MONITOR_PERFORMANCE,
    "Measure and print performance statistics (e.g., runtime or flops)",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);

  /* initialize (deactivated) performance counters */
  rhea_amr_perfmon_init (0, 0);
}

/******************************************************************************
 * Monitoring
 *****************************************************************************/

/* perfomance monitor tags and names */
typedef enum
{
  RHEA_AMR_PERFMON_AMR_FLAG,
  RHEA_AMR_PERFMON_AMR_COARSEN,
  RHEA_AMR_PERFMON_AMR_REFINE,
  RHEA_AMR_PERFMON_AMR_BALANCE,
  RHEA_AMR_PERFMON_AMR_PARTITION,
  RHEA_AMR_PERFMON_AMR_DATA_INIT,
  RHEA_AMR_PERFMON_AMR_DATA_PROJECT,
  RHEA_AMR_PERFMON_AMR_DATA_PARTITION,
  RHEA_AMR_PERFMON_AMR_DATA_FINALIZE,
  RHEA_AMR_PERFMON_AMR_TOTAL,
  RHEA_AMR_PERFMON_N
}
rhea_amr_perfmon_idx_t;

static const char  *rhea_amr_perfmon_name[RHEA_AMR_PERFMON_N] =
{
  "AMR: Flag",
  "AMR: Coarsen",
  "AMR: Refine",
  "AMR: Balance",
  "AMR: Partition",
  "AMR: Data initialize",
  "AMR: Data project",
  "AMR: Data partition",
  "AMR: Data finalize",
  "AMR (total)"
};
ymir_perf_counter_t rhea_amr_perfmon[RHEA_AMR_PERFMON_N];

void
rhea_amr_perfmon_init (const int activate, const int skip_if_active)
{
  const int           active = activate && rhea_amr_monitor_performance;

  ymir_perf_counter_init_all_ext (rhea_amr_perfmon, rhea_amr_perfmon_name,
                                  RHEA_AMR_PERFMON_N, active, skip_if_active);
}

void
rhea_amr_perfmon_print (sc_MPI_Comm mpicomm,
                        const int print_wtime,
                        const int print_n_calls,
                        const int print_flops)
{
  const int           active = rhea_amr_monitor_performance;
  const int           print = (print_wtime || print_n_calls || print_flops);
  int                 n_stats = RHEA_AMR_PERFMON_N * YMIR_PERF_COUNTER_N_STATS;
  sc_statinfo_t       stats[n_stats];
  char                stats_name[n_stats][YMIR_PERF_COUNTER_NAME_SIZE];

  /* exit if nothing to do */
  if (!active || !print) {
    return;
  }

  /* gather performance statistics */
  n_stats = ymir_perf_counter_gather_stats (
      rhea_amr_perfmon, RHEA_AMR_PERFMON_N, stats, stats_name,
      mpicomm, print_wtime, print_n_calls, print_flops);

  /* print performance statistics */
  ymir_perf_counter_print_stats (stats, n_stats, "AMR");
}

/******************************************************************************
 * Initial AMR for p4est
 *****************************************************************************/

/* types of initial refinement for a p4est mesh */
typedef enum
{
  RHEA_AMR_INIT_REFINE_UNKNOWN = -2,
  RHEA_AMR_INIT_REFINE_NONE    = -1,
  RHEA_AMR_INIT_REFINE_UNIFORM =  0, /* default, does not have a refine fn. */
  RHEA_AMR_INIT_REFINE_HALF,
  RHEA_AMR_INIT_REFINE_DEPTH,
  RHEA_AMR_INIT_REFINE_N
}
rhea_amr_init_refine_t;

/* name for types of initial refinement */
const char         *rhea_amr_init_refine_name_list[RHEA_AMR_INIT_REFINE_N] =
{
  "uniform",
  "half",
  "depth"
};

static rhea_amr_init_refine_t
rhea_amr_init_refine_decode (const char * name)
{
  const int           n_types = (int) RHEA_AMR_INIT_REFINE_N;
  int                 typeid;
  rhea_amr_init_refine_t type;

  /* exit if nothing to do */
  if (name == NULL) {
    return RHEA_AMR_INIT_REFINE_NONE;
  }

  /* match name to a type */
  type = RHEA_AMR_INIT_REFINE_UNKNOWN;
  for (typeid = 0; typeid < n_types; typeid++) {
    if (!strcmp (name, rhea_amr_init_refine_name_list[typeid])) {
      type = (rhea_amr_init_refine_t) typeid;
    }
  }

  /* return refinement type */
  return type;
}

int
rhea_amr_init_refine (p4est_t *p4est,
                      int level_min,
                      int level_max,
                      rhea_domain_options_t *domain_options)
{
  /* arguments for refinement */
  const char         *refine_name = rhea_amr_init_refine_name;
  rhea_amr_init_refine_t type;
  p4est_refine_t      refine_fn = NULL;
  void               *refine_data = NULL;
  int                 recursively = 0;
  const p4est_init_t    init_fn = rhea_p4est_init_fn;
  const p4est_replace_t replace_fn = RHEA_AMR_P4EST_REPLACE_FN;
  /* arguments for partitioning */
  const int           partition_for_coarsening = 1;
  p4est_weight_t      partition_weight_fn = RHEA_AMR_P4EST_PARTITION_WEIGHT_FN;

  /* get type of refinement from name; update refinement min/max levels */
  type = rhea_amr_init_refine_decode (refine_name);
  if (0 <= rhea_amr_init_refine_level_min) {
    level_min = SC_MAX (level_min, rhea_amr_init_refine_level_min);
  }
  if (0 <= rhea_amr_init_refine_level_max) {
    level_max = SC_MIN (level_max, rhea_amr_init_refine_level_max);
  }
  RHEA_ASSERT (level_min <= level_max);

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (
      __func__, "refine=\"%s\", type=%i, level_min=%i, level_max=%i",
      refine_name, type, level_min, level_max);

  /* set refinement function and data */
  switch (type) {
  case RHEA_AMR_INIT_REFINE_UNKNOWN:
    RHEA_GLOBAL_LERROR ("Unknown refinement type");
    break;
  case RHEA_AMR_INIT_REFINE_NONE:
    /* no function necessary */
    break;
  case RHEA_AMR_INIT_REFINE_UNIFORM:
    refine_data = &level_min;
    refine_fn = rhea_amr_refine_to_level_fn;
    recursively = 1;
    break;
  case RHEA_AMR_INIT_REFINE_HALF:
    refine_data = &level_min;
    refine_fn = rhea_amr_refine_half_fn;
    recursively = 1;
    break;
  case RHEA_AMR_INIT_REFINE_DEPTH:
    {
      double             *depth_m = NULL;
      int                 count, k;
      const double        rmin = domain_options->radius_min;
      const double        rmax = domain_options->radius_max;
      double              d, depth_rel, p4est_zmax;
      int                *depth;
      rhea_amr_refine_depth_data_t *data;

      /* convert input string to double values */
      count = ymir_options_convert_string_to_double (
          rhea_amr_init_refine_depth_m, &depth_m);

      /* skip refinement if depths were not provided */
      if (count == 0) {
        type = RHEA_AMR_INIT_REFINE_NONE;
        break;
      }

      /* check values */
      for (k = 0; k < count; k++) {
        RHEA_CHECK_ABORT (0.0 <= depth_m[k], "Refinement depth is negative");
        if (0 < k) {
          RHEA_CHECK_ABORT (depth_m[k] <= depth_m[k-1],
                            "Refinement depths are not ascending");
        }
      }

      /* set z-coordinate of the domain's top */
      switch (domain_options->shape) {
      case RHEA_DOMAIN_CUBE:
      case RHEA_DOMAIN_BOX:
      case RHEA_DOMAIN_SHELL:
      case RHEA_DOMAIN_CUBE_SPHERICAL:
        p4est_zmax = NAN;
        refine_fn = rhea_amr_refine_depth_fn;
        break;
      case RHEA_DOMAIN_BOX_SPHERICAL:
        p4est_zmax = 2.0;
        refine_fn = rhea_amr_refine_depth_box_fn;
        break;
      default: /* unknown domain shape */
        RHEA_ABORT_NOT_REACHED ();
      }

      /* transform dimensional depths to p4est quadrant coordinates */
      depth = RHEA_ALLOC (int, count);
      for (k = 0; k < count; k++) {
        d = rhea_domain_depth_m_to_depth (depth_m[k], domain_options);
        depth_rel = d / (rmax - rmin);
        RHEA_ASSERT (0.0 <= depth_rel && depth_rel <= 1.0);
        depth[k] = (int) round (depth_rel * (double) P4EST_ROOT_LEN);
      }
      YMIR_FREE (depth_m); /* was allocated in ymir */

      /* set refinement data */
      data = RHEA_ALLOC (rhea_amr_refine_depth_data_t, 1);
      data->depth = depth;
      data->count = count;
      data->level_min = level_min;
      data->p4est_zmax = p4est_zmax;
      refine_data = data;
      recursively = 1;
    }
    break;
  default: /* unknown refinement type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* perform refinement */
  if (refine_fn != NULL) {
    void               *user_pointer = p4est->user_pointer;

    /* refine */
    p4est->user_pointer = (void *) refine_data;
    p4est_refine_ext (p4est, recursively, level_max, refine_fn, init_fn,
                      replace_fn);
    p4est->user_pointer = user_pointer;

    /* partition and balance */
    p4est_partition_ext (p4est, partition_for_coarsening, partition_weight_fn);
    p4est_balance (p4est, P4EST_CONNECT_FULL, init_fn);
    p4est_partition_ext (p4est, partition_for_coarsening, partition_weight_fn);
  }

  /* destroy refinement data */
  switch (type) {
  case RHEA_AMR_INIT_REFINE_NONE:
  case RHEA_AMR_INIT_REFINE_UNIFORM:
  case RHEA_AMR_INIT_REFINE_HALF:
    break;
  case RHEA_AMR_INIT_REFINE_DEPTH:
    {
      rhea_amr_refine_depth_data_t *data = refine_data;

      RHEA_FREE (data->depth);
      RHEA_FREE (data);
    }
    break;
  default: /* unknown refinement type */
    RHEA_ABORT_NOT_REACHED ();
  }

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);

  /* return whether refinement was performed */
  return (refine_fn != NULL);
}

/******************************************************************************
 * AMR for ymir
 *****************************************************************************/

int
rhea_amr_flag_is_valid (const rhea_amr_flag_t amr_flag)
{
  switch (amr_flag) {
  case RHEA_AMR_FLAG_INIT:
    return 0;
  case RHEA_AMR_FLAG_COARSEN:
  case RHEA_AMR_FLAG_NO_CHANGE:
  case RHEA_AMR_FLAG_REFINE:
    return 1;
  default: /* unknown flag */
    RHEA_ABORT_NOT_REACHED ();
  }
}

static int
rhea_amr_coarsen_via_flag_fn (p4est_t *p4est, p4est_topidx_t tree,
                              p4est_quadrant_t *quadrants[])
{
  int                 k;

  SC_CHECK_ABORT (p4est_quadrant_is_familypv (quadrants), "Coarsen invocation");

  for (k = 0; k < P4EST_CHILDREN; k++) {
    rhea_p4est_quadrant_data_t *d = quadrants[k]->p.user_data;

    /* if at least one child is not flagged for coarsening */
    if (d->amr_flag != RHEA_AMR_FLAG_COARSEN) {
      return 0;
    }
  }

  /* if all of the children are flagged for coarsening */
  return 1;
}

static int
rhea_amr_refine_via_flag_fn (p4est_t *p4est, p4est_topidx_t tree,
                             p4est_quadrant_t *quadrant)
{
  rhea_p4est_quadrant_data_t *d = quadrant->p.user_data;

  /* if this quadrant is flagged for refinement */
  return (d->amr_flag == RHEA_AMR_FLAG_REFINE);
}

static p4est_gloidx_t
rhea_amr_get_global_n_elems (p4est_t *p4est)
{
  return p4est->global_num_quadrants;
}

static int
rhea_amr_get_level_max (p4est_t *p4est)
{
  sc_MPI_Comm         mpicomm = p4est->mpicomm;
  int                 mpiret;
  int                 level_max_loc;
  int                 level_max_glo;
  p4est_topidx_t      ti;

  /* find processor local max level among all trees */
  level_max_loc = 0;
  for (ti = p4est->first_local_tree; ti <= p4est->last_local_tree; ++ti) {
    p4est_tree_t       *tree = p4est_tree_array_index (p4est->trees, ti);

    level_max_loc = SC_MAX (level_max_loc, (int) tree->maxlevel);
  }

  /* get processor-global max level */
  level_max_glo = 0;
  mpiret = sc_MPI_Allreduce (&level_max_loc, &level_max_glo, 1, sc_MPI_INT,
                             sc_MPI_MAX, mpicomm); SC_CHECK_MPI (mpiret);

  return level_max_glo;
}

static int
rhea_amr_get_level_histogram (p4est_t *p4est,
                              int level_histogram[P4EST_QMAXLEVEL+1])
{
  sc_MPI_Comm         mpicomm = p4est->mpicomm;
  int                 mpiret;
  int                 level_histogram_loc[P4EST_QMAXLEVEL + 1];
  int                 l;
  p4est_topidx_t      ti;
  size_t              tqi;

  /* set counts to zero */
  for (l = 0; l <= P4EST_QMAXLEVEL; l++) {
    level_histogram_loc[l] = 0;
  }

  /* count levels on this processor */
  for (ti = p4est->first_local_tree; ti <= p4est->last_local_tree; ++ti) {
    p4est_tree_t       *t = p4est_tree_array_index (p4est->trees, ti);
    sc_array_t         *tquadrants = &(t->quadrants);

    for (tqi = 0; tqi < tquadrants->elem_count; ++tqi) {
      p4est_quadrant_t   *q = p4est_quadrant_array_index (tquadrants, tqi);

      RHEA_ASSERT (0 <= q->level && q->level <= P4EST_QMAXLEVEL);
      level_histogram_loc[q->level]++;
    }
  }

  /* sum levels from all processors */
  mpiret = sc_MPI_Allreduce (level_histogram_loc, level_histogram,
                             P4EST_QMAXLEVEL + 1, sc_MPI_INT, sc_MPI_SUM,
                             mpicomm); SC_CHECK_MPI (mpiret);

  /* find max level */
  for (l = P4EST_QMAXLEVEL; 0 <= l; l--) {
    if (0 < level_histogram[l]) {
      break;
    }
  }
  return l;
}

int
rhea_amr (p4est_t *p4est,
          const int n_cycles,
          const int n_cycles_uniform,
          const double flagged_elements_thresh_begin,
          const double flagged_elements_thresh_cycle,
          rhea_amr_flag_elements_fn_t flag_elements_fn,
          void *flag_elements_data,
          rhea_amr_data_initialize_fn_t data_initialize_fn,
          rhea_amr_data_finalize_fn_t data_finalize_fn,
          rhea_amr_data_project_fn_t data_project_fn,
          rhea_amr_data_partition_fn_t data_partition_fn,
          void *data)
{
  const double        n_elements_maxd = rhea_amr_n_elements_max;
  double              n_elements_estd;
  int                 reached_n_elements_max = 0;
  const p4est_gloidx_t  n_elements_begin = rhea_amr_get_global_n_elems (p4est);
  p4est_gloidx_t      n_elements_curr = n_elements_begin;
  p4est_gloidx_t      n_elements_prev, n_elements_coar, n_elements_refn;
  const int           iter_max = SC_MAX (1, n_cycles);
  const int           iter_uniform_max = SC_MAX (0, n_cycles_uniform);
  int                 iter, iter_uniform;
  p4est_gloidx_t      n_flagged_coar, n_flagged_refn;
  double              flagged_thresh, flagged_rel;
//double              coarsened_rel;
  double              refined_rel, changed_rel;
  int                 data_altered = 0;
  /* arguments for coarsening/refinement */
  const int           coarsen_recursively = 0;
  const int           refine_recursively = 0;
  p4est_coarsen_t     coarsen_fn = rhea_amr_coarsen_via_flag_fn;
  p4est_refine_t      refine_fn = rhea_amr_refine_via_flag_fn;
  p4est_init_t        init_fn = rhea_p4est_init_fn;
  /* arguments for partitioning */
  const int           partition_for_coarsening = 1;
  p4est_weight_t      partition_weight_fn = RHEA_AMR_P4EST_PARTITION_WEIGHT_FN;

  //TODO set flag `uniform` for if a uniform refinement function is detected;
  //     then avoid p4est calls to balance/partition

  /* start performance monitors */
  ymir_perf_counter_start (&rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_TOTAL]);

  RHEA_GLOBAL_INFOF_FN_BEGIN (
      __func__,
      "#cycles=%i, #elements_flagged_threshold_begin=%g, "
      "#elements_flagged_threshold_cycle=%g",
      n_cycles, flagged_elements_thresh_begin, flagged_elements_thresh_cycle);

  /* check input */
  RHEA_ASSERT (flag_elements_fn != NULL);

  for (iter = 0; iter < iter_max; iter++) { /* BEGIN: AMR iterations */
    RHEA_GLOBAL_INFOF_FN_TAG (__func__, "cycle=%i, #elements=%lli",
                              iter, (long long int) n_elements_curr);

    /*
     * Flag Elements
     */

    /* call function that flags elements for coarsening/refinement */
    ymir_perf_counter_start (&rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_FLAG]);
    flagged_rel = flag_elements_fn (p4est, flag_elements_data,
                                    &n_flagged_coar, &n_flagged_refn);
    ymir_perf_counter_stop_add (&rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_FLAG]);
    RHEA_GLOBAL_INFOF_FN_TAG (__func__, "cycle=%i, #elements_flagged_rel=%g, "
                              "n_flagged_coarsen=%lli, n_flagged_refine=%lli",
                              iter, flagged_rel,
                              (long long int) n_flagged_coar,
                              (long long int) n_flagged_refn);

    /* set threshold for #elements flagged for coarsening/refinement */
    if (isfinite (flagged_elements_thresh_begin) &&
        0.0 <= flagged_elements_thresh_begin) {
      flagged_thresh = flagged_elements_thresh_begin;
    }
    else {
      flagged_thresh = 0.0;
    }
    if (0 < iter && isfinite (flagged_elements_thresh_cycle) &&
        0.0 <= flagged_elements_thresh_cycle) {
      flagged_thresh = flagged_elements_thresh_cycle;
    }

    /* stop AMR if not enough elements were flagged */
    if (flagged_rel <= flagged_thresh) {
      RHEA_GLOBAL_INFOF_FN_TAG (
          __func__,
          "cycle=%i, stop_reason=\"#elements_flagged_rel %g <= threshold %g\"",
          iter, flagged_rel, flagged_thresh);
      break;
    }

    /*
     * Coarsen & Refine Elements
     */

    /* initialize data at first iteration */
    if (0 == iter && data_initialize_fn != NULL) {
      ymir_perf_counter_start (
          &rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_DATA_INIT]);
      data_initialize_fn (p4est, data);
      data_altered = 1;
      ymir_perf_counter_stop_add (
          &rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_DATA_INIT]);
    }

    /* update previous element count */
    n_elements_prev = n_elements_curr;

    /* coarsen p4est */
    ymir_perf_counter_start (&rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_COARSEN]);
    p4est_coarsen (p4est, coarsen_recursively, coarsen_fn, init_fn);
    ymir_perf_counter_stop_add (
        &rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_COARSEN]);
    n_elements_coar = rhea_amr_get_global_n_elems (p4est);
    RHEA_ASSERT (0 <= n_elements_prev - n_elements_coar);
  //coarsened_rel = (double) (n_elements_prev - n_elements_coar) /
  //                (double) n_elements_prev;

    /* check if max #elements is estimated to be reached */
    n_elements_estd =
      (double) (n_elements_coar + (P4EST_CHILDREN - 1) * n_flagged_refn);
    RHEA_GLOBAL_INFOF_FN_TAG (
        __func__, "cycle=%i, #elements_estimate=%g, #elements_max=%g",
        iter, n_elements_estd, n_elements_maxd);
    if (n_elements_maxd < n_elements_estd) {
      reached_n_elements_max = 1;
    }

    /* refine p4est */
    if (!reached_n_elements_max) {
      ymir_perf_counter_start (&rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_REFINE]);
      p4est_refine (p4est, refine_recursively, refine_fn, init_fn);
      ymir_perf_counter_stop_add (
          &rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_REFINE]);
    }
    n_elements_refn = rhea_amr_get_global_n_elems (p4est);
    RHEA_ASSERT (0 <= n_elements_refn - n_elements_coar);
    refined_rel = (double) (n_elements_refn - n_elements_coar) /
                  (double) n_elements_coar;

    /* balance p4est */
    ymir_perf_counter_start (&rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_BALANCE]);
    p4est_balance (p4est, P4EST_CONNECT_FULL, init_fn);
    ymir_perf_counter_stop_add (
        &rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_BALANCE]);
    n_elements_curr = rhea_amr_get_global_n_elems (p4est);
    RHEA_ASSERT (0 <= n_elements_curr - n_elements_refn);

    /* project data onto adapted mesh */
    if (data_project_fn != NULL) {
      ymir_perf_counter_start (
          &rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_DATA_PROJECT]);
      data_project_fn (p4est, data);
      ymir_perf_counter_stop_add (
          &rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_DATA_PROJECT]);
    }

    /*
     * Partition
     */

    /* partition p4est */
    ymir_perf_counter_start (&rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_PARTITION]);
    p4est_partition_ext (p4est, partition_for_coarsening, partition_weight_fn);
    ymir_perf_counter_stop_add (
        &rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_PARTITION]);

    /* partition data */
    if (data_partition_fn != NULL) {
      ymir_perf_counter_start (
          &rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_DATA_PARTITION]);
      data_partition_fn (p4est, data);
      ymir_perf_counter_stop_add (
          &rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_DATA_PARTITION]);
    }

    /*
     * Check Element Count Changes
     */

    /* stop AMR if max #elements was estimated to be reached */
    if (reached_n_elements_max) {
      RHEA_GLOBAL_INFOF_FN_TAG (
          __func__,
          "cycle=%i, stop_reason=\"#elements est. %g > #elements max %g\"",
          iter, n_elements_estd, n_elements_maxd);
      iter++; /* adjust #cycles */
      break;
    }

    /* stop AMR if not enough elements were actually changed due to
     * coarsening-balancing tradeoff, i.e., the following is low:
     *   (n_elements_refn - n_elements_coar) / n_elements_prev
     */
    changed_rel = fabs ((double) (n_elements_curr - n_elements_prev)) /
                  (double) n_elements_prev;
    if (changed_rel <= flagged_thresh && 2.0*refined_rel <= flagged_thresh) {
      RHEA_GLOBAL_INFOF_FN_TAG (
          __func__, "cycle=%i, stop_reason=\"rel. #elements changed %g and "
          "2x rel. #elements refined %g <= flagged threshold %g\"",
          iter, changed_rel, refined_rel, flagged_thresh);
      iter++; /* adjust #cycles */
      break;
    }
  } /* END: AMR iterations */

  /* print message if stopped due to reaching max #cycles */
  if (iter_max == iter) {
    RHEA_GLOBAL_INFOF_FN_TAG (
        __func__, "cycle=%i, stop_reason=\"max #cycles reached\"",
        iter - 1);
  }

  /* perform uniform refinement in addition to AMR */
  if (data_altered && 0 < iter_uniform_max) {
    RHEA_GLOBAL_INFOF_FN_TAG (__func__, "#cycles_uniform=%i",
                              iter_uniform_max);

    for (iter_uniform = 0; iter_uniform < iter_uniform_max; iter_uniform++) {
      RHEA_GLOBAL_INFOF_FN_TAG (
          __func__, "cycle_uniform=%i, #elements=%lli",
          iter_uniform, (long long int) n_elements_curr);

      /* check if max #elements is estimated to be reached */
      n_elements_estd = (double) (P4EST_CHILDREN * n_elements_curr);
      RHEA_GLOBAL_INFOF_FN_TAG (
          __func__, "cycle_uniform=%i, #elements_estimate=%g, #elements_max=%g",
          iter_uniform, n_elements_estd, n_elements_maxd);
      if (n_elements_maxd < n_elements_estd) {
        RHEA_GLOBAL_INFOF_FN_TAG (
            __func__,
            "cycle=%i, stop_reason=\"#elements est. %g > #elements max %g\"",
            iter_uniform, n_elements_estd, n_elements_maxd);
        break;
      }

      /* refine p4est */
      ymir_perf_counter_start (&rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_REFINE]);
      p4est_refine (p4est, refine_recursively, rhea_amr_refine_all_fn, init_fn);
      ymir_perf_counter_stop_add (
          &rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_REFINE]);

      /* skip balancing of p4est mesh */

      /* project data onto adapted mesh */
      if (data_project_fn != NULL) {
        ymir_perf_counter_start (
            &rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_DATA_PROJECT]);
        data_project_fn (p4est, data);
        ymir_perf_counter_stop_add (
            &rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_DATA_PROJECT]);
      }

      /* skip partitioning of p4est mesh */

      /* partition data */
      if (data_partition_fn != NULL) {
        ymir_perf_counter_start (
            &rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_DATA_PARTITION]);
        data_partition_fn (p4est, data);
        ymir_perf_counter_stop_add (
            &rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_DATA_PARTITION]);
      }

      /* udpate element counts */
      n_elements_prev = n_elements_curr;
      n_elements_curr = rhea_amr_get_global_n_elems (p4est);
    }
  }
  else {
    iter_uniform = 0;
  }

  /* finalize data if AMR was performed */
  if (data_altered && data_finalize_fn != NULL) {
    ymir_perf_counter_start (
        &rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_DATA_FINALIZE]);
    data_finalize_fn (p4est, data);
    ymir_perf_counter_stop_add (
        &rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_DATA_FINALIZE]);
  }

  /* print statistics */
  RHEA_GLOBAL_INFOF_FN_TAG (
      __func__, "#cycles=%i, #cycles_uniform=%i, "
      "#elements_beginning=%lli, #elements_end=%lli, "
      "change=%+.1f%%", iter, iter_uniform,
      (long long int) n_elements_begin, (long long int) n_elements_curr,
      (double) (n_elements_curr - n_elements_begin) /
      (double) n_elements_begin * 100.0);
  if (0 < iter || 0 < iter_uniform) {
    int                 level_max = -1;

    if (rhea_amr_log_level) {
      int                 level_histogram[P4EST_QMAXLEVEL + 1];
      int                 l;

      level_max = rhea_amr_get_level_histogram (p4est, level_histogram);
      RHEA_GLOBAL_INFOF ("<%s_mesh_level_histogram>\n", __func__);
      for (l = 0; l <= level_max; l++) {
        RHEA_GLOBAL_INFOF ("level %2i ... %9i elements\n",
                           l, level_histogram[l]);
      }
      RHEA_GLOBAL_INFOF ("</%s_mesh_level_histogram>\n", __func__);
    }
    else if (rhea_amr_log_level_max) {
      level_max = rhea_amr_get_level_max (p4est);
    }

    if (0 <= level_max) {
      RHEA_GLOBAL_INFOF_FN_TAG (__func__, "mesh_level_max=%i", level_max);
    }
  }

  RHEA_GLOBAL_INFO_FN_END (__func__);

  /* stop performance monitors */
  ymir_perf_counter_stop_add (&rhea_amr_perfmon[RHEA_AMR_PERFMON_AMR_TOTAL]);

  /* return number of performed AMR iterations */
  return iter;
}

/******************************************************************************
 * Generic Callback Functions for Coarsening/Refinement
 *****************************************************************************/

int
rhea_amr_coarsen_all_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                         p4est_quadrant_t * quadrant)
{
  return (0 < quadrant->level);
}

int
rhea_amr_refine_all_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                        p4est_quadrant_t * quadrant)
{
  return 1;
}

int
rhea_amr_coarsen_to_level_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                              p4est_quadrant_t * quadrant)
{
  int                *level_max = p4est->user_pointer;

  return (0 <= *level_max && *level_max < quadrant->level);
}

int
rhea_amr_refine_to_level_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                             p4est_quadrant_t * quadrant)
{
  int                *level_min = p4est->user_pointer;

  return (quadrant->level < *level_min);
}

int
rhea_amr_coarsen_half_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                          p4est_quadrant_t * quadrant)
{
  return (0 < quadrant->level) &&
         !rhea_amr_refine_half_fn (p4est, which_tree, quadrant);
}

int
rhea_amr_refine_half_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                         p4est_quadrant_t * quadrant)
{
  const int          *level_min = p4est->user_pointer;

  return (quadrant->level < *level_min) || (P4EST_ROOT_LEN / 2 <= quadrant->z);
}

int
rhea_amr_refine_depth_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                          p4est_quadrant_t * quadrant)
{
  rhea_amr_refine_depth_data_t  *d = p4est->user_pointer;
  const int          *depth = d->depth;
  const int           quad_depth = P4EST_ROOT_LEN - quadrant->z;
  int                 k;

  if ((int) quadrant->level < d->level_min) {
    return 1;
  }

  for (k = 0; k < d->count; k++) {
    if (quad_depth <= depth[k] && (int) quadrant->level < d->level_min+k+1) {
      return 1;
    }
  }

  return 0;
}

int
rhea_amr_refine_depth_box_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                              p4est_quadrant_t * quadrant)
{
  rhea_amr_refine_depth_data_t  *d = p4est->user_pointer;
  const int          *depth = d->depth;
  const double        zmax = d->p4est_zmax;
  const double        zscale = (double) P4EST_ROOT_LEN;
  double              xyz[3];
  int                 quad_depth;
  int                 k;

  RHEA_ASSERT (isfinite (zmax) && 0.0 <= zmax);

  if ((int) quadrant->level < d->level_min) {
    return 1;
  }

  p4est_qcoord_to_vertex (p4est->connectivity, which_tree, quadrant->x,
                          quadrant->y, quadrant->z, xyz);
  quad_depth = (int) round (zscale * (zmax - xyz[P4EST_DIM-1]));

  for (k = 0; k < d->count; k++) {
    if (quad_depth <= depth[k] && (int) quadrant->level < d->level_min+k+1) {
      return 1;
    }
  }

  return 0;
}

/******************************************************************************
 * Generic Flagging for Coarsening/Refinement
 *****************************************************************************/

double
rhea_amr_get_global_num_flagged (const p4est_locidx_t n_flagged_coar_loc,
                                 const p4est_locidx_t n_flagged_refn_loc,
                                 p4est_t *p4est,
                                 p4est_gloidx_t *n_flagged_coar_glo,
                                 p4est_gloidx_t *n_flagged_refn_glo)
{
  sc_MPI_Comm         mpicomm = p4est->mpicomm;
  int                 mpiret;
  int64_t             n_loc[2], n_glo[2];

  /* sum local contributions to get the global number of flagged quadrants */
  n_loc[0] = (int64_t) n_flagged_coar_loc;
  n_loc[1] = (int64_t) n_flagged_refn_loc;
  mpiret = sc_MPI_Allreduce (n_loc, n_glo, 2, MPI_INT64_T, sc_MPI_SUM,
                             mpicomm); SC_CHECK_MPI (mpiret);
  //TODO `sc_MPI_INT64_T` does not exist

  /* return relative number of flagged quadrants */
  if (n_flagged_coar_glo != NULL) {
    *n_flagged_coar_glo = n_glo[0];
  }
  if (n_flagged_refn_glo != NULL) {
    *n_flagged_refn_glo = n_glo[1];
  }
  return (double) (n_glo[0] + n_glo[1]) / (double) p4est->global_num_quadrants;
}

double
rhea_amr_flag_coarsen_half_fn (p4est_t *p4est, void *data,
                               p4est_gloidx_t *n_flagged_coar,
                               p4est_gloidx_t *n_flagged_refn)
{
  void               *user_pointer = p4est->user_pointer;
  p4est_locidx_t      n_flagged_loc;
  p4est_topidx_t      ti;
  size_t              tqi;

  p4est->user_pointer = data;

  /* flag quadrants */
  n_flagged_loc = 0;
  for (ti = p4est->first_local_tree; ti <= p4est->last_local_tree; ++ti) {
    p4est_tree_t       *t = p4est_tree_array_index (p4est->trees, ti);
    sc_array_t         *tquadrants = &(t->quadrants);

    for (tqi = 0; tqi < tquadrants->elem_count; ++tqi) {
      p4est_quadrant_t   *q = p4est_quadrant_array_index (tquadrants, tqi);
      rhea_p4est_quadrant_data_t *qd = q->p.user_data;

      if (rhea_amr_coarsen_half_fn (p4est, ti, q)) {
        qd->amr_flag = RHEA_AMR_FLAG_COARSEN;
        n_flagged_loc++;
      }
      else {
        qd->amr_flag = RHEA_AMR_FLAG_NO_CHANGE;
      }
    }
  }

  p4est->user_pointer = user_pointer;

  /* return absolute & relative numbers of flagged quadrants */
  return rhea_amr_get_global_num_flagged (n_flagged_loc, 0, p4est,
                                          n_flagged_coar, n_flagged_refn);
}

double
rhea_amr_flag_refine_half_fn (p4est_t *p4est, void *data,
                              p4est_gloidx_t *n_flagged_coar,
                              p4est_gloidx_t *n_flagged_refn)
{
  void               *user_pointer = p4est->user_pointer;
  p4est_locidx_t      n_flagged_loc;
  p4est_topidx_t      ti;
  size_t              tqi;

  p4est->user_pointer = data;

  /* flag quadrants */
  n_flagged_loc = 0;
  for (ti = p4est->first_local_tree; ti <= p4est->last_local_tree; ++ti) {
    p4est_tree_t       *t = p4est_tree_array_index (p4est->trees, ti);
    sc_array_t         *tquadrants = &(t->quadrants);

    for (tqi = 0; tqi < tquadrants->elem_count; ++tqi) {
      p4est_quadrant_t   *q = p4est_quadrant_array_index (tquadrants, tqi);
      rhea_p4est_quadrant_data_t *qd = q->p.user_data;

      if (rhea_amr_refine_half_fn (p4est, ti, q)) {
        qd->amr_flag = RHEA_AMR_FLAG_REFINE;
        n_flagged_loc++;
      }
      else {
        qd->amr_flag = RHEA_AMR_FLAG_NO_CHANGE;
      }
    }
  }

  p4est->user_pointer = user_pointer;

  /* return absolute & relative numbers of flagged quadrants */
  return rhea_amr_get_global_num_flagged (0, n_flagged_loc, p4est,
                                          n_flagged_coar, n_flagged_refn);
}

double
rhea_amr_flag_coarsen_to_level_fn (p4est_t *p4est, void *data,
                                   p4est_gloidx_t *n_flagged_coar,
                                   p4est_gloidx_t *n_flagged_refn)
{
  void               *user_pointer = p4est->user_pointer;
  p4est_locidx_t      n_flagged_loc;
  p4est_topidx_t      ti;
  size_t              tqi;

  p4est->user_pointer = data;

  /* flag quadrants */
  n_flagged_loc = 0;
  for (ti = p4est->first_local_tree; ti <= p4est->last_local_tree; ++ti) {
    p4est_tree_t       *t = p4est_tree_array_index (p4est->trees, ti);
    sc_array_t         *tquadrants = &(t->quadrants);

    for (tqi = 0; tqi < tquadrants->elem_count; ++tqi) {
      p4est_quadrant_t   *q = p4est_quadrant_array_index (tquadrants, tqi);
      rhea_p4est_quadrant_data_t *qd = q->p.user_data;

      if (rhea_amr_coarsen_to_level_fn (p4est, ti, q)) {
        qd->amr_flag = RHEA_AMR_FLAG_COARSEN;
        n_flagged_loc++;
      }
      else {
        qd->amr_flag = RHEA_AMR_FLAG_NO_CHANGE;
      }
    }
  }

  p4est->user_pointer = user_pointer;

  /* return absolute & relative numbers of flagged quadrants */
  return rhea_amr_get_global_num_flagged (n_flagged_loc, 0, p4est,
                                          n_flagged_coar, n_flagged_refn);
}

double
rhea_amr_flag_refine_to_level_fn (p4est_t *p4est, void *data,
                                  p4est_gloidx_t *n_flagged_coar,
                                  p4est_gloidx_t *n_flagged_refn)
{
  void               *user_pointer = p4est->user_pointer;
  p4est_locidx_t      n_flagged_loc;
  p4est_topidx_t      ti;
  size_t              tqi;

  p4est->user_pointer = data;

  /* flag quadrants */
  n_flagged_loc = 0;
  for (ti = p4est->first_local_tree; ti <= p4est->last_local_tree; ++ti) {
    p4est_tree_t       *t = p4est_tree_array_index (p4est->trees, ti);
    sc_array_t         *tquadrants = &(t->quadrants);

    for (tqi = 0; tqi < tquadrants->elem_count; ++tqi) {
      p4est_quadrant_t   *q = p4est_quadrant_array_index (tquadrants, tqi);
      rhea_p4est_quadrant_data_t *qd = q->p.user_data;

      if (rhea_amr_refine_to_level_fn (p4est, ti, q)) {
        qd->amr_flag = RHEA_AMR_FLAG_REFINE;
        n_flagged_loc++;
      }
      else {
        qd->amr_flag = RHEA_AMR_FLAG_NO_CHANGE;
      }
    }
  }

  p4est->user_pointer = user_pointer;

  /* return absolute & relative numbers of flagged quadrants */
  return rhea_amr_get_global_num_flagged (0, n_flagged_loc, p4est,
                                          n_flagged_coar, n_flagged_refn);
}
