/*
 */

#include <rhea_amr.h>
#include <rhea_base.h>
#include <rhea_discretization.h>
#include <p8est_bits.h>
#include <p8est_extended.h>

#define RHEA_AMR_P4EST_REPLACE_FN NULL
#define RHEA_AMR_P4EST_PARTITION_WEIGHT_FN NULL

/* default options */
#define RHEA_AMR_DEFAULT_INIT_REFINE_NAME "uniform"
#define RHEA_AMR_DEFAULT_INIT_REFINE_DEPTH_M NULL

char               *rhea_amr_init_refine_name =
  RHEA_AMR_DEFAULT_INIT_REFINE_NAME;
char               *rhea_amr_init_refine_depth_m =
  RHEA_AMR_DEFAULT_INIT_REFINE_DEPTH_M;

void
rhea_amr_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "AMR";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_S, "init-refine", '\0',
    &(rhea_amr_init_refine_name), RHEA_AMR_DEFAULT_INIT_REFINE_NAME,
    "Refinement of new/initial mesh: uniform, half, radial",
  YMIR_OPTIONS_S, "init-refine-depth", '\0',
    &(rhea_amr_init_refine_depth_m), RHEA_AMR_DEFAULT_INIT_REFINE_DEPTH_M,
    "Refinement 'depth': Sorted list of depths [m] (put shallowest at last)",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);
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
                      const int level_min,
                      const int level_max,
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

  /* get type of refinement from name */
  type = rhea_amr_init_refine_decode (refine_name);

  /* set refinement function and data */
  switch (type) {
  case RHEA_AMR_INIT_REFINE_UNKNOWN:
    RHEA_GLOBAL_LERROR ("Unknown refinement type");
    break;
  case RHEA_AMR_INIT_REFINE_NONE:
  case RHEA_AMR_INIT_REFINE_UNIFORM:
    /* no function necessary */
    break;
  case RHEA_AMR_INIT_REFINE_HALF:
    refine_fn = rhea_amr_refine_half_fn;
    break;
  case RHEA_AMR_INIT_REFINE_DEPTH:
    {
      double             *depth_m = NULL;
      int                 count, k;
      const double        rmin = domain_options->radius_min;
      const double        rmax = domain_options->radius_max;
      double              d, depth_rel;
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
      refine_data = data;
      refine_fn = rhea_amr_refine_depth_fn;
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

  /* return whether refinement was performed */
  return (refine_fn != NULL);
}

/******************************************************************************
 * AMR for ymir
 *****************************************************************************/

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

int
rhea_amr (p4est_t *p4est,
          const double n_flagged_elements_tol,
          const int amr_recursive_count,
          const double n_flagged_elements_recursive_tol,
          rhea_amr_flag_elements_fn_t flag_elements_fn,
          void *flag_elements_data,
          rhea_amr_data_initialize_fn_t data_initialize_fn,
          rhea_amr_data_finalize_fn_t data_finalize_fn,
          rhea_amr_data_project_fn_t data_project_fn,
          rhea_amr_data_partition_fn_t data_partition_fn,
          void *data)
{
  const char         *this_fn_name = "rhea_amr";
  /* arguments for coarsening/refinement */
  const int           iter_max = SC_MAX (1, amr_recursive_count);
  const int           coarsen_recursively = 0;
  const int           refine_recursively = 0;
  p4est_coarsen_t     coarsen_fn = rhea_amr_coarsen_via_flag_fn;
  p4est_refine_t      refine_fn = rhea_amr_refine_via_flag_fn;
  p4est_init_t        init_fn = rhea_p4est_init_fn;
  /* arguments for partitioning */
  const int           partition_for_coarsening = 1;
  p4est_weight_t      partition_weight_fn = RHEA_AMR_P4EST_PARTITION_WEIGHT_FN;

  double              n_flagged_rel;
  int                 iter;

  RHEA_GLOBAL_INFOF ("Into %s (#elements flagged tol. %g, recursive %i, "
                     "#elements flagged recursive tol. %g)\n",
                     this_fn_name, n_flagged_elements_tol, amr_recursive_count,
                     n_flagged_elements_recursive_tol);

  /* check input */
  RHEA_ASSERT (flag_elements_fn != NULL);

  for (iter = 0; iter < iter_max; iter++) { /* BEGIN: AMR iter */
    /*
     * Flag Elements
     */

    /* call function that flags elements for coarsening/refinement */
    n_flagged_rel = flag_elements_fn (p4est, flag_elements_data);

    /* print statistics */
    //TODO

    /* stop AMR if not enough elements were flagged */
    if (0 == iter) {
      if (isfinite (n_flagged_elements_tol) &&
          n_flagged_rel < n_flagged_elements_tol) {
        RHEA_GLOBAL_INFOF (
            "%s -- iter %i: Exit AMR, "
            "rel. #elements flagged %g, tolerance %g\n",
            this_fn_name, iter, n_flagged_rel, n_flagged_elements_tol);
        break;
      }
    }
    else {
      if (isfinite (n_flagged_elements_recursive_tol) &&
          n_flagged_rel < n_flagged_elements_recursive_tol) {
        RHEA_GLOBAL_INFOF (
            "%s -- iter %i: Stop AMR recursion, "
            "rel. #elements flagged %g, tolerance %g\n",
            this_fn_name, iter, n_flagged_rel,
            n_flagged_elements_recursive_tol);
        break;
      }
    }

    /*
     * Coarsen & Refine Elements
     */

    /* initialize data at first iteration */
    if (0 == iter && data_initialize_fn != NULL) {
      data_initialize_fn (p4est, data);
    }

    /* coarsen p4est */
    p4est_coarsen (p4est, coarsen_recursively, coarsen_fn, init_fn);

    /* print statistics */
    //TODO

    /* refine p4est */
    p4est_refine (p4est, refine_recursively, refine_fn, init_fn);

    /* print statistics */
    //TODO

    /* balance p4est */
    p4est_balance (p4est, P4EST_CONNECT_FULL, init_fn);

    /* print statistics */
    //TODO

    /* project data onto adapted mesh */
    if (data_project_fn != NULL) {
      data_project_fn (p4est, data);
    }

    /*
     * Partition
     */

    /* partition p4est */
    p4est_partition_ext (p4est, partition_for_coarsening, partition_weight_fn);

    /* print statistics */
    //TODO

    /* partition data */
    if (data_partition_fn != NULL) {
      data_partition_fn (p4est, data);
    }
  } /* END: AMR iter */

  /* finalize data if AMR was performed */
  if (0 < iter && data_finalize_fn != NULL) {
    data_finalize_fn (p4est, data);
  }

  /* print statistics */
  //TODO

  RHEA_GLOBAL_INFOF ("Done %s\n", this_fn_name);

  /* return number of performed AMR iterations */
  return iter;
}

/******************************************************************************
 * Generic Callback Functions for Coarsening/Refinement
 *****************************************************************************/

int
rhea_amr_refine_half_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                         p4est_quadrant_t * quadrant)
{
  if (P4EST_ROOT_LEN / 2 <= quadrant->z) {
    return 1;
  }
  else {
    return 0;
  }
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

/******************************************************************************
 * Generic Flagging for Coarsening/Refinement
 *****************************************************************************/

double
rhea_amr_flag_refine_half_fn (p4est_t *p4est, void *data)
{
  const p4est_locidx_t  n_quadrants = p4est->local_num_quadrants;
  p4est_locidx_t      n_flagged;
  p4est_topidx_t      ti;
  size_t              tqi;

  /* flag quadrants */
  n_flagged = 0;
  for (ti = p4est->first_local_tree; ti <= p4est->last_local_tree; ++ti) {
    p4est_tree_t       *t = p4est_tree_array_index (p4est->trees, ti);
    sc_array_t         *tquadrants = &(t->quadrants);

    for (tqi = 0; tqi < tquadrants->elem_count; ++tqi) {
      p4est_quadrant_t   *q = p4est_quadrant_array_index (tquadrants, tqi);
      rhea_p4est_quadrant_data_t *d = q->p.user_data;

      d->amr_flag = rhea_amr_refine_half_fn (p4est, ti, q);
      RHEA_ASSERT (d->amr_flag == RHEA_AMR_FLAG_REFINE ||
                   d->amr_flag == RHEA_AMR_FLAG_NO_CHANGE);
      if (d->amr_flag == RHEA_AMR_FLAG_REFINE) {
        n_flagged++;
      }
    }
  }

  /* return relative number of flagged quadrants */
  return ((double) n_flagged) / ((double) n_quadrants);
}
