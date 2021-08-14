#include <rhea_inversion_qoi.h>
#include <rhea_base.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_INVERSION_QOI_DEFAULT_STRESS_TYPE_LIST NULL

/* global options */
rhea_inversion_qoi_options_t rhea_inversion_qoi_options_global;

void
rhea_inversion_qoi_add_options (rhea_inversion_qoi_options_t *inv_qoi_options,
                                ymir_options_t *opt_sup)
{
  const char         *opt_prefix = "QuantitiesOfInterest";
  ymir_options_t     *opt = ymir_options_new ();
  rhea_inversion_qoi_options_t *inv_qoi_opt;

  /* set options storage */
  if (inv_qoi_options != NULL) {
    /* choose provided options */
    inv_qoi_opt = inv_qoi_options;
  }
  else {
    /* choose global options */
    inv_qoi_opt = &rhea_inversion_qoi_options_global;
  }

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_S, "stress-type-list", '\0',
    &(inv_qoi_opt->stress_type_list),
    RHEA_INVERSION_QOI_DEFAULT_STRESS_TYPE_LIST,
    "List of types of stress QOIs",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);
}

/******************************************************************************
 * Quantities Of Interest (QOI)
 *****************************************************************************/

/* inversion quantities of interest */
struct rhea_inversion_qoi
{
  /* stress inside mantle */
  rhea_inversion_obs_stress_t *stress_qoi_type;
  int                          stress_qoi_n;
  int                          stress_qoi_weak_n;

  /* Stokes problem (not owned) */
  rhea_stokes_problem_t  *stokes_problem;

  /* options (not owned) */
  rhea_weakzone_options_t  *weak_options;
};

rhea_inversion_qoi_t *
rhea_inversion_qoi_new (rhea_stokes_problem_t *stokes_problem,
                        rhea_inversion_qoi_options_t *inv_qoi_options)
{
  rhea_inversion_qoi_options_t *inv_qoi_opt;
  rhea_inversion_qoi_t         *inv_qoi;

  /* set options storage */
  if (inv_qoi_options != NULL) {
    /* choose provided options */
    inv_qoi_opt = inv_qoi_options;
  }
  else {
    /* choose global options */
    inv_qoi_opt = &rhea_inversion_qoi_options_global;
  }

  /* initialize QOI */
  inv_qoi = RHEA_ALLOC (rhea_inversion_qoi_t, 1);
  inv_qoi->stokes_problem = stokes_problem;
  inv_qoi->weak_options =
    rhea_stokes_problem_get_weakzone_options (stokes_problem);

  /* set stress QOI from options */
  if (NULL != inv_qoi_opt->stress_type_list &&
      rhea_weakzone_exists (inv_qoi->weak_options)) {
    double             *list = NULL;
    int                 n_entries, i;

    n_entries = ymir_options_convert_string_to_double (
        inv_qoi_opt->stress_type_list, &list);
    if (0 < n_entries) {
      inv_qoi->stress_qoi_weak_n = rhea_weakzone_get_total_n_labels (
          inv_qoi->weak_options);
      inv_qoi->stress_qoi_n      = n_entries;
      inv_qoi->stress_qoi_type   = RHEA_ALLOC (rhea_inversion_obs_stress_t,
                                               n_entries);
      for (i = 0; i < n_entries; i++) {
        inv_qoi->stress_qoi_type[i] =
          (rhea_inversion_obs_stress_t) ((int) list[i]);
      }
    }
    YMIR_FREE (list); /* was allocated in ymir */
  }
  else {
    inv_qoi->stress_qoi_type   = NULL;
    inv_qoi->stress_qoi_n      = 0;
    inv_qoi->stress_qoi_weak_n = 0;
  }

  /* return QOI */
  return inv_qoi;
}

void
rhea_inversion_qoi_destroy (rhea_inversion_qoi_t *inv_qoi)
{
  if (NULL != inv_qoi->stress_qoi_type) {
    RHEA_FREE (inv_qoi->stress_qoi_type);
  }
  RHEA_FREE (inv_qoi);
}

static int
rhea_inversion_qoi_stress_exists (rhea_inversion_qoi_t *inv_qoi)
{
  return (
      0 < inv_qoi->stress_qoi_n &&
      0 < inv_qoi->stress_qoi_weak_n &&
      inv_qoi->stress_qoi_type != NULL
  );
}

#if 0 //TODO unused
static int
rhea_inversion_qoi_exists (rhea_inversion_qoi_t *inv_qoi)
{
  return (
      rhea_inversion_qoi_stress_exists (inv_qoi)
  );
}
#endif

/******************************************************************************
 * Mapping QOI to Observation Operators
 *****************************************************************************/

static int
_stress_get_reduced_index (const int stress_qoi_index,
                           const int weakzone_label,
                           rhea_inversion_qoi_t *inv_qoi)
{
  int                 weakzone_index;

  weakzone_index = rhea_weakzone_lookup_index_from_label (
      weakzone_label, inv_qoi->weak_options);
  return stress_qoi_index * inv_qoi->stress_qoi_weak_n + weakzone_index;
}

static void
_stress_split_reduced_index (int *stress_qoi_index,
                             int *weakzone_label,
                             const int reduced_index,
                             rhea_inversion_qoi_t *inv_qoi)
{
  int                 weakzone_index;

  *stress_qoi_index = reduced_index / inv_qoi->stress_qoi_weak_n;
  weakzone_index    = reduced_index % inv_qoi->stress_qoi_weak_n;
  *weakzone_label   = (int) rhea_weakzone_lookup_label_from_index (
      weakzone_index, inv_qoi->weak_options);
}

static void
_stress_create_obs (rhea_inversion_obs_stress_t *stress_obs_type,
                    ymir_vec_t **stress_obs_weight,
                    rhea_inversion_qoi_t *inv_qoi,
                    const int stress_qoi_index,
                    const int weakzone_label)
{
  ymir_mesh_t        *ymir_mesh;

  /* exif if nothing to do */
  if (!rhea_inversion_qoi_stress_exists (inv_qoi)) {
    return;
  }

  /* check input */
  RHEA_ASSERT (0 <= stress_qoi_index && stress_qoi_index < inv_qoi->stress_qoi_n);

  /* set type of observation operator */
  *stress_obs_type = inv_qoi->stress_qoi_type[stress_qoi_index];

  /* create weight of observation operator */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (inv_qoi->stokes_problem);
  *stress_obs_weight = rhea_weakzone_new (ymir_mesh);
  rhea_weakzone_compute_indicator (
      *stress_obs_weight, weakzone_label, inv_qoi->weak_options);
}

void
rhea_inversion_qoi_stress_create_obs (
                                  rhea_inversion_obs_stress_t *stress_obs_type,
                                  ymir_vec_t **stress_obs_weight,
                                  rhea_inversion_qoi_t *inv_qoi,
                                  const int reduced_index)
{
  int                 stress_qoi_index, weakzone_label;

  _stress_split_reduced_index (&stress_qoi_index, &weakzone_label,
                               reduced_index, inv_qoi);
  _stress_create_obs (stress_obs_type, stress_obs_weight, inv_qoi,
                      stress_qoi_index, weakzone_label);
}

void
rhea_inversion_qoi_stress_clear_obs (
                                  ymir_vec_t *stress_obs_weight,
                                  rhea_inversion_qoi_t *inv_qoi)
{
  /* exif if nothing to do */
  if (!rhea_inversion_qoi_stress_exists (inv_qoi)) {
    return;
  }

  rhea_weakzone_destroy (stress_obs_weight);
}

/******************************************************************************
 * QOI Vector
 *****************************************************************************/

#if 0 //TODO unused
ymir_vec_t *
rhea_inversion_qoi_vec_new (rhea_inversion_qoi_t *inv_qoi)
{
  /* exif if nothing to do */
  if (!rhea_inversion_qoi_exists (inv_qoi)) {
    return NULL;
  }

  return ymir_vec_new_meshfree (rhea_inversion_qoi_get_count_all (inv_qoi));
}

void
rhea_inversion_qoi_vec_destroy (ymir_vec_t *vec)
{
  ymir_vec_destroy (vec);
}

int
rhea_inversion_qoi_vec_check_type (ymir_vec_t *vec,
                                   rhea_inversion_qoi_t *inv_qoi)
{
  return (
      ymir_vec_is_meshfree (vec) &&
      vec->n_meshfree == rhea_inversion_qoi_get_count_all (inv_qoi)
  );
}

int
rhea_inversion_qoi_vec_is_valid (ymir_vec_t *vec,
                                 sc_MPI_Comm mpicomm,
                                 rhea_inversion_qoi_t *inv_qoi)
{
  const int           n_qoi = rhea_inversion_qoi_get_count_all (inv_qoi);
  double             *v = vec->meshfree->e[0];
  int                 idx;

  /* check input */
  RHEA_ASSERT (rhea_inversion_qoi_vec_check_type (vec, inv_qoi));

  /* inspect processor-local values */
  for (idx = 0; idx < n_qoi; idx++) {
    if (!isfinite (v[idx])) {
      return 0;
    }
  }

  /* inspect values across processors (using a ring topology) */
  {
    int                 mpisize, mpirank, mpiret;
    int                 mpirank_prev, mpirank_next;
    const int           THIS_MPI_TAG = 1;
    sc_MPI_Request      request_recv, request_send;
    sc_MPI_Status       status;
    double             *v_prev;

    /* get parallel environment */
    mpiret = sc_MPI_Comm_size (mpicomm, &mpisize); SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);

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
    v_prev = RHEA_ALLOC (double, n_qoi);
    mpiret = sc_MPI_Irecv (v_prev, n_qoi, sc_MPI_DOUBLE, mpirank_prev,
                           THIS_MPI_TAG, mpicomm, &request_recv);
    SC_CHECK_MPI (mpiret);

    /* send values to next neighbor */
    mpiret = sc_MPI_Isend (v, n_qoi, sc_MPI_DOUBLE, mpirank_next,
                           THIS_MPI_TAG, mpicomm, &request_send);
    SC_CHECK_MPI (mpiret);

    /* wait until receive is complete */
    mpiret = sc_MPI_Waitall (1, &request_recv, &status);
    SC_CHECK_MPI (mpiret);

    /* compare values of previous neighbor */
    for (idx = 0; idx < n_qoi; idx++) {
      const double        rel_err = fabs (v[idx] - v_prev[idx]) / fabs (v[idx]);

      if (SC_1000_EPS < rel_err) {
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
rhea_inversion_qoi_vec_print (ymir_vec_t *qoi_vec,
                              rhea_inversion_qoi_t *inv_qoi)
{
  const int           n_qoi = rhea_inversion_qoi_get_count_all (inv_qoi);
  const double       *qoi = qoi_vec->meshfree->e[0];
  int                 idx;

  /* check input */
  RHEA_ASSERT (rhea_inversion_qoi_vec_check_type (qoi_vec, inv_qoi));

  /* print each QOI */
  for (idx = 0; idx < n_qoi; idx++) {
    RHEA_GLOBAL_INFOF ("qoi# %3i: %g\n", idx, qoi[idx]);
  }
}

void
rhea_inversion_qoi_vec_reduced_copy (
                                  ymir_vec_t *vec,
                                  sc_dmatrix_t *vec_reduced,
                                  rhea_inversion_qoi_t *inv_qoi,
                                  const rhea_inversion_qoi_class_t class_id)
{
  const double       *r = vec_reduced->e[0];
  double             *v = vec->meshfree->e[0];
  int                 offset = 0;
  int                 nrows  = 0;
  int                 ncols  = 0;
  int                 i, j, idx_reduced;

  /* get offset and size */
  switch (class_id) {
  case RHEA_INVERSION_QOI_STRESS:
    offset = 0;
    nrows  = inv_qoi->stress_qoi_n;
    ncols  = inv_qoi->stress_qoi_weak_n;
    break;
  default: /* unknown QOI class */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* check input */
  RHEA_ASSERT (rhea_inversion_qoi_vec_check_type (vec, inv_qoi));
  RHEA_ASSERT (vec_reduced->m == nrows);
  RHEA_ASSERT (vec_reduced->n == ncols);

  /* fill entries */
  idx_reduced = 0;
  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      RHEA_ASSERT (idx_reduced < vec_reduced->m * vec_reduced->n);
      v[offset + i*ncols + j] = r[idx_reduced];
      idx_reduced++;
    }
  }
}
#endif

sc_dmatrix_t *
rhea_inversion_qoi_stress_vec_reduced_new (rhea_inversion_qoi_t *inv_qoi)
{
  const int           n_qoi = inv_qoi->stress_qoi_n * inv_qoi->stress_qoi_weak_n;

  /* exif if nothing to do */
  if (!rhea_inversion_qoi_stress_exists (inv_qoi)) {
    return NULL;
  }

  /* create reduced QOI vector */
  return sc_dmatrix_new (n_qoi, 1);
}

void
rhea_inversion_qoi_stress_vec_reduced_destroy (sc_dmatrix_t *vec_reduced)
{
  sc_dmatrix_destroy (vec_reduced);
}

int
rhea_inversion_qoi_stress_vec_set_entry (
                                  sc_dmatrix_t *vec_reduced,
                                  const double val,
                                  rhea_inversion_qoi_t *inv_qoi,
                                  const int stress_qoi_index,
                                  const int weakzone_label)
{
  const int           idx_reduced =
    _stress_get_reduced_index (stress_qoi_index, weakzone_label, inv_qoi);

  /* set value and return index */
  vec_reduced->e[idx_reduced][0] = val;
  return idx_reduced;
}

sc_dmatrix_t *
rhea_inversion_qoi_stress_to_param_matrix_new (
                                  const int n_parameters,
                                  rhea_inversion_qoi_t *inv_qoi)
{
  const int           n_qoi = inv_qoi->stress_qoi_n * inv_qoi->stress_qoi_weak_n;

  /* exif if nothing to do */
  if (!rhea_inversion_qoi_stress_exists (inv_qoi)) {
    return NULL;
  }

  /* create matrix mapping QOI to parameters */
  return sc_dmatrix_new (n_parameters, n_qoi);
}

void
rhea_inversion_qoi_stress_to_param_matrix_destroy (sc_dmatrix_t *mat)
{
  sc_dmatrix_destroy (mat);
}

int
rhea_inversion_qoi_stress_to_param_matrix_set_column (
                                  sc_dmatrix_t *mat,
                                  sc_dmatrix_t *col,
                                  rhea_inversion_qoi_t *inv_qoi,
                                  const int stress_qoi_index,
                                  const int weakzone_label)
{
  const int           idx_reduced =
    _stress_get_reduced_index (stress_qoi_index, weakzone_label, inv_qoi);
  int                 idx_param;

  /* check input */
  RHEA_ASSERT (1 == col->n);
  RHEA_ASSERT (mat->m == col->m);
  RHEA_ASSERT (mat->n == rhea_inversion_qoi_get_count_all (inv_qoi));
  RHEA_ASSERT (0 <= stress_qoi_index && stress_qoi_index < inv_qoi->stress_qoi_n);
  RHEA_ASSERT (0 <= idx_reduced && idx_reduced < mat->n);

  /* fill one column */
  for (idx_param = 0; idx_param < mat->m; idx_param++) {
    mat->e[idx_param][idx_reduced] = col->e[idx_param][0];
  }

  /* return index */
  return idx_reduced;
}

/******************************************************************************
 * Data Access
 *****************************************************************************/

int
rhea_inversion_qoi_get_count_all (rhea_inversion_qoi_t *inv_qoi)
{
  return inv_qoi->stress_qoi_n * inv_qoi->stress_qoi_weak_n;
}
