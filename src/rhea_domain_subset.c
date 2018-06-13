#include <rhea_domain_subset.h>
#include <rhea_base.h>
#include <rhea_point_in_polygon.h>
#include <rhea_io_std.h>
#include <ymir_mesh.h>

/******************************************************************************
 * Column Subset
 *****************************************************************************/

int
rhea_domain_subset_point_is_in_column (const double x,
                                       const double y,
                                       const double z,
                                       rhea_domain_subset_column_t *column,
                                       rhea_domain_options_t *domain_options)
{
  const float         rmin = column->radius_min;
  const float         rmax = column->radius_max;
  float               test_r;
  int                 is_inside;

  /* check input */
  RHEA_ASSERT (column->radius_min <= column->radius_max);

  /* calculate radius */
  test_r = (float) rhea_domain_compute_radius (x, y, z, domain_options);

  /* check radial coordinate first, then lateral coordinates */
  if (rmin <= test_r && test_r <= rmax) { /* if radial coordinate matches */
    const float        *_sc_restrict vx = column->polygon_vertices_x;
    const float        *_sc_restrict vy = column->polygon_vertices_y;
    const size_t        nv = column->polygon_n_vertices;
    double              tmp_x, tmp_y;
    float               test_x, test_y;

    /* calculate (2-dim) lateral test coordinates */
    rhea_domain_extract_lateral (&tmp_x, &tmp_y, x, y, z,
                                 column->polygon_coord_type, domain_options);
    test_x = (float) tmp_x;
    test_y = (float) tmp_y;

    /* check if lateral coordinates are inside polygon */
    is_inside = rhea_point_in_polygon_is_inside (test_x, test_y, vx, vy, nv);
  }
  else { /* if radial coordinate is out of bounds */
    is_inside = 0;
  }

  /* return if point is inside column */
  return is_inside;
}

/******************************************************************************
 * Column Filter
 *****************************************************************************/

static ymir_gloidx_t
rhea_domain_subset_get_n_points_global (ymir_locidx_t n_points_in_subset_local,
                                        sc_MPI_Comm mpicomm)
{
  int                 mpiret;
  int64_t             n_loc = (int64_t) n_points_in_subset_local;
  int64_t             n_glo;

  mpiret = sc_MPI_Allreduce (&n_loc, &n_glo, 1, MPI_INT64_T, sc_MPI_SUM,
                             mpicomm); SC_CHECK_MPI (mpiret);
  return (ymir_gloidx_t) n_glo;
}

/** Filter parameters */
typedef struct rhea_domain_subset_filter_params
{
  /* domain and subset parameters */
  rhea_domain_options_t  *domain_options;
  rhea_domain_subset_column_t  *column;

  /* vector parameters */
  int                 n_fields;

  /* status */
  ymir_locidx_t       n_points_in_subset_local;
}
rhea_domain_subset_filter_params_t;

/** Filters at a node. */
static void
rhea_domain_subset_apply_filter_node_fn (double *v, double x, double y,
                                         double z, ymir_locidx_t nodeid,
                                         void *data)
{
  rhea_domain_subset_filter_params_t *params = data;

  if (!rhea_domain_subset_point_is_in_column (x, y, z, params->column,
                                              params->domain_options)) {
    const int           n_fields = params->n_fields;
    int                 fieldid;

    /* set outside values to zero */
    for (fieldid = 0; fieldid < n_fields; fieldid++) {
      v[fieldid] = 0.0;
    }
  }
  else {
    params->n_points_in_subset_local++;
  }
}

ymir_locidx_t
rhea_domain_subset_apply_filter (ymir_vec_t *vec,
                                 ymir_gloidx_t *n_points_in_subset_global,
                                 rhea_domain_options_t *domain_options,
                                 rhea_domain_subset_column_t *column)
{
  rhea_domain_subset_filter_params_t  params;

  /* set/init parameters */
  params.domain_options = domain_options;
  params.column = column;
  params.n_points_in_subset_local = 0;

  /* apply filter */
  if (ymir_vec_is_cvec (vec)) {
    params.n_fields = vec->ncfields;
    ymir_cvec_set_function (vec, rhea_domain_subset_apply_filter_node_fn,
                            &params);
  }
  else if (ymir_vec_is_dvec (vec)) {
    params.n_fields = vec->ndfields;
    ymir_dvec_set_function (vec, rhea_domain_subset_apply_filter_node_fn,
                            &params);
  }
  else {
    RHEA_ABORT_NOT_REACHED ();
  }

  /* count global #nodes inside subset */
  if (n_points_in_subset_global != NULL) {
    *n_points_in_subset_global = rhea_domain_subset_get_n_points_global (
        params.n_points_in_subset_local,
        ymir_mesh_get_MPI_Comm (ymir_vec_get_mesh (vec)));
  }

  /* return #nodes inside subset */
  return params.n_points_in_subset_local;
}

/******************************************************************************
 * Column I/O
 *****************************************************************************/

/** I/O parameters */
typedef struct rhea_domain_subset_io_params
{
  /* domain and subset parameters */
  rhea_domain_options_t  *domain_options;
  rhea_domain_subset_column_t  *column;
  rhea_domain_coordinate_type_t coord_type;

  /* vector parameters */
  int                 n_fields;

  /* coordinates & values */
  ymir_locidx_t       n_points_in_subset_local;
  ymir_locidx_t       pointid;
  float              *coordinate_value;
}
rhea_domain_subset_io_params_t;

/** Gathers values of one node. */
static void
rhea_domain_subset_io_gather_node_fn (double *v, double x, double y, double z,
                                      ymir_locidx_t nodeid, void *data)
{
  rhea_domain_subset_io_params_t *params = data;

  if (rhea_domain_subset_point_is_in_column (x, y, z, params->column,
                                             params->domain_options)) {
    const int           n_fields = params->n_fields;
    const int           n_entries_per_line = 3 + n_fields;
    const size_t        coord_idx = n_entries_per_line*params->pointid;
    const size_t        value_idx = coord_idx + 3;
    float              *coord = &(params->coordinate_value[coord_idx]);
    float              *value = &(params->coordinate_value[value_idx]);
    double              c0, c1, c2;
    int                 fieldid;

    /* set coordinates */
    rhea_domain_convert_coordinates (&c0, &c1, &c2, x, y, z,
                                     params->coord_type,
                                     params->domain_options);
    coord[0] = (float) c0;
    coord[1] = (float) c1;
    coord[2] = (float) c2;

    /* set values */
    for (fieldid = 0; fieldid < n_fields; fieldid++) {
      value[fieldid] = (float) v[fieldid];
    }

    /* increment array index */
    params->pointid++;
  }
}

ymir_locidx_t
rhea_domain_subset_write_txt (const char *file_path_base,
                              ymir_gloidx_t *n_points_in_subset_global,
                              ymir_vec_t *vec,
                              rhea_domain_options_t *domain_options,
                              rhea_domain_subset_column_t *column)
{
  rhea_domain_subset_io_params_t  params;
  ymir_vec_t         *dummy_vec = ymir_vec_template (vec);

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "path=\"%s\"", file_path_base);

  /* set parameters */
  params.domain_options = domain_options;
  params.column = column;
  params.coord_type = RHEA_DOMAIN_COORDINATE_SPHERICAL_GEO;
  params.n_points_in_subset_local = rhea_domain_subset_apply_filter (
      dummy_vec, NULL, domain_options, column);
  ymir_vec_destroy (dummy_vec);

  /* gather & write values */
  if (0 < params.n_points_in_subset_local) {
    ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vec);
    const int           mpisize = ymir_mesh_get_MPI_Comm_size (ymir_mesh);
    const int           mpirank = ymir_mesh_get_MPI_Comm_rank (ymir_mesh);
    const int           n_digits =
                          1 + (int) floor (log ((double) mpisize) / M_LN10);
    int                 n_entries_per_line;
    size_t              n_entries;
    char                path[BUFSIZ];

    RHEA_INFOF ("%s n_points_in_subset_local=%i\n", __func__,
                (int) params.n_points_in_subset_local);

    /* gather values */
    params.pointid = 0;
    if (ymir_vec_is_cvec (vec)) {
      params.n_fields = vec->ncfields;
      n_entries_per_line = 3 + params.n_fields;
      n_entries = (size_t) (n_entries_per_line*params.n_points_in_subset_local);
      params.coordinate_value = RHEA_ALLOC (float, n_entries);
      ymir_cvec_set_function (vec, rhea_domain_subset_io_gather_node_fn,
                              &params);
    }
    else if (ymir_vec_is_dvec (vec)) {
      params.n_fields = vec->ndfields;
      n_entries_per_line = 3 + params.n_fields;
      n_entries = (size_t) (n_entries_per_line*params.n_points_in_subset_local);
      params.coordinate_value = RHEA_ALLOC (float, n_entries);
      ymir_dvec_set_function (vec, rhea_domain_subset_io_gather_node_fn,
                              &params);
    }
    else {
      RHEA_ABORT_NOT_REACHED ();
    }

    /* write file */
    snprintf (path, BUFSIZ, "%s_%0*i.txt", file_path_base, n_digits, mpirank);
    rhea_io_std_write_float_to_txt (path, params.coordinate_value, n_entries,
                                    n_entries_per_line);

    /* destroy */
    RHEA_FREE (params.coordinate_value);
  }

  /* count global #nodes inside subset */
  if (n_points_in_subset_global != NULL) {
    *n_points_in_subset_global = rhea_domain_subset_get_n_points_global (
        params.n_points_in_subset_local,
        ymir_mesh_get_MPI_Comm (ymir_vec_get_mesh (vec)));
    RHEA_GLOBAL_INFOF_FN_END (__func__, "n_points_in_subset_global=%lli",
                              (long long int) *n_points_in_subset_global);
  }
  else {
    RHEA_GLOBAL_INFO_FN_END (__func__);
  }

  /* return #nodes inside subset */
  return params.n_points_in_subset_local;
}

/******************************************************************************
 * Specific Domain Subsets
 *****************************************************************************/

/**
 * Aleutian subduction:
 * - east longitude: 170E .. 220E
 * - latitude:       45N .. 65N (colatitude: 45..25)
 * - depth:          200km .. 20km
 */
static const size_t aleutian_n_vertices = 5;
static const float  aleutian_vertices_x[] = {170.0, 220.0, 220.0, 170.0, 170.0};
static const float  aleutian_vertices_y[] = { 45.0,  45.0,  25.0,  25.0,  45.0};
static const float  aleutian_depth_top_m = 20.0e3;
static const float  aleutian_depth_bottom_m = 200.0e3;

static rhea_domain_subset_column_t *
rhea_domain_subset_column_new_aleutian (rhea_domain_options_t *domain_options)
{
  rhea_domain_subset_column_t *column;

  /* create new column */
  column = RHEA_ALLOC (rhea_domain_subset_column_t, 1);

  /* set lateral parameters of column */
  column->polygon_vertices_x = RHEA_ALLOC (float, aleutian_n_vertices);
  column->polygon_vertices_y = RHEA_ALLOC (float, aleutian_n_vertices);
  memcpy (column->polygon_vertices_x, aleutian_vertices_x,
          aleutian_n_vertices * sizeof (float));
  memcpy (column->polygon_vertices_y, aleutian_vertices_y,
          aleutian_n_vertices * sizeof (float));
  column->polygon_n_vertices = aleutian_n_vertices;
  column->polygon_coord_type = RHEA_DOMAIN_COORDINATE_SPHERICAL_GEO_DIM;

  /* set radial parameters of column */
  column->radius_min = domain_options->radius_max -
                       aleutian_depth_bottom_m / domain_options->radius_max_m;
  column->radius_max = domain_options->radius_max -
                       aleutian_depth_top_m / domain_options->radius_max_m;

  return column;
}

static void
rhea_domain_subset_column_destroy_aleutian (rhea_domain_subset_column_t *column)
{
  RHEA_FREE (column->polygon_vertices_x);
  RHEA_FREE (column->polygon_vertices_y);
  RHEA_FREE (column);
}

ymir_locidx_t
rhea_domain_subset_apply_filter_aleutian (
                                      ymir_vec_t *vec,
                                      ymir_gloidx_t *n_points_in_subset_global,
                                      rhea_domain_options_t *domain_options)
{
  rhea_domain_subset_column_t *column;
  ymir_locidx_t       n_points_in_subset;

  /* create column */
  column = rhea_domain_subset_column_new_aleutian (domain_options);

  /* apply filter */
  n_points_in_subset = rhea_domain_subset_apply_filter (
      vec, n_points_in_subset_global, domain_options, column);

  /* destroy column */
  rhea_domain_subset_column_destroy_aleutian (column);

  /* return #nodes inside subset */
  return n_points_in_subset;
}

ymir_locidx_t
rhea_domain_subset_write_txt_aleutian (const char *file_path_base,
                                       ymir_vec_t *vec,
                                       ymir_gloidx_t *n_points_in_subset_global,
                                       rhea_domain_options_t *domain_options)
{
  rhea_domain_subset_column_t *column;
  ymir_locidx_t       n_points_in_subset;

  /* create column */
  column = rhea_domain_subset_column_new_aleutian (domain_options);

  /* write file */
  n_points_in_subset = rhea_domain_subset_write_txt (
      file_path_base, n_points_in_subset_global, vec, domain_options, column);

  /* destroy column */
  rhea_domain_subset_column_destroy_aleutian (column);

  /* return #nodes inside subset */
  return n_points_in_subset;
}
