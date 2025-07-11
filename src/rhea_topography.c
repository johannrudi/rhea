/*
 */

#include <rhea_topography.h>
#include <rhea_base.h>
#include <rhea_domain.h>
#include <rhea_io_mpi.h>
#include <rhea_io_std.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_TOPOGRAPHY_DEFAULT_TYPE_NAME "NONE"
#define RHEA_TOPOGRAPHY_DEFAULT_POINTS_FILE_PATH_BIN NULL
#define RHEA_TOPOGRAPHY_DEFAULT_POINTS_FILE_PATH_TXT NULL
#define RHEA_TOPOGRAPHY_DEFAULT_DISPLACEMENTS_FILE_PATH_BIN NULL
#define RHEA_TOPOGRAPHY_DEFAULT_DISPLACEMENTS_FILE_PATH_TXT NULL
#define RHEA_TOPOGRAPHY_DEFAULT_LABELS_FILE_PATH_BIN NULL
#define RHEA_TOPOGRAPHY_DEFAULT_LABELS_FILE_PATH_TXT NULL
#define RHEA_TOPOGRAPHY_DEFAULT_N_POINTS (5000000)
#define RHEA_TOPOGRAPHY_DEFAULT_WRITE_POINTS_FILE_PATH_BIN NULL
#define RHEA_TOPOGRAPHY_DEFAULT_WRITE_LABELS_FILE_PATH_BIN NULL
#define RHEA_TOPOGRAPHY_DEFAULT_WRITE_DISPLACEMENTS_FILE_PATH_BIN NULL

/* initialize options */
char               *rhea_topography_type_name =
  RHEA_TOPOGRAPHY_DEFAULT_TYPE_NAME;
char               *rhea_topography_points_file_path_bin =
  RHEA_TOPOGRAPHY_DEFAULT_POINTS_FILE_PATH_BIN;
char               *rhea_topography_points_file_path_txt =
  RHEA_TOPOGRAPHY_DEFAULT_POINTS_FILE_PATH_TXT;
char               *rhea_topography_displacements_file_path_bin =
  RHEA_TOPOGRAPHY_DEFAULT_DISPLACEMENTS_FILE_PATH_BIN;
char               *rhea_topography_displacements_file_path_txt =
  RHEA_TOPOGRAPHY_DEFAULT_DISPLACEMENTS_FILE_PATH_TXT;
char               *rhea_topography_labels_file_path_bin =
  RHEA_TOPOGRAPHY_DEFAULT_LABELS_FILE_PATH_BIN;
char               *rhea_topography_labels_file_path_txt =
  RHEA_TOPOGRAPHY_DEFAULT_LABELS_FILE_PATH_TXT;
int                 rhea_topography_n_points =
  RHEA_TOPOGRAPHY_DEFAULT_N_POINTS;
char               *rhea_topography_write_points_file_path_bin =
  RHEA_TOPOGRAPHY_DEFAULT_WRITE_POINTS_FILE_PATH_BIN;
char               *rhea_topography_write_labels_file_path_bin =
  RHEA_TOPOGRAPHY_DEFAULT_WRITE_LABELS_FILE_PATH_BIN;
char               *rhea_topography_write_displacements_file_path_bin =
  RHEA_TOPOGRAPHY_DEFAULT_WRITE_DISPLACEMENTS_FILE_PATH_BIN;

void
rhea_topography_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Topography";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_S, "type", '\0',
    &(rhea_topography_type_name), RHEA_TOPOGRAPHY_DEFAULT_TYPE_NAME,
    "Topography type name: data_points_displacements, "
    "data_points_displacements_labels",

  YMIR_OPTIONS_S, "points-file-path-bin", '\0',
    &(rhea_topography_points_file_path_bin),
    RHEA_TOPOGRAPHY_DEFAULT_POINTS_FILE_PATH_BIN,
    "Path to a binary file with (x,y,z) coordinates of topography points",
  YMIR_OPTIONS_S, "points-file-path-txt", '\0',
    &(rhea_topography_points_file_path_txt),
    RHEA_TOPOGRAPHY_DEFAULT_POINTS_FILE_PATH_TXT,
    "Path to a text file with (x,y,z) coordinates of topography points",
  YMIR_OPTIONS_S, "displacements-file-path-bin", '\0',
    &(rhea_topography_displacements_file_path_bin),
    RHEA_TOPOGRAPHY_DEFAULT_DISPLACEMENTS_FILE_PATH_BIN,
    "Path to a binary file with displacements of topography points",
  YMIR_OPTIONS_S, "displacements-file-path-txt", '\0',
    &(rhea_topography_displacements_file_path_txt),
    RHEA_TOPOGRAPHY_DEFAULT_DISPLACEMENTS_FILE_PATH_TXT,
    "Path to a text file with displacements of topography points",
  YMIR_OPTIONS_S, "labels-file-path-bin", '\0',
    &(rhea_topography_labels_file_path_bin),
    RHEA_TOPOGRAPHY_DEFAULT_LABELS_FILE_PATH_BIN,
    "Path to a binary file with labels of topography points",
  YMIR_OPTIONS_S, "labels-file-path-txt", '\0',
    &(rhea_topography_labels_file_path_txt),
    RHEA_TOPOGRAPHY_DEFAULT_LABELS_FILE_PATH_TXT,
    "Path to a text file with labels of topography points",

  YMIR_OPTIONS_I, "num-points", '\0',
    &(rhea_topography_n_points), RHEA_TOPOGRAPHY_DEFAULT_N_POINTS,
    "Number of points that are imported",

  YMIR_OPTIONS_S, "write-points-file-path-bin", '\0',
    &(rhea_topography_write_points_file_path_bin),
    RHEA_TOPOGRAPHY_DEFAULT_WRITE_POINTS_FILE_PATH_BIN,
    "Output path for a binary file with coordinates of topography points",
  YMIR_OPTIONS_S, "write-displacements-file-path-bin", '\0',
    &(rhea_topography_write_displacements_file_path_bin),
    RHEA_TOPOGRAPHY_DEFAULT_WRITE_DISPLACEMENTS_FILE_PATH_BIN,
    "Output path for a binary file with displacements of topography points",
  YMIR_OPTIONS_S, "write-labels-file-path-bin", '\0',
    &(rhea_topography_write_labels_file_path_bin),
    RHEA_TOPOGRAPHY_DEFAULT_WRITE_LABELS_FILE_PATH_BIN,
    "Output path for a binary file with labels of topography points",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);
}

void
rhea_topography_process_options (rhea_topography_options_t *opt,
                               rhea_domain_options_t *domain_options)
{
  /* set topography type */
  if (strcmp (rhea_topography_type_name, "NONE") == 0) {
    opt->type = RHEA_TOPOGRAPHY_NONE;
  }
  else if (strcmp (rhea_topography_type_name,
                   "data_points_displacements") == 0) {
    opt->type = RHEA_TOPOGRAPHY_DATA_POINTS_DISPLS;
  }
  else if (strcmp (rhea_topography_type_name,
                   "data_points_displacements_labels") == 0)
  {
    opt->type = RHEA_TOPOGRAPHY_DATA_POINTS_DISPLS_LABELS;
  }
  else { /* unknown type name */
    RHEA_ABORT ("Unknown topography type name");
  }

  /* set paths to binary & text files */
  opt->points_file_path_bin = rhea_topography_points_file_path_bin;
  opt->points_file_path_txt = rhea_topography_points_file_path_txt;
  opt->labels_file_path_bin = rhea_topography_labels_file_path_bin;
  opt->labels_file_path_txt = rhea_topography_labels_file_path_txt;
  opt->displacements_file_path_bin =
    rhea_topography_displacements_file_path_bin;
  opt->displacements_file_path_txt =
    rhea_topography_displacements_file_path_txt;

  /* set number of points */
  opt->n_points = rhea_topography_n_points;

  /* set output paths */
  opt->write_points_file_path_bin = rhea_topography_write_points_file_path_bin;
  opt->write_labels_file_path_bin = rhea_topography_write_labels_file_path_bin;
  opt->write_displacements_file_path_bin =
    rhea_topography_write_displacements_file_path_bin;

  /* init data */
  opt->pointcloud = NULL;

  /* set dependent options */
  opt->domain_options = domain_options;
}

int
rhea_topography_exists (rhea_topography_options_t *opt)
{
  switch (opt->type) {
  case RHEA_TOPOGRAPHY_NONE:
    return 0;
  case RHEA_TOPOGRAPHY_DATA_POINTS_DISPLS:
  case RHEA_TOPOGRAPHY_DATA_POINTS_DISPLS_LABELS:
    return 1;
  default: /* unknown topography type */
    RHEA_ABORT_NOT_REACHED ();
    return 0;
  }
}

/******************************************************************************
 * Topography Vector
 *****************************************************************************/

ymir_vec_t *
rhea_topography_new (ymir_mesh_t *ymir_mesh)
{
  return ymir_face_cvec_new (ymir_mesh, RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
}

void
rhea_topography_destroy (ymir_vec_t *topography)
{
  ymir_vec_destroy (topography);
}

int
rhea_topography_check_vec_type (ymir_vec_t *vec)
{
  return (
      ymir_vec_is_cvec (vec) &&
      vec->ncfields == 1 &&
      vec->meshnum == RHEA_DOMAIN_BOUNDARY_FACE_TOP
  );
}

int
rhea_topography_is_valid (ymir_vec_t *vec)
{
  return sc_dmatrix_is_valid (vec->dataown) && sc_dmatrix_is_valid (vec->coff);
}

/******************************************************************************
 * Data
 *****************************************************************************/

void
rhea_topography_data_create (rhea_topography_options_t *opt,
                             sc_MPI_Comm mpicomm)
{
  rhea_domain_options_t *domain_options = opt->domain_options;
  int                 n_read;
  int                 create_coordinates;
  int                 create_displacements;
  int                 create_labels;
  double             *coordinates = NULL;
  double             *displacements = NULL;
  int                *labels = NULL;

  /* init flags */
  if (opt->points_file_path_bin != NULL ||
      opt->points_file_path_txt != NULL) {
    create_coordinates = 1;
  }
  else {
    create_coordinates = 0;
  }
  if (opt->displacements_file_path_bin != NULL ||
      opt->displacements_file_path_txt != NULL) {
    create_displacements = 1;
  }
  else {
    create_displacements = 0;
  }
  if (opt->labels_file_path_bin != NULL ||
      opt->labels_file_path_txt != NULL) {
    create_labels = 1;
  }
  else {
    create_labels = 0;
  }

  /* decide what data needs to be created */
  switch (opt->type) {
  case RHEA_TOPOGRAPHY_NONE:
    return;
  case RHEA_TOPOGRAPHY_DATA_POINTS_DISPLS:
    RHEA_CHECK_ABORT (create_coordinates,
                      "file path missing: topography coordinates");
    RHEA_CHECK_ABORT (create_displacements,
                      "file path missing: topography displacements");
    create_labels = 0;
    break;
  case RHEA_TOPOGRAPHY_DATA_POINTS_DISPLS_LABELS:
    RHEA_CHECK_ABORT (create_coordinates,
                      "file path missing: topography coordinates");
    RHEA_CHECK_ABORT (create_displacements,
                      "file path missing: topography displacements");
    RHEA_CHECK_ABORT (create_labels,
                      "file path missing: topography labels");
    break;
  default: /* unknown topography type */
    RHEA_ABORT_NOT_REACHED ();
  }

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "type=%i", opt->type);

  /* check input */
  RHEA_ASSERT (0 < opt->n_points);

  /* create, read, and write coordinates (updates #points in options) */
  if (create_coordinates) {
    const int           n_entries = 3 * opt->n_points;
    const char         *file_path_bin = opt->points_file_path_bin;
    const char         *file_path_txt = opt->points_file_path_txt;
    const char         *write_file_path_bin = opt->write_points_file_path_bin;

    coordinates = RHEA_ALLOC (double, n_entries);

    /* read from file */
    if (file_path_bin != NULL) { /* if read from binary file */
      n_read = rhea_io_mpi_read_broadcast_double (
          coordinates, n_entries, file_path_bin, NULL /* path txt */, mpicomm);
      RHEA_ASSERT (n_read == n_entries);
    }
    else if (file_path_txt != NULL) { /* if read from text file */
      n_read = rhea_io_mpi_read_broadcast_double (
          coordinates, 0 /* #entries */, write_file_path_bin, file_path_txt,
          mpicomm);
      RHEA_ASSERT (0 < n_read && n_read <= n_entries);
      RHEA_ASSERT (!(n_read % 3));
      opt->n_points = n_read / 3;
    }
    else { /* otherwise no reading possible */
      RHEA_ABORT_NOT_REACHED ();
    }

    RHEA_GLOBAL_INFOF ("%s Number of topography points=%i\n", __func__,
                       opt->n_points);
  }

  /* create, read, and write displacements */
  if (create_displacements) {
    const int           n_entries = opt->n_points;
    const char         *file_path_bin = opt->displacements_file_path_bin;
    const char         *file_path_txt = opt->displacements_file_path_txt;
    const char         *write_file_path_bin =
                          opt->write_displacements_file_path_bin;

    displacements = RHEA_ALLOC (double, n_entries);

    /* read from file */
    if (file_path_bin != NULL) { /* if read from binary file */
      n_read = rhea_io_mpi_read_broadcast_double (
          displacements, n_entries, file_path_bin, NULL /* path txt */,
          mpicomm);
      RHEA_ASSERT (n_read == n_entries);
    }
    else if (file_path_txt != NULL) { /* if read from text file */
      n_read = rhea_io_mpi_read_broadcast_double (
          displacements, 0 /* #entries */, write_file_path_bin, file_path_txt,
          mpicomm);
      RHEA_ASSERT (0 < n_read && n_read <= n_entries);
    }
    else { /* otherwise no reading possible */
      RHEA_ABORT_NOT_REACHED ();
    }
  }

  /* create, read, and write labels */
  if (create_labels) {
    const int           n_entries = opt->n_points;
    const char         *file_path_bin = opt->labels_file_path_bin;
    const char         *file_path_txt = opt->labels_file_path_txt;
    const char         *write_file_path_bin = opt->write_labels_file_path_bin;

    labels = RHEA_ALLOC (int, n_entries);

    /* read from file */
    if (file_path_bin != NULL) { /* if read from binary file */
      n_read = rhea_io_mpi_read_broadcast_int (
          labels, n_entries, file_path_bin, NULL /* path txt */, mpicomm);
      RHEA_ASSERT (n_read == n_entries);
    }
    else if (file_path_txt != NULL) { /* if read from text file */
      n_read = rhea_io_mpi_read_broadcast_int (
          labels, 0 /* #entries */, write_file_path_bin, file_path_txt,
          mpicomm);
      RHEA_ASSERT (0 < n_read && n_read <= n_entries);
    }
    else { /* otherwise no reading possible */
      RHEA_ABORT_NOT_REACHED ();
    }
  }

  /* create point cloud */
  opt->pointcloud = rhea_pointcloud_topography_new (
      domain_options->x_min, domain_options->x_max,
      domain_options->y_min, domain_options->y_max,
      domain_options->z_min, domain_options->z_max,
      coordinates, opt->n_points);
  if (create_displacements) {
    rhea_pointcloud_topography_set_displacements (opt->pointcloud,
                                                  displacements);
  }
  if (create_labels) {
    rhea_pointcloud_topography_set_labels (opt->pointcloud, labels);
  }

  /* destroy */
  RHEA_FREE (coordinates);
  if (create_displacements) {
    RHEA_FREE (displacements);
  }
  if (create_labels) {
    RHEA_FREE (labels);
  }

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

void
rhea_topography_data_clear (rhea_topography_options_t *opt)
{
  if (opt->pointcloud != NULL) {
    rhea_pointcloud_topography_destroy (opt->pointcloud);
  }
}

/******************************************************************************
 * Topography Computation
 *****************************************************************************/

static double
rhea_topography_displacement_integrate (
                                      const double *_sc_restrict nearest_dist,
                                      const double *_sc_restrict nearest_displ,
                                      const int n_nearest)
{
  double              displ;
  int                 k;

  /* compute average */ //TODO better compute via interpolation
  displ = 0.0;
  for (k = 0; k < n_nearest; k++) {
    displ += nearest_displ[k];
  }
  displ *= 1.0 / (double) n_nearest;

  return displ;
}

#define RHEA_TOPOGRAPHY_DISPLACEMENT_N_NEAREST (1) //TODO increase

double
rhea_topography_displacement_node (int *nearest_label,
                                   double x, double y, double z,
                                   rhea_topography_options_t *opt)
{
  double              nearest_dist[RHEA_TOPOGRAPHY_DISPLACEMENT_N_NEAREST];
  double              nearest_displ[RHEA_TOPOGRAPHY_DISPLACEMENT_N_NEAREST];
  double              displ;
  int                 n_found;

  /* check input */
  RHEA_ASSERT (opt->pointcloud != NULL);

  /* project target coordinates onto the surface */
  rhea_domain_project_to_surface (&x, &y, &z, opt->domain_options);

  /* find `n` nearest points from the cloud */
  {
    const double        pt[3] = {x, y, z};
    const int           n_nearest = RHEA_TOPOGRAPHY_DISPLACEMENT_N_NEAREST;
#ifdef RHEA_ENABLE_DEBUG
    int                 k;
#endif

    switch (opt->type) {
    case RHEA_TOPOGRAPHY_DATA_POINTS_DISPLS:
      n_found = rhea_pointcloud_topography_find_n_nearest (
          nearest_dist, NULL /* coord's */, nearest_displ, NULL /* labels */,
          n_nearest, opt->pointcloud, pt);
      break;
    case RHEA_TOPOGRAPHY_DATA_POINTS_DISPLS_LABELS:
      RHEA_ASSERT (nearest_label != NULL);
      n_found = rhea_pointcloud_topography_find_n_nearest (
          nearest_dist, NULL /* coord's */, nearest_displ, nearest_label,
          n_nearest, opt->pointcloud, pt);
      break;
    default: /* unknown topography type */
      RHEA_ABORT_NOT_REACHED ();
    }
    RHEA_ASSERT (n_nearest == n_found);
#ifdef RHEA_ENABLE_DEBUG
    for (k = 0; k < n_found; k++) {
      RHEA_ASSERT (isfinite (nearest_displ[k]));
    }
#endif
  }

  /* calculate displacement for target coordinates */
  displ = rhea_topography_displacement_integrate (nearest_dist, nearest_displ,
                                                  n_found);
  RHEA_ASSERT (isfinite (displ));
  return displ;
}

static void
rhea_topography_displacement_node_vol_fn (double *displ, double x, double y,
                                          double z, ymir_locidx_t nid,
                                          void *data)
{
  rhea_topography_options_t *opt = data;
  int                 label[RHEA_TOPOGRAPHY_DISPLACEMENT_N_NEAREST];

  *displ = rhea_topography_displacement_node (label, x, y, z, opt);
}

static void
rhea_topography_displacement_node_face_fn (double *displ, double x, double y,
                                           double z, double nx, double ny,
                                           double nz, ymir_topidx_t face,
                                           ymir_locidx_t node_id, void *data)
{
  rhea_topography_options_t *opt = data;
  int                 label[RHEA_TOPOGRAPHY_DISPLACEMENT_N_NEAREST];

  *displ = rhea_topography_displacement_node (label, x, y, z, opt);
}

void
rhea_topography_displacement_vec (ymir_vec_t *displacement,
                                  rhea_topography_options_t *opt)
{
  int                 is_cvec = 0;
  int                 is_dvec = 0;

  /* check input */
  RHEA_ASSERT (opt->type == RHEA_TOPOGRAPHY_NONE || displacement != NULL);

  /* set vector type */
  if (displacement != NULL) {
    is_cvec = ymir_vec_is_cvec (displacement);
    is_dvec = ymir_vec_is_dvec (displacement);
  }

  /* set values of displacement vector */
  switch (opt->type) {
  case RHEA_TOPOGRAPHY_NONE:
    if (displacement != NULL) {
      ymir_vec_set_value (displacement, NAN);
    }
    break;

  case RHEA_TOPOGRAPHY_DATA_POINTS_DISPLS:
  case RHEA_TOPOGRAPHY_DATA_POINTS_DISPLS_LABELS:
    if (!ymir_vec_is_face_vec (displacement)) { /* if volume vector */
      if (is_cvec) {
        ymir_cvec_set_function (
            displacement, rhea_topography_displacement_node_vol_fn, opt);
      }
      else if (is_dvec) {
        ymir_dvec_set_function (
            displacement, rhea_topography_displacement_node_vol_fn, opt);
      }
      else {
        RHEA_ABORT_NOT_REACHED ();
      }
    }
    else { /* if face vector */
      if (is_cvec) {
        ymir_face_cvec_set_function (
            displacement, rhea_topography_displacement_node_face_fn, opt);
      }
      else if (is_dvec) {
        ymir_face_dvec_set_function (
            displacement, rhea_topography_displacement_node_face_fn, opt);
      }
      else {
        RHEA_ABORT_NOT_REACHED ();
      }
    }
    RHEA_ASSERT (!is_cvec || (sc_dmatrix_is_valid (displacement->dataown) &&
                              sc_dmatrix_is_valid (displacement->coff)));
    RHEA_ASSERT (!is_dvec || sc_dmatrix_is_valid (displacement->dataown));
    break;

  default: /* unknown topography type */
    RHEA_ABORT_NOT_REACHED ();
  }
}
