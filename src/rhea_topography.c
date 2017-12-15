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
  const char         *this_fn_name = "rhea_topography_data_create";
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
  default: /* unknown weak zone type */
    RHEA_ABORT_NOT_REACHED ();
  }

  RHEA_GLOBAL_INFOF ("Into %s\n", this_fn_name);

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

    RHEA_GLOBAL_INFOF ("%s: Number of topography points: %i\n", this_fn_name,
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

  RHEA_GLOBAL_INFOF ("Done %s\n", this_fn_name);
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

//TODO
