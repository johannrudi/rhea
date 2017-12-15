/*
 */

#include <rhea_weakzone.h>
#include <rhea_base.h>
#include <rhea_io_mpi.h>
#include <rhea_io_std.h>
#include <ymir_vec_getset.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_WEAKZONE_DEFAULT_TYPE_NAME "NONE"
#define RHEA_WEAKZONE_DEFAULT_THICKNESS_M (20.0e3)
#define RHEA_WEAKZONE_DEFAULT_THICKNESS_CONST_M (5.0e3)
#define RHEA_WEAKZONE_DEFAULT_WEAK_FACTOR_INTERIOR (1.0e-5)
#define RHEA_WEAKZONE_DEFAULT_POINTS_FILE_PATH_BIN NULL
#define RHEA_WEAKZONE_DEFAULT_POINTS_FILE_PATH_TXT NULL
#define RHEA_WEAKZONE_DEFAULT_LABELS_FILE_PATH_BIN NULL
#define RHEA_WEAKZONE_DEFAULT_LABELS_FILE_PATH_TXT NULL
#define RHEA_WEAKZONE_DEFAULT_FACTORS_FILE_PATH_BIN NULL
#define RHEA_WEAKZONE_DEFAULT_FACTORS_FILE_PATH_TXT NULL
#define RHEA_WEAKZONE_DEFAULT_N_POINTS (5000000)
#define RHEA_WEAKZONE_DEFAULT_WRITE_POINTS_FILE_PATH_BIN NULL
#define RHEA_WEAKZONE_DEFAULT_WRITE_POINTS_FILE_PATH_TXT NULL
#define RHEA_WEAKZONE_DEFAULT_WRITE_LABELS_FILE_PATH_BIN NULL
#define RHEA_WEAKZONE_DEFAULT_WRITE_FACTORS_FILE_PATH_BIN NULL

/* initialize options */
char               *rhea_weakzone_type_name = RHEA_WEAKZONE_DEFAULT_TYPE_NAME;
double              rhea_weakzone_thickness_m =
  RHEA_WEAKZONE_DEFAULT_THICKNESS_M;
double              rhea_weakzone_thickness_const_m =
  RHEA_WEAKZONE_DEFAULT_THICKNESS_CONST_M;
double              rhea_weakzone_weak_factor_interior =
  RHEA_WEAKZONE_DEFAULT_WEAK_FACTOR_INTERIOR;
char               *rhea_weakzone_points_file_path_bin =
  RHEA_WEAKZONE_DEFAULT_POINTS_FILE_PATH_BIN;
char               *rhea_weakzone_points_file_path_txt =
  RHEA_WEAKZONE_DEFAULT_POINTS_FILE_PATH_TXT;
char               *rhea_weakzone_labels_file_path_bin =
  RHEA_WEAKZONE_DEFAULT_LABELS_FILE_PATH_BIN;
char               *rhea_weakzone_labels_file_path_txt =
  RHEA_WEAKZONE_DEFAULT_LABELS_FILE_PATH_TXT;
char               *rhea_weakzone_factors_file_path_bin =
  RHEA_WEAKZONE_DEFAULT_FACTORS_FILE_PATH_BIN;
char               *rhea_weakzone_factors_file_path_txt =
  RHEA_WEAKZONE_DEFAULT_FACTORS_FILE_PATH_TXT;
int                 rhea_weakzone_n_points =
  RHEA_WEAKZONE_DEFAULT_N_POINTS;
char               *rhea_weakzone_write_points_file_path_bin =
  RHEA_WEAKZONE_DEFAULT_WRITE_POINTS_FILE_PATH_BIN;
char               *rhea_weakzone_write_points_file_path_txt =
  RHEA_WEAKZONE_DEFAULT_WRITE_POINTS_FILE_PATH_TXT;
char               *rhea_weakzone_write_labels_file_path_bin =
  RHEA_WEAKZONE_DEFAULT_WRITE_LABELS_FILE_PATH_BIN;
char               *rhea_weakzone_write_factors_file_path_bin =
  RHEA_WEAKZONE_DEFAULT_WRITE_FACTORS_FILE_PATH_BIN;

void
rhea_weakzone_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Weakzone";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_S, "type", '\0',
    &(rhea_weakzone_type_name), RHEA_WEAKZONE_DEFAULT_TYPE_NAME,
    "Weak zone type name: data_points, data_points_labels, "
    "data_points_labels_factors",

  YMIR_OPTIONS_D, "thickness", '\0',
    &(rhea_weakzone_thickness_m), RHEA_WEAKZONE_DEFAULT_THICKNESS_M,
    "Width about center of weak zone [m]",
  YMIR_OPTIONS_D, "thickness-const", '\0',
    &(rhea_weakzone_thickness_const_m), RHEA_WEAKZONE_DEFAULT_THICKNESS_CONST_M,
    "Width of smoothing of edges of weak zone [m]",
  YMIR_OPTIONS_D, "weak-factor-interior", '\0',
    &(rhea_weakzone_weak_factor_interior),
    RHEA_WEAKZONE_DEFAULT_WEAK_FACTOR_INTERIOR,
    "Min weak zone factor, which is assumed in the weak zone's interior",

  YMIR_OPTIONS_S, "points-file-path-bin", '\0',
    &(rhea_weakzone_points_file_path_bin),
    RHEA_WEAKZONE_DEFAULT_POINTS_FILE_PATH_BIN,
    "Path to a binary file with (x,y,z) coordinates of weak zone points",
  YMIR_OPTIONS_S, "points-file-path-txt", '\0',
    &(rhea_weakzone_points_file_path_txt),
    RHEA_WEAKZONE_DEFAULT_POINTS_FILE_PATH_TXT,
    "Path to a text file with (x,y,z) coordinates of weak zone points",
  YMIR_OPTIONS_S, "labels-file-path-bin", '\0',
    &(rhea_weakzone_labels_file_path_bin),
    RHEA_WEAKZONE_DEFAULT_LABELS_FILE_PATH_BIN,
    "Path to a binary file with labels of weak zone points",
  YMIR_OPTIONS_S, "labels-file-path-txt", '\0',
    &(rhea_weakzone_labels_file_path_txt),
    RHEA_WEAKZONE_DEFAULT_LABELS_FILE_PATH_TXT,
    "Path to a text file with labels of weak zone points",
  YMIR_OPTIONS_S, "factors-file-path-bin", '\0',
    &(rhea_weakzone_factors_file_path_bin),
    RHEA_WEAKZONE_DEFAULT_FACTORS_FILE_PATH_BIN,
    "Path to a binary file with factors of weak zone points",
  YMIR_OPTIONS_S, "factors-file-path-txt", '\0',
    &(rhea_weakzone_factors_file_path_txt),
    RHEA_WEAKZONE_DEFAULT_FACTORS_FILE_PATH_TXT,
    "Path to a text file with factors of weak zone points",

  YMIR_OPTIONS_I, "num-points", '\0',
    &(rhea_weakzone_n_points), RHEA_WEAKZONE_DEFAULT_N_POINTS,
    "Number of points that are imported",

  YMIR_OPTIONS_S, "write-points-file-path-bin", '\0',
    &(rhea_weakzone_write_points_file_path_bin),
    RHEA_WEAKZONE_DEFAULT_WRITE_POINTS_FILE_PATH_BIN,
    "Output path for a binary file with coordinates of weak zone points",
  YMIR_OPTIONS_S, "write-points-file-path-txt", '\0',
    &(rhea_weakzone_write_points_file_path_txt),
    RHEA_WEAKZONE_DEFAULT_WRITE_POINTS_FILE_PATH_TXT,
    "Output path for a text file with coordinates of weak zone points",
  YMIR_OPTIONS_S, "write-labels-file-path-bin", '\0',
    &(rhea_weakzone_write_labels_file_path_bin),
    RHEA_WEAKZONE_DEFAULT_WRITE_LABELS_FILE_PATH_BIN,
    "Output path for a binary file with labels of weak zone points",
  YMIR_OPTIONS_S, "write-factors-file-path-bin", '\0',
    &(rhea_weakzone_write_factors_file_path_bin),
    RHEA_WEAKZONE_DEFAULT_WRITE_FACTORS_FILE_PATH_BIN,
    "Output path for a binary file with factors of weak zone points",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);
}

void
rhea_weakzone_process_options (rhea_weakzone_options_t *opt,
                               rhea_domain_options_t *domain_options)
{
  /* set weak zone type */
  if (strcmp (rhea_weakzone_type_name, "NONE") == 0) {
    opt->type = RHEA_WEAKZONE_NONE;
  }
  else if (strcmp (rhea_weakzone_type_name, "data_points") == 0) {
    opt->type = RHEA_WEAKZONE_DATA_POINTS;
  }
  else if (strcmp (rhea_weakzone_type_name, "data_points_labels") == 0) {
    opt->type = RHEA_WEAKZONE_DATA_POINTS_LABELS;
  }
  else if (strcmp (rhea_weakzone_type_name, "data_points_labels_factors") == 0)
  {
    opt->type = RHEA_WEAKZONE_DATA_POINTS_LABELS_FACTORS;
  }
  else { /* unknown type name */
    RHEA_ABORT ("Unknown weak zone type name");
  }

  /* set (nondimensional) parameters for smoothed weak zones */
  opt->thickness = rhea_weakzone_thickness_m / domain_options->radius_max_m;
  opt->thickness_const = rhea_weakzone_thickness_const_m /
                         domain_options->radius_max_m;
  opt->weak_factor_interior = rhea_weakzone_weak_factor_interior;

  /* set paths to binary & text files */
  opt->points_file_path_bin = rhea_weakzone_points_file_path_bin;
  opt->points_file_path_txt = rhea_weakzone_points_file_path_txt;
  opt->labels_file_path_bin = rhea_weakzone_labels_file_path_bin;
  opt->labels_file_path_txt = rhea_weakzone_labels_file_path_txt;
  opt->factors_file_path_bin = rhea_weakzone_factors_file_path_bin;
  opt->factors_file_path_txt = rhea_weakzone_factors_file_path_txt;

  /* set number of points */
  opt->n_points = rhea_weakzone_n_points;

  /* set output paths */
  opt->write_points_file_path_bin = rhea_weakzone_write_points_file_path_bin;
  opt->write_points_file_path_txt = rhea_weakzone_write_points_file_path_txt;
  opt->write_labels_file_path_bin = rhea_weakzone_write_labels_file_path_bin;
  opt->write_factors_file_path_bin = rhea_weakzone_write_factors_file_path_bin;

  /* init data */
  opt->pointcloud = NULL;

  /* set dependent options */
  opt->domain_options = domain_options;
}

int
rhea_weakzone_exists (rhea_weakzone_options_t *opt)
{
  switch (opt->type) {
  case RHEA_WEAKZONE_NONE:
    return 0;
  case RHEA_WEAKZONE_DATA_POINTS:
  case RHEA_WEAKZONE_DATA_POINTS_LABELS:
  case RHEA_WEAKZONE_DATA_POINTS_LABELS_FACTORS:
    return 1;
  default: /* unknown weak zone type */
    RHEA_ABORT_NOT_REACHED ();
    return 0;
  }
}

/******************************************************************************
 * Weak Zone Vector
 *****************************************************************************/

ymir_vec_t *
rhea_weakzone_new (ymir_mesh_t *ymir_mesh)
{
  return ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
}

void
rhea_weakzone_destroy (ymir_vec_t *weakzone)
{
  ymir_vec_destroy (weakzone);
}

int
rhea_weakzone_check_vec_type (ymir_vec_t *vec)
{
  return (
      ymir_vec_is_dvec (vec) &&
      vec->ndfields == 1 &&
      vec->node_type == YMIR_GAUSS_NODE
  );
}

int
rhea_weakzone_is_valid (ymir_vec_t *vec)
{
  return (
      sc_dmatrix_is_valid (vec->dataown) &&
      0.0 < ymir_dvec_min_global (vec) && ymir_dvec_max_global (vec) <= 1.0
  );
}

/******************************************************************************
 * Data
 *****************************************************************************/

void
rhea_weakzone_data_create (rhea_weakzone_options_t *opt, sc_MPI_Comm mpicomm)
{
  const char         *this_fn_name = "rhea_weakzone_data_create";
  rhea_domain_options_t *domain_options = opt->domain_options;
  int                 n_read;
  int                 create_coordinates;
  int                 create_labels;
  int                 create_factors;
  double             *coordinates = NULL;
  int                *labels = NULL;
  double             *factors = NULL;

  /* init flags */
  if (opt->points_file_path_bin != NULL ||
      opt->points_file_path_txt != NULL) {
    create_coordinates = 1;
  }
  else {
    create_coordinates = 0;
  }
  if (opt->labels_file_path_bin != NULL ||
      opt->labels_file_path_txt != NULL) {
    create_labels = 1;
  }
  else {
    create_labels = 0;
  }
  if (opt->factors_file_path_bin != NULL ||
      opt->factors_file_path_txt != NULL) {
    create_factors = 1;
  }
  else {
    create_factors = 0;
  }

  /* decide what data needs to be created */
  switch (opt->type) {
  case RHEA_WEAKZONE_NONE:
    return;
  case RHEA_WEAKZONE_DATA_POINTS:
    RHEA_CHECK_ABORT (create_coordinates,
                      "file path missing: weak zone coordinates");
    create_labels = 0;
    create_factors = 0;
    break;
  case RHEA_WEAKZONE_DATA_POINTS_LABELS:
    RHEA_CHECK_ABORT (create_coordinates,
                      "file path missing: weak zone coordinates");
    RHEA_CHECK_ABORT (create_labels,
                      "file path missing: weak zone labels");
    create_factors = 0;
    break;
  case RHEA_WEAKZONE_DATA_POINTS_LABELS_FACTORS:
    RHEA_CHECK_ABORT (create_coordinates,
                      "file path missing: weak zone coordinates");
    RHEA_CHECK_ABORT (create_labels,
                      "file path missing: weak zone labels");
    RHEA_CHECK_ABORT (create_factors,
                      "file path missing: weak zone factors");
    break;
  default: /* unknown weak zone type */
    RHEA_ABORT_NOT_REACHED ();
  }

  RHEA_GLOBAL_INFOF ("Into %s (type %i)\n", this_fn_name, opt->type);

  /* check input */
  RHEA_ASSERT (0 < opt->n_points);

  /* create, read, and write coordinates (updates #points in options) */
  if (create_coordinates) {
    const int           n_entries = 3 * opt->n_points;
    const char         *file_path_bin = opt->points_file_path_bin;
    const char         *file_path_txt = opt->points_file_path_txt;
    const char         *write_file_path_bin = opt->write_points_file_path_bin;
    const char         *write_file_path_txt = opt->write_points_file_path_txt;

    coordinates = RHEA_ALLOC (double, n_entries);

    /* read from file */
    if (file_path_bin != NULL) { /* if read from binary file */
      n_read = rhea_io_mpi_read_broadcast_double (
          coordinates, n_entries, file_path_bin, NULL /* path txt */, mpicomm);
      RHEA_ASSERT (n_read == n_entries);
    }
    else if (file_path_txt != NULL) { /* if read from text file */
      n_read = rhea_io_mpi_read_broadcast_double (
          coordinates, 0 /* #entries */, NULL /* path bin */, file_path_txt,
          mpicomm);
      RHEA_ASSERT (0 < n_read && n_read <= n_entries);
      RHEA_ASSERT (!(n_read % 3));
      opt->n_points = n_read / 3;
    }
    else { /* otherwise no reading possible */
      RHEA_ABORT_NOT_REACHED ();
    }

    RHEA_GLOBAL_INFOF ("%s: Number of weak zone points: %i\n", this_fn_name,
                       opt->n_points);

    /* write to file */
    if (write_file_path_bin != NULL || write_file_path_txt != NULL) {
      int                 mpirank, mpiret;

      /* get parallel environment */
      mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);

      /* write coordinates to binary file */
      if (mpirank == 0 && write_file_path_bin != NULL) {
        rhea_io_std_write_double (write_file_path_bin, coordinates,
                                  (size_t) 3 * opt->n_points);
      }

      /* write coordinates to text file */
      if (mpirank == 0 && write_file_path_txt != NULL) {
        rhea_io_std_write_double_to_txt (write_file_path_txt, coordinates,
                                         (size_t) 3 * opt->n_points,
                                         3 /* entries per line */);
      }
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

  /* create, read, and write factors */
  if (create_factors) {
    const int           n_entries = opt->n_points;
    const char         *file_path_bin = opt->factors_file_path_bin;
    const char         *file_path_txt = opt->factors_file_path_txt;
    const char         *write_file_path_bin = opt->write_factors_file_path_bin;

    factors = RHEA_ALLOC (double, n_entries);

    /* read from file */
    if (file_path_bin != NULL) { /* if read from binary file */
      n_read = rhea_io_mpi_read_broadcast_double (
          factors, n_entries, file_path_bin, NULL /* path txt */, mpicomm);
      RHEA_ASSERT (n_read == n_entries);
    }
    else if (file_path_txt != NULL) { /* if read from text file */
      n_read = rhea_io_mpi_read_broadcast_double (
          factors, 0 /* #entries */, write_file_path_bin, file_path_txt,
          mpicomm);
      RHEA_ASSERT (0 < n_read && n_read <= n_entries);
    }
    else { /* otherwise no reading possible */
      RHEA_ABORT_NOT_REACHED ();
    }
  }

  /* create point cloud */
  opt->pointcloud = rhea_pointcloud_weakzone_new (
      domain_options->x_min, domain_options->x_max,
      domain_options->y_min, domain_options->y_max,
      domain_options->z_min, domain_options->z_max,
      coordinates, opt->n_points);
  if (create_labels) {
    rhea_pointcloud_weakzone_set_labels (opt->pointcloud, labels);
  }
  if (create_factors) {
    rhea_pointcloud_weakzone_set_factors (opt->pointcloud, factors);
  }

  /* destroy */
  RHEA_FREE (coordinates);
  if (create_labels) {
    RHEA_FREE (labels);
  }
  if (create_factors) {
    RHEA_FREE (factors);
  }

  RHEA_GLOBAL_INFOF ("Done %s (type %i)\n", this_fn_name, opt->type);
}

void
rhea_weakzone_data_clear (rhea_weakzone_options_t *opt)
{
  if (opt->pointcloud != NULL) {
    rhea_pointcloud_weakzone_destroy (opt->pointcloud);
  }
}

/******************************************************************************
 * Weak Zone Computation
 *****************************************************************************/

/**
 * Computes the distance to the weak zone surface (e.g., shortest distance to
 * point cloud).
 */
double
rhea_weakzone_dist_node (int *nearest_label, double *nearest_factor,
                         const double x, const double y, const double z,
                         rhea_weakzone_options_t *opt)
{
  const double        pt[3] = {x, y, z};
  double              dist;

  switch (opt->type) {
  case RHEA_WEAKZONE_DATA_POINTS:
    RHEA_ASSERT (opt->pointcloud != NULL);
    dist = rhea_pointcloud_weakzone_find_nearest (NULL, NULL, NULL,
                                                  opt->pointcloud, pt);
    break;
  case RHEA_WEAKZONE_DATA_POINTS_LABELS:
    RHEA_ASSERT (opt->pointcloud != NULL);
    RHEA_ASSERT (nearest_label != NULL);
    dist = rhea_pointcloud_weakzone_find_nearest (NULL, NULL, nearest_label,
                                                  opt->pointcloud, pt);
    break;
  case RHEA_WEAKZONE_DATA_POINTS_LABELS_FACTORS:
    RHEA_ASSERT (opt->pointcloud != NULL);
    RHEA_ASSERT (nearest_label != NULL);
    RHEA_ASSERT (nearest_factor != NULL);
    dist = rhea_pointcloud_weakzone_find_nearest (NULL, nearest_factor,
                                                  nearest_label,
                                                  opt->pointcloud, pt);
    RHEA_ASSERT (isfinite (*nearest_factor));
    RHEA_ASSERT (0.0 < *nearest_factor && *nearest_factor <= 1.0);
    break;
  default: /* unknown weak zone type */
    RHEA_ABORT_NOT_REACHED ();
  }
  RHEA_ASSERT (isfinite (dist));
  RHEA_ASSERT (0.0 <= dist);

  return dist;
}

/**
 * Computes the weak zone factor depending on the distance to a weak zone.
 * The edges of the weak zone are smoothed by a Gaussian.
 *
 *   1 - (1 - factor_interior) * exp ( - dist^2 / (2 * (0.5*thickness)^2) )
 */
double
rhea_weakzone_factor_node (const double distance,
                           const double thickness,
                           const double thickness_const,
                           const double factor_interior)
{
  const double        d = distance - 0.5 * thickness_const;
  const double        std_dev = 0.5 * (thickness - thickness_const);
  double              factor;

  /* check input */
  RHEA_ASSERT (isfinite (distance));
  RHEA_ASSERT (0.0 <= distance);
  RHEA_ASSERT (thickness_const <= thickness);
  RHEA_ASSERT (isfinite (factor_interior));
  RHEA_ASSERT (0.0 < factor_interior && factor_interior <= 1.0);

  if (d <= 0.0) { /* if inside (=constant) zone */
    factor = factor_interior;
  }
  else { /* otherwise in smoothed zone */
    factor = 1.0 - (1.0 - factor_interior) * exp (-d*d / (2.0*std_dev*std_dev));
  }
  RHEA_ASSERT (isfinite (factor));
  RHEA_ASSERT (0.0 < factor && factor <= 1.0);

  /* return weak factor */
  return factor;
}

double
rhea_weakzone_factor_deriv_node (const double distance,
                                 const double thickness,
                                 const double thickness_const,
                                 const double factor_interior)
{
  const double        d = distance - 0.5 * thickness_const;
  const double        std_dev = 0.5 * (thickness - thickness_const);
  double              factor_deriv;

  /* check input */
  RHEA_ASSERT (isfinite (distance));
  RHEA_ASSERT (0.0 <= distance);
  RHEA_ASSERT (thickness_const <= thickness);
  RHEA_ASSERT (isfinite (factor_interior));
  RHEA_ASSERT (0.0 < factor_interior && factor_interior <= 1.0);

  if (d <= 0.0) { /* if inside (=constant) zone */
    factor_deriv = 0.0;
  }
  else { /* otherwise in smoothed zone */
    factor_deriv = (1.0 - factor_interior) * (d / (std_dev*std_dev)) *
                   exp (-d*d / (2.0*std_dev*std_dev));
  }
  RHEA_ASSERT (isfinite (factor_deriv));
  RHEA_ASSERT (0.0 <= factor_deriv);

  /* return weak factor */
  return factor_deriv;
}

/**
 * Computes the value of the weak zone at given coordinates.
 */
static double
rhea_weakzone_node (const double x, const double y, const double z,
                    rhea_weakzone_options_t *opt)
{
  const double        thickness = opt->thickness;
  const double        thickness_const = opt->thickness_const;
  double              distance;
  int                 label; //TODO unused yet
  double              factor_interior;

  switch (opt->type) {
  case RHEA_WEAKZONE_DATA_POINTS:
  case RHEA_WEAKZONE_DATA_POINTS_LABELS:
    distance = rhea_weakzone_dist_node (&label, NULL, x, y, z, opt);
    factor_interior = opt->weak_factor_interior;
    break;
  case RHEA_WEAKZONE_DATA_POINTS_LABELS_FACTORS:
    distance = rhea_weakzone_dist_node (&label, &factor_interior, x, y, z, opt);
    break;
  default: /* unknown weak zone type */
    RHEA_ABORT_NOT_REACHED ();
  }

  return rhea_weakzone_factor_node (distance, thickness, thickness_const,
                                    factor_interior);
}

static void
rhea_weakzone_node_fn (double *w, double x, double y, double z,
                       ymir_locidx_t nid, void *data)
{
  rhea_weakzone_options_t *opt = data;

  *w = rhea_weakzone_node (x, y, z, opt);
}

void
rhea_weakzone_compute_elem (double *_sc_restrict weak_elem,
                            const double *_sc_restrict x,
                            const double *_sc_restrict y,
                            const double *_sc_restrict z,
                            const int n_nodes,
                            rhea_weakzone_options_t *opt)
{
  int                 nodeid;

  for (nodeid = 0; nodeid < n_nodes; nodeid++) {
    weak_elem[nodeid] = rhea_weakzone_node (x[nodeid], y[nodeid], z[nodeid],
                                            opt);
  }
}

void
rhea_weakzone_compute (ymir_vec_t *weakzone, void *data)
{
  rhea_weakzone_options_t *opt = data;

  /* check input */
  RHEA_ASSERT (opt->type == RHEA_WEAKZONE_NONE || weakzone != NULL);
  RHEA_ASSERT (opt->type == RHEA_WEAKZONE_NONE ||
               rhea_weakzone_check_vec_type (weakzone));

  switch (opt->type) {
  case RHEA_WEAKZONE_NONE:
    if (weakzone != NULL) {
      ymir_vec_set_value (weakzone, 1.0);
    }
    break;
  case RHEA_WEAKZONE_DATA_POINTS:
  case RHEA_WEAKZONE_DATA_POINTS_LABELS:
  case RHEA_WEAKZONE_DATA_POINTS_LABELS_FACTORS:
    ymir_dvec_set_function (weakzone, rhea_weakzone_node_fn, opt);
    break;
  default: /* unknown weak zone type */
    RHEA_ABORT_NOT_REACHED ();
  }
  RHEA_ASSERT (rhea_weakzone_is_valid (weakzone));
}

static void
rhea_weakzone_dist_node_fn (double *dist, double x, double y, double z,
                            ymir_locidx_t nid, void *data)
{
  rhea_weakzone_options_t *opt = data;
  int                 label;
  double              factor_interior;

  *dist = rhea_weakzone_dist_node (&label, &factor_interior, x, y, z, opt);
}

void
rhea_weakzone_compute_distance (ymir_vec_t *distance,
                                rhea_weakzone_options_t *opt)
{
  /* check input */
  RHEA_ASSERT (opt->type == RHEA_WEAKZONE_NONE || distance != NULL);
  RHEA_ASSERT (opt->type == RHEA_WEAKZONE_NONE ||
               rhea_weakzone_check_vec_type (distance));

  switch (opt->type) {
  case RHEA_WEAKZONE_NONE:
    if (distance != NULL) {
      ymir_vec_set_value (distance, NAN);
    }
    break;
  case RHEA_WEAKZONE_DATA_POINTS:
  case RHEA_WEAKZONE_DATA_POINTS_LABELS:
  case RHEA_WEAKZONE_DATA_POINTS_LABELS_FACTORS:
    ymir_dvec_set_function (distance, rhea_weakzone_dist_node_fn, opt);
    break;
  default: /* unknown weak zone type */
    RHEA_ABORT_NOT_REACHED ();
  }
  RHEA_ASSERT (rhea_weakzone_is_valid (distance));
}

/******************************************************************************
 * Get & Set Values
 *****************************************************************************/

double *
rhea_weakzone_get_elem_gauss (sc_dmatrix_t *weak_el_mat, ymir_vec_t *weak_vec,
                              const ymir_locidx_t elid)
{
#ifdef RHEA_ENABLE_DEBUG
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (weak_vec);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  /* check input */
  RHEA_ASSERT (rhea_weakzone_check_vec_type (weak_vec));
  RHEA_ASSERT (weak_el_mat->m == n_nodes_per_el);
  RHEA_ASSERT (weak_el_mat->n == 1);
#endif

  ymir_dvec_get_elem (weak_vec, weak_el_mat, YMIR_STRIDE_NODE, elid, YMIR_READ);
  return weak_el_mat->e[0];
}

void
rhea_weakzone_set_elem_gauss (ymir_vec_t *weak_vec, sc_dmatrix_t *weak_el_mat,
                              const ymir_locidx_t elid)
{
#ifdef RHEA_ENABLE_DEBUG
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (weak_vec);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  /* check input */
  RHEA_ASSERT (rhea_weakzone_check_vec_type (weak_vec));
  RHEA_ASSERT (weak_el_mat->m == n_nodes_per_el);
  RHEA_ASSERT (weak_el_mat->n == 1);
#endif

  ymir_dvec_set_elem (weak_vec, weak_el_mat, YMIR_STRIDE_NODE, elid, YMIR_SET);
}
