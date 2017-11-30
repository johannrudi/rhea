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
#define RHEA_WEAKZONE_DEFAULT_N_POINTS (5000000)
#define RHEA_WEAKZONE_DEFAULT_WRITE_POINTS_FILE_PATH_BIN NULL
#define RHEA_WEAKZONE_DEFAULT_WRITE_POINTS_FILE_PATH_TXT NULL

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
int                 rhea_weakzone_n_points =
  RHEA_WEAKZONE_DEFAULT_N_POINTS;
char               *rhea_weakzone_write_points_file_path_bin =
  RHEA_WEAKZONE_DEFAULT_WRITE_POINTS_FILE_PATH_BIN;
char               *rhea_weakzone_write_points_file_path_txt =
  RHEA_WEAKZONE_DEFAULT_WRITE_POINTS_FILE_PATH_TXT;

void
rhea_weakzone_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Weakzone";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_S, "type", '\0',
    &(rhea_weakzone_type_name), RHEA_WEAKZONE_DEFAULT_TYPE_NAME,
    "Weak zone type name: 'DATA'",

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

  YMIR_OPTIONS_I, "num-points", '\0',
    &(rhea_weakzone_n_points), RHEA_WEAKZONE_DEFAULT_N_POINTS,
    "Number of points that are imported",

  YMIR_OPTIONS_S, "write-points-file-path-bin", '\0',
    &(rhea_weakzone_write_points_file_path_bin),
    RHEA_WEAKZONE_DEFAULT_WRITE_POINTS_FILE_PATH_BIN,
    "Output path for a binary file with weak zone points",

  YMIR_OPTIONS_S, "write-points-file-path-txt", '\0',
    &(rhea_weakzone_write_points_file_path_txt),
    RHEA_WEAKZONE_DEFAULT_WRITE_POINTS_FILE_PATH_TXT,
    "Output path for a text file with weak zone points",

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
  else { /* unknown type name */
    RHEA_ABORT ("Unknown weak zone type name");
  }

  /* set (nondimensional) parameters for smoothed weak zones */
  opt->thickness = rhea_weakzone_thickness_m / domain_options->radius_max_m;
  opt->thickness_const = rhea_weakzone_thickness_const_m /
                         domain_options->radius_max_m;
  opt->weak_factor_interior = rhea_weakzone_weak_factor_interior;

  /* set paths to binary & text files of (x,y,z) coordinates */
  opt->points_file_path_bin = rhea_weakzone_points_file_path_bin;
  opt->points_file_path_txt = rhea_weakzone_points_file_path_txt;

  /* set number of points */
  opt->n_points = rhea_weakzone_n_points;

  /* set output paths */
  opt->write_points_file_path_bin = rhea_weakzone_write_points_file_path_bin;
  opt->write_points_file_path_txt = rhea_weakzone_write_points_file_path_txt;

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

/******************************************************************************
 * Data
 *****************************************************************************/

void
rhea_weakzone_data_create (rhea_weakzone_options_t *opt, sc_MPI_Comm mpicomm)
{
  const char         *this_fn_name = "rhea_weakzone_data_create";
  rhea_domain_options_t *domain_options = opt->domain_options;
  const char         *file_path_bin = opt->points_file_path_bin;
  const char         *file_path_txt = opt->points_file_path_txt;
  const char         *write_file_path_bin = opt->write_points_file_path_bin;
  const char         *write_file_path_txt = opt->write_points_file_path_txt;
  const int           n_entries = 3 * opt->n_points;
  int                 n_read;

  int                 create_coordinates = 0;
  double             *coordinates = NULL;

//TODO
//int                 create_labels = 0;
//int                *labels = NULL;

//TODO
//int                 create_factors = 0;
//double             *factors = NULL;

  /* decide what data needs to be created */
  switch (opt->type) {
  case RHEA_WEAKZONE_NONE:
    return;
  case RHEA_WEAKZONE_DATA_POINTS:
    create_coordinates = 1;
    break;
  default: /* unknown weak zone type */
    RHEA_ABORT_NOT_REACHED ();
  }

  RHEA_GLOBAL_INFOF ("Into %s (type %i)\n", this_fn_name, opt->type);

  /* check input */
  RHEA_ASSERT (0 < opt->n_points);

  /* create coordinates */
  if (create_coordinates) {
    coordinates = RHEA_ALLOC (double, n_entries);

    /* read coordinates of weak zone points */
    if (file_path_bin != NULL) { /* if read from binary file */
      n_read = rhea_io_mpi_read_broadcast_double (
          coordinates, n_entries, file_path_bin, NULL /* path txt */, mpicomm);
      RHEA_ASSERT (n_entries == n_read);
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

    /* write weak zone data */
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

  /* create point cloud */
  opt->pointcloud = rhea_pointcloud_weakzone_new (
      domain_options->x_min, domain_options->x_max,
      domain_options->y_min, domain_options->y_max,
      domain_options->z_min, domain_options->z_max,
      coordinates, opt->n_points);

  /* destroy */
  RHEA_FREE (coordinates);

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
static double
rhea_weakzone_dist_node (const double x, const double y, const double z,
                         rhea_weakzone_options_t *opt)
{
  const double        pt[3] = {x, y, z};
  double              dist;

  switch (opt->type) {
  case RHEA_WEAKZONE_DATA_POINTS:
  case RHEA_WEAKZONE_DATA_POINTS_LABELS:
  case RHEA_WEAKZONE_DATA_POINTS_LABELS_FACTORS:
    RHEA_ASSERT (opt->pointcloud != NULL);
    dist = rhea_pointcloud_weakzone_find_nearest (NULL, NULL, NULL,
                                                  opt->pointcloud, pt);
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
 *   1 - (1 - weak_factor_interior) * exp ( - dist^2 / (2 * (0.5*thickness)^2) )
 */
static double
rhea_weakzone_factor_node (const double distance,
                           const double thickness,
                           const double thickness_const,
                           const double weak_factor_interior)
{
  const double        d = distance - 0.5 * thickness_const;
  const double        std_dev = 0.5 * (thickness - thickness_const);

  /* check input */
  RHEA_ASSERT (isfinite (distance));
  RHEA_ASSERT (0.0 <= distance);
  RHEA_ASSERT (thickness_const <= thickness);

  if (d <= 0.0) {
    /* return value inside zone with constant weak factor */
    return weak_factor_interior;
  }
  else {
    /* return smoothed weak zone */
    return 1.0 - (1.0 - weak_factor_interior) *
                 exp (-d*d / (2.0 * std_dev*std_dev));
  }
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
  const double        weak_factor_interior = opt->weak_factor_interior;
  double              distance;

  distance = rhea_weakzone_dist_node (x, y, z, opt);
  return rhea_weakzone_factor_node (distance, thickness, thickness_const,
                                    weak_factor_interior);
}

static void
rhea_weakzone_node_fn (double *w, double x, double y, double z,
                       ymir_locidx_t nid, void *data)
{
  rhea_weakzone_options_t *opt = data;

  *w = rhea_weakzone_node (x, y, z, opt);
}

void
rhea_weakzone_compute (ymir_vec_t *weakzone, rhea_weakzone_options_t *opt)
{
  /* check input */
  RHEA_ASSERT (opt->type == RHEA_WEAKZONE_NONE || weakzone != NULL);

  switch (opt->type) {
  case RHEA_WEAKZONE_NONE:
    if (weakzone != NULL) {
      ymir_vec_set_value (weakzone, 1.0);
    }
    break;
  case RHEA_WEAKZONE_DATA_POINTS:
    ymir_dvec_set_function (weakzone, rhea_weakzone_node_fn, opt);
    break;
  default: /* unknown weak zone type */
    RHEA_ABORT_NOT_REACHED ();
  }
}

static void
rhea_weakzone_dist_node_fn (double *dist, double x, double y, double z,
                            ymir_locidx_t nid, void *data)
{
  rhea_weakzone_options_t *opt = data;

  *dist = rhea_weakzone_dist_node (x, y, z, opt);
}

void
rhea_weakzone_compute_distance (ymir_vec_t *distance,
                                rhea_weakzone_options_t *opt)
{
  /* check input */
  RHEA_ASSERT (opt->type == RHEA_WEAKZONE_NONE || distance != NULL);

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
}
