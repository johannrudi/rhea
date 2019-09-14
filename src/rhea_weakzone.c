#include <rhea_weakzone.h>
#include <rhea_weakzone_label.h>
#include <rhea_base.h>
#include <rhea_io_mpi.h>
#include <rhea_io_std.h>
#include <ymir_vec_getset.h>
#include <ymir_perf_counter.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_WEAKZONE_DEFAULT_TYPE_NAME "NONE"
#define RHEA_WEAKZONE_DEFAULT_POINTS_FILE_PATH_BIN NULL
#define RHEA_WEAKZONE_DEFAULT_POINTS_FILE_PATH_TXT NULL
#define RHEA_WEAKZONE_DEFAULT_LABELS_FILE_PATH_BIN NULL
#define RHEA_WEAKZONE_DEFAULT_LABELS_FILE_PATH_TXT NULL
#define RHEA_WEAKZONE_DEFAULT_FACTORS_FILE_PATH_BIN NULL
#define RHEA_WEAKZONE_DEFAULT_FACTORS_FILE_PATH_TXT NULL
#define RHEA_WEAKZONE_DEFAULT_N_POINTS (100000000)
#define RHEA_WEAKZONE_DEFAULT_WRITE_POINTS_FILE_PATH_BIN NULL
#define RHEA_WEAKZONE_DEFAULT_WRITE_POINTS_FILE_PATH_TXT NULL
#define RHEA_WEAKZONE_DEFAULT_WRITE_LABELS_FILE_PATH_BIN NULL
#define RHEA_WEAKZONE_DEFAULT_WRITE_FACTORS_FILE_PATH_BIN NULL
#define RHEA_WEAKZONE_DEFAULT_THICKNESS_M (20.0e3)
#define RHEA_WEAKZONE_DEFAULT_THICKNESS_CLASS_SLAB_M (NAN)
#define RHEA_WEAKZONE_DEFAULT_THICKNESS_CLASS_RIDGE_M (NAN)
#define RHEA_WEAKZONE_DEFAULT_THICKNESS_CLASS_FRACTURE_M (NAN)
#define RHEA_WEAKZONE_DEFAULT_THICKNESS_CONST_M (5.0e3)
#define RHEA_WEAKZONE_DEFAULT_THICKNESS_CONST_CLASS_SLAB_M (NAN)
#define RHEA_WEAKZONE_DEFAULT_THICKNESS_CONST_CLASS_RIDGE_M (NAN)
#define RHEA_WEAKZONE_DEFAULT_THICKNESS_CONST_CLASS_FRACTURE_M (NAN)
#define RHEA_WEAKZONE_DEFAULT_WEAK_FACTOR_INTERIOR (1.0e-5)
#define RHEA_WEAKZONE_DEFAULT_WEAK_FACTOR_INTERIOR_CLASS_SLAB (NAN)
#define RHEA_WEAKZONE_DEFAULT_WEAK_FACTOR_INTERIOR_CLASS_RIDGE (NAN)
#define RHEA_WEAKZONE_DEFAULT_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE (NAN)
#define RHEA_WEAKZONE_DEFAULT_WEAK_FACTOR_INTERIOR_EARTH_FILE_PATH_TXT NULL
#define RHEA_WEAKZONE_DEFAULT_MONITOR_PERFORMANCE (0)

/* initialize options */
char               *rhea_weakzone_type_name = RHEA_WEAKZONE_DEFAULT_TYPE_NAME;
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
double              rhea_weakzone_thickness_m =
  RHEA_WEAKZONE_DEFAULT_THICKNESS_M;
double              rhea_weakzone_thickness_class_slab_m =
  RHEA_WEAKZONE_DEFAULT_THICKNESS_CLASS_SLAB_M;
double              rhea_weakzone_thickness_class_ridge_m =
  RHEA_WEAKZONE_DEFAULT_THICKNESS_CLASS_RIDGE_M;
double              rhea_weakzone_thickness_class_fracture_m =
  RHEA_WEAKZONE_DEFAULT_THICKNESS_CLASS_FRACTURE_M;
double              rhea_weakzone_thickness_const_m =
  RHEA_WEAKZONE_DEFAULT_THICKNESS_CONST_M;
double              rhea_weakzone_thickness_const_class_slab_m =
  RHEA_WEAKZONE_DEFAULT_THICKNESS_CONST_CLASS_SLAB_M;
double              rhea_weakzone_thickness_const_class_ridge_m =
  RHEA_WEAKZONE_DEFAULT_THICKNESS_CONST_CLASS_RIDGE_M;
double              rhea_weakzone_thickness_const_class_fracture_m =
  RHEA_WEAKZONE_DEFAULT_THICKNESS_CONST_CLASS_FRACTURE_M;
double              rhea_weakzone_weak_factor_interior =
  RHEA_WEAKZONE_DEFAULT_WEAK_FACTOR_INTERIOR;
double              rhea_weakzone_weak_factor_interior_class_slab =
  RHEA_WEAKZONE_DEFAULT_WEAK_FACTOR_INTERIOR_CLASS_SLAB;
double              rhea_weakzone_weak_factor_interior_class_ridge =
  RHEA_WEAKZONE_DEFAULT_WEAK_FACTOR_INTERIOR_CLASS_RIDGE;
double              rhea_weakzone_weak_factor_interior_class_fracture =
  RHEA_WEAKZONE_DEFAULT_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE;
char               *rhea_weakzone_weak_factor_interior_earth_file_path_txt =
  RHEA_WEAKZONE_DEFAULT_WEAK_FACTOR_INTERIOR_EARTH_FILE_PATH_TXT;
int                 rhea_weakzone_monitor_performance =
                      RHEA_WEAKZONE_DEFAULT_MONITOR_PERFORMANCE;

void
rhea_weakzone_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Weakzone";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_S, "type", '\0',
    &(rhea_weakzone_type_name), RHEA_WEAKZONE_DEFAULT_TYPE_NAME,
    "Weak zone type name: depth, data_points, data_points_labels, "
    "data_points_labels_factors",

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

  YMIR_OPTIONS_D, "thickness", '\0',
    &(rhea_weakzone_thickness_m), RHEA_WEAKZONE_DEFAULT_THICKNESS_M,
    "Width about center of weak zones [m]",
  YMIR_OPTIONS_D, "thickness-class-slab", '\0',
    &(rhea_weakzone_thickness_class_slab_m),
    RHEA_WEAKZONE_DEFAULT_THICKNESS_CLASS_SLAB_M,
    "Width about center of slabs [m]",
  YMIR_OPTIONS_D, "thickness-class-ridge", '\0',
    &(rhea_weakzone_thickness_class_ridge_m),
    RHEA_WEAKZONE_DEFAULT_THICKNESS_CLASS_RIDGE_M,
    "Width about center of ridges [m]",
  YMIR_OPTIONS_D, "thickness-class-fracture", '\0',
    &(rhea_weakzone_thickness_class_fracture_m),
    RHEA_WEAKZONE_DEFAULT_THICKNESS_CLASS_FRACTURE_M,
    "Width about center of fractures [m]",

  YMIR_OPTIONS_D, "thickness-const", '\0',
    &(rhea_weakzone_thickness_const_m), RHEA_WEAKZONE_DEFAULT_THICKNESS_CONST_M,
    "Width of interior with min factor in weak zones [m]",
  YMIR_OPTIONS_D, "thickness-const-class-slab", '\0',
    &(rhea_weakzone_thickness_const_class_slab_m),
    RHEA_WEAKZONE_DEFAULT_THICKNESS_CONST_CLASS_SLAB_M,
    "Width of interior with min factor in slabs [m]",
  YMIR_OPTIONS_D, "thickness-const-class-ridge", '\0',
    &(rhea_weakzone_thickness_const_class_ridge_m),
    RHEA_WEAKZONE_DEFAULT_THICKNESS_CONST_CLASS_RIDGE_M,
    "Width of interior with min factor in ridges [m]",
  YMIR_OPTIONS_D, "thickness-const-class-fracture", '\0',
    &(rhea_weakzone_thickness_const_class_fracture_m),
    RHEA_WEAKZONE_DEFAULT_THICKNESS_CONST_CLASS_FRACTURE_M,
    "Width of interior with min factor in fractures [m]",

  YMIR_OPTIONS_D, "weak-factor-interior", '\0',
    &(rhea_weakzone_weak_factor_interior),
    RHEA_WEAKZONE_DEFAULT_WEAK_FACTOR_INTERIOR,
    "Min weak zone factor, which is assumed in the weak zone's interior",
  YMIR_OPTIONS_D, "weak-factor-interior-class-slab", '\0',
    &(rhea_weakzone_weak_factor_interior_class_slab),
    RHEA_WEAKZONE_DEFAULT_WEAK_FACTOR_INTERIOR_CLASS_SLAB,
    "Min weak zone factor of slabs",
  YMIR_OPTIONS_D, "weak-factor-interior-class-ridge", '\0',
    &(rhea_weakzone_weak_factor_interior_class_ridge),
    RHEA_WEAKZONE_DEFAULT_WEAK_FACTOR_INTERIOR_CLASS_RIDGE,
    "Min weak zone factor of ridges",
  YMIR_OPTIONS_D, "weak-factor-interior-class-fracture", '\0',
    &(rhea_weakzone_weak_factor_interior_class_fracture),
    RHEA_WEAKZONE_DEFAULT_WEAK_FACTOR_INTERIOR_CLASS_FRACTURE,
    "Min weak zone factor of fractures",

  YMIR_OPTIONS_S, "weak-factor-interior-earth-file-path-txt", '\0',
    &(rhea_weakzone_weak_factor_interior_earth_file_path_txt),
    RHEA_WEAKZONE_DEFAULT_WEAK_FACTOR_INTERIOR_EARTH_FILE_PATH_TXT,
    "Path to a text file with min weak zone factors for earth "
    "(needs to have 120 lines)",

  YMIR_OPTIONS_B, "monitor-performance", '\0',
    &(rhea_weakzone_monitor_performance),
    RHEA_WEAKZONE_DEFAULT_MONITOR_PERFORMANCE,
    "Measure and print performance statistics (e.g., runtime or flops)",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);

  /* initialize (deactivated) performance counters */
  rhea_weakzone_perfmon_init (0, 0);
}

void
rhea_weakzone_process_options (rhea_weakzone_options_t *opt,
                               rhea_domain_options_t *domain_options)
{
  /* set weak zone type */
  if (strcmp (rhea_weakzone_type_name, "NONE") == 0) {
    opt->type = RHEA_WEAKZONE_NONE;
  }
  else if (strcmp (rhea_weakzone_type_name, "depth") == 0) {
    opt->type = RHEA_WEAKZONE_DEPTH;
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

  /* set (nondimensional) parameters for weak zone geometry */
  RHEA_ASSERT (0.0 < domain_options->radius_max_m);
  RHEA_ASSERT (isfinite (rhea_weakzone_thickness_m));
  opt->thickness = rhea_weakzone_thickness_m / domain_options->radius_max_m;
  opt->thickness_class_slab     = opt->thickness;
  opt->thickness_class_ridge    = opt->thickness;
  opt->thickness_class_fracture = opt->thickness;
  if (isfinite (rhea_weakzone_thickness_class_slab_m)) {
    opt->thickness_class_slab = rhea_weakzone_thickness_class_slab_m /
                                domain_options->radius_max_m;
  }
  if (isfinite (rhea_weakzone_thickness_class_ridge_m)) {
    opt->thickness_class_ridge = rhea_weakzone_thickness_class_ridge_m /
                                 domain_options->radius_max_m;
  }
  if (isfinite (rhea_weakzone_thickness_class_fracture_m)) {
    opt->thickness_class_fracture = rhea_weakzone_thickness_class_fracture_m /
                                    domain_options->radius_max_m;
  }
  RHEA_ASSERT (isfinite (rhea_weakzone_thickness_const_m));
  opt->thickness_const = rhea_weakzone_thickness_const_m /
                         domain_options->radius_max_m;
  opt->thickness_const_class_slab     = opt->thickness_const;
  opt->thickness_const_class_ridge    = opt->thickness_const;
  opt->thickness_const_class_fracture = opt->thickness_const;
  if (isfinite (rhea_weakzone_thickness_const_class_slab_m)) {
    opt->thickness_const_class_slab =
      rhea_weakzone_thickness_const_class_slab_m /
      domain_options->radius_max_m;
  }
  if (isfinite (rhea_weakzone_thickness_const_class_ridge_m)) {
    opt->thickness_const_class_ridge =
      rhea_weakzone_thickness_const_class_ridge_m /
      domain_options->radius_max_m;
  }
  if (isfinite (rhea_weakzone_thickness_const_class_fracture_m)) {
    opt->thickness_const_class_fracture =
      rhea_weakzone_thickness_const_class_fracture_m /
      domain_options->radius_max_m;
  }

  /* set weak zone factors */
  RHEA_ASSERT (isfinite (rhea_weakzone_weak_factor_interior));
  opt->weak_factor_interior = rhea_weakzone_weak_factor_interior;
  opt->weak_factor_interior_class_slab     = opt->weak_factor_interior;
  opt->weak_factor_interior_class_ridge    = opt->weak_factor_interior;
  opt->weak_factor_interior_class_fracture = opt->weak_factor_interior;
  if (isfinite (rhea_weakzone_weak_factor_interior_class_slab)) {
    opt->weak_factor_interior_class_slab =
      rhea_weakzone_weak_factor_interior_class_slab;
  }
  if (isfinite (rhea_weakzone_weak_factor_interior_class_ridge)) {
    opt->weak_factor_interior_class_ridge =
      rhea_weakzone_weak_factor_interior_class_ridge;
  }
  if (isfinite (rhea_weakzone_weak_factor_interior_class_fracture)) {
    opt->weak_factor_interior_class_fracture =
      rhea_weakzone_weak_factor_interior_class_fracture;
  }
  opt->weak_factor_interior_earth = NULL;

  /* init data */
  opt->pointcloud = NULL;

  /* init statistics */
  opt->stats_radius_min = NAN;

  opt->stats_thickness_max = opt->thickness + opt->thickness_const;
  opt->stats_thickness_max = SC_MAX (
      opt->stats_thickness_max,
      opt->thickness_class_slab + opt->thickness_const_class_slab);
  opt->stats_thickness_max = SC_MAX (
      opt->stats_thickness_max,
      opt->thickness_class_ridge + opt->thickness_const_class_ridge);
  opt->stats_thickness_max = SC_MAX (
      opt->stats_thickness_max,
      opt->thickness_class_fracture + opt->thickness_const_class_fracture);

  opt->stats_factor_interior_min = opt->weak_factor_interior;
  opt->stats_factor_interior_min = SC_MIN (
      opt->stats_factor_interior_min,
      opt->weak_factor_interior_class_slab);
  opt->stats_factor_interior_min = SC_MIN (
      opt->stats_factor_interior_min,
      opt->weak_factor_interior_class_ridge);
  opt->stats_factor_interior_min = SC_MIN (
      opt->stats_factor_interior_min,
      opt->weak_factor_interior_class_fracture);

  /* set dependent options */
  opt->domain_options = domain_options;
}

int
rhea_weakzone_exists (rhea_weakzone_options_t *opt)
{
  switch (opt->type) {
  case RHEA_WEAKZONE_NONE:
    return 0;
  case RHEA_WEAKZONE_DEPTH:
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
 * Monitoring
 *****************************************************************************/

/* perfomance monitor tags and names */
typedef enum
{
  RHEA_WEAKZONE_PERFMON_CREATE_COORDINATES,
  RHEA_WEAKZONE_PERFMON_CREATE_LABELS,
  RHEA_WEAKZONE_PERFMON_CREATE_FACTORS,
  RHEA_WEAKZONE_PERFMON_CREATE_POINTCLOUD,
  RHEA_WEAKZONE_PERFMON_DIST_NODE,
  RHEA_WEAKZONE_PERFMON_N
}
rhea_weakzone_perfmon_idx_t;

static const char  *rhea_weakzone_perfmon_name[RHEA_WEAKZONE_PERFMON_N] =
{
  "Create coordinates",
  "Create labels",
  "Create factors",
  "Create point cloud",
  "Compute distance"
};
ymir_perf_counter_t rhea_weakzone_perfmon[RHEA_WEAKZONE_PERFMON_N];

void
rhea_weakzone_perfmon_init (const int activate, const int skip_if_active)
{
  const int           active = activate && rhea_weakzone_monitor_performance;

  ymir_perf_counter_init_all_ext (rhea_weakzone_perfmon,
                                  rhea_weakzone_perfmon_name,
                                  RHEA_WEAKZONE_PERFMON_N,
                                  active, skip_if_active);
}

void
rhea_weakzone_perfmon_print (sc_MPI_Comm mpicomm,
                        const int print_wtime,
                        const int print_n_calls,
                        const int print_flops)
{
  const int           active = rhea_weakzone_monitor_performance;
  const int           print = (print_wtime || print_n_calls || print_flops);
  int                 n_stats = RHEA_WEAKZONE_PERFMON_N *
                                YMIR_PERF_COUNTER_N_STATS;
  sc_statinfo_t       stats[n_stats];
  char                stats_name[n_stats][YMIR_PERF_COUNTER_NAME_SIZE];

  /* exit if nothing to do */
  if (!active || !print) {
    return;
  }

  /* gather performance statistics */
  n_stats = ymir_perf_counter_gather_stats (
      rhea_weakzone_perfmon, RHEA_WEAKZONE_PERFMON_N, stats, stats_name,
      mpicomm, print_wtime, print_n_calls, print_flops);

  /* print performance statistics */
  ymir_perf_counter_print_stats (stats, n_stats, "Weak Zone");
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
  case RHEA_WEAKZONE_DEPTH:
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

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "type=%i", opt->type);

  /* check input */
  RHEA_ASSERT (0 < opt->n_points);

  /* create, read, and write coordinates (updates #points in options) */
  if (create_coordinates) {
    const int           n_entries = 3 * opt->n_points;
    const char         *file_path_bin = opt->points_file_path_bin;
    const char         *file_path_txt = opt->points_file_path_txt;
    const char         *write_file_path_bin = opt->write_points_file_path_bin;
    const char         *write_file_path_txt = opt->write_points_file_path_txt;
    int                 idx;

    /* start performance monitors */
    ymir_perf_counter_start (
        &rhea_weakzone_perfmon[RHEA_WEAKZONE_PERFMON_CREATE_COORDINATES]);

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

    RHEA_GLOBAL_INFOF ("%s Number of weak zone points=%i\n", __func__,
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

    /* find min radius */
    opt->stats_radius_min = opt->domain_options->radius_max;
    for (idx = 0; idx < 3*opt->n_points; idx += 3) {
      const double        r = rhea_domain_compute_radius (coordinates[idx    ],
                                                          coordinates[idx + 1],
                                                          coordinates[idx + 2],
                                                          opt->domain_options);

      opt->stats_radius_min = SC_MIN (r, opt->stats_radius_min);
    }

    /* stop performance monitors */
    ymir_perf_counter_stop_add (
        &rhea_weakzone_perfmon[RHEA_WEAKZONE_PERFMON_CREATE_COORDINATES]);
  }

  /* create, read, and write labels */
  if (create_labels) {
    const int           n_entries = opt->n_points;
    const char         *file_path_bin = opt->labels_file_path_bin;
    const char         *file_path_txt = opt->labels_file_path_txt;
    const char         *write_file_path_bin = opt->write_labels_file_path_bin;

    /* start performance monitors */
    ymir_perf_counter_start (
        &rhea_weakzone_perfmon[RHEA_WEAKZONE_PERFMON_CREATE_LABELS]);

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

    RHEA_GLOBAL_INFOF ("%s Number of labels=%i\n", __func__, n_read);

    /* stop performance monitors */
    ymir_perf_counter_stop_add (
        &rhea_weakzone_perfmon[RHEA_WEAKZONE_PERFMON_CREATE_LABELS]);
  }

  /* create, read, and write factors */
  if (create_factors) {
    const int           n_entries = opt->n_points;
    const char         *file_path_bin = opt->factors_file_path_bin;
    const char         *file_path_txt = opt->factors_file_path_txt;
    const char         *write_file_path_bin = opt->write_factors_file_path_bin;
    int                 idx;

    /* start performance monitors */
    ymir_perf_counter_start (
        &rhea_weakzone_perfmon[RHEA_WEAKZONE_PERFMON_CREATE_FACTORS]);

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

    RHEA_GLOBAL_INFOF ("%s Number of factors=%i\n", __func__, n_read);

    /* find min factor */
    for (idx = 0; idx < n_entries; idx++) {
      opt->stats_factor_interior_min = SC_MIN (opt->stats_factor_interior_min,
                                               factors[idx]);
    }

    /* stop performance monitors */
    ymir_perf_counter_stop_add (
        &rhea_weakzone_perfmon[RHEA_WEAKZONE_PERFMON_CREATE_FACTORS]);
  }

  /* create point cloud */
  ymir_perf_counter_start (
      &rhea_weakzone_perfmon[RHEA_WEAKZONE_PERFMON_CREATE_POINTCLOUD]);
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
  ymir_perf_counter_stop_add (
      &rhea_weakzone_perfmon[RHEA_WEAKZONE_PERFMON_CREATE_POINTCLOUD]);

  /* destroy */
  RHEA_FREE (coordinates);
  if (create_labels) {
    RHEA_FREE (labels);
  }
  if (create_factors) {
    RHEA_FREE (factors);
  }

  /* read weak zone factors for each label */
  if (create_labels &&
      rhea_weakzone_weak_factor_interior_earth_file_path_txt != NULL) {
    const char         *file_path_txt =
      rhea_weakzone_weak_factor_interior_earth_file_path_txt;
    const int           n_entries = RHEA_WEAKZONE_LABEL_EARTH_N;
    int                 n_read;
    int                 idx;

    /* start performance monitors */
    ymir_perf_counter_start (
        &rhea_weakzone_perfmon[RHEA_WEAKZONE_PERFMON_CREATE_LABELS]);

    RHEA_ASSERT (opt->weak_factor_interior_earth == NULL);
    opt->weak_factor_interior_earth = RHEA_ALLOC (double, n_entries);

    n_read = rhea_io_mpi_read_broadcast_double (
        opt->weak_factor_interior_earth, n_entries, NULL /* path bin */,
        file_path_txt, mpicomm);
    RHEA_ASSERT (n_read == n_entries);

    RHEA_GLOBAL_INFOF ("%s Number of distinct labels=%i\n", __func__,
                       n_read);

    /* find min factor */
    for (idx = 0; idx < n_entries; idx++) {
      opt->stats_factor_interior_min = SC_MIN (
          opt->stats_factor_interior_min,
          opt->weak_factor_interior_earth[idx]);
    }

    /* stop performance monitors */
    ymir_perf_counter_stop_add (
        &rhea_weakzone_perfmon[RHEA_WEAKZONE_PERFMON_CREATE_LABELS]);
  }

  /* print statistics */
  RHEA_GLOBAL_INFOF ("%s Weak zone statistics: Min radius=%g\n",
                     __func__, opt->stats_radius_min);
  RHEA_GLOBAL_INFOF ("%s Weak zone statistics: Max thickness=%g\n",
                     __func__, opt->stats_thickness_max);
  RHEA_GLOBAL_INFOF ("%s Weak zone statistics: Min factor interior=%g\n",
                     __func__, opt->stats_factor_interior_min);

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

void
rhea_weakzone_data_clear (rhea_weakzone_options_t *opt)
{
  if (opt->pointcloud != NULL) {
    rhea_pointcloud_weakzone_destroy (opt->pointcloud);
  }
  if (opt->weak_factor_interior_earth != NULL) {
    RHEA_FREE (opt->weak_factor_interior_earth);
  }
}

/******************************************************************************
 * Weak Zone Computation
 *****************************************************************************/

double
rhea_weakzone_lookup_thickness (const int label, rhea_weakzone_options_t *opt)
{
  const rhea_weakzone_label_t class_id =
    rhea_weakzone_label_get_class ((rhea_weakzone_label_t) label);

  /* check input */
  RHEA_ASSERT (rhea_weakzone_label_is_valid_int (label));

  /* return thickness depending on label class */
  switch (class_id) {
  case RHEA_WEAKZONE_LABEL_CLASS_NONE:
    return opt->thickness;
  case RHEA_WEAKZONE_LABEL_CLASS_SLAB:
    return opt->thickness_class_slab;
  case RHEA_WEAKZONE_LABEL_CLASS_RIDGE:
    return opt->thickness_class_ridge;
  case RHEA_WEAKZONE_LABEL_CLASS_FRACTURE:
    return opt->thickness_class_fracture;
  default:
    RHEA_ABORT_NOT_REACHED ();
    return NAN;
  }
}

double
rhea_weakzone_lookup_thickness_const (const int label,
                                      rhea_weakzone_options_t *opt)
{
  const rhea_weakzone_label_t class_id =
    rhea_weakzone_label_get_class ((rhea_weakzone_label_t) label);

  /* check input */
  RHEA_ASSERT (rhea_weakzone_label_is_valid_int (label));

  /* return thickness depending on label class */
  switch (class_id) {
  case RHEA_WEAKZONE_LABEL_CLASS_NONE:
    return opt->thickness_const;
  case RHEA_WEAKZONE_LABEL_CLASS_SLAB:
    return opt->thickness_const_class_slab;
  case RHEA_WEAKZONE_LABEL_CLASS_RIDGE:
    return opt->thickness_const_class_ridge;
  case RHEA_WEAKZONE_LABEL_CLASS_FRACTURE:
    return opt->thickness_const_class_fracture;
  default:
    RHEA_ABORT_NOT_REACHED ();
    return NAN;
  }
}

double
rhea_weakzone_lookup_factor_interior (const int label,
                                      rhea_weakzone_options_t *opt)
{
  switch (opt->type) {
  case RHEA_WEAKZONE_DEPTH: /* single factor for all weak zones */
  case RHEA_WEAKZONE_DATA_POINTS:
    return opt->weak_factor_interior;

  case RHEA_WEAKZONE_DATA_POINTS_LABELS: /* different label dependent factors */
  case RHEA_WEAKZONE_DATA_POINTS_LABELS_FACTORS:
    RHEA_ASSERT (rhea_weakzone_label_is_valid_int (label));
    if (rhea_weakzone_label_is_class ((rhea_weakzone_label_t) label) ||
        opt->weak_factor_interior_earth == NULL) { /* if based on class */
      const rhea_weakzone_label_t class_id =
        rhea_weakzone_label_get_class ((rhea_weakzone_label_t) label);

      switch (class_id) {
      case RHEA_WEAKZONE_LABEL_CLASS_NONE:
        return opt->weak_factor_interior;
      case RHEA_WEAKZONE_LABEL_CLASS_SLAB:
        return opt->weak_factor_interior_class_slab;
      case RHEA_WEAKZONE_LABEL_CLASS_RIDGE:
        return opt->weak_factor_interior_class_ridge;
      case RHEA_WEAKZONE_LABEL_CLASS_FRACTURE:
        return opt->weak_factor_interior_class_fracture;
      default: /* unknown class */
        RHEA_ABORT_NOT_REACHED ();
        return NAN;
      }
    }
    else if (opt->weak_factor_interior_earth != NULL) { /* if vals individual */
      const int           idx =
        rhea_weakzone_label_earth_get_idx ((rhea_weakzone_label_t) label);

      RHEA_ASSERT (0 <= idx && idx < RHEA_WEAKZONE_LABEL_EARTH_N);
      return opt->weak_factor_interior_earth[idx];
    }
    else { /* otherwise use generic value */
      return opt->weak_factor_interior;
    }

  default: /* unknown weak zone type */
    RHEA_ABORT_NOT_REACHED ();
    return NAN;
  }
}

/**
 * Checks whether coordinates can expected to be far from weak zones s.t.
 *
 *   (1 - factor_interior) * exp ( - dist^2 / (2 * (0.5*thickness)^2) ) < eps
 */
static int
rhea_weakzone_dist_node_is_large (const double x, const double y,
                                  const double z, rhea_weakzone_options_t *opt)
{
  const double        eps = 1.0e-6;
  double              d, std_dev, factor_interior, thresh;

  /* exit if statistics not available */
  if (!isfinite (opt->stats_radius_min) ||
      !isfinite (opt->stats_thickness_max) ||
      !isfinite (opt->stats_factor_interior_min)) {
    return 0;
  }

  /* compute distance^2 */
  d = opt->stats_radius_min -
      rhea_domain_compute_radius (x, y, z, opt->domain_options);
  if (d <= 0.0) { /* if radius_min <= r(x,y,z) */
    return 0;
  }

  /* set square of max width */
  std_dev = 0.5 * opt->stats_thickness_max;

  /* set min factor */
  factor_interior = opt->stats_factor_interior_min;

  /* calculate threshold for distance^2 and return */
  thresh = - (2.0*std_dev*std_dev) * log (eps / (1.0 - factor_interior));
  RHEA_ASSERT (0.0 <= thresh);
  return (thresh < d*d);
}

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

  /* start performance monitors */
  ymir_perf_counter_start (
      &rhea_weakzone_perfmon[RHEA_WEAKZONE_PERFMON_DIST_NODE]);

  /* exit if distance is expected to be large */
  if (rhea_weakzone_dist_node_is_large (x, y, z, opt)) {
    if (nearest_label != NULL) {
      *nearest_label= RHEA_WEAKZONE_LABEL_CLASS_NONE;
    }
    if (nearest_factor != NULL) {
      *nearest_factor = RHEA_WEAKZONE_NEUTRAL_VALUE;
    }
    return opt->domain_options->depth;
  }

  /* compute distance */
  switch (opt->type) {
  case RHEA_WEAKZONE_DEPTH:
    dist = opt->domain_options->radius_max -
           rhea_domain_compute_radius (x, y, z, opt->domain_options);
    dist = SC_MAX (0.0, dist);
    break;
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

  /* stop performance monitors */
  ymir_perf_counter_stop_add (
      &rhea_weakzone_perfmon[RHEA_WEAKZONE_PERFMON_DIST_NODE]);

  return dist;
}

double
rhea_weakzone_indicator_node (const double distance,
                              const double thickness,
                              const double thickness_const)
{
  const double        d = distance - 0.5 * thickness_const;
  const double        std_dev = 0.5 * (thickness - thickness_const);
  double              indicator;

  /* check input */
  RHEA_ASSERT (isfinite (distance));
  RHEA_ASSERT (isfinite (thickness));
  RHEA_ASSERT (isfinite (thickness_const));
  RHEA_ASSERT (0.0 <= distance);
  RHEA_ASSERT (thickness_const <= thickness);

  if (d <= 0.0) { /* if inside (=constant) zone */
    indicator = 1.0;
  }
  else { /* otherwise in smoothed zone */
    indicator = exp (-d*d / (2.0*std_dev*std_dev));
  }
  RHEA_ASSERT (isfinite (indicator));
  RHEA_ASSERT (0.0 <= indicator && indicator <= 1.0);

  /* return weak zone indicator */
  return indicator;
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
  double              factor;

  /* check input */
  RHEA_ASSERT (isfinite (factor_interior));
  RHEA_ASSERT (0.0 < factor_interior && factor_interior <= 1.0);

  factor = 1.0 - (1.0 - factor_interior) *
           rhea_weakzone_indicator_node (distance, thickness, thickness_const);
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
  int                 label;
  double              distance;
  double              thickness, thickness_const, factor_interior;

  /* compute distance to a surface/manifold describing the weak zone */
  switch (opt->type) {
  case RHEA_WEAKZONE_DEPTH:
  case RHEA_WEAKZONE_DATA_POINTS:
    label = RHEA_WEAKZONE_LABEL_CLASS_NONE;
    distance = rhea_weakzone_dist_node (NULL, NULL, x, y, z, opt);
    thickness = rhea_weakzone_lookup_thickness (label, opt);
    thickness_const = rhea_weakzone_lookup_thickness_const (label, opt);
    factor_interior = rhea_weakzone_lookup_factor_interior (label, opt);
    break;
  case RHEA_WEAKZONE_DATA_POINTS_LABELS:
    label = RHEA_WEAKZONE_LABEL_UNKNOWN;
    distance = rhea_weakzone_dist_node (&label, NULL, x, y, z, opt);
    RHEA_ASSERT (rhea_weakzone_label_is_valid_int (label));
    thickness = rhea_weakzone_lookup_thickness (label, opt);
    thickness_const = rhea_weakzone_lookup_thickness_const (label, opt);
    factor_interior = rhea_weakzone_lookup_factor_interior (label, opt);
    break;
  case RHEA_WEAKZONE_DATA_POINTS_LABELS_FACTORS:
    label = RHEA_WEAKZONE_LABEL_UNKNOWN;
    distance = rhea_weakzone_dist_node (&label, &factor_interior, x, y, z, opt);
    RHEA_ASSERT (rhea_weakzone_label_is_valid_int (label));
    thickness = rhea_weakzone_lookup_thickness (label, opt);
    thickness_const = rhea_weakzone_lookup_thickness_const (label, opt);
    if (!isfinite (factor_interior)) { /* if factor is invalid */
      factor_interior = rhea_weakzone_lookup_factor_interior (label, opt);
    }
    break;
  default: /* unknown weak zone type */
    RHEA_ABORT_NOT_REACHED ();
  }
  RHEA_ASSERT (isfinite (distance));
  RHEA_ASSERT (isfinite (factor_interior));

  /* compute weak zone value */
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
  case RHEA_WEAKZONE_DEPTH:
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

  *dist = rhea_weakzone_dist_node (&label, NULL, x, y, z, opt);
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
  case RHEA_WEAKZONE_DEPTH:
  case RHEA_WEAKZONE_DATA_POINTS:
  case RHEA_WEAKZONE_DATA_POINTS_LABELS:
  case RHEA_WEAKZONE_DATA_POINTS_LABELS_FACTORS:
    ymir_dvec_set_function (distance, rhea_weakzone_dist_node_fn, opt);
    RHEA_ASSERT (rhea_weakzone_is_valid (distance));
    break;
  default: /* unknown weak zone type */
    RHEA_ABORT_NOT_REACHED ();
  }
}

typedef struct rhea_weakzone_indicator_node_data
{
  rhea_weakzone_options_t  *weak_options;
  rhea_weakzone_label_t     label_filter;
}
rhea_weakzone_indicator_node_data_t;

static void
rhea_weakzone_indicator_node_fn (double *indicator, double x, double y,
                                 double z, ymir_locidx_t nid, void *data)
{
  rhea_weakzone_indicator_node_data_t  *d = data;
  rhea_weakzone_options_t     *opt = d->weak_options;
  const rhea_weakzone_label_t label_filter = d->label_filter;
  rhea_weakzone_label_t       label;
  int                 label_int;
  double              distance, thickness, thickness_const, indc;

  /* get distance and weak zone parameters */
  switch (opt->type) {
  case RHEA_WEAKZONE_DEPTH:
  case RHEA_WEAKZONE_DATA_POINTS:
    label = RHEA_WEAKZONE_LABEL_CLASS_NONE;
    distance = rhea_weakzone_dist_node (NULL, NULL, x, y, z, opt);
    thickness = rhea_weakzone_lookup_thickness (label, opt);
    thickness_const = rhea_weakzone_lookup_thickness_const (label, opt);
    break;
  case RHEA_WEAKZONE_DATA_POINTS_LABELS:
  case RHEA_WEAKZONE_DATA_POINTS_LABELS_FACTORS:
    label_int = RHEA_WEAKZONE_LABEL_UNKNOWN;
    distance = rhea_weakzone_dist_node (&label_int, NULL, x, y, z, opt);
    RHEA_ASSERT (rhea_weakzone_label_is_valid_int (label_int));
    label = (rhea_weakzone_label_t) label_int;
    thickness = rhea_weakzone_lookup_thickness (label, opt);
    thickness_const = rhea_weakzone_lookup_thickness_const (label, opt);
    break;
  default: /* unknown weak zone type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* compute indicator value */
  indc = rhea_weakzone_indicator_node (distance, thickness, thickness_const);

  /* apply filter with respect to label */
  if (rhea_weakzone_label_is_class (label_filter)) {
    const rhea_weakzone_label_t class_id =
      rhea_weakzone_label_get_class ((rhea_weakzone_label_t) label);

    if (label_filter == class_id) *indicator = indc;
    else                          *indicator = 0.0;
  }
  else {
    if (label_filter == label) *indicator = indc;
    else                       *indicator = 0.0;
  }
}

void
rhea_weakzone_compute_indicator (ymir_vec_t *indicator,
                                 const int label_filter,
                                 rhea_weakzone_options_t *opt)
{
  rhea_weakzone_indicator_node_data_t data;

  /* check input */
  RHEA_ASSERT (opt->type == RHEA_WEAKZONE_NONE || indicator != NULL);
  RHEA_ASSERT (opt->type == RHEA_WEAKZONE_NONE ||
               rhea_weakzone_check_vec_type (indicator));
  RHEA_ASSERT (rhea_weakzone_label_is_valid_int (label_filter));

  /* set data */
  data.weak_options = opt;
  data.label_filter = (rhea_weakzone_label_t) label_filter;

  /* compute indicator */
  switch (opt->type) {
  case RHEA_WEAKZONE_NONE:
    if (indicator != NULL) {
      ymir_vec_set_value (indicator, 0.0);
    }
    break;
  case RHEA_WEAKZONE_DEPTH:
  case RHEA_WEAKZONE_DATA_POINTS:
  case RHEA_WEAKZONE_DATA_POINTS_LABELS:
  case RHEA_WEAKZONE_DATA_POINTS_LABELS_FACTORS:
    ymir_dvec_set_function (indicator, rhea_weakzone_indicator_node_fn, &data);
    break;
  default: /* unknown weak zone type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* check output */
  RHEA_ASSERT (sc_dmatrix_is_valid (indicator->dataown) &&
               0.0 <= ymir_dvec_min_global (indicator) &&
               ymir_dvec_max_global (indicator) <= 1.0);
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
