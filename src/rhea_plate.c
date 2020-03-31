#include <rhea_plate.h>
#include <rhea_base.h>
#include <rhea_point_in_polygon.h>
#include <rhea_io_mpi.h>
#include <rhea_velocity.h>
#include <ymir_velocity_vec.h>
#include <ymir_mass_vec.h>
#include <ymir_perf_counter.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_PLATE_DEFAULT_POLYGON_N_POLYGONS (0)
#define RHEA_PLATE_DEFAULT_POLYGON_VERTICES_FILE_PATH_TXT NULL
#define RHEA_PLATE_DEFAULT_POLYGON_N_VERTICES_TOTAL (0)
#define RHEA_PLATE_DEFAULT_POLYGON_VERTICES_COARSE_FILE_PATH_TXT NULL
#define RHEA_PLATE_DEFAULT_POLYGON_N_VERTICES_COARSE_TOTAL (0)
#define RHEA_PLATE_DEFAULT_POLYGON_X_MIN (0.0)
#define RHEA_PLATE_DEFAULT_POLYGON_X_MAX (360.0)
#define RHEA_PLATE_DEFAULT_POLYGON_Y_MIN (0.0)
#define RHEA_PLATE_DEFAULT_POLYGON_Y_MAX (180.0)
#define RHEA_PLATE_DEFAULT_XSECTION_BOUNDARY_LON_LIST NULL
#define RHEA_PLATE_DEFAULT_XSECTION_VELOCITY_MM_YR_LIST NULL
#define RHEA_PLATE_DEFAULT_XSECTION_SHRINK_FACTOR (NAN)
#define RHEA_PLATE_DEFAULT_MONITOR_PERFORMANCE (0)

/* initialize options */
int                 rhea_plate_polygon_n_polygons =
  RHEA_PLATE_DEFAULT_POLYGON_N_POLYGONS;
char               *rhea_plate_polygon_vertices_file_path_txt =
  RHEA_PLATE_DEFAULT_POLYGON_VERTICES_FILE_PATH_TXT;
int                 rhea_plate_polygon_n_vertices_total =
  RHEA_PLATE_DEFAULT_POLYGON_N_VERTICES_TOTAL;
char               *rhea_plate_polygon_vertices_coarse_file_path_txt =
  RHEA_PLATE_DEFAULT_POLYGON_VERTICES_COARSE_FILE_PATH_TXT;
int                 rhea_plate_polygon_n_vertices_coarse_total =
  RHEA_PLATE_DEFAULT_POLYGON_N_VERTICES_COARSE_TOTAL;
double              rhea_plate_polygon_x_min =
  RHEA_PLATE_DEFAULT_POLYGON_X_MIN;
double              rhea_plate_polygon_x_max =
  RHEA_PLATE_DEFAULT_POLYGON_X_MAX;
double              rhea_plate_polygon_y_min =
  RHEA_PLATE_DEFAULT_POLYGON_Y_MIN;
double              rhea_plate_polygon_y_max =
  RHEA_PLATE_DEFAULT_POLYGON_Y_MAX;
char               *rhea_plate_xsection_boundary_lon_list =
  RHEA_PLATE_DEFAULT_XSECTION_BOUNDARY_LON_LIST;
char               *rhea_plate_xsection_velocity_mm_yr_list =
  RHEA_PLATE_DEFAULT_XSECTION_VELOCITY_MM_YR_LIST;
double              rhea_plate_xsection_shrink_factor =
  RHEA_PLATE_DEFAULT_XSECTION_SHRINK_FACTOR;
int                 rhea_plate_monitor_performance =
                      RHEA_PLATE_DEFAULT_MONITOR_PERFORMANCE;

void
rhea_plate_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Plate";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  /*** polygon based (2D) plates ***/

  YMIR_OPTIONS_I, "num-polygons", '\0',
    &(rhea_plate_polygon_n_polygons), RHEA_PLATE_DEFAULT_POLYGON_N_POLYGONS,
    "Number of polygons, each of which represents one plate",

  YMIR_OPTIONS_S, "polygon-vertices-file-path-txt", '\0',
    &(rhea_plate_polygon_vertices_file_path_txt),
    RHEA_PLATE_DEFAULT_POLYGON_VERTICES_FILE_PATH_TXT,
    "Path to a text file with (lon,lat) vertices of plate polygons",
  YMIR_OPTIONS_I, "polygon-vertices-num-total", '\0',
    &(rhea_plate_polygon_n_vertices_total),
    RHEA_PLATE_DEFAULT_POLYGON_N_VERTICES_TOTAL,
    "Total number of vertices of all combined plate polygons",

  YMIR_OPTIONS_S, "polygon-vertices-coarse-file-path-txt", '\0',
    &(rhea_plate_polygon_vertices_coarse_file_path_txt),
    RHEA_PLATE_DEFAULT_POLYGON_VERTICES_COARSE_FILE_PATH_TXT,
    "Path to a text file with vertices to coarse polygon containers",
  YMIR_OPTIONS_I, "polygon-vertices-coarse-num-total", '\0',
    &(rhea_plate_polygon_n_vertices_coarse_total),
    RHEA_PLATE_DEFAULT_POLYGON_N_VERTICES_COARSE_TOTAL,
    "Total number of vertices of all combined coarse polygon containers",

  YMIR_OPTIONS_D, "polygon-x-min", '\0',
    &(rhea_plate_polygon_x_min), RHEA_PLATE_DEFAULT_POLYGON_X_MIN,
    "Polygon vertices: Min x",
  YMIR_OPTIONS_D, "polygon-x-max", '\0',
    &(rhea_plate_polygon_x_max), RHEA_PLATE_DEFAULT_POLYGON_X_MAX,
    "Polygon vertices: Max x",
  YMIR_OPTIONS_D, "polygon-y-min", '\0',
    &(rhea_plate_polygon_y_min), RHEA_PLATE_DEFAULT_POLYGON_Y_MIN,
    "Polygon vertices: Min y",
  YMIR_OPTIONS_D, "polygon-y-max", '\0',
    &(rhea_plate_polygon_y_max), RHEA_PLATE_DEFAULT_POLYGON_Y_MAX,
    "Polygon vertices: Max y",

  /*** interval based (1D) plates ***/

  YMIR_OPTIONS_S, "cross-section-plate-boundary-lon-list", '\0',
    &(rhea_plate_xsection_boundary_lon_list),
    RHEA_PLATE_DEFAULT_XSECTION_BOUNDARY_LON_LIST,
    "Cross section domain: Plate boundaries (comma separated list of lon's)",
  YMIR_OPTIONS_S, "cross-section-plate-velocity-mm-yr-list", '\0',
    &(rhea_plate_xsection_velocity_mm_yr_list),
    RHEA_PLATE_DEFAULT_XSECTION_VELOCITY_MM_YR_LIST,
    "Cross section domain: Plate velocities [mm/yr] (comma separated list)",
  YMIR_OPTIONS_D, "cross-section-plate-shrink-factor", '\0',
    &(rhea_plate_xsection_shrink_factor),
    RHEA_PLATE_DEFAULT_XSECTION_SHRINK_FACTOR,
    "Cross section domain: Shrink plates by this factor in (0,1)",

  /*** other options ***/

  YMIR_OPTIONS_B, "monitor-performance", '\0',
    &(rhea_plate_monitor_performance), RHEA_PLATE_DEFAULT_MONITOR_PERFORMANCE,
    "Measure and print performance statistics (e.g., runtime or flops)",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);

  /* initialize (deactivated) performance counters */
  rhea_plate_perfmon_init (0, 0);
}

void
rhea_plate_process_options (rhea_plate_options_t *opt,
                            rhea_domain_options_t *domain_options,
                            rhea_temperature_options_t *temp_options)
{
  /* exit if nothing to do */
  if (opt == NULL) {
    return;
  }

  /* set options for cross sectional domain */
  opt->xsection_boundary_lon_list = rhea_plate_xsection_boundary_lon_list;
  opt->xsection_boundary = NULL;
  opt->xsection_tangential_velocity_mm_yr_list =
    rhea_plate_xsection_velocity_mm_yr_list;
  opt->xsection_tangential_velocity = NULL;
  opt->xsection_n_intervals = 0;
  opt->xsection_shrink_factor = rhea_plate_xsection_shrink_factor;

  /* set number of polygons */
  if (0 < rhea_plate_polygon_n_polygons) { /* if set from option */
    opt->n_polygons = rhea_plate_polygon_n_polygons;
  }
  else { /* otherwise set to default value */
    switch (domain_options->shape) {
    case RHEA_DOMAIN_CUBE:
      opt->n_polygons = RHEA_PLATE_CUBE_N;
      break;
    case RHEA_DOMAIN_BOX:
      opt->n_polygons = 0;
      break;
    case RHEA_DOMAIN_SHELL:
      opt->n_polygons = RHEA_PLATE_EARTH_N;
      break;
    case RHEA_DOMAIN_CUBE_SPHERICAL:
      opt->n_polygons = 0;
      break;
    case RHEA_DOMAIN_BOX_SPHERICAL:
      opt->n_polygons = 0;
      break;
    default: /* unknown domain shape */
      RHEA_ABORT_NOT_REACHED ();
    }
  }

  /* set paths to text files */
  opt->vertices_file_path_txt = rhea_plate_polygon_vertices_file_path_txt;
  opt->n_vertices_total = rhea_plate_polygon_n_vertices_total;
  opt->vertices_coarse_container_file_path_txt =
    rhea_plate_polygon_vertices_coarse_file_path_txt;
  opt->n_vertices_coarse_total = rhea_plate_polygon_n_vertices_coarse_total;

  /* initialize storage variables */
  opt->vertices_x = NULL;
  opt->vertices_y = NULL;
  opt->n_vertices = NULL;
  opt->vertices_coarse_container_x = NULL;
  opt->vertices_coarse_container_y = NULL;
  opt->n_vertices_coarse_container = NULL;
  opt->translation_x = NULL;
  opt->translation_y = NULL;

  /* initialize range of values */
  opt->x_min = (float) rhea_plate_polygon_x_min;
  opt->x_max = (float) rhea_plate_polygon_x_max;
  opt->y_min = (float) rhea_plate_polygon_y_min;
  opt->y_max = (float) rhea_plate_polygon_y_max;

  /* initialize plate velocities */
  opt->angular_velocity = NULL;

  /* set dependent options */
  opt->domain_options = domain_options;
  opt->temp_options = temp_options;
}

/******************************************************************************
 * Monitoring
 *****************************************************************************/

/* perfomance monitor tags and names */
typedef enum
{
  RHEA_PLATE_PERFMON_DATA_CREATE,
  RHEA_PLATE_PERFMON_POINT_INSIDE_POLYGON,
  RHEA_PLATE_PERFMON_SET_LABEL,
  RHEA_PLATE_PERFMON_APPLY_FILTER,
  RHEA_PLATE_PERFMON_N
}
rhea_plate_perfmon_idx_t;

static const char  *rhea_plate_perfmon_name[RHEA_PLATE_PERFMON_N] =
{
  "Create data",
  "Point inside polygon test",
  "Set plate label",
  "Apply plate filter"
};
ymir_perf_counter_t rhea_plate_perfmon[RHEA_PLATE_PERFMON_N];

void
rhea_plate_perfmon_init (const int activate, const int skip_if_active)
{
  const int           active = activate && rhea_plate_monitor_performance;

  ymir_perf_counter_init_all_ext (rhea_plate_perfmon,
                                  rhea_plate_perfmon_name,
                                  RHEA_PLATE_PERFMON_N,
                                  active, skip_if_active);
}

void
rhea_plate_perfmon_print (sc_MPI_Comm mpicomm,
                          const int print_wtime,
                          const int print_n_calls,
                          const int print_flops)
{
  const int           active = rhea_plate_monitor_performance;
  const int           print = (print_wtime || print_n_calls || print_flops);
  int                 n_stats = RHEA_PLATE_PERFMON_N *
                                YMIR_PERF_COUNTER_N_STATS;
  sc_statinfo_t       stats[n_stats];
  char                stats_name[n_stats][YMIR_PERF_COUNTER_NAME_SIZE];

  /* exit if nothing to do */
  if (!active || !print) {
    return;
  }

  /* gather performance statistics */
  n_stats = ymir_perf_counter_gather_stats (
      rhea_plate_perfmon, RHEA_PLATE_PERFMON_N, stats, stats_name,
      mpicomm, print_wtime, print_n_calls, print_flops);

  /* print performance statistics */
  ymir_perf_counter_print_stats (stats, n_stats, "Plate");
}

/******************************************************************************
 * Plate Boundary Data
 *****************************************************************************/

static int
rhea_plate_polygon_is_separator (const float vert_x, const float vert_y)
{
  return (!isfinite (vert_x) && !isfinite (vert_y));
}

static int
rhea_plate_polygon_get_n_vertices (const float *vertices_all,
                                   const int n_vertices_total)
{
  int                 k;

  for (k = 0; k < n_vertices_total; k++) {
    const float         vx = vertices_all[2*k    ];
    const float         vy = vertices_all[2*k + 1];

    /* stop loop if at the beginning of a new polygon */
    if (rhea_plate_polygon_is_separator (vx, vy)) {
      break;
    }
  }

  return k;
}

static void
rhea_plate_process_vertices_cube (float **vertices_x,
                                  float **vertices_y,
                                  size_t *n_vertices,
                                  const int n_plates,
                                  const float *vertices_all,
                                  const int n_vertices_total,
                                  const float x_min, const float x_max,
                                  const float y_min, const float y_max)
{
  int                 idx_curr = 0;
  int                 n_total = n_vertices_total;
  int                 pid;
  float              *vert_x, *vert_y;
  size_t              n_vert, vid;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (__func__, "n_plates=%i, n_vertices_total=%i",
                                 n_plates, n_vertices_total);

  /* skip first */
  if (rhea_plate_polygon_is_separator (vertices_all[0], vertices_all[1])) {
    idx_curr++;
    n_total--;
  }

  for (pid = 0; pid < n_plates; pid++) {
    RHEA_ASSERT (idx_curr < n_vertices_total);
    RHEA_ASSERT (0 < n_total);

    /* get #vertices for current plate `pid` */
    n_vert = n_vertices[pid] = (size_t) rhea_plate_polygon_get_n_vertices (
        vertices_all + 2*idx_curr, n_total);
    RHEA_GLOBAL_VERBOSEF_FN_TAG (__func__, "plate_idx=%i, n_vertices=%i",
                                 pid, (int) n_vert);
    RHEA_ASSERT (0 < n_vert);

    /* set vertex coordinates */
    vert_x = vertices_x[pid] = RHEA_ALLOC (float, n_vert);
    vert_y = vertices_y[pid] = RHEA_ALLOC (float, n_vert);
    for (vid = 0; vid < n_vert; vid++) {
      vert_x[vid] = vertices_all[2*(idx_curr + vid)    ];
      vert_y[vid] = vertices_all[2*(idx_curr + vid) + 1];
      RHEA_ASSERT (isfinite (vert_x[vid]));
      RHEA_ASSERT (isfinite (vert_y[vid]));

      vert_x[vid] = (vert_x[vid] - x_min) / (x_max - x_min);
      vert_y[vid] = (vert_y[vid] - y_min) / (y_max - y_min);
      RHEA_ASSERT (0.0 <= vert_x[vid] && vert_x[vid] <= 1.0);
      RHEA_ASSERT (0.0 <= vert_y[vid] && vert_y[vid] <= 1.0);
    }

    /* increment vertex counter */
    idx_curr += (n_vert + 1);
    n_total -= (n_vert + 1);
  }

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static void
rhea_plate_modify_vertices_shell_init (int bins_x[36],
                                       int bins_y[180],
                                       const float *vertices_x,
                                       const float *vertices_y,
                                       const size_t n_vertices)
{
  size_t              vid;
  int                 bx, by;

  /* init bins to zero */
  memset (bins_x, 0, 36 * sizeof (int));
  memset (bins_y, 0, 180 * sizeof (int));

  /* add vertex coordinates to bins */
  for (vid = 0; vid < n_vertices; vid++) {
    bx = SC_MIN (SC_MAX (0, (int) floor (vertices_x[vid]/10.0)), 35);
    by = SC_MIN (SC_MAX (0, (int) floor (vertices_y[vid])), 179);
    bins_x[bx]++;
    bins_y[by]++;
  }
}

static void
rhea_plate_modify_vertices_shell_set (float *translation_x,
                                      float *translation_y,
                                      float *insert_y,
                                      const int bins_x[36],
                                      const int bins_y[180])
{
  int                 bin_east, bin_west, bin_south, bin_north;
  int                 i;

  /* find east & west bounds */
  bin_east = -1;
  for (i = 0; i <= 35; i++) { /* loop east (counter clock-wise) */
    if (0 < bins_x[i]) {
      bin_east = i;
    }
    else if (0 <= bin_east) {
      break;
    }
  }
  bin_west = -1;
  for (i = 35; 0 <= i; i--) { /* loop west (clock-wise) */
    if (0 < bins_x[i]) {
      bin_west = i;
    }
    else if (0 <= bin_west) {
      break;
    }
  }

  /* find north & south bounds */
  for (i = 0; i <= 179; i++) { /* loop south (down) */
    if (0 < bins_y[i]) {
      bin_north = i;
      break;
    }
  }
  for (i = 179; 0 <= i; i--) { /* loop north (up) */
    if (0 < bins_y[i]) {
      bin_south = i;
      break;
    }
  }
  if (0 == bin_west && bin_east == 35) { /* if whole circle */
    if (bin_north < 90 && bin_south < 90) { /* if polygon in north hemisph. */
      bin_north = 0;
    }
    else if (90 < bin_north && 90 < bin_south) { /* if polygon in south hemi. */
      bin_south = 179;
    }
    else {
      RHEA_ABORT_NOT_REACHED ();
    }
  }

  /* set translations */
  if (bin_west <= bin_east && 0 < bin_north && bin_south < 179) { /* no wrap */
    *translation_x = 0.0;
    *translation_y = 0.0;
    *insert_y = NAN;
  }
  else if (0 < bin_north && bin_south < 179) { /* set translation in west dir */
    RHEA_ASSERT (bin_east < bin_west);
    *translation_x = -10.0 * (float) bin_west;
    *translation_y = 0.0;
    *insert_y = NAN;
  }
  else if (0 == bin_north) { /* insert north pole (0,0) and (360,0) */
    RHEA_ASSERT (bin_south < 179);
    *translation_x = 0.0;
    *translation_y = 0.0;
    *insert_y = 0.0;
  }
  else if (bin_south == 179) { /* insert south pole (0,180) and (360,180) */
    RHEA_ASSERT (0 < bin_north);
    *translation_x = 0.0;
    *translation_y = 0.0;
    *insert_y = 180.0;
  }
  else {
    RHEA_ABORT_NOT_REACHED ();
  }

  RHEA_GLOBAL_VERBOSEF_FN_TAG (
      __func__, "west=%i, east=%i, north=%i, south=%i, "
      "translation_x=%g, translation_y=%g, insert_y=%g",
      10*bin_west, 10*bin_east, bin_north, bin_south,
      *translation_x, *translation_y, *insert_y);
}

static void
rhea_plate_modify_vertices_shell_translate (float *vertices_x,
                                            float *vertices_y,
                                            const size_t n_vertices,
                                            const float translation_x,
                                            const float translation_y)
{
  size_t              vid;

  /* apply translation into east-west direction */
  if (0.0 < fabs (translation_x)) {
    for (vid = 0; vid < n_vertices; vid++) {
      if (0.0 <= (vertices_x[vid] + translation_x)) {
        vertices_x[vid] += translation_x;
      }
      else {
        vertices_x[vid] += translation_x + 360.0;
      }
      RHEA_ASSERT (0.0 <= vertices_x[vid] && vertices_x[vid] <= 360.0);
    }
  }

  /* apply translation into south-north direction */
  if (0.0 < fabs (translation_y)) {
    for (vid = 0; vid < n_vertices; vid++) {
      vertices_y[vid] += translation_y;
      RHEA_ASSERT (0.0 <= vertices_y[vid] && vertices_y[vid] <= 180.0);
    }
  }
}

static size_t
rhea_plate_modify_vertices_shell_insert (float **vertices_x,
                                         float **vertices_y,
                                         const size_t n_vertices,
                                         const float insert_y)
{
  float              *vertices_old_x = *vertices_x;
  float              *vertices_old_y = *vertices_y;
  float               x, y, prev_x, prev_y;
  size_t              vid;

  /* create larger arrays */
  *vertices_x = RHEA_ALLOC (float, n_vertices + 4);
  *vertices_y = RHEA_ALLOC (float, n_vertices + 4);

  /* copy vertices until wrapping occurs */
  prev_x = vertices_old_x[0];
  prev_y = vertices_old_y[0];
  for (vid = 0; vid < n_vertices; vid++) {
    x = vertices_old_x[vid];
    y = vertices_old_y[vid];
    if (fabs (x - prev_x) <= 10.0) {
      (*vertices_x)[vid] = x;
      (*vertices_y)[vid] = y;
      prev_x = x;
      prev_y = y;
    }
    else {
      break;
    }
  }

  /* insert new vertex coordinates */
  if (x <= prev_x) { /* if wraps counter-clockwise */
    (*vertices_x)[vid  ] = 360.0;
    (*vertices_y)[vid  ] = prev_y;
    (*vertices_x)[vid+1] = 360.0;
    (*vertices_y)[vid+1] = insert_y;
    (*vertices_x)[vid+2] = 0.0;
    (*vertices_y)[vid+2] = insert_y;
    (*vertices_x)[vid+3] = 0.0;
    (*vertices_y)[vid+3] = prev_y;
  }
  else { /* otherwise wraps clockwise */
    (*vertices_x)[vid  ] = 0.0;
    (*vertices_y)[vid  ] = prev_y;
    (*vertices_x)[vid+1] = 0.0;
    (*vertices_y)[vid+1] = insert_y;
    (*vertices_x)[vid+2] = 360.0;
    (*vertices_y)[vid+2] = insert_y;
    (*vertices_x)[vid+3] = 360.0;
    (*vertices_y)[vid+3] = prev_y;
  }

  /* continue copying vertices */
  prev_x = x;
  prev_y = y;
  for (/* resume */; vid < n_vertices; vid++) {
    x = vertices_old_x[vid];
    y = vertices_old_y[vid];
    RHEA_ASSERT (fabs (x - prev_x) <= 10.0);
    (*vertices_x)[vid+4] = x;
    (*vertices_y)[vid+4] = y;
    prev_x = x;
    prev_y = y;
  }

  /* destroy old arrays */
  RHEA_FREE (vertices_old_x);
  RHEA_FREE (vertices_old_y);

  /* return new #vertices */
  return n_vertices + 4;
}

static void
rhea_plate_modify_vertices_shell (float *translation_x,
                                  float *translation_y,
                                  float **vertices_x,
                                  float **vertices_y,
                                  size_t *n_vertices,
                                  const int n_plates)
{
  int                 bins_x[36], bins_y[180];
  int                 pid;
  float               insert_y;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (__func__, "n_plates=%i", n_plates);

  for (pid = 0; pid < n_plates; pid++) {
    /* populate bins */
    rhea_plate_modify_vertices_shell_init (
        bins_x, bins_y, vertices_x[pid], vertices_y[pid], n_vertices[pid]);

    /* set translations */
    rhea_plate_modify_vertices_shell_set (
        &translation_x[pid], &translation_y[pid], &insert_y, bins_x, bins_y);

    /* apply translations */
    rhea_plate_modify_vertices_shell_translate (
        vertices_x[pid], vertices_y[pid], n_vertices[pid],
        translation_x[pid], translation_y[pid]);

    /* insert vertices */
    if (isfinite (insert_y)) {
      n_vertices[pid] = rhea_plate_modify_vertices_shell_insert (
          &vertices_x[pid], &vertices_y[pid], n_vertices[pid], insert_y);
    }
  }

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static void
rhea_plate_process_vertices_shell (float **vertices_x,
                                   float **vertices_y,
                                   size_t *n_vertices,
                                   const int n_plates,
                                   const float *vertices_all,
                                   const int n_vertices_total,
                                   const float x_min, const float x_max,
                                   const float y_min, const float y_max)
{
  int                 idx_curr = 0;
  int                 n_total = n_vertices_total;
  int                 pid;
  float              *vert_x, *vert_y;
  size_t              n_vert, vid;
#ifdef RHEA_ENABLE_DEBUG
  float               vert_x_min, vert_x_max;
  float               vert_y_min, vert_y_max;
#endif

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (__func__, "n_plates=%i, n_vertices_total=%i",
                                 n_plates, n_vertices_total);

  /* skip first */
  if (rhea_plate_polygon_is_separator (vertices_all[0], vertices_all[1])) {
    idx_curr++;
    n_total--;
  }

  for (pid = 0; pid < n_plates; pid++) {
    RHEA_ASSERT (idx_curr < n_vertices_total);
    RHEA_ASSERT (0 < n_total);

    /* get #vertices for current plate `pid` */
    n_vert = n_vertices[pid] = (size_t) rhea_plate_polygon_get_n_vertices (
        vertices_all + 2*idx_curr, n_total);
    RHEA_GLOBAL_VERBOSEF_FN_TAG (__func__, "plate_idx=%i, n_vertices=%i",
                                 pid, (int) n_vert);
    RHEA_ASSERT (0 < n_vert);

    /* set vertex coordinates */
    vert_x = vertices_x[pid] = RHEA_ALLOC (float, n_vert);
    vert_y = vertices_y[pid] = RHEA_ALLOC (float, n_vert);
#ifdef RHEA_ENABLE_DEBUG
    vert_x_min = 360.0; vert_x_max = 0.0;
    vert_y_min = 180.0; vert_y_max = 0.0;
#endif
    for (vid = 0; vid < n_vert; vid++) {
      vert_x[vid] = vertices_all[2*(idx_curr + vid)    ];
      vert_y[vid] = vertices_all[2*(idx_curr + vid) + 1];
      RHEA_ASSERT (isfinite (vert_x[vid]));
      RHEA_ASSERT (isfinite (vert_y[vid]));

      vert_x[vid] = 360.0 * (vert_x[vid] - x_min) / (x_max - x_min);
      vert_y[vid] = 180.0 * (vert_y[vid] - y_min) / (y_max - y_min);
      RHEA_ASSERT (-1.0e8*SC_EPS < vert_x[vid] &&
                   vert_x[vid] < 360.0 + 1.0e8*SC_EPS);
      RHEA_ASSERT (-1.0e8*SC_EPS < vert_y[vid] &&
                   vert_y[vid] < 180.0 + 1.0e8*SC_EPS);
      vert_x[vid] = SC_MIN (SC_MAX (0.0, vert_x[vid]), 360.0);
      vert_y[vid] = SC_MIN (SC_MAX (0.0, vert_y[vid]), 180.0);
#ifdef RHEA_ENABLE_DEBUG
      vert_x_min = SC_MIN (vert_x[vid], vert_x_min);
      vert_x_max = SC_MAX (vert_x[vid], vert_x_max);
      vert_y_min = SC_MIN (vert_y[vid], vert_y_min);
      vert_y_max = SC_MAX (vert_y[vid], vert_y_max);
#endif
    }
#ifdef RHEA_ENABLE_DEBUG
    RHEA_GLOBAL_VERBOSEF_FN_TAG (
        __func__, "plate_idx=%i, "
        "vert_x_min=%g, vert_x_max=%g, vert_y_min=%g, vert_y_max=%g",
        pid, vert_x_min, vert_x_max, vert_y_min, vert_y_max);
#endif

    /* increment vertex counter */
    idx_curr += (n_vert + 1);
    n_total -= (n_vert + 1);
  }

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

/******************************************************************************
 * Plate Velocity Data
 *****************************************************************************/

/* Cartesian angular velocities.  For each plate we have:
 *   [0] x
 *   [1] y
 *   [2] z
 */
static const double rhea_plate_velocity_cube[3*RHEA_PLATE_CUBE_N] =
{
  /*RHEA_PLATE_CUBE_NW*/ 0.0, 0.0,  1.0,
  /*RHEA_PLATE_CUBE_NE*/ 0.0, 0.0, -1.0,
  /*RHEA_PLATE_CUBE_SE*/ 0.0, 0.0,  1.0,
  /*RHEA_PLATE_CUBE_SW*/ 0.0, 0.0, -1.0
};

/* Spherical angular velocities from NNR-MORVEL56.  For each plate we have:
 *   [0] longitude        -180 <   phi  <= +180
 *   [1] latitude         - 90 <= theta <= + 90
 *   [2] angular velocity    0 <= omega
 */
static const double rhea_plate_velocity_earth_morvel56[3*RHEA_PLATE_EARTH_N] =
{
  /*RHEA_PLATE_EARTH_AM*/  63.1704, -122.8242, 0.2973,
  /*RHEA_PLATE_EARTH_AN*/  65.4235, -118.1053, 0.2500,
  /*RHEA_PLATE_EARTH_AR*/  48.8807, -  8.4909, 0.5588,
  /*RHEA_PLATE_EARTH_AU*/  33.8612,   37.9414, 0.6316,
  /*RHEA_PLATE_EARTH_CP*/  44.4352,   23.0880, 0.6080,
  /*RHEA_PLATE_EARTH_CA*/  35.1956, - 92.6236, 0.2862,
  /*RHEA_PLATE_EARTH_CO*/  26.9346, -124.3074, 1.1978,
  /*RHEA_PLATE_EARTH_EU*/  48.8509, -106.5007, 0.2227,
  /*RHEA_PLATE_EARTH_IN*/  50.3722, -  3.2898, 0.5438,
  /*RHEA_PLATE_EARTH_JF*/ -38.3086,   60.0379, 0.9513,
  /*RHEA_PLATE_EARTH_LW*/  51.8860, - 69.5195, 0.2856,
  /*RHEA_PLATE_EARTH_MQ*/  49.1891,   11.0524, 1.1440,
  /*RHEA_PLATE_EARTH_NZ*/  46.2348, -101.0564, 0.6957,
  /*RHEA_PLATE_EARTH_NA*/ - 4.8548, - 80.6447, 0.2087,
///*RHEA_PLATE_EARTH_NU*/  47.6763, - 68.4377, 0.2921,
  /*RHEA_PLATE_EARTH_PA*/ -63.5756,  114.6975, 0.6509,
  /*RHEA_PLATE_EARTH_PS*/ -46.0242, - 31.3615, 0.9098,
  /*RHEA_PLATE_EARTH_RI*/  20.2450, -107.2861, 4.5359,
  /*RHEA_PLATE_EARTH_SW*/ -29.9420, - 36.8671, 1.3616,
  /*RHEA_PLATE_EARTH_SC*/  22.5244, -106.1485, 0.1464,
  /*RHEA_PLATE_EARTH_SM*/  49.9506, - 84.5154, 0.3393,
  /*RHEA_PLATE_EARTH_SA*/ -22.6179, -112.8327, 0.1090,
  /*RHEA_PLATE_EARTH_SU*/  50.0558, - 95.0218, 0.3368,
  /*RHEA_PLATE_EARTH_SR*/ -32.4957, -111.3224, 0.1072,
  /*RHEA_PLATE_EARTH_YZ*/  63.0285, -116.6180, 0.3335,

  /*RHEA_PLATE_EARTH_AS*/  19.4251,  122.8665, 0.1239,
  /*RHEA_PLATE_EARTH_AP*/ - 6.5763, - 83.9776, 0.4881,
  /*RHEA_PLATE_EARTH_AT*/  40.1121,   26.6585, 1.2105,
  /*RHEA_PLATE_EARTH_BR*/ -63.7420,  142.0636, 0.4898,
  /*RHEA_PLATE_EARTH_BS*/ - 1.4855,  121.6413, 2.4753,
  /*RHEA_PLATE_EARTH_BH*/ -39.9983,  100.4994, 0.7988,
  /*RHEA_PLATE_EARTH_BU*/ - 6.1254, - 78.1008, 2.2287,
  /*RHEA_PLATE_EARTH_CL*/ -72.7849,   72.0525, 0.6066,
///*RHEA_PLATE_EARTH_CR*/ -20.3985,  170.5303, 3.9232,
  /*RHEA_PLATE_EARTH_EA*/  24.9729,   67.5269, 11.334,
///*RHEA_PLATE_EARTH_FT*/ -16.3322,  178.0679, 5.1006,
///*RHEA_PLATE_EARTH_GP*/   2.5287,   81.1806, 5.4868,
///*RHEA_PLATE_EARTH_JZ*/  34.2507,   70.7429, 22.368,
  /*RHEA_PLATE_EARTH_KE*/  39.9929,    6.4584, 2.3474,
///*RHEA_PLATE_EARTH_MN*/ - 3.6699,  150.2676, 51.569,
///*RHEA_PLATE_EARTH_MO*/  14.2480,   92.6656, 0.7742,
  /*RHEA_PLATE_EARTH_MA*/  11.0533,  137.8404, 1.3061,
  /*RHEA_PLATE_EARTH_MS*/   2.1477, - 56.0916, 3.5655,
  /*RHEA_PLATE_EARTH_NH*/   0.5684, -  6.6018, 2.4688,
///*RHEA_PLATE_EARTH_NI*/ - 3.2883, -174.4882, 3.3136,
  /*RHEA_PLATE_EARTH_ND*/  17.7331, -122.6815, 0.1162,
  /*RHEA_PLATE_EARTH_NB*/ -45.0406,  127.6370, 0.8563,
  /*RHEA_PLATE_EARTH_OK*/  30.3022, - 92.2813, 0.2290,
  /*RHEA_PLATE_EARTH_ON*/  36.1163,  137.9182, 2.5391,
  /*RHEA_PLATE_EARTH_PM*/  31.3510, -113.9038, 0.3171,
///*RHEA_PLATE_EARTH_SL*/  50.7058, -143.4675, 0.2677,
  /*RHEA_PLATE_EARTH_SS*/ - 2.8685,  130.6236, 1.7029,
  /*RHEA_PLATE_EARTH_SB*/   6.8767, - 31.8883, 8.1107,
  /*RHEA_PLATE_EARTH_TI*/ - 4.4363,  113.4976, 1.8639,
  /*RHEA_PLATE_EARTH_TO*/  25.8737,    4.4767, 8.9417,
  /*RHEA_PLATE_EARTH_WL*/   0.1050,  128.5186, 1.7444
};

static void
rhea_plate_velocity_earth_apply_dim_scal (double rot_spherical[3],
                                          const int remove_dim,
                                          rhea_plate_options_t *opt)
{
  const double  sec_per_yr = RHEA_TEMPERATURE_SECONDS_PER_YEAR;
  const double  radius_m = opt->domain_options->radius_max_m;
  const double  thermal_diffus_m2_s =
                  opt->temp_options->thermal_diffusivity_m2_s;

  if (!remove_dim) {
    RHEA_ASSERT (0.0 <= rot_spherical[0]);
    RHEA_ASSERT (0.0 <= rot_spherical[1] && rot_spherical[1] <= 2.0*M_PI);
    RHEA_ASSERT (0.0 <= rot_spherical[2] && rot_spherical[2] <= M_PI);
    rot_spherical[0] /= M_PI / 180.0;
    rot_spherical[0] /= radius_m / (thermal_diffus_m2_s/radius_m);
    rot_spherical[0] /= 1.0 / (1.0e6*sec_per_yr);
    rot_spherical[1] = rot_spherical[1]/M_PI * 180.0 - 180.0;
    rot_spherical[2] = -rot_spherical[2]/M_PI * 180.0 + 90.0;
  }
  else {
    RHEA_ASSERT (0.0 <= rot_spherical[0]);
    RHEA_ASSERT (-180.0 <= rot_spherical[1] && rot_spherical[1] <= 180.0);
    RHEA_ASSERT (- 90.0 <= rot_spherical[2] && rot_spherical[2] <=  90.0);
    rot_spherical[0] *= M_PI / 180.0;
    rot_spherical[0] *= radius_m / (thermal_diffus_m2_s/radius_m);
    rot_spherical[0] *= 1.0 / (1.0e6*sec_per_yr);
    rot_spherical[1] = (rot_spherical[1] + 180.0)/180.0 * M_PI;
    rot_spherical[2] = -(rot_spherical[2] - 90.0)/180.0 * M_PI;
  }
}

/**
 * Assigns angular velocity of plate `plate_label`.
 */
static void
rhea_plate_velocity_lookup_angular_velocity (double rot_axis[3],
                                             const int plate_label,
                                             rhea_plate_options_t *opt)
{
  /* check input */
  RHEA_ASSERT (RHEA_PLATE_NONE < plate_label);

  switch (opt->domain_options->shape) {
  case RHEA_DOMAIN_CUBE:
    /* get values */
    RHEA_ASSERT (plate_label < RHEA_PLATE_CUBE_N);
    rot_axis[0] = rhea_plate_velocity_cube[3*plate_label    ];
    rot_axis[1] = rhea_plate_velocity_cube[3*plate_label + 1];
    rot_axis[2] = rhea_plate_velocity_cube[3*plate_label + 2];
    break;

  case RHEA_DOMAIN_SHELL:
    {
      double        rot_sp[3];

      /* get values (in spherical coordinates) */
      RHEA_ASSERT (plate_label < RHEA_PLATE_EARTH_N);
      rot_sp[0] = rhea_plate_velocity_earth_morvel56[3*plate_label + 2];
      rot_sp[1] = rhea_plate_velocity_earth_morvel56[3*plate_label + 1];
      rot_sp[2] = rhea_plate_velocity_earth_morvel56[3*plate_label    ];

      /* non-dimensionalize */
      rhea_plate_velocity_earth_apply_dim_scal (
          rot_sp, 1 /* remove dim. */, opt);

      /* convert into Cartesian coordinates */
      RHEA_ASSERT (0.0 <= rot_sp[0]);
      RHEA_ASSERT (0.0 <= rot_sp[1] && rot_sp[1] <= 2.0*M_PI);
      RHEA_ASSERT (0.0 <= rot_sp[2] && rot_sp[2] <= M_PI);
      rot_axis[0] = rot_sp[0] * sin (rot_sp[2]) * cos (rot_sp[1]);
      rot_axis[1] = rot_sp[0] * sin (rot_sp[2]) * sin (rot_sp[1]);
      rot_axis[2] = rot_sp[0] * cos (rot_sp[2]);
    }
    break;

  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }
}

static void
rhea_plate_velocity_xsection_apply_dim_scal (double *velocity,
                                             const int remove_dim,
                                             rhea_plate_options_t *opt)
{
  /* check input */
  RHEA_ASSERT (isfinite (*velocity));

  if (!remove_dim) { /* multiply by dimensional scaling */
    *velocity *= opt->temp_options->thermal_diffusivity_m2_s /
                 opt->domain_options->radius_max_m;
  }
  else { /* otherwise remove dimensional scaling */
    *velocity /= opt->temp_options->thermal_diffusivity_m2_s /
                 opt->domain_options->radius_max_m;
  }
}

/******************************************************************************
 * Plate Data
 *****************************************************************************/

static int
rhea_plate_data_exists (rhea_plate_options_t *opt)
{
  /* exit if nothing to do */
  if (opt == NULL) {
    return 0;
  }

  /* return if data was allocated */
  switch (opt->domain_options->shape) {
  case RHEA_DOMAIN_CUBE:
  case RHEA_DOMAIN_SHELL:
  case RHEA_DOMAIN_CUBE_SPHERICAL:
    return (opt->vertices_x != NULL && opt->vertices_y != NULL &&
            opt->n_vertices != NULL && opt->angular_velocity != NULL);
  case RHEA_DOMAIN_BOX:
  case RHEA_DOMAIN_BOX_SPHERICAL:
    return (opt->xsection_boundary != NULL &&
            opt->xsection_tangential_velocity != NULL);
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }
}

int
rhea_plate_data_create (rhea_plate_options_t *opt, sc_MPI_Comm mpicomm)
{
  int                 n_plates = rhea_plate_get_n_plates (opt);
  int                 n_vertices_total, n_vertices_coarse_total;
  float              *vertices_all;
  int                 success_xsection = 0;
  int                 success_vert = 0;
  int                 success_vert_coarse = 0;

  /* check input */
  RHEA_ASSERT (opt != NULL);

  /* exit if nothing to do */
  if (rhea_plate_data_exists (opt)) {
    return 0;
  }

  /* start performance monitors */
  ymir_perf_counter_start (
      &rhea_plate_perfmon[RHEA_PLATE_PERFMON_DATA_CREATE]);

  /*
   * Interval Based (1D) Plates
   */

  /* fill plate boundary list for cross sectional domain */
  if (NULL != opt->xsection_boundary_lon_list &&
      NULL != opt->xsection_tangential_velocity_mm_yr_list) {
    double             *list;
    int                 n_entries, pid;

    /* get plate boundary values from string */
    list = NULL;
    n_entries = ymir_options_convert_string_to_double (
        opt->xsection_boundary_lon_list, &list);
    success_xsection += n_entries;
    opt->xsection_n_intervals = n_entries - 1;

    /* set plate boundaries; convert degrees to radians */
    opt->xsection_boundary = RHEA_ALLOC (float, 2 * opt->xsection_n_intervals);
    opt->xsection_boundary[0] = (float) list[0] * M_PI/180.0;
    for (pid = 1; pid < n_entries-1; pid++) {
      opt->xsection_boundary[2*pid-1] = (float) list[pid] * M_PI/180.0;
      opt->xsection_boundary[2*pid  ] = (float) list[pid] * M_PI/180.0;
    }
    opt->xsection_boundary[2*(n_entries-1) - 1] = (float) list[n_entries-1] *
                                                  M_PI/180.0;
    YMIR_FREE (list); /* was allocated in ymir */

    /* shrink plate sizes */
    if (isfinite (opt->xsection_shrink_factor) &&
        0.0 < opt->xsection_shrink_factor) {
      double              length, shift;

      for (pid = 0; pid < opt->xsection_n_intervals; pid++) {
        length = opt->xsection_boundary[2*pid+1] -
                 opt->xsection_boundary[2*pid];
        shift = 0.5 * length * (1.0 - opt->xsection_shrink_factor);
        RHEA_ASSERT (0.0 < length);
        RHEA_ASSERT (0.0 < shift);
        opt->xsection_boundary[2*pid  ] += shift;
        opt->xsection_boundary[2*pid+1] -= shift;
      }
    }

    /* get plate velocity values from string */
    list = NULL;
    n_entries = ymir_options_convert_string_to_double (
        opt->xsection_tangential_velocity_mm_yr_list, &list);
    success_xsection += n_entries;
    RHEA_CHECK_ABORT (opt->xsection_n_intervals == n_entries,
                      "Mismatch of given #plates and #velocities");

    /* set plate velocities; remove dimensional scaling */
    opt->xsection_tangential_velocity =
      RHEA_ALLOC (double, opt->xsection_n_intervals);
    for (pid = 0; pid < opt->xsection_n_intervals; pid++) {
      opt->xsection_tangential_velocity[pid] =
          list[pid] / (1000.0 * RHEA_TEMPERATURE_SECONDS_PER_YEAR);
      rhea_plate_velocity_xsection_apply_dim_scal (
          &(opt->xsection_tangential_velocity[pid]), 1 /* remove_dim */, opt);
    }
    YMIR_FREE (list); /* was allocated in ymir */

    /* print plate data */
#ifdef RHEA_ENABLE_DEBUG
    for (pid = 0; pid < opt->xsection_n_intervals; pid++) {
      RHEA_GLOBAL_VERBOSEF_FN_TAG (
          __func__,
          "xsection_plate_idx=%i, boundaries=(%g,%g), tangential_velocity=%g",
          pid, opt->xsection_boundary[2*pid], opt->xsection_boundary[2*pid+1],
          opt->xsection_tangential_velocity[pid]);
      RHEA_ASSERT (opt->xsection_boundary[2*pid] <
                   opt->xsection_boundary[2*pid+1]);
    }
#endif
  }
  else {
    opt->xsection_boundary = NULL;
    opt->xsection_tangential_velocity = NULL;
    opt->xsection_n_intervals = 0;
  }

  /*
   * Polygon Based (2D) Plates
   */

  /* correct #polygons if no vertices are given and exit */
  if (opt->vertices_file_path_txt == NULL) {
    opt->n_polygons = 0;

    /* stop performance monitors */
    ymir_perf_counter_stop_add (
        &rhea_plate_perfmon[RHEA_PLATE_PERFMON_DATA_CREATE]);

    return (success_xsection + success_vert + success_vert_coarse);
  }

  /* get parameters from options */
  n_vertices_total = opt->n_vertices_total;
  n_vertices_coarse_total = opt->n_vertices_coarse_total;

  /* read (fine) polygon vertices (if they do not exist) */
  if (opt->vertices_file_path_txt != NULL &&
      opt->vertices_x == NULL && opt->vertices_y == NULL &&
      opt->n_vertices == NULL) {
    /* read file */
    vertices_all = RHEA_ALLOC (float, 2*n_vertices_total);
    success_vert = rhea_io_mpi_read_broadcast_float (
        vertices_all, 2*n_vertices_total, NULL /* path bin */,
        opt->vertices_file_path_txt, mpicomm);
    RHEA_ASSERT (success_vert == 2*n_vertices_total);

    /* process vertices */
    if (success_vert) {
      opt->vertices_x = RHEA_ALLOC (float *, n_plates);
      opt->vertices_y = RHEA_ALLOC (float *, n_plates);
      opt->n_vertices = RHEA_ALLOC (size_t, n_plates);

      switch (opt->domain_options->shape) {
      case RHEA_DOMAIN_CUBE:
        rhea_plate_process_vertices_cube (
            opt->vertices_x, opt->vertices_y, opt->n_vertices, n_plates,
            vertices_all, n_vertices_total,
            opt->x_min, opt->x_max, opt->y_min, opt->y_max);
        break;
      case RHEA_DOMAIN_SHELL:
        rhea_plate_process_vertices_shell (
            opt->vertices_x, opt->vertices_y, opt->n_vertices, n_plates,
            vertices_all, n_vertices_total,
            opt->x_min, opt->x_max, opt->y_min, opt->y_max);
        break;
      default: /* unknown domain shape */
        RHEA_ABORT_NOT_REACHED ();
      }
    }

    /* destroy */
    RHEA_FREE (vertices_all);
  }
  else {
    success_vert = 1;
  }

  /* read coarse container vertices (if they do not exist) */
  if (opt->vertices_coarse_container_file_path_txt != NULL &&
      opt->vertices_coarse_container_x == NULL &&
      opt->vertices_coarse_container_y == NULL &&
      opt->n_vertices_coarse_container == NULL) {
    /* read file */
    vertices_all = RHEA_ALLOC (float, 2*n_vertices_coarse_total);
    success_vert_coarse = rhea_io_mpi_read_broadcast_float (
        vertices_all, 2*n_vertices_coarse_total, NULL /* path bin */,
        opt->vertices_coarse_container_file_path_txt, mpicomm);
    RHEA_ASSERT (success_vert == 2*n_vertices_coarse_total);

    /* process vertices */
    if (success_vert_coarse) {
      opt->vertices_coarse_container_x = RHEA_ALLOC (float *, n_plates);
      opt->vertices_coarse_container_y = RHEA_ALLOC (float *, n_plates);
      opt->n_vertices_coarse_container = RHEA_ALLOC (size_t, n_plates);

      switch (opt->domain_options->shape) {
      case RHEA_DOMAIN_CUBE:
        rhea_plate_process_vertices_cube (
            opt->vertices_coarse_container_x, opt->vertices_coarse_container_y,
            opt->n_vertices_coarse_container, n_plates,
            vertices_all, n_vertices_total,
            opt->x_min, opt->x_max, opt->y_min, opt->y_max);
        break;
      case RHEA_DOMAIN_SHELL:
        rhea_plate_process_vertices_shell (
            opt->vertices_coarse_container_x, opt->vertices_coarse_container_y,
            opt->n_vertices_coarse_container, n_plates,
            vertices_all, n_vertices_total,
            opt->x_min, opt->x_max, opt->y_min, opt->y_max);
        break;
      default: /* unknown domain shape */
        RHEA_ABORT_NOT_REACHED ();
      }
    }

    /* destroy */
    RHEA_FREE (vertices_all);
  }
  else {
    success_vert_coarse = 1;
  }

  /* set translations */
  if (opt->vertices_x != NULL && opt->vertices_y != NULL &&
      opt->n_vertices != NULL) { /* if vertices exist */
    switch (opt->domain_options->shape) {
    case RHEA_DOMAIN_CUBE:
      break;
    case RHEA_DOMAIN_SHELL:
      opt->translation_x = RHEA_ALLOC (float, n_plates);
      opt->translation_y = RHEA_ALLOC (float, n_plates);
      if (opt->vertices_coarse_container_x != NULL &&
          opt->vertices_coarse_container_y != NULL &&
          opt->n_vertices_coarse_container != NULL) { /* use coarse vertices */
        RHEA_ABORT_NOT_REACHED (); //TODO remove caorse containers
      }
      else { /* use vertices of (fine) polygon */
        rhea_plate_modify_vertices_shell (
            opt->translation_x, opt->translation_y,
            opt->vertices_x, opt->vertices_y, opt->n_vertices, n_plates);
      }
      break;
    default: /* unknown domain shape */
      RHEA_ABORT_NOT_REACHED ();
    }
  }

  /* create velocities (if they do not exist) */
  if (opt->angular_velocity == NULL) {
    int                 pid;
    double              rot_axis[3];

    opt->angular_velocity = RHEA_ALLOC (double, n_plates*3);
    for (pid = 0; pid < n_plates; pid++) { /* loop over all plates */
      rhea_plate_velocity_lookup_angular_velocity (rot_axis, pid, opt);
      opt->angular_velocity[3*pid    ] = rot_axis[0];
      opt->angular_velocity[3*pid + 1] = rot_axis[1];
      opt->angular_velocity[3*pid + 2] = rot_axis[2];
    }
  }

  /* stop performance monitors */
  ymir_perf_counter_stop_add (
      &rhea_plate_perfmon[RHEA_PLATE_PERFMON_DATA_CREATE]);

  return (success_xsection + success_vert + success_vert_coarse);
}

void
rhea_plate_data_clear (rhea_plate_options_t *opt)
{
  int                 pid;

  /* exit if nothing to do */
  if (!rhea_plate_data_exists (opt)) {
    return;
  }

  /*
   * Interval Based (1D) Plates
   */

  /* destroy plate boundary data for cross sectional domain */
  if (NULL != opt->xsection_boundary) {
    RHEA_FREE (opt->xsection_boundary);
    opt->xsection_boundary = NULL;
    opt->xsection_n_intervals = 0;
  }
  if (NULL != opt->xsection_tangential_velocity) {
    RHEA_FREE (opt->xsection_tangential_velocity);
    opt->xsection_tangential_velocity = NULL;
  }

  /*
   * Polygon Based (2D) Plates
   */

  /* destroy polygon vertices */
  if (opt->vertices_x != NULL && opt->vertices_y != NULL &&
      opt->n_vertices != NULL) {
    for (pid = 0; pid < rhea_plate_get_n_plates (opt); pid++) {
      RHEA_FREE (opt->vertices_x[pid]);
      RHEA_FREE (opt->vertices_y[pid]);
    }

    RHEA_FREE (opt->vertices_x);
    RHEA_FREE (opt->vertices_y);
    RHEA_FREE (opt->n_vertices);
    opt->vertices_x = NULL;
    opt->vertices_y = NULL;
    opt->n_vertices = NULL;
  }

  /* destroy container vertices */
  if (opt->vertices_coarse_container_x != NULL &&
      opt->vertices_coarse_container_y != NULL &&
      opt->n_vertices_coarse_container != NULL) {
    for (pid = 0; pid < rhea_plate_get_n_plates (opt); pid++) {
      RHEA_FREE (opt->vertices_coarse_container_x[pid]);
      RHEA_FREE (opt->vertices_coarse_container_y[pid]);
    }

    RHEA_FREE (opt->vertices_coarse_container_x);
    RHEA_FREE (opt->vertices_coarse_container_y);
    RHEA_FREE (opt->n_vertices_coarse_container);
    opt->vertices_coarse_container_x = NULL;
    opt->vertices_coarse_container_y = NULL;
    opt->n_vertices_coarse_container = NULL;
  }

  /* destroy translation */
  if (opt->translation_x != NULL) {
    RHEA_FREE (opt->translation_x);
    opt->translation_x = NULL;
  }
  if (opt->translation_y != NULL) {
    RHEA_FREE (opt->translation_y);
    opt->translation_y = NULL;
  }

  /* destroy velocities */
  if (opt->angular_velocity != NULL) {
    RHEA_FREE (opt->angular_velocity);
  }
}

int
rhea_plate_get_n_plates (rhea_plate_options_t *opt)
{
  /* exit if nothing to do */
  if (opt == NULL) {
    return -1;
  }

  switch (opt->domain_options->shape) {
  case RHEA_DOMAIN_CUBE:
  case RHEA_DOMAIN_SHELL:
  case RHEA_DOMAIN_CUBE_SPHERICAL:
    /* return number of polygons */
    return opt->n_polygons;
  case RHEA_DOMAIN_BOX:
  case RHEA_DOMAIN_BOX_SPHERICAL:
    /* return number of intercals */
    return opt->xsection_n_intervals;
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }
}

/******************************************************************************
 * Plate Retrieval
 *****************************************************************************/

/**
 * Checks whether the given (x,y) coordinates are inside a polygon.
 */
static int
rhea_plate_point_inside_polygon_test (
                            const float test_x,
                            const float test_y,
                            const float *_sc_restrict vert_x,
                            const float *_sc_restrict vert_y,
                            const size_t n_vert,
                            const float *_sc_restrict vert_coarse_container_x,
                            const float *_sc_restrict vert_coarse_container_y,
                            const size_t n_vert_coarse_container)
{
  int                 is_inside = 1;

  /* check input */
  RHEA_ASSERT (vert_x != NULL);
  RHEA_ASSERT (vert_y != NULL);
  RHEA_ASSERT (0 < n_vert);
  RHEA_ASSERT (0 >= n_vert_coarse_container || vert_coarse_container_x != NULL);
  RHEA_ASSERT (0 >= n_vert_coarse_container || vert_coarse_container_y != NULL);

  /* start performance monitors */
  ymir_perf_counter_start (
      &rhea_plate_perfmon[RHEA_PLATE_PERFMON_POINT_INSIDE_POLYGON]);

  /* check coarse container of plate */
  if (0 < n_vert_coarse_container) {
    is_inside = rhea_point_in_polygon_is_inside (
        test_x, test_y, vert_coarse_container_x, vert_coarse_container_y,
        n_vert_coarse_container);
  }

  /* check fine-resolution boundary of plate */
  if (is_inside) {
    is_inside = rhea_point_in_polygon_is_inside (test_x, test_y,
                                                 vert_x, vert_y, n_vert);
  }
  else {
    is_inside = 0;
  }

  /* stop performance monitors */
  ymir_perf_counter_stop_add (
      &rhea_plate_perfmon[RHEA_PLATE_PERFMON_POINT_INSIDE_POLYGON]);

  /* return if test coordinate is inside polygon */
  return is_inside;
}

static int
rhea_plate_is_inside_polygon (const double x, const double y, const double z,
                              const int plate_label, rhea_plate_options_t *opt)
{
  const float        *vx = opt->vertices_x[plate_label];
  const float        *vy = opt->vertices_y[plate_label];
  const size_t        nv = opt->n_vertices[plate_label];
  float              *vcx = NULL, *vcy = NULL;
  size_t              nvc = 0;
  double              tmp_x, tmp_y;
  float               test_x, test_y;

  /* check input */
  RHEA_ASSERT (opt->vertices_x != NULL);
  RHEA_ASSERT (opt->vertices_y != NULL);
  RHEA_ASSERT (opt->n_vertices != NULL);
  RHEA_ASSERT (plate_label < opt->n_polygons);

  /* calculate test coordinates */
  rhea_domain_extract_lateral (&tmp_x, &tmp_y, x, y, z,
                               RHEA_DOMAIN_COORDINATE_SPHERICAL_GEO_DIM,
                               opt->domain_options);
  test_x = (float) tmp_x;
  test_y = (float) tmp_y;

  /* apply translation into east-west direction */
  if (opt->translation_x != NULL) {
    const float         translation_x = opt->translation_x[plate_label];

    if (0.0 <= (test_x + translation_x)) {
      test_x += translation_x;
    }
    else {
      test_x += translation_x + 360.0;
    }
  }
  RHEA_ASSERT (0.0 <= test_x && test_x <= 360.0);

  /* apply translation into south-north direction */
  if (opt->translation_y != NULL) {
    const float         translation_y = opt->translation_y[plate_label];

    test_y += translation_y;
  }
  RHEA_ASSERT (0.0 <= test_y);

  /* get vertices of coarse container (if they exist) */
  if (opt->vertices_coarse_container_x != NULL &&
      opt->vertices_coarse_container_y != NULL &&
      opt->n_vertices_coarse_container != NULL) {
    vcx = opt->vertices_coarse_container_x[plate_label];
    vcy = opt->vertices_coarse_container_y[plate_label];
    nvc = opt->n_vertices_coarse_container[plate_label];
  }

  /* return if inside polygon */
  return rhea_plate_point_inside_polygon_test (test_x, test_y, vx, vy, nv,
                                               vcx, vcy, nvc);
}

static int
rhea_plate_is_inside_interval (const double x, const double y, const double z,
                               const int plate_label, rhea_plate_options_t *opt)
{
  const float         boundary_begin = opt->xsection_boundary[2*plate_label  ];
  const float         boundary_end   = opt->xsection_boundary[2*plate_label+1];
  double              tmp_x, tmp_y;

  /* check input */
  RHEA_ASSERT (opt->xsection_boundary != NULL);
  RHEA_ASSERT (plate_label < opt->xsection_n_intervals);

  /* calculate test coordinates */
  rhea_domain_extract_lateral (&tmp_x, &tmp_y, x, y, z,
                               RHEA_DOMAIN_COORDINATE_SPHERICAL_GEO,
                               opt->domain_options);

  /* shift by min longitude */
  switch (opt->domain_options->shape) {
  case RHEA_DOMAIN_BOX:
    break;
  case RHEA_DOMAIN_BOX_SPHERICAL:
    tmp_x -= opt->domain_options->lon_min;
    break;
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* return if inside interval */
  return (boundary_begin <= tmp_x && tmp_x <= boundary_end);
}

int
rhea_plate_is_inside (const double x, const double y, const double z,
                      const int plate_label, rhea_plate_options_t *opt)
{
  switch (opt->domain_options->shape) {
  case RHEA_DOMAIN_CUBE:
  case RHEA_DOMAIN_SHELL:
  case RHEA_DOMAIN_CUBE_SPHERICAL:
    return rhea_plate_is_inside_polygon (x, y, z, plate_label, opt);
  case RHEA_DOMAIN_BOX:
  case RHEA_DOMAIN_BOX_SPHERICAL:
    return rhea_plate_is_inside_interval (x, y, z, plate_label, opt);
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }
}

/**
 * Gets the plate label at the given (x,y,z) coordinates.
 */
static int
rhea_plate_get_label (const double x, const double y, const double z,
                      rhea_plate_options_t *opt)
{
  const int           n_plates = rhea_plate_get_n_plates (opt);
  int                 pid;

  /* find and return plate label */
  for (pid = 0; pid < n_plates; pid++) { /* loop over all plates */
    if (rhea_plate_is_inside (x, y, z, pid, opt)) {
      return pid;
    }
  }

  /* return that test coordinates are outside of any plate */
  return RHEA_PLATE_NONE;
}

/* data for plate retrieval */
typedef struct rhea_plate_node_fn_data
{
  rhea_plate_options_t *opt;
  int                 plate_label;
  int                 n_fields;
  ymir_vec_t         *filter;
}
rhea_plate_node_fn_data_t;

/**
 * Sets plate label at a single node.
 */
static void
rhea_plate_set_label_node_fn (double *v, double x, double y, double z,
                              ymir_locidx_t nodeid, void *data)
{
  rhea_plate_node_fn_data_t *d = data;
  const int           n_fields = d->n_fields;
  int                 fieldid;
  int                 plate_label;

  plate_label = rhea_plate_get_label (x, y, z, d->opt);
  for (fieldid = 0; fieldid < n_fields; fieldid++) {
    v[fieldid] = (double) plate_label;
  }
}

static void
rhea_plate_set_label_face_node_fn (double *v, double x, double y, double z,
                                   double nx, double ny, double nz,
                                   ymir_topidx_t face, ymir_locidx_t nodeid,
                                   void *data)
{
  rhea_plate_node_fn_data_t *d = data;
  const int           n_fields = d->n_fields;
  int                 fieldid;
  int                 plate_label;

  plate_label = rhea_plate_get_label (x, y, z, d->opt);
  for (fieldid = 0; fieldid < n_fields; fieldid++) {
    v[fieldid] = (double) plate_label;
  }
}

void
rhea_plate_set_label_vec (ymir_vec_t *vec, rhea_plate_options_t *opt)
{
  rhea_plate_node_fn_data_t  data;

  /* exit if nothing to do */
  if (!rhea_plate_data_exists (opt)) {
    ymir_vec_set_value (vec, (double) RHEA_PLATE_NONE);
    return;
  }

  /* start performance monitors */
  ymir_perf_counter_start (
      &rhea_plate_perfmon[RHEA_PLATE_PERFMON_SET_LABEL]);

  /* set data parameters */
  data.opt = opt;
  data.plate_label = RHEA_PLATE_NONE;

  if (!ymir_vec_is_face_vec (vec)) { /* if volume vector */
    /* set labels at continuous GLL nodes */
    if (ymir_vec_has_cvec (vec)) {
      data.n_fields = vec->ncfields;
      ymir_cvec_set_function (vec, rhea_plate_set_label_node_fn, &data);
    }

    /* set labels at discontinuous GLL/Gauss nodes */
    if (ymir_vec_has_dvec (vec)) {
      data.n_fields = vec->ndfields;
      ymir_dvec_set_function (vec, rhea_plate_set_label_node_fn, &data);
    }

    /* note: processing of evec's is not implemented */
    if (ymir_vec_has_evec (vec)) {
      RHEA_ABORT_NOT_REACHED ();
    }
  }
  else { /* if face vector */
    /* set labels at continuous GLL nodes */
    if (ymir_vec_has_cvec (vec)) {
      data.n_fields = vec->ncfields;
      ymir_face_cvec_set_function (vec, rhea_plate_set_label_face_node_fn,
                                   &data);
    }

    /* set labels at discontinuous GLL/Gauss nodes */
    if (ymir_vec_has_dvec (vec)) {
      data.n_fields = vec->ndfields;
      ymir_face_dvec_set_function (vec, rhea_plate_set_label_face_node_fn,
                                   &data);
    }

    /* note: processing of evec's is not implemented */
    if (ymir_vec_has_evec (vec)) {
      RHEA_ABORT_NOT_REACHED ();
    }
  }

  /* stop performance monitors */
  ymir_perf_counter_stop_add (
      &rhea_plate_perfmon[RHEA_PLATE_PERFMON_SET_LABEL]);
}

void
rhea_plate_set_weight_vec (ymir_vec_t *vec,
                           rhea_plate_area_to_weight_fn_t area_to_weight_fn,
                           rhea_plate_options_t *opt)
{
  ymir_vec_t         *plate_filter, *plate_filter_mass;
  double              total_area, plate_area, weight;
  const int           n_plates = rhea_plate_get_n_plates (opt);
  int                 pid;

  /* exit if nothing to do */
  if (!rhea_plate_data_exists (opt)) {
    return;
  }

  /* check input */
  RHEA_ASSERT (ymir_vec_is_face_vec (vec) &&
               vec->meshnum == RHEA_DOMAIN_BOUNDARY_FACE_TOP);

  /* initialize */
  if ( (ymir_vec_is_cvec (vec) && vec->ncfields == 1) ||
       (ymir_vec_is_dvec (vec) && vec->ndfields == 1) ) {
    plate_filter      = ymir_vec_template (vec);
    plate_filter_mass = ymir_vec_template (vec);
  }
  else {
    RHEA_ABORT_NOT_REACHED ();
  }
  ymir_vec_set_zero (vec);

  /* compute total area */
  ymir_vec_set_value (plate_filter, 1.0);
  ymir_mass_apply (plate_filter, plate_filter_mass);
  total_area = ymir_vec_innerprod (plate_filter_mass, plate_filter);

  /* add each plate with its weight to the combined filter */
  for (pid = 0; pid < n_plates; pid++) {
    /* filter plate */
    ymir_vec_set_value (plate_filter, 1.0);
    rhea_plate_apply_filter_vec (plate_filter, pid, opt);

    /* compute plate area */
    ymir_mass_apply (plate_filter, plate_filter_mass);
    plate_area = ymir_vec_innerprod (plate_filter_mass, plate_filter);

    /* add plate's weight */
    if (NULL != area_to_weight_fn) {
      weight = area_to_weight_fn (plate_area, total_area);
    }
    else {
      weight = total_area/plate_area;
    }
    ymir_vec_add (weight, plate_filter, vec);
  }

  /* destroy */
  ymir_vec_destroy (plate_filter);
  ymir_vec_destroy (plate_filter_mass);
}

/**
 * Filters values at a single node.
 */
static void
rhea_plate_apply_filter_node_fn (double *v, double x, double y, double z,
                                 ymir_locidx_t nodeid, void *data)
{
  rhea_plate_node_fn_data_t *d = data;

  if (!rhea_plate_is_inside (x, y, z, d->plate_label, d->opt)) { /* if out */
    const int           n_fields = d->n_fields;
    int                 fieldid;

    /* set values to zero outside of plate */
    for (fieldid = 0; fieldid < n_fields; fieldid++) {
      v[fieldid] = 0.0;
    }
  }
}

static void
rhea_plate_apply_filter_face_node_fn (double *v, double x, double y, double z,
                                      double nx, double ny, double nz,
                                      ymir_topidx_t face, ymir_locidx_t nodeid,
                                      void *data)
{
  rhea_plate_node_fn_data_t *d = data;

  if (!rhea_plate_is_inside (x, y, z, d->plate_label, d->opt)) { /* if out */
    const int           n_fields = d->n_fields;
    int                 fieldid;

    /* set values to zero outside of plate */
    for (fieldid = 0; fieldid < n_fields; fieldid++) {
      v[fieldid] = 0.0;
    }
  }
}

void
rhea_plate_apply_filter_vec (ymir_vec_t *vec, const int plate_label,
                             rhea_plate_options_t *opt)
{
  rhea_plate_node_fn_data_t  data;

  /* exit if nothing to do */
  if (!rhea_plate_data_exists (opt)) {
    return;
  }

  /* start performance monitors */
  ymir_perf_counter_start (
      &rhea_plate_perfmon[RHEA_PLATE_PERFMON_APPLY_FILTER]);

  /* set data parameters */
  data.opt = opt;
  data.plate_label = plate_label;

  if (!ymir_vec_is_face_vec (vec)) { /* if volume vector */
    /* filter vector at continuous GLL nodes */
    if (ymir_vec_has_cvec (vec)) {
      data.n_fields = vec->ncfields;
      ymir_cvec_set_function (vec, rhea_plate_apply_filter_node_fn, &data);
    }

    /* filter vector at discontinuous GLL/Gauss nodes */
    if (ymir_vec_has_dvec (vec)) {
      data.n_fields = vec->ndfields;
      ymir_dvec_set_function (vec, rhea_plate_apply_filter_node_fn, &data);
    }

    /* note: processing of evec's is not implemented */
    if (ymir_vec_has_evec (vec)) {
      RHEA_ABORT_NOT_REACHED ();
    }
  }
  else { /* if face vector */
    /* filter vector at continuous GLL nodes */
    if (ymir_vec_has_cvec (vec)) {
      data.n_fields = vec->ncfields;
      ymir_face_cvec_set_function (vec, rhea_plate_apply_filter_face_node_fn,
                                   &data);
    }

    /* filter vector at discontinuous GLL/Gauss nodes */
    if (ymir_vec_has_dvec (vec)) {
      data.n_fields = vec->ndfields;
      ymir_face_dvec_set_function (vec, rhea_plate_apply_filter_face_node_fn,
                                   &data);
    }

    /* note: processing of evec's is not implemented */
    if (ymir_vec_has_evec (vec)) {
      RHEA_ABORT_NOT_REACHED ();
    }
  }

  /* stop performance monitors */
  ymir_perf_counter_stop_add (
      &rhea_plate_perfmon[RHEA_PLATE_PERFMON_APPLY_FILTER]);
}

/**
 * Filters values at a single node.
 */
static void
rhea_plate_apply_filter_all_node_fn (double *v, double x, double y, double z,
                                     ymir_locidx_t nodeid, void *data)
{
  rhea_plate_node_fn_data_t *d = data;
  const double        filter_val = *ymir_cvec_index (d->filter, nodeid, 0);

  if (filter_val < 0.0) {
    const int           n_fields = d->n_fields;
    int                 fieldid;

    /* set values to zero outside of plate */
    for (fieldid = 0; fieldid < n_fields; fieldid++) {
      v[fieldid] = 0.0;
    }
  }
}

static void
rhea_plate_apply_filter_all_face_node_fn (double *v, double x, double y,
                                          double z, double nx, double ny,
                                          double nz, ymir_topidx_t face,
                                          ymir_locidx_t nodeid, void *data)
{
  rhea_plate_node_fn_data_t *d = data;
  const double        filter_val = *ymir_cvec_index (d->filter, nodeid, 0);

  if (filter_val < 0.0) {
    const int           n_fields = d->n_fields;
    int                 fieldid;

    /* set values to zero outside of plate */
    for (fieldid = 0; fieldid < n_fields; fieldid++) {
      v[fieldid] = 0.0;
    }
  }
}

void
rhea_plate_apply_filter_all_vec (ymir_vec_t *vec, rhea_plate_options_t *opt)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vec);
  rhea_plate_node_fn_data_t  data;

  /* exit if nothing to do */
  if (!rhea_plate_data_exists (opt)) {
    return;
  }

  if (!ymir_vec_is_face_vec (vec)) { /* if volume vector */
    /* filter vector at continuous GLL nodes */
    if (ymir_vec_has_cvec (vec)) {
      data.n_fields = vec->ncfields;
      data.filter = ymir_cvec_new (ymir_mesh, 1);
      rhea_plate_set_label_vec (data.filter, opt);
      ymir_cvec_set_function (vec, rhea_plate_apply_filter_all_node_fn, &data);
      ymir_vec_destroy (data.filter);
    }

    /* note: processing of dvec's is not implemented */
    if (ymir_vec_has_dvec (vec)) {
      RHEA_ABORT_NOT_REACHED ();
    }

    /* note: processing of evec's is not implemented */
    if (ymir_vec_has_evec (vec)) {
      RHEA_ABORT_NOT_REACHED ();
    }
  }
  else { /* if face vector */
    /* filter vector at continuous GLL nodes */
    if (ymir_vec_has_cvec (vec)) {
      data.n_fields = vec->ncfields;
      data.filter = ymir_face_cvec_new (ymir_mesh, vec->meshnum, 1);
      rhea_plate_set_label_vec (data.filter, opt);
      ymir_face_cvec_set_function (
          vec, rhea_plate_apply_filter_all_face_node_fn, &data);
      ymir_vec_destroy (data.filter);
    }

    /* note: processing of dvec's is not implemented */
    if (ymir_vec_has_dvec (vec)) {
      RHEA_ABORT_NOT_REACHED ();
    }

    /* note: processing of evec's is not implemented */
    if (ymir_vec_has_evec (vec)) {
      RHEA_ABORT_NOT_REACHED ();
    }
  }
}

/******************************************************************************
 * Plate Velocity Fields
 *****************************************************************************/

void
rhea_plate_velocity_generate_from_rotation (ymir_vec_t *vel,
                                            const double rot_axis[3],
                                            rhea_plate_options_t *opt)
{
  ymir_velocity_vec_generate_rotation (vel, opt->domain_options->center,
                                       rot_axis);
}

static void
rhea_plate_velocity_generate_polygon (ymir_vec_t *vel, const int plate_label,
                                      rhea_plate_options_t *opt)
{
  double              rot_axis[3];

  /* check input */
  RHEA_ASSERT (RHEA_PLATE_NONE <= plate_label &&
               plate_label < rhea_plate_get_n_plates (opt));

  /* exit if nothing to do */
  if (!rhea_plate_data_exists (opt)) {
    ymir_vec_set_value (vel, NAN);
    return;
  }
  else if (RHEA_PLATE_NONE == plate_label) {
    ymir_vec_set_value (vel, 0.0);
    return;
  }

  /* get angular velocity */
  rot_axis[0] = opt->angular_velocity[3*plate_label    ];
  rot_axis[1] = opt->angular_velocity[3*plate_label + 1];
  rot_axis[2] = opt->angular_velocity[3*plate_label + 2];

  /* generate velocity field and apply filter for the plate's interior */
  rhea_plate_velocity_generate_from_rotation (vel, rot_axis, opt);
  rhea_plate_apply_filter_vec (vel, plate_label, opt);
}

static void
rhea_plate_set_unit_tangential_component_fn (double *vec,
                                             double X, double Y, double Z,
                                             double nx, double ny, double nz,
                                             ymir_topidx_t face,
                                             ymir_locidx_t node_id,
                                             void *data)
{
  rhea_plate_options_t *opt = data;
  double              magn;

  RHEA_ASSERT (opt->domain_options->box_subdivision_y == 1);

  /* rotate surface normal vector by 90 degrees clockwise in (x,z)-plane:
   *          [ 0 0 +1]   [nx]
   *   tang = [ 0 1  0] * [ny]
   *          [-1 0  0]   [nz]
   */
  vec[0] = +nz;
  vec[1] = 0.0; /* neutralize y-component */
  vec[2] = -nx;

  /* normalize */
  magn = sqrt (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
  vec[0] *= 1.0/magn;
  vec[1] *= 1.0/magn;
  vec[2] *= 1.0/magn;
}

static void
rhea_plate_velocity_generate_interval (ymir_vec_t *vel, const int plate_label,
                                       rhea_plate_options_t *opt)
{
  /* check input */
  RHEA_ASSERT (RHEA_PLATE_NONE <= plate_label &&
               plate_label < rhea_plate_get_n_plates (opt));
  RHEA_ASSERT (rhea_velocity_surface_check_vec_type (vel));

  /* exit if nothing to do */
  if (!rhea_plate_data_exists (opt)) {
    ymir_vec_set_value (vel, NAN);
    return;
  }
  else if (RHEA_PLATE_NONE == plate_label) {
    ymir_vec_set_value (vel, 0.0);
    return;
  }

  /* set unit tangential component at the surface (clockwise) */
  ymir_face_cvec_set_function (
      vel, rhea_plate_set_unit_tangential_component_fn, opt);

  /* scale by plate velocity; apply filter for the plate's interior */
  ymir_vec_scale (opt->xsection_tangential_velocity[plate_label], vel);
  rhea_plate_apply_filter_vec (vel, plate_label, opt);
}

void
rhea_plate_velocity_generate (ymir_vec_t *vel, const int plate_label,
                              rhea_plate_options_t *opt)
{
  switch (opt->domain_options->shape) {
  case RHEA_DOMAIN_CUBE:
  case RHEA_DOMAIN_SHELL:
  case RHEA_DOMAIN_CUBE_SPHERICAL:
    return rhea_plate_velocity_generate_polygon (vel, plate_label, opt);
  case RHEA_DOMAIN_BOX:
  case RHEA_DOMAIN_BOX_SPHERICAL:
    return rhea_plate_velocity_generate_interval (vel, plate_label, opt);
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }
}

void
rhea_plate_velocity_generate_all (ymir_vec_t *vel, rhea_plate_options_t *opt)
{
  ymir_vec_t         *plate_vel = ymir_vec_template (vel);
  int                 pid;
#ifdef RHEA_ENABLE_DEBUG
  double              mean_vel_magn;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);
#endif

  /* check input */
  RHEA_ASSERT (opt != NULL);

  /* initialize output */
  ymir_vec_set_value (vel, 0.0);

  /* get velocities of all plates */
  for (pid = 0; pid < rhea_plate_get_n_plates (opt); pid++) {
    rhea_plate_velocity_generate (plate_vel, pid, opt);
    ymir_vec_add (1.0, plate_vel, vel);
#ifdef RHEA_ENABLE_DEBUG
    rhea_velocity_convert_to_dimensional_mm_yr (
        plate_vel, opt->domain_options, opt->temp_options);
    mean_vel_magn = rhea_plate_velocity_get_mean_magnitude (
        plate_vel, pid, 0 /* !projection */, opt);
    RHEA_GLOBAL_VERBOSEF_FN_TAG (
        __func__, "plate_idx=%i, mean_velocity_magn=\"%g mm/yr\"",
        pid, mean_vel_magn);
#endif
  }

  /* destroy */
  ymir_vec_destroy (plate_vel);

#ifdef RHEA_ENABLE_DEBUG
  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
#endif
}

void
rhea_plate_velocity_evaluate_rotation (double rot_axis[3],
                                       ymir_vec_t *vel,
                                       const int plate_label,
                                       const int project_out_mean_rot,
                                       rhea_plate_options_t *opt)
{
  const double       *center = opt->domain_options->center;
  double              moment_of_inertia[3];
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vel);
  ymir_vec_t         *vel_surf;

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (vel) ||
               rhea_velocity_surface_check_vec_type (vel));
  RHEA_ASSERT (opt != NULL);
  RHEA_ASSERT (RHEA_PLATE_NONE == plate_label || rhea_plate_data_exists (opt));

  /* create surface velocity */
  if (!ymir_vec_is_face_vec (vel)) { /* if volume vector */
    vel_surf = rhea_velocity_surface_new_from_vol (vel);
  }
  else { /* if face vector */
    vel_surf = rhea_velocity_surface_new (ymir_mesh);
    ymir_vec_copy (vel, vel_surf);
  }

  /* project out mean rotation (no net rotation) */
  if (project_out_mean_rot) {
    rhea_plate_velocity_project_out_mean_rotation (vel_surf, opt);
  }

  /* filter plate */
  if (RHEA_PLATE_NONE != plate_label) {
    const ymir_topidx_t face_surf = RHEA_DOMAIN_BOUNDARY_FACE_TOP;
    ymir_vec_t         *density;

    /* filter plate velocity */
    rhea_plate_apply_filter_vec (vel_surf, plate_label, opt);

    /* calculate the moment of inertia for the plate */
    density = ymir_face_dvec_new (ymir_mesh, face_surf, 1, YMIR_GAUSS_NODE);
    ymir_vec_set_value (density, 1.0);
    rhea_plate_apply_filter_vec (density, plate_label, opt);
    ymir_velocity_vec_compute_moment_of_inertia (moment_of_inertia, center,
                                                 density, NULL /* unused */,
                                                 -1 /* unused */);
    ymir_vec_destroy (density);
  }
  else {
    moment_of_inertia[0] = opt->domain_options->moment_of_inertia_surface[0];
    moment_of_inertia[1] = opt->domain_options->moment_of_inertia_surface[1];
    moment_of_inertia[2] = opt->domain_options->moment_of_inertia_surface[2];
  }

  /* compute the mean rotation */
  ymir_velocity_vec_compute_mean_rotation (
      rot_axis, vel_surf, center, moment_of_inertia, 0 /* !residual_space */);

  /* destroy */
  rhea_velocity_surface_destroy (vel_surf);

  /* check output */
  RHEA_ASSERT (isfinite (rot_axis[0]));
  RHEA_ASSERT (isfinite (rot_axis[1]));
  RHEA_ASSERT (isfinite (rot_axis[2]));
}

static double
rhea_plate_velocity_compute_mean_magnitude_internal (
                                                  ymir_vec_t *vel_surf,
                                                  const int plate_label,
                                                  ymir_vec_t *magn_surf,
                                                  ymir_vec_t *filter_surf,
                                                  ymir_vec_t *filter_surf_mass,
                                                  rhea_plate_options_t *opt)
{
  double              filter_area;

  /* check input */
  RHEA_ASSERT (rhea_velocity_surface_check_vec_type (vel_surf));
  RHEA_ASSERT (ymir_vec_is_cvec (magn_surf));
  RHEA_ASSERT (ymir_vec_is_cvec (filter_surf));
  RHEA_ASSERT (ymir_vec_is_cvec (filter_surf_mass));
  RHEA_ASSERT (ymir_vec_is_face_vec (magn_surf));
  RHEA_ASSERT (ymir_vec_is_face_vec (filter_surf));
  RHEA_ASSERT (ymir_vec_is_face_vec (filter_surf_mass));
  RHEA_ASSERT (RHEA_DOMAIN_BOUNDARY_FACE_TOP == magn_surf->meshnum);
  RHEA_ASSERT (RHEA_DOMAIN_BOUNDARY_FACE_TOP == filter_surf->meshnum);
  RHEA_ASSERT (RHEA_DOMAIN_BOUNDARY_FACE_TOP == filter_surf_mass->meshnum);
  RHEA_ASSERT (1 == magn_surf->ncfields);
  RHEA_ASSERT (1 == filter_surf->ncfields);
  RHEA_ASSERT (1 == filter_surf_mass->ncfields);
  RHEA_ASSERT (opt != NULL);

  /* generate filter for plate */
  ymir_vec_set_value (filter_surf, 1.0);
  rhea_plate_apply_filter_vec (filter_surf, plate_label, opt);

  /* calculate area */
  ymir_mass_apply (filter_surf, filter_surf_mass);
  filter_area = ymir_vec_innerprod (filter_surf_mass, filter_surf);

  /* compute pointwise magnitude */
  ymir_cvec_innerprod_pointwise (magn_surf, vel_surf, vel_surf);
  ymir_cvec_sqrt (magn_surf, magn_surf);

  /* compute mean velocity */
  ymir_vec_multiply_in (filter_surf, magn_surf);
  return ymir_vec_innerprod (filter_surf_mass, magn_surf) / filter_area;
}

double
rhea_plate_velocity_get_mean_magnitude (ymir_vec_t *vel,
                                        const int plate_label,
                                        const int project_out_mean_rot,
                                        rhea_plate_options_t *opt)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vel);
  const ymir_topidx_t face_surf = RHEA_DOMAIN_BOUNDARY_FACE_TOP;
  ymir_vec_t         *vel_surf;
  ymir_vec_t         *magn_surf, *filter_surf, *filter_surf_mass;
  double              mean_vel_magn;

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (vel) ||
               rhea_velocity_surface_check_vec_type (vel));

  /* exit if nothing to do */
  if (!rhea_plate_data_exists (opt)) {
    return NAN;
  }

  /* create work variables */
  magn_surf = ymir_face_cvec_new (ymir_mesh, face_surf, 1);
  filter_surf = ymir_face_cvec_new (ymir_mesh, face_surf, 1);
  filter_surf_mass = ymir_face_cvec_new (ymir_mesh, face_surf, 1);

  /* create surface velocity */
  if (!ymir_vec_is_face_vec (vel)) { /* if volume vector */
    vel_surf = rhea_velocity_surface_new_from_vol (vel);
  }
  else { /* if face vector */
    vel_surf = rhea_velocity_surface_new (ymir_mesh);
    ymir_vec_copy (vel, vel_surf);
  }

  /* project out mean rotation (no net rotation) */
  if (project_out_mean_rot) {
    rhea_plate_velocity_project_out_mean_rotation (vel_surf, opt);
  }

  /* compute mean velocity */
  mean_vel_magn = rhea_plate_velocity_compute_mean_magnitude_internal (
      vel_surf, plate_label, magn_surf, filter_surf, filter_surf_mass, opt);

  /* destroy */
  rhea_velocity_surface_destroy (vel_surf);
  ymir_vec_destroy (magn_surf);
  ymir_vec_destroy (filter_surf);
  ymir_vec_destroy (filter_surf_mass);

  /* return mean velocity */
  return mean_vel_magn;
}

void
rhea_plate_velocity_get_mean_magnitude_all (double *mean_vel_magn,
                                            ymir_vec_t *vel,
                                            const int *plate_label,
                                            const int project_out_mean_rot,
                                            rhea_plate_options_t *opt)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vel);
  const ymir_topidx_t face_surf = RHEA_DOMAIN_BOUNDARY_FACE_TOP;
  ymir_vec_t         *vel_surf;
  ymir_vec_t         *magn_surf, *filter_surf, *filter_surf_mass;
  const int           n_plates = rhea_plate_get_n_plates (opt);
  int                 pid;

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (vel) ||
               rhea_velocity_surface_check_vec_type (vel));

  /* exit if nothing to do */
  if (!rhea_plate_data_exists (opt)) {
    for (pid = 0; pid < n_plates; pid++) {
      mean_vel_magn[pid] = NAN;
    }
    return;
  }

  /* create work variables */
  magn_surf = ymir_face_cvec_new (ymir_mesh, face_surf, 1);
  filter_surf = ymir_face_cvec_new (ymir_mesh, face_surf, 1);
  filter_surf_mass = ymir_face_cvec_new (ymir_mesh, face_surf, 1);

  /* create surface velocity */
  if (!ymir_vec_is_face_vec (vel)) { /* if volume vector */
    vel_surf = rhea_velocity_surface_new_from_vol (vel);
  }
  else { /* if face vector */
    vel_surf = rhea_velocity_surface_new (ymir_mesh);
    ymir_vec_copy (vel, vel_surf);
  }

  /* project out mean rotation (no net rotation) */
  if (project_out_mean_rot) {
    rhea_plate_velocity_project_out_mean_rotation (vel_surf, opt);
  }

  /* compute mean velocities */
  for (pid = 0; pid < n_plates; pid++) {
    if (plate_label != NULL && plate_label[pid] <= 0) { /* if skip this plate */
      mean_vel_magn[pid] = NAN;
      continue;
    }
    mean_vel_magn[pid] = rhea_plate_velocity_compute_mean_magnitude_internal (
        vel_surf, pid, magn_surf, filter_surf, filter_surf_mass, opt);
  }

  /* destroy */
  rhea_velocity_surface_destroy (vel_surf);
  ymir_vec_destroy (magn_surf);
  ymir_vec_destroy (filter_surf);
  ymir_vec_destroy (filter_surf_mass);
}

void
rhea_plate_velocity_project_out_mean_rotation (ymir_vec_t *vel,
                                               rhea_plate_options_t *opt)
{
  double              rot_center[3], moment_of_inertia[3];

  /* exit if nothing to do */
  if (!rhea_plate_data_exists (opt)) {
    return;
  }

  /* get center from options */
  rot_center[0] = opt->domain_options->center[0];
  rot_center[1] = opt->domain_options->center[1];
  rot_center[2] = opt->domain_options->center[2];

  /* get moment of inertia from options */
  if (!ymir_vec_is_face_vec (vel)) { /* if volume vector */
    moment_of_inertia[0] = opt->domain_options->moment_of_inertia[0];
    moment_of_inertia[1] = opt->domain_options->moment_of_inertia[1];
    moment_of_inertia[2] = opt->domain_options->moment_of_inertia[2];
  }
  else { /* if face vector */
    moment_of_inertia[0] = opt->domain_options->moment_of_inertia_surface[0];
    moment_of_inertia[1] = opt->domain_options->moment_of_inertia_surface[1];
    moment_of_inertia[2] = opt->domain_options->moment_of_inertia_surface[2];
  }

  /* project out mean rotation */
  ymir_velocity_vec_project_out_mean_rotation (
      vel, rot_center, moment_of_inertia, 0 /* !residual_space */);
}
