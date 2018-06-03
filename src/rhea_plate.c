#include <rhea_plate.h>
#include <rhea_base.h>
#include <rhea_point_in_polygon.h>
#include <rhea_io_mpi.h>
#include <rhea_velocity.h>
#include <ymir_velocity_vec.h>
#include <ymir_mass_vec.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_PLATE_DEFAULT_POLYGON_VERTICES_FILE_PATH_TXT NULL
#define RHEA_PLATE_DEFAULT_POLYGON_N_VERTICES_TOTAL (0)
#define RHEA_PLATE_DEFAULT_POLYGON_VERTICES_COARSE_FILE_PATH_TXT NULL
#define RHEA_PLATE_DEFAULT_POLYGON_N_VERTICES_COARSE_TOTAL (0)
#define RHEA_PLATE_DEFAULT_POLYGON_X_MIN (0.0)
#define RHEA_PLATE_DEFAULT_POLYGON_X_MAX (360.0)
#define RHEA_PLATE_DEFAULT_POLYGON_Y_MIN (0.0)
#define RHEA_PLATE_DEFAULT_POLYGON_Y_MAX (180.0)

/* initialize options */
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

void
rhea_plate_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Plate";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

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

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);
}

void
rhea_plate_process_options (rhea_plate_options_t *opt,
                            rhea_domain_options_t *domain_options)
{
  /* exit if nothing to do */
  if (opt == NULL) {
    return;
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

  /* set dependent options */
  opt->domain_options = domain_options;
}

/******************************************************************************
 * Plate Boundary Data
 *****************************************************************************/

static int
rhea_plate_get_n_plates (rhea_plate_options_t *opt)
{
  switch (opt->domain_options->shape) {
  case RHEA_DOMAIN_CUBE:  return RHEA_PLATE_CUBE_N;
  case RHEA_DOMAIN_SHELL: return RHEA_PLATE_EARTH_N;
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }
}

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
    RHEA_GLOBAL_VERBOSEF (" %s Polygon of plate %i has %i vertices\n",
                          __func__, pid, (int) n_vert);
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
rhea_plate_polygon_translation_shell_init (int bins_x[36],
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
rhea_plate_polygon_translation_shell_set (float *translation_x,
                                          float *translation_y,
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
    if (bin_north < 90 && bin_south < 90) { /* if polygon in south hemisph. */
      bin_north = 0;
    }
    else if (90 < bin_north && 90 < bin_south) { /* if polygon in north hemi. */
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
  }
  else if (0 < bin_north && bin_south < 179) { /* set translation in west dir */
    RHEA_ASSERT (bin_east < bin_west);
    *translation_x = -10.0 * (double) bin_west;
    *translation_y = 0.0;
  }
  else if (0 == bin_north) { /* set translation in south direction */
    RHEA_ASSERT (bin_south <= 8);
    *translation_x = 0.0;
  //*translation_y = +10.0 * (double) bin_south;
    *translation_y = 0.0; //TODO add pts to data: (0,0) and (360,0)
  }
  else if (bin_south == 179) { /* set translation in south direction */
    RHEA_ASSERT (9 <= bin_north);
    *translation_x = 0.0;
  //*translation_y = -10.0 * (double) bin_north;
    *translation_y = 0.0; //TODO add pts to data: (0,180) and (360,180)
  }
  else {
    RHEA_ABORT_NOT_REACHED ();
  }

  RHEA_GLOBAL_VERBOSEF (
      " %s Polygon bounds: west %i, east %i, north %i, south %i; "
      "translation: x %g, y %g\n",
      __func__, 10*bin_west, 10*bin_east, bin_north, bin_south,
      *translation_x, *translation_y);
}

static void
rhea_plate_polygon_translation_shell_apply (float *vertices_x,
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

static void
rhea_plate_translate_vertices_shell (float *translation_x,
                                     float *translation_y,
                                     float **vertices_x,
                                     float **vertices_y,
                                     const size_t *n_vertices,
                                     const int n_plates)
{
  int                 bins_x[36], bins_y[180];
  int                 pid;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (__func__, "n_plates=%i", n_plates);

  for (pid = 0; pid < n_plates; pid++) {
    /* populate bins */
    rhea_plate_polygon_translation_shell_init (
        bins_x, bins_y, vertices_x[pid], vertices_y[pid], n_vertices[pid]);

    /* set translations */
    rhea_plate_polygon_translation_shell_set (
        &translation_x[pid], &translation_y[pid], bins_x, bins_y);

    /* apply translations */
    rhea_plate_polygon_translation_shell_apply (
        vertices_x[pid], vertices_y[pid], n_vertices[pid],
        translation_x[pid], translation_y[pid]);
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
    RHEA_GLOBAL_VERBOSEF (" %s Polygon of plate %i has %i vertices\n",
                          __func__, pid, (int) n_vert);
    RHEA_ASSERT (0 < n_vert);

    /* set vertex coordinates */
    vert_x = vertices_x[pid] = RHEA_ALLOC (float, n_vert);
    vert_y = vertices_y[pid] = RHEA_ALLOC (float, n_vert);
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
    }

    /* increment vertex counter */
    idx_curr += (n_vert + 1);
    n_total -= (n_vert + 1);
  }

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

int
rhea_plate_data_create (rhea_plate_options_t *opt, sc_MPI_Comm mpicomm)
{
  int                 n_plates;
  int                 n_vertices_total, n_vertices_coarse_total;
  float              *vertices_all;
  int                 success_vert = 0;
  int                 success_vert_coarse = 0;

  /* exit if nothing to do */
  if (opt == NULL) {
    return 1;
  }

  /* get parameters from options */
  n_plates = rhea_plate_get_n_plates (opt);
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
        vertices_all, 2*n_vertices_total, NULL /* path bin */,
        opt->vertices_coarse_container_file_path_txt, mpicomm);

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
  switch (opt->domain_options->shape) {
  case RHEA_DOMAIN_CUBE:
    break;
  case RHEA_DOMAIN_SHELL:
    opt->translation_x = RHEA_ALLOC (float, n_plates);
    opt->translation_y = RHEA_ALLOC (float, n_plates);
    if (opt->vertices_coarse_container_x != NULL &&
        opt->vertices_coarse_container_y != NULL &&
        opt->n_vertices_coarse_container != NULL) { /* use coarse vertices */
      int                 pid;

      rhea_plate_translate_vertices_shell (
          opt->translation_x, opt->translation_y,
          opt->vertices_coarse_container_x, opt->vertices_coarse_container_y,
          opt->n_vertices_coarse_container, n_plates);
      for (pid = 0; pid < n_plates; pid++) {
        rhea_plate_polygon_translation_shell_apply (
            opt->vertices_x[pid], opt->vertices_y[pid], opt->n_vertices[pid],
            opt->translation_x[pid], opt->translation_y[pid]);
      }
    }
    else { /* use vertices of (fine) polygon */
      rhea_plate_translate_vertices_shell (
          opt->translation_x, opt->translation_y,
          opt->vertices_x, opt->vertices_y, opt->n_vertices, n_plates);
    }
    break;
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }

  return (success_vert + success_vert_coarse);
}

void
rhea_plate_data_clear (rhea_plate_options_t *opt)
{
  int                 pid;

  /* exit if nothing to do */
  if (opt == NULL) {
    return;
  }

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
}

static int
rhea_plate_data_exitst (rhea_plate_options_t *opt)
{
  return (opt->vertices_x != NULL && opt->vertices_y != NULL &&
          opt->n_vertices != NULL);
}

/******************************************************************************
 * Plate Retrieval
 *****************************************************************************/

/**
 * Checks whether the given (x,y) coordinates are inside a polygon.
 */
static int
rhea_plate_is_inside_polygon (const float test_x,
                              const float test_y,
                              const float *_sc_restrict vert_x,
                              const float *_sc_restrict vert_y,
                              const size_t n_vert,
                              const float *_sc_restrict vert_coarse_container_x,
                              const float *_sc_restrict vert_coarse_container_y,
                              const size_t n_vert_coarse_container)
{
  /* check input */
  RHEA_ASSERT (vert_x != NULL);
  RHEA_ASSERT (vert_y != NULL);
  RHEA_ASSERT (0 < n_vert);
  RHEA_ASSERT (0 >= n_vert_coarse_container || vert_coarse_container_x != NULL);
  RHEA_ASSERT (0 >= n_vert_coarse_container || vert_coarse_container_y != NULL);

  /* check coarse container of plate */
  if (0 < n_vert_coarse_container) {
    int                 is_inside;

    /* run point-in-polygon test on coarse container */
    is_inside = rhea_point_in_polygon_is_inside (
        test_x, test_y, vert_coarse_container_x, vert_coarse_container_y,
        n_vert_coarse_container);

    /* exit if point is outside of container */
    if (!is_inside) {
      return 0;
    }
  }

  /* check fine-resolution boundary of plate */
  return rhea_point_in_polygon_is_inside (test_x, test_y, vert_x, vert_y,
                                          n_vert);
}

/**
 * Checks whether the given (x,y,z) coordinates are inside a specific plate.
 */
static int
rhea_plate_is_inside (const double x, const double y, const double z,
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

  return rhea_plate_is_inside_polygon (test_x, test_y, vx, vy, nv,
                                       vcx, vcy, nvc);
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

  /* check input */
  RHEA_ASSERT (opt != NULL);

  /* exit if nothing to do */
  if (!rhea_plate_data_exitst (opt)) {
    ymir_vec_set_value (vec, (double) RHEA_PLATE_NONE);
    return;
  }

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

  /* check input */
  RHEA_ASSERT (opt != NULL);

  /* exit if nothing to do */
  if (!rhea_plate_data_exitst (opt)) {
    return;
  }

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

  /* check input */
  RHEA_ASSERT (opt != NULL);

  /* exit if nothing to do */
  if (!rhea_plate_data_exitst (opt)) {
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
      data.filter = ymir_face_cvec_new (ymir_mesh, 1, vec->meshnum);
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
 * Plate Velocities
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
  const double  thermal_diffus_m2_s = 1.0e-6; //TODO take from temp options
  const double  radius_m = opt->domain_options->radius_max_m;

  if (!remove_dim) {
    RHEA_ASSERT (0.0 <= rot_spherical[0]);
    RHEA_ASSERT (0.0 <= rot_spherical[1] && rot_spherical[1] <= 2.0*M_PI);
    RHEA_ASSERT (0.0 <= rot_spherical[2] && rot_spherical[2] <= M_PI);
    rot_spherical[0] = rot_spherical[0]/M_PI * 180.0;
    rot_spherical[0] /= radius_m / (thermal_diffus_m2_s/radius_m);
    rot_spherical[0] /= 1.0 / (1.0e6*sec_per_yr);
    rot_spherical[1] = rot_spherical[1]/M_PI * 180.0 - 180.0;
    rot_spherical[2] = -rot_spherical[2]/M_PI * 180.0 + 90.0;
  }
  else {
    RHEA_ASSERT (0.0 <= rot_spherical[0]);
    RHEA_ASSERT (-180.0 <= rot_spherical[1] && rot_spherical[1] <= 180.0);
    RHEA_ASSERT (- 90.0 <= rot_spherical[2] && rot_spherical[2] <=  90.0);
    rot_spherical[0] = rot_spherical[0]/180.0 * M_PI;
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
rhea_plate_velocity_get_angular_velocity (double rot_axis[3],
                                          const int plate_label,
                                          rhea_plate_options_t *opt)
{
  /* check input */
  RHEA_ASSERT (RHEA_PLATE_NONE < plate_label &&
               plate_label < rhea_plate_get_n_plates (opt));

  switch (opt->domain_options->shape) {
  case RHEA_DOMAIN_CUBE:
    /* get values */
    rot_axis[0] = rhea_plate_velocity_cube[3*plate_label    ];
    rot_axis[1] = rhea_plate_velocity_cube[3*plate_label + 1];
    rot_axis[2] = rhea_plate_velocity_cube[3*plate_label + 2];
    break;

  case RHEA_DOMAIN_SHELL:
    {
      double        rot_sp[3];

      /* get values */
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

void
rhea_plate_velocity_generate (ymir_vec_t *vel, const int plate_label,
                              rhea_plate_options_t *opt)
{
  double              rot_center[3], rot_axis[3];

  /* check input */
  RHEA_ASSERT (opt != NULL);
  RHEA_ASSERT (RHEA_PLATE_NONE <= plate_label &&
               plate_label < rhea_plate_get_n_plates (opt));

  /* exit if nothing to do */
  if (!rhea_plate_data_exitst (opt)) {
    ymir_vec_set_value (vel, NAN);
    return;
  }
  else if (RHEA_PLATE_NONE == plate_label ) {
    ymir_vec_set_value (vel, 0.0);
    return;
  }

  /* get center from options */
  rot_center[0] = opt->domain_options->center[0];
  rot_center[1] = opt->domain_options->center[1];
  rot_center[2] = opt->domain_options->center[2];

  /* assign angular velocity */
  rhea_plate_velocity_get_angular_velocity (rot_axis, plate_label, opt);

  /* generate velocity field and apply filter for the plate's interior */
  ymir_velocity_vec_generate_rotation (vel, rot_center, rot_axis);
  rhea_plate_apply_filter_vec (vel, plate_label, opt);
}

void
rhea_plate_velocity_generate_all (ymir_vec_t *vel, rhea_plate_options_t *opt)
{
  ymir_vec_t         *plate_vel = ymir_vec_template (vel);
  int                 pid;

  /* check input */
  RHEA_ASSERT (opt != NULL);

  /* initialize output */
  ymir_vec_set_value (vel, 0.0);

  /* get velocities of all plates */
  for (pid = 0; pid < rhea_plate_get_n_plates (opt); pid++) {
    rhea_plate_velocity_generate (plate_vel, pid, opt);
    ymir_vec_add (1.0, plate_vel, vel);
  }

  /* destroy */
  ymir_vec_destroy (plate_vel);
}

double
rhea_plate_velocity_get_mean_magnitude (ymir_vec_t *vel,
                                        const int plate_label,
                                        rhea_plate_options_t *opt)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vel);
  const ymir_topidx_t face_surf = RHEA_DOMAIN_BOUNDARY_FACE_TOP;
  ymir_vec_t         *vel_surf;
  ymir_vec_t         *magn_surf, *filter_surf, *tmp_mass;
  double              vel_mean;

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (vel) ||
               rhea_velocity_surface_check_vec_type (vel));
  RHEA_ASSERT (opt != NULL);

  /* exit if nothing to do */
  if (!rhea_plate_data_exitst (opt)) {
    return NAN;
  }

  /* create surface velocity without net rotations */
  if (!ymir_vec_is_face_vec (vel)) { /* if volume vector */
    vel_surf = rhea_velocity_surface_new_from_vol (vel);
  }
  else { /* if face vector */
    vel_surf = rhea_velocity_surface_new (ymir_mesh);
    ymir_vec_copy (vel, vel_surf);
  }
  rhea_plate_velocity_project_out_mean_rotation (vel_surf, opt);

  /* compute pointwise magnitude */
  magn_surf = ymir_face_cvec_new (ymir_mesh, 1, face_surf);
  ymir_cvec_innerprod_pointwise (magn_surf, vel_surf, vel_surf);
  ymir_cvec_sqrt (magn_surf, magn_surf);
  rhea_velocity_surface_destroy (vel_surf);

  /* generate filter for plate; apply to velocity vector */
  filter_surf = ymir_face_cvec_new (ymir_mesh, 1, face_surf);
  ymir_vec_set_value (filter_surf, 1.0);
  rhea_plate_apply_filter_vec (filter_surf, plate_label, opt);
  ymir_vec_multiply_in1 (filter_surf, magn_surf);

  /* compute mean velocity */
  tmp_mass = ymir_face_cvec_new (ymir_mesh, 1, face_surf);
  ymir_mass_apply (magn_surf, tmp_mass);
  vel_mean = ymir_vec_innerprod (tmp_mass, magn_surf);
  ymir_mass_apply (filter_surf, tmp_mass);
  vel_mean /= ymir_vec_innerprod (tmp_mass, filter_surf);

  /* destroy */
  ymir_vec_destroy (magn_surf);
  ymir_vec_destroy (filter_surf);
  ymir_vec_destroy (tmp_mass);

  /* return mean velocity */
  return vel_mean;
}

void
rhea_plate_velocity_project_out_mean_rotation (ymir_vec_t *vel,
                                               rhea_plate_options_t *opt)
{
  double              rot_center[3], moment_of_inertia[3];

  /* check input */
  RHEA_ASSERT (opt != NULL);

  /* exit if nothing to do */
  if (!rhea_plate_data_exitst (opt)) {
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
      vel, rot_center, moment_of_inertia, 0 /*residual_space*/);
}
