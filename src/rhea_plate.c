#include <rhea_plate.h>
#include <rhea_base.h>
#include <rhea_point_in_polygon.h>
#include <rhea_io_mpi.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_PLATE_DEFAULT_POLYGON_VERTICES_FILE_PATH_TXT NULL
#define RHEA_PLATE_DEFAULT_POLYGON_N_VERTICES_TOTAL (0)
#define RHEA_PLATE_DEFAULT_POLYGON_VERTICES_COARSE_FILE_PATH_TXT NULL
#define RHEA_PLATE_DEFAULT_POLYGON_N_VERTICES_COARSE_TOTAL (0)

/* initialize options */
char               *rhea_plate_polygon_vertices_file_path_txt =
  RHEA_PLATE_DEFAULT_POLYGON_VERTICES_FILE_PATH_TXT;
int                 rhea_plate_polygon_n_vertices_total =
  RHEA_PLATE_DEFAULT_POLYGON_N_VERTICES_TOTAL;
char               *rhea_plate_polygon_vertices_coarse_file_path_txt =
  RHEA_PLATE_DEFAULT_POLYGON_VERTICES_COARSE_FILE_PATH_TXT;
int                 rhea_plate_polygon_n_vertices_coarse_total =
  RHEA_PLATE_DEFAULT_POLYGON_N_VERTICES_COARSE_TOTAL;

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
  opt->translation_x = NULL;
  opt->translation_y = NULL;
  opt->vertices_coarse_container_x = NULL;
  opt->vertices_coarse_container_y = NULL;
  opt->n_vertices_coarse_container = NULL;
  opt->translation_coarse_container_x = NULL;
  opt->translation_coarse_container_y = NULL;

  /* set dependent options */
  opt->domain_options = domain_options;
}

/******************************************************************************
 * Data
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
                                  const int n_vertices_total)
{
  int                 idx_curr = 0;
  int                 n_total = n_vertices_total;
  int                 pid;
  float              *vert_x, *vert_y;
  size_t              n_vert, vid;

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
    RHEA_GLOBAL_VERBOSEF ("%s: Polygon of plate %i has %i vertices\n",
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
    }

    /* increment vertex counter */
    idx_curr += (n_vert + 1);
    n_total -= (n_vert + 1);
  }
}

int
rhea_plate_data_create (rhea_plate_options_t *opt, sc_MPI_Comm mpicomm)
{
  int                 n_plates;
  int                 n_vertices_total, n_vertices_coarse_total;
  float              *vertices_all;
  int                 success = 0;

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
    opt->vertices_x = RHEA_ALLOC (float *, n_plates);
    opt->vertices_y = RHEA_ALLOC (float *, n_plates);
    opt->n_vertices = RHEA_ALLOC (size_t, n_plates);

    /* read file */
    vertices_all = RHEA_ALLOC (float, 2*n_vertices_total);
    success += rhea_io_mpi_read_broadcast_float (
        vertices_all, 2*n_vertices_total, NULL /* path bin */,
        opt->vertices_file_path_txt, mpicomm);

    /* process vertices */
    switch (opt->domain_options->shape) {
    case RHEA_DOMAIN_CUBE:
      rhea_plate_process_vertices_cube (
        opt->vertices_x, opt->vertices_y, opt->n_vertices, n_plates,
        vertices_all, n_vertices_total);
      break;
    case RHEA_DOMAIN_SHELL:
      //TODO
      RHEA_ABORT_NOT_REACHED ();
      break;
    default: /* unknown domain shape */
      RHEA_ABORT_NOT_REACHED ();
    }

    /* destroy */
    RHEA_FREE (vertices_all);
  }
  else {
    success++;
  }

  /* read coarse container vertices (if they do not exist) */
  if (opt->vertices_coarse_container_file_path_txt != NULL &&
      opt->vertices_coarse_container_x == NULL &&
      opt->vertices_coarse_container_y == NULL &&
      opt->n_vertices_coarse_container == NULL) {
    opt->vertices_coarse_container_x = RHEA_ALLOC (float *, n_plates);
    opt->vertices_coarse_container_y = RHEA_ALLOC (float *, n_plates);
    opt->n_vertices_coarse_container = RHEA_ALLOC (size_t, n_plates);

    /* read file */
    vertices_all = RHEA_ALLOC (float, 2*n_vertices_coarse_total);
    success += rhea_io_mpi_read_broadcast_float (
        vertices_all, 2*n_vertices_total, NULL /* path bin */,
        opt->vertices_coarse_container_file_path_txt, mpicomm);

    /* process vertices */
    switch (opt->domain_options->shape) {
    case RHEA_DOMAIN_CUBE:
      rhea_plate_process_vertices_cube (
        opt->vertices_coarse_container_x, opt->vertices_coarse_container_y,
        opt->n_vertices_coarse_container, n_plates,
        vertices_all, n_vertices_total);
      break;
    case RHEA_DOMAIN_SHELL:
      //TODO
      RHEA_ABORT_NOT_REACHED ();
      break;
    default: /* unknown domain shape */
      RHEA_ABORT_NOT_REACHED ();
    }

    /* destroy */
    RHEA_FREE (vertices_all);
  }
  else {
    success++;
  }

  return success;
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
  if (opt->translation_x != NULL && opt->translation_y != NULL) {
    RHEA_FREE (opt->translation_x);
    RHEA_FREE (opt->translation_y);
    opt->translation_x = NULL;
    opt->translation_y = NULL;
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
  if (opt->translation_coarse_container_x != NULL &&
      opt->translation_coarse_container_y != NULL) {
    RHEA_FREE (opt->translation_coarse_container_x);
    RHEA_FREE (opt->translation_coarse_container_y);
    opt->translation_coarse_container_x = NULL;
    opt->translation_coarse_container_y = NULL;
  }
}

/******************************************************************************
 * Plate Retrieval
 *****************************************************************************/

static int
rhea_plate_is_inside (const float test_x,
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

int
rhea_plate_get_label (const double x, const double y, const double z,
                      rhea_plate_options_t *opt)
{
  const int           n_plates = rhea_plate_get_n_plates (opt);
  int                 pid;
  float               test_x, test_y;
  float              *vx, *vy, *vcx = NULL, *vcy = NULL;
  size_t              nv, nvc = 0;

  /* check input */
  RHEA_ASSERT (opt->vertices_x != NULL);
  RHEA_ASSERT (opt->vertices_y != NULL);
  RHEA_ASSERT (opt->n_vertices != NULL);

  /* transform test coordinates */
  switch (opt->domain_options->shape) {
  case RHEA_DOMAIN_CUBE:
    //TODO apply a coordinate mapping from `rhea_domain` (maybe implement that)
    test_x = (float) x;
    test_y = (float) y;
    break;
  case RHEA_DOMAIN_SHELL:
    //TODO convert to spherical and apply shift
    //test_x = (float) x;
    //test_y = (float) y;
    break;
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* find and return plate label */
  for (pid = 0; pid < n_plates; pid++) { /* loop over all plates */
    vx = opt->vertices_x[pid];
    vy = opt->vertices_y[pid];
    nv = opt->n_vertices[pid];
    if (opt->vertices_coarse_container_x != NULL &&
        opt->vertices_coarse_container_y != NULL &&
        opt->n_vertices_coarse_container != NULL) {
      vcx = opt->vertices_coarse_container_x[pid];
      vcy = opt->vertices_coarse_container_y[pid];
      nvc = opt->n_vertices_coarse_container[pid];
    }

    if (rhea_plate_is_inside (test_x, test_y, vx, vy, nv, vcx, vcy, nvc)) {
      /* return that test coordinates are inside current plate */
      return pid;
    }
  }

  /* return that test coordinates are outside of any plate */
  return RHEA_PLATE_NONE;
}

/*TODO delete
const float **plate_vert_x;
const float **plate_vert_y;
const size_t *plate_n_vert;
const float **plate_vert_coarse_container_x;
const float **plate_vert_coarse_container_y;
const size_t *plate_n_vert_coarse_container;

rhea_plate_cube_label_t
rhea_plate_cube_get_label (const double test_x, const double test_y,
                           rhea_plate_options_t *opt)
{
  const int           n_plates = RHEA_PLATE_CUBE_N;

  return (rhea_plate_cube_label_t) rhea_plate_get_label (
      (float) test_x, (float) test_y, plate_vert_x, plate_vert_y, plate_n_vert,
      plate_vert_coarse_container_x, plate_vert_coarse_container_y,
      plate_n_vert_coarse_container, n_plates);
}

rhea_plate_earth_label_t
rhea_plate_earth_get_label (const double test_x, const double test_y,
                            rhea_plate_options_t *opt)
{
  const int           n_plates = RHEA_PLATE_EARTH_N;

  return (rhea_plate_earth_label_t) rhea_plate_get_label (
      (float) test_x, (float) test_y, plate_vert_x, plate_vert_y, plate_n_vert,
      plate_vert_coarse_container_x, plate_vert_coarse_container_y,
      plate_n_vert_coarse_container, n_plates);
}
*/
