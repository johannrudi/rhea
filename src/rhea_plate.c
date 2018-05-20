#include <rhea_plate.h>
#include <rhea_base.h>
#include <rhea_point_in_polygon.h>
#include <rhea_io_mpi.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_PLATE_DEFAULT_POLYGON_VERTICES_FILE_PATH_TXT NULL
#define RHEA_PLATE_DEFAULT_POLYGON_VERTICES_COARSE_FILE_PATH_TXT NULL
#define RHEA_PLATE_DEFAULT_POLYGON_N_VERTICES_TOTAL (0)

/* initialize options */
char               *rhea_plate_polygon_vertices_file_path_txt =
  RHEA_PLATE_DEFAULT_POLYGON_VERTICES_FILE_PATH_TXT;
char               *rhea_plate_polygon_vertices_coarse_file_path_txt =
  RHEA_PLATE_DEFAULT_POLYGON_VERTICES_COARSE_FILE_PATH_TXT;
int                 rhea_plate_polygon_n_vertices_total =
  RHEA_PLATE_DEFAULT_POLYGON_N_VERTICES_TOTAL;

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

  YMIR_OPTIONS_S, "polygon-vertices-coarse-file-path-txt", '\0',
    &(rhea_plate_polygon_vertices_coarse_file_path_txt),
    RHEA_PLATE_DEFAULT_POLYGON_VERTICES_COARSE_FILE_PATH_TXT,
    "Path to a text file with vertices to coarse polygon containers",

  YMIR_OPTIONS_I, "polygon-vertices-num-total", '\0',
    &(rhea_plate_polygon_n_vertices_total),
    RHEA_PLATE_DEFAULT_POLYGON_N_VERTICES_TOTAL,
    "Total number of vertices of all combined plate polygons",

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
  opt->vertices_coarse_container_file_path_txt =
    rhea_plate_polygon_vertices_coarse_file_path_txt;
  opt->n_vertices_total = rhea_plate_polygon_n_vertices_total;

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

int
rhea_plate_data_create (rhea_plate_options_t *opt, sc_MPI_Comm mpicomm)
{
  int                 n_vertices_total;
  int                 n_plates, k;
  float              *vert = NULL;
  int                 success = 0;

  /* exit if nothing to do */
  if (opt == NULL) {
    return 1;
  }

  /* get parameters from options */
  n_vertices_total = opt->n_vertices_total;
  n_plates = rhea_plate_get_n_plates (opt);

  /* read (fine) polygon vertices (if they do not exist) */
  if (opt->vertices_file_path_txt != NULL &&
      opt->vertices_x == NULL && opt->vertices_y == NULL &&
      opt->n_vertices == NULL) {
    opt->vertices_x = RHEA_ALLOC (float *, n_plates);
    opt->vertices_y = RHEA_ALLOC (float *, n_plates);
    opt->n_vertices = RHEA_ALLOC (size_t, n_plates);

    /* read file */
    vert = RHEA_ALLOC (float, 2*n_vertices_total);
    success += rhea_io_mpi_read_broadcast_float (
        vert, n_vertices_total, NULL /* path bin */,
        opt->vertices_file_path_txt, mpicomm);

    /* fill vertex arrays */
    for (k = 0; k < n_plates; k++) {
      //TODO
    }
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
    if (vert == NULL) {
      vert = RHEA_ALLOC (float, 2*n_vertices_total);
    }

    /* fill vertex arrays */
    for (k = 0; k < n_plates; k++) {
    }
  }
  else {
    success++;
  }

  /* destroy */
  if (vert != NULL) {
    RHEA_FREE (vert);
  }

  return 1;
}

void
rhea_plate_data_clear (rhea_plate_options_t *opt)
{
  int                 n_plates, k;

  /* exit if nothing to do */
  if (opt == NULL) {
    return;
  }

  /* get parameters from options */
  n_plates = rhea_plate_get_n_plates (opt);

  /* destroy polygon vertices */
  if (opt->vertices_x != NULL && opt->vertices_y != NULL &&
      opt->n_vertices != NULL) {
    RHEA_FREE (opt->vertices_x);
    RHEA_FREE (opt->vertices_y);
    RHEA_FREE (opt->n_vertices);
    opt->vertices_x = NULL;
    opt->vertices_y = NULL;
    opt->n_vertices = NULL;
    opt->translation_x = NULL;
    opt->translation_y = NULL;
  }

  /* destroy container vertices */
  if (opt->vertices_coarse_container_x != NULL &&
      opt->vertices_coarse_container_y != NULL &&
      opt->n_vertices_coarse_container != NULL) {
    RHEA_FREE (opt->vertices_coarse_container_x);
    RHEA_FREE (opt->vertices_coarse_container_y);
    RHEA_FREE (opt->n_vertices_coarse_container);
    opt->vertices_coarse_container_x = NULL;
    opt->vertices_coarse_container_y = NULL;
    opt->n_vertices_coarse_container = NULL;
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
  int                 k;
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
  for (k = 0; k < n_plates; k++) { /* loop over all plates */
    vx = opt->vertices_x[k];
    vy = opt->vertices_y[k];
    nv = opt->n_vertices[k];
    if (opt->vertices_coarse_container_x != NULL &&
        opt->vertices_coarse_container_y != NULL &&
        opt->n_vertices_coarse_container != NULL) {
      vcx = opt->vertices_coarse_container_x[k];
      vcy = opt->vertices_coarse_container_y[k];
      nvc = opt->n_vertices_coarse_container[k];
    }

    if (rhea_plate_is_inside (test_x, test_y, vx, vy, nv, vcx, vcy, nvc)) {
      /* return that test coordinates are inside current plate */
      return k;
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
