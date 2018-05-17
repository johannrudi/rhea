#include <rhea_plate.h>
#include <rhea_base.h>
#include <rhea_point_in_polygon.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_PLATE_DEFAULT_POLYGONS_FILE_PATH_TXT NULL

/* initialize options */
char               *rhea_plate_polygons_file_path_txt =
  RHEA_PLATE_DEFAULT_POLYGONS_FILE_PATH_TXT;

void
rhea_plate_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Plate";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_S, "polygons-file-path-txt", '\0',
    &(rhea_plate_polygons_file_path_txt),
    RHEA_PLATE_DEFAULT_POLYGONS_FILE_PATH_TXT,
    "Path to a text file with (lon,lat) coordinates of plate polygons",

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

  /* set paths to binary & text files */
  opt->polygons_file_path_txt = rhea_plate_polygons_file_path_txt;

  /* set dependent options */
  opt->domain_options = domain_options;
}

/******************************************************************************
 * TODO
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

static int
rhea_plate_get_label (const float test_x,
                      const float test_y,
                      const float **plate_vert_x,
                      const float **plate_vert_y,
                      const size_t *plate_n_vert,
                      const float **plate_vert_coarse_container_x,
                      const float **plate_vert_coarse_container_y,
                      const size_t *plate_n_vert_coarse_container,
                      const int n_plates)
{
  int                 k;

  for (k = 0; k < n_plates; k++) { /* loop over all plates */
    const float        *_sc_restrict vx = plate_vert_x[k];
    const float        *_sc_restrict vy = plate_vert_y[k];
    const size_t        nv = plate_n_vert[k];
    const float        *_sc_restrict vcx = plate_vert_coarse_container_x[k];
    const float        *_sc_restrict vcy = plate_vert_coarse_container_y[k];
    const size_t        nvc = plate_n_vert_coarse_container[k];

    if (rhea_plate_is_inside (test_x, test_y, vx, vy, nv, vcx, vcy, nvc)) {
      /* return that test coordinates are inside current plate */
      return k;
    }
  }

  /* return that test coordinates are outside of any plate */
  return RHEA_PLATE_NONE;
}

//TODO read from options
const float **plate_vert_x;
const float **plate_vert_y;
const size_t *plate_n_vert;
const float **plate_vert_coarse_container_x;
const float **plate_vert_coarse_container_y;
const size_t *plate_n_vert_coarse_container;

rhea_plate_cube_label_t
rhea_plate_cube_get_label (const double test_x, const double test_y)
{
  const int           n_plates = RHEA_PLATE_CUBE_N;

  return (rhea_plate_cube_label_t) rhea_plate_get_label (
      (float) test_x, (float) test_y, plate_vert_x, plate_vert_y, plate_n_vert,
      plate_vert_coarse_container_x, plate_vert_coarse_container_y,
      plate_n_vert_coarse_container, n_plates);
}

rhea_plate_earth_label_t
rhea_plate_earth_get_label (const double test_x, const double test_y)
{
  const int           n_plates = RHEA_PLATE_EARTH_N;

  return (rhea_plate_earth_label_t) rhea_plate_get_label (
      (float) test_x, (float) test_y, plate_vert_x, plate_vert_y, plate_n_vert,
      plate_vert_coarse_container_x, plate_vert_coarse_container_y,
      plate_n_vert_coarse_container, n_plates);
}
