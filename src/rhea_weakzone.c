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
#define RHEA_WEAKZONE_DEFAULT_POINTS_FILE_PATH_BIN NULL
#define RHEA_WEAKZONE_DEFAULT_POINTS_FILE_PATH_TXT NULL
#define RHEA_WEAKZONE_DEFAULT_N_POINTS (5000000)
#define RHEA_WEAKZONE_DEFAULT_WRITE_POINTS_FILE_PATH_BIN NULL
#define RHEA_WEAKZONE_DEFAULT_WRITE_POINTS_FILE_PATH_TXT NULL

/* initialize options */
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
rhea_weakzone_process_options (rhea_weakzone_options_t *opt)
{
  /* set weak zone type */
  //TODO
//if (strcmp (rhea_viscosity_model_name, "UWYL") == 0) {
//  opt->model = RHEA_VISCOSITY_MODEL_UWYL;
//}
//else if (strcmp (rhea_viscosity_model_name, "UWYL_LADD_UCUT") == 0) {
//  opt->model = RHEA_VISCOSITY_MODEL_UWYL_LADD_UCUT;
//}
//else if (strcmp (rhea_viscosity_model_name, "UWYL_LADD_USHIFT") == 0) {
//  opt->model = RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT;
//}
//else { /* unknown model name */
//  RHEA_ABORT ("Unknown viscosity model name");
//}

  /* set paths to binary & text files of (x,y,z) coordinates */
  opt->points_file_path_bin = rhea_weakzone_points_file_path_bin;
  opt->points_file_path_txt = rhea_weakzone_points_file_path_txt;

  /* set number of points */
  opt->n_points = rhea_weakzone_n_points;

  /* set output paths */
  opt->write_points_file_path_bin = rhea_weakzone_write_points_file_path_bin;
  opt->write_points_file_path_txt = rhea_weakzone_write_points_file_path_txt;

  /* init data */
  opt->coordinates = NULL;
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
 * Get & Set Functions
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
  const char         *file_path_bin = opt->points_file_path_bin;
  const char         *file_path_txt = opt->points_file_path_txt;
  const char         *write_file_path_bin = opt->write_points_file_path_bin;
  const char         *write_file_path_txt = opt->write_points_file_path_txt;
  const int           n_entries = 3 * opt->n_points;
  int                 n_read;
  double             *coordinates;

  RHEA_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (0 < opt->n_points);

  /* create coordinates array */
  opt->coordinates = coordinates = RHEA_ALLOC (double, n_entries);

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

  RHEA_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

void
rhea_weakzone_data_clear (rhea_weakzone_options_t *opt)
{
  /* destroy */
  if (opt->coordinates != NULL) {
    RHEA_FREE (opt->coordinates);
  }
}
