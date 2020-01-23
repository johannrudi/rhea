/*
 */

#include <rhea_composition.h>
#include <rhea_base.h>
#include <rhea_io_mpi.h>
#include <ymir_vec_getset.h>
#include <ymir_comm.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_COMPOSITION_DEFAULT_TYPE_NAME "NONE"
#define RHEA_COMPOSITION_DEFAULT_DATA_FILE_PATH_BIN NULL
#define RHEA_COMPOSITION_DEFAULT_DATA_FILE_PATH_TXT NULL
#define RHEA_COMPOSITION_DEFAULT_WRITE_DATA_FILE_PATH_BIN NULL
/* initialize options */
char               *rhea_composition_type_name =
  RHEA_COMPOSITION_DEFAULT_TYPE_NAME;
char               *rhea_composition_data_file_path_bin =
  RHEA_COMPOSITION_DEFAULT_DATA_FILE_PATH_BIN;
char               *rhea_composition_data_file_path_txt =
  RHEA_COMPOSITION_DEFAULT_DATA_FILE_PATH_TXT;
char               *rhea_composition_write_data_file_path_bin =
  RHEA_COMPOSITION_DEFAULT_WRITE_DATA_FILE_PATH_BIN;

void
rhea_composition_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Composition";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

    YMIR_OPTIONS_S, "type", '\0',
	&(rhea_composition_type_name), RHEA_COMPOSITION_DEFAULT_TYPE_NAME,
	"Type of composition: NONE, integer index, density field",

	YMIR_OPTIONS_S, "data-file-path-bin", '\0',
	&(rhea_composition_data_file_path_bin),
	RHEA_COMPOSITION_DEFAULT_DATA_FILE_PATH_BIN,
	"Path to a binary file that contains a temperature field",
	YMIR_OPTIONS_S, "data-file-path-txt", '\0',
	&(rhea_composition_data_file_path_txt),
	RHEA_COMPOSITION_DEFAULT_DATA_FILE_PATH_TXT,
	"Path to a text file that contains a temperature field",
	YMIR_OPTIONS_S, "write-data-file-path-bin", '\0',
	&(rhea_composition_write_data_file_path_bin),
	RHEA_COMPOSITION_DEFAULT_WRITE_DATA_FILE_PATH_BIN,
	"Output path for a binary file containing a temperature field",

  YMIR_OPTIONS_END_OF_LIST);

  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);

}

ymir_vec_t *
rhea_composition_new (ymir_mesh_t *ymir_mesh)
{
  return ymir_cvec_new (ymir_mesh, 1);
}

void
rhea_composition_destroy (ymir_vec_t *composition)
{
  ymir_vec_destroy (composition);
}

int
rhea_composition_check_vec_type (ymir_vec_t *vec)
{
  return (
      ymir_vec_is_cvec (vec) &&
      vec->ncfields == 1 &&
      vec->node_type == YMIR_GLL_NODE
  );
}

MPI_Offset *
rhea_composition_segment_offset_create (ymir_vec_t *vec)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vec);
  const mangll_locidx_t *n_nodes = ymir_mesh->cnodes->Ngo;
  const int           n_fields = vec->ncfields;
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  int                 mpisize, mpiret;
  MPI_Offset         *segment_offset;
  int                 r;

  /* get parallel environment */
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize); SC_CHECK_MPI (mpiret);

  /* create segment offsets */
  segment_offset = RHEA_ALLOC (MPI_Offset, mpisize + 1);
  segment_offset[0] = 0;
  for (r = 0; r < mpisize; r++) {
    segment_offset[r+1] = segment_offset[r] +
                          (MPI_Offset) (n_fields * n_nodes[r]);
  }

  return segment_offset;
}

int
rhea_composition_segment_size_get (ymir_vec_t *vec)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vec);
  const mangll_locidx_t *n_nodes = ymir_mesh->cnodes->Ngo;
  const int           n_fields = vec->ncfields;
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  int                 mpirank, mpiret;

  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);
  return (int) (n_fields * n_nodes[mpirank]);
}

static void
rhea_composition_read_dens (ymir_vec_t *composition)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (composition);
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  int                 mpisize, mpirank, mpiret;
  MPI_Offset         *segment_offset;
  char               *file_path_bin;
  char               *file_path_txt;
  double             *comp_data = ymir_cvec_index (composition, 0, 0);

  /* set file paths */
  if (rhea_composition_data_file_path_bin != NULL) {
    file_path_bin = rhea_composition_data_file_path_bin;
    file_path_txt = NULL;
  }
  else if (rhea_composition_data_file_path_txt != NULL) {
    file_path_bin = rhea_composition_write_data_file_path_bin;
    file_path_txt = rhea_composition_data_file_path_txt;
  }
  else { /* unknown file path */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* get parallel environment */
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize); SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);

  /* create segment offsets */
  segment_offset = rhea_composition_segment_offset_create (composition);

  /* read composition */
  if (file_path_txt != NULL || segment_offset[mpisize] < INT_MAX) {
    int                *segment_offset_int = RHEA_ALLOC (int, mpisize + 1);
    int                 r;

    for (r = 0; r <= mpisize; r++) {
      segment_offset_int[r] = (int) segment_offset[r];
    }
    rhea_io_mpi_read_scatter_double (comp_data, segment_offset_int,
                                     file_path_bin, file_path_txt, mpicomm);
    RHEA_FREE (segment_offset_int);
  }
  else {
    const int           segment_size = (int) ymir_mesh->cnodes->Ngo[mpirank];
    int                 n_read;

    RHEA_ASSERT (file_path_bin != NULL);
    n_read = rhea_io_mpi_read_segment_double (
        comp_data, segment_offset[mpirank], segment_size,
        file_path_bin, mpicomm);
    if (!(n_read == segment_size)) {
      RHEA_LERRORF ("%s: Mismatch of #entries read: "
                    "#requested %lli, #read %lli, file path %s\n",
                    __func__, (long long int) segment_size,
                    (long long int) n_read, file_path_bin);
    }
  }

  /* destroy */
  RHEA_FREE (segment_offset);

  /* communicate shared node values */
  ymir_vec_share_owned (composition);
}

int
rhea_composition_write (char *file_path_bin,
                        ymir_vec_t *composition,
                        sc_MPI_Comm mpicomm)
{
  MPI_Offset          segment_offset;
  int                 segment_size;
  int                 n_written;
  double             *data;

  /* check input */
  RHEA_ASSERT (rhea_composition_check_vec_type (composition));
  RHEA_ASSERT (sc_dmatrix_is_valid (composition->dataown));

  /* set parameters for I/O */
  data = ymir_cvec_index (composition, 0, 0);
  segment_offset = rhea_composition_segment_offset_get (composition);
  segment_size = rhea_composition_segment_size_get (composition);

  /* write */
  n_written = rhea_io_mpi_write_segment_double (
        file_path_bin, data, segment_offset, segment_size, mpicomm);
  RHEA_ASSERT (n_written == segment_size);

  /* return whether writing was success */
  return (n_written == segment_size);
}
