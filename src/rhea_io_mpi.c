/*
 */

#include <rhea_io_mpi.h>
#include <rhea_io_std.h>
#include <rhea_base.h>

int
rhea_io_mpi_read_broadcast_double (double *values_all,
                                   int n_entries,
                                   const char *file_path_bin,
                                   const char *file_path_txt,
                                   sc_MPI_Comm mpicomm)
{
  const char         *this_fn_name = "rhea_io_mpi_read_broadcast_double";
  int                 mpirank, mpiret;
  size_t              total_size = (size_t) n_entries;

  if (file_path_txt == NULL) { /* if read from binary file */
    RHEA_GLOBAL_INFOF ("Into %s (%s, #entries %i)\n", this_fn_name,
                       file_path_bin, n_entries);
  }
  else {
    RHEA_GLOBAL_INFOF ("Into %s (bin: %s, txt: %s)\n", this_fn_name,
                       file_path_bin, file_path_txt);
  }

  /* check input */
  RHEA_ASSERT (file_path_txt != NULL ||
               (file_path_bin != NULL && 0 < n_entries));

  /* get parallel environment */
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);

  /* read file with one processor */
  if (mpirank == 0) {
    if (file_path_txt == NULL) { /* if read from binary file */
      rhea_io_std_read_double (values_all, total_size, file_path_bin);
    }
    else { /* otherwise read from text file */
      total_size = rhea_io_std_read_double_from_txt (values_all, total_size,
                                                     file_path_txt);
      if (file_path_bin != NULL) {
        rhea_io_std_write_double (file_path_bin, values_all, total_size);
      }
    }
  }

  /* broadcast values to all processors */
  if (n_entries <= 0) { /* if #entries is not known */
    n_entries = (int) total_size;
    mpiret = sc_MPI_Bcast (&n_entries, 1, sc_MPI_INT, 0, mpicomm);
    SC_CHECK_MPI (mpiret);
  }
  mpiret = sc_MPI_Bcast (values_all, n_entries, sc_MPI_DOUBLE, 0, mpicomm);
  SC_CHECK_MPI (mpiret);

  if (file_path_txt == NULL) { /* if read from binary file */
    RHEA_GLOBAL_INFOF ("Done %s (%s, #entries %i)\n", this_fn_name,
                       file_path_bin, n_entries);
  }
  else {
    RHEA_GLOBAL_INFOF ("Done %s (bin: %s, txt: %s)\n", this_fn_name,
                       file_path_bin, file_path_txt);
  }

  /* return number of entries read */
  return n_entries;
}

void
rhea_io_mpi_read_scatter_double (double *values_segment,
                                 const int *segment_offset,
                                 const char *file_path_bin,
                                 const char *file_path_txt,
                                 sc_MPI_Comm mpicomm)
{
  const char         *this_fn_name = "rhea_io_mpi_read_scatter_double";
  int                 mpisize, mpirank, mpiret;
  double             *values_all;
  int                *segment_size;
  int                 r;

  /* get parallel environment */
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize); SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);

  if (file_path_txt == NULL) { /* if read from binary file */
    RHEA_GLOBAL_INFOF ("Into %s (%s, #entries %i)\n", this_fn_name,
                       file_path_bin, segment_offset[mpisize]);
  }
  else {
    RHEA_GLOBAL_INFOF ("Into %s (bin: %s, txt: %s)\n", this_fn_name,
                       file_path_bin, file_path_txt);
  }

  /* check input */
  RHEA_ASSERT (file_path_txt != NULL || file_path_bin != NULL);
  RHEA_ASSERT (0 < segment_offset[mpisize]);

  /* read file with one processor */
  if (mpirank == 0) {
    size_t              total_size = (size_t) segment_offset[mpisize];

    values_all = RHEA_ALLOC (double, total_size);
    if (file_path_txt == NULL) { /* if read from binary file */
      rhea_io_std_read_double (values_all, total_size, file_path_bin);
    }
    else { /* otherwise read from text file */
      total_size = rhea_io_std_read_double_from_txt (values_all, total_size,
                                                     file_path_txt);
      if (file_path_bin != NULL) {
        rhea_io_std_write_double (file_path_bin, values_all, total_size);
      }
    }
  }
  else {
    values_all = NULL;
  }

  /* create segment size array */
  segment_size = RHEA_ALLOC (int, mpisize);
  for (r = 0; r < mpisize; r++) {
    segment_size[r] = segment_offset[r+1] - segment_offset[r];
  }

  /* scatter segments to corresponding processors */
  mpiret = MPI_Scatterv (
      values_all, segment_size, segment_offset, sc_MPI_DOUBLE,
      values_segment, segment_size[mpirank], sc_MPI_DOUBLE, 0, mpicomm);
  SC_CHECK_MPI (mpiret);
  //TODO `sc_MPI_Scatterv` is not defined

  /* destroy */
  if (mpirank == 0) {
    RHEA_FREE (values_all);
  }
  RHEA_FREE (segment_size);

  if (file_path_txt == NULL) { /* if read from binary file */
    RHEA_GLOBAL_INFOF ("Done %s (%s)\n", this_fn_name, file_path_bin);
  }
  else {
    RHEA_GLOBAL_INFOF ("Done %s (bin: %s, txt: %s)\n", this_fn_name,
                       file_path_bin, file_path_txt);
  }
}
