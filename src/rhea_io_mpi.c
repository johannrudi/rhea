/*
 */

#include <rhea_io_mpi.h>
#include <rhea_io_std.h>
#include <rhea_base.h>

void
rhea_io_mpi_read_bcast_double (double *values, const size_t n_entries,
                               const char *file_path_bin,
                               const char *file_path_txt,
                               sc_MPI_Comm mpicomm)
{
  const char         *this_fn_name = "rhea_io_mpi_read_bcast_double";
  int                 mpirank, mpiret;

  if (file_path_txt == NULL) { /* if read from binary file */
    RHEA_GLOBAL_INFOF ("Into %s (%s)\n", this_fn_name, file_path_bin);
  }
  else {
    RHEA_GLOBAL_INFOF ("Into %s (bin: %s, txt: %s)\n", this_fn_name,
                       file_path_bin, file_path_txt);
  }

  /* check input */
  RHEA_ASSERT (file_path_bin != NULL || file_path_txt != NULL);

  /* get parallel environment */
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);

  /* read from file with one processor */
  if (mpirank == 0) {
    if (file_path_txt == NULL) { /* if read from binary file */
      rhea_io_std_read_double (values, n_entries, file_path_bin);
    }
    else { /* otherwise read from text file */
      rhea_io_std_read_double_from_txt (values, n_entries, file_path_txt);
      if (file_path_bin != NULL) {
        rhea_io_std_write_double (file_path_bin, values, n_entries);
      }
    }
  }

  /* broadcast values to all processors */
  mpiret = sc_MPI_Bcast (values, (int) n_entries, sc_MPI_DOUBLE, 0, mpicomm);
  SC_CHECK_MPI (mpiret);

  if (file_path_txt == NULL) { /* if read from binary file */
    RHEA_GLOBAL_INFOF ("Done %s (%s)\n", this_fn_name, file_path_bin);
  }
  else {
    RHEA_GLOBAL_INFOF ("Done %s (bin: %s, txt: %s)\n", this_fn_name,
                       file_path_bin, file_path_txt);
  }
}
