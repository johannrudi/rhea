/*
 */

#include <rhea_io_mpi.h>
#include <rhea_io_std.h>
#include <rhea_base.h>

void
rhea_io_mpi_read_bcast_double (double *values, const size_t n_entries,
                               const char *file_path, sc_MPI_Comm mpicomm)
{
  const char         *this_fn_name = "rhea_io_mpi_read_bcast_double";
  int                 mpirank, mpiret;

  RHEA_GLOBAL_INFOF ("Into %s (%s)\n", this_fn_name, file_path);

  /* get parallel environment */
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);

  /* read from file with one processor */
  if (mpirank == 0) {
    rhea_io_std_read_double (values, n_entries, file_path);
  }

  /* broadcast values to all processors */
  mpiret = sc_MPI_Bcast (values, (int) n_entries, sc_MPI_DOUBLE, 0, mpicomm);
  SC_CHECK_MPI (mpiret);

  RHEA_GLOBAL_INFOF ("Done %s (%s)\n", this_fn_name, file_path);
}
