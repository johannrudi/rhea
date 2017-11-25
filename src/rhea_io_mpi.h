/*
 */

#ifndef RHEA_IO_MPI_H
#define RHEA_IO_MPI_H

#include <sc.h>

/**
 * Reads double values from a binary file with one processor (MPI rank), then
 * broadcasts to all processors.
 */
void                rhea_io_mpi_read_bcast_double (double *values,
                                                   const size_t n_entries,
                                                   const char *file_path_bin,
                                                   const char *file_path_txt,
                                                   sc_MPI_Comm mpicomm);

#endif /* RHEA_IO_MPI_H */
