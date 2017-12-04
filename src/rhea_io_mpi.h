/*
 */

#ifndef RHEA_IO_MPI_H
#define RHEA_IO_MPI_H

#include <sc.h>

/**
 * Reads double values from a binary or text file with one processor (i.e., MPI
 * rank), then broadcasts to all processors.
 */
int                 rhea_io_mpi_read_broadcast_double (
                                                    double *values_all,
                                                    int n_entries,
                                                    const char *file_path_bin,
                                                    const char *file_path_txt,
                                                    sc_MPI_Comm mpicomm);

/**
 * Reads int values from a binary or text file with one processor (i.e., MPI
 * rank), then broadcasts to all processors.
 */
int                 rhea_io_mpi_read_broadcast_int (int *values_all,
                                                    int n_entries,
                                                    const char *file_path_bin,
                                                    const char *file_path_txt,
                                                    sc_MPI_Comm mpicomm);

/**
 * Reads double values from a binary or text file with one processor (i.e., MPI
 * rank), then scatters corresponding segments to all processors.
 */
void                rhea_io_mpi_read_scatter_double (double *values_segment,
                                                     const int *segment_offset,
                                                     const char *file_path_bin,
                                                     const char *file_path_txt,
                                                     sc_MPI_Comm mpicomm);

/**
 * Gathers double values from all processors, then writes to a text file with
 * one processor (i.e., MPI rank).
 */
void                rhea_io_mpi_gather_write_double_to_txt (
                                                  const char *file_path_txt,
                                                  const double *values_segment,
                                                  const int *segment_offset,
                                                  int n_entries_per_line,
                                                  sc_MPI_Comm mpicomm);

#endif /* RHEA_IO_MPI_H */
