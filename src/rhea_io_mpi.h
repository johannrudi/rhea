/* RHEA_IO_MPI  Reads & writes files that are distributed in parallel. */

#ifndef RHEA_IO_MPI_H
#define RHEA_IO_MPI_H

#include <sc.h>
#include <ymir_options.h>

/******************************************************************************
 * Options & Monitoring
 *****************************************************************************/

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_io_mpi_add_options (ymir_options_t * opt_sup);

/**
 * Initializes performance counters.
 */
void                rhea_io_mpi_perfmon_init (const int activate,
                                              const int skip_if_active);

/**
 * Prints statistics collected by performance monitors.
 */
void                rhea_io_mpi_perfmon_print (sc_MPI_Comm mpicomm,
                                               const int print_wtime,
                                               const int print_n_calls,
                                               const int print_flops);

/******************************************************************************
 * Read
 *****************************************************************************/

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
 * Reads float values from a binary or text file with one processor (i.e., MPI
 * rank), then broadcasts to all processors.
 */
int                 rhea_io_mpi_read_broadcast_float (
                                                    float *values_all,
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
 * Reads segments of double values from a binary file collectively wity all
 * processors.
 */
int                 rhea_io_mpi_read_segment_double (double *values_segment,
                                                     MPI_Offset segment_offset,
                                                     int segment_size,
                                                     char *file_path_bin,
                                                     sc_MPI_Comm mpicomm);

/******************************************************************************
 * Write
 *****************************************************************************/

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

/**
 * Writes segments of double values into a binary file collectively wity all
 * processors.
 */
int                 rhea_io_mpi_write_segment_double (char *file_path_bin,
                                                      double *values_segment,
                                                      MPI_Offset segment_offset,
                                                      int segment_size,
                                                      sc_MPI_Comm mpicomm);

#endif /* RHEA_IO_MPI_H */
