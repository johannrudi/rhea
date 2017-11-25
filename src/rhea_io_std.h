/*
 */

#ifndef RHEA_IO_STD_H
#define RHEA_IO_STD_H

#include <sc.h>

/**
 * Reads double values from a binary file.
 */
size_t              rhea_io_std_read_double (double *values,
                                             const size_t n_entries,
                                             const char *file_path);

/**
 * Writes double values to a binary file.
 */
size_t              rhea_io_std_write_double (const char *file_path,
                                              const double *values,
                                              const size_t n_entries);

#endif /* RHEA_IO_STD_H */
