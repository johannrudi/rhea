/*
 */

#ifndef RHEA_IO_STD_H
#define RHEA_IO_STD_H

#include <sc.h>

/******************************************************************************
 * Binary Files
 *****************************************************************************/

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

/******************************************************************************
 * Text Files
 *****************************************************************************/

/**
 * Reads double values from a text file.
 */
size_t              rhea_io_std_read_double_from_txt (double *values,
                                                      const size_t n_entries,
                                                      const char *file_path);

/**
 * Writes double values to a text file.
 */
size_t              rhea_io_std_write_double_to_txt (
                                          const char *file_path,
                                          const char *value_whitespace_format,
                                          const int n_entries_per_line,
                                          const double *values,
                                          const size_t n_entries);

#endif /* RHEA_IO_STD_H */
