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
 * Issues fscanf command with appropriate format string and stores items into
 * `data`.
 *
 * \param [out] data      Data that is read from file
 * \param [in]  params    User parameters
 * \param [in]  file_ptr  File pointer
 * \param [in]  n_read    Number of items read so far
 * \return                Status returned from fscanf
 */
typedef int       (*rhea_io_std_fscanf_fn_t) (void *data, void *params,
                                              FILE *file_ptr,
                                              const size_t n_read);

/**
 * Issues fprintf command with appropriate format string and stores items from
 * `data`.
 *
 * \param [in] file_ptr   File pointer
 * \param [in] data       Data that is read from file
 * \param [in] params     User parameters
 * \param [in] n_written  Number of items read so far
 * \return                Status returned from fprintf
 */
typedef int       (*rhea_io_std_fprintf_fn_t) (FILE *file_ptr, const void *data,
                                               void *params,
                                               const size_t n_written);

/**
 * Reads generic data from a text file.
 */
size_t              rhea_io_std_read_txt (void *data,
                                          rhea_io_std_fscanf_fn_t fscanf_fn,
                                          void *fscanf_params,
                                          const char *file_path);

/**
 * Writes generic data to a text file.
 */
size_t              rhea_io_std_write_txt (const char *file_path,
                                           const void *data,
                                           rhea_io_std_fprintf_fn_t fprintf_fn,
                                           void *fprintf_params);

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
                                          const double *values,
                                          const size_t n_entries,
                                          int n_entries_per_line);

#endif /* RHEA_IO_STD_H */
