/*
  This file is part of the ymir Library.
  ymir is a C library for modeling ice sheets

  Copyright (C) 2010, 2011 Carsten Burstedde, Toby Isaac, Georg Stadler,
                           Lucas Wilcox.

  The ymir Library is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The ymir Library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the ymir Library.  If not, see <http://www.gnu.org/licenses/>.

  ---

  This is the slabs example for global instantaneous mantle flow with plates.

*/

#ifndef SLABS_IO_H
#define SLABS_IO_H

#include <ymir_mesh.h>
#include <slabs_base.h>

/* switch for using MPI I/O functions */
//#define SLABS_IO_USE_MPI_IO

/* types of coordinates */
typedef enum
{
  SL_CARTESIAN_COORDINATE,
  SL_SPHERICAL_COORDINATE_MATH_CONV,
  SL_SPHERICAL_COORDINATE_GEO_CONV
}
slabs_io_coordinate_type_t;

/**
 *
 */
void
slabs_io_write_node_coordinates_to_textfile (const char *textfile,
                                             ymir_mesh_t *mesh,
                                             slabs_node_type_t node_type,
                                             slabs_io_coordinate_type_t
                                               coord_type);

/**
 *
 */
void
slabs_io_write_cvec_to_textfile (const char *textfile, ymir_vec_t *vec,
                                 slabs_io_coordinate_type_t coord_type);

/**
 *
 */
void
slabs_io_write_dvec_to_textfile (const char *textfile, ymir_vec_t *vec,
                                 slabs_io_coordinate_type_t coord_type);

/**
 *
 */
void
slabs_io_write_double_vec_to_textfile (const char *textfile, const double *vec,
                                       const int n_values, const int n_fields,
                                       const int write_id);

/**
 * Reads, then scatters a vector.
 */
void
slabs_io_read_scatter_vector (const char *binaryfile, double *vec,
                              const ymir_gloidx_t *segment_offset,
                              const ymir_locidx_t *segment_size,
                              MPI_Comm mpicomm);

/**
 * Gathers, then writes a vector.
 */
void
slabs_io_gather_write_vector (const char *binaryfile, const double *vec,
                              const ymir_gloidx_t *segment_offset,
                              const ymir_locidx_t *segment_size,
                              MPI_Comm mpicomm);

/**
 *
 */
void
slabs_io_read_temperature (sc_dmatrix_t *temperature,
                           mangll_t *mangll, mangll_cnodes_t *cnodes,
#ifndef SLABS_IO_USE_MPI_IO
                           const char *filename_txt, const char *filename_bin
#else
                           char *filename_txt, char *filename_bin
#endif
                           );

/**
 *
 */
double *
slabs_io_read_weakzone (MPI_Comm mpicomm,
#ifndef SLABS_IO_USE_MPI_IO
                        const char *filename_txt, const char *filename_bin,
#else
                        char *filename_txt, char *filename_bin,
#endif
                        int *n_points_txt, const int n_points_bin);

/**
 * Initializes performance counters.
 */
void
slabs_io_perf_counter_init (const int active);

/**
 * Gathers global statistics of performance counters.
 */
void
slabs_io_perf_counter_gather (MPI_Comm mpicomm, const int print_wtime,
                              const int print_n_calls, const int print_flops);

/**
 * Prints statistics of performance counters.
 */
void
slabs_io_perf_counter_print ();

#endif /* SLABS_IO_H */
