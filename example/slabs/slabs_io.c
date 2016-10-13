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

#include <slabs_io.h>
#include <slabs_physics.h>
#include <ymir_perf_counter.h>

/* number to string conversion for writing node coordinates and vector values */
#define SL_IO_WRITE_COORD_NODE_ID "%12lli"
#define SL_IO_WRITE_COORD_PRECISION "%+1.15e"
#define SL_IO_WRITE_VAL_PRECISION "%+1.15e"

/* number of characters for temperature import from text file */
#define SL_IO_TEMP_IMPORT_N_CHARS_INDEX 12
#define SL_IO_TEMP_IMPORT_N_CHARS_VALUE 22

/* parameters for parallel I/O */
#define SLABS_IO_MAX_SIZE_1_PROC (5000000)
#define SLABS_IO_NUM_PROCS_PER_IO_GROUP (512)

/* performance counters */
typedef enum
{
  SLABS_IO_COUNTER_TEMPERATURE,
  SLABS_IO_COUNTER_WEAKZONE,
  SLABS_IO_COUNTER_READ_SCATTER_VECTOR,
  SLABS_IO_COUNTER_GATHER_WRITE_VECTOR,
  SLABS_IO_COUNTER_N
}
slabs_io_counter_idx_t;
ymir_perf_counter_t slabs_io_counter[SLABS_IO_COUNTER_N];
const char         *slabs_io_counter_name[SLABS_IO_COUNTER_N] =
{
  "IO temperature",
  "IO weakzone",
  "IO read & scatter vector",
  "IO gather & write vector"
};
sc_statinfo_t       slabs_io_stats[
                      SLABS_IO_COUNTER_N * YMIR_PERF_COUNTER_N_STATS];
char                slabs_io_stats_name[
                      SLABS_IO_COUNTER_N * YMIR_PERF_COUNTER_N_STATS][
                      YMIR_PERF_COUNTER_NAME_SIZE];
int                 slabs_io_n_stats;

/* helper structure for writing vectors to file */
typedef struct slabs_io_write_vec_to_txt
{
  FILE               *fp;
  mangll_cnodes_t    *cnodes;
  mangll_gloidx_t     dnode_offset;
  int                 n_fields;
  slabs_io_coordinate_type_t  coord_type;
  int                 write_coord_only;
}
slabs_io_write_vec_to_txt_t;

/**
 * Calculates the number of digits for MPI ranks, given an MPI size.
 */
static int
slabs_io_set_num_digits_for_mpiranks (const int mpisize)
{
  int                 n_digits;

  if (mpisize <= 10000) {
    n_digits = 4;
  }
  else if (mpisize <= 1000000) {
    n_digits = 6;
  }
  else {
    n_digits = 9; /* hopefully this is enough */
    YMIR_ASSERT (mpisize <= 1000000000);
  }

  return n_digits;
}

/*
 * Opens a file.
 */
static FILE *
slabs_io_file_open (const char *filename, const char *access_type,
                    const char *this_fn_name)
{
  FILE               *fp;

  /* open file */
  fp = fopen (filename, access_type);

  /* check status */
  if (!fp) {
    if (this_fn_name != NULL) {
      YMIR_ABORTF ("%s: Could not open file \"%s\"\n", this_fn_name, filename);
    }
    else {
      YMIR_ABORTF ("Could not open file \"%s\"\n", filename);
    }
  }

  /* return file pointer */
  return fp;
}

/*
 * Closes a file.
 */
static void
slabs_io_file_close (FILE *fp, const char *filename, const char *this_fn_name)
{
  int                 status;

  /* close file */
  status = fclose (fp);

  /* check status */
  if (status != 0) {
    if (this_fn_name != NULL && filename != NULL) {
      YMIR_LERRORF ("%s: Error closing file \"%s\"\n", this_fn_name, filename);
    }
    else if (filename != NULL) {
      YMIR_LERRORF ("Error closing file \"%s\"\n", filename);
    }
    else {
      YMIR_LERROR ("Error closing file\n");
    }
  }
}

/**
 *
 */
static void
slabs_io_write_vector_to_textfile_fn (double *val, double x, double y, double z,
                                      ymir_locidx_t nid, void *data)
{
  slabs_io_write_vec_to_txt_t  *params = (slabs_io_write_vec_to_txt_t *) data;
  long long int       node_gloidx;
  double              coord1, coord2, coord3;
  int                 fieldid;

  /* set global node index */
  if (params->cnodes != NULL) { /* if continuous node */
    if (nid < params->cnodes->owncount) { /* if node is owned by MPI process */
      node_gloidx = (long long int) mangll_cnodes_global_index (params->cnodes,
                                                                nid);
    }
    else { /* if node is not owned by MPI process */
      node_gloidx = -1;
    }
  }
  else { /* if discontinuous node */
    node_gloidx = (long long int) params->dnode_offset + nid;
  }

  /* write node coordinates and vector value at this node */
  if (0 <= node_gloidx) { /* if node should be written to file */
    /* set coordinates in the right format */
    switch (params->coord_type) {
    case SL_CARTESIAN_COORDINATE:
      coord1 = x;
      coord2 = y;
      coord3 = z;
      break;

    case SL_SPHERICAL_COORDINATE_MATH_CONV:
      slabs_convert_cartesian_to_spherical_math_conv (x, y, z, &coord1,
                                                      &coord2, &coord3);
      break;

    case SL_SPHERICAL_COORDINATE_GEO_CONV:
      slabs_convert_cartesian_to_spherical_geo_conv (x, y, z, &coord1,
                                                     &coord2, &coord3);
      break;

    default: /* unknown coordinate type */
      YMIR_ABORT_NOT_REACHED ();
    }

    /* write node coordinates */
    fprintf (params->fp,
             SL_IO_WRITE_COORD_NODE_ID " " SL_IO_WRITE_COORD_PRECISION " "
             SL_IO_WRITE_COORD_PRECISION " " SL_IO_WRITE_COORD_PRECISION,
             node_gloidx, coord1, coord2, coord3);

    if (params->write_coord_only) {
      /* write only node coordinates */
      fprintf (params->fp, "\n");
    }
    else {
      /* write values of the vector */
      for (fieldid = 0; fieldid < params->n_fields; fieldid++) {
        fprintf (params->fp, " " SL_IO_WRITE_VAL_PRECISION, val[fieldid]);
      }
      fprintf (params->fp, "\n");
    }
  }
}

/**
 *
 */
void
slabs_io_write_node_coordinates_to_textfile (const char *textfile,
                                             ymir_mesh_t *mesh,
                                             slabs_node_type_t node_type,
                                             slabs_io_coordinate_type_t
                                               coord_type)
{
  const char         *this_fn_name =
                        "slabs_io_write_node_coordinates_to_textfile";
  const int           mpisize = mesh->ma->mpisize;
  const int           mpirank = mesh->ma->mpirank;
  const int           N = ymir_n (mesh->ma->N);
  const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);
  int                 empty = 1;
  ymir_vec_t         *dummyvec;
  char                filename[BUFSIZ];
  FILE               *fp;
  slabs_io_write_vec_to_txt_t  params;

  YMIR_GLOBAL_INFOF ("Into %s: Write node coordinates to file \"%s*.txt\"\n",
                     this_fn_name, textfile);

  /* create a dummy vector with the desired node type */
  switch (node_type) {
  case SL_GLL_CONTINUOUS_NODE:
    if (0 < mesh->cnodes->owncount) {
      empty = 0;
      dummyvec = ymir_cvec_new (mesh, 1);
      params.cnodes = mesh->cnodes;
    }
    break;

  case SL_GLL_DISCONTINUOUS_NODE:
    if (0 < mesh->ma->mesh->K) {
      empty = 0;
      dummyvec = ymir_dvec_new (mesh, 1, YMIR_GLL_NODE);
      params.cnodes = NULL;
      params.dnode_offset = mesh->ma->mesh->RtoGEO[mpirank] * n_nodes_per_el;
    }
    break;

  case SL_GAUSS_NODE:
    if (0 < mesh->ma->mesh->K) {
      empty = 0;
      dummyvec = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
      params.cnodes = NULL;
      params.dnode_offset = mesh->ma->mesh->RtoGEO[mpirank] * n_nodes_per_el;
    }
    break;

  default: /* unknown node type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* do nothing if processor is empty */
  if (empty) {
    YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
    return;
  }

  /* append MPI rank to filename */
  if (1 < mpisize) {
    int                 n_digits;

    /* set number of digits for mpiranks */
    n_digits = slabs_io_set_num_digits_for_mpiranks (mpisize);

    snprintf (filename, BUFSIZ, "%s.%0*i.txt", textfile, n_digits, mpirank);
  }
  else {
    snprintf (filename, BUFSIZ, "%s.txt", textfile);
  }

  /* open file */
  fp = slabs_io_file_open (filename, "w", this_fn_name);

  /* write coordinates */
  params.fp = fp;
  params.coord_type = coord_type;
  params.write_coord_only = 1;
  switch (node_type) {
  case SL_GLL_CONTINUOUS_NODE:
    ymir_cvec_set_function (dummyvec, slabs_io_write_vector_to_textfile_fn,
                            &params);
    break;

  case SL_GLL_DISCONTINUOUS_NODE:
  case SL_GAUSS_NODE:
    ymir_dvec_set_function (dummyvec, slabs_io_write_vector_to_textfile_fn,
                            &params);
    break;

  default: /* unknown node type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* close file */
  slabs_io_file_close (fp, filename, this_fn_name);

  /* destroy */
  ymir_vec_destroy (dummyvec);

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

/**
 *
 */
void
slabs_io_write_cvec_to_textfile (const char *textfile, ymir_vec_t *vec,
                                 slabs_io_coordinate_type_t coord_type)
{
  const char         *this_fn_name = "slabs_io_write_cvec_to_textfile";
  ymir_mesh_t        *mesh = vec->mesh;
  const int           mpisize = mesh->ma->mpisize;
  const int           mpirank = mesh->ma->mpirank;
  char                filename[BUFSIZ];
  FILE               *fp;
  slabs_io_write_vec_to_txt_t  params;

  YMIR_GLOBAL_INFOF ("Into %s: Write continuous node vector to "
                     "file \"%s*.txt\"\n", this_fn_name, textfile);

  /* check input parameters */
  YMIR_ASSERT_HAS_CVEC (vec);

  /* do nothing if processor is empty */
  if (vec->owncount <= 0) {
    YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
    return;
  }

  /* append MPI rank to filename */
  if (1 < mpisize) {
    int                 n_digits;

    /* set number of digits for mpiranks */
    n_digits = slabs_io_set_num_digits_for_mpiranks (mpisize);

    snprintf (filename, BUFSIZ, "%s.%0*i.txt", textfile, n_digits, mpirank);
  }
  else {
    snprintf (filename, BUFSIZ, "%s.txt", textfile);
  }

  /* open file */
  fp = slabs_io_file_open (filename, "w", this_fn_name);

  /* write vector coordinates and values */
  params.fp = fp;
  params.cnodes = mesh->cnodes;
  params.n_fields = vec->ncfields;
  params.coord_type = coord_type;
  params.write_coord_only = 0;
  ymir_cvec_set_function (vec, slabs_io_write_vector_to_textfile_fn, &params);

  /* close file */
  slabs_io_file_close (fp, filename, this_fn_name);

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

/**
 *
 */
void
slabs_io_write_dvec_to_textfile (const char *textfile, ymir_vec_t *vec,
                                 slabs_io_coordinate_type_t coord_type)
{
  const char         *this_fn_name = "slabs_io_write_dvec_to_textfile";
  ymir_mesh_t        *mesh = vec->mesh;
  const int           mpisize = mesh->ma->mpisize;
  const int           mpirank = mesh->ma->mpirank;
  const int           N = ymir_n (mesh->ma->N);
  const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);
  const mangll_gloidx_t  dnode_offset = mesh->ma->mesh->RtoGEO[mpirank] *
                                        n_nodes_per_el;
  char                filename[BUFSIZ];
  FILE               *fp;
  slabs_io_write_vec_to_txt_t  params;

  YMIR_GLOBAL_INFOF ("Into %s: Write discontinuous node vector to "
                     "file \"%s*.txt\"\n", this_fn_name, textfile);

  /* check input parameters */
  YMIR_ASSERT_HAS_DVEC (vec);
/* do nothing if processor is empty */
  if (vec->K <= 0) {
    YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
    return;
  }

  /* append MPI rank to filename */
  if (1 < mpisize) {
    int                 n_digits;

    /* set number of digits for mpiranks */
    n_digits = slabs_io_set_num_digits_for_mpiranks (mpisize);

    snprintf (filename, BUFSIZ, "%s.%0*i.txt", textfile, n_digits, mpirank);
  }
  else {
    snprintf (filename, BUFSIZ, "%s.txt", textfile);
  }

  /* open file */
  fp = slabs_io_file_open (filename, "w", this_fn_name);

  /* write vector coordinates and values */
  params.fp = fp;
  params.cnodes = NULL;
  params.dnode_offset = dnode_offset;
  params.n_fields = vec->ndfields;
  params.coord_type = coord_type;
  params.write_coord_only = 0;
  ymir_dvec_set_function (vec, slabs_io_write_vector_to_textfile_fn, &params);

  /* close file */
  slabs_io_file_close (fp, filename, this_fn_name);

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

/**
 *
 */
void
slabs_io_write_double_vec_to_textfile (const char *textfile, const double *vec,
                                       const int n_values, const int n_fields,
                                       const int write_id)
{
  const char         *this_fn_name = "slabs_io_write_double_vec_to_textfile";
  char                filename[BUFSIZ];
  FILE               *fp;
  int                 i, j;

  YMIR_INFOF ("Into %s: Write double vector to file \"%s.txt\"\n",
              this_fn_name, textfile);

  /* set filename */
  snprintf (filename, BUFSIZ, "%s.txt", textfile);

  /* open file */
  fp = slabs_io_file_open (filename, "w", this_fn_name);

  /* write values */
  for (i = 0; i < n_values; i++) {
    if (write_id) {
      fprintf (fp, "%12i ", i);
    }
    for (j = 0; j < n_fields; j++) {
      fprintf (fp, SL_IO_WRITE_VAL_PRECISION " ", vec[n_fields*i + j]);
    }
    fprintf (fp, "\n");
  }

  /* close file */
  slabs_io_file_close (fp, filename, this_fn_name);

  YMIR_INFOF ("Done %s\n", this_fn_name);
}

#if !defined(SLABS_IO_USE_MPI_IO)

/**
 * Reads double vector from binary file serially.
 */
static int
slabs_io_read_vector (const char *binaryfile, double *vec,
                      const size_t size)
{
  const char         *this_fn_name = "slabs_io_read_vector";
  FILE               *fp;
  size_t              size_read;

  YMIR_INFOF ("Into %s: Read vector from binary file \"%s\"\n", this_fn_name,
              binaryfile);

  /* open file */
  fp = slabs_io_file_open (binaryfile, "rb", this_fn_name);

  /* read file */
  size_read = fread (vec, sizeof (double), size, fp);

  /* check size */
  if (size_read != size) {
    YMIR_LERRORF ("%s: Mismatch of read size for file %s. "
                  "Requested %lli, actually read %lli\n",
                  this_fn_name, binaryfile, (long long int) size,
                  (long long int) size_read);
  }

  /* close file */
  slabs_io_file_close (fp, binaryfile, this_fn_name);

  YMIR_INFOF ("Done %s\n", this_fn_name);

  return size_read;
}

/**
 * Writes double vector to binary file serially.
 */
static void
slabs_io_write_vector (const char *binaryfile, const double *vec,
                       const size_t size)
{
  const char         *this_fn_name = "slabs_io_write_vector";
  FILE               *fp;
  size_t              size_written;

  YMIR_INFOF ("Into %s: Write vector to binary file \"%s\"\n", this_fn_name,
              binaryfile);

  /* open file */
  fp = slabs_io_file_open (binaryfile, "wb", this_fn_name);

  /* write file */
  size_written = fwrite (vec, sizeof (double), size, fp);

  /* check size */
  if (size_written != size) {
    YMIR_LERRORF ("%s: Mismatch of write size for file %s. "
                  "Requested %lli, actually written %lli\n",
                  this_fn_name, binaryfile, (long long int) size,
                  (long long int) size_written);
  }

  /* close file */
  slabs_io_file_close (fp, binaryfile, this_fn_name);

  YMIR_INFOF ("Done %s\n", this_fn_name);
}

/**
 * Reads, then broadcasts a vector.
 */
static void
slabs_io_read_bcast_vector_short (const char *binaryfile, double *vec,
                                  const int size, MPI_Comm mpicomm)
{
  const char         *this_fn_name = "slabs_io_read_bcast_vector_short";
  int                 mpirank, mpiret;

  YMIR_GLOBAL_INFOF ("Into %s: Read & bcast vector from binary file \"%s\"\n",
                     this_fn_name, binaryfile);

  /* get parallel environment */
  mpiret = MPI_Comm_rank (mpicomm, &mpirank); YMIR_CHECK_MPI (mpiret);

  /* read file with one processor */
  if (mpirank == 0) {
    slabs_io_read_vector (binaryfile, vec, (size_t) size);
  }

  /* broadcast data */
  mpiret = MPI_Bcast (vec, size, MPI_DOUBLE, 0, mpicomm);
  YMIR_CHECK_MPI (mpiret);

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

/**
 * Reads, then scatters a vector.
 */
static void
slabs_io_read_scatter_vector_short (const char *binaryfile, double *vec,
                                    const int *segment_offset,
                                    const int *segment_size,
                                    MPI_Comm mpicomm)
{
  const char         *this_fn_name = "slabs_io_read_scatter_vector_short";
  int                 mpisize, mpirank, mpiret;
  double             *vec_global;

  YMIR_GLOBAL_INFOF ("Into %s: Read & scatter vector from binary file \"%s\"\n",
                     this_fn_name, binaryfile);

  /* get parallel environment */
  mpiret = MPI_Comm_size (mpicomm, &mpisize); YMIR_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &mpirank); YMIR_CHECK_MPI (mpiret);

  /* read file with one processor */
  if (mpirank == 0) {
    vec_global = YMIR_ALLOC (double, segment_offset[mpisize]);
    slabs_io_read_vector (binaryfile, vec_global,
                          (size_t) segment_offset[mpisize]);
  }
  else {
    vec_global = NULL;
  }
  mpiret = MPI_Barrier (mpicomm); YMIR_CHECK_MPI (mpiret);

  /* scatter data */
  mpiret = MPI_Scatterv (vec_global, segment_size, segment_offset, MPI_DOUBLE,
                         vec, segment_size[mpirank], MPI_DOUBLE, 0, mpicomm);
  YMIR_CHECK_MPI (mpiret);

  /* destroy */
  if (mpirank == 0) {
    YMIR_FREE (vec_global);
  }

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

void
slabs_io_read_scatter_vector (const char *binaryfile, double *vec,
                              const ymir_gloidx_t *segment_offset,
                              const ymir_locidx_t *segment_size,
                              MPI_Comm mpicomm)
{
  const char         *this_fn_name = "slabs_io_read_scatter_vector";
  int                 mpisize, mpirank, mpiret;
  MPI_Comm            mpicomm_io, mpicomm_scatter;
  int                 n_groups, group_id;
  int                 group_size, group_rank_first, group_rank_last;
  int                 remaining_procs;
  int                *group_segment_offset, *group_segment_size;
  int                 vec_group_size;
  double             *vec_group;
  int                 this_proc_size;
  int                 r;

  /* start performance counters */
  ymir_perf_counter_start (
      &slabs_io_counter[SLABS_IO_COUNTER_READ_SCATTER_VECTOR]);

  /* get parallel environment */
  mpiret = MPI_Comm_size (mpicomm, &mpisize); YMIR_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &mpirank); YMIR_CHECK_MPI (mpiret);

  /* use simple version for short vectors */
  if (segment_offset[mpisize] < SLABS_IO_MAX_SIZE_1_PROC) {
    int                *offset = YMIR_ALLOC (int, mpisize + 1);
    int                *size = YMIR_ALLOC (int, mpisize);

    for (r = 0; r < mpisize; r++) {
      offset[r] = (int) segment_offset[r];
      size[r] = (int) segment_size[r];
    }
    offset[mpisize] = (int) segment_offset[mpisize];

    slabs_io_read_scatter_vector_short (binaryfile, vec, offset, size, mpicomm);

    YMIR_FREE (offset);
    YMIR_FREE (size);

    /* stop performance counters */
    ymir_perf_counter_stop_add (
        &slabs_io_counter[SLABS_IO_COUNTER_READ_SCATTER_VECTOR]);

    return;
  }

  YMIR_GLOBAL_INFOF ("Into %s: Read & scatter vector from binary file \"%s\"\n",
                     this_fn_name, binaryfile);

  /* determine groups for IO */
  group_size = SC_MIN (SLABS_IO_NUM_PROCS_PER_IO_GROUP, mpisize);
  n_groups = mpisize / group_size;
  group_id = mpirank / group_size;
  remaining_procs = mpisize - n_groups * group_size;
  /* apply correction for remaining processors */
  group_rank_first = group_id * group_size;
  if (0 < remaining_procs && (n_groups - 1) <= group_id) {
    /* add remaining ranks to last group, if mpisize/group_size has remainder */
    group_id = n_groups - 1;
    group_rank_first = group_id * group_size;
    group_size += remaining_procs;
  }
  group_rank_last = group_rank_first + group_size - 1;
  YMIR_ASSERT (0 <= group_rank_first && group_rank_last < mpisize);

  YMIR_GLOBAL_INFOF ("%s: Use %i groups (~%i ranks per group) to read "
                     "vector of global size %lli\n",
                     this_fn_name, n_groups, group_size,
                     (long long int) segment_offset[mpisize]);

  /* set offsets and sizes for this group */
  group_segment_offset = YMIR_ALLOC (int, group_size + 1);
  group_segment_size = YMIR_ALLOC (int, group_size);
  group_segment_offset[0] = 0;
  for (r = 0; r < group_size; r++) {
    group_segment_offset[r+1] = (int) (segment_offset[group_rank_first+r+1] -
                                       segment_offset[group_rank_first]);
    group_segment_size[r] = group_segment_offset[r+1] - group_segment_offset[r];
    YMIR_ASSERT (group_segment_size[r] == segment_size[group_rank_first+r]);
  }

  /* allocate storage for vector segment of this group */
  if (mpirank == group_rank_first) {
    vec_group_size = group_segment_offset[group_size];
    vec_group = YMIR_ALLOC (double, vec_group_size);
  }
  else {
    vec_group_size = 0;
    vec_group = NULL;
  }

  /* create MPI communicator for IO */
  mpiret = MPI_Comm_split (mpicomm,
                           (mpirank == group_rank_first ? 0 : MPI_UNDEFINED),
                           mpirank, &mpicomm_io);
  YMIR_CHECK_MPI (mpiret);

  /* read vector from file with one processor per group */
  if (mpirank == group_rank_first) {
    MPI_File            fp;
    MPI_Offset          mpi_offset = segment_offset[mpirank] * sizeof (double);
    MPI_Status          mpistatus;
    int                 count;

    YMIR_ASSERT (mpicomm_io != MPI_COMM_NULL);

    mpiret = MPI_File_open (mpicomm_io, binaryfile, MPI_MODE_RDONLY,
                            MPI_INFO_NULL, &fp); YMIR_CHECK_MPI (mpiret);

    mpiret = MPI_File_read_at (fp, mpi_offset, vec_group, vec_group_size,
                               MPI_DOUBLE, &mpistatus); YMIR_CHECK_MPI (mpiret);
    mpiret = MPI_Get_count (&mpistatus, MPI_DOUBLE, &count);
    YMIR_CHECK_MPI (mpiret);
    YMIR_CHECK_ABORT (count == vec_group_size, "read count");

    mpiret = MPI_File_close (&fp); YMIR_CHECK_MPI (mpiret);
  }

  /* create MPI communicator for scattering */
  mpiret = MPI_Comm_split (mpicomm, group_id, mpirank, &mpicomm_scatter);
  YMIR_CHECK_MPI (mpiret);

  /* scatter vector segments within groups */
  this_proc_size = (int) segment_size[mpirank];
  mpiret = MPI_Scatterv (vec_group, group_segment_size, group_segment_offset,
                         MPI_DOUBLE, vec, this_proc_size, MPI_DOUBLE,
                         0, mpicomm_scatter); YMIR_CHECK_MPI (mpiret);

  /* destroy */
  YMIR_FREE (group_segment_offset);
  YMIR_FREE (group_segment_size);
  if (mpirank == group_rank_first) {
    YMIR_FREE (vec_group);
  }
  mpiret = MPI_Comm_free (&mpicomm_scatter); YMIR_CHECK_MPI (mpiret);
  if (mpicomm_io != MPI_COMM_NULL) {
    mpiret = MPI_Comm_free (&mpicomm_io); YMIR_CHECK_MPI (mpiret);
  }

  /* stop performance counters */
  ymir_perf_counter_stop_add (
      &slabs_io_counter[SLABS_IO_COUNTER_READ_SCATTER_VECTOR]);

  YMIR_GLOBAL_INFOF (
      "Done %s---time %g sec\n", this_fn_name,
      slabs_io_counter[SLABS_IO_COUNTER_READ_SCATTER_VECTOR].wtime_intv);
}

/**
 * Gathers, then writes a vector.
 */
static void
slabs_io_gather_write_vector_short (const char *binaryfile, const double *vec,
                                    const int *segment_offset,
                                    const int *segment_size,
                                    MPI_Comm mpicomm)
{
  const char         *this_fn_name = "slabs_io_gather_write_vector_short";
  int                 mpisize, mpirank, mpiret;
  double             *vec_global;

  YMIR_GLOBAL_INFOF ("Into %s: Gather & write vector to binary file \"%s\"\n",
                     this_fn_name, binaryfile);

  /* get parallel environment */
  mpiret = MPI_Comm_size (mpicomm, &mpisize); YMIR_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &mpirank); YMIR_CHECK_MPI (mpiret);

  /* gather data */
  if (mpirank == 0) {
    vec_global = YMIR_ALLOC (double, segment_offset[mpisize]);
  }
  else {
    vec_global = NULL;
  }
  mpiret = MPI_Gatherv (vec, segment_size[mpirank], MPI_DOUBLE,
                        vec_global, segment_size, segment_offset, MPI_DOUBLE,
                        0, mpicomm); YMIR_CHECK_MPI (mpiret);

  /* write file with one processor */
  if (mpirank == 0) {
    slabs_io_write_vector (binaryfile, vec_global,
                           (size_t) segment_offset[mpisize]);
    YMIR_FREE (vec_global);
  }

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

/* TODO this function is a quick hack to get io working on BGQ,
 * see original below */
void
slabs_io_gather_write_vector (const char *binaryfile, const double *vec,
                              const ymir_gloidx_t *segment_offset,
                              const ymir_locidx_t *segment_size,
                              MPI_Comm mpicomm)
{
  const char         *this_fn_name = "slabs_io_gather_write_vector";
  int                 mpisize, mpirank;
  int                 n_groups, group_id;
  int                 group_size;
  int                 i;

  /* start performance counters */
  ymir_perf_counter_start (
      &slabs_io_counter[SLABS_IO_COUNTER_GATHER_WRITE_VECTOR]);

  YMIR_GLOBAL_INFOF ("Into %s: Gather & write vector to binary file \"%s\"\n",
                     this_fn_name, binaryfile);

  /* get parallel environment */
  MPI_Comm_size (mpicomm, &mpisize);
  MPI_Comm_rank (mpicomm, &mpirank);

  /* set segment offset and size for this processor */
  long int this_proc_offset = segment_offset[mpirank];
  int      this_proc_size = segment_size[mpirank];

  /* determine groups for IO */
  group_size = SC_MIN (SLABS_IO_NUM_PROCS_PER_IO_GROUP, mpisize);
  n_groups = mpisize / group_size;
  group_id = mpirank / group_size;
  if (mpisize - n_groups * group_size > 0)
  {
      printf("ON BG/Q NO REMAINING PROCESS SHOULD EXIST!!\n");
      fflush(stdout);
      MPI_Abort(MPI_COMM_WORLD, -1);
  }

  YMIR_GLOBAL_INFOF ("%s: Use %i groups (~%i ranks per group) to write "
                     "vector of global size %lli - NEW IMPLEMENTATION\n",
                     this_fn_name, n_groups, group_size,
                     (long long int) segment_offset[mpisize]);
  YMIR_GLOBAL_INFO ("THIS IS THE NEW IMPLEMENTATION!!!\n");

  MPI_Comm scomm;
  MPI_Comm_split(mpicomm, group_id, mpirank, &scomm );
  if (scomm == MPI_COMM_NULL)
     printf("RANK %d, SCOM IS NULL! \n", mpirank);

  int local_rank = -1;
  MPI_Comm_rank(scomm, &local_rank);

  double coltime = MPI_Wtime();
  int *local_size = NULL;
  local_size = (int*) calloc(group_size, sizeof(int));
  local_size[local_rank] = this_proc_size;

  //TODO fix: this is known already, no communication required
  if (local_rank == 0)
     MPI_Reduce( MPI_IN_PLACE, local_size, group_size, MPI_INT, MPI_SUM, 0, scomm);
  else
     MPI_Reduce( local_size, local_size, group_size, MPI_INT, MPI_SUM, 0, scomm);

  int *displs = NULL;
  double *abkxa_global_vec = NULL; //TODO terrible variable name
  int global_size = 0;
  if (local_rank == 0) {
      for(i = 0; i < group_size; i++) {
          global_size += local_size[i];
      }

      displs = (int*) calloc(group_size, sizeof(int));
      for(i = 1; i < group_size; i++) {
           displs[i] = displs[i-1] + local_size[i-1];
      }

      abkxa_global_vec = (double*) malloc(global_size * sizeof(double));
      if(abkxa_global_vec == NULL) {
        printf("ERROR cannot allocate internal vector!!!\n");
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, -1);
      }
  }
  MPI_Gatherv(vec, this_proc_size, MPI_DOUBLE, abkxa_global_vec, local_size, displs, MPI_DOUBLE, 0, scomm);
  if(local_size != NULL) free(displs);
  if(local_size != NULL) free(local_size);
  if(scomm != MPI_COMM_NULL) MPI_Comm_free(&scomm);

  coltime = MPI_Wtime() - coltime;
//if (local_rank == 0)
//    printf("Rank %d, Write comm time: %.2f\n", mpirank, coltime);
      //TODO these are a loooot of prints, not for large runs!
  int leader_color = MPI_UNDEFINED;
  if(local_rank == 0)
    leader_color = 0;

//  MPI_Barrier(mpicomm);
  MPI_Comm inter_comm;
  MPI_Comm_split(mpicomm, leader_color, mpirank, &inter_comm);
//  MPI_Barrier(mpicomm);

  double wtime = MPI_Wtime();
  if(local_rank == 0) {
      MPI_File file;
      MPI_Status status;
      MPI_File_open(inter_comm, binaryfile, MPI_MODE_CREATE|MPI_MODE_WRONLY,
                                MPI_INFO_NULL, &file);
      MPI_Offset offset = sizeof(double)*this_proc_offset;
      MPI_File_write_at(file, offset, abkxa_global_vec, global_size, MPI_DOUBLE, &status);
      MPI_File_close(&file);
  }

  wtime = MPI_Wtime() - wtime;
//if (local_rank == 0)
//    printf("Rank: %d, write file time: %.2f\n", mpirank, wtime);
      //TODO these are a loooot of prints, not for large runs!
//  MPI_Barrier(mpicomm);
  if(inter_comm != MPI_COMM_NULL) MPI_Comm_free(&inter_comm);
  if(local_rank == 0) { free(abkxa_global_vec); abkxa_global_vec = NULL; }

  /* stop performance counters */
  ymir_perf_counter_stop_add (
      &slabs_io_counter[SLABS_IO_COUNTER_GATHER_WRITE_VECTOR]);

  YMIR_GLOBAL_INFOF (
      "Done %s---time %g sec\n", this_fn_name,
      slabs_io_counter[SLABS_IO_COUNTER_GATHER_WRITE_VECTOR].wtime_intv);
}

//TODO below is the actual non-MPI-IO version, this code needs to be cleaned up
#if 0
void
slabs_io_gather_write_vector (const char *binaryfile, const double *vec,
                              const ymir_gloidx_t *segment_offset,
                              const ymir_locidx_t *segment_size,
                              MPI_Comm mpicomm)
{
  const char         *this_fn_name = "slabs_io_gather_write_vector";
  ymir_gloidx_t       this_proc_offset;
  ymir_locidx_t       this_proc_size;
  int                 mpisize, mpirank, mpiret;
  int                 n_groups, group_id;
  int                 group_size, remaining_procs;
  int                 group_rank_first, group_rank_last;
  MPI_Status          mpistatus;
  FILE               *fp;
  int                 start_io;
  int                 ret;
  size_t              size_written;

  /* start performance counters */
  ymir_perf_counter_start (
      &slabs_io_counter[SLABS_IO_COUNTER_GATHER_WRITE_VECTOR]);

  /* get parallel environment */
  mpiret = MPI_Comm_size (mpicomm, &mpisize); YMIR_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &mpirank); YMIR_CHECK_MPI (mpiret);

  /* use simple version for short vectors */
  if (segment_offset[mpisize] < SLABS_IO_MAX_SIZE_1_PROC) {
    int                *offset = YMIR_ALLOC (int, mpisize + 1);
    int                *size = YMIR_ALLOC (int, mpisize);
    int                 r;

    for (r = 0; r < mpisize; r++) {
      offset[r] = (int) segment_offset[r];
      size[r] = (int) segment_size[r];
    }
    offset[mpisize] = (int) segment_offset[mpisize];

    slabs_io_gather_write_vector_short (binaryfile, vec, offset, size, mpicomm);

    YMIR_FREE (offset);
    YMIR_FREE (size);

    /* stop performance counters */
    ymir_perf_counter_stop_add (
        &slabs_io_counter[SLABS_IO_COUNTER_GATHER_WRITE_VECTOR]);

    return;
  }

  YMIR_GLOBAL_INFOF ("Into %s: Gather & write vector to binary file \"%s\"\n",
                     this_fn_name, binaryfile);

  /* set segment offset and size for this processor */
  this_proc_offset = segment_offset[mpirank];
  this_proc_size = segment_size[mpirank];

  /* determine groups for IO */
  group_size = SLABS_IO_NUM_PROCS_PER_IO_GROUP;
  n_groups = mpisize / group_size;
  group_id = mpirank / group_size;
  remaining_procs = mpisize - n_groups * group_size;

  group_rank_first = group_id * group_size;
  if (0 < remaining_procs && (n_groups - 1) <= group_id) {
    /* add remaining ranks to last group, if mpisize/group_size has remainder */
    group_id = n_groups - 1;
    group_rank_first = group_id * group_size;
    group_size += remaining_procs;
  }
  group_rank_last = group_rank_first + group_size - 1;
  YMIR_ASSERT (0 <= group_rank_first && group_rank_last < mpisize);

  YMIR_GLOBAL_INFOF ("%s: Use %i groups (~%i ranks per group) to write "
                     "vector of global size %lli\n",
                     this_fn_name, n_groups, group_size,
                     (long long int) segment_offset[mpisize]);

  /* open/create file with root rank */
  if (0 == mpirank) {
    start_io = 1;
    fp = slabs_io_file_open (binaryfile, "wb", this_fn_name);
  }
  mpiret = MPI_Barrier (mpicomm); YMIR_CHECK_MPI (mpiret);

  /* open file in update mode on all other ranks */
  if (0 < mpirank) {
    if (group_rank_first == mpirank) { /* if rank is first in group */
      start_io = 1;
      fp = slabs_io_file_open (binaryfile, "rb+", this_fn_name);
    }
    else { /* if rank is not first in group */
      start_io = 0;

      /* wait for sequential synchronization within group */
      mpiret = MPI_Recv (&start_io, 1, MPI_INT, mpirank - 1, group_id,
                         mpicomm, &mpistatus); YMIR_CHECK_MPI (mpiret);
      fp = slabs_io_file_open (binaryfile, "rb+", this_fn_name);
    }
  }

  /* move file pointer to the beginning of this processor's segment */
  if (0 < mpirank) {
    ret = fseek (fp, this_proc_offset * sizeof (double), SEEK_SET);
    YMIR_CHECK_ABORT (ret == 0, "seek for this processor's segment");
  }

  /* write own segment of the vector */
  size_written = fwrite (vec, sizeof (double), (size_t) this_proc_size, fp);
  YMIR_CHECK_ABORT (size_written == (size_t) this_proc_size, "fwrite vector");

  /* close file */
  ret = fflush (fp); YMIR_CHECK_ABORT (ret == 0, "file flush");
  slabs_io_file_close (fp, binaryfile, this_fn_name);
  fp = NULL;

  /* initiate sequential synchronization within group */
  if (mpirank < group_rank_last) {
    mpiret = MPI_Send (&start_io, 1, MPI_INT, mpirank + 1, group_id, mpicomm);
    YMIR_CHECK_MPI (mpiret);
  }

  /* stop performance counters */
  ymir_perf_counter_stop_add (
      &slabs_io_counter[SLABS_IO_COUNTER_GATHER_WRITE_VECTOR]);

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}
#endif

#else /* if defined SLABS_IO_USE_MPI_IO */

/**
 *
 */
static MPI_File *
slabs_io_mpi_read_vector_begin (char *binaryfile, double *vec_data,
                                int vec_size, MPI_Comm mpicomm)
{
  const char         *this_fn_name = "slabs_io_mpi_read_vector_begin";
  MPI_File           *fh = YMIR_ALLOC (MPI_File, 1);
  int                 mpiret;

  YMIR_GLOBAL_INFOF ("%s: Begin reading vector from binary file \"%s\"\n",
                     this_fn_name, binaryfile);

  /* open file */
  mpiret = MPI_File_open (mpicomm, binaryfile, MPI_MODE_RDONLY, MPI_INFO_NULL,
                          fh);
  YMIR_CHECK_MPI (mpiret);

  /* begin reading file */
  mpiret = MPI_File_read_all_begin (*fh, vec_data, vec_size, MPI_DOUBLE);
  YMIR_CHECK_MPI (mpiret);

  /* return file handle */
  return fh;
}

/**
 *
 */
static int
slabs_io_mpi_read_vector_end (char *binaryfile, double *vec_data,
                              int vec_size, MPI_File *fh)
{
  const char         *this_fn_name = "slabs_io_mpi_read_vector_end";
  MPI_Status          status;
  int                 mpiret;
  int                 count;

  /* end reading file */
  mpiret = MPI_File_read_all_end (*fh, vec_data, &status);
  YMIR_CHECK_MPI (mpiret);

  /* check number of read entries */
  mpiret = MPI_Get_count (&status, MPI_DOUBLE, &count);
  YMIR_CHECK_MPI (mpiret);

  /* close file */
  mpiret = MPI_File_close (fh); YMIR_CHECK_MPI (mpiret);
  YMIR_FREE (fh);

  YMIR_GLOBAL_INFOF ("%s: End reading vector from binary file \"%s\"\n",
                     this_fn_name, binaryfile);

  /* return number of read entries */
  return count;
}

/**
 *
 */
static void
slabs_io_mpi_read_vector_bcast (char *binaryfile, double *vec_data,
                                int vec_size, MPI_Comm mpicomm)
{
  const char         *this_fn_name = "slabs_io_mpi_read_vector_bcast";
  int                 mpirank;
  int                 mpiret;

  YMIR_GLOBAL_INFOF ("%s: Begin reading vector from binary file \"%s\"\n",
                     this_fn_name, binaryfile);

  /* get parallel environment */
  mpiret = MPI_Comm_rank (mpicomm, &mpirank); YMIR_CHECK_MPI (mpiret);

  /* read file with one processor */
  if (mpirank == 0) {
    MPI_File            fh_bin;
    MPI_Status          status;
    int                 count;

    /* open binary file (file access independently from process to process) */
    mpiret = MPI_File_open (MPI_COMM_SELF, binaryfile, MPI_MODE_RDONLY,
                            MPI_INFO_NULL, &fh_bin);
    YMIR_CHECK_MPI (mpiret);

    /* write vector to file */
    mpiret = MPI_File_read (fh_bin, vec_data, vec_size, MPI_DOUBLE, &status);
    YMIR_CHECK_MPI (mpiret);

    /* check number of read entries */
    mpiret = MPI_Get_count (&status, MPI_DOUBLE, &count);
    YMIR_CHECK_MPI (mpiret);
    if (count != vec_size) {
      YMIR_LERRORF ("%s: Wrong count %i for file read. Should be %i\n",
                    this_fn_name, count, vec_size);
    }

    /* close file */
    mpiret = MPI_File_close (&fh_bin);
    YMIR_CHECK_MPI (mpiret);
  }

  /* broadcast data */
  mpiret = MPI_Bcast (vec_data, vec_size, MPI_DOUBLE, 0, mpicomm);
  YMIR_CHECK_MPI (mpiret);

  YMIR_GLOBAL_INFOF ("%s: End reading vector from binary file \"%s\"\n",
                     this_fn_name, binaryfile);
}

/**
 * Reads double vector from binary file in parallel.
 */
static int
slabs_io_mpi_read_vector_segment (char *binaryfile, const int offset,
                                  double *vec_data, int vec_size,
                                  MPI_Comm mpicomm)
{
  const char         *this_fn_name = "slabs_io_mpi_read_vector_segment";
  MPI_File            fh;
  MPI_Status          status;
  int                 mpiret;
  int                 count;

  YMIR_GLOBAL_INFOF ("%s: Begin reading vector from binary file \"%s\"\n",
                     this_fn_name, binaryfile);

  /* open file */
  mpiret = MPI_File_open (mpicomm, binaryfile, MPI_MODE_RDONLY, MPI_INFO_NULL,
                          &fh);
  YMIR_CHECK_MPI (mpiret);

  /* read file */
  mpiret = MPI_File_read_at_all (fh, offset * sizeof (double),
                                 vec_data, vec_size,
                                 MPI_DOUBLE, &status);
  YMIR_CHECK_MPI (mpiret);

  /* check number of read entries */
  mpiret = MPI_Get_count (&status, MPI_DOUBLE, &count);
  YMIR_CHECK_MPI (mpiret);

  //YMIR_TRACEF ("%s: Vector with %i of %i elements (offset %i) read from "
  //             "binary file %s\n",
  //             this_fn_name, count, vec_size, offset, binaryfile);

  /* close file */
  mpiret = MPI_File_close (&fh); YMIR_CHECK_MPI (mpiret);

  YMIR_GLOBAL_INFOF ("%s: End reading vector from binary file \"%s\"\n",
                     this_fn_name, binaryfile);

  /* return number of read entries */
  return count;
}

#endif /* if defined SLABS_IO_USE_MPI_IO */

#if !defined(SLABS_IO_USE_MPI_IO)

/**
 * Copies temperature vector from text file to binary file.
 */
static long long int
slabs_io_copy_temperature_from_txt_to_bin (const char *textfile,
                                           const char *binaryfile)
{
  const char         *this_fn_name =
                        "slabs_io_copy_temperature_from_txt_to_bin";
  FILE               *fp_txt, *fp_bin;
  int                 status;
  long long int       node_count;
  double              temp;

  YMIR_INFOF ("Into %s: Copy temperature vector from text file \"%s\" "
              "to binary file \"%s\"\n", this_fn_name, textfile, binaryfile);

  /* open files */
  fp_txt = slabs_io_file_open (textfile, "r", this_fn_name);
  fp_bin = slabs_io_file_open (binaryfile, "wb", this_fn_name);

  /* copy temperature from text file to binary file */
  while (!feof (fp_txt)) { /* while not at the end of file */
    status = fscanf (fp_txt, "%lli %lf", &node_count, &temp);
    if (status == 2) {
      fwrite (&temp, sizeof (double), 1, fp_bin);
    }
  }

  /* close files */
  slabs_io_file_close (fp_txt, textfile, this_fn_name);
  slabs_io_file_close (fp_bin, binaryfile, this_fn_name);

  YMIR_INFOF ("Done %s\n", this_fn_name);

  /* return the number of read nodes */
  return node_count;
}

#else /* if SLABS_IO_USE_MPI_IO */

/**
 * Copies temperature vector from text file to binary file using MPI I/O.
 */
static void
slabs_io_mpi_copy_temperature_from_txt_to_bin (const char *textfile,
                                               const long long int n_values,
                                               char *binaryfile)
{
  const char         *this_fn_name =
                        "slabs_io_mpi_copy_temperature_from_txt_to_bin";
  FILE               *fp_txt;
  MPI_File            fh_bin;
  MPI_Status          status;
  int                 mpiret;
  int                 count;
  long long int       node, i;
  double             *temp;
  double             *temp_pos;

  YMIR_INFOF ("Into %s: Copy temperature vector from text file \"%s\" "
              "to binary file \"%s\"\n", this_fn_name, textfile, binaryfile);

  /* create temperature vector */
  temp = YMIR_ALLOC (double, n_values);

  /* open text file */
  fp_txt = slabs_io_file_open (textfile, "r", this_fn_name);

  /* copy temperature from text file into vector */
  temp_pos = temp;
  for (i = 0; i < n_values; i++) {
    YMIR_ASSERT (!feof (fp_txt));

    fscanf (fp_txt, "%lli %lf", &node, temp_pos);
    YMIR_ASSERT (node == i);

    temp_pos++;
  }

  /* close text file */
  slabs_io_file_close (fp_txt, textfile, this_fn_name);

  /* open binary file (file access independently from process to process) */
  mpiret = MPI_File_open (MPI_COMM_SELF, binaryfile,
                          MPI_MODE_WRONLY | MPI_MODE_CREATE,
                          MPI_INFO_NULL, &fh_bin);
  YMIR_CHECK_MPI (mpiret);

  /* write vector to file */
  mpiret = MPI_File_write (fh_bin, temp, n_values, MPI_DOUBLE, &status);
  YMIR_CHECK_MPI (mpiret);

  /* check number of read entries */
  mpiret = MPI_Get_count (&status, MPI_DOUBLE, &count);
  YMIR_CHECK_MPI (mpiret);
  if (count != n_values) {
    YMIR_LERRORF ("%s: Wrong count %i on file write. Should be %lli\n",
                  this_fn_name, count, n_values);
  }

  /* close file */
  mpiret = MPI_File_close (&fh_bin);
  YMIR_CHECK_MPI (mpiret);

  YMIR_INFOF ("Done %s\n", this_fn_name);

  /* destroy */
  YMIR_FREE (temp);
}

#endif /* end SLABS_IO_USE_MPI_IO */

/**
 *
 */
#if !defined(SLABS_IO_USE_MPI_IO)
void
slabs_io_read_temperature (sc_dmatrix_t *temperature,
                           mangll_t *mangll, mangll_cnodes_t *cnodes,
                           const char *filename_txt, const char *filename_bin)
{
  MPI_Comm            mpicomm = mangll->mpicomm;
  const int           mpisize = mangll->mpisize;
  const int           mpirank = mangll->mpirank;
  int                 mpiret;

  mangll_locidx_t    *n_cnodes = cnodes->Ngo;
  const ymir_locidx_t  n_local_cnodes = cnodes->Ncn;
  const int           n_fields = temperature->n;
  double             *temp_data = temperature->e[0];

  ymir_locidx_t      *n_local_owned_cnodes;
  ymir_gloidx_t      *global_offset;
  int                 r;

  /* check input */
  YMIR_ASSERT (temperature->m == (sc_bint_t) n_local_cnodes);

  /* start performance counters */
  ymir_perf_counter_start (&slabs_io_counter[SLABS_IO_COUNTER_TEMPERATURE]);

  /* initialize counters for continuous nodes */
  n_local_owned_cnodes = YMIR_ALLOC (ymir_locidx_t, mpisize);
  global_offset = YMIR_ALLOC (ymir_gloidx_t, mpisize + 1);
  global_offset[0] = 0;

  /* count continuous nodes */
  for (r = 0; r < mpisize; r++) {
    n_local_owned_cnodes[r] = (ymir_locidx_t) n_cnodes[r];
    global_offset[r + 1] = global_offset[r] + (ymir_gloidx_t) n_cnodes[r];
  }

  /* write temperature data from text file to binary file */
  if (filename_txt != NULL) {
    if (mpirank == 0) {
#ifdef YMIR_DEBUG
      long long int       n_nodes_read =
#endif
      slabs_io_copy_temperature_from_txt_to_bin (filename_txt, filename_bin);
      YMIR_ASSERT (n_nodes_read == global_offset[mpisize]);
    }
    mpiret = MPI_Barrier (mpicomm); YMIR_CHECK_MPI (mpiret);
  }

  /* read temperature data from binary file */
  slabs_io_read_scatter_vector (filename_bin, temp_data, global_offset,
                                n_local_owned_cnodes, mpicomm);

  /* communicate shared node values */
  {
    sc_array_t         *temp_array;

    temp_array = sc_array_new_data (temp_data,
                                    (size_t) n_fields * sizeof (double),
                                    (size_t) n_local_cnodes);
    mangll_cnodes_share_owned (temp_array, cnodes);
    sc_array_destroy (temp_array);
  }

  /* destroy */
  YMIR_FREE (global_offset);
  YMIR_FREE (n_local_owned_cnodes);

  /* stop performance counters */
  ymir_perf_counter_stop_add (&slabs_io_counter[SLABS_IO_COUNTER_TEMPERATURE]);
}
#else /* if SLABS_IO_USE_MPI_IO */
void
slabs_io_read_temperature (sc_dmatrix_t *temperature,
                           mangll_t *mangll, mangll_cnodes_t *cnodes,
                           char *filename_txt, char *filename_bin)
{
  MPI_Comm            mpicomm = mangll->mpicomm;
  const int           mpisize = mangll->mpisize;
  const int           mpirank = mangll->mpirank;
  int                 mpiret;

  const int           offset = (int) cnodes->offset;
  const int           n_local_owned_cnodes = (int) cnodes->owncount;
  mangll_locidx_t    *n_cnodes = cnodes->Ngo;
  const int           n_local_cnodes = (int) cnodes->Ncn;
  long long int       n_global_cnodes;
  const int           n_fields = temperature->n;
  double             *temp_data = temperature->e[0];
  sc_array_t         *temp_array;
  int                 p;

  /* check input */
  YMIR_ASSERT (temperature->m == n_local_cnodes);

  /* start performance counters */
  ymir_perf_counter_start (&slabs_io_counter[SLABS_IO_COUNTER_TEMPERATURE]);

  /* initialize counters for continuous nodes */
  n_global_cnodes = 0;

  /* count continuous nodes */
  for (p = 0; p < mpisize; p++) {
    /* update number of global continuous nodes */
    n_global_cnodes += (long long int) n_cnodes[p];
  }
  YMIR_ASSERT (n_global_cnodes <= INT_MAX);

  /* write temperature data from text file to binary file */
  if (filename_txt != NULL) {
    if (mpirank == 0) {
      slabs_io_mpi_copy_temperature_from_txt_to_bin (
          filename_txt, n_global_cnodes, filename_bin);
    }
    mpiret = MPI_Barrier (mpicomm); YMIR_CHECK_MPI (mpiret);
  }

  /* read temperature data from binary file */
  slabs_io_mpi_read_vector_segment (
      filename_bin, offset, temp_data, n_local_owned_cnodes, mpicomm);

  /* create array view of temperature field */
  temp_array = sc_array_new_data (temp_data,
                                  (size_t) n_fields * sizeof (double),
                                  (size_t) n_local_cnodes);

  /* communicate shared node values */
  mangll_cnodes_share_owned (temp_array, cnodes);

  /* destroy */
  sc_array_destroy (temp_array);

  /* stop performance counters */
  ymir_perf_counter_stop_add (&slabs_io_counter[SLABS_IO_COUNTER_TEMPERATURE]);
}
#endif /* end SLABS_IO_USE_MPI_IO */

/**
 * Copies weak zone point cloud from text file to binary file.
 */
static int
slabs_io_copy_weakzone_from_txt_to_bin (const char *textfile,
                                        const char *binaryfile)
{
  const char         *this_fn_name = "slabs_io_copy_weakzone_from_txt_to_bin";
  FILE               *fp_txt, *fp_bin;
  int                 status;
  int                 point_count = 0;
  double              coord[3];

  YMIR_INFOF ("Into %s: Copy weak zone point cloud "
              "from text file \"%s\" to binary file \"%s\"\n",
              this_fn_name, textfile, binaryfile);

  /* open files */
  fp_txt = slabs_io_file_open (textfile, "r", this_fn_name);
  fp_bin = slabs_io_file_open (binaryfile, "wb", this_fn_name);

  /* copy weak zone coordinates from text file to binary file */
  while (!feof (fp_txt)) { /* while not at the end of file */
    status = fscanf (fp_txt, "%lf %lf %lf", &coord[0], &coord[1], &coord[2]);
    if (status == 3) {
      fwrite (&coord, sizeof (double), 3, fp_bin);
      point_count++;
    }
  }

  YMIR_INFOF ("%s: %i points copied\n", this_fn_name, point_count);

  /* close files */
  slabs_io_file_close (fp_txt, textfile, this_fn_name);
  slabs_io_file_close (fp_bin, binaryfile, this_fn_name);

  YMIR_INFOF ("Done %s\n", this_fn_name);

  /* return the number of points in the cloud */
  return point_count;
}

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
                        int *n_points_txt, const int n_points_bin)
{
  int                 mpirank, mpiret;
  int                 n_points;
  double             *point;

  /* check input */
  YMIR_ASSERT (filename_bin != NULL);

  /* start performance counters */
  ymir_perf_counter_start (&slabs_io_counter[SLABS_IO_COUNTER_WEAKZONE]);

  /* get parallel environment */
  mpiret = MPI_Comm_rank (mpicomm, &mpirank); YMIR_CHECK_MPI (mpiret);

  /* write weak zone point coordinates from text file to binary file */
  if (filename_txt != NULL) {
    if (mpirank == 0) {
      n_points = slabs_io_copy_weakzone_from_txt_to_bin (filename_txt,
                                                         filename_bin);
      YMIR_ASSERT (n_points_bin <= 0 || n_points_bin == n_points);
    }
    mpiret = MPI_Barrier (mpicomm); YMIR_CHECK_MPI (mpiret);

    /* communicate number of points in txt file */
    mpiret = MPI_Bcast (&n_points, 1, MPI_INT, 0, mpicomm);
    YMIR_CHECK_MPI (mpiret);

    /* set output: #points in txt file */
    if (n_points_txt != NULL) {
      *n_points_txt = n_points;
    }
  }
  else {
    n_points = n_points_bin;

    /* set output: #points in txt file */
    if (n_points_txt != NULL) {
      *n_points_txt = -1;
    }
  }

  /* initialize storage for point coordinates */
  YMIR_ASSERT (0 < n_points);
  point = YMIR_ALLOC (double, 3 * n_points);

  /* read weak zone point coordinates from binary file on all processors */
#ifndef SLABS_IO_USE_MPI_IO
  slabs_io_read_bcast_vector_short (filename_bin, point, 3 * n_points, mpicomm);
#else
  slabs_io_mpi_read_vector_segment (filename_bin, 0, point, 3 * n_points,
                                    mpicomm);
#endif

  /* stop performance counters */
  ymir_perf_counter_stop_add (&slabs_io_counter[SLABS_IO_COUNTER_WEAKZONE]);

  /* return point coordinates */
  return point;
}

void
slabs_io_perf_counter_init (const int active)
{
  ymir_perf_counter_init_all (slabs_io_counter, slabs_io_counter_name,
                              SLABS_IO_COUNTER_N, active);
}

void
slabs_io_perf_counter_gather (MPI_Comm mpicomm, const int print_wtime,
                              const int print_n_calls, const int print_flops)
{
  slabs_io_n_stats = ymir_perf_counter_gather_stats (
        slabs_io_counter, SLABS_IO_COUNTER_N,
        slabs_io_stats, slabs_io_stats_name, mpicomm,
        print_wtime, print_n_calls, print_flops);
}

void
slabs_io_perf_counter_print ()
{
  ymir_perf_counter_print_stats (slabs_io_stats, slabs_io_n_stats, "Slabs IO");
}

/**
 *TODO delete
 */
#if 0
MPI_File           *slabs_io_read_weakzone_fh = NULL;
double             *slabs_io_read_weakzone_point = NULL;

/**
 *
 */
int
slabs_io_read_weakzone_begin (MPI_Comm mpicomm,
                              slabs_physics_options_t *physics_options)
{
  int                 mpirank;
  int                 mpiret;
  MPI_File           *fh;

  const slabs_weakzone_t  weakzone_type = physics_options->weakzone_type;
  char               *filename_txt =
                        physics_options->weakzone_import_filename_txt;
  char               *filename_bin =
                        physics_options->weakzone_import_filename_bin;
  int                 n_points =
                        physics_options->weakzone_import_n_points;
  double             *point;

  /* start performance counters */
  if (slabs_io_counter[SLABS_IO_COUNTER_WEAKZONE].active) {
    ymir_perf_counter_start (&slabs_io_counter[SLABS_IO_COUNTER_WEAKZONE]);
  }

  /* get parallel environment */
  mpiret = MPI_Comm_rank (mpicomm, &mpirank); YMIR_CHECK_MPI (mpiret);

  /* write weak zone point coordinates from text file to binary file */
  if (   mpirank == 0
      && weakzone_type == SL_WEAKZONE_IMPORT_FILE && filename_txt != NULL ) {
    n_points = slabs_io_copy_weakzone_from_txt_to_bin (filename_txt,
                                                       filename_bin);
    physics_options->weakzone_import_filename_txt = NULL;
    YMIR_ASSERT (physics_options->weakzone_import_n_points <= 0 ||
                 physics_options->weakzone_import_n_points == n_points);
  }
  mpiret = MPI_Barrier (mpicomm); YMIR_CHECK_MPI (mpiret);

  /* initialize storage for point coordinates */
  YMIR_ASSERT (0 < n_points);
  point = YMIR_ALLOC (double, 3 * n_points);

  /* begin reading weak zone point coordinates from binary file on all procs */
  fh = slabs_io_read_vector_begin (filename_bin, point, 3 * n_points, mpicomm);

  /* set global variables */
  slabs_io_read_weakzone_fh = fh;
  slabs_io_read_weakzone_point = point;

  /* stop performance counters */
  if (slabs_io_counter[SLABS_IO_COUNTER_WEAKZONE].active) {
    ymir_perf_counter_stop_add (&slabs_io_counter[SLABS_IO_COUNTER_WEAKZONE]);
  }

  /* return number of point coordinates */
  return n_points;
}

/**
 *
 */
double *
slabs_io_read_weakzone_end (slabs_physics_options_t *physics_options)
{
  char               *filename_bin =
                        physics_options->weakzone_import_filename_bin;
  int                 n_points =
                        physics_options->weakzone_import_n_points;
  MPI_File           *fh = slabs_io_read_weakzone_fh;
  double             *point = slabs_io_read_weakzone_point;

  /* start performance counters */
  if (slabs_io_counter[SLABS_IO_COUNTER_WEAKZONE].active) {
    ymir_perf_counter_start (&slabs_io_counter[SLABS_IO_COUNTER_WEAKZONE]);
  }

  /* end reading weak zone point coordinates from binary file on all procs */
  slabs_io_read_vector_end (filename_bin, point, 3 * n_points, fh);

  /* stop performance counters */
  if (slabs_io_counter[SLABS_IO_COUNTER_WEAKZONE].active) {
    ymir_perf_counter_stop_add (&slabs_io_counter[SLABS_IO_COUNTER_WEAKZONE]);
  }

  /* return point coordinates */
  return point;
}
#endif
