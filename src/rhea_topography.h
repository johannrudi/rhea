/*
 */

#ifndef RHEA_TOPOGRAPHY_H
#define RHEA_TOPOGRAPHY_H

#include <rhea_domain.h>
#include <rhea_pointcloud_adaptor.h>
#include <ymir_vec_ops.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* options for topography */
typedef struct rhea_topography_options
{
  /* binary/text files with coordinates, labels, displ. of topography points */
  char               *points_file_path_bin;
  char               *points_file_path_txt;
  char               *labels_file_path_bin;
  char               *labels_file_path_txt;
  char               *displacements_file_path_bin;
  char               *displacements_file_path_txt;

  /* number of points */
  int                 n_points;

  /* output */
  char               *write_points_file_path_bin;
  char               *write_labels_file_path_bin;
  char               *write_displacements_file_path_bin;

  /* data */
  rhea_pointcloud_topography_t *pointcloud;

  /* options & properties of the computational domain */
  rhea_domain_options_t *domain_options;
}
rhea_topography_options_t;

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_topography_add_options (ymir_options_t * opt_sup);

/**
 * Processes options and stores them.
 */
void                rhea_topography_process_options (
                                        rhea_topography_options_t *opt,
                                        rhea_domain_options_t *domain_options);

/******************************************************************************
 * Topography Vector
 *****************************************************************************/

/**
 * Creates a new topography vector.
 */
ymir_vec_t         *rhea_topography_new (ymir_mesh_t *ymir_mesh);

/**
 * Destroys a topography vector.
 */
void                rhea_topography_destroy (ymir_vec_t *topography);

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_topography_check_vec_type (ymir_vec_t *vec);

/**
 * Checks entries of a vector.
 */
int                 rhea_topography_is_valid (ymir_vec_t *vec);

/******************************************************************************
 * Data
 *****************************************************************************/

/**
 * Allocates and performs the setup of topography data.
 */
void                rhea_topography_data_create (rhea_topography_options_t *opt,
                                                 sc_MPI_Comm mpicomm);

/**
 * Clears storage of topography data.
 */
void                rhea_topography_data_clear (rhea_topography_options_t *opt);

/******************************************************************************
 * Topography Computation
 *****************************************************************************/

//TODO

#endif /* RHEA_TOPOGRAPHY_H */
