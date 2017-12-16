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

/* types of topography */
typedef enum
{
  RHEA_TOPOGRAPHY_NONE = -1,                /* topography does not exist */
  RHEA_TOPOGRAPHY_DATA_POINTS_DISPLS = 0,   /* points (coord, displ) */
  RHEA_TOPOGRAPHY_DATA_POINTS_DISPLS_LABELS /* points (coord, displ, labels) */
}
rhea_topography_t;

/* options for topography */
typedef struct rhea_topography_options
{
  /* type of topography */
  rhea_topography_t   type;

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

/**
 * Checks whether a topography exists for the given set of options.
 */
int                 rhea_topography_exists (rhea_topography_options_t *opt);

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

/**
 * Calculates topography displacement for a single point with coordinates
 * (x,y,z).
 */
double              rhea_topography_displacement_node (
                                              int *nearest_label,
                                              double x, double y, double z,
                                              rhea_topography_options_t *opt);

/**
 * Fills values of a vector with displacement calculations.
 */
void                rhea_topography_displacement_vec (
                                              ymir_vec_t *displacement,
                                              rhea_topography_options_t *opt);

#endif /* RHEA_TOPOGRAPHY_H */
