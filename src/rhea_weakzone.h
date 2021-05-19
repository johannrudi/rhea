/** RHEA_WEAKZONE
 *
 * Weak zones to decouple plate boundaries.
 */

#ifndef RHEA_WEAKZONE_H
#define RHEA_WEAKZONE_H

#include <rhea_domain.h>
#include <rhea_pointcloud_adaptor.h>
#include <rhea_weakzone_label.h>
#include <ymir_vec_ops.h>

/* constant: neutral/default value for weak zone (i.e., no weakening) */
#define RHEA_WEAKZONE_NEUTRAL_VALUE (1.0)

/******************************************************************************
 * Options
 *****************************************************************************/

/* types of weak zones */
typedef enum
{
  RHEA_WEAKZONE_NONE = -1,                  /* weak zone does not exist */
  RHEA_WEAKZONE_DEPTH = 0,                  /* horizontal slab along surface */
  RHEA_WEAKZONE_DATA_POINTS,                /* points (coordinates) */
  RHEA_WEAKZONE_DATA_POINTS_LABELS,         /* pt's (coord, labels) */
  RHEA_WEAKZONE_DATA_POINTS_LABELS_FACTORS, /* pt's (coord, labels, factors) */
}
rhea_weakzone_t;

/* options for weak zones */
typedef struct rhea_weakzone_options
{
  /* type of weak zone */
  rhea_weakzone_t     type;

  /* weak zone geometry */
  double              thickness;
  double              thickness_class_slab;
  double              thickness_class_ridge;
  double              thickness_class_fracture;
  double              thickness_const;
  double              thickness_const_class_slab;
  double              thickness_const_class_ridge;
  double              thickness_const_class_fracture;

  /* max weakening in the interior of weak zones;
   * factors correspond to labels in `rhea_weakzone_label.h` */
  double              weak_factor_interior;
  double              weak_factor_interior_class_slab;
  double              weak_factor_interior_class_ridge;
  double              weak_factor_interior_class_fracture;
  double             *weak_factor_interior_label;

  /* binary/text files with coordinates, labels, factors of weak zone points */
  char               *points_file_path_bin;
  char               *points_file_path_txt;
  char               *labels_file_path_bin;
  char               *labels_file_path_txt;
  char               *factors_file_path_bin;
  char               *factors_file_path_txt;

  /* number of points and labels per class */
  int                 n_points;
  int                 n_labels[RHEA_WEAKZONE_LABEL_CLASS_N];

  /* output */
  char               *write_points_file_path_bin;
  char               *write_points_file_path_txt;
  char               *write_labels_file_path_bin;
  char               *write_factors_file_path_bin;

  /* data */
  rhea_pointcloud_weakzone_t *pointcloud;

  /* statistics */
  double              stats_radius_min;
  double              stats_thickness_max;
  double              stats_factor_interior_min;

  /* options & properties of the computational domain */
  rhea_domain_options_t *domain_options;
}
rhea_weakzone_options_t;

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_weakzone_add_options (ymir_options_t * opt_sup);

/**
 * Processes options and stores them.
 */
void                rhea_weakzone_process_options (
                                        rhea_weakzone_options_t *opt,
                                        rhea_domain_options_t *domain_options);

/**
 * Checks whether a weak zone exists for the given set of options.
 */
int                 rhea_weakzone_exists (rhea_weakzone_options_t *opt);

/******************************************************************************
 * Monitoring
 *****************************************************************************/

/**
 * Initializes performance counters.
 */
void                rhea_weakzone_perfmon_init (const int activate,
                                                const int skip_if_active);

/**
 * Prints statistics collected by performance monitors.
 */
void                rhea_weakzone_perfmon_print (sc_MPI_Comm mpicomm,
                                                 const int print_wtime,
                                                 const int print_n_calls,
                                                 const int print_flops);

/******************************************************************************
 * Weak Zone Vector
 *****************************************************************************/

/**
 * Creates a new weak zone vector.
 */
ymir_vec_t         *rhea_weakzone_new (ymir_mesh_t *ymir_mesh);
ymir_vec_t         *rhea_weakzone_normal_new (ymir_mesh_t *ymir_mesh);

/**
 * Destroys a weak zone vector.
 */
void                rhea_weakzone_destroy (ymir_vec_t *weakzone);
void                rhea_weakzone_normal_destroy (ymir_vec_t *weakzone_normal);

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_weakzone_check_vec_type (ymir_vec_t *vec);
int                 rhea_weakzone_normal_check_vec_type (ymir_vec_t *vec);

/**
 * Checks entries of a vector.
 */
int                 rhea_weakzone_is_valid (ymir_vec_t *vec);

/******************************************************************************
 * Data
 *****************************************************************************/

/**
 * Allocates and performs the setup of data that is required to compute the
 * weak zone field.
 */
void                rhea_weakzone_data_create (rhea_weakzone_options_t *opt,
                                               sc_MPI_Comm mpicomm);

/**
 * Clears storage of weak zone data.
 */
void                rhea_weakzone_data_clear (rhea_weakzone_options_t *opt);

/******************************************************************************
 * Weak Zone Computation
 *****************************************************************************/

/**
 * Computes a custom weak zone.  Callback function for weak zone computation.
 */
typedef void      (*rhea_weakzone_compute_fn_t) (ymir_vec_t *weakzone,
                                                 void *data);

/**
 * Computes the weak zone at each node of the mesh.
 */
void                rhea_weakzone_compute (ymir_vec_t *weakzone, void *data);

/**
 * Computes the distances to the weak zone surfaces at each node of the mesh.
 */
void                rhea_weakzone_compute_distance (
                                                ymir_vec_t *distance,
                                                rhea_weakzone_options_t *opt);

/**
 * Computes the (smooth) indicator for weak zones that is 1 inside of weak
 * zones and goes to zero far away from weak zones.
 */
void                rhea_weakzone_compute_indicator (
                                                ymir_vec_t *indicator,
                                                const int label_filter,
                                                rhea_weakzone_options_t *opt);

/**
 * Computes the normal to the weak zone surfaces at each node of the mesh.
 */
void                rhea_weakzone_compute_normal (
                                                ymir_vec_t *weakzone_normal,
                                                rhea_weakzone_options_t *opt);

/**
 * Computes the weak zone distance/factor/factor derivative at one node.
 */
double              rhea_weakzone_dist_node (int *nearest_label,
                                             double *nearest_factor,
                                             const double x, const double y,
                                             const double z,
                                             rhea_weakzone_options_t *opt);

double              rhea_weakzone_factor_node (const double distance,
                                               const double thickness,
                                               const double thickness_const,
                                               const double factor_interior);

double              rhea_weakzone_factor_deriv_node (
                                                 const double distance,
                                                 const double thickness,
                                                 const double thickness_const,
                                                 const double factor_interior);

/**
 * Computes the weak zone factor of one element.
 */
void                rhea_weakzone_compute_elem (double *_sc_restrict weak_elem,
                                                const double *_sc_restrict x,
                                                const double *_sc_restrict y,
                                                const double *_sc_restrict z,
                                                const int n_nodes,
                                                rhea_weakzone_options_t *opt);

/******************************************************************************
 * Get & Set Parameters
 *****************************************************************************/

/**
 * Returns the total number of labels.
 */
int                 rhea_weakzone_get_total_n_labels (
                                                rhea_weakzone_options_t *opt);

/**
 * Looks up weak zone parameters.
 */
double              rhea_weakzone_lookup_thickness (
                                              const int label,
                                              rhea_weakzone_options_t *opt);

double              rhea_weakzone_lookup_thickness_const (
                                              const int label,
                                              rhea_weakzone_options_t *opt);

double              rhea_weakzone_lookup_factor_interior (
                                              const int label,
                                              rhea_weakzone_options_t *opt);

int                 rhea_weakzone_lookup_index_from_label (
                                              const rhea_weakzone_label_t label,
                                              rhea_weakzone_options_t *opt);

rhea_weakzone_label_t rhea_weakzone_lookup_label_from_index (
                                              const int idx,
                                              rhea_weakzone_options_t *opt);

/******************************************************************************
 * Get & Set Values
 *****************************************************************************/

/**
 * Gets the weak zone of one element at Gauss nodes.
 */
double             *rhea_weakzone_get_elem_gauss (sc_dmatrix_t *weak_el_mat,
                                                  ymir_vec_t *weak_vec,
                                                  const ymir_locidx_t elid);

/**
 * Sets the weak zone of one element at Gauss nodes.
 */
void                rhea_weakzone_set_elem_gauss (ymir_vec_t *weak_vec,
                                                  sc_dmatrix_t *weak_el_mat,
                                                  const ymir_locidx_t elid);

#endif /* RHEA_WEAKZONE_H */
