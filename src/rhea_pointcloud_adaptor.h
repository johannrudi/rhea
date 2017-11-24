/*
*/

#ifndef RHEA_POINTCLOUD_ADAPTOR_H
#define RHEA_POINTCLOUD_ADAPTOR_H

#include <sc.h>

/* tell the C++ compiler to use the C style name mangling */
SC_EXTERN_C_BEGIN;

/******************************************************************************
 * Weak Zone Point Cloud
 *****************************************************************************/

/* weak zone point cloud (opaque) */
typedef struct rhea_pointcloud_weakzone rhea_pointcloud_weakzone_t;

/**
 * Creates a new weak zone point cloud object.
 */
rhea_pointcloud_weakzone_t *rhea_pointcloud_weakzone_new (
                                        const double *point_coordinates,
                                        const size_t n_points);

/**
 * Destroys a weak zone point cloud object.
 */
void                rhea_pointcloud_weakzone_destroy (
                                        rhea_pointcloud_weakzone_t *ptcl_weak);

/**
 * Sets coordinates of the points and adjusts the size of the cloud.
 */
void                rhea_pointcloud_weakzone_set_coordinates (
                                        rhea_pointcloud_weakzone_t *ptcl_weak,
                                        const double *point_coordinates,
                                        const size_t n_points);

/**
 * Sets factors corresponding to the points.
 * Assumes that the number of factors equals the point count.
 */
void                rhea_pointcloud_weakzone_set_factors (
                                        rhea_pointcloud_weakzone_t *ptcl_weak,
                                        const double *factors);

/**
 * Sets labels corresponding to the points.
 * Assumes that the number of factors equals the point count.
 */
void                rhea_pointcloud_weakzone_set_labels (
                                        rhea_pointcloud_weakzone_t *ptcl_weak,
                                        const int *labels);

/**
 * Finds point of the cloud that is nearest to a target point. Returns distance
 * to nearest point.
 */
double              rhea_pointcloud_weakzone_find_nearest (
                                        double *nearest_coordinates,
                                        double *nearest_factor,
                                        int *nearest_label,
                                        rhea_pointcloud_weakzone_t *ptcl_weak,
                                        const double *target_coordinates);

/**
 * Finds multiple points of the cloud that are nearest to a target point.
 */
int                 rhea_pointcloud_weakzone_find_n_nearest (
                                        double *nearest_dist,
                                        double *nearest_coordinates,
                                        double *nearest_factor,
                                        int *nearest_label,
                                        const int n_nearest,
                                        rhea_pointcloud_weakzone_t *ptcl_weak,
                                        const double *target_coordinates);

/******************************************************************************
 * Topography Point Cloud
 *****************************************************************************/

/* topography point cloud (opaque) */
typedef struct rhea_pointcloud_topography rhea_pointcloud_topography_t;

//TODO rhea_pointcloud_topography_new
//TODO rhea_pointcloud_topography_destroy
//TODO rhea_pointcloud_topography_set_coordinates (2D?)
//TODO rhea_pointcloud_topography_set_displacements
//TODO rhea_pointcloud_topography_set_labels




//TODO deprecated below:

/******************************************************************************
 * Point Cloud
 *****************************************************************************/

typedef struct rhea_pointcloud_Cloud rhea_pointcloud_Cloud_t;

rhea_pointcloud_Cloud_t *rhea_pointcloud_Cloud_new ();

void                rhea_pointcloud_Cloud_destroy (
                                          rhea_pointcloud_Cloud_t *cloud);

void                rhea_pointcloud_Cloud_set_points (
                                          rhea_pointcloud_Cloud_t *cloud,
                                          const double *point,
                                          const size_t n_points);

/******************************************************************************
 * KD-Tree
 *****************************************************************************/

typedef struct rhea_pointcloud_KDTree rhea_pointcloud_KDTree_t;

rhea_pointcloud_KDTree_t *rhea_pointcloud_KDTree_new (
                                          const rhea_pointcloud_Cloud_t *cloud);

void                rhea_pointcloud_KDTree_destroy (
                                          rhea_pointcloud_KDTree_t *tree);

double              rhea_pointcloud_KDTree_find_shortest_distance_single (
                                          rhea_pointcloud_KDTree_t *tree,
                                          const double *target_pt);

void                rhea_pointcloud_KDTree_find_shortest_distance_multi (
                                          rhea_pointcloud_KDTree_t *tree,
                                          double *dist,
                                          const double *target_x,
                                          const double *target_y,
                                          const double *target_z,
                                          const unsigned int n_target_points);

SC_EXTERN_C_END;

#endif /* RHEA_POINTCLOUD_ADAPTOR_H */
