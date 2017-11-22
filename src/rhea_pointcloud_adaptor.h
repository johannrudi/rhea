/*
*/

#ifndef RHEA_POINTCLOUD_ADAPTOR_H
#define RHEA_POINTCLOUD_ADAPTOR_H

#include <sc.h>

/* tell the C++ compiler to use the C style name mangling */
SC_EXTERN_C_BEGIN;

/******************************************************************************
 * Point Cloud
 *****************************************************************************/

typedef struct rhea_pointcloud_Cloud rhea_pointcloud_Cloud;

rhea_pointcloud_Cloud *rhea_pointcloud_Cloud_new ();

void                rhea_pointcloud_Cloud_destroy (
                                          rhea_pointcloud_Cloud *cloud);

void                rhea_pointcloud_Cloud_set_points (
                                          rhea_pointcloud_Cloud *cloud,
                                          const double *point,
                                          const size_t n_points);

/******************************************************************************
 * KD-Tree
 *****************************************************************************/

typedef struct rhea_pointcloud_KDTree rhea_pointcloud_KDTree;

rhea_pointcloud_KDTree *rhea_pointcloud_KDTree_new (
                                          const rhea_pointcloud_Cloud *cloud);

void                rhea_pointcloud_KDTree_destroy (
                                          rhea_pointcloud_KDTree *tree);

double              rhea_pointcloud_KDTree_find_shortest_distance_single (
                                          rhea_pointcloud_KDTree *tree,
                                          const double *target_pt);

void                rhea_pointcloud_KDTree_find_shortest_distance_multi (
                                          rhea_pointcloud_KDTree *tree,
                                          double *dist,
                                          const double *target_x,
                                          const double *target_y,
                                          const double *target_z,
                                          const unsigned int n_target_points);

SC_EXTERN_C_END;

#endif /* RHEA_POINTCLOUD_ADAPTOR_H */
