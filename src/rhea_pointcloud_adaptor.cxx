/*
*/

#include <rhea_pointcloud_adaptor.h>
#include <rhea_pointcloud.hpp>

SC_EXTERN_C_BEGIN;

/******************************************************************************
 * Point Cloud
 *****************************************************************************/

rhea_pointcloud_Cloud *
rhea_pointcloud_Cloud_new ()
{
  return new rhea_pointcloud_Cloud ();
}

void
rhea_pointcloud_Cloud_destroy (rhea_pointcloud_Cloud *cloud)
{
  delete cloud;
}

void
rhea_pointcloud_Cloud_set_points (rhea_pointcloud_Cloud *cloud,
                                  const double *point, const size_t n_points)
{
  cloud->set_points (point, n_points);
}

/******************************************************************************
 * KD-Tree
 *****************************************************************************/

rhea_pointcloud_KDTree *
rhea_pointcloud_KDTree_new (const rhea_pointcloud_Cloud *cloud)
{
  return new rhea_pointcloud_KDTree (*cloud);
}

void
rhea_pointcloud_KDTree_destroy (rhea_pointcloud_KDTree *tree)
{
  delete tree;
}

double
rhea_pointcloud_KDTree_find_shortest_distance_single (
                                          rhea_pointcloud_KDTree *tree,
                                          const double *target_pt)
{
  return tree->find_shortest_distance_single (target_pt);
}

void
rhea_pointcloud_KDTree_find_shortest_distance_multi (
                                          rhea_pointcloud_KDTree *tree,
                                          double *dist,
                                          const double *target_x,
                                          const double *target_y,
                                          const double *target_z,
                                          const unsigned int n_target_points)
{
  tree->find_shortest_distance_multi (dist, target_x, target_y, target_z,
                                      n_target_points);
}

SC_EXTERN_C_END;
