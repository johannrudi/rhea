/*
*/

#include <rhea_pointcloud_adaptor.h>
#include <rhea_pointcloud.hpp>
#include <rhea_base.h>

SC_EXTERN_C_BEGIN;

/******************************************************************************
 * Weak Zone Point Cloud
 *****************************************************************************/

/* weak zone point cloud */
struct rhea_pointcloud_weakzone
{
  rhea_pointcloud_Cloud  *cloud;
  rhea_pointcloud_KDTree *tree;
};

rhea_pointcloud_weakzone_t *
rhea_pointcloud_weakzone_new (const double *point_coordinates,
                              const int *point_labels,
                              const size_t n_points)
{
  rhea_pointcloud_weakzone_t *ptcl_weak;

  /* create weak zone object */
  ptcl_weak = RHEA_ALLOC (rhea_pointcloud_weakzone_t, 1);

  /* create point cloud and fill with points */
  ptcl_weak->cloud = new rhea_pointcloud_Cloud ();
  ptcl_weak->cloud->set_points (point_coordinates, n_points);

  /* create kd-tree corresponding to point cloud */
  ptcl_weak->tree = new rhea_pointcloud_KDTree (*(ptcl_weak->cloud));

  return ptcl_weak;
}

void
rhea_pointcloud_weakzone_destroy (rhea_pointcloud_weakzone_t *ptcl_weak)
{
  delete ptcl_weak->tree;
  delete ptcl_weak->cloud;
  RHEA_FREE (ptcl_weak);
}




//TODO deprecated below:

/******************************************************************************
 * Point Cloud
 *****************************************************************************/

rhea_pointcloud_Cloud_t *
rhea_pointcloud_Cloud_new ()
{
  return new rhea_pointcloud_Cloud ();
}

void
rhea_pointcloud_Cloud_destroy (rhea_pointcloud_Cloud_t *cloud)
{
  delete cloud;
}

void
rhea_pointcloud_Cloud_set_points (rhea_pointcloud_Cloud_t *cloud,
                                  const double *point, const size_t n_points)
{
  cloud->set_points (point, n_points);
}

/******************************************************************************
 * KD-Tree
 *****************************************************************************/

rhea_pointcloud_KDTree_t *
rhea_pointcloud_KDTree_new (const rhea_pointcloud_Cloud_t *cloud)
{
  return new rhea_pointcloud_KDTree (*cloud);
}

void
rhea_pointcloud_KDTree_destroy (rhea_pointcloud_KDTree_t *tree)
{
  delete tree;
}

double
rhea_pointcloud_KDTree_find_shortest_distance_single (
                                          rhea_pointcloud_KDTree_t *tree,
                                          const double *target_pt)
{
  return tree->find_shortest_distance (target_pt);
}

void
rhea_pointcloud_KDTree_find_shortest_distance_multi (
                                          rhea_pointcloud_KDTree_t *tree,
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
