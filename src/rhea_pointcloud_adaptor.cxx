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
rhea_pointcloud_weakzone_new (const double x_min, const double x_max,
                              const double y_min, const double y_max,
                              const double z_min, const double z_max,
                              const double *point_coordinates,
                              const size_t n_points)
{
  rhea_pointcloud_weakzone_t *ptcl_weak;

  /* create weak zone object */
  ptcl_weak = RHEA_ALLOC (rhea_pointcloud_weakzone_t, 1);

  /* create point cloud */
  ptcl_weak->cloud = new rhea_pointcloud_Cloud ();

  /* set the bounding box containing all points */
  ptcl_weak->cloud->set_bounding_box (x_min, x_max, y_min, y_max, z_min, z_max,
                                      1 /* add a margin */);

  /* fill the cloud with points */
  ptcl_weak->cloud->set_point_coordinates_all (point_coordinates, n_points);

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

void
rhea_pointcloud_weakzone_set_coordinates (rhea_pointcloud_weakzone_t *ptcl_weak,
                                          const double *point_coordinates,
                                          const size_t n_points)
{
  ptcl_weak->cloud->set_point_coordinates_all (point_coordinates, n_points);
  ptcl_weak->tree->setup ();
}

void
rhea_pointcloud_weakzone_set_factors (rhea_pointcloud_weakzone_t *ptcl_weak,
                                      const double *factors)
{
  ptcl_weak->cloud->set_point_value_all (factors);
}

void
rhea_pointcloud_weakzone_set_labels (rhea_pointcloud_weakzone_t *ptcl_weak,
                                     const int *labels)
{
  ptcl_weak->cloud->set_point_label_all (labels);
}

double
rhea_pointcloud_weakzone_find_nearest (double *nearest_coordinates,
                                       double *nearest_factor,
                                       int *nearest_label,
                                       rhea_pointcloud_weakzone_t *ptcl_weak,
                                       const double *target_coordinates)
{
  return ptcl_weak->tree->find_nearest (nearest_coordinates, nearest_factor,
                                        nearest_label, target_coordinates);
}

int
rhea_pointcloud_weakzone_find_n_nearest (double *nearest_dist,
                                         double *nearest_coordinates,
                                         double *nearest_factor,
                                         int *nearest_label,
                                         const int n_nearest,
                                         rhea_pointcloud_weakzone_t *ptcl_weak,
                                         const double *target_coordinates)
{
  return (int) ptcl_weak->tree->find_n_nearest (
      nearest_dist, nearest_coordinates, nearest_factor, nearest_label,
      (size_t) n_nearest, target_coordinates);
}

/******************************************************************************
 * Topography Point Cloud
 *****************************************************************************/

/* weak zone point cloud */
struct rhea_pointcloud_topography
{
  rhea_pointcloud_Cloud  *cloud;
  rhea_pointcloud_KDTree *tree;
};

rhea_pointcloud_topography_t *
rhea_pointcloud_topography_new (const double x_min, const double x_max,
                                const double y_min, const double y_max,
                                const double z_min, const double z_max,
                                const double *point_coordinates,
                                const size_t n_points)
{
  rhea_pointcloud_topography_t *ptcl_topo;

  /* create weak zone object */
  ptcl_topo = RHEA_ALLOC (rhea_pointcloud_topography_t, 1);

  /* create point cloud */
  ptcl_topo->cloud = new rhea_pointcloud_Cloud ();

  /* set the bounding box containing all points */
  ptcl_topo->cloud->set_bounding_box (x_min, x_max, y_min, y_max, z_min, z_max,
                                      1 /* add a margin */);

  /* fill the cloud with points */
  ptcl_topo->cloud->set_point_coordinates_all (point_coordinates, n_points);

  /* create kd-tree corresponding to point cloud */
  ptcl_topo->tree = new rhea_pointcloud_KDTree (*(ptcl_topo->cloud));

  return ptcl_topo;
}

void
rhea_pointcloud_topography_destroy (rhea_pointcloud_topography_t *ptcl_topo)
{
  delete ptcl_topo->tree;
  delete ptcl_topo->cloud;
  RHEA_FREE (ptcl_topo);
}

void
rhea_pointcloud_topography_set_displacements (
                                      rhea_pointcloud_topography_t *ptcl_topo,
                                      const double *factors)
{
  ptcl_topo->cloud->set_point_value_all (factors);
}

void
rhea_pointcloud_topography_set_labels (rhea_pointcloud_topography_t *ptcl_topo,
                                       const int *labels)
{
  ptcl_topo->cloud->set_point_label_all (labels);
}

int
rhea_pointcloud_topography_find_n_nearest (
                                      double *nearest_dist,
                                      double *nearest_coordinates,
                                      double *nearest_factor,
                                      int *nearest_label,
                                      const int n_nearest,
                                      rhea_pointcloud_topography_t *ptcl_topo,
                                      const double *target_coordinates)
{
  return (int) ptcl_topo->tree->find_n_nearest (
      nearest_dist, nearest_coordinates, nearest_factor, nearest_label,
      (size_t) n_nearest, target_coordinates);
}

SC_EXTERN_C_END;
