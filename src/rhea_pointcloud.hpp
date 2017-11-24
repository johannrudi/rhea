/*
*/

#ifndef RHEA_POINTCLOUD_HPP
#define RHEA_POINTCLOUD_HPP

#include "nanoflann/include/nanoflann.hpp"
#include <cstdlib>

/******************************************************************************
 * Point Cloud
 *****************************************************************************/

/* spatial dimension */
#define RHEA_POINTCLOUD_DIM 3

/* margin for the limits of the bounding box containing the point cloud */
#define RHEA_POINTCLOUD_BOX_MARGIN 0.01

/* undefined label of a point */
#define RHEA_POINTCLOUD_LABEL_UNDEF (-404)

/**
 * rhea_pointcloud_point: One single point
 *
 *  x, y, z --- coordinates
 *  value   --- can be the weak zone factor or topography displacement
 *  label   --- a label to distinguish geological entities, e.g., plates
 */
struct rhea_pointcloud_point
{
  double              x, y, z;
  double              value;
  int                 label;
};

/**
 * rhea_pointcloud_Cloud: Collection of points that constitute the point cloud.
 */
class rhea_pointcloud_Cloud
{
  private:

  /* point data */
  std::vector<rhea_pointcloud_point>  point;

  /* bounding box size */
  double              bbox_x_min, bbox_x_max;
  double              bbox_y_min, bbox_y_max;
  double              bbox_z_min, bbox_z_max;

  public:

  /**
   * Constructs a new object (or an instance of the class).
   */
  rhea_pointcloud_Cloud ()
  {
    const double        xyz_limit = 1.0 + RHEA_POINTCLOUD_BOX_MARGIN;

    /* initialize limits of the bounding box */
    bbox_x_min = -xyz_limit;
    bbox_x_max = +xyz_limit;

    bbox_y_min = -xyz_limit;
    bbox_y_max = +xyz_limit;

    bbox_z_min = -xyz_limit;
    bbox_z_max = +xyz_limit;
  }

  /**
   * Destroys an object.
   */
  ~rhea_pointcloud_Cloud () {}

  /**
   * Sets limits of the bounding box.
   */
  void
  set_bounding_box (const double x_min, const double x_max,
                    const double y_min, const double y_max,
                    const double z_min, const double z_max,
                    const int add_margin)
  {
    /* set limits */
    bbox_x_min = x_min;
    bbox_x_max = x_max;

    bbox_y_min = y_min;
    bbox_y_max = y_max;

    bbox_z_min = z_min;
    bbox_z_max = z_max;

    /* add a margin to the limits */
    if (add_margin) {
      const double        margin = RHEA_POINTCLOUD_BOX_MARGIN;
      double              length;

      length = fabs (x_max - x_min);
      bbox_x_min -= length * margin;
      bbox_x_max += length * margin;

      length = fabs (y_max - y_min);
      bbox_y_min -= length * margin;
      bbox_y_max += length * margin;

      length = fabs (z_max - z_min);
      bbox_z_min -= length * margin;
      bbox_z_max += length * margin;
    }
  }

  /**
   * Gets (x,y,z) coordinates of one point.
   */
  inline void
  get_point_coordinates (double *coordinates, const size_t idx)
  const
  {
    coordinates[0] = point[idx].x;
    coordinates[1] = point[idx].y;
    coordinates[2] = point[idx].z;
  }

  /**
   * Sets (x,y,z) coordinates of one point.
   */
  inline void
  set_point_coordinates (const size_t idx, const double x, const double y,
                         const double z)
  {
    point[idx].x = x;
    point[idx].y = y;
    point[idx].z = z;
  }

  /**
   * Gets value of one point.
   */
  inline double
  get_point_value (const size_t idx)
  const
  {
    return point[idx].value;
  }

  /**
   * Sets value of one point.
   */
  inline void
  set_point_value (const size_t idx, const double value)
  {
    point[idx].value = value;
  }

  /**
   * Gets label of one point.
   */
  inline int
  get_point_label (const size_t idx)
  const
  {
    return point[idx].label;
  }

  /**
   * Sets label of one point.
   */
  inline void
  set_point_label (const size_t idx, const int label)
  {
    point[idx].label = label;
  }

  /**
   * Sets (x,y,z) coordinates and data of one point.
   */
  inline void
  set_point (const size_t idx, const double x, const double y, const double z,
             const double value, const int label)
  {
    set_point_coordinates (idx, x, y, z);
    set_point_value (idx, value);
    set_point_label (idx, label);
  }

  /**
   * Sets coordinates of all points.
   */
  void
  set_point_coordinates_all (const double *xyz, const size_t n_points);

  /**
   * Sets values of all points.
   */
  void
  set_point_value_all (const double *value);

  /**
   * Sets labels of all points.
   */
  void
  set_point_label_all (const int *label);

  /**
   * Sets coordinates of all points.
   */
  void
  set_point_all (const double *xyz, const double *value, const int *label,
                 const size_t n_points);

  /**
   * Adds a new point to the end of the cloud array.
   */
  inline void
  push_back_point (rhea_pointcloud_point pt)
  {
    point.push_back (pt);
  }

  /**
   * Resizes the point cloud.
   */
  void
  resize (const size_t n_points)
  {
    point.resize (n_points);
  }

  /**
   * Reserves a size for point cloud.
   */
  void
  reserve (const size_t n_points)
  {
    point.reserve (n_points);
  }

  /**
   * Gets size of point cloud.
   */
  inline size_t
  size ()
  const
  {
    return point.size ();
  }

  /**
   * Reads points from a text file and fills them into the cloud.
   */
  void
  read_points_text_file (const char *filename, const size_t estimated_n_points);

  /**
   * Prints all points in the cloud.
   */
  void
  print_points ()
  const;

  /* -------------------- Nanoflann Interface -------------------- */

  /**
   * Returns the number of data points.
   */
  inline size_t
  kdtree_get_point_count ()
  const
  {
    return point.size ();
  }

  /**
   * Computes the square distance between the vector `pt[0:size-1]` and the
   * point associated to index `idx`.
   */
  inline double
  kdtree_distance (const double *pt, const size_t idx, size_t size)
  const
  {
    const double      d0 = pt[0] - point[idx].x;
    const double      d1 = pt[1] - point[idx].y;
    const double      d2 = pt[2] - point[idx].z;

    return d0*d0 + d1*d1 + d2*d2;
  }

  /**
   * Returns the dim'th component of the point with index `idx`.
   * Note: Since this is inlined and the `dim` argument is typically an
   *       immediate value, the "if/else's" are actually solved at compile time.
   */
  inline double
  kdtree_get_pt (const size_t idx, int dim)
  const
  {
    if (dim == 0) {
      return point[idx].x;
    }
    else if (dim == 1) {
      return point[idx].y;
    }
    else {
      return point[idx].z;
    }
  }

  /**
   * Gets the bounding box for the point cloud.
   */
  template <class BBOX>
  bool
  kdtree_get_bbox (BBOX &bbox)
  const
  {
    /* set limits for each dimension */
    bbox[0].low  = bbox_x_min;
    bbox[0].high = bbox_x_max;

    bbox[1].low  = bbox_y_min;
    bbox[1].high = bbox_y_max;

    bbox[2].low  = bbox_z_min;
    bbox[2].high = bbox_z_max;

    return true;
  }
}; /* class rhea_pointcloud_Cloud */

/******************************************************************************
 * KD-Tree
 *****************************************************************************/

/* max number of nodes that a leaf of the kd-tree can have */
#define RHEA_POINTCLOUD_KDTREE_LEAF_MAX_SIZE 10

/**
 * rhea_pointcloud_KDTree: Kd-tree topology for fast computation of the shortest
 * distance to any point in the cloud.
 */
class rhea_pointcloud_KDTree
{
  private:

  /* point cloud */
  const rhea_pointcloud_Cloud  &cloud;

  /* kd-tree */
  nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<double, rhea_pointcloud_Cloud>,
      rhea_pointcloud_Cloud, RHEA_POINTCLOUD_DIM
  >  *tree;

  public:

  /**
   * Constructs a new object.
   */
  rhea_pointcloud_KDTree (const rhea_pointcloud_Cloud &cloud_in) :
    cloud (cloud_in)
  {
    const size_t        leaf_max_size = RHEA_POINTCLOUD_KDTREE_LEAF_MAX_SIZE;

    /* create new tree */
    tree = new nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, rhea_pointcloud_Cloud>,
        rhea_pointcloud_Cloud, RHEA_POINTCLOUD_DIM
    > (
        RHEA_POINTCLOUD_DIM, cloud,
        nanoflann::KDTreeSingleIndexAdaptorParams (leaf_max_size)
    );

    /* setup tree */
    setup ();
  }

  /**
   * Destroys an object.
   */
  ~rhea_pointcloud_KDTree ()
  {
    delete tree;
  }

  /**
   * Builds kd-tree index.
   */
  void
  setup ()
  {
    tree->buildIndex ();
  }

  /**
   * Finds point of the cloud that is nearest to a target point.
   * Returns Euclidean distance to nearest point, or NAN if no point was found.
   */
  double
  find_nearest (double *nearest_coordinates, double *nearest_value,
                int *nearest_label, const double *target_coordinates)
  const
  {
    size_t              nearest_idx;
    double              nearest_dist_sq;
    nanoflann::KNNResultSet<double> result_set (1 /* #nearest */);
    size_t              n_results;

    /* perform nearest neighbor search */
    n_results = tree->knnSearch (target_coordinates, 1 /* #nearest */,
                                 &nearest_idx, &nearest_dist_sq);

    /* set output data and return Euclidean distance */
    if (n_results) {
      if (nearest_coordinates != NULL) {
        cloud.get_point_coordinates (nearest_coordinates, nearest_idx);
      }
      if (nearest_value != NULL) {
        *nearest_value = cloud.get_point_value (nearest_idx);
      }
      if (nearest_label != NULL) {
        *nearest_label = cloud.get_point_label (nearest_idx);
      }
      return sqrt (nearest_dist_sq);
    }
    else {
      if (nearest_coordinates != NULL) {
        nearest_coordinates[0] = NAN;
        nearest_coordinates[1] = NAN;
        nearest_coordinates[2] = NAN;
      }
      if (nearest_value != NULL) {
        *nearest_value = NAN;
      }
      if (nearest_label != NULL) {
        *nearest_label = RHEA_POINTCLOUD_LABEL_UNDEF;
      }
      return NAN;
    }
  }

  /**
   * Finds multiple points of the cloud that are nearest to a target point.
   * Returns number of nearest points found.
   */
  size_t
  find_n_nearest (double *nearest_dist, double *nearest_coordinates,
                  double *nearest_value, int *nearest_label,
                  const size_t n_nearest, const double *target_coordinates)
  const
  {
    std::vector<size_t> nearest_idx (n_nearest);
    std::vector<double> nearest_dist_sq (n_nearest);
    size_t              n_results;

    /* perform nearest neighbor search */
    n_results = tree->knnSearch (target_coordinates, n_nearest,
                                 &nearest_idx[0], &nearest_dist_sq[0]);

    /* set output data */
    if (nearest_dist != NULL) {
      for (size_t k = 0; k < n_results; k++) {
        nearest_dist[k] = sqrt (nearest_dist_sq[k]);
      }
      for (size_t k = n_results; k < n_nearest; k++) {
        nearest_dist[k] = NAN;
      }
    }
    if (nearest_coordinates != NULL) {
      for (size_t k = 0; k < n_results; k++) {
        cloud.get_point_coordinates (&nearest_coordinates[3*k], nearest_idx[k]);
      }
      for (size_t k = n_results; k < n_nearest; k++) {
        nearest_coordinates[3*k    ] = NAN;
        nearest_coordinates[3*k + 1] = NAN;
        nearest_coordinates[3*k + 2] = NAN;
      }
    }
    if (nearest_value != NULL) {
      for (size_t k = 0; k < n_results; k++) {
        nearest_value[k] = cloud.get_point_value (nearest_idx[k]);
      }
      for (size_t k = n_results; k < n_nearest; k++) {
        nearest_value[k] = NAN;
      }
    }
    if (nearest_label != NULL) {
      for (size_t k = 0; k < n_results; k++) {
        nearest_label[k] = cloud.get_point_label (nearest_idx[k]);
      }
      for (size_t k = n_results; k < n_nearest; k++) {
        nearest_label[k] = RHEA_POINTCLOUD_LABEL_UNDEF;
      }
    }

    /* return number of nearest points found */
    return n_results;
  }

  /**
   * Finds the shortest distance from the target point to any point of the
   * cloud.
   */
  inline double
  find_shortest_distance (const double *target_pt)
  const
  {
    size_t              return_index;
    nanoflann::KNNResultSet<double>  result_set (1);
    double              dist_sq;

    /* perform nearest neighbor search */
    result_set.init (&return_index, &dist_sq);
    tree->findNeighbors (result_set, target_pt, nanoflann::SearchParams ());

    /* return Euclidean distance */
    return sqrt (dist_sq);
  }

  /**
   * Finds the shortest distance from a set of target points to any point of
   * the cloud.
   */
  inline void
  find_shortest_distance_multi (double *dist,
                                const double *target_x,
                                const double *target_y,
                                const double *target_z,
                                const unsigned int n_target_points)
  const
  {
    size_t              return_index;
    nanoflann::KNNResultSet<double>  result_set (1);
    nanoflann::SearchParams  search_params;
    double              pt[3];
    double              dist_sq;

    /* perform nearest neighbor search for all target points */
    for (size_t i = 0; i < n_target_points; i++) {
      pt[0] = target_x[i];
      pt[1] = target_y[i];
      pt[2] = target_z[i];

      /* initialize result set object */
      result_set.init (&return_index, &dist_sq);

      /* find nearest neighbor */
      tree->findNeighbors (result_set, pt, search_params);

      /* set Euclidean distance */
      dist[i] = sqrt (dist_sq);
    }
  }
}; /* class rhea_pointcloud_KDTree */

#endif /* RHEA_POINTCLOUD_HPP */
