/*
*/

#ifndef RHEA_POINTCLOUD_HPP
#define RHEA_POINTCLOUD_HPP

#include "nanoflann/include/nanoflann.hpp"
#include <cstdlib>

/* set dimension */
#define RHEA_POINTCLOUD_DIM 3

/* set bottom and top radius of spherical shell */
#define RHEA_POINTCLOUD_SHELL_RADIUS_BOTTOM 0.55
#define RHEA_POINTCLOUD_SHELL_RADIUS_TOP 1.0

/* set margin factor for point cloud bounding box */
#define RHEA_POINTCLOUD_MARGIN 1.01

/* set max number of nodes that a leaf of the kd-tree can have */
#define RHEA_POINTCLOUD_KDTREE_LEAF_MAX_SIZE 10

/******************************************************************************
 * Point Cloud
 *****************************************************************************/

/**
 * Structure defining a single point.
 */
struct rhea_pointcloud_Point
{
  double              x, y, z;
};

/**
 * Dataset class that defines a point cloud in a spherical shell domain.
 */
class rhea_pointcloud_Cloud
{
  private:

  /* point cloud data */
  std::vector<rhea_pointcloud_Point>  point;

  public:

  /**
   * Constructs a new rhea_pointcloud_Cloud object.
   */
  rhea_pointcloud_Cloud () {}

  /**
   * Destroys a rhea_pointcloud_Cloud object.
   */
  ~rhea_pointcloud_Cloud () {}

  /* -------------------- Nanoflann Interface -------------------- */

  /**
   * Returns the number of data points.
   */
  inline size_t
  kdtree_get_point_count () const
  {
    return point.size ();
  }

  /**
   * Returns the square distance between the vector `pt[0:size-1]` and the data
   * point with index `idx` stored in the class.
   */
  inline double
  kdtree_distance (const double *pt, const size_t idx, size_t size)
  const
  {
    const double      d0 = pt[0] - point[idx].x;
    const double      d1 = pt[1] - point[idx].y;
    const double      d2 = pt[2] - point[idx].z;

    return d0 * d0 + d1 * d1 + d2 * d2;
  }

  /**
   * Returns the dim'th component of the point with index `idx`.
   * Note: Since this is inlined and the `dim` argument is typically an
   * immediate value, the "if/else's" are actually solved at compile time.
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
   * Sets the bounding box for the point cloud.
   */
  template <class BBOX>
  bool
  kdtree_get_bbox (BBOX &bb)
  const
  {
    /* 0th dimension limits */
    bb[0].low  = - RHEA_POINTCLOUD_MARGIN
                 * RHEA_POINTCLOUD_SHELL_RADIUS_TOP;
    bb[0].high = + RHEA_POINTCLOUD_MARGIN
                 * RHEA_POINTCLOUD_SHELL_RADIUS_TOP;

    /* 1st dimension limits */
    bb[1].low  = - RHEA_POINTCLOUD_MARGIN
                 * RHEA_POINTCLOUD_SHELL_RADIUS_TOP;
    bb[1].high = + RHEA_POINTCLOUD_MARGIN
                 * RHEA_POINTCLOUD_SHELL_RADIUS_TOP;

    /* 2nd dimension limits */
    bb[2].low  = - RHEA_POINTCLOUD_MARGIN
                 * RHEA_POINTCLOUD_SHELL_RADIUS_TOP;
    bb[2].high = + RHEA_POINTCLOUD_MARGIN
                 * RHEA_POINTCLOUD_SHELL_RADIUS_TOP;

    return true;
  }

  /* -------------------- Functions -------------------- */

  /**
   * Get coordinates of one point.
   */
  inline void
  get_point (double *pt, const size_t idx)
  const
  {
    pt[0] = point[idx].x;
    pt[1] = point[idx].y;
    pt[2] = point[idx].z;
  }

  /**
   * Set coordinates of one point.
   */
  inline void
  set_point (const size_t idx, const double x, const double y, const double z)
  {
    point[idx].x = x;
    point[idx].y = y;
    point[idx].z = z;
  }

  /**
   * Set coordinates of all points.
   */
  void
  set_points (const double *coord, const size_t n_points);

  /**
   * Add new point to the end of the cloud array.
   */
  inline void
  push_back_point (rhea_pointcloud_Point pt)
  {
    point.push_back (pt);
  }

  /**
   * Resize point cloud.
   */
  void
  resize (const size_t n_points)
  {
    point.resize (n_points);
  }

  /**
   * Reserve a size for point cloud.
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
  size () const
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
  print_points () const;
}; /* rhea_pointcloud_Cloud */

/******************************************************************************
 * KD-Tree
 *****************************************************************************/

/**
 * Kd-tree for computing the shortest distance to any point in the cloud.
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
   * Constructs a new rhea_pointcloud_KDTree object.
   */
  rhea_pointcloud_KDTree (const rhea_pointcloud_Cloud &cloud_in) :
    cloud (cloud_in)
  {
    /* create new tree */
    tree = new nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, rhea_pointcloud_Cloud>,
        rhea_pointcloud_Cloud, RHEA_POINTCLOUD_DIM
    > (
        RHEA_POINTCLOUD_DIM, cloud,
        nanoflann::KDTreeSingleIndexAdaptorParams (
            RHEA_POINTCLOUD_KDTREE_LEAF_MAX_SIZE)
    );

    /* build kd-tree index */
    tree->buildIndex ();
  }

  /**
   * Destroys a rhea_pointcloud_KDTree object.
   */
  ~rhea_pointcloud_KDTree ()
  {
    delete tree;
  }

  /**
   * Finds the shortest distance from the target point to a point of the cloud.
   */
  inline double
  find_shortest_distance_single (const double *target_pt)
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
   * Finds the shortest distance from a set of target points to a point of
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
};

#endif /* RHEA_POINTCLOUD_HPP */
