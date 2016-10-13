/*
  This file is part of the ymir Library.
  ymir is a C library for modeling ice sheets

  Copyright (C) 2010, 2011 Carsten Burstedde, Toby Isaac, Georg Stadler,
                           Lucas Wilcox.

  The ymir Library is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The ymir Library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the ymir Library.  If not, see <http://www.gnu.org/licenses/>.

  ---

  This is the slabs example for global instantaneous mantle flow with plates.

*/

#ifndef  WEAKZONE_HH
#define  WEAKZONE_HH

#include <nanoflann.hh>

#include <cstdlib>
#include <iostream>
#include <fstream>

//#define NDEBUG /* if set then asserts are ignored */
#include <assert.h>

/* set dimension */
#define WEAKZONE_POINTCLOUD_DIM 3

/* set bottom and top radius of spherical shell */
#define WEAKZONE_POINTCLOUD_SHELL_RADIUS_BOTTOM 0.55
#define WEAKZONE_POINTCLOUD_SHELL_RADIUS_TOP 1.0

/* set margin factor for point cloud bounding box */
#define WEAKZONE_POINTCLOUD_MARGIN 1.01

/* set max number of nodes that a leaf of the kd-tree can have */
#define WEAKZONE_KDTREE_LEAF_MAX_SIZE 10

/*****************************************************************************
 * Weakzone Point Cloud
 ****************************************************************************/

/**
 * Structure defining a single point.
 */
struct WeakzonePoint
{
  double              x, y, z;
};

/**
 * Dataset class that defines a point cloud in a spherical shell domain.
 */
class WeakzonePointCloud
{
  private:

  /* point cloud data */
  std::vector<WeakzonePoint>  point;

  public:

  /**
   * Constructs a new WeakzonePointCloud object.
   */
  WeakzonePointCloud () {}

  /**
   * Destroys a WeakzonePointCloud object.
   */
  ~WeakzonePointCloud () {}

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
   * Returns the distance between the vector `pt[0:size-1]` and the data point
   * with index `idx` stored in the class.
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
    bb[0].low  = - WEAKZONE_POINTCLOUD_MARGIN
                 * WEAKZONE_POINTCLOUD_SHELL_RADIUS_TOP;
    bb[0].high = + WEAKZONE_POINTCLOUD_MARGIN
                 * WEAKZONE_POINTCLOUD_SHELL_RADIUS_TOP;

    /* 1th dimension limits */
    bb[1].low  = - WEAKZONE_POINTCLOUD_MARGIN
                 * WEAKZONE_POINTCLOUD_SHELL_RADIUS_TOP;
    bb[1].high = + WEAKZONE_POINTCLOUD_MARGIN
                 * WEAKZONE_POINTCLOUD_SHELL_RADIUS_TOP;

    /* 2th dimension limits */
    bb[2].low  = - WEAKZONE_POINTCLOUD_MARGIN
                 * WEAKZONE_POINTCLOUD_SHELL_RADIUS_TOP;
    bb[2].high = + WEAKZONE_POINTCLOUD_MARGIN
                 * WEAKZONE_POINTCLOUD_SHELL_RADIUS_TOP;

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
  push_back_point (WeakzonePoint pt)
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
   * Reads points from a text file and fills then into the cloud.
   */
  void
  read_points_text_file (const char *filename, const size_t estimated_n_points);

  /**
   * Prints all points in the cloud.
   */
  void
  print_points () const;
}; /* WeakzonePointCloud */

/*****************************************************************************
 * Weakzone KD-Tree
 ****************************************************************************/

/**
 * Kd-tree for computing the shortest distance to weak zones.
 */
class WeakzoneKDTree
{
  private:

  /* point cloud */
  const WeakzonePointCloud  &cloud;

  /* kd-tree */
  nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<double, WeakzonePointCloud>,
      WeakzonePointCloud,
      WEAKZONE_POINTCLOUD_DIM
  >  *tree;

  public:

  /**
   * Constructs a new WeakzoneKDTree object.
   */
  WeakzoneKDTree (const WeakzonePointCloud &point_cloud) : cloud(point_cloud)
  {
    /* create new tree */
    tree = new nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, WeakzonePointCloud>,
        WeakzonePointCloud,
        WEAKZONE_POINTCLOUD_DIM
    > (
        WEAKZONE_POINTCLOUD_DIM,
        cloud,
        nanoflann::KDTreeSingleIndexAdaptorParams (
            WEAKZONE_KDTREE_LEAF_MAX_SIZE
        )
    );

    /* build kd-tree index */
    tree->buildIndex ();
  }

  /**
   * Destroys a WeakzoneKDTree object.
   */
  ~WeakzoneKDTree ()
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
    double              dist_sqr;
    nanoflann::KNNResultSet<double>  result_set (1);

    /* perform nearest neighbor search */
    result_set.init (&return_index, &dist_sqr);
    tree->findNeighbors (result_set, target_pt, nanoflann::SearchParams ());

    /* return Euclidean distance */
    return sqrt (dist_sqr);
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

    /* perform nearest neighbor search for all target points */
    for (size_t i = 0; i < n_target_points; i++) {
      pt[0] = target_x[i];
      pt[1] = target_y[i];
      pt[2] = target_z[i];

      /* initialize result set object */
      result_set.init (&return_index, &dist[i]);

      /* find nearest neighbor */
      tree->findNeighbors (result_set, pt, search_params);

      /* set distance */
      dist[i] = sqrt (dist[i]);
    }
  }
};

#endif /* WEAKZONE_HH */
