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

#include <weakzone.hh>
#include <weakzone_wrapper.h>

SC_EXTERN_C_BEGIN;

/*****************************************************************************
 * Weakzone Point Cloud
 ****************************************************************************/

WeakzonePointCloud *
WeakzonePointCloud_new ()
{
  return new WeakzonePointCloud ();
}

void
WeakzonePointCloud_destroy (WeakzonePointCloud *cloud)
{
  delete cloud;
}

void
WeakzonePointCloud_set_points (WeakzonePointCloud *cloud,
                               const double *point, const size_t n_points)
{
  cloud->set_points (point, n_points);
}

/*****************************************************************************
 * Weakzone KD-Tree
 ****************************************************************************/

WeakzoneKDTree *
WeakzoneKDTree_new (const WeakzonePointCloud *point_cloud)
{
  return new WeakzoneKDTree (*point_cloud);
}

void
WeakzoneKDTree_destroy (WeakzoneKDTree *tree)
{
  delete tree;
}

double
WeakzoneKDTree_find_shortest_distance_single (WeakzoneKDTree *tree,
                                              const double *target_pt)
{
  return tree->find_shortest_distance_single (target_pt);
}

void
WeakzoneKDTree_find_shortest_distance_multi (WeakzoneKDTree *tree,
                                             double *dist,
                                             const double *target_x,
                                             const double *target_y,
                                             const double *target_z,
                                             const unsigned int
                                               n_target_points)
{
  tree->find_shortest_distance_multi (dist, target_x, target_y, target_z,
                                      n_target_points);
}

SC_EXTERN_C_END;

