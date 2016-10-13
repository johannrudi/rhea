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

void
WeakzonePointCloud::set_points (const double *coord, const size_t n_points)
{
  /* resize cloud to number of points */
  resize (n_points);

  for (size_t i = 0; i < n_points; i++) { /* loop over all points */
    assert (
        (- WEAKZONE_POINTCLOUD_MARGIN * WEAKZONE_POINTCLOUD_SHELL_RADIUS_TOP)
        <= coord[3 * i]
        &&
        coord[3 * i]
        <= (WEAKZONE_POINTCLOUD_MARGIN * WEAKZONE_POINTCLOUD_SHELL_RADIUS_TOP)
    );
    assert (
        (- WEAKZONE_POINTCLOUD_MARGIN * WEAKZONE_POINTCLOUD_SHELL_RADIUS_TOP)
        <= coord[3 * i + 1]
        &&
        coord[3 * i + 1]
        <= (WEAKZONE_POINTCLOUD_MARGIN * WEAKZONE_POINTCLOUD_SHELL_RADIUS_TOP)
    );
    assert (
        (- WEAKZONE_POINTCLOUD_MARGIN * WEAKZONE_POINTCLOUD_SHELL_RADIUS_TOP)
        <= coord[3 * i + 2]
        &&
        coord[3 * i + 2]
        <= (WEAKZONE_POINTCLOUD_MARGIN * WEAKZONE_POINTCLOUD_SHELL_RADIUS_TOP)
    );

    /* set coordinates of a point */
    set_point (i, coord[3 * i], coord[3 * i + 1], coord[3 * i + 2]);
  }
}

void
WeakzonePointCloud::read_points_text_file (const char *filename,
                                           const size_t estimated_n_points = 10)
{
  const char          this_fn_name[] =
                        "[WeakzonePointCloud::read_points_text_file]";
  size_t              n_points = 0;
  std::ifstream       file;
  WeakzonePoint       pt;

  /* set estimated cloud size */
  reserve (estimated_n_points);

  /* open file */
  file.open (filename);
  if (!file.is_open ()) { /* if not opened */
    std::cout << this_fn_name << " Error: " << "Could not open file "
              << filename << std::endl;
    return;
  }

  /* read points */
  while (!file.eof ()) {
    /* read x,y,z coordinates */
    if ( (file >> pt.x) && (file >> pt.y) && (file >> pt.z) ) {
      /* store point in point cloud */
      push_back_point (pt);

      /* count points */
      n_points++;
    }
  }

  /* close file */
  file.close ();
}

void
WeakzonePointCloud::print_points () const
{
  for (size_t i = 0; i < size (); i++) {
    WeakzonePoint     pt = point[i];

    printf ("point %09lu: %+1.8f, %+1.8f, %+1.8f\n", (long unsigned int) i,
            pt.x, pt.y, pt.z);
  }
}

