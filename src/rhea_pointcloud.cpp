/*
*/

#include <rhea_pointcloud.hpp>
#include <rhea_base.h>

#include <iostream>
#include <fstream>
#ifdef RHEA_ENABLE_DEBUG
#include <assert.h>
#endif

void
rhea_pointcloud_Cloud::set_points (const double *coordinates,
                                   const size_t n_points)
{
#ifdef RHEA_ENABLE_DEBUG
  const double        xyz_min = - RHEA_POINTCLOUD_MARGIN *
                                  RHEA_POINTCLOUD_SHELL_RADIUS_TOP;
  const double        xyz_max = + RHEA_POINTCLOUD_MARGIN *
                                  RHEA_POINTCLOUD_SHELL_RADIUS_TOP;
#endif

  /* resize cloud to number of points */
  resize (n_points);

  for (size_t i = 0; i < n_points; i++) { /* loop over all points */
    const double        x = coordinates[3*i    ];
    const double        y = coordinates[3*i + 1];
    const double        z = coordinates[3*i + 2];

#ifdef RHEA_ENABLE_DEBUG
    assert (xyz_min <= x && x <= xyz_max);
    assert (xyz_min <= y && y <= xyz_max);
    assert (xyz_min <= z && z <= xyz_max);
#endif

    /* set coordinates of a point */
    set_point (i, x, y, z);
  }
}

void
rhea_pointcloud_Cloud::read_points_text_file (
                                          const char *filename,
                                          const size_t estimated_n_points = 10)
{
  const char          this_fn_name[] =
                        "[rhea_pointcloud_Cloud::read_points_text_file]";
  size_t              n_points = 0;
  std::ifstream       file;
  rhea_pointcloud_Point pt;

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
rhea_pointcloud_Cloud::print_points () const
{
  for (size_t i = 0; i < size (); i++) {
    rhea_pointcloud_Point pt = point[i];

    printf ("point %09lu: %+.8f, %+.8f, %+.8f\n", (long unsigned int) i,
            pt.x, pt.y, pt.z);
  }
}
