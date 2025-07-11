/*
*/

#include <rhea_pointcloud.hpp>
#include <rhea_base.h>

#include <iostream>
#include <fstream>

void
rhea_pointcloud_Cloud::set_point_coordinates_all (const double *xyz,
                                                  const size_t n_points)
{
  /* resize cloud to number of points */
  resize (n_points);

  for (size_t i = 0; i < n_points; i++) { /* loop over all points */
    const double        x = xyz[3*i    ];
    const double        y = xyz[3*i + 1];
    const double        z = xyz[3*i + 2];

    /* check input */
    RHEA_ASSERT (bbox_x_min <= x && x <= bbox_x_max);
    RHEA_ASSERT (bbox_y_min <= y && y <= bbox_y_max);
    RHEA_ASSERT (bbox_z_min <= z && z <= bbox_z_max);

    /* set coordinates of this point */
    set_point_coordinates (i, x, y, z);
  }
}

void
rhea_pointcloud_Cloud::set_point_value_all (const double *value)
{
  for (size_t i = 0; i < size (); i++) { /* loop over all points */
    set_point_value (i, value[i]);
  }
}

void
rhea_pointcloud_Cloud::set_point_label_all (const int *label)
{
  for (size_t i = 0; i < size (); i++) { /* loop over all points */
    set_point_label (i, label[i]);
  }
}

void
rhea_pointcloud_Cloud::set_point_all (const double *xyz, const double *value,
                                      const int *label, const size_t n_points)
{
  set_point_coordinates_all (xyz, n_points);
  set_point_value_all (value);
  set_point_label_all (label);
}

void
rhea_pointcloud_Cloud::print_points ()
const
{
  for (size_t i = 0; i < size (); i++) {
    rhea_pointcloud_point pt = point[i];

    printf ("point %09lu: %+.8f,%+.8f,%+.8f, value %g, label %i\n",
            (long unsigned int) i, pt.x, pt.y, pt.z, pt.value, pt.label);
  }
}
