/* EXAMPLE_SHARE_IO  Shared functions for rhea examples.
 *
 * Author:             Johann Rudi <johann@ices.utexas.edu>
 */

#ifndef EXAMPLE_SHARE_IO_H
#define EXAMPLE_SHARE_IO_H

#include <rhea_stokes_problem.h>

/**
 * Writes velocity and normal stress at surface to text file.
 */
void                example_share_io_write_solution_surf_txt (
                                      const char *file_path_txt,
                                      rhea_domain_coordinate_type_t coord_type,
                                      ymir_vec_t *sol_vel_press,
                                      rhea_stokes_problem_t *stokes_problem);

#endif /* EXAMPLE_SHARE_IO_H */
