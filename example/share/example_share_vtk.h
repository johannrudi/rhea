/* EXAMPLE_SHARE_VTK  Shared functions for rhea examples.
 *
 * Author:             Johann Rudi <johann@ices.utexas.edu>
 */

#ifndef EXAMPLE_SHARE_VTK_H
#define EXAMPLE_SHARE_VTK_H

#include <rhea_stokes_problem.h>
#include <rhea_plate.h>

/**
 * Writes VTK of input data.
 */
void                example_share_vtk_write_input_data (
                                      const char *vtk_write_input_path,
                                      rhea_stokes_problem_t *stokes_problem,
                                      rhea_plate_options_t *plate_options);

/**
 * Writes VTK of solution.
 */
void                example_share_vtk_write_solution (
                                      const char *vtk_write_solution_path,
                                      ymir_vec_t *sol_vel_press,
                                      rhea_stokes_problem_t *stokes_problem);

#endif /* EXAMPLE_SHARE_VTK_H */
