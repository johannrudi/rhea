/**
 * Shared functions for rhea examples.
 *
 ******************************************************************************
 * Author:             Johann Rudi <johann@ices.utexas.edu>
 *****************************************************************************/

#ifndef EXAMPLE_SHARE_VTK_H
#define EXAMPLE_SHARE_VTK_H

#include <rhea_stokes_problem.h>

/**
 * Writes VTK of input data.
 */
void                example_share_vtk_write_input_data (
                                      const char *vtk_write_input_path,
                                      ymir_mesh_t *ymir_mesh,
                                      rhea_temperature_options_t *temp_options,
                                      rhea_viscosity_options_t *visc_options);

/**
 * Writes VTK of solution.
 */
void                example_share_vtk_write_solution (
                                      const char *vtk_write_solution_path,
                                      ymir_vec_t *sol_vel_press,
                                      rhea_stokes_problem_t *stokes_problem);

#endif /* EXAMPLE_SHARE_VTK_H */
