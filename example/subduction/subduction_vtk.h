#ifndef SUBDUCTION_VTK_H
#define SUBDUCTION_VTK_H

#include <subduction_physics.h>
#include <subduction_test.h>
#include <subduction_options.h>
#include <subduction_postp.h>
#include <subduction.h>
#include <ymir_vtk.h>
#include <rhea_vtk.h>

void subduction_add_vtk_options (ymir_options_t * opt);

/* Write vtk of input data. */
void
subd_write_input (ymir_mesh_t *ymir_mesh,
                   rhea_stokes_problem_t *stokes_problem,
                   rhea_temperature_options_t *temp_options,
                   ymir_vec_t *temperature,
                   ymir_vec_t *visc_TI_svisc,
                   ymir_vec_t *visc_TI_rotate,
                   const char *vtk_write_input_path);

void
subd_write_input_basic (rhea_stokes_problem_t *stokes_problem,
                       rhea_temperature_options_t *temp_options,
                       const char *vtk_write_input_path);


void
subd_vtk_write (rhea_stokes_problem_t * stokes_problem,
                subd_options_t * subd_options);

void
subd_vtk_write_solution (ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         rhea_stokes_problem_t *stokes_problem,
                         ymir_vec_t *sol_vel_press);

void
subd_vtk_write_test (ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         rhea_stokes_problem_t *stokes_problem,
                         ymir_vec_t *sol_vel_press,
                         subd_options_t *subd_options);

void
subd_vtk_write_postp (ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         rhea_stokes_problem_t *stokes_problem,
                         ymir_vec_t *sol_vel_press,
                         subd_options_t *subd_options);

void
subd_vtk_write_stress (ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         rhea_stokes_problem_t *stokes_problem,
                         ymir_vec_t *sol_vel_press,
                         subd_options_t *subd_options);

void
subd_vtk_write_freesurface (ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         ymir_vec_t *sol_vel_press);

void
subd_vtk_write_neumann (ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         rhea_stokes_problem_t *stokes_problem,
                         ymir_vec_t *sol_vel_press);


#endif /*SUBDUCTION_VTK_H*/
