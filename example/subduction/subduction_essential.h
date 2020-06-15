#ifndef SUBDUCTION_ESSENTIAL_H
#define SUBDUCTION_ESSENTIAL_H

#include <example_share_mesh.h>
#include <example_share_stokes.h>
#include <rhea_stokes_problem_amr.h>
#include <subduction.h>
#include <subduction_physics.h>
#include <subduction_test.h>
#include <subduction_TI.h>
#include <subduction_vtk.h>
#include <subduction_adjoint.h>

void
subd_setup_mesh (p4est_t **p4est,
                    ymir_mesh_t **ymir_mesh,
                    ymir_pressure_elem_t **press_elem,
                    MPI_Comm mpicomm,
                    rhea_domain_options_t *domain_options,
                    rhea_topography_options_t *topo_options,
                    rhea_discretization_options_t *discr_options,
                    subd_options_t *subd_options);

void
subd_setup_stokes (rhea_stokes_problem_t **stokes_problem,
                    p4est_t     *p4est,
                    ymir_mesh_t **ymir_mesh,
                    ymir_pressure_elem_t **press_elem,
                    rhea_discretization_options_t *discr_options,
                    rhea_domain_options_t *domain_options,
                    rhea_temperature_options_t *temp_options,
                    rhea_weakzone_options_t *weak_options,
                    rhea_viscosity_options_t *visc_options,
					rhea_composition_options_t *comp_options,
                    subd_options_t *subd_options,
                    const char *vtk_write_input_path,
                    const char *bin_solver_path,
                    const char *vtk_solver_path);


void
subd_setup_clear_all (rhea_stokes_problem_t *stokes_problem,
                         p4est_t *p4est,
                         ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         rhea_temperature_options_t *temp_options,
                         rhea_viscosity_options_t *visc_options,
                         rhea_weakzone_options_t *weak_options,
                         rhea_topography_options_t *topo_options,
                         rhea_plate_options_t *plate_options,
                         rhea_discretization_options_t *discr_options);

void
subd_run_solver (ymir_vec_t *sol_vel_press,
                    ymir_mesh_t *ymir_mesh,
                    ymir_pressure_elem_t *press_elem,
                    rhea_stokes_problem_t *stokes_problem,
                    const int nonzero_initial_guess,
                    const int iter_max, const double rel_tol);

#endif /*SUBDUCTION_ESSENTIAL_H*/
