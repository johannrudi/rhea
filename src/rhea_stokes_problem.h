/** RHEA_STOKES_PROBLEM
 *
 * Provides the Stokes solver.
 */

#ifndef RHEA_STOKES_PROBLEM_H
#define RHEA_STOKES_PROBLEM_H

#include <rhea_domain.h>
#include <rhea_discretization.h>
#include <rhea_temperature.h>
#include <rhea_composition.h>
#include <rhea_plate.h>
#include <rhea_weakzone.h>
#include <rhea_viscosity.h>
#include <rhea_velocity.h>
#include <rhea_pressure.h>
#include <rhea_velocity_pressure.h>
#include <ymir_stokes_op.h>

/* types of linearization */
typedef enum
{
  RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_PICARD,
  RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_REGULAR,
  RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL,
  RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_PRIMALDUAL_SYMM,
  RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_DEV1, //TODO change name
  RHEA_STOKES_PROBLEM_NONLINEAR_LINEARIZATION_NEWTON_DEV2  //TODO change name
}
rhea_stokes_problem_nonlinear_linearization_t;

/******************************************************************************
 * Options
 *****************************************************************************/

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_stokes_problem_add_options (ymir_options_t * opt_sup);

/**
 * Processes options and stores them.
 */
void                rhea_stokes_problem_process_options ();

/******************************************************************************
 * Monitoring
 *****************************************************************************/

/**
 * Initializes performance counters.
 */
void                rhea_stokes_problem_perfmon_init (const int activate,
                                                      const int skip_if_active);

/**
 * Prints statistics collected by performance monitors.
 */
void                rhea_stokes_problem_perfmon_print (sc_MPI_Comm mpicomm,
                                                       const int print_wtime,
                                                       const int print_n_calls,
                                                       const int print_flops);

/******************************************************************************
 * Stokes Problem
 *****************************************************************************/

/* Stokes problem (opaque) */
typedef struct rhea_stokes_problem rhea_stokes_problem_t;

/**
 * Creates/destroys a Stokes problem.
 */
rhea_stokes_problem_t *rhea_stokes_problem_new (
                                    ymir_mesh_t *ymir_mesh,
                                    ymir_pressure_elem_t *press_elem,
                                    ymir_vec_t *temperature,
                                    ymir_vec_t *composition,
                                    rhea_domain_options_t *domain_options,
                                    rhea_temperature_options_t *temp_options,
                                    rhea_composition_options_t *comp_options,
                                    rhea_weakzone_options_t *weak_options,
                                    rhea_viscosity_options_t *visc_options);

void                rhea_stokes_problem_destroy (
                                    rhea_stokes_problem_t *stokes_problem);

/**
 * Sets up/clears objects that have dependencies on the mesh.
 */
void                rhea_stokes_problem_create_mesh_dependencies (
                                    rhea_stokes_problem_t *stokes_problem,
                                    ymir_mesh_t *ymir_mesh,
                                    ymir_pressure_elem_t *press_elem);

void                rhea_stokes_problem_clear_mesh_dependencies (
                                    rhea_stokes_problem_t *stokes_problem);

/******************************************************************************
 * Stokes Solver
 *****************************************************************************/

/**
 * Sets up the solver of a Stokes problem.
 */
void                rhea_stokes_problem_setup_solver (
                                    rhea_stokes_problem_t *stokes_problem);

/**
 * Solves a Stokes problem.
 */
int                 rhea_stokes_problem_solve (
                                    ymir_vec_t **sol_vel_press,
                                    const int nonzero_initial_guess,
                                    const int iter_max,
                                    const double rtol,
                                    rhea_stokes_problem_t *stokes_problem);

int                 rhea_stokes_problem_has_converged (
                                    const int stop_reason,
                                    rhea_stokes_problem_t *stokes_problem);

/******************************************************************************
 * Stokes Solver, Extended
 *****************************************************************************/

/**
 * Computes the viscous stress coefficient of a Stokes problem.
 */
void                rhea_stokes_problem_compute_coefficient (
                                    rhea_stokes_problem_t *stokes_problem,
                                    ymir_vec_t *velocity_pressure,
                                    const int nonlinear_init);

/**
 * Computes and updates the viscous stress coefficient of a Stokes problem.
 */
void                rhea_stokes_problem_compute_and_update_coefficient (
                                    rhea_stokes_problem_t *stokes_problem,
                                    ymir_vec_t *velocity_pressure,
                                    const int nonlinear_init);

/**
 * Sets up the solver of a Stokes problem.
 */
void                rhea_stokes_problem_setup_solver_ext (
                                    rhea_stokes_problem_t *stokes_problem,
                                    const int num_krylov_solvers);

/**
 * Updates the solver of a Stokes problem.
 */
void                rhea_stokes_problem_update_solver (
                                    rhea_stokes_problem_t *stokes_problem,
                                    const int update_coeff,
                                    ymir_vec_t *vel_press,
                                    const int override_rhs,
                                    ymir_vec_t *rhs_vel_press,
                                    ymir_vec_t *rhs_vel,
                                    ymir_vec_t *vel_nonzero_dirichlet,
                                    ymir_vec_t **rhs_vel_face_nonzero_neumann);

/**
 * Solves a Stokes problem.
 */
int                 rhea_stokes_problem_solve_ext (
                                    ymir_vec_t **sol_vel_press,
                                    const int nonzero_initial_guess,
                                    const int iter_max,
                                    const double rtol,
                                    rhea_stokes_problem_t *stokes_problem,
                                    const int resume,
                                    const int force_linear_solve,
                                    const int krylov_solver_idx,
                                    int *num_iterations,
                                    double *residual_reduction,
                                    int *mesh_modified_by_solver);

int                 rhea_stokes_problem_has_converged_ext (
                                    const int stop_reason,
                                    rhea_stokes_problem_t *stokes_problem,
                                    const int force_linear_solve);

/**
 * Sets data to enable solver AMR/grid continuation.
 */
void                rhea_stokes_problem_set_solver_amr (
                                rhea_stokes_problem_t *stokes_problem,
                                p4est_t *p4est,
                                rhea_discretization_options_t *discr_options);

/**
 * Sets output path for binary output the iterations of a nonlinear solve.
 */
void                rhea_stokes_problem_set_solver_bin_output (
                                    rhea_stokes_problem_t *stokes_problem,
                                    char *bin_path);

/**
 * Sets output path for vtk output the iterations of a nonlinear solve.
 */
void                rhea_stokes_problem_set_solver_vtk_output (
                                    rhea_stokes_problem_t *stokes_problem,
                                    char *vtk_path);

/******************************************************************************
 * Data Access
 *****************************************************************************/

/**
 * Accesses data of a Stokes problem.
 */
ymir_mesh_t        *rhea_stokes_problem_get_ymir_mesh (
                                    rhea_stokes_problem_t *stokes_problem);
ymir_pressure_elem_t *rhea_stokes_problem_get_press_elem (
                                    rhea_stokes_problem_t *stokes_problem);

void                rhea_stokes_problem_set_temperature (
                                    rhea_stokes_problem_t *stokes_problem,
                                    ymir_vec_t *temperature);
ymir_vec_t         *rhea_stokes_problem_get_temperature (
                                    rhea_stokes_problem_t *stokes_problem);
void                rhea_stokes_problem_remove_temperature (
                                    rhea_stokes_problem_t *stokes_problem);

void                rhea_stokes_problem_set_composition (
                                    rhea_stokes_problem_t *stokes_problem,
                                    ymir_vec_t *composition_density,
                                    ymir_vec_t *composition_viscosity);
ymir_vec_t         *rhea_stokes_problem_get_composition_density (
                                    rhea_stokes_problem_t *stokes_problem);
void                rhea_stokes_problem_remove_composition (
                                    rhea_stokes_problem_t *stokes_problem);

void                rhea_stokes_problem_set_velocity_pressure (
                                    rhea_stokes_problem_t *stokes_problem,
                                    ymir_vec_t *velocity_pressure);
ymir_vec_t         *rhea_stokes_problem_get_velocity_pressure (
                                    rhea_stokes_problem_t *stokes_problem);
void                rhea_stokes_problem_remove_velocity_pressure (
                                    rhea_stokes_problem_t *stokes_problem);

int                 rhea_stokes_problem_is_nonlinear (
                                    rhea_stokes_problem_t *stokes_problem);

rhea_domain_options_t      *rhea_stokes_problem_get_domain_options (
                                    rhea_stokes_problem_t *stokes_problem);
void                        rhea_stokes_problem_set_domain_options (
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_domain_options_t *domain_options);
rhea_temperature_options_t *rhea_stokes_problem_get_temperature_options (
                                    rhea_stokes_problem_t *stokes_problem);
void                        rhea_stokes_problem_set_temperature_options (
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_temperature_options_t *temp_options);
rhea_composition_options_t *rhea_stokes_problem_get_composition_options (
                                    rhea_stokes_problem_t *stokes_problem);
rhea_plate_options_t       *rhea_stokes_problem_get_plate_options (
                                    rhea_stokes_problem_t *stokes_problem);
void                        rhea_stokes_problem_set_plate_options (
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_plate_options_t *plate_options);
rhea_weakzone_options_t    *rhea_stokes_problem_get_weakzone_options (
                                    rhea_stokes_problem_t *stokes_problem);
void                        rhea_stokes_problem_set_weakzone_options (
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_weakzone_options_t *weak_options);
rhea_viscosity_options_t   *rhea_stokes_problem_get_viscosity_options (
                                    rhea_stokes_problem_t *stokes_problem);
void                        rhea_stokes_problem_set_viscosity_options (
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_viscosity_options_t *visc_options);

void                rhea_stokes_problem_copy_viscosity (
                                    ymir_vec_t *viscosity,
                                    rhea_stokes_problem_t *stokes_problem);
void                rhea_stokes_problem_copy_marker (
                                    ymir_vec_t *marker,
                                    rhea_stokes_problem_t *stokes_problem);

ymir_vec_t         *rhea_stokes_problem_get_weakzone (
                                    rhea_stokes_problem_t *stokes_problem);
ymir_vec_t         *rhea_stokes_problem_get_rhs_vel_press (
                                    rhea_stokes_problem_t *stokes_problem);//XI
void                rhea_stokes_problem_set_rhs_vel_press (
                                    rhea_stokes_problem_t *stokes_problem,
                                    ymir_vec_t *rhs_vel_press);//XI
ymir_vec_t         *rhea_stokes_problem_get_rhs_vel (
                                    rhea_stokes_problem_t *stokes_problem);
ymir_vec_t         *rhea_stokes_problem_get_rhs_vel_nonzero_dirichlet (
                                    rhea_stokes_problem_t *stokes_problem);
ymir_vec_t        **rhea_stokes_problem_get_rhs_vel_nonzero_neumann (
                                    rhea_stokes_problem_t *stokes_problem);
ymir_vec_t         *rhea_stokes_problem_get_rhs_vel_nonzero_neumann_surface (
                                    rhea_stokes_problem_t *stokes_problem);

void                rhea_stokes_problem_set_weakzone_compute_fn (
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_weakzone_compute_fn_t fn,
                                    void *data);
void                rhea_stokes_problem_set_viscosity_compute_fn (
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_viscosity_compute_fn_t fn,
                                    void *data);
void                rhea_stokes_problem_set_rhs_vel_compute_fn (
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_velocity_rhs_compute_fn_t fn,
                                    void *data);
void                rhea_stokes_problem_set_rhs_vel_nonzero_dir_compute_fn (
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_velocity_rhs_nz_dir_compute_fn_t fn,
                                    void *data);
void                rhea_stokes_problem_set_rhs_vel_nonzero_neu_compute_fn (
                                    rhea_stokes_problem_t *stokes_problem,
                                    rhea_velocity_rhs_nz_neu_compute_fn_t fn,
                                    void *data);

void                rhea_stokes_problem_weakzone_compute (
                                    ymir_vec_t *weakzone,
                                    rhea_stokes_problem_t *stokes_problem);
void                rhea_stokes_problem_viscosity_compute (
                                    ymir_vec_t *viscosity,
                                    ymir_vec_t *proj_scal,
                                    ymir_vec_t *marker,
                                    ymir_vec_t *temperature,
                                    ymir_vec_t *weakzone,
                                    ymir_vec_t *velocity,
                                    rhea_stokes_problem_t *stokes_problem);

ymir_stokes_op_t   *rhea_stokes_problem_get_stokes_op (
                                    rhea_stokes_problem_t *stokes_problem);

int                 rhea_stokes_problem_get_num_krylov_solvers (
                                    rhea_stokes_problem_t *stokes_problem);

/******************************************************************************
 * I/O
 *****************************************************************************/

/**
 * Writes a Stokes problem to disk.
 */
int                 rhea_stokes_problem_write (
                                    char *base_path_bin,
                                    rhea_stokes_problem_t *stokes_problem);

/******************************************************************************
 * Vector Operations
 *****************************************************************************/

/**
 * Sets velocity components on the boundary, which are constrained by Dirichlet
 * boundary conditions.
 */
int                 rhea_stokes_problem_velocity_set_boundary (
                                        ymir_vec_t *velocity,
                                        ymir_vec_t *velocity_on_boundary,
                                        rhea_stokes_problem_t *stokes_problem);

int                 rhea_stokes_problem_velocity_set_boundary_zero (
                                        ymir_vec_t *velocity,
                                        rhea_stokes_problem_t *stokes_problem);

int                 rhea_stokes_problem_velocity_enforce_boundary_conditions (
                                        ymir_vec_t *velocity,
                                        rhea_stokes_problem_t *stokes_problem);

int                 rhea_stokes_problem_enforce_boundary_conditions (
                                        ymir_vec_t *velocity_pressure,
                                        rhea_stokes_problem_t *stokes_problem);

/**
 * Computes mean rotation of the velocity.
 */
int                 rhea_stokes_problem_velocity_compute_mean_rotation (
                                    double mean_rot_axis[3],
                                    ymir_vec_t *velocity,
                                    rhea_stokes_problem_t *stokes_problem);

/**
 * Projects out mean rotation of the velocity.
 */
int                 rhea_stokes_problem_velocity_project_out_mean_rotation (
                                        ymir_vec_t *velocity,
                                        const int residual_space,
                                        rhea_stokes_problem_t *stokes_problem);

/**
 * Projects out null spaces.
 */
int                 rhea_stokes_problem_project_out_nullspace (
                                        ymir_vec_t *vel_press,
                                        rhea_stokes_problem_t *stokes_problem);

/**
 * Computes stress tensor.
 */
int                 rhea_stokes_problem_stress_compute (
                                    ymir_vec_t *stress,
                                    ymir_vec_t *vel_press,
                                    rhea_stokes_problem_t *stokes_problem,
                                    ymir_stress_op_t *override_stress_op,
                                    int linearized_visc_coeff,
                                    const int skip_pressure);

/**
 * Computes normal stress at the surface of the domain.
 */
int                 rhea_stokes_problem_stress_compute_normal_at_surface (
                                    ymir_vec_t *stress_norm_surf,
                                    ymir_vec_t *vel_press,
                                    rhea_stokes_problem_t *stokes_problem);

/**
 * Computes divergence of strain rate multiplied by viscous stress coefficient.
 */
int                 rhea_stokes_problem_stress_div_compute (
                                    ymir_vec_t *velocity,
                                    ymir_vec_t *stress_like_tensor,
                                    rhea_stokes_problem_t *stokes_problem,
                                    ymir_stress_op_t *override_stress_op,
                                    int linearized_visc_coeff);

#endif /* RHEA_STOKES_PROBLEM_H */
