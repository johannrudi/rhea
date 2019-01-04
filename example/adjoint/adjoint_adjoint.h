#ifndef ADJOINT_ADJOINT_H
#define ADJOINT_ADJOINT_H

#include <adjoint_math.h>
#include <adjoint_options.h>
#include <adjoint_io.h>
#include <adjoint_vtk.h>
#include <adjoint_essential.h>
#include <rhea.h>
#include <ymir_vec_ops.h>
#include <ymir_vec_getset.h>
#include <ymir_stress_op.h>
#include <ymir_velocity_vec.h>
#include <ymir_comm.h>
#include <ymir_stokes_vec.h>
#include <ymir_derivative_elem.h>
#include <ymir_stokes_pc.h>
#include <ymir_interp_vec.h>
#include <ymir_mass_vec.h>

typedef struct adjoint_problem
{
 ymir_vec_t           *msol;
 ymir_vec_t           *hessian;
 ymir_vec_t           *neg_grad;
 double               objective;

 ymir_vec_t           *sol_vel_press;
 ymir_vec_t           *usol;
 ymir_vec_t           *vsol;
 ymir_vec_t           **Husol;
 ymir_vec_t           **Hvsol;
 ymir_vec_t           **coeff;

 rhea_stokes_problem_t  *stokes_problem;
 p4est_t                *p4est;
 ymir_mesh_t            *ymir_mesh;
 ymir_pressure_elem_t   *press_elem;

 rhea_discretization_options_t  *discr_options;
 rhea_temperature_options_t     *temp_options;
 subd_options_t                 *subd_options;

 int                        solver_iter_max;
 double                     solver_rel_tol;

 FILE     *fp_record;
}
adjoint_problem_t;

/*
 * called in main, setup adjoint and newton problem
 */
void
adjoint_stokes_new (rhea_stokes_problem_t **stokes_problem,
                    ymir_mesh_t **ymir_mesh,
                    ymir_pressure_elem_t **press_elem,
                    rhea_domain_options_t *domain_options,
                    rhea_temperature_options_t *temp_options,
                    rhea_weakzone_options_t *weak_options,
                    rhea_viscosity_options_t *visc_options,
                    subd_options_t *subd_options);

adjoint_problem_t *
adjoint_problem_new (ymir_vec_t *solution,
                     rhea_stokes_problem_t *stokes_problem,
                     p4est_t *p4est, ymir_mesh_t *ymir_mesh,
                     ymir_pressure_elem_t *press_elem,
                     rhea_discretization_options_t *discr_options,
                     rhea_temperature_options_t *temp_options,
                     rhea_newton_options_t *newton_options,
                     subd_options_t *subd_options,
                     int     solver_iter_max,
                     double  solver_rel_tol);

void
adjoint_setup_newton (rhea_newton_problem_t **newton_problem,
                      adjoint_problem_t *adjoint_problem);


void
adjoint_destroy_newton (rhea_newton_problem_t *newton_problem);

/********************************************
 * user defined function used in rhea_newton
 ********************************************/
int
adjoint_get_field (int *field_nums, int n_components, int i);

void
adjoint_setup_stencil  (ymir_vec_t * stencil,
                        subd_options_t * subd_options,
                        subd_adjoint_stencil_field_t stencil_field);

subd_adjoint_visc_parameter_t
adjoint_get_parameter (int *parameter_nums, int n_components, int i);

void
adjoint_compute_visc_grad_m (ymir_vec_t  *visc_grad,
                         rhea_stokes_problem_t *stokes_problem,
                         subd_options_t  *subd_options,
                         subd_adjoint_visc_parameter_t parameter);

void
adjoint_run_solver (ymir_vec_t *usol, adjoint_problem_t *adjoint_problem);

void
adjoint_stokes_update_init (rhea_stokes_problem_t *stokes_problem,
                          p4est_t     *p4est,
                          ymir_mesh_t *ymir_mesh,
                          ymir_pressure_elem_t *press_elem,
                          rhea_discretization_options_t *discr_options,
                          rhea_temperature_options_t *temp_options,
                          subd_options_t *subd_options);

void
adjoint_solve_init (adjoint_problem_t *adjoint_problem);

/******************* newton function ***********************************/
void
adjoint_data_init (ymir_vec_t *solution,
                         void *data);

void
adjoint_stokes_update_forward (rhea_stokes_problem_t *stokes_problem,
                              subd_options_t *subd_options);

void
adjoint_solve_forward (adjoint_problem_t *adjoint_problem);

void
adjoint_coeff_update (adjoint_problem_t *adjoint_problem);

/******************* newton function ***********************************/
void
adjoint_update_operator_fn (ymir_vec_t *solution, void *data);

/******************* newton function ***********************************/
double
adjoint_evaluate_objective (ymir_vec_t *solution,
                          void *data);

void
adjoint_stokes_update_adjoint (ymir_vec_t *usol,
                              rhea_stokes_problem_t *stokes_problem,
                              subd_options_t *subd_options);

void
adjoint_solve_adjoint (adjoint_problem_t *adjoint_problem);

void
adjoint_update_coeff (adjoint_problem_t *adjoint_problem);

void
adjoint_compute_neg_grad (adjoint_problem_t *adjoint_problem);

/******************* newton function ***********************************/
void
adjoint_solve_negative_gradient (ymir_vec_t *neg_grad,
                                   ymir_vec_t *solution, void *data);

/******************* newton function ***********************************/
double
adjoint_compute_gradient_norm (ymir_vec_t *neg_grad,
                                    void *data, double *norm);

void
adjoint_stokes_update_hessian_forward (adjoint_problem_t *adjoint_problem, int i);

void
adjoint_solve_hessian_forward (adjoint_problem_t *adjoint_problem, int i);

void
adjoint_stokes_update_hessian_adjoint (ymir_vec_t *Husol,
                                      rhea_stokes_problem_t *stokes_problem,
                                      subd_options_t  *subd_options);

void
adjoint_solve_hessian_adjoint (adjoint_problem_t *adjoint_problem, int i);

void
adjoint_compute_hessian (ymir_vec_t *hessian,
                         adjoint_problem_t *adjoint_problem);

void
adjoint_compute_step (ymir_vec_t *step,
                         adjoint_problem_t *adjoint_problem);

/******************* newton function ***********************************/
int
adjoint_solve_hessian_system_fn (ymir_vec_t *step, ymir_vec_t *neg_grad,
                                const int lin_iter_max, const double lin_res_norm_rtol,
                                const int nonzero_initial_guess, void *data,
                                int *lin_iter_count);

/******************* newton function ***********************************/
void
adjoint_apply_hessian_fn (ymir_vec_t *hessian_vec, ymir_vec_t *dir_vec,
                         void *data);

/******************* newton function ***********************************/
void
adjoint_update_hessian_fn (ymir_vec_t *solution, ymir_vec_t *step,
                          const double step_length, void *data);

/******************* newton function ***********************************/
void
adjoint_output_prestep (ymir_vec_t *solution, const int iter, void *data);

double
subd_adjoint_gradient (ymir_vec_t *edot0, ymir_vec_t *edot1,
                            ymir_vec_t *temp, ymir_vec_t *visc);

#endif
