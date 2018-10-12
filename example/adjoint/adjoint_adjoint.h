#ifndef ADJOINT_ADJOINT_H
#define ADJOINT_ADJOINT_H

#include <adjoint_options.h>
#include <adjoint_io.h>
//#include <adjoint_essential.h>
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

void
adjoint_data_initialize (ymir_vec_t *solution,
                         void *data);

void adjoint_compute_gradient_norm (ymir_vec_t *neg_gradient,
                                    void *data, double *norm);
void
adjoint_compute_negative_gradient (ymir_vec_t *neg_gradient,
                                   ymir_vec_t *solution, void *data);

double
adjoint_evaluate_objective (ymir_vec_t *solution,
                          void *data);

typedef struct adjoint_problem
{
 ymir_vec_t *solution;

 ymir_vec_t           *sol_vel_press;
 ymir_vec_t           *usol;
 ymir_vec_t           *vsol;
 ymir_vec_t           *Husol;
 ymir_vec_t           *Hvsol;

 // initial viscosity, use parameter or ymir_vec_t visc?

 rhea_stokes_problem_t  *stokes_problem;
 p4est_t                *p4est;
 ymir_mesh_t            *ymir_mesh;
 ymir_pressure_elem_t   *press_elem;

 rhea_discretization_options_t  *discr_options;
 rhea_temperature_options_t     *temp_options;
 subd_options_t                 *subd_options;
 const char                     *vtk_write_input_path;

 int                        solver_iter_max;
 double                     solver_rel_tol;
}
adjoint_problem_t;

void
adjoint_setup_adjoint_problem (adjoint_problem_t *adjoint_problem,
                              rhea_stokes_problem_t *stokes_problem,
                              p4est_t *p4est, ymir_mesh_t *ymir_mesh,
                              ymir_pressure_elem_t *press_elem,
                              rhea_discretization_options_t *discr_options,
                              rhea_temperature_options_t *temp_options,
                              subd_options_t *subd_options,
                              const char *vtk_write_input_path,
                              int     solver_iter_max,
                              double  solver_rel_tol);

void
adjoint_setup_newton (rhea_newton_problem_t **newton_problem,
                      adjoint_problem_t *adjoint_problem);

void
subd_adjoint_stencil_visc  (ymir_vec_t * stencil_visc,
                         subd_options_t * subd_options);

void
subd_adjoint_rhs_hessian_forward (ymir_vec_t *rhs_vel_press,
                                    rhea_stokes_problem_t *stokes_problem,
                                    subd_options_t  *subd_opt);

double *
adjoint_vec_gradient (ymir_vec_t *usol, ymir_vec_t *vsol);

double
subd_adjoint_gradient (ymir_vec_t *edot0, ymir_vec_t *edot1,
                            ymir_vec_t *temp, ymir_vec_t *visc);


#endif
