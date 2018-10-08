#ifndef SUBDUCTION_ADJOINT_H
#define SUBDUCTION_ADJOINT_H

#include <subduction_options.h>
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
subd_adjoint_stencil_visc  (ymir_vec_t * stencil_visc,
                         subd_options_t * subd_options);

void
subd_adjoint_rhs_hessian_forward (ymir_vec_t *rhs_vel_press,
                                    rhea_stokes_problem_t *stokes_problem,
                                    subd_options_t  *subd_opt);

double
subd_adjoint_gradient (ymir_vec_t *edot0, ymir_vec_t *edot1,
                            ymir_vec_t *temp, ymir_vec_t *visc);


#endif
