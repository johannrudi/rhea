#ifndef RHEA_STOKES_PROBLEM_AMR_H
#define RHEA_STOKES_PROBLEM_AMR_H

#include <rhea_stokes_problem.h>
#include <rhea_discretization.h>

int                 rhea_stokes_problem_amr (
                                  rhea_stokes_problem_t *stokes_problem,
                                  ymir_vec_t **velocity_pressure,
                                  p4est_t *p4est,
                                  rhea_discretization_options_t *discr_options);

#endif /* RHEA_STOKES_PROBLEM_AMR_H */
