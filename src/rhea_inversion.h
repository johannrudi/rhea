/** RHEA_INVERSION
 *
 * Inversion for Stokes parameters.
 */

#ifndef RHEA_INVERSION_H
#define RHEA_INVERSION_H

#include <rhea_stokes_problem.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_inversion_add_options (ymir_options_t * opt_sup);

/******************************************************************************
 * Inverse Problem
 *****************************************************************************/

/* inverse problem (opaque) */
typedef struct rhea_inversion_problem rhea_inversion_problem_t;

/**
 * Creates/destroys an inverse problem.
 */
rhea_inversion_problem_t *rhea_inversion_new (
                            rhea_stokes_problem_t *stokes_problem);

void                      rhea_inversion_destroy (
                            rhea_inversion_problem_t *inv_problem);

#endif /* RHEA_INVERSION_H */
