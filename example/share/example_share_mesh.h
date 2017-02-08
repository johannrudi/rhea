/**
 * Shared functions for rhea examples.
 *
 ******************************************************************************
 * Author:             Johann Rudi <johann@ices.utexas.edu>
 *****************************************************************************/

#ifndef EXAMPLE_SHARE_MESH_H
#define EXAMPLE_SHARE_MESH_H

#include <rhea_domain.h>
#include <rhea_discretization.h>
#include <ymir_pressure_elem.h>

/**
 * Creates a new mesh.
 */
void                example_share_mesh_new (
                                p4est_t **p4est,
                                ymir_mesh_t **ymir_mesh,
                                ymir_pressure_elem_t **press_elem,
                                MPI_Comm mpicomm,
                                rhea_domain_options_t *domain_options,
                                rhea_discretization_options_t *discr_options);

/**
 * Destroys a mesh.
 */
void                example_share_mesh_destroy (
                                ymir_mesh_t *ymir_mesh,
                                ymir_pressure_elem_t *press_elem,
                                p4est_t *p4est,
                                rhea_discretization_options_t *discr_options);

#endif /* EXAMPLE_SHARE_MESH_H */
