/* EXAMPLE_SHARE_MESH  Shared functions for rhea examples.
 *
 * Author:             Johann Rudi <johann@ices.utexas.edu>
 */

#ifndef EXAMPLE_SHARE_MESH_H
#define EXAMPLE_SHARE_MESH_H

#include <rhea_domain.h>
#include <rhea_topography.h>
#include <rhea_discretization.h>
#include <ymir_pressure_elem.h>

/**
 * Creates a new mesh.
 */
void                example_share_mesh_new (
                                p4est_t **p4est,
                                ymir_mesh_t **ymir_mesh,
                                ymir_pressure_elem_t **press_elem,
                                sc_MPI_Comm mpicomm,
                                rhea_domain_options_t *domain_options,
                                rhea_topography_options_t *topo_options,
                                rhea_discretization_options_t *discr_options,
                                const int performance_monitor_index_mesh);

/**
 * Destroys a mesh.
 */
void                example_share_mesh_destroy (
                                ymir_mesh_t *ymir_mesh,
                                ymir_pressure_elem_t *press_elem,
                                p4est_t *p4est,
                                rhea_topography_options_t *topo_options,
                                rhea_discretization_options_t *discr_options);

#endif /* EXAMPLE_SHARE_MESH_H */
