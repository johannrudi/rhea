/**
 * Shared functions for rhea examples.
 *
 ******************************************************************************
 * Author:             Johann Rudi <johann@ices.utexas.edu>
 *****************************************************************************/

#include <example_share_mesh.h>
#include <rhea_base.h>

void
example_share_mesh_new (p4est_t **p4est,
                        ymir_mesh_t **ymir_mesh,
                        ymir_pressure_elem_t **press_elem,
                        sc_MPI_Comm mpicomm,
                        rhea_domain_options_t *domain_options,
                        rhea_discretization_options_t *discr_options)
{
  const char         *this_fn_name = "example_share_mesh_new";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* create p4est */
  *p4est = rhea_discretization_p4est_new (mpicomm, discr_options,
                                          domain_options);

  /* set up boundary, store in `discr_options` */
  rhea_discretization_boundary_create (discr_options, *p4est, domain_options);

  /* create ymir mesh and pressure element */
  rhea_discretization_ymir_mesh_new_from_p4est (ymir_mesh, press_elem, *p4est,
                                                discr_options);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

void
example_share_mesh_destroy (ymir_mesh_t *ymir_mesh,
                            ymir_pressure_elem_t *press_elem,
                            p4est_t *p4est,
                            rhea_discretization_options_t *discr_options)
{
  const char         *this_fn_name = "example_share_mesh_destroy";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* destroy mesh */
  rhea_discretization_ymir_mesh_destroy (ymir_mesh, press_elem);
  rhea_discretization_p4est_destroy (p4est);

  /* destroy boundary */
  rhea_discretization_boundary_clear (discr_options);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}
