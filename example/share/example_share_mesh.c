#include <example_share_mesh.h>
#include <rhea.h>

void
example_share_mesh_new (p4est_t **p4est,
                        ymir_mesh_t **ymir_mesh,
                        ymir_pressure_elem_t **press_elem,
                        sc_MPI_Comm mpicomm,
                        rhea_domain_options_t *domain_options,
                        rhea_topography_options_t *topo_options,
                        rhea_discretization_options_t *discr_options,
                        const int performance_monitor_index_mesh)
{
  rhea_performance_monitor_start_barrier (performance_monitor_index_mesh);
  RHEA_GLOBAL_PRODUCTION_FN_BEGIN (__func__);

  /* set up data */
  rhea_topography_data_create (topo_options, mpicomm);

  /* create p4est */
  *p4est = rhea_discretization_p4est_new (mpicomm, discr_options,
                                          domain_options);

  /* set up boundary, store in `discr_options` */
  rhea_discretization_boundary_create (discr_options, *p4est, domain_options);

  /* create ymir mesh and pressure element */
  rhea_discretization_ymir_mesh_new_from_p4est (ymir_mesh, press_elem, *p4est,
                                                discr_options);

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);
  rhea_performance_monitor_stop_add (performance_monitor_index_mesh);
}

void
example_share_mesh_destroy (ymir_mesh_t *ymir_mesh,
                            ymir_pressure_elem_t *press_elem,
                            p4est_t *p4est,
                            rhea_topography_options_t *topo_options,
                            rhea_discretization_options_t *discr_options)
{
  RHEA_GLOBAL_PRODUCTION_FN_BEGIN (__func__);

  /* destroy mesh */
  rhea_discretization_ymir_mesh_destroy (ymir_mesh, press_elem);
  rhea_discretization_p4est_destroy (p4est);

  /* destroy boundary */
  rhea_discretization_boundary_clear (discr_options);

  /* destroy data */
  rhea_topography_data_clear (topo_options);

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);
}
