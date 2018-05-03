/*
 */

#include <rhea_vis.h>
#include <rhea_base.h>
#include <ymir_elem_iter.h>
#include <ymir_vec_optimized.h>
#ifdef RHEA_USE_CATALYST
# include <rhea_vis_adaptor.h>
#endif

void
rhea_vis_initialize (const char *catalyst_scripts[],
                     const int n_catalyst_scripts)
{
#ifdef RHEA_USE_CATALYST
  RHEA_GLOBAL_VERBOSEF ("Into %s\n", __func__);
  rhea_vis_adaptor_initialize (catalyst_scripts, n_catalyst_scripts);
  RHEA_GLOBAL_VERBOSEF ("Done %s\n", __func__);
#endif
}

void
rhea_vis_finalize ()
{
#ifdef RHEA_USE_CATALYST
  RHEA_GLOBAL_VERBOSEF ("Into %s\n", __func__);
  rhea_vis_adaptor_finalize ();
  RHEA_GLOBAL_VERBOSEF ("Done %s\n", __func__);
#endif
}

void
rhea_vis_process_primary (ymir_vec_t *velocity, ymir_vec_t *pressure)
{
#ifdef RHEA_USE_CATALYST
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (velocity);
  const ymir_locidx_t n_elements = ymir_mesh_get_num_elems_loc (ymir_mesh);
  const int           order = ymir_mesh_get_order (ymir_mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (
                                                                  ymir_mesh);
  unsigned int        n_nodes = n_elements * n_nodes_per_el;
  const int           node_map[27] = { 0,  2,  8,  6, // bottom vertices
                                      18, 20, 26, 24, // top vertices
                                       1,  5,  7,  3, // bottom edge centers
                                      19, 23, 25, 21, // top edge centers
                                       9, 11, 17, 15, // mid edge centers
                                      12, 14, 10, 16, // mid faces
                                       4, 22,         // bottom & top faces
                                      13};            // volume center
  unsigned int       *elem_to_coord_idx;
  double             *coordinates;
  double             *velocity_data;
  double             *pressure_data;
  ymir_locidx_t       elid;
  int                 i;

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", __func__);

  /* check input */
  RHEA_ASSERT (order == 2);

  /* create work variables */
  elem_to_coord_idx = RHEA_ALLOC (unsigned int, n_nodes);
  coordinates = RHEA_ALLOC (double, n_nodes * 3);
  velocity_data = RHEA_ALLOC (double, n_nodes * 3);
  pressure_data = RHEA_ALLOC (double, n_nodes);

  /* set coordinates */
  {
    for (elid = 0; elid < n_elements; elid++) {
      const double *x = ymir_mesh_get_elem_coord_x (elid, ymir_mesh);
      const double *y = ymir_mesh_get_elem_coord_y (elid, ymir_mesh);
      const double *z = ymir_mesh_get_elem_coord_z (elid, ymir_mesh);

      for (i = 0; i < n_nodes_per_el; i++) {
        const int           nodeid = node_map[i];
        const ymir_locidx_t idx = n_nodes_per_el*elid + i;

        elem_to_coord_idx[idx] = idx;
        coordinates[3*idx    ] = x[nodeid];
        coordinates[3*idx + 1] = y[nodeid];
        coordinates[3*idx + 2] = z[nodeid];
      }
    }
  }

  /* set velocity and pressure */
  {
    const double       *_sc_restrict Binv = ymir_mesh->ma->refel->Brinv->e[0];
    const double       *_sc_restrict Binvt = ymir_mesh->ma->refel->Brinvt->e[0];
    sc_dmatrix_t       *elmat_vel = sc_dmatrix_new (n_nodes_per_el, 3);
    sc_dmatrix_t       *elmat_press_gauss = sc_dmatrix_new (1, n_nodes_per_el);
    sc_dmatrix_t       *elmat_press_gll = sc_dmatrix_new (1, n_nodes_per_el);
    double             *_sc_restrict U = elmat_vel->e[0];
    double             *_sc_restrict P_gauss = elmat_press_gauss->e[0];
    double             *_sc_restrict P_gll = elmat_press_gll->e[0];
    ymir_elem_iter_t    iter;

    ymir_elem_iter_init_dirty (&iter, 1 /* ompsize */, 0 /* ompthread */,
                               velocity, &velocity);
    for (iter.elid = 0; iter.elid < n_elements; iter.elid++) {
      for (i = 0; i < n_nodes_per_el; i++) {
        const int           nodeid = node_map[i];
        const ymir_locidx_t idx = n_nodes_per_el*iter.elid + i;

        ymir_vel_get_elem_optimized (elmat_vel, velocity, &iter);
        ymir_press_get_elem_optimized (elmat_press_gauss, pressure, iter.elid,
                                       0 /* ompthread */);
        mangll_tensor_optimized_interp_scal (P_gll, Binv, Binvt, P_gauss,
                                             order, 0 /* !transpose */);

        velocity_data[3*idx    ] = U[3*nodeid    ];
        velocity_data[3*idx + 1] = U[3*nodeid + 1];
        velocity_data[3*idx + 2] = U[3*nodeid + 2];
        pressure_data[idx] = P_gll[nodeid];
      }
    }

    sc_dmatrix_destroy (elmat_vel);
    sc_dmatrix_destroy (elmat_press_gauss);
    sc_dmatrix_destroy (elmat_press_gll);
  }

  /* process visualization */
  rhea_vis_adaptor_process (order, n_nodes, coordinates,
                            n_elements, elem_to_coord_idx,
                            velocity_data, pressure_data);

  /* destroy */
  RHEA_FREE (elem_to_coord_idx);
  RHEA_FREE (coordinates);
  RHEA_FREE (velocity_data);
  RHEA_FREE (pressure_data);

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", __func__);
#endif
}
