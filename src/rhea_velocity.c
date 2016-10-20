/*
 */

#include <rhea_velocity.h>
#include <rhea_base.h>
#include <ymir_vec_getset.h>
#include <ymir_derivative_elem.h>

void
rhea_velocity_get_elem_gll (sc_dmatrix_t *vel_el_mat,
                            ymir_vec_t *vel_vec,
                            const ymir_locidx_t elid)
{
#ifdef RHEA_ENABLE_DEBUG
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (vel_vec);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  /* check input */
  YMIR_ASSERT_IS_CVEC (vel_vec);
  RHEA_ASSERT (vel_vec->ncfields == 3);
  RHEA_ASSERT (vel_el_mat->m == n_nodes_per_el);
  RHEA_ASSERT (vel_el_mat->n == 3);
#endif

  /* interpolate from continuous nodes to discontinuous GLL nodes */
  ymir_cvec_get_elem_interp (vel_vec, vel_el_mat, YMIR_STRIDE_NODE,
                             elid, YMIR_GLL_NODE, YMIR_READ);
}

void
rhea_velocity_compute_strain_rate_2inv_elem_gauss (
                                        sc_dmatrix_t *strain_rate_2inv_el_mat,
                                        sc_dmatrix_t *vel_el_mat,
                                        ymir_vec_t *vel_vec,
                                        const ymir_locidx_t elid,
                                        sc_dmatrix_t *tmp_grad_vel,
                                        sc_dmatrix_t *tmp_dvel,
                                        sc_dmatrix_t *tmp_vel)
{
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (vel_vec);
  const int           order = ymir_mesh_get_order (mesh);
  const int           n_nodes = ymir_mesh_get_num_nodes_per_elem (mesh);
  const double       *_sc_restrict Dr = mesh->drst->e[0];
  const double       *_sc_restrict Br = mesh->brst->e[0];
  const double       *_sc_restrict rx = mesh->drdx->e[elid];
  const double       *_sc_restrict sx = mesh->dsdx->e[elid];
  const double       *_sc_restrict tx = mesh->dtdx->e[elid];
  const double       *_sc_restrict ry = mesh->drdy->e[elid];
  const double       *_sc_restrict sy = mesh->dsdy->e[elid];
  const double       *_sc_restrict ty = mesh->dtdy->e[elid];
  const double       *_sc_restrict rz = mesh->drdz->e[elid];
  const double       *_sc_restrict sz = mesh->dsdz->e[elid];
  const double       *_sc_restrict tz = mesh->dtdz->e[elid];

  /* check input */
  RHEA_ASSERT (strain_rate_2inv_el_mat->m == n_nodes);
  RHEA_ASSERT (strain_rate_2inv_el_mat->n == 1);
  RHEA_ASSERT (vel_el_mat->m == n_nodes);
  RHEA_ASSERT (vel_el_mat->n == 3);
  RHEA_ASSERT (tmp_grad_vel->m == n_nodes);
  RHEA_ASSERT (tmp_grad_vel->n == 9);
  RHEA_ASSERT (tmp_dvel->m == n_nodes);
  RHEA_ASSERT (tmp_dvel->n == 3);
  RHEA_ASSERT (tmp_vel->m == n_nodes);
  RHEA_ASSERT (tmp_vel->n == 3);

  /* get velocity at GLL nodes */
  rhea_velocity_get_elem_gll (vel_el_mat, vel_vec, elid);

  /* compute the gradient of the velocity at Gauss nodes */
  ymir_derivative_elem_grad (
      order, 3, Dr, Br, rx, sx, tx, ry, sy, ty, rz, sz, tz,
      vel_el_mat, tmp_grad_vel, tmp_dvel, NULL, tmp_vel, 0 /* !transpose */);

  /* compute the 2nd invariant of the strain rate (also at Gauss nodes) */
  {
    const double       *_sc_restrict grad_vel_data = tmp_grad_vel->e[0];
    double             *_sc_restrict sr2_data = strain_rate_2inv_el_mat->e[0];
    int                 nodeid, nodeid_grad = 0;
    double              v, s;

    for (nodeid = 0; nodeid < n_nodes; nodeid++, nodeid_grad += 9) {
      const double       *_sc_restrict g = &grad_vel_data[nodeid_grad];

      /* compute:
       *   sr2 = sqrt ( 1/2 * grad_sym (vel) : grad_sym (vel) )
       * where
       *   grad_sym (vel) = 1/2 * (grad (vel) + grad (vel)^T)
       */
      v  = g[0] * g[0];
      s  = g[1] + g[3]; v += 0.5 * s * s;
      s  = g[2] + g[6]; v += 0.5 * s * s;
      v += g[4] * g[4];
      s  = g[5] + g[7]; v += 0.5 * s * s;
      v += g[8] * g[8];
      sr2_data[nodeid] = sqrt (0.5 * v);
    }
  }
}
