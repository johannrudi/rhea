/*
 */

#include <rhea_velocity.h>
#include <rhea_base.h>
#include <ymir_vec_getset.h>
#include <ymir_velocity_elem.h>

ymir_vec_t *
rhea_velocity_new (ymir_mesh_t *ymir_mesh)
{
  return ymir_cvec_new (ymir_mesh, 3);
}

void
rhea_velocity_destroy (ymir_vec_t *velocity)
{
  ymir_vec_destroy (velocity);
}

int
rhea_velocity_check_vec_type (ymir_vec_t *vec)
{
  return (
      ymir_vec_is_cvec (vec) &&
      vec->ncfields == 3 &&
      vec->node_type == YMIR_GLL_NODE
  );
}

int
rhea_velocity_is_valid (ymir_vec_t *vec)
{
  return sc_dmatrix_is_valid (vec->dataown) && sc_dmatrix_is_valid (vec->coff);
}

/******************************************************************************
 * Get & Set Functions
 *****************************************************************************/

void
rhea_velocity_get_elem_gll (sc_dmatrix_t *vel_el_mat,
                            ymir_vec_t *vel_vec,
                            const ymir_locidx_t elid)
{
#ifdef RHEA_ENABLE_DEBUG
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (vel_vec);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (vel_vec));
  RHEA_ASSERT (vel_el_mat->m == n_nodes_per_el);
  RHEA_ASSERT (vel_el_mat->n == 3);
#endif

  /* interpolate from continuous nodes to discontinuous GLL nodes */
  ymir_cvec_get_elem_interp (vel_vec, vel_el_mat, YMIR_STRIDE_NODE,
                             elid, YMIR_GLL_NODE, YMIR_READ);
}
