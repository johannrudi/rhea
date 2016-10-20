/*
 */

#include <rhea_temperature.h>
#include <rhea_base.h>
#include <ymir_vec_getset.h>

/**
 * Restrict temperature to its valid range [0,1].
 */
static void
rhea_temperature_bound_range_data (double *_sc_restrict temp_data,
                                   const sc_bint_t size)
{
  sc_bint_t           i;

  for (i = 0; i < size; i++) {
    RHEA_ASSERT (isfinite (temp_data[i]));
    temp_data[i] = SC_MAX (0.0, SC_MIN (temp_data[i], 1.0));
  }
}

void
rhea_temperature_get_elem_gauss (sc_dmatrix_t *temp_el_mat,
                                 ymir_vec_t *temp_vec,
                                 const ymir_locidx_t elid)
{
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (temp_vec);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  /* check input */
  YMIR_ASSERT_IS_CVEC (temp_vec);
  RHEA_ASSERT (temp_vec->ncfields == 1);
  RHEA_ASSERT (temp_el_mat->m == n_nodes_per_el);
  RHEA_ASSERT (temp_el_mat->n == 1);

  /* interpolate from continuous nodes to (discontinuous) Gauss nodes */
  ymir_cvec_get_elem_interp (temp_vec, temp_el_mat, YMIR_STRIDE_NODE, elid,
                             YMIR_GAUSS_NODE, YMIR_READ);

  /* restrict temperature to valid range to account for interpolation erros */
  rhea_temperature_bound_range_data (temp_el_mat->e[0], n_nodes_per_el);
}
