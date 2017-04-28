/*
 */

#include <rhea_strainrate.h>
#include <rhea_base.h>
#include <rhea_velocity.h>
#include <ymir_velocity_elem.h>

ymir_vec_t *
rhea_strainrate_new (ymir_mesh_t *ymir_mesh)
{
  return ymir_dvec_new (ymir_mesh, 6, YMIR_GAUSS_NODE);
}

void
rhea_strainrate_destroy (ymir_vec_t *strainrate)
{
  ymir_vec_destroy (strainrate);
}

int
rhea_strainrate_check_vec_type (ymir_vec_t *vec)
{
  return (
      ymir_vec_is_dvec (vec) &&
      vec->ndfields == 6 &&
      vec->node_type == YMIR_GAUSS_NODE
  );
}

int
rhea_strainrate_is_valid (ymir_vec_t *vec)
{
  return sc_dmatrix_is_valid (vec->dataown);
}

ymir_vec_t *
rhea_strainrate_2inv_new (ymir_mesh_t *ymir_mesh)
{
  return ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
}

void
rhea_strainrate_2inv_destroy (ymir_vec_t *strainrate_2inv)
{
  ymir_vec_destroy (strainrate_2inv);
}

int
rhea_strainrate_2inv_check_vec_type (ymir_vec_t *vec)
{
  return (
      ymir_vec_is_dvec (vec) &&
      vec->ndfields == 1 &&
      vec->node_type == YMIR_GAUSS_NODE
  );
}

int
rhea_strainrate_2inv_is_valid (ymir_vec_t *vec)
{
  return (
      sc_dmatrix_is_valid (vec->dataown) &&
      0.0 <= ymir_dvec_min_global (vec)
  );
}

void
rhea_strainrate_compute_sqrt_of_2inv (ymir_vec_t *strainrate_sqrt_2inv,
                                      ymir_vec_t *velocity)
{
  ymir_mesh_t           *ymir_mesh = ymir_vec_get_mesh (velocity);
  ymir_velocity_elem_t  *vel_elem;

  /* check input */
  RHEA_ASSERT (rhea_strainrate_2inv_check_vec_type (strainrate_sqrt_2inv));
  RHEA_ASSERT (rhea_velocity_check_vec_type (velocity));

  /* create work variables */
  vel_elem = ymir_velocity_elem_new (ymir_mesh->ma->N, ymir_mesh->ma->ompsize);

  /* compute the square root of the second invariant of the strain rate */
  ymir_second_invariant_vec (velocity, strainrate_sqrt_2inv, vel_elem);
  ymir_dvec_sqrt (strainrate_sqrt_2inv, strainrate_sqrt_2inv);

  /* destroy */
  ymir_velocity_elem_destroy (vel_elem);
}

void
rhea_strainrate_compute_sqrt_of_2inv_elem_gauss (
                                        sc_dmatrix_t *strainrate_2inv_el_mat,
                                        sc_dmatrix_t *vel_el_mat,
                                        ymir_vec_t *vel_vec,
                                        const ymir_locidx_t elid,
                                        sc_dmatrix_t *tmp_grad_vel,
                                        sc_dmatrix_t *tmp_dvel,
                                        sc_dmatrix_t *tmp_vel)
{
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (vel_vec);

  /* get velocity at GLL nodes */
  rhea_velocity_get_elem_gll (vel_el_mat, vel_vec, elid);

  /* compute the 2nd invariant of the strain rate at Gauss nodes */
  ymir_velocity_elem_compute_strain_rate_2inv_gauss (
      strainrate_2inv_el_mat, vel_el_mat, mesh, elid,
      tmp_grad_vel, tmp_dvel, tmp_vel);
}
