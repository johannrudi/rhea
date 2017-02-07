/*
 */

#include <rhea_stress.h>
#include <rhea_base.h>
#include <rhea_strainrate.h>
#include <rhea_viscosity.h>

ymir_vec_t *
rhea_stress_2inv_new (ymir_mesh_t *ymir_mesh)
{
  return ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
}

void
rhea_stress_2inv_destroy (ymir_vec_t *stress_2inv)
{
  ymir_vec_destroy (stress_2inv);
}

int
rhea_stress_2inv_check_vec_type (ymir_vec_t *vec)
{
  return (
      ymir_vec_is_dvec (vec) &&
      vec->ndfields == 1 &&
      vec->node_type == YMIR_GAUSS_NODE
  );
}

int
rhea_stress_2inv_is_valid (ymir_vec_t *vec)
{
  return (
      sc_dmatrix_is_valid (vec->dataown) &&
      0.0 <= ymir_dvec_min_global (vec)
  );
}

void
rhea_stress_compute_viscstress_sqrt_of_2inv (ymir_vec_t *viscstress_sqrt_2inv,
                                             ymir_vec_t *strainrate_sqrt_2inv,
                                             ymir_vec_t *viscosity)
{
  RHEA_ASSERT (rhea_stress_2inv_check_vec_type (viscstress_sqrt_2inv));
  RHEA_ASSERT (rhea_strainrate_2inv_check_vec_type (viscstress_sqrt_2inv));
  RHEA_ASSERT (rhea_viscosity_check_vec_type (viscstress_sqrt_2inv));

  ymir_dvec_copy (strainrate_sqrt_2inv, viscstress_sqrt_2inv);
  ymir_dvec_multiply_in (viscosity, viscstress_sqrt_2inv);
  ymir_dvec_scale (2.0, viscstress_sqrt_2inv);
}
