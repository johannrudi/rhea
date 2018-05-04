/*
 */

#include <rhea_pressure.h>
#include <rhea_base.h>
#include <ymir_pressure_vec.h>
#include <ymir_interp_vec.h>

/******************************************************************************
 * Pressure Vector
 *****************************************************************************/

ymir_vec_t *
rhea_pressure_new (ymir_mesh_t *ymir_mesh, ymir_pressure_elem_t *press_elem)
{
  return ymir_pressure_vec_new (ymir_mesh, press_elem);
}

void
rhea_pressure_destroy (ymir_vec_t *pressure)
{
  ymir_vec_destroy (pressure);
}

static double
rhea_pressure_get_dim_scal (rhea_domain_options_t *domain_options,
                            rhea_temperature_options_t *temp_options,
                            rhea_viscosity_options_t *visc_options)
{
  return visc_options->representative_Pas *
         temp_options->thermal_diffusivity_m2_s /
         (domain_options->radius_max_m * domain_options->radius_max_m);
}

void
rhea_pressure_convert_to_dimensional_Pa (
                                      ymir_vec_t * pressure,
                                      rhea_domain_options_t *domain_options,
                                      rhea_temperature_options_t *temp_options,
                                      rhea_viscosity_options_t *visc_options)
{
  ymir_vec_scale (rhea_pressure_get_dim_scal (domain_options, temp_options,
                                              visc_options),
                  pressure);
}

int
rhea_pressure_check_vec_type (ymir_vec_t *vec, ymir_pressure_elem_t *press_elem)
{
  if (press_elem != NULL) {
    return ymir_pressure_vec_is_press_vec (vec, press_elem);
  }
  else {
    return ymir_pressure_vec_is_press_vec_simple (vec);
  }
}

int
rhea_pressure_is_valid (ymir_vec_t *vec)
{
  return sc_dmatrix_is_valid (vec->dataown);
}

MPI_Offset *
rhea_pressure_segment_offset_create (ymir_vec_t *vec)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vec);
  const mangll_gloidx_t *n_elem_offset = ymir_mesh->ma->mesh->RtoGEO;
  const int           n_fields = vec->nefields;
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  int                 mpisize, mpiret;
  MPI_Offset         *segment_offset;
  int                 r;

  /* get parallel environment */
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize); SC_CHECK_MPI (mpiret);

  /* create segment offsets */
  segment_offset = RHEA_ALLOC (MPI_Offset, mpisize + 1);
  segment_offset[0] = 0;
  for (r = 0; r <= mpisize; r++) {
    segment_offset[r] = (MPI_Offset) (n_fields * n_elem_offset[r]);
  }

  return segment_offset;
}

MPI_Offset
rhea_pressure_segment_offset_get (ymir_vec_t *vec)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vec);
  const mangll_gloidx_t *n_elem_offset = ymir_mesh->ma->mesh->RtoGEO;
  const int           n_fields = vec->nefields;
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  int                 mpirank, mpiret;

  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);
  return (MPI_Offset) (n_fields * n_elem_offset[mpirank]);
}

int
rhea_pressure_segment_size_get (ymir_vec_t *vec)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vec);
  const mangll_gloidx_t *n_elem_offset = ymir_mesh->ma->mesh->RtoGEO;
  const int           n_fields = vec->nefields;
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  int                 mpirank, mpiret;

  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);
  return (int) (n_fields * (n_elem_offset[mpirank+1] - n_elem_offset[mpirank]));
}

double
rhea_pressure_compute_mean (ymir_vec_t *pressure,
                            ymir_pressure_elem_t *press_elem)
{
  /* check input */
  RHEA_ASSERT (rhea_pressure_check_vec_type (pressure, press_elem));
  RHEA_ASSERT (rhea_pressure_is_valid (pressure));

  return ymir_pressure_vec_mean (pressure, press_elem, 0);
}

/******************************************************************************
 * Statistics
 *****************************************************************************/

void
rhea_pressure_stats_get_global (double *abs_min_Pa, double *abs_max_Pa,
                                double *mean_Pa, ymir_vec_t *pressure,
                                ymir_pressure_elem_t *press_elem,
                                rhea_domain_options_t *domain_options,
                                rhea_temperature_options_t *temp_options,
                                rhea_viscosity_options_t *visc_options)
{
  const double        dim_scal = rhea_pressure_get_dim_scal (
                          domain_options, temp_options, visc_options);
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (pressure);
  ymir_vec_t         *press_abs = ymir_dvec_new (ymir_mesh, 1,
                                                 YMIR_GAUSS_NODE);

  /* check input */
  RHEA_ASSERT (rhea_pressure_check_vec_type (pressure, press_elem));
  RHEA_ASSERT (rhea_pressure_is_valid (pressure));

  /* set abs values of pressure */
  ymir_interp_vec (pressure, press_abs);
  ymir_vec_fabs (press_abs, press_abs);

  /* find global values */
  if (abs_min_Pa != NULL) {
    *abs_min_Pa = dim_scal * ymir_vec_min_global (press_abs);
  }
  if (abs_max_Pa != NULL) {
    *abs_max_Pa = dim_scal * ymir_vec_max_global (press_abs);
  }
  if (mean_Pa != NULL) {
    *mean_Pa = dim_scal * rhea_pressure_compute_mean (pressure, press_elem);
  }

  /* destroy */
  ymir_vec_destroy (press_abs);
}
