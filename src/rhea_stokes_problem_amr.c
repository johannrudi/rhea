/*
 */

#include <rhea_stokes_problem_amr.h>
#include <rhea_base.h>
#include <rhea_amr.h>
#include <rhea_discretization.h>
#include <rhea_temperature.h>
#include <rhea_velocity.h>
#include <rhea_pressure.h>
#include <rhea_velocity_pressure.h>
#include <rhea_viscosity.h>
#include <ymir_interp_vec.h>
#include <ymir_hmg_intergrid_h.h>
#include <ymir_mass_vec.h>
#include <ymir_pressure_vec.h>
#include <mangll_fields.h>

/******************************************************************************
 * AMR Data
 *****************************************************************************/

typedef struct rhea_stokes_problem_amr_data
{
  /* options (not owned) */
  rhea_discretization_options_t  *discr_options;
  rhea_viscosity_options_t       *visc_options;

  /* mesh (owned) */
  ymir_mesh_t        *ymir_mesh;
  ymir_pressure_elem_t *press_elem;

  /* fields of a Stokes problem (owned) */
  ymir_vec_t         *temperature;
  ymir_vec_t         *velocity;
  ymir_vec_t         *pressure;

  /* AMR iterations */
  int                 first_amr;
  mangll_t           *mangll_original;
  mangll_t           *mangll_adapted;
  mangll_t           *mangll_partitioned;

  /* buffers for projection & partitioning */
  sc_dmatrix_t       *temperature_original;
  sc_dmatrix_t       *temperature_adapted;

  sc_dmatrix_t       *velocity_original;
  sc_dmatrix_t       *velocity_adapted;

  sc_dmatrix_t       *pressure_original;
  sc_dmatrix_t       *pressure_adapted;
}
rhea_stokes_problem_amr_data_t;

static void
rhea_stokes_problem_amr_data_fields_clear (
                                      rhea_stokes_problem_amr_data_t *amr_data)
{
  if (amr_data->temperature != NULL) {
    ymir_vec_destroy (amr_data->temperature);
  }
  if (amr_data->velocity != NULL) {
    ymir_vec_destroy (amr_data->velocity);
  }
  if (amr_data->pressure != NULL) {
    ymir_vec_destroy (amr_data->pressure);
  }

  amr_data->temperature = NULL;
  amr_data->velocity = NULL;
  amr_data->pressure = NULL;
}

static void
rhea_stokes_problem_amr_data_buffer_clear (
                                      rhea_stokes_problem_amr_data_t *amr_data)
{
  if (amr_data->temperature_original != NULL) {
    sc_dmatrix_destroy (amr_data->temperature_original);
  }
  if (amr_data->temperature_adapted != NULL) {
    sc_dmatrix_destroy (amr_data->temperature_adapted);
  }

  if (amr_data->velocity_original != NULL) {
    sc_dmatrix_destroy (amr_data->velocity_original);
  }
  if (amr_data->velocity_adapted != NULL) {
    sc_dmatrix_destroy (amr_data->velocity_adapted);
  }

  if (amr_data->pressure_original != NULL) {
    sc_dmatrix_destroy (amr_data->pressure_original);
  }
  if (amr_data->pressure_adapted != NULL) {
    sc_dmatrix_destroy (amr_data->pressure_adapted);
  }

  amr_data->temperature_original = NULL;
  amr_data->temperature_adapted = NULL;
  amr_data->velocity_original = NULL;
  amr_data->velocity_adapted = NULL;
  amr_data->pressure_original = NULL;
  amr_data->pressure_adapted = NULL;
}

rhea_stokes_problem_amr_data_t *
rhea_stokes_problem_amr_data_new (ymir_mesh_t *ymir_mesh,
                                  ymir_pressure_elem_t *press_elem,
                                  ymir_vec_t *temperature,
                                  ymir_vec_t *velocity_pressure,
                                  rhea_discretization_options_t *discr_options,
                                  rhea_viscosity_options_t *visc_options)
{
  rhea_stokes_problem_amr_data_t *amr_data;

  /* check input */
  RHEA_ASSERT (rhea_temperature_check_vec_type (temperature));
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (velocity_pressure));

  /* create */
  amr_data = RHEA_ALLOC (rhea_stokes_problem_amr_data_t, 1);

  /* set data */
  amr_data->ymir_mesh = ymir_mesh;
  amr_data->press_elem = press_elem;

  amr_data->discr_options = discr_options;
  amr_data->visc_options = visc_options;

  if (temperature != NULL) {
    amr_data->temperature = rhea_temperature_new (ymir_mesh);
    ymir_vec_copy (temperature, amr_data->temperature);
  }
  else {
    amr_data->temperature = NULL;
  }
  if (velocity_pressure != NULL) {
    amr_data->velocity = rhea_velocity_new (ymir_mesh);
    amr_data->pressure = rhea_pressure_new (ymir_mesh, press_elem);
    rhea_velocity_pressure_copy_components (
        amr_data->velocity, amr_data->pressure, velocity_pressure, press_elem);
  }
  else {
    amr_data->velocity = NULL;
    amr_data->pressure = NULL;
  }

  amr_data->first_amr = -1;
  amr_data->mangll_original = NULL;
  amr_data->mangll_adapted = NULL;
  amr_data->mangll_partitioned = NULL;

  /* clear buffers */
  rhea_stokes_problem_amr_data_buffer_clear (amr_data);

  return amr_data;
}

void
rhea_stokes_problem_amr_data_destroy (rhea_stokes_problem_amr_data_t *amr_data)
{
  rhea_stokes_problem_amr_data_fields_clear (amr_data);
  rhea_stokes_problem_amr_data_buffer_clear (amr_data);
  RHEA_FREE (amr_data);
}

/******************************************************************************
 * AMR for Single Field
 *****************************************************************************/

static sc_dmatrix_t *
rhea_stokes_problem_amr_init_field (ymir_vec_t * vector, const int n_fields)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vector);
  const mangll_locidx_t n_elements = ymir_mesh->ma->mesh->K;
  const int           n_nodes_per_el = ymir_mesh->ma->Np;
  sc_dmatrix_t       *buffer_mat;
  ymir_vec_t         *buffer_view;

  /* check input */
  RHEA_ASSERT (!ymir_vec_is_cvec (vector) || n_fields == vector->ncfields);
  RHEA_ASSERT (!ymir_vec_is_dvec (vector) || n_fields == vector->ndfields);
  RHEA_ASSERT (!ymir_vec_is_evec (vector) || n_fields == 1);

  /* create buffer and view onto that buffer */
  buffer_mat = sc_dmatrix_new (n_elements, n_nodes_per_el);
  buffer_view = ymir_dvec_new_data (ymir_mesh, n_fields, YMIR_GAUSS_NODE,
                                    buffer_mat);

  /* project onto Gauss nodes */
  ymir_interp_vec (vector, buffer_view);

  /* destroy */
  ymir_vec_destroy (buffer_view);

  /* return buffer */
  return buffer_mat;
}

static sc_dmatrix_t *
rhea_stokes_problem_amr_project_field (sc_dmatrix_t * buffer_original,
                                       const int n_fields,
                                       mangll_t * mangll_original,
                                       mangll_t * mangll_adapted)
{
  /* adapted mesh */
  const mangll_locidx_t  n_elements = mangll_adapted->mesh->K;
  const int           n_nodes_per_el = mangll_adapted->Np;
  sc_dmatrix_t       *buffer_adapted;

  /* create buffer for adapted data */
  buffer_adapted = sc_dmatrix_new (n_elements, n_nodes_per_el);

  /* project original data */
  ymir_hmg_intergrid_h_project_gauss (
      buffer_adapted, mangll_adapted,
      buffer_original, mangll_original, n_fields);

  /* destroy buffer for original data */
  sc_dmatrix_destroy (buffer_original);

  /* return buffer for adapted data */
  return buffer_adapted;
}

static sc_dmatrix_t *
rhea_stokes_problem_amr_partition_field (sc_dmatrix_t * buffer_adapted,
                                         const int n_fields,
                                         mangll_t * mangll_adapted,
                                         mangll_t * mangll_partitioned)
{
  /* adapted mesh */
  mangll_gloidx_t    *RtoGEO_adapted = mangll_adapted->mesh->RtoGEO;
  MPI_Comm            mpicomm = mangll_adapted->mpicomm;
  const int           mpisize = mangll_adapted->mpisize;
  const int           mpirank = mangll_adapted->mpirank;
  /* partitoned mesh */
  mangll_gloidx_t    *RtoGEO_partitioned = mangll_partitioned->mesh->RtoGEO;
  const mangll_locidx_t  n_elements = mangll_partitioned->mesh->K;
  const int           n_nodes_per_el = mangll_partitioned->Np;
  sc_dmatrix_t       *buffer_partitioned;

  /* create buffer for partitioned data */
  buffer_partitioned = sc_dmatrix_new (n_elements, n_nodes_per_el);

  /* partition adapted data */
  mangll_field_partition (n_fields, n_nodes_per_el, mpirank, mpisize, mpicomm,
                          RtoGEO_adapted, buffer_adapted,
                          RtoGEO_partitioned, buffer_partitioned);

  /* destroy buffer for adapted data */
  sc_dmatrix_destroy (buffer_adapted);

  /* return buffer for partitioned data */
  return buffer_partitioned;
}

/******************************************************************************
 * AMR Callback Functions
 *****************************************************************************/

static void
rhea_stokes_problem_amr_data_initialize_fn (p4est_t *p4est, void *data)
{
  rhea_stokes_problem_amr_data_t *d = data;
  const int           has_temp = (d->temperature != NULL);
  const int           has_vel = (d->velocity != NULL);
  const int           has_press = (d->pressure != NULL);

  /* check input */
  RHEA_ASSERT (d->mangll_original == NULL);
  RHEA_ASSERT (d->mangll_adapted == NULL);
  RHEA_ASSERT (d->mangll_partitioned == NULL);

  /*
   * Fields
   */

  /* project fields onto Gauss nodes, store in buffers `*_original` */
  if (has_temp) {
    RHEA_ASSERT (rhea_temperature_check_vec_type (d->temperature));
    RHEA_ASSERT (d->temperature_original == NULL);
    d->temperature_original = rhea_stokes_problem_amr_init_field (
        d->temperature, 1 /* #fields */);
  }
  if (has_vel) {
    RHEA_ASSERT (rhea_velocity_check_vec_type (d->velocity));
    RHEA_ASSERT (d->velocity_original == NULL);
    d->velocity_original = rhea_stokes_problem_amr_init_field (
        d->velocity, 3 /* #fields */);
  }
  if (has_press) {
    RHEA_ASSERT (rhea_pressure_check_vec_type (d->pressure, d->press_elem));
    RHEA_ASSERT (d->pressure_original == NULL);
    d->pressure_original = rhea_stokes_problem_amr_init_field (
        d->pressure, 1 /* #fields */);
  }

  /* destroy fields */
  rhea_stokes_problem_amr_data_fields_clear (d);

  /*
   * Mesh
   */

  /* set AMR parameters */
  d->first_amr = 1;
  d->mangll_original = d->ymir_mesh->ma;
  d->mangll_adapted = NULL;
  d->mangll_partitioned = NULL;

  /* destroy mesh partially (keep only the mangll object) */
  rhea_discretization_mangll_continuous_destroy (
      NULL, d->ymir_mesh->cnodes);
  ymir_mesh_destroy (d->ymir_mesh);
  ymir_pressure_elem_destroy (d->press_elem);
  d->ymir_mesh = NULL;
  d->press_elem = NULL;
}

static void
rhea_stokes_problem_amr_data_finalize_fn (p4est_t *p4est, void *data)
{
  rhea_stokes_problem_amr_data_t *d = data;
  const int           has_temp = (d->temperature_original != NULL);
  const int           has_vel = (d->velocity_original != NULL);
  const int           has_press = (d->pressure_original != NULL);
  ymir_vec_t         *temperature_view = NULL;
  ymir_vec_t         *velocity_view = NULL;
  ymir_vec_t         *pressure_view = NULL;
  ymir_vec_t         *temperature_mass;
  ymir_vec_t         *velocity_mass;
  ymir_vec_t         *pressure_lump;

  /* check input */
  RHEA_ASSERT (d->mangll_original != NULL);
  RHEA_ASSERT (d->mangll_adapted == NULL);
  RHEA_ASSERT (d->mangll_partitioned == NULL);
  RHEA_ASSERT (d->ymir_mesh == NULL);
  RHEA_ASSERT (d->press_elem == NULL);

  /*
   * Mesh
   */

  /* destroy discontinuous mangll */
  RHEA_ASSERT (d->first_amr == 0);
  rhea_discretization_mangll_discontinuous_destroy (d->mangll_original);
  d->mangll_original = NULL;

  /* create ymir mesh */
  rhea_discretization_ymir_mesh_new_from_p4est (
      &(d->ymir_mesh), &(d->press_elem), p4est, d->discr_options);

  /*
   * Fields
   */

  /* create views onto buffers and fields */
  if (has_temp) {
    RHEA_ASSERT (d->temperature == NULL);
    temperature_view = ymir_dvec_new_data (d->ymir_mesh, 1, YMIR_GAUSS_NODE,
                                           d->temperature_original);
    temperature_mass = rhea_temperature_new (d->ymir_mesh);
    d->temperature = rhea_temperature_new (d->ymir_mesh);
  }
  if (has_vel) {
    RHEA_ASSERT (d->velocity == NULL);
    velocity_view = ymir_dvec_new_data (d->ymir_mesh, 3, YMIR_GAUSS_NODE,
                                        d->velocity_original);
    velocity_mass = rhea_velocity_new (d->ymir_mesh);
    d->velocity = rhea_velocity_new (d->ymir_mesh);
  }
  if (has_press) {
    RHEA_ASSERT (d->pressure == NULL);
    pressure_view = ymir_dvec_new_data (d->ymir_mesh, 1, YMIR_GAUSS_NODE,
                                        d->pressure_original);
    pressure_lump = rhea_pressure_new (d->ymir_mesh, d->press_elem);
    d->pressure = rhea_pressure_new (d->ymir_mesh, d->press_elem);
  }

  /* project back fields from Gauss nodes, stored in buffers `*_original` */
  if (has_temp) {
    ymir_mass_apply_gauss (temperature_view);
    ymir_interp_vec (temperature_view, temperature_mass);
    ymir_mass_invert (temperature_mass, d->temperature);
    ymir_vec_destroy (temperature_mass);
    ymir_vec_destroy (temperature_view);
  }
  if (has_vel) {
    ymir_mass_apply_gauss (velocity_view);
    ymir_interp_vec (velocity_view, velocity_mass);
    ymir_mass_invert (velocity_mass, d->velocity);
    ymir_vec_destroy (velocity_mass);
    ymir_vec_destroy (velocity_view);
  }
  if (has_press) {
    ymir_pressure_vec_lump_mass (pressure_lump, d->press_elem);

    ymir_mass_apply_gauss (pressure_view);
    ymir_interp_vec (pressure_view, d->pressure);
    ymir_vec_divide_in (pressure_lump, d->pressure);
    ymir_vec_destroy (pressure_view);
    ymir_vec_destroy (pressure_lump);
  }

  /* destroy buffers */
  rhea_stokes_problem_amr_data_buffer_clear (d);
}

static void
rhea_stokes_problem_amr_data_project_fn (p4est_t *p4est, void *data)
{
  rhea_stokes_problem_amr_data_t *d = data;

  /* check input */
  RHEA_ASSERT (d->mangll_original != NULL);
  RHEA_ASSERT (d->mangll_adapted == NULL);
  RHEA_ASSERT (d->mangll_partitioned == NULL);
  RHEA_ASSERT (d->first_amr == 0 || d->first_amr == 1);

  /* create adapted mangll object */
  d->mangll_adapted = rhea_discretization_mangll_interpolation_new (
      p4est, d->discr_options->order);
  RHEA_ASSERT (d->mangll_adapted != NULL);
  RHEA_ASSERT (d->mangll_adapted->N == d->mangll_original->N);

  /* perform projections */
  {
    const int           has_temp = (d->temperature_original != NULL);
    const int           has_vel = (d->velocity_original != NULL);
    const int           has_press = (d->pressure_original != NULL);

    /* project from original to adapted */
    if (has_temp) {
      RHEA_ASSERT (d->temperature_adapted == NULL);
      d->temperature_adapted = rhea_stokes_problem_amr_project_field (
          d->temperature_original, 1 /* #fields */,
          d->mangll_original, d->mangll_adapted);
      d->temperature_original = NULL;

      /* bound values */
      //TODO
    }
    if (has_vel) {
      RHEA_ASSERT (d->velocity_adapted == NULL);
      d->velocity_adapted = rhea_stokes_problem_amr_project_field (
          d->velocity_original, 3 /* #fields */,
          d->mangll_original, d->mangll_adapted);
      d->velocity_original = NULL;
    }
    if (has_press) {
      RHEA_ASSERT (d->pressure_adapted == NULL);
      d->pressure_adapted = rhea_stokes_problem_amr_project_field (
          d->pressure_original, 1 /* #fields */,
          d->mangll_original, d->mangll_adapted);
      d->pressure_original = NULL;
    }
  }

  /* destroy original mangll */
  if (d->first_amr) {
    rhea_discretization_mangll_continuous_destroy (d->mangll_original, NULL);
    d->first_amr = 0;
  }
  else {
    rhea_discretization_mangll_discontinuous_destroy (d->mangll_original);
  }
  d->mangll_original = NULL;
}

static void
rhea_stokes_problem_amr_data_partition_fn (p4est_t *p4est, void *data)
{
  rhea_stokes_problem_amr_data_t *d = data;

  /* check input */
  RHEA_ASSERT (d->mangll_original == NULL);
  RHEA_ASSERT (d->mangll_adapted != NULL);
  RHEA_ASSERT (d->mangll_partitioned == NULL);

  /* create partitioned mangll object */
  d->mangll_partitioned = rhea_discretization_mangll_discontinuous_new (
      p4est, d->discr_options->order,
      d->discr_options->X_fn, d->discr_options->X_data);
  RHEA_ASSERT (d->mangll_partitioned != NULL);
  RHEA_ASSERT (d->mangll_partitioned->N == d->mangll_adapted->N);

  /* perform partitioning */
  {
    const int           has_temp = (d->temperature_adapted != NULL);
    const int           has_vel = (d->velocity_adapted != NULL);
    const int           has_press = (d->pressure_adapted != NULL);

    /* partition from adapted to original */
    if (has_temp) {
      RHEA_ASSERT (d->temperature_original == NULL);
      d->temperature_original = rhea_stokes_problem_amr_partition_field (
          d->temperature_adapted, 1 /* #fields*/,
          d->mangll_adapted, d->mangll_partitioned);
      d->temperature_adapted = NULL;
    }
    if (has_vel) {
      RHEA_ASSERT (d->velocity_original == NULL);
      d->velocity_original = rhea_stokes_problem_amr_partition_field (
          d->velocity_adapted, 3 /* #fields*/,
          d->mangll_adapted, d->mangll_partitioned);
      d->velocity_adapted = NULL;
    }
    if (has_press) {
      RHEA_ASSERT (d->pressure_original == NULL);
      d->pressure_original = rhea_stokes_problem_amr_partition_field (
          d->pressure_adapted, 1 /* #fields*/,
          d->mangll_adapted, d->mangll_partitioned);
      d->pressure_adapted = NULL;
    }
  }

  /* destroy adapted but unpartitioned mangll */
  rhea_discretization_mangll_interpolation_destroy (d->mangll_adapted);
  d->mangll_adapted = NULL;

  /* reassign original mangll for next projection */
  d->mangll_original = d->mangll_partitioned;
  d->mangll_partitioned = NULL;
}

/******************************************************************************
 * AMR Main Functions
 *****************************************************************************/

#if 0
int
rhea_stokes_problem_amr (p4est_t *p4est)
{
  rhea_stokes_problem_amr_data_t  data;
  int                 amr_iter;

  amr_iter = rhea_amr (p4est,
                       rhea_stokes_problem_amr_data_initialize_fn,
                       rhea_stokes_problem_amr_data_finalize_fn,
                       rhea_stokes_problem_amr_data_project_fn,
                       rhea_stokes_problem_amr_data_partition_fn,
                       &data);

  return amr_iter;
}
#endif
