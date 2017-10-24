/*
 */

#include <rhea_stokes_problem_amr.h>
#include <rhea_discretization.h>
#include <rhea_amr.h>
#include <rhea_base.h>

typedef struct rhea_stokes_problem_amr_data
{
  ymir_mesh_t        *ymir_mesh;
  ymir_pressure_elem_t *press_elem;
  rhea_discretization_options_t *discr_options;

  /* AMR iterations */
  int                 first_amr;
  mangll_t           *mangll_original;
  mangll_t           *mangll_adapted;
  mangll_t           *mangll_partitioned;

  /* fields */
  sc_dmatrix_t       *temperature_original;
  sc_dmatrix_t       *temperature_adapted;

  sc_dmatrix_t       *weakzone_original;
  sc_dmatrix_t       *weakzone_adapted;

  sc_dmatrix_t       *velocity_original;
  sc_dmatrix_t       *velocity_adapted;

  sc_dmatrix_t       *pressure_original;
  sc_dmatrix_t       *pressure_adapted;
}
rhea_stokes_problem_amr_data_t;

static void
rhea_stokes_problem_amr_data_initialize_fn (p4est_t *p4est, void *data)
{
  rhea_stokes_problem_amr_data_t *d = data;
  mangll_t           *mangll = NULL; //TODO

  /* set fields */
  //TODO

  /* set AMR parameters */
  d->first_amr = 1;
  d->mangll_original = mangll;
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

  /* check input */
  RHEA_ASSERT (d->mangll_original != NULL);
  RHEA_ASSERT (d->mangll_adapted == NULL);
  RHEA_ASSERT (d->mangll_partitioned == NULL);

  /* destroy discontinuous mangll */
  RHEA_ASSERT (!d->first_amr);
  rhea_discretization_mangll_discontinuous_destroy (d->mangll_original);
  d->mangll_original = NULL;

  /* create ymir objects */
  rhea_discretization_ymir_mesh_new_from_p4est (
      &(d->ymir_mesh), &(d->press_elem), p4est, d->discr_options);
}

static void
rhea_stokes_problem_amr_data_project_fn (p4est_t *p4est, void *data)
{
  rhea_stokes_problem_amr_data_t *d = data;

  /* check input */
  RHEA_ASSERT (d->mangll_original != NULL);
  RHEA_ASSERT (d->mangll_adapted == NULL);
  RHEA_ASSERT (d->mangll_partitioned == NULL);

  /* create adapted mangll object */
  d->mangll_adapted = rhea_discretization_mangll_interpolation_new (
      p4est, d->discr_options->order);
  RHEA_ASSERT (d->mangll_adapted != NULL);

  /* project fields */
  //TODO

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

  /* partition fields */
  //TODO

  /* destroy adapted but unpartitioned mangll */
  rhea_discretization_mangll_interpolation_destroy (d->mangll_adapted);
  d->mangll_adapted = NULL;

  /* reassign original mangll for next projection */
  d->mangll_original = d->mangll_partitioned;
  d->mangll_partitioned = NULL;
}

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
