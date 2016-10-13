/*
  This file is part of the ymir Library.
  ymir is a C library for modeling ice sheets

  Copyright (C) 2010, 2011 Carsten Burstedde, Toby Isaac, Georg Stadler,
                           Lucas Wilcox.

  The ymir Library is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The ymir Library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the ymir Library.  If not, see <http://www.gnu.org/licenses/>.

  ---

  This is the slabs example for global instantaneous mantle flow with plates.

*/

#include <slabs_base.h>
#include <slabs_stokes_state.h>
#include <mangll_fields.h>
#include <ymir_gmg_intergrid_extended.h>
#include <ymir_mass_vec.h>
#include <ymir_pressure_vec.h>
#include <ymir_stokes_vec.h>
#include <ymir_interp_vec.h>
#include <slabs_io.h>
#include <slabs_discretization.h>
#include <p8est_extended.h>

/**
 *
 */
slabs_stokes_state_t *
slabs_stokes_state_new (p8est_t *p8est)
{
  slabs_stokes_state_t  *state;

  /* create new Stokes state */
  state = YMIR_ALLOC (slabs_stokes_state_t, 1);

  /* initialize p8est mesh */
  state->p8est = p8est;

  /* initialize state variables */
  state->temperature = NULL;
  state->weakzone = NULL;
  state->velocity = NULL;
  state->pressure = NULL;
  state->dual_tensor = NULL;
  state->temp_vec = NULL;
  state->weak_vec = NULL;
  state->vel_press_vec = NULL;
  state->dual_vec = NULL;

  /* initialize weak zone function */
  state->set_weakzone_fn = NULL;
  state->set_weakzone_data = NULL;

  /* initialize buffer variables */
  state->velocity_bc = NULL;
  state->vel_bc_vec = NULL;
  state->amr_buffer_coarse_temp = NULL;
  state->amr_buffer_coarse_vel = NULL;
  state->amr_buffer_coarse_press = NULL;
  state->amr_buffer_coarse_dual = NULL;
  state->amr_buffer_fine_temp = NULL;
  state->amr_buffer_fine_vel = NULL;
  state->amr_buffer_fine_press = NULL;
  state->amr_buffer_fine_dual = NULL;

  /* return Stokes state */
  return state;
}

/**
 *
 */
sc_dmatrix_t *
slabs_stokes_state_init_temp (slabs_stokes_state_t *state,
                              mangll_cnodes_t *cnodes)
{
  const mangll_locidx_t  n_nodes = cnodes->Ncn;

  /* check state */
  YMIR_ASSERT (state->temperature == NULL);

  /* create temperature field */
  state->temperature = sc_dmatrix_new (n_nodes, 1);

  /* return temperature field */
  return state->temperature;
}

/**
 *
 */
ymir_vec_t *
slabs_stokes_state_init_temp_vec (slabs_stokes_state_t *state,
                                  ymir_mesh_t *mesh)
{
  /* check state */
  YMIR_ASSERT (state->temperature != NULL);
  YMIR_ASSERT (state->temp_vec == NULL);

  /* create vector view for temperature */
  state->temp_vec = ymir_cvec_new_data (mesh, 1, state->temperature);

  /* return temperature vector */
  return state->temp_vec;
}

/**
 *
 */
void
slabs_stokes_state_clear_temp (slabs_stokes_state_t *state)
{
  if (state->temperature != NULL) {
    sc_dmatrix_destroy (state->temperature);
    state->temperature = NULL;
  }
  if (state->temp_vec != NULL) {
    ymir_vec_destroy (state->temp_vec);
    state->temp_vec = NULL;
  }
}

/**
 *
 */
ymir_vec_t *
slabs_stokes_state_init_weakzone (slabs_stokes_state_t *state,
                                  ymir_mesh_t *mesh,
                                  slabs_stokes_state_set_weakzone_t set_fn,
                                  void *set_data)
{
  const mangll_locidx_t  n_elements = mesh->cnodes->K;
  const int           N = ymir_n (mesh->ma->N);
  const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);

  /* check input */
  YMIR_ASSERT (set_fn != NULL);

  /* check state */
  YMIR_ASSERT (state->weak_vec == NULL);
  YMIR_ASSERT (state->weakzone == NULL);

  /* create weak zone vector */
  state->weak_vec = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);

  /* create view onto velocity vector */
  state->weakzone = sc_dmatrix_new_view (n_elements, n_nodes_per_el,
                                         state->weak_vec->dvec);

  /* set weak zone generator function */
  state->set_weakzone_fn = set_fn;
  state->set_weakzone_data = set_data;

  /* return weak zone vector */
  return state->weak_vec;
}

/**
 *
 */
void
slabs_stokes_state_clear_weakzone (slabs_stokes_state_t *state)
{
  YMIR_ASSERT (state->weakzone != NULL);
  YMIR_ASSERT (state->weak_vec != NULL);

  sc_dmatrix_destroy (state->weakzone);
  state->weakzone = NULL;

  ymir_vec_destroy (state->weak_vec);
  state->weak_vec = NULL;
}

/**
 *
 */
static void
slabs_stokes_state_compute_weakzone (slabs_stokes_state_t *state)
{
  /* check state */
  YMIR_ASSERT (state->weak_vec != NULL);
  YMIR_ASSERT (state->weakzone != NULL);
  YMIR_ASSERT (state->set_weakzone_fn != NULL);

  /* compute weak zone factors */
  state->set_weakzone_fn (state->weak_vec, state->set_weakzone_data);
}

/**
 *
 */
ymir_vec_t *
slabs_stokes_state_init_vel_press (slabs_stokes_state_t *state,
                                   ymir_mesh_t *mesh,
                                   ymir_pressure_elem_t *press_elem)
{
  const mangll_locidx_t  n_cnodes = mesh->cnodes->Ncn;
  const mangll_locidx_t  n_elements = mesh->cnodes->K;

  /* check state */
  YMIR_ASSERT (state->vel_bc_vec == NULL);
  YMIR_ASSERT (state->velocity_bc == NULL);
  YMIR_ASSERT (state->vel_press_vec == NULL);
  YMIR_ASSERT (state->velocity == NULL);
  YMIR_ASSERT (state->pressure == NULL);
  //YMIR_ASSERT (1 < ymir_n (mesh->ma->N));

  /* create (velocity,pressure) pair vector */
  state->vel_press_vec = ymir_stokes_vec_new (mesh, press_elem);

  /* create views onto velocity and pressure vectors */
  switch (press_elem->space) {
  case YMIR_PRESSURE_SPACE_POLY:
  case YMIR_PRESSURE_SPACE_TENS:
    /* create view onto velocity vector */
    YMIR_ASSERT (state->vel_press_vec->ncfields == 3);
    state->velocity = sc_dmatrix_new_view (n_cnodes, 3,
                                           state->vel_press_vec->cvec);

    /* create view onto pressure vector
     * Warning: Requires discontinuous pressure space, otherwise pressure
     *          unknowns are stored in cvec */
    YMIR_ASSERT (state->vel_press_vec->nefields == press_elem->Np);
    state->pressure = sc_dmatrix_new_view (n_elements, press_elem->Np,
                                           state->vel_press_vec->evec);
    break;

  case YMIR_PRESSURE_SPACE_STAB:
    state->velocity = NULL;
    state->pressure = NULL;
    break;

  default: /* unknown pressure space */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* create velocity vector */
  state->vel_bc_vec = ymir_cvec_new (mesh, 3);

  /* create view onto velocity vector */
  state->velocity_bc = sc_dmatrix_new_view (n_cnodes,
                                            state->vel_bc_vec->ncfields,
                                            state->vel_bc_vec->cvec);

  /* return (velocity,pressure) vector */
  return state->vel_press_vec;
}

/**
 *
 */
void
slabs_stokes_state_clear_vel_press (slabs_stokes_state_t *state)
{
  YMIR_ASSERT (state->velocity_bc != NULL);
  YMIR_ASSERT (state->vel_bc_vec != NULL);
  YMIR_ASSERT (state->vel_press_vec != NULL);

  /* destroy velocity field */
  sc_dmatrix_destroy (state->velocity_bc);
  state->velocity_bc = NULL;
  ymir_vec_destroy (state->vel_bc_vec);
  state->vel_bc_vec = NULL;

  /* destroy (velocity,pressure) pair fields */
  if (state->velocity != NULL) {
    sc_dmatrix_destroy (state->velocity);
    state->velocity = NULL;
  }
  if (state->pressure != NULL) {
    sc_dmatrix_destroy (state->pressure);
    state->pressure = NULL;
  }
  ymir_vec_destroy (state->vel_press_vec);
  state->vel_press_vec = NULL;
}

/**
 *
 */
ymir_vec_t *
slabs_stokes_state_init_dual_tensor (slabs_stokes_state_t *state,
                                     ymir_mesh_t *mesh)
{
  const mangll_locidx_t  n_elements = mesh->cnodes->K;
  const int           N = ymir_n (mesh->ma->N);
  const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);

  /* check state */
  YMIR_ASSERT (state->dual_vec == NULL);
  YMIR_ASSERT (state->dual_tensor == NULL);

  /* create ymir vector for dual tensor */
  state->dual_vec = ymir_dvec_new (mesh, 6, YMIR_GAUSS_NODE);

  /* create view onto dual tensor */
  state->dual_tensor = sc_dmatrix_new_view (n_elements,
      n_nodes_per_el * state->dual_vec->ndfields,
      state->dual_vec->dvec);

  /* return dual tensor */
  return state->dual_vec;
}

/**
 *
 */
void
slabs_stokes_state_clear_dual_tensor (slabs_stokes_state_t *state)
{
  YMIR_ASSERT (state->dual_tensor != NULL);
  YMIR_ASSERT (state->dual_vec != NULL);

  sc_dmatrix_destroy (state->dual_tensor);
  state->dual_tensor = NULL;

  ymir_vec_destroy (state->dual_vec);
  state->dual_vec = NULL;
}

/**
 *
 */
static void
slabs_stokes_state_clear (slabs_stokes_state_t *state)
{
  /* destroy temperature field */
  slabs_stokes_state_clear_temp (state);

  /* destroy weak zone */
  if (state->weakzone != NULL && state->weak_vec != NULL) {
    slabs_stokes_state_clear_weakzone (state);
  }
#ifdef YMIR_DEBUG
  else {
    YMIR_ASSERT (state->weakzone == NULL);
    YMIR_ASSERT (state->weak_vec == NULL);
  }
#endif

  /* destroy velocity and pressure fields */
  if (state->vel_bc_vec != NULL && state->vel_press_vec != NULL) {
    slabs_stokes_state_clear_vel_press (state);
  }
#ifdef YMIR_DEBUG
  else {
    YMIR_ASSERT (state->velocity_bc == NULL);
    YMIR_ASSERT (state->vel_bc_vec == NULL);
    YMIR_ASSERT (state->velocity == NULL);
    YMIR_ASSERT (state->pressure == NULL);
    YMIR_ASSERT (state->vel_press_vec == NULL);
  }
#endif

  /* destroy dual tensor */
  if (state->dual_tensor != NULL && state->dual_vec != NULL) {
    slabs_stokes_state_clear_dual_tensor (state);
  }
#ifdef YMIR_DEBUG
  else {
    YMIR_ASSERT (state->dual_tensor == NULL);
    YMIR_ASSERT (state->dual_vec == NULL);
  }
#endif
}

/**
 *
 */
static void
slabs_stokes_state_clear_amr_buffer (slabs_stokes_state_t *state)
{
  if (state->amr_buffer_coarse_temp != NULL) {
    sc_dmatrix_destroy (state->amr_buffer_coarse_temp);
    state->amr_buffer_coarse_temp = NULL;
  }
  if (state->amr_buffer_coarse_vel != NULL) {
    sc_dmatrix_destroy (state->amr_buffer_coarse_vel);
    state->amr_buffer_coarse_vel = NULL;
  }
  if (state->amr_buffer_coarse_press != NULL) {
    sc_dmatrix_destroy (state->amr_buffer_coarse_press);
    state->amr_buffer_coarse_press = NULL;
  }
  if (state->amr_buffer_coarse_dual != NULL) {
    sc_dmatrix_destroy (state->amr_buffer_coarse_dual);
    state->amr_buffer_coarse_dual = NULL;
  }

  if (state->amr_buffer_fine_temp != NULL) {
    sc_dmatrix_destroy (state->amr_buffer_fine_temp);
    state->amr_buffer_fine_temp = NULL;
  }
  if (state->amr_buffer_fine_vel != NULL) {
    sc_dmatrix_destroy (state->amr_buffer_fine_vel);
    state->amr_buffer_fine_vel = NULL;
  }
  if (state->amr_buffer_fine_press != NULL) {
    sc_dmatrix_destroy (state->amr_buffer_fine_press);
    state->amr_buffer_fine_press = NULL;
  }
  if (state->amr_buffer_fine_dual != NULL) {
    sc_dmatrix_destroy (state->amr_buffer_fine_dual);
    state->amr_buffer_fine_dual = NULL;
  }
}

/**
 *
 */
void
slabs_stokes_state_destroy (slabs_stokes_state_t *state)
{
  /* stop execution if nothing to do */
  if (state == NULL) {
    return;
  }

  /* destroy state variables */
  slabs_stokes_state_clear (state);

  /* destroy AMR buffer */
  slabs_stokes_state_clear_amr_buffer (state);

  /* destroy state */
  YMIR_FREE (state);
}

/**
 *
 */
void
slabs_stokes_state_save (slabs_stokes_state_t *state, const char *filepath)
{
  const char         *this_fn_name = "slabs_stokes_state_save";
  MPI_Comm            mpicomm;
  int                 mpisize;
  ymir_mesh_t        *mesh;
  ymir_locidx_t      *size_loc;
  ymir_gloidx_t      *offset_glo;
  char                path[BUFSIZ];
  int                 r;

  /* check input */
  YMIR_ASSERT (state->p8est != NULL);
  YMIR_ASSERT (state->temp_vec != NULL);

  YMIR_GLOBAL_INFOF ("Into %s: \"%s\"\n", this_fn_name, filepath);

  /* save p4est */
  snprintf (path, BUFSIZ, "%s_mesh", filepath);
  //p8est_io_write (path, state->p8est, 0 /* do not write data */);
  //TODO: add that function

  /* get mesh & parallel environment */
  mesh = state->temp_vec->mesh;
  mpisize = mesh->ma->mpisize;
  mpicomm = mesh->ma->mpicomm;

  /* get number of local cnodes */
  size_loc = YMIR_ALLOC (ymir_locidx_t, mpisize);
  memcpy (size_loc, mesh->cnodes->Ngo, mpisize * sizeof (mangll_locidx_t));

  /* compute global cnode offsets */
  offset_glo = YMIR_ALLOC (ymir_gloidx_t, mpisize + 1);
  offset_glo[0] = 0;
  for (r = 0; r < mpisize; r++) { /* loop over all processors */
    offset_glo[r + 1] = offset_glo[r] + (ymir_gloidx_t) size_loc[r];
  }

  /* save temperature field */
  snprintf (path, BUFSIZ, "%s_temperature.bin", filepath);
  slabs_io_gather_write_vector (path, state->temperature->e[0],
                                offset_glo, size_loc, mpicomm);

  /* TODO weak zone? */

  /* save velocity and pressure fields */
  if (state->velocity != NULL && state->pressure != NULL) {
    int                 n_fields;
    p4est_gloidx_t     *quadrant_offset = state->p8est->global_first_quadrant;

    /* set local sizes and offsets */
    n_fields = state->velocity->n;
    for (r = 0; r < mpisize; r++) { /* loop over all processors */
      size_loc[r] *= n_fields;
      offset_glo[r] *= n_fields;
    }
    offset_glo[mpisize] *= n_fields;

    /* save velocity */
    snprintf (path, BUFSIZ, "%s_velocity.bin", filepath);
    slabs_io_gather_write_vector (path, state->velocity->e[0],
                                  offset_glo, size_loc, mpicomm);

    /* set local sizes and offsets */
    n_fields = state->pressure->n;
    offset_glo[0] = 0;
    for (r = 0; r < mpisize; r++) { /* loop over all processors */
      offset_glo[r+1] = (ymir_gloidx_t) n_fields * quadrant_offset[r+1];
      size_loc[r] = (ymir_locidx_t) (offset_glo[r+1] - offset_glo[r]);
    }

    /* save pressure */
    snprintf (path, BUFSIZ, "%s_pressure.bin", filepath);
    slabs_io_gather_write_vector (path, state->pressure->e[0],
                                  offset_glo, size_loc, mpicomm);
  }

  /* save dual tensor */
  if (state->dual_tensor != NULL) {
    int                 n_fields = state->dual_tensor->n;
    p4est_gloidx_t     *quadrant_offset = state->p8est->global_first_quadrant;

    /* set local sizes and offsets */
    offset_glo[0] = 0;
    for (r = 0; r < mpisize; r++) { /* loop over all processors */
      offset_glo[r+1] = (ymir_gloidx_t) n_fields * quadrant_offset[r+1];
      size_loc[r] = (ymir_locidx_t) (offset_glo[r+1] - offset_glo[r]);
    }

    /* save dual tensor */
    snprintf (path, BUFSIZ, "%s_dualtensor.bin", filepath);
    slabs_io_gather_write_vector (path, state->dual_tensor->e[0],
                                  offset_glo, size_loc, mpicomm);
  }

  /* destroy */
  YMIR_FREE (size_loc);
  YMIR_FREE (offset_glo);

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

/**
 *
 */
slabs_stokes_state_t *
slabs_stokes_state_load_init (p8est_t ** p8est, const char *filepath_p4est,
                              MPI_Comm mpicomm, const int inspect_p4est)
{
  slabs_stokes_state_t *state;

  /* load p4est */
  //*p8est = p8est_io_read (filepath_p4est, mpicomm,
  //                        0 /* zero data size */, 0 /* do not read data */);
  //TODO: add that function

  /* initialize user data */
  slabs_discr_p8est_init_data (*p8est, inspect_p4est);

  /* partition p4est mesh for multigrid coarsening */
  p8est_partition_ext (*p8est, 1 /* for coarsening */, NULL);

  /* create new Stokes state */
  state = slabs_stokes_state_new (*p8est);

  /* return Stokes state */
  return state;
}

void
slabs_stokes_state_load_temp (slabs_stokes_state_t *state, ymir_mesh_t *mesh,
                              const char *filepath_temperature_txt,
                              const char *filepath_temperature_bin)
{
  mangll_t           *mangll = mesh->ma;
  mangll_cnodes_t    *cnodes = mesh->cnodes;

  /* check input */
  YMIR_ASSERT (state->temperature == NULL);

  /* initialize temperature in Stokes state */
  slabs_stokes_state_init_temp (state, cnodes);

  /* read temperature from file into Stokes state */
  slabs_io_read_temperature (state->temperature, mangll, cnodes,
                             filepath_temperature_txt,
                             filepath_temperature_bin);

  /* initialize temperature vector in Stokes state */
  slabs_stokes_state_init_temp_vec (state, mesh);
}

void
slabs_stokes_state_load_vel_press (slabs_stokes_state_t *state,
                                   ymir_mesh_t *mesh,
                                   ymir_pressure_elem_t *press_elem,
                                   const char *filepath_velocity,
                                   const char *filepath_pressure)
{
  MPI_Comm            mpicomm;
  int                 mpisize;
  p4est_gloidx_t     *quadrant_offset = state->p8est->global_first_quadrant;
  ymir_locidx_t      *size_loc;
  ymir_gloidx_t      *offset_glo;
  int                 n_fields;
  int                 r;

  /* check input */
  YMIR_ASSERT (state->velocity == NULL);
  YMIR_ASSERT (state->pressure == NULL);
  YMIR_ASSERT (state->vel_press_vec == NULL);

  /* initialize velocity and pressure in Stokes state */
  slabs_stokes_state_init_vel_press (state, mesh, press_elem);

  /* get mesh & parallel environment */
  mpisize = mesh->ma->mpisize;
  mpicomm = mesh->ma->mpicomm;

  /* get number of local cnodes */
  size_loc = YMIR_ALLOC (ymir_locidx_t, mpisize);
  memcpy (size_loc, mesh->cnodes->Ngo, mpisize * sizeof (mangll_locidx_t));

  /* set local sizes and offsets for velocity */
  n_fields = state->velocity->n;
  offset_glo = YMIR_ALLOC (ymir_gloidx_t, mpisize + 1);
  offset_glo[0] = 0;
  for (r = 0; r < mpisize; r++) { /* loop over all processors */
    size_loc[r] *= n_fields;
    offset_glo[r + 1] = offset_glo[r] + (ymir_gloidx_t) size_loc[r];
  }

  /* read velocity */
  slabs_io_read_scatter_vector (filepath_velocity, state->velocity->e[0],
                                offset_glo, size_loc, mpicomm);

  /* communicate shared node values */
  {
    sc_array_t         *vel_array;

    vel_array = sc_array_new_data (state->velocity->e[0],
                                   (size_t) n_fields * sizeof (double),
                                   (size_t) mesh->cnodes->Ncn);
    mangll_cnodes_share_owned (vel_array, mesh->cnodes);
    sc_array_destroy (vel_array);
  }

  /* set local sizes and offsets for pressure */
  n_fields = state->pressure->n;
  offset_glo[0] = 0;
  for (r = 0; r < mpisize; r++) { /* loop over all processors */
    offset_glo[r+1] = (ymir_gloidx_t) n_fields * quadrant_offset[r+1];
    size_loc[r] = (ymir_locidx_t) (offset_glo[r+1] - offset_glo[r]);
  }

  /* read pressure */
  slabs_io_read_scatter_vector (filepath_pressure, state->pressure->e[0],
                                offset_glo, size_loc, mpicomm);

  /* destroy */
  YMIR_FREE (size_loc);
  YMIR_FREE (offset_glo);
}

/**
 *
 */
void
slabs_stokes_state_amr_prepare (slabs_stokes_state_t *state,
                                ymir_mesh_t *mesh)
{
  const char         *this_fn_name = "slabs_stokes_state_amr_prepare";
  const mangll_locidx_t  n_elements = mesh->cnodes->K;
  const int           N = ymir_n (mesh->ma->N);
  const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);
  ymir_vec_t         *temp_vec;
  ymir_vec_t         *vel_press_vec = state->vel_press_vec;
  ymir_vec_t         *dual_vec = state->dual_vec;

  ymir_vec_t         *vel_vec;
  ymir_vec_t         *press_vec;
  ymir_vec_t         *temp_buffer_vec;
  ymir_vec_t         *vel_buffer_vec;
  ymir_vec_t         *press_buffer_vec;
  ymir_vec_t         *dual_buffer_vec;
  int                 process_vel_press;
  int                 process_dual;

  YMIR_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  /* check state */
  YMIR_ASSERT (state->temperature != NULL);
  YMIR_ASSERT (state->amr_buffer_coarse_temp == NULL);
  YMIR_ASSERT (state->amr_buffer_coarse_vel == NULL);
  YMIR_ASSERT (state->amr_buffer_coarse_press == NULL);
  YMIR_ASSERT (state->amr_buffer_coarse_dual == NULL);

  /* decide whether to process velocity and pressure fields */
  if (vel_press_vec != NULL) {
    process_vel_press = 1;
  }
  else {
    process_vel_press = 0;
  }

  /* decide whether to process dual tensor */
  if (dual_vec != NULL) {
    process_dual = 1;
  }
  else {
    process_dual = 0;
  }

  /*
   * Setup source variables
   */

  /* create vector view for temperature if it does not exist already */
  if (state->temp_vec == NULL) {
    temp_vec = slabs_stokes_state_init_temp_vec (state, mesh);
  }
  else {
    temp_vec = state->temp_vec;
  }

  /* create vector views onto the velocity and pressure components */
  if (process_vel_press) {
    slabs_stokes_vec_get_components_view (&vel_vec, &press_vec, vel_press_vec);
  }

  /*
   * Setup destination variables
   */

  /* create buffer for discontinuous fields on GLL nodes */
  state->amr_buffer_coarse_temp = sc_dmatrix_new (n_elements, n_nodes_per_el);
  if (process_vel_press) {
    state->amr_buffer_coarse_vel = sc_dmatrix_new (n_elements,
                                                   3 * n_nodes_per_el);
    state->amr_buffer_coarse_press = sc_dmatrix_new (n_elements,
                                                     n_nodes_per_el);
  }
  if (process_dual) {
    state->amr_buffer_coarse_dual = sc_dmatrix_new (n_elements,
                                                    6 * n_nodes_per_el);
  }

  /* create vector views of buffer fields */
  temp_buffer_vec = ymir_dvec_new_data (mesh, 1, YMIR_GLL_NODE,
                                        state->amr_buffer_coarse_temp);
  if (process_vel_press) {
    vel_buffer_vec = ymir_dvec_new_data (mesh, 3, YMIR_GLL_NODE,
                                         state->amr_buffer_coarse_vel);
    press_buffer_vec = ymir_dvec_new_data (mesh, 1, YMIR_GAUSS_NODE,
                                           state->amr_buffer_coarse_press);
  }
  if (process_dual) {
    dual_buffer_vec = ymir_dvec_new_data (mesh, 6, YMIR_GAUSS_NODE,
                                          state->amr_buffer_coarse_dual);
  }

  /*
   * Interpolate from source to destination variables
   */

  /* interpolate/copy state fields onto discontinuous nodes */
  ymir_interp_vec (temp_vec, temp_buffer_vec);
  if (process_vel_press) {
    ymir_interp_vec (vel_vec, vel_buffer_vec);
    ymir_interp_vec (press_vec, press_buffer_vec);
  }
  if (process_dual) {
    ymir_dvec_copy (dual_vec, dual_buffer_vec);
  }

  /* destroy work variables */
  ymir_vec_destroy (temp_buffer_vec);
  if (process_vel_press) {
    ymir_vec_destroy (vel_vec);
    ymir_vec_destroy (press_vec);
    ymir_vec_destroy (vel_buffer_vec);
    ymir_vec_destroy (press_buffer_vec);
  }
  if (process_dual) {
    ymir_vec_destroy (dual_buffer_vec);
  }

  /* destroy state variables that are not necessary for AMR */
  slabs_stokes_state_clear (state);

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

/**
 *
 */
void
slabs_stokes_state_amr_finalize (slabs_stokes_state_t *state,
                                 ymir_mesh_t *mesh,
                                 ymir_pressure_elem_t *press_elem)
{
  const char         *this_fn_name = "slabs_stokes_state_amr_finalize";
  mangll_cnodes_t    *cnodes = mesh->cnodes;
  ymir_vec_t         *temp_vec;
  ymir_vec_t         *vel_press_vec;
  ymir_vec_t         *dual_vec;

  ymir_vec_t         *temp_buffer_vec, *temp_buffer_vec_mass;
  ymir_vec_t         *vel_buffer_vec = NULL, *vel_buffer_vec_mass = NULL;
  ymir_vec_t         *press_buffer_vec = NULL, *press_buffer_vec_mass = NULL;
  ymir_vec_t         *dual_buffer_vec = NULL;
  ymir_vec_t         *temp_mass;
  ymir_vec_t         *vel_vec, *vel_mass = NULL;
  ymir_vec_t         *press_vec, *press_lump = NULL;
  int                 process_vel_press;
  int                 process_dual;

  YMIR_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  /* check state */
  YMIR_ASSERT (state->amr_buffer_coarse_temp != NULL);

  /* decide whether to process velocity and pressure fields */
  if (   state->amr_buffer_coarse_vel != NULL
      && state->amr_buffer_coarse_press != NULL ) {
    process_vel_press = 1;
  }
  else if (   state->amr_buffer_coarse_vel == NULL
           && state->amr_buffer_coarse_press == NULL ) {
    process_vel_press = 0;
  }
  else {
    YMIR_ABORT_NOT_REACHED ();
  }

  /* decide whether to process dual tensor */
  if (state->amr_buffer_coarse_dual != NULL) {
    process_dual = 1;
  }
  else {
    process_dual = 0;
  }

  /*
   * Setup source variables
   */

  /* create vector views of buffer fields */
  temp_buffer_vec = ymir_dvec_new_data (mesh, 1, YMIR_GLL_NODE,
                                        state->amr_buffer_coarse_temp);
  if (process_vel_press) {
    vel_buffer_vec = ymir_dvec_new_data (mesh, 3, YMIR_GLL_NODE,
                                         state->amr_buffer_coarse_vel);
    press_buffer_vec = ymir_dvec_new_data (mesh, 1, YMIR_GAUSS_NODE,
                                           state->amr_buffer_coarse_press);
  }
  if (process_dual) {
    dual_buffer_vec = ymir_dvec_new_data (mesh, 6, YMIR_GAUSS_NODE,
                                          state->amr_buffer_coarse_dual);
  }

  /* apply mass matrix before interpolating to native nodes */
  temp_buffer_vec_mass = ymir_dvec_new (mesh, 1, YMIR_GLL_NODE);
  ymir_mass_apply (temp_buffer_vec, temp_buffer_vec_mass);
  if (process_vel_press) {
    vel_buffer_vec_mass = ymir_dvec_new (mesh, 3, YMIR_GLL_NODE);
    ymir_mass_apply (vel_buffer_vec, vel_buffer_vec_mass);
    press_buffer_vec_mass = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
    ymir_mass_apply (press_buffer_vec, press_buffer_vec_mass);
  }
  /* Note: dual does not require mass matrix or inv. mass matrix */

  /*
   * Setup destination variables
   */

  /* initialize temperature field */
  slabs_stokes_state_init_temp (state, cnodes);
  temp_vec = slabs_stokes_state_init_temp_vec (state, mesh);

  /* initialize velocity and pressure fields */
  if (process_vel_press) {
    slabs_stokes_state_init_vel_press (state, mesh, press_elem);
    vel_press_vec = state->vel_press_vec;
    slabs_stokes_vec_get_components_view (&vel_vec, &press_vec, vel_press_vec);
  }

  /* initialize dual tensor */
  if (process_dual) {
    slabs_stokes_state_init_dual_tensor (state, mesh);
    dual_vec = state->dual_vec;
  }

  /* create additional temporary vectors */
  temp_mass = ymir_vec_template (temp_vec);
  if (process_vel_press) {
    vel_mass = ymir_vec_template (vel_vec);
    press_lump = ymir_pressure_vec_new (mesh, press_elem);
    ymir_pressure_vec_lump_mass (press_lump, press_elem);
  }

  /*
   * Interpolate from source to destination variables
   */

  /* interpolate/copy state fields on discont. nodes onto their native nodes */
  ymir_interp_vec (temp_buffer_vec_mass, temp_mass);
  if (process_vel_press) {
    ymir_interp_vec (vel_buffer_vec_mass, vel_mass);
    ymir_interp_vec (press_buffer_vec_mass, press_vec);
  }
  if (process_dual) {
    ymir_dvec_copy (dual_buffer_vec, dual_vec);
  }

  /* apply inverse mass matrix after interpolating to native nodes */
  ymir_mass_invert (temp_mass, temp_vec);
  if (process_vel_press) {
    ymir_mass_invert (vel_mass, vel_vec);
    ymir_vec_divide_in (press_lump, press_vec);
  }
  /* Note: dual does not require mass matrix or inv. mass matrix */

  /* destroy work variables */
  ymir_vec_destroy (temp_buffer_vec_mass);
  ymir_vec_destroy (temp_buffer_vec);
  ymir_vec_destroy (temp_mass);
  if (process_vel_press) {
    ymir_vec_destroy (vel_buffer_vec_mass);
    ymir_vec_destroy (vel_buffer_vec);
    ymir_vec_destroy (press_buffer_vec_mass);
    ymir_vec_destroy (press_buffer_vec);
    ymir_vec_destroy (vel_vec);
    ymir_vec_destroy (vel_mass);
    ymir_vec_destroy (press_vec);
    ymir_vec_destroy (press_lump);
  }
  if (process_dual) {
    ymir_vec_destroy (dual_buffer_vec);
  }

  /* destroy AMR buffer */
  slabs_stokes_state_clear_amr_buffer (state);

  /*
   * Postprocessing
   */

  /* compute weak zone factor */
  if (state->set_weakzone_fn != NULL) {
    slabs_stokes_state_init_weakzone (state, mesh, state->set_weakzone_fn,
                                      state->set_weakzone_data);
    slabs_stokes_state_compute_weakzone (state);
  }

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

/**
 *
 */
void
slabs_stokes_state_project (slabs_stokes_state_t *state,
                            mangll_t *mangll_source, mangll_t *mangll_dest)
{
  const char         *this_fn_name = "slabs_stokes_state_project";
  const mangll_locidx_t  n_elements_dest = mangll_dest->mesh->K;
  const int           N = ymir_n (mangll_source->N);
  const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);

  int                 process_vel_press;
  int                 process_dual;
  sc_dmatrix_t       *Irm_gauss, *Irp_gauss;

  YMIR_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  /* check input */
  YMIR_ASSERT (mangll_source->N == mangll_dest->N);

  /* check state */
  YMIR_ASSERT (state->amr_buffer_coarse_temp != NULL);
  YMIR_ASSERT (state->amr_buffer_fine_temp == NULL);
  YMIR_ASSERT (state->amr_buffer_fine_vel == NULL);
  YMIR_ASSERT (state->amr_buffer_fine_press == NULL);
  YMIR_ASSERT (state->amr_buffer_fine_dual == NULL);

  /* decide whether to process velocity and pressure fields */
  if (   state->amr_buffer_coarse_vel != NULL
      && state->amr_buffer_coarse_press != NULL ) {
    process_vel_press = 1;
  }
  else if (   state->amr_buffer_coarse_vel == NULL
           && state->amr_buffer_coarse_press == NULL ) {
    process_vel_press = 0;
  }
  else {
    YMIR_ABORT_NOT_REACHED ();
  }

  /* decide whether to process dual tensor */
  if (state->amr_buffer_coarse_dual != NULL) {
    process_dual = 1;
  }
  else {
    process_dual = 0;
  }

  /* create interpolation matrices for Gauss nodes */
  if (process_vel_press || process_dual) {
    ymir_gmg_intergrid_gauss_to_gauss_interp_matrix (
        &Irm_gauss, &Irp_gauss, mangll_source);
  }

  /* setup destination fields: create buffer for discont. fields on GLL nodes */
  state->amr_buffer_fine_temp = sc_dmatrix_new (n_elements_dest,
                                                n_nodes_per_el);
  if (process_vel_press) {
    state->amr_buffer_fine_vel = sc_dmatrix_new (n_elements_dest,
                                                 3 * n_nodes_per_el);
    state->amr_buffer_fine_press = sc_dmatrix_new (n_elements_dest,
                                                   n_nodes_per_el);
  }
  if (process_dual) {
    state->amr_buffer_fine_dual = sc_dmatrix_new (n_elements_dest,
                                                  6 * n_nodes_per_el);
  }

  /* project all state fields */
  ymir_gmg_intergrid_project_OLD (
      1, mangll_source, state->amr_buffer_coarse_temp,
      mangll_dest, state->amr_buffer_fine_temp,
      mangll_dest->Irm, mangll_dest->Irp);
  if (process_vel_press) {
    ymir_gmg_intergrid_project_OLD (
        3, mangll_source, state->amr_buffer_coarse_vel,
        mangll_dest, state->amr_buffer_fine_vel,
        mangll_dest->Irm, mangll_dest->Irp);
    ymir_gmg_intergrid_project_OLD (
        1, mangll_source, state->amr_buffer_coarse_press,
        mangll_dest, state->amr_buffer_fine_press, Irm_gauss, Irp_gauss);
  }
  if (process_dual) {
    ymir_gmg_intergrid_project_OLD (
        6, mangll_source, state->amr_buffer_coarse_dual,
        mangll_dest, state->amr_buffer_fine_dual, Irm_gauss, Irp_gauss);
  }

  //TODO interpolation of coarsend continuous fields
  // (above works only for discont.)

  /* destroy buffer for source fields */
  sc_dmatrix_destroy (state->amr_buffer_coarse_temp);
  state->amr_buffer_coarse_temp = NULL;
  if (process_vel_press) {
    sc_dmatrix_destroy (state->amr_buffer_coarse_vel);
    state->amr_buffer_coarse_vel = NULL;
    sc_dmatrix_destroy (state->amr_buffer_coarse_press);
    state->amr_buffer_coarse_press = NULL;
  }
  if (process_dual) {
    sc_dmatrix_destroy (state->amr_buffer_coarse_dual);
    state->amr_buffer_coarse_dual = NULL;
  }
  if (process_vel_press || process_dual) {
    sc_dmatrix_destroy (Irm_gauss);
    sc_dmatrix_destroy (Irp_gauss);
  }

  /* bound temperature to valid interval [0, 1] */
  slabs_matrix_bound_values (state->amr_buffer_fine_temp, 0.0, 1.0);

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

/**
 *
 */
void
slabs_stokes_state_partition (slabs_stokes_state_t *state,
                              mangll_t *mangll_source, mangll_t *mangll_dest)
{
  const char         *this_fn_name = "slabs_stokes_state_partition";
  MPI_Comm            mpicomm = mangll_source->mpicomm;
  const int           mpisize = mangll_source->mpisize;
  const int           mpirank = mangll_source->mpirank;

  mangll_gloidx_t    *RtoGEO_source = mangll_source->mesh->RtoGEO;
  mangll_gloidx_t    *RtoGEO_dest = mangll_dest->mesh->RtoGEO;
  const mangll_locidx_t  n_elements_dest = mangll_dest->mesh->K;
  const int           N = ymir_n (mangll_source->N);
  const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);

  int                 process_vel_press;
  int                 process_dual;

  YMIR_GLOBAL_INFOF ("Into %s\n", this_fn_name);

  /* check input */
  YMIR_ASSERT (mangll_source->N == mangll_dest->N);

  /* check state */
  YMIR_ASSERT (state->amr_buffer_coarse_temp == NULL);
  YMIR_ASSERT (state->amr_buffer_coarse_vel == NULL);
  YMIR_ASSERT (state->amr_buffer_coarse_press == NULL);
  YMIR_ASSERT (state->amr_buffer_coarse_dual == NULL);
  YMIR_ASSERT (state->amr_buffer_fine_temp != NULL);

  /* decide whether to process velocity and pressure fields */
  if (   state->amr_buffer_fine_vel != NULL
      && state->amr_buffer_fine_press != NULL ) {
    process_vel_press = 1;
  }
  else if (   state->amr_buffer_fine_vel == NULL
           && state->amr_buffer_fine_press == NULL ) {
    process_vel_press = 0;
  }
  else {
    YMIR_ABORT_NOT_REACHED ();
  }

  /* decide whether to process dual tensor */
  if (state->amr_buffer_fine_dual != NULL) {
    process_dual = 1;
  }
  else {
    process_dual = 0;
  }

  /* setup destination fields: create buffer for discont. fields on GLL nodes */
  state->amr_buffer_coarse_temp = sc_dmatrix_new (n_elements_dest,
                                                  n_nodes_per_el);
  if (process_vel_press) {
    state->amr_buffer_coarse_vel = sc_dmatrix_new (n_elements_dest,
                                                   3 * n_nodes_per_el);
    state->amr_buffer_coarse_press = sc_dmatrix_new (n_elements_dest,
                                                     n_nodes_per_el);
  }
  if (process_dual) {
    state->amr_buffer_coarse_dual = sc_dmatrix_new (n_elements_dest,
                                                    6 * n_nodes_per_el);
  }

  /* partition all state fields */
  mangll_field_partition (1, n_nodes_per_el, mpirank, mpisize, mpicomm,
                          RtoGEO_source, state->amr_buffer_fine_temp,
                          RtoGEO_dest, state->amr_buffer_coarse_temp);
  if (process_vel_press) {
    mangll_field_partition (3, n_nodes_per_el, mpirank, mpisize, mpicomm,
                            RtoGEO_source, state->amr_buffer_fine_vel,
                            RtoGEO_dest, state->amr_buffer_coarse_vel);
    mangll_field_partition (1, n_nodes_per_el, mpirank, mpisize, mpicomm,
                            RtoGEO_source, state->amr_buffer_fine_press,
                            RtoGEO_dest, state->amr_buffer_coarse_press);
  }
  if (process_dual) {
    mangll_field_partition (6, n_nodes_per_el, mpirank, mpisize, mpicomm,
                            RtoGEO_source, state->amr_buffer_fine_dual,
                            RtoGEO_dest, state->amr_buffer_coarse_dual);
  }

  /* destroy buffer for source fields */
  sc_dmatrix_destroy (state->amr_buffer_fine_temp);
  state->amr_buffer_fine_temp = NULL;
  if (process_vel_press) {
    sc_dmatrix_destroy (state->amr_buffer_fine_vel);
    state->amr_buffer_fine_vel = NULL;
    sc_dmatrix_destroy (state->amr_buffer_fine_press);
    state->amr_buffer_fine_press = NULL;
  }
  if (process_dual) {
    sc_dmatrix_destroy (state->amr_buffer_fine_dual);
    state->amr_buffer_fine_dual = NULL;
  }

  YMIR_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

