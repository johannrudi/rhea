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

#ifndef SLABS_STOKES_STATE_H
#define SLABS_STOKES_STATE_H

#include <p8est.h>
#include <ymir_mesh.h>
#include <ymir_pressure_elem.h>
#include <ymir_vec_ops.h>

/* prototype for functions that set weak zone distance values */
typedef void (*slabs_stokes_state_set_weakzone_t) (ymir_vec_t *weak_vec,
                                                   void *);

/* fields that describe the state of a Stokes problem */
typedef struct slabs_stokes_state
{
  /* temperature, size: [number of continuous nodes] x [1] */
  sc_dmatrix_t       *temperature;

  /* weak zone, size: [number of elements] x [number of nodes per element] */
  sc_dmatrix_t       *weakzone;

  /* velocity, size: [number of continuous nodes] x [3] */
  sc_dmatrix_t       *velocity;

  /* pressure, size: [number of elements] x [number of nodes per element] */
  sc_dmatrix_t       *pressure;

  /* dual tensor on Gauss nodes (for primal-dual nonlinear solver),
   * size: [number of elements] x [number of nodes per element * 6] */
  sc_dmatrix_t       *dual_tensor;

  /* the variables above are associated with the following vectors */
  ymir_vec_t         *temp_vec;
  ymir_vec_t         *weak_vec;
  ymir_vec_t         *vel_press_vec;
  ymir_vec_t         *dual_vec;

  /* discretization */
  p8est_t            *p8est;

  /* function that sets weak zone distance values and data */
  slabs_stokes_state_set_weakzone_t  set_weakzone_fn;
  void               *set_weakzone_data;

  /* buffer for velocity with enforced (Dirichlet) boundary conditions */
  sc_dmatrix_t       *velocity_bc;
  ymir_vec_t         *vel_bc_vec;

  /* buffers for AMR,
   * size: [number of elements] x [number of nodes per element * num fields] */
  sc_dmatrix_t       *amr_buffer_coarse_temp;
  sc_dmatrix_t       *amr_buffer_coarse_vel;
  sc_dmatrix_t       *amr_buffer_coarse_press;
  sc_dmatrix_t       *amr_buffer_coarse_dual;
  sc_dmatrix_t       *amr_buffer_fine_temp;
  sc_dmatrix_t       *amr_buffer_fine_vel;
  sc_dmatrix_t       *amr_buffer_fine_press;
  sc_dmatrix_t       *amr_buffer_fine_dual;
}
slabs_stokes_state_t;

/**
 *
 */
slabs_stokes_state_t *
slabs_stokes_state_new (p8est_t *p8est);

/**
 *
 */
sc_dmatrix_t *
slabs_stokes_state_init_temp (slabs_stokes_state_t *state,
                              mangll_cnodes_t *cnodes);

/**
 *
 */
ymir_vec_t *
slabs_stokes_state_init_temp_vec (slabs_stokes_state_t *state,
                                  ymir_mesh_t *mesh);

/**
 *
 */
void
slabs_stokes_state_clear_temp (slabs_stokes_state_t *state);

/**
 *
 */
ymir_vec_t *
slabs_stokes_state_init_weakzone (slabs_stokes_state_t *state,
                                  ymir_mesh_t *mesh,
                                  slabs_stokes_state_set_weakzone_t set_fn,
                                  void *set_data);

/**
 *
 */
void
slabs_stokes_state_clear_weakzone (slabs_stokes_state_t *state);

/**
 *
 */
ymir_vec_t *
slabs_stokes_state_init_vel_press (slabs_stokes_state_t *state,
                                   ymir_mesh_t *mesh,
                                   ymir_pressure_elem_t *press_elem);

/**
 *
 */
void
slabs_stokes_state_clear_vel_press (slabs_stokes_state_t *state);

/**
 *
 */
ymir_vec_t *
slabs_stokes_state_init_dual_tensor (slabs_stokes_state_t *state,
                                     ymir_mesh_t *mesh);

/**
 *
 */
void
slabs_stokes_state_clear_dual_tensor (slabs_stokes_state_t *state);

/**
 *
 */
void
slabs_stokes_state_destroy (slabs_stokes_state_t *state);

/**
 *
 */
void
slabs_stokes_state_save (slabs_stokes_state_t *state, const char *filepath);

/**
 *
 */
slabs_stokes_state_t *
slabs_stokes_state_load_init (p8est_t ** p8est, const char *filepath_p4est,
                              MPI_Comm mpicomm, const int inspect_p4est);

void
slabs_stokes_state_load_temp (slabs_stokes_state_t *state, ymir_mesh_t *mesh,
                              const char *filepath_temperature_txt,
                              const char *filepath_temperature_bin);

void
slabs_stokes_state_load_vel_press (slabs_stokes_state_t *state,
                                   ymir_mesh_t *mesh,
                                   ymir_pressure_elem_t *press_elem,
                                   const char *filepath_velocity,
                                   const char *filepath_pressure);

/**
 *
 */
//void
//slabs_stokes_state_get (ymir_cvec_t *temp_vec,
//                        ymir_cvec_t *vel_vec,
//                        ymir_evec_t *press_vec,
//                        slabs_stokes_state_t *state);

/**
 *
 */
//void
//slabs_stokes_state_set (slabs_stokes_state_t *state,
//                        ymir_cvec_t *temp_vec,
//                        ymir_cvec_t *vel_vec,
//                        ymir_evec_t *press_vec);

/**
 *
 */
void
slabs_stokes_state_amr_prepare (slabs_stokes_state_t *state,
                                ymir_mesh_t *mesh);

/**
 *
 */
void
slabs_stokes_state_amr_finalize (slabs_stokes_state_t *state,
                                 ymir_mesh_t *mesh,
                                 ymir_pressure_elem_t *press_elem);

/**
 *
 */
void
slabs_stokes_state_project (slabs_stokes_state_t *state,
                            mangll_t *mangll_source, mangll_t *mangll_dest);

/**
 *
 */
void
slabs_stokes_state_partition (slabs_stokes_state_t *state,
                              mangll_t *mangll_source, mangll_t *mangll_dest);

#endif /* SLABS_STOKES_STATE_H */
