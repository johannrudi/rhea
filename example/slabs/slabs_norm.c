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

#include <slabs_norm.h>
#include <ymir_mass_vec.h>
#include <ymir_velocity_vec.h>
#include <ymir_pressure_vec.h>
#include <ymir_stokes_vec.h>
#include <ymir_comm.h>
#include <slabs_physics.h>
#include <slabs_discretization.h>

/**
 *
 */
static void
slabs_norm_vec_L2_innerprod_comp (double *innerprod_cvec,
                                  double *innerprod_dvec,
                                  double *innerprod_evec,
                                  ymir_vec_t *vec_l, ymir_vec_t *vec_r)
{
  ymir_vec_t         *Mvec_r = ymir_vec_template (vec_r);

  /* apply mass matrix */
  ymir_mass_apply (vec_r, Mvec_r);

  /* compute inner product for each component */
  if (innerprod_cvec != NULL) {
    *innerprod_cvec = ymir_cvec_innerprod (vec_l, Mvec_r);
  }
  if (innerprod_dvec != NULL) {
    *innerprod_dvec = ymir_dvec_innerprod (vec_l, Mvec_r);
  }
  if (innerprod_evec != NULL) {
    *innerprod_evec = ymir_evec_innerprod (vec_l, Mvec_r);
  }

  /* destroy */
  ymir_vec_destroy (Mvec_r);
}

/**
 *
 */
static void
slabs_norm_fnc_L2_innerprod_comp (double *innerprod_cvec,
                                  double *innerprod_dvec,
                                  double *innerprod_evec,
                                  ymir_vec_t *vec_l, ymir_vec_t *vec_r,
                                  ymir_pressure_elem_t *press_elem)
{
  ymir_mesh_t        *mesh = vec_l->mesh;

  /* compute inner product for continuous fields */
  if (innerprod_cvec != NULL && 0 < vec_l->ncfields) {
    ymir_cvec_t        *clump = ymir_cvec_new (mesh, vec_l->ncfields);
    ymir_cvec_t        *cvec_l_lump = ymir_cvec_new (mesh, vec_l->ncfields);

    /* divide in lumped mass matrix */
    ymir_cvec_copy (vec_l, cvec_l_lump);
    ymir_mass_lump (clump);
    ymir_cvec_divide_in (clump, cvec_l_lump);

    /* compute inner product */
    *innerprod_cvec = ymir_cvec_innerprod (cvec_l_lump, vec_r);

    /* destroy */
    ymir_vec_destroy (clump);
    ymir_vec_destroy (cvec_l_lump);
  }
  else if (innerprod_cvec != NULL) {
    *innerprod_cvec = 0.0;
  }

  /* compute inner product for discontinuous fields */
  if (innerprod_dvec != NULL && 0 < vec_l->ndfields) {
    ymir_dvec_t        *dlump = ymir_dvec_new (mesh, vec_l->ndfields,
                                               vec_l->node_type);
    ymir_dvec_t        *dvec_l_lump = ymir_dvec_new (mesh, vec_l->ndfields,
                                                     vec_l->node_type);

    /* divide in lumped mass matrix */
    ymir_dvec_copy (vec_l, dvec_l_lump);
    ymir_mass_lump (dlump);
    ymir_dvec_divide_in (dlump, dvec_l_lump);

    /* compute inner product */
    *innerprod_dvec = ymir_dvec_innerprod (dvec_l_lump, vec_r);

    /* destroy */
    ymir_vec_destroy (dlump);
    ymir_vec_destroy (dvec_l_lump);
  }
  else if (innerprod_dvec != NULL) {
    *innerprod_dvec = 0.0;
  }

  /* compute inner product for element fields */
  if (innerprod_evec != NULL && 0 < vec_l->nefields && press_elem != NULL) {
    ymir_vec_t         *elump = ymir_pressure_vec_new (mesh, press_elem);
    ymir_vec_t         *evec_l_lump = ymir_pressure_vec_new (mesh, press_elem);

    /* divide in lumped mass matrix */
    ymir_evec_copy (vec_l, evec_l_lump);
    ymir_pressure_vec_lump_mass (elump, press_elem);
    ymir_evec_divide_in (elump, evec_l_lump);

    /* compute inner product */
    *innerprod_evec = ymir_evec_innerprod (evec_l_lump, vec_r);

    /* destroy */
    ymir_vec_destroy (elump);
    ymir_vec_destroy (evec_l_lump);
  }
  else if (innerprod_evec != NULL) {
    YMIR_ASSERT (0 == vec_l->nefields && press_elem == NULL);
    *innerprod_evec = 0.0;
  }
}

/**
 *
 */
static void
slabs_norm_Hminus1_L2_innerprod_comp (double *innerprod_vel,
                                      double *innerprod_press,
                                      ymir_vec_t *up_l, ymir_vec_t *up_r,
                                      ymir_pressure_elem_t *press_elem,
                                      ymir_Hminus1_norm_op_t *norm_op)
{
  /* compute H^-1 inner product of velocity */
  *innerprod_vel = ymir_Hminus1_norm_innerprod (up_l, up_r, norm_op);

  /* compute L^2 inner product of pressure */
  slabs_norm_fnc_L2_innerprod_comp (NULL, NULL, innerprod_press, up_l, up_r,
                                    press_elem);
}

/**
 *
 */
void
slabs_norm_innerprod_comp (double *innerprod_vel, double *innerprod_press,
                           ymir_vec_t *up_l, ymir_vec_t *up_r,
                           ymir_pressure_elem_t *press_elem,
                           slabs_norm_type_t norm_type,
                           ymir_Hminus1_norm_op_t *norm_op)
{
  YMIR_ASSERT (innerprod_vel != NULL);
  YMIR_ASSERT (innerprod_press != NULL);

  switch (norm_type) {
  case SL_NORM_VEC_L2: /* l^2 inner product of two vectors */
    *innerprod_vel = ymir_cvec_innerprod (up_l, up_r);
    *innerprod_press = ymir_evec_innerprod (up_l, up_r);
    break;

  case SL_NORM_VEC_L2_MASS_WEIGHTED: /* L^2 i.p. by weighting with mass mat */
    slabs_norm_vec_L2_innerprod_comp (innerprod_vel, NULL, innerprod_press,
                                      up_l, up_r);
    break;

  case SL_NORM_FNC_L2: /* L^2 i.p. of two vectors which are mass weighted */
    slabs_norm_fnc_L2_innerprod_comp (innerprod_vel, NULL, innerprod_press,
                                      up_l, up_r, press_elem);
    break;

  case SL_NORM_FNC_HMINUS1: /* H^-1 inner prod of two continuous vectors */
    YMIR_ASSERT (norm_op != NULL);
    *innerprod_vel = ymir_Hminus1_norm_innerprod (up_l, up_r, norm_op);
    *innerprod_press = 0.0;
    break;

  case SL_NORM_FNC_HMINUS1_L2: /* (H^-1,L^2) i.p. of two vectors */
    YMIR_ASSERT (norm_op != NULL);
    slabs_norm_Hminus1_L2_innerprod_comp (innerprod_vel, innerprod_press,
                                          up_l, up_r, press_elem, norm_op);
    break;

  default: /* unknown norm type */
    YMIR_ABORT_NOT_REACHED ();
  }
}

/**
 *
 */
double
slabs_norm (double *norm_vel, double *norm_press,
            ymir_vec_t *up, ymir_pressure_elem_t *press_elem,
            slabs_norm_type_t norm_type, ymir_Hminus1_norm_op_t *norm_op)
{
  double              norm, nrm_vel, nrm_press;

  /* compute norm */
  slabs_norm_innerprod_comp (&nrm_vel, &nrm_press, up, up, press_elem,
                             norm_type, norm_op);
  norm = sqrt (nrm_vel + nrm_press);

  /* set norm of components */
  if (norm_vel != NULL) {
    *norm_vel = sqrt (nrm_vel);
  }
  if (norm_press != NULL) {
    *norm_press = sqrt (nrm_press);
  }

  /* return combined norm */
  return norm;
}

/**
 *
 */
void
slabs_norm_symtensor_innerprod_frobenius_elem (const sc_dmatrix_t *tensor_l,
                                               const sc_dmatrix_t *tensor_r,
                                               sc_dmatrix_t *innerprod)
{
  const ymir_locidx_t  n_nodes = tensor_l->m;
  ymir_locidx_t       nodeid;
  double             *innerprod_data = innerprod->e[0];

  YMIR_ASSERT (tensor_l->m == tensor_r->m && tensor_l->m == innerprod->m);
  YMIR_ASSERT (tensor_l->n == 6);
  YMIR_ASSERT (tensor_r->n == 6);
  YMIR_ASSERT (innerprod->n == 1);

  for (nodeid = 0; nodeid < n_nodes; nodeid++) { /* loop over all nodes */
    const double        *tensor_l_node = tensor_l->e[0] + 6 * nodeid;
    const double        *tensor_r_node = tensor_r->e[0] + 6 * nodeid;

    /* compute inner product for this node */
    innerprod_data[nodeid] =   tensor_l_node[0] * tensor_r_node[0]
                             + 2.0 * tensor_l_node[1] * tensor_r_node[1]
                             + 2.0 * tensor_l_node[2] * tensor_r_node[2]
                             + tensor_l_node[3] * tensor_r_node[3]
                             + 2.0 * tensor_l_node[4] * tensor_r_node[4]
                             + tensor_l_node[5] * tensor_r_node[5];
  }
}

/**
 *
 */
void
slabs_norm_symtensor_frobenius_elem (const sc_dmatrix_t *tensor,
                                     sc_dmatrix_t *norm)
{
  const ymir_locidx_t  n_nodes = tensor->m;
  ymir_locidx_t       nodeid;
  double             *norm_data = norm->e[0];

  YMIR_ASSERT (tensor->m == norm->m);
  YMIR_ASSERT (tensor->n == 6);
  YMIR_ASSERT (norm->n == 1);

  for (nodeid = 0; nodeid < n_nodes; nodeid++) { /* loop over all nodes */
    const double        *tensor_node = tensor->e[0] + 6 * nodeid;

    /* compute norm for this node */
    norm_data[nodeid] = sqrt (  tensor_node[0] * tensor_node[0]
                              + 2.0 * tensor_node[1] * tensor_node[1]
                              + 2.0 * tensor_node[2] * tensor_node[2]
                              + tensor_node[3] * tensor_node[3]
                              + 2.0 * tensor_node[4] * tensor_node[4]
                              + tensor_node[5] * tensor_node[5] );
  }
}

/**
 *
 */
double
slabs_norm_symtensor_frobenius_min (const ymir_dvec_t *tensor)
{
  MPI_Comm            mpicomm = tensor->mesh->ma->mpicomm;
  int                 mpiret;

  const ymir_locidx_t  n_nodes = tensor->dvec->m * tensor->dvec->n
                                 / tensor->ndfields;
  ymir_locidx_t       nodeid;
  double              norm;
  double              min_loc = DBL_MAX;
  double              min_glo;

  YMIR_ASSERT (tensor->ndfields == 6);

  for (nodeid = 0; nodeid < n_nodes; nodeid++) { /* loop over all nodes */
    const double        *tensor_node = tensor->dvec->e[0] + 6 * nodeid;

    /* compute norm for this node */
    norm = sqrt(  tensor_node[0] * tensor_node[0]
                + 2.0 * tensor_node[1] * tensor_node[1]
                + 2.0 * tensor_node[2] * tensor_node[2]
                + tensor_node[3] * tensor_node[3]
                + 2.0 * tensor_node[4] * tensor_node[4]
                + tensor_node[5] * tensor_node[5] );

    /* update processor-local min */
    min_loc = SC_MIN (min_loc, norm);
  }

  /* communicate global min */
  mpiret = MPI_Allreduce (&min_loc, &min_glo, 1, MPI_DOUBLE, MPI_MIN, mpicomm);
  YMIR_CHECK_MPI (mpiret);

  /* return global min Frobenius-norm of symmetric tensor */
  return min_glo;
}

/**
 *
 */
double
slabs_norm_symtensor_frobenius_max (const ymir_dvec_t *tensor)
{
  MPI_Comm            mpicomm = tensor->mesh->ma->mpicomm;
  int                 mpiret;

  const ymir_locidx_t  n_nodes = tensor->dvec->m * tensor->dvec->n
                                 / tensor->ndfields;
  ymir_locidx_t       nodeid;
  double              norm;
  double              max_loc = 0.0;
  double              max_glo;

  YMIR_ASSERT (tensor->ndfields == 6);

  for (nodeid = 0; nodeid < n_nodes; nodeid++) { /* loop over all nodes */
    const double        *tensor_node = tensor->dvec->e[0] + 6 * nodeid;

    /* compute norm for this node */
    norm = sqrt(  tensor_node[0] * tensor_node[0]
                + 2.0 * tensor_node[1] * tensor_node[1]
                + 2.0 * tensor_node[2] * tensor_node[2]
                + tensor_node[3] * tensor_node[3]
                + 2.0 * tensor_node[4] * tensor_node[4]
                + tensor_node[5] * tensor_node[5] );

    /* update processor-local max */
    max_loc = SC_MAX (max_loc, norm);
  }

  /* communicate global max */
  mpiret = MPI_Allreduce (&max_loc, &max_glo, 1, MPI_DOUBLE, MPI_MIN, mpicomm);
  YMIR_CHECK_MPI (mpiret);

  /* return max Frobenius-norm of symmetric tensor */
  return max_glo;
}

/**
 *
 */
void
slabs_norm_symtensor_normalize_threshold_frobenius_elem (sc_dmatrix_t *tensor,
                                                         const double
                                                         threshold)
{
  const ymir_locidx_t  n_nodes = tensor->m;
  ymir_locidx_t       nodeid;
  int                 fieldid;
  double              norm;

  YMIR_ASSERT (tensor->n == 6);

  for (nodeid = 0; nodeid < n_nodes; nodeid++) { /* loop over all nodes */
    double             *tensor_node = tensor->e[0] + 6 * nodeid;

    /* compute norm for this node */
    norm =   tensor_node[0] * tensor_node[0]
           + 2.0 * tensor_node[1] * tensor_node[1]
           + 2.0 * tensor_node[2] * tensor_node[2]
           + tensor_node[3] * tensor_node[3]
           + 2.0 * tensor_node[4] * tensor_node[4]
           + tensor_node[5] * tensor_node[5];
    norm = sqrt (norm);

    /* normalize all fields at this node */
    if (threshold < norm) {
      for (fieldid = 0; fieldid < 6; fieldid++) { /* loop over all fields */
        tensor_node[fieldid] /= norm;
      }
    }
  }
}

/**
 *
 */
void
slabs_norm_symtensor_normalize_frobenius_elem (sc_dmatrix_t *tensor)
{
  slabs_norm_symtensor_normalize_threshold_frobenius_elem (tensor, 0.0);
}

/**
 *
 */
ymir_gloidx_t
slabs_norm_symtensor_normalize_threshold_frobenius (ymir_dvec_t *tensor,
                                                    const double threshold)
{
  MPI_Comm            mpicomm = tensor->mesh->ma->mpicomm;
  int                 mpiret;

  const ymir_locidx_t  n_nodes = tensor->dvec->m * tensor->dvec->n
                                 / tensor->ndfields;
  ymir_locidx_t       nodeid;
  int                 fieldid;
  double              norm;
  ymir_gloidx_t       n_scaled_nodes_loc = 0;
  ymir_gloidx_t       n_scaled_nodes_glo;

  YMIR_ASSERT (tensor->ndfields == 6);

  for (nodeid = 0; nodeid < n_nodes; nodeid++) { /* loop over all nodes */
    double             *tensor_node = tensor->dvec->e[0] + 6 * nodeid;

    /* compute norm for this node */
    norm =   tensor_node[0] * tensor_node[0]
           + 2.0 * tensor_node[1] * tensor_node[1]
           + 2.0 * tensor_node[2] * tensor_node[2]
           + tensor_node[3] * tensor_node[3]
           + 2.0 * tensor_node[4] * tensor_node[4]
           + tensor_node[5] * tensor_node[5];
    norm = sqrt (norm);

    /* normalize all fields at this node */
    if (threshold < norm) {
      for (fieldid = 0; fieldid < 6; fieldid++) { /* loop over all fields */
        tensor_node[fieldid] /= norm;
      }
      n_scaled_nodes_loc++;
    }
  }

  /* get processor-global number of scaled nodes */
  mpiret = MPI_Allreduce (&n_scaled_nodes_loc, &n_scaled_nodes_glo,
                          1, MPI_INT64_T, MPI_SUM, mpicomm);
  YMIR_CHECK_MPI (mpiret);

  /* return number of scaled nodes */
  return n_scaled_nodes_glo;
}

/**
 *
 */
ymir_gloidx_t
slabs_norm_symtensor_normalize_frobenius (ymir_dvec_t *tensor)
{
  return slabs_norm_symtensor_normalize_threshold_frobenius (tensor, 0.0);
}

/**
 *
 */
void
slabs_norm_symtensor_scale_threshold_frobenius_elem (sc_dmatrix_t *tensor,
                                                     sc_dmatrix_t *scaling,
                                                     const double threshold)
{
  const ymir_locidx_t  n_nodes = tensor->m;
  const double       *scaling_data = scaling->e[0];
  ymir_locidx_t       nodeid;
  int                 fieldid;
  double              norm;

  YMIR_ASSERT (tensor->n == 6);
  YMIR_ASSERT (scaling->n == 1);
  YMIR_ASSERT (tensor->m == scaling->m);

  for (nodeid = 0; nodeid < n_nodes; nodeid++) { /* loop over all nodes */
    double             *tensor_node = tensor->e[0] + 6 * nodeid;

    /* compute norm for this node */
    norm =   tensor_node[0] * tensor_node[0]
           + 2.0 * tensor_node[1] * tensor_node[1]
           + 2.0 * tensor_node[2] * tensor_node[2]
           + tensor_node[3] * tensor_node[3]
           + 2.0 * tensor_node[4] * tensor_node[4]
           + tensor_node[5] * tensor_node[5];
    norm = sqrt (norm);

    /* normalize all fields at this node */
    if ((threshold * scaling_data[nodeid]) < norm) {
      for (fieldid = 0; fieldid < 6; fieldid++) { /* loop over all fields */
        tensor_node[fieldid] *= scaling_data[nodeid] / norm;
      }
    }
  }
}

/**
 *
 */
void
slabs_norm_symtensor_scale_frobenius_elem (sc_dmatrix_t *tensor,
                                           sc_dmatrix_t *scaling)
{
  slabs_norm_symtensor_scale_threshold_frobenius_elem (tensor, scaling, 0.0);
}

/**
 *
 */
ymir_gloidx_t
slabs_norm_symtensor_scale_threshold_frobenius (ymir_dvec_t *tensor,
                                                const ymir_dvec_t *scaling,
                                                const double threshold)
{
  MPI_Comm            mpicomm = tensor->mesh->ma->mpicomm;
  int                 mpiret;

  const ymir_locidx_t  n_nodes = tensor->dvec->m * tensor->dvec->n
                                 / tensor->ndfields;
  const double       *scaling_data = scaling->dvec->e[0];
  ymir_locidx_t       nodeid;
  int                 fieldid;
  double              norm;
  ymir_gloidx_t       n_scaled_nodes_loc = 0;
  ymir_gloidx_t       n_scaled_nodes_glo;

  YMIR_ASSERT (tensor->ndfields == 6);
  YMIR_ASSERT (scaling->ndfields == 1);
  YMIR_ASSERT ((scaling->dvec->m * scaling->dvec->n) == n_nodes);

  for (nodeid = 0; nodeid < n_nodes; nodeid++) { /* loop over all nodes */
    double             *tensor_node = tensor->dvec->e[0] + 6 * nodeid;

    /* compute norm for this node */
    norm =   tensor_node[0] * tensor_node[0]
           + 2.0 * tensor_node[1] * tensor_node[1]
           + 2.0 * tensor_node[2] * tensor_node[2]
           + tensor_node[3] * tensor_node[3]
           + 2.0 * tensor_node[4] * tensor_node[4]
           + tensor_node[5] * tensor_node[5];
    norm = sqrt (norm);

    /* scale all fields at this node */
    if ((threshold * scaling_data[nodeid]) < norm) {
      for (fieldid = 0; fieldid < 6; fieldid++) { /* loop over all fields */
        tensor_node[fieldid] *= scaling_data[nodeid] / norm;
      }
      n_scaled_nodes_loc++;
    }
  }

  /* get processor-global number of scaled nodes */
  mpiret = MPI_Allreduce (&n_scaled_nodes_loc, &n_scaled_nodes_glo,
                          1, MPI_INT64_T, MPI_SUM, mpicomm);
  YMIR_CHECK_MPI (mpiret);

  /* return number of scaled nodes */
  return n_scaled_nodes_glo;
}

/**
 *
 */
ymir_gloidx_t
slabs_norm_symtensor_scale_frobenius (ymir_dvec_t *tensor,
                                      const ymir_dvec_t *scale)
{
  return slabs_norm_symtensor_scale_threshold_frobenius (tensor, scale, 0.0);
}

/**
 * MINRES: preconditioner norm
 *   || r ||_{P^-1} = sqrt(r^T * P * r)
 *
 * GMRES with left block preconditioning: preconditioned residual
 *   || P * r ||_{l^2}
 *
 * GMRES with right or symmetric block preconditioning: unprecond. residual
 *   || r ||_{l^2}
 */
double
slabs_norm_of_residual (double *norm_vel, double *norm_press,
                        ymir_vec_t *residual_up,
                        ymir_stokes_op_t *stokes_op,
                        ymir_stokes_pc_t *stokes_pc,
                        slabs_krylov_type_t krylov_type,
                        slabs_norm_type_t norm_type,
                        ymir_Hminus1_norm_op_t *norm_op)
{
  ymir_pressure_elem_t  *press_elem = stokes_op->press_elem;

  ymir_vec_t         *up_l, *up_r;
  double              norm_res, nrm_vel, nrm_press;

  /* apply preconditioner */
  switch (krylov_type) {
  case SL_KRYLOV_MINRES:
    YMIR_ASSERT (stokes_pc != NULL);
    up_r = ymir_vec_template (residual_up);
    YMIR_ABORT_NOT_REACHED (); //TODO following line calls wrong function
    ymir_nlstokes_pc_apply (residual_up, up_r, stokes_pc);
    up_l = residual_up;
    break;

  case SL_KRYLOV_GMRES:
    if (ymir_stokes_pc_block_type == YMIR_STOKES_PC_BLOCK_L) {
      /* if left preconditioning */
      YMIR_ASSERT (stokes_pc != NULL);
      YMIR_ASSERT (norm_type != SL_NORM_FNC_HMINUS1_L2);
      up_r = ymir_vec_template (residual_up);
      if (ymir_stress_op_is_nl (stokes_pc->stokes_op->stress_op)) {
        ymir_nlstokes_pc_apply_lower (residual_up, up_r, stokes_pc);
      }
      else {
        ymir_stokes_pc_apply_lower (residual_up, up_r, stokes_pc);
      }
      up_l = up_r;
    }
    else { /* if right or symmetric preconditioning */
      up_r = residual_up;
      up_l = residual_up;
    }
    break;

  default: /* unknown Krylov method */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* compute norm of residual */
  slabs_norm_innerprod_comp (&nrm_vel, &nrm_press, up_l, up_r, press_elem,
                             norm_type, norm_op);
  norm_res = sqrt (nrm_vel + nrm_press);
  if (norm_vel != NULL) {
    *norm_vel = sqrt (nrm_vel);
  }
  if (norm_press != NULL) {
    *norm_press = sqrt (nrm_press);
  }

  /* destroy */
  switch (krylov_type) {
  case SL_KRYLOV_MINRES:
    ymir_vec_destroy (up_r);
    break;

  case SL_KRYLOV_GMRES:
    if (ymir_stokes_pc_block_type == YMIR_STOKES_PC_BLOCK_L) {
      ymir_vec_destroy (up_r);
    }
    break;

  default: /* unknown Krylov method */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* return norm of residual */
  return norm_res;
}

/**
 *
 */
double
slabs_norm_of_residual_l2 (ymir_vec_t *residual_up,
                           ymir_stokes_op_t *stokes_op,
                           ymir_stokes_pc_t *stokes_pc,
                           slabs_krylov_type_t krylov_type)
{
  return slabs_norm_of_residual (NULL, NULL, residual_up, stokes_op, stokes_pc,
                                 krylov_type, SL_NORM_VEC_L2, NULL);
}

/**
 *
 */
double
slabs_norm_compute_residual (ymir_vec_t *residual_up,
                             double *norm_vel, double *norm_press,
                             ymir_vec_t *up,
                             ymir_cvec_t *rhs_u_point,
                             ymir_stokes_op_t *stokes_op,
                             ymir_stokes_pc_t *stokes_pc,
                             slabs_krylov_type_t krylov_type,
                             slabs_norm_type_t norm_type,
                             ymir_Hminus1_norm_op_t *norm_op)
{
  ymir_mesh_t        *mesh = rhs_u_point->mesh;
  ymir_pressure_elem_t  *press_elem = stokes_op->press_elem;

  ymir_vec_t         *rhs;
  double              norm_res;

  /* construct the right-hand side */
  rhs = ymir_stokes_vec_new (mesh, press_elem);
  ymir_stokes_pc_construct_rhs (rhs, rhs_u_point, NULL, NULL,
                                1 /* incompressible */, stokes_op, 0);
  YMIR_ASSERT (sc_dmatrix_is_valid (rhs->dataown));
  YMIR_ASSERT (sc_dmatrix_is_valid (rhs->coff));
  YMIR_ASSERT (ymir_vec_is_not_dirty (rhs));

  /* determine the residual */
  if (up != NULL) {
    YMIR_ASSERT (residual_up != NULL);
    YMIR_ASSERT (sc_dmatrix_is_valid (up->dataown));
    YMIR_ASSERT (sc_dmatrix_is_valid (up->coff));
    YMIR_ASSERT (ymir_vec_is_not_dirty (up));

    /* compute residual
     *   r = b - K * x
     * where
     *   r --- residual
     *   b --- right-hand side
     *   K --- Stokes operator
     *   x --- velocity-pressure vector */
    ymir_stokes_pc_apply_stokes_op (up, residual_up, stokes_op, 0, 0);
    ymir_vec_add (-1.0, rhs, residual_up);
    ymir_vec_scale (-1.0, residual_up);
  }
  else {
    YMIR_ASSERT (residual_up == NULL);

    /* set residual to be the right-hand side */
    residual_up = rhs;
  }
  YMIR_ASSERT (sc_dmatrix_is_valid (residual_up->dataown));
  YMIR_ASSERT (sc_dmatrix_is_valid (residual_up->coff));
  YMIR_ASSERT (slabs_vec_is_finite (residual_up));

  /* compute norm of residual */
  norm_res = slabs_norm_of_residual (norm_vel, norm_press, residual_up,
                                     stokes_op, stokes_pc, krylov_type,
                                     norm_type, norm_op);

  /* destroy vectors */
  ymir_vec_destroy (rhs);

  /* return norm of residual */
  return norm_res;
}

/**
 *
 */
double
slabs_norm_compute_residual_l2 (ymir_vec_t *residual_up,
                                ymir_vec_t *up,
                                ymir_cvec_t *rhs_u_point,
                                ymir_stokes_op_t *stokes_op,
                                ymir_stokes_pc_t *stokes_pc,
                                slabs_krylov_type_t krylov_type)
{
  return slabs_norm_compute_residual (residual_up, NULL, NULL, up, rhs_u_point,
                                      stokes_op, stokes_pc, krylov_type,
                                      SL_NORM_VEC_L2, NULL);
}

/**
 * Computes the residual for the dual equation.
 */
double
slabs_norm_compute_residual_dual (ymir_dvec_t *residual_dual,
                                  ymir_vec_t *up,
                                  ymir_dvec_t *dual,
                                  ymir_stokes_op_t *stokes_op,
                                  const slabs_nl_solver_primaldual_type_t
                                    nl_solver_primaldual_type,
                                  const slabs_krylov_type_t krylov_type)
{
  double              norm_res;

  switch (nl_solver_primaldual_type) {
  case SL_NL_SOLVER_PRIMALDUAL_NONE:
    norm_res = 0.0;
    break;

  case SL_NL_SOLVER_PRIMALDUAL_NORMSTRAIN:
    /* compute residual:
     *
     *   residual_dual = grad_sym (u) - || grad_sym (u) ||_{F} * dual
     */
    {
      ymir_mesh_t        *mesh = up->mesh;
      mangll_t           *mangll = mesh->ma;
      const mangll_locidx_t  n_elements = mangll->mesh->K;
      const int           N = ymir_n (mangll->N);
      const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);
      ymir_cvec_t        *u;
      sc_dmatrix_t       *vel_el_mat, *dual_el_mat;
      sc_dmatrix_t       *grad_vel_el_mat, *tmp_dvel, *tmp;
      sc_dmatrix_t       *strain_frob_el_mat;
      sc_dmatrix_t       *residual_dual_el_mat;
      mangll_locidx_t     elid;

      /* create views */
      slabs_stokes_vec_get_components_view (&u, NULL, up);

      /* create work variables */
      vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
      dual_el_mat = sc_dmatrix_new (n_nodes_per_el, 6);
      grad_vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 9);
      tmp_dvel = sc_dmatrix_new (n_nodes_per_el, 3);
      tmp = sc_dmatrix_new (n_nodes_per_el, 3);
      strain_frob_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
      residual_dual_el_mat = sc_dmatrix_new (n_nodes_per_el, 6);

      for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
        /* get velocity field of this element at GLL nodes */
        ymir_cvec_get_elem_interp (u, vel_el_mat, YMIR_STRIDE_NODE, elid,
                                   YMIR_GLL_NODE, YMIR_READ);

        /* get dual tensor for this element */
        ymir_dvec_get_elem (dual, dual_el_mat, YMIR_STRIDE_NODE, elid,
                            YMIR_COPY);

        /* compute 2nd invariant of the strain rate at Gauss nodes
         * (store in `strain_frob_el_mat`) */
        slabs_second_invariant_elem (vel_el_mat, strain_frob_el_mat, mangll,
                                     elid, grad_vel_el_mat, tmp_dvel, tmp,
                                     SL_GAUSS_NODE);

        /* transform 2nd invarian to Frobenius-norm of the strain rate tensor */
        sc_dmatrix_scale (2.0, strain_frob_el_mat);
        sc_dmatrix_sqrt (strain_frob_el_mat, strain_frob_el_mat);

        /* multiply in Frobenius-norm to dual */
        slabs_matrix_multiply_in_1d (strain_frob_el_mat, dual_el_mat);

        /* compute strain rate tensor for this element
         * (store in `residual_dual_el_mat`) */
        slabs_strain_rate_tensor_elem (residual_dual_el_mat, mangll,
                                       grad_vel_el_mat);

        /* compute difference */
        sc_dmatrix_add (-1.0, dual_el_mat, residual_dual_el_mat);

        /* set values of residual vector */
        ymir_dvec_set_elem (residual_dual, residual_dual_el_mat,
                            YMIR_STRIDE_NODE, elid, YMIR_SET);
      }

      /* compute L^2-norm of residual */
      slabs_norm_vec_L2_innerprod_comp (NULL, &norm_res, NULL,
                                        residual_dual, residual_dual);
      norm_res = sqrt (norm_res);

      /* destroy work variables */
      sc_dmatrix_destroy (vel_el_mat);
      sc_dmatrix_destroy (dual_el_mat);
      sc_dmatrix_destroy (grad_vel_el_mat);
      sc_dmatrix_destroy (tmp_dvel);
      sc_dmatrix_destroy (tmp);
      sc_dmatrix_destroy (strain_frob_el_mat);
      sc_dmatrix_destroy (residual_dual_el_mat);
      ymir_vec_destroy (u);
    }
    break;

  case SL_NL_SOLVER_PRIMALDUAL_STRESS:
    /* compute residual:
     *
     *   residual_dual = grad_sym (u) - 1 / (2 * viscosity) * dual
     */
    {
      ymir_dvec_t        *viscosity = stokes_op->stress_op->viscosity;
      ymir_cvec_t        *u;
      ymir_dvec_t        *tmp_6d = ymir_vec_template (dual);

      /* create view */
      slabs_stokes_vec_get_components_view (&u, NULL, up);

      /* compute the strain rate of the velocity iterate */
      ymir_velocity_strain_rate (u, residual_dual, 0);

      /* get the dual iterate */
      ymir_dvec_copy (dual, tmp_6d);

      /* divide in viscosity */
      ymir_dvec_divide_in1 (viscosity, tmp_6d);

      /* set difference */
      ymir_dvec_add (-1.0, tmp_6d, residual_dual);

      /* compute L^2-norm of residual */
      slabs_norm_vec_L2_innerprod_comp (NULL, &norm_res, NULL,
                                        residual_dual, residual_dual);
      norm_res = sqrt (norm_res);

      /* destroy */
      ymir_vec_destroy (tmp_6d);
      ymir_vec_destroy (u);
    }
    break;

  default: /* unknown primal-dual type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* return norm of dual residual */
  return norm_res;
}

/**
 *
 */
void
slabs_norm_weight_residual (ymir_vec_t *residual_up,
                            ymir_pressure_elem_t *press_elem)
{
  ymir_mesh_t        *mesh = residual_up->mesh;
  ymir_vec_t         *res_u, *res_p;
  ymir_vec_t         *vel_lump = ymir_cvec_new (mesh, 3);
  ymir_vec_t         *press_lump = ymir_pressure_vec_new (mesh, press_elem);

  /* get vector components */
  slabs_stokes_vec_get_components (&res_u, &res_p, residual_up, press_elem);

  /* create lumped mass matrices */
  ymir_mass_lump (vel_lump);
  ymir_pressure_vec_lump_mass (press_lump, press_elem);

  /* scale by inverse of lumped mass matrix */
  ymir_vec_divide_in (vel_lump, res_u);
  ymir_vec_divide_in (press_lump, res_p);
  slabs_stokes_vec_set_components (residual_up, res_u, res_p, press_elem);

  /* destroy */
  ymir_vec_destroy (vel_lump);
  ymir_vec_destroy (press_lump);
  ymir_vec_destroy (res_u);
  ymir_vec_destroy (res_p);
}

/******************************************************************************
 * Tests
 *****************************************************************************/

/* parameter list for test of Frobenius-norm */
typedef struct slabs_norm_test_frobenius_set_data
{
  double              cx, cy, cz;
}
slabs_norm_test_frobenius_set_data_t;

/**
 *
 */
static void
slabs_norm_test_frobenius_set_vel_fn (double *val, double x, double y, double z,
                                      ymir_locidx_t nid, void *data)
{
  slabs_norm_test_frobenius_set_data_t  *d =
    (slabs_norm_test_frobenius_set_data_t *) data;
  const double        cx = d->cx;
  const double        cy = d->cy;
  const double        cz = d->cz;

  val[0] = exp (M_PI * x);
  val[1] = cos (cx * M_PI * x) * cos (cy * M_PI * y) * cos (cz * M_PI * z);
  val[2] = sin (cx * M_PI * x) * sin (cy * M_PI * y) * sin (cz * M_PI * z);
}

/**
 * Tests computation of Frobenius-norm of the strain rate tensor.
 */
void
slabs_norm_test_frobenius_of_strain_rate (ymir_mesh_t *mesh)
{
  const char         *this_fn_name = "slabs_norm_test_frobenius_of_strain_rate";
  ymir_cvec_t        *vel_vec;
  ymir_dvec_t        *frob_IIe_direct = ymir_dvec_new (mesh, 1,
                                                       YMIR_GAUSS_NODE);
  ymir_dvec_t        *frob_IIe_elem = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
  ymir_dvec_t        *frob_strain_elem = ymir_dvec_new (mesh, 1,
                                                        YMIR_GAUSS_NODE);
  slabs_norm_test_frobenius_set_data_t  data;
  double              norm_frob;
  double              abs_error_IIe_direct;
  double              abs_error_IIe_elem;

  /* create velocity vector */
  data.cx = 2.0;
  data.cy = 4.0;
  data.cz = 1.0;
  vel_vec = ymir_cvec_new (mesh, 3);
  ymir_cvec_set_function (vel_vec, slabs_norm_test_frobenius_set_vel_fn, &data);

  /* compute Frobenius-norm via 2nd invariant */
  {
    ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (
                                          mesh->ma->N, mesh->ma->ompsize);

    ymir_second_invariant_vec (vel_vec, frob_IIe_direct, vel_elem);
    ymir_dvec_scale (2.0, frob_IIe_direct);
    ymir_dvec_sqrt (frob_IIe_direct, frob_IIe_direct);
    ymir_velocity_elem_destroy (vel_elem);
  }

  /* compute Frobenius-norm via 2nd invariant and from strain rate tensor
   * element wise */
  {
    mangll_t           *mangll = mesh->ma;
    const mangll_locidx_t  n_elements = mangll->mesh->K;
    const int           N = ymir_n (mangll->N);
    const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);
    sc_dmatrix_t       *vel_el_mat;
    sc_dmatrix_t       *grad_vel_el_mat, *tmp_dvel, *tmp_3d;
    sc_dmatrix_t       *IIe_el_mat;
    sc_dmatrix_t       *strain_vel_el_mat;
    sc_dmatrix_t       *frob_IIe_el_mat;
    sc_dmatrix_t       *frob_strain_el_mat;
    mangll_locidx_t     elid;

    /* create work variables */
    vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
    grad_vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 9);
    tmp_dvel = sc_dmatrix_new (n_nodes_per_el, 3);
    tmp_3d = sc_dmatrix_new (n_nodes_per_el, 3);
    IIe_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
    strain_vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 6);
    frob_IIe_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
    frob_strain_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);

    for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
      /* get velocity of this element at GLL nodes */
      ymir_cvec_get_elem_interp (vel_vec, vel_el_mat, YMIR_STRIDE_NODE, elid,
                                 YMIR_GLL_NODE, YMIR_READ);

      /* compute 2nd invariant for this element */
      slabs_second_invariant_elem (vel_el_mat, IIe_el_mat, mangll, elid,
                                   grad_vel_el_mat, tmp_dvel, tmp_3d,
                                   SL_GAUSS_NODE);

      /* compute strain rate tensor for this element */
      slabs_strain_rate_tensor_elem (strain_vel_el_mat, mangll,
                                     grad_vel_el_mat);

      /* compute Frobenius-norm via 2nd invariant */
      sc_dmatrix_copy (IIe_el_mat, frob_IIe_el_mat);
      sc_dmatrix_scale (2.0, frob_IIe_el_mat);
      sc_dmatrix_sqrt (frob_IIe_el_mat, frob_IIe_el_mat);

      /* compute Frobenius-norm of strain rate tensor */
      slabs_norm_symtensor_frobenius_elem (strain_vel_el_mat,
                                           frob_strain_el_mat);

      /* set values of Frobenius-norm */
      ymir_dvec_set_elem (frob_IIe_elem, frob_IIe_el_mat,
                          YMIR_STRIDE_NODE, elid, YMIR_SET);
      ymir_dvec_set_elem (frob_strain_elem, frob_strain_el_mat,
                          YMIR_STRIDE_NODE, elid, YMIR_SET);
    }

    /* destroy work variables */
    sc_dmatrix_destroy (vel_el_mat);
    sc_dmatrix_destroy (grad_vel_el_mat);
    sc_dmatrix_destroy (tmp_dvel);
    sc_dmatrix_destroy (tmp_3d);
    sc_dmatrix_destroy (IIe_el_mat);
    sc_dmatrix_destroy (strain_vel_el_mat);
    sc_dmatrix_destroy (frob_IIe_el_mat);
    sc_dmatrix_destroy (frob_strain_el_mat);
  }

  /* check error */
  ymir_dvec_set_zero (frob_strain_elem);
  norm_frob = ymir_dvec_norm (frob_strain_elem);
  ymir_vec_add (-1.0, frob_strain_elem, frob_IIe_direct);
  ymir_vec_add (-1.0, frob_strain_elem, frob_IIe_elem);
  abs_error_IIe_direct = ymir_dvec_norm (frob_IIe_direct);
  abs_error_IIe_elem = ymir_dvec_norm (frob_IIe_elem);

  /* print error */
  if (0.0 < norm_frob) {
    YMIR_GLOBAL_VERBOSEF ("%s: IIe direct: abs error %1.3e, rel error %1.3e\n",
                          this_fn_name, abs_error_IIe_direct,
                          abs_error_IIe_direct / norm_frob);
    YMIR_GLOBAL_VERBOSEF ("%s: IIe elem: abs error %1.3e, rel error %1.3e\n",
                          this_fn_name, abs_error_IIe_elem,
                          abs_error_IIe_elem / norm_frob);
  }
  else {
    YMIR_GLOBAL_VERBOSEF ("%s: IIe direct: abs error %1.3e\n",
                          this_fn_name, abs_error_IIe_direct);
    YMIR_GLOBAL_VERBOSEF ("%s: IIe elem: abs error %1.3e\n",
                          this_fn_name, abs_error_IIe_elem);
  }

  /* destroy */
  ymir_vec_destroy (vel_vec);
  ymir_vec_destroy (frob_IIe_direct);
  ymir_vec_destroy (frob_IIe_elem);
  ymir_vec_destroy (frob_strain_elem);
}

/* parameter list for test of H^-1 norm operator */
typedef struct slabs_norm_Hminus1_stiff_test_sol_data
{
  double              cx, cy, cz;
}
slabs_norm_Hminus1_stiff_test_sol_data_t;

/**
 *
 */
static void
slabs_norm_Hminus1_stiff_test_sol_fn (double *sol, double x, double y,
                                      double z, ymir_locidx_t nid, void *data)
{
  slabs_norm_Hminus1_stiff_test_sol_data_t
    *d = (slabs_norm_Hminus1_stiff_test_sol_data_t *) data;

  *sol = sin (M_PI * d->cx * x) * sin (M_PI * d->cy * y)
       * sin (M_PI * d->cz * z);
}

static void
slabs_norm_Hminus1_stiff_mass_test_sol_fn (double *sol, double x, double y,
                                           double z, ymir_locidx_t nid,
                                           void *data)
{
  slabs_norm_Hminus1_stiff_test_sol_data_t
    *d = (slabs_norm_Hminus1_stiff_test_sol_data_t *) data;

  *sol = cos (M_PI * d->cx * x) * cos (M_PI * d->cy * y)
       * cos (M_PI * d->cz * z);
}

/**
 * Tests the Operator of the H^-1 norm by solving a Poisson problem and
 * comparing the error to the exact solution.
 */
void
slabs_norm_Hminus1_stiff_test (ymir_Hminus1_norm_op_t *norm_op,
                               ymir_mesh_t *mesh,
                               slabs_domain_shape_t domain_shape)
{
  const char         *this_fn_name = "slabs_norm_Hminus1_stiff_test";
  const double        mass_scaling = norm_op->mass_scaling;
  ymir_vec_t         *sol, *sol_exact;
  ymir_vec_t         *rhs, *rhs_point;
  double              cx, cy, cz;
  slabs_norm_Hminus1_stiff_test_sol_data_t  data;
  double              rel_error;

  /* test is only implemented for the cube domain */
  YMIR_ASSERT (domain_shape == SL_DOMAIN_CUBE);

  /* initialize vectors */
  sol = ymir_cvec_new (mesh, 1);
  sol_exact = ymir_cvec_new (mesh, 1);
  rhs = ymir_cvec_new (mesh, 1);
  rhs_point = ymir_cvec_new (mesh, 1);

  /* generate random numbers */
  cx = 1.0; //round ( ((double) rand () / (double) RAND_MAX) * 2.0 );
  cy = 2.0; //round ( ((double) rand () / (double) RAND_MAX) * 3.0 );
  cz = 3.0; //round ( ((double) rand () / (double) RAND_MAX) * 4.0 );
  data.cx = cx;
  data.cy = cy;
  data.cz = cz;

  /* set exact solution */
  if (0 < mass_scaling) {
    /* u(x) = sin(PI * cx * x) + sin(PI * cy * y) + sin(PI * cz * z) */
    ymir_cvec_set_function (
        sol_exact, slabs_norm_Hminus1_stiff_mass_test_sol_fn, &data);
  }
  else {
    /* u(x) = cos(PI * cx * x) + cos(PI * cy * y) + cos(PI * cz * z) */
    ymir_cvec_set_function (
        sol_exact, slabs_norm_Hminus1_stiff_test_sol_fn, &data);
  }

  /* set point values of right-hand side:
   *
   *   f(x) = (PI^2 * (cx^2 + cy^2 + cz^2) + mass_scaling) * u(x)
   */
  ymir_cvec_copy (sol_exact, rhs_point);
  ymir_cvec_scale (M_PI * M_PI * (cx * cx + cy * cy + cz * cz) + mass_scaling,
                   rhs_point);

  /* compute right-hand side for algebraic system */
  ymir_mass_apply (rhs_point, rhs);

#ifdef YMIR_PETSC
  {
    KSP                 ksp;

    /* get KSP */
    if (0 < mass_scaling) {
      ksp = norm_op->ksp;
    }
    else {
      ksp = norm_op->stiff_pc->ksp;
    }

    /* solve Poisson problem */
    ymir_petsc_ksp_solve (sol, rhs, ksp, 0 /* zero initial guess */,
                          YMIR_MESH_PETSCLAYOUT_NONE, NULL);
  }
#endif

  /* compute relative error */
  ymir_cvec_add (-1.0, sol_exact, sol);
  rel_error = ymir_cvec_norm (sol) / ymir_cvec_norm (sol_exact);

  /* output of results */
  YMIR_GLOBAL_VERBOSEF ("%s: (cx,cy,cz) = (%g,%g,%g), mass_scaling %g, "
                        "rel solution error %1.3e\n",
                        this_fn_name, cx, cy, cz, mass_scaling, rel_error);

  /* destroy */
  ymir_vec_destroy (sol);
  ymir_vec_destroy (sol_exact);
  ymir_vec_destroy (rhs);
  ymir_vec_destroy (rhs_point);
}

/**
 *
 */
static int
slabs_norm_test_refine_all_fn (p8est_t *p8est, p4est_topidx_t tree,
                               p8est_quadrant_t *quadrant)
{
  return 1;
}

#if 0
static int
slabs_norm_test_refine_corner_fn (p8est_t *p8est, p4est_topidx_t tree,
                                  p8est_quadrant_t *quadrant)
{
  double              coord[P8EST_DIM];

  /* get coordinates of quadrant */
  p8est_qcoord_to_vertex (p8est->connectivity, tree, quadrant->x,
                          quadrant->y, quadrant->z, coord);

  /* refine in corner */
  if ( (0.4 < coord[0] && 0.4 < coord[1] && 0.4 < coord[2]) ||
       (0.4 > coord[0] && 0.4 > coord[1] && 0.4 > coord[2]) ) {
    return 1;
  }
  else {
    return 0;
  }
}
#endif

/* parameter list for testing mesh-independence */
typedef struct slabs_norm_test_mesh_independence_data
{
  double              cx, cy, cz;
}
slabs_norm_test_mesh_independence_data_t;

/**
 *
 */
static void
slabs_norm_test_mesh_independence_fn (double *val, double x, double y, double z,
                                      ymir_locidx_t nid, void *data)
{
  slabs_norm_test_mesh_independence_data_t
    *d = (slabs_norm_test_mesh_independence_data_t *) data;

  if (x < 0.5) {
    *val = 0.0;
  }
  else {
    *val = exp (d->cx * x) * sin (M_PI * d->cy * y) * sin (M_PI * d->cz * z);
  }
}

/**
 * Tests norms for mesh-independence.
 */
void
slabs_norm_test_mesh_independence (MPI_Comm mpicomm,
                                   const slabs_domain_shape_t domain_shape,
                                   const int minlevel, const int maxlevel,
                                   char *refine, const int order,
                                   double mass_scaling,
                                   const int max_refine_steps)
{
  const char         *this_fn_name = "slabs_norm_test_mesh_independence";

  slabs_physics_options_t  physics_options;
  slabs_discr_options_t  discr_options;
  p8est_t            *p8est;
  mangll_t           *mangll;
  mangll_cnodes_t    *cnodes;
  ymir_mesh_t        *mesh;
  ymir_pressure_elem_t  *press_elem;
  ymir_Hminus1_norm_op_t  *norm_op;
  ymir_cvec_t        *sol_exact, *rhs, *rhs_point;
  ymir_cvec_t        *guess, *guess_point;
  ymir_cvec_t        *residual;
  double              cx, cy, cz;
  slabs_norm_Hminus1_stiff_test_sol_data_t  data_op;
  slabs_norm_test_mesh_independence_data_t  data;
  int                 refine_step;
  sc_dmatrix_t       *results;

  /* test is only implemented for the cube domain */
  YMIR_ASSERT (domain_shape == SL_DOMAIN_CUBE);

  /* set required options */
  physics_options.domain_shape = domain_shape;
  discr_options.minlevel = minlevel;
  discr_options.maxlevel = maxlevel;
  discr_options.refine = refine;
  discr_options.X_fn = slabs_discr_identity_X;
  discr_options.order = order;

  /* create p4est and boundary variables */
  p8est = slabs_discr_p8est_new (mpicomm, &physics_options, &discr_options);
  slabs_discr_options_set_boundary (&discr_options, p8est, &physics_options);

  /* create vector to store results */
  results = sc_dmatrix_new (max_refine_steps, 3);

  for (refine_step = 0; refine_step < max_refine_steps; refine_step++) {
    /* refine p4est */
    // 1) refine all
    p8est_refine (p8est, 0, slabs_norm_test_refine_all_fn, NULL);
    // 2) refine in corner
    //p8est_refine (p8est, 0, slabs_norm_test_refine_corner_fn, NULL);
    //p8est_balance (p8est, P8EST_CONNECT_FULL, NULL);

    /* create mesh */
    slabs_discr_mangll_and_cnodes_new (&mangll, &cnodes, p8est, &discr_options);
    slabs_discr_ymir_new (&mesh, &press_elem, mangll, cnodes, &discr_options);

    /* create H^-1 norm operator */
    norm_op = ymir_Hminus1_norm_op_new (mesh, mass_scaling);
    norm_op->mass_scaling = 0.0;

    /* initialize vectors */
    sol_exact = ymir_cvec_new (mesh, 1);
    rhs = ymir_cvec_new (mesh, 1);
    rhs_point = ymir_cvec_new (mesh, 1);
    guess = ymir_cvec_new (mesh, 1);
    guess_point = ymir_cvec_new (mesh, 1);
    residual = ymir_cvec_new (mesh, 1);

    /* generate random numbers */
    cx = 2.0; //round ( ((double) rand () / (double) RAND_MAX) * 2.0 );
    cy = 3.0; //round ( ((double) rand () / (double) RAND_MAX) * 3.0 );
    cz = 4.0; //round ( ((double) rand () / (double) RAND_MAX) * 4.0 );
    data_op.cx = cx;
    data_op.cy = cy;
    data_op.cz = cz;

    /* set exact solution:
     *
     *   u(x) = cos(PI * cx * x) + cos(PI * cy * y) + cos(PI * cz * z)
     */
    ymir_cvec_set_function (
        sol_exact, slabs_norm_Hminus1_stiff_mass_test_sol_fn, &data_op);

    /* set point values of right-hand side:
     *
     *   f(x) = (PI^2 * (cx^2 + cy^2 + cz^2) + mass_scaling) * u(x)
     */
    ymir_cvec_copy (sol_exact, rhs_point);
    ymir_cvec_scale (
        M_PI*M_PI * (cx*cx + cy*cy + cz*cz) + norm_op->mass_scaling, rhs_point);

    /* compute right-hand side for algebraic system */
    ymir_mass_apply (rhs_point, rhs);

    /* set guess */
    data.cx = cx;
    data.cy = cy;
    data.cz = cz;
    ymir_cvec_set_function (guess_point, slabs_norm_test_mesh_independence_fn,
                            &data);

    /* compute residual */
    ymir_cvec_copy (rhs, residual);
    ymir_Hminus1_norm_stiff_mass_apply (guess_point, guess,//TODO outdated
                                        &(norm_op->mass_scaling));
    ymir_cvec_add (-1.0, guess, residual);
    norm_op->mass_scaling = mass_scaling;

    /* [0] compute l^2 norm */
    results->e[refine_step][0] = slabs_norm (NULL, NULL, residual, NULL,
                                             SL_NORM_VEC_L2, NULL);

    /* [1] compute L^2 norm */
    results->e[refine_step][1] = slabs_norm (NULL, NULL, residual, NULL,
                                             SL_NORM_FNC_L2, NULL);

    /* [2] compute H^-1 norm */
    results->e[refine_step][2] = slabs_norm (NULL, NULL, residual, NULL,
                                             SL_NORM_FNC_HMINUS1, norm_op);

    /* destroy */
    ymir_vec_destroy (sol_exact);
    ymir_vec_destroy (rhs);
    ymir_vec_destroy (rhs_point);
    ymir_vec_destroy (guess);
    ymir_vec_destroy (guess_point);
    ymir_vec_destroy (residual);
    ymir_Hminus1_norm_op_destroy (norm_op);
    slabs_discr_ymir_mangll_destroy (mesh, press_elem);
  }

  /* destroy p4est and options */
  slabs_discr_p8est_destroy (p8est);
  slabs_discr_options_clear_boundary (&discr_options);

  /* output of results */
  for (refine_step = 0; refine_step < max_refine_steps; refine_step++) {
    YMIR_GLOBAL_VERBOSEF ("%s: level %i, "
                          "l^2 %1.3e, L^2 %1.3e, H^-1 comp %1.3e\n",
                          this_fn_name, minlevel + refine_step + 1,
                          results->e[refine_step][0],
                          results->e[refine_step][1],
                          results->e[refine_step][2]);
  }

  /* destroy results */
  sc_dmatrix_destroy (results);
}

