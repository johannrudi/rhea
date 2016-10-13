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
#include <mangll_tensor.h>
#include <ymir_velocity_vec.h>
#include <ymir_pressure_vec.h>
#include <ymir_stokes_vec.h>
#ifdef YMIR_PETSC
#include <petscsys.h>
#endif

/* initialize processor-local memory usage */
unsigned long slabs_initial_local_memory_usage = 0;

ymir_topidx_t *
slabs_cube_tree_to_bf (p8est_connectivity_t *conn)
{
  ymir_topidx_t      *tree_to_bf;
  ymir_topidx_t       ntrees = conn->num_trees;
  ymir_topidx_t       ti;
  int                 fi;

  tree_to_bf = YMIR_ALLOC (ymir_topidx_t, 6 * ntrees);

  for (ti = 0; ti < ntrees; ti++) { /* loop over all trees */
    /* set sides */
    for (fi = 0; fi < 4; fi++) { /* loop over all side faces */
      if (conn->tree_to_tree[6 * ti + fi] == ti) {
        switch (fi) {
        case 0:
          tree_to_bf[6 * ti + fi] = SL_SIDE1;
          break;
        case 1:
          tree_to_bf[6 * ti + fi] = SL_SIDE2;
          break;
        case 2:
          tree_to_bf[6 * ti + fi] = SL_SIDE3;
          break;
        case 3:
          tree_to_bf[6 * ti + fi] = SL_SIDE4;
          break;
        default: /* unknown face type */
          YMIR_ABORT_NOT_REACHED ();
        }
      }
      else {
        tree_to_bf[6 * ti + fi] = SL_NONE;
      }
    }

    /* set top and bottom */
    tree_to_bf[6 * ti + 4] = SL_BASE;
    tree_to_bf[6 * ti + 5] = SL_TOP;
  }

  return tree_to_bf;
}

ymir_topidx_t *
slabs_shell_tree_to_bf (p8est_connectivity_t *conn)
{
  ymir_topidx_t      *tree_to_bf;
  ymir_topidx_t       ntrees = conn->num_trees;
  ymir_topidx_t       ti;
  int                 fi;

  tree_to_bf = YMIR_ALLOC (ymir_topidx_t, 6 * ntrees);

  for (ti = 0; ti < ntrees; ti++) { /* loop over all trees */
    /* set sides */
    for (fi = 0; fi < 4; fi++) { /* loop over all side faces */
      tree_to_bf[6 * ti + fi] = SL_NONE;
    }

    /* set top and bottom */
    tree_to_bf[6 * ti + 4] = SL_BASE;
    tree_to_bf[6 * ti + 5] = SL_TOP;
  }

  return tree_to_bf;
}

/**
 *
 */
ymir_dir_code_t
slabs_cube_vel_dir_fn (double X, double Y, double Z,
                       double nx, double ny, double nz,
                       ymir_topidx_t face, ymir_locidx_t node_id,
                       void *dir_data)
{
  int                 bc_type = *((int *) dir_data);

  switch (bc_type) {
  case SL_VEL_BC_DIRICHLET_ALL:
    /* set Dirichlet in all directions for all points */
    return YMIR_VEL_DIRICHLET_ALL;
    break;

  case SL_VEL_BC_DIRICHLET_NORM:
    /* set Dirichlet in normal direction for all points */
    return YMIR_VEL_DIRICHLET_NORM;
    break;

  case SL_VEL_BC_DIR_SIDES_NEU_TB:
    if ((face != SL_BASE) && (face != SL_TOP)) {
      /* set Dirichlet in all directions on side faces */
      return YMIR_VEL_DIRICHLET_ALL;
    }
    else {
      /* bottom and top faces will get Neumann BC's */
      return YMIR_VEL_DIRICHLET_NONE;
    }
    break;

  case SL_VEL_BC_DIR_NORM_SIDES_DIR_ALL_TB:
    if ((face != SL_BASE) && (face != SL_TOP)) {
      /* set Dirichlet in normal direction on side faces */
      return YMIR_VEL_DIRICHLET_NORM;
    }
    else {
      /* set Dirichlet in all directions on bottom and top faces */
      return YMIR_VEL_DIRICHLET_ALL;
    }
    break;

  case SL_VEL_BC_NEUMANN_ALL:
    /* set Neumann for all points */
    return YMIR_VEL_DIRICHLET_NONE;
    break;

  default: /* unknown boundary condition */
    YMIR_ABORT_NOT_REACHED ();
  }
}

/**
 *
 */
ymir_dir_code_t
slabs_shell_vel_dir_fn (double X, double Y, double Z,
                        double nx, double ny, double nz,
                        ymir_topidx_t face, ymir_locidx_t node_id,
                        void *dir_data)
{
  int                 bc_type = *((int *) dir_data);

  switch (bc_type) {
  case SL_VEL_BC_DIRICHLET_ALL:
    /* set Dirichlet BC in all directions for all points */
    return YMIR_VEL_DIRICHLET_ALL;
    break;

  case SL_VEL_BC_DIRICHLET_NORM:
  case SL_VEL_BC_DIRICHLET_NORM_FIXDOF:
    /* set Dirichlet BC in normal direction for all points */
    return YMIR_VEL_DIRICHLET_NORM;
    break;

  case SL_VEL_BC_DIRICHLET_ALL_INNER:
    if (face == SL_TOP) {
      /* set Dirichlet BC in normal direction on outer sphere of a shell */
      return YMIR_VEL_DIRICHLET_NORM;
    }
    else {
      /* set Dirichlet BC in all direction on inner sphere of a shell */
      return YMIR_VEL_DIRICHLET_ALL;
    }
    break;

  case SL_VEL_BC_NEUMANN_ALL:
    /* set Neumann for all points */
    return YMIR_VEL_DIRICHLET_NONE;
    break;

  default: /* unknown boundary condition */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* set one point to all Dirichlet */
/*  if (face == SL_BASE && fabs (nz - 1.0) < SC_EPS) {
    YMIR_INFOF ("Set boundary node %d with coord. (%g %g %g) and normal (%g %g %g) to all Dirichlet\n", node_id, X, Y, Z, nx, ny, nz);
    return YMIR_VEL_DIRICHLET_ALL;
  }*/

  /* set another point to all Dirichlet */
/*  if (face == SL_BASE && fabs (ny - 1.0) < SC_EPS) {
    YMIR_INFOF ("Set boundary node %d with coord. (%g %g %g) and normal (%g %g %g) to all Dirichlet\n", node_id, X, Y, Z, nx, ny, nz);
    return YMIR_VEL_DIRICHLET_ALL;
  }*/
}

/**
 * called by `slabs_shell_vel_dir_restrict_dof` below
 */
static void
slabs_shell_vel_dir_restrict_dof_fn (double *dummy,
                                     double X, double Y, double Z,
                                     double nx, double ny, double nz,
                                     ymir_topidx_t face, ymir_locidx_t node_id,
                                     void *data)
{
  ymir_vel_dir_t     *vel_dir = (ymir_vel_dir_t *) data;

  /* set one point to all Dirichlet */
  if (face == SL_BASE && fabs (nz - 1.0) < SC_EPS) {
    ymir_mesh_t        *mesh = vel_dir->mesh;
    ymir_face_mesh_t   *fmesh = &mesh->fmeshes[face];
    ymir_locidx_t       vcnid = fmesh->CntoVCn[node_id];
    double             *vn = ymir_cvec_index (vel_dir->nvec, vcnid, 0);

    /* set Dirichlet in all directions */
    vel_dir->rank[vcnid] = YMIR_VEL_DIRICHLET_ALL;
    vn[0] = 0;
    vn[1] = 0;
    vn[2] = 0;

    YMIR_INFOF ("Set boundary node %d with coord. (%g %g %g) "
                "and normal (%g %g %g) to all Dirichlet\n",
                node_id, X, Y, Z, nx, ny, nz);
  }

  /* set another point to Dirichlet in normal and one tangential direction */
  if (face == SL_BASE && fabs (ny - 1.0) < SC_EPS) {
    ymir_mesh_t        *mesh = vel_dir->mesh;
    ymir_face_mesh_t   *fmesh = &mesh->fmeshes[face];
    ymir_locidx_t       vcnid = fmesh->CntoVCn[node_id];
    double             *vn = ymir_cvec_index (vel_dir->nvec, vcnid, 0);

    /* set the way the vector `vn` is interpreted */
    vel_dir->rank[vcnid] = YMIR_VEL_DIRICHLET_TANG;

    /* set Dirichlet in y- and z-direction */
    vn[0] = 1;
    vn[1] = 0;
    vn[2] = 0;

    YMIR_INFOF ("Set boundary node %d with coord. (%g %g %g) "
                "and normal (%g %g %g) to Dirichlet in normal "
                "and one tangential direction\n",
                node_id, X, Y, Z, nx, ny, nz);
  }

}

void
slabs_shell_vel_dir_restrict_dof (ymir_vel_dir_t * vel_dir)
{
  ymir_mesh_t        *mesh = vel_dir->mesh;
  ymir_topidx_t       nface = mesh->num_face_meshes;
  ymir_topidx_t       fm;

  for (fm = 0; fm < nface; fm++) {      // loop over all faces
    ymir_cvec_t        *dummy_face_cvec;

    /* create dummy vector for current face */
    dummy_face_cvec = ymir_face_cvec_new (mesh, fm, 1);

    /* set additional Dirichlet boundary conditions for current face */
    ymir_face_cvec_set_function (dummy_face_cvec,
                                 slabs_shell_vel_dir_restrict_dof_fn, vel_dir);

    /* destroy vector */
    ymir_vec_destroy (dummy_face_cvec);
  }
}

/**
 * Gets maximum absolute value of all entries of a matrix.
 */
double
slabs_matrix_compute_abs_max (const sc_dmatrix_t *mat)
{
  const sc_bint_t     totalsize = mat->m * mat->n;
  const double       *data = mat->e[0];
  sc_bint_t           i;
  double              max = 0.0;

  /* find max */
  for (i = 0; i < totalsize; ++i) {
    max = SC_MAX (fabs (data[i]), max);
  }

  /* return max */
  return max;
}

/**
 * Gets minimum and maximum absolute value of all entries of a matrix.
 */
void
slabs_matrix_compute_abs_min_max (double *min, double *max,
                                  const sc_dmatrix_t *mat)
{
  const sc_bint_t     totalsize = mat->m * mat->n;
  const double       *data = mat->e[0];
  sc_bint_t           i;

  /* initialize min and max */
  *min = DBL_MAX;
  *max = 0.0;

  /* find min and max */
  for (i = 0; i < totalsize; ++i) {
    *min = SC_MIN (fabs (data[i]), *min);
    *max = SC_MAX (fabs (data[i]), *max);
  }
}

/**
 * Bounds all entries of a matrix to a minimum and a maximum.
 */
void
slabs_matrix_bound_values (sc_dmatrix_t *mat, const double min,
                           const double max)
{
  const int           totalsize = mat->m * mat->n;
  double             *data = mat->e[0];
  int                 i;

  for (i = 0; i < totalsize; i++) {
    data[i] = SC_MAX (data[i], min);
    data[i] = SC_MIN (data[i], max);
  }
}

/**
 * Rounds all entries of a matrix.
 */
void
slabs_matrix_round (sc_dmatrix_t *mat)
{
  const int           totalsize = mat->m * mat->n;
  double             *data = mat->e[0];
  int                 i;

  for (i = 0; i < totalsize; i++) {
    data[i] = round (data[i]);
  }
}

/**
 * Applies floor to all entries of a matrix.
 */
void
slabs_matrix_floor (sc_dmatrix_t *mat)
{
  const int           totalsize = mat->m * mat->n;
  double             *data = mat->e[0];
  int                 i;

  for (i = 0; i < totalsize; i++) {
    data[i] = floor (data[i]);
  }
}

/**
 * Applies ceil to all entries of a matrix.
 */
void
slabs_matrix_ceil (sc_dmatrix_t *mat)
{
  const int           totalsize = mat->m * mat->n;
  double             *data = mat->e[0];
  int                 i;

  for (i = 0; i < totalsize; i++) {
    data[i] = ceil (data[i]);
  }
}

/**
 * Perform element-wise multiplication, Y[i,:] := Y[i,:] .* X[i,0], for all i.
 */
void
slabs_matrix_multiply_in_1d (const sc_dmatrix_t *X, sc_dmatrix_t *Y)
{
  const sc_bint_t     nrows = Y->m;
  const sc_bint_t     ncols = Y->n;
  const double       *Xdata = X->e[0];
  double             *Ydata = Y->e[0];
  sc_bint_t           i, j;

  SC_ASSERT (X->m == Y->m && X->n == 1);

  for (i = 0; i < nrows; ++i) {
    for (j = 0; j < ncols; ++j) {
      Ydata[ncols * i + j] *= Xdata[i];
    }
  }
}

/**
 * Perform element-wise division, Y[i,:] := Y[i,:] ./ X[i,0], for all i.
 */
void
slabs_matrix_divide_in_1d (const sc_dmatrix_t *X, sc_dmatrix_t *Y)
{
  const sc_bint_t     nrows = Y->m;
  const sc_bint_t     ncols = Y->n;
  const double       *Xdata = X->e[0];
  double             *Ydata = Y->e[0];
  sc_bint_t           i, j;

  SC_ASSERT (X->m == Y->m && X->n == 1);

  for (i = 0; i < nrows; ++i) {
    for (j = 0; j < ncols; ++j) {
      Ydata[ncols * i + j] /= Xdata[i];
    }
  }
}

/**
 *
 */
 double
slabs_elem_volume (mangll_t *mangll, mangll_locidx_t elid)
{
  const int           N = ymir_n (mangll->N);
  const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);
  double             *mass = mangll->MassMatrix->e[0];
  double             *detJ = mangll->J->e[elid];
  int                 nodeid;
  double              volume = 0.0;

  /* compute volume of this element: 1^T * M * 1 */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    volume += mass[nodeid] * detJ[nodeid];
  }

  /* return volume of this element */
  return volume;
}

/**
 *
 */
 void
slabs_elem_get_gauss_coordinates (double *x, double *y, double *z,
                                  const mangll_locidx_t elid,
                                  mangll_t *mangll, double *tmp_el)
{
  const int           N = ymir_n (mangll->N);
  double             *ma_x = mangll->X->e[elid];
  double             *ma_y = mangll->Y->e[elid];
  double             *ma_z = mangll->Z->e[elid];
  double             *Br_data = mangll->Br->e[0];

  /* compute coordinates of Gauss nodes */
  MANGLL_TENSOR_IIAX_APPLY_ELEM (N + 1, Br_data, ma_x, x);
  MANGLL_TENSOR_IAIX_APPLY_ELEM (N + 1, Br_data, x, tmp_el);
  MANGLL_TENSOR_AIIX_APPLY_ELEM (N + 1, Br_data, tmp_el, x);

  MANGLL_TENSOR_IIAX_APPLY_ELEM (N + 1, Br_data, ma_y, y);
  MANGLL_TENSOR_IAIX_APPLY_ELEM (N + 1, Br_data, y, tmp_el);
  MANGLL_TENSOR_AIIX_APPLY_ELEM (N + 1, Br_data, tmp_el, y);

  MANGLL_TENSOR_IIAX_APPLY_ELEM (N + 1, Br_data, ma_z, z);
  MANGLL_TENSOR_IAIX_APPLY_ELEM (N + 1, Br_data, z, tmp_el);
  MANGLL_TENSOR_AIIX_APPLY_ELEM (N + 1, Br_data, tmp_el, z);
}

/**
 *
 */
void
slabs_interp_gll_to_gauss (sc_dmatrix_t *mat_gll, sc_dmatrix_t *mat_gauss,
                           mangll_t *mangll)
{
  const mangll_locidx_t  n_elements = mangll->mesh->K;
  const int           N = ymir_n (mangll->N);
  const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);
  const int           n_fields = mat_gll->n / n_nodes_per_el;
  double             *Br_data = mangll->refel->Br->e[0];
  double             *mat_gll_el;
  double             *mat_gauss_el;
  double             *tmp = YMIR_ALLOC (double, mat_gll->n);
  mangll_locidx_t     elid;

  YMIR_ASSERT (mat_gll->m == n_elements);
  YMIR_ASSERT (mat_gll->m == mat_gauss->m && mat_gll->n == mat_gauss->n);
  YMIR_ASSERT (sc_dmatrix_is_valid (mat_gll));

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    mat_gll_el = mat_gll->e[elid];
    mat_gauss_el = mat_gauss->e[elid];

    mangll_tensor_IIAx_apply_elem (N + 1, n_fields, Br_data,
                                   mat_gll_el, mat_gauss_el);
    mangll_tensor_IAIx_apply_elem (N + 1, n_fields, Br_data, mat_gauss_el, tmp);
    mangll_tensor_AIIx_apply_elem (N + 1, n_fields, Br_data, tmp, mat_gauss_el);
  }

  YMIR_ASSERT (sc_dmatrix_is_valid (mat_gauss));

  /* destroy */
  YMIR_FREE (tmp);
}

/**
 *
 */
void
slabs_interp_gauss_to_gll (sc_dmatrix_t *mat_gauss, sc_dmatrix_t *mat_gll,
                           mangll_t *mangll)
{
  const mangll_locidx_t  n_elements = mangll->mesh->K;
  const int           N = ymir_n (mangll->N);
  const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);
  const int           n_fields = mat_gauss->n / n_nodes_per_el;
  double             *Br_data = mangll->refel->Brinv->e[0];
  double             *mat_gll_el;
  double             *mat_gauss_el;
  double             *tmp = YMIR_ALLOC (double, mat_gll->n);
  mangll_locidx_t     elid;

  YMIR_ASSERT (mat_gauss->m == n_elements);
  YMIR_ASSERT (mat_gll->m == mat_gauss->m && mat_gll->n == mat_gauss->n);
  YMIR_ASSERT (sc_dmatrix_is_valid (mat_gauss));

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    mat_gauss_el = mat_gauss->e[elid];
    mat_gll_el = mat_gll->e[elid];

    mangll_tensor_IIAx_apply_elem (N + 1, n_fields, Br_data,
                                   mat_gauss_el, mat_gll_el);
    mangll_tensor_IAIx_apply_elem (N + 1, n_fields, Br_data, mat_gll_el, tmp);
    mangll_tensor_AIIx_apply_elem (N + 1, n_fields, Br_data, tmp, mat_gll_el);
  }

  YMIR_ASSERT (sc_dmatrix_is_valid (mat_gll));

  /* destroy */
  YMIR_FREE (tmp);
}

/**
 *
 */
static inline void
slabs_deriv_elem_geom_transf (double *_sc_restrict du,
                              const int order, const int n_fields,
                              double *_sc_restrict grad_u,
                              double *_sc_restrict ref_to_x,
                              double *_sc_restrict ref_to_y,
                              double *_sc_restrict ref_to_z)
{
  const int           N = ymir_n (order);
  const int           n_nodes = (N + 1) * (N + 1) * (N + 1);
  double              du_val;
  int                 nodeid, fieldid;

  for (nodeid = 0; nodeid < n_nodes; nodeid++) { /* loop over all nodes */
    for (fieldid = 0; fieldid < n_fields; fieldid++) { /* loop over all
                                                        * fields */
      du_val = *(du++);

      *(grad_u++) += ref_to_x[nodeid] * du_val;
      *(grad_u++) += ref_to_y[nodeid] * du_val;
      *(grad_u++) += ref_to_z[nodeid] * du_val;
    }
  }
}

/**
 * Computes the gradient on an element by the matrix vector product:
 *
 *   --                 --   --             --
 *   | dr/dx ds/dx ds/dx |   | I (x) I (x) D |
 *   | dr/dy ds/dy dt/dy | * | I (x) D (x) I | * u TODO might be wrong
 *   | dr/dz ds/dz dt/dz |   | D (x) I (x) I |
 *   --                 --   --             --
 *
 * Note: It is assumed that the result `grad_u` was initialized to zero before.
 */
static inline void
slabs_deriv_gll_to_gll_elem_grad (double *_sc_restrict u,
                                  const int order, int n_fields,
                                  double *_sc_restrict grad_u,
                                  double *_sc_restrict Dr,
                                  double *_sc_restrict rx,
                                  double *_sc_restrict sx,
                                  double *_sc_restrict tx,
                                  double *_sc_restrict ry,
                                  double *_sc_restrict sy,
                                  double *_sc_restrict ty,
                                  double *_sc_restrict rz,
                                  double *_sc_restrict sz,
                                  double *_sc_restrict tz,
                                  double *_sc_restrict tmp_du)
{
  const int           N = ymir_n (order);

  /*
   * [dr/dx dr/dy dr/dz] * ( [I (x) I (x) D] * u )^T
   */

  /* compute derivative in r-direction in reference domain */
  mangll_tensor_IIAx_apply_elem (N + 1, n_fields, Dr, u, tmp_du);

  /* apply chain rule to get gradient in physical domain */
  slabs_deriv_elem_geom_transf (tmp_du, order, n_fields, grad_u, rx, ry, rz);

  /*
   * [ds/dx ds/dy ds/dz] * ( [I (x) D (x) I] * u )^T
   */

  /* compute derivative in r-direction in reference domain */
  mangll_tensor_IAIx_apply_elem (N + 1, n_fields, Dr, u, tmp_du);

  /* apply chain rule to get gradient in physical domain */
  slabs_deriv_elem_geom_transf (tmp_du, order, n_fields, grad_u, sx, sy, sz);

  /*
   * [dt/dx dt/dy dt/dz] * ( [D (x) I (x) I] * u )^T
   */

  /* compute derivative in r-direction in reference domain */
  mangll_tensor_AIIx_apply_elem (N + 1, n_fields, Dr, u, tmp_du);

  /* apply chain rule to get gradient in physical domain */
  slabs_deriv_elem_geom_transf (tmp_du, order, n_fields, grad_u, tx, ty, tz);
}

/**
 *
 */
 void
slabs_gradient_gll_to_gll_elem (sc_dmatrix_t *u, sc_dmatrix_t *grad_u,
                                mangll_t *mangll, const mangll_locidx_t elid,
                                sc_dmatrix_t *tmp_du)
{
  const int           n_fields = u->n;
  const int           N = ymir_n (mangll->N);
#ifdef YMIR_DEBUG
  const int           n_nodes = (N + 1) * (N + 1) * (N + 1);
#endif
  double             *Dr_data = mangll->refel->Drll->e[0];
  double             *rx_data = mangll->rx->e[elid];
  double             *sx_data = mangll->sx->e[elid];
  double             *tx_data = mangll->tx->e[elid];
  double             *ry_data = mangll->ry->e[elid];
  double             *sy_data = mangll->sy->e[elid];
  double             *ty_data = mangll->ty->e[elid];
  double             *rz_data = mangll->rz->e[elid];
  double             *sz_data = mangll->sz->e[elid];
  double             *tz_data = mangll->tz->e[elid];

  /* check input parameters */
  YMIR_ASSERT (u->m == n_nodes);
  YMIR_ASSERT (grad_u->m == n_nodes && grad_u->n == (3 * n_fields));
  YMIR_ASSERT (tmp_du->m == n_nodes && tmp_du->n == n_fields);

  /* initialize gradient to zero */
  sc_dmatrix_set_zero (grad_u);

  /* compute derivative on GLL nodes for this element */
  slabs_deriv_gll_to_gll_elem_grad (u->e[0], N, n_fields, grad_u->e[0], Dr_data,
                                    rx_data, sx_data, tx_data,
                                    ry_data, sy_data, ty_data,
                                    rz_data, sz_data, tz_data,
                                    tmp_du->e[0]);
}

/**
 *
 */
static inline void
slabs_deriv_gll_to_gauss_elem_dr (double *_sc_restrict u,
                                  const int order, int n_fields,
                                  double *_sc_restrict Dr,
                                  double *_sc_restrict Br,
                                  double *_sc_restrict du,
                                  double *_sc_restrict tmp)
{
  const int           N = ymir_n (order);

  mangll_tensor_IIAx_apply_elem (N + 1, n_fields, Dr, u, du);
  mangll_tensor_IAIx_apply_elem (N + 1, n_fields, Br, du, tmp);
  mangll_tensor_AIIx_apply_elem (N + 1, n_fields, Br, tmp, du);
}

static inline void
slabs_deriv_gll_to_gauss_elem_ds (double *_sc_restrict u,
                                  const int order, int n_fields,
                                  double *_sc_restrict Dr,
                                  double *_sc_restrict Br,
                                  double *_sc_restrict du,
                                  double *_sc_restrict tmp)
{
  const int           N = ymir_n (order);

  mangll_tensor_IIAx_apply_elem (N + 1, n_fields, Br, u, du);
  mangll_tensor_IAIx_apply_elem (N + 1, n_fields, Dr, du, tmp);
  mangll_tensor_AIIx_apply_elem (N + 1, n_fields, Br, tmp, du);
}

static inline void
slabs_deriv_gll_to_gauss_elem_dt (double *_sc_restrict u,
                                  const int order, int n_fields,
                                  double *_sc_restrict Dr,
                                  double *_sc_restrict Br,
                                  double *_sc_restrict du,
                                  double *_sc_restrict tmp)
{
  const int           N = ymir_n (order);

  mangll_tensor_IIAx_apply_elem (N + 1, n_fields, Br, u, du);
  mangll_tensor_IAIx_apply_elem (N + 1, n_fields, Br, du, tmp);
  mangll_tensor_AIIx_apply_elem (N + 1, n_fields, Dr, tmp, du);
}

/**
 * Computes the gradient on an element by the matrix vector product:
 *
 *   --                 --   --             --
 *   | dr/dx ds/dx ds/dx |   | B (x) B (x) D |
 *   | dr/dy ds/dy dt/dy | * | B (x) D (x) B | * u
 *   | dr/dz ds/dz dt/dz |   | D (x) B (x) B |
 *   --                 --   --             --
 *
 * Note: It is assumed that the result `grad_u` was initialized to zero before.
 */
static inline void
slabs_deriv_gll_to_gauss_elem_grad (double *_sc_restrict u,
                                    const int order, int n_fields,
                                    double *_sc_restrict grad_u,
                                    double *_sc_restrict Dr,
                                    double *_sc_restrict Br,
                                    double *_sc_restrict rx,
                                    double *_sc_restrict sx,
                                    double *_sc_restrict tx,
                                    double *_sc_restrict ry,
                                    double *_sc_restrict sy,
                                    double *_sc_restrict ty,
                                    double *_sc_restrict rz,
                                    double *_sc_restrict sz,
                                    double *_sc_restrict tz,
                                    double *_sc_restrict tmp_du,
                                    double *_sc_restrict tmp)
{
  /*
   * [dr/dx dr/dy dr/dz] * ( [B (x) B (x) D] * u )^T
   */

  /* compute derivative in r-direction in reference domain */
  slabs_deriv_gll_to_gauss_elem_dr (u, order, n_fields, Dr, Br, tmp_du, tmp);

  /* apply chain rule to get gradient in physical domain */
  slabs_deriv_elem_geom_transf (tmp_du, order, n_fields, grad_u, rx, ry, rz);

  /*
   * [ds/dx ds/dy ds/dz] * ( [B (x) D (x) B] * u )^T
   */

  /* compute derivative in r-direction in reference domain */
  slabs_deriv_gll_to_gauss_elem_ds (u, order, n_fields, Dr, Br, tmp_du, tmp);

  /* apply chain rule to get gradient in physical domain */
  slabs_deriv_elem_geom_transf (tmp_du, order, n_fields, grad_u, sx, sy, sz);

  /*
   * [dt/dx dt/dy dt/dz] * ( [D (x) B (x) B] * u )^T
   */

  /* compute derivative in r-direction in reference domain */
  slabs_deriv_gll_to_gauss_elem_dt (u, order, n_fields, Dr, Br, tmp_du, tmp);

  /* apply chain rule to get gradient in physical domain */
  slabs_deriv_elem_geom_transf (tmp_du, order, n_fields, grad_u, tx, ty, tz);
}

/**
 *
 */
 void
slabs_gradient_gll_to_gauss_elem (sc_dmatrix_t *u, sc_dmatrix_t *grad_u,
                                  mangll_t *mangll, mangll_locidx_t elid,
                                  sc_dmatrix_t *tmp_du, sc_dmatrix_t *tmp)
{
  const int           n_fields = u->n;
  const int           N = ymir_n (mangll->N);
#ifdef YMIR_DEBUG
  const int           Np = (N + 1) * (N + 1) * (N + 1);
#endif
  double             *Dr_data = mangll->refel->Dr->e[0];
  double             *Br_data = mangll->refel->Br->e[0];
  double             *rx_data = mangll->rx->e[elid];
  double             *sx_data = mangll->sx->e[elid];
  double             *tx_data = mangll->tx->e[elid];
  double             *ry_data = mangll->ry->e[elid];
  double             *sy_data = mangll->sy->e[elid];
  double             *ty_data = mangll->ty->e[elid];
  double             *rz_data = mangll->rz->e[elid];
  double             *sz_data = mangll->sz->e[elid];
  double             *tz_data = mangll->tz->e[elid];

  /* check input parameters */
  YMIR_ASSERT (N > 1);
  YMIR_ASSERT (u->m == Np);
  YMIR_ASSERT (grad_u->m == Np && grad_u->n == (3 * n_fields));
  YMIR_ASSERT (tmp_du->m == Np && tmp_du->n == n_fields);
  YMIR_ASSERT (tmp->m == Np && tmp->n == n_fields);

  /* initialize gradient to zero */
  sc_dmatrix_set_zero (grad_u);

  /* compute derivative on Gauss nodes for this element */
  slabs_deriv_gll_to_gauss_elem_grad (u->e[0], N, n_fields, grad_u->e[0],
                                      Dr_data, Br_data,
                                      rx_data, sx_data, tx_data,
                                      ry_data, sy_data, ty_data,
                                      rz_data, sz_data, tz_data,
                                      tmp_du->e[0], tmp->e[0]);
}

/**
 *
 */
void
slabs_gradient_gll_to_gauss (sc_dmatrix_t *u, sc_dmatrix_t *grad_u,
                             mangll_t *mangll)
{
  const int           n_fields = u->n;
  mangll_locidx_t     n_elements = mangll->mesh->K;
  const int           N = ymir_n (mangll->N);
  const int           Np = (N + 1) * (N + 1) * (N + 1);
  double             *Dr_data = mangll->refel->Dr->e[0];
  double             *Br_data = mangll->refel->Br->e[0];
  sc_dmatrix_t       *rx = mangll->rx;
  sc_dmatrix_t       *sx = mangll->sx;
  sc_dmatrix_t       *tx = mangll->tx;
  sc_dmatrix_t       *ry = mangll->ry;
  sc_dmatrix_t       *sy = mangll->sy;
  sc_dmatrix_t       *ty = mangll->ty;
  sc_dmatrix_t       *rz = mangll->rz;
  sc_dmatrix_t       *sz = mangll->sz;
  sc_dmatrix_t       *tz = mangll->tz;
  double             *rx_data, *sx_data, *tx_data;
  double             *ry_data, *sy_data, *ty_data;
  double             *rz_data, *sz_data, *tz_data;

  double             *tmp_du = YMIR_ALLOC (double, Np * n_fields);
  double             *tmp = YMIR_ALLOC (double, Np * n_fields);
  double             *u_data, *grad_u_data;
  mangll_locidx_t     elid;

  /* check input parameters */
  YMIR_ASSERT (N > 1);
  YMIR_ASSERT (u->m == (n_elements * Np));
  YMIR_ASSERT (grad_u->m == (n_elements * Np) && grad_u->n == (3 * n_fields));

  /* initialize gradient to zero */
  sc_dmatrix_set_zero (grad_u);

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    /* get pointer to input and output vectors */
    u_data = u->e[Np * elid];
    grad_u_data = grad_u->e[Np * elid];

    /* get factors pertainting to geometric transformation */
    rx_data = rx->e[elid];
    sx_data = sx->e[elid];
    tx_data = tx->e[elid];
    ry_data = ry->e[elid];
    sy_data = sy->e[elid];
    ty_data = ty->e[elid];
    rz_data = rz->e[elid];
    sz_data = sz->e[elid];
    tz_data = tz->e[elid];

    /* compute derivative on Gauss nodes for this element */
    slabs_deriv_gll_to_gauss_elem_grad (u_data, N, n_fields, grad_u_data,
                                        Dr_data, Br_data,
                                        rx_data, sx_data, tx_data,
                                        ry_data, sy_data, ty_data,
                                        rz_data, sz_data, tz_data,
                                        tmp_du, tmp);
  }

  /* destroy */
  YMIR_FREE (tmp_du);
  YMIR_FREE (tmp);
}

/**
 *
 */
 void
slabs_second_invariant_elem (sc_dmatrix_t *u, sc_dmatrix_t *IIe,
                             mangll_t *mangll, mangll_locidx_t elid,
                             sc_dmatrix_t *tmp_grad_u,
                             sc_dmatrix_t *tmp_du,
                             sc_dmatrix_t *tmp,
                             slabs_node_type_t node_type)
{
  const int           N = ymir_n (mangll->N);
  const int           n_nodes = (N + 1) * (N + 1) * (N + 1);
  int                 nodeid;
  double              IIe_val;
  double             *_sc_restrict IIe_data = IIe->e[0];

  /* check input parameters */
  YMIR_ASSERT (u->m == n_nodes && u->n == 3);
  YMIR_ASSERT (IIe->m * IIe->n == n_nodes);
  YMIR_ASSERT (   node_type != SL_GAUSS_NODE
               || (node_type == SL_GAUSS_NODE && tmp != NULL) );

  /* compute gradient for this element */
  if (node_type == SL_GAUSS_NODE) {
    slabs_gradient_gll_to_gauss_elem (u, tmp_grad_u, mangll, elid, tmp_du, tmp);
  }
  else {
    slabs_gradient_gll_to_gll_elem (u, tmp_grad_u, mangll, elid, tmp_du);
  }

  /* compute 2nd invariant of the strain rate */
  for (nodeid = 0; nodeid < n_nodes; nodeid++) { /* loop over all nodes */
    double             *_sc_restrict grad_u_data = tmp_grad_u->e[0] + 9*nodeid;

    /* compute 2nd invariant for this node */
    IIe_val = 0.0;
    IIe_val += SC_SQR (grad_u_data[0]);
    IIe_val += SC_SQR (grad_u_data[1] + grad_u_data[3]) * 0.5;
    IIe_val += SC_SQR (grad_u_data[2] + grad_u_data[6]) * 0.5;
    IIe_val += SC_SQR (grad_u_data[4]);
    IIe_val += SC_SQR (grad_u_data[5] + grad_u_data[7]) * 0.5;
    IIe_val += SC_SQR (grad_u_data[8]);

    /* set 2nd invariant */
    IIe_data[nodeid] = 0.5 * IIe_val;
  }
}

/**
 *
 */
void
slabs_second_invariant (sc_dmatrix_t *u, sc_dmatrix_t *IIe,
                        mangll_t *mangll, slabs_node_type_t node_type)
{
  const int           n_fields = u->n;
  mangll_locidx_t     n_elements = mangll->mesh->K;
  const int           N = ymir_n (mangll->N);
  const int           n_nodes = (N + 1) * (N + 1) * (N + 1);

  sc_dmatrix_t       *tmp_grad_u = sc_dmatrix_new (n_nodes, 3 * n_fields);
  sc_dmatrix_t       *tmp_du = sc_dmatrix_new (n_nodes, n_fields);
  sc_dmatrix_t       *tmp;
  mangll_locidx_t     elid;

  /* check input parameters */
  YMIR_ASSERT (n_fields == 3);
  YMIR_ASSERT (u->m == (n_elements * n_nodes));
  YMIR_ASSERT ((IIe->m * IIe->n) == (n_elements * n_nodes));

  /* create additional variables for computations on Gauss nodes */
  if (node_type == SL_GAUSS_NODE) {
    tmp = sc_dmatrix_new (n_nodes, n_fields);
  }
  else {
    tmp = NULL;
  }

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    sc_dmatrix_t       *u_el = sc_dmatrix_new_view_offset (n_nodes * elid,
                                                           n_nodes, n_fields,
                                                           u);
    sc_dmatrix_t       *IIe_el = sc_dmatrix_new_view_offset (n_nodes * elid,
                                                             n_nodes, IIe->n,
                                                             IIe);

    /* compute 2nd invariant for this element */
    if (node_type == SL_GAUSS_NODE) {
      /* compute 2nd invariant on Gauss nodes */
      slabs_second_invariant_elem (u_el, IIe_el, mangll, elid, tmp_grad_u,
                                   tmp_du, tmp, SL_GAUSS_NODE);
    }
    else {
      /* compute 2nd invariant on GLL nodes */
      slabs_second_invariant_elem (u_el, IIe_el, mangll, elid, tmp_grad_u,
                                   tmp_du, tmp, SL_GLL_DISCONTINUOUS_NODE);
    }

    /* destroy */
    sc_dmatrix_destroy (u_el);
    sc_dmatrix_destroy (IIe_el);
  }

  /* destroy */
  sc_dmatrix_destroy (tmp_grad_u);
  sc_dmatrix_destroy (tmp_du);
  if (node_type == SL_GAUSS_NODE) {
    sc_dmatrix_destroy (tmp);
  }
}

/**
 *
 */
 void
slabs_strain_rate_tensor_elem (sc_dmatrix_t *strain_rate_tensor,
                               mangll_t *mangll,
                               sc_dmatrix_t *grad_u)
{
  const int           N = ymir_n (mangll->N);
  const int           n_nodes = (N + 1) * (N + 1) * (N + 1);
  int                 nodeid;

  /* check input parameters */
  YMIR_ASSERT (grad_u->m == n_nodes && grad_u->n == 9);
  YMIR_ASSERT (strain_rate_tensor->m == n_nodes && strain_rate_tensor->n == 6);

  /* compute the strain rate tensor */
  for (nodeid = 0; nodeid < n_nodes; nodeid++) { /* loop over all nodes */
    double             *_sc_restrict grad_u_data = grad_u->e[0] + 9 * nodeid;
    double             *_sc_restrict strain_rate_data =
                                       strain_rate_tensor->e[0] + 6 * nodeid;

    /* compute strain rate tensor for this node */
    strain_rate_data[0] = grad_u_data[0];
    strain_rate_data[1] = 0.5 * (grad_u_data[1] + grad_u_data[3]);
    strain_rate_data[2] = 0.5 * (grad_u_data[2] + grad_u_data[6]);
    strain_rate_data[3] = grad_u_data[4];
    strain_rate_data[4] = 0.5 * (grad_u_data[5] + grad_u_data[7]);
    strain_rate_data[5] = grad_u_data[8];
  }
}

/**
 *
 */
 int
slabs_strain_rate_tensor_idx (const int row, const int col)
{
  return   ( 3*SC_MIN (row, col) - (SC_MIN (row, col) - 1)*SC_MIN (row, col)/2 )
         + ( SC_MAX (row, col) - SC_MIN (row, col) );
}

/**
 *
 */
static void
slabs_vec_is_finite_fn (double *val, double x, double y, double z,
                        ymir_locidx_t nid, void *data)
{
  int                *finite = (int *) data;

  if (!isfinite (*val)) {
    *finite = 0;
    //TODO remove
    //YMIR_INFOF ("###DEV### nid %i, (x,y,z) = (%.3f,%.3f,%.3f), val %g "
    //            "not finite\n",
    //            nid, x, y, z, *val);
  }
}

static void
slabs_vec_is_finite_elem_fn (double *val, double *x, double *y, double *z,
                             ymir_locidx_t elem_id, void *data)
{
  int                *finite = (int *) data;

  if (!isfinite (*val)) {
    *finite = 0;
    //TODO remove
    //YMIR_INFOF ("###DEV### elid %i, (x,y,z) = (%.2f,%.2f,%.2f), val %g "
    //            "not finite\n",
    //            elem_id, *x, *y, *z, *val);
  }
  else {
    //TODO remove
    //YMIR_INFOF ("###DEV### elid %i, (x,y,z) = (%.2f,%.2f,%.2f), val %g\n",
    //            elem_id, *x, *y, *z, *val);
  }
}

/**
 *
 */
int
slabs_vec_is_finite (ymir_vec_t *vec)
{
  int                 finite_cvec = 1;
  int                 finite_dvec = 1;
  int                 finite_evec = 1;

  ymir_vec_set_function (vec,
                         slabs_vec_is_finite_fn, &finite_cvec,
                         slabs_vec_is_finite_fn, &finite_dvec,
                         slabs_vec_is_finite_elem_fn, &finite_evec);

  return finite_cvec && finite_dvec && finite_evec;
}

/**
 *
 */
void
slabs_cvec_set_min_value (ymir_cvec_t *cvec, const double min)
{
  const int           totalsize = cvec->cvec->m * cvec->cvec->n;
  double             *data = cvec->cvec->e[0];
  int                 i;

  for (i = 0; i < totalsize; i++) {
    if (isfinite (data[i])) {
      data[i] = SC_MIN (data[i], min);
    }
    else {
      data[i] = min;
    }
  }
}

/**
 *
 */
void
slabs_cvec_bound_values (ymir_vec_t *cvec, const double min,
                         const double max)
{
  YMIR_ASSERT_HAS_CVEC (cvec);
  slabs_matrix_bound_values (cvec->cvec, min, max);
}

void
slabs_dvec_bound_values (ymir_vec_t *dvec, const double min,
                         const double max)
{
  YMIR_ASSERT_HAS_DVEC (dvec);
  slabs_matrix_bound_values (dvec->dvec, min, max);
}

/**
 *
 */
void
slabs_cvec_round (ymir_vec_t *cvec)
{
  YMIR_ASSERT_HAS_CVEC (cvec);
  slabs_matrix_round (cvec->cvec);
}

void
slabs_dvec_round (ymir_vec_t *dvec)
{
  YMIR_ASSERT_HAS_DVEC (dvec);
  slabs_matrix_round (dvec->dvec);
}

/**
 *
 */
void
slabs_cvec_floor (ymir_vec_t *cvec)
{
  YMIR_ASSERT_HAS_CVEC (cvec);
  slabs_matrix_floor (cvec->cvec);
}

void
slabs_dvec_floor (ymir_vec_t *dvec)
{
  YMIR_ASSERT_HAS_DVEC (dvec);
  slabs_matrix_floor (dvec->dvec);
}

/**
 *
 */
void
slabs_cvec_ceil (ymir_vec_t *cvec)
{
  YMIR_ASSERT_HAS_CVEC (cvec);
  slabs_matrix_ceil (cvec->cvec);
}

void
slabs_dvec_ceil (ymir_vec_t *dvec)
{
  YMIR_ASSERT_HAS_DVEC (dvec);
  slabs_matrix_ceil (dvec->dvec);
}

/**
 *
 */
void
slabs_cvec_compute_magnitude (const ymir_vec_t *cvec, ymir_vec_t *magn)
{
  ymir_face_mesh_t   *fmesh = &cvec->mesh->fmeshes[cvec->meshnum];
  const ymir_locidx_t n_cnodes = fmesh->Ncn;
  const int           n_fields = cvec->ncfields;
  const double       *cvec_data = cvec->cvec->e[0];
  double             *magn_data = magn->cvec->e[0];
  ymir_locidx_t       cnid;
  int                 fieldid;
  double              sum;

  YMIR_ASSERT (&magn->mesh->fmeshes[magn->meshnum] == fmesh);
  YMIR_ASSERT (magn->ncfields == 1);

  if (YMIR_CVEC_STRIDE == YMIR_STRIDE_NODE) {
    for (cnid = 0; cnid < n_cnodes; cnid++) {
      sum = 0.0;
      for (fieldid = 0; fieldid < n_fields; fieldid++) {
        sum += SC_SQR (cvec_data[n_fields * cnid + fieldid]);
      }
      magn_data[cnid] = sqrt (sum);
    }
  }
  else {
    YMIR_ABORT_NOT_REACHED ();
  }
}

/**
 *
 */
void
slabs_dvec_compute_magnitude (const ymir_vec_t *dvec, ymir_vec_t *magn)
{
  const ymir_locidx_t n_elements = dvec->K;
  const int           n_nodes_per_el = dvec->Np;
  const int           n_fields = dvec->ndfields;
  const double       *dvec_data = dvec->dvec->e[0];
  double             *magn_data = magn->dvec->e[0];
  ymir_locidx_t       i;
  int                 fieldid;
  double              sum;

  YMIR_ASSERT (magn->K == n_elements);
  YMIR_ASSERT (magn->Np == n_nodes_per_el);
  YMIR_ASSERT (magn->ndfields == 1);

  if (YMIR_DVEC_STRIDE == YMIR_STRIDE_NODE) {
    for (i = 0; i < (n_elements * n_nodes_per_el); i++) {
      sum = 0.0;
      for (fieldid = 0; fieldid < n_fields; fieldid++) {
        sum += SC_SQR (dvec_data[n_fields * i + fieldid]);
      }
      magn_data[i] = sqrt (sum);
    }
  }
  else {
    YMIR_ABORT_NOT_REACHED ();
  }
}

/**
 *
 */
ymir_dvec_t *
slabs_dvec_new_from_element_data (sc_dmatrix_t *data,
                                  ymir_mesh_t *mesh,
                                  ymir_node_type_t node_type)
{
  sc_dmatrix_t       *eldata;
  ymir_locidx_t       elid;
  int                 nodeid;
  ymir_dvec_t        *dvec;

  /* create ymir dvec */
  dvec = ymir_dvec_new (mesh, 1, node_type);

  /* check input parameters */
  YMIR_ASSERT (data->m > 0 && data->n == dvec->K);

  /* create vector to store element values */
  eldata = sc_dmatrix_new (1, dvec->Np);

  for (elid = 0; elid < dvec->K; elid++) { /* loop over all elem's */
    /* set values for this element */
    ymir_dvec_get_elem (dvec, eldata, YMIR_STRIDE_COMP, elid, YMIR_WRITE);
    for (nodeid = 0; nodeid < dvec->Np; nodeid++) { /* loop over all nodes */
      eldata->e[0][nodeid] = data->e[0][elid];
    }
    ymir_dvec_set_elem (dvec, eldata, YMIR_STRIDE_COMP, elid, YMIR_SET);
  }

  /* destroy */
  sc_dmatrix_destroy (eldata);

  /* return ymir dvec */
  return dvec;
}

/**
 *
 */
void
slabs_stokes_vec_get_components_view (ymir_vec_t **u, ymir_vec_t **p,
                                      ymir_vec_t *up)
{
  ymir_mesh_t        *mesh = up->mesh;

  YMIR_ASSERT (1 < ymir_n (mesh->ma->N));

  /* create view of velocity vector */
  if (u != NULL) {
    *u = ymir_cvec_new_data (mesh, up->ncfields, up->cvec);
  }

  /* create view of pressure vector */
  if (p != NULL) {
    *p = ymir_evec_new_data (mesh, up->nefields, up->e_to_d_fn, up->e_to_d_data,
                             up->evec);
  }
}

/**
 *
 */
void
slabs_stokes_vec_get_components (ymir_vec_t **u, ymir_vec_t **p,
                                 ymir_vec_t *up,
                                 ymir_pressure_elem_t *press_elem)
{
  if (press_elem->space == YMIR_PRESSURE_SPACE_STAB) {
    ymir_mesh_t        *mesh = up->mesh;

    if (u != NULL && p != NULL) {
      *u = ymir_cvec_new (mesh, 3);
      *p = ymir_pressure_vec_new (mesh, press_elem);
      ymir_stokes_vec_get_components (up, *u, *p, press_elem);
    }
    else if (u != NULL) {
      *u = ymir_cvec_new (mesh, 3);
      ymir_stokes_vec_get_velocity (up, *u, press_elem);
    }
    else if (p != NULL) {
      *p = ymir_pressure_vec_new (mesh, press_elem);
      ymir_stokes_vec_get_pressure (up, *p, press_elem);
    }
  }
  else {
    slabs_stokes_vec_get_components_view (u, p, up);
  }
}

/**
 *
 */
void
slabs_stokes_vec_get_components_update (ymir_vec_t *u, ymir_vec_t *p,
                                        ymir_vec_t *up,
                                        ymir_pressure_elem_t *press_elem)
{
  if (press_elem->space == YMIR_PRESSURE_SPACE_STAB) {
    if (u != NULL && p != NULL) {
      ymir_stokes_vec_get_components (up, u, p, press_elem);
    }
    else if (u != NULL) {
      ymir_stokes_vec_get_velocity (up, u, press_elem);
    }
    else if (p != NULL) {
      ymir_stokes_vec_get_pressure (up, p, press_elem);
    }
  }
  /* otherwise nothing needs to be done due to view onto vector data */
}

/**
 *
 */
void
slabs_stokes_vec_set_components (ymir_vec_t *up,
                                 ymir_vec_t *u, ymir_vec_t *p,
                                 ymir_pressure_elem_t *press_elem)
{
  if (press_elem->space == YMIR_PRESSURE_SPACE_STAB) {
    if (u != NULL && p != NULL) {
      ymir_stokes_vec_set_components (u, p, up, press_elem);
    }
    else if (u != NULL) {
      ymir_stokes_vec_set_velocity (u, up, press_elem);
    }
    else if (p != NULL) {
      ymir_stokes_vec_set_pressure (p, up, press_elem);
    }
  }
  /* otherwise nothing needs to be done due to view onto vector data */
}

//TODO deprecated
#if 0
/* timing results for linear and nonlinear simulation */
typedef struct slabs_runtime
{
  double              setup;
  double              setup_amr;
  double              setup_stokes;

  double              solver;
  double              lin_solver_setup_stokes_op;
  double              lin_solver_setup_stokes_pc;
  double              nl_solver_setup_stokes_op;
  double              nl_solver_setup_stokes_pc;
  double              nl_solver_lin_solves;
}
slabs_runtime_t;

/**
 * Prints runtimes.
 */
void
slabs_print_timing (slabs_runtime_t *runtime,
                    const int nl_solver_type,
                    const int mpisize,
                    const long long int n_unknowns_total)
{
  const double        timing_norm = (double) (mpisize * n_unknowns_total);
  const double        total = runtime->setup + runtime->solver;

  /* print table header */
  YMIR_GLOBAL_INFO  ("===================== Runtime =====================\n");
  YMIR_GLOBAL_INFOF ("  N = #unknowns = %lli (at solution)\n",
                     n_unknowns_total);
  YMIR_GLOBAL_INFOF ("  p = #MPI procs = %i\n", mpisize);
  YMIR_GLOBAL_INFO  ("---------------------------------------------------\n");
  YMIR_GLOBAL_INFO  ("  Task                  time (sec)   time/p/N (sec)\n");
  YMIR_GLOBAL_INFO  ("---------------------------------------------------\n");

  /* print setup runtime */
  YMIR_GLOBAL_INFOF ("  Setup                 %1.3e\n", runtime->setup);
  YMIR_GLOBAL_INFOF ("  - init AMR            %1.3e\n", runtime->setup_amr);
  YMIR_GLOBAL_INFOF ("  - create Stokes prob  %1.3e\n", runtime->setup_stokes);

  /* print solver runtime */
  if (nl_solver_type == SL_NL_SOLVER_NONE) { /* if linear solve */
    YMIR_GLOBAL_INFOF ("  Linear solver setup   %1.3e    %1.3e\n",
                       runtime->lin_solver_setup_stokes_op +
                       runtime->lin_solver_setup_stokes_pc,
                       ( runtime->lin_solver_setup_stokes_op +
                         runtime->lin_solver_setup_stokes_pc ) / timing_norm);
    YMIR_GLOBAL_INFOF ("  - setup Stokes op     %1.3e    %1.3e\n",
                       runtime->lin_solver_setup_stokes_op,
                       runtime->lin_solver_setup_stokes_op / timing_norm);
    YMIR_GLOBAL_INFOF ("  - setup Stokes pc     %1.3e    %1.3e\n",
                       runtime->lin_solver_setup_stokes_pc,
                       runtime->lin_solver_setup_stokes_pc / timing_norm);
    YMIR_GLOBAL_INFOF ("  Linear solver         %1.3e    %1.3e\n",
                       runtime->solver, runtime->solver / timing_norm);
  }
  else { /* if nonlinear solve */
    YMIR_GLOBAL_INFOF ("  Nonlinear solver      %1.3e    %1.3e\n",
                       runtime->solver, runtime->solver / timing_norm);
    YMIR_GLOBAL_INFOF ("  - setup Stokes op     %1.3e    %1.3e\n",
                       runtime->nl_solver_setup_stokes_op,
                       runtime->nl_solver_setup_stokes_op / timing_norm);
    YMIR_GLOBAL_INFOF ("  - setup Stokes pc     %1.3e    %1.3e\n",
                       runtime->nl_solver_setup_stokes_pc,
                       runtime->nl_solver_setup_stokes_pc / timing_norm);
    YMIR_GLOBAL_INFOF ("  - linear solves       %1.3e    %1.3e\n",
                       runtime->nl_solver_lin_solves,
                       runtime->nl_solver_lin_solves / timing_norm);
  }

  /* print table footer */
  YMIR_GLOBAL_INFO  ("---------------------------------------------------\n");
  YMIR_GLOBAL_INFOF ("  Total time            %1.3e    %1.3e\n",
                     total, total / timing_norm);
  YMIR_GLOBAL_INFO  ("===================================================\n");
}
#endif

/******************************************************************************
 * Tests
 *****************************************************************************/

/* parameter list for test of 2nd invariant */
typedef struct slabs_base_test_strain_set_data
{
  double              cx, cy, cz;
}
slabs_base_test_strain_set_data_t;

/**
 *
 */
static void
slabs_base_test_strain_set_vel_fn (double *val, double x, double y, double z,
                                   ymir_locidx_t nid, void *data)
{
  slabs_base_test_strain_set_data_t  *d =
    (slabs_base_test_strain_set_data_t *) data;
  const double        cx = d->cx;
  const double        cy = d->cy;
  const double        cz = d->cz;

  val[0] = exp (M_PI * x);
  val[1] = cos (cx * M_PI * x) * cos (cy * M_PI * y) * cos (cz * M_PI * z);
  val[2] = sin (cx * M_PI * x) * sin (cy * M_PI * y) * sin (cz * M_PI * z);
}

/**
 *
 */
void
slabs_base_test_second_invariant (ymir_mesh_t *mesh)
{
  const char         *this_fn_name = "slabs_base_test_second_invariant";
  ymir_cvec_t        *vel_vec;
  ymir_dvec_t        *IIe_direct = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
  ymir_dvec_t        *IIe_elem = ymir_dvec_new (mesh, 1, YMIR_GAUSS_NODE);
  slabs_base_test_strain_set_data_t  data;
  double              norm_IIe;
  double              abs_error;

  /* create velocity vector */
  data.cx = 2.0;
  data.cy = 4.0;
  data.cz = 1.0;
  vel_vec = ymir_cvec_new (mesh, 3);
  ymir_cvec_set_function (vel_vec, slabs_base_test_strain_set_vel_fn, &data);

  /* compute second invariant of the strain rate with vector function */
  {
    ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (
                                          mesh->ma->N, mesh->ma->ompsize);

    ymir_second_invariant_vec (vel_vec, IIe_direct, vel_elem);
    ymir_velocity_elem_destroy (vel_elem);
  }

  /* compute second invariant of the strain rate element wise */
  {
    mangll_t           *mangll = mesh->ma;
    const mangll_locidx_t  n_elements = mangll->mesh->K;
    const int           N = ymir_n (mangll->N);
    const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);
    sc_dmatrix_t       *vel_el_mat;
    sc_dmatrix_t       *grad_vel_el_mat, *tmp_dvel, *tmp_3d;
    sc_dmatrix_t       *IIe_el_mat;
    mangll_locidx_t     elid;

    /* create work variables */
    vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
    grad_vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 9);
    tmp_dvel = sc_dmatrix_new (n_nodes_per_el, 3);
    tmp_3d = sc_dmatrix_new (n_nodes_per_el, 3);
    IIe_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);

    for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
      /* get velocity of this element at GLL nodes */
      ymir_cvec_get_elem_interp (vel_vec, vel_el_mat, YMIR_STRIDE_NODE, elid,
                                 YMIR_GLL_NODE, YMIR_READ);

      /* compute 2nd invariant for this element */
      slabs_second_invariant_elem (vel_el_mat, IIe_el_mat, mangll, elid,
                                   grad_vel_el_mat, tmp_dvel, tmp_3d,
                                   SL_GAUSS_NODE);

      /* set values of 2nd invariant */
      ymir_dvec_set_elem (IIe_elem, IIe_el_mat, YMIR_STRIDE_NODE, elid,
                          YMIR_SET);
    }

    /* destroy work variables */
    sc_dmatrix_destroy (vel_el_mat);
    sc_dmatrix_destroy (grad_vel_el_mat);
    sc_dmatrix_destroy (tmp_dvel);
    sc_dmatrix_destroy (tmp_3d);
    sc_dmatrix_destroy (IIe_el_mat);
  }

  /* check error */
  norm_IIe = ymir_dvec_norm (IIe_direct);
  ymir_vec_add (-1.0, IIe_direct, IIe_elem);
  abs_error = ymir_dvec_norm (IIe_elem);

  /* print error */
  if (0.0 < norm_IIe) {
    YMIR_GLOBAL_VERBOSEF ("%s: abs error %1.3e, rel error %1.3e\n",
                          this_fn_name, abs_error, abs_error / norm_IIe);
  }
  else {
    YMIR_GLOBAL_VERBOSEF ("%s: abs error %1.3e\n", this_fn_name, abs_error);
  }

  /* destroy */
  ymir_vec_destroy (vel_vec);
  ymir_vec_destroy (IIe_direct);
  ymir_vec_destroy (IIe_elem);
}

/**
 *
 */
void
slabs_base_test_strain_rate_tensor (ymir_mesh_t *mesh)
{
  const char         *this_fn_name = "slabs_base_test_strain_rate_tensor";
  ymir_cvec_t        *vel_vec;
  ymir_dvec_t        *strain_rate_direct = ymir_dvec_new (mesh, 6,
                                                          YMIR_GAUSS_NODE);
  ymir_dvec_t        *strain_rate_elem = ymir_dvec_new (mesh, 6,
                                                        YMIR_GAUSS_NODE);
  slabs_base_test_strain_set_data_t  data;
  double              norm_strain_rate;
  double              abs_error;

  /* create velocity vector */
  data.cx = 2.0;
  data.cy = 4.0;
  data.cz = 1.0;
  vel_vec = ymir_cvec_new (mesh, 3);
  ymir_cvec_set_function (vel_vec, slabs_base_test_strain_set_vel_fn, &data);

  /* compute strain rate tensor with vector function */
  ymir_velocity_strain_rate (vel_vec, strain_rate_direct, 0);

  /* compute strain rate tensor element wise */
  {
    mangll_t           *mangll = mesh->ma;
    const mangll_locidx_t  n_elements = mangll->mesh->K;
    const int           N = ymir_n (mangll->N);
    const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);
    sc_dmatrix_t       *vel_el_mat;
    sc_dmatrix_t       *grad_vel_el_mat, *tmp_dvel, *tmp_3d;
    sc_dmatrix_t       *IIe_el_mat;
    sc_dmatrix_t       *strain_vel_el_mat;
    mangll_locidx_t     elid;

    /* create work variables */
    vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 3);
    grad_vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 9);
    tmp_dvel = sc_dmatrix_new (n_nodes_per_el, 3);
    tmp_3d = sc_dmatrix_new (n_nodes_per_el, 3);
    IIe_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
    strain_vel_el_mat = sc_dmatrix_new (n_nodes_per_el, 6);

    for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
      /* get velocity of this element at GLL nodes */
      ymir_cvec_get_elem_interp (vel_vec, vel_el_mat, YMIR_STRIDE_NODE, elid,
                                 YMIR_GLL_NODE, YMIR_READ);

      /* compute strain rate tensor for this element */
      slabs_second_invariant_elem (vel_el_mat, IIe_el_mat, mangll, elid,
                                   grad_vel_el_mat, tmp_dvel, tmp_3d,
                                   SL_GAUSS_NODE);
      slabs_strain_rate_tensor_elem (strain_vel_el_mat, mangll,
                                     grad_vel_el_mat);

      /* set values of strain rate tensor */
      ymir_dvec_set_elem (strain_rate_elem, strain_vel_el_mat,
                          YMIR_STRIDE_NODE, elid, YMIR_SET);
    }

    /* destroy work variables */
    sc_dmatrix_destroy (vel_el_mat);
    sc_dmatrix_destroy (grad_vel_el_mat);
    sc_dmatrix_destroy (tmp_dvel);
    sc_dmatrix_destroy (tmp_3d);
    sc_dmatrix_destroy (IIe_el_mat);
    sc_dmatrix_destroy (strain_vel_el_mat);
  }

  /* check error */
  norm_strain_rate = ymir_dvec_norm (strain_rate_direct);
  ymir_vec_add (-1.0, strain_rate_direct, strain_rate_elem);
  abs_error = ymir_dvec_norm (strain_rate_elem);

  /* print error */
  if (0.0 < norm_strain_rate) {
    YMIR_GLOBAL_VERBOSEF ("%s: abs error %1.3e, rel error %1.3e\n",
                          this_fn_name, abs_error,
                          abs_error / norm_strain_rate);
  }
  else {
    YMIR_GLOBAL_VERBOSEF ("%s: abs error %1.3e\n", this_fn_name, abs_error);
  }

  /* destroy */
  ymir_vec_destroy (vel_vec);
  ymir_vec_destroy (strain_rate_direct);
  ymir_vec_destroy (strain_rate_elem);
}

