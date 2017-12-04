/***********************************************************
 * Post-processing for 2nd invariant stress, traction, .etc.
 ***********************************************************/

#include <slabtest_stress.h>
#include <rhea.h>
#include <ymir_derivative_elem.h>
#include <ymir_vec_ops.h>
#include <ymir_vec_getset.h>
#include <ymir_velocity_vec.h>
#include <ymir_stokes_pc.h>
#include <ymir_stokes_vec.h>
#include <ymir_interp_vec.h>
#include <ymir_mass_vec.h>

/* compute 9 components of the stress tensor at each element*/
void
slabs_stress_elem (sc_dmatrix_t * in, sc_dmatrix_t * out,
                   sc_dmatrix_t * visc, ymir_velocity_elem_t * vel_elem,
                   double *_sc_restrict rxd, double *_sc_restrict sxd,
                   double *_sc_restrict txd, double *_sc_restrict ryd,
                   double *_sc_restrict syd, double *_sc_restrict tyd,
                   double *_sc_restrict rzd, double *_sc_restrict szd,
                   double *_sc_restrict tzd, sc_dmatrix_t * drst, sc_dmatrix_t * brst)
{
  const int           N = ymir_n (vel_elem->N);
  const int           Np = ymir_n (vel_elem->Np);
  int                 gp, i, j, k, l;
  double              S[9];
  double             *_sc_restrict viscd  = visc->e[0];

  sc_dmatrix_t       *tempvec1 = vel_elem->tempvec1;
  sc_dmatrix_t       *tempvec2 = vel_elem->tempvec2;
  sc_dmatrix_t       *tempvec3 = vel_elem->tempvec3;
  sc_dmatrix_t       *temptens = vel_elem->temptens1;
  ymir_derivative_elem_grad (N, 3, drst->e[0], brst->e[0],
                             rxd, sxd, txd, ryd, syd, tyd,
                             rzd, szd, tzd, in, temptens, tempvec1, tempvec2,
                             tempvec3, 0);

  /* create stress from duvw/dxyz * viscosity */
  for (gp = 0; gp < Np; gp++) {
    double               *_sc_restrict temptensd = temptens->e[0] + 9 * gp;
    double               *_sc_restrict outd  = out->e[0] + 9 * gp;

    S[0] = temptensd[0];
    S[1] = (temptensd[1] + temptensd[3]) * (1. / 2.);
    S[2] = (temptensd[2] + temptensd[6]) * (1. / 2.);
    S[3] = S[1];
    S[4] = temptensd[4];
    S[5] = (temptensd[5] + temptensd[7]) * (1. / 2.);
    S[6] = S[2];
    S[7] = S[5];
    S[8] = temptensd[8];

    /* compute linear combination of viscous stress and 3x3 tensor */
    for (i = 0; i < 9; i++)  {
      outd[i] = 2.0 * viscd[gp] * S[i];
    }
  }
}

/* compute 9 components of the stress tensor*/
void
slabs_stress (ymir_cvec_t * vel, ymir_dvec_t * tau,
              ymir_dvec_t * visc, ymir_velocity_elem_t * vel_elem)
{
  ymir_mesh_t         *mesh = vel->mesh;
  ymir_locidx_t       elid;
  const int           N  = ymir_n (mesh->cnodes->N);
  const int           Np = ymir_np (mesh->cnodes->N);
  const int           K  = mesh->cnodes->K;
  sc_dmatrix_t       *elemin = sc_dmatrix_new (1, 3 * Np);
  sc_dmatrix_t       *elemout = sc_dmatrix_new (1,9 * Np);
  sc_dmatrix_t       *elemvisc = sc_dmatrix_new (1, Np);
  const char          *this_fn_name = "slabs_stress";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  ymir_dvec_set_zero (tau);

  for (elid = 0; elid < K; elid++)  {
    double             *_sc_restrict rxd = mesh->drdx->e[elid];
    double             *_sc_restrict sxd = mesh->dsdx->e[elid];
    double             *_sc_restrict txd = mesh->dtdx->e[elid];
    double             *_sc_restrict ryd = mesh->drdy->e[elid];
    double             *_sc_restrict syd = mesh->dsdy->e[elid];
    double             *_sc_restrict tyd = mesh->dtdy->e[elid];
    double             *_sc_restrict rzd = mesh->drdz->e[elid];
    double             *_sc_restrict szd = mesh->dsdz->e[elid];
    double             *_sc_restrict tzd = mesh->dtdz->e[elid];
    double             *_sc_restrict Jdetd = mesh->Jdet->e[elid];

    ymir_cvec_get_elem_interp (vel, elemin, YMIR_STRIDE_NODE, elid,
                               YMIR_GLL_NODE, YMIR_COPY);
    ymir_dvec_get_elem_interp (visc, elemvisc, YMIR_STRIDE_COMP, elid,
                              YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (tau, elemout, YMIR_STRIDE_NODE, elid,
                             YMIR_GAUSS_NODE, YMIR_WRITE);

    slabs_stress_elem (elemin, elemout,
                       elemvisc,
                       vel_elem, rxd, sxd, txd, ryd,
                       syd, tyd, rzd, szd, tzd, mesh->drst, mesh->brst);

    ymir_dvec_set_elem_interp (tau, elemout, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);

    ymir_read_view_release (elemvisc);
  }

  sc_dmatrix_destroy (elemin);
  sc_dmatrix_destroy (elemout);
  sc_dmatrix_destroy (elemvisc);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}


/* compute traction as well as normal/shear stress at each element
 * the procedure is computing strain_rate tensor, then traction on the plane,
 * then normal and shear stress*/
void
slabs_normal_shear_stress_elem (sc_dmatrix_t * in, sc_dmatrix_t * out1, sc_dmatrix_t * out2,
                            sc_dmatrix_t * out3, double * n_dir,
                  sc_dmatrix_t * visc_mat, ymir_velocity_elem_t * vel_elem,
                  double *_sc_restrict rxd, double *_sc_restrict sxd,
                  double *_sc_restrict txd, double *_sc_restrict ryd,
                  double *_sc_restrict syd, double *_sc_restrict tyd,
                  double *_sc_restrict rzd, double *_sc_restrict szd,
                  double *_sc_restrict tzd, sc_dmatrix_t * drst, sc_dmatrix_t * brst)
{
  const int           N = ymir_n (vel_elem->N);
  int                 gp, i, j, k, l;
  double              tempmatd[9];
  double              tempvec[3];
  double              temp;
  double             *_sc_restrict viscd = visc_mat->e[0];
  sc_dmatrix_t       *tempvec1 = vel_elem->tempvec1;
  sc_dmatrix_t       *tempvec2 = vel_elem->tempvec2;
  sc_dmatrix_t       *tempvec3 = vel_elem->tempvec3;
  sc_dmatrix_t       *temptens = vel_elem->temptens1;

  ymir_derivative_elem_grad (N, 3, drst->e[0], brst->e[0],
                             rxd, sxd, txd, ryd, syd, tyd,
                             rzd, szd, tzd, in, temptens, tempvec1, tempvec2,
                             tempvec3, 0);

  /* create stress from duvw/dxyz * viscosity */
  for (gp = 0, k = 0; k < N + 1; k++) {
    for (j = 0; j < N + 1; j++) {
      for (i = 0; i < N + 1; i++, gp++) {
        double               *_sc_restrict temptensd = temptens->e[0] + 9 * gp;
        double               *_sc_restrict normal = out1->e[0] + gp;
        double               *_sc_restrict shear  = out2->e[0] + gp;
        double               *_sc_restrict trac   = out3->e[0] + 3 * gp;

        tempmatd[0] = temptensd[0];
        tempmatd[1] = (temptensd[1] + temptensd[3]) * (1. / 2.);
        tempmatd[2] = (temptensd[2] + temptensd[6]) * (1. / 2.);
        tempmatd[3] = tempmatd[1];
        tempmatd[4] = temptensd[4];
        tempmatd[5] = (temptensd[5] + temptensd[7]) * (1. / 2.);
        tempmatd[6] = tempmatd[2];
        tempmatd[7] = tempmatd[5];
        tempmatd[8] = temptensd[8];
        *(trac++) = tempvec[0] = ( tempmatd[0] * n_dir[0]
                                 + tempmatd[1] * n_dir[1]
                                 + tempmatd[2] * n_dir[2]) * 2.0 * viscd[gp];
        *(trac++) = tempvec[1] = ( tempmatd[3] * n_dir[0]
                                 + tempmatd[4] * n_dir[1]
                                 + tempmatd[5] * n_dir[2]) * 2.0 * viscd[gp];
        *(trac++) = tempvec[2] = ( tempmatd[6] * n_dir[0]
                                 + tempmatd[7] * n_dir[1]
                                 + tempmatd[8] * n_dir[2]) * 2.0 * viscd[gp];
        *normal = temp = tempvec[0] * n_dir[0] + tempvec[1] * n_dir[1] + tempvec[2] * n_dir[2];
        *shear  = sqrt( tempvec[0] * tempvec[0] + tempvec[1] * tempvec[1] + tempvec[2] * tempvec[2]
                        - temp * temp);
      }
    }
  }
}

/* compute traction as well as normal and shear stress*/
void
slabs_normal_shear_stress (ymir_cvec_t * vel, ymir_dvec_t * n_tau, ymir_dvec_t * s_tau,
                        ymir_dvec_t * traction, double * n_dir,
                       ymir_dvec_t * visc, ymir_velocity_elem_t * vel_elem)
{
  ymir_mesh_t         *mesh = vel->mesh;
  ymir_locidx_t       elid;
  const int           N  = ymir_n (mesh->cnodes->N);
  const int           Np = ymir_np (mesh->cnodes->N);
  const int           K  = mesh->cnodes->K;
  sc_dmatrix_t       *elemin = sc_dmatrix_new (1, 3 * Np);
  sc_dmatrix_t       *elemout1 = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemout2 = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemout3 = sc_dmatrix_new (1, 3 * Np);
  sc_dmatrix_t       *elemvisc = sc_dmatrix_new (1, Np);

  ymir_dvec_set_zero (traction);
  ymir_dvec_set_zero (n_tau);
  ymir_dvec_set_zero (s_tau);

  for (elid = 0; elid < K; elid++)  {
    double             *_sc_restrict rxd = mesh->drdx->e[elid];
    double             *_sc_restrict sxd = mesh->dsdx->e[elid];
    double             *_sc_restrict txd = mesh->dtdx->e[elid];
    double             *_sc_restrict ryd = mesh->drdy->e[elid];
    double             *_sc_restrict syd = mesh->dsdy->e[elid];
    double             *_sc_restrict tyd = mesh->dtdy->e[elid];
    double             *_sc_restrict rzd = mesh->drdz->e[elid];
    double             *_sc_restrict szd = mesh->dsdz->e[elid];
    double             *_sc_restrict tzd = mesh->dtdz->e[elid];
    double             *_sc_restrict Jdetd = mesh->Jdet->e[elid];

    ymir_cvec_get_elem_interp (vel, elemin, YMIR_STRIDE_NODE, elid,
                               YMIR_GLL_NODE, YMIR_COPY);
    ymir_dvec_get_elem_interp (visc, elemvisc, YMIR_STRIDE_COMP, elid,
                              YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (n_tau, elemout1, YMIR_STRIDE_COMP, elid,
                             YMIR_GAUSS_NODE, YMIR_WRITE);
    ymir_dvec_get_elem_interp (s_tau, elemout2, YMIR_STRIDE_COMP, elid,
                             YMIR_GAUSS_NODE, YMIR_WRITE);

    slabs_normal_shear_stress_elem (elemin, elemout1, elemout2, elemout3, n_dir, elemvisc,
                            vel_elem, rxd, sxd, txd, ryd,
                            syd, tyd, rzd, szd, tzd, mesh->drst, mesh->brst);

    ymir_dvec_set_elem_interp (n_tau, elemout1, YMIR_STRIDE_COMP, elid,
                             YMIR_GAUSS_NODE, YMIR_SET);
    ymir_dvec_set_elem_interp (s_tau, elemout2, YMIR_STRIDE_COMP, elid,
                             YMIR_GAUSS_NODE, YMIR_SET);
    ymir_dvec_set_elem_interp (traction, elemout3, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);

    ymir_read_view_release (elemvisc);
  }

  sc_dmatrix_destroy (elemin);
  sc_dmatrix_destroy (elemout1);
  sc_dmatrix_destroy (elemout2);
  sc_dmatrix_destroy (elemout3);
  sc_dmatrix_destroy (elemvisc);
}

/* compute traction only, can be expanded to compute normal and shear stress.
 * The function basically duplicate slabs_stress,
 * It is an alternative approach that takes advantage of an existing
 * subroutine ymir_velocity_strain_rate and directly compute traction on each node*/
void
slabs_traction (ymir_cvec_t * vel, ymir_dvec_t *traction,
                  double *n_dir, ymir_dvec_t *visc)
{
  ymir_mesh_t        *mesh = vel->mesh;
  const int           N = ymir_n (mesh->cnodes->N);
  const int           Np = (N + 1)*(N+1)*(N+1);
  ymir_locidx_t       K = mesh->cnodes->K;
  ymir_dvec_t         *tau_tensor = ymir_dvec_new (mesh, 6,
                                                     YMIR_GAUSS_NODE);
  ymir_locidx_t       elid;
  sc_dmatrix_t       *an = sc_dmatrix_new (0, 0);
  sc_dmatrix_t       *bn = sc_dmatrix_new (0, 0);
  int                 i, j;

  ymir_velocity_strain_rate (vel, tau_tensor, 0);
  ymir_dvec_multiply_in1 (visc, tau_tensor);
  for (elid = 0; elid < K; elid++) {
    for (j = 0; j < Np; j++) {
      double              val = 0.;
      ymir_dvec_get_node (tau_tensor, an, elid, j, YMIR_READ);
      ymir_dvec_get_node (traction,   bn, elid, j, YMIR_WRITE);
      bn->e[0][0] = an->e[0][0] * n_dir[0] + an->e[0][1] * n_dir[1] + an->e[0][2] * n_dir[2];
      bn->e[0][1] = an->e[0][1] * n_dir[0] + an->e[0][3] * n_dir[1] + an->e[0][4] * n_dir[2];
      bn->e[0][2] = an->e[0][2] * n_dir[0] + an->e[0][4] * n_dir[1] + an->e[0][5] * n_dir[2];
      ymir_read_view_release (an);
      ymir_dvec_set_node (traction, bn, elid, j, YMIR_SET);
    }
  }
  sc_dmatrix_destroy (an);
  sc_dmatrix_destroy (bn);
}


static void
slabs_physics_normal_boundary_stress_fn (double *stress_norm,
                                         double x, double y, double z,
                                         double nx, double ny, double nz,
                                         ymir_topidx_t face,
                                         ymir_locidx_t node_id,
                                         void *data)
{
  ymir_vec_t         *vec_bndr = (ymir_vec_t *) data;
  double             *v = ymir_cvec_index (vec_bndr, node_id, 0);

  RHEA_ASSERT (vec_bndr->ncfields == 3);

  /* compute inner product with boundary outer normal vector */
  *stress_norm = nx * v[0] + ny * v[1] + nz * v[2];
}

/* normal stress on the surface */
void
slabs_physics_compute_normal_boundary_stress (ymir_vec_t *stress_bndr_norm,
                                              ymir_vec_t *up,
                                              ymir_vec_t *rhs_u_point,
                                              ymir_stokes_op_t *stokes_op)
{
  ymir_mesh_t        *mesh = up->mesh;
  ymir_pressure_elem_t  *press_elem = stokes_op->press_elem;
  ymir_stress_op_t   *stress_op = stokes_op->stress_op;
  const int           skip_dir = stress_op->skip_dir;
  const ymir_topidx_t face_id = stress_bndr_norm->meshnum;

  ymir_vec_t         *rhs = ymir_stokes_vec_new (mesh, press_elem);
  ymir_vec_t         *residual_up = ymir_stokes_vec_new (mesh, press_elem);
  ymir_vec_t         *residual_u = ymir_cvec_new (mesh, 3);
  ymir_vec_t         *residual_bndr = ymir_face_cvec_new (mesh, face_id, 3);
  ymir_vec_t         *mass_lump_boundary;

  /* check input */
  YMIR_ASSERT_IS_CVEC (stress_bndr_norm);
  RHEA_ASSERT (stress_bndr_norm->ncfields == 1);
  RHEA_ASSERT (ymir_stokes_vec_is_stokes_vec (up));
  YMIR_ASSERT_IS_CVEC (rhs_u_point);
  RHEA_ASSERT (ymir_vec_is_not_dirty (up));
  RHEA_ASSERT (ymir_vec_is_not_dirty (rhs_u_point));

  /* construct the right-hand side */
  ymir_stokes_pc_construct_rhs (rhs, rhs_u_point, NULL, NULL,
                                1 /* incompressible */, stokes_op, 0);
  RHEA_ASSERT (sc_dmatrix_is_valid (rhs->dataown));
  RHEA_ASSERT (sc_dmatrix_is_valid (rhs->coff));
  RHEA_ASSERT (ymir_vec_is_not_dirty (rhs));

  /* turn off boundary constraints */
  stress_op->skip_dir = 1;

  /* compute (unconstrained) residual
   *   r_mom  = a * u + b^t * p - f
   *   r_mass = b * u
   */
  ymir_stokes_pc_apply_stokes_op (up, residual_up, stokes_op, 0, 0);
  ymir_vec_add (-1.0, rhs, residual_up);
  RHEA_ASSERT (sc_dmatrix_is_valid (residual_up->dataown));
  RHEA_ASSERT (sc_dmatrix_is_valid (residual_up->coff));
  ymir_vec_destroy (rhs);

  /* restore boundary constraints */
  stress_op->skip_dir = skip_dir;

  /* get the velocity component of the residual */
  ymir_stokes_vec_get_velocity (residual_up, residual_u, stokes_op->press_elem);

  /* interpolate residual onto boundary */
  ymir_interp_vec (residual_u, residual_bndr);
  RHEA_ASSERT (sc_dmatrix_is_valid (residual_bndr->dataown));
  RHEA_ASSERT (sc_dmatrix_is_valid (residual_bndr->coff));
  ymir_vec_destroy (residual_up);
  ymir_vec_destroy (residual_u);

 /* get the normal part of the residual */
  ymir_face_cvec_set_function (stress_bndr_norm,
                               slabs_physics_normal_boundary_stress_fn,
                               residual_bndr);
  ymir_vec_destroy (residual_bndr);

  /* invert mass matrix on boundary */
  mass_lump_boundary = ymir_face_cvec_new (mesh, face_id, 1);
  ymir_mass_lump (mass_lump_boundary);
  ymir_vec_divide_in (mass_lump_boundary, stress_bndr_norm);
  ymir_vec_destroy (mass_lump_boundary);
}


