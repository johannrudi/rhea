/* Collide Example:
 *
 * Cartesian domain.  A constant viscosity is used with discontinuous jumps.
 * The viscosity is defined with the coordinate.  For example, in the lower
 * mantle, visc=100, and in the upper mantle, the viscosity is 1 except a weak
 * zone in the center where the vsic is 0.1.  The temperature should also be
 * defined with the location.  Same thing for the boundary condition.  The
 * normal velocity boundary condition is prefixed with the location.  For
 * example, on the left side, there can be inflow at the upper half and outflow
 * in the lower half.
 */

#include <rhea.h>
#include <ymir_velocity_vec.h>
#include <ymir_derivative_elem.h>
#include <ymir_vec_ops.h>
#include <ymir_stokes_pc.h>
#include <ymir_stokes_vec.h>
#include <ymir_stress_op.h>
#include <ymir_pressure_vec.h>


/* basic constants */
#define COLLIDE_SEC_PER_YEAR (31557600.0)    /* seconds in a year (365.25*24*3600) */
#define COLLIDE_EARTH_RADIUS (6371.0e3)      /* mean radius of the Earth [m]       */
#define COLLIDE_UPPER_MANTLE_DEPTH (660.0e3) /* approx. depth of upper mantle [m]  */
#define COLLIDE_THERM_DIFFUS (1.0e-6)        /* thermal diffusivity [m^2 / s]      */
#define COLLIDE_TEMP_DIFF (1400.0)           /* temperature difference [K]         */
#define COLLIDE_VISC_REP (1.0e20)            /* representative viscosity [Pa s]    */


/* enumerator for domain shapes */
typedef enum
{
  COLLIDE_VEL_DIR_BC_WALLSLIDE,
  COLLIDE_VEL_DIR_BC_INOUTFLOW_SIN,
  COLLIDE_VEL_DIR_BC_INOUTFLOW_TANH
}
collide_vel_dir_bc_t;

typedef enum
{
  COLLIDE_X_FUNCTION_IDENTITY,
  COLLIDE_X_FUNCTION_SINE,
  COLLIDE_X_FUNCTION_RAMP,
  COLLIDE_X_FUNCTION_PROFILE
}
collide_x_func_t;

/* options of collide example */
typedef struct collide_options
{
  collide_vel_dir_bc_t vel_dir_bc;
  collide_x_func_t     x_func;
  double              wallslide_vel;
  double              flow_scale;
  double              uwkzone_loc_upper;
  double              uwkzone_loc_lower;
}
collide_options_t;

/**
 * Geometry transformation.
 */

static void
collide_X_fn_identity (mangll_tag_t tag, mangll_locidx_t np,
                       const double *_sc_restrict EX,
                       const double *_sc_restrict EY,
                       const double *_sc_restrict EZ,
                       double *_sc_restrict X,
                       double *_sc_restrict Y,
                       double *_sc_restrict Z, void *data)
{
  mangll_locidx_t     il;

  for (il = 0; il < np; ++il) {
    X[il] = EX[il];
    Y[il] = EY[il];
    Z[il] = EZ[il];
  }
}

static void
collide_X_fn_sine (mangll_tag_t tag, mangll_locidx_t np,
                       const double *_sc_restrict EX,
                       const double *_sc_restrict EY,
                       const double *_sc_restrict EZ,
                       double *_sc_restrict X,
                       double *_sc_restrict Y,
                       double *_sc_restrict Z, void *data)
{
  mangll_locidx_t     il;
  const double        slope = 0.5;

  for (il = 0; il < np; ++il) {
    X[il] = EX[il];
    Y[il] = EY[il];
    Z[il] = EZ[il] * (1 + slope * sin(M_PI * EY[il]));
  }
}


static void
collide_X_fn_ramp (mangll_tag_t tag, mangll_locidx_t np,
                   const double *_sc_restrict EX,
                   const double *_sc_restrict EY,
                   const double *_sc_restrict EZ,
                   double *_sc_restrict X,
                   double *_sc_restrict Y,
                   double *_sc_restrict Z, void *data)
{
  const double        z_length = 1.0;
  const double        slope = 0.5;
  mangll_locidx_t     il;

  for (il = 0; il < np; ++il) {
    X[il] = EX[il];
    Y[il] = EY[il];
    Z[il] = EZ[il] * (z_length + slope * EY[il]);
  }
}

static void
collide_X_fn_profile (mangll_tag_t tag, mangll_locidx_t np,
                   const double *_sc_restrict EX,
                   const double *_sc_restrict EY,
                   const double *_sc_restrict EZ,
                   double *_sc_restrict X,
                   double *_sc_restrict Y,
                   double *_sc_restrict Z, void *data)
{
  mangll_locidx_t     il;
  FILE *fp;
  char input_s[180];
  double temp;

  fp = fopen("profile.topo","r");
  fgets(input_s,200,fp);
  for (il = 0; il < np; ++il) {
    fscanf(fp,"%lf\n",&temp);
    X[il] = EX[il];
    Y[il] = EY[il];
    Z[il] = EZ[il]*temp;
  }
}

///* compute traction. It is an alternative approach that takes advantage of an existing
//   subroutine ymir_velocity_strain_rate and directly compute traction on each node*/
//collide_traction (ymir_cvec_t * vel, ymir_dvec_t *traction,
//                  double *n_dir, ymir_dvec_t *visc)
//{
//  ymir_mesh_t        *mesh = vel->mesh;
//  const int           N = ymir_n (mesh->cnodes->N);
//  const int           Np = (N + 1)*(N+1)*(N+1);
//  ymir_locidx_t       K = mesh->cnodes->K;
//  ymir_dvec_t         *tau_tensor = ymir_dvec_new (mesh, 6,
//                                                     YMIR_GAUSS_NODE);
//
//
//  ymir_locidx_t       elid;
//  sc_dmatrix_t       *an = sc_dmatrix_new (0, 0);
//  sc_dmatrix_t       *bn = sc_dmatrix_new (0, 0);
//  int                 i, j;
//
//  ymir_velocity_strain_rate (vel, tau_tensor, 0);
//  ymir_dvec_multiply_in1 (visc, tau_tensor);
//  for (elid = 0; elid < K; elid++) {
//    for (j = 0; j < Np; j++) {
//      double              val = 0.;
//
//      ymir_dvec_get_node (tau_tensor, an, elid, j, YMIR_READ);
//      ymir_dvec_get_node (traction,   bn, elid, j, YMIR_WRITE);
//      bn->e[0][0] = an->e[0][0] * n_dir[0] + an->e[0][1] * n_dir[1] + an->e[0][2] * n_dir[2];
//      bn->e[0][1] = an->e[0][1] * n_dir[0] + an->e[0][3] * n_dir[1] + an->e[0][4] * n_dir[2];
//      bn->e[0][2] = an->e[0][2] * n_dir[0] + an->e[0][4] * n_dir[1] + an->e[0][5] * n_dir[2];
//      ymir_read_view_release (an);
//      ymir_dvec_set_node (traction, bn, elid, j, YMIR_SET);
//    }
//  }
//  sc_dmatrix_destroy (an);
//  sc_dmatrix_destroy (bn);
//
//}
//
///* compute traction as well as normal and shear stress using
// * exsiting subroutine (above) collide_traction
// * and directly projecting the traction to the nodes*/
//void
//collide_normal_stress (ymir_cvec_t * vel, ymir_dvec_t * n_tau, ymir_dvec_t * s_tau,
//                      ymir_dvec_t * traction,  double *n_dir,  ymir_dvec_t * visc)
//{
//  ymir_mesh_t        *mesh = vel->mesh;
//  const int           N = ymir_n (mesh->cnodes->N);
//  const int           Np = (N + 1)*(N+1)*(N+1);
//  ymir_locidx_t       K = mesh->cnodes->K;
//
//  ymir_locidx_t       elid;
//  sc_dmatrix_t       *an = sc_dmatrix_new (0, 0);
//  sc_dmatrix_t       *bn = sc_dmatrix_new (0, 0);
//  sc_dmatrix_t       *cn = sc_dmatrix_new (0, 0);
//  int                 i, j;
//
//  collide_traction (vel, traction, n_dir, visc);
//  for (elid = 0; elid < K; elid++) {
//    for (j = 0; j < Np; j++) {
//      ymir_dvec_get_node (traction, an, elid, j, YMIR_READ);
//      ymir_dvec_get_node (n_tau,    bn, elid, j, YMIR_WRITE);
//      ymir_dvec_get_node (s_tau,    cn, elid, j, YMIR_WRITE);
//      bn->e[0][0] = an->e[0][0] * n_dir[0] + an->e[0][1] * n_dir[1] + an->e[0][2] * n_dir[2];
//      cn->e[0][0] = sqrt(
//                          (   (an->e[0][0]*an->e[0][0])
//                            + (an->e[0][1]*an->e[0][1])
//                            + (an->e[0][2]*an->e[0][2])  )
//                            - (bn->e[0][0]*bn->e[0][0])
//                        );
//      ymir_read_view_release (an);
//      ymir_dvec_set_node (n_tau, bn, elid, j, YMIR_SET);
//      ymir_dvec_set_node (s_tau, cn, elid, j, YMIR_SET);
//    }
//  }
//  sc_dmatrix_destroy (an);
//  sc_dmatrix_destroy (bn);
//}



///* compute traction at each element*/
//void
//collide_traction_elem (sc_dmatrix_t * in, sc_dmatrix_t * out,
//                  double * n_dir,
//                  sc_dmatrix_t * visc_mat, ymir_velocity_elem_t * vel_elem,
//                  double *_sc_restrict rxd, double *_sc_restrict sxd,
//                  double *_sc_restrict txd, double *_sc_restrict ryd,
//                  double *_sc_restrict syd, double *_sc_restrict tyd,
//                  double *_sc_restrict rzd, double *_sc_restrict szd,
//                  double *_sc_restrict tzd, sc_dmatrix_t * drst, sc_dmatrix_t * brst)
//{
//  const int           N = ymir_n (vel_elem->N);
//  int                 gp, i, j, k, l;
//  double              tempmatd[9];
//  double             *_sc_restrict viscd = visc_mat->e[0];
//  double              factor;
//  sc_dmatrix_t       *tempvec1 = vel_elem->tempvec1;
//  sc_dmatrix_t       *tempvec2 = vel_elem->tempvec2;
//  sc_dmatrix_t       *tempvec3 = vel_elem->tempvec3;
//  sc_dmatrix_t       *temptens = vel_elem->temptens1;
//
//  ymir_derivative_elem_grad (N, 3, drst->e[0], brst->e[0],
//                             rxd, sxd, txd, ryd, syd, tyd,
//                             rzd, szd, tzd, in, temptens, tempvec1, tempvec2,
//                             tempvec3, 0);
//
//  /* create stress from duvw/dxyz * viscosity, multiply by gauss
//  * weights and Jdet */
//  for (gp = 0, k = 0; k < N + 1; k++) {
//    for (j = 0; j < N + 1; j++) {
//      for (i = 0; i < N + 1; i++, gp++) {
//        double               *_sc_restrict temptensd = temptens->e[0] + 9 * gp;
//        double               *_sc_restrict tempvec = out->e[0] + 3 * gp;
//
//        tempmatd[0] = temptensd[0];
//        tempmatd[1] = (temptensd[1] + temptensd[3]) * (1. / 2.);
//        tempmatd[2] = (temptensd[2] + temptensd[6]) * (1. / 2.);
//        tempmatd[3] = tempmatd[1];
//        tempmatd[4] = temptensd[4];
//        tempmatd[5] = (temptensd[5] + temptensd[7]) * (1. / 2.);
//        tempmatd[6] = tempmatd[2];
//        tempmatd[7] = tempmatd[5];
//        tempmatd[8] = temptensd[8];
//        *(tempvec++) = (tempmatd[0]*n_dir[0] + tempmatd[1]*n_dir[1] + tempmatd[2]*n_dir[2])*viscd[gp];
//        *(tempvec++) = (tempmatd[3]*n_dir[0] + tempmatd[4]*n_dir[1] + tempmatd[5]*n_dir[2])*viscd[gp];
//        *(tempvec++) = (tempmatd[6]*n_dir[0] + tempmatd[7]*n_dir[1] + tempmatd[8]*n_dir[2])*viscd[gp];
//      }
//    }
//  }
//
//}
//
///* compute traction*/
//void
//collide_traction (ymir_cvec_t * vel, ymir_dvec_t * traction, double * n_dir,
//                 ymir_dvec_t * visc, ymir_velocity_elem_t * vel_elem)
//{
//  ymir_mesh_t         *mesh = vel->mesh;
//  ymir_dvec_t         *tau_tensor = ymir_dvec_new (mesh, 6,
//                                                     YMIR_GAUSS_NODE);
//  ymir_locidx_t       elid;
//  const int           N  = ymir_n (mesh->cnodes->N);
//  const int           Np = ymir_np (mesh->cnodes->N);
//  const int           K  = mesh->cnodes->K;
//  sc_dmatrix_t       *elemin = sc_dmatrix_new (1, 3 * Np);
//  sc_dmatrix_t       *elemout = sc_dmatrix_new (1, 3 * Np);
//  sc_dmatrix_t       *elemvisc = sc_dmatrix_new (1, Np);
//
//  ymir_dvec_set_zero (traction);
//
//    for (elid = 0; elid < K; elid++)  {
//
//      double             *_sc_restrict rxd = mesh->drdx->e[elid];
//      double             *_sc_restrict sxd = mesh->dsdx->e[elid];
//      double             *_sc_restrict txd = mesh->dtdx->e[elid];
//      double             *_sc_restrict ryd = mesh->drdy->e[elid];
//      double             *_sc_restrict syd = mesh->dsdy->e[elid];
//      double             *_sc_restrict tyd = mesh->dtdy->e[elid];
//      double             *_sc_restrict rzd = mesh->drdz->e[elid];
//      double             *_sc_restrict szd = mesh->dsdz->e[elid];
//      double             *_sc_restrict tzd = mesh->dtdz->e[elid];
//      double             *_sc_restrict Jdetd = mesh->Jdet->e[elid];
//
//      ymir_cvec_get_elem_interp (vel, elemin, YMIR_STRIDE_NODE, elid,
//                                 YMIR_GLL_NODE, YMIR_COPY);
//      ymir_dvec_get_elem_interp (visc, elemvisc, YMIR_STRIDE_COMP, elid,
//                                YMIR_GAUSS_NODE, YMIR_READ);
//      collide_traction_elem (elemin, elemout, n_dir, elemvisc, vel_elem, rxd, sxd, txd, ryd,
//                        syd, tyd, rzd, szd, tzd, mesh->drst, mesh->brst);
//      ymir_dvec_set_elem_interp (traction, elemout, YMIR_STRIDE_NODE, elid,
//                                 YMIR_GAUSS_NODE, YMIR_SET);
//      ymir_read_view_release (elemvisc);
//    }
//
//  sc_dmatrix_destroy (elemin);
//  sc_dmatrix_destroy (elemout);
//  sc_dmatrix_destroy (elemvisc);
//}
//

/* compute traction as well as normal/shear stress at each element*/
void
collide_normal_stress_elem (sc_dmatrix_t * in, sc_dmatrix_t * out1, sc_dmatrix_t * out2,
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

  /* create stress from duvw/dxyz * viscosity, multiply by gauss
  * weights and Jdet */
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
                                 + tempmatd[2] * n_dir[2]) * viscd[gp];
        *(trac++) = tempvec[1] = ( tempmatd[3] * n_dir[0]
                                 + tempmatd[4] * n_dir[1]
                                 + tempmatd[5] * n_dir[2]) * viscd[gp];
        *(trac++) = tempvec[2] = ( tempmatd[6] * n_dir[0]
                                 + tempmatd[7] * n_dir[1]
                                 + tempmatd[8] * n_dir[2]) * viscd[gp];
        *normal = temp = tempvec[0] * n_dir[0] + tempvec[1] * n_dir[1] + tempvec[2] * n_dir[2];
        *shear  = sqrt( tempvec[0] * tempvec[0] + tempvec[1] * tempvec[1] + tempvec[2] * tempvec[2]
                        - temp * temp);

      }
    }
  }

}

/* compute traction as well as normal and shear stress*/
void
collide_normal_stress (ymir_cvec_t * vel, ymir_dvec_t * n_tau, ymir_dvec_t * s_tau,
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

      collide_normal_stress_elem (elemin, elemout1, elemout2, elemout3, n_dir, elemvisc,
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



/**
 *
 */
static void
collide_physics_normal_boundary_stress_fn (double *stress_norm,
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

void
collide_physics_compute_normal_boundary_stress (ymir_vec_t *stress_bndr_norm,
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
  ymir_stokes_op_construct_rhs_ext (rhs_u_point, NULL, NULL, rhs,
                                    1 /* incompressible */, stokes_op);
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
                               collide_physics_normal_boundary_stress_fn,
                               residual_bndr);
  ymir_vec_destroy (residual_bndr);

  /* invert mass matrix on boundary */
  mass_lump_boundary = ymir_face_cvec_new (mesh, face_id, 1);
  ymir_mass_lump (mass_lump_boundary);
  ymir_vec_divide_in (mass_lump_boundary, stress_bndr_norm);
  ymir_vec_destroy (mass_lump_boundary);
}

/* TODO compute integral of stress on specified area*/
/*
void collide_physics_compute_integral (ymir_vec_t *stress) {
  ymir_vec_t *mass = ymir_vec_template (stress);
  ymir_vec_t *ones = ymir_vec_template (stress);
  double temp_stress;

  ymir_mass_apply (stress,mass);
  ymir_vec_set_value (ones, 1.0);
  temp_stress = ymir_vec_innerprod (ones, mass);
  YMIR_GLOBAL_INFOF ("integral of stress in weakzone: %1.3e\n",temp_stress);

}
*/

/**
 * Sets up the mesh.
 */
static void
collide_setup_mesh (p4est_t **p4est,
                    ymir_mesh_t **ymir_mesh,
                    ymir_pressure_elem_t **press_elem,
                    MPI_Comm mpicomm,
                    rhea_domain_options_t *domain_options,
                    rhea_discretization_options_t *discr_options,
                    collide_options_t *collide_options)

{
  const char         *this_fn_name = "collide_setup_mesh";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* set custom X-function */
  switch (collide_options->x_func) {
    case COLLIDE_X_FUNCTION_IDENTITY:
      rhea_discretization_set_user_X_fn (discr_options,
                                     collide_X_fn_identity, NULL);
      break;

    case COLLIDE_X_FUNCTION_SINE:
      rhea_discretization_set_user_X_fn (discr_options,
                                       collide_X_fn_sine, NULL);
      break;

    case COLLIDE_X_FUNCTION_RAMP:
      rhea_discretization_set_user_X_fn (discr_options,
                                       collide_X_fn_ramp, NULL);
      break;

    case COLLIDE_X_FUNCTION_PROFILE:
      rhea_discretization_set_user_X_fn (discr_options,
                                       collide_X_fn_profile, NULL);
      break;

    default:
      RHEA_ABORT_NOT_REACHED ();
    break;
  }

  /* create p4est */
  *p4est = rhea_discretization_p4est_new (mpicomm, discr_options,
                                          domain_options);

  /* set up boundary, store in `discr_options` */
  rhea_discretization_options_set_boundary (discr_options, *p4est,
                                            domain_options);

  /* create ymir mesh and pressure element */
  rhea_discretization_ymir_mesh_new_from_p4est (ymir_mesh, press_elem, *p4est,
                                                discr_options);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/* Dirichlet all and tangential on one side: A test case 'wallslide'*/

static ymir_dir_code_t
collide_set_vel_dir_freeslip (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
    return YMIR_VEL_DIRICHLET_NORM;
}

static ymir_dir_code_t
collide_set_vel_dir_wallslide (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
  if (face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE3 && 0.5 < Z) {
    return YMIR_VEL_DIRICHLET_ALL;
  }
  else {
    return YMIR_VEL_DIRICHLET_NORM;
  }
}

void
collide_set_rhs_vel_nonzero_dir_wallslide (
    double *vel, double x, double y, double z,
    ymir_locidx_t nid, void *data)
{
  collide_options_t  *collide_options = data;
  const double        wallslide_vel = collide_options->wallslide_vel;

  if (fabs (y) < SC_1000_EPS && 0.5 < z) {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = wallslide_vel * sin (2.0 * M_PI * z);
  }
  else {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }
}

/*
 * Dirichlet all on one side of the domain: SIDE3 (y=0).
 */
static ymir_dir_code_t
collide_set_vel_dir_inoutflow (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
  if (face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE3) {
    return YMIR_VEL_DIRICHLET_ALL;
  }
  else if (face == RHEA_DOMAIN_BOUNDARY_FACE_BASE) {
     return YMIR_VEL_DIRICHLET_NONE;
  }
  else {
    return YMIR_VEL_DIRICHLET_NORM;
  }
}

/**
 * In-out flow sine velocity on one side of the domain.
 */
void
collide_set_rhs_vel_nonzero_dir_inoutflow_sin (
    double *vel, double x, double y, double z,
    ymir_locidx_t nid, void *data)
{
  collide_options_t  *collide_options = data;
  const double        flow_scale = collide_options->flow_scale;

  if (fabs (y) < SC_1000_EPS) {
    vel[0] = 0.0;
    vel[1] = -flow_scale * sin (2.0 * M_PI * z);
    vel[2] = 0.0;
  }
  else {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }
}

/**
 * In-out flow tanh velocity on one side of the domain.
 */
void
collide_set_rhs_vel_nonzero_dir_inoutflow_tanh (
    double *vel, double x, double y, double z,
    ymir_locidx_t nid, void *data)
{
  collide_options_t  *collide_options = data;
  const double        flow_scale = collide_options->flow_scale;
  const double        zU = collide_options->uwkzone_loc_upper;
  const double        zL = collide_options->uwkzone_loc_lower;
  const double        a = zU-zL, b = 1.0-a, c = 0.5*(zU+zL);
  const double        shape = 5.0, scaling = 0.5*(b+a)*flow_scale, shift = 0.5*(b-a)*flow_scale;
  double txL = shape*(z-zL),txU = shape*(z-zU);

  if (fabs (y) < SC_1000_EPS) {
    vel[0] = 0.0;
    vel[2] = 0.0;


    if (z<=c)
      vel[1] = shift + scaling *
             ( (exp (txL) - exp (-txL)) /
             (exp (txL) + exp (-txL)) );
    else
      vel[1] = shift - scaling *
             ( (exp (txU) - exp (-txU)) /
             (exp (txU) + exp (-txU)) );
  }
  else {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }
}

/**
 * Sets up a linear Stokes problem.
 */
static void
collide_setup_stokes (rhea_stokes_problem_t **stokes_problem,
                      ymir_mesh_t *ymir_mesh,
                      ymir_pressure_elem_t *press_elem,
                      rhea_domain_options_t *domain_options,
                      rhea_temperature_options_t *temp_options,
                      rhea_viscosity_options_t *visc_options,
                      collide_options_t *collide_options,
                      const char *vtk_write_input_path)
{
  const char         *this_fn_name = "collide_setup_stokes";
  ymir_vec_t         *temperature, *weakzone;
  ymir_vec_t         *rhs_vel_nonzero_dirichlet = NULL;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* compute temperature */
  temperature = rhea_temperature_new (ymir_mesh);
  //TODO xi, is the following computation of temperature what you want?
  rhea_temperature_compute (temperature, temp_options);

  /* compute weak zone */
  weakzone = rhea_weakzone_new (ymir_mesh);
  //TODO xi, here you need to provide the weak zone factor
  ymir_vec_set_value (weakzone, 1.0);

  /* set velocity boundary conditions & nonzero Dirichlet values */
  if (domain_options->velocity_bc_type == RHEA_DOMAIN_VELOCITY_BC_USER) {
    switch (collide_options->vel_dir_bc) {
      case COLLIDE_VEL_DIR_BC_WALLSLIDE:
        rhea_domain_set_user_velocity_dirichlet_bc (
            collide_set_vel_dir_wallslide, NULL /* no data necessary */,
            0 /* TODO don't need this flag */);
        rhs_vel_nonzero_dirichlet = rhea_velocity_new (ymir_mesh);
        ymir_cvec_set_function (rhs_vel_nonzero_dirichlet,
                                collide_set_rhs_vel_nonzero_dir_wallslide,
                                collide_options);
        break;

      case COLLIDE_VEL_DIR_BC_INOUTFLOW_SIN:
        rhea_domain_set_user_velocity_dirichlet_bc (
            collide_set_vel_dir_inoutflow, NULL /* no data necessary */,
            0 /* TODO don't need this flag */);
        rhs_vel_nonzero_dirichlet = rhea_velocity_new (ymir_mesh);
        ymir_cvec_set_function (rhs_vel_nonzero_dirichlet,
                                collide_set_rhs_vel_nonzero_dir_inoutflow_sin,
                                collide_options);
        break;

      case COLLIDE_VEL_DIR_BC_INOUTFLOW_TANH:
        rhea_domain_set_user_velocity_dirichlet_bc (
            collide_set_vel_dir_inoutflow, NULL /* no data necessary */,
            0 /* TODO don't need this flag */);
        rhs_vel_nonzero_dirichlet = rhea_velocity_new (ymir_mesh);
        ymir_cvec_set_function (rhs_vel_nonzero_dirichlet,
                                collide_set_rhs_vel_nonzero_dir_inoutflow_tanh,
                                collide_options);
        break;
      default: /* BC not set */
        RHEA_ABORT_NOT_REACHED ();
    }
  }
  /* write vtk of input data */ //TODO better move this into main fnc
  if (vtk_write_input_path != NULL) {
    const rhea_viscosity_t  viscosity_type = visc_options->type;
    ymir_vec_t         *background_temp = rhea_temperature_new (ymir_mesh);
    ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *rhs_vel = rhea_velocity_new (ymir_mesh);

    rhea_temperature_background_compute (background_temp, temp_options);
    if (viscosity_type == RHEA_VISCOSITY_NONLINEAR) {
      visc_options->type = RHEA_VISCOSITY_LINEAR;
    }
    rhea_viscosity_compute (viscosity,
                            NULL /* nl. Stokes output */,
                            NULL /* nl. Stokes output */,
                            NULL /* nl. Stokes output */,
                            temperature, weakzone,
                            NULL /* nl. Stokes input */,
                            visc_options);
    if (viscosity_type == RHEA_VISCOSITY_NONLINEAR) {
      visc_options->type = viscosity_type;
    }
    rhea_temperature_compute_rhs_vel (rhs_vel, temperature, temp_options);

    rhea_vtk_write_input_data (vtk_write_input_path, temperature,
                               background_temp, weakzone, viscosity, rhs_vel);
/*
    if (rhs_vel_nonzero_dirichlet != NULL) {
      char            path[BUFSIZ];

      snprintf (path, BUFSIZ, "%s_vel_nonzero_dirichlet", vtk_write_input_path);
      ymir_vtk_write (ymir_mesh, path, rhs_vel_nonzero_dirichlet,
                      "vel_nonzero_dirichlet", NULL);
    }
*/
    rhea_temperature_destroy (background_temp);
    rhea_viscosity_destroy (viscosity);
    rhea_velocity_destroy (rhs_vel);
  }

  /* create Stokes problem */
  *stokes_problem = rhea_stokes_problem_new (
      temperature, weakzone, rhs_vel_nonzero_dirichlet,
      ymir_mesh, press_elem,
      domain_options, temp_options, visc_options);
  rhea_stokes_problem_setup_solver (*stokes_problem);

  /* destroy */
  rhea_temperature_destroy (temperature);
  rhea_temperature_destroy (weakzone);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Cleans up Stokes problem and mesh.
 */
static void
collide_setup_clear_all (rhea_stokes_problem_t *stokes_problem,
                         p4est_t *p4est,
                         ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         rhea_discretization_options_t *discr_options)
{
  const char         *this_fn_name = "collide_setup_clear_all";
  ymir_vec_t         *rhs_vel_nonzero_dirichlet = NULL;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* destroy Stokes problem */
  rhs_vel_nonzero_dirichlet =
    rhea_stokes_problem_get_rhs_vel_nonzero_dirichlet (stokes_problem);
  if (rhs_vel_nonzero_dirichlet != NULL) {
    rhea_velocity_destroy (rhs_vel_nonzero_dirichlet);
  }
  rhea_stokes_problem_destroy (stokes_problem);

  /* destroy mesh */
  rhea_discretization_ymir_mesh_destroy (ymir_mesh, press_elem);
  rhea_discretization_p4est_destroy (p4est);

  /* destroy (some) options */
  rhea_discretization_options_clear (discr_options);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Runs Stokes solver.
 */
static void
collide_run_solver (ymir_vec_t *sol_vel_press,
                    ymir_mesh_t *ymir_mesh,
                    ymir_pressure_elem_t *press_elem,
                    rhea_stokes_problem_t *stokes_problem,
                    const int iter_max, const double rel_tol)
{
  const char         *this_fn_name = "collide_run_solver";
  ymir_vec_t         *rhs_vel_nonzero_dirichlet;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* run solver */
  rhea_stokes_problem_solve (sol_vel_press, iter_max, rel_tol, stokes_problem);

  /* add nonzero Dirichlet values of the velocity to the solution */
  rhs_vel_nonzero_dirichlet =
    rhea_stokes_problem_get_rhs_vel_nonzero_dirichlet (stokes_problem);
  if (rhs_vel_nonzero_dirichlet != NULL) {
    ymir_vec_t         *sol_vel = rhea_velocity_new (ymir_mesh);

    ymir_stokes_vec_get_velocity (sol_vel_press, sol_vel, press_elem);
    ymir_vec_add (1.0, rhs_vel_nonzero_dirichlet, sol_vel);
    ymir_stokes_vec_set_velocity (sol_vel, sol_vel_press, press_elem);
    rhea_velocity_destroy (sol_vel);
  }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/* dimensionalize the output*/
void
collide_transform_to_dimensional_stress (ymir_vec_t * stress)
{
  const double        dim_scal = COLLIDE_VISC_REP * COLLIDE_THERM_DIFFUS /
                                 (COLLIDE_EARTH_RADIUS * COLLIDE_EARTH_RADIUS);

  ymir_vec_scale (dim_scal, stress);
}

void
collide_transform_to_dimensional_temperature (ymir_vec_t * temp)
{
  YMIR_ABORT_NOT_REACHED (); //TODO implement this
}

void
collide_transform_to_dimensional_viscosity (ymir_vec_t * visc)
{
  ymir_vec_scale (COLLIDE_VISC_REP, visc);
}

void
collide_transform_to_dimensional_velocity (ymir_vec_t * vel)
{
  const double        dim_scal = COLLIDE_THERM_DIFFUS / COLLIDE_EARTH_RADIUS;

  ymir_vec_scale (dim_scal, vel);
}

void
collide_transform_to_dimensional_strain_rate (ymir_vec_t * strain_rate)
{
  const double        dim_scal = COLLIDE_THERM_DIFFUS /
                                 (COLLIDE_EARTH_RADIUS * COLLIDE_EARTH_RADIUS);

  ymir_vec_scale (dim_scal, strain_rate);
}

void
collide_transform_to_dimensional_pressure (ymir_vec_t * press)
{
  const double        dim_scal = COLLIDE_VISC_REP * COLLIDE_THERM_DIFFUS /
                                 (COLLIDE_EARTH_RADIUS * COLLIDE_EARTH_RADIUS);

  ymir_vec_scale (dim_scal, press);
}


/**
 * Runs the program.
 */
int
main (int argc, char **argv)
{
  const char         *this_fn_name = "collide:main";
  /* MPI */
  MPI_Comm            mpicomm = MPI_COMM_WORLD;
  int                 mpisize, mpirank, ompsize;
  int                 mpiret;
  /* options */
  ymir_options_t     *opt;
  rhea_domain_options_t         domain_options;
  rhea_temperature_options_t    temp_options;
  rhea_viscosity_options_t      visc_options;
  rhea_discretization_options_t discr_options;
  /* collide options */
  int                 vel_dir_bc;
  double              wallslide_vel;
  double              flow_scale;
  double              uwkzone_loc_upper;
  double              uwkzone_loc_lower;
  int                 x_func;
  collide_options_t   collide_options;

  /* options local to this function */
  int                 production_run;
  int                 solver_iter_max;
  double              solver_rel_tol;
  char               *vtk_write_input_path;
  char               *vtk_write_solution_path;
  /* mesh */
  p4est_t            *p4est;
  ymir_mesh_t        *ymir_mesh;
  ymir_pressure_elem_t  *press_elem;
  /* Stokes */
  rhea_stokes_problem_t *stokes_problem;
  ymir_vec_t         *sol_vel_press;

  /*
   * Initialize Libraries
   */

  /* initialize rhea and sub-packages */
  rhea_initialize (argc, argv, mpicomm);

  /* get parallel environment */
  mpiret = MPI_Comm_size (mpicomm, &mpisize); YMIR_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &mpirank); YMIR_CHECK_MPI (mpiret);

#ifdef RHEA_ENABLE_OPENMP
  ompsize = omp_get_max_threads ();
#else
  ompsize = 1;
#endif

  /*
   * Define & Parse Options
   */

  opt = ymir_options_global_new (argv[0] /* program path */);

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  /* basic options */
  YMIR_OPTIONS_CALLBACK, "help", 'h', 0 /* no callback fn args */,
    ymir_options_print_usage_and_exit_fn, NULL /* no arg usage */,
    "Print usage and exit",
  YMIR_OPTIONS_INIFILE, "options-file", 'f',
    ".ini file with option values",

  /* velocity Dirichlet BC's */
  YMIR_OPTIONS_I, "velocity-dirichlet-bc", '\0',
    &vel_dir_bc, COLLIDE_VEL_DIR_BC_WALLSLIDE,
    "Velocity Dirichlet boundary condition",
  YMIR_OPTIONS_D, "wallslide-velocity", '\0',
    &wallslide_vel, 1.0,
    "Tangential velocity",
  YMIR_OPTIONS_D, "flow-scaling", '\0',
    &flow_scale, 1.0,
    "scaling of velocity BC.",
  /* surface location  */
  YMIR_OPTIONS_I, "bound-x-function", '\0',
    &x_func, COLLIDE_X_FUNCTION_IDENTITY,
    "boundary location: surface topography",
 /* weakzone upper and lower bounds*/
  YMIR_OPTIONS_D, "user-weakzone-location-upper",'\0',
    &(uwkzone_loc_upper),
    0.75,
    "user defined weakzone: upper bound",
  YMIR_OPTIONS_D, "user-weakzone-location-lower",'\0',
    &(uwkzone_loc_lower),
    0.25,
    "user defined weakzone: lower bound",

  /* solver options */
  YMIR_OPTIONS_I, "solver-iter-max", '\0',
    &solver_iter_max, 100,
    "Maximum number of iterations for Stokes solver",
  YMIR_OPTIONS_D, "solver-rel-tol", '\0',
    &solver_rel_tol, 1.0e-6,
    "Relative tolerance for Stokes solver",

  /* performance & monitoring options */
  YMIR_OPTIONS_B, "production-run", '\0',
    &(production_run), 0,
    "Execute as a production run (to reduce some overhead and checks)",

  /* vtk output options */
  YMIR_OPTIONS_S, "vtk-write-input-path", '\0',
    &(vtk_write_input_path), NULL,
    "File path for vtk files for the input of the Stokes problem",
  YMIR_OPTIONS_S, "vtk-write-solution-path", '\0',
    &(vtk_write_solution_path), NULL,
    "File path for vtk files for the solution of the Stokes problem",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add sub-options */
  rhea_add_options_all (opt);
  ymir_options_add_suboptions_solver_stokes (opt);

  /* parse options */
  {
    int                 optret;

    optret = ymir_options_parse (SC_LP_INFO, opt, argc, argv);
    if (optret < 0) { /* if parsing was not successful */
      ymir_options_print_usage (SC_LP_INFO, opt, NULL /* args usage */);
      RHEA_GLOBAL_INFO ("Option parsing failed\n");
      exit (0);
    }
  }

  /*
   * Process Collide Options
   */

  collide_options.vel_dir_bc = (collide_vel_dir_bc_t) vel_dir_bc;
  collide_options.wallslide_vel = wallslide_vel;
  collide_options.flow_scale = flow_scale;
  collide_options.uwkzone_loc_upper = uwkzone_loc_upper;
  collide_options.uwkzone_loc_lower = uwkzone_loc_lower;

  collide_options.x_func = (collide_x_func_t) x_func;
  /*
   * Initialize Main Program
   */

  RHEA_GLOBAL_PRODUCTIONF (
      "Into %s (production %i)\n", this_fn_name, production_run);
  RHEA_GLOBAL_PRODUCTIONF (
      "Parallel environment: MPI size %i, OpenMP size %i\n", mpisize, ompsize);
  ymir_set_up (argc, argv, mpicomm, production_run);

  /* print & process options */
  ymir_options_print_summary (SC_LP_INFO, opt);
  rhea_process_options_all (&domain_options, &temp_options,
                            &visc_options, &discr_options);

  /*
   * Setup Mesh
   */

  collide_setup_mesh (&p4est, &ymir_mesh, &press_elem, mpicomm,
                      &domain_options, &discr_options, &collide_options);

  /*
   * Setup Stokes Problem
   */

  collide_setup_stokes (&stokes_problem, ymir_mesh, press_elem,
                        &domain_options, &temp_options, &visc_options,
                        &collide_options, vtk_write_input_path);

  /*
   * Solve Stokes Problem
   */

  /* initialize solution vector */
  sol_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);

  /* run solver */
  collide_run_solver (sol_vel_press, ymir_mesh, press_elem, stokes_problem,
                      solver_iter_max, solver_rel_tol);

  /* write vtk of solution */
  if (vtk_write_solution_path != NULL) {
    ymir_vec_t         *velocity = rhea_velocity_new (ymir_mesh);
    ymir_vec_t         *pressure = rhea_pressure_new (ymir_mesh, press_elem);
    ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);

    ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (
                                                    ymir_mesh->ma->N, ymir_mesh->ma->ompsize);
    ymir_vec_t         *edotII = ymir_dvec_new (ymir_mesh, 1,
                                                     YMIR_GAUSS_NODE);
    ymir_vec_t         *tauII = ymir_dvec_new (ymir_mesh, 1,
                                                     YMIR_GAUSS_NODE);

    ymir_vec_t         *surf_normal_stress = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    ymir_vec_t         *surf_tauII = ymir_face_dvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1,
                                                     YMIR_GAUSS_NODE);
    ymir_vec_t         *vert_n_tau = ymir_dvec_new (ymir_mesh, 1,
                                                     YMIR_GAUSS_NODE);
    ymir_vec_t         *vert_s_tau = ymir_dvec_new (ymir_mesh, 1,
                                                     YMIR_GAUSS_NODE);
    ymir_vec_t         *hori_n_tau = ymir_dvec_new (ymir_mesh, 1,
                                                     YMIR_GAUSS_NODE);
    ymir_vec_t         *hori_s_tau = ymir_dvec_new (ymir_mesh, 1,
                                                     YMIR_GAUSS_NODE);
    ymir_vec_t         *vert_traction = ymir_dvec_new (ymir_mesh, 3,
                                                     YMIR_GAUSS_NODE);
    ymir_vec_t         *hori_traction = ymir_dvec_new (ymir_mesh, 3,
                                                     YMIR_GAUSS_NODE);

    ymir_stokes_vec_get_components (sol_vel_press, velocity, pressure,
                                    press_elem);
    rhea_stokes_problem_get_viscosity (viscosity, stokes_problem);

    rhea_vtk_write_solution (vtk_write_solution_path, velocity, pressure,
                             viscosity);


    /* compute 2nd invariant of the strain rate */
    ymir_second_invariant_vec (velocity, edotII, vel_elem);
    ymir_vec_sqrt (edotII, edotII);

    /* compute 2nd invariant of deviatoric stress tau = 2* (2nd invariant of strain_rate * viscosity )
      and its projection on the surface */
    ymir_vec_copy (edotII,tauII)
    ymir_vec_multiply_in (viscosity, tauII);
    ymir_vec_scale (2.0, tauII);
//    ymir_interp_vec (tauII, surf_tauII);


/* compute surface normal stress sigma */
    collide_physics_compute_normal_boundary_stress (
                   surf_normal_stress, sol_vel_press,
                   rhea_stokes_problem_get_rhs_vel (stokes_problem),
                   rhea_stokes_problem_get_stokes_op (stokes_problem));

/* TODO unsuccessfull, dimensionalization for output */
//
//    collide_transform_to_dimensional_velocity (velocity);
//    collide_transform_to_dimensional_pressure (pressure);
//    collide_transform_to_dimensional_viscosity (viscosity);
//    collide_transform_to_dimensional_stress (surf_normal_stress);


    double vert_n_dir[3] = {0.0, 0.0, 1.0};
    collide_normal_stress (velocity, vert_n_tau, vert_s_tau, vert_traction, vert_n_dir,
                            viscosity, vel_elem);
    double hori_n_dir[3] = {0.0, 1.0, 0.0};
    collide_normal_stress (velocity, hori_n_tau, hori_s_tau, hori_traction, hori_n_dir,
                            viscosity, vel_elem);

    if (edotII != NULL) {
      char            path[BUFSIZ];

      snprintf (path, BUFSIZ, "%s_stress", vtk_write_solution_path);
      ymir_vtk_write (ymir_mesh, path,
                      edotII, "edotII",
                      tauII, "tauII",
                      vert_traction, "vertical tau_vec",
                      vert_n_tau, "vertical normal tau",
                      vert_s_tau, "vertical shear tau",
                      hori_traction, "horizontal tau_vec",
                      hori_n_tau, "horizontal normal tau",
                      hori_s_tau, "horizontal shear tau",
//                      surf_tauII, "surf_tauII",
                      surf_normal_stress, "surf_normal_stress",
                      NULL);
    }

    /* destroy */
    ymir_vec_destroy (edotII);
    ymir_vec_destroy (tauII);
    ymir_vec_destroy (vert_traction);
    ymir_vec_destroy (vert_n_tau);
    ymir_vec_destroy (vert_s_tau);
    ymir_vec_destroy (hori_traction);
    ymir_vec_destroy (hori_n_tau);
    ymir_vec_destroy (hori_s_tau);
    ymir_vec_destroy (surf_tauII);
    ymir_vec_destroy (surf_normal_stress);
    ymir_velocity_elem_destroy (vel_elem);
    rhea_velocity_destroy (velocity);
    rhea_pressure_destroy (pressure);
    rhea_viscosity_destroy (viscosity);
  }

  /* destroy */
  rhea_velocity_pressure_destroy (sol_vel_press);

  /*
   * Finalize
   */

  /* destroy Stokes problem and mesh */
  collide_setup_clear_all (stokes_problem, p4est, ymir_mesh, press_elem,
                           &discr_options);

  /* destroy options */
  ymir_options_global_destroy ();

  /* print that this function is ending */
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

  /* finalize rhea */
  rhea_finalize ();

  return 0;
}
