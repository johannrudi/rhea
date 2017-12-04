#include <rhea.h>
#include <ymir_velocity_vec.h>
#include <ymir_derivative_elem.h>
#include <ymir_vec_ops.h>
#include <ymir_stokes_pc.h>
#include <ymir_stokes_vec.h>
#include <ymir_stress_op.h>
#include <ymir_pressure_vec.h>

/* weak zone 2plates_poly2 */
#define SLABS_WEAKZONE_2PLATES_SUBDU_LONGITUDE (-100.0)      /* rhea1:  0.13 */
#define SLABS_WEAKZONE_2PLATES_SUBDU_DIP_ANGLE (-1.0)        /* rhea1:   N/A */
#define SLABS_WEAKZONE_2PLATES_SUBDU_DEPTH (80.0e3)          /* rhea1: 50 km */
#define SLABS_WEAKZONE_2PLATES_SUBDU_WIDTH (-1.0)            /* rhea1:   N/A */
#define SLABS_WEAKZONE_2PLATES_SUBDU_THICKNESS (20.0e3)      /* rhea1: 20 km */
#define SLABS_WEAKZONE_2PLATES_SUBDU_THICKNESS_CONST (5.0e3) /* rhea1: 10 km */
#define SLABS_WEAKZONE_2PLATES_SUBDU_WEAK_FACTOR (1.0e-5)    /* rhea1:  1e-5 */
#define SLABS_WEAKZONE_2PLATES_RIDGE_DEPTH (30.0e3)          /* rhea1: 30 km */
#define SLABS_WEAKZONE_2PLATES_RIDGE_WIDTH (30.0e3)          /* rhea1: 30 km */
#define SLABS_WEAKZONE_2PLATES_RIDGE_SMOOTHWIDTH (10.0e3)    /* rhea1:  5 km */
#define SLABS_WEAKZONE_2PLATES_RIDGE_WEAK_FACTOR (1.0e-5)    /* rhea1:  1e-5 */


/* compute traction as well as normal/shear stress at each element*/
void
slabs_postp_weakzone_traction_elem (sc_dmatrix_t * in, sc_dmatrix_t * out,
                                    sc_dmatrix_t * visc, sc_dmatrix_t * svisc,
                                    sc_dmatrix_t * TItens, ymir_velocity_elem_t * vel_elem,
                                    double *_sc_restrict rxd, double *_sc_restrict sxd,
                                    double *_sc_restrict txd, double *_sc_restrict ryd,
                                    double *_sc_restrict syd, double *_sc_restrict tyd,
                                    double *_sc_restrict rzd, double *_sc_restrict szd,
                                    double *_sc_restrict tzd, sc_dmatrix_t * drst, sc_dmatrix_t * brst)
{
  const int           N = ymir_n (vel_elem->N);
  const int           Np = ymir_n (vel_elem->Np);
  int                 gp, i, j, k, l;
  double              S[9], E[9];
  double              cs, ct;
  double             *_sc_restrict viscd  = visc->e[0];
  double             *_sc_restrict sviscd = svisc->e[0];
  double             *_sc_restrict outd   = out->e[0];

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
    double               outv = .0;
    double               *_sc_restrict temptensd = temptens->e[0] + 9 * gp;
    double               *_sc_restrict T   = TItens->e[0] + 9 * gp;

    S[0] = temptensd[0];
    S[1] = (temptensd[1] + temptensd[3]) * (1. / 2.);
    S[2] = (temptensd[2] + temptensd[6]) * (1. / 2.);
    S[3] = S[1];
    S[4] = temptensd[4];
    S[5] = (temptensd[5] + temptensd[7]) * (1. / 2.);
    S[6] = S[2];
    S[7] = S[5];
    S[8] = temptensd[8];

    for (k = 0; k < 9; k++) {
      E[k] = 0.0;
    }
    for (k=0; k<3; k++) {
      for (l=0; l<3; l++) {
        E[0] += T[    k] * S[3*k + l] * T[3*l    ];
        E[1] += T[    k] * S[3*k + l] * T[3*l + 1];
        E[2] += T[    k] * S[3*k + l] * T[3*l + 2];
        E[4] += T[3 + k] * S[3*k + l] * T[3*l + 1];
        E[5] += T[3 + k] * S[3*k + l] * T[3*l + 2];
        E[8] += T[6 + k] * S[3*k + l] * T[3*l + 2];
      }
    }
    E[3] = E[1];
    E[6] = E[2];
    E[7] = E[5];

    cs = viscd[gp] + sviscd[gp];
    ct = viscd[gp] - sviscd[gp];

    /* compute linear combination of viscous stress and 3x3 tensor */
    S[0] = cs * S[0] + ct * E[0];
    S[1] = cs * S[1] + ct * E[1];
    S[2] = cs * S[2] + ct * E[2];
    S[3] = cs * S[3] + ct * E[3];
    S[4] = cs * S[4] + ct * E[4];
    S[5] = cs * S[5] + ct * E[5];
    S[6] = cs * S[6] + ct * E[6];
    S[7] = cs * S[7] + ct * E[7];
    S[8] = cs * S[8] + ct * E[8];

    for (k = 0; k < 9; k++) {
      outv += SC_SQR (S[k]);
    }
    outd[gp] = sqrt (0.5 * outv);
  }
}

/* compute traction as well as normal and shear stress*/
void
slabs_postp_weakzone_traction (ymir_cvec_t * vel, ymir_dvec_t * traction,
                               ymir_dvec_t * visc, ymir_dvec_t *svisc,
                               ymir_dvec_t * TItens, ymir_velocity_elem_t * vel_elem)
{
  ymir_mesh_t         *mesh = vel->mesh;
  ymir_locidx_t       elid;
  const int           N  = ymir_n (mesh->cnodes->N);
  const int           Np = ymir_np (mesh->cnodes->N);
  const int           K  = mesh->cnodes->K;
  sc_dmatrix_t       *elemin = sc_dmatrix_new (1, 3 * Np);
  sc_dmatrix_t       *elemout = sc_dmatrix_new (1, Np); // TODO
  sc_dmatrix_t       *elemvisc = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemsvisc = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemTItens = sc_dmatrix_new (1, 9 * Np);

  ymir_dvec_set_zero (traction);

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
    ymir_dvec_get_elem_interp (svisc, elemsvisc, YMIR_STRIDE_COMP, elid,
                              YMIR_GAUSS_NODE, YMIR_READ);
    ymir_dvec_get_elem_interp (TItens, elemTItens, YMIR_STRIDE_COMP, elid,
                              YMIR_GAUSS_NODE, YMIR_READ);


    slabs_postp_weakzone_traction_elem (elemin, elemout,
                                        elemvisc, elemsvisc, elemTItens,
                                        vel_elem, rxd, sxd, txd, ryd,
                                        syd, tyd, rzd, szd, tzd, mesh->drst, mesh->brst);

    ymir_dvec_set_elem_interp (traction, elemout, YMIR_STRIDE_COMP, elid,
                               YMIR_GAUSS_NODE, YMIR_SET);

    ymir_read_view_release (elemvisc);
    ymir_read_view_release (elemsvisc);
    ymir_read_view_release (elemTItens);
  }

  sc_dmatrix_destroy (elemin);
  sc_dmatrix_destroy (elemout);
  sc_dmatrix_destroy (elemvisc);
  sc_dmatrix_destroy (elemsvisc);
  sc_dmatrix_destroy (elemTItens);
}


