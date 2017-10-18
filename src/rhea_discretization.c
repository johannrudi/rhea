/*
 */

#include <rhea_discretization.h>
#include <rhea_base.h>
#include <rhea_amr.h>
#include <p8est_extended.h>
#include <mangll_p8est.h>
#include <ymir_stress_pc.h>
#include <ymir_gmg.h>

/* default options */
#define RHEA_DISCRETIZATION_DEFAULT_ORDER (2)
#define RHEA_DISCRETIZATION_DEFAULT_LEVEL_MIN (1)
#define RHEA_DISCRETIZATION_DEFAULT_LEVEL_MAX (20)
#define RHEA_DISCRETIZATION_DEFAULT_REFINEMENT_TYPE "uniform"

int                 rhea_discretization_order =
  RHEA_DISCRETIZATION_DEFAULT_ORDER;
int                 rhea_discretization_level_min =
  RHEA_DISCRETIZATION_DEFAULT_LEVEL_MIN;
int                 rhea_discretization_level_max =
  RHEA_DISCRETIZATION_DEFAULT_LEVEL_MAX;
char               *rhea_discretization_refinement_type =
  RHEA_DISCRETIZATION_DEFAULT_REFINEMENT_TYPE;

void
rhea_discretization_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Discretization";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_I, "order", '\0',
    &(rhea_discretization_order), RHEA_DISCRETIZATION_DEFAULT_ORDER,
    "Order of finite element discretization",

  YMIR_OPTIONS_I, "level-min", '\0',
    &(rhea_discretization_level_min), RHEA_DISCRETIZATION_DEFAULT_LEVEL_MIN,
    "Minumum level of mesh refinement",
  YMIR_OPTIONS_I, "level-max", '\0',
    &(rhea_discretization_level_max), RHEA_DISCRETIZATION_DEFAULT_LEVEL_MAX,
    "Maximum level of mesh refinement",

  YMIR_OPTIONS_S, "refinement-type", '\0',
    &(rhea_discretization_refinement_type),
    RHEA_DISCRETIZATION_DEFAULT_REFINEMENT_TYPE,
    "Init refinement type: uniform, half",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);
}

/**
 * Identity geometry transformation.
 */
static void
rhea_discretization_X_fn_identity (mangll_tag_t tag, mangll_locidx_t np,
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

#define RHEA_DISCR_SHELL_X_TAG_TOP 3
#define RHEA_DISCR_SHELL_X_TAG_LEFT 2
#define RHEA_DISCR_SHELL_X_TAG_RIGHT 0
#define RHEA_DISCR_SHELL_X_TAG_BOTTOM 1
#define RHEA_DISCR_SHELL_X_TAG_BACK 4
#define RHEA_DISCR_SHELL_X_TAG_FRONT 5

/**
 * Geometry transformation for 24-octree hollow sphere.
 */
static void
rhea_discretization_X_fn_shell (mangll_tag_t tag, mangll_locidx_t np,
                                const double *_sc_restrict EX,
                                const double *_sc_restrict EY,
                                const double *_sc_restrict EZ,
                                double *_sc_restrict X,
                                double *_sc_restrict Y,
                                double *_sc_restrict Z, void *data)
{
  rhea_domain_options_t  *domain_options = data;
  const double        R1 = domain_options->radius_min;
  const double        R2 = domain_options->radius_max;
  const double        R2byR1 = R2 / R1;
  const double        R1sqrbyR2 = R1 * R1 / R2;

  mangll_locidx_t     il;
  double              x, y, R, q;

  RHEA_ASSERT (0 <= tag && tag < 24);
#ifdef RHEA_ENABLE_DEBUG
  for (il = 0; il < np; ++il) {
    RHEA_ASSERT (EX[il] < 1.0 + SC_1000_EPS && EX[il] > -1.0 - SC_1000_EPS);
    RHEA_ASSERT (EY[il] < 1.0 + SC_1000_EPS && EY[il] > -1.0 - SC_1000_EPS);
    RHEA_ASSERT (EZ[il] < 2.0 + SC_1000_EPS && EZ[il] > +1.0 - SC_1000_EPS);
  }
#endif

  switch (tag / 4) {
  case RHEA_DISCR_SHELL_X_TAG_TOP:  /* top (+z) */
    for (il = 0; il < np; ++il) {
      /* transform abc[0] and y in-place for nicer grading */
      x = tan (EX[il] * M_PI_4);
      y = tan (EY[il] * M_PI_4);

      /* compute transformation ingredients */
      R = R1sqrbyR2 * pow (R2byR1, EZ[il]);
      q = R / sqrt (x * x + y * y + 1.);

      /* assign correct coordinates based on patch id */
      X[il] = +q * y;
      Y[il] = -q * x;
      Z[il] = +q;
    }
    break;
  case RHEA_DISCR_SHELL_X_TAG_LEFT:  /* left (-y) */
    for (il = 0; il < np; ++il) {
      x = tan (EX[il] * M_PI_4);
      y = tan (EY[il] * M_PI_4);
      R = R1sqrbyR2 * pow (R2byR1, EZ[il]);
      q = R / sqrt (x * x + y * y + 1.);
      X[il] = -q;
      Y[il] = -q * x;
      Z[il] = +q * y;
    }
    break;
  case RHEA_DISCR_SHELL_X_TAG_BOTTOM:  /* bottom (-z) */
    for (il = 0; il < np; ++il) {
      x = tan (EX[il] * M_PI_4);
      y = tan (EY[il] * M_PI_4);
      R = R1sqrbyR2 * pow (R2byR1, EZ[il]);
      q = R / sqrt (x * x + y * y + 1.);
      X[il] = -q * y;
      Y[il] = -q * x;
      Z[il] = -q;
    }
    break;
  case RHEA_DISCR_SHELL_X_TAG_RIGHT:  /* right (+y) */
    for (il = 0; il < np; ++il) {
      x = tan (EX[il] * M_PI_4);
      y = tan (EY[il] * M_PI_4);
      R = R1sqrbyR2 * pow (R2byR1, EZ[il]);
      q = R / sqrt (x * x + y * y + 1.);
      X[il] = +q;
      Y[il] = -q * x;
      Z[il] = -q * y;
    }
    break;
  case RHEA_DISCR_SHELL_X_TAG_BACK:  /* back (-x) */
    for (il = 0; il < np; ++il) {
      x = tan (EX[il] * M_PI_4);
      y = tan (EY[il] * M_PI_4);
      R = R1sqrbyR2 * pow (R2byR1, EZ[il]);
      q = R / sqrt (x * x + y * y + 1.);
      X[il] = -q * x;
      Y[il] = +q;
      Z[il] = +q * y;
    }
    break;
  case RHEA_DISCR_SHELL_X_TAG_FRONT:  /* front (+x) */
    for (il = 0; il < np; ++il) {
      x = tan (EX[il] * M_PI_4);
      y = tan (EY[il] * M_PI_4);
      R = R1sqrbyR2 * pow (R2byR1, EZ[il]);
      q = R / sqrt (x * x + y * y + 1.);
      X[il] = +q * x;
      Y[il] = -q;
      Z[il] = +q * y;
    }
    break;
  default:
    RHEA_ABORT_NOT_REACHED ();
  }
}

/**
 * Gets tag corresponding to a part of the spherical shell domain from element
 * coordinates:
 *   If (EX,EY) \in [-1,+1) x [-1,+1), then return `top` tag.
 *   If (EX,EY) \in [-3,-1) x [-1,+1), then return `back` tag.
 *   If (EX,EY) \in [+1,+3) x [-1,+1), then return `front` tag.
 *   If (EX,EY) \in [-1,+1) x [-3,-1), then return `left` tag.
 *   If (EX,EY) \in [-1,+1) x [+1,+3), then return `right` tag.
 */
static mangll_tag_t
rhea_discretization_X_fn_box_spherical_get_tag (const double EX,
                                                const double EY)
{
  const mangll_tag_t  n_tags_per_part = 4;

  if ( (-1.0 <= EX && EX < 1.0) && (-1.0 <= EY && EY < 1.0) ) {
    /* return `top` if (EX,EY) \in [-1,+1] x [-1,+1] */
    return RHEA_DISCR_SHELL_X_TAG_TOP * n_tags_per_part;
  }
  else if ( (-3.0 <= EX && EX < -1.0) && (-1.0 <= EY && EY < 1.0) ) {
    /* return `back` if (EX,EY) \in [-3,-1) x [-1,+1] */
    return RHEA_DISCR_SHELL_X_TAG_BACK * n_tags_per_part;
  }
  else if ( (1.0 <= EX && EX < 3.0) && (-1.0 <= EY && EY < 1.0) ) {
    /* return `front` if (EX,EY) \in (+1,+3] x [-1,+1] */
    return RHEA_DISCR_SHELL_X_TAG_FRONT * n_tags_per_part;
  }
  else if ( (-1.0 <= EX && EX < 1.0) && (-3.0 <= EY && EY < -1.0) ) {
    /* return `left` if (EX,EY) \in [-1,+1] x [-3,-1) */
    return RHEA_DISCR_SHELL_X_TAG_LEFT * n_tags_per_part;
  }
  else if ( (-1.0 <= EX && EX < 1.0) && (1.0 <= EY && EY < 3.0) ) {
    /* return `right` if (EX,EY) \in [-1,+1] x (+1,+3] */
    return RHEA_DISCR_SHELL_X_TAG_RIGHT * n_tags_per_part;
  }
  else { /* unknown part of spherical shell domain */
    RHEA_ABORT_NOT_REACHED ();
  }
}

/**
 * Geometry transformation for chunk and slice of hollow sphere.
 */
static void
rhea_discretization_X_fn_box_spherical (mangll_tag_t tag, mangll_locidx_t np,
                                        const double *_sc_restrict EX,
                                        const double *_sc_restrict EY,
                                        const double *_sc_restrict EZ,
                                        double *_sc_restrict X,
                                        double *_sc_restrict Y,
                                        double *_sc_restrict Z, void *data)
{
  const mangll_tag_t  n_tags_per_part = 4;
  mangll_tag_t        shell_tag;
  double             *E_replace;
  mangll_locidx_t     il;

  /* set tag for shell X-function */
  shell_tag = rhea_discretization_X_fn_box_spherical_get_tag (EX[0], EY[0]);

  /* init replacement for a component of element coordinates */
  if (shell_tag/n_tags_per_part != RHEA_DISCR_SHELL_X_TAG_TOP) {
    E_replace = RHEA_ALLOC (double, np);
  }

  /* call shell X-function */
  switch (shell_tag/n_tags_per_part) {
  case RHEA_DISCR_SHELL_X_TAG_TOP:  /* top */
    rhea_discretization_X_fn_shell (shell_tag, np, EX, EY, EZ, X, Y, Z, data);
    break;
  case RHEA_DISCR_SHELL_X_TAG_BACK:  /* back */
    for (il = 0; il < np; ++il) { /* shift EX from [-3,-1] to [-1,+1] */
      E_replace[il] = EX[il] + 2.0;
    }
    rhea_discretization_X_fn_shell (shell_tag, np, E_replace, EY, EZ, X, Y, Z,
                                    data);
    break;
  case RHEA_DISCR_SHELL_X_TAG_FRONT:  /* front */
    for (il = 0; il < np; ++il) { /* shift EX from [+1,+3] to [-1,+1] */
      E_replace[il] = EX[il] - 2.0;
    }
    rhea_discretization_X_fn_shell (shell_tag, np, E_replace, EY, EZ, X, Y, Z,
                                    data);
    break;
  case RHEA_DISCR_SHELL_X_TAG_LEFT:  /* left */
    for (il = 0; il < np; ++il) { /* shift EY from [-3,-1] to [-1,+1] */
      E_replace[il] = EY[il] + 2.0;
    }
    rhea_discretization_X_fn_shell (shell_tag, np, EX, E_replace, EZ, X, Y, Z,
                                    data);
    break;
  case RHEA_DISCR_SHELL_X_TAG_RIGHT:  /* right */
    for (il = 0; il < np; ++il) { /* shift EY from [+1,+3] to [-1,+1] */
      E_replace[il] = EY[il] - 2.0;
    }
    rhea_discretization_X_fn_shell (shell_tag, np, EX, E_replace, EZ, X, Y, Z,
                                    data);
    break;
  default: /* unknown part of spherical shell domain */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* destroy */
  if (shell_tag/n_tags_per_part != RHEA_DISCR_SHELL_X_TAG_TOP) {
    RHEA_FREE (E_replace);
  }
}

void
rhea_discretization_process_options (rhea_discretization_options_t *opt,
                                     rhea_domain_options_t *domain_options)
{
  /* set discretization order and mesh refinement levels */
  opt->order = rhea_discretization_order;
  opt->level_min = rhea_discretization_level_min;
  opt->level_max = rhea_discretization_level_max;

  /* set undefined boundary information */
  opt->boundary = NULL;

  /* set reference-to-physical transformation function */
  switch (domain_options->shape) {
  case RHEA_DOMAIN_CUBE:
  case RHEA_DOMAIN_BOX:
    opt->X_fn = rhea_discretization_X_fn_identity;
    opt->X_data = NULL;
    break;
  case RHEA_DOMAIN_SHELL:
    opt->X_fn = rhea_discretization_X_fn_shell;
    opt->X_data = domain_options;
    break;
  case RHEA_DOMAIN_CUBE_SPHERICAL:
  case RHEA_DOMAIN_BOX_SPHERICAL:
    opt->X_fn = rhea_discretization_X_fn_box_spherical;
    opt->X_data = domain_options;
    break;
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }
}

void
rhea_discretization_set_user_X_fn (rhea_discretization_options_t *opt,
                                   mangll_X_t X_fn, void *X_data)
{
  opt->X_fn = X_fn;
  opt->X_data = X_data;
}

void
rhea_discretization_options_set_boundary (rhea_discretization_options_t *opt,
                                          p4est_t *p4est,
                                          rhea_domain_options_t *domain_options)
{
  opt->boundary = rhea_domain_boundary_new (p4est, domain_options);

  /* initialize geometric multigrid */
  ymir_gmg_init_mesh (p4est, opt->boundary->e_to_fm_fn,
                      opt->boundary->tree_to_bf);
  ymir_gmg_init_bc (rhea_domain_create_velocity_dirichlet_bc, domain_options,
                    NULL /* Robin BC's */, NULL /* Robin data */);
}

void
rhea_discretization_options_clear (rhea_discretization_options_t *opt)
{
  rhea_domain_boundary_destroy (opt->boundary);
  opt->boundary = NULL;
}

p4est_t *
rhea_discretization_p4est_new (MPI_Comm mpicomm,
                               rhea_discretization_options_t *opt,
                               rhea_domain_options_t *domain_options)
{
  /* arguments for creating a new p4est object */
  const int           n_quadrants_init = 0;
  int                 level_init = opt->level_min;
  const int           fill_uniformly = 1;
  size_t              data_size = 0;
  p4est_init_t        init_fn = NULL;
  /* arguments for first refinement */
  const char         *refine = rhea_discretization_refinement_type;
  p4est_refine_t      refine_fn = NULL;
  void               *refine_data = NULL;
  const int           refine_recursively = 1;
  const int           refine_level_max = opt->level_min;
  const int           partition_for_coarsening = 1;
  p4est_weight_t      partition_weight_fn = NULL;
  p4est_replace_t     replace_fn = NULL;
  /* p4est objects */
  p4est_connectivity_t *conn;
  p4est_t            *p4est;

  /*
   * Create Connectivity
   */

  switch (domain_options->shape) {
  case RHEA_DOMAIN_CUBE:
    conn = p8est_connectivity_new_unitcube ();
    break;

  case RHEA_DOMAIN_BOX:
    RHEA_ASSERT (isfinite (domain_options->x_min));
    RHEA_ASSERT (isfinite (domain_options->y_min));
    RHEA_ASSERT (isfinite (domain_options->z_min));
    RHEA_ASSERT (isfinite (domain_options->x_max));
    RHEA_ASSERT (isfinite (domain_options->y_max));
    RHEA_ASSERT (isfinite (domain_options->z_max));
    RHEA_ASSERT (domain_options->x_min < domain_options->x_max);
    RHEA_ASSERT (domain_options->y_min < domain_options->y_max);
    RHEA_ASSERT (domain_options->z_min < domain_options->z_max);
    RHEA_ASSERT (0 < domain_options->box_subdivision_x);
    RHEA_ASSERT (0 < domain_options->box_subdivision_y);
    RHEA_ASSERT (0 < domain_options->box_subdivision_z);

    conn = p8est_connectivity_new_brick (domain_options->box_subdivision_x,
                                         domain_options->box_subdivision_y,
                                         domain_options->box_subdivision_z,
                                         0, 0, 0);
    /* scale & shift reference coordinates of vertices such that
     *   x component: [0, subdivision_x] -> [x_min, x_max]
     *   y component: [0, subdivision_y] -> [y_min, y_max]
     *   z component: [0, subdivision_z] -> [z_min, z_max]
     */
    {
      const double        xmin = domain_options->x_min;
      const double        ymin = domain_options->y_min;
      const double        zmin = domain_options->z_min;
      const double        xmax = domain_options->x_max;
      const double        ymax = domain_options->y_max;
      const double        zmax = domain_options->z_max;
      const double        subx = (double) domain_options->box_subdivision_x;
      const double        suby = (double) domain_options->box_subdivision_y;
      const double        subz = (double) domain_options->box_subdivision_z;
      double             *vertices = conn->vertices;
      p4est_topidx_t      vi;

      for (vi = 0; vi < conn->num_vertices; vi++) { /* loop over all vertices */
        vertices[3*vi    ] = (xmax - xmin) * vertices[3*vi    ] / subx + xmin;
        vertices[3*vi + 1] = (ymax - ymin) * vertices[3*vi + 1] / suby + ymin;
        vertices[3*vi + 2] = (zmax - zmin) * vertices[3*vi + 2] / subz + zmin;
      }
    }
    break;

  case RHEA_DOMAIN_SHELL:
    conn = p8est_connectivity_new_shell ();
    break;

  case RHEA_DOMAIN_CUBE_SPHERICAL:
    conn = p8est_connectivity_new_unitcube ();
    /* shift coordinates of vertices, so (modified) shell_X fnc can be used:
     *   x component: [0, 1] -> [-0.5, 0.5]
     *   y component: [0, 1] -> [-0.5, 0.5]
     *   z component: [0, 1] -> [ 1  , 2  ]
     */
    {
      double         *vertices = conn->vertices;
      p4est_topidx_t  vi;

      for (vi = 0; vi < conn->num_vertices; vi++) { /* loop over all vertices */
        vertices[3*vi    ] -= 0.5;
        vertices[3*vi + 1] -= 0.5;
        vertices[3*vi + 2] += 1.0;
      }
    }
    break;

  case RHEA_DOMAIN_BOX_SPHERICAL:
    RHEA_ASSERT (0 < domain_options->box_subdivision_x);
    RHEA_ASSERT (0 < domain_options->box_subdivision_y);
    RHEA_ASSERT (0 < domain_options->box_subdivision_z);

    conn = p8est_connectivity_new_brick (domain_options->box_subdivision_x,
                                         domain_options->box_subdivision_y,
                                         domain_options->box_subdivision_z,
                                         0, 0, 0);
    /* shift coordinates of vertices, so (modified) shell_X fnc can be used:
     *   x component: [0, subdiv_x] -> 1/subdiv_z * [-subdiv_x/2, subdiv_x/2]
     *   y component: [0, subdiv_y] -> 1/subdiv_z * [-subdiv_y/2, subdiv_y/2]
     *   z component: [0, subdiv_z] ->              [ 1         , 2         ]
     */
    {
      const double        subx = (double) domain_options->box_subdivision_x;
      const double        suby = (double) domain_options->box_subdivision_y;
      const double        subz = (double) domain_options->box_subdivision_z;
      double             *vertices = conn->vertices;
      p4est_topidx_t      vi;

      for (vi = 0; vi < conn->num_vertices; vi++) { /* loop over all vertices */
        vertices[3*vi    ] = (vertices[3*vi    ] - 0.5 * subx) / subz;
        vertices[3*vi + 1] = (vertices[3*vi + 1] - 0.5 * suby) / subz;
        vertices[3*vi + 2] = (vertices[3*vi + 2] / subz) + 1.0;
      }
    }
    break;

  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }

  /*
   * Set Up Initial Refinement
   */

  /* set refine function; reduce min level for new p4est to allow additional
   * refinement */
  if (refine != NULL) {
    if (!strcmp (refine, "uniform")) {
      /* uniform refinement is the default action */
    }
    else if (!strcmp (refine, "half")) {
      refine_fn = rhea_amr_p4est_refine_half;
      level_init = level_init - 1;
    }
  //else if (!strcmp (refine, "evenodd")) {
  //  refine_fn = mangll_p8est_refine_evenodd;
  //  level_init = level_init - 1;
  //}
  //else if (!strcmp (refine, "oddeven")) {
  //  refine_fn = mangll_p8est_refine_oddeven;
  //  level_init = level_init - 1;
  //}
  //else if (!strcmp (refine, "origin")) {
  //  refine_fn = mangll_p8est_refine_origin;
  //  level_init = 0;
  //}
  //else if (!strcmp (refine, "fractal")) {
  //  refine_fn = mangll_p8est_refine_fractal;
  //  level_init = level_init - 4;
  //}
  //else if (!strcmp (refine, "finercenter2")) {
  //  refine_fn = mangll_p8est_refine_finercenter2;
  //  level_init = level_init - 1;
  //}
  //else if (!strcmp (refine, "finercenter4")) {
  //  refine_fn = mangll_p8est_refine_finercenter4;
  //  level_init = level_init - 1;
  //}
    else {
      RHEA_GLOBAL_LERROR ("Unknown refinement type");
      return NULL;
    }
  }

  /* restrict min level */
  level_init = SC_MAX (0, level_init);
  if (opt->level_min <= opt->level_max) {
    level_init = SC_MIN (level_init, opt->level_max);
  }

  /*
   * Create p4est
   */

  /* create new p4est */
  p4est = p4est_new_ext (mpicomm, conn, n_quadrants_init, level_init,
                         fill_uniformly, data_size, init_fn, refine_data);

  /* refine */
  if (refine_fn != NULL) {
    p4est_refine_ext (p4est, refine_recursively, refine_level_max, refine_fn,
                      init_fn, replace_fn);
    p4est_partition_ext (p4est, partition_for_coarsening, partition_weight_fn);
    p4est_balance (p4est, P4EST_CONNECT_FULL, init_fn);
    p4est_partition_ext (p4est, partition_for_coarsening, partition_weight_fn);
  }

  /* initialize user data */
  //slabs_discr_p8est_init_data (p4est, discr_options->inspect_p4est);
  //TODO

  /* return p4est */
  return p4est;
}

void
rhea_discretization_p4est_destroy (p4est_t *p4est)
{
  p4est_connectivity_t *conn = p4est->connectivity;

  //slabs_discr_p8est_clear_data (p4est);
  //TODO
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);
}

void
rhea_discretization_mangll_and_cnodes_new (mangll_t **mangll,
                                           mangll_cnodes_t **cnodes,
                                           p4est_t *p4est,
                                           rhea_discretization_options_t *opt)
{
  MPI_Comm            mpicomm = p4est->mpicomm;
  const int           order = opt->order;
  const mangll_refel_quadrature_type_t  quad_type = MANGLL_REFEL_QUAD_GAUSS;

  p4est_ghost_t      *ghost;
  mangll_mesh_t      *mangll_mesh;

  /* create p4est ghost */
  if (cnodes != NULL) {
    ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
  }
  else {
    ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  }

  /* create mangll mesh structure */
  mangll_mesh = mangll_p8est_mesh_new_full (p4est, ghost);

  /* assign mapping to physical space */
  mangll_mesh->X_fn = opt->X_fn;
  mangll_mesh->X_data = opt->X_data;

  RHEA_ASSERT (mangll != NULL);
  if (cnodes != NULL) { /* if mesh with continuous nodes should be created */
    /* create continuous node structure */
    *cnodes = mangll_p8est_cnodes_new (p4est, ghost, order);

    /* create mangll structure with continuous geometry */
    *mangll = mangll_new_ext (mpicomm, order, quad_type, mangll_mesh, *cnodes);
  }
  else {
    /* create mangll structure with discontinuous geometry */
    *mangll = mangll_new_ext (mpicomm, order, quad_type, mangll_mesh, NULL);
  }

  /* destroy ghost */
  p4est_ghost_destroy (ghost);
}

void
rhea_discretization_ymir_mesh_new_from_mangll (
                                          ymir_mesh_t **ymir_mesh,
                                          ymir_pressure_elem_t **press_elem,
                                          mangll_t *mangll,
                                          mangll_cnodes_t *cnodes,
                                          rhea_discretization_options_t *opt)
{
  const int           skip_face_prealloc = 1;
  const int           skip_diag_prealloc = ymir_stress_pc_gmg;
  rhea_domain_boundary_t  *boundary = opt->boundary;

  /* check_input */
  RHEA_ASSERT (opt->boundary != NULL);

  /* create ymir mesh */
  *ymir_mesh = ymir_mesh_new_ext (mangll, cnodes, boundary->e_to_fm_fn,
                                  boundary->tree_to_bf, skip_face_prealloc,
                                  skip_diag_prealloc);

  /* create pressure element */
  if (press_elem != NULL) {
    *press_elem = ymir_pressure_elem_new (mangll->refel, mangll->ompsize);
  }
}

void
rhea_discretization_ymir_mesh_new_from_p4est (
                                          ymir_mesh_t **ymir_mesh,
                                          ymir_pressure_elem_t **press_elem,
                                          p4est_t *p4est,
                                          rhea_discretization_options_t *opt)
{
  mangll_t           *mangll;
  mangll_cnodes_t    *cnodes;

  /* create mangll & cnodes */
  rhea_discretization_mangll_and_cnodes_new (&mangll, &cnodes, p4est, opt);

  /* create ymir mesh & pressure element */
  rhea_discretization_ymir_mesh_new_from_mangll (ymir_mesh, press_elem,
                                                 mangll, cnodes, opt);
}

/**
 * Destroys ymir mesh and corresponding mangll structures.
 */
void
rhea_discretization_ymir_mesh_destroy (ymir_mesh_t *ymir_mesh,
                                       ymir_pressure_elem_t *press_elem)
{
  /* destroy pressure element */
  if (press_elem != NULL) {
    ymir_pressure_elem_destroy (press_elem);
  }

  /* destroy mangll, cnodes, and ymir_mesh */
  if (ymir_mesh != NULL) {
    mangll_destroy (ymir_mesh->ma);
    mangll_p8est_cnodes_destroy (ymir_mesh->cnodes);
    ymir_mesh_destroy (ymir_mesh);
  }
}
