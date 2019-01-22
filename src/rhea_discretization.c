#include <rhea_discretization.h>
#include <rhea_base.h>
#include <rhea_amr.h>
#include <rhea_io_mpi.h>
#include <p8est_extended.h>
#include <mangll_p8est.h>
#include <ymir_stress_pc.h>
#include <ymir_gmg.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_DISCRETIZATION_DEFAULT_ORDER (2)
#define RHEA_DISCRETIZATION_DEFAULT_LEVEL_MIN (1)
#define RHEA_DISCRETIZATION_DEFAULT_LEVEL_MAX (18)
#define RHEA_DISCRETIZATION_DEFAULT_P4EST_FILE_PATH NULL

int                 rhea_discretization_order =
  RHEA_DISCRETIZATION_DEFAULT_ORDER;
int                 rhea_discretization_level_min =
  RHEA_DISCRETIZATION_DEFAULT_LEVEL_MIN;
int                 rhea_discretization_level_max =
  RHEA_DISCRETIZATION_DEFAULT_LEVEL_MAX;
char               *rhea_discretization_p4est_file_path =
  RHEA_DISCRETIZATION_DEFAULT_P4EST_FILE_PATH;

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

  YMIR_OPTIONS_S, "p4est-file-path", '\0',
    &(rhea_discretization_p4est_file_path),
    RHEA_DISCRETIZATION_DEFAULT_P4EST_FILE_PATH,
    "Read p4est from disk at this file path",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add sub-options */
  rhea_amr_add_options (opt);

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
    RHEA_ASSERT (-1.0 - SC_1000_EPS < EX[il] && EX[il] < 1.0 + SC_1000_EPS);
    RHEA_ASSERT (-1.0 - SC_1000_EPS < EY[il] && EY[il] < 1.0 + SC_1000_EPS);
    RHEA_ASSERT (+1.0 - SC_1000_EPS < EZ[il] && EZ[il] < 2.0 + SC_1000_EPS);
  }
#endif

  switch (tag / 4) {
  case RHEA_DISCR_SHELL_X_TAG_TOP:  /* top (+z) */
    for (il = 0; il < np; ++il) {
      /* transform in-place for nicer grading */
      x = tan (EX[il] * M_PI_4);
      y = tan (EY[il] * M_PI_4);

      /* compute transformation components */
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

static double
rhea_discretization_X_fn_box_spherical_remove_mean (double *E_no_mean,
                                                    const double *E,
                                                    const mangll_locidx_t np)
{
  double              mean = 0.0;
  mangll_locidx_t     il;

  /* compute mean */
  for (il = 0; il < np; ++il) {
    mean += E[il];
  }
  mean *= 1.0 / (double) np;

  /* remove mean */
  for (il = 0; il < np; ++il) {
    E_no_mean[il] = E[il] - mean;
  }

  return mean;
}

/**
 * Geometry transformation for subsections of the full shell domain.
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
  rhea_domain_options_t  *domain_options = data;
  const int           distortion_corr =
                        domain_options->box_spherical_distortion_corr;
  const mangll_tag_t  n_tags_per_part = 4;
  mangll_tag_t        shell_tag;
  double             *E_corr, *E_shift;
  mangll_locidx_t     il;

  /* set tag for shell X-function */
  shell_tag = rhea_discretization_X_fn_box_spherical_get_tag (EX[0], EY[0]);

  /* create work variables */
  if (distortion_corr) {
    E_corr = RHEA_ALLOC (double, np);
  }
  else {
    E_corr = NULL;
  }
  if (shell_tag/n_tags_per_part != RHEA_DISCR_SHELL_X_TAG_TOP) {
    E_shift = RHEA_ALLOC (double, np);
  }
  else {
    E_shift = NULL;
  }

  /* call shell X-function */
  switch (shell_tag/n_tags_per_part) {
  case RHEA_DISCR_SHELL_X_TAG_TOP:  /* top */
    if (!distortion_corr) {
      rhea_discretization_X_fn_shell (shell_tag, np, EX, EY, EZ, X, Y, Z, data);
    }
    else {
      rhea_discretization_X_fn_box_spherical_remove_mean (E_corr, EY, np);
      rhea_discretization_X_fn_shell (shell_tag, np, EX, E_corr, EZ, X, Y, Z,
                                      data);
      rhea_discretization_X_fn_shell (shell_tag, np, EX, EY, EZ, X, E_corr, Z,
                                      data);
    }
    break;
  case RHEA_DISCR_SHELL_X_TAG_BACK:  /* back */
    for (il = 0; il < np; ++il) { /* shift EX from [-3,-1] to [-1,+1] */
      E_shift[il] = EX[il] + 2.0;
    }
    if (!distortion_corr) {
      rhea_discretization_X_fn_shell (shell_tag, np, E_shift, EY, EZ, X, Y, Z,
                                      data);
    }
    else {
      rhea_discretization_X_fn_box_spherical_remove_mean (E_corr, EY, np);
      rhea_discretization_X_fn_shell (shell_tag, np, E_shift, E_corr, EZ,
                                      X, Y, Z, data);
      rhea_discretization_X_fn_shell (shell_tag, np, E_shift, EY, EZ,
                                      X, E_corr, Z, data);
    }
    break;
  case RHEA_DISCR_SHELL_X_TAG_FRONT:  /* front */
    for (il = 0; il < np; ++il) { /* shift EX from [+1,+3] to [-1,+1] */
      E_shift[il] = EX[il] - 2.0;
    }
    if (!distortion_corr) {
      rhea_discretization_X_fn_shell (shell_tag, np, E_shift, EY, EZ, X, Y, Z,
                                      data);
    }
    else {
      rhea_discretization_X_fn_box_spherical_remove_mean (E_corr, EY, np);
      rhea_discretization_X_fn_shell (shell_tag, np, E_shift, E_corr, EZ,
                                      X, Y, Z, data);
      rhea_discretization_X_fn_shell (shell_tag, np, E_shift, EY, EZ,
                                      X, E_corr, Z, data);
    }
    break;
  case RHEA_DISCR_SHELL_X_TAG_LEFT:  /* left */
    for (il = 0; il < np; ++il) { /* shift EY from [-3,-1] to [-1,+1] */
      E_shift[il] = EY[il] + 2.0;
    }
    if (!distortion_corr) {
      rhea_discretization_X_fn_shell (shell_tag, np, EX, E_shift, EZ, X, Y, Z,
                                      data);
    }
    else {
      rhea_discretization_X_fn_box_spherical_remove_mean (E_corr, E_shift, np);
      rhea_discretization_X_fn_shell (shell_tag, np, EX, E_corr, EZ,
                                      X, Y, Z, data);
      rhea_discretization_X_fn_shell (shell_tag, np, EX, E_shift, EZ,
                                      X, E_corr, Z, data);
    }
    break;
  case RHEA_DISCR_SHELL_X_TAG_RIGHT:  /* right */
    for (il = 0; il < np; ++il) { /* shift EY from [+1,+3] to [-1,+1] */
      E_shift[il] = EY[il] - 2.0;
    }
    if (!distortion_corr) {
      rhea_discretization_X_fn_shell (shell_tag, np, EX, E_shift, EZ, X, Y, Z,
                                      data);
    }
    else {
      rhea_discretization_X_fn_box_spherical_remove_mean (E_corr, E_shift, np);
      rhea_discretization_X_fn_shell (shell_tag, np, EX, E_corr, EZ,
                                      X, Y, Z, data);
      rhea_discretization_X_fn_shell (shell_tag, np, EX, E_shift, EZ,
                                      X, E_corr, Z, data);
    }
    break;
  default: /* unknown part of shell domain */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* destroy */
  if (distortion_corr) {
    RHEA_FREE (E_corr);
  }
  if (shell_tag/n_tags_per_part != RHEA_DISCR_SHELL_X_TAG_TOP) {
    RHEA_FREE (E_shift);
  }
}

/**
 * Geometry transformation functions with added topography.
 */
static void
rhea_discretization_X_fn_identity_topo (mangll_tag_t tag, mangll_locidx_t np,
                                        const double *_sc_restrict EX,
                                        const double *_sc_restrict EY,
                                        const double *_sc_restrict EZ,
                                        double *_sc_restrict X,
                                        double *_sc_restrict Y,
                                        double *_sc_restrict Z, void *data)
{
  rhea_topography_options_t *topo_options = data;
  rhea_domain_options_t     *domain_options = topo_options->domain_options;
  const double        z_min = domain_options->z_min;
  const double        z_max = domain_options->z_max;
  const double        z_diff = z_max - z_min;
  double              displ, scaling;
  mangll_locidx_t     il;

  rhea_discretization_X_fn_identity (tag, np, EX, EY, EZ, X, Y, Z, NULL);

  /* scale each node such that
   *   z <- z + displ * (z - z_min)/(z_max - z_min) */
  for (il = 0; il < np; ++il) {
    if (0.0 < Z[il]) {
      displ = rhea_topography_displacement_node (NULL /* label */, X[il], Y[il],
                                                 Z[il], topo_options);
      scaling = 1.0 + displ * (Z[il] - z_min)/z_diff / Z[il];
      Z[il] *= scaling;
    }
  }
}

static void
rhea_discretization_X_fn_shell_topo (mangll_tag_t tag, mangll_locidx_t np,
                                     const double *_sc_restrict EX,
                                     const double *_sc_restrict EY,
                                     const double *_sc_restrict EZ,
                                     double *_sc_restrict X,
                                     double *_sc_restrict Y,
                                     double *_sc_restrict Z, void *data)
{
  rhea_topography_options_t *topo_options = data;
  rhea_domain_options_t     *domain_options = topo_options->domain_options;
  const double        radius_min = domain_options->radius_min;
  const double        radius_max = domain_options->radius_max;
  const double        radius_diff = radius_max - radius_min;
  double              radius, displ, scaling;
  mangll_locidx_t     il;

  rhea_discretization_X_fn_shell (tag, np, EX, EY, EZ, X, Y, Z,
                                  domain_options);

  /* scale each node such that
   *   radius <- radius + displ * (radius - radius_min)/(radius_max - radius_min) */
  for (il = 0; il < np; ++il) {
    radius = rhea_domain_compute_radius (X[il], Y[il], Z[il], domain_options);
    if (0.0 < radius) {
      displ = rhea_topography_displacement_node (NULL /* label */, X[il], Y[il],
                                                 Z[il], topo_options);
      scaling = 1.0 + displ * (radius - radius_min)/radius_diff / radius;
      X[il] *= scaling;
      Y[il] *= scaling;
      Z[il] *= scaling;
    }
  }
}

static void
rhea_discretization_X_fn_box_spherical_topo (mangll_tag_t tag,
                                             mangll_locidx_t np,
                                             const double *_sc_restrict EX,
                                             const double *_sc_restrict EY,
                                             const double *_sc_restrict EZ,
                                             double *_sc_restrict X,
                                             double *_sc_restrict Y,
                                             double *_sc_restrict Z,
                                             void *data)
{
  rhea_topography_options_t *topo_options = data;
  rhea_domain_options_t     *domain_options = topo_options->domain_options;
  const double        radius_min = domain_options->radius_min;
  const double        radius_max = domain_options->radius_max;
  const double        radius_diff = radius_max - radius_min;
  double              radius, displ, scaling;
  mangll_locidx_t     il;

  rhea_discretization_X_fn_box_spherical (tag, np, EX, EY, EZ, X, Y, Z,
                                          domain_options);

  /* scale each node such that
   *   radius <- radius + displ * (radius - radius_min)/(radius_max - radius_min) */
  for (il = 0; il < np; ++il) {
    radius = rhea_domain_compute_radius (X[il], Y[il], Z[il], domain_options);
    if (0.0 < radius) {
      displ = rhea_topography_displacement_node (NULL /* label */, X[il], Y[il],
                                                 Z[il], topo_options);
      scaling = 1.0 + displ * (radius - radius_min)/radius_diff / radius;
      X[il] *= scaling;
      Y[il] *= scaling;
      Z[il] *= scaling;
    }
  }
}

void
rhea_discretization_process_options (rhea_discretization_options_t *opt,
                                     rhea_domain_options_t *domain_options,
                                     rhea_topography_options_t *topo_options)
{
  int                 topo_exists;
  mangll_X_t          X_fn;
  void               *X_data;

  /* set discretization order and mesh refinement levels */
  opt->order = rhea_discretization_order;
  opt->level_min = rhea_discretization_level_min;
  opt->level_max = rhea_discretization_level_max;

  /* set p4est file path */
  opt->p4est_file_path = rhea_discretization_p4est_file_path;

  /* set undefined boundary information */
  opt->boundary = NULL;

  /* check whether topography exists */
  if (topo_options != NULL) {
    topo_exists = rhea_topography_exists (topo_options);
  }
  else {
    topo_exists = 0;
  }

  /* set reference-to-physical transformation function */
  switch (domain_options->shape) {
  case RHEA_DOMAIN_CUBE:
  case RHEA_DOMAIN_BOX:
    if (!topo_exists) {
      X_fn = rhea_discretization_X_fn_identity;
      X_data = NULL;
    }
    else {
      X_fn = rhea_discretization_X_fn_identity_topo;
      X_data = topo_options;
    }
    break;
  case RHEA_DOMAIN_SHELL:
    if (!topo_exists) {
      X_fn = rhea_discretization_X_fn_shell;
      X_data = domain_options;
    }
    else {
      X_fn = rhea_discretization_X_fn_shell_topo;
      X_data = topo_options;
    }
    break;
  case RHEA_DOMAIN_CUBE_SPHERICAL:
  case RHEA_DOMAIN_BOX_SPHERICAL:
    if (!topo_exists) {
      X_fn = rhea_discretization_X_fn_box_spherical;
      X_data = domain_options;
    }
    else {
      X_fn = rhea_discretization_X_fn_box_spherical_topo;
      X_data = topo_options;
    }
    break;
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }
  rhea_discretization_set_X_fn (opt, X_fn, X_data);
}

void
rhea_discretization_set_X_fn (rhea_discretization_options_t *opt,
                              mangll_X_t X_fn, void *X_data)
{
  opt->X_fn = X_fn;
  opt->X_data = X_data;
}

/******************************************************************************
 * Boundary
 *****************************************************************************/

void
rhea_discretization_boundary_create (rhea_discretization_options_t *opt,
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
rhea_discretization_boundary_clear (rhea_discretization_options_t *opt)
{
  rhea_domain_boundary_destroy (opt->boundary);
  opt->boundary = NULL;
}

/******************************************************************************
 * Constructor/Destructor for p4est
 *****************************************************************************/

void
rhea_p4est_init_fn (p4est_t *p4est, p4est_topidx_t tree,
                    p4est_quadrant_t *quadrant)
{
  rhea_p4est_quadrant_data_t *d = quadrant->p.user_data;

  d->amr_flag = RHEA_P4EST_AMR_FLAG_INIT;
}

p4est_t *
rhea_discretization_p4est_new (sc_MPI_Comm mpicomm,
                               rhea_discretization_options_t *opt,
                               rhea_domain_options_t *domain_options)
{
  /* arguments for creating a new p4est object */
  const int           level_min = opt->level_min;
  const int           level_max = opt->level_max;
  const int           n_quadrants_init = 0;
  const int           fill_uniformly = 1;
  const size_t        data_size = sizeof (rhea_p4est_quadrant_data_t);
  const p4est_init_t  init_fn = rhea_p4est_init_fn;
  void               *user_pointer = NULL;
  int                 subdivision_x, subdivision_y, subdivision_z;
  /* p4est objects */
  p4est_connectivity_t *conn;
  p4est_t            *p4est;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /*
   * Read p4est from file
   */

  if (opt->p4est_file_path != NULL) {
    /* read p4est */
    p4est = p4est_load_ext (opt->p4est_file_path, mpicomm,
                            data_size, 0 /* do not read */,
                            0 /* !partition, i.e., load partition dependent */,
                            1 /* broadcast head */, user_pointer, &conn);

    /* initialize user data */
    p4est_reset_data (p4est, data_size, init_fn, user_pointer);

    /* partition for multigrid coarsening */
    p4est_partition_ext (p4est, 1 /* for coarsening */, NULL);

    RHEA_GLOBAL_VERBOSE_FN_END (__func__);

    /* return p4est */
    return p4est;
  }

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

    subdivision_x = domain_options->box_subdivision_x;
    subdivision_y = domain_options->box_subdivision_y;
    subdivision_z = domain_options->box_subdivision_z;
    conn = p8est_connectivity_new_brick (subdivision_x, subdivision_y,
                                         subdivision_z, 0, 0, 0);
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
      const double        subx = (double) subdivision_x;
      const double        suby = (double) subdivision_y;
      const double        subz = (double) subdivision_z;
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

    subdivision_x = domain_options->box_subdivision_y; /* flip x and y */
    subdivision_y = domain_options->box_subdivision_x; /* flip x and y */
    subdivision_z = domain_options->box_subdivision_z;
    conn = p8est_connectivity_new_brick (subdivision_x, subdivision_y,
                                         subdivision_z, 0, 0, 0);
    /* shift coordinates of vertices, so (modified) shell_X fnc can be used:
     *   x component: [0, subdiv_x] -> subdiv_x/subdiv_z * [-0.5, 0.5]
     *   y component: [0, subdiv_y] -> subdiv_y/subdiv_z * [-0.5, 0.5]
     *   z component: [0, subdiv_z] ->                     [ 1  , 2  ]
     */
    {
      const double        subx = (double) subdivision_x;
      const double        suby = (double) subdivision_y;
      const double        subz = (double) subdivision_z;
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
   * Create p4est
   */

  /* create new p4est */
  p4est = p4est_new_ext (mpicomm, conn, n_quadrants_init, level_min,
                         fill_uniformly, data_size, init_fn, user_pointer);

  /* refine */
  rhea_amr_init_refine (p4est, level_min, level_max, domain_options);

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);

  /* return p4est */
  return p4est;
}

void
rhea_discretization_p4est_destroy (p4est_t *p4est)
{
  p4est_connectivity_t *conn = p4est->connectivity;

  /* destroy p4est */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);
}

/******************************************************************************
 * Constructor/Destructor for mangll
 *****************************************************************************/

void
rhea_discretization_mangll_continuous_new (mangll_t **mangll,
                                           mangll_cnodes_t **cnodes,
                                           p4est_t *p4est,
                                           rhea_discretization_options_t *opt)
{
  sc_MPI_Comm         mpicomm = p4est->mpicomm;
  const int           order = opt->order;
  const mangll_refel_quadrature_type_t  quad_type = MANGLL_REFEL_QUAD_GAUSS;

  p4est_ghost_t      *ghost;
  mangll_mesh_t      *mangll_mesh;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* create p4est ghost */
  if (cnodes != NULL) {
    ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
  }
  else {
    ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  }

  /* create mangll mesh */
  mangll_mesh = mangll_p8est_mesh_new_full (p4est, ghost);

  /* assign mapping to physical space */
  mangll_mesh->X_fn = opt->X_fn;
  mangll_mesh->X_data = opt->X_data;

  RHEA_ASSERT (mangll != NULL);
  if (cnodes != NULL) { /* if mesh with continuous nodes should be created */
    /* create continuous node */
    *cnodes = mangll_p8est_cnodes_new (p4est, ghost, order);

    /* create mangll object with continuous geometry */
    *mangll = mangll_new_ext (mpicomm, order, quad_type, mangll_mesh, *cnodes);
  }
  else {
    /* create mangll object with discontinuous geometry */
    *mangll = mangll_new_ext (mpicomm, order, quad_type, mangll_mesh, NULL);
  }

  /* destroy ghost */
  p4est_ghost_destroy (ghost);

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

void
rhea_discretization_mangll_continuous_destroy (mangll_t *mangll,
                                               mangll_cnodes_t *cnodes)
{
  if (mangll != NULL) {
    mangll_destroy (mangll);
  }
  if (cnodes != NULL) {
    mangll_p8est_cnodes_destroy (cnodes);
  }
}

mangll_t *
rhea_discretization_mangll_discontinuous_new (
                                            p4est_t *p4est,
                                            rhea_discretization_options_t *opt)
{
  mangll_t           *mangll;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  rhea_discretization_mangll_continuous_new (&mangll, NULL, p4est, opt);

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);

  return mangll;
}

void
rhea_discretization_mangll_discontinuous_destroy (mangll_t *mangll)
{
  mangll_destroy (mangll);
}

/******************************************************************************
 * Constructor/Destructor for ymir
 *****************************************************************************/

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

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

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

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
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
  rhea_discretization_mangll_continuous_new (&mangll, &cnodes, p4est, opt);

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
    rhea_discretization_mangll_continuous_destroy (ymir_mesh->ma,
                                                   ymir_mesh->cnodes);
    ymir_mesh_destroy (ymir_mesh);
  }
}

/******************************************************************************
 * Coordinates
 *****************************************************************************/

static void
rhea_discretization_set_cont_coordinates (
                                  ymir_vec_t *coordinates,
                                  rhea_domain_coordinate_type_t coord_type,
                                  rhea_domain_options_t *domain_options)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (coordinates);
  ymir_topidx_t       meshid = coordinates->meshnum;
  const ymir_locidx_t n_nodes = ymir_mesh->fmeshes[meshid].Ncn;
  ymir_locidx_t       nodeid;
  ymir_vec_t         *x_vec, *y_vec, *z_vec;

  /* exit if nothing to do */
  if (n_nodes <= 0) {
    return;
  }

  /* get cartesian coordinates */
  x_vec = ymir_face_cvec_new (ymir_mesh, meshid, 1);
  y_vec = ymir_face_cvec_new (ymir_mesh, meshid, 1);
  z_vec = ymir_face_cvec_new (ymir_mesh, meshid, 1);
  ymir_vec_get_coords (x_vec, y_vec, z_vec);

  /* set coordinates */
  if (YMIR_CVEC_STRIDE == YMIR_STRIDE_NODE) {
    const double       *_sc_restrict x = ymir_cvec_index (x_vec, 0, 0);
    const double       *_sc_restrict y = ymir_cvec_index (y_vec, 0, 0);
    const double       *_sc_restrict z = ymir_cvec_index (z_vec, 0, 0);
    double             *_sc_restrict coord =
                          ymir_cvec_index (coordinates, 0, 0);

    for (nodeid = 0; nodeid < n_nodes; nodeid++) {
      rhea_domain_convert_coordinates (&coord[3*nodeid    ],
                                       &coord[3*nodeid + 1],
                                       &coord[3*nodeid + 2],
                                       x[nodeid], y[nodeid], z[nodeid],
                                       coord_type, domain_options);
    }
  }
  else {
    RHEA_ABORT_NOT_REACHED ();
  }

  /* destroy */
  ymir_vec_destroy (x_vec);
  ymir_vec_destroy (y_vec);
  ymir_vec_destroy (z_vec);
}

void
rhea_discretization_write_cont_coordinates (
                                  const char *file_path_txt,
                                  ymir_mesh_t *ymir_mesh,
                                  ymir_topidx_t meshid,
                                  rhea_domain_coordinate_type_t coord_type,
                                  rhea_domain_options_t *domain_options)
{
  ymir_face_mesh_t   *face_mesh = &(ymir_mesh->fmeshes[meshid]);
  const ymir_locidx_t *n_nodes = face_mesh->Ngo;
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  const int           mpisize = ymir_mesh_get_MPI_Comm_size (ymir_mesh);
  int                *segment_offset;
  int                 r;
  ymir_vec_t         *coordinates;
  double             *coord_data;

  /* create segment offsets */
  segment_offset = RHEA_ALLOC (int, mpisize + 1);
  segment_offset[0] = 0;
  for (r = 0; r < mpisize; r++) {
    RHEA_ASSERT ((segment_offset[r] + 3*n_nodes[r]) <= INT_MAX);
    segment_offset[r+1] = segment_offset[r] + (int) 3*n_nodes[r];
  }

  /* set coordinates */
  coordinates = ymir_face_cvec_new (ymir_mesh, meshid, 3);
  rhea_discretization_set_cont_coordinates (coordinates, coord_type,
                                            domain_options);
  if (0 < face_mesh->Ncn) { /* if nodes exist on this rank */
    coord_data = ymir_cvec_index (coordinates, 0, 0);
  }
  else { /* otherwise this rank is empty */
    coord_data = NULL;
  }

  /* write coordiantes */
  rhea_io_mpi_gather_write_double_to_txt (file_path_txt, coord_data,
                                          segment_offset, 3, mpicomm);

  /* destroy */
  ymir_vec_destroy (coordinates);
  RHEA_FREE (segment_offset);
}

void
rhea_discretization_write_cont_coordinates_volume (
                                  const char *file_path_txt,
                                  ymir_mesh_t *ymir_mesh,
                                  rhea_domain_coordinate_type_t coord_type,
                                  rhea_domain_options_t *domain_options)
{
  rhea_discretization_write_cont_coordinates (file_path_txt, ymir_mesh,
                                              YMIR_VOL_MESH, coord_type,
                                              domain_options);
}

void
rhea_discretization_write_cont_coordinates_surface (
                                  const char *file_path_txt,
                                  ymir_mesh_t *ymir_mesh,
                                  rhea_domain_coordinate_type_t coord_type,
                                  rhea_domain_options_t *domain_options)
{
  rhea_discretization_write_cont_coordinates (file_path_txt, ymir_mesh,
                                              RHEA_DOMAIN_BOUNDARY_FACE_TOP,
                                              coord_type, domain_options);
}
