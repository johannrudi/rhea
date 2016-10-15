/*
 */

#include <rhea_viscosity.h>
#include <rhea_base.h>
#include <ymir_vec.h>

/**
 * Calculates linear (i.e., temperature dependent) viscosity term from
 * Arrhenius relationship.
 *
 *   visc (T) = exp (E * (0.5 - T))
 *
 *   E --- activation energy > 0
 *   T --- temperature in [0,1]
 */
static double
rhea_viscosity_linear_arrhenius (const double activation_energy,
                                 const double temp)
{
  return exp (activation_energy * (0.5 - temp));
}

/**
 * Computes linear viscosity at one node.  It incorporates temperature
 * dependence and weak zones.
 */
static double
rhea_viscosity_linear_node (const double temp,
                            const double weak,
                            const double scaling,
                            const double activation_energy,
                            rhea_viscosity_options_t *opt,
                            const int restrict_to_bounds)
{
  const rhea_viscosity_model_t  model = opt->model;
  const double        visc_min = opt->min;
  const double        visc_max = opt->max;

  double              visc_node;

  /* compute viscosity from Arrhenius relationship */
  visc_node = scaling * rhea_viscosity_linear_arrhenius (activation_energy,
                                                         temp);

  /* apply weak zone and restrict to bounds */
  switch (model) {
  case RHEA_VISCOSITY_MODEL_UWYL:
    /* restrict viscosity to upper bound */
    if (restrict_to_bounds && 0.0 < visc_max && visc_max < visc_node) {
      visc_node = visc_max;
    }

    /* multiply by weak zone */
    visc_node *= weak;

    /* restrict viscosity to lower bound */
    if (restrict_to_bounds && 0.0 < visc_min && visc_node < visc_min) {
      visc_node = visc_min;
    }
    break;

  case RHEA_VISCOSITY_MODEL_UWYL_LREG:
  case RHEA_VISCOSITY_MODEL_UWYL_SHIFT_LREG:
    /* restrict viscosity to upper bound */
    if (restrict_to_bounds && 0.0 < visc_max && visc_max < visc_node) {
      visc_node = visc_max;
    }

    /* multiply by weak zone */
    visc_node *= weak;

    /* restrict viscosity to lower bound */
    if (restrict_to_bounds && 0.0 < visc_min) {
      visc_node += visc_min;
    }
    break;

  default: /* unknown viscosity model */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* return viscosity */
  return visc_node;
}

/**
 * Computes linear viscosity in an element.
 */
static void
rhea_viscosity_linear_elem (sc_dmatrix_t *visc_el_mat,
                            const double *_sc_restrict x,
                            const double *_sc_restrict y,
                            const double *_sc_restrict z,
                            const int *Vmask,
                            const sc_dmatrix_t *temp_el_mat,
                            const sc_dmatrix_t *weak_el_mat,
                            rhea_viscosity_options_t *opt,
                            const int restrict_to_bounds)
{
  rhea_domain_options_t  *domain_options = opt->domain_options;
  const double        lm_um_interface_radius =
                        domain_options->lm_um_interface_radius;
  const double        transition_width =
                        domain_options->lm_um_interface_smooth_transition_width;
  const double        um_scaling = opt->upper_mantle_scaling;
  const double        um_activ_energy = opt->upper_mantle_activation_energy;
  const double        lm_scaling = opt->lower_mantle_scaling;
  const double        lm_activ_energy = opt->lower_mantle_activation_energy;
  const int           n_nodes_per_el = visc_el_mat->m;
  const double       *temp_el_data = temp_el_mat->e[0];
  const double       *weak_el_data = ( weak_el_mat != NULL ?
                                       weak_el_mat->e[0] : NULL );
  double             *visc_el_data = visc_el_mat->e[0];
  double              scaling, activ_energy;
  int                 nodeid;

  /* check input */
  RHEA_ASSERT (visc_el_mat->n == 1);
  RHEA_ASSERT (temp_el_mat != NULL);
  RHEA_ASSERT (temp_el_mat->m == visc_el_mat->m);
  RHEA_ASSERT (temp_el_mat->n == visc_el_mat->n);
  RHEA_ASSERT (sc_dmatrix_is_valid (temp_el_mat));
  RHEA_ASSERT (weak_el_mat == NULL || weak_el_mat->m == visc_el_mat->m);
  RHEA_ASSERT (weak_el_mat == NULL || weak_el_mat->n == visc_el_mat->n);
  RHEA_ASSERT (weak_el_mat == NULL || sc_dmatrix_is_valid (weak_el_mat));

  /* set parameters depending on location in lower or upper mantle */
  if (rhea_domain_elem_is_in_upper_mantle (x, y, z, Vmask, domain_options)) {
    scaling = um_scaling;
    activ_energy = um_activ_energy;
  }
  else {
    scaling = lm_scaling;
    activ_energy = lm_activ_energy;
  }

  /* compute viscosity at each node of this element */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    const double        r = rhea_domain_compute_radius (x[nodeid], y[nodeid],
                                                        z[nodeid],
                                                        domain_options);
    const double        temp = temp_el_data[nodeid];
    const double        weak = ( weak_el_mat != NULL ?
                                 weak_el_data[nodeid] : 1.0 );

    /* check temperature for valid range [0,1] */
    RHEA_ASSERT (0.0 <= temp && temp <= 1.0);
    /* check weak zone for valid range (0,1] */
    RHEA_ASSERT (weak_el_mat == NULL || (0.0 < weak && weak <= 1.0));

    /* compute viscosity */
    if (lm_um_interface_radius <= 0.0 || transition_width <= 0.0 ||
        transition_width < fabs (r - lm_um_interface_radius)) {
      /* compute viscosity "sufficiently far" from LM/UM interface or with a
       * discontinuous LM/UM interface */
      visc_el_data[nodeid] = rhea_viscosity_linear_node (
          temp, weak, scaling, activ_energy, opt, restrict_to_bounds);
    }
    else { /* if close to LM/UM interface and must apply smoothing */
      const double        c =
        (transition_width + (r - lm_um_interface_radius)) /
        (2.0 * transition_width);
      double              visc_lm, visc_um;

      /* compute viscosity with smooth transition at LM/UM interface (via a
       * convex combination) */
      visc_lm = rhea_viscosity_linear_node (
          temp, weak, lm_scaling, lm_activ_energy, opt, restrict_to_bounds);
      visc_um = rhea_viscosity_linear_node (
          temp, weak, um_scaling, um_activ_energy, opt, restrict_to_bounds);
      visc_el_data[nodeid] = (1.0 - c) * visc_lm + c * visc_um;
    }

    /* check viscosity for `nan`, `inf`, and positivity */
    RHEA_ASSERT (isfinite (visc_el_data[nodeid]));
    RHEA_ASSERT (0.0 < visc_el_data[nodeid]);
  }
}

/**
 * Computes linear viscosity.
 */
static void
rhea_viscosity_linear_vec (ymir_vec_t *visc_vec,
                           ymir_vec_t *temp_vec,
                           ymir_vec_t *weak_vec,
                           rhea_viscosity_options_t *opt)
{
  const int           restrict_to_bounds = 1;
//TODO provide functions for these:
  ymir_mesh_t        *mesh = visc_vec->mesh;
  mangll_t           *mangll = mesh->ma;
  const ymir_locidx_t  n_elements = mesh->cnodes->K;
  const int           N = ymir_n (mangll->N);
  const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);
  const int          *Vmask = mangll->refel->Vmask;

  sc_dmatrix_t       *temp_el_mat;
  sc_dmatrix_t       *weak_el_mat;
  sc_dmatrix_t       *visc_el_mat;
  double             *x, *y, *z, *tmp_el;
  ymir_locidx_t       elid;

  /* check input */
  YMIR_ASSERT_IS_CVEC (temp_vec);
  YMIR_ASSERT_IS_DVEC (weak_vec);
  YMIR_ASSERT_IS_DVEC (visc_vec);

  /* create work variables */
  temp_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  weak_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  visc_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  x = RHEA_ALLOC (double, n_nodes_per_el);
  y = RHEA_ALLOC (double, n_nodes_per_el);
  z = RHEA_ALLOC (double, n_nodes_per_el);
  tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    /* get coordinates at Gauss nodes */
    //TODO
    //slabs_elem_get_gauss_coordinates (x, y, z, elid, mangll, tmp_el);

    /* get temperature field at Gauss nodes */
    ymir_cvec_get_elem_interp (temp_vec, temp_el_mat, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_READ);
    //TODO
    //slabs_matrix_bound_values (temp_el_mat, 0.0, 1.0);

    /* get weak zone */
    ymir_dvec_get_elem (weak_vec, weak_el_mat, YMIR_STRIDE_NODE, elid,
                        YMIR_READ);

    /* compute linear viscosity */
    rhea_viscosity_linear_elem (visc_el_mat, x, y, z, Vmask, temp_el_mat,
                                weak_el_mat, opt, restrict_to_bounds);

    /* set viscosity */
    ymir_dvec_set_elem (visc_vec, visc_el_mat, YMIR_STRIDE_NODE, elid,
                        YMIR_SET);
  }

  /* destroy */
  sc_dmatrix_destroy (temp_el_mat);
  sc_dmatrix_destroy (weak_el_mat);
  sc_dmatrix_destroy (visc_el_mat);
  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);
}
