/*
 */

#include <rhea_viscosity.h>
#include <rhea_base.h>

/**
 * Calculates viscosity for temperature dependent (only) viscosity law.
 *
 *   visc (T) = exp (E * (0.5 - T))
 */
 double
slabs_visc_temp_fn (const double visc_temp_decay, const double temp)
{
  return exp (visc_temp_decay * (0.5 - temp));
}

/**
 * Computes temperature dependent viscosity at a node.
 */
static inline double
slabs_visc_temp_node (const double temp, const double weak,
                      const double scaling, const double visc_temp_decay,
                      slabs_physics_options_t *physics_options,
                      const int restrict_to_bounds)
{
  const slabs_viscosity_model_t  viscosity_model =
                        physics_options->viscosity_model_type;
  const double        visc_min = physics_options->viscosity_min;
  const double        visc_max = physics_options->viscosity_max;
  const double        visc_temp_max = physics_options->viscosity_temp_max;

  double              visc_node;

  /* compute viscosity */
  visc_node = scaling * slabs_visc_temp_fn (visc_temp_decay, temp);

  /* apply weak zone and restrict to bounds */
  switch (viscosity_model) {
  case SL_VISCOSITY_MODEL_WYUL:
  case SL_VISCOSITY_MODEL_UWYUL:
    /* restrict viscosity to upper bound for temp. dep. viscosity */
    if (   viscosity_model == SL_VISCOSITY_MODEL_UWYUL && restrict_to_bounds
        && 0.0 < visc_temp_max && visc_temp_max < visc_node ) {
      visc_node = visc_temp_max;
    }

    /* multiply by weak zone */
    visc_node *= weak;

    /* restrict viscosity to upper bound */
    if ( restrict_to_bounds && 0.0 < visc_max && visc_max < visc_node ) {
      visc_node = visc_max;
    }

    /* restrict viscosity to lower bound */
    if ( restrict_to_bounds && 0.0 < visc_min && visc_node < visc_min ) {
      visc_node = visc_min;
    }
    break;

  case SL_VISCOSITY_MODEL_UWYL:
  case SL_VISCOSITY_MODEL_UYWL:
  case SL_VISCOSITY_MODEL_UYWL_SHIFT:
  case SL_VISCOSITY_MODEL_UWL_IIE_REG:
    /* restrict viscosity to upper bound */
    if ( restrict_to_bounds && 0.0 < visc_max && visc_max < visc_node ) {
      visc_node = visc_max;
    }

    /* multiply by weak zone */
    visc_node *= weak;

    /* restrict viscosity to lower bound */
    if ( restrict_to_bounds && 0.0 < visc_min && visc_node < visc_min ) {
      visc_node = visc_min;
    }
    break;

  case SL_VISCOSITY_MODEL_UWYL_LREG:
  case SL_VISCOSITY_MODEL_UWYL_SHIFT_LREG:
    /* restrict viscosity to upper bound */
    if ( restrict_to_bounds && 0.0 < visc_max && visc_max < visc_node ) {
      visc_node = visc_max;
    }

    /* multiply by weak zone */
    visc_node *= weak;

    /* restrict viscosity to lower bound */
    if ( restrict_to_bounds && 0.0 < visc_min ) {
      visc_node += visc_min;
    }
    break;

  default: /* unknown viscosity model type */
    YMIR_ABORT_NOT_REACHED ();
  }

  /* return temperature dependent viscosity */
  return visc_node;
}

/**
 * Computes temperature dependent viscosity in an element.
 */
void
slabs_visc_temp_elem (sc_dmatrix_t *visc_el_mat,
                      const double *x, const double *y, const double *z,
                      const int *Vmask,
                      const sc_dmatrix_t *temp_el_mat,
                      const sc_dmatrix_t *weak_el_mat,
                      slabs_physics_options_t *physics_options,
                      const int restrict_to_bounds)
{
  const double        upper_mantle_radius =
                        physics_options->viscosity_upper_mantle_radius;
  const double        transition_zone =
                        physics_options->viscosity_lower_upper_transition_zone;
  const int           n_nodes_per_el = visc_el_mat->m;
  const double       *temp_el_data = temp_el_mat->e[0];
  double             *visc_el_data = visc_el_mat->e[0];
  double              scaling, visc_temp_decay;
  int                 nodeid;

  /* check input parameters */
  YMIR_ASSERT (visc_el_mat->n == 1);
  YMIR_ASSERT (temp_el_mat != NULL);
  YMIR_ASSERT (temp_el_mat->m == n_nodes_per_el);
  YMIR_ASSERT (temp_el_mat->n == visc_el_mat->n);
  YMIR_ASSERT (weak_el_mat == NULL || weak_el_mat->m == n_nodes_per_el);
  YMIR_ASSERT (weak_el_mat == NULL || weak_el_mat->n == visc_el_mat->n);
  YMIR_ASSERT (sc_dmatrix_is_valid (temp_el_mat));

  /* set parameters depending on location in lower or upper mantle */
  if (slabs_physics_elem_in_upper_mantle (x, y, z, Vmask, physics_options)) {
    /* set upper mantle parameters */
    scaling = physics_options->viscosity_scaling;
    visc_temp_decay = physics_options->viscosity_temp_decay;
  }
  else {
    /* set lower mantle parameters */
    scaling = physics_options->viscosity_lower_mantle_scaling;
    visc_temp_decay = physics_options->viscosity_lower_mantle_temp_decay;
  }

  /* compute viscosity in this element */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) { /* loop over all
                                                         * nodes */
    const double        r =
      slabs_compute_radius (x[nodeid], y[nodeid], z[nodeid], physics_options);
    const double        temp = temp_el_data[nodeid];
    double              weak;

    /* check temperature for valid range [0,1] */
    YMIR_ASSERT (0.0 <= temp && temp <= 1.0);

    /* set weak zone */
    if (weak_el_mat != NULL) {
      weak = weak_el_mat->e[nodeid][0];
      YMIR_ASSERT (isfinite (weak));
      YMIR_ASSERT (0.0 < weak && weak <= 1.0);
    }
    else {
      weak = 1.0;
    }

    /* compute viscosity */
    if (upper_mantle_radius <= 0.0 || transition_zone <= 0.0 ||
        transition_zone < fabs (r - upper_mantle_radius)) {
      /* compute viscosity disregarding LM/UM transition zone */
      visc_el_data[nodeid] = slabs_visc_temp_node (temp, weak, scaling,
                                                   visc_temp_decay,
                                                   physics_options,
                                                   restrict_to_bounds);
    }
    else {
      const double        factor =
        (transition_zone + (r - upper_mantle_radius)) / (2.0*transition_zone);
      double              visc_lm, visc_um;

      /* compute viscosity with LM/UM transition zone (convex combination) */
      visc_lm = slabs_visc_temp_node (
          temp, weak, physics_options->viscosity_lower_mantle_scaling,
          physics_options->viscosity_lower_mantle_temp_decay,
          physics_options, restrict_to_bounds);
      visc_um = slabs_visc_temp_node (
          temp, weak, physics_options->viscosity_scaling,
          physics_options->viscosity_temp_decay,
          physics_options, restrict_to_bounds);
      visc_el_data[nodeid] = (1.0 - factor) * visc_lm + factor * visc_um;
    }

    /* check viscosity for `nan`, `inf`, and positivity */
    YMIR_ASSERT (isfinite (visc_el_data[nodeid]));
    YMIR_ASSERT (0.0 < visc_el_data[nodeid]);
  }
}

/**
 * Computes the temperature dependent viscosity vector.
 */
static void
slabs_viscosity_linear (ymir_dvec_t *viscosity,
                        slabs_stokes_state_t *state,
                        slabs_physics_options_t *physics_options)
{
  ymir_cvec_t        *temp_vec = state->temp_vec;
  ymir_dvec_t        *weak_vec = state->weak_vec;
  ymir_mesh_t        *mesh = viscosity->mesh;
  mangll_t           *mangll = mesh->ma;
  const mangll_locidx_t  n_elements = mesh->cnodes->K;
  const int           N = ymir_n (mangll->N);
  const int           n_nodes_per_el = (N + 1) * (N + 1) * (N + 1);

  sc_dmatrix_t       *temp_el_mat;
  sc_dmatrix_t       *weak_el_mat;
  sc_dmatrix_t       *visc_el_mat;
  double             *x, *y, *z, *tmp_el;
  mangll_locidx_t     elid;

  /* check input parameters */
  YMIR_ASSERT (state->temperature != NULL);
  YMIR_ASSERT (state->temp_vec != NULL);
  YMIR_ASSERT (state->weakzone != NULL);
  YMIR_ASSERT (state->weak_vec != NULL);

  /* create work variables */
  temp_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  weak_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  visc_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  x = YMIR_ALLOC (double, n_nodes_per_el);
  y = YMIR_ALLOC (double, n_nodes_per_el);
  z = YMIR_ALLOC (double, n_nodes_per_el);
  tmp_el = YMIR_ALLOC (double, n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    /* get coordinates of this element at Gauss nodes */
    slabs_elem_get_gauss_coordinates (x, y, z, elid, mangll, tmp_el);

    /* get temperature field of this element from state at Gauss nodes */
    ymir_cvec_get_elem_interp (temp_vec, temp_el_mat, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_READ);
    slabs_matrix_bound_values (temp_el_mat, 0.0, 1.0);

    /* get weak zone of this element */
    ymir_dvec_get_elem (weak_vec, weak_el_mat, YMIR_STRIDE_NODE, elid,
                        YMIR_READ);

    /* compute temperature dependent viscosity (restrict to bounds) */
    slabs_visc_temp_elem (visc_el_mat, x, y, z, mangll->refel->Vmask,
                          temp_el_mat, weak_el_mat, physics_options, 1);

    /* set viscosity of this element */
    ymir_dvec_set_elem (viscosity, visc_el_mat, YMIR_STRIDE_NODE, elid,
                        YMIR_SET);
  }

  /* destroy */
  sc_dmatrix_destroy (temp_el_mat);
  sc_dmatrix_destroy (weak_el_mat);
  sc_dmatrix_destroy (visc_el_mat);
  YMIR_FREE (x);
  YMIR_FREE (y);
  YMIR_FREE (z);
  YMIR_FREE (tmp_el);
}
