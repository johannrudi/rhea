
/**************************************
 * Viscosity Computation
 *************************************/
static void
slabs_viscosity_elem (double *_sc_restrict visc_elem,
                        const double *_sc_restrict x,
                        const double *_sc_restrict y,
                        const double *_sc_restrict z,
                        const double *_sc_restrict weak_elem,
                        const int n_nodes_per_el,
                        slabs_options_t *slabs_options)
{
  int                 nodeid;
  double              z_mid = slabs_options->slabs_visc_options->z_lith;

  /* compute viscosity in this element */
  for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
    if (z[nodeid] >= z_mid)  {
      visc_elem[nodeid] = slabs_options->slabs_visc_options->visc_lith;
    }
    else {
      visc_elem[nodeid] = slabs_options->slabs_visc_options->visc_asthen;
    }
    visc_elem[nodeid] *= weak_elem[nodeid];

    /* check viscosity for `nan`, `inf`, and positivity */
    RHEA_ASSERT (isfinite (visc_elem[nodeid]));
    RHEA_ASSERT (0.0 < visc_elem[nodeid]);
  }
}

static void
slabs_viscosity_compute (ymir_vec_t *viscosity,
                         ymir_vec_t *rank1_tensor_scal,
                         ymir_vec_t *bounds_marker,
                         ymir_vec_t *yielding_marker,
                         ymir_vec_t *temperature,
                         ymir_vec_t *weakzone,
                         ymir_vec_t *velocity,
                         rhea_viscosity_options_t *viscosity_options,
                         void *data)
{
  slabs_options_t  *slabs_options = data;
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (viscosity);
  const ymir_locidx_t  n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);
  mangll_t           *mangll = mesh->ma;
  const int           N = ymir_n (mangll->N);

  sc_dmatrix_t       *visc_el_mat, *weak_el_mat;
  double             *x, *y, *z, *tmp_el,*visc_el_data, *weak_el_data;
  ymir_locidx_t       elid;

  /* create work variables */
  weak_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  visc_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  visc_el_data = visc_el_mat->e[0];
  x = RHEA_ALLOC (double, n_nodes_per_el);
  y = RHEA_ALLOC (double, n_nodes_per_el);
  z = RHEA_ALLOC (double, n_nodes_per_el);
  tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) {
    /* get coordinates of this element at Gauss nodes */
    ymir_mesh_get_elem_coord_gauss (x, y, z, elid, mesh, tmp_el);

    weak_el_data = rhea_viscosity_get_elem_gauss (weak_el_mat, weakzone, elid);
    /* compute user defined weak zone viscosity*/
    slabs_viscosity_elem (visc_el_data, x, y, z, weak_el_data, n_nodes_per_el,
                          slabs_options);
    /* set viscosity of this element */
    rhea_viscosity_set_elem_gauss (viscosity, visc_el_mat, elid);
  }

  /* destroy */
  sc_dmatrix_destroy (weak_el_mat);
  sc_dmatrix_destroy (visc_el_mat);
  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);
}

