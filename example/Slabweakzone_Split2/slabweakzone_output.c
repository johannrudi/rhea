/**************************************
 * Output Functions
 *************************************/

/* Write vtk of input data. */
slabs_write_input (ymir_mesh_t *ymir_mesh,
                   rhea_stokes_problem_t *stokes_problem,
                   rhea_temperature_options_t *temp_options,
                   ymir_vec_t *temperature,
                   ymir_vec_t *weakzone,
                   ymir_vec_t *visc_TI_svisc,
                   ymir_vec_t *visc_TI_rotate,
                   const char *vtk_write_input_path)
{
  const char         *this_fn_name = "slabs_write_input";
  ymir_vec_t         *background_temp = rhea_temperature_new (ymir_mesh);
  ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);
  ymir_vec_t         *rhs_vel, *rhs_vel_nonzero_dirichlet;
  char                path[BUFSIZ];

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);
  rhs_vel = rhea_stokes_problem_get_rhs_vel (stokes_problem);

  rhea_temperature_background_compute (background_temp, temp_options);


  rhea_vtk_write_input_data (vtk_write_input_path, temperature,
                             background_temp, weakzone, viscosity, NULL,
                             rhs_vel);

  rhea_temperature_destroy (background_temp);
  rhea_viscosity_destroy (viscosity);

  /* output transversely isotropy parameters*/
  {
    if (visc_TI_svisc != NULL && visc_TI_rotate != NULL) {
      ymir_vec_t         *shear_visc = ymir_vec_clone (visc_TI_svisc);
      ymir_vec_scale (0.5, shear_visc);
      ymir_vec_t *angle = ymir_vec_clone (visc_TI_rotate);
      ymir_vec_scale (180.0/M_PI, angle);
      snprintf (path, BUFSIZ, "%s_anisotropic_viscosity", vtk_write_input_path);
      ymir_vtk_write (ymir_mesh, path,
                      shear_visc, "shear_viscosity",
                      angle, "rotation",
                      NULL);
      rhea_viscosity_destroy (shear_visc);
      rhea_viscosity_destroy (angle);
    }
  }

 // {
 //   rhs_vel_nonzero_dirichlet = rhea_stokes_problem_get_rhs_vel_nonzero_dirichlet (stokes_problem);
 //   if (rhs_vel_nonzero_dirichlet != NULL)
 //     snprintf (path, BUFSIZ, "%s_nonzero_dirichlet", vtk_write_input_path);
 //     ymir_vtk_write (ymir_mesh, path,
 //                     rhs_vel_nonzero_dirichlet, "velocity B.C.",
 //                     NULL);
 // }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

static int
slabs_output_pressure (const char *filepath, ymir_vec_t *pressure)
{
  const char *this_fn_name = "slabs_output_pressure";

  ymir_mesh_t        *mesh = ymir_vec_get_mesh (pressure);
  mangll_cnodes_t    *cnodes = mesh->cnodes;
  const int           N = cnodes->N;
  int                 Np;
  ymir_locidx_t       n_elements;
  ymir_locidx_t       elid;
  ymir_face_mesh_t   *fmesh;
  ymir_locidx_t       Ntotal;
  int                 i, j, k;
  double             *x, *y, *z, *tmp_el;
  ymir_topidx_t       fm = pressure->meshnum;
  int                 mpirank = mesh->ma->mpirank;
  sc_dmatrix_t       *elem;
  double             *elemd;
  FILE               *prsfile;
  char                prsfilename[BUFSIZ];


  if (fm == YMIR_VOL_MESH) {
    Np = (N + 1) * (N + 1) * (N + 1);
    n_elements = cnodes->K;
  }
  else {
    fmesh = &(mesh->fmeshes[fm]);
    Np = (N + 1) * (N + 1);
    n_elements = fmesh->K;
  }
  Ntotal = Np * n_elements;

  /* Have each proc write to its own file */
  if (fm == YMIR_VOL_MESH) {
    snprintf (prsfilename, BUFSIZ, "%s_pressure_%04d", filepath, mpirank);
  }
  else {
    snprintf (prsfilename, BUFSIZ, "%s_pressure_%04d.face%d", filepath, mpirank,
              (int) fm);
  }
  prsfile = fopen (prsfilename, "w");
  if (prsfile == NULL) {
    YMIR_LERRORF ("Could not open %s for output!\n", prsfilename);
    return -1;
  }
  x = RHEA_ALLOC (double, Np);
  y = RHEA_ALLOC (double, Np);
  z = RHEA_ALLOC (double, Np);
  tmp_el = RHEA_ALLOC (double, Np);

  elem = sc_dmatrix_new (Np, 1);
  for (elid = 0; elid < n_elements; elid++) {
    if (pressure->nefields)
      ymir_vec_get_elem_interp (pressure, elem, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_COPY);
    else if (pressure->ndfields)
      ymir_dvec_get_elem (pressure, elem, YMIR_STRIDE_NODE, elid,
                          YMIR_COPY);

    ymir_mesh_get_elem_coord_gauss (x, y, z, elid, mesh, tmp_el);

    elemd = elem->e[0];
    for (i = 0; i < Np; ++i)  {
      fprintf (prsfile, "%24.16e    %24.16e    %24.16e    %24.16e\n", x[i],
             y[i], z[i], elemd[i]);
    }
  }
  if (fclose (prsfile)) {
    YMIR_LERROR ("ymir_vtk: Error closing footer\n");
    return -1;
  }
  prsfile = NULL;

  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);
  sc_dmatrix_destroy (elem);

  return 0;
}
