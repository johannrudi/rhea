#include <subduction_io.h>

char               *txt_write_visual_mesh_path;
char               *txt_write_cart_vol_coord_path;
char               *txt_write_cart_surf_coord_path;

char               *txt_write_temp_path;
char               *txt_write_surfvelo_path;
char               *txt_write_surfveloZ_path;
char               *txt_write_surfvelo2_path;
char               *txt_write_velo_path;
char               *txt_write_topo_path;
char               *txt_write_topodispl_path;
char               *txt_write_edotsym_path;

char               *txt_read_surfvelo_path;
char               *txt_read_edotsym_path;

void subduction_add_txt_options (ymir_options_t * opt)
{
  ymir_options_addv (opt,

  /* vtk output options */
  YMIR_OPTIONS_S, "txt-write-cartesian-volume-coordinates-path", '\0',
    &(txt_write_cart_vol_coord_path), NULL,
    "File path for write temperature txt",

  YMIR_OPTIONS_S, "txt-write-cartesian-surface-coordinates-path", '\0',
    &(txt_write_cart_surf_coord_path), NULL,
    "File path for write temperature txt",

  YMIR_OPTIONS_S, "txt-write-visual-mesh-path", '\0',
    &(txt_write_visual_mesh_path), NULL,
    "File path for write mesh txt for purpose of visualization",

  YMIR_OPTIONS_S, "txt-write-temp-path", '\0',
    &(txt_write_temp_path), NULL,
    "File path for write temperature txt",

  YMIR_OPTIONS_S, "txt-write-surfvelo-path", '\0',
    &(txt_write_surfvelo_path), NULL,
    "File path for write surface velocity txt",

  YMIR_OPTIONS_S, "txt-write-surfveloZ-path", '\0',
    &(txt_write_surfveloZ_path), NULL,
    "File path for write surface velocity txt",

  YMIR_OPTIONS_S, "txt-write-surfvelo2-path", '\0',
    &(txt_write_surfvelo2_path), NULL,
    "File path for write surface velocity txt",

  YMIR_OPTIONS_S, "txt-write-velo-path", '\0',
    &(txt_write_velo_path), NULL,
    "File path for write velocity txt",

  YMIR_OPTIONS_S, "txt-write-topo-path", '\0',
    &(txt_write_topo_path), NULL,
    "File path for topography",

  YMIR_OPTIONS_S, "txt-write-topodispl-path", '\0',
    &(txt_write_topodispl_path), NULL,
    "File path for topography displacement",

  YMIR_OPTIONS_S, "txt-write-edotsym-path", '\0',
    &(txt_write_edotsym_path), NULL,
    "File path for symmetric strain rate, a vector of 6",

  YMIR_OPTIONS_S, "txt-read-surfvelo-path", '\0',
    &(txt_read_surfvelo_path), NULL,
    "File path for read surface velocity txt",

  YMIR_OPTIONS_S, "txt-read-edotsym-path", '\0',
    &(txt_read_edotsym_path), NULL,
    "File path for read surface velocity txt",

   YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */
}

void
subd_txt_io (rhea_stokes_problem_t *stokes_problem, subd_options_t *subd_options)
{
  ymir_mesh_t     *ymir_mesh
                    = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  const char      *this_fn_name = "subd_io_txt";
  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* write Cartesian volume coordinates */
  if (txt_write_cart_vol_coord_path != NULL) {
    rhea_discretization_write_cont_coordinates_volume (
        txt_write_cart_vol_coord_path, ymir_mesh,
        RHEA_DOMAIN_COORDINATE_CARTESIAN,
//        RHEA_DOMAIN_COORDINATE_SPHERICAL_GEO,
        subd_options->domain_options);
  }

  /* write Cartesian volume coordinates */
  if (txt_write_cart_surf_coord_path != NULL) {
    rhea_discretization_write_cont_coordinates_surface (
        txt_write_cart_surf_coord_path, ymir_mesh,
        RHEA_DOMAIN_COORDINATE_CARTESIAN,
        subd_options->domain_options);
  }

  /* write mesh (coordinates and connectivity)
   * for the purpose of visualization*/
  if (txt_write_visual_mesh_path != NULL) {
    subd_visual_mesh_write (ymir_mesh, txt_write_visual_mesh_path);
  }


  /* write temperature */
  if (txt_write_temp_path != NULL) {
    ymir_vec_t *temperature;

    temperature = rhea_stokes_problem_get_temperature (stokes_problem);
    subd_facevec_write (temperature, txt_write_temp_path);
  }

  /* write temperature */
  if (txt_write_velo_path != NULL) {
    ymir_vec_t              *velocity = rhea_velocity_new (ymir_mesh);

    ymir_stokes_vec_get_velocity (
                rhea_stokes_problem_get_velocity_pressure (stokes_problem),
                velocity,
                rhea_stokes_problem_get_press_elem (stokes_problem));

    subd_facevec_write (velocity, txt_write_velo_path);
    rhea_velocity_destroy (velocity);
  }

  /* write surface velocity to be read in as Neumann B.C. in adjoint equations */
  if (txt_write_surfvelo_path != NULL) {
    ymir_vec_t              *velocity = rhea_velocity_new (ymir_mesh);
    ymir_vec_t              *surf_velocity = ymir_face_cvec_new (ymir_mesh,
                                              RHEA_DOMAIN_BOUNDARY_FACE_TOP, 3);
    ymir_stokes_vec_get_velocity (
                rhea_stokes_problem_get_velocity_pressure (stokes_problem),
                velocity,
                rhea_stokes_problem_get_press_elem (stokes_problem));

    ymir_interp_vec (velocity, surf_velocity);
    subd_facevec_write (surf_velocity, txt_write_surfvelo_path);

    rhea_velocity_destroy (velocity);
    ymir_vec_destroy (surf_velocity);
  }

  /* write surface velocity to be read in as Neumann B.C. in adjoint equations */
  if (txt_write_surfveloZ_path != NULL) {
    ymir_vec_t              *velocity = rhea_velocity_new (ymir_mesh);
    ymir_vec_t              *surf_velocity = ymir_face_cvec_new (ymir_mesh,
                                              RHEA_DOMAIN_BOUNDARY_FACE_TOP, 3);
    ymir_vec_t              *surf_velocityZ = ymir_face_cvec_new (ymir_mesh,
                                              RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    ymir_vec_t          *surf_velo2 = ymir_vec_template (surf_velocity);

    ymir_vec_set_zero (surf_velo2);

    sc_dmatrix_t  *surf_velocityZ_mat = surf_velocityZ->cvec;

    ymir_stokes_vec_get_velocity (
                rhea_stokes_problem_get_velocity_pressure (stokes_problem),
                velocity,
                rhea_stokes_problem_get_press_elem (stokes_problem));

    ymir_interp_vec (velocity, surf_velocity);
    ymir_cvec_get_comp (surf_velocity, surf_velocityZ_mat, 2, YMIR_COPY);

    ymir_cvec_set_comp (surf_velo2, surf_velocityZ_mat, 2, YMIR_SET);

    subd_facevec_write (surf_velocityZ, txt_write_surfveloZ_path);
    subd_facevec_write (surf_velo2, txt_write_surfvelo2_path);

    rhea_velocity_destroy (velocity);
    ymir_vec_destroy (surf_velocity);
    ymir_vec_destroy (surf_velocityZ);
    ymir_vec_destroy (surf_velo2);
  }


  if (txt_write_topo_path != NULL) {
    ymir_vec_t        *topography = ymir_face_cvec_new (ymir_mesh,
                            RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    subd_postp_topography (topography, stokes_problem, subd_options);

    subd_facevec_write (topography, txt_write_topo_path);

    ymir_vec_destroy (topography);
  }

  if (txt_write_topodispl_path != NULL) {
    ymir_vec_t        *topo_displ = ymir_face_cvec_new (ymir_mesh,
                            RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    double            r_max = subd_options->domain_options->radius_max;

    subd_postp_topography (topo_displ, stokes_problem, subd_options);
    ymir_vec_shift (-1.0, topo_displ);

    subd_facevec_write (topo_displ, txt_write_topodispl_path);

    ymir_vec_destroy (topo_displ);
  }

  if (txt_write_edotsym_path != NULL)  {
    ymir_vec_t              *velocity = rhea_velocity_new (ymir_mesh);
    ymir_vec_t *edotsym = ymir_dvec_new (ymir_mesh, 6, YMIR_GAUSS_NODE);
    ymir_stokes_vec_get_velocity (
                rhea_stokes_problem_get_velocity_pressure (stokes_problem),
                velocity,
                rhea_stokes_problem_get_press_elem (stokes_problem));
    ymir_velocity_strain_rate (velocity, edotsym, 0);
    subd_vol_dvec_write (edotsym, txt_write_edotsym_path);
    rhea_velocity_destroy (velocity);
    ymir_vec_destroy (edotsym);
  }

  if (txt_read_edotsym_path != NULL) {
    ymir_vec_t              *velocity = rhea_velocity_new (ymir_mesh);
    ymir_vec_t              *viscosity = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t            *visc_stencil = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t              *temperature;
    ymir_vec_t *edotsym0 = ymir_dvec_new (ymir_mesh, 6, YMIR_GAUSS_NODE);
    ymir_vec_t *edotsym1 = ymir_dvec_new (ymir_mesh, 6, YMIR_GAUSS_NODE);
    double  gradient;

    temperature = rhea_stokes_problem_get_temperature (stokes_problem);
    subd_vol_dvec_read (edotsym0, txt_read_edotsym_path);

    ymir_stokes_vec_get_velocity (
                rhea_stokes_problem_get_velocity_pressure (stokes_problem),
                velocity,
                rhea_stokes_problem_get_press_elem (stokes_problem));
    ymir_velocity_strain_rate (velocity, edotsym1, 0);

    rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

//    gradient = subd_postp_adjoint_gradient (edotsym0, edotsym1, temperature, viscosity);
     subd_adjoint_stencil_visc (visc_stencil, subd_options);
     ymir_vec_multiply_in (visc_stencil,viscosity);
    gradient = subd_adjoint_gradient (edotsym0, edotsym1, temperature, viscosity);

RHEA_GLOBAL_PRODUCTIONF ("the gradient of adjoint is %+.16e\n", gradient);

    rhea_velocity_destroy (velocity);
    rhea_viscosity_destroy (viscosity);
    rhea_viscosity_destroy (visc_stencil);
    ymir_vec_destroy (edotsym0);
    ymir_vec_destroy (edotsym1);
  }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**************************************
 * Write txt file Functions
 *************************************/
void
subd_visual_mesh_write (ymir_mesh_t * mesh, const char *file_path_txt)
{
  mangll_cnodes_t    *cnodes = mesh->cnodes;
  const int           N = cnodes->N;
  const int           Nrp = N + 1;
  int                 N3;
  int                 Np;
  ymir_locidx_t       K;
  ymir_locidx_t       il, ik;
  ymir_locidx_t       Ntotal;
  ymir_locidx_t       Ncells;
  int                 i, j, k;
  double             *Xd, *Yd, *Zd;
  const int           dim = 3;
  int                 Nvertices;
  int                 mpirank = mesh->ma->mpirank;
  int                 numProcs = mesh->ma->mpisize;
  sc_dmatrix_t       *elem;
  double             *elemd;
  int                 p;
  FILE                *nodefile;
  FILE                *connectfile;
  char                filename[BUFSIZ];

  N3 = N * N * N;
  Np = (N + 1) * (N + 1) * (N + 1);
  K = cnodes->K;
  Xd = mesh->ma->X->e[0];
  Yd = mesh->ma->Y->e[0];
  Zd = mesh->ma->Z->e[0];
  Nvertices = 8;
  Ntotal = Np * K;
  Ncells = N3 * K;

  snprintf (filename, BUFSIZ, "%s_nodes", file_path_txt);
  nodefile = fopen (filename, "w");
  for (il = 0; il < Ntotal; ++il) {
    fprintf (nodefile, "%16.8e %16.8e %16.8e\n", Xd[il],
             Yd[il], Zd[il]);
  }
  fclose (nodefile);

  snprintf (filename, BUFSIZ, "%s_connectivity", file_path_txt);
  connectfile = fopen (filename, "w");
  for (ik = 0; ik < K; ++ik) {
    for (k = 0; k < N; ++k) {
      for (j = 0; j < N; ++j) {
        for (i = 0; i < N; ++i) {
          fprintf (connectfile,
                   "%lld %lld %lld %lld %lld %lld %lld %lld\n",
                   (long long) ik * Np + k * Nrp * Nrp + j * Nrp + i,
                   (long long) ik * Np + k * Nrp * Nrp + j * Nrp + (i + 1),
                   (long long) ik * Np + k * Nrp * Nrp + (j + 1) * Nrp +
                   (i + 1),
                   (long long) ik * Np + k * Nrp * Nrp + (j + 1) * Nrp + i,
                   (long long) ik * Np + (k + 1) * Nrp * Nrp + j * Nrp + i,
                   (long long) ik * Np + (k + 1) * Nrp * Nrp + j * Nrp +
                   (i + 1),
                   (long long) ik * Np + (k + 1) * Nrp * Nrp + (j +
                                                                1) * Nrp +
                   (i + 1),
                   (long long) ik * Np + (k + 1) * Nrp * Nrp + (j +
                                                                1) * Nrp +
                   i);
        }
      }
    }
  }
  fclose (connectfile);
}

void
subd_gauss_coord_pressure_write (ymir_vec_t *pressure, const char *file_path_txt)
{
  const char *this_fn_name = "subd_pressure_gauss_write";

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
    snprintf (prsfilename, BUFSIZ, "%s_%04d", file_path_txt, mpirank);
  }
  else {
    snprintf (prsfilename, BUFSIZ, "%s_%04d.face%d", file_path_txt, mpirank,
              (int) fm);
  }
  prsfile = fopen (prsfilename, "w");
  if (prsfile == NULL) {
    YMIR_LERRORF ("Could not open %s for output!\n", prsfilename);
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
  }
  prsfile = NULL;

  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);
  sc_dmatrix_destroy (elem);
}

void
subd_coord_topo_write (ymir_vec_t *topography,
                        const char *file_path_txt)
{
  ymir_mesh_t           *ymir_mesh = ymir_vec_get_mesh (topography);
  mangll_cnodes_t       *cnodes = ymir_mesh->cnodes;
  const int             N = cnodes->N;
  int                   Np = (N + 1) * (N + 1);
  int                   mpirank = ymir_mesh->ma->mpirank;

  ymir_locidx_t         fm = topography->meshnum;
  ymir_face_mesh_t      *fmesh = &(ymir_mesh->fmeshes[fm]);
  ymir_locidx_t         n_elements, elid, Ntotal;
  int                   i, j, k;
  double                *Xd, *Yd, *Zd;
  double                *elemd;
  sc_dmatrix_t          *elem;
  FILE                  *outfile;
  char                  outfilename[BUFSIZ];
  const char            *this_fn_name = "subd_coord_topo_wirte";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  n_elements = fmesh->K;
  Ntotal = Np * n_elements;
  Xd = fmesh->X->e[0];
  Yd = fmesh->Y->e[0];
  Zd = fmesh->Z->e[0];

  snprintf (outfilename, BUFSIZ, "%s%04d.face%d", file_path_txt, mpirank,
            (int) fm);

  outfile = fopen (outfilename, "w");
  if (outfile == NULL) {
    YMIR_LERRORF ("Could not open %s for output!\n", outfilename);
    return;
  }

  elem = sc_dmatrix_new (Np, 1);
  elemd = elem->e[0];
  for (elid = 0; elid < n_elements; elid++) {
    ymir_cvec_get_elem_interp (topography, elem, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_COPY); //TODO: GLL or GAUSS?

    for (i = 0; i < Np; ++i)  {
      j = elid * Np + i;
      fprintf (outfile, "%24.16e    %24.16e    %24.16e    %24.16e\n", Xd[j],
             Yd[j], Zd[j], elemd[i]);
    }
  }
  if (fclose (outfile)) {
    YMIR_LERROR ("main: Error closing footer\n");
    return;
  }

  sc_dmatrix_destroy (elem);
  ymir_vec_destroy (topography);
}

/*write out volume ymir_dvec in txt*/
void
subd_vol_dvec_write (ymir_vec_t *volvec, const char *file_path_txt)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (volvec);
  const ymir_locidx_t K = volvec->K;
  const int           Np = volvec->Np;
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  int                 mpisize, mpiret;
  int                 nd = volvec->ndfields;

  double             *volvec_data = ymir_dvec_index (volvec, 0, 0, 0);
  int                *segment_offset;
  int                 r;
  char                * this_fn_name = "subd_vol_dvec_write";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* get parallel environment */
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize); SC_CHECK_MPI (mpiret);

  /* create segment offsets */
  segment_offset = RHEA_ALLOC (int, mpisize + 1);
  segment_offset[0] = 0;
  for (r = 0; r < mpisize; r++) {
    RHEA_ASSERT ((segment_offset[r] + nd * Np * K) <= INT_MAX);
    segment_offset[r+1] = segment_offset[r] + (int) nd * Np * K;
  }

  /* write volvec */
  rhea_io_mpi_gather_write_double_to_txt (file_path_txt, volvec_data,
                                          segment_offset, nd, mpicomm);

  /* destroy */
  RHEA_FREE (segment_offset);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/*write out face mesh ymir_vec in txt. This is more general, can replace subd_volvec_write*/
void
subd_facevec_write (ymir_vec_t *facevec, const char *file_path_txt)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (facevec);
  ymir_topidx_t       meshid = facevec->meshnum;
  ymir_face_mesh_t   *face_mesh = &(ymir_mesh->fmeshes[meshid]);
  const ymir_locidx_t *n_nodes = face_mesh->Ngo;
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  int                 mpisize = ymir_mesh_get_MPI_Comm_size (ymir_mesh);

  int                 nc = facevec->ncfields;
  double             *facevec_data;
  int                *segment_offset;
  int                 r;
  char                * this_fn_name = "subd_facevec_write";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* create segment offsets */
  segment_offset = RHEA_ALLOC (int, mpisize + 1);
  segment_offset[0] = 0;
  for (r = 0; r < mpisize; r++) {
    RHEA_ASSERT ((segment_offset[r] + (int) nc*n_nodes[r]) <= INT_MAX);
    segment_offset[r+1] = segment_offset[r] + (int) nc*n_nodes[r];
  }

  if (face_mesh->Ncn > 0) { /* if nodes exist on this rank */
    facevec_data = ymir_cvec_index (facevec, 0, 0);
  }
  else { /* otherwise this rank is empty */
    facevec_data = NULL;
  }

  /* write facevec */
  rhea_io_mpi_gather_write_double_to_txt (file_path_txt, facevec_data,
                                          segment_offset, nc, mpicomm);

  /* destroy */
  RHEA_FREE (segment_offset);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**************************************
 * Read txt file Functions
 *************************************/

/*read volume ymir_cvec from txt*/
void
subd_vol_cvec_read (ymir_vec_t *volvec, const char *file_path_txt)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (volvec);
  const mangll_locidx_t *n_nodes = ymir_mesh->cnodes->Ngo;
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  int                 mpisize, mpiret;
  int                 nc = volvec->ncfields;

  double             *volvec_data = ymir_cvec_index (volvec, 0, 0);
  int                *segment_offset;
  int                 r;

  /* get parallel environment */
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize); SC_CHECK_MPI (mpiret);

  /* create segment offsets */
  segment_offset = RHEA_ALLOC (int, mpisize + 1);
  segment_offset[0] = 0;
  for (r = 0; r < mpisize; r++) {
    RHEA_ASSERT ((segment_offset[r] + n_nodes[r]) <= INT_MAX);
    segment_offset[r+1] = segment_offset[r] + (int) nc * n_nodes[r];
  }

  /* read volvec */
  rhea_io_mpi_read_scatter_double (volvec_data, segment_offset,
                                     NULL, file_path_txt, mpicomm);
  /* destroy */
  RHEA_FREE (segment_offset);
}

void
subd_vol_dvec_read (ymir_vec_t *volvec, const char *file_path_txt)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (volvec);
  int                 nd = volvec->ndfields;
  ymir_topidx_t       meshid = volvec->meshnum;
  const ymir_locidx_t K = volvec->K;
  const int           Np = volvec->Np;
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  int                 mpisize, mpiret;

  double             *volvec_data = ymir_dvec_index (volvec, 0, 0, 0);
  int                *segment_offset;
  int                 r;
  char                * this_fn_name = "subd_vol_dvec_read";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* get parallel environment */
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize); SC_CHECK_MPI (mpiret);

  /* create segment offsets */
  segment_offset = RHEA_ALLOC (int, mpisize + 1);
  segment_offset[0] = 0;
  for (r = 0; r < mpisize; r++) {
    RHEA_ASSERT ((segment_offset[r] + (int) nd*Np*K) <= INT_MAX);
    segment_offset[r+1] = segment_offset[r] + (int) nd * Np * K;
  }

  rhea_io_mpi_read_scatter_double (volvec_data, segment_offset,
                                     NULL, file_path_txt, mpicomm);

  RHEA_FREE (segment_offset);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}


void
subd_facevec_read (ymir_vec_t *facevec, const char *file_path_txt)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (facevec);
  int                 nc = facevec->ncfields;
  ymir_topidx_t       meshid = facevec->meshnum;
  ymir_face_mesh_t   *face_mesh = &(ymir_mesh->fmeshes[meshid]);
  const ymir_locidx_t *n_nodes = face_mesh->Ngo;
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  int                 mpisize, mpiret;

  double             *neum_data;
  int                *segment_offset;
  int                 r;
  char                * this_fn_name = "subd_facevec_read";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* get parallel environment */
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize); SC_CHECK_MPI (mpiret);

  /* create segment offsets */
  segment_offset = RHEA_ALLOC (int, mpisize + 1);
  segment_offset[0] = 0;
  for (r = 0; r < mpisize; r++) {
    RHEA_ASSERT ((segment_offset[r] + (int) nc*n_nodes[r]) <= INT_MAX);
    segment_offset[r+1] = segment_offset[r] + (int) nc*n_nodes[r];
  }

  /* read facevec */
  if (face_mesh->Ncn > 0) { /* if nodes exist on this rank */
    neum_data = ymir_cvec_index (facevec, 0, 0);
  } else {
    neum_data = NULL;
  }

  rhea_io_mpi_read_scatter_double (neum_data, segment_offset,
                                     NULL, file_path_txt, mpicomm);

  RHEA_FREE (segment_offset);

  /* communicate shared node values */
  ymir_vec_share_owned (facevec);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}


/*
#if 0
if (vtk_write_io_path != NULL) {
    mangll_cnodes_t       *cnodes = ymir_mesh->cnodes;
    const int             N = cnodes->N;
    int                   Np, mpirank;
    ymir_locidx_t         n_elements, elid, Ntotal, fm;
    int                   i, j, k;
    double                *x, *y, *z, *tmp_el, *elemd;
    sc_dmatrix_t          *elem;
    FILE                  *outfile;
      char                  outfilename[BUFSIZ];

      ymir_vec_t            *viscosity = rhea_viscosity_new (ymir_mesh);

      RHEA_GLOBAL_PRODUCTIONF ("In %s: Start vtk_write_io\n", this_fn_name);

      rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

      mpirank = ymir_mesh->ma->mpirank;

      Np = (N + 1) * (N + 1) * (N + 1);
      n_elements = cnodes->K;
      Ntotal = Np * n_elements;

      snprintf (outfilename, BUFSIZ, "%s_visc_%04d", vtk_write_io_path, mpirank);

      outfile = fopen (outfilename, "w");
      if (outfile == NULL) {
        YMIR_LERRORF ("Could not open %s for output!\n", outfilename);
        return -1;
      }

      x = RHEA_ALLOC (double, Np);
      y = RHEA_ALLOC (double, Np);
      z = RHEA_ALLOC (double, Np);
      tmp_el = RHEA_ALLOC (double, Np);

      elem = sc_dmatrix_new (Np, 1);
      elemd = elem->e[0];
      for (elid = 0; elid < n_elements; elid++) {
        elemd =  rhea_viscosity_get_elem_gauss (elem, viscosity, elid);

        ymir_mesh_get_elem_coord_gauss (x, y, z, elid, ymir_mesh, tmp_el);

        for (i = 0; i < Np; ++i)  {
          fprintf (outfile, "%24.16e    %24.16e    %24.16e    %24.16e\n", x[i],
                 y[i], z[i], elemd[i]);
        }
      }
      if (fclose (outfile)) {
        YMIR_LERROR ("main: Error closing footer\n");
        return -1;
      }

      RHEA_FREE (x);
      RHEA_FREE (y);
      RHEA_FREE (z);
      RHEA_FREE (tmp_el);
      sc_dmatrix_destroy (elem);
      ymir_vec_destroy (viscosity);

    }

    if (vtk_write_mpiio_path != NULL) {
      mangll_cnodes_t       *cnodes = ymir_mesh->cnodes;
      const int             N = cnodes->N;
      int                   Np, mpirank, mpisize;
      ymir_locidx_t         n_elements, elid, Ntotal, fm;
      int                   i, j, k;
      double                *x, *y, *z, *tmp_el, *elemd;
      double                *Ax, *Ay, *Az, *Avisc;
      double                *Px, *Py, *Pz, *Pvisc;
      sc_dmatrix_t          *elem;
      FILE                  *outfile;
      char                  outfilename[BUFSIZ];

      ymir_vec_t            *viscosity = rhea_viscosity_new (ymir_mesh);

      RHEA_GLOBAL_PRODUCTIONF ("In %s: Start vtk_write_mpiio\n", this_fn_name);

      rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

      mpirank = ymir_mesh->ma->mpirank;
      mpisize = ymir_mesh->ma->mpisize;

      Np = (N + 1) * (N + 1) * (N + 1);
      n_elements = cnodes->K;
      Ntotal = Np * n_elements;

      snprintf (outfilename, BUFSIZ, "%s_visc_merged", vtk_write_mpiio_path);
      if (mpirank == 0) {
        outfile = fopen (outfilename, "w");
        if (outfile == NULL) {
          YMIR_LERRORF ("Could not open %s for output!\n", outfilename);
          return -1;
        }
        Ax = RHEA_ALLOC (double, mpisize * n_elements * Np);
        Ay = RHEA_ALLOC (double, mpisize * n_elements * Np);
        Az = RHEA_ALLOC (double, mpisize * n_elements * Np);
        Avisc = RHEA_ALLOC (double, mpisize * n_elements * Np);
      }

      Px = RHEA_ALLOC (double, n_elements * Np);
      Py = RHEA_ALLOC (double, n_elements * Np);
      Pz = RHEA_ALLOC (double, n_elements * Np);
      Pvisc = RHEA_ALLOC (double, n_elements * Np);

      x = RHEA_ALLOC (double, Np);
      y = RHEA_ALLOC (double, Np);
      z = RHEA_ALLOC (double, Np);
      tmp_el = RHEA_ALLOC (double, Np);

      elem = sc_dmatrix_new (Np, 1);
      elemd = elem->e[0];
      for (elid = 0; elid < n_elements; elid++) {
        elemd =  rhea_viscosity_get_elem_gauss (elem, viscosity, elid);

        ymir_mesh_get_elem_coord_gauss (x, y, z, elid, ymir_mesh, tmp_el);

        for (i = 0; i < Np; ++i)  {
          Px[elid * Np + i] = x[i];
          Py[elid * Np + i] = y[i];
          Pz[elid * Np + i] = z[i];
          Pvisc[elid * Np + i] = elemd[i];
        }
      }
      MPI_Gather (Px, n_elements*Np, MPI_DOUBLE,
                  Ax, n_elements*Np, MPI_DOUBLE,
                  0, MPI_COMM_WORLD);
      MPI_Gather (Py, n_elements*Np, MPI_DOUBLE,
                  Ay, n_elements*Np, MPI_DOUBLE,
                  0, MPI_COMM_WORLD);
      MPI_Gather (Pz, n_elements*Np, MPI_DOUBLE,
                  Az, n_elements*Np, MPI_DOUBLE,
                  0, MPI_COMM_WORLD);
      MPI_Gather (Pvisc, n_elements*Np, MPI_DOUBLE,
                  Avisc, n_elements*Np, MPI_DOUBLE,
                  0, MPI_COMM_WORLD);

      if (mpirank == 0) {
        for (i = 0; i < mpisize*n_elements*Np; i++) {
          fprintf (outfile, "%24.16e    %24.16e    %24.16e    %24.16e\n", Ax[i],
                   Ay[i], Az[i], Avisc[i]);
        }
        if (fclose (outfile)) {
          YMIR_LERROR ("main: Error closing footer\n");
          return -1;
        }
        RHEA_FREE (Ax);
        RHEA_FREE (Ay);
        RHEA_FREE (Az);
        RHEA_FREE (Avisc);
      }

      char            path[BUFSIZ];

      snprintf (path, BUFSIZ, "%s_original", vtk_write_mpiio_path);
      ymir_vtk_write (ymir_mesh, path,
                      viscosity, "original_visc",
                      NULL);

      RHEA_FREE (x);
      RHEA_FREE (y);
      RHEA_FREE (z);
      RHEA_FREE (tmp_el);
      RHEA_FREE (Px);
      RHEA_FREE (Py);
      RHEA_FREE (Pz);
      RHEA_FREE (Pvisc);
      sc_dmatrix_destroy (elem);
      ymir_vec_destroy (viscosity);
    }


    if (vtk_read_io_path != NULL)  {
      mangll_cnodes_t       *cnodes = ymir_mesh->cnodes;
      const int             N = cnodes->N;
      int                   Np, mpirank;
      ymir_locidx_t         n_elements, elid, Ntotal, fm;
      int                   i, j, k;
      double                *elemd;
      double                tmp;
      sc_dmatrix_t          *elem;
      FILE                  *infile;
      char                  infilename[BUFSIZ];
      ymir_vec_t            *readin_visc = rhea_viscosity_new (ymir_mesh);

      RHEA_GLOBAL_PRODUCTIONF ("In %s: Start vtk_read_io\n", this_fn_name);

      mpirank = ymir_mesh->ma->mpirank;

      Np = (N + 1) * (N + 1) * (N + 1);
      n_elements = cnodes->K;
      Ntotal = Np * n_elements;

      snprintf (infilename, BUFSIZ, "%s_visc_%04d", vtk_read_io_path, mpirank);

      infile = fopen (infilename, "r");
      if (infile == NULL) {
        YMIR_LERRORF ("Could not open %s for reading!\n", infilename);
        return -1;
      }

      elem = sc_dmatrix_new (Np, 1);
      elemd = elem->e[0];
      for (elid = 0; elid < n_elements; elid++) {
        for (i = 0; i < Np; ++i)  {
          fscanf (infile, "%lf %lf %lf %lf\n", &tmp,
                 &tmp, &tmp, &elemd[i]);
        }
      rhea_viscosity_set_elem_gauss (readin_visc, elem, elid);

    }
    if (fclose (infile)) {
      YMIR_LERROR ("main: Error closing footer\n");
      return -1;
    }


    sc_dmatrix_destroy (elem);
    ymir_vec_destroy (readin_visc);
  }

  if (vtk_read_mpiio_path != NULL)  {
    mangll_cnodes_t       *cnodes = ymir_mesh->cnodes;
    const int             N = cnodes->N;
    int                   Np, mpirank;
    ymir_locidx_t         n_elements, elid, Ntotal, fm;
    int                   i, j, k;
    double                *elemd;
    double                *Avisc, *Pvisc;
    double                tmp, norm_diff;
    sc_dmatrix_t          *elem;
    FILE                  *infile;
    char                  infilename[BUFSIZ];
    ymir_vec_t            *readin_visc = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t            *viscosity = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t            *diff_visc = rhea_viscosity_new (ymir_mesh);

    RHEA_GLOBAL_PRODUCTIONF ("In %s: Start vtk_read_mpiio\n", this_fn_name);

    rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

    mpirank = ymir_mesh->ma->mpirank;

    Np = (N + 1) * (N + 1) * (N + 1);
    n_elements = cnodes->K;
    Ntotal = Np * n_elements;

    snprintf (infilename, BUFSIZ, "%s_visc_merged", vtk_read_mpiio_path);
    if (mpirank == 0) {
      infile = fopen (infilename, "r");
      if (infile == NULL) {
        YMIR_LERRORF ("Could not open %s for reading!\n", infilename);
        return -1;
      }
      Avisc = RHEA_ALLOC (double, mpisize * n_elements * Np);
      for (i = 0; i < mpisize*n_elements*Np; ++i)  {
        fscanf (infile, "%lf %lf %lf %lf\n", &tmp,
               &tmp, &tmp, &Avisc[i]);
      }
    }

    Pvisc = RHEA_ALLOC (double, n_elements * Np);
    MPI_Scatter (Avisc, n_elements*Np, MPI_DOUBLE,
                 Pvisc, n_elements*Np, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    elem = sc_dmatrix_new (Np, 1);
    elemd = elem->e[0];
    for (elid = 0; elid < n_elements; elid++) {
      for (i = 0; i < Np; ++i)  {
        elemd[i] = Pvisc[elid * Np + i];
      }
      rhea_viscosity_set_elem_gauss (readin_visc, elem, elid);
    }

    if (mpirank == 0) {
      if (fclose (infile)) {
        YMIR_LERROR ("main: Error closing footer\n");
        return -1;
      }
      RHEA_FREE (Avisc);
    }

    ymir_vec_copy (readin_visc, diff_visc);
    ymir_vec_add (-1.0, viscosity, diff_visc);
    norm_diff = ymir_vec_norm (diff_visc);
    RHEA_GLOBAL_INFOF ("mpiio_read: norm_diff=%f\n", norm_diff);

    char            path[BUFSIZ];

    snprintf (path, BUFSIZ, "%s_readin", vtk_read_mpiio_path);
    ymir_vtk_write (ymir_mesh, path,
                    readin_visc, "readin_visc",
                    diff_visc, "diff_visc",
                    NULL);


    RHEA_FREE (Pvisc);
    sc_dmatrix_destroy (elem);
    ymir_vec_destroy (readin_visc);
    ymir_vec_destroy (viscosity);
    ymir_vec_destroy (diff_visc);
  }


  if (vtk_read_ioface_path != NULL) {
    ymir_vec_t            *surf_normal_stress = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    mangll_cnodes_t       *cnodes = ymir_mesh->cnodes;
    const int             N = cnodes->N;
    int                   Np;
    int                   mpirank = ymir_mesh->ma->mpirank;
    ymir_locidx_t         n_elements, elid, Ntotal;
    ymir_locidx_t         fm = surf_normal_stress->meshnum;
    ymir_face_mesh_t      *fmesh = &(ymir_mesh->fmeshes[fm]);
    int                   i, j, k;
    double                *elemd;
    double                tmp;
    sc_dmatrix_t          *elem;
    FILE                  *infile;
    char                  infilename[BUFSIZ];


    RHEA_GLOBAL_PRODUCTIONF ("In %s: Start vtk_read_ioface\n", this_fn_name);

    Np = (N + 1) * (N + 1);
    n_elements = fmesh->K;
    Ntotal = Np * n_elements;

    snprintf (infilename, BUFSIZ, "%s_nstress_%04d.face%d", vtk_read_ioface_path, mpirank,
              (int) fm);

    infile = fopen (infilename, "r");
    if (infile == NULL) {
      YMIR_LERRORF ("Could not open %s for output!\n", infilename);
      return -1;
    }

    elem = sc_dmatrix_new (Np, 1);
    elemd = elem->e[0];
    for (elid = 0; elid < n_elements; elid++) {
      for (i = 0; i < Np; ++i) {
        fscanf(infile, "%lf %lf %lf %lf\n", &tmp, &tmp, &tmp, &elemd[i]);
      }
      ymir_cvec_set_elem_interp (surf_normal_stress, elem, YMIR_STRIDE_NODE, elid,
                                 YMIR_GAUSS_NODE, YMIR_SET);
    }
    if (fclose (infile)) {
      YMIR_LERROR ("main: Error closing footer\n");
      return -1;
    }

    char        path[BUFSIZ];
    snprintf (path, BUFSIZ, "%s_readin", vtk_read_ioface_path);
    ymir_vtk_write (ymir_mesh, path,
                    surf_normal_stress, "readin_surf_normal_stress",
                    NULL);

    sc_dmatrix_destroy (elem);
    ymir_vec_destroy (surf_normal_stress);
  }

#endif
*/
