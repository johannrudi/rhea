/*
 */

#include <rhea_vtk.h>
#include <rhea_base.h>
#include <rhea_temperature.h>
#include <rhea_weakzone.h>
#include <rhea_viscosity.h>
#include <rhea_velocity.h>
#include <rhea_pressure.h>
#include <ymir_vtk.h>

void
rhea_vtk_write_input_data (const char *filepath,
                           ymir_vec_t *temperature,
                           ymir_vec_t *background_temp,
                           ymir_vec_t *weakzone,
                           ymir_vec_t *viscosity,
                           ymir_vec_t *rhs_vel)
{
  const char         *this_fn_name = "rhea_vtk_write_input_data";
  const int           in_temp = (temperature != NULL ? 1 : 0);
  const int           in_back = (background_temp != NULL ? 1 : 0);
  const int           in_weak = (weakzone != NULL ? 1 : 0);
  const int           in_visc = (viscosity != NULL ? 1 : 0);
  ymir_mesh_t        *mesh;

  RHEA_GLOBAL_INFOF ("Into %s: Write VTK file \"%s\"\n",
                     this_fn_name, filepath);

  /* check input */
  RHEA_ASSERT (in_temp || in_back || in_weak || in_visc || rhs_vel != NULL);

  /* get ymir mesh */
  if (in_temp) {
    mesh = ymir_vec_get_mesh (temperature);
  }
  else if (in_back) {
    mesh = ymir_vec_get_mesh (background_temp);
  }
  else if (in_weak) {
    mesh = ymir_vec_get_mesh (weakzone);
  }
  else if (in_visc) {
    mesh = ymir_vec_get_mesh (viscosity);
  }
  else {
    mesh = ymir_vec_get_mesh (rhs_vel);
  }

  /* create vectors that were not given as input */
  if (!in_temp) {
    temperature = rhea_temperature_new (mesh);
    ymir_vec_set_value (temperature, -1.0);
  }
  if (!in_back) {
    background_temp = rhea_temperature_new (mesh);
    ymir_vec_set_value (background_temp, -1.0);
  }
  if (!in_weak) {
    weakzone = rhea_weakzone_new (mesh);
    ymir_vec_set_value (weakzone, -1.0);
  }
  if (!in_visc) {
    viscosity = rhea_viscosity_new (mesh);
    ymir_vec_set_value (viscosity, -1.0);
  }

  /* write vtk file */
  ymir_vtk_write (mesh, filepath,
                  temperature, "temperature",
                  background_temp, "background_temp",
                  weakzone, "weakzone",
                  viscosity, "viscosity",
                  rhs_vel, "rhs_vel",
                  NULL);

  /* destroy */
  if (!in_temp) {
    rhea_temperature_destroy (temperature);
  }
  if (!in_back) {
    rhea_temperature_destroy (background_temp);
  }
  if (!in_weak) {
    rhea_weakzone_destroy (weakzone);
  }
  if (!in_visc) {
    rhea_viscosity_destroy (viscosity);
  }

  RHEA_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

void
rhea_vtk_write_solution (const char *filepath,
                         ymir_vec_t *velocity,
                         ymir_vec_t *pressure,
                         ymir_vec_t *viscosity)
{
  const char         *this_fn_name = "rhea_vtk_write_solution";
  const int           in_vel = (velocity != NULL ? 1 : 0);
  const int           in_press = (pressure != NULL ? 1 : 0);
  ymir_mesh_t        *mesh;

  RHEA_GLOBAL_INFOF ("Into %s: Write VTK file \"%s\"\n",
                     this_fn_name, filepath);

  /* check input */
  RHEA_ASSERT (in_vel || in_press || viscosity != NULL);

  /* get ymir mesh */
  if (in_vel) {
    mesh = ymir_vec_get_mesh (velocity);
  }
  else if (in_press) {
    mesh = ymir_vec_get_mesh (pressure);
  }
  else {
    mesh = ymir_vec_get_mesh (viscosity);
  }

  /* create vectors that were not given as input */
  if (!in_vel) {
    velocity = rhea_velocity_new (mesh);
    ymir_vec_set_value (velocity, 0.0);
  }
  if (!in_press) {
    pressure = ymir_dvec_new (mesh, 1, YMIR_GLL_NODE); /* since no press_elem */
    ymir_vec_set_value (pressure, 0.0);
  }

  /* write vtk file */
  ymir_vtk_write (mesh, filepath,
                  velocity, "velocity",
                  pressure, "pressure",
                  viscosity, "viscosity",
                  NULL);

  /* destroy */
  if (!in_vel) {
    rhea_velocity_destroy (velocity);
  }
  if (!in_press) {
    ymir_vec_destroy (pressure);
  }

  RHEA_GLOBAL_INFOF ("Done %s\n", this_fn_name);
}

int
rhea_output_pressure(const char *filepath, ymir_vec_t *pressure)
{
  const char *this_fn_name = "rhea_output_pressure";

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

  elem = sc_dmatrix_new (Np, 3);
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

  return 0;
}
