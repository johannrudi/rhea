/** IO
 *
 * Runs I/O tests.
 *
 ******************************************************************************
 * Author:             Johann Rudi <johann@ices.utexas.edu>
 *****************************************************************************/

#include <rhea.h>
#include <example_share_mesh.h>
#include <example_share_stokes.h>
#include <example_share_vtk.h>
#include <ymir_vtk.h> //###DEV###

/******************************************************************************
 * Main Program
 *****************************************************************************/

/**
 * Runs the program.
 */
int
main (int argc, char **argv)
{
  static const char   func_name[] = "io:main";
  /* parallel environment */
  MPI_Comm            mpicomm = sc_MPI_COMM_WORLD;
  int                 mpisize, mpirank, ompsize;
  /* options */
  ymir_options_t               *opt;
  rhea_domain_options_t         domain_options;
  rhea_temperature_options_t    temp_options;
  rhea_weakzone_options_t       weak_options;
  rhea_viscosity_options_t      visc_options;
  rhea_topography_options_t     topo_options;
  rhea_discretization_options_t discr_options;
  rhea_newton_options_t         newton_options;
  /* options local to this program */
  char               *write_vol_coord_file_path_txt;
  char               *write_surf_coord_file_path_txt;
  char               *vtk_input_path;
  /* mesh */
  p4est_t            *p4est;
  ymir_mesh_t        *ymir_mesh;
  ymir_pressure_elem_t  *press_elem;
  /* Stokes */
  rhea_stokes_problem_t *stokes_problem;

  /*
   * Initialize Program
   */

  /* begin program initialization */
  rhea_init_begin (&mpisize, &mpirank, &ompsize, argc, argv, mpicomm);

  /* create options */
  opt = ymir_options_global_new (argv[0] /* program path */);
  rhea_add_options_base (opt);

  /* add options of this program */
  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  /* coordinates output */
  YMIR_OPTIONS_S, "write-volume-coordinates-file-path-txt", '\0',
    &(write_vol_coord_file_path_txt), NULL,
    "Output path for text file with coordinates of continuous nodes in volume",
  YMIR_OPTIONS_S, "write-surface-coordinates-file-path-txt", '\0',
    &(write_surf_coord_file_path_txt), NULL,
    "Output path for text file with coordinates of continuous nodes at surface",

  /* vtk output */
  YMIR_OPTIONS_S, "vtk-write-input-path", '\0',
    &(vtk_input_path), NULL,
    "VTK file path for the input of the Stokes problem",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add sub-options */
  rhea_add_options_all (opt);
  ymir_options_add_suboptions_solver_stokes (opt);

  /* end program initialization */
  rhea_init_end (opt);

  /*
   * Print Environment and Options
   */

  RHEA_GLOBAL_PRODUCTIONF (
      "Into %s (production %i)\n", func_name, rhea_production_run_get ());
  RHEA_GLOBAL_PRODUCTIONF (
      "Parallel environment: MPI size %i, OpenMP size %i\n", mpisize, ompsize);

  /* print & process options */
  ymir_options_print_summary (SC_LP_INFO, opt);
  rhea_process_options_all (&domain_options, &temp_options, &weak_options,
                            &visc_options, &topo_options, &discr_options,
                            &newton_options);

  /*
   * Setup Mesh
   */

  example_share_mesh_new (&p4est, &ymir_mesh, &press_elem, mpicomm,
                          &domain_options, &topo_options, &discr_options);

  /* write coordinates */
  if (write_vol_coord_file_path_txt != NULL) {
    rhea_discretization_write_cont_coordinates_volume (
        write_vol_coord_file_path_txt, ymir_mesh,
        RHEA_DISCRETIZATION_COORDINATE_SPHERICAL_GEO);
  }
  if (write_surf_coord_file_path_txt != NULL) {
    rhea_discretization_write_cont_coordinates_surface (
        write_surf_coord_file_path_txt, ymir_mesh,
        RHEA_DISCRETIZATION_COORDINATE_SPHERICAL_GEO);
  }

  /*
   * Setup Stokes Problem
   */

  example_share_stokes_new (&stokes_problem, &ymir_mesh, &press_elem,
                            &temp_options, &weak_options, &visc_options,
                            &newton_options, NULL, NULL, NULL);

  /* write vtk of input data */
  example_share_vtk_write_input_data (vtk_input_path, stokes_problem,
                                      &temp_options, &visc_options);

  /* write vtk of weak zone data */
  if (rhea_weakzone_exists (&weak_options)) {
    ymir_vec_t         *distance;
    char                path[BUFSIZ];

    distance = rhea_weakzone_new (ymir_mesh);
    rhea_weakzone_compute_distance (distance, &weak_options);

    snprintf (path, BUFSIZ, "%s_weakzone", vtk_input_path);
    ymir_vtk_write (ymir_mesh, path, distance, "distance", NULL);

    rhea_weakzone_destroy (distance);
  }

  /* write vtk of topography data */
  if (rhea_topography_exists (&topo_options)) {
    ymir_vec_t         *displacement, *displacement_surf;
    char                path[BUFSIZ];

    displacement = ymir_dvec_new (ymir_mesh, 1, YMIR_GLL_NODE);
    displacement_surf = rhea_topography_new (ymir_mesh);
    rhea_topography_displacement_vec (displacement, &topo_options);
    rhea_topography_displacement_vec (displacement_surf, &topo_options);

    snprintf (path, BUFSIZ, "%s_topography", vtk_input_path);
    ymir_vtk_write (ymir_mesh, path,
                    displacement, "displacement",
                    displacement_surf, "displacement_surf", NULL);

    ymir_vec_destroy (displacement);
    rhea_topography_destroy (displacement_surf);
  }

  /*
   * Clear Stokes Problem & Mesh
   */

  /* destroy Stokes problem */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);
  example_share_stokes_destroy (stokes_problem, &temp_options, &weak_options,
                                &visc_options);

  /* destroy mesh */
  example_share_mesh_destroy (ymir_mesh, press_elem, p4est, &topo_options,
                              &discr_options);

  /*
   * Finalize
   */

  /* destroy options */
  ymir_options_global_destroy ();

  /* print that this function is ending */
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", func_name);

  /* finalize rhea */
  rhea_finalize ();

  return 0;
}
