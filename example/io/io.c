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
#include <ymir_vtk.h>

/******************************************************************************
 * Monitoring
 *****************************************************************************/

/* perfomance monitor tags and names */
typedef enum
{
  RHEA_MAIN_PERFMON_SETUP_MESH,
  RHEA_MAIN_PERFMON_SETUP_STOKES,
  RHEA_MAIN_PERFMON_TOTAL,
  RHEA_MAIN_PERFMON_N
}
rhea_main_performance_monitor_idx_t;

static const char  *rhea_main_performance_monitor_name[RHEA_MAIN_PERFMON_N] =
{
  "Setup Mesh",
  "Setup Stokes",
  "Total"
};

/******************************************************************************
 * Functions
 *****************************************************************************/

/**
 * Filters out a specific region and applies that filter to a vector.
 *
 *  east longitude: 170E .. 220E
 *  latitude:       45N .. 65N (colatitude: 45..25)
 *  depth:          200km .. 20km
 */
static const size_t aleutian_n_vertices = 5;
static const float  aleutian_vertices_x[] = {170.0, 220.0, 220.0, 170.0, 170.0};
static const float  aleutian_vertices_y[] = { 45.0,  45.0,  25.0,  25.0,  45.0};
static const float  aleutian_depth_top = 20.0e3;
static const float  aleutian_depth_bottom = 200.0e3;

/* Filter data */
typedef struct rhea_io_filter_node_fn_data
{
  rhea_domain_column_t   *column;
  rhea_domain_options_t  *domain_options;
  int                     n_fields;
}
rhea_io_filter_node_fn_data_t;

/** Filters at a node. */
static void
rhea_io_apply_filter_aleutian_node_fn (double *v, double x, double y, double z,
                                       ymir_locidx_t nodeid, void *data)
{
  rhea_io_filter_node_fn_data_t  *d = data;

  if (!rhea_domain_point_is_in_column (x, y, z, d->column, d->domain_options)) {
    const int           n_fields = d->n_fields;
    int                 fieldid;

    /* set outside values to zero */
    for (fieldid = 0; fieldid < n_fields; fieldid++) {
      v[fieldid] = 0.0;
    }
  }
}

/** Filters a vector. */
static void
rhea_io_apply_filter_aleutian (ymir_vec_t *vec,
                               rhea_domain_options_t *domain_options)
{
  rhea_domain_column_t  column;
  rhea_io_filter_node_fn_data_t data;

  /* set lateral parameters of column */
  column.polygon_vertices_x = RHEA_ALLOC (float, aleutian_n_vertices);
  column.polygon_vertices_y = RHEA_ALLOC (float, aleutian_n_vertices);
  memcpy (column.polygon_vertices_x, aleutian_vertices_x,
          aleutian_n_vertices * sizeof (float));
  memcpy (column.polygon_vertices_y, aleutian_vertices_y,
          aleutian_n_vertices * sizeof (float));
  column.polygon_n_vertices = aleutian_n_vertices;
  column.polygon_coord_type = RHEA_DOMAIN_COORDINATE_SPHERICAL_GEO_DIM;

  /* set radial parameters of column */
  column.radius_min = domain_options->radius_max -
                      aleutian_depth_bottom / domain_options->radius_max_m;
  column.radius_max = domain_options->radius_max -
                      aleutian_depth_top / domain_options->radius_max_m;

  /* set filter data */
  data.column = &column;
  data.domain_options = domain_options;

  /* apply filter */
  if (ymir_vec_is_cvec (vec)) {
    data.n_fields = vec->ncfields;
    ymir_cvec_set_function (vec, rhea_io_apply_filter_aleutian_node_fn, &data);
  }
  else if (ymir_vec_is_dvec (vec)) {
    data.n_fields = vec->ndfields;
    ymir_dvec_set_function (vec, rhea_io_apply_filter_aleutian_node_fn, &data);
  }
  else {
    RHEA_ABORT_NOT_REACHED ();
  }

  /* destroy */
  RHEA_FREE (column.polygon_vertices_x);
  RHEA_FREE (column.polygon_vertices_y);
}

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
  rhea_plate_options_t          plate_options;
  rhea_weakzone_options_t       weak_options;
  rhea_topography_options_t     topo_options;
  rhea_viscosity_options_t      visc_options;
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

  /* initialize performance monitors */
  rhea_performance_monitor_init (rhea_main_performance_monitor_name,
                                 RHEA_MAIN_PERFMON_N);

  /* start performance monitors */
  rhea_performance_monitor_start_barrier (RHEA_MAIN_PERFMON_TOTAL);

  /*
   * Print Environment and Options
   */

  RHEA_GLOBAL_PRODUCTIONF (
      "Into %s (production %i)\n", func_name, rhea_production_run_get ());
  RHEA_GLOBAL_PRODUCTIONF (
      "Parallel environment: MPI size %i, OpenMP size %i\n", mpisize, ompsize);

  /* print & process options */
  ymir_options_print_summary (SC_LP_INFO, opt);
  rhea_process_options_all (&domain_options, &temp_options, &plate_options,
                            &weak_options, &topo_options, &visc_options,
                            &discr_options, &newton_options);

  /*
   * Setup Mesh
   */

  example_share_mesh_new (&p4est, &ymir_mesh, &press_elem, mpicomm,
                          &domain_options, &topo_options, &discr_options,
                          RHEA_MAIN_PERFMON_SETUP_MESH);

  /* write coordinates */
  if (write_vol_coord_file_path_txt != NULL) {
    rhea_discretization_write_cont_coordinates_volume (
        write_vol_coord_file_path_txt, ymir_mesh,
        RHEA_DOMAIN_COORDINATE_SPHERICAL_GEO, &domain_options);
  }
  if (write_surf_coord_file_path_txt != NULL) {
    rhea_discretization_write_cont_coordinates_surface (
        write_surf_coord_file_path_txt, ymir_mesh,
        RHEA_DOMAIN_COORDINATE_SPHERICAL_GEO, &domain_options);
  }

  /*
   * Setup Stokes Problem
   */

  example_share_stokes_new (&stokes_problem, &ymir_mesh, &press_elem,
                            &temp_options, &plate_options,
                            &weak_options, &visc_options,
                            &newton_options, NULL, NULL,
                            RHEA_MAIN_PERFMON_SETUP_MESH,
                            RHEA_MAIN_PERFMON_SETUP_STOKES, NULL, NULL);

  /* write vtk of input data */
  example_share_vtk_write_input_data (vtk_input_path, stokes_problem,
                                      &plate_options);

  /* write vtk of weak zone data */
  if (vtk_input_path != NULL && rhea_weakzone_exists (&weak_options)) {
    ymir_vec_t         *distance;
    char                path[BUFSIZ];

    distance = rhea_weakzone_new (ymir_mesh);
    rhea_weakzone_compute_distance (distance, &weak_options);

    snprintf (path, BUFSIZ, "%s_weakzone", vtk_input_path);
    ymir_vtk_write (ymir_mesh, path, distance, "distance", NULL);

    rhea_weakzone_destroy (distance);
  }

  /* write vtk of topography data */
  if (vtk_input_path != NULL && rhea_topography_exists (&topo_options)) {
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

  /* write vtk of Aleutian */
  if (vtk_input_path != NULL) {
    ymir_vec_t         *aleutian_filter;
    char                path[BUFSIZ];

    aleutian_filter = ymir_dvec_new (ymir_mesh, 1, YMIR_GLL_NODE);
    ymir_vec_set_value (aleutian_filter, 1.0);
    rhea_io_apply_filter_aleutian (aleutian_filter, &domain_options);

    snprintf (path, BUFSIZ, "%s_aleutian", vtk_input_path);
    ymir_vtk_write (ymir_mesh, path,
                    aleutian_filter, "aleutian_filter", NULL);

    ymir_vec_destroy (aleutian_filter);
  }

  /*
   * Clear Stokes Problem & Mesh
   */

  /* destroy Stokes problem */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);
  example_share_stokes_destroy (stokes_problem, &temp_options, &plate_options,
                                &weak_options, &visc_options);

  /* destroy mesh */
  example_share_mesh_destroy (ymir_mesh, press_elem, p4est, &topo_options,
                              &discr_options);

  /*
   * Finalize
   */

  /* stop performance monitors */
  rhea_performance_monitor_stop_add_barrier (RHEA_MAIN_PERFMON_TOTAL);

  /* print performance statistics */
  rhea_performance_monitor_print (func_name,
                                  1 /* print wtime */,
                                  0 /* print #calls */,
                                  0 /* print flops */,
                                  0 /* print ymir */);
  rhea_performance_monitor_finalize ();

  /* destroy options */
  ymir_options_global_destroy ();

  /* print that this function is ending */
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", func_name);

  /* finalize rhea */
  rhea_finalize ();

  return 0;
}
