/** YIELDING
 *
 * Runs challenging nonlinear Stokes problems that involve yielding.
 *
 ******************************************************************************
 * Author:             Johann Rudi <johann@ices.utexas.edu>
 *****************************************************************************/

#include <rhea.h>
#include <example_share_mesh.h>
#include <example_share_stokes.h>
#include <example_share_vtk.h>

/******************************************************************************
 * Monitoring
 *****************************************************************************/

/* perfomance monitor tags and names */
typedef enum
{
  RHEA_MAIN_PERFMON_SETUP_MESH,
  RHEA_MAIN_PERFMON_SETUP_STOKES,
  RHEA_MAIN_PERFMON_SETUP_SOLVER,
  RHEA_MAIN_PERFMON_SOLVE,
  RHEA_MAIN_PERFMON_TOTAL,
  RHEA_MAIN_PERFMON_N
}
rhea_main_performance_monitor_idx_t;

static const char  *rhea_main_performance_monitor_name[RHEA_MAIN_PERFMON_N] =
{
  "Setup Mesh",
  "Setup Stokes",
  "Setup Solver",
  "Solve",
  "Total"
};

/******************************************************************************
 * Boundary Conditions
 *****************************************************************************/

/* types of boundary conditions */
typedef enum
{
  YIELDING_NZ_DIR_BC_NONE,
  YIELDING_NZ_DIR_BC_INFLOW_LR_CONST  /* const inflow from left & right side */
}
yielding_nonzero_dir_bc_t;

/* data for boundary conditions */
typedef struct yielding_nonzero_dir_bc_data
{
  yielding_nonzero_dir_bc_t type;
  double              inflow_velocity_magn;
}
yielding_nonzero_dir_bc_data_t;

static void
yielding_set_rhs_vel_nonzero_dir_fn (ymir_vec_t *vel_nonzero_dirichlet,
                                     void *data)
{
  yielding_nonzero_dir_bc_data_t *d = data;

  switch (d->type) {
  case YIELDING_NZ_DIR_BC_NONE:
    ymir_vec_set_zero (vel_nonzero_dirichlet);
    break;

  case YIELDING_NZ_DIR_BC_INFLOW_LR_CONST:
    {
      const double    vel = -d->inflow_velocity_magn;
      double normal_flow_vel[RHEA_DOMAIN_BOUNDARY_FACE_N];

      normal_flow_vel[RHEA_DOMAIN_BOUNDARY_FACE_BASE ] = NAN; /* -z */
      normal_flow_vel[RHEA_DOMAIN_BOUNDARY_FACE_TOP  ] = NAN; /* +z */
      normal_flow_vel[RHEA_DOMAIN_BOUNDARY_FACE_SIDE1] = vel; /* -x */
      normal_flow_vel[RHEA_DOMAIN_BOUNDARY_FACE_SIDE2] = vel; /* +x */
      normal_flow_vel[RHEA_DOMAIN_BOUNDARY_FACE_SIDE3] = NAN; /* -y */
      normal_flow_vel[RHEA_DOMAIN_BOUNDARY_FACE_SIDE4] = NAN; /* +y */
      rhea_velocity_nonzero_boundary_set_face_normals_fn (vel_nonzero_dirichlet,
                                                          normal_flow_vel);
    }
    break;

  default: /* unkown boundary condition */
    RHEA_ABORT_NOT_REACHED ();
  }
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
  static const char   func_name[] = "yielding:main";
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
  rhea_composition_options_t    comp_options;
  rhea_discretization_options_t discr_options;
  /* options local to this program */
  yielding_nonzero_dir_bc_t       nonzero_dir_bc_type;
  yielding_nonzero_dir_bc_data_t  nonzero_dir_bc_data;
  double              nonzero_dir_bc_inflow_velocity_magn;
  int                 solver_iter_max;
  double              solver_rel_tol;
  char               *velocity_file_path_bin;
  char               *pressure_file_path_bin;
  char               *bin_solver_path;
  char               *vtk_input_path;
  char               *vtk_solution_path;
  char               *vtk_solver_path;
#ifdef RHEA_USE_CATALYST
  char               *vis_catalyst_script;
#endif
  /* mesh */
  p4est_t              *p4est;
  ymir_mesh_t          *ymir_mesh;
  ymir_pressure_elem_t *press_elem;
  /* Stokes */
  rhea_stokes_problem_t  *stokes_problem;
  ymir_vec_t         *sol_vel_press;
  int                 nonzero_inital_guess;

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

  /* boundary conditions options */
  YMIR_OPTIONS_I, "nonzero-dirichlet-bc-type", '\0',
    &nonzero_dir_bc_type, YIELDING_NZ_DIR_BC_NONE,
    "Type of nonzero Dirichlet boundary conditions",
  YMIR_OPTIONS_D, "nonzero-dirichlet-bc-inflow-velocity-magn", '\0',
    &nonzero_dir_bc_inflow_velocity_magn, 1.0,
    "Nonzero Dir. BC: Magnitude of inflow velocity",

  /* solver options */
  YMIR_OPTIONS_I, "solver-iter-max", '\0',
    &solver_iter_max, 100,
    "Maximum number of iterations for Stokes solver",
  YMIR_OPTIONS_D, "solver-rel-tol", '\0',
    &solver_rel_tol, 1.0e-6,
    "Relative tolerance for Stokes solver",

  /* velocity & pressure (initial guess) intput */
  YMIR_OPTIONS_S, "velocity-file-path", '\0',
    &(velocity_file_path_bin), NULL,
    "Path to a binary file that contains a velocity field",
  YMIR_OPTIONS_S, "pressure-file-path", '\0',
    &(pressure_file_path_bin), NULL,
    "Path to a binary file that contains a pressure field",

  /* binary file output */
  YMIR_OPTIONS_S, "bin-write-solver-path", '\0',
    &(bin_solver_path), NULL,
    "Bin file path for solver internals (e.g., iterations of Newton's method)",

  /* vtk file output */
  YMIR_OPTIONS_S, "vtk-write-input-path", '\0',
    &(vtk_input_path), NULL,
    "VTK file path for the input of the Stokes problem",
  YMIR_OPTIONS_S, "vtk-write-solution-path", '\0',
    &(vtk_solution_path), NULL,
    "VTK file path for the solution of the Stokes problem",
  YMIR_OPTIONS_S, "vtk-write-solver-path", '\0',
    &(vtk_solver_path), NULL,
    "VTK file path for solver internals (e.g., iterations of Newton's method)",

  /* visualization */
#ifdef RHEA_USE_CATALYST
  YMIR_OPTIONS_S, "vis-catalyst-script", '\0',
    &(vis_catalyst_script), NULL,
    "Python script for ParaView-Catalyst",
#endif

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
  rhea_process_options_all (&domain_options, &temp_options, &comp_options,
                            &plate_options, &weak_options, &topo_options,
                            &visc_options, &discr_options);

  /*
   * Setup Mesh
   */

  example_share_mesh_new (&p4est, &ymir_mesh, &press_elem, mpicomm,
                          &domain_options, &topo_options, &discr_options,
                          RHEA_MAIN_PERFMON_SETUP_MESH);

  /*
   * Setup Stokes Problem
   */

  example_share_stokes_new (&stokes_problem, &ymir_mesh, &press_elem,
                            &temp_options, &comp_options, &plate_options,
                            &weak_options, &visc_options,
                            p4est, &discr_options,
                            RHEA_MAIN_PERFMON_SETUP_MESH,
                            RHEA_MAIN_PERFMON_SETUP_STOKES,
                            bin_solver_path, vtk_solver_path);

  /* set nonzero Dirichlet BC's */
  if (YIELDING_NZ_DIR_BC_NONE != nonzero_dir_bc_type) {
    /* set params from options */
    nonzero_dir_bc_data.type = nonzero_dir_bc_type;
    nonzero_dir_bc_data.inflow_velocity_magn =
      nonzero_dir_bc_inflow_velocity_magn;

    /* set callback function for nonzero velocity */
    rhea_stokes_problem_set_rhs_vel_nonzero_dir_compute_fn (
        stokes_problem, yielding_set_rhs_vel_nonzero_dir_fn,
        &nonzero_dir_bc_data);
  }

  /* write vtk of input data */
  example_share_vtk_write_input_data (vtk_input_path, stokes_problem,
                                      &plate_options);

  /*
   * Solve Stokes Problem
   */

  /* setup solver */
  rhea_performance_monitor_start_barrier (RHEA_MAIN_PERFMON_SETUP_SOLVER);
  rhea_stokes_problem_setup_solver (stokes_problem);
  rhea_performance_monitor_stop_add (RHEA_MAIN_PERFMON_SETUP_SOLVER);

  /* initialize solution vector */
  sol_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);
  if (velocity_file_path_bin != NULL && pressure_file_path_bin != NULL) {
    int                 read_success;

    read_success = rhea_velocity_pressure_read (
        sol_vel_press, velocity_file_path_bin, pressure_file_path_bin,
        press_elem, mpicomm);
    if (!read_success) {
      RHEA_LERRORF (
        "%s: Failed reading velocity & pressure from files: \"%s\", \"%s\"\n",
        func_name, velocity_file_path_bin, pressure_file_path_bin);
    }
    nonzero_inital_guess = (read_success != 0);
  }
  else {
    nonzero_inital_guess = 0;
  }

  /* run solver */
  rhea_performance_monitor_start_barrier (RHEA_MAIN_PERFMON_SOLVE);
  rhea_stokes_problem_solve (&sol_vel_press, nonzero_inital_guess,
                             solver_iter_max, solver_rel_tol, stokes_problem);
  rhea_performance_monitor_stop_add (RHEA_MAIN_PERFMON_SOLVE);

  /* write vtk of solution */
  example_share_vtk_write_solution (vtk_solution_path, sol_vel_press,
                                    stokes_problem);

#ifdef RHEA_USE_CATALYST
  if (vis_catalyst_script != NULL) {
    /* initialize visualization */
    rhea_vis_initialize ((const char **) &vis_catalyst_script, 1);

    /* process visualization */
    {
      ymir_vec_t         *velocity, *pressure;

      velocity = rhea_velocity_new (ymir_mesh);
      pressure = rhea_pressure_new (ymir_mesh, press_elem);
      rhea_velocity_pressure_copy_components (velocity, pressure, sol_vel_press,
                                              press_elem);
      rhea_vis_process_primary (velocity, pressure);
      rhea_velocity_destroy (velocity);
      rhea_pressure_destroy (pressure);
    }

    /* finalize visualization */
    rhea_vis_finalize ();
  }
#endif

  /*
   * Clear Stokes Problem & Mesh
   */

  /* destroy solution */
  rhea_velocity_pressure_destroy (sol_vel_press);

  /* destroy Stokes problem */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);
  example_share_stokes_destroy (stokes_problem, &temp_options, &comp_options,
                                &plate_options, &weak_options, &visc_options);

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
                                  RHEA_PERFMON_PRINT_WTIME_ALL,
                                  RHEA_PERFMON_PRINT_NCALLS_ESSENTIAL,
                                  RHEA_PERFMON_PRINT_FLOPS_NONE,
                                  RHEA_PERFMON_PRINT_YMIR_NONE);
  rhea_performance_monitor_finalize ();

  /* destroy options */
  ymir_options_global_destroy ();

  /* print that this function is ending */
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", func_name);

  /* finalize rhea */
  rhea_finalize ();

  return 0;
}
