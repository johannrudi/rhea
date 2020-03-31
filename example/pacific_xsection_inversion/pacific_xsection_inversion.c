/** PACIFIC_XSECTION_INVERSION
 *
 * Runs an inverse problem solver for a mantle convection model that is
 * supported on a cross section of Earth across the Pacific.
 */

#include <rhea.h>
#include <example_share_mesh.h>
#include <example_share_stokes.h>
#include <example_share_vtk.h>
#include <ymir_interp_vec.h>

/******************************************************************************
 * Monitoring
 *****************************************************************************/

/* perfomance monitor tags and names */
typedef enum
{
  RHEA_MAIN_PERFMON_SETUP_MESH,
  RHEA_MAIN_PERFMON_SETUP_STOKES,
  RHEA_MAIN_PERFMON_SOLVE_INVERSION,
  RHEA_MAIN_PERFMON_TOTAL,
  RHEA_MAIN_PERFMON_N
}
rhea_main_performance_monitor_idx_t;

static const char  *rhea_main_performance_monitor_name[RHEA_MAIN_PERFMON_N] =
{
  "Setup Mesh",
  "Setup Stokes",
  "Solve Inversion",
  "Total"
};

/******************************************************************************
 * Boundary Conditions
 *****************************************************************************/

/* types of boundary conditions */
typedef enum
{
  /* default: set from domain options */
  PAXSEC_VEL_BC_DIR_DEFAULT,

  /* Dirichlet in normal direction & zero velocity of South American plate */
  PAXSEC_VEL_BC_DIR_NORMAL_ZERO_1ST_PLATE,

  /* Dirichlet in normal direction & fixed velocity of South American plate */
  PAXSEC_VEL_BC_DIR_NORMAL_FIX_1ST_PLATE
}
paxsec_vel_bc_t;

/* data for boundary conditions */
typedef struct paxsec_vel_bc_data
{
  paxsec_vel_bc_t     type;
  rhea_plate_options_t *plate_options;
}
paxsec_vel_bc_data_t;

/**
 * Sets type of Dirichlet boundary conditions for velocity.
 */
static ymir_dir_code_t
paxsec_vel_bc_dir_fn (double x, double y, double z,
                      double nx, double ny, double nz,
                      ymir_topidx_t face, ymir_locidx_t node_id,
                      void *data)
{
  paxsec_vel_bc_data_t *d = data;

  switch (d->type) {
  case PAXSEC_VEL_BC_DIR_DEFAULT:
    /* set Dirichlet in normal direction for all points */
    return YMIR_VEL_DIRICHLET_NORM;
  case PAXSEC_VEL_BC_DIR_NORMAL_ZERO_1ST_PLATE:
  case PAXSEC_VEL_BC_DIR_NORMAL_FIX_1ST_PLATE:
    {
      const int     plate_label = 0; /* index of South American plate */

      if ( face == RHEA_DOMAIN_BOUNDARY_FACE_TOP &&
           rhea_plate_is_inside (x, y, z, plate_label, d->plate_options) ) {
        /* set Dirichlet in all directions */
        return YMIR_VEL_DIRICHLET_ALL;
      }
      else {
        /* set Dirichlet in normal direction */
        return YMIR_VEL_DIRICHLET_NORM;
      }
    }
  default: /* unknown boundary condition */
    RHEA_ABORT_NOT_REACHED ();
  }
}

/**
 * Sets values for Dirichlet boundary conditions for velocity.
 */
static void
paxsec_vel_bc_dir_nonzero_fn (ymir_vec_t *vel_nonzero_dirichlet, void *data)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vel_nonzero_dirichlet);
  paxsec_vel_bc_data_t *d = data;

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (vel_nonzero_dirichlet));

  /* initialize */
  ymir_vec_set_zero (vel_nonzero_dirichlet);

  /* set nonzero boundary values */
  switch (d->type) {
  case PAXSEC_VEL_BC_DIR_DEFAULT:
  case PAXSEC_VEL_BC_DIR_NORMAL_ZERO_1ST_PLATE:
    /* nothing to do */
    break;
  case PAXSEC_VEL_BC_DIR_NORMAL_FIX_1ST_PLATE:
    {
      const ymir_topidx_t face = RHEA_DOMAIN_BOUNDARY_FACE_TOP;
      const int           plate_label = 0; /* index of South American plate */
      ymir_vec_t         *vel_vol = rhea_velocity_new (ymir_mesh);
      ymir_vec_t         *vel_face = ymir_face_cvec_new (ymir_mesh, face, 3);

      /* set values of face vector */
      rhea_plate_velocity_generate (vel_face, plate_label, d->plate_options);

      /* interpolate face vector to volume vector */
      ymir_interp_vec (vel_face, vel_vol);
      ymir_vec_add (1.0, vel_vol, vel_nonzero_dirichlet);

      /* destroy */
      ymir_vec_destroy (vel_face);
      rhea_velocity_destroy (vel_vol);
    }
    break;
  default: /* unknown boundary condition */
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
  static const char   func_name[] = "pacific_xsection_inversion:main";
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
  /* options local to this program */
  int                 vel_bc_type;
  paxsec_vel_bc_data_t  vel_bc_data;
  char               *bin_solver_path;
  char               *txt_inv_solver_path;
  char               *vtk_input_path;
  char               *vtk_solution_path;
  char               *vtk_solver_path;
  char               *vtk_inv_solver_vol_path;
  char               *vtk_inv_solver_surf_path;
  char               *vtk_inv_param_derivative_path;
  /* mesh */
  p4est_t                *p4est;
  ymir_mesh_t            *ymir_mesh;
  ymir_pressure_elem_t   *press_elem;
  /* Stokes */
  rhea_stokes_problem_t  *stokes_problem;
  /* Inversion */
  rhea_inversion_problem_t *inv_problem;

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

  /* solver options */
  YMIR_OPTIONS_I, "velocity-boundary-condition", '\0',
    &(vel_bc_type), PAXSEC_VEL_BC_DIR_DEFAULT,
    "Type of velocity boundary condition",

  /* binary file output */
  YMIR_OPTIONS_S, "bin-write-solver-path", '\0',
    &(bin_solver_path), NULL,
    "Bin file path for solver internals (e.g., iterations of Newton's method)",

  /* text file output */
  YMIR_OPTIONS_S, "txt-write-inversion-solver-path", '\0',
    &(txt_inv_solver_path), NULL,
    "Text file path for solver internals of the inversion",

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
  YMIR_OPTIONS_S, "vtk-write-inversion-solver-volume-path", '\0',
    &(vtk_inv_solver_vol_path), NULL,
    "VTK file path for solver internals of the inversion (volume fields).",
  YMIR_OPTIONS_S, "vtk-write-inversion-solver-surface-path", '\0',
    &(vtk_inv_solver_surf_path), NULL,
    "VTK file path for solver internals of the inversion (surface fields).",
  YMIR_OPTIONS_S, "vtk-write-inversion-param-derivative-path", '\0',
    &(vtk_inv_param_derivative_path), NULL,
    "VTK file path for testing derivatives of the viscosity w.r.t. parameters.",

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
                            &discr_options);

  /*
   * Setup Mesh
   */

  example_share_mesh_new (&p4est, &ymir_mesh, &press_elem, mpicomm,
                          &domain_options, &topo_options, &discr_options,
                          RHEA_MAIN_PERFMON_SETUP_MESH);

  /*
   * Setup Stokes Problem
   */

  /* set boundary conditions */
  switch (vel_bc_type) {
  case PAXSEC_VEL_BC_DIR_DEFAULT:
    break;
  case PAXSEC_VEL_BC_DIR_NORMAL_ZERO_1ST_PLATE:
  case PAXSEC_VEL_BC_DIR_NORMAL_FIX_1ST_PLATE:
    domain_options.velocity_bc_type = RHEA_DOMAIN_VELOCITY_BC_USER;
    rhea_plate_data_create (&plate_options, ymir_mesh_get_MPI_Comm (ymir_mesh));
    vel_bc_data.type = vel_bc_type;
    vel_bc_data.plate_options = &plate_options;
    rhea_domain_set_user_velocity_dirichlet_bc (
        paxsec_vel_bc_dir_fn, &vel_bc_data, 0 /* zero Dir (nonzero -> RHS) */);
    break;
  default: /* unknown boundary condition */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* setup stokes problem */
  example_share_stokes_new (&stokes_problem, &ymir_mesh, &press_elem,
                            &temp_options, &plate_options,
                            &weak_options, &visc_options,
                            p4est, &discr_options,
                            RHEA_MAIN_PERFMON_SETUP_MESH,
                            RHEA_MAIN_PERFMON_SETUP_STOKES,
                            bin_solver_path, vtk_solver_path);

  /* set callback function for nonzero Dirichlet boundary conditions */
  switch (vel_bc_type) {
  case PAXSEC_VEL_BC_DIR_DEFAULT:
  case PAXSEC_VEL_BC_DIR_NORMAL_ZERO_1ST_PLATE:
    break;
  case PAXSEC_VEL_BC_DIR_NORMAL_FIX_1ST_PLATE:
    rhea_stokes_problem_set_rhs_vel_nonzero_dir_compute_fn (
        stokes_problem, paxsec_vel_bc_dir_nonzero_fn, &vel_bc_data);
    break;
  default: /* unknown boundary condition */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* write vtk of input data */
  example_share_vtk_write_input_data (vtk_input_path, stokes_problem,
                                      &plate_options);

  /* setup Stokes solver */
  rhea_performance_monitor_start_barrier (RHEA_MAIN_PERFMON_SETUP_STOKES);
  rhea_stokes_problem_setup_solver_ext (stokes_problem,
                                        RHEA_INVERSION_KRYLOV_SOLVER_N);
  rhea_performance_monitor_stop_add (RHEA_MAIN_PERFMON_SETUP_STOKES);

  /*
   * Solve Inverse Problem
   */

  /* setup inversion solver */
  inv_problem = rhea_inversion_new (stokes_problem);
  rhea_inversion_set_txt_output (inv_problem, txt_inv_solver_path);
  rhea_inversion_set_vtk_output (inv_problem, vtk_inv_solver_vol_path,
                                 vtk_inv_solver_surf_path);

  /* run inversion solver */
  rhea_performance_monitor_start_barrier (RHEA_MAIN_PERFMON_SOLVE_INVERSION);
  rhea_inversion_solve (inv_problem,
                        0 /* use model param's as initial guess */,
                        NULL /* no parameter vector */);
  rhea_performance_monitor_stop_add_barrier (RHEA_MAIN_PERFMON_SOLVE_INVERSION);

  /* destroy */
  rhea_inversion_destroy (inv_problem);

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
