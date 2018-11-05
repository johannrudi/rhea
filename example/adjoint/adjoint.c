/* Slabs_Weakzone Example:
 *
 * Cartesian domain.
*/


#include <example_share_mesh.h>
#include <example_share_stokes.h>
#include <example_share_vtk.h>
#include <adjoint_options.h>
#include <adjoint_geometry.h>
#include <adjoint_physics.h>
#include <adjoint_test.h>
#include <adjoint_TI.h>
#include <adjoint_postp.h>
#include <adjoint_vtk.h>
#include <adjoint_io.h>
#include <adjoint_essential.h>
#include <ymir_pressure_vec.h>
#include <ymir_vtk.h>
#include <ymir.h>

/******************************************************************************
 *  Monitoring
 ******************************************************************************/

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
 *  Main Program
 ******************************************************************************/
/**
 * Runs the program.
 */
int
main (int argc, char **argv)
{
  const char         *this_fn_name = "subduction:main";

  /* MPI */
  MPI_Comm            mpicomm = MPI_COMM_WORLD;
  int                 mpisize, mpirank, ompsize;

  /* options */
  ymir_options_t     *opt;
  rhea_domain_options_t         domain_options;
  rhea_temperature_options_t    temp_options;
  rhea_weakzone_options_t       weak_options;
  rhea_viscosity_options_t      visc_options;
  rhea_topography_options_t     topo_options;
  rhea_plate_options_t     plate_options;
  rhea_discretization_options_t discr_options;
  rhea_newton_options_t         newton_options;

  /* subd options */
  subd_para_options_t     subd_para_options;
  subd_temp_options_t     subd_temp_options;
  subd_visc_options_t     subd_visc_options;
  subd_weak_options_t     subd_weak_options;
  subd_surf_options_t     subd_surf_options;
  subd_velbc_options_t    subd_velbc_options;
  subd_test_options_t     subd_test_options;
  subd_adjoint_options_t     subd_adjoint_options;
  subd_options_t          subd_options;

  /* options local to this function */
  int                 nonzero_initial_guess = 0;
  int                 solver_iter_max;
  double              solver_rel_tol;

  char               *vtk_write_input_path;
  /* mesh */
  p4est_t            *p4est;
  ymir_mesh_t        *ymir_mesh;
  ymir_pressure_elem_t  *press_elem;

  /* Stokes */
  rhea_stokes_problem_t *stokes_problem;

  /* Newton & adjoint */
  ymir_vec_t            *solution;
  rhea_newton_problem_t *newton_problem;
  adjoint_problem_t     *adjoint_problem;

  /*
   * Initialize Libraries
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
  YMIR_OPTIONS_I, "solver-iter-max", '\0',
    &solver_iter_max, 100,
    "Maximum number of iterations for Stokes solver",
  YMIR_OPTIONS_D, "solver-rel-tol", '\0',
    &solver_rel_tol, 1.0e-6,
    "Relative tolerance for Stokes solver",

  YMIR_OPTIONS_S, "vtk-write-input-path", '\0',
    &(vtk_write_input_path), NULL,
    "File path for vtk files for the input of the Stokes problem",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add sub-options */
  subduction_add_options (opt);
  subduction_add_vtk_options (opt);
  subduction_add_txt_options (opt);

  rhea_add_options_all (opt);
  ymir_options_add_suboptions_solver_stokes (opt);


  /* end program initialization */
  rhea_init_end (opt);

  /*
   * Print Environment and Options
   */
  RHEA_GLOBAL_PRODUCTIONF (
      "Into %s (production %i)\n", this_fn_name, rhea_production_run_get ());
  RHEA_GLOBAL_PRODUCTIONF (
      "Parallel environment: MPI size %i, OpenMP size %i\n", mpisize, ompsize);

  /* Print & Process Options */
  ymir_options_print_summary (SC_LP_INFO, opt);

  /*process rhea options */
  rhea_process_options_all (&domain_options, &temp_options, &plate_options, &weak_options,
                            &topo_options, &visc_options, &discr_options,
                            &newton_options);

  /*process subduction options */
  subduction_process_options (&subd_options, &domain_options, &subd_para_options,
                              &subd_temp_options, &subd_visc_options,
                              &subd_weak_options, &subd_surf_options,
                              &subd_velbc_options, &subd_test_options,
                              &subd_adjoint_options);


  /*
   * Setup Mesh
   */
//  example_share_mesh_new (&p4est, &ymir_mesh, &press_elem, mpicomm,
//                          &domain_options, &topo_options, &discr_options,
//                          RHEA_MAIN_PERFMON_SETUP_MESH);  //TODO this function causes errors

  subd_setup_mesh (&p4est, &ymir_mesh, &press_elem, mpicomm,
                      &domain_options, &topo_options, &discr_options, &subd_options);

 /*
  * Setup Stokes Problem
  */
  adjoint_stokes_new (&stokes_problem, &ymir_mesh, &press_elem,
                      &domain_options, &temp_options, &weak_options,
                      &visc_options, &subd_options);

  /*
   * Solve Adjonit Problem
   */


  adjoint_problem = adjoint_problem_new (stokes_problem, p4est, ymir_mesh, press_elem,
                          &discr_options, &temp_options, &subd_options, vtk_write_input_path,
                          solver_iter_max, solver_rel_tol);

  adjoint_setup_newton (&newton_problem, adjoint_problem);
//  adjoint_stokes_update (&stokes_problem, p4est, &ymir_mesh, &press_elem,
//                        &discr_options, &temp_options, &subd_options,
//                        vtk_write_input_path);

  newton_options.status_verbosity = 1;
  newton_options.nonzero_initial_guess = 1;
  rhea_newton_problem_set_check_gradient (1, newton_problem);

  /* initialize solution vector */
  solution = ymir_vec_new_meshfree (1);

  rhea_newton_solve (&solution, newton_problem, &newton_options);

//  rhea_stokes_problem_set_velocity_pressure (stokes_problem, adjoint_problem->sol_vel_press);
//  subd_vtk_write (stokes_problem, &subd_options);
//  subd_txt_io (stokes_problem, &subd_options);

  /* destroy */
  ymir_vec_destroy (solution);
   /*
    * Finalize
    */
  /* destroy Stokes problem and mesh */
  adjoint_setup_clear_all (stokes_problem, newton_problem,
                        p4est, ymir_mesh, press_elem,
                        &temp_options, &visc_options, &weak_options,
                        &topo_options, &plate_options, &discr_options);
  /* destroy options */
  ymir_options_global_destroy ();

  /* print that this function is ending */
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

  /* finalize rhea */
  rhea_finalize ();

  return 0;
}
