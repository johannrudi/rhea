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
#include <rhea_vis.h> //TODO put into rhea.h
//#include <ymir_vtk.h> //###DEV###

/******************************************************************************
 * Main Program
 *****************************************************************************/

/**
 * Runs the program.
 */
int
main (int argc, char **argv)
{
  const char         *this_fn_name = "io:main";
  /* parallel environment */
  MPI_Comm            mpicomm = MPI_COMM_WORLD;
  int                 mpisize, mpirank, ompsize;
  /* options */
  ymir_options_t               *opt;
  rhea_domain_options_t         domain_options;
  rhea_temperature_options_t    temp_options;
  rhea_weakzone_options_t       weak_options;
  rhea_viscosity_options_t      visc_options;
  rhea_discretization_options_t discr_options;
  rhea_newton_options_t         newton_options;
  /* options local to this program */
  char               *vtk_write_input_path;
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

  /* vtk output options */
  YMIR_OPTIONS_S, "vtk-write-input-path", '\0',
    &(vtk_write_input_path), NULL,
    "File path for vtk files for the input of the Stokes problem",

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
      "Into %s (production %i)\n", this_fn_name, rhea_get_production_run ());
  RHEA_GLOBAL_PRODUCTIONF (
      "Parallel environment: MPI size %i, OpenMP size %i\n", mpisize, ompsize);

  /* print & process options */
  ymir_options_print_summary (SC_LP_INFO, opt);
  rhea_process_options_all (&domain_options, &temp_options, &weak_options,
                            &visc_options, &discr_options,
                            &newton_options);

  /*
   * Setup Mesh
   */

  example_share_mesh_new (&p4est, &ymir_mesh, &press_elem, mpicomm,
                          &domain_options, &discr_options);

  /*
   * Setup Stokes Problem
   */

  example_share_stokes_new (&stokes_problem, ymir_mesh, press_elem,
                            &domain_options, &temp_options, &visc_options,
                            &newton_options, NULL);

  /* write vtk of input data */
  example_share_vtk_write_input_data (vtk_write_input_path, stokes_problem,
                                      &temp_options, &visc_options);

  /*
   * Finalize
   */

  /* destroy Stokes problem */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);
  example_share_stokes_destroy (stokes_problem);

  /* destroy mesh */
  example_share_mesh_destroy (ymir_mesh, press_elem, p4est, &discr_options);

  /* destroy options */
  ymir_options_global_destroy ();

  /* print that this function is ending */
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

  /* finalize rhea */
  rhea_finalize ();

  return 0;
}
