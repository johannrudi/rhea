/** NEWTON_POLYNOMIAL
 *
 * Runs Newton solver for a problem based on finding the minimial distance to a
 * quadratic polynomial.
 *
 ******************************************************************************
 *
 * GIVEN: Coefficients c0, c1, c2 forming a quadratic polynomial
 *
 *   f(x) := c0 + c1*x + c2*x^2,  x \in \R
 *
 * GIVEN DATA: 2-dimensional Cartesian coordinates
 *
 *   d := [a , b]^T,  d \in \R^2
 *
 * WANT: Find a minimizer `x` for the objective functional:
 *
 *   J(x) := 1/2 * || F(x) - d ||^2
 *
 * where
 *
 *   F(x) := [x , f(x)]^T,  F(x) \in \R^2
 *
 ******************************************************************************
 *
 * GRADIENT:
 *
 *   g(y) := delta_x [J(x)] (y) := [ d J(x + eps*y) / d eps ]_{eps=0}
 *                               = (F(x) - d) * [1 , f'(x)] * y
 *                               = (x - a)*y + f'(x)*(f(x) - b)*y
 *
 * yields
 *
 *   g = (x - a) + f'(x)*(f(x) - b)
 *
 ******************************************************************************
 *
 * HESSIAN:
 *
 *   H(y) := delta_x [g(x)] (y) := [ d g(x + eps*y) / d eps ]_{eps=0}
 *
 * yields
 *
 *   H = 1 + f''(x)*(f(x) - b) + f'(x)*f'(x)
 *
 ******************************************************************************
 * Author:             Johann Rudi <johann@ices.utexas.edu>
 *****************************************************************************/

#include <newton_polynomial_base.h>
#include <rhea.h>
#include <ymir_vtk.h>

/******************************************************************************
 * Main Program
 *****************************************************************************/

/**
 * Sets up the mesh.
 */
static void
newton_polynomial_setup_mesh (p4est_t **p4est,
                              ymir_mesh_t **ymir_mesh,
                              MPI_Comm mpicomm,
                              rhea_domain_options_t *domain_options,
                              rhea_discretization_options_t *discr_options)
{
  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", __func__);

  /* create p4est */
  *p4est = rhea_discretization_p4est_new (mpicomm, discr_options,
                                          domain_options);

  /* set up boundary, store in `discr_options` */
  rhea_discretization_boundary_create (discr_options, *p4est, domain_options);

  /* create ymir mesh and pressure element */
  rhea_discretization_ymir_mesh_new_from_p4est (ymir_mesh, NULL, *p4est,
                                                discr_options);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", __func__);
}

/**
 * Callback function to set values to coordinates.
 */
static void
newton_polynomial_set_x_coord_fn (double *out, double x, double y, double z,
                                  ymir_locidx_t nodeid, void *data)
{
  *out = x;
}

static void
newton_polynomial_set_y_coord_fn (double *out, double x, double y, double z,
                                  ymir_locidx_t nodeid, void *data)
{
  *out = y;
}

/**
 * Sets up a Newton problem.
 */
static void
newton_polynomial_setup_newton (rhea_newton_problem_t **nl_problem,
                                ymir_mesh_t *ymir_mesh)
{
  ymir_vec_t         *tmp_vec = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
  newton_polynomial_problem_t *poly_problem;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", __func__);

  /* initialize data for Newton problem */
  poly_problem = RHEA_ALLOC (newton_polynomial_problem_t, 1);
  poly_problem->sol_x = ymir_vec_template (tmp_vec);

  /* set coefficient of polynomial function */
  {
    const double        start_node = 0.0;
    const double        start_val = 1.0;
    const double        start_deriv = 0.1;
    const double        end_node = 1.0;
    const double        end_val = 0.0;
    double             *coeff;

    coeff = newton_polynomial_new_p2_coeff_hermite (
        start_node, start_val, start_deriv, end_node, end_val);

    poly_problem->coeff[0] = coeff[0];
    poly_problem->coeff[1] = coeff[1];
    poly_problem->coeff[2] = coeff[2];

    RHEA_GLOBAL_INFOF ("%s: coefficient of polynomial (%g, %g, %g)\n",
                       __func__, coeff[0], coeff[1], coeff[2]);

    RHEA_FREE (coeff);
  }

  /* set coordiantes as data */
  poly_problem->data_a = ymir_vec_template (tmp_vec);
  poly_problem->data_b = ymir_vec_template (tmp_vec);
  ymir_dvec_set_function (
      poly_problem->data_a, newton_polynomial_set_x_coord_fn, NULL);
  ymir_dvec_set_function (
      poly_problem->data_b, newton_polynomial_set_y_coord_fn, NULL);

  /* create Newton problem */
  {
    ymir_vec_t         *step_vec = ymir_vec_template (tmp_vec);
    ymir_vec_t         *neg_gradient_vec = ymir_vec_template (tmp_vec);

    *nl_problem = rhea_newton_problem_new (
        newton_polynomial_compute_negative_gradient,
        newton_polynomial_compute_gradient_norm,
        0 /* no multi-component norms */,
        newton_polynomial_solve_hessian_system);
    rhea_newton_problem_set_vectors (
        neg_gradient_vec, step_vec, *nl_problem);
    rhea_newton_problem_set_data_fn (
        poly_problem, newton_polynomial_data_init, NULL /* no clear fnc. */,
        *nl_problem);
    rhea_newton_problem_set_evaluate_objective_fn (
        newton_polynomial_evaluate_objective, 0 /* no multi-component obj */,
        rhea_error_stats_default_fn, *nl_problem);
    rhea_newton_problem_set_apply_hessian_fn (
        newton_polynomial_apply_hessian, *nl_problem);
    rhea_newton_problem_set_update_fn (
        NULL /* updating linearized operator is not necessary */,
        newton_polynomial_update_hessian,
        NULL /* updating RHS of linearized system is not necessary */,
        *nl_problem);
  }
#ifdef RHEA_ENABLE_DEBUG
  rhea_newton_problem_set_checks (1 /* grad */, 1 /* Hessian */, *nl_problem);
#endif

  /* destroy */
  ymir_vec_destroy (tmp_vec);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", __func__);
}

/**
 * Cleans up Newton problem and mesh.
 */
static void
newton_polynomial_setup_clear_all (rhea_newton_problem_t *nl_problem,
                                   p4est_t *p4est,
                                   ymir_mesh_t *ymir_mesh,
                                   rhea_discretization_options_t *discr_options)
{
  newton_polynomial_problem_t *poly_problem;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", __func__);

  /* destroy data of Newton problem */
  poly_problem =
    (newton_polynomial_problem_t *) rhea_newton_problem_get_data (nl_problem);
  ymir_vec_destroy (poly_problem->sol_x);
  ymir_vec_destroy (poly_problem->data_a);
  ymir_vec_destroy (poly_problem->data_b);
  RHEA_FREE (poly_problem);

  /* destroy Newton problem */
  ymir_vec_destroy (rhea_newton_problem_get_neg_gradient_vec (nl_problem));
  ymir_vec_destroy (rhea_newton_problem_get_step_vec (nl_problem));
  rhea_newton_problem_destroy (nl_problem);

  /* destroy mesh */
  rhea_discretization_ymir_mesh_destroy (ymir_mesh, NULL);
  rhea_discretization_p4est_destroy (p4est);

  /* destroy boundary data */
  rhea_discretization_boundary_clear (discr_options);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", __func__);
}

/**
 * Runs the program.
 */
int
main (int argc, char **argv)
{
  static const char   func_name[] = "newton_polynomial:main";
  /* MPI */
  MPI_Comm            mpicomm = sc_MPI_COMM_WORLD;
  int                 mpisize, mpirank, ompsize;
  /* options */
  ymir_options_t     *opt;
  rhea_domain_options_t         domain_options;
  rhea_discretization_options_t discr_options;
  rhea_newton_options_t         newton_options;
  /* options local to this function */
  char               *vtk_solution_path;
  /* mesh */
  p4est_t            *p4est;
  ymir_mesh_t        *ymir_mesh;
  /* Newton */
  rhea_newton_problem_t *nl_problem;
  ymir_vec_t         *solution;

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
  YMIR_OPTIONS_S, "vtk-write-solution-path", '\0',
    &(vtk_solution_path), NULL,
    "VTK file path for the solution of the Newton solver",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add sub-options */
  rhea_add_options_newton (opt);

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
  rhea_process_options_newton (&domain_options, &discr_options,
                               &newton_options);

  /*
   * Setup Mesh
   */

  newton_polynomial_setup_mesh (&p4est, &ymir_mesh, mpicomm, &domain_options,
                                &discr_options);

  /*
   * Setup Newton Problem
   */

  newton_polynomial_setup_newton (&nl_problem, ymir_mesh);

  /*
   * Solve Newton Problem
   */

  /* create solution vector */
  solution = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);

  /* set initial guess */
  ymir_vec_set_value (solution, 0.5);

  /* run Newton solver */
  newton_options.nonzero_initial_guess = 1;
  newton_options.status_verbosity = 2;
  rhea_newton_solve (&solution, nl_problem, &newton_options);

  /*
   * Output of Solution
   */

  if (vtk_solution_path != NULL) {
    newton_polynomial_problem_t *poly_problem =
      (newton_polynomial_problem_t *) rhea_newton_problem_get_data (nl_problem);
    ymir_vec_t         *dist = ymir_vec_template (solution);

    /* compute distance to polynomial */
    newton_polynomial_compute_distance (dist, solution, poly_problem);

    /* write vtk file */
    ymir_vtk_write (ymir_mesh, vtk_solution_path,
                    poly_problem->data_a, "data_a",
                    poly_problem->data_b, "data_b",
                    solution, "solution",
                    dist, "distance", NULL);

    /* destroy */
    ymir_vec_destroy (dist);
  }

  /*
   * Clear Newton Problem & Mesh
   */

  /* destroy */
  ymir_vec_destroy (solution);

  /* destroy Newton problem and mesh */
  newton_polynomial_setup_clear_all (nl_problem, p4est, ymir_mesh,
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
