/** NEWTON_POLYNOMIAL_MESHFREE
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

/******************************************************************************
 * Main Program
 *****************************************************************************/

/**
 * Runs the program.
 */
int
main (int argc, char **argv)
{
  static const char   func_name[] = "newton_polynomial_meshfree:main";
  /* MPI */
  MPI_Comm            mpicomm = sc_MPI_COMM_WORLD;
  int                 mpisize, mpirank, ompsize;
  /* options */
  ymir_options_t     *opt;
  rhea_newton_options_t newton_options;
  /* Newton */
  newton_polynomial_problem_t *poly_problem;
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
  rhea_newton_add_options (opt);

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
  rhea_newton_process_options (&newton_options);

  /*
   * Setup Problem
   */

  /* create solution vector */
  solution = ymir_vec_new_meshfree (1);

  /* initialize data for Newton problem */
  poly_problem = RHEA_ALLOC (newton_polynomial_problem_t, 1);
  poly_problem->sol_x = ymir_vec_template (solution);

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
                       func_name, coeff[0], coeff[1], coeff[2]);

    RHEA_FREE (coeff);
  }

  /* set data of misfit */
  poly_problem->data_a = ymir_vec_template (solution);
  poly_problem->data_b = ymir_vec_template (solution);
  ymir_vec_set_value (poly_problem->data_a, 0.0);
  ymir_vec_set_value (poly_problem->data_b, 0.0);

  /* create Newton problem */
  {
    ymir_vec_t         *step_vec = ymir_vec_template (solution);
    ymir_vec_t         *neg_gradient_vec = ymir_vec_template (solution);

    nl_problem = rhea_newton_problem_new (
        newton_polynomial_compute_negative_gradient,
        newton_polynomial_solve_hessian_system);
    rhea_newton_problem_set_vectors (
        neg_gradient_vec, step_vec, nl_problem);
    rhea_newton_problem_set_data_fn (
        poly_problem, newton_polynomial_data_init, NULL /* no clear fnc. */,
        nl_problem);
    rhea_newton_problem_set_conv_criterion_fn (
        RHEA_NEWTON_CONV_CRITERION_OBJECTIVE,
        newton_polynomial_evaluate_objective,
        newton_polynomial_compute_gradient_norm,
        0 /* no multi-component norms */, nl_problem);
    rhea_newton_problem_set_apply_hessian_fn (
        newton_polynomial_apply_hessian, nl_problem);
    rhea_newton_problem_set_update_fn (
        NULL /* updating linearized operator is not necessary */,
        newton_polynomial_update_hessian,
        NULL /* updating RHS of linearized system is not necessary */,
        nl_problem);
  }
#ifdef RHEA_ENABLE_DEBUG
  rhea_newton_problem_set_checks (1 /* grad */, 1 /* Hessian */, nl_problem);
#endif

  /*
   * Solve
   */

  /* set initial guess */
  ymir_vec_set_value (solution, 0.5);

  /* run Newton solver */
  newton_options.nonzero_initial_guess = 1;
  newton_options.status_verbosity = 2;
  rhea_newton_solve (&solution, nl_problem, &newton_options);

  /*
   * Clear Problem
   */

  /* destroy */
  ymir_vec_destroy (solution);

  /* destroy data of Newton problem */
  ymir_vec_destroy (poly_problem->sol_x);
  ymir_vec_destroy (poly_problem->data_a);
  ymir_vec_destroy (poly_problem->data_b);
  RHEA_FREE (poly_problem);

  /* destroy Newton problem */
  ymir_vec_destroy (rhea_newton_problem_get_neg_gradient_vec (nl_problem));
  ymir_vec_destroy (rhea_newton_problem_get_step_vec (nl_problem));
  rhea_newton_problem_destroy (nl_problem);

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
