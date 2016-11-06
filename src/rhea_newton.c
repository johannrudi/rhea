/* Newton residual */
typedef struct rhea_newton_residual
{
  ymir_vec_t         *residual;

  double              res_norm_rtol;

  double              res_norm;
  double             *res_norm_comp;
  double              res_norm_prev;
  double             *res_norm_prev_comp;
  double              res_norm_init;
}
rhea_newton_residual_t;

/* Newton step */
typedef struct rhea_newton_step
{
  int                 iter;

  ymir_vec_t         *step;
  double              step_length;

  int                 lin_iter_max;
  int                 lin_iter_count;
  double              lin_res_norm_rtol;
  double              lin_res_norm_reduction;
  double              lin_convergence;
}
rhea_newton_step_t;

void
rhea_newton_solve (ymir_vec_t *solution)
{
  const char         *this_fn_name = "rhea_newton_solve";
  const int           print_summary = 1; //TODO

  const int           iter_start = ;
  const int           iter_max = ;
  int                 iter;
  rhea_newton_step_t *step;
  rhea_newton_residual_t *residual;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /*
   * Initialize
   */
  {
    //TODO
  }

  /*
   * Iterations Loop
   */

  for (iter = iter_start; iter <= iter_max; iter++) { /* BEGIN: Newton iter */

    /*
     * Pre-Step Output
     */
    {
      //TODO
    }

    /*
     * Newton Step
     */
    {
      /* check for loop exit conditions */
      if (0.0 < res_norm_rtol && res_norm < (res_norm_rtol * res_norm_init)) {
        RHEA_GLOBAL_PRODUCTIONF (
            "%s: Nonlinear problem converged to rtol (%.3e)\n",
            this_fn_name, res_norm_rtol);
        break;
      }
      if (iter_max <= iter) {
        RHEA_GLOBAL_PRODUCTIONF (
            "%s: Maximum number of nonlinear iterations reached (%d)\n",
            this_fn_name, iter_max);
        break;
      }

      /* calculate the accuracy for the linear solve (or "forcing") */
      rhea_newton_set_step_accuracy (step, residual, opt);

      /* solve for Newton step `step` */
      rhea_newton_solve_step (
          /* out: */ step,
          /* in:  */ residual, nl_problem, opt);

      /* perform line search for `step_length` (updates the Jacobian) */
      rhea_newton_search_step_length (
          /* out: */ solution, step, residual,
          /* in:  */ nl_problem, opt);
    }

    /*
     * Post-Step Output
     */
    {
      double              lin_conv;

      /* calculate convergence of the linear solver `rtol^(1/#iter)` */
      lin_conv = exp ( log (lin_res_norm_reduction) / ((double) lin_n_iter) );

      RHEA_GLOBAL_INFOF (
          "%s: Newton iter %d, linear #iter %d, step length %g\n",
          rhea_newton_solve, iter, lin_n_iter, step_length);

      RHEA_GLOBAL_INFOF (
          "%s:   Linear solver residual reduction %.3e (prescribed rtol %.3e), "
          "convergence (rtol^(1/#iter)) %.3e\n",
          rhea_newton_solve, lin_res_norm_reduction, lin_rtol, lin_conv);

      if (print_summary) {
        //TODO print summary at the end
        RHEA_GLOBAL_INFOF (
            "%s: Newton step summary: [%2d, %4d, %.6f, %.3e, %.3e, %.3e]\n",
            rhea_newton_solve, iter, lin_n_iter, step_length,
            lin_res_norm_reduction, lin_rtol, lin_conv);
      }
    }

    /*
     * Setup Next Step
     */
    {
      //TODO
    }

  } /* END: Newton iter */

  /*
   * Finalize
   */
  {
    //TODO
  }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}
