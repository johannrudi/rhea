/*
 */

#include <rhea_newton.h>
#include <rhea_base.h>

/* Newton residual */
typedef struct rhea_newton_residual
{
  ymir_vec_t         *residual_vec;

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

  ymir_vec_t         *step_vec;
  double              step_length;

  int                 lin_iter_max;
  int                 lin_iter_count;
  double              lin_res_norm_rtol;
  double              lin_res_norm_reduction;
  double              lin_convergence;
}
rhea_newton_step_t;

/**
 * Calculates the accuracy for a linear solve within a nonlinear iteration
 * (i.e., the Krylov relative tolerance or "forcing").
 */
static void
rhea_newton_set_step_accuracy (step, residual, opt)
{
  const char         *this_fn_name = "rhea_newton_set_step_accuracy";
  const int           iter = step->iter;
  double              lin_rtol;

  if (iter <= opt->adaptive_lin_rtol_init_n_iter) { /* if set init rtol */
    lin_rtol = opt->adaptive_lin_rtol_init;

    RHEA_GLOBAL_INFOF ("%s: Newton iter %d, use initial linear rtol %.3e\n",
                       this_fn_name, iter, lin_rtol);
  }
  else { /* otherwise set adaptive rtol */
    const double        lin_rtol_max = opt->adaptive_lin_rtol_max;
    const double        exponent = opt->adaptive_lin_rtol_exponent;
    const double        res_norm = residual->res_norm;
    const double        res_norm_prev = residual->res_norm_prev;

    const int           prog_reduction_active =
                          (0 < opt->adaptive_lin_rtol_progressive_n_iter);
    double              prog_reduction = 1.0;

    const int           lin_rtol_min_active =
                          opt->adaptive_lin_rtol_min_active;
    const double        lin_rtol_min_thres =
                          opt->adaptive_lin_rtol_min_threshold;
    double              lin_rtol_min = 0.0;

    /* set the relative tolerance for the linear solver of inexact Newton
     * according to "Choice 2" in [Eisenstat, Walker, 1996]:
     *
     *   rtol_max * (residual / prev residual)^exponent
     */
    lin_rtol = lin_rtol_max * pow (res_norm/res_norm_prev, exponent);

    /* set progressive reduction factor:
     *
     *   (exp(-(#iter / prog #iter)^2) + eps/rtol) /
     *   (1 + eps/rtol)
     */
    if (prog_reduction_active) {
      const double        iter2 = (double) iter * iter;
      const double        prog_n_iter2 = (double)
                            opt->adaptive_lin_rtol_progressive_n_iter *
                            opt->adaptive_lin_rtol_progressive_n_iter;
      const double        eps = SC_1000_EPS;
      const double        rt = SC_MIN (lin_rtol, lin_rtol_max);

      prog_reduction = (exp (-iter2/prog_n_iter2) + eps/rt) / (1.0 + eps/rt);
    }

    /* set min rtol, which serves as a saveguard to avoid oversolving
     *
     *   rtol_max * ( (prev reduction + residual / prev residual)/2 )^exponent
     */
    if (lin_rtol_min_active) {
      const double        reduction_prev = step->lin_res_norm_reduction;
      double              reduction_avg;

      reduction_avg = 0.5 * ( SC_MIN (reduction_prev, 1.0) +
                              res_norm/res_norm_prev );
      lin_rtol_min = lin_rtol_max * pow (reduction_avg, exponent);
    }

    RHEA_GLOBAL_INFOF (
        "%s: Newton iter %d, candidate for adaptive linear rtol %.3e, "
        "max %.3e, min %.3e, min threshold %.3e, progressive reduction %.3e\n",
        this_fn_name, iter, lin_rtol, lin_rtol_max, lin_rtol_min,
        lin_rtol_min_thresh, prog_reduction);

    /* deactivate min rtol below a given threshold */
    if (lin_rtol_min <= lin_rtol_min_thresh) {
      lin_rtol_min = 0.0;
    }

    /* set adaptive rtol */
    lin_rtol = SC_MIN (lin_rtol, lin_rtol_max);
    lin_rtol *= prog_reduction;
    lin_rtol = SC_MAX (lin_rtol_min, lin_rtol);

    RHEA_GLOBAL_INFOF ("%s: Newton iter %d, set adaptive linear rtol %.3e\n",
                       this_fn_name, iter, lin_rtol);
  }

  /* store the relative tolerance */
  step->lin_res_norm_rtol = lin_rtol;
}

/**
 * Checks whether to continue the loop of multiple runs of the linear solver.
 */
static int
rhea_newton_solve_step_continue_linear_solver_loop (
                                          const int lin_iter_count,
                                          const int lin_iter_max,
                                          const double lin_res_norm_rtol,
                                          const double lin_res_norm_reduction)
{
  const double        thresh = 0.9;
  int                 initial_run, iterations_available, missing_accuracy;

  initial_run = (0 < lin_iter_count);
  iterations_available = (
    (double) lin_iter_count < thres * (double) lin_iter_max
  );
  missing_accuracy = (
    lin_res_norm_rtol / lin_res_norm_reduction < thresh
  );

  return (initial_run || (iterations_available && missing_accuracy));
}

/**
 * Solves for Newton step.
 */
static void
rhea_newton_solve_step (
    /* out: */ step,
    /* in:  */ residual, nl_problem, opt);
{
  const char         *this_fn_name = "rhea_newton_solve_step";
  const int           lin_iter_max = ;
  int                 lin_iter_count = 0;
  const double        lin_res_norm_rtol = ;
  double              lin_res_norm_reduction = 1.0;

  ymir_vec_t         *res_vec = residual->residual_vec;
  ymir_vec_t         *lin_res_vec = ymir_vec_template (res_vec);
  ymir_vec_t         *step_vec = step->step_vec;

  double              lin_res_norm_init, lin_res_norm;
  double              rt, at;
  int                 it_max, it_count = 0;
  int                 gmres_n_vecs;

  /* compute l^2 norm of residual before running the linear solver */
  lin_res_norm_init = ymir_vec_norm (res_vec);

  /* solve linear system */
  while (
    rhea_newton_solve_step_continue_linear_solver_loop (
        lin_iter_count, lin_iter_max,
        lin_res_norm_rtol, lin_res_norm_reduction)
  ) {
    const int           nonzero_initial_guess = (0 < lin_iter_count);

    //TODO:
    RHEA_GLOBAL_INFOF ("Krylov step %i: l^2 norm of initial Krylov residual "
                       "%1.3e\n", k, lin_res_norm_init);

    /* run linear solver (use zero initial guess) */
    it_count = -1;
    ymir_stokes_pc_solve (
        res_vec, step_vec, nl_problem->stokes_pc, nonzero_initial_guess,
        rt, at, it_max, gmres_n_vecs, &it_count);
    RHEA_ASSERT (-1 < it_count);
    lin_iter_count += it_count;

    /* compute l^2 norm of the residual of the linear system */
    ymir_nlstokes_op_apply (step_vec, lin_res_vec, nl_problem->stokes_op);
    ymir_vec_add (-1.0, res_vec, lin_res_vec);
    lin_res_norm = ymir_vec_norm (lin_res_vec);

    //TODO:
    RHEA_GLOBAL_INFOF ("Krylov step %i: l^2 norm of Krylov residual %1.3e\n",
                       k, lin_res_norm);

    /* update residual reduction */
    lin_res_norm_reduction = lin_res_norm / lin_res_norm_init;

    //TODO:
    RHEA_GLOBAL_INFOF ("Krylov step %i: %i Krylov iter for rtol %1.3e "
                       "(prescribed rtol %1.3e)\n",
                       k, lin_iter_count, lin_res_norm_reduction,
                       lin_res_norm_rtol);

    /* update Krylov rtol and maxiter */
    rt *= lin_res_norm_rtol / lin_res_norm_reduction;
    it_max = lin_iter_max - lin_iter_count;
    k++;
  }

  /* destroy */
  ymir_vec_destroy (lin_res_vec);

  /* set residual reduction and return the number of Krylov iterations */
  step->lin_res_iter_count; = lin_res_iter_count;
  step->lin_res_norm_reduction = lin_res_norm_reduction;
}

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
