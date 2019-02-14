/*
 */

#include <rhea_newton_check.h>
#include <rhea_base.h>

#define RHEA_NEWTON_CHECK_N_TRIALS (5)

void
rhea_newton_check_gradient (ymir_vec_t *solution,
                            rhea_newton_problem_t *nl_problem)
{
  ymir_vec_t         *neg_gradient_vec =
                        rhea_newton_problem_get_neg_gradient_vec (nl_problem);
  ymir_vec_t         *step_vec = rhea_newton_problem_get_step_vec (nl_problem);

  ymir_vec_t         *sol_vec, *dir_vec, *perturb_vec;
  double              grad_dir_ref, grad_dir_chk;
  double              obj_val, perturb_obj_val;
  double              abs_error, rel_error;
  sc_dmatrix_t       *result = sc_dmatrix_new (RHEA_NEWTON_CHECK_N_TRIALS, 5);
  const int           exp_incr = 2;
  const int           exp_min = 2;
  int                 n;

  /* check input */
  RHEA_ASSERT (NULL != rhea_newton_problem_get_neg_gradient_vec (nl_problem));
  RHEA_ASSERT (NULL != rhea_newton_problem_get_step_vec (nl_problem));

  /* exit if check cannot be executed */
  if (!rhea_newton_problem_evaluate_objective_exists (nl_problem)) {
    RHEA_GLOBAL_VERBOSEF_FN_TAG (
        __func__, "stop_reason=\"%s\"",
        "Cannot execute because objective functional is not provided");
    return;
  }

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);

  /* create work vectors */
  sol_vec = ymir_vec_template (step_vec);
  dir_vec = ymir_vec_template (step_vec);
  perturb_vec = ymir_vec_template (step_vec);

  /* set position and direction vectors (possibly adjust scale of direction) */
  if (solution == NULL) {
    ymir_vec_set_zero (sol_vec);
    ymir_vec_set_random (dir_vec);
  }
  else {
    ymir_vec_copy (solution, sol_vec);
    ymir_vec_set_random (dir_vec);
    ymir_vec_multiply_in (sol_vec, dir_vec);
  }

  /* compute reference derivative */
  grad_dir_ref = ymir_vec_innerprod (neg_gradient_vec, dir_vec);
  grad_dir_ref *= -1.0;

  /* compare with finite difference derivative */
  rhea_newton_problem_update_operator (sol_vec, nl_problem);
  obj_val = rhea_newton_problem_evaluate_objective (sol_vec, nl_problem);
  for (n = 0; n < RHEA_NEWTON_CHECK_N_TRIALS; n++) {
    const int           exp = exp_min + exp_incr * n;
    const double        eps = pow (10.0, (double) -exp);

    ymir_vec_copy (sol_vec, perturb_vec);
    ymir_vec_add (eps, dir_vec, perturb_vec);
    rhea_newton_problem_update_operator (perturb_vec, nl_problem);
    perturb_obj_val = rhea_newton_problem_evaluate_objective (perturb_vec,
                                                              nl_problem);
    grad_dir_chk = (perturb_obj_val - obj_val) / eps;

    abs_error = fabs (grad_dir_ref - grad_dir_chk);
    rel_error = abs_error / fabs (grad_dir_ref);

    result->e[n][0] = eps;
    result->e[n][1] = abs_error;
    result->e[n][2] = rel_error;
    result->e[n][3] = grad_dir_ref;
    result->e[n][4] = grad_dir_chk;
  }

  /* restore */
  rhea_newton_problem_update_operator (solution, nl_problem);

  /* print results */
  RHEA_GLOBAL_INFO ("========================================\n");
  RHEA_GLOBAL_INFOF ("%s\n", __func__);
  RHEA_GLOBAL_INFO ("----------------------------------------\n");
  for (n = 0; n < RHEA_NEWTON_CHECK_N_TRIALS; n++) {
    RHEA_GLOBAL_INFOF (
        "eps=%.1e ; error abs=%.1e, rel=%.1e ; "
        "grad^T*dir ref=%.12e, chk=%.12e\n",
        result->e[n][0], result->e[n][1], result->e[n][2],
        result->e[n][3], result->e[n][4]);
  }
  RHEA_GLOBAL_INFO ("========================================\n");

  /* destroy */
  ymir_vec_destroy (sol_vec);
  ymir_vec_destroy (dir_vec);
  ymir_vec_destroy (perturb_vec);
  sc_dmatrix_destroy (result);

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

void
rhea_newton_check_hessian (ymir_vec_t *solution,
                           rhea_newton_problem_t *nl_problem)
{
  const int           check_gradient =
                        rhea_newton_problem_get_check_gradient (nl_problem);
  ymir_vec_t         *neg_gradient_vec =
                        rhea_newton_problem_get_neg_gradient_vec (nl_problem);
  ymir_vec_t         *step_vec = rhea_newton_problem_get_step_vec (nl_problem);

  ymir_vec_t         *sol_vec, *dir_vec, *perturb_vec;
  ymir_vec_t         *neg_grad_vec, *H_dir_ref, *H_dir_chk;
  double              abs_error, rel_error, norm_ref;
  sc_dmatrix_t       *result = sc_dmatrix_new (RHEA_NEWTON_CHECK_N_TRIALS, 3);
  const int           exp_incr = 2;
  const int           exp_min = 2;
  int                 n;

  /* check input */
  RHEA_ASSERT (NULL != solution);
  RHEA_ASSERT (NULL != rhea_newton_problem_get_neg_gradient_vec (nl_problem));
  RHEA_ASSERT (NULL != rhea_newton_problem_get_step_vec (nl_problem));
  RHEA_ASSERT (rhea_newton_problem_compute_neg_gradient_exists (nl_problem));

  /* exit if check cannot be executed */
  if (!rhea_newton_problem_apply_hessian_exists (nl_problem)) {
    RHEA_GLOBAL_VERBOSEF_FN_TAG (
        __func__, "stop_reason=\"%s\"",
        "Cannot execute because Hessian-apply function is not provided");
    return;
  }

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);

  /* create work vectors */
  sol_vec = solution;
  dir_vec = ymir_vec_template (step_vec);
  perturb_vec = ymir_vec_template (step_vec);
  neg_grad_vec = ymir_vec_template (neg_gradient_vec);
  H_dir_ref = ymir_vec_template (neg_gradient_vec);
  H_dir_chk = ymir_vec_template (neg_gradient_vec);

  /* set direction vector (scale similarly as solution vector) */
  ymir_vec_set_random (dir_vec);
  ymir_vec_multiply_in (sol_vec, dir_vec);

  /* compute reference Hessian-vector apply */
  rhea_newton_problem_apply_hessian (H_dir_ref, dir_vec, nl_problem);
  norm_ref = ymir_vec_norm (H_dir_ref);

  /* setup finite difference Hessian */
  rhea_newton_problem_set_check_gradient (0, nl_problem);
  rhea_newton_problem_update_operator (sol_vec, nl_problem);
  rhea_newton_problem_compute_neg_gradient (neg_grad_vec, sol_vec,
                                            nl_problem);

  /* compare result of reference Hessian-vector application with
   * finite difference Hessian */
  for (n = 0; n < RHEA_NEWTON_CHECK_N_TRIALS; n++) {
    const int           exp = exp_min + exp_incr * n;
    const double        eps = pow (10.0, (double) -exp);

    ymir_vec_copy (sol_vec, perturb_vec);
    ymir_vec_add (eps, dir_vec, perturb_vec);
    rhea_newton_problem_update_operator (perturb_vec, nl_problem);
    rhea_newton_problem_compute_neg_gradient (H_dir_chk, perturb_vec,
                                              nl_problem);
    ymir_vec_scale (-1.0, H_dir_chk);
    ymir_vec_add (1.0, neg_grad_vec, H_dir_chk);
    ymir_vec_scale (1.0/eps, H_dir_chk);

    ymir_vec_add (-1.0, H_dir_ref, H_dir_chk);
    abs_error = ymir_vec_norm (H_dir_chk);
    rel_error = abs_error / norm_ref;

    result->e[n][0] = eps;
    result->e[n][1] = abs_error;
    result->e[n][2] = rel_error;
  }

  /* restore */
  rhea_newton_problem_set_check_gradient (check_gradient, nl_problem);
  rhea_newton_problem_update_operator (solution, nl_problem);

  /* print results */
  RHEA_GLOBAL_INFO ("========================================\n");
  RHEA_GLOBAL_INFOF ("%s\n", __func__);
  RHEA_GLOBAL_INFO ("----------------------------------------\n");
  for (n = 0; n < RHEA_NEWTON_CHECK_N_TRIALS; n++) {
    RHEA_GLOBAL_INFOF (
        "eps=%.1e ; error abs=%.1e, rel=%.1e\n",
        result->e[n][0], result->e[n][1], result->e[n][2]);
  }
  RHEA_GLOBAL_INFO ("========================================\n");

  /* destroy */
  ymir_vec_destroy (dir_vec);
  ymir_vec_destroy (perturb_vec);
  ymir_vec_destroy (neg_grad_vec);
  ymir_vec_destroy (H_dir_ref);
  ymir_vec_destroy (H_dir_chk);
  sc_dmatrix_destroy (result);

  RHEA_GLOBAL_INFO_FN_END (__func__);
}
