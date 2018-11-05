/*
 */

#include <rhea_newton_check.h>
#include <rhea_base.h>

void
rhea_newton_check_gradient (ymir_vec_t *solution,
                            rhea_newton_problem_t *nl_problem)
{
  ymir_vec_t         *neg_gradient_vec =
                        rhea_newton_problem_get_neg_gradient_vec (nl_problem);
  ymir_vec_t         *step_vec = rhea_newton_problem_get_step_vec (nl_problem);

RHEA_GLOBAL_PRODUCTIONF ("Hello%d\n",10);
RHEA_ASSERT (step_vec != NULL);
RHEA_ASSERT (step_vec->meshnum = YMIR_VEC_MESHFREE);
RHEA_GLOBAL_PRODUCTIONF ("m=%d n=%d\n",step_vec->meshfree->m, step_vec->meshfree->n);
RHEA_ASSERT (step_vec->meshfree->m > 0 && step_vec->meshfree->n > 0);
RHEA_GLOBAL_PRODUCTIONF ("step_vec=%f\n",step_vec->meshfree->e[0][0]);

  ymir_vec_t         *sol_vec, *dir_vec, *perturb_vec;
  double              grad_dir_ref, grad_dir_chk;
  double              obj_val, perturb_obj_val;
  double              abs_error, rel_error;
  int                 n;

  /* check input */
  RHEA_ASSERT (NULL != rhea_newton_problem_get_neg_gradient_vec (nl_problem));
  RHEA_ASSERT (NULL != rhea_newton_problem_get_step_vec (nl_problem));

  /* return if check cannot be executed */
  if (!rhea_newton_problem_evaluate_objective_exists (nl_problem)) {
    RHEA_GLOBAL_VERBOSEF (
        "%s Cannot execute because objective functional is not provided.\n",
        __func__);
    return;
  }

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);

  /* create work vectors */
  sol_vec = ymir_vec_new_meshfree (1);
  dir_vec = ymir_vec_new_meshfree (1);
  perturb_vec = ymir_vec_new_meshfree (1);
//  sol_vec = ymir_vec_template (step_vec);
//RHEA_GLOBAL_PRODUCTIONF ("Hello%d\n",0);
//  dir_vec = ymir_vec_template (step_vec);
//RHEA_GLOBAL_PRODUCTIONF ("Hello%d\n",1);
//  perturb_vec = ymir_vec_template (step_vec);
//RHEA_GLOBAL_PRODUCTIONF ("Hello%d\n",2);
//exit clearly here
  /* set position and direction vectors (possibly adjust scale of direction) */
  if (solution == NULL) {
    ymir_vec_set_zero (sol_vec);
    ymir_vec_set_random (dir_vec);
  }
  else {
RHEA_GLOBAL_PRODUCTIONF ("solution!=NULL%d\n",1);
//exit clearly here
    ymir_vec_copy (solution, sol_vec); //problem is here
    ymir_vec_set_random (dir_vec);
RHEA_GLOBAL_PRODUCTIONF ("solution!=NULL%d\n",2);
    ymir_vec_multiply_in (sol_vec, dir_vec);
RHEA_GLOBAL_PRODUCTIONF ("solution!=NULL%d\n",3);
  }
  /* compute reference derivative */
  grad_dir_ref = ymir_vec_innerprod (neg_gradient_vec, dir_vec);
  grad_dir_ref *= -1.0;

  /* compare with finite difference derivative */
  rhea_newton_problem_update_operator (sol_vec, nl_problem);
  obj_val = rhea_newton_problem_evaluate_objective (sol_vec, nl_problem);
  for (n = 2; n <= 10; n += 2) {
    const double        eps = pow (10.0, (double) -n);

    ymir_vec_copy (sol_vec, perturb_vec);
    ymir_vec_add (eps, dir_vec, perturb_vec);
    rhea_newton_problem_update_operator (perturb_vec, nl_problem);
    perturb_obj_val = rhea_newton_problem_evaluate_objective (perturb_vec,
                                                              nl_problem);
    grad_dir_chk = (perturb_obj_val - obj_val) / eps;

    abs_error = fabs (grad_dir_ref - grad_dir_chk);
    rel_error = abs_error / fabs (grad_dir_ref);

    RHEA_GLOBAL_INFOF ("<%s_result eps=%.1e, error abs=%.1e, rel=%.1e "
                       "(grad^T*dir ref=%.12e, chk=%.12e) />\n",
                       __func__, eps, abs_error, rel_error,
                       grad_dir_ref, grad_dir_chk);
  }

  /* restore */
  rhea_newton_problem_update_operator (solution, nl_problem);

  /* destroy */
  ymir_vec_destroy (sol_vec);
  ymir_vec_destroy (dir_vec);
  ymir_vec_destroy (perturb_vec);

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
  int                 n;

  /* check input */
  RHEA_ASSERT (NULL != solution);
  RHEA_ASSERT (NULL != rhea_newton_problem_get_neg_gradient_vec (nl_problem));
  RHEA_ASSERT (NULL != rhea_newton_problem_get_step_vec (nl_problem));
  RHEA_ASSERT (rhea_newton_problem_compute_neg_gradient_exists (nl_problem));

  /* return if check cannot be executed */
  if (!rhea_newton_problem_apply_hessian_exists (nl_problem)) {
    RHEA_GLOBAL_VERBOSEF (
        "%s Cannot execute because Hessian-apply function is not provided.\n",
        __func__);
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
  for (n = 2; n <= 10; n += 2) {
    const double        eps = pow (10.0, (double) -n);

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

    RHEA_GLOBAL_INFOF ("<%s_result eps=%.1e, error abs=%.1e, rel=%.1e />\n",
                       __func__, eps, abs_error, rel_error);
  }

  /* restore */
  rhea_newton_problem_set_check_gradient (check_gradient, nl_problem);
  rhea_newton_problem_update_operator (solution, nl_problem);

  /* destroy */
  ymir_vec_destroy (dir_vec);
  ymir_vec_destroy (perturb_vec);
  ymir_vec_destroy (neg_grad_vec);
  ymir_vec_destroy (H_dir_ref);
  ymir_vec_destroy (H_dir_chk);

  RHEA_GLOBAL_INFO_FN_END (__func__);
}
