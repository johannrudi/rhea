/*
 */

#include <rhea_newton_check.h>
#include <rhea_base.h>

void
rhea_newton_check_gradient (ymir_vec_t *solution,
                            rhea_newton_problem_t *nl_problem)
{
  const char         *this_fn_name = "rhea_newton_check_gradient";
  ymir_vec_t         *neg_gradient_vec =
                        rhea_newton_problem_get_neg_gradient_vec (nl_problem);
  ymir_vec_t         *step_vec = rhea_newton_problem_get_step_vec (nl_problem);

  ymir_vec_t         *sol_vec, *dir_vec, *perturb_vec;
  double              grad_dir_ref, grad_dir_chk;
  double              obj_val, perturb_obj_val;
  double              abs_error, rel_error;
  int                 n;

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (NULL != rhea_newton_problem_get_neg_gradient_vec (nl_problem));
  RHEA_ASSERT (NULL != rhea_newton_problem_get_step_vec (nl_problem));

  /* return if check cannot be executed */
  if (!rhea_newton_problem_evaluate_objective_exists (nl_problem)) {
    RHEA_GLOBAL_VERBOSEF (
        "%s: Cannot execute because objective functional is not provided.\n",
        this_fn_name);
    return;
  }

  /* create work vectors */
  sol_vec = ymir_vec_template (step_vec);
  dir_vec = ymir_vec_template (step_vec);
  perturb_vec = ymir_vec_template (step_vec);

  /* set position and direction vectors */
  if (solution == NULL) {
    ymir_vec_set_zero (sol_vec);
  }
  else {
    ymir_vec_copy (solution, sol_vec);
  }
#ifdef YMIR_WITH_PETSC
  ymir_petsc_vec_set_random (dir_vec, YMIR_MESH_PETSCLAYOUT_NONE);
#else
  RHEA_ABORT_NOT_REACHED ();
#endif

  /* compute reference derivative */
  grad_dir_ref = ymir_vec_innerprod (neg_gradient_vec, dir_vec);
  grad_dir_ref *= -1.0;

  /* compare with finite difference derivative */
  obj_val = rhea_newton_problem_evaluate_objective (sol_vec, nl_problem);
  for (n = 2; n <= 6; n++) {
    const double        eps = pow (10.0, (double) -n);

    ymir_vec_copy (sol_vec, perturb_vec);
    ymir_vec_add (eps, dir_vec, perturb_vec);
    perturb_obj_val = rhea_newton_problem_evaluate_objective (perturb_vec,
                                                              nl_problem);
    grad_dir_chk = (perturb_obj_val - obj_val) / eps;

    abs_error = fabs (grad_dir_ref - grad_dir_chk);
    rel_error = abs_error / fabs (grad_dir_ref);

    RHEA_GLOBAL_INFOF ("%s: eps %.1e, error abs %.1e, rel %.1e "
                       "(grad^T*dir ref %.12e, chk %.12e)\n",
                       this_fn_name, eps, abs_error, rel_error,
                       grad_dir_ref, grad_dir_chk);
  }

  /* destroy */
  ymir_vec_destroy (sol_vec);
  ymir_vec_destroy (dir_vec);
  ymir_vec_destroy (perturb_vec);

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
}

void
rhea_newton_check_hessian (ymir_vec_t *solution,
                           rhea_newton_problem_t *nl_problem)
{
  const char         *this_fn_name = "rhea_newton_check_hessian";
  const int           check_gradient =
                        rhea_newton_problem_get_check_gradient (nl_problem);
  ymir_vec_t         *neg_gradient_vec =
                        rhea_newton_problem_get_neg_gradient_vec (nl_problem);
  ymir_vec_t         *step_vec = rhea_newton_problem_get_step_vec (nl_problem);

  ymir_vec_t         *sol_vec, *dir_vec, *perturb_vec;
  ymir_vec_t         *H_dir_ref, *H_dir_chk;
  double              abs_error, rel_error;
  int                 n;

  RHEA_GLOBAL_VERBOSEF ("Into %s\n", this_fn_name);

  /* check input */
  RHEA_ASSERT (NULL != solution);
  RHEA_ASSERT (NULL != rhea_newton_problem_get_neg_gradient_vec (nl_problem));
  RHEA_ASSERT (NULL != rhea_newton_problem_get_step_vec (nl_problem));
  RHEA_ASSERT (rhea_newton_problem_compute_neg_gradient_exists (nl_problem));

  /* return if check cannot be executed */
  if (!rhea_newton_problem_apply_hessian_exists (nl_problem)) {
    RHEA_GLOBAL_VERBOSEF (
        "%s: Cannot execute because Hessian-apply function is not provided.\n",
        this_fn_name);
    return;
  }

  /* create work vectors */
  sol_vec = solution;
  dir_vec = ymir_vec_template (step_vec);
  perturb_vec = ymir_vec_template (step_vec);
  H_dir_ref = ymir_vec_template (step_vec);
  H_dir_chk = ymir_vec_template (step_vec);

  /* set position and direction vectors */
#ifdef YMIR_WITH_PETSC
  ymir_petsc_vec_set_random (dir_vec, YMIR_MESH_PETSCLAYOUT_NONE);
#else
  RHEA_ABORT_NOT_REACHED ();
#endif

  /* compute reference Hessian-vector apply */
  rhea_newton_problem_apply_hessian (H_dir_ref, dir_vec, nl_problem);

  /* compare with finite difference Hessian */
  rhea_newton_problem_set_check_gradient (0, nl_problem);
  for (n = 2; n <= 6; n++) {
    const double        eps = pow (10.0, (double) -n);

    ymir_vec_copy (sol_vec, perturb_vec);
    ymir_vec_add (eps, dir_vec, perturb_vec);
    rhea_newton_problem_compute_neg_gradient (H_dir_chk, perturb_vec,
                                              nl_problem);
    ymir_vec_scale (-1.0, H_dir_chk);
    ymir_vec_add (1.0, neg_gradient_vec, H_dir_chk);
    ymir_vec_scale (1.0/eps, H_dir_chk);

    ymir_vec_add (-1.0, H_dir_ref, H_dir_chk);
    abs_error = ymir_vec_norm (H_dir_chk);
    rel_error = abs_error / ymir_vec_norm (H_dir_ref);

    RHEA_GLOBAL_INFOF ("%s: eps %.1e, error abs %.1e, rel %.1e\n",
                       this_fn_name, eps, abs_error, rel_error);
  }
  rhea_newton_problem_set_check_gradient (check_gradient, nl_problem);

  /* destroy */
  ymir_vec_destroy (dir_vec);
  ymir_vec_destroy (perturb_vec);
  ymir_vec_destroy (H_dir_ref);
  ymir_vec_destroy (H_dir_chk);

  RHEA_GLOBAL_VERBOSEF ("Done %s\n", this_fn_name);
}
