/*
 */

#include <rhea_newton_check.h>
#include <rhea_base.h>
#include <ymir_vec_getset.h>
#include <ymir_mass_vec.h>
#include <ymir_pressure_vec.h>

/* number of reductions of the perturbation for finite differences */
#define RHEA_NEWTON_CHECK_N_TRIALS (6)

static void
rhea_newton_check_scale_cvec (ymir_vec_t *scale, ymir_vec_t *source)
{
  ymir_face_mesh_t   *face_mesh = ymir_vec_get_face_mesh (source);
  const ymir_locidx_t n_nodes = face_mesh->Ncn;
  const int           n_fields = source->ncfields;
  ymir_locidx_t       nodeid;
  int                 fieldid;

  for (nodeid = 0; nodeid < n_nodes; nodeid++) {
    const double       *src = ymir_cvec_index (source, nodeid, 0);
    double             *scl = ymir_cvec_index (scale, nodeid, 0);
    double              magn = 0.0;

    for (fieldid = 0; fieldid < n_fields; fieldid++) {
      magn += src[fieldid] * src[fieldid];
    }
    magn = sqrt (magn);

    if (1.0 < magn) {
      for (fieldid = 0; fieldid < n_fields; fieldid++) {
        scl[fieldid] *= magn;
      }
    }
  }
}

static void
rhea_newton_check_scale_dvec (ymir_vec_t *scale, ymir_vec_t *source)
{
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (source);
  const ymir_locidx_t n_elements = mesh->cnodes->K;
  const int           n_nodes = ymir_np (mesh->cnodes->N);
  const int           n_fields = source->ndfields;
  ymir_locidx_t       elid;
  int                 nodeid, fieldid;

  for (elid = 0; elid < n_elements; elid++) {
    for (nodeid = 0; nodeid < n_nodes; nodeid++) {
      const double       *src = ymir_dvec_index (source, elid, nodeid, 0);
      double             *scl = ymir_dvec_index (scale, elid, nodeid, 0);
      double              magn = 0.0;

      for (fieldid = 0; fieldid < n_fields; fieldid++) {
        magn += src[fieldid] * src[fieldid];
      }
      magn = sqrt (magn);

      if (1.0 < magn) {
        for (fieldid = 0; fieldid < n_fields; fieldid++) {
          scl[fieldid] *= magn;
        }
      }
    }
  }
}

static void
rhea_newton_check_scale_evec (ymir_vec_t *scale, ymir_vec_t *source)
{
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (source);
  const ymir_locidx_t n_elements = mesh->cnodes->K;
  const int           n_nodes = source->nefields;
  ymir_locidx_t       elid;
  int                 nodeid;
  double              magn;

  for (elid = 0; elid < n_elements; elid++) {
    const double       *src = ymir_evec_index (source, elid, 0);
    double             *scl = ymir_evec_index (scale, elid, 0);

    for (nodeid = 0; nodeid < n_nodes; nodeid++) {
      magn = fabs (src[nodeid]);
      if (1.0 < magn) {
        scl[nodeid] *= magn;
      }
    }
  }
}

static void
rhea_newton_check_scale_meshfree (ymir_vec_t *scale, ymir_vec_t *source)
{
  const int           n_entries = source->n_meshfree;
  int                 i;
  const double       *src = source->meshfree->e[0];
  double             *scl = scale->meshfree->e[0];
  double              magn;

  for (i = 0; i < n_entries; i++) {
    magn = fabs (src[i]);
    if (1.0 < magn) {
      scl[i] *= magn;
    }
    scl[i] *= 0.5;
  }
}

static void
rhea_newton_check_set_perturb_vec (ymir_vec_t *perturb_vec, ymir_vec_t *sol_vec,
                                   sc_MPI_Comm mpicomm)
{
  /* set to random values in [0,1] */
  ymir_vec_set_random_uniform (perturb_vec);

  /* communicate values in case of meshfree vectors */
  if (perturb_vec->n_meshfree && mpicomm != MPI_COMM_NULL) {
    ymir_vec_meshfree_sync (perturb_vec, 0 /* mpirank_master */, mpicomm);
  }
  else if (perturb_vec->n_meshfree) {
    RHEA_GLOBAL_INFOF (
        "%s: Warning: meshfree direction is not equal across mpiranks; "
        "provide MPI communicator to Newton problem.\n", __func__);
  }

  /* exit if no solution vector given */
  if (sol_vec == NULL) {
    return;
  }

  /* scale w.r.t. magnitude of solution vector */
  if (sol_vec->ncfields) {
    rhea_newton_check_scale_cvec (perturb_vec, sol_vec);
  }
  if (sol_vec->ndfields) {
    rhea_newton_check_scale_dvec (perturb_vec, sol_vec);
  }
  if (sol_vec->nefields) {
    rhea_newton_check_scale_evec (perturb_vec, sol_vec);
  }
  if (sol_vec->n_meshfree) {
    rhea_newton_check_scale_meshfree (perturb_vec, sol_vec);
  }
}

void
rhea_newton_check_gradient (ymir_vec_t *solution,
                            ymir_vec_t *neg_gradient,
                            int n_trials,
                            ymir_vec_t **perturb_vec,
                            int n_perturb_vecs,
                            const int iter,
                            rhea_newton_problem_t *nl_problem)
{
  ymir_vec_t         *sol_vec, *perturb_sol_vec;
  const int           create_perturb_vec = (n_perturb_vecs <= 0);
  double              grad_dir_ref, grad_dir_fd;
  double              obj_val, perturb_obj_val;
  int                 obj_err_exists;
  double              obj_err_mean, obj_err_stddev;
  double              abs_error, rel_error;
  sc_dmatrix_t       *result;
  const int           exp_incr = 2;
  const int           exp_min = 0;
  int                 vidx, tidx;

  /* check input */
  RHEA_ASSERT (NULL != neg_gradient);
  RHEA_ASSERT (NULL != rhea_newton_problem_get_step_vec (nl_problem));

  /* exit if nothing to do */
  if (n_trials <= 0) {
    return;
  }
  /* if needs setting default #trials */
  if (1 == n_trials) {
    n_trials = RHEA_NEWTON_CHECK_N_TRIALS;
  }

  /* exit if check cannot be executed */
  if (!rhea_newton_problem_evaluate_objective_exists (nl_problem)) {
    RHEA_GLOBAL_VERBOSEF_FN_TAG (
        __func__, "stop_reason=\"%s\"",
        "Cannot execute because objective functional is not provided");
    return;
  }

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "newton_iter=%i", iter);

  /* create work vectors */
  if (solution == NULL) {
    sol_vec = ymir_vec_template (rhea_newton_problem_get_step_vec (nl_problem));
    ymir_vec_set_zero (sol_vec);
  }
  else {
    sol_vec = ymir_vec_template (solution);
    ymir_vec_copy (solution, sol_vec);
  }
  perturb_sol_vec = ymir_vec_template (sol_vec);
  if (create_perturb_vec) {
    n_perturb_vecs = 1;
    perturb_vec = RHEA_ALLOC (ymir_vec_t *, 1);
    perturb_vec[0] = ymir_vec_template (sol_vec);
    rhea_newton_check_set_perturb_vec (
        perturb_vec[0], solution,
        rhea_newton_problem_get_mpicomm (nl_problem));
  }
  result = sc_dmatrix_new (n_trials, 7);

  /* set up the finite difference derivative */
  obj_val = rhea_newton_problem_evaluate_objective (solution, nl_problem,
                                                    NULL /* components */);
  obj_err_exists = rhea_newton_problem_evaluate_objective_err (
      &obj_err_mean, &obj_err_stddev, nl_problem);

  for (vidx = 0; vidx < n_perturb_vecs; vidx++) {
    /* compute the reference directional derivative */
    grad_dir_ref = ymir_vec_innerprod (neg_gradient, perturb_vec[vidx]);
    grad_dir_ref *= -1.0;

    /* compare reference with finite difference derivative */
    for (tidx = 0; tidx < n_trials; tidx++) {
      const int           exp = exp_min + exp_incr * tidx;
      const double        eps = pow (10.0, (double) -exp);

      ymir_vec_copy (sol_vec, perturb_sol_vec);
      ymir_vec_add (eps, perturb_vec[vidx], perturb_sol_vec);
      rhea_newton_problem_update_operator (perturb_sol_vec, nl_problem);
      perturb_obj_val = rhea_newton_problem_evaluate_objective (
          perturb_sol_vec, nl_problem, NULL /* components */);
      obj_err_exists = rhea_newton_problem_evaluate_objective_err (
          &obj_err_mean, &obj_err_stddev, nl_problem);
      grad_dir_fd = (perturb_obj_val - obj_val) / eps;

      abs_error = fabs (grad_dir_ref - grad_dir_fd);
      rel_error = abs_error / fabs (grad_dir_ref);

      result->e[tidx][0] = eps;
      result->e[tidx][1] = abs_error;
      result->e[tidx][2] = rel_error;
      result->e[tidx][3] = grad_dir_ref;
      result->e[tidx][4] = grad_dir_fd;
      result->e[tidx][5] = obj_err_mean;
      result->e[tidx][6] = obj_err_stddev;
    }

    /* print results */
    RHEA_GLOBAL_INFO ("========================================\n");
    RHEA_GLOBAL_INFOF ("%s, newton_iter=%i, perturb_vec_idx=%i\n",
                       __func__, iter, vidx);
    RHEA_GLOBAL_INFO ("----------------------------------------\n");
    for (tidx = 0; tidx < n_trials; tidx++) {
      if (obj_err_exists) {
        RHEA_GLOBAL_INFOF (
            "F.D. eps=%.1e ; error abs=%.1e, rel=%.1e ; "
            "grad^T*dir ref=%.12e, F.D.=%.12e ; "
            "obj eval err (%.1e,%.1e)\n",
            result->e[tidx][0], result->e[tidx][1], result->e[tidx][2],
            result->e[tidx][3], result->e[tidx][4],
            result->e[tidx][5], result->e[tidx][6]);
      }
      else {
        RHEA_GLOBAL_INFOF (
            "F.D. eps=%.1e ; error abs=%.1e, rel=%.1e ; "
            "grad^T*dir ref=%.12e, F.D.=%.12e\n",
            result->e[tidx][0], result->e[tidx][1], result->e[tidx][2],
            result->e[tidx][3], result->e[tidx][4]);
      }
    }
    RHEA_GLOBAL_INFO ("========================================\n");
  }

  /* restore operator */
  rhea_newton_problem_update_operator (solution, nl_problem);

  /* destroy */
  ymir_vec_destroy (sol_vec);
  ymir_vec_destroy (perturb_sol_vec);
  if (create_perturb_vec) {
    ymir_vec_destroy (perturb_vec[0]);
    RHEA_FREE (perturb_vec);
  }
  sc_dmatrix_destroy (result);

  RHEA_GLOBAL_INFOF_FN_END (__func__, "newton_iter=%i", iter);
}

static void
rhea_newton_check_hessian_compute_error (double *abs_error,
                                         double *rel_error,
                                         ymir_vec_t *hess_dir_ref,
                                         ymir_vec_t *hess_dir_fd)
{
  int                 n_error_components = 0;

  /* check input */
  RHEA_ASSERT (hess_dir_ref->ncfields == hess_dir_fd->ncfields);
  RHEA_ASSERT (hess_dir_ref->ndfields == hess_dir_fd->ndfields);
  RHEA_ASSERT (hess_dir_ref->nefields == hess_dir_fd->nefields);
  RHEA_ASSERT (hess_dir_ref->n_meshfree == hess_dir_fd->n_meshfree);

  /* initialize output */
  *abs_error = 0.0;
  *rel_error = 0.0;

  /* compute the L2-norm of the error */
  if (hess_dir_ref->ncfields ||
      hess_dir_ref->ndfields ||
      hess_dir_ref->nefields) {
    ymir_vec_t         *tmp = ymir_vec_template (hess_dir_ref);
    ymir_vec_t         *inv_mass_lump = ymir_vec_template (hess_dir_ref);
    double              norm_ref, abs_err;

    /* construct the inverse of the lumped mass matrix */
    ymir_mass_lump (inv_mass_lump);
    if (hess_dir_ref->nefields) {
      ymir_pressure_vec_lump_mass (inv_mass_lump, NULL);
    }
    ymir_vec_fabs (inv_mass_lump, inv_mass_lump);
    ymir_vec_reciprocal (inv_mass_lump);

    /* compute norm of reference */
    ymir_vec_copy (hess_dir_ref, tmp);
    ymir_vec_multiply_in (inv_mass_lump, tmp);
    norm_ref = sqrt (ymir_vec_innerprod (hess_dir_ref, tmp));

    /* compute error */
    ymir_vec_add (-1.0, hess_dir_ref, hess_dir_fd);

    /* compute norm of error */
    ymir_vec_copy (hess_dir_fd, tmp);
    ymir_vec_multiply_in (inv_mass_lump, tmp);
    abs_err = sqrt (ymir_vec_innerprod (hess_dir_fd, tmp));
    *abs_error += abs_err;
    if (SC_EPS * abs_err < norm_ref) {
      *rel_error += abs_err / norm_ref;
    }
    else {
      *rel_error += 1.0;
    }
    n_error_components++;

    /* destroy */
    ymir_vec_destroy (tmp);
    ymir_vec_destroy (inv_mass_lump);
  }

  /* compute the root mean square (RMS) of the error */
  if (hess_dir_ref->n_meshfree) {
    const int           n_entries = hess_dir_ref->n_meshfree;
    int                 i;
    const double       *rf = hess_dir_ref->meshfree->e[0];
    const double       *fd = hess_dir_fd->meshfree->e[0];
    double              diff_sq;
    double              abs_err_sq = 0.0;
    double              rel_err_sq = 0.0;
    int                 n_entries_nonzero = 0;

    RHEA_GLOBAL_VERBOSE ("========================================\n");
    RHEA_GLOBAL_VERBOSEF ("%s\n", __func__);
    RHEA_GLOBAL_VERBOSE ("----------------------------------------\n");
    for (i = 0; i < n_entries; i++) {
      if (isfinite (rf[i]) && isfinite (fd[i]) &&
          (DBL_MIN < fabs (rf[i]) || DBL_MIN < fabs (fd[i]))) {
        diff_sq = (rf[i] - fd[i]) * (rf[i] - fd[i]);
        abs_err_sq += diff_sq;
        if (SC_EPS*SC_EPS * diff_sq < (rf[i]*rf[i])) {
          rel_err_sq += diff_sq / (rf[i]*rf[i]);
        }
        else {
          rel_err_sq += diff_sq;
        }
        n_entries_nonzero++;
        RHEA_GLOBAL_VERBOSEF (
            "param# %3i: value ref %+.6e, fd %+.6e ; error abs %.3e, rel %.3e\n",
            i, rf[i], fd[i], sqrt (diff_sq), sqrt (diff_sq/(rf[i]*rf[i])));
      }
    }
    RHEA_GLOBAL_VERBOSE ("========================================\n");

    *abs_error += sqrt (abs_err_sq / (double) n_entries_nonzero);
    *rel_error += sqrt (rel_err_sq / (double) n_entries_nonzero);
    n_error_components++;
  }

  /* finalize relative error */
  if (1 < n_error_components) {
    *rel_error /= sqrt ((double) n_error_components);
  }
}

void
rhea_newton_check_hessian (ymir_vec_t *solution,
                           ymir_vec_t *neg_gradient,
                           int n_trials,
                           const int iter,
                           rhea_newton_problem_t *nl_problem)
{
  const int           check_gradient =
                        rhea_newton_problem_get_check_gradient (nl_problem);
  ymir_vec_t         *perturb_sol_vec, *perturb_vec;
  ymir_vec_t         *hess_dir_ref, *hess_dir_fd;
  double              abs_error, rel_error;
  sc_dmatrix_t       *result;
  const int           exp_incr = 2;
  const int           exp_min = 0;
  int                 n;

  /* check input */
  RHEA_ASSERT (NULL != solution);
  RHEA_ASSERT (NULL != neg_gradient);
  RHEA_ASSERT (rhea_newton_problem_compute_neg_gradient_exists (nl_problem));

  /* exit if nothing to do */
  if (n_trials <= 0) {
    return;
  }
  /* if needs setting default #trials */
  if (1 == n_trials) {
    n_trials = RHEA_NEWTON_CHECK_N_TRIALS;
  }

  /* exit if check cannot be executed */
  if (!rhea_newton_problem_apply_hessian_exists (nl_problem)) {
    RHEA_GLOBAL_VERBOSEF_FN_TAG (
        __func__, "stop_reason=\"%s\"",
        "Cannot execute because Hessian-apply function is not provided");
    return;
  }

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "newton_iter=%i", iter);

  /* create work vectors */
  perturb_sol_vec = ymir_vec_template (solution);
  perturb_vec = ymir_vec_template (solution);
  hess_dir_ref = ymir_vec_template (neg_gradient);
  hess_dir_fd = ymir_vec_template (neg_gradient);
  result = sc_dmatrix_new (n_trials, 3);

  /* set direction vector */
  rhea_newton_check_set_perturb_vec (
      perturb_vec, solution, rhea_newton_problem_get_mpicomm (nl_problem));

  /* compute the reference directional derivative */
  rhea_newton_problem_apply_hessian (hess_dir_ref, perturb_vec, nl_problem);

  /* set up the finite difference derivative */
  rhea_newton_problem_set_check_gradient (0, nl_problem);

  /* compare reference with finite difference derivative */
  for (n = 0; n < n_trials; n++) {
    const int           exp = exp_min + exp_incr * n;
    const double        eps = pow (10.0, (double) -exp);

    ymir_vec_copy (solution, perturb_sol_vec);
    ymir_vec_add (eps, perturb_vec, perturb_sol_vec);
    rhea_newton_problem_update_operator (perturb_sol_vec, nl_problem);
    rhea_newton_problem_compute_neg_gradient (hess_dir_fd, perturb_sol_vec,
                                              nl_problem);
    ymir_vec_scale (-1.0, hess_dir_fd);
    ymir_vec_add (1.0, neg_gradient, hess_dir_fd);
    ymir_vec_scale (1.0/eps, hess_dir_fd);

    rhea_newton_check_hessian_compute_error (&abs_error, &rel_error,
                                             hess_dir_ref, hess_dir_fd);

    result->e[n][0] = eps;
    result->e[n][1] = abs_error;
    result->e[n][2] = rel_error;
  }

  /* restore */
  rhea_newton_problem_set_check_gradient (check_gradient, nl_problem);
  rhea_newton_problem_update_operator (solution, nl_problem);

  /* print results */
  RHEA_GLOBAL_INFO ("========================================\n");
  RHEA_GLOBAL_INFOF ("%s, newton_iter=%i\n", __func__, iter);
  RHEA_GLOBAL_INFO ("----------------------------------------\n");
  for (n = 0; n < n_trials; n++) {
    RHEA_GLOBAL_INFOF (
        "F.D. eps=%.1e ; error abs=%.1e, rel=%.1e\n",
        result->e[n][0], result->e[n][1], result->e[n][2]);
  }
  RHEA_GLOBAL_INFO ("========================================\n");

  /* destroy */
  ymir_vec_destroy (perturb_sol_vec);
  ymir_vec_destroy (perturb_vec);
  ymir_vec_destroy (hess_dir_ref);
  ymir_vec_destroy (hess_dir_fd);
  sc_dmatrix_destroy (result);

  RHEA_GLOBAL_INFOF_FN_END (__func__, "newton_iter=%i", iter);
}
