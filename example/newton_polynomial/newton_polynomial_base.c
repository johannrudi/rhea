/*
 */

#include <newton_polynomial_base.h>
#include <rhea_base.h>

double *
newton_polynomial_new_p2_coeff_hermite (const double start_node,
                                        const double start_val,
                                        const double start_deriv,
                                        const double end_node,
                                        const double end_val)
{
  double             *coeff;

  /* allocate coefficients */
  coeff = RHEA_ALLOC (double, 3);

  /* calculate coefficients via Hermite interpolation */
  coeff[2] = ((end_val - start_val) / (end_node - start_node) - start_deriv) /
             (end_node - start_node);
  coeff[1] = start_deriv - 2.0 * coeff[2] * start_node;
  coeff[0] = start_val - (start_deriv * start_node) +
             (coeff[2] * start_node * start_node);

  /* return coefficients */
  return coeff;
}

double
newton_polynomial_eval_p2 (const double x, const double coeff[3])
{
  return coeff[0] + coeff[1]*x + coeff[2]*x*x;
}

double
newton_polynomial_eval_p2d (const double x, const double coeff[3])
{
  return coeff[1] + 2.0*coeff[2]*x;
}

double
newton_polynomial_eval_p2dd (const double x, const double coeff[3])
{
  return 2.0*coeff[2];
}

void
newton_polynomial_evaluate (ymir_vec_t *f_vec, ymir_vec_t *x_vec,
                            const double coeff[3], const int derivative)
{
  const int           has_shared = (x_vec->coff != NULL ? 1 : 0);

  const sc_dmatrix_t *x_owned_mat = x_vec->dataown;
  const double       *x_owned_data = x_owned_mat->e[0];
  const sc_dmatrix_t *x_shared_mat = (has_shared ? x_vec->coff : NULL);
  const double       *x_shared_data = (has_shared ? x_shared_mat->e[0] : NULL);

  sc_dmatrix_t       *f_owned_mat = f_vec->dataown;
  double             *f_owned_data = f_owned_mat->e[0];
  sc_dmatrix_t       *f_shared_mat = (has_shared ? f_vec->coff : NULL);
  double             *f_shared_data = (has_shared ? f_shared_mat->e[0] : NULL);

  const sc_bint_t     size_owned = x_owned_mat->m * x_owned_mat->n;
  const sc_bint_t     size_shared = (has_shared ?
                                     x_shared_mat->m * x_shared_mat->n : 0);
  sc_bint_t           i;

  /* check input */
  RHEA_ASSERT (sc_dmatrix_is_valid (x_vec->dataown));
  RHEA_ASSERT (!has_shared || sc_dmatrix_is_valid (x_vec->coff));
  RHEA_ASSERT (!has_shared || f_vec->coff != NULL);
  RHEA_ASSERT ( (f_owned_mat->m * f_owned_mat->n) == size_owned );
  RHEA_ASSERT ( !has_shared ||
                (f_shared_mat->m * f_shared_mat->n) == size_shared );

  /* evaluate function for processor-owned portion of the vector */
  switch (derivative) {
  case 0:
    for (i = 0; i < size_owned; i++) {
      f_owned_data[i] = newton_polynomial_eval_p2 (x_owned_data[i], coeff);
    }
    break;
  case 1:
    for (i = 0; i < size_owned; i++) {
      f_owned_data[i] = newton_polynomial_eval_p2d (x_owned_data[i], coeff);
    }
    break;
  case 2:
    for (i = 0; i < size_owned; i++) {
      f_owned_data[i] = newton_polynomial_eval_p2dd (x_owned_data[i], coeff);
    }
    break;
  default:
    RHEA_ABORT_NOT_REACHED ();
  }

  /* exit if no shared portion */
  if (!has_shared) {
    return;
  }

  /* evaluate function for processor-shared portion of the vector */
  switch (derivative) {
  case 0:
    for (i = 0; i < size_shared; i++) {
      f_shared_data[i] = newton_polynomial_eval_p2 (x_shared_data[i], coeff);
    }
    break;
  case 1:
    for (i = 0; i < size_shared; i++) {
      f_shared_data[i] = newton_polynomial_eval_p2d (x_shared_data[i], coeff);
    }
    break;
  case 2:
    for (i = 0; i < size_shared; i++) {
      f_shared_data[i] = newton_polynomial_eval_p2dd (x_shared_data[i], coeff);
    }
    break;
  default:
    RHEA_ABORT_NOT_REACHED ();
  }

  /* check output */
  RHEA_ASSERT (sc_dmatrix_is_valid (f_vec->dataown));
  RHEA_ASSERT (!has_shared || sc_dmatrix_is_valid (f_vec->coff));
}

void
newton_polynomial_compute_misfit (ymir_vec_t *misfit_x,
                                  ymir_vec_t *misfit_y,
                                  ymir_vec_t *sol_x,
                                  newton_polynomial_problem_t *poly_problem)
{
  const double       *coeff = poly_problem->coeff;
  ymir_vec_t         *data_a = poly_problem->data_a;
  ymir_vec_t         *data_b = poly_problem->data_b;

  /* check output */
  RHEA_ASSERT (sc_dmatrix_is_valid (data_a->dataown));
  RHEA_ASSERT (!data_a->coff || sc_dmatrix_is_valid (data_a->coff));
  RHEA_ASSERT (sc_dmatrix_is_valid (data_b->dataown));
  RHEA_ASSERT (!data_b->coff || sc_dmatrix_is_valid (data_b->coff));

  /* evaluate polynomial */
  ymir_vec_copy (sol_x, misfit_x);
  newton_polynomial_evaluate (misfit_y, misfit_x, coeff, 0 /* !deriv */);

  /* subtract data */
  ymir_vec_add (-1.0, data_a, misfit_x);
  ymir_vec_add (-1.0, data_b, misfit_y);

  /* check output */
  RHEA_ASSERT (sc_dmatrix_is_valid (misfit_x->dataown));
  RHEA_ASSERT (!misfit_x->coff || sc_dmatrix_is_valid (misfit_x->coff));
  RHEA_ASSERT (sc_dmatrix_is_valid (misfit_y->dataown));
  RHEA_ASSERT (!misfit_y->coff || sc_dmatrix_is_valid (misfit_x->coff));
}

void
newton_polynomial_compute_distance (ymir_vec_t *dist, ymir_vec_t *x_vec,
                                    newton_polynomial_problem_t *poly_problem)
{
  const int           has_shared = (x_vec->coff != NULL ? 1 : 0);
  ymir_vec_t         *misfit_x = ymir_vec_template (x_vec);
  ymir_vec_t         *misfit_y = ymir_vec_template (x_vec);

  const sc_dmatrix_t *misfit_x_owned_mat = misfit_x->dataown;
  const double       *misfit_x_owned_data = misfit_x_owned_mat->e[0];
  const sc_dmatrix_t *misfit_x_shared_mat =
                        (has_shared ? misfit_x->coff : NULL);
  const double       *misfit_x_shared_data =
                        (has_shared ? misfit_x_shared_mat->e[0] : NULL);

  const sc_dmatrix_t *misfit_y_owned_mat = misfit_y->dataown;
  const double       *misfit_y_owned_data = misfit_y_owned_mat->e[0];
  const sc_dmatrix_t *misfit_y_shared_mat =
                        (has_shared ? misfit_y->coff : NULL);
  const double       *misfit_y_shared_data =
                        (has_shared ? misfit_y_shared_mat->e[0] : NULL);

  sc_dmatrix_t       *dist_owned_mat = dist->dataown;
  double             *dist_owned_data = dist_owned_mat->e[0];
  sc_dmatrix_t       *dist_shared_mat = (has_shared ? dist->coff : NULL);
  double             *dist_shared_data =
                        (has_shared ? dist_shared_mat->e[0] : NULL);

  const sc_bint_t     size_owned = dist_owned_mat->m * dist_owned_mat->n;
  const sc_bint_t     size_shared = (has_shared ?
                                     dist_shared_mat->m * dist_shared_mat->n :
                                     0);
  sc_bint_t           i;

  /* compute misfit */
  newton_polynomial_compute_misfit (misfit_x, misfit_y, x_vec, poly_problem);

  /* compute distance for owned nodes */
  for (i = 0; i < size_owned; i++) {
    const double        dx = misfit_x_owned_data[i];
    const double        dy = misfit_y_owned_data[i];

    dist_owned_data[i] = sqrt (dx*dx + dy*dy);
  }

  /* compute distance for shared nodes */
  if (has_shared) {
    for (i = 0; i < size_shared; i++) {
      const double        dx = misfit_x_shared_data[i];
      const double        dy = misfit_y_shared_data[i];

      dist_shared_data[i] = sqrt (dx*dx + dy*dy);
    }
  }

  /* destroy */
  ymir_vec_destroy (misfit_x);
  ymir_vec_destroy (misfit_y);
}

/******************************************************************************
 * Callback Functions for Newton Solver
 *****************************************************************************/

void
newton_polynomial_data_init (ymir_vec_t *solution, void *data)
{
  newton_polynomial_problem_t  *poly_problem = data;

  /* setup Hessian (see also: newton_polynomial_update_hessian (...)) */
  if (solution != NULL) { /* if nonzero initial guess */
    ymir_vec_copy (solution, poly_problem->sol_x);
  }
  else { /* otherwise assume zero initial guess */
    ymir_vec_set_value (poly_problem->sol_x, 0.0);
  }
}

double
newton_polynomial_evaluate_objective (ymir_vec_t *solution, void *data)
{
  newton_polynomial_problem_t  *poly_problem = data;
  ymir_vec_t         *misfit_x = ymir_vec_template (solution);
  ymir_vec_t         *misfit_y = ymir_vec_template (solution);
  double              obj_val;

  /* compute misfit */
  newton_polynomial_compute_misfit (misfit_x, misfit_y, solution, poly_problem);

  /* evaluate objective functional */
  obj_val = ymir_vec_innerprod (misfit_x, misfit_x) +
            ymir_vec_innerprod (misfit_y, misfit_y);
  obj_val *= 0.5;

  /* destroy */
  ymir_vec_destroy (misfit_x);
  ymir_vec_destroy (misfit_y);

  /* return value of objective functional */
  return obj_val;
}

void
newton_polynomial_compute_negative_gradient (ymir_vec_t *neg_gradient,
                                             ymir_vec_t *solution, void *data)
{
  newton_polynomial_problem_t  *poly_problem = data;
  const double       *coeff = poly_problem->coeff;
  ymir_vec_t         *misfit_x = ymir_vec_template (poly_problem->data_a);
  ymir_vec_t         *misfit_y = ymir_vec_template (poly_problem->data_b);
  ymir_vec_t         *sol;

  /* check for NULL input (mandatory) */
  if (solution == NULL) {
    sol = ymir_vec_template (poly_problem->sol_x);
    ymir_vec_set_zero (sol);
  }
  else {
    sol = solution;
  }

  /* compute misfit */
  newton_polynomial_compute_misfit (misfit_x, misfit_y, sol, poly_problem);

  /* evaluate derivative of polynomial */
  newton_polynomial_evaluate (neg_gradient, sol, coeff, 1 /* deriv */);

  /* compute negative gradient */
  ymir_vec_multiply_in (misfit_y, neg_gradient);
  ymir_vec_add (1.0, misfit_x, neg_gradient);
  ymir_vec_scale (-1.0, neg_gradient);

  /* destroy */
  ymir_vec_destroy (misfit_x);
  ymir_vec_destroy (misfit_y);
  if (solution == NULL) {
    ymir_vec_destroy (sol);
  }

  /* check output */
  RHEA_ASSERT (sc_dmatrix_is_valid (neg_gradient->dataown));
  RHEA_ASSERT (!neg_gradient->coff || sc_dmatrix_is_valid (neg_gradient->coff));
}

double
newton_polynomial_compute_gradient_norm (ymir_vec_t *neg_gradient, void *data,
                                         double *norm_comp)
{
  RHEA_ASSERT (norm_comp == NULL);
  return ymir_vec_norm (neg_gradient);
}

/**
 * Applies the Hessian operator or its inverse.
 */
static void
newton_polynomial_apply_hessian_internal (
                                    ymir_vec_t *out, ymir_vec_t *in,
                                    newton_polynomial_problem_t *poly_problem,
                                    const int inverse)
{
  const double       *coeff = poly_problem->coeff;
  ymir_vec_t         *sol_x = poly_problem->sol_x;
  ymir_vec_t         *misfit_x = ymir_vec_template (sol_x);
  ymir_vec_t         *misfit_y = ymir_vec_template (sol_x);
  ymir_vec_t         *df = misfit_x;

  /* evaluate 2nd derivative of polynomial */
  newton_polynomial_evaluate (out, sol_x, coeff, 2 /* 2nd deriv */);

  /* multiply element-wise by misfit */
  newton_polynomial_compute_misfit (misfit_x, misfit_y, sol_x, poly_problem);
  ymir_vec_multiply_in (misfit_y, out);

  /* add square of 1st derivative of polynomial */
  newton_polynomial_evaluate (df, sol_x, coeff, 1 /* 1st deriv */);
  ymir_vec_multiply_in (df, df);
  ymir_vec_add (1.0, df, out);

  /* add one */
  ymir_vec_shift (1.0, out);

  /* invert Hessian element-wise (stored in `out`) */
  if (inverse) {
    ymir_vec_reciprocal (out);
  }

  /* multiply by input vector element-wise */
  ymir_vec_multiply_in (in, out);

  /* destroy */
  ymir_vec_destroy (misfit_x);
  ymir_vec_destroy (misfit_y);
}

void
newton_polynomial_apply_hessian (ymir_vec_t *out, ymir_vec_t *in, void *data)
{
  newton_polynomial_problem_t  *poly_problem = data;
  const int           inverse = 0;

  newton_polynomial_apply_hessian_internal (out, in, poly_problem, inverse);
}

int
newton_polynomial_solve_hessian_system (ymir_vec_t *step,
                                        ymir_vec_t *neg_gradient,
                                        const int lin_iter_max,
                                        const double lin_res_norm_rtol,
                                        const int nonzero_initial_guess,
                                        void *data,
                                        int *lin_iter_count)
{
  newton_polynomial_problem_t  *poly_problem = data;
  const int           inverse = 1;

  newton_polynomial_apply_hessian_internal (step, neg_gradient, poly_problem,
                                            inverse);

  /* return iteraton count and "stopping" reason */
  if (lin_iter_count != NULL) {
    *lin_iter_count = 1;
  }
  return 1;
}

void
newton_polynomial_update_hessian (ymir_vec_t *solution, ymir_vec_t *step_vec,
                                  const double step_length, void *data)
{
  newton_polynomial_problem_t  *poly_problem = data;

  ymir_vec_copy (solution, poly_problem->sol_x);
}
