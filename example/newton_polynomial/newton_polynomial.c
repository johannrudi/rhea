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
 *   d := [a , b],  d \in \R^2
 *
 * WANT: Find a minimizer `x` for the objective functional:
 *
 *   J(x) := 1/2 * || F(x) - d ||^2
 *
 * where
 *
 *   F(x) := [x , f(x)],  F(x) \in \R^2
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

#include <rhea.h>
#include <ymir_vtk.h>

/**
 * Computes interpolating quadratic polynomial via Hermite interpolation.
 */
static double *
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

/**
 * Evaluates quadratic polynomial.
 */
static double
newton_polynomial_eval_p2 (const double x, const double coeff[3])
{
  return coeff[0] + coeff[1]*x + coeff[2]*x*x;
}

/**
 * Evaluates 1st derivative of quadratic polynomial.
 */
static double
newton_polynomial_eval_p2d (const double x, const double coeff[3])
{
  return coeff[1] + 2.0*coeff[2]*x;
}

/**
 * Evaluates 2nd derivative of quadratic polynomial.
 */
static double
newton_polynomial_eval_p2dd (const double x, const double coeff[3])
{
  return 2.0*coeff[2];
}

/**
 * Evaluates quadratic polynomial element-wise over a vector as input values.
 */
static void
newton_polynomial_evaluate (ymir_vec_t *f_vec, ymir_vec_t *x_vec,
                            const double coeff[3], const int derivative)
{
  const int           has_shared = (x_vec->coff != NULL ? 1 : 0);

  const sc_dmatrix_t *x_owned_mat = x_vec->dataown;
  const double       *x_owned_data = x_owned_mat->e[0];
  const sc_dmatrix_t *x_shared_mat = (has_shared ? x_vec->coff : NULL);
  const double       *x_shared_data = (has_shared ? x_shared_mat->e[0] : NULL);
  const sc_bint_t     size_owned = x_owned_mat->m * x_owned_mat->n;
  const sc_bint_t     size_shared = (has_shared ?
                                     x_shared_mat->m * x_shared_mat->n : 0);

  sc_dmatrix_t       *f_owned_mat = f_vec->dataown;
  double             *f_owned_data = f_owned_mat->e[0];
  sc_dmatrix_t       *f_shared_mat = (has_shared ? f_vec->coff : NULL);
  double             *f_shared_data = (has_shared ? f_shared_mat->e[0] : NULL);
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

/* Object containing information about the problem */
typedef struct newton_polynomial_problem
{
  double              coeff[3];
  ymir_vec_t         *data_a;
  ymir_vec_t         *data_b;
  ymir_vec_t         *sol_x;
}
newton_polynomial_problem_t;

/**
 * Computes the misfit between the data and a trial solution vector.
 */
static void
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
  newton_polynomial_evaluate (misfit_y, sol_x, coeff, 0 /* !deriv */);

  /* compute misfit */
  ymir_vec_copy (sol_x, misfit_x);
  ymir_vec_add (-1.0, data_a, misfit_x);
  ymir_vec_add (-1.0, data_b, misfit_y);

  /* check output */
  RHEA_ASSERT (sc_dmatrix_is_valid (misfit_x->dataown));
  RHEA_ASSERT (!misfit_x->coff || sc_dmatrix_is_valid (misfit_x->coff));
  RHEA_ASSERT (sc_dmatrix_is_valid (misfit_y->dataown));
  RHEA_ASSERT (!misfit_y->coff || sc_dmatrix_is_valid (misfit_x->coff));
}

/**
 * Evaluates the objective functional.
 * (Callback function for Newton's method.)
 */
static double
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

/**
 * Computes the negative gradient of the objective functional.
 * (Callback function for Newton's method.)
 */
static void
newton_polynomial_compute_negative_gradient (ymir_vec_t *neg_gradient,
                                             ymir_vec_t *solution, void *data)
{
  newton_polynomial_problem_t  *poly_problem = data;
  const double       *coeff = poly_problem->coeff;
  ymir_vec_t         *misfit_x = ymir_vec_template (solution);
  ymir_vec_t         *misfit_y = ymir_vec_template (solution);

  /* evaluate derivative of polynomial */
  newton_polynomial_evaluate (neg_gradient, solution, coeff, 1 /* deriv */);

  /* compute misfit */
  newton_polynomial_compute_misfit (misfit_x, misfit_y, solution, poly_problem);

  /* compute negative gradient */
  ymir_vec_multiply_in (misfit_y, neg_gradient);
  ymir_vec_add (1.0, misfit_x, neg_gradient);
  ymir_vec_scale (-1.0, neg_gradient);

  /* destroy */
  ymir_vec_destroy (misfit_x);
  ymir_vec_destroy (misfit_y);

  /* check output */
  RHEA_ASSERT (sc_dmatrix_is_valid (neg_gradient->dataown));
  RHEA_ASSERT (!neg_gradient->coff || sc_dmatrix_is_valid (neg_gradient->coff));
}

/**
 * Computes the norm of the gradient.
 * (Callback function for Newton's method.)
 */
static double
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

/**
 * Applies Hessian operator.
 * (Callback function for Newton's method.)
 */
static void
newton_polynomial_apply_hessian (ymir_vec_t *out, ymir_vec_t *in, void *data)
{
  newton_polynomial_problem_t  *poly_problem = data;
  const int           inverse = 0;

  newton_polynomial_apply_hessian_internal (out, in, poly_problem, inverse);
}

/**
 * Applies inverse of the Hessian operator.
 * (Callback function for Newton's method.)
 */
static int
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

/**
 * Updates the nonlinear operator, given a trial solution vector (nothing to
 * do in this case).
 * (Callback function for Newton's method.)
 */
static void
newton_polynomial_update_operator (ymir_vec_t *solution, void *data)
{
  return;
}

/**
 * Updates the Hessian operator, given a trial solution vector.
 * (Callback function for Newton's method.)
 */
static void
newton_polynomial_update_hessian (ymir_vec_t *solution, void *data)
{
  newton_polynomial_problem_t  *poly_problem = data;

  ymir_vec_copy (solution, poly_problem->sol_x);
}

#if 0
/* TODO outdated
 * J(x) := 1/2 * || [x ; f(x)] - [a ; b] ||^2
 *
 * A(x) := [x ; f(x)]
 *
 * d/dx A(x) * y = [diag(1) ; diag(f'(x))] * y
 */

typedef struct newton_polynomial_problem
{
  double              coeff[3];
  ymir_vec_t         *x_vec;
  ymir_vec_t         *xy_rhs_vec;
}
newton_polynomial_problem_t;

static void
newton_polynomial_apply_op (ymir_vec_t *out, ymir_vec_t *in,
                            newton_polynomial_problem_t *poly_problem)
{
  const double       *coeff = poly_problem->coeff;
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (in);
  const ymir_locidx_t  n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  sc_dmatrix_t       *in_el_mat, *out_el_mat;
  double             *in_el_data, *out_el_data;
  ymir_locidx_t       elid;
  int                 nodeid;

  /* check input */
//TODO
//RHEA_ASSERT (visc_vec != NULL);

  /* create work variables */
  in_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  out_el_mat = sc_dmatrix_new (n_nodes_per_el, 2);

  in_el_data = NULL;
  out_el_data = out_el_mat->e[0];

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    /* get input */
    ymir_dvec_get_elem (in, in_el_mat, YMIR_STRIDE_NODE, elid, YMIR_READ);
    in_el_data = in_el_mat->e[0];

    /* apply operator at each node */
    for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
      const double        in_val = in_el_data[nodeid];
      double              f_in, eval_p2d;

      f_in = newton_polynomial_eval_p2 (in_val, coeff);

      out_el_data[2*nodeid    ] = in_val;
      out_el_data[2*nodeid + 1] = f_in;
    }

    /* set output */
    ymir_dvec_set_elem (out, out_el_mat, YMIR_STRIDE_NODE, elid, YMIR_SET);
  }

  /* destroy */
  sc_dmatrix_destroy (in_el_mat);
  sc_dmatrix_destroy (out_el_mat);
}

static void
newton_polynomial_apply_jac (ymir_vec_t *out, ymir_vec_t *in,
                             newton_polynomial_problem_t *poly_problem,
                             const int inverse)
{
  const double       *coeff = poly_problem->coeff;
  ymir_vec_t         *x_vec = poly_problem->x_vec;
  ymir_vec_t         *rhs_vec = poly_problem->xy_rhs_vec;
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (in);
  const ymir_locidx_t  n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  sc_dmatrix_t       *in_el_mat, *out_el_mat;
  double             *in_el_data, *out_el_data;
  sc_dmatrix_t       *x_el_mat, *rhs_el_mat;
  double             *x_el_data, *rhs_el_data;
  ymir_locidx_t       elid;
  int                 nodeid;

  /* check input */
//TODO
//RHEA_ASSERT (visc_vec != NULL);

  /* create work variables */
  in_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  out_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  x_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  rhs_el_mat = sc_dmatrix_new (n_nodes_per_el, 2);

  in_el_data = NULL;
  out_el_data = out_el_mat->e[0];
  x_el_data = NULL;
  rhs_el_data = NULL;

  for (elid = 0; elid < n_elements; elid++) { /* loop over all elements */
    /* get input */
    ymir_dvec_get_elem (in, in_el_mat, YMIR_STRIDE_NODE, elid, YMIR_READ);
    in_el_data = in_el_mat->e[0];

    /* get current solution and right-hand side */
    ymir_dvec_get_elem (x_vec, x_el_mat, YMIR_STRIDE_NODE, elid, YMIR_READ);
    x_el_data = x_el_mat->e[0];
    ymir_dvec_get_elem (rhs_vec, rhs_el_mat, YMIR_STRIDE_NODE, elid, YMIR_READ);
    rhs_el_data = rhs_el_mat->e[0];

    /* apply operator/Jacobian at each node */
    for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
      const double        in_val = in_el_data[nodeid];
      const double        x_val = x_el_data[nodeid];
      double              f_in, df_x, jac_val;

      f_in = newton_polynomial_eval_p2 (in_val, coeff);
      df_x = newton_polynomial_eval_p2d (x_val, coeff);

      jac_val = (rhs_el_data[2*nodeid    ] - x_val) +
                (rhs_el_data[2*nodeid + 1] - f_in) * df_x;

      switch (inverse) {
      case 0:
        out_el_data[nodeid] = in_val * jac_val;
        break;
      case 1:
        out_el_data[nodeid] = in_val / jac_val;
        break;
      default:
        RHEA_ABORT_NOT_REACHED ();
    }

    /* set output */
    ymir_dvec_set_elem (out, out_el_mat, YMIR_STRIDE_NODE, elid, YMIR_SET);
  }

  /* destroy */
  sc_dmatrix_destroy (in_el_mat);
  sc_dmatrix_destroy (out_el_mat);
  sc_dmatrix_destroy (x_el_mat);
  sc_dmatrix_destroy (rhs_el_mat);
}

void
newton_polynomial_apply_operator (ymir_vec_t *out, ymir_vec_t *in, void *data)
{
  newton_polynomial_problem_t  *poly_problem = data;

  newton_polynomial_apply_op (out, in, poly_problem);
}

void
newton_polynomial_apply_jacobian (ymir_vec_t *out, ymir_vec_t *in, void *data)
{
  newton_polynomial_problem_t  *poly_problem = data;
  const int           inverse = 0;

  newton_polynomial_apply_jac (out, in, poly_problem, inverse);
}

void
newton_polynomial_update_operator (ymir_vec_t *solution, void *data)
{
  return;
}

void
newton_polynomial_update_jacobian (ymir_vec_t *solution, void *data)
{
  newton_polynomial_problem_t  *poly_problem = data;
  ymir_vec_t         *x_vec = poly_problem->x_vec;

  ymir_vec_copy (solution, x_vec);
}

void
newton_polynomial_compute_residual (ymir_vec_t *residual, ymir_vec_t *solution,
                                    void *data)
{
  newton_polynomial_problem_t  *poly_problem = data;
  ymir_vec_t         *xy_rhs_vec = poly_problem->xy_rhs_vec;

  newton_polynomial_apply_op (residual, solution, poly_problem);
  ymir_vec_add (-1.0, xy_rhs_vec, residual);
}

double
newton_polynomial_compute_residual_norm (ymir_vec_t *residual, void *data,
                                         double *res_norm_comp)
{
  return ymir_vec_norm (residual);
}

int
newton_polynomial_solve_linearized (ymir_vec_t *step,
                                    ymir_vec_t *residual,
                                    const int lin_iter_max,
                                    const double lin_res_norm_rtol,
                                    const int nonzero_initial_guess,
                                    void *data,
                                    int *lin_iter_count)
{
  newton_polynomial_problem_t  *poly_problem = data;
  const int           inverse = 1;

  newton_polynomial_apply_jac (step, residual, poly_problem, inverse);

  /* return iteraton count and "stopping" reason */
  if (lin_iter_count != NULL) {
    *lin_iter_count = 1;
  }
  return 1;
}
#endif

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
  const char         *this_fn_name = "newton_polynomial_setup_mesh";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* create p4est */
  *p4est = rhea_discretization_p4est_new (mpicomm, discr_options,
                                          domain_options);

  /* set up boundary, store in `discr_options` */
  rhea_discretization_options_set_boundary (discr_options, *p4est,
                                            domain_options);

  /* create ymir mesh and pressure element */
  rhea_discretization_ymir_mesh_new_from_p4est (ymir_mesh, NULL, *p4est,
                                                discr_options);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
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
  const char         *this_fn_name = "newton_polynomial_setup_newton";
  ymir_vec_t         *tmp_vec = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
  newton_polynomial_problem_t *poly_problem;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* initialize data for Newton problem */
  poly_problem = RHEA_ALLOC (newton_polynomial_problem_t, 1);
  poly_problem->sol_x = ymir_vec_template (tmp_vec);

  /* set coefficient of polynomial function */
  {
    const double        start_node = 0.0;
    const double        start_val = 1.0;
    const double        start_deriv = 0.2;
    const double        end_node = 1.0;
    const double        end_val = 0.0;
    double             *coeff;

    coeff = newton_polynomial_new_p2_coeff_hermite (
        start_node, start_val, start_deriv, end_node, end_val);

    poly_problem->coeff[0] = coeff[0];
    poly_problem->coeff[1] = coeff[1];
    poly_problem->coeff[2] = coeff[2];

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
        neg_gradient_vec, step_vec,
        newton_polynomial_evaluate_objective,
        newton_polynomial_compute_negative_gradient,
        newton_polynomial_compute_gradient_norm,
        newton_polynomial_apply_hessian,
        newton_polynomial_solve_hessian_system,
        newton_polynomial_update_operator,
        newton_polynomial_update_hessian,
        0 /* no multi-component norms */,
        poly_problem);
  }

  /* destroy */
  ymir_vec_destroy (tmp_vec);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
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
  const char         *this_fn_name = "newton_polynomial_setup_clear_all";
  newton_polynomial_problem_t *poly_problem;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

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

  /* destroy (some) options */
  rhea_discretization_options_clear (discr_options);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Runs the program.
 */
int
main (int argc, char **argv)
{
  const char         *this_fn_name = "newton_polynomial:main";
  /* MPI */
  MPI_Comm            mpicomm = MPI_COMM_WORLD;
  int                 mpisize, mpirank, ompsize;
  int                 mpiret;
  /* options */
  ymir_options_t     *opt;
  rhea_domain_options_t         domain_options;
  rhea_discretization_options_t discr_options;
  rhea_newton_options_t         newton_options;
  /* options local to this function */
  int                 production_run;
  char               *vtk_write_solution_path;
  /* mesh */
  p4est_t            *p4est;
  ymir_mesh_t        *ymir_mesh;
  /* Newton */
  rhea_newton_problem_t *nl_problem;
  ymir_vec_t         *solution;

  /*
   * Initialize Libraries
   */

  /* initialize rhea and sub-packages */
  rhea_initialize (argc, argv, mpicomm);

  /* get parallel environment */
  mpiret = MPI_Comm_size (mpicomm, &mpisize); YMIR_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &mpirank); YMIR_CHECK_MPI (mpiret);

#ifdef RHEA_ENABLE_OPENMP
  ompsize = omp_get_max_threads ();
#else
  ompsize = 1;
#endif

  /*
   * Define & Parse Options
   */

  opt = ymir_options_global_new (argv[0] /* program path */);

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  /* basic options */
  YMIR_OPTIONS_CALLBACK, "help", 'h', 0 /* no callback fn args */,
    ymir_options_print_usage_and_exit_fn, NULL /* no arg usage */,
    "Print usage and exit",
  YMIR_OPTIONS_INIFILE, "options-file", 'f',
    ".ini file with option values",

  /* performance & monitoring options */
  YMIR_OPTIONS_B, "production-run", '\0',
    &(production_run), 0,
    "Execute as a production run (to reduce some overhead and checks)",

  /* vtk output options */
  YMIR_OPTIONS_S, "vtk-write-solution-path", '\0',
    &(vtk_write_solution_path), NULL,
    "File path for vtk files for the solution of the Stokes problem",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add sub-options */
  rhea_add_options_newton (opt);
//ymir_options_add_suboptions_solver_stokes (opt); //TODO del

  /* parse options */
  {
    int                 optret;

    optret = ymir_options_parse (SC_LP_INFO, opt, argc, argv);
    if (optret < 0) { /* if parsing was not successful */
      ymir_options_print_usage (SC_LP_INFO, opt, NULL /* args usage */);
      RHEA_GLOBAL_INFO ("Option parsing failed\n");
      exit (0);
    }
  }

  /*
   * Initialize Main Program
   */

  RHEA_GLOBAL_PRODUCTIONF (
      "Into %s (production %i)\n", this_fn_name, production_run);
  RHEA_GLOBAL_PRODUCTIONF (
      "Parallel environment: MPI size %i, OpenMP size %i\n", mpisize, ompsize);
  ymir_set_up (argc, argv, mpicomm, production_run);

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

  solution = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
  ymir_vec_set_zero (solution);
//TODO move following to rhea_newton?
  newton_polynomial_update_hessian (
      solution, rhea_newton_problem_get_data (nl_problem));
  rhea_newton_solve (solution, nl_problem, &newton_options);

  /* write vtk of solution */
  if (vtk_write_solution_path != NULL) {
    newton_polynomial_problem_t *poly_problem =
      (newton_polynomial_problem_t *) rhea_newton_problem_get_data (nl_problem);

    ymir_vtk_write (ymir_mesh, vtk_write_solution_path,
                    poly_problem->data_a, "data_a",
                    poly_problem->data_b, "data_b",
                    solution, "solution", NULL);
  }

  /* destroy */
  ymir_vec_destroy (solution);

  /*
   * Finalize
   */

  /* destroy Newton problem and mesh */
  newton_polynomial_setup_clear_all (nl_problem, p4est, ymir_mesh,
                                     &discr_options);

  /* destroy options */
  ymir_options_global_destroy ();

  /* print that this function is ending */
  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

  /* finalize rhea */
  rhea_finalize ();

  return 0;
}
