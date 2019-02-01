#include <adjoint_adjoint.h>

/*********************************************************************************************
 * Sets up a linear Stokes problem.
 *********************************************************************************************/
void
adjoint_stokes_new (rhea_stokes_problem_t **stokes_problem,
                    ymir_mesh_t **ymir_mesh,
                    ymir_pressure_elem_t **press_elem,
                    rhea_domain_options_t *domain_options,
                    rhea_temperature_options_t *temp_options,
                    rhea_weakzone_options_t *weak_options,
                    rhea_viscosity_options_t *visc_options,
                    subd_options_t *subd_options)
{
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (*ymir_mesh);
  ymir_vec_t         *temperature;

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);

  /* set up data */
  rhea_weakzone_data_create (weak_options, mpicomm);

  /* compute temperature */
  temperature = rhea_temperature_new (*ymir_mesh);
  subd_compute_temperature (temperature, temp_options, subd_options);

  subd_set_velocity_dirichlet_bc (domain_options, subd_options);

  /* create Stokes problem */
  *stokes_problem = rhea_stokes_problem_new (
      *ymir_mesh, *press_elem, temperature, domain_options, temp_options,
      weak_options, visc_options);

//  subd_compute_weakzone (stokes_problem, subd_options); //see whether it is necessary for different cases

  /* set visc function for user defined viscosity   */
  if (subd_options->visc_options->type != SUBD_VISC_RHEA) {
    subd_viscosity_set_function (*stokes_problem, subd_options);
  }


  RHEA_GLOBAL_INFO_FN_END (__func__);
}

/********************************************************************************************
 * create adjoint_problem
 ********************************************************************************************/
adjoint_problem_t *
adjoint_problem_new (ymir_vec_t *solution,
                     rhea_stokes_problem_t *stokes_problem,
                     p4est_t *p4est, ymir_mesh_t *ymir_mesh,
                     ymir_pressure_elem_t *press_elem,
                     rhea_discretization_options_t *discr_options,
                     rhea_temperature_options_t *temp_options,
                     subd_options_t *subd_options,
                     int       solver_iter_max,
                     double    solver_rel_tol)
{
  adjoint_problem_t *adjoint_problem = RHEA_ALLOC (adjoint_problem_t, 1);
  int         i,  n = subd_options->adjoint_options->n_components;

  adjoint_problem->stokes_problem = stokes_problem;
  adjoint_problem->p4est = p4est;
  adjoint_problem->ymir_mesh = ymir_mesh;
  adjoint_problem->press_elem = press_elem;
  adjoint_problem->discr_options = discr_options;
  adjoint_problem->temp_options = temp_options;
  adjoint_problem->subd_options = subd_options;

  adjoint_problem->solver_iter_max = solver_iter_max;
  adjoint_problem->solver_rel_tol = solver_rel_tol;

  adjoint_problem->msol = ymir_vec_template (solution);
  adjoint_problem->neg_grad = ymir_vec_template (solution);
  adjoint_problem->hessian = ymir_vec_new_meshfree ( n *  n);

  adjoint_problem->sol_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);
  adjoint_problem->usol = rhea_velocity_new (ymir_mesh);
  adjoint_problem->vsol = rhea_velocity_new (ymir_mesh);
  adjoint_problem->Husol = YMIR_ALLOC (ymir_vec_t *,  n);
  adjoint_problem->Hvsol = YMIR_ALLOC (ymir_vec_t *,  n);
  adjoint_problem->coeff = YMIR_ALLOC (ymir_vec_t *,  n);
  for (i = 0; i <  n; i++)  {
    adjoint_problem->Husol[i] = rhea_velocity_new (ymir_mesh);
    adjoint_problem->Hvsol[i] = rhea_velocity_new (ymir_mesh);
    adjoint_problem->coeff[i] = rhea_viscosity_new (ymir_mesh);
  }
  subd_adjoint_options_t *adjoint_options = subd_options->adjoint_options;
  char    *record_path = adjoint_options->txt_write_record_path;
  if (record_path) {
    adjoint_problem->fp_record = fopen (record_path, "w");
  }

  return (adjoint_problem);
}

void
adjoint_problem_destroy (void *data)
{
  adjoint_problem_t *adjoint_problem = data;
  FILE  *fp_record = adjoint_problem->fp_record;
  int    i,  n = adjoint_problem->subd_options->adjoint_options->n_components;

  RHEA_ASSERT (adjoint_problem->msol != NULL);
  RHEA_ASSERT (adjoint_problem->usol != NULL);
  RHEA_ASSERT (adjoint_problem->vsol != NULL);
  RHEA_ASSERT (adjoint_problem->Husol != NULL);
  RHEA_ASSERT (adjoint_problem->Hvsol != NULL);
  RHEA_ASSERT (adjoint_problem->sol_vel_press != NULL);

  if (fp_record) { //TODO
    fclose (fp_record);
  }

  {
    const char *solution_path =
        adjoint_problem->subd_options->adjoint_options->vtk_write_solution_path;
    if (solution_path) {
      rhea_stokes_problem_t *stokes_problem = adjoint_problem->stokes_problem;
      ymir_mesh_t *ymir_mesh = adjoint_problem->ymir_mesh;
      ymir_pressure_elem_t *press_elem = adjoint_problem->press_elem;
      ymir_vec_t *sol_vel_press = adjoint_problem->sol_vel_press;

        char   path[BUFSIZ];
        snprintf (path, BUFSIZ, "%s_last", solution_path);
        subd_vtk_write_solution (ymir_mesh, press_elem, stokes_problem, sol_vel_press);
    }
  }

  ymir_vec_destroy (adjoint_problem->msol);
  ymir_vec_destroy (adjoint_problem->hessian);
  ymir_vec_destroy (adjoint_problem->neg_grad);
  rhea_velocity_destroy (adjoint_problem->usol);
  rhea_velocity_destroy (adjoint_problem->vsol);
  rhea_velocity_pressure_destroy (adjoint_problem->sol_vel_press);

  for (i = 0; i <  n; i++)  {
    rhea_velocity_destroy (adjoint_problem->Husol[i]);
    rhea_velocity_destroy (adjoint_problem->Hvsol[i]);
    rhea_viscosity_destroy (adjoint_problem->coeff[i]);
  }
  YMIR_FREE (adjoint_problem->Husol);
  YMIR_FREE (adjoint_problem->Hvsol);
  YMIR_FREE (adjoint_problem->coeff);

  RHEA_FREE (adjoint_problem);

}

/********************************************************************************************
 * create newton_problem
 ********************************************************************************************/
void
adjoint_setup_newton (rhea_newton_problem_t **newton_problem,
                      adjoint_problem_t *adjoint_problem)
{
  ymir_vec_t *sol_vec = adjoint_problem->msol;

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);

  /*create newton problem */
  {
    ymir_vec_t    *neg_grad_vec = ymir_vec_template (sol_vec);
    ymir_vec_t    *step_vec = ymir_vec_template (sol_vec);

    *newton_problem = rhea_newton_problem_new (
        adjoint_solve_negative_gradient,
        adjoint_compute_gradient_norm,
        0 /* no multi-component norms */,
        adjoint_solve_hessian_system_fn);

    rhea_newton_problem_set_vectors (
        neg_grad_vec, step_vec, *newton_problem);

    rhea_newton_problem_set_data_fn (
        adjoint_problem, adjoint_data_init, adjoint_problem_destroy,
        *newton_problem);

    rhea_newton_problem_set_evaluate_objective_fn (
        adjoint_evaluate_objective,
        *newton_problem);

    rhea_newton_problem_set_apply_hessian_fn (
       /* adjoint_apply_hessian_fn */ NULL, *newton_problem);

    rhea_newton_problem_set_update_fn (
        adjoint_update_operator_fn,
        adjoint_update_hessian_fn,
        NULL, /* adjoint_modify_hessian_system_fn, is not necessary */
        *newton_problem);

    rhea_newton_problem_set_output_fn (
        adjoint_output_prestep,
        *newton_problem);
  }

#ifdef RHEA_ENABLE_DEBUG
  rhea_newton_problem_set_checks (1 /* grad */, 1 /* Hessian */, *newton_problem);
#endif

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

void
adjoint_destroy_newton (rhea_newton_problem_t *newton_problem)
{
  /* destroy Newton problem */
  ymir_vec_destroy (rhea_newton_problem_get_neg_gradient_vec (newton_problem));
  ymir_vec_destroy (rhea_newton_problem_get_step_vec (newton_problem));
  rhea_newton_problem_destroy (newton_problem);

}
/********************************************************************************************
 * user defined function used in newton solver
 ********************************************************************************************/
/****************helper functions for Stencil***************************/
int
adjoint_get_field (int *field_nums, int n_components, int i)
{
  int digit, field_indicator;
  subd_adjoint_stencil_field_t field;

  digit = n_components - i;
  field_indicator = *field_nums / int_pow (100, digit - 1);

RHEA_GLOBAL_PRODUCTIONF ("TESTTESTTEST: i=%d digit=%d, stencil_indicator=%d, stencil_nums=%d, critical=%d\n", i, digit, field_indicator, *field_nums, int_pow (100, digit - 1));
  field = (subd_adjoint_stencil_field_t) field_indicator;
  *field_nums -= field_indicator * int_pow (100, digit - 1);

  return (field);
}

static void
subd_adjoint_stencil_visc_rhea_um_elem (double *_sc_restrict visc_elem,
                                            const double *_sc_restrict x,
                                            const double *_sc_restrict y,
                                            const double *_sc_restrict z,
                                            const int n_nodes,
                                            const int *Vmask,
                                            subd_options_t *subd_options)
{
  int                 nodeid;
  double              value = subd_options->adjoint_options->stencil_options->value;
  int                 m;
  rhea_domain_options_t  *domain_options = subd_options->domain_options;
  int                 is_in_upper_mantle;

  /* set flag if location is in lower or upper mantle */
  is_in_upper_mantle = rhea_domain_elem_is_in_upper_mantle (x, y, z, Vmask,
                                    domain_options);

  /* compute viscosity in this element */
  if (is_in_upper_mantle)
    for (nodeid = 0; nodeid < n_nodes; nodeid++) {
        visc_elem[nodeid] = value;
    }
  else
    for (nodeid = 0; nodeid < n_nodes; nodeid++) {
        visc_elem[nodeid] = .0;
    }

    /* check viscosity for `nan`, `inf`, and positivity */
    RHEA_ASSERT (isfinite (visc_elem[nodeid]));
}

static void
subd_adjoint_stencil_visc_rhea_lm_elem (double *_sc_restrict visc_elem,
                                            const double *_sc_restrict x,
                                            const double *_sc_restrict y,
                                            const double *_sc_restrict z,
                                            const int n_nodes,
                                            const int *Vmask,
                                            subd_options_t *subd_options)
{
  int                 nodeid;
  double              value = subd_options->adjoint_options->stencil_options->value;
  int                 m;
  rhea_domain_options_t  *domain_options = subd_options->domain_options;
  int                 is_in_upper_mantle;

  /* set flag if location is in lower or upper mantle */
  is_in_upper_mantle = rhea_domain_elem_is_in_upper_mantle (x, y, z, Vmask,
                                    domain_options);

  /* compute viscosity in this element */
  if (is_in_upper_mantle)
    for (nodeid = 0; nodeid < n_nodes; nodeid++) {
        visc_elem[nodeid] = .0;
    }
  else
    for (nodeid = 0; nodeid < n_nodes; nodeid++) {
        visc_elem[nodeid] = value;
    }

    /* check viscosity for `nan`, `inf`, and positivity */
    RHEA_ASSERT (isfinite (visc_elem[nodeid]));
}

static void
subd_adjoint_stencil_visc_custom_layers_um_elem (double *_sc_restrict visc_elem,
                                            const double *_sc_restrict x,
                                            const double *_sc_restrict y,
                                            const double *_sc_restrict z,
                                            const int n_nodes,
                                            subd_options_t *subd_options)
{
  int                 nodeid;
  double              z_mid = subd_options->visc_options->z_lith;
  double              value = subd_options->adjoint_options->stencil_options->value;
  int                 m;

  /* compute viscosity in this element */
    for (nodeid = 0; nodeid < n_nodes; nodeid++) {
      if (z[nodeid] > z_mid)  {
        visc_elem[nodeid] = value;
      }
      else {
        visc_elem[nodeid] = .0;
      }
    }
    /* check viscosity for `nan`, `inf`, and positivity */
    RHEA_ASSERT (isfinite (visc_elem[nodeid]));
}

void
adjoint_setup_stencil  (ymir_vec_t *stencil,
                        subd_options_t *subd_options,
                        subd_adjoint_stencil_field_t stencil_field)
{
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (stencil);
  const ymir_locidx_t  n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);
  const int           *Vmask = ymir_mesh_get_vertex_indices (mesh);

  mangll_t           *mangll = mesh->ma;
  const int           N = ymir_n (mangll->N);

  sc_dmatrix_t       *stencil_el_mat;
  double             *x, *y, *z, *tmp_el,*stencil_el_data;
  ymir_locidx_t       elid;

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);

  /* create work variables */
  stencil_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  stencil_el_data = stencil_el_mat->e[0];
  x = RHEA_ALLOC (double, n_nodes_per_el);
  y = RHEA_ALLOC (double, n_nodes_per_el);
  z = RHEA_ALLOC (double, n_nodes_per_el);
  tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

    /* compute stencilosity stencil*/
  for (elid = 0; elid < n_elements; elid++) {
    /* get coordinates of this element at Gauss nodes */
    ymir_mesh_get_elem_coord_gauss (x, y, z, elid, mesh, tmp_el);

    switch (stencil_field) {
      case SUBD_ADJOINT_STENCIL_VISC_RHEA_UM:
        subd_adjoint_stencil_visc_rhea_um_elem (stencil_el_data, x, y, z, n_nodes_per_el,
                                          Vmask, subd_options);
      break;

      case SUBD_ADJOINT_STENCIL_VISC_RHEA_LM:
        subd_adjoint_stencil_visc_rhea_lm_elem (stencil_el_data, x, y, z, n_nodes_per_el,
                                          Vmask, subd_options);
      break;

      case SUBD_ADJOINT_STENCIL_VISC_CUSTOM_LAYERS_UM:
        subd_adjoint_stencil_visc_custom_layers_um_elem (stencil_el_data, x, y, z, n_nodes_per_el,
                                         subd_options);
      break;

      default:
      RHEA_ABORT_NOT_REACHED ();
    }

    /* set viscosity of this element */
    rhea_viscosity_set_elem_gauss (stencil, stencil_el_mat, elid);
  }

  /* destroy */
  sc_dmatrix_destroy (stencil_el_mat);
  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

/**************** viscosity gradient **************************************/

subd_adjoint_visc_parameter_t
adjoint_get_parameter (int *parameter_nums, int n_components, int i)
{
  int digit, para_indicator;
  subd_adjoint_visc_parameter_t parameter;

  digit = n_components - i;
  para_indicator = *parameter_nums / int_pow (10, digit - 1);
RHEA_GLOBAL_PRODUCTIONF ("TESTTESTTEST: i=%d digit=%d, para_indicator=%d, parameter_nums=%d, critical=%d\n", i, digit, para_indicator, *parameter_nums, int_pow (10, digit - 1));

  parameter = (subd_adjoint_visc_parameter_t) para_indicator;
  *parameter_nums -= para_indicator * int_pow (10, digit - 1);
  return (parameter);
}

static void
adjoint_compute_visc_grad_m_prefactor (ymir_vec_t   *visc_grad)
{
  ymir_vec_set_value (visc_grad, 1.0);
  return;
}

static void
adjoint_compute_visc_grad_m_activation_energy (ymir_vec_t *grad,
                      ymir_vec_t *visc, ymir_vec_t *temp)
{
  ymir_mesh_t *ymir_mesh = ymir_vec_get_mesh (grad);
  const int           N  = ymir_n (ymir_mesh->cnodes->N);
  const int           Np = ymir_np (ymir_mesh->cnodes->N);
  const int           K  = ymir_mesh->cnodes->K;
  sc_dmatrix_t       *elemtemp = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemvisc = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemgrad = sc_dmatrix_new (1, Np);

  int elid, gp;
  double neutral_temp = RHEA_TEMPERATURE_NEUTRAL_VALUE;

  /*visc * (RHEA_TEMPERATURE_NEUTRAL_VALUE - temp)*/
  for (elid = 0; elid < K; elid++) {
    ymir_cvec_get_elem_interp (temp, elemtemp, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_COPY);
    ymir_dvec_get_elem (visc, elemvisc, YMIR_STRIDE_NODE, elid,
                               YMIR_COPY);
    ymir_dvec_get_elem (grad, elemgrad, YMIR_STRIDE_NODE, elid,
                               YMIR_WRITE);
    double     *_sc_restrict tempd = elemtemp->e[0];
    double     *_sc_restrict viscd = elemvisc->e[0];
    double     *_sc_restrict gradd = elemgrad->e[0];
    for (gp = 0; gp < Np; gp++) {
      gradd[gp] = viscd[gp] * (neutral_temp - tempd[gp]);
    }
    ymir_dvec_set_elem (grad, elemgrad, YMIR_STRIDE_NODE, elid,
                        YMIR_SET);
  }

  sc_dmatrix_destroy (elemtemp);
  sc_dmatrix_destroy (elemvisc);
  sc_dmatrix_destroy (elemgrad);
}

void
adjoint_compute_visc_grad_m (ymir_vec_t  *visc_grad,
                         rhea_stokes_problem_t *stokes_problem,
                         subd_options_t     *subd_options,
                         subd_adjoint_visc_parameter_t parameter)
{
  ymir_mesh_t     *ymir_mesh = ymir_vec_get_mesh (visc_grad);
  ymir_vec_t      *visc = rhea_viscosity_new (ymir_mesh);
  ymir_vec_t      *temp;

  rhea_stokes_problem_copy_viscosity (visc, stokes_problem);
  temp = rhea_stokes_problem_get_temperature (stokes_problem);

  switch (parameter) {
  case (SUBD_ADJOINT_VISC_PRE):
    adjoint_compute_visc_grad_m_prefactor (visc_grad);
    break;
  case (SUBD_ADJOINT_VISC_E):
    adjoint_compute_visc_grad_m_activation_energy (visc_grad, visc, temp);
    break;
  default:
    RHEA_ABORT_NOT_REACHED ();
  }

  if (subd_options->adjoint_options->use_exponent) {
    ymir_vec_multiply_in (visc, visc_grad);
  }

  rhea_viscosity_destroy (visc);
}

/*******************Stokes solves*********************************************/
void
adjoint_run_solver (ymir_vec_t *usol,
                    adjoint_problem_t *adjoint_problem)
{
  ymir_vec_t    *sol_vel_press = adjoint_problem->sol_vel_press;
  rhea_stokes_problem_t *stokes_problem = adjoint_problem->stokes_problem;
  ymir_mesh_t   *ymir_mesh = adjoint_problem->ymir_mesh;
  ymir_pressure_elem_t *press_elem = adjoint_problem->press_elem;
  int iter_max = adjoint_problem->solver_iter_max;
  double rel_tol = adjoint_problem->solver_rel_tol;
  const int nonzero_initial_guess = 0;

  ymir_vec_set_zero (sol_vel_press);
  subd_run_solver (sol_vel_press, ymir_mesh, press_elem, stokes_problem,
                   nonzero_initial_guess, iter_max, rel_tol);

  ymir_upvec_get_u (sol_vel_press, usol, YMIR_COPY);
}

/**********************Initialize parameter and solve****************/

static void
adjoint_parameter_init (ymir_vec_t *solution, ymir_vec_t *msol,
                            subd_options_t *subd_options)
{
  subd_adjoint_options_t *adjoint_options = subd_options->adjoint_options;
  subd_adjoint_visc_field_t visc_field;
  int            n = adjoint_options->n_components;
  int           field_nums = adjoint_options->visc_fields;
  int           use_exponent = adjoint_options->use_exponent;
  double        m_para;
  int           i;

  if (solution == NULL) {
    RHEA_GLOBAL_PRODUCTIONF ("set newton_options.nonzero_initial_guess as %d\n", 1);
    return;
  }

  for (i = 0; i <  n; i++)  {
    visc_field = (subd_adjoint_visc_field_t ) adjoint_get_field (&field_nums,  n, i);
    switch (visc_field) {
    case SUBD_ADJOINT_FIELD_VISC_RHEA_UM:
      m_para = subd_options->rhea_visc_options->upper_mantle_scaling;
    break;

    case SUBD_ADJOINT_FIELD_VISC_RHEA_LM:
      m_para = subd_options->rhea_visc_options->lower_mantle_scaling;
    break;

    case SUBD_ADJOINT_FIELD_VISC_RHEA_UM_E:
      m_para = subd_options->rhea_visc_options->upper_mantle_arrhenius_activation_energy;
      break;

    case SUBD_ADJOINT_FIELD_VISC_RHEA_LM_E:
      m_para = subd_options->rhea_visc_options->lower_mantle_arrhenius_activation_energy;
      break;

    case SUBD_ADJOINT_FIELD_VISC_CUSTOM_LAYERS_UM:
      m_para = subd_options->visc_options->visc_lith;

    default:
      RHEA_ABORT_NOT_REACHED ();
    }
    if (!use_exponent)
      solution->meshfree->e[0][i] = m_para;
    else
      solution->meshfree->e[0][i] = log(m_para);
  }

}

void
adjoint_stokes_update_init (rhea_stokes_problem_t *stokes_problem,
                          p4est_t     *p4est,
                          ymir_mesh_t *ymir_mesh,
                          ymir_pressure_elem_t *press_elem,
                          rhea_discretization_options_t *discr_options,
                          rhea_temperature_options_t *temp_options,
                          subd_options_t *subd_options)

{
  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);

  subd_options->rhs_type = SUBD_RHS_DENSITY;
  subd_compute_rhs_vel (stokes_problem, subd_options);

  /* about amr */
  {
    rhea_stokes_problem_set_solver_amr (stokes_problem, p4est, discr_options);
    /* perform initial AMR */
    if (p4est != NULL && discr_options != NULL) {
      rhea_stokes_problem_init_amr (stokes_problem, p4est, discr_options);

      /* retrieve adapted mesh */
      ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
      press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);
    }
  }

  /* set up Stokes solver */
  rhea_stokes_problem_setup_solver (stokes_problem);

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

void
adjoint_solve_init (adjoint_problem_t *adjoint_problem)
{
  rhea_stokes_problem_t *stokes_problem = adjoint_problem->stokes_problem;
  p4est_t   *p4est = adjoint_problem->p4est;
  ymir_mesh_t   *ymir_mesh = adjoint_problem->ymir_mesh;
  ymir_pressure_elem_t *press_elem = adjoint_problem->press_elem;
  rhea_discretization_options_t *discr_options = adjoint_problem->discr_options;
  rhea_temperature_options_t  *temp_options = adjoint_problem->temp_options;
  subd_options_t *subd_options = adjoint_problem->subd_options;
  ymir_vec_t *usol = adjoint_problem->usol;

  adjoint_stokes_update_init (stokes_problem, p4est, ymir_mesh, press_elem,
                              discr_options, temp_options, subd_options);
  /* solve Stokes problem*/
  adjoint_run_solver (usol, adjoint_problem);

//  {
//    char *vtk_write_input_path = "/scratch/02600/xiliu/rhea/Adjoint/Test_newton/Test_invertE/Test4/test/input";
//    subd_write_input_basic (stokes_problem, temp_options,
//                            vtk_write_input_path);
//  }
}

void
adjoint_data_init (ymir_vec_t *solution,
                         void *data)
{
  adjoint_problem_t *adjoint_problem = data;
  subd_options_t  *subd_options = adjoint_problem->subd_options;
  ymir_vec_t *msol = adjoint_problem->msol;

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);

  adjoint_parameter_init (solution, msol, subd_options);

  adjoint_solve_init (adjoint_problem);
  adjoint_coeff_update (adjoint_problem);

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

/********************used in newton iterations*********************************/
static void
adjoint_parameter_update (ymir_vec_t *solution, subd_options_t *subd_options)
{
  subd_adjoint_options_t *adjoint_options = subd_options->adjoint_options;
  subd_adjoint_visc_field_t visc_field;
  int            n = adjoint_options->n_components;
  int           field_nums = adjoint_options->visc_fields;
  int           use_exponent = adjoint_options->use_exponent;
  double        m_para;
  int           i;

  for (i = 0; i <  n; i++)  {
    if (!use_exponent)
      m_para = solution->meshfree->e[0][i];
    else
      m_para = exp(solution->meshfree->e[0][i]);

    visc_field = (subd_adjoint_visc_field_t ) adjoint_get_field (&field_nums,  n, i);
    switch (visc_field) {
    case SUBD_ADJOINT_FIELD_VISC_RHEA_UM:
      subd_options->rhea_visc_options->upper_mantle_scaling = m_para;
    break;

    case SUBD_ADJOINT_FIELD_VISC_RHEA_LM:
      subd_options->rhea_visc_options->lower_mantle_scaling = m_para;
    break;

    case SUBD_ADJOINT_FIELD_VISC_RHEA_UM_E:
      subd_options->rhea_visc_options->upper_mantle_arrhenius_activation_energy = m_para;
      break;

    case SUBD_ADJOINT_FIELD_VISC_RHEA_LM_E:
      subd_options->rhea_visc_options->lower_mantle_arrhenius_activation_energy = m_para;
      break;

    case SUBD_ADJOINT_FIELD_VISC_CUSTOM_LAYERS_UM:
      subd_options->visc_options->visc_lith = m_para;

    default:
      RHEA_ABORT_NOT_REACHED ();
    }
  }
}

void
adjoint_stokes_update_forward (rhea_stokes_problem_t *stokes_problem,
                              subd_options_t *subd_options)
{
  ymir_vec_t      *coeff;
  ymir_stokes_op_t *stokes_op;
  ymir_stress_op_t *stress_op;

  ymir_vec_t *rhs_vel;
  ymir_vec_t *rhs_vel_press;

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);

  rhea_stokes_problem_compute_coefficient (stokes_problem, 0 /* !init */);
  coeff = rhea_stokes_problem_get_coeff (stokes_problem);
  stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
  stress_op = stokes_op->stress_op;
  ymir_stress_op_set_coeff_scal (stress_op, coeff);

  subd_options->rhs_type = SUBD_RHS_DENSITY;
  subd_compute_rhs_vel (stokes_problem, subd_options);

  rhs_vel_press = rhea_stokes_problem_get_rhs_vel_press (stokes_problem);
  rhs_vel = rhea_stokes_problem_get_rhs_vel (stokes_problem);
  rhea_stokes_problem_update_solver (
      stokes_problem, NULL,
      1 /* new_rhs */,
      NULL,
      rhs_vel, //rhs_vel: volume forcing is zero
      NULL, // from forward Stokes
      NULL); //rhs_vel_nonzero_dirichlet is zero

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

void
adjoint_solve_forward (adjoint_problem_t *adjoint_problem)
{
  rhea_stokes_problem_t *stokes_problem = adjoint_problem->stokes_problem;
  subd_options_t  *subd_options = adjoint_problem->subd_options;
  ymir_vec_t      *usol = adjoint_problem->usol;

  adjoint_stokes_update_forward (stokes_problem, subd_options);

  adjoint_run_solver (usol, adjoint_problem);
}

void
adjoint_coeff_update (adjoint_problem_t *adjoint_problem)
{
  ymir_mesh_t *ymir_mesh = adjoint_problem->ymir_mesh;
  rhea_stokes_problem_t *stokes_problem = adjoint_problem->stokes_problem;
  subd_options_t  *subd_options = adjoint_problem->subd_options;

  ymir_dvec_t *coeff;
  ymir_dvec_t *stencil = rhea_viscosity_new (ymir_mesh);
  subd_adjoint_options_t *adjoint_options = subd_options->adjoint_options;
  int  n = adjoint_options->n_components;
  int parameter_nums = adjoint_options->visc_parameters;
  subd_adjoint_visc_parameter_t parameter;
  int field_nums = adjoint_options->stencil_options->fields;
  subd_adjoint_stencil_field_t stencil_field;
  int i;

  for (i = 0; i <  n; i++)  {
    coeff = adjoint_problem->coeff[i];
    parameter = adjoint_get_parameter (&parameter_nums,  n, i);
    adjoint_compute_visc_grad_m (coeff, stokes_problem, subd_options, parameter);
    stencil_field = (subd_adjoint_stencil_field_t) adjoint_get_field (&field_nums,  n, i);
RHEA_GLOBAL_PRODUCTIONF ("i=%d: stencil_field=%d\n", i, (int) stencil_field);
    adjoint_setup_stencil (stencil, subd_options, stencil_field);
    ymir_vec_multiply_in (stencil, coeff);
    ymir_vec_scale (2.0, coeff);
  }

  rhea_viscosity_destroy (stencil);
}

/******************* newton function ***********************************/
void
adjoint_update_operator_fn (ymir_vec_t *solution, void *data)
{
  adjoint_problem_t *adjoint_problem = data;
  subd_options_t  *subd_options = adjoint_problem->subd_options;

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);

  adjoint_parameter_update (solution, subd_options);
  adjoint_solve_forward (adjoint_problem);

  adjoint_coeff_update (adjoint_problem);

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

/******************* newton function ***********************************/
double
adjoint_evaluate_objective (ymir_vec_t *solution, void *data)
{
  adjoint_problem_t *adjoint_problem = data;
  double        obj_val;
  ymir_vec_t    *usol = adjoint_problem->usol;
  ymir_mesh_t   *ymir_mesh = adjoint_problem->ymir_mesh;
  ymir_vec_t    *usol_surf = ymir_face_cvec_new (ymir_mesh,
                          RHEA_DOMAIN_BOUNDARY_FACE_TOP, 3);
  ymir_vec_t    *usol_surfZ = ymir_face_cvec_new (ymir_mesh,
                          RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
  sc_dmatrix_t  *usol_surfZ_mat = usol_surfZ->cvec;
  ymir_cvec_t  *mass = ymir_vec_template (usol_surfZ);

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);

  if (solution == NULL) {
    RHEA_GLOBAL_PRODUCTIONF ("set newton_options.nonzero_initial_guess as %d\n", 1);
    return 0;
  }

  ymir_interp_vec (usol, usol_surf);
  ymir_cvec_get_comp (usol_surf, usol_surfZ_mat, 2, YMIR_COPY);

  ymir_mass_apply (usol_surfZ, mass);
  obj_val = ymir_vec_innerprod (usol_surfZ, mass);
  obj_val *= 0.5;
  adjoint_problem->objective = obj_val;


  RHEA_GLOBAL_PRODUCTIONF ("objective is %.12e\n", obj_val);
  ymir_vec_destroy (usol_surf);
  ymir_vec_destroy (usol_surfZ);
  ymir_vec_destroy (mass);

  RHEA_GLOBAL_INFO_FN_END (__func__);

  return (obj_val);
}

void
adjoint_stokes_update_adjoint (ymir_vec_t *usol,
                              rhea_stokes_problem_t *stokes_problem,
                              subd_options_t *subd_options)
{
  ymir_vec_t **rhs_vel_neumann;
  ymir_vec_t *rhs_vel_press;
  ymir_stokes_op_t    *stokes_op;

  subd_options->velbc_options->velbc_nonzero_neu = SUBD_VELBC_NEU_VEC;
  subd_options->data = usol;
  subd_compute_rhs_velbc_neumann (stokes_problem, subd_options);

  rhs_vel_press = rhea_stokes_problem_get_rhs_vel_press (stokes_problem);
  rhs_vel_neumann = rhea_stokes_problem_get_rhs_vel_neumann (stokes_problem);
  stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);

  ymir_stokes_pc_construct_rhs (
      rhs_vel_press, //out and set to stokes
      NULL, //rhs_vel: volume forcing is zero
      rhs_vel_neumann, // from forward Stokes
      NULL, //rhs_vel_nonzero_dirichlet is zero
      1, //incompressible
      stokes_op,
      0);

  subd_options->data = NULL;
}

void
adjoint_solve_adjoint (adjoint_problem_t *adjoint_problem)
{
  ymir_vec_t  *usol = adjoint_problem->usol;
  rhea_stokes_problem_t *stokes_problem = adjoint_problem->stokes_problem;
  subd_options_t *subd_options = adjoint_problem->subd_options;
  ymir_vec_t  *vsol = adjoint_problem->vsol;

  adjoint_stokes_update_adjoint (usol, stokes_problem, subd_options);

  adjoint_run_solver (vsol, adjoint_problem);
}

static void
adjoint_strainrate_dotprod (ymir_vec_t *dotprod,
                    ymir_vec_t *usol, ymir_vec_t *vsol)
{
  ymir_mesh_t *ymir_mesh = usol->mesh;
  ymir_vec_t *edotu = ymir_dvec_new (ymir_mesh, 6, YMIR_GAUSS_NODE);
  ymir_vec_t *edotv = ymir_dvec_new (ymir_mesh, 6, YMIR_GAUSS_NODE);


  ymir_velocity_strain_rate (usol, edotu, 0);
  ymir_velocity_strain_rate (vsol, edotv, 0);

  ymir_velocity_symtens_dotprod (edotu, edotv, dotprod);

  ymir_vec_destroy (edotu);
  ymir_vec_destroy (edotv);
}

void
adjoint_compute_neg_grad (adjoint_problem_t *adjoint_problem)
{
  ymir_vec_t  *usol = adjoint_problem->usol;
  ymir_vec_t  *vsol = adjoint_problem->vsol;
  ymir_vec_t  *neg_grad = adjoint_problem->neg_grad;
  ymir_mesh_t *ymir_mesh = usol->mesh;
  ymir_vec_t *dotprod = ymir_dvec_new(ymir_mesh, 1, YMIR_GAUSS_NODE);
  ymir_dvec_t *dotprod_mass = ymir_vec_template (dotprod);
  ymir_vec_t *coeff;
  int i,  n = adjoint_problem->subd_options->adjoint_options->n_components;

  ymir_vec_set_zero (neg_grad);
  for (i = 0; i <  n; i++)  {
    coeff = adjoint_problem->coeff[i];
    adjoint_strainrate_dotprod (dotprod, usol, vsol);
    ymir_mass_apply (dotprod, dotprod_mass);
    neg_grad->meshfree->e[0][i] = -1.0 * ymir_vec_innerprod (coeff, dotprod_mass);
  }

  ymir_vec_destroy (dotprod);
  ymir_vec_destroy (dotprod_mass);
}

/******************* newton function ***********************************/
void
adjoint_solve_negative_gradient (ymir_vec_t *neg_grad,
                                   ymir_vec_t *solution, void *data)
{
  adjoint_problem_t *adjoint_problem = data;
  rhea_stokes_problem_t *stokes_problem = adjoint_problem->stokes_problem;
  subd_options_t *subd_options = adjoint_problem->subd_options;
  int n = subd_options->adjoint_options->n_components;

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);

  if (solution == NULL) {
    RHEA_GLOBAL_PRODUCTIONF ("set newton_options.nonzero_initial_guess as %d\n", 1);
    return;
  }

  adjoint_solve_adjoint (adjoint_problem);

  adjoint_compute_neg_grad (adjoint_problem);
  ymir_vec_copy (adjoint_problem->neg_grad, neg_grad);

  for (int i = 0; i < n; i++) {
    RHEA_GLOBAL_PRODUCTIONF ("negative_gradient%d = %.12e\n", i, neg_grad->meshfree->e[0][i]);
  }

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

/******************* newton function ***********************************/
double
adjoint_compute_gradient_norm (ymir_vec_t *neg_grad,
                               void *data, double *norm_comp)
{
  double            norm;
  adjoint_problem_t *adjoint_problem = data;

  norm = ymir_vec_norm (neg_grad);

  return (norm);
}

void
adjoint_stokes_update_hessian_forward (adjoint_problem_t *adjoint_problem, int i)
{
  subd_options_t        *subd_options = adjoint_problem->subd_options;
  ymir_mesh_t           *ymir_mesh = adjoint_problem->ymir_mesh;
  rhea_stokes_problem_t *stokes_problem = adjoint_problem->stokes_problem;
  ymir_vec_t            *usol = adjoint_problem->usol;
  ymir_vec_t            *coeff;
  ymir_vec_t            *rhs_vel_press;

  ymir_stokes_op_t      *stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
  ymir_stress_op_t      *stress_op;
  ymir_vel_dir_t        *vel_dir = rhea_stokes_problem_get_vel_dir (stokes_problem);

  ymir_vec_t            *rhs_vel = rhea_velocity_new (ymir_mesh);

  subd_adjoint_options_t *adjoint_options = subd_options->adjoint_options;

  coeff = adjoint_problem->coeff[i];
  stress_op = ymir_stress_op_new (coeff, vel_dir, NULL, NULL, NULL);
  ymir_stress_pc_apply_stress_op (usol, rhs_vel, stress_op, 0, 0);
  ymir_vec_scale (-1.0, rhs_vel);

  rhs_vel_press = rhea_stokes_problem_get_rhs_vel_press (stokes_problem);
  ymir_vec_set_zero (rhs_vel_press);
  ymir_upvec_set_u (rhs_vel_press, rhs_vel, YMIR_SET);

  ymir_stress_op_destroy (stress_op);
  rhea_velocity_destroy (rhs_vel);
}

void
adjoint_solve_hessian_forward (adjoint_problem_t *adjoint_problem, int i)
{
  ymir_vec_t *Husol = adjoint_problem->Husol[i];

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);
  RHEA_GLOBAL_PRODUCTIONF ("component_id = %d\n", i);

  adjoint_stokes_update_hessian_forward (adjoint_problem, i);
  adjoint_run_solver (Husol, adjoint_problem);

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

void
adjoint_stokes_update_hessian_adjoint (ymir_vec_t *Husol,
                                      rhea_stokes_problem_t *stokes_problem,
                                      subd_options_t  *subd_options)
{
  ymir_vec_t    *rhs_vel_press;
  ymir_vec_t    **rhs_vel_neumann;
  ymir_stokes_op_t  *stokes_op;

  subd_options->rhs_type = SUBD_RHS_ADJOINT_NONE;
  subd_compute_rhs_vel (stokes_problem, subd_options);

  subd_options->velbc_options->velbc_nonzero_neu = SUBD_VELBC_NEU_VEC;
  subd_options->data = Husol;
  subd_compute_rhs_velbc_neumann (stokes_problem, subd_options);

  rhs_vel_press = rhea_stokes_problem_get_rhs_vel_press (stokes_problem);
  rhs_vel_neumann = rhea_stokes_problem_get_rhs_vel_neumann (stokes_problem);
  stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);

  ymir_stokes_pc_construct_rhs (
    rhs_vel_press, //out and set to stokes
    NULL, //rhs_vel: volume forcing is zero
    rhs_vel_neumann, // from stokes
    NULL, //rhs_vel_nonzero_dirichlet is zero
    1, //incompressible
    stokes_op,
    0);
}

void
adjoint_solve_hessian_adjoint (adjoint_problem_t *adjoint_problem, int i)
{
  rhea_stokes_problem_t *stokes_problem = adjoint_problem->stokes_problem;
  subd_options_t *subd_options = adjoint_problem->subd_options;

  ymir_vec_t *Husol = adjoint_problem->Husol[i];
  ymir_vec_t *Hvsol = adjoint_problem->Hvsol[i];

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);
  RHEA_GLOBAL_PRODUCTIONF ("component_id = %d\n", i);

  adjoint_stokes_update_hessian_adjoint (Husol, stokes_problem, subd_options);
  adjoint_run_solver (Hvsol, adjoint_problem);

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

void
adjoint_compute_hessian (ymir_vec_t *hessian,
                            adjoint_problem_t *adjoint_problem)
{
  ymir_vec_t *usol = adjoint_problem->usol;
  ymir_mesh_t *ymir_mesh = usol->mesh;
  subd_adjoint_options_t *adjoint_options = adjoint_problem->subd_options->adjoint_options;
  int n = adjoint_options->n_components;

  ymir_vec_t  **coeff = adjoint_problem->coeff;
  ymir_vec_t  **Hvsol = adjoint_problem->Hvsol;

  ymir_dvec_t *dotprod = ymir_dvec_new(ymir_mesh, 1, YMIR_GAUSS_NODE);
  ymir_vec_t  **dotprod_mass = YMIR_ALLOC (ymir_vec_t *, n);

  int i, j, k = 0;

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);

  for (j = 0; j < n; j++)  {
    dotprod_mass[j] = ymir_dvec_new(ymir_mesh, 1, YMIR_GAUSS_NODE);
    adjoint_strainrate_dotprod (dotprod, usol, Hvsol[j]);
    ymir_mass_apply (dotprod, dotprod_mass[j]);
  }

  for (i = 0; i < n; i++)  {
    for (j = 0; j < n; j++) {
      hessian->meshfree->e[0][k] = ymir_vec_innerprod (coeff[i], dotprod_mass[j]);
      k++;
    }
  }

  for (i = 0; i < n; i++) {
    RHEA_GLOBAL_PRODUCTIONF ("diagonal hessian%d = %.12e\n", i, hessian->meshfree->e[0][(i+1)*(i+1) - 1]);
  }
  ymir_vec_destroy (dotprod);
  for (i = 0; i < n; i++) {
    ymir_vec_destroy (dotprod_mass[i]);
  }
  YMIR_FREE (dotprod_mass);

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

void
adjoint_compute_step (ymir_vec_t *step,
                        adjoint_problem_t *adjoint_problem)
{
  int n = adjoint_problem->subd_options->adjoint_options->n_components;
  double *hessian_data = adjoint_problem->hessian->meshfree->e[0];
  double *neg_grad_data = adjoint_problem->neg_grad->meshfree->e[0];
  double *step_data = step->meshfree->e[0];
  double **hessian_matrix;
  int k = 0;

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);

  hessian_matrix = allocate_square_matrix (n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      hessian_matrix[i][j] = hessian_data[k];
      k++;
    }
  }

  solve_x_from_Ax_b (step_data, hessian_matrix, neg_grad_data, n);
  destroy_square_matrix (hessian_matrix, n);

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

/******************* newton function ***********************************/
int
adjoint_solve_hessian_system_fn (ymir_vec_t *step, ymir_vec_t *neg_grad,
                                const int lin_iter_max, const double lin_res_norm_rtol,
                                const int nonzero_initial_guess, void *data,
                                int *lin_iter_count)
{
  adjoint_problem_t *adjoint_problem = data;
  rhea_stokes_problem_t *stokes_problem = adjoint_problem->stokes_problem;
  subd_options_t *subd_options = adjoint_problem->subd_options;
  ymir_vec_t *hessian = adjoint_problem->hessian;
  int n = adjoint_problem->subd_options->adjoint_options->n_components;

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);

   /*solve for Hessian forward*/
  for (int i = 0; i < n; i++)  {
    adjoint_solve_hessian_forward (adjoint_problem, i);
    adjoint_solve_hessian_adjoint (adjoint_problem, i);
  }
  ymir_vec_set_zero (hessian);
  adjoint_compute_hessian (hessian, adjoint_problem);

  adjoint_compute_step (step, adjoint_problem);

  /* return iteraton count and "stopping" reason */
  if (lin_iter_count != NULL) {
    *lin_iter_count = 1;
  }

  RHEA_GLOBAL_INFO_FN_END (__func__);

  return 1;
}

/******************* newton function ***********************************/
void
adjoint_apply_hessian_fn (ymir_vec_t *hessian_vec, ymir_vec_t *dir_vec,
                         void *data)
{
   adjoint_problem_t *adjoint_problem = data;

   ymir_vec_copy (adjoint_problem->hessian, hessian_vec);
   ymir_vec_multiply_in (dir_vec, hessian_vec);
}

/******************* newton function ***********************************/
void
adjoint_update_hessian_fn (ymir_vec_t *solution, ymir_vec_t *step,
                          const double step_length, void *data)
{
  adjoint_problem_t   *adjoint_problem = (adjoint_problem_t *) data;
  int  n = adjoint_problem->subd_options->adjoint_options->n_components;
  int  use_exponent = adjoint_problem->subd_options->adjoint_options->use_exponent;


  ymir_vec_copy (solution, adjoint_problem->msol);
  for (int i = 0; i < n; i++)  {
    if (use_exponent)
      RHEA_GLOBAL_PRODUCTIONF ("solution%d = %.12e\n", i, exp(solution->meshfree->e[0][i]));
    else
      RHEA_GLOBAL_PRODUCTIONF ("solution%d = %.12e\n", i, solution->meshfree->e[0][i]);

  }
}

void
adjoint_output_prestep (ymir_vec_t *solution, const int iter, void *data)
{
  adjoint_problem_t *adjoint_problem = data;
  rhea_stokes_problem_t *stokes_problem = adjoint_problem->stokes_problem;
  subd_options_t    *subd_options = adjoint_problem->subd_options;
  ymir_mesh_t *ymir_mesh = adjoint_problem->ymir_mesh;
  ymir_pressure_elem_t *press_elem = adjoint_problem->press_elem;
  ymir_vec_t *sol_vel_press = adjoint_problem->sol_vel_press;
  const char *solution_path = subd_options->adjoint_options->vtk_write_solution_path;

  double      objective = adjoint_problem->objective;
  ymir_vec_t *neg_grad = adjoint_problem->neg_grad;
  ymir_vec_t *hessian  = adjoint_problem->hessian;
  FILE *fp = adjoint_problem->fp_record;
  int         use_exponent = subd_options->adjoint_options->use_exponent;
  int         n = subd_options->adjoint_options->n_components;
  double      m_para[n];


  for (int i = 0; i < n; i++) {
    if (!use_exponent)
      m_para[i] = solution->meshfree->e[0][i];
    else
      m_para[i] = exp(solution->meshfree->e[0][i]);
  }

  if (iter == 0) {
    if (solution_path) {
      char   path[BUFSIZ];
      snprintf (path, BUFSIZ, "%s_init", solution_path);
      subd_vtk_write_solution (ymir_mesh, press_elem, stokes_problem, sol_vel_press);
    }
     if (fp) {
       fprintf (fp, "%4s %16s %16s %16s %16s %16s %16s %16s %.16s %.16s\n", "iter", "parameter1", "parameter2", "objective", "gradient1", "gradient2", "hessian00", "hessian01", "hessian10", "hessian11");
       fprintf (fp, "%4d %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n", iter, m_para[0], m_para[1], objective, -1.0 * neg_grad->meshfree->e[0][0], -1.0 * neg_grad->meshfree->e[0][1], hessian->meshfree->e[0][0], hessian->meshfree->e[0][1], hessian->meshfree->e[0][2],hessian->meshfree->e[0][3]);
     }
  }
  else {
    if (fp) {
       fprintf (fp, "%4d %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n", iter, m_para[0], m_para[1], objective, -1.0 * neg_grad->meshfree->e[0][0], -1.0 * neg_grad->meshfree->e[0][1], hessian->meshfree->e[0][0], hessian->meshfree->e[0][1], hessian->meshfree->e[0][2],hessian->meshfree->e[0][3]);
    }
  }

}

//TODO obsolete
double
subd_adjoint_gradient (ymir_vec_t *edot0, ymir_vec_t *edot1,
                            ymir_vec_t *temp, ymir_vec_t *visc)
{
  ymir_mesh_t *ymir_mesh = ymir_vec_get_mesh (edot0);
  ymir_dvec_t *dotprod = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
  ymir_dvec_t *integrant = ymir_vec_template (dotprod);
  ymir_dvec_t *vec_e = ymir_vec_template (dotprod);
  ymir_dvec_t *masse = ymir_vec_template (dotprod);

  const int           N  = ymir_n (ymir_mesh->cnodes->N);
  const int           Np = ymir_np (ymir_mesh->cnodes->N);
  const int           K  = ymir_mesh->cnodes->K;
  sc_dmatrix_t       *elemtemp = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemvisc = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemdotprod = sc_dmatrix_new (1, Np);
  sc_dmatrix_t       *elemintegrant = sc_dmatrix_new (1, Np);

  int     elid, gp;
  double integral = .0;

  ymir_velocity_symtens_dotprod (edot0, edot1, dotprod);

  for (elid = 0; elid < K; elid++) {
    ymir_cvec_get_elem_interp (temp, elemtemp, YMIR_STRIDE_NODE, elid,
                               YMIR_GAUSS_NODE, YMIR_COPY);
    ymir_dvec_get_elem (visc, elemvisc, YMIR_STRIDE_NODE, elid,
                               YMIR_COPY);
    ymir_dvec_get_elem (dotprod, elemdotprod, YMIR_STRIDE_NODE, elid,
                               YMIR_COPY);
    ymir_dvec_get_elem (integrant, elemintegrant, YMIR_STRIDE_NODE, elid,
                               YMIR_WRITE);
    double     *_sc_restrict tempd = elemtemp->e[0];
    double     *_sc_restrict viscd = elemvisc->e[0];
    double     *_sc_restrict dotprodd = elemdotprod->e[0];
    double     *_sc_restrict integrantd = elemintegrant->e[0];


    for (gp = 0; gp < Np; gp++) {
      integrantd[gp] = 2.0 * viscd[gp] * dotprodd[gp];

    }
    ymir_dvec_set_elem (integrant, elemintegrant, YMIR_STRIDE_NODE, elid,
                        YMIR_SET);
  }
  ymir_vec_set_value (vec_e, 1.0);
  ymir_mass_apply (vec_e, masse);
  integral = ymir_vec_innerprod (integrant, masse);
//  integral = 2.0 * ymir_vec_innerprod (dotprod, masse);

  sc_dmatrix_destroy (elemtemp);
  sc_dmatrix_destroy (elemvisc);
  sc_dmatrix_destroy (elemdotprod);
  sc_dmatrix_destroy (elemintegrant);
  ymir_vec_destroy (dotprod);
  ymir_vec_destroy (integrant);
  ymir_vec_destroy (vec_e);
  ymir_vec_destroy (masse);

  return (integral);
}


