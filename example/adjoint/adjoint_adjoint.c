#include <adjoint_adjoint.h>

void
adjoint_data_initialize (ymir_vec_t *solution,
                         void *data)
{
  const char    *this_fn_name = "adjoint_data_initialize";
  adjoint_problem_t *adjoint_problem = data;
  rhea_stokes_problem_t *stokes_problem = adjoint_problem->stokes_problem;
  p4est_t   *p4est = adjoint_problem->p4est;
  ymir_mesh_t   *ymir_mesh = adjoint_problem->ymir_mesh;
  ymir_pressure_elem_t *press_elem = adjoint_problem->press_elem;
  rhea_discretization_options_t *discr_options = adjoint_problem->discr_options;
  rhea_temperature_options_t  *temp_options = adjoint_problem->temp_options;
  subd_options_t  *subd_options = adjoint_problem->subd_options;
  const char *vtk_write_input_path = adjoint_problem->vtk_write_input_path;

  ymir_vec_t *sol_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);
  ymir_vec_t *usol = rhea_velocity_new (ymir_mesh);

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);
  int nonzero_initial_guess = 0;
  int solver_iter_max = adjoint_problem->solver_iter_max;
  double solver_rel_tol = adjoint_problem->solver_rel_tol;

  //TODO: modify subd_options to recompute viscosity, also reassign rhs_vel, rhs_nonzero Dirichlet and Neumann
  adjoint_step_initialize_stokes (&stokes_problem, p4est, &ymir_mesh, &press_elem,
                    discr_options, temp_options, subd_options,
                    vtk_write_input_path);

  /* solve Stokes problem*/
  subd_run_solver (sol_vel_press, ymir_mesh, press_elem, stokes_problem,
                   nonzero_initial_guess, solver_iter_max, solver_rel_tol);

  ymir_upvec_get_u (sol_vel_press, usol, YMIR_COPY);
  adjoint_problem->usol = usol;
  ymir_vec_destroy (sol_vel_press);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

double
adjoint_evaluate_objective (ymir_vec_t *solution, void *data)
{
  const char    *this_fn_name = "adjoint_evaluate_objective";
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

  ymir_interp_vec (usol, usol_surf);
  ymir_cvec_get_comp (usol_surf, usol_surfZ_mat, 2, YMIR_COPY);

  ymir_mass_apply (usol_surfZ, mass);
  obj_val = ymir_vec_innerprod (usol_surfZ, mass);

  ymir_vec_destroy (usol_surf);
  ymir_vec_destroy (usol_surfZ);
  ymir_vec_destroy (mass);

  return (obj_val);
}
void adjoint_compute_gradient_norm (ymir_vec_t *neg_gradient,
                                    void *data, double *norm)
{
  adjoint_problem_t *adjoint_problem = data;
  adjoint_compute_negative_gradient (neg_gradient, NULL, adjoint_problem);
  norm = neg_gradient->meshfree->e[0];
}

void
adjoint_compute_negative_gradient (ymir_vec_t *neg_gradient,
                                   ymir_vec_t *solution, void *data)
{
  const char      *this_fn_name = "adjoint_compute_negative_gradient";
  adjoint_problem_t *adjoint_problem = data;
  ymir_vec_t  *usol = adjoint_problem->usol;
  subd_options_t *subd_options = adjoint_problem->subd_options;
  rhea_stokes_problem_t *stokes_problem = adjoint_problem->stokes_problem;
  ymir_stokes_op_t    *stokes_op;
  ymir_vec_t **rhs_vel_neumann;
  ymir_vec_t *rhs_vel_press;

  ymir_mesh_t   *ymir_mesh = adjoint_problem->ymir_mesh;
  ymir_pressure_elem_t *press_elem = adjoint_problem->press_elem;

  ymir_vec_t *sol_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);
  ymir_vec_t *vsol = rhea_velocity_new (ymir_mesh);

  int nonzero_initial_guess = 0;
  int solver_iter_max = adjoint_problem->solver_iter_max;
  double solver_rel_tol = adjoint_problem->solver_rel_tol;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  subd_options->rhs_type = (subd_rhs_type_t) 2;
  subd_options->velbc_options->velbc_nonzero_neu = (subd_velbc_nonzero_neu_t) 3;
  subd_options->velo = adjoint_problem->usol;

  subd_compute_rhs_vel (stokes_problem, subd_options);
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

  /* solve Stokes problem*/
  subd_run_solver (sol_vel_press, ymir_mesh, press_elem, stokes_problem,
                   nonzero_initial_guess, solver_iter_max, solver_rel_tol);

  ymir_upvec_get_u (sol_vel_press, vsol, YMIR_COPY);
  adjoint_problem->vsol = vsol;
  ymir_vec_destroy (sol_vel_press);
  neg_gradient->meshfree->e[0] = adjoint_vec_gradient (usol, vsol);
exit(10);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

void
adjoint_setup_adjoint_problem (adjoint_problem_t *adjoint_problem,
                              rhea_stokes_problem_t *stokes_problem,
                              p4est_t *p4est, ymir_mesh_t *ymir_mesh,
                              ymir_pressure_elem_t *press_elem,
                              rhea_discretization_options_t *discr_options,
                              rhea_temperature_options_t *temp_options,
                              subd_options_t *subd_options,
                              const char *vtk_write_input_path,
                              int       solver_iter_max,
                              double    solver_rel_tol)
{
  adjoint_problem->stokes_problem = stokes_problem;
  adjoint_problem->p4est = p4est;
  adjoint_problem->ymir_mesh = ymir_mesh;
  adjoint_problem->press_elem = press_elem;
  adjoint_problem->discr_options = discr_options;
  adjoint_problem->temp_options = temp_options;
  adjoint_problem->subd_options = subd_options;
  adjoint_problem->vtk_write_input_path = vtk_write_input_path;

  adjoint_problem->solver_iter_max = solver_iter_max;
  adjoint_problem->solver_rel_tol = solver_rel_tol;

  adjoint_problem->vsol = NULL;
  adjoint_problem->Husol = NULL ;
  adjoint_problem->Hvsol = NULL;
  adjoint_problem->sol_vel_press = NULL;
}

void
adjoint_setup_newton (rhea_newton_problem_t **newton_problem,
                      adjoint_problem_t *adjoint_problem)
{
  char      *this_fn_name = "adjoint_setup_newton";

  RHEA_GLOBAL_PRODUCTIONF ("into %s\n", this_fn_name);

  /*create newton problem */
  {
    ymir_vec_t    *neg_gradient_vec = ymir_vec_new_meshfree (1);
    ymir_vec_t    *step_vec = ymir_vec_template (neg_gradient_vec);

    *newton_problem = rhea_newton_problem_new (
        adjoint_compute_negative_gradient,
//        NULL,
//        adjoint_solve_hessian_system_fn);
        NULL);
    rhea_newton_problem_set_vectors (
        neg_gradient_vec, step_vec, *newton_problem);
    rhea_newton_problem_set_data_fn (
        adjoint_problem, adjoint_data_initialize, NULL /* no clear fnc. */,
        *newton_problem);
    rhea_newton_problem_set_conv_criterion_fn (
        RHEA_NEWTON_CONV_CRITERION_GRADIENT_NORM,
        adjoint_evaluate_objective,
//        NULL, // Is that OK?
        adjoint_compute_gradient_norm,
        0 /* no multi-component norms */, *newton_problem);
//    rhea_newton_problem_set_apply_hessian_fn (
//        adjoint_apply_hessian_fn, *newton_problem);
//    rhea_newton_problem_set_update_fn (
//        adjoint_update_operator_fn /* updating linearized operator is not necessary */,
//        adjoint_update_hessian_fn,
//        adjoint_modify_hessian_system_fn /* updating RHS of linearized system is not necessary */,
//        *newton_problem);
  }
#ifdef RHEA_ENABLE_DEBUG
  rhea_newton_problem_set_checks (1 /* grad */, 0 /* Hessian */, *nl_problem);
#endif


  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", __func__);
}



/**************************************************************
 * Adjoint related
***************************************************************/
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
  double              stencil_um = subd_options->adjoint_options->stencil_options->value;
  double              stencil_lm = .0;
  int                 m;

  /* compute viscosity in this element */
    for (nodeid = 0; nodeid < n_nodes; nodeid++) {
      if (z[nodeid] > z_mid)  {
        visc_elem[nodeid] = stencil_um;
      }
      else {
        visc_elem[nodeid] = stencil_lm;
      }
    }
    /* check viscosity for `nan`, `inf`, and positivity */
    RHEA_ASSERT (isfinite (visc_elem[nodeid]));
}

void
subd_adjoint_stencil_visc  (ymir_vec_t *stencil_visc,
                         subd_options_t *subd_options)
{
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (stencil_visc);
  const ymir_locidx_t  n_elements = ymir_mesh_get_num_elems_loc (mesh);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);
  mangll_t           *mangll = mesh->ma;
  const int           N = ymir_n (mangll->N);
  char                *this_fn_name = "subd_adjoint_visc_stencil";

  sc_dmatrix_t       *stencil_el_mat;
  double             *x, *y, *z, *tmp_el,*stencil_el_data;
  ymir_locidx_t       elid;

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);
  /* create work variables */
  stencil_el_mat = sc_dmatrix_new (n_nodes_per_el, 1);
  stencil_el_data = stencil_el_mat->e[0];
  x = RHEA_ALLOC (double, n_nodes_per_el);
  y = RHEA_ALLOC (double, n_nodes_per_el);
  z = RHEA_ALLOC (double, n_nodes_per_el);
  tmp_el = RHEA_ALLOC (double, n_nodes_per_el);

  for (elid = 0; elid < n_elements; elid++) {
    /* get coordinates of this element at Gauss nodes */
    ymir_mesh_get_elem_coord_gauss (x, y, z, elid, mesh, tmp_el);

    /* compute stencilosity stencil*/ //TODO: later on add visc_rhea, and other viscosity
    if (subd_options->adjoint_options->stencil_options->type == SUBD_ADJOINT_STENCIL_VISC_CUSTOM) {
      switch (subd_options->adjoint_options->stencil_options->custom_type) {
        case SUBD_ADJOINT_STENCIL_VISC_CUSTOM_LAYERS_UM:
          subd_adjoint_stencil_visc_custom_layers_um_elem (stencil_el_data, x, y, z, n_nodes_per_el,
                                           subd_options);
        break;

        default:
        RHEA_ABORT_NOT_REACHED ();
      }
    }
    /* set viscosity of this element */
    rhea_viscosity_set_elem_gauss (stencil_visc, stencil_el_mat, elid);
  }

  /* destroy */
  sc_dmatrix_destroy (stencil_el_mat);
  RHEA_FREE (x);
  RHEA_FREE (y);
  RHEA_FREE (z);
  RHEA_FREE (tmp_el);
}

void
subd_adjoint_rhs_hessian_forward (ymir_vec_t *rhs_vel_press,
                                    rhea_stokes_problem_t *stokes_problem,
                                    subd_options_t  *subd_opt)
{
  rhea_domain_options_t        *domain_opt = subd_opt->domain_options;
  ymir_mesh_t           *mesh = rhs_vel_press->mesh;
  ymir_vec_t            *stencil_visc = rhea_viscosity_new (mesh);
  ymir_vec_t            *viscosity = rhea_viscosity_new (mesh);
  ymir_stress_op_t      *stress_op;
  ymir_vel_dir_t        *vel_dir = rhea_stokes_problem_get_vel_dir (stokes_problem);

  ymir_vec_t            *velo = rhea_velocity_new (mesh);
  ymir_vec_t            *rhs_vel = rhea_velocity_new (mesh);
  char                  *this_fn_name = "subd_adjoint_rhs_hessian_forward";
  char                  *filepath = "/scratch/02600/xiliu/rhea/Subduction/Test_adjoint/Test_loop/topo/invert_um/Invert2/iter1_forward/velo";

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);
  subd_facevec_read (velo, filepath);
  rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);
  subd_adjoint_stencil_visc (stencil_visc, subd_opt);
  ymir_vec_multiply_in (stencil_visc,viscosity);

  stress_op = ymir_stress_op_new (viscosity, vel_dir, NULL, velo, NULL);
  ymir_stress_pc_apply_stress_op (velo, rhs_vel, stress_op, 0, 0);
  ymir_vec_scale (-1.0, rhs_vel);

  ymir_vec_set_zero (rhs_vel_press);
  ymir_upvec_set_u (rhs_vel_press, rhs_vel, YMIR_SET);

  ymir_stress_op_destroy (stress_op);
  rhea_viscosity_destroy (stencil_visc);
  rhea_viscosity_destroy (viscosity);
  rhea_velocity_destroy (velo);
  rhea_velocity_destroy (rhs_vel);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

double *
adjoint_vec_gradient (ymir_vec_t *usol, ymir_vec_t *vsol)
{
  ymir_mesh_t *ymir_mesh = usol->mesh;
  ymir_vec_t *edotu = ymir_dvec_new (ymir_mesh, 6, YMIR_GAUSS_NODE);
  ymir_vec_t *edotv = ymir_dvec_new (ymir_mesh, 6, YMIR_GAUSS_NODE);
  ymir_vec_t *dotprod = ymir_dvec_new(ymir_mesh, 1, YMIR_GAUSS_NODE);
  ymir_dvec_t *vec_e = ymir_vec_template (dotprod);
  ymir_dvec_t *masse = ymir_vec_template (dotprod);
  double  *integral = RHEA_ALLOC (double, 1);

  ymir_velocity_strain_rate (usol, edotu, 0);
  ymir_velocity_strain_rate (vsol, edotv, 0);
  ymir_velocity_symtens_dotprod (edotu, edotv, dotprod);

  ymir_vec_set_value (vec_e, 1.0);
  ymir_mass_apply (vec_e, masse);
  integral[0] = ymir_vec_innerprod (dotprod, masse);

  return (integral);
}

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
