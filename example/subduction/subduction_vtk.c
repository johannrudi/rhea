#include <subduction_vtk.h>

char               *vtk_write_solution_path;
char               *vtk_write_stress_path;
char               *vtk_write_freesurface_path;
char               *vtk_write_postp_path;
char               *vtk_write_test_path;
char               *vtk_write_neumann_path;

void subduction_add_vtk_options (ymir_options_t * opt)
{
  ymir_options_addv (opt,

  /* vtk output options */
  YMIR_OPTIONS_S, "vtk-write-solution-path", '\0',
    &(vtk_write_solution_path), NULL,
    "File path for vtk files for the solution of the Stokes problem",

  YMIR_OPTIONS_S, "vtk-write-stress-path", '\0',
    &(vtk_write_stress_path), NULL,
    "File path for vtk files for stress results",

  YMIR_OPTIONS_S, "vtk-write-freesurface-path", '\0',
    &(vtk_write_freesurface_path), NULL,
    "File path for vtk files for freesurface results",

  YMIR_OPTIONS_S, "vtk-write-postp-path", '\0',
    &(vtk_write_postp_path), NULL,
    "File path for vtk files of post processing data",

  YMIR_OPTIONS_S, "vtk-write-test-path", '\0',
    &(vtk_write_test_path), NULL,
    "File path for vtk files of test: stress_op or manufactured",

  YMIR_OPTIONS_S, "vtk-write-neumann-path", '\0',
    &(vtk_write_neumann_path), NULL,
    "File path for the neumann",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

}

/* Write vtk of input data. */
void
subd_write_input_basic (rhea_stokes_problem_t *stokes_problem,
                       rhea_temperature_options_t *temp_options,
                       rhea_viscosity_options_t *visc_options,
                       const char *vtk_write_input_path)
{
  const char         *this_fn_name = "subd_write_input_basic";
  ymir_mesh_t        *ymir_mesh;
  ymir_vec_t         *background_temp, *viscosity, *bounds_marker;
  ymir_vec_t         *temperature, *weakzone, *rhs_vel;
  char                path[BUFSIZ];

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* get mesh */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);

  temperature = rhea_stokes_problem_get_temperature (stokes_problem);
  weakzone = rhea_stokes_problem_get_weakzone (stokes_problem);
  rhs_vel = rhea_stokes_problem_get_rhs_vel (stokes_problem);

  viscosity = rhea_viscosity_new (ymir_mesh);
  bounds_marker = rhea_viscosity_new (ymir_mesh);
  switch (visc_options->type) {
  case RHEA_VISCOSITY_LINEAR:
    rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);
    break;
  case RHEA_VISCOSITY_NONLINEAR:
    rhea_viscosity_compute_nonlinear_init (viscosity,
                                           NULL /* nl. Stokes output */,
                                           bounds_marker,
                                           temperature, weakzone,
                                           visc_options);
    break;
  default: /* unknown viscosity type */
    RHEA_ABORT_NOT_REACHED ();
  }


  background_temp = rhea_temperature_new (ymir_mesh);
  rhea_temperature_background_compute (background_temp, temp_options);

  RHEA_ASSERT (rhs_vel != NULL);
  RHEA_ASSERT (rhs_vel->meshnum == YMIR_VOL_MESH);

  rhea_vtk_write_input_data (vtk_write_input_path, temperature,
                             background_temp, weakzone, viscosity, bounds_marker,
                             NULL, rhs_vel);

  rhea_temperature_destroy (background_temp);
  rhea_viscosity_destroy (viscosity);
  rhea_viscosity_destroy (bounds_marker);

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}


/* Write vtk of input data. */
void
subd_write_input (ymir_mesh_t *ymir_mesh,
                   rhea_stokes_problem_t *stokes_problem,
                   rhea_temperature_options_t *temp_options,
                   ymir_vec_t *temperature,
                   ymir_vec_t *visc_TI_svisc,
                   ymir_vec_t *visc_TI_rotate,
                   const char *vtk_write_input_path)
{
  const char         *this_fn_name = "subd_write_input";
  ymir_vec_t         *background_temp = rhea_temperature_new (ymir_mesh);
  ymir_vec_t         *weakzone, *rhs_vel, *rhs_vel_nonzero_dirichlet;
  char                path[BUFSIZ];

  RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  weakzone = rhea_stokes_problem_get_weakzone (stokes_problem);
  rhs_vel = rhea_stokes_problem_get_rhs_vel (stokes_problem);
  rhea_temperature_background_compute (background_temp, temp_options);

  RHEA_ASSERT (rhs_vel != NULL);
  RHEA_ASSERT (rhs_vel->meshnum == YMIR_VOL_MESH);

      ymir_vtk_write (ymir_mesh, vtk_write_input_path,
                    temperature, "temperature",
                    rhs_vel, "rhs_vel",
                    background_temp, "background_temp",
                    weakzone, "weakzone",
                    NULL);


  rhea_temperature_destroy (background_temp);

  /* output transversely isotropy parameters*/
  {
    if (visc_TI_svisc != NULL && visc_TI_rotate != NULL) {
      ymir_vec_t         *shear_visc = ymir_vec_clone (visc_TI_svisc);
      ymir_vec_scale (0.5, shear_visc);
      ymir_vec_t *angle = ymir_vec_clone (visc_TI_rotate);
      ymir_vec_scale (180.0/M_PI, angle);
      snprintf (path, BUFSIZ, "%s_anisotropic_viscosity", vtk_write_input_path);
      ymir_vtk_write (ymir_mesh, path,
                      shear_visc, "shear_viscosity",
                      angle, "rotation",
                      NULL);
      rhea_viscosity_destroy (shear_visc);
      rhea_viscosity_destroy (angle);
    }
  }

 // {
 //   rhs_vel_nonzero_dirichlet = rhea_stokes_problem_get_rhs_vel_nonzero_dirichlet (stokes_problem);
 //   if (rhs_vel_nonzero_dirichlet != NULL)
 //     snprintf (path, BUFSIZ, "%s_nonzero_dirichlet", vtk_write_input_path);
 //     ymir_vtk_write (ymir_mesh, path,
 //                     rhs_vel_nonzero_dirichlet, "velocity B.C.",
 //                     NULL);
 // }

  RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

void
subd_vtk_write (rhea_stokes_problem_t *stokes_problem,
                subd_options_t *subd_options)
{
  const char * this_fn_name = "subd_vtk_write";
  ymir_mesh_t         *ymir_mesh
                          = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  ymir_pressure_elem_t   *press_elem
                          = rhea_stokes_problem_get_press_elem (stokes_problem);
  ymir_vec_t          *sol_vel_press
                          = rhea_stokes_problem_get_velocity_pressure (stokes_problem);

  subd_vtk_write_solution (ymir_mesh, press_elem, stokes_problem, sol_vel_press);

  subd_vtk_write_test (ymir_mesh, press_elem, stokes_problem, sol_vel_press, subd_options);

  subd_vtk_write_postp (ymir_mesh, press_elem, stokes_problem, sol_vel_press, subd_options);

  subd_vtk_write_stress (ymir_mesh, press_elem, stokes_problem, sol_vel_press, subd_options);

  subd_vtk_write_freesurface (ymir_mesh, press_elem, sol_vel_press);

  subd_vtk_write_neumann (ymir_mesh, press_elem, stokes_problem, sol_vel_press);

}

void
subd_vtk_write_solution (ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         rhea_stokes_problem_t *stokes_problem,
                         ymir_vec_t *sol_vel_press)
{
  if (vtk_write_solution_path != NULL)  {
  /* write vtk of solution */
    ymir_vec_t         *velocity = rhea_velocity_new (ymir_mesh);
    ymir_vec_t         *pressure = rhea_pressure_new (ymir_mesh, press_elem);
    ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);

    ymir_stokes_vec_get_components (sol_vel_press, velocity, pressure,
                                    press_elem);
    rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

  //DEPRECATED
  //rhea_vtk_write_solution (vtk_write_solution_path, velocity, pressure,
  //                         viscosity, NULL, NAN);

    rhea_pressure_destroy (pressure);
    rhea_viscosity_destroy (viscosity);
    rhea_velocity_destroy (velocity);
  }
}

void
subd_vtk_write_test (ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         rhea_stokes_problem_t *stokes_problem,
                         ymir_vec_t *sol_vel_press,
                         subd_options_t *subd_options)
{
  const char      *this_fn_name = "subd_vtk_write_test";

  if (vtk_write_test_path != NULL)  {
    if(subd_options->test_options->test_manufactured != SUBD_TEST_MANUFACTURED_NONE) {
      ymir_vec_t         *vel_ref = rhea_velocity_new (ymir_mesh);
      ymir_vec_t         *vel_chk = rhea_velocity_new (ymir_mesh);
      ymir_vec_t         *pres_ref = ymir_cvec_new (ymir_mesh, 1);
      ymir_vec_t         *pres_chk = rhea_pressure_new (ymir_mesh, press_elem);
      double              vel_abs_err, vel_rel_err;
      ymir_vec_t         *vel_err = rhea_velocity_new (ymir_mesh);
      ymir_vec_t         *pres_err = ymir_cvec_new (ymir_mesh, 1);
      const               subd_test_manufactured_t
                          test_type = subd_options->test_options->test_manufactured;
      ymir_stokes_op_t   *stokes_op;
      ymir_stress_op_t   *stress_op;

      stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
      stress_op = stokes_op->stress_op;
      ymir_vec_set_value (pres_ref, .0);

      /* compute velocity fields */
      switch (test_type) {
        case SUBD_TEST_MANUFACTURED_SINCOS1_ISO:
        case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT90:
        case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT45:
        case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT60:
        case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT60_VISCEXP60:
          /* compute reference velocity field (output) */
          ymir_cvec_set_function (vel_ref, subd_test_sincos1_vel_in_fn,
                                  NULL);
          break;

        case SUBD_TEST_MANUFACTURED_POLY1_TIROT90:
        case SUBD_TEST_MANUFACTURED_POLY1_TIROT90_VISCEXP:
          /* compute reference velocity field (output) */
          ymir_cvec_set_function (vel_ref, subd_test_poly1_vel_in_fn,
                                  NULL);

          break;

       default:
          RHEA_ABORT_NOT_REACHED ();
      }

      ymir_stokes_vec_get_components (sol_vel_press, vel_chk, pres_chk,
                                        press_elem);
      subd_test_manufactured_compute_vel_err (&vel_abs_err, &vel_rel_err,
                                                 vel_err, vel_ref, vel_chk,
                                                 stress_op);

      RHEA_GLOBAL_INFOF ("manufactured solution test (test type %i): abs error %1.3e rel error %1.3e\n",
                          test_type, vel_abs_err, vel_rel_err);


      char      path[BUFSIZ];
      snprintf (path, BUFSIZ, "%s_manufactured", vtk_write_test_path);
      ymir_vtk_write (ymir_mesh, path,
                      vel_ref, "vel_reference",
                      vel_chk, "vel_check",
                      vel_err, "vel_error",
                      pres_ref, "pressure_reference",
                      pres_chk, "pressure_check",
                      NULL);
      rhea_velocity_destroy (vel_ref);
      rhea_velocity_destroy (vel_chk);
      rhea_velocity_destroy (vel_err);
      rhea_pressure_destroy (pres_chk);
      ymir_vec_destroy (pres_ref);
      ymir_vec_destroy (pres_err);
    }

    if (subd_options->test_options->test_stress_comp != SUBD_TEST_MANUFACTURED_NONE)
    {
      RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

      const               subd_test_manufactured_t
                          test_type = subd_options->test_options->test_stress_comp;
      ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (ymir_mesh->ma->N,
                                                                ymir_mesh->ma->ompsize);
      ymir_vec_t         *velocity = rhea_velocity_new (ymir_mesh);
      ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);
      ymir_vec_t         *shear_visc = rhea_viscosity_new (ymir_mesh);
      ymir_vec_t         *TI_tensor = ymir_dvec_new (ymir_mesh, 9, YMIR_GAUSS_NODE);
      ymir_vec_t         *trac_n = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
      ymir_vec_t         *trac_s = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
      ymir_vec_t         *normal_force = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
      ymir_vec_t         *shear_force = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
      ymir_vec_t         *traction_ref = rhea_velocity_new (ymir_mesh);
      ymir_vec_t         *stress = ymir_dvec_new (ymir_mesh, 3, YMIR_GAUSS_NODE);
      ymir_vec_t         *stress_ref = rhea_velocity_new (ymir_mesh);

      ymir_stokes_vec_get_velocity (sol_vel_press, velocity,
                                      press_elem);
      rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

      if (subd_options->visc_options->anisotropy_type == SUBD_VISC_TRANSVERSELY_ISOTROPY)  {
        ymir_stokes_op_t      *stokes_op;
        ymir_stress_op_t      *stress_op;

        /* get the viscous stress operator */
        stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
        stress_op = stokes_op->stress_op;

        /* copy shear viscosity */
        subd_stress_op_copy_shear_visc (shear_visc, stress_op);
        subd_stress_op_copy_TI_tensor (TI_tensor, stress_op);

        subd_manufactured_strainvec_coupling_compute (velocity, trac_n, trac_s,
                                               shear_visc, viscosity);

        subd_manufactured_stressvec_coupling_compute (velocity, viscosity,
                                              shear_visc, TI_tensor,
                                              normal_force, shear_force,
                                              vel_elem, subd_options);

        subd_stressvec_TI (velocity, stress, viscosity, shear_visc, TI_tensor, vel_elem);
        switch (test_type) {
          case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT90:
            ymir_cvec_set_function (traction_ref, subd_test_sincos1_TIrot90_traction_fn,
                                      NULL);
            ymir_cvec_set_function (stress_ref, subd_test_sincos1_TIrot90_stress_fn,
                                      NULL);
          case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT45:
            ymir_cvec_set_function (traction_ref, subd_test_sincos1_TIrot45_traction_fn,
                                      NULL);
            ymir_cvec_set_function (stress_ref, subd_test_sincos1_TIrot45_stress_fn,
                                      NULL);
          case SUBD_TEST_MANUFACTURED_SINCOS1_TIROT60:
            ymir_cvec_set_function (traction_ref, subd_test_sincos1_TIrot60_traction_fn,
                                      NULL);
            ymir_cvec_set_function (stress_ref, subd_test_sincos1_TIrot60_stress_fn,
                                      NULL);
            break;

          default:
            RHEA_ABORT_NOT_REACHED ();
        }
      }
      char      path[BUFSIZ];
      snprintf (path, BUFSIZ, "%s_stress_component", vtk_write_test_path);
      ymir_vtk_write (ymir_mesh, path,
                      stress_ref, "stress_reference (yy,zz,yz)",
                      stress, "stress_check",
                      traction_ref, "traction_reference (0,fn,fs)",
                      normal_force, "normal_force (from stress)",
                      shear_force, "shear_force (from stress)",
                      trac_n, "normal_force (from strain)",
                      trac_s, "shear_force (from strain)",
                      NULL);

      ymir_velocity_elem_destroy (vel_elem);
      rhea_velocity_destroy (velocity);
      rhea_viscosity_destroy (viscosity);
      rhea_viscosity_destroy (shear_visc);
      ymir_vec_destroy (TI_tensor);
      ymir_vec_destroy (trac_n);
      ymir_vec_destroy (trac_s);
      ymir_vec_destroy (normal_force);
      ymir_vec_destroy (shear_force);
      ymir_vec_destroy (traction_ref);
      ymir_vec_destroy (stress);
      ymir_vec_destroy (stress_ref);
    } /* end of stress component*/
    /*stress operator test, TODO*/
  }
}

void
subd_vtk_write_postp (ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         rhea_stokes_problem_t *stokes_problem,
                         ymir_vec_t *sol_vel_press,
                         subd_options_t *subd_options)
{
  const char *this_fn_name = "subd_vtk_write_postp";

     /* compute and output analysis of stress */
  if (vtk_write_postp_path != NULL)  {
    RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

    ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (ymir_mesh->ma->N,
                                                              ymir_mesh->ma->ompsize);
    ymir_vec_t         *velocity = rhea_velocity_new (ymir_mesh);
    ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *shear_visc = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *TI_tensor = ymir_dvec_new (ymir_mesh, 9, YMIR_GAUSS_NODE);
    ymir_vec_t         *traction = ymir_dvec_new (ymir_mesh, 3, YMIR_GAUSS_NODE);
    ymir_vec_t         *normal_force = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
    ymir_vec_t         *shear_force = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);

    ymir_stokes_vec_get_velocity (sol_vel_press, velocity,
                                    press_elem);
    rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

    if (subd_options->visc_options->anisotropy_type == SUBD_VISC_TRANSVERSELY_ISOTROPY)  {
      ymir_stokes_op_t      *stokes_op;
      ymir_stress_op_t      *stress_op;

      /* get the viscous stress operator */
      stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
      stress_op = stokes_op->stress_op;

      /* copy shear viscosity */
      subd_stress_op_copy_shear_visc (shear_visc, stress_op);
      subd_stress_op_copy_TI_tensor (TI_tensor, stress_op);

      subd_postp_weakzone_traction_compute (velocity, traction,
                                             shear_visc, viscosity,
                                             subd_options);

      subd_postp_weakzone_coupling_compute (velocity, viscosity,
                                            shear_visc, TI_tensor,
                                            normal_force, shear_force,
                                            vel_elem, subd_options);
    }
    else  {
      subd_postp_weakzone_traction_compute (velocity, traction,
                                             viscosity, viscosity,
                                             subd_options);

      subd_postp_weakzone_coupling_compute (velocity, viscosity,
                                             shear_visc, TI_tensor,
                                             normal_force, shear_force,
                                             vel_elem, subd_options);
   }
   rhea_viscosity_destroy (shear_visc);
   ymir_vec_destroy (TI_tensor);

   {
     char            path[BUFSIZ];

     snprintf (path, BUFSIZ, "%s", vtk_write_postp_path);
     ymir_vtk_write (ymir_mesh, path,
                     normal_force, "normal_force from sigma",
                     shear_force, "shear_force from sigma",
                     traction, "trac_n from edot",
                     NULL);
   }

    /* destroy */
    rhea_viscosity_destroy (viscosity);
    rhea_velocity_destroy (velocity);
    ymir_vec_destroy (normal_force);
    ymir_vec_destroy (shear_force);
    ymir_vec_destroy (traction);
    ymir_velocity_elem_destroy (vel_elem);

    RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
  }
}

void
subd_vtk_write_stress (ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         rhea_stokes_problem_t *stokes_problem,
                         ymir_vec_t *sol_vel_press,
                         subd_options_t *subd_options)
{
  const char        *this_fn_name = "subd_vtk_write_stress";
  /* compute and output second invariant strain_rate, stress, and surface normal stress  */
  if (vtk_write_stress_path != NULL)  {
    RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

    ymir_vec_t         *velocity = rhea_velocity_new (ymir_mesh);
    ymir_vec_t         *viscosity = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *n_stress = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *s_stress = rhea_viscosity_new (ymir_mesh);
    ymir_vec_t         *surf_ntau = ymir_face_dvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1,
                                                     YMIR_GAUSS_NODE);
    ymir_vec_t         *surf_stau = ymir_face_dvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1,
                                                     YMIR_GAUSS_NODE);
    ymir_vec_t         *traction = ymir_dvec_new (ymir_mesh, 3, YMIR_GAUSS_NODE);
    double             n_dir[3] = {.0, .0, 1.0};

    ymir_velocity_elem_t  *vel_elem = ymir_velocity_elem_new (ymir_mesh->ma->N,
                                                              ymir_mesh->ma->ompsize);
    ymir_vec_t            *edotII = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
    ymir_vec_t            *tauII = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
    ymir_vec_t            *topography = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);

    ymir_stokes_vec_get_velocity (sol_vel_press, velocity,
                                    press_elem);
    rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);

    /* compute 2nd invariant of the strain rate */
    ymir_second_invariant_vec (velocity, edotII, vel_elem);
    ymir_vec_sqrt (edotII, edotII);

    /* compute 2nd invariant of deviatoric stress tau = 2* (2nd invariant of strain_rate * viscosity )
      and its projection on the surface */
    if (subd_options->visc_options->anisotropy_type == SUBD_VISC_TRANSVERSELY_ISOTROPY)  {
      ymir_stokes_op_t      *stokes_op;
      ymir_stress_op_t      *stress_op;
      ymir_vec_t            *shear_visc = rhea_viscosity_new (ymir_mesh);
      ymir_vec_t            *TI_tensor = ymir_dvec_new (ymir_mesh, 9,
                                                        YMIR_GAUSS_NODE);

      /* get the viscous stress operator */
      stokes_op = rhea_stokes_problem_get_stokes_op (stokes_problem);
      stress_op = stokes_op->stress_op;

      /* copy shear viscosity */
      subd_stress_op_copy_shear_visc (shear_visc, stress_op);
      subd_stress_op_copy_TI_tensor (TI_tensor, stress_op);

      subd_2inv_stress_TI (velocity, tauII,
                              viscosity, shear_visc, TI_tensor, vel_elem);

      rhea_viscosity_destroy (shear_visc);
      ymir_vec_destroy (TI_tensor);
    }
    else  {
      ymir_vec_copy (edotII, tauII)
      ymir_vec_multiply_in (viscosity, tauII);
      ymir_vec_scale (2.0, tauII);
    }

    subd_postp_topography (topography, stokes_problem, subd_options);

    subd_normal_stress (velocity, n_stress, s_stress, traction, n_dir, viscosity, vel_elem);
    ymir_interp_vec (n_stress, surf_ntau);
    ymir_interp_vec (s_stress, surf_stau);
    {
      char            path[BUFSIZ];

      snprintf (path, BUFSIZ, "%s", vtk_write_stress_path);
      ymir_vtk_write (ymir_mesh, path,
                      edotII, "edotII",
                      tauII, "tauII",
                      n_stress, "normal_stress",
                      s_stress, "shear_stress",
                      surf_ntau, "surf_normal_tau",
                      surf_stau, "surf_shear_tau",
                      topography, "topography",
                      NULL);
    }

    /* destroy */
    rhea_viscosity_destroy (viscosity);
    rhea_velocity_destroy (velocity);
    ymir_vec_destroy (n_stress);
    ymir_vec_destroy (s_stress);
    ymir_vec_destroy (surf_ntau);
    ymir_vec_destroy (surf_stau);
    ymir_vec_destroy (traction);
    ymir_vec_destroy (edotII);
    ymir_vec_destroy (tauII);
    ymir_vec_destroy (topography);
    ymir_velocity_elem_destroy (vel_elem);

    RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
  }
}

void
subd_vtk_write_freesurface (ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         ymir_vec_t *sol_vel_press)
{
  const char      *this_fn_name = "subd_vtk_write_freesurface";

 /* compute and output second invariant strain_rate, stress, and surface normal stress  */
  if (vtk_write_freesurface_path != NULL)  {
    RHEA_GLOBAL_PRODUCTIONF ("In %s: Start vtk_write_freesurface\n", this_fn_name);

    ymir_vec_t         *velocity = rhea_velocity_new (ymir_mesh);
    ymir_vec_t         *surf_normal_velocity = ymir_face_cvec_new (ymir_mesh,
                                                     RHEA_DOMAIN_BOUNDARY_FACE_TOP, 3);

    ymir_stokes_vec_get_velocity (sol_vel_press, velocity,
                                    press_elem);

    ymir_interp_vec (velocity, surf_normal_velocity);

    {
      char            path[BUFSIZ];

      snprintf (path, BUFSIZ, "%s", vtk_write_freesurface_path);
      ymir_vtk_write (ymir_mesh, path,
                      surf_normal_velocity, "surf_normal_velocity",
                      NULL);
    }

    /* destroy */
    rhea_velocity_destroy (velocity);
    ymir_vec_destroy (surf_normal_velocity);
    RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
  }
}

void
subd_vtk_write_neumann (ymir_mesh_t *ymir_mesh,
                         ymir_pressure_elem_t *press_elem,
                         rhea_stokes_problem_t *stokes_problem,
                         ymir_vec_t *sol_vel_press)
{
  const char        *this_fn_name = "subd_vtk_write_neumann";
  if (vtk_write_neumann_path != NULL)  {
    RHEA_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

    ymir_vec_t         *surf_neumann;
    ymir_vec_t      *velocity = rhea_velocity_new (ymir_mesh);
    ymir_vec_t      *surf_velocity = ymir_face_cvec_new (ymir_mesh,
                                                             RHEA_DOMAIN_BOUNDARY_FACE_TOP, 3);

    ymir_stokes_vec_get_velocity (sol_vel_press, velocity,
                                            press_elem);
    ymir_interp_vec (velocity, surf_velocity);

    surf_neumann = rhea_stokes_problem_get_rhs_vel_nonzero_neumann_surface (stokes_problem);

    {
      char            path[BUFSIZ];

      snprintf (path, BUFSIZ, "%s", vtk_write_neumann_path);
      if (surf_neumann == NULL) {
        ymir_vtk_write (ymir_mesh, path,
                        surf_velocity, "surf_velocity",
                        NULL);
      }
      else  {
        ymir_vtk_write (ymir_mesh, path,
                        surf_velocity, "surf_velocity",
                        surf_neumann, "surf_neumann",
                        NULL);
      }
    }

    rhea_velocity_destroy (velocity);
    ymir_vec_destroy (surf_velocity);
    /* destroy */
    RHEA_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
  }
}
