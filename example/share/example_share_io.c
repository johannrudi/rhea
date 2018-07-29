#include <example_share_io.h>
#include <rhea_base.h>
#include <rhea_velocity.h>
#include <rhea_velocity_pressure.h>
#include <rhea_stress.h>
#include <rhea_io_std.h>

void
example_share_io_write_solution_surf_txt (
                                      const char *file_path_txt,
                                      rhea_domain_coordinate_type_t coord_type,
                                      ymir_vec_t *sol_vel_press,
                                      rhea_stokes_problem_t *stokes_problem)
{
  const ymir_topidx_t meshid = RHEA_DOMAIN_BOUNDARY_FACE_TOP;
  const int           n_fields_out = 3 + 3 + 1;
  rhea_domain_options_t      *domain_options;
  rhea_temperature_options_t *temp_options;
  rhea_viscosity_options_t   *visc_options;
  ymir_mesh_t          *ymir_mesh;
  ymir_pressure_elem_t *press_elem;
  ymir_vec_t         *velocity;
  ymir_vec_t         *velocity_surf, *stress_norm_surf;
  ymir_vec_t         *out_surf;
  double             *out_data;

  /* exit if nothing to do */
  if (file_path_txt == NULL) {
    return;
  }

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (__func__, "path=\"%s\"", file_path_txt);

  /* get options */
  domain_options = rhea_stokes_problem_get_domain_options (stokes_problem);
  temp_options = rhea_stokes_problem_get_temperature_options (stokes_problem);
  visc_options = rhea_stokes_problem_get_viscosity_options (stokes_problem);

  /* get mesh */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);

  /* create & get velocity */
  velocity = rhea_velocity_new (ymir_mesh);
  rhea_velocity_pressure_copy_components (velocity, NULL, sol_vel_press,
                                          press_elem);

  /* create surface variables */
  velocity_surf = rhea_velocity_surface_new (ymir_mesh);
  stress_norm_surf = rhea_stress_surface_new (ymir_mesh);

  /* get surface fields */
  rhea_velocity_surface_interpolate (velocity_surf, velocity);
  rhea_stokes_problem_stress_compute_normal_at_surface (stress_norm_surf,
                                                        sol_vel_press,
                                                        stokes_problem);

  /* convert to physical dimensions */
  rhea_velocity_convert_to_dimensional_mm_yr (velocity_surf, domain_options,
                                              temp_options);
  rhea_stress_convert_to_dimensional_Pa (stress_norm_surf, domain_options,
                                         temp_options, visc_options);

  /* create & set output vector */
  out_surf = ymir_face_cvec_new (ymir_mesh, meshid, n_fields_out);
  {
    ymir_vec_t         *x_vec, *y_vec, *z_vec;

    /* get cartesian coordinates */
    x_vec = ymir_face_cvec_new (ymir_mesh, meshid, 1);
    y_vec = ymir_face_cvec_new (ymir_mesh, meshid, 1);
    z_vec = ymir_face_cvec_new (ymir_mesh, meshid, 1);
    ymir_vec_get_coords (x_vec, y_vec, z_vec);

    /* set output entries */
    if (YMIR_CVEC_STRIDE == YMIR_STRIDE_NODE) {
      const ymir_locidx_t n_nodes = ymir_mesh->fmeshes[meshid].Ncn;
      ymir_locidx_t       nodeid;
      const double       *_sc_restrict x = ymir_cvec_index (x_vec, 0, 0);
      const double       *_sc_restrict y = ymir_cvec_index (y_vec, 0, 0);
      const double       *_sc_restrict z = ymir_cvec_index (z_vec, 0, 0);
      const double       *_sc_restrict v = ymir_cvec_index (velocity_surf,
                                                            0, 0);
      const double       *_sc_restrict s = ymir_cvec_index (stress_norm_surf,
                                                            0, 0);
      double             *_sc_restrict out = ymir_cvec_index (out_surf, 0, 0);

      for (nodeid = 0; nodeid < n_nodes; nodeid++) {
        /* set coordinates */
        rhea_domain_convert_coordinates (&out[n_fields_out*nodeid    ],
                                         &out[n_fields_out*nodeid + 1],
                                         &out[n_fields_out*nodeid + 2],
                                         x[nodeid], y[nodeid], z[nodeid],
                                         coord_type, domain_options);

        /* set velocity */
        out[n_fields_out*nodeid + 3] = v[3*nodeid    ];
        out[n_fields_out*nodeid + 4] = v[3*nodeid + 1];
        out[n_fields_out*nodeid + 5] = v[3*nodeid + 2];

        /* set normal stress */
        out[n_fields_out*nodeid + 6] = s[nodeid];
      }
    }
    else {
      RHEA_ABORT_NOT_REACHED ();
    }

    /* destroy */
    ymir_vec_destroy (x_vec);
    ymir_vec_destroy (y_vec);
    ymir_vec_destroy (z_vec);
  }

  /* write txt */
  {
    ymir_face_mesh_t   *face_mesh = &(ymir_mesh->fmeshes[meshid]);
    const ymir_locidx_t n_nodes = face_mesh->Ncn;
    const int           mpisize = ymir_mesh_get_MPI_Comm_size (ymir_mesh);
    const int           mpirank = ymir_mesh_get_MPI_Comm_rank (ymir_mesh);
    const int           n_digits =
                          1 + (int) floor (log ((double) mpisize) / M_LN10);
    char                path[BUFSIZ];

    /* get pointer to data */
    if (0 < n_nodes) { /* if nodes exist on this rank */
      out_data = ymir_cvec_index (out_surf, 0, 0);
    }
    else { /* otherwise this rank is empty */
      out_data = NULL;
    }

    /* write */
    snprintf (path, BUFSIZ, "%s_%0*i.txt", file_path_txt, n_digits, mpirank);
    rhea_io_std_write_double_to_txt (path, out_data, n_fields_out*n_nodes,
                                     n_fields_out);
  }

  /* destroy */
  ymir_vec_destroy (out_surf);
  rhea_velocity_destroy (velocity);
  rhea_velocity_surface_destroy (velocity_surf);
  rhea_stress_surface_destroy (stress_norm_surf);

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}
