#include <rhea_inversion_obs_viscosity.h>
#include <rhea_base.h>
#include <ymir_mass_vec.h>

rhea_domain_subset_column_t **
rhea_inversion_obs_viscosity_new (int *n_columns,
                                  double **value,
                                  double **weight,
                                  const rhea_inversion_obs_viscosity_t obs_type,
                                  char *value_Pas_list,
                                  char *stddev_list,
                                  rhea_plate_options_t *plate_options,
                                  rhea_viscosity_options_t *visc_options)
{
  const double        visc_dim_Pas =
                        rhea_viscosity_get_dim_Pas (visc_options);
  int                 n_entries, cid;
  rhea_domain_subset_column_t **column;

  /* check input */
  RHEA_ASSERT (NULL != visc_options);

  /* create columns */
  switch (obs_type) {
  case RHEA_INVERSION_OBS_VISCOSITY_NONE:
    *n_columns = 0;
    column = NULL;
    break;
  case RHEA_INVERSION_OBS_VISCOSITY_AVERAGE_REGION:
    RHEA_ABORT_NOT_REACHED (); //TODO
    break;
  case RHEA_INVERSION_OBS_VISCOSITY_AVERAGE_UNDER_PLATES:
    RHEA_ASSERT (NULL != plate_options);
    {
      rhea_domain_options_t  *domain_options = plate_options->domain_options;
      const double        rmin = domain_options->lm_um_interface_radius;
      const double        rmax = domain_options->radius_max;
      int                 n_vert;

      /* set number of columns (=plates) */
      *n_columns = rhea_plate_get_n_plates (plate_options);
      if (*n_columns <= 0) {
        column = NULL;
        break;
      }

      /* create columns */
      column = RHEA_ALLOC (rhea_domain_subset_column_t *, *n_columns);
      for (cid = 0; cid < *n_columns; cid++) {
        column[cid] = RHEA_ALLOC (rhea_domain_subset_column_t, 1);

        if (0 < plate_options->n_polygons &&
            0 < plate_options->n_vertices[cid]) {
          RHEA_ASSERT (*n_columns == plate_options->n_polygons);
          column[cid]->polygon_n_vertices = n_vert
                                          = plate_options->n_vertices[cid];
          column[cid]->polygon_vertices_x = RHEA_ALLOC (float, n_vert);
          column[cid]->polygon_vertices_y = RHEA_ALLOC (float, n_vert);
          memcpy (column[cid]->polygon_vertices_x,
                  plate_options->vertices_x[cid], n_vert * sizeof (float));
          memcpy (column[cid]->polygon_vertices_y,
                  plate_options->vertices_y[cid], n_vert * sizeof (float));
          column[cid]->polygon_coord_type =
            RHEA_DOMAIN_COORDINATE_SPHERICAL_GEO_DIM;
        }

        if (0 < plate_options->xsection_n_intervals) {
          RHEA_ASSERT (*n_columns == plate_options->xsection_n_intervals);
          column[cid]->xsection_boundary[0] =
            plate_options->xsection_boundary[2*cid  ];
          column[cid]->xsection_boundary[1] =
            plate_options->xsection_boundary[2*cid+1];
        }

        column[cid]->radius_min = rmin;
        column[cid]->radius_max = 0.25*rmin + 0.75*rmax; /* min+3/4*(max-min) */

        column[cid]->volume = NAN; /* will be computed later */
      }
    }
    break;
  default: /* unknown type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* create observatianal values and weights */
  if (0 < *n_columns) {
    /* get observational values and standard deviations */
    n_entries = ymir_options_convert_string_to_double (value_Pas_list, value);
    RHEA_ASSERT (n_entries == *n_columns);
    n_entries = ymir_options_convert_string_to_double (stddev_list, weight);
    RHEA_ASSERT (n_entries == *n_columns);

    /* convert entries */
    for (cid = 0; cid < *n_columns; cid++) {
      (*value)[cid] *= 1.0 / visc_dim_Pas;
      (*weight)[cid] = 1.0 / ((*weight)[cid] * (*weight)[cid]);
    }
  }
  else {
    *value  = NULL;
    *weight = NULL;
  }

  return column;
}

void
rhea_inversion_obs_viscosity_destroy (rhea_domain_subset_column_t **column,
                                      const int n_columns,
                                      double *value,
                                      double *weight)
{
  int                 cid;

  if (0 < n_columns) {
    for (cid = 0; cid < n_columns; cid++) {
      RHEA_FREE (column[cid]);
    }
    RHEA_FREE (column);
    YMIR_FREE (value);
    YMIR_FREE (weight);
  }
}

static double
_get_volume (rhea_domain_subset_column_t *column,
             ymir_mesh_t *ymir_mesh,
             rhea_domain_options_t *domain_options)
{
  ymir_vec_t         *unit, *unit_mass;

  /* return volume if known */
  if (isfinite (column->volume)) {
    return column->volume;
  }

  /* compute volume */
  RHEA_ASSERT (NULL != ymir_mesh);
  RHEA_ASSERT (NULL != domain_options);

  unit      = rhea_viscosity_new (ymir_mesh);
  unit_mass = rhea_viscosity_new (ymir_mesh);

  ymir_vec_set_value (unit, 1.0);
  ymir_mass_apply (unit, unit_mass);
  rhea_domain_subset_apply_filter (unit, NULL, domain_options, column);
  column->volume = ymir_vec_innerprod (unit, unit_mass);

  rhea_viscosity_destroy (unit);
  rhea_viscosity_destroy (unit_mass);

  RHEA_ASSERT (isfinite (column->volume));
  RHEA_ASSERT (0.0 < column->volume);
  return column->volume;
}

static void
_log_fn (double *v, double x, double y, double z, ymir_locidx_t nodeid,
         void *data)
{
  *v = log (*v);
}

static double
rhea_inversion_obs_viscosity_misfit_log_avg_diff (
                                ymir_vec_t *forward_vel,
                                const rhea_inversion_obs_viscosity_t obs_type,
                                rhea_domain_subset_column_t *column,
                                const double obs_val,
                                rhea_stokes_problem_t *stokes_problem,
                                ymir_vec_t *viscosity)
{
  rhea_viscosity_options_t *visc_options =
    rhea_stokes_problem_get_viscosity_options (stokes_problem);
  rhea_domain_options_t *domain_options =
    rhea_stokes_problem_get_domain_options (stokes_problem);
  ymir_mesh_t        *ymir_mesh =
                        rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  ymir_vec_t         *temperature =
                        rhea_stokes_problem_get_temperature (stokes_problem);
  ymir_vec_t         *weakzone =
                        rhea_stokes_problem_get_weakzone (stokes_problem);
  ymir_vec_t         *visc, *marker, *unit_mass;
  const double        volume = _get_volume (column, ymir_mesh, domain_options);
  double              log_visc_int;

  /* create work vectors */
  visc      = rhea_viscosity_new (ymir_mesh);
  marker    = rhea_viscosity_new (ymir_mesh);
  unit_mass = rhea_viscosity_new (ymir_mesh);

  /* get viscosity */
  if (NULL != viscosity) {
    ymir_vec_copy (viscosity, visc);
  }
  else {
    rhea_viscosity_compute (
        /* out: */ visc, NULL, marker,
        /* in:  */ temperature, weakzone, forward_vel, visc_options);
  }

  /* compute log of viscosity */
  ymir_dvec_set_function (visc, _log_fn, NULL);

  /* filter viscosity in a region */
  rhea_domain_subset_apply_filter (visc, NULL, domain_options, column);

  /* compute integral of log viscosity */
  ymir_vec_set_value (marker, 1.0);
  ymir_mass_apply (marker, unit_mass);
  log_visc_int = ymir_vec_innerprod (visc, unit_mass);

  /* destroy */
  rhea_viscosity_destroy (visc);
  rhea_viscosity_destroy (marker);
  rhea_viscosity_destroy (unit_mass);

  /* return distance of viscosities: log(visc_avg) - log(\int visc/volume) */
  RHEA_ASSERT (isfinite (log_visc_int));
  RHEA_ASSERT (isfinite (volume));
  return log (obs_val) - log_visc_int/volume;
}

double
rhea_inversion_obs_viscosity_misfit (
                                ymir_vec_t *forward_vel,
                                const rhea_inversion_obs_viscosity_t obs_type,
                                const int n_columns,
                                rhea_domain_subset_column_t **column,
                                const double *obs_val,
                                const double *weight,
                                rhea_stokes_problem_t *stokes_problem)
{
  rhea_viscosity_options_t *visc_options =
    rhea_stokes_problem_get_viscosity_options (stokes_problem);
  ymir_mesh_t        *ymir_mesh =
                        rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  ymir_vec_t         *temperature =
                        rhea_stokes_problem_get_temperature (stokes_problem);
  ymir_vec_t         *weakzone =
                        rhea_stokes_problem_get_weakzone (stokes_problem);
  ymir_vec_t         *viscosity, *marker;
  double              misfit, diff;
  int                 cid;

  /* return zero if nothing to do */
  if (RHEA_INVERSION_OBS_VISCOSITY_NONE == obs_type) {
    return 0.0;
  }

  /* compute viscosity */
  viscosity = rhea_viscosity_new (ymir_mesh);
  marker    = rhea_viscosity_new (ymir_mesh);
  rhea_viscosity_compute (
      /* out: */ viscosity, NULL, marker,
      /* in:  */ temperature, weakzone, forward_vel, visc_options);
  rhea_viscosity_destroy (marker);

  /* compute sum of weighted data misfits */
  misfit = 0.0;
  for (cid = 0; cid < n_columns; cid++) {
    /* compute distance between average model viscosity and data */
    diff = rhea_inversion_obs_viscosity_misfit_log_avg_diff (
        forward_vel, obs_type, column[cid], obs_val[cid], stokes_problem,
        viscosity);
    misfit += diff * weight[cid] * diff;
  }
  rhea_viscosity_destroy (viscosity);

  /* return misfit term of objective functional */
  return misfit;
}

void
rhea_inversion_obs_viscosity_adjoint_rhs (
                                ymir_vec_t *rhs_vel_mass,
                                ymir_vec_t *forward_vel,
                                const rhea_inversion_obs_viscosity_t obs_type)
{

  /* return zero if nothing to do */
  if (RHEA_INVERSION_OBS_VISCOSITY_NONE == obs_type) {
    ymir_vec_set_zero (rhs_vel_mass);
    return;
  }

}

void
rhea_inversion_obs_viscosity_incremental_adjoint_rhs (
                                ymir_vec_t *rhs_vel_mass,
                                ymir_vec_t *forward_vel,
                                const rhea_inversion_obs_viscosity_t obs_type)
{

  /* return zero if nothing to do */
  if (RHEA_INVERSION_OBS_VISCOSITY_NONE == obs_type) {
    ymir_vec_set_zero (rhs_vel_mass);
    return;
  }

}
