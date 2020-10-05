#include <rhea_inversion_obs_viscosity.h>
#include <rhea_strainrate.h>
#include <rhea_base.h>
#include <ymir_mass_vec.h>
#include <ymir_stress_pc.h>

rhea_domain_subset_column_t **
rhea_inversion_obs_viscosity_new (int *n_columns,
                                  double **value,
                                  double **weight,
                                  const rhea_inversion_obs_viscosity_t obs_type,
                                  char *value_Pas_list,
                                  char *stddev_rel_list,
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

      /* set number of columns (= #plates) */
      *n_columns = rhea_plate_get_n_plates (plate_options);
      if (*n_columns <= 0) {
        column = NULL;
        break;
      }

      /* create columns */
      column = RHEA_ALLOC (rhea_domain_subset_column_t *, *n_columns);
      for (cid = 0; cid < *n_columns; cid++) {
        column[cid] = RHEA_ALLOC (rhea_domain_subset_column_t, 1);

        /* copy polygon vertices if they exist */
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
        else {
          column[cid]->polygon_n_vertices = 0;
          column[cid]->polygon_vertices_x = NULL;
          column[cid]->polygon_vertices_y = NULL;
        }

        /* copy interval boundaries if they exist */
        if (0 < plate_options->xsection_n_intervals) {
          RHEA_ASSERT (*n_columns == plate_options->xsection_n_intervals);
          column[cid]->xsection_boundary[0] =
            plate_options->xsection_boundary[2*cid  ];
          column[cid]->xsection_boundary[1] =
            plate_options->xsection_boundary[2*cid+1];
        }
        else {
          column[cid]->xsection_boundary[0] = NAN;
          column[cid]->xsection_boundary[1] = NAN;
        }

        /* set radii of column */
        column[cid]->radius_min = rmin;
        column[cid]->radius_max = 0.25*rmin + 0.75*rmax; /* min+3/4*(max-min) */

        /* init volume of column */
        column[cid]->volume = NAN; /* will be computed later */
      }
    }
    break;
  default: /* unknown type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* create observatianal values and weights */
  if (0 < *n_columns) {
    double             *value_Pas  = NULL;
    double             *stddev_rel = NULL;

    /* get observational values and standard deviations */
    n_entries = ymir_options_convert_string_to_double (value_Pas_list,
                                                       &value_Pas);
    RHEA_ASSERT (n_entries == *n_columns);
    n_entries = ymir_options_convert_string_to_double (stddev_rel_list,
                                                       &stddev_rel);
    RHEA_ASSERT (n_entries == *n_columns);

    /* convert entries */
    *value  = RHEA_ALLOC (double, *n_columns);
    *weight = RHEA_ALLOC (double, *n_columns);
    for (cid = 0; cid < *n_columns; cid++) {
      const double        val = value_Pas[cid] / visc_dim_Pas;
      const double        std = log (stddev_rel[cid]);

      (*value)[cid]  = val;
      (*weight)[cid] = 1.0 / std;
    }
    YMIR_FREE (value_Pas);  /* was allocated in ymir */
    YMIR_FREE (stddev_rel); /* was allocated in ymir */
  }
  else {
    *value  = NULL;
    *weight = NULL;
  }

  /* print values and weights */
  if (0 < *n_columns) {
    RHEA_GLOBAL_INFO ("========================================\n");
    RHEA_GLOBAL_INFOF ("%s: value, value [Pas], weight, stddev\n", __func__);
    RHEA_GLOBAL_INFO ("----------------------------------------\n");
    for (cid = 0; cid < *n_columns; cid++) {
      RHEA_GLOBAL_INFOF ("column_idx %3i: %.3e, %.3e, %.3e, %g\n",
                         cid, (*value)[cid], visc_dim_Pas*(*value)[cid],
                         (*weight)[cid], exp (1.0/(*weight)[cid]));
    }
    RHEA_GLOBAL_INFO ("========================================\n");
  }

  return column;
}

void
rhea_inversion_obs_viscosity_destroy (rhea_domain_subset_column_t **column,
                                      const int n_columns,
                                      double *value,
                                      double *weight)
{
  if (0 < n_columns) {
    int                 cid;

    for (cid = 0; cid < n_columns; cid++) {
      if (NULL != column[cid]->polygon_vertices_x) {
        RHEA_FREE (column[cid]->polygon_vertices_x);
      }
      if (NULL != column[cid]->polygon_vertices_y) {
        RHEA_FREE (column[cid]->polygon_vertices_y);
      }
      RHEA_FREE (column[cid]);
    }
    RHEA_FREE (column);
    RHEA_FREE (value);
    RHEA_FREE (weight);
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

//TODO implement without callback fn. in ymir_vec_ops
static void
_log_fn (double *v, double x, double y, double z, ymir_locidx_t nodeid,
         void *data)
{
  *v = log (*v);
}

static double
rhea_inversion_obs_viscosity_compute_log_avg (
                                ymir_vec_t *forward_vel,
                                rhea_domain_subset_column_t *column,
                                rhea_stokes_problem_t *stokes_problem,
                                ymir_vec_t *viscosity)
{
  rhea_domain_options_t *domain_options =
    rhea_stokes_problem_get_domain_options (stokes_problem);
  ymir_mesh_t        *ymir_mesh =
                        rhea_stokes_problem_get_ymir_mesh (stokes_problem);
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
    ymir_vec_t         *temperature =
                          rhea_stokes_problem_get_temperature (stokes_problem);
    ymir_vec_t         *weakzone =
                          rhea_stokes_problem_get_weakzone (stokes_problem);

    rhea_stokes_problem_viscosity_compute (
        /* out: */ visc, NULL, marker,
        /* in:  */ temperature, weakzone, forward_vel, stokes_problem);
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
  return log_visc_int/volume;
}

/**
 * Compute the difference between the avarage viscosity and the corresponding
 * value from data:
 *   lob(visc(vel)) - log(visc_obs)
 */
static double
rhea_inversion_obs_viscosity_diff (
                                ymir_vec_t *forward_vel,
                                rhea_domain_subset_column_t *column,
                                const double obs_val,
                                rhea_stokes_problem_t *stokes_problem,
                                ymir_vec_t *viscosity)
{
  double              log_visc_avg;

  log_visc_avg = rhea_inversion_obs_viscosity_compute_log_avg (
      forward_vel, column, stokes_problem, viscosity);
  return log_visc_avg - log (obs_val);
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
  ymir_mesh_t        *ymir_mesh =
                        rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  ymir_vec_t         *temperature =
                        rhea_stokes_problem_get_temperature (stokes_problem);
  ymir_vec_t         *weakzone =
                        rhea_stokes_problem_get_weakzone (stokes_problem);
  ymir_vec_t         *viscosity, *marker;
  double              misfit_norm_sq, diff;
  int                 cid;

  /* return zero if nothing to do */
  if (RHEA_INVERSION_OBS_VISCOSITY_NONE == obs_type) {
    return 0.0;
  }

  /* compute viscosity */
  viscosity = rhea_viscosity_new (ymir_mesh);
  marker    = rhea_viscosity_new (ymir_mesh);
  rhea_stokes_problem_viscosity_compute (
      /* out: */ viscosity, NULL, marker,
      /* in:  */ temperature, weakzone, forward_vel, stokes_problem);
  rhea_viscosity_destroy (marker);

  /* compute sum of weighted data misfits */
  RHEA_GLOBAL_VERBOSE ("========================================\n");
  RHEA_GLOBAL_VERBOSEF ("%s: weight, log(obs_val), diff\n", __func__);
  RHEA_GLOBAL_VERBOSE ("----------------------------------------\n");
  misfit_norm_sq = 0.0;
  for (cid = 0; cid < n_columns; cid++) {
    /* compute distance between average model viscosity and data */
    diff = rhea_inversion_obs_viscosity_diff (
        forward_vel, column[cid], obs_val[cid], stokes_problem, viscosity);
    misfit_norm_sq += weight[cid]*diff * weight[cid]*diff;
    RHEA_GLOBAL_VERBOSEF ("column_idx %3i: %.3e, %.3e, %.3e\n",
                          cid, weight[cid], log (obs_val[cid]), diff);
  }
  rhea_viscosity_destroy (viscosity);
  RHEA_GLOBAL_VERBOSE ("========================================\n");

  /* return misfit term of objective functional */
  return 0.5*misfit_norm_sq;
}

void
rhea_inversion_obs_viscosity_add_adjoint_rhs (
                                ymir_vec_t *rhs_vel_mass,
                                ymir_vec_t *forward_vel,
                                const rhea_inversion_obs_viscosity_t obs_type,
                                const int n_columns,
                                rhea_domain_subset_column_t **column,
                                const double *obs_val,
                                const double *weight,
                                rhea_stokes_problem_t *stokes_problem)
{
  rhea_domain_options_t *domain_options =
    rhea_stokes_problem_get_domain_options (stokes_problem);
  ymir_mesh_t        *ymir_mesh =
                        rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  ymir_vec_t         *temperature =
                        rhea_stokes_problem_get_temperature (stokes_problem);
  ymir_vec_t         *weakzone =
                        rhea_stokes_problem_get_weakzone (stokes_problem);
  ymir_vec_t         *viscosity, *proj_scal, *marker, *coeff;
  ymir_vec_t         *strainrate_2inv, *op_out_vel;
  ymir_vel_dir_t     *vel_dir;
  ymir_stress_op_t   *stress_op;
  double              diff, volume;
  int                 cid;

  /* return if nothing to do */
  if (RHEA_INVERSION_OBS_VISCOSITY_NONE == obs_type) {
    return;
  }

  /* create work vectors */
  viscosity       = rhea_viscosity_new (ymir_mesh);
  proj_scal       = rhea_viscosity_new (ymir_mesh);
  marker = coeff  = rhea_viscosity_new (ymir_mesh);
  strainrate_2inv = rhea_strainrate_2inv_new (ymir_mesh);
  op_out_vel      = rhea_velocity_new (ymir_mesh);

  /* compute viscosity */
  rhea_stokes_problem_viscosity_compute (
      /* out: */ viscosity, proj_scal, marker,
      /* in:  */ temperature, weakzone, forward_vel, stokes_problem);

  /* compute 2nd invariant of the strain rate*/
  rhea_strainrate_compute_sqrt_of_2inv (strainrate_2inv, forward_vel);

  /* create stress operator */
  ymir_vec_set_value (coeff, 1.0);
  vel_dir = rhea_domain_create_velocity_dirichlet_bc (
      ymir_mesh, NULL /* dirscal */, domain_options);
  stress_op = ymir_stress_op_new_ext (
      coeff, vel_dir,
      NULL /* Robin BC's */,
      NULL /* deprecated */,
      NULL /* deprecated */,
      domain_options->center, domain_options->moment_of_inertia);

  /* compute sum of right-hand sides */
  for (cid = 0; cid < n_columns; cid++) {
    /* init coefficient: proj_scal/eII^2 */
    ymir_vec_copy (proj_scal, coeff);
    ymir_vec_divide_in (strainrate_2inv, coeff);
    ymir_vec_divide_in (strainrate_2inv, coeff);

    /* filter coefficient */
    rhea_domain_subset_apply_filter (coeff, NULL, domain_options, column[cid]);

    /* compute distance between average model viscosity and data */
    diff = rhea_inversion_obs_viscosity_diff (
        forward_vel, column[cid], obs_val[cid], stokes_problem, viscosity);

    /* scale coefficient */
    volume  = _get_volume (column[cid], ymir_mesh, domain_options);
    ymir_vec_scale (2.0/2.0 * weight[cid]*diff * weight[cid] / (2.0*volume),
                    coeff);

    /* set the coefficient of the viscous stress operator */
    ymir_stress_op_set_coeff_scal (stress_op, coeff);

    /* apply viscous stress operator to forward velocity */
    ymir_stress_pc_apply_stress_op (forward_vel, op_out_vel, stress_op,
                                    0 /* !linearized */, 0 /* !dirty */);

    /* add output to right-hand side (change sign to obtain RHS) */
    ymir_vec_add (-1.0, op_out_vel, rhs_vel_mass);
  }

  /* destroy */
  ymir_stress_op_destroy (stress_op);
  ymir_vel_dir_destroy (vel_dir);
  rhea_viscosity_destroy (viscosity);
  rhea_viscosity_destroy (proj_scal);
  rhea_viscosity_destroy (coeff);
  rhea_strainrate_2inv_destroy (strainrate_2inv);
  rhea_velocity_destroy (op_out_vel);
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

  //TODO
}
