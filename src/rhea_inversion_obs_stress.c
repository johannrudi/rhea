#include <rhea_inversion_obs_stress.h>
#include <rhea_base.h>
#include <rhea_stress.h>
#include <rhea_velocity.h>
#include <ymir_mass_vec.h>
#include <ymir_stokes_op.h>
#include <ymir_stress_op_optimized.h>

static int
rhea_inversion_obs_stress_weight_check_vec_type (ymir_vec_t *vec)
{
  if (1 == vec->ndfields) {
    return (
        ymir_vec_is_dvec (vec) &&
        vec->node_type == YMIR_GAUSS_NODE
    );
  }
  else {
    return rhea_stress_check_vec_type (vec);
  }
}

static int
rhea_inversion_obs_stress_weight_is_valid (ymir_vec_t *vec)
{
  return rhea_stress_is_valid (vec);
}

void
rhea_inversion_obs_stress_diff (ymir_vec_t *misfit_stress,
                                ymir_vec_t *forward_vel_press,
                                ymir_vec_t *obs_stress,
                                ymir_vec_t *weight,
                                rhea_stokes_problem_t *stokes_problem)
{
  /* check input */
  RHEA_ASSERT (rhea_stress_check_vec_type (misfit_stress));
  RHEA_ASSERT (forward_vel_press == NULL ||
               rhea_velocity_pressure_check_vec_type (forward_vel_press));
  RHEA_ASSERT (forward_vel_press == NULL ||
               rhea_velocity_pressure_is_valid (forward_vel_press));
  RHEA_ASSERT (obs_stress == NULL || rhea_stress_check_vec_type (obs_stress));
  RHEA_ASSERT (obs_stress == NULL || rhea_stress_is_valid (obs_stress));
  RHEA_ASSERT (weight == NULL ||
               rhea_inversion_obs_stress_weight_check_vec_type (weight));
  RHEA_ASSERT (weight == NULL ||
               rhea_inversion_obs_stress_weight_is_valid (weight));

  /* compute stress */
  if (forward_vel_press != NULL) {
    rhea_stokes_problem_stress_compute (misfit_stress, forward_vel_press,
                                        stokes_problem);
  }
  else {
    ymir_vec_set_zero (misfit_stress);
  }

  /* compute difference between forward and observed stress */
  if (obs_stress != NULL) {
    ymir_vec_add (-1.0, obs_stress, misfit_stress);
  }

  /* apply weights */
  if (weight != NULL) {
    if (1 == weight->ndfields) {
      ymir_vec_multiply_in1 (weight, misfit_stress);
    }
    else {
      RHEA_ASSERT (6 == weight->ndfields);
      ymir_vec_multiply_in (weight, misfit_stress);
    }
  }
}

static void
rhea_inversion_obs_stress_diff_adjoint (
                                  ymir_vec_t *output_vel,
                                  ymir_vec_t *stress,
                                  ymir_vec_t *forward_vel_press,
                                  ymir_vec_t *weight,
                                  const rhea_inversion_obs_stress_t obs_type,
                                  rhea_stokes_problem_t *stokes_problem)
{
  /* apply weights */
  if (weight != NULL) {
    if (1 == weight->ndfields) {
      ymir_vec_multiply_in1 (weight, stress);
    }
    else {
      RHEA_ASSERT (6 == weight->ndfields);
      ymir_vec_multiply_in (weight, stress);
    }
  }

  /* compute output velocity from stress */
  //TODO output_vel = 2 * mu'(forward_vel) * stress - pressure
}

static double
_compute_volume_gauss (ymir_vec_t *vec)
{
  ymir_vec_t         *unit;
  double              volume;

  /* check input */
  RHEA_ASSERT (vec->node_type == YMIR_GAUSS_NODE);

  /* setup mass weighted unit vector */
  unit = ymir_vec_template (vec);
  ymir_vec_set_value (unit, 1.0);
  ymir_mass_apply_gauss (unit);

  /* compute inner product to integrate */
  volume = ymir_vec_innerprod (vec, unit);
  ymir_vec_destroy (unit);

  return volume;
}

static void
_extract_single_field (ymir_vec_t *field, ymir_vec_t *vec, const int fieldid)
{
  const ymir_locidx_t n_elements     = vec->K;
  const int           n_nodes_per_el = vec->Np;
  ymir_locidx_t       elid;
  int                 nodeid;

  /* check input */
  RHEA_ASSERT (1 == field->ndfields);

  /* extract a single field */
  for (elid = 0; elid < n_elements; elid++) {
    for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
      const double *V = ymir_dvec_index (vec, elid, nodeid, fieldid);
      double       *F = ymir_dvec_index (field, elid, nodeid, 0);

      F[0] = V[0];
    }
  }
}

static void
_compute_n_volumes_gauss (double *volume, ymir_vec_t *vec)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vec);
  ymir_vec_t         *unit, *field;
  const int           n_fields = vec->ndfields;
  int                 fieldid;

  /* check input */
  RHEA_ASSERT (vec->node_type == YMIR_GAUSS_NODE);

  /* setup mass weighted unit vector */
  unit = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
  ymir_vec_set_value (unit, 1.0);
  ymir_mass_apply_gauss (unit);

  /* compute inner product to integrate */
  field = ymir_dvec_new (ymir_mesh, 1, YMIR_GAUSS_NODE);
  for (fieldid = 0; fieldid < n_fields; fieldid++) {
    _extract_single_field (field, vec, fieldid);
    volume[fieldid] = ymir_vec_innerprod (field, unit);
  }
  ymir_vec_destroy (field);
  ymir_vec_destroy (unit);
}

double
rhea_inversion_obs_stress_misfit (double *misfit_comp,
                                  ymir_vec_t *forward_vel_press,
                                  ymir_vec_t *obs_stress,
                                  ymir_vec_t *weight,
                                  const rhea_inversion_obs_stress_t obs_type,
                                  const rhea_weakzone_label_t weak_label,
                                  rhea_stokes_problem_t *stokes_problem)
{
  ymir_mesh_t        *ymir_mesh;
  ymir_vec_t         *misfit_stress;
  double              misfit;

  /* return if nothing to do */
  if (RHEA_INVERSION_OBS_STRESS_NONE == obs_type) {
    return 0.0;
  }

  /* get ymir mesh */
  if (forward_vel_press != NULL) {
    ymir_mesh = ymir_vec_get_mesh (forward_vel_press);
  }
  else if (obs_stress != NULL) {
    ymir_mesh = ymir_vec_get_mesh (obs_stress);
  }
  else {
    RHEA_ABORT_NOT_REACHED ();
  }

  /* create work vectors */
  misfit_stress = rhea_stress_new (ymir_mesh);

  /* compute misfit vector */
  rhea_inversion_obs_stress_diff (
      misfit_stress, forward_vel_press, obs_stress, weight, stokes_problem);

  /* compute misfit value(s) */
  switch (obs_type) {
  case RHEA_INVERSION_OBS_STRESS_VOLUME:
    {
      ymir_vec_t         *misfit_mass = rhea_stress_new (ymir_mesh);

      /* compute squared L2-norm */
      ymir_mass_apply (misfit_stress, misfit_mass);
      misfit = 0.5 * ymir_vec_innerprod (misfit_stress, misfit_mass);
      rhea_stress_destroy (misfit_mass);

      if (misfit_comp != NULL) {
        misfit_comp[0] = misfit;
      }
    }
    break;
  case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_NORMAL:
    RHEA_ASSERT (rhea_weakzone_label_is_valid (weak_label));
    RHEA_ASSERT (weight != NULL);
    RHEA_ASSERT (1 == weight->ndfields);
    /* note: assume weak zone indicator is implicitly given by `weight` */
    {
      rhea_weakzone_options_t *weak_options =
        rhea_stokes_problem_get_weakzone_options (stokes_problem);
      ymir_vec_t         *weak_normal, *stress_normal_normal;
      double              volume = _compute_volume_gauss (weight);

      /* compute normal direction at plate boundaries */
      weak_normal = rhea_weakzone_normal_new (ymir_mesh);
      rhea_weakzone_compute_normal (weak_normal, weak_options);

      /* compute normal component of the normal stress */
      stress_normal_normal = rhea_stress_normal_new (ymir_mesh);
      rhea_stress_normal_compute_normal (stress_normal_normal,
                                         misfit_stress, weak_normal);
      rhea_weakzone_normal_destroy (weak_normal);

      /* compute integral */
      misfit = _compute_volume_gauss (stress_normal_normal) / volume;
      rhea_stress_normal_destroy (stress_normal_normal);

      if (misfit_comp != NULL) {
        misfit_comp[0] = misfit;
      }
    }
    break;
  case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL:
    RHEA_ASSERT (rhea_weakzone_label_is_valid (weak_label));
    RHEA_ASSERT (weight != NULL);
    RHEA_ASSERT (1 == weight->ndfields);
    {
      rhea_weakzone_options_t *weak_options =
        rhea_stokes_problem_get_weakzone_options (stokes_problem);
      ymir_vec_t         *weak_normal, *stress_normal_tangential;
      double              volume = _compute_volume_gauss (weight);
      double              comp[3];
      int                 compid;

      /* compute normal direction at plate boundaries */
      weak_normal = rhea_weakzone_normal_new (ymir_mesh);
      rhea_weakzone_compute_normal (weak_normal, weak_options);

      /* compute tangential component of the normal stress */
      stress_normal_tangential = rhea_stress_tangential_new (ymir_mesh);
      rhea_stress_normal_compute_tangential (stress_normal_tangential,
                                             misfit_stress, weak_normal);
      rhea_weakzone_normal_destroy (weak_normal);

      /* compute integral for each component */
      _compute_n_volumes_gauss (comp, stress_normal_tangential);
      rhea_stress_tangential_destroy (stress_normal_tangential);
      for (compid = 0; compid < 3; compid++) {
        comp[compid] *= 1.0 / volume;
      }

      misfit = sqrt (comp[0]*comp[0] + comp[1]*comp[1] + comp[2]*comp[2]);
      if (misfit_comp != NULL) {
        for (compid = 0; compid < 3; compid++) {
          misfit_comp[compid] = comp[compid];
        }
      }
    }
    break;
  default: /* unknown observation type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* destroy */
  rhea_stress_destroy (misfit_stress);

  return misfit;
}

void
rhea_inversion_obs_stress_add_adjoint_rhs (
                                  ymir_vec_t *rhs_vel_mass,
                                  ymir_vec_t *forward_vel_press,
                                  ymir_vec_t *obs_stress,
                                  ymir_vec_t *weight,
                                  const rhea_inversion_obs_stress_t obs_type,
                                  rhea_stokes_problem_t *stokes_problem)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (rhs_vel_mass);
  ymir_vec_t         *misfit, *misfit_mass;

  /* check input */
  RHEA_ASSERT (rhea_velocity_check_vec_type (rhs_vel_mass));
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (forward_vel_press));
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (forward_vel_press));
  RHEA_ASSERT (obs_stress == NULL || rhea_stress_check_vec_type (obs_stress));
  RHEA_ASSERT (obs_stress == NULL || rhea_stress_is_valid (obs_stress));
  RHEA_ASSERT (weight == NULL ||
               rhea_inversion_obs_stress_weight_check_vec_type (weight));
  RHEA_ASSERT (weight == NULL ||
               rhea_inversion_obs_stress_weight_is_valid (weight));

  /* create work vectors */
  misfit      = rhea_stress_new (ymir_mesh);
  misfit_mass = rhea_stress_new (ymir_mesh);

#if 0
  /* compute misfit vector */
  rhea_inversion_obs_stress_diff (
      misfit, forward_vel, obs_stress, weight, obs_type, stokes_problem);

  /* apply mass matrix */
  ymir_mass_apply (misfit, misfit_mass);

  /* apply adjoint of misfit vector operators */
  rhea_inversion_obs_stress_diff_adjoint (
      misfit_mass, misfitf_mass, weight_surf, obs_type, domain_options);

  /* add output to right-hand side (change sign to obtain RHS) */
  ymir_vec_add (-1.0, misfit_vol_mass, rhs_vel_mass);
#endif

  /* destroy */
  rhea_stress_destroy (misfit);
  rhea_stress_destroy (misfit_mass);

}
