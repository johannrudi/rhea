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
                                        stokes_problem, 0 /* !linearized */);
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
    if (weight->ndfields == 1) {
      ymir_vec_multiply_in1 (weight, misfit_stress);
    }
    else {
      RHEA_ASSERT (weight->ndfields == misfit_stress->ndfields);
      ymir_vec_multiply_in (weight, misfit_stress);
    }
  }
}

/**
 * Compute the adjoint of the d/(d vel) derivative of the observation operator:
 *   [d/(d vel) ObsOp(stress(vel,press)]^T * weight
 */
static void
rhea_inversion_obs_stress_diff_adjoint (
                                  ymir_vec_t *velocity,
                                  ymir_vec_t *strain_rate,
                                  ymir_vec_t *weight,
                                  rhea_stokes_problem_t *stokes_problem)
{
  /* apply weights */
  if (weight != NULL) {
    if (weight->ndfields == 1) {
      ymir_vec_multiply_in1 (weight, strain_rate);
    }
    else {
      RHEA_ASSERT (weight->ndfields == strain_rate->ndfields);
      ymir_vec_multiply_in (weight, strain_rate);
    }
  }

  /* compute divergence of viscous coefficient multiplied by strain rate */
  rhea_stokes_problem_stress_div_compute (
      velocity, strain_rate, stokes_problem, 1 /* linearized */);
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
  RHEA_ASSERT (0 <= fieldid && fieldid < vec->ndfields);

  /* extract a single field */
  for (elid = 0; elid < n_elements; elid++) {
    for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
      const double *V = ymir_dvec_index (vec, elid, nodeid, fieldid);
      double       *F = ymir_dvec_index (field, elid, nodeid, 0);

      F[0] = V[0];
    }
  }
}

//TODO deprecated
#if 0
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
#endif

double
rhea_inversion_obs_stress_misfit (ymir_vec_t *forward_vel_press,
                                  ymir_vec_t *obs_stress,
                                  ymir_vec_t *weight,
                                  const rhea_inversion_obs_stress_t obs_type,
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
      ymir_vec_copy (misfit_stress, misfit_mass);
      ymir_mass_apply_gauss (misfit_mass);
      misfit = 0.5 * ymir_vec_innerprod (misfit_stress, misfit_mass);
      rhea_stress_destroy (misfit_mass);
    }
    break;
  case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_NORMAL:
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
    }
    break;
  case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_0:
  case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_1:
  case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_2:
    RHEA_ASSERT (weight != NULL);
    RHEA_ASSERT (1 == weight->ndfields);
    /* note: assume weak zone indicator is implicitly given by `weight` */
    {
      rhea_weakzone_options_t *weak_options =
        rhea_stokes_problem_get_weakzone_options (stokes_problem);
      ymir_vec_t         *weak_normal, *stress_normal_tangential, *stress_comp;
      double              volume = _compute_volume_gauss (weight);

      /* compute normal direction at plate boundaries */
      weak_normal = rhea_weakzone_normal_new (ymir_mesh);
      rhea_weakzone_compute_normal (weak_normal, weak_options);

      /* compute tangential component of the normal stress */
      stress_normal_tangential = rhea_stress_tangential_new (ymir_mesh);
      rhea_stress_normal_compute_tangential (stress_normal_tangential,
                                             misfit_stress, weak_normal);
      rhea_weakzone_normal_destroy (weak_normal);

      /* extract component */
      stress_comp = rhea_stress_normal_new (ymir_mesh);
      switch (obs_type) {
      case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_0:
        _extract_single_field (stress_comp, stress_normal_tangential, 0);
        break;
      case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_1:
        _extract_single_field (stress_comp, stress_normal_tangential, 1);
        break;
      case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_2:
        _extract_single_field (stress_comp, stress_normal_tangential, 2);
        break;
      default: /* unknown observation type */
        RHEA_ABORT_NOT_REACHED ();
      }
      rhea_stress_tangential_destroy (stress_normal_tangential);

      /* compute integral */
      misfit = _compute_volume_gauss (stress_comp) / volume;
      rhea_stress_normal_destroy (stress_comp);
    }
    break;
  default: /* unknown observation type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* destroy */
  rhea_stress_destroy (misfit_stress);

  return misfit;
}

static void
_set_normal_normal_tensor (ymir_vec_t *tensor, ymir_vec_t *normal)
{
  const ymir_locidx_t n_elements     = normal->K;
  const int           n_nodes_per_el = normal->Np;
  ymir_locidx_t       elid;
  int                 nodeid;

  /* check input */
  RHEA_ASSERT (6 == tensor->ndfields);
  RHEA_ASSERT (3 == normal->ndfields);

  /* extract a single field */
  for (elid = 0; elid < n_elements; elid++) {
    for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
      const double *N = ymir_dvec_index (normal, elid, nodeid, 0);
      double       *T = ymir_dvec_index (tensor, elid, nodeid, 0);

      /* compute T = N * N^T */
      T[0] = N[0]*N[0];  T[1] = N[0]*N[1];  T[2] = N[0]*N[2];
                         T[3] = N[1]*N[1];  T[4] = N[1]*N[2];
                                            T[5] = N[2]*N[2];
    }
  }
}

static void
_set_tangential_normal_tensor (ymir_vec_t *tensor, ymir_vec_t *normal,
                               const int fieldid)
{
  const ymir_locidx_t n_elements     = normal->K;
  const int           n_nodes_per_el = normal->Np;
  ymir_locidx_t       elid;
  int                 nodeid;
  double              E[3], V[3];

  /* check input */
  RHEA_ASSERT (9 == tensor->ndfields);
  RHEA_ASSERT (3 == normal->ndfields);
  RHEA_ASSERT (0 <= fieldid && fieldid < 3);

  /* set unit normal vector */
  E[0] = 0.0;
  E[1] = 0.0;
  E[2] = 0.0;
  E[fieldid] = 1.0;

  /* extract a single field */
  for (elid = 0; elid < n_elements; elid++) {
    for (nodeid = 0; nodeid < n_nodes_per_el; nodeid++) {
      const double *N = ymir_dvec_index (normal, elid, nodeid, 0);
      double       *T = ymir_dvec_index (tensor, elid, nodeid, 0);

      /* compute T = (I - N * N^T) * (E * N^T) = E * N^T - N_i * N * N^T */
      /*                                       = (E - N_i * N) * N^T */
      V[0] = E[0] - N[fieldid] * N[0];
      V[1] = E[1] - N[fieldid] * N[1];
      V[2] = E[2] - N[fieldid] * N[2];
      T[0] = V[0]*N[0];  T[1] = V[0]*N[1];  T[2] = V[0]*N[2];
      T[3] = V[1]*N[0];  T[4] = V[1]*N[1];  T[5] = V[1]*N[2];
      T[6] = V[2]*N[0];  T[7] = V[2]*N[1];  T[8] = V[2]*N[2];
    }
  }
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
  ymir_vec_t         *misfit_mass, *rhs_add;

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
  switch (obs_type) {
  case RHEA_INVERSION_OBS_STRESS_VOLUME:
  case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_NORMAL:
    misfit_mass = rhea_stress_new (ymir_mesh);
    break;
  case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_0:
  case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_1:
  case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_2:
    misfit_mass = rhea_stress_nonsymmetric_new (ymir_mesh);
    break;
  default: /* unknown observation type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* compute adjoint of misfit operator */
  switch (obs_type) {
  case RHEA_INVERSION_OBS_STRESS_VOLUME:
    /* compute misfit vector */
    rhea_inversion_obs_stress_diff (
        misfit_mass, forward_vel_press, obs_stress, weight, stokes_problem);

    /* apply mass matrix */
    ymir_mass_apply_gauss (misfit_mass);
    break;
  case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_NORMAL:
    RHEA_ASSERT (weight != NULL);
    RHEA_ASSERT (1 == weight->ndfields);
    /* note: assume weak zone indicator is implicitly given by `weight` */
    {
      rhea_weakzone_options_t *weak_options =
        rhea_stokes_problem_get_weakzone_options (stokes_problem);
      ymir_vec_t         *weak_normal;
      double              volume = _compute_volume_gauss (weight);

      /* compute normal direction at plate boundaries */
      weak_normal = rhea_weakzone_normal_new (ymir_mesh);
      rhea_weakzone_compute_normal (weak_normal, weak_options);

      /* set (normal,normal) tensor (like rank-1 outer product) */
      _set_normal_normal_tensor (misfit_mass, weak_normal);

      /* scale by volume */
      ymir_vec_scale (1.0/volume, misfit_mass);

      /* apply mass matrix */
      ymir_mass_apply_gauss (misfit_mass);
    }
    break;
  case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_0:
  case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_1:
  case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_2:
    RHEA_ASSERT (weight != NULL);
    RHEA_ASSERT (1 == weight->ndfields);
    /* note: assume weak zone indicator is implicitly given by `weight` */
    {
      rhea_weakzone_options_t *weak_options =
        rhea_stokes_problem_get_weakzone_options (stokes_problem);
      ymir_vec_t         *weak_normal;
      double              volume = _compute_volume_gauss (weight);

      /* compute normal direction at plate boundaries */
      weak_normal = rhea_weakzone_normal_new (ymir_mesh);
      rhea_weakzone_compute_normal (weak_normal, weak_options);

      /* set (tangential,normal) tensor */
      switch (obs_type) {
      case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_0:
        _set_tangential_normal_tensor (misfit_mass, weak_normal, 0);
        break;
      case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_1:
        _set_tangential_normal_tensor (misfit_mass, weak_normal, 1);
        break;
      case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_2:
        _set_tangential_normal_tensor (misfit_mass, weak_normal, 2);
        break;
      default: /* unknown observation type */
        RHEA_ABORT_NOT_REACHED ();
      }

      /* scale by volume */
      ymir_vec_scale (1.0/volume, misfit_mass);

      /* apply mass matrix */
      ymir_mass_apply_gauss (misfit_mass);
    }
    break;
  default: /* unknown observation type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* apply adjoint of difference operator */
  rhs_add = rhea_velocity_new (ymir_mesh);
  rhea_inversion_obs_stress_diff_adjoint (
      rhs_add, misfit_mass, weight, stokes_problem);

  /* destroy */
  switch (obs_type) {
  case RHEA_INVERSION_OBS_STRESS_VOLUME:
  case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_NORMAL:
    rhea_stress_destroy (misfit_mass);
    break;
  case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_0:
  case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_1:
  case RHEA_INVERSION_OBS_STRESS_QOI_PLATE_BOUNDARY_TANGENTIAL_2:
    rhea_stress_nonsymmetric_destroy (misfit_mass);
    break;
  default: /* unknown observation type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* add output to right-hand side (change sign to obtain RHS) */
  ymir_vec_add (-1.0, rhs_add, rhs_vel_mass);
  rhea_velocity_destroy (rhs_add);
}
