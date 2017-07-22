/*
 */

#include <rhea_stokes_norm.h>
#include <rhea_base.h>
#include <rhea_velocity_pressure.h>
#include <ymir_mass_vec.h>

static void
rhea_stokes_norm_l2_innerprod (double *innerprod_vel,
                               double *innerprod_press,
                               ymir_vec_t *arg_left,
                               ymir_vec_t *arg_right,
                               ymir_pressure_elem_t *press_elem)
{
  ymir_vec_t         *vel_left, *press_left;
  ymir_vec_t         *vel_right, *press_right;

  /* create component vectors */
  rhea_velocity_pressure_create_components (&vel_left, &press_left, arg_left,
                                            press_elem);
  rhea_velocity_pressure_create_components (&vel_right, &press_right, arg_right,
                                            press_elem);

  /* compute inner product for each component */
  if (innerprod_vel != NULL) {
    *innerprod_vel = ymir_vec_innerprod (vel_left, vel_right);
  }
  if (innerprod_press != NULL) {
    *innerprod_press = ymir_vec_innerprod (press_left, press_right);
  }

  /* destroy */
  ymir_vec_destroy (vel_left);
  ymir_vec_destroy (press_left);
  ymir_vec_destroy (vel_right);
  ymir_vec_destroy (press_right);
}

static void
rhea_stokes_norm_L2_primal_innerprod (double *innerprod_vel,
                                      double *innerprod_press,
                                      ymir_vec_t *arg_left,
                                      ymir_vec_t *arg_right,
                                      ymir_pressure_elem_t *press_elem)
{
  ymir_vec_t         *arg_right_mass = ymir_vec_template (arg_right);

  /* apply mass matrix */
  ymir_mass_apply (arg_right, arg_right_mass);

  /* compute inner product */
  rhea_stokes_norm_l2_innerprod(innerprod_vel, innerprod_press,
                                arg_left, arg_right_mass, press_elem);

  /* destroy */
  ymir_vec_destroy (arg_right_mass);
}

static void
rhea_stokes_norm_L2_dual_innerprod (double *innerprod_vel,
                                    double *innerprod_press,
                                    ymir_vec_t *arg_left,
                                    ymir_vec_t *arg_right,
                                    ymir_pressure_elem_t *press_elem)
{
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (arg_left);
  ymir_vec_t         *vel_left, *press_left;
  ymir_vec_t         *vel_right, *press_right;
  ymir_vec_t         *vel_right_lump, *press_right_lump;

  /* create component vectors */
  rhea_velocity_pressure_create_components (&vel_left, &press_left, arg_left,
                                            press_elem);
  rhea_velocity_pressure_create_components (&vel_right, &press_right, arg_right,
                                            press_elem);

  /* invert (lumped) mass matrices */
  if (innerprod_vel != NULL) {
    vel_right_lump = ymir_cvec_new (mesh, vel_right->ncfields);
    ymir_mass_lump (vel_right_lump);
    ymir_cvec_fabs (vel_right_lump, vel_right_lump);
    ymir_cvec_reciprocal (vel_right_lump);
    ymir_cvec_multiply_in (vel_right, vel_right_lump);
  }
  else {
    vel_right_lump = NULL;
  }
  if (innerprod_press != NULL) {
    press_right_lump = ymir_pressure_vec_new (mesh, press_elem);
    ymir_pressure_vec_lump_mass (press_right_lump, press_elem);
    ymir_vec_fabs (press_right_lump, press_right_lump);
    ymir_vec_reciprocal (press_right_lump);
    ymir_vec_multiply_in (press_right, press_right_lump);
  }
  else {
    press_right_lump = NULL;
  }

  /* compute inner product for each component */
  if (innerprod_vel != NULL) {
    *innerprod_vel = ymir_vec_innerprod (vel_left, vel_right_lump);
    ymir_vec_destroy (vel_right_lump);
  }
  if (innerprod_press != NULL) {
    *innerprod_press = ymir_vec_innerprod (press_left, press_right_lump);
    ymir_vec_destroy (press_right_lump);
  }

  /* destroy */
  ymir_vec_destroy (vel_left);
  ymir_vec_destroy (press_left);
  ymir_vec_destroy (vel_right);
  ymir_vec_destroy (press_right);
}

static void
rhea_stokes_norm_Hminus1_L2_innerprod (double *innerprod_vel,
                                       double *innerprod_press,
                                       ymir_vec_t *arg_left,
                                       ymir_vec_t *arg_right,
                                       ymir_pressure_elem_t *press_elem,
                                       ymir_Hminus1_norm_op_t *norm_op)
{
  ymir_vec_t         *vel_left, *vel_right;

  /* create component vectors */
  rhea_velocity_pressure_create_components (&vel_left, NULL, arg_left,
                                            press_elem);
  rhea_velocity_pressure_create_components (&vel_right, NULL, arg_right,
                                            press_elem);

  /* compute H^-1 inner product of velocity */
  if (innerprod_vel != NULL) {
    *innerprod_vel = ymir_Hminus1_norm_innerprod (vel_left, vel_right, norm_op);
  }

  /* compute L^2 inner product of pressure */
  if (innerprod_press != NULL) {
    rhea_stokes_norm_L2_primal_innerprod (NULL, innerprod_press,
                                          arg_left, arg_right, press_elem);
  }

  /* destroy */
  ymir_vec_destroy (vel_left);
  ymir_vec_destroy (vel_right);
}

void
rhea_stokes_norm_innerprod (double *innerprod_vel,
                            double *innerprod_press,
                            ymir_vec_t *arg_left,
                            ymir_vec_t *arg_right,
                            rhea_stokes_norm_type_t norm_type,
                            ymir_Hminus1_norm_op_t *norm_op,
                            ymir_pressure_elem_t *press_elem)
{
  RHEA_ASSERT (innerprod_vel != NULL);
  RHEA_ASSERT (innerprod_press != NULL);

  switch (norm_type) {
  case RHEA_STOKES_NORM_L2_VEC_SP:
    rhea_stokes_norm_l2_innerprod (innerprod_vel, innerprod_press,
                                   arg_left, arg_right, press_elem);
    break;

  case RHEA_STOKES_NORM_L2_PRIMAL:
    rhea_stokes_norm_L2_primal_innerprod (innerprod_vel, innerprod_press,
                                          arg_left, arg_right, press_elem);
    break;

  case RHEA_STOKES_NORM_L2_DUAL:
    rhea_stokes_norm_L2_dual_innerprod (innerprod_vel, innerprod_press,
                                        arg_left, arg_right, press_elem);
    break;

  case RHEA_STOKES_NORM_HMINUS1_L2:
    RHEA_ASSERT (norm_op != NULL);
    rhea_stokes_norm_Hminus1_L2_innerprod (innerprod_vel, innerprod_press,
                                           arg_left, arg_right,
                                           press_elem, norm_op);
    break;

  default: /* unknown norm type */
    RHEA_ABORT_NOT_REACHED ();
  }
}

double
rhea_stokes_norm_compute (double *norm_vel,
                          double *norm_press,
                          ymir_vec_t *vec,
                          rhea_stokes_norm_type_t norm_type,
                          ymir_Hminus1_norm_op_t *norm_op,
                          ymir_pressure_elem_t *press_elem)
{
  double              ip_vel, ip_press;

  /* compute inner product */
  rhea_stokes_norm_innerprod (&ip_vel, &ip_press, vec, vec,
                              norm_type, norm_op, press_elem);

  /* compute norms */
  if (norm_vel != NULL) {
    *norm_vel = sqrt (ip_vel);
  }
  if (norm_press != NULL) {
    *norm_press = sqrt (ip_press);
  }
  return sqrt (ip_vel + ip_press);
}
