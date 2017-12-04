/**************************************************************
 * Non-zero Dirichlet boundary conditions
***************************************************************/

/* Dirichlet all */
static ymir_dir_code_t
slabs_set_vel_dir_all (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
     return YMIR_VEL_DIRICHLET_ALL;
}

static ymir_dir_code_t
slabs_set_vel_dir_all_2D (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
  if (face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE3 ||
      face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE4 ||
      face == RHEA_DOMAIN_BOUNDARY_FACE_BASE  ||
      face == RHEA_DOMAIN_BOUNDARY_FACE_TOP) {
      return YMIR_VEL_DIRICHLET_ALL;
  }
  else
    return YMIR_VEL_DIRICHLET_NORM;
}

/* free-slip */
static ymir_dir_code_t
slabs_set_vel_dir_freeslip (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
    return YMIR_VEL_DIRICHLET_NORM;
}

/* Dirichlet all on one side of the domain: SIDE3 (y=0).*/
static ymir_dir_code_t
slabs_set_vel_dir_inoutflow (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
  if (face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE3) {
    return YMIR_VEL_DIRICHLET_ALL;
  }
  else {
    return YMIR_VEL_DIRICHLET_NORM;
  }
}

/* Dirichlet all on both left and right sides of the domain: SIDE3 (y=0), and SIDE4 (y=ymax).*/
static ymir_dir_code_t
slabs_set_vel_dir_inoutflow_double (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
  if (face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE3 ||
      face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE4) {
    return YMIR_VEL_DIRICHLET_ALL;
  }
  else {
    return YMIR_VEL_DIRICHLET_NORM;
  }
}

/* Dirichlet all on one side of the domain: SIDE3 (y=0), and Neumann at the base*/
static ymir_dir_code_t
slabs_set_vel_dir_inoutflow_basefree (
    double X, double Y, double Z,
    double nx, double ny, double nz,
    ymir_topidx_t face, ymir_locidx_t node_id,
    void *data)
{
  if (face == RHEA_DOMAIN_BOUNDARY_FACE_SIDE3) {
    return YMIR_VEL_DIRICHLET_ALL;
  }
  else if (face == RHEA_DOMAIN_BOUNDARY_FACE_BASE) {
     return YMIR_VEL_DIRICHLET_NONE;
  }
  else {
    return YMIR_VEL_DIRICHLET_NORM;
  }
}

/* In-out flow sine velocity on one side of the domain. */
void
slabs_set_rhs_vel_nonzero_dir_inoutflow_sin (
    double *vel, double x, double y, double z,
    ymir_locidx_t nid, void *data)
{
  slabs_options_t  *slabs_options = data;
  const double     flow_scale = slabs_options->slabs_velbc_options->flow_scale;

  if (y < SC_1000_EPS) {
    vel[0] = 0.0;
    vel[1] = -flow_scale * sin (2.0 * M_PI * z);
    vel[2] = 0.0;
  }
  else {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }
}

/* In-out flow tanh velocity on one side of the domain. 2 layers */
void
slabs_set_rhs_vel_nonzero_dir_inoutflow_tanh_2layer (double *vel, double x, double y, double z,
                                              ymir_locidx_t nid, void *data)
{
  slabs_options_t  *slabs_options = data;
  slabs_velbc_options_t *velbc_options = slabs_options->slabs_velbc_options;
  const double        z_max = slabs_options->slabs_domain_options->z_max;
  const double        flow_scale = velbc_options->flow_scale;
  const double        zM = velbc_options->vel_dir_bc_middle;
  const double        a = z_max-zM, b = z_max-a;
  const double        shape = 2.0 * M_PI, scaling = 0.5*flow_scale, shift = 0.5*(b-a)*flow_scale;
  double txM = shape * (z-zM);

  if (y < SC_1000_EPS) {
    vel[0] = 0.0;
    vel[2] = 0.0;
    vel[1] = shift + scaling *
             ( (exp (txM) - exp (-txM)) /
             (exp (txM) + exp (-txM)) );
  }
  else {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }
}

/* In-out flow tanh velocity on one side of the domain.*/
void
slabs_set_rhs_vel_nonzero_dir_inoutflow_tanh_3layer (double *vel, double x, double y, double z,
                                              ymir_locidx_t nid, void *data)
{
  slabs_options_t  *slabs_options = data;
  slabs_velbc_options_t *velbc_options = slabs_options->slabs_velbc_options;
  const double        z_max = slabs_options->slabs_domain_options->z_max;
  const double        flow_scale = velbc_options->flow_scale;
  const double        zU = velbc_options->vel_dir_bc_upper;
  const double        zL = velbc_options->vel_dir_bc_lower;
  const double        a = zU-zL, b = z_max-a, c = 0.5*(zU+zL);
  const double        shape = 2.0 * M_PI, scaling = 0.5*(b+a)*flow_scale, shift = 0.5*(b-a)*flow_scale;
  double txL = shape*(z-zL),txU = shape*(z-zU);

  if (y < SC_1000_EPS) {
    vel[0] = 0.0;
    vel[2] = 0.0;
    if (z<=c)
      vel[1] = shift + scaling *
             ( (exp (txL) - exp (-txL)) /
             (exp (txL) + exp (-txL)) );
    else
      vel[1] = shift - scaling *
             ( (exp (txU) - exp (-txU)) /
             (exp (txU) + exp (-txU)) );
  }
  else {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }
}

/* In-out flow tanh velocity on both left and right sides of the domain.*/
void
slabs_set_rhs_vel_nonzero_dir_inoutflow_double_tanh_3layer (
    double *vel, double x, double y, double z,
    ymir_locidx_t nid, void *data)
{
  slabs_options_t  *slabs_options = data;
  slabs_velbc_options_t *velbc_options = slabs_options->slabs_velbc_options;
  const double        z_max = slabs_options->slabs_domain_options->z_max;
  const double        y_max = slabs_options->slabs_domain_options->y_max;
  const double        flow_scale = velbc_options->flow_scale;
  const double        zU = velbc_options->vel_dir_bc_upper;
  const double        zL = velbc_options->vel_dir_bc_lower;
  const double        a = zU-zL, b = z_max-a, c = 0.5*(zU+zL);
  const double        shape = 2.0 * M_PI, scaling = 0.5*(b+a)*flow_scale, shift = 0.5*(b-a)*flow_scale;
  double txL = shape*(z-zL),txU = shape*(z-zU);

  if (y < SC_1000_EPS) {
    vel[0] = 0.0;
    vel[2] = 0.0;
    if (z<=c)
      vel[1] = shift + scaling *
               ( (exp (txL) - exp (-txL)) /
               (exp (txL) + exp (-txL)) );
    else
      vel[1] = shift - scaling *
               ( (exp (txU) - exp (-txU)) /
               (exp (txU) + exp (-txU)) );
  }
  else if ((y_max - y) < SC_1000_EPS) {
    vel[0] = 0.0;
    vel[2] = 0.0;
    if (z<=c)
      vel[1] = -shift - scaling *
             ( (exp (txL) - exp (-txL)) /
             (exp (txL) + exp (-txL)) );
    else
      vel[1] = -shift + scaling *
             ( (exp (txU) - exp (-txU)) /
             (exp (txU) + exp (-txU)) );
  }
  else {
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
  }
}

void
slabs_vel_nonzero_dirichlet_compute ( ymir_vec_t * rhs_vel_nonzero_dirichlet,
                                      slabs_options_t * slabs_options)
{
  switch (slabs_options->slabs_velbc_options->vel_dir_bc) {
    case SLABS_VEL_DIR_BC_INOUTFLOW_SIN:
      rhea_domain_set_user_velocity_dirichlet_bc (
          slabs_set_vel_dir_inoutflow, NULL /* no data necessary */,
          0 /* TODO don't need this flag */);
      ymir_cvec_set_function (rhs_vel_nonzero_dirichlet,
                              slabs_set_rhs_vel_nonzero_dir_inoutflow_sin,
                              slabs_options);
      break;

    case SLABS_VEL_DIR_BC_INOUTFLOW_TANH_TWOLAYER:
      rhea_domain_set_user_velocity_dirichlet_bc (
          slabs_set_vel_dir_inoutflow, NULL /* no data necessary */,
          0 /* TODO don't need this flag */);
      ymir_cvec_set_function (rhs_vel_nonzero_dirichlet,
                              slabs_set_rhs_vel_nonzero_dir_inoutflow_tanh_2layer,
                              slabs_options);
      break;

    case SLABS_VEL_DIR_BC_INOUTFLOW_TANH_THREELAYER:
      rhea_domain_set_user_velocity_dirichlet_bc (
          slabs_set_vel_dir_inoutflow, NULL /* no data necessary */,
          0 /* TODO don't need this flag */);
      ymir_cvec_set_function (rhs_vel_nonzero_dirichlet,
                              slabs_set_rhs_vel_nonzero_dir_inoutflow_tanh_3layer,
                              slabs_options);
      break;

   default: /* BC not set */
      RHEA_ABORT_NOT_REACHED ();
  }
}


