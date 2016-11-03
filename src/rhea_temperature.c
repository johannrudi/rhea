/*
 */

#include <rhea_temperature.h>
#include <rhea_base.h>
#include <ymir_vec_getset.h>

/* constant: seconds in a year (= 365.25*24*3600) */
#define RHEA_TEMPERATURE_SECONDS_PER_YEAR (31557600.0)

/* default options */
#define RHEA_TEMPERATURE_DEFAULT_TYPE_NAME "none"
#define RHEA_TEMPERATURE_DEFAULT_COLD_PLATE_AGE_YR (50.0e6)        /* [yr] */
#define RHEA_TEMPERATURE_DEFAULT_THERMAL_DIFFUSIVITY_M2_S (1.0e-6) /* [m^2/s] */
#define RHEA_TEMPERATURE_DEFAULT_RHS_SCALING (1.0)
#define RHEA_TEMPERATURE_SINKER_DEFAULT_ACTIVE (0)
#define RHEA_TEMPERATURE_SINKER_DEFAULT_CENTER_X (0.5)
#define RHEA_TEMPERATURE_SINKER_DEFAULT_CENTER_Y (0.5)
#define RHEA_TEMPERATURE_SINKER_DEFAULT_CENTER_Z (0.5)
#define RHEA_TEMPERATURE_SINKER_DEFAULT_DECAY (100.0)
#define RHEA_TEMPERATURE_SINKER_DEFAULT_WIDTH (0.1)
#define RHEA_TEMPERATURE_SINKER_DEFAULT_SCALING (1.0)

/* initialize options */
char               *rhea_temperature_type_name =
  RHEA_TEMPERATURE_DEFAULT_TYPE_NAME;
double              rhea_temperature_cold_plate_age_yr =
  RHEA_TEMPERATURE_DEFAULT_COLD_PLATE_AGE_YR;
double              rhea_temperature_thermal_diffusivity_m2_s =
  RHEA_TEMPERATURE_DEFAULT_THERMAL_DIFFUSIVITY_M2_S;
double              rhea_temperature_rhs_scaling =
  RHEA_TEMPERATURE_DEFAULT_RHS_SCALING;
int                 rhea_temperature_sinker_active =
  RHEA_TEMPERATURE_SINKER_DEFAULT_ACTIVE;
double              rhea_temperature_sinker_center_x =
  RHEA_TEMPERATURE_SINKER_DEFAULT_CENTER_X;
double              rhea_temperature_sinker_center_y =
  RHEA_TEMPERATURE_SINKER_DEFAULT_CENTER_Y;
double              rhea_temperature_sinker_center_z =
  RHEA_TEMPERATURE_SINKER_DEFAULT_CENTER_Z;
double              rhea_temperature_sinker_decay =
  RHEA_TEMPERATURE_SINKER_DEFAULT_DECAY;
double              rhea_temperature_sinker_width =
  RHEA_TEMPERATURE_SINKER_DEFAULT_WIDTH;
double              rhea_temperature_sinker_scaling =
  RHEA_TEMPERATURE_SINKER_DEFAULT_SCALING;

void
rhea_temperature_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Temperature";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_S, "type", '\0',
    &(rhea_temperature_type_name), RHEA_TEMPERATURE_DEFAULT_TYPE_NAME,
    "Type of temperature: none, import, cold_plate, 2plates_poly2",

  YMIR_OPTIONS_D, "cold-plate-model-plate-age", '\0',
    &(rhea_temperature_cold_plate_age_yr),
    RHEA_TEMPERATURE_DEFAULT_COLD_PLATE_AGE_YR,
    "Cold plate model: plate age [yr] for half-space cooling model",

  YMIR_OPTIONS_D, "thermal-diffusivity", '\0',
    &(rhea_temperature_thermal_diffusivity_m2_s),
    RHEA_TEMPERATURE_DEFAULT_THERMAL_DIFFUSIVITY_M2_S,
    "Thermal diffusivity [m^2 / s]",

  YMIR_OPTIONS_D, "right-hand-side-scaling", '\0',
    &(rhea_temperature_rhs_scaling), RHEA_TEMPERATURE_DEFAULT_RHS_SCALING,
    "Scaling factor for velocity right-hand side from temperature",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);
}

void
rhea_temperature_add_options_sinker (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "TemperatureSinker";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_B, "active", '\0',
    &(rhea_temperature_sinker_active), RHEA_TEMPERATURE_SINKER_DEFAULT_ACTIVE,
    "Activate sinker",

  YMIR_OPTIONS_D, "center-x", '\0',
    &(rhea_temperature_sinker_center_x),
    RHEA_TEMPERATURE_SINKER_DEFAULT_CENTER_X,
    "Center of Sinker: x-coordinate",
  YMIR_OPTIONS_D, "center-y", '\0',
    &(rhea_temperature_sinker_center_y),
    RHEA_TEMPERATURE_SINKER_DEFAULT_CENTER_Y,
    "Center of Sinker: y-coordinate",
  YMIR_OPTIONS_D, "center-z", '\0',
    &(rhea_temperature_sinker_center_z),
    RHEA_TEMPERATURE_SINKER_DEFAULT_CENTER_Z,
    "Center of Sinker: z-coordinate",

  YMIR_OPTIONS_D, "decay", '\0',
    &(rhea_temperature_sinker_decay), RHEA_TEMPERATURE_SINKER_DEFAULT_DECAY,
    "Decay of temperature",
  YMIR_OPTIONS_D, "width", '\0',
    &(rhea_temperature_sinker_width), RHEA_TEMPERATURE_SINKER_DEFAULT_WIDTH,
    "Width of sinker where the max temperature is reached",

  YMIR_OPTIONS_D, "scaling", '\0',
    &(rhea_temperature_sinker_scaling), RHEA_TEMPERATURE_SINKER_DEFAULT_SCALING,
    "Scaling factor",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);
}

void
rhea_temperature_process_options (rhea_temperature_options_t *opt,
                                  rhea_domain_options_t *domain_options)
{
  /* set temperature type */
  if (strcmp (rhea_temperature_type_name, "none") == 0) {
    opt->type = RHEA_TEMPERATURE_NONE;
  }
  else if (strcmp (rhea_temperature_type_name, "import") == 0) {
    opt->type = RHEA_TEMPERATURE_IMPORT;
  }
  else if (strcmp (rhea_temperature_type_name, "cold_plate") == 0) {
    opt->type = RHEA_TEMPERATURE_COLD_PLATE;
  }
  else if (strcmp (rhea_temperature_type_name, "2plates_poly2") == 0) {
    opt->type = RHEA_TEMPERATURE_2PLATES_POLY2;
  }
  else { /* unknown temperature type */
    RHEA_ABORT ("Unknown temperature type");
  }

  /* set options for cold plate model */
  opt->cold_plate_model_plate_age_yr = rhea_temperature_cold_plate_age_yr;

  /* set thermal constants */
  opt->thermal_diffusivity_m2_s = rhea_temperature_thermal_diffusivity_m2_s;

  /* set sinker options */
  opt->sinker_active = rhea_temperature_sinker_active;
  opt->sinker_center_x = rhea_temperature_sinker_center_x;
  opt->sinker_center_y = rhea_temperature_sinker_center_y;
  opt->sinker_center_z = rhea_temperature_sinker_center_z;
  opt->sinker_decay = rhea_temperature_sinker_decay;
  opt->sinker_width = rhea_temperature_sinker_width;
  opt->sinker_scaling = rhea_temperature_sinker_scaling;

  /* set right-hand side options */
  opt->rhs_scaling = rhea_temperature_rhs_scaling;

  /* store domain options */
  opt->domain_options = domain_options;
}

ymir_vec_t *
rhea_temperature_new (ymir_mesh_t *ymir_mesh)
{
  return ymir_cvec_new (ymir_mesh, 1);
}

void
rhea_temperature_destroy (ymir_vec_t *temperature)
{
  ymir_vec_destroy (temperature);
}

int
rhea_temperature_check_vec_type (ymir_vec_t *vec)
{
  return (
      ymir_vec_is_cvec (vec) &&
      vec->ncfields == 1 &&
      vec->node_type == YMIR_GLL_NODE
  );
}

int
rhea_temperature_is_valid (ymir_vec_t *vec)
{
  return (
      sc_dmatrix_is_valid (vec->dataown) && sc_dmatrix_is_valid (vec->coff) &&
      0.0 <= ymir_cvec_min_global (vec) && ymir_cvec_max_global (vec) <= 1.0
  );
}

/******************************************************************************
 * Temperature Computation
 *****************************************************************************/

/**
 * Computes the nondimensional scaling factor that is necessary to calculate
 * the temperature distribution with the half-space cooling model (HSCM).
 *
 *   scaling = R / (2 * sqrt(kappa * t))
 *
 * where
 *   R     --- earth radius [m]
 *   kappa --- thermal diffusivity [m^2 / s]
 *   t     --- plate age [yr]
 */
static double
rhea_temperature_cold_plate_hscm_plate_age_to_scaling (
                                              const double radius_max_m,
                                              const double plate_age_yr,
                                              const double thermal_diffus_m2_s)
{
  const double        plate_age_s =
                        plate_age_yr * RHEA_TEMPERATURE_SECONDS_PER_YEAR;

  return radius_max_m / (2.0 * sqrt (plate_age_s * thermal_diffus_m2_s));
}

/**
 * Computes radially dependent temperature distribution according to the
 * half-space cooling model.
 *
 *   T(r) = erf( c * (r_t - r)),  for 0 <= r <= r_t
 *
 * where
 *   r_t --- nondimensional earth radius (usually 1)
 *   c   --- appropriate scaling factor
 *
 * Range(T) = [0, 1]
 */
static double
rhea_temperature_cold_plate_hscm (const double radius,
                                  const double radius_max,
                                  const double radius_max_m,
                                  const double plate_age_yr,
                                  const double thermal_diffus_m2_s)
{
  double              scaling_nondim;

  scaling_nondim = rhea_temperature_cold_plate_hscm_plate_age_to_scaling (
      radius_max_m, plate_age_yr, thermal_diffus_m2_s);

  return erf (scaling_nondim * (radius_max - radius));
}

/**
 * Computes the temperature of one sinker that is shaped like a Gaussian:
 *
 *   T(x) = c * ( 1.0 - exp(-decay * max(0, ||center - x|| - width/2)^2) )
 *
 * where
 *   x      --- 3D coordinates (x,y,z)
 *   center --- coordinates of sinker center
 *   decay  --- decay
 *   width  --- width of min temperature around the center
 *   c      --- scaling factor (typically 1)
 */
static double
rhea_temperature_sinker (const double x, const double y, const double z,
                         const double center_x, const double center_y,
                         const double center_z, const double decay,
                         const double width, const double scaling)
{
  double              val;

  /* compute distance from center */
  val = sqrt ( (center_x - x) * (center_x - x) +
               (center_y - y) * (center_y - y) +
               (center_z - z) * (center_z - z) );

  /* subtract bubble width */
  val = SC_MAX (0.0, val - width / 2.0);

  /* compute bubble temperature */
  val = 1.0 - exp (-decay * val * val);
  val *= scaling;

  /* return temperature after bounding to valid interval [0,1] */
  return SC_MAX (0.0, SC_MIN (val, 1.0));
}

/**
 * Computes the temperature at one node.
 */
static double
rhea_temperature_node (const double x, const double y, const double z,
                       rhea_temperature_options_t *opt)
{
  rhea_domain_options_t  *domain_options = opt->domain_options;
  double              temp;

  switch (opt->type) {
  case RHEA_TEMPERATURE_NONE:
    temp = RHEA_TEMPERATURE_DEFAULT_VALUE;
    break;
  case RHEA_TEMPERATURE_COLD_PLATE:
    {
      const double        radius_max = domain_options->radius_max;
      const double        radius_max_m = domain_options->radius_max_m;
      const double        plate_age_yr = opt->cold_plate_model_plate_age_yr;
      const double        thermal_diffus_m2_s = opt->thermal_diffusivity_m2_s;
      double              radius;

      radius = rhea_domain_compute_radius (x, y, z, domain_options);
      temp = rhea_temperature_cold_plate_hscm (
          radius, radius_max, radius_max_m, plate_age_yr, thermal_diffus_m2_s);
    }
    break;
//case RHEA_TEMPERATURE_2PLATES_POLY2:
//  break;
  default: /* unknown temperature type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* multiply in sinker temperature */
  if (opt->sinker_active) {
    temp *= rhea_temperature_sinker (
        x, y, z, opt->sinker_center_x, opt->sinker_center_y,
        opt->sinker_center_z, opt->sinker_decay, opt->sinker_width,
        opt->sinker_scaling);
  }

  /* check background temperature for `nan` and `inf` */
  RHEA_ASSERT (isfinite (temp));

  /* bound background temperature to valid interval */
  temp = SC_MAX (0.0, SC_MIN (temp, 1.0));

  /* return background temperature */
  return temp;
}

/**
 * Callback function to compute the temperature.
 */
static void
rhea_temperature_compute_fn (double *temp, double x, double y, double z,
                             ymir_locidx_t nodeid, void *data)
{
  rhea_temperature_options_t  *opt = data;

  *temp = rhea_temperature_node (x, y, z, opt);
}

void
rhea_temperature_compute (ymir_vec_t *temperature,
                          rhea_temperature_options_t *opt)
{
  ymir_cvec_set_function (temperature, rhea_temperature_compute_fn, opt);
}

/******************************************************************************
 * Background Temperature Computation
 *****************************************************************************/

/**
 * Computes the background temperature at one node.
 */
static double
rhea_temperature_background_node (const double x, const double y,
                                  const double z,
                                  rhea_temperature_options_t *opt)
{
  rhea_domain_options_t  *domain_options = opt->domain_options;
  double              back_temp;

  switch (opt->type) {
  case RHEA_TEMPERATURE_NONE:
    back_temp = RHEA_TEMPERATURE_DEFAULT_VALUE;
    break;
  case RHEA_TEMPERATURE_COLD_PLATE:
    {
      const double        radius_max = domain_options->radius_max;
      const double        radius_max_m = domain_options->radius_max_m;
      const double        plate_age_yr = opt->cold_plate_model_plate_age_yr;
      const double        thermal_diffus_m2_s = opt->thermal_diffusivity_m2_s;
      double              radius;

      radius = rhea_domain_compute_radius (x, y, z, domain_options);
      back_temp = rhea_temperature_cold_plate_hscm (
          radius, radius_max, radius_max_m, plate_age_yr, thermal_diffus_m2_s);
    }
    break;
//case RHEA_TEMPERATURE_2PLATES_POLY2:
//  break;
  default: /* unknown temperature type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* check background temperature for `nan` and `inf` */
  RHEA_ASSERT (isfinite (back_temp));

  /* bound background temperature to valid interval */
  back_temp = SC_MAX (0.0, SC_MIN (back_temp, 1.0));

  /* return background temperature */
  return back_temp;
}

/**
 * Callback function to compute the background temperature.
 */
static void
rhea_temperature_background_compute_fn (double *temp, double x, double y,
                                        double z, ymir_locidx_t nodeid,
                                        void *data)
{
  rhea_temperature_options_t  *opt = data;

  *temp = rhea_temperature_background_node (x, y, z, opt);
}

void
rhea_temperature_background_compute (ymir_vec_t *back_temperature,
                                     rhea_temperature_options_t *opt)
{
  ymir_cvec_set_function (back_temperature,
                          rhea_temperature_background_compute_fn, opt);
}

/******************************************************************************
 * Get & Set Functions
 *****************************************************************************/

/**
 * Restrict temperature to its valid range [0,1].
 */
static void
rhea_temperature_bound_range_data (double *_sc_restrict temp_data,
                                   const sc_bint_t size)
{
  sc_bint_t           i;

  for (i = 0; i < size; i++) {
    RHEA_ASSERT (isfinite (temp_data[i]));
    temp_data[i] = SC_MAX (0.0, SC_MIN (temp_data[i], 1.0));
  }
}

double *
rhea_temperature_get_elem_gauss (sc_dmatrix_t *temp_el_mat,
                                 ymir_vec_t *temp_vec,
                                 const ymir_locidx_t elid)
{
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (temp_vec);
  const int           n_nodes_per_el = ymir_mesh_get_num_nodes_per_elem (mesh);

  /* check input */
  YMIR_ASSERT_IS_CVEC (temp_vec);
  RHEA_ASSERT (temp_vec->ncfields == 1);
  RHEA_ASSERT (temp_el_mat->m == n_nodes_per_el);
  RHEA_ASSERT (temp_el_mat->n == 1);

  /* interpolate from continuous nodes to (discontinuous) Gauss nodes */
  ymir_cvec_get_elem_interp (temp_vec, temp_el_mat, YMIR_STRIDE_NODE, elid,
                             YMIR_GAUSS_NODE, YMIR_READ);

  /* restrict temperature to valid range to account for interpolation erros */
  rhea_temperature_bound_range_data (temp_el_mat->e[0], n_nodes_per_el);

  return temp_el_mat->e[0];
}

/******************************************************************************
 * Right-Hand Side Computation
 *****************************************************************************/

/**
 * Computes velocity right-hand side at one node from temperature and
 * background temperature.
 */
static void
rhea_temperature_rhs_vel_node (double *rhs, const double x, const double y,
                               const double z, const double temp,
                               const double back_temp,
                               rhea_domain_options_t *domain_options)
{
  const double        scaling = rhea_temperature_rhs_scaling;

  switch (domain_options->shape) {
  case RHEA_DOMAIN_CUBE:
  case RHEA_DOMAIN_BOX:
    /**
     * Computes right-hand side from temperature and background temperature
     * on the cube domain:
     *
     *   f(x) = e_z * (T - T_0)
     *
     * where `e_z` is unit vector in z-direction, `T` is temperature,
     * `T_0` is background temperature.
     */
    rhs[0] = 0.0;
    rhs[1] = 0.0;
    rhs[2] = scaling * (temp - back_temp);
    break;
  case RHEA_DOMAIN_SHELL:
  case RHEA_DOMAIN_CUBE_SPHERICAL:
  case RHEA_DOMAIN_BOX_SPHERICAL:
    /**
     * Computes right-hand side from temperature and background temperature
     * on the shell domain:
     *
     *   f(x) = e_r * (T - T_0)
     *
     * where `e_r` is normalized spherical position vector, `T` is temperature,
     * `T_0` is background temperature.
     */
    {
      const double        radius = rhea_domain_compute_radius (x, y, z,
                                                               domain_options);

      rhs[0] = scaling * x / radius * (temp - back_temp);
      rhs[1] = scaling * y / radius * (temp - back_temp);
      rhs[2] = scaling * z / radius * (temp - back_temp);
    }
    break;
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }
}

/* data for callback function to compute the velocity right-hand side */
typedef struct rhea_temperature_compute_rhs_vel_fn_data
{
  ymir_vec_t          *temperature;
  rhea_temperature_options_t  *temp_options;
}
rhea_temperature_compute_rhs_vel_fn_data_t;

/**
 * Callback function to compute the velocity right-hand side.
 */
static void
rhea_temperature_compute_rhs_vel_fn (double *rhs, double x, double y, double z,
                                     ymir_locidx_t nodeid, void *data)
{
  rhea_temperature_compute_rhs_vel_fn_data_t  *d = data;
  rhea_temperature_options_t  *opt = d->temp_options;
  rhea_domain_options_t  *domain_options = opt->domain_options;
  const double        temp = *ymir_cvec_index (d->temperature, nodeid, 0);
  double              back_temp;

  /* compute background temperature */
  back_temp = rhea_temperature_background_node (x, y, z, opt);

  /* compute right-hand side */
  rhea_temperature_rhs_vel_node (rhs, x, y, z, temp, back_temp, domain_options);
}

void
rhea_temperature_compute_rhs_vel (ymir_vec_t *rhs_vel,
                                  ymir_vec_t *temperature,
                                  rhea_temperature_options_t *opt)
{
  rhea_temperature_compute_rhs_vel_fn_data_t  data;

  /* set right-hand side */
  data.temperature = temperature;
  data.temp_options = opt;
  ymir_cvec_set_function (rhs_vel, rhea_temperature_compute_rhs_vel_fn, &data);
}
