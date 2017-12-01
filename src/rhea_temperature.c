/*
 */

#include <rhea_temperature.h>
#include <rhea_base.h>
#include <rhea_io_mpi.h>
#include <ymir_vec_getset.h>
#include <ymir_comm.h>

/* constant: seconds in a year (= 365.25*24*3600) */
#define RHEA_TEMPERATURE_SECONDS_PER_YEAR (31557600.0)

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_TEMPERATURE_DEFAULT_TYPE_NAME "NONE"
#define RHEA_TEMPERATURE_DEFAULT_SCALE (1.0)
#define RHEA_TEMPERATURE_DEFAULT_SHIFT (0.0)
#define RHEA_TEMPERATURE_DEFAULT_COLD_PLATE_AGE_YR (50.0e6)        /* [yr] */
#define RHEA_TEMPERATURE_DEFAULT_THERMAL_DIFFUSIVITY_M2_S (1.0e-6) /* [m^2/s] */
#define RHEA_TEMPERATURE_DEFAULT_RHS_SCALING (1.0)
#define RHEA_TEMPERATURE_DEFAULT_DATA_FILE_PATH_BIN NULL
#define RHEA_TEMPERATURE_DEFAULT_DATA_FILE_PATH_TXT NULL
#define RHEA_TEMPERATURE_DEFAULT_WRITE_DATA_FILE_PATH_BIN NULL

#define RHEA_TEMPERATURE_SINKER_DEFAULT_ACTIVE (0)
#define RHEA_TEMPERATURE_SINKER_DEFAULT_RANDOM_COUNT (0)
#define RHEA_TEMPERATURE_SINKER_DEFAULT_DECAY (100.0)
#define RHEA_TEMPERATURE_SINKER_DEFAULT_WIDTH (0.1)
#define RHEA_TEMPERATURE_SINKER_DEFAULT_SCALING (1.0)
#define RHEA_TEMPERATURE_SINKER_DEFAULT_CENTER_X (0.5)
#define RHEA_TEMPERATURE_SINKER_DEFAULT_CENTER_Y (0.5)
#define RHEA_TEMPERATURE_SINKER_DEFAULT_CENTER_Z (0.5)
#define RHEA_TEMPERATURE_SINKER_DEFAULT_DILATATION (1.0)
#define RHEA_TEMPERATURE_SINKER_DEFAULT_TRANSLATION_X (0.0)
#define RHEA_TEMPERATURE_SINKER_DEFAULT_TRANSLATION_Y (0.0)
#define RHEA_TEMPERATURE_SINKER_DEFAULT_TRANSLATION_Z (0.0)

#define RHEA_TEMPERATURE_PLUME_DEFAULT_ACTIVE (0)
#define RHEA_TEMPERATURE_PLUME_DEFAULT_RANDOM_COUNT (0)
#define RHEA_TEMPERATURE_PLUME_DEFAULT_DECAY (100.0)
#define RHEA_TEMPERATURE_PLUME_DEFAULT_WIDTH (0.1)
#define RHEA_TEMPERATURE_PLUME_DEFAULT_SCALING (1.0)
#define RHEA_TEMPERATURE_PLUME_DEFAULT_CENTER_X (0.5)
#define RHEA_TEMPERATURE_PLUME_DEFAULT_CENTER_Y (0.5)
#define RHEA_TEMPERATURE_PLUME_DEFAULT_CENTER_Z (0.5)
#define RHEA_TEMPERATURE_PLUME_DEFAULT_DILATATION (1.0)
#define RHEA_TEMPERATURE_PLUME_DEFAULT_TRANSLATION_X (0.0)
#define RHEA_TEMPERATURE_PLUME_DEFAULT_TRANSLATION_Y (0.0)
#define RHEA_TEMPERATURE_PLUME_DEFAULT_TRANSLATION_Z (0.0)

/* initialize options */
char               *rhea_temperature_type_name =
  RHEA_TEMPERATURE_DEFAULT_TYPE_NAME;
double              rhea_temperature_scale = RHEA_TEMPERATURE_DEFAULT_SCALE;
double              rhea_temperature_shift = RHEA_TEMPERATURE_DEFAULT_SHIFT;
double              rhea_temperature_cold_plate_age_yr =
  RHEA_TEMPERATURE_DEFAULT_COLD_PLATE_AGE_YR;
double              rhea_temperature_thermal_diffusivity_m2_s =
  RHEA_TEMPERATURE_DEFAULT_THERMAL_DIFFUSIVITY_M2_S;
double              rhea_temperature_rhs_scaling =
  RHEA_TEMPERATURE_DEFAULT_RHS_SCALING;
char               *rhea_temperature_data_file_path_bin =
  RHEA_TEMPERATURE_DEFAULT_DATA_FILE_PATH_BIN;
char               *rhea_temperature_data_file_path_txt =
  RHEA_TEMPERATURE_DEFAULT_DATA_FILE_PATH_TXT;
char               *rhea_temperature_write_data_file_path_bin =
  RHEA_TEMPERATURE_DEFAULT_WRITE_DATA_FILE_PATH_BIN;

int                 rhea_temperature_sinker_active =
  RHEA_TEMPERATURE_SINKER_DEFAULT_ACTIVE;
int                 rhea_temperature_sinker_random_count =
  RHEA_TEMPERATURE_SINKER_DEFAULT_RANDOM_COUNT;
double              rhea_temperature_sinker_decay =
  RHEA_TEMPERATURE_SINKER_DEFAULT_DECAY;
double              rhea_temperature_sinker_width =
  RHEA_TEMPERATURE_SINKER_DEFAULT_WIDTH;
double              rhea_temperature_sinker_scaling =
  RHEA_TEMPERATURE_SINKER_DEFAULT_SCALING;
double              rhea_temperature_sinker_center_x =
  RHEA_TEMPERATURE_SINKER_DEFAULT_CENTER_X;
double              rhea_temperature_sinker_center_y =
  RHEA_TEMPERATURE_SINKER_DEFAULT_CENTER_Y;
double              rhea_temperature_sinker_center_z =
  RHEA_TEMPERATURE_SINKER_DEFAULT_CENTER_Z;
double              rhea_temperature_sinker_dilatation =
  RHEA_TEMPERATURE_SINKER_DEFAULT_DILATATION;
double              rhea_temperature_sinker_translation_x =
  RHEA_TEMPERATURE_SINKER_DEFAULT_TRANSLATION_X;
double              rhea_temperature_sinker_translation_y =
  RHEA_TEMPERATURE_SINKER_DEFAULT_TRANSLATION_Y;
double              rhea_temperature_sinker_translation_z =
  RHEA_TEMPERATURE_SINKER_DEFAULT_TRANSLATION_Z;

int                 rhea_temperature_plume_active =
  RHEA_TEMPERATURE_PLUME_DEFAULT_ACTIVE;
int                 rhea_temperature_plume_random_count =
  RHEA_TEMPERATURE_PLUME_DEFAULT_RANDOM_COUNT;
double              rhea_temperature_plume_decay =
  RHEA_TEMPERATURE_PLUME_DEFAULT_DECAY;
double              rhea_temperature_plume_width =
  RHEA_TEMPERATURE_PLUME_DEFAULT_WIDTH;
double              rhea_temperature_plume_scaling =
  RHEA_TEMPERATURE_PLUME_DEFAULT_SCALING;
double              rhea_temperature_plume_center_x =
  RHEA_TEMPERATURE_PLUME_DEFAULT_CENTER_X;
double              rhea_temperature_plume_center_y =
  RHEA_TEMPERATURE_PLUME_DEFAULT_CENTER_Y;
double              rhea_temperature_plume_center_z =
  RHEA_TEMPERATURE_PLUME_DEFAULT_CENTER_Z;
double              rhea_temperature_plume_dilatation =
  RHEA_TEMPERATURE_PLUME_DEFAULT_DILATATION;
double              rhea_temperature_plume_translation_x =
  RHEA_TEMPERATURE_PLUME_DEFAULT_TRANSLATION_X;
double              rhea_temperature_plume_translation_y =
  RHEA_TEMPERATURE_PLUME_DEFAULT_TRANSLATION_Y;
double              rhea_temperature_plume_translation_z =
  RHEA_TEMPERATURE_PLUME_DEFAULT_TRANSLATION_Z;

void
rhea_temperature_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Temperature";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_S, "type", '\0',
    &(rhea_temperature_type_name), RHEA_TEMPERATURE_DEFAULT_TYPE_NAME,
    "Type of temperature: NONE, data, cold_plate, 2plates_poly2",

  YMIR_OPTIONS_D, "scale", '\0',
    &(rhea_temperature_scale), RHEA_TEMPERATURE_DEFAULT_SCALE,
    "Scaling factor multiplied to temperature",
  YMIR_OPTIONS_D, "shift", '\0',
    &(rhea_temperature_shift), RHEA_TEMPERATURE_DEFAULT_SHIFT,
    "Shift added to temperature",

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

  YMIR_OPTIONS_S, "data-file-path-bin", '\0',
    &(rhea_temperature_data_file_path_bin),
    RHEA_TEMPERATURE_DEFAULT_DATA_FILE_PATH_BIN,
    "Path to a binary file that contains a temperature field",
  YMIR_OPTIONS_S, "data-file-path-txt", '\0',
    &(rhea_temperature_data_file_path_txt),
    RHEA_TEMPERATURE_DEFAULT_DATA_FILE_PATH_TXT,
    "Path to a text file that contains a temperature field",
  YMIR_OPTIONS_S, "write-data-file-path-bin", '\0',
    &(rhea_temperature_write_data_file_path_bin),
    RHEA_TEMPERATURE_DEFAULT_WRITE_DATA_FILE_PATH_BIN,
    "Output path for a binary file containing a temperature field",

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
  YMIR_OPTIONS_I, "random-count", '\0',
    &(rhea_temperature_sinker_random_count),
    RHEA_TEMPERATURE_SINKER_DEFAULT_RANDOM_COUNT,
    "Number of sinkers to place randomly inside domain",

  YMIR_OPTIONS_D, "decay", '\0',
    &(rhea_temperature_sinker_decay), RHEA_TEMPERATURE_SINKER_DEFAULT_DECAY,
    "Decay of temperature",
  YMIR_OPTIONS_D, "width", '\0',
    &(rhea_temperature_sinker_width), RHEA_TEMPERATURE_SINKER_DEFAULT_WIDTH,
    "Width of inner sphere with maximum value",
  YMIR_OPTIONS_D, "scaling", '\0',
    &(rhea_temperature_sinker_scaling), RHEA_TEMPERATURE_SINKER_DEFAULT_SCALING,
    "Scaling factor",

  YMIR_OPTIONS_D, "center-x", '\0',
    &(rhea_temperature_sinker_center_x),
    RHEA_TEMPERATURE_SINKER_DEFAULT_CENTER_X,
    "Center of sigle sinker: x-coordinate",
  YMIR_OPTIONS_D, "center-y", '\0',
    &(rhea_temperature_sinker_center_y),
    RHEA_TEMPERATURE_SINKER_DEFAULT_CENTER_Y,
    "Center of sigle sinker: y-coordinate",
  YMIR_OPTIONS_D, "center-z", '\0',
    &(rhea_temperature_sinker_center_z),
    RHEA_TEMPERATURE_SINKER_DEFAULT_CENTER_Z,
    "Center of sigle sinker: z-coordinate",

  YMIR_OPTIONS_D, "dilatation", '\0',
    &(rhea_temperature_sinker_dilatation),
    RHEA_TEMPERATURE_SINKER_DEFAULT_DILATATION,
    "Move random center points inward or outward relative to domain center",
  YMIR_OPTIONS_D, "translation-x", '\0',
    &(rhea_temperature_sinker_translation_x),
    RHEA_TEMPERATURE_SINKER_DEFAULT_TRANSLATION_X,
    "Move random center points in x-direction",
  YMIR_OPTIONS_D, "translation-y", '\0',
    &(rhea_temperature_sinker_translation_y),
    RHEA_TEMPERATURE_SINKER_DEFAULT_TRANSLATION_Y,
    "Move random center points in y-direction",
  YMIR_OPTIONS_D, "translation-z", '\0',
    &(rhea_temperature_sinker_translation_z),
    RHEA_TEMPERATURE_SINKER_DEFAULT_TRANSLATION_Z,
    "Move random center points in z-direction",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);
}

void
rhea_temperature_add_options_plume (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "TemperaturePlume";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_B, "active", '\0',
    &(rhea_temperature_plume_active), RHEA_TEMPERATURE_PLUME_DEFAULT_ACTIVE,
    "Activate plume",
  YMIR_OPTIONS_I, "random-count", '\0',
    &(rhea_temperature_plume_random_count),
    RHEA_TEMPERATURE_PLUME_DEFAULT_RANDOM_COUNT,
    "Number of plumes to place randomly inside domain",

  YMIR_OPTIONS_D, "decay", '\0',
    &(rhea_temperature_plume_decay), RHEA_TEMPERATURE_PLUME_DEFAULT_DECAY,
    "Decay of temperature",
  YMIR_OPTIONS_D, "width", '\0',
    &(rhea_temperature_plume_width), RHEA_TEMPERATURE_PLUME_DEFAULT_WIDTH,
    "Width of inner sphere with minimum value",
  YMIR_OPTIONS_D, "scaling", '\0',
    &(rhea_temperature_plume_scaling), RHEA_TEMPERATURE_PLUME_DEFAULT_SCALING,
    "Scaling factor",

  YMIR_OPTIONS_D, "center-x", '\0',
    &(rhea_temperature_plume_center_x),
    RHEA_TEMPERATURE_PLUME_DEFAULT_CENTER_X,
    "Center of single plume: x-coordinate",
  YMIR_OPTIONS_D, "center-y", '\0',
    &(rhea_temperature_plume_center_y),
    RHEA_TEMPERATURE_PLUME_DEFAULT_CENTER_Y,
    "Center of single plume: y-coordinate",
  YMIR_OPTIONS_D, "center-z", '\0',
    &(rhea_temperature_plume_center_z),
    RHEA_TEMPERATURE_PLUME_DEFAULT_CENTER_Z,
    "Center of single plume: z-coordinate",

  YMIR_OPTIONS_D, "dilatation", '\0',
    &(rhea_temperature_plume_dilatation),
    RHEA_TEMPERATURE_PLUME_DEFAULT_DILATATION,
    "Move random center points inward or outward relative to domain center",
  YMIR_OPTIONS_D, "translation-x", '\0',
    &(rhea_temperature_plume_translation_x),
    RHEA_TEMPERATURE_PLUME_DEFAULT_TRANSLATION_X,
    "Move random center points in x-direction",
  YMIR_OPTIONS_D, "translation-y", '\0',
    &(rhea_temperature_plume_translation_y),
    RHEA_TEMPERATURE_PLUME_DEFAULT_TRANSLATION_Y,
    "Move random center points in y-direction",
  YMIR_OPTIONS_D, "translation-z", '\0',
    &(rhea_temperature_plume_translation_z),
    RHEA_TEMPERATURE_PLUME_DEFAULT_TRANSLATION_Z,
    "Move random center points in z-direction",

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
  if (strcmp (rhea_temperature_type_name, "NONE") == 0) {
    opt->type = RHEA_TEMPERATURE_NONE;
  }
  else if (strcmp (rhea_temperature_type_name, "data") == 0) {
    opt->type = RHEA_TEMPERATURE_DATA;
  }
  else if (strcmp (rhea_temperature_type_name, "cold_plate") == 0) {
    opt->type = RHEA_TEMPERATURE_COLD_PLATE;
  }
  else { /* unknown temperature type */
    RHEA_ABORT ("Unknown temperature type");
  }

  /* set scaling and shifting values */
  opt->scale = rhea_temperature_scale;
  opt->shift = rhea_temperature_shift;

  /* set options for cold plate model */
  opt->cold_plate_model_plate_age_yr = rhea_temperature_cold_plate_age_yr;

  /* set thermal constants */
  opt->thermal_diffusivity_m2_s = rhea_temperature_thermal_diffusivity_m2_s;

  /* set sinker options */
  opt->sinker_active = rhea_temperature_sinker_active;
  opt->sinker_random_count = rhea_temperature_sinker_random_count;
  opt->sinker_decay = rhea_temperature_sinker_decay;
  opt->sinker_width = rhea_temperature_sinker_width;
  opt->sinker_scaling = rhea_temperature_sinker_scaling;
  opt->sinker_center_x = rhea_temperature_sinker_center_x;
  opt->sinker_center_y = rhea_temperature_sinker_center_y;
  opt->sinker_center_z = rhea_temperature_sinker_center_z;
  opt->sinker_dilatation = rhea_temperature_sinker_dilatation;
  opt->sinker_translation_x = rhea_temperature_sinker_translation_x;
  opt->sinker_translation_y = rhea_temperature_sinker_translation_y;
  opt->sinker_translation_z = rhea_temperature_sinker_translation_z;

  /* set plume options */
  opt->plume_active = rhea_temperature_plume_active;
  opt->plume_random_count = rhea_temperature_plume_random_count;
  opt->plume_decay = rhea_temperature_plume_decay;
  opt->plume_width = rhea_temperature_plume_width;
  opt->plume_scaling = rhea_temperature_plume_scaling;
  opt->plume_center_x = rhea_temperature_plume_center_x;
  opt->plume_center_y = rhea_temperature_plume_center_y;
  opt->plume_center_z = rhea_temperature_plume_center_z;
  opt->plume_dilatation = rhea_temperature_plume_dilatation;
  opt->plume_translation_x = rhea_temperature_plume_translation_x;
  opt->plume_translation_y = rhea_temperature_plume_translation_y;
  opt->plume_translation_z = rhea_temperature_plume_translation_z;

  /* set right-hand side options */
  opt->rhs_scaling = rhea_temperature_rhs_scaling;

  /* store domain options */
  opt->domain_options = domain_options;
}

/******************************************************************************
 * Vector
 *****************************************************************************/

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
 * Get & Set Values
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
  RHEA_ASSERT (rhea_temperature_check_vec_type (temp_vec));
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
 * Fundamental Formulas
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
 * Computes a smooth spherical indicator with values in the range [0,1] and
 * that is shaped like a Gaussian:
 *
 *   Indc(x) = exp ( -decay * max (0, ||center - x|| - width/2)^2 )
 *
 * where
 *   x      --- 3D coordinates (x,y,z)
 *   center --- coordinates of sphere center
 *   decay  --- exponential decay
 *   width  --- width of sphere with max indicator = 1
 */
static double
rhea_temperature_indicator (const double x,
                            const double y,
                            const double z,
                            const double center_x,
                            const double center_y,
                            const double center_z,
                            const double decay,
                            const double width)
{
  double              val;

  /* compute distance from center */
  val = sqrt ( (center_x - x) * (center_x - x) +
               (center_y - y) * (center_y - y) +
               (center_z - z) * (center_z - z) );

  /* subtract sphere width from distance */
  val = SC_MAX (0.0, val - width / 2.0);

  /* compute Gaussian indicator */
  val = exp (-decay * val * val);
  RHEA_ASSERT (0.0 <= val && val <= 1.0);

  /* return indicator */
  return val;
}

/* uniformly random numbers that define points in a unit cube */
const double rhea_temperature_random_point_x[64] = {
  0.88438893850986044, 0.14846502294340846, 0.44565621143957523,
  0.30385162277765243, 0.79848582952154912, 0.23687978322605474,
  0.97370490267498400, 0.86009887685416064, 0.98523663941483031,
  0.72034320600074808, 0.88763686588960411, 0.99217530249530528,
  0.90134812303451162, 0.10843642785551477, 0.56714436189799555,
  0.66251609843508175, 0.96199376691617278, 0.69631446341112724,
  0.89003623322821313, 0.11394864199518917, 0.65199716721700318,
  0.28182023756086683, 0.75591412090796117, 0.11393063867842590,
  0.05064649767054341, 0.63020512499974968, 0.60315632043989653,
  0.03542346068278846, 0.10804603661555889, 0.55114033580441879,
  0.87223550392159788, 0.45964201011136652, 0.45187460904250087,
  0.74090530289358758, 0.42719355110128809, 0.36835084162880050,
  0.82905590877880064, 0.65050763538006051, 0.87757393430765984,
  0.17991487850424015, 0.58109322757849147, 0.86462201253355300,
  0.52892246632842854, 0.54327994624985110, 0.44454216287347137,
  0.62945048609947929, 0.34530786221516552, 0.95381302519792688,
  0.77502781398196741, 0.29553419616487600, 0.78485459109314959,
  0.32102321708895110, 0.57068285048187828, 0.69913355930120900,
  0.44621561169122159, 0.67537531841659471, 0.74719694413071300,
  0.13183066547551237, 0.14573209949629140, 0.82232622218530615,
  0.49263859846178115, 0.53852557446902116, 0.03642551552548012,
  0.36503262530574743
};
const double rhea_temperature_random_point_y[64] = {
  0.43898966680559026, 0.61981592726087100, 0.84399951147493613,
  0.48329456773428559, 0.98748757518935293, 0.70223663266789038,
  0.97230555568845700, 0.40188339800808559, 0.55947740598727547,
  0.48403851681234178, 0.19873675009892122, 0.40235161693594812,
  0.99538177651788629, 0.03611403079535602, 0.96196465780532600,
  0.52331331841200956, 0.54020403751352009, 0.51971615792087811,
  0.33020224251402053, 0.31092271345360578, 0.06616012459832676,
  0.88006625972605623, 0.60329637581534845, 0.97856388515997794,
  0.46620183379063596, 0.23029916377767379, 0.59987909583285700,
  0.51381483503827652, 0.45987562894231393, 0.80540431700799764,
  0.05219215242890096, 0.95853359612900957, 0.33342818646086647,
  0.50679452405550252, 0.16869027871227116, 0.94181779523988429,
  0.62659097495527594, 0.72662953230651417, 0.01436214394598578,
  0.92629426819795000, 0.63715122252522471, 0.05595252940702089,
  0.69435056062244160, 0.70252027215346380, 0.08539783081264496,
  0.79617906649845815, 0.94681666716159452, 0.07359563580819994,
  0.91418782124377351, 0.15184572231448457, 0.27083150215131735,
  0.82956180439916749, 0.57182963622884753, 0.79625794328157706,
  0.46566240836511064, 0.90366452636353989, 0.26051150708354609,
  0.12350083231980058, 0.58504361526871551, 0.72290297326846520,
  0.65488289844966852, 0.28220516455822231, 0.32624457306031240,
  0.30914961889591730
};
const double rhea_temperature_random_point_z[64] = {
  0.78172261291716627, 0.26062367939442066, 0.19620491873125200,
  0.33781204387102370, 0.15904755472239451, 0.37547166269754739,
  0.64369804942140474, 0.63193079797745710, 0.93359191584220900,
  0.63903109268123026, 0.39536627395521018, 0.65885648107293115,
  0.65316328109469035, 0.61809123904697372, 0.74610546941888767,
  0.25989428148439742, 0.03027016341976874, 0.05903055566448412,
  0.22970119787111953, 0.22843232213573739, 0.27543137986091693,
  0.44433035869923465, 0.78326593740403627, 0.84859667547954900,
  0.32565327494803353, 0.57988497328128900, 0.44842788038391523,
  0.40773016177263700, 0.45088275367583042, 0.70085007278606060,
  0.21968131140534608, 0.79004535697943823, 0.05909531268008061,
  0.19992541313935464, 0.75169459085468726, 0.01717254091071729,
  0.53874651515721916, 0.09448856392057847, 0.29430262684623720,
  0.06818043682703453, 0.65126926476494229, 0.81685517377385020,
  0.21240489628183400, 0.95643455726046889, 0.05734014754125749,
  0.69119133285814349, 0.52019031821147721, 0.20703194807307400,
  0.78255064771029959, 0.84791052231378217, 0.22781070481613375,
  0.82218219460097586, 0.28601827249938028, 0.44158905615122279,
  0.27903918727461485, 0.90852589884837465, 0.68963784007410900,
  0.19090285324449119, 0.07336169008046678, 0.92585803801748578,
  0.89012347796168689, 0.97595751787996410, 0.97301362389282042,
  0.12091238458062825
};

/**
 * Computes an indicator value as a superposition of 1,...,64 randomly placed
 * Gaussian spheres.
 */
static double
rhea_temperature_indicator_random (const int n_spheres,
                                   const double x,
                                   const double y,
                                   const double z,
                                   const double decay,
                                   const double width,
                                   const double dilatation,
                                   const double translation_x,
                                   const double translation_y,
                                   const double translation_z)
{
  double              center_x, center_y, center_z, val;
  int                 i;

  /* check input */
  RHEA_ASSERT (0.0 < dilatation);
  RHEA_ASSERT (0 < n_spheres && n_spheres <= 64);

  /* compute temperature of each bubble; combine multiplicatively */
  val = 0.0;
  for (i = 0; i < n_spheres; i++) { /* loop over all spheres */
    /* get center point */
    center_x = rhea_temperature_random_point_x[i];
    center_y = rhea_temperature_random_point_y[i];
    center_z = rhea_temperature_random_point_z[i];

    /* dilate and translate center point (assumes cube domain) */
    center_x = dilatation * (center_x - 0.5) + 0.5 + translation_x;
    center_y = dilatation * (center_y - 0.5) + 0.5 + translation_y;
    center_z = dilatation * (center_z - 0.5) + 0.5 + translation_z;

    /* add indicator value */
    val += rhea_temperature_indicator (x, y, z, center_x, center_y, center_z,
                                       decay, width);
  }

  /* return anomaly factor that will be multiplied by the temperature */
  return SC_MAX (0.0, SC_MIN (val, 1.0));
}

/**
 * Computes the temperature anomaly of one (cold) sinker, with values in the
 * range [0,1], that is shaped like a Gaussian:
 *
 *   dT(x) = c * ( 1.0 - exp(-decay * max(0, ||center - x|| - width/2)^2) )
 *
 * where
 *   x      --- 3D coordinates (x,y,z)
 *   center --- coordinates of sinker center
 *   decay  --- exponential decay
 *   width  --- width of min temperature around the center
 *   c      --- scaling factor (typically 1)
 */
static double
rhea_temperature_sinker (const double x,
                         const double y,
                         const double z,
                         const double center_x,
                         const double center_y,
                         const double center_z,
                         const double decay,
                         const double width,
                         const double scaling)
{
  double              val;

  /* get indicator value */
  val = rhea_temperature_indicator (x, y, z, center_x, center_y, center_z,
                                    decay, width);

  /* transform indicator to sinker anomaly */
  val = scaling * (1.0 - val);

  /* return anomaly factor that will be multiplied by the temperature */
  return SC_MAX (0.0, SC_MIN (val, 1.0));
}

/**
 * Computes the temperature anomaly of 1,...,64 randomly placed sinkers.
 */
static double
rhea_temperature_sinker_random (const int n_sinkers,
                                const double x,
                                const double y,
                                const double z,
                                const double decay,
                                const double width,
                                const double scaling,
                                const double dilatation,
                                const double translation_x,
                                const double translation_y,
                                const double translation_z)
{
  double              val;

  /* get indicator value */
  val = rhea_temperature_indicator_random (n_sinkers, x, y, z, decay, width,
                                           dilatation, translation_x,
                                           translation_y, translation_z);

  /* transform indicator to sinker anomaly */
  val = scaling * (1.0 - val);
  RHEA_ASSERT (isfinite (val));

  /* return anomaly factor that will be multiplied by the temperature */
  return SC_MAX (0.0, SC_MIN (val, 1.0));
}

/**
 * Computes the temperature anomaly of one (hot) plume, with values in the
 * range [1,2], that is shaped like a Gaussian:
 *
 *   dT(x) = 1 + c * exp(-decay * max(0, ||center - x|| - width/2)^2)
 *
 * where
 *   x      --- 3D coordinates (x,y,z)
 *   center --- coordinates of plume center
 *   decay  --- exponential decay
 *   width  --- width of max temperature around the center
 *   c      --- scaling factor (typically 1)
 */
static double
rhea_temperature_plume (const double x,
                        const double y,
                        const double z,
                        const double center_x,
                        const double center_y,
                        const double center_z,
                        const double decay,
                        const double width,
                        const double scaling)
{
  double              val;

  /* get indicator value */
  val = rhea_temperature_indicator (x, y, z, center_x, center_y, center_z,
                                    decay, width);

  /* transform indicator to plume anomaly */
  val = 1.0 + scaling * val;

  /* return anomaly factor that will be multiplied by the temperature */
  return SC_MAX (1.0, SC_MIN (val, 2.0));
}

/**
 * Computes the temperature anomaly of 1,...,64 randomly placed plumes.
 */
static double
rhea_temperature_plume_random (const int n_plumes,
                               const double x,
                               const double y,
                               const double z,
                               const double decay,
                               const double width,
                               const double scaling,
                               const double dilatation,
                               const double translation_x,
                               const double translation_y,
                               const double translation_z)
{
  double              val;

  /* get indicator value */
  val = rhea_temperature_indicator_random (n_plumes, x, y, z, decay, width,
                                           dilatation, translation_x,
                                           translation_y, translation_z);

  /* transform indicator to plume anomaly */
  val = 1.0 + scaling * val;

  /* return anomaly factor that will be multiplied by the temperature */
  return SC_MAX (1.0, SC_MIN (val, 2.0));
}

/******************************************************************************
 * Temperature Computation
 *****************************************************************************/

/**
 * Sets multiplier for sinker anomalies and plume anomalies.
 */
static double
rhea_temperature_anomaly_multiplier (const double x, const double y,
                                     const double z,
                                     rhea_temperature_options_t *opt)
{
  double              mult = 1.0;

  /* multiply in sinker anomaly */
  if (opt->sinker_active) {
    const double        decay = opt->sinker_decay;
    const double        width = opt->sinker_width;
    const double        scaling = opt->sinker_scaling;

    if (0 < opt->sinker_random_count) { /* if multiple randomly placed */
      const double        dilat = opt->sinker_dilatation;
      const double        transl_x = opt->sinker_translation_x;
      const double        transl_y = opt->sinker_translation_y;
      const double        transl_z = opt->sinker_translation_z;

      mult *= rhea_temperature_sinker_random (
          opt->sinker_random_count, x, y, z,
          decay, width, scaling, dilat, transl_x, transl_y, transl_z);
    }
    else { /* if single sinker */
      const double        center_x = opt->sinker_center_x;
      const double        center_y = opt->sinker_center_y;
      const double        center_z = opt->sinker_center_z;

      mult *= rhea_temperature_sinker (
          x, y, z, center_x, center_y, center_z, decay, width, scaling);
    }
  }

  /* multiply in plume anomaly */
  if (opt->plume_active) {
    const double        decay = opt->plume_decay;
    const double        width = opt->plume_width;
    const double        scaling = opt->plume_scaling;

    if (0 < opt->plume_random_count) { /* if multiple randomly placed */
      const double        dilat = opt->plume_dilatation;
      const double        transl_x = opt->plume_translation_x;
      const double        transl_y = opt->plume_translation_y;
      const double        transl_z = opt->plume_translation_z;

      mult *= rhea_temperature_plume_random (
          opt->plume_random_count, x, y, z,
          decay, width, scaling, dilat, transl_x, transl_y, transl_z);
    }
    else { /* if single plume */
      const double        center_x = opt->plume_center_x;
      const double        center_y = opt->plume_center_y;
      const double        center_z = opt->plume_center_z;

      mult *= rhea_temperature_plume (
          x, y, z, center_x, center_y, center_z, decay, width, scaling);
    }
  }

  return mult;
}

/**
 * Callback function to multiply in anomalies.
 */
static void
rhea_temperature_anomaly_multiply_in_fn (double *temp, double x, double y,
                                         double z, ymir_locidx_t nodeid,
                                         void *data)
{
  rhea_temperature_options_t  *opt = data;

  *temp *= rhea_temperature_anomaly_multiplier (x, y, z, opt);
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
    temp = RHEA_TEMPERATURE_NEUTRAL_VALUE;
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
  default: /* unknown temperature type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* multiply in sinker/plume anomalies */
  temp *= rhea_temperature_anomaly_multiplier (x, y, z, opt);

  /* scale and shift */
  if (isfinite (opt->scale)) {
    temp *= opt->scale;
  }
  if (isfinite (opt->shift)) {
    temp += opt->shift;
  }

  /* check temperature for `nan` and `inf` */
  RHEA_ASSERT (isfinite (temp));

  /* bound temperature to valid interval */
  temp = SC_MAX (0.0, SC_MIN (temp, 1.0));

  /* return temperature */
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

static void
rhea_temperature_read (ymir_vec_t *temperature)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (temperature);
  const mangll_locidx_t *n_nodes = ymir_mesh->cnodes->Ngo;
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  int                 mpisize, mpiret;

  char               *file_path_bin;
  char               *file_path_txt;
  double             *temp_data = ymir_cvec_index (temperature, 0, 0);
  int                *segment_offset;
  int                 r;

  /* set file paths */
  if (rhea_temperature_data_file_path_bin != NULL) {
    file_path_bin = rhea_temperature_data_file_path_bin;
    file_path_txt = NULL;
  }
  else if (rhea_temperature_data_file_path_txt != NULL) {
    file_path_bin = rhea_temperature_write_data_file_path_bin;
    file_path_txt = rhea_temperature_data_file_path_txt;
  }
  else { /* unknown file path */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* get parallel environment */
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize); SC_CHECK_MPI (mpiret);

  /* create segment offsets */
  segment_offset = RHEA_ALLOC (int, mpisize + 1);
  segment_offset[0] = 0;
  for (r = 0; r < mpisize; r++) {
    RHEA_ASSERT ((segment_offset[r] + n_nodes[r]) <= INT_MAX);
    segment_offset[r+1] = segment_offset[r] + (int) n_nodes[r];
  }

  /* read temperature */
  rhea_io_mpi_read_scatter_double (temp_data, segment_offset,
                                   file_path_bin, file_path_txt, mpicomm);
  RHEA_FREE (segment_offset);

  /* communicate shared node values */
  ymir_vec_share_owned (temperature);
}

void
rhea_temperature_compute (ymir_vec_t *temperature,
                          rhea_temperature_options_t *opt)
{
  switch (opt->type) {
  case RHEA_TEMPERATURE_DATA:
    /* read temperature */
    rhea_temperature_read (temperature);

    /* multiply in sinker/plume anomalies */
    ymir_cvec_set_function (temperature,
                            rhea_temperature_anomaly_multiply_in_fn, opt);

    /* scale and shift */
    if (isfinite (opt->scale)) {
      ymir_vec_scale (opt->scale, temperature);
    }
    if (isfinite (opt->shift)) {
      ymir_vec_shift (opt->shift, temperature);
    }

    /* bound temperature to valid interval */
    //TODO create separate fnc, use in amr
    ymir_vec_bound_min (temperature, 0.0);
    ymir_vec_bound_max (temperature, 1.0);
    break;
  case RHEA_TEMPERATURE_NONE:
  case RHEA_TEMPERATURE_COLD_PLATE:
    ymir_cvec_set_function (temperature, rhea_temperature_compute_fn, opt);
    break;
  default: /* unknown temperature type */
    RHEA_ABORT_NOT_REACHED ();
  }
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
    back_temp = RHEA_TEMPERATURE_NEUTRAL_VALUE;
    break;
  case RHEA_TEMPERATURE_DATA:
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
  default: /* unknown temperature type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* scale and shift */
  if (isfinite (opt->scale)) {
    back_temp *= opt->scale;
  }
  if (isfinite (opt->shift)) {
    back_temp += opt->shift;
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
