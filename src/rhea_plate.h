/* RHEA_PLATE  Assigns plate labels and retrieves them from coordinates. */

#ifndef RHEA_PLATE_H
#define RHEA_PLATE_H

#include <rhea_domain.h>
#include <rhea_temperature.h>

/******************************************************************************
 * Plate Labels
 *****************************************************************************/

/* generic label for plate "none" */
#define RHEA_PLATE_NONE (-1)

/* plates of cube domain */
typedef enum
{
  RHEA_PLATE_CUBE_NONE = RHEA_PLATE_NONE, /* "none" (must come first) */
  RHEA_PLATE_CUBE_NW = 0,                 /* North-West (must start at zero) */
  RHEA_PLATE_CUBE_NE,                     /* North-East */
  RHEA_PLATE_CUBE_SE,                     /* South-East */
  RHEA_PLATE_CUBE_SW,                     /* South-West */
  RHEA_PLATE_CUBE_N                       /* (number of plates) */
}
rhea_plate_cube_label_t;

/* plates of earth domain */
typedef enum
{
  RHEA_PLATE_EARTH_NONE = RHEA_PLATE_NONE,  /* "none" (must come first) */

  /* MORVEL(25) (25 plates) */
  RHEA_PLATE_EARTH_AM = 0,  /* Amur (must start at zero) */
  RHEA_PLATE_EARTH_AN,      /* Antarctic */
  RHEA_PLATE_EARTH_AR,      /* Arabia */
  RHEA_PLATE_EARTH_AU,      /* Australia */
  RHEA_PLATE_EARTH_CP,      /* Capricorn */
  RHEA_PLATE_EARTH_CA,      /* Caribbean */
  RHEA_PLATE_EARTH_CO,      /* Cocos */
  RHEA_PLATE_EARTH_EU,      /* Eurasia */
  RHEA_PLATE_EARTH_IN,      /* India */
  RHEA_PLATE_EARTH_JF,      /* Juan de Fuca */
  RHEA_PLATE_EARTH_LW,      /* Lwandle */
  RHEA_PLATE_EARTH_MQ,      /* Macquarie */
  RHEA_PLATE_EARTH_NZ,      /* Nazca */
  RHEA_PLATE_EARTH_NA,      /* North America */
//RHEA_PLATE_EARTH_NU,      /* Nubia */
  RHEA_PLATE_EARTH_PA,      /* Pacific */
  RHEA_PLATE_EARTH_PS,      /* Philippine Sea */
  RHEA_PLATE_EARTH_RI,      /* Rivera */
  RHEA_PLATE_EARTH_SW,      /* Sandwich */
  RHEA_PLATE_EARTH_SC,      /* Scotia */
  RHEA_PLATE_EARTH_SM,      /* Somalia */
  RHEA_PLATE_EARTH_SA,      /* South America */
  RHEA_PLATE_EARTH_SU,      /* Sundaland */
  RHEA_PLATE_EARTH_SR,      /* Sur */
  RHEA_PLATE_EARTH_YZ,      /* Yangtze */

  /* Bird, 2003 (31 plates) */
  RHEA_PLATE_EARTH_AS,      /* Aegean Sea */
  RHEA_PLATE_EARTH_AP,      /* Altiplano */
  RHEA_PLATE_EARTH_AT,      /* Anatolia */
  RHEA_PLATE_EARTH_BR,      /* Bahnoral Reef */
  RHEA_PLATE_EARTH_BS,      /* Banda Sea */
  RHEA_PLATE_EARTH_BH,      /* Birds Head */
  RHEA_PLATE_EARTH_BU,      /* Burma */
  RHEA_PLATE_EARTH_CL,      /* Caroline */
//RHEA_PLATE_EARTH_CR,      /* Conway Reef */
  RHEA_PLATE_EARTH_EA,      /* Easter */
//RHEA_PLATE_EARTH_FT,      /* Futuna */
//RHEA_PLATE_EARTH_GP,      /* Galapagos */
//RHEA_PLATE_EARTH_JZ,      /* Juan Fernandez */
  RHEA_PLATE_EARTH_KE,      /* Kermadec */
//RHEA_PLATE_EARTH_MN,      /* Manus */
//RHEA_PLATE_EARTH_MO,      /* Maoke */
  RHEA_PLATE_EARTH_MA,      /* Mariana */
  RHEA_PLATE_EARTH_MS,      /* Molucca Sea */
  RHEA_PLATE_EARTH_NH,      /* New Hebrides */
//RHEA_PLATE_EARTH_NI,      /* Niuafo'ou */
  RHEA_PLATE_EARTH_ND,      /* North Andes */
  RHEA_PLATE_EARTH_NB,      /* North Bismarck */
  RHEA_PLATE_EARTH_OK,      /* Okhotsk */
  RHEA_PLATE_EARTH_ON,      /* Okinawa */
  RHEA_PLATE_EARTH_PM,      /* Panama */
//RHEA_PLATE_EARTH_SL,      /* Shetland */
  RHEA_PLATE_EARTH_SS,      /* Solomon Sea */
  RHEA_PLATE_EARTH_SB,      /* South Bismarck */
  RHEA_PLATE_EARTH_TI,      /* Timor */
  RHEA_PLATE_EARTH_TO,      /* Tonga */
  RHEA_PLATE_EARTH_WL,      /* Woodlark */

  RHEA_PLATE_EARTH_N        /* (number of plates) */
}
rhea_plate_earth_label_t;

#define RHEA_PLATE_EARTH_MORVEL25_BEGIN RHEA_PLATE_EARTH_AM
#define RHEA_PLATE_EARTH_MORVEL25_END   RHEA_PLATE_EARTH_YZ
#define RHEA_PLATE_EARTH_BIRD2003_BEGIN RHEA_PLATE_EARTH_AS
#define RHEA_PLATE_EARTH_BIRD2003_END   RHEA_PLATE_EARTH_WL

/******************************************************************************
 * Options & Monitoring
 *****************************************************************************/

/* options for plates */
typedef struct rhea_plate_options
{
  /* number of polygons, each of which represents one plate */
  int                 n_polygons;

  /* text files with vertices of plate polygons */
  char               *vertices_file_path_txt;
  int                 n_vertices_total;
  char               *vertices_coarse_container_file_path_txt;
  int                 n_vertices_coarse_total;

  /* storage of polygon vertices */
  float             **vertices_x;
  float             **vertices_y;
  size_t             *n_vertices;

  /* storage of vertices of coarse polygon containers */
  float             **vertices_coarse_container_x;
  float             **vertices_coarse_container_y;
  size_t             *n_vertices_coarse_container;

  /* translation of polygon & container vertices */
  float              *translation_x;
  float              *translation_y;

  /* range of values for the vertex coordinates */
  float               x_min;
  float               x_max;
  float               y_min;
  float               y_max;

  /* plate velocities (Euler poles) */
  double             *angular_velocity;

  /* cross sectional domain: plate boundaries and (tangential) velocities */
  char               *xsection_boundary_lon_list;
  float              *xsection_boundary;
  char               *xsection_tangential_velocity_mm_yr_list;
  double             *xsection_tangential_velocity;
  int                 xsection_n_intervals;
  double              xsection_shrink_factor;

  /* options (not owned) */
  rhea_domain_options_t      *domain_options;
  rhea_temperature_options_t *temp_options;
}
rhea_plate_options_t;

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_plate_add_options (ymir_options_t * opt_sup);

/**
 * Processes options and stores them.
 */
void                rhea_plate_process_options (
                                    rhea_plate_options_t *opt,
                                    rhea_domain_options_t *domain_options,
                                    rhea_temperature_options_t *temp_options);

/**
 * Initializes performance counters.
 */
void                rhea_plate_perfmon_init (const int activate,
                                             const int skip_if_active);

/**
 * Prints statistics collected by performance monitors.
 */
void                rhea_plate_perfmon_print (sc_MPI_Comm mpicomm,
                                              const int print_wtime,
                                              const int print_n_calls,
                                              const int print_flops);

/******************************************************************************
 * Plate Data
 *****************************************************************************/

/**
 * Allocates and performs the setup of data that is required for plate
 * retrieval.
 */
int                 rhea_plate_data_create (rhea_plate_options_t *opt,
                                            sc_MPI_Comm mpicomm);

/**
 * Clears storage of plate data.
 */
void                rhea_plate_data_clear (rhea_plate_options_t *opt);

/**
 * Returns the number of plates.
 */
int                 rhea_plate_get_n_plates (rhea_plate_options_t *opt);

/******************************************************************************
 * Plate Retrieval
 *****************************************************************************/

/**
 * Checks whether the given (x,y,z) coordinates are inside a specific plate.
 */
int                 rhea_plate_is_inside (const double x,
                                          const double y,
                                          const double z,
                                          const int plate_label,
                                          rhea_plate_options_t *opt);

/**
 * Sets plate labels at all entries of a vector.
 */
void                rhea_plate_set_label_vec (ymir_vec_t *vec,
                                              rhea_plate_options_t *opt);

/**
 * Sets weights of each plate as the inverse relative plate area:
 *   weight = rhea_plate_set_weight_fn_t (total_area / plate_area)
 *
 * Sets values outside of all plates to zero.
 */
typedef double    (*rhea_plate_area_to_weight_fn_t) (double plate_area,
                                                     double total_area);

void                rhea_plate_set_weight_vec (
                             ymir_vec_t *vec,
                             rhea_plate_area_to_weight_fn_t area_to_weight_fn,
                             rhea_plate_options_t *opt);

/**
 * Filters values of a vector inside a plate.  Sets values outside this plate
 * to zero.
 */
void                rhea_plate_apply_filter_vec (ymir_vec_t *vec,
                                                 const int plate_label,
                                                 rhea_plate_options_t *opt);

/**
 * Filters values of a vector inside any plate.  Sets values outside of all
 * plates to zero.
 */
void                rhea_plate_apply_filter_all_vec (ymir_vec_t *vec,
                                                     rhea_plate_options_t *opt);

/******************************************************************************
 * Plate Velocities
 *****************************************************************************/

/**
 * Generates plate velocity at all coordinates in the plate's interior from
 * a rotational axis.
 */
void                rhea_plate_velocity_generate_from_rotation (
                                                  ymir_vec_t *vel,
                                                  const double rot_axis[3],
                                                  rhea_plate_options_t *opt);

/**
 * Generates plate velocity at all coordinates in the plate's interior from
 * (Euler pole) data.  Velocity is zero outside of the plate.
 */
void                rhea_plate_velocity_generate (ymir_vec_t *vel,
                                                  const int plate_label,
                                                  rhea_plate_options_t *opt);

/**
 * Generates velocities of all plates at all coordinates from (Euler pole)
 * data.  Velocity is zero if outside of any plate.
 */
void                rhea_plate_velocity_generate_all (
                                                  ymir_vec_t *vel,
                                                  rhea_plate_options_t *opt);

/**
 * Computes the rotational axis of a plate from the given velocity `vel`.
 */
void                rhea_plate_velocity_evaluate_rotation (
                                                double rot_axis[3],
                                                ymir_vec_t *vel,
                                                const int plate_label,
                                                const int project_out_mean_rot,
                                                rhea_plate_options_t *opt);

/**
 * Computes the mean velocity of a plate from the given velocity `vel`.
 */
double              rhea_plate_velocity_get_mean_magnitude (
                                                ymir_vec_t *vel,
                                                const int plate_label,
                                                const int project_out_mean_rot,
                                                rhea_plate_options_t *opt);

/**
 * Computes mean velocities of all plates from the given velocity `vel`.
 */
void                rhea_plate_velocity_get_mean_magnitude_all (
                                                double *mean_vel_magn,
                                                ymir_vec_t *vel,
                                                const int *plate_label,
                                                const int project_out_mean_rot,
                                                rhea_plate_options_t *opt);

/**
 * Removes mean rotation (or net rotation) from a velocity field.
 */
void                rhea_plate_velocity_project_out_mean_rotation (
                                                    ymir_vec_t *vel,
                                                    rhea_plate_options_t *opt);

#endif /* RHEA_PLATE_H */
