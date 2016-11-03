/*
 */

#ifndef RHEA_TEMPERATURE_H
#define RHEA_TEMPERATURE_H

#include <rhea_domain.h>
#include <ymir_vec_ops.h>

/* constant: default value for temperature (gives const viscosity) */
#define RHEA_TEMPERATURE_DEFAULT_VALUE (0.5)

/* enumerator for types of temperature */
typedef enum
{
  RHEA_TEMPERATURE_NONE,
  RHEA_TEMPERATURE_COLD_PLATE,
  RHEA_TEMPERATURE_2PLATES_POLY2,
  RHEA_TEMPERATURE_IMPORT
}
rhea_temperature_t;

/* options of the mantle's temperature */
typedef struct rhea_temperature_options
{
  /* type of the temperature */
  rhea_temperature_t  type;

  /* data imported from file */
  char               *import_path_txt;
  char               *import_path_bin;

  /* cold plate model */
  double              cold_plate_model_plate_age_yr;

  /* thermal constants */
  double              thermal_diffusivity_m2_s;

  /* sinker */
  int                 sinker_active;
  double              sinker_center_x;
  double              sinker_center_y;
  double              sinker_center_z;
  double              sinker_decay;
  double              sinker_width;
  double              sinker_scaling;

  /* buoyancy right-hand side derived from temperature */
  double              rhs_scaling;

  /* options & properties of the computational domain */
  rhea_domain_options_t  *domain_options;
}
rhea_temperature_options_t;

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_temperature_add_options (ymir_options_t * opt_sup);

/**
 * Processes options and stores them.
 */
void                rhea_temperature_process_options (
                                        rhea_temperature_options_t *opt,
                                        rhea_domain_options_t *domain_options);

/**
 * Creates a new temperature vector.
 */
ymir_vec_t         *rhea_temperature_new (ymir_mesh_t *ymir_mesh);

/**
 * Destroys a temperature vector.
 */
void                rhea_temperature_destroy (ymir_vec_t *temperature);

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_temperature_check_vec_type (ymir_vec_t *vec);

/**
 * Checks entries of a vector.
 */
int                 rhea_temperature_is_valid (ymir_vec_t *vec);

/**
 * Computes the temperature.
 */
void                rhea_temperature_compute (ymir_vec_t *temperature,
                                              rhea_temperature_options_t *opt);

/**
 * Computes the background temperature.
 */
void                rhea_temperature_background_compute (
                                              ymir_vec_t *back_temperature,
                                              rhea_temperature_options_t *opt);

/**
 * Gets the temperature of one element at Gauss nodes.
 */
double             *rhea_temperature_get_elem_gauss (sc_dmatrix_t *temp_el_mat,
                                                     ymir_vec_t *temp_vec,
                                                     const ymir_locidx_t elid);

/**
 * Computes velocity right-hand side in (primal) function space, given a
 * temperature vector.
 */
void                rhea_temperature_compute_rhs_vel (
                                              ymir_vec_t *rhs_vel,
                                              ymir_vec_t *temperature,
                                              rhea_temperature_options_t *opt);

#endif /* RHEA_TEMPERATURE_H */
