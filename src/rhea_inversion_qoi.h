/** RHEA_INVERSION_QOI
 *
 * Quantities of interest (QOI) of inverse problems for mantle convection.
 */

#ifndef RHEA_INVERSION_QOI_H
#define RHEA_INVERSION_QOI_H

#include <rhea_inversion_obs_stress.h>

/* types of QOI classes */
typedef enum
{
  RHEA_INVERSION_QOI_STRESS
}
rhea_inversion_qoi_class_t;

/******************************************************************************
 * Options
 *****************************************************************************/

/* QOI options regarding the inverse problem */
typedef struct rhea_inversion_qoi_options
{
  char               *stress_type_list;
  char               *stress_activate_weakzone_label_file_path_txt;
}
rhea_inversion_qoi_options_t;

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_inversion_qoi_add_options (
                              rhea_inversion_qoi_options_t *inv_qoi_options,
                              ymir_options_t *opt_sup);

/******************************************************************************
 * Quantities Of Interest (QOI)
 *****************************************************************************/

/* quantities of interest (opaque) */
typedef struct rhea_inversion_qoi rhea_inversion_qoi_t;

/**
 * Creates a new set of quantities of interest.
 */
rhea_inversion_qoi_t *rhea_inversion_qoi_new (
                            rhea_stokes_problem_t *stokes_problem,
                            rhea_inversion_qoi_options_t *inv_qoi_options);

/**
 * Destroys a set of quantities of interest.
 */
void                rhea_inversion_qoi_destroy (
                            rhea_inversion_qoi_t *inv_qoi);

/**
 * Checks if any QOI are set.
 */
int                 rhea_inversion_qoi_exists (
                            rhea_inversion_qoi_t *inv_qoi);

/******************************************************************************
 * Mapping QOI to Observation Operators
 *****************************************************************************/

void                rhea_inversion_qoi_stress_create_obs (
                                  rhea_inversion_obs_stress_t *stress_obs_type,
                                  ymir_vec_t **stress_obs_weight,
                                  rhea_inversion_qoi_t *inv_qoi,
                                  const int reduced_index);

void                rhea_inversion_qoi_stress_clear_obs (
                                  ymir_vec_t *stress_obs_weight,
                                  rhea_inversion_qoi_t *inv_qoi);

/******************************************************************************
 * QOI Vector
 *****************************************************************************/

#if 0
/**
 * Creates/destroys a QOI vector.
 */
ymir_vec_t         *rhea_inversion_qoi_vec_new (rhea_inversion_qoi_t *qoi);

void                rhea_inversion_qoi_vec_destroy (ymir_vec_t *vec);

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_inversion_qoi_vec_check_type (
                                            ymir_vec_t *vec,
                                            rhea_inversion_qoi_t *inv_qoi);

/**
 * Checks entries of a vector.
 */
int                 rhea_inversion_qoi_vec_is_valid (
                                            ymir_vec_t *vec,
                                            sc_MPI_Comm mpicomm,
                                            rhea_inversion_qoi_t *inv_qoi);

/**
 * Prints the values of a QOI vector.
 */
void                rhea_inversion_qoi_vec_print (
                                            ymir_vec_t *vec,
                                            rhea_inversion_qoi_t *inv_qoi);

/**
 * Copy entries from a reduced to a full QOI vector.
 */
void                rhea_inversion_qoi_vec_reduced_copy (
                                  ymir_vec_t *vec,
                                  sc_dmatrix_t *vec_reduced,
                                  rhea_inversion_qoi_t *inv_qoi,
                                  const rhea_inversion_qoi_class_t class_id);
#endif

/**
 * Creates/destroys a reduced QOI vector.
 */
sc_dmatrix_t       *rhea_inversion_qoi_stress_vec_reduced_new (
                                  rhea_inversion_qoi_t *inv_qoi);

void                rhea_inversion_qoi_stress_vec_reduced_destroy (
                                  sc_dmatrix_t *vec_reduced);

/**
 * Sets an entry of a QOI vector.
 * Returns index of which entry was set.
 */
int                 rhea_inversion_qoi_stress_vec_set_entry (
                                  sc_dmatrix_t *vec_reduced,
                                  const double val,
                                  rhea_inversion_qoi_t *inv_qoi,
                                  const int stress_qoi_index,
                                  const int weakzone_label);

/**
 * Creates/destroys a matrix mapping a reduced QOI vector to a reduced
 * parameter space.
 */
sc_dmatrix_t       *rhea_inversion_qoi_stress_to_param_matrix_new (
                                  const int n_parameters,
                                  rhea_inversion_qoi_t *inv_qoi);

void                rhea_inversion_qoi_stress_to_param_matrix_destroy (
                                  sc_dmatrix_t *mat);

/**
 * Sets a column of a QOI-to-parameter matrix.
 * Returns index of which column was set.
 */
int                 rhea_inversion_qoi_stress_to_param_matrix_set_column (
                                  sc_dmatrix_t *mat,
                                  sc_dmatrix_t *param_vec_reduced,
                                  rhea_inversion_qoi_t *inv_qoi,
                                  const int stress_qoi_index,
                                  const int weakzone_label);

/******************************************************************************
 * Data Access
 *****************************************************************************/

/**
 * Returns the number of QOI.
 */
int                 rhea_inversion_qoi_get_count_all (
                                    rhea_inversion_qoi_t *inv_qoi);

#endif /* RHEA_INVERSION_QOI_H */
