/*
 */

#ifndef RHEA_WEAKZONE_H
#define RHEA_WEAKZONE_H

#include <ymir_vec_ops.h>

/* constant: neutral/default value for weak zone (i.e., no weakening) */
#define RHEA_WEAKZONE_NEUTRAL_VALUE (1.0)

/******************************************************************************
 * Options
 *****************************************************************************/

/* options for weak zones */
typedef struct rhea_weakzone_options
{
  /* binary/text file with (x,y,z) coordinates of weak zone points */
  char               *points_file_path_bin;
  char               *points_file_path_txt;
}
rhea_weakzone_options_t;

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_weakzone_add_options (ymir_options_t * opt_sup);

/**
 * Processes options and stores them.
 */
void                rhea_weakzone_process_options (
                                                rhea_weakzone_options_t *opt);

/******************************************************************************
 * Weak Zone Vector
 *****************************************************************************/

/**
 * Creates a new weak zone vector.
 */
ymir_vec_t         *rhea_weakzone_new (ymir_mesh_t *ymir_mesh);

/**
 * Destroys a weak zone vector.
 */
void                rhea_weakzone_destroy (ymir_vec_t *weakzone);

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_weakzone_check_vec_type (ymir_vec_t *vec);

/**
 * Checks entries of a vector.
 */
int                 rhea_weakzone_is_valid (ymir_vec_t *vec);

/******************************************************************************
 * Get & Set Functions
 *****************************************************************************/

/**
 * Gets the weak zone of one element at Gauss nodes.
 */
double             *rhea_weakzone_get_elem_gauss (sc_dmatrix_t *weak_el_mat,
                                                  ymir_vec_t *weak_vec,
                                                  const ymir_locidx_t elid);

/**
 * Sets the weak zone of one element at Gauss nodes.
 */
void                rhea_weakzone_set_elem_gauss (ymir_vec_t *weak_vec,
                                                  sc_dmatrix_t *weak_el_mat,
                                                  const ymir_locidx_t elid);

#endif /* RHEA_WEAKZONE_H */
