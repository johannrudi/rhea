/*
 */

#ifndef RHEA_COMPOSITION_H
#define RHEA_COMPOSITION_H

#include <rhea_domain.h>
#include <ymir_vec_ops.h>

/* enumerator for types of composition field */
typedef enum
{
  RHEA_COMPOSITION_NONE, /* no input for composition */
  RHEA_COMPOSITION_INT, /* use integers to define composition */
  RHEA_COMPOSITION_DENS /* continuous field such as density */
}
rhea_composition_t;

/******************************************************************************
 * Options
 *****************************************************************************/

/* options of the mantle's composition */
typedef struct rhea_composition_options
{
  /* type of composition field */
  rhea_composition_t type;

  /* data imported from file */
  char               *import_path_txt;
  char               *import_path_bin;

  /* options & properties of the computational domain */
  rhea_domain_options_t  *domain_options;

}
rhea_composition_options_t;

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_composition_add_options (ymir_options_t * opt_sup);

/**
 * Creates a new composition vector.
 */
ymir_vec_t         *rhea_composition_new (ymir_mesh_t *ymir_mesh);

/**
 * Destroys a composition vector.
 */
void                rhea_composition_destroy (ymir_vec_t *compositino);

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_composition_check_vec_type (ymir_vec_t *vec);

/**
 * Gets rank-global offsets or rank-local sizes of a distributed vector for
 * each MPI-rank.
 */
MPI_Offset         *rhea_composition_segment_offset_create (ymir_vec_t *vec);
int                 rhea_composition_segment_size_get (ymir_vec_t *vec);

/**
 * Write out composition
 */
void				rhea_composition_write (char *file_path_bin,
										ymir_vec_t *composition, sc_MPI_Comm mpicomm);

#endif /* RHEA_COMPOSITION_H */
