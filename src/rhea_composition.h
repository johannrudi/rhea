/*
 */

#ifndef RHEA_COMPOSITION_H
#define RHEA_COMPOSITION_H

#include <rhea_domain.h>
#include <ymir_vec_ops.h>

#define RHEA_COMPOSITION_MAXIMUM_FLAVOR (100)

/* enumerator for types of composition field */
typedef enum
{
  RHEA_COMPOSITION_NONE, /* no input for composition */
  RHEA_COMPOSITION_FLAVOR, /* use integers to define composition */
  RHEA_COMPOSITION_DENSITY /* continuous field such as density */
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

  /* number of flavor */
  int				 n_flavor;

  /* density anomaly for each flavor */
  double			 density_flavor[RHEA_COMPOSITION_MAXIMUM_FLAVOR];
  /* viscosity effect for each flavor */
  double			 visc_flavor[RHEA_COMPOSITION_MAXIMUM_FLAVOR];

  /* data imported from file */
  char               *import_path_txt;
  char               *import_path_bin;

  /* buoyancy right-hand side derived from composition */
  double              rhs_scaling_comp;

  /* options & properties of the computational domain */
  rhea_domain_options_t  *domain_options;

}
rhea_composition_options_t;

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_composition_add_options (ymir_options_t * opt_sup);

/**
 * Processes options and stores them.
 */
void                rhea_composition_process_options (
                                        rhea_composition_options_t *opt,
                                        rhea_domain_options_t *domain_options);

/**
 * Returns wheter a composition exists.
 */
int                 rhea_composition_exists (rhea_composition_options_t *opt);

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
MPI_Offset          rhea_composition_segment_offset_get (ymir_vec_t *vec);
int                 rhea_composition_segment_size_get (ymir_vec_t *vec);

/**
 * convert composition to density and viscosity
 */
void rhea_composition_convert (ymir_vec_t *composition,
							ymir_vec_t **comp_density, ymir_vec_t **comp_visc,
							ymir_mesh_t *ymir_mesh,
							rhea_composition_options_t *opt);

/**
 * Read composition field
 */
void				rhea_composition_read (ymir_vec_t *composition,
										rhea_composition_options_t *opt);

/**
 * Write out composition
 */
int				rhea_composition_write (char *file_path_bin,
										ymir_vec_t *composition, sc_MPI_Comm mpicomm);
/**
 * Computes and adds velocity right-hand side in (primal) function space, given
 * a compositional density vector.
 */
void rhea_composition_add_rhs_vel (ymir_vec_t *rhs_vel,
                                   ymir_vec_t *comp_density,
                                   rhea_composition_options_t *opt);

#endif /* RHEA_COMPOSITION_H */
