#ifndef RHEA_DISCRETIZATION_H
#define RHEA_DISCRETIZATION_H

#include <rhea_domain.h>
#include <ymir_pressure_elem.h>

/* options for the discretization */
typedef struct rhea_discretization_options
{
  int                 order;
  int8_t              level_min;
  int8_t              level_max;

  /* boundary information */
  rhea_domain_boundary_t  *boundary;

  /* F.E. transformation function from reference space to physical space */
  mangll_X_t          X_fn;
  void               *X_data;
}
rhea_discretization_options_t;

/* global options */
extern int          rhea_discretization_level_min;

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_discretization_add_options (ymir_options_t * opt_sup);

/**
 * Processes options and stores them.
 */
void                rhea_discretization_process_options (
                                        rhea_discretization_options_t *opt,
                                        rhea_domain_options_t *domain_options);

/**
 * Sets boundary information in the discretization options object.
 */
void                rhea_discretization_options_set_boundary (
                                        rhea_discretization_options_t *opt,
                                        p4est_t *p4est,
                                        rhea_domain_options_t *domain_options);

/**
 * Clears discretization options object.
 */
void                rhea_discretization_options_clear (
                                        rhea_discretization_options_t *opt);

/**
 * Creates new p4est.
 */
p4est_t            *rhea_discretization_p4est_new (
                                        MPI_Comm mpicomm,
                                        rhea_discretization_options_t *opt,
                                        rhea_domain_options_t *domain_options);

/**
 * Creates new mangll and cnodes objects.
 */
void                rhea_discretization_mangll_and_cnodes_new (
                                        mangll_t **mangll,
                                        mangll_cnodes_t **cnodes,
                                        p4est_t *p4est,
                                        rhea_discretization_options_t *opt);

/**
 * Creates new ymir mesh and pressure element structures.
 */
void                rhea_discretization_ymir_new (
                                        ymir_mesh_t **mesh,
                                        ymir_pressure_elem_t **press_elem,
                                        mangll_t *mangll,
                                        mangll_cnodes_t *cnodes,
                                        rhea_discretization_options_t *opt);

#endif /* RHEA_DISCRETIZATION_H */
