#ifndef RHEA_DISCRETIZATION_H
#define RHEA_DISCRETIZATION_H

#include <rhea_domain.h>
#include <ymir_pressure_elem.h>

/* options for the discretization */
typedef struct rhea_discretization_options
{
  int                 order;
  int                 level_min;
  int                 level_max;

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
 * Sets an X-function for geometry transformation of the mesh.
 */
void                rhea_discretization_set_user_X_fn (
                                        rhea_discretization_options_t *opt,
                                        mangll_X_t X_fn, void *X_data);

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

/******************************************************************************
 * Constructor/Destructor for p4est
 *****************************************************************************/

p4est_t            *rhea_discretization_p4est_new (
                                        MPI_Comm mpicomm,
                                        rhea_discretization_options_t *opt,
                                        rhea_domain_options_t *domain_options);

void                rhea_discretization_p4est_destroy (p4est_t *p4est);

/* user data of a p4est quadrant */
typedef struct
{
  int                 amr_flag;
}
rhea_p4est_quadrant_data_t;

/**
 * Callback function for initializing new p4est quadrants.
 */
void                rhea_p4est_init_fn (p4est_t *p4est, p4est_topidx_t tree,
                                        p4est_quadrant_t *quadrant);

/******************************************************************************
 * Constructor/Destructor for mangll
 *****************************************************************************/

void                rhea_discretization_mangll_continuous_new (
                                        mangll_t **mangll,
                                        mangll_cnodes_t **cnodes,
                                        p4est_t *p4est,
                                        rhea_discretization_options_t *opt);

void                rhea_discretization_mangll_continuous_destroy (
                                        mangll_t *mangll,
                                        mangll_cnodes_t *cnodes);

mangll_t           *rhea_discretization_mangll_discontinuous_new (
                                        p4est_t *p4est, const int order,
                                        mangll_X_t X_fn, void *X_data);

void                rhea_discretization_mangll_discontinuous_destroy (
                                        mangll_t *mangll);

mangll_t           *rhea_discretization_mangll_interpolation_new (
                                        p4est_t *p4est, const int order);

void                rhea_discretization_mangll_interpolation_destroy (
                                        mangll_t *mangll);

/******************************************************************************
 * Constructor/Destructor for ymir
 *****************************************************************************/

void                rhea_discretization_ymir_mesh_new_from_mangll (
                                          ymir_mesh_t **ymir_mesh,
                                          ymir_pressure_elem_t **press_elem,
                                          mangll_t *mangll,
                                          mangll_cnodes_t *cnodes,
                                          rhea_discretization_options_t *opt);

void                rhea_discretization_ymir_mesh_new_from_p4est (
                                          ymir_mesh_t **ymir_mesh,
                                          ymir_pressure_elem_t **press_elem,
                                          p4est_t *p4est,
                                          rhea_discretization_options_t *opt);

void                rhea_discretization_ymir_mesh_destroy (
                                          ymir_mesh_t *ymir_mesh,
                                          ymir_pressure_elem_t *press_elem);

#endif /* RHEA_DISCRETIZATION_H */
