/* RHEA_STRESS  Stress tensor and its second invariant. */

#ifndef RHEA_STRESS_H
#define RHEA_STRESS_H

#include <rhea_domain.h>
#include <rhea_temperature.h>
#include <rhea_viscosity.h>

/******************************************************************************
 * Stress Vector
 *****************************************************************************/

/**
 * Creates a new (second invariant of) stress tensor.
 */
ymir_vec_t         *rhea_stress_new (ymir_mesh_t *ymir_mesh);

ymir_vec_t         *rhea_stress_2inv_new (ymir_mesh_t *ymir_mesh);

/**
 * Destroys a (second invariant of) stress tensor.
 */
void                rhea_stress_destroy (ymir_vec_t *stress);

void                rhea_stress_2inv_destroy (ymir_vec_t *stress_2inv);

/**
 * Converts entries of a nondimensional stress vector into dimensional
 * values:
 *
 *   [Pa]
 */
void                rhea_stress_convert_to_dimensional_Pa (
                                      ymir_vec_t * stress,
                                      rhea_domain_options_t *domain_options,
                                      rhea_temperature_options_t *temp_options,
                                      rhea_viscosity_options_t *visc_options);

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_stress_check_vec_type (ymir_vec_t *vec);

int                 rhea_stress_2inv_check_vec_type (ymir_vec_t *vec);

/**
 * Checks entries of a vector.
 */
int                 rhea_stress_is_valid (ymir_vec_t *vec);

int                 rhea_stress_2inv_is_valid (ymir_vec_t *vec);

/**
 * Computes the viscous stress tensor.
 */
void                rhea_stress_compute_viscstress (ymir_vec_t *viscstress,
                                                    ymir_vec_t *strainrate,
                                                    ymir_vec_t *viscosity);

/**
 * Computes the square root of the second invariant of the viscous stress
 * tensor.
 */
void                rhea_stress_compute_viscstress_sqrt_of_2inv (
                                             ymir_vec_t *viscstress_sqrt_2inv,
                                             ymir_vec_t *strainrate_sqrt_2inv,
                                             ymir_vec_t *viscosity);

/**
 * Computes the norm of a stress tensor.
 */
double              rhea_stress_compute_norm (ymir_vec_t *stress);

/**
 * Separates the 6 components of a (symmetric) stress tensor into diagonal and
 * off-diagonal components.
 */
void                rhea_stress_separate_diag_offdiag (
                                                    ymir_vec_t *stress_diag,
                                                    ymir_vec_t *stress_offdiag,
                                                    ymir_vec_t *stress);

/******************************************************************************
 * Stress Surface Vector
 *****************************************************************************/

/**
 * Creates a new stress vector at surface.
 */
ymir_vec_t         *rhea_stress_surface_new (ymir_mesh_t *ymir_mesh);

/**
 * Destroys a stress vector at surface.
 */
void                rhea_stress_surface_destroy (ymir_vec_t *vel_surf);

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_stress_surface_check_vec_type (ymir_vec_t *vec);

/**
 * Checks entries of a vector.
 */
int                 rhea_stress_surface_is_valid (ymir_vec_t *vec);

/**
 * Extracts the normal component of the stress at the surface from the residual
 * of the Stokes momentum equation.
 */
void                rhea_stress_surface_extract_from_residual (
                                                ymir_vec_t *stress_norm_surf,
                                                ymir_vec_t *residual_mom);

/******************************************************************************
 * Statistics
 *****************************************************************************/

/**
 * Computes global strain rate statistics.
 */
void                rhea_stress_stats_get_global (
                              double *min_Pa, double *max_Pa, double *mean_Pa,
                              ymir_vec_t *velocity, ymir_vec_t *viscosity,
                              rhea_domain_options_t *domain_options,
                              rhea_temperature_options_t *temp_options,
                              rhea_viscosity_options_t *visc_options);

void                rhea_stress_surface_stats_get_global (
                              double *min_Pa, double *max_Pa, double *mean_Pa,
                              ymir_vec_t *stress_norm_surf,
                              rhea_domain_options_t *domain_options,
                              rhea_temperature_options_t *temp_options,
                              rhea_viscosity_options_t *visc_options);

#endif /* RHEA_STRESS_H */
