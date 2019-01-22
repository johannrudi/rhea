#ifndef SUBDUCTION_PHYSICS_H
#define SUBDUCTION_PHYSICS_H

#include <subduction_options.h>
#include <subduction_geometry.h>
#include <subduction_options.h>
#include <subduction_io.h>
#include <subduction_test.h>
#include <rhea.h>
#include <ymir_vec_getset.h>
#include <ymir_vec_ops.h>
#include <ymir_comm.h>


double
subd_temperature_brick_2plates_poly2 (double r, double lon,
                                       subd_options_t * subd_options);

void
subd_temperature_set_fn (double *temp, double x, double y, double z,
                          ymir_locidx_t nid, void *data);

void
subd_poly2_temperature_compute (ymir_vec_t *temperature,
                           subd_options_t *subd_options);

double
subd_temperature_cold_plate_hscm (const double radius,
                                  const double radius_max,
                                  const double radius_max_m,
                                  const double plate_age_yr,
                                  const double thermal_diffus_m2_s);

void
subd_custom_temperature_compute (ymir_vec_t *temp,
                                    subd_options_t *subd_options);

void
subd_compute_temperature (ymir_vec_t *temp,
                          rhea_temperature_options_t *temp_options,
                          subd_options_t *subd_options);

double
subd_weakzone_factor_fn (const double distance,
                          const double thickness,
                          const double thickness_const,
                          const double weak_factor);

/* Computes distance to weak zone between plates. */
double
subd_weakzone_subduct_dist_2plates_poly2 (double r, double lon,
                                         double *poly2_coeff,
                                         double start_node,
                                         double start_val,
                                         double start_deriv,
                                         double end_node,
                                         double end_val);

/* Computes distance to weak zone at ridge. */
double
subd_weakzone_ridge_dist_2plates_poly2 (double r, double lon,
                                         double end_node,
                                         double end_val);

double
subd_weakzone_brick_2plates_poly2 (double r, double lon,
                              subd_options_t * subd_options);

double
subd_weakzone_node (const double x, const double y, const double z,
                     subd_options_t *subd_options);

/* compute weak zone factor of an element */
void
subd_weakzone_elem (double *_sc_restrict weak_elem,
                     const double *x,
                     const double *y,
                     const double *z,
                     const int n_nodes_per_el,
                     subd_options_t *subd_options);

void
subd_poly2_weakzone_compute (ymir_dvec_t *weakzone, void *data);


void
subd_compute_weakzone (rhea_stokes_problem_t *stokes_problem,
                        subd_options_t *subd_options);

void
subd_layers_viscosity_old_elem (double *_sc_restrict visc_elem,
                            const double *_sc_restrict x,
                            const double *_sc_restrict y,
                            const double *_sc_restrict z,
                            const double *_sc_restrict weak_elem,
                            const int n_nodes_per_el,
                            subd_options_t *subd_options);

void
subd_custom_visc (ymir_vec_t *viscosity,
                         ymir_vec_t *rank1_tensor_scal,
                         ymir_vec_t *bounds_marker,
                         ymir_vec_t *yielding_marker,
                         ymir_vec_t *temperature,
                         ymir_vec_t *weakzone,
                         ymir_vec_t *velocity,
                         void *data);


void
subd_viscosity_set_function (rhea_stokes_problem_t *stokes_problem,
                  subd_options_t *subd_options);

void
subd_compute_rhs_vel (rhea_stokes_problem_t *stokes_problem,
                      void *data);

void
subd_set_velocity_dirichlet_bc (rhea_domain_options_t *domain_options,
                                subd_options_t *subd_options);
void
subd_compute_rhs_velbc_dirichlet (rhea_stokes_problem_t *stokes_problem,
                              subd_options_t *subd_options);

void
subd_compute_rhs_velbc_neumann (rhea_stokes_problem_t *stokes_problem,
                              subd_options_t *subd_options);


#endif
