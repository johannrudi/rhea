/** RHEA_VISCOSITY
 *
 * Viscosity modeling mantle convection.
 */

#ifndef RHEA_VISCOSITY_H
#define RHEA_VISCOSITY_H

#include <rhea_domain.h>
#include <ymir_vec_ops.h>

/* constant: neutral/default value for viscosity */
#define RHEA_VISCOSITY_NEUTRAL_VALUE (1.0)

/* types of viscosities */
typedef enum
{
  RHEA_VISCOSITY_LINEAR,    /* linear rheology */
  RHEA_VISCOSITY_NONLINEAR  /* nonlinear rheology */
}
rhea_viscosity_t;

/* types of linear components in the viscosity */
typedef enum
{
  RHEA_VISCOSITY_LINEAR_CONST,        /* constant linear viscosity */
  RHEA_VISCOSITY_LINEAR_TEMPREVERSE,  /* reverse tempererature: (1 - T) */
  RHEA_VISCOSITY_LINEAR_ARRHENIUS     /* Arrhenius relationship */
}
rhea_viscosity_linear_t;

/* types of nonlinear components in the viscosity */
typedef enum
{
  RHEA_VISCOSITY_NONLINEAR_SRW,     /* strain rate weakening */
  RHEA_VISCOSITY_NONLINEAR_YLD,     /* plastic yielding */
  RHEA_VISCOSITY_NONLINEAR_SRW_YLD  /* strain rate weakening & yielding */
}
rhea_viscosity_nonlinear_t;

/* types of initial viscosities for nonlinear solvers */
typedef enum
{
  RHEA_VISCOSITY_NONLINEAR_INIT_DEFAULT,
  RHEA_VISCOSITY_NONLINEAR_INIT_LIN,
  RHEA_VISCOSITY_NONLINEAR_INIT_LIN_RESCALE_UM
}
rhea_viscosity_nonlinear_init_t;

/* types of viscosity models that determine how the components are combined */
typedef enum
{
  /* (1) upper bound, (2) weak zone, (3) yielding, (4) lower bound;
   * viscosity bounds via cut-off */
  RHEA_VISCOSITY_MODEL_UWYL,

  /* (1) upper bound, (2) weak zone, (3) yielding, (4) lower bound;
   * upper viscosity bound via cut-off, lower viscosity bound via addition */
  RHEA_VISCOSITY_MODEL_UWYL_LADD_UCUT,

  /* (1) upper bound + shift, (2) weak zone, (3) yielding, (4) lower bound;
   * upper viscosity bound via shift, lower viscosity bound via addition */
  RHEA_VISCOSITY_MODEL_UWYL_LADD_USHIFT
}
rhea_viscosity_model_t;

/******************************************************************************
 * Options
 *****************************************************************************/

/* options of the mantle's viscosity */
typedef struct rhea_viscosity_options
{
  /* types of linear/nonlinear viscosity components and its composition */
  rhea_viscosity_t                 type;
  rhea_viscosity_linear_t          type_linear;
  rhea_viscosity_nonlinear_t       type_nonlinear;
  rhea_viscosity_nonlinear_init_t  type_nonlinear_init;
  rhea_viscosity_model_t           model;

  /* viscosity constants */
  double              representative_Pas;

  /* lower and upper bounds for the viscosity */
  double              min;
  double              max;

  /* scaling factors */
  double              upper_mantle_scaling;
  double              lower_mantle_scaling;

  /* activation energy in Arrhenius relationship */
  double              upper_mantle_arrhenius_activation_energy;
  double              lower_mantle_arrhenius_activation_energy;

  /* stress exponent that governs strain rate weakening (aka. `n`) */
  double              stress_exponent;

  /* parameters for plastic yielding */
  double              yield_strength;

  /* regularization for projector of nonlinear viscosity */
  double              nonlinear_projector_regularization;

  /* options & properties of the computational domain */
  rhea_domain_options_t  *domain_options;
}
rhea_viscosity_options_t;

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_viscosity_add_options (ymir_options_t * opt_sup);

/**
 * Processes options and stores them.
 */
void                rhea_viscosity_process_options (
                                        rhea_viscosity_options_t *opt,
                                        rhea_domain_options_t *domain_options);

/******************************************************************************
 * Viscosity Vector
 *****************************************************************************/

/**
 * Creates a new viscosity vector.
 */
ymir_vec_t         *rhea_viscosity_new (ymir_mesh_t *ymir_mesh);

/**
 * Destroys a viscosity vector.
 */
void                rhea_viscosity_destroy (ymir_vec_t *viscosity);

/**
 * Converts entries of a nondimensional viscosity vector into dimensional
 * values:
 *
 *   [Pa*s]
 */
void                rhea_viscosity_convert_to_dimensional_Pas (
                                                ymir_vec_t * viscosity,
                                                rhea_viscosity_options_t *opt);

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_viscosity_check_vec_type (ymir_vec_t *vec);

/**
 * Checks entries of a vector.
 */
int                 rhea_viscosity_is_valid (ymir_vec_t *vec);

/******************************************************************************
 * Viscosity Surface Vector
 *****************************************************************************/

/**
 * Creates a new viscosity vector at surface.
 */
ymir_vec_t         *rhea_viscosity_surface_new (ymir_mesh_t *ymir_mesh);

/**
 * Destroys a viscosity vector at surface.
 */
void                rhea_viscosity_surface_destroy (ymir_vec_t *visc_surf);

/**
 * Checks whether a vector is of the right type.
 */
int                 rhea_viscosity_surface_check_vec_type (ymir_vec_t *vec);

/**
 * Checks entries of a vector.
 */
int                 rhea_viscosity_surface_is_valid (ymir_vec_t *vec);

/**
 * Creates a new viscosity vector at surface with values interpolated from a
 * viscosity volume vector.
 */
ymir_vec_t         *rhea_viscosity_surface_new_from_vol (ymir_vec_t *visc_vol);

/**
 * Interpolates viscosity from volume to surface.
 */
void                rhea_viscosity_surface_interpolate (
                                                  ymir_vec_t *visc_surf,
                                                  ymir_vec_t *visc_vol,
                                                  const double visc_surf_min,
                                                  const double visc_surf_max);

/******************************************************************************
 * Viscosity Computation
 *****************************************************************************/

/**
 * Computes a custom viscosity.  Callback function for viscosity computation.
 *
 * \param [out] viscosity       Viscosity (Gauss quadrature nodes)
 * \param [out] proj_scal       Scaling for rank-1 fourth-order tensor
 *                                (Gauss quadrature nodes, may be NULL)
 * \param [out] marker          Marker for different physics in the effective
 *                              viscosity
 *                                (Gauss quadrature nodes, may be NULL)
 * \param [in] temperature      Temperature (GLL nodes, may be NULL)
 * \param [in] weakzone         Weak zone factor (Gauss nodes, may be NULL)
 * \param [in] velocity         Velocity (GLL nodes, may be NULL)
 * \param [in] opt              Viscosity options
 * \param [in] data             User data
 */
typedef void      (*rhea_viscosity_compute_fn_t) (ymir_vec_t *viscosity,
                                                  ymir_vec_t *proj_scal,
                                                  ymir_vec_t *marker,
                                                  ymir_vec_t *temperature,
                                                  ymir_vec_t *weakzone,
                                                  ymir_vec_t *velocity,
                                                  void *data);

/**
 * Computes the viscosity.
 *
 * \param [out] viscosity       Viscosity (Gauss quadrature nodes)
 * \param [out] proj_scal       Scaling for rank-1 fourth-order tensor
 *                                (Gauss quadrature nodes, may be NULL)
 * \param [out] marker          Marker for different physics in the effective
 *                              viscosity
 *                                (Gauss quadrature nodes, may be NULL)
 * \param [in] temperature      Temperature (GLL nodes, may be NULL)
 * \param [in] weakzone         Weak zone factor (Gauss nodes, may be NULL)
 * \param [in] velocity         Velocity (GLL nodes, may be NULL)
 * \param [in] opt              Viscosity options
 * \param [in] data             User data
 */
void                rhea_viscosity_compute (ymir_vec_t *viscosity,
                                            ymir_vec_t *proj_scal,
                                            ymir_vec_t *marker,
                                            ymir_vec_t *temperature,
                                            ymir_vec_t *weakzone,
                                            ymir_vec_t *velocity,
                                            void *data);

/**
 * Computes the viscosity to initialize a nonlinear solver with zero velocity.
 *
 * \param [out] viscosity       Viscosity (Gauss quadrature nodes)
 * \param [out] proj_scal       Scaling for rank-1 fourth-order tensor
 *                                (Gauss quadrature nodes, may be NULL)
 * \param [out] marker          Marker for different physics in the effective
 *                              viscosity
 *                                (Gauss quadrature nodes, may be NULL)
 * \param [in] temperature      Temperature (GLL nodes, may be NULL)
 * \param [in] weakzone         Weak zone factor (Gauss nodes, may be NULL)
 * \param [in] opt              Viscosity options
 */
void                rhea_viscosity_compute_nonlinear_init (
                                               ymir_vec_t *viscosity,
                                               ymir_vec_t *proj_scal,
                                               ymir_vec_t *marker,
                                               ymir_vec_t *temperature,
                                               ymir_vec_t *weakzone,
                                               void *data);

/**
 * Computes the viscosity of one element.
 */
void                rhea_viscosity_compute_elem (
                                double *_sc_restrict visc_elem,
                                double *_sc_restrict proj_scal_elem,
                                double *_sc_restrict marker_elem,
                                const double *_sc_restrict temp_elem,
                                const double *_sc_restrict weak_elem,
                                const double *_sc_restrict strt_sqrt_2inv_elem,
                                const double *_sc_restrict x,
                                const double *_sc_restrict y,
                                const double *_sc_restrict z,
                                const int n_nodes,
                                const int *_sc_restrict Vmask,
                                rhea_viscosity_options_t *opt);

void                rhea_viscosity_compute_nonlinear_init_elem (
                                double *_sc_restrict visc_elem,
                                double *_sc_restrict proj_scal_elem,
                                double *_sc_restrict marker_elem,
                                const double *_sc_restrict temp_elem,
                                const double *_sc_restrict weak_elem,
                                const double *_sc_restrict x,
                                const double *_sc_restrict y,
                                const double *_sc_restrict z,
                                const int n_nodes,
                                const int *_sc_restrict Vmask,
                                rhea_viscosity_options_t *opt);

/******************************************************************************
 * Properties
 *****************************************************************************/

/**
 * Returns whether restriction to min/max viscosity bound is active.
 */
int                 rhea_viscosity_restrict_min (rhea_viscosity_options_t *opt);
int                 rhea_viscosity_restrict_max (rhea_viscosity_options_t *opt);

/**
 * Gets the scaling factor.
 */
double              rhea_viscosity_get_scaling (rhea_viscosity_options_t *opt,
                                                const int is_in_upper_mantle,
                                                const int restrict_min);

/**
 * Returns whether temperature dependence via Arrhenius relationship is enabled.
 */
int                 rhea_viscosity_has_arrhenius (
                                                rhea_viscosity_options_t *opt);

/**
 * Returns whether strain rate weakening physics are enabled.
 */
int                 rhea_viscosity_has_strain_rate_weakening (
                                                rhea_viscosity_options_t *opt);

/**
 * Gets the exponent for the sqrt of the 2nd invariant of the strain rate
 * tensor: `1/n`, where `n` is the stress exponent.
 */
double              rhea_viscosity_get_strain_rate_weakening_exp (
                                                rhea_viscosity_options_t *opt);

/**
 * Returns whether yielding physics are enabled.
 */
int                 rhea_viscosity_has_yielding (rhea_viscosity_options_t *opt);

/**
 * Gets the yield strength.
 */
double              rhea_viscosity_get_yield_strength (
                                                rhea_viscosity_options_t *opt);

/**
 * Returns whether regularization for the projector of the nonlinear viscosity
 * is enabled.
 */
int                 rhea_viscosity_has_nonlinear_projector_regularization (
                                                rhea_viscosity_options_t *opt);

/**
 * Gets the regularization value for the projector of the nonlinear viscosity.
 */
double              rhea_viscosity_get_nonlinear_projector_regularization (
                                                rhea_viscosity_options_t *opt);

/**
 * Returns the value by which the whole viscosity is shifted.
 */
double              rhea_viscosity_get_visc_shift (
                                               rhea_viscosity_options_t *opt);

/**
 * Returns the value by which the projection part of the (linearized) viscosity
 * coefficient is shifted.
 */
double              rhea_viscosity_get_visc_shift_proj (
                                               rhea_viscosity_options_t *opt);

/******************************************************************************
 * Get & Set Values
 *****************************************************************************/

/**
 * Gets the viscosity of one element at Gauss nodes.
 */
double             *rhea_viscosity_get_elem_gauss (sc_dmatrix_t *visc_el_mat,
                                                   ymir_vec_t *visc_vec,
                                                   const ymir_locidx_t elid);

/**
 * Sets the viscosity of one element at Gauss nodes.
 */
void                rhea_viscosity_set_elem_gauss (ymir_vec_t *visc_vec,
                                                   sc_dmatrix_t *visc_el_mat,
                                                   const ymir_locidx_t elid);

/**
 * Gets the rank-1 4th-order tensor of one element at Gauss nodes.
 */
double             *rhea_viscosity_proj_scal_get_elem_gauss (
                                                sc_dmatrix_t *proj_scal_el_mat,
                                                ymir_vec_t *proj_scal_vec,
                                                const ymir_locidx_t elid);

/**
 * Sets the rank-1 4th-order tensor of one element at Gauss nodes.
 */
void                rhea_viscosity_proj_scal_set_elem_gauss (
                                                ymir_vec_t *proj_scal_vec,
                                                sc_dmatrix_t *proj_scal_el_mat,
                                                const ymir_locidx_t elid);

/**
 * Gets a marker of one element at Gauss nodes.
 */
double             *rhea_viscosity_marker_get_elem_gauss (
                                                  sc_dmatrix_t *marker_el_mat,
                                                  ymir_vec_t *marker_vec,
                                                  const ymir_locidx_t elid);

/**
 * Sets a marker of one element at Gauss nodes.
 */
void                rhea_viscosity_marker_set_elem_gauss (
                                                  ymir_vec_t *marker_vec,
                                                  sc_dmatrix_t *marker_el_mat,
                                                  const ymir_locidx_t elid);

/******************************************************************************
 * Markers
 *****************************************************************************/

/**
 * Filters a vector where the min/max viscosity bound is active, otherwise sets
 * the values at nodes away from min/max bound to zero.
 */
void                rhea_viscosity_marker_filter_min (ymir_vec_t *vec,
                                                      ymir_vec_t *marker,
                                                      const int invert_filter);

void                rhea_viscosity_marker_filter_max (ymir_vec_t *vec,
                                                      ymir_vec_t *marker,
                                                      const int invert_filter);

/**
 * Filters a vector where yielding occurs, and otherwise sets the values at
 * nodes without yielding to zero.
 */
void                rhea_viscosity_marker_filter_yielding (
                                                   ymir_vec_t *vec,
                                                   ymir_vec_t *marker,
                                                   const int invert_filter);

/**
 * Computes the total volume where each marker is active.
 */
void                rhea_viscosity_marker_get_volume (
                                                  double *vol_min,
                                                  double *vol_max,
                                                  double *vol_yielding,
                                                  ymir_vec_t *marker);

/**
 * Computes the volume where bounds are active.
 */
void                rhea_viscosity_marker_get_bounds_volume (
                                                  double *vol_min,
                                                  double *vol_max,
                                                  ymir_vec_t *marker);

/**
 * Computes the volume where yielding is active.
 */
double              rhea_viscosity_marker_get_yielding_volume (
                                                  ymir_vec_t *yielding_marker);

/******************************************************************************
 * Filters
 *****************************************************************************/

/**
 * Computes the volume of a filter.  A filter is understood as a vector with
 * ones where the filter is active and zeros otherwise.
 */
double              rhea_viscosity_filter_compute_volume (ymir_vec_t *filter);

/**
 * Sets a filter for the upper/lower mantle.
 */
void                rhea_viscosity_stats_filter_upper_mantle (
                                        ymir_vec_t *filter,
                                        rhea_domain_options_t *domain_options);

void                rhea_viscosity_stats_filter_lower_mantle (
                                        ymir_vec_t *filter,
                                        rhea_domain_options_t *domain_options);

/**
 * Sets a filter that is active at high (i.e., lithospheric) viscosity values.
 */
void                rhea_viscosity_stats_filter_lithosphere (
                                                ymir_vec_t *filter,
                                                ymir_vec_t *viscosity,
                                                rhea_viscosity_options_t *opt,
                                                double threshold);

void                rhea_viscosity_stats_filter_lithosphere_surf (
                                                ymir_vec_t *filter_surf,
                                                ymir_vec_t *viscosity_vol,
                                                rhea_viscosity_options_t *opt,
                                                double threshold);

/**
 * Sets a filter in the upper mantle that is active away from the lithospere.
 */
void                rhea_viscosity_stats_filter_asthenosphere (
                                                ymir_vec_t *filter,
                                                ymir_vec_t *viscosity,
                                                rhea_viscosity_options_t *opt,
                                                double threshold);

/******************************************************************************
 * Statistics
 *****************************************************************************/

/**
 * Computes global viscosity statistics.
 */
void                rhea_viscosity_stats_get_global (
                                                double *min_Pas,
                                                double *max_Pas,
                                                double *mean_Pas,
                                                ymir_vec_t *viscosity,
                                                rhea_viscosity_options_t *opt);

void                rhea_viscosity_stats_get_regional (
                                                double *upper_mantle_mean_Pas,
                                                double *lower_mantle_mean_Pas,
                                                double *lith_mean_Pas,
                                                double *asth_mean_Pas,
                                                ymir_vec_t *viscosity,
                                                rhea_viscosity_options_t *opt);

/**
 * Computes the (approx.) volume of the lithosphere.
 */
double              rhea_viscosity_stats_get_lithosphere_volume (
                                                ymir_vec_t *viscosity,
                                                rhea_viscosity_options_t *opt);

/**
 * Computes the (approx.) volume of the asthenosphere.
 */
double              rhea_viscosity_stats_get_asthenosphere_volume (
                                                ymir_vec_t *viscosity,
                                                rhea_viscosity_options_t *opt);

#endif /* RHEA_VISCOSITY_H */
