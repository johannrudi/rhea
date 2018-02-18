#ifndef RHEA_AMR_H
#define RHEA_AMR_H

#include <rhea_domain.h>

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_amr_add_options (ymir_options_t * opt_sup);

/******************************************************************************
 * Initial AMR for p4est
 *****************************************************************************/

/**
 * Performs initial refinement of a p4est mesh.
 */
int                 rhea_amr_init_refine (p4est_t *p4est,
                                          int level_min,
                                          int level_max,
                                          rhea_domain_options_t
                                            *domain_options);

/******************************************************************************
 * AMR for ymir
 *****************************************************************************/

/* flags for AMR */
typedef enum
{
  RHEA_AMR_FLAG_COARSEN   = -1,
  RHEA_AMR_FLAG_NO_CHANGE =  0,
  RHEA_AMR_FLAG_REFINE    =  1
}
rhea_amr_flag_t;

/**
 * Flags elements of a p4est mesh for coarsening/refinement.
 *
 * \param [in] p4est    p4est mesh
 * \param [in/out] data User data
 * \return              Relative #elements flagged for coarsening/refinement
 */
typedef double    (*rhea_amr_flag_elements_fn_t) (p4est_t *p4est, void *data);

/**
 * Initializes data before mesh is altered through coarsening/refinement.
 *
 * \param [in] p4est    p4est mesh before AMR
 * \param [in/out] data User data
 */
typedef void      (*rhea_amr_data_initialize_fn_t) (p4est_t *p4est, void *data);

/**
 * Finalizes data after AMR is finished.
 *
 * \param [in] p4est    p4est mesh after AMR
 * \param [in/out] data User data
 */
typedef void      (*rhea_amr_data_finalize_fn_t) (p4est_t *p4est, void *data);

/**
 * Projects data to an adapted (i.e., coarsened/refined) mesh.
 *
 * \param [in] p4est    Adapted p4est mesh
 * \param [in/out] data User data
 */
typedef void      (*rhea_amr_data_project_fn_t) (p4est_t *p4est, void *data);

/**
 * Partitions data from one mesh to another.
 *
 * \param [in] p4est    Adapted & partitioned p4est mesh
 * \param [in/out] data User data
 */
typedef void      (*rhea_amr_data_partition_fn_t) (p4est_t *p4est, void *data);

/**
 * Performs adaptive coarsening and refinement of a p4est mesh.
 */
int                 rhea_amr (p4est_t *p4est,
                              const double n_flagged_elements_tol,
                              const int amr_recursive_count,
                              const double n_flagged_elements_recursive_tol,
                              rhea_amr_flag_elements_fn_t flag_elements_fn,
                              void *flag_elements_data,
                              rhea_amr_data_initialize_fn_t data_initialize_fn,
                              rhea_amr_data_finalize_fn_t data_finalize_fn,
                              rhea_amr_data_project_fn_t data_project_fn,
                              rhea_amr_data_partition_fn_t data_partition_fn,
                              void *data);

/******************************************************************************
 * Generic Callback Functions for Coarsening/Refinement
 *****************************************************************************/

int                 rhea_amr_coarsen_all_fn (p4est_t * p4est,
                                             p4est_topidx_t which_tree,
                                             p4est_quadrant_t * quadrant);
int                 rhea_amr_refine_all_fn (p4est_t * p4est,
                                            p4est_topidx_t which_tree,
                                            p4est_quadrant_t * quadrant);

int                 rhea_amr_refine_to_min_level_fn (
                                                  p4est_t * p4est,
                                                  p4est_topidx_t which_tree,
                                                  p4est_quadrant_t * quadrant);

int                 rhea_amr_coarsen_half_fn (p4est_t * p4est,
                                              p4est_topidx_t which_tree,
                                              p4est_quadrant_t * quadrant);
int                 rhea_amr_refine_half_fn (p4est_t * p4est,
                                             p4est_topidx_t which_tree,
                                             p4est_quadrant_t * quadrant);

typedef struct rhea_amr_refine_depth_data
{
  int                *depth;
  int                 count;
  int                 level_min;
  double              p4est_zmax;
}
rhea_amr_refine_depth_data_t;

int                 rhea_amr_refine_depth_fn (p4est_t * p4est,
                                              p4est_topidx_t which_tree,
                                              p4est_quadrant_t * quadrant);

int                 rhea_amr_refine_depth_box_fn (p4est_t * p4est,
                                                  p4est_topidx_t which_tree,
                                                  p4est_quadrant_t * quadrant);

/******************************************************************************
 * Generic Flagging for Coarsening/Refinement
 *****************************************************************************/

/**
 * Sums local contributions to get the global number of flagged quadrants.
 * Returns the relative number of flagged quadrants.
 */
double              rhea_amr_get_relative_global_num_flagged (
                                            const p4est_locidx_t n_flagged_loc,
                                            p4est_t *p4est);

double              rhea_amr_flag_coarsen_half_fn (p4est_t *p4est, void *data);
double              rhea_amr_flag_refine_half_fn (p4est_t *p4est, void *data);

#endif /* RHEA_AMR_H */
