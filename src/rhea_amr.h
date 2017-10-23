#ifndef RHEA_AMR_H
#define RHEA_AMR_H

#include <rhea_domain.h>

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_amr_add_options (ymir_options_t * opt_sup);

/**
 * Performs initial refinement of a p4est mesh.
 */
int                 rhea_amr_init_refine (p4est_t *p4est,
                                          const int level_min,
                                          const int level_max,
                                          rhea_domain_options_t
                                            *domain_options);

/**
 * Marks elements of a p4est mesh for coarsening/refinement.
 *
 * \param [in] p4est    p4est mesh
 * \param [in/out] data User data
 * \return              Relative #elements marked for coarsening/refinement
 */
typedef double    (*rhea_amr_mark_elements_fn_t) (p4est_t *p4est, void *data);

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
                              const double n_marked_elements_tol,
                              const int amr_recursive_count,
                              const double n_marked_elements_recursive_tol,
                              rhea_amr_mark_elements_fn_t mark_elements_fn,
                              void *mark_elements_data,
                              rhea_amr_data_initialize_fn_t data_initialize_fn,
                              rhea_amr_data_finalize_fn_t data_finalize_fn,
                              rhea_amr_data_project_fn_t data_project_fn,
                              rhea_amr_data_partition_fn_t data_partition_fn,
                              void *data);

#endif /* RHEA_AMR_H */
