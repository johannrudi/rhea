/* RHEA_DOMAIN_SUBSET  Describes geometric subsets of a domain. */

#ifndef RHEA_DOMAIN_SUBSET_H
#define RHEA_DOMAIN_SUBSET_H

#include <rhea_domain.h>

/******************************************************************************
 * Column Subset
 *****************************************************************************/

/** Properties of a column inside a domain. */
typedef struct rhea_domain_subset_column
{
  /* base & top of column are defined by a polygon */
  float                        *polygon_vertices_x;
  float                        *polygon_vertices_y;
  size_t                        polygon_n_vertices;
  rhea_domain_coordinate_type_t polygon_coord_type;

  /* cross sectional domain: column boundaries */
  float               xsection_boundary[2];

  /* column boundaries w.r.t. radius/depth */
  float               radius_min;
  float               radius_max;

  /* volume of column */
  float               volume;
}
rhea_domain_subset_column_t;

/**
 * Checks whether the given (x,y,z) coordinates are inside a column within the
 * domain.
 */
int                 rhea_domain_subset_point_is_in_column (
                                        const double x,
                                        const double y,
                                        const double z,
                                        rhea_domain_subset_column_t *column,
                                        rhea_domain_options_t *domain_options);

/******************************************************************************
 * Column Filter
 *****************************************************************************/

/**
 * Filters a vector according to a column subset.
 */
ymir_locidx_t       rhea_domain_subset_apply_filter (
                                      ymir_vec_t *vec,
                                      ymir_gloidx_t *n_points_in_subset_global,
                                      rhea_domain_options_t *domain_options,
                                      rhea_domain_subset_column_t *column);

/******************************************************************************
 * Column I/O
 *****************************************************************************/

/**
 * Writes coordinates & values of a vector that lie within a column subset to a
 * text file.
 */
ymir_locidx_t       rhea_domain_subset_write_txt (
                                      const char *file_path_base,
                                      ymir_gloidx_t *n_points_in_subset_global,
                                      ymir_vec_t *vec,
                                      rhea_domain_options_t *domain_options,
                                      rhea_domain_subset_column_t *column);

/******************************************************************************
 * Specific Domain Subsets
 *****************************************************************************/

/**
 * Filters a vector around Aleutian subduction.
 */
ymir_locidx_t       rhea_domain_subset_apply_filter_aleutian (
                                      ymir_vec_t *vec,
                                      ymir_gloidx_t *n_points_in_subset_global,
                                      rhea_domain_options_t *domain_options);

/**
 * Writes coordinates & values of a vector around Aleutian subduction.
 */
ymir_locidx_t       rhea_domain_subset_write_txt_aleutian (
                                      const char *file_path_base,
                                      ymir_vec_t *vec,
                                      ymir_gloidx_t *n_points_in_subset_global,
                                      rhea_domain_options_t *domain_options);

#endif /* RHEA_DOMAIN_SUBSET_H */
