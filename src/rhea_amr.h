#ifndef RHEA_AMR_H
#define RHEA_AMR_H

#include <rhea_domain.h>

/* types of refinement when creating a new p4est mesh */
typedef enum
{
  RHEA_AMR_P4EST_REFINE_UNKNOWN = -2,
  RHEA_AMR_P4EST_REFINE_NONE    = -1,
  RHEA_AMR_P4EST_REFINE_UNIFORM =  0, /* default, does not have a refine fn. */
  RHEA_AMR_P4EST_REFINE_HALF,
  RHEA_AMR_P4EST_REFINE_DEPTH,
  RHEA_AMR_P4EST_REFINE_N
}
rhea_amr_p4est_refine_t;

rhea_amr_p4est_refine_t rhea_amr_p4est_refine_decode (const char * name);

int                 rhea_amr_p4est_refine_half (p4est_t * p4est,
                                                p4est_topidx_t which_tree,
                                                p4est_quadrant_t * quadrant);

typedef struct rhea_amr_p4est_refine_depth_data
{
  int                *depth;
  int                 count;
  int                 level_min;
}
rhea_amr_p4est_refine_depth_data_t;

int                 rhea_amr_p4est_refine_depth (p4est_t * p4est,
                                                 p4est_topidx_t which_tree,
                                                 p4est_quadrant_t * quadrant);

#endif /* RHEA_AMR_H */
