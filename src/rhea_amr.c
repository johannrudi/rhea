/*
 */

#include <rhea_amr.h>

const char         *rhea_amr_p4est_refine_name[RHEA_AMR_P4EST_REFINE_N] =
{
  "uniform",
  "half",
  "depth"
};

rhea_amr_p4est_refine_t
rhea_amr_p4est_refine_decode (const char * name)
{
  const int           n_types = (int) RHEA_AMR_P4EST_REFINE_N;
  int                 typeid;
  rhea_amr_p4est_refine_t type;

  /* exit if nothing to do */
  if (name == NULL) {
    return RHEA_AMR_P4EST_REFINE_NONE;
  }

  /* match name to a type */
  type = RHEA_AMR_P4EST_REFINE_UNKNOWN;
  for (typeid = 0; typeid < n_types; typeid++) {
    if (!strcmp (name, rhea_amr_p4est_refine_name[typeid])) {
      type = (rhea_amr_p4est_refine_t) typeid;
    }
  }

  /* return refinement type */
  return type;
}

int
rhea_amr_p4est_refine_half (p4est_t * p4est, p4est_topidx_t which_tree,
                            p4est_quadrant_t * quadrant)
{
  if (P4EST_ROOT_LEN / 2 <= quadrant->z) {
    return 1;
  }
  else {
    return 0;
  }
}

int
rhea_amr_p4est_refine_depth (p4est_t * p4est, p4est_topidx_t which_tree,
                              p4est_quadrant_t * quadrant)
{
  rhea_amr_p4est_refine_depth_data_t  *d = p4est->user_pointer;
  const int          *depth = d->depth;
  const int           quad_depth = P4EST_ROOT_LEN - quadrant->z;
  int                 k;

  if ((int) quadrant->level < d->level_min) {
    return 1;
  }

  for (k = 0; k < d->count; k++) {
    if (quad_depth <= depth[k] && (int) quadrant->level < d->level_min+k+1) {
      return 1;
    }
  }

  return 0;
}
