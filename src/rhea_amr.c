/*
 */

#include <rhea_amr.h>

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
