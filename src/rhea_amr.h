#ifndef RHEA_AMR_H
#define RHEA_AMR_H

#include <rhea_domain.h>

int                 rhea_amr_p4est_refine_half (p4est_t * p4est,
                                                p4est_topidx_t which_tree,
                                                p4est_quadrant_t * quadrant);

#endif /* RHEA_AMR_H */
