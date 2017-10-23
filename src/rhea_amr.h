#ifndef RHEA_AMR_H
#define RHEA_AMR_H

#include <rhea_domain.h>

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_amr_add_options (ymir_options_t * opt_sup);

/**
 * Performs initial refinement for a p4est mesh.
 */
int                 rhea_amr_p4est_refine (p4est_t *p4est,
                                           const int level_min,
                                           const int level_max,
                                           rhea_domain_options_t
                                             *domain_options);

#endif /* RHEA_AMR_H */
