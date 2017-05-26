/*
*/

#ifndef RHEA_VIS_H
#define RHEA_VIS_H

#include <ymir_vec_ops.h>

void                rhea_vis_initialize (const char *catalyst_scripts[],
                                         const int n_catalyst_scripts);

void                rhea_vis_finalize ();

void                rhea_vis_process_primary (ymir_vec_t *velocity,
                                              ymir_vec_t *pressure);

#endif /* RHEA_VIS_H */
