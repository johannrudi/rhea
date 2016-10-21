#ifndef RHEA_DISCRETIZATION_H
#define RHEA_DISCRETIZATION_H

#include <ymir_options.h>

/* global options */
extern int          rhea_discretization_level_min;

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_discretization_add_options (ymir_options_t * opt_sup);

#endif /* RHEA_DISCRETIZATION_H */
