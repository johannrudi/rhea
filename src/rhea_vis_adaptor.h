/*
*/

#ifndef RHEA_VIS_ADAPTOR_H
#define RHEA_VIS_ADAPTOR_H

#ifdef RHEA_USE_CATALYST
#include <sc.h>

/* tell the C++ compiler to use the C style name mangling */
SC_EXTERN_C_BEGIN;

void                rhea_vis_adaptor_initialize (const char *catalyst_scripts[],
                                                 const int n_catalyst_scripts);

void                rhea_vis_adaptor_finalize ();

void                rhea_vis_adaptor_process (unsigned int order,
                                              unsigned int nPoints,
                                              const double *pointCoords,
                                              unsigned int nCells,
                                              unsigned int *cellToCoordIdx,
                                              const double *velocityData,
                                              const double *pressureData);

SC_EXTERN_C_END;
#endif /* RHEA_USE_CATALYST */

#endif /* RHEA_VIS_ADAPTOR_H */
