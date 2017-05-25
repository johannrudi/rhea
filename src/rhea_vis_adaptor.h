/*
*/

#ifndef RHEA_VIS_ADAPTOR_HH
#define RHEA_VIS_ADAPTOR_HH

#include <sc.h>

/* tell the C++ compiler to use the C style name mangling */
SC_EXTERN_C_BEGIN;

void CatalystInitialize (int numScripts, char* scripts[]);

void CatalystFinalize ();

void CatalystCoProcess (unsigned int numberOfPoints, double* pointsData,
                        unsigned int numberOfCells, unsigned int* cellsData,
                        int order,
                        double* velocityData, double* pressureData);
                        //double time, unsigned int timeStep, int lastTimeStep);

SC_EXTERN_C_END;

#endif /* RHEA_VIS_ADAPTOR_HH */
