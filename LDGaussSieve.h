#ifndef LDGAUSSSIEVE_H
#define LDGAUSSSIEVE_H 

#include <NTL/quad_float.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>
#include "LinAlg.h"
#include "Lattice.h"


long* LDGaussSieve(t_BaseGS* base, long MaxVectors,long Psize,double aalpha, double bbeta, double goal);

#endif