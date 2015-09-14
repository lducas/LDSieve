#ifndef GAUSSSIEVE_H
#define GAUSSSIEVE_H 

#include <NTL/quad_float.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>
#include "LinAlg.h"
#include "Lattice.h"


long* GaussSieve(t_BaseGS* base, long MaxVectors);

#endif