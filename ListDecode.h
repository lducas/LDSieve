#ifndef LISTDECODING_H
#define LISTDECODING_H 
#include <string.h>
#include <limits.h>
#include <iomanip>
#include <fstream>

#include "Lattice.h"
#include "LinAlg.h"

#define LD_NumberOfBlock 3


void LD_Setup(int n, int PPsize);
void LD_Search(long* I, long* t,double alpha);

// long* LD_Search(long* t,double alpha,long min_out, long max_out);

// long* LD_Search_Slow(long* t,double alpha);
#endif