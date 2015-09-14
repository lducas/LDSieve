#ifndef LATTICE_H
#define LATTICE_H 

#include <NTL/quad_float.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>
#include "LinAlg.h"

typedef struct {
  long dim;
  long* B;
  double* Binv;
  double* muinv;
  double* Bstar;
  double* mu;
  double* r;
  double D;
} t_BaseGS;

t_BaseGS* new_BaseGS(mat_ZZ B_ZZ);

void spherical_sampling(int n,double *res);
void gen_rand_lattice(mat_ZZ& B,long n);
double* B_To_Bstar(t_BaseGS* base, long* x);
long* Babai_Bstar_To_B(t_BaseGS* base, double *x);
double* Canon_To_Bstar(t_BaseGS* base, double* x);
long* B_To_Canon(t_BaseGS* base, long* x);
void B_To_Canon(t_BaseGS* base, long* res,long* x);
long* Canon_To_B(t_BaseGS* base, long* x);
double* B_To_Bstar(t_BaseGS* base, long* x);
double* Bstar_To_Canon(t_BaseGS* base, double* x);
t_BaseGS* randomized_BaseGS(mat_ZZ B_ZZ);
void gaussian_sampling(t_BaseGS* base, double *x,long* z,double sigma);

float Predict_attack(t_BaseGS* base, int i1, int i2);

#endif