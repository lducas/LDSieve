#ifndef LINALG_H
#define LINALG_H 

#include <NTL/quad_float.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>


using namespace NTL;
using namespace std;

double square_norm_rel(long n, double* x,double* r);
double square_norm(long n, double* x);
double square_norm(long n, long* x);

void MatVecProd(long n, double* m, double* v,double* res);
void MatVecProd(long n, double* m, long* v,double* res);
void IntMatVecProd(long n, long* m, long* v,long* res);
double* To_mat_double(mat_RR& M);
double* To_vec_double(vec_RR& v);
void conv_mat(mat_RR& res,mat_ZZ& B);
void conv_mat(mat_ZZ& res,mat_RR& B);
long* To_mat_long(mat_ZZ& M);
void cerr_vec(long n,double* v);
void cerr_vec(long n,long* v);
void cerr_vec(long n,long long* v);  

#endif
