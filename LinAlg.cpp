#include <string.h>
#include <limits.h>
#include <iomanip>
#include <fstream>

#include <NTL/quad_float.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>


using namespace NTL;
using namespace std;

double square_norm_rel(long n, double* x,double* r){
  int i;
  double d=0;
  for(i=0;i<n;i++){
    d+= x[i]*x[i]*r[i];
  }
  return d;
}
double square_norm(long n, long* x){
  int i;
  double d=0;
  for(i=0;i<n;i++){
    d+= ((double) x[i])*((double) x[i]);
  }
  return d;
}
double square_norm(long n, double* x){
  int i;
  double d=0;
  for(i=0;i<n;i++){
    d+= x[i]*x[i];
  }
  return d;
}


void MatVecProd(long n, double* m, double* v,double* res){
  int i,j;
  for(i=0;i<n;i++){
    res[i]=0;
    for(j=0;j<n;j++){
      res[i] += v[j]*m[n*j+i];
    }
  }
}
void MatVecProd(long n, double* m, long* v,double* res){
  int i,j;
  for(i=0;i<n;i++){
    res[i]=0;
    for(j=0;j<n;j++){
      res[i] += v[j]*m[n*j+i];
    }
  }
}
void IntMatVecProd(long n, long* m, long* v,long* res){
  int i,j;
  for(i=0;i<n;i++){
    res[i]=0;
    for(j=0;j<n;j++){
      res[i] += v[j]*m[n*j+i];
    }
  }
}

double* To_mat_double(mat_RR& M){
  long n = M.NumRows();
  long m = M.NumCols(); long i,j; double f; double *res = new double[n*m];
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      conv(f,M(i+1,j+1));
      res[m*i+j] = f;
    }
  }
  return res;
}

double* To_vec_double(vec_RR& v){
  long n = v.length();
  long i,j; double f; double *res = new double[n];
  for(i=0;i<n;i++){
      conv(f,v(i+1));
      res[i] = f;
    }
  return res;
}

void conv_mat(mat_RR& res,mat_ZZ& B){
  long n = B.NumRows();
  long m = B.NumCols(); long i,j;
  res.SetDims(n,m);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      conv(res(i+1,j+1),B(i+1,j+1));
    }
  }
}

void conv_mat(mat_ZZ& res,mat_RR& B){
  long n = B.NumRows();
  long m = B.NumCols(); long i,j;
  res.SetDims(n,m);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      conv(res(i+1,j+1),floor(B(i+1,j+1)+.5) );
    }
  }
}
long* To_mat_long(mat_ZZ& M){
  long n = M.NumRows(); long m = M.NumCols(); long i,j; long f;
  ZZ z,max= to_ZZ(LONG_MAX);
  long *res = new long[n*m];
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      z = M(i+1,j+1);
      if (abs(z)>= max){
	cerr << "Conversetion to Long Failed : Matrix too Large" << endl;
	exit(0);
      }
      conv(f,z);
      res[m*i+j] = f;
    }
  }
  return res;
}

void cerr_vec(long n,double* v){
  cerr << "[";
  for (int i=0;i<n;i++)
    cerr << v[i] << " ";
  cerr << "]" << endl;
}


void cerr_vec(long n,long* v){
  cerr << "[";
  for (int i=0;i<n;i++)
    cerr << v[i] << " ";
  cerr << "]" << endl;
}


void cerr_vec(long n,long long* v){
  cerr << "[";
  for (int i=0;i<n;i++)
    cerr << v[i] << " ";
  cerr << "]" << endl;
}

