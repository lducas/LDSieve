#include <string.h>
#include <limits.h>
#include <iomanip>
#include <fstream>

#include "Lattice.h"
#include "LinAlg.h"
#include "ListDecode.h"
#include <algorithm>

using namespace NTL;
using namespace std;

#define m LD_NumberOfBlock


void quickSort(long* si, double* val, int debut, int fin){
	////////cerr << "Sort " << debut << " : " << fin << endl;
	long gauche = debut-1;    long droite = fin+1;    long sw;
	const double pivot = val[si[(gauche+droite)/2]];

	if(debut >= fin) 
		return;
	while(1)
	{
		do droite--; while(val[si[droite]] < pivot);
		do gauche++; while(val[si[gauche]] > pivot);

		if(gauche < droite){
			sw = si[droite];
			si[droite]=si[gauche];
			si[gauche]=sw;
		}
		else break;
	}
	quickSort(si,val, debut, droite);
	quickSort(si,val, droite+1, fin);
}

int b,nn;
long N,Psize;
double* Polytope;

long* TmpRes;
long TmpResCount;
double* tt;

double* TmpInnerProd;
long long* Sort;

double* sums,*minsums;

void improve_poly(double delta){
	double Norm;
	for (int i = 0; i < Psize; ++i)
	{
		for (int j = 0; j < Psize; ++j)
		{
			double scal = 0;
			Norm = sqrt(square_norm(b,&Polytope[i*b]));
			for (int k = 0; k < b; ++k)
				scal += Polytope[i*b +k] * Polytope[j*b+k];
			scal /=Norm;
			if (scal<0)
				scal = 0;

			for (int k = 0; k < b; ++k)
				Polytope[i*b +k] -= scal * delta * Polytope[j*b+k];
		}
		Norm = sqrt(square_norm(b,&Polytope[i*b]));
		for (int k = 0; k < b; ++k)
					Polytope[i*b +k] /= Norm;

	}
}


void LD_Setup(int n, int PPsize){
	nn = n;

	Psize = PPsize;
	if (n % m){
		cerr << "n must be a multiple of m"  << endl;
		exit(0);
	}
	b = nn/m;
	Polytope = new double[b*Psize];
	N = 1;
	for (int i = 0; i < m; ++i)
		N*= Psize;
	for (int i = 0; i < Psize; ++i)
	{
		spherical_sampling(b,&Polytope[b*i]);
	}

	cerr << "Improving Polytope" << endl;
	for (int i = 0; i < 100; ++i)
		improve_poly(.1);
	cerr << "Done" << endl;

	TmpRes = new long[N];
	//////cerr << "TmpRes Ok" << endl;
	tt = new double[n];
	TmpInnerProd = new double[m*Psize];
	Sort = new long long[m*Psize];
	//////cerr << "TmpOrderInnerProd Ok" << endl;
	sums = new double[m];
	minsums = new double[m];
}




#define deux28 (((long long) 2) << 28)

void LD_Search(long* I, long* t,double alpha){
//	return LD_Search_Slow(t,alpha);

	double Norm = sqrt(square_norm(nn,t));
	////cerr << Norm << endl;
	for (int i = 0; i < nn; ++i)
		tt[i] = ((double) t[i])/ (Norm *sqrt(m)) ;

	double* val = TmpInnerProd;

	for (int i = 0; i < m; ++i){
		for (int p = 0; p < Psize; ++p){
			val[i*Psize + p] = 0.;
			for (int k = 0; k < b; ++k){
				val[i*Psize + p] += tt[b*i + k] * Polytope[b*p + k];
			}
		}
	}

	// cerr_vec(Psize,&val[0*Psize]);

	// compact results and indices in 1 value

	for (int i = 0; i < m; ++i){
		for (int p = 0; p < Psize; ++p){
			Sort[i*Psize + p] = ((long) ceil(.5 + val[i*Psize + p] * deux28)) << 32;
			Sort[i*Psize + p] += p;
		}
		sort(&Sort[i*Psize], &Sort[(i+1)*Psize], greater<long>());
	}	

	// split the results
	for (int i = 0; i < m; ++i){
		for (int p = 0; p < Psize; ++p){
			val[i*Psize + p] = ((double) (Sort[i*Psize + p] >> 32)) / ((double) deux28);
			Sort[i*Psize + p] = (Sort[i*Psize + p] << 32) >> 32 ;
		}
	}	

//	cerr_vec(Psize,val);

	// cerr_vec(Psize,&Sort[0*Psize]);
	// cerr_vec(Psize,&val[0*Psize]);

	TmpResCount = 0;

	/// Setting Up Kannan-Style Enumeration
	int res, i = m-1;
	double* v0 = val;
	double* v1 = &val[1*Psize];
	double* v2 = &val[2*Psize];

	double sum = 0;
	sums[0] = val[0];

	for (int i = 0; i < m; ++i){
		minsums[i] = alpha;
		for (int j = i+1; j < m; ++j)
		{
			minsums[i] -= val[j*Psize + 0];
		}
	}

	for (int j0 = 0; (j0 < Psize) && (v0[j0] >= minsums[0]); ++j0){
		sums[0] = val[0*Psize + j0];
		for (int j1 = 0; (j1 < Psize) && (v1[j1] + sums[0] >= minsums[1]); ++j1){
			sums[1] = sums[0] + val[1*Psize + j1];
			for (int j2 = 0; (j2 < Psize) && (v2[j2] + sums[1] >= minsums[2]); ++j2){
				res = 	Psize*Psize*Sort[0*Psize + j0] 
							+ Psize*Sort[1*Psize + j1] 
					        +       Sort[2*Psize + j2];
			    // cerr << Sort[0*Psize + j0] << "," 
			    // 	 << Sort[1*Psize + j1] << ","
			    // 	 << Sort[2*Psize + j2] << "  ";

				TmpResCount++;
				I[TmpResCount] = res;			
			}
		}
	}
	I[0] = TmpResCount;
}



// long* LD_Search_Slow(long* t,double alpha){
// 	TmpResCount =0;
// 	double* v = new double[nn];
// 	double sum;
// 	double Norm = sqrt(square_norm(nn,t));
// 	for (int i = 0; i < nn; ++i){
// 		tt[i] = ((double) t[i])/ (Norm *sqrt(m)) ;
// 	}
// 	for (int a1 = 0; a1 < Psize; ++a1)
// 	{
// 		memcpy(&(v[0*b]),&Polytope[b*a1],b*sizeof(double));
// 		for (int a2 = 0; a2 < Psize; ++a2)
// 		{
// 			memcpy(&(v[1*b]),&Polytope[b*a2],b*sizeof(double));
// 			for (int a3 = 0; a3 < Psize; ++a3)
// 			{
// 				memcpy(&(v[2*b]),&Polytope[b*a3],b*sizeof(double));
// 				sum = 0;
// 				for (int i = 0; i < nn; ++i)
// 					sum += tt[i]*v[i];
// 				//cerr << sum << ":" << TmpResCount;
// 				if (sum > alpha){
// 					TmpResCount++;
// 					TmpRes[TmpResCount] = a1*Psize*Psize + a2*Psize + a3;						
// 				}
// 			}
// 		}
// 	}
// 	TmpRes[0] = TmpResCount;
// 	cerr << TmpResCount << endl;
// 	return TmpRes;
// }
