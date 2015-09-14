#include <string.h>
#include <limits.h>
#include <iomanip>
#include <fstream>
#include "enum.h"
#include "Lattice.h"
#include "LinAlg.h"

using namespace NTL;
using namespace std;

long VectorsIndex;
long* Vectors;
double* Norms;
long* UpdateList;
long UpdateIndex;
long* temp_reduce;
long n;

double ShortestNorm;
long ShortestIndex;
double sigma;
double* zero;
long counter;


void GS_Setup(t_BaseGS* base, long MaxVectors){
	n = base->dim;
	Vectors = new long[n*MaxVectors];
	Norms  = new double[MaxVectors];
	UpdateList = new long[MaxVectors];
	temp_reduce = new long[n];
	zero = new double[n];
	memset(zero,0.0,n*sizeof(double));
	VectorsIndex = 0;
	UpdateIndex = 0;
	ShortestNorm = 1000.*n*base->r[0];
	sigma = .5 * sqrt(base->r[0]);
	counter = 0;
}


void PushUpdate(long x){
////	cerr << "Pushing " << x << endl;
	for (int i = 0; i < UpdateIndex; ++i)
	{
		if (UpdateList[i]==x)
			return;
	}
	UpdateList[UpdateIndex]= x;
	UpdateIndex++;
//	cerr << "Pushed " << endl;
} 

long PopUpdate(){
//	cerr << "Poping" << UpdateList[UpdateIndex-1] << endl;

	UpdateIndex--;
	return UpdateList[UpdateIndex];
}

int reduce(t_BaseGS* base,long a, long b){
//	cerr << a << " " << b<< endl;	
	counter++;
	for (int i = 0; i < n; ++i){
		temp_reduce[i] = Vectors[a*n+i]-Vectors[b*n+i];
	}
	double Norm = square_norm(n,temp_reduce);
	//cerr << Norm << endl;
	if (Norm ==0.){
		return 0;
	}

	if (Norm < Norms[a]){
		memcpy(&Vectors[a*n],temp_reduce,n*sizeof(long));
		Norms[a] = Norm;
		PushUpdate(a);
		if (Norm<ShortestNorm){
			ShortestNorm = Norm;
			ShortestIndex = a;
		}
		return 1;
	}
	if (Norm < Norms[b]){
		memcpy(&Vectors[b*n],temp_reduce,n*sizeof(long));
		Norms[b] = Norm;
		PushUpdate(b);
		if (Norm<ShortestNorm){
			ShortestNorm = Norm;
			ShortestIndex = a;
		}
		return 0;
	}
	return 0;
}

long* GaussSieve(t_BaseGS* base, long MaxVectors){
	GS_Setup(base,MaxVectors);
	int a;
	double Norm;
	long* z = new long[n];
	int aaa = 0;
	int bbb = 0;
	while(true){
		aaa++;
		if (!(aaa % 1000))
			cerr << VectorsIndex << "   " << counter / 1000000 <<  "M  " << sqrt(ShortestNorm) << endl;
		if (!UpdateIndex){
			if (VectorsIndex==MaxVectors){
				cerr << "Done after that many reduce trials" << counter << endl;
							return &Vectors[ShortestIndex*n];
						}
		//	cerr << "sampling" << endl;
			gaussian_sampling(base,zero,z ,sigma);
			B_To_Canon(base,&Vectors[VectorsIndex*n],z);
			Norm = square_norm(n,&Vectors[VectorsIndex*n]);
			Norms[VectorsIndex] = Norm;

			PushUpdate(VectorsIndex);
			VectorsIndex++;
		}
		a = PopUpdate();
		for (int b = 0; b < VectorsIndex; ++b)
		{
			bbb = 0;
			if (a!=b)
				bbb = reduce(base,a,b);
			if (bbb)
				break;
		}
	}
}
