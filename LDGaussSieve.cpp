#include <string.h>
#include <limits.h>
#include <iomanip>
#include <fstream>

#include "Lattice.h"
#include "LinAlg.h"
#include "ListDecode.h"

#include "Buckets.cpp"
#include <algorithm>


// Uncomment the following to run standard GaussSieve
//#define NO_BUCKETING

using namespace NTL;
using namespace std;

long VectorsIndex;
long* Vectors;
long TimeStamp;
long** UpdateList;
long UpdateIndex;
long* temp_reduce;
long* temp_neg;
long n;
long* Isearch;
long* Iadd;
long* Irem;

double ShortestNorm;
long* Shortest;
double sigma;
double* zero;
long counter;

t_BaseGS* base_glob;

double alpha,beta;

long number_of_buckets;

//typedef unordered_multimap<long,long*> t_HT;

//t_HT HashTable;


double Current_ATime, Tmp_ATime;
void Init_ATime() { Current_ATime = GetTime(); }

void Elapsed_ATime()
{
  Tmp_ATime = GetTime()-Current_ATime;
  cout << Tmp_ATime;
}

void Update_ATime()
{
  Elapsed_ATime();
  Init_ATime();
}


uint64_t get_cycles(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}




inline long long_square_norm(long* x){
  long d=0;
  for(int i=0;i<n;i++){
    d+= x[i]*x[i];
  }
  return d;
}




void LDGS_Setup(t_BaseGS* base, long MaxVectors, long Psize){
	n = base->dim;
	Vectors = new long[(n+3)*MaxVectors];
	UpdateList = new long*[MaxVectors];
	temp_reduce = new long[n];
	temp_neg = new long[n];
	zero = new double[n];	
	memset(Vectors,0,(n+3)*MaxVectors*sizeof(long));
	memset(zero,0.0,n*sizeof(double));
	VectorsIndex = 0;
	UpdateIndex = 0;
	ShortestNorm = 1000.*n*base->r[0];
	sigma = .2 * sqrt(base->r[n/2]);
	counter = 0;
	number_of_buckets =0;
	base_glob = base;

#ifndef NO_BUCKETING
	number_of_buckets = 1;
	for (int i = 0; i < LD_NumberOfBlock; ++i)
	{
		number_of_buckets*=Psize;
	}

	Isearch = new long[number_of_buckets];
	Irem = new long[number_of_buckets];
	Iadd = new long[number_of_buckets];
	LD_Setup(n,Psize);
	Buckets_Setup(number_of_buckets);
#endif

}

inline void PushUpdate(long* x){
	if (!x[n+2]){
		x[n+2] =1;
		UpdateList[UpdateIndex]= x;
		UpdateIndex++;
		}
} 


inline long* PopUpdate(){
	UpdateIndex--;
	UpdateList[UpdateIndex][n+2]=0;
	return UpdateList[UpdateIndex];
}


#ifndef NO_BUCKETING
inline void AddToTable(long* a,double delta){
	LD_Search(Iadd,a, delta);
	for (int i = 1; i <= Iadd[0]; ++i)
	{
		Add_To_Bucket(a,Iadd[i]);
	}
//	cerr << "Done Add" << endl;
}

inline void RemFromTable(long* a,double delta){
	LD_Search(Irem,a, delta);
	for (int i = 1; i <= Irem[0]; ++i)
	{ 
		Rem_From_Bucket(a,Irem[i]);	
	}
//	cerr << "Done Rem" << endl;
}
#endif



int doreduce(long* a, long* b){

	counter++;
	long ip = 0;
	for (int i = 0; i < n; ++i){
		ip += a[i]*b[i];
	}
	int sign = (ip>0? -1:1);
	long Norm = a[n] + b[n] + sign* 2*ip;

//	cerr << ip << " " << a[n] + b[n] << " " << Norm << endl;
	if (Norm == 0){
		if (a==b)
			return 0;

		#ifndef NO_BUCKETING
			RemFromTable(b,alpha);
		#endif
		long* z = new long[n];
		gaussian_sampling(base_glob,zero,z ,sigma);
		B_To_Canon(base_glob,b,z);
		Norm = long_square_norm(b);
		b[n] = Norm;
		#ifndef NO_BUCKETING
			AddToTable(b,alpha);
		#endif		
		PushUpdate(b);
		delete z;
		return 2;
	}

	if (Norm < b[n]){
		#ifndef NO_BUCKETING
			RemFromTable(b,alpha);
		#endif
		for (int i = 0; i < n; ++i)
		{
			b[i] += sign*a[i];
		}
		b[n] = Norm;
		#ifndef NO_BUCKETING
			AddToTable(b,alpha);
		#endif		
		PushUpdate(b);
		if (Norm<ShortestNorm){
			ShortestNorm = Norm;
			Shortest = b;
		}
		return 2;
	}

	if (Norm < a[n]){
		#ifndef NO_BUCKETING
			RemFromTable(a,alpha);
		#endif				
		for (int i = 0; i < n; ++i)
		{
			a[i] += sign*b[i];
		}
		a[n] = Norm;
		#ifndef NO_BUCKETING		
			AddToTable(a,alpha);
		#endif		
		PushUpdate(a);
		if (Norm<ShortestNorm){
			ShortestNorm = Norm;
			Shortest = a;
		}
		return 1;
	}

	return 0;
}

inline int reduce(long* a, long* b){

	counter++;
	long ip = 0;
	for (int i = 0; i < n; ++i){
		ip += a[i]*b[i];
	}
	int sign = (ip>0? -1:1);
	long Norm = a[n] + b[n] + sign* 2*ip;
	if (Norm < b[n] || Norm<a[n])
		return doreduce(a,b);
	return 0;
}



#ifndef NO_BUCKETING
long SearchAndReduce(long* a,double delta){

	long sign = 1;

	LD_Search(Isearch,a, delta );
	int res;
	long* b;

	retry_negated:
	TimeStamp ++;
	a[n+1] = TimeStamp;

	for (int i = 1; i <= Isearch[0]; ++i){
		long** bucket = Get_Bucket(Isearch[i]);
		for (long j =0; bucket[j]; j++) 
		{
			b = bucket[j];
			if ( TimeStamp > b[n+1])
			{ 
				b[n+1] = TimeStamp;
				res = reduce(a,b);
				if (res){
					if (res == 1)
						return 1;
					if (res == 2)
						j--;
				}
			}
		}
	}

	if (sign==-1)
		return 0;

	sign = -1;
    for (int i = 0; i < n; ++i)
	  temp_neg[i] = -a[i];
	LD_Search(Isearch, temp_neg, delta );
	goto retry_negated;
}
#endif

long* LDGaussSieve(t_BaseGS* base, long MaxVectors,long Psize,double aalpha, double bbeta, double goal){
	LDGS_Setup(base,MaxVectors,Psize);

	long* a;
	double Norm;
	long* z = new long[n];
	int aaa = 0;
	int bbb = 0;
	alpha = aalpha;
	beta = bbeta;


	Init_ATime();
	unsigned long cycle_start = get_cycles();
	while(true){
		aaa++;
		if (!(aaa % 50000))
			cerr << "Vecs : " << VectorsIndex << "\t" << "Store / Buc : " << TotalSize << "/" 
				 << number_of_buckets << " Reds :" << ceil(counter / 100000.)/10. <<  "M\t shortest :" << sqrt(ShortestNorm) << endl;


		if ((VectorsIndex==MaxVectors) || (sqrt(ShortestNorm) <=goal)){
			cerr << "Sieving Done !" << endl;
			if (sqrt(ShortestNorm) <=goal){
				cout << (get_cycles() - cycle_start)/1000000  << " ";
				Elapsed_ATime();
				cout << " " 
					 << " "  << n
					 << " " << Psize << " " << aalpha << " " << bbeta 
					 << " " <<  VectorsIndex 
					 << " " << TotalSize
					 << " " << counter << "\n";
				}
			cerr << "Vecs : " << VectorsIndex << "\t" << "Store / Buc : " << TotalSize << "/" 
				 << number_of_buckets << "Reds :" << ceil(counter / 100000.)/10. <<  "M\t shortest :" << sqrt(ShortestNorm) << endl;
			return Shortest;
			}

		if (!UpdateIndex){
			gaussian_sampling(base,zero,z ,sigma);
			a = &Vectors[VectorsIndex*(n+3)];
			B_To_Canon(base,a,z);
			Norm = long_square_norm(a);
			a[n] = Norm;

#ifndef NO_BUCKETING
			AddToTable(a,alpha);
#endif		
			PushUpdate(a);
			VectorsIndex++;
		}
		a = PopUpdate();
#ifndef NO_BUCKETING
		SearchAndReduce(a,beta);
#else		
		for (int i = 0; i < VectorsIndex; ++i)
		{
			bbb = 0;
			long* b = &Vectors[(n+3)*i];
			if (reduce(a,b) == 1)
				break;
		}




#endif
	}
}
