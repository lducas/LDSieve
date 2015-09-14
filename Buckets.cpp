#include <string.h>
#include <limits.h>
#include <iomanip>
#include <fstream>
#include <iostream>
#include "Buckets.h"

using namespace std;

long nb_buckets;
long max_Bucket_size;
long* Buckets_size;
long*** Buckets;
long** tmp_res;
long TotalSize;


// Buckets[i][j] is the i-th entry of bucket j (not the opposite)
// That way, increasing the bucket size is simply done by allocating a new
// Bucket[max+1] of a fixed size.

void Buckets_Setup(long nb_bucketss){
	TotalSize = 0;
	nb_buckets = nb_bucketss;
	max_Bucket_size = IncreaseSizeStep;
	Buckets_size = new long[nb_buckets];
	memset(Buckets_size,0,nb_buckets*sizeof(long));
	tmp_res = new long*[MaxMaxBucketsSize+1];
	Buckets = new long**[MaxMaxBucketsSize];
	for (int i = 0; i < max_Bucket_size; ++i){
		Buckets[i] = new long*[nb_buckets];
	}

} 

void Buckets_Increase_Size(){
	for (int i = max_Bucket_size; i < max_Bucket_size+IncreaseSizeStep; ++i)
	{
		Buckets[i] = new long*[nb_buckets];
	}
	max_Bucket_size+=IncreaseSizeStep;
}

inline void Add_To_Bucket(long* a, long bu){
	if (Buckets_size[bu] == max_Bucket_size)
		Buckets_Increase_Size();
	Buckets[Buckets_size[bu]][bu] = a;
	Buckets_size[bu]++;
	TotalSize++;
}

inline void Rem_From_Bucket(long* a, long bu){
	for (int i = 0; i < Buckets_size[bu]; ++i){
		if (Buckets[i][bu]==a){
			Buckets_size[bu]--;
			Buckets[i][bu] = Buckets[Buckets_size[bu]][bu];
			TotalSize--;
			return;
		}
	}
}

inline long** Get_Bucket(long bu){
//	cerr << bu;
	for (int i = 0; i < Buckets_size[bu]; ++i){
//		cerr << " " << i;
		tmp_res[i] = Buckets[i][bu];
	}
//	cerr << endl;
	tmp_res[Buckets_size[bu]] = NULL;
	return tmp_res;
}


 