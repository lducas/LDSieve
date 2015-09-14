#ifndef BUCKETS_H
#define BUCKETS_H 

#include <string.h>
#include <limits.h>
#include <iomanip>
#include <fstream>


#define MaxMaxBucketsSize 10000000
#define IncreaseSizeStep 1


void Buckets_Setup(long nb_bucketss);
void Buckets_Increase_Size();
inline void Add_To_Bucket(long* a, long bu);
inline void Rem_From_Bucket(long* a, long bu);
inline long** Get_Bucket(long bu);





#endif