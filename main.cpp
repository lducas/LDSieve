#include <string.h>
#include <limits.h>
#include <iomanip>
#include <fstream>

#include <NTL/quad_float.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>

#include "LDGaussSieve.h"
#include "LinAlg.h"
#include "Lattice.h"
#include "math.h"
#include "ListDecode.h"

using namespace NTL;
using namespace std;

double Current_Time, Tmp_Time;
void Init_Time() { Current_Time = GetTime(); }
void Display_Time(double t)
{
  cerr << "\t Time  = " << t << "sec" << endl;
}

void Elapsed_Time()
{
  Tmp_Time = GetTime()-Current_Time;
  Display_Time(Tmp_Time);
}

void Update_Time()
{
  Elapsed_Time();
  Init_Time();
}

int main(int argc, char *argv[]){
  int n;
  mat_ZZ B_ZZ;
  double ratio_len;
  double e,e2;

  RR::SetPrecision(150);

  if (argc!=6){
    cerr << "Usage : \n cat filename | ldsieve N g C a b " << endl;
    cerr << "filename : basis as a square matrix (in text NTL format)" << endl;
    cerr << "N : maximal number of vectors" << endl;
    cerr << "g : square of goal norm" << endl;
    cerr << "C : Size of each component of the code" << endl;
    cerr << "a : storing sphere-cap parameter alpha " << endl;
    cerr << "b : search sphere-cap parameter beta " << endl;
    exit(1);
  }


  double max_vecs   = atof(argv[1]);
  double goal       = sqrt(atof(argv[2]));
  long Psize        = atoi(argv[3]);
  double alpha      = atof(argv[4]);
  double beta       = atof(argv[5]);

  cerr << "Reading matrix from stdin ... ";
  cin >> B_ZZ;
  cerr << "Done";

  n = B_ZZ.NumCols();

  cerr << "dim         = " << n << endl;
  cerr << "goal norm   = " << goal << endl;
  cerr << "RPC size C  = " << Psize << endl;
  cerr << "alpha       = " << alpha << endl;
  cerr << "beta        = " << beta << endl;
  cerr << "RPC #blocs m = " << LD_NumberOfBlock << "  (Edit ListDecode.* to change) " << endl;

  Init_Time();
  t_BaseGS* base = new_BaseGS(B_ZZ);  
  Elapsed_Time();

  long* v;
  float norm;
  
  Init_Time();
  v = LDGaussSieve(base,max_vecs, Psize, alpha,beta,goal); 
  Elapsed_Time();
  if (!v)
  {
    cerr << "Failed : no result" << endl;
    exit(0);
  }
//  cerr_vec(n,v);
  norm = sqrt(square_norm(n,v));
  cerr << "Shortest vector Found " << norm << "(goal : " << goal << ")" << endl;
  if (norm < 1.001*goal){
    cerr << "Success " << endl;
  }
  exit(0);
}
