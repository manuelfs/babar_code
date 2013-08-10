#include "babar_code/PlotsThesis/PlotUtils.cc"
#include "TH1F.h"
#include "TChain.h"
#include "TString.h"
#include "TRandom3.h"

using namespace std;
using std::cout;
using std::endl;

#define Nmodes 868

void prd_efficiency(int seed=0){
  TRandom3 rand(seed); // 1 makes the first seed constant. 0 is a time dependent seed.
  double fmode[Nmodes], eff[Nmodes], norm=0;
  for(int mode=0; mode<Nmodes; mode++){
    fmode[mode] = rand.Uniform(1);
    norm += fmode[mode];
    eff[mode] = rand.Uniform(1);
  }
  for(int mode=0; mode<Nmodes; mode++) fmode /= norm;
}
