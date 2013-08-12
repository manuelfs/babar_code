#include "babar_code/PlotsThesis/PlotUtils.cc"
#include "TH1F.h"
#include "TChain.h"
#include "TString.h"
#include "TRandom3.h"
#include <fstream>

using namespace std;
using std::cout;
using std::endl;

int findMode(double x, double cfmode[], int Nmodes) {
  int minMode=0, maxMode=Nmodes;
  int midMode = (maxMode-minMode)/2+minMode, rep=0;
  while(maxMode>minMode+1 && rep <20){
    // cout<<x<<" is between cfmode["<<minMode<<"] = "<<cfmode[minMode]
    // 	<<" and cfmode["<<maxMode<<"] = "<<cfmode[maxMode]
    // 	<<", with cfmode["<<midMode<<"] = "<<cfmode[midMode]<<endl;
    if(x<cfmode[midMode]) maxMode = midMode;
    else minMode = midMode;
    midMode = (maxMode-minMode)/2+minMode;
    rep++;
  }
  return minMode;
}

void prd_efficiency(long Nevents=1200, long Nrep=10, double effmin=0.02, double effmax=0.08, 
		    int seed=1){
  TString filename = "babar_code/prd/fraction_eff.txt";
  ifstream infile(filename);
  TRandom3 rand(seed); // 1: first seed constant - 0: time dependent seed.
  double fmode[1000], cfmode[1000], eff[1000], f_eff, norm=0, AverEff=0;
  int Nmodes=0, imode;
  while(infile>> f_eff){
    eff[Nmodes] = rand.Uniform(effmin, effmax);
    fmode[Nmodes] = f_eff/eff[Nmodes];
    norm += fmode[Nmodes];
    Nmodes++;
  }
  cfmode[0] = 0;
  for(int mode=0; mode<Nmodes; mode++) {
    fmode[mode] /= norm;
    cfmode[mode+1] = fmode[mode] + cfmode[mode];
    AverEff += fmode[mode]*eff[mode];
  }
  double Ntotal = Nevents/AverEff;
  double mean = 0, rms = 0;
  for(int rep=0; rep<Nrep; rep++){
    double Npass=0;
    for(int eve=0; eve < Ntotal; eve++) {
      imode = findMode(rand.Uniform(1), cfmode, Nmodes);
      if(rand.Uniform(1) < eff[imode]) Npass++;
    }
    double effrep = Npass/Ntotal;
    mean += effrep; rms += effrep*effrep;
  }
  double DNrep = (double)Nrep;
  mean /= DNrep; rms = sqrt((rms - mean*mean*DNrep)/(DNrep-1));

  int digits = 4;
  cout<<"Paper: ("<<RoundNumber(AverEff*100,digits)<<" +- "<<RoundNumber(sqrt(Nevents)*100,digits,Ntotal)
      <<")%  -  Simulation: ("<<RoundNumber(mean*100, digits)<<" +- "<<RoundNumber(rms*100, digits)
      <<")%"<<endl;
}

void calc_fractions(){
  TChain chain("ntp1");
  chain.Add("AWG82/ntuples/small/RAll_RunAll.root");

  int candBMode;
  chain.SetBranchAddress("candBMode",&candBMode);

  double events[10000];
  for(int eve=0; eve < 10000; eve++) events[eve] = 0;
  for(int eve=0; eve < chain.GetEntries(); eve++){
    chain.GetEvent(eve);
    events[candBMode-10000]++;
  }

  TString filename = "babar_code/prd/fraction_eff.txt";
  ofstream outfile(filename);
  double total = 0;
  int Nmodes = 0;
  for(int eve=0; eve < 10000; eve++) total += events[eve];
  for(int eve=0; eve < 10000; eve++) {
    if(events[eve] > 0) {
      outfile << RoundNumber(events[eve]*100, 5, total) << " ";
      Nmodes++;
      if(Nmodes%10 == 0) outfile<<endl;
    }
  }
  cout<<"Printed "<<Nmodes<<" in "<<filename<<endl;
  outfile.close();
}

