#include "TChain.h"
#include "TString.h"
#include "TRandom3.h"
#include "TFile.h"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using std::cout;
using std::endl;

TString round(double n, int e, double d=1.);

void GenerateBkg(int nFiles){
  //TString genName  = "AWG82/ntuples/small/FitRAllWe1x100_RunAll.root";
  TString genName  = "AWG82/ntuples/small/FitRAllNewx100_RunAll.root";
  TString CombContfileName = "AWG82/systematics/BF/CombContYieldsStat.txt";
  fstream CombContFile; CombContFile.open(CombContfileName,fstream::out);
  double Bkg[8][2], NBkg[8][2];
  for(int cand=0; cand<8; cand++) 
    for(int i=0; i<2; i++) {Bkg[cand][i]=0;NBkg[cand][i]=0;}
  float weight;
  int candType, MCType;
  TChain genChain("ntp1");
  genChain.Add(genName);
  genChain.SetBranchAddress("weight",&weight);
  genChain.SetBranchAddress("candType",&candType);
  genChain.SetBranchAddress("MCType",&MCType);
  for(int entry=0; entry<genChain.GetEntries(); entry++){
    genChain.GetEvent(entry);
    for(int cand=0; cand<8; cand++){
      if(MCType==0 && candType==cand+1) {Bkg[cand][0] += weight; NBkg[cand][0]++;}
      if(MCType<0 && candType==cand+1)  {Bkg[cand][1] += weight; NBkg[cand][1]++;}
    }
  }
   
  TRandom3 rand(0);
  double meanBkg[8][2], rmsBkg[8][2];
  for(int cand=0; cand<8; cand++){
    for(int i=0; i<2; i++){
      meanBkg[cand][i] = 0; rmsBkg[cand][i] = 0;
    }
  }
  for(int file=0; file<nFiles; file++){
    if(file%2000==0) cout<<"Doing iteration "<<file<<" of "<<nFiles<<endl;
    for(int i=0; i<2; i++){
      for(int cand=0; cand<8; cand++){
	double randYield = rand.PoissonD(NBkg[cand][i])*Bkg[cand][i]/NBkg[cand][i];
	CombContFile<<RoundNumber(randYield,2)<<endl;
	meanBkg[cand][i] += randYield;
	rmsBkg[cand][i] += (randYield*randYield);
      }
      CombContFile<<endl;
    }
    CombContFile<<endl;
  }

  for(int i=0; i<2; i++){
    for(int cand=0; cand<8; cand++){
      meanBkg[cand][i] /= (double)nFiles;
      rmsBkg[cand][i] = sqrt((rmsBkg[cand][i] - meanBkg[cand][i]*meanBkg[cand][i]*nFiles)/(double)(nFiles-1));
      cout<<"Bkg "<<cand<<" is     \t"<<round(meanBkg[cand][i],1)<<" +- "<<round(rmsBkg[cand][i],1)<<",\t a "<<
	round(rmsBkg[cand][i]/meanBkg[cand][i]*100,2)<<" % error"<<endl;
    }
    cout<<endl;
  }

  cout<<"Written "<<CombContfileName<<endl;
  CombContFile.close();
}

TString round(double n, int e, double d){
  if(d==0) return " - ";
  double neg = 1; if(n*d<0) neg = -1;
  double b = (int)(neg*n/d*pow(10.,(double)e)+0.5);
  b /= pow(10.,(double)e)*neg;
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(!result.Contains(".") && e != 0) result += ".";
  
  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<e-afterdot.Length(); i++)
    result += "0";
  return result;
}
