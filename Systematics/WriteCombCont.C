#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

double getWentries(TTree *tree, int iMCType, int icandType, TString Wname="weight");
TString round(double n, int e, double d=1.);

void WriteCombCont(int nFiles){
  double nBkg[2][8], nBkg2[2][8];
  for(int mc=0; mc<2; mc++){
    for(int cand=0; cand<8; cand++){
      nBkg[mc][cand] = 0;
      nBkg2[mc][cand] = 0;
    }
  }
  TChain *tree;
  TString CombContfileName = "AWG82/systematics/BF/CombContYields.txt";
  fstream CombContFile; CombContFile.open(CombContfileName,fstream::out);
  for(int file=0; file<=nFiles; file++){
    if(file%10==0) cout<<"Doing iteration "<<file<<" of "<<nFiles<<endl;
    TString ntupleName = "AWG82/systematics/BF/FitRoot/Fit"; ntupleName += file; ntupleName += "RAll_RunAll.root";
    tree = new TChain("ntp1"); tree->Add(ntupleName);
    for(int mc=0; mc>-2; mc--){
      for(int cand=1; cand<=8; cand++)	{
	double Entries = getWentries(tree,mc,cand);
	CombContFile<<Entries<<endl;
	nBkg[mc+1][cand-1] += Entries;
	nBkg2[mc+1][cand-1] += (Entries*Entries);
      }
      CombContFile<<endl;
    }
    CombContFile<<endl;
    tree->Delete();
  }
  nFiles++; //To calculate mean and RMS
  for(int mc=0; mc<2; mc++){
    for(int cand=0; cand<8; cand++){
      nBkg[mc][cand] /= (double)nFiles;
      nBkg2[mc][cand] = sqrt((nBkg2[mc][cand] - nBkg[mc][cand]*nBkg[mc][cand]*nFiles)/(double)(nFiles-1));
      cout<<"Bkg "<<cand<<" is     \t"<<round(nBkg[mc][cand],1)<<" +- "
	  <<round(nBkg2[mc][cand],1)<<",\t a "<<
	round(nBkg2[mc][cand]/nBkg[mc][cand]*100,2)<<" % error"<<endl;
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
double getWentries(TTree *tree, int iMCType, int icandType, TString Wname){
  double Wentries = 0; float weight=0;
  int MCType, candType;
  if(tree==0) return 0;
  tree->SetBranchAddress(Wname,&weight);
  tree->SetBranchAddress("MCType",&MCType);
  tree->SetBranchAddress("candType",&candType);
  for(int entry=0; entry<tree->GetEntries(); entry++){
    tree->GetEvent(entry);
    if(MCType==-2) MCType = -1;
    if(candType==icandType && MCType==iMCType) Wentries += weight;
  }
  return Wentries;
}










