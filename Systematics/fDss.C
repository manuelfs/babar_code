#include "TString.h"
#include "TChain.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

TString round(double n, int e, double d=1.);

void fDss(int nFiles, int doW = 1){

  TString folder = "AWG82/systematics/BF/FitRoot/Fit"; 
  TString channels[] = {"$D^0$","$D^{*0}$","$D^+$","$D^{*+}$","$D$","$D^{*}$"};
  
  float weight;
  int candType, MCType;
  double mean[6], rms[6];
  for(int i=0; i<6; i++){mean[i] = 0; rms[i] = 0;}
  for(int file=0; file<nFiles; file++){
    TString genName = folder; genName += file; genName += "RAll_RunAll.root";
    TChain *genChain = new TChain("ntp1");
    genChain->Add(genName);
    double totW[2][6], totW2[2][6];
    for(int i=0; i<4; i++){
      for(int j=0; j<2; j++){
	totW[j][i] = 0; totW2[j][i] = 0; 
      }
    }
    genChain->SetBranchAddress("weight",&weight);
    genChain->SetBranchAddress("candType",&candType);
    genChain->SetBranchAddress("MCType",&MCType);
    for(int entry=0; entry<genChain->GetEntries(); entry++){
      genChain->GetEvent(entry);
      if(doW==0) weight = 1;

      if(MCType==14) {
	totW[candType>4][(candType-1)%4] += weight;
	totW2[candType>4][(candType-1)%4] += weight*weight;
      }
    }
    for(int chan = 0; chan<2; chan++){
      for(int j=0; j<2; j++){
	totW[j][chan+4] = totW[j][chan]+totW[j][chan+2];
	totW2[j][chan+4] = totW2[j][chan]+totW2[j][chan+2];
      }
    }
    for(int i=0; i<6; i++){
      double n = totW[0][i];
      double N = totW[1][i];
      double ratio = -1;
      if(N!=0) {
	ratio = n/N;
      }
      mean[i] += ratio;
      rms[i] += ratio*ratio;
    }
    genChain->Delete();
  }
  for(int i=0; i<6; i++){
    mean[i] /= (double)nFiles;
    rms[i] = sqrt((rms[i] - mean[i]*mean[i]*nFiles)/(nFiles-1));
    cout<<"Ratio "<<channels[i]<<" is     \t"<<round(mean[i],4)<<" +- "<<round(rms[i],4)<<",\t a "<<
      round(rms[i]/mean[i]*100,2)<<" % error"<<endl;
  }
  cout<<endl;
  cout<<endl;
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
