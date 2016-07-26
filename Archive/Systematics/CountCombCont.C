#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

TString round(double n, int e, double d=1.);

void CountCombCont(TString CombContfileName = "AWG82/systematics/BF/CombContYields.txt", int nFiles = 1000){
  double nBkg[2][8], nBkg2[2][8];
  for(int mc=0; mc<2; mc++){
    for(int cand=0; cand<8; cand++){
      nBkg[mc][cand] = 0;
      nBkg2[mc][cand] = 0;
    }
  }
  fstream CombContFile; CombContFile.open(CombContfileName,fstream::in);
  for(int file=0; file<=nFiles; file++){
    //if(file%200==0) cout<<"Doing iteration "<<file<<" of "<<nFiles<<endl;
    for(int mc=0; mc>-2; mc--){
      for(int cand=1; cand<=8; cand++)	{
	double Entries = 0;
	CombContFile>>Entries;
	nBkg[mc+1][cand-1] += Entries;
	nBkg2[mc+1][cand-1] += (Entries*Entries);
      }
    }
  }
  nFiles++; //To calculate mean and RMS
  for(int mc=1; mc>=0; mc--){
    for(int cand=0; cand<8; cand++){
      nBkg[mc][cand] /= (double)nFiles;
      nBkg2[mc][cand] = sqrt((nBkg2[mc][cand] - nBkg[mc][cand]*nBkg[mc][cand]*nFiles)/(double)(nFiles-1));
      cout<<"Yield "<<51+cand+(1-mc)*10<<" is "<<round(nBkg[mc][cand],1)<<" +- "
	  <<round(nBkg2[mc][cand],1)<<",\t a "<<
	round(nBkg2[mc][cand]/nBkg[mc][cand]*100,2)<<" % error"<<endl;
    }
    cout<<endl;
  }
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









