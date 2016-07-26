#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

TString round(double n, int e, double d=1.);

void CountCombCont(TString CombContfileName = "keys/text/vary_Const.txt", 
		   int index = 5, int nFiles = 1000){
  double nBkg[70], nBkg2[70];
  for(int mc=0; mc<70; mc++){
    nBkg[mc] = 0;
    nBkg2[mc] = 0;
  }
  int component = 0; double constraint;
  fstream CombContFile; CombContFile.open(CombContfileName,fstream::in);
  for(int file=0; file<=nFiles; file++){
    //if(file%200==0) cout<<"Doing iteration "<<file<<" of "<<nFiles<<endl;
    while(component != 40){
      CombContFile>>component>>constraint;
      nBkg[component] += constraint;
      nBkg2[component] += (constraint*constraint);
      //cout<<component<<"\t"<<constraint<<endl;
    }
    component = 0;//cout<<endl;
  }
  nFiles++; //To calculate mean and RMS
  nBkg[index] /= (double)nFiles;
  nBkg2[index] = sqrt((nBkg2[index] - nBkg[index]*nBkg[index]*nFiles)/(double)(nFiles-1));
  cout<<"Constraint "<<index<<" is     \t"<<round(nBkg[index],4)<<" +- "
      <<round(nBkg2[index],4)<<",\t a "<<
    round(nBkg2[index]/nBkg[index]*100,2)<<" % error"<<endl;
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









