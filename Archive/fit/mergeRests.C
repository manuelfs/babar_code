#include "TString.h"
#include "TChain.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void mergeRests(){
  int nFiles[] = {4,14};
  TString relName[] = {"R24","R26"};
  for(int rel=0; rel<2; rel++){
    TString stemName = "AWG82/ntuples/small/"; stemName += relName[rel]; 
    TString fName=stemName; fName+="*_"; fName+=nFiles[rel]; fName+="_RunAll.root";
    TString merName=stemName; merName+="_RunAll.root";
    cout<<"Merging "<<fName<<" into "<<merName<<endl;
    TChain c("ntp1");
    c.Add(fName);
    c.Merge(merName);
    stemName += "Rests_";
    for(int file=1; file<=nFiles[rel]; file=file+2){
      TString f1Name=stemName; f1Name+=file; f1Name+="_"; f1Name+=nFiles[rel]; f1Name+="_RunAll.root";
      TString f2Name=stemName; f2Name+=(file+1); f2Name+="_"; f2Name+=nFiles[rel]; f2Name+="_RunAll.root";
      merName=stemName; merName+=(file+1)/2; merName+="_"; merName+=nFiles[rel]/2; merName+="_RunAll.root";
      cout<<"Merging "<<f1Name<<" and "<<f2Name<<" into "<<merName<<endl;
      TChain c12("ntp1");
      c12.Add(f1Name);c12.Add(f2Name);
      c12.Merge(merName);
    }
  }
  TChain cAll("ntp1");
  cAll.Add("AWG82/ntuples/small/R2*_*4_RunAll.root");
  cAll.Merge("AWG82/ntuples/small/RAll_RunAll.root");
}

