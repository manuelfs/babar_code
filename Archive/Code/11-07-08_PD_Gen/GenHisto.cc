//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: GenHisto.cc,v 1.2 2011/02/22 01:01:50 manuelf Exp $
//
// Description:
//      GenHisto - Creates histograms with the generics distributions to
//                 have "perfect pdf"
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      10/03/08 manuelf -- Created
//------------------------------------------------------------------------

// Creates the 20 pdfSamples and 12 dssSamples cutting in MCType and candType

#include "TROOT.h"
#include "TChain.h"
#include "TH2F.h"
#include "TFile.h"
#include "DonutUtils/KeysUtils.cc"
#include <fstream>
#include <iostream>
using std::cout;
using std::endl;

void saveFile(TChain *inputTree, Int_t thecut, int nbins);

int main(int argc, char *argv[]){
  if (argc < 2 || argc > 3) {
    cout << "USAGE: GenHisto sample=RAll [nbins=250]" << endl;
    return 1;
  }

  TString sample = "RAll";
  if (argc>1) sample = argv[1];
  int nbins = 250;
  if (argc>2) {TString temp = argv[2]; nbins = temp.Atoi();}

  TString genFile = "AWG82/ntuples/small/Fit"; genFile += sample; genFile +="_RunAll.root";
  TChain *tree = new TChain("ntp1");
  tree->Add(genFile);

  for (int i = 1 ; i <= 68 ; i ++) {
    if(i==33)i=41; if(i==49)i=51; if(i==59)i=61; 
    saveFile(tree,i, nbins);
  } 
  cout << "Done with all samples..." << endl;
  
  return 1;
}

void saveFile(TChain *inputTree, Int_t thecut, int nbins){
  bool doPrint = true;
  int mbins = 0, pbins = 0;
  if(thecut<=8)                            {mbins = 50;  pbins = 20;}
  if(thecut>=9&&thecut<=13 || thecut ==15) {mbins = 250; pbins = 22;}
  if(thecut==14 || thecut==16)             {mbins = 100; pbins = 20;}
  if(thecut>=17&&thecut<=20)               {mbins = 40;  pbins = 20;}
  if(thecut>=20&&thecut<=24)               {mbins = 200; pbins = 22;}
  if(thecut>=25&&thecut<=32)               {mbins = 80;  pbins = 20;}
  if(thecut>=40&&thecut<70)                {mbins = 30;  pbins = 20;}
  mbins = nbins; pbins = 100;

  char buf[64];
  sprintf(buf,"keys/root/GenHisto/pdfKeys_%d_Fit.root",thecut);
  TFile *f = new TFile(buf,"RECREATE");
  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4;
  TH2F h2("h2","KEYS",mbins,m2min,m2max,pbins,plmin,plmax);  

  Float_t weight, candM2, candPstarLep;
  Int_t MCType, candType;
  inputTree->SetBranchAddress("weight",&weight);
  inputTree->SetBranchAddress("MCType",&MCType);
  inputTree->SetBranchAddress("candM2",&candM2);
  inputTree->SetBranchAddress("candPstarLep",&candPstarLep);
  inputTree->SetBranchAddress("candType",&candType);
  double nl = 0;
  Int_t nevt = (Int_t)inputTree -> GetEntries();
  for (int j = 0 ; j < nevt ; j ++) {
    inputTree -> GetEvent(j);
    if(!isSample(thecut, MCType, candType)) continue;

    h2.Fill(candM2,candPstarLep,weight);
    nl+=weight;
  }
  if(doPrint) {
    cout<<thecut<<"\t"<<RoundNumber(nl,1);
    cout<<endl;
  }
  for(int i=1; i<mbins+1; i++)
    for(int j=1; j<pbins+1; j++)
      if(h2.GetBinContent(i,j)<0.000001) h2.SetBinContent(i,j,0.000001);
  h2.Write();
  f->Close(); f->Delete();
}
