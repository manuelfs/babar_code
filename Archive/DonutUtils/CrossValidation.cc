//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: CrossValidation.cc,v 1.5 2012/03/04 00:33:02 manuelf Exp $
//
// Description:
//      CrossValidation - Calculates the crossvalidation algorithme on 
//      several KEYS fits
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      10/02/12 manuelf -- Cleaned up and made it 2D by looping over the
//               smoothing of p*l as well
//      09/12/08 manuelf -- Created
//------------------------------------------------------------------------

#include "TString.h"
#include "TFile.h"
#include "TH2F.h"
#include "TChain.h"
#include "TStyle.h"
#include "TSystem.h"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooAddPdf.hh"
#include "RooFitCore/RooAbsReal.hh"
#include "RooFitCore/RooDataSet.hh"
#include "DonutUtils/RooNDKeysPdfDonut.hh"
#include "DonutUtils/weightManager.hh"
#include "DonutUtils/KeysUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using namespace RooFit;
using std::cout;
using std::endl;

int main(int argc, char *argv[]){

  if (argc < 2 || argc > 11 ) {
    cout << "USAGE: CrossValidation sample [mTimes=20] [mMin=30] [mMax=200] [pTimes=5] [pMin=100] [pMax=300]"
	 << "[doRot=noRot] [sample=RAll] [nSigma=3]" << endl;
    return 1;
  }
  TString sam = argv[1]; int Sam = sam.Atoi();
  int mTimes = 20;
  if(argc>2)  {TString mTimes_s = argv[2]; mTimes = mTimes_s.Atoi();}
  TString mMin_s = "30";
  if(argc>3)  mMin_s = argv[3]; 
  double mMin = mMin_s.Atof()/100.;
  TString mMax_s = "200";
  if(argc>4)  mMax_s = argv[4]; 
  double mMax = mMax_s.Atof()/100.;
  int pTimes = 5;
  if(argc>5)  {TString pTimes_s = argv[5]; pTimes = pTimes_s.Atoi();}
  TString pMin_s = "100";
  if(argc>6)  pMin_s = argv[6]; 
  double pMin = pMin_s.Atof()/100.;
  TString pMax_s = "300";
  if(argc>7)  pMax_s = argv[7]; 
  double pMax = pMax_s.Atof()/100.;
  bool doRot = false;
  if(argc>8)  {TString t_s = argv[8]; if(t_s == "Rot") doRot = true;}
  if(doRot) cout<<"Rotating data"<<endl;
  TString Sample = "RAll";
  if (argc>9) Sample = argv[9];
  double nSigma = 3.;
  if(argc>10)  {TString sTemp = argv[10]; nSigma = sTemp.Atof();}
  if(mTimes == 1) mTimes++; if(pTimes == 1) pTimes++; 

  TString adaptive = "a";
  TString Addf2 = "nof2";
  TString Base = "_"; Base += sam; Base += "_"; Base += mMin_s; Base += "to"; Base += mMax_s;
  Base += "_"; Base += pMin_s; Base += "to"; Base += pMax_s;

  TString inputfile = nameData(Sam,Sample);
  TChain tree("ntp1");
  tree.Add(inputfile);
  float candM2,candPstarLep, weight=-1;
  tree.SetBranchAddress("candM2",&candM2);
  tree.SetBranchAddress("candPstarLep",&candPstarLep);
  tree.SetBranchAddress("weight",&weight);

  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4;
  double entries = tree.GetEntries();
  RooRealVar mmiss2("candM2","candM2",m2min,m2max);
  RooRealVar pstarl("candPstarLep","candPstarLep",plmin,plmax);
  RooRealVar totWeight("totWeight","totWeight",0.,100.);
  RooArgSet myVars(mmiss2,pstarl,totWeight);
  RooDataSet  data("data","data",myVars);
  for (int evt = 0 ; evt < entries; evt ++) {
    tree.GetEvent(evt);
    mmiss2.setVal(candM2);
    pstarl.setVal(candPstarLep);
    totWeight.setVal(weight);
    if(candM2<m2min||candM2>m2max || candPstarLep<plmin||candPstarLep>plmax || weight<0||weight>100) continue;
    data.add(RooArgSet(mmiss2,pstarl,totWeight));
  }
  data.setWeightVar(totWeight);

  TString hname = "AWG82/results/keys/root/CrossVali/hCrossKeys";hname += Base; hname+=".root";
  TFile* hfile = new TFile(hname,"RECREATE"); 
  hfile->Close(); hfile->Delete();
  double bin_m = (mMax-mMin)/((double)(mTimes-1)), minCross = 9999.0, minSmooth_m = mMin, cv = 0;
  double bin_p = (pMax-pMin)/((double)(pTimes-1)), smoo_p = pMin, minSmooth_p = pMin;
  TH2F h2("cross","Cross Validation", mTimes, mMin-bin_m/2, mMax+bin_m/2, pTimes, pMin-bin_p/2, pMax+bin_p/2);
  for(int rep_p = 0; rep_p < pTimes; rep_p++){
    double smoo_m = mMin;
    for(int rep_m = 0; rep_m < mTimes; rep_m++){
      //RooNDKeysPdfDonut DPpdf("DPpdf","DPpdf",RooArgList(mmiss2,pstarl),data,adaptive,smoo_m,smoo_p*smoo_m,3.,doRot);
      RooNDKeysPdfDonut DPpdf("DPpdf","DPpdf",RooArgList(mmiss2,pstarl),data,adaptive,smoo_m,smoo_p,3.,doRot);
      if(Addf2=="withf2") cv = DPpdf.CrossVali(true,nSigma);
      else cv = DPpdf.CrossVali(false,nSigma);
      h2.SetBinContent(rep_m+1,rep_p+1,cv);
      hfile = new TFile(hname,"RECREATE"); 
      h2.Write();
      hfile->Close(); hfile->Delete();
      if(cv<minCross) {
	minCross=cv; minSmooth_m=smoo_m; minSmooth_p=smoo_p;
	cout<<rep_m+1<<"/"<< mTimes <<", "<< rep_p+1<<"/"<< pTimes <<":\t CV "<<RoundNumber(cv*100,3)<<
	  " at mmiss smoothing "<<RoundNumber(smoo_m,3)<<" and scaleP "<<RoundNumber(smoo_p,3)<<endl; 
      }
      smoo_m += bin_m;
    }
    cout<<"Finished "<<rep_p+1<<"/"<< pTimes <<":\t CV "<<RoundNumber(cv*100,2)<<
      " at mmiss smoothing "<<RoundNumber(smoo_m-bin_m,3)<<" and scaleP "<<RoundNumber(smoo_p,3)<<endl; 
    smoo_p += bin_p;
  }

  cout<<"Saved "<<hname<<endl;
  cout<<endl<<endl<<"Minimum CrossV "<<RoundNumber(minCross*100,3)<<" at smoo_m = "<<RoundNumber(minSmooth_m,3)
      <<" and smoo_p = "<<RoundNumber(minSmooth_p,3)<<endl<<endl; 
  return 1;
}

