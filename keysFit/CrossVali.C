#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TStyle.h"
#include "TSystem.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void crossVali(TString sam="0", TString adaptive = "a", int Times = 20,
	       double wMin = 0.3, double wMax = 3., TString Addf2="withf2") {
//   gSystem->Load("libHtml");
//   gSystem->Load("libMinuit");
//   gSystem->Load("libRooFitCore.so");
//   gSystem->Load("libRooFitModels.so");
//   using namespace RooFit;
  time_t start,end;
  time (&start);
  double dif;

  TString Base = "ND"; Base += adaptive; Base += "_"; Base += sam; Base += "_"; Base += Addf2;

  TString inputfile = "fitSamples/pdfSample"; inputfile += sam; inputfile += ".root";
  cout << "File = " << inputfile << endl;	
  TChain c("ntp1");
  c.Add(inputfile);

  Float_t candM2,candPstarLep;
  c.SetBranchAddress("candM2",&candM2);
  c.SetBranchAddress("candPstarLep",&candPstarLep);
  TCanvas mm("mm","Cross validation output",1200,800);
  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4;
  double entries = c.GetEntries();
  TString hname = "AWG82/results/keys/root/CrossVali/hCross";hname += Base; hname+=".root";
  TFile* hfile = new TFile(hname,"RECREATE"); 

  TH1F h1("cross","Cross Validation",Times, wMin, wMax);

  for(int rep = 0; rep < Times; rep++){
//     RooRealVar mmiss2("candM2","candM2",m2min,m2max);
//     RooRealVar pstarl("candPstarLep","candPstarLep",plmin,plmax);
//     RooArgSet myVars(mmiss2,pstarl);
//     RooDataSet  data("data","data",myVars);
//     if(rep%10==0) cout<<"Doing repetition "<<rep<<" of "<< Times <<endl; 
//     for (int evt = 0 ; evt < entries; evt ++) {
//       c.GetEvent(evt);
//       mmiss2.setVal(candM2);
//       pstarl.setVal(candPstarLep);
//       data.add(RooArgSet(mmiss2,pstarl));
//     }
    double smoo = wMin + rep/(double)(Times-1)*(wMax-wMin);
//     RooNDKeysPdf DPpdf("DPpdf","DPpdf",RooArgList(mmiss2,pstarl),data,adaptive,smoo,4);
    time (&end);dif = difftime (end,start);
    if(rep%50==0) cout<<dif<<" seconds after finding the KEYS function"<<endl;
    time (&start);

    double cv=0;
//     if(Addf2=="withf2") cv = DPpdf.CrossVali(true);
//     else cv = DPpdf.CrossVali(false);
    h1.SetBinContent(rep+1,cv);
    time (&end);dif = difftime (end,start);
    if(Times==1) cout<<dif<<" seconds after Cross validation"<<endl;
    time (&start);
  }
  h1.Draw();
  mm.SaveAs("testCross.eps");
  h1.Write();
//   hfile->Close();
  hfile->Delete();
  cout<<hname<<" saved"<<endl; 
  return;
} 

