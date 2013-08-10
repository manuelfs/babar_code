//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: FitGenSample.cc,v 1.4 2011/02/22 01:01:50 manuelf Exp $
//
// Description:
//      FitGenSample - Fits a generic sample from the truth-matched and
//      non-truth-matched components in cocktail
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      10/02/10 manuelf -- Created
//------------------------------------------------------------------------

#include "TString.h"
#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TTree.h"
#include "TSystem.h"
#include "TIterator.h"
#include "TRandom3.h"
#include "RooFitCore/RooFitResult.hh"
#include "RooFitCore/RooArgList.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooFormulaVar.hh"
#include "RooFitCore/RooAddPdf.hh"
#include "RooFitCore/RooAbsReal.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooDataHist.hh"
#include "RooFitCore/RooHistPdf.hh"
#include "RooFitCore/RooPlot.hh"
#include "RooFitCore/RooNLLVar.hh"
#include "RooFitCore/RooMinuit.hh"
#include "RooFitCore/RooGlobalFunc.hh"
#include "DonutUtils/KeysUtils.cc"
#include <fstream>
#include <iostream>

using namespace RooFit;
using namespace std;
using std::cout;
using std::endl;

int main(int argc, char *argv[]){
  if (argc < 3 || argc > 7 ) {
    cout << "USAGE: FitGenSample HistoTru HistoNonTru [dataFolder = Rest2noSep] [Times=1] [HistoTau=""] [fixTau]" << endl;
    return 0;
  }
  // Setting the input parameters
  TString hnameTru = argv[1];
  TString hnameNonTru = argv[2];
  TString dataFolder = "Rest2noSep";
  if(argc>3) dataFolder  = argv[3]; 
  int Times = 1;
  if(argc>4)  {TString sTimes = argv[4]; Times = sTimes.Atoi();}
  TString hnameTau = "";
  if(argc>5) hnameTau = argv[5];
  double fixTau = -1;
  if(argc>6) {TString temp = argv[6]; fixTau = temp.Atof();}
  
  TString sam, smoo, scaleP, Base;
  Int_t Sam, nM2bin, nPlbin;
  ParseName(hnameTru, sam, smoo, scaleP, Base, Sam, nM2bin, nPlbin);
  double entries = 0.;
  TString inputfile = nameData(Sam,dataFolder);
  TChain *treeData = new TChain("ntp1");
  treeData->Add(inputfile);
  TString Totname = "TotIntegral"; 
  TString Totvari = "candPstarLep>>"; Totvari+=Totname; 
  TH1F *TotHisto = new TH1F(Totname, "", 10, 0, 2.4);
  if(treeData) treeData->Draw(Totvari,"weight");
  entries = TotHisto->Integral();
  cout<<"Fitting "<<inputfile<<endl;

  TString hnameFit = "AWG82/results/keys/root/";
  if(Times==1) hnameFit += "GenFit/pdfKeys_"; 
   else hnameFit += "TimesGenFit/hTimes_"; 
  hnameFit += Sam; hnameFit += "_Fit.root";
  TFile* hfileFit = new TFile(hnameFit,"RECREATE"); 
  hfileFit->Close(); hfileFit->Delete();
  for(int rep = 0; rep < Times; rep++){
    if(rep%50==0) cout<<"Repetition "<<rep+1<<" of "<< Times <<endl; 
    // Setting the fitting parameters
    bool Histo3 = false;
    if(hnameTau != "") Histo3 = true;
    TString h2Name = "h2";
    if(Times>1){ h2Name = "hTimes"; h2Name += rep;}
    RooRealVar mmiss2("candM2","candM2",-4,12);
    RooRealVar pstarl("candPstarLep","candPstarLep",0.,2.4);
    RooRealVar weight("weight","weight",0.,100.);
    RooDataSet data("data","data",treeData,RooArgSet(mmiss2,pstarl,weight),"","weight");
    TFile hfileTru(hnameTru); 
    TH2F *h2Tru = (TH2F *)(hfileTru.Get(h2Name));    
    if(h2Tru==0){
      cout<<h2Name<<" does not exit in "<<hnameTru<<endl;
      return 0;
    }
    h2Tru->SetDirectory(0);
    RooDataHist RDTru("RDTru", "RDTru", RooArgList(mmiss2,pstarl), h2Tru);
    RooHistPdf pdfTru("pdfTru","pdfTru",RooArgList(mmiss2,pstarl),RDTru,2);
    TFile hfileNonTru(hnameNonTru); hfileNonTru.cd();
    TH2F *h2NonTru = (TH2F *)(hfileNonTru.Get(h2Name));    
    if(h2NonTru==0){
      cout<<h2Name<<" does not exit in "<<hnameNonTru<<endl;
      return 0;
    }
    h2NonTru->SetDirectory(0);
    RooDataHist RDNonTru("RDNonTru", "RDNonTru", RooArgList(mmiss2,pstarl), h2NonTru);
    RooHistPdf pdfNonTru("pdfNonTru","pdfNonTru",RooArgList(mmiss2,pstarl),RDNonTru,2);
    RooDataHist *RDTau=0;
    RooHistPdf *pdfTau=0;
    TH2F *h2Tau=0;
    if(Histo3){
      TFile hfileTau(hnameTau); 
      h2Tau = (TH2F *)(hfileTau.Get(h2Name));    
      if(h2Tau==0){
	cout<<h2Name<<" does not exit in "<<hnameTau<<endl;
	return 0;
      }
      h2Tau->SetDirectory(0);
      RDTau = new RooDataHist("RDTau", "RDTau", RooArgList(mmiss2,pstarl), h2Tau);
      pdfTau = new RooHistPdf("pdfTau","pdfTau",RooArgList(mmiss2,pstarl),*RDTau,2);
    }

    RooRealVar *fNonTru = new RooRealVar("fNonTru","fNonTru",0.2,0,1.);
    RooRealVar *fTau=0;
    RooFormulaVar *fTru;
    if(Histo3) {
      fTau = new RooRealVar("fTau","fTau",0.1,0,1.);
      fTru = new RooFormulaVar("fTru","fTru","1-fTau-fNonTru",RooArgList(*fTau,*fNonTru));
    } else fTru = new RooFormulaVar("fTru","fTru","1-fNonTru",RooArgList(*fNonTru));
    RooAddPdf *pdfSum;
    if(!Histo3) pdfSum = new RooAddPdf("pdfSum","pdfSum",RooArgList(pdfTru,pdfNonTru),RooArgList(*fTru,*fNonTru));
    else pdfSum = new RooAddPdf("pdfSum","pdfSum",RooArgList(pdfTru,pdfNonTru,*pdfTau),RooArgList(*fTru,*fNonTru,*fTau));

    // Fitting
    RooNLLVar nll("nll","nll",*pdfSum,data,kFALSE);
    RooMinuit min(nll);
    if(Times>1){min.setPrintLevel(-1); min.setVerbose(kFALSE);}
    min.setStrategy(2);
    min.hesse(); min.migrad();

    if(Times==1){
      fstream os;
      TString textName = "keys/text/GenFit/Fractions_"; textName += sam; textName += ".txt";
      os.open(textName,fstream::out);
      if(Histo3){
	RooFitResult *fitres = min.save();
	RooArgList lis = fitres->randomizePars();
	TIterator *iter= lis.createIterator();
	RooRealVar *par(0);
	for(int i=0; i<10000; i++){
	  fitres->randomizePars();
	  iter= lis.createIterator();
	  while((0 != (par= (RooRealVar*)iter->Next()))) {
	    double genR = par->getVal();
	    if(genR>=1) genR=1-1e-6;
	    if(genR<=0) genR=1e-6;
	    os<<genR<<"\t";
	  }
	  os<<endl;
	}
	delete iter;
      } else {
	double mean = fNonTru->getVal(), rms = fNonTru->getError();
	TRandom3 rand(0);
	for(int i=0; i<10000; i++) {
	  double genR = rand.Gaus(mean,rms);
	  while(genR>=1 || genR<=0) genR = rand.Gaus(mean,rms);
	  os<<genR<<endl;    
	}
      }
    }

    // Adding histograms and saving
    double ffitTru = fTru->getVal(), ffitNonTru = fNonTru->getVal();
    double scaleNonTru = h2Tru->Integral()/h2NonTru->Integral()*ffitNonTru/ffitTru;
    if(Histo3){
      double ffitTau = fTau->getVal(); 
      if(fixTau>0) ffitTau = fixTau*(ffitTru+ffitNonTru)/(entries-fixTau);
      double scaleTau = h2Tru->Integral()/h2Tau->Integral()*ffitTau/ffitTru;
      h2Tru->Add(h2Tau, scaleTau);
    } 
    h2Tru->Add(h2NonTru, scaleNonTru);
    hfileFit = new TFile(hnameFit,"UPDATE"); 
    hfileFit->cd();
    h2Tru->Write(h2Name);
    hfileFit->Close(); hfileFit->Delete();
    delete pdfSum;
    delete pdfTau;
    //RDTau->Delete();fTau->Delete();fTru->Delete();
    //h2Tau->Delete();h2Tru->Delete();h2NonTru->Delete();
  }
  cout<<"Saved the sum in "<<hnameFit<<endl;

  // Plotting
  if(Times==1) PlotEps(hnameFit, dataFolder, "GenFit", "1", "EPS");
}

//   RooRealVar nTru("nTru","nTru",entries*0.7,0,entries);
//   RooRealVar *nNonTru = new RooRealVar("nNonTru","nNonTru",entries*0.2,0,entries);
//   RooRealVar *nTau=0;
//   if(Histo3) {
//     nTau = new RooRealVar("nTau","nTau",entries*0.1,0,entries);
//   }
//   RooAddPdf *pdfSum;
//   if(!Histo3) pdfSum = new RooAddPdf("pdfSum","pdfSum",RooArgList(pdfTru,pdfNonTru),RooArgList(nTru,*nNonTru));
//   else pdfSum = new RooAddPdf("pdfSum","pdfSum",RooArgList(pdfTru,pdfNonTru,*pdfTau),RooArgList(nTru,*nNonTru,*nTau));

//   // Fitting
//   RooNLLVar nll("nll","nll",*pdfSum,data,kTRUE);
//   RooMinuit min(nll);
//   min.setStrategy(2);
//   min.hesse(); min.migrad();
//   min.hesse(); min.minos();

//   // Adding histograms and saving
//   double nfitTru = nTru.getVal(), nfitNonTru = nNonTru->getVal();
//   double scaleNonTru = h2Tru->Integral()/h2NonTru->Integral()*nfitNonTru/nfitTru;
//   if(Histo3){
//     double nfitTau = nTau->getVal(); //nfitTau = entries - nfitTru - nfitNonTru;
//     if(fixTau>0) nfitTau = fixTau*(nfitTru+nfitNonTru)/(entries-fixTau);
//     double nTotal = nfitTau+nfitNonTru+nfitTru;
//     cout<<"Sample "<<Sam<<" uses fNon = "<<nfitNonTru/nTotal<<", and fTau = "<<nfitTau/nTotal<<endl;
//     double scaleTau = h2Tru->Integral()/h2Tau->Integral()*nfitTau/nfitTru;
//     h2Tru->Add(h2Tau, scaleTau);
//   } else {
//     double nTotal = nfitNonTru+nfitTru;
//     cout<<"Sample "<<Sam<<" uses fNon = "<<nfitNonTru/nTotal<<endl;
//   }
//   h2Tru->Add(h2NonTru, scaleNonTru);
//   TString hnameFit = "AWG82/results/keys/root/GenFit/pdfKeys_"; hnameFit += Sam; hnameFit += "_Fit.root";
//   TFile* hfileFit = new TFile(hnameFit,"RECREATE"); 
//   hfileFit->cd();
//   h2Tru->Write("h2");
//   hfileFit->Close(); hfileFit->Delete();
//   cout<<"Saved the sum in "<<hnameFit<<endl;

//   // Plotting
//   PlotEps(hnameFit, dataFolder, "GenFit", weightName, "EPS");
// }


