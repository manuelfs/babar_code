//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DonutFitFinal.cc,v 1.3 2011/02/22 01:01:50 manuelf Exp $
//
// Description:
//      DonutFitFinal - Global fit to B->D(*)taunu with histograms
//                      Includes code to run toys
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      10/03/18 manuelf -- Adapted from Mazur's DonutFitWrapperKeys.cc 
//------------------------------------------------------------------------

#include "TFile.h"
#include "TChain.h"
#include "TH2.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooDataHist.hh"
#include "RooFitCore/RooHistPdf.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooFitResult.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooAddPdf.hh"
#include "RooFitCore/RooStringVar.hh"
#include "RooFitCore/RooCategory.hh"
#include "RooFitCore/RooSimultaneous.hh"
#include "RooFitCore/RooNLLVar.hh"
#include "RooFitCore/RooFormulaVar.hh"
#include "RooFitCore/RooMinuit.hh"
#include "DonutUtils/KeysUtils.cc"
#include <iomanip>
#include <iostream>

using std::cout;
using std::endl;
using std::fixed;
using std::setprecision;
using namespace RooFit;

void initializeYields();
void initializePdfs(int rep=1);
void readConfigFile(TString configFile);
Bool_t parseToBool(RooStringVar &rsv, Bool_t &val);
void readDataset();
void generateToy();
Bool_t isFloated(int channel);
int doChannel(int channel);
void printAsLatex();
void printAsText();
void printYieldsChannel(ostream& os, int channel, double &dF, double &dT);
void printCorrelMat(ostream& os, int isDc);


// Fitter options
Bool_t FixYield[70];
TString rootfile, latexfile, plotfile, textfile, pullfile, constfile, inifile, TypeFit;

// Various
double m2min, m2max, plmin, plmax, truYield[70], truYieldFix[70], nData[8];
int IndfYield[70], typeToy=0, nhTimes[70], floatDinDss=1, IsoConst=0, RatioFit=0, isData=0;
int IndPdf[2][9] = {{1,9,5,13,17,41,51,61,37},{21,25,29,45,55,65,33,-1,-1}};
TString candLab[8] = {"D0","Ds0","Dc","Dsc","DssD0","DssDs0","DssDc","DssDsc"};
TChain *treeData;

// Fit variables
RooRealVar *mmiss2, *pstarl, *totWeight;
RooCategory *candType0,*candTypec;
RooAbsReal *Yield[70];
RooRealVar *fYield[70], *RatioTau[5];
RooDataHist* RDpdf[70];
RooHistPdf* pdf[70];
RooAddPdf *pdfsum[8];
RooSimultaneous *sim0pdf, *simcpdf;
RooDataSet *data0,*datac;
RooNLLVar *nll0, *nllc;
RooMinuit *min0, *minc;
RooFitResult *fitres0,*fitresc;
RooDataSet *protoDS[70];

int main(int argc, char *argv[]){
  if (argc < 2 || argc > 6) {
    cout << "USAGE: DonutFitFinal configfile [typeToy=0] [Times=1] [Tag=R24Rests_1_2] [TagOut=""]" << endl;
    cout << "typeToy: -1 Write tables, 0 Normal fit, 1 Standard, 2 Bootstrap"<<
      ", 3 D** constraint, 4 Constraints Error, 5 KEYS Error" << endl;
    return 0;
  }

  TString configName = argv[1];
  if(argc>2)  {TString sTemp = argv[2];  typeToy = sTemp.Atoi();}
  int Times = 1;
  if(argc>3)  {TString sTimes = argv[3]; Times = sTimes.Atoi();}
  TString Tag = "R24Rests_1_2";
  if(argc>4)  Tag = argv[4];
  TString TagOut = "";
  if(argc>5)  TagOut = argv[5];

  if(typeToy == 0) Times=1;
  cout<<"Reading "<<configName<<endl;
  readConfigFile(configName);
  if(rootfile==""){cout<<configName<<" not properly read"<<endl; return 0;}
  if(!TypeFit.Contains("D0") && !TypeFit.Contains("Dc") && !TypeFit.Contains("Iso") && !TypeFit.Contains("RatioFit")) {
    cout<<"TypeFit has to be (at least) one of: D0, Dc, Iso, RatioFit"<<endl;
    return 0;
  }
  rootfile   = "AWG82/ntuples/small/Fit"; rootfile += Tag; rootfile += "_RunAll.root";
  latexfile  = "FitAll/fits/TableFinal"; latexfile += Tag; latexfile += TagOut; latexfile += ".tex";
  plotfile   = "FitAll/fits/PlotsFinal"; plotfile += Tag; plotfile += TagOut; plotfile += ".ps";
  textfile   = "FitAll/fits/TextFinal"; textfile += Tag; textfile += TagOut; textfile += ".txt";
  pullfile   = "FitAll/fits/Pulls"; pullfile += typeToy; pullfile += Tag; 
  pullfile += TagOut; pullfile += ".root";
  TString textName = "keys/text/vary_RestsConst.txt";
  TString DssConstName = "keys/text/DssConst.txt";
  double EffRatio[] = {0.747, 0.476, 0.768, 0.440}, DssConst[2][4];

  if(TypeFit.Contains("Iso")) IsoConst = 1;          if(TypeFit.Contains("IsoTau")) IsoConst = 2;
  if(TypeFit.Contains("fixDinDss")) floatDinDss = 0; if(TypeFit.Contains("RatioFit")) RatioFit=1;
  if(rootfile.Contains("Data")) isData = 1;

  treeData = new TChain("ntp1");
  treeData->Add(rootfile);  
  initializeYields();
  if(typeToy == -1){
    readDataset();
    cout<<constfile<<" and "<<inifile<<" saved"<<endl;
    return 1;
  }
  if(typeToy<2 || typeToy >2) readDataset();
  // Initializing histograms for toys and error measurement
  TH1F *hPull[70], *hYield[70], *hRatio[5];
  if(Times>1){
    for(int chan=0; chan<4; chan++){
      TString hRatioName = "Ratio"; hRatioName += chan+1;
      double nTau = truYield[IndPdf[0][0]+chan], nL = truYield[IndPdf[0][1]+chan];
      if(chan%2==1) {
	nTau += truYield[IndPdf[0][2]+chan-1];
	nL   += truYield[IndPdf[0][3]+chan-1];
      } else {
	nTau += truYield[IndPdf[0][2]+chan+1];
	nL   += truYield[IndPdf[0][3]+chan+1];
      }
      double Ratio = nTau/nL/EffRatio[chan]*2;
      hRatio[chan] = new TH1F(hRatioName,hRatioName,300,0.6*Ratio,1.5*Ratio);
    }
    for(int chan=0; chan<8; chan++){
      if(!doChannel(chan)) continue;
      for(int i=0; i<9; i++){
	if(IndPdf[chan>3][i]<0) continue;
	int index = IndPdf[chan>3][i]+chan-4*(chan>3);
	if(isFloated(index) && FixYield[index]==0){
	  TString hPullName = "Pull"; hPullName+=index;
	  hPull[index] = new TH1F(hPullName,hPullName,200,-4,4);
	  TString hYieldName = "YieldFit"; hYieldName+=index;
	  hYield[index] = new TH1F(hYieldName,hYieldName,300,0.6*truYield[index],1.5*truYield[index]);
	}
      }
    }
  }
  fstream varyCons, fileDssConst;
  if(typeToy == 3) {
    fileDssConst.open(DssConstName,fstream::in);
    for(int chan=0; chan<4; chan++){
      fileDssConst>>DssConst[0][chan];
      DssConst[1][chan] = fYield[IndPdf[0][4]+chan]->getVal();
    }
  }
  if(typeToy == 4) varyCons.open(textName,fstream::in);
  if(typeToy == 5) for(int i=0; i<70; i++) nhTimes[i] = 0;
  cout<<"Reading pdfs and initializing Fit with typeToy "<<typeToy<<endl;
  if(typeToy!=5) initializePdfs();

  // Loop for toys
  for(int rep = 0; rep < Times; rep++){
    cout<<"Repetition "<<rep+1<<" of "<< Times <<endl; 
    if(nll0)nll0->Delete();if(min0)min0->Delete();
    if(typeToy==1) generateToy();
    if(typeToy==2) readDataset();
    if(typeToy == 3) {
      TRandom3 rand(0); 
      for(int chan=0; chan<4; chan++){
	fYield[IndPdf[0][4]+chan]->setVal(rand.Gaus(DssConst[1][chan],DssConst[1][chan]*DssConst[0][chan]));
      }
    }
    if(typeToy == 4) {
      int findex;
      double fraction;
      for(int i=0; i<9; i++){
	for(int chan=0; chan<8; chan++){
	  if(IndPdf[chan>3][i]<0) continue;
	  int index = IndPdf[chan>3][i]+chan-4*(chan>3);
	  if(!isFloated(index)) {
	    varyCons>>findex>>fraction;
	    if(index==findex) fYield[findex]->setVal(fraction);
	    else {
	      cout<<textName<<" does not have the appropriate format"<<endl;
	      return 0;
	    }
	  }
	}
      }
    }
    if(typeToy==5) initializePdfs(rep);
    Bool_t doExt = kTRUE;
    if(TypeFit.Contains("D0") || IsoConst){
      if(IsoConst==0) cout<<endl<<"Fitting the neutral modes"<<endl<<endl;
      else cout<<endl<<"Fitting the neutral and charged modes with Isospin constrained"<<endl<<endl;
      nll0 = new RooNLLVar("nll0","nll0",*sim0pdf,*data0,doExt,0,0,1,kFALSE,kFALSE);
      min0 = new RooMinuit(*nll0); 
      min0->setPrintLevel(-1); min0->setVerbose(kFALSE);
      min0->setStrategy(2);  min0->optimizeConst(kTRUE);
      min0 -> hesse();       min0 -> migrad();
      if(Times>1){
	for(int chan=0; chan<8; chan++){
	  if(!doChannel(chan) || (chan%4)>1) continue;
	  for(int i=0; i<9; i++){
	    if(IndPdf[chan>3][i]<0) continue;
	    int index = IndPdf[chan>3][i]+chan-4*(chan>3);
	    if(isFloated(index) && FixYield[index]==0){
	      double FittedY = Yield[index]->getVal();
	      double pull;
	      if(RatioFit && index<=4){
		pull = (RatioTau[index]->getVal()-truYieldFix[index]/truYieldFix[index+8])/
		  (static_cast <RooRealVar*>(RatioTau[index]))->getError();
	      } else pull = (FittedY-truYieldFix[index])/(static_cast <RooRealVar*>(Yield[index]))->getError();
	      hPull[index]->Fill(pull);
	      hYield[index]->Fill(FittedY);
	    }
	  }
	}
      }
      fitres0 = min0-> save();
    }
    if(nllc)nllc->Delete();if(minc)minc->Delete();
    if(TypeFit.Contains("Dc")) {
      cout<<endl<<"Fitting the charged modes"<<endl<<endl;
      nllc = new RooNLLVar("nllc","nllc",*simcpdf,*datac,doExt);
      minc = new RooMinuit(*nllc); 
      minc->setPrintLevel(-1); minc->setVerbose(kFALSE);
      minc->setStrategy(2);  minc->optimizeConst(kTRUE);
      minc -> hesse();       minc -> migrad();
      if(Times>1){
	for(int chan=0; chan<8; chan++){
	  if(!doChannel(chan) || (chan%4)<2) continue;
	  for(int i=0; i<9; i++){
	    if(IndPdf[chan>3][i]<0) continue;
	    int index = IndPdf[chan>3][i]+chan-4*(chan>3);
	    if(isFloated(index) && FixYield[index]==0){
	      double FittedY = Yield[index]->getVal();
	      double pull;
	      if(RatioFit && index<=4){
		pull = (RatioTau[index]->getVal()-truYieldFix[index]/truYieldFix[index+8])/
		  (static_cast <RooRealVar*>(RatioTau[index]))->getError();
	      } else pull = (FittedY-truYieldFix[index])/(static_cast <RooRealVar*>(Yield[index]))->getError();
	      hPull[index]->Fill(pull);
	      hYield[index]->Fill(FittedY);
	    }
	  }
	}
      } 
      fitresc = minc-> save();
    }
    if(Times>1){ // Filling the histograms
      TFile fPull(pullfile,"RECREATE"); 
      for(int chan=0; chan<4; chan++) {
	double nTau = Yield[IndPdf[0][0]+chan]->getVal(), nL = Yield[IndPdf[0][1]+chan]->getVal();
	if(chan%2==1) {
	  nTau += Yield[IndPdf[0][2]+chan-1]->getVal();
	  nL   += Yield[IndPdf[0][3]+chan-1]->getVal();
	} else {
	  nTau += Yield[IndPdf[0][2]+chan+1]->getVal();
	  nL   += Yield[IndPdf[0][3]+chan+1]->getVal();
	}
	double Ratio = nTau/nL/EffRatio[chan]*2;
	hRatio[chan]->Fill(Ratio);
	hRatio[chan]->Write();
      }
      for(int chan=0; chan<8; chan++){
	if(!doChannel(chan)) continue;
	for(int i=0; i<9; i++){
	  if(IndPdf[chan>3][i]<0) continue;
	  int index = IndPdf[chan>3][i]+chan-4*(chan>3);
	  if(isFloated(index) && FixYield[index]==0){
	    hPull[index]->Write();
	    hYield[index]->Write();
	  }
	}
      }
      fPull.Close();
    }
  }
  if(Times==1){    
    printAsLatex();
    printAsText();
  }
  cout<<endl<<endl;
  RooFitResult *fitres[2] = {fitres0, fitresc};
  for(int chan=0; chan<4; chan++){
    double errorRatio, Ratio, errorTau;
    int index = IndPdf[0][0]+chan, indexL = IndPdf[0][1]+chan;
    if(Times==1){      
      double yieldL = Yield[indexL]->getVal() + Yield[indexL+5-2*(chan%2)]->getVal();
      double errorL = (static_cast <RooRealVar*>(Yield[indexL]))->getError();
      if(chan%2==1) errorL = sqrt(pow(errorL,2)+pow((static_cast <RooRealVar*>(Yield[indexL+3]))->getError(),2));
      else errorL *= yieldL/Yield[indexL]->getVal();
      if(RatioFit) {
	Ratio = RatioTau[index]->getVal()/EffRatio[chan]*2;
	errorRatio = (static_cast <RooRealVar*>(RatioTau[index]))->getError()/EffRatio[chan]*2;
	errorTau = sqrt(pow(yieldL*errorRatio,2)+pow(RatioTau[index]->getVal()*errorL,2)+
			fitres[index>2]->correlation(*Yield[indexL],*RatioTau[index])*errorRatio*errorL*Ratio*yieldL);
      } else {
	double yieldTau = Yield[index]->getVal() + Yield[index+5-2*(chan%2)]->getVal();
	errorTau = (static_cast <RooRealVar*>(Yield[index]))->getError()*yieldTau/Yield[index]->getVal();
	Ratio = yieldTau/yieldL/EffRatio[chan]*2;
	errorRatio = Ratio*sqrt(pow(errorTau/yieldTau,2) + pow(errorL/yieldL,2));
      }
    } else {
      Ratio = hRatio[chan]->GetMean();
      errorRatio = hRatio[chan]->GetRMS();
    }
    cout<<"Ratio "<<index<<": "<<RoundNumber(Ratio,4)<<" +- "<< RoundNumber(errorRatio,4)<<
      ", a "<<RoundNumber(errorRatio*100,1,Ratio)<<" % error"<<endl;
  }
  cout<<endl<<endl;
  delete mmiss2; delete pstarl;
  for (int i = 0 ; i < 70 ; i ++) {
    if(pdf[i])    pdf[i]   ->Delete(); 
    if(RDpdf[i])  RDpdf[i] ->Delete(); 
    if(Yield[i])  Yield[i] ->Delete(); 
    if(fYield[i]) fYield[i]->Delete(); 
  }
  treeData->Delete();
}


// Yields are initialized
void initializeYields() {
  m2min = -4; m2max = 12; plmin = 0; plmax = 2.4;
  mmiss2 = new RooRealVar("candM2","m_{miss}^{2}",m2min,m2max,"GeV^{2}/c^{4}");
  pstarl = new RooRealVar("candPstarLep","lepton momentum",plmin,plmax,"GeV/c");

  for(int chan=0; chan<8; chan++){
    if(!doChannel(chan)) continue;
    for(int i=0; i<9; i++){
      if(IndPdf[chan>3][i]<0) continue;
      int index = IndPdf[chan>3][i]+chan-4*(chan>3);
      if(isFloated(index)){
	TString yieldName = "Yield"; yieldName += index;
	if(!(RatioFit && index<=4)) Yield[index]  = new RooRealVar(yieldName,yieldName,3000,0,50000);
      } else {
	TString fyieldName = "fYield"; fyieldName += index;
	fYield[index] = new RooRealVar(fyieldName,fyieldName,1.0);
	if(IsoConst){
	  if(index>=5 && index<=6 || index==14) IndfYield[index] = index-5+2*(index%2);
	  else if(index>=7 && index<=8 || index==16) IndfYield[index] = index-5+2*(index%2)-2;
	  else if(index== 3 || index== 4 || index==11 || index==12 || index==15 || 
		  index==23 || index==24 || floatDinDss&&index>=27&&index<=28) IndfYield[index] = index-2;
	  else if(index>=17 && index<=18) IndfYield[index] = index+4;
	  else if(index>=19 && index<=20) IndfYield[index] = index+2;
	  else if(index>=25 && index<=26) IndfYield[index] = index-16+(index%2)*4;
	  else if(index>=27 && index<=28) IndfYield[index] = index-18+(index%2)*4;
	  else if(index>=29 && index<=30){ 
	    if(floatDinDss) IndfYield[index] = index-4;
	    else IndfYield[index] = index-20-(index%2);
	  } else if(index>=31 && index<=32){ 
	    if(floatDinDss) IndfYield[index] = index-6;
	    else IndfYield[index] = index-22-(index%2);
	  }
	  else if(index>=37 && index<=38) IndfYield[index] = index-4;
	  else if(index>=39 && index<=40) IndfYield[index] = index-6;
	} else {
	  if(index>=5 && index<=8 || index==14 || index==16) IndfYield[index] = index-5+2*(index%2);
	  else if(index>=17 && index<=20) IndfYield[index] = index+4;
	  else if(index>=25 && index<=28) IndfYield[index] = index-16+(index%2)*4;
	  else if(index>=29 && index<=32){ 
	    if(floatDinDss) IndfYield[index] = index-4;
	    else IndfYield[index] = index-20-(index%2==0);
	  }
	  else if(index>=37 && index<=40) IndfYield[index] = index-4;
	} //IsoConst
      }
    }
  }
  for(int i=0; i<9; i++){
    for(int chan=0; chan<8; chan++){
      if(!doChannel(chan)) continue;
      if(IndPdf[chan>3][i]<0) continue;
      int index = IndPdf[chan>3][i]+chan-4*(chan>3);
      TString yieldName = "Yield"; yieldName += index;
      if(!isFloated(index)){
	TString formula = "Yield"; formula += IndfYield[index]; formula += "*";
	formula += "fYield"; formula += index;
	Yield[index]  = new RooFormulaVar(yieldName,yieldName,formula,
					  RooArgList(*Yield[IndfYield[index]],*fYield[index]));
      } else if(RatioFit && index<=4){
	TString ratioName = "RatioTau"; ratioName += index;
	RatioTau[index] = new RooRealVar(ratioName,ratioName,0.1,0.0,0.5);
	TString formula = "Yield"; formula += index+8; formula += "*";
	formula += ratioName;
	Yield[index]  = new RooFormulaVar(yieldName,yieldName,formula,
					  RooArgList(*Yield[index+8],*RatioTau[index]));
      }
    }
  }
  candType0 = new RooCategory("candType0","candType0");
  candTypec = new RooCategory("candTypec","candTypec");
  for(int chan=0; chan<8; chan++) {
    if(!doChannel(chan)) continue;
    if(chan<2 || chan>=4&&chan<6 || IsoConst) candType0 -> defineType(candLab[chan],chan+1);
    else candTypec -> defineType(candLab[chan],chan+1);
  }
}

// PDFs are initialized
void initializePdfs(int rep) {
  fstream ResWFile; ResWFile.open("../DonutUtils/ResolutionWidths.txt",fstream::in);
  double ResWidth[4];
  for(int chan=0; chan<4; chan++) ResWFile >> ResWidth[chan];
  for(int chan=0; chan<8; chan++){
    if(!doChannel(chan)) continue;
    for(int i=0; i<9; i++){
      if(IndPdf[chan>3][i]<0) continue;
      int index = IndPdf[chan>3][i]+chan-4*(chan>3);
      if(pdf[index]){
	pdf[index]->Delete(); RDpdf[index]->Delete();
      }
      TString hname = "AWG82/results/keys/root/Fit/pdfKeys_"; 
      if(typeToy==5) hname = "AWG82/results/keys/root/Times/hTimes_"; 
      hname += index; hname += "_Fit.root";
      TFile hfile(hname); 
      TString h2Name = "h2"; 
      if(typeToy==5) {
	h2Name = "hTimes"; h2Name += nhTimes[index];
      }
      TH2F *h2 = (TH2F *)(hfile.Get(h2Name));
      if(h2==0){
	if(typeToy==5){
	  h2 = (TH2F *)(hfile.Get("hTimes0"));
	  nhTimes[index] = 1;
	} else {
	  cout<<hname<<" does not contain "<<h2Name<<endl;
	  return;
	}
      } else nhTimes[index]++;
      if(TypeFit.Contains("Conv") && index>=21&&index<=24) ConvKeys(h2, ResWidth[index-21]);
      h2->SetDirectory(0);
      TString dhName = "Rdh"; dhName += index;
      RDpdf[index] = new RooDataHist(dhName, dhName, RooArgList(*mmiss2,*pstarl), h2);
      TString pdfName = "pdf"; pdfName += index;
      pdf[index] = new RooHistPdf(pdfName,pdfName,RooArgList(*mmiss2,*pstarl),*RDpdf[index],2);
      h2->Delete();
    }
  }

  // PDF sums
  RooArgList pdfs[8],pars[8];
  for(int chan=0; chan<8; chan++){
    if(!doChannel(chan)) continue;
    for(int i=0; i<9; i++){
      if(IndPdf[chan>3][i]<0) continue;
      int index = IndPdf[chan>3][i]+chan-4*(chan>3);
      if(isFloated(index)) {
	if(RatioFit && index<=4) (static_cast <RooRealVar*>(RatioTau[index]))->setConstant(FixYield[index]);
	else (static_cast <RooRealVar*>(Yield[index]))->setConstant(FixYield[index]);
      }
      pdfs[chan].add(*pdf[index]);
      pars[chan].add(*Yield[index]);
    }
  }

  for(int chan=0; chan<8; chan++) {
    if(!doChannel(chan)) continue;
    TString pdfName = "pdfsum"; pdfName += chan;
    if(pdfsum[chan]) pdfsum[chan]->Delete();
    pdfsum[chan] = new RooAddPdf(pdfName,pdfName,pdfs[chan],pars[chan]);
  }

  if(sim0pdf) sim0pdf->Delete(); if(simcpdf) simcpdf->Delete();
  sim0pdf = new RooSimultaneous("sim0pdf","sim0pdf",*candType0);
  simcpdf = new RooSimultaneous("simcpdf","simcpdf",*candTypec);
  for(int chan=0; chan<8; chan++){
    if(!doChannel(chan)) continue;
    if(chan<2 || chan>=4&&chan<6 || IsoConst) sim0pdf -> addPdf(*pdfsum[chan],candLab[chan]);
    else simcpdf -> addPdf(*pdfsum[chan],candLab[chan]);
  }
}

Bool_t parseToBool(RooStringVar &rsv, Bool_t &val){
  TString ValString = rsv.getVal();
  if(ValString == "1") val = true;
  else val = false;
  return kTRUE;
}

void readDataset(){

  // Initialized the true yields and the tree with its corresponding weights
  for(int i=0; i<70; i++) truYield[i] = 0;
  Int_t MCType2,candType2;
  float candM2,candPstarLep, weight=-1;
  treeData->SetBranchAddress("MCType",&MCType2);
  treeData->SetBranchAddress("candType",&candType2);
  treeData->SetBranchAddress("candM2",&candM2);
  treeData->SetBranchAddress("candPstarLep",&candPstarLep);
  treeData->SetBranchAddress("weight",&weight);

  totWeight = new RooRealVar("totWeight","totWeight",0.,100.);
  RooArgSet myVars0(*mmiss2,*pstarl,*candType0,*totWeight);
  RooArgSet myVarsc(*mmiss2,*pstarl,*candTypec,*totWeight);
  if(data0) data0->Delete();if(datac) datac->Delete();
  data0 = new RooDataSet("data0","data0",myVars0);
  datac = new RooDataSet("datac","datac",myVarsc);
  double entries = treeData->GetEntries();
  TRandom3 rand(0); 
  cout<<"DonutFitFinal::readDataset - Looping over events in "<<rootfile<<" and assigning weights "<<endl;
  for (int evt = 0 ; evt < entries; evt ++) {
    if(typeToy==2) {
      treeData->GetEvent((int)rand.Uniform(entries));
    } else treeData->GetEvent(evt);
    if(!doChannel(candType2-1)) continue;
    if(candM2<m2min||candM2>m2max || candPstarLep<plmin||candPstarLep>plmax || weight<0||weight>100) continue;
    mmiss2->setVal(candM2);
    pstarl->setVal(candPstarLep);
    totWeight->setVal(weight);
    if(candType2<3 || candType2>4&&candType2<7 || IsoConst){
      candType0->setIndex(candType2);
      data0->add(RooArgSet(*mmiss2,*pstarl,*candType0,*totWeight));
    } else {
      candTypec->setIndex(candType2);
      datac->add(RooArgSet(*mmiss2,*pstarl,*candTypec,*totWeight));
    }
    for(int i=0; i<70; i++) if(isSample(i,MCType2,candType2)) truYield[i] += weight;
    nData[candType2-1] += weight;
  }
  data0->setWeightVar(*totWeight);
  datac->setWeightVar(*totWeight);

  // Setting the initial (actual) values
  if(typeToy==1) {
    for(int i=0; i<70; i++) {
      truYield[i] = (int)truYield[i];     // You integers as Poisson means
      truYieldFix[i] = (int)truYield[i];  // truYieldFix has a fixed number of events, not the mean?
    }
  }
  if(typeToy>1) for(int i=0; i<70; i++) truYieldFix[i] = truYield[i];
  fstream ini, os;
  if(typeToy==-1) {
    ini.open(inifile,fstream::out);
    os.open(constfile,fstream::out);
  } else {
    ini.open(inifile,fstream::in);
    os.open(constfile,fstream::in);
  }
  TString buffer;
  for(int i=0; i<9; i++){
    for(int chan=0; chan<8; chan++){
      if(!doChannel(chan)) continue;
      if(IndPdf[chan>3][i]<0) continue;
      int index = IndPdf[chan>3][i]+chan-4*(chan>3);
      if(!isFloated(index)) {
	double fraction = truYield[index]/truYield[IndfYield[index]];
	if(typeToy==-1){
	  double n = truYield[index], d = truYield[IndfYield[index]];
	  double error = sqrt(n/d/d+n*n/d/d/d);
	  os<<index<<"\t"<<RoundNumber(fraction,4)<<"  +-  "<<RoundNumber(error,4)<<endl;
	} else if(typeToy!=1) os>>index>>fraction>>buffer>>buffer;
	fYield[index]->setVal(fraction);
      } else {
	double iniYield = truYield[index];
	if(typeToy==-1) ini<<index<<"\t"<<iniYield<<endl;
	else ini>>index>>iniYield;
	if(isData) truYield[index] = iniYield;
	if(rootfile.Contains("_4_Run") || rootfile.Contains("_14_Run")) iniYield /= 2.;
	if(RatioFit && index<=4) (static_cast <RooRealVar*>(RatioTau[index]))->setVal(truYield[index]/truYield[index+8]);
	else (static_cast <RooRealVar*>(Yield[index]))->setVal(iniYield);
      }
    }
  }
  if(isData){
    for(int i=0; i<9; i++){
      for(int chan=0; chan<8; chan++){
	if(!doChannel(chan)) continue;
	if(IndPdf[chan>3][i]<0) continue;
	int index = IndPdf[chan>3][i]+chan-4*(chan>3);
	if(!isFloated(index)) {
	  truYield[index] = fYield[index]->getVal()*truYield[IndfYield[index]];
	}
      }
    }
  }
}

void generateToy(){
  TRandom3 rand(0); 
  RooArgSet myVars(*mmiss2,*pstarl);
  RooArgSet myVars0(*mmiss2,*pstarl,*candType0);
  RooArgSet myVarsc(*mmiss2,*pstarl,*candTypec);
  if(data0) data0->Delete();if(datac) datac->Delete();
  data0 = new RooDataSet("data0","data0",myVars0);
  datac = new RooDataSet("datac","datac",myVarsc);
  for(int chan=0; chan<8; chan++){
    if(!doChannel(chan)) continue;
    for(int i=0; i<9; i++){
      if(IndPdf[chan>3][i]<0) continue;
      int index = IndPdf[chan>3][i]+chan-4*(chan>3);
      truYield[index] = rand.Poisson(truYieldFix[index]);
      int nEvents = (int)truYield[index];
      if(nEvents==0) continue;
      if(protoDS[index]) protoDS[index]->Delete();
      TString proName = "proto"; proName += index;
      if(chan<2 || chan>3&&chan<6 || IsoConst) {
	protoDS[index] = new RooDataSet(proName, proName, *candType0);
	candType0->setIndex(chan+1);
	for(int evt=0; evt<nEvents; evt++) protoDS[index]->add(*candType0);
      }else {
	protoDS[index] = new RooDataSet(proName, proName, *candTypec);
	candTypec->setIndex(chan+1);
	for(int evt=0; evt<nEvents; evt++) protoDS[index]->add(*candTypec);
      }
      if(chan<2 || chan>3&&chan<6) data0->append(*(pdf[index]->generate(myVars,*protoDS[index],false)));
      else  datac->append(*(pdf[index]->generate(myVars,*protoDS[index],false)));
    }
  }
}

Bool_t isFloated(int i){
  if(IsoConst){
    if(IsoConst==2){
      if(i>=1&&i<=2 || i>=9&&i<=13 || i==15 || i>=21&&i<=24 || floatDinDss==1&&i>=25&&i<=28 || i>=41&&i<=68) return true;
      else return false;
    } else {
      if(i>=1&&i<=2 || i>=9&&i<=10 || i==13 || i>=21&&i<=22 || floatDinDss==1&&i>=25&&i<=26 || i>=41&&i<=68) return true;
      else return false;
    }
  }
  if(i>=1&&i<=4 || i>=9&&i<=13 || i==15 || i>=21&&i<=24 || floatDinDss==1&&i>=25&&i<=28 || i>=33&&i<=36 || i>=41&&i<=68) return true;
  return false;
}
int doChannel(int chan){
  if(IsoConst) return 3;
  if(TypeFit.Contains("D0") && (chan==0 || chan==1 || chan==4 || chan==5)) return 1; 
  if(TypeFit.Contains("Dc") && (chan==2 || chan==3 || chan==6 || chan==7)) return 2; 
  return 0;
}

void readConfigFile(TString configFile){
  for(int i=0; i<70; i++) FixYield[i] = false;
  RooStringVar rsvrootfile   ("rootfile",   "rootfile", "");
  RooStringVar rsvlatexfile  ("latexfile",  "latexfile", "");
  RooStringVar rsvplotfile   ("plotfile",   "plotfile", "");
  RooStringVar rsvtextfile   ("textfile",   "textfile", "");
  RooStringVar rsvpullfile   ("pullfile",   "pullfile", "");
  RooStringVar rsvconstfile  ("constfile",  "constfile","");
  RooStringVar rsvinifile  ("inifile",  "inifile","");
  RooStringVar rsvTypeFit   ("TypeFit",   "TypeFit", "");
  RooStringVar rsvconstDssD0 ("ConstDssD0",     "ConstDssD0", "0");
  RooStringVar rsvconstDssDs0("ConstDssDs0",    "ConstDssDs0","0");
  RooStringVar rsvconstDs0D0 ("ConstDs0D0",     "ConstDs0D0", "0");
  RooStringVar rsvconstDs0Ds0("ConstDs0Ds0",    "ConstDs0Ds0","0");
  RooStringVar rsvconstD0D0  ("ConstD0D0",      "ConstD0D0",  "0");
  RooStringVar rsvconstTauD0 ("ConstTauD0",     "ConstTauD0", "0");
  RooStringVar rsvconstTauDs0("ConstTauDs0",    "ConstTauDs0","0");
  RooStringVar rsvconstCombD0("ConstCombD0",    "ConstCombD0","1");
  RooStringVar rsvconstXFD0  ("ConstXFD0",      "ConstXFD0",  "1");
  RooStringVar rsvconstCombDs0("ConstCombDs0",  "ConstCombDs0","1");
  RooStringVar rsvconstXFDs0 ("ConstXFDs0",     "ConstXFDs0", "1");
  RooStringVar rsvconstCombDssD0("ConstCombDssD0","ConstCombDssD0", "0");
  RooStringVar rsvconstXFDssD0  ("ConstXFDssD0",  "ConstXFDssD0", "1");
  RooStringVar rsvconstCombDssDs0("ConstCombDssDs0", "ConstCombDssDs0", "0");
  RooStringVar rsvconstXFDssDs0 ("ConstXFDssDs0", "ConstXFDssDs0", "1");
  RooStringVar rsvconstContD0("ConstContD0",    "ConstContD0","1");
  RooStringVar rsvconstContDs0("ConstContDs0",  "ConstContDs0","1");
  RooStringVar rsvconstContDssD0("ConstContDssD0","ConstContDssD0", "1");
  RooStringVar rsvconstContDssDs0("ConstContDssDs0", "ConstContDssDs0", "1");
  RooStringVar rsvconstDssDpipi ("ConstDssDpipi",     "ConstDssDpipi", "0");
  RooStringVar rsvconstDssDspipi("ConstDssDspipi",    "ConstDssDspipi","0");
  RooArgSet fitterPars(rsvconstDssD0,rsvconstDssDs0,rsvconstDs0D0,rsvconstDs0Ds0,rsvconstD0D0,rsvconstTauD0,rsvconstTauDs0,
		       rsvconstCombD0,rsvconstXFD0);
  fitterPars.add(rsvconstCombDs0);    
  fitterPars.add(rsvconstXFDs0 );     
  fitterPars.add(rsvconstCombDssD0);  
  fitterPars.add(rsvconstXFDssD0  );  
  fitterPars.add(rsvconstCombDssDs0); 
  fitterPars.add(rsvconstXFDssDs0 );  
  fitterPars.add(rsvconstContD0);    
  fitterPars.add(rsvconstContDs0);    
  fitterPars.add(rsvconstContDssD0);  
  fitterPars.add(rsvconstContDssDs0); 
  fitterPars.add(rsvconstDssDpipi);  
  fitterPars.add(rsvconstDssDspipi); 
  fitterPars.add(rsvrootfile );       
  fitterPars.add(rsvlatexfile );      
  fitterPars.add(rsvplotfile );      
  fitterPars.add(rsvtextfile );      
  fitterPars.add(rsvpullfile );      
  fitterPars.add(rsvconstfile );      
  fitterPars.add(rsvinifile );      
  fitterPars.add(rsvTypeFit );      
  fitterPars.readFromFile(configFile);
  parseToBool(rsvconstDssD0,    FixYield[21]);  parseToBool(rsvconstDssD0,    FixYield[23]); 
  parseToBool(rsvconstDssDs0,   FixYield[22]);  parseToBool(rsvconstDssDs0,   FixYield[24]); 
  parseToBool(rsvconstDs0D0,    FixYield[13]);  parseToBool(rsvconstDs0D0,    FixYield[15]); 
  parseToBool(rsvconstDs0Ds0,   FixYield[10]);  parseToBool(rsvconstDs0Ds0,   FixYield[12]); 
  parseToBool(rsvconstD0D0,     FixYield[9]);   parseToBool(rsvconstD0D0,     FixYield[11]); 
  parseToBool(rsvconstTauD0,    FixYield[1]);  	parseToBool(rsvconstTauD0,    FixYield[3]);  
  parseToBool(rsvconstTauDs0,   FixYield[2]);  	parseToBool(rsvconstTauDs0,   FixYield[4]);  
  parseToBool(rsvconstCombD0    ,FixYield[51]); parseToBool(rsvconstCombD0    ,FixYield[53]);
  parseToBool(rsvconstXFD0      ,FixYield[41]); parseToBool(rsvconstXFD0      ,FixYield[43]);
  parseToBool(rsvconstCombDs0   ,FixYield[52]); parseToBool(rsvconstCombDs0   ,FixYield[54]);
  parseToBool(rsvconstXFDs0     ,FixYield[42]); parseToBool(rsvconstXFDs0     ,FixYield[44]);
  parseToBool(rsvconstCombDssD0 ,FixYield[55]); parseToBool(rsvconstCombDssD0 ,FixYield[57]);
  parseToBool(rsvconstXFDssD0   ,FixYield[45]); parseToBool(rsvconstXFDssD0   ,FixYield[47]);
  parseToBool(rsvconstCombDssDs0,FixYield[56]); parseToBool(rsvconstCombDssDs0,FixYield[58]);
  parseToBool(rsvconstXFDssDs0  ,FixYield[46]); parseToBool(rsvconstXFDssDs0  ,FixYield[48]);  
  parseToBool(rsvconstContD0    ,FixYield[61]); parseToBool(rsvconstContD0    ,FixYield[63]);
  parseToBool(rsvconstContDs0   ,FixYield[62]); parseToBool(rsvconstContDs0   ,FixYield[64]);
  parseToBool(rsvconstContDssD0 ,FixYield[65]); parseToBool(rsvconstContDssD0 ,FixYield[67]);
  parseToBool(rsvconstContDssDs0,FixYield[66]); parseToBool(rsvconstContDssDs0,FixYield[68]);
  parseToBool(rsvconstDssDpipi,  FixYield[33]); parseToBool(rsvconstDssDpipi,  FixYield[35]); 
  parseToBool(rsvconstDssDspipi, FixYield[34]); parseToBool(rsvconstDssDspipi, FixYield[36]); 
  rootfile = rsvrootfile.getVal();
  latexfile = rsvlatexfile.getVal();
  plotfile = rsvplotfile.getVal();
  textfile = rsvtextfile.getVal();
  pullfile = rsvpullfile.getVal();
  constfile = rsvconstfile.getVal();
  inifile = rsvinifile.getVal();
  TypeFit = rsvTypeFit.getVal();
 }

void printAsLatex(){

  double totFit=0, totTru=0;
  TString TorEshort = "true", TorElong = " true";
  if(isData){TorEshort = "exp"; TorElong = " expected";}
  fstream os;
  os.open(latexfile,fstream::out);
  RooNLLVar *nll[2] = {nll0, nllc};
  if(IsoConst) nll[1] = nll0;
  os << "\\documentclass[6pt]{article}\n\\usepackage{rotating}\\usepackage{lscape}\n\\usepackage{amsmath}\n\\usepackage{colortbl}"<<endl<<
    "\\usepackage[left=2cm,top=1cm,right=3cm,nohead,nofoot]{geometry}\n\\pagestyle{empty}\n\\input{babarsym.tex}"<<endl<<
    "\\newcommand{\\Red} [1]  {\\textcolor{red}{#1}}\n\\begin{document} "<<endl;
  for(int i=0; i<2; i++){
    if(IsoConst==0){
      if(!TypeFit.Contains("D0") && i==0) continue;
      if(!TypeFit.Contains("Dc") && i==1) continue;
    }
    double likeli = nll[i]->getVal();
    totFit=0; totTru=0;
    os << "\\begin{tabular}{l l l r @{ $\\pm$ } l l}\\hline\\hline" << endl <<
      "Channel & Component & $N_\\mathrm{"<<TorEshort<<"}$ & \\multicolumn{2}{c}{$N_\\mathrm{fit}$} & Pull\\\\ \\hline" << endl;
    printYieldsChannel(os, 2*i, totFit, totTru);
    printYieldsChannel(os, 2*i + 1, totFit, totTru);
    os << "\\hline" << endl << "\\end{tabular}" << endl << 
      "\\begin{tabular}{l l l r @{ $\\pm$ } l l}\\hline\\hline" << endl <<
      "Channel & Component & $N_\\mathrm{"<<TorEshort<<"}$ & \\multicolumn{2}{c}{$N_\\mathrm{fit}$} & Pull\\\\ \\hline" << endl;
    printYieldsChannel(os, 2*i + 4, totFit, totTru);
    printYieldsChannel(os, 2*i + 5, totFit, totTru);
    os <<  "\\hline" << endl << "\\end{tabular}"<< endl;
    os <<" \\\\ \\\\"<< endl<<endl<<"\\begin{center}Total events: "<<RoundNumber(totTru,0)<<TorElong<<" and "<<
      RoundNumber(totFit,0)<<" fitted  -  The total likelihood is "<<RoundNumber(likeli,1)<<" \\end{center} \\vspace{0.8in}"<<endl<<endl;
    cout<<"Likelihood "<<i<<" is "<<likeli<<endl;
  }
    os << endl << "\\newpage"<< endl<< endl;
  if(TypeFit.Contains("D0") && floatDinDss) printCorrelMat(os,0);
  if(TypeFit.Contains("Dc") && floatDinDss) printCorrelMat(os,1);
  os << endl<< "\\end{document}"<< endl;
  cout << endl << endl << "LaTeX Fit Results printed in " << latexfile << endl << endl;
}

void printYieldsChannel(ostream& os, int channel, double &dF, double &dT){
  double totFit=0, totTru=0;
  int digits = 0, isD0=0;
  if(channel<2 || channel>3&&channel<6) isD0 = 1;
  TString sep = " & ", Labels[10], chanLabels[9];
  chanLabels[1]  = "\\Dstarz";            chanLabels[0]  = "\\Dz";
  chanLabels[3]  = "\\Dstarp";            chanLabels[2]  = "\\Dp";
  chanLabels[5]  = "$\\Dstarz\\piz$";     chanLabels[4]  = "$\\Dz\\piz$";
  chanLabels[7]  = "$\\Dstarp\\piz$";     chanLabels[6]  = "$\\Dp\\piz$";
  Labels[0]  = "$\\Dstarz\\taum\\nutb$";  Labels[1]  = "$\\Dstarz\\ellm\\nulb$";
  Labels[2]  = "$\\Dz\\taum\\nutb$";      Labels[3]  = "$\\Dz\\ellm\\nulb$";
  Labels[4]  = "$D^{**}\\ellm\\nulb$";    Labels[5]  = "Charge XF";
  Labels[6]  = "Comb. BG";                Labels[7]  = "Continuum";
  Labels[8]  = "$\\DorDsz\\ellm\\nulb$";  Labels[9]  = "$D^{(*)}\\pi\\pi\\ellm\\nulb$";
  if(!isD0){
    Labels[0]  = "$\\Dstarp\\taum\\nutb$";  Labels[1]  = "$\\Dstarp\\ellm\\nulb$";
    Labels[2]  = "$\\Dp\\taum\\nutb$";      Labels[3]  = "$\\Dp\\ellm\\nulb$";
    Labels[8]  = "$\\DorDsp\\ellm\\nulb$";
  }
  int IndLab[2][9] = {{2,3,0,1,4,5,6,7,9},{4,1,3,5,6,7,9,-1,-1}};
  for(int i=0; i<9; i++){
    if(channel>3 && IndPdf[1][i]<0) continue;
    int index = IndPdf[0][i]+channel;
    if(channel>3) index = IndPdf[1][i]+channel-4;
    totFit += Yield[index]->getVal();
    totTru += truYield[index];
  }
  os << chanLabels[channel];
  for(int i=0; i<9; i++){
    if(i==7){
      RooDataSet *rds;
      if(isD0){
	TString candReduce = "candType0==candType0::"; candReduce += candLab[channel];
	rds = (RooDataSet*)data0->reduce(candReduce);
      }else{
	TString candReduce = "candTypec==candTypec::"; candReduce += candLab[channel];
	rds = (RooDataSet*)datac->reduce(candReduce);
      }
      RooNLLVar nllRed("nllRed","nllRed",*pdfsum[channel],*rds,kTRUE);
      os <<"$"<< RoundNumber(nllRed.getVal(),1)<<"_{LH}$ ";
      delete rds;
    }
    if(!(IndLab[channel>3][i]<0)) {
      int index = IndPdf[0][i]+channel, indexLab = i;
      if(channel==0 || channel==2) indexLab = IndLab[0][i];
      if(channel>3) {
	index = IndPdf[1][i]+channel-4;
	indexLab = IndLab[1][i];
      }
      if(i==2) os<<"Total ";
      if(i==3){os<<"$"<<RoundNumber(totTru,digits); if(isData) os<<"_{exp}$ "; else os<<"_{true}$ ";}
      if(i==4) os<<"$"<<RoundNumber(totFit,digits)<<"_{fit}$ ";
      if(i==5 && isData) os<<"$"<<RoundNumber(nData[channel],digits)<<"_{data}$ ";
      os << sep << Labels[indexLab] << sep << RoundNumber(truYield[index],digits) << sep << RoundNumber(Yield[index]->getVal(),digits) << sep;
      double yError=-99;
      if(isFloated(index)) {  
	if(RatioFit && index<=4) {
	  double Dval=Yield[index+8]->getVal(), Derror=(static_cast <RooRealVar*>(Yield[index+8]))->getError();
	  double Rval=RatioTau[index]->getVal(), Rerror=(static_cast <RooRealVar*>(RatioTau[index]))->getError();
	  yError = sqrt(pow(Dval*Rerror,2)+pow(Rval*Derror,2));
	}else yError = (static_cast <RooRealVar*>(Yield[index]))->getError();
	os<<RoundNumber(yError,digits);
      }else os << "---";
      os << sep;
      if(isFloated(index)){
	if(Yield[index]->isConstant()) os << "fixed";
	else os << RoundNumber(Yield[index]->getVal()-truYield[index],2,yError);
      }else os << "---";
    }
    os << " \\\\"<<endl ;
  }
  os << " \\hline" << endl;
  dF += totFit; dT += totTru;
  return;
}

void printCorrelMat(ostream& os, int isDc){
  TString sep = " & ",redContentPre = "\\Red{",redContentPost = "}", zOrp = "z";
  if(isDc)zOrp = "p";
  int IndVar[11] = {1,13,9,2,10,21,25,55,22,26,56};
  TString topLabels[11],sideLabels[11];
  RooFitResult *fitres[2] = {fitres0, fitresc};
  RooAbsReal *Vars[11];
  for(int i=0; i<11; i++) Vars[i] = Yield[IndVar[i]+2*isDc];
  if(RatioFit){ 
    Vars[0] = RatioTau[1+2*isDc];
    Vars[3] = RatioTau[2+2*isDc];
  }
  sideLabels[0]  = "$\\Dz\\tau\\nu\\Rightarrow\\Dz$";
  sideLabels[1]  = "$\\Dstarz\\ell\\nu\\Rightarrow\\Dz$";
  sideLabels[2]  = "$\\Dz\\ell\\nu\\Rightarrow\\Dz$";
  sideLabels[3]  = "$\\Dstarz\\tau\\nu\\Rightarrow\\Dstarz$";
  sideLabels[4]  = "$\\Dstarz\\ell\\nu\\Rightarrow\\Dstarz$";
  sideLabels[5]  = "$D^{**}(\\ell/\\tau)\\nu\\Rightarrow\\Dz\\piz$";
  sideLabels[6]  = "$\\Dstarz(\\ell/\\tau)\\nu\\Rightarrow\\Dz\\piz$";
  sideLabels[7]  = "Comb. $\\Dz\\piz$";
  sideLabels[8]  = "$D^{**}(\\ell/\\tau)\\nu\\Rightarrow\\Dstarz\\piz$";
  sideLabels[9]  = "$\\Dstarz(\\ell/\\tau)\\nu\\Rightarrow\\Dstarz\\piz$";
  sideLabels[10] = "Comb. $\\Dstarz\\piz$";
  topLabels[0]  = "$\\Dz\\tau$";      topLabels[1]  = "$\\Dstarz\\ell$";
  topLabels[2]  = "$\\Dz\\ell$";      topLabels[3]  = "$\\Dstarz\\tau$";
  topLabels[4]  = "$\\Dstarz\\ell$";  topLabels[5]  = "$D^{**}$";
  topLabels[6]  = "$D^{*}$";          topLabels[7]  = "Comb.";
  topLabels[8]  = "$D^{**}$";         topLabels[9]  = "$D^{*}$";
  topLabels[10] = "Comb.";
  if(isDc){
    sideLabels[0]  = "$\\Dp\\tau\\nu\\Rightarrow\\Dp$";
    sideLabels[1]  = "$\\Dstarp\\ell\\nu\\Rightarrow\\Dp$";
    sideLabels[2]  = "$\\Dp\\ell\\nu\\Rightarrow\\Dp$";
    sideLabels[3]  = "$\\Dstarp\\tau\\nu\\Rightarrow\\Dstarp$";
    sideLabels[4]  = "$\\Dstarp\\ell\\nu\\Rightarrow\\Dstarp$";
    sideLabels[5]  = "$D^{**}(\\ell/\\tau)\\nu\\Rightarrow\\Dp\\piz$";
    sideLabels[6]  = "$\\Dstarp(\\ell/\\tau)\\nu\\Rightarrow\\Dp\\piz$";
    sideLabels[7]  = "Comb. $\\Dp\\piz$";
    sideLabels[8]  = "$D^{**}(\\ell/\\tau)\\nu\\Rightarrow\\Dstarp\\piz$";
    sideLabels[9]  = "$\\Dstarp(\\ell/\\tau)\\nu\\Rightarrow\\Dstarp\\piz$";
    sideLabels[10] = "Comb. $\\Dstarp\\piz$";
    topLabels[0]  = "$\\Dp\\tau$";      topLabels[1]  = "$\\Dstarp\\ell$";
    topLabels[2]  = "$\\Dp\\ell$";      topLabels[3]  = "$\\Dstarp\\tau$";
    topLabels[4]  = "$\\Dstarp\\ell$";  
  }
  os<<"\\begin{center}\\begin{tabular}{| l | rrr|rr | rrr|rrr |}\\cline{2-12}"<<endl;
  os<<"\\multicolumn{1}{l|}{} & \\multicolumn{3}{c|}{\\D"<<zOrp<<"} & \\multicolumn{2}{c|}{\\Dstar"<<zOrp<<"} & "<<
    "\\multicolumn{3}{c|}{$\\Dstar"<<zOrp<<"\\piz$} & \\multicolumn{3}{c|}{$\\D"<<zOrp<<"\\piz$}\\\\ \\cline{2-12}"<<
    endl<<"\\multicolumn{1}{l|}{}";
  for(int i=0; i<11; i++) os<<sep<<topLabels[i];
  os<<"\\\\ \\hline"<<endl;
  for(int row=0; row<11; row++){
    os<<sideLabels[row];
    for(int col=0; col<11; col++){
      os<<sep;
      bool red = false; double cor = 0.0, den = 1;
      if (!Vars[row]->isConstant() && !Vars[col]->isConstant())
	cor = fitres[isDc]->correlation(*Vars[row],*Vars[col]);
       else den = 0;
      if (row != col && fabs(cor) > 0.2) red = true;
      if (red) os << redContentPre;
      os << RoundNumber(cor,2,den);
      if (red) os << redContentPost;
    }
    os<<" \\\\"<<endl;
    if(row==2 || row==4 || row==6) os<<"\\hline "; 
  }
  os <<"\\hline \\end{tabular} \\end{center} "<<endl;

}
void printAsText(){
  fstream os;
  os.open(textfile,fstream::out);
  int IndChan[8] = {0,4,1,5,2,6,3,7}, digits = 1;
  for(int ichan=0; ichan<8; ichan++){
    int chan = IndChan[ichan];
    if(!doChannel(chan)) continue;
    for(int i=0; i<9; i++){
      int index = IndPdf[0][i]+chan;
      if(chan>3) index = IndPdf[1][i]+chan-4;
      if(IndPdf[chan>3][i]<0) {
	os<<"----------"<<endl;
	break;
      }else{
	os<<"Yield["<<index<<"]\t"<<RoundNumber(truYield[index],digits)<<"\t"<<
	  RoundNumber(Yield[index]->getVal(),digits)<<"\t+/-\t";
	if(isFloated(index)) {
	  double errorY = (static_cast <RooRealVar*>(Yield[index]))->getError();
	  os<<RoundNumber(errorY,digits)<<"\t";
	  if(errorY>0) os<<RoundNumber((Yield[index]->getVal()-truYield[index])/errorY,2);
	  else os<<"fixed";
	  os<<endl;
	} else os<<"-\t-"<<endl;
      }
    }
    os<<endl;
  }
  os.close();
  cout << endl << "Text Fit Results printed in " << textfile << endl<<endl;
}


