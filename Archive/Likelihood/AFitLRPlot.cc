/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: AFit
 * Author: A Bevan (a.j.bevan@qmul.ac.uk)
 *
 * Copyright (C) 2007 QMUL
 *****************************************************************************/

// -- CLASS DESCRIPTION [PDF] --
//
// Implementation of a class to produce plots of the ratio of signal to
// background likelihoods obtained in a likelihood fit to a data sample.
//
// To use this class, one has to provide it with
//  i)   signal PDF (RooAbsPdf*)
//  ii)  background PDF (RooAbsPdf*)
//  iii) variables used in the likelihood (RooArgSet*)
//  iv)  the data fitted (RooDataSet*)
//  v)   the fitted signal yield
//  vi)  the fitted background yield
//
// The output of this class is an eps file and a ROOT file.  The eps file
// shows the likelihood ratio plot of data, with an overlay of the
// expectation of background (RED) and of the background+signal (GREEN)
// based on the specified fit yields.  The ROOT file contains a TCanvas
// called can2 that was used to make the afforementioned eps file.  In addition
// to this canvas, there are a number of other objects that can be used
// for debugging any problems that may be encountered.
//

#include "AFit/AFitLRPlot.hh"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
 
#ifdef __BABARROOFIT__
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooAbsPdf.hh"
#include "RooFitCore/RooPlot.hh"
#include "RooFitCore/RooTreeDataStore.hh"
#else
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooAbsPdf.h"
#include "RooPlot.h"
#include "RooTreeDataStore.h"
#endif                                           

ClassImp(AFitLRPlot)
;

AFitLRPlot::AFitLRPlot() :
  _sigPDF(0)
  , _bgPDF(0)
  , _sigLike(0)
  , _bgLike(0)
  , _sGen(0)
  , _bgGen(0)
  , _dataCut("1")
  , _nToGen(100000) 
{
  setOutputFile("test.root");
  // ctor
}

AFitLRPlot::~AFitLRPlot() {}
 
RooPlot * AFitLRPlot::makePlot(RooAbsPdf * theSigPDF, RooAbsPdf * theBgPDF, RooArgSet *genVars, RooDataSet * theRefData,
			       Double_t nSig, Double_t nBG) {
  _sigPDF = theSigPDF;
  _bgPDF = theBgPDF;

  _sigLike  =new RooFormulaVar("sigLike", "signal hypothesis",     "@0", * _sigPDF);
  _bgLike   =new RooFormulaVar("bgLike",  "background hypothesis", "@0", * _bgPDF);
  _SovBLike =new RooFormulaVar("SovBLike",  "s/(S+B) likelihood ratio", "@0/(@0+@1)", RooArgList(*_sigLike, *_bgLike) );

  cout << "\tGenerating signal" << endl;
  _sGen = _sigPDF->generate(*genVars, _nToGen);
  _sGen->addColumn( *_sigLike );
  _sGen->addColumn( *_bgLike );
  _sGen->addColumn( *_SovBLike );

  cout << "\tGenerating background" << endl;
  _bgGen = _bgPDF->generate(*genVars, _nToGen);
  _bgGen->addColumn( *_sigLike );
  _bgGen->addColumn( *_bgLike );
  _bgGen->addColumn( *_SovBLike );

  _theData = (RooDataSet*)theRefData->reduce(_dataCut);
  _theData->addColumn( *_sigLike );
  _theData->addColumn( *_bgLike );
  _theData->addColumn( *_SovBLike );
  Int_t nTot = _theData->numEntries();

  cout << ClassName() << ": Finished generating signal and background distributions" << endl;
  /*
  TTree * treeSig = (TTree*)(_sGen->tree().Clone("signalTree"));
  TTree * treeBG  = (TTree*)(_bgGen->tree().Clone("backgroundTree"));
  TTree * treeDat = (TTree*)(_theData->tree().Clone("dataTree"));
  */

  TTree * treeSig = (TTree*)(((RooTreeDataStore*)_sGen->store())->tree().Clone("data") );
  TTree * treeBG = (TTree*)(((RooTreeDataStore*)_bgGen->store())->tree().Clone("data") );
  TTree * treeDat = (TTree*)(((RooTreeDataStore*)_theData->store())->tree().Clone("data") );

  double normSiga = treeSig->GetMaximum("sigLike");
  double normSigb = treeSig->GetMaximum("bgLike");

  double normBGa = treeBG->GetMaximum("sigLike");
  double normBGb = treeBG->GetMaximum("bgLike");

  cout << ClassName() << ": " <<normSiga << " " << normSigb << endl;
  cout << ClassName() << ": " <<normBGa << " " << normBGb << endl;

  TCanvas can("can", "");
  can.Divide(2,2);
  can.cd(1);
  Double_t range_LH[2] = {0.0, 1.0};
  Int_t nBins = 100;
  Double_t normalisationScale = (Double_t)nTot/(Double_t)_nToGen;
  TH1F hSig("hSig", "SovB", nBins, range_LH[0], range_LH[1]);
  TH1F hBG("hBG", "SovB", nBins, range_LH[0], range_LH[1]);
  TH1F hDat("hDat", "SovB", nBins, range_LH[0], range_LH[1]);
  hSig.SetLineColor(kBlue);
  hBG.SetLineColor(kRed);
  hDat.SetLineColor(kBlack);
  treeDat->Draw("SovBLike>>hDat");
  treeBG->Draw("SovBLike>>hBG");
  treeSig->Draw("SovBLike>>hSig");
  hBG.Scale(normalisationScale);
  hSig.Scale(normalisationScale);
  hSig.Draw();
  hBG.Draw("same");
  hDat.Draw("same");

  can.cd(2);
  TH1F hDatBG("hDatBG", "BG likelihood for data", nBins, range_LH[0], range_LH[1]);
  TH1F hDatSig("hDatSig", "signal likelihood for data", nBins, range_LH[0], range_LH[1]);
  hDatBG.SetLineColor(kRed);
  hDatSig.SetLineColor(kBlue);
  treeDat->Draw("bgLike>>hDatBG");
  treeDat->Draw("sigLike>>hDatSig");
  hDatBG.Draw();
  hDatSig.Draw("same");

  can.cd(3);
  TH1F hSigBG("hSigBG", "BG likelihood for signal", nBins, range_LH[0], range_LH[1]);
  TH1F hSigSig("hSigSig", "signal likelihood for signal", nBins, range_LH[0], range_LH[1]);
  hSigBG.SetLineColor(kRed);
  hSigSig.SetLineColor(kBlue);
  treeSig->Draw("bgLike>>hSigBG");
  treeSig->Draw("sigLike>>hSigSig");
  hSigBG.Draw();
  hSigSig.Draw("same");

  can.cd(4);
  TH1F hBGBG("hBGBG", "BG likelihood for background", nBins, range_LH[0], range_LH[1]);
  TH1F hBGSig("hBGSig", "signal likelihood for background", nBins, range_LH[0], range_LH[1]);
  hBGBG.SetLineColor(kRed);
  hBGSig.SetLineColor(kBlue);
  treeBG->Draw("bgLike>>hBGBG");
  treeBG->Draw("sigLike>>hBGSig");
  hBGBG.Draw();
  hBGSig.Draw("same");

  /*
   * Having made the debugging plots above, make the final LR ratio plot
   */
  TCanvas can2("can2", "");
  TH1F hSig2("hSig2", "SovB", nBins, range_LH[0], range_LH[1]);
  TH1F hBG2("hBG2", "SovB", nBins, range_LH[0], range_LH[1]);
  TH1F hTOT2("hTOT2", "SovB", nBins, range_LH[0], range_LH[1]);
  TH1F hDat2("hDat2", "SovB", nBins, range_LH[0], range_LH[1]);
  hSig2.SetLineColor(kBlue);
  hBG2.SetLineColor(kRed);
  hBG2.SetFillColor(kRed);
  hTOT2.SetLineColor(kGreen);
  hTOT2.SetFillColor(kGreen);
  hDat2.SetLineColor(kBlack);
  treeDat->Draw("SovBLike>>hDat2");
  treeBG->Draw("SovBLike>>hBG2");
  treeSig->Draw("SovBLike>>hSig2");
  hBG2.Scale(nBG/_nToGen);
  hSig2.Scale(nSig/_nToGen);
  hTOT2.Add(&hBG2, 1.0);
  hTOT2.Add(&hSig2, 1.0);
  hTOT2.SetTitle("");
  hTOT2.SetStats(0);
  hTOT2.GetXaxis()->SetTitle("Likelihood Ratio");
  hTOT2.GetYaxis()->SetTitle("Number of Entries");
  hTOT2.Draw();
  hBG2.Draw("same");
  hDat2.Draw("esame");

  TString epsFName = _outputFile;
  epsFName.ReplaceAll(".root", ".eps");
  can2.Print(epsFName);
  hTOT2.SetTitle("SovB");

  /*
   * Output histograms and cancases to a ROOT file
   */
  cout << ClassName() << ": Writing output to the following ROOT file: " << _outputFile << endl;
  TFile fo(_outputFile, "RECREATE");
  fo.cd();

  can.Write();
  can2.Write();
  hSig.Write();
  hBG.Write();
  hDat.Write();
  hDatBG.Write();
  hDatSig.Write();
  hSigBG.Write();
  hSigSig.Write();
  hBGBG.Write();
  hBGSig.Write();
  hSig2.Write();
  hBG2.Write();
  hTOT2.Write();
  hDat2.Write();

  fo.Write();
  fo.Close();

  delete _sigLike;
  delete _bgLike;
  delete _SovBLike;
  return 0;
} 
