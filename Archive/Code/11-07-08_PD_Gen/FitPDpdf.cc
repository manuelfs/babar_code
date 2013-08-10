//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: FitPdfSampleKEYS.cc,v 1.6 2011/02/22 01:01:50 manuelf Exp $
//
// Description:
//      FitPDpdf - Fits the PD spectrum
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      11/06/29 manuelf -- Adapted from FitPdfSampleKEYS.cc
//------------------------------------------------------------------------

#include "TString.h"
#include "TFile.h"
#include "TCut.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooAddPdf.hh"
#include "RooFitCore/RooAbsReal.hh"
#include "RooFitCore/RooDataSet.hh"
#include "DonutUtils/RooNDKeysPdfDonut.hh"
#include "DonutUtils/KeysUtils.cc"
#include "DonutUtils/Styles.cc"
#include <fstream>
#include <iostream>
#include <ctime>

using namespace std;
using namespace RooFit;
using std::cout;
using std::endl;

void Stitch(TH1F *hResult, TH1F *hTail, double stitch, int nbins, double entries[]);
void SaveHisto(TH1F *h2, TString hname, int Sam=0);

int main(int argc, char *argv[]){

  time_t start,end;
  time (&start);
  double dif;

  if (argc < 3 || argc > 7 ) {
    cout << "USAGE: FitPDpdf sample smoo1 [smoo2=-1] [stitch=1.7] [nbins=300] [dataFolder=RAll]"<< endl;
    return 0;
  }
  TString sam = argv[1]; int Sam = sam.Atoi();
  TString smoo1_s = argv[2]; double smoo1 = smoo1_s.Atof()/100.;
  TString scaleP_s = "200";
  TString smoo2_s = "-1"; double smoo2 = -1;
  if(argc>3)  {smoo2_s = argv[3]; smoo2 = smoo2_s.Atof()/100.;}
  double stitch = 1.7;
  if(argc>4)  {TString sTemp = argv[4]; stitch = sTemp.Atof();}
  int nbins = 300;
  if(argc>5)  {TString sTemp = argv[5]; nbins = sTemp.Atoi();}
  TString dataFolder = "RAll";
  if(argc>6)  dataFolder = argv[6]; 

  double m2cut[2];
  TString m2CutName = "babar_code/Systematics/MmissCut.txt";
  fstream m2File; m2File.open(m2CutName,fstream::in);
  m2File>>m2cut[0]>>m2cut[1];
  Styles style; style.setPadsStyle(1); style.applyStyle();
  double pdmin = 0, pdmax = 2.;
  bool doStitch = false; if(smoo2>0) doStitch = true;
  if(doStitch){
    int binCut = (int)(stitch/pdmax*(double)nbins)+1;
    stitch = (binCut-1)*pdmax/nbins;
  }
  TString adaptive = "a";
  TString Base =  "_"; Base += sam; Base += "_"; Base += smoo1_s; 
  if(doStitch) {Base+="and"; Base+=smoo2_s;}

  double entriesCut[]={0,0,0,0};
  float candPstarD, candM2, weight=-1;
  TString inputfile = nameData(Sam,dataFolder);
  TChain *tree = new TChain("ntp1");
  tree->Add(inputfile);
  tree->SetBranchAddress("candM2",&candM2);
  tree->SetBranchAddress("candPstarD",&candPstarD);
  tree->SetBranchAddress("weight",&weight);

  double totMCB = 0, totuds = 0, totccbar = 0, totdata = 0, totOffdata = 0;
  getNumberB("AWG82/ntuples/small/RAll_RunAll.root", "All", totMCB, totdata, totuds, totccbar, totOffdata);
  //double wMC = totdata/totMCB; 
  TRandom3 rand; 
  double entries = tree->GetEntries();
  TH1F *h2 = new TH1F("h2","KEYS",nbins,pdmin,pdmax);
  TString hnameFit = "keys/root/PDfitFinal/pdfKeys_"; hnameFit += sam; hnameFit += "_Fit.root";
  TString hname = "keys/root/fitPdf/hKeys"; hname += Base; hname += ".root";
  TFile* hfile = new TFile(hname,"RECREATE"); 
  hfile->Close(); hfile->Delete();
  int Times = 1;
  for(int rep = 0; rep < Times; rep++){
    for(int iCut=0; iCut<4; iCut++) entriesCut[iCut] = 0;
    RooRealVar pstard("candPstarD","candPstarD",pdmin,pdmax);
    RooRealVar totWeight("totWeight","totWeight",0.,100.);
    RooArgSet myVars(pstard,totWeight);
    RooDataSet  data("data","data",myVars);
    if(rep%50==0) cout<<"Repetition "<<rep+1<<" of "<< Times <<" - Smoothings: "<<smoo1_s<<"/"<<smoo2_s<<endl; 
    for (int evt = 0 ; evt < entries; evt ++) {
      tree->GetEvent(evt);
      if(candPstarD<pdmin||candPstarD>pdmax || weight<0||weight>100 ||candM2<=m2cut[0] ||candM2>m2cut[1]) continue;
      pstard.setVal(candPstarD);
      totWeight.setVal(weight);
      data.add(RooArgSet(pstard,totWeight));
      if(candPstarD>=stitch) entriesCut[0] += weight;
      else entriesCut[1] += weight;
      entriesCut[3] += weight;
    }
    data.setWeightVar(totWeight);
    RooNDKeysPdfDonut KEYS("KEYS","KEYS",RooArgList(pstard),data,adaptive,smoo1,3.);
    KEYS.fillHistogram(h2,RooArgList(pstard));
    time (&end);dif = difftime (end,start);
    if(rep%50==0) cout<<dif<<" seconds after making histogram"<<endl;
    time (&start);
    if(doStitch) {
      if(rep%50==0) cout<<"Stitching smoothings "<<smoo2_s<<" at p*D "<<stitch<<endl;
      RooNDKeysPdfDonut KEYS3("KEYS3","KEYS",RooArgList(pstard),data,adaptive,smoo2,3.);
      TH1F *hsti = new TH1F("hsti3","Histo to be stitched",nbins,pdmin,pdmax);
      KEYS3.fillHistogram(hsti,RooArgList(pstard));
      Stitch(h2, hsti, stitch, nbins,entriesCut);
      hsti->Delete();
    }
    hfile = new TFile(hname,"UPDATE"); 
    h2->Write();
    hfile->Close(); hfile->Delete();
    if(Times==1) SaveHisto(h2, hnameFit,Sam);
    if(Times>1) h2->Delete();
  }
  cout<<hname<<" saved"<<endl; 
  if(Times>1) {return Times;}

  // Plotting
  int nbinsPlot = 30;
  TCanvas can("can","KEYS fits");
  TH1F *dOff = new TH1F("dOff","",nbinsPlot,pdmin, pdmax);
  TString cutTree = "(candM2>"; cutTree += m2cut[0]; cutTree += "&&candM2<="; cutTree += m2cut[1]; 
  cutTree += ")*weight";
  tree->Draw("candPstarD>>dOff",cutTree);
  dOff->Sumw2(); 
  style.setMarkers(dOff, 0.6, 20);
  dOff->SetMaximum(dOff->GetMaximum()*1.16);
  dOff->SetMinimum(0);
  dOff->Draw();
  style.setTitles(dOff,"p*_{D} (GeV)","Events/(5 MeV)");
  h2->Scale(dOff->Integral()/h2->Integral()*(double)nbins/(double)nbinsPlot);
  h2->SetLineColor(4);
  h2->Draw("c same");

  double legW = 0.2, legH = 0.17;
  double legX = style.PadLeftMargin+legW+0.05, legY = 1-style.PadTopMargin-0.01;
  TLegend legOff(legX-legW, legY-legH, legX, legY);
  legOff.SetTextSize(style.TitleSize); legOff.SetFillColor(0); legOff.SetTextFont(style.nFont);
  legOff.SetBorderSize(0);
  TString label = "MC";
  legOff.AddEntry(dOff,label);
  legOff.AddEntry(h2,"KEYS");
  legOff.Draw();

  TString plotName = "keys/eps/fitPdf/PDFit_"; plotName+=sam; plotName+=".eps";
  can.SaveAs(plotName);
  return 1; 

}

void Stitch(TH1F *hResult, TH1F *hTail, double stitch, int nbins, double entries[]){

  int binRange = nbins/20;
  double pdmin = 0, pdmax = 2.;
  int binCut = (int)(stitch/(pdmax-pdmin)*(double)nbins+pdmin)+1;
  if(binCut+binRange>nbins) binRange = nbins-binCut-1;
  double factor = hResult->Integral(1,binCut-1)/hTail->Integral(binCut,nbins)*entries[0]/(entries[3]-entries[0]);
  for(int j=binCut-binRange+1; j<nbins+1; j++){
    double val=0;
    if(j>=binCut+binRange) {
      val = hTail->GetBinContent(j)*factor;
    }else  {
      double weight = (j-binCut+binRange)/(double)(2*binRange);
      val = (1.-weight)*hResult->GetBinContent(j) + weight*hTail->GetBinContent(j)*factor;
    }
    hResult->SetBinContent(j,val);
  }

  return;
}

// For combinatoric and crossfeed the same fit is used several times, 
// so several copies are saved
void SaveHisto(TH1F *h2, TString hname, int Sam){
  TFile* hfileFit = new TFile(hname,"RECREATE"); 
  hfileFit->cd();
  h2->Write();
  hfileFit->Close(); hfileFit->Delete();
  if(Sam==700){ // I might use it for pdf 69 and 70
    hname.ReplaceAll("48","55");
    SaveHisto(h2, hname);
    hname.ReplaceAll("55","57");
    SaveHisto(h2, hname);
  }
}

