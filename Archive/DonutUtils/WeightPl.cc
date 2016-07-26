//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: WeightPl.cc,v 1.4 2012/03/04 00:33:04 manuelf Exp $
//
// Description:
//      WeightPl - Fits a KEYS to the p*l spectrum of Offpeak and
//      continuum MC, and saves the ratio in a text file
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      10/03/03 manuelf -- Created
//------------------------------------------------------------------------

#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCut.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "RooFitCore/RooArgList.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooFormulaVar.hh"
#include "RooFitCore/RooAddPdf.hh"
#include "RooFitCore/RooAbsReal.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooDataHist.hh"
#include "RooFitCore/RooHistPdf.hh"
#include "DonutUtils/RooNDKeysPdfDonut.hh"
#include "DonutUtils/KeysUtils.cc"
#include "DonutUtils/weightManager.cc"
#include "DonutUtils/cuts.cc"
#include "DonutUtils/Styles.cc"
#include <fstream>
#include <iostream>
#include <ctime>

using namespace RooFit;
using namespace std;
using std::cout;
using std::endl;

int main(int argc, char *argv[]){
  time_t start,end;
  time (&start);
  double dif;
  if (argc < 3 || argc > 5 ) {
    cout << "USAGE: WeightPl smoothing Channel [factor=10] [extraCuts=candQ2>4]" << endl;
    return 0;
  }
  // Setting the input parameters
  TString smoo_s = argv[1]; double smoo = smoo_s.Atof()/100.;
  TString Channel_s = argv[2]; //int Channel = Channel_s.Atoi();
  double factor = 10;
  if(argc>3) {TString temp_s = argv[3]; factor = temp_s.Atof();} 
  TString cuts_s = "candQ2>4&&candMES>5.27";
  if (argc>4) cuts_s = argv[4];
  cuts_s.ReplaceAll("XX","&&");

  Styles style3; 
  style3.setPadsStyle(3); style3.applyStyle();

  //TString titles[] = {"D^{0}","D*^{0}","D^{+}","D*^{+}"};
  TString OffFile = "AWG82/ntuples/small/OffPeak_RunAll.root";
  TString udsFile = "AWG82/ntuples/small/uds_RunAll.root";
  TString ccbarFile = "AWG82/ntuples/small/ccbar_RunAll.root";
  double OffEntries = 0, contiEntries = 0, wConti=0;
  double wuds = 0.042759, wccbar = 0.051672;
  cuts_s += "&&candType=="; cuts_s += Channel_s;
  TCut cuts = basic+M2P; cuts += cuts_s;
  TChain OffChain("ntp1"), udsChain("ntp1"), ccbarChain("ntp1");
  OffChain.Add(OffFile);
  udsChain.Add(udsFile);
  ccbarChain.Add(ccbarFile);
  TTree *Off = OffChain.CopyTree(cuts);
  TTree *uds = udsChain.CopyTree(cuts);
  TTree *ccbar = ccbarChain.CopyTree(cuts);
  TTree *conti = new TTree("ntp1","ntp1");
  Float_t candPstarLep, weight=-1;
  uds->SetBranchAddress("candPstarLep",&candPstarLep);
  ccbar->SetBranchAddress("candPstarLep",&candPstarLep);
  conti->Branch("candPstarLep",&candPstarLep,"candPstarLep/F");
  conti->Branch("weight",&weight,"weight/F");
  weight = wuds;
  for (int evt = 0 ; evt < uds->GetEntries()/factor; evt ++) {
    uds->GetEvent(evt);
    conti->Fill();
    wConti += weight;
  }
  weight = wccbar;
  for (int evt = 0 ; evt < ccbar->GetEntries()/factor; evt ++) {
    ccbar->GetEvent(evt);
    conti->Fill();
    wConti += weight;
  }
  contiEntries = uds->GetEntries()/factor+ccbar->GetEntries()/factor;
  OffEntries = Off->GetEntries();

  // Setting the fitting parameters
  RooRealVar pstarl("candPstarLep","candPstarLep",0,2.4);
  RooRealVar RRVweight("weight","weight",0.,1.);
  RooDataSet OffData("OffData","OffData",Off,RooArgSet(pstarl));
  RooDataSet contiData("contiData","contiData",conti,RooArgSet(pstarl,RRVweight),"","weight");
  int nbins = 1000;
  TH1F hOff("hOff","KEYS fit to OffPeak",nbins,0,2.4);
  TH1F hconti("hconti","KEYS fit to cc+uds",nbins,0,2.4);
  TH1F hrat("hrat","",nbins,0,2.4);

  // Fitting
  time (&end);dif = difftime (end,start);
  cout<<dif<<" sec after loading. Fitting "<<OffEntries<<" of OffPeak and "<<contiEntries
      <<" of continuum with smoothing "<<smoo<<endl;
  time (&start);
  RooNDKeysPdfDonut OffKEYS("OffKEYS","OffKEYS",RooArgList(pstarl),OffData,"a",smoo,3.);
  time (&end);dif = difftime (end,start);
  cout<<dif<<" sec after OffPeak KEYing. Doing continuum "<<endl;
  time (&start);
  RooNDKeysPdfDonut contiKEYS("contiKEYS","contiKEYS",RooArgList(pstarl),contiData,"a",smoo*0.7,3.);
  time (&end);dif = difftime (end,start);
  cout<<dif<<" sec after continuum MC KEYing. Doing histograms "<<endl;
  time (&start);
  OffKEYS.fillHistogram(&hOff,RooArgList(pstarl));
  contiKEYS.fillHistogram(&hconti,RooArgList(pstarl));
  time (&end);dif = difftime (end,start);
  cout<<dif<<" sec after histograms "<<endl;
  time (&start);

  // Saving
  TString ratioName = "keys/weights/ContiWeightRatio_"; ratioName+=Channel_s; ratioName+="_t";
  ratioName+=smoo_s; ratioName+=".root";
  TFile ratioFile(ratioName,"RECREATE");
  if(factor>1)  hrat.Divide(&hOff,&hconti,hconti.Integral(),hOff.Integral());
  else {
    cout<<"Including normalization in the re-weighting. Off: "<<OffEntries<<", Conti: "<<wConti<<endl;
    hrat.Divide(&hOff,&hconti,hconti.Integral()*OffEntries,hOff.Integral()*wConti);
  }
  for(int i=1; i<nbins+1; i++) if(hrat.GetBinContent(i)>4)hrat.SetBinContent(i,4);
  hrat.Write(); hOff.Write(); hconti.Write(); ratioFile.Close();
  TString textName = "keys/weights/ContiWeightText_"; textName+=Channel_s; textName+="_";
  textName+=smoo_s; textName+=".txt";
  fstream textFile;
  textFile.open(textName,fstream::out);
  double width = 2.4/(double)nbins;
  double plBin = width/2.;
  for(int i=1; i<nbins+1; i++){
    textFile<<plBin<<"\t"<<hrat.GetBinContent(i)<<endl;
    plBin += width;
  }
  textFile.close();
  cout<<"Saved "<<textName<<endl;

  // Plotting
  TCanvas can("can","KEYS fits");
  can.Divide(3); TPad *cPad = (TPad *)can.cd(1);
  TH1F *dOff = new TH1F("dOff","",40,0,2.4);
  Off->Draw("candPstarLep>>dOff");
  dOff->Sumw2(); 
  style3.setMarkers(dOff, 0.6, 20);
  dOff->SetMaximum(dOff->GetMaximum()*1.16);
  dOff->Draw();
  style3.fixYAxis(dOff,cPad);
  style3.setTitles(dOff,"p*_{l} (GeV)","Events/(60 MeV)","a)");
  hOff.Scale(dOff->Integral()/hOff.Integral()*1000./40.);
  hOff.SetLineColor(4);
  hOff.Draw("c same");

  double legW = 0.5, legH = 0.17;
  double legX = 1-style3.PadRightMargin-0.01, legY = 1-style3.PadTopMargin-0.01;
  TLegend legOff(legX-legW, legY-legH, legX, legY);
  legOff.SetTextSize(style3.TitleSize); legOff.SetFillColor(0); legOff.SetTextFont(style3.nFont);
  legOff.SetBorderSize(0);
  TString label = "Data ("; label += (int)OffEntries; label += ")";
  legOff.AddEntry(dOff,label);
  legOff.AddEntry(&hOff,"KEYS");
  legOff.Draw();

  cPad = (TPad *)can.cd(2);
  TH1F *dconti = new TH1F("dconti","",40,0,2.4);
  conti->Draw("candPstarLep>>dconti","weight");
  dconti->Scale(2/(wuds+wccbar)); dconti->Sumw2();  
  style3.setMarkers(dconti, 0.6, 20);
  style3.fixYAxis(dconti,cPad);
  dconti->SetMaximum(dconti->GetMaximum()*1.16);
  dconti->Draw(); dconti->SetMinimum(0); 
  style3.setTitles(dconti,"p*_{l} (GeV)","Events/(60 MeV)","b)");
  hconti.Scale(dconti->Integral()/hconti.Integral()*1000./40.);
  hconti.SetLineColor(4);
  hconti.Draw("c same");
  TLegend legconti(legX-legW, legY-legH, legX, legY);
  legconti.SetTextSize(style3.TitleSize); legconti.SetFillColor(0); 
  legconti.SetTextFont(style3.nFont);  legconti.SetBorderSize(0);
  label = "MC ("; label += (int)contiEntries; label += ")";
  legconti.AddEntry(dconti,label);
  legconti.AddEntry(&hconti,"KEYS");
  legconti.Draw();

  cPad = (TPad *)can.cd(3);
  cPad->SetLeftMargin(0.11);
  TLine lin; lin.SetLineColor(4);
  hrat.SetMaximum(1.5);
  hrat.Draw();
  style3.setTitles(&hrat,"p*_{l} (GeV)","","","c)");
  lin.DrawLine(0,1,2.4,1);hrat.Draw("same");
  TString plotName = "keys/weights/ContiWeightPl_"; plotName+=Channel_s; plotName+="_";
  plotName+=smoo_s; plotName+=".eps";
  can.SaveAs(plotName);

  delete conti; 
  return 1;
}


