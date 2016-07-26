#include "TH1F.h"
#include "TCut.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TString.h"
#include "mycode/cuts.hh"
#include <fstream>
#include <iostream>
using std::cout;
using std::endl;

void PlotdataR22R24(){
  TChain r22("ntp1");
  TChain r24("ntp1");
  r22.Add("AWG82/ntuples/small/Adddata30_Run1.root");
  r22.Add("AWG82/ntuples/small/Adddata30_Run2.root");
  r24.Add("AWG82/ntuples/small/Adddata_R24_Run1.root");
  r24.Add("AWG82/ntuples/small/Adddata_R24_Run2.root");
  r24.SetLineColor(2);

  TString var[4] = {"candM2>>h(40,-1.5,3)","candM2>>h(40,-0.6,1.5)",
		    "candPstarLep>>h(25,0,2.4)","mm2pi0>>h(35,-3,8)"};
  TCut cuts[4] = {"candQ2<5","candQ2<5"+Mva,"candQ2<5"+Mva,dssMvaAll};
  TString titles[4] = {"q^{2} < 5 GeV^{2}","q^{2} < 5 GeV^{2}, BDT cuts",
		       "q^{2} < 5 GeV^{2}, BDT cuts","D^{(*)}#pi^{0} sample cuts"};
  TH1F *h[2][4]; TLatex label; TLegend *leg[4];
  label.SetTextSize(0.07);
  TCanvas c("R22_R24","Comparison R22/R24",1000,800);
  c.Divide(2,2,0.001,0.001);
  gStyle->SetOptStat(0);
  for(int i=0; i<4; i++){
    c.cd(i+1);
    double n22 = r22.Draw(var[i],cuts[i]);
    TH1F *hTemp = (TH1F*)gDirectory->Get("h");
    h[0][i] = (TH1F*)hTemp->Clone("r22_"+i);
    double n24 = r24.Draw(var[i],cuts[i]);
    hTemp = (TH1F*)gDirectory->Get("h");
    h[1][i] = (TH1F*)hTemp->Clone("r24_"+i);
    h[0][i]->SetTitle("");
    double maxi = h[0][i]->GetMaximum();
    if(maxi<h[1][i]->GetMaximum()) maxi = h[1][i]->GetMaximum();
    h[0][i]->SetMaximum(maxi*1.1);
    h[0][i]->SetXTitle("m^{2}_{miss} [Gev^{2}]");
    if(i==2) h[0][i]->SetXTitle("p*_{l} [Gev]");
    h[0][i]->Draw();
    h[1][i]->Draw("same");
    if(i==2) leg[i] = new TLegend(0.1,.74,0.4,0.9);
    else leg[i] = new TLegend(0.6,.74,0.9,0.9);
    leg[i]->SetFillColor(0);
    TString tag = "R22 ("; tag += (int)n22; tag += ")";
    leg[i]->AddEntry(h[0][i],tag);
    tag = "R24 ("; tag += (int)(n24); tag += ")";
    leg[i]->AddEntry(h[1][i],tag);
    leg[i]->Draw();
    label.SetNDC(kTRUE);
    label.DrawLatex(0.13,0.92,titles[i]);  
  }
  c.SaveAs("babar_code/R24/dataR22_R24.eps");
}

