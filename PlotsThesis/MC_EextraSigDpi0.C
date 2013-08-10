#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TChain.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "babar_code/Styles/Styles.cc"
#include "../DonutUtils/cuts.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void MC_EextraSigDpi0(){

  Styles style2; style2.setPadsStyle(1); style2.applyStyle();

  TChain tree("ntp1");
  tree.Add("AWG82/ntuples/small/RAll_RunAll.root");

  TCanvas can("can","Mmiss resolution comparison");
  //can.Divide(2,1);
  TPad *cPad = (TPad *)can.cd(0);
  int colors[] = {4,2,4,28};
  TH1F* h[4]; double n[4];
  for(int i=0; i<4; i++) {
    TString hname = "h"; hname += i;
    h[i] = new TH1F(hname,"",26, 0, 1.3);
    h[i]->SetLineWidth(2);
    h[i]->SetLineColor(colors[i]);
  }
  MvaAll *= "weight"; dssMvaAll *= "weight";
  n[0] = tree.Draw("candEExtra>>h0",MvaAll);
  n[1] = tree.Draw("candEExtra>>h1",dssMvaAll);
  h[0]->Scale(1000/h[0]->Integral());
  h[1]->Scale(1000/h[1]->Integral());
  //h[0]->SetMaximum(h[0]->GetMaximum()*1.2);
  h[0]->Draw();
  style2.fixYAxis(h[0],cPad);
  style2.setTitles(h[0],"E_{Extra} (GeV)","Entries/(50 MeV)");
  h[1]->Draw("same");
  double legW = 0.2, legH = 0.17;
  double legX = 1-style2.PadRightMargin-0.05, legY = 1-style2.PadTopMargin-0.06;
  TLegend leg(legX-legW, legY-legH, legX, legY);
  leg.SetTextSize(style2.LabelSize); leg.SetFillColor(0); leg.SetTextFont(style2.nFont);
  leg.SetBorderSize(0);
  TString label = "Signal sample";
  leg.AddEntry(h[0],label);
  label = "D#pi^{0} sample"; 
  leg.AddEntry(h[1],label);
  leg.Draw();

//   cPad = (TPad *)can.cd(2);
//   n[2] = tree.Draw("candEExtra>>h2","(candPMiss>0.2&&candQ2>4&&(MCType==2||MCType==4)&&candType==2||(MCType==8||MCType==10)&&candType==4)*weight");
//   n[3] = tree.Draw("candEExtra>>h3","(candPMiss>0.2&&candQ2>4&&MCType>12)*weight");
//   h[2]->Scale(1000/h[2]->Integral());
//   h[3]->Scale(1000/h[3]->Integral());
//   h[2]->SetMaximum(h[2]->GetMaximum()*1.2);
//   h[2]->Draw();
//   style2.fixYAxis(h[2],cPad);
//   style2.setTitles(h[2],"E_{Extra} (GeV)","Entries/(50 MeV)","b)");
//   h[3]->Draw("same");
//   legW = 0.22;
//   TLegend leg2(legX-legW, legY-legH, legX, legY);
//   leg2.SetTextSize(style2.LabelSize); leg2.SetFillColor(0); leg2.SetTextFont(style2.nFont);
//   leg2.SetBorderSize(0);
//   label = "D* l #nu"; 
//   leg2.AddEntry(h[2],label);
//   label = "D** l #nu"; 
//   leg2.AddEntry(h[3],label);
//   leg2.Draw();

  TString pName = "public_html/MC_EextraSigDpi0.eps"; 
  can.SaveAs(pName);
  for(int i=0; i<4; i++) h[i]->Delete();

}

