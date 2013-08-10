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
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

TString RoundNumber(double n, int e, double d=1);

void EvtSel_MmissRes(){

  Styles style2; style2.setPadsStyle(2); style2.applyStyle();

  TChain tree("ntp1");
  tree.Add("AWG82/ntuples/small/RAll_RunAll.root");

  TCanvas can("can","Mmiss resolution comparison");
  can.Divide(2,1);
  TPad *cPad = (TPad *)can.cd(1);
  TH1F* h[4]; double n[4];
  for(int i=0; i<4; i++) {
    TString hname = "h"; hname += i;
    h[i] = new TH1F(hname,"",40,-0.5,0.5);
    h[i]->SetLineWidth(2);
    if(i%2==0) h[i]->SetLineColor(4);
    else h[i]->SetLineColor(2);
  }
  TString cuts = "(MCType==1||MCType==3||MCType==7||MCType==9)";
  n[0] = tree.Draw("candM2-candM2Tru>>h0",cuts);
  n[1] = tree.Draw("candM2NF-candM2Tru>>h1",cuts);
  h[0]->Draw();
  style2.fixYAxis(h[0],cPad);
  style2.setTitles(h[0],"m^{2}_{miss}-m^{2}_{miss,true} (GeV^{2})","Events/(0.025 GeV^{2})","a)");
  //h[1]->Scale(h[0]->Integral()/h[1]->Integral());
  h[1]->Draw("same"); h[0]->Draw("same");
  double legW = 0.26, legH = 0.17;
  double legX = 1-style2.PadRightMargin-0.01, legY = 1-style2.PadTopMargin-0.01;
  TLegend leg(legX-legW, legY-legH, legX, legY);
  leg.SetTextSize(style2.LabelSize); leg.SetFillColor(0); leg.SetTextFont(style2.nFont);
  leg.SetBorderSize(0);
  TString label = "Re-fitted";
  leg.AddEntry(h[0],label);
  label = "No fit"; 
  leg.AddEntry(h[1],label);
  leg.Draw();

  can.cd(2);
  cuts = "MCType==5||MCType==11";
  n[2] = tree.Draw("candM2-candM2Tru>>h2",cuts);
  n[3] = tree.Draw("candM2NF-candM2Tru>>h3",cuts);
  h[2]->Draw();
  style2.setTitles(h[2],"m^{2}_{miss}-m^{2}_{miss,true} (GeV^{2})","Events/(0.025 GeV^{2})","b)");
  //h[3]->Scale(h[2]->Integral()/h[3]->Integral());
  h[3]->Draw("same");h[2]->Draw("same");
  TLegend leg2(legX-legW, legY-legH, legX, legY);
  leg2.SetTextSize(style2.LabelSize); leg2.SetFillColor(0); leg2.SetTextFont(style2.nFont);
  leg2.SetBorderSize(0);
  label = "Re-fitted"; 
  leg2.AddEntry(h[2],label);
  label = "No fit"; 
  leg2.AddEntry(h[3],label);
  leg2.Draw();

  TString pName = "public_html/EvtSel_MmissRes.eps"; 
  can.SaveAs(pName);
  for(int i=0; i<4; i++) h[i]->Delete();

}

TString RoundNumber(double n, int e, double d){
  if(d==0) return " - ";
  double neg = 1; if(n*d<0) neg = -1;
  double b = (int)(neg*n/d*pow(10.,(double)e)+0.5);
  b /= pow(10.,(double)e)*neg;
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(!result.Contains(".") && e != 0) result += ".";
  
  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<e-afterdot.Length(); i++)
    result += "0";
  return result;
}
