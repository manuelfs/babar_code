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

void Francon(){

  Styles style2; style2.setPadsStyle(2); style2.applyStyle();

  TChain tree("ntp1");
  tree.Add("AWG82/ntuples/small/FitRAll_RunAll.root");
  TChain data("ntp1");
  data.Add("AWG82/ntuples/small/FitData_RunAll.root");

  TCanvas can("can","Mmiss resolution comparison");
  can.Divide(2,1);
  TPad *cPad = (TPad *)can.cd(1);
  TH1F* h[4]; double n[4];
  for(int i=0; i<4; i++) {
    TString hname = "h"; hname += i;
    h[i] = new TH1F(hname,"",50,-2.5,7.5);
    h[i]->SetLineWidth(2);
    if(i%2==1) h[i]->SetLineColor(4);
    else h[i]->Sumw2();
  }
  TString cuts = "(MCType==1||MCType==3||MCType==7||MCType==9)";
  n[0] = data.Draw("candM2>>h0","candType>4&&candIsMu==0");
  n[1] = tree.Draw("candM2>>h1","weight*(candType>4&&candIsMu==0)");
  h[0]->Draw();
  style2.fixYAxis(h[0],cPad);
  style2.setTitles(h[0],"m^{2}_{miss} (GeV^{2})","Events/(0.2 GeV^{2})","e");
  h[1]->Scale(h[0]->Integral()/h[1]->Integral());
  h[1]->Draw("same");
  double legW = 0.26, legH = 0.17;
  double legX = 1-style2.PadRightMargin-0.01, legY = 1-style2.PadTopMargin-0.01;
  TLegend leg(legX-legW, legY-legH, legX, legY);
  leg.SetTextSize(style2.LabelSize); leg.SetFillColor(0); leg.SetTextFont(style2.nFont);
  leg.SetBorderSize(0);
  TString label = "Data";
  leg.AddEntry(h[0],label);
  label = "MC"; 
  leg.AddEntry(h[1],label);
  leg.Draw();

  can.cd(2);
  cuts = "MCType==5||MCType==11";
  n[2] = data.Draw("candM2>>h2","candType>4&&candIsMu==1");
  n[3] = tree.Draw("candM2>>h3","weight*(candType>4&&candIsMu==1)");
  h[2]->Draw();
  style2.setTitles(h[2],"m^{2} (GeV^{2})","Events/(0.2 GeV^{2})","#mu");
  h[3]->Scale(h[2]->Integral()/h[3]->Integral());
  h[3]->Draw("same");
  TLegend leg2(legX-legW, legY-legH, legX, legY);
  leg2.SetTextSize(style2.LabelSize); leg2.SetFillColor(0); leg2.SetTextFont(style2.nFont);
  leg2.SetBorderSize(0);
  label = "Data"; 
  leg2.AddEntry(h[2],label);
  label = "MC"; 
  leg2.AddEntry(h[3],label);
  leg2.Draw();

  TString pName = "public_html/Francon.eps"; 
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
