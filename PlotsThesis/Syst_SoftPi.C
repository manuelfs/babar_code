#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCut.h"
#include "babar_code/Styles/Styles.cc"
#include "../DonutUtils/cuts.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

TString RoundNumber(double n, int e, double d=1);

void Syst_SoftPi(){

  Styles style2; style2.setPadsStyle(2); style2.applyStyle();

  double legW = 0.38, legH = 0.17;
  double legX = 1-style2.PadRightMargin-0.01, legY = 1-style2.PadTopMargin-0.01;
  TLegend leg(legX-legW, legY-legH, legX, legY);
  TChain RAll("ntp1");
  RAll.Add("AWG82/ntuples/small/FitRAll_RunAll.root");

  int nbins[] = {20,20}, colors[] = {4, 1}; 
  double minX[] = {0, 0}, maxX[] = {0.35, 0.35};
  TString Variables[] = {"candPisoftP","candPisoftP"}, TagLeft[] = {"a)","b)"}, xTitle[] = {"p_{#pi^{#pm}} (GeV)","p_{#pi^{0}} (GeV)"};
  TString Units[] = {" MeV", " MeV"};
  double factor[] = {1000, 1000}; int digits[] = {0,0};
  TCut CutPad[2] = {"candDstarType==1","candDstarType==2"}, CutHist[2] = {"(MCType==8||MCType==10)","MCType==12"};
  //TCut CutPad[2] = {"candDstarType==1","candDstarType==2"}, CutHist[2] = {"(MCType==2||MCType==4)","MCType==6"};
  TCut CutAll = "candType==4"; 
  TCanvas can("can","Soft Pion momenta");
  can.Divide(2,1);
  TH1F* h[2][2]; double n[2][2];
  for(int pad=1; pad>=0; pad--){
    can.cd(pad+1);
    for(int hist=0; hist<2; hist++) {
      TString hname = "h"; hname += pad; hname += hist;
      h[pad][hist] = new TH1F(hname,"",nbins[pad],minX[pad],maxX[pad]);
      h[pad][hist]->SetLineWidth(2);
      h[pad][hist]->SetLineColor(colors[hist]);
      TString vari = Variables[pad]; vari += ">>"; vari += hname;
      TCut Cuts = CutAll; Cuts += CutPad[pad]; Cuts += CutHist[hist];
      RAll.Draw(vari,Cuts);
      n[pad][hist] = h[pad][hist]->Integral(); 
    }
    h[pad][0]->SetMaximum(h[pad][0]->GetMaximum()*1.3);
    h[pad][0]->Draw();
    TString yTitle = "Events/("; yTitle += RoundNumber(maxX[pad]-minX[pad],digits[pad],(double)nbins[pad]/factor[pad]); 
    yTitle += Units[pad]; yTitle += ")";
    style2.setTitles(h[pad][0],xTitle[pad],yTitle,TagLeft[pad]);
    h[pad][1]->Scale(n[pad][0]/n[pad][1]);
    h[pad][1]->Draw("same");
    if(pad==0){
      leg.SetTextSize(style2.LabelSize); leg.SetFillColor(0); leg.SetTextFont(style2.nFont);
      leg.SetBorderSize(0);
      TString label = "D*l#nu ("; label += (int)n[pad][0]; label += ")";
      leg.AddEntry(h[pad][0],label);
      label = "D*#tau#nu ("; label += (int)n[pad][1]; label += ")";
      leg.AddEntry(h[pad][1],label);
      leg.Draw();
    }
  }

  TString pName = "public_html/Syst_SoftPi.eps"; 
  can.SaveAs(pName);
  for(int pad=0; pad<2; pad++)
    for(int hist=0; hist<2; hist++) 
      h[pad][hist]->Delete();

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
