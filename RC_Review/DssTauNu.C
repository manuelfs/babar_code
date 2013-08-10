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

void DssTauNu(int isDpi0=1){

  Styles style2; style2.setPadsStyle(2); style2.applyStyle();

  double legW = 0.37, legH = 0.17;
  double legX = 1-style2.PadRightMargin-0.01, legY = 1-style2.PadTopMargin-0.02;
  TChain RAll("ntp1");
  RAll.Add("AWG82/ntuples/small/FitRAllNewx100_RunAll.root");

  TString label, Decay[] = {"D**l#nu (", "D**#tau#nu (", "D(*)#tau#nu ("};
  int nbins[] = {25,15}, colors[] = {4, 1, 2}, nHist = 2; 
  double minX[] = {-2.5, 0}, maxX[] = {7.5, 2.4}, maxH[] = {0,0};
  TString Variables[] = {"candM2","candPstarLep"}, TagLeft[] = {"",""}, xTitle[] = {"m^{2}_{miss} (GeV^{2})","p*_{l} (GeV)"};
  TString Units[] = {" GeV^{2}", " MeV"};
  double factor[] = {1, 1000}; int digits[] = {2,0};
  TString CutHist[3] = {"MCTaumode<0&&MCType==14","MCTaumode>=0&&MCType==14","(MCType==5||MCType==6||MCType==11||MCType==12)"};
  TString CutAll = "weight*(candType>4&&"; 
  if(isDpi0==0) {
    CutAll = "weight*(candType<=4&&"; 
    minX[0] = -1; maxX[0] = 9.5; legH = 0.25;
    nHist = 3; nbins[1] = 20;
  }
  TLegend leg(legX-legW, legY-legH, legX, legY);
  TCanvas can("can","D**taunu background");
  can.Divide(2,1);
  TH1F* h[2][3]; double n[2][3];
  for(int pad=1; pad>=0; pad--){
    can.cd(pad+1);
    for(int hist=0; hist<nHist; hist++) {
      TString hname = "h"; hname += pad; hname += hist;
      h[pad][hist] = new TH1F(hname,"",nbins[pad],minX[pad],maxX[pad]);
      h[pad][hist]->SetLineWidth(2);
      h[pad][hist]->SetLineColor(colors[hist]);
      TString Cuts = CutAll; Cuts += CutHist[hist]; Cuts += ")";
      RAll.Project(hname,Variables[pad],Cuts);
      n[pad][hist] = h[pad][hist]->Integral(); 
      h[pad][hist]->Scale(1000/n[pad][hist]);
      if(h[pad][hist]->GetMaximum() > maxH[pad]) maxH[pad] = h[pad][hist]->GetMaximum();
    }
  }
  for(int pad=1; pad>=0; pad--){
    can.cd(pad+1);
    h[pad][0]->SetMaximum(maxH[pad]*1.1);
    h[pad][0]->Draw();
    TString yTitle = "Events/("; yTitle += RoundNumber(maxX[pad]-minX[pad],digits[pad],(double)nbins[pad]/factor[pad]); 
    yTitle += Units[pad]; yTitle += ")";
    style2.setTitles(h[pad][0],xTitle[pad],yTitle,TagLeft[pad]);
    h[pad][1]->Draw("same");
    if(isDpi0==0) h[pad][2]->Draw("same");
    if(pad==0){
      leg.SetTextSize(style2.LabelSize); leg.SetFillColor(0); leg.SetTextFont(style2.nFont);
      leg.SetBorderSize(0);
      for(int hist=0; hist<nHist; hist++) {
	if(hist==2) n[pad][hist] *= 1.342;
	label = Decay[hist]; label += (int)n[pad][hist]; label += ")";
	leg.AddEntry(h[pad][hist],label);
      }
      leg.Draw();
    }
  }
  TString pName = "public_html/DssTauNu"; pName += isDpi0; pName += ".eps"; 
  can.SaveAs(pName);
  for(int pad=0; pad<2; pad++)
    for(int hist=0; hist<nHist; hist++) 
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
