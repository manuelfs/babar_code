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

void Syst_Dsspi(){

  Styles style2; style2.setPadsStyle(3); style2.CanvasH = 250;
style2.applyStyle();

  double legW = 0.55, legH = 0.15;
  double legX = 1-style2.PadRightMargin-0.02, legY = 1-style2.PadTopMargin-0.02;
  TLegend *leg[2];
  for(int ileg=0; ileg<2; ileg++) leg[ileg] = new TLegend(legX-legW, legY-legH, legX, legY);
  TChain RAll("ntp1");
  RAll.Add("AWG82/ntuples/small/RAll_RunAll.root");

  int nbins[] = {28,28,27}, colors[] = {4, 1}; 
  double minX[] = {-1, -2, 5.205}, maxX[] = {9, 7, 5.29};
  TString Variables[] = {"candM2","mm2pi0","candMES"}, TagLeft[] = {"a)","b)","c)"};
  TString xTitle[] = {"m^{2}_{miss} (GeV^{2})","m^{2}_{miss} (GeV^{2})","m_{ES} (GeV)"};
  TString Units[] = {" GeV^{2}", " GeV^{2}", " MeV"};
  double factor[] = {1, 1, 1000}; int digits[] = {2,2,0};
  TCut baseCut[] = {MvaAll, dssMvaAll, dssMvaAll}, Cuts;
  TString addCut[] = {"MCType==14&&MCPions==10","MCType==14&&MCPions==1"};

  TCanvas can("can","D** decays");
  can.Divide(3,1);
  TH1F* h[3][2]; double n[3][2];
  for(int pad=0; pad<3; pad++){
    can.cd(pad+1);
    for(int hist=0; hist<2; hist++) {
      TString hname = "h"; hname += pad; hname += hist;
      h[pad][hist] = new TH1F(hname,"",nbins[pad],minX[pad],maxX[pad]);
      h[pad][hist]->SetLineWidth(2);
      h[pad][hist]->SetLineColor(colors[hist]);
      TString vari = Variables[pad]; vari += ">>"; vari += hname;
      Cuts = baseCut[pad]; Cuts += addCut[hist];
      if(pad<2) Cuts += "candMES>5.27";
      RAll.Draw(vari,Cuts);
      n[pad][hist] = h[pad][hist]->Integral(); 
    }
    h[pad][0]->SetMaximum(h[pad][0]->GetMaximum()*1.1);
    h[pad][0]->Draw();
    TString yTitle = "Events/("; yTitle += RoundNumber(maxX[pad]-minX[pad],digits[pad],(double)nbins[pad]/factor[pad]); 
    yTitle += Units[pad]; yTitle += ")";
    style2.setTitles(h[pad][0],xTitle[pad],yTitle,TagLeft[pad]);
    h[pad][1]->Scale(n[pad][0]/n[pad][1]);
    h[pad][1]->Draw("same");
    if(pad<2){
      leg[pad]->SetTextSize(0.06); leg[pad]->SetFillColor(0); leg[pad]->SetTextFont(style2.nFont);
      leg[pad]->SetBorderSize(0);
      TString label = "D** #rightarrow D #pi^{0} ("; label += (int)n[pad][0]; label += ")";
      leg[pad]->AddEntry(h[pad][0],label);
      label = "D** #rightarrow D #pi^{#pm} ("; label += (int)n[pad][1]; label += ")";
      leg[pad]->AddEntry(h[pad][1],label);
      leg[pad]->Draw();
    }
  }

  TString pName = "public_html/Syst_Dsspi.eps";  
  can.SaveAs(pName);
  for(int pad=0; pad<3; pad++)
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
