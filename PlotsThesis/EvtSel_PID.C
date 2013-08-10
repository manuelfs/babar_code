#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
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

void EvtSel_PID(TString isMu = "1"){

  Styles style2; style2.setPadsStyle(2); style2.applyStyle();

  TChain newPID("ntp1"), oldPID("ntp1");
  oldPID.Add("AWG82/ntuples/PIDR26/Merged/SP-1235-BSemiExclAdd-Run4-R26a-v04.root");
  newPID.Add("AWG82/ntuples/R26/Merged/SP-1235-BSemiExclAdd-Run4-R26a-v04.root");

  TCanvas can("can","PID comparison");
  can.Divide(2,1);
  can.cd(1);
  TH1F* h[4]; double n[4];
  for(int i=0; i<4; i++) {
    TString hname = "h"; hname += i;
    h[i] = new TH1F(hname,"",24,0,3.6);
    h[i]->SetLineWidth(2);
    if(i%2==0) h[i]->SetLineColor(4);
    else h[i]->SetLineColor(2);
  }
  TString cuts = "candIsMu=="; cuts += isMu; cuts += "&&candLepTru==1&&MCType>0";
  n[0] = newPID.Draw("candPLep>>h0",cuts);
  n[1] = oldPID.Draw("candPLep>>h1",cuts);
  h[0]->SetMaximum(h[0]->GetMaximum()*1.2);
  h[0]->Draw();
  style2.setTitles(h[0],"p_{l} (GeV)","Events/(150 MeV)","a)");
  h[1]->Draw("same"); h[0]->Draw("same");
  double legW = 0.43, legH = 0.17;
  double legX = 1-style2.PadRightMargin-0.01, legY = 1-style2.PadTopMargin-0.01;
  TLegend leg(legX-legW, legY-legH, legX, legY);
  leg.SetTextSize(style2.LabelSize); leg.SetFillColor(0); leg.SetTextFont(style2.nFont);
  leg.SetBorderSize(0);
  TString label = "New PID ("; label += (int)n[0]; label += ")";
  leg.AddEntry(h[0],label);
  label = "Old PID ("; label += (int)n[1]; label += ")";
  leg.AddEntry(h[1],label);
  leg.Draw();

  can.cd(2);
  cuts = "candIsMu=="; cuts += isMu; cuts += "&&candLepTru==0&&MCType>0";
  n[2] = newPID.Draw("candPLep>>h2",cuts);
  n[3] = oldPID.Draw("candPLep>>h3",cuts);
  h[2]->SetMaximum(h[2]->GetMaximum()*1.2);
  h[2]->Draw();
  style2.setTitles(h[2],"p_{l} (GeV)","Events/(150 MeV)","b)");
  h[3]->Draw("same");h[2]->Draw("same");
  TLegend leg2(legX-legW, legY-legH, legX, legY);
  leg2.SetTextSize(style2.LabelSize); leg2.SetFillColor(0); leg2.SetTextFont(style2.nFont);
  leg2.SetBorderSize(0);
  label = "New PID ("; label += (int)n[2]; label += ")";
  leg2.AddEntry(h[2],label);
  label = "Old PID ("; label += (int)n[3]; label += ")";
  leg2.AddEntry(h[3],label);
  leg2.Draw();

  TString pName = "public_html/EvtSel_PID"; pName += isMu; pName += ".eps"; 
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
