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

void PID(TString Variable = "candPLep"){

  Styles style2; style2.setPadsStyle(2); style2.applyStyle();

  TChain newPID("ntp1");
  newPID.Add("AWG82/ntuples/small/RAll_RunAll.root");

  TCanvas can("can","PID comparison");
  can.Divide(2,1);
  can.cd(1);
  TH1F* h[4]; double n[4];
  double pmax = 3.6; if(Variable=="candPstarLep") pmax = 2.4;
  for(int i=0; i<4; i++) {
    TString hname = "h"; hname += i;
    h[i] = new TH1F(hname,"",24,0,pmax);
    h[i]->SetLineWidth(2);
    if(i%2==0) h[i]->SetLineColor(4);
  }
  TString cutsTau = "candPLep>0.4&&candIsMu==0&&(MCType==5||MCType==6||MCType==11||MCType==12)";
  TString cutsLep = "candPLep>0.4&&candIsMu==0&&(MCType==1||MCType==2||MCType==7||MCType==8)";
  TString var = Variable; var+=">>h0"; n[0] = newPID.Draw(var,cutsTau);
  var = Variable; var+=">>h1"; n[1] = newPID.Draw(var,cutsLep);
  h[0]->Scale(1000/h[0]->Integral());
  h[1]->Scale(1000/h[1]->Integral());
  h[0]->SetMaximum(h[0]->GetMaximum()*1.2);
  h[0]->Draw();
  TString xtitle = "p_{l} (GeV)"; if(Variable=="candPstarLep") xtitle = "p*_{l} (GeV)";
  style2.setTitles(h[0],xtitle,"Events/(150 MeV)","e");
  h[1]->Draw("same");
  double legW = 0.43, legH = 0.17;
  double legX = 1-style2.PadRightMargin-0.01, legY = 1-style2.PadTopMargin-0.01;
  TLegend leg(legX-legW, legY-legH, legX, legY);
  leg.SetTextSize(style2.LabelSize); leg.SetFillColor(0); leg.SetTextFont(style2.nFont);
  leg.SetBorderSize(0);
  TString label = "D#tau#nu ("; label += (int)n[0]; label += ")";
  leg.AddEntry(h[0],label);
  label = "Dl#nu ("; label += (int)n[1]; label += ")";
  leg.AddEntry(h[1],label);
  leg.Draw();

  can.cd(2);
  cutsTau = "candIsMu==1&&(MCType==5||MCType==6||MCType==11||MCType==12)";
  cutsLep = "candIsMu==1&&(MCType==3||MCType==4||MCType==9||MCType==10)";
  var = Variable; var+=">>h2"; n[2] = newPID.Draw(var,cutsTau);
  var = Variable; var+=">>h3"; n[3] = newPID.Draw(var,cutsLep);
  h[2]->Scale(1000/h[2]->Integral());
  h[3]->Scale(1000/h[3]->Integral());
  h[2]->SetMaximum(h[2]->GetMaximum()*1.2);
  h[2]->Draw();
  style2.setTitles(h[2],xtitle,"Events/(150 MeV)","#mu");
  h[3]->Draw("same");
  TLegend leg2(legX-legW, legY-legH, legX, legY);
  leg2.SetTextSize(style2.LabelSize); leg2.SetFillColor(0); leg2.SetTextFont(style2.nFont);
  leg2.SetBorderSize(0);
  label = "D#tau#nu ("; label += (int)n[2]; label += ")";
  leg2.AddEntry(h[2],label);
  label = "Dl#nu ("; label += (int)n[3]; label += ")";
  leg2.AddEntry(h[3],label);
  leg2.Draw();

  TString pName = "public_html/Jul6_PLep.eps"; 
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
