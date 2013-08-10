#include "TString.h"
#include "TFile.h"
#include "TLine.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH2F.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

TString RoundNumber(double n, int e, double d=1);

void plotCrossV(TString hname){
  double min = 999;
  TCanvas c("c", "Cross-validation values");
  
  TFile hfile(hname); hfile.cd();
  TH2F *cross = (TH2F *)gDirectory->Get("cross");
  double ini_m = cross->GetXaxis()->GetXmin(), fin_m = cross->GetXaxis()->GetXmax();
  double ini_p = cross->GetYaxis()->GetXmin(), fin_p = cross->GetYaxis()->GetXmax();
  int minBin = cross->GetMinimumBin(), mbin = cross->GetNbinsX(), pbin = cross->GetNbinsY();
  double bw_m = (fin_m-ini_m)/(double)(mbin), bw_p = (fin_p-ini_p)/(double)(pbin);
  ini_m += bw_m/2, fin_m -= bw_m/2; // I had the edges, and I wanted the centers
  ini_p += bw_p/2, fin_p -= bw_p/2;
  int minBin_m = minBin % (mbin+2), minBin_p = (minBin / (mbin+2));
  double min_m = ini_m + bw_m*(double)(minBin_m-1);
  double min_p = ini_p + bw_p*(double)(minBin_p-1);

  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  cross->SetXTitle("m_{miss} smoothing");
  cross->SetYTitle("p*_{l} smoothing factor");
  TString title = "Cross-validation minimum at s_{m} = "; title += RoundNumber(min_m,2);
  title += " and s_{p} = "; title += RoundNumber(min_p,2);
  cross->SetTitle(title);
  cross->Draw("cont4z");

  double len_m = (fin_m-ini_m)/20.;
  double len_p = (fin_p-ini_p)/20.;
  TLine line;
  line.SetLineColor(1);
  line.DrawLine(min_m-len_m, min_p, min_m+len_m, min_p);
  line.DrawLine(min_m, min_p-len_p, min_m, min_p+len_p);

  hname.Remove(0,hname.Last('/')+1);
  hname.Remove(hname.First('.'),hname.Length());
  TString epsName = "babar_code/keysFit/eps/"; epsName += hname; epsName+=".eps";
  c.SaveAs(epsName);
}

TString RoundNumber(double n, int e, double d){
  if(d==0) return " - ";
  double b = (int)(n/d*pow(10.,(double)e)+0.5);
  b /= pow(10.,(double)e);
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(!result.Contains(".") && e != 0) result += ".";
  
  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<e-afterdot.Length(); i++)
    result += "0";
  return result;
}
