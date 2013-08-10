#include "TH1F.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TString.h"
#include <fstream>
#include <iostream>
using std::cout;
using std::endl;

TString RoundNumber(double n, int e, double d=1);

void PlotErrorRatio(TString fname){
  TFile f(fname);
  TCanvas c("c","Error on the Ratios",900,600);
  c.Divide(2,2);
  TLatex label;
  label.SetTextSize(0.065);label.SetNDC(kTRUE);
  TString titles[4] = {"D^{0}","D*^{0}","D^{+}","D*^{+}",};
  gStyle->SetOptStat(1110);
  TH1F *h[4];
  for(int i=0; i<4; i++){
    c.cd(i+1);
    TString hname = "Ratio"; hname += i;
    h[i] = (TH1F *)f.Get(hname);
    double mean = h[i]->GetMean(), rms = h[i]->GetRMS();
    int nBins = h[i]->GetNbinsX(), ibin = 1, inibin, finbin;
    while(h[i]->GetBinContent(ibin) == 0) ibin++;
    inibin = ibin-1; ibin = nBins;
    while(h[i]->GetBinContent(ibin) == 0) ibin--;
    finbin = ibin+1; 
    h[i]->GetXaxis()->SetRange(inibin,finbin);
    h[i]->SetTitle(""); h[i]->Draw();
    TString Title = titles[i]; Title+=" #tau #nu / "; Title += titles[i];
    Title += " l #nu: Error ";
    Title += RoundNumber(rms*100,1,mean); Title += " %";
    label.DrawLatex(0.11,0.92,Title);
  }
  c.SaveAs("babar_code/keysFit/ErrorRatio.eps");
  f.Close();
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
