#include "TString.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TSystem.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

TString RoundNumber(double n, int e, double d=1);

void AddPulls(TString nPull = "1", int beg=1, int end=4){
  TString hName = "Pull"; hName += nPull;
  TString sumName = "hSum"; sumName += nPull;
  TH1F *hPull[4], *hSum = new TH1F(sumName,"",200,-4,4);
  for(int i=beg; i<=end; i++){
    TString fName = "FitAll/fits/Pulls"; fName += i; fName += "RAll.root";
    TFile f(fName); f.cd();
    hPull[i-1] = (TH1F*)gDirectory->Get(hName);
    hSum->Add(hPull[i-1]);
    f.Close();
  }
  TCanvas c;
  TF1 fGaus("Gaussian","gaus",-4,4);
  hSum->Fit(&fGaus,"L Q");
  double Mean = fGaus.GetParameter(1), eMean = fGaus.GetParError(1);
  double Sigma = fGaus.GetParameter(2), eSigma = fGaus.GetParError(2);
  cout<<hName<<" has "<<hSum->GetEntries()<<" entries: \t Mean is "<<RoundNumber(Mean,3)<<" +- "<<
    RoundNumber(eMean,3)<<" and Sigma "<<RoundNumber(Sigma,3)<<" +- "<<RoundNumber(eSigma,3)<<endl;
   hSum->Draw();
   c.SaveAs("canvas.eps");
   hSum->Delete();
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



