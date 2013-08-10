#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TString.h"
#include "TRandom3.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void Diff_Ratio(double mean=15){
  TRandom3 rand(0);
  gStyle->SetOptFit(1);
  TCanvas c("c","Difference Vs Ratio",600,220);
  c.Divide(3);
  int min = (int)(mean-3*sqrt(mean));
  int max = (int)(mean+3*sqrt(mean));
  if(min<0) min=0;
  TH1F hPoi("Poisson","",(int)(max-min),min,max);
  TH1F hDiff("Difference","",(int)(3*mean),-1.5*mean,1.5*mean);
  TH1F hRat("Ratio","",(int)(10*sqrt(mean)),0,3);
  hPoi.SetLineColor(4);hDiff.SetLineColor(4);hRat.SetLineColor(4);
  hPoi.SetLineWidth(2);hDiff.SetLineWidth(2);hRat.SetLineWidth(2);
  for(int i=0; i<1000000; i++){
    double r1 = rand.Poisson(mean), r2 = rand.Poisson(mean);
    hPoi.Fill(r1);hPoi.Fill(r2);
    hDiff.Fill(r1-r2);
    if(r2) hRat.Fill(r1/r2);
    else hRat.Fill(4);
  }
  TF1 fDiff("fDiff","gaus",-1.5*mean,1.5*mean);
  TF1 fRat("fRat","gaus",0,3);
  c.cd(1);hDiff.Fit(&fDiff,"L Q");
  hRat.Fit(&fRat,"L Q");
  c.cd(1); hPoi.Draw();
  c.cd(2); hDiff.Draw();
  c.cd(3); hRat.Draw();
  TString epsName = "public_html/DiffRat_"; epsName += (int)mean; epsName += ".eps";
  c.SaveAs(epsName);
}
