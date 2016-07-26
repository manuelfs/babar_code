#include "TString.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;



void softPi0_Weights(){

  TCanvas c;
  int nbins = 13;
  double minX = 0.2, maxX = 6.7;
  TH1D *H_pi0Cor = new TH1D("H_pi0Cor","H_pi0Cor",nbins,minX,maxX);
  H_pi0Cor->SetBinContent(0,0.9210917);
  H_pi0Cor->SetBinContent(1,0.94641);
  H_pi0Cor->SetBinContent(2,0.9539677);
  H_pi0Cor->SetBinContent(3,0.9525576);
  H_pi0Cor->SetBinContent(4,0.9574195);
  H_pi0Cor->SetBinContent(5,0.9683621);
  H_pi0Cor->SetBinContent(6,0.9791807);
  H_pi0Cor->SetBinContent(7,0.9807274);
  H_pi0Cor->SetBinContent(8,0.9867322);
  H_pi0Cor->SetBinContent(9,0.9927288);
  H_pi0Cor->SetBinContent(10,1.002873);
  H_pi0Cor->SetBinContent(11,1.002569);
  H_pi0Cor->SetBinContent(12,1.014027);
  H_pi0Cor->SetBinContent(13,1.02584);
  H_pi0Cor->SetBinContent(14,1.072968);
  H_pi0Cor->SetBinError(0,0.008606788);
  H_pi0Cor->SetBinError(1,0.001542882);
  H_pi0Cor->SetBinError(2,0.001353615);
  H_pi0Cor->SetBinError(3,0.001474407);
  H_pi0Cor->SetBinError(4,0.001673845);
  H_pi0Cor->SetBinError(5,0.0019275);
  H_pi0Cor->SetBinError(6,0.002296575);
  H_pi0Cor->SetBinError(7,0.002842807);
  H_pi0Cor->SetBinError(8,0.003600341);
  H_pi0Cor->SetBinError(9,0.004615527);
  H_pi0Cor->SetBinError(10,0.006091183);
  H_pi0Cor->SetBinError(11,0.008257723);
  H_pi0Cor->SetBinError(12,0.01205016);
  H_pi0Cor->SetBinError(13,0.02023907);
  H_pi0Cor->SetBinError(14,0.04756652);
  H_pi0Cor->SetEntries(226053.7);
  H_pi0Cor->Draw();

  TF1 fFit("fFit","[0]*x+[1]",minX,maxX);
  H_pi0Cor->Fit(&fFit,"N");
  fFit.SetLineWidth(2); fFit.SetLineColor(4);

  c.SaveAs("babar_code/Systematics/softPi0_Fit.eps");
   
}


