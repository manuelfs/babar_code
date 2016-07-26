#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TArrow.h"
#include "TChain.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFile.h"
#include "TH2F.h"
#include "babar_code/Styles/Styles.cc"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void Syst_RDBootstrap(){

  Styles style; style.setPadsStyle(); style.applyStyle();

  TCanvas can("can","RD from bootstrap");
  TFile PullFile("FitAll/Errors/Pulls5DataNewx100_RunAllIso.root"); PullFile.cd();
  TH2F *hFile;
  hFile = (TH2F *)(PullFile.Get("Correl0"));
  double limitsR[2][2] = {{0.38,0.53},{0.305,0.37}}, nomR[] = {0.444, 0.3332};
  int minbin[2], maxbin[2];
  for(int pad=0; pad<2; pad++){
    int nbins = hFile->GetNbinsX(); 
    double minH = hFile->GetXaxis()->GetBinLowEdge(1), maxH = hFile->GetXaxis()->GetBinLowEdge(nbins+1);
    if(pad==1) {
      nbins = hFile->GetNbinsY(); 
      minH = hFile->GetYaxis()->GetBinLowEdge(1); maxH = hFile->GetYaxis()->GetBinLowEdge(nbins+1);
    }
    double binW = (maxH-minH)/(double)nbins;
    int binShift = (int)((nomR[pad]-hFile->GetMean(pad+1))/binW+0.5);
    minbin[pad] = (int)((limitsR[pad][0]-minH)/binW-binShift+1);
    maxbin[pad] = (int)((limitsR[pad][1]-minH)/binW-binShift+1);
    cout<<"binShift "<<binShift<<", minbin "<<minbin[pad]<<", maxbin "<<maxbin[pad]<<", binW "<<binW<<endl;
  }
  TH2F hRatio("hRatio","",maxbin[0]-minbin[0]+1,limitsR[0][0],limitsR[0][1],
	      maxbin[1]-minbin[1]+1,limitsR[1][0],limitsR[1][1]);
  for(int xbin=minbin[0]; xbin<=maxbin[0]; xbin++)
    for(int ybin=minbin[1]; ybin<=maxbin[1]; ybin++)
      hRatio.SetBinContent(xbin-minbin[0]+1, ybin-minbin[1]+1, hFile->GetBinContent(xbin,ybin));
  hRatio.SetXTitle("R(D)"); hRatio.SetYTitle("R(D*)"); 
  hRatio.Draw("cont0");
  TLatex label; label.SetTextSize(0.068); label.SetTextFont(132); label.SetTextAlign(19);
  TString RMSx = "RMS(R(D))   = "; RMSx += RoundNumber(hFile->GetRMS(1),3);
  TString RMSy = "RMS(R(D*)) = "; RMSy += RoundNumber(hFile->GetRMS(2),3);
  TString rho  = "#rho = "; rho += RoundNumber(hFile->GetCorrelationFactor(),2);
  label.DrawLatex(limitsR[0][1]-0.045, limitsR[1][1]-0.006, RMSx);
  label.DrawLatex(limitsR[0][1]-0.045, limitsR[1][1]-0.0125, RMSy);
  label.DrawLatex(limitsR[0][1]-0.045, limitsR[1][1]-0.019, rho);

  TString pName = "public_html/Syst_RDBootstrap.eps"; 
  can.SaveAs(pName);
  hFile->Delete();
}

