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
#include "TH1F.h"
#include "TF1.h"
#include "babar_code/Styles/Styles.cc"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;


void Plot_Final_Gaussian(TString typeError){

  Styles style; style.setPadsStyle(-8); 
  style.CanvasH = 350;
  style.applyStyle();

  TCanvas can("can","RD from bootstrap");
  can.Divide(2,1); TPad *cPad;
  TString yTitleBeg = "Entries/(";
  TString xTitle[] = {"R(D)", "R(D*)"};
  TString PadLabel[] = {"a)", "b)"};
  TString PullName = "FitAll/Errors/Pulls"; PullName += typeError;
  PullName += "DataNewx100_RunAllIso.root";
  TFile PullFile(PullName); PullFile.cd();
  TH1F* h[2], *hRatio[2];
  TF1 *fGaus[2]; 
  for(int pad=0; pad<2; pad++){
    cPad = (TPad *)can.cd(pad+1);
    TString hName = "Ratio"; hName += pad+1;
    h[pad] = (TH1F *)(PullFile.Get(hName));
    double mean = h[pad]->GetMean(), rms = h[pad]->GetRMS();
    double minX = mean-3.5*rms, maxX = mean+3.5*rms;
    int nbins = h[pad]->GetNbinsX(); 
    double minH = h[pad]->GetBinLowEdge(1), maxH = h[pad]->GetBinLowEdge(nbins+1);
    double binW = (maxH-minH)/(double)nbins;
    int minbin = (int)((minX-minH)/binW+1);
    int maxbin = (int)((maxX-minH)/binW+1);

    hName += "Subset";
    hRatio[pad] = new TH1F(hName,"",maxbin-minbin+1,minX,maxX);
    for(int bin=minbin; bin<=maxbin; bin++){
      hRatio[pad]->SetBinContent(bin-minbin+1, h[pad]->GetBinContent(bin));
    }
    hRatio[pad]->GetXaxis()->CenterTitle(true);
    hRatio[pad]->SetLabelOffset(0.01,"Y");
    hRatio[pad]->Draw();
    hName = "Gaussian"; hName += pad+1;
    fGaus[pad] = new TF1(hName,"gaus",mean-4*rms, mean+4*rms);
    //fGaus[pad]->SetLineColor(2);
    fGaus[pad]->SetLineColor(4);
    hRatio[pad]->Fit(fGaus[pad],"L Q");
    double Sigma = fGaus[pad]->GetParameter(2);

    TString yTitle = yTitleBeg; yTitle += RoundNumber(maxX-minX,4,
						      (double)(maxbin-minbin+1)); yTitle+=")";
//     TString rLabel = "#splitline{RMS = "; rLabel += RoundNumber(rms,4);
//     rLabel += "}{   #sigma_{fit}  = "; rLabel += RoundNumber(Sigma,4); rLabel += "}";
    TString rLabel = "#sigma_{fit}  = "; rLabel += RoundNumber(Sigma,4); 
    style.setTitles(hRatio[pad],xTitle[pad],yTitle);
  }

  TString pName = "public_html/Stability_Gaussian_";
  if(typeError=="3") pName += "fDss";
  else pName += "PDFstat";
  pName += ".eps"; 
  can.SaveAs(pName);
  for(int pad=0; pad<2; pad++){
      h[pad]->Delete();
      hRatio[pad]->Delete();
      fGaus[pad]->Delete();
  }
}

void Final_Gaussian(){
  Plot_Final_Gaussian("3");
}
