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
#include "babar_code/Styles/Styles.cc"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void Syst_Nopi0Yields(){

  Styles style; style.setPadsStyle(4); style.applyStyle();

  TCanvas can("can","Yields and no Dpi0 fit");
  can.Divide(2,2); TPad *cPad;
  TArrow arrow; arrow.SetLineWidth(2);
  double limits[4][2] = {{0,0.7},{0.1,0.7},{0.35,0.9},{0.3,0.8}};
  TString yTitleBeg = "Entries/(";
  TString xTitle = "Signal yield";
  TString PadLabel[] = {"a)", "b)","c)", "d)"};
  TString rLabel[] = {"D^{0}","D*^{0}","D^{+}","D*^{+}"};
  double YieldNoDpi0[2][4] = {{217.6, 272.5, 131.8, 167.9},{201.1, 311.4, 131.8, 167.9}}; // The second is with 0.11 convolution
  TFile PullFile("FitAll/Errors/Pulls3Data.root"); PullFile.cd();
  TH1F* h[4], *hRatio[4];
  for(int pad=0; pad<4; pad++){
    cPad = (TPad *)can.cd(pad+1);
    TString hName = "YieldFit"; hName += pad+1;
    h[pad] = (TH1F *)(PullFile.Get(hName));
    int nbins = h[pad]->GetNbinsX();
    int minbin = (int)(nbins*limits[pad][0]), maxbin = (int)(nbins*limits[pad][1]);
    double minx = h[pad]->GetXaxis()->GetBinLowEdge(minbin);
    double maxx = h[pad]->GetXaxis()->GetBinLowEdge(maxbin+1);
    hName += "Copy";
    hRatio[pad] = new TH1F(hName,"",maxbin-minbin+1,minx,maxx);
    for(int bin=minbin; bin<=maxbin; bin++){
      hRatio[pad]->SetBinContent(bin-minbin+1, h[pad]->GetBinContent(bin));
    }
    hRatio[pad]->Draw();
    TString yTitle = yTitleBeg; yTitle += RoundNumber(maxx-minx,1,(double)(maxbin-minbin+1)); yTitle+=")";
    style.setTitles(hRatio[pad],xTitle,yTitle,PadLabel[pad],rLabel[pad]);
    arrow.SetLineColor(2); arrow.SetFillColor(2); 
    arrow.DrawArrow(YieldNoDpi0[0][pad],hRatio[pad]->GetMaximum()*0.8,YieldNoDpi0[0][pad],0.2, 0.02);
    if(pad<2){
      arrow.SetLineColor(4); arrow.SetFillColor(4); 
      arrow.DrawArrow(YieldNoDpi0[1][pad],hRatio[pad]->GetMaximum()*0.8,YieldNoDpi0[1][pad],0.2, 0.02);
    }
  }

  TString pName = "public_html/Syst_Nopi0Yields.eps"; 
  can.SaveAs(pName);
  for(int pad=0; pad<4; pad++){
      h[pad]->Delete();
      hRatio[pad]->Delete();
  }
}

