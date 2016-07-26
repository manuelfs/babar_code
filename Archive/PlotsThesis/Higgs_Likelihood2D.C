#include "babar_code/PlotsThesis/PlotUtils.cc"
#include "babar_code/Styles/Styles.cc"
#include "TString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMatrixT.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include <fstream>
#include <iostream>

#define Nsig 9
using namespace std;
using std::cout;
using std::endl;

void Higgs_Likelihood2D(){

  Styles style; style.setPadsStyle(1); //style.PadLeftMargin
  style.applyStyle();
  TCanvas can("can","Likelihood scan");
  TPad *cPad = (TPad *)can.cd(0);
  cPad->SetLeftMargin(0.12);
  TString xTitle[] = {"Yield(B#rightarrowD#tau#nu)", "Yield(B#rightarrowD*#tau#nu)"};

  TString textName = "FitAll/fits/TextFinalIsoDataNe2x100.txt";
  double Yield[2][70], Error[70], RD[2][2] = {{0.4441, 0.3332}, {0.297, 0.252}};
  ReadFitFile(textName, Yield, Error);

  TFile fLikeli("FitAll/Errors/Pulls21DataNewx100_RunAllIso_61x61.root"); fLikeli.cd();
  TH2F *hLikeli;
  hLikeli = (TH2F *)(fLikeli.Get("Likeli2D_0"));
  double minLikeli = hLikeli->GetMinimum();
  int nBins[] = {hLikeli->GetNbinsX(), hLikeli->GetNbinsY()};
  double Limit[2][2] = {{hLikeli->GetXaxis()->GetBinLowEdge(1), hLikeli->GetXaxis()->GetBinLowEdge(nBins[0]+1)},
			{hLikeli->GetYaxis()->GetBinLowEdge(1), hLikeli->GetYaxis()->GetBinLowEdge(nBins[1]+1)}};

  TH2F hFinal("hFinal","",nBins[0],Limit[0][0],Limit[0][1], nBins[1], Limit[1][0],Limit[1][1]);
  for(int bin=1; bin<=nBins[0]; bin++)
    for(int ybin=1; ybin<=nBins[1]; ybin++)
	hFinal.SetBinContent(bin,ybin, hLikeli->GetBinContent(bin,ybin)-minLikeli);

  double zCont[Nsig];
  int colors[] = {kBlue+3, kBlue+2, kBlue+1, kBlue-3, kBlue-4, kBlue-7, kBlue-9, kBlue-10, 0};
  gStyle->SetPalette(Nsig, colors);
  for(int ns=0; ns<Nsig; ns++) zCont[ns] = ns*ns/2.;
  //for(int ns=Nsig; ns>=1; ns--) zCont[ns-1] = ns*ns/2.;
  hFinal.SetContour(Nsig,zCont);
  hFinal.SetXTitle(xTitle[0]); hFinal.SetYTitle(xTitle[1]); 
  hFinal.Draw("cont4");

  double padL[2][2] = {{cPad->GetLeftMargin(),   1-cPad->GetRightMargin()}, // Margins are reversed for x-y!
		       {cPad->GetBottomMargin(), 1-cPad->GetTopMargin()}};
  double SMyield[2], frac=0.018;
  TLine line; line.SetLineWidth(2); line.SetLineColor(1);
  for(int chan=0; chan<2; chan++){
    double range = Limit[chan][1]-Limit[chan][0];
    SMyield[chan] = Yield[0][chan+1]*RD[1][chan]/RD[0][chan];
    SMyield[chan] = padL[chan][0] + (padL[chan][1]-padL[chan][0])*(SMyield[chan]-Limit[chan][0])/range;
  }
  line.DrawLineNDC(SMyield[0]-frac, SMyield[1], SMyield[0]+frac, SMyield[1]);
  line.DrawLineNDC(SMyield[0], SMyield[1]-2*frac, SMyield[0], SMyield[1]+2*frac);

  TLatex latex; latex.SetNDC(kTRUE); latex.SetTextAlign(33); latex.SetTextSize(style.TextSize);
  latex.DrawLatex(SMyield[0]-frac,SMyield[1]-frac,"SM");

  line.SetLineColor(0);
  for(int chan=0; chan<2; chan++){
    double range = Limit[chan][1]-Limit[chan][0];
    SMyield[chan] = Yield[0][chan+1];
    SMyield[chan] = padL[chan][0] + (padL[chan][1]-padL[chan][0])*(SMyield[chan]-Limit[chan][0])/range;
  }
  line.DrawLineNDC(SMyield[0]-frac, SMyield[1], SMyield[0]+frac, SMyield[1]);
  line.DrawLineNDC(SMyield[0], SMyield[1]-2*frac, SMyield[0], SMyield[1]+2*frac);


  TH1F *histo[Nsig];
  double legW = 0.06, legH = 0.05*Nsig;
  double legX = 1-style.PadRightMargin-0.02, legY = 1-style.PadTopMargin-0.02;
  TLegend *leg = new TLegend(legX-legW, legY-legH, legX, legY);
  leg->SetTextSize(0.06); leg->SetFillColor(0); leg->SetTextFont(style.nFont);
  for(int ileg=0; ileg<Nsig-1; ileg++) {
    TString label = "histo"; label += ileg+1; 
    histo[ileg] = new TH1F(label,"histo",10,0,10);
    histo[ileg]->SetLineColor(colors[ileg]);histo[ileg]->SetFillColor(colors[ileg]);
    label = " "; label += ileg+1; label += "#sigma";
    leg->AddEntry(histo[ileg],label);
  }
  leg->Draw();



  TString pName = "public_html/Higgs_Likelihood2D.eps"; 
  can.SaveAs(pName);
  fLikeli.Close();
}


