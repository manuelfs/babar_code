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
#include "TLegend.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void Final_Likelihood(){

  Styles style; style.setPadsStyle(-8); 
  style.CanvasH = 350; 
  style.PadRightMargin = 0.049; style.PadLeftMargin = 0.184; style.yTitleOffset = 0.92;
  style.applyStyle();
  TCanvas can("can","Likelihood scan");
  can.Divide(2,1); TPad *cPad;
  //TString xTitle[] = {"Yield(B#rightarrowD#tau#nu)", "Yield(B#rightarrowD*#tau#nu)"};
  TString xTitle[] = {"a", "b"};

  TString textName = "FitAll/fits/TextFinalIsoDataNe2x100.txt";
  double Yield[2][70], Error[70];
  ReadFitFile(textName, Yield, Error);
  double YieldFactor[] = {(Yield[0][1]+Yield[0][3]+Yield[0][6]+Yield[0][8])/Yield[0][1],
			   (Yield[0][2]+Yield[0][4]+Yield[0][5]+Yield[0][7])/Yield[0][2]};
  TH1F *hLikeli[2], *hLikeli2[2], *hParabola[2];
  int nBins = 200;
  TLine line; line.SetLineStyle(2);

  TFile fLikeli("FitAll/Errors/Pulls20DataNewx100_RunAllIso.root"); fLikeli.cd();
  for(int pad=0; pad<2; pad++){
    cPad = (TPad *)can.cd(pad+1);
    TString hName = "Likelihood"; hName += pad+1; 
    hLikeli2[pad] = (TH1F *)(fLikeli.Get(hName));
    int BinsLikeli = hLikeli2[pad]->GetNbinsX();
    double minLikeli = hLikeli2[pad]->GetMinimum();
    hName += "Expanded";
    hLikeli[pad] = new TH1F(hName,"", BinsLikeli, YieldFactor[pad]*hLikeli2[pad]->GetXaxis()->GetBinLowEdge(1),
			    YieldFactor[pad]*hLikeli2[pad]->GetXaxis()->GetBinLowEdge(BinsLikeli+1));
    for(int bin=1; bin<=BinsLikeli; bin++)
      hLikeli[pad]->SetBinContent(bin,hLikeli2[pad]->GetBinContent(bin)-minLikeli);
    hName += "Parab";
    double mu = YieldFactor[pad]*Yield[0][pad+1], sigma = YieldFactor[pad]*Error[pad+1], nSig = 3.6;;
    if(pad==0) mu *= 1.01; // Needed because the likelihood scan is for a slightly different fit
    double minX = mu-nSig*sigma, maxX = mu+nSig*sigma;
    hParabola[pad] = new TH1F(hName,"",nBins, minX, maxX);
    for(int bin=1; bin<=nBins; bin++){
      double x = hParabola[pad]->GetBinCenter(bin);
      hParabola[pad]->SetBinContent(bin,pow((x-mu)/sigma,2)/2.);
    }
    hParabola[pad]->SetLineWidth(3); hLikeli[pad]->SetLineWidth(3); 
    hParabola[pad]->SetLineColor(2); hLikeli[pad]->SetLineColor(4); 
    hParabola[pad]->SetLineStyle(3);
    hParabola[pad]->SetLabelOffset(0.01,"Y");
    hParabola[pad]->GetXaxis()->CenterTitle(true);
    hParabola[pad]->Draw("c");
    for(int ns=1; ns<=3; ns++){
      double Y = ns*ns/2.;
      line.DrawLine(minX, Y, maxX, Y);  
      line.DrawLine(mu-ns*sigma, 0, mu-ns*sigma, Y); line.DrawLine(mu+ns*sigma, 0, mu+ns*sigma, Y);
    }
    style.setTitles(hParabola[pad],xTitle[pad],"n");//"n^{2}_{#sigma}/2");
    hLikeli[pad]->Draw("c same");

  }
  double legW = 0.38, legH = 0.24;
  double legX = 1-style.PadRightMargin-0.3, legY = 1-style.PadTopMargin-0.01;
  TLegend leg(legX-legW, legY-legH, legX, legY);
  leg.SetTextSize(style.LabelSize); leg.SetFillColor(0); leg.SetTextFont(style.nFont);
  leg.SetBorderSize(0);
  leg.AddEntry(hLikeli[0],"-Log(Likelihood)");
  leg.AddEntry(hParabola[0],"Parabola");
  can.cd(1); leg.Draw();

  TString pName = "public_html/Stability_Likelihood.eps"; 
  can.SaveAs(pName);
  for(int pad=0; pad<2; pad++){
      hParabola[pad]->Delete();
      hLikeli[pad]->Delete();
  }
}


