#include <iostream>

#include "TH1D.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TPad.h"
#include "TMarker.h"

#include "styles.hpp"
#include "plot_results.hpp"

using namespace std;

int main(){
  TString s_babarST="BABAR (ST)", s_babarHT="BABAR (HT)", s_belleST="BELLE (ST)", s_belleHT="BELLE (HT)";
  TString s_lhcb="LHCb";
  //vector<TString> allNames({s_babarST, s_babarHT, s_belleST, s_belleHT, s_lhcb});
  vector<TString> allNames({s_babarHT, s_belleHT, s_lhcb});

  vector<PadResults> pads;

  //////////////////// Taunu ////////////////////////
  vector<Results> resTaunu;
  resTaunu.push_back(Results(s_babarST, 1.7, {0.8}, {0.2}));
  resTaunu.push_back(Results(s_babarHT, 1.83, {0.53, 0.49}, {0.24}));
  resTaunu.push_back(Results(s_belleST, 1.25, {0.28}, {0.27}));
  resTaunu.push_back(Results(s_belleHT, 0.72, {0.27, 0.25}, {0.11}));
  Results resTaunuSM("SM", 0.75, {0.10, 0.05});
  Results resTaunuAverage("Average", 1.06, {0.19});
  // pads.push_back(PadResults("#it{B}(B^{+}#rightarrow #tau^{+}#nu_{#tau}) [10^{-4}]", 0., 2.8, 
  // 			    resTaunu, resTaunuSM, resTaunuAverage));

  //////////////////// RD ////////////////////////
  vector<Results> resRD;
  resRD.push_back(Results(s_babarHT, 0.440, {0.058}, {0.042}));
  resRD.push_back(Results(s_belleHT, 0.375, {0.064}, {0.026}));
  Results resRDSM("SM", 0.297, {0.017});
  Results resRDAverage("Average", 0.391, {0.041},{0.028});
  pads.push_back(PadResults("R(D)", 0.24, 0.54, resRD, resRDSM, resRDAverage));

  //////////////////// RDs ////////////////////////
  vector<Results> resRDs;
  resRDs.push_back(Results(s_babarHT, 0.332, {0.024}, {0.018}));
  resRDs.push_back(Results(s_belleHT, 0.293, {0.038}, {0.015}));
  resRDs.push_back(Results(s_lhcb, 0.336, {0.027}, {0.030}));
  Results resRDsSM("SM", 0.252, {0.003});
  Results resRDsAverage("Average", 0.322, {0.018},{0.012});
  pads.push_back(PadResults("R(D#lower[-.1]{*})", 0.22, 0.41, resRDs, resRDsSM, resRDsAverage));


  //////////////////////////// Making plots ////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  //// Setting style
  //int colorSM = kGray+2, colorAverage = kRed;
  int colorSM = kGreen+2, colorAverage = kRed+1, colorRes = 1;
  float sideTextSize = 0.08, axisTextSize = 0.08;
  float lMargin = 0.02, rMargin = 0.02;
  float bMargin = 0.17, tMargin = 0.03;
  gStyle->SetOptStat(0);
  gStyle->SetPadBottomMargin(bMargin);
  gStyle->SetPadTopMargin(tMargin);
  gStyle->SetPadLeftMargin(lMargin);
  gStyle->SetPadRightMargin(rMargin);
  gStyle->SetTitleSize(axisTextSize,"x");     // Set the 2 axes title size
  gStyle->SetLabelSize(axisTextSize,"x");     // Set the 2 axes label size
  gStyle->SetNdivisions(506, "x");   
  gStyle->SetTitleOffset(1,"x");     
  gStyle->SetTitleFont(132,"xyz");          // Set the all 2 axes title font
  gStyle->SetLabelFont(132,"xyz");          // Set the all 2 axes label font
  gStyle->SetTextFont(132);                // Set global text font
  gStyle->SetPadTickX(1);             // Ticks at the top
  gStyle->SetPadTickY(0);             // No ticks at the right


  //// Creating canvas
  float nPads = pads.size(); 
  float sideW = 230, padW = 500;
  float canW = sideW + nPads*padW;
  TCanvas can("can","", canW, 500);


  //// Plotting results names
  float nRes = allNames.size();
  float resultH = (1-bMargin-tMargin)/(nRes);

  TLatex label; label.SetNDC(kTRUE); label.SetTextFont(132);
  label.SetTextSize(sideTextSize); label.SetTextAlign(12);
  for(size_t ind=0; ind<nRes; ind++){
    label.DrawLatex(0.01, bMargin + resultH*(nRes-ind-0.5), allNames[ind]);
  }

  TString hname;
  vector<TPad*> Pads;
  vector<TH1D*> hRD;
  TBox box; box.SetLineColor(0);
  TLine line; 
  TMarker marker; marker.SetMarkerStyle(8); marker.SetMarkerSize(1.6);
  float errH = 0.02;
  for(size_t pad=0; pad<nPads; pad++){
    //// Drawing pads
    can.cd(0);
    hname = "Pad"; hname += pad;
    Pads.push_back(new TPad(hname,"",(sideW+pad*padW)/canW, 0, (sideW+(pad+1)*padW)/canW, 1));
    Pads[pad]->Draw(); Pads[pad]->cd();

    //// Drawing base histogram
    hname = "hRD"; hname += pad;
    hRD.push_back(new TH1D(hname, "", 10, pads[pad].minX, pads[pad].maxX));
    hRD[pad]->GetYaxis()->SetNdivisions(0);
    hRD[pad]->SetMinimum(0);
    hRD[pad]->SetMaximum(1);
    hRD[pad]->GetXaxis()->CenterTitle(true);
    hRD[pad]->SetXTitle(pads[pad].title);
    hRD[pad]->Draw();

    //// Drawing SM
    box.SetFillColor(colorSM); box.SetFillColorAlpha(colorSM, 0.3);
    box.DrawBox(pads[pad].sm.value-pads[pad].sm.errDown(), 0,
		pads[pad].sm.value+pads[pad].sm.errUp(), 1);
    line.SetLineColor(colorSM); line.SetLineWidth(2);
    line.DrawLine(pads[pad].sm.value, 0, pads[pad].sm.value, 1);

    //// Drawing Average
    box.SetFillColor(colorAverage); box.SetFillColorAlpha(colorAverage, 0.3);
    box.DrawBox(pads[pad].average.value-pads[pad].average.errDown(), 0,
		pads[pad].average.value+pads[pad].average.errUp(), 1);
    line.SetLineColor(colorAverage); line.SetLineWidth(2);
    line.DrawLine(pads[pad].average.value, 0, pads[pad].average.value, 1);

    //// Drawing measuremnts
    line.SetLineColor(colorRes); line.SetLineWidth(2);
    for(auto &result : pads[pad].vresults){
      float resY=0;
      for(size_t ind=0; ind<nRes; ind++) 
	if(allNames[ind]==result.name) resY = (nRes-ind-0.5)/nRes;
      line.DrawLine(result.value - result.errDown(), resY, result.value + result.errUp(), resY);
      line.DrawLine(result.value - result.errDown(), resY-errH, result.value - result.errDown(), resY+errH);
      line.DrawLine(result.value + result.errUp(), resY-errH, result.value + result.errUp(), resY+errH);
      line.DrawLine(result.value - result.statDown(), resY-0.7*errH, result.value - result.statDown(), resY+0.7*errH);
      line.DrawLine(result.value + result.statUp(), resY-0.7*errH, result.value + result.statUp(), resY+0.7*errH);
      marker.DrawMarker(result.value, resY);
    }

   hRD[pad]->Draw("same axes");
  }

  TString plotName = "plots/results.pdf";
  can.SaveAs(plotName);

  return 1;
}

Results::Results(TString iname, float ivalue, vector<float> istatError, vector<float> isystError):
  name(iname),
  value(ivalue),
  statError(istatError),
  systError(isystError){
  if(statError.size()==1) statError.push_back(statError[0]);
  if(systError.size()==1) systError.push_back(systError[0]);
  }
  
PadResults::PadResults(TString ititle, float iminX, float imaxX, std::vector<Results> ivresults, 
		       Results ism, Results iaverage):
  title(ititle),
  minX(iminX),
  maxX(imaxX),
  vresults(ivresults),
  sm(ism),
  average(iaverage){

  }
