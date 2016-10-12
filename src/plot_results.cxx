#include <iostream>
#include <getopt.h>

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
#include "keys_utils.hpp"

using namespace std;

namespace{
  bool printResults = true;
  enum whichPlots {bdxtaunu, btaunu, both};
  int whichPlot = bdxtaunu;
}

int main(int argc, char *argv[]){
  GetOptions(argc, argv);

  TString s_babarST="BABAR (ST)", s_babarHT="BABAR (HT)", s_belleST="Belle (ST)", s_belleHT="Belle (HT)";
  TString s_lhcb="LHCb";
  vector<TString> allNames({s_babarHT, s_belleHT, s_lhcb});
  if(whichPlot == btaunu) allNames = vector<TString>({s_babarST, s_babarHT, s_belleST, s_belleHT});
  if(whichPlot == both) allNames = vector<TString>({s_babarST, s_babarHT, s_belleST, s_belleHT, s_lhcb});
  vector<PadResults> pads;

  //////////////////// Taunu ////////////////////////
  vector<Results> resTaunu;
  resTaunu.push_back(Results(s_babarST, 1.7, {0.8}, {0.2}));
  resTaunu.push_back(Results(s_babarHT, 1.83, {0.53, 0.49}, {0.24}));
  resTaunu.push_back(Results(s_belleST, 1.25, {0.28}, {0.27}));
  resTaunu.push_back(Results(s_belleHT, 0.72, {0.27, 0.25}, {0.11}));
  Results resTaunuSM("SM", 0.75, {0.10, 0.05});
  Results resTaunuAverage("Average", 1.06, {0.19});
  float maxTaunu = (printResults?4.4:2.8);
  if(whichPlot == btaunu || whichPlot == both) 
    pads.push_back(PadResults("#it{B}(B^{-}#rightarrow #tau^{-}#nu_{#tau}) [10^{-4}]", 0., maxTaunu, 
			      resTaunu, resTaunuSM, resTaunuAverage));

  //////////////////// RD ////////////////////////
  vector<Results> resRD;
  resRD.push_back(Results(s_babarHT, 0.440, {0.058}, {0.042}));
  resRD.push_back(Results(s_belleHT, 0.375, {0.064}, {0.026}));
  Results resRDSM("SM", 0.297, {0.017});
  Results resRDAverage("Average", 0.391, {0.041},{0.028});
  float maxRD = (printResults?0.82:0.54);
  if(whichPlot == bdxtaunu || whichPlot == both) 
    pads.push_back(PadResults("R(D)", 0.24, maxRD, resRD, resRDSM, resRDAverage));

  //////////////////// RDs ////////////////////////
  vector<Results> resRDs;
  resRDs.push_back(Results(s_babarHT, 0.332, {0.024}, {0.018}));
  resRDs.push_back(Results(s_belleHT, 0.293, {0.038}, {0.015}));
  resRDs.push_back(Results(s_lhcb, 0.336, {0.027}, {0.030}));
  Results resRDsSM("SM", 0.252, {0.003});
  Results resRDsAverage("Average", 0.322, {0.018},{0.012});
  float maxRDs = (printResults?0.56:0.4);
  if(whichPlot == bdxtaunu || whichPlot == both) 
    pads.push_back(PadResults("R(D#lower[-.1]{*})", 0.22, maxRDs, resRDs, resRDsSM, resRDsAverage));


  //////////////////////////// Making plots ////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  //// Setting style
  //int colorSM = kGray+2, colorAverage = kRed;
  int colorSM = kGreen+2, colorAverage = kRed+1, colorRes = 1;
  float sideTextSize = 0.188, axisTextSize = 0.08;
  float lMargin = 0.02, rMargin = 0.02;
  float bMargin = 0.2, tMargin = 0.03;
  gStyle->SetOptStat(0);
  gStyle->SetPadBottomMargin(bMargin);
  gStyle->SetPadTopMargin(tMargin);
  gStyle->SetPadLeftMargin(lMargin);
  gStyle->SetPadRightMargin(rMargin);
  gStyle->SetTitleSize(axisTextSize,"x");     // Set the 2 axes title size
  gStyle->SetLabelSize(axisTextSize,"x");     // Set the 2 axes label size
  gStyle->SetNdivisions(506, "x");   
  gStyle->SetTitleOffset(1.17,"x");     
  gStyle->SetTitleFont(132,"xyz");          // Set the all 2 axes title font
  gStyle->SetLabelFont(132,"xyz");          // Set the all 2 axes label font
  gStyle->SetTextFont(132);                // Set global text font
  gStyle->SetPadTickX(1);             // Ticks at the top
  gStyle->SetPadTickY(0);             // No ticks at the right


  //// Creating canvas
  float nPads = pads.size(); 
  float sideW = 105, padW = (printResults?300:250);
  float canW = sideW + nPads*padW;
  TCanvas can("can","", canW, 250);


  //// Plotting results names
  float nRes = allNames.size();
  float resultH = (1-bMargin-tMargin)/(nRes);
  TPad sidePad("sidePad", "", 0, 0, sideW/canW, 1); 
  sidePad.Draw(); sidePad.cd();

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
  TMarker marker; marker.SetMarkerStyle(8); marker.SetMarkerSize(0.8);
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
    hRD[pad]->Draw("");

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
    line.SetLineColor(colorRes); 
    for(auto &result : pads[pad].vresults){
      float resY=0;
      for(size_t ind=0; ind<nRes; ind++) 
	if(allNames[ind]==result.name) resY = (nRes-ind-0.5)/nRes;
      line.SetLineWidth(1);
      line.DrawLine(result.value - result.errDown(), resY, result.value + result.errUp(), resY);
      line.SetLineWidth(1);
      line.DrawLine(result.value - result.errDown(), resY-errH, result.value - result.errDown(), resY+errH);
      line.DrawLine(result.value + result.errUp(), resY-errH, result.value + result.errUp(), resY+errH);
      line.DrawLine(result.value - result.statDown(), resY-errH, result.value - result.statDown(), resY+errH);
      line.DrawLine(result.value + result.statUp(), resY-errH, result.value + result.statUp(), resY+errH);
      marker.DrawMarker(result.value, resY);

      int digits = 3;
      if(!pads[pad].title.Contains("D")){
	if(result.name == s_babarST) digits = 1;
	else digits = 2;
      }
      TString s_res = RoundNumber(result.value,digits)+" ";
      if(result.statUp() == result.statDown()) s_res += "#pm "+RoundNumber(result.statUp(),digits)+" ";
      else s_res += "^{+"+RoundNumber(result.statUp(),digits)+"}_{-"+RoundNumber(result.statDown(),digits)+"} ";
      if(result.systUp() == result.systDown()) s_res += "#pm "+RoundNumber(result.systUp(),digits);
      else s_res += "^{+"+RoundNumber(result.systUp(),digits)+"}_{-"+RoundNumber(result.systDown(),digits)+"}";
      label.SetTextSize(axisTextSize/1.2); label.SetTextAlign(32); label.SetNDC(kFALSE); 
      if(printResults) label.DrawLatex(pads[pad].maxX-0.03*(pads[pad].maxX-pads[pad].minX), resY, s_res);
    }

   hRD[pad]->Draw("same axis");
  }

  TString plotName = "plots/results";
  if(!printResults) plotName += "_notext";
  if(whichPlot == bdxtaunu || whichPlot == both) plotName += "_dxtaunu";
  if(whichPlot == btaunu || whichPlot == both) plotName += "_taunu";
  plotName += ".pdf";
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


void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"notext", no_argument, 0, 'n'},   
      {"plots", required_argument, 0, 'p'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "p:n", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'n':
      printResults = false;
      break;
    case 'p':
      whichPlot = atoi(optarg);
      break;
    case 0:
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
