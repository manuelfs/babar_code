#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TBox.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TChain.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFile.h"
#include "TMath.h"
#include "babar_code/Styles/Styles.cc"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void ResultsRD(){

  Styles style; style.setPadsStyle(2); style.applyStyle();
  gStyle->SetPadLeftMargin(0.3); gStyle->SetPadTopMargin(0.12);

  double minRD[2] = {99., 99.}, maxRD[2] = {-99., -99.}, chi2[] = {0,0}, mean[] = {0,0}, RMS[] = {0,0}; 
  double nomiRD[2][2] = {{0.4402,0.0715},{0.332,0.0299}};
  double SMRD[2][2] = {{0.297, 0.017}, {0.252, 0.003}}, maxHisto = 1.;
  int ndof[2] = {-1,-1};
  int colors[] = {28,4,28,28,4};
  int nFiles = 5; TString fileTag[] = {"Belle 2007", "BaBar 2008", "Belle 2009", "Belle 2010", "BaBar 2012"};
  double RD[5][2][2] = {{{0.42, 0.13},{0.436, 0.119}}, {{0.416, 0.13},{0.297, 0.06}}, {{0.592, 0.163},{0.475, 0.104}}, 
			{{0.345,0.114},{0.426,0.085}}, {{0.4402,0.0715},{0.3316,0.0299}} };
  for(int file=0; file<nFiles; file++){
    for(int cand=0; cand<2; cand++){
      if(RD[file][cand][0]-RD[file][cand][1] < minRD[cand%2]) minRD[cand%2] = RD[file][cand][0]-RD[file][cand][1];
      if(RD[file][cand][0]+RD[file][cand][1] > maxRD[cand%2]) maxRD[cand%2] = RD[file][cand][0]+RD[file][cand][1];
      if(file<4 && (file>0||cand==1)){
	double Val = RD[file][cand][0], eVal = RD[file][cand][1];
	mean[cand] += Val/eVal/eVal;
	RMS[cand] += 1/eVal/eVal;
      }
    }
  }
  for(int cand=0; cand<2; cand++){
    mean[cand] /= RMS[cand]; RMS[cand] = sqrt(1/RMS[cand]);
    nomiRD[cand][0] = mean[cand]; nomiRD[cand][1] = RMS[cand]; 
  }
  nomiRD[0][0] = 0.4184; nomiRD[0][1] = 0.0810; // Found with babar_code/PlotsThesis/AverErrors.C
  nomiRD[1][0] = 0.3534; nomiRD[1][1] = 0.0474;
  //maxRD[0] = 0.78; maxRD[1] = 0.57; 
  TCanvas can("can","RD results");
  can.Divide(2,1); TPad *cPad;
  TString xTitle[] = {"R(D)", "R(D*)"};
  TH1F* hRD[2];
  TBox box; box.SetFillStyle(3002); box.SetLineColor(0);box.SetFillColor(34);
  TLatex label; label.SetTextSize(0.08); label.SetTextFont(style.nFont); label.SetTextAlign(12);
  TMarker mark; mark.SetMarkerStyle(8); mark.SetMarkerSize(1.1);
  TLine lin; lin.SetLineStyle(1); lin.SetLineWidth(3);
  for(int pad=0; pad<2; pad++){
    cPad = (TPad *)can.cd(pad+1);

    double delta = maxRD[pad] - minRD[pad];
    TString hName = "Ratio"; hName += pad+1;
    hRD[pad] = new TH1F(hName, "", 10, minRD[pad] - 0.1*delta, maxRD[pad] + 0.1*delta);
    hRD[pad]->GetYaxis()->SetNdivisions(0);
    hRD[pad]->SetMaximum(maxHisto);
    hRD[pad]->Draw();
    cout<<RoundNumber(nomiRD[pad][0],3)<<" +- "<<RoundNumber(nomiRD[pad][1],3)<<endl;
    box.SetFillColor(33); box.DrawBox(nomiRD[pad][0]-nomiRD[pad][1], 0, nomiRD[pad][0]+nomiRD[pad][1], maxHisto);
    lin.SetLineColor(33); lin.DrawLine(nomiRD[pad][0], 0, nomiRD[pad][0], maxHisto);
    box.SetFillColor(46);box.DrawBox(SMRD[pad][0]-SMRD[pad][1], 0, SMRD[pad][0]+SMRD[pad][1], maxHisto);
    lin.SetLineColor(46); lin.DrawLine(SMRD[pad][0], 0, SMRD[pad][0], maxHisto);
    for(int file=0; file<nFiles; file++){
      if(pad==0 && file==0) continue;
      int nBeg = 0; if(pad==1) nBeg = 1;
      double halfH = 1/(double)nFiles/2.;
      label.SetTextFont(22); label.SetTextColor(colors[file]); label.SetTextAlign(12);
      label.SetTextSize(0.06); label.DrawLatex(minRD[pad]-0.6*delta, 1-halfH*(2*file+1)+0.075, fileTag[file]);
      int digits = 2; if(file==4) digits = 3;
      TString Error = RoundNumber(RD[file][pad][0],digits); Error += " #pm ";
      Error += RoundNumber(RD[file][pad][1],digits); 
      label.SetTextSize(0.06); label.DrawLatex(minRD[pad]-0.6*delta, 1-halfH*(2*file+1)-0.0, Error);
      if(file==4){
	label.SetTextColor(2); 
	//label.DrawLatex(minRD[pad]+0.6*delta, 1-halfH*(2*file+1)+0.02, "Preliminary");
      }
      for(int cand=nBeg; cand<4; cand+=4){
	double height = 1 - (halfH*(2*file+1) - halfH/3.2);
	//double height = 1 - (halfH*(2*file+1) - halfH/3.2);
	//if(cand>1) height = 1 - (halfH*(2*file+1) + halfH/3.2);
	mark.SetMarkerColor(colors[file]); lin.SetLineColor(colors[file]);
	if(file==4){
	  mark.SetMarkerColor(2); lin.SetLineColor(2);
	}
	mark.DrawMarker(RD[file][cand][0], height);
	lin.DrawLine(RD[file][cand][0]-RD[file][cand][1], height, RD[file][cand][0]+RD[file][cand][1], height);
	chi2[pad] += pow((nomiRD[pad][0]-RD[file][cand][0])/RD[file][cand][1],2); ndof[pad]++;
      }
    }
    TString rightLabel = "#splitline{#chi^{2}: "; rightLabel += RoundNumber(chi2[pad],1); rightLabel += "/";
    rightLabel += ndof[pad]; rightLabel += " = "; rightLabel += RoundNumber(chi2[pad],2,ndof[pad]);
    rightLabel += "}{Prob. = "; rightLabel += RoundNumber(TMath::Prob(chi2[pad],ndof[pad])*100,2); rightLabel += "%}";
    style.setTitles(hRD[pad],xTitle[pad]);//,"","",rightLabel);
    label.SetTextSize(0.07); label.SetTextColor(1);label.SetTextAlign(23);
    label.DrawLatex(SMRD[pad][0], 1.12, "SM");
    label.DrawLatex(nomiRD[pad][0], 1.12, "Aver.");
  }

  TString pName = "public_html/ResultsRD.eps"; 
  can.SaveAs(pName);
  for(int pad=0; pad<2; pad++){
      hRD[pad]->Delete();
  }
}

