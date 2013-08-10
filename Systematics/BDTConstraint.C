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
#include "TH1F.h"
#include "babar_code/Styles/Styles.cc"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;
#define nFiles 6

void BDTConstraint(){

  Styles style; style.setPadsStyle(2); style.applyStyle();
  //gStyle->SetPadLeftMargin(0.3); gStyle->SetPadTopMargin(0.12);

  int colors[] = {2,4,28,9, 46, 38};
  double RD[nFiles][4][2];
  double minRD[2] = {99., 99.}, maxRD[2] = {-99., -99.}; 

  TString RDTag[] = {"D^{0}", "D*^{0}", "D^{+}", "D*^{+}"};
  double nomiRD[2][2] = {{0.8364, 0.0209}, {0.9613, 0.0423}};
  TString fileTag[] = {"0.1", "0.2", "0.3", "0.4", "0.5", "0.6"};

  TString fileName = "babar_code/Systematics/Text/BDTConstraint.txt", text; 
  fstream textFile; textFile.open(fileName,fstream::in);
  for(int file=0; file<nFiles; file++){
    for(int cand=0; cand<2; cand++){
      textFile>>RD[file][cand][0]>>text>>RD[file][cand][1]; 
      //cout<<fileTag[file]<<"\t"<<RD[file][cand][0]<<" #pm "<<RD[file][cand][1]<<endl;
    }
  }


  double minBDT = 0, maxBDT = 0.7;
  TCanvas can("can","RD results");
  can.Divide(2,1); TPad *cPad;
  TString yTitle[] = {"f_{BDT}^{0}/f_{Nom}^{0}", "f_{BDT}^{+}/f_{Nom}^{+}"};
  TH1F* hRD[2];
  TBox box; box.SetFillStyle(3002); box.SetLineColor(0);box.SetFillColor(34);
  TLatex label; label.SetTextSize(0.08); label.SetTextFont(132); label.SetTextAlign(12);
  TMarker mark; mark.SetMarkerStyle(8); mark.SetMarkerSize(1.05);
  TLine lin; lin.SetLineStyle(1); lin.SetLineWidth(2);
  double minRD2[2] = {0.85, 0.85}, maxRD2[2] = {1.15, 1.15}; 
  for(int pad=0; pad<2; pad++){
    cPad = (TPad *)can.cd(pad+1);

    double delta = maxRD[pad] - minRD[pad];
    TString hName = "Ratio"; hName += pad+1;
    hRD[pad] = new TH1F(hName, "", 10, minBDT, maxBDT);
    hRD[pad]->SetMaximum(maxRD[pad] + 0.1*delta);
    hRD[pad]->SetMinimum(minRD[pad] - 0.1*delta);
    hRD[pad]->SetMaximum(maxRD2[pad]);
    hRD[pad]->SetMinimum(minRD2[pad]);
    hRD[pad]->Draw();
    lin.SetLineColor(colors[4+pad]); lin.DrawLine(minBDT, 1, maxBDT, 1);
    double NomiError = nomiRD[pad][1]/nomiRD[pad][0];
    box.SetFillColor(colors[4+pad]); box.DrawBox(minBDT, 1-NomiError, maxBDT, 1+NomiError);
    for(int file=0; file<nFiles; file++){
      double height = RD[file][pad][0]/nomiRD[pad][0], dH = RD[file][pad][1]/nomiRD[pad][0], BDT = fileTag[file].Atof();
      mark.SetMarkerColor(colors[pad]); lin.SetLineColor(colors[pad]);
      mark.DrawMarker(BDT, height);
      double maxHeight = height+dH; if(maxHeight>maxRD2[pad]) maxHeight=maxRD2[pad];
      double minHeight = height-dH; if(minHeight<minRD2[pad]) minHeight=minRD2[pad];
      lin.DrawLine(BDT, minHeight, BDT, maxHeight);
      if(file==0){
	double delta = 0.1; 
	height = minRD2[pad] + (maxRD2[pad]-minRD2[pad])*(0.93-0.1*(pad>1));
	mark.DrawMarker(25+delta, height);
	label.SetTextSize(0.075); label.SetTextColor(colors[pad]);label.SetTextAlign(12);
	label.DrawLatex(34+delta, height*1.005, RDTag[pad]);	  
      }
      style.setTitles(hRD[pad], "BDT cut", yTitle[pad]);
    }
  }

  TString pName = "public_html/BDTConstraint.eps"; 
  can.SaveAs(pName);
  for(int pad=0; pad<2; pad++){
      hRD[pad]->Delete();
  }
}

