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
#define nFiles 3
#define nBDT 9

void Final_Chi2(int Variable = 2, int Channel = 4){
  //Variable: 0 mmiss, 1 pstarl, 2 combined
  //Channel: 0 B+, 2 B-, 4 Iso

  Styles style; style.setPadsStyle(2); style.applyStyle();
  //gStyle->SetPadLeftMargin(0.3); gStyle->SetPadTopMargin(0.12);

  int colors[] = {2,4,28,9, 46, 38};
  double pChi2[nFiles][nBDT][6][3],ExtremChi2[nBDT][6][3][2];
  double minpChi2[2] = {99., 99.}, maxpChi2[2] = {-99., -99.}; 
  for(int iBDT=0; iBDT<nBDT; iBDT++){
    for(int cand=0; cand<6; cand++){
      for(int type=0; type<3; type++){
	ExtremChi2[iBDT][cand][type][0] = 100.;
	ExtremChi2[iBDT][cand][type][1] = 0.;
      }
    }
  }
  TString folder = "babar_code/Systematics/Text/Chi2_"; 
  TString nameFiles[] = {"05_05", "10_10", "15_15"};
  TString nameBDT[] = {"030", "050", "075", "100", "120", "150", "200", "250", "300"};

  TString pChi2Tag[] = {"D^{0}", "D*^{0}", "D^{+}", "D*^{+}"};

  for(int file=0; file<nFiles; file++){
    TString fileName = folder; fileName += nameFiles[file]; fileName += ".txt";
    fstream textFile; textFile.open(fileName,fstream::in);
    TString text = ""; 
    for(int iBDT=0; iBDT<nBDT; iBDT++){
      for(int cand=0; cand<4; cand++){
	while(text != "M2") textFile>>text; 
	textFile>>text; pChi2[file][iBDT][cand][0] = text.Atof();
	textFile>>text; textFile>>text; pChi2[file][iBDT][cand][1] = text.Atof();
	textFile>>text; textFile>>text; pChi2[file][iBDT][cand][2] = text.Atof();
	if(cand>1){
	  textFile>>text; textFile>>text; pChi2[file][iBDT][cand+2][0] = text.Atof();
	  textFile>>text; textFile>>text; pChi2[file][iBDT][cand+2][1] = text.Atof();
	  textFile>>text; textFile>>text; pChi2[file][iBDT][cand+2][2] = text.Atof();
	}
	//cout<<pChi2[file][iBDT][cand][0]<<" and total "<<pChi2[file][iBDT][cand][2]<<endl;
      }
	for(int i=0; i<9; i++) textFile>>text; 
	//cout<<endl;
    }
  }
  double minBDT = 0, maxBDT = 330;
  //gPad->SetLogy(1);
  TCanvas can("can","pChi2 results");
  can.Divide(2,1); TPad *cPad;
  TString yTitle[] = {"p(#chi^{2}) (%)", "p(#chi^{2}) (%)"};
  TH1F* hpChi2[2];
  TBox box; box.SetFillStyle(3002); box.SetLineColor(0);box.SetFillColor(34);
  TLatex label; label.SetTextSize(0.08); label.SetTextFont(132); label.SetTextAlign(12);
  TMarker mark; mark.SetMarkerStyle(8); mark.SetMarkerSize(0.92);
  TLine lin; lin.SetLineStyle(1); lin.SetLineWidth(2);
  double minpChi22[2] = {0.00001, 0.001}, maxpChi22[2] = {200, 200}; 
  for(int pad=0; pad<2; pad++){
    cPad = (TPad *)can.cd(pad+1);
    cPad->SetLogy(1);
    double delta = maxpChi2[pad] - minpChi2[pad];
    TString hName = "Ratio"; hName += pad+1;
    hpChi2[pad] = new TH1F(hName, "", 10, minBDT, maxBDT);
    //hpChi2[pad]->GetYaxis()->SetNdivisions(0);
    hpChi2[pad]->SetMaximum(maxpChi2[pad] + 0.1*delta);
    hpChi2[pad]->SetMinimum(minpChi2[pad] - 0.1*delta);
    hpChi2[pad]->SetMaximum(maxpChi22[pad]);
    hpChi2[pad]->SetMinimum(minpChi22[pad]);
    hpChi2[pad]->Draw();

    double NomiMin=100, NomiMax=0, val;
    for(int file=0; file<3; file++){
      val = pChi2[file][3][pad+Channel][2]/100;
      if(NomiMin > val) NomiMin = val;
      if(NomiMax < val) NomiMax = val;
    }
    box.SetFillColor(colors[pad]); box.DrawBox(minBDT, NomiMin, maxBDT, NomiMax);

    for(int iBDT=0; iBDT<nBDT; iBDT++){
      if(nameBDT[iBDT]=="080") continue;
      int nBeg = 0; if(pad==1) nBeg = 1;
      for(int cand=nBeg; cand<2; cand+=2){
	double BDT = nameBDT[iBDT].Atof();
	for(int file=0; file<3; file++){
	  double height = pChi2[file][iBDT][cand+Channel][2]/100.;	
	  //if(cand>1) BDT += 8;
	  if(height<ExtremChi2[iBDT][cand][Variable][0]) ExtremChi2[iBDT][cand][Variable][0] = height;
	  if(height>ExtremChi2[iBDT][cand][Variable][1]) ExtremChi2[iBDT][cand][Variable][1] = height;
	  mark.SetMarkerColor(colors[cand]); lin.SetLineColor(colors[cand]);
	  if(file>0 || 1){
	    mark.DrawMarker(BDT, height);
	  }
	  if(file==0 && 0){
	    double delta = 25;
	    height = minpChi22[pad] + (maxpChi22[pad]-minpChi22[pad])*(0.93-0.1*(cand>1));
	    mark.DrawMarker(25+delta, height);
	    label.SetTextSize(0.075); label.SetTextColor(colors[cand]);label.SetTextAlign(12);
	    label.DrawLatex(34+delta, height*1.005, pChi2Tag[cand]);
	  }

	}
	lin.DrawLine(BDT, ExtremChi2[iBDT][cand][Variable][0], BDT, ExtremChi2[iBDT][cand][Variable][1]);
	style.setTitles(hpChi2[pad], "% of the nominal sample", yTitle[pad]);
      }
    }
  }

  TString pName = "public_html/Final_Chi2.eps"; 
  can.SaveAs(pName);
  for(int pad=0; pad<2; pad++){
      hpChi2[pad]->Delete();
  }
}

