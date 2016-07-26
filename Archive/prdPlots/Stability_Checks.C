#include "TString.h"
#include "TMath.h"
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
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void Stability_Checks(){

  double SysRD[6][8][2];
  double RD[6][8][2], minRD[2] = {0.05, 0.17}, maxRD[2] = {1.1, 0.67}, chi2[] = {0,0,0,0}; 
  double nomiRD[2][2] = {{0.4441, 0.0576}, {0.3332, 0.0237}};
  double SMRD[2][2] = {{0.297, 0.017}, {0.252, 0.003}};
  int ndof[] = {-1,-1, -1,-1};
  TString folder = "FitAll/fits/TextFinalDataNe2x100_Run";
  int nFiles = 4; TString fileTag[] = {"Runs 123", "Run 4", "Run 5", "Run 6", "Electron", "Muon"};
  TString nameFiles[] = {"123", "4", "5", "6","e", "mu"};
  for(int file=0; file<nFiles; file++){
    TString fileName = folder; fileName += nameFiles[file]; fileName += ".txt";
    fstream textFile; textFile.open(fileName,fstream::in);
    TString text = "";
    while(!text.Contains("R(D")) textFile>>text;
    for(int cand=0; cand<4; cand++){
      if(cand>0) textFile>>text;
      textFile>>text; RD[file][cand][0] = text.Atof();
      textFile>>text; textFile>>text; text.ReplaceAll(",",""); RD[file][cand][1] = text.Atof();
      for(int i=0; i<12; i++) textFile>>text; 
      //cout<<RD[file][cand][0]<<" #pm "<<RD[file][cand][1]<<endl;
      SysRD[file][cand][1] = 0;
    }
    //cout<<endl;
  }
  for(int file=0; file<2; file++){
    TString fileName = "FitAll/fits/TextFinal"; fileName += nameFiles[file+4]; 
    fileName += "DataNe2x100_RunAll.txt";
    fstream textFile; textFile.open(fileName,fstream::in);
    TString text = "";
    while(!text.Contains("R(D")) textFile>>text;
    for(int cand=4; cand<8; cand++){
      if(cand>4) textFile>>text;
      textFile>>text; RD[file][cand][0] = text.Atof();
      textFile>>text; textFile>>text; text.ReplaceAll(",",""); RD[file][cand][1] = text.Atof();
      for(int i=0; i<12; i++) textFile>>text; 
      //cout<<RD[file][cand][0]<<" #pm "<<RD[file][cand][1]<<endl;
    }
    //cout<<endl;
  }
  for(int file=0; file<2; file++){
    TString fileName = "FitAll/fits/TextFinal"; fileName += nameFiles[file+4]; 
    fileName += "DataNe6x100_RunAll.txt";
    //cout<<"Opening "<<fileName<<endl;
    fstream textFile; textFile.open(fileName,fstream::in);
    TString text = "";
    while(!text.Contains("R(D")) textFile>>text;
    for(int cand=4; cand<8; cand++){
      if(cand>4) textFile>>text;
      textFile>>text; SysRD[file][cand][0] = text.Atof();
      textFile>>text; textFile>>text; text.ReplaceAll(",",""); SysRD[file][cand][1] = text.Atof();
      for(int i=0; i<12; i++) textFile>>text; 
      //cout<<"VarRD "<<SysRD[file][cand][0]<<", NomRD "<<RD[file][cand][0]<<", errRD "<<RD[file][cand][1];
      SysRD[file][cand][1] = fabs(SysRD[file][cand][0]-RD[file][cand][0]);
      //cout<<", totRD "<<RD[file][cand][1]<<endl;

    }
  }

  gStyle->SetOptStat(0);              // No Stats box
  gStyle->SetPadTickX(0);             // No ticks at the right
  gStyle->SetPadTickY(0);             // No ticks at the top
  gStyle->SetNdivisions(503, "xy");   // 5 primary ticks and 4 secondary ticks
  gStyle->SetTitleFont(132,"xy");          // Set the all 2 axes title font
  gStyle->SetLabelFont(132,"xy");          // Set the all 2 axes label font
  gStyle->SetTextFont(132);                // Set global text font
  gStyle->SetTitleSize(0.16,"xy");     // Set the 2 axes title size
  gStyle->SetLabelSize(0.16,"xy");     // Set the 2 axes label size
  double LMargin = 0.17, TopF = 0.42, BMargin = 0.3;
  TCanvas can("can","RD per Run", 800, 450); can.cd();
  TPad *cPad[4];
  cPad[0] = new TPad("pad0","", LMargin, TopF, LMargin + (1-LMargin)/2., 1);
  cPad[1] = new TPad("pad1","", LMargin + (1-LMargin)/2., TopF, 1, 1);
  cPad[2] = new TPad("pad2","", LMargin, 0, LMargin + (1-LMargin)/2., TopF);
  cPad[3] = new TPad("pad3","", LMargin + (1-LMargin)/2., 0, 1, TopF);
  for(int pad=0; pad<4; pad++){
    cPad[pad]->SetTopMargin(0); cPad[pad]->SetRightMargin(0);  cPad[pad]->SetLeftMargin(0); 
    if(pad<2) cPad[pad]->SetBottomMargin(0); 
    else      cPad[pad]->SetBottomMargin(BMargin); 
  }
  TString xTitle[] = {"R(D)", "R(D*)"};
  TH1F* hRD[4];
  int colors[] = {4, 2, kGreen+2, kSpring-4, kRed+2, kRed-9}, markStyle[] = {8, 21};
  TBox box; box.SetLineColor(0);box.SetFillColor(8);
  TLatex label; label.SetTextSize(0.08); label.SetTextFont(132); label.SetTextAlign(12); label.SetNDC(kTRUE);
  TMarker mark; mark.SetMarkerStyle(8); mark.SetMarkerSize(1.1);
  TLine lin; lin.SetLineStyle(1); lin.SetLineWidth(2);
  //cout<<"Starting pad for"<<endl;
  for(int ipad=0; ipad<4; ipad++){
    can.cd();
    cPad[ipad]->Draw(); cPad[ipad]->cd();
    int pad = ipad%2;

    TString hName = "Ratio"; hName += ipad+1;
    hRD[ipad] = new TH1F(hName, "", 10, minRD[pad], maxRD[pad]);
    hRD[ipad]->GetYaxis()->SetNdivisions(0);
    hRD[ipad]->SetMaximum(1);
    if(ipad>1) {
      hRD[ipad]->SetXTitle(xTitle[pad]); 
      hRD[ipad]->GetXaxis()->CenterTitle(true);
      nFiles = 2;
    } else nFiles = 4;
    hRD[ipad]->Draw();
    box.SetFillColor(colors[3]); box.DrawBox(nomiRD[pad][0]-nomiRD[pad][1], 0.003, nomiRD[pad][0]+nomiRD[pad][1], 0.997);
    box.SetFillColor(colors[5]);box.DrawBox(SMRD[pad][0]-SMRD[pad][1], 0.003, SMRD[pad][0]+SMRD[pad][1], 0.997);
    lin.SetLineColor(colors[2]);lin.DrawLine(nomiRD[pad][0], 0.003, nomiRD[pad][0], 0.997);
    lin.SetLineColor(colors[4]);lin.DrawLine(SMRD[pad][0], 0.003, SMRD[pad][0], 0.997);
    for(int file=0; file<nFiles; file++){
      int nBeg = 0; if(pad==1) nBeg = 1;
      double halfH = 1/(double)nFiles/2.;
      for(int cand=nBeg; cand<4; cand+=2){
	double height = 1 - (halfH*(2*file+1) - halfH/3.2);
	if(cand>1) height = 1 - (halfH*(2*file+1) + halfH/3.2);
	mark.SetMarkerStyle(markStyle[cand>1]);
	mark.SetMarkerColor(colors[cand>1]); lin.SetLineColor(colors[cand>1]);
	double val = RD[file][cand+4*(ipad>1)][0], dH = 0.015;
	double eStat = RD[file][cand+4*(ipad>1)][1], eSys = sqrt(pow(SysRD[file][cand+4*(ipad>1)][1],2)+pow(eStat,2));
	if(ipad>1) dH *= 2;
	mark.DrawMarker(val, height);
	lin.DrawLine(val-eSys, height, val+eSys, height);
	lin.DrawLine(val-eSys, height-dH, val-eSys, height+dH); lin.DrawLine(val+eSys, height-dH, val+eSys, height+dH);
	if(ipad>1)
	  lin.DrawLine(val-eStat, height-dH, val-eStat, height+dH); lin.DrawLine(val+eStat, height-dH, val+eStat, height+dH);
	chi2[ipad] += pow((nomiRD[pad][0]-val)/eSys,2); ndof[ipad]++;
      }
    }
    TString rightLabel = "#splitline{#chi^{2}: "; rightLabel += RoundNumber(chi2[ipad],1); rightLabel += "/";
    rightLabel += ndof[ipad]; rightLabel += "}{p = "; rightLabel += RoundNumber(TMath::Prob(chi2[ipad],ndof[ipad])*100,1); 
    rightLabel += "%}";
    label.SetTextAlign(33); label.SetTextSize(0.105); 
    if(ipad>1) label.SetTextSize(0.14); 
    label.DrawLatex(0.98, 0.95, rightLabel);
    if(ipad==1){
      double x = -0.08, y = -0.5;
      label.SetTextAlign(12); label.SetNDC(kFALSE);

      box.SetFillColor(colors[5]);box.DrawBox(0.55+x, 1.03+y, 0.575+x, 1.08+y);
      label.SetTextSize(0.07); label.DrawLatex(0.59+x, 1.045+y, "SM");
      box.SetFillColor(colors[3]);box.DrawBox(0.55+x, 0.93+y, 0.575+x, 0.98+y);
      label.SetTextSize(0.07); label.DrawLatex(0.59+x, 0.94+y, "All data");

      mark.SetMarkerColor(colors[0]); mark.SetMarkerStyle(markStyle[0]);
      mark.DrawMarker(0.5625+x, 0.83+y);   label.DrawLatex(0.59+x, 0.83+y, "Bm");
      mark.SetMarkerColor(colors[1]); mark.SetMarkerStyle(markStyle[1]);
      mark.DrawMarker(0.5625+x, 0.715+y);  label.DrawLatex(0.59+x, 0.715+y, "B0");
      label.SetNDC(kTRUE);
    }
    hRD[ipad]->Draw("axis same");

  }

  can.cd(); label.SetTextAlign(12); label.SetTextSize(0.07); 
  double dH = (1-TopF)/4., BotF = TopF*BMargin;
  for(int file=0; file<4; file++)
      label.DrawLatex(0.015, 1-dH*(file+0.5), fileTag[file]);
  dH = (TopF-BotF)/2.;
  for(int file=0; file<4; file++)
      label.DrawLatex(0.015, TopF-dH*(file+0.5), fileTag[file+4]);


  TString pName = "public_html/Stability_Checks.eps"; 
  can.SaveAs(pName);
  for(int pad=0; pad<4; pad++){
      hRD[pad]->Delete();
  }
}

