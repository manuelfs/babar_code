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
#include "babar_code/Styles/Styles.cc"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void Final_EmuRD(int isIso=0, TString NomFit = "Ne2", TString VarFit = "Ne6"){

  Styles style; style.setPadsStyle(2); style.TextSize = 0.072; style.applyStyle();
  gStyle->SetPadLeftMargin(0.23);

  double SysRD[6][4][2];
  double RD[6][4][2], minRD[2] = {99., 99.}, maxRD[2] = {-99., -99.}, chi2[] = {0,0}; 
  double nomiRD[2][2] = {{0.4441, 0.058}, {0.3331, 0.024}};
  double SMRD[2][2] = {{0.297, 0.017}, {0.252, 0.003}};
  int ndof[2] = {-1,-1}, nCand = 4;
  if(isIso) nCand=2;
  TString folder = "FitAll/fits/TextFinal";
  int nFiles = 2; TString fileTag[] = {"Electron", "Muon"};
  TString nameFiles[] = {"e", "mu"};
  for(int file=0; file<nFiles; file++){
    TString fileName = folder; if(isIso) fileName += "Iso"; fileName += nameFiles[file]; 
    fileName += "Data"; fileName += NomFit; fileName += "x100_RunAll.txt";
    fstream textFile; textFile.open(fileName,fstream::in);
    TString text = "";
    while(!text.Contains("R(D")) textFile>>text;
    for(int cand=0; cand<nCand; cand++){
      if(cand>0) textFile>>text;
      textFile>>text; RD[file][cand][0] = text.Atof();
      textFile>>text; textFile>>text; text.ReplaceAll(",",""); RD[file][cand][1] = text.Atof();
      for(int i=0; i<12; i++) textFile>>text; 
      //for(int i=0; i<4; i++) textFile>>text; 
      if(RD[file][cand][0]-RD[file][cand][1] < minRD[cand%2]) minRD[cand%2] = RD[file][cand][0]-RD[file][cand][1];
      if(RD[file][cand][0]+RD[file][cand][1] > maxRD[cand%2]) maxRD[cand%2] = RD[file][cand][0]+RD[file][cand][1];
      //cout<<RD[file][cand][0]<<" #pm "<<RD[file][cand][1]<<endl;
    }
    //cout<<endl;
  }
  for(int file=0; file<nFiles; file++){
    TString fileName = folder; if(isIso) fileName += "Iso"; fileName += nameFiles[file]; 
    fileName += "Data"; fileName += VarFit; fileName += "x100_RunAll.txt";
    fstream textFile; textFile.open(fileName,fstream::in);
    TString text = "";
    while(!text.Contains("R(D")) textFile>>text;
    for(int cand=0; cand<nCand; cand++){
      if(cand>0) textFile>>text;
      textFile>>text; SysRD[file][cand][0] = text.Atof();
      textFile>>text; textFile>>text; text.ReplaceAll(",",""); SysRD[file][cand][1] = text.Atof();
      for(int i=0; i<12; i++) textFile>>text; 
      //cout<<"VarRD "<<SysRD[file][cand][0]<<", NomRD "<<RD[file][cand][0]<<", errRD "<<RD[file][cand][1];
      SysRD[file][cand][1] = fabs(SysRD[file][cand][0]-RD[file][cand][0]);
      //cout<<", totRD "<<RD[file][cand][1]<<endl;

    }
  }
  maxRD[0] = 0.86; maxRD[1] = 0.61; 
  TCanvas can("can","RD per Run");
  can.Divide(2,1); TPad *cPad;
  TString xTitle[] = {"R(D)", "R(D*)"};
  TH1F* hRD[2];
  int colors[] = {4, 2};
  TBox box; box.SetFillStyle(3002); box.SetLineColor(0);box.SetFillColor(8);
  TLatex label; label.SetTextSize(0.08); label.SetTextFont(style.nFont); label.SetTextAlign(12);
  TMarker mark; mark.SetMarkerStyle(8); mark.SetMarkerSize(0.75);
  TLine lin; lin.SetLineStyle(1); lin.SetLineWidth(2);
  for(int pad=0; pad<2; pad++){
    cPad = (TPad *)can.cd(pad+1);

    double delta = maxRD[pad] - minRD[pad];
    TString hName = "Ratio"; hName += pad+1;
    hRD[pad] = new TH1F(hName, "", 10, minRD[pad] - 0.1*delta, maxRD[pad] + 0.1*delta);
    hRD[pad]->GetYaxis()->SetNdivisions(0);
    hRD[pad]->SetMaximum(1);
    hRD[pad]->Draw();
    box.SetFillColor(8); box.DrawBox(nomiRD[pad][0]-nomiRD[pad][1], 0, nomiRD[pad][0]+nomiRD[pad][1], 1);
    box.SetFillColor(46);box.DrawBox(SMRD[pad][0]-SMRD[pad][1], 0, SMRD[pad][0]+SMRD[pad][1], 1);
    lin.SetLineColor(8);lin.DrawLine(nomiRD[pad][0], 0, nomiRD[pad][0], 1);
    lin.SetLineColor(46);lin.DrawLine(SMRD[pad][0], 0, SMRD[pad][0], 1);
    for(int file=0; file<nFiles; file++){
      int nBeg = 0; if(pad==1) nBeg = 1;
      double halfH = 1/(double)nFiles/2.;
      label.DrawLatex(minRD[pad]-0.45*delta, 1-halfH*(2*file+1), fileTag[file]);
      for(int cand=nBeg; cand<nCand; cand+=2){
	double height = 1 - (halfH*(2*file+1) - halfH/3.2), eStat = RD[file][cand][1], dH = 0.01;
	double eSys = sqrt(pow(SysRD[file][cand][1],2)+pow(eStat,2)), val = RD[file][cand][0];
	if(cand>1) height = 1 - (halfH*(2*file+1) + halfH/3.2);
	mark.SetMarkerColor(colors[cand>1]); lin.SetLineColor(colors[cand>1]);
	mark.DrawMarker(RD[file][cand][0], height);
	lin.DrawLine(val-eSys, height, val+eSys, height);
	lin.DrawLine(val-eSys, height-dH, val-eSys, height+dH); lin.DrawLine(val+eSys, height-dH, val+eSys, height+dH);
	lin.DrawLine(val-eStat, height-dH, val-eStat, height+dH); lin.DrawLine(val+eStat, height-dH, val+eStat, height+dH);
	chi2[pad] += pow((nomiRD[pad][0]-val)/eSys,2); ndof[pad]++;
      }
    }
    double pChi2 = TMath::Prob(chi2[pad],ndof[pad])*100;
    TString rightLabel = "#splitline{#chi^{2}: "; rightLabel += RoundNumber(chi2[pad],1); rightLabel += "/";
    rightLabel += ndof[pad]; rightLabel += " = "; rightLabel += RoundNumber(chi2[pad],2,ndof[pad]);
    rightLabel += "}{Prob. = "; rightLabel += RoundNumber(pChi2,2); rightLabel += "%}";
    style.setTitles(hRD[pad],xTitle[pad],"","",rightLabel);
    cout<<"Chi2 probability for "<<xTitle[pad]<<" is "<<RoundNumber(pChi2,2)<<endl;
    if(pad==1){
      double x = -0.08, y = -0.38;
      box.SetFillColor(46);box.DrawBox(0.55+x, 0.98+y, 0.575+x, 1.03+y);
      label.SetTextSize(0.07); label.DrawLatex(0.59+x, 1.005+y, "SM");
      box.SetFillColor(8);box.DrawBox(0.55+x, 0.89+y, 0.575+x, 0.94+y);
      label.SetTextSize(0.07); label.DrawLatex(0.59+x, 0.915+y, "All data");
      mark.SetMarkerColor(colors[0]); lin.SetLineColor(colors[0]);
      mark.DrawMarker(0.5625+x, 0.81+y);
      lin.DrawLine(0.55+x, 0.81+y, 0.575+x, 0.81+y);
      label.DrawLatex(0.59+x, 0.825+y, "B^{+}");
      mark.SetMarkerColor(colors[1]); lin.SetLineColor(colors[1]);
      mark.DrawMarker(0.5625+x, 0.715+y);
      lin.DrawLine(0.55+x, 0.715+y, 0.575+x, 0.715+y);
      label.DrawLatex(0.59+x, 0.73+y, "B^{0}");
    }
  }

  TString pName = "public_html/Final_EmuRD.eps"; 
  if(isIso) pName.ReplaceAll(".eps","Iso.eps");
  can.SaveAs(pName);
  for(int pad=0; pad<2; pad++){
      hRD[pad]->Delete();
  }
}

