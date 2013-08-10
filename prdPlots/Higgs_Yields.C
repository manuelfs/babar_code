//------------------------------------------------------------------------
// File and Version Information:
//      $Id: RateCalc.hh,v 1.1 2012/03/04 00:33:02 manuelf Exp $
//
// Description:
//      PlotRate - Plots the B->D(*)TauNu rates based on
//                 hep-ph 1203.2654
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      12/03/17 manuelf -- Created 
//------------------------------------------------------------------------

#include "babar_code/NewPhysics/RateCalc.cc"
#include "babar_code/FF/BToDtaunu.cc"
#include "babar_code/FF/BToDstaunu.cc"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include "TMath.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TBox.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include "TLegend.h"
#include "TSpline.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"
#include <fstream>
#include <iostream>

#define nPads 4

using namespace std;
using std::cout;
using std::endl;

void Higgs_Yields(){

  TBox box; box.SetLineColor(10);box.SetFillColor(10);
  int nRows = 2, nCols = 2, nFont = 132, Colors[2][2] = {{kRed+2, kRed-9}, {kAzure-1, kAzure-9}};
  double Range[2][2][2] = {{{97, 111}, {97, 111}}, {{165, 635}, {755, 1055}}};
  double dRows = (double)nRows, dCols = (double) nCols;
  double bMargin = 0.14, padH = (1-bMargin)/dRows, padW = 1/dCols, LeftMargin = 0.23;
  gStyle->SetOptStat(0);              // No Stats box
  gStyle->SetPadTickX(0);             // No ticks at the right
  gStyle->SetPadTickY(0);             // No ticks at the top
  gStyle->SetNdivisions(503, "xy");   // 5 primary ticks and 4 secondary ticks
  gStyle->SetTitleFont(nFont,"xy");          // Set the all 2 axes title font
  gStyle->SetLabelFont(nFont,"xy");          // Set the all 2 axes label font
  gStyle->SetTextFont(nFont);                // Set global text font
  gStyle->SetTitleSize(0.16,"xy");     // Set the 2 axes title size
  gStyle->SetLabelSize(0.16,"xy");     // Set the 2 axes label size
  TCanvas can("can","RD per Run", 800, 450); can.cd();
  TPad *Pads[nPads][2];

  TString yTitle[] = {"E", "Y"}, HigName, hName;
  TSpline3 *MeasuredRD[4], *SplineRD[2], *SplineErrRD[2];
  double RD[2][130], errorRD[2][130], List_tBmH[130], RelerrRD[2][130];
  float weight;
  int candType, MCType;
  double totW[2][6], totW2[2][6], totN[2][2][6], nSM[2]={1,1}, errorSM[2]={1,1};

  int whichB = 2, dHig = 5, file = 0;
  int iniHig = 0, iHig = iniHig, finHig = 100, bin, nFiles = (finHig-iniHig)/dHig+1, nBins = 1000;
  double higRD[4][30];
  TString folder = "FitAll/fits/TextFinalIsoDataHi2x", text = ""; 

  for(int iHig=iniHig; iHig<=finHig; iHig += dHig){
    HigName = "AWG82/ntuples/small/FitRAllHigx"; 
    if(iHig<100) HigName += "0";
    if(iHig<10) HigName += "0";
    HigName += iHig; HigName += "_RunAll.root";
    TChain genChain("ntp1"); genChain.Add(HigName);
    genChain.SetBranchAddress("weight",&weight);
    genChain.SetBranchAddress("candType",&candType);
    genChain.SetBranchAddress("MCType",&MCType);
    List_tBmH[file] = (double)iHig/100.;

    for(int i=0; i<4; i++){
      for(int j=0; j<2; j++){
	totW[j][i] = 0; totW2[j][i] = 0; 
	for(int k=0; k<2; k++) totN[k][j][i] = 0;
      }
    }
    for(int entry=0; entry<genChain.GetEntries(); entry++){
      genChain.GetEvent(entry);
      if(candType<=2 && MCType==5) {
	totW[0][0] += weight; totW2[0][0] += weight*weight;
	if(candType==1) totN[0][0][0] += weight;
	else totN[1][0][0] += weight;
      }
      if(candType<=2 && MCType==6) {
	totW[0][1] += weight; totW2[0][1] += weight*weight;
	if(candType==2) totN[0][0][1] += weight;
	else totN[1][0][1] += weight;
      }
      if(candType>=3 && candType<=4 && MCType==11){
	totW[0][2] += weight; totW2[0][2] += weight*weight;
	if(candType==3) totN[0][0][2] += weight;
	else totN[1][0][2] += weight;
      }
      if(candType>=3 && candType<=4 && MCType==12){
	totW[0][3] += weight; totW2[0][3] += weight*weight;
	if(candType==4) totN[0][0][3] += weight;
	else totN[1][0][3] += weight;
      }

    }
    double DiffError;
    for(int isDs=0; isDs<2; isDs++) {
      if(iHig==0){
	if(whichB==0){
	  nSM[isDs] = totW[0][isDs]; 
	  errorSM[isDs] = totW2[0][isDs]; 
	} else if(whichB==1){
	  nSM[isDs] = totW[0][isDs+2]; 
	  errorSM[isDs] = totW2[0][isDs+2]; 
	} else {
	  nSM[isDs] = totW[0][isDs]+totW[0][isDs+2];
	  errorSM[isDs] = totW2[0][isDs]+totW2[0][isDs+2]; 
	}
      }
      if(whichB==0){
	RD[isDs][file] = totW[0][isDs]/nSM[isDs]*100;
	DiffError = totW2[0][isDs]/pow(totW[0][isDs],2)-errorSM[isDs]/pow(nSM[isDs],2);
	RelerrRD[isDs][file] = sqrt(totW2[0][isDs]/pow(totW[0][isDs],2));
      } else if(whichB==1){
	RD[isDs][file] = totW[0][isDs+2]/nSM[isDs]*100;
	DiffError = totW2[0][isDs+2]/pow(totW[0][isDs+2],2)-errorSM[isDs]/pow(nSM[isDs],2);
	RelerrRD[isDs][file] = sqrt(totW2[0][isDs+2]/pow(totW[0][isDs+2],2));
      } else {
	RD[isDs][file] = (totW[0][isDs]+totW[0][isDs+2])/nSM[isDs]*100;
	DiffError = (totW2[0][isDs]+totW2[0][isDs+2])/(pow(totW[0][isDs]+totW[0][isDs+2],2))
	  -errorSM[isDs]/pow(nSM[isDs],2);
	RelerrRD[isDs][file] = sqrt((totW2[0][isDs]+totW2[0][isDs+2])/
				    (pow(totW[0][isDs]+totW[0][isDs+2],2)));
      }
      if(DiffError<0) DiffError=0;
      errorRD[isDs][file] = sqrt(DiffError)*RD[isDs][file];
      //cout<<RoundNumber(List_tBmH[file],2)<<": Y/Y_SM = "<<RoundNumber(RD[isDs][file],2)
      //<<" +- "<<RoundNumber(errorRD[isDs][file],2)<<" \t";
    }
    //cout<<endl;
    file++;
  }
  for(int isDs=0; isDs<2; isDs++) {
    can.cd(isDs+1);
    HigName = "Spline"; HigName += isDs; 
    SplineRD[isDs] = new TSpline3(HigName,List_tBmH,RD[isDs],nFiles,"ble1",0); 
    HigName += "Error";
    SplineErrRD[isDs] = new TSpline3(HigName,List_tBmH,errorRD[isDs],nFiles,"ble1",0);
  }


  for(file=0; file<nFiles; file++){
    TString fileName = folder; 
    if(iHig<100) fileName += "0";
    if(iHig<10) fileName += "0";
    fileName += iHig; fileName += ".txt";
    fstream textFile; textFile.open(fileName,fstream::in);
    bin = 0;
    while(!text.Contains("R(D")) {
      textFile>>text; 
      bin++;
      if(bin>1000) {cout<<"R(D) not found in "<<fileName<<endl; return;}
    }
    for(int cand=0; cand<2; cand++){
      if(cand>0) textFile>>text;
      textFile>>text; higRD[cand][file] = text.Atof();
      textFile>>text; textFile>>text; text.ReplaceAll(",",""); higRD[cand+2][file] = text.Atof();
      for(int i=0; i<5; i++) textFile>>text; 
      textFile>>text; higRD[cand][file] = text.Atof();
      textFile>>text; textFile>>text; text.ReplaceAll(",",""); higRD[cand+2][file] = text.Atof();
      for(int i=0; i<4; i++) textFile>>text; 
      //cout<<higRD[cand][file]<<" #pm "<<higRD[cand+2][file]<<" - a "
      //<<RoundNumber(higRD[cand+2][file]*100,1,higRD[cand][file])<<"% error"<<endl;
      List_tBmH[file] = (double)iHig/100.;
    }
    //cout<<endl;
    iHig += dHig;
  }
  for(int isDs=0; isDs<4; isDs++){
    text = "SplineRD"; text += isDs;
    MeasuredRD[isDs] = new TSpline3(text,List_tBmH,higRD[isDs],nFiles,"ble1",0); 
  }

  TH1F *Histo[2][2][2];
  for(int his=0; his<2; his++) {
    for(int col=0; col<nCols; col++){
      for(int isDs=0; isDs<2; isDs++) {
	HigName = "Histo"; HigName += col; HigName += isDs; HigName += his;
	Histo[col][isDs][his] = new TH1F(HigName,"",nBins,0,1);
	Histo[col][isDs][his]->SetLineWidth(3); 
      }
    }
  }

  for(int col=0; col<nCols; col++){
    for(int row=0; row<nRows; row++){
      can.cd(0);
      int pad = nRows*col+row;
      double TextSize = 0.11;
      double RightMargin = 0.03, TopMargin=0, BottomMargin=bMargin/(bMargin+padH)*row;
      double PadXY[2][2] = {{padW*col, padW*(col+1)},{(padH+bMargin)*(dRows-1-row), bMargin+padH*(2-row)}};
      if(row==0) TextSize *= (bMargin+padH)/padH;
      hName = "Pad0_"; hName += pad;
      //cout<<hName<<": "<<PadXY[0][0]<<", "<<PadXY[1][0]<<", "<<PadXY[0][1]<<", "<<PadXY[1][1]<<endl;
      Pads[pad][0] = new TPad(hName,"",PadXY[0][0], PadXY[1][0], PadXY[0][1], PadXY[1][1]);
      Pads[pad][0]->SetLeftMargin(LeftMargin);     Pads[pad][0]->SetRightMargin(RightMargin); 
      Pads[pad][0]->SetBottomMargin(BottomMargin); Pads[pad][0]->SetTopMargin(TopMargin); 
      Pads[pad][0]->Draw(); Pads[pad][0]->cd();

      for(int bin=1; bin<=nBins; bin++){
	double tBmH = Histo[col][row][0]->GetBinCenter(bin), valEffi = MeasuredRD[row]->Eval(tBmH);
	double errorEffi = MeasuredRD[row+2]->Eval(tBmH);
	if(col==0){
	  valEffi = SplineRD[row]->Eval(tBmH);
	  errorEffi = SplineErrRD[row]->Eval(tBmH);
	}
	for(int his=0; his<2; his++) {
	  Histo[col][row][his]->SetBinContent(bin, valEffi);
	  if(his==0) Histo[col][row][his]->SetBinError(bin, errorEffi);
	  else       Histo[col][row][his]->SetBinError(bin,0);
	}
      }
      Histo[col][row][0]->SetTitleOffset(1.1,"x");         
      Histo[col][row][0]->SetTitleSize(TextSize,"xy");      // Set the 2 axes title size
      Histo[col][row][0]->SetLabelSize(TextSize,"xy");      // Set the 2 axes label size
      Histo[col][row][0]->SetTitleFont(nFont,"xy");         // Set the all 2 axes title font
      Histo[col][row][0]->SetLabelFont(nFont,"xy");         // Set the all 2 axes label font
      Histo[col][row][0]->SetNdivisions(804, "xy");         // 5 primary ticks and 4 secondary ticks
      Histo[col][row][0]->SetLabelOffset(0.01,"Y");
      Histo[col][row][0]->GetXaxis()->CenterTitle(true);
      if(row==1) Histo[col][row][0]->SetXTitle("B");
      Histo[col][row][0]->SetLineColor(Colors[col][1]); Histo[col][row][0]->SetFillColor(Colors[col][1]);
      Histo[col][row][0]->SetMinimum(Range[col][row][0]); Histo[col][row][0]->SetMaximum(Range[col][row][1]); 
      Histo[col][row][0]->Draw("e3");


      Histo[col][row][1]->SetLineColor(Colors[col][0]); 
      Histo[col][row][1]->Draw("c same");    
    }
    can.cd(0);
    TLatex label; label.SetTextFont(nFont); label.SetNDC(kTRUE);
    label.SetTextAlign(23); label.SetTextAngle(90); label.SetTextSize(0.07);
    label.DrawLatex(0.01+col/dCols, 0.58, yTitle[col]);
  }

  TString epsName = "public_html/Higgs_Yields.eps";
  can.SaveAs(epsName);

  for(int row=0; row<nRows; row++){
    if(SplineRD[row]) SplineRD[row]->Delete();
    if(SplineErrRD[row]) SplineErrRD[row]->Delete();
    if(MeasuredRD[row]) MeasuredRD[row]->Delete();
    if(MeasuredRD[row+2]) MeasuredRD[row+2]->Delete();
    for(int col=0; col<nCols; col++){
      for(int his=0; his<2; his++) 
	if(Histo[col][row][his]) Histo[col][row][his]->Delete();
    }
  }
}




