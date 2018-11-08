//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: PRL_Plot.cc,v 1.2 2012/08/23 02:22:18 manuelf Exp $
//
// Description:
//      PlotFit - Plots the results of the final fit for the PRL
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      12/02/17 manuelf -- Adapted from PlotFit.cc
//------------------------------------------------------------------------
//////////////////////////////
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TText.h"
#include <fstream>
#include <iostream>

/////////////////////////////
#include "styles.hpp"
#include "keys_utils.hpp"

namespace{
  bool texLabels=true;
};

using namespace TMath;
using namespace std;

void BABARLabel(Double_t xpos= 0.73, Double_t ypos= 0.85, Double_t scale= 1.0) {
  TText *babar = new TText();
  babar->SetNDC(kTRUE);
  babar->SetTextFont(32);
  babar->SetTextSize(0.10*scale);
  babar->DrawText(xpos,ypos,"B");
  babar->SetTextSize(0.075*scale);
  babar->DrawText(xpos+0.042*scale,ypos,"A");
  babar->SetTextSize(0.10*scale);
  babar->DrawText(xpos+0.078*scale,ypos,"B");
  babar->SetTextSize(0.075*scale);
  babar->DrawText(xpos+0.120*scale,ypos,"AR");
  delete babar;
}

void PlotFinalFit2(TString textName, TTree *tree, int nbins = 25, double minx=-0.4, double maxx=1.4, double minY=0, double maxY=2.4,
		  TString typePlot="histo", TString var="candM2", TString sample="sig", 
		  double maxFactor = 1., TString psEPS = "EPS", double maxFactorPl = 1., int nbinsPl = 25);

int main(int argc, char *argv[]){
  if (argc < 2 || argc > 14 ) {
    cout << "USAGE: ./run/prd_fit.exe textFile [sample=sig] [Subtract=] [var=candM2] [nbins=25] "<< 
      "[minX=-0.5] [maxX=1.55] [typePlot=curve] [maxFactor=1] [minY=0] [maxY=2.4] [maxFactorPl=1] [nbinsPl]" << endl;
    return 0;
  }
  // Setting the input parameters
  TString textName = argv[1];
  TString sample = "sig";
  if (argc>2) sample = argv[2];
  TString Subtract = "";
  if (argc>3) Subtract = argv[3];
  TString var = "candM2";
  if (argc>4) var = argv[4];
  int nbins = 25;
  if(argc>5) {TString temp_s = argv[5]; nbins = temp_s.Atoi();} 
  double min = -0.5;
  if(argc>6) {TString temp_s = argv[6]; min = temp_s.Atof();} 
  double max = 1.55;
  if(argc>7) {TString temp_s = argv[7]; max = temp_s.Atof();} 
  TString typePlot = "curve";
  if(argc>8) typePlot = argv[8];
  double maxFactor = 1;
  if(argc>9) {TString temp_s = argv[9]; maxFactor = temp_s.Atof();} 
  double minY = 0, maxY = 2.4;
  if(var=="candPstarLep"){minY = -4; maxY = 12;}
  if(argc>10) {TString temp_s = argv[10]; minY = temp_s.Atof();} 
  if(argc>11) {TString temp_s = argv[11]; maxY = temp_s.Atof();} 
  double maxFactorPl = 1;
  if(argc>12) {TString temp_s = argv[12]; maxFactorPl = temp_s.Atof();}  
  int nbinsPl = nbins;
  if(argc>13) {TString temp_s = argv[13]; nbinsPl = temp_s.Atoi();} 

  TString TagSample = textName;
  TagSample.Remove(0,TagSample.Index("Final")+5); TagSample.ReplaceAll(".txt",""); TagSample.ReplaceAll("IsoData","Data");
  TagSample.ReplaceAll("Ne2","New"); 
  TagSample.ReplaceAll("mES","Sig"); TagSample.ReplaceAll("on2","onD"); TagSample.ReplaceAll("o22","onD");
  TagSample.ReplaceAll("Si2","Ori"); TagSample.ReplaceAll("Sid","Ori"); TagSample.ReplaceAll("Or2","Ori"); 
  TagSample.ReplaceAll("Df2","Dfi"); TagSample.ReplaceAll("Ds2","Dsf"); 

  TString ntufolder = "/cms2r0/slac/small/";
  TString genFile = ntufolder+"Fit"; genFile += TagSample; genFile += "_RunAll.root";
  if(textName.Contains("Data")) {
    //if(!TagSample.Contains("Sigx")) genFile = ntufolder+"FitData_RunAll.root";
    if(textName.Contains("1234")) genFile = ntufolder+"FitData_Run1234.root";
    if(textName.Contains("56")) genFile = ntufolder+"FitData_Run56.root";
    if(textName.Contains("eData")) genFile = ntufolder+"FiteDataSidx100_RunAll.root";
    if(textName.Contains("muData")) genFile = ntufolder+"FitmuDataSidx100_RunAll.root";
  }
  if(textName.Contains("Dpipi")) {
    genFile = ntufolder+"FitDataSidx100_RunAll.root";
    //typePlot += "Dpipi";
  }
  TChain *gen = new TChain("ntp1");
  gen->Add(genFile);
  PlotFinalFit2(textName, gen, nbins, min, max, minY, maxY, typePlot, var, sample, maxFactor, Subtract, maxFactorPl, nbinsPl);

  delete gen;
  return 1;
}

void PlotFinalFit2(TString textName, TTree *tree, int nbins, double minx, double maxx, double minY, double maxY,
		  TString typePlot, TString var, TString sample, 
		  double maxFactor, TString Subtract, double maxFactorPl, int nbinsPl){

  if(maxFactorPl!=1){cout<<"Changed maxFactorPl to "<<maxFactorPl<<endl;} // Had to use it to compile
  int truInd=0, isDss = 0, doSubtract = 0, isIso = 0, isSum = 0, nPads = 4;
  if(sample.Contains("dss")) isDss=1;
  if(typePlot.Contains("Tru")) truInd=8;
  if(typePlot.Contains("Sum")) {isSum = 1; nPads  = 1;}
  if(typePlot.Contains("Iso")) {isIso = 1; nPads /= 2;}

  double Yields[16][9], hIntegral[8][2], iYields[70][2];
  for(int pad=0; pad<8; pad++) for(int mc=0; mc<2; mc++) hIntegral[pad][mc] = 0;
  TString buffer, Yname;
  fstream textFile; textFile.open(textName,fstream::in);
  int begChan = 0, endChan = 4, IndText[8] = {0,4,1,5,2,6,3,7};
  for(int i=2*begChan; i<2*endChan; i++){
    int index = IndText[i];
    for(int j=0; j<9; j++){
      textFile>>Yname>>Yields[index+8][j]>>Yields[index][j]>>buffer>>buffer>>buffer;
      //cout<<"True: "<<Yields[index+8][j]<<"\tFit: "<<Yields[index][j]<<"\tChannel "<<index<<endl;
      Yname.ReplaceAll("Yield[",""); Yname.ReplaceAll("]","");
      iYields[Yname.Atoi()][1] = Yields[index+8][j];
      iYields[Yname.Atoi()][0] = Yields[index][j];
      if(index>3 && j==6){
	textFile>>buffer;
	break;
      }
    }
    //cout<<endl;
  }
  int nm2bin = 0, nplbin = 0, Projbins = nbins, ProjbinsPl = nbinsPl;
  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4, insetPlmax=1.99;
  fstream ResWFile; ResWFile.open("txt/resolution_widths.txt",fstream::in);
  double ResWidth[4];
  for(int chan=0; chan<4; chan++) ResWFile >> ResWidth[chan];

  TString TagSample = textName;
  TagSample.Remove(0,TagSample.Index("Final")+5); TagSample.ReplaceAll(".txt","");
  if(TagSample.Contains("x")) {
    if(TagSample.Contains("Data")) TagSample.Remove(0,TagSample.Index("Data")+4); 
    if(TagSample.Contains("RAll")) TagSample.Remove(0,TagSample.Index("RAll")+4); 
    TagSample.Remove(TagSample.First('x')+4,TagSample.Sizeof()+1);
    TagSample.ReplaceAll("Ne2","New"); 
    TagSample.ReplaceAll("mES","Sig"); TagSample.ReplaceAll("on2","onD"); TagSample.ReplaceAll("o22","onD");
    TagSample.ReplaceAll("Si2","Ori"); TagSample.ReplaceAll("Sid","Ori"); TagSample.ReplaceAll("Or2","Ori"); 
    TagSample.ReplaceAll("Df2","Dfi"); TagSample.ReplaceAll("Ds2","Dsf"); TagSample.ReplaceAll("IsoData","Data");
  }
  if(typePlot.Contains("Dpipi")) TagSample = "Sigx100";
  if(textName.Contains("eData")) TagSample = "eSigx100";
  if(textName.Contains("muData")) TagSample = "muSigx100";
  TH2F *pdf[70];
  for (int i = 1 ; i <= 68 ; i ++) {
    //if(i==37)i=41; 
    if(i==49)i=51; if(i==59)i=61; 
    TString hname = "/cms2r0/slac/keys/root/Fit/pdfKeys_"; 
    if(TagSample.Contains("x")) {hname = "/cms2r0/slac/keys/root/fit"; hname += TagSample; hname += "/pdfKeys_"; }
    hname += i; hname += "_Fit.root";
    TFile hfile(hname); 
    TString pdfName = "pdf"; pdfName += i;
    pdf[i] = static_cast<TH2F *>((hfile.Get("h2"))->Clone(pdfName));
    if(typePlot.Contains("Conv") && i>=21&&i<=24) ConvKeys(pdf[i],ResWidth[i-21]);
    pdf[i]->SetDirectory(0);
    if(i==1) {nm2bin = pdf[1]->GetNbinsX(); nplbin = pdf[1]->GetNbinsY();}
    //if(nm2bin!=pdf[i]->GetNbinsX() || nplbin!=pdf[i]->GetNbinsY()) return can;
  }
  int nxbins = nm2bin; double ProjMin=minx, ProjMax=maxx; 
  TString otherVar = "candPstarLep";
  if(var=="candPstarLep") {
    nxbins = nplbin;
    otherVar = "candM2";
  }
  if(typePlot.Contains("curve")) {
    Projbins = nxbins; ProjbinsPl = nplbin;
    if(var=="candPstarLep") {ProjMin=0; ProjMax=2.4;}
    else {ProjMin=-4; ProjMax=12;}
  }
  int Indpdf[8][9] = {{61,51,41,37,13,17,9,5,1},   {62,52,42,38,10,18,14,2,6},
		      {63,53,43,39,15,19,11,7,3},  {64,54,44,40,12,20,16,4,8},
		      {65,55,45,29,25,21,33,-1,-1},{66,56,46,30,26,22,34,-1,-1},
		      {67,57,47,31,27,23,35,-1,-1},{68,58,48,32,28,24,36,-1,-1}};
  int colDsl=kBlue-8, colDl=kCyan+1, colBkg=kYellow-7;
  int Colors[2][9] = {{colBkg,colBkg,colBkg,colBkg,colDsl,28,colDl,kGreen-7,kRed-7},
  		      {colBkg,colBkg,colBkg,colDl,colDsl,28,0,0,0}};

  // int colDsl=kBlue-8, colDl=kCyan+1, colBkg=kGreen-9;
  // int Colors[2][9] = {{colBkg,colBkg,colBkg,colBkg,colDsl,28,colDl,kRed-7, kOrange-4},
  // 		      {colBkg,colBkg,colBkg,colDl,colDsl,28,0,0,0}};
  TString channelTitle[4]={"D0","D*0","D+","D*+"};
  //TString PadLabel[] = {"a", "b", "(c)", "(d)", "e)", "f)", "g)", "h)"};
  TH1F *htree[8], pdfProj[8][9], *pdfProjCopy[8][9];
  double chi2[8][2] = {{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}}; 
  int ndof[8][2] = {{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}};
  int IndPad[2][2][4] = {{{1,2,3,4},{1,2,1,2}},{{1,3,5,7},{1,3,1,3}}};
  TH1F hSum("hSum","",Projbins,ProjMin,ProjMax), hSumPl("hSumPl","",ProjbinsPl,0,insetPlmax);
  hSum.SetFillColor(0);
  hSum.SetLineColor(0);
  TString hname;
  TLatex label; label.SetNDC(kTRUE); label.SetTextFont(132);

  gStyle->SetHatchesLineWidth(1);
  int fillDl = 3654, fillDsl = 3645;
  styles style; 
  style.LabelSize = 0.126; style.nDivisions = 708;
  style.setDefaultStyle();
  gStyle->SetNdivisions(505, "y");
  double canH = (isIso ? 900:1200);
  if(typePlot.Contains("Ds")) canH = 1100;
  if(typePlot.Contains("Inset")) canH = 800;
  TCanvas can("can","Final fit", 800, canH);
  TPad *Pads[4];
  //gPad->SetTicks(1,1);

  int maxProj = 8; if(isDss) maxProj = 6;
  double bMargin = (isIso?0.12:0.065);
  double tMargin = 0.04;
  double padH = (1-bMargin-tMargin)/nPads;
  double LeftMargin = 0.14, RightMargin = 0.1, TopMargin=0, BottomMargin=0;
  //if(var=="candM2"){LeftMargin = 0.16; RightMargin = 0.08;}
  LeftMargin = 0.18; RightMargin = 0.065;
  for(int pad=begChan; pad<endChan; pad++){
    if(isSum) IndPad[0][isIso][pad] = 1;
    double PadY[2] = {bMargin+padH*(3-pad), bMargin+padH*(4-pad)};
    if(isIso){
      PadY[0] = bMargin+padH*(1-pad%2);
      PadY[1] = bMargin+padH*(2-pad%2);
    }
    if(pad==endChan-1){
      PadY[0] = 0;
      BottomMargin=bMargin/(bMargin+padH);
    } else {
      if(typePlot.Contains("Ds")) PadY[0] += 0.05;
    }
    hname = "Pad0_"; hname += pad;
    //cout<<"Creating pad "<<pad<<" y1,y2 "<<PadY[0]<<","<<PadY[1]<<endl;
    //if(!isIso || pad>=2) {
    can.cd();
      Pads[pad] = new TPad(hname,"",0, PadY[0], 1, PadY[1]);
      Pads[pad]->SetLeftMargin(LeftMargin);     Pads[pad]->SetRightMargin(RightMargin); 
      Pads[pad]->SetBottomMargin(BottomMargin); Pads[pad]->SetTopMargin(TopMargin); 
      Pads[pad]->Draw(); Pads[pad]->cd();
      //}

    int type = 0;
    if(!isDss) {if(pad%2==1) type=1;}
    else type=2;
    int channel = pad+4*isDss;
    int candType = pad+1; if(isDss)candType+=4;
    hname = "hData"; hname += pad; 
    TString plotvar = var; plotvar+=">>"; plotvar+=hname; 
    htree[pad] = new TH1F(hname,"",nbins,minx,maxx);
    htree[pad]->Sumw2();
    TString totCuts = "(candType=="; totCuts+=candType;  
    totCuts+="&&"; totCuts+=otherVar; totCuts+=">=";
    totCuts+=minY; totCuts+="&&"; totCuts+=otherVar; totCuts+="<="; totCuts+=maxY; 
    totCuts+=")*weight";
    tree->Project(hname,var,totCuts);

    //// Looping over PDFs
    for(int i=0; i<9; i++) {
      if(Indpdf[channel][i]<0) break;
      TString SubComp = ""; SubComp += i;
      TString pdfName = "pdfProj"; pdfName += pad; pdfName += i;
      int Bins = Projbins, BinsPl = ProjbinsPl;
      double MinX = ProjMin, MaxX = ProjMax;
      if(Subtract.Contains(SubComp)){
	Bins = nbins; BinsPl = nbinsPl;
	MinX=minx; MaxX=maxx;
      }
      double scale = iYields[Indpdf[channel][i]][0]*static_cast<double>(Bins)/
	pdf[Indpdf[channel][i]]->Integral()/static_cast<double>(nbins)*(maxx-minx)/(MaxX-MinX);
      //cout<<"Yield "<<Yields[channel+truInd][IndYield[type][i]]<<" for "<< channel+truInd<<" and "<<
      //IndYield[type][i] <<", scale "<<scale<<", Bins "<<Bins<<
      //", Integral "<< pdf[Indpdf[channel][i]]->Integral()<<endl;
      pdfProj[pad][i] = binHisto(pdf[Indpdf[channel][i]],Bins,MinX,MaxX,minY,maxY,pdfName,var);
      pdfProj[pad][i].Scale(scale);
      if(Subtract.Contains(SubComp)){
	for(int bin = 1; bin<=nbins; bin++) pdfProj[pad][i].SetBinError(bin,0);
	htree[pad]->Add(&pdfProj[pad][i],-1);
	pdfProj[pad][i].Scale(0);
	doSubtract = 1;
      } else {
	hSum.Add(&pdfProj[pad][i]);
	pdfProj[pad][i].Add(&hSum,&hSum,1,0);
      }
      pdfProj[pad][i].SetFillColor(Colors[isDss][i]);
      pdfProj[pad][i].SetLineColor(Colors[isDss][i]);
      if(i==1) {pdfProj[pad][i].SetLineColor(1);pdfProj[pad][i].SetLineStyle(2);}
      else pdfProj[pad][i].SetLineWidth(0);
    } // Loop over PDFs

    style.setMarkers(htree[pad], 0.7, 20);
    if(pad >1 && isIso){
      vector<TString> IsoTags = {"D","D*"};
      if(texLabels) {
	IsoTags = vector<TString>({"Dl","D#lower[.4]{^{*}}l"});
	if(var=="candM2"){
	  if(maxFactor>5) IsoTags = vector<TString>({"(b)","(e)"});
	  else IsoTags = vector<TString>({"(a)","(d)"});
	} else IsoTags = vector<TString>({"(c)","(f)"});
      }
      channelTitle[pad] = IsoTags[pad-2];
      int iPad = pad-2; 
      htree[pad]->Add(htree[iPad]);
      for(int i=maxProj; i>=0; i--) pdfProj[pad][i].Add(&pdfProj[iPad][i]);
    }
    if(pad == 3 && isSum){
      channelTitle[pad] = "All";
      for(int iPad = 0; iPad<3; iPad++){
	htree[pad]->Add(htree[iPad]);
	for(int i=maxProj; i>=0; i--) pdfProj[pad][i].Add(&pdfProj[iPad][i]);
      }
    }
    if(typePlot.Contains("histo")){
      int iniBin = static_cast<int>((1-minx)*static_cast<double>(nbins)/(maxx-minx));
      int finBin = static_cast<int>((2-minx)*static_cast<double>(nbins)/(maxx-minx));
      //cout<<"iniBin "<<iniBin<<" and finBin "<<finBin<<endl;
      for(int bin=1; bin<=nbins; bin++){
	hIntegral[pad][0] += htree[pad]->GetBinContent(bin);
	hIntegral[pad][1] += pdfProj[pad][maxProj].GetBinContent(bin);
	double uncert = htree[pad]->GetBinError(bin);
	if(uncert<sqrt(8)) continue;
	double diff = htree[pad]->GetBinContent(bin)-pdfProj[pad][maxProj].GetBinContent(bin);
	if(uncert>0){
	  chi2[pad][0] += pow(diff/uncert,2);
	  ndof[pad][0]++;
	  if(var=="candM2"&&bin>=iniBin && bin<=finBin){
	    //cout<<"diff "<<diff<<" and uncer "<<uncert<<endl;
	    chi2[pad][1] += pow(diff/uncert,2);
	    ndof[pad][1]++;
	  }
	}
      }
      double nFit = hIntegral[pad][1], nData = hIntegral[pad][0], eData = sqrt(hIntegral[pad][0]); if(nData<=0) nData = 1;
      double Ratio = nData/nFit;
      cout<<"Data: "<<RoundNumber(nData,0)<<", Fit: "<<RoundNumber(nFit,0)<<"\t  -  Ratio: "<< RoundNumber(Ratio,5) 
	  << " +- "<< RoundNumber(eData/nFit,5)<<"\t That's "<< RoundNumber(Ratio-1,2,eData/nFit)<<" sigma away"<<endl;

    }
    float maxi = htree[pad]->GetMaximum(); 
    if(pdfProj[pad][maxProj].GetMaximum()>maxi) maxi = pdfProj[pad][maxProj].GetMaximum();
    htree[pad]->SetMaximum(1.13*maxi/maxFactor);
    if(var=="candM2" && maxFactor>5 && pad%2==1) htree[pad]->SetMaximum(1.13*maxi/maxFactor/3.);
    //if(pad==3)  htree[pad]->SetMaximum(1.13*maxi/maxFactor/3.15);
    if(doSubtract==0) htree[pad]->SetMinimum(0);
    else htree[pad]->SetMinimum(htree[pad]->GetMinimum()*1.6);
    if(typePlot.Contains("Log")) {
      htree[pad]->SetMinimum(0.9);
      Pads[pad]->SetLogy(1);
    }
    if(typePlot.Contains("pull") && ((pad>1 && isIso) || (isIso==0&&(isSum==0 || (isSum==1&&pad==3))))){
      for(int bin=1; bin<=nbins; bin++){
	double BinErr = htree[pad]->GetBinError(bin);
	if(BinErr) htree[pad]->SetBinContent(bin,(htree[pad]->GetBinContent(bin)-pdfProj[pad][maxProj].GetBinContent(bin))/BinErr);
	else htree[pad]->SetBinContent(bin,0);
	htree[pad]->SetBinError(bin,0);
      }
      htree[pad]->SetMaximum(4);
      htree[pad]->SetMinimum(-4);
      
    }
    if(pad==3) htree[pad]->SetLabelSize(style.LabelSize*padH/(padH+bMargin),"xy");
    else htree[pad]->SetLabelSize(style.LabelSize,"xy");
    if(typePlot.Contains("Ds") && !typePlot.Contains("Inset")) htree[pad]->SetLabelSize(style.LabelSize/1.8,"xy");

    htree[pad]->SetLabelOffset(0.01,"Y");
    if(!isIso || pad>=2) htree[pad]->Draw("e1");
    if(pad<3 && typePlot.Contains("Sum")){
      chi2[3][0] += chi2[pad][0]; ndof[3][0] += ndof[pad][0];
      chi2[7][0] += chi2[pad+4][0]; ndof[7][0] += ndof[pad+4][0];}
    if((var=="candM2" && (minY>0||maxY<2.4)) || ((var=="candPstarLep") && (minY>-4||maxY<12))){
      TString nameVari = "m^{2}_{miss}", otherVari = "p*_{l}";
      double minLabel = plmin, maxLabel = plmax;
      if(var=="candPstarLep"){
	nameVari = "p*_{l}"; otherVari = "m^{2}_{miss}"; 
	minLabel = m2min; maxLabel = m2max;
      }
    }
    TString optionPlot = "lf2 same";
    for(int i=maxProj; i>=0; i--) {
      TString SubComp = ""; SubComp += i;
      if(!typePlot.Contains("pull") && !Subtract.Contains(SubComp)) {
	if(!isIso || pad>=2) pdfProj[pad][i].Draw(optionPlot);
	if(Colors[isDss][i] == colDl || Colors[isDss][i] == colDsl){
	  pdfProjCopy[pad][i] = static_cast<TH1F*>(pdfProj[pad][i].Clone()); pdfProjCopy[pad][i]->SetFillColor(1); 
	  if(Colors[isDss][i] == colDsl) pdfProjCopy[pad][i]->SetFillStyle(fillDsl); 
	  else pdfProjCopy[pad][i]->SetFillStyle(fillDl); 
	  //pdfProjCopy[pad][i]->Draw(optionPlot);	  
	  if(!isIso || pad>=2) pdfProjCopy[pad][i]->Draw("h same");	  
	}
      }
    }
    hSum.Scale(0); hSumPl.Scale(0);
    if(!typePlot.Contains("pull")){
      hSum.Draw("same");
      htree[pad]->Draw("e1 sameaxis");
      htree[pad]->Draw("e0 e1 same");
    }
    style.setTitles(htree[pad],"");
  } // Pad for

  // Plotting labels
  TString ytitle = "Events/(", xTitle = "m2miss";
  if(texLabels) xTitle = "m^{2}_{miss} [GeV^{2}]";
  if(var=="candM2") {
    ytitle += RoundNumber((maxx-minx),2,static_cast<double>(nbins)); ytitle += " GeV^{2})";
  } else {
    //ytitle += RoundNumber((maxx-minx)*1000,0,static_cast<double>(nbins)); ytitle += " MeV)";
    ytitle += RoundNumber((maxx-minx),1,static_cast<double>(nbins)); ytitle += " GeV)";
    xTitle = "p";
    if(texLabels) xTitle = "E#lower[.6]{^{*}}#kern[-.9]{#lower[-.0]{_{l}}} [GeV]"; 
  }
  can.cd(0);
  label.SetTextSize(0.056); 
  label.SetTextAngle(90); label.SetTextAlign(23);
  if(typePlot.Contains("Ds")) {
    if(!typePlot.Contains("Inset"))label.DrawLatex(0.004,0.38,ytitle);
    else  label.DrawLatex(0.05,0.33,ytitle);
  } else label.DrawLatex(0.004,0.55,ytitle);
  label.SetTextAngle(0);  label.SetTextAlign(21);
  if(typePlot.Contains("Ds") && !typePlot.Contains("Inset")) label.DrawLatex(0.55,0.05,xTitle);
  else label.DrawLatex(0.55,0.02,xTitle);

  TBox box; box.SetFillStyle(1001); box.SetLineColor(10); box.SetFillColor(10);
  label.SetTextAlign(13); label.SetTextSize(style.LabelSize/3); 
  if(isIso) label.SetTextSize(style.LabelSize/2.4);
  double labLeft = 0.035;
  if(var=="candPstarLep" && isDss) labLeft = 0.6;
  //if(var=="candPstarLep" && typePlot.Contains("Sig")) labLeft = 0.65;
  for(int pad=0; pad<4; pad++){
    if(!typePlot.Contains("Ds")){
      if(!isIso){
        if(pad<3){
          box.DrawBox(LeftMargin-0.03, bMargin+padH*(3-pad), LeftMargin-0.002, bMargin+padH*(3-pad)+0.02);
          label.DrawLatex(LeftMargin-0.025, bMargin+padH*(3-pad)+0.007,"0");
        }
        label.DrawLatex(LeftMargin+labLeft, bMargin+padH*(4-pad)-0.03, channelTitle[pad]);
      } else if(pad>=2){
        if(pad<3){
          box.DrawBox(LeftMargin-0.07, bMargin+padH*(1-pad%2), LeftMargin-0.002, bMargin+padH*(1-pad%2)+0.02);
          label.DrawLatex(LeftMargin-0.032, bMargin+padH*(1-pad%2)+0.018,"0");
        }
        label.DrawLatex(LeftMargin+labLeft, bMargin+padH*(2-pad%2)-0.03, channelTitle[pad]);
      }

    }
  }

  int legPad = 3, nEntries=0;
  vector<TString> legLabel = {"Bkg", "Bkg", "Bkg", "Bkg", "Dssl", "Dsl", "Dl", "Dstau", "Dtau"};
  if(texLabels) legLabel = vector<TString>({"Bkg.", "Bkg.", "Bkg.", "Bkg.", "B #rightarrow D#lower[.4]{^{**}}(l/#tau)#nu",
					   "B #rightarrow D#lower[.4]{^{*}}l#nu", "B #rightarrow Dl#nu",
	"B #rightarrow D#lower[.4]{^{*}}#tau#nu","B #rightarrow D#tau#nu"});
  double legXY[2][2] = {{0.56, 0.77}, {bMargin+padH*(3-legPad+0.19), bMargin+padH*(3-legPad+0.87)}};
  if(isDss) {
    legXY[1][0] = bMargin+padH*(3-legPad+0.52);
    legLabel[3] = "Dl"; legLabel[4] = "Dsl"; legLabel[5] = "Dssl"; 
    maxProj--;
  }
  TLegend leg(legXY[0][0], legXY[1][0], legXY[0][1], legXY[1][1]);
  leg.SetTextSize(style.LabelSize/3.2); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont); leg.SetBorderSize(0);
  if(isIso) leg.SetTextSize(style.LabelSize/2.5);
  leg.AddEntry(htree[0], "Data"); nEntries++;
  for(int ipdf=maxProj; ipdf>0; ipdf--) {
    if(ipdf==2 || (ipdf==3&&isDss==0)) continue;
    int iipdf = ipdf;
    if(isDss==0){
      if(ipdf==4) iipdf = 5;
      if(ipdf==5) iipdf = 4;
    } 
    TString option = "f";
    if(pdfProj[0][iipdf].GetLineStyle() == 2) option = "fl";
    leg.AddEntry(&pdfProj[0][iipdf], legLabel[ipdf],option);
    nEntries++;
  }
  if(typePlot.Contains("Leg")) {
    can.cd(0);
    leg.Draw();
    double iDl = 3, iDsl = 2;
    if(isDss) iDl = 1; 
    fillDl -= 300; fillDsl -= 300; 
    box.SetFillStyle(fillDl); box.SetLineColor(0); box.SetFillColor(1);
    double legW = legXY[0][1]-legXY[0][0], legH = (legXY[1][1]-legXY[1][0])/static_cast<double>(nEntries);
    box.DrawBox(legXY[0][0]+legW*0.04, legXY[1][0]+legH*(iDl+0.15), 
     		legXY[0][0]+legW*0.21, legXY[1][0]+legH*(iDl+0.85));
    box.SetFillStyle(fillDsl); 
    box.DrawBox(legXY[0][0]+legW*0.04, legXY[1][0]+legH*(iDsl+0.15), 
     		legXY[0][0]+legW*0.21, legXY[1][0]+legH*(iDsl+0.85));
    box.SetFillColor(0); box.SetFillStyle(0); box.SetLineWidth(1); box.SetLineColor(1); 
    for(int row=0; row<nEntries-1; row++)
      box.DrawBox(legXY[0][0]+legW*0.04, legXY[1][0]+legH*(row+0.15), 
    		  legXY[0][0]+legW*0.21, legXY[1][0]+legH*(row+0.85));
  }
  
  double cutsH = 0.51, cutsX = 0.9;
  label.SetTextAlign(33); label.SetTextSize(style.LabelSize/2.5); 
  //label.DrawLatex(cutsX, cutsH, "q^{2} > 4 GeV^{2}");
  //if(minY>=1) label.DrawLatex(cutsX, cutsH-0.05, "m^{2}_{miss} > "+RoundNumber(minY,0)+" GeV^{2}");
  if(minY>=1) {
    label.DrawLatex(cutsX, 1-tMargin-0.08, "m^{2}_{miss} > "+RoundNumber(minY,0)+" GeV^{2}");
    label.DrawLatex(cutsX, cutsH, "m^{2}_{miss} > "+RoundNumber(minY,0)+" GeV^{2}");
  }

  BABARLabel(0.75,0.9,0.62);
  if(typePlot.Contains("Ds") && !typePlot.Contains("Inset")) BABARLabel(0.75,0.48,0.62);
  TString plotName = "plots/PRL"; plotName+=sample; plotName+=var; plotName+="_"; 
  plotName+=typePlot; plotName+=".pdf";
  can.SaveAs(plotName);
  return;
}



