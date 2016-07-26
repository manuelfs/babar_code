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
#include "DonutUtils/Styles.cc"
/////////////////////////////
#include "TString.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLatex.h"
#include "TLine.h"
#include "DonutUtils/KeysUtils.cc"
#include <fstream>
#include <iostream>

using namespace TMath;
using namespace std;
using std::cout;
using std::endl;

void PlotFinalFit2(TString textName, TTree *tree, int nbins = 25, double minx=-0.4, double maxx=1.4, double minY=0, double maxY=2.4,
		  TString typePlot="histo", TString var="candM2", TString sample="sig", 
		  double maxFactor = 1., TString psEPS = "EPS", double maxFactorPl = 1., int nbinsPl = 25);

int main(int argc, char *argv[]){
  if (argc < 2 || argc > 14 ) {
    cout << "USAGE: PlotFit textFile [sample=sig] [Subtract=] [var=candM2] [nbins=25] "<< 
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
  TString genFile = "AWG82/ntuples/small/Fit"; genFile += TagSample; genFile += "_RunAll.root";
  if(textName.Contains("Data")) {
    //if(!TagSample.Contains("Sigx")) genFile = "AWG82/ntuples/small/FitData_RunAll.root";
    if(textName.Contains("1234")) genFile = "AWG82/ntuples/small/FitData_Run1234.root";
    if(textName.Contains("56")) genFile = "AWG82/ntuples/small/FitData_Run56.root";
    if(textName.Contains("eData")) genFile = "AWG82/ntuples/small/FiteDataSidx100_RunAll.root";
    if(textName.Contains("muData")) genFile = "AWG82/ntuples/small/FitmuDataSidx100_RunAll.root";
  }
  if(textName.Contains("Dpipi")) {
    genFile = "AWG82/ntuples/small/FitDataSidx100_RunAll.root";
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

  int truInd=0, isDss = 0, doSubtract = 0, doPl = 0, isIso = 0, isSum = 0, nPads = 4;
  if(sample.Contains("dss")) isDss=1;
  if(typePlot.Contains("Tru")) truInd=8;
  if(typePlot.Contains("Sum")) {isSum = 1; nPads  = 1;}
  if(typePlot.Contains("Pl"))  {doPl = 1;  nPads *= 2;}
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
  fstream ResWFile; ResWFile.open("../DonutUtils/ResolutionWidths.txt",fstream::in);
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
    TString hname = "keys/root/Fit/pdfKeys_"; 
    if(TagSample.Contains("x")) {hname = "keys/root/fit"; hname += TagSample; hname += "/pdfKeys_"; }
    hname += i; hname += "_Fit.root";
    TFile hfile(hname); 
    TString pdfName = "pdf"; pdfName += i;
    pdf[i] = (TH2F *)(hfile.Get("h2"))->Clone(pdfName);    
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
  //int Colors[7] = {5,28,4,2,3,38,43};
  int Colors[2][9] = {{10,10,10,0,kBlue-9,kRed+4,kYellow-4,kGreen-7,kRed-7},{5,5,5,2,4,28,38,0,0}};
  TString channelTitle[4]={"D^{0}","D*^{0}","D^{+}","D*^{+}"};
  TString PadLabel[] = {"a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)"};
  TH1F *htree[8], pdfProj[8][9], *pdfProjCopy[8][9];
  double chi2[8][2] = {{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}}; 
  int ndof[8][2] = {{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}};
  int IndPad[2][2][4] = {{{1,2,3,4},{1,2,1,2}},{{1,3,5,7},{1,3,1,3}}};
  TH1F hSum("hSum","",Projbins,ProjMin,ProjMax), hSumPl("hSumPl","",ProjbinsPl,0,insetPlmax);
  hSum.SetFillColor(0);
  hSum.SetLineColor(0);
  TString hname;
  TLine line; line.SetLineWidth(1);
  TLatex label; label.SetNDC(kTRUE); label.SetTextFont(132);

  Styles style; style.setPadsStyle(4); 
  style.PadTopMargin = 0; style.PadBottomMargin = 0; style.PadRightMargin = 0.01; 
  style.yTitleOffset = 0.015; style.LabelSize = 0.08; style.nDivisions = 308;
  style.applyStyle();
  gStyle->SetNdivisions(205, "y");   
  TCanvas can("can","Final fit",800,1500);
  TPad *mainpad = new TPad("mainpad","",0,0.03,1.0,1.0), *inset[4];
  mainpad->Draw();
  mainpad->cd();
  mainpad->Divide(1,4,0,0);
  gPad->SetTicks(1,1);


  for(int pad=begChan; pad<endChan; pad++){
    if(isSum) IndPad[doPl][isIso][pad] = 1;
    mainpad->cd(pad+1);
    if(pad==3) gPad->SetBottomMargin(0.08);
    int type = 0;
    if(!isDss) {if(pad%2==1) type=1;}
    else type=2;
    int channel = pad+4*isDss;
    int candType = pad+1; if(isDss)candType+=4;
    hname = "hData"; hname += pad; 
    TString plotvar = var; plotvar+=">>"; plotvar+=hname; 
    htree[pad] = new TH1F(hname,hname,nbins,minx,maxx);
    htree[pad]->Sumw2();
    TString totCuts = "(candType=="; totCuts+=candType;  
    if(doPl==0){
      totCuts+="&&"; totCuts+=otherVar; totCuts+=">=";
      totCuts+=minY; totCuts+="&&"; totCuts+=otherVar; totCuts+="<="; totCuts+=maxY; 
    }
    totCuts+=")*weight";
    tree->Project(hname,var,totCuts);
    if(doPl){
      hname += "Pl"; 
      plotvar ="candPstarLep>>"; plotvar+=hname; 
      htree[pad+4] = new TH1F(hname,hname,nbinsPl,0,insetPlmax);
      htree[pad+4]->Sumw2();
      totCuts = "(candType=="; totCuts+=candType;  totCuts+="&&"; totCuts+="candM2>=";
      totCuts+=minY; totCuts+="&&"; totCuts+="candM2<="; totCuts+=maxY; totCuts+=")*weight";
      tree->Project(hname,"candPstarLep",totCuts);
    }
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
      double scale = iYields[Indpdf[channel][i]][0]*(double)Bins/
	pdf[Indpdf[channel][i]]->Integral()/(double)nbins*(maxx-minx)/(MaxX-MinX);
      //cout<<"Yield "<<Yields[channel+truInd][IndYield[type][i]]<<" for "<< channel+truInd<<" and "<<
      //IndYield[type][i] <<", scale "<<scale<<", Bins "<<Bins<<
      //", Integral "<< pdf[Indpdf[channel][i]]->Integral()<<endl;
      if(doPl) pdfProj[pad][i] = binHisto(pdf[Indpdf[channel][i]],Bins,MinX,MaxX,0,insetPlmax,pdfName,var);
      else pdfProj[pad][i] = binHisto(pdf[Indpdf[channel][i]],Bins,MinX,MaxX,minY,maxY,pdfName,var);
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
      if(doPl){
	pdfName += "Pl";
	scale = iYields[Indpdf[channel][i]][0]*(double)BinsPl/
	  pdf[Indpdf[channel][i]]->Integral()/(double)nbinsPl;
	pdfProj[pad+4][i] = binHisto(pdf[Indpdf[channel][i]],BinsPl,0,insetPlmax,minY,maxY,pdfName,"candPstarLep");
	pdfProj[pad+4][i].Scale(scale);
	if(Subtract.Contains(SubComp)){
	  for(int bin = 1; bin<=nbinsPl; bin++) pdfProj[pad+4][i].SetBinError(bin,0);
	  htree[pad+4]->Add(&pdfProj[pad+4][i],-1);
	  pdfProj[pad+4][i].Scale(0);
	} else {
	  hSumPl.Add(&pdfProj[pad+4][i]);
	  pdfProj[pad+4][i].Add(&hSumPl,&hSumPl,1,0);
	}
	pdfProj[pad+4][i].SetFillColor(Colors[isDss][i]);
	pdfProj[pad+4][i].SetLineColor(Colors[isDss][i]);
	if(i==1) {pdfProj[pad+4][i].SetLineColor(1);pdfProj[pad+4][i].SetLineStyle(2);}
      }
    }

    int maxProj = 8; if(isDss) maxProj = 6;
    style.setMarkers(htree[pad], 0.6, 20);
    if(pad >1 && isIso){
      TString IsoTags[] = {"D","D*"};
      channelTitle[pad] = IsoTags[pad-2];
      int iPad = pad-2; 
      htree[pad]->Add(htree[iPad]);
      for(int i=maxProj; i>=0; i--) pdfProj[pad][i].Add(&pdfProj[iPad][i]);
      if(doPl){
	htree[pad+4]->Add(htree[iPad+4]);
	for(int i=maxProj; i>=0; i--) pdfProj[pad+4][i].Add(&pdfProj[iPad+4][i]);
      }
    }
    if(pad == 3 && isSum){
      channelTitle[pad] = "All";
      for(int iPad = 0; iPad<3; iPad++){
	htree[pad]->Add(htree[iPad]);
	for(int i=maxProj; i>=0; i--) pdfProj[pad][i].Add(&pdfProj[iPad][i]);
	if(doPl){
	  htree[pad+4]->Add(htree[iPad+4]);
	  for(int i=maxProj; i>=0; i--) pdfProj[pad+4][i].Add(&pdfProj[iPad+4][i]);
	}
      }
    }
    if(typePlot.Contains("histo")){
      int iniBin = (int)((1-minx)*(double)nbins/(maxx-minx)), finBin = (int)((2-minx)*(double)nbins/(maxx-minx));
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

      if(doPl){
	for(int bin=1; bin<=nbinsPl; bin++){
	  hIntegral[pad+4][0] += htree[pad+4]->GetBinContent(bin);
	  hIntegral[pad+4][1] += pdfProj[pad+4][maxProj].GetBinContent(bin);
	  double uncert = htree[pad+4]->GetBinError(bin);
	  if(uncert<sqrt(8)) continue;
	  double diff = htree[pad+4]->GetBinContent(bin)-pdfProj[pad+4][maxProj].GetBinContent(bin);
	  if(uncert>0){
	    chi2[pad+4][0] += pow(diff/uncert,2);
	    ndof[pad+4][0]++;
	    if(var=="candM2"&&bin>=iniBin && bin<=finBin){
	      chi2[pad+4][1] += pow(diff/uncert,2);
	      ndof[pad+4][1]++;
	    }
	  }
	}
      }
    }
    float maxi = htree[pad]->GetMaximum(); 
    if(pdfProj[pad][maxProj].GetMaximum()>maxi) maxi = pdfProj[pad][maxProj].GetMaximum();
    htree[pad]->SetMaximum(1.13*maxi/maxFactor);
    if(var=="candM2" && maxFactor>5 && pad%2==1) htree[pad]->SetMaximum(1.13*maxi/maxFactor/3.);
    if(pad==3)  htree[pad]->SetMaximum(1.13*maxi/maxFactor/3.15);
    if(doSubtract==0) htree[pad]->SetMinimum(0);
    else htree[pad]->SetMinimum(htree[pad]->GetMinimum()*1.6);
    if(typePlot.Contains("pull") && (pad>1 && isIso || isIso==0&&(isSum==0||isSum==1&&pad==3))){
      for(int bin=1; bin<=nbins; bin++){
	double BinErr = htree[pad]->GetBinError(bin);
	if(BinErr) htree[pad]->SetBinContent(bin,(htree[pad]->GetBinContent(bin)-pdfProj[pad][maxProj].GetBinContent(bin))/BinErr);
	else htree[pad]->SetBinContent(bin,0);
	htree[pad]->SetBinError(bin,0);
      }
      htree[pad]->SetMaximum(4);
      htree[pad]->SetMinimum(-4);
      
    }
    htree[pad]->Draw("e1");
    if(doPl){
      hname = "inset"; hname += pad;
      inset[pad] = new TPad(hname,"",0.57,0.34,0.99,0.98);
      inset[pad]->SetRightMargin(0); inset[pad]->SetTopMargin(0);
      inset[pad]->SetLeftMargin(0.15); inset[pad]->SetBottomMargin(0.34);
      inset[pad]->Draw(); inset[pad]->cd();

      //can.cd(IndPad[doPl][isIso][pad]+1);
      maxi = htree[pad+4]->GetMaximum(); 
      if(pdfProj[pad+4][maxProj].GetMaximum()>maxi) maxi = pdfProj[pad+4][maxProj].GetMaximum();
      htree[pad+4]->SetMaximum(1.13*maxi/maxFactorPl);
      if(doSubtract==0) htree[pad+4]->SetMinimum(0);
      else htree[pad+4]->SetMinimum(htree[pad+4]->GetMinimum()*1.6);
      if(typePlot.Contains("pull") && (pad>1 && isIso || isIso==0&&(isSum==0||isSum==1&&pad==3))){
	for(int bin=1; bin<=nbinsPl; bin++){
	  double BinErr = htree[pad+4]->GetBinError(bin);
	  if(BinErr) htree[pad+4]->SetBinContent(bin,(htree[pad+4]->GetBinContent(bin)-pdfProj[pad+4][maxProj].GetBinContent(bin))/BinErr);
	  else htree[pad+4]->SetBinContent(bin,0);
	  htree[pad+4]->SetBinError(bin,0);
	}
	htree[pad+4]->SetMaximum(4);
	htree[pad+4]->SetMinimum(-4);
      
      }
      htree[pad+4]->GetXaxis()->SetNdivisions(305); htree[pad+4]->GetYaxis()->SetNdivisions(203);
      htree[pad+4]->SetTitleOffset(1.02,"x");     
      htree[pad+4]->SetLabelSize(0.13); htree[pad+4]->SetLabelSize(0.13,"Y"); htree[pad+4]->SetTitleSize(0.15);
      htree[pad+4]->Draw("e1");
      mainpad->cd(pad+1);
    }
    if(pad<3 && typePlot.Contains("Sum")){
      chi2[3][0] += chi2[pad][0]; ndof[3][0] += ndof[pad][0];
      chi2[7][0] += chi2[pad+4][0]; ndof[7][0] += ndof[pad+4][0];}
    if(var=="candM2"&&(minY>0||maxY<2.4) || (doPl||var=="candPstarLep")&&(minY>-4||maxY<12)){
      TString nameVari = "m^{2}_{miss}", otherVari = "p*_{l}";
      double minLabel = plmin, maxLabel = plmax;
      if(var=="candPstarLep" || doPl){
	nameVari = "p*_{l}"; otherVari = "m^{2}_{miss}"; 
	minLabel = m2min; maxLabel = m2max;
      }
    }
    TString optionPlot = "lf2 same";
    for(int i=maxProj; i>=0; i--) {
      TString SubComp = ""; SubComp += i;
      if(!typePlot.Contains("pull") && !Subtract.Contains(SubComp)) {
	pdfProj[pad][i].Draw(optionPlot);
	if(i==maxProj-2 || i==maxProj-4){
	  pdfProjCopy[pad][i] = (TH1F*)pdfProj[pad][i].Clone(); pdfProjCopy[pad][i]->SetFillColor(1); 
	  if(i==maxProj-2) pdfProjCopy[pad][i]->SetFillStyle(3004); 
	  else pdfProjCopy[pad][i]->SetFillStyle(3005); 
	  pdfProjCopy[pad][i]->Draw(optionPlot);	  
	}
      }
    }
    hSum.Scale(0); hSumPl.Scale(0);
    if(!typePlot.Contains("pull")){
      hSum.Draw("same");
      htree[pad]->Draw("e1 sameaxis");
      htree[pad]->Draw("e0 e1 same");
    }
    if(doPl){
      inset[pad]->cd(); 
      for(int i=maxProj; i>=0; i--) {
	TString SubComp = ""; SubComp += i;
	if(!typePlot.Contains("pull") && !Subtract.Contains(SubComp)) {
	  pdfProj[pad+4][i].Draw(optionPlot);
	  if(i==maxProj-2 || i==maxProj-4){
	    pdfProjCopy[pad+4][i] = (TH1F*)pdfProj[pad+4][i].Clone(); pdfProjCopy[pad+4][i]->SetFillColor(1); 
	    if(i==maxProj-2) pdfProjCopy[pad+4][i]->SetFillStyle(3004); 
	    else pdfProjCopy[pad+4][i]->SetFillStyle(3005); 
	    pdfProjCopy[pad+4][i]->Draw(optionPlot);	  
	}
      }
      }
      if(!typePlot.Contains("pull")){
	hSum.Draw("same");
	htree[pad+4]->Draw("e1 sameaxis");
	htree[pad+4]->Draw("e0 e1 same");
      }
      style.setTitles(htree[pad+4],"|#font[22]{p}*_{l}| (GeV)");
      if(typePlot.Contains("psfrag")) style.setTitles(htree[pad+4],"p");
    }
    double padLabelX = 0.27;
    if(typePlot.Contains("psfrag")) {
      padLabelX = 0.245;
    }
    mainpad->cd(pad+1);
    label.SetTextSize(0.12); label.SetTextAlign(33); label.DrawLatex(padLabelX,0.9,PadLabel[pad]);

    style.setTitles(htree[pad],"");
  }  
  TString ytitle = "Events/(", xTitle = "m^{2}_{miss} (GeV^{2})";
  double xTitleY = 0.029;
  if(var=="candM2") {
    ytitle += RoundNumber((maxx-minx),2,(double)nbins); ytitle += " GeV^{2})    [Events/(100 MeV)  in insets]";
  } else {
    ytitle += RoundNumber((maxx-minx)*1000,0,(double)nbins); ytitle += " MeV)";
    xTitle = "p*_{l} (GeV)";
  }
  mainpad->cd();
  label.SetTextSize(0.053); 
  if(typePlot.Contains("psfrag")) {
    xTitle = "mmiss_gev";
    label.SetTextSize(0.058); 
    xTitleY = 0.023;
  }
  label.SetTextAngle(90); label.DrawLatex(0.033,0.81,ytitle);
  can.cd();
  label.SetTextSize(0.049); 
  label.SetTextAngle(0); label.SetTextAlign(33); label.DrawLatex(0.9,xTitleY,xTitle);
  TString plotName = "keys/eps/FinalFit/PRL"; plotName+=sample; plotName+=var; plotName+="_"; 
  TBox box; box.SetFillStyle(1001); box.SetLineColor(10);box.SetFillColor(10);
  box.DrawBox(0.158,0.75, 0.176,0.77);  box.DrawBox(0.158,0.51, 0.176,0.53);  box.DrawBox(0.158,0.26, 0.176,0.28);
  label.SetTextSize(0.038); 
  label.DrawLatex(0.174,0.767,"0");
  label.DrawLatex(0.174,0.522,"0");
  label.DrawLatex(0.174,0.28,"0");
  plotName+=typePlot; plotName+=".eps";
  can.SaveAs(plotName);
  return;
}



