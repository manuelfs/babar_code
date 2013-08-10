//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: PlotFit.cc,v 1.3 2011/02/22 01:01:50 manuelf Exp $
//
// Description:
//      PlotPD - Plots the results of the fit to PD
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      11/06/29 manuelf -- Adapted from PlotFit.cc
//------------------------------------------------------------------------
//////////////////////////////
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "DonutUtils/Styles.cc"
/////////////////////////////
#include "TString.h"
#include "TChain.h"
#include "TCanvas.h"
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
  TagSample.Remove(0,TagSample.Index("Final")+5); TagSample.ReplaceAll(".txt","");
  TString genFile = "AWG82/ntuples/small/Fit"; genFile += TagSample; genFile += "_RunAll.root";
  if(textName.Contains("Data")) {
    genFile = "AWG82/ntuples/small/FitData_RunAll.root";
    if(textName.Contains("1234")) genFile = "AWG82/ntuples/small/FitData_Run1234.root";
    if(textName.Contains("56")) genFile = "AWG82/ntuples/small/FitData_Run56.root";
    if(textName.Contains("eData")) genFile = "AWG82/ntuples/small/FiteData_RunAll.root";
    if(textName.Contains("muData")) genFile = "AWG82/ntuples/small/FitmuData_RunAll.root";
  }
  TChain *gen = new TChain("ntp1");
  gen->Add(genFile);
  if(Subtract.Contains("ps")){
    TString plotfile = "keys/eps/FinalFit/FinalFit_"; plotfile+=typePlot; plotfile+=".ps";
    cout << endl << "Plots "<<typePlot<<" of the fits saved in " << plotfile << endl;
    TCanvas canvas; // Needed to open and close the ps file
    canvas.Print(plotfile+"[");
    PlotFinalFit2(textName, gen, 30, -0.4, 1.4, 0, 2.4, typePlot, "candM2", "sig", 1., plotfile);
    PlotFinalFit2(textName, gen, 34, -0.9, 2.9, 0, 2.4, typePlot, "candM2", "sig", 1., plotfile);
    PlotFinalFit2(textName, gen, 34, -0.9, 2.9, 0, 2.4, typePlot, "candM2", "sig", 11., plotfile);
    PlotFinalFit2(textName, gen, 34, -1.8, 5.8, 0, 2.4, typePlot, "candM2", "sig", 1., plotfile);
    PlotFinalFit2(textName, gen, 34, -1.8, 5.8, 0, 2.4, typePlot, "candM2", "sig", 20., plotfile);
    PlotFinalFit2(textName, gen, 34, -2., 9.75, 0, 2.4, typePlot, "candM2", "sig", 1., plotfile);
    PlotFinalFit2(textName, gen, 34, -2., 9.75, 0, 2.4, typePlot, "candM2", "sig", 25., plotfile);
    PlotFinalFit2(textName, gen, 22, -0.9, 2.1, -4, 12, typePlot, "candM2", "dss", 1., plotfile);
    PlotFinalFit2(textName, gen, 28, -1.8, 4.2, -4, 12, typePlot, "candM2", "dss", 1., plotfile);
    PlotFinalFit2(textName, gen, 30, -3.6, 8.4, -4, 12, typePlot, "candM2", "dss", 1., plotfile);
    PlotFinalFit2(textName, gen, 25, 0., 2.4,   0, 2.4, typePlot, "candPstarLep", "sig", 1., plotfile);
    PlotFinalFit2(textName, gen, 15, 0., 2.4,   -4, 12, typePlot, "candPstarLep", "dss", 1., plotfile);
    canvas.Print(plotfile+"]");
  } else PlotFinalFit2(textName, gen, nbins, min, max, minY, maxY, typePlot, var, sample, maxFactor, Subtract, maxFactorPl, nbinsPl);

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

  Styles style; style.setPadsStyle(nPads); style.applyStyle();

  double Yields[16][9];
  TString buffer;
  fstream textFile; textFile.open(textName,fstream::in);
  int begChan = 0, endChan = 4;
  if(typePlot.Contains("D0")&&!typePlot.Contains("Dc")) endChan = 2; 
  if(typePlot.Contains("Dc")&&!typePlot.Contains("D0")) begChan = 2;
  for(int index=0; index<4; index++){
    for(int j=0; j<8; j++){
      textFile>>buffer>>Yields[index+8][j]>>Yields[index][j]>>buffer>>buffer>>buffer;
      //cout<<"True: "<<Yields[index+8][j]<<"\tFit: "<<Yields[index][j]<<"\tChannel "<<index<<endl;
    }
    //cout<<endl;
  }
  double m2cut[2];
  TString m2CutName = "babar_code/Systematics/MmissCut.txt";
  fstream m2File; m2File.open(m2CutName,fstream::in);
  m2File>>m2cut[0]>>m2cut[1];
  TLine line; line.SetLineWidth(1);
  TCanvas can("can","Final fit");
  if(nPads>1) can.Divide(2,nPads/2);
  int nm2bin = 0, nplbin = 0, Projbins = nbins, ProjbinsPl = nbinsPl;
  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4, pdmin = 0, pdmax = 2.;
  TH1F *pdf[70];
  for (int i = 1 ; i <= 64 ; i ++) {
    if(i==21)i=41; 
    if(i==45)i=51; if(i==55)i=61; 
    TString hname = "keys/root/PDfitFinal/pdfKeys_"; hname += i; hname += "_Fit.root";
    TFile hfile(hname); 
    TString pdfName = "pdf"; pdfName += i;
    pdf[i] = (TH1F *)(hfile.Get("h2"))->Clone(pdfName);    
    pdf[i]->SetDirectory(0);
    if(i==1) {nm2bin = pdf[1]->GetNbinsX();}
  }
  int nxbins = nm2bin; double ProjMin=minx, ProjMax=maxx; 
  TString otherVar = "candPstarLep";
  if(typePlot.Contains("curve")) {
    Projbins = nxbins; ProjbinsPl = nplbin;
    ProjMin=pdmin; ProjMax=pdmax;
  }
  int Indpdf[8][9] = {{61,51,41,17,13,9,5,1},   {62,52,42,18,14,10,6,2},
		      {63,53,43,19,15,11,7,3},  {64,54,44,20,16,12,8,4},
		      {65,55,45,29,25,21,-1,-1},{66,56,46,30,26,22,-1,-1},
		      {67,57,47,31,27,23,-1,-1},{68,58,48,32,28,24,-1,-1}};
  int IndYield[3][8] = {{7,6,5,4,3,1,2,0},{7,6,5,4,3,1,2,0},{5,4,3,2,1,0,-1,-1}};
  int Colors[7] = {5,28,4,2,3,38,43};
  TString legNames[7] = {"Bkg.","D**l#nu","D*l#nu","Dl#nu","D*#tau#nu","D#tau#nu","D#pi#pil#nu"};
  int IndColor[3][8] = {{0,0,0,1,2,3,4,5},{0,0,0,1,3,2,5,4},{0,0,0,3,2,1,-1,-1}};
  TString channelTitle[4]={"D^{0}","D*^{0}","D^{+}","D*^{+}"};
  TString PadLabel[] = {"a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)"};
  TH1F *htree[8], pdfProj[8][9];
  TLegend *leg[4];
  double chi2[8][2] = {{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}}; 
  int ndof[8][2] = {{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}};
  int IndPad[2][2][4] = {{{1,2,3,4},{1,2,1,2}},{{1,3,5,7},{1,3,1,3}}};
  TH1F hSum("hSum","",Projbins,ProjMin,ProjMax), hSumPl("hSumPl","",ProjbinsPl,0,2.4);
  hSum.SetFillColor(0);
  hSum.SetLineColor(0);
  for(int pad=begChan; pad<endChan; pad++){
    if(isSum) IndPad[doPl][isIso][pad] = 1;
    can.cd(IndPad[doPl][isIso][pad]);
    int type = 0;
    if(!isDss) {if(pad%2==1) type=1;}
    else type=2;
    int channel = pad+4*isDss;
    int candType = pad+1; if(isDss)candType+=4;
    TString hname = "hData"; hname += pad; 
    TString plotvar = var; plotvar+=">>"; plotvar+=hname; 
    htree[pad] = new TH1F(hname,hname,nbins,minx,maxx);
    htree[pad]->Sumw2();
    TString totCuts = "(candType=="; totCuts+=candType; totCuts+="&&candM2>";
    totCuts+=m2cut[0]; totCuts+="&&candM2<="; totCuts+=m2cut[1]; totCuts+=")*weight";
    tree->Project(hname,var,totCuts);
    for(int i=0; i<8; i++) {
      if(Indpdf[channel][i]<0) break;
      TString SubComp = ""; SubComp += i;
      TString pdfName = "pdfProj"; pdfName += pad; pdfName += i;
      int Bins = Projbins, BinsPl = ProjbinsPl;
      double MinX = ProjMin, MaxX = ProjMax;
      if(Subtract.Contains(SubComp)){
	Bins = nbins; BinsPl = nbinsPl;
	MinX=minx; MaxX=maxx;
      }
    
      double scale = Yields[channel+truInd][IndYield[type][i]]*(double)Bins/
	pdf[Indpdf[channel][i]]->Integral()/(double)nbins*(maxx-minx)/(MaxX-MinX);
      //cout<<"Yield "<<Yields[channel+truInd][IndYield[type][i]]<<" for "<< channel+truInd<<" and "<<
      //IndYield[type][i] <<", scale "<<scale<<", Bins "<<Bins<<
      //", Integral "<< pdf[Indpdf[channel][i]]->Integral()<<endl;
      pdfProj[pad][i] = *pdf[Indpdf[channel][i]];
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
      pdfProj[pad][i].SetFillColor(Colors[IndColor[type][i]]);
      pdfProj[pad][i].SetLineColor(Colors[IndColor[type][i]]);
      if(i==1) {pdfProj[pad][i].SetLineColor(1);pdfProj[pad][i].SetLineStyle(2);}
    }

    TString ytitle = "Events/(", xTitle = "p*_{D} (Gev)";
    ytitle += RoundNumber((maxx-minx)*1000,0,(double)nbins); ytitle += " MeV)";
    int maxProj = 7; if(isDss) maxProj = 6;
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
      if(doPl){
	for(int bin=1; bin<=nbinsPl; bin++){
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
    if(var=="candM2" && maxFactor>10 && pad%2==1) htree[pad]->SetMaximum(1.13*maxi/maxFactor/2);
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
    if(typePlot.Contains("pull")){
      ytitle = "Pull";
      htree[pad]->SetMarkerSize(0.75);
      htree[pad]->Draw("p");
      line.SetLineStyle(1); line.SetLineColor(4); line.DrawLine(minx,0,maxx,0);
      line.SetLineStyle(2); line.SetLineColor(2); line.DrawLine(minx,2.,maxx,2.); line.DrawLine(minx,-2.,maxx,-2.);
    }
    if(doPl){
      can.cd(IndPad[doPl][isIso][pad]+1);
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
      htree[pad+4]->Draw("e1");
      can.cd(IndPad[doPl][isIso][pad]);
    }
    double legW = 0.18, legH = 0.41;
    double legX = 1-style.PadRightMargin-0.01, legY = 1-style.PadTopMargin-0.01;
    ndof[pad][0] -= 3; ndof[pad+4][0] -= 3; 
    if(pad%2==1) {ndof[pad][0]++; ndof[pad+4][0]++;}
    if(var=="candM2"){
      if(maxx<2) ndof[pad][0]++;
      if(minx>1.2) ndof[pad][0]++;
    } 
    if(var=="candPstarLep"){
      if(minY>=1) ndof[pad][0]++;
    } 
    if(typePlot.Contains("m1")){ndof[pad][0]--; ndof[pad][1]--; ndof[pad+4][0]--; ndof[pad+4][1]--;}
    if(typePlot.Contains("m2")){ndof[pad][0]-=2; ndof[pad][1]-=2; ndof[pad+4][0]-=2; ndof[pad+4][1]-=2;}
    if(pad<3 && typePlot.Contains("Sum")){
      chi2[3][0] += chi2[pad][0]; ndof[3][0] += ndof[pad][0];
      chi2[7][0] += chi2[pad+4][0]; ndof[7][0] += ndof[pad+4][0];}
    TString rightLabel = "";
    if(var=="candM2"&&(minY>0||maxY<2.4) || (doPl||var=="candPstarLep")&&(minY>-4||maxY<12)){
      TString units="GeV^{2})", otherUnits="GeV", nameVari = "m^{2}_{miss}", otherVari = "p*_{l}";
      double minLabel = plmin, maxLabel = plmax;
      if(var=="candPstarLep" || doPl){
	units="GeV)"; nameVari = "p*_{l}"; otherVari = "m^{2}_{miss}"; otherUnits="GeV^{2}";
	minLabel = m2min; maxLabel = m2max;
      }
      rightLabel = RoundNumber(minY,1); 
      if(maxY > maxLabel) {rightLabel+=" "; rightLabel += otherUnits;}
      rightLabel+=" #leq "; 
      if(minY <= minLabel) rightLabel=otherVari; 
      else rightLabel+=otherVari; 
      if(maxY <= maxLabel) {
	rightLabel+=" < "; rightLabel += RoundNumber(maxY,1); rightLabel += " "; rightLabel += otherUnits;
      }
    }
    if(typePlot.Contains("histo")){
      rightLabel = "#splitline{#chi^{2}: "; rightLabel += RoundNumber(chi2[pad][0],1); rightLabel += "/";
      rightLabel += ndof[pad][0]; rightLabel += " = "; rightLabel += RoundNumber(chi2[pad][0],2,ndof[pad][0]);
      rightLabel += "}{Prob. = "; rightLabel += RoundNumber(Prob(chi2[pad][0],ndof[pad][0])*100,2); rightLabel += "%}";
    }
     if(channel<4){
      if(var=="candPstarLep")  {
	legX = style.PadLeftMargin+legW+0.03; 
	TString tempS = rightLabel;
	rightLabel = channelTitle[pad];
	channelTitle[pad] = tempS;
      }
    }else legH = 0.32;
    leg[pad] = new TLegend(legX-legW, legY-legH, legX, legY);
    leg[pad]->SetTextSize(style.LabelSize); leg[pad]->SetFillColor(0); leg[pad]->SetTextFont(style.nFont);
    leg[pad]->SetBorderSize(0);
    TString optionPlot = "lf2 same";
    if(!typePlot.Contains("curve")) optionPlot = "same";
    for(int i=maxProj; i>=0; i--) {
      TString SubComp = ""; SubComp += i;
      if(!typePlot.Contains("pull") && !Subtract.Contains(SubComp)) pdfProj[pad][i].Draw(optionPlot);
      if((i>2 || i==1) && (Indpdf[channel][i]<33||Indpdf[channel][i]>40||typePlot.Contains("Dpipi")))leg[pad]->AddEntry(&pdfProj[pad][i],legNames[IndColor[type][i]]);
    }
    TString pos = "right"; if(!isDss && var=="candPstarLep") pos="left";
    hSum.Scale(0); hSumPl.Scale(0);
    if(!typePlot.Contains("pull")){
      hSum.Draw("same");
      htree[pad]->Draw("e1 sameaxis");
      htree[pad]->Draw("e0 e1 same");
    }
    int didLeg = 0;
    if(typePlot.Contains("Leg") && (IndPad[doPl][isIso][pad]==begChan+1 || typePlot.Contains("Sum"))) {
      leg[pad]->Draw();
      didLeg = 1;
    }
    TString safeRLabel = rightLabel;
    can.cd(IndPad[doPl][isIso][pad]);
    rightLabel = safeRLabel;
    if(didLeg || typePlot.Contains("pull") ) rightLabel = "";	
    if(didLeg) channelTitle[pad] += "                     ";
    TString RealLabel = rightLabel; RealLabel += ",  ";
    if(rightLabel=="" || doPl) RealLabel = channelTitle[pad];
    else RealLabel += channelTitle[pad];
    //RealLabel = rightLabel;
    style.setTitles(htree[pad],xTitle,ytitle,PadLabel[IndPad[doPl][isIso][pad]-1],RealLabel);
    //style.setTitles(htree[pad],xTitle,ytitle);
  }
  if(Subtract.Contains(".ps")) {
    can.Print(Subtract);
  } else {
    TString plotName = "keys/eps/FinalFit/"; plotName+=sample; plotName+=var; plotName+="_"; 
    plotName += RoundNumber(m2cut[0]*100,0);
    plotName += "_"; plotName += RoundNumber(m2cut[1]*100,0); plotName+=".eps";
    can.SaveAs(plotName);
  }
  //for(int i=0; i<70;i++) if(pdf[i]) pdf[i]->Delete();
  return;
}



