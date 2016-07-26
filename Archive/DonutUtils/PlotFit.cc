//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: PlotFit.cc,v 1.5 2012/08/23 02:22:18 manuelf Exp $
//
// Description:
//      PlotFit - Plots the results of the final fit
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      10/03/05 manuelf -- Created
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
  TagSample.Remove(0,TagSample.Index("Final")+5); TagSample.ReplaceAll(".txt",""); TagSample.ReplaceAll("IsoData","Data");
  TagSample.ReplaceAll("Ne2","New"); TagSample.ReplaceAll("Hi2","Hig"); 
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
    if(textName.Contains("Hi")) genFile = "AWG82/ntuples/small/FitDataNewx100_RunAll.root";
    if(textName.Contains("Dsstau")) genFile = "AWG82/ntuples/small/FitDataNewx100_RunAll.root";
  }
  if(textName.Contains("Dpip")) {
    genFile = "AWG82/ntuples/small/FitDataNewx100_RunAll.root";
    //typePlot += "Dpipi";
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

  double Yields[16][9], hIntegral[8][2];
  for(int pad=0; pad<8; pad++) for(int mc=0; mc<2; mc++) hIntegral[pad][mc] = 0;
  TString buffer;
  fstream textFile; textFile.open(textName,fstream::in);
  int begChan = 0, endChan = 4, IndText[8] = {0,4,1,5,2,6,3,7};
  if(typePlot.Contains("D0")&&!typePlot.Contains("Dc")) endChan = 2; 
  if(typePlot.Contains("Dc")&&!typePlot.Contains("D0")) begChan = 2;
  for(int i=2*begChan; i<2*endChan; i++){
    int index = IndText[i];
    for(int j=0; j<9; j++){
      textFile>>buffer>>Yields[index+8][j]>>Yields[index][j]>>buffer>>buffer>>buffer;
      //cout<<"True: "<<Yields[index+8][j]<<"\tFit: "<<Yields[index][j]<<"\tChannel "<<index<<endl;
      if(index>3 && j==6){
	textFile>>buffer;
	break;
      }
    }
    //cout<<endl;
  }
  TLine line; line.SetLineWidth(1);
  TCanvas can("can","Final fit");
  if(nPads>1) can.Divide(2,nPads/2);
  int nm2bin = 0, nplbin = 0, Projbins = nbins, ProjbinsPl = nbinsPl;
  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4;
  fstream ResWFile; ResWFile.open("../DonutUtils/ResolutionWidths.txt",fstream::in);
  double ResWidth[4];
  for(int chan=0; chan<4; chan++) ResWFile >> ResWidth[chan];

  TString TagSample = textName;
  TagSample.Remove(0,TagSample.Index("Final")+5); TagSample.ReplaceAll(".txt","");
  if(TagSample.Contains("x")) {
    if(TagSample.Contains("Data")) TagSample.Remove(0,TagSample.Index("Data")+4); 
    if(TagSample.Contains("RAll")) TagSample.Remove(0,TagSample.Index("RAll")+4); 
    TagSample.Remove(TagSample.First('x')+4,TagSample.Sizeof()+1);
    TagSample.ReplaceAll("Ne2","New"); TagSample.ReplaceAll("Hi2","Hig"); 
    TagSample.ReplaceAll("mES","Sig"); TagSample.ReplaceAll("on2","onD"); TagSample.ReplaceAll("o22","onD");
    TagSample.ReplaceAll("Si2","Ori"); TagSample.ReplaceAll("Sid","Ori"); TagSample.ReplaceAll("Or2","Ori"); 
    TagSample.ReplaceAll("Df2","Dfi"); TagSample.ReplaceAll("Ds2","Dsf"); TagSample.ReplaceAll("IsoData","Data");
  }
  if(typePlot.Contains("Dsstau")) TagSample = "Dsstaux100";
  if(typePlot.Contains("Dpipi")) TagSample = "Newx100";
  if(textName.Contains("eData")) TagSample = "eNewx100";
  if(textName.Contains("muData")) TagSample = "muNewx100";
  TH2F *pdf[70];
  for (int i = 1 ; i <= 68 ; i ++) {
    //if(i==37)i=41; 
    if(i==49)i=51; if(i==59)i=61; 
    TString hname = "keys/root/Fit/pdfKeys_"; 
    if(TagSample.Contains("x")) {hname = "keys/root/fit"; hname += TagSample; hname += "/pdfKeys_"; }
    if(TagSample.Contains("Dpip")) {hname = "keys/root/fitNewx100/pdfKeys_"; }
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
  int Indpdf[8][9] = {{61,51,41,17,37,13,9,5,1},   {62,52,42,18,38,14,10,6,2},
		      {63,53,43,19,39,15,11,7,3},  {64,54,44,20,40,16,12,8,4},
		      {65,55,45,29,25,21,33,-1,-1},{66,56,46,30,26,22,34,-1,-1},
		      {67,57,47,31,27,23,35,-1,-1},{68,58,48,32,28,24,36,-1,-1}};
  int IndYield[3][9] = {{7,6,5,4,8,3,1,2,0},{7,6,5,4,8,3,1,2,0},{5,4,3,2,1,0,6,-1,-1}};
  int Colors[7] = {5,28,4,2,3,38,43};
  TString legNames[7] = {"Bkg.","D**l#nu","D*l#nu","Dl#nu","D*#tau#nu","D#tau#nu","D#pi#pil#nu"};
  int IndColor[3][9] = {{0,0,0,1,6,2,3,4,5},{0,0,0,1,6,3,2,5,4},{0,0,0,3,2,1,6,-1,-1}};
  TString channelTitle[4]={"D^{0}","D*^{0}","D^{+}","D*^{+}"};
  //TString PadLabel[] = {"a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)"};
  TString PadLabel[] = {"D0", "D0", "Ds0", "Ds0", "Dp", "Dp", "Dsp", "Dsp"};
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
      htree[pad+4] = new TH1F(hname,hname,nbinsPl,0,2.4);
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
      double scale = Yields[channel+truInd][IndYield[type][i]]*(double)Bins/
	pdf[Indpdf[channel][i]]->Integral()/(double)nbins*(maxx-minx)/(MaxX-MinX);
      //cout<<"Yield "<<Yields[channel+truInd][IndYield[type][i]]<<" for "<< channel+truInd<<" and "<<
      //IndYield[type][i] <<", scale "<<scale<<", Bins "<<Bins<<
      //", Integral "<< pdf[Indpdf[channel][i]]->Integral()<<endl;
      if(doPl) pdfProj[pad][i] = binHisto(pdf[Indpdf[channel][i]],Bins,MinX,MaxX,0,2.4,pdfName,var);
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
      pdfProj[pad][i].SetFillColor(Colors[IndColor[type][i]]);
      pdfProj[pad][i].SetLineColor(Colors[IndColor[type][i]]);
      if(i==1) {pdfProj[pad][i].SetLineColor(1);pdfProj[pad][i].SetLineStyle(2);}
      if(doPl){
	pdfName += "Pl";
	scale = Yields[channel+truInd][IndYield[type][i]]*(double)BinsPl/
	  pdf[Indpdf[channel][i]]->Integral()/(double)nbinsPl;
	pdfProj[pad+4][i] = binHisto(pdf[Indpdf[channel][i]],BinsPl,0,2.4,minY,maxY,pdfName,"candPstarLep");
	pdfProj[pad+4][i].Scale(scale);
	if(Subtract.Contains(SubComp)){
	  for(int bin = 1; bin<=nbinsPl; bin++) pdfProj[pad+4][i].SetBinError(bin,0);
	  htree[pad+4]->Add(&pdfProj[pad+4][i],-1);
	  pdfProj[pad+4][i].Scale(0);
	} else {
	  hSumPl.Add(&pdfProj[pad+4][i]);
	  pdfProj[pad+4][i].Add(&hSumPl,&hSumPl,1,0);
	}
	pdfProj[pad+4][i].SetFillColor(Colors[IndColor[type][i]]);
	pdfProj[pad+4][i].SetLineColor(Colors[IndColor[type][i]]);
	if(i==1) {pdfProj[pad+4][i].SetLineColor(1);pdfProj[pad+4][i].SetLineStyle(2);}
      }
    }

    TString ytitle = "Events/(", xTitle = "m^{2}_{miss} (GeV^{2})", ytitle2 = "Events/(", xTitle2 = "x";//"p*_{l} (GeV)";
    if(var=="candM2") {
      ytitle += RoundNumber((maxx-minx),2,(double)nbins); ytitle += " GeV^{2})";
    } else {
      ytitle += RoundNumber((maxx-minx)*1000,0,(double)nbins); ytitle += " MeV)";
      xTitle = "x";//""p*_{l} (GeV)";
    }
    if(doPl){
      ytitle2 += RoundNumber(2400,0,(double)nbinsPl); ytitle2 += " MeV)";
      style.setMarkers(htree[pad+4], 0.6, 20);
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
      //double nFit = hIntegral[pad][1], nData = hIntegral[pad][0], eData = sqrt(hIntegral[pad][0]); 
      //if(nData<=0) nData = 1;
      //double Ratio = nData/nFit;
      //cout<<"Data: "<<RoundNumber(nData,0)<<", Fit: "<<RoundNumber(nFit,0)<<"\t  -  Ratio: "<< RoundNumber(Ratio,5) 
      //<< " +- "<< RoundNumber(eData/nFit,5)<<"\t That's "<< RoundNumber(Ratio-1,2,eData/nFit)<<" sigma away"<<endl;

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
      if(typePlot.Contains("pull")){
	ytitle2 = "Pull";
	htree[pad+4]->SetMarkerSize(0.75);
	htree[pad+4]->Draw("p");
	line.SetLineStyle(1); line.SetLineColor(4); line.DrawLine(0,0,2.4,0);
	line.SetLineStyle(2); line.SetLineColor(2); line.DrawLine(0,2.,2.4,2.); line.DrawLine(0,-2.,2.4,-2.);
      }
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
      double pM2 = Prob(chi2[pad][0],ndof[pad][0]), pPl = Prob(chi2[pad+4][0],ndof[pad+4][0]);
      double pBoth = Prob(chi2[pad][0]+chi2[pad+4][0], ndof[pad][0]+ndof[pad+4][0]); 
      rightLabel = "#splitline{#chi^{2}: "; rightLabel += RoundNumber(chi2[pad][0],1); rightLabel += "/";
      rightLabel += ndof[pad][0]; rightLabel += " = "; rightLabel += RoundNumber(chi2[pad][0],2,ndof[pad][0]);
      rightLabel += "}{Prob. = "; rightLabel += RoundNumber(pM2*100,2); rightLabel += "%}";
      cout<<channelTitle[pad]<<" p(Chi2):\t M2 "  <<RoundNumber(pM2*10000,2)
      <<"\t Pl "  <<RoundNumber(pPl*10000,2)<<"\t Total "  <<RoundNumber(pBoth*10000,2);
      if(pad>1){
	pM2 = Prob(chi2[pad-2][0]+chi2[pad][0],ndof[pad-2][0]+ndof[pad][0]); 
	pPl = Prob(chi2[pad-2+4][0]+chi2[pad+4][0],ndof[pad-2+4][0]+ndof[pad+4][0]);
	pBoth = Prob(chi2[pad-2][0]+chi2[pad][0]+chi2[pad-2+4][0]+chi2[pad+4][0],
		     ndof[pad-2][0]+ndof[pad][0]+ndof[pad-2+4][0]+ndof[pad+4][0]); 
	cout<<"\t M2Tot "  <<RoundNumber(pM2*10000,2)<<"\t PlTot "  <<RoundNumber(pPl*10000,2)
	<<"\t  Combined "  <<RoundNumber(pBoth*10000,2);
      }
      cout<<endl;
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
      if((i>2 || i==1) && (Indpdf[channel][i]<33||Indpdf[channel][i]>40||
			   typePlot.Contains("Dpipi")))leg[pad]->AddEntry(&pdfProj[pad][i],legNames[IndColor[type][i]]);
    }
    //TString pos = "right"; if(!isDss && var=="candPstarLep") pos="left";
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
    if(doPl){
      can.cd(IndPad[doPl][isIso][pad]+1);
      if(typePlot.Contains("histo")){
	rightLabel = "#splitline{#chi^{2}: "; rightLabel += RoundNumber(chi2[pad+4][0],1); rightLabel += "/";
	rightLabel += ndof[pad+4][0]; rightLabel += " = "; rightLabel += RoundNumber(chi2[pad+4][0],2,ndof[pad+4][0]);
	rightLabel += "}{Prob. = "; rightLabel += RoundNumber(Prob(chi2[pad+4][0],ndof[pad+4][0])*100,2); rightLabel += "%}";
      }
      if(typePlot.Contains("pull")) rightLabel = "";				     
      for(int i=maxProj; i>=0; i--) {
	TString SubComp = ""; SubComp += i;
	if(!typePlot.Contains("pull") && !Subtract.Contains(SubComp)) pdfProj[pad+4][i].Draw(optionPlot);
      }
      //TString pos = "right"; if(!isDss && var=="candPstarLep") pos="left";
      if(!typePlot.Contains("pull")){
	hSum.Draw("same");
	htree[pad+4]->Draw("e1 sameaxis");
	htree[pad+4]->Draw("e0 e1 same");
      }
      TString RealLabel = rightLabel; RealLabel += ",  ";
      if(rightLabel=="") RealLabel = channelTitle[pad];
      else RealLabel += channelTitle[pad];
      RealLabel = rightLabel; //For thesis
      if(style.isThesis=="yes") {
	if(isDss && !typePlot.Contains("histo")){
	  TString dssLabel = "x"; dssLabel += pad;
	  style.setTitles(htree[pad+4],xTitle2,ytitle2,RealLabel,dssLabel);
	}else
	  style.setTitles(htree[pad+4],xTitle2,ytitle2,PadLabel[IndPad[doPl][isIso][pad]],RealLabel);
      } else style.setTitles(htree[pad+4],xTitle2,ytitle2,channelTitle[pad],rightLabel);
    }
    can.cd(IndPad[doPl][isIso][pad]);
    rightLabel = safeRLabel;
    if(didLeg || typePlot.Contains("pull") ) rightLabel = "";	
    if(didLeg) channelTitle[pad] += "                     ";
    TString RealLabel = rightLabel; RealLabel += ",  ";
    if(rightLabel=="" || doPl) RealLabel = channelTitle[pad];
    else RealLabel += channelTitle[pad];
    RealLabel = rightLabel; //For thesis
    if(doPl && !typePlot.Contains("histo")) RealLabel = "";
    if(style.isThesis=="yes") style.setTitles(htree[pad],xTitle,ytitle,PadLabel[IndPad[doPl][isIso][pad]-1],RealLabel);
    else style.setTitles(htree[pad],xTitle,ytitle,channelTitle[pad],safeRLabel);
    //style.setTitles(htree[pad],xTitle,ytitle);
  }
  if(Subtract.Contains(".ps")) {
    can.Print(Subtract);
  } else {
    TString plotName = "keys/eps/FinalFit/"; plotName+=sample; plotName+=var; plotName+="_"; 
    plotName+=typePlot; plotName+=".eps";
    can.SaveAs(plotName);
  }
//   double totChi2[2] = {0,0}; int totndof[2] = {0,0};
//   for(int i=begChan; i<endChan;i++) {
//     leg[i]->Delete();
//     htree[i]->Delete();
//     cout<<"Total: "<<RoundNumber(chi2[i][0],2)<<" chi2 for "<<ndof[i][0]<<" dof -> "<<RoundNumber(chi2[i][0],2,ndof[i][0])
// 	<<"\t\tSubset: "<<RoundNumber(chi2[i][1],2)<<" chi2 for "<<ndof[i][1]<<" dof -> "
// 	<<RoundNumber(chi2[i][1],2,ndof[i][1]);
//     if(typePlot.Contains("m1")) cout<<"/"<<RoundNumber(chi2[i][1],2,ndof[i][1]+1);
//     cout<<endl;
//     totChi2[0] += chi2[i][0];totChi2[1] += chi2[i][1];
//     totndof[0] += ndof[i][0]; totndof[1] += ndof[i][1]; 
//   }
//   cout<<endl<<"Total: "<<RoundNumber(totChi2[0],2)<<" chi2 for "<<totndof[0]<<" dof -> "<<RoundNumber(totChi2[0],2,totndof[0])
//       <<"\t\tSubset: "<<RoundNumber(totChi2[1],2)<<" chi2 for "<<totndof[1]<<" dof: "
//       <<RoundNumber(totChi2[1],2,totndof[1]);
//     if(typePlot.Contains("m1")) cout<<"/"<<RoundNumber(totChi2[1],2,totndof[1]+4);
//     cout<<endl;
  //for(int i=0; i<70;i++) if(pdf[i]) pdf[i]->Delete();
  return;
}



