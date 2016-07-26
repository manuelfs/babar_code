//------------------------------------------------------------------------
// File and Version Information:
//      $Id: SidebandmES.cc,v 1.2 2012/08/23 02:22:18 manuelf Exp $
//
// Description:
//      plotVariable - Plots variables using the results of the fit
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      11/08/31 manuelf -- Created from mycode/plotVariable.cc
//------------------------------------------------------------------------

#include "TString.h"
#include "TPad.h"
#include "TCut.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TArrow.h"
#include "TH1F.h"
#include "TF1.h"
#include "THStack.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TTree.h"
#include "DonutUtils/KeysUtils.cc"
#include "DonutUtils/cuts.cc"
#include "DonutUtils/Styles.cc"
#include <fstream>
#include <iostream>

using namespace TMath;
using namespace std;
using std::cout;
using std::endl;
#define nBDT 9

int main(int argc, char *argv[]){
  if (argc < 2 || argc > 15 ) {
    cout << "USAGE: SidebandmES typeSample [Variable=candMES] [nbins=20] [minX=5.2] [maxX=5.26] "<< 
      "[extraCuts=Q2 MVA] [PosLeg=no] [doScale=no] [weightFile=wTotal] [subtract=356789] [emu=both] [textName=Ne2] "<<
      "[ runs=All] [tagName] " << endl;
    return 0;
  }

  TString typeSample = argv[1];
  TString Variable = "candMES"; 
  if(argc>2) Variable = argv[2]; 
  int nbins = 20;
  if(argc>3) {TString temp_s = argv[3]; nbins = temp_s.Atoi();} 
  double minX = 5.2;
  if(argc>4) {TString temp_s = argv[4]; minX = temp_s.Atof();} 
  double maxX = 5.26;
  if(argc>5) {TString temp_s = argv[5]; maxX = temp_s.Atof();} 
  TString extraCuts = "candQ2>4XXcandMvaDl>0.05";
  if (argc>6) extraCuts = argv[6];
  extraCuts.ReplaceAll("XX","&&");
  TString PosLeg = "no";
  if (argc>7) PosLeg = argv[7];
  TString doScale = "no";
  if (argc>8) doScale = argv[8];
  TString weightName = "babar_code/Reweight/wTotal.txt";
  if (argc>9) weightName = argv[9];
  TString subtract = "356789";
  if (argc>10) subtract = argv[10];
  TString emu = "";
  if (argc>11) emu = argv[11];
  TString textName = "FitAll/fits/TextFinalDataNe2x100.txt";
  if (argc>12) textName = argv[12];
  TString runs = "All";
  if (argc>13) runs = argv[13];
  TString tagName = "";
  if (argc>14) tagName = argv[14];

  Styles style4; style4.setPadsStyle(4); 

  int nSamples = 10, isDss = 0; bool isOff = false; 
  if(typeSample.Contains("Off")) {nSamples = 2; isOff = true;}
  if(typeSample.Contains("dss")) isDss = 1;

 
  TCut cuts = CSq2; if(typeSample=="pl") cuts = CSpl; 
  if(typeSample.Contains("MVA")) cuts = MvaAll;
  if(typeSample.Contains("noCut")) cuts = PMiss+M2P;//+!MvaAll; 
  if(typeSample.Contains("Reg2")) cuts = PMiss+M2P+lowQ2+Mva; 
  if(typeSample.Contains("noMVA")) cuts = PMiss+M2P+!Mva;//+!MvaAll; 
  if(typeSample.Contains("Reg251")) cuts = PMiss+M2P+lowQ2+Mva51; 
  if(isDss) cuts = dssMvaAll; if(isDss && typeSample.Contains("noCut")) cuts = dss+!MvaAll;
  if(isDss && typeSample.Contains("MVA")) cuts = dss+dssMva+!MvaAll;
  if(isDss && typeSample.Contains("Dl")) cuts = dss+dssMvaDl+cosT+!MvaAll;
  if(isDss && typeSample.Contains("Comb")) cuts = dss+dssMvaComb+cosT+!MvaAll;
  if(isDss && typeSample.Contains("CosT")) cuts = dss+cosT+!MvaAll;
  if(isDss && typeSample.Contains("Maz")) cuts = dsseeAll;
  if(isDss && typeSample.Contains("51")) cuts = dssMvaAll51;
  if(isOff) cuts = MvaAll;    if(isOff && typeSample.Contains("noCut")) cuts = PMiss+M2P;  
  if(typeSample.Contains("signal")) {
    cuts = MvaAll; 
    if(typeSample.Contains("Mazur")) cuts = PMiss+M2P+ee;
    if(typeSample.Contains("noCut")) cuts = PMiss+M2P; 
    if(typeSample == "signal51") cuts = MvaAll51; 
  }
  if(tagName.Contains("Maz")) cuts = PMiss+M2P; 
  TString rangeCuts = Variable; rangeCuts += ">="; rangeCuts += minX; rangeCuts += "&&"; 
  rangeCuts += Variable;  rangeCuts += "<="; rangeCuts += maxX;
  cuts += extraCuts; cuts += rangeCuts;
  extraCuts.ReplaceAll("candMvaDl","BDT"); extraCuts.ReplaceAll("candQ2","q^{2}"); 
  extraCuts.ReplaceAll("candType","type"); extraCuts.ReplaceAll("&&"," & ");
  extraCuts.ReplaceAll("candMvaDss","BDT_");extraCuts.ReplaceAll("candEExtra","E_{Extra}");
  extraCuts.ReplaceAll("abs(candCosT)","cos(#theta_{T})");
  extraCuts.ReplaceAll("eextrapi0","E_{Extra}");extraCuts.ReplaceAll("ppi0","p_{#pi^{0}}");
  extraCuts.ReplaceAll("mpi0","m_{#pi^{0}}");
  double gEntries=0,uEntries=0,cEntries=0,dEntries=0,nSample[10];
  TTree *gen=0;

  TString tupleFolder = "AWG82/ntuples/small/";
  if(tagName!="") tupleFolder = "AWG82/ntuples/Newsmall/";
  tagName += ".root";
  TString dataName = tupleFolder; dataName += "Data_RunAll"; dataName += tagName;
  TString genName = tupleFolder; genName += "RAll_RunAll"; genName += tagName;
  TString udsName = tupleFolder; udsName += "uds_RunAll"; udsName += tagName;
  TString ccbarName = tupleFolder; ccbarName += "ccbar_RunAll"; ccbarName += tagName;
  TString DpipiName = tupleFolder; DpipiName += "GenDss_RunAll"; DpipiName += tagName;
  if(typeSample.Contains("R24")) genName.ReplaceAll("RAll_","R24_");
  if(typeSample.Contains("R26")) genName.ReplaceAll("RAll_","R26_");
  if(isOff)  dataName = "AWG82/ntuples/small/OffPeak_RunAll.root";
  if(typeSample.Contains("signal")) dataName = genName;
  TTree *data = WeightedTree(dataName, dEntries, "1",0,cuts, runs);
  if(!isOff) gen = WeightedTree(genName, gEntries, weightName,0,cuts,runs);
  TTree *uds = WeightedTree(udsName, uEntries, weightName,-1,cuts);
  TTree *ccbar = WeightedTree(ccbarName, cEntries, weightName,-1,cuts);
  TTree *treeDpipi = WeightedTree(DpipiName, cEntries, weightName,0,cuts);

  if(tagName.Contains("Maz")) runs = "1234";
  double totDssDpipi = 7049000+7036000;
  double totMCB = 0, totuds = 0, totccbar = 0, totdata = 0, totOffdata = 0;
  getNumberB(genName, runs, totMCB, totdata, totuds, totccbar, totOffdata);
  double wuds = totMCB/totuds*2.09/1.05;              // 0.862068
  double wccbar = totMCB/totccbar*1.3/1.05;           // 1.041766
  double wMC = totdata/totMCB;                        // 0.496529
  double wDpipi = 0.01/(totDssDpipi/totMCB);
  if(isOff){ 
    wuds   = totOffdata/totuds*2.09*1000;             // 0.042759
    wccbar = totOffdata/totccbar*1.3*1000;            // 0.051672
    wMC = 1;
  }
  wccbar = wuds;                                      // weightManager converts ccbar into uds
  if(tagName.Contains("Maz")){wuds = 0; wccbar = 0;}

  TString xtitle = Variable, Units = " MeV";
  TLatex *label = new TLatex(); label->SetNDC(kTRUE); label->SetTextFont(style4.nFont);
  TLine line; line.SetLineStyle(2); line.SetLineColor(28); line.SetLineWidth(2);
  TArrow arrow; arrow.SetLineColor(28); arrow.SetFillColor(28); arrow.SetLineWidth(2);
  float nSamInt[4][9], totYields[]={0,0,0}; 
  TH1F *hdata[4], *hStack[4][10];
  THStack *hs[4];
  for (int i=0; i<4; i++) {TString hsName = "hs"; hsName+=i; hs[i] = new THStack(hsName,"");}
  gROOT->SetStyle("Plain");gStyle->SetOptStat(0);
  TCanvas c("dataMC","data Vs MC",style4.CanvasW,style4.CanvasH);
  c.Divide(2,2);
  double leftMargin = 0.15, rightMargin = 0.1;
  TString samCut[2][7] = {{"MCType==0","MCType>6&&MCType<13","MCType>=13&&MCType<=14","(MCType==1||MCType==3)",
			   "(MCType==2||MCType==4)","MCType==5","MCType==6"},
			  {"MCType==0","MCType>0&&MCType<7","MCType>=13&&MCType<=14","(MCType==7||MCType==9)",
			   "(MCType==8||MCType==10)","MCType==11","MCType==12"}};
  int indSam[2][10] = {{0,1,2,3,4,5,6,7,8,9},{0,1,2,3,5,6,7,8,4,9}};
  int IndText[8] = {0,4,1,5,2,6,3,7}, IndScale[3][7]={{6,5,4,1,3,0,2},{6,5,4,3,1,2,0},{4,3,0,2,1,2,1}};
  double Yields[16][9];
  TString buffer;
  for(int i=0; i<16; i++) for(int j=0; j<9; j++)Yields[i][j] = 0;
  TString BDTNames[] = {"030", "050", "075", "100", "120", "150", "200", "250", "300", "400", "500"};
  float BDTCuts[nBDT][4] = {{0.7207, 0.6110, 0.6000, 0.5850},     //0.3x 
			    {0.6300, 0.5650, 0.4740, 0.5250},     //0.5x
			    {0.5413, 0.5358, 0.4320, 0.4870},     //0.75x
			    {0.4580, 0.4800, 0.3600, 0.4100},     //1.0x
			    {0.3810, 0.4546, 0.3080, 0.3659},     //1.2x
			    {0.3065, 0.4133, 0.2530, 0.3171},     //1.5x
			    {0.1709, 0.3440, 0.1787, 0.2650},     //2.0x
			    {0.0990, 0.2804, 0.1275, 0.2240},     //2.5x
			    {0.0638, 0.2350, 0.0889, 0.1940}};    //3.0x


  double genYield[4], datYield[4], subYield[4];
  int isDpipi = 0; if(textName.Contains("Dpipi")) isDpipi = 1;
  for(int iBDT = 0; iBDT < nBDT; iBDT++){
    //for(int iBDT = 3; iBDT < 4; iBDT++){
    cout<<endl<<"Ratio for "<<BDTNames[iBDT]<<"% sample"<<endl;
    textName = "FitAll/fits/TextFinalDataNe2x"; textName += BDTNames[iBDT]; textName += ".txt";
    //textName.ReplaceAll("x300","x250");
    fstream textFile; textFile.open(textName,fstream::in);
    for(int i=0; i<16; i++) for(int j=0; j<9; j++)Yields[i][j] = 0;
    for(int i=0; i<8; i++){
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
    TString mESCorrName = "babar_code/Reweight/wBkg_CombAllmESx"; mESCorrName += BDTNames[iBDT];
    mESCorrName += emu;

    fstream mESCorr; mESCorr.open(mESCorrName,fstream::out);
    for(int i=0; i<4; i++){
      c.cd(i+1);
      TPad *p1 = new TPad("p1","",0,0.26,1,0.98);
      p1->Draw(); p1->cd(); 
      TString totCut = "candType=="; totCut += i+1; 
      if(emu=="e")  totCut += "&&candIsMu==0";
      if(emu=="mu") totCut += "&&candIsMu==1";
      totCut += "&&candMvaDl>"; totCut += BDTCuts[iBDT][i];
      TString hname = "data"; hname += i;
      TString vari = Variable; vari += ">>"; vari += hname;//vari += "*10.5794/10.5394>>";
      //float xbins[] = {0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.15, 1.3, 1.5, 2.4};
      hdata[i] = new TH1F(hname,"",nbins,minX,maxX);
      hdata[i]->Sumw2();
      data->Project(hname,Variable,totCut);
      if(typeSample.Contains("signal")) hdata[i]->Scale(wMC);
      dEntries = hdata[i]->Integral();
      if(i==0) {
	style4.fixYAxis(hdata[i], p1);
	if(style4.nFont==62) leftMargin *= 1.1;
      }
      leftMargin = style4.PadLeftMargin; 
      p1->SetLeftMargin(leftMargin); p1->SetRightMargin(rightMargin); 
      if(style4.isThesis == "yes") {p1->SetTopMargin(0.01);}

      TString ContiCut = "("; ContiCut += totCut; ContiCut += ")*weight";
      hname = "ccbar"; hname += i;
      vari = Variable; vari.ReplaceAll("candM2NF","candM2");
      vari += ">>"; vari += hname;
      hStack[i][0]  = new TH1F(hname,"",nbins,minX,maxX);
      hStack[i][0]->Sumw2();
      ccbar->Project(hname,Variable,ContiCut);
      hStack[i][0]->Scale(wccbar);
      hname = "uds"; hname += i;
      vari = Variable;  vari.ReplaceAll("candM2NF","candM2");
      vari += ">>"; vari += hname;
      hStack[i][1]  = new TH1F(hname,"",nbins,minX,maxX);
      hStack[i][1]->Sumw2();
      nSample[1] = uds->Project(hname,Variable,ContiCut);
      hStack[i][1]->Scale(wuds);
      nSample[0] = hStack[i][0]->Integral(); 
      nSample[1] = hStack[i][1]->Integral(); 
      //     totccbar = hStack[i][0]->Integral(); 
      //     totuds = hStack[i][1]->Integral(); 
      //     hStack[i][0]->Scale(Yields[i+4*isDss][7]/(totccbar+totuds));
      //     hStack[i][1]->Scale(Yields[i+4*isDss][7]/(totccbar+totuds));
      //     nSample[0] = Yields[i+4*isDss][7]/(totccbar+totuds)*totccbar;
      //     nSample[1] = Yields[i+4*isDss][7]/(totccbar+totuds)*totuds;
      hname = "Dpipi"; hname += i;
      vari = Variable;  vari.ReplaceAll("candM2NF","candM2");
      vari += ">>"; vari += hname;
      hStack[i][9]  = new TH1F(hname,"",nbins,minX,maxX);
      hStack[i][9]->Sumw2();
      ContiCut = "("; ContiCut += totCut; ContiCut += "&&MCType==14"; ContiCut += ")*weight";
      nSample[9] = treeDpipi->Project(hname,Variable,ContiCut);
      hStack[i][9]->Scale(wDpipi);
      int iDpipi = 8; if(isDss) iDpipi = 6;
      double scale = Yields[i+4*isDss][iDpipi];
      if(scale>0) scale /= Yields[8+i+4*isDss][iDpipi];
      hStack[i][9]->Scale(scale);
      nSample[9] = hStack[i][9]->Integral(); 
      //cout<<"Scale "<<nSample[9]*wMC<<", true "<<Yields[8+i+4*isDss][iDpipi]<<" and ntuple "<<nSample[9]/wDpipi<<endl;

      int a = 0; if(i>1) a=1;
      gEntries = 0; double gSub=0;
      for(int sam2 = 0; sam2<nSamples; sam2++){
	int sam = indSam[isDss][sam2];
	if(i%2==0){
	  if(sam==7) sam = 8;
	  else if(sam==8) sam = 7;
	}
	if(sam>1&&sam<9){
	  TString genCut = "("; genCut += totCut; genCut += "&&"; 
	  genCut += samCut[a][sam-2]; genCut += ")*weight";
	  hname = "gen"; hname += i+10*sam;
	  vari = Variable; vari += ">>"; vari += hname; 
	  hStack[i][sam]  = new TH1F(hname,"",nbins,minX,maxX);
	  hStack[i][sam]->Sumw2();
	  nSample[sam] = gen->Project(hname,Variable,genCut);
	  int iIndScale = i%2; if(isDss) iIndScale = 2;
	  double scale = Yields[i+4*isDss][IndScale[iIndScale][sam-2]];
	  if(scale>0) scale /= Yields[8+i+4*isDss][IndScale[iIndScale][sam-2]];
	  hStack[i][sam]->Scale(scale);
	  //hStack[i][sam]->Scale(Yields[i+4*isDss][IndScale[iIndScale][sam-2]]/hStack[i][sam]->Integral());
	  //cout<<Yields[i][IndScale[i%2][sam-4]]<<"/"<<Yields[8+i][IndScale[i%2][sam-2]]<<"\t For index "<<IndScale[i%2][sam-2]<<endl;
	}
	TString Sam=""; Sam+=sam;
	double samInt = hStack[i][sam]->Integral();
	nSamInt[i][sam] = samInt;
	nSample[sam] = samInt*wMC;
	//nSample[sam] = samInt;
	if(subtract.Contains(Sam)) {
	  dEntries -= nSample[sam];
	  gSub += nSample[sam];
	} else {
	  gEntries += samInt;
	  if(typeSample.Contains("signal") && sam<2) dEntries += nSample[sam];
	}
      }
      int isCreated = 0, SubisCreated = 0;
      for(int sam2 = 0; sam2<nSamples; sam2++){
	int sam = indSam[isDss][sam2];
	if(i%2==0){
	  if(sam==7) sam = 8;
	  else if(sam==8) sam = 7;
	}
	if(nSample[sam]) {
	  TString Sam=""; Sam+=sam;
	  if(!subtract.Contains(Sam)){
	    if(isCreated==0 && nSample[sam]>0) {
	      isCreated++;
	    }
	  } else {
	    if(SubisCreated==0 && nSample[sam]>0) {
	      SubisCreated++;
	    }
	  }
	}
      }
      if(isCreated==0){
	cout<<"No events pass the cuts in channel "<<i<<endl;
	continue;
      }
      //cout<<"Cont: "<<nSample[0]+nSample[1]<<endl;
      genYield[i] = gEntries; datYield[i] = dEntries; subYield[i] = gSub;
      if(i>=2 && typeSample.Contains("Sum")){
	genYield[i] += genYield[i-2]; genYield[i-2] = genYield[i];
	datYield[i] += datYield[i-2]; datYield[i-2] = datYield[i];
	subYield[i] += subYield[i-2]; subYield[i-2] = subYield[i];
      }
      // Deleting histograms
      hdata[i]->Delete();
      hStack[i][0]->Delete();
      hStack[i][1]->Delete();
      hStack[i][9]->Delete();
      for(int sam2 = 0; sam2<nSamples; sam2++){
	int sam = indSam[isDss][sam2];
	if(sam>1&&sam<9&&hStack[i][sam]) hStack[i][sam]->Delete();
      }
    }
    for(int i=0; i<4; i++){
      gEntries = genYield[i]; dEntries = datYield[i]; double gSub = subYield[i];
      double errRatio = sqrt(gEntries/pow(dEntries,2)+pow(gEntries,2)/pow(dEntries,4)*(dEntries+gSub+gSub*wMC));
      totYields[0] += gEntries; totYields[1] += dEntries; totYields[2] += gSub;
      cout<<"MC: "<<RoundNumber(gEntries*wMC,0)<<", data: "<<RoundNumber(dEntries,1)<<"\t  -  Ratio: "<< 
	RoundNumber(gEntries*wMC,5,dEntries) << " +- "<< RoundNumber(errRatio*wMC,5)<<endl;
      mESCorr<<i<<"   "<<RoundNumber(gEntries*wMC,5,dEntries)<<endl;
    }
    mESCorr.close();
  }

//   double errRatio = sqrt(totYields[0]/pow(totYields[1],2)+
// 			 pow(totYields[0],2)/pow(totYields[1],4)*(totYields[1]+2*totYields[2]));
//   double Ratio = totYields[0]*wMC/totYields[1];
//   cout<<"MC: "<<RoundNumber(totYields[0]*wMC,0)<<", data: "<<RoundNumber(totYields[1],0)<<"\t  -  Ratio: "<< 
//     RoundNumber(Ratio,5) << " +- "<< RoundNumber(errRatio*wMC,5)<<"\t That's "<<
//     RoundNumber(Ratio-1,2,errRatio*wMC)<<" sigma away"<<endl;
//   cout<<"Dss:\t";
//   for(int chan=0;chan<4;chan++) cout<<RoundNumber(nSamInt[chan][4],1)<<", ";
//   cout<<endl<<"Ds:\t";
//   for(int chan=0;chan<4;chan++) cout<<RoundNumber(nSamInt[chan][6],1)<<", ";
//   cout<<endl<<"D:\t";
//   for(int chan=0;chan<4;chan++) cout<<RoundNumber(nSamInt[chan][5],1)<<", ";
//   cout<<endl;

  return 1; 
}

