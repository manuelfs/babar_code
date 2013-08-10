//------------------------------------------------------------------------
// File and Version Information:
//      $Id: OptiBDT.cc,v 1.3 2012/03/04 00:33:03 manuelf Exp $
//
// Description:
//      OptiBDT - Plots purity and significance for different BDT cuts
//                It can also just print the yields for the BDT cuts
//                by including "Yield" in typeSample
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      10/09/26 manuelf -- Created from babar_code/CutOpti/*
//------------------------------------------------------------------------

#include "TString.h"
#include "TPad.h"
#include "TCut.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TTree.h"
#include "DonutUtils/KeysUtils.cc"
#include "DonutUtils/cuts.cc"
#include "DonutUtils/Styles.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

int main(int argc, char *argv[]){
  if (argc < 5 || argc > 7 ) {
    cout << "USAGE: OptiBDT typeSample nbins(candType) minX(Dl) maxX(Comb) "<< 
      "[fileTag=RAll] [extraCuts=1]" << endl;
    return 0;
  }

  Styles style2; style2.setPadsStyle(2); style2.applyStyle();
  TString typeSample = argv[1];
  TString temp_s = argv[2]; int nbins = temp_s.Atoi(); 
  temp_s = argv[3]; double minX = temp_s.Atof(); 
  temp_s = argv[4]; double maxX = temp_s.Atof(); 
  TString fileTag = "RAll";
  if (argc>5) fileTag = argv[5];
  TString extraCuts = "";
  if (argc>6) extraCuts = argv[6];
  extraCuts.ReplaceAll("XX","&&");
  TString weightName = "babar_code/Reweight/wTotal.txt";

  int isOld = 0, isDss = 0, is2D = 0;
  double gEntries=0, uEntries=0, cEntries=0, Yields[6];
  TString tupleFolder = "AWG82/ntuples/small/";
  if(fileTag.Contains("old")) {fileTag.ReplaceAll("old",""); tupleFolder = "AWG82/ntuples/Oldsmall/"; isOld=1;}
  TString genName   = tupleFolder; genName   += fileTag; genName += "_RunAll.root";
  TString udsName   = tupleFolder; udsName   += "uds_RunAll.root"; 
  TString ccbarName = tupleFolder; ccbarName += "ccbar_RunAll.root"; 

  TCut cuts = PMiss+Q2+M2P+M2Signal+MEScut;
  if(isDss) {
    cuts = dss+cosT+Q2;
    if(isOld==0) cuts += Mpi0;
  }
  if(extraCuts != "") cuts += extraCuts;
  TTree *gen = WeightedTree(genName, gEntries, weightName,0, cuts);
  TTree *uds = WeightedTree(udsName, uEntries, weightName,-1, cuts);
  TTree *ccbar = WeightedTree(ccbarName, cEntries, weightName,-1, cuts);

  double totMCB = 0, totuds = 0, totccbar = 0, totdata = 0, totOffdata = 0;
  getNumberB(genName, "All", totMCB, totdata, totuds, totccbar, totOffdata);
  double wuds = totMCB/totuds*2.09/1.05;                      // 0.862068
  double wccbar = totMCB/totccbar*1.3/1.05;                   // 1.041766
  double wMC = totdata/totMCB;                                // 0.496529

  if(typeSample.Contains("NoSca")) wMC = 1;

  if(typeSample.Contains("dss")) isDss = 1;
  if(!typeSample.Contains("Comb") && !typeSample.Contains("Dl")) is2D = 1;
  TString sigCuts[] = {"MCType==5", "MCType==6", "MCType==11", "MCType==12"};
  TString NosigCuts[] = {"MCType!=6", "MCType!=5", "MCType!=12", "MCType!=11"};
  TString samCut[2][5] = {{"(MCType==5||MCType==6)","MCType>0&&MCType<5","MCType>=13&&MCType<=14","MCType>6&&MCType<13","MCType==0"},
			  {"(MCType==11||MCType==12)","MCType>6&&MCType<11","MCType>=13&&MCType<=14","MCType>0&&MCType<7","MCType==0"}};
  int indSam[2][6] = {{0,1,2,3,4,5},{2,1,0,3,4,5}};
  TString legTag[] = {"Signal", "Norm", "D**", "Cross", "Comb", "Cont"};
  TString hInt_s = "hIntegral";
  TH1F *hIntegral = new TH1F(hInt_s,"",5,-40,40);;
  if(typeSample.Contains("Yield")) {
    TString sCuts = "candType=="; sCuts += nbins; sCuts += "&&candMva"; 
    if(isDss) sCuts += "Dss";
    sCuts += "Dl>"; sCuts += minX;  sCuts += "&&candMva"; 
    if(isDss) sCuts += "Dss";
    sCuts += "Comb>"; sCuts += maxX; 
    for(int sam = 0; sam < 5; sam++){
      TString SamCut = sCuts; SamCut += "&&"; 
      if(!isDss && sam==0 && typeSample.Contains("sep"))SamCut += sigCuts[nbins-1];
      else SamCut += samCut[nbins>2][indSam[isDss][sam]];
      TCut totCut = "1"; totCut += SamCut; totCut *= "weight";
      gen->Project(hInt_s,"candM2",totCut);
      Yields[sam] = hIntegral->Integral()*wMC;
    }
    TCut totCut = "1"; totCut += sCuts; totCut *= "weight";
    uds->Project(hInt_s,"candM2",totCut);
    Yields[5] = hIntegral->Integral()*wuds*wMC;
    ccbar->Project(hInt_s,"candM2",totCut);
    Yields[5] += hIntegral->Integral()*wccbar*wMC;
    double total = 0;
    for(int sam = 0; sam < 6; sam++) total += Yields[sam];
    cout<<"S/(S+B): "<<RoundNumber(Yields[0],2,total)<<", S/sqrt(S+B): "<<RoundNumber(Yields[0],2,sqrt(total))<<"\t";
    if(isDss) Yields[1] += Yields[2];
    for(int sam = 0; sam < 6; sam++) {
      if(isDss && sam == 2) continue;
      cout<<legTag[indSam[isDss][sam]]<<": "<<RoundNumber(Yields[sam],0)<<". ";
    }
    cout<<endl;
    return 1;
  }

  cout<<"Using "<<genName<<endl;
  TH1F *hBDT1D[4];
  TH2F *hBDT2D[4];
  for(int i=0; i<4; i++){
    TString hname = "BDT"; hname += i;
    if(is2D) hBDT2D[i] = new TH2F(hname,"", nbins, minX, maxX, nbins, minX, maxX);
    else hBDT1D[i] = new TH1F(hname,"", nbins, minX, maxX);
  }
  TString Dnames[] = {"D^{0}","D*^{0}","D^{+}","D*^{+}"};
  TString EexCut[] = {"0.2","0.2","0.15","0.3"};
  TLatex *label = new TLatex();
  label->SetNDC(kTRUE);
  label->SetTextSize(0.05);
  int nybins = 1;
  double binWidth = (maxX-minX)/((double)nbins), maxi=-99, OldSigni[4];
  if(is2D) nybins = nbins;
  double legW = 0.3, legH = 0.19;
  double legX = style2.PadLeftMargin+legW+0.03, legY = style2.PadBottomMargin+legH+0.085;
  TLegend *leg[2];
  leg[0] = new TLegend(legX-legW, legY-legH, legX, legY);
  leg[1] = new TLegend(legX-legW, legY-legH, legX, legY);
  TCanvas c("cBDT","Optimization BDT");
  c.Divide(2,1);
  for(int chan = 0; chan < 4; chan++){
    TPad *cPad = (TPad *)c.cd(chan/2+1);
    cout<<"Calculating "<<Dnames[chan]<<endl;
    
    TString sCuts = "candType=="; sCuts+=chan+1; sCuts+="&&candEExtra<"; sCuts+=EexCut[chan]; 
    TCut totCut = "1"; totCut += sCuts; 
    if(!isDss && typeSample.Contains("sep")) totCut += NosigCuts[chan];
    totCut *= "weight";
    TCut sigCut = "1"; sigCut += sCuts; 
    if(!isDss && typeSample.Contains("sep")) sigCut += sigCuts[chan];
    else sigCut += samCut[chan>1][indSam[isDss][0]];
    sigCut *= "weight";
    TString vari = "candM2>>"; vari += hInt_s;
    gen->Draw(vari,totCut);
    double total = hIntegral->Integral()*wMC;
    uds->Draw(vari,totCut);
    total += hIntegral->Integral()*wuds*wMC;
    ccbar->Draw(vari,totCut);
    total += hIntegral->Integral()*wuds*wMC;
    gen->Draw(vari,sigCut);
    double signal = hIntegral->Integral()*wMC;
    double signi = 0;
    if(total) OldSigni[chan] = sqrt(signal*signal/total);
    cout<<"Old significance was "<<RoundNumber(OldSigni[chan],2)<<endl;
    double hiSig = 0, hiDl=-99, hiComb=-99;
    for(int xbin = 1; xbin < nbins+1; xbin++){
      if(xbin%10==0)cout<<"Bin "<<xbin<<" of "<<nbins<<endl;
      for(int ybin = 1; ybin < nybins+1; ybin++){
	double Dlcut = ((double)xbin-0.5)*binWidth+minX;
	double Combcut = ((double)ybin-0.5)*binWidth+minX;
	if(typeSample.Contains("Dl")) Combcut=-99.;
	if(typeSample.Contains("Comb")) {Combcut=Dlcut; Dlcut=-99.;}
	sCuts = "candType=="; sCuts += chan+1; sCuts += "&&candMva"; 
	if(isDss) sCuts += "Dss";
	sCuts += "Dl>"; sCuts += Dlcut;  sCuts += "&&candMva"; 
	if(isDss) sCuts += "Dss";
	sCuts += "Comb>"; sCuts += Combcut; 
	totCut = "1"; totCut += sCuts; 
	if(!isDss && typeSample.Contains("sep")) totCut += NosigCuts[chan];
	totCut *= "weight";
	sigCut = "1"; sigCut += sCuts; 
	if(!isDss && typeSample.Contains("sep")) sigCut += sigCuts[chan];
	else sigCut += samCut[chan>1][indSam[isDss][0]];
	sigCut *= "weight";
	gen->Draw(vari,totCut);
	total = hIntegral->Integral()*wMC;
	uds->Draw(vari,totCut);
	total += hIntegral->Integral()*wuds*wMC;
	ccbar->Draw(vari,totCut);
	total += hIntegral->Integral()*wuds*wMC;
	gen->Draw(vari,sigCut);
	signal = hIntegral->Integral()*wMC;
	signi = 0;
	if(total) signi = sqrt(signal*signal/total);
	if(signi==0) signi = -0.5;
	if(is2D) hBDT2D[chan]->SetBinContent(xbin,ybin,signi);
	else hBDT1D[chan]->SetBinContent(xbin,signi);
	if(signi>hiSig){
	  hiSig=signi;
	  hiDl = Dlcut; hiComb = Combcut;
	}
      }
    }
    TString hTitle = Dnames[chan]; hTitle += " - S/#surd(S+B)_{max} is "; hTitle += RoundNumber(hiSig,2);
    TLine line;  line.SetLineStyle(2); //line.SetLineWidth(2);
    if(is2D){
      hTitle += ": BDT_{Dl} > "; hTitle += RoundNumber(hiDl,2); hTitle += " & BDT_{Comb} > ";
      hTitle += RoundNumber(hiComb,2);
      hBDT2D[chan]->SetXTitle("Combinatoric BDT cut");
      hBDT2D[chan]->SetYTitle("Dlnu BDT cut");
      hBDT2D[chan]->GetYaxis()->SetLabelSize(0.057);
      hBDT2D[chan]->GetXaxis()->SetLabelSize(0.055);
      hBDT2D[chan]->GetYaxis()->SetNdivisions(6+100*2);
      hBDT2D[chan]->GetXaxis()->SetNdivisions(7+100*2);
      hBDT2D[chan]->GetXaxis()->SetTitleSize(0.05);
      hBDT2D[chan]->GetYaxis()->SetTitleSize(0.05);
      hBDT2D[chan]->Draw("cont4z");
    } else {
      TString xtitle = "Signal BDT cut";
      if(!typeSample.Contains("Dl")) xtitle = "Combinatoric BDT cut";
      Dnames[chan] += " channel";
      hBDT1D[chan]->SetLineWidth(2);
      if(chan%2 == 0) {
	hBDT1D[chan]->SetLineColor(2);
	hBDT1D[chan]->SetMinimum(0);
	style2.setTitles(hBDT1D[chan],xtitle, "Significance");
	maxi = hBDT1D[chan]->GetMaximum();
	leg[chan/2]->AddEntry(hBDT1D[chan],Dnames[chan]);
      } else {
	if(maxi<hBDT1D[chan]->GetMaximum()) maxi = hBDT1D[chan]->GetMaximum();
	hBDT1D[chan-1]->SetMaximum(maxi*1.1);
	hBDT1D[chan-1]->Draw("hist c");
	style2.fixYAxis(hBDT1D[chan-1],cPad);
	line.SetLineColor(2); line.DrawLine(minX,OldSigni[chan-1],maxX,OldSigni[chan-1]);
	hBDT1D[chan]->SetLineColor(4);
	hBDT1D[chan]->Draw("hist c same");
	line.SetLineColor(4); line.DrawLine(minX,OldSigni[chan],maxX,OldSigni[chan]);
	leg[chan/2]->AddEntry(hBDT1D[chan],Dnames[chan]);
	leg[chan/2]->SetTextSize(style2.LabelSize); leg[chan/2]->SetFillColor(0); 
	leg[chan/2]->SetTextFont(style2.nFont);  leg[chan/2]->SetBorderSize(0);
	leg[chan/2]->Draw();
      }
    }
    //label->DrawLatex(0.09,0.93,hTitle); 
    cout<<hTitle<<endl;
  }

  TString fileName = "keys/eps/OptiBDT/"; fileName += typeSample; fileName += "_"; 
  if(isOld) fileName += "Old";
  fileName += fileTag; fileName += ".eps";
  if(style2.isThesis == "yes")  fileName = "public_html/EvtSel_OptiSigni.eps"; 
  c.SaveAs(fileName);

  hIntegral->Delete(); label->Delete();
  for(int i=0; i<4; i++){
    if(is2D) hBDT2D[i]->Delete();
    else hBDT1D[i]->Delete();
  }
  return 1; 
}

