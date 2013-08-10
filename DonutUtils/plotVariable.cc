//------------------------------------------------------------------------
// File and Version Information:
//      $Id: plotVariable.cc,v 1.2 2012/08/23 02:22:18 manuelf Exp $
//
// Description:
//      plotVariable - Plots variables using the results of the fit
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      11/06/10 manuelf -- Created from mycode/CSample.cc
//      12/03/03 manuelf -- Adapted it to use the exact files employed in
//                          fit to avoid weight problems
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

int main(int argc, char *argv[]){
  if (argc < 2 || argc > 15 ) {
    cout << "USAGE: plotVariable typeSample [Variable=candM2] [nbins=20] [minX=-0.5] [maxX=1.5] "<< 
      "[extraCuts=1] [PosLeg=r] [maxY=1.1] [weightFile=wTotal] [subtract] [textName=Iter2] "<<
      "[emu=both] [ runs=All] [tagName=Newx100] " << endl;
    return 0;
  }

  TString typeSample = argv[1];
  TString Variable = "candM2"; 
  if(argc>2) Variable = argv[2]; 
  int nbins = 20;
  if(argc>3) {TString temp_s = argv[3]; nbins = temp_s.Atoi();} 
  double minX = -0.5;
  if(argc>4) {TString temp_s = argv[4]; minX = temp_s.Atof();} 
  double maxX = 1.5;
  if(argc>5) {TString temp_s = argv[5]; maxX = temp_s.Atof();} 
  TString extraCuts = "1";
  if (argc>6) extraCuts = argv[6];
  extraCuts.ReplaceAll("XX","&&");
  TString PosLeg = "r";
  if (argc>7) PosLeg = argv[7];
  double maxY = 1.1;
  if (argc>8) {TString temp_s = argv[8]; maxY = temp_s.Atof();} 
  TString weightName = "babar_code/Reweight/wTotal.txt";
  if (argc>9) weightName = argv[9];
  TString subtract = "";
  if (argc>10) subtract = argv[10];
  TString textName = "FitAll/fits/TextFinalDataNe2x100.txt";
  if (argc>11) textName = argv[11];
  TString emu = "both";
  if (argc>12) emu = argv[12];
  TString runs = "All";
  if (argc>13) runs = argv[13];
  TString tagName = "Newx100";
  if (argc>14) tagName = argv[14];

  TString doScale = "no";
  int nSamples = 10, isDss = 0; 
  if(typeSample.Contains("dss")) isDss = 1;
  if(extraCuts=="") extraCuts = "1";

  double gEntries=0,dEntries=0,nSample[10];
  TString tupleFolder = "AWG82/ntuples/small/Fit";
  TString dataName = tupleFolder; dataName += "Data"; dataName += tagName; dataName += "_RunAll.root";
  TString genName = tupleFolder; genName += "RAll"; genName += tagName; genName += "_RunAll.root";
  if(typeSample.Contains("signal")) dataName = genName;
  TChain *data = new TChain("ntp1"); data->Add(dataName);
  TChain *gen = new TChain("ntp1"); gen->Add(genName);


  fstream ResWFile; ResWFile.open("../DonutUtils/ResolutionWidths.txt",fstream::in);
  double ResWidth[4];
  for(int chan=0; chan<4; chan++) ResWFile >> ResWidth[chan];
  TString legTag[] = {"Cont (", "Cont (", "Bkg (", "Cross (", "D**l#nu (", "Dl#nu (", "D*l#nu (", "D#tau#nu (", "D*#tau#nu (", "D#pi#pil#nu ("};
  TString titles[] = {"D^{0}: MC/data = ","D*^{0}: MC/data = ","D^{+}: MC/data = ","D*^{+}: MC/data = "};
  TString xtitle = Variable, Units = " MeV";
  if(Variable=="candM2"||Variable == "mm2pi0"||Variable=="candM2NF"||Variable=="candQ2") Units = " GeV^{2}";
  if(Variable=="candCosT"||Variable.Contains("candTagNeutMult")||Variable.Contains("candTagChargedMult")
     ||Variable.Contains("Mva")||Variable.Contains("candW")) Units = "";
  if(Variable=="candM2" || Variable == "mm2pi0") xtitle = "m^{2}_{miss} (GeV^{2})"; 
  if(Variable=="candPstarLep") xtitle = "p*_{l} (GeV)"; 
  if(Variable=="candPstarD") xtitle = "p*_{D} (GeV)"; 
  if(Variable=="candW") xtitle = "w"; 
  if(Variable=="candM2NF") xtitle = "m^{2}_{miss,no fit} (GeV^{2})"; 
  if(Variable=="candQ2") xtitle = "q^{2} (GeV^{2})"; if(Variable=="candMES") xtitle = "m_{ES} (GeV)"; 
  if(Variable=="candEExtra" || Variable=="eextrapi0") xtitle = "E_{Extra} (GeV)"; 
  if(Variable=="candMvaDl") xtitle = "Signal BDT"; if(Variable=="candMvaDssDl") xtitle = "Semileptonic BDT";
  if(Variable=="candMvaDssComb") xtitle = "Comb. BDT"; if(Variable=="candDmass") xtitle = "";
  if(Variable=="candBTagDmass") xtitle = "";
  if(Variable=="candBTagDeltam") xtitle = "";if(Variable=="candDeltam") xtitle = "";
  if(Variable=="candDeltaE") xtitle = "#DeltaE (GeV)";if(Variable=="candCosT") xtitle = "cos(#theta_{T})";
  if(Variable.Contains("candTagChargedMult")) xtitle = "Tag charged mult.";
  if(Variable.Contains("candTagNeutMult")) xtitle = "Neutral mult."; if(Variable=="mpi0") xtitle = "m(#pi^{0}) (GeV)";
  if(Variable=="dmpi0") xtitle = "m(D#pi^{0})-m(D) (GeV)";
  if(Variable=="ppi0") xtitle = "p(#pi^{0}) (GeV)"; if(Variable=="e1pi0") xtitle = "max(E_{#gamma,#pi^{0}}) (GeV)";
  if(Variable=="pmisspi0" || Variable=="candPMiss") xtitle = "p_{miss} (GeV)";
  if(isDss){
    if(Variable=="candM2") Variable = "mm2pi0";
    if(Variable=="candMvaDl") Variable = "candMvaDssDl";
    if(Variable=="candMvaComb") Variable = "candMvaDssComb";
  }


  Styles style4; style4.setPadsStyle(4);
  if(typeSample.Contains("Iso")) {
    style4.CanvasH /= 2;
    titles[2] = "D: MC/data = ";
    titles[3] = "D*: MC/data = ";
  }

  int doRoot = 0;
  if(typeSample.Contains("Root")) doRoot = 1;
  TH1F hYields("hYields","",40,0,40);
  tagName.ReplaceAll(".root",""); 
  TString histosName = "keys/eps/CSample/root/"; histosName += typeSample;  histosName += tagName; 
  histosName += "_"; histosName += Variable;  histosName += "_"; histosName += runs;  histosName += "_"; 
  histosName += emu;  histosName += ".root";
  TFile *fHistos=0; if(doRoot) fHistos = new TFile(histosName,"RECREATE");

  TLatex *label = new TLatex(); label->SetNDC(kTRUE); label->SetTextFont(style4.nFont);
  TLine line; line.SetLineStyle(2); line.SetLineColor(28); line.SetLineWidth(2);
  TArrow arrow; arrow.SetLineColor(28); arrow.SetFillColor(28); arrow.SetLineWidth(2);
  float nTau[4], nSamInt[4][9], totYields[]={0,0,0}; 
  TH1F *hdata[4], *hdata2[4], *hDiff[4], *hStack[4][10], *hMCsum[4], *hMCsub[4]={0,0,0,0}, *hTau[4];
  THStack *hs[4];
  for (int i=0; i<4; i++) {TString hsName = "hs"; hsName+=i; hs[i] = new THStack(hsName,"");}
  gROOT->SetStyle("Plain");gStyle->SetOptStat(0);
  TCanvas c("dataMC","data Vs MC",style4.CanvasW,style4.CanvasH);
  if(!typeSample.Contains("Iso")) c.Divide(2,2);
  else c.Divide(2,1);
  double leftMargin = 0.15, rightMargin = 0.1;
  if(style4.isThesis == "yes") rightMargin = style4.PadRightMargin;
  int colors[] = {5,5,5,5,28,2,9,38,3,43,1};
  if(isDss){colors[7] = 9; colors[8] = 9;}
  TString samCut[2][7] = {{"MCType==0","MCType>6&&MCType<13","MCType>=13&&MCType<=14","(MCType==1||MCType==3)",
			   "(MCType==2||MCType==4)","MCType==5","MCType==6"},
			  {"MCType==0","MCType>0&&MCType<7","MCType>=13&&MCType<=14","(MCType==7||MCType==9)",
			   "(MCType==8||MCType==10)","MCType==11","MCType==12"}};
  int indSam[2][10] = {{0,1,2,3,4,5,6,7,8,9},{0,1,2,3,5,6,7,8,4,9}};
  int IndText[8] = {0,4,1,5,2,6,3,7}, IndScale[3][7]={{6,5,4,1,3,0,2},{6,5,4,3,1,2,0},{4,3,0,2,1,2,1}};
  double Yields[16][9];
  TString buffer;
  if(typeSample.Contains("Dpipi")) textName = "FitAll/fits/TextFinalDpipiDataDss.txt";
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
  int isDpipi = 0, pad; if(textName.Contains("Dpipi")) isDpipi = 1;
  for(int i=0; i<4; i++){
    pad = i+1; if(typeSample.Contains("Iso") && i>1) pad -= 2;
    c.cd(pad);
    TPad *p1 = new TPad("p1","",0,0.26,1,0.98);
    p1->Draw(); p1->cd(); 
    TString totCut = extraCuts; totCut += "&&candType=="; totCut += i+1; 
    if(emu=="e")  totCut += "&&candIsMu==0";
    if(emu=="mu") totCut += "&&candIsMu==1";
    TString hname = "Data"; hname += i;
    TString vari = Variable; vari += ">>"; vari += hname;
    hdata[i] = new TH1F(hname,"",nbins,minX,maxX);
    hdata[i]->Sumw2();
    data->Project(hname,Variable,totCut);
    dEntries = hdata[i]->Integral();
    hYields.SetBinContent(i*10+9, dEntries);
    if(i==0) {
      style4.fixYAxis(hdata[i], p1);
      if(style4.nFont==62) leftMargin *= 1.1;
    }
    leftMargin = style4.PadLeftMargin; 
    p1->SetLeftMargin(leftMargin); p1->SetRightMargin(rightMargin); 
    if(style4.isThesis == "yes") {p1->SetTopMargin(0.01);}

    TString ContiCut = "(MCType==-1&&"; ContiCut += totCut; ContiCut += ")*weight";
    hname = "MC"; hname += (i+1)*10;
    vari = Variable; 
    vari += ">>"; vari += hname;
    hStack[i][0]  = new TH1F(hname,"",nbins,minX,maxX);
    hStack[i][0]->Sumw2();
    gen->Project(hname,Variable,ContiCut);
    ContiCut = "(MCType==-2&&"; ContiCut += totCut; ContiCut += ")*weight";
    hname = "MC"; hname += (i+1)*10+1;
    vari = Variable;  
    vari += ">>"; vari += hname;
    hStack[i][1]  = new TH1F(hname,"",nbins,minX,maxX);
    hStack[i][1]->Sumw2();
    nSample[1] = gen->Project(hname,Variable,ContiCut);
    nSample[0] = hStack[i][0]->Integral(); 
    nSample[1] = hStack[i][1]->Integral(); 
    hname = "Dpipi"; hname += i;
    vari = Variable;  vari.ReplaceAll("candM2NF","candM2");
    vari += ">>"; vari += hname;
    hStack[i][9]  = new TH1F(hname,"",nbins,minX,maxX);
    hStack[i][9]->Sumw2();
    ContiCut = "("; ContiCut += totCut; ContiCut += "&&MCType==16"; ContiCut += ")*weight";
    nSample[9] = gen->Project(hname,Variable,ContiCut);
    int iDpipi = 8; if(isDss) iDpipi = 6;
    double scale = Yields[i+4*isDss][iDpipi];
    if(scale>0) scale /= Yields[8+i+4*isDss][iDpipi];
    hStack[i][9]->Scale(scale);
    nSample[9] = hStack[i][9]->Integral(); 

    int a = 0; if(i>1) a=1;
    TString TauCut = "("; TauCut += totCut; TauCut += "&&MCType=="; 
    TauCut += i+5+a*4; TauCut += ")*weight";
    TString hTauName = "hTau"; hTauName += i;
    hTau[i] = new TH1F(hTauName,"",nbins,minX,maxX);
    vari = Variable; vari += ">>"; vari += hTauName;
    gen->Project(hTauName,Variable,TauCut);
    nTau[i] = hTau[i]->Integral();
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
	hname = "MC"; hname += (i+1)*10+sam;
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
      if(typeSample.Contains("Conv")&&(i%2==0&&sam==5||i%2==1&&sam==6 || sam==4)){
	hname = "Sample"; hname += 10*i+sam;
	hStack[i][sam] = GaussConv(hStack[i][sam],ResWidth[i],hname,8,nbins,minX,maxX);
      }
      hStack[i][sam]->SetFillColor(colors[sam]);
      hStack[i][sam]->SetLineColor(colors[sam]);
      if(sam==2) {
	hStack[i][sam]->SetLineColor(colors[10]);
	hStack[i][sam]->SetLineStyle(2);
      }
      double samInt = hStack[i][sam]->Integral();
      nSamInt[i][sam] = samInt;
      nSample[sam] = samInt;
      hYields.SetBinContent(i*10+sam, nSample[sam]);
      //cout<<"sam "<<sam<<": "<<nSamInt[i][sam]<<" - "<<nSample[sam]<<endl;
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
	  if(doScale=="yes") hStack[i][sam]->Scale(dEntries/gEntries);
	  if(typeSample.Contains("signal") && sam<2) hdata[i]->Add(hStack[i][sam]);
	  if(nSample[sam]) hs[i]->Add(hStack[i][sam],"hist");
	  if(isCreated>0) hMCsum[i]->Add(hStack[i][sam]);
	  if(isCreated==0 && nSample[sam]>0) {
	    TString sumName = "MCsum"; sumName += i+10*sam;
	    hMCsum[i] = (TH1F*)hStack[i][sam]->Clone(sumName);
	    isCreated++;
	  }
	} else {
	  if(SubisCreated>0) hMCsub[i]->Add(hStack[i][sam]);
	  if(SubisCreated==0 && nSample[sam]>0) {
	    TString subName = "MCsub"; subName += i+10*sam;
	    hMCsub[i] = (TH1F*)hStack[i][sam]->Clone(subName);
	    SubisCreated++;
	  }
	}
	if(i>1 && typeSample.Contains("Iso")) {
	  hStack[i][sam]->Add(hStack[i-2][sam]);
	}
      }
      if(doRoot) {fHistos->cd(); hStack[i][sam]->Write();}
    }
    if(i>1 && typeSample.Contains("Iso")) {
      hdata[i]->Add(hdata[i-2]);
      hMCsum[i]->Add(hMCsum[i-2]);
    }
    if(isCreated==0){
      cout<<"No events pass the cuts in channel "<<i<<endl;
      continue;
    }
    if(hMCsub[i]) {
      hdata[i]->Add(hMCsub[i],-1);
      if(hdata[i]->GetMinimum()>0) hdata[i]->SetMinimum(0);
    } else hdata[i]->SetMinimum(0);
    if(doRoot) {fHistos->cd(); hdata[i]->Write();}
    float maxi = hdata[i]->GetMaximum();
    if(hs[i]->GetMaximum()>maxi) maxi = hs[i]->GetMaximum();
    hdata[i]->SetMaximum(1.15*maxi);
    style4.setMarkers(hdata[i], 0.6, 20);
    hdata[i]->GetYaxis()->SetLabelFont(style4.nFont);    hdata[i]->GetXaxis()->SetLabelFont(style4.nFont);
    hs[i]->Draw("");
    hs[i]->GetYaxis()->SetLabelFont(style4.nFont);       hs[i]->GetXaxis()->SetLabelFont(style4.nFont);
    hs[i]->GetYaxis()->SetLabelSize(style4.LabelSize*1.15);   hs[i]->GetXaxis()->SetLabelSize(style4.LabelSize*1.15);
    hs[i]->GetYaxis()->SetNdivisions(style4.nDivisions); hs[i]->GetXaxis()->SetNdivisions(style4.nDivisions);

    hdata[i]->GetYaxis()->SetLabelFont(style4.nFont);       hdata[i]->GetXaxis()->SetLabelFont(style4.nFont);
    hdata[i]->GetYaxis()->SetLabelSize(style4.LabelSize*1.15);   hdata[i]->GetXaxis()->SetLabelSize(style4.LabelSize*1.15);
    hdata[i]->GetYaxis()->SetNdivisions(style4.nDivisions); hdata[i]->GetXaxis()->SetNdivisions(style4.nDivisions);

    hdata[i]->SetMaximum(hdata[i]->GetMaximum()*maxY);
    if(Variable.Contains("Q2") && hdata[i]->GetMinimum() < 0) hdata[i]->SetMinimum(-4);
    hdata[i]->Draw("a");
    hs[i]->Draw("same");
    hdata[i]->Draw("same");
    hdata[i]->Draw("axis same");
    if(Variable.Contains("candMva")){
      double BDTcut = BDTCuts[i];
      if(Variable == "candMvaDssDl") BDTcut = -0.45;
      if(Variable == "candMvaDssComb") BDTcut = -0.35;
      if(Variable == "candMvaDssComb" && i == 2) BDTcut = -0.3;
      line.DrawLine(BDTcut,0, BDTcut,maxi*1.15);
      arrow.DrawArrow(BDTcut,maxi*1.05,BDTcut+0.1,maxi*1.05,0.02,"|>");
    }
    TString padLeg = ""; padLeg += i;
    if(PosLeg!="no" && (style4.isThesis!="yes" || style4.isThesis=="yes"&&PosLeg.Contains(padLeg))){
      TLegend *leg;
      double legW = 0.27, legH = 0.57;
      double legX = 1-rightMargin, legY = 0.9;
      if(style4.isThesis=="yes") {
	legX-=0.012; legW=0.17; legH=0.59;
	legY = 0.97;
      }
      if(isDss) legH -= 0.15;
      if(isDpipi) legH += 0.07;
      if(style4.nFont==62) legW += 0.02;
      if(dEntries<10000) legW -= 0.02;
      if(PosLeg.Contains("l")) {
	legX = leftMargin+legW;
	if(style4.isThesis=="yes") legX+=0.04;
      }  
      leg = new TLegend(legX-legW, legY-legH, legX, legY);
      leg->SetTextFont(style4.nFont);
      leg->SetTextSize(style4.LabelSize*1.1);
      leg->SetFillColor(0); 
      if(style4.isThesis=="yes") {leg->SetBorderSize(0);leg->SetTextSize(style4.LabelSize*1.2);}
      TString dleg = "data ("; if(typeSample.Contains("signal")) dleg = "MC (";
      dleg += RoundNumber(dEntries,0); dleg += ")";
      if(style4.isThesis=="yes"){dleg = "data";}
      leg->AddEntry(hdata[i],dleg);
      int nhFirst = 8; if(isDpipi) nhFirst = 9;
      for(int nh2=nhFirst; nh2>=0; nh2--){
	int nh = indSam[isDss][nh2];
	if(i%2==0 && isDss==0){
	  if(nh==7) nh = 8;
	  else if(nh==8) nh = 7;
	}
	  //cout<<i<<": nh "<<nh<<" has "<<RoundNumber(nSample[nh],0)<<" entries"<<endl;
	if(nh==0 || nh==1 || nh==3) nSample[2] += nSample[3]+nSample[1]+nSample[0];
	else if(isDss &&(nh==7 || nh==8)) nSample[nh-2] += nSample[nh];
	else{
	  dleg = legTag[nh]; 
	  TString Sam=""; Sam+=nh;
	  if(subtract.Contains(Sam))dleg += "-";
	  dleg+=RoundNumber(nSample[nh],0); dleg+=")";
	  if(style4.isThesis=="yes"){dleg = legTag[nh]; dleg.ReplaceAll(" (","");}
	  leg->AddEntry(hStack[i][nh], dleg);
	}
      }
      leg->Draw();
    }

    Float_t min = hdata[i]->GetXaxis()->GetBinLowEdge(hdata[i]->GetXaxis()->GetFirst());
    Float_t max = hdata[i]->GetXaxis()->GetBinLowEdge(hdata[i]->GetXaxis()->GetLast()+1);
    hdata2[i] = (TH1F *)hdata[i]->Clone();
    hDiff[i]  = (TH1F *)hdata[i]->Clone();
    //for(int bin = 1; bin<hMCsum[i]->GetNbinsX()+1; bin++) hMCsum[i]->SetBinError(bin,0);
    //hdata2[i]->Add(hMCsum[i],-1);
    hdata2[i]->Divide(hMCsum[i]);
    hDiff[i]->Add(hMCsum[i],-1);
    int ndof=-1; double chi2=0, sErr=0, nErr=0;
    for(int bin = 1; bin<hdata2[i]->GetNbinsX()+1; bin++){
      double eBin = hdata2[i]->GetBinError(bin);
      if(eBin>0 && eBin<0.4){sErr += eBin; nErr++;}
    }
    double avg3Error = (double)((int)((3.*sErr/nErr+0.05)*10.))/10.;
    for(int bin = 1; bin<hDiff[i]->GetNbinsX()+1; bin++){
      double vRat = hdata2[i]->GetBinContent(bin);
      if(vRat>1.+avg3Error || vRat<1.-avg3Error) hdata2[i]->SetBinContent(bin,-10.);
      double vBin = hDiff[i]->GetBinContent(bin), eBin = hDiff[i]->GetBinError(bin);
      //double nBin = hdata[i]->GetBinContent(bin)+hMCsum[i]->GetBinContent(bin);
      //cout<<ndof<<": vBin "<<vBin<<", eBin "<<eBin<<", nBin "<<nBin<<"  -  chi2 "<<chi2<<endl;
      if(vBin != 0 && eBin >= sqrt(8)){
	chi2 += pow((vBin-1.)/eBin,2);
	ndof++;
      }
    }
    c.cd(pad);
    TPad *p2 = new TPad("p2","",0,0,1,0.2);
    p2->SetLeftMargin(leftMargin); p2->SetRightMargin(rightMargin); 
    p2->Draw(); p2->cd();
    if(hdata2[i]->GetMaximum()>-1.6) hdata2[i]->SetMaximum(1.+avg3Error);
    if(hdata2[i]->GetMinimum()<90.5) hdata2[i]->SetMinimum(1.-avg3Error);
    hdata2[i]->GetYaxis()->SetNdivisions(2+100*2);
    hdata2[i]->SetTitle("");
    hdata2[i]->SetXTitle("");
    hdata2[i]->GetYaxis()->SetLabelSize(style4.LabelSize*4.1);
    hdata2[i]->GetXaxis()->SetLabelSize(0.01);
    hdata2[i]->SetMarkerSize(.5);
    hdata2[i]->Draw("E0 P");
    if(doScale=="yes") ndof--;
    double ChiProb = Prob(chi2,ndof), xLabel = 1-rightMargin+0.01;
    TLine l;
    l.SetLineColor(2);
    l.DrawLine(min, 1.0, max, 1.0);
    c.cd(pad);
    label->SetTextAlign(33);
    TString rightLabel = "";
    if(style4.isThesis != "yes" || typeSample.Contains("Chi")) {
      rightLabel = "#splitline{#chi^{2}: "; rightLabel += RoundNumber(chi2,1); rightLabel += "/";
      rightLabel += ndof; rightLabel += " = "; rightLabel += RoundNumber(chi2,2,ndof);
      rightLabel += "}{Prob. = "; rightLabel += RoundNumber(100*ChiProb,1); rightLabel += "%}";

      label->SetTextSize(0.065); 
      if(PosLeg.Contains("l")) {label->SetTextAlign(13);label->DrawLatex(0.18,0.94,rightLabel);}
      else label->DrawLatex(xLabel-0.02,0.94,rightLabel);
      label->SetTextAlign(33);
    }
    label->SetTextSize(style4.TitleSize);
    if(Variable=="candPstarD"){xtitle = "i"; xtitle += i;}
    label->DrawLatex(xLabel,0.27,xtitle);
    double yFactor = 1; int digits = 2;
    if(Units.Contains("MeV")){ yFactor = 1000; digits = 0;}
    if((maxX-minX)*yFactor/((double)nbins)<1 && digits<1) digits = 1;
    TString ytitle = "Events/("; ytitle += RoundNumber((maxX-minX)*yFactor,digits,(double)nbins);
    ytitle += Units; ytitle += ")";
    if(Variable.Contains("Mult")) ytitle = "Events";
    label->SetTextAngle(90);
    label->DrawLatex(0.01,0.92,ytitle);
    label->SetTextAngle(0);label->SetTextAlign(11);
    double errRatio = sqrt(gEntries/pow(dEntries,2)+pow(gEntries,2)/pow(dEntries,4)*(dEntries+2*gSub));
    TString RatioTitle = titles[i]; RatioTitle += RoundNumber(gEntries,2,dEntries);
    RatioTitle+=" #pm "; RatioTitle += RoundNumber(errRatio,2);
    if(emu=="e")  RatioTitle+=",  Electrons";
    if(emu=="mu") RatioTitle+=",  Muons";
    if(typeSample.Contains("signal")) RatioTitle.ReplaceAll("data","MC");
    label->SetTextSize(style4.TextSize*0.92);
    if(style4.isThesis != "yes") label->DrawLatex(leftMargin,0.925,RatioTitle);
    else {
      RatioTitle = titles[i]; RatioTitle.ReplaceAll(": MC/data = ","");
      label->SetTextAlign(13); 
      double xLabel = leftMargin + 0.035;
      if(PosLeg.Contains("l")) {label->SetTextAlign(33); xLabel = 1-rightMargin - 0.02;}
      label->DrawLatex(xLabel,0.93,RatioTitle);
    }
    totYields[0] += gEntries; totYields[1] += dEntries; totYields[2] += gSub;
     cout<<"MC: "<<RoundNumber(gEntries,0)<<", data: "<<RoundNumber(dEntries,0)<<"\t  -  Ratio: "<< 
       RoundNumber(gEntries,5,dEntries) << " +- "<< RoundNumber(errRatio,5)<<endl;
    TString labRun = "";
    if(i==2){
      labRun = runs; labRun += " runs, ";labRun += nbins;labRun += " bins";
    }
    if(i==3){
      labRun = "Signal: ";
      for(int tau=0; tau<4; tau++){labRun += RoundNumber(nTau[tau],0); if(tau<3)labRun += ", ";}
    }
    label->SetTextSize(0.055);
    //if(style4.isThesis != "yes") label->DrawLatex(0.02,0.22,labRun);
  }
  if(style4.isThesis != "yes") {
    c.cd(0);
    TString Scaled_s = "not scaled"; 
    if(doScale=="yes")  Scaled_s = "scaled";
    label->SetTextSize(0.031);label->SetTextAlign(21);
    //label->DrawLatex(0.5,0.46,Scaled_s);
    //label->DrawLatex(0.5,0.490,"Areas");
  }

  if(doRoot) {fHistos->cd(); hYields.Write(); fHistos->Close(); fHistos->Delete();}
  double errRatio = sqrt(totYields[0]/pow(totYields[1],2)+
			 pow(totYields[0],2)/pow(totYields[1],4)*(totYields[1]+2*totYields[2]));
  double Ratio = totYields[0]/totYields[1];
  cout<<"MC: "<<RoundNumber(totYields[0],0)<<", data: "<<RoundNumber(totYields[1],0)<<"\t  -  Ratio: "<< 
    RoundNumber(Ratio,5) << " +- "<< RoundNumber(errRatio,5)<<"\t That's "<<
    RoundNumber(Ratio-1,2,errRatio)<<" sigma away"<<endl;
//   cout<<"Dss:\t";
//   for(int chan=0;chan<4;chan++) cout<<RoundNumber(nSamInt[chan][4],1)<<", ";
//   cout<<endl<<"Ds:\t";
//   for(int chan=0;chan<4;chan++) cout<<RoundNumber(nSamInt[chan][6],1)<<", ";
//   cout<<endl<<"D:\t";
//   for(int chan=0;chan<4;chan++) cout<<RoundNumber(nSamInt[chan][5],1)<<", ";
//   cout<<endl;

  tagName.ReplaceAll(".root",""); Variable.ReplaceAll("-","_");
  TString fileName = "keys/eps/CSample/"; fileName += typeSample;  fileName += tagName; fileName += "_";
  fileName += Variable;  fileName += "_"; fileName += runs;  fileName += "_"; fileName += emu;  fileName += ".eps";
  c.SaveAs(fileName);
  return 1; 
}

