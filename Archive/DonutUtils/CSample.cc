//------------------------------------------------------------------------
// File and Version Information:
//      $Id: CSample.cc,v 1.4 2012/08/23 02:22:17 manuelf Exp $
//
// Description:
//      CSample - Plots control samples
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      10/10/24 manuelf -- Updated chi^2 and set standard plot format
//      10/07/09 manuelf -- Created from mycode/CSampleBF.C
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
  if (argc < 2 || argc > 14 ) {
    cout << "USAGE: CSample typeSample [Variable=candM2] [nbins=20] [minX=-0.5] [maxX=1.5] "<< 
      "[extraCuts=1] [PosLeg=r] [doScale=yes] [weightFile=wTotal] [subtract] "<<
      "[emu=both] [ runs=All] [tagName] " << endl;
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
  TString extraCuts = "";
  if (argc>6) extraCuts = argv[6];
  extraCuts.ReplaceAll("XX","&&");
  TString PosLeg = "r";
  if (argc>7) PosLeg = argv[7];
  TString doScale = "yes";
  if (argc>8) doScale = argv[8];
  TString weightName = "babar_code/Reweight/wTotal.txt";
  if (argc>9) weightName = argv[9];
  TString subtract = "";
  if (argc>10) subtract = argv[10];
  TString emu = "both";
  if (argc>11) emu = argv[11];
  TString runs = "All";
  if (argc>12) runs = argv[12];
  TString tagName = "";
  if (argc>13) tagName = argv[13];

  int nSamples = 7, isDss = 0; bool isOff = false; 
  if(typeSample.Contains("Off")) {nSamples = 2; isOff = true;}
  if(typeSample.Contains("dss")) isDss = 1;

  TCut cuts = CSq2; if(typeSample=="pl") cuts = CSpl; 
  if(typeSample.Contains("MVA")) cuts = MvaAll;
  if(typeSample.Contains("noCut")) cuts = PMiss+M2P+MEScut;//+!MvaAll; 
  if(typeSample.Contains("noMES")) cuts = PMiss+M2P;//+!MvaAll; 
  if(typeSample.Contains("Reg2")) cuts = PMiss+M2P+lowQ2+Mva+MEScut; 
  if(typeSample.Contains("noMva")) cuts = PMiss+M2P+!MvaAll; 
  if(typeSample.Contains("noMVA")) cuts = PMiss+M2P+!Mva+MEScut;//+!MvaAll; 
  if(typeSample.Contains("Reg251")) cuts = PMiss+M2P+lowQ2+Mva51+MEScut; 
  if(isDss) cuts = dssMvaAll; if(isDss && typeSample.Contains("noCut")) cuts = dss+!MvaAll+MEScut;
  if(isDss && typeSample.Contains("MVA")) cuts = dss+dssMva+!MvaAll+MEScut;
  if(isDss && typeSample.Contains("Dl")) cuts = dss+dssMvaDl+cosT+!MvaAll+MEScut;
  if(isDss && typeSample.Contains("Comb")) cuts = dss+dssMvaComb+cosT+!MvaAll+MEScut;
  if(isDss && typeSample.Contains("CosT")) cuts = dss+cosT+!MvaAll+MEScut;
  if(isDss && typeSample.Contains("Maz")) cuts = dsseeAll;
  if(isDss && typeSample.Contains("51")) cuts = dssMvaAll51;
  if(isOff) cuts = MvaAll;    if(isOff && typeSample.Contains("noCut")) cuts = PMiss+M2P;  
  if(typeSample.Contains("signal")) {
    cuts = MvaAll; 
    if(typeSample.Contains("Mazur")) cuts = PMiss+M2P+ee;
    if(typeSample.Contains("noCut")) cuts = PMiss+M2P+MEScut; 
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
  double gEntries=0,uEntries=0,cEntries=0,dEntries=0,nSample[7];
  TChain *gen=0;

  fstream ResWFile; ResWFile.open("../DonutUtils/ResolutionWidths.txt",fstream::in);
  double ResWidth[4];
  for(int chan=0; chan<4; chan++) ResWFile >> ResWidth[chan];
  TString tupleFolder = "AWG82/ntuples/small/";  
  if(tagName!="") tupleFolder = "AWG82/ntuples/Newsmall/";
  tagName += ".root";
  TString dataName = tupleFolder; dataName += "Data_RunAll"; dataName += tagName;
  TString genName = tupleFolder; genName += "RAll_RunAll"; genName += tagName;
  TString udsName = tupleFolder; udsName += "uds_RunAll"; udsName += tagName;
  TString ccbarName = tupleFolder; ccbarName += "ccbar_RunAll"; ccbarName += tagName;
  if(typeSample.Contains("R24")) genName.ReplaceAll("RAll_","R24_");
  if(typeSample.Contains("R26")) genName.ReplaceAll("RAll_","R26_");
  if(isOff)  dataName = "AWG82/ntuples/small/OffPeak_RunAll.root";
  if(typeSample.Contains("signal")) dataName = genName;
  TString rootFile = "AWG82/ntuples/temp/data"; rootFile += typeSample; rootFile += "temp.root";
  WeightedTree3(dataName, dEntries, "1",0,cuts,"All", rootFile);
  TChain *data = new TChain("ntp1"); data->Add(rootFile);
  //TTree *data = WeightedTree(dataName, dEntries, "1",0,cuts, runs);
  //TTree *uds = WeightedTree(udsName, uEntries, weightName,0,cuts);
  //TTree *ccbar = WeightedTree(ccbarName, cEntries, weightName,0,cuts);
  rootFile = "AWG82/ntuples/temp/ccbar"; rootFile += typeSample; rootFile += "temp.root";
  WeightedTree3(ccbarName, cEntries, weightName,0,cuts,"All", rootFile);
  TChain *ccbar = new TChain("ntp1"); ccbar->Add(rootFile);
  rootFile = "AWG82/ntuples/temp/uds"; rootFile += typeSample; rootFile += "temp.root";
  WeightedTree3(udsName, uEntries, weightName,0,cuts,"All", rootFile);
  TChain *uds = new TChain("ntp1"); uds->Add(rootFile);
  if(!isOff) {
    rootFile = "AWG82/ntuples/temp/BB"; rootFile += typeSample; rootFile += "temp.root";
    WeightedTree3(genName, uEntries, weightName,0,cuts,"All", rootFile);
    gen = new TChain("ntp1"); gen->Add(rootFile);
    //gen = WeightedTree(genName, gEntries, weightName,0,cuts,runs);
  }
  
  if(tagName.Contains("Maz")) runs = "1234";
  double totMCB = 0, totuds = 0, totccbar = 0, totdata = 0, totOffdata = 0;
  getNumberB(genName, runs, totMCB, totdata, totuds, totccbar, totOffdata);
  double wuds = totMCB/totuds*2.09/1.05;              // 3.83785
  double wccbar = totMCB/totccbar*1.3/1.05;           // 
  double wMC = totdata/totMCB;                        // 4.63785
  if(isOff){ 
    wuds   = totOffdata/totuds*2.09*1000;             // 0.042759
    wccbar = totOffdata/totccbar*1.3*1000;            // 0.051672
    wMC = 1;
  }

  wccbar = wuds;                                      // weightManager converts ccbar into uds
  if(tagName.Contains("Maz")){wuds = 0; wccbar = 0;}

  TString samCut[2][5] = {{"MCType==0","MCType>6&&MCType<13","MCType>=13&&MCType<=14","(MCType==1||MCType==3||MCType==5)",
			   "(MCType==2||MCType==4||MCType==6)"},
			  {"MCType==0","MCType>0&&MCType<7","MCType>=13&&MCType<=14","(MCType==7||MCType==9||MCType==11)",
			   "(MCType==8||MCType==10||MCType==12)"}};
  int indSam[2][7] = {{0,1,2,3,4,5,6},{0,1,2,3,5,6,4}};
  TString legTag[7] = {"Cont (", "Cont (", "BB (", "Cross (", "D**l#nu (", "Dl#nu (", "D*l#nu ("};
  TString titles[] = {"D^{0}: MC/data = ","D*^{0}: MC/data = ","D^{+}: MC/data = ","D*^{+}: MC/data = "};
  TString xtitle = Variable, Units = " MeV";
  if(Variable=="candM2"||Variable == "mm2pi0"||Variable=="candM2NF"||Variable=="candQ2") Units = " GeV^{2}";
  if(Variable=="candCosT"||Variable.Contains("candTagNeutMult")||Variable.Contains("candTagChargedMult")
     ||Variable.Contains("Mva")) Units = "";
  if(Variable=="candM2" || Variable == "mm2pi0") xtitle = "m^{2}_{miss} (GeV^{2})"; 
  if(Variable=="candPstarLep") xtitle = "p*_{l} (GeV)"; 
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
  if(Variable=="ppi0") xtitle = "p(#pi^{0}) (GeV)"; if(Variable=="e1pi0") xtitle = "max(E_{#gamma,#pi^{0}}) (GeV)";
  if(Variable=="pmisspi0" || Variable=="candPMiss") xtitle = "p_{miss} (GeV)";
  if(isDss){
    if(Variable=="candM2") Variable = "mm2pi0";
    if(Variable=="candMvaDl") Variable = "candMvaDssDl";
    if(Variable=="candMvaComb") Variable = "candMvaDssComb";
  }


  int nPads = 4, pad; if(typeSample.Contains("Sum")) nPads = 1; 
  Styles style4; style4.setPadsStyle(nPads); 
  if(typeSample.Contains("Sum")) {
    style4.CanvasW = 380;
    style4.TitleSize = 0.066;
    titles[0] = "All channels: MC/data = ";
    samCut[0][1] = "(candType<3&&MCType>6&&MCType<13||candType>2&&MCType>0&&MCType<7)";
    samCut[0][3] = "(candType<3&&(MCType==1||MCType==3||MCType==5)||candType>2&&(MCType==7||MCType==9||MCType==11))";
    samCut[0][4] = "(candType<3&&(MCType==2||MCType==4||MCType==6)||candType>2&&(MCType==8||MCType==10||MCType==12))";
  }
  
  TLatex *label = new TLatex(); label->SetNDC(kTRUE); label->SetTextFont(style4.nFont);
  TLine line; line.SetLineStyle(2); line.SetLineColor(28); line.SetLineWidth(2);
  TArrow arrow; arrow.SetLineColor(28); arrow.SetFillColor(28); arrow.SetLineWidth(2);
  int colors[8] = {5,5,3,3,28,2,9,1};
  float nTau[4], nSamInt[4][7], totYields[]={0,0,0}; 
  if(isOff){colors[1] = 2; colors[7] = 5;}
  TH1F *hdata[4], *hdata2[4], *hDiff[4], *hStack[4][7], *hMCsum[4], *hMCsub[4]={0,0,0,0}, *hTau[4];


  int doRoot = 0;
  if(typeSample.Contains("Root")) doRoot = 1;
  TH1F hYields("hYields","",40,0,40);
  tagName.ReplaceAll(".root",""); 
  TString histosName = "keys/eps/CSample/root/"; histosName += typeSample;  histosName += tagName; 
  histosName += "_"; histosName += Variable;  histosName += "_"; histosName += runs;  histosName += "_"; 
  histosName += emu;  histosName += ".root";
  TFile *fHistos=0; if(doRoot) fHistos = new TFile(histosName,"RECREATE");

  THStack *hs[4];
  for (int i=0; i<4; i++) {TString hsName = "hs"; hsName+=i; hs[i] = new THStack(hsName,"");}
  gROOT->SetStyle("Plain");gStyle->SetOptStat(0);
  TCanvas c("dataMC","data Vs MC",style4.CanvasW,style4.CanvasH);
  if(!typeSample.Contains("Sum")) c.Divide(2,2);
  double leftMargin = 0.15, rightMargin = 0.1;
  if(style4.isThesis == "yes") rightMargin = style4.PadRightMargin;
  // Temporary for MVA fit
  // ====================================================================================
  TH1F *hFit[4]; TF1 *fFit[4];
  TCanvas cFit("FitCanvas","Linear MVA fit");
  cFit.Divide(2,2);
  // ====================================================================================
  for(int i=0; i<nPads; i++){
    pad = i+1; if(typeSample.Contains("Sum")) pad = 0;
    c.cd(pad);
    TPad *p1 = new TPad("p1","",0,.26,1,0.98);
    p1->Draw(); p1->cd(); 
    TString totCut = "candType=="; totCut += i+1; 
    if(typeSample.Contains("Sum")) totCut = "1";
    if(emu=="e")  totCut += "&&candIsMu==0";
    if(emu=="mu") totCut += "&&candIsMu==1";
    TString hname = "Data"; hname += i;
    TString vari = Variable; vari += ">>"; vari += hname;
    hdata[i] = new TH1F(hname,"",nbins,minX,maxX);
    hdata[i]->Sumw2();
    //cout<<i<<" Before projecting data: "<<hname<<", "<<Variable<<", "<<totCut<<endl;
    //cout<<"Entries data "<<data->GetEntries()<<endl;
    data->Project(hname,Variable,totCut);
    //cout<<"Done projecting data"<<endl;
    if(typeSample.Contains("signal")) hdata[i]->Scale(wMC);
    dEntries = hdata[i]->Integral();
    hYields.SetBinContent(i*10+9, dEntries);

    if(i==0) {
      style4.fixYAxis(hdata[i], p1);
      if(style4.nFont==62) leftMargin *= 1.1;
    }
    leftMargin = style4.PadLeftMargin; 
    p1->SetLeftMargin(leftMargin); p1->SetRightMargin(rightMargin);

    TString ContiCut = "("; ContiCut += totCut; ContiCut += ")*weight";
    hname = "MC"; hname += (i+1)*10;
    vari = Variable; vari.ReplaceAll("candM2NF","candM2");
    vari += ">>"; vari += hname;
    hStack[i][0]  = new TH1F(hname,"",nbins,minX,maxX);
    hStack[i][0]->Sumw2();
    nSample[0] = ccbar->Project(hname,Variable,ContiCut);
    hStack[i][0]->Scale(wccbar);
    hname = "MC"; hname += (i+1)*10+1;
    vari = Variable;  vari.ReplaceAll("candM2NF","candM2");
    vari += ">>"; vari += hname;
    hStack[i][1]  = new TH1F(hname,"",nbins,minX,maxX);
    hStack[i][1]->Sumw2();
    nSample[1] = uds->Project(hname,Variable,ContiCut);
    hStack[i][1]->Scale(wuds);
    int a = 0; if(i>1) a=1;
    TString TauCut = "("; TauCut += totCut; TauCut += "&&MCType=="; 
    TauCut += i+5+a*4; TauCut += ")*weight";
    if(!isOff){
      TString hTauName = "hTau"; hTauName += i;
      hTau[i] = new TH1F(hTauName,"",nbins,minX,maxX);
      vari = Variable; vari += ">>"; vari += hTauName;
      gen->Project(hTauName,Variable,TauCut);
      nTau[i] = hTau[i]->Integral()*wMC;
    } else nTau[i] = 0.;
    gEntries = 0; double gSub=0;
    for(int sam2 = 0; sam2<nSamples; sam2++){
      int sam = indSam[isDss][sam2];
      if(sam>1){
	TString genCut = "("; genCut += totCut; genCut += "&&"; 
 	genCut += samCut[a][sam-2]; genCut += ")*weight";
	//if(sam-2==3) genCut += "*0.444/0.30";
	//if(sam-2==4) genCut += "*0.333/0.25";
	hname = "MC"; hname += (i+1)*10+sam;
	vari = Variable; vari += ">>"; vari += hname; 
	hStack[i][sam]  = new TH1F(hname,"",nbins,minX,maxX);
	hStack[i][sam]->Sumw2();
	nSample[sam] = gen->Project(hname,Variable,genCut);
      }
      TString Sam=""; Sam+=sam;
      //      if(typeSample.Contains("Conv")&&!subtract.Contains(Sam)){
      if(typeSample.Contains("Conv")&&(i%2==0&&sam==5||i%2==1&&sam==6 || sam==4)){
	hname = "Sample"; hname += 10*i+sam;
	hStack[i][sam] = GaussConv(hStack[i][sam],ResWidth[i],hname,8,nbins,minX,maxX);
      }
      hStack[i][sam]->SetFillColor(colors[sam]);
      hStack[i][sam]->SetLineColor(colors[sam]);
      if(sam==0 || sam==2) {
	hStack[i][sam]->SetLineColor(colors[7]);
	if(!isOff)hStack[i][sam]->SetLineStyle(2);
      }
      double samInt = hStack[i][sam]->Integral();
      nSamInt[i][sam] = samInt;
      nSample[sam] = samInt*wMC;
      hYields.SetBinContent(i*10+sam, nSample[sam]);
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
      if(doRoot) {fHistos->cd(); hStack[i][sam]->Write();}
      if(nSample[sam]) {
	TString Sam=""; Sam+=sam;
	if(!subtract.Contains(Sam)){
	  if(doScale=="yes") hStack[i][sam]->Scale(dEntries/gEntries);
	  else hStack[i][sam]->Scale(wMC);
	  //if(i>1) hStack[i][sam]->Add(hStack[i-2][sam]);  //Prov
	  if(typeSample.Contains("signal") && sam<2) hdata[i]->Add(hStack[i][sam]);
	  if(nSample[sam]) hs[i]->Add(hStack[i][sam],"hist");
	  if(isCreated>0) hMCsum[i]->Add(hStack[i][sam]);
	  if(isCreated==0 && nSample[sam]>0) {
	    TString sumName = "MCsum"; sumName += i+10*sam;
	    hMCsum[i] = (TH1F*)hStack[i][sam]->Clone(sumName);
	    isCreated++;
	  }
	} else {
	  hStack[i][sam]->Scale(wMC);
	  if(SubisCreated>0) hMCsub[i]->Add(hStack[i][sam]);
	  if(SubisCreated==0 && nSample[sam]>0) {
	    TString subName = "MCsub"; subName += i+10*sam;
	    hMCsub[i] = (TH1F*)hStack[i][sam]->Clone(subName);
	    SubisCreated++;
	  }
	}
      }
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
    //if(i>1) hdata[i]->Add(hdata[i-2]);  //Prov
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
    if(PosLeg!="no"){
      TLegend *leg;
      double legW = 0.28, legH = 0.52;
      if(isOff){
	legW = 0.3; legH = 0.25;
      }
      if(style4.nFont==62) legW += 0.02;
      if(dEntries<10000) legW -= 0.02;
      if(PosLeg=="l") leg = new TLegend(leftMargin, 0.9-legH, leftMargin+legW, 0.9);
      else            leg = new TLegend(1-rightMargin-legW, 0.9-legH, 1-rightMargin, 0.9);
      leg->SetTextFont(style4.nFont);
      leg->SetTextSize(style4.LabelSize*1.1);
      leg->SetFillColor(0);
      TString dleg = "data ("; if(typeSample.Contains("signal")) dleg = "MC (";
      dleg += RoundNumber(dEntries,0); dleg += ")";
      leg->AddEntry(hdata[i],dleg);
      if(!isOff){
	for(int nh2=6; nh2>=0; nh2--){
	  int nh = indSam[isDss][nh2];
	  if(nh==1 || nh==3) nSample[nh-1] += nSample[nh];
	  else{
	    dleg = legTag[nh]; 
	    TString Sam=""; Sam+=nh;
	    if(subtract.Contains(Sam))dleg += "-";
	    dleg+=RoundNumber(nSample[nh],0); dleg+=")";
	    leg->AddEntry(hStack[i][nh], dleg);
	  }
	}
      } else {
	dleg = "uds ("; dleg+=RoundNumber(nSample[1],0); dleg+=")";
	leg->AddEntry(hStack[i][1], dleg);
	dleg = "ccbar ("; dleg+=RoundNumber(nSample[0],0); dleg+=")";
	leg->AddEntry(hStack[i][0], dleg);
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
    // Temporary for MVA fit
    // ====================================================================================
    if(Variable=="candMvaDl"){
      //cFit.cd(i+1);	
      double finVal[] = {0.9, 0.75, 0.85, 0.75};
      int iniBin = (int)((BDTCuts[i]-0.2-minX)*(double)nbins/(maxX-minX))+1;
      int finBin = (int)((finVal[i]-minX)*(double)nbins/(maxX-minX))+1;
      TString FitName = "hFit"; FitName += i;
      double RealIni = hdata[i]->GetXaxis()->GetBinLowEdge(iniBin);
      double RealFin = hdata[i]->GetXaxis()->GetBinLowEdge(finBin+1);
      int Fitbins = finBin-iniBin+1;
      hFit[i] = new TH1F(FitName,"",Fitbins,RealIni,RealFin);
      for(int bin = 1; bin<=Fitbins; bin++){
	hFit[i]->SetBinContent(bin,hdata2[i]->GetBinContent(bin+iniBin-1));
	hFit[i]->SetBinError(bin,hdata2[i]->GetBinError(bin+iniBin-1));
      }
      TString fName = "function"; fName += i;
      fFit[i] = new TF1(fName,"[0]*x+[1]",RealIni,RealFin);
      hFit[i]->Fit(fFit[i],"N Q");
      fFit[i]->SetLineWidth(2); fFit[i]->SetLineColor(4);
      //hFit[i]->Draw();
      cout<<RoundNumber(fFit[i]->GetParameter(0),2)<<" +- "<<RoundNumber(fFit[i]->GetParError(0),2)<<endl;
    }
    // ====================================================================================
    int dof=0; double chi2=0, sErr=0, nErr=0, sRat=0;
    for(int bin = 1; bin<hdata2[i]->GetNbinsX()+1; bin++){
      double eBin = hdata2[i]->GetBinError(bin);
      double rBin = hdata2[i]->GetBinContent(bin);      
      if(eBin>0 && eBin<0.4){
	sErr += eBin; nErr++; sRat += rBin;
      }
    }
    //double maxRat = (sRat+4*sErr)/nErr; 
    //double minRat = (sRat-4*sErr)/nErr; 
    double maxRat = (nErr+4*sErr)/nErr; 
    double minRat = (nErr-4*sErr)/nErr; 
    maxRat = (double)((int)((maxRat+0.05)*10.))/10.;
    minRat = (double)((int)((minRat+0.05)*10.))/10.;
    if(maxRat<=1) maxRat=1.1;
    if(minRat>=1) minRat=0.9;
    if(maxRat>1.9) maxRat=1.9;
    if(minRat<0.1) minRat=0.1;
    //maxRat=1.2; minRat=0.8;   //Prov
    //if(i==3){maxRat=1.4; minRat=0.8;}   //Prov
    for(int bin = 1; bin<hDiff[i]->GetNbinsX()+1; bin++){
      double vRat = hdata2[i]->GetBinContent(bin);
      if(vRat>maxRat || vRat<minRat) hdata2[i]->SetBinContent(bin,-10.);
      double vBin = hDiff[i]->GetBinContent(bin), eBin = hDiff[i]->GetBinError(bin);
      double nBin = hdata[i]->GetBinContent(bin)+hMCsum[i]->GetBinContent(bin);
      //cout<<dof<<": vBin "<<vBin<<", eBin "<<eBin<<", nBin "<<nBin<<"  -  chi2 "<<chi2<<endl;
      if(vBin != 0 && eBin != 0 && nBin>10){
	chi2 += pow((vBin-1.)/eBin,2);
	dof++;
      }
    }
    c.cd(pad);
    TPad *p2 = new TPad("p2","",0,0,1,0.2);
    p2->SetLeftMargin(leftMargin);  p2->SetRightMargin(rightMargin);
    p2->Draw(); p2->cd();
    if(hdata2[i]->GetMaximum()>-1.6) hdata2[i]->SetMaximum(maxRat);
    if(hdata2[i]->GetMinimum()<90.5) hdata2[i]->SetMinimum(minRat);
    hdata2[i]->GetYaxis()->SetNdivisions(2+100*2);
    hdata2[i]->SetTitle("");
    hdata2[i]->SetXTitle("");
    hdata2[i]->GetYaxis()->SetLabelSize(style4.LabelSize*4.1);
    hdata2[i]->GetXaxis()->SetLabelSize(0.01);
    hdata2[i]->SetMarkerSize(.5);
    hdata2[i]->Draw("E0 P");
    if(Variable=="candMvaDl" && style4.isThesis != "yes") fFit[i]->Draw("same");
    if(doScale=="yes") dof--;
    if(typeSample.Contains("Conv")) dof--;
    double ChiProb = Prob(chi2,dof);
    TLine l;
    l.SetLineColor(2);
    l.DrawLine(min, 1.0, max, 1.0);
    c.cd(pad);
    if(style4.isThesis != "yes") {
      TString ChiProb_s = RoundNumber(100*ChiProb,1); ChiProb_s += "%";
      label->SetTextSize(0.05);
      label->DrawLatex(0.91,0.164,"#tilde#chi^{2} =");
      label->DrawLatex(0.908,0.115,RoundNumber(chi2,2,(double)dof));
      label->DrawLatex(0.91,0.069,"p =");
      label->DrawLatex(0.908,0.02,ChiProb_s);
    }

    label->SetTextSize(style4.TitleSize);label->SetTextAlign(33);
    double xLabel = 1-rightMargin+0.01; 
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
    digits = 2; if(typeSample.Contains("Sum")) digits = 3; 
    TString RatioTitle = titles[i]; RatioTitle += RoundNumber(gEntries*wMC,digits,dEntries);
    RatioTitle+=" #pm "; RatioTitle += RoundNumber(errRatio*wMC,digits);
    if(emu=="e")  RatioTitle+=",  Electrons";
    if(emu=="mu") RatioTitle+=",  Muons";
    if(typeSample.Contains("signal")) RatioTitle.ReplaceAll("data","MC");
    label->SetTextSize(style4.TextSize*0.92);
    label->DrawLatex(leftMargin,0.925,RatioTitle);
    totYields[0] += gEntries; totYields[1] += dEntries; totYields[2] += gSub;
     cout<<"MC: "<<RoundNumber(gEntries*wMC,0)<<", data: "<<RoundNumber(dEntries,0)<<"\t  -  Ratio: "<< 
       RoundNumber(gEntries*wMC,5,dEntries) << " +- "<< RoundNumber(errRatio*wMC,5)<<endl;
    TString labRun = "";
    if(i<2 && extraCuts!=""){
      if(extraCuts.Length()<40) {if(i==0) labRun = extraCuts;}
      else {
	labRun = extraCuts;
	if(i==0) {labRun.Remove(labRun.Last('&'),labRun.Sizeof());}
	else {labRun.Remove(0,labRun.Last('&')+1);}
      }
    }
    if(i==2){
      labRun = runs; labRun += " runs, ";labRun += nbins;labRun += " bins";
    }
    if(i==3){
      labRun = "Signal: ";
      for(int tau=0; tau<4; tau++){labRun += RoundNumber(nTau[tau],0); if(tau<3)labRun += ", ";}
    }
    if(Variable=="candMvaDl") {
      labRun = "Slope: "; labRun += RoundNumber(fFit[i]->GetParameter(0),2);
      labRun += " #pm "; labRun += RoundNumber(fFit[i]->GetParError(0),2);
    }
    label->SetTextSize(0.055);
    if(style4.isThesis != "yes") label->DrawLatex(0.02,0.22,labRun);
  }
  if(style4.isThesis != "yes") {
    c.cd(0);
    TString Scaled_s = "not scaled"; 
    if(doScale=="yes")  Scaled_s = "scaled";
    label->SetTextSize(0.031);label->SetTextAlign(21);
    label->DrawLatex(0.5,0.46,Scaled_s);
    label->DrawLatex(0.5,0.490,"Areas");
  }
  if(doRoot) {fHistos->cd(); hYields.Write(); fHistos->Close(); fHistos->Delete();}
  double errRatio = sqrt(totYields[0]/pow(totYields[1],2)+
			 pow(totYields[0],2)/pow(totYields[1],4)*(totYields[1]+2*totYields[2]));
  double Ratio = totYields[0]*wMC/totYields[1];
  cout<<"MC: "<<RoundNumber(totYields[0]*wMC,0)<<", data: "<<RoundNumber(totYields[1],0)<<"\t  -  Ratio: "<< 
    RoundNumber(Ratio,5) << " +- "<< RoundNumber(errRatio*wMC,5)<<"\t That's "<<
    RoundNumber(Ratio-1,2,errRatio*wMC)<<" sigma away"<<endl;
//   cout<<"Dss:\t";
//   for(int chan=0;chan<4;chan++) cout<<RoundNumber(nSamInt[chan][4],1)<<", ";
//   cout<<endl<<"Ds:\t";
//   for(int chan=0;chan<4;chan++) cout<<RoundNumber(nSamInt[chan][6],1)<<", ";
//   cout<<endl<<"D:\t";
//   for(int chan=0;chan<4;chan++) cout<<RoundNumber(nSamInt[chan][5],1)<<", ";
//   cout<<endl;

  Variable.ReplaceAll("-","_");
  TString fileName = "keys/eps/CSample/"; fileName += typeSample;  fileName += tagName; fileName += "_";
  fileName += Variable;  fileName += "_"; fileName += runs;  fileName += "_"; fileName += emu;  fileName += ".eps";
  c.SaveAs(fileName);
//   cFit.SaveAs("FitMVA.eps");
  return 1; 
}

