#include "TString.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TSystem.h"
#include "RooFitCore/RooGlobalFunc.hh"
#include "RooFitCore/RooAbsData.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooPlot.hh"
#include "DonutUtils/RooDonutSemilep.hh"
#include "DonutUtils/RooDonutCB.hh"
#include "DonutUtils/weightManager.hh"
#include "DonutUtils/KeysUtils.cc"
#include <fstream>
#include <iostream>

using namespace RooFit;
using namespace std;
using std::cout;
using std::endl;

int main(int argc, char *argv[]){
  if (argc < 2 || argc > 3 ) {
    cout << "USAGE: MazurPlot sample [weightFile]" << endl;
    return 0;
  }

  TString Sam = argv[1];
  Int_t sam = Sam.Atoi();
  TString weightName = "wFF";
  if(sam==6 || sam==7 || sam==16 || sam==17 || (sam>=20&&sam<32)) weightName = "wMix";
  if (argc>2) weightName = argv[2];

  bool isSL = true;
  if(sam==0||sam==1||sam==2||sam==10||sam==11||sam==12||sam==20||
     sam==23||sam==26||sam==29) isSL = false;
  TString fileName="babar_code/fit/ParFit";fileName+=sam;fileName+=".txt";
  fstream parFile;
  parFile.open(fileName,fstream::in);
  // Reading file with parameter range
  Double_t lowpeak, initpeak, highpeak, lowmuL, initmuL, highmuL, lowmuPeak, initmuPeak, highmuPeak;
  Double_t lowmuR, initmuR, highmuR, lownuL, initnuL, highnuL, lownuR, initnuR, highnuR;  
  Double_t lowbc1, initbc1, highbc1, lowbo1, initbo1, highbo1, lowbs1, initbs1, highbs1;
  Double_t lowbc2, initbc2, highbc2, lowbo2, initbo2, highbo2, lowbs2, initbs2, highbs2;
  Double_t lowfl, initfl, highfl, lowfh, initfh, highfh, lowsh1, initsh1, highsh1;
  Double_t lowsh2, initsh2, highsh2, lowsl1, initsl1, highsl1;
  Double_t lowsl2, initsl2, highsl2;
  Double_t lowbfrac, initbfrac, highbfrac, lowbmean, initbmean, highbmean;  
  Double_t lowbsigr, initbsigr, highbsigr, lowbsigl, initbsigl, highbsigl;  
  Double_t lowbpow, initbpow, highbpow, lowbsigpow,initbsigpow, highbsigpow;
  TString pName; int Nfit;
  parFile>>pName>>Nfit;
  if(Nfit!=sam){cout<<"The file "<<fileName<<" is not for sample "<<sam<<endl; return 0;}
  parFile>>pName>>lowpeak>>initpeak>>highpeak;
  parFile>>pName>>lowmuL>>initmuL>>highmuL;
  parFile>>pName>>lowmuPeak>>initmuPeak>>highmuPeak;
  parFile>>pName>>lowmuR>>initmuR>>highmuR;
  parFile>>pName>>lownuL>>initnuL>>highnuL;
  parFile>>pName>>lownuR>>initnuR>>highnuR;
  parFile>>pName>>lowbo1>>initbo1>>highbo1;
  parFile>>pName>>lowbs1>>initbs1>>highbs1;
  parFile>>pName>>lowbc1>>initbc1>>highbc1;
  parFile>>pName>>lowsl1>>initsl1>>highsl1;
  parFile>>pName>>lowsh1>>initsh1>>highsh1;
  parFile>>pName>>lowbo2>>initbo2>>highbo2;
  parFile>>pName>>lowbs2>>initbs2>>highbs2;
  parFile>>pName>>lowbc2>>initbc2>>highbc2;
  parFile>>pName>>lowsl2>>initsl2>>highsl2;
  parFile>>pName>>lowsh2>>initsh2>>highsh2;
  parFile>>pName>>lowfl>>initfl>>highfl;
  cout<<"The file "<<fileName<<" opened"<<endl;
  parFile>>pName>>lowfh>>initfh>>highfh;
  if(!isSL){
      parFile>>pName>>lowbfrac>>initbfrac>>highbfrac;
      parFile>>pName>>lowbmean>>initbmean>>highbmean;
      parFile>>pName>>lowbsigr>>initbsigr>>highbsigr;
      parFile>>pName>>lowbsigl>>initbsigl>>highbsigl;
      parFile>>pName>>lowbpow>>initbpow>>highbpow;
      parFile>>pName>>lowbsigpow>>initbsigpow>>highbsigpow;
  }

  RooRealVar mmiss2("candM2","candM2",-4,12);
  RooRealVar pstarl("candPstarLep","candPstarLep",0.,2.4);
  RooRealVar totWeight("totWeight","totWeight",0.,100.);

  RooRealVar peak("peak","peak"        ,initpeak);
  RooRealVar muL("muL","muL"           ,initmuL);
  RooRealVar muPeak("muPeak","muPeak"  ,initmuPeak);
  RooRealVar muR("muR","muR"           ,initmuR);
  RooRealVar nuL("nuL","nuL"           ,initnuL);
  RooRealVar nuR("nuR","nuR"           ,initnuR);
  RooRealVar bo1("bo1","bias1 offset"  ,initbo1);
  RooRealVar bs1("bs1","bias1 slope"   ,initbs1);
  RooRealVar bc1("bc1","bias1 curv"    ,initbc1);
  RooRealVar bo2("bo2","bias2 offset"  ,initbo2);
  RooRealVar bs2("bs2","bias2 slope"   ,initbs2);
  RooRealVar bc2("bc2","bias2 curv"    ,initbc2);
  RooRealVar sl1("sl1","sig low 1"     ,initsl1);
  RooRealVar sh1("sh1","sig hi 1"      ,initsh1);
  RooRealVar sl2("sl2","sig low 2"     ,initsl2);
  RooRealVar sh2("sh2","sig hi 2"      ,initsh2);
  RooRealVar fl("fl","frac low"        ,initfl);
  RooRealVar fh("fh","frac hi"         ,initfh);
  RooRealVar bfrac("bfrac","bfrac"      ,initbfrac);
  RooRealVar bmean("bmean","bmean"      ,initbmean);
  RooRealVar bsigl("bsigl","bsigl"      ,initbsigl);
  RooRealVar bsigr("bsigr","bsigr"      ,initbsigr);
  RooRealVar bpow("bpow","bpow"         ,initbpow);
  RooRealVar bsigpow("bsigpow","bsigpow",initbsigpow);

  RooDonutSemilep SL("rds","rds",mmiss2,pstarl,peak,muL,muPeak,muR,nuL,nuR,
		      bo1,bs1,bc1,bo2,bs2,bc2,sl1,sh1,sl2,sh2,fl,fh);
  RooDonutCB CB("rdc","rdc",mmiss2,pstarl,peak,muL,muPeak,muR,nuL,nuR,
		 bo1,bs1,bc1,bo2,bs2,bc2,
		 sl1,sh1,sl2,sh2,bfrac,bmean,bsigr,bsigl,bpow,bsigpow,fl,fh);

  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4;
  double xlow = m2min,xhigh = m2max, ylow = plmin,yhigh = plmax;
  Int_t nbinx = 80, nbiny = 80, nM2bin = 400, nPlbin = 240;
  setBins(sam, xlow, xhigh, nbinx, nbiny);
  double entries = 0.;
  Int_t isCocktail = 0;
  if(Sam<32) isCocktail = 1;
  TString inputfile = nameData(sam);
  TTree *treeData = WeightedTree(inputfile, entries, weightName, isCocktail);
  RooDataSet data("data","data",RooArgSet(mmiss2,pstarl,totWeight));
  Float_t weight=0, candM2, candPstarLep;
  treeData->SetBranchAddress("weight",&weight);
  treeData->SetBranchAddress("candM2",&candM2);
  treeData->SetBranchAddress("candPstarLep",&candPstarLep);
  for (int evt = 0 ; evt < treeData->GetEntries() ; evt ++) {
    treeData->GetEvent(evt);
    mmiss2.setVal(candM2);
    pstarl.setVal(candPstarLep);
    totWeight.setVal(weight);
    data.add(RooArgSet(mmiss2,pstarl,totWeight));
  }
  data.setWeightVar(totWeight);
  TCanvas RooCanvas("RooCanvas","Roofit plot of p*l");
  pstarl.setBins(nbiny);
  RooPlot *pf = pstarl.frame();
  pf->SetName("pstarhist");
  pf->SetTitle("Roofit plot of p*l");
  data.plotOn(pf,DataError(RooAbsData::SumW2));
  CB.plotOn(pf);
  pf->Draw();
  RooCanvas.SaveAs("RoofitPl.eps");
  pstarl.setVal(1.9); mmiss2.setVal(0);
  cout<<"CB at mm=0 and pl=1.9 is "<<CB.evaluate()<<" and bfrac "<<bfrac.getVal()<<
    ", muL "<<nuL.getVal()<<endl;
  
  TString M2titles[] = {"0 < p*_{l} < 1 GeV","1 < p*_{l} < 1.4 GeV","1.4 < p*_{l} < 1.8 GeV",
		      "1.8 < p*_{l} < 2.4 GeV","0 < p*_{l} < 2.4 GeV"};
  TCut M2cuts[] = {"candPstarLep<1","candPstarLep>1&&candPstarLep<1.4",
		      "candPstarLep>1.4&&candPstarLep<1.8","candPstarLep>1.8&&candPstarLep<2.4", ""};
  TCut Plcuts[] = {"candM2<1","candM2>=1",""};
  TH1F *hm2[5], *m2[5], *hpl[3], *pl[3];
  TString M2names[5], Plnames[3];
  for(int i=0;i<5;i++){
    M2names[i] = "hm2_"; M2names[i] += i;
    hm2[i] = new TH1F(M2names[i],"",nM2bin,xlow,xhigh); 
    hm2[i]->SetLineColor(4);
    hm2[i]->SetLineWidth(1);
    TString hdname = "dm2"; hdname += i;
    TString vari = "candM2>>"; vari+=hdname; vari+="("; vari+= nbinx; vari+=",";vari+= xlow; 
    vari+=",";vari+= xhigh; vari+=")";
    M2cuts[i] *= "weight";
    treeData->Draw(vari,M2cuts[i]);
    m2[i] = (TH1F*)gDirectory->Get(hdname);
    m2[i]->SetXTitle("m^{2}_{miss} [GeV^{2}]");
    formatHisto(m2[i]);
    if(i<3) {
      Plnames[i] = "hpl_"; Plnames[i] += i;
      hpl[i] = new TH1F(Plnames[i],"",nPlbin,ylow,yhigh); 
      hpl[i]->SetLineColor(4);
      hpl[i]->SetLineWidth(1);
      TString hdname = "pl"; hdname += i;
      TString vari = "candPstarLep>>"; vari+=hdname; vari+="("; vari+= nbiny; vari+=",";vari+= ylow; 
      vari+=",";vari+= yhigh; vari+=")";
      Plcuts[i] *= "weight";
      treeData->Draw(vari,Plcuts[i]);
      pl[i] = (TH1F*)gDirectory->Get(hdname);
      pl[i]->SetXTitle("p*_{l} [GeV]");
      formatHisto(pl[i]);
    }
  }
  TLatex *label = new TLatex(); label->SetNDC(kTRUE); label->SetTextSize(0.055);
  TCanvas all6("All5","Mazur ansatz for fit "+Sam,1700,1800);
  all6.Divide(2,3,0.001,0.001);
  all6.cd(1);gPad->SetTopMargin(0.004);gPad->SetRightMargin(0.003);gPad->SetBottomMargin(0.11);
  gPad->SetFillColor(10);gPad->SetFrameFillColor(10);
  gStyle->SetOptStat(0);
  double hIntegral = 0;
  for(double mm=xlow; mm<xhigh; mm+=(xhigh-xlow)/400.){
    double limits[] = {0.,1.,1.4,1.8,2.4};
    for(int i=0; i<4; i++){
      double val = 0;
      for(double pl=limits[i]; pl<limits[i+1]; pl+=0.03){
	mmiss2.setVal(mm); pstarl.setVal(pl);
	if(isSL) val += 0;//SL.evaluate();
	else val += CB.evaluate();
      }
      hm2[i]->Fill(mm,val);
      hm2[4]->Fill(mm,val);
      hIntegral += val/(double)nM2bin;
    }
  }
  for(int i=0;i<5;i++){
    hm2[i]->Scale(entries/(double)nbinx/hIntegral);
    float maxHisto = hm2[i]->GetMaximum();
    if(m2[i]->GetMaximum()<maxHisto) m2[i]->SetMaximum(maxHisto*1.05);
  }
  hIntegral = 0;
  for(double pl=0.005; pl<2.4; pl+=0.01){
    double val = 0;
    for(double mm=-3; mm<5; mm+=0.04){
      mmiss2.setVal(mm); pstarl.setVal(pl);
      if(isSL) val += 0;//SL.evaluate();
      else val += CB.evaluate();
    }
    hpl[2]->Fill(pl,val);
    hIntegral += val/(double)nPlbin;
  }
  hpl[2]->Scale(entries/(double)nbiny/hIntegral);
  float maxHisto = hpl[2]->GetMaximum();
  if(pl[2]->GetMaximum()<maxHisto) pl[2]->SetMaximum(maxHisto*1.05);
  //cout<<"entries "<<entries<<", hIntegral "<<hIntegral<<endl;
  all6.cd(1);gPad->SetTopMargin(0.004);gPad->SetRightMargin(0.003);gPad->SetBottomMargin(0.11);
  gPad->SetFillColor(10);gPad->SetFrameFillColor(10);
  m2[4]->Draw("e0");  m2[4]->Draw("e1 same"); hm2[4]->Draw("c same");label->DrawLatex(0.15,0.92,"Mazur");
  all6.cd(2);gPad->SetTopMargin(0.004);gPad->SetRightMargin(0.003);gPad->SetBottomMargin(0.11);
  gPad->SetFillColor(10);gPad->SetFrameFillColor(10);
  pl[2]->Draw("e0");  pl[2]->Draw("e1 same"); hpl[2]->Draw("c same");
  for(int i=0; i<4; i++){
    all6.cd(i+3);gPad->SetTopMargin(0.003);gPad->SetRightMargin(0.003);gPad->SetBottomMargin(0.11);
    gPad->SetFillColor(10);gPad->SetFrameFillColor(10);
    m2[i]->Draw("e0"); label->DrawLatex(0.68,0.92,M2titles[i]);
    m2[i]->Draw("e1 same"); hm2[i]->Draw("c same"); 
  }
  TString cName = "babar_code/fit/eps/mmiss";cName+=sam;cName+=".eps";
  all6.SaveAs(cName);
  for(int i=0; i<5; i++)hm2[i]->Delete();


}



