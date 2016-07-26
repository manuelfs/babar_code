#include "TString.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TSystem.h"
#include "DonutUtils/MazurPlot.hh"
#include "DonutUtils/weightManager.hh"
#include "DonutUtils/KeysUtils.cc"
#include <fstream>
#include <iostream>

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

  RooDonutSemilep SL(initpeak, initmuL, initmuPeak,initmuR, initnuL, initnuR,initbo1,initbs1,initbc1,initbo2, 
		     initbs2, initbc2,initsl1,initsh1, initsl2,initsh2, initfl, initfh);
  RooDonutCB CB(initpeak, initmuL, initmuPeak,initmuR, initnuL, initnuR,initbo1,initbs1,initbc1,initbo2, 
		initbs2, initbc2,initsl1,initsh1, initsl2,initsh2, initbfrac, initbmean,
		initbsigr, initbsigl, initbpow,initbpow,initfl, initfh);

  double entries = 0.;
  Int_t isCocktail = 0;
  if(Sam<32) isCocktail = 1;
  TString inputfile = nameData(sam);
  TTree *treeData = WeightedTree(inputfile, entries, weightName, isCocktail);
  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4;
  double xlow = m2min,xhigh = m2max, ylow = plmin,yhigh = plmax;
  Int_t nbinx = 80, nbiny = 80, nM2bin = 400, nPlbin = 240;
  setBins(sam, xlow, xhigh, nbinx, nbiny);
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
	if(isSL) val += SL.evaluate(pl,mm);
	else val += CB.evaluate(pl,mm);
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
    for(double mm=-3; mm<4; mm+=0.04){
      if(isSL) val += SL.evaluate(pl,mm);
      else val += CB.evaluate(pl,mm);
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

Double_t RooDonutSemilep::evaluate(double pstarl, double mm2) 
{
  Double_t norm = getNorm(pstarl);
  Double_t arg1 = mm2 - getMean1(pstarl);
  Double_t arg2 = mm2 - getMean2(pstarl);
  Double_t sigma1 = getSigma1(pstarl);
  Double_t sigma2 = getSigma2(pstarl);
  Double_t frac = getFrac(pstarl);
  if (frac < 0) frac = 0;
  if (frac > 1) frac = 1;
  Double_t eval = norm * (frac*exp(-0.5*arg1*arg1/(sigma1*sigma1)) +
			  (1.0-frac)*exp(-0.5*arg2*arg2/(sigma2*sigma2)));
  return eval;
}

Double_t RooDonutSemilep::getNorm(Double_t y) const
{
  Double_t arg,mu,nu;
  if (y > peak) {
    arg = y - peak;
    mu = muPeak + (muR-muPeak)*(y-peak)/(2.4-peak);
  } else {
    arg = peak - y;
    mu = muL + (muPeak-muL)*(y)/(peak);
  }
  nu = nuL + (nuR-nuL)*(y)/2.4;
  arg /= mu;

  return exp(-1.*pow(arg,nu));
}

Double_t RooDonutSemilep::getMean1(Double_t y) const
{
  return biasOffset1 + biasSlope1*y + biasCurv1*y*y;
}

Double_t RooDonutSemilep::getMean2(Double_t y) const
{
  return biasOffset2 + biasSlope2*y + biasCurv2*y*y;
}

Double_t RooDonutSemilep::getSigma1(Double_t y) const
{
  return sigLow1 + (sigHigh1-sigLow1)*(y/2.4);
}

Double_t RooDonutSemilep::getSigma2(Double_t y) const
{
  return sigLow2 + (sigHigh2-sigLow2)*(y/2.4);
}

Double_t RooDonutSemilep::getFrac(Double_t y) const
{
  return fracLow + (fracHigh-fracLow)*(y/2.4);
}

RooDonutSemilep::RooDonutSemilep(double _peak, double _muL, double _muPeak,
				 double _muR, double _nuL, double _nuR,
				 double _biasOffset1, double _biasSlope1, double _biasCurv1,
				 double _biasOffset2, double _biasSlope2, double _biasCurv2,
				 double _sigLow1, double _sigHigh1, double _sigLow2,
				 double _sigHigh2, double _fracLow, double _fracHigh){
  peak       =_peak        ;
  muL        =_muL         ;
  muPeak     =_muPeak      ;
  muR        =_muR         ;
  nuL        =_nuL         ;
  nuR        =_nuR         ;
  biasOffset1=_biasOffset1 ;
  biasSlope1 =_biasSlope1  ;
  biasCurv1  =_biasCurv1   ;
  biasOffset2=_biasOffset2 ;
  biasSlope2 =_biasSlope2  ;
  biasCurv2  =_biasCurv2   ;
  sigLow1    =_sigLow1     ;
  sigHigh1   =_sigHigh1    ;
  sigLow2    =_sigLow2     ;
  sigHigh2   =_sigHigh2    ;
  fracLow    =_fracLow     ;
  fracHigh   =_fracHigh    ;
}

Double_t RooDonutCB::evaluate(double pstarl, double mm2) 
{
  Double_t norm = getNorm(pstarl);
  Double_t arg1 = mm2 - getMean1(pstarl);
  Double_t arg2 = mm2 - getMean2(pstarl);
  Double_t sigma1 = getSigma1(pstarl);
  Double_t sigma2 = getSigma2(pstarl);
  Double_t frac = getFrac(pstarl);
  if (frac < 0) frac = 0;
  if (frac > 1) frac = 1;
  Double_t gaus  = exp(-0.5*arg1*arg1/(sigma1*sigma1));
  Double_t gaus2 = exp(-0.5*arg2*arg2/(sigma2*sigma2));
  Double_t eval = norm * (frac*gaus + (1.0-frac)*gaus2);

  Double_t bf2 = bfrac * pow((2.4-pstarl)/(2.2),bpow);
  if (bf2 > 1.0) bf2 = 1.0;
  eval *= (1.0-bf2);
  if (mm2>bmean)
    eval += bf2*norm*exp(-0.5*(bmean-mm2)*(bmean-mm2)/(getBsigr(pstarl)*getBsigr(pstarl)));
  else
    eval += bf2*norm*exp(-0.5*(bmean-mm2)*(bmean-mm2)/(bsig2*bsig2));

  return eval;
}

Double_t RooDonutCB::getNorm(Double_t pstarl) const
{
  Double_t arg,mu,nu;
  if (pstarl > peak) {
    arg = pstarl - peak;
    mu = muPeak + (muR-muPeak)*(pstarl-peak)/(2.4-peak);
  } else {
    arg = peak - pstarl;
    mu = muL + (muPeak-muL)*(pstarl)/(peak);
  }
  nu = nuL + (nuR-nuL)*(pstarl)/2.4;
  arg /= mu;

  return exp(-1.*pow(arg,nu));
}

Double_t RooDonutCB::getMean1(Double_t y) const
{
  return biasOffset1 + biasSlope1*y + biasCurv1*y*y;
}

Double_t RooDonutCB::getMean2(Double_t y) const
{
  return biasOffset2 + biasSlope2*y + biasCurv2*y*y;
}

Double_t RooDonutCB::getSigma1(Double_t y) const
{
  return sigLow1 + (sigHigh1-sigLow1)*(y/2.4);
}

Double_t RooDonutCB::getSigma2(Double_t y) const
{
  return sigLow2 + (sigHigh2-sigLow2)*(y/2.4);
}

Double_t RooDonutCB::getFrac(Double_t y) const
{
  return fracLow + (fracHigh-fracLow)*(y/2.4);
}

Double_t RooDonutCB::getBsigr(Double_t y) const
{
  return bsig1 * pow(1.0-y/2.4,bsigpow);
}

RooDonutCB::RooDonutCB(double &_peak, double &_muL, double &_muPeak,
		       double &_muR, double &_nuL, double &_nuR,
		       double &_biasOffset1, double &_biasSlope1, double &_biasCurv1,
		       double &_biasOffset2, double &_biasSlope2, double &_biasCurv2,
		       double &_sigLow1, double &_sigHigh1, double &_sigLow2,
		       double &_sigHigh2, double &_bfrac, double &_bmean,
		       double &_bsig1, double &_bsig2, double &_bpow,
		       double &_bsigpow, double &_fracLow, double &_fracHigh){
   peak        = _peak        ;
   muL         = _muL         ;
   muPeak      = _muPeak      ;
   muR         = _muR         ;
   nuL         = _nuL         ;
   nuR         = _nuR         ;
   biasOffset1 = _biasOffset1 ;
   biasSlope1  = _biasSlope1  ;
   biasCurv1   = _biasCurv1   ;
   biasOffset2 = _biasOffset2 ;
   biasSlope2  = _biasSlope2  ;
   biasCurv2   = _biasCurv2   ;
   sigLow1     = _sigLow1     ;
   sigHigh1    = _sigHigh1    ;
   sigLow2     = _sigLow2     ;
   sigHigh2    = _sigHigh2    ;
   bfrac       = _bfrac       ;
   bmean       = _bmean       ;
   bsig1       = _bsig1       ;
   bsig2       = _bsig2       ;
   bpow        = _bpow        ;
   bsigpow     = _bsigpow     ;
   fracLow     = _fracLow     ;
   fracHigh    = _fracHigh    ; 
}



