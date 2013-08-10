//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: ResolutionDss.cc,v 1.2 2012/03/04 00:33:03 manuelf Exp $
//
// Description:
//      ResolutionDss - Finds the Gaussian to be convolved with the
//                      cocktail to get Generics in D**
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      10/12/15 manuelf -- Created from ResolutionFit
//------------------------------------------------------------------------

#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TChain.h"
#include "TChain.h"
#include "TSystem.h"
#include "RooFitCore/RooMinuit.hh"
#include "RooFitCore/RooArgList.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooFormulaVar.hh"
#include "RooFitCore/RooAddPdf.hh"
#include "RooFitCore/RooAbsReal.hh"
#include "RooFitCore/RooDataHist.hh"
#include "RooFitCore/RooHistPdf.hh"
#include "RooFitCore/RooChi2Var.hh"
#include "RooFitCore/RooPlot.hh"
#include "RooFitCore/RooAbsData.hh"
#include "RooFitCore/RooGlobalFunc.hh"
#include "RooFitCore/RooNumConvPdf.hh"
#include "RooFitModels/RooGaussian.hh"
#include "RooFitModels/RooBifurGauss.hh"
#include "DonutUtils/KeysUtils.cc"
#include "DonutUtils/Styles.cc"
#include "DonutUtils/cuts.cc"
#include <fstream>
#include <iostream>

using namespace RooFit;
using namespace std;
using std::cout;
using std::endl;

TH1F *_hGlobal;
Double_t HistoConv(Double_t *x, Double_t *par);
TH1F *BGaussConv(TH1F *h, float Lwidth, float Rwidth, TString hname, int PointsBin, int Cnbins, double CminX, double CmaxX,
		 double Goffset = 0);
TH1F *XConv(TH1F *h, double offset, double expo, double start, TString hname, int PointsBin, 
	    int Cnbins, double CminX, double CmaxX);

int main(int argc, char *argv[]){
  if (argc < 2 || argc > 5 ) {
    cout << "USAGE: ResolutionDss typeFit=Convo [Offset] [Expo] [Start] "<< endl;
    return 0;
  }
  // Setting the input parameters
  TString typeFit = "Convo";
  double Offset = 1;
  if (argc>2) {TString temp_s = argv[2]; Offset = temp_s.Atof();}
  double Expo = 2;
  if (argc>3) {TString temp_s = argv[3]; Expo = temp_s.Atof();}
  double Start = 0;
  if (argc>4) {TString temp_s = argv[4]; Start = temp_s.Atof();}

  Styles style4; 
  style4.setPadsStyle(1); style4.PadTopMargin = 0.11; style4.applyStyle();
  TString tupleFolder = "AWG82/ntuples/small/Fit";
  TString DssName = tupleFolder; DssName += "cocktail_RunAll.root"; 
  //TString DssName = tupleFolder; DssName += "Dpipi_RunAll.root"; 
  TString genName = tupleFolder; genName += "RAll_RunAll.root"; 
  TChain gen("ntp1"), Dss("ntp1");
  gen.Add(genName);
  Dss.Add(DssName);

  int nbins = 65, nbinsConv = 400;  // Enough bins to plot it smoothly (option "c")
  float minX = -0.6, maxX = 5;
  float parFit[4][4][2];
  TF1 *hConvolution[4];
  for(int chan=0; chan<4; chan++) {
    TString fName = "HistoConvo"; fName += chan;
    hConvolution[chan] = new TF1(fName,HistoConv,minX,maxX,4);
    hConvolution[chan]->SetParameter(0,1); 
    hConvolution[chan]->SetParameter(1,0.03); hConvolution[chan]->FixParameter(2,1.04); 
    hConvolution[chan]->FixParameter(3,-0.028);
    hConvolution[chan]->SetParLimits(1,0.01,0.2); hConvolution[chan]->SetParLimits(2,1.01,1.2);hConvolution[chan]->FixParameter(2,1.04); 
//     hConvolution[chan] = new TF1(fName,HistoConv,minX,maxX,4);
//     hConvolution[chan]->SetParameter(0,1); 
//     hConvolution[chan]->SetParameter(1,0.2); hConvolution[chan]->SetParameter(2,0.1); 
//     hConvolution[chan]->SetParameter(3,-0.1);
//     hConvolution[chan]->SetParLimits(1,0.001,1); hConvolution[chan]->SetParLimits(2,0.005,1);
  }
  TLatex *label = new TLatex(); label->SetNDC(kTRUE);
  TString titles[] = {"D^{0}","D*^{0}","D^{+}","D*^{+}"};
  TString Variable = "candM2"; 
  TH1F *hGen[4], *hDss[4], *hConvDlnu[4],*hEx[4][4];
  TCanvas can("Resolution_fit","D(*)lnu Resolution fit");
  //can.Divide(2,2); 
  for(int i=0; i<1; i++){
    //can.cd(i+1);
    //TString totCut = "(MCType==14&&candType=="; totCut += i+5; totCut += ")*weight";
    TString totCut = "(MCType==14)*weight";
    TString hname = "Gen"; hname += i;
    TString vari = Variable; vari += ">>"; vari += hname;
    hGen[i] = new TH1F(hname,"",nbins,minX,maxX);
    hGen[i]->Sumw2();
    gen.Draw(vari,totCut);

    hname = "Dss"; hname += i;
    vari = Variable; vari += ">>"; vari += hname;
    hDss[i]  = new TH1F(hname,"",nbins,minX,maxX);
    hDss[i]->Sumw2();
    Dss.Draw(vari,totCut);
    hDss[i]->Scale(hGen[i]->Integral()/hDss[i]->Integral());
    hDss[i]->SetLineColor(28); 

    style4.setMarkers(hGen[i], 0.6, 20);
    float maxi = hGen[i]->GetMaximum();
    if(hDss[i]->GetMaximum()>maxi) maxi = hDss[i]->GetMaximum();
    hGen[i]->SetMaximum(1.12*maxi);
    hGen[i]->Draw();
    hDss[i]->Draw("hist same");
    TString ytitle = "Events/("; ytitle += RoundNumber((maxX-minX),2,(double)nbins); ytitle += " GeV^{2})";
    style4.setTitles(hGen[i],"m^{2}_{miss} (GeV^{2})",ytitle);
    _hGlobal = (TH1F*)hDss[i]->Clone("hGlobal");
    hGen[i]->Fit(hConvolution[i], "N E Q B"); //N: No plot, E: Minos errors, Q: Not verbose, B: Bounds (par limits)
    for(int par=0; par<4; par++){
      parFit[i][par][0] = hConvolution[i]->GetParameter(par);
      parFit[i][par][1] = hConvolution[i]->GetParError(par);
    }
    hConvolution[i]->SetLineColor(7);
    hConvolution[i]->Draw("c same");
    //double Widths[3][3] = {{0.16, 0.0001, 0.13}, {0.16, 0.01, 0.13},  {0.16, 0.1, 0.13}};
    double Widths[3][3] = {{Offset, Expo, Start}, {0.03, 1.04, -0.028},  {0.03, 1.04, -0.029}};
    int Colors[3] = {2,3,9};
    for(int ex=0; ex<1; ex++) {
      TString exName = "ex"; exName+=i; exName+=ex;
//       hEx[i][ex] = BGaussConv(hDss[i], Widths[ex][0], Widths[ex][1], exName, 2, nbinsConv, minX, maxX, Widths[ex][2]);
//       hEx[i][ex]->SetLineColor(Colors[ex]); 
//       hEx[i][ex]->Scale(hGen[i]->Integral()/hEx[i][ex]->Integral()*(double)nbinsConv/(double)nbins);
//       hEx[i][ex]->Draw("c same");
       hEx[i][ex] = XConv(hDss[i], Widths[ex][0], Widths[ex][1], Widths[ex][2], exName, 2, nbinsConv, minX, maxX);
       hEx[i][ex]->SetLineColor(Colors[ex]); 
       hEx[i][ex]->Scale(hGen[i]->Integral()/hEx[i][ex]->Integral()*(double)nbinsConv/(double)nbins);
       hEx[i][ex]->Draw("c same");
    }
    TString ConvName = "ConvDlnu"; ConvName += i;
    hConvDlnu[i]  = BGaussConv(hDss[i], parFit[i][1][0], parFit[i][2][0], ConvName, 2, nbinsConv, minX, maxX, parFit[i][3][0]);
    hConvDlnu[i]  = XConv(hDss[i], parFit[i][1][0], parFit[i][2][0], parFit[i][3][0], ConvName, 3, nbinsConv, minX, maxX);
    hConvDlnu[i]->SetLineColor(4); hConvDlnu[i]->SetLineWidth(2);
    hConvDlnu[i]->Scale(hGen[i]->Integral()/hConvDlnu[i]->Integral()*(double)nbinsConv/(double)nbins);
    //hConvDlnu[i]->Draw("c same");
    hGen[i]->Draw("same");
    TString PlotTitle = titles[i]; PlotTitle += ": #sigma_{fit} = ("; PlotTitle += RoundNumber(1000*parFit[i][1][0],0);
    PlotTitle += " #pm "; PlotTitle += RoundNumber(1000*parFit[i][1][1],0); PlotTitle += ") #upoint 10^{-3} GeV^{2} "; 
    double legW = 0.155, legH = 0.225;
    double legX = 1-style4.PadRightMargin-0.01, legY = 1-style4.PadTopMargin-0.01;
    TLegend *leg = new TLegend(legX-legW, legY-legH, legX, legY);
    leg->SetTextSize(style4.TitleSize); leg->SetFillColor(0);leg->SetBorderSize(0); leg->SetTextFont(style4.nFont);
    leg->AddEntry(hGen[i],"Gen");
    leg->AddEntry(hDss[i],"Coc");
    leg->AddEntry(hConvDlnu[i],"Fit");
    leg->Draw();
    can.cd(i+1);
    label->DrawLatex(style4.PadLeftMargin+0.03,0.924,PlotTitle);
    
    int nDig = 3;
//     cout<<endl<<titles[i]<<":\t"<<RoundNumber(1000*parFit[i][1][0],1)<<" +- "<<RoundNumber(1000*parFit[i][1][1],1);
//     cout<<" \t "<<RoundNumber(1000*parFit[i][2][0],1)<<" +- "<<RoundNumber(1000*parFit[i][2][1],1);
    cout<<endl<<titles[i]<<":\t"<<RoundNumber(parFit[i][1][0],nDig)<<" +- "<<RoundNumber(parFit[i][1][1],nDig);
    cout<<" \t "<<RoundNumber(parFit[i][2][0],nDig)<<" +- "<<RoundNumber(parFit[i][2][1],nDig);
    cout<<" \t "<<RoundNumber(parFit[i][3][0],nDig)<<" +- "<<RoundNumber(parFit[i][3][1],nDig)<<endl;
  }
  TString plotName = "keys/eps/Resolution/Dss"; plotName += typeFit; plotName+=".eps";
  can.SaveAs(plotName);

  return 1;
}

Double_t IntX(double offset, double expo, double start, double minX, double maxX) {
  if(minX<start) minX = start; if(maxX<start) maxX = start;
  if((start+offset) < 0 || (minX+offset) < 0 || (maxX+offset) < 0) return 0;  
  return pow(start+offset,expo-1) * (pow(minX+offset,1-expo) - pow(maxX+offset,1-expo));
}

//par[0] the normalization, par[1] is the offset, par[2] the exponent, par[3] the start of the function
Double_t XHistoConv(Double_t *x, Double_t *par){

  double offset = par[1], expo = par[2], start = par[3];

  int nbins = _hGlobal->GetNbinsX();
  double minX = _hGlobal->GetXaxis()->GetBinLowEdge(_hGlobal->GetXaxis()->GetFirst());
  double maxX = _hGlobal->GetXaxis()->GetBinLowEdge(_hGlobal->GetXaxis()->GetLast()+1);
  double wBin = (maxX-minX)/(double)nbins;
  double minG = x[0]-maxX, maxG = x[0]-minX;
  if(minG<start) minG=start;
  if(minG>=maxG) return 0;

  double AreaG = IntX(offset, expo, start, minG, maxG), valH = 0;
  int iniBin = (int)((x[0]-maxG-minX)/wBin)+1, finBin = (int)((x[0]-minG-minX)/wBin)+1;
  for(int binH=iniBin; binH<=finBin; binH++){
    double maxH = x[0]-(double)(binH-1)*wBin-minX, minH = x[0]-(double)(binH)*wBin-minX;
    if(minH<minG) minH=minG; if(maxH>maxG) maxH=maxG; 
    valH += IntX(offset, expo, start, minH, maxH)*_hGlobal->GetBinContent(binH);
  }

  return par[0]*valH/AreaG;
}

TH1F *XConv(TH1F *h, double offset, double expo, double start, TString hname, int PointsBin, 
	    int Cnbins, double CminX, double CmaxX) {

  //_hGlobal = h;
  int Hnbins = h->GetNbinsX();
  double HminX = h->GetXaxis()->GetBinLowEdge(h->GetXaxis()->GetFirst());
  double HmaxX = h->GetXaxis()->GetBinLowEdge(h->GetXaxis()->GetLast()+1);
  TH1F *hConv;
  if(Cnbins<0) {
    hConv = (TH1F *)h->Clone(hname);
    Cnbins = Hnbins; CminX = HminX; CmaxX = HmaxX;
  } else hConv = new TH1F(hname, hname, Cnbins, CminX, CmaxX);
  if(CminX<HminX || CmaxX>HmaxX){
    cout<<"Range of convolved histogram has to be a subset of the original"<<endl;
    return 0;
  }
  double CwBin = (CmaxX-CminX)/(double)Cnbins;
  TF1 g("g",XHistoConv,CminX,CmaxX,4);
  g.SetParameter(0,1);
  g.SetParameter(1,offset);
  g.SetParameter(2,expo);
  g.SetParameter(3,start);

  double val[2000], dx = CwBin/((double)(PointsBin+1));
  for(int bin=0; bin<Hnbins; bin++) val[bin] = h->GetBinContent(bin+1);
  for(int binC=1; binC<=Cnbins; binC++) {
    double valC = 0;
    for(int point=0; point<PointsBin; point++){
      double x = CminX+(double)(binC-1)*CwBin+(double)(point+1)*dx;
      valC += g.Eval(x)/((double)PointsBin);
    }
    hConv->SetBinContent(binC,valC);
  }
  return hConv;
}

double IntBG(double mean, double Lw, double Rw, double minX, double maxX){
  if(minX < mean){
    if(maxX < mean) return IntG(mean, Lw, minX, maxX);
    else return (IntG(mean, Lw, minX, mean) + IntG(mean, Rw, mean, maxX)*Rw/Lw);
  } else return IntG(mean, Rw, minX, maxX)*Rw/Lw;
}

//par[0] the normalization, par[1/2] is the left/right width, and par[3] the offset
Double_t HistoConv(Double_t *x, Double_t *par){
  int nbins = _hGlobal->GetNbinsX();

  double minX = _hGlobal->GetXaxis()->GetBinLowEdge(_hGlobal->GetXaxis()->GetFirst());
  double maxX = _hGlobal->GetXaxis()->GetBinLowEdge(_hGlobal->GetXaxis()->GetLast()+1);
  double nSigma = 6., wBin = (maxX-minX)/(double)nbins;
  double Lwidth = par[1], Rwidth = par[2], Gmean = par[3];
  x[0] += Gmean;

  double minG = x[0]-(double)nSigma*Lwidth, maxG = x[0]+(double)nSigma*Rwidth;
  if(minG<minX) {minG=minX;} if(maxG>maxX) {maxG=maxX;} 
  double AreaG = IntBG(x[0],Lwidth,Rwidth,minG,maxG), valH = 0;
  int iniBin = (int)((minG-minX)/wBin)+1, finBin = (int)((maxG-minX)/wBin)+1;
  for(int binH=iniBin; binH<=finBin; binH++){
    double minH = (double)(binH-1)*wBin+minX, maxH = (double)(binH)*wBin+minX;
    if(minH<minG) minH=minG; if(maxH>maxG) maxH=maxG; 
    valH += IntBG(x[0],Lwidth,Rwidth,minH,maxH)*_hGlobal->GetBinContent(binH);
  }

  return par[0]*valH/AreaG;
}

TH1F *BGaussConv(TH1F *h, float Lwidth, float Rwidth, TString hname, int PointsBin, int Cnbins, double CminX, double CmaxX,
		double Goffset) {
  int Hnbins = h->GetNbinsX();
  double HminX = h->GetXaxis()->GetBinLowEdge(h->GetXaxis()->GetFirst());
  double HmaxX = h->GetXaxis()->GetBinLowEdge(h->GetXaxis()->GetLast()+1);
  TH1F *hConv;
  if(Cnbins<0) {
    hConv = (TH1F *)h->Clone(hname);
    Cnbins = Hnbins; CminX = HminX; CmaxX = HmaxX;
  } else hConv = new TH1F(hname, hname, Cnbins, CminX, CmaxX);
  if(CminX<HminX || CmaxX>HmaxX){
    cout<<"Range of convolved histogram has to be a subset of the original"<<endl;
    return 0;
  }
  double nSigma = 6.;
  double CwBin = (CmaxX-CminX)/(double)Cnbins, HwBin = (HmaxX-HminX)/(double)Hnbins;

  double val[2000], dx = CwBin/((double)(PointsBin+1)), totAreaG = IntBG(0,Lwidth,Rwidth,-nSigma*Lwidth,nSigma*Rwidth);
  for(int bin=0; bin<Hnbins; bin++) val[bin] = h->GetBinContent(bin+1);
  for(int binC=1; binC<=Cnbins; binC++) {
    double valC = 0;
    for(int point=0; point<PointsBin; point++){
      double x = CminX+(double)(binC-1)*CwBin+(double)(point+1)*dx+Goffset;
      double minG = x-(double)nSigma*Lwidth, maxG = x+(double)nSigma*Rwidth;
      double AreaG = totAreaG, valH = 0;
      bool lessRange = false;
      if(minG<HminX) {minG=HminX; lessRange = true;} if(maxG>HmaxX) {maxG=HmaxX; lessRange = true;} 
      if(lessRange) AreaG = IntBG(x,Lwidth,Rwidth,minG,maxG);
      int iniBin = (int)((minG-HminX)/HwBin)+1, finBin = (int)((maxG-HminX)/HwBin)+1;
      for(int binH=iniBin; binH<=finBin; binH++){
	double minH = (double)(binH-1)*HwBin+HminX, maxH = (double)(binH)*HwBin+HminX;
	if(minH<minG) minH=minG; if(maxH>maxG) maxH=maxG; 
	valH += IntBG(x, Lwidth, Rwidth, minH, maxH)*val[binH-1];
      }
      valC += valH/(AreaG*(double)PointsBin);
    }
    hConv->SetBinContent(binC,valC);
  }
  return hConv;
}


