//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: ResolutionFit.cc,v 1.3 2012/03/04 00:33:03 manuelf Exp $
//
// Description:
//      ResolutionFit - Finds the resolution of the mmiss peaks in data and MC
//                      It can either fit a Gaussian+Bifurcated gaussian, or
//                      fit the width of the gaussian to be convolved with the
//                      MC to get data
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      10/10/12 manuelf -- Created from WeightPl
//------------------------------------------------------------------------

#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TChain.h"
#include "TTree.h"
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

int main(int argc, char *argv[]){
  if (argc < 2 || argc > 8 ) {
    cout << "USAGE: ResolutionFit typeFit=Convo [minQ2=-5] [maxQ2=4] [minEex=0] [maxEex=8] [extraCuts] "<<
      "[weightBB]" << endl;
    return 0;
  }
  // Setting the input parameters
  TString typeFit = "Convo";
  if (argc>1) typeFit = argv[1];
  TString minQ2 = "-5";
  if (argc>2) minQ2 = argv[2];
  TString maxQ2 = "4";
  if (argc>3) maxQ2 = argv[3];
  TString minEex = "0";
  if (argc>4) minEex = argv[4];
  TString maxEex = "8";
  if (argc>5) maxEex = argv[5];
  TString cut_s = "1";
  if (argc>6) cut_s = argv[6];
  cut_s.ReplaceAll("XX","&&");
  TString weightName = "babar_code/Reweight/wTotal.txt";
  if (argc>7) weightName = argv[7];

  Styles style; style.setPadsStyle(4); style.PadTopMargin = 0.11; 
  style.applyStyle();
  TCanvas can("Resolution_fit","D(*)lnu Resolution fit");
  can.Divide(2,2); 
  if(style.isThesis=="yes") {
    if(typeFit.Contains("Convo")) style.setPadsStyle(1); 
    else if(typeFit.Contains("Roo")) style.setPadsStyle(2); 
    else {style.setPadsStyle(4); style.PadTopMargin = 0.11; }
  } else {style.setPadsStyle(4); style.PadTopMargin = 0.11; }
  style.applyStyle();
  TCanvas canD("Resolution_fitD","D(*)lnu Resolution fit");
  if(style.isThesis=="yes") {
    if(typeFit.Contains("Roo")) canD.Divide(2); 
  } else canD.Divide(2,2);

  cut_s += "&&candQ2>"; cut_s += minQ2; cut_s += "&&candQ2<="; cut_s += maxQ2; 
  cut_s += "&&candEExtra>="; cut_s += minEex.Atof()/10.; cut_s += "&&candEExtra<="; cut_s += maxEex.Atof()/10.; 
  TCut cuts = PMiss+M2P; cuts += cut_s; 
  if(typeFit.Contains("BDT")) cuts += Mva;
  double gEntries=0,uEntries=0,cEntries=0,dEntries=0;
  TString tupleFolder = "AWG82/ntuples/small/";
  TString dataName = tupleFolder; dataName += "Data_RunAll.root"; 
  TString genName = tupleFolder; genName += "RAll_RunAll.root"; 
  TString udsName = tupleFolder; udsName += "uds_RunAll.root"; 
  TString ccbarName = tupleFolder; ccbarName += "ccbar_RunAll.root"; 
  if(typeFit.Contains("R24")) genName.ReplaceAll("RAll_","R24_");
  TTree *gen = WeightedTree(genName, gEntries, weightName,0,cuts);
  TTree *data = WeightedTree(dataName, dEntries, "1",0,cuts);
  TTree *uds = WeightedTree(udsName, uEntries, weightName,-1,cuts);
  TTree *ccbar = WeightedTree(ccbarName, cEntries, weightName,-1,cuts);
  double totMCB = 0, totuds = 0, totccbar = 0, totdata = 0, totOffdata = 0;
  getNumberB(genName, "All", totMCB, totdata, totuds, totccbar, totOffdata);
  double wuds = totMCB/totuds*2.09/1.05;              // 0.862068
  double wccbar = totMCB/totccbar*1.3/1.05;           // 1.041766
  double wMC = totdata/totMCB;                        // 0.496529
  wccbar = wuds;                                      // weightManager converts ccbar into uds

  int nbins = 35, nbinsConv = 200;  // Enough bins to plot it smoothly (option "c")
  float minX = -0.3, maxX = 0.4;
  float valResW[] = {0.025, 0.017, 0.019, 0.023,0.025, 0.017, 0.019, 0.023}, errResW[8];
  float iniPar[] = {0, 0.04, 0, 0.14, 0.2};
  float minPar[] = {-0.2, 0.01, -0.2, 0.01, 0.01};
  TF1 *hConvolution[4];
  for(int chan=0; chan<4; chan++) {
    TString fName = "HistoConvo"; fName += chan;
    hConvolution[chan] = new TF1(fName,HistoConv,minX,maxX,2);
    hConvolution[chan]->SetParameter(0,0.025); hConvolution[chan]->SetParameter(1,1); 
    hConvolution[chan]->SetParLimits(0,0,0.1);
  }
  TLatex *label = new TLatex(); label->SetNDC(kTRUE); label->SetTextAlign(13);
  TString titles[] = {"D^{0}","D*^{0}","D^{+}","D*^{+}"};
  TString ParNames[] = {"Gmean", "Gsig", "BGmean", "BGsigL", "BGsigR", "NG", "NBG"}; 
  RooRealVar *ParFit[4][7][2], *ResWidth[4], Mmiss("Mmiss","m_{miss}^{2}",minX,maxX,"GeV^{2}");
  RooRealVar ResMean("ResMean","ResMean",0); ResMean.setConstant(0);
  RooGaussian *ResGauss[4];
  RooGaussian *Gauss[4][2];
  RooBifurGauss *BGauss[4][2];
  for(int chan=0; chan<4; chan++){
    TString name;
    for(int sam=0; sam<2; sam++){
      for(int par=0; par<7; par++){
	name = ParNames[par]; name += 100*(sam+1)+10*(chan+1)+par;
	if(par<5) ParFit[chan][par][sam] = new RooRealVar(name, name, iniPar[par], minPar[par], 0.6);
	else ParFit[chan][par][sam] = new RooRealVar(name, name, dEntries/8, 0, dEntries);
      }
      name = "Gauss"; name += 10*(sam+1)+chan+1;
      Gauss[chan][sam] = new RooGaussian(name, name, Mmiss, *ParFit[chan][0][sam], *ParFit[chan][1][sam]);
      name = "BGauss"; name += 10*(sam+1)+chan+1;
      BGauss[chan][sam] = new RooBifurGauss(name, name, Mmiss, 
					    *ParFit[chan][2][sam], *ParFit[chan][3][sam], *ParFit[chan][4][sam]);
    }
    name = "ResWidth"; name += (chan+1);
    ResWidth[chan] = new RooRealVar(name, name, valResW[chan]);
    name = "ResGauss"; name += (chan+1);
    ResGauss[chan] = new RooGaussian(name, name, Mmiss, ResMean, *ResWidth[chan]);
  }
  TString Variable = "candM2"; if(typeFit.Contains("NF")) Variable = "candM2NF";
  TString samCut[4] = {"(MCType==1||MCType==3)", "(MCType==2||MCType==4)",
		       "(MCType==7||MCType==9)", "(MCType==8||MCType==10)"};
  TH1F *hdata[4], *hDlnu[4], *hBkg[4], *hConti[4], *hConvDlnu[4];
  for(int i=0; i<4; i++){
    TPad *cPad = (TPad *)can.cd(i+1);
    TString totCut = "candType=="; totCut += i+1; 
    TString hname = "data"; hname += i;
    TString vari = Variable; vari += ">>"; vari += hname;
    hdata[i] = new TH1F(hname,"",nbins,minX,maxX);
    data->Draw(vari,totCut);
    hdata[i]->Sumw2();

    TString DlnuCut = "("; DlnuCut += totCut; DlnuCut += "&&"; DlnuCut += samCut[i]; 
    DlnuCut += ")*weight";
    hname = "Dlnu"; hname += i;
    vari = Variable; vari += ">>"; vari += hname;
    hDlnu[i]  = new TH1F(hname,"",nbins,minX,maxX);
    hDlnu[i]->Sumw2();
    gen->Draw(vari,DlnuCut);
    hDlnu[i]->Scale(wMC);
    TString BkgCut = "("; BkgCut += totCut; BkgCut += "&&!"; 
    BkgCut += samCut[i]; BkgCut += ")*weight";
    hname = "Bkg"; hname += i;
    vari = Variable; vari += ">>"; vari += hname;
    hBkg[i]  = new TH1F(hname,"",nbins,minX,maxX);
    hBkg[i]->Sumw2();
    gen->Draw(vari,BkgCut);
    hBkg[i]->Scale(wMC);

    TString ContiCut = "("; ContiCut += totCut; ContiCut += ")*weight";
    hname = "Conti"; hname += i;
    vari = "candM2"; vari += ">>"; vari += hname;
    hConti[i]  = new TH1F(hname,"",nbins,minX,maxX);
    hConti[i]->Sumw2();
    ccbar->Draw(vari,ContiCut);
    hConti[i]->Scale(wccbar*wMC);
    hBkg[i]->Add(hConti[i]);
    vari += "uds"; hname += "uds";
    hConti[i] = new TH1F(hname,"",nbins,minX,maxX);
    uds->Draw(vari,ContiCut);
    hConti[i]->Scale(wuds*wMC);
    hBkg[i]->Add(hConti[i]);

    if(typeFit.Contains("Roo")) hdata[i]->Add(hBkg[i],-1);
    else hDlnu[i]->Add(hBkg[i],1);

    TString PlotTitle = "";
    if(typeFit.Contains("Roo")){
      RooArgList pdfList(*Gauss[i][1], *BGauss[i][1]), parList(*ParFit[i][5][1],*ParFit[i][6][1]);
      RooAddPdf pdf("pdfSum","pdfSum",pdfList,parList);
      RooDataHist rooData(typeFit,typeFit,RooArgList(Mmiss),hDlnu[i]);
      RooChi2Var Chi2("Chi2","Chi2", pdf, rooData, kTRUE);
      RooMinuit rooMin(Chi2);
      rooMin.setPrintLevel(-1); rooMin.setVerbose(kFALSE);
      rooMin.setStrategy(2);    rooMin.optimizeConst(kTRUE);
      rooMin. hesse();          rooMin.migrad();
      for(int par=2; par<5; par++){
	ParFit[i][par][0]->setVal(ParFit[i][par][1]->getVal());
	ParFit[i][par][0]->setConstant();
      }
      RooArgList pdfListD(*Gauss[i][0], *BGauss[i][0]), parListD(*ParFit[i][5][0],*ParFit[i][6][0]);
      RooAddPdf pdfD("pdfSumD","pdfSumD",pdfListD,parListD);
      RooDataHist rooDataD("rooDataD","rooDataD",RooArgList(Mmiss),hdata[i]);
      RooChi2Var Chi2D("Chi2D","Chi2D", pdfD, rooDataD, kTRUE);
      RooMinuit rooMinD(Chi2D);
      rooMinD.setPrintLevel(-1); rooMinD.setVerbose(kFALSE);
      rooMinD.setStrategy(2);    rooMinD.optimizeConst(kTRUE);
      rooMinD. hesse();          rooMinD.migrad();

      if(style.isThesis=="yes") {
	if(i==3){
	  canD.cd(1);
	  RooPlot *plotM2 = Mmiss.frame();
	  rooData.plotOn(plotM2,DataError(RooAbsData::SumW2));
	  plotM2->getAttMarker()->SetMarkerSize(0.6);
	  pdf.plotOn(plotM2,Components(*BGauss[i][1]),LineStyle(kDashed));
	  plotM2->getAttLine()->SetLineColor(2);
	  plotM2->getAttLine()->SetLineWidth(2);
	  pdf.plotOn(plotM2);
	  plotM2->getAttLine()->SetLineWidth(2);
	  plotM2->SetTitle(""); 
	  plotM2->Draw();
	  label->DrawLatex(style.PadLeftMargin+0.04,1-style.PadTopMargin-0.03,"a)");  

	  canD.cd(2);
	  RooPlot *plotM2D = Mmiss.frame();
	  rooDataD.plotOn(plotM2D,DataError(RooAbsData::SumW2));
	  plotM2D->getAttMarker()->SetMarkerSize(0.6);
	  pdfD.plotOn(plotM2D,Components(*BGauss[i][0]),LineStyle(kDashed));
	  plotM2D->getAttLine()->SetLineColor(2);
	  plotM2D->getAttLine()->SetLineWidth(2);
	  pdfD.plotOn(plotM2D);
	  plotM2D->getAttLine()->SetLineWidth(2);
	  plotM2D->SetTitle(""); 
	  plotM2D->SetMinimum(0);
	  plotM2D->Draw();
	  label->DrawLatex(style.PadLeftMargin+0.04,1-style.PadTopMargin-0.03,"b)");  
	}
      } else {
	RooPlot *plotM2 = Mmiss.frame();
	rooData.plotOn(plotM2,DataError(RooAbsData::SumW2));
	pdf.plotOn(plotM2);
	plotM2->SetTitle("");
	plotM2->Draw();

	canD.cd(i+1);
	RooPlot *plotM2D = Mmiss.frame();
	rooDataD.plotOn(plotM2D,DataError(RooAbsData::SumW2));
	pdfD.plotOn(plotM2D);
	plotM2D->SetTitle("");
	plotM2D->Draw();
	PlotTitle = titles[i];
	label->DrawLatex(0.12,0.928,PlotTitle);
      }
    } else {
      if(style.isThesis=="yes" && i==3) cPad = (TPad *)canD.cd();
      hDlnu[i]->SetLineColor(28); 
      hDlnu[i]->Scale(hdata[i]->Integral()/hDlnu[i]->Integral());

      style.setMarkers(hdata[i], 0.6, 20);
      float maxi = hdata[i]->GetMaximum();
      if(hDlnu[i]->GetMaximum()>maxi) maxi = hDlnu[i]->GetMaximum();
      hdata[i]->SetMaximum(1.12*maxi);
      hdata[i]->Draw();
      style.fixYAxis(hdata[i],cPad);
      hDlnu[i]->Draw("hist same");
      TString ytitle = "Events/("; ytitle += RoundNumber((maxX-minX),2,(double)nbins); ytitle += " GeV^{2})";
      style.setTitles(hdata[i],"m^{2}_{miss} (GeV^{2})",ytitle);
      if(typeFit.Contains("Example")){
	double legW = 0.35, legH = 0.3;
	double legX = 1-style.PadRightMargin-0.01, legY = 1-style.PadTopMargin-0.01;
	TLegend *leg = new TLegend(legX-legW, legY-legH, legX, legY);
	leg->SetTextSize(style.LabelSize); leg->SetFillColor(0); leg->SetTextFont(style.nFont);
	leg->AddEntry(hdata[i],"Data");
	leg->AddEntry(hDlnu[i],"MC");
	int nExamples = 3, colors[] = {2, 3, 4, 12, 10};
	TH1F *hExample[5];
	PlotTitle = titles[i];
	double Res = 0.01, dRes = 0.020;
	for(int ex=0; ex<nExamples; ex++){
	  TString ExName = "ExName"; ExName += 10*i+ex;
	  hExample[ex] = GaussConv(hDlnu[i], Res, ExName, 2, nbinsConv, minX, maxX);
	  hExample[ex]->SetLineColor(colors[ex]); hExample[ex]->SetLineWidth(2);
	  hExample[ex]->Scale(hdata[i]->Integral()/hExample[ex]->Integral()*(double)nbinsConv/(double)nbins);
	  hExample[ex]->Draw("c same");
	  TString legTag = "#sigma = "; legTag += RoundNumber(Res,3); legTag += " GeV^{2}"; 
	  leg->AddEntry(hExample[ex],legTag);
	  Res += dRes;
	}
	hdata[i]->Draw("same");
 	leg->Draw();
     } else {       
	_hGlobal = (TH1F*)hDlnu[i]->Clone("hGlobal");
	hdata[i]->Fit(hConvolution[i], "N E Q B"); //N: No plot, E: Minos errors, Q: Not verbose, B: Bounds (par limits)
	valResW[i] = hConvolution[i]->GetParameter(0); errResW[i] = hConvolution[i]->GetParError(0);
	TString ConvName = "ConvDlnu"; ConvName += i;
	if(typeFit.Contains("Bad")){
	  RooDataHist rooDlnu("RooDlnu","RooDlnu",RooArgList(Mmiss),hDlnu[i]);
	  RooHistPdf pdfDlnu("pdfDlnu","pdfDlnu", RooArgList(Mmiss),rooDlnu,2);
	  RooNumConvPdf convDlnu("convDlnu","convDlnu",Mmiss,pdfDlnu,*ResGauss[i]);
	  hConvDlnu[i]  = new TH1F(ConvName,"",nbinsConv,minX,maxX);
	  convDlnu.fillHistogram(hConvDlnu[i],Mmiss);
	}else hConvDlnu[i]  = GaussConv(hDlnu[i], valResW[i], ConvName, 2, nbinsConv, minX, maxX);

	hConvDlnu[i]->SetLineColor(4); hConvDlnu[i]->SetLineWidth(2);
	hConvDlnu[i]->Scale(hdata[i]->Integral()/hConvDlnu[i]->Integral()*(double)nbinsConv/(double)nbins);
	hConvDlnu[i]->Draw("c same");
	hdata[i]->Draw("same");
	//PlotTitle = titles[i]; PlotTitle += ": #sigma_{fit} = ("; PlotTitle += RoundNumber(1000*sqrt(valResW[i]),0);
	//PlotTitle += " #pm "; PlotTitle += RoundNumber(1000*errResW[i]/sqrt(valResW[i]),0); PlotTitle += " MeV)^{2} "; 
	PlotTitle = titles[i]; PlotTitle += ": #sigma_{fit} = ("; PlotTitle += RoundNumber(1000*valResW[i],0);
	PlotTitle += " #pm "; PlotTitle += RoundNumber(1000*errResW[i],0); PlotTitle += ") #upoint 10^{-3} GeV^{2} "; 
	double legW = 0.15, legH = 0.24;
	double legX = 1-style.PadRightMargin-0.01, legY = 1-style.PadTopMargin-0.04;
	TLegend *leg = new TLegend(legX-legW, legY-legH, legX, legY);
	leg->SetTextSize(style.TitleSize); leg->SetFillColor(0);leg->SetBorderSize(0); leg->SetTextFont(style.nFont);
	leg->AddEntry(hdata[i],"Data");
	leg->AddEntry(hDlnu[i],"MC");
	leg->AddEntry(hConvDlnu[i],"Fit");
	leg->Draw();
      }
    }
    can.cd(i+1);
    if(style.isThesis=="no") label->DrawLatex(style.PadLeftMargin+0.03,0.924,PlotTitle);
  }
  cout<<endl;
  if(typeFit.Contains("Roo")){
    for(int sam=0; sam<2; sam++){
      for(int par=0; par<7; par++) cout<<ParNames[par] << "\t\t";
      cout<<endl;
      for(int chan=0; chan<4; chan++){
	for(int par=0; par<7; par++){
	  double factor = 1000.; int digits = 1;
	  if(par>4) {factor = 1; digits = 0;}
	  cout<<RoundNumber(ParFit[chan][par][sam]->getVal()*factor, digits) << " +- ";
	  cout<<RoundNumber(ParFit[chan][par][sam]->getError()*factor, digits) << "\t";
	  if(par==1) {
	    valResW[chan+sam*4] = ParFit[chan][par][sam]->getVal();
	    errResW[chan+sam*4] = ParFit[chan][par][sam]->getError();
	  }
	}
	cout<<endl;
      }
      cout<<endl;
    }
  } else if(typeFit.Contains("Convo")){
    for(int chan=0; chan<4; chan++){
      cout<<titles[chan]<<":\t"<<RoundNumber(1000*valResW[chan],1)<<" +- "<<RoundNumber(1000*errResW[chan],1)<<endl;
    }
  }
  double ve=0, e=0, ee=0, veR=0, eR=0, eeR=0;
  for(int chan=0; chan<4; chan++){
    ve += valResW[chan]/errResW[chan]/errResW[chan];
    e  += 1/errResW[chan];
    ee += 1/errResW[chan]/errResW[chan];
    if(typeFit.Contains("Roo")){
      veR += valResW[chan+4]/errResW[chan+4]/errResW[chan+4];
      eR  += 1/errResW[chan+4];
      eeR += 1/errResW[chan+4]/errResW[chan+4];
    }
  }
  if(!typeFit.Contains("Roo")){
    cout<<endl<<"Aver.:\t"<<RoundNumber(1000*ve,1,ee)<<" +- "<<RoundNumber(1000*1/sqrt(ee),1)<<endl<<endl;
    can.cd(0);
    TString labAve = "#sigma_{fit,av} = ("; labAve += RoundNumber(1000*ve,1,ee); labAve += " #pm ";
    labAve += RoundNumber(1000*1/sqrt(ee),1); labAve += ") #upoint 10^{-3} GeV^{2} ";
    label->SetTextSize(0.036);
    if(style.isThesis=="no") label->DrawLatex(0.53,0.02,labAve);    
  } else {
    double Wdata = ve/ee, WMC = veR/eeR;
    double Unsim = sqrt(Wdata*Wdata-WMC*WMC);
    double eUnsim = sqrt(Wdata*Wdata/ee+WMC*WMC/eeR)/Unsim;
    cout<<endl<<"Width data:\t"<<RoundNumber(1000*Wdata,1)<<" +- "<<RoundNumber(1000*1/sqrt(ee),1)<<endl;
    cout<<"Width MC:\t"<<RoundNumber(1000*WMC,1)<<" +- "<<RoundNumber(1000*1/sqrt(eeR),1)<<endl;
    cout<<"Width Unsim:\t"<<RoundNumber(1000*Unsim,1)<<" +- "<<RoundNumber(1000*eUnsim,1)<<endl<<endl;
  }
  TString plotName = "keys/eps/Resolution/"; plotName += typeFit; plotName+="_";
  plotName += minQ2; plotName+="-"; plotName += maxQ2; plotName+="_";
  plotName += minEex; plotName+="-"; plotName += maxEex; plotName+=".eps";
  can.SaveAs(plotName);
  plotName.ReplaceAll("Resolution/","Resolution/Data");
  if(typeFit.Contains("Roo") || style.isThesis=="yes") canD.SaveAs(plotName);

//   if(_hGlobal) _hGlobal->Delete();
//   for(int chan=0; chan<4; chan++){
//     if(hConvolution[chan]) hConvolution[chan]->Delete();
//     if(hdata[chan]    ) hdata[chan]    ->Delete();
//     if(hDlnu[chan]    ) hDlnu[chan]    ->Delete();
//     if(hConvDlnu[chan]) hConvDlnu[chan]->Delete();
//     if(hBkg[chan]     ) hBkg[chan]     ->Delete();
//     if(hConti[chan]   ) hConti[chan]   ->Delete();
//     if(Gauss[chan]    ) Gauss[chan]    ->Delete();
//     if(BGauss[chan]   ) BGauss[chan]   ->Delete();
//     for(int par=0; par<7; par++){
//       if(ParFit[chan][par]) ParFit[chan][par]->Delete();
//     }
//   }
  return 1;
}

//par[0] is the width and par[1] the normalization
Double_t HistoConv(Double_t *x, Double_t *par){
  int nbins = _hGlobal->GetNbinsX();
  double minX = _hGlobal->GetXaxis()->GetBinLowEdge(_hGlobal->GetXaxis()->GetFirst());
  double maxX = _hGlobal->GetXaxis()->GetBinLowEdge(_hGlobal->GetXaxis()->GetLast()+1);
  double nSigma = 4., wBin = (maxX-minX)/(double)nbins;
  double width = par[0];

  double minG = x[0]-(double)nSigma*width, maxG = x[0]+(double)nSigma*width;
  if(minG<minX) {minG=minX;} if(maxG>maxX) {maxG=maxX;} 
  double AreaG = IntG(x[0],width,minG,maxG), valH = 0;
  int iniBin = (int)((minG-minX)/wBin)+1, finBin = (int)((maxG-minX)/wBin)+1;
  for(int binH=iniBin; binH<=finBin; binH++){
    double minH = (double)(binH-1)*wBin+minX, maxH = (double)(binH)*wBin+minX;
    if(minH<minG) minH=minG; if(maxH>maxG) maxH=maxG; 
    valH += IntG(x[0], width, minH, maxH)*_hGlobal->GetBinContent(binH);
  }

  return par[1]*valH/AreaG;
}


