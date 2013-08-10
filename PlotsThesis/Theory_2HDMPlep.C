#include "TPad.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "babar_code/Styles/Styles.cc"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace TMath;
using namespace std;
using std::cout;
using std::endl;

TH1F SplitBin(TH1F *histo, TString hName);

void Theory_2HDMPlep(int PDFindex=5){

  int nHis = 4, nPads = 3;
  Styles style; style.setPadsStyle(nPads); 
  style.PadLeftMargin = 0.151; style.yTitleOffset = 1.06; 
  style.applyStyle();
  TCanvas can("can","2HDM PDFs");
  can.Divide(nPads,1); 


  double legW = 0.48, legH = 0.28;
  double legX = 1-style.PadRightMargin-0.01-legW, legY = 1-style.PadTopMargin-0.01;
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(style.LabelSize); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont); leg.SetBorderSize(0);


  TString Name, Variable[] = {"truePssLep", "truePstarLep", "truePLep"}, legName[4];
  TString xTitle[] = {"c", "d", "a"};
  TString yTitle[] = {"t3", "t2", "t1"};

  TString HigTag[] = {"000", "030", "050", "100"};
  TString Cuts = "MCType=="; Cuts += PDFindex;
  int nBins = 60, color[] = {2,4,8,1};
  //double maxH[] = {0,0,0}, tBmH[4], limX[3][2] = {{3, 12}, {0, 3.2}, {0, 2.2}};
  double maxH[] = {0,0,0}, tBmH[4], limX[3][2] = {{0, 1}, {0, 2.2}, {0, 3}};
  
  TH1F Histo[3][4], *hProv[3][4];
  for(int his=0; his<nHis; his++){
    TString cName = "AWG82/generator/4M"; cName += HigTag[his]; cName += "Dxtaunu.root";
    //TString cName = "AWG82/generator/Archive/4MEvtGen_Hig"; cName += HigTag[his]; cName += "Dxtaunu.root";
    TChain chain("GqaMCAnalysis/ntp1");
    chain.Add(cName);
    tBmH[his] = HigTag[his].Atof()/100.;
    legName[his] = "t#beta/m_{H} = "; legName[his] += RoundNumber(tBmH[his],1);
    if(fabs(tBmH[his])<1e-6) legName[his] = "SM"; 
    for(int pad=0; pad<nPads; pad++){
      can.cd(pad+1);
      Name = "Histo"; Name += pad; Name += his;
      hProv[pad][his] = new TH1F(Name, "",nBins,limX[pad][0],limX[pad][1]);
      TString vari = Variable[pad];
      chain.Project(Name, vari, Cuts);
      Name += "Split";
      Histo[pad][his] = SplitBin(hProv[pad][his], Name);
      Name += "Split";
      Histo[pad][his] = SplitBin(&Histo[pad][his], Name);
      Name += "Split";
      Histo[pad][his] = SplitBin(&Histo[pad][his], Name);
      Histo[pad][his].SetLineWidth(2);
      Histo[pad][his].SetLineColor(color[his]);
      Histo[pad][his].Scale(nBins*8/Histo[pad][his].Integral()/(limX[pad][1]-limX[pad][0]));
      if(maxH[pad]<Histo[pad][his].GetMaximum()) maxH[pad]=Histo[pad][his].GetMaximum();
      if(pad==0) leg.AddEntry(&Histo[pad][his], legName[his]);
    }
  }
    
  for(int pad=0; pad<nPads; pad++){
    can.cd(pad+1);
    for(int his=0; his<nHis; his++){
      if(his==0) {
	Histo[pad][his].SetMinimum(0);
	Histo[pad][his].SetMaximum(maxH[pad]*1.05);
	Histo[pad][his].Draw("c");
	style.setTitles(&Histo[pad][his],xTitle[pad],yTitle[pad]);
      }else Histo[pad][his].Draw("c same");
    }
    if(pad==2) leg.Draw();
  }
  TString epsName = "public_html/Theory_2HDMPlep_"; epsName += PDFindex; epsName += ".eps";
  can.SaveAs(epsName);
  for(int pad=0; pad<nPads; pad++)
    for(int his=0; his<nHis; his++)
      hProv[pad][his]->Delete();
}


TH1F SplitBin(TH1F *histo, TString hName){
  int nBins = histo->GetNbinsX();
  double minX = histo->GetXaxis()->GetBinLowEdge(1), maxX = histo->GetXaxis()->GetBinLowEdge(nBins+1);
  double prevBin = 0, thisBin = 0;

  TH1F hResult(hName,"",nBins*2,minX,maxX);
  for(int bin=1; bin<=nBins; bin++){
    thisBin = histo->GetBinContent(bin);
    if(bin>1) hResult.SetBinContent(2*bin-2,prevBin*2/3. + thisBin/3.);
    hResult.SetBinContent(2*bin-1,prevBin/3. + thisBin*2/3.);
    prevBin = thisBin;
  }
  hResult.SetBinContent(2*nBins, thisBin*2/3.);

  return hResult;
}
