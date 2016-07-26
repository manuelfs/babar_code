//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KeysPlot.cc,v 1.7 2012/03/04 00:33:03 manuelf Exp $
//
// Description:
//      KeysPlot - Plots mmiss and pl projections of a KEYS fit given
//      a histogram
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      10/02/12 manuelf -- Code moved to KeysUtils.cc
//      09/12/15 manuelf -- Created
//------------------------------------------------------------------------

#include "TString.h"
#include "DonutUtils/KeysUtils.cc"
#include "DonutUtils/Styles.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void PlotSingle(TString hname, TString dataFolder, int nbins = -1, double minM2 = -99., double maxM2 = -99., 
		double maxFactor = 1.,  double intFactor = 1., TString typePlot = "curve");

int main(int argc, char *argv[]){
  if (argc < 2 || argc > 10 ) {
    cout << "USAGE: KeysPlot HistoName [dataFolder=RAll] [typePlot=curve] [nbins=-1] [minM2=-99] "<< 
      "[maxM2=-99] [maxFactor=1] [intFactor=1] [weightFile]" << endl;
    return 0;
  }

  TString hname = argv[1];
  TString dataFolder = "RAll";
  if(argc>2) dataFolder = argv[2]; 
  TString typePlot = "curve";
  if(argc>3) typePlot = argv[3]; 
  int nbins = -1;;
  if(argc>4) {TString temp_s = argv[4]; nbins = temp_s.Atoi();} 
  double minM2 = -99.;
  if(argc>5) {TString temp_s = argv[5]; minM2 = temp_s.Atof();} 
  double maxM2 = -99.;
  if(argc>6) {TString temp_s = argv[6]; maxM2 = temp_s.Atof();} 
  double maxFactor = 1.;
  if(argc>7) {TString temp_s = argv[7]; maxFactor = temp_s.Atof();} 
  double intFactor = 1.;
  if(argc>8) {TString temp_s = argv[8]; intFactor = temp_s.Atof();} 
  TString sam, smoo, scaleP, Base;
  Int_t Sam, nM2bin, nPlbin;
  ParseName(hname, sam, smoo, scaleP, Base, Sam, nM2bin, nPlbin);
  TString weightName = "babar_code/Reweight/wTotal.txt";
  if (argc>9) weightName = argv[9];

  if(typePlot.Contains("pl") || typePlot.Contains("m2")) 
    PlotSingle(hname,dataFolder,nbins, minM2, maxM2,maxFactor,intFactor,typePlot);
  else PlotEps(hname, dataFolder, "Plot", weightName, "EPScomb", maxFactor, nbins, minM2, maxM2, intFactor, typePlot);
  
  return 1; 
}
  
void PlotSingle(TString hname, TString dataFolder, int nbins, double minX, double maxX, 
		double maxFactor,  double intFactor, TString typePlot) {
  
  TString sam, smoo, scaleP, Base;
  Int_t Sam, nM2bin, nPlbin, isPl=0;
  ParseName(hname, sam, smoo, scaleP, Base, Sam, nM2bin, nPlbin);
  if(typePlot.Contains("pl")) isPl=1;
  TString Variable = "candM2", otherVariable = "candPstarLep"; 
  if(isPl) {
    Variable="candPstarLep";
    otherVariable = "candM2";
  }
  Styles style1; 
  style1.setPadsStyle(1); style1.applyStyle();
  TString inputfile = nameData(Sam, dataFolder);
  TChain *treeData = new TChain("ntp1");
  treeData->Add(inputfile);
  TFile hfile(hname); hfile.cd();
  TH2F *h2 = (TH2F *)gDirectory->Get("h2");
  nM2bin = h2->GetNbinsX(); nPlbin = h2->GetNbinsY(); 
  TCanvas mm("mm","KEYS fits to mmiss-pstarl");
  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4;
  double minK = minX, maxK = maxX, minCut = maxFactor, maxCut = intFactor;
  int nbinsK = nbins;
  if(typePlot.Contains("curve")){
    if(isPl){
      nbinsK = nPlbin;
      minK = plmin; maxK = plmax;
    } else {
      nbinsK = nM2bin;
      minK = m2min; maxK = m2max;
    }
  }

  TString hdname = "hData"; 
  TH1F* hData = new TH1F("hData","",nbins,minX,maxX);
  TString vari = Variable; vari+=">>"; vari+=hdname; 
  TString Cuts = "("; Cuts+=otherVariable; Cuts+=">="; Cuts+=minCut; Cuts+="&&";
  Cuts+=otherVariable; Cuts+="<"; Cuts+=maxCut; Cuts+=")*weight";
  if(treeData) treeData->Draw(vari,Cuts);
  style1.setMarkers(hData, 0.7, 20);
  hData->SetMinimum(0);

  double entries = getWentries(treeData);
  if(entries<0.00001) entries = 1;
  double scale = entries*nbinsK*(maxX-minX)/nbins/(maxK-minK)/h2->Integral();

  TH1F hKeys = binHisto(h2,nbinsK,minK,maxK,minCut,maxCut,"hKeys",Variable);
  hKeys.Scale(scale);
  hKeys.SetLineColor(4);
  hKeys.SetLineWidth(1);
  hKeys.SetAxisRange(minX,maxX);
  float maxHisto = hKeys.GetMaximum();
  if(treeData){
    if(typePlot == "diff") hData->Add(&hKeys,-1);
    else {
      if(hData->GetMaximum()>maxHisto) maxHisto = hData->GetMaximum();
      if(maxHisto<1) maxHisto = 1.;
      hData->SetMaximum(maxHisto*1.2);
    }
  }

  TString units="GeV^{2})", otherUnits="GeV", nameVari = "m^{2}_{miss}", otherVari = "p*_{l}";
  if(isPl){
    units="GeV)"; nameVari = "p*_{l}"; otherVari = "m^{2}_{miss}"; otherUnits="GeV^{2}";
  }
  TString ytitle = "Entries/("; ytitle += RoundNumber((maxX-minX),2,(double)nbins); 
  ytitle += " "; ytitle += units;
  TString xtitle = nameVari; xtitle+=" ("; xtitle+=units;
  TString sCuts = RoundNumber(minCut,1); sCuts+=" #leq "; 
  if(fabs(minCut)<0.001) sCuts=otherVari; else sCuts+=otherVari; 
  sCuts+=" < "; sCuts += RoundNumber(maxCut,1); sCuts += " "; sCuts += otherUnits;
  if(isPl){
    if(minCut<=m2min && maxCut>=m2max) sCuts = "";
  } else {
    if(minCut<=plmin && maxCut>=plmax) sCuts = "";
  }
  TPad *cPad = (TPad *)mm.cd(0);
  TString optionPlot = "c same";
  if(!typePlot.Contains("curve")) optionPlot = "same";
  TLine lin; lin.SetLineColor(4);
  if(treeData){
    hData->Draw("e0");  
    style1.fixYAxis(hData,cPad);
    if(typePlot.Contains("diff"))lin.DrawLine(minX,0,maxX,0);
    hData->Draw("e1 same"); 
    if(!typePlot.Contains("diff")) hKeys.Draw(optionPlot);
    style1.setTitles(hData,xtitle,ytitle,"",sCuts);
  }else hKeys.Draw("c");
  TString epsname = "keys/eps/Plot/Single_"; epsname += dataFolder; 
  epsname += Base; epsname += ".eps";
  mm.SaveAs(epsname);
 
  h2->Delete(); hData->Delete();
}
