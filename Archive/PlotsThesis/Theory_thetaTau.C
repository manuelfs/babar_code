#include "babar_code/NewPhysics/RateCalc.cc"
#include "babar_code/FF/BToDtaunu.cc"
#include "babar_code/FF/BToDstaunu.cc"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include "babar_code/Styles/Styles.cc"
#include "TMath.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TBox.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include "TLegend.h"
#include "TSpline.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout; 
using std::endl;

void Theory_thetaTau(int isSM=1, double tBmH1=0.3, double tBmH2=0.5, double tBmH3=1, int doScale=1){
  int nHis = 3;
  Styles style; style.setPadsStyle(2); style.applyStyle();
  BToDtaunu Dtaunu; BToDstaunu Dstaunu;


  int nBins = 1000, color[] = {8,4,2,1};
  double limTL[] = {-1, 1}, ml[] = {Dtaunu.mE,Dtaunu.mMu,Dtaunu.mTau,Dtaunu.mTau}, intTL[2][4], meanTL[2][4];
  TString yTitle[] = {"D0", "D1"}, hName;
//   TString yTitle[] = {"d#Gamma(Dl#nu)/dcos(#theta_{#tau}) (10^{16}GeV^{-1})", 
// 		      "d#Gamma(Dl#nu)/(dcos(#theta_{#tau})) (10^{16}GeV^{-1})"}, hName;
  TString xTitle = "cos#theta_{#tau}", legName[] = {"De", "Dmu", "Dtau", "#tau"};
  TString decName[] = {"Dtaunu ","D*taunu"}, PanelTag[] = {"a)", "b)"};
  double tL, tBmH[] = {0, tBmH1, tBmH2, tBmH3}, maxH[] = {0,0};
  TCanvas can; can.Divide(2,1);
  TH1F *hTL[2][4];
  double legW = 0.27, legH = 0.31;
  double legX = style.PadLeftMargin+0.05, legY = style.PadBottomMargin+0.04;
  if(isSM==0) {
    ml[0] = ml[2]; ml[1] = ml[2]; 
    color[0] = 2; color[2] = 8;
    nHis = 4;
    legH = 0.34;
  }
  TLegend leg(legX, legY, legX+legW, legY+legH);
  leg.SetTextSize(style.LabelSize); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont); leg.SetBorderSize(0);
  for(int isDs=0; isDs<=1; isDs++){
    can.cd(isDs+1);
    for(int histo=0; histo<nHis; histo++){
      hName = "hTL"; hName += isDs; hName += histo;
      hTL[isDs][histo] = new TH1F(hName,"",nBins,limTL[0], limTL[1]);
      hTL[isDs][histo]->SetLineWidth(2);
      hTL[isDs][histo]->SetLineColor(color[histo]);
      if(isSM==0) {
	if(isDs==0) Dtaunu._gSR = -Dtaunu.mb_quark*pow(tBmH[histo],2);
	else        Dstaunu._gSR = -Dstaunu.mb_quark*pow(tBmH[histo],2);
	legName[histo] = "t#beta/m_{H} = "; legName[histo] += RoundNumber(tBmH[histo],1);
	legName[histo] += " GeV^{-1}";
	if(fabs(tBmH[histo])<1e-6) legName[histo] = "SM"; 
	if(ml[histo]==Dtaunu.mE) legName[histo] += ", e";
	if(ml[histo]==Dtaunu.mMu) legName[histo] += ", #mu"; 
      }
      for(int bin=1; bin<=nBins; bin++){
	tL = hTL[isDs][histo]->GetBinCenter(bin);
	if(isDs==0) hTL[isDs][histo]->SetBinContent(bin,1e16*Dtaunu.ComputetL(acos(tL),1,ml[histo]));
	else hTL[isDs][histo]->SetBinContent(bin,1e16*Dstaunu.ComputetL(tL,1,ml[histo]));
	//else hTL[isDs][histo]->SetBinContent(bin,1e16*Dstaunu.Compute(tBmH1,tL,1,ml[histo]));
      }
      meanTL[isDs][histo] = hTL[isDs][histo]->GetMean();
      intTL[isDs][histo] = hTL[isDs][histo]->Integral();
      if(isDs==0) leg.AddEntry(hTL[isDs][histo], legName[histo]);
      if(doScale) {
	hTL[isDs][histo]->Scale(nBins/hTL[isDs][histo]->Integral()/(limTL[1]-limTL[0]));
	yTitle[isDs].ReplaceAll(" (10^{16}GeV^{-1})","");
      }
      if(maxH[isDs]<hTL[isDs][histo]->GetMaximum()) maxH[isDs]=hTL[isDs][histo]->GetMaximum();
    }
  }
  cout<<endl<<"tanBeta/mH: \t";
  for(int histo=0; histo<nHis; histo++)cout<<RoundNumber(tBmH[histo],1)<<"\t";
  for(int isDs=0; isDs<=1; isDs++){
    can.cd(isDs+1);
    for(int histo=0; histo<nHis; histo++){
      if(histo==0){
	hTL[isDs][histo]->SetMinimum(0); hTL[isDs][histo]->SetMaximum(maxH[isDs]*1.17);
	hTL[isDs][histo]->Draw("c");
	style.setTitles(hTL[isDs][histo],xTitle,yTitle[isDs],PanelTag[isDs]);
      } else hTL[isDs][histo]->Draw("c same");
    }
    if(isDs==1) leg.Draw();

    cout<<endl<<decName[isDs]<<" ->   mean tL:  ";
    for(int histo=0; histo<nHis; histo++) cout<<RoundNumber(meanTL[isDs][histo],2)<<"\t";
    cout<<"  Rate/RateSM:  ";
    for(int histo=0; histo<nHis; histo++) cout<<RoundNumber(intTL[isDs][histo],2,intTL[isDs][0])<<"\t";

  }
  cout<<endl<<endl;
  TString epsName = "public_html/Theory_thetaTau_"; epsName += isSM; epsName += ".eps";
  can.SaveAs(epsName);
  for(int isDs=0; isDs<=1; isDs++)
    for(int histo=0; histo<nHis; histo++) 
      hTL[isDs][histo]->Delete();
}


