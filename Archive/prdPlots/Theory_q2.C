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

void Theory_q2(int isSM = 0, int doScale=1, double tBmH1=0.3, double tBmH2=0.5, double tBmH3=1){
  int nHis = 3;
  Styles style; style.setPadsStyle(-8); 
  style.CanvasH = 350;
  style.applyStyle();
  BToDtaunu Dtaunu; BToDstaunu Dstaunu;
  RateCalc RDs;
  RDs.ReadFF("babar_code/NewPhysics/FFinputs_HFAG11");


  int nBins = 1000, color[] = {8,2,4,1}, lineStyle[] = {3,2,1,4};
  double limQ2[] = {0, 12.}, ml[] = {Dtaunu.mE,Dtaunu.mMu,Dtaunu.mTau,Dtaunu.mTau}, intTL[2][4], meanTL[2][4];
//   TString yTitle[] = {"d#Gamma(Dlnu)/dq^{2} (10^{16}GeV^{-1})", 
// 		      "d#Gamma(D*lnu)/dq^{2} (10^{16}GeV^{-1})"}, hName, PanelTag[] = {"a)", "b)"};
  TString yTitle[] = {"x", "x"}, hName, decName[] = {"Dtaunu ","D*taunu"}, PanelTag[] = {"(a)", "(b)"};
  //TString xTitle = "q^{2} (GeV^{2})", legName[] = {"De", "Dmu", "Dtau", "#tau"};
  TString xTitle = "q", legName[] = {"De", "Dmu", "Dtau", "#tau"};
  double q2, tBmH[] = {0, tBmH1, tBmH2, tBmH3}, maxH[] = {0,0}, SMRD[] = {0.297, 0.252};
  TCanvas can; can.Divide(2,1);
  TH1F *hQ2[2][4];
  double legW = 0.35, legH = 0.4;
  double legX = 1-style.PadRightMargin-0.02, legY = 1-style.PadTopMargin-0.01;
  if(isSM<1) {
    limQ2[0] = 3;
    ml[0] = ml[2]; ml[1] = ml[2]; 
    color[0] = 4; color[2] = 8; 
    nHis = 4+isSM;
    legX = 0.5; legY = 0.94; 
    lineStyle[0] = 1; lineStyle[1] = 3; lineStyle[2] = 2; lineStyle[3] = 4; 
  } 
  TLegend leg(legX-legW, legY-legH, legX, legY);
  leg.SetTextSize(style.LabelSize*0.88); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont); leg.SetBorderSize(0);

  ///////////////////////////
  //double mB  = 5.2795, mDs = 2.01025; 
  //double mB2 = mB*mB, mDs2 = mDs*mDs; 
  ////////////////////////////////
  for(int isDs=0; isDs<=1; isDs++){
    can.cd(isDs+1);
    for(int histo=0; histo<nHis; histo++){
      hName = "hQ2"; hName += isDs; hName += histo;
      hQ2[isDs][histo] = new TH1F(hName,"",nBins,limQ2[0], limQ2[1]);
      hQ2[isDs][histo]->SetLineStyle(lineStyle[histo]);
      hQ2[isDs][histo]->SetLineWidth(3);
      hQ2[isDs][histo]->SetLineColor(color[histo]);
      if(isSM<1) {
	if(isDs==0) Dtaunu._gSR = -Dtaunu.mb_quark*pow(tBmH[histo],2);
	else        Dstaunu._gSR = -Dstaunu.mb_quark*pow(tBmH[histo],2);
	legName[histo] = "t#beta/m_{H} = "; legName[histo] += RoundNumber(tBmH[histo],1);
	legName[histo] += " GeV^{-1}";
	if(fabs(tBmH[histo])<1e-6) legName[histo] = "SM"; 
      }
      for(int bin=1; bin<=nBins; bin++){
	q2 = hQ2[isDs][histo]->GetBinCenter(bin);
	if(isDs==0) hQ2[isDs][histo]->SetBinContent(bin,1e16*Dtaunu.Compute(q2,1,ml[histo]));
	else hQ2[isDs][histo]->SetBinContent(bin,1e16*Dstaunu.Compute(q2,1,ml[histo]));
	if(histo==3) {
	  RDs.gsR = -RDs.mb*pow(tBmH[histo],2);
	  hQ2[isDs][histo]->SetBinContent(bin,1e16*RDs.Gamma_q2(isDs,q2,ml[histo]));
	}
	//////////// pD plot /////////////
	//double squ = mB2*mB2+mDs2*mDs2+q2*q2-2*(mB2*mDs2+mDs2*q2+q2*mB2); 
	//if(squ<0) squ=0;
	//hQ2[isDs][histo]->SetBinContent(bin,sqrt(squ)/(2*mB));
      }
      meanTL[isDs][histo] = hQ2[isDs][histo]->GetMean();
      intTL[isDs][histo] = hQ2[isDs][histo]->Integral();
      if(isDs==0) leg.AddEntry(hQ2[isDs][histo], legName[histo]);
      if(doScale) {
	hQ2[isDs][histo]->Scale(1000/hQ2[isDs][histo]->Integral());
	yTitle[isDs].ReplaceAll(" (10^{16}GeV^{-1})","");
      }
      if(maxH[isDs]<hQ2[isDs][histo]->GetMaximum()) maxH[isDs]=hQ2[isDs][histo]->GetMaximum();
    }
  }
  //maxH[0] = 0.18; maxH[1] = 0.19;
  cout<<endl<<"tanBeta/mH: \t";
  for(int histo=0; histo<nHis; histo++)cout<<RoundNumber(tBmH[histo],1)<<"\t";
  for(int isDs=0; isDs<=1; isDs++){
    can.cd(isDs+1);
    for(int histo=0; histo<nHis; histo++){
      if(histo==0){
	hQ2[isDs][histo]->SetMaximum(maxH[isDs]*1.15);
	hQ2[isDs][histo]->SetLabelOffset(0.01,"Y");
	hQ2[isDs][histo]->GetXaxis()->CenterTitle(true);
	if(isSM<1) {
	  //PanelTag[isDs] = "";
	  //hQ2[isDs][histo]->SetMaximum(0.198);
	}
	hQ2[isDs][histo]->Draw("c");
	//if(isSM<0) PanelTag[isDs] = "";
	style.setTitles(hQ2[isDs][histo],xTitle,yTitle[isDs], PanelTag[isDs]);
      } else hQ2[isDs][histo]->Draw("c same");
    }
    if(isDs==0 && isSM==1) leg.Draw();
    cout<<endl<<decName[isDs]<<" ->   mean tL:  ";
    for(int histo=0; histo<nHis; histo++) cout<<RoundNumber(meanTL[isDs][histo],2)<<"\t";
    cout<<"  Rate/RateSM:  ";
    for(int histo=0; histo<nHis; histo++) cout<<RoundNumber(intTL[isDs][histo]*SMRD[isDs],3,intTL[isDs][0])<<"\t";
  }
  cout<<endl<<endl;
  TString epsName = "public_html/Theory_q2_"; epsName += isSM; epsName += ".eps";
  can.SaveAs(epsName);

  if(isSM==0){
    style.CanvasH = 90;
    style.applyStyle();
    TCanvas can2; can2.cd();
    legX = 0.078; legY = 0.89; legW = 0.2; legH = 0.86; double shift = 0.47;
    TLegend leg2(legX, legY-legH, legX+legW, legY);
    leg2.SetTextSize(0.5); leg2.SetFillColor(0); leg2.SetTextFont(style.nFont); leg2.SetBorderSize(0);
    TLegend leg3(legX+shift, legY-legH, legX+legW+shift, legY);
    leg3.SetTextSize(0.5); leg3.SetFillColor(0); leg3.SetTextFont(style.nFont); leg3.SetBorderSize(0);
    TString legName2[] = {"a", "b", "c", "d"};
    for(int histo=0; histo<2; histo++) leg2.AddEntry(hQ2[0][histo], legName2[histo]);
    for(int histo=2; histo<4; histo++) leg3.AddEntry(hQ2[0][histo], legName2[histo]);
    leg2.Draw(); leg3.Draw();

    epsName = "public_html/Legend_tBmH.eps";
    can2.SaveAs(epsName);
  }
  for(int isDs=0; isDs<=1; isDs++)
    for(int histo=0; histo<nHis; histo++) 
      hQ2[isDs][histo]->Delete();
}

