//------------------------------------------------------------------------
// File and Version Information:
//      $Id: RateCalc.hh,v 1.1 2012/03/04 00:33:02 manuelf Exp $
//
// Description:
//      PlotRate - Plots the B->D(*)TauNu rates based on
//                 hep-ph 1203.2654
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      12/03/17 manuelf -- Created 
//------------------------------------------------------------------------

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

void PlotAngles(double tBmH1, double tBmH2, double tBmH3, double q2, 
		double thetaL, double thetaV=0, double chi=0, int isDgamma =0, int doScale=1){
  int isSM = 0, nHis =4;
  Styles style; style.setPadsStyle(2); style.applyStyle();
  BToDtaunu Dtaunu; BToDstaunu Dstaunu;

  double legW = 0.2, legH = 0.27;
  double legX = style.PadLeftMargin+0.27, legY = 1-style.PadTopMargin-0.52;

  int nBins = 1000, color[] = {2,4,28,1};
  double pi = 3.1415927;
  double limTL[] = {0, pi}, ml[] = {Dtaunu.mE,Dtaunu.mMu,Dtaunu.mTau,Dtaunu.mTau}, intTL[2][4], meanTL[2][4];
  TString yTitle[] = {"d#Gamma(D*#tau#nu)/d#theta_{V} (10^{16}GeV^{-1})", 
		      "d#Gamma(D*#tau#nu)/d#chi (10^{16}GeV^{-1})"}, hName;
  TString xTitle[] = {"#theta_{V}","#chi"}, legName[] = {"e", "#mu", "#tau", "#tau"};
  TString decName[] = {"Dtaunu ","D*taunu"};
  double tL, tBmH[] = {0, tBmH1, tBmH2, tBmH3}, maxH[] = {0,0};
  TCanvas can; can.Divide(2,1);
  TH1F *hTL[2][4];
  if(isSM==0) {
    ml[0] = ml[2]; 
    ml[1] = ml[2]; 
  }
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(style.LabelSize); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont); leg.SetBorderSize(0);
  for(int isDs=0; isDs<=1; isDs++){
    can.cd(isDs+1);
    for(int histo=0; histo<nHis; histo++){
      hName = "hTL"; hName += isDs; hName += histo;
      hTL[isDs][histo] = new TH1F(hName,"",nBins,limTL[0], limTL[1]);
      hTL[isDs][histo]->SetLineWidth(1);
      hTL[isDs][histo]->SetLineColor(color[histo]);
      if(isSM==0) {
	Dstaunu._gSR = -Dstaunu.mb_quark*pow(tBmH[histo],2);
	legName[histo] = "t#beta/m_{H} = "; legName[histo] += RoundNumber(tBmH[histo],1);
	if(fabs(tBmH[histo])<1e-6) legName[histo] = "SM"; 
	if(ml[histo]==Dtaunu.mE) legName[histo] += ", e";
	if(ml[histo]==Dtaunu.mMu) legName[histo] += ", #mu";
      }
      for(int bin=1; bin<=nBins; bin++){
	tL = hTL[isDs][histo]->GetBinCenter(bin);
	if(q2>0){
// 	  if(bin==500 || bin ==600){
// 	    cout<<"q2 "<<q2<<", thetaL "<<thetaL<<", thetaV "<<thetaV<<", chi "<<chi<<", tL "<<tL
// 		<<". Compute "<<Dstaunu.Compute(q2,cos(thetaL),cos(tL),chi,isDgamma,true,1,ml[histo])<<endl;
// 	  }
	  if(isDs==0) hTL[isDs][histo]->SetBinContent(bin,1e16*Dstaunu.Compute(q2,cos(thetaL),cos(tL),chi,isDgamma,true,1,ml[histo]));
	  else hTL[isDs][histo]->SetBinContent(bin,1e16*Dstaunu.Compute(q2,cos(thetaL),cos(thetaV),tL,isDgamma,true,1,ml[histo]));
	}else{
	  double q2min = ml[histo]*ml[histo], q2max = pow(Dstaunu._mB-Dstaunu._mDs,2);
	  int nPoints = 1000;
	  double dq2 = (q2max-q2min)/(double)nPoints, IntegratedTL=0, valq2=q2min, valRate;
	  //if(bin==1)cout<<"q2min "<<q2min<<", q2max "<<q2max<<", dq2 "<<dq2<<endl;
	  int isCLN=1; //if(histo==1) isCLN=0;
	  for(int iq2=0; iq2<=nPoints; iq2++){
	    if(isDs==0) valRate = Dstaunu.Compute(valq2,cos(thetaL),cos(tL),chi,isDgamma,true,isCLN,ml[histo]);
	    else        valRate = Dstaunu.Compute(valq2,cos(thetaL),cos(thetaV),tL,isDgamma,true,isCLN,ml[histo]);
	    if(iq2==0 || iq2==nPoints) valRate /= 2.;
	    IntegratedTL += valRate;
	    valq2 += dq2; 
	  }
	  hTL[isDs][histo]->SetBinContent(bin,IntegratedTL);
	}
      }
      meanTL[isDs][histo] = hTL[isDs][histo]->GetMean();
      intTL[isDs][histo] = hTL[isDs][histo]->Integral();
      if(isDs==0) leg.AddEntry(hTL[isDs][histo], legName[histo]);
      if(doScale) {
	hTL[isDs][histo]->Scale(100/hTL[isDs][histo]->Integral());
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
	TString sTheta = "q^{2} = "; 
	if(q2>0)sTheta += RoundNumber(q2,1); else sTheta += "All";
	sTheta += ", #theta_{L} = "; sTheta += RoundNumber(thetaL,1);
	if(isDs==0) {
	  sTheta += ", #chi = "; sTheta += RoundNumber(chi,1);
	} else {
	  sTheta += ", #theta_{V} = "; sTheta += RoundNumber(thetaV,1);
	}
	hTL[isDs][histo]->SetMinimum(0); hTL[isDs][histo]->SetMaximum(maxH[isDs]*1.17);
	hTL[isDs][histo]->Draw("c");
	style.setTitles(hTL[isDs][histo],xTitle[isDs],yTitle[isDs],sTheta);
      } else hTL[isDs][histo]->Draw("c same");
    }
    if(isDs==1) leg.Draw();

    cout<<endl<<decName[isDs]<<" ->   mean tL:  ";
    for(int histo=0; histo<nHis; histo++) cout<<RoundNumber(meanTL[isDs][histo],2)<<"\t";
    cout<<"  Rate/RateSM:  ";
    for(int histo=0; histo<nHis; histo++) cout<<RoundNumber(intTL[isDs][histo],2,intTL[isDs][0])<<"\t";

  }
  cout<<endl<<endl;
  TString epsName = "public_html/RD_Angles_"; epsName += doScale; epsName += ".eps";
  can.SaveAs(epsName);
  for(int isDs=0; isDs<=1; isDs++)
    for(int histo=0; histo<nHis; histo++) 
      hTL[isDs][histo]->Delete();
}

void PlottL(int isSM, double tBmH1, double tBmH2, double tBmH3, double q2, 
	    double thetaV=0.5, double chi=0, int isDgamma =0, int doScale=1){
  int nHis = 3;
  Styles style; style.setPadsStyle(2); style.applyStyle();
  BToDtaunu Dtaunu; BToDstaunu Dstaunu;

  double legW = 0.2, legH = 0.27;
  double legX = style.PadLeftMargin+0.37, legY = 1-style.PadTopMargin-0.52;

  int nBins = 1000, color[] = {28,4,2,1};
  double limTL[] = {-1, 1}, ml[] = {Dtaunu.mE,Dtaunu.mMu,Dtaunu.mTau,Dtaunu.mTau}, intTL[2][4], meanTL[2][4];
  TString yTitle[] = {"d#Gamma(D#tau#nu)/dcos(#theta_{#tau}) (10^{16}GeV^{-1})", 
		      "d#Gamma(D*#tau#nu)/(dcos(#theta_{#tau})d#theta_{V}d#chi) (10^{16}GeV^{-1})"}, hName;
  TString xTitle = "cos(#theta_{#tau})", legName[] = {"De", "Dmu", "Dtau", "#tau"};
  TString decName[] = {"Dtaunu ","D*taunu"};
  double tL, tBmH[] = {0, tBmH1, tBmH2, tBmH3}, maxH[] = {0,0};
  TCanvas can; can.Divide(2,1);
  TH1F *hTL[2][4];
  if(isSM==0) {
    ml[0] = ml[2]; ml[1] = ml[2]; 
    color[0] = 2; color[2] = 28;
    nHis = 4;
  }
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(style.LabelSize); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont); leg.SetBorderSize(0);
  for(int isDs=0; isDs<=1; isDs++){
    can.cd(isDs+1);
    for(int histo=0; histo<nHis; histo++){
      hName = "hTL"; hName += isDs; hName += histo;
      hTL[isDs][histo] = new TH1F(hName,"",nBins,limTL[0], limTL[1]);
      hTL[isDs][histo]->SetLineWidth(1);
      hTL[isDs][histo]->SetLineColor(color[histo]);
      if(isSM==0) {
	if(isDs==0) Dtaunu._gSR = -Dtaunu.mb_quark*pow(tBmH[histo],2);
	else        Dstaunu._gSR = -Dstaunu.mb_quark*pow(tBmH[histo],2);
	legName[histo] = "t#beta/m_{H} = "; legName[histo] += RoundNumber(tBmH[histo],1);
	if(fabs(tBmH[histo])<1e-6) legName[histo] = "SM"; 
	if(ml[histo]==Dtaunu.mE) legName[histo] += ", e";
	if(ml[histo]==Dtaunu.mMu) legName[histo] += ", #mu";
      }
      for(int bin=1; bin<=nBins; bin++){
	tL = hTL[isDs][histo]->GetBinCenter(bin);
	if(q2>0){
	  if(isDs==0) hTL[isDs][histo]->SetBinContent(bin,1e16*Dtaunu.Compute(q2,tL,1,ml[histo]));
	  else hTL[isDs][histo]->SetBinContent(bin,1e16*Dstaunu.Compute(q2,cos(tL),cos(thetaV),chi,isDgamma,true,1,ml[histo]));
	}else{
	  yTitle[isDs].ReplaceAll("d#theta_{V}d#chi","");
	  double q2min = ml[histo]*ml[histo], q2max = pow(Dtaunu._mB-Dtaunu._mD,2);
	  if(isDs) q2max = pow(Dstaunu._mB-Dstaunu._mDs,2);
	  int nPoints = 1000;
	  double dq2 = (q2max-q2min)/(double)nPoints, IntegratedTL=0, valq2=q2min, valRate;
	  //if(bin==1)cout<<"q2min "<<q2min<<", q2max "<<q2max<<", dq2 "<<dq2<<endl;
	  int isCLN=1; //if(histo==1) isCLN=0;
	  //int isBm=1;  if(histo==2) isBm=0;
	  //Dtaunu.SetMasses(isBm); Dstaunu.SetMasses(isBm);
	  for(int iq2=0; iq2<=nPoints; iq2++){
	    if(isDs==0) valRate = Dtaunu.Compute(valq2,acos(tL),isCLN,ml[histo]);
	    else        valRate = Dstaunu.Compute(valq2,tL,cos(thetaV),chi,isDgamma,true,isCLN,ml[histo]);
	    if(iq2==0 || iq2==nPoints) valRate /= 2.;
	    IntegratedTL += valRate;
	    valq2 += dq2; 
	  }
	  hTL[isDs][histo]->SetBinContent(bin,IntegratedTL);
	}
      }
      meanTL[isDs][histo] = hTL[isDs][histo]->GetMean();
      intTL[isDs][histo] = hTL[isDs][histo]->Integral();
      if(isDs==0) leg.AddEntry(hTL[isDs][histo], legName[histo]);
      if(doScale) {
	hTL[isDs][histo]->Scale(100/hTL[isDs][histo]->Integral());
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
	TString sTheta = "q^{2} = "; 
	if(q2>0)sTheta += RoundNumber(q2,1); else sTheta += "All";
	if(isDs) {
	  sTheta += ", #theta_{V} = "; sTheta += RoundNumber(thetaV,1);
	  sTheta += ", #chi = "; sTheta += RoundNumber(chi,1);
	}
	hTL[isDs][histo]->SetMinimum(0); hTL[isDs][histo]->SetMaximum(maxH[isDs]*1.17);
	hTL[isDs][histo]->Draw("c");
	style.setTitles(hTL[isDs][histo],xTitle,yTitle[isDs],sTheta);
      } else hTL[isDs][histo]->Draw("c same");
    }
    if(isDs==1) leg.Draw();

    cout<<endl<<decName[isDs]<<" ->   mean tL:  ";
    for(int histo=0; histo<nHis; histo++) cout<<RoundNumber(meanTL[isDs][histo],2)<<"\t";
    cout<<"  Rate/RateSM:  ";
    for(int histo=0; histo<nHis; histo++) cout<<RoundNumber(intTL[isDs][histo],2,intTL[isDs][0])<<"\t";

  }
  cout<<endl<<endl;
  TString epsName = "public_html/RD_tL_"; epsName += doScale; epsName += ".eps";
  can.SaveAs(epsName);
  for(int isDs=0; isDs<=1; isDs++)
    for(int histo=0; histo<nHis; histo++) 
      hTL[isDs][histo]->Delete();
}

void PlotQ2(double tBmH1, double tBmH2, double tBmH3, double theta, 
	    double ctv=0.5, double chi=0, int isDgamma =0, int doScale=1){
  int isSM = 0, nHis =4;
  Styles style; style.setPadsStyle(2); style.applyStyle();
  BToDtaunu Dtaunu; BToDstaunu Dstaunu;

  double legW = 0.2, legH = 0.27;
  double legX = style.PadLeftMargin+0.27, legY = 1-style.PadTopMargin-0.52;

  int nBins = 1000, color[] = {2,4,28,1};
  double limQ2[] = {0, 12.}, ml[] = {Dtaunu.mTau,Dtaunu.mTau,Dtaunu.mTau,Dtaunu.mTau}, intQ2[2][4], meanQ2[2][4];
  TString yTitle[] = {"d#Gamma(D#tau#nu)/dq^{2} (10^{16}GeV^{-1})", 
		      "d#Gamma(D*#tau#nu)/dq^{2} (10^{16}GeV^{-1})"}, hName;
  TString xTitle = "q^{2} (GeV^{2})", legName[] = {"e", "#mu", "#tau", "#tau"};
  TString decName[] = {"Dtaunu ","D*taunu"};
  double q2, tBmH[] = {0, tBmH1, tBmH2, tBmH3}, maxH[] = {0,0};
  TCanvas can; can.Divide(2,1);
  TH1F *hQ2[2][4];
  if(isSM==0) {
    limQ2[0] = 3;
    ml[0] = ml[2]; ml[1] = ml[2]; 
  }
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(style.LabelSize); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont); leg.SetBorderSize(0);
  for(int isDs=0; isDs<=1; isDs++){
    can.cd(isDs+1);
    for(int histo=0; histo<nHis; histo++){
      hName = "hQ2"; hName += isDs; hName += histo;
      hQ2[isDs][histo] = new TH1F(hName,"",nBins,limQ2[0], limQ2[1]);
      hQ2[isDs][histo]->SetLineWidth(1);
      hQ2[isDs][histo]->SetLineColor(color[histo]);
      if(isSM==0) {
	if(isDs==0) Dtaunu._gSR = -Dtaunu.mb_quark*pow(tBmH[histo],2);
	else        Dstaunu._gSR = -Dstaunu.mb_quark*pow(tBmH[histo],2);
	legName[histo] = "t#beta/m_{H} = "; legName[histo] += RoundNumber(tBmH[histo],1);
	if(fabs(tBmH[histo])<1e-6) legName[histo] = "SM"; 
      }
      for(int bin=1; bin<=nBins; bin++){
	q2 = hQ2[isDs][histo]->GetBinCenter(bin);
	if(isDs==0) hQ2[isDs][histo]->SetBinContent(bin,1e16*Dtaunu.Compute(q2,theta,1,ml[histo]));
	else hQ2[isDs][histo]->SetBinContent(bin,1e16*Dstaunu.Compute(q2,cos(theta),ctv,chi,isDgamma,true,1,ml[histo]));
	//if(histo==1 && isDs==0) hQ2[isDs][histo]->SetBinContent(bin,1e16*RDs.GammaD_q2(q2,ml[histo]));
      }
      meanQ2[isDs][histo] = hQ2[isDs][histo]->GetMean();
      intQ2[isDs][histo] = hQ2[isDs][histo]->Integral();
      if(isDs==0) leg.AddEntry(hQ2[isDs][histo], legName[histo]);
      if(doScale) {
	hQ2[isDs][histo]->Scale(100/hQ2[isDs][histo]->Integral());
	yTitle[isDs].ReplaceAll(" (10^{16}GeV^{-1})","");
      }
      if(maxH[isDs]<hQ2[isDs][histo]->GetMaximum()) maxH[isDs]=hQ2[isDs][histo]->GetMaximum();
    }
  }
  cout<<endl<<"tanBeta/mH: \t";
  for(int histo=0; histo<nHis; histo++)cout<<RoundNumber(tBmH[histo],1)<<"\t";
  for(int isDs=0; isDs<=1; isDs++){
    can.cd(isDs+1);
    for(int histo=0; histo<nHis; histo++){
      if(histo==0){
	TString sTheta = "#theta_{L} = "; sTheta += RoundNumber(theta,1);
	if(isDs) {
	  sTheta += ", cos(#theta_{V}) = "; sTheta += RoundNumber(ctv,1);
	  sTheta += ", #chi = "; sTheta += RoundNumber(chi,1);
	}
	hQ2[isDs][histo]->SetMaximum(maxH[isDs]*1.15);
	hQ2[isDs][histo]->Draw("c");
	style.setTitles(hQ2[isDs][histo],xTitle,yTitle[isDs],sTheta);
      } else hQ2[isDs][histo]->Draw("c same");
    }
    if(isDs==1) leg.Draw();

    cout<<endl<<decName[isDs]<<" ->   mean q2:  ";
    for(int histo=0; histo<nHis; histo++) cout<<RoundNumber(meanQ2[isDs][histo],2)<<"\t";
    cout<<"  Rate/RateSM:  ";
    for(int histo=0; histo<nHis; histo++) cout<<RoundNumber(intQ2[isDs][histo],2,intQ2[isDs][0])<<"\t";

  }
  cout<<endl<<endl;
  TString epsName = "public_html/RD_Q2tL_"; epsName += doScale; epsName += ".eps";
  can.SaveAs(epsName);
  for(int isDs=0; isDs<=1; isDs++)
    for(int histo=0; histo<nHis; histo++) 
      hQ2[isDs][histo]->Delete();
}



void PlotEffi(int whichB = 2, int dHig = 5){
  Styles style; style.setPadsStyle(2); style.applyStyle();

  TCanvas can; can.Divide(2,1);
  TH1F *Histo[2][2];
//   TString yTitle[] = {"#epsilon(D#tau#nu)/#epsilon_{SM} (%)", 
// 		      "#epsilon(D*#tau#nu)/#epsilon_{SM} (%)"}, HigName;
  TString yTitle[] = {"t0", "t1"}, HigName;

  TSpline3 *SplineRD[2], *SplineErrRD[2];
  double RD[2][130], errorRD[2][130], List_tBmH[130], RelerrRD[2][130];
  float weight;
  int candType, MCType;
  double totW[2][6], totW2[2][6], totN[2][2][6], nSM[2]={1,1}, errorSM[2]={1,1};
  int iniHig = 0, finHig = 100, file=0, nFiles = (finHig-iniHig)/dHig+1, nBins = 1000;
  for(int his=0; his<2; his++) {
    for(int isDs=0; isDs<2; isDs++) {
      HigName = "Histo"; HigName += isDs; HigName += his;
      Histo[isDs][his] = new TH1F(HigName,"",nBins,0,1);
      Histo[isDs][his]->SetLineWidth(2); 
    }
  }
  for(int iHig=iniHig; iHig<=finHig; iHig += dHig){
    HigName = "AWG82/ntuples/small/FitRAllHigx"; 
    if(iHig<100) HigName += "0";
    if(iHig<10) HigName += "0";
    HigName += iHig; HigName += "_RunAll.root";
    TChain genChain("ntp1"); genChain.Add(HigName);
    genChain.SetBranchAddress("weight",&weight);
    genChain.SetBranchAddress("candType",&candType);
    genChain.SetBranchAddress("MCType",&MCType);
    List_tBmH[file] = (double)iHig/100.;

    for(int i=0; i<4; i++){
      for(int j=0; j<2; j++){
	totW[j][i] = 0; totW2[j][i] = 0; 
	for(int k=0; k<2; k++) totN[k][j][i] = 0;
      }
    }
    for(int entry=0; entry<genChain.GetEntries(); entry++){
      genChain.GetEvent(entry);
      if(candType<=2 && MCType==5) {
	totW[0][0] += weight; totW2[0][0] += weight*weight;
	if(candType==1) totN[0][0][0] += weight;
	else totN[1][0][0] += weight;
      }
      if(candType<=2 && MCType==6) {
	totW[0][1] += weight; totW2[0][1] += weight*weight;
	if(candType==2) totN[0][0][1] += weight;
	else totN[1][0][1] += weight;
      }
      if(candType>=3 && candType<=4 && MCType==11){
	totW[0][2] += weight; totW2[0][2] += weight*weight;
	if(candType==3) totN[0][0][2] += weight;
	else totN[1][0][2] += weight;
      }
      if(candType>=3 && candType<=4 && MCType==12){
	totW[0][3] += weight; totW2[0][3] += weight*weight;
	if(candType==4) totN[0][0][3] += weight;
	else totN[1][0][3] += weight;
      }

    }
    double DiffError;
    for(int isDs=0; isDs<2; isDs++) {
      if(iHig==0){
	if(whichB==0){
	  nSM[isDs] = totW[0][isDs]; 
	  errorSM[isDs] = totW2[0][isDs]; 
	} else if(whichB==1){
	  nSM[isDs] = totW[0][isDs+2]; 
	  errorSM[isDs] = totW2[0][isDs+2]; 
	} else {
	  nSM[isDs] = totW[0][isDs]+totW[0][isDs+2];
	  errorSM[isDs] = totW2[0][isDs]+totW2[0][isDs+2]; 
	}
      }
      if(whichB==0){
	RD[isDs][file] = totW[0][isDs]/nSM[isDs]*100;
	DiffError = totW2[0][isDs]/pow(totW[0][isDs],2)-errorSM[isDs]/pow(nSM[isDs],2);
	RelerrRD[isDs][file] = sqrt(totW2[0][isDs]/pow(totW[0][isDs],2));
      } else if(whichB==1){
	RD[isDs][file] = totW[0][isDs+2]/nSM[isDs]*100;
	DiffError = totW2[0][isDs+2]/pow(totW[0][isDs+2],2)-errorSM[isDs]/pow(nSM[isDs],2);
	RelerrRD[isDs][file] = sqrt(totW2[0][isDs+2]/pow(totW[0][isDs+2],2));
      } else {
	RD[isDs][file] = (totW[0][isDs]+totW[0][isDs+2])/nSM[isDs]*100;
	DiffError = (totW2[0][isDs]+totW2[0][isDs+2])/(pow(totW[0][isDs]+totW[0][isDs+2],2))
	  -errorSM[isDs]/pow(nSM[isDs],2);
	RelerrRD[isDs][file] = sqrt((totW2[0][isDs]+totW2[0][isDs+2])/
				    (pow(totW[0][isDs]+totW[0][isDs+2],2)));
      }
      if(DiffError<0) DiffError=0;
      errorRD[isDs][file] = sqrt(DiffError)*RD[isDs][file];
      //cout<<RoundNumber(List_tBmH[file],2)<<": Y/Y_SM = "<<RoundNumber(RD[isDs][file],2)
      //<<" +- "<<RoundNumber(errorRD[isDs][file],2)<<" \t";
    }
    //cout<<endl;
    file++;
  }
  for(int isDs=0; isDs<2; isDs++) {
    can.cd(isDs+1);
    HigName = "Spline"; HigName += isDs; 
    SplineRD[isDs] = new TSpline3(HigName,List_tBmH,RD[isDs],nFiles,"ble1",0); 
    HigName += "Error";
    SplineErrRD[isDs] = new TSpline3(HigName,List_tBmH,errorRD[isDs],nFiles,"ble1",0);
    //cout<<SplineRD[isDs]->Eval(0.01)<<" +- "<<SplineErrRD[isDs]->Eval(0.01)<<endl;
    for(int bin=1; bin<=nBins; bin++){
      double tBmH = Histo[isDs][0]->GetBinCenter(bin), valEffi = SplineRD[isDs]->Eval(tBmH);
      double errorEffi = SplineErrRD[isDs]->Eval(tBmH);
      for(int his=0; his<2; his++) {
	Histo[isDs][his]->SetBinContent(bin, valEffi);
	if(his==0) Histo[isDs][his]->SetBinError(bin, errorEffi);
	else       Histo[isDs][his]->SetBinError(bin,0);
      }
    }
    Histo[isDs][0]->SetLineColor(2); Histo[isDs][0]->SetFillColor(2);Histo[isDs][0]->SetFillStyle(3002);
    Histo[isDs][0]->Draw("e3");
    style.setTitles(Histo[isDs][0],"tan#beta/m_{H} (GeV^{-1})",yTitle[isDs]); 
    Histo[isDs][1]->Draw("c same");    
  }
  TString epsName = "public_html/Higgs_Efficiency"; epsName+=whichB; epsName+=".eps";
  can.SaveAs(epsName);

  if(dHig==5){
    cout<<"double dEffError[2][] = {";
    for(int isDs=0; isDs<2; isDs++){
      cout<<"{";
      for(file=0; file<nFiles; file++){
	if(file<nFiles-1) cout<<RoundNumber(errorRD[isDs][file]*100,2,RD[isDs][file])<<", ";
	else   cout<<RoundNumber(errorRD[isDs][file]*100,2,RD[isDs][file]);
      }
      if(isDs==0) cout<<"}, "<<endl;
      else cout<<"}";
    }
    cout<<"};"<<endl<<endl;


    cout<<"double RelError[2][] = {";
    for(int isDs=0; isDs<2; isDs++){
      cout<<"{";
      for(file=0; file<nFiles; file++){
	if(file<nFiles-1) cout<<RoundNumber(RelerrRD[isDs][file]*100,2)<<", ";
	else   cout<<RoundNumber(RelerrRD[isDs][file]*100,2);
      }
      if(isDs==0) cout<<"}, "<<endl;
      else cout<<"}";
    }
    cout<<"};"<<endl;


  }
  for(int isDs=0; isDs<2; isDs++){
    SplineRD[isDs]->Delete();
    SplineErrRD[isDs]->Delete();
    for(int his=0; his<2; his++) 
      if(Histo[isDs][his]) Histo[isDs][his]->Delete();
  }
}


void PlotYields(){
  Styles style; style.setPadsStyle(2); style.applyStyle();

  TCanvas can; can.Divide(2,1);
  //TString yTitle[] = {"Yield(B#rightarrowD#tau#nu)", "Yield(B#rightarrowD*#tau#nu)"}, HigName;
  TString yTitle[] = {"Y1", "Y2"}, HigName;

  TSpline3 *MeasuredRD[4];
  int iniHig = 0, iHig = iniHig, finHig = 100, dHig = 5, bin, nFiles = (finHig-iniHig)/dHig+1, nBins = 1000;
  double higRD[4][30], List_tBmH[30];
  TString folder = "FitAll/fits/TextFinalIsoDataHi2x", text = ""; 
  for(int file=0; file<nFiles; file++){
    TString fileName = folder; 
    if(iHig<100) fileName += "0";
    if(iHig<10) fileName += "0";
    fileName += iHig; fileName += ".txt";
    fstream textFile; textFile.open(fileName,fstream::in);
    bin = 0;
    while(!text.Contains("R(D")) {
      textFile>>text; 
      bin++;
      if(bin>1000) {cout<<"R(D) not found in "<<fileName<<endl; return;}
    }
    for(int cand=0; cand<2; cand++){
      if(cand>0) textFile>>text;
      textFile>>text; higRD[cand][file] = text.Atof();
      textFile>>text; textFile>>text; text.ReplaceAll(",",""); higRD[cand+2][file] = text.Atof();
      for(int i=0; i<5; i++) textFile>>text; 
      textFile>>text; higRD[cand][file] = text.Atof();
      textFile>>text; textFile>>text; text.ReplaceAll(",",""); higRD[cand+2][file] = text.Atof();
      for(int i=0; i<4; i++) textFile>>text; 
      //cout<<higRD[cand][file]<<" #pm "<<higRD[cand+2][file]<<" - a "
      //<<RoundNumber(higRD[cand+2][file]*100,1,higRD[cand][file])<<"% error"<<endl;
      List_tBmH[file] = (double)iHig/100.;
    }
    //cout<<endl;
    iHig += dHig;
  }
  for(int isDs=0; isDs<4; isDs++){
    text = "SplineRD"; text += isDs;
    MeasuredRD[isDs] = new TSpline3(text,List_tBmH,higRD[isDs],nFiles,"ble1",0); 
  }

  TH1F *Histo[2][2];
  for(int his=0; his<2; his++) {
    for(int isDs=0; isDs<2; isDs++) {
      HigName = "Histo"; HigName += isDs; HigName += his;
      Histo[isDs][his] = new TH1F(HigName,"",nBins,0,1);
      Histo[isDs][his]->SetLineWidth(2); 
    }
  }

  for(int isDs=0; isDs<2; isDs++) {
    can.cd(isDs+1);
    for(int bin=1; bin<=nBins; bin++){
      double tBmH = Histo[isDs][0]->GetBinCenter(bin), valEffi = MeasuredRD[isDs]->Eval(tBmH);
      double errorEffi = MeasuredRD[isDs+2]->Eval(tBmH);
      for(int his=0; his<2; his++) {
	Histo[isDs][his]->SetBinContent(bin, valEffi);
	if(his==0) Histo[isDs][his]->SetBinError(bin, errorEffi);
	else       Histo[isDs][his]->SetBinError(bin,0);
      }
    }
    Histo[isDs][0]->SetLineColor(4); Histo[isDs][0]->SetFillColor(4);Histo[isDs][0]->SetFillStyle(3002);
    Histo[isDs][0]->Draw("e3");
    style.setTitles(Histo[isDs][0],"tan#beta/m_{H} (GeV^{-1})",yTitle[isDs]); 
    Histo[isDs][1]->SetLineColor(4); 
    Histo[isDs][1]->Draw("c same");    
  }
  TString epsName = "public_html/Higgs_Yields.eps";
  can.SaveAs(epsName);

  for(int isDs=0; isDs<2; isDs++){
    MeasuredRD[isDs]->Delete();
    MeasuredRD[isDs+2]->Delete();
    for(int his=0; his<2; his++) 
      if(Histo[isDs][his]) Histo[isDs][his]->Delete();
  }
}



void PlotRate(int isDs = 1, double gsmax = 4){

  RateCalc RDs;
  RDs.ReadFF("babar_code/NewPhysics/FFinputs_HFAG11");

  Styles style; 
  style.setPadsStyle(1); style.CanvasW = 350; style.CanvasH = 225; 
  style.TitleSize = 0.08; style.xTitleOffset = 0.898; style.yTitleOffset = 0.6; 
  style.PadBottomMargin = 0.18;
  style.applyStyle();
  double legW = 0.38, legH = 0.17;
  double legX = style.PadLeftMargin+0.01, legY = 1-style.PadTopMargin-0.01;
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(style.LabelSize*0.94); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont); leg.SetBorderSize(0);

  int nBins = 300;
  double minY[] = {0, 0.};
  //TString paper[] = {"Eq. 9: Kamenik, Mescia (2008)", 
  TString paper[] = {"Eq. 4: Fajfer, Kamenik, Nisandzic (2012)", 
		     "Eq. 11: Fajfer, Kamenik, Nisandzic (2012)"};
  TString yTitle[] = {"R(D)", "R(D*)"};
  TString xTitle[] = {"-(g_{SR}+g_{SL}) (GeV^{-1})", "-(g_{SR}-g_{SL}) (GeV^{-1})"};
  double gs = 0, dgs = gsmax/(double)nBins, PRDrate;
  TCanvas can; can.cd();
  TH1F hRDs("RDs","",nBins,0,gsmax);
  TH1F hRDs_PRD("RDs_PRD","",nBins,0,gsmax); hRDs_PRD.SetLineColor(2);
  for(int bin=1; bin<=nBins; bin++){
    RDs.gsR = gs;
    hRDs.SetBinContent(bin,RDs.RRate(isDs));
    if(isDs) PRDrate = 0.252*(1 + 0.12*RDs.mTau*gs + 0.044*pow(RDs.mTau*gs,2));
    else     PRDrate = 0.296*(1 + 1.5 *RDs.mTau*gs + 1.1  *pow(RDs.mTau*gs,2));
    hRDs_PRD.SetBinContent(bin, PRDrate);
    gs -= dgs;
  }
  style.setTitles(&hRDs,xTitle[isDs],yTitle[isDs]); 
  hRDs.SetMinimum(minY[isDs]);
  hRDs.Draw("c");
  hRDs_PRD.Draw("c same");
  leg.AddEntry(&hRDs_PRD, paper[isDs]);
  leg.AddEntry(&hRDs,"Our calculation");
  leg.Draw();
  TString epsName = "public_html/RD_Implementation_"; epsName += isDs; epsName += ".eps";
  can.SaveAs(epsName);
}


void PlottL(int isSM, double tBmH1=0.3, double tBmH2=0.5, double tBmH3=1, int doScale=1){
  int nHis = 3;
  Styles style; style.setPadsStyle(2); style.applyStyle();
  BToDtaunu Dtaunu; BToDstaunu Dstaunu;


  int nBins = 1000, color[] = {28,4,2,1};
  double limTL[] = {-1, 1}, ml[] = {Dtaunu.mE,Dtaunu.mMu,Dtaunu.mTau,Dtaunu.mTau}, intTL[2][4], meanTL[2][4];
  TString yTitle[] = {"D0", "D1"}, hName;
//   TString yTitle[] = {"d#Gamma(Dl#nu)/dcos(#theta_{#tau}) (10^{16}GeV^{-1})", 
// 		      "d#Gamma(Dl#nu)/(dcos(#theta_{#tau})) (10^{16}GeV^{-1})"}, hName;
  TString xTitle = "cos(#theta_{#tau})", legName[] = {"De", "Dmu", "Dtau", "#tau"};
  TString decName[] = {"Dtaunu ","D*taunu"}, PanelTag[] = {"a)", "b)"};
  double tL, tBmH[] = {0, tBmH1, tBmH2, tBmH3}, maxH[] = {0,0};
  TCanvas can; can.Divide(2,1);
  TH1F *hTL[2][4];
  if(isSM==0) {
    ml[0] = ml[2]; ml[1] = ml[2]; 
    color[0] = 2; color[2] = 28;
    nHis = 4;
  }
  double legW = 0.23, legH = 0.27;
  double legX = 1-style.PadRightMargin-0.01, legY = 1-style.PadTopMargin-0.01;
  TLegend leg(legX-legW, legY-legH, legX, legY);
  leg.SetTextSize(style.LabelSize); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont); leg.SetBorderSize(0);
  for(int isDs=0; isDs<=1; isDs++){
    can.cd(isDs+1);
    for(int histo=0; histo<nHis; histo++){
      hName = "hTL"; hName += isDs; hName += histo;
      hTL[isDs][histo] = new TH1F(hName,"",nBins,limTL[0], limTL[1]);
      hTL[isDs][histo]->SetLineWidth(1);
      hTL[isDs][histo]->SetLineColor(color[histo]);
      if(isSM==0) {
	if(isDs==0) Dtaunu._gSR = -Dtaunu.mb_quark*pow(tBmH[histo],2);
	else        Dstaunu._gSR = -Dstaunu.mb_quark*pow(tBmH[histo],2);
	legName[histo] = "t#beta/m_{H} = "; legName[histo] += RoundNumber(tBmH[histo],1);
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
  TString epsName = "public_html/RD_tL_"; epsName += doScale; epsName += ".eps";
  can.SaveAs(epsName);
  for(int isDs=0; isDs<=1; isDs++)
    for(int histo=0; histo<nHis; histo++) 
      hTL[isDs][histo]->Delete();
}


void PlotQ2(int isSM = 1, int doScale=0, double tBmH1=0.3, double tBmH2=0.5, double tBmH3=1){
  int nHis = 3;
  Styles style; style.setPadsStyle(2); style.applyStyle();
  BToDtaunu Dtaunu; BToDstaunu Dstaunu;
  RateCalc RDs;
  RDs.ReadFF("babar_code/NewPhysics/FFinputs_HFAG11");


  int nBins = 1000, color[] = {28,4,2,1};
  double limQ2[] = {0, 12.}, ml[] = {Dtaunu.mE,Dtaunu.mMu,Dtaunu.mTau,Dtaunu.mTau};
//   TString yTitle[] = {"d#Gamma(Dlnu)/dq^{2} (10^{16}GeV^{-1})", 
// 		      "d#Gamma(D*lnu)/dq^{2} (10^{16}GeV^{-1})"}, hName, PanelTag[] = {"a)", "b)"};
  TString yTitle[] = {"D0", "D1"}, hName, PanelTag[] = {"a)", "b)"};
  TString xTitle = "q^{2} (GeV^{2})", legName[] = {"De", "Dmu", "Dtau", "#tau"};
  double q2, tBmH[] = {0, tBmH1, tBmH2, tBmH3}, maxH[] = {0,0};
  TCanvas can; can.Divide(2,1);
  TH1F *hQ2[2][4];
  if(isSM==0) {
    limQ2[0] = 3;
    ml[0] = ml[2]; ml[1] = ml[2]; 
  }
  double legW = 0.23, legH = 0.27;
  double legX = 1-style.PadRightMargin-0.01, legY = 1-style.PadTopMargin-0.01;
  TLegend leg(legX-legW, legY-legH, legX, legY);
  leg.SetTextSize(style.LabelSize); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont); leg.SetBorderSize(0);
  for(int isDs=0; isDs<=1; isDs++){
    can.cd(isDs+1);
    for(int histo=0; histo<nHis; histo++){
      hName = "hQ2"; hName += isDs; hName += histo;
      hQ2[isDs][histo] = new TH1F(hName,"",nBins,limQ2[0], limQ2[1]);
      hQ2[isDs][histo]->SetLineWidth(1);
      hQ2[isDs][histo]->SetLineColor(color[histo]);
      if(isSM==0) {
	if(isDs==0) Dtaunu._gSR = -Dtaunu.mb_quark*pow(tBmH[histo],2);
	else        Dstaunu._gSR = -Dstaunu.mb_quark*pow(tBmH[histo],2);
	legName[histo] = "t#beta/m_{H} = "; legName[histo] += RoundNumber(tBmH[histo],1);
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
      }
      if(isDs==0) leg.AddEntry(hQ2[isDs][histo], legName[histo]);
      if(doScale) {
	hQ2[isDs][histo]->Scale(100/hQ2[isDs][histo]->Integral());
	yTitle[isDs].ReplaceAll(" (10^{16}GeV^{-1})","");
      }
      if(maxH[isDs]<hQ2[isDs][histo]->GetMaximum()) maxH[isDs]=hQ2[isDs][histo]->GetMaximum();
    }
  }
  for(int isDs=0; isDs<=1; isDs++){
    can.cd(isDs+1);
    for(int histo=0; histo<nHis; histo++){
      if(histo==0){
	hQ2[isDs][histo]->SetMaximum(maxH[isDs]*1.1);
	hQ2[isDs][histo]->Draw("c");
	style.setTitles(hQ2[isDs][histo],xTitle,yTitle[isDs], PanelTag[isDs]);
      } else hQ2[isDs][histo]->Draw("c same");
    }
    if(isDs==0) leg.Draw();
  }
  TString epsName = "public_html/RD_Q2_"; epsName += doScale; epsName += ".eps";
  can.SaveAs(epsName);
  for(int isDs=0; isDs<=1; isDs++)
    for(int histo=0; histo<nHis; histo++) 
      hQ2[isDs][histo]->Delete();
}

