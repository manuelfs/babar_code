#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TBox.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TChain.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFile.h"
#include "TMath.h"

#include "styles.hpp"
#include "keys_utils.hpp" 
#include <fstream>
#include <iostream>

using namespace std;

void CalcRatio(double BTau[3], double Bl[2], double RD[3]);
void averageN(double RD[][3], int nRD, double rho, double averRD[]);
void PrintRD(double RD[3], int isDs, int digits=3);
double totError(double RD[]);

int main(){
  int nFiles = 5, AverAll = 0, Verbose = 1;
  
  if(nFiles<5) AverAll = 0;
  styles style; style.setPadsStyle(-8); 
  style.CanvasH = 360; style.PadTopMargin = 0; style.PadLeftMargin = 0; style.PadRightMargin = 0.02;
  style.LabelSize *= 0.92;
  style.setDefaultStyle();

  // Setting values of measured RD(*)
  double minRD[2] = {0.16, 0.21}, maxRD[2] = {0.85, 0.62}, nomiRD[2][2]; 
  double SMRD[2][2] = {{0.297, 0.017}, {0.252, 0.003}}, maxHisto = 1., frac = maxHisto/40., frac2 = 0.8*frac;
  double BelleRho = 1, BaBeRho = 0.25;
  double B0Dslnu[] = {4.59, 0.26}, BmDslnu[2], BmDlnu[] = {2.26, 0.11}; // Dslnu from Belle, Dlnu from PDG2010
  double RD[2][5][3];
  TString SepLine = "==========================================================================";
  for(int iE=0; iE<2; iE++) BmDslnu[iE] = B0Dslnu[iE]*1.638e-12/1.525e-12; // Ratio of B0/B+ lifetimes

  // Belle 2007
  double BF_B0DsmTauNu_07[] = {2.02, 0.39, 0.37}; // BF(B0->D*-TauNu)
  CalcRatio(BF_B0DsmTauNu_07, B0Dslnu, RD[1][0]);

  // Belle 2009
  double RD_09[2][3]  = {{0.70, 0.19, 0.10}, {0.48, 0.21, 0.06}};  // R(D0)  and R(D+)
  double RDs_09[2][3] = {{0.47, 0.11, 0.07}, {0.48, 0.13, 0.05}};  // R(D*0) and R(D*+)
  averageN(RD_09,  2, BelleRho, RD[0][1]);
  averageN(RDs_09, 2, BelleRho, RD[1][1]);

  // Belle 2010
  double BF_BmD0TauNu_10[]  = {0.77, 0.22, 0.12}; // BF(B-->D0TauNu)
  double BF_BmDs0TauNu_10[] = {2.12, 0.28, 0.29}; // BF(B-->D*0TauNu)
  CalcRatio(BF_BmD0TauNu_10,  BmDlnu,  RD[0][2]);
  CalcRatio(BF_BmDs0TauNu_10, BmDslnu, RD[1][2]);

  // BaBar 2008
  RD[0][3][0] = 0.416;  RD[0][3][1] = 0.117;  RD[0][3][2] = 0.052; 
  RD[1][3][0] = 0.297;  RD[1][3][1] = 0.056;  RD[1][3][2] = 0.018; 

  // BaBar 2012
  RD[0][4][0] = 0.4402; RD[0][4][1] = 0.0577; RD[0][4][2] = 0.0423; 
  RD[1][4][0] = 0.3316; RD[1][4][1] = 0.0236; RD[1][4][2] = 0.0183; 

  // Belle averages
  double AverBelle[2][2][3];
  averageN(RD[0]+1, 2, BelleRho, AverBelle[0][0]);
  averageN(RD[1],   3, BelleRho, AverBelle[1][0]);
  if(Verbose){
    cout<<endl<<"Belle measurements"<<endl<<SepLine<<endl;
    for(int mea=0; mea<3; mea++)
      for(int isDs=0; isDs<2; isDs++) PrintRD(RD[isDs][mea], isDs, 2);
    cout<<endl<<"BaBar measurement"<<endl<<SepLine<<endl;
    for(int isDs=0; isDs<2; isDs++) PrintRD(RD[isDs][3], isDs, 2);
    cout<<endl<<"Belle averages"<<endl<<SepLine<<endl;
    for(int isDs=0; isDs<2; isDs++) PrintRD(AverBelle[isDs][0], isDs);
  }

  // BaBe averages
  double AverBaBe[2][2][3];
  if(Verbose) cout<<endl<<"BaBe 2008 averages"<<endl<<SepLine<<endl;
  for(int isDs=0; isDs<2; isDs++) {
    for(int iRD=0; iRD<3; iRD++) AverBelle[isDs][1][iRD] = RD[isDs][3][iRD];
    averageN(AverBelle[isDs], 2, BaBeRho, AverBaBe[0][isDs]);
    if(Verbose) PrintRD(AverBaBe[0][isDs], isDs, 4);
  }
  if(Verbose) cout<<endl<<"BaBe 2012 averages"<<endl<<SepLine<<endl;
  for(int isDs=0; isDs<2; isDs++) {
    for(int iRD=0; iRD<3; iRD++) AverBelle[isDs][1][iRD] = RD[isDs][4][iRD];
    averageN(AverBelle[isDs], 2, BaBeRho, AverBaBe[1][isDs]);
    if(Verbose) PrintRD(AverBaBe[1][isDs], isDs, 4);
  }

  double nSigma[2];
  for(int isDs=0; isDs<2; isDs++) {
    nomiRD[isDs][0] = AverBaBe[AverAll][isDs][0]; 
    nomiRD[isDs][1] = totError(AverBaBe[AverAll][isDs]);
    nSigma[isDs] = pow(nomiRD[isDs][0] - SMRD[isDs][0],2)/(pow(nomiRD[isDs][1],2)+pow(SMRD[isDs][1],2));
  }
  double pValue = TMath::Prob(nSigma[0]+nSigma[1],2);
  cout<<endl<<"R(D) is "<<RoundNumber(sqrt(nSigma[0]),1)<<" sigma away, and R(D*) is "<<RoundNumber(sqrt(nSigma[1]),1)
      <<" sigma away. A total of "<< RoundNumber(sqrt(2)*TMath::ErfInverse(1 - pValue),1)<<" sigma."<<endl<<endl;


  int colors[] = {28,28,28,4,4, kGray+2, kGray, kRed+2, kRed-9}, iMeas[] = {0,3,1,2,4};
  TString fileTag[] = {"Belle 2007", "Belle 2009", "Belle 2010", "BaBar 2008", "BaBar 2012"};
  double LeftMargin = 0.21, padW = (1-LeftMargin)/2., padH = (1-style.PadBottomMargin)/static_cast<double>(nFiles);
  TCanvas can("can","RD results");
  TPad *Pads[2];
  TString xTitle[] = {"R(D)", "R(D*)"};
  TH1F* hRD[2];
  TBox box; box.SetLineColor(0);box.SetFillColor(34);
  TLatex label; label.SetTextSize(0.08); label.SetTextFont(style.nFont); label.SetTextAlign(12);
  TMarker mark; mark.SetMarkerStyle(8); mark.SetMarkerSize(1.1);
  TLine lin; lin.SetLineStyle(1); lin.SetLineWidth(2);
  for(int pad=0; pad<2; pad++){
    can.cd(0);
    double PadXY[2][2] = {{LeftMargin+padW*pad, LeftMargin+padW*(pad+1)},{0, 1}};
    TString hName = "Pad_"; hName += pad;
    Pads[pad] = new TPad(hName,"",PadXY[0][0], PadXY[1][0], PadXY[0][1], PadXY[1][1]);
    Pads[pad]->Draw(); Pads[pad]->cd();

    hName = "Ratio"; hName += pad+1;
    hRD[pad] = new TH1F(hName, "", 10, minRD[pad], maxRD[pad]);
    hRD[pad]->GetYaxis()->SetNdivisions(0);
    hRD[pad]->SetMaximum(maxHisto);
    hRD[pad]->GetXaxis()->CenterTitle(true);
    hRD[pad]->Draw();
    box.SetFillColor(colors[6]); box.DrawBox(nomiRD[pad][0]-nomiRD[pad][1], 0, nomiRD[pad][0]+nomiRD[pad][1], maxHisto);
    lin.SetLineColor(colors[5]); lin.DrawLine(nomiRD[pad][0], 0, nomiRD[pad][0], maxHisto);
    box.SetFillColor(colors[8]); box.DrawBox(SMRD[pad][0]-SMRD[pad][1], 0, SMRD[pad][0]+SMRD[pad][1], maxHisto);
    lin.SetLineColor(colors[7]); lin.DrawLine(SMRD[pad][0], 0, SMRD[pad][0], maxHisto);
    for(int ifile=0; ifile<nFiles; ifile++) {
      int file = iMeas[ifile];
      if(pad==0 && file==0) continue;
      int nBeg = 0; if(pad==1) nBeg = 1;
      double halfH = 1/static_cast<double>(nFiles)/2.;
      if(file==4){
	label.SetTextColor(2); 
      }
      for(int cand=nBeg; cand<4; cand+=4){
	double height = 1 - (halfH*(2*ifile+1) - halfH/3.2);
	mark.SetMarkerColor(colors[file]); lin.SetLineColor(colors[file]);
	if(file==4){
	  mark.SetMarkerColor(2); lin.SetLineColor(2);
	}
	mark.DrawMarker(RD[cand][file][0], height);
	double minTot = RD[cand][file][0]-totError(RD[cand][file]), maxTot = RD[cand][file][0]+totError(RD[cand][file]);
	double minSta = RD[cand][file][0]-RD[cand][file][1], maxSta = RD[cand][file][0]+RD[cand][file][1];
	lin.SetLineWidth(1);
	lin.DrawLine(minTot, height+frac, minTot, height-frac); lin.DrawLine(maxTot, height+frac, maxTot, height-frac);
	lin.DrawLine(minSta, height+frac2, minSta, height-frac2); lin.DrawLine(maxSta, height+frac2, maxSta, height-frac2);
	lin.SetLineWidth(2); lin.DrawLine(minTot, height, maxTot, height);
      }
    }
    hRD[pad]->Draw("axis same");
    style.setTitles(hRD[pad],xTitle[pad]);
    label.SetTextSize(0.07); label.SetTextColor(1);label.SetTextAlign(23);
  }

  can.cd(0);
  for(int ifile=0; ifile<nFiles; ifile++) {
    int file = iMeas[ifile];
    label.SetTextFont(22); label.SetTextColor(colors[file]); label.SetTextAlign(12);
    label.SetTextSize(0.08); label.DrawLatex(0.01, style.PadBottomMargin+0.02+padH*(nFiles-ifile-0.5), fileTag[file]);
  }
  TString pName = "plots/old_RDx.pdf"; 
  can.SaveAs(pName);
  for(int pad=0; pad<2; pad++){
      hRD[pad]->Delete();
  }

  return 0;
}


double totError(double RD[]){
  return sqrt(pow(RD[1],2)+pow(RD[2],2));
}

// Average of nRD measurements, all with a correlation of rho. 
// Uses Bob Kowalewski "Methodology and code used in HFAG semileptonic averages" (2005)
void averageN(double RD[][3], int nRD, double rho, double averRD[]){
  TMatrixT<double> StatE(nRD,nRD), SystE(nRD,nRD);
  
  for(int iR=0; iR<nRD; iR++){
    for(int iRD=0; iRD<3; iRD++) averRD[iRD] = 0;
    for(int jR=0; jR<nRD; jR++){
      if(iR==jR){
	StatE(iR, jR) = pow(RD[iR][1],2);
	SystE(iR, jR) = pow(RD[iR][2],2);
      } else {
	StatE(iR, jR) = 0;
	SystE(iR, jR) = rho*RD[iR][2]*RD[(iR+1)%nRD][2];
      }
    }
  }

  TMatrixT<double> TotE(StatE); TotE += SystE; TotE.Invert();
  TMatrixT<double> WStatW(nRD,nRD), WSystW(nRD,nRD);
  WStatW.Mult(TotE, StatE); WStatW *= TotE;
  WSystW.Mult(TotE, SystE); WSystW *= TotE;
  double sumW = 0;
  for(int iR=0; iR<nRD; iR++){
    for(int jR=0; jR<nRD; jR++){
      averRD[0] += TotE(iR,jR)*RD[iR][0];
      averRD[1] += WStatW(iR,jR);
      averRD[2] += WSystW(iR,jR);
      sumW += TotE(iR,jR);
    }
  }
  for(int iRD=0; iRD<3; iRD++){
    averRD[iRD] /= sumW;
    if(iRD>0) {
      averRD[iRD] /= sumW;
      averRD[iRD] = sqrt(averRD[iRD]);
    }
  }
}

void CalcRatio(double BTau[3], double Bl[2], double RD[3]){
  RD[0] = BTau[0]/Bl[0];
  RD[1] = BTau[1]/Bl[0];
  RD[2] = RD[0]*sqrt(pow(BTau[2]/BTau[0],2)+pow(Bl[1]/Bl[0],2));
}

void PrintRD(double RD[3], int isDs, int digits){
  TString RDNames[] = {"R(D)  = ", "R(D*) = "};
  cout<<RDNames[isDs]<<RoundNumber(RD[0],digits)<<" +- "<<RoundNumber(RD[1],digits)<<" +- "<<
    RoundNumber(RD[2],digits)<<" => Total error is "<<RoundNumber(totError(RD),digits)<<", a "<<
    RoundNumber(totError(RD)*100,1,RD[0])<<"% error"<<endl;
}

