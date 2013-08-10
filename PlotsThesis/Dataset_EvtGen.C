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

void Dataset_EvtGen(int isTheta=0){
  Styles style; style.setPadsStyle(2); 
  style.yTitleOffset = 1.4; style.PadLeftMargin = 0.182;
  style.applyStyle();
  BToDtaunu Dtaunu; BToDstaunu Dstaunu;

  double legW = 0.22, legH = 0.25;
  double legX = style.PadLeftMargin+0.03, legY = 1-style.PadTopMargin-0.01;
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(style.LabelSize); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont); leg.SetBorderSize(0);

  TChain mc("GqaMCAnalysis/ntp1");
  mc.Add("AWG82/generator/Archive/OneMEvtGen_BpDxtaunu.root");

  int nBins[] = {100,100,100}, color[] = {1,4,2,1}, nHis=3, lineWidth[] = {2,1,2}, isCLN=0, iLeg[] = {2,1,0};
  int digits[] = {2, 2};
  TString cuts[] = {"(MCType==5||MCType==11)","(MCType==6||MCType==12)"};
  double limQ2[] = {3, 12.}, ml[] = {Dtaunu.mTau,Dtaunu.mTau,Dtaunu.mTau,Dtaunu.mTau}, rmsQ2[2][4], meanQ2[2][4];
  TString yTitle[] = {" GeV^{2})", ")"}, hName;
  TString xTitle = "q^{2} (GeV^{2})", legName[] = {"MC", "HQET", "ISGW2"}, rightTag[] = {"a)", "b)"};;
  TString decName[] = {"Dtaunu ","D*taunu"}, varName[] = {"trueQ2","trueCTL"};
  double q2, maxH[] = {0,0}, factor[] = {1, 1};
  TCanvas can; can.Divide(2,1);
  TH1F *hQ2[2][3];
  if(isTheta){
    limQ2[0] = -1; limQ2[1] = 1; 
    xTitle = "cos(#theta_{#tau})";
  }
  for(int isDs=0; isDs<=1; isDs++){
    can.cd(isDs+1);
    for(int histo=0; histo<nHis; histo++){
      hName = "hQ2"; hName += isDs; hName += histo;
      hQ2[isDs][histo] = new TH1F(hName,"",nBins[histo],limQ2[0], limQ2[1]);
      hQ2[isDs][histo]->SetLineWidth(lineWidth[histo]);
      hQ2[isDs][histo]->SetLineColor(color[histo]);
      if(histo==0){
	if(isTheta) {
	  if(isDs==0) varName[isTheta] = "-cos(trueCTL)";
	  else        varName[isTheta] = "-trueCTL";
	}
	mc.Project(hName,varName[isTheta],cuts[isDs]);
	hQ2[isDs][histo]->Sumw2();
	style.setMarkers(hQ2[isDs][histo], 0.6, 20);
      } else {
	isCLN = 2-histo;
	for(int bin=1; bin<=nBins[histo]; bin++){
	  q2 = hQ2[isDs][histo]->GetBinCenter(bin);
	  if(isTheta==0){
	    if(isDs==0) hQ2[isDs][histo]->SetBinContent(bin,1e16*Dtaunu.Compute(q2,isCLN,ml[histo]));
	    else        hQ2[isDs][histo]->SetBinContent(bin,1e16*Dstaunu.Compute(q2,isCLN,ml[histo]));
	  } else {
	    if(isDs==0) hQ2[isDs][histo]->SetBinContent(bin,Dtaunu .ComputetL(acos(q2),isCLN,ml[histo]));
	    else        hQ2[isDs][histo]->SetBinContent(bin,Dstaunu.ComputetL(q2,isCLN,ml[histo]));
	  }
	}
	hQ2[isDs][histo]->Scale(hQ2[isDs][0]->Integral()/hQ2[isDs][histo]->Integral()*
				(double)nBins[histo]/(double)nBins[0]);
      }
      if(maxH[isDs]<hQ2[isDs][histo]->GetMaximum()) maxH[isDs]=hQ2[isDs][histo]->GetMaximum();
      meanQ2[isDs][histo] = hQ2[isDs][histo]->GetMean();
      rmsQ2[isDs][histo] = hQ2[isDs][histo]->GetRMS();
    }
  }
  for(int histo=0; histo<nHis; histo++) leg.AddEntry(hQ2[0][iLeg[histo]], legName[iLeg[histo]]);
  for(int isDs=0; isDs<2; isDs++){
    can.cd(isDs+1);
    for(int histo=0; histo<nHis; histo++){
      if(histo==0){
	hQ2[isDs][histo]->SetMaximum(maxH[isDs]*1.15); hQ2[isDs][histo]->SetMinimum(0);
	if(isTheta) hQ2[isDs][histo]->SetMaximum(maxH[isDs]*1.2); 
	hQ2[isDs][histo]->Draw();
	TString Title = "Events/("; 
	Title += RoundNumber((limQ2[1]-limQ2[0])*factor[isDs], digits[isDs], nBins[isDs]);
	Title += yTitle[isTheta];
	style.setTitles(hQ2[isDs][histo],xTitle,Title,"",rightTag[isDs]);
      } else hQ2[isDs][histo]->Draw("c same");
    }
    if(isDs==0 && isTheta==1) leg.Draw();
    if(isDs==1 && isTheta==0) leg.Draw();

    cout<<endl<<decName[isDs]<<" ->   mean q2:  ";
    for(int histo=0; histo<nHis; histo++) cout<<RoundNumber(meanQ2[isDs][histo],3)<<"\t";
    cout<<"  RMS:  ";
    for(int histo=0; histo<nHis; histo++) cout<<RoundNumber(rmsQ2[isDs][histo],3)<<"\t";

  }
  cout<<endl<<endl;
  TString epsName = "public_html/Dataset_EvtGen_Q2.eps";
  if(isTheta) epsName = "public_html/Dataset_EvtGen_thetaL.eps";
  can.SaveAs(epsName);
  for(int isDs=0; isDs<=1; isDs++)
    for(int histo=0; histo<nHis; histo++) 
      hQ2[isDs][histo]->Delete();
}

void PlotAngles(int isDgamma =0, double thetaL=2.5, double thetaV=1.5, double Chi=2){
  Styles style; style.setPadsStyle(2); style.yTitleOffset += 0.05; style.applyStyle();
  BToDtaunu Dtaunu; BToDstaunu Dstaunu;

  double legW = 0.22, legH = 0.25;
  double legX = style.PadLeftMargin+0.05, legY = style.PadBottomMargin+0.04;
  TLegend leg(legX, legY, legX+legW, legY+legH);
  leg.SetTextSize(style.LabelSize); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont); leg.SetBorderSize(0);

  TChain mc("GqaMCAnalysis/ntp1");
  mc.Add("AWG82/generator/Archive/4MEvtGen_B0Dxtaunu.root");
  mc.Add("AWG82/generator/Archive/4MEvtGen_BpDxtaunu.root");
  mc.Add("AWG82/generator/Archive/OneMEvtGen_B0Dxtaunu.root");
  mc.Add("AWG82/generator/Archive/OneMEvtGen_BpDxtaunu.root");
  mc.Add("AWG82/generator/Archive/5M1EvtGen_DsDxtaunu.root");
  mc.Add("AWG82/generator/Archive/5M2EvtGen_DsDxtaunu.root");
  mc.Add("AWG82/generator/Archive/5M2EvtGen_DsDxtaunu.root");
  mc.Add("AWG82/generator/Archive/2M1EvtGen_DsDxtaunu.root");
  mc.Add("AWG82/generator/Archive/2M2EvtGen_DsDxtaunu.root");
  mc.Add("AWG82/generator/Archive/2M3EvtGen_DsDxtaunu.root");
  mc.Add("AWG82/generator/Archive/2M4EvtGen_DsDxtaunu.root");
  mc.Add("AWG82/generator/Archive/2M5EvtGen_DsDxtaunu.root");

  double PI   = 3.14159265, dTheta = 0.1;
  int nBins[] = {50,150,150}, color[] = {1,4,2,1}, nHis=3, lineWidth[] = {2,1,2}, isCLN=0;
  int nPoints = 1000, iLeg[] = {0,1,2}, digits[] = {2, 2};
  TString cuts[] = {"(MCType==6||MCType==12)&&isDgamma==","(MCType==6||MCType==12)&&isDgamma=="};
  cuts[0] += isDgamma; cuts[1] += isDgamma;
  for(int isDs=0; isDs<=1; isDs++){
      cuts[isDs] += "&&trueCTL<"; cuts[isDs] += cos(thetaL-dTheta);
      cuts[isDs] += "&&trueCTL>"; cuts[isDs] += cos(thetaL+dTheta);
    if(isDs){
      cuts[isDs] += "&&trueCTV<"; cuts[isDs] += cos(thetaV-dTheta);
      cuts[isDs] += "&&trueCTV>"; cuts[isDs] += cos(thetaV+dTheta);
    } else {
      cuts[isDs] += "&&trueChi>"; cuts[isDs] += Chi-dTheta;
      cuts[isDs] += "&&trueChi<"; cuts[isDs] += Chi+dTheta;
    }
  }
  double limQ2[2][2] = {{-1, 1},{0,2*PI}}, ml[] = {Dtaunu.mTau,Dtaunu.mTau,Dtaunu.mTau,Dtaunu.mTau};
  double rmsQ2[2][4], meanQ2[2][4];
  TString yTitle[]= {")",")"}, hName;
  TString xTitle[] = {"cos(#theta_{V})","#chi"}, legName[] = {"MC", "HQET","ISGW2"};
  TString decName[] = {"Dtaunu ","D*taunu"}, varName[] = {"trueCTV","trueChi"};
  double theta, maxH[] = {0,0}, factor[] = {1, 1};
  TCanvas can; can.Divide(2,1);
  TH1F *hQ2[2][3];
  for(int isDs=0; isDs<=1; isDs++){
    can.cd(isDs+1);
    for(int histo=0; histo<nHis; histo++){
      hName = "hQ2"; hName += isDs; hName += histo;
      hQ2[isDs][histo] = new TH1F(hName,"",nBins[histo],limQ2[isDs][0], limQ2[isDs][1]);
      hQ2[isDs][histo]->SetLineWidth(lineWidth[histo]);
      hQ2[isDs][histo]->SetLineColor(color[histo]);
      if(histo==0){
	mc.Project(hName,varName[isDs],cuts[isDs]);
	hQ2[isDs][histo]->Sumw2();
	style.setMarkers(hQ2[isDs][histo], 0.6, 20);
      } else {
	isCLN = 2-histo;
	for(int bin=1; bin<=nBins[histo]; bin++){
	  theta = hQ2[isDs][histo]->GetBinCenter(bin);
	  double q2min = ml[histo]*ml[histo], q2max = pow(Dstaunu._mB-Dstaunu._mDs,2);
	  double dq2 = (q2max-q2min)/(double)nPoints, IntegratedTL=0, valq2=q2min, valRate;
	  for(int iq2=0; iq2<=nPoints; iq2++){
	    if(isDs==0) valRate = Dstaunu.Compute(valq2,-cos(thetaL),theta,Chi,isDgamma,true,isCLN,ml[histo]);
	    else        valRate = Dstaunu.Compute(valq2,-cos(thetaL),cos(thetaV),theta,isDgamma,true,isCLN,ml[histo]);
	    if(iq2==0 || iq2==nPoints) valRate /= 2.;
	    IntegratedTL += valRate;
	    valq2 += dq2; 
	  }
	  hQ2[isDs][histo]->SetBinContent(bin,IntegratedTL);
	}
	hQ2[isDs][histo]->Scale(hQ2[isDs][0]->Integral()/hQ2[isDs][histo]->Integral()*
				(double)nBins[histo]/(double)nBins[0]);
      }
      if(maxH[isDs]<hQ2[isDs][histo]->GetMaximum()) maxH[isDs]=hQ2[isDs][histo]->GetMaximum();
      meanQ2[isDs][histo] = hQ2[isDs][histo]->GetMean();
      rmsQ2[isDs][histo] = hQ2[isDs][histo]->GetRMS();
    }
  }
  for(int histo=0; histo<nHis; histo++) leg.AddEntry(hQ2[0][iLeg[histo]], legName[iLeg[histo]]);
  for(int isDs=0; isDs<2; isDs++){
    can.cd(isDs+1);
    for(int histo=0; histo<nHis; histo++){
      if(histo==0){
	TString sTheta = "#theta_{#tau} = "; sTheta += RoundNumber(thetaL,1); 
	if(isDs==1) {
	  sTheta += ", #theta_{V} = "; sTheta += RoundNumber(thetaV,1);
	} else {
	  sTheta += ", #chi = "; sTheta += RoundNumber(Chi,1);
	}
	hQ2[isDs][histo]->SetMaximum(maxH[isDs]*1.3); hQ2[isDs][histo]->SetMinimum(0);
	hQ2[isDs][histo]->Draw();
	TString Title = "Events/("; 
	Title += RoundNumber((limQ2[1]-limQ2[0])*factor[isDs], digits[isDs], nBins[isDs]);
	Title += yTitle[isDs];
	style.setTitles(hQ2[isDs][histo],xTitle[isDs],Title,sTheta);
      } else hQ2[isDs][histo]->Draw("c same");
    }
    if(isDs==1) leg.Draw();

    cout<<endl<<decName[isDs]<<" ->   mean theta:  ";
    for(int histo=0; histo<nHis; histo++) cout<<RoundNumber(meanQ2[isDs][histo],3)<<"\t";
    cout<<"  RMS:  ";
    for(int histo=0; histo<nHis; histo++) cout<<RoundNumber(rmsQ2[isDs][histo],3)<<"\t";

  }
  cout<<endl<<endl;
  TString epsName = "public_html/Dataset_EvtGen_Angles_"; epsName += isDgamma; epsName+=".eps";
  can.SaveAs(epsName);
  for(int isDs=0; isDs<=1; isDs++)
    for(int histo=0; histo<nHis; histo++) 
      hQ2[isDs][histo]->Delete();
}

