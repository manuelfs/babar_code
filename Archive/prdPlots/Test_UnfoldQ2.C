#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TFile.h"
#include "THStack.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLegend.h"
#include "babar_code/FF/BToDtaunu.cc"
#include "babar_code/FF/BToDstaunu.cc"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

#define nPads 6
#define nType 2
#define nFont 132
using namespace TMath;
using namespace std;
using std::cout;
using std::endl;

void formatHisto(TH1D *histo, int color, double TextSize);

void Test_UnfoldQ2(int Type = 0){
  TString TypeName[] = {"Efficiency", "Unfold"};
  if(Type==-1) {
    for(int type=0; type<nType; type++)
      cout<<type<<" "<<TypeName[type]<<", "; 
    cout<<endl; 
    return;
  }
  BToDtaunu Dtaunu; BToDstaunu Dstaunu;
  int nBins = 18, nValues = 50, nBinsRes = 80;
  double limQ2[] = {3.5, 12.5}, limRes = 2, ml[] = {Dtaunu.mTau, Dtaunu.mTau, Dtaunu.mMu};
  if(Type==0){limQ2[0] = 4; nBins = 17;}

  TString mcCuts[3][6] = {{"weight*((MCType==5&&candType==1||MCType==11&&candType==3)&&candM2>1.5)",
			   "weight*((MCType==6&&candType==2||MCType==12&&candType==4)&&candM2>1.5)" ,
			   "weight*((MCType==5&&candType==1||MCType==11&&candType==3))",
			   "weight*((MCType==6&&candType==2||MCType==12&&candType==4))",
			   "weight*(((MCType==1||MCType==3)&&candType==1||(MCType==7||MCType==9)&&candType==3))",
			   "weight*(((MCType==2||MCType==4)&&candType==2||(MCType==8||MCType==10)&&candType==4))"},
			  {"weight*((MCType==11&&candType==3)&&candM2>1.5)",
			   "weight*((MCType==12&&candType==4)&&candM2>1.5)" ,
			   "weight*((MCType==11&&candType==3))",
			   "weight*((MCType==12&&candType==4))",
			   "weight*(((MCType==7||MCType==9)&&candType==3))",
			   "weight*(((MCType==8||MCType==10)&&candType==4))"},
			  {"weight*((MCType==5&&candType==1)&&candM2>1.5)",
			   "weight*((MCType==6&&candType==2)&&candM2>1.5)" ,
			   "weight*((MCType==5&&candType==1))",
			   "weight*((MCType==6&&candType==2))",
			   "weight*(((MCType==1||MCType==3)&&candType==1))",
			   "weight*(((MCType==2||MCType==4)&&candType==2))"}};

  TString PadLabel[] = {"a", "b", "x", "d", "e", "f", "(g)", "(h)", "(i)"};
  TString Units[] = {" GeV^{2})", "0.35 GeV^{2})"}, hName;
  TString xTitle[] = {"Q", "R"};
  TBox box; box.SetLineColor(10);box.SetFillColor(10);

  gStyle->SetOptStat(0);
  int nRows = 2, nCols = 3;
  int Colors[] = {8,4,2}, ColorP[] = {4,1,2,8};
  TLatex label; label.SetTextFont(132); label.SetNDC(kTRUE);
  TCanvas can("dataMC","data Vs MC",700,150*nRows); 
  TPad *Pads[nPads][2];
  TH1D *hbini[nPads], *hxini[nPads][3], *hQ2[nPads][3], *hEff[nPads][3], *hRes[nPads], *hRes2[nPads];
  TH1D *hP[2][4], hTemp("hTemp","",100,0,3.4);
  TChain MC("ntp1");
  MC.Add("AWG82/ntuples/small/FitRAllNewx100_RunAll.root");

  TLegend *leg1[3];
  double scaleEff = 1;
  double dRows = (double)nRows, dCols = (double) nCols;
  double bMargin = 0.12, padH = (1-bMargin)/dRows, padW = 1/dCols, LeftMargin = 0.2;
  for(int col=0; col<nCols; col++){
    for(int row=0; row<nRows; row++){
      can.cd(0);
      int pad = nRows*col+row;
      double TextSize = 0.08;
      double RightMargin = 0.04, TopMargin=0, BottomMargin=bMargin/(bMargin+padH)*row;
      double PadXY[2][2] = {{padW*col, padW*(col+1)},{(padH+bMargin)*(dRows-1-row), bMargin+padH*(2-row)}};
      if(row==0) TextSize *= (bMargin+padH)/padH;
      hName = "Pad0_"; hName += pad;
      //cout<<hName<<": "<<PadXY[0][0]<<", "<<PadXY[1][0]<<", "<<PadXY[0][1]<<", "<<PadXY[1][1]<<endl;
      Pads[pad][0] = new TPad(hName,"",PadXY[0][0], PadXY[1][0], PadXY[0][1], PadXY[1][1]);
      Pads[pad][0]->SetLeftMargin(LeftMargin);     Pads[pad][0]->SetRightMargin(RightMargin); 
      Pads[pad][0]->SetBottomMargin(BottomMargin); Pads[pad][0]->SetTopMargin(TopMargin); 
      Pads[pad][0]->Draw(); Pads[pad][0]->cd();

      if(Type==0){
	for(int isNorm=2; isNorm>=0; isNorm--){
	  hName = "xini"; hName += pad; hName += isNorm;
	  hxini[pad][isNorm] = new TH1D(hName,"",nBins, limQ2[0], limQ2[1]);
	  formatHisto(hxini[pad][isNorm], Colors[isNorm], TextSize);
	  MC.Project(hName,"candQ2",mcCuts[col][row+2*isNorm]);
	  if(row==1) hxini[pad][isNorm]->SetXTitle(xTitle[0]);

	  //cout<<pad<<" Calculating q2"<<isNorm<<endl;
	  hName = "hQ2"; hName += pad; hName += isNorm;
	  hQ2[pad][isNorm] = new TH1D(hName,"",nBins, limQ2[0], limQ2[1]);
	  hQ2[pad][isNorm]->SetLineColor(Colors[isNorm]);
	  for(int bin=1; bin<=nBins; bin++){
	    double q2 = hQ2[pad][isNorm]->GetBinLowEdge(bin), val=0;
	    double dq2 = (hQ2[pad][isNorm]->GetBinLowEdge(bin)-q2)/(double)nValues;
	    q2 += dq2/2.;
	    for(int ival=0; ival<nValues; ival++) {
	      if(row==0) val += 1e16*Dtaunu.Compute(q2,1, ml[isNorm]);
	      else       val += 1e16*Dstaunu.Compute(q2,1, ml[isNorm]);
	      q2 += dq2;
	    }
	    hQ2[pad][isNorm]->SetBinContent(bin,val); hQ2[pad][isNorm]->SetBinError(bin,0);
	  }
	  hEff[pad][isNorm] = (TH1D*)hxini[pad][isNorm]->Clone(); 
	  hEff[pad][isNorm]->Divide(hQ2[pad][isNorm]);
	  if(isNorm==2) scaleEff = 1/hEff[pad][isNorm]->GetMaximum();
	  if(isNorm==1) scaleEff = hEff[pad][2]->Integral(nBins/3,nBins/3*2)/
	    hEff[pad][1]->Integral(nBins/3,nBins/3*2);
	  hEff[pad][isNorm]->Scale(scaleEff);
	  hEff[pad][isNorm]->SetMinimum(0);
	  if(row==1) hEff[pad][isNorm]->SetMaximum(1.72);
	}
	//cout<<pad<<" Drawing"<<endl;
	hEff[pad][1]->Draw();
	if(col==0) hEff[pad][0]->Draw("same"); 
	hEff[pad][2]->Draw("same"); hEff[pad][1]->Draw("same");
	//hQ2[pad][2]->Draw(); hQ2[pad][1]->Draw("same"); 

	// Labels
	label.SetTextAngle(0); label.SetTextSize(TextSize/0.93); 
	TString chanLabel = PadLabel[pad]; 
	double Xchan = 0.92, Ychan = 0.89;
	if(Type>2)  Xchan = LeftMargin+0.06;
	label.SetTextAlign(22); label.DrawLatex(Xchan, Ychan, chanLabel);
      }
      if(Type==1){
	int isNorm=1;
	switch(col){
	case 0:
	  hName = "xini"; hName += pad; hName += isNorm;
	  hxini[pad][isNorm] = new TH1D(hName,"",nBins, limQ2[0], limQ2[1]);
	  formatHisto(hxini[pad][isNorm], Colors[isNorm], TextSize);
	  MC.Project(hName,"trueQ2",mcCuts[col][row+2*isNorm]);
	  if(row==1) hxini[pad][isNorm]->SetXTitle(xTitle[0]);
	  hName = "bini"; hName += pad; hName += isNorm;
	  hbini[pad] = new TH1D(hName,"",nBins, limQ2[0], limQ2[1]);
	  formatHisto(hbini[pad], 1, TextSize);
	  MC.Project(hName,"candQ2",mcCuts[col][row+2*isNorm]);
	  hxini[pad][isNorm]->SetMarkerSize(0.); 
	  hxini[pad][isNorm]->SetMaximum(hxini[pad][isNorm]->GetMaximum()*1.2);
	  hxini[pad][isNorm]->Draw("hist");
	  hbini[pad]->Draw("same");
	  if(row==0){
	    double legXY[2][2] = {{0.78, 0.98}, {0.67, 0.95}};
	    leg1[col] = new TLegend(legXY[0][0], legXY[1][0], legXY[0][1], legXY[1][1]);
	    leg1[col]->SetTextSize(TextSize); leg1[col]->SetFillStyle(0); 
	    leg1[col]->SetTextFont(132); leg1[col]->SetBorderSize(0);
	    leg1[col]->AddEntry(hxini[pad][isNorm],"t");
	    leg1[col]->AddEntry(hbini[pad],"d");
	    leg1[col]->Draw();
	  }
	  break;
	case 1:
	  isNorm=1;
	  hName = "Res"; hName += row;
	  hRes[row] = new TH1D(hName,"",nBinsRes, -limRes, limRes);
	  formatHisto(hRes[row], Colors[isNorm], TextSize);
	  MC.Project(hName,"candQ2-trueQ2", mcCuts[0][row+2*isNorm]);
	  isNorm=2;
	  hName = "Res2"; hName += row;
	  hRes2[row] = new TH1D(hName,"",nBinsRes, -limRes, limRes);
	  formatHisto(hRes2[row], Colors[isNorm], TextSize);
	  hRes[row]->SetMarkerSize(0.3); hRes2[row]->SetMarkerSize(0.); 
	  MC.Project(hName,"candQ2-trueQ2", mcCuts[0][row+2*isNorm]);
	  hRes2[row]->Scale(hRes[row]->Integral()/hRes2[row]->Integral());
	  hRes2[row]->SetMaximum(hRes2[row]->GetMaximum()*1.15);
	  hRes2[row]->Draw("hist"); hRes[row]->Draw("same");
	  if(row==1) hRes2[row]->SetXTitle(xTitle[1]);
	  hName = "#splitline{#mu_{sig} = "; 
	  hName += RoundNumber(hRes[row]->GetMean(),2); hName += "}{#mu_{nor} = ";
	  hName += RoundNumber(hRes2[row]->GetMean(),2); hName += "}";
	  label.SetTextAlign(33); label.SetTextAngle(0); label.SetTextSize(TextSize);
	  label.DrawLatex(1-RightMargin-0.05, 1-0.05, hName);
	  hName = "#splitline{#sigma_{sig} = ";
	  hName += RoundNumber(hRes[row]->GetRMS(),2); hName += "}{#sigma_{nor} = ";
	  hName += RoundNumber(hRes2[row]->GetRMS(),2); hName += "}";
	  label.DrawLatex(1-RightMargin-0.05, 1-0.37-0.06*(row==0), hName);
	  if(row==0){
	    double legXY[2][2] = {{LeftMargin+0.05, LeftMargin+0.25}, {0.67, 0.95}};
	    leg1[col] = new TLegend(legXY[0][0], legXY[1][0], legXY[0][1], legXY[1][1]);
	    leg1[col]->SetTextSize(TextSize); leg1[col]->SetFillStyle(0); 
	    leg1[col]->SetTextFont(132); leg1[col]->SetBorderSize(0);
	    leg1[col]->AddEntry(hRes2[row],"Norm.");
	    leg1[col]->AddEntry(hRes[row],"Signal");
	    leg1[col]->Draw();
	  }
	  break;
	case 2:
	  double maxH = 0;
	  for(isNorm=0; isNorm<3; isNorm++){
	    hName = "hP"; hName += pad;  hName += isNorm;
	    hP[pad][isNorm] = new TH1D(hName,"",nBins, limQ2[0], limQ2[1]);
	    formatHisto(hP[pad][isNorm], ColorP[isNorm], TextSize);
	    double q2 = limQ2[0], dq2 = (limQ2[1]-limQ2[0])/(double)nBins;
	    for(int bin=1; bin<=nBins; bin++){
	      hName = "weight*(trueQ2>="; hName += q2;
	      hName += "&&trueQ2<"; hName += q2+dq2; hName += "&&";
	      TString Cuts = mcCuts[0][2+row+2*(isNorm>1)];
	      Cuts.ReplaceAll("weight*(", hName);
	      TString var = "candPLep"; if(isNorm%2==1) var = "candPstarD";
	      double entries = MC.Project("hTemp",var, Cuts);
	      //cout<<hTemp.GetMean()<<", "<<Cuts<<endl;
	      hP[pad][isNorm]->SetBinContent(bin, hTemp.GetMean());
	      if(entries>1) hP[pad][isNorm]->SetBinError(bin, hTemp.GetRMS()/sqrt(entries));
	      else hP[pad][isNorm]->SetBinError(bin, hTemp.GetMean());
	      
	      q2 += dq2;
	    }
	    if(hP[pad][isNorm]->GetMaximum()>maxH) maxH = hP[pad][isNorm]->GetMaximum();
	  }
	  for(isNorm=0; isNorm<3; isNorm++){
	    if(isNorm==0) {
	      hP[pad][isNorm]->SetMaximum(maxH*1.17);
	      hP[pad][isNorm]->Draw();
	      if(row==1) hP[pad][isNorm]->SetXTitle(xTitle[0]);
	    } else{
	      //hP[pad][isNorm]->Scale(hP[pad][0]->Integral()/hP[pad][isNorm]->Integral());
	      hP[pad][isNorm]->Draw("same");
	    }
	  }
	  if(row==1){
	    double legXY[2][2] = {{LeftMargin+0.05, LeftMargin+0.25}, {0.25, 0.5}};
	    leg1[col] = new TLegend(legXY[0][0], legXY[1][0], legXY[0][1], legXY[1][1]);
	    leg1[col]->SetTextSize(TextSize); leg1[col]->SetFillStyle(0); 
	    leg1[col]->SetTextFont(132); leg1[col]->SetBorderSize(0);
	    leg1[col]->AddEntry(hP[pad][0],"ls");
	    leg1[col]->AddEntry(hP[pad][2],"ln");
	    leg1[col]->AddEntry(hP[pad][1],"D");
	    leg1[col]->Draw();
	  }

	  break;
	}
      }
    }
    can.cd(0);
    hName = "Normalized efficiency/("; 
    if(Type==1) hName = "Events/("; 
    if(Type==1&&col==1) hName += RoundNumber(2*limRes,2,nBinsRes);
    else hName += RoundNumber((limQ2[1]-limQ2[0]),2,nBins);
    hName += Units[0];
    if(Type==1&&col==2) hName = "Momentum (GeV)";
    label.SetTextAlign(23); label.SetTextAngle(90); label.SetTextSize(0.055);
    label.DrawLatex(0.01+col/dCols, 0.56, hName);
    if(Type==0 || Type==1) {
      box.DrawBox(padW*(LeftMargin+col)-0.02, bMargin+padH, padW*(LeftMargin+col)-0.002, bMargin+padH+0.02);
      label.SetTextAlign(32); label.SetTextAngle(0); label.SetTextSize(0.044);
      label.DrawLatex(padW*(LeftMargin+col)-0.0025, bMargin+padH+0.003, "0");
    }

  }



  TString pName = "public_html/Test_Q2_"; pName += TypeName[Type]; pName += ".eps"; 
  can.SaveAs(pName);
  for(int row=0; row<nRows; row++){
    if(Type==1){
      if(hRes[row]) hRes[row]->Delete();
      if(hRes2[row]) hRes2[row]->Delete();
    }
    for(int col=0; col<nCols; col++){
      int pad = nRows*col+row;
      //if(Type==1 && hbini[pad]) hbini[pad]->Delete();
      for(int isNorm=0; isNorm<3; isNorm++){
	if(Type==0){
	  if(hEff[pad][isNorm]) hEff[pad][isNorm]->Delete();
	  if(hQ2[pad][isNorm]) hQ2[pad][isNorm]->Delete();
	  if(hxini[pad][isNorm]) hxini[pad][isNorm]->Delete();
	}
      }
    }
  }
}


void formatHisto(TH1D *histo, int color, double TextSize){
  histo->Sumw2();
  histo->SetMarkerStyle(20); histo->SetMarkerSize(0.5); 
  histo->SetLineColor(color); histo->SetMarkerColor(color);
  
  histo->SetTitleSize(TextSize,"xy");      // Set the 2 axes title size
  histo->SetLabelSize(TextSize,"xy");      // Set the 2 axes label size
  histo->SetTitleFont(nFont,"xy");         // Set the all 2 axes title font
  histo->SetLabelFont(nFont,"xy");         // Set the all 2 axes label font
  histo->SetNdivisions(904, "xy");         // 5 primary ticks and 4 secondary ticks
  histo->SetLabelOffset(0.01,"Y");
  histo->GetXaxis()->CenterTitle(true);

}
