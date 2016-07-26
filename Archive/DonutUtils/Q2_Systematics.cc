#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "DonutUtils/RooNDKeysPdfDonut.hh"
#include "DonutUtils/KeysUtils.cc"
#include "RooFitCore/RooArgList.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooRealVar.hh"
#include <fstream>
#include <iostream>

#define nCuts 3
#define nFont 132
using namespace TMath;
using namespace std;
using std::cout;
using std::endl;

void formatHisto(TH1F *histo, int color, double TextSize);

int main(int argc, char *argv[]){
  double smoothing = 1.2;
  if(argc>1) {TString smoo_s = argv[1]; smoothing = smoo_s.Atof();}
  int nBins = 16;
  if(argc>2) {TString smoo_s = argv[2]; nBins = smoo_s.Atoi();}

  TString Folder = "AWG82/ntuples/small/FitRAll", treeName[] = {"Newx100", "Newx100", "Dpipix100", "Newx100"};
  TChain *MC[nCuts]; TTree *cPDF[nCuts][2];
  
  gStyle->SetOptStat(0); gROOT->SetStyle("Plain"); gStyle->SetOptStat(0); 
  int nRows = 2, nCols = 2, nBinsPDF = 30, Colors[] = {1,4,2,8};
  double dRows = (double)nRows, dCols = (double) nCols;
  double bMargin = 0.12, padH = (1-bMargin)/dRows, padW = 1/dCols, LeftMargin = 0.2;
  double limQ2[] = {4, 12}, TextSize = 0.08;
  TCanvas can("dataMC","data Vs MC",233*nCols,150*nRows); 
  TPad *Pads[nCuts][2];
  TString hName;
  TH1F *hMC[nCuts][2], *hPDF[nCuts][2];
  TString mcCuts[nCuts][2] = {{"weight*(MCType==14&&MCTaumode<0&&(candType==1||candType==3)&&candM2>1.5)",
			       "weight*(MCType==14&&MCTaumode<0&&(candType==2||candType==4)&&candM2>1.5)"},
			      {"weight*(MCType==14&&MCTaumode>=0&&(candType==1||candType==3)&&candM2>1.5)",
			       "weight*(MCType==14&&MCTaumode>=0&&(candType==2||candType==4)&&candM2>1.5)"},
 			      {"weight*(MCType==16&&(candType==1||candType==3)&&candM2>1.5)",
 			       "weight*(MCType==16&&(candType==2||candType==4)&&candM2>1.5)"}};
//  			      {"(MCType==13&&(candType==1||candType==3)&&candM2>1.5)",
// 				  "(MCType==13&&(candType==2||candType==4)&&candM2>1.5)"}};



  RooRealVar RRVq2("candQ2","candQ2",4,12.5);
  RooRealVar RRVweight("weight","weight",0.,100.);
  TString textName = "public_html/Uncertainties_BkgDss.txt";
  fstream textFile; textFile.open(textName,fstream::out);

  for(int cut=0; cut<nCuts; cut++){
    hName = Folder; hName += treeName[cut]; hName += "_RunAll.root";
    MC[cut] = new TChain("ntp1");
    MC[cut]->Add(hName);
    for(int row=0; row<nRows; row++){
      hName = "hMC"; hName += cut; hName += row; 
      hMC[cut][row] = new TH1F(hName,"", nBins, limQ2[0], limQ2[1]);
      formatHisto(hMC[cut][row], Colors[cut], TextSize);
      MC[cut]->Project(hName,"candQ2",mcCuts[cut][row]);
      if(cut>0) hMC[cut][row]->Scale(hMC[0][row]->Integral()/hMC[cut][row]->Integral());

      hName = "hPDF"; hName += cut; hName += row; 
      hPDF[cut][row] = new TH1F(hName,"", nBinsPDF*nBins, limQ2[0], limQ2[1]);
      hPDF[cut][row]->SetDirectory(0);
      cPDF[cut][row] = MC[cut]->CopyTree(mcCuts[cut][row]);
      RooDataSet cutData("cutData","",cPDF[cut][row],RooArgSet(RRVq2,RRVweight),"","weight");
      //cout<<"Doing cut "<<cut<<" with bins "<<hPDF[cut][row]->GetBinX()<<endl;
      RooNDKeysPdfDonut DssPdf("DssPdf","",RooArgList(RRVq2),cutData,"a",smoothing,3);
      DssPdf.fillHistogram(hPDF[cut][row],RooArgList(RRVq2));
      hPDF[cut][row]->Scale(hMC[0][row]->Integral()/hPDF[cut][row]->Integral()*(double)(nBins*nBinsPDF)/(double)nBins);
      hPDF[cut][row]->SetLineColor(Colors[cut]);
    }
  }

  for(int col=0; col<nCols; col++){
    for(int row=0; row<nRows; row++){
      can.cd(0);
      int pad = nRows*col+row;
      double RightMargin = 0.04, TopMargin=0, BottomMargin=bMargin/(bMargin+padH)*row;
      double PadXY[2][2] = {{padW*col, padW*(col+1)},{(padH+bMargin)*(dRows-1-row), bMargin+padH*(2-row)}};
      if(row==0) TextSize *= (bMargin+padH)/padH;
      hName = "Pad0_"; hName += pad;
      //cout<<hName<<": "<<PadXY[0][0]<<", "<<PadXY[1][0]<<", "<<PadXY[0][1]<<", "<<PadXY[1][1]<<endl;
      Pads[pad][0] = new TPad(hName,"",PadXY[0][0], PadXY[1][0], PadXY[0][1], PadXY[1][1]);
      Pads[pad][0]->SetLeftMargin(LeftMargin);     Pads[pad][0]->SetRightMargin(RightMargin); 
      Pads[pad][0]->SetBottomMargin(BottomMargin); Pads[pad][0]->SetTopMargin(TopMargin); 
      Pads[pad][0]->Draw(); Pads[pad][0]->cd();

      hMC[0][row]->SetMinimum(0);
      hMC[0][row]->SetMaximum(hMC[0][row]->GetMaximum()*1.3);
      hMC[0][row]->Draw("hist");
      hMC[col+1][row]->Draw("e0 same");
      hPDF[0][row]->Draw("c same");
      hPDF[col+1][row]->Draw("c same");
   }
  }
  for(int row=0; row<nRows; row++){
    double q2 = limQ2[0], dq2 = hMC[0][0]->GetXaxis()->GetBinWidth(1);
    for(int bin=1; bin<=nBins; bin++){
      textFile<<RoundNumber(q2,1)<<"-"<<RoundNumber(q2+dq2,1)<<"  \t"<<RoundNumber(hMC[0][row]->GetBinContent(bin),2);
      textFile<<"\t"<<RoundNumber(hMC[0][row]->GetBinError(bin),2)<<" \t\t";
      for(int cut=1; cut<nCuts; cut++){
	double val=0;
	for(int binP=1; binP<=nBinsPDF; binP++) val += hPDF[cut][row]->GetBinContent(nBinsPDF*(bin-1)+binP);
	textFile<<RoundNumber(val,2,(double)nBinsPDF)<<"\t"<<RoundNumber(hMC[cut][row]->GetBinError(bin),2)<<" \t\t";
      }
      textFile<<endl;
      q2 += dq2;
    }
    textFile<<endl;
  }

  textFile.close();
  cout<<"Written "<<textName<<endl;
  TString pName = "public_html/Q2_Systematics.eps"; 
  can.SaveAs(pName);
  for(int row=0; row<nRows; row++){
    for(int cut=0; cut<nCuts; cut++){
      if(hMC[cut][row]) hMC[cut][row]->Delete();
      //if(hPDF[cut][row]) hMC[cut][row]->Delete();
    }
  }
}

void formatHisto(TH1F *histo, int color, double TextSize){
  histo->Sumw2();
  histo->SetMarkerStyle(20); histo->SetMarkerSize(0.6); 
  histo->SetLineColor(color); histo->SetMarkerColor(color);
  
  histo->SetTitleSize(TextSize,"xy");      // Set the 2 axes title size
  histo->SetLabelSize(TextSize,"xy");      // Set the 2 axes label size
  histo->SetTitleFont(nFont,"xy");         // Set the all 2 axes title font
  histo->SetLabelFont(nFont,"xy");         // Set the all 2 axes label font
  histo->SetNdivisions(904, "xy");         // 5 primary ticks and 4 secondary ticks
  histo->SetLabelOffset(0.01,"Y");
  histo->GetXaxis()->CenterTitle(true);

}
