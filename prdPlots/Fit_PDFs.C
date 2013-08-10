#include "TString.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TFile.h"
#include "TBox.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TTree.h"
#include <fstream>
#include <iostream>

#define nPads 10
using namespace TMath;
using namespace std;
using std::cout;
using std::endl;

void Fit_PDFs(int oneCol = 1){
  TString Folder = "keys/eps/CSample/root/";
  TString PadLabel[] = {"(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"};
  TString fName[] = {"EpsKeys_1_Fit", "EpsKeys_9_Fit", "EpsKeys_13_Fit", 
		     "EpsKeys_17_Fit", "EpsKeys_51_Fit", };
  TString histoName[2][4] = {{"m24", "hm2_4", "hm2_4sigma", "hm2_4sigmasigma"},
			    {"pl2", "hpl_2", "hpl_2sigma", "hpl_2sigmasigma"}};
  TString histoOption[4] = {"e0", "c", "e4", "e4"};
  TString Units[] = {"0.25 GeV^{2})", "50 MeV)"};
  int canW = 700;
  double TextSize = 0.11, labelSize = 0.15;
  if(oneCol) {
    canW = 500; TextSize = 0.13, labelSize = 0.16;
  }
  TCanvas can("dataMC","data Vs MC",canW,140*5);
  TPad *Pads[nPads][2];
  TH1F *hMC[nPads][4], *hTemp;
  TString hName;

  int nRows = 5, nCols = 2, nFont = 132, nPDFs = 4;
  int Colors[] = {1,4,38,9};
  TLatex label; label.SetTextSize(labelSize); label.SetTextFont(132); 
  label.SetTextAlign(33);label.SetNDC(kTRUE);

  double dRows = (double)nRows;
  double bMargin = 0.06, top = 0., padH = (1-top-bMargin)/dRows, padW = 0.5, yOffset = 0.5;
  double LeftMargin = 0.19, RightMargin = 0.03, BottomMargin = 0, TopMargin = 0;
  if(oneCol) {
    LeftMargin = 0.25; bMargin = 0.07; yOffset = 0.55;
  }
  for(int row=0; row<nRows; row++){
    TString histosName = Folder; histosName += fName[row]; histosName += ".root";
    TFile fHistos(histosName);
    for(int col=0; col<nCols; col++){
      int pad = col+nCols*row;
      //cout<<"Doing pad "<<pad<<": File "<<histosName<<endl;
      for(int pdf=0; pdf<nPDFs; pdf++){
	fHistos.cd(); 
	//cout<<"Doing pdf "<<pdf<<": histo "<<histoName[col][pdf]<<endl;
	hTemp = (TH1F *)fHistos.Get(histoName[col][pdf]);
	hName = "MC_";  hName += pad; hName += pdf;
	if(hTemp) hMC[pad][pdf]  = (TH1F*)hTemp->Clone(hName);
	else {cout<<row<<", "<<col<<": PDF "<<pdf<<endl; return;}
	hMC[pad][pdf]->SetDirectory(0);
	if(pdf>1) hMC[pad][pdf]->SetFillColor(Colors[pdf]);
	if(pdf==1) {
	  hMC[pad][pdf]->SetLineWidth(1);
	  hMC[pad][pdf]->SetLineColor(Colors[pdf]);
	  int nBins = hMC[pad][pdf]->GetNbinsX();
	  for(int bin=1; bin<=nBins; bin++) {
	    hMC[pad][pdf]->SetBinError(bin, 0);
	  }
	}
      }

      double PadXY[2][2] = {{padW*col, padW*(col+1)},{padH*(dRows-1-row)+bMargin, padH*(dRows-row)+bMargin}};
      if(row==0) {
	PadXY[1][1] = 1; TopMargin = top/(top+padH);
      }
      if(row==nRows-1) {
	PadXY[1][0] = 0; BottomMargin = bMargin/(bMargin+padH);
	if(col==0) TextSize *= padH/(bMargin+padH);
	yOffset *= (bMargin+padH)/padH;
	label.SetTextSize(labelSize*padH/(bMargin+padH));
      }
      hName = "Pad0_"; hName += pad;
      //cout<<hName<<": "<<PadXY[0][0]<<", "<<PadXY[1][0]+padH*fRatio<<", "<<PadXY[0][1]<<", "<<PadXY[1][1]<<endl;
      Pads[pad][0] = new TPad(hName,"",PadXY[0][0], PadXY[1][0], PadXY[0][1], PadXY[1][1]);
      Pads[pad][0]->SetLeftMargin(LeftMargin);      Pads[pad][0]->SetRightMargin(RightMargin); 
      Pads[pad][0]->SetBottomMargin(BottomMargin);  Pads[pad][0]->SetTopMargin(TopMargin);
      can.cd();
      Pads[pad][0]->Draw(); Pads[pad][0]->cd();

      hMC[pad][0]->SetMarkerStyle(20); hMC[pad][0]->SetMarkerSize(0.28);
      hMC[pad][0]->SetTitleSize(TextSize*1.2,"xy");  // Set the 2 axes title size
      hMC[pad][0]->SetLabelSize(TextSize,"xy");      // Set the 2 axes label size
      hMC[pad][0]->SetTitleFont(nFont,"xy");         // Set the all 2 axes title font
      hMC[pad][0]->SetLabelFont(nFont,"xy");         // Set the all 2 axes label font
      hMC[pad][0]->SetNdivisions(504, "xy");         // 5 primary ticks and 4 secondary ticks
      hMC[pad][0]->SetTitleOffset(yOffset,"y");         // Set y offset
      hName = "m"; if(col) hName = "p"; if(row<nRows-1) hName = "";
      hMC[pad][0]->GetXaxis()->CenterTitle(true);
      hMC[pad][0]->SetXTitle(hName);
      hMC[pad][0]->SetYTitle("");

      double maxH = hMC[pad][0]->GetMaximum();
      //if(maxH<hMC[pad][nPDFs-1]->GetMaximum()) maxH = hMC[pad][nPDFs-1]->GetMaximum();
      if(col==0) hMC[pad][0]->SetMaximum(maxH*1.01);
      else if(row==1 || row==3)  hMC[pad][0]->SetMaximum(maxH*1.2);
      else hMC[pad][0]->SetMaximum(maxH*1.14);
      hMC[pad][0]->Draw("");
      for(int pdf=nPDFs-1; pdf>=0; pdf--) {
	histoOption[pdf] += " same";
	hMC[pad][pdf]->Draw(histoOption[pdf]);
      }
      //hMC[pad][0]->Draw("axis same");
      //hMC[pad][0]->Draw("e0 same");
      label.DrawLatex(0.93, 0.89, PadLabel[row]);

    }
  }
  can.cd(0);
  label.SetTextAngle(90); label.SetTextAlign(23); label.SetTextSize(labelSize*0.025/0.15);
  if(oneCol) label.SetTextSize(labelSize*0.035/0.15);
  hName = "Events/("; hName += Units[0];
  label.DrawLatex(0.008, 0.53, hName);
  hName = "Events/("; hName += Units[1];
  label.DrawLatex(0.515, 0.53, hName);


  TBox box; box.SetFillStyle(1001); box.SetLineColor(10); box.SetFillColor(10);
  label.SetTextAngle(0); label.SetTextAlign(13); label.SetTextSize(labelSize*0.035/0.15/1.19); 
  for(int row=0; row<nRows-1; row++){
    for(int col=0; col<nCols; col++){
      box.DrawBox(padW*(LeftMargin+col)-0.03,  bMargin+padH*(dRows-row-1), 
		  padW*(LeftMargin+col)-0.002, bMargin+padH*(dRows-row-1)+0.02);
      label.DrawLatex(padW*(LeftMargin+col)-0.017, bMargin+padH*(dRows-row-1)+0.009,"0");
    }
  }

  TString pName = "public_html/Fit_PDFs.eps"; 
  can.SaveAs(pName);

  for(int row=0; row<nRows; row++){
    for(int col=0; col<nCols; col++){
      //int pad = col+3*row;
      //for(int isRat=0; isRat<2; isRat++) if(Pads[pad][isRat]) Pads[pad][isRat]->Delete();
    }
  }
  
}

