#include "TString.h"
#include "TPad.h"
#include "TF1.h"
#include "TCut.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TFile.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TBox.h"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

#define nPads 9
#define nType 3
using namespace TMath;
using namespace std;
using std::cout;
using std::endl;

void Control_Samples(int Type=2){
  TString TypeName[] = {"Control_Samples", "Test_Q2_Sideband2", "Higgs_Q2Bkg"};
  if(Type==-1) {
    for(int type=0; type<nType; type++)
      cout<<type<<" "<<TypeName[type]<<", "; 
    cout<<endl; 
    return;
  }
  TString textName = "public_html/Uncertainties_Bkg.txt";
  fstream textFile; if(Type==2) textFile.open(textName,fstream::out);
  TString Folder = "keys/eps/CSample/root/";
  TString PadLabel[] = {"(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"};
  TString fName[nType][6] = {{"Pl_OffNoW", "Pl_HiEexNoW", "Pl_HiEexW", 
			      "Eex_mESSB", "mES_HiEex",   "mES_LoEex"},
			     //	{"Q2_SB_050", "Q2_SB_mES", "Q2_SB_Eex", 
			     // "Q2_SB_050", "Q2_SB_mES", "Q2_SB_Eex"}};
			     {"Q2_SB2_Eex", "Q2_SB2_mES", "Q2_SB2_100", 
			      "Q2_SB2_Eex", "Q2_SB2_mES", "Q2_SB2_100"},
			     {"Q2",         "Q2_SB2_Eex", "Q2_SB2_100", 
			      "Q2",         "Q2_SB2_Eex", "Q2_SB2_100"}};
  TString Units[] = {"40 MeV)", "40 MeV)", "40 MeV)", "100 MeV)", "6 MeV)", "6 MeV)"};
  int nRows=2, nCols = 3, nFont = 132, nPDFs = 7, colDsl=kBlue-8, colDl=kCyan+1, colBkg=kYellow-7, nBinsError=500;
  int Colors[] = {kYellow+1,kYellow+1,colBkg,colBkg,28,colDl,colDsl,kRed-7,kGreen-7,1};
  TCanvas can("dataMC","data Vs MC",700,180*nRows);
  TPad *Pads[nPads][2];
  TH1F *hData[nPads], *hYields[nPads], *hRatio[nPads], *hMC[nPads][10], *hTemp, *hCopy[nPads][10];
  TH1F *hFit[4]; TF1 *fFit[2];
  TString hName;

  gStyle->SetHatchesLineWidth(1);
  int fillDl = 3654, fillDsl = 3645;
  double TextSize = 0.11, TextSizeR = 0.165;
  double LeftMargin = 0.18, RightMargin = 0.03;
  TLatex label; label.SetTextSize(0.15); label.SetTextFont(132); label.SetTextAlign(33);label.SetNDC(kTRUE);

  double dRows = (double)nRows, dCols = (double) nCols, fRatio = 0.4;
  double padH = 1/dRows, padW = 1/dCols, fSameArea=1;
  for(int row=0; row<nRows; row++){
    for(int col=0; col<nCols; col++){
      int pad = col+3*row;
      int candi = 2+row;
      if(pad==0 && Type==0) nPDFs = 2;
      if(col==0 && Type==2) nPDFs = 9;
      else nPDFs = 7;
      //cout<<"Doing pad "<<pad<<": File "<<histosName<<endl;
      TString histosName = Folder; histosName += fName[Type][pad]; histosName += ".root";
      TFile fHistos(histosName);
      fHistos.cd(); can.cd();
      hName = "Data0"; 
      if(fName[Type][pad].Contains("Res")) hName = "Data3"; 
      if(Type==1 || Type==2&&col>0){hName = "Data"; hName += row;}
      else if(Type==2){hName = "Data"; hName += candi;}
      hTemp = (TH1F *)fHistos.Get(hName);
      hName = "Data_"; hName += pad;
      if(hTemp) hData[pad]  = (TH1F*)hTemp->Clone(hName);
      else continue;
      hData[pad]->SetDirectory(0);
      if(Type==1 || Type==2&&col>0){
	hName = "Data"; hName += row+2; 
	hTemp = (TH1F *)fHistos.Get(hName);
	hData[pad]->Add(hTemp);
      }
      hName = "hYields";  hName += pad;
      hTemp = (TH1F *)fHistos.Get("hYields");
      if(hTemp) hYields[pad]  = (TH1F*)hTemp->Clone(hName);
      else continue;
      if(Type==1 || Type==2&&col>0){
	fSameArea=0;
	for(int pdf=0; pdf<=nPDFs; pdf++)
	  fSameArea += hYields[pad]->GetBinContent(pdf+10*row)+hYields[pad]->GetBinContent(pdf+10*row+20);
	fSameArea = (hYields[pad]->GetBinContent(9+10*row)+hYields[pad]->GetBinContent(9+10*row+20))/fSameArea;
      }
      for(int pdf=0; pdf<nPDFs; pdf++){
	int ipdf = pdf;
	if(pdf==7) ipdf = 8;
	if(pdf==8) ipdf = 7;
	hName = "MC1"; hName += ipdf;
	if(Type==1 || Type==2&&col>0){hName = "MC"; hName += row+1; hName += ipdf;}
	else if(Type==2){hName = "MC"; hName += candi+1; hName += ipdf;}
	int shift = 0;
	if(fName[Type][pad].Contains("Res")){
	  hName.ReplaceAll("MC1","MC4");
	  if(fName[Type][pad].Contains("Conv") && (ipdf==4 || ipdf==6)) hName.ReplaceAll("MC4","Sample3");
	  shift = 30;
	  //cout<<hYields[pad]->GetBinContent(ipdf+shift)<<", "<<endl;
	}
	hTemp = (TH1F *)fHistos.Get(hName);
	hName = "MC_";  hName += pad; hName += ipdf;
	if(hTemp) hMC[pad][pdf]  = (TH1F*)hTemp->Clone(hName);
	else {cout<<row<<", "<<col<<": PDF "<<ipdf<<endl; return;}
	hMC[pad][pdf]->SetDirectory(0);
	double yield = hYields[pad]->GetBinContent(ipdf+shift);
	if(Type==1 || Type==2&&col>0){
	  hName = "MC"; hName += row+3; hName += ipdf;
	  hTemp = (TH1F *)fHistos.Get(hName);
	  hMC[pad][pdf]->Add(hTemp);
	  yield  = hYields[pad]->GetBinContent(ipdf+shift+10*row);
	  yield += hYields[pad]->GetBinContent(ipdf+shift+10*row+20);
	  yield *= fSameArea;
	} else if(Type==2) {
	  yield = hYields[pad]->GetBinContent(ipdf+10*row);
	  yield += hYields[pad]->GetBinContent(ipdf+10*row+20);
	}
	hMC[pad][pdf]->Scale(yield/hMC[pad][pdf]->Integral());
	hMC[pad][pdf]->SetFillColor(Colors[ipdf]);
	hMC[pad][pdf]->SetLineColor(Colors[ipdf]); hMC[pad][pdf]->SetLineWidth(0);
	int nBins = hMC[pad][pdf]->GetNbinsX();
	for(int bin=1; bin<=nBins; bin++) {
	  hMC[pad][pdf]->SetBinError(bin, 0);
	}
	if(pdf>0) hMC[pad][pdf]->Add(hMC[pad][pdf-1]);
      }
      double PadXY[2][2] = {{padW*col, padW*(col+1)},{padH*(dRows-1-row), padH*(dRows-row)}};
      hName = "Pad0_"; hName += pad;
      //cout<<hName<<": "<<PadXY[0][0]<<", "<<PadXY[1][0]+padH*fRatio<<", "<<PadXY[0][1]<<", "<<PadXY[1][1]<<endl;
      Pads[pad][0] = new TPad(hName,"",PadXY[0][0], PadXY[1][0]+padH*fRatio, PadXY[0][1], PadXY[1][1]);
      Pads[pad][0]->SetLeftMargin(LeftMargin); Pads[pad][0]->SetRightMargin(RightMargin); 
      Pads[pad][0]->SetBottomMargin(0);        Pads[pad][0]->SetTopMargin(0.05); 
      Pads[pad][0]->Draw(); Pads[pad][0]->cd();

      hData[pad]->SetMarkerStyle(20); hData[pad]->SetMarkerSize(0.28);
      hData[pad]->SetTitleSize(TextSize*1.1,"xy");  // Set the 2 axes title size
      hData[pad]->SetLabelSize(TextSize,"xy");      // Set the 2 axes label size
      hData[pad]->SetTitleFont(nFont,"xy");         // Set the all 2 axes title font
      hData[pad]->SetLabelFont(nFont,"xy");         // Set the all 2 axes label font
      hData[pad]->SetNdivisions(504, "xy");         // 5 primary ticks and 4 secondary ticks
      hData[pad]->SetTitleOffset(0.73,"y");         // Set y offset
      if(Type==1 && col>=0 || Type==2) Units[pad] = "0.50 GeV^{2})";
      hName = "Events/("; hName += Units[pad];
      hData[pad]->SetYTitle(hName);

      double maxH = hData[pad]->GetMaximum();
      if(maxH<hMC[pad][nPDFs-1]->GetMaximum()) maxH = hMC[pad][nPDFs-1]->GetMaximum();
      hData[pad]->SetMaximum(maxH*1.15);
      hData[pad]->SetLabelOffset(0.01,"Y");
      hData[pad]->GetXaxis()->CenterTitle(true);
      hData[pad]->Draw("");
      for(int pdf=nPDFs-1; pdf>=0; pdf--) {
	hMC[pad][pdf]->Draw("h same");
	if(Colors[pdf] == colDl || Colors[pdf] == colDsl){
	  hCopy[pad][pdf] = (TH1F*)hMC[pad][pdf]->Clone(); hCopy[pad][pdf]->SetFillColor(1); 
	  hCopy[pad][pdf]->SetDirectory(0);
	  if(Colors[pdf] == colDsl) hCopy[pad][pdf]->SetFillStyle(fillDsl); 
	  else hCopy[pad][pdf]->SetFillStyle(fillDl); 
	  hCopy[pad][pdf]->Draw("h same");	  
	}
      } 
      hData[pad]->Draw("axis same");
      hData[pad]->Draw("e0 same");
      double shift = 0;
      if(Type==2){
	PadLabel[pad].ReplaceAll("(",""); PadLabel[pad].ReplaceAll(")","");
	shift = 0.03;
      }
      if(Type==0&&(col==0 || row==1) || Type==1&&(col!=0||row!=1)|| Type==2) 
	label.DrawLatex(LeftMargin+0.105, 0.88-shift, PadLabel[pad]);
      else label.DrawLatex(0.91, 0.88-shift, PadLabel[pad]);

      /////////////////   MC/Data   //////////////////////////////////
      can.cd();
      hName = "Pad1_"; hName += pad;
      //cout<<hName<<": "<<PadXY[0][0]<<", "<<PadXY[1][0]<<", "<<PadXY[0][1]<<", "<<PadXY[1][0]+padH*fRatio<<endl;
      Pads[pad][1] = new TPad(hName,"",PadXY[0][0], PadXY[1][0], PadXY[0][1], PadXY[1][0]+padH*fRatio);
      Pads[pad][1]->SetLeftMargin(LeftMargin); Pads[pad][1]->SetRightMargin(RightMargin); 
      Pads[pad][1]->SetBottomMargin(0.45);      Pads[pad][1]->SetTopMargin(0); 
      Pads[pad][1]->Draw(); Pads[pad][1]->cd();
      hName = "Clone_"; hName += pad;
      hRatio[pad]  = (TH1F*)hData[pad]->Clone(hName);
      hRatio[pad]->SetDirectory(0);
      hRatio[pad]->Divide(hMC[pad][nPDFs-1]);
      int nBins = hRatio[pad]->GetNbinsX();
      double minY = 0.66, maxY = 1.34;
      if(Type>0 && col==2){
	minY = 0.61, maxY = 1.39;
	//hRatio[pad]->GetYaxis()->SetNdivisions(202,false);
      }      
      if(Type==2&&col==2){
	double finVal[] = {12, 11}, minX = hRatio[pad]->GetBinLowEdge(1), maxX = hRatio[pad]->GetBinLowEdge(nBins+1);
	int iniBin = 1, finBin = (int)((finVal[row]-minX)*(double)nBins/(maxX-minX));
	TString FitName = "hFit"; FitName += pad;
	double RealIni = hRatio[pad]->GetXaxis()->GetBinLowEdge(iniBin);
	double RealFin = hRatio[pad]->GetXaxis()->GetBinLowEdge(finBin+1);
	int Fitbins = finBin-iniBin+1;
	hFit[row] = new TH1F(FitName,"",Fitbins,RealIni,RealFin);
	for(int bin = 1; bin<=Fitbins; bin++){
	  hFit[row]->SetBinContent(bin,hRatio[pad]->GetBinContent(bin+iniBin-1));
	  hFit[row]->SetBinError(bin,hRatio[pad]->GetBinError(bin+iniBin-1));
	}
	TString fName = "function"; fName += row;
	//fFit[row] = new TF1(fName,"[5]*x*x*x*x*x+[4]*x*x*x*x+[3]*x*x*x+[2]*x*x+[1]*x+[0]",RealIni,RealFin);
	fFit[row] = new TF1(fName,"[4]*x*x*x*x+[3]*x*x*x+[2]*x*x+[1]*x+[0]",RealIni,RealFin);
	//fFit[row] = new TF1(fName,"[3]*x*x*x+[2]*x*x+[1]*x+[0]",RealIni,RealFin);
	hFit[row]->Fit(fFit[row],"N Q");
	fFit[row]->SetLineWidth(1); fFit[row]->SetLineColor(4);
	double q2 = minX, dq2 = hRatio[pad]->GetXaxis()->GetBinWidth(1);
	for(int bin=1; bin<=nBins; bin++){
	  textFile<<RoundNumber(q2,1)<<"-"<<RoundNumber(q2+dq2,1)<<"  \t"<<RoundNumber(fFit[row]->Integral(q2,q2+dq2),2,dq2);
	  textFile<<"\t"<<RoundNumber(hRatio[pad]->GetBinError(bin),2)<<endl;
	  q2 += dq2;
	}
	textFile<<endl;
	FitName = "hFitError"; FitName += pad;
	hFit[row+2] = new TH1F(FitName,"",nBinsError,RealIni,RealFin);
	hFit[row+2]->SetDirectory(0);
	for(int bin=1; bin<=nBinsError; bin++){
	  q2 = hFit[row+2]->GetBinCenter(bin);
	  int binH = (int)((q2-minX)/dq2+1);
	  double errBin = hRatio[pad]->GetBinError(binH), val = fFit[row]->Eval(q2);
	  double loError = val-errBin, hiError = val+errBin;
	  if(loError > 1-errBin) loError = 1-errBin;
	  if(hiError < 1+errBin) hiError = 1+errBin;
	  val = (hiError+loError)/2.;
	  hFit[row+2]->SetBinContent(bin, val);
	  hFit[row+2]->SetBinError(bin, hiError-val);
	}
	hFit[row+2]->SetFillColor(kBlue-10); hFit[row+2]->SetLineColor(kBlue-9); 
      }
      for(int bin=1; bin<=nBins; bin++) {
	double val = hRatio[pad]->GetBinContent(bin);
	if(val>maxY || val<minY) {hRatio[pad]->SetBinContent(bin, -99);hRatio[pad]->SetBinError(bin, 0);}
      }

      hName = "x"; if(Type==0) hName += pad; 
      hRatio[pad]->SetXTitle(hName);
      hRatio[pad]->SetTitleSize(TextSizeR*1.15,"xy");      // Set the 2 axes title size
      hRatio[pad]->SetLabelSize(TextSizeR,"xy");      // Set the 2 axes label size
      hRatio[pad]->SetYTitle("data/MC");
      hRatio[pad]->GetXaxis()->SetNdivisions(505);         // 5 primary ticks and 4 secondary ticks
      //if(fName[Type][pad].Contains("mES")) hRatio[pad]->GetXaxis()->SetNdivisions(1006);
      //hRatio[pad]->GetYaxis()->SetNdivisions(201);         // 5 primary ticks and 4 secondary ticks
      hRatio[pad]->SetTitleOffset(0.48,"y");         // Set y offset
      hRatio[pad]->SetMinimum(minY); hRatio[pad]->SetMaximum(maxY); 
      hRatio[pad]->GetYaxis()->SetRangeUser(minY,maxY); 
      hRatio[pad]->Draw("E0 P");
      if(Type==2&&col==2){
	hFit[row+2]->Draw("e3 same");
	fFit[row]->Draw("same");
      }
      TLine line; line.SetLineColor(2);
      line.DrawLine(hRatio[pad]->GetXaxis()->GetBinLowEdge(1), 1.0, 
		    hRatio[pad]->GetXaxis()->GetBinLowEdge(nBins+1), 1.0);


      hRatio[pad]->Draw("E0 P same");
      hRatio[pad]->Draw("axis same");
    }
  }
  can.cd(0);
  TBox box; box.SetFillStyle(1001); box.SetLineColor(10); box.SetFillColor(10);
  label.SetTextAngle(0); label.SetTextAlign(13); label.SetTextSize(TextSize/3.5); 
  for(int row=0; row<nRows; row++){
    for(int col=0; col<nCols; col++){
      box.DrawBox(padW*(LeftMargin+col)-0.03,  padH*(dRows-row-1+fRatio), 
		  padW*(LeftMargin+col)-0.002, padH*(dRows-row-1+fRatio)+0.02);
      label.DrawLatex(padW*(LeftMargin+col)-0.01, padH*(dRows-row-1+fRatio)+0.014,"0");
    }
  }

  if(Type<2){
    int nEntries=0;
    TString legLabel[] = {"ccbar", "uds", "BB", "xfeed", "Dssl", "Dsl", "Dl"};
    double legXY[2][2] = {{padW*0.62, padW*0.87}, {padH*1.55, padH*1.95}};
    if(Type==1) {
      legXY[0][0] = padW*0.21; legXY[0][1] = padW*0.43; 
      legXY[1][0] = padH*0.55; legXY[1][1] = padH*0.95; }
    TLegend leg(legXY[0][0], legXY[1][0], legXY[0][1], legXY[1][1]);
    leg.SetTextSize(0.04); leg.SetFillStyle(0); leg.SetTextFont(132); leg.SetBorderSize(0);
    for(int pdf=nPDFs-1; pdf>0; pdf--) {
      if(pdf==3) continue;
      int ipdf = pdf;
      //if(pdf==5) ipdf = 6;
      //if(pdf==6) ipdf = 5;
      leg.AddEntry(hMC[1][ipdf], legLabel[pdf]);
      nEntries++;
    }
    double iDl = 4, iDsl = 3;
    leg.Draw();
    fillDl -= 400; fillDsl -= 400; 
    box.SetFillStyle(fillDl); box.SetLineColor(0); box.SetFillColor(1);
    double legW = legXY[0][1]-legXY[0][0], legH = (legXY[1][1]-legXY[1][0])/(double)nEntries;
    box.DrawBox(legXY[0][0]+legW*0.04, legXY[1][0]+legH*(iDl+0.15), 
		legXY[0][0]+legW*0.21, legXY[1][0]+legH*(iDl+0.85));
    box.SetFillStyle(fillDsl); 
    box.DrawBox(legXY[0][0]+legW*0.04, legXY[1][0]+legH*(iDsl+0.15), 
		legXY[0][0]+legW*0.21, legXY[1][0]+legH*(iDsl+0.85));


    box.SetFillColor(0); box.SetFillStyle(0); box.SetLineWidth(1); box.SetLineColor(1); 
    for(int row=0; row<nEntries; row++)
      box.DrawBox(legXY[0][0]+legW*0.04, legXY[1][0]+legH*(row+0.15), 
		  legXY[0][0]+legW*0.21, legXY[1][0]+legH*(row+0.85));
  } else {
    textFile.close();
    cout<<"Written "<<textName<<endl;
  }
  TString pName = "public_html/"; pName += TypeName[Type]; pName += ".eps"; 
  can.SaveAs(pName);

  for(int row=0; row<nRows; row++){
    for(int col=0; col<nCols; col++){
      //int pad = col+3*row;
      //for(int isRat=0; isRat<2; isRat++) if(Pads[pad][isRat]) Pads[pad][isRat]->Delete();
    }
  }
  
}

