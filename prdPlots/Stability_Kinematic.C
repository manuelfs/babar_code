#include "TMath.h"
#include "TString.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TH1F.h"
#include "TFile.h"
#include "THStack.h"
#include "TChain.h"
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

#define nPads 9
#define nType 10
#define nFont 132
using namespace TMath;
using namespace std;
using std::cout;
using std::endl;

void formatHisto(TH1F *histo, int color, double TextSize);

void Stability_Kinematic(int Type = 9){
  TString TypeName[] = {"EExtra", "Kinematic", "Q2", "Q2_Bkg", "Q2_BDT", "Q2_BDTBkg", 
			"Q2_BDTBkg2", "Q2_Norm", "Q2_Efficiency", "Q2_NormEffDss"};
  if(Type==-1) {
    for(int type=0; type<nType; type++)
      cout<<type<<" "<<TypeName[type]<<", "; 
    cout<<endl; 
    return;
  }
  int doRoot=0; if(Type==2) doRoot = 1;
  TString histosName = "keys/eps/CSample/root/Q2_Folded.root";
  TFile *fHistos2=0; if(doRoot) fHistos2 = new TFile(histosName,"RECREATE");

  TString Folder = "keys/eps/CSample/root/";
  //TString PadLabel[] = {"(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"};
  TString PadLabel[] = {"a", "b", "x", "d", "e", "f", "(g)", "(h)", "(i)"};
  TString fName[nType][3] = {{"EExtra",       "EExtra",    "EExtra_Sub"},
			     {"MES",          "MES_Sub",   "Pl_Sub"},
			     {"Q2_Sub",       "Q2_Sub030", "Q2_Sub045"},
			     {"Q2",           "Q2_030",    "Q2_045"},
			     {"Q2_SubBDT050", "Q2_Sub",    "Q2_SubBDT150"},
			     {"Q2_BDT050",    "Q2_BDT100", "Q2_BDT150"},
			     {"Q2_BDT030",    "Q2_BDT050x2", "Q2_BDT075"},
			     {"Q2_Sub",       "Q2",        "Q2_Norm"},
			     {"Q2_Unfolded",  "Q2_Unfolded", "Q2_Unfolded"},
			     {"Q2_Norm",  "Q2_Norm",  "Q2_Norm"}};
  TString Units[] = {"1 MeV)", "1 MeV)", "100 MeV)", "50 MeV)", "0.50 GeV^{2})", "0.35 GeV^{2})", "0.25 GeV^{2})"}, hName;
  TString xTitle[] = {"M", "M", "P", "E", "Q"};

  int nPDFs[nType][3] = {{9,9,2}, {9,2,2}, {2,2,2}, {9,9,9}, {2,2,2}, 
			 {9,9,9}, {9,9,9}, {2,9,9}, {1,1,1}, {9,9,9}};
  int nRows = 2, nCols = 3, colBkg=kYellow-7, nBinsE=17;
  int Colors[] = {colBkg,colBkg,colBkg,colBkg,28,kCyan+1,kBlue-8,kRed-7,kGreen-7}, ColorsE[] = {8,4,2};
  TLatex label; label.SetTextFont(132); label.SetNDC(kTRUE);
  TCanvas can("dataMC","data Vs MC",700,150*nRows); 
  TPad *Pads[nPads][2];
  TH1F *hData[nPads], *hYields[nPads], *hMC[nPads][10], *hTemp, *hCopy[nPads][10], *hDssPDF[3][2];
  TH1F *hQ2[nPads][3], *hEff[nPads][3];
  double limQ2[] = {4, 12.5};

  TChain MC("ntp1");
  MC.Add("AWG82/ntuples/small/FitRAllNewx100_RunAll.root");
  TString mcCuts[6] = {"weight*((MCType==5&&candType==1||MCType==11&&candType==3)&&candM2>1.5)",
		       "weight*((MCType==6&&candType==2||MCType==12&&candType==4)&&candM2>1.5)" ,
		       "weight*((MCType==5&&candType==1||MCType==11&&candType==3))",
		       "weight*((MCType==6&&candType==2||MCType==12&&candType==4))",
		       "weight*(((MCType==1||MCType==3)&&candType==1||(MCType==7||MCType==9)&&candType==3))",
		       "weight*(((MCType==2||MCType==4)&&candType==2||(MCType==8||MCType==10)&&candType==4))"};

  if(Type==9){ // Reading D**(l/tau)nu histograms
    TString textName = "public_html/Uncertainties_BkgDss.txt";
    fstream textFile; textFile.open(textName,fstream::in);
    for(int row=0; row<nRows; row++){
      for(int dss=0; dss<3; dss++){
	hName = "hDssPDF"; hName += dss; hName += row;
	hDssPDF[dss][row] = new TH1F(hName, "", nBinsE, limQ2[0], limQ2[1]);
	hDssPDF[dss][row]->SetDirectory(0);
	hDssPDF[dss][row]->SetLineWidth(2); 
      }
      for(int bin=1; bin<nBinsE; bin++){
	textFile>>hName;
	double val[2];
	for(int dss=0; dss<3; dss++){
	  textFile>>val[0]>>val[1];
	  hDssPDF[dss][row]->SetBinContent(bin, val[0]);
	  //hDssPDF[dss][row]->SetBinError(bin, val[1]);
	  //hDssPDF[dss][row]->SetBinError(bin, 0);
	}
      }
    }
  }
  gStyle->SetOptStat(0);
  gStyle->SetHatchesLineWidth(1);
  int fillDl = 3654, fillDsl = 3645;
  TBox box; box.SetLineColor(10);box.SetFillColor(10);
  double dRows = (double)nRows, dCols = (double) nCols;
  double bMargin = 0.12, padH = (1-bMargin)/dRows, padW = 1/dCols, LeftMargin = 0.2;
  for(int col=0; col<nCols; col++){
    TString histosName = Folder; histosName += fName[Type][col]; histosName += ".root";
    TFile fHistos(histosName);
    //cout<<"Opened "<<histosName<<endl;
    for(int row=0; row<nRows; row++){
      int onePDF = 0; if(Type==2 || Type==4 || Type==7&&col==0) onePDF = 1;
      int pad = nRows*col+row;
      int candi = 2+row;
      if(nPDFs[Type][col]==9 && col==0 && Type==0) candi = row;
      fHistos.cd(); can.cd(0);
      hName = "Data"; hName += candi;
      if(Type == 8) {hName = "UnfoldedData_"; hName += pad;}
      hTemp = (TH1F *)fHistos.Get(hName);
      hName = "Data_"; hName += pad;
      if(hTemp) hData[pad]  = (TH1F*)hTemp->Clone(hName);
      else continue;
      hData[pad]->SetDirectory(0);
      int nBins = hData[pad]->GetNbinsX();
      hName = "hYields";  hName += pad;
      hTemp = (TH1F *)fHistos.Get("hYields");
      if(hTemp && Type<8 || Type==9) hYields[pad]  = (TH1F*)hTemp->Clone(hName);
      else if(Type<8) continue; 
      for(int pdf=0; pdf<nPDFs[Type][col]; pdf++){
	int ipdf = pdf+(pdf==(nPDFs[Type][col]-2))+(nPDFs[Type][col]==2)*7;
	if(pdf==(nPDFs[Type][col]-1)) ipdf = 7;
	if(onePDF){
	  if(pdf==0&&row==0) ipdf=7;
	  if(pdf==1&&row==0) ipdf=8;
	}
	if(Type==8&&row==0) ipdf = 7;
	if(Type==8&&row==1) ipdf = 8;
	hName = "MC"; hName += candi+1; hName += ipdf;
	double yield = 0;
	if(Type<8 || Type==9){
	  if(!(col==1&&Type==0)) yield += hYields[pad]->GetBinContent(ipdf+10*row);
	  if(!(col==0&&Type==0)) yield += hYields[pad]->GetBinContent(ipdf+10*row+20);
	} else {
	  yield = 1;
	  hName = "UnfoldedMC_"; hName += pad;
	}
	if(yield<=0){
	  cout<<pdf<<": ipdf "<<ipdf<<", yield "<<yield<<", row "<<row<<", col "<<col<<endl;
	  hMC[pad][pdf]  = (TH1F*)hMC[pad][pdf-1]->Clone(hName);
	} else {
	  hTemp = (TH1F *)fHistos.Get(hName);
	  hName = "MC_";  hName += pad; hName += pdf;
	  if(hTemp) hMC[pad][pdf]  = (TH1F*)hTemp->Clone(hName);
	  else {cout<<row<<", "<<col<<": PDF "<<pdf<<endl; return;}
	  if(hMC[pad][pdf]->Integral()<=0) {
	    cout<<pdf<<": ipdf "<<ipdf<<", yield "<<yield<<", row "<<row<<", col "<<col<<", candi "<<candi<<endl;
	    hName = "MC"; hName += candi+1-2; hName += ipdf;
	    hTemp = (TH1F *)fHistos.Get(hName);
	    hMC[pad][pdf]->Add(hTemp);
	  }
	}
	hMC[pad][pdf]->SetDirectory(0);
	if(Type<8 || Type==9) hMC[pad][pdf]->Scale(yield/hMC[pad][pdf]->Integral());
	hMC[pad][pdf]->SetFillColor(Colors[ipdf]);
	hMC[pad][pdf]->SetLineColor(Colors[ipdf]); hMC[pad][pdf]->SetLineWidth(0);
	if(pdf==1){
	  hMC[pad][pdf]->SetLineWidth(1); hMC[pad][pdf]->SetLineColor(1);
	  hMC[pad][pdf]->SetLineStyle(2);
	}
	if(pdf>0) {
	  if(onePDF==0) hMC[pad][pdf]->Add(hMC[pad][pdf-1]);
	  else hData[pad]->Add(hMC[pad][pdf],-1);
	}
      }
      if(onePDF) nPDFs[Type][col] = 1;
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

      hData[pad]->SetMarkerStyle(20); hData[pad]->SetMarkerSize(0.5);
      hData[pad]->SetTitleSize(TextSize,"xy");  // Set the 2 axes title size
      hData[pad]->SetLabelSize(TextSize,"xy");      // Set the 2 axes label size
      hData[pad]->SetTitleFont(nFont,"xy");         // Set the all 2 axes title font
      hData[pad]->SetLabelFont(nFont,"xy");         // Set the all 2 axes label font
      hData[pad]->SetNdivisions(504, "xy");         // 5 primary ticks and 4 secondary ticks
      hData[pad]->SetLabelOffset(0.01,"Y");
      hData[pad]->GetXaxis()->CenterTitle(true);
      if(row==1) {
	if(Type==0) hName = xTitle[3];
	else if(Type==1) hName = xTitle[col];
	else hName = xTitle[4];
      } else hName = "";
      hData[pad]->SetXTitle(hName);

      // Chi2 calculation
      double chi2 = 0; int ndof = -2, digits = 1;
      for(int bin=1; bin<=nBins; bin++) {
	double vBin = hData[pad]->GetBinContent(bin)-hMC[pad][nPDFs[Type][col]-1]->GetBinContent(bin);
	double eBin = sqrt(pow(hData[pad]->GetBinError(bin),2)+
			   pow(hMC[pad][nPDFs[Type][col]-1]->GetBinError(bin),2));
	if(vBin != 0 && hData[pad]->GetBinError(bin) >= sqrt(8)){
	  chi2 += pow(vBin/eBin,2);
	  ndof++;
	}
      }
      double maxH = hData[pad]->GetMaximum();
      if(maxH<hMC[pad][nPDFs[Type][col]-1]->GetMaximum()) maxH = hMC[pad][nPDFs[Type][col]-1]->GetMaximum();
      hData[pad]->SetMaximum(maxH*1.28);
      if(Type>1) hData[pad]->SetMaximum(maxH*1.6);
      if(Type==0) hData[pad]->SetMaximum(maxH*1.2);
      if(Type==9) hData[pad]->SetMaximum(maxH*1.39);
      if(nPDFs[Type][col]==1) {
	double minY = 1.4*hData[pad]->GetMinimum();
	if(minY>-maxH/12.) minY = -maxH/12.;
	hData[pad]->SetMinimum(minY);
	if(row==0) hData[pad]->SetMaximum(maxH*1.7);
      }
      if(Type==2||Type==8) hData[pad]->SetMaximum(80);
      if(Type==2||Type==8 && row==0) hData[pad]->SetMinimum(-8);
      if(nPDFs[Type][col]==2) hData[pad]->SetMinimum(-maxH/12.);
      if(nPDFs[Type][col]==2 && Type==0) hData[pad]->SetMinimum(-maxH/20.);
      hData[pad]->Draw("");
      for(int pdf=nPDFs[Type][col]-1; pdf>=0; pdf--) {
	for(int bin=1; bin<=nBins; bin++) {
	  hMC[pad][pdf]->SetBinError(bin, 0);
	}
	hMC[pad][pdf]->Draw("h same");
	if(pdf==5 || pdf==6){
	  hName = "Copy"; hName += pad; hName += pdf;
	  hCopy[pad][pdf] = (TH1F*)hMC[pad][pdf]->Clone(hName); hCopy[pad][pdf]->SetFillColor(1); 
	  hCopy[pad][pdf]->SetDirectory(0);
	  if(pdf==5) hCopy[pad][pdf]->SetFillStyle(fillDl); 
	  else hCopy[pad][pdf]->SetFillStyle(fillDsl); 
	  hCopy[pad][pdf]->Draw("h same");
	}
      }
      hData[pad]->Draw("axis same");
      hData[pad]->Draw("e0 same");
      if(doRoot) {
	fHistos2->cd(); 
	hData[pad]->Write(); hMC[pad][nPDFs[Type][col]-1]->Write(); 
      }

      if(Type==9 && col==1){
	BToDtaunu Dtaunu; BToDstaunu Dstaunu;
	double ml[] = {Dtaunu.mTau, Dtaunu.mTau, Dtaunu.mMu}, scaleEff = 1;
	int nValues = 50, markerStyle[] = {8,4,21};
	for(int isNorm=2; isNorm>=0; isNorm--){
	  hName = "Eff"; hName += pad; hName += isNorm;
	  hEff[pad][isNorm] = new TH1F(hName,"",nBinsE, limQ2[0], limQ2[1]);
	  formatHisto(hEff[pad][isNorm], ColorsE[isNorm], TextSize);
	  MC.Project(hName,"candQ2",mcCuts[row+2*isNorm]);
	  if(row==1) hEff[pad][isNorm]->SetXTitle(xTitle[4]);
	  hEff[pad][isNorm]->SetDirectory(0);
	  hEff[pad][isNorm]->SetMarkerStyle(markerStyle[isNorm]);

	  //cout<<pad<<" Calculating q2"<<isNorm<<endl;
	  hName = "hQ2"; hName += pad; hName += isNorm;
	  hQ2[pad][isNorm] = new TH1F(hName,"",nBinsE, limQ2[0], limQ2[1]);
	  hQ2[pad][isNorm]->SetLineColor(Colors[isNorm]);
	  for(int bin=1; bin<=nBinsE; bin++){
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
	  hEff[pad][isNorm]->Divide(hQ2[pad][isNorm]);
	  if(isNorm==2) scaleEff = 1/hEff[pad][isNorm]->GetMaximum();
	  if(isNorm==1) scaleEff = hEff[pad][2]->Integral(nBinsE/3,nBinsE/3*2)/
	    hEff[pad][1]->Integral(nBinsE/3,nBinsE/3*2);
	  hEff[pad][isNorm]->Scale(scaleEff);
	  hEff[pad][isNorm]->SetMinimum(0);
	  if(row==1) hEff[pad][isNorm]->SetMaximum(1.72);
	}
	//cout<<pad<<" Drawing "<<scaleEff<<endl;
	hEff[pad][1]->Draw();
	hEff[pad][0]->Draw("same"); hEff[pad][2]->Draw("same"); hEff[pad][1]->Draw("same");

      }
      if(Type==9 && col==2){    
	ColorsE[0] = 4; ColorsE[1] = 2; ColorsE[2] = 8; 
	int lineStyle[] = {1,3,2};
	for(int dss=0; dss<3; dss++){
	  formatHisto(hDssPDF[dss][row], ColorsE[dss], TextSize);
	  hDssPDF[dss][row]->SetLineStyle(lineStyle[dss]); hDssPDF[dss][row]->SetMarkerSize(0);
	  hDssPDF[dss][row]->Scale(100/hDssPDF[dss][row]->Integral());
	  if(dss==0){ 
	    if(row==0) hDssPDF[dss][row]->SetMaximum(hDssPDF[dss][row]->GetMaximum()*1.65);
	    else hDssPDF[dss][row]->SetMaximum(hDssPDF[dss][row]->GetMaximum()*1.4);
	    hDssPDF[dss][row]->Draw("hist");
	    if(row==1)hDssPDF [dss][row]->SetXTitle(xTitle[4]);
	  } else hDssPDF[dss][row]->Draw("hist same");
	}
      }


      // Labels
      if((Type==2||Type==3||Type==8) && col==2 && row==0) digits = 4;
      if((Type==4||Type==5) && col==0 && row==0) digits = 2;
      label.SetTextAngle(0); label.SetTextSize(TextSize/0.93); 
      TString chi2Label = "#splitline{#chi^{2}: ", chanLabel = PadLabel[pad]; 
      chi2Label += RoundNumber(chi2,1); chi2Label += "/"; chi2Label += ndof; chi2Label += "}{p = "; 
      chi2Label += RoundNumber(TMath::Prob(chi2,ndof)*100,digits); 
      chi2Label += "%}";
      if(nPDFs[Type][col]==9 && Type<2) chi2Label = "";
      if(Type>1) {
	chi2Label.ReplaceAll("#splitline{","");
	chi2Label.ReplaceAll("}{",", ");
	chi2Label.ReplaceAll("%}","%");
      }
      double Xchi2 = LeftMargin+0.04, Ychi2 = 0.94, Xchan = 0.92, Ychan = 0.89;
      if(Type<2){
	Xchi2 = 0.93; label.SetTextAlign(33);
      } else label.SetTextAlign(13); 
      if(!(Type==9&&col>0)) label.DrawLatex(Xchi2, Ychi2, chi2Label);

      if(Type==0 && (col==0&&row==0 || col==2)) Xchan = LeftMargin+0.3;
      else if(Type==0) Xchan = 0.92;
      else if(Type<2)  Xchan = LeftMargin+0.06;
      label.SetTextAlign(22); label.DrawLatex(Xchan, Ychan, chanLabel);
      if(onePDF) nPDFs[Type][col] = 2;
    } 
    can.cd(0);
    hName = "Events/("; 
    if(Type==8) hName = "Weighted events/("; 
    if(Type==9&&col==1) hName = "Normalized efficiency/("; 
    if(Type==0) hName += Units[3];
    else if(Type==1) hName += Units[col];
    else if(Type<6 || Type==7&&col<2 || Type==8 || Type==9&&col>=1) hName += Units[4];
    else if(Type==7&&col==2 || Type==9&&col==0) hName += Units[6];
    else hName += Units[5];
    label.SetTextAlign(23); label.SetTextAngle(90); label.SetTextSize(0.055);
    label.DrawLatex(0.01+col/dCols, 0.56, hName);
    if(nPDFs[Type][col]==9) {
      box.DrawBox(padW*(LeftMargin+col)-0.02, bMargin+padH, padW*(LeftMargin+col)-0.002, bMargin+padH+0.02);
      label.SetTextAlign(32); label.SetTextAngle(0); label.SetTextSize(0.044);
      label.DrawLatex(padW*(LeftMargin+col)-0.0025, bMargin+padH+0.003, "0");
    }

  }

  if(doRoot) {fHistos2->Close(); fHistos2->Delete();}
  int nEntries=0;
  TString legLabel[] = {"Bkg", "Bkg", "Bkg", "Bkg", "Dssl", "Dsl", "Dl", "Dstau", "Dtau"};
  double legXY[2][2] = {{padW*0.62, padW*0.87}, {bMargin+padH*1.28, bMargin+padH*1.95}};
  TLegend leg(legXY[0][0], legXY[1][0], legXY[0][1], legXY[1][1]);
  if(Type==0) {
    can.cd(0);
    leg.SetTextSize(0.05); leg.SetFillColor(0); leg.SetTextFont(132); leg.SetBorderSize(0);
    for(int pdf=8; pdf>1; pdf--) {
      if(pdf==3) continue;
      int ipdf = pdf;
      if(pdf==5) ipdf = 6;
      if(pdf==6) ipdf = 5;
      leg.AddEntry(hMC[0][ipdf], legLabel[pdf]);
      nEntries++;
    }
    double iDl = 3, iDsl = 2;
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
  }
  double leg9XY[2][2] = {{padW*1.26, padW*1.7}, {bMargin+padH*1.05, bMargin+padH*1.42}};
  TLegend leg9(leg9XY[0][0], leg9XY[1][0], leg9XY[0][1], leg9XY[1][1]);
  //double legDssXY[2][2] = {{padW*2.45, padW*2.7}, {bMargin+padH*1.05, bMargin+padH*1.42}};
  double legDssXY[2][2] = {{padW*2.24, padW*2.55}, {bMargin+padH*1.6, bMargin+padH*1.95}};
  TLegend legDss(legDssXY[0][0], legDssXY[1][0], legDssXY[0][1], legDssXY[1][1]);
  if(Type==9){
    can.cd(0);
    legLabel[0] = "Sigm2"; legLabel[1] = "Sig"; legLabel[2] = "Norm";
    leg9.SetTextSize(0.042); leg9.SetFillStyle(0); leg9.SetTextFont(132); leg9.SetBorderSize(0);
    for(int isNorm=2; isNorm>=0; isNorm--)
      leg9.AddEntry(hEff[2][isNorm], legLabel[isNorm]);
    leg9.Draw();

    legLabel[0] = "Dss"; legLabel[1] = "Dsstau"; legLabel[2] = "Dpipi";
    legDss.SetTextSize(0.042); legDss.SetFillStyle(0); legDss.SetTextFont(132); legDss.SetBorderSize(0);
    for(int dss=2; dss>=0; dss--)
      legDss.AddEntry(hDssPDF[dss][0], legLabel[dss]);
    legDss.Draw();
  }


  TString pName = "public_html/Stability_"; pName += TypeName[Type]; pName += ".eps"; 
  if(Type==2) pName.ReplaceAll("Stability","Higgs");
  if(Type>2 &&Type<8) pName.ReplaceAll("Stability","Test");
  if(Type==2 || Type==8) pName.ReplaceAll("Stability","Higgs");
  can.SaveAs(pName);
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
