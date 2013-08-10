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
#include "TH1F.h"
#include "babar_code/Styles/Styles.cc"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;
#define nFiles 9

void RD_BDT(TString FileTag = "mES", bool isError=false, int isIso=0){

  Styles style; style.setPadsStyle(2); style.applyStyle();
  //gStyle->SetPadLeftMargin(0.3); gStyle->SetPadTopMargin(0.12);

  int colors[] = {2,4,28,9, 46, 38};
  double RD[nFiles][4][2];
  double minRD[2] = {99., 99.}, maxRD[2] = {-99., -99.}; 

  double nomiRD[2][2] = {{0.456, 0.053}, {0.325, 0.022}};
  TString folder = "FitAll/fits/TextFinalRAll"; if(isIso==1) folder = "FitAll/fits/TextFinalIsoData"; 
  folder+=FileTag; folder+="x";
  TString nameFiles[] = {"030", "050", "075", "100", "120", "150", "200", "250", "300", "400"};

  TString RDTag[] = {"R(D^{0})", "R(D*^{0})", "R(D^{+})", "R(D*^{+})"};

  if(FileTag.Contains("DD")) folder.ReplaceAll(FileTag,"Df2");
  for(int file=0; file<nFiles; file++){
    TString fileName = folder; fileName += nameFiles[file]; fileName += ".txt";
    fstream textFile; textFile.open(fileName,fstream::in);
    TString text = ""; int times = 0;
    while(!text.Contains("R(D")) {
      textFile>>text; 
      times++;
      if(times>1000) {cout<<"R(D) not found in "<<fileName<<endl; break;}
    }
    int endChan = 4; if(isIso==1) endChan = 2;
    for(int cand=0; cand<endChan; cand++){
      if(cand>0) textFile>>text;
      textFile>>text; RD[file][cand][0] = text.Atof();
      textFile>>text; textFile>>text; text.ReplaceAll(",",""); RD[file][cand][1] = text.Atof();
      for(int i=0; i<12; i++) textFile>>text; 
      if(RD[file][cand][0]-RD[file][cand][1] < minRD[cand%2]) minRD[cand%2] = RD[file][cand][0]-RD[file][cand][1];
      if(RD[file][cand][0]+RD[file][cand][1] > maxRD[cand%2]) maxRD[cand%2] = RD[file][cand][0]+RD[file][cand][1];
      //cout<<RD[file][cand][0]<<" #pm "<<RD[file][cand][1]<<endl;
    }
    //cout<<endl;
  }
  if(FileTag.Contains("DD")) {
    folder.ReplaceAll("Df2","Ds2");
    for(int file=0; file<nFiles; file++){
      TString fileName = folder; fileName += nameFiles[file]; fileName += ".txt";
      fstream textFile; textFile.open(fileName,fstream::in);
      TString text = ""; int times = 0;
      while(!text.Contains("R(D")) {
	textFile>>text; 
	times++;
	if(times>1000) {cout<<"R(D) not found in "<<fileName<<endl; break;}
      }
      int endChan = 4; if(isIso==1) endChan = 2;
      for(int cand=0; cand<endChan; cand++){
	if(cand>0) textFile>>text;
	textFile>>text; if(cand%2==0) RD[file][cand][0] = text.Atof();
	textFile>>text; textFile>>text; text.ReplaceAll(",",""); if(cand%2==0) RD[file][cand][1] = text.Atof();
	for(int i=0; i<12; i++) textFile>>text; 
	if(RD[file][cand][0]-RD[file][cand][1] < minRD[cand%2]) minRD[cand%2] = RD[file][cand][0]-RD[file][cand][1];
	if(RD[file][cand][0]+RD[file][cand][1] > maxRD[cand%2]) maxRD[cand%2] = RD[file][cand][0]+RD[file][cand][1];
      }
    }
  }

//   if(maxRD[1] > maxRD[0]) maxRD[0] = maxRD[1];
//   if(minRD[1] < minRD[0]) minRD[0] = minRD[1];

  TString fileName = folder; fileName += "100.txt"; if(isIso==0) fileName.ReplaceAll("Data","IsoData");
  fstream textFile; textFile.open(fileName,fstream::in);
  TString text = ""; int times = 0;
  while(!text.Contains("R(D")) {
    textFile>>text; 
    times++;
    if(times>1000) {cout<<"R(D) not found in "<<fileName<<endl; break;}
  }
  for(int cand=0; cand<2; cand++){
    if(cand>0) textFile>>text;
    textFile>>text; nomiRD[cand][0] = text.Atof();
    textFile>>text; textFile>>text; text.ReplaceAll(",",""); nomiRD[cand][1] = text.Atof();
    for(int i=0; i<12; i++) textFile>>text; 
  }
  //nomiRD[0][0] = 0.415; 
  double minBDT = 0, maxBDT = 330;
  TCanvas can("can","RD results");
  can.Divide(2,1); TPad *cPad;
  TString yTitle[] = {"R(D)", "R(D*)"};
  TH1F* hRD[2];
  TBox box; box.SetFillStyle(3002); box.SetLineColor(0);box.SetFillColor(34);
  TLatex label; label.SetTextSize(0.08); label.SetTextFont(132); label.SetTextAlign(12);
  TMarker mark; mark.SetMarkerStyle(8); mark.SetMarkerSize(1.05);
  TLine lin; lin.SetLineStyle(1); lin.SetLineWidth(2);
  double minRD2[2] = {0.19, 0.21}, maxRD2[2] = {0.76, 0.48}; 
  if(isError) {
    minRD2[0] = 8; minRD2[1] = 5; maxRD2[0] = 32; maxRD2[1] = 21;
    yTitle[0] = "Error (%)"; yTitle[1] = "Error (%)"; 
    if(isIso!=0){maxRD2[0] = 21; maxRD2[1] = 12;}
  } else if(isIso==1) {
    minRD2[0] = 0.24; minRD2[1] = 0.24; maxRD2[0] = 0.59; maxRD2[1] = 0.39;
  }
  for(int pad=0; pad<2; pad++){
    cPad = (TPad *)can.cd(pad+1);

    double delta = maxRD[pad] - minRD[pad];
    TString hName = "Ratio"; hName += pad+1;
    hRD[pad] = new TH1F(hName, "", 10, minBDT, maxBDT);
    //hRD[pad]->GetYaxis()->SetNdivisions(0);
    hRD[pad]->SetMaximum(maxRD[pad] + 0.1*delta);
    hRD[pad]->SetMinimum(minRD[pad] - 0.1*delta);
    hRD[pad]->SetMaximum(maxRD2[pad]);
    hRD[pad]->SetMinimum(minRD2[pad]);
    hRD[pad]->Draw();
    if(!isError) {box.SetFillColor(colors[4+pad]); box.DrawBox(minBDT, nomiRD[pad][0]-nomiRD[pad][1], maxBDT, nomiRD[pad][0]+nomiRD[pad][1]);}
    if(!isError) {lin.SetLineColor(colors[4+pad]); lin.DrawLine(minBDT, nomiRD[pad][0], maxBDT, nomiRD[pad][0]);}
    for(int file=0; file<nFiles; file++){
      int nBeg = 0; if(isIso==3) nBeg=2;
      if(pad==1) nBeg++;
//       label.SetTextFont(22); label.SetTextColor(colors[file]); label.SetTextAlign(12);
//       label.SetTextSize(0.072); label.DrawLatex(minRD[pad]-0.6*delta, 1-halfH*(2*file+1)+0.085, fileTag[file]);
      int endChan = 4; if(isIso==1 || isIso==2) endChan = 2;
      for(int cand=nBeg; cand<endChan; cand+=2){
	double height = RD[file][cand][0], dH = RD[file][cand][1], BDT = nameFiles[file].Atof()-4;
	if(cand>1) BDT += 8;
	//double height = 1 - (halfH*(2*file+1) - halfH/3.2);
	//if(cand>1) height = 1 - (halfH*(2*file+1) + halfH/3.2);
	mark.SetMarkerColor(colors[cand]); lin.SetLineColor(colors[cand]);
	if(isError) height = RD[file][cand][1]/RD[file][cand][0]*100;
	if(file>0 || 1){
	  mark.DrawMarker(BDT, height);
	  if(!isError) lin.DrawLine(BDT, height-dH, BDT, height+dH);
	}
	if(file==0){
	  double delta = 25; if(isError) delta = 70;
	  height = minRD2[pad] + (maxRD2[pad]-minRD2[pad])*(0.93-0.1*(cand>1));
	  mark.DrawMarker(25+delta, height);
	  label.SetTextSize(0.075); label.SetTextColor(colors[cand]);label.SetTextAlign(12);
	  label.DrawLatex(34+delta, height*1.005, RDTag[cand]);
	  
	}
      }
      style.setTitles(hRD[pad], "% of the nominal sample", yTitle[pad]);
    }
//     TString rightLabel = "#splitline{#chi^{2}: "; rightLabel += RoundNumber(chi2[pad],1); rightLabel += "/";
//     rightLabel += ndof[pad]; rightLabel += " = "; rightLabel += RoundNumber(chi2[pad],2,ndof[pad]);
//     rightLabel += "}{Prob. = "; rightLabel += RoundNumber(TMath::Prob(chi2[pad],ndof[pad])*100,2); rightLabel += "%}";
//     style.setTitles(hRD[pad],xTitle[pad]);//,"","",rightLabel);
//     label.SetTextSize(0.07); label.SetTextColor(1);label.SetTextAlign(23);
//     label.DrawLatex(SMRD[pad][0], 1.12, "SM");
//     label.DrawLatex(nomiRD[pad][0], 1.12, "Aver.");
  }

  TString pName = "public_html/RD_BDT_"; pName+=FileTag; if(isIso) pName+="Iso"; if(isError) pName+="_Err"; 
  pName+=".eps"; 
  can.SaveAs(pName);
  for(int pad=0; pad<2; pad++){
      hRD[pad]->Delete();
  }
}

