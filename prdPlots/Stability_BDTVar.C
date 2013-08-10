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

void Stability_BDTVar(TString FileTag = "Ne2", bool isError=false, int isIso=0){

  Styles style; style.setPadsStyle(-8); 
  style.CanvasH = 350; 
  style.applyStyle();

  int colors[] = {4,4,2,2, kGreen+2, kSpring-4};
  double RD[nFiles][4][2];
  double minRD[2] = {99., 99.}, maxRD[2] = {-99., -99.}; 

  double nomiRD[2][2] = {{0.456, 0.053}, {0.325, 0.022}};
  TString folder = "FitAll/fits/TextFinalData"; if(isIso==1) folder = "FitAll/fits/TextFinalIsoData"; 
  folder+=FileTag; folder+="x";
  TString nameFiles[] = {"030", "050", "075", "100", "120", "150", "200", "250", "300", "400"};

  TString RDTag[] = {"Bm", "R*0", "B0", "R*p"};
  if(isIso){RDTag[0] = "R(D)"; RDTag[1] = "R(D*)"; }

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
  double minBDT = 0, maxBDT = 330;
  TCanvas can("can","RD results");
  can.Divide(2,1); TPad *cPad;
  TString yTitle[] = {"R(D)", "R(D*)"};
  TH1F* hRD[2];
  TBox box; box.SetLineColor(0);box.SetFillColor(34);
  TLatex label; label.SetTextSize(0.08); label.SetTextFont(132); label.SetTextAlign(12);
  TMarker mark; mark.SetMarkerStyle(8); mark.SetMarkerSize(1.03);
  TLine lin; lin.SetLineStyle(1); lin.SetLineWidth(2);
  double minRD2[2] = {0.19, 0.23}, maxRD2[2] = {0.75, 0.51}; 
  if(isError) {
    minRD2[0] = 0.04; minRD2[1] = 0.01; maxRD2[0] = 0.145; maxRD2[1] = 0.075;
    if(isIso!=0){maxRD2[0] = 0.099; maxRD2[1] = 0.043;}
    yTitle[0] = "Uncertainty R(D)"; yTitle[1] = "Uncertainty R(D*)"; 
  } else if(isIso==1) {
    minRD2[0] = 0.27; minRD2[1] = 0.26; maxRD2[0] = 0.62; maxRD2[1] = 0.42;
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
    hRD[pad]->SetLabelOffset(0.01,"Y");
    hRD[pad]->GetXaxis()->CenterTitle(true);
    hRD[pad]->Draw();
    box.SetFillColor(colors[5]); 
    box.DrawBox(minBDT, nomiRD[pad][0]-nomiRD[pad][1], maxBDT, nomiRD[pad][0]+nomiRD[pad][1]);
    lin.SetLineColor(colors[4]); lin.SetLineWidth(3);
    lin.DrawLine(minBDT, nomiRD[pad][0], maxBDT, nomiRD[pad][0]);lin.SetLineWidth(2);
    for(int file=0; file<nFiles; file++){
      int nBeg = 0; if(isIso==3) nBeg=2;
      if(pad==1) nBeg++;

      int endChan = 4; if(isIso==1 || isIso==2) endChan = 2;
      for(int cand=nBeg; cand<endChan; cand+=2){
	double height = RD[file][cand][0], dH = RD[file][cand][1], BDT = nameFiles[file].Atof()-4;
	if(cand>1) BDT += 8;
	mark.SetMarkerColor(colors[cand]); lin.SetLineColor(colors[cand]);
	if(isError) height = RD[file][cand][1];
	if(file>0 || 1){
	  if(cand<2) mark.SetMarkerStyle(8);
	  else       mark.SetMarkerStyle(21);
	  lin.DrawLine(BDT, height-dH, BDT, height+dH);
	  double wid = 3;
	  lin.DrawLine(BDT-wid,  height+dH, BDT+wid,  height+dH);
	  lin.DrawLine(BDT-wid,  height-dH, BDT+wid,  height-dH);
	  mark.DrawMarker(BDT, height);
	}
	if(file==0 && cand%2==0){
	  double delta = 10; 
	  height = minRD2[pad] + (maxRD2[pad]-minRD2[pad])*(0.91-0.14*(cand>1));
	  mark.DrawMarker(25+delta, height);
	  label.SetTextSize(0.075); label.SetTextColor(colors[cand]);label.SetTextAlign(12);
	  label.DrawLatex(34+delta, height*1.0, RDTag[cand]);
	  
	}
      }
      style.setTitles(hRD[pad], "% of the nominal sample", yTitle[pad]);
    }
    hRD[pad]->Draw("axis same");
  }

  TString pName = "public_html/Stability_BDTVar"; //pName+=FileTag; 
  if(isIso) pName+="Iso"; if(isError) pName+="_Err"; 
  pName+=".eps"; 
  can.SaveAs(pName);
  for(int pad=0; pad<2; pad++){
      hRD[pad]->Delete();
  }
}

