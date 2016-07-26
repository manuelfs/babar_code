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
#include "TF1.h"
#include "babar_code/Styles/Styles.cc"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;
#define nFiles 3

void mESCorr_Bin(bool isError=false){

  Styles style; style.setPadsStyle(2); style.applyStyle();
  //gStyle->SetPadLeftMargin(0.3); gStyle->SetPadTopMargin(0.12);

  int colors[] = {2,4,28,9, 46, 38};
  double RD[nFiles][4][2], Yield[nFiles][4][3];
  double minRD[2] = {99., 99.}, maxRD[2] = {-99., -99.}; 

  TString RDTag[] = {"D^{0}", "D*^{0}", "D^{+}", "D*^{+}"};
  double nomiRD[2][2] = {{1, 0.053}, {1, 0.022}};
  //TString fileTag[] = {"50", "67", "80", "100", "120", "150", "200", "250", "300"};
  TString fileTag[] = {"100", "200", "300"};

  TString fileName = "../DonutUtils/mESCorrections.txt"; 
  fstream textFile; textFile.open(fileName,fstream::in);
  for(int file=0; file<nFiles; file++){
    TString text = ""; int times = 0;
    while(!text.Contains("sample")) {
      textFile>>text; 
      times++;
      if(times>1000) {cout<<"R(D) not found in "<<fileName<<endl; break;}
    }
    for(int cand=0; cand<4; cand++){
      textFile>>text>>text; 
      Yield[file][cand][0] = text.Atof();
      textFile>>text>>text; 
      Yield[file][cand][1] = text.Atof();
      for(int i=0; i<2; i++) textFile>>text; 
      textFile>>text; RD[file][cand][0] = text.Atof();
      textFile>>text; textFile>>text; RD[file][cand][1] = text.Atof();
      if(file>0) {
	double n1 = Yield[file][cand][0], n0 = Yield[file-1][cand][0];
	double d1 = Yield[file][cand][1], d0 = Yield[file-1][cand][1];
	double e1 = RD[file][cand][1]*d1*d1/n1, e0 = RD[file-1][cand][1]*d0*d0/n0;
	Yield[file][cand][2] = (n1-n0)/(d1-d0)/(d1-d0)*sqrt(e1*e1+e0*e0);
	RD[file][cand][0] = (n1-n0)/(d1-d0);
      } else Yield[file][cand][2] = RD[file][cand][1];
      if(RD[file][cand][0]-RD[file][cand][1] < minRD[cand%2]) minRD[cand%2] = RD[file][cand][0]-RD[file][cand][1];
      if(RD[file][cand][0]+RD[file][cand][1] > maxRD[cand%2]) maxRD[cand%2] = RD[file][cand][0]+RD[file][cand][1];
      //cout<<fileTag[file]<<"\t"<<RD[file][cand][0]<<" #pm "<<RD[file][cand][1]<<endl;
    }
    //cout<<endl;
  }
//   if(maxRD[1] > maxRD[0]) maxRD[0] = maxRD[1];
//   if(minRD[1] < minRD[0]) minRD[0] = minRD[1];

  fileName = "babar_code/Systematics/Text/mESFits.txt"; 
  fstream mESFile; mESFile.open(fileName,fstream::out);
  double shift = 100, parFits[4][2];
  TH1F h("h","",3,50-shift,350-shift);
  TF1* line = new TF1("line","[0]+[1]*x",50-shift,350-shift);
  line->SetParameter(0,0.85); line->SetParameter(1,0.0002); 
  for(int cand=0; cand<4; cand++){
    for(int file=0; file<nFiles; file++){
      h.SetBinContent(file+1,RD[file][cand][0]); h.SetBinError(file+1,RD[file][cand][1]); 
    }
    h.Fit("line");
    parFits[cand][0] = line->GetParameter(0); parFits[cand][1] = line->GetParameter(1); 
    mESFile<<RDTag[cand]<<"\t"<<RoundNumber(line->GetParameter(0),4)<<"\t"<<RoundNumber(line->GetParameter(1)*10000,4)<<endl;
  }
  mESFile.close();
  cout<<"Saved fits in "<<fileName<<endl;
  double minBDT = 0, maxBDT = 350;
  TCanvas can("can","RD results");
  can.Divide(2,1); TPad *cPad;
  TString yTitle[] = {"R(D)", "R(D*)"};
  TH1F* hRD[2];
  TBox box; box.SetFillStyle(3002); box.SetLineColor(0);box.SetFillColor(34);
  TLatex label; label.SetTextSize(0.08); label.SetTextFont(132); label.SetTextAlign(12);
  TMarker mark; mark.SetMarkerStyle(8); mark.SetMarkerSize(1.05);
  TLine lin; lin.SetLineStyle(1); lin.SetLineWidth(2);
  double minRD2[2] = {0.75, 0.75}, maxRD2[2] = {1.1, 1.1}; 
  if(isError) {
    minRD2[0] = 0; minRD2[1] = 0; maxRD2[0] = 41; maxRD2[1] = 25;
    yTitle[0] = "Error (%)"; yTitle[1] = "Error (%)"; 
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
    //if(!isError) {box.SetFillColor(colors[4+pad]); box.DrawBox(minBDT, nomiRD[pad][0]-nomiRD[pad][1], maxBDT, nomiRD[pad][0]+nomiRD[pad][1]);}
    if(!isError) {lin.SetLineColor(colors[4+pad]); lin.DrawLine(minBDT, nomiRD[pad][0], maxBDT, nomiRD[pad][0]);}
    for(int file=0; file<nFiles; file++){
      if(fileTag[file]=="67" || fileTag[file]=="80" || fileTag[file]=="120") continue;
      if(fileTag[file].Contains("50")) continue;
      int nBeg = 0; if(pad==1) nBeg = 1;
      for(int cand=nBeg; cand<4; cand+=2){
	double height = RD[file][cand][0], dH = Yield[file][cand][2], BDT = fileTag[file].Atof()-4;
	if(cand>1) BDT += 8;
	//double height = 1 - (halfH*(2*file+1) - halfH/3.2);
	//if(cand>1) height = 1 - (halfH*(2*file+1) + halfH/3.2);
	mark.SetMarkerColor(colors[cand]); lin.SetLineColor(colors[cand]);
	if(isError) height = RD[file][cand][1]/RD[file][cand][0]*100;
	mark.DrawMarker(BDT, height);
	double maxHeight = height+dH; if(maxHeight>maxRD2[pad]) maxHeight=maxRD2[pad];
	if(!isError) lin.DrawLine(BDT, height-dH, BDT, maxHeight);
	if(file==1){
	  double delta = 250; if(isError) delta = 100;
	  height = minRD2[pad] + (maxRD2[pad]-minRD2[pad])*(0.93-0.1*(cand>1));
	  mark.DrawMarker(25+delta, height);
	  label.SetTextSize(0.075); label.SetTextColor(colors[cand]);label.SetTextAlign(12);
	  label.DrawLatex(34+delta, height*1.005, RDTag[cand]);
	  
	}
      }
      style.setTitles(hRD[pad], "% of the nominal sample", "Sideband correction");
    }
  }
  for(int cand=0; cand<4; cand++){
    can.cd((cand%2)+1);
    lin.SetLineColor(colors[cand]);
    lin.DrawLine(0, parFits[cand][0]+parFits[cand][1]*(0-shift), 350, parFits[cand][0]+parFits[cand][1]*(350-shift));
  }

  TString pName = "public_html/mESCorr_Bin.eps"; 
  can.SaveAs(pName);
  for(int pad=0; pad<2; pad++){
      hRD[pad]->Delete();
  }
}

