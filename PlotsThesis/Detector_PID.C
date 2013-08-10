#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include "babar_code/Styles/Styles.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

#define nbinsMax 60

void readPIDTable(TString textName, double Tracks[3][nbinsMax], double Pranges[3][nbinsMax], int &nPbins);

void Detector_PID(int isMu=0, int isMC=0, int isRun6=1){

  TString Folder = "/afs/slac.stanford.edu/g/babar/physicstools/pid/tables/";
  TString Run = "all-r24c/";
  if(isRun6) Run = "run6-r24c/";
  TString MC = Folder; MC += "SP-r24c/asData/";
  TString PIDName[2][3] = {{"electrons+/e.emcrad.KM.Loose","electrons+/pi.ksTau31.KM.Loose", "electrons+/k.dstar.KM.Loose"},
			   {"muons+/mu.mumug2.BDT.VeryLoose", "muons+/pi.tau31.BDT.VeryLoose", "muons+/k.Dstar.BDT.VeryLoose"}};

  int nPbins[6];
  double Pranges[6][3][nbinsMax], Tracks[6][3][nbinsMax];
  for(int tab=0; tab<3; tab++){
    TString TabName = Folder; TabName += Run; TabName += PIDName[isMu][tab];
    cout<<"Doing "<< TabName<<endl;
    readPIDTable(TabName, Tracks[tab], Pranges[tab], nPbins[tab]);
    if(isMC) TabName.ReplaceAll(Folder,MC);
    else TabName.ReplaceAll("+/","-/");
    cout<<"Doing "<< TabName<<endl;
    readPIDTable(TabName, Tracks[tab+3], Pranges[tab+3], nPbins[tab+3]);
  }


  Styles style; style.setPadsStyle(3); style.applyStyle();
  TCanvas can; can.Divide(3,1);
  int colors[] = {2,4};
  TLegend *leg[3];
  TString legName[] = {"e^{+}", "#pi^{+}", "K^{+}"};
  double Hranges[3][nbinsMax];
  TH1F *hPID[3][2];
  for(int tab=0; tab<3; tab++){
    can.cd(tab+1);
    double legW = 0.17, legH = 0.19;
    double legX = 1-style.PadRightMargin-0.01, legY = 1-style.PadTopMargin-0.03;
    if(tab==0) legY -= 0.2;
    leg[tab] = new TLegend(legX-legW, legY-legH, legX, legY);
    leg[tab]->SetTextSize(style.LabelSize*1.1); leg[tab]->SetFillColor(0); leg[tab]->SetTextFont(style.nFont);
    leg[tab]->SetBorderSize(0);
    for(int p=0; p<=nPbins[tab]; p++) {
      Hranges[tab][p] = Pranges[tab][0][p];
      //cout<<Pranges[tab][0][p]<<"  ";
    }
    Hranges[tab][nPbins[tab]] = 10; 
    //cout<<endl;
    for(int cha=0; cha<2; cha++){
      TString hname = "hPID"; hname += tab; hname += cha;
      hPID[tab][cha] = new TH1F(hname,"",nPbins[tab],Hranges[tab]);
      for(int p=0; p<nPbins[tab]; p++){
	//double Error = 0, N = Tracks[tab][0][p]+Tracks[tab+3][0][p], n = Tracks[tab][1][p]+Tracks[tab+3][1][p];
	double Error = 0, N = Tracks[tab+3*cha][0][p], n = Tracks[tab+3*cha][1][p];
	if(N) {
	  Error = sqrt(n*N*N-n*n*N)/N/N;
	  hPID[tab][cha]->SetBinContent(p+1,n/N*100);
	}
	hPID[tab][cha]->SetBinError(p+1,Error*100);
      }
      style.setMarkers(hPID[tab][cha], 0.5, 20);
      hPID[tab][cha]->SetMarkerColor(colors[cha]);
      hPID[tab][cha]->SetLineColor(colors[cha]);
      hPID[tab][cha]->SetAxisRange(0,3.9);
      
      TString legLabel = legName[tab];
      if(isMu) legLabel.ReplaceAll("e","#mu");
      if(cha) legLabel.ReplaceAll("+","-");
      leg[tab]->AddEntry(hPID[tab][cha],legLabel);
      if(cha==0){
	if(tab==0) hPID[tab][cha]->SetMinimum(50);
	hPID[tab][cha]->Draw();
	style.setTitles(hPID[tab][cha],"x","Efficiency (%)");
      } else {
	hPID[tab][cha]->Draw("same");
	leg[tab]->Draw();
      }
    }
  }
  TString epsName = "public_html/Detector_PID_";
  if(isMu) epsName += "mu"; else epsName += "e";
  epsName += ".eps";
  can.SaveAs(epsName);
  for(int tab=0; tab<3; tab++) for(int cha=0; cha<2; cha++) hPID[tab][cha]->Delete(); 

}

void readPIDTable(TString textName, double Tracks[3][nbinsMax], double Pranges[3][nbinsMax], int &nPbins){
  fstream textFile; textFile.open(textName,fstream::in);

  nPbins=0;
  double pmin, pmax, phimin, phimax, thetamin, thetamax, eff, effError, iniTracks, finTracks;
  int nLines = 0; 
  for(int p=0; p<nbinsMax; p++) 
    for(int i=0; i<3; i++) {Tracks[i][p]=0; Pranges[i][p]=-1;}
  while(textFile && nLines<3000) {
    textFile>>pmin>>pmax>>thetamin>>thetamax>>phimin>> phimax>> eff>> effError>>finTracks>>iniTracks;
//     cout<<"  "<<pmin<<"  "<<pmax<<"  "<<thetamin<<"  "<<thetamax<<"  "<<phimin<<"  "<< phimax<<"  "
// 	<< eff<<"  "<< effError<<"  "<<finTracks<<"  "<<iniTracks<<endl;
    //if(iniTracks && finTracks/iniTracks>.9 && pmin>1.8 && pmin<3.3) cout<<pmin<<"  "<<pmax<<"  "<<thetamin<<"  "<<thetamax<<"  "<<finTracks/iniTracks<<endl;
    if(thetamin<22 || thetamax>138) continue;
    for(int p=0; p<=nPbins; p++){
      if(pmin==Pranges[0][p]) {Tracks[0][p]+=iniTracks; Tracks[1][p]+=finTracks; break;}
      if(p==nPbins) {
	Tracks[0][p]+=iniTracks; Tracks[1][p]+=finTracks;
	Pranges[0][p] = pmin; Pranges[1][p] = pmax; nPbins++; 
	break;}
    }
    nLines++;
  }
}


