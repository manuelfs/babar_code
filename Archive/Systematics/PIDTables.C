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

#define nbinsMax 50

void textPIDTable(int Tracks[4][2][nbinsMax], double Pranges[4][2][nbinsMax], int nPbins[4], int iTable);
void plotPIDTable(int Tracks[4][2][nbinsMax], double Pranges[4][2][nbinsMax], int nPbins[4]);
void readPITTable(TString textName, int Tracks[2][nbinsMax], double Pranges[2][nbinsMax], int &nPbins);
void printPIDTable(int Tracks[2][nbinsMax], double Pranges[2][nbinsMax], int nPbins);

void PIDTables(int iTable=-1){

  if(iTable<0 || iTable>3){
    cout<<"iTable: 0 Standard e, 1 Standard mu, 2 Beamspot e, 3 Beamspot mu"<<endl;
    return;
  }

//   TString PIDName[4][2] = {{"/afs/slac.stanford.edu/g/babar/physicstools/pid/tables/SP-r24c/asData/run6-r24c/electrons+/e.emcrad.KM.Loose",
// 			    "/afs/slac.stanford.edu/g/babar/physicstools/pid/tables/run6-r24c/electrons+/e.emcrad.KM.Loose"},
// 			   {"/afs/slac.stanford.edu/g/babar/physicstools/pid/tables/SP-r24c/asData/run6-r24c/muons+/mu.mumug2.BDT.VeryLoose",
// 			    "/afs/slac.stanford.edu/g/babar/physicstools/pid/tables/run6-r24c/muons+/mu.mumug2.BDT.VeryLoose"},
  TString PIDName[4][2] = {{"/afs/slac.stanford.edu/g/babar/physicstools/pid/tables/SP-r24c/asData/all-r24c/electrons+/e.emcrad.KM.Loose",
			    "/afs/slac.stanford.edu/g/babar/physicstools/pid/tables/all-r24c/electrons+/e.emcrad.KM.Loose"},
			   {"/afs/slac.stanford.edu/g/babar/physicstools/pid/tables/SP-r24c/asData/all-r24c/muons+/mu.mumug2.BDT.VeryLoose",
			    "/afs/slac.stanford.edu/g/babar/physicstools/pid/tables/all-r24c/muons+/mu.mumug2.BDT.VeryLoose"},
			   {"/nfs/farm/babar/AWG100/PID/BeamSpotTest/tables/SP-r25/asData/modified/electrons+/e.emcrad.KM.Loose",
			    "/nfs/farm/babar/AWG100/PID/BeamSpotTest/tables/SP-r25/asData/standard/electrons+/e.emcrad.KM.Loose"},
			   {"/nfs/farm/babar/AWG100/PID/BeamSpotTest/tables/SP-r25/asData/modified/muons+/mu.mumug2.BDT.VeryLoose",
			    "/nfs/farm/babar/AWG100/PID/BeamSpotTest/tables/SP-r25/asData/standard/muons+/mu.mumug2.BDT.VeryLoose"}};

  int nPbins[4], Tracks[4][2][nbinsMax];
  double Pranges[4][2][nbinsMax];
  for(int tab=0; tab<2; tab++){
    TString TabName = PIDName[iTable][tab];
    cout<<endl<<"Doing "<< TabName<<endl<<endl;
    readPITTable(TabName, Tracks[tab], Pranges[tab], nPbins[tab]);
    //printPIDTable(Tracks[tab], Pranges[tab], nPbins[tab]);
    TabName.ReplaceAll("+/","-/");
    readPITTable(TabName, Tracks[tab+2], Pranges[tab+2], nPbins[tab+2]);
  }
  for(int tab=0; tab<2; tab++)
    for(int p=0; p<nPbins[0]; p++) 
      for(int i=0; i<2; i++) 
	Tracks[tab][i][p] =Tracks[tab+2][i][p];

  textPIDTable(Tracks, Pranges, nPbins, iTable);

  plotPIDTable(Tracks, Pranges, nPbins);

}

void textPIDTable(int Tracks[4][2][nbinsMax], double Pranges[4][2][nbinsMax], int nPbins[4], int iTable){

  TString outName = "babar_code/Systematics/Text/RatioPID"; outName += iTable; outName += ".txt";
  fstream outFile; outFile.open(outName,fstream::out);

  for(int p=0; p<nPbins[0]; p++) {
    double Ratio = 0, iniM = Tracks[0][0][p], finM = Tracks[0][1][p], iniS = Tracks[1][0][p], finS = Tracks[1][1][p];
    if(finM*iniS) Ratio = iniM*finS/finM/iniS;
    outFile<<RoundNumber(Pranges[0][0][p],1)<<" - "<<RoundNumber(Pranges[0][1][p],1)
	   <<"\t"<<RoundNumber(Ratio,5)<<endl;
  }
  cout<<"Written "<<outName<<endl;
}

void readPITTable(TString textName, int Tracks[2][nbinsMax], double Pranges[2][nbinsMax], int &nPbins){
  fstream textFile; textFile.open(textName,fstream::in);

  nPbins=0;
  double pmin, pmax, phimin, phimax, thetamin, thetamax, eff, effError;
  int nLines = 0, iniTracks, finTracks; 
  for(int p=0; p<nbinsMax; p++) 
    for(int i=0; i<2; i++) {Tracks[i][p]=0; Pranges[i][p]=-1;}
  while(textFile && nLines<3000) {
    textFile>>pmin>>pmax>>thetamin>>thetamax>>phimin>> phimax>> eff>> effError>>finTracks>>iniTracks;
//       cout<<"  "<<pmin<<"  "<<pmax<<"  "<<thetamin<<"  "<<thetamax<<"  "<<phimin<<"  "<< phimax<<"  "
//       	  << eff<<"  "<< effError<<"  "<<finTracks<<"  "<<iniTracks<<endl;
    if(thetamin<10 || thetamin>145) continue;
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

void plotPIDTable(int Tracks[4][2][nbinsMax], double Pranges[4][2][nbinsMax], int nPbins[4]){

  Styles style; style.setPadsStyle(1); style.applyStyle();
  TCanvas can;
  int colors[] = {2,4};
  double Hranges[nbinsMax];
  for(int p=0; p<nPbins[0]; p++) Hranges[p] = Pranges[0][0][p];
  Hranges[nPbins[0]] = 6.2;
  TH1F *hPID[2];
  for(int tab=0; tab<2; tab++){
    TString hname = "hPID"; hname += tab;
    hPID[tab] = new TH1F(hname,"",nPbins[0],Hranges);
    for(int p=0; p<nPbins[0]; p++){
      double Error = 0, N = Tracks[tab][0][p]+Tracks[tab+2][0][p], n = Tracks[tab][1][p]+Tracks[tab+2][1][p];
      if(N) {
	Error = sqrt(n*N*N-n*n*N)/N/N;
	hPID[tab]->SetBinContent(p+1,n/N);
      }
      hPID[tab]->SetBinError(p+1,Error);
      hPID[tab]->SetMarkerColor(colors[tab]);
      hPID[tab]->SetLineColor(colors[tab]);
    }
    //hPID[tab]->SetMaximum(0.86);
    hPID[tab]->SetMinimum(0.7);
  }
  hPID[0]->Draw();
  hPID[1]->Draw("same");
  style.setTitles(hPID[0],"p_{l} (GeV)","Efficiency");
  can.SaveAs("public_html/Jul6_PIDBeamspot.eps");
  hPID[0]->Delete();  hPID[1]->Delete();
}

void printPIDTable(int Tracks[2][nbinsMax], double Pranges[2][nbinsMax], int nPbins){
  double Perror[nbinsMax];
  for(int p=0; p<nPbins; p++) {
    double n = Tracks[1][p], N = Tracks[0][p];
    if(N) Perror[p] = sqrt(n*N*N-n*n*N)/N/N;
    else Perror[p] = 999;
    cout<<RoundNumber(Pranges[0][p],1)<<"-"<<RoundNumber(Pranges[1][p],1)
	<<":\t "<<RoundNumber(n,3,N)<<" +- "
	<<RoundNumber(Perror[p],3)<<", a "<<RoundNumber(Perror[p]*N*100,1,n)<<" % error"
	<<"\t\t"<<n<<" of "<<N<<" tracks reconstructed"<<endl;
  }
}
