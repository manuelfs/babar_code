//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: FitPdfSampleKEYS.cc,v 1.7 2012/03/04 00:33:02 manuelf Exp $
//
// Description:
//      FitPdfSampleKEYS - Fits one sample using the KEYS function
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      10/02/12 manuelf -- Plotting code moved to KeysUtils.cc
//      10/02/10 manuelf -- Adapted to the new numbering of the samples
//      09/12/08 manuelf -- Created
//------------------------------------------------------------------------

#include "TString.h"
#include "TFile.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooAddPdf.hh"
#include "RooFitCore/RooAbsReal.hh"
#include "RooFitCore/RooDataSet.hh"
#include "DonutUtils/RooNDKeysPdfDonut.hh"
#include "DonutUtils/KeysUtils.cc"
#include <fstream>
#include <iostream>
#include <ctime>

using namespace std;
using namespace RooFit;
using std::cout;
using std::endl;

void Stitch(TH2F *hResult, TH2F *hTail, double *stitches, int nM2bin, int nPlbin, double entries[], int stitchP=1);
void SaveHisto(TH2F *h2, TString hname, int Sam=0);

int main(int argc, char *argv[]){

  time_t start,end;
  time (&start);
  double dif;

  if (argc < 3 || argc > 18 ) {
    cout << "USAGE: FitPdfSampleKEYS sample smoothing [scaleP=200] [nM2bin=1000] [nPlbin=300] "
	 <<"[smooth2=50] [scaleP2=200] [stitch1=0.35] [smooth3=50] [scaleP3=200] [stitch2=1.5] [doRot=Rot] [dataFolder=RAll] "
	 << "[smooth4=50] [scaleP4=200] [stitch4=-4] [Times=1]" << endl;
    return 0;
  }
  TString sam = argv[1]; int Sam = sam.Atoi();
  TString smoo_s = argv[2]; double smoo = smoo_s.Atof()/100.;
  TString scaleP_s = "200";
  if(argc>3) scaleP_s = argv[3]; 
  double scaleP = scaleP_s.Atof()/100.;
  int nM2bin = 1000;
  if(argc>4)  {TString sM2bin = argv[4]; nM2bin = sM2bin.Atoi();}
  int nPlbin = 300;
  if(argc>5)  {TString sPlbin = argv[5]; nPlbin = sPlbin.Atoi();}
  TString smoo2_s = "50";
  if(argc>6)  smoo2_s = argv[6]; 
  double smoo2 = smoo2_s.Atof()/100.;
  TString scaleP2_s = "200";
  if(argc>7) scaleP2_s = argv[7]; 
  double scaleP2 = scaleP2_s.Atof()/100.;
  double stitches[] = {0.35, 1.5, -4., -4.};
  if(argc>8)  {TString sTemp = argv[8]; stitches[0] = sTemp.Atof();}
  TString smoo3_s = "50";
  if(argc>9)  smoo3_s = argv[9]; 
  double smoo3 = smoo3_s.Atof()/100.;
  TString scaleP3_s = "200";
  if(argc>10) scaleP3_s = argv[10]; 
  double scaleP3 = scaleP3_s.Atof()/100.;
  if(argc>11)  {TString sTemp = argv[11]; stitches[1] = sTemp.Atof();}
  bool doRot = true;
  if(argc>12)  {TString t_s = argv[12]; if(t_s == "noRot") doRot = false;}
  if(doRot) cout<<"Rotating data"<<endl;
  TString dataFolder = "RAll";
  if(argc>13) dataFolder = argv[13]; 
  TString smoo4_s = "50";
  if(argc>14)  smoo4_s = argv[14]; 
  double smoo4 = smoo4_s.Atof()/100.;
  TString scaleP4_s = "200";
  if(argc>15) scaleP4_s = argv[15]; 
  double scaleP4 = scaleP4_s.Atof()/100.;
  if(argc>16)  {TString sTemp = argv[16]; stitches[2] = sTemp.Atof();}
  int Times = 1, IsCocktail = 0;
  if(argc>17)  {TString sTimes = argv[17]; Times = sTimes.Atoi();}
  int isSyst = 0; if(Times<0){isSyst = 1; Times *= -1;}
  if(dataFolder.Contains("cocktail")) IsCocktail = 1;

  bool doStitch = false; if(stitches[0]>0 || stitches[1]>0) doStitch = true;
  if(stitches[0]>0){
    int binCut = (int)(stitches[0]/2.4*(double)nPlbin)+1;
    stitches[0] = (binCut-1)*2.4/nPlbin;
  }
  if(stitches[1]>0){
    int binCut = (int)(stitches[1]/2.4*(double)nPlbin)+1;
    stitches[1] = (binCut-1)*2.4/nPlbin;
  }
  if(stitches[2]>-4){
    int binCut = (int)((stitches[2]+4)/16.*(double)nM2bin)+1;
    stitches[2] = (binCut-1)*16./nM2bin-4;
    //binCut = (int)((stitches[2]+4.4)/16.*(double)nM2bin)+1;
    //if(stitches[2]<1.8 && Sam<21) binCut = (int)(5.8/16.*(double)nM2bin)+1;  // If the stitching is for mmiss < 1.8 GeV^2,
    stitches[3] = (binCut-1)*16./nM2bin-4;                         // I still normalize the tail with only mmiss > 1.8
  }
  TString adaptive = "a";
  TString Base =  "_"; Base += sam; Base += "_"; 
  Base += nM2bin*1000+nPlbin; Base += "_"; Base += smoo_s; 
  if(doStitch) {Base+="and"; Base+=smoo2_s;}
  Base += "_"; Base += scaleP_s; if(doStitch) {Base+="and"; Base+=scaleP2_s;}

  double entriesCut[]={0,0,0,0};
  float candM2,candPstarLep, weight=-1;
  TString inputfile = nameData(Sam,dataFolder);
  TChain *tree = new TChain("ntp1");
  if(isSyst==0){
    tree->Add(inputfile);
    tree->SetBranchAddress("candM2",&candM2);
    tree->SetBranchAddress("candPstarLep",&candPstarLep);
    tree->SetBranchAddress("weight",&weight);
  }

  double totMCB = 0, totuds = 0, totccbar = 0, totdata = 0, totOffdata = 0;
  getNumberB("AWG82/ntuples/small/RAll_RunAll.root", "All", totMCB, totdata, totuds, totccbar, totOffdata);
  double wMC = totdata/totMCB; 
  TRandom3 rand; 
  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4;
  double entries = tree->GetEntries();
  TH2F *h2 = new TH2F("h2","KEYS",nM2bin,m2min,m2max,nPlbin,plmin,plmax);
  TString hnameFit = "keys/root/fitFinal/pdfKeys_";
  hnameFit += sam; hnameFit += "_Fit.root";
  TString hname = "keys/root/";
  if(Times==1) {
    hname += "fitPdf/hKeys"; hname += Base; hname += ".root";
  } else {
    if(isSyst){
      hname = "AWG82/systematics/BF/KeysTimes/hTimes_"; hname += sam; hname += "_Fit.root";
    } else {
      hname += "Times/hTimes_"; hname += sam; hname += "_Fit.root";
    }
  }
  TFile* hfile = new TFile(hname,"RECREATE"); 
  hfile->Close(); hfile->Delete();
  for(int rep = 0; rep < Times; rep++){
    if(isSyst) {
      tree->Delete();
      tree = new TChain("ntp1");
      inputfile = "AWG82/systematics/BF/fitSamples/Iter"; inputfile += rep; 
      inputfile += "/pdfSample"; inputfile += Sam; inputfile += ".root";
      tree->Add(inputfile);
      tree->SetBranchAddress("candM2",&candM2);
      tree->SetBranchAddress("candPstarLep",&candPstarLep);
      tree->SetBranchAddress("weight",&weight);
      entries = tree->GetEntries();
    }
    for(int iCut=0; iCut<4; iCut++) entriesCut[iCut] = 0;
    RooRealVar mmiss2("candM2","candM2",m2min,m2max);
    RooRealVar pstarl("candPstarLep","candPstarLep",plmin,plmax);
    RooRealVar totWeight("totWeight","totWeight",0.,100.);
    RooArgSet myVars(mmiss2,pstarl,totWeight);
    RooDataSet  data("data","data",myVars);
    if(rep%50==0) cout<<"Repetition "<<rep+1<<" of "<< Times <<" - Smoothings: "<<smoo_s<<"/"<<scaleP_s<<endl; 
    for (int evt = 0 ; evt < entries; evt ++) {
      if(Times>1 && isSyst==0) {
	tree->GetEvent((int)rand.Uniform(entries));
      } else tree->GetEvent(evt);
      if(isSyst) weight /= wMC;
      if(candM2<m2min||candM2>m2max || candPstarLep<plmin||candPstarLep>plmax || weight<0||weight>100) continue;
      mmiss2.setVal(candM2);
      pstarl.setVal(candPstarLep);
      totWeight.setVal(weight);
      data.add(RooArgSet(mmiss2,pstarl,totWeight));
      if(candPstarLep<stitches[0]) entriesCut[0] += weight;
      if(candPstarLep>stitches[1]) entriesCut[1] += weight;
      if(candM2>stitches[3]) entriesCut[2] += weight;
      entriesCut[3] += weight;
    }
    data.setWeightVar(totWeight);
    RooNDKeysPdfDonut KEYS("KEYS","KEYS",RooArgList(mmiss2,pstarl),data,adaptive,smoo,smoo*scaleP,3.,doRot);
    time (&end);dif = difftime (end,start);
    if(rep%50==0) cout<<dif<<" sec after first KEYS function with smoothings "<<smoo_s<<"/"<<scaleP_s<<endl;
    time (&start);
    TString hNameTimes = "hTimes"; hNameTimes += rep;
    if(Times>1) h2 = new TH2F(hNameTimes,"KEYS",nM2bin,m2min,m2max,nPlbin,plmin,plmax);
    KEYS.fillHistogram(h2,RooArgList(mmiss2,pstarl));
    time (&end);dif = difftime (end,start);
    if(rep%50==0) cout<<dif<<" seconds after storing the first KEYS function"<<endl;
    time (&start);
    if(doStitch) {
      if(stitches[2]>-4) {
	if(rep%50==0) cout<<"Stitching smoothings "<<smoo4_s<<"/"<<scaleP4_s<<" at mmiss "<<stitches[2]<<endl;
	RooNDKeysPdfDonut KEYS3("KEYS3","KEYS",RooArgList(mmiss2,pstarl),data,adaptive,smoo4,smoo4*scaleP4,3.,doRot);
	TH2F *hsti = new TH2F("hsti3","Histo to be stitched",nM2bin,m2min,m2max,nPlbin,plmin,plmax);
	KEYS3.fillHistogram(hsti,RooArgList(mmiss2,pstarl));
	Stitch(h2, hsti, stitches, nM2bin, nPlbin,entriesCut, 4);
	hsti->Delete();
      }
      if(rep%50==0) cout<<"Stitching smoothings "<<smoo2_s<<"/"<<scaleP2_s<<" at p*l "<<stitches[0]<<" and "
			<<smoo3_s<<"/"<<scaleP3_s<<" at p*l "<<stitches[1]<<endl;
      if(stitches[0]>0) {
	RooNDKeysPdfDonut KEYS2("KEYS2","KEYS",RooArgList(mmiss2,pstarl),data,adaptive,smoo2,smoo2*scaleP2,3.,doRot);
	TH2F *hsti = new TH2F("hsti","Histo to be stitched",nM2bin,m2min,m2max,nPlbin,plmin,plmax);
	KEYS2.fillHistogram(hsti,RooArgList(mmiss2,pstarl));
	Stitch(h2, hsti, stitches, nM2bin, nPlbin,entriesCut,2);
	hsti->Delete();
	time (&end);dif = difftime (end,start);
	if(rep%50==0) cout<<dif<<" seconds after the first Pl stitches "<<endl;
	time (&start);
      }
      if(stitches[1]>0) {
	RooNDKeysPdfDonut KEYS2("KEYS2","KEYS",RooArgList(mmiss2,pstarl),data,adaptive,smoo3,smoo3*scaleP3,3.,doRot);
	TH2F *hsti = new TH2F("hsti","Histo to be stitched",nM2bin,m2min,m2max,nPlbin,plmin,plmax);
	KEYS2.fillHistogram(hsti,RooArgList(mmiss2,pstarl));
	Stitch(h2, hsti, stitches, nM2bin, nPlbin,entriesCut,3);
	hsti->Delete();
	time (&end);dif = difftime (end,start);
	if(rep%50==0) cout<<dif<<" seconds after the second Pl stitches "<<endl;
	time (&start);
      }
    }
    hfile = new TFile(hname,"UPDATE"); 
    h2->Write();
    hfile->Close(); hfile->Delete();
    if(Times==1) SaveHisto(h2, hnameFit,Sam);
    time (&end);dif = difftime (end,start);
    if(rep%50==0) cout<<dif<<" seconds after making histogram"<<endl;
    time (&start);
    if(Times>1) h2->Delete();
  }
  cout<<hname<<" saved"<<endl; 
  if(Times>1) {return Times;}

  PlotEps(hname, dataFolder, "fitPdf", "", "EPS");

  return 1; 

}

void Stitch(TH2F *hResult, TH2F *hTail, double *stitches, int nM2bin, int nPlbin, double entries[], int stitchP){

  int binRange = nPlbin/20;
  if(stitches[0]>0 && stitchP==2){
    int binCut = (int)(stitches[0]/2.4*(double)nPlbin+0.5)+1;
    if(binRange+1>binCut) binRange = binCut-1;
    double factor = hResult->Integral(1,nM2bin,binCut+1,nPlbin)/hTail->Integral(1,nM2bin,1,binCut)*entries[0]/(entries[3]-entries[0]);
    for(int i=1; i<nM2bin+1; i++){
      for(int j=1; j<binCut+binRange+1; j++){
	double val=0;
	if(j<=binCut-binRange) {
	  val = hTail->GetBinContent(i,j)*factor;
	}else  {
	  double weight = (j-binCut+binRange)/(double)(2*binRange);
	  val = weight*hResult->GetBinContent(i,j) + (1.-weight)*hTail->GetBinContent(i,j)*factor;
	}
	hResult->SetBinContent(i,j,val);
      }
    }
  }
  if(stitches[1]>0 && stitchP==3){
    int binCut = (int)(stitches[1]/2.4*(double)nPlbin+0.5)+1;
    if(binCut+binRange>nPlbin) binRange = nPlbin-binCut-1;
    double factor = hResult->Integral(1,nM2bin,1,binCut-1)/hTail->Integral(1,nM2bin,binCut,nPlbin)*entries[1]/(entries[3]-entries[1]);
//     cout<<"binCut "<<binCut<<" and binRange "<<binRange<<". factor "<<factor<<endl;
//     cout<<"resInt "<<hResult->Integral(1,nM2bin,1,binCut-1)<<", tailInt "<<hTail->Integral(1,nM2bin,binCut,nPlbin)
// 	<<" - ratio "<<hTail->Integral(1,nM2bin,binCut,nPlbin)/hResult->Integral(1,nM2bin,1,binCut-1)<<endl;
//     cout<<"entries[3]-entries[1] "<<entries[3]-entries[1]<<", entries[1] "<<entries[1]<<" - ratio "<<
//    entries[1]/(entries[3]-entries[1])<<endl;
    for(int i=1; i<nM2bin+1; i++){
      for(int j=binCut-binRange+1; j<nPlbin+1; j++){
	double val=0;
	if(j>=binCut+binRange) {
	  val = hTail->GetBinContent(i,j)*factor;
	}else  {
	  double weight = (j-binCut+binRange)/(double)(2*binRange);
	  val = (1.-weight)*hResult->GetBinContent(i,j) + weight*hTail->GetBinContent(i,j)*factor;
	}
	hResult->SetBinContent(i,j,val);
      }
    }
  }
  binRange = nM2bin/40;
  if(stitchP==4){
    int binCut = (int)((stitches[2]+4)/16.*(double)nM2bin+0.5)+1;
    int binCut3 = (int)((stitches[3]+4)/16.*(double)nM2bin+0.5)+1;
    if(binCut+binRange>nM2bin) binRange = nM2bin-binCut-1;
    double factor = hResult->Integral(1,binCut3-1,1,nPlbin)/hTail->Integral(binCut3,nM2bin,1,nPlbin)*entries[2]/(entries[3]-entries[2]);
    //factor *= 1.2;
//     cout<<"binCut "<<binCut3<<" and binRange "<<binRange<<". factor "<<factor<<endl;
//     cout<<"resInt "<<hResult->Integral(1,binCut3-1,1,nPlbin)<<", tailInt "<<hTail->Integral(binCut3,nM2bin,1,nPlbin)
// 	<<". entries[2] "<<entries[2]<<". entries[3] "<<entries[3]<<endl;
    for(int i=1; i<nPlbin+1; i++){
      for(int j=binCut-binRange+1; j<nM2bin+1; j++){
	double val=0;
	if(j>=binCut+binRange) {
	  val = hTail->GetBinContent(j,i)*factor;
	}else  {
	  double weight = (j-binCut+binRange)/(double)(2*binRange);
	  val = (1.-weight)*hResult->GetBinContent(j,i) + weight*hTail->GetBinContent(j,i)*factor;
	}
	hResult->SetBinContent(j,i,val);
      }
    }
  }
  return;
}

// For combinatoric and crossfeed the same fit is used several times, 
// so several copies are saved
void SaveHisto(TH2F *h2, TString hname, int Sam){
  TFile* hfileFit = new TFile(hname,"RECREATE"); 
  hfileFit->cd();
  h2->Write();
  hfileFit->Close(); hfileFit->Delete();
  if(Sam==700){ // I might use it for pdf 69 and 70
    hname.ReplaceAll("48","55");
    SaveHisto(h2, hname);
    hname.ReplaceAll("55","57");
    SaveHisto(h2, hname);
  }
}

