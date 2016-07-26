//----------------------------------------------------------------------------
// keys_utils: utilities for fitting with the keys functions
//
// Description:
//      UTILITIES FOR FITTING WITH THE KEYS
//      ParseName    - Parses the name of a histogram, finding sample number
//                     number of bins and smoothing
//      nameData     - Returns the filename of a pdf sample
//      FormatHisto  - Format for histograms with data points
//      setBins      - Returns plotting range and number of bins
//      WeightedTree - Returns a tree with a branch "weight"
//      RoundNumber  - Converts a number to a string with certain decimals
//      binHisto     - Projects a 2D histogram with certain binning
//      PlotFinalFit - Reads the results of a final fit a plots them
//      PlotEps      - Produces a ps/eps plot from a histogram file
//      IntG         - Integral of a gaussian in a range
//      GaussConv    - Convolution of a histogram and a gaussian
//      getNumberB   - Number of B (and uds,cc) in the files
//      Choleski     - Choleski decomposition
//      ReadFitFile  - Reads a text file with the results of a fit
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      12/04/25 manuelf -- Added ReadFitFile
//      12/03/24 manuelf -- Added Choleski
//      10/10/28 manuelf -- Added getNumberB
//      10/10/15 manuelf -- Added IntG and GaussConv
//      10/03/10 manuelf -- Added PlotFinalFit
//      10/03/02 manuelf -- Added binHisto
//      10/02/12 manuelf -- Added plotting code, PlotEps
//      10/01/18 manuelf -- Created
//----------------------------------------------------------------------------


#include <cmath>
#include <fstream>
#include <iostream>

#include "TString.h"
#include "TFile.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TMath.h"
#include "TROOT.h"
#include "TMatrixT.h"

#include "styles.hpp"
#include "keys_utils.hpp"

using namespace TMath;
using namespace std;

double IntG(double mean, double sigma, double minX, double maxX){
  return (erf((maxX-mean)/sigma/sqrt(2.))-erf((minX-mean)/sigma/sqrt(2.)))/2.;
}

double getWentries(TTree *tree, TString Wname){
  double Wentries = 0; float weight=0;
  if(tree==0) return 0;
  tree->SetBranchAddress(Wname,&weight);
  for(int entry=0; entry<tree->GetEntries(); entry++){
    tree->GetEvent(entry);
    Wentries += weight;
  }
  return Wentries;
}

void ParseName(TString hname, TString &sam, TString &smoo, TString &scaleP, TString &Base, 
	       int &Sam, int &nM2bin, int &nPlbin){
  sam = hname; sam.Remove(0,sam.Last('/')+1); sam.Remove(0,sam.First('_')+1);
  sam.Remove(sam.First('_'),sam.Length());
  Sam = sam.Atoi();

  scaleP = hname; scaleP.Remove(0,scaleP.Last('_')+1);
  scaleP.Remove(scaleP.First('.'),scaleP.Length());

  smoo = hname;smoo.Remove(0,smoo.First('_')+1);smoo.Remove(0,smoo.First('_')+1);
  smoo.Remove(0,smoo.First('_')+1); smoo.Remove(smoo.First('.'),smoo.Length());

  TString binString = hname; binString.Remove(binString.Last('_'),binString.Length());
  binString.Remove(binString.Last('_'),binString.Length());
  binString.Remove(0,binString.Last('_')+1);
  nM2bin = binString.Atoi(); nPlbin = nM2bin%1000;
  nM2bin = nM2bin/1000;

  Base = hname; Int_t slashpos = Base.Last('s');
  if (slashpos>=0) {
    Base.Remove(0,slashpos+1); Base.Remove(Base.First('.'),Base.Length());
  } else { cout<<hname<<" is incorrect"<<endl; return;}
  return;
}

TString nameData(int Sam, TString folder){
  TString inputfile = "fitSamples/"; 
  if(folder != "") {inputfile += folder; inputfile += "/";}
  inputfile += "pdfSample";inputfile += Sam; inputfile += ".root";

  return inputfile;
}

// Format for histograms with data points
void formatHisto(TH1F *h){
  h->SetMinimum(0);
  h->SetTitle("");
  h->SetMarkerStyle(20);
  h->SetMarkerSize(.6);
  h->GetYaxis()->SetLabelSize(0.062);
  h->GetXaxis()->SetLabelSize(0.055);
  h->GetYaxis()->SetNdivisions(6+100*2);
  h->GetXaxis()->SetNdivisions(7+100*2);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(1.);
  return;
}

// Plotting range and number of bins
void setBins(int Sam, double &xlow, double &xhigh, int &nbinx, int &nbiny, int IsCocktail){
  if(IsCocktail!=1){cout<<"IsCocktail is "<<IsCocktail<<endl;} // Needed to use it to compile
  if (Sam>=1 && Sam<=8) {
    xlow = -1; xhigh = 10; 
    if(Sam<6) {nbinx = 40; nbiny = 40;}
    else {nbinx = 30; nbiny = 30;}
  }
  if (Sam>=9 && Sam<=12) {
    xlow = -0.4; xhigh = 0.5; 
    nbinx = 50; nbiny = 50;
  }
  if(Sam==13 || Sam==15){
    xlow = -0.6; xhigh = 3;
    nbinx = 50; nbiny = 50;
  }
  if(Sam==14 || Sam==16){
    xlow = -2; xhigh = 1.5;
    nbinx = 35; nbiny = 35;
  }
  if(Sam>=17 && Sam <=20){
    xlow = -1; xhigh = 10;
    nbinx = 40; nbiny = 40;
  }
  if(Sam>=21 && Sam <=24){
    xlow = -1.; xhigh = 6;
    nbinx = 40; nbiny = 40;
  }
  if (Sam>=25 && Sam<=28) {
    xhigh = 7; 
    nbinx = 40; nbiny = 40;
  }
  if (Sam>=29 && Sam<=32) {
    xhigh = 7; 
    nbinx = 35; nbiny = 35;
  }
  if (Sam>=33 && Sam<=36) {
    xlow = -1.5; xhigh = 7.5; 
    nbinx = 30; nbiny = 30;
  }
  if (Sam>=37 && Sam<=40) {
    xlow = -1.5; xhigh = 8.5; 
    nbinx = 25; nbiny = 25;
  }
  if (Sam >= 41 && Sam<=44) {
    xlow = -1.5; xhigh = 10;
    nbinx = 40; nbiny = 40;
  }
  if (Sam >= 45 && Sam<=50) {
    xlow = -3; xhigh = 10;
    nbinx = 35; nbiny = 35;
  }
  if (Sam >= 51 && Sam<=60) {
    xlow = -1.5; xhigh = 10;
    nbinx = 40; nbiny = 40;
  }
  if (Sam >= 61 && Sam<=70) {
    xlow = -1.5; xhigh = 10;
    nbinx = 35; nbiny = 35;
  }
  if(Sam==8 || Sam==16 || Sam==32 || Sam==64) { 
    nbinx = 15; nbiny = 15;
  }
  return;
}


bool isSample(int sam, int MCType, int candType){
  if (sam==1  && (candType==1 && MCType==5 )) return true;
  if (sam==2  && (candType==2 && MCType==6 )) return true;
  if (sam==3  && (candType==3 && MCType==11)) return true;
  if (sam==4  && (candType==4 && MCType==12)) return true;
  if (sam==5  && (candType==1 && MCType==6)) return true;
  if (sam==6  && (candType==2 && MCType==5)) return true;
  if (sam==7  && (candType==3 && MCType==12)) return true;
  if (sam==8  && (candType==4 && MCType==11)) return true;
  if (sam==9  && (candType==1 && (MCType==1||MCType==3))) return true;
  if (sam==10 && (candType==2 && (MCType==2||MCType==4))) return true;
  if (sam==11 && (candType==3 && (MCType==7||MCType==9))) return true;
  if (sam==12 && (candType==4 && (MCType==8||MCType==10))) return true;
  if (sam==13 && (candType==1 && (MCType==2||MCType==4) )) return true;
  if (sam==14 && (candType==2 && (MCType==1||MCType==3))) return true;
  if (sam==15 && (candType==3 && (MCType==8||MCType==10))) return true;
  if (sam==16 && (candType==4 && (MCType==7||MCType==9))) return true;
  if (sam==17 && (candType==1 && (MCType>=13&&MCType<=14))) return true;
  if (sam==18 && (candType==2 && (MCType>=13&&MCType<=14))) return true;
  if (sam==19 && (candType==3 && (MCType>=13&&MCType<=14))) return true;
  if (sam==20 && (candType==4 && (MCType>=13&&MCType<=14))) return true;

  if (sam==21 && (candType==5&&(MCType>=13&&MCType<=14))) return true;
  if (sam==22 && (candType==6&&(MCType>=13&&MCType<=14))) return true;
  if (sam==23 && (candType==7&&(MCType>=13&&MCType<=14))) return true;
  if (sam==24 && (candType==8&&(MCType>=13&&MCType<=14))) return true;
  if (sam==25 && (candType==5&&(MCType==2 || MCType== 4 || MCType== 6))) return true;
  if (sam==26 && (candType==6&&(MCType==2 || MCType== 4 || MCType== 6))) return true;
  if (sam==27 && (candType==7&&(MCType==8 || MCType==10 || MCType==12))) return true;
  if (sam==28 && (candType==8&&(MCType==8 || MCType==10 || MCType==12))) return true;
  if (sam==29 && (candType==5&&(MCType==1 || MCType== 3 || MCType== 5))) return true;
  if (sam==30 && (candType==6&&(MCType==1 || MCType== 3 || MCType== 5))) return true;
  if (sam==31 && (candType==7&&(MCType==7 || MCType== 9 || MCType==11))) return true;
  if (sam==32 && (candType==8&&(MCType==7 || MCType== 9 || MCType==11))) return true;
  ///////////////// Dpipilnu /////////////////////////////////////////
  if (sam==33 && (candType==5&&(MCType==16))) return true;
  if (sam==34 && (candType==6&&(MCType==16))) return true;
  if (sam==35 && (candType==7&&(MCType==16))) return true;
  if (sam==36 && (candType==8&&(MCType==16))) return true;
  if (sam==37 && (candType==1&&(MCType==16))) return true;
  if (sam==38 && (candType==2&&(MCType==16))) return true;
  if (sam==39 && (candType==3&&(MCType==16))) return true;
  if (sam==40 && (candType==4&&(MCType==16))) return true;
  if (sam==49 && (candType>4&&(MCType==16))) return true;
  if (sam==50 && (candType<5&&(MCType==16))) return true;
  ///////////////// Dpipilnu /////////////////////////////////////////

  if (sam==41 && (candType==1&&MCType>6&&MCType<13)) return true;
  if (sam==42 && (candType==2&&MCType>6&&MCType<13)) return true;
  if (sam==43 && (candType==3&&MCType>0&&MCType<7)) return true;
  if (sam==44 && (candType==4&&MCType>0&&MCType<7)) return true;
  if (sam==45 && (candType==5&&MCType>6&&MCType<13)) return true;
  if (sam==46 && (candType==6&&MCType>6&&MCType<13)) return true;
  if (sam==47 && (candType==7&&MCType>0&&MCType<7)) return true;
  if (sam==48 && (candType==8&&MCType>0&&MCType<7)) return true;

  if (sam==51 && (candType==1&&MCType==0)) return true;
  if (sam==52 && (candType==2&&MCType==0)) return true;
  if (sam==53 && (candType==3&&MCType==0)) return true;
  if (sam==54 && (candType==4&&MCType==0)) return true;
  if (sam==55 && (candType==5&&MCType==0)) return true;
  if (sam==56 && (candType==6&&MCType==0)) return true;
  if (sam==57 && (candType==7&&MCType==0)) return true;
  if (sam==58 && (candType==8&&MCType==0)) return true;

  if (sam==61 && (candType==1&&MCType<0)) return true;
  if (sam==62 && (candType==2&&MCType<0)) return true;
  if (sam==63 && (candType==3&&MCType<0)) return true;
  if (sam==64 && (candType==4&&MCType<0)) return true;
  if (sam==65 && (candType==5&&MCType<0)) return true;
  if (sam==66 && (candType==6&&MCType<0)) return true;
  if (sam==67 && (candType==7&&MCType<0)) return true;
  if (sam==68 && (candType==8&&MCType<0)) return true;
  if (sam==69 && ((candType==1||candType==3)&&MCType<0)) return true;
  if (sam==70 && ((candType==2||candType==4)&&MCType<0)) return true;

  return false;
}

TString RoundNumber(double num, int decimals, double denom){
  if(denom==0 || !isfinite(num) || !isfinite(denom)) return " - ";
  double neg = 1; if(num*denom<0) neg = -1;
  num /= neg*denom; num += 0.5*pow(10.,-decimals);
  if(abs(num) > 1e16) return "-";
  long num_int = static_cast<long>(num);
  long num_dec = static_cast<long>((1+num-num_int)*pow(10.,decimals));
  TString s_dec = ""; s_dec += num_dec; s_dec.Remove(0,1);
  TString result="";
  if(neg<0) result+="-";
  result+= num_int;
  if(decimals>0) {
    result+="."; result+=s_dec;
  }

  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<decimals-afterdot.Length(); i++)
    result += "0";
  if(result.Length()>15) cout<<"num "<<num<<", denom "<<denom<<"  ---->  "<<result<<endl;
  return result;
}

void getNumberB(TString genName, TString runs, double &totMCB, double &totdata, 
		double &totuds, double &totccbar, double &totOffdata){
  double NBR24[2][6] = {{  34878000, 101690000,  56035000, 166784000, 215168000, 130336000},
			{  34941000, 104188000,  57888000, 169801000, 215953000, 135224000}};
  double NBR26[2][6] = {{  78999000, 235345000, 120771000, 389670000, 509712000, 290169000},
			{  78560000, 246656000, 122374000, 383657000, 526675000, 291862000}};
  double nccbar[] =     {  55254000, 164146000,  88321000, 267308000, 343667000, 208664000};
  double nuds[] =       { 160514000, 451636000, 275869000, 421599000, 553604000, 327032000};
  double ndata[] =      { 22556256.9, 68438426.0, 35763257.9, 111429669.4, 147620363.4, 85194672.2};
  double lumOffdata[] = {2621.575, 7030.748, 2495.880, 10228.272, 14546.087, 7887.3};
  double fracRest24 = 2/3., fracRest26 = 1.;
  totMCB = 0; totuds = 0; totccbar = 0; totdata = 0; totOffdata = 0;
  if(genName.Contains("R24Rest")) fracRest24 = 1/3.;
  if(genName.Contains("R26Rest")) fracRest26 = 1/7.;
  for(int i=0; i<6; i++){
    TString Run = ""; Run += i+1;
    if(runs=="All" || runs.Contains(Run)){
      totdata += ndata[i];			                     //  471.003M events
      if(genName.Contains("RAll") || genName.Contains("R24")) 
	totMCB += (NBR24[0][i]+NBR24[1][i])*fracRest24;              //  948.591M events
      if(genName.Contains("RAll") || genName.Contains("R26")) 
	totMCB += (NBR26[0][i]+NBR26[1][i])*fracRest26;              // 3274.450M events
    }
    totOffdata += lumOffdata[i];                                     //   44.810k pb
    totuds += nuds[i];				                     // 2190.254M events
    totccbar += nccbar[i];			                     // 1127.360M events
  }
}

TH1F binHisto(TH2F *h2, int nbins, double xlow, double xhigh, double miny, double maxy, 
	      TString hname, TString var){
  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4;
  int nm2bin = h2->GetNbinsX(); 
  int nplbin = h2->GetNbinsY(); 
  double h2xmin = m2min, h2xmax = m2max, h2ymin = plmin, h2ymax = plmax;
  int nh2xbin = nm2bin, nh2ybin = nplbin;
  if(var=="candPstarLep"){
    h2xmin = plmin, h2xmax = plmax; 
    h2ymin = m2min, h2ymax = m2max; 
    nh2xbin = nplbin; nh2ybin = nm2bin;}
  double h2width = (h2xmax-h2xmin)/static_cast<double>(nh2xbin);
  double h1width = (xhigh-xlow)/static_cast<double>(nbins);
  double h2yini = (miny-h2ymin)*static_cast<double>(nh2ybin)/(h2ymax-h2ymin)+1;
  double h2yfin = (maxy-h2ymin)*static_cast<double>(nh2ybin)/(h2ymax-h2ymin)+1;
  int h2yini_i = static_cast<int>(h2yini), h2yfin_i = static_cast<int>(h2yfin); 
  double iniyfrac = h2yini_i+1 - h2yini, finyfrac = h2yfin - h2yfin_i;
  if(h2yini_i<1) h2yini_i = 1;	        
  if(h2yfin_i>nh2ybin) h2yfin_i = nh2ybin;

  TH1F h1(hname,hname,nbins, xlow, xhigh);
  //cout<<endl<<"Binning from "<<xlow<<" to "<<xhigh<<" with "<<nbins<<" bins. Initial bins: "<<nh2xbin<<endl;
  for(int bin = 1; bin <= nbins; bin++){
    double binvar = 0, tempvar = 0;
    double binmin = xlow+static_cast<double>((bin-1))*h1width;
    double h2xini = (binmin-h2xmin)/h2width+1, h2xfin = (binmin+h1width-h2xmin)/h2width+1;
    int h2xini_i = static_cast<int>(h2xini), h2xfin_i = static_cast<int>(h2xfin); 
    double inifrac = h2xini_i+1 - h2xini, finfrac = h2xfin - h2xfin_i;
    if(h2xini_i<1) h2xini_i = 1;	        
    if(h2xfin_i>nh2xbin) h2xfin_i = nh2xbin;
    if(h2xini_i==h2xfin_i){inifrac = h1width/h2width; finfrac = 1.;}
    //cout<<"h2xini_i: "<<h2xini_i<<"\th2xfin_i: "<<h2xfin_i<<"\tinifrac "<<inifrac<<",\t";
    //cout<<"h2yini_i: "<<h2yini_i<<"\th2yfin_i: "<<h2yfin_i<<"\tiniyfrac "<<iniyfrac<<endl;
    for(int h2xbin = h2xini_i; h2xbin<h2xfin_i+1; h2xbin++){
      for(int h2ybin = h2yini_i; h2ybin < h2yfin_i+1; h2ybin++){
	if(var=="candPstarLep") tempvar = h2->GetBinContent(h2ybin, h2xbin);
	else tempvar = h2->GetBinContent(h2xbin, h2ybin);
	if(h2xbin == h2xini_i) tempvar *= inifrac;
	if(h2xbin == h2xfin_i) tempvar *= finfrac;
	if(h2ybin == h2yini_i) tempvar *= iniyfrac;
	if(h2ybin == h2yfin_i) tempvar *= finyfrac;
	binvar += tempvar;
      }
    }
    h1.SetBinContent(bin,binvar);
    h1.SetBinError(bin,0);
  }

  return h1;
}

void FillHisto(TH1F *h, TTree *tree, TString varName, int cand){
  int entries = tree->GetEntries(), candType;
  float varVal, weight;
  tree->SetBranchAddress(varName,&varVal);
  tree->SetBranchAddress("candType",&candType);
  tree->SetBranchAddress("weight",&weight);
  for (int evt = 0 ; evt < entries; evt ++) {
    tree->GetEvent(evt);
    if(candType==cand) h->Fill(varVal,weight);
  }
}
void PlotFinalFit(TString textName, TTree *tree, int nbins, double minx, double maxx,
		  TString typePlot, TString var, TString sample, 
		  double maxFactor, TString psEPS){
  double Yields[16][8];
  TString buffer;
  fstream textFile; textFile.open(textName,fstream::in);
  int begChan = 0, endChan = 4, IndText[8] = {1,5,0,4,3,7,2,6};
  if(typePlot.Contains("D0")&&!typePlot.Contains("Dc")) endChan = 2; 
  if(typePlot.Contains("Dc")&&!typePlot.Contains("D0")) begChan = 2;
  for(int i=2*begChan; i<2*endChan; i++){
    int index = IndText[i];
    for(int j=0; j<8; j++){
      textFile>>buffer>>Yields[index+8][j]>>Yields[index][j]>>buffer>>buffer>>buffer;
      cout<<"True: "<<Yields[index+8][j]<<"\tFit: "<<Yields[index][j]<<"\tChannel "<<index<<endl;
      if(index>3 && j==4){
	textFile>>buffer;
	break;
      }
    }
    cout<<endl;
  }
  gROOT->SetStyle("Plain");gStyle->SetOptStat(0);
  TCanvas can("can","Final fit",1200,800);
  can.Divide(2,2);
  int nm2bin = 0, nplbin = 0, Projbins = nbins;
  TH2F *pdf[70];
  for (int i = 1 ; i <= 58 ; i ++) {
    if(i==29)i=31; if(i==39)i=41; if(i==49)i=51; 
    TString hname = "/cms2r0/slac/keys/root/Fit/pdfKeys_"; hname += i; hname += "_Fit.root";
    TFile hfile(hname); 
    TString pdfName = "pdf"; pdfName += i;
    pdf[i] = static_cast<TH2F *>((hfile.Get("h2"))->Clone(pdfName));
    pdf[i]->SetDirectory(0);
    if(i==0) {nm2bin = pdf[0]->GetNbinsX(); nplbin = pdf[0]->GetNbinsY();}
    //if(nm2bin!=pdf[i]->GetNbinsX() || nplbin!=pdf[i]->GetNbinsY()) return can;
  }
  int nxbins = nm2bin; double miny=0, maxy=2.4;
  TString xTitle = "m^{2}_{miss} (Gev^{2})";
  if(var=="candPstarLep") {
    nxbins = nplbin; miny=-4; maxy=12;
    xTitle = "p*_{l} (Gev)";}
  if(typePlot.Contains("curve") || typePlot.Contains("All")) Projbins = nxbins;
  int Indpdf[8][8] = {{51,41,31,17,9,13,5,1},{52,42,32,18,10,14,6,2},
		      {53,43,33,19,11,15,7,3}, {54,44,34,20,12,16,8,4},
		      {55,45,35,25,21,-1,-1,-1},{56,46,36,26,22,-1,-1,-1},
		      {57,47,37,27,23,-1,-1,-1},{58,48,38,28,24,-1,-1,-1}};
  int Colors[7] = {5,14,28,2,4,3,6};
  TString legNames[7] = {"Cont.","Bkg.","D**l#nu","Dl#nu","D*l#nu","D*#tau#nu","D#tau#nu"};
  int IndColor[3][8] = {{0,1,1,2,3,4,5,6},{0,1,1,3,2,4,6,5},{0,1,1,4,2,-1,-1,-1}};
  int IndYield[3][8] = {{7,6,5,4,1,3,2,0},{7,6,5,4,1,3,2,0},{4,3,2,1,0,-1,-1,-1}};
  TString channelTitle[4]={"D*^{0} channel","D^{0} channel","D*^{+} channel","D^{+} channel"};
  int truInd=0;
  if(typePlot.Contains("Tru")) truInd=8;
  TH1F *htree[4], pdfProj[4][8];
  TLegend *leg[4];
  TLatex label; label.SetTextSize(0.075); label.SetNDC(kTRUE);
  double weird = Yields[7][0];
  for(int pad=begChan; pad<endChan; pad++){
    can.cd(pad+1);
    int type = 0;
    if(sample!="dss") {if(pad%2==1) type=1;}
    else type=2;
    int channel = pad+4*(sample=="dss");
    if(channel<4){
      if(var=="candM2") leg[pad] = new TLegend(.71,.57,0.9,0.9);
      else leg[pad] = new TLegend(0.1,.57,.29,.9);
    }else leg[pad] = new TLegend(0.71,0.7,0.9,0.9);
    leg[pad]->SetFillColor(kWhite);
    leg[pad]->SetTextSize(0.05);
    TString hname = "hData"; hname += pad; 
    TString plotvar = var; plotvar+=">>"; plotvar+=hname; //plotvar+="(";plotvar+=nbins;
    //plotvar+=",";plotvar+=minx;plotvar+=",";plotvar+=maxx;plotvar+=")";
    int candType = pad; if(pad%2==0)candType+=2; if(sample=="dss")candType+=4;
    TString cuts = "(candType=="; cuts+=candType; cuts+=")*weight";

    htree[pad] = new TH1F(hname,hname,nbins,minx,maxx);
    FillHisto(htree[pad],tree,var,candType);
//     tree->Project(hname,var,cuts); // Element Yields[7][0] changes value!!
    Yields[7][0] = weird;
    formatHisto(htree[pad]);
    htree[pad]->GetXaxis()->SetTitle(xTitle);
    htree[pad]->SetMaximum(htree[pad]->GetMaximum()*1.15/maxFactor);
    htree[pad]->Draw("e1");
    for(int i=0; i<8; i++) {
      TString pdfName = "pdfProj"; pdfName += pad; pdfName += i;
      double scale = Yields[channel+truInd][IndYield[type][i]]*static_cast<double>(Projbins)/
	pdf[Indpdf[channel][i]]->Integral()/static_cast<double>(nbins);
      //cout<<"Yield "<<Yields[channel+truInd][IndYield[type][i]]<<" for "<< channel+truInd<<" and "<<
      //IndYield[type][i] <<", scale "<<scale<<endl;
      pdfProj[pad][i] = binHisto(pdf[Indpdf[channel][i]],Projbins,minx,maxx,miny,maxy,pdfName,var);
      pdfProj[pad][i].Scale(scale);
      if(i>0) pdfProj[pad][i].Add(&pdfProj[pad][i-1]);
      pdfProj[pad][i].SetFillColor(Colors[IndColor[type][i]]);
      pdfProj[pad][i].SetLineColor(Colors[IndColor[type][i]]);
      if(i==0) pdfProj[pad][i].SetLineColor(1);
      if(sample=="dss" && i==4) break;
    }
    TString optionPlot = "c same";
    if(!typePlot.Contains("curve")) optionPlot = "same";
    for(int i=7; i>=0; i--) {
      if(sample=="dss" && i==7) i=4;
      pdfProj[pad][i].Draw(optionPlot);
      if(i>0)leg[pad]->AddEntry(&pdfProj[pad][i],legNames[IndColor[type][i]]);
    }
    label.DrawLatex(0.1,0.92,channelTitle[pad]);
    TString pos = "right"; if(sample!="dss" && var=="candPstarLep") pos="left";
    leg[pad]->Draw();
    htree[pad]->Draw("e0 same");
  }
  if(psEPS.Contains(".ps")) {
    can.Print(psEPS);
  }
  if(psEPS.Contains("EPS")) {
    TString plotName = "plots/"; plotName+=sample; plotName+=var; plotName+="_"; 
    plotName+=typePlot; plotName+=".pdf";
    can.SaveAs(plotName);
  }
  cout<<"KeysUtils::PlotFinalFit - Deleting leg, htree, pdf"<<endl;
  for(int i=begChan; i<endChan;i++) {
    leg[i]->Delete();
    htree[i]->Delete();
  }
  for(int i=0; i<70;i++) if(pdf[i]) pdf[i]->Delete();
  return;
}

TH1F *GaussConv(TH1F *h, float width, TString hname, int PointsBin, int Cnbins, double CminX, double CmaxX,
		double Goffset) {
  int Hnbins = h->GetNbinsX();
  double HminX = h->GetXaxis()->GetBinLowEdge(h->GetXaxis()->GetFirst());
  double HmaxX = h->GetXaxis()->GetBinLowEdge(h->GetXaxis()->GetLast()+1);
  TH1F *hConv;
  if(Cnbins<0) {
    hConv = static_cast<TH1F *>(h->Clone(hname));
    Cnbins = Hnbins; CminX = HminX; CmaxX = HmaxX;
  } else hConv = new TH1F(hname, hname, Cnbins, CminX, CmaxX);
  if(CminX<HminX || CmaxX>HmaxX){
    cout<<"Range of convolved histogram has to be a subset of the original"<<endl;
    return 0;
  }
  double nSigma = 3.;
  double CwBin = (CmaxX-CminX)/static_cast<double>(Cnbins), HwBin = (HmaxX-HminX)/static_cast<double>(Hnbins);

  double val[2000], dx = CwBin/(static_cast<double>((PointsBin+1))), totAreaG = IntG(0,1,-nSigma,nSigma);
  for(int bin=0; bin<Hnbins; bin++) val[bin] = h->GetBinContent(bin+1);
  for(int binC=1; binC<=Cnbins; binC++) {
    double valC = 0;
    for(int point=0; point<PointsBin; point++){
      double x = CminX+static_cast<double>((binC-1))*CwBin+static_cast<double>((point+1))*dx+Goffset;
      double minG = x-static_cast<double>(nSigma)*width, maxG = x+static_cast<double>(nSigma)*width;
      double AreaG = totAreaG, valH = 0;
      bool lessRange = false;
      if(minG<HminX) {minG=HminX; lessRange = true;} if(maxG>HmaxX) {maxG=HmaxX; lessRange = true;} 
      if(lessRange) AreaG = IntG(x,width,minG,maxG);
      int iniBin = static_cast<int>(((minG-HminX)/HwBin))+1, finBin = static_cast<int>(((maxG-HminX)/HwBin))+1;
      for(int binH=iniBin; binH<=finBin; binH++){
	double minH = static_cast<double>((binH-1))*HwBin+HminX, maxH = static_cast<double>((binH))*HwBin+HminX;
	if(minH<minG) minH=minG; if(maxH>maxG) maxH=maxG; 
	valH += IntG(x, width, minH, maxH)*val[binH-1];
      }
      valC += valH/(AreaG*static_cast<double>(PointsBin));
    }
    hConv->SetBinContent(binC,valC);
  }
  return hConv;
}

void ConvKeys(TH2F *h2, double width){
  int nM2bin = h2->GetNbinsX(), nPlbin = h2->GetNbinsY(); 
  double minX = h2->GetXaxis()->GetBinLowEdge(h2->GetXaxis()->GetFirst());
  double maxX = h2->GetXaxis()->GetBinLowEdge(h2->GetXaxis()->GetLast()+1);

  int PointsBin = 1;
  double nSigma = 3., val[2000];
  double wBin = (maxX-minX)/static_cast<double>(nM2bin);
  double dx = wBin/(static_cast<double>((PointsBin+1))), totAreaG = IntG(0,1,-nSigma,nSigma);


  for(int bin=1; bin<=nPlbin; bin++){
    for(int mbin=1; mbin<=nM2bin; mbin++) val[mbin-1] = h2->GetBinContent(mbin,bin);
    for(int mbin=1; mbin<=nM2bin; mbin++) {
      double valC = 0;
      for(int point=0; point<PointsBin; point++){
	double x = minX+static_cast<double>((mbin-1))*wBin+static_cast<double>((point+1))*dx;
	double minG = x-static_cast<double>(nSigma)*width, maxG = x+static_cast<double>(nSigma)*width;
	double AreaG = totAreaG, valH = 0;
	bool lessRange = false;
	if(minG<minX) {minG=minX; lessRange = true;} if(maxG>maxX) {maxG=maxX; lessRange = true;} 
	if(lessRange) AreaG = IntG(x,width,minG,maxG);
	int iniBin = static_cast<int>(((minG-minX)/wBin))+1, finBin = static_cast<int>(((maxG-minX)/wBin))+1;
	for(int binH=iniBin; binH<=finBin; binH++){
	  double minH = static_cast<double>((binH-1))*wBin+minX, maxH = static_cast<double>((binH))*wBin+minX;
	  if(minH<minG) minH=minG; if(maxH>maxG) maxH=maxG; 
	  valH += IntG(x, width, minH, maxH)*val[binH-1];
	}
	valC += valH/(AreaG*static_cast<double>(PointsBin));
      }
      h2->SetBinContent(mbin,bin,valC);
    }
  }
}

void PlotEps(TString hname, TString dataFolder, TString epsFolder, TString weightName, TString psEPS,
	     double maxFactor, int nbins, double minM2, double maxM2, 
	     double intFactor, TString typePlot) {
  if(weightName=="") cout<<"Needed to use it for compiling"<<endl;
  TString sam, smoo, scaleP, Base;
  Int_t Sam, nM2bin, nPlbin;
  ParseName(hname, sam, smoo, scaleP, Base, Sam, nM2bin, nPlbin);
  bool doPS = false, doEPS = false;
  if(psEPS.Contains("EPS")) doEPS = true;
  if(psEPS.Contains("ps")) doPS = true;

  styles style6; 
  //style6.setPadsStyle(6); 
  style6.setDefaultStyle();
  TString inputfile = nameData(Sam, dataFolder);
  TChain *treeData = new TChain("ntp1");
  treeData->Add(inputfile);
  TFile hfile(hname); hfile.cd();
  TH2F *h2 = static_cast<TH2F *>(gDirectory->Get("h2"));
  nM2bin = h2->GetNbinsX(); nPlbin = h2->GetNbinsY(); 
  if(typePlot.Contains("Conv")) ConvKeys(h2,intFactor);
  intFactor = 1;
  TCanvas mm("mm","KEYS fits to mmiss-pstarl");
  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4;
  double xlow = m2min,xhigh = m2max, ylow = plmin,yhigh = plmax;
  Int_t nbinx = 80,nbiny = 80;
  setBins(Sam, xlow, xhigh, nbinx, nbiny);
  if(nbins>0){
    nbinx = nbins;
    nbiny = nbins;
  }
  if(minM2>-50) xlow = minM2; if(maxM2>-50) xhigh = maxM2; 
  if(!typePlot.Contains("curve") && !typePlot.Contains("All")) {
    nM2bin = nbinx; nPlbin = nbiny;
    m2min = xlow; m2max = xhigh;
    plmin = ylow; plmax = yhigh;
  }

  TString M2titles[] = {"p*_{l} < 1 GeV","1 #leq p*_{l} < 1.4 GeV","1.4 #leq p*_{l} < 1.8 GeV",
		      "1.8 #leq p*_{l} < 2.4 GeV","0 < p*_{l} < 2.4 GeV"};
  TCut M2cuts[] = {"candPstarLep<1","candPstarLep>=1&&candPstarLep<1.4",
		      "candPstarLep>=1.4&&candPstarLep<1.8","candPstarLep>=1.8&&candPstarLep<2.4", ""};
  TCut Plcuts[] = {"candM2<1","candM2>=1",""};
  double PlLimits[5][2] = {{0,1},{1,1.4},{1.4,1.8},{1.8,2.4},{0,2.4}};
  double M2Limits[3][2] = {{-4., 1.},{1.,12.},{-4.,12.}};

  TH1F hm2[5], *m2[5], hpl[3], *pl[3];
  for(int i=0;i<5;i++){
    TString hm2Name = "hm2"; hm2Name += i;
    hm2[i] = binHisto(h2,nM2bin,m2min,m2max,PlLimits[i][0],PlLimits[i][1],hm2Name,"candM2"); 
    if(i==2) {
      TString hplName = "hpl"; hplName += i;
      hpl[i] = binHisto(h2,nPlbin,plmin,plmax,M2Limits[i][0],M2Limits[i][1],hplName,"candPstarLep"); 
    }
  }
  double entries = getWentries(treeData);
  if(entries<0.00001) entries = 1;
  double scaleM2 = entries*nM2bin*(xhigh-xlow)/nbinx/(m2max-m2min)/h2->Integral()*intFactor;
  double scalePl = entries*nPlbin/nbiny/h2->Integral()*intFactor;
  for(int i=0;i<5;i++){
    TString hdname = "dm2"; hdname += i;
    TString vari = "candM2>>"; vari+=hdname; 
    M2cuts[i] *= "weight";
    m2[i] = new TH1F(hdname, "", nbinx, xlow, xhigh);
    m2[i]->Sumw2();
    if(treeData) treeData->Draw(vari,M2cuts[i]);
    style6.setMarkers(m2[i], 0.7, 20);
    hm2[i].Scale(scaleM2);
    hm2[i].SetLineColor(4);
    hm2[i].SetLineWidth(1);
    hm2[i].SetAxisRange(xlow,xhigh);
    float maxHisto = hm2[i].GetMaximum();
    if(treeData){
      if(typePlot == "diff") m2[i]->Add(&hm2[i],-1);
      else {
	if(m2[i]->GetMaximum()>maxHisto) maxHisto = m2[i]->GetMaximum();
	if(maxHisto<1) maxHisto = 1.;
	m2[i]->SetMaximum(maxHisto/maxFactor*1.2);
      }
    }
  }
  for(int i=2;i<3;i++){
    TString hdname = "pl"; hdname += i;
    TString vari = "candPstarLep>>"; vari+=hdname; 
    Plcuts[i] *= "weight";
    pl[i] = new TH1F(hdname, "", nbiny, ylow, yhigh);
    pl[i]->Sumw2();
    if(treeData) treeData->Draw(vari,Plcuts[i]);
    style6.setMarkers(pl[i], 0.7, 20);
    hpl[i].Scale(scalePl);
    hpl[i].SetLineColor(4);
    hpl[i].SetLineWidth(1);
    hpl[i].SetAxisRange(ylow,yhigh);
    float maxHisto = hpl[i].GetMaximum();
    if(treeData){
      if(typePlot == "diff") pl[i]->Add(&hpl[i],-1);
      else {
	if(pl[i]->GetMaximum()>maxHisto) maxHisto = pl[i]->GetMaximum();
	pl[i]->SetMaximum(maxHisto*1.2);
      }
    }
  }
  TString optionPlot = "c same";
  TLine lin; lin.SetLineColor(4);
  if(!typePlot.Contains("curve")) optionPlot = "same";
  if(doPS){
    TString psname = "plots/"; psname += epsFolder; psname += "/psKeys"; 
    psname += Base; psname += ".ps";
    mm.Print(psname+"[");
    if(treeData){m2[4]->Draw("e0"); m2[4]->Draw("e1 same"); hm2[4].Draw(optionPlot); 
    }else hm2[4].Draw("c"); 
    mm.Print(psname);
    if(treeData){pl[2]->Draw("e0"); pl[2]->Draw("e1 same"); hpl[2].Draw(optionPlot); 
    }else hpl[2].Draw("c"); 
    mm.Print(psname);
    for(int i=0;i<4;i++){
      if(treeData){m2[i]->Draw("e0"); m2[i]->Draw("e1 same"); hm2[i].Draw(optionPlot); 
      }else hm2[i].Draw("c"); 
      mm.Print(psname);
    }
    for(int i=0;i<2;i++){
      if(treeData){pl[i]->Draw("e0"); pl[i]->Draw("e1 same"); hpl[i].Draw(optionPlot); 
      }else hpl[i].Draw("c"); 
      mm.Print(psname);
    }
    mm.Print(psname+"]");
  }
  
  if(doEPS){
    TString rightLabel[] = {"",""};
    if(typePlot.Contains("histo")){
      double chi2[2] = {0,0}; int ndof[2] = {-1,-1};
      for(int bin=1; bin<=nbinx; bin++){
	if(m2[4]->GetBinContent(bin)<10) continue;
	double diff = m2[4]->GetBinContent(bin)-hm2[4].GetBinContent(bin);
	double uncert = m2[4]->GetBinError(bin);
	if(uncert>0){
	  chi2[0] += pow(diff/uncert,2);
	  ndof[0]++;
	}
      }
      for(int bin=1; bin<=nbiny; bin++){
	if(pl[2]->GetBinContent(bin)<10) continue;
	double diff = pl[2]->GetBinContent(bin)-hpl[2].GetBinContent(bin);
	double uncert = pl[2]->GetBinError(bin);
	if(uncert>0){
	  chi2[1] += pow(diff/uncert,2);
	  ndof[1]++;
	}
      }
      rightLabel[0] = "#chi^{2} = "; rightLabel[0] += RoundNumber(chi2[0],1); rightLabel[0] += "/";
      rightLabel[0] += ndof[0]; rightLabel[0] += "; p = "; rightLabel[0] += RoundNumber(TMath::Prob(chi2[0],ndof[0])*100,1); rightLabel[0] += "%";    
      rightLabel[1] = "#chi^{2} = "; rightLabel[1] += RoundNumber(chi2[1],1); rightLabel[1] += "/";
      rightLabel[1] += ndof[1]; rightLabel[1] += "; p = "; rightLabel[1] += RoundNumber(TMath::Prob(chi2[1],ndof[1])*100,1); rightLabel[1] += "%";    
    }
    TString Mytitle = "Entries/("; Mytitle += RoundNumber(xhigh-xlow,2,static_cast<double>(nbinx)); Mytitle += " GeV^{2})";
    TString Pytitle = "Entries/("; Pytitle += RoundNumber((yhigh-ylow)*1000,0,static_cast<double>(nbiny)); Pytitle += " MeV)";
    TCanvas all6("all6","KEYS fits to mmiss-pstarl");
    all6.Divide(2,3);
    all6.cd(1);
    if(treeData){
      m2[4]->Draw("e0");  if(typePlot.Contains("diff"))lin.DrawLine(xlow,0,xhigh,0);
      m2[4]->Draw("e1 same"); 
      if(!typePlot.Contains("diff")) hm2[4].Draw(optionPlot);
      style6.setTitles(m2[4],"m^{2}_{miss} (GeV^{2})",Mytitle,"",rightLabel[0]);
    }else hm2[4].Draw("c");
    double legW = 0.32, legH = 0.2;
    double legX = 1-style6.PadRightMargin-0.01, legY = 1-style6.PadTopMargin-0.008;
    TLegend leg(legX-legW, legY-legH, legX, legY);
    leg.SetFillColor(0); leg.SetTextSize(style6.TitleSize); leg.SetBorderSize(0);
    leg.SetTextFont(style6.nFont); 

    // if(style6.isThesis=="no") leg.AddEntry(m2[4],dataFolder);
    // else leg.AddEntry(m2[4],"MC");
    leg.AddEntry(m2[4],"MC");

    leg.AddEntry(&hm2[4],"KEYS p.d.f.");
    if(rightLabel[0]=="") leg.Draw();
    all6.cd(2);
    if(treeData){
      pl[2]->Draw("e0");  if(typePlot == "diff")lin.DrawLine(ylow,0,yhigh,0);
      pl[2]->Draw("e1 same"); if(typePlot != "diff")hpl[2].Draw(optionPlot);
      style6.setTitles(pl[2],"p*_{l} (GeV)",Pytitle,"",rightLabel[1]);
    }else hpl[2].Draw("c");
    for(int i=0;i<4;i++){
      all6.cd(i+3);
      if(treeData){
	m2[i]->Draw("e0"); if(typePlot == "diff")lin.DrawLine(xlow,0,xhigh,0);
	m2[i]->Draw("e1 same"); if(typePlot != "diff")hm2[i].Draw(optionPlot); 
	style6.setTitles(m2[i],"m^{2}_{miss} (GeV^{2})",Mytitle,"",M2titles[i]);
      }else hm2[i].Draw("c"); 
    }
    TString epsname = "plots/"; epsname += epsFolder; epsname += "/"; epsname += dataFolder; 
    epsname += "epsKeys"; epsname += Base; epsname += ".pdf";
    all6.SaveAs(epsname);
  }

  h2->Delete();
  return; 

}

// Cholesky-Banachiewicz algorithm
TMatrixT<double> Choleski(TMatrixT<double> A, int nRows){
  TMatrixT<double> L(nRows,nRows);
  
  for(int row=0; row<nRows; row++)
    for(int col=0; col<nRows; col++)
      L(row,col) = 0;
  
  for(int row=0; row<nRows; row++){
    for(int col=0; col<(row+1); col++){
      if(col==row) {
	for(int irow=0; irow<col; irow++) L(row,col) += L(col,irow)*L(col,irow);
	L(row,col) = sqrt(A(row,col) - L(row,col));
      } else {
	for(int irow=0; irow<col; irow++) L(row,col) += L(row,irow)*L(col,irow);
	L(row,col) = (A(row,col) - L(row,col))/L(col,col);
      }
    }
  }
  return L;
}


// Reads a text file with the results of a fit
void ReadFitFile(TString textName, double Yield[2][70], double Error[70]){
  TString buffer, sError;
  fstream textFile; textFile.open(textName,fstream::in);
  for(int i=0; i<8; i++){
    for(int j=0; j<9; j++){
      textFile>>buffer;
      buffer.ReplaceAll("Yield[",""); buffer.ReplaceAll("]",""); 
      int index = buffer.Atoi();
      textFile>>Yield[1][index]>>Yield[0][index]>>buffer>>sError>>buffer;
      if(sError.Atof()>0) Error[index] = sError.Atof();
      else Error[index] = 0;
      //cout<<"True: "<<Yield[1][index]<<"\tFit: "<<Yield[0][index]<<" +- "<<Error[index]<<"\tIndex "<<index<<endl;
      if(i%2==1 && j==6){
	textFile>>buffer;
	break;
      }
    }
    //cout<<endl;
  }

}

