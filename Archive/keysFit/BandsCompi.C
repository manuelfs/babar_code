#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TObject.h"
#include "TKey.h"
#include "TCollection.h"
#include "TStyle.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;
using std::cout;
using std::endl;

void bands(TString hname){
  TString sam = hname; sam.Remove(0,sam.First('_')+1);
  sam.Remove(sam.First('_'),sam.Length());
  TString smoo = hname;smoo.Remove(0,smoo.Last('_')+1);
  smoo.Remove(smoo.First('.'),smoo.Length());
  TString binString = hname; binString.Remove(binString.Last('_'),binString.Length());
  binString.Remove(0,binString.Last('_')+1);
  int nM2bin = binString.Atoi(); int nPlbin = nM2bin%1000;
  nM2bin = nM2bin/1000;
  TString basename(hname);
  Int_t slashpos = basename.Last('D');
  if (slashpos>=0) {
    basename.Remove(0,slashpos+1); basename.Remove(basename.First('.'),basename.Length());
  } else { cout<<hname<<" is incorrect"<<endl; return;}

  TCanvas mm("Bands","Bands");
  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4, tot=-999.;
  double PlLimits[] = {0, 1, 1.4, 1.8, 2.4};
  double M2Limits[] = {0., 5., 16.};
  int binlim[5]; int plbinlim[3];
  for(int i=0;i<5;i++) {
    PlLimits[i] = PlLimits[i]/2.4*(double)nPlbin;
    binlim[i] = (int)PlLimits[i];
    if(i<3){
      M2Limits[i] = M2Limits[i]/16.*(double)nM2bin;
      plbinlim[i] = (int)M2Limits[i];
    }
  }
  double M[5][1001], M2[5][1001];
  double p[3][301], p2[3][301];
  for(int k=1; k<nM2bin+1; k++){
    for(int i=0; i<5; i++){
      M[i][k] = 0; M2[i][k] = 0;
      if(k<nPlbin+1 && i<3) {p[i][k] = 0; p2[i][k] = 0;}
    }
  }

  TH2F *h2i; int nHisto=0;
  TFile f(hname);
  TIter next(f.GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)next())) {
    TObject* obj= key->ReadObj();
    if(nHisto%20==0) cout<<"Doing "<<nHisto<<" of some"<<endl;
    h2i = (TH2F *)gDirectory->Get(obj->GetName());
    for(int k=1; k<nM2bin+1; k++){
      double sumPlep[] = {0,0,0,0,0};
      for(int j=1; j<nPlbin+1; j++){
	double val = h2i->GetBinContent(k,j);
	for(int s=0; s<4; s++){
	  if(j>binlim[s] && j<binlim[s+1]+1) {
	    M[s][k] += val;
	    sumPlep[s] += val;
	  }
	  if(j==binlim[s+1]) M2[s][k] += sumPlep[s]*sumPlep[s];
	}
	M[4][k] += val;
	sumPlep[4] += val;
      }
      M2[4][k] += sumPlep[4]*sumPlep[4];
    }
    for(int k=1; k<nPlbin+1; k++){
      double sumM2[] = {0,0,0};
      for(int j=1; j<nM2bin+1; j++){
	double valp = h2i->GetBinContent(j,k);
	for(int s=0; s<2; s++){
	  if(j>plbinlim[s] && j<plbinlim[s+1]+1) {
	    p[s][k] += valp;
	    sumM2[s] += valp;
	  }
	  if(j==plbinlim[s+1]) p2[s][k] += sumM2[s]*sumM2[s];
	}
	p[2][k] += valp;
	sumM2[2] += valp;
      }
      p2[2][k] += sumM2[2]*sumM2[2];
    }
    if(nHisto==0) tot = h2i->Integral();
    nHisto++;
    if(nHisto==10) break;
  }
  f.Close();

  TString inputfile = "fitSamples/pdfSample"; inputfile += sam; inputfile += ".root";
  cout << "File = " << inputfile << endl;	
  TChain c("ntp1");
  c.Add(inputfile);
  double entries = c.GetEntries();
  Float_t xlow,xhigh;
  Int_t nbinx,nbiny,Sam = sam.Atoi();
  xlow = m2min;
  xhigh = m2max;
  nbinx = 80;
  nbiny = 80;
  if (Sam==0 || Sam==2 || Sam==10 || Sam==12 || Sam == 20 || Sam == 23 || Sam == 26 || Sam == 29) {
    xlow = -2; xhigh = 4;
    if (Sam > 12) {nbinx = 40; nbiny = 40;}
  }
  else if (Sam==1 || Sam==11) {
    xlow = -2; xhigh = 6;
  }
  if (Sam==6 || Sam==7 || Sam==16 || Sam==17) {
    nbinx = 40; nbiny = 40;
  }
  if (Sam==8 || Sam==18) {
    xhigh = 4; nbinx = 40; nbiny = 40;
  }
  if (Sam==9 || Sam==19) {
    nbinx = 40; nbiny = 40;
  }
  if (Sam==21 || Sam==22 || Sam==24 || Sam==25 || Sam==27 || Sam==28 || Sam==30 || Sam==31) {
    xhigh = 8; nbinx = 40; nbiny = 20;
  }
  if (Sam > 31) {
    nbinx = 40; nbiny = 20;
  }

  TString M2titles[] = {"0 < p*_{l} < 1 GeV","1 < p*_{l} < 1.4 GeV","1.4 < p*_{l} < 1.8 GeV",
		      "1.8 < p*_{l} < 2.4 GeV","0 < p*_{l} < 2.4 GeV"};
  TString Pltitles[] = {"-4 < m^{2}_{miss} < 1 GeV^{2}","1 < m^{2}_{miss} < 12 GeV^{2}",
			"-4 < m^{2}_{miss} < 12 GeV^{2}"};
  TString M2cuts[] = {"candPstarLep<1","candPstarLep>1&&candPstarLep<1.4",
		      "candPstarLep>1.4&&candPstarLep<1.8","candPstarLep>1.8&&candPstarLep<2.4", ""};
  TString Plcuts[] = {"candM2<1.","candM2>=1.",""};

  TString psname = "AWG82/results/keys/Bands/epsBandKeys"; psname+=basename; psname += ".ps";
  mm.Print(psname+"[");
  TH1F *hm2[5], *hm2Error[5], *m2[5], *hpl[3], *pl[3], *hplError[3];
  TString M2names[5], Plnames[3];
  for(int i=0;i<5;i++){
    M2names[i] = "hm2_"; M2names[i] += i;
    hm2[i] = new TH1F(M2names[i],M2titles[i],nM2bin,m2min,m2max); 
    M2names[i] += "Max";
    hm2Error[i] = new TH1F(M2names[i],M2titles[i],nM2bin,m2min,m2max); 
    M2names[i] += "Min";
    if(i<3) {
      Plnames[i] = "hpl_"; Plnames[i] += i;
      hpl[i] = new TH1F(Plnames[i],Pltitles[i],nPlbin,plmin,plmax); 
      Plnames[i] += "Max";
      hplError[i] = new TH1F(Plnames[i],Pltitles[i],nPlbin,plmin,plmax); 
    }
  }    
  for(int i=0;i<5;i++){
    TString hname = "m2"; hname += i;
    TString vari = "candM2>>"; vari+=hname; vari+="("; vari+= nbinx; vari+=",";vari+= xlow; 
    vari+=",";vari+= xhigh; vari+=")";
    c.Draw(vari,M2cuts[i]);
    m2[i] = (TH1F*)gDirectory->Get(hname);
    m2[i]->SetXTitle("m^{2}_{miss} [GeV^{2}]");
    m2[i]->SetTitle(M2titles[i]);
    m2[i]->Sumw2();
    m2[i]->SetMarkerStyle(20);
    m2[i]->SetMarkerSize(1);
    gStyle->SetOptStat(0);
    if(i<4){
      for(int j=1; j<nM2bin+1; j++){
	double mean, Square, variance, rms;
	double factor = entries*nM2bin*(xhigh-xlow)/nbinx/(m2max-m2min)/tot;
	mean = M[i][j]*factor/nHisto;
	Square = M2[i][j]*factor*factor;
	variance = Square/(nHisto-1)-mean*mean*nHisto/(nHisto-1);
	rms = 0;
	if(variance>=0) rms = sqrt(variance);
	else cout<<"Bin "<<j<<" is wrong. Square "<<Square<<" and mean "<<mean<<endl;
	hm2Error[i]->SetBinContent(j,mean);
	hm2Error[i]->SetBinError(j,rms);
	hm2[i]->SetBinContent(j,mean);
	//if(j%50==0)cout<<"mean "<<mean<<", rms "<<rms<<" at "<<j<<". Square "<<Square<<", mean "<<mean<<endl;
      }
    } 
    hm2[i]->SetLineColor(4);
    hm2[i]->SetLineWidth(2);
    hm2Error[i]->SetLineColor(9);
    hm2Error[i]->SetFillColor(38);
    hm2Error[i]->SetLineWidth(2);
    if(i<4) {
      hm2[4]->Add(hm2[i]);
      hm2Error[4]->Add(hm2Error[i]);
    }
  }
  for(int i=0;i<3;i++){
    TString hname = "pl"; hname += i;
    TString vari = "candPstarLep>>"; vari+=hname; vari+="("; vari+= nbiny; vari+=",0,2.4)";
    c.Draw(vari,Plcuts[i]);
    pl[i] = (TH1F*)gDirectory->Get(hname);
    pl[i]->SetXTitle("p*_{l} [GeV]");
    pl[i]->SetTitle(Pltitles[i]);
    pl[i]->Sumw2();
    pl[i]->SetMarkerStyle(20);
    pl[i]->SetMarkerSize(1);
    gStyle->SetOptStat(0);
    if(i<2){
      for(int j=1; j<nPlbin+1; j++){
	double mean, Square, variance, rms;
	double factor = entries*nPlbin/nbiny/tot;
	mean = p[i][j]*factor/nHisto;
	Square = p2[i][j]*factor*factor;
	variance = Square/(nHisto-1)-mean*mean*nHisto/(nHisto-1);
	rms = 0;
	if(variance>=0) rms = sqrt(variance);
	else cout<<"Bin "<<j<<" is wrong. Square "<<Square<<" and mean "<<mean<<endl;
	hplError[i]->SetBinContent(j,mean);
	hplError[i]->SetBinError(j,rms);
	hpl[i]->SetBinContent(j,mean);
      }
    }
    hpl[i]->SetLineColor(4);
    hpl[i]->SetLineWidth(2);
    hplError[i]->SetLineColor(9);
    hplError[i]->SetFillColor(38);
    hplError[i]->SetLineWidth(2);
    if(i<2) {
      hpl[2]->Add(hpl[i]);
      hplError[2]->Add(hplError[i]);
    }
  }
  m2[4]->Draw("e1"); hm2Error[4]->Draw("e4 same"); 
  m2[4]->Draw("e1 same"); hm2[4]->Draw("c same"); mm.Print(psname);
  pl[2]->Draw("e1"); hplError[2]->Draw("e4 same"); 
  pl[2]->Draw("e1 same"); hpl[2]->Draw("c same"); mm.Print(psname);
  for(int i=0;i<4;i++){
    m2[i]->Draw("e1");
    hm2Error[i]->Draw("e4 same"); 
    m2[i]->Draw("e1 same"); hm2[i]->Draw("c same"); 
    mm.Print(psname);
  }
  for(int i=0;i<2;i++){
    pl[i]->Draw("e1");
    hplError[i]->Draw("e4 same"); 
    pl[i]->Draw("e1 same"); hpl[i]->Draw("c same"); 
    mm.Print(psname);
  }
  mm.Print(psname+"]");
  return; 
}
