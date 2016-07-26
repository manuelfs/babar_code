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
#include "TLatex.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;
using std::cout;
using std::endl;

void formatHisto(TH1F *h);

void bands(TString hname){
  TString sam = hname; sam.Remove(0,sam.First('_')+1);
  sam.Remove(sam.First('_'),sam.Length()); int Sam = sam.Atoi();
  TString smoo = hname;smoo.Remove(0,smoo.Last('_')+1);
  smoo.Remove(smoo.First('.'),smoo.Length());
  TString binString = hname; binString.Remove(binString.Last('_'),binString.Length());
  binString.Remove(0,binString.Last('_')+1);
  int nM2bin = binString.Atoi(); int nPlbin = nM2bin%1000;
  nM2bin = nM2bin/1000;
  TString Base(hname);
  Int_t slashpos = Base.Last('D');
  if (slashpos>=0) {
    Base.Remove(0,slashpos+1); Base.Remove(Base.First('.'),Base.Length());
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
  vector<Double_t> M[5][1001];
  vector<Double_t> p[3][301];

  TH2F *h2i; int nHisto=0;
  TFile f(hname);
  TIter next(f.GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)next())) {
    TObject* obj= key->ReadObj();
    if(nHisto%20==0) cout<<"Doing "<<nHisto<<" of some"<<endl;
    h2i = (TH2F *)gDirectory->Get(obj->GetName());
    for(int k=1; k<nM2bin+1; k++){
      for(int j=1; j<nPlbin+1; j++){
	double val = h2i->GetBinContent(k,j);
	if(j==1) M[4][k].push_back(val);
	else M[4][k][nHisto] += val;
	if(k==1) p[2][j].push_back(val);
	else p[2][j][nHisto] += val;

	for(int s=0; s<4; s++){
	  if(j>binlim[s] && j<binlim[s+1]+1) {
	    if(j==binlim[s]+1) M[s][k].push_back(val);
	    else M[s][k][nHisto] += val;
	  }
	}
	for(int s=0; s<2; s++){
	  if(k>plbinlim[s] && k<plbinlim[s+1]+1) {
	    if(k==plbinlim[s]+1) p[s][j].push_back(val);
	    else p[s][j][nHisto] += val;
	  }
	}
      }
    }
    if(nHisto==0) tot = h2i->Integral();
    nHisto++;
    //if(nHisto==20) break;
  }
  f.Close();

  TString inputfile = "fitSamples/"; 
  if(Sam<20) inputfile+="pdfSample";
  else if(Sam<32){
    inputfile+="dssSample";
    sam = ""; sam+=Sam-20;
  }
  inputfile += sam; inputfile += ".root";
  cout << "File = " << inputfile << endl;	
  TChain c("ntp1");
  c.Add(inputfile);
  double entries = c.GetEntries();
  Float_t maxM2[]={0,0,0,0,0};
  double xlow = m2min,xhigh = m2max;
  Int_t nbinx = 80,nbiny = 80;
  if (Sam==0 || Sam==2 || Sam==10 || Sam==12 || Sam == 20 || Sam == 23 || Sam == 26 || Sam == 29) {
    xlow = -2; xhigh = 4;
    if (Sam > 12) {nbinx = 40; nbiny = 40;}
  } else if (Sam==1 || Sam==11) { xlow = -2; xhigh = 6;}
  if (Sam==6 || Sam==7 || Sam==16 || Sam==17) {nbinx = 40; nbiny = 40;}
  if (Sam==8 || Sam==18) {xhigh = 4; nbinx = 40; nbiny = 40;}
  if (Sam==9 || Sam==19) {nbinx = 40; nbiny = 40;}
  if (Sam==21 || Sam==22 || Sam==24 || Sam==25 || Sam==27 || Sam==28 || Sam==30 || Sam==31) {
    xhigh = 8; nbinx = 40; nbiny = 20;}
  if (Sam > 31) {nbinx = 40; nbiny = 20;}

  TString M2titles[] = {"0 < p*_{l} < 1 GeV","1 < p*_{l} < 1.4 GeV","1.4 < p*_{l} < 1.8 GeV",
		      "1.8 < p*_{l} < 2.4 GeV","0 < p*_{l} < 2.4 GeV"};
  TString Pltitles[] = {"-4 < m^{2}_{miss} < 1 GeV^{2}","1 < m^{2}_{miss} < 12 GeV^{2}",
			"-4 < m^{2}_{miss} < 12 GeV^{2}"};
  TString M2cuts[] = {"candPstarLep<1","candPstarLep>1&&candPstarLep<1.4",
		      "candPstarLep>1.4&&candPstarLep<1.8","candPstarLep>1.8&&candPstarLep<2.4", ""};
  TString Plcuts[] = {"candM2<1.","candM2>=1.",""};

  TString psname = "AWG82/results/keys/Bands/epsBandKeys"; psname+=Base; psname += ".ps";
  mm.Print(psname+"[");
  TH1F *hm2[5], *hm2Error[5], *hm2Error2[5], *m2[5], *hpl[3], *pl[3], *hplError[3], *hplError2[3];
  TString M2names[5], Plnames[3];
  for(int i=0;i<5;i++){
    M2names[i] = "hm2_"; M2names[i] += i;
    hm2[i] = new TH1F(M2names[i],M2titles[i],nM2bin,m2min,m2max); 
    M2names[i] += "sigma";
    hm2Error[i] = new TH1F(M2names[i],M2titles[i],nM2bin,m2min,m2max); 
    M2names[i] += "sigma";
    hm2Error2[i] = new TH1F(M2names[i],M2titles[i],nM2bin,m2min,m2max); 
    if(i<3) {
      Plnames[i] = "hpl_"; Plnames[i] += i;
      hpl[i] = new TH1F(Plnames[i],Pltitles[i],nPlbin,plmin,plmax); 
      Plnames[i] += "sigma";
      hplError[i] = new TH1F(Plnames[i],Pltitles[i],nPlbin,plmin,plmax); 
      Plnames[i] += "sigma";
      hplError2[i] = new TH1F(Plnames[i],Pltitles[i],nPlbin,plmin,plmax); 
    }
  }    
  for(int i=0;i<5;i++){
    TString hname = "m2"; hname += i;
    TString vari = "candM2>>"; vari+=hname; vari+="("; vari+= nbinx; vari+=",";vari+= xlow; 
    vari+=",";vari+= xhigh; vari+=")";
    c.Draw(vari,M2cuts[i]);
    m2[i] = (TH1F*)gDirectory->Get(hname);
    m2[i]->SetXTitle("m^{2}_{miss} [GeV^{2}]");
    formatHisto(m2[i]);
    gStyle->SetOptStat(0);
    for(int j=1; j<nM2bin+1; j++){
      double mean=0;
      double factor = entries*nM2bin*(xhigh-xlow)/nbinx/(m2max-m2min)/tot;


      for(int his=0; his<nHisto; his++) {
	mean += M[i][j][his];
      }
      mean = mean/nHisto;
      sort(M[i][j].begin(),M[i][j].end());
      double imean = 0.;
      for(int his=0; his<nHisto-1; his++) {
	double v0 = M[i][j][his], v1 = M[i][j][his+1];
	if(mean>v0 && mean<=v1){
	  if(v0==v1) imean = (double)(his+1);
	  else imean = (double)his+(mean-v0)/(v1-v0);
	  break;
	}
      }
      imean = (double)(nHisto-1)/2;
      double isigmas[] = {imean-0.4772*(double)nHisto, imean-0.3413*(double)nHisto, 
			  imean+0.3413*(double)nHisto, imean+0.4772*(double)nHisto, };
      if(isigmas[0]<0){
	isigmas[3] -= isigmas[0];
	isigmas[0] = 0.;
      }
      if(isigmas[1]<0){
	isigmas[2] -= isigmas[1];
	isigmas[1] = 0.;
      }
      if(isigmas[3]>(double)(nHisto-1)){
	isigmas[0] -= isigmas[3]-nHisto+1;
	isigmas[3] = (double)(nHisto-1);
      }
      if(isigmas[2]>(double)(nHisto-1)){
	isigmas[1] -= isigmas[2]-nHisto+1;
	isigmas[2] = (double)(nHisto-1);
      }
      double valsigmas[4];
      for(int sig=0; sig<4; sig++){
	int isig = (int)isigmas[sig];
	double w = isigmas[sig] - (double)isig;
	if(w==0) valsigmas[sig] = M[i][j][isig];
	else valsigmas[sig] = M[i][j][isig]*(1-w)+M[i][j][isig+1]*w;
	valsigmas[sig] *= factor;
      }
      mean *= factor;
      hm2Error[i]->SetBinContent(j,(valsigmas[2]+valsigmas[1])/2);
      hm2Error[i]->SetBinError(j,(valsigmas[2]-valsigmas[1])/2);
      hm2Error2[i]->SetBinContent(j,(valsigmas[3]+valsigmas[0])/2);
      hm2Error2[i]->SetBinError(j,(valsigmas[3]-valsigmas[0])/2);
      hm2[i]->SetBinContent(j,mean);
      if(maxM2[i]<valsigmas[3])maxM2[i]=valsigmas[3];
    }
    hm2[i]->SetLineColor(4);
    hm2[i]->SetLineWidth(1);
    hm2Error[i]->SetFillColor(38);
    hm2Error2[i]->SetFillColor(9);
    if(maxM2[i]>m2[i]->GetMaximum()) m2[i]->SetMaximum(maxM2[i]*1.05);
  }
  for(int i=0;i<3;i++){
    TString hname = "pl"; hname += i;
    TString vari = "candPstarLep>>"; vari+=hname; vari+="("; vari+= nbiny; vari+=",0,2.4)";
    c.Draw(vari,Plcuts[i]);
    pl[i] = (TH1F*)gDirectory->Get(hname);
    pl[i]->SetXTitle("p*_{l} [GeV]");
    formatHisto(pl[i]);
    gStyle->SetOptStat(0);
    for(int j=1; j<nPlbin+1; j++){
      double mean=0;
      double factor = entries*nPlbin/nbiny/tot;
      for(int his=0; his<nHisto; his++) {
	mean += p[i][j][his];
      }
      mean = mean/nHisto;
      sort(p[i][j].begin(),p[i][j].end());
      double imean = 0.;
      for(int his=0; his<nHisto-1; his++) {
	double v0 = p[i][j][his], v1 = p[i][j][his+1];
	if(mean>v0 && mean<=v1){
	  if(v0==v1) imean = (double)(his+1);
	  else imean = (double)his+(mean-v0)/(v1-v0);
	  break;
	}
      }
      imean = (double)(nHisto-1)/2;
      double isigmas[] = {imean-0.4772*(double)nHisto, imean-0.3413*(double)nHisto, 
			  imean+0.3413*(double)nHisto, imean+0.4772*(double)nHisto, };
      if(isigmas[0]<0){
	isigmas[3] -= isigmas[0];
	isigmas[0] = 0.;
      }
      if(isigmas[1]<0){
	isigmas[2] -= isigmas[1];
	isigmas[1] = 0.;
      }
      if(isigmas[3]>(double)(nHisto-1)){
	isigmas[0] -= isigmas[3]-nHisto+1;
	isigmas[3] = (double)(nHisto-1);
      }
      if(isigmas[2]>(double)(nHisto-1)){
	isigmas[1] -= isigmas[2]-nHisto+1;
	isigmas[2] = (double)(nHisto-1);
      }
      double valsigmas[4];
      for(int sig=0; sig<4; sig++){
	int isig = (int)isigmas[sig];
	double w = isigmas[sig] - (double)isig;
	if(w==0) valsigmas[sig] = p[i][j][isig];
	else valsigmas[sig] = p[i][j][isig]*(1-w)+p[i][j][isig+1]*w;
	valsigmas[sig] *= factor;
      }
      mean *= factor;
      if(j>25&&j<50 && i==-1) {
	cout<<"mean: "<<mean<<", imean: "<<imean<<", pl: "<<(double)j/(double)nPlbin*2.4<<endl;
	cout<<"isigmas[1]: "<<isigmas[1]<<", valsigmas[1]: "<<valsigmas[1]
	    <<", isigmas[2]: "<<isigmas[2]<<", valsigmas[2]: "<<valsigmas[2]<<endl;
	cout<<"isigmas[0]: "<<isigmas[0]<<", valsigmas[0]: "<<valsigmas[0]
	    <<", isigmas[3]: "<<isigmas[3]<<", valsigmas[3]: "<<valsigmas[3]<<endl;
      }
      hplError[i]->SetBinContent(j,(valsigmas[2]+valsigmas[1])/2);
      hplError[i]->SetBinError(j,(valsigmas[2]-valsigmas[1])/2);
      hplError2[i]->SetBinContent(j,(valsigmas[3]+valsigmas[0])/2);
      hplError2[i]->SetBinError(j,(valsigmas[3]-valsigmas[0])/2);
      hpl[i]->SetBinContent(j,mean);
    }
    hpl[i]->SetLineColor(4);
    hpl[i]->SetLineWidth(1);
    hplError[i]->SetFillColor(38);
    hplError2[i]->SetFillColor(9);
  }
  TLatex *label = new TLatex(); label->SetNDC(kTRUE); label->SetTextSize(0.055);
  m2[4]->Draw("e1"); hm2Error2[4]->Draw("e4 same"); hm2Error[4]->Draw("e4 same"); 
  m2[4]->Draw("e0 same"); m2[4]->Draw("e1 same"); hm2[4]->Draw("c same"); mm.Print(psname);
  pl[2]->Draw("e1"); hplError2[2]->Draw("e4 same"); hplError[2]->Draw("e4 same"); 
  pl[2]->Draw("e0 same"); pl[2]->Draw("e1 same"); hpl[2]->Draw("c same"); mm.Print(psname);
  for(int i=0;i<4;i++){
    m2[i]->Draw("e1");
    hm2Error2[i]->Draw("e4 same"); hm2Error[i]->Draw("e4 same"); 
    m2[i]->Draw("e0 same"); m2[i]->Draw("e1 same"); hm2[i]->Draw("c same"); 
    mm.Print(psname);
  }
  for(int i=0;i<2;i++){
    pl[i]->Draw("e1");
    hplError2[i]->Draw("e4 same"); hplError[i]->Draw("e4 same"); 
    pl[i]->Draw("e0 same"); pl[i]->Draw("e1 same"); hpl[i]->Draw("c same"); 
    mm.Print(psname);
  }
  mm.Print(psname+"]");
  
  TString sHisto = ""; sHisto+=nHisto; sHisto += " fits";
  bool doSigma2 = true;
  TCanvas all6("all6","KEYS fits to mmiss-pstarl",1700,1800);
  all6.cd();
  all6.Divide(2,3,0.001,0.001);
  all6.cd(1);gPad->SetTopMargin(0.004);gPad->SetRightMargin(0.003);gPad->SetBottomMargin(0.11);
  gPad->SetFillColor(10);gPad->SetFrameFillColor(10);
  m2[4]->Draw("e0"); if(doSigma2) hm2Error2[4]->Draw("e4 same");  
  hm2Error[4]->Draw("e4 same"); m2[4]->Draw("e0 same"); m2[4]->Draw("e1 same"); hm2[4]->Draw("c same");
  all6.cd(2);gPad->SetTopMargin(0.004);gPad->SetRightMargin(0.003);gPad->SetBottomMargin(0.11);
  gPad->SetFillColor(10);gPad->SetFrameFillColor(10);
  pl[2]->Draw("e0"); if(doSigma2) hplError2[2]->Draw("e4 same");  
  hplError[2]->Draw("e4 same"); pl[2]->Draw("e0 same"); pl[2]->Draw("e1 same"); hpl[2]->Draw("c same");
  for(int i=0;i<4;i++){
    all6.cd(i+3);gPad->SetTopMargin(0.003);gPad->SetRightMargin(0.003);gPad->SetBottomMargin(0.11);
    gPad->SetFillColor(10);gPad->SetFrameFillColor(10);
    m2[i]->Draw("e0");  if(doSigma2) hm2Error2[i]->Draw("e4 same");
    hm2Error[i]->Draw("e4 same"); label->DrawLatex(0.69,0.92,M2titles[i]);
    m2[i]->Draw("e0 same"); m2[i]->Draw("e1 same"); hm2[i]->Draw("c same"); 
    if(i==2) label->DrawLatex(0.01,0.01,sHisto);
  }
  TString epsname = "AWG82/results/keys/eps/EpsKeys";
  epsname += Base; epsname += ".eps";
  all6.SaveAs(epsname);

  return; 
}

void formatHisto(TH1F *h){
  h->SetTitle("");
  h->Sumw2();
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

