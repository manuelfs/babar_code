//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: BootstrapBands.cc,v 1.5 2012/08/23 02:22:17 manuelf Exp $
//
// Description:
//      BootstrapBands - Plots the 1 and 2 sigma envelopes of the bootstrapped
//      fits to a sample
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      09/12/18 manuelf -- Created
//------------------------------------------------------------------------

#include "TCanvas.h"
#include "TChain.h"
#include "TCollection.h"
#include "TCut.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TKey.h"
#include "TLatex.h"
#include "TObject.h"
#include "TString.h"
#include "TStyle.h"
#include "DonutUtils/Styles.cc"
#include "DonutUtils/KeysUtils.cc"
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;
using std::cout;
using std::endl;

int main(int argc, char *argv[]){
  
  if (argc < 2 || argc > 8 ) {
    cout << "USAGE: BootstrapBands Filename [nbins=-1] [MinM2=-99] [MaxM2=-99] [maxFactor=1] [nPads=6]"
	 <<" [DataFolder=RAllNewx100]" << endl;
    return 0;
  }
  TString hname = argv[1]; 
  int nbins = -1;;
  if(argc>2) {TString temp_s = argv[2]; nbins = temp_s.Atoi();} 
  double MinM2 = -99.;
  if(argc>3) {TString temp_s = argv[3]; MinM2 = temp_s.Atof();} 
  double MaxM2 = -99.;
  if(argc>4) {TString temp_s = argv[4]; MaxM2 = temp_s.Atof();} 
  double maxFactor = 1.;
  if(argc>5) {TString temp_s = argv[5]; maxFactor = temp_s.Atof();} 
  int nPads = 6;
  if(argc>6) {TString temp_s = argv[6]; nPads = temp_s.Atoi();} 
  TString DataFolder = "RAllNewx100";
  if (argc>7) DataFolder = argv[7];
  TString sam, smoo, scaleP, Base;
  Int_t Sam, nM2bin, nPlbin;
  ParseName(hname, sam, smoo, scaleP, Base, Sam, nM2bin, nPlbin);

  Styles style6; style6.setPadsStyle(nPads); style6.applyStyle();
  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4, hIntegral=-999.;
  //TString singleName = "keys/root/fitNewx100/pdfKeys"; singleName+=Base; singleName+=".root";
  TString singleName = "keys/root/fitNoConv/pdfKeys"; singleName+=Base; singleName+=".root";
  //TString singleName = "keys/root/Fit/pdfKeys"; singleName+=Base; singleName+=".root"; //Prov
  TFile singleFile(singleName); TH2F* hSingle = 0;
  if(!singleFile.IsZombie()) {
    hSingle = (TH2F*)gDirectory->Get("h2");
    hIntegral = hSingle->Integral();
  } else cout<<singleName<<" not found. Using the average of the fits"<<endl;
  nM2bin = hSingle->GetNbinsX(); nPlbin = hSingle->GetNbinsY();   
  TCanvas mm("Bands","Bands");
  double PlLimits[] = {0, 1, 1.4, 1.8, 2.4};
  double M2Limits[] = {0., 5., 16.};
  int binlim[6]; int plbinlim[4];
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
    if(nHisto%50==0) cout<<"Doing "<<nHisto<<" of some number of histograms"<<endl;
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
    if(nHisto==0) hIntegral = h2i->Integral();
    nHisto++;
    //if(nHisto==100) break;
  }
  f.Close();

  TString inputfile = nameData(Sam, DataFolder);
  TChain *treeData = new TChain("ntp1");
  treeData->Add(inputfile);
  double entries = getWentries(treeData);
  if(entries<0.00001) entries = 1;

  Float_t maxM2[]={0,0,0,0,0};
  double xlow = m2min,xhigh = m2max;
  Int_t nbinx = 80,nbiny = 80;
  if(nbins>0){
    nbinx=nbins; nbiny=nbins; xlow=MinM2; xhigh=MaxM2;
  } else setBins(Sam, xlow, xhigh, nbinx, nbiny);

  TString M2titles[] = {"p*_{l} < 1 GeV","1 #leq p*_{l} < 1.4 GeV","1.4 #leq p*_{l} < 1.8 GeV",
		      "1.8 #leq p*_{l} < 2.4 GeV","0 < p*_{l} < 2.4 GeV"};
  TString Pltitles[] = {"-4 < m^{2}_{miss} < 1 GeV^{2}","1 < m^{2}_{miss} < 12 GeV^{2}",
			"-4 < m^{2}_{miss} < 12 GeV^{2}"};
  TCut M2cuts[] = {"candPstarLep<1","candPstarLep>=1&&candPstarLep<1.4",
		      "candPstarLep>=1.4&&candPstarLep<1.8","candPstarLep>=1.8&&candPstarLep<2.4", ""};
  TCut Plcuts[] = {"candM2<1.","candM2>=1.",""};

  TString psname = "AWG82/results/keys/Bands/epsBandKeys"; psname+=Base; psname += ".ps";
  //mm.Print(psname+"[");
  //TString singleName = "AWG82/results/keys/root/single/hKeysND"; singleName+=Base; singleName+=".root";
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
    TString hdname = "m2"; hdname += i;
    m2[i] = new TH1F(hdname, "", nbinx, xlow, xhigh);
    m2[i]->Sumw2();
    TString vari = "candM2>>"; vari+=hdname; 
    TCut totCut = M2cuts[i]; totCut *= "weight";
    treeData->Draw(vari,totCut);
    style6.setMarkers(m2[i], 0.7, 20);
    m2[i]->SetXTitle("m^{2}_{miss} [GeV^{2}]");
    for(int j=1; j<nM2bin+1; j++){
      double mean=0;
      double factor = entries*nM2bin*(xhigh-xlow)/nbinx/(m2max-m2min)/hIntegral;
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
      if(hSingle){
	double binVal = 0;
	if(i==4){binlim[4]=0; binlim[5]=nPlbin;}
	for(int binp = binlim[i]+1; binp < binlim[i+1]+1; binp++){
	  binVal+=hSingle->GetBinContent(j,binp)*entries*nM2bin*(xhigh-xlow)/nbinx/(m2max-m2min)/hIntegral;
	}
	hm2[i]->SetBinContent(j,binVal);
      } else hm2[i]->SetBinContent(j,mean);
      if(maxM2[i]<valsigmas[3])maxM2[i]=valsigmas[3];
    }
    hm2[i]->SetLineColor(4); //hm2[i]->SetLineColor(2); //Prov
    hm2[i]->SetLineWidth(1);
    hm2Error[i]->SetFillColor(38);
    hm2Error2[i]->SetFillColor(9);
    double maxi = m2[i]->GetMaximum();
    if(maxM2[i]>maxi) maxi = maxM2[i];
    if(1 > m2[i]->GetMaximum()) m2[i]->SetMaximum(1.15);
    m2[i]->SetMaximum(maxi*1.1/maxFactor);
  }
  for(int i=2;i<3;i++){
    TString hdname = "pl"; hdname += i;
    pl[i] = new TH1F(hdname, "", nbiny, 0, 2.4);
    pl[i]->Sumw2();
    TString vari = "candPstarLep>>"; vari+=hdname; 
    TCut totCut = Plcuts[i]; totCut *= "weight";
    treeData->Draw(vari,totCut);
    style6.setMarkers(pl[i], 0.7, 20);
    for(int j=1; j<nPlbin+1; j++){
      double mean=0;
      double factor = entries*nPlbin/nbiny/hIntegral;
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
      if(hSingle){
	double binVal = 0;
	if(i==2){plbinlim[i]=0; plbinlim[i+1]=nM2bin;}
	for(int binp = plbinlim[i]+1; binp < plbinlim[i+1]+1; binp++){
	  binVal+=hSingle->GetBinContent(binp,j)*entries*nPlbin/nbiny/hIntegral;
	}
	hpl[i]->SetBinContent(j,binVal);
      } else hpl[i]->SetBinContent(j,mean);
    }
    hpl[i]->SetLineColor(4);
    hpl[i]->SetLineWidth(1);
    hplError[i]->SetFillColor(38);
    hplError2[i]->SetFillColor(9);
  }
  TLatex *label = new TLatex(); label->SetNDC(kTRUE); label->SetTextSize(0.055);
//   m2[4]->Draw("e1"); hm2Error2[4]->Draw("e4 same"); hm2Error[4]->Draw("e4 same"); 
//   m2[4]->Draw("e0 same"); m2[4]->Draw("e1 same"); hm2[4]->Draw("c same"); //mm.Print(psname);
//   pl[2]->Draw("e1"); hplError2[2]->Draw("e4 same"); hplError[2]->Draw("e4 same"); 
//   pl[2]->Draw("e0 same"); pl[2]->Draw("e1 same"); hpl[2]->Draw("c same"); //mm.Print(psname);
//   for(int i=0;i<4;i++){
//     m2[i]->Draw("e1");
//     hm2Error2[i]->Draw("e4 same"); hm2Error[i]->Draw("e4 same"); 
//     m2[i]->Draw("e0 same"); m2[i]->Draw("e1 same"); hm2[i]->Draw("c same"); 
//     //mm.Print(psname);
//   }
//   for(int i=0;i<2;i++){
//     pl[i]->Draw("e1");
//     hplError2[i]->Draw("e4 same"); hplError[i]->Draw("e4 same"); 
//     pl[i]->Draw("e0 same"); pl[i]->Draw("e1 same"); hpl[i]->Draw("c same"); 
//     //mm.Print(psname);
//   }
  //mm.Print(psname+"]");
  
  TString sHisto = ""; sHisto+=nHisto; sHisto += " fits";
  TString Mytitle = "Entries/("; Mytitle += RoundNumber(xhigh-xlow,2,(double)nbinx); Mytitle += " GeV^{2})";
  TString Pytitle = "Entries/("; Pytitle += RoundNumber(2.4*1000,0,(double)nbiny); Pytitle += " MeV)";
  bool doSigma2 = true;
  TCanvas all6("all6","KEYS fits to mmiss-pstarl");
  if(nPads==6){
    all6.Divide(2,3);  all6.cd(1);
  }
  m2[4]->Draw("e0"); if(doSigma2) hm2Error2[4]->Draw("e4 same");  
  hm2Error[4]->Draw("e4 same"); m2[4]->Draw("e0 same"); m2[4]->Draw("e1 same"); hm2[4]->Draw("c same");
  style6.setTitles(m2[4],"m^{2}_{miss} (GeV^{2})",Mytitle);
  if(nPads==6){
    all6.cd(2);
    pl[2]->Draw("e0"); if(doSigma2) hplError2[2]->Draw("e4 same");  
    hplError[2]->Draw("e4 same"); pl[2]->Draw("e0 same"); pl[2]->Draw("e1 same"); hpl[2]->Draw("c same");
    style6.setTitles(pl[2],"p*_{l} (GeV)",Pytitle);
    for(int i=0;i<4;i++){
      all6.cd(i+3);
      m2[i]->Draw("e0");  if(doSigma2) hm2Error2[i]->Draw("e4 same");
      hm2Error[i]->Draw("e4 same"); 
      m2[i]->Draw("e0 same"); m2[i]->Draw("e1 same"); hm2[i]->Draw("c same"); 
      style6.setTitles(m2[i],"m^{2}_{miss} (GeV^{2})",Mytitle,"",M2titles[i]);
      if(i==2 && style6.isThesis=="no") label->DrawLatex(0.01,0.01,sHisto);
    }
  }

  TString epsname = "AWG82/results/keys/eps/Bands/EpsKeys";
  epsname += Base; epsname += ".eps";
  all6.SaveAs(epsname);
  cout<<nHisto<<" histograms used to plot bands on "<<epsname<<endl;

  TString histosName = "keys/eps/CSample/root/EpsKeys"; histosName += Base; histosName += ".root";
  TFile fHistos(histosName,"RECREATE");  fHistos.cd(); 
  m2[4]->Write(); hm2Error2[4]->Write(); hm2Error[4]->Write(); hm2[4]->Write(); 
  pl[2]->Write(); hplError2[2]->Write(); hplError[2]->Write(); hpl[2]->Write(); 
  fHistos.Close();
  cout<<"Histograms saved in "<<histosName<<endl;

  treeData->Delete(); 
  return 1; 
}


