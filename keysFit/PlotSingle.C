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
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void formatHisto(TH1F *h);

void PlotSingle(TString hname){
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

  TString Base = "_"; Base += sam; Base += "_"; 
  Base += binString; Base += "_"; Base += smoo; 

  TFile hfile(hname); 
  TH2F *h2 = (TH2F *)gDirectory->Get("h2");
  TString inputfile = "fitSamples/pdfSample"; inputfile += sam; inputfile += ".root";
  cout << "File = " << inputfile << endl;	
  TChain c("ntp1");
  c.Add(inputfile);
  TCanvas mm("mm","KEYS fits to mmiss-pstarl",1200,800);
  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4;
  double xlow = m2min,xhigh = m2max, ylow = plmin,yhigh = plmax;
  Int_t nbinx = 80,nbiny = 80, Sam = sam.Atoi();
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
  if(Sam==0){xlow=-1; xhigh=1;}
  double entries = c.GetEntries();
  Double_t hIntegral=h2->Integral();

  TString M2titles[] = {"0 < p*_{l} < 1 GeV","1 < p*_{l} < 1.4 GeV","1.4 < p*_{l} < 1.8 GeV",
		      "1.8 < p*_{l} < 2.4 GeV","0 < p*_{l} < 2.4 GeV"};
  TString Pltitles[] = {"-4 < m^{2}_{miss} < 1 GeV^{2}","1 < m^{2}_{miss} < 12 GeV^{2}",
			"-4 < m^{2}_{miss} < 12 GeV^{2}"};
  TCut M2cuts[] = {"candPstarLep<1","candPstarLep>1&&candPstarLep<1.4",
		      "candPstarLep>1.4&&candPstarLep<1.8","candPstarLep>1.8&&candPstarLep<2.4", ""};
  TCut Plcuts[] = {"candM2<1","candM2>=1",""};
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

  TString psname = "AWG82/results/keys/eps/Stitch/epsKeys"; psname += Base; psname += ".ps";
  mm.Print(psname+"[");
  TH1F *hm2[5], *m2[5], *hpl[3], *pl[3];
  TString M2names[5], Plnames[3];
  for(int i=0;i<5;i++){
    M2names[i] = "hm2_"; M2names[i] += i;
    hm2[i] = new TH1F(M2names[i],M2titles[i],nM2bin,m2min,m2max); 
    if(i<3) {
      Plnames[i] = "hpl_"; Plnames[i] += i;
      hpl[i] = new TH1F(Plnames[i],Pltitles[i],nPlbin,plmin,plmax); 
    }
  }
  
  hIntegral = h2->Integral();
  for(int i=0;i<5;i++){
    TString hdname = "dm2"; hdname += i;
    TString vari = "candM2>>"; vari+=hdname; vari+="("; vari+= nbinx; vari+=",";vari+= xlow; 
    vari+=",";vari+= xhigh; vari+=")";
    //M2cuts[i] += "weight<100";M2cuts[i] *= "wFF";
    c.Draw(vari,M2cuts[i]);
    m2[i] = (TH1F*)gDirectory->Get(hdname);
    m2[i]->SetXTitle("m^{2}_{miss} [GeV^{2}]");
    formatHisto(m2[i]);
    gStyle->SetOptStat(0);
    if(i<4){
      for(int j=1; j<nM2bin+1; j++){
	double binVal = 0;
	for(int binp = binlim[i]+1; binp < binlim[i+1]+1; binp++){
	  binVal+=h2->GetBinContent(j,binp)*entries*nM2bin*(xhigh-xlow)/nbinx/(m2max-m2min)/hIntegral;
	}
	hm2[i]->SetBinContent(j,binVal);
      }
    } 
    hm2[i]->SetLineColor(4);
    hm2[i]->SetLineWidth(1);
    if(i<4) hm2[4]->Add(hm2[i]);
  }
  for(int i=0;i<3;i++){
    TString hdname = "pl"; hdname += i;
    TString vari = "candPstarLep>>"; vari+=hdname; vari+="("; vari+= nbiny; vari+=",";vari+= ylow; 
    vari+=",";vari+= yhigh; vari+=")";
    cout<<vari<<"  "<<endl;
    //Plcuts[i] += "weight<100";Plcuts[i] *= "wFF";
    c.Draw(vari,Plcuts[i]);
    pl[i] = (TH1F*)gDirectory->Get(hdname);
    pl[i]->SetXTitle("p*_{l} [GeV]");
    formatHisto(pl[i]);
    gStyle->SetOptStat(0);
    if(i<2){
      for(int j=1; j<nPlbin+1; j++){
	double binVal = 0;
	for(int binp = plbinlim[i]+1; binp < plbinlim[i+1]+1; binp++){
	  binVal += h2->GetBinContent(binp,j)*entries*nPlbin/nbiny/hIntegral;	
	}
	hpl[i]->SetBinContent(j,binVal);
      }
    }
    hpl[i]->SetLineColor(4);
    hpl[i]->SetLineWidth(1);
    if(i<2) hpl[2]->Add(hpl[i]);
  }
  cout<<"End for"<<endl;
  m2[4]->Draw("e0"); m2[4]->Draw("e1 same"); hm2[4]->Draw("c same"); mm.Print(psname);
  pl[2]->Draw("e0"); pl[2]->Draw("e1 same"); hpl[2]->Draw("c same"); mm.Print(psname);
  for(int i=0;i<4;i++){
    m2[i]->Draw("e0"); m2[i]->Draw("e1 same"); hm2[i]->Draw("c same"); mm.Print(psname);
  }
  for(int i=0;i<2;i++){
    pl[i]->Draw("e0"); pl[i]->Draw("e1 same"); hpl[i]->Draw("c same"); mm.Print(psname);
  }
  mm.Print(psname+"]");

  TLatex *label = new TLatex(); label->SetNDC(kTRUE); label->SetTextSize(0.055);

  TCanvas all6("all6","KEYS fits to mmiss-pstarl",1700,1800);
  all6.cd();
  all6.Divide(2,3,0.001,0.001);
  all6.cd(1);gPad->SetTopMargin(0.004);gPad->SetRightMargin(0.003);gPad->SetBottomMargin(0.11);
  gPad->SetFillColor(10);gPad->SetFrameFillColor(10);
  m2[4]->Draw("e0");  m2[4]->Draw("e1 same"); hm2[4]->Draw("c same");
  all6.cd(2);gPad->SetTopMargin(0.004);gPad->SetRightMargin(0.003);gPad->SetBottomMargin(0.11);
  gPad->SetFillColor(10);gPad->SetFrameFillColor(10);
  pl[2]->Draw("e0");  pl[2]->Draw("e1 same"); hpl[2]->Draw("c same");
  for(int i=0;i<4;i++){
    all6.cd(i+3);gPad->SetTopMargin(0.003);gPad->SetRightMargin(0.003);gPad->SetBottomMargin(0.11);
    gPad->SetFillColor(10);gPad->SetFrameFillColor(10);
    m2[i]->Draw("e0"); label->DrawLatex(0.69,0.92,M2titles[i]);
    m2[i]->Draw("e1 same"); hm2[i]->Draw("c same"); 
  }
  TString epsname = "AWG82/results/keys/eps/Stitch/EpsKeys";
  epsname += Base; epsname += ".eps";
  all6.SaveAs(epsname);

  h2->Delete();
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
