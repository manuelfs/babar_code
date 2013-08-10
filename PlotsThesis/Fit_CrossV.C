#include "TPad.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "babar_code/Styles/Styles.cc"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

TH1F binHisto(TH2F *h2, int nbins, double xlow, double xhigh, double miny, double maxy, 
	      TString hname, TString var = "candM2");
double getWentries(TTree *tree, TString Wname="weight");

//bsub -q xlong -o AWG82/logfiles/results/keys/CrossV/logCrossV_29.log CrossValidation 29 50 25 250 6 65 240 noRot RAll
//bsub -q long -o AWG82/logfiles/results/keys/FitPdf/logFit_29.log FitPdfSampleKEYS 29 118 114 1000 300 100 30 0.3 50 50 2.0 noRot RAll 150 50 1.6
void Fit_CrossV(){

  Styles style; style.setPadsStyle(2); style.applyStyle();
  double LeftMargin = style.PadLeftMargin;
  style.PadLeftMargin *= 1.1;

  TCanvas can("can","Cross Validation values");
  can.Divide(2,1); TPad *cPad = (TPad *)can.cd(1);
  TFile crossFile("babar_code/PlotsThesis/root/hCrossKeys_29_25to250_65to240.root");crossFile.cd();
  TH2F *hCross = (TH2F *)(gDirectory->Get("cross"))->Clone("hCross");
  hCross->SetDirectory(0);
  crossFile.Close();
  TH1D *h1Cross[3];
  int binP[] = {1,3,6}, colors[] = {2,4,8};
  double minC = 1e9, maxC = -1e9;
  double legW = 0.24, legH = 0.28;
  double legX = 1-style.PadRightMargin-0.02, legY = 1-style.PadTopMargin-0.02;
  TLegend leg(legX-legW, legY-legH, legX, legY);
  leg.SetTextSize(style.LabelSize); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont);  leg.SetBorderSize(0);
  TString legNames[] = {"#rho_{p} = 0.65", "#rho_{p} = 1.35", "#rho_{p} = 2.40"};
  for(int i=0; i<3; i++){
    TString h1Name = "h1"; h1Name += i;
    h1Cross[i] = (TH1D*)hCross->ProjectionX(h1Name, binP[i],binP[i]);
    h1Cross[i]->SetTitleSize(style.TitleSize,"xy");
    h1Cross[i]->SetLabelSize(style.LabelSize,"xy");
    h1Cross[i]->SetTitleOffset(style.xTitleOffset,"x");
    h1Cross[i]->SetTitleOffset(style.yTitleOffset*1.1,"y");
    cPad->SetLeftMargin(style.PadLeftMargin);
    h1Cross[i]->GetXaxis()->SetTitleFont(style.nFont);
    h1Cross[i]->SetLabelFont(style.nFont,"xy");
    h1Cross[i]->SetNdivisions(style.nDivisions, "xy");
    h1Cross[i]->SetLineColor(colors[i]);
    h1Cross[i]->SetLineWidth(2);
    if(h1Cross[i]->GetMinimum()<minC) minC = h1Cross[i]->GetMinimum();
    if(h1Cross[i]->GetMaximum()>maxC) maxC = h1Cross[i]->GetMaximum();
    if(i==0) h1Cross[i]->Draw("c");
    else h1Cross[i]->Draw("c same");
    leg.AddEntry(h1Name, legNames[i]);
  }
  leg.Draw();
  h1Cross[0]->SetMaximum(maxC+0.1*(maxC-minC));
  h1Cross[0]->SetMinimum(minC-0.1*(maxC-minC));
  style.setTitles(h1Cross[0],"#rho_{m}","CV output","a)");

  style.PadLeftMargin = LeftMargin;
  can.cd(2);
  TChain *treeData = new TChain("ntp1");
  treeData->Add("fitSamples/RAll/pdfSample29.root");
  TH2F *h2[3]; double integral[3];
  TString hNames[] = {"40", "120", "200"};
  for(int i=0; i<3; i++){
    TString name = "babar_code/PlotsThesis/root/pdfKeys_29_Fit_"; name += hNames[i]; name += ".root";
    TFile hfile(name); hfile.cd();
    TString h2Name = "h2"; h2Name += i;
    h2[i] = (TH2F *)(gDirectory->Get("h2"))->Clone(h2Name);
    h2[i]->SetDirectory(0);
    hfile.Close();
    integral[i] = h2[i]->Integral();
  }
  Int_t nM2bin = h2[0]->GetNbinsX(); 
  int nbins = 35;
  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4;
  double entries = getWentries(treeData), minX = -4, maxX = 7;
  if(entries<0.00001) entries = 1;

  TH1F* hData = new TH1F("hData","",nbins,minX,maxX);
  TString vari = "candM2>>hData"; 
  if(treeData) treeData->Draw(vari,"weight");
  style.setMarkers(hData, 0.5, 20);
  hData->SetMinimum(0);
  TString ytitle = "Entries/("; ytitle += RoundNumber((maxX-minX),2,(double)nbins); ytitle += " GeV^{2})"; 
  TString xtitle = "m^{2}_{miss} (GeV^{2})"; 

  TH1F hKeys[3]; TString leg2Names[] = {"#rho_{m} = 0.4", "#rho_{m} = 1.2", "#rho_{m} = 2.0"};
  TLegend leg2(legX-legW, legY-legH, legX, legY);
  leg2.SetTextSize(style.LabelSize); leg2.SetFillColor(0); 
  leg2.SetTextFont(style.nFont);  leg2.SetBorderSize(0);
  for(int i=0; i<3; i++){
    double scale = entries*nM2bin*(maxX-minX)/nbins/(m2max-m2min)/integral[i];
    TString hName = "h"; hName += i;
    hKeys[i]= binHisto(h2[i],nM2bin,m2min,m2max,plmin,plmax,hName,"candM2");
    hKeys[i].Scale(scale);
    hKeys[i].SetLineColor(colors[i]);
    hKeys[i].SetLineWidth(2);
    hKeys[i].SetAxisRange(minX,maxX);
    leg2.AddEntry(&hKeys[i], leg2Names[i]);
  }
  if(treeData){
    hData->Draw("e0");  
    for(int i=0; i<3; i++) hKeys[i].Draw("c same");
    style.setTitles(hData,xtitle,ytitle,"b)");
    hData->Draw("e1 same"); 
    leg2.Draw();
  }
  TString pName = "public_html/Fit_CrossV.eps"; 
  can.SaveAs(pName);

  for(int i=0; i<3; i++) {
    h1Cross[i]->Delete();
    h2[i]->Delete();
  }
  hData->Delete(); treeData->Delete();
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
  double h2width = (h2xmax-h2xmin)/(double)nh2xbin;
  double h1width = (xhigh-xlow)/(double)nbins;
  double h2yini = (miny-h2ymin)*(double)nh2ybin/(h2ymax-h2ymin)+1;
  double h2yfin = (maxy-h2ymin)*(double)nh2ybin/(h2ymax-h2ymin)+1;
  int h2yini_i = (int)h2yini, h2yfin_i = (int)h2yfin; 
  double iniyfrac = h2yini_i+1 - h2yini, finyfrac = h2yfin - h2yfin_i;
  if(h2yini_i<1) h2yini_i = 1;	        
  if(h2yfin_i>nh2ybin) h2yfin_i = nh2ybin;

  TH1F h1(hname,hname,nbins, xlow, xhigh);
  for(int bin = 1; bin <= nbins; bin++){
    double binvar = 0, tempvar = 0;
    double binmin = xlow+(double)(bin-1)*h1width;
    double h2xini = (binmin-h2xmin)/h2width+1, h2xfin = (binmin+h1width-h2xmin)/h2width+1;
    int h2xini_i = (int)h2xini, h2xfin_i = (int)h2xfin; 
    double inifrac = h2xini_i+1 - h2xini, finfrac = h2xfin - h2xfin_i;
    if(h2xini_i<1) h2xini_i = 1;	        
    if(h2xfin_i>nh2xbin) h2xfin_i = nh2xbin;
    if(h2xini_i==h2xfin_i){inifrac = h1width/h2width; finfrac = 1.;}
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




