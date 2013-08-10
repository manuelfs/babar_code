#include "TPad.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "babar_code/Styles/Styles.cc"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace TMath;
using namespace std;
using std::cout;
using std::endl;

int minEvents = 8, minBins = 4;
double calcChi(TH1F *hData, TH1F *hPdf, int &ndof);
double FillChi(TH1F *hChi, TString Variable, TString sample, TString samplePDF, int nbins);
TH1F Histos(TH1F **hData, int pdf, TString Variable, TString sample, TString samplePDF, int nbins);
TH1F binHisto(TH2F *h2, int nbins, double xlow, double xhigh, double miny, double maxy, 
	      TString hname, TString var = "candM2");
double getWentries(TTree *tree, TString Wname="weight");
void setBins(int Sam, double &xlow, double &xhigh);

void Fit_Chi2(TString sample, int Bin1 = 15, int Bin2 = 90, int plotPDF = 18){

  Styles style; style.setPadsStyle(4); 
  style.PadRightMargin = 0.029; style.PadLeftMargin = 0.151; style.yTitleOffset = 1.06; 
  style.applyStyle();
  TCanvas can("can","Chi squares of the PDFs");
  can.Divide(2,2); 

  double legW = 0.42, legH = 0.22;
  double legX = 1-style.PadRightMargin-0.02, legY = 1-style.PadTopMargin-0.02;
  if(sample=="R26") {legX = style.PadLeftMargin+legW+0.02; legY = 1-style.PadTopMargin-0.09;}
  TLegend *leg[2];
  for(int iLeg=0; iLeg<2; iLeg++){
    leg[iLeg] = new TLegend(legX-legW, legY-legH, legX, legY);
    leg[iLeg]->SetTextSize(style.LabelSize); leg[iLeg]->SetFillStyle(0); 
    leg[iLeg]->SetTextFont(style.nFont);  leg[iLeg]->SetBorderSize(0);
  }
  TH1F hPdf[4]; TH1F *hData[4], *hChi[4][2];
  TString lTag[] = {"a)", "b)", "c)", "d)"};
  TString sTag[] = {sample,"R24"};
  TString vTag[] = {"candM2","candPstarLep"};
  TString legTag[] = {"m^{2}_{miss} (#bar{p} = ", "p*_{l} (#bar{p} = "};
  int Bins[] = {Bin1, Bin2}, colors[] = {4,2};
  double minX = 0, maxX = 2.4;
  setBins(plotPDF, minX,maxX);
  for(int pad=0; pad<2; pad++){
    TPad *cPad = (TPad *)can.cd(pad+1+(pad>1)*2);
    hPdf[pad] = Histos(&hData[pad], plotPDF, "candM2", sTag[pad>1], "R26", Bins[pad%2]);
    style.setMarkers(hData[pad], 0.5, 20);
    hData[pad]->Draw("e0");  
    hData[pad]->Draw("e1 same"); 
    hPdf[pad].Draw("same");
    int ndof = -1;
    double chi2 = calcChi(hData[pad], &hPdf[pad], ndof);
    TString pLabel = "#chi^{2} = "; pLabel += RoundNumber(chi2,1); pLabel += "/";
    pLabel += ndof; pLabel += "; p = "; pLabel += RoundNumber(Prob(chi2,ndof)*100,1); pLabel += "%";    
    TString yTitle = "Entries/("; yTitle += RoundNumber(maxX-minX,2,Bins[pad%2]); yTitle += " GeV^{2})";
    style.setTitles(hData[pad], "m^{2}_{miss} (GeV^{2})",yTitle,lTag[pad],pLabel);

    cPad = (TPad *)can.cd(pad+3+(pad>1)*2);
    for(int iVari = 0; iVari<2; iVari++){
      cout<<endl<<"Doing "<<Bins[pad%2]<<" bins and variable "<< vTag[iVari]<<endl;
      TString nameChi = "Chi"; nameChi += pad; nameChi += iVari;
      hChi[pad][iVari] = new TH1F(nameChi,"",10,0,100);
      double meanChi = FillChi(hChi[pad][iVari], vTag[iVari], sTag[pad>1], "R26", Bins[pad%2]);
      hChi[pad][iVari]->SetLineColor(colors[iVari]);
      TString legLabel = legTag[iVari]; legLabel += RoundNumber(meanChi,0); legLabel += "%)";
      leg[pad]->AddEntry(hChi[pad][iVari],legLabel);
      if(iVari==1){
	double maxi = hChi[pad][0]->GetMaximum();
	if(hChi[pad][1]->GetMaximum()>maxi) maxi = hChi[pad][1]->GetMaximum();
	if(sTag[pad>1]=="R24") maxi *= 1.33;
	hChi[pad][0]->SetMaximum(maxi*1.08);
	hChi[pad][0]->SetLineColor(4);
	hChi[pad][0]->Draw();
	style.setTitles(hChi[pad][0], "p-Value (%)","Entries/(10%)",lTag[pad+2]);
	hChi[pad][1]->SetLineColor(2);
	hChi[pad][1]->Draw("same");
	leg[pad]->Draw();
      }
    }  
  }

  TString pName = "public_html/Fit_Chi2"; pName += sample; pName += ".eps";
  can.SaveAs(pName);
  for(int pad=0; pad<2; pad++){
    hData[pad]->Delete();
    for(int iVari = 0; iVari<2; iVari++) hChi[pad][iVari]->Delete();
  }
}

double calcChi(TH1F *hData, TH1F *hPdf, int &ndof){
  double chi2 = 0; ndof = -1;
  for(int bin=1; bin<=hData->GetNbinsX(); bin++){
      if(hData->GetBinContent(bin)<minEvents) continue;
      double diff = hData->GetBinContent(bin)-hPdf->GetBinContent(bin);
      double uncert = hData->GetBinError(bin);
      if(uncert>0){
	chi2 += pow(diff/uncert,2);
	ndof++;
      }
    }
  return chi2;
}

double FillChi(TH1F *hChi, TString Variable, TString sample, TString samplePDF, int nbins){
  double meanChi = 0; int nChi = 0;
  for(int pdf=1; pdf<=68; pdf++){
    if(pdf==33) pdf = 41; if(pdf==49) pdf = 51; if(pdf==59) pdf = 61; 
    if(sample=="R24" && pdf>60) continue;
    TH1F* hData;
    TH1F hPdf = Histos(&hData, pdf, Variable, sample, samplePDF, nbins);
    int ndof = -1;
    double chi2 = calcChi(hData, &hPdf, ndof);
    TString pLabel = "PDF "; pLabel += pdf; pLabel += ": chi^2 = "; pLabel += RoundNumber(chi2,1); pLabel += "/";
    pLabel += ndof; pLabel += "; p = "; pLabel += RoundNumber(Prob(chi2,ndof)*100,1); pLabel += "%";    
    if(ndof>=(minBins-1)) {
      double p = Prob(chi2,ndof)*100;
      hChi->Fill(p);
      if(p<=10) cout<<pLabel<<endl;
      meanChi += p; nChi++;
    }
    hData->Delete();
  }
  hChi->SetMinimum(0); hChi->SetLineWidth(2); 
  return meanChi/(double)nChi;
}

TH1F Histos(TH1F **hData, int pdf, TString Variable, TString sample, TString samplePDF, int nbins){
  double minX = 0, maxX = 2.4, minOther = -4, maxOther = 12;
  if(Variable=="candM2") {
    minOther = 0; maxOther = 2.4;
    setBins(pdf, minX,maxX);
  }
  TString nameData = "fitSamples/"; nameData += sample; nameData += "/pdfSample";
  nameData += pdf; nameData += ".root";
  TChain treeData("ntp1");
  treeData.Add(nameData);
  double entries = getWentries(&treeData);
  if(entries<0.00001) entries = 1;
  TString namePdf = "keys/Archive/11-04-19/root/FitCV_"; namePdf += samplePDF;
  namePdf += "/pdfKeys_"; namePdf += pdf; namePdf += "_Fit.root";
  TFile fPdf(namePdf); fPdf.cd();
  TH2F *hPdf2 = (TH2F *)(gDirectory->Get("h2"))->Clone("hPdf2");
  hPdf2->SetDirectory(0);
  fPdf.Close();
  double scale = entries/hPdf2->Integral();
  TString nData = "hData"; nData+=pdf; nData+=sample; nData+= samplePDF; nData+=nbins; nData+=(int)&nbins;
  *hData = new TH1F(nData,"",nbins,minX,maxX);
  (*hData)->Sumw2();
  (*hData)->SetMinimum(0);
  TString vari = Variable; vari += ">>"; vari += nData; 
  treeData.Draw(vari,"weight");
  TString nPdf = "hPdf"; nPdf += nData;
  TH1F hPdf = binHisto(hPdf2,nbins,minX,maxX,minOther,maxOther,nPdf,Variable);
  hPdf.SetLineColor(4); hPdf.SetLineWidth(2); hPdf.Scale(scale);

  hPdf2->Delete(); 
  return hPdf;
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

void setBins(int Sam, double &xlow, double &xhigh){
  xlow = -4; xhigh = 12;
  if (Sam>=1 && Sam<=8) {
    xlow = -1; xhigh = 10; 
  }
  if (Sam>=9 && Sam<=12) {
    xlow = -0.4; xhigh = 0.5; 
  }
  if(Sam==13 || Sam==15){
    xlow = -0.6; xhigh = 3;
  }
  if(Sam==14 || Sam==16){
    xlow = -2; xhigh = 1.5;
  }
  if(Sam>=17 && Sam <=20){
    xlow = -1; xhigh = 10;
  }
  if(Sam>=21 && Sam <=24){
    xlow = -1.; xhigh = 6;
  }
  if (Sam>=25 && Sam<=28) {
    xhigh = 7; 
  }
  if (Sam>=29 && Sam<=36) {
    xhigh = 7; 
  }
  if (Sam >= 41 && Sam<=44) {
    xlow = -1.5; xhigh = 10;
  }
  if (Sam >= 45 && Sam<=50) {
    xlow = -3; xhigh = 10;
  }
  if (Sam >= 51 && Sam<=60) {
    xlow = -1.5; xhigh = 10;
  }
  if (Sam >= 61 && Sam<=70) {
    xlow = -1.5; xhigh = 10;
  }
  return;
}




