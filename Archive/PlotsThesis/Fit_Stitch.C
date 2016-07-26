#include "TPad.h"
#include "TString.h"
#include "TFile.h"
#include "TArrow.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TText.h"
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

void Fit_Stitch(){

  Styles style; style.setPadsStyle(4); style.applyStyle();

  TCanvas can("can","Stitching regions");
  can.Divide(2,2); 

  TPad *cPad = (TPad *)can.cd(4);
  TChain *treeData = new TChain("ntp1");
  treeData->Add("fitSamples/RAll/pdfSample9.root");
  TH2F *h2[3]; double integral[3];
  TString hNames[] = {"_NoStitchNew", "New"};
  for(int i=0; i<2; i++){
    TString name = "babar_code/PlotsThesis/root/pdfKeys_9_Fit"; name += hNames[i]; name += ".root";
    TFile hfile(name); hfile.cd();
    TString h2Name = "h2"; h2Name += i;
    h2[i] = (TH2F *)(gDirectory->Get("h2"))->Clone(h2Name);
    h2[i]->SetDirectory(0);
    hfile.Close();
    integral[i] = h2[i]->Integral();
  }
  int colors[] = {2,4,8};
  double legW = 0.38, legH = 0.3;
  double legX = 1-style.PadRightMargin-0.02, legY = 1-style.PadTopMargin-0.02;
  Int_t nM2bin = h2[0]->GetNbinsX(), nPlbin = h2[0]->GetNbinsY(); 
  int nbins = 50;
  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4, minX = -2.5, maxX = 6;
  double entries = getWentries(treeData);
  if(entries<0.00001) entries = 1;

  TH1F* hData = new TH1F("hData","",nbins,minX,maxX);
  TString vari = "candM2>>hData"; 
  if(treeData) treeData->Draw(vari,"weight");
  style.setMarkers(hData, 0.5, 20);
  hData->SetMinimum(0);
  TString ytitle = "Entries/("; ytitle += RoundNumber((maxX-minX),2,(double)nbins); ytitle += " GeV^{2})"; 
  TString xtitle = "m^{2}_{miss} (GeV^{2})"; 

  TH1F hKeys[3]; TString leg2Names[] = {"No stitching", "With stitching", "#rho_{m} = 1.46"};
  TLegend leg2(legX-legW, legY-legH, legX, legY);
  leg2.SetTextSize(style.LabelSize); leg2.SetFillColor(0); 
  leg2.SetTextFont(style.nFont);  leg2.SetBorderSize(0);
  leg2.AddEntry(hData, "MC");
  double Mins[] = {0.9,minX};
  for(int i=0; i<2; i++){
    double scale = entries*nM2bin*(maxX-minX)/nbins/(m2max-m2min)/integral[i];
    TString hName = "h"; hName += i;
    hKeys[i]= binHisto(h2[i],nM2bin,m2min,m2max,plmin,plmax,hName,"candM2");
    hKeys[i].Scale(scale);
    hKeys[i].SetLineColor(colors[i]);
    hKeys[i].SetLineWidth(2);
    hKeys[i].SetAxisRange(Mins[i],maxX);
    leg2.AddEntry(&hKeys[i], leg2Names[i]);
  }
  double maxi = hKeys[1].GetMaximum()*1.1;
  hData->SetMaximum(maxi);
  hData->Draw("e0");  
  hKeys[1].Draw("c same");
  hData->Draw("e1 same"); 
  style.fixYAxis(hData,cPad); style.applyStyle();
  style.setTitles(hData,xtitle,ytitle,"d)");
  TBox box; box.SetFillStyle(0); box.SetLineColor(1);box.SetFillColor(1);
  box.DrawBox(minX-0.5,-1000,maxX+0.3,800);
  TLine line; line.SetLineColor(1); line.SetLineWidth(2);
  line.DrawLine(1.1,0,1.1,maxi);
  TLatex label; label.SetTextAlign(12);
  TArrow arrow; arrow.SetFillColor(1);
  label.DrawLatex(1.5,10000,"m_{s} #pm 0.2 GeV^{2}");
  arrow.DrawArrow(0.9,10000.,1.3,10000.,0.0075,"<|>");

  cPad = (TPad *)can.cd(3);
  TH1F *hData2 = (TH1F*)hData->Clone("hData2");
  TH1F* hKeys2[2];
  for(int i=0; i<2; i++) hKeys2[i] = (TH1F*)hKeys[i].Clone();  
  hData2->SetMaximum(hData->GetMaximum()/212.);
  hData2->Draw("e0");  
  for(int i=0; i<2; i++) hKeys2[i]->Draw("c same");
  hData2->Draw("e1 same"); 
  leg2.Draw();
  style.fixYAxis(hData2,cPad); style.applyStyle();
  style.setTitles(hData,xtitle,ytitle,"c)");
  line.DrawLine(1.1,0,1.1,hData->GetMaximum()/212.);

  cPad = (TPad *)can.cd(2);
  TH1F* hData3 = new TH1F("hData3","",nbins,plmin,plmax);
  vari = "candPstarLep>>hData3"; 
  if(treeData) treeData->Draw(vari,"weight");
  style.setMarkers(hData3, 0.5, 20);
  hData->SetMinimum(0);
  ytitle = "Entries/("; ytitle += RoundNumber((plmax-plmin)*1000,0,(double)nbins); ytitle += " MeV)"; 
  xtitle = "p*_{l} (GeV)"; 
  TH1F hPl[3]; 
  double PMins[] = {1.94,0};
  for(int i=0; i<2; i++){
    double scale = entries*nPlbin/nbins/integral[i];
    TString hName = "hPl"; hName += i;
    hPl[i]= binHisto(h2[i],nPlbin,plmin,plmax,m2min,m2max,hName,"candPstarLep");
    hPl[i].Scale(scale);
    hPl[i].SetLineColor(colors[i]);
    hPl[i].SetLineWidth(2);
    hPl[i].SetAxisRange(PMins[i],2.4);
  }
  hData3->SetMinimum(0);
  hData3->SetMaximum(hData3->GetMaximum()*1.1);
  hData3->Draw("e0");  
  for(int i=0; i<2; i++) hPl[i].Draw("c same");
  hData3->Draw("e1 same"); 
  style.fixYAxis(hData3,cPad); style.applyStyle();
  style.setTitles(hData3,xtitle,ytitle,"b)");
  line.SetLineStyle(2); 
  line.DrawLine(0.8,0,0.8,hData3->GetMaximum());
  line.SetLineStyle(1); 
  line.DrawLine(2,0,2,hData3->GetMaximum());
  label.SetTextAlign(32);label.DrawLatex(0.72,650,"p_{s}^{l}");
  label.SetTextAlign(32);label.DrawLatex(0.72,560,"#pm 60 MeV");
  arrow.DrawArrow(0.74,600.,0.86,600.,0.008,"<|>");
  label.DrawLatex(1.92,95,"p_{s}^{h} #pm 60 MeV");
  arrow.DrawArrow(1.94,100.,2.06,100.,0.008,"<|>");

  cPad = (TPad *)can.cd(1);
  gStyle->SetPalette(1);
  TH2F *hRegions = new TH2F("hRegions","",50,minX,maxX,50,plmin,plmax);
  treeData->Draw("candPstarLep:candM2>>hRegions","","cont4");
  hRegions->SetTitleSize(style.TitleSize,"xyz");
  hRegions->SetLabelSize(style.LabelSize,"xyz");
  hRegions->SetTitleOffset(style.xTitleOffset,"x");
  hRegions->SetTitleOffset(0.7,"y");
  hRegions->GetXaxis()->SetTitleFont(style.nFont);
  hRegions->SetLabelFont(style.nFont,"xyz");
  hRegions->SetNdivisions(style.nDivisions, "xyz");
  hRegions->SetXTitle("m^{2}_{miss} (GeV^{2})");
  hRegions->SetYTitle("p*_{l} (GeV)");
  cPad->SetLeftMargin(0.12);
  double hTop = cPad->GetBottomMargin()+(1-cPad->GetTopMargin()-cPad->GetBottomMargin())*2/2.4;
  double hBot = cPad->GetBottomMargin()+(1-cPad->GetTopMargin()-cPad->GetBottomMargin())*0.8/2.4;
  double hLef = cPad->GetLeftMargin()+(1-cPad->GetRightMargin()-cPad->GetLeftMargin())*(1.1-minX)/(maxX-minX);
  line.DrawLineNDC(cPad->GetLeftMargin(),hTop,1-cPad->GetRightMargin(),hTop);
  line.DrawLineNDC(hLef,hTop,hLef,hBot);
  line.SetLineStyle(2); 
  line.DrawLineNDC(cPad->GetLeftMargin(),hBot,1-cPad->GetRightMargin(),hBot);
  label.SetNDC(kTRUE); label.SetTextAlign(13); label.SetTextColor(0);
  label.DrawLatex(style.PadLeftMargin+0.04,1-style.PadTopMargin-0.03,"a)");  
  label.DrawLatex(0.14,0.65,"Region 1");
  label.DrawLatex(0.72,0.65,"Region 2");
  label.DrawLatex(0.72,0.928,"Region 3");
  label.DrawLatex(0.72,0.33,"Region 4");

  can.cd(0);
  arrow.DrawArrow(0.56,0.09,0.5,0.13,0.01,"|>");
  label.SetTextSize(0.032);
  label.SetTextColor(1); label.DrawLatex(0.525,0.145,"x200");  

  TString pName = "public_html/Fit_Stitch.eps"; 
  can.SaveAs(pName);

  for(int i=0; i<2; i++) {
    h2[i]->Delete();
    hKeys2[i]->Delete();
  }
  hData->Delete(); hData2->Delete(); hData3->Delete(); treeData->Delete(); hRegions->Delete();
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




