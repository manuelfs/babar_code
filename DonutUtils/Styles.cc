#ifndef STYLES_HH
#define STYLES_HH

//----------------------------------------------------------------------------
// File and Version Information:
//      $Id: Styles.cc,v 1.2 2012/03/04 00:33:03 manuelf Exp $
//
// Description:
//      Styles - Class to set default styles, depending on the number of pads
//  
//  
// 
// Author List:
//      Jochen Dingfelder                         SLAC
//      Wells Wulsin                              Stanford University
//      Manuel Franco Sevilla                     Stanford University
//
// History:
//      10/10/22  manuelf -- Adapted to Donut's needs
//      10/08/09  wulsin  -- Copied from SetStyles.C.  Made into a proper class.  
//----------------------------------------------------------------------------

#if !defined(__CINT__)

#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooStringVar.hh"
#include "TCanvas.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TLatex.h"
#include <iostream>
#include <iomanip>

using std::cin;
using std::cout;
using std::endl;
#endif

class Styles {
public: 
  Styles();  
  void setPadsStyle(int numberPads = 1);
  void testGlobalStyle(bool fixY = true, float scale = 1000.); 
  void setGlobalStyle();
  void applyStyle();
  void printValues();
  void fixYAxis(TH1 *h, TPad *pad);
  void styleHist(TH1 *h, Int_t color = 1, Int_t fillstyle = 0,
		 Int_t symbol = 8,Double_t size = 0.7, Int_t width = 1);
  void setMarkers(TH1 *h, float Msize=0.6, int Mstyle=20) ;
  void setTitles(TH1 *h, TString xTitle="", TString yTitle="", TString Left="", TString Right="");
  void setTitleSizes(TH1 *h,  float size, float lsize, int font=62, 
		     float xoff=1., float yoff=1., int divisions=405);
  int nFont, nPads, nDivisions;
  int CanvasW, CanvasH;
  float TextSize, TitleSize, LabelSize, xTitleOffset, yTitleOffset;
  float PadRightMargin, PadBottomMargin, PadLeftMargin, PadTopMargin;
  TString isThesis;

};

Styles::Styles() {
  RooRealVar   rrvnFont      ("nFont"	   ,"nFont"      ,0);
  RooRealVar   rrvnDivisions ("nDivisions" ,"nDivisions" ,0);
  RooStringVar rrvisThesis   ("isThesis"   ,"nDivisions" ,"yes");
  RooArgSet AllPars(rrvnFont, rrvnDivisions, rrvisThesis);
  AllPars.readFromFile("../DonutUtils/ConfigurationPlots.txt",0,"General");

  nFont      = (int)rrvnFont.getVal();   
  nDivisions = (int)rrvnDivisions.getVal();
  isThesis   = rrvisThesis.getVal();
  nPads      = 1;
}

// Default styles for each possible number of pads
// The values of PadLeftMargin and yTitleOffset are set in between the optimal
// values for 4 and 3 digits in the y axis
void Styles::setPadsStyle(int numberPads) {
  nPads = numberPads;
  TString Section = "Pads_"; Section += nPads;

  RooRealVar rrvCanvasW	       ("CanvasW"	 ,"CanvasW"	    ,0);
  RooRealVar rrvCanvasH	       ("CanvasH"	 ,"CanvasH"	    ,0);
  RooRealVar rrvTextSize       ("TextSize"	 ,"TextSize"	    ,0);          
  RooRealVar rrvTitleSize      ("TitleSize"      ,"TitleSize"	    ,0);  
  RooRealVar rrvLabelSize      ("LabelSize"      ,"LabelSize"	    ,0);   
  RooRealVar rrvPadRightMargin ("PadRightMargin" ,"PadRightMargin"  ,0); 
  RooRealVar rrvPadTopMargin   ("PadTopMargin"   ,"PadTopMargin"    ,0); 
  RooRealVar rrvPadBottomMargin("PadBottomMargin","PadBottomMargin" ,0); 
  RooRealVar rrvxTitleOffset   ("xTitleOffset"   ,"xTitleOffset"    ,0); 
  RooRealVar rrvPadLeftMargin  ("PadLeftMargin"  ,"PadLeftMargin"   ,0); 
  RooRealVar rrvyTitleOffset   ("yTitleOffset"   ,"yTitleOffset"    ,0);  
  
  RooArgSet AllPars(rrvCanvasW, rrvCanvasH, rrvTextSize, rrvTitleSize, rrvLabelSize, 
		    rrvPadRightMargin, rrvPadTopMargin, rrvPadBottomMargin, rrvxTitleOffset);
  AllPars.add(rrvPadLeftMargin);
  AllPars.add(rrvyTitleOffset);  
  AllPars.readFromFile("../DonutUtils/ConfigurationPlots.txt",0,Section, true);

  CanvasW         = (int)rrvCanvasW   .getVal();   
  CanvasH         = (int)rrvCanvasH   .getVal();   
  TextSize        = rrvTextSize	      .getVal();   
  TitleSize       = rrvTitleSize      .getVal();
  LabelSize       = rrvLabelSize      .getVal();
  PadRightMargin  = rrvPadRightMargin .getVal();
  PadTopMargin    = rrvPadTopMargin   .getVal();
  PadBottomMargin = rrvPadBottomMargin.getVal();
  xTitleOffset    = rrvxTitleOffset   .getVal();
  PadLeftMargin   = rrvPadLeftMargin  .getVal();
  yTitleOffset    = rrvyTitleOffset   .getVal();
}

// ----------------------------------------------------------------------
void Styles::fixYAxis(TH1 *h, TPad *pad){
  float maxi = h->GetMaximum()*1.15;
  int digits = (int)(log(maxi)/log(10.)+0.001)+1;
  if(digits<2) digits = 2;

  TString Section = "Pads_"; Section += nPads;
  Section += "_Digits_"; Section += digits;
  RooRealVar rrvPadLeftMargin  ("PadLeftMargin"  ,"PadLeftMargin"   ,0); 
  RooRealVar rrvyTitleOffset   ("yTitleOffset"   ,"yTitleOffset"    ,0);  
  RooArgSet AllPars(rrvPadLeftMargin, rrvyTitleOffset);
  AllPars.readFromFile("../DonutUtils/ConfigurationPlots.txt",0,Section, true);

  PadLeftMargin   = rrvPadLeftMargin  .getVal();
  yTitleOffset    = rrvyTitleOffset   .getVal();

  h->SetTitleOffset(yTitleOffset,"y");
  pad->SetLeftMargin(PadLeftMargin);
}

// Test the global style settings for a generic histogram.  
void Styles::testGlobalStyle(bool fixY, float scale) {
  
  setPadsStyle(nPads); setGlobalStyle(); applyStyle();
  
  TH1* h = new TH1F("h", "h", 50, 0, 50);
  TH1* hc[6];
  for (int i=1; i<=50; i++) {
    double value = scale*exp(-0.5*pow(((i-25.)/5.),2));  // Gaussian shape
    h->SetBinContent(i, value);
  }

  TCanvas c;
  if(nPads == 2) c.Divide(2);
  if(nPads == 3) c.Divide(3);
  if(nPads == 4) c.Divide(2,2);
  if(nPads == 6) c.Divide(3,2);
  TPad *cPad = (TPad *)c.cd(1); h->Draw();
  if(fixY) fixYAxis(h,cPad);
  setTitles(h, "D^{(*)0/+} channels", "xlabel^{2}_{miss} (GeV^{2})", "Events/(10 MeV^{2})");
  float scales[] = {0.1, 10, 0.01};
  for(int pads = 2; pads<=4; pads++){
    if(nPads>=pads){
      cPad = (TPad *)c.cd(pads); 
      hc[pads-2] = (TH1F*)h->Clone();
      hc[pads-2]->Scale(scales[pads-2]); 
      if(fixY) fixYAxis(hc[pads-2],cPad);
      hc[pads-2]->Draw();
      setTitles(hc[pads-2], "D^{(*)0/+} channels", "xlabel^{2}_{miss} (GeV^{2})", "Events/(1000 MeV^{2})");
    }
  }
  TString epsName = "babar_code/Styles/Plot_"; epsName += nPads; epsName += "Pads.eps";
  c.Print(epsName);
  
}

void Styles::applyStyle() {
  setGlobalStyle();
  gStyle->SetCanvasDefW(CanvasW);
  gStyle->SetCanvasDefH(CanvasH);
  gStyle->SetTextSize(TextSize);            // Set global text size
  gStyle->SetTitleSize(TitleSize,"xy");     // Set the 2 axes title size
  gStyle->SetLabelSize(LabelSize,"xy");     // Set the 2 axes label size

  gStyle->SetTitleOffset(xTitleOffset,"x");     
  gStyle->SetTitleOffset(yTitleOffset,"y");     
  gStyle->SetPadRightMargin (PadRightMargin);    
  gStyle->SetPadBottomMargin(PadBottomMargin); 
  gStyle->SetPadTopMargin(PadTopMargin); 
  gStyle->SetPadLeftMargin  (PadLeftMargin); 
  gStyle->SetNdivisions(nDivisions, "xy");   // 5 primary ticks and 4 secondary ticks

  gStyle->SetTitleFont(nFont,"xy");          // Set the all 2 axes title font
  gStyle->SetLabelFont(nFont,"xy");          // Set the all 2 axes label font
  gStyle->SetTextFont(nFont);                // Set global text font
}

// Set default styles globally.   
void Styles::setGlobalStyle() {
  gROOT->SetStyle("Plain");           // The standard style for BaBar plots
  gStyle->SetPalette(1);              // Decent colors for 2D plots
  gStyle->SetOptStat(0);              // No Stats box
  gStyle->SetPadTickX(0);             // No ticks at the right
  gStyle->SetPadTickY(0);             // No ticks at the top
}

void Styles::printValues() {
  cout<<"nFont           = " << nFont           << endl;
  cout<<"nPads           = " << nPads           << endl;
  cout<<"nDivisions      = " << nDivisions      << endl;
  cout<<"CanvasW         = " << CanvasW         << endl;   
  cout<<"CanvasH         = " << CanvasH         << endl;   
  cout<<"TextSize        = " << TextSize        << endl;   
  cout<<"TitleSize       = " << TitleSize       << endl;
  cout<<"LabelSize       = " << LabelSize       << endl;
  cout<<"PadRightMargin  = " << PadRightMargin  << endl;
  cout<<"PadTopMargin    = " << PadTopMargin    << endl;
  cout<<"PadBottomMargin = " << PadBottomMargin << endl;
  cout<<"xTitleOffset    = " << xTitleOffset    << endl;
  cout<<"PadLeftMargin   = " << PadLeftMargin   << endl;
  cout<<"yTitleOffset    = " << yTitleOffset    << endl;

}

// ----------------------------------------------------------------------
void Styles::setMarkers(TH1 *h, float Msize, int Mstyle) {
  h->SetMarkerStyle(Mstyle);
  h->SetMarkerSize(Msize);
}

// ----------------------------------------------------------------------
void Styles::setTitles(TH1 *h, TString xTitle, TString yTitle, TString Left, TString Right) {
  if (0==h) {
    cout << " Histogram not defined" << endl;
  } else {
    h->SetTitle(""); h->SetXTitle(xTitle); h->SetYTitle(yTitle);
    TLatex label; label.SetNDC(kTRUE);
    label.SetTextAlign(13);
    label.DrawLatex(PadLeftMargin+0.04,1-PadTopMargin-0.03,Left);  
    label.SetTextAlign(33);
    label.DrawLatex(1-PadRightMargin-0.02,1-PadTopMargin-0.03,Right);  
 }
}

// ----------------------------------------------------------------------
void Styles::styleHist(TH1 *h, Int_t color, Int_t fillstyle,
		       Int_t symbol, Double_t size, Int_t width) {
  h->SetLineColor(color);   
  h->SetLineWidth(width);
  h->SetMarkerColor(color); 
  h->SetMarkerStyle(symbol);  
  h->SetMarkerSize(size); 
  h->SetStats(kFALSE); 
  h->SetFillStyle(fillstyle); 
  h->SetFillColor(color);
}

// ----------------------------------------------------------------------
void Styles::setTitleSizes(TH1 *h,  float size, float lsize, int font,
			   float xoff, float yoff, int divisions) {
  if (0==h) {
    cout << " Histogram not defined" << endl;
  } else {
    h->SetTitleOffset(xoff, "x");      h->SetTitleOffset(yoff, "y");
    h->SetTitleSize(size,   "x");      h->SetTitleSize(size,   "y");
    h->SetLabelSize(lsize,  "x");      h->SetLabelSize(lsize,  "y");
    h->SetLabelFont(font,   "x");      h->SetLabelFont(font,   "y");
    h->GetXaxis()->SetTitleFont(font); h->GetYaxis()->SetTitleFont(font);
    h->SetNdivisions(divisions,"X");   h->SetNdivisions(divisions,   "Y");
  }
}


#endif	/* STYLES_HH */
