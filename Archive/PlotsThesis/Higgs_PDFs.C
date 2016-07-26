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

TH1F binHisto(TH2F *h2, int nbins, double xlow, double xhigh, double miny, double maxy, 
	      TString hname, TString var = "candM2");

TH1F SplitBin(TH1F *histo, TString hName);

void Higgs_PDFs(int PDFindex=1, int nBins = 50, TString option = "c"){

  Styles style; style.setPadsStyle(2); 
  style.PadRightMargin = 0.029; style.PadLeftMargin = 0.151; style.yTitleOffset = 1.06; 
  style.applyStyle();
  TCanvas can("can","2HDM PDFs");
  can.Divide(2,1); 


  double legW = 0.47, legH = 0.3;
  double legX = 1-style.PadRightMargin-0.01-legW, legY = 1-style.PadTopMargin-0.01;
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(style.LabelSize); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont); leg.SetBorderSize(0);


  TString Name, HigTag[] = {"000", "030", "050", "100"}, Variable[] = {"candM2","candPstarLep"}, legName[4];
  TString xTitle[] = {"m_{miss}^{2} (GeV^{2})", "p"}, yTitle[] = {"Density/(GeV^{2})","Density/(GeV)"};
  int nHis = 4, nBinPDF[2], color[] = {2,4,8,1};
  double maxH[] = {0,0}, tBmH[4], limX[2][2] = {{-1,10},{0,2.2}};
  TH2F *pdf[4];
  TH1F Histo[2][4];
  for(int his=0; his<nHis; his++){
    Name = "keys/root/fitHigx"; Name += HigTag[his]; Name += "/pdfKeys_";
    Name += PDFindex; Name += "_Fit.root";
    TFile hfile(Name); 
    TString pdfName = "pdf"; pdfName += his;
    pdf[his] = (TH2F *)(hfile.Get("h2"))->Clone(pdfName);    
    pdf[his]->SetDirectory(0);
    if(his==0) {nBinPDF[0] = pdf[his]->GetNbinsX(); nBinPDF[1] = pdf[his]->GetNbinsY();}
    tBmH[his] = HigTag[his].Atof()/100.;
    legName[his] = "t#beta/m_{H} = "; legName[his] += RoundNumber(tBmH[his],1);
    legName[his] += " GeV^{-1}";
    if(fabs(tBmH[his])<1e-6) legName[his] = "SM"; 
  }
  nBinPDF[0] = nBins; nBinPDF[1] = nBins;
  for(int pad=0; pad<2; pad++){
    can.cd(pad+1);
    for(int his=0; his<nHis; his++){
      Name = "Histo"; Name += pad; Name += his;
      Histo[pad][his] = binHisto(pdf[his],nBinPDF[pad],limX[pad][0],limX[pad][1],
				 limX[(pad+1)%2][0],limX[(pad+1)%2][1],Name,Variable[pad]);
//       Name += "Split";
//       Histo[pad][his] = SplitBin(&Histo[pad][his], Name);
//       Name += "Split";
//       Histo[pad][his] = SplitBin(&Histo[pad][his], Name);
      Histo[pad][his].SetLineWidth(2);
      Histo[pad][his].SetLineColor(color[his]);
      Histo[pad][his].Scale(nBinPDF[pad]/Histo[pad][his].Integral()/(limX[pad][1]-limX[pad][0]));
      if(maxH[pad]<Histo[pad][his].GetMaximum()) maxH[pad]=Histo[pad][his].GetMaximum();
      if(pad==0) leg.AddEntry(&Histo[pad][his], legName[his]);
    }
  }
  //maxH[0] = 0.24; maxH[1] = 1.65;
  for(int pad=0; pad<2; pad++){
    can.cd(pad+1);
    for(int his=0; his<nHis; his++){
      if(his==0) {
	Histo[pad][his].SetMaximum(maxH[pad]*1.05);
	Histo[pad][his].Draw(option);
	style.setTitles(&Histo[pad][his],xTitle[pad],yTitle[pad]);
      }else {
	TString soption = option; soption += " same";
	Histo[pad][his].Draw(soption);
      }
    }
    if(pad==1) leg.Draw();
  }
  TString epsName = "public_html/Higgs_PDFs_"; epsName += PDFindex; epsName += ".eps";
  can.SaveAs(epsName);

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

  TH1F h1(hname,"",nbins, xlow, xhigh);
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


TH1F SplitBin(TH1F *histo, TString hName){
  int nBins = histo->GetNbinsX();
  double minX = histo->GetXaxis()->GetBinLowEdge(1), maxX = histo->GetXaxis()->GetBinLowEdge(nBins+1);
  double prevBin = 0, thisBin = 0;

  TH1F hResult(hName,"",nBins*2,minX,maxX);
  for(int bin=1; bin<=nBins; bin++){
    thisBin = histo->GetBinContent(bin);
    if(bin>1) hResult.SetBinContent(2*bin-2,prevBin*2/3. + thisBin/3.);
    hResult.SetBinContent(2*bin-1,prevBin/3. + thisBin*2/3.);
    prevBin = thisBin;
  }
  hResult.SetBinContent(2*nBins, thisBin*2/3.);

  return hResult;
}
