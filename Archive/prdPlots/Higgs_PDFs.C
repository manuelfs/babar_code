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


TH1F SplitBin(TH1F *histo, TString hName);

void Higgs_PDFs(int PDFindex=1, int nBins = 25, TString option = "c"){

  Styles style; style.setPadsStyle(-8); 
  style.CanvasH = 350; style.PadRightMargin += 0.003;
  //style.PadRightMargin = 0.029; style.PadLeftMargin = 0.151; style.yTitleOffset = 1.06; 
  style.applyStyle();
  TCanvas can("can","2HDM PDFs");
  can.Divide(2,1); 


  TString Name, HigTag[] = {"000", "030", "050", "100"}, Variable[] = {"candM2","candPstarLep"};
  TString xTitle[] = {"m", "p"}, yTitle[] = {"d","x"};
  int nHis = 4, nBinPDF[2], color[] = {4,2,8,1}, lineStyle[] = {1,3,2,4};
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
      Histo[pad][his].SetLineWidth(3);
      Histo[pad][his].SetLineStyle(lineStyle[his]);
      Histo[pad][his].SetLineColor(color[his]);
      Histo[pad][his].Scale(10*nBinPDF[pad]/Histo[pad][his].Integral()/(limX[pad][1]-limX[pad][0]));
      if(maxH[pad]<Histo[pad][his].GetMaximum()) maxH[pad]=Histo[pad][his].GetMaximum();
    }
  }
  //maxH[0] = 0.24; maxH[1] = 1.65;
  for(int pad=0; pad<2; pad++){
    can.cd(pad+1);
    for(int his=0; his<nHis; his++){
      if(his==0) {
	Histo[pad][his].SetMaximum(maxH[pad]*1.05);
	Histo[pad][his].SetLabelOffset(0.01,"Y");
	Histo[pad][his].GetXaxis()->CenterTitle(true);
	Histo[pad][his].Draw(option);
	style.setTitles(&Histo[pad][his],xTitle[pad],yTitle[pad]);
      }else {
	TString soption = option; soption += " same";
	Histo[pad][his].Draw(soption);
      }
    }
  }
  TString epsName = "public_html/Higgs_PDFs_"; epsName += PDFindex; epsName += ".eps";
  can.SaveAs(epsName);

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
