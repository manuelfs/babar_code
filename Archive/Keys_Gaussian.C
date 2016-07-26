#include "TPad.h"
#include "TString.h"
#include "TF1.h"
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

#define nPoints 4
using namespace std;
using std::cout;
using std::endl;

void Keys_Gaussian(double rms=0.5){

  Styles style; style.setPadsStyle(1); style.applyStyle();
  TCanvas can("can","Keys",500,400);

  double val[2][4] = {{1.5, 5.5, 6.5, 7.5},{1, 2, 4, 3}};
  double minX = 0, maxX = 10;
  int nBins = 10;
  TString func = "", gau = "*TMath::Gaus(x,";

  TH1F histo("histo","",nBins,minX,maxX);
  for(int poi=0; poi<nPoints; poi++){
    histo.Fill(val[0][poi], val[1][poi]);
    func += "+"; func += val[1][poi]; func += gau; func += val[0][poi]; 
    func += ","; func += rms;func += ",1)"; 
  }
  histo.Sumw2();
  style.setMarkers(&histo, 1.2, 20);
  style.setTitles(&histo,"","");
  TF1 Function("f",func,0,10);
  Function.SetLineColor(2);
  histo.Draw("e1");
  Function.Draw("same");

  TString pName = "public_html/Keys_Gaussian.eps"; 
  can.SaveAs(pName);
}
