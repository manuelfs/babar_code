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

void THDM_Comp(TString HigTag = "100", int PDFindex=5, int isPLep = 0){

  int nHis = 2, nPads = 3;
  Styles style; style.setPadsStyle(nPads); 
  //style.PadRightMargin = 0.029; style.PadLeftMargin = 0.151; style.yTitleOffset = 1.06; 
  style.applyStyle();
  TCanvas can("can","2HDM PDFs");
  can.Divide(nPads,1); 


  double legW = 0.3, legH = 0.18;
  double legX = 1-style.PadRightMargin-0.01-legW, legY = 1-style.PadTopMargin-0.01;
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(style.LabelSize); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont); leg.SetBorderSize(0);


  TString Name, Variable[] = {"trueQ2", "trueCTL", "candPstarLep"}, legName[4];
  TString xTitle[] = {"q_{true}^{2} (GeV^{2})", "cos(#theta_{#tau,true})", "p*_{l} (GeV)"};
  TString yTitle[] = {"Density/(GeV^{2})","Density","Density/(GeV)"};

  TString Cuts = "weight*(MCTaumode>0&&(MCType=="; Cuts += PDFindex; Cuts += "||MCType=="; Cuts += PDFindex+6; Cuts += "))";
  TString Folder[] = {"AWG82/ntuples/Arxivsmall/FitRAllHigx",  "AWG82/ntuples/small/FitRAllHigx"};
  TString EndName[] = {"_RunAll.root", "_RunAll.root"};
  TString LegName[] = {"ME", "ME + p_{l}"};
  int nBins = 20, color[] = {4,1,28};
  if(HigTag=="000") color[1] = 2;
  //double maxH[] = {0,0,0}, tBmH[4], limX[3][2] = {{3, 12}, {0, 3.2}, {0, 2.2}};
  double maxH[] = {0,0,0}, tBmH[4], limX[3][2] = {{3, 12}, {-1, 1}, {0, 2.2}};
  if(isPLep){
    Variable[2] = "truePLep";
    limX[2][1] = 3;
    xTitle[2] = "p_{l,true} (GeV)";
  }
  
  TH1F *Histo[3][4];
  for(int his=0; his<nHis; his++){
    TString cName = Folder[his]; cName += HigTag; cName += EndName[his];
    //cout<<"Reading "<<cName<<endl;
    TChain chain("ntp1");
    chain.Add(cName);
    tBmH[his] = HigTag.Atof()/100.;
    legName[his] = LegName[his]; 
    for(int pad=0; pad<nPads; pad++){
      can.cd(pad+1);
      Name = "Histo"; Name += pad; Name += his;
      Histo[pad][his] = new TH1F(Name, "",nBins,limX[pad][0],limX[pad][1]);
      TString vari = Variable[pad];
      //if(pad==1 && PDFindex==6) vari = "acos(trueCTL)";
      if(pad==1 && PDFindex==5) vari = "cos(trueCTL)";
      chain.Project(Name, vari, Cuts);
      Histo[pad][his]->SetLineWidth(2);
      Histo[pad][his]->SetLineColor(color[his]);
      Histo[pad][his]->Scale(nBins/Histo[pad][his]->Integral()/(limX[pad][1]-limX[pad][0]));
      if(maxH[pad]<Histo[pad][his]->GetMaximum()) maxH[pad]=Histo[pad][his]->GetMaximum();
      if(pad==0) leg.AddEntry(Histo[pad][his], legName[his]);
    }
  }
    
  for(int pad=0; pad<nPads; pad++){
    can.cd(pad+1);
    for(int his=0; his<nHis; his++){
      if(his==0) {
	Histo[pad][his]->SetMinimum(0);
	Histo[pad][his]->SetMaximum(maxH[pad]*1.05);
	Histo[pad][his]->Draw("");
	style.setTitles(Histo[pad][his],xTitle[pad],yTitle[pad]);
      }else Histo[pad][his]->Draw("same");
    }
    if(pad==2) leg.Draw();
  }
  TString epsName = "public_html/THDM_Comp_"; epsName += HigTag; epsName += "_";
  epsName += PDFindex; epsName += ".eps";
  can.SaveAs(epsName);

  for(int pad=0; pad<nPads; pad++)
    for(int his=0; his<nHis; his++)
      if(Histo[pad][his]) Histo[pad][his]->Delete();
}



