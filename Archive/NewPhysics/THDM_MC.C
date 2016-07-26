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

void THDM_MC(int PDFindex=5, int isPLep = 1){

  Styles style; style.setPadsStyle(2); 
  style.PadRightMargin = 0.029; style.PadLeftMargin = 0.151; style.yTitleOffset = 1.06; 
  style.applyStyle();
  TCanvas can("can","2HDM PDFs");
  can.Divide(2,1); 


  double legW = 0.32, legH = 0.27;
  double legX = 1-style.PadRightMargin-0.01-legW, legY = 1-style.PadTopMargin-0.01;
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(style.LabelSize); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont); leg.SetBorderSize(0);


  TString Name, HigTag[] = {"000", "030", "050", "100"}, Variable[] = {"candM2Tru","truePstarLep"}, legName[4];
  TString xTitle[] = {"m_{miss,true}^{2} (GeV^{2})", "p*_{l,true} (GeV)"}, yTitle[] = {"Density/(GeV^{2})","Density/(GeV)"};
  TString Cuts = "MCType=="; Cuts += PDFindex;
  int nHis = 4, nBinPDF[] = {100,100}, color[] = {2,4,28,1};
  double maxH[] = {0,0}, tBmH[4], limX[2][2] = {{-1,10},{0,2.2}};
  if(isPLep){
    Variable[1] = "truePLep";
    limX[1][1] = 3;
    xTitle[1] = "p_{l,true} (GeV)";
  }
  
  TH1F *Histo[2][4];
  for(int his=0; his<nHis; his++){
    TString cName = "AWG82/generator/4M"; cName += HigTag[his]; cName += "Dxtaunu.root";
    //TString cName = "AWG82/generator/Archive/4MEvtGen_Hig"; cName += HigTag[his]; cName += "Dxtaunu.root";
    TChain chain("GqaMCAnalysis/ntp1");
    chain.Add(cName);
    tBmH[his] = HigTag[his].Atof()/100.;
    legName[his] = "t#beta/m_{H} = "; legName[his] += RoundNumber(tBmH[his],1);
    if(fabs(tBmH[his])<1e-6) legName[his] = "SM"; 
    for(int pad=0; pad<2; pad++){
      can.cd(pad+1);
      Name = "Histo"; Name += pad; Name += his;
      Histo[pad][his] = new TH1F(Name, "",nBinPDF[pad],limX[pad][0],limX[pad][1]);
      chain.Project(Name, Variable[pad], Cuts);
      Histo[pad][his]->SetLineWidth(2);
      Histo[pad][his]->SetLineColor(color[his]);
      Histo[pad][his]->Scale(nBinPDF[pad]/Histo[pad][his]->Integral()/(limX[pad][1]-limX[pad][0]));
      if(maxH[pad]<Histo[pad][his]->GetMaximum()) maxH[pad]=Histo[pad][his]->GetMaximum();
      if(pad==0) leg.AddEntry(Histo[pad][his], legName[his]);
    }
  }
    
  for(int pad=0; pad<2; pad++){
    can.cd(pad+1);
    for(int his=0; his<nHis; his++){
      if(his==0) {
	Histo[pad][his]->SetMaximum(maxH[pad]*1.05);
	Histo[pad][his]->Draw("c");
	style.setTitles(Histo[pad][his],xTitle[pad],yTitle[pad]);
      }else Histo[pad][his]->Draw("c same");
    }
    if(pad==1) leg.Draw();
  }
  TString epsName = "public_html/THDM_MC_"; epsName += PDFindex; epsName += ".eps";
  can.SaveAs(epsName);
  for(int pad=0; pad<2; pad++)
    for(int his=0; his<nHis; his++)
      Histo[pad][his]->Delete();
}



