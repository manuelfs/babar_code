#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TArrow.h"
#include "TChain.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCut.h"
#include "../DonutUtils/cuts.cc"
#include "babar_code/Styles/Styles.cc"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void Syst_MES_Eex(){

  Styles style; style.setPadsStyle(3); style.applyStyle();

  TChain gen("ntp1");
  gen.Add("AWG82/ntuples/small/FitRAllNewx100_RunAll.root");

  double nTotal, yield[4], maxi;
  TString cuts[] = {"(candType<5&&(MCType>=1&&MCType<=4||MCType>=7&&MCType<=10)&&candM2<2.5)*weight",
		    "(candType<5&&(MCType>=1&&MCType<=4||MCType>=7&&MCType<=10)&&candM2>=2.5)*weight",
		    "(candType<5&&MCType==0)*weight",
		    "(candType<5&&MCType<0)*weight"};


  TH1F *hCount = new TH1F("hCount","",100,-4,12);
  for(int i=0; i<4; i++) {
    gen.Draw("candM2>>hCount",cuts[i]);
    yield[i] = hCount->Integral();
  }
  gen.Draw("candM2>>hCount","(candType<5)*weight");
  nTotal = hCount->Integral();
  //nTotal += yield[3];
  hCount->Delete();

  TCanvas can("can","BDT variables");
  can.Divide(3,1); TPad *cPad;
//   TString Variable[] = {"candMES","candEExtra","candDmass","candDeltam",
// 			"candBTagDmass","candBTagDeltam","candTagChargedMult-0.5","candCosT"};
//   int bins[] = {15,10,40,22,24,22,7,40}, colors[] = {4,2,1,8};
//   double xrange[8][2] = {{5.27,5.291}, {0,0.5}, {1.94,2.06}, {0.13,0.152}, 
// 			 {1.841,1.889}, {0.13,0.152}, {1.5,8.5}, {-1,1}};
  TString Variable[] = {"candMES","candEExtra","candPstarLep","candDmass","candDeltam",
			"candBTagDmass","candBTagDeltam","candTagChargedMult-0.5","candCosT"};
  int bins[] = {15,10,15,22,24,22,7,40}, colors[] = {4,2,1,8};
  double xrange[8][2] = {{5.27,5.291}, {0,0.5}, {0,2.4}, {0.13,0.152}, 
			 {1.841,1.889}, {0.13,0.152}, {1.5,8.5}, {-1,1}};
  TString yTitleEnd[] = {" Mev)", " MeV)", " MeV)", "1 MeV)", "2 MeV)", "1 MeV)", "",""};
  TString xTitle[] = {"m_{ES} (Gev)", "E_{Extra} (Gev)", "p*_{l} (Gev)", "", 
		      "", "", "B_{tag} charged multiplicity", "cos(#Delta#theta_{T})"};
  TString labels[] = {"Norm. m_{miss}^{2}<2.5", "Norm. m_{miss}^{2}>2.5", "BB bkg.","Cont."};
  double legW = 0.55, legH = 0.21;
  double legX = 1-style.PadRightMargin-0.02, legY = 1-style.PadTopMargin-0.02;
  TLegend leg(legX-legW, legY-legH, legX, legY);
  leg.SetTextSize(style.LabelSize*0.95); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont);  leg.SetBorderSize(0);
  TH1F* h[8][4];
  for(int pad=0; pad<3; pad++){
    cPad = (TPad *)can.cd(pad+1);
    maxi = -99;
    for(int i=0; i<3; i++) {
      TString hname = "h"; hname += pad; hname += i;
      h[pad][i] = new TH1F(hname,"",bins[pad],xrange[pad][0],xrange[pad][1]);
      h[pad][i]->SetLineWidth(2);  h[pad][i]->SetLineColor(colors[i]);
      TString vari = Variable[pad]; vari += ">>"; vari += hname;
      gen.Draw(vari,cuts[i]);
      if(pad==0) {
	h[pad][i]->Integral();
	//labels[i] += RoundNumber(yield[i]*100,0,nTotal); labels[i] += "%)";
	leg.AddEntry(h[pad][i],labels[i]);
      }
      h[pad][i]->Scale(1000/h[pad][i]->Integral());
      if(h[pad][i]->GetMaximum()>maxi) maxi = h[pad][i]->GetMaximum();
    }
    h[pad][0]->SetMaximum(maxi*1.1);
    h[pad][0]->SetMinimum(0);
    h[pad][0]->Draw("");
    if(pad==1) leg.Draw();
    double factor = 1000; int digits = 0;
    TString yTitle = "Entries/("; yTitle += RoundNumber((xrange[pad][1]-xrange[pad][0])*factor,digits,(double)bins[pad]);
    yTitle += yTitleEnd[pad];
    style.setTitles(h[pad][0],xTitle[pad],yTitle);
    //h[pad][3]->Draw("same");
    h[pad][2]->Draw("same");h[pad][1]->Draw("same");h[pad][0]->Draw("same");
  }

  TString pName = "public_html/Syst_MES_Eex.eps"; 
  can.SaveAs(pName);
  for(int pad=0; pad<3; pad++){
    for(int i=0; i<3; i++){
      h[pad][i]->Delete();
    }
  }
}

