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

void Intro_M2Pl(){

  Styles style; style.setPadsStyle(2); style.applyStyle();

  TChain gen("ntp1");
  gen.Add("AWG82/ntuples/small/FitRAllNewx100_RunAll.root");

  double nTotal, yield[4], maxi;
  TString cuts[] = {//"(candType==1&&MCType==5||candType==2&&MCType==6||candType==3&&MCType==11||candType==4&&MCType==12)*weight",
		    "(candType<5&&(MCType>=5&&MCType<=6||MCType>=11&&MCType<=12))*weight",
		    "(candType<5&&(MCType>0&&MCType<5||MCType>6&&MCType<11))*weight",
		    //"(candType==1&&(MCType==1||MCType==3)||candType==2&&(MCType==2||MCType==4)||candType==3&&(MCType==7||MCType==9)||candType==4&&(MCType==8||MCType==10))*weight",
		    "(candType<5&&MCType==14)*weight",
		    "(candType<5&&MCType<=0)*weight"};
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
  can.Divide(2,1); TPad *cPad;
  TString Variable[] = {"candM2","candPstarLep","candDmass","candDeltam",
			"candBTagDmass","candBTagDeltam","candTagChargedMult-0.5","candCosT"};
  int bins[] = {15,15,40,22,24,22,7,40}, colors[] = {4,2,1,8};
  double xrange[8][2] = {{-2,10}, {0,2.4}, {1.94,2.06}, {0.13,0.152}, 
			 {1.841,1.889}, {0.13,0.152}, {1.5,8.5}, {-1,1}};
  TString yTitleEnd[] = {" Gev^{2})", " MeV)", "3 MeV)", "1 MeV)", "2 MeV)", "1 MeV)", "",""};
  TString xTitle[] = {"m^{2}_{miss} (Gev^{2})", "p*_{l} (Gev)", "", "", 
		      "", "", "B_{tag} charged multiplicity", "cos(#Delta#theta_{T})"};
  TString labels[] = {"Signal (", "Normalization (", "D**l#nu (", "Comb. + Cont. ("};
  TString PadLabel[] = {"a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)"};
  double legW = 0.52, legH = 0.36;
  double legX = 1-style.PadRightMargin-0.02, legY = 1-style.PadTopMargin-0.02;
  TLegend leg(legX-legW, legY-legH, legX, legY);
  leg.SetTextSize(style.LabelSize); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont);  leg.SetBorderSize(0);
  TH1F* h[8][4];
  for(int pad=0; pad<2; pad++){
    cPad = (TPad *)can.cd(pad+1);
    maxi = -99;
    for(int i=0; i<4; i++) {
      TString hname = "h"; hname += pad; hname += i;
      h[pad][i] = new TH1F(hname,"",bins[pad],xrange[pad][0],xrange[pad][1]);
      h[pad][i]->SetLineWidth(2);  h[pad][i]->SetLineColor(colors[i]);
      TString vari = Variable[pad]; vari += ">>"; vari += hname;
      gen.Draw(vari,cuts[i]);
      if(pad==0) {
	h[pad][i]->Integral();
	labels[i] += RoundNumber(yield[i]*100,0,nTotal); labels[i] += "%)";
	leg.AddEntry(h[pad][i],labels[i]);
      }
      h[pad][i]->Scale(1000/h[pad][i]->Integral());
      if(h[pad][i]->GetMaximum()>maxi) maxi = h[pad][i]->GetMaximum();
    }
    if(pad==0) maxi *= 0.5;
    h[pad][0]->SetMaximum(maxi*1.1);
    h[pad][0]->SetMinimum(0);
    h[pad][0]->Draw("c");
    if(pad==0) {
      style.fixYAxis(h[pad][0],cPad); style.applyStyle();
      leg.Draw();
    } else {
      h[pad][0]->SetTitleOffset(style.yTitleOffset,"y");
      cPad->SetLeftMargin(style.PadLeftMargin);
    }
    double factor = 1; int digits = 2;
    if(pad==1){factor = 1000; digits = 0;}
    TString yTitle = "Entries/("; yTitle += RoundNumber((xrange[pad][1]-xrange[pad][0])*factor,digits,(double)bins[pad]);
    yTitle += yTitleEnd[pad];
    style.setTitles(h[pad][0],xTitle[pad],yTitle,PadLabel[pad]);
    h[pad][3]->Draw("c same");h[pad][2]->Draw("c same");h[pad][1]->Draw("c same");h[pad][0]->Draw("c same");
    //h[pad][3]->Draw("same");h[pad][2]->Draw("same");h[pad][1]->Draw("same");h[pad][0]->Draw("same");
  }

  TString pName = "public_html/Intro_M2Pl.eps"; 
  can.SaveAs(pName);
  for(int pad=0; pad<2; pad++){
    for(int i=0; i<4; i++){
      h[pad][i]->Delete();
    }
  }
}

