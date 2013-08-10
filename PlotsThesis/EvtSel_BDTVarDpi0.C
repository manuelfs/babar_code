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
#include "babar_code/Styles/Styles.cc"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void EvtSel_BDTVarDpi0(){

  Styles style; style.setPadsStyle(8); style.applyStyle();

  TString NameTrees[3] = {"AWG82/ntuples/small/RAll_RunAll.root", 
			  "AWG82/ntuples/small/uds_RunAll.root", "AWG82/ntuples/small/ccbar_RunAll.root"};
  TChain gen("ntp1"), cont("ntp1");
  gen.Add(NameTrees[0]);
  for(int t=1; t<3; t++) cont.Add(NameTrees[t]);
  double totMCB = 0, totuds = 0, totccbar = 0, totdata = 0, totOffdata = 0;
  getNumberB(NameTrees[0], "All", totMCB, totdata, totuds, totccbar, totOffdata);
  double wuds = totMCB/totuds*2.09/1.05;     


  double nTotal, yield[4], maxi;
  TString cuts[] = {"(abs(candCosT)<.8&&mpi0>.12&&mpi0<.15&&pmisspi0>.2&&candQ2>4&&candType==2&&(MCType>12))*weight",
		    "(abs(candCosT)<.8&&mpi0>.12&&mpi0<.15&&pmisspi0>.2&&candQ2>4&&candType==2&&(MCType>0&&MCType<13))*weight",
		    "(abs(candCosT)<.8&&mpi0>.12&&mpi0<.15&&pmisspi0>.2&&candQ2>4&&candType==2&&(MCType==0))*weight",
		    "(abs(candCosT)<.8&&mpi0>.12&&mpi0<.15&&pmisspi0>.2&&candQ2>4&&candType==2)*weight"};
  TH1F *hCount = new TH1F("hCount","",100,-4,12);
  for(int i=0; i<3; i++) {
    gen.Draw("mm2pi0>>hCount",cuts[i]);
    yield[i] = hCount->Integral();
  }
  cont.Draw("mm2pi0>>hCount",cuts[3]);
  yield[3] = hCount->Integral()*wuds;
  gen.Draw("mm2pi0>>hCount",cuts[3]);
  nTotal = hCount->Integral();
  nTotal += yield[3];
  hCount->Delete();

  TCanvas can("can","BDT variables Dpi0");
  can.Divide(2,4); TPad *cPad;
  TString Variable[] = {"eextrapi0","candDeltaE","candDmass","candDeltam",
			"ppi0","dmpi0","e1pi0","candCosT"};
  int bins[] = {25,24,30,22,30,30,24,20}, colors[] = {4,2,1,8};
  double xrange[8][2] = {{-0.5,2}, {-0.072,0.072}, {1.94,2.06}, {0.13,0.152}, 
			 {0,1.5}, {0,1.5}, {0,1.2}, {-0.8,0.8}};
  TString yTitleEnd[] = {"100 MeV)", "6 MeV)", "4 MeV)", "1 MeV)", "50 MeV)", "50 MeV)", "50 MeV)",""};
  TString xTitle[] = {"E_{Extra} (GeV)", "#DeltaE (GeV)", "", "", 
		      "p(#pi^{0}) (GeV)", "", "max(E_{#gamma}) (GeV)", "cos(#Delta#theta_{T})"};
  TString labels[] = {"D** l #nu (", "Semilep. bkg. (", "Comb. bkg. (", "Cont. bkg. ("};
  TString PadLabel[] = {"a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)"};
  double legW = 0.5, legH = 0.4;
  double legX = 1-style.PadRightMargin-0.02, legY = 1-style.PadTopMargin-0.02;
  TLegend leg(legX-legW, legY-legH, legX, legY);
  leg.SetTextSize(style.LabelSize); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont);  leg.SetBorderSize(0);
  TH1F* h[8][4];
  for(int pad=0; pad<8; pad++){
    cPad = (TPad *)can.cd(pad+1);
    maxi = -99;
    for(int i=0; i<4; i++) {
      TString hname = "h"; hname += pad; hname += i;
      h[pad][i] = new TH1F(hname,"",bins[pad],xrange[pad][0],xrange[pad][1]);
      h[pad][i]->SetLineWidth(2);  h[pad][i]->SetLineColor(colors[i]);
      TString vari = Variable[pad]; vari += ">>"; vari += hname;
      if(i<3) gen.Draw(vari,cuts[i]);
      else cont.Draw(vari,cuts[i]);
      if(pad==0) {
	h[pad][i]->Integral();
	labels[i] += RoundNumber(yield[i]*100,0,nTotal); labels[i] += "%)";
	leg.AddEntry(h[pad][i],labels[i]);
      }
      h[pad][i]->Scale(1000/h[pad][i]->Integral());
      if(h[pad][i]->GetMaximum()>maxi) maxi = h[pad][i]->GetMaximum();
    }
    if(pad==4 || pad==6) maxi *= 1.1;
    h[pad][0]->SetMaximum(maxi*1.1);
    h[pad][0]->SetMinimum(0);
    h[pad][0]->Draw();
    if(pad==0) {
      style.fixYAxis(h[pad][0],cPad); style.applyStyle();
      leg.Draw();
    } else {
      h[pad][0]->SetTitleOffset(style.yTitleOffset,"y");
      cPad->SetLeftMargin(style.PadLeftMargin);
    }
    TString yTitle = "Entries/(";
    if(pad<7) yTitle += yTitleEnd[pad];
    else yTitle = "Entries/(0.05)";
    style.setTitles(h[pad][0],xTitle[pad],yTitle,PadLabel[pad]);
    h[pad][3]->Draw("same");h[pad][2]->Draw("same");h[pad][1]->Draw("same");h[pad][0]->Draw("same");
  }

  TString pName = "public_html/EvtSel_BDTVarDpi0.eps"; 
  can.SaveAs(pName);
  for(int pad=0; pad<8; pad++){
    for(int i=0; i<4; i++){
      h[pad][i]->Delete();
    }
  }
}

