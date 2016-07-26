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

void EvtSel_BDTVar(){

  Styles style; style.setPadsStyle(-8); style.applyStyle();

  TString NameTrees[3] = {"AWG82/ntuples/small/RAll_RunAll.root", 
  //TString NameTrees[3] = {"AWG82/ntuples/small/RAll_Run6.root", 
			  "AWG82/ntuples/small/uds_RunAll.root", "AWG82/ntuples/small/ccbar_RunAll.root"};
  TChain gen("ntp1"), cont("ntp1");
  gen.Add(NameTrees[0]);
  for(int t=1; t<3; t++) cont.Add(NameTrees[t]);
  double totMCB = 0, totuds = 0, totccbar = 0, totdata = 0, totOffdata = 0;
  getNumberB(NameTrees[0], "All", totMCB, totdata, totuds, totccbar, totOffdata);
  double wuds = totMCB/totuds*2.09/1.05;     

  double nTotal, yield[4], maxi;
  TString cuts[] = {"(candPMiss>0.2&&candQ2>4&&candType==2&&(MCType==2||MCType==4||MCType==6))*weight",
		    "(candPMiss>0.2&&candQ2>4&&candType==2&&(MCType==1||MCType==3||MCType==5||MCType>6))*weight",
		    "(candPMiss>0.2&&candQ2>4&&candType==2&&(MCType==0))*weight",
		    "(candPMiss>0.2&&candQ2>4&&candType==2)*weight"};
//   TString cuts[] = {"(candMvaDl>0.48&&candPMiss>0.2&&candQ2>4&&candType==2&&(MCType==2||MCType==4||MCType==6))*weight",
// 		    //"(candPMiss>0.2&&candQ2>4&&candType==2&&(MCType==1||MCType==3||MCType==5||MCType>6))*weight",
// 		    "(candMvaDl>0.48&&candPMiss>0.2&&candQ2>4&&candType==2&&(MCType>12))*weight",
// 		    "(candMvaDl>0.48&&candPMiss>0.2&&candQ2>4&&candType==2&&(MCType==0))*weight",
// 		    "(candMvaDl>0.48&&candPMiss>0.2&&candQ2>4&&candType==2)*weight"};
  TH1F *hCount = new TH1F("hCount","",100,-4,12);
  for(int i=0; i<3; i++) {
    gen.Draw("candM2>>hCount",cuts[i]);
    yield[i] = hCount->Integral();
  }
  cont.Draw("candM2>>hCount",cuts[3]);
  yield[3] = hCount->Integral()*wuds;
  gen.Draw("candM2>>hCount",cuts[3]);
  nTotal = hCount->Integral();
  nTotal += yield[3];
  hCount->Delete();

  TCanvas can("can","BDT variables");
  can.Divide(2,4); TPad *cPad;
  TString Variable[] = {"candEExtra","candDeltaE","candDmass","candDeltam",
			"candBTagDmass","candBTagDeltam","candTagChargedMult-0.5","candCosT"};
  int bins[] = {32,24,40,22,24,22,7,40}, colors[] = {4,2,8,1}, lineStyle[] = {1,2,4,3};
  double xrange[8][2] = {{0,1.6}, {-0.072,0.072}, {1.94,2.06}, {0.13,0.152}, 
			 {1.841,1.889}, {0.13,0.152}, {1.5,8.5}, {-1,1}};
//   int bins[] = {10,24,40,22,20,15,7,40}, colors[] = {4,2,1,8};
//   double xrange[8][2] = {{0,0.5}, {-0.072,0.072}, {1.94,2.06}, {0.13,0.152}, 
// 			 {1.841,1.889}, {0.13,0.152}, {1.5,8.5}, {-1,1}};
  TString yTitleEnd[] = {"50 MeV)", "6 MeV)", "3 MeV)", "1 MeV)", "2 MeV)", "1 MeV)", "",""};
  TString xTitle[] = {"E_{Extra} (GeV)", "#DeltaE (GeV)", "", "", 
		      "", "", "B_{tag} charged multiplicity", "cos(#Delta#theta_{T})"};
  //TString labels[] = {"Signal + norm. (", "Semilep. bkg. (", "Comb. bkg. (", "Cont. bkg. ("};
  TString labels[] = {"a", "b", "c", "d"};
  TString PadLabel[] = {"(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"};
  double legW = 0.52, legH = 0.44;
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
      h[pad][i]->SetLineWidth(3);  h[pad][i]->SetLineColor(colors[i]);
      h[pad][i]->SetLineStyle(lineStyle[i]);
      TString vari = Variable[pad]; vari += ">>"; vari += hname;
      if(i<3) gen.Draw(vari,cuts[i]);
      else cont.Draw(vari,cuts[i]);
      if(pad==0) {
	h[pad][i]->Integral();
	//labels[i] += RoundNumber(yield[i]*100,0,nTotal); labels[i] += "%)";
	leg.AddEntry(h[pad][i],labels[i]);
      }
      h[pad][i]->Scale(1000/h[pad][i]->Integral());
      if(h[pad][i]->GetMaximum()>maxi) maxi = h[pad][i]->GetMaximum();
    }
    if(pad==0) maxi *= 1.1;
    if(pad==6) h[pad][0]->SetNdivisions(207,"x");
    h[pad][0]->SetLabelOffset(0.01,"Y");
    h[pad][0]->SetMaximum(maxi*1.1);
    h[pad][0]->SetMinimum(0);
    h[pad][0]->GetXaxis()->CenterTitle(true);
    h[pad][0]->GetYaxis()->CenterTitle(true);
    h[pad][0]->Draw();
//     if(pad==0) {
//       style.fixYAxis(h[pad][0],cPad); style.applyStyle();
//       leg.Draw();
//     } else {
//       h[pad][0]->SetTitleOffset(style.yTitleOffset,"y");
//       cPad->SetLeftMargin(style.PadLeftMargin);
//     }
    if(pad==0){
      leg.Draw();
    }
    TString yTitle = "Entries/(";
    if(pad<6) yTitle += yTitleEnd[pad];
    else if(pad==6) yTitle = "Entries";
    else yTitle = "Entries/(0.05)";
    xTitle[pad] = "t"; xTitle[pad] += pad;
    style.setTitles(h[pad][0],xTitle[pad],yTitle,PadLabel[pad]);
    h[pad][3]->Draw("same");h[pad][2]->Draw("same");h[pad][1]->Draw("same");h[pad][0]->Draw("same");
  }

  TString pName = "public_html/EvtSel_BDTVar.eps"; 
  can.SaveAs(pName);
  for(int pad=0; pad<8; pad++){
    for(int i=0; i<4; i++){
      h[pad][i]->Delete();
    }
  }
}

