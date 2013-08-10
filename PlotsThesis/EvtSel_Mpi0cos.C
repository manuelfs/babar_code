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

void EvtSel_Mpi0cos(){

  Styles style2; style2.setPadsStyle(2); style2.applyStyle();

  TString NameTrees[3] = {"AWG82/ntuples/small/RAll_RunAll.root", 
			  "AWG82/ntuples/small/uds_RunAll.root", "AWG82/ntuples/small/ccbar_RunAll.root"};
  TChain gen("ntp1"), cont("ntp1");
  gen.Add(NameTrees[0]);
  for(int t=1; t<3; t++) cont.Add(NameTrees[t]);
  double totMCB = 0, totuds = 0, totccbar = 0, totdata = 0, totOffdata = 0;
  getNumberB(NameTrees[0], "All", totMCB, totdata, totuds, totccbar, totOffdata);
  double wuds = totMCB/totuds*2.09/1.05;     

  TH1F *hCount = new TH1F("hCount","",100,-4,12);
  gen.Draw("mm2pi0>>hCount","weight");
  double nTotal = hCount->Integral();
  cont.Draw("mm2pi0>>hCount","weight");
  nTotal += hCount->Integral()*wuds;

  TLine line; line.SetLineStyle(2); line.SetLineColor(28); line.SetLineWidth(2);
  TArrow arrow; arrow.SetLineColor(28); arrow.SetFillColor(28); arrow.SetLineWidth(2);
  TCanvas can("can","Mpi0 and costhetaT cuts");
  can.Divide(2,1); TPad *cPad = (TPad *)can.cd(1);
  int bins[] = {50,40}, colors[2][4] = {{2,4,1,3},{2,4,1,1}};
  double xrange[2][2] = {{0.09,0.165},{-1,1}}, yield[2][4], maxi[] = {-99,-99};
  TString Variable[] = {"mpi0","candCosT"};
  TString labels[2][4] = {{"D** l #nu (", "B#bar{B} Bkg. (","Cont. Bkg. (",""},
 			  {"D** l #nu (", "BB bkg. (","Cont. bkg. (",""}};
  TString cuts[2][4] = {{"(mm2pi0>-4&&mm2pi0<12&&MCType>12)*weight", "(mm2pi0>-4&&mm2pi0<12&&MCType<13)*weight", 
			 "(mm2pi0>-4&&mm2pi0<12)*weight", ""},
			{"(mm2pi0>-4&&mm2pi0<12&&MCType>12)*weight", "(mm2pi0>-4&&mm2pi0<12&&MCType<13)*weight", 
			 "(mm2pi0>-4&&mm2pi0<12)*weight", ""}};
  double legW = 0.4, legH = 0.225;
  double legX = 0.65, legY = 1-style2.PadTopMargin-0.02;
  TLegend *leg[2];
  leg[0] = new TLegend(legX-legW, legY-legH, legX, legY);
  leg[1] = new TLegend(legX-legW, legY-legH, legX, legY);
  TH1F* h[2][4];
  for(int pad=0; pad<2; pad++){
    leg[pad]->SetTextSize(style2.LabelSize); leg[pad]->SetFillColor(0); 
    leg[pad]->SetTextFont(style2.nFont);  leg[pad]->SetBorderSize(0);
    for(int i=0; i<3; i++) {
      TString hname = "h"; hname += pad; hname += i;
      h[pad][i] = new TH1F(hname,"",bins[pad],xrange[pad][0],xrange[pad][1]);
      h[pad][i]->SetLineWidth(2);  h[pad][i]->SetLineColor(colors[pad][i]);
      TString vari = Variable[pad]; vari += ">>"; vari += hname;
      if(i<2){
	gen.Draw(vari,cuts[pad][i]);
	yield[pad][i] = h[pad][i]->Integral();
      } else {
	cont.Draw(vari,cuts[pad][i]);
	yield[pad][i] = h[pad][i]->Integral()*wuds;
      }
      h[pad][i]->Scale(1000/h[pad][i]->Integral());
      if(h[pad][i]->GetMaximum()>maxi[pad]) maxi[pad] = h[pad][i]->GetMaximum();
      labels[pad][i] += RoundNumber(yield[pad][i]*100,0,nTotal); labels[pad][i] += "%)";
      leg[pad]->AddEntry(h[pad][i],labels[pad][i]);
    }
    h[pad][0]->SetMaximum(maxi[pad]*1.13);
    h[pad][0]->SetMinimum(0);
  }
  cPad = (TPad *)can.cd(2);
  h[0][0]->Draw();
  style2.fixYAxis(h[0][0],cPad);
  style2.setTitles(h[0][0],"m(#gamma#gamma) (GeV)","Entries/(1.5 MeV)","b)");
  h[0][1]->Draw("same");h[0][2]->Draw("same");
  double Lh = 0.93, Ah = 0.1;
  line.DrawLine(0.12,h[0][0]->GetMinimum(), 0.12,maxi[0]*Lh);
  line.DrawLine(0.15,h[0][0]->GetMinimum(), 0.15,maxi[0]*Lh);
  arrow.DrawArrow(0.12,maxi[0]*(Lh-Ah),0.126,maxi[0]*(Lh-Ah),0.01,"|>");
  arrow.DrawArrow(0.15,maxi[0]*(Lh-Ah),0.144,maxi[0]*(Lh-Ah),0.01,"|>");


  cPad = (TPad *)can.cd(1);
  h[1][0]->Draw();
  style2.fixYAxis(h[1][0],cPad);
  style2.setTitles(h[1][0],"cos(#Delta#theta_{T})","Entries/(0.05)","a)");
  h[1][1]->Draw("same"); h[1][2]->Draw("same"); 
  leg[1]->Draw();
  Lh = 0.4;
  line.DrawLine(-0.8,h[1][0]->GetMinimum(), -0.8,maxi[1]*Lh);
  line.DrawLine(0.8,h[1][0]->GetMinimum(), 0.8,maxi[1]*Lh);
  arrow.DrawArrow(-0.8,maxi[1]*(Lh-Ah),-0.6,maxi[1]*(Lh-Ah),0.01,"|>");
  arrow.DrawArrow(0.8,maxi[1]*(Lh-Ah),0.6,maxi[1]*(Lh-Ah),0.01,"|>");

  TString pName = "public_html/EvtSel_Mpi0cos.eps"; 
  can.SaveAs(pName);
  for(int pad=0; pad<2; pad++){
    leg[pad]->Delete();
    for(int i=0; i<3; i++){
      h[pad][i]->Delete();
    }
  }
  hCount->Delete();
}

