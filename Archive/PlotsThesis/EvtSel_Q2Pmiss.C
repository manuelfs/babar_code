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

void EvtSel_Q2Pmiss(){

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
  gen.Draw("candM2>>hCount","weight");
  double nTotal = hCount->Integral();
  cont.Draw("candM2>>hCount","weight");
  nTotal += hCount->Integral()*wuds;

  TLine line; line.SetLineStyle(2); line.SetLineColor(28); line.SetLineWidth(2);
  TArrow arrow; arrow.SetLineColor(28); arrow.SetFillColor(28); arrow.SetLineWidth(2);
  TCanvas can("can","Pmiss and q2 cuts");
  can.Divide(2,1); TPad *cPad = (TPad *)can.cd(1);
  int bins[] = {42,40}, colors[2][4] = {{8,4,1,3},{8,2,4,1}};
  double xrange[2][2] = {{0,4.2},{-3,13}}, yield[2][4], maxi[] = {-99,-99};
  TString Variable[] = {"candPMiss","candQ2"};
  TString labels[2][4] = {{"Signal (", "Normaliz. (","Had. bkg. (",""},
 			  {"Signal (", "D l #nu (", "D* l #nu (", "Bkg. ("}};
//   TString labels[2][4] = {{"Signal", "Normaliz.","Had. Bkg.",""},
// 			  {"Signal", "D l #nu", "D* l #nu", "Bkg."}};
  TString cuts[2][4] = {{"(candType<3&&MCType>4&&MCType<7||candType>2&&MCType>10&&MCType<13)*weight",
			 "(candType<3&&MCType>0&&MCType<5||candType>2&&MCType>6&&MCType<11)*weight",
			 "(MCType==0&&MCCombmode==12)*weight", ""},
			{"(candType<3&&MCType>4&&MCType<7||candType>2&&MCType>10&&MCType<13)*weight",
			 "(candType<3&&(MCType==1||MCType==3)||candType>2&&(MCType==7||MCType==9))*weight",
			 "(candType<3&&(MCType==2||MCType==4)||candType>2&&(MCType==8||MCType==10))*weight",
			 "(!(candType<3&&MCType>0&&MCType<7||candType>2&&MCType>6&&MCType<13))*weight"}};
  double legW = 0.4, legH = 0.225;
  double legX = 1-style2.PadRightMargin-0.02, legY = 1-style2.PadTopMargin-0.02;
  TLegend *leg[2];
  leg[0] = new TLegend(legX-legW, legY-legH, legX, legY);
  legW = 0.24; legH = 0.285; legX = 0.47;
  leg[1] = new TLegend(legX-legW, legY-legH, legX, legY);
  TH1F* h[2][4];
  for(int pad=0; pad<2; pad++){
    leg[pad]->SetTextSize(style2.LabelSize); leg[pad]->SetFillColor(0); 
    leg[pad]->SetTextFont(style2.nFont);  leg[pad]->SetBorderSize(0);
    for(int i=0; i<4; i++) {
      if(pad==0 && i==3) continue;
      TString hname = "h"; hname += pad; hname += i;
      h[pad][i] = new TH1F(hname,"",bins[pad],xrange[pad][0],xrange[pad][1]);
      h[pad][i]->SetLineWidth(2);  h[pad][i]->SetLineColor(colors[pad][i]);
      TString vari = Variable[pad]; vari += ">>"; vari += hname;
      gen.Draw(vari,cuts[pad][i]);
      if(i==3){
	hname = "hCont"; hname += pad; hname += i;
	TH1F *hCont = new TH1F(hname,"",bins[pad],xrange[pad][0],xrange[pad][1]);
	TString vari = Variable[pad]; vari += ">>"; vari += hname;
	cont.Draw(vari,cuts[pad][i]);
	hCont->Scale(wuds);
	h[pad][i]->Add(hCont);
	hCont->Delete();
      }
      yield[pad][i] = h[pad][i]->Integral();
      h[pad][i]->Scale(1000/h[pad][i]->Integral());
      if(h[pad][i]->GetMaximum()>maxi[pad]) maxi[pad] = h[pad][i]->GetMaximum();
      labels[pad][i] += RoundNumber(yield[pad][i]*100,0,nTotal); labels[pad][i] += "%)";
      leg[pad]->AddEntry(h[pad][i],labels[pad][i]);
    }
    h[pad][0]->SetMaximum(maxi[pad]*1.22);
  }
  h[0][0]->Draw();
  style2.fixYAxis(h[0][0],cPad);
  style2.setTitles(h[0][0],"|p_{miss}| (GeV)","Entries/(100 MeV)","a)");
  h[0][1]->Draw("same");h[0][2]->Draw("same");
  leg[0]->Draw();
  line.DrawLine(0.2,h[0][0]->GetMinimum(), 0.2,maxi[0]/1.45);
  arrow.DrawArrow(0.2,maxi[0]/1.65,0.5,maxi[0]/1.65,0.01,"|>");


  cPad = (TPad *)can.cd(2);
  h[1][0]->Draw();
  style2.fixYAxis(h[1][0],cPad);
  style2.setTitles(h[1][0],"q^{2} (GeV^{2})","Entries/(0.4 GeV^{2})","b)");
  h[1][1]->Draw("same"); h[1][2]->Draw("same"); h[1][3]->Draw("same");
  leg[1]->Draw();
  line.DrawLine(4,h[1][0]->GetMinimum(), 4,maxi[1]/1.45);
  arrow.DrawArrow(4,maxi[1]/1.65,5.4,maxi[1]/1.65,0.01,"|>");

  TString pName = "public_html/EvtSel_Q2Pmiss.eps"; 
  can.SaveAs(pName);
  for(int pad=0; pad<2; pad++){
    leg[pad]->Delete();
    for(int i=0; i<4; i++){
      if(pad==0 && i==3) continue;
      h[pad][i]->Delete();
    }
  }
  hCount->Delete();
}

