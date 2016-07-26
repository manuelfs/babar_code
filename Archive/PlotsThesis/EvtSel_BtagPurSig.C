#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TPad.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "babar_code/Styles/Styles.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void EvtSel_BtagPurSig(){

  Styles style2; style2.setPadsStyle(2); style2.applyStyle();
  fstream optimi;
  optimi.open("babar_code/txt/Newoptim.txt",fstream::in);

  double values[9][1000], cut1, cut2, puri, signi;
  double S, Norm, Dss, Cross, Comb;
  int row=0;
  TString rub;
  while(optimi && row<1000){
    optimi>>cut1>>cut2>>rub>>puri>>signi>>S>>Norm>>Dss>>Cross>>Comb;
    values[0][row] = cut1*100; values[1][row] = cut2;
    values[2][row] = (S+Norm)/(S+Norm+Dss+Cross+Comb); 
    values[3][row] = (S+Norm)/sqrt(S+Norm+Dss+Cross+Comb); 
    values[4][row] = S; values[5][row] = Norm; 
    values[6][row] = Dss; values[7][row] = Cross; 
    values[8][row] = Comb; 
    row++;
    //cout<<row<<": "<<cut1<<"\t"<<cut2<<"\t"<<puri<<"\t"<<signi<<endl;
  }

  int npoints = 100, ngraphs = 3;
  int col[8] = {1,4,2,4,6,28,9,13};
  TString names[] = {"No neutrals cut","< 1 non true neutrals",
	       "< 0.5 non true neutrals"};
  TGraph *g[2][3];
  TCanvas c("Graph","Purity and significance for Btag truth-mathched cuts");
  c.Divide(2,1);
  
  TLine line[2];  
  for(int i=0; i<2; i++){line[i].SetLineStyle(2);  line[i].SetLineColor(28);} 
  TLatex label; label.SetNDC(kTRUE); label.SetTextAlign(13);
  TString Ytitle[] = {"S/(S+B)", "S/#surd(S+B)"};
  TString Tag[] = {"a)", "b)"};
  double legW = 0.55, legH = 0.22;
  double legX = 1-style2.PadRightMargin-0.01, legY = style2.PadBottomMargin+legH+0.05;
  TLegend leg(legX-legW, legY-legH, legX, legY);
  leg.SetTextSize(style2.LabelSize); leg.SetFillColor(0); leg.SetTextFont(style2.nFont);
  leg.SetBorderSize(0);
  for(int pad=0; pad<2; pad++){
    c.cd(pad+1);
    for(int i=0; i<ngraphs; i++){
      g[pad][i] = new TGraph(npoints,values[0]+npoints*i,values[2+pad]+npoints*i);
      g[pad][i]->SetFillColor(0);
      g[pad][i]->SetLineColor(col[i]);
      g[pad][i]->SetLineWidth(2);
      if(i==0){
	g[pad][i]->GetXaxis()->SetLimits(0,100.);
	if(pad==0) g[pad][i]->SetMinimum(0);
	g[pad][i]->SetTitle("");
	g[pad][i]->GetXaxis()->SetTitle("Cut on R_{tc} (%)");
	g[pad][i]->GetYaxis()->SetTitle(Ytitle[pad]);
	g[pad][i]->Draw("ACL");
      }else{
	g[pad][i]->Draw("CL");   
	if(i==2) label.DrawLatex(style2.PadLeftMargin+0.04,1-style2.PadTopMargin-0.03,Tag[pad]);  
      }
      if(pad==0) leg.AddEntry(g[pad][i],names[i]);
    }
    if(pad==0) leg.Draw();
    if(pad==1) line[pad].DrawLine(30,g[pad][0]->GetYaxis()->GetXmin(),
		       30,g[pad][0]->GetYaxis()->GetXmax());
    else line[pad].DrawLine(30,0,30,g[pad][0]->GetYaxis()->GetXmax());
  }
  c.SaveAs("public_html/EvtSel_BtagPurSig.eps");
  for(int pad=0; pad<2; pad++)
    for(int i=0; i<ngraphs; i++)g[pad][i]->Delete();
 
}


