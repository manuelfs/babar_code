#include "TString.h"
#include "TFile.h"
#include "TLine.h"
#include "TCanvas.h"
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

void EvtSel_BtagYield(){

  Styles style; style.setPadsStyle(1); style.applyStyle();
  fstream optimi;
  optimi.open("babar_code/txt/Newoptim.txt",fstream::in);

  double values[9][1000], cut1, cut2, puri, signi;
  double S, Norm, Dss, Cross, Comb, f=3;
  int row=0;
  TString rub;
  while(optimi && row<1000){
    optimi>>cut1>>cut2>>rub>>puri>>signi>>S>>Norm>>Dss>>Cross>>Comb;
    values[0][row] = cut1*100; values[1][row] = cut2;
    values[2][row] = (S+Norm)/(S+Norm+Dss+Cross+Comb); 
    values[3][row] = (S+Norm)/sqrt(S+Norm+Dss+Cross+Comb); 
    values[4][row] = (S+Norm)/f; values[5][row] = Norm; 
    values[6][row] = Dss/f; values[7][row] = Cross/f; 
    values[8][row] = Comb/f; 
    row++;
    //cout<<row<<": "<<cut1<<"\t"<<cut2<<"\t"<<puri<<"\t"<<signi<<endl;
  }

  int npoints = 100, ngraphs = 5;
  int vari[] = {8, 4, 8, 6, 7};
  int col[8] = {4,8,2,1,4,28,9,13};
  TString names[] = {"rubbish","Signal + norm.","Combinatoric","D** l #nu","Crossfeed"};

  TGraph *g[5];
  TCanvas c("Graph","Purity and significance for Btag truth-mathched cuts");
  
  double legW = 0.25, legH = 0.3;
  double legX = 1-style.PadRightMargin-0.02, legY = 1-style.PadTopMargin-0.04;
  TLegend leg(legX-legW, legY-legH, legX, legY);
  leg.SetTextSize(style.TitleSize); leg.SetFillColor(0); leg.SetTextFont(style.nFont);
  leg.SetBorderSize(0);

  for(int i=0; i<ngraphs; i++){
    g[i] = new TGraph(npoints,values[0],values[vari[i]]);
    g[i]->SetFillColor(0);
    g[i]->SetLineColor(col[i]);
    g[i]->SetLineWidth(2);
    if(i==0){
      g[i]->GetXaxis()->SetLimits(0,100.);
      g[i]->GetXaxis()->SetTitle("Cut on R_{tc} (%)");
      g[i]->GetYaxis()->SetTitle("Yields");
      g[i]->SetTitle("");
      g[i]->Draw("ALP");
    }else{
      g[i]->Draw("LP");   
      leg.AddEntry(g[i],names[i]);
    }
  }
  leg.Draw();
  TLine line;  line.SetLineStyle(2);  line.SetLineColor(28); line.SetLineWidth(2);
  line.DrawLine(30,0,30,values[vari[0]][0]*1.1);
  c.SaveAs("public_html/EvtSel_BtagYield.eps");
}


