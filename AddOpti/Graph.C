#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraph.h"
#include "TChain.h"
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

void Graph(){
  fstream optimi;
  optimi.open("babar_code/txt/Newoptim.txt",fstream::in);

  double values[9][1000], cut1, cut2, puri, signi;
  double S, Norm, Dss, Cross, Comb;
  int row=0;
  TString rub;
  while(optimi && row<1000){
    optimi>>cut1>>cut2>>rub>>puri>>signi>>S>>Norm>>Dss>>Cross>>Comb;
    values[0][row] = cut1*100; values[1][row] = cut2;
    values[2][row] = (S+Norm)/sqrt(S+Norm+Dss+Cross+Comb); 
    values[3][row] = signi; 
    values[4][row] = S; values[5][row] = Norm; 
    values[6][row] = Dss; values[7][row] = Cross; 
    values[8][row] = Comb; 
    row++;
    //cout<<row<<": "<<cut1<<"\t"<<cut2<<"\t"<<puri<<"\t"<<signi<<endl;
  }

  int npoints = 100, ngraphs = 3;
  int col[8] = {1,3,46,4,6,28,9,13};
  TString names[] = {"No neutrals cut","<2 non-truthmatched neutrals",
	       "<1 non-truthmatched neutrals"};
  TGraph *g[3];
  TCanvas c("Graph","Graph",800,550);
  TLegend leg(.12,.75,.4,.88);
  leg.SetFillColor(0);
  leg.SetTextSize(0.03);
  leg.SetBorderSize(0);
  for(int i=0; i<ngraphs; i++){
    g[i] = new TGraph(npoints,values[0]+npoints*i,values[2]+npoints*i);
//     g[i]->SetMarkerSize(.7);
//     g[i]->SetMarkerStyle(20);
//     g[i]->SetMarkerColor(col[i]);
    g[i]->SetFillColor(0);
    g[i]->SetLineColor(col[i]);
    g[i]->SetLineWidth(2);
    if(i==0){
      g[i]->GetXaxis()->SetLimits(0,100.);
      //g[i]->GetYaxis()->SetLimits(0,g[i]->GetMaximum()*1.2);
      g[i]->SetMinimum(0);
      g[i]->GetXaxis()->SetTitle("Cut on % of events with all tracks truthmatched");
      g[i]->SetTitle("Purity for different cuts");
      g[i]->Draw("CALP");
    }else{
      g[i]->Draw("CLP");   
    }
    leg.AddEntry(g[i],names[i]);
  }
//   TLine line;
//   line.SetLineStyle(2);
//   TLatex tag;
//   tag.SetTextSize(.04);
//   line.SetLineColor(4); tag.SetTextColor(4);
//   line.DrawLine(0,.427,100,.427);
//   tag.DrawText(10,.41,"BSemiExcl, old PID");
  //leg.Draw();
  c.SaveAs("babar_code/eps/Newpurity.eps");
}


