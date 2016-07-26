//---------------------------------------------------------------------------------
// Description:
//      Plots the significance for different cuts in the % of tracks truthmatched
//      for each Breco mode
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      09/05/23 manuelf -- Clean up from old GraphSigni.C
//---------------------------------------------------------------------------------

void GraphSigni(){

  double values[2][9][1000], cut1, cut2, puri, signi;
  double S, Norm, Dss, Cross, Comb;
  int row=0;
  TString rub;
  fstream optimiB0;
  optimiB0.open("babar_code/txt/Newoptim.txt",fstream::in);
  while(optimiB0 && row<1000){
    optimiB0>>cut1>>cut2>>rub>>puri>>signi>>S>>Norm>>Dss>>Cross>>Comb;
    values[0][0][row] = cut1*100; values[0][1][row] = cut2*100;
    values[0][2][row] = puri; values[0][3][row] = (S+Norm)/(S+Norm+Dss+Cross+Comb); 
    values[0][4][row] = (double)S; values[0][5][row] = (double)Norm; 
    values[0][6][row] = (double)Dss; values[0][7][row] = (double)Cross; 
    values[0][8][row] = (double)Comb; 
    row++;
  }

  int npoints = 100, ngraphs = 1;
  int col[8] = {1,1,1,1,6,28,9,13};
  TString names[] = {"B^{0}","B^{+}"};
  TGraph *g[4];
  TCanvas c("Graph","Graph",1000,500);
  for(int i=0; i<ngraphs; i++){
    g[i] = new TGraph(npoints,values[i][0],values[i][3]);
    g[i]->SetFillColor(0);
    g[i]->SetLineColor(col[i]);
    g[i]->SetLineWidth(2.9);
    g[i]->GetXaxis()->SetLimits(0,100.);
    g[i]->GetXaxis()->SetTitle("Cut on % of events with all tracks truthmatched");
    g[i]->SetTitle("Signal sample: Significance for "+names[i]);
    g[i]->Draw("ALP");
  }
  c.SaveAs("babar_code/eps/NewSigni.eps");
}


