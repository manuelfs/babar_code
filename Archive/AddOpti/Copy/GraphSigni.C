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

void GraphSigni(TString filename="optimB0truth12346_yieldRest6_cuts"){

  double values[2][9][1000], cut1, cut2, puri, signi;
  int S, Norm, Dss, Cross, Comb, row=0;
  TString rub;
  fstream optimiB0;
  optimiB0.open("babar_code/txt/"+filename+".txt",fstream::in);
  while(optimiB0 && row<1000){
    optimiB0>>cut1>>cut2>>rub>>puri>>signi>>S>>Norm>>Dss>>Cross>>Comb;
    values[0][0][row] = cut1*100; values[0][1][row] = cut2*100;
    values[0][2][row] = puri; values[0][3][row] = signi; 
    values[0][4][row] = (double)S; values[0][5][row] = (double)Norm; 
    values[0][6][row] = (double)Dss; values[0][7][row] = (double)Cross; 
    values[0][8][row] = (double)Comb; 
    row++;
  }
  row=0;
  filename.ReplaceAll("B0","Bp");
  fstream optimiBp;
  optimiBp.open("babar_code/txt/"+filename+".txt",fstream::in);
  while(optimiBp && row<1000){
    optimiBp>>cut1>>cut2>>rub>>puri>>signi>>S>>Norm>>Dss>>Cross>>Comb;
    values[1][0][row] = cut1*100; values[1][1][row] = cut2*100;
    values[1][2][row] = puri; values[1][3][row] = signi; 
    values[1][4][row] = (double)S; values[1][5][row] = (double)Norm; 
    values[1][6][row] = (double)Dss; values[1][7][row] = (double)Cross; 
    values[1][8][row] = (double)Comb; 
    row++;
  }

  TLine line;
  line.SetLineStyle(2);
  TLatex tag;
  tag.SetTextSize(.04);
  int npoints = 100, ngraphs = 2;
  int col[8] = {1,1,1,1,6,28,9,13};
  TString names[] = {"B^{0}","B^{+}"};
  TGraph *g[4];
  TCanvas c("Graph","Graph",1000,500);
  c.Divide(2,1);
  for(int i=0; i<ngraphs; i++){
    c.cd(i+1);
    g[i] = new TGraph(npoints,values[i][0],values[i][3]);
    g[i]->SetFillColor(0);
    g[i]->SetLineColor(col[i]);
    g[i]->SetLineWidth(2.9);
    g[i]->GetXaxis()->SetLimits(0,100.);
    g[i]->GetXaxis()->SetTitle("Cut on % of events with all tracks truthmatched");
    g[i]->SetTitle("Signal sample: Significance for "+names[i]);
    g[i]->Draw("ALP");
    line.SetLineColor(28); 
    line.DrawLine(30,0,30,1.07*values[i][3][0]);
  }
  c.SaveAs("babar_code/eps/Signi_"+filename+".eps");
}


