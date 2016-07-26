//---------------------------------------------------------------------------------
// Description:
//      Plots the yields for different cuts in the % of tracks truthmatched
//      for each Breco mode
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      09/05/23 manuelf -- Clean up from old GraphYields.C
//---------------------------------------------------------------------------------

void GraphYields(TString filename="optimB0truth12346_yieldRest6_cuts"){

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

  int npoints = 100, ngraphs = 5;
  TLine line;
  line.SetLineStyle(2);
  TLatex tag;
  tag.SetTextSize(.04);
  int vari[] = {8, 4, 6, 7, 8};
  int col[] = {4,1,3,2,4,6,28,9,13};
  TString names[] = {"rubbish","Signal","D**","Crossfeed","Combinatoric"};
  TString titles[] = {"B^{0}","B^{+}"};
  TGraph *g[5][2];
  TCanvas c("Graph","Graph",1000,500);
  //TLegend leg(.69,.67,.85,.85);
  TLegend leg(.625,.73,.9,.9);
  leg.SetFillColor(0);
  leg.SetTextSize(0.035);
  //leg.SetBorderSize(0);
  c.Divide(2,1);
  for(int j=0; j<2; j++) {
    c.cd(j+1);
    if(j==1) vari[0] = 4;
    for(int i=0; i<ngraphs; i++){
      g[i][j] = new TGraph(npoints,values[j][0],values[j][vari[i]]);
      g[i][j]->SetFillColor(0);
      g[i][j]->SetLineColor(col[i]);
      g[i][j]->SetLineWidth(2.9);
      if(i==0){
	g[i][j]->GetXaxis()->SetLimits(0,100.);
	g[i][j]->GetXaxis()->SetTitle("Cut on % of events with all tracks truthmatched");
	g[i][j]->SetTitle("Yields for "+titles[j]);
	g[i][j]->Draw("ALP");
      }else{
	g[i][j]->Draw("LP");   
	if(j==0) leg.AddEntry(g[i][j],names[i]);
      }
    }
    leg.Draw();
    line.SetLineColor(28); 
    line.DrawLine(30,0,30,1.07*values[j][vari[0]][0]);
  }
  c.SaveAs("babar_code/eps/Yields_"+filename+".eps");
}


