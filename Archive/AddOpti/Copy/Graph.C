void Graph(){
  fstream optimi;
  optimi.open("modes/txt/TableRun24_optimization.txt",fstream::in);

  double values[9][1000], cut1, cut2, puri, signi;
  int S, Norm, Dss, Cross, Comb, row=0;
  TString rub;
  while(optimi && row<1000){
    optimi>>cut1>>cut2>>rub>>puri>>signi>>S>>Norm>>Dss>>Cross>>Comb;
    values[0][row] = cut1*100; values[1][row] = cut2*100;
    values[2][row] = puri; values[3][row] = signi; 
    values[4][row] = (double)S; values[5][row] = (double)Norm; 
    values[6][row] = (double)Dss; values[7][row] = (double)Cross; 
    values[8][row] = (double)Comb; 
    row++;
  }

  int npoints = 100, ngraphs = 3;
  TLine line;
  line.SetLineStyle(2);
  TLatex tag;
  tag.SetTextSize(.04);
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
    g[i]->SetLineWidth(2.9);
    if(i==0){
      g[i]->GetXaxis()->SetLimits(0,100.);
      g[i]->GetXaxis()->SetTitle("Cut on % of events with all tracks truthmatched");
      g[i]->SetTitle("Purity for different cuts");
      g[i]->Draw("ALP");
    }else{
      g[i]->Draw("LP");   
    }
    leg.AddEntry(g[i],names[i]);
  }
  line.SetLineColor(4); tag.SetTextColor(4);
  line.DrawLine(0,.427,100,.427);
  tag.DrawText(10,.41,"BSemiExcl, old PID");
  leg.Draw();
  c.SaveAs("modes/eps/purity.eps");
}


