void GraphComb(TString Run="24"){
  fstream optimi;
  optimi.open("modes/txt/TableRun"+Run+"_opti_comb.txt",fstream::in);

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

  TString tagRun = "Runs 2 and 4";
  if(Run=="3") tagRun = "Run 3";
  int npoints = 100, ngraphs = 7;
  TLine line;
  line.SetLineStyle(2);
  TLatex tag;
  tag.SetTextSize(.04);
  int vari[] = {4, 3, 4, 5, 6, 7, 8, 9};
  int col[] = {4,1,3,2,4,6,28,9,13};
  TString names[] = {"Max plot", "D_{s} #rightarrow #tau","2 leptons","D #rightarrow e",
		     "D #rightarrow #mu","K_{L}","Else"};
  TGraph *g[8];
  TCanvas c("Graph","Graph",800,550);
  //TLegend leg(.71,.63,.86,.85);
  TLegend leg(.75,.64,.9,.9);
  leg.SetFillColor(0);
  leg.SetTextSize(0.035);
  //leg.SetBorderSize(0);
  for(int i=0; i<ngraphs; i++){
    g[i] = new TGraph(npoints,values[0],values[vari[i]]);
    g[i]->SetFillColor(0);
    g[i]->SetLineColor(col[i]);
    g[i]->SetLineWidth(2.9);
    if(i==0){
      g[i]->GetXaxis()->SetLimits(0,100.);
      g[i]->GetXaxis()->SetTitle("Cut on % of events with all tracks truthmatched");
      g[i]->SetTitle("Combinatoric bkg for different cuts: "+tagRun);
      g[i]->Draw("ALP");
    }else{
      g[i]->Draw("LP");   
      leg.AddEntry(g[i],names[i]);
    }
  }
  line.SetLineColor(4); tag.SetTextColor(4);
  //line.DrawLine(0,.427,100,.427);
  //tag.DrawText(10,.41,"BSemiExcl, old PID");
  leg.Draw();
  c.SaveAs("modes/eps/comb"+Run+".eps");
}


