void GraphSigni(TString Run="24"){

  double values[4][9][1000], cut1, cut2, puri, signi;
  int S, Norm, Dss, Cross, Comb, row=0;
  TString rub;
  fstream optimi;
  optimi.open("modes/txt/TableRun"+Run+"_optimization.txt",fstream::in);
  while(optimi && row<1000){
    optimi>>cut1>>cut2>>rub>>puri>>signi>>S>>Norm>>Dss>>Cross>>Comb;
    values[0][0][row] = cut1*100; values[0][1][row] = cut2*100;
    values[0][2][row] = puri; values[0][3][row] = signi; 
    values[0][4][row] = (double)S; values[0][5][row] = (double)Norm; 
    values[0][6][row] = (double)Dss; values[0][7][row] = (double)Cross; 
    values[0][8][row] = (double)Comb; 
    row++;
  }
  row=0;
  fstream fDss;
  fDss.open("modes/txt/Dss"+Run+".txt",fstream::in);
  while(fDss && row<1000){
    fDss>>cut1>>cut2>>rub>>puri>>signi>>S>>Norm>>Dss>>Cross>>Comb;
    values[1][0][row] = cut1*100; values[1][1][row] = cut2*100;
    values[1][2][row] = puri; values[1][3][row] = signi; 
    values[1][4][row] = (double)S; values[1][5][row] = (double)Norm; 
    values[1][6][row] = (double)Dss; values[1][7][row] = (double)Cross; 
    values[1][8][row] = (double)Comb; 
    row++;
  }
  row=0;
  fstream Cr;
  Cr.open("modes/txt/Cr"+Run+".txt",fstream::in);
  while(Cr && row<1000){
    Cr>>cut1>>cut2>>rub>>puri>>signi>>S>>Norm>>Dss>>Cross>>Comb;
    values[2][0][row] = cut1*100; values[2][1][row] = cut2*100;
    values[2][2][row] = puri; values[2][3][row] = signi; 
    values[2][4][row] = (double)S; values[2][5][row] = (double)Norm; 
    values[2][6][row] = (double)Dss; values[2][7][row] = (double)Cross; 
    values[2][8][row] = (double)Comb; 
    row++;
  }
  row=0;
  fstream Co;
  Co.open("modes/txt/Co"+Run+".txt",fstream::in);
  while(Co && row<1000){
    Co>>cut1>>cut2>>rub>>puri>>signi>>S>>Norm>>Dss>>Cross>>Comb;
    values[3][0][row] = cut1*100; values[3][1][row] = cut2*100;
    values[3][2][row] = puri; values[3][3][row] = signi; 
    values[3][4][row] = (double)S; values[3][5][row] = (double)Norm; 
    values[3][6][row] = (double)Dss; values[3][7][row] = (double)Cross; 
    values[3][8][row] = (double)Comb; 
    row++;
  }

  TLine line;
  line.SetLineStyle(2);
  TLatex tag;
  tag.SetTextSize(.04);
  int npoints = 100, ngraphs = 4;
  //double AddSigni[] = {297.66, 640.34, 732.3, 538.3};
  double AddSigni[4];
  for(int i=0;i<4;i++)AddSigni[i] = values[i][3][0];
  double AddSigni3[] = {53.49, 113.88,132., 100.86};  
  double noAddSigni[] = {296.2, 413.4, 567.3, 481.6};
  double noAddSigni3[] = {46.08, 66.15, 94.80, 80.45};
  int col[8] = {1,1,1,1,6,28,9,13};
  TString tagRun = "runs 2 and 4";
  if(Run=="3") {
    tagRun = "run 3";
    for(int i=0; i<4; i++){
      AddSigni[i] = AddSigni3[i];
      noAddSigni[i] = noAddSigni3[i];
    }
  }
  TString names[] = {"D**+Cross+Comb","D**","Cross","Comb"};
  TGraph *g[4];
  TCanvas c("Graph","Graph",1000,550);
  c.Divide(2,2);
  for(int i=0; i<ngraphs; i++){
    c.cd(i+1);
    g[i] = new TGraph(npoints,values[i][0],values[i][3]);
    g[i]->SetFillColor(0);
    g[i]->SetLineColor(col[i]);
    g[i]->SetLineWidth(2.9);
    g[i]->GetXaxis()->SetLimits(0,100.);
    g[i]->GetXaxis()->SetTitle("Cut on % of events with all tracks truthmatched");
    g[i]->SetTitle("Significance in "+tagRun+" for different cuts: B="+names[i]);
    g[i]->Draw("ALP");
    line.SetLineColor(2); tag.SetTextColor(2);
    line.DrawLine(0,AddSigni[i],100,AddSigni[i]);
    tag.DrawText(34,AddSigni[i]*1.03,"BSemiExclAdd, new PID");
    line.SetLineColor(4); tag.SetTextColor(4);
    tag.DrawText(38,noAddSigni[i]*.925,"BSemiExcl, old PID");
    line.DrawLine(0,noAddSigni[i],100,noAddSigni[i]);
  }
  c.SaveAs("modes/eps/significance"+Run+".eps");
}


