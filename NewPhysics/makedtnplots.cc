void readGsfile(int gs, double &rd, double &rderr, double &rds, double &rdserr){
  TString Gsname = "babar_code/NewPhysics/CLN_gs/r_"; Gsname += gs;
  fstream Gsfile; Gsfile.open(Gsname,fstream::in);
  Gsfile>>rd>>rderr>>rds>>rdserr;
}

void readRD(double gs){

  double RD[2][4];
  int gs_main = (int)(gs*1000);
  double frac = gs*1000-gs_main;
  if(gs_main>=0 && gs_main<2000) {
    readGsfile(gs_main,RD[0][0],RD[0][1],RD[0][2],RD[0][3]);
    if(frac>0.001 && gs_main+1<2000){
      readGsfile(gs_main+1,RD[1][0],RD[1][1],RD[1][2],RD[1][3]);
      for(int i=0; i<4; i++)RD[0][i] = (1-frac)*RD[0][i] + frac*RD[1][i];
    }
  }
      

  cout<<RD[0][0]<<", "<<RD[0][1]<<", "<<RD[0][2]<<", "<<RD[0][3]<<", "<<frac<<endl;

}


void makedtnplots() {
  TFile f("gs_cln.root");
  TH1F *rd = (TH1F*) f.Get("rd");
  TH1F *rds = (TH1F*) f.Get("rds");
  TH1F *errd = (TH1F*) f.Get("errd");
  TH1F *errds = (TH1F*) f.Get("errds");
  TH1F *rdnoerr = (TH1F*) f.Get("rdnoerr");
  TH1F *rdsnoerr = (TH1F*) f.Get("rdsnoerr");
  TH1F *reld = (TH1F*) f.Get("reld");
  TH1F *relds = (TH1F*) f.Get("relds");
  TH1F *pol = (TH1F*) f.Get("pol");
  TH1F *fbD = (TH1F*) f.Get("fbD");
  TH1F *fbDs = (TH1F*) f.Get("fbDs");
  TH1F *correl = (TH1F*) f.Get("correl");

  TCanvas *c1 = new TCanvas;
//   rd->Draw("e3");
//   c1->SaveAs("br_dtn.eps");
//   rds->Draw("e3");
//   c1->SaveAs("br_dstn.eps");
//   pol->Draw("e3");
//   c1->SaveAs("longpol.eps");
//   fbD->Draw("e3");
//   c1->SaveAs("afb_d.eps");
//   fbDs->Draw("e3");
//   c1->SaveAs("afb_ds.eps");
//   correl->Draw("hist");
//   c1->SaveAs("correl.eps");

  TF1 *fun = new TF1("fun","pol2");
  FILE *tf = fopen("babar_code/NewPhysics/dtn_gs","wt");
  fun->SetLineColor(kRed);
  rdnoerr->SetLineWidth(8);
  rdnoerr->Fit("fun", "Q B");
  fprintf(tf,"RD %f %f %f\n",fun->GetParameter(0),
	  fun->GetParameter(1),fun->GetParameter(2));
  rdnoerr->Draw();
  c1->SaveAs("br_dtn_param.eps");

  errd->SetLineWidth(8);
  errd->Fit("fun", "Q B");
  fprintf(tf,"errD %f %f %f\n",fun->GetParameter(0),
	  fun->GetParameter(1),fun->GetParameter(2));
  errd->Draw();
  c1->SaveAs("err_dtn_param.eps");

  rdsnoerr->SetLineWidth(8);
  rdsnoerr->Fit("fun", "Q B");
  fprintf(tf,"RDs %f %f %f\n",fun->GetParameter(0),
	  fun->GetParameter(1),fun->GetParameter(2));
  rdsnoerr->Draw();
  c1->SaveAs("br_dstn_param.eps");

  errds->SetLineWidth(8);
  errds->Fit("fun", "Q B");
  fprintf(tf,"errDs %f %f %f\n",fun->GetParameter(0),
	  fun->GetParameter(1),fun->GetParameter(2));
  errds->Draw();
  c1->SaveAs("err_dstn_param.eps");

  fclose(tf);
  gApplication->Terminate();
}
