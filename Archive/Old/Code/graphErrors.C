{

    TString names[] = {"Results_New_skim_cos","Results_Old_skim_cos"};
//    TString names[] = {"Results_New_skim_cospl","Results_Old_skim_cospl"};
//    TString names[] = {"Results_New_skim_pl1","Results_Old_skim_pl1"};

TLegend leg(.5,.7,.88,.88);
leg.SetFillColor(0);
leg.SetTextSize(0.025);
leg.SetBorderSize(0);
int col[8] = {1,2,3,4,6,28,9,13};
TGraphErrors *g[6];
for(int j=0; j<2; j++){
  fstream textfile;
  textfile.open("Run1comp/"+names[j]+".txt",fstream::in);
  double purity[10], signal[10], bkg[10], signalerr[10], bkgerr[10], purerr[10];
  int number;
  TString rubbish;
  for(int i=0; i<7;i++)
    textfile>>rubbish;
  for(int i=0; i<10;i++){
    textfile>>number>>purity[i]>>signal[i]>>signalerr[i]>>bkg[i]>>bkgerr[i];
    purerr[i] = bkg[i]/pow(signal[i]+bkg[i],2)*signalerr[i]+signal[i]/pow(signal[i]+bkg[i],2)*bkgerr[i];
    //purity[i] = signal[i]/(sqrt(signal[i]+bkg[i]));
  }
  g[j] = new TGraphErrors(10,purity,signal,purerr,signalerr);
  g[j]->SetMarkerSize(.7);
  g[j]->SetMarkerStyle(20);
  g[j]->SetMarkerColor(col[j]);
  g[j]->SetFillColor(0);
  g[j]->SetLineColor(col[j]);
  if(j==0){
    g[j]->GetXaxis()->SetTitle("Purity");
     g[j]->GetXaxis()->SetLimits(0,1.);
    g[j]->GetYaxis()->SetLimits(0.,21000.);
    //    g[j]->GetYaxis()->SetTitle("Signal");
    g[j]->SetTitle("Signal Vs Purity");
    g[j]->Draw("ALP");
  }else
    g[j]->Draw("LP");   
  g[0]->GetYaxis()->SetLimits(0.,22000.);
  leg.AddEntry(g[j],names[j]);

}
 leg.Draw();
textfile.close();


}
