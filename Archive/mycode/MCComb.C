void MCComb(TString Dmode = "D0", TString when = "after", bool doScale = false){
  int n = 3;
  // Loading files
  TChain ch[4];
  int colors[]={1,2,4,8};
  for(int i=0; i<4; i++)   ch[i] = new TChain("ntp1");
  for(int i=0; i<4; i++)   ch[i].SetLineColor(colors[i]);

  ch[0].Add("AWG82/ntuples/NN/SP-1235-BSemiExcl-Run2-R22d-v08/*");
  ch[0].Add("AWG82/ntuples/NN/SP-1235-BSemiExcl-Run4-R22d-v08/*");
  ch[0].Add("AWG82/ntuples/NN/SP-1237-BSemiExcl-Run2-R22d-v08/*");
  ch[0].Add("AWG82/ntuples/NN/SP-1237-BSemiExcl-Run4-R22d-v08/*");

  ch[1].Add("AWG82/ntuples/BDT_deltaE/SP-1235-BSemiExcl-Run2-R22d-v08/*");
  ch[1].Add("AWG82/ntuples/BDT_deltaE/SP-1235-BSemiExcl-Run4-R22d-v08/*");
  ch[1].Add("AWG82/ntuples/BDT_deltaE/SP-1237-BSemiExcl-Run2-R22d-v08/*");
  ch[1].Add("AWG82/ntuples/BDT_deltaE/SP-1237-BSemiExcl-Run4-R22d-v08/*");

  ch[2].Add("AWG82/ntuples/NN_KM_deltaE/SP-1235-BSemiExcl-Run2-R22d-v08/*");
  ch[2].Add("AWG82/ntuples/NN_KM_deltaE/SP-1235-BSemiExcl-Run4-R22d-v08/*");
  ch[2].Add("AWG82/ntuples/NN_KM_deltaE/SP-1237-BSemiExcl-Run2-R22d-v08/*");
  ch[2].Add("AWG82/ntuples/NN_KM_deltaE/SP-1237-BSemiExcl-Run4-R22d-v08/*");


  TCut MMiss = "candPMiss>.2";
  TCut Q2 = "candQ2>4";
  TCut ee1 = "candType<3&&candEExtra<.2";
  TCut ee2 = "candType==3&&candEExtra<.15";
  TCut ee3 = "candType==4&&candEExtra<.3";
  TCut ee = (ee1||ee2||ee3);
  TCut basic = MMiss+Q2+ee;

  TH1F *histo[6];
  TString BinName[]={"D_{s} #rightarrow #tau","#tau no D_{s}","ee","#mu #mu","e #mu","D #rightarrow e",
  "D #rightarrow #mu","e no D","#mu no D",">3 l","K_{L}","Else"};
  TString names[]={"NN Base (","BDT (","eKMLoose ("};
  TCut cuts;
  TString Dname;
  if(Dmode == "D0"){
    cuts = "candType==1";
    Dname = "D^{0}";
  } else if(Dmode == "Dstar0"){
    cuts = "candType==2";
    Dname = "D*^{0}";
  } else if(Dmode == "Dplus"){
    cuts = "candType==3";
    Dname = "D^{+}";
  } else {
    cuts = "candType==4";
    Dname = "D*^{+}";
  } 
  cuts = "MCType==0&&candIsMu==0";
  if(when == "after")
    cuts += basic;
  TCanvas c("Yields","Yields",600,400);
  gStyle->SetOptStat(0);
  TLegend *leg;
  if(when == "after" && (Dmode == "D0" || Dmode == "Dstar0" || Dmode=="Dstar"))
    leg= new TLegend(0.63,1.-.087*n,0.9,0.9);
  else
    leg = new TLegend(0.1,1.-.087*n,0.4,0.9);
  leg->SetFillColor(0);
  float maxi = 0;
  int nJPsi=0;
  TCut cut3 = cuts;
  for(int i=0; i<n;i++){
    if(i==2){
      cuts = cuts+"MCCombmode<3"||cuts+"MCCombmode>=3&&candBMode>11000";
    }
    ch[i].Draw("MCCombmode>>h(12,1,13)",cuts);
    if(i==2){
      ch[i].GetEntries();
      nJPsi = ch[i].GetEntries(cut3+"MCCombmode>=3&&candBMode<11000");
}
    TH1F *h =  (TH1F*)gDirectory->Get("h");
    histo[i] = (TH1F*)h->Clone("histo"+i);
    histo[i]->Fill(12,nJPsi);
    cout<<"Entries "<<histo[i]->GetEntries()<<" and nJPsi "<<nJPsi<<endl;
    int entries = histo[i]->GetEntries()+nJPsi-1;
    if(doScale)
      histo[i]->Scale(100./entries);
    if(histo[i]->GetMaximum()>maxi) maxi = histo[i]->GetMaximum();
    names[i] += entries; names[i] += " entries)";
 }
  for(int i=0; i<n;i++){
    histo[i]->SetLineColor(colors[i]);
    histo[i]->SetLineWidth(2);
    if(i==0){
      histo[i]->SetTitle("Background types for the electron mode");
      //      histo[i]->SetTitle("Background types for the "+Dname+" "+when+" cuts");
      histo[i]->SetMaximum(1.05*maxi);
      histo[i]->Draw();
      for(int j=1;j<13;j++)
	histo[i]->GetXaxis()->SetBinLabel(j,BinName[j-1]);
    }else
      histo[i]->Draw("same");
    leg->AddEntry(histo[i],names[i]);
  }
  leg->Draw();
  if(doScale)
    c.SaveAs("mycode/eps/MCType/MCComb_Scaled_"+Dmode+"_"+when+"Cuts.eps");
  else
    c.SaveAs("mycode/eps/MCType/MCCombAlle_"+Dmode+"_"+when+"Cuts.eps");
}


