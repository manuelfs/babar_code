void SinglePlot(TString Dmode = "D0", TString when = "after", bool doScale = false){
  int n = 2;
  // Loading files
  TChain ch[3];
  for(int i=0; i<3; i++)   ch[i] = new TChain("ntp1");

  ch[0].Add("AWG82/ntuples/NN/SP-1235-BSemiExcl-Run2-R22d-v08/*");
  ch[0].Add("AWG82/ntuples/NN/SP-1235-BSemiExcl-Run4-R22d-v08/*");
  ch[0].Add("AWG82/ntuples/NN/SP-1237-BSemiExcl-Run2-R22d-v08/*");
  ch[0].Add("AWG82/ntuples/NN/SP-1237-BSemiExcl-Run4-R22d-v08/*");

  ch[1].Add("AWG82/ntuples/NN_Tweaked/SP-1235-BSemiExcl-Run2-R22d-v08/*");
  ch[1].Add("AWG82/ntuples/NN_Tweaked/SP-1235-BSemiExcl-Run4-R22d-v08/*");
  ch[1].Add("AWG82/ntuples/NN_Tweaked/SP-1237-BSemiExcl-Run2-R22d-v08/*");
  ch[1].Add("AWG82/ntuples/NN_Tweaked/SP-1237-BSemiExcl-Run4-R22d-v08/*");

//   ch[1].Add("AWG82/ntuples/BDT_deltaE/SP-1235-BSemiExcl-Run2-R22d-v08/*");
//   ch[1].Add("AWG82/ntuples/BDT_deltaE/SP-1235-BSemiExcl-Run4-R22d-v08/*");
//   ch[1].Add("AWG82/ntuples/BDT_deltaE/SP-1237-BSemiExcl-Run2-R22d-v08/*");
//   ch[1].Add("AWG82/ntuples/BDT_deltaE/SP-1237-BSemiExcl-Run4-R22d-v08/*");

//   ch[2].Add("AWG82/ntuples/NN_KM_deltaE/SP-1235-BSemiExcl-Run2-R22d-v08/*");
//   ch[2].Add("AWG82/ntuples/NN_KM_deltaE/SP-1235-BSemiExcl-Run4-R22d-v08/*");
//   ch[2].Add("AWG82/ntuples/NN_KM_deltaE/SP-1237-BSemiExcl-Run2-R22d-v08/*");
//   ch[2].Add("AWG82/ntuples/NN_KM_deltaE/SP-1237-BSemiExcl-Run4-R22d-v08/*");


  TCut MMiss = "candPMiss>.2";
  TCut Q2 = "candQ2>4";
  TCut ee1 = "candType<3&&candEExtra<.2";
  TCut ee2 = "candType==3&&candEExtra<.15";
  TCut ee3 = "candType==4&&candEExtra<.3";
  TCut ee = (ee1||ee2||ee3);
  TCut basic = MMiss+Q2+ee;

  TH1F *histo[2];
  TString names[]={"NN Base (","PID Tweak (","eKMLoose ("};
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
  cuts = "candIsMu==0";
  if(when == "after")
    cuts += basic;
  int colors[]={1,2,4,8};
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
    ch[i].Draw("candPLep>>h(50,0,3.6)",cuts);
    TH1F *h =  (TH1F*)gDirectory->Get("h");
    histo[i] = (TH1F*)h->Clone("histo"+i);
    cout<<"Entries "<<histo[i]->GetEntries()<<endl;
    int entries = histo[i]->GetEntries();
    if(doScale)
      histo[i]->Scale(100./entries);
    if(histo[i]->GetMaximum()>maxi) maxi = histo[i]->GetMaximum();
    names[i] += entries; names[i] += " entries)";
  }
  for(int i=0; i<n;i++){
    histo[i]->SetLineColor(colors[i]);
    histo[i]->SetLineWidth(2);
    if(i==0){
      histo[i]->SetTitle("Electron momentum in the lab frame");
      histo[i]->SetXTitle("p_{#mu} [GeV]");
      histo[i]->SetMaximum(1.05*maxi);
      histo[i]->Draw();
    }else
      histo[i]->Draw("same");
    leg->AddEntry(histo[i],names[i]);
  }
  leg->Draw();
  if(doScale)
    c.SaveAs("mycode/eps/MCType/MCComb_Scaled_"+Dmode+"_"+when+"Cuts.eps");
  else
    c.SaveAs("mycode/eps/candPLepeTwe_"+when+"Cuts.eps");
}


