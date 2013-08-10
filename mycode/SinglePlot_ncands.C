void SinglePlot_ncands(){
  int n = 3;
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

  ch[2].Add("AWG82/ntuples/Add_BDT_KM/SP-1235-BSemiExclAdd-Run2-R22f-v01/*");
  ch[2].Add("AWG82/ntuples/Add_BDT_KM/SP-1237-BSemiExclAdd-Run2-R22f-v01/*");
  ch[2].Add("AWG82/ntuples/Add_BDT_KM/SP-1235-BSemiExclAdd-Run4-R22f-v01/*");
  ch[2].Add("AWG82/ntuples/Add_BDT_KM/SP-1237-BSemiExclAdd-Run4-R22f-v01/*");

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
  TString names[]={"Base (","Tracks matched (","Tracks and neutrals ("};
  TCut cuts[4];
  for(int i=0; i<n; i++) 
    cuts[i] = basic;
  cuts[1] += "candLepTru==1&&candBntCha==0&&candDntCha==0";
  cuts[2] += "candLepTru==1&&candBntCha==0&&candDntCha==0&&candBntNeu==0&&candDntNeu==0";
  int colors[]={1,2,4,8};
  TCanvas c("ncands","ncands",600,400);
  gStyle->SetOptStat(1111);
  TLegend *leg;
  leg= new TLegend(0.35,1.-.087*n,0.7,0.9);
  leg->SetFillColor(0);
  float maxi = 0;
  int nJPsi=0;
  for(int i=0; i<n;i++){
    ch[2].Draw("ncands>>base(18,0,18)",cuts[i]);
    TH1F *h =  (TH1F*)gDirectory->Get("base");
    histo[i] = (TH1F*)h->Clone("Base"+i);
    cout<<"Entries "<<histo[i]->GetEntries()<<endl;
    int entries = histo[i]->GetEntries();
    if(histo[i]->GetMaximum()>maxi) maxi = histo[i]->GetMaximum();
    names[i] += entries; names[i] += " entries)";
  }
  for(int i=0; i<n;i++){
    histo[i]->SetLineColor(colors[i]);
    histo[i]->SetLineWidth(2);
    if(i==0){
      histo[i]->SetTitle("Number of candidates at best B selection");
      histo[i]->SetXTitle("");
      histo[i]->SetMaximum(1.05*maxi);
      histo[i]->Draw();
    }else
      histo[i]->Draw("same");
    leg->AddEntry(histo[i],names[i]);
  }
  leg->Draw();
  c.SaveAs("mycode/eps/ncandsAdd.eps");
}


