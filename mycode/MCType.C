void MCType(TString Dmode = "D0", TString when = "after", bool doScale = false){

  TCut MMiss = "candPMiss>.2";
  TCut Q2 = "candQ2>4";
  TCut ee1 = "candType<3&&candEExtra<.2";
  TCut ee2 = "candType==3&&candEExtra<.15";
  TCut ee3 = "candType==4&&candEExtra<.3";
  TCut ee = (ee1||ee2||ee3);
  TCut basic = MMiss+Q2+ee;
  // Loading files
  int n = 2;
  TChain ch[10];
  int colors[]={1,2,4,8,3,28,6,7,10,12};
  for(int i=0; i<10; i++)  ch[i] = new TChain("ntp1");
  for(int i=0; i<10; i++)   ch[i].SetLineColor(colors[i]);

  ch[0].Add("AWG82/ntuples/NN/Merged/*Run2*");
  ch[0].Add("AWG82/ntuples/NN/Merged/*Run4*");
  ch[0].Add("AWG82/ntuples/NN/Merged/*Run3*");

  ch[1].Add("AWG82/ntuples/Add_BDT_KM/Merged/*Run2*");
  ch[1].Add("AWG82/ntuples/Add_BDT_KM/Merged/*Run4*");
  ch[1].Add("AWG82/ntuples/Add_BDT_KM/Merged/*Run3*");

  ch[2].Add("AWG82/ntuples/data/Merged/*Run2*");
  ch[2].Add("AWG82/ntuples/data/Merged/*Run4*");
  ch[2].Add("AWG82/ntuples/data/Merged/*Run3*");

  ch[2].Add("AWG82/ntuples/Add_BDT_KM_Truth55/Merged/*Run2*");
  ch[2].Add("AWG82/ntuples/Add_BDT_KM_Truth55/Merged/*Run4*");
  ch[2].Add("AWG82/ntuples/Add_BDT_KM_Truth55/Merged/*Run3*");

  ch[3].Add("AWG82/ntuples/Add_BDT_KM_Q2/Merged/*Run2*");
  ch[3].Add("AWG82/ntuples/Add_BDT_KM_Q2/Merged/*Run4*");
  ch[3].Add("AWG82/ntuples/Add_BDT_KM_Q2/Merged/*Run3*");

//   ch[3].Add("AWG82/ntuples/Add_NN_LH_Truth55/Merged/*Run2*");
//   ch[3].Add("AWG82/ntuples/Add_NN_LH_Truth55/Merged/*Run4*");
//   ch[3].Add("AWG82/ntuples/Add_NN_LH_Truth55/Merged/*Run3*");

//   ch[4].Add("AWG82/ntuples/NN_KM/Merged/*Run2*");
//   ch[4].Add("AWG82/ntuples/NN_KM/Merged/*Run4*");
//   ch[4].Add("AWG82/ntuples/NN_KM/Merged/*Run3*");

//   ch[5].Add("AWG82/ntuples/BDT_LH/Merged/*Run2*");
//   ch[5].Add("AWG82/ntuples/BDT_LH/Merged/*Run4*");
//   ch[5].Add("AWG82/ntuples/BDT_LH/Merged/*Run3*");

//   ch[2].Add("AWG82/ntuples/Add_BDT_KM/SP-1235-BSemiExclAdd-Run2-R22f-v01/*");
//   ch[2].Add("AWG82/ntuples/Add_BDT_KM/SP-1237-BSemiExclAdd-Run2-R22f-v01/*");
//   ch[2].Add("AWG82/ntuples/Add_BDT_KM/SP-1235-BSemiExclAdd-Run4-R22f-v01/*");
//   ch[2].Add("AWG82/ntuples/Add_BDT_KM/SP-1237-BSemiExclAdd-Run4-R22f-v01/*");
//   ch[2].Add("AWG82/ntuples/Add_BDT_KM/SP-1235-BSemiExclAdd-Run3-R22f-v01/*");
//   ch[2].Add("AWG82/ntuples/Add_BDT_KM/SP-1237-BSemiExclAdd-Run3-R22f-v01/*");



  TH1F *histo[6];
  TString BinName[]={"D^{0}e","D*^{0}e","D^{0}#mu","D*^{0}#mu",
		     "D^{0}#tau","D*^{0}#tau","D^{+}e","D*^{+}e","D^{+}#mu",
		     "D*^{+}#mu","D^{+}#tau","D*^{+}#tau","D^{(}*^{)}#pi l",
		     "D**l","Bkg"};
  TString names[]={"NN Base (","BDT (","Add: restricted ("};
  TCut cuts="";
  TString Dname="";
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
//   cuts += "candIsMu==1";
  if(when == "after")
    cuts += basic;
  int colors[]={1,2,4,8};
  TCanvas c("Yields","Yields",600,400);
  gStyle->SetOptStat(0);
  TLegend *leg;
  if(when == "after" && (Dmode == "D0" || Dmode == "Dstar0"))
    leg= new TLegend(0.61,1.-.089*n,0.9,0.9);
  else
    leg = new TLegend(0.1,1.-.089*n,0.43,0.9);
  leg->SetFillColor(0);
  float maxi = 0, first;
  for(int i=0; i<n;i++){
    ch[i].Draw("MCType>>h",cuts);
    histo[i] = (TH1F*)h->Clone("histo"+i);
    int entries = histo[i]->GetEntries();
    if(doScale)
      histo[i]->Scale(100./entries);
    if(histo[i]->GetMaximum()>maxi) maxi = histo[i]->GetMaximum();
    first = histo[i]->GetBinContent(1);
    for(int j=2;j<16;j++)
      histo[i]->SetBinContent(j-1,histo[i]->GetBinContent(j));
    histo[i]->SetBinContent(15,first);
    names[i] += entries; names[i] += " entries)";
 }
  for(int i=0; i<n;i++){
    histo[i]->SetLineColor(colors[i]);
    histo[i]->SetLineWidth(2);
    if(i==0){
      histo[i]->SetTitle("MC types for the "+Dname+" "+when+" cuts");
      histo[i]->SetMaximum(1.05*maxi);
      histo[i]->Draw();
      for(int j=1;j<16;j++)
	histo[i]->GetXaxis()->SetBinLabel(j,BinName[j-1]);
    }else
      histo[i]->Draw("same");
    leg->AddEntry(histo[i],names[i]);
  }
  leg->Draw();
  if(doScale)
    c.SaveAs("mycode/eps/MCType/MCType_Scaled_"+Dmode+"_"+when+"Cuts.eps");
  else
    c.SaveAs("mycode/eps/MCType/MCType2_"+Dmode+"_"+when+"Cuts.eps");
}


