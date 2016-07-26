void MCplots(){
  // Loading files
  TChain ch("ntp1");
  ch.Add("AWG82/ntuples/SP-1235-BSemiExcl-Run2-R22d-v06/*");
  ch.Add("AWG82/ntuples/SP-1235-BSemiExcl-Run4-R22d-v06/*");
  ch.Add("AWG82/ntuples/SP-1237-BSemiExcl-Run2-R22d-v06/*");
  ch.Add("AWG82/ntuples/SP-1237-BSemiExcl-Run4-R22d-v06/*");
  TCut MMiss = "candPMiss>.2";
  TCut Q2 = "candQ2>4";
  TCut ee1 = "candType<3&&candEExtra<.2";
  TCut ee2 = "candType==3&&candEExtra<.15";
  TCut ee3 = "candType==4&&candEExtra<.3";
  TCut ee = (ee1||ee2||ee3);

  TH1F *histo[6];
  TString BinName[]={"D^{0}e","D*^{0}e","D^{0}#mu","D*^{0}#mu",
		     "D^{0}#tau","D*^{0}#tau","D^{+}e","D*^{+}e","D^{+}#mu",
		     "D*^{0}#mu","D^{+}#tau","D*^{+}#tau","D^{(}*^{)}#pi l",
		     "D**l","Bkg"};
  TCut D0cut = "candType<3"+MMiss+Q2+ee;
  TCut D0pcut = "candType==4&&candDstarType==1";
  TCut Dpluscut = "candType==3||(candType==4&&candDstarType==5)"+MMiss+Q2+ee;

  int n = 4;
//   TCut cuts[]={D0cut+"candDType<5",D0cut+"candDType==5",
// 	       D0cut+"candDType==8"};
//   TString names[]={"Mazur (","K_{S}^{0} #pi^{+} #pi^{-} #pi^{0} (","#pi^{+} #pi^{-} ("};

//   TCut cuts[]={D0cut+"candDType<5",D0cut+"candDType==7",
// 	       D0cut+"candDType==9"};
//   TString names[]={"Mazur (","K_{S}^{0} #pi^{0} (","K^{+} K^{-} ("};
  TCut cuts[]={Dpluscut+"(candDType<4||candDType==6)",Dpluscut+"candDType==4",
	       Dpluscut+"candDType==5",Dpluscut+"candDType==7"};
  TString names[]={"Mazur (","K_{S}^{0} #pi^{+} #pi^{0} (","K_{S}^{0} #pi^{+} #pi^{+} #pi^{-} (","K_{S}^{0} K^{+} ("};


  int colors[]={1,2,4,8};
  TCanvas c("Yields","Yields",600,400);
  gStyle->SetOptStat(0);
  TLegend *leg= new TLegend(0.1,1.-.089*n,0.42,0.9);
//   TLegend *leg= new TLegend(0.62,1.-.09*n,0.9,0.9);
  leg->SetFillColor(0);
  float maxi = 0, first;
  for(int i=0; i<n;i++){
    ch.Draw("MCType>>h",cuts[i]);
    histo[i] = (TH1F*)h->Clone("histo"+i);
    int entries = histo[i]->GetEntries();
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
      histo[i]->SetTitle("% MC types for the different D^{+} modes");
      histo[i]->SetMaximum(1.05*maxi);
      histo[i]->Draw();
      for(int j=1;j<16;j++)
	histo[i]->GetXaxis()->SetBinLabel(j,BinName[j-1]);
    }else
      histo[i]->Draw("same");
    leg->AddEntry(histo[i],names[i]);
  }
  leg->Draw();
  c.SaveAs("mycode/eps/MCType_Dplus.eps");
}


