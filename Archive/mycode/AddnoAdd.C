TString round(double n, double d);

void AddnoAdd(int nvari=0, TString type = ""){
  if(nvari<0||nvari>4){
    cout<<"Variable must be between 0 and 4. Exiting"<<endl;
    return;
  }
  TCut acc = "candPMiss>.2";
  TCut mes = "candMES>5.27";
  TCut bbb1 = "candExtraTracks==0&&candThetaLep>.4&&candThetaLep<2.6";
  TCut ee1 = "candType<3&&candEExtra<.2";
  TCut ee2 = "candType==3&&candEExtra<.15";
  TCut ee3 = "candType==4&&candEExtra<.3";
  TCut ee = (ee1||ee2||ee3);
  TCut basic = acc+bbb1+mes+ee;

  TChain a("ntp1");
  //a.Add("AWG82/ntuples/Add_BDT_KM/Merged/*Run*");
  a.Add("AWG82/ntuples/Add_NN_LH/Merged/*Run*");
  TChain b("ntp1");
  b.Add("AWG82/ntuples/NN/Merged/*Run*");
  TChain s("ntp1");
  s.Add("AWG82/ntuples/semilep/Merged/*");
  a.SetLineColor(1);
  b.SetLineColor(2);
  s.SetLineColor(4);
  TCut tot = "";
  if(type == "ex")
    tot = basic;

  TString vari[] = {"candM2","candPstarLep","candQ2","candEExtra","candMES"};
  TString limits[] = {"(25,-1,2)","(15,0,2.5)","(15,-2,14)","(25,0,2)","(15,5.27,5.29)",
		      "(25,-1,2)","(15,0,2.5)","(15,-1,5)","(25,0,2)","(15,5.27,5.29)"};
  TString titles[] = {"D^{0}: Add/noAdd = ","D*^{0}: Add/noAdd = ","D^{+}: Add/noAdd = ","D*^{+}: Add/noAdd = "};
  TString xtitles[] = {"m^{2}_{miss} [GeV^{2}]","p*_{l} [GeV]","q^{2} [GeV^{2}]","E_{Extra} [GeV]","m_{ES} [GeV]"};
  TCut cuts[] = {"(MCType==1||MCType==3)","(MCType==2||MCType==4)","(MCType==7||MCType==9)","(MCType==8||MCType==10)"};
  TH1F *Add[4];
  TH1F *noAdd[4];
  TH1F *Add2[4];
  TH1F *noAdd2[4];
  TH1F *Sig[4];
  TH1F *Sig2[4];
  TCanvas c("AddNOADD","Add Vs NOADD",800,600);
  c.Divide(2,2);

  for(int i=0; i<4; i++){
    c.cd(i+1);
    TPad *p1 = new TPad("p1","",0,.2,1,1);
    p1->Draw();
    p1->cd();
    TString cand = "candType=="; cand += i+1; //cand+="&&candEExtra<.3&&candPLep>1.5";
    TCut totcand = cuts[i]; totcand += cand; totcand += tot;

    a.Draw(vari[nvari]+">>ADD"+limits[nvari],totcand);
    TH1F *h = (TH1F*)gDirectory->Get("ADD");
    Add[i] = (TH1F*)h->Clone("Add"+i);
    //Add[i]->Sumw2();
    Add[i]->SetMarkerStyle(20);
    Add[i]->SetMarkerSize(.6);
    Add[i]->GetXaxis()->SetTitleOffset(.7);
    Add[i]->GetXaxis()->SetTitleSize(0.05);
    Add[i]->SetXTitle(xtitles[nvari]);
    gStyle->SetOptStat(0);

    b.Draw(vari[nvari]+">>NOADD"+limits[nvari],totcand);
    h = (TH1F*)gDirectory->Get("NOADD");
    noAdd[i] = (TH1F*)h->Clone("MonteCarlo"+i);
    s.Draw(vari[nvari]+">>Sig"+limits[nvari],totcand);
    h = (TH1F*)gDirectory->Get("Sig");
    Sig[i] = (TH1F*)h->Clone("Signal"+i);

    double nAdd = (double)Add[i]->GetEntries(); double nnoAdd = (double)noAdd[i]->GetEntries();
    double nSig = (double)Sig[i]->GetEntries();
    Add[i]->SetTitle(titles[i]+round(nAdd,nnoAdd));
    cout<<titles[i]<<": "<<nAdd<<" Vs "<<nnoAdd<<endl;
    if(nnoAdd)
      noAdd[i]->Scale(nAdd/(nnoAdd));
    if(nSig)
      Sig[i]->Scale(nAdd/(nSig));

    float maxi = Add[i]->GetMaximum();
    if(noAdd[i]->GetMaximum()>maxi) maxi = noAdd[i]->GetMaximum();
    if(Sig[i]->GetMaximum()>maxi) maxi = Sig[i]->GetMaximum();
    Add[i]->SetMaximum(1.15*maxi);
    Add[i]->GetYaxis()->SetLabelSize(0.062);
    Add[i]->GetXaxis()->SetLabelSize(0.06);
    Add[i]->GetYaxis()->SetNdivisions(6+100*2);
    Add[i]->GetXaxis()->SetNdivisions(7+100*2);
    Add[i]->Draw();
    noAdd[i]->Draw("same");
    Sig[i]->Draw("same");
    TLegend *leg= new TLegend(0.58,.7,0.9,0.9);
    leg->SetTextSize(0.046);
    //leg->SetBorderSize(0);
    leg->SetFillColor(0);
    TString dleg = "Add ("; dleg+=(int)nAdd;dleg+=" entries)";
    leg->AddEntry(Add[i],dleg);
    TString mleg = "noAdd ("; mleg+=(int)nnoAdd;mleg+=" entries)";
    leg->AddEntry(noAdd[i],mleg);
    TString sleg = "Signal ("; sleg+=(int)nSig;sleg+=" entries)";
    leg->AddEntry(Sig[i],sleg);
    leg->Draw();
    c.cd(i+1);
    TPad *p2 = new TPad("p2","",0,0,1,.2);
    p2->Draw();
    p2->cd();
    
    Float_t min = Add[i]->GetXaxis()->GetBinLowEdge(Add[i]->GetXaxis()->GetFirst());
    Float_t max = Add[i]->GetXaxis()->GetBinLowEdge(Add[i]->GetXaxis()->GetLast()+1);
    Add2[i] = (TH1F *)Add[i]->Clone();
    Add2[i]->Sumw2();
    noAdd2[i] = (TH1F *)noAdd[i]->Clone();
    noAdd2[i]->Sumw2();
    Add2[i]->Divide(noAdd2[i]);
    TLine l;
    l.SetLineColor(2);
    Add2[i]->GetYaxis()->SetNdivisions(3+100*2);
    Add2[i]->SetTitle("");
    Add2[i]->SetXTitle("");
    Add2[i]->SetMaximum(1.75);
    Add2[i]->SetMinimum(0.25);
    Add2[i]->GetYaxis()->SetLabelSize(0.2);
    Add2[i]->GetXaxis()->SetLabelSize(0.01);
    Add2[i]->SetMarkerSize(.5);
    Add2[i]->Draw("E0 P");
    l.DrawLine(min, 1.0, max, 1.0);
  }
  c.SaveAs("mycode/PlotsAdd/Add"+type+vari[nvari]+".eps");
}


TString round(double n, double d){
  if(d==0) return " - ";
  double b = ((int)(n/d *100+.5));
  b/=100.;
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(b<10){
    if(result.Length() == 1)
      result += ".00";
    if(result.Length() == 3)
      result += "0";
  }
  return result;
}


