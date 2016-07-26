#include "mycode/cuts.hh"

TString round(double n, double d);

void CSampleconti(int nvari=0, TString type = "pl", TString run = "56", TString cuts="Mva"){
  TString isMaz = "";
  if(nvari<0||nvari>8){
    cout<<"Variable must be between 0 and 8. Exiting"<<endl;
    return;
  }
  basic = PMiss;
  if(cuts=="Mva") basic += Mva;
  else if (nvari!=3) basic += ee;
  if(cuts.Contains("cosT"))
    basic += cosT;

  if(cuts.Contains("All"))
    dss += mpi0+dssee;
  else if(cuts.Contains("Mva"))
    dss += dssMva;
  if(cuts.Contains("cosT"))
    dss += cosT;


  TChain uds("ntp1");
  TChain ccbar("ntp1");
  uds.Add("AWG82/ntuples/small/uds30_Run12356.root");
  ccbar.Add("AWG82/ntuples/small/ccbar30_Run12356.root");
  TString dfile = "AWG82/ntuples/small/Adddata30_Run"; dfile += run; dfile +=".root";
  TString mfile = "AWG82/ntuples/small/Add_Truth30Rest_Run"; mfile += run; mfile +=".root";
  TChain d("ntp1");
  d.Add(dfile);
  double fracRest = 1.;
  if(mfile.Contains("Rest")) fracRest = 2/3.;
  TChain m("ntp1");
  m.SetLineColor(2);
  m.Add(mfile);
  cout<<"Data: "<<dfile<<" with entries: "<<d.GetEntries()<<endl;   
  cout<<"MC: "<<mfile<<" with entries: "<<m.GetEntries()<<endl;   

  TCut tot;
  int nlim = 0;
  if(type == "pl")
    tot = basic+"candPstarLep>1.5";
  else if(type == "q2"){
    tot = basic+"candQ2<5";
    nlim = 9;
  } else if(type == "dss"){
    tot = dss;
    nlim = 18;
  }else
    tot = basic+"candPstarLep>1.5&&candQ2>4"; 
   
  TCut CombCut[] = {"MCType==0||MCType>6","MCType>=0&&MCType<7||MCType>12"};
  TCut DssCombCut = "MCType==0";
  TCut NormCut[] = {"MCType>0&&MCType<5","MCType>6&&MCType<11"};
  TCut DssDlCut = "MCType>0&&MCType<13";
  TCut SigCut[] = {"MCType>4&&MCType<7","MCType>10&&MCType<13"};
  TCut DssSigCut = "MCType>12";

  TString vari[] = {"candM2","candPstarLep","candQ2","candEExtra","candMES","candMvaDl",
		    "candRejectedTracks","candRejectedPhotons","abs(candCosT)"};
  if(type=="dss"){
    vari[0] = "mm2pi0";
    vari[3] = "eextrapi0";
  }
  TString Namevari[] = {"M2","PstarLep","Q2","EExtra","MES","MvaDl",
		    "RejectedTracks","RejectedPhotons","CosT"};
  if(nvari==3) vari[3] = "eextrapi0";
  TString limits[] = {"(15,-.5,1.5)","(15,1.5,2.5)","(15,-2,14)","(25,0,2)","(15,5.27,5.29)","(15,-0.8,0.8)","(5,0,5)","(12,0,12)","(15,0,1)",
		      "(15,-.5,1.5)","(15,0,2.5)","(15,-1,5)","(25,0,.5)","(15,5.27,5.29)","(15,-0.8,0.8)","(5,0,5)","(12,0,12)","(15,0,1)",
		      "(20,-2.5,5)","(15,0,2.4)","(15,-1,5)","(15,0,.5)","(15,5.27,5.29)","(15,-0.,0.8)","(5,0,5)","(12,0,12)","(15,0,1)"};
  TString titles[] = {"D^{0}: MC/data = ","D*^{0}: MC/data = ","D^{+}: MC/data = ","D*^{+}: MC/data = "};
  TString xtitles[] = {"m^{2}_{miss} [GeV^{2}]","p*_{l} [GeV]","q^{2} [GeV^{2}]",
		       "E_{Extra} [GeV]","m_{ES} [GeV]","Dlnu MVA","Rejected tracks", "Rejected photons", "abs[cos(#theta_{T})]"};

  //double ratios[] = {3.31255, 3.06587, 2.82047, 3.03722, 3.32314, 1.60632};
  double nccbar[] = {58900000., 168844000., 83974000., 0., 366758000., 104778000.}; // Run 4 is 252830000
  double nuds[] = {47180000., 130858000., 66892000., 0., 317846000., 84414000}; // Run 4 is 213380000
  double MCBs[] = {36968000+37200000., 103498000+103124000., 50556000+49766000., 167994000+167466000., 244322000+244812000., 68148000+68016000.};
  double dataBs[] = {22389980.4, 67394307.5, 35569248.8, 110449802.7, 147190396.5, 84767412.6};
  double nMCB = 0, ndataBs = 0;
  double totuds = 0, totccbar = 0;
  for(int i=1; i<7; i++){
    totuds += nuds[i-1];
    totccbar += nccbar[i-1];
    TString irun = ""; irun += i;
    if(run.Contains(irun)){
      nMCB += MCBs[i-1]; ndataBs += dataBs[i-1];
    }
  }
  double wuds = nMCB/totuds*2.09/1.05*fracRest; 
  double wccbar = nMCB/totccbar*1.3/1.05*fracRest; 
  cout<<"Weight uds: "<<wuds<<"; Weight ccbar: "<<wccbar<<endl;
  double MCDlumi = 1;
  if(ndataBs) MCDlumi = nMCB/ndataBs;

  TLatex *label = new TLatex();
  label->SetNDC(kTRUE);
  TH1F *data[4];
  TH1F *mc[4];  
  TH1F *hStack[4][4];
  TH1F *huds[4];
  TH1F *hccbar[4];
  TH1F *data2[4];
  TH1F *mc2[4];
  TH1F *mcDs[4];
  THStack hs[4];
  int colors[] = {5,34,30,2}; 
  TString legTag[] = {"Cont (","Bkg (","Dl#nu (","D#tau#nu (","Cont (","Comb (","Dl#nu (","D** ("};
  int legN[4];
  for (int i=0; i<4; i++) hs[i] = new THStack("hs"+i,"Stack for MC");
  TCanvas c("dataMC","data Vs MC",800,600);
  c.Divide(2,2);

  for(int i=0; i<4; i++){
    c.cd(i+1);
    TPad *p1 = new TPad("p1","",0,.25,1,1);
    p1->Draw();
    p1->cd();
    TString cand = "candType=="; cand+=i+1;
    TCut totcand = tot; totcand += cand;
    d.Draw(vari[nvari]+">>dat"+limits[nvari+nlim],totcand);
    TH1F *h = (TH1F*)gDirectory->Get("dat");
    data[i] = (TH1F*)h->Clone("data"+i);
    data[i]->Sumw2();
    data[i]->SetMarkerStyle(20);
    data[i]->SetMarkerSize(.6);
    gStyle->SetOptStat(0);

    m.Draw(vari[nvari]+">>MC"+limits[nvari+nlim],totcand);
    h = (TH1F*)gDirectory->Get("MC");
    mc[i] = (TH1F*)h->Clone("MonteCarlo"+i);
    int a = 0; if(i>1) a=1;
    if(type=="dss")
      m.Draw(vari[nvari]+">>MCbkg1"+limits[nvari+nlim],totcand+DssCombCut);
    else 
      m.Draw(vari[nvari]+">>MCbkg1"+limits[nvari+nlim],totcand+CombCut[a]);
    h = (TH1F*)gDirectory->Get("MCbkg1");
    hStack[1][i] = (TH1F*)h->Clone("Bkg1"+i);
    if(type=="dss")
      m.Draw(vari[nvari]+">>MCbkg2"+limits[nvari+nlim],totcand+DssDlCut);
    else 
      m.Draw(vari[nvari]+">>MCbkg2"+limits[nvari+nlim],totcand+NormCut[a]);
    h = (TH1F*)gDirectory->Get("MCbkg2");
    hStack[2][i] = (TH1F*)h->Clone("Bkg2"+i);
    if(type=="dss")
      m.Draw(vari[nvari]+">>MCsig"+limits[nvari+nlim],totcand+DssSigCut);
    else 
      m.Draw(vari[nvari]+">>MCsig"+limits[nvari+nlim],totcand+SigCut[a]);
    h = (TH1F*)gDirectory->Get("MCsig");
    hStack[3][i] = (TH1F*)h->Clone("Sig"+i);
    uds.Draw(vari[nvari]+">>uds"+limits[nvari+nlim],totcand);
    h = (TH1F*)gDirectory->Get("uds");
    huds[i] = (TH1F*)h->Clone("Uds"+i);
    ccbar.Draw(vari[nvari]+">>ccbar"+limits[nvari+nlim],totcand);
    h = (TH1F*)gDirectory->Get("ccbar");
    hccbar[i] = (TH1F*)h->Clone("Ccbar"+i);

    mc[i]->Add(huds[i],wuds);
    mc[i]->Add(hccbar[i],wccbar);
    double ndata = (double)data[i]->GetEntries(); double nmc = (double)mc[i]->GetEntries();
    data[i]->SetTitle("");
    cout<<titles[i]<<": "<<ndata<<" Vs "<<nmc<<"; uds "<<huds[i]->GetEntries();
    cout<<" Vs "<<hccbar[i]->GetEntries()<<endl;
    if(nmc)
      mc[i]->Scale(ndata/nmc);

    legN[0] = wuds*huds[i]->GetEntries()+wccbar*hccbar[i]->GetEntries();
    huds[i]->Scale(wuds);
    huds[i]->SetLineStyle(2);
    hccbar[i]->Scale(wccbar);
    huds[i]->SetFillColor(colors[0]);
    hccbar[i]->SetFillColor(colors[0]);
    hccbar[i]->SetLineColor(colors[0]);
    hStack[0][i] = (TH1F*)huds[i]->Clone("Cont"+i);
    hStack[0][i]->Add(hccbar[i]);
    huds[i]->Scale(ndata/nmc);
    hccbar[i]->Scale(ndata/nmc);
    hs[i].Add(huds[i]);
    hs[i].Add(hccbar[i]);
    for(int nh=1; nh<4; nh++){
      legN[nh] = hStack[nh][i]->GetEntries();
      hStack[nh][i]->Scale(ndata/nmc);
      hStack[nh][i]->SetFillColor(colors[nh]);
      hStack[nh][i]->SetLineColor(colors[nh]);
      hs[i].Add(hStack[nh][i]);
    }

    float maxi = data[i]->GetMaximum();
    if(mc[i]->GetMaximum()>maxi) maxi = mc[i]->GetMaximum();
    data[i]->SetMaximum(1.15*maxi);
    data[i]->SetMinimum(0);
    hs[i].Draw("");
    hs[i].GetYaxis()->SetLabelSize(0.062);
    hs[i].GetXaxis()->SetLabelSize(0.08);
    hs[i].GetYaxis()->SetNdivisions(6+100*2);
    hs[i].GetXaxis()->SetNdivisions(7+100*2);
    data[i]->Draw("a");
    mc[i]->Draw("same");
    hs[i].Draw("same");
    data[i]->Draw("same");
    TString RatioTitle = titles[i]; RatioTitle += round(nmc/MCDlumi/fracRest,ndata);
    label->SetTextSize(0.085);
    label->DrawLatex(0.02,0.92,RatioTitle);
    TLegend *leg;
    if(nvari==8){
      if(type=="dss")
	leg = new TLegend(0.1,.7,0.5,0.9);
      else
	leg = new TLegend(0.5,.1,0.9,0.3);
    }
    else
      leg = new TLegend(0.68,.56,0.9,0.9);
    leg->SetTextSize(0.06);
    //leg->SetBorderSize(0);
    leg->SetFillColor(0);
    TString dleg = "data ("; dleg+=(int)ndata;dleg+=")";
    leg->AddEntry(data[i],dleg);
    int offset = 0; if(type=="dss") offset = 4;
    for(nh=3; nh>=0; nh--){
      TString mleg = legTag[nh+offset]; mleg+=legN[nh];mleg+=")";
      leg->AddEntry(hStack[nh][i],mleg);
    }
    leg->Draw();
    c.cd(i+1);
    TPad *p2 = new TPad("p2","",0,0,1,.2);
    p2->Draw();
    p2->cd();
    
    Float_t min = data[i]->GetXaxis()->GetBinLowEdge(data[i]->GetXaxis()->GetFirst());
    Float_t max = data[i]->GetXaxis()->GetBinLowEdge(data[i]->GetXaxis()->GetLast()+1);
    data2[i] = (TH1F *)data[i]->Clone();
    mc2[i] = (TH1F *)mc[i]->Clone();
    mc2[i]->Sumw2();
    data2[i]->Divide(mc2[i]);
    TLine l;
    l.SetLineColor(2);
    data2[i]->GetYaxis()->SetNdivisions(3+100*2);
    data2[i]->SetTitle("");
    data2[i]->SetXTitle("");
    data2[i]->SetMaximum(1.5);
    data2[i]->SetMinimum(0.5);
    data2[i]->GetYaxis()->SetLabelSize(0.25);
    data2[i]->GetXaxis()->SetLabelSize(0.01);
    data2[i]->SetMarkerSize(.5);
    data2[i]->Draw("E0 P");
    l.DrawLine(min, 1.0, max, 1.0);
    c.cd(i+1);
    label->SetTextSize(0.06);
    label->DrawLatex(0.68,0.21,xtitles[nvari]);
    
  }
  c.SaveAs("mycode/PlotsCS/conti"+isMaz+type+Namevari[nvari]+"Run"+run+cuts+".eps");
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


