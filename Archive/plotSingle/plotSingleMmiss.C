//---------------------------------------------------------------------------------
// Description:
//      Plots the mmiss resolution
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      09/05/31 manuelf -- Creation
//---------------------------------------------------------------------------------
TString round(double n, double d);

void plotSingle(TString Tag="Signal"){
  TChain b("ntp1");
   b.Add("AWG82/ntuples/Dss/Add_Truth30_Run5.root");
   b.Add("AWG82/ntuples/Dss/Add_Truth30_Run6.root");
   //  b.Add("AWG82/ntuples/Dss/uds30_Run6.root");
   //  b.Add("AWG82/ntuples/Dss/ccbar30_Run6.root");

  TCut cosT = "abs(candCosT)<.9";
  TCut dssacc = "mm2pi0>-4&&mm2pi0<12&&candPstarLep>0&&candPstarLep<2.4&&candMES>5.2&&candMES<5.3&&pmisspi0>.2";
  TCut Mva2 = "(candMvaDssComb>0.2&&candMvaDssDl>0.1)";
  TCut dssMva = Mva2+dssacc;

  TH1F *mmiss[4];
  TString Dnames[] = {"D^{0}","D*^{0}","D^{+}","D*^{+}"};
  TLatex *label = new TLatex();
  label->SetNDC(kTRUE);
  label->SetTextSize(0.06);
  TCanvas c("cMvaDl","Optimization MvaDl",800,600);
  c.Divide(2,2);
  gStyle->SetOptStat(1111110);
  for(int i=0; i<4; i++){
    c.cd(i+1);
    TString cuts = "(candType<3&&MCType>0&&MCType<7||candType>2&&MCType>6&&MCType<13)&&candType=="; cuts += i+1;
    TCut tot = dssMva; tot += cuts;
    b.Draw("mm2pi0>>h(60,-4.,8.)",tot);
    TH1F *h = (TH1F*)gDirectory->Get("h");
    mmiss[i] = (TH1F*)h->Clone("cand"+i);
    mmiss[i]->SetLineWidth(2);
    mmiss[i]->SetTitle("");
    mmiss[i]->SetXTitle("m^{2}_{miss} [GeV^{2}]");
    mmiss[i]->GetXaxis()->SetTitleSize(0.048);
    mmiss[i]->GetXaxis()->SetTitleOffset(.9);
    mmiss[i]->Draw();
    TString tag = Dnames[i]; tag+= ": "; tag+= Tag; 
    label->DrawLatex(.11,0.93,tag);
  }
  c.SaveAs("babar_code/eps/DssMmiss_"+Tag+".eps");
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
