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

void plotMmiss(bool Truth = false, TString cut = "0"){
  TChain b("ntp1");
  b.Add("AWG82/ntuples/Add_Truth30/Merged/*Run5*");
  b.Add("AWG82/ntuples/Add_Truth30/Merged/*Run6*");

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
    TString cuts = "candMvaDl>0."; cuts += cut; 
    cuts += "&&candType==", cuts += i+1;
    b.Draw("candM2Tru-candM2>>h(60,-1.,1.)",cuts);
    TH1F *h = (TH1F*)gDirectory->Get("h");
    mmiss[i] = (TH1F*)h->Clone("cand"+i);
    mmiss[i]->SetLineWidth(2);
    mmiss[i]->SetTitle("");
    mmiss[i]->SetXTitle("m^{2}_{miss,true} - m^{2}_{miss,reco}  [GeV^{2}]");
    mmiss[i]->GetXaxis()->SetTitleSize(0.048);
    mmiss[i]->GetXaxis()->SetTitleOffset(.9);
    mmiss[i]->Draw();
    TString tag = Dnames[i]; tag+=": mean "; tag+=round(mmiss[i]->GetMean(),1.);
    tag+=", RMS "; tag+=round(mmiss[i]->GetRMS(),1.);
    label->DrawLatex(.11,0.93,tag);
  }
  c.SaveAs("babar_code/eps/Mmiss_"+cut+".eps");
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
