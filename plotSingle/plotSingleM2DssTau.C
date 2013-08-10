//---------------------------------------------------------------------------------
// Description:
//      Plots the a single variable for the 4 channels
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      09/06/30 manuelf -- Upgrade
//      09/05/31 manuelf -- Creation
//---------------------------------------------------------------------------------

#include "mycode/cuts.hh"
TString round(double n, int e, double d=1.);

void plotSingleM2DssTau(TString Tag="Dl"){

  TString Dssfile = "AWG82/ntuples/small/semFull_RunAll.root";
  TChain Dss("ntp1");
  Dss.Add(Dssfile);
  TChain sig("ntp1");
  sig.Add("AWG82/ntuples/small/sigFull_RunAll.root");
  sig.SetLineColor(2);

  TCut modeCut[2];
  TString decay[2];
  TString limits[2];
  modeCut[0] = "(MCType==5||MCType==6)";
  modeCut[1] = "(MCType==11||MCType==12)";
  decay[0] = "D** #tau #nu";
  decay[1] = "D** #tau #nu";
  limits[0] = "(50,-0.5,9)";
  limits[1] = "(50,-0.5,9)";
  TString Dnames[] = {"D^{0}","D*^{0}","D^{+}","D*^{+}"};
  TH1F *mmiss[4][3];
  TLatex *label = new TLatex();
  label->SetNDC(kTRUE);
  label->SetTextSize(0.08);
  TCanvas c("cMvaDl","Plot for the different channels",1000,600);
  c.Divide(2,2);
  gStyle->SetOptStat(0);
  int a=0;
  for(int i=0; i<4; i++){
    a = i%2; 
    c.cd(i+1);
    TString cuts = "candType=="; cuts += i+1; 
    TCut tot1 = "(MCType>12&&MCTaumode>-1)"; tot1+=MvaAll; tot1 += cuts;  tot1*="wBF";
    TCut tot2 = modeCut[i>1]; tot2+=MvaAll; tot2 += cuts;  tot2*="wBF";
    TCut tot3 = "(MCType>12&&MCTaumode==-1)"; tot3+=MvaAll; tot3 += cuts;  tot3*="wBF";
    TString vari = "candM2>>h"; vari+=limits[a]; 
    if(Tag.Contains("pi0")){
      vari = "mm2pi0>>h"; vari+=limits[a]; 
    }
    Dss.Draw(vari,tot1);
    TH1F *h = (TH1F*)gDirectory->Get("h");
    mmiss[i][0] = (TH1F*)h->Clone("cand"+i);
    sig.Draw(vari,tot2);
    TH1F *h = (TH1F*)gDirectory->Get("h");
    mmiss[i][1] = (TH1F*)h->Clone("candMva"+i);
    Dss.Draw(vari,tot3);
    TH1F *h = (TH1F*)gDirectory->Get("h");
    mmiss[i][2] = (TH1F*)h->Clone("candMva"+i);
    mmiss[i][2]->SetLineColor(4);
    double ngen = mmiss[i][0]->Integral();
    double nsig = mmiss[i][1]->Integral();
    double nold = mmiss[i][2]->Integral();
    if(nsig) {
      if(ngen) mmiss[i][1]->Scale(ngen/nsig);
      else if(nold) mmiss[i][2]->Scale(nsig/nold);
    }      
    if(nold && ngen) mmiss[i][2]->Scale(ngen/nold);

    double maxi = mmiss[i][0]->GetMaximum();
    if(maxi<mmiss[i][1]->GetMaximum()) maxi = mmiss[i][1]->GetMaximum();
    if(maxi<mmiss[i][2]->GetMaximum()) maxi = mmiss[i][2]->GetMaximum();

    mmiss[i][0]->SetMaximum(1.15*maxi);
    mmiss[i][0]->GetXaxis()->SetTitleOffset(.7);
    mmiss[i][0]->GetXaxis()->SetTitleSize(0.059);
    mmiss[i][0]->GetYaxis()->SetLabelSize(0.062);
    mmiss[i][0]->GetXaxis()->SetLabelSize(0.056);
    mmiss[i][0]->GetYaxis()->SetNdivisions(5+100*2);
    mmiss[i][0]->GetXaxis()->SetNdivisions(6+100*2);
    mmiss[i][0]->SetLineWidth(2);
    mmiss[i][0]->SetTitle("");
    mmiss[i][0]->SetXTitle("m^{2}_{miss} [GeV^{2}]");
    mmiss[i][0]->Draw();
    mmiss[i][1]->SetLineWidth(2);
    mmiss[i][1]->SetLineColor(2);
    mmiss[i][2]->Draw("same");
    mmiss[i][0]->Draw("same");
    mmiss[i][1]->Draw("same");
    TString tag = Dnames[i]; if(Tag.Contains("pi0")) tag += "#pi^{0}";
    tag+= " channel"; 
    label->DrawLatex(.11,0.92,tag);
    TLegend *leg;
    if(Tag=="pi0Ds0"||Tag=="pi0D0")
      leg = new TLegend(0.1,.66,0.43,0.9);
    else
      leg = new TLegend(0.53,.62,0.9,0.9);
    leg->SetTextSize(0.068);
    leg->SetFillColor(0);
    TString legtag = "D** #tau #nu ("; legtag+=round(ngen,0); legtag+=")";
    leg->AddEntry(mmiss[i][0],legtag);
    legtag = "D^{(*)} #tau #nu ("; legtag+=round(nsig,0); legtag+=")";
    leg->AddEntry(mmiss[i][1],legtag);
    legtag = "D** l #nu ("; legtag+=round(nold,0); legtag+=")";
    leg->AddEntry(mmiss[i][2],legtag);
    leg->Draw();
  }
  c.SaveAs("babar_code/eps/DssTauMmiss_"+Tag+".eps");
}


TString round(double n, int e, double d){
  if(d==0) return " - ";
  double b = (int)(n/d*pow(10,e)+0.5);
  b /= pow(10,e);
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(!result.Contains(".")&&e!=0) result += ".";
  
  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<e-afterdot.Length(); i++)
    result += "0";
  return result;
}
