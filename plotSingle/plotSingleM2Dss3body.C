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

void plotSingleM2Dss3body(TString Tag="Dl"){

  TString Dssfile = "AWG82/ntuples/small/semFull_RunAll.root";
  TChain Dss("ntp1");
  Dss.Add(Dssfile);
  TString Dss3file = "AWG82/ntuples/small/Dss3body_RunAll.root";
  TChain Dss3("ntp1");
  Dss3.Add(Dss3file);
  
  TString lund[] = {"10423","10413","425","415"};
  TString name[] = {"D_{1} (D#pi#pi) l #nu (","D_{1} (D*#pi) l #nu (","D_{2} (D^{(}*^{)}#pi) l #nu ("};
  if(Tag.Contains("D1p")){
    lund[0] = "20413";lund[1] = "20423";lund[2] = "10411";lund[3] = "10421";
    name[0] = "D_{1}' (D^{(}*^{)}#pi#pi) l #nu ("; name[1] = "D_{1}' (D*#pi) l #nu (";
    name[2] = "D_{0} (D#pi) l #nu (";
  }
  TCut modeCut[2];
  TString decay[2];
  TString limits[2];
  modeCut[0] = "(MCType==5||MCType==6)";
  modeCut[1] = "(MCType==11||MCType==12)";
  decay[0] = "D** #tau #nu";
  decay[1] = "D** #tau #nu";
  limits[0] = "(35,-0.5,8)";
  limits[1] = "(35,-0.5,8)";
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
    TString cut1 = "((abs(MCD)=="; cut1+=lund[0]; cut1+="||abs(MCD)=="; cut1+=lund[1]; cut1+=")&&MCTaumode==-1)";
    TString cut2 = "((abs(MCD)=="; cut2+=lund[2]; cut2+="||abs(MCD)=="; cut2+=lund[3]; cut2+=")&&MCTaumode==-1)";
    TCut tot1 = cut1; tot1+=MvaAll; tot1 += cuts;  tot1*="wBF";
    TCut tot2 = cut2; tot2+=MvaAll; tot2 += cuts;  tot2*="wBF";
    TString vari = "candM2>>h"; vari+=limits[a]; 
    if(Tag.Contains("pi0")){
      vari = "mm2pi0>>h"; vari+=limits[a]; 
    }
    Dss3.Draw(vari,tot1);
    TH1F *h = (TH1F*)gDirectory->Get("h");
    mmiss[i][0] = (TH1F*)h->Clone("cand"+i);
    Dss.Draw(vari,tot1);
    TH1F *h = (TH1F*)gDirectory->Get("h");
    mmiss[i][1] = (TH1F*)h->Clone("candMva"+i);
    Dss.Draw(vari,tot2);
    TH1F *h = (TH1F*)gDirectory->Get("h");
    mmiss[i][2] = (TH1F*)h->Clone("candMva"+i);
    mmiss[i][2]->SetLineColor(4);
    double ngen = mmiss[i][0]->Integral();
    double nsig = mmiss[i][1]->Integral();
    double nold = mmiss[i][2]->Integral();
    if(nsig!=0) {
      if(ngen!=0) mmiss[i][1]->Scale(ngen/nsig);
      else if(nold!=0) mmiss[i][2]->Scale(nsig/nold);
    }      
    if(nold!=0 && ngen!=0) mmiss[i][2]->Scale(ngen/nold);

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
      leg = new TLegend(0.45,.62,0.9,0.9);
    leg->SetTextSize(0.068);
    leg->SetFillColor(0);
    TString legtag = name[0]; legtag+=round(ngen,0); legtag+=")";
    leg->AddEntry(mmiss[i][0],legtag);
    legtag = name[1]; legtag+=round(nsig,0); legtag+=")";
    leg->AddEntry(mmiss[i][1],legtag);
    legtag = name[2]; legtag+=round(nold,0); legtag+=")";
    leg->AddEntry(mmiss[i][2],legtag);
    leg->Draw();
  }
  c.SaveAs("babar_code/eps/Dss3bodyMmiss_"+Tag+".eps");
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
