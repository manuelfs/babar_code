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
TString round(double n, double d);

void plotSingleDss(TString Tag="D1", TString run = "1234"){

  TString mfile = "AWG82/ntuples/small/Add_Truth30Rest_Run"; mfile += run; mfile +=".root";
  TChain b("ntp1");
  b.Add(mfile);

  TCut type = "MCType == 14 && abs(MCD) == 10423"; TCut type2 = "MCType == 14 && abs(MCD) == 10413"; 
  TString typeTitle = "D_{1}";
  if (Tag=="D2"){
    type = "MCType == 14 && abs(MCD) == 425"; type2 = "MCType == 14 && abs(MCD) == 415"; 
    typeTitle = "D_{2}";
  }
  if (Tag=="D0"){
    type = "MCType == 14 && abs(MCD) == 10421"; type2 = "MCType == 14 && abs(MCD) == 10411"; 
    typeTitle = "D_{0}";
  }
  if (Tag=="D1p"){
    type = "MCType == 14 && abs(MCD) == 20423"; type2 = "MCType == 14 && abs(MCD) == 20413"; 
    typeTitle = "D_{1p}";
  }
  if (Tag=="NR"){
    type = "MCType == 13 && (MCDssmode < 3 || MCDssmode > 5)"; type2 = "MCType == 14 && (MCDssmode >= 3 && MCDssmode <= 5)"; 
    typeTitle = "NR";
  }
  TH1F *mmiss[4][2];
  TString Dnames[] = {"D^{0}","D*^{0}","D^{+}","D*^{+}"};
  TLatex *label = new TLatex();
  label->SetNDC(kTRUE);
  label->SetTextSize(0.06);
  TCanvas c("cMvaDl","Plot for the different channels",1000,600);
  c.Divide(2,2);
  gStyle->SetOptStat(0);
  for(int i=0; i<4; i++){
    c.cd(i+1);
    TString cuts = "candType=="; cuts += i+1; 
    TCut tot = dss+dssMva+type; tot += cuts;
    b.Draw("mm2pi0>>h(20,-2.5,5.)",tot);
    TH1F *h = (TH1F*)gDirectory->Get("h");
    mmiss[i][0] = (TH1F*)h->Clone("cand"+i);
    tot = dss+dssMva+type2; tot += cuts;
    b.Draw("mm2pi0>>h(20,-2.5,5.)",tot);
    TH1F *h = (TH1F*)gDirectory->Get("h");
    mmiss[i][1] = (TH1F*)h->Clone("candMva"+i);
    double maxi = mmiss[i][0]->GetMaximum();
    if(maxi<mmiss[i][1]->GetMaximum()) maxi = mmiss[i][1]->GetMaximum();

    mmiss[i][0]->SetMaximum(1.15*maxi);
    mmiss[i][0]->GetXaxis()->SetTitleOffset(.68);
    mmiss[i][0]->GetXaxis()->SetTitleSize(0.06);
    mmiss[i][0]->GetYaxis()->SetLabelSize(0.062);
    mmiss[i][0]->GetXaxis()->SetLabelSize(0.06);
    mmiss[i][0]->GetYaxis()->SetNdivisions(6+100*2);
    mmiss[i][0]->GetXaxis()->SetNdivisions(7+100*2);
    mmiss[i][0]->SetLineWidth(2);
    mmiss[i][0]->SetTitle("");
    mmiss[i][0]->SetXTitle("m^{2}_{miss} [GeV^{2}]");
    mmiss[i][0]->Draw();
    mmiss[i][1]->SetLineWidth(2);
    mmiss[i][1]->SetLineColor(2);
    mmiss[i][1]->Draw("same");
    TString tag = Dnames[i]; tag+= " channel: "; tag+= typeTitle; 
    label->DrawLatex(.11,0.93,tag);
    TLegend *leg;
    if(Tag=="Comb")
      leg = new TLegend(0.1,.73,0.33,0.9);
    else
      leg = new TLegend(0.67,.73,0.9,0.9);
    leg->SetTextSize(0.058);
    leg->SetFillColor(0);
    TString legtag = typeTitle; legtag+="^{0} ("; legtag+=(int)mmiss[i][0]->GetEntries(); legtag+=")";
    leg->AddEntry(mmiss[i][0],legtag);
    legtag = typeTitle; legtag+="^{+} ("; legtag+=(int)mmiss[i][1]->GetEntries(); legtag+=")";
    leg->AddEntry(mmiss[i][1],legtag);
    leg->Draw();
  }
  c.SaveAs("babar_code/eps/"+Tag+"Mmiss_"+run+".eps");
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
