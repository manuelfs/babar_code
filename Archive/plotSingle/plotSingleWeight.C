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

void plotSingleWeight(TString Tag="Q2"){

  //TString mfile = "AWG82/ntuples/small/Add_Truth30Rest_Run123456.root";
  TString mfile = "AWG82/ntuples/small/s*";
  TChain gen("ntp1");
  gen.Add(mfile);

  //TString title = "True q^{2} [GeV^{2}]";
  TString title = "True Vs reco. q^{2} [GeV^{2}]";
  TString limits[] = {"(50,0,12.5)","(50,0,12.5)"};
  if(Tag=="PLep"){
    limits[0] = "(50,0,4)"; 
    limits[1] = "(50,0,3)";
    title = "True p_{l} [GeV]";
  }

  TH1F *mmiss[4][2];
  TString Dnames[] = {"D l #nu","D* l #nu","D #tau #nu","D* #tau #nu"};
  TCut modeCut[] = {"(MCType==1||MCType==3||MCType==7||MCType==9)","(MCType==2||MCType==4||MCType==8||MCType==10)",
		    "(MCType==5||MCType==11)","(MCType==6||MCType==12)"};
  TCut weird = "wFF<100";
  TLatex *label = new TLatex();
  label->SetNDC(kTRUE);
  label->SetTextSize(0.08);
  TCanvas c("cMvaDl","Plot for the different channels",1000,600);
  c.Divide(2,2);
  gStyle->SetOptStat(0);
  int a=0;
  for(int i=0; i<4; i++){
    if(i>1) a = 1; 
    TString vari = "trueQ2>>h"; 
    if(Tag=="PLep"){
      vari = "truePLep>>h";
    }
    vari += limits[a];
    c.cd(i+1);
    TCut tot = modeCut[i]+weird; 
    tot *= "wFF";
    gen.Draw(vari,tot);
    TH1F *h = (TH1F*)gDirectory->Get("h");
    mmiss[i][0] = (TH1F*)h->Clone("cand"+i);
    vari = "candQ2>>h"; vari += limits[a];
    gen.Draw(vari,tot);
    h = (TH1F*)gDirectory->Get("h");
    if(h) mmiss[i][1] = (TH1F*)h->Clone("candMva"+i);
    else{ 
      cout<<"No re-weighted histogram"<<endl;
      continue;
    }
    double ngen = mmiss[i][0]->Integral();
    double nsig = mmiss[i][1]->Integral();
    double maxi = mmiss[i][0]->GetMaximum();
    if(maxi<mmiss[i][1]->GetMaximum()) maxi = mmiss[i][1]->GetMaximum();

    mmiss[i][0]->SetMaximum(1.15*maxi);
    mmiss[i][0]->GetXaxis()->SetTitleOffset(.8);
    mmiss[i][0]->GetXaxis()->SetTitleSize(0.06);
    mmiss[i][0]->GetYaxis()->SetLabelSize(0.06);
    mmiss[i][0]->GetXaxis()->SetLabelSize(0.06);
    mmiss[i][0]->GetYaxis()->SetNdivisions(5+100*2);
    mmiss[i][0]->GetXaxis()->SetNdivisions(6+100*2);
    mmiss[i][0]->SetLineWidth(2);
    mmiss[i][0]->SetTitle("");
    mmiss[i][0]->SetXTitle("");
    mmiss[i][0]->Draw();
    mmiss[i][1]->SetLineWidth(2);
    mmiss[i][1]->SetLineColor(2);
    mmiss[i][1]->Draw("same");
    TString tag = Dnames[i]; tag += ": "; tag += title;
    label->DrawLatex(.12,0.93,tag);
    TLegend *leg;
    if(i>1&&Tag=="Q2")
      leg = new TLegend(0.1,.73,0.43,0.9);
    else
      leg = new TLegend(0.57,.73,0.9,0.9);
    leg->SetTextSize(0.058);
    leg->SetFillColor(0);
    //TString legtag = "ISGW2 ("; legtag+=round(ngen,0); legtag+=")";
    TString legtag = "True ("; legtag+=round(ngen,0); legtag+=")";
    leg->AddEntry(mmiss[i][0],legtag);
    //legtag = "CLN ("; legtag+=round(nsig,0); legtag+=")";
    legtag = "Reco. ("; legtag+=round(nsig,0); legtag+=")";
    leg->AddEntry(mmiss[i][1],legtag);
    leg->Draw();
  }
  c.SaveAs("babar_code/eps/"+Tag+"FFWeightRec.eps");
}


TString round(double n, int e, double d){
  if(d==0) return " - ";
  double b = (int)(n/d*pow(10,e)+0.5);
  b /= pow(10,e);
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(b<10 && e==2){
    if(result.Length() == 1)
      result += ".00";
    if(result.Length() == 3)
      result += "0";
  }
  if(e==1 && !result.Contains(".")) result += ".0";
  return result;
}
