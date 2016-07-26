//---------------------------------------------------------------------------------
// Description:
//      Plots the D** FF weighted distributions for the 4 D**
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      09/11/22 manuelf -- Upgrade
//      09/05/31 manuelf -- Creation
//---------------------------------------------------------------------------------

#include "mycode/cuts.hh"
TString round(double n, int e, double d=1.);

void plotSingleDssFFw(TString Tag="Mass"){

  TString mfile = "AWG82/ntuples/small/semFull_RunAll.root";
  TChain semil("ntp1");
  semil.Add(mfile);

  TString limits[2], xTitle, vari;
  if(Tag=="Mass"){
    limits[0] = "(70,1.9,5)"; 
    limits[1] = "(70,2.3,2.7)";
    xTitle = "Mass [GeV]";
    vari = "trueDmass";
  } else if(Tag=="Q2"){
    limits[0] = "(50,0,9)"; 
    limits[1] = "(50,0,9)";
    xTitle = "q^{2} [GeV^2]";
    vari = "trueQ2";
  } else if(Tag=="w"){
    limits[0] = "(50,1,1.6)"; 
    limits[1] = "(50,1,1.6)";
    xTitle = "w";
    vari = "(5.279*5.279+trueDmass*trueDmass-trueQ2)/(2*5.279*trueDmass)";
  } else if(Tag=="Plep"){
    limits[0] = "(50,0,3.5)"; 
    limits[1] = "(50,0,3.5)";
    xTitle = "p_{l,lab} [GeV]";
    vari = "truePLep";
  } else {
    cout<<"Tag has to be on of {Mass, Q2, w, Plep}"<<endl;
    return;
  }

  TH1F *histo[4][2];
  TString Dnames[] = {"D*_{0} l #nu","D'_{1} l #nu","D_{1} l #nu","D*_{2} l #nu"};
  TCut modeCut[] = {"(abs(MCD)==10421||abs(MCD)==10411)","(abs(MCD)==20423||abs(MCD)==20413)",
		    "(abs(MCD)==10423||abs(MCD)==10413)", "(abs(MCD)==425||abs(MCD)==415)"};
  TCut weird = "MCTaumode==-1&&wFF<100";
  TLatex *label = new TLatex();
  label->SetNDC(kTRUE);
  label->SetTextSize(0.08);
  TCanvas c("canvas","Comparison of D** FF weights",1000,600);
  c.Divide(2,2);
  gStyle->SetOptStat(0);
  int a=0;
  for(int i=0; i<4; i++){
    c.cd(i+1);
    TString totvari = vari; totvari += ">>h"; totvari += limits[i>1]; 
    TCut tot = modeCut[i]+weird; 
    semil.Draw(totvari,tot);
    TH1F *h = (TH1F*)gDirectory->Get("h");
    histo[i][0] = (TH1F*)h->Clone("noWeight"+i);
    tot *= "wFF";
    semil.Draw(totvari,tot);
    h = (TH1F*)gDirectory->Get("h");
    if(h) histo[i][1] = (TH1F*)h->Clone("Weight"+i);
    else{ 
      cout<<"No re-weighted histogram"<<endl;
      continue;
    }
    double ngen = histo[i][0]->Integral();
    double nsig = histo[i][1]->Integral();
    double maxi = histo[i][0]->GetMaximum();
    if(maxi<histo[i][1]->GetMaximum()) maxi = histo[i][1]->GetMaximum();

    histo[i][0]->SetMaximum(1.15*maxi);
    histo[i][0]->GetXaxis()->SetTitleOffset(.8);
    histo[i][0]->GetXaxis()->SetTitleSize(0.06);
    histo[i][0]->GetYaxis()->SetLabelSize(0.06);
    histo[i][0]->GetXaxis()->SetLabelSize(0.06);
    histo[i][0]->GetYaxis()->SetNdivisions(5+100*2);
    histo[i][0]->GetXaxis()->SetNdivisions(6+100*2);
    histo[i][0]->SetLineWidth(2);
    histo[i][0]->SetTitle("");
    histo[i][0]->SetXTitle(xTitle);
    histo[i][0]->Draw();
    histo[i][1]->SetLineWidth(2);
    histo[i][1]->SetLineColor(2);
    histo[i][1]->Draw("same");
    TString tag = Dnames[i]; 
    label->DrawLatex(.12,0.93,tag);
    TLegend *leg;
    if(0)
      leg = new TLegend(0.1,.73,0.43,0.9);
    else
      leg = new TLegend(0.59,.73,0.9,0.9);
    leg->SetTextSize(0.058);
    leg->SetFillColor(0);
    TString legtag = "ISGW2 ("; legtag+=round(ngen,0); legtag+=")";
    leg->AddEntry(histo[i][0],legtag);
    legtag = "LLSW ("; legtag+=round(nsig,0); legtag+=")";
    leg->AddEntry(histo[i][1],legtag);
    leg->Draw();
  }
  c.SaveAs("babar_code/eps/DssFFWeight_"+Tag+".eps");
}

TString round(double n, int e, double d){
  if(d==0) return " - ";
  double b = (int)(n/d*pow(10,e)+0.5);
  b /= pow(10,e);
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(!result.Contains(".") && e != 0) result += ".";
  
  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<e-afterdot.Length(); i++)
    result += "0";
  return result;
}

