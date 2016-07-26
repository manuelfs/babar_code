#include "TCut.h"
#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <iostream>
using std::cout;
using std::endl;
TString round(double n, int d);

void optieex(){
  TChain b("ntp1");
  b.Add("AWG82/ntuples/mva/Add_BDT_KMAll.root");

  TH1F *hEExtra[4];
  for(int i=0; i<4; i++){
    TString hname = "EExtra"; hname += i;
    hEExtra[i] = new TH1F(hname,"Optimization EExtra",60,0,0.6);
  }

  TCut sigCut[] = {"candLepTru==1&&(MCType==5||MCType==6)","candLepTru==1&&MCType==6",
		   "candLepTru==1&&(MCType==11||MCType==12)","candLepTru==1&&MCType==12"};
  TCut normCut[] = {"candLepTru==1&&(MCType>0&&MCType<5)","candLepTru==1&&(MCType==2||MCType==4)",
		    "candLepTru==1&&(MCType>6&&MCType<11)","candLepTru==1&&(MCType==8||MCType==10)"};
  TCut crossCut[] = {"MCType>6&&MCType<13","MCType>6&&MCType<13",
		     "MCType>0&&MCType<7","MCType>0&&MCType<7"};
  TString Dnames[] = {"D^{0}","D*^{0}","D^{+}","D*^{+}"};
  TLatex *label = new TLatex();
  label->SetNDC(kTRUE);
  label->SetTextSize(0.06);
  Long64_t signal, norm, cross, dss, comb;

  TCanvas c("cEExtra","Optimization EExtra",800,600);
  c.Divide(2,2);
  for(int j=1; j<5; j++){
    c.cd(j);
    cout<<"Calculating "<<Dnames[j-1]<<endl;
    double hiSig = 0, hiEEx=0;
    for(int i=1; i<61; i++){
      if(i%16==0)cout<<"Bin "<<i<<" of 60"<<endl;
      double Ecut = ((double)i)/100.;
      TString Cut = "candType=="; Cut+=j; Cut+="&&candM2>1&&candEExtra<"; Cut += Ecut;
      Long64_t total = b.GetEntries(Cut);
      TCut sc = sigCut[j-1]; sc += Cut.Data(); 
      signal = b.GetEntries(sc);
      TCut nc = normCut[j-1]; nc += Cut.Data();
      norm = b.GetEntries(nc);
      TCut cc = crossCut[j-1]; cc += Cut.Data(); 
      cross = b.GetEntries(cc);
      TCut dssCut = "MCType>12"; dssCut += Cut.Data();
      dss = b.GetEntries(dssCut);
      comb = total-signal-norm-cross-dss;

      double den = signal+cross+dss+comb;
      double signi = 0;
      if(den) signi = signal*signal/den;
      hEExtra[j-1]->SetBinContent(i,signi);
      if(signi>hiSig){
	hiSig=signi;
	hiEEx = Ecut;
      }
    }
    gStyle->SetOptStat(0);
    TString hTitle = "Highest S^{2}/(S+B) is "; hTitle += round(hiSig,0);
    hTitle += ", for E_{Extra} > "; hTitle += round(hiEEx,2); 
    hEExtra[j-1]->SetTitle(hTitle);
    hEExtra[j-1]->SetXTitle("E_{Extra} cut [GeV]");
    hEExtra[j-1]->SetMinimum(0);
    hEExtra[j-1]->SetLineWidth(2);
    hEExtra[j-1]->Draw();
    label->DrawLatex(.82,0.93,Dnames[j-1]);
  }
  c.SaveAs("babar_code/eps/EExtraOpti.eps");

}

TString round(double n, int d){
  double b = ((int)(n*pow(10,d)+.5));
  b/=pow(10,d);
  TString result; result+= b;
  result.ReplaceAll(" ","");
  return result;
}
