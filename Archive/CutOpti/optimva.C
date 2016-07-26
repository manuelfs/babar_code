//---------------------------------------------------------------------------------
// Description:
//      Optimization of the Dln and Comb MVA cuts
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      09/05/15 manuelf -- Adaptation from optieex.C
//---------------------------------------------------------------------------------

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

void optimva(){
  TChain b("ntp1");
  b.Add("AWG82/ntuples/mva/AddMC56.root");

  int xbin = 30, ybin = 30;
  double xlow = -0.5, ylow = -0.3, xhi = 0.5, yhi = 0.4;
  TH2F *hMVA[4];
  for(int i=0; i<4; i++){
    TString hname = "MVA"; hname += i;
    hMVA[i] = new TH2F(hname,"Optimization MVA",xbin,xlow,xhi,ybin,ylow,yhi);
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
  TCanvas c("cMVA","Optimization MVA",800,600);
  c.Divide(2,2);
  for(int j=1; j<5; j++){
    c.cd(j);
    cout<<"Calculating "<<Dnames[j-1]<<endl;
    double hiSig = 0, hiDl=0, hiComb=0;
    for(int i=1; i<ybin+1; i++){
      for(int m=1; m<xbin+1; m++){
	if(i%5==0 && m==1)cout<<"Bin "<<i<<" of "<<ybin<<endl;
	double Dlcut = ((double)i)*(yhi-ylow)/((double)ybin)+ylow;
	double Combcut = ((double)m)*(xhi-xlow)/((double)xbin)+xlow;
	TString Cut = "candType=="; Cut+=j; Cut+="&&candM2>1&&candMvaDl>"; Cut += Dlcut;
	Cut += "&&candMvaComb>"; Cut += Combcut;
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
	hMVA[j-1]->SetBinContent(m,i,signi);
	if(signi>hiSig){
	  hiSig=signi;
	  hiDl = Dlcut; hiComb = Combcut;
	}
      }
    }
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    TString hTitle = "Highest S^{2}/(S+B) is "; hTitle += round(hiSig,0);
    hTitle += ", for Dlnu > "; hTitle += round(hiDl,2); hTitle += " and Comb > ";
    hTitle += round(hiComb,2);
    hMVA[j-1]->SetTitle(hTitle);
    hMVA[j-1]->SetXTitle("Combinatoric MVA cut");
    hMVA[j-1]->SetYTitle("Dlnu MVA cut");
    hMVA[j-1]->Draw("cont4z");
    label->DrawLatex(.9,0.93,Dnames[j-1]); 
  }
  c.SaveAs("babar_code/eps/MVAOpti.eps");
}


TString round(double n, int d){
  double b = ((int)(n*pow(10,d)+.5));
  b/=pow(10,d);
  TString result; result+= b;
  result.ReplaceAll(" ","");
  return result;
}
