//---------------------------------------------------------------------------------
// Description:
//      Optimization of the Dln and Comb MVA cuts for the D** sample
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      09/05/29 manuelf -- Adaptation from optimva.C
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

void optiDssConti(){
  TString mfile = "AWG82/ntuples/small/Add_R24Rest_RunAll.root";
  TChain b("ntp1");
  b.Add(mfile);
  TChain uds("ntp1");
  TChain ccbar("ntp1");
  uds.Add("AWG82/ntuples/small/uds30_Run12356.root");
  ccbar.Add("AWG82/ntuples/small/ccbar30_Run12356.root");

  double fracRest = 1.;
  if(mfile.Contains("Rest")) fracRest = 2/3.;
  double nccbar[] = {58900000., 168844000., 83974000., 0., 366758000., 104778000.}; // Run 4 is 252830000
  double nuds[] = {47180000., 130858000., 66892000., 0., 317846000., 84414000}; // Run 4 is 213380000
  double NBpRun[2][6] = {{36968000, 103498000, 50556000, 167994000, 244322000,  68148000},
                         {34878000, 101690000, 56035000, 166784000, 215168000, 130336000}};
  double NB0Run[2][6] = {{37200000, 103124000, 49766000, 167466000, 244812000,  68016000},
                         {34941000, 104188000, 57888000, 169801000, 215953000, 135224000}};
  double nMCB = 0, totuds = 0, totccbar = 0;
  for(int i=1; i<7; i++){
    totuds += nuds[i-1];
    totccbar += nccbar[i-1];
    nMCB += NBpRun[1][i-1]+NB0Run[1][i-1]; 
  }
  double wuds = nMCB/totuds*2.09/1.05*fracRest; 
  double wccbar = nMCB/totccbar*1.3/1.05*fracRest; 

  int xbin = 40, ybin = 40;
  double xlow = -0.8, ylow = -0.8, xhi = 0.1, yhi = 0.1;
  TH2F *hMVA[4];
  for(int i=0; i<4; i++){
    TString hname = "MVA"; hname += i;
    hMVA[i] = new TH2F(hname,"Optimization MVA",xbin,xlow,xhi,ybin,ylow,yhi);
  }

  TCut sigCut = "MCType>12";
  TString Dnames[] = {"D^{0}","D*^{0}","D^{+}","D*^{+}"};
  TLatex *label = new TLatex();
  label->SetNDC(kTRUE);
  label->SetTextSize(0.06);
  double signal;
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
	TString Cut = "eextrapi0<.5&&ppi0>.4&&mpi0>.125&&mpi0<.145&&abs(candCosT)<.8&&candType=="; 
	Cut+=j; Cut+="&&candMvaDssDl>"; Cut += Dlcut;
	Cut += "&&candMvaDssComb>"; Cut += Combcut;
	double total = b.GetEntries(Cut);
	double Nuds = uds.GetEntries(Cut)*wuds;
	double Nccbar = ccbar.GetEntries(Cut)*wccbar;
	total += Nuds; total += Nccbar;
	TCut sc = sigCut; sc += Cut.Data(); 
	signal = b.GetEntries(sc);

	double signi = 0;
	if(total) signi = signal/sqrt(total);
	hMVA[j-1]->SetBinContent(m,i,signi);
	if(signi>hiSig){
	  hiSig=signi;
	  hiDl = Dlcut; hiComb = Combcut;
	}
      }
    }
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    TString hTitle = "Highest S/sqr(S+B) is "; hTitle += round(hiSig,2);
    hTitle += ", for Dlnu > "; hTitle += round(hiDl,2); hTitle += " and Comb > ";
    hTitle += round(hiComb,2);
    hMVA[j-1]->SetTitle(hTitle);
    hMVA[j-1]->SetXTitle("Combinatoric MVA cut");
    hMVA[j-1]->SetYTitle("Dlnu MVA cut");
    hMVA[j-1]->Draw("cont4z");
    label->DrawLatex(.9,0.93,Dnames[j-1]); 
  }
  c.SaveAs("babar_code/eps/contiDssOpti2.eps");
}


TString round(double n, int d){
  double b = ((int)(n*pow(10,d)+.5));
  b/=pow(10,d);
  TString result; result+= b;
  result.ReplaceAll(" ","");
  return result;
}
