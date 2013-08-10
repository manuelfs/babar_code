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
#include "TH2F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <iostream>
using std::cout;
using std::endl;
TString round(double n, int e, double d=1);

void optiMva2(){
  TChain b("ntp1");
  b.Add("AWG82/ntuples/small/Add_R24Rest_RunAll.root");
  TChain uds("ntp1");
  TChain ccbar("ntp1");
  uds.Add("AWG82/ntuples/small/udsRests_2_3_RunAll.root");
  uds.Add("AWG82/ntuples/small/udsRests_3_3_RunAll.root");
  ccbar.Add("AWG82/ntuples/small/ccbarRests_2_3_RunAll.root");
  ccbar.Add("AWG82/ntuples/small/ccbarRests_3_3_RunAll.root");
  double NBpRun[2][6] = {{34878000, 101690000, 56035000, 166784000, 215168000, 130336000},
			 {44172000, 130056000, 67039000, 216196000, 291011000, 160817000}};
  double NB0Run[2][6] = {{34941000, 104188000, 57888000, 169801000, 215953000, 135224000},
			 {43741000, 128709000, 67786000, 213996000, 245874000, 124453000}};
  double nccbar[] = {55254000, 164146000, 88321000, 267308000, 343667000, 208664000};
  double nuds[] = {160514000, 451636000, 275869000, 421599000, 553604000, 327032000};
  double nMCB = 0, totuds = 0, totccbar = 0, fracRest = 1/3.*3/2;
  for(int i=1; i<7; i++){
    totuds += nuds[i-1];
    totccbar += nccbar[i-1];
    nMCB += NBpRun[0][i-1]+NB0Run[0][i-1]; 
  }
  double wuds = nMCB/totuds*2.09/1.05*fracRest; 
  double wccbar = nMCB/totccbar*1.3/1.05*fracRest; 

  int xbin = 40, ybin = 40;
  double xlow = -0.2, ylow = -0.1, xhi = 0.6, yhi = 0.6;
  TH2F *hMVA[4];
  for(int i=0; i<4; i++){
    TString hname = "MVA"; hname += i;
    hMVA[i] = new TH2F(hname,"Optimization MVA",xbin,xlow,xhi,ybin,ylow,yhi);
  }

  TCut sigCut[] = {"MCType==5","MCType==6", "MCType==11","MCType==12"};
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
	TString Cut = "candType=="; Cut+=j; Cut+="&&candPMiss>.2&&candQ2>4&&candM2>";
	if(j%2==1)Cut+="2";
	else Cut +="1";
	Cut+="&&candMvaDl>"; Cut += Dlcut;
	Cut += "&&candMvaComb>"; Cut += Combcut;
	double total = b.GetEntries(Cut);
	TCut sc = sigCut[j-1]; sc += Cut.Data(); 
	signal = b.GetEntries(sc);
// 	if(j<3) {
// 	  total /= 1.04974787;
// 	  signal /= 1.04974787;
// 	} else {
// 	  total /= 1.0710204897;
// 	  signal /= 1.0710204897;
// 	}

	total += uds.GetEntries(Cut)*wuds+ccbar.GetEntries(Cut)*wccbar;
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
  c.SaveAs("babar_code/eps/MVA2Opti.eps");
}


TString round(double n, int e, double d){
  if(d==0) return " - ";
  double neg = 1; if(n*d<0) neg = -1;
  double b = (int)(neg*n/d*pow(10.,(double)e)+0.5);
  b /= pow(10.,(double)e)*neg;
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(!result.Contains(".") && e != 0) result += ".";
  
  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<e-afterdot.Length(); i++)
    result += "0";
  return result;
}
