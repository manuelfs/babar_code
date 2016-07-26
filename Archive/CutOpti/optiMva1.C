//---------------------------------------------------------------------------------
// Description:
//      Optimization of the Dln cuts
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
TString round(double n, int e, double d=1);

void optiMva1(){

  TString weightName = "babar_code/Reweight/wFFBF.txt";
  TString weightConti = "babar_code/Reweight/wPl.txt";

  TString folder = "small/", tag="Add_R24", runs = "All", cuts = "candPMiss>.2&&candQ2>4&&candM2>1.5";
  double gEntries=0,uEntries=0,cEntries=0;
  if(tag.Contains("old")) {tag.ReplaceAll("old",""); folder = "Oldsmall/";}
  TString genName = "AWG82/ntuples/"; genName += folder; genName += tag; genName += "_RunAll.root";
  TString udsName = "AWG82/ntuples/"; udsName += folder; udsName += "uds_RunAll.root";
  TString ccbarName = "AWG82/ntuples/"; ccbarName += folder; ccbarName += "ccbar_RunAll.root";
  TTree *gen = WeightedTree(genName, gEntries, weightName,0,cuts,runs);
  TTree *uds = WeightedTree(udsName, uEntries, weightConti,-1,cuts);
  TTree *ccbar = WeightedTree(ccbarName, cEntries, weightConti,-1,cuts);

  double NBpRun[2][6] = {{34878000, 101690000, 56035000, 166784000, 215168000, 130336000},
			 {44172000, 130056000, 67039000, 216196000, 291011000, 160817000}};
  double NB0Run[2][6] = {{34941000, 104188000, 57888000, 169801000, 215953000, 135224000},
			 {43741000, 128709000, 67786000, 213996000, 245874000, 124453000}};
  double nccbar[] = {55254000, 164146000, 88321000, 267308000, 343667000, 208664000};
  double nuds[] = {160514000, 451636000, 275869000, 421599000, 553604000, 327032000};
  double nMCB = 0, totuds = 0, totccbar = 0, fracRest = 6/3.*3/2;
  for(int i=1; i<7; i++){
    totuds += nuds[i-1];
    totccbar += nccbar[i-1];
    nMCB += NBpRun[0][i-1]+NB0Run[0][i-1]; 
  }
  double wuds = nMCB/totuds*2.09/1.05*fracRest; 
  double wccbar = nMCB/totccbar*1.3/1.05*fracRest; 

  int xbin = 80;
  double xlow = -0.1, xhi = 0.7;
  TH1F *hMvaDl[4], hIntegral;
  for(int i=0; i<4; i++){
    TString hname = "MvaDl"; hname += i;
    hMvaDl[i] = new TH1F(hname,"Optimization MvaDl",xbin,xlow,xhi);
  }
  TString hInt = "hIntegral";
  hIntegral = new TH1F(hInt,"",5,1.5,12);
  TCut sigCut[] = {"MCType==5","MCType==6","MCType==11","MCType==12"};
  TString Dnames[] = {"D^{0}","D*^{0}","D^{+}","D*^{+}"};
  TLatex *label = new TLatex();
  label->SetNDC(kTRUE);
  label->SetTextSize(0.06);
  double signal;

  TCanvas c("cMvaDl","Optimization MvaDl",800,600);
  c.Divide(2,2);
  for(int j=1; j<5; j++){
    c.cd(j);
    cout<<"Calculating "<<Dnames[j-1]<<endl;
    double hiSig = 0, hiDl=-99;
    for(int i=1; i<xbin+1; i++){
      if(i%16==0)cout<<"Bin "<<i<<" of "<<xbin<<endl;
      double Dlcut = ((double)i)*(xhi-xlow)/((double)xbin)+xlow;
      TString Cut = "candType=="; Cut+=j; Cut+="&&candMvaDl>"; Cut += Dlcut;
      gen->Project(hInt,"candM2",Cut);
      double total = hIntegral->Integral();
      uds->Project(hInt,"candM2",Cut);
      total += hIntegral->Integral()*wuds;
      ccbar->Project(hInt,"candM2",Cut);
      total += hIntegral->Integral()*wccbar;
      TCut sc = sigCut[j-1]; sc += Cut.Data(); 
      gen->Project(hInt,"candM2",sc);
      signal = hIntegral->Integral();
      double signi = 0;
      if(total) signi = sqrt(signal*signal/total);
      hMvaDl[j-1]->SetBinContent(i,signi);
      if(signi>hiSig){
	hiSig=signi;
	hiDl = Dlcut;
      }
    }
    gStyle->SetOptStat(0);
    TString hTitle = "Highest S/sqr(S+B) is "; hTitle += round(hiSig,2);
    hTitle += ", for Dlnu > "; hTitle += round(hiDl,2); 
    hMvaDl[j-1]->SetTitle(hTitle);
    hMvaDl[j-1]->SetXTitle("Dlnu MVA cut");
    hMvaDl[j-1]->SetMinimum(0);
    hMvaDl[j-1]->SetLineWidth(2);
    hMvaDl[j-1]->Draw();
    label->DrawLatex(.82,0.93,Dnames[j-1]);
  }
  c.SaveAs("babar_code/eps/MvaDlOptiFinal_new.eps");

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

