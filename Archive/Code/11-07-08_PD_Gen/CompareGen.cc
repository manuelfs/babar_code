//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: CompareGen.cc,v 1.2 2011/02/22 01:01:50 manuelf Exp $
//
// Description:
//      CompareGen - Compares the distributions of cocktail and generics
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      10/01/20 manuelf -- Created
//------------------------------------------------------------------------

#include "TString.h"
#include "TFile.h"
#include "TCut.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TSystem.h"
#include "DonutUtils/cuts.cc"
#include "DonutUtils/weightManager.hh"
#include "DonutUtils/KeysUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

int main(int argc, char *argv[]){
  if (argc < 2 || argc > 4 ) {
    cout << "USAGE: CompareGen Type [mcGen=Rest2] [weightFile]" << endl;
    return 0;
  }

  TString Tag = argv[1];
  TString mcGen = "Rest2";
  if (argc>2) mcGen = argv[2];
  TString weightName = "babar_code/Reweight/wFFBF.txt";
  if(Tag.Contains("pi0") || Tag.Contains("Dss")) 
    weightName = "babar_code/Reweight/wFFBFSig.txt";
  if (argc>3) weightName = argv[3];
  TCut Btag = "(candBMode>11000&&candBMode<11400||candBMode>12000&&candBMode<12500)";
  //MvaAll += Btag;

  if(Tag.Contains("pi0")) MvaAll = dssMvaAll;
  double genEntries=0, sigEntries=0, oldEntries=0;
  TString genfile = "AWG82/ntuples/small/Gen_Add"; genfile+=mcGen;genfile+="_RunAll.root";
  TTree *gen = WeightedTree(genfile, genEntries, "1",0,MvaAll);
  TTree *sig = WeightedTree("AWG82/ntuples/small/*Full*All.root", sigEntries, weightName,1,MvaAll);
  sig->SetLineColor(2);
  TChain oldChain("ntp1");
  TString oldfile = "AWG82/ntuples/small/BDT_KM_RunAll.root";
  TTree *old = WeightedTree(oldfile, oldEntries, "1",0,MvaAll);
  old->SetLineColor(4);


  TCut modeCut[2];
  TString decay[2];
  TString Titles[2];
  TString limits[2];
  if(Tag=="D"){
    modeCut[0] = "(MCType==1||MCType==3)&&candType==1";
    modeCut[1] = "(MCType==7||MCType==9)&&candType==3";
    Titles[0] = "D^{0} channel: ";
    Titles[1] = "D^{+} channel: ";
    decay[0] = "D^{0} (e/#mu) #nu";
    decay[1] = "D^{+} (e/#mu) #nu";
    limits[0] = "candM2>>h(50,-0.6,1.)";
    limits[1] = "candPstarLep>>h(50,0,2.4)";
  }else if(Tag=="Ds"){
    modeCut[0] = "(MCType==2||MCType==4)&&candType==2";
    modeCut[1] = "(MCType==8||MCType==10)&&candType==4";
    Titles[0] = "D*^{0} channel: ";
    Titles[1] = "D*^{+} channel: ";
    decay[0] = "D*^{0} (e/#mu) #nu";
    decay[1] = "D*^{+} (e/#mu) #nu";
    limits[0] = "candM2>>h(50,-0.6,1)";
    limits[1] = "candPstarLep>>h(50,0,2.4)";
  }else if(Tag=="Dfu"){
    modeCut[0] = "(MCType==1||MCType==3)&&candType==2";
    modeCut[1] = "(MCType==7||MCType==9)&&candType==4";
    Titles[0] = "D*^{0} channel: ";
    Titles[1] = "D*^{+} channel: ";
    decay[0] = "D^{0} (e/#mu) #nu";
    decay[1] = "D^{+} (e/#mu) #nu";
    limits[0] = "candM2>>h(50,-3.,3.)";
    limits[1] = "candPstarLep>>h(50,0,2.4)";
  }else if(Tag=="Dsfd"){
    modeCut[0] = "(MCType==2||MCType==4)&&candType==1";
    modeCut[1] = "(MCType==8||MCType==10)&&candType==3";
    Titles[0] = "D^{0} channel: ";
    Titles[1] = "D^{+} channel: ";
    decay[0] = "D*^{0} (e/#mu) #nu";
    decay[1] = "D*^{+} (e/#mu) #nu";
    limits[0] = "candM2>>h(50,-0.5,2.2)";
    limits[1] = "candPstarLep>>h(50,0,2.4)";
  }else if(Tag=="Dtau"){
    modeCut[0] = "(MCType==5)&&candType==1";
    modeCut[1] = "(MCType==11)&&candType==3";
    Titles[0] = "D^{0} channel: ";
    Titles[1] = "D^{+} channel: ";
    decay[0] = "D^{0} #tau #nu";
    decay[1] = "D^{+} #tau #nu";
    limits[0] = "candM2>>h(50,-0.5,10)";
    limits[1] = "candPstarLep>>h(50,0,2.4)";
  }else if(Tag=="Dstau"){
    modeCut[0] = "(MCType==6)&&candType==2";
    modeCut[1] = "(MCType==12)&&candType==4";
    Titles[0] = "D*^{0} channel: ";
    Titles[1] = "D*^{+} channel: ";
    decay[0] = "D*^{0} #tau #nu";
    decay[1] = "D*^{+} #tau #nu";
    limits[0] = "candM2>>h(50,-0.5,10)";
    limits[1] = "candPstarLep>>h(50,0,2.4)";
  }else if(Tag=="Dstaufd"){
    modeCut[0] = "(MCType==6)&&candType==1";
    modeCut[1] = "(MCType==12)&&candType==3";
    Titles[0] = "D^{0} channel: ";
    Titles[1] = "D^{+} channel: ";
    decay[0] = "D*^{0} #tau #nu";
    decay[1] = "D*^{+} #tau #nu";
    limits[0] = "candM2>>h(50,-0.5,10)";
    limits[1] = "candPstarLep>>h(50,0,2.4)";
  }else if(Tag=="Dssfdd"){
    modeCut[0] = "(MCType>=13&&MCType<=14&&candType==1)";
    modeCut[1] = "(MCType>=13&&MCType<=14)&&candType==3";
    Titles[0] = "D^{0} channel: ";
    Titles[1] = "D^{+} channel: ";
    decay[0] = "D** (e/#mu/#tau) #nu + NR";
    decay[1] = "D** (e/#mu/#tau) #nu + NR";
    limits[0] = "candM2>>h(50,-1,8)";
    limits[1] = "candPstarLep>>h(50,0,2.4)";
  }else if(Tag=="Dssfd"){
    modeCut[0] = "(MCType>=13&&MCType<=14&&candType==2)";
    modeCut[1] = "(MCType>=13&&MCType<=14)&&candType==4";
    Titles[0] = "D*^{0} channel: ";
    Titles[1] = "D*^{+} channel: ";
    decay[0] = "D** (e/#mu/#tau) #nu + NR";
    decay[1] = "D** (e/#mu/#tau) #nu + NR";
    limits[0] = "candM2>>h(50,-1,7)";
    limits[1] = "candPstarLep>>h(50,0,2.4)";
  }else if(Tag=="pi0DssD"){
    modeCut[0] = "(MCType>=13&&MCType<=14)&&candType==1";
    modeCut[1] = "(MCType>=13&&MCType<=14)&&candType==3";
    Titles[0] = "D^{0}#pi^{0} channel: ";
    Titles[1] = "D^{+}#pi^{0} channel: ";
    decay[0] = "D** (e/#mu/#tau) #nu + NR";
    decay[1] = "D** (e/#mu/#tau) #nu + NR";
    limits[0] = "candM2>>h(40,-0.5,2)";
    limits[1] = "candPstarLep>>h(35,0,2.4)";
  }else if(Tag=="pi0DssDs"){
    modeCut[0] = "(MCType>=13&&MCType<=14)&&candType==2";
    modeCut[1] = "(MCType>=13&&MCType<=14)&&candType==4";
    Titles[0] = "D*^{0}#pi^{0} channel: ";
    Titles[1] = "D*^{+}#pi^{0} channel: ";
    decay[0] = "D** (e/#mu/#tau) #nu + NR";
    decay[1] = "D** (e/#mu/#tau) #nu + NR";
    limits[0] = "candM2>>h(50,-0.5,2)";
    limits[1] = "candPstarLep>>h(50,0,2.4)";
  }else if(Tag=="pi0Dfu"){
    modeCut[0] = "(MCType==1||MCType==3)&&candType==1";
    modeCut[1] = "(MCType==7||MCType==9)&&candType==3";
    Titles[0] = "D^{0}#pi^{0} channel: ";
    Titles[1] = "D^{+}#pi^{0} channel: ";
    decay[0] = "D^{0} (e/#mu) #nu";
    decay[1] = "D^{+} (e/#mu) #nu";
    limits[0] = "candM2>>h(50,-4.,2.)";
    limits[1] = "candPstarLep>>h(50,0,2.4)";
  }else if(Tag=="pi0Dsfu"){
    modeCut[0] = "(MCType==2||MCType==4)&&candType==2";
    modeCut[1] = "(MCType==8||MCType==10)&&candType==4";
    Titles[0] = "D*^{0}#pi^{0} channel: ";
    Titles[1] = "D*^{+}#pi^{0} channel: ";
    decay[0] = "D*^{0} (e/#mu) #nu";
    decay[1] = "D*^{+} (e/#mu) #nu";
    limits[0] = "candM2>>h(50,-4,2.)";
    limits[1] = "candPstarLep>>h(50,0,2.4)";
  }else {
    cout<<"Select one: D, Ds, Dfu, Dsfd, Dtau, Dstau, Dstaufd, Dssfdd, Dssfd, pi0DssD, pi0DssDs, pi0Dfu, pi0Dsfu"<<endl;
    return 0;
  }
  TH1F *mmiss[4][3];
  TLatex *label = new TLatex();
  label->SetNDC(kTRUE);
  label->SetTextSize(0.08);
  TCanvas c("cMvaDl","Plot for the different channels",1000,600);
  c.Divide(2,2);
  gStyle->SetOptStat(0);c.SetFillColor(10);c.SetFrameFillColor(10);
  int a=0;
  for(int i=0; i<4; i++){
    a = i%2; 
    c.cd(i+1);
    gPad->SetFillColor(10);gPad->SetFrameFillColor(10);
    TCut tot = modeCut[i>1]; tot *= "weight";
    TString vari = limits[a]; 
    if(Tag.Contains("pi0") && a==0){
      vari.Replace(vari.First("candM2"),6,"mm2pi0"); 
    }
    cout<<"Doing "<<vari<<endl;
    TCut tot2 = tot; //tot2+="MCTaumode==-1"; 
    gen->Draw(vari,tot);
    TH1F *h = (TH1F*)gDirectory->Get("h");
    mmiss[i][0] = (TH1F*)h->Clone("cand"+i);
    sig->Draw(vari,tot);
    h = (TH1F*)gDirectory->Get("h");
    mmiss[i][1] = (TH1F*)h->Clone("candMva"+i);
    old->Draw(vari,tot);
    h = (TH1F*)gDirectory->Get("h");
    mmiss[i][2] = (TH1F*)h->Clone("candMva"+i);
    double ngen = mmiss[i][0]->Integral();
    double nsig = mmiss[i][1]->Integral();
    double nold = mmiss[i][2]->Integral();
    if(nsig) {
      if(ngen) mmiss[i][1]->Scale(ngen/nsig);
      else if(nold) mmiss[i][2]->Scale(nsig/nold);
    }      
    if(nold && ngen) mmiss[i][2]->Scale(ngen/nold);

    formatHisto(mmiss[i][0]);
    mmiss[i][0]->SetXTitle("m^{2}_{miss} [GeV^{2}]");
    if(a) mmiss[i][0]->SetXTitle("p*_l [GeV]");
    double maxi = mmiss[i][0]->GetMaximum();
    if(maxi<mmiss[i][1]->GetMaximum()) maxi = mmiss[i][1]->GetMaximum();
    if(maxi<mmiss[i][2]->GetMaximum()) maxi = mmiss[i][2]->GetMaximum();
    mmiss[i][0]->SetMaximum(1.15*maxi);
    mmiss[i][0]->Draw();
    mmiss[i][1]->SetLineWidth(2);
    mmiss[i][1]->SetLineColor(2);
    mmiss[i][2]->Draw("same");
    mmiss[i][0]->Draw("same");
    mmiss[i][1]->Draw("same");
    TString tag = Titles[i>1]; tag+= decay[i>1]; 
    label->DrawLatex(.11,0.92,tag);
    TLegend *leg;
    if(Tag=="pi0Dsfu"||Tag=="pi0Dfu")
      leg = new TLegend(0.1,.66,0.43,0.9);
    else
      leg = new TLegend(0.57,.66,0.9,0.9);
    leg->SetTextSize(0.058);
    leg->SetFillColor(0);
    TString legtag = "BrecoAdd ("; legtag+=RoundNumber(ngen,0); legtag+=")";
    leg->AddEntry(mmiss[i][0],legtag);
    legtag = "Breco ("; legtag+=RoundNumber(nold,0); legtag+=")";
    leg->AddEntry(mmiss[i][2],legtag);
    legtag = "Cocktail ("; legtag+=RoundNumber(nsig,0); legtag+=")";
    leg->AddEntry(mmiss[i][1],legtag);
    if(a==0)leg->Draw();
  }
  c.SaveAs("babar_code/eps/"+Tag+"compareGen"+mcGen+".eps");

  return 1; 

}

