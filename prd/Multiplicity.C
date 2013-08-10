#include "TString.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "DonutUtils/cuts.cc"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;


void Multiplicity(){
  int candBnCharged = 0, candBnNeutral = 0, candBMode = 0, mult = 0;
  int nMult[10000][3];
  for(int dec=0; dec<10000; dec++)
    for(int var=0; var<3; var++) nMult[dec][var] = -1;

  TH1F *histo[3][2];
  TH2F histo2D("Cha_Neu","",12,0,12,12,0,12);
  TString hName, Bname[] = {"B0","Bp"}, Varname[] = {"Mult_", "Char_", "Neu_"};
  for(int var=0; var<3; var++)
    for(int his=0; his<2; his++) {
      hName = Varname[var]; hName += Bname[his];
      histo[var][his] = new TH1F(hName, "", 12, 0, 12);
    }
  int isBp;

  TChain tree("ntp1");
  tree.Add("AWG82/ntuples/small/RAll_RunAll.root");
  tree.SetBranchAddress("candBnCharged",&candBnCharged);
  tree.SetBranchAddress("candBnNeutral",&candBnNeutral);
  tree.SetBranchAddress("candBMode",&candBMode);
  long Nmax = tree.GetEntries();
  //Nmax = 100000;
  for(int entry=0; entry<Nmax; entry++){
    tree.GetEvent(entry);
    mult = candBnCharged + candBnNeutral;
    //cout<<"Mode "<<candBMode<<" has "<<mult<<" multiplicity"<<endl;
    nMult[candBMode-10000][0] = mult;
    nMult[candBMode-10000][1] = candBnCharged;
    nMult[candBMode-10000][2] = candBnNeutral;
      if(candBMode <12000 || (candBMode>=14000 && candBMode<16000)) isBp = 1;
      else isBp = 0;
      histo[0][isBp]->Fill(mult);
      histo[1][isBp]->Fill(candBnCharged);
      histo[2][isBp]->Fill(candBnNeutral);
histo2D.Fill(candBnCharged,candBnNeutral);
  }


  // double mean[2] = {0, 0};
  // int nMean[2] = {0, 0}, isBp;
  // for(int dec=0; dec<10000; dec++) {
  //   for(int var=0; var<3; var++) {
  //     if(nMult[dec][var]<0) continue;
  //     if(dec <2000 || (dec>=4000 && dec<6000)) isBp = 1;
  //     else isBp = 0;
  //     histo[var][isBp]->Fill(nMult[dec][var]);
  //     histo2D.Fill(nMult[dec][1],nMult[dec][2]);
  //   }
  // }
  TCanvas can("can","Multiplicities",800,1200); can.Divide(2,3);
  
  for(int var=0; var<3; var++)
    for(int his=0; his<2; his++)  {
      can.cd(his+2*var+1);
      histo[var][his]->Draw();
    }
  histo2D.Draw("");
  can.SaveAs("multi.gif");
  for(int var=0; var<3; var++)
    for(int his=0; his<2; his++)  histo[var][his]->Delete();
}




