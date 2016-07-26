//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: WeightVariable.cc,v 1.1 2010/10/21 18:25:43 manuelf Exp $
//
// Description:
//      WeightVariable - Re-weights a generic variable in bins. By defaults
//                       it calculates weights after BDT cuts
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      10/08/12 manuelf -- Created
//------------------------------------------------------------------------

#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TTree.h"
#include "TSystem.h"
#include "DonutUtils/cuts.cc"
#include "DonutUtils/KeysUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

int main(int argc, char *argv[]){
  if (argc < 5 || argc > 8 ) {
    cout << "USAGE: WeightVariable Variable nbins minX maxX [subBkg_AddD=yes_YES] [extraCuts] [wFileBB]" << endl;
    return 0;
  }
  // Setting the input parameters
  TString Variable = argv[1]; 
  TString temp_s = argv[2]; int nbins = temp_s.Atoi(); 
  temp_s = argv[3]; float minX = temp_s.Atof(); 
  temp_s = argv[4]; float maxX = temp_s.Atof(); 
  TString subBkg = "yes_YES";
  if (argc>5) subBkg = argv[5];
  TString extraCuts = "";
  if (argc>6) extraCuts = argv[6];
  extraCuts.ReplaceAll("XX","&&");
  TString weightName = "babar_code/Reweight/wCombCont.txt";
  if (argc>7) weightName = argv[7];

  TString textName = "../DonutUtils/Weights_"; textName+=Variable; textName+=".txt";
  TCut cuts = PMiss+M2P+Mva+!MvaAll; 
  if(Variable.Contains("noMVA")){
    Variable.ReplaceAll("noMVA","");
    cuts = PMiss+M2P+!Mva+!MvaAll; 
  }
  cuts += extraCuts;
  double gEntries=0,uEntries=0,cEntries=0,dEntries=0;
  TString genName = "AWG82/ntuples/small/RAll_RunAll.root";
  TTree *trees[4];
  trees[0] = WeightedTree("AWG82/ntuples/small/Data_RunAll.root", dEntries, "1",0,cuts,"All","datatree");
  trees[1] = WeightedTree(genName, gEntries, weightName,0,cuts,"All","BBtree");
  trees[2] = WeightedTree("AWG82/ntuples/small/uds_RunAll.root", uEntries, weightName,-1,cuts,"All","udstree");
  trees[3] = WeightedTree("AWG82/ntuples/small/ccbar_RunAll.root", cEntries, weightName,-1,cuts,"All","ccbartree");
  TString sNames[4] = {"Data", "BB", "uds", "ccbar"};
  double totMCB = 0, totuds = 0, totccbar = 0, totdata = 0, totOffdata = 0;
  getNumberB(genName, "All", totMCB, totdata, totuds, totccbar, totOffdata);
  double wuds = totMCB/totuds*2.09/1.05;              // 0.862068
  double wccbar = totMCB/totccbar*1.3/1.05;           // 1.041766
  double wMC = totdata/totMCB;                        // 0.496529
  wccbar = wuds;

  gROOT->SetStyle("Plain");gStyle->SetOptStat(10);
  TCanvas can("can","Distributions",1000,700);
  can.Divide(4,4); 
  TH1F* histo[4][5];
  int Nchan = 4;
  if(subBkg.Contains("YES")) Nchan = 2;
  for(int chan = 0; chan < Nchan; chan++){
    for(int i = 0; i < 4; i++){
      TString hname = sNames[i]; hname += chan;
      histo[chan][i] = new TH1F(hname, "", nbins, minX, maxX);
      TString vari = Variable; vari += ">>"; vari += hname;
      TString totCut = "((candType=="; totCut += chan+1; 
      if(subBkg.Contains("YES")) {totCut += "||candType=="; totCut += chan+3;}
      if(subBkg.Contains("yes") && i==1) totCut += ")&&MCType>0";
      else totCut += ")";
      totCut += ")*weight";
      can.cd(i+4*chan+1); trees[i]->Draw(vari,totCut);
    }
    if(subBkg.Contains("yes")) {
      TString hname = "CombSub"; hname += chan;
      histo[chan][4] = new TH1F(hname, "", nbins, minX, maxX);
      TString vari = Variable; vari += ">>"; vari += hname;
      TString totCut = "((candType=="; totCut += chan+1; 
      if(subBkg.Contains("YES")) {totCut += "||candType=="; totCut += chan+3;}
      totCut += ")&&MCType==0)*weight";
      trees[1]->Draw(vari,totCut);
      histo[chan][0]->Add(histo[chan][4],-wMC);
      histo[chan][0]->Add(histo[chan][2],-wuds*wMC);
      histo[chan][0]->Add(histo[chan][3],-wccbar*wMC);
    } else {
      histo[chan][1]->Add(histo[chan][2],wuds);
      histo[chan][1]->Add(histo[chan][3],wccbar);
    }
    histo[chan][0]->Divide(histo[chan][0],histo[chan][1],1,wMC);
    if(subBkg.Contains("YES")) histo[chan+2][0] = (TH1F*)histo[chan][0]->Clone();
  }
  can.SaveAs("WeightVariable.eps");
  fstream textFile;
  textFile.open(textName,fstream::out);
  textFile<<nbins<<"\t"<<minX<<"\t"<<maxX<<endl<<endl;
  double width = (maxX-minX)/(double)nbins;
  double variBin = minX;
  for(int i=1; i<nbins+1; i++){
    textFile<<variBin;
    for(int chan=0; chan<4; chan++)textFile<<"  \t"<<RoundNumber(histo[chan][0]->GetBinContent(i),5);
    textFile<<endl;
    variBin += width;
  }
  textFile.close();
  cout<<"Saved "<<textName<<endl;

  return 1;
}


