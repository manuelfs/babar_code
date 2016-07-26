//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: WeightVariable.cc,v 1.4 2012/08/23 02:22:18 manuelf Exp $
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
  if (argc < 2 || argc > 6 ) {
    cout << "USAGE: WeightVariable Variable [BinName] [subBkg_AddD=yes_NO] [extraCuts] [wFileBB]" << endl;
    return 0;
  }
  // Setting the input parameters
  TString Variable = argv[1]; 
  TString BinName = "babar_code/WeightsVari/BinName.txt"; 
  if (argc>2) BinName = argv[2];
  TString subBkg = "yes_NO";
  if (argc>3) subBkg = argv[3];
  TString extraCuts = "";
  if (argc>4) extraCuts = argv[4];
  extraCuts.ReplaceAll("XX","&&");
  TString weightName = "babar_code/Reweight/wCombCont.txt";
  if (argc>5) weightName = argv[5];

  fstream BinFile;
  BinFile.open(BinName,fstream::in);
  int nbins; double binning[20];
  BinFile >> nbins;
  for(int bin=0; bin<=nbins; bin++) BinFile >> binning[bin];
  TString textName = "../DonutUtils/Weights_"; textName+=Variable; textName+=".txt";
  TCut cuts = PMiss+M2P+lowQ2+MEScut+Mva+!MvaAll; 
  if(Variable.Contains("noMVA")){
    Variable.ReplaceAll("noMVA","");
    cuts = PMiss+M2P+lowQ2+MEScut+!Mva+!MvaAll; 
  }
  if(Variable.Contains("All")){
    Variable.ReplaceAll("All","");
    cuts = PMiss+M2P+lowQ2+MEScut+!MvaAll; 
  }
  cuts += extraCuts;
  double gEntries=0,uEntries=0,cEntries=0,dEntries=0;
  TString genName = "AWG82/ntuples/small/RAll_RunAll.root";
  TString tempFolder = "AWG82/ntuples/temp/"; 
  TString NameTree[4], TagTree[] = {"datatree", "BBtree", "udstree", "ccbartree"};
  for(int tre=0; tre<4; tre++) {
    NameTree[tre] = tempFolder;
    NameTree[tre] += TagTree[tre];
    NameTree[tre] += ".root";
  }
  WeightedTree3("AWG82/ntuples/small/Data_RunAll.root", dEntries, "1",0,cuts,"All",NameTree[0]);
  WeightedTree3(genName, gEntries, weightName,0,cuts,"All",NameTree[1]);
  WeightedTree3("AWG82/ntuples/small/uds_RunAll.root", uEntries, weightName,-1,cuts,"All",NameTree[2]);
  WeightedTree3("AWG82/ntuples/small/ccbar_RunAll.root", cEntries, weightName,-1,cuts,"All",NameTree[3]);
  TChain *trees[4];
  for(int tre=0; tre<4; tre++) {
    trees[tre] = new TChain("ntp1");
    trees[tre]->Add(NameTree[tre]); 
  }

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
      histo[chan][i] = new TH1F(hname, "", nbins, binning);
      TString vari = Variable; vari += ">>"; vari += hname;
      TString totCut = "((candType=="; totCut += chan+1; 
      if(subBkg.Contains("YES")) {totCut += "||candType=="; totCut += chan+3;}
      if(subBkg.Contains("yes") && i==1) totCut += ")&&MCType>0&&MCType<13";
      //if(subBkg.Contains("yes") && i==1) totCut += ")&&MCType>0&&MCType<15"; // Including the D** in the weights
      else totCut += ")";
      totCut += ")*weight";
      can.cd(i+4*chan+1); trees[i]->Draw(vari,totCut);
    }
    if(subBkg.Contains("yes")) {
      TString hname = "CombSub"; hname += chan;
      histo[chan][4] = new TH1F(hname, "", nbins, binning);
      TString vari = Variable; vari += ">>"; vari += hname;
      TString totCut = "((candType=="; totCut += chan+1; 
      if(subBkg.Contains("YES")) {totCut += "||candType=="; totCut += chan+3;}
      totCut += ")&&(MCType==0||MCType>=13&&MCType<=14))*weight";
      //totCut += ")&&(MCType==0))*weight"; // Including the D** in the weights
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
  textFile<<nbins<<endl<<endl;
  for(int i=1; i<nbins+1; i++){
    textFile<<binning[i-1];
    for(int chan=0; chan<4; chan++)textFile<<"  \t"<<RoundNumber(histo[chan][0]->GetBinContent(i),5);
    textFile<<endl;
  }
  textFile<<binning[nbins]<<endl;
  textFile.close();
  cout<<"Saved "<<textName<<endl;

  return 1;
}


