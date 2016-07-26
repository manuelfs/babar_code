#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace TMath;
using namespace std;
using std::cout;
using std::endl;

void Weight_HiggsPsl(){

  int nHig = 21, nBins = 50;
  double pMin = 0, pMax = 2, nEntries[2][2];
  TString treeName, folder = "babar_code/Reweight/HiggsPsl/wHiggsPsl_";
  TChain *Higtree[21];
  TH1F *Histo[2][2];
  for(int his=0; his<2; his++){
    for(int isDs=0; isDs<2; isDs++){
      treeName = "Histo"; treeName += his; treeName += isDs;
      Histo[his][isDs] = new TH1F(treeName,"",nBins,pMin,pMax);
    }
  }

  for(int hig=0; hig<nHig; hig++){

    int tBmH = hig*5;
    TString HigTag = "";
    if(tBmH<10) HigTag += "0";
    if(tBmH<100) HigTag += "0";
    HigTag += tBmH;

    Higtree[hig] = new TChain("GqaMCAnalysis/ntp1");
    treeName = "AWG82/generator/4M";
    treeName += HigTag; treeName += "Dxtaunu.root";
    Higtree[hig]->Add(treeName);

    for(int isDs=0; isDs<2; isDs++){
      treeName = "Histo"; treeName += hig>0; treeName += isDs;
      TString Cuts = "MCType==5";
      if(isDs) Cuts = "MCType==6";
      nEntries[hig>0][isDs] = Higtree[hig]->Project(treeName, "truePstarLep",Cuts);

      TString textName = folder; textName += isDs;
      textName += "_"; textName += HigTag; textName += ".txt";
      fstream textFile;
      textFile.open(textName,fstream::out);
      for(int bin=1; bin<=nBins; bin++){
	double nSM = Histo[0][isDs]->GetBinContent(bin), ratio = 1;
	if(nSM) ratio = Histo[hig>0][isDs]->GetBinContent(bin)/nSM
	  *(nEntries[0][isDs]/nEntries[hig>0][isDs]);
	textFile<<RoundNumber(ratio,4)<<endl;
      }
      textFile.close();
    }
  }

  cout<<"Written "<<nHig<<" values of tanB/mH to "<<folder<<endl;
  for(int his=0; his<2; his++) 
    for(int isDs=0; isDs<2; isDs++) Histo[his][isDs]->Delete();
}

