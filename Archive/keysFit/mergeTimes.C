#include "TH2F.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TChain.h"
#include "TKey.h"
#include "TString.h"
#include <fstream>
#include <iostream>
using std::cout;
using std::endl;
void mergeFiles(TString f1Name, TString f2Name, int doCount=0);

void mergeTimes(TString fromFolder, int doCount=1, TString toFolder = "keys/root/Times/"){
  TString f1Name, f2Name;
  for(int rep = 1; rep<=68; rep++){
    if(rep>32&&rep<41 || rep>48&&rep<51 || rep>58&&rep<61) continue;
    f1Name = toFolder; f2Name = fromFolder; 
    if(doCount) f1Name = fromFolder; 
    f1Name += "hTimes_"; f2Name += "hTimes_"; 
    f1Name += rep; f2Name += rep;
    f1Name += "_Fit.root"; f2Name += "_Fit.root";
    if(FILE *fCheck1 = fopen(f1Name,"r")){
      fclose(fCheck1); 
      if(FILE *fCheck2 = fopen(f2Name,"r")) {
	fclose(fCheck2);
	mergeFiles(f1Name, f2Name, doCount);
	//cout<<f1Name<<" does exist"<<endl;
      } //else cout<<f2Name<<" does not exist"<<endl;
    } //else cout<<f1Name<<" does not exist"<<endl;
    if(rep==20 || rep==32) cout<<endl;
  }
}

void mergeFiles(TString f1Name, TString f2Name, int doCount){
  TString fOption = "UPDATE";
  if(doCount) fOption = "";
  TFile f1(f1Name, fOption);

  TIter iter(f1.GetListOfKeys());
  TKey *key;
  int nHisto = 0;
  while ( (key = (TKey*)iter())) nHisto++;
  f1.cd();
  int nf1 = nHisto;
  if(doCount)  cout<<"There are "<<nf1<<" histograms in "<<f1Name<<endl;
  else {
    TFile f2(f2Name);
    TIter iter2(f2.GetListOfKeys());
    while ( (key = (TKey*)iter2())) {
      if(nHisto>=1000) break;
      TH2F *h2 = (TH2F*)key->ReadObj();
      TString hname = "hTimes"; hname += nHisto;
      f1.cd();
      h2->Write(hname); h2->Delete();
      nHisto++;
      f2.cd();
    }
    cout<<"Written "<<nHisto-nf1<<" histograms from "<<f2Name<<endl;
    cout<<f1Name<<": "<<nf1<<" and "<< nHisto<<endl;
    f2.Close();
  }
  f1.Close();  
}
