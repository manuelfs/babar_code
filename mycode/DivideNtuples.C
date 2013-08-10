#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"
#include "TTree.h"
#include <iostream>
#include <string>

#define nFiles 20

void DivideNtuples(const char* name, int nRests = 7) {

  TString MergedDir = name; 
  if(!MergedDir.Contains("Merged")) {cout<<"It doesn't have Merged"<<endl; return;}
   
  TChain *ch[nFiles];
  for(int i=0; i<nFiles; i++)  ch[i] = new TChain("ntp1");
  TChain *chRests[nFiles][nFiles];
  TString MdirRests[nFiles];
  MergedDir.Remove(MergedDir.Last('M'),MergedDir.Sizeof());
  for(int i=0; i<nRests; i++){
    MdirRests[i] = MergedDir; MdirRests[i] += "/MergedRests_";MdirRests[i] += i+1;
    MdirRests[i] += "_";MdirRests[i] += nRests; MdirRests[i] += "/";
    gSystem->mkdir(MdirRests[i]);
  }
  void *dir = gSystem->OpenDirectory(gSystem->ExpandPathName(name));
  TString subdir;
  if (dir) {
    int ndir = 0;
    while ((subdir = gSystem->GetDirEntry(dir)) && ndir < 40) {
      if (subdir=="") break;
      if (subdir.Contains("..") || subdir == ".") continue;
      if (!subdir.Contains("1235")) continue;
      if (!subdir.Contains("Run6-")) continue;
      TString fname = name + subdir; 
      subdir.Remove(subdir.Last('.'),subdir.Sizeof());
      TString outfileRests[nFiles];
      for(int i=0; i<nRests; i++){
	outfileRests[i] = MdirRests[i];	outfileRests[i] += subdir; 
	outfileRests[i] += "Rests_"; outfileRests[i] += i+1; outfileRests[i] += "_"; 
	outfileRests[i] += nRests; outfileRests[i] += ".root";
      }
      ch[ndir]->Add(fname);
      Int_t lowerID;
      TBranch *b_lowerID;
      ch[ndir]->SetBranchAddress("lowerID", &lowerID, &b_lowerID);
      int entries = ch[ndir]->GetEntries(); 
      int minR = 0, maxR = nRests;
      for(int i=minR; i<maxR; i++) chRests[i][ndir] = (TChain *)ch[ndir]->CloneTree(0);
      cout<<"Dividing  "<<fname<<" with "<<entries<<" entries"<<endl;
      int iRests = 0;
      for (int i=0; i<entries; i++){
	//if(i%1000==0) cout<<"Entry "<<i<<" of "<<entries<<endl;
	ch[ndir]->GetEntry(i);
	iRests = abs(lowerID % nRests);
	if(iRests>=minR && iRests<maxR) chRests[iRests][ndir]->Fill();
      }
      for(int i=minR; i<maxR; i++){
	TFile *Restfile = new TFile(outfileRests[i],"RECREATE");
	Restfile->cd();
	chRests[i][ndir]->Write();
	delete Restfile;
      }
      ndir++;
   }
  } else cout<<"manuelf: The directory does not exist"<<endl;
}   






