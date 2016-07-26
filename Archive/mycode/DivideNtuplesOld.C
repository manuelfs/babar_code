#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"
#include "TTree.h"
#include <iostream>
#include <string>

#define nFiles 15

void DivideNtuplesOld(const char* name, int nRests = 2) {

  TString MergedDir = name; 
  if(!MergedDir.Contains("Merged")) {cout<<"It doesn't have Merged"<<endl; return;}
  int fMVA = 3;
   
  TChain *ch[nFiles];
  for(int i=0; i<nFiles; i++)  ch[i] = new TChain("ntp1");
  TChain *chRest[nFiles];
  TChain *chMVA[nFiles];
  TChain *chRest2[nFiles];
//   TChain *chMVAs[10][20];
//   TChain *chRests[10][20];
//   TChain *chRests2[10][20];
  MergedDir.Remove(MergedDir.Last('M'),MergedDir.Sizeof());
  TString MdirRest = MergedDir; MdirRest += "/MergedRest/";
  gSystem->mkdir(MdirRest);
  TString MdirMVA = MergedDir; MdirMVA += "/MergedMVA/";
  gSystem->mkdir(MdirMVA);
  TString MdirRest2 = MergedDir; MdirRest2 += "/MergedRest2/";
  gSystem->mkdir(MdirRest2);
//   TString MdirRests[10];
//   TString MdirRests2[10];
//   TString MdirMVAs[10];
//   for(int i=0; i<nRests; i++){
//     MdirRests[i] = MergedDir; MdirRests[i] += "/MergedRests_";MdirRests[i] += i;
//     MdirRests[i] += "_";MdirRests[i] += nRests;MdirRests[i] += "/";
//     gSystem->mkdir(MdirRests[i]);
//     MdirRests2[i] = MergedDir; MdirRests2[i] += "/MergedRests2_";MdirRests2[i] += i;
//     MdirRests2[i] += "_";MdirRests2[i] += nRests;MdirRests2[i] += "/";
//     gSystem->mkdir(MdirRests2[i]);
//     MdirMVAs[i] = MergedDir; MdirMVAs[i] += "/MergedMVAs_";MdirMVAs[i] += i;
//     MdirMVAs[i] += "_";MdirMVAs[i] += nRests;MdirMVAs[i] += "/";
//     gSystem->mkdir(MdirMVAs[i]);
//   }

  void *dir = gSystem->OpenDirectory(gSystem->ExpandPathName(name));
  TString subdir;
  if (dir) {
    int ndir = 0;
    while ((subdir = gSystem->GetDirEntry(dir)) && ndir < nFiles) {
      if (subdir=="") break;
      if (subdir.Contains("..") || subdir == ".") continue;
      //if (subdir.Contains("Run2") || subdir.Contains("Run4") || subdir.Contains("Run5")) continue;
      if (!subdir.Contains("5-BSemiExclAdd-Run4")) continue;
      TString fname = name + subdir; 
      TString outfileRest = MdirRest;
      subdir.Remove(subdir.Last('.'),subdir.Sizeof());
      outfileRest += subdir; outfileRest += "Rest.root";
      TString outfileMVA = MdirMVA;
      outfileMVA += subdir; outfileMVA += "MVA.root";
      TString outfileRest2 = MdirRest2;
      outfileRest2 += subdir; outfileRest2 += "Rest2.root";
//       TString outfileRests[10];
//       TString outfileRests2[10];
//       TString outfileMVAs[10];
//       for(int i=0; i<nRests; i++){
// 	outfileRests[i] = MdirRests[i];	outfileRests[i] += subdir; 
// 	outfileRests[i] += "Rests_"; outfileRests[i] += i; outfileRests[i] += "_"; 
// 	outfileRests[i] += nRests; outfileRests[i] += ".root";
// 	outfileRests2[i] = MdirRests2[i];	outfileRests2[i] += subdir; 
// 	outfileRests2[i] += "Rests2_"; outfileRests2[i] += i; outfileRests2[i] += "_"; 
// 	outfileRests2[i] += nRests; outfileRests2[i] += ".root";
// 	outfileMVAs[i] = MdirMVAs[i];	outfileMVAs[i] += subdir; 
// 	outfileMVAs[i] += "MVAs_"; outfileMVAs[i] += i; outfileMVAs[i] += "_"; 
// 	outfileMVAs[i] += nRests; outfileMVAs[i] += ".root";
//       }
      ch[ndir]->Add(fname);
      Int_t lowerID;
      TBranch *b_lowerID;
      ch[ndir]->SetBranchAddress("lowerID", &lowerID, &b_lowerID);
      int entries = ch[ndir]->GetEntries(); 
      chRest[ndir] = (TChain *)ch[ndir]->CloneTree(0);
//       chMVA[ndir] = (TChain *)ch[ndir]->CloneTree(0);
      chRest2[ndir] = (TChain *)ch[ndir]->CloneTree(0);
      //entries = 100;
//       for(int i=0; i<nRests; i++) chRests[i][ndir] = (TChain *)ch[ndir]->CloneTree(0);
//       for(int i=0; i<nRests; i++) chRests2[i][ndir] = (TChain *)ch[ndir]->CloneTree(0);
//       for(int i=0; i<nRests; i++) chMVAs[i][ndir] = (TChain *)ch[ndir]->CloneTree(0);
      cout<<"Dividing  "<<fname<<" with "<<entries<<" entries"<<endl;
      int iRests = 0, iRests2 = 0, iMVAs = 0;
      for (int i=0; i<entries; i++){
	//if(i%1000==0) cout<<"Entry "<<i<<" of "<<entries<<endl;
	ch[ndir]->GetEntry(i);
	if((lowerID % fMVA) == 0){
// 	  chMVA[ndir]->Fill();
// 	  chMVAs[iMVAs][ndir]->Fill();
	  iMVAs++;
	  if(iMVAs==nRests) iMVAs=0;
	}else if(abs(lowerID % fMVA) == 1){
	  chRest2[ndir]->Fill();
// 	  chRests2[iRests2][ndir]->Fill();
	  iRests2++;
	  if(iRests2==nRests) iRests2=0;
	}else{
	  chRest[ndir]->Fill();
// 	  chRests[iRests][ndir]->Fill();
	  iRests++;
	  if(iRests==nRests) iRests=0;
	}
      }
      TFile *Restfile = new TFile(outfileRest,"RECREATE");
      Restfile->cd();
      chRest[ndir]->Write();
      delete Restfile;
//       TFile *MVAfile = new TFile(outfileMVA,"RECREATE");
//       MVAfile->cd();
//       chMVA[ndir]->Write();
//       delete MVAfile;
      TFile *Rest2file = new TFile(outfileRest2,"RECREATE");
      Rest2file->cd();
      chRest2[ndir]->Write();
      delete Rest2file;
//       for(int i=0; i<nRests; i++){
// 	TFile *Restfile = new TFile(outfileRests[i],"RECREATE");
// 	Restfile->cd();
// 	chRests[i][ndir]->Write();
// 	delete Restfile;
// 	Restfile = new TFile(outfileRests2[i],"RECREATE");
// 	Restfile->cd();
// 	chRests2[i][ndir]->Write();
// 	delete Restfile;
// 	Restfile = new TFile(outfileMVAs[i],"RECREATE");
// 	Restfile->cd();
// 	chMVAs[i][ndir]->Write();
// 	delete Restfile;
//       }
      ndir++;
   }
  } else cout<<"manuelf: The directory does not exist"<<endl;
}   






