#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "iostream.h"

void mergeNtuples(const char* name) {
   
  TChain ch[12];
  for(int i=0; i<12; i++)  ch[i] = new TChain("ntp1");
  TString Mdir = name; Mdir += "/Merged/";
  gSystem->mkdir(Mdir);

  void *dir = gSystem->OpenDirectory(gSystem->ExpandPathName(name));
  int nfiles = 0;
  TString subdir;
  if (dir) {
    int ndir = 0;
    while ((subdir = gSystem->GetDirEntry(dir)) && ndir < 1) {
      if (subdir=="") break;
      if (subdir.Contains(".") || subdir.Contains("Merged")) continue;
      //if (!subdir.Contains("Run6")) continue;
      //if (subdir.Contains("1235") || subdir.Contains("Run1") || subdir.Contains("Run2") || subdir.Contains("Run3")) continue;
      TString basename(subdir);
      Int_t slashpos = basename.Last('/');
      if (slashpos>=0) {
	basename.Remove(0,slashpos+1);      // Remove what comes before the /
      } else {
	basename = "File"; basename += ndir;
      }
      TString outfile = Mdir;
      outfile += subdir; outfile += ".root";
      subdir = name + subdir; subdir += "/*";  
      ch[ndir].Add(subdir);

      cout<<"Merging "<<subdir<<endl;
      ch[ndir].Merge(outfile);
      ndir++;
    }
  } else cout<<"manuelf: The directory does not exist"<<endl;
}   






