#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void GetErrors(TString folder="FitAll/Errors/", int nError = -1){

  TString ErrorName[] = {"fD**", "Constraints statistical", "PDF statistical", 
			 "Bkg. corrections", "D(*)lnu FF", "D**lnu BF", "Bkg. stat", 
			 "Crossfeed corrections", "Iso crossfeed"};
  TString Tag[] = {"R(D0):  ", "R(D*0): ", "R(D+):  ", "R(D*+): ", "R(D):   ", "R(D*):  ", 
		   "Correlation R(D0)-R(D*0) is ", "Correlation R(D+)-R(D*+) is ", 
		   "Correlation R(D)-R(D*) is   "};
  
  TH2F* h2;
  TH1F* histo;
  TString hName, nameFile;
  int begFile = 3, endFile = 11;
  if(nError>=3 && nError<=11){begFile = nError; endFile = nError;}
  for(int file=begFile; file<=endFile; file++){
    nameFile = folder; nameFile += "/Pulls"; nameFile += file; nameFile += "DataNewx100_RunAll.root";
    if(FILE *fCheck = fopen(nameFile,"r")) {
      fclose(fCheck);
      cout<<endl<<endl<<file<<" - "<<ErrorName[file-3]<<endl;
      cout<<"============================================================"<<endl;
      TFile PullFile(nameFile); PullFile.cd();
      for(int pad=0; pad<4; pad++){
	hName = "Ratio"; hName += pad+1;
	histo = (TH1F *)(PullFile.Get(hName));
	double rms = histo->GetRMS();
	double mean = histo->GetMean();
	cout<<Tag[pad]<<RoundNumber(mean,4)<<" +- "<<RoundNumber(rms,4)<<", a "<<
	  RoundNumber(100*rms,2,mean)<<" % error with "<<histo->GetEntries()<<" entries"<<endl;
	histo->Delete();
      }
      for(int pad=0; pad<2; pad++){
	hName = "Correl"; hName += pad;
	h2 = (TH2F *)(PullFile.Get(hName));
	cout<<Tag[pad+6]<<RoundNumber(h2->GetCorrelationFactor(),2)<<endl;
	h2->Delete();
      }
      PullFile.Close();
    }
    nameFile = folder; nameFile += "/Pulls"; nameFile += file; nameFile += "DataNewx100_RunAllIso.root";
    if(FILE *fCheck = fopen(nameFile,"r")) {
      fclose(fCheck);
      if(file==11){
	cout<<endl<<endl<<file<<" - "<<ErrorName[file-3]<<endl;
	cout<<"============================================================"<<endl;
      } else cout<<endl;
      TFile PullFile(nameFile); PullFile.cd();
      for(int pad=0; pad<2; pad++){
	hName = "Ratio"; hName += pad+1;
	histo = (TH1F *)(PullFile.Get(hName));
	double rms = histo->GetRMS();
	double mean = histo->GetMean();
	cout<<Tag[pad+4]<<RoundNumber(mean,4)<<" +- "<<RoundNumber(rms,4)<<", a "<<
	  RoundNumber(100*rms,2,mean)<<" % error with "<<histo->GetEntries()<<" entries"<<endl;
	histo->Delete();
	if(pad%2==1){
	  hName = "Correl"; hName += pad/2;
	  h2 = (TH2F *)(PullFile.Get(hName));
	  cout<<Tag[pad/2+8]<<RoundNumber(h2->GetCorrelationFactor(),2)<<endl;
	  h2->Delete();
	}
      }
      PullFile.Close();
    }
  }
}

