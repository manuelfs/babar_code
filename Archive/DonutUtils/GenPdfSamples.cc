//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: GenPdfSamples.cc,v 1.9 2012/03/04 00:33:02 manuelf Exp $
//
// Description:
//      GenPdfSamples - Creates the samples to be fitted 
//
// Author List:
//      Michael Mazur (original)                  UC Santa Barbara
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      10/10/21 manuelf -- New numbering and some clean up.
//               Merging with GenCombSamples.cc
//      10/02/11 manuelf -- Deleted the duplicaded "dss" variables and 
//               adapted to the new numbering
//      09/05/05 manuelf -- Adapted to new code
//------------------------------------------------------------------------

#include "TROOT.h"
#include "TCut.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "DonutUtils/KeysUtils.cc"
#include <fstream>
#include <iostream>
using std::cout;
using std::endl;

Int_t MCType,candType,candIsMu;
Float_t candM2,candPstarD,candPstarLep,weight;
double totMCB = 0, totuds = 0, totccbar = 0, totdata = 0, totOffdata = 0;

void fillTree(TTree *inputTree, TTree *outputTree, Int_t thecut, Int_t dss, TString sepTru_Print);

int main(int argc, char *argv[])
{
  if (argc < 2 || argc > 3) {
    cout << "USAGE: GenPdfSamples ntupleTag=RAll [Print_SCALE = yes_YES]" << endl;
    return 1;
  }

  TString ntupleTag = argv[1];
  TString Print_SCALE = "yes_YES";
  if(argc>2)  Print_SCALE = argv[2]; 

  TString ntupleName = "AWG82/ntuples/small/Fit"; ntupleName += ntupleTag; ntupleName += "_RunAll.root";
  TChain *chain = new TChain("ntp1");
  chain->Add(ntupleName);
  cout<<endl<<"Running over "<<ntupleName<<" with "<<chain->GetEntries()<<" entries"<<endl<<endl;
  getNumberB(ntupleName, "All", totMCB, totdata, totuds, totccbar, totOffdata);
  if(ntupleTag.Contains("Dpipi") || ntupleTag.Contains("GenDss")) 
    {totMCB = 1; totdata = 1; totuds = 1;}

  chain -> SetBranchAddress("candPstarLep",  &candPstarLep);            
  chain -> SetBranchAddress("candPstarD",    &candPstarD);            
  chain -> SetBranchAddress("candM2",        &candM2);                  
  chain -> SetBranchAddress("candType",      &candType);                
  chain -> SetBranchAddress("candIsMu",      &candIsMu);                
  chain -> SetBranchAddress("MCType",        &MCType);                  
  chain -> SetBranchAddress("weight",        &weight);                       

  TFile *f;
  TString pdfFolder = "fitSamples/"; pdfFolder += ntupleTag; 
  gSystem->mkdir(pdfFolder);
  int iBeg = 1, iEnd = 70;
  if(ntupleTag.Contains("GenDss")) {iBeg = 33; iEnd = 40;}
  for (int i = iBeg ; i <= iEnd ; i ++) {
    if(i==51) i=55; if(i==59) i=61;   //Skipping non-existent numbers, and the signal combinatoric
    int isDss = 0;
    if(i>=21 && i<=36 || i>=45 && i<=49 || i>=55 && i<=58 || i>=65 && i<=68) isDss = 1;
    TString pdfName = pdfFolder; pdfName += "/pdfSample";pdfName += i; pdfName += ".root";
    f = new TFile(pdfName,"RECREATE");
    TTree *t = new TTree("ntp1","cands");
    fillTree(chain,t,i,isDss, Print_SCALE);
    t->Write();
    t->ResetBranchAddresses();
    f->Write();
    f->Close();
  } 

  if(ntupleTag.Contains("GenDss")) return 0;

  // Adding the MVA events that were not used in the MVA training
  if(Print_SCALE.Contains("MVA")){
    TString MVAName = "AWG82/ntuples/small/FitR24MVA_RunAll.root";
    cout<<endl<<"Adding "<<MVAName<<endl;
    chain->Add(MVAName);
  }
  for (int i = 51 ; i <= 54 ; i ++) {
    TString pdfName = pdfFolder; pdfName += "/pdfSample";pdfName += i; pdfName += ".root";
    f = new TFile(pdfName,"RECREATE");
    TTree *t = new TTree("ntp1","cands");
    fillTree(chain,t,i,0, Print_SCALE);
    t->Write();
    t->ResetBranchAddresses();
    f->Write();
    f->Close();
  } 

  cout << "Done with all samples..." << endl;
  return 0;
}

void fillTree(TTree *inputTree, TTree *outputTree, Int_t thecut, Int_t isDss, TString Print_SCALE){
  bool doPrint_SCALE = true;
  if(Print_SCALE.Contains("no")) doPrint_SCALE = false;
  Float_t myCandM2;

  outputTree -> Branch("MCType",&MCType,"MCType/I");
  outputTree -> Branch("candPstarLep",&candPstarLep,"candPstarLep/F");
  outputTree -> Branch("candPstarD",&candPstarD,"candPstarD/F");
  outputTree -> Branch("candM2",&myCandM2,"candM2/F");
  outputTree -> Branch("candType",&candType,"candType/I");
  outputTree -> Branch("weight",&weight,"weight/F");

  int doEmu = 2;
  if(Print_SCALE.Contains("doE")) doEmu = 1;
  if(Print_SCALE.Contains("doMu")) doEmu = 0;
  double totDss = 7049000+7036000, totWeight=0, totEvents=0;
  double xData=totMCB/totdata;
  if(thecut>=33&&thecut<=40) xData *= (0.01/(totDss/totdata));
  if(thecut>=51&&thecut<=54 && Print_SCALE.Contains("MVA")) xData+=1;
  if(thecut>=61&&thecut<=70) xData*=totuds/totMCB;
  bool Sep12 = false;
  if(thecut>=17 && thecut<=24) Sep12 = true;
  Int_t n1 = 0, n2 = 0, nevt = (Int_t)inputTree -> GetEntries();
  for (int j = 0 ; j < nevt ; j ++) {
    inputTree -> GetEvent(j);
    if(candIsMu == doEmu) continue;
    if(thecut>=35&&thecut<=36 || thecut>=39&&thecut<=40) candType += 2;
    if(!isSample(thecut, MCType, candType)) continue;
    totWeight += weight; totEvents++;
  }
  for (int j = 0 ; j < nevt ; j ++) {
    inputTree -> GetEvent(j);
    if(candIsMu == doEmu) continue;

    if(thecut>=35&&thecut<=36 || thecut>=39&&thecut<=40) candType += 2;
    if(!isSample(thecut, MCType, candType)) continue;
    if(Print_SCALE.Contains("YES")) weight *= (totEvents/totWeight);
    myCandM2 = candM2;
    outputTree -> Fill();
    if(!Sep12) n1++;
    else {
      if(thecut>=17&&thecut<25){
	if(MCType==13) n1++;
	else n2++;
      }
      if(thecut>=24&&thecut<29){
	if(MCType%2==1) n1++;
	else n2++;
      }
    }
  }

  if(doPrint_SCALE) {
    cout<<thecut<<"\t       "<<n1;
    if(Sep12) cout<<"/"<<n2;
    cout<<"\t&\t"<<RoundNumber(n1,0,xData);
    if(Sep12) cout<<"/"<<RoundNumber(n2,0,xData);
    cout<<endl;
  }
}
