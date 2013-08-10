#include "babar_code/FF/bBToDtaunu.cc"
#include "babar_code/FF/bBToDstaunu.cc"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include "TMath.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TFile.h"
#include "TBox.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void TreeWeight(long nEntries=1000, TString tag = "tau"){

  BToDtaunu Dtaunu[3][4]; BToDstaunu Dstaunu[3][4];

  //TTree tree("ntp1","cands");
  TChain mc("GqaMCAnalysis/ntp1");
  if(!tag.Contains("lnu")){
    mc.Add("AWG82/generator/4MEvtGen_BpDxtaunu.root");
    mc.Add("AWG82/generator/OneMEvtGen_BpDxtaunu.root");
    //mc.Add("AWG82/generator/4MEvtGen_B0Dxtaunu.root");
    //mc.Add("AWG82/generator/OneMEvtGen_B0Dxtaunu.root");
  }else mc.Add("AWG82/generator/EvtGen_lnuDxtaunu.root");
  //mc.Add("AWG82/generator/5M1EvtGen_DsDxtaunu.root");

  int nHiggs = 2, nBad = 0;
  float wHig[4], tBmH[] = {0, 1, 0.50, 0.25}; 
  double PI   = 3.14159265, sumW[4][2][3], ml[] = {BToDtaunu::mE,BToDtaunu::mMu,BToDtaunu::mTau};
  long nW[4][3];
  TString Types[] = {"Dtaunu", "D*->Dpi", "D*->Dgamma"};
  for(int iHig=0; iHig<nHiggs; iHig++){
    for(int type=0; type<3; type++){
      if(type==2) {
	Dtaunu[type][iHig].SetTanBetamH(tBmH[iHig]);
	Dstaunu[type][iHig].SetTanBetamH(tBmH[iHig]);
      }
      Dtaunu[type][iHig].SetNorm(ml[type]);
      Dstaunu[type][iHig].SetNorm(ml[type]);
      nW[iHig][type] = 0;
      for(int rms=0; rms<2; rms++)
	sumW[iHig][rms][type] = 0;
    }
  }
  //float wCLN=-99., wH025=-99., wH050=-99., wH100=-99., tBmH[] = {0, 0.25, 0.50, 1.}; 
  float trueCTL=-99., trueCTV=-99., trueChi=-99., trueQ2=-99.; // manuelf
  int MCType=-1, isDgamma=-1;
//   tree.Branch("trueCTL",&trueCTL,"trueCTL/F");
//   tree.Branch("trueCTV",&trueCTV,"trueCTV/F");
//   tree.Branch("trueChi",&trueChi,"trueChi/F");
//   tree.Branch("trueQ2" ,&trueQ2 ,"trueQ2/F" );
//   tree.Branch("MCType" ,&MCType ,"MCType/I" );
//   tree.Branch("isDgamma" ,&isDgamma ,"isDgamma/I" );
//   tree.Branch("wCLN" , &wHig[0] ,"wCLN/F" );
//   tree.Branch("wH100" ,&wHig[1] ,"wH100/F" );
//   tree.Branch("wH025" ,&wHig[3] ,"wH025/F" );
//   tree.Branch("wH050" ,&wHig[2] ,"wH050/F" );
  mc.SetBranchAddress("trueCTL",&trueCTL);
  mc.SetBranchAddress("trueCTV",&trueCTV);
  mc.SetBranchAddress("trueChi",&trueChi);
  mc.SetBranchAddress("trueQ2", &trueQ2);
  mc.SetBranchAddress("MCType", &MCType);
  mc.SetBranchAddress("isDgamma", &isDgamma);

  if(nEntries<0 || nEntries > mc.GetEntries()) nEntries = mc.GetEntries();
  for(int entry=0; entry<nEntries; entry++){
    if(entry%(nEntries/5)==0) cout<<"Done entry "<<entry<<" of "<<nEntries<<endl;
    mc.GetEvent(entry);
    bool keep = true;
    int iMCType = MCType; if(iMCType>6) iMCType -= 6;
    for(int iHig=0; iHig<nHiggs; iHig++){
      if(iMCType==5) {
	wHig[iHig] = Dtaunu[2][iHig].FromSP8ToThisModel(trueQ2, PI-trueCTL, BToDtaunu::mTau);
      } else if(iMCType==1) {
	wHig[iHig] = Dtaunu[0][iHig].FromSP8ToThisModel(trueQ2, PI-trueCTL, BToDtaunu::mE);
      } else if(iMCType==3) {
	wHig[iHig] = Dtaunu[1][iHig].FromSP8ToThisModel(trueQ2, PI-trueCTL, BToDtaunu::mMu);
      } else if(iMCType==6) {
	if(!tag.Contains("gam")) 
	  wHig[iHig] = Dstaunu[2][iHig].FromSP8ToThisModel(trueQ2,-trueCTL,trueCTV,trueChi,isDgamma,true,BToDtaunu::mTau);
	else wHig[iHig] = Dstaunu[2][iHig].FromSP8ToThisModel(trueQ2,-trueCTL,trueCTV,trueChi,0,true,BToDtaunu::mTau);
      } else if(iMCType==2) {
	if(!tag.Contains("gam")) 
	  wHig[iHig] = Dstaunu[0][iHig].FromSP8ToThisModel(trueQ2,-trueCTL,trueCTV,trueChi,isDgamma,true,BToDtaunu::mE);
	else wHig[iHig] = Dstaunu[0][iHig].FromSP8ToThisModel(trueQ2,-trueCTL,trueCTV,trueChi,0,true,BToDtaunu::mE);
      } else if(iMCType==4) {
	if(!tag.Contains("gam")) 
	  wHig[iHig] = Dstaunu[1][iHig].FromSP8ToThisModel(trueQ2,-trueCTL,trueCTV,trueChi,isDgamma,true,BToDtaunu::mMu);
	else wHig[iHig] = Dstaunu[1][iHig].FromSP8ToThisModel(trueQ2,-trueCTL,trueCTV,trueChi,0,true,BToDtaunu::mMu);
      } else keep = false;
      if(wHig[iHig]!=wHig[iHig]) keep = false;
    }
    if(keep) {
      //tree.Fill();
      for(int iHig=0; iHig<nHiggs; iHig++){
	if(iMCType==1 ||iMCType==3 ||iMCType==5){
	  sumW[iHig][0][0] += wHig[iHig]; sumW[iHig][1][0] += wHig[iHig]*wHig[iHig]; 
	  nW[iHig][0]++;
	} else {
	  sumW[iHig][0][1+isDgamma] += wHig[iHig]; sumW[iHig][1][1+isDgamma] += wHig[iHig]*wHig[iHig]; 
	  nW[iHig][1+isDgamma]++;
	}
      }
    }else nBad++;
  }
//   TString treeName = "AWG82/generator/TreeWeight"; treeName+=tag; treeName+=".root";
//   TFile f(treeName,"RECREATE");
//   f.cd();
//   tree.Write();
//   f.Write();
//   f.Close();
//   cout<<endl<<"Written "<<treeName<<endl<<endl;
  if(nBad>0) cout<<"There were "<<nBad<<" bad entries."<<endl;
  cout<<endl<<"Entries: \t   ";
  for(int type=0; type<3; type++){
    cout<<Types[type]<<":  "<<nW[0][type]<<"\t";
    if(type==1)cout<<"\t";
  }
  cout<<endl;
  for(int iHig=0; iHig<nHiggs; iHig++){
    cout<<"For tBmH = "<<RoundNumber(tBmH[iHig],2)<<":   ";
    for(int type=0; type<3; type++){
      double N = (double)nW[iHig][type];
      double mean = sumW[iHig][0][type]/N;
      //double rms = sqrt(fabs(sumW[iHig][1][type]/pow(sumW[iHig][0][type],2)-1./N))*mean;
      double rms = sqrt(fabs((sumW[iHig][1][type]*N-pow(sumW[iHig][0][type],2)))/pow(N,3));
      cout<<RoundNumber(mean*100-100,3)<<" +- "<<RoundNumber(rms*100,3)<<" \t";
    }
    cout<<endl;
  }
}

