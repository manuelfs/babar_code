//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: EventMixer.cc,v 1.7 2012/03/04 00:33:02 manuelf Exp $
//
// Description:
//      EvenMixer - Prepares a file for the global fit, storing the right
//      mmiss and candtype
//
// Author List:
//      Michael Mazur (original)                  UC Santa Barbara
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      10/11/04 manuelf -- Based it on KeysUtils
//      10/01/18 manuelf -- Adapted to the new code
//------------------------------------------------------------------------

#include "TCut.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "DonutUtils/cuts.cc"
#include "DonutUtils/KeysUtils.cc"
#include <fstream>
#include <iostream>

using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  if (argc < 2 || argc > 5) {
    cout << "USAGE: EventMixer RootFile [scale2Data_cont_Dpipi=yes_YES] [emu=both] [weightFile]" << endl;
    return 1;
  }

  TString RootFile = argv[1];
  TString scale2Data = "yes";
  if (argc>2) scale2Data = argv[2];
  TString emu = "";
  if (argc>3) emu = argv[3];
  TString weightName = "babar_code/Reweight/wTotal.txt";
  if (argc>4) weightName = argv[4];

  int nIsmu = 1;
  if(emu!="") nIsmu = 2;
  float candM2, candPstarD,candMES, candPstarLep,candPLep, mm2pi0, weight=-1, MCRD = 0.3, MCRDs = 0.25;
  float candPisoftP=-99.;
  int MCType, candIsMu, candType, MCTaumode, MCD, candDstarType,MCPions;
  TTree *outputTree[2];
  for(int isMu = 0; isMu<nIsmu; isMu++){
    outputTree[isMu] = new TTree("ntp1","cands");
    outputTree[isMu] -> Branch("weight",&weight,"weight/F");
    outputTree[isMu] -> Branch("candM2",&candM2,"candM2/F");
    outputTree[isMu] -> Branch("candMES",&candMES,"candMES/F");
    outputTree[isMu] -> Branch("candPstarLep",&candPstarLep,"candPstarLep/F");
    outputTree[isMu] -> Branch("candPLep",&candPLep,"candPLep/F");
    outputTree[isMu] -> Branch("candPstarD",&candPstarD,"candPstarD/F");
    outputTree[isMu] -> Branch("candType",&candType,"candType/I");
    outputTree[isMu] -> Branch("candDstarType",&candDstarType,"candDstarType/I");
    outputTree[isMu] -> Branch("candPisoftP",&candPisoftP,"candPisoftP/F");
    outputTree[isMu] -> Branch("MCType",&MCType,"MCType/I");
    outputTree[isMu] -> Branch("MCD",&MCD,"MCD/I");
    outputTree[isMu] -> Branch("MCPions",&MCPions,"MCPions/I");
    outputTree[isMu] -> Branch("MCRD",&MCRD,"MCRD/F");
    outputTree[isMu] -> Branch("MCRDs",&MCRDs,"MCRDs/F");
    outputTree[isMu] -> Branch("MCTaumode",&MCTaumode,"MCTaumode/I");
    outputTree[isMu] -> Branch("candIsMu",&candIsMu,"candIsMu/I");
  }
  TString Runs = RootFile;
  Runs.Remove(0,Runs.Last('_')+1); Runs.ReplaceAll("Run",""); Runs.ReplaceAll(".root","");
  double totMCB = 0, totuds = 0, totccbar = 0, totdata = 0, totOffdata = 0;
  double totDpipi = 6874000+6799000, totDss = 7049000+7036000;
  getNumberB(RootFile, Runs, totMCB, totdata, totuds, totccbar, totOffdata);
  double wuds = totMCB/totuds*2.09/1.05;     
  double wMC = totdata/totMCB; if(!scale2Data.Contains("yes")) wMC = 1;
  if(RootFile.Contains("_4_Run") || RootFile.Contains("_14_Run")) wuds /= 2.;
  TString NameTrees[4] = {RootFile, "AWG82/ntuples/small/uds_RunAll.root", "AWG82/ntuples/small/ccbar_RunAll.root",
			  "AWG82/ntuples/small/GenDss_RunAll.root"};
  TCut Cuts[2] = {MvaAll, dssMvaAll};
  int nTrees = 3;
  if(scale2Data.Contains("Dpipi")) nTrees = 4;
  for(int cut=0; cut<2; cut++){
    Cuts[cut] += MEScut;
    for(int tree=0; tree<nTrees; tree++){
      if((RootFile.Contains("Data") ||scale2Data.Contains("NO")) && tree>0) continue;
      double dEntries;
      int iEntries, isCocktail=0;
      if(tree>0) isCocktail=-1;
      TTree *inputTree = WeightedTree(NameTrees[tree], dEntries, weightName,isCocktail,Cuts[cut]);
      inputTree -> SetBranchAddress("weight",&weight);
      inputTree -> SetBranchAddress("mm2pi0",&mm2pi0);
      inputTree -> SetBranchAddress("candM2",&candM2);
      inputTree -> SetBranchAddress("candMES",&candMES);
      inputTree -> SetBranchAddress("candPLep",&candPLep);
      inputTree -> SetBranchAddress("candPstarLep",&candPstarLep);
      inputTree -> SetBranchAddress("candPstarD",&candPstarD);
      inputTree -> SetBranchAddress("candPisoftP",&candPisoftP);
      inputTree -> SetBranchAddress("candType",&candType);
      inputTree -> SetBranchAddress("candDstarType",&candDstarType);
      inputTree -> SetBranchAddress("MCType",&MCType);
      inputTree -> SetBranchAddress("MCTaumode",&MCTaumode);
      inputTree -> SetBranchAddress("MCD",&MCD);
      inputTree -> SetBranchAddress("MCPions",&MCPions);
      inputTree -> SetBranchAddress("candIsMu",&candIsMu);
      iEntries = inputTree->GetEntries();
      for(int entry=0; entry<iEntries; entry++){
	inputTree -> GetEvent(entry);
	if(!NameTrees[tree].Contains("GenDss")) weight *= wMC;
	else {
	  if(MCType!=14) continue;
	  if(scale2Data.Contains("yes")) {
	    if(MCType==14) weight *= 0.01/(totDss/totdata);
	    else if(MCType==15) weight *= 0.01/(totDpipi*0.25/totdata);
	    else if(MCType==16) weight *= 0.01/(totDpipi*0.75/totdata);
	    else cout<<"Error. In GenDss we have MCType "<<MCType<<endl;
	  } 
	  if(tree==0 && candType>2) candType -= 2;
	  MCType = 16; 
	}
	if(scale2Data.Contains("Dsspi") && MCPions==1) MCType = 16;
	if(cut==1) {
	  candM2 = mm2pi0;
	  candType += 4;
	}
	if(tree==1 || tree==2) weight *= wuds;
	if(RootFile.Contains("Data")) weight = 1.;
	int emuTree = candIsMu;
	if(nIsmu==1) emuTree = 0;
	outputTree[emuTree] -> Fill();
      }
      inputTree->Delete();
    }
  }

  TString nameFolder = RootFile;
  RootFile.Remove(0,RootFile.Last('/')+1);
  nameFolder.Remove(nameFolder.Last('/')+1,nameFolder.Length());
  nameFolder += "Fit"; 
  nameFolder.ReplaceAll("Defsmall","small");
  TString emu_s[] = {"e","mu"}, nameFile;
  for(int isMu = 0; isMu<nIsmu; isMu++){
    nameFile = nameFolder; if(emu!="") nameFile += emu_s[isMu];
    nameFile += RootFile;
    TFile f(nameFile,"RECREATE");
    f.cd();
    outputTree[isMu]->Write();
    f.Write();
    f.Close();
  }
  cout<<"Written "<<nameFile<<endl;
}

