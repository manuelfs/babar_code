//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: MultiEventMixer.cc,v 1.2 2012/08/23 02:22:18 manuelf Exp $
//
// Description:
//      ultiEvenMixer - Prepares multiple file for the global fit, storing the right
//      mmiss and candtype
//
// Author List:
//      Michael Mazur (original)                  UC Santa Barbara
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      10/09/01 manuelf -- Adapted EvenMixer.cc
//------------------------------------------------------------------------

#include "TCut.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "DonutUtils/cuts.cc"
#include "DonutUtils/KeysUtils.cc"
#include <fstream>
#include <iostream>
#include <ctime>

using std::cout;
using std::endl;
#define nBDT 9


int main(int argc, char *argv[])
{
  if (argc < 2 || argc > 7) {
    cout << "USAGE: MultiEventMixer RootFile [NameTag=""] [scale2Data_cont_Dpipi_mES=yes_YES_mES] "
	 << "[tBmH=0] [emu=both] [weightFile]" << endl;
    return 1;
  }

  float candM2=0, candPstarD=0,candMES=0, candPstarLep=0,candPLep=0, mm2pi0=0, candMvaDl=0, weight=-1;
  float candMvaDssDl=0, candMvaDssComb=0, MCRD=0, MCRDs=0, candQ2=0, candEExtra=0, candPMiss=0, candThetaLep=0;
  float candDmass=0, candDeltam=0, trueCTL=0, trueCTV=0, trueChi=0, trueQ2=0, candM2Tru=0, truePLep=0;
  float candW=0, mD[] = {1.865, 2.007, 1.869, 2.010}, mB = 5.279;
  int MCType=0, candIsMu=0, candType=0, MCTaumode=0, MCD=0, candDstarType=0,MCPions=0, MCSubmode=0;

//   double entr;
//   TTree *tre = WeightedTree2("AWG82/ntuples/small/Data_RunAll.root", entr, "babar_code/Reweight/wTotal.txt", 0);
//   tre->SetBranchAddress("MCType",&MCType);
//   //tre->SetBranchAddress("candPstarLep",&pl);
//   for(int i=0; i<5; i++){
//     tre->GetEvent(i);
//     cout<<i<<" "<<MCType<<endl;
//   }
//   return 1;

  time_t start,end; double dif;
  time (&start);
  TString RootFile = argv[1];
  TString NameTag = "";
  if (argc>2) NameTag = argv[2];
  TString scale2Data = "yes_YES_mES";
  if (argc>3) scale2Data = argv[3];
  TString tBmH = "";
  if (argc>4) tBmH = argv[4];
  TString emu = "";
  if (argc>5) emu = argv[5];
  TString weightName = "babar_code/Reweight/wTotal.txt";
  if (argc>6) weightName = argv[6];

  int nIsmu = 1, isHiggs = 0;
  double wHiggsPsl[51][2];
  if(emu!="") nIsmu = 2; 
  if(tBmH != "") {
    isHiggs = 1;
    for(int isDs=0; isDs<2; isDs++){
      TString PslName = "babar_code/Reweight/HiggsPsl/wHiggsPsl_"; PslName += isDs;
      PslName += "_"; PslName += tBmH; PslName += ".txt";
      fstream textFile;
      textFile.open(PslName,fstream::in);
      for(int bin=1; bin<=50; bin++) textFile>>wHiggsPsl[bin][isDs];
      textFile.close();
    }
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
  TString nameFolder = RootFile;
  RootFile.Remove(0,RootFile.Last('/')+1);
  nameFolder.Remove(nameFolder.Last('/')+1,nameFolder.Length());
  nameFolder += "Fit"; 
  TString emu_s[] = {"e","mu"};

  TString BDTNames[] = {"030", "050", "075", "100", "120", "150", "200", "250", "300"};
  float BDTCuts[nBDT][4] = {{0.7207, 0.6110, 0.6000, 0.5850},     //0.3x 
			    {0.6300, 0.5650, 0.4740, 0.5250},     //0.5x
			    {0.5413, 0.5358, 0.4320, 0.4870},     //0.75x
			    {0.4580, 0.4800, 0.3600, 0.4100},     //1.0x
			    {0.3810, 0.4546, 0.3080, 0.3659},     //1.2x
			    {0.3065, 0.4133, 0.2530, 0.3171},     //1.5x
			    {0.1709, 0.3440, 0.1787, 0.2650},     //2.0x
			    {0.0990, 0.2804, 0.1275, 0.2240},     //2.5x
			    {0.0638, 0.2350, 0.0889, 0.1940}};    //3.0x
//   float BDTCuts[nBDT][4] = {{0.7207, 0.4800, 0.6000, 0.4100},     //0.3x 
// 			    {0.6300, 0.4800, 0.4740, 0.4100},     //0.5x
// 			    {0.5413, 0.4800, 0.4320, 0.4100},     //0.75x
// 			    {0.4580, 0.4800, 0.3600, 0.4100},     //1.0x
// 			    {0.3810, 0.4800, 0.3080, 0.4100},     //1.2x
// 			    {0.3065, 0.4800, 0.2530, 0.4100},     //1.5x
// 			    {0.1709, 0.4800, 0.1787, 0.4100},     //2.0x
// 			    {0.0990, 0.4800, 0.1275, 0.4100},     //2.5x
// 			    {0.0638, 0.4800, 0.0889, 0.4100}};    //3.0x
//   float BDTCuts[nBDT][4] = {{0.4580, 0.6110, 0.3600, 0.5850},     //0.3x 
// 			    {0.4580, 0.5650, 0.3600, 0.5250},     //0.5x
// 			    {0.4580, 0.5358, 0.3600, 0.4870},     //0.75x
// 			    {0.4580, 0.4800, 0.3600, 0.4100},     //1.0x
// 			    {0.4580, 0.4546, 0.3600, 0.3659},     //1.2x
// 			    {0.4580, 0.4133, 0.3600, 0.3171},     //1.5x
// 			    {0.4580, 0.3440, 0.3600, 0.2650},     //2.0x
// 			    {0.4580, 0.2804, 0.3600, 0.2240},     //2.5x
// 			    {0.4580, 0.2350, 0.3600, 0.1940}};    //3.0x
  float DssBDTCuts[2][4] = {{-0.45, -0.45, -0.45, -0.45}, {-0.35, -0.35, -0.3, -0.35}};
  TCut Cuts[2] = {basic, dss};
  Cuts[0] += M2P; Cuts[0] += "candMvaDl>0.09";
  Cuts[1] +=  cosT; Cuts[1] +=  Q2; Cuts[1] += Mpi0; Cuts[1] += "candMvaDssDl>-0.45"; Cuts[1] += "candMvaDssComb>-0.4";
  if(isHiggs){
    Cuts[0] = MvaAll;
    Cuts[1] = dssMvaAll;
  }
  int nTrees = 3;
  if(scale2Data.Contains("Dpipi")) nTrees = 4;

  float mESWeight[nBDT][4][2];
  TString mESBaseName = "babar_code/Reweight/wBkg_CombAllmESx", dummy;
  TTree *outputTree[nBDT][2];
  int begBDT = 3, endBDT = 4;
  if(scale2Data.Contains("All")) {begBDT = 0; endBDT = nBDT;}
  for(int iBDT = begBDT; iBDT < endBDT; iBDT++){
    for(int isMu = 0; isMu<nIsmu; isMu++){
      outputTree[iBDT][isMu] = new TTree("ntp1","cands");
      outputTree[iBDT][isMu] -> Branch("weight",&weight,"weight/F");
      outputTree[iBDT][isMu] -> Branch("candM2",&candM2,"candM2/F");
      outputTree[iBDT][isMu] -> Branch("candMES",&candMES,"candMES/F");
      outputTree[iBDT][isMu] -> Branch("candPstarLep",&candPstarLep,"candPstarLep/F");
      outputTree[iBDT][isMu] -> Branch("candPLep",&candPLep,"candPLep/F");
      outputTree[iBDT][isMu] -> Branch("candMvaDl",&candMvaDl,"candMvaDl/F");
      outputTree[iBDT][isMu] -> Branch("candPstarD",&candPstarD,"candPstarD/F");
      outputTree[iBDT][isMu] -> Branch("candType",&candType,"candType/I");
      outputTree[iBDT][isMu] -> Branch("candDstarType",&candDstarType,"candDstarType/I");
      outputTree[iBDT][isMu] -> Branch("MCType",&MCType,"MCType/I");
      outputTree[iBDT][isMu] -> Branch("MCD",&MCD,"MCD/I");
      outputTree[iBDT][isMu] -> Branch("MCSubmode",&MCSubmode,"MCSubmode/I");
      outputTree[iBDT][isMu] -> Branch("MCPions",&MCPions,"MCPions/I");
      outputTree[iBDT][isMu] -> Branch("MCRD",&MCRD,"MCRD/F");
      outputTree[iBDT][isMu] -> Branch("MCRDs",&MCRDs,"MCRDs/F");
      outputTree[iBDT][isMu] -> Branch("MCTaumode",&MCTaumode,"MCTaumode/I");
      outputTree[iBDT][isMu] -> Branch("candIsMu",&candIsMu,"candIsMu/I");
      outputTree[iBDT][isMu] -> Branch("candQ2",&candQ2,"candQ2/F");
      outputTree[iBDT][isMu] -> Branch("candW",&candW,"candW/F");
      outputTree[iBDT][isMu] -> Branch("candDmass",&candDmass,"candDmass/F");
      outputTree[iBDT][isMu] -> Branch("candDeltam",&candDeltam,"candDeltam/F");
      outputTree[iBDT][isMu] -> Branch("candEExtra",&candEExtra,"candEExtra/F");
      outputTree[iBDT][isMu] -> Branch("candPMiss",&candPMiss,"candPMiss/F");
      outputTree[iBDT][isMu] -> Branch("candThetaLep",&candThetaLep,"candThetaLep/F");
      outputTree[iBDT][isMu] -> Branch("trueCTL",&trueCTL,"trueCTL/F");
      outputTree[iBDT][isMu] -> Branch("trueCTV",&trueCTV,"trueCTV/F");
      outputTree[iBDT][isMu] -> Branch("trueChi",&trueChi,"trueChi/F");
      outputTree[iBDT][isMu] -> Branch("trueQ2" ,&trueQ2 ,"trueQ2/F" );
      outputTree[iBDT][isMu] -> Branch("candM2Tru",&candM2Tru,"candM2Tru/F");
      outputTree[iBDT][isMu] -> Branch("truePLep",&truePLep,"truePLep/F");

      TString mESName = mESBaseName; mESName += BDTNames[iBDT]; 
      if(nIsmu==2) mESName += emu_s[isMu];
      fstream mESFile; mESFile.open(mESName,fstream::in);
      for(int cand = 0; cand<4; cand++) mESFile>>dummy>>mESWeight[iBDT][cand][isMu];
      mESFile.close();
    }
  }
  for(int cut=0; cut<2; cut++){
    Cuts[cut] += MEScut;
    for(int tree=0; tree<nTrees; tree++){
      if((RootFile.Contains("Data") ||scale2Data.Contains("NO")) && tree>0) continue;
      double dEntries, gSR=0;
      int iEntries, isCocktail=0;
      if(tree>0) isCocktail=-1;
      if(isHiggs) gSR = -4.2*pow(tBmH.Atof()/100.,2);
      TString treeName = "tree"; treeName += cut; treeName += tree; 
      //TTree *inputTree =  new TTree(treeName,"ntp1");
      TString rootFile = "AWG82/ntuples/temp/"; rootFile += NameTag; rootFile += scale2Data; 
      rootFile += tBmH; rootFile += "temp.root";
      WeightedTree3(NameTrees[tree], dEntries, weightName,gSR,Cuts[cut],"All", rootFile);
      TChain *inputTree = new TChain("ntp1");
      inputTree->Add(rootFile);
      time (&end);dif = difftime (end,start);
      //cout<<dif<<" sec for "<<NameTrees[tree]<<endl;
      time (&start);
      //cout<<"Loaded tree "<< NameTrees[tree]<<", dEntries "<< dEntries<<endl;
      inputTree -> SetBranchAddress("weight",&weight);
      inputTree -> SetBranchAddress("mm2pi0",&mm2pi0);
      inputTree -> SetBranchAddress("candM2",&candM2);
      inputTree -> SetBranchAddress("candMES",&candMES);
      inputTree -> SetBranchAddress("candPLep",&candPLep);
      inputTree -> SetBranchAddress("candPstarLep",&candPstarLep);
      inputTree -> SetBranchAddress("candPstarD",&candPstarD);
      inputTree -> SetBranchAddress("candType",&candType);
      inputTree -> SetBranchAddress("candDstarType",&candDstarType);
      inputTree -> SetBranchAddress("candDmass",&candDmass);
      inputTree -> SetBranchAddress("candDeltam",&candDeltam);
      inputTree -> SetBranchAddress("MCType",&MCType);
      inputTree -> SetBranchAddress("MCSubmode",&MCSubmode);
      inputTree -> SetBranchAddress("MCTaumode",&MCTaumode);
      inputTree -> SetBranchAddress("MCD",&MCD);
      inputTree -> SetBranchAddress("MCPions",&MCPions);
      inputTree -> SetBranchAddress("candIsMu",&candIsMu);
      inputTree -> SetBranchAddress("candMvaDl",&candMvaDl);
      inputTree -> SetBranchAddress("candMvaDssDl",&candMvaDssDl);
      inputTree -> SetBranchAddress("candMvaDssComb",&candMvaDssComb);
      inputTree -> SetBranchAddress("candQ2",&candQ2);
      inputTree -> SetBranchAddress("candEExtra",&candEExtra);
      inputTree -> SetBranchAddress("trueCTL",&trueCTL);
      inputTree -> SetBranchAddress("trueCTV",&trueCTV);
      inputTree -> SetBranchAddress("trueChi",&trueChi);
      inputTree -> SetBranchAddress("trueQ2", &trueQ2);
      inputTree -> SetBranchAddress("candThetaLep",&candThetaLep);
      inputTree -> SetBranchAddress("candPMiss",&candPMiss);
      inputTree -> SetBranchAddress("candM2Tru",&candM2Tru);
      inputTree -> SetBranchAddress("truePLep",&truePLep);
      iEntries = inputTree->GetEntries();
      for(int entry=0; entry<iEntries; entry++){
	inputTree -> GetEvent(entry);

	int icandType = candType-1;
	MCRD = 0.3; MCRDs = 0.25;
	//MCRD = 0.449; MCRDs = 0.323;
	if(!NameTrees[tree].Contains("GenDss")) weight *= wMC;
	else {
	  if(MCType!=14) continue;
	  if(!RootFile.Contains("Data")  && scale2Data.Contains("yes")) {
	    if(MCType==14) weight *= 0.01/(totDss/totdata);
	    else if(MCType==15) weight *= 0.01/(totDpipi*0.25/totdata);
	    else if(MCType==16) weight *= 0.01/(totDpipi*0.75/totdata);
	    else cout<<"Error. In GenDss we have MCType "<<MCType<<endl;
	  } 
	  if(tree==0 && candType>2) candType -= 2;
	  MCType = 16; 
	}
	if(cut==1) {
	  candM2 = mm2pi0;
	  candType += 4;
	}
	if(isHiggs && MCTaumode>0 && fabs(candPLep-truePLep)<0.4){
	  int isDs = -1;
	  if(MCType==5 || MCType==11) isDs = 0;
	  if(MCType==6 || MCType==12) isDs = 1;
	  int iBin = (int)(candPstarLep/2*50.); // The re-weighting is done with 50 bins in the [0, 2] GeV range
	  if(isDs>=0 && iBin>=0 && iBin<50) weight *= wHiggsPsl[iBin][isDs];
	}
	if(!RootFile.Contains("Data")  && scale2Data.Contains("Dsstau") && MCType==14 && MCTaumode>=0) 
	  weight *= 1.50; // 50% systematic on R(D**)
	if(!RootFile.Contains("Data")  && scale2Data.Contains("Dsspi") && MCPions==1) 
	  weight *= 1.05; // 5% systematic on pi/pi0 efficiencies
	if(!RootFile.Contains("Data")  && scale2Data.Contains("M2tail") && candM2>2.5 && candType<=4 && (MCType>=1&&MCType<=4 || MCType>=7&&MCType<=10)) 
	  weight *= 1.05; // 5% systematic on the mmiss tail of normalization
	candW = (mB*mB+mD[icandType]*mD[icandType]-candQ2)/(2*mB*mD[icandType]);
	if(tree==1 || tree==2) weight *= wuds;
	if(RootFile.Contains("Data")) weight = 1.;
	int emuTree = candIsMu;
	if(nIsmu==1) emuTree = 0;
	double fixWeight = weight;
	for(int iBDT = begBDT; iBDT < endBDT; iBDT++) {
	  if(cut==0&&candMvaDl>BDTCuts[iBDT][icandType] || 
	     cut==1&&candMvaDl<=BDTCuts[iBDT][icandType]&&candMvaDssDl>DssBDTCuts[0][icandType]&&candMvaDssComb>DssBDTCuts[1][icandType]) {
	    if((MCType<=0||MCType>=13) && cut==0 && scale2Data.Contains("mES") && !RootFile.Contains("Data")) {
	      if(NameTag.Contains("on") && icandType%2 == 1) weight = fixWeight/mESWeight[4][icandType][emuTree];
	      else                       weight = fixWeight/mESWeight[iBDT][icandType][emuTree];
	    }
	    outputTree[iBDT][emuTree]->Fill();
	  }
	}
      }
      inputTree->Delete();
    }
  }

  for(int iBDT = begBDT; iBDT < endBDT; iBDT++){
    for(int isMu = 0; isMu<nIsmu; isMu++){
      TString BDTTag = NameTag; BDTTag += "x"; 
      if(isHiggs==0) BDTTag += BDTNames[iBDT]; 
      else BDTTag += tBmH;
      BDTTag += "_";
      TString nameFile = nameFolder; if(emu!="") nameFile += emu_s[isMu];
      nameFile += RootFile;
      nameFile.ReplaceAll("_", BDTTag);
      TFile f(nameFile,"RECREATE");
      f.cd();
      outputTree[iBDT][isMu]->Write();
      f.Write();
      f.Close();
      outputTree[iBDT][isMu]->Delete();
      cout<<"Written "<<nameFile<<endl;
    }
  }
}

