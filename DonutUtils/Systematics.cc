//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: Systematics.cc,v 1.2 2012/08/23 02:22:18 manuelf Exp $
//
// Description:
//      Systematics - Creates 1000 small ntuples varying the B combinatoric
//      and D**lnu BF within errors, as well as the D(*)lnu FF.
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      11/06/02 manuelf -- Created
//------------------------------------------------------------------------

#include "TCut.h"
#include "TString.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TFile.h"
#include "TMatrixTSym.h"
#include "TF3.h"
#include "DonutUtils/cuts.cc"
#include "DonutUtils/KeysUtils.cc"
#include <fstream>
#include <iostream>

using std::cout;
using std::endl;

double SBmES_Corr[2][2] = {{1,1},{0.0189, 0.0339}};
void EventMixRAll(TString wFile, TString ntupleName);
void GenPdfRAll(TString ntupleName, int number);
void fillTree(TTree *inputTree, TTree *outputTree, Int_t thecut, Int_t dss);

int main(int argc, char *argv[]){

  if (argc < 1 || argc > 3) {
    cout << "USAGE: Systematics nBegin nFiles"<< endl;
    return 1;
  }
  TString temp = argv[1]; int nBegin = temp.Atoi();
  temp = argv[2];         int nFiles = temp.Atoi();

  TString GenFile  = "AWG82/ntuples/small/RAll_RunAll.root";
  double totMCB = 0, totuds = 0, totccbar = 0, totdata = 0, totOffdata = 0;
  getNumberB(GenFile, "All", totMCB, totdata, totuds, totccbar, totOffdata);
  double wMC = totdata/totMCB; 

  // Branching fractions and errors from SP8 and PDG 2010
  int nBF = 49;
  TString BFNames[] = {"weightD10", "weightD20", "weightD00", "weightD1prime0", 
		       "weightD1p", "weightD2p", "weightD0p", "weightD1primep", "wDoubleSL", 
		       "wDsTaunu", "wDsmunu ", "wDsEtalnu", "wDsEtaprimelnu", "wDsphilnu", 
		       "wDzbDs  ", "wDzbDsstar", "wDstarzbDs", "wDstarzbDsstar", "wDstarstarDs", 
		       "wDpDzbKz", "wDzbDstarpKz", "wDstarzbDpKz", "wDstarzbDstarpKz", "wDzDzbKp", 
		       "wDstarzbDstarzKp", "wDstarzDzbKp", "wDstarmDpKp", "wDzba1p ", "wDstarzba1p", 
		       "wDzbRhop", "wDstarzbRhop", 
		       "wDmDs   ", "wDstarmDs", "wDmDsstar", "wDstarmDsstar", "wDmDzKp ", "wDmDstarzKp", 
		       "wDstarmDstarzKp", "wDmDpKz  ", "wDstarmDpKz", "wDstarmDstarpKz", "wDstarzDzbKz", 
		       "wDstarpDstarm", "wDpDstarm", "wDpDm   ", "wDma1p   ", "wDstarma1p", "wDmRhop ", "wDstarmRhop"};
  double BF[] = {43, 41, 41, 45, 40, 38, 38, 41, 100,
		 560, 58, 267, 99, 249, 
		 100, 76, 82, 171, 270, 14, 52, 6.3, 78, 21, 53, 47, 15, 40, 190, 134, 98,
		 72, 80, 74, 177, 17, 46, 118, 8.5, 65, 78, 18.5, 8.2, 6.1, 2.11, 60, 130, 76, 68};
  double BFSP[]={56, 30, 49, 90, 52, 23, 45, 83, 100,
		 640, 61.6, 307, 106, 242, 
		 129, 111, 124, 278, 136, 17, 49, 31, 100, 15, 70, 47, 20, 88.7, 159.7, 134, 98,
		 90, 126, 90, 240, 17, 49, 100, 15, 65, 70, 20, 8.3, 6.7, 2.7, 83.4, 120, 77, 68};
  double ErrorsBF[]= {3, 3, 8, 9, 5, 6, 8, 12, 10,
		      40, 4, 29, 23, 11, 
		      17, 16, 17, 24, 120, 14,12, 1.7, 26, 2.6, 16, 10, 4, 40, 50, 18, 17, 
		      8, 11, 16, 14, 40, 10, 20, 8.5, 16, 11, 18.5, 0.9, 1.5, 0.31, 33, 27, 13, 9};
  
  TString FFNames[] = {"rhoD2", "Delta", "HQET_R0", "HQET_rho2", "HQET_R1", "HQET_R2"};
  double FF[] =       {1.186, 1, 1.14, 1.207, 1.401, 0.854};    // From HFAG 2011, Tanaka 2010, and 1203.2654
  double ErrorsFF[] = {0.055, 1, 0.07, 0.028, 0.033, 0.020};

  double valFF[3][3] = {{FF[3], FF[4], FF[5]}, {ErrorsFF[3], ErrorsFF[4], ErrorsFF[5]},{0.575, -0.697, -0.872}};
  TMatrixT<double> CovFF(3,3); //CovFF is the covariance matrix
  for(int iFF=0; iFF<3; iFF++) {
    CovFF(iFF,iFF) = pow(valFF[1][iFF],2);
    CovFF(iFF,(iFF+1)%3) = valFF[2][iFF]*valFF[1][iFF]*valFF[1][(iFF+1)%3];
    CovFF((iFF+1)%3,iFF) = CovFF(iFF,(iFF+1)%3);
  }
  CovFF = Choleski(CovFF,3);
  TRandom3 rand(0); // 1 makes the first seed constant. 0 is a time dependent seed.



  TString TotalNames[] = {"FFWeights = ", "DlnFF = ", "BFWeights = ", 
			  "PlWeightsD  = ", 
			  "variWeights = babar_code/WeightsVari/Weights_", "BkgWeights  = "}; 


  //CSample noCut candEExtra 9 1.2 2.4 "candQ2>4" r no babar_code/Reweight/wPl.txt SUB
  double ContEntries[] = {13913, 8985, 11530, 1289}; // Continuum in 1.2<EExtra<2.4 and q2>4 (SUB = 23456)
  double dEntries = 53064; // Data entries minus semileptonic (SUB = 3456)
  //SidebandmES noCutSum candMES 20 5.2 5.26 "candQ2>4XXcandMvaDl>0.05XXcandM2>2" no no babar_code/Reweight/wTotal.txt 356789
  double dmES_SB[] = {2974.5, 1056.6}; // Data minus semileptonic in mES SB
  //CSample MVA candMES 9 5.2 5.26 "candQ2>4XXcandM2>2" r no babar_code/Reweight/wPl.txt 23456
  double ContmES_SB[] = {1190.57, 504.29, 568.20, 60.62}; // Continuum in 5.2<mES<5.26 and mmiss>2 

  fstream PlFileD;
  Float_t PlBins[1000], wPlD[4][1000], ErrorPl[] = {0.02968, 0.03609, 0.03110, 0.08859};
  PlFileD.open("babar_code/Reweight/wPlConti",fstream::in);
  for(int cand=0; cand<4; cand++){
    for(int i=0; i<1000; i++){
      PlFileD>>PlBins[i]>>wPlD[cand][i];
    }
  }
  PlFileD.close();
  TString folder = "AWG82/systematics/BF/Text/";

  // ====================================  Iterations  ==========================================
  for(int file=nBegin; file<=nFiles; file++){
    if(file%10==0) cout<<"Doing iteration "<<file<<" of "<<nFiles<<endl;

    // ==============  Variations on the background BF  ===================
    TString BFfileName = folder; BFfileName += "wAllBF_"; BFfileName += file;
    //BFfileName = "babar_code/Reweight/wAllBF.Nov";   // Provisional if default values are to be written
    fstream BFFile; BFFile.open(BFfileName,fstream::out);
    for(int w=0; w<nBF; w++){
      double BFGaus = rand.Gaus(BF[w], ErrorsBF[w]);
      //BFGaus = BF[w];   // Provisional if default values are to be written
      if(BFGaus<0) BFGaus = 0;
      BFFile<<BFNames[w]<<"\t= "<<RoundNumber(BFGaus,4,BFSP[w])<<endl;
      if(w==8 || w==13 || w==30) BFFile<<endl;
    }
    BFFile.close();
    //continue;   // Provisional if default values are to be written

    // ==============  Variations on the D(*)lnu FF  ===================
    TString FFfileName = folder; FFfileName += "wAllFF_"; FFfileName += file;
    fstream FFFile; FFFile.open(FFfileName,fstream::out);
    double RandFFGaus[6], fZ[3];
    for(int chan=0; chan<3; chan++) fZ[chan] = rand.Gaus(0, 1);
    for(int w=0; w<6; w++){
      if(w<3) RandFFGaus[w] = rand.Gaus(FF[w], ErrorsFF[w]);
      else {
	RandFFGaus[w] = FF[w];
	for(int col=0; col<3; col++) RandFFGaus[w] += fZ[col]*CovFF(w-3,col);
      }
      if(RandFFGaus[w]<0 && w!=1) RandFFGaus[w] = 0;
      FFFile<<FFNames[w]<<"\t= "<<RoundNumber(RandFFGaus[w],4)<<endl;
    }
    FFFile.close();

    // ===========  Writing the names of the re-weighting files  ================
    TString ContfileName = folder; ContfileName += "wPlConti_"; ContfileName += file; 
    TString CombfileName = folder; CombfileName += "wComb_"; CombfileName += file; 
    TString TotalfileName = folder; TotalfileName += "wTotal_"; TotalfileName += file; TotalfileName += ".txt";
    fstream TotalFile; TotalFile.open(TotalfileName,fstream::out);
    for(int w=0; w<5; w++){
      TotalFile<<TotalNames[w];
      if(w<=1) TotalFile<<FFfileName;
      if(w==2) TotalFile<<BFfileName;
      if(w==3) TotalFile<<ContfileName;
      TotalFile<<endl;
    }

    // ============  Variations on the Continuum Pl weights  =================
    double dEntriesSub = dEntries, wCandPl[4];
    fstream ContFile; ContFile.open(ContfileName,fstream::out);
    for(int cand=0; cand<4; cand++){
      wCandPl[cand] = rand.Gaus(1,ErrorPl[cand]); if(wCandPl<=0) wCandPl[cand] = 0.001;
      for(int i=0; i<1000; i++){
	ContFile<<PlBins[i]<<"\t"<<wPlD[cand][i]*wCandPl[cand]<<endl;
      }
      ContFile<<endl<<endl;
      dEntriesSub -= ContEntries[cand]*wCandPl[cand];
    }
    ContFile.close();
    // =========  Variations on the EExtra and mES corrections of the BB bkg  ==============
    double BkgEntries = 0, nBB_SBmES_D[] = {0,0};
    double wCombGaus = rand.Gaus(1,0.0186); if(wCombGaus<=0) wCombGaus = 0.;
    // c.CopyTree(basic+"MCType==0&&candEExtra>1.2&&candEExtra<2.4&&candMES>5.27&&candQ2>4")
    WeightedTree("AWG82/systematics/BB_SBEExtra.root", BkgEntries, TotalfileName,0,"1","All","BkgTree");
    fstream CombFile; CombFile.open(CombfileName,fstream::out);
    for(int w=0; w<4; w++){
      CombFile<<w<<"\t"<<RoundNumber(BkgEntries*wMC*wCombGaus,4,dEntriesSub)<<endl;
    }
    CombFile.close();
    TotalFile<<TotalNames[5]<<CombfileName;
    TotalFile.close();

    WeightedTree("AWG82/systematics/BB_SBmES_D.root", nBB_SBmES_D[0], TotalfileName,0,"1","All","DBkgTree");
    WeightedTree("AWG82/systematics/BB_SBmES_Ds.root", nBB_SBmES_D[1], TotalfileName,0,"1","All","DsBkgTree");
    for(int w=0; w<2; w++) {
      double totCont = ContmES_SB[w]*wCandPl[w] + ContmES_SB[w+2]*wCandPl[w+2];
      SBmES_Corr[0][w] = (totCont+nBB_SBmES_D[w]*wMC)/dmES_SB[w]*rand.Gaus(1,SBmES_Corr[1][w]);
      //cout<<SBmES_Corr[0][w]<<" \t nBB_SBmES_Dendl "<<nBB_SBmES_D[w]*wMC<<"\t totCont "<<totCont<<endl;
    }

    TString ntupleName = "AWG82/systematics/BF/FitRoot/Fit"; ntupleName += file; ntupleName += "RAll_RunAll.root";
    EventMixRAll(TotalfileName, ntupleName);
    GenPdfRAll(ntupleName, file);
  }

  cout<<endl<<"Made files from "<<nBegin<<" to "<<nFiles<<" at "<<folder<<endl<<endl;
  return 1;
}

void GenPdfRAll(TString ntupleName, int number){
  TChain *chain = new TChain("ntp1");
  chain->Add(ntupleName);

  TFile *f;
  TString pdfFolder = "AWG82/systematics/BF/fitSamples/Iter"; pdfFolder += number; 
  gSystem->mkdir(pdfFolder);
  for (int i = 1 ; i <= 58 ; i ++) {
    if(i==33) i=41; if(i==49) i=51;   //Skipping non-existent numbers
    int isDss = 1;
    TString pdfName = pdfFolder; pdfName += "/pdfSample";pdfName += i; pdfName += ".root";
    f = new TFile(pdfName,"RECREATE");
    TTree *t = new TTree("ntp1","cands");
    fillTree(chain,t,i,isDss);
    t->Write();
    t->ResetBranchAddresses();
    f->Write();
    f->Close();
    f->Delete();
    //t->Delete();  // It crashes if I try to delete it, but I think this might be leaking
  } 
  chain->Delete();
}

void fillTree(TTree *inputTree, TTree *outputTree, Int_t thecut, Int_t isDss){
  Float_t myCandM2;
  Int_t MCType,candType;
  Float_t candM2,candPstarLep,weight;

  inputTree -> SetBranchAddress("candPstarLep",  &candPstarLep);            
  inputTree -> SetBranchAddress("candM2",        &candM2);                  
  inputTree -> SetBranchAddress("candType",      &candType);                
  inputTree -> SetBranchAddress("MCType",        &MCType);                  
  inputTree -> SetBranchAddress("weight",        &weight);                       

  outputTree -> Branch("MCType",&MCType,"MCType/I");
  outputTree -> Branch("candPstarLep",&candPstarLep,"candPstarLep/F");
  outputTree -> Branch("candM2",&myCandM2,"candM2/F");
  outputTree -> Branch("candType",&candType,"candType/I");
  outputTree -> Branch("weight",&weight,"weight/F");

  bool Sep12 = false;
  Int_t n1 = 0, n2 = 0, nevt = (Int_t)inputTree -> GetEntries();
  double totWeight=0, totEvents=0;
  for (int j = 0 ; j < nevt ; j ++) {
    inputTree -> GetEvent(j);
    if(thecut>=35&&thecut<=36 || thecut>=39&&thecut<=40) candType += 2;
    if(!isSample(thecut, MCType, candType)) continue;
    totWeight += weight; totEvents++;
  }
  for (int j = 0 ; j < nevt ; j ++) {
    inputTree -> GetEvent(j);

    if(!isSample(thecut, MCType, candType)) continue;
    weight *= (totEvents/totWeight);
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
}

void EventMixRAll(TString weightName, TString ntupleName){
  TString RootFile = "AWG82/ntuples/small/RAll_RunAll.root";

  float candM2,candPstarLep, mm2pi0, weight=-1, MCRD, MCRDs;
  int MCType, candIsMu, candType, MCTaumode, MCD;
  TTree *outputTree = new TTree("ntp1","cands");
  outputTree -> Branch("weight",&weight,"weight/F");
  outputTree -> Branch("candM2",&candM2,"candM2/F");
  outputTree -> Branch("candPstarLep",&candPstarLep,"candPstarLep/F");
  outputTree -> Branch("candType",&candType,"candType/I");
  outputTree -> Branch("MCType",&MCType,"MCType/I");
  outputTree -> Branch("MCD",&MCD,"MCD/I");
  outputTree -> Branch("MCTaumode",&MCTaumode,"MCTaumode/I");
  outputTree -> Branch("candIsMu",&candIsMu,"candIsMu/I");
  outputTree -> Branch("MCRD",&MCRD,"MCRD/F");
  outputTree -> Branch("MCRDs",&MCRDs,"MCRDs/F");

  double totMCB = 0, totuds = 0, totccbar = 0, totdata = 0, totOffdata = 0;
  getNumberB(RootFile, "All", totMCB, totdata, totuds, totccbar, totOffdata);
  double wuds = totMCB/totuds*2.09/1.05;     
  double wMC = totdata/totMCB; 
  TString NameTrees[] = {RootFile, "AWG82/ntuples/small/uds_RunAll.root", "AWG82/ntuples/small/ccbar_RunAll.root"};
  TCut Cuts[2] = {MvaAll+MEScut, dssMvaAll+MEScut};
  for(int cut=0; cut<2; cut++){
    for(int tree=0; tree<3; tree++){
      double dEntries;
      int iEntries, isCocktail=0;
      if(tree>0) isCocktail=-1;
      TTree *inputTree = WeightedTree(NameTrees[tree], dEntries, weightName,isCocktail,Cuts[cut]);
      inputTree -> SetBranchStatus("*",0);  // It started crashing when I loaded all branches
      inputTree -> SetBranchStatus("weight",1);
      inputTree -> SetBranchStatus("mm2pi0",1);
      inputTree -> SetBranchStatus("candM2",1);
      inputTree -> SetBranchStatus("candPstarLep",1);
      inputTree -> SetBranchStatus("candType",1);
      inputTree -> SetBranchStatus("MCType",1);
      inputTree -> SetBranchStatus("MCTaumode",1);
      inputTree -> SetBranchStatus("MCD",1);
      inputTree -> SetBranchStatus("candIsMu",1);
      inputTree -> SetBranchAddress("weight",&weight);
      inputTree -> SetBranchAddress("mm2pi0",&mm2pi0);
      inputTree -> SetBranchAddress("candM2",&candM2);
      inputTree -> SetBranchAddress("candPstarLep",&candPstarLep);
      inputTree -> SetBranchAddress("candType",&candType);
      inputTree -> SetBranchAddress("MCType",&MCType);
      inputTree -> SetBranchAddress("MCTaumode",&MCTaumode);
      inputTree -> SetBranchAddress("MCD",&MCD);
      inputTree -> SetBranchAddress("candIsMu",&candIsMu);
      iEntries = inputTree->GetEntries();
      for(int entry=0; entry<iEntries; entry++){
	inputTree -> GetEvent(entry);
	MCRD = 0.3; MCRDs = 0.25;
	weight *= wMC;
	if(cut==1) {
	  candM2 = mm2pi0;
	  candType += 4;
	}
	if(tree==1 || tree==2) weight *= wuds;
	if((MCType<=0 || MCType>=13) && cut==0) weight /= SBmES_Corr[0][candType%2==0];
	outputTree -> Fill();
      }
      inputTree->Delete();
    }
  }

  TFile f(ntupleName,"RECREATE");
  f.cd();
  outputTree->Write();
  f.Write();
  f.Close();
  outputTree->Delete();
}

