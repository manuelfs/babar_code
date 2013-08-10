//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: SmallNtuple.cc,v 1.5 2010/05/11 19:00:18 manuelf Exp $
//
// Description:
//      SmallNtuple - Creates a flat ntuple, recalculating the BDT and
//      weights
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      09/07/18 manuelf -- Created
//------------------------------------------------------------------------

#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include "TMVA/Reader.h"
#include "DonutUtils/weightManager.hh"
#include "DonutUtils/PidCorrectMesMean.hh"
using std::cout;
using std::endl;

Long64_t LoadTree(Long64_t entry, TChain* fChain);

int main(int argc, char *argv[])
{
  if (argc < 3 || argc > 5 ) {
    cout << "USAGE: Folder RunNumber [tagRest] [weightFile]" << endl;
    return 1;
  }

  TChain c("ntp1");

  TString folder = argv[1];
  TString Run = argv[2];
  TString tagRest;
  if(argc ==4) tagRest = argv[3];
  else tagRest = "";
  TString weightName = "babar_code/Reweight/wFF.txt";
  if (argc>4) weightName = argv[4];

  Int_t nFiles = 0;
   if (folder != "") {
     TString directory = "AWG82/ntuples/";
     directory += folder;
     TString filename;
     if(Run == "All"){
       TString dirs = "/Merged"; dirs+=tagRest; dirs+="/*";
       filename=directory;
       filename += dirs;
       nFiles += c.Add(filename);
       cout<<"Added "<<nFiles<<" files from "<<filename<<" with entries: "<<c.GetEntries()<<endl;   
     }else{
       for(int i=1; i<7; i++){
	 TString irun = ""; irun += i;
	 if(Run.Contains(irun)){
	   TString dirs = "/Merged"; dirs+=tagRest; dirs+="/*Run"; dirs += irun; dirs+="*";
	   filename=directory;
	   filename += dirs;
	   nFiles += c.Add(filename);
	   cout<<"Added "<<nFiles<<" files from "<<filename<<" with entries: "<<c.GetEntries()<<endl;   
	 }
       }
     }
   } else {
     cout<<"Specify a folder name"<<endl;
     return 1;
   }
   bool mEScorrection = false;
   //if(folder.Contains("data") && (Run.Contains("5") || Run.Contains("6"))) mEScorrection = true;
   if(mEScorrection){
     cout<<"Doing the mES correction"<<endl;
     PidCorrectMesMean::initialize("../DonutUtils/mES_means_Runs123456_final.txt");
   }

  fstream textFile;
  textFile.open("babar_code/modesCocktail.txt",fstream::in);
  int modeMap[87][54], BMode, inCocktail, Dmode, Xmode;
  for(int i=0; i<87; i++)
    for(int j=0; j<54; j++)
      modeMap[i][j] = 0;
  while(textFile){
    textFile>>BMode;
    Dmode = (BMode/100)%100; Xmode = BMode%100;
    modeMap[Dmode][Xmode] = 1;
  }


  // Declaration of leave types
  Int_t runnum,platform,partition,upperID,lowerID,isSP6,MCType,MCSubmode,MCCombmode;
  Int_t MCPions,MCTaumode,MCBrem,MCScatter,MCD,MC2Body,MCDoubleSL,MCCombB,MCCombDs,nTrueP;
  Int_t nTrueKL,nTrueNu,isBzero,MCDssmode,MCUnsim,trueLepMother,MCComblong,MCComblongD;
  Float_t truePLep,truePD,trueDmass,trueCTL,trueCTV,trueChi,trueQ2,trueTauFlight;
  Float_t truePPi0,trueDssPPi0,truePLep1,truePLep2,truePLep3;
  Int_t trueLepCharge,ntrueLep,trueMotherLep1,trueMotherLep2,trueMotherLep3;
  Int_t trueLundLep1,trueLundLep2,trueLundLep3,trueTagPi,trueTagK,trueTagKs,trueTagPi0;
  Int_t trueTagE,trueTagMu,trueTagDs,trueTagD,trueTagYPi,trueTagYK,trueTagYKs,trueTagYPi0;
  Int_t trueTagElse,nUpsilon,nups1,nups2,nups3,nups4,ncands,ncands1,ncands2,ncands3,ncands4;
  Int_t candType,candDstarType,candDType,candIsMixed,candExtraTracks,candRejectedTracks,candRejectedPhotons;
  Int_t candIsMu,candLepCharge,candDLund,candLepTru,candDTru,candBTru,candFitStat,candFitNdof;
  Float_t eventM2,candDmass,candDeltam,candFitChi2,candFitProb,candPisoftP,candPisoftTheta;
  Float_t candPisoftPhi,candPisoftE1,candPisoftE2,candPisoftM2,candKsMass,candPi0Mass,candPi0P;
  Float_t candPLep,candThetaLep,candPstarLep,candPD,candPstarD,candEstarD,candBremE,candBremAngle;
  Float_t candEExtra,candUsedEnergy;
  Int_t candIsBrem,candExtraPhotons,candUsedPhotons,candBMode,candBntCha,candBntNeu;
  Int_t candDntCha,candDntNeu,candBnCharged,candBnNeutral,candDnCharged,candDnNeutral;
  Float_t candMExtra,candMTot,candMES,candMES_orig,candDeltaE,candBTagDeltaP,candBTagDeltaPx;
  Float_t candBTagDeltaPy,candBTagDeltaPz,candBTagDeltaE,candBTagYMass,candBTagDmass;
  Float_t candBTagDeltam,candDDeltaP,candDDeltaPx,candDDeltaPy,candDDeltaPz,candDDeltaE;
  Float_t candLepDeltaP,candLepDeltaPx,candLepDeltaPy,candLepDeltaPz,candLepDeltaE,candNuDeltaP;
  Float_t candNuDeltaPx,candNuDeltaPy,candNuDeltaPz,candNuDeltaE,candBDPx,candBDPy,candBDPz;
  Float_t candBDE,candMissPx,candMissPy,candMissPz,candMissE,candIntPur;
  Int_t  candTagLeptons;
  Int_t  candTagChargedMult,candTagNeutMult;
  Float_t emcfwd,emcbwd,trkfwd,trkbwd,candM2,candQ2,candPMiss,candEMiss,candThetaMiss;
  Float_t candPhiMiss,candThetaDL,candHelicity,candM2Tru,candM2TruBTag,candM2TruD,candM2TruLep;
  Float_t candMvaBestB,candMvaComb,candMvaDl,candCosT;
  Float_t         candMvaBestBL[390];   //[ncands]
  Int_t           candBModeL[390];   //[ncands]
  Float_t         candEExtraL[390];   //[ncands]
   Int_t           candBntChaL[390];   //[ncands]
   Int_t           candBntNeuL[390];   //[ncands]
   Int_t           candDntChaL[390];   //[ncands]
   Int_t           candDntNeuL[390];   //[ncands]
   Int_t           candTypeL[390];   //[ncands]
   Int_t           candLepTruL[390];   //[ncands]
   Float_t         candDeltaEL[390];   //[ncands]
   Int_t           npi0;
   Float_t         q2pi0[920];   //[npi0]
   Float_t         mpi0[920];   //[npi0]
   Float_t         ppi0[920];   //[npi0]
   Float_t         dmpi0[920];   //[npi0]
   Float_t         mm2pi0[920];   //[npi0]
   Float_t         pmisspi0[920];   //[npi0]
   Float_t         eextrapi0[920];   //[npi0]
   Float_t         e1pi0[920];   //[npi0]
   Float_t         e2pi0[920];   //[npi0]
   Int_t           bestepi0[920];   //[npi0]
   Int_t           bestmasspi0[920];   //[npi0]
   Int_t           ndpi0;
   Float_t         dmdpi0[68160];   //[ndpi0]
   Float_t         mm2dpi0[68160];   //[ndpi0]
   Float_t         pmissdpi0[68160];   //[ndpi0]
   Float_t         eextradpi0[68160];   //[ndpi0]
   Float_t         p1dpi0[68160];   //[ndpi0]
   Float_t         p2dpi0[68160];   //[ndpi0]
   Float_t         m1dpi0[68160];   //[ndpi0]
   Float_t         m2dpi0[68160];   //[ndpi0]
   Int_t           bestdpi0[68160];   //[ndpi0]
   Int_t           nGTLHard;
   Int_t           nGTLSoft;
   Int_t           nOtherTracks;
   Int_t           nCandK;
   Int_t           nCandPi;
   Float_t         candPKaons[20];   //[nCandK]
   Float_t         candPPions[30];   //[nCandPi]
   Float_t trkWeight,candKSdxy,candKSTheta,candKSPT,candPPi0,kaonWeight,pionWeight;
   Float_t elecWeight,muonWeight;

   c.SetBranchAddress("runnum", &runnum);                      
   c.SetBranchAddress("platform", &platform);                  
   c.SetBranchAddress("partition", &partition);                
   c.SetBranchAddress("upperID", &upperID);                    
   c.SetBranchAddress("lowerID", &lowerID);                    
   c.SetBranchAddress("isSP6", &isSP6);                        
   c.SetBranchAddress("MCType", &MCType);                      
   c.SetBranchAddress("MCSubmode", &MCSubmode);                
   c.SetBranchAddress("MCCombmode", &MCCombmode);              
   c.SetBranchAddress("MCPions", &MCPions);                    
   c.SetBranchAddress("MCTaumode", &MCTaumode);                
   c.SetBranchAddress("MCBrem", &MCBrem);                      
   c.SetBranchAddress("MCScatter", &MCScatter);                
   c.SetBranchAddress("MCD", &MCD);                            
   c.SetBranchAddress("MC2Body", &MC2Body);                    
   c.SetBranchAddress("MCDoubleSL", &MCDoubleSL);              
   c.SetBranchAddress("MCComblong", &MCComblong);                    
   c.SetBranchAddress("MCComblongD", &MCComblongD);                    
   c.SetBranchAddress("MCCombB", &MCCombB);                    
   c.SetBranchAddress("MCCombDs", &MCCombDs);                  
   c.SetBranchAddress("nTrueP", &nTrueP);                      
   c.SetBranchAddress("nTrueKL", &nTrueKL);                    
   c.SetBranchAddress("nTrueNu", &nTrueNu);                    
   c.SetBranchAddress("isBzero", &isBzero);                    
   c.SetBranchAddress("MCDssmode", &MCDssmode);                
   c.SetBranchAddress("MCUnsim", &MCUnsim);                    
   c.SetBranchAddress("truePLep", &truePLep);                  
   c.SetBranchAddress("trueLepMother", &trueLepMother);        
   c.SetBranchAddress("truePD", &truePD);                      
   c.SetBranchAddress("trueDmass", &trueDmass);                
   c.SetBranchAddress("trueCTL", &trueCTL);                    
   c.SetBranchAddress("trueCTV", &trueCTV);                    
   c.SetBranchAddress("trueChi", &trueChi);                    
   c.SetBranchAddress("trueQ2", &trueQ2);                      
   c.SetBranchAddress("trueLepCharge", &trueLepCharge);        
   c.SetBranchAddress("trueTauFlight", &trueTauFlight);        
   c.SetBranchAddress("truePPi0", &truePPi0);                  
   c.SetBranchAddress("trueDssPPi0", &trueDssPPi0);            
   c.SetBranchAddress("ntrueLep", &ntrueLep);                  
   c.SetBranchAddress("trueMotherLep1", &trueMotherLep1);      
   c.SetBranchAddress("trueMotherLep2", &trueMotherLep2);      
   c.SetBranchAddress("trueMotherLep3", &trueMotherLep3);      
   c.SetBranchAddress("trueLundLep1", &trueLundLep1);          
   c.SetBranchAddress("trueLundLep2", &trueLundLep2);          
   c.SetBranchAddress("trueLundLep3", &trueLundLep3);         
   c.SetBranchAddress("truePLep1", &truePLep1);                
   c.SetBranchAddress("truePLep2", &truePLep2);                
   c.SetBranchAddress("truePLep3", &truePLep3);                
   c.SetBranchAddress("trueTagPi", &trueTagPi);                
   c.SetBranchAddress("trueTagK", &trueTagK);                  
   c.SetBranchAddress("trueTagKs", &trueTagKs);                
   c.SetBranchAddress("trueTagPi0", &trueTagPi0);              
   c.SetBranchAddress("trueTagE", &trueTagE);                  
   c.SetBranchAddress("trueTagMu", &trueTagMu);                
   c.SetBranchAddress("trueTagDs", &trueTagDs);                
   c.SetBranchAddress("trueTagD", &trueTagD);                  
   c.SetBranchAddress("trueTagYPi", &trueTagYPi);              
   c.SetBranchAddress("trueTagYK", &trueTagYK);                
   c.SetBranchAddress("trueTagYKs", &trueTagYKs);              
   c.SetBranchAddress("trueTagYPi0", &trueTagYPi0);            
   c.SetBranchAddress("trueTagElse", &trueTagElse);            
   c.SetBranchAddress("nUpsilon", &nUpsilon);                  
   c.SetBranchAddress("nups1", &nups1);                        
   c.SetBranchAddress("nups2", &nups2);                        
   c.SetBranchAddress("nups3", &nups3);                        
   c.SetBranchAddress("nups4", &nups4);                        
   c.SetBranchAddress("ncands", &ncands);                      
   c.SetBranchAddress("ncands1", &ncands1);                    
   c.SetBranchAddress("ncands2", &ncands2);                    
   c.SetBranchAddress("ncands3", &ncands3);                    
   c.SetBranchAddress("ncands4", &ncands4);                    
   c.SetBranchAddress("candType", &candType);                  
   c.SetBranchAddress("candDstarType", &candDstarType);      
   c.SetBranchAddress("candDType", &candDType);                
   c.SetBranchAddress("candIsMixed", &candIsMixed);            
   c.SetBranchAddress("candExtraTracks", &candExtraTracks);   
   c.SetBranchAddress("candRejectedTracks", &candRejectedTracks);
   c.SetBranchAddress("candRejectedPhotons", &candRejectedPhotons);
   c.SetBranchAddress("eventM2", &eventM2);                    
   c.SetBranchAddress("candIsMu", &candIsMu);                  
   c.SetBranchAddress("candLepCharge", &candLepCharge);       
   c.SetBranchAddress("candDLund", &candDLund);                
   c.SetBranchAddress("candLepTru", &candLepTru);              
   c.SetBranchAddress("candDTru", &candDTru);                  
   c.SetBranchAddress("candBTru", &candBTru);                  
   c.SetBranchAddress("candDmass", &candDmass);                
   c.SetBranchAddress("candDeltam", &candDeltam);              
   c.SetBranchAddress("candFitStat", &candFitStat);            
   c.SetBranchAddress("candFitNdof", &candFitNdof);            
   c.SetBranchAddress("candFitChi2", &candFitChi2);            
   c.SetBranchAddress("candFitProb", &candFitProb);               
   c.SetBranchAddress("candPisoftP", &candPisoftP);               
   c.SetBranchAddress("candPisoftTheta", &candPisoftTheta);       
   c.SetBranchAddress("candPisoftPhi", &candPisoftPhi);           
   c.SetBranchAddress("candPisoftE1", &candPisoftE1);             
   c.SetBranchAddress("candPisoftE2", &candPisoftE2);             
   c.SetBranchAddress("candPisoftM2", &candPisoftM2);             
   c.SetBranchAddress("candKsMass", &candKsMass);                 
   c.SetBranchAddress("candPi0Mass", &candPi0Mass);               
   c.SetBranchAddress("candPi0P", &candPi0P);                     
   c.SetBranchAddress("candPLep", &candPLep);                     
   c.SetBranchAddress("candThetaLep", &candThetaLep);             
   c.SetBranchAddress("candPstarLep", &candPstarLep);             
   c.SetBranchAddress("candPD", &candPD);                         
   c.SetBranchAddress("candPstarD", &candPstarD);                 
   c.SetBranchAddress("candEstarD", &candEstarD);                 
   c.SetBranchAddress("candIsBrem", &candIsBrem);                 
   c.SetBranchAddress("candBremE", &candBremE);                   
   c.SetBranchAddress("candBremAngle", &candBremAngle);           
   c.SetBranchAddress("candEExtra", &candEExtra);                 
   c.SetBranchAddress("candExtraPhotons", &candExtraPhotons);     
   c.SetBranchAddress("candUsedEnergy", &candUsedEnergy);         
   c.SetBranchAddress("candUsedPhotons", &candUsedPhotons);       
   c.SetBranchAddress("candMExtra", &candMExtra);                 
   c.SetBranchAddress("candMTot", &candMTot);                     
   c.SetBranchAddress("candBMode", &candBMode);                   
   c.SetBranchAddress("candMES", &candMES);                       
   c.SetBranchAddress("candMES_orig", &candMES_orig);             
   c.SetBranchAddress("candDeltaE", &candDeltaE);                 
   c.SetBranchAddress("candBTagDeltaP", &candBTagDeltaP);         
   c.SetBranchAddress("candBTagDeltaPx", &candBTagDeltaPx);       
   c.SetBranchAddress("candBTagDeltaPy", &candBTagDeltaPy);       
   c.SetBranchAddress("candBTagDeltaPz", &candBTagDeltaPz);       
   c.SetBranchAddress("candBTagDeltaE", &candBTagDeltaE);         
   c.SetBranchAddress("candBTagYMass", &candBTagYMass);           
   c.SetBranchAddress("candBTagDmass", &candBTagDmass);           
   c.SetBranchAddress("candBTagDeltam", &candBTagDeltam);         
   c.SetBranchAddress("candBntCha", &candBntCha);                 
   c.SetBranchAddress("candBntNeu", &candBntNeu);                 
   c.SetBranchAddress("candDntCha", &candDntCha);                 
   c.SetBranchAddress("candDntNeu", &candDntNeu);                 
   c.SetBranchAddress("candBnCharged", &candBnCharged);           
   c.SetBranchAddress("candBnNeutral", &candBnNeutral);           
   c.SetBranchAddress("candDnCharged", &candDnCharged);           
   c.SetBranchAddress("candDnNeutral", &candDnNeutral);           
   c.SetBranchAddress("candDDeltaP", &candDDeltaP);               
   c.SetBranchAddress("candDDeltaPx", &candDDeltaPx);             
   c.SetBranchAddress("candDDeltaPy", &candDDeltaPy);             
   c.SetBranchAddress("candDDeltaPz", &candDDeltaPz);             
   c.SetBranchAddress("candDDeltaE", &candDDeltaE);               
   c.SetBranchAddress("candLepDeltaP", &candLepDeltaP);           
   c.SetBranchAddress("candLepDeltaPx", &candLepDeltaPx);         
   c.SetBranchAddress("candLepDeltaPy", &candLepDeltaPy);         
   c.SetBranchAddress("candLepDeltaPz", &candLepDeltaPz);         
   c.SetBranchAddress("candLepDeltaE", &candLepDeltaE);           
   c.SetBranchAddress("candNuDeltaP", &candNuDeltaP);             
   c.SetBranchAddress("candNuDeltaPx", &candNuDeltaPx);           
   c.SetBranchAddress("candNuDeltaPy", &candNuDeltaPy);           
   c.SetBranchAddress("candNuDeltaPz", &candNuDeltaPz);           
   c.SetBranchAddress("candNuDeltaE", &candNuDeltaE);             
   c.SetBranchAddress("candBDPx", &candBDPx);                     
   c.SetBranchAddress("candBDPy", &candBDPy);                     
   c.SetBranchAddress("candBDPz", &candBDPz);                     
   c.SetBranchAddress("candBDE", &candBDE);                       
   c.SetBranchAddress("candMissPx", &candMissPx);                 
   c.SetBranchAddress("candMissPy", &candMissPy);                 
   c.SetBranchAddress("candMissPz", &candMissPz);                 
   c.SetBranchAddress("candMissE", &candMissE);                   
   c.SetBranchAddress("candIntPur", &candIntPur);                 
   c.SetBranchAddress("candTagLeptons", &candTagLeptons);         
   c.SetBranchAddress("candTagChargedMult", &candTagChargedMult); 
   c.SetBranchAddress("candTagNeutMult", &candTagNeutMult);       
   c.SetBranchAddress("emcfwd", &emcfwd);                         
   c.SetBranchAddress("emcbwd", &emcbwd);                         
   c.SetBranchAddress("trkfwd", &trkfwd);                         
   c.SetBranchAddress("trkbwd", &trkbwd);                         
   c.SetBranchAddress("candM2", &candM2);                         
   c.SetBranchAddress("candQ2", &candQ2);                         
   c.SetBranchAddress("candPMiss", &candPMiss);                   
   c.SetBranchAddress("candEMiss", &candEMiss);                   
   c.SetBranchAddress("candThetaMiss", &candThetaMiss);           
   c.SetBranchAddress("candPhiMiss", &candPhiMiss);               
   c.SetBranchAddress("candThetaDL", &candThetaDL);               
   c.SetBranchAddress("candHelicity", &candHelicity);             
   c.SetBranchAddress("candM2Tru", &candM2Tru);                   
   c.SetBranchAddress("candM2TruBTag", &candM2TruBTag);           
   c.SetBranchAddress("candM2TruD", &candM2TruD);                 
   c.SetBranchAddress("candM2TruLep", &candM2TruLep);             
   c.SetBranchAddress("candMvaBestB", &candMvaBestB);             
   c.SetBranchAddress("candMvaComb", &candMvaComb);               
   c.SetBranchAddress("candMvaDl", &candMvaDl);                   
   c.SetBranchAddress("candCosT", &candCosT);                     
   c.SetBranchAddress("candMvaBestBL", candMvaBestBL);            
   c.SetBranchAddress("candBModeL", candBModeL);                  
   c.SetBranchAddress("candEExtraL", candEExtraL);                
   c.SetBranchAddress("candBntChaL", candBntChaL);                
   c.SetBranchAddress("candBntNeuL", candBntNeuL);                
   c.SetBranchAddress("candDntChaL", candDntChaL);                
   c.SetBranchAddress("candDntNeuL", candDntNeuL);                
   c.SetBranchAddress("candTypeL", candTypeL);                    
   c.SetBranchAddress("candLepTruL", candLepTruL);                
   c.SetBranchAddress("candDeltaEL", candDeltaEL);                
   c.SetBranchAddress("npi0", &npi0);                             
   c.SetBranchAddress("q2pi0", q2pi0);                             
   c.SetBranchAddress("mpi0", mpi0);                             
   c.SetBranchAddress("ppi0", ppi0);                             
   c.SetBranchAddress("dmpi0", dmpi0);                            
   c.SetBranchAddress("mm2pi0", mm2pi0);                          
   c.SetBranchAddress("pmisspi0", pmisspi0);                      
   c.SetBranchAddress("eextrapi0", eextrapi0);                    
   c.SetBranchAddress("e1pi0", e1pi0);                            
   c.SetBranchAddress("e2pi0", e2pi0);                            
   c.SetBranchAddress("bestepi0", bestepi0);                      
   c.SetBranchAddress("bestmasspi0", bestmasspi0);                   
   c.SetBranchAddress("ndpi0", &ndpi0);                         
   c.SetBranchAddress("dmdpi0", dmdpi0);                        
   c.SetBranchAddress("mm2dpi0", mm2dpi0);                      
   c.SetBranchAddress("pmissdpi0", pmissdpi0);                  
   c.SetBranchAddress("eextradpi0", eextradpi0);                
   c.SetBranchAddress("p1dpi0", p1dpi0);                        
   c.SetBranchAddress("p2dpi0", p2dpi0);                        
   c.SetBranchAddress("m1dpi0", m1dpi0);                        
   c.SetBranchAddress("m2dpi0", m2dpi0);                        
   c.SetBranchAddress("bestdpi0", bestdpi0);                    
   c.SetBranchAddress("nGTLHard", &nGTLHard);                   
   c.SetBranchAddress("nGTLSoft", &nGTLSoft);                   
   c.SetBranchAddress("nOtherTracks", &nOtherTracks);           
   c.SetBranchAddress("trkWeight", &trkWeight);                 
   c.SetBranchAddress("candKSdxy", &candKSdxy);                 
   c.SetBranchAddress("candKSTheta", &candKSTheta);             
   c.SetBranchAddress("candKSPT", &candKSPT);                   
   c.SetBranchAddress("nCandK", &nCandK);                       
   c.SetBranchAddress("nCandPi", &nCandPi);                     
   c.SetBranchAddress("candPKaons", candPKaons);                
   c.SetBranchAddress("candPPions", candPPions);                
   c.SetBranchAddress("candPPi0", &candPPi0);                   
   c.SetBranchAddress("kaonWeight", &kaonWeight);               
   c.SetBranchAddress("pionWeight", &pionWeight);               
   c.SetBranchAddress("elecWeight", &elecWeight);               
   c.SetBranchAddress("muonWeight", &muonWeight);                 

   Float_t impi0, ie1pi0, ie2pi0, iq2pi0, idmpi0, ieextrapi0, ippi0, imm2pi0, ipmisspi0, candMvaDssComb, candMvaDssDl;
   TString method = "BDT"; 
   float MES;
   TString basename = folder; basename += tagRest;
   TString ftree = "AWG82/ntuples/small/"; ftree += basename;
   ftree += "_Run"; ftree += Run; ftree += ".root";
   TFile f(ftree,"RECREATE");
   f.cd();

   int ibestepi0=-99, weighted = 0;
   float candMvaComb2=-99, candMvaDl2 = -99., Dldiff=0., Combdiff=0., wBF=1.0, wComb=1.0, wFF=1.0;
   cout<<"Starting tree"<<endl;
   TTree *Mva = new TTree("ntp1","ntp1");
   Mva->Branch("npi0",&npi0,"npi0/I");
   Mva->Branch("bestepi0",&ibestepi0,"bestepi0/I");
   Mva->Branch("q2pi0",&iq2pi0,"q2pi0/F");
   Mva->Branch("mpi0",&impi0,"mpi0/F");
   Mva->Branch("e1pi0",&ie1pi0,"e1pi0/F");
   Mva->Branch("e2pi0",&ie2pi0,"e2pi0/F");
   Mva->Branch("dmpi0",&idmpi0,"dmpi0/F");
   Mva->Branch("eextrapi0",&ieextrapi0,"eextrapi0/F");
   Mva->Branch("ppi0",&ippi0,"ppi0/F");
   Mva->Branch("mm2pi0",&imm2pi0,"mm2pi0/F");
   Mva->Branch("pmisspi0",&ipmisspi0,"pmisspi0/F");
   Mva->Branch("candEExtra",&candEExtra,"candEExtra/F");
   Mva->Branch("candMES",&MES,"candMES/F");
   Mva->Branch("candDeltaE",&candDeltaE,"candDeltaE/F");
   Mva->Branch("candDmass",&candDmass,"candDmass/F");
   Mva->Branch("candDeltam",&candDeltam,"candDeltam/F");
   Mva->Branch("candBTagDeltam",&candBTagDeltam,"candBTagDeltam/F");
   Mva->Branch("candBTagDmass",&candBTagDmass,"candBTagDmass/F");
   Mva->Branch("candRejectedTracks",&candRejectedTracks,"candRejectedTracks/I");
   Mva->Branch("candRejectedPhotons",&candRejectedPhotons,"candRejectedPhotons/I");
   Mva->Branch("candTagChargedMult",&candTagChargedMult,"candTagChargedMult/I");
   Mva->Branch("candTagNeutMult",&candTagNeutMult,"candTagNeutMult/I");
   Mva->Branch("candPMiss",&candPMiss,"candPMiss/F");
   Mva->Branch("candQ2",&candQ2,"candQ2/F");
   Mva->Branch("candM2",&candM2,"candM2/F");
   Mva->Branch("candCosT",&candCosT,"candCosT/F");
   Mva->Branch("candPstarLep",&candPstarLep,"candPstarLep/F");
   Mva->Branch("candLepTru",&candLepTru,"candLepTru/I");
   Mva->Branch("candIsMu",&candIsMu,"candIsMu/I");
   Mva->Branch("candBMode",&candBMode,"candBMode/I");
   Mva->Branch("candThetaLep",&candThetaLep,"candThetaLep/F");
   Mva->Branch("candFitProb",&candFitProb,"candFitProb/F");
   Mva->Branch("candExtraTracks",&candExtraTracks,"candExtraTracks/I");
   Mva->Branch("candBTru",&candBTru,"candBTru/I");

   Mva->Branch("trueLepMother",&trueLepMother,"trueLepMother/I");
   Mva->Branch("ntrueLep",&ntrueLep,"ntrueLep/I");
   Mva->Branch("trueMotherLep1",&trueMotherLep1,"trueMotherLep1/I");
   Mva->Branch("trueMotherLep2",&trueMotherLep2,"trueMotherLep2/I");
   Mva->Branch("trueMotherLep3",&trueMotherLep3,"trueMotherLep3/I");
   Mva->Branch("trueLundLep1",&trueLundLep1,"trueLundLep1/I");
   Mva->Branch("trueLundLep2",&trueLundLep2,"trueLundLep2/I");
   Mva->Branch("trueLundLep3",&trueLundLep3,"trueLundLep3/I");

   Mva->Branch("isSP6",&isSP6,"isSP6/I");
   Mva->Branch("MCType",&MCType,"MCType/I");
   Mva->Branch("MCSubmode",&MCSubmode,"MCSubmode/I");
   Mva->Branch("MCTaumode",&MCTaumode,"MCTaumode/I");
   Mva->Branch("MCCombmode",&MCCombmode,"MCCombmode/I");
   Mva->Branch("MCDssmode",&MCDssmode,"MCDssmode/I");
   Mva->Branch("MCD",&MCD,"MCD/I");
   Mva->Branch("MCPions",&MCPions,"MCPions/I");
   Mva->Branch("MCComblong",&MCComblong,"MCComblong/I");
   Mva->Branch("MCComblongD",&MCComblongD,"MCComblongD/I");
   Mva->Branch("MCCombB",&MCCombB,"MCCombB/I");
   Mva->Branch("MCCombDs",&MCCombDs,"MCCombDs/I");
   Mva->Branch("isBzero",&isBzero,"isBzero/I");
   Mva->Branch("MCUnsim",&MCUnsim,"MCUnsim/I");
   Mva->Branch("trueDmass",&trueDmass,"trueDmass/F");
   Mva->Branch("truePLep",&truePLep,"truePLep/F");
   Mva->Branch("truePD",&truePD,"truePD/F");
   Mva->Branch("truePPi0",&truePPi0,"truePPi0/F");
   Mva->Branch("trueDssPPi0",&trueDssPPi0,"trueDssPPi0/F");
   Mva->Branch("trueCTL",&trueCTL,"trueCTL/F");
   Mva->Branch("trueCTV",&trueCTV,"trueCTV/F");
   Mva->Branch("trueChi",&trueChi,"trueChi/F");
   Mva->Branch("trueQ2",&trueQ2,"trueQ2/F");
   Mva->Branch("trueLepCharge",&trueLepCharge,"trueLepCharge/I");
   Mva->Branch("trueTagDs",&trueTagDs,"trueTagDs/I");
   Mva->Branch("trueTagD",&trueTagD,"trueTagD/I");
   Mva->Branch("MCDoubleSL",&MCDoubleSL,"MCDoubleSL/I");

   Mva->Branch("candBntCha",&candBntCha,"candBntCha/I");
   Mva->Branch("candBntNeu",&candBntNeu,"candBntNeu/I");
   Mva->Branch("candDntCha",&candDntCha,"candDntCha/I");
   Mva->Branch("candDntNeu",&candDntNeu,"candDntNeu/I");
   Mva->Branch("candBnCharged",&candBnCharged,"candBnCharged/I");
   Mva->Branch("candBnNeutral",&candBnNeutral,"candBnNeutral/I");
   Mva->Branch("candDnCharged",&candDnCharged,"candDnCharged/I");
   Mva->Branch("candDnNeutral",&candDnNeutral,"candDnNeutral/I");
   Mva->Branch("candType",&candType,"candType/I");
   Mva->Branch("inCocktail",&inCocktail,"inCocktail/I");
   Mva->Branch("candDstarType",&candDstarType,"candDstarType/I");
   Mva->Branch("candDType",&candDType,"candDType/I");
   Mva->Branch("candMvaComb",&candMvaComb2,"candMvaComb/F");
   Mva->Branch("candMvaDl",&candMvaDl2,"candMvaDl/F");
   Mva->Branch("candMvaDssComb",&candMvaDssComb,"candMvaDssComb/F");
   Mva->Branch("candMvaDssDl",&candMvaDssDl,"candMvaDssDl/F");
   Mva->Branch("Dldiff",&Dldiff,"Dldiff/F");
   Mva->Branch("Combdiff",&Combdiff,"Combdiff/F");
   Mva->Branch("wBF",&wBF,"wBF/F");
   Mva->Branch("wComb",&wComb,"wComb/F");
   Mva->Branch("wFF",&wFF,"wFF/F");
   Mva->Branch("weighted",&weighted,"weighted/I");


  TMVA::Reader * readCombD0;
  TMVA::Reader * readCombDs0;
  TMVA::Reader * readCombDp;
  TMVA::Reader * readCombDsp;
  TMVA::Reader * readDlD0;
  TMVA::Reader * readDlDs0;
  TMVA::Reader * readDlDp; 
  TMVA::Reader * readDlDsp;

  TMVA::Reader * readDssCombD0;
  TMVA::Reader * readDssCombDs0;
  TMVA::Reader * readDssCombDp;
  TMVA::Reader * readDssCombDsp;
  TMVA::Reader * readDssDlD0;
  TMVA::Reader * readDssDlDs0;
  TMVA::Reader * readDssDlDp; 
  TMVA::Reader * readDssDlDsp;

   readCombD0 = new TMVA::Reader("V:Silent");
   readCombDs0 = new TMVA::Reader("!V:Silent");
   readCombDp = new TMVA::Reader("!V:Silent");
   readCombDsp = new TMVA::Reader("!V:Silent");
   readDlD0 = new TMVA::Reader("V:Silent");
   readDlDs0 = new TMVA::Reader("V:Silent");
   readDlDp = new TMVA::Reader("V:Silent");
   readDlDsp = new TMVA::Reader("V:Silent");



  readCombD0->AddVariable("candEExtra", &candEExtra);
  readCombD0->AddVariable("candMES", &MES);
  readCombD0->AddVariable("candDmass", &candDmass);
  readCombD0->AddVariable("candTagChargedMult", &candTagChargedMult);
  readCombD0->AddVariable("candBTagDeltam", &candBTagDeltam);
  readCombD0->AddVariable("candBTagDmass", &candBTagDmass);
  readCombD0->AddVariable("candDeltaE", &candDeltaE);
  readCombD0->AddVariable("candCosT", &candCosT);

  TString file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaCombD0_BDT.weights.txt";
  readCombD0->BookMVA(method,file);

  readDlD0->AddVariable("candEExtra", &candEExtra);
  readDlD0->AddVariable("candMES", &MES);
  readDlD0->AddVariable("candDmass", &candDmass);
  readDlD0->AddVariable("candTagChargedMult", &candTagChargedMult);
  readDlD0->AddVariable("candBTagDeltam", &candBTagDeltam);
  readDlD0->AddVariable("candBTagDmass", &candBTagDmass);
  readDlD0->AddVariable("candDeltaE", &candDeltaE);
  readDlD0->AddVariable("candCosT", &candCosT);

  file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDlD0_BDT.weights.txt";
  readDlD0->BookMVA(method,file);

  readCombDs0->AddVariable("candEExtra", &candEExtra);
  readCombDs0->AddVariable("candMES", &MES);
  readCombDs0->AddVariable("candDmass", &candDmass);
  readCombDs0->AddVariable("candDeltam", &candDeltam);
  readCombDs0->AddVariable("candTagChargedMult", &candTagChargedMult);
  readCombDs0->AddVariable("candBTagDeltam", &candBTagDeltam);
  readCombDs0->AddVariable("candBTagDmass", &candBTagDmass);
  readCombDs0->AddVariable("candDeltaE", &candDeltaE);
  readCombDs0->AddVariable("candCosT", &candCosT);

  file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaCombDs0_BDT.weights.txt";
  readCombDs0->BookMVA(method,file);

  readDlDs0->AddVariable("candEExtra", &candEExtra);
  readDlDs0->AddVariable("candMES", &MES);
  readDlDs0->AddVariable("candDmass", &candDmass);
  readDlDs0->AddVariable("candDeltam", &candDeltam);
  readDlDs0->AddVariable("candTagChargedMult", &candTagChargedMult);
  readDlDs0->AddVariable("candBTagDeltam", &candBTagDeltam);
  readDlDs0->AddVariable("candBTagDmass", &candBTagDmass);
  readDlDs0->AddVariable("candDeltaE", &candDeltaE);
  readDlDs0->AddVariable("candCosT", &candCosT);

  file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDlDs0_BDT.weights.txt";
  readDlDs0->BookMVA(method,file);

  readCombDp->AddVariable("candEExtra", &candEExtra);
  readCombDp->AddVariable("candMES", &MES);
  readCombDp->AddVariable("candDmass", &candDmass);
  readCombDp->AddVariable("candTagChargedMult", &candTagChargedMult);
  readCombDp->AddVariable("candBTagDeltam", &candBTagDeltam);
  readCombDp->AddVariable("candBTagDmass", &candBTagDmass);
  readCombDp->AddVariable("candDeltaE", &candDeltaE);
  readCombDp->AddVariable("candCosT", &candCosT);

  file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaCombDp_BDT.weights.txt";
  readCombDp->BookMVA(method,file);

  readDlDp->AddVariable("candEExtra", &candEExtra);
  readDlDp->AddVariable("candMES", &MES);
  readDlDp->AddVariable("candDmass", &candDmass);
  readDlDp->AddVariable("candTagChargedMult", &candTagChargedMult);
  readDlDp->AddVariable("candBTagDeltam", &candBTagDeltam);
  readDlDp->AddVariable("candBTagDmass", &candBTagDmass);
  readDlDp->AddVariable("candDeltaE", &candDeltaE);
  readDlDp->AddVariable("candCosT", &candCosT);

  file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDlDp_BDT.weights.txt";
  readDlDp->BookMVA(method,file);

  readCombDsp->AddVariable("candEExtra", &candEExtra);
  readCombDsp->AddVariable("candMES", &MES);
  readCombDsp->AddVariable("candDmass", &candDmass);
  readCombDsp->AddVariable("candDeltam", &candDeltam);
  readCombDsp->AddVariable("candTagChargedMult", &candTagChargedMult);
  readCombDsp->AddVariable("candBTagDeltam", &candBTagDeltam);
  readCombDsp->AddVariable("candBTagDmass", &candBTagDmass);
  readCombDsp->AddVariable("candDeltaE", &candDeltaE);
  readCombDsp->AddVariable("candCosT", &candCosT);

  file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaCombDsp_BDT.weights.txt";
  readCombDsp->BookMVA(method,file);

  readDlDsp->AddVariable("candEExtra", &candEExtra);
  readDlDsp->AddVariable("candMES", &MES);
  readDlDsp->AddVariable("candDmass", &candDmass);
  readDlDsp->AddVariable("candDeltam", &candDeltam);
  readDlDsp->AddVariable("candTagChargedMult", &candTagChargedMult);
  readDlDsp->AddVariable("candBTagDeltam", &candBTagDeltam);
  readDlDsp->AddVariable("candBTagDmass", &candBTagDmass);
  readDlDsp->AddVariable("candDeltaE", &candDeltaE);
  readDlDsp->AddVariable("candCosT", &candCosT);

  file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDlDsp_BDT.weights.txt";
  readDlDsp->BookMVA(method,file);

  //===================== D** MVAs ===============================================
   readDssCombD0 = new TMVA::Reader("V:Silent");
   readDssDlD0 = new TMVA::Reader("V:Silent");
   readDssCombDs0 = new TMVA::Reader("V:Silent");
   readDssDlDs0 = new TMVA::Reader("V:Silent");
   readDssCombDp = new TMVA::Reader("V:Silent");
   readDssDlDp = new TMVA::Reader("V:Silent");
   readDssCombDsp = new TMVA::Reader("V:Silent");
   readDssDlDsp = new TMVA::Reader("V:Silent");

   readDssCombD0->AddVariable("mpi0", &impi0);
   readDssCombD0->AddVariable("candDmass", &candDmass);
   readDssCombD0->AddVariable("dmpi0", &idmpi0);
   readDssCombD0->AddVariable("eextrapi0", &ieextrapi0);
   readDssCombD0->AddVariable("ppi0", &ippi0);
   readDssCombD0->AddVariable("e1pi0", &ie1pi0);
   readDssCombD0->AddVariable("candCosT", &candCosT);
   readDssCombD0->AddVariable("candMES", &MES);
   readDssCombD0->AddVariable("candDeltaE", &candDeltaE);
   file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDssCombD0_BDT.weights.txt";
   readDssCombD0->BookMVA(method,file);

   readDssDlD0->AddVariable("mpi0", &impi0);
   readDssDlD0->AddVariable("candDmass", &candDmass);
   readDssDlD0->AddVariable("dmpi0", &idmpi0);
   readDssDlD0->AddVariable("eextrapi0", &ieextrapi0);
   readDssDlD0->AddVariable("ppi0", &ippi0);
   readDssDlD0->AddVariable("e1pi0", &ie1pi0);
   readDssDlD0->AddVariable("candCosT", &candCosT);
   readDssDlD0->AddVariable("candMES", &MES);
   readDssDlD0->AddVariable("candDeltaE", &candDeltaE);
   file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDssDlD0_BDT.weights.txt";
   readDssDlD0->BookMVA(method,file);

   readDssCombDs0->AddVariable("mpi0", &impi0);
   readDssCombDs0->AddVariable("candDmass", &candDmass);
   readDssCombDs0->AddVariable("dmpi0", &idmpi0);
   readDssCombDs0->AddVariable("eextrapi0", &ieextrapi0);
   readDssCombDs0->AddVariable("ppi0", &ippi0);
   readDssCombDs0->AddVariable("e1pi0", &ie1pi0);
   readDssCombDs0->AddVariable("candCosT", &candCosT);
   readDssCombDs0->AddVariable("candDeltam", &candDmass);
   readDssCombDs0->AddVariable("candMES", &MES);
   readDssCombDs0->AddVariable("candDeltaE", &candDeltaE);
   file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDssCombDs0_BDT.weights.txt";
   readDssCombDs0->BookMVA(method,file);

   readDssDlDs0->AddVariable("mpi0", &impi0);
   readDssDlDs0->AddVariable("candDmass", &candDmass);
   readDssDlDs0->AddVariable("dmpi0", &idmpi0);
   readDssDlDs0->AddVariable("eextrapi0", &ieextrapi0);
   readDssDlDs0->AddVariable("ppi0", &ippi0);
   readDssDlDs0->AddVariable("e1pi0", &ie1pi0);
   readDssDlDs0->AddVariable("candCosT", &candCosT);
   readDssDlDs0->AddVariable("candDeltam", &candDmass);
   readDssDlDs0->AddVariable("candMES", &MES);
   readDssDlDs0->AddVariable("candDeltaE", &candDeltaE);
   file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDssDlDs0_BDT.weights.txt";
   readDssDlDs0->BookMVA(method,file);

   readDssCombDp->AddVariable("mpi0", &impi0);
   readDssCombDp->AddVariable("candDmass", &candDmass);
   readDssCombDp->AddVariable("dmpi0", &idmpi0);
   readDssCombDp->AddVariable("eextrapi0", &ieextrapi0);
   readDssCombDp->AddVariable("ppi0", &ippi0);
   readDssCombDp->AddVariable("e1pi0", &ie1pi0);
   readDssCombDp->AddVariable("candCosT", &candCosT);
   readDssCombDp->AddVariable("candMES", &MES);
   readDssCombDp->AddVariable("candDeltaE", &candDeltaE);
   file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDssCombDp_BDT.weights.txt";
   readDssCombDp->BookMVA(method,file);

   readDssDlDp->AddVariable("mpi0", &impi0);
   readDssDlDp->AddVariable("candDmass", &candDmass);
   readDssDlDp->AddVariable("dmpi0", &idmpi0);
   readDssDlDp->AddVariable("eextrapi0", &ieextrapi0);
   readDssDlDp->AddVariable("ppi0", &ippi0);
   readDssDlDp->AddVariable("e1pi0", &ie1pi0);
   readDssDlDp->AddVariable("candCosT", &candCosT);
   readDssDlDp->AddVariable("candMES", &MES);
   readDssDlDp->AddVariable("candDeltaE", &candDeltaE);
   file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDssDlDp_BDT.weights.txt";
   readDssDlDp->BookMVA(method,file);

   readDssCombDsp->AddVariable("mpi0", &impi0);
   readDssCombDsp->AddVariable("candDmass", &candDmass);
   readDssCombDsp->AddVariable("dmpi0", &idmpi0);
   readDssCombDsp->AddVariable("eextrapi0", &ieextrapi0);
   readDssCombDsp->AddVariable("ppi0", &ippi0);
   readDssCombDsp->AddVariable("e1pi0", &ie1pi0);
   readDssCombDsp->AddVariable("candCosT", &candCosT);
   readDssCombDsp->AddVariable("candDeltam", &candDmass);
   readDssCombDsp->AddVariable("candMES", &MES);
   readDssCombDsp->AddVariable("candDeltaE", &candDeltaE);
   file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDssCombDsp_BDT.weights.txt";
   readDssCombDsp->BookMVA(method,file);

   readDssDlDsp->AddVariable("mpi0", &impi0);
   readDssDlDsp->AddVariable("candDmass", &candDmass);
   readDssDlDsp->AddVariable("dmpi0", &idmpi0);
   readDssDlDsp->AddVariable("eextrapi0", &ieextrapi0);
   readDssDlDsp->AddVariable("ppi0", &ippi0);
   readDssDlDsp->AddVariable("e1pi0", &ie1pi0);
   readDssDlDsp->AddVariable("candCosT", &candCosT);
   readDssDlDsp->AddVariable("candDeltam", &candDmass);
   readDssDlDsp->AddVariable("candMES", &MES);
   readDssDlDsp->AddVariable("candDeltaE", &candDeltaE);
   file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDssDlDsp_BDT.weights.txt";
   readDssDlDsp->BookMVA(method,file);

   int IsCocktail = 0;
   if(folder.Contains("cocktail")) IsCocktail = 1;
   WeightManager myWM("babar_code/Reweight/wFFBF.txt",IsCocktail);
   WeightManager myWMFF(weightName,IsCocktail);
   Long64_t nentries = c.GetEntries();
   //nentries = 100;

   for (Long64_t jentry=0; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry, &c);
      if (ientry < 0) break;
      c.GetEntry(jentry);
      if (jentry%10000 == 0) cout<<"Entry "<<jentry<<" of "<<nentries<<endl;
      
      if(mEScorrection){
	MES = PidCorrectMesMean::correctMes(candMES,runnum);	
      } else MES = candMES;
      if (npi0==0){
	iq2pi0 = impi0 = ie1pi0 = ie2pi0 = idmpi0 = ieextrapi0 = ippi0 = imm2pi0 = ipmisspi0 = -99.;
	candMvaDssComb = candMvaDssDl = -99.;
	ibestepi0 = 0;
      } else {
	int i=0; ibestepi0 = 1;
	while(i<npi0) {
	  if(bestepi0[i] == 1) break;
	  i++;
	}
	//if(bestepi0[i] != 1 || npi0==0) continue;
	iq2pi0 = q2pi0[i];
	impi0 = mpi0[i];
	ie1pi0 = e1pi0[i];
	ie2pi0 = e2pi0[i];
	idmpi0 = dmpi0[i];
	ieextrapi0 = eextrapi0[i];
	ippi0 = ppi0[i];
	imm2pi0 = mm2pi0[i];
	ipmisspi0 = pmisspi0[i];
	if(candType==1){
	  candMvaDssComb = readDssCombD0->EvaluateMVA(method);
	  candMvaDssDl = readDssDlD0->EvaluateMVA(method);
	}
	if(candType==2){
	  candMvaDssComb = readDssCombDs0->EvaluateMVA(method);
	  candMvaDssDl = readDssDlDs0->EvaluateMVA(method);
	}
	if(candType==3){
	  candMvaDssComb = readDssCombDp->EvaluateMVA(method);
	  candMvaDssDl = readDssDlDp->EvaluateMVA(method);
	}
	if(candType==4){
	  candMvaDssComb = readDssCombDsp->EvaluateMVA(method);
	  candMvaDssDl = readDssDlDsp->EvaluateMVA(method);
	}
      }

      if(candType==1){
	candMvaComb2 = readCombD0->EvaluateMVA(method);
	candMvaDl2 = readDlD0->EvaluateMVA(method);
      }
      if(candType==2){
	candMvaComb2 = readCombDs0->EvaluateMVA(method);
	candMvaDl2 = readDlDs0->EvaluateMVA(method);
     }
      if(candType==3){
	candMvaComb2 = readCombDp->EvaluateMVA(method);
	candMvaDl2 = readDlDp->EvaluateMVA(method);
      }
      if(candType==4){
	candMvaComb2 = readCombDsp->EvaluateMVA(method);
	candMvaDl2 = readDlDsp->EvaluateMVA(method);
      }
      Dldiff = candMvaDl - candMvaDl2;
      Combdiff = candMvaComb - candMvaComb2;
      int Type = candType;
      if(folder.Contains("uds")) MCType = -2;
      if(folder.Contains("ccbar")) MCType = -1;
      if(imm2pi0>-4&&imm2pi0<12&&candPstarLep>0&&candPstarLep<2.4&&candMES>5.2&&candMES<5.3&&ipmisspi0>.2&&
	 ibestepi0==1&&fabs(candCosT)<.8&&(candMvaDssComb>0.13&&candMvaDssDl>0.15&&candType==1)||
	 (candMvaDssComb>0.13&&candMvaDssDl>0.15&&candType==2)||
	 (candMvaDssComb>0.15&&candMvaDssDl>0.11&&candType==3)||
	 (candMvaDssComb>0.1&&candMvaDssDl>0.17&&candType==4)) Type += 4;
      wBF = myWM.getEventWeight(Type,candDstarType,MCType,MCSubmode,MCDssmode,MCD,MCPions,
			       isBzero,isSP6,MCTaumode,candPstarLep,trueDmass,trueCTL,trueCTV,trueChi,trueQ2,
			       trueLepCharge,candM2);
      wFF = myWMFF.getEventWeight(candType,candDstarType,MCType,MCSubmode,MCDssmode,MCD,MCPions,
			       isBzero,isSP6,MCTaumode,candPstarLep,trueDmass,trueCTL,trueCTV,trueChi,trueQ2,
			       trueLepCharge,candM2);
      wComb = myWM.getCombWeight(MCCombB,MCCombDs,MCDoubleSL,candLepTru);
      Dmode = (candBMode/100)%100; Xmode = candBMode%100;
      inCocktail = modeMap[Dmode][Xmode];
      if(wComb<0){
	weighted = 0;
	wComb = fabs(wComb);
      } else weighted = 1;

      Mva->Fill();
   }
   Mva->Write();
   cout<<"Writing the file "<<ftree<<endl;
   f.Close();

   return 1;


}


Long64_t LoadTree(Long64_t entry, TChain* fChain)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   return centry;
}
