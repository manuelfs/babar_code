// $Id: weightManager.cc,v 1.10 2012/08/23 02:22:18 manuelf Exp $
//---------------------------------------------------------------------------------
// Description:
//      Macro that takes care of the re-weighting
//
// Author List:
//      Michael Mazur                             UC Santa Barbara
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      09/07/09 manuelf -- Adaptation from Mazur's weightManager.cc
//---------------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "DonutUtils/weightManager.hh"
#include "DonutUtils/BToDstaunu.hh"
#include "DonutUtils/BToDtaunu.hh"

#include "TROOT.h"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooStringVar.hh"
#include "XslFFReweighting/XSLBToDSSlnu_LLSW.hh"
#include "XslFFReweighting/XSLEvtFFWeight.hh"

#include <fstream>
#include <iostream>
using std::cout;
using std::endl;

Float_t WeightManager::getVariWeight(TString Vari, Float_t value, Int_t candType){
  if(!doVariWeights) return 1.0;
  for(int vari=0; vari<nVari; vari++){
     if(Variables[vari] == Vari){
       if(value < ranges[vari][0] || value > ranges[vari][nbins[vari]]) return 1.0;
       for(int bin=0; bin<nbins[vari]; bin++){
	 if(value >= ranges[vari][bin] && value < ranges[vari][bin+1]) return variWeights[vari][candType-1][bin];
       }
     }
  }
  return 1.0;
}

Float_t WeightManager::getBDTBias(){
  if(!doVariWeights) return 1.0;
  return 1/BDTBias;
}

Float_t WeightManager::getEventWeight(Int_t candType, Int_t candDstarType, Int_t MCType,Int_t MCSubmode,
				      Int_t MCDssmode,Int_t MCD,Int_t MCPions, Int_t isBzero, Int_t isSP6,
				      Int_t MCTaumode, Float_t candPstarLep,Float_t trueDmass,
				      Float_t ctl,Float_t ctv,Float_t chi,Float_t q2,Int_t lcharge,
				      Float_t candM2)
{
  Float_t w = 1.0, pi=3.14159265;
  if(MCType>=0){
    if (doBkgWeights)  w /= wBkg[candType-1]; //Btag efficiency taken from the combinatorial
    if (doFFWeights && MCType==14)  w *= DssFFWeight(MCD,ctl,q2,trueDmass);
    if (doFFWeights)    w *= FFWeight(MCType,MCSubmode,-ctl,ctv,chi,q2,lcharge);
    if (doBFWeights)    w *= BFWeight(candType, MCType,MCSubmode,MCD,MCPions,MCTaumode);
    if (doDlnFFWeights) w *= DlnFFWeight(MCType,pi-ctl,q2);
    if (doSigMC)        w *= SigMCWeight(candType,MCType,MCSubmode,MCDssmode,isBzero,isSP6,MCTaumode);
  } else {
    if (doPlWeights)    w *= PlWeight(candPstarLep, candType);
    if (MCType==-1)     w *= ccbarTouds;
    if (doBkgWeights)   w *= BkgWeight(MCType, candType);
  }
  if(isZero) w = 0;
  
  return w;
}

Float_t WeightManager::BkgWeight(Int_t MCType, Int_t candType){
  double w = 1.;
  if(candType>=1 && candType<=4){
    if(MCType==0) w /= wBkg[candType-1];
    //w /= wBkgmES[candType-1];  // For now, the mES sideband bias is only corrected in MultiEventMixer
    return w;
  }
  cout<<"BkgWeight: Bad candType: "<<candType<<endl; return 1.;
}

Float_t WeightManager::PlWeight(Float_t candPstarLep, Int_t candType){
  int bin = (int)(candPstarLep/2.4*1000);
  if(bin<0 || bin>=1000) {cout<<"Bad p*l: "<<candPstarLep<<endl; return 1.;}
  if(candType>=1 && candType<=4) return wPlD[candType-1][bin];
  cout<<"PlWeight: Bad candType: "<<candType<<endl; return 1.;
}

Float_t WeightManager::BFWeight(Int_t candType, Int_t MCType,Int_t MCSubmode,Int_t MCD,Int_t MCPions,Int_t MCTaumode){
  Float_t wBF = 1.0;
  Float_t RD0 = 0.313, RDc = 0.338, RDs0 = 0.259, RDsc = 0.281, RDss0 = 0.251, RDssc = 0.272;
  Float_t RD = 0.30, RDs = 0.25, RDss = 0.18; //R(Dtau/Dl), R(D*tau/D*l) and R(D**tau/D**l)
  //Float_t RD = 0.449, RDs = 0.323, RDss = 0.18; //R(Dtau/Dl), R(D*tau/D*l) and R(D**tau/D**l)

  Int_t lund = abs(MCD);
  if (MCType == 13) {
    if (lund == 411 && MCPions == 1)       wBF = wDpi;
    else if (lund == 421 && MCPions == 10) wBF = wD0pi0;
    else if (lund == 413 && MCPions == 1)  wBF = wDstarpi;
    else if (lund == 423 && MCPions == 10) wBF = wDstar0pi0;
    else if (lund == 411 && MCPions == 10) wBF = wDpi0;
    else if (lund == 421 && MCPions == 1)  wBF = wD0pi;
    else if (lund == 413 && MCPions == 10) wBF = wDstarpi0;
    else if (lund == 423 && MCPions == 1)  wBF = wDstar0pi;
    else if (MCPions == 2 || MCPions == 11 || MCPions == 20)  wBF = 1.;
    else cout << "cannot reweight nonresonant D** mode with D = " << MCD << " and pi = " << MCPions <<
	   endl << "taking weight = 1.0" << endl;
  }
  if (MCType == 14) {
    switch (lund) {
    case 10423: wBF = wD10; break;
    case 425  : wBF = wD20; break;
    case 10421: wBF = wD00; break;
    case 20423: wBF = wD1prime0; break;
    case 10413: wBF = wD1p; break;
    case 415  : wBF = wD2p; break;
    case 10411: wBF = wD0p; break;
    case 20413: wBF = wD1primep; break;
    default: cout << "cannot reweight resonant D** mode " << MCD << endl << "taking weight = 1.0" << endl;
    }
    if(MCTaumode > -1){
      if(lund%100 > 20) wBF *= RDss/RDss0;
      else              wBF *= RDss/RDssc;
    }
  }
  // In R18 Generics, there's 40 NR, 203 D**l and 66 D**tau in B0
  // and 38 NR, 225 D**l and 66 D**tau in B+
  if(isCocktail==1 && MCType>=13&&MCType<=14){
    double tot;
    if (MCType == 13) {
      int charge = (lund%100)+MCPions;
      if(charge>15&&charge<25) tot = 270/243.; // charge of D(*)+pi0 is 21(23)
      else tot = 270/263.;
      if (lund == 411 && MCPions == 1)       wBF *= 19./60*tot;
      else if (lund == 421 && MCPions == 10) wBF *= 10./30*tot;
      else if (lund == 413 && MCPions == 1)  wBF *=  6./20*tot;
      else if (lund == 423 && MCPions == 10) wBF *=  3./10*tot;
      else if (lund == 411 && MCPions == 10) wBF *= 10./30*tot;
      else if (lund == 421 && MCPions == 1)  wBF *= 20./60*tot;
      else if (lund == 413 && MCPions == 10) wBF *=  3./10*tot;
      else if (lund == 423 && MCPions == 1)  wBF *=  7./20*tot;
    }
    if(MCTaumode>-1) {
      tot = 150/66.;
      switch (lund) {
      case 10423: wBF *= 13./56*tot; break;
      case 425  : wBF *= 20./37*tot; break;
      case 10421: wBF *= 13./20*tot; break;
      case 20423: wBF *= 20./37*tot; break;
      case 10413: wBF *= 13./56*tot; break;
      case 415  : wBF *= 20./37*tot; break;
      case 10411: wBF *= 13./20*tot; break;
      case 20413: wBF *= 20./37*tot; break;
      }
    } else { 
      if((lund%100)>15) tot = 270/263.;
      else tot = 270/243.;
      switch (lund) {
      case 10423: wBF *= 56./56*tot;  break;
      case 425  : wBF *= 30./37*tot;  break;
      case 10421: wBF *= 49./20*tot;  break;
      case 20423: wBF *= 90./37*tot;  break;
      case 10413: wBF *= 52./56*tot;  break;
      case 415  : wBF *= 23./37*tot;  break;
      case 10411: wBF *= 45./20*tot;  break;
      case 20413: wBF *= 83./37*tot;  break;
      }
    }
  }
//   if (MCType == 1 || MCType == 3 || MCType == 5) wBF = wD0*BFIterRatio[0];
//   if (MCType == 2 || MCType == 4 || MCType == 6) wBF = wDs0*BFIterRatio[1];
//   if (MCType == 7 || MCType == 9 || MCType == 11) wBF = wDp*BFIterRatio[2];
//   if (MCType == 8 || MCType == 10 || MCType == 12) wBF = wDsp*BFIterRatio[3];
//   if (MCType == 5)   wBF *= RD/RD0*BFIterRatio[4];
//   if (MCType == 6)   wBF *= RDs/RDs0*BFIterRatio[5];
//   if (MCType == 11)  wBF *= RD/RDc*BFIterRatio[4];
//   if (MCType == 12)  wBF *= RDs/RDsc*BFIterRatio[5];
//   if ((MCType == 2 || MCType == 4 || MCType == 6)&&candType==1) wBF *= BFIterRatio[6];
//   if ((MCType == 8 || MCType == 10 || MCType == 12)&&candType==3) wBF *= BFIterRatio[7];

  if (MCType == 1 || MCType == 3 || MCType == 5) wBF = wD0;
  if (MCType == 2 || MCType == 4 || MCType == 6) wBF = wDs0;
  if (MCType == 7 || MCType == 9 || MCType == 11) wBF = wDp;
  if (MCType == 8 || MCType == 10 || MCType == 12) wBF = wDsp;
  if (MCType == 5)   wBF *= RD/RD0;
  if (MCType == 6)   wBF *= RDs/RDs0;
  if (MCType == 11)  wBF *= RD/RDc;
  if (MCType == 12)  wBF *= RDs/RDsc;
  if (MCPions == 2 || MCPions == 11 || MCPions == 20 || MCPions == 100 || MCPions==0&&MCType==14)  wBF = 1.;
  return wBF;
}

Float_t WeightManager::DssFFWeight(Int_t MCD, Float_t ctl, Float_t q2, Float_t mDss){
  Float_t w = 1.0;
  int dss;
  MCD = abs(MCD);
  if(MCD==20423      || MCD==20413){ dss = -1; }    // D1'
  else if(MCD==10421 || MCD==10411){ dss = 0;  }     // D0
  else if(MCD==10423 || MCD==10413){ dss = 1;  }     // D1
  else if(MCD==425   || MCD==415){ dss = 2;;   }      // D2
  else return w;

  XSLBToDSSlnu_LLSW ffw(mB,mDss,q2,ctl,dss,-1.5,1);
  w *= ffw.FromSP8ToThisModel();
  return w;
}

Float_t WeightManager::FFWeight(Int_t MCType,Int_t MCSubmode,Float_t ctl,Float_t ctv,Float_t chi,
				Float_t q2,Int_t lcharge) {
  double w = 1;
  int isDgamma = 0; if(MCSubmode>=200&&MCType<7) isDgamma=1;
  BToDstaunu ffw(HQET_rho2, HQET_R1, HQET_R2, HQET_R0, 0);
  if (MCType==8 || MCType==10 || MCType==12) ffw.SetMasses(0);
  if (MCType == 6 || MCType == 12) {
    ffw._gSR = gSR;
    w = ffw.FromSP8ToThisModel(q2,ctl,ctv,chi,isDgamma,lcharge==1,ffw.mTau);
  }else if(MCType == 4 || MCType == 10) w = ffw.FromSP8ToThisModel(q2,ctl,ctv,chi,isDgamma,lcharge==1,ffw.mMu);
  else if(MCType == 2 || MCType == 8)  w = ffw.FromSP8ToThisModel(q2,ctl,ctv,chi,isDgamma,lcharge==1,ffw.mE);

  if (isnan(w)) return 0;
  return w;
}

Float_t WeightManager::DlnFFWeight(Int_t MCType,Float_t ctl,Float_t q2){
  double w = 1;
  BToDtaunu ffw(rhoD2, Delta, 0);
  if (MCType==7 || MCType==9 || MCType==11) ffw.SetMasses(0);
  if (MCType==1 || MCType==7)          w = ffw.FromSP8ToThisModel(q2, ctl, ffw.mE);
  else if(MCType == 3 || MCType == 9)  w = ffw.FromSP8ToThisModel(q2, ctl, ffw.mMu);
  else if(MCType == 5 || MCType == 11){
    ffw._gSR = gSR;
    w = ffw.FromSP8ToThisModel(q2, ctl, ffw.mTau);
  }
  if (isnan(w)) return 0;
  return w;
}

Float_t WeightManager::SigMCWeight(Int_t candType, Int_t MCType,Int_t MCSubmode,Int_t MCDssmode,
				   Int_t isBzero,Int_t isSP6,Int_t MCTaumode){
  Float_t RD0 = 0.313, RDc = 0.338, RDs0 = 0.259, RDsc = 0.281, RDss0 = 0.251, RDssc = 0.272; //SP8 values
  Float_t w = 1.0, TauToL = (0.1778+0.1731)/2; // B(tau->lnunu)=0.176 in SP8

  if (sigMode == 1 && isCocktail) {
    if(MCType==14 && MCTaumode>-1) {  //Mixing 3.4M of D**taunu and 10.1M of D**lnu 
      if(isBzero)       w = RDssc*TauToL*10131/3354; 
      else              w = RDss0*TauToL*10137/3354;
    }
    if(candType>4) {
      if(MCType==5)       w = RD0*TauToL*4.474/4.492;
      if(MCType==11)      w = RDc*TauToL*4.492/4.486;
      if(MCType==6)       w = RDs0*TauToL*4.477/4.488; 
      if(MCType==12){
	if(MCSubmode<200) w = RDsc*TauToL*4.475/4.492*0.683; // B(D*+ -> D0pi) = 68.3% in SP8
	else              w = RDsc*TauToL*4.475/4.486*0.306; // B(D*+ -> Dpi0) = 30.6% in SP8
      }
    }
  }
  return w;
}

Float_t WeightManager::getCombWeight(Int_t MCCombB,Int_t MCCombDs,Int_t MCDoubleSL,Int_t candType){
  if (!doBFWeights) return 1.;
  Float_t wcomb = 1.0;
  if(doBkgWeights) wcomb *= BkgWeight(0, candType);
  if (MCDoubleSL==1) return wDoubleSL*wcomb;
  //if (MCCombB == 0 && candTruLep == 0) return wFakeLep*wcomb;
  if (MCCombB == 0) return wUnknown*wcomb;


  if (MCCombDs) switch (MCCombDs) {
  case 1: wcomb *= wDsTaunu; break;
  case 2: wcomb *= wDsmunu; break;
  case 3: wcomb *= wUnknown; break;
  case 4: wcomb *= wDsphilnu; break;
  case 5: wcomb *= wDsEtalnu; break;
  case 6: wcomb *= wDsEtaprimelnu; break;
  }

  switch (MCDoubleSL) {
  case 10: wcomb *= wDzba1p; break;
  case 20: wcomb *= wDma1p; break;
  case 30: wcomb *= wDstarzba1p; break;
  case 40: wcomb *= wDstarma1p; break;
  }
  switch (MCCombB) {
  case 71: wcomb *= wDzbDs; break;
  case 81: wcomb *= wDzbDsstar; break;
  case 73: wcomb *= wDstarzbDs; break;
  case 83: wcomb *= wDstarzbDsstar; break;
  case 75: 
  case 85: wcomb *= wDstarstarDs; break;
  case 421: wcomb *= wDpDzbKz; break;
  case 441: wcomb *= wDzbDstarpKz; break;
  case 443: wcomb *= wDstarzbDstarpKz; break;
  case 432: wcomb *= wDstarzbDpKz; break;
  case 311: wcomb *= wDzDzbKp; break;
  case 333: wcomb *= wDstarzbDstarzKp; break;
  case 331: wcomb *= wDstarzDzbKp; break;
  case 342: wcomb *= wDstarmDpKp; break;
  case 210: wcomb *= wDzbRhop; break;
  case 230: wcomb *= wDstarzbRhop; break;

  case 72: wcomb *= wDmDs; break;
  case 74: wcomb *= wDstarmDs; break;
  case 82: wcomb *= wDmDsstar; break;
  case 84: wcomb *= wDstarmDsstar; break;
  case 321: wcomb *= wDmDzKp; break;
  case 332: wcomb *= wDmDstarzKp; break;
  case 343: wcomb *= wDstarmDstarzKp; break;
  case 422: wcomb *= wDpDmKz; break;
  case 442: wcomb *= wDstarmDpKz; break;
  case 444: wcomb *= wDstarmDstarpKz; break;
  case 431: wcomb *= wDstarzDzbKz; break;
  case 44: wcomb *= wDstarpDstarm; break;
  case 42: wcomb *= wDpDstarm; break;
  case 22: wcomb *= wDpDm; break;
  case 220: wcomb *= wDmRhop; break;
  case 240: wcomb *= wDstarmRhop; break;
  default: wcomb *= wUnknown;
  }
  return wcomb;
}

WeightManager::WeightManager(const char *weightFile, int IsCocktail, Bool_t verb){
  doVariWeights  = kFALSE;
  doBkgWeights   = kFALSE;
  doPlWeights    = kFALSE;
  doBFWeights    = kFALSE;
  doFFWeights    = kFALSE;
  doDlnFFWeights = kFALSE;
  doPi0Weights   = kFALSE;
  doSigMC        = kFALSE;
  verbose        = verb;
  if(IsCocktail > 0)  if(verbose) cout<<"Doing Cocktail. Using SP4 FF parameters and BF"<<endl;
  isCocktail = IsCocktail;
  isZero = false;

  sigMode = 0;

  for(int cand=0; cand<4; cand++) wBkg[cand]=1;
  BDTBias = 1;
  mD0 = 1.865; mDp = 1.869; mDs0 = 2.007; mDsp = 2.010;
  mB = 5.279;
  double totuds = 2.190254e+09, totccbar = 1.12736e+09;
  ccbarTouds = totuds/totccbar/2.09*1.3; 
  if (weightFile == 0) {
    if(verbose) cout << "running without reweighting" << endl;
  } else if (*weightFile == '0') {
      isZero = true;
  } else if (weightFile) {
    if(verbose) cout << "reading weight file " << weightFile << endl;
    RooStringVar variweights("variWeights","variWeights","");// Re-weighting for BDT
    RooStringVar bkgweights ("BkgWeights","BkgWeights","");  // Normalization correction for comb+cont
    RooStringVar plweightsD ("PlWeightsD","PlWeightsD","");  // Lepton momentum re-weighting for continuum
    RooStringVar bfweights  ("BFWeights" ,"BFWeights","");   // Branching Fractions
    RooStringVar ffweights  ("FFWeights" ,"FFWeights","");   // FF of D*lnu: r1, r2, rho^2
    RooStringVar dlnff      ("DlnFF"     ,"DlnFF"    ,"");   // FF of Dlnu: rho^2
    RooStringVar pi0weights ("Pi0Weights","Pi0Weights","");
    RooStringVar sigmc      ("SigMC"     ,"SigMC"    ,"");
    RooArgSet weightpars(variweights,bkgweights,plweightsD,bfweights,ffweights,dlnff,pi0weights,sigmc);
    if (weightpars.readFromFile(weightFile)) {
      cout << "Error! could not read a valid reweighting configuration from the file " << weightFile << endl;
      cout << "No reweighting will be done" << endl;
    } else {
      initializeVariWeights(variweights.getVal());
      initializePlWeights(plweightsD.getVal());
      initializeBkgWeights(bkgweights.getVal());
      if(IsCocktail>=0){
	initializeBFWeights(bfweights.getVal());
	initializeFFWeights(ffweights.getVal());
	initializeDlnFFWeights(dlnff.getVal());
	initializePi0Weights(pi0weights.getVal());
	initializeSigMC(sigmc.getVal());
      }
    }
  } else cout<<"Weight file is not valid"<<endl;
}
void WeightManager::Print(){
  cout<<"BF "<<doBFWeights<<", D* FF "<<doFFWeights<<", D FF "<<doDlnFFWeights<<", SigMC "<<doSigMC<<endl;
  cout<<"Is cocktail? "<<isCocktail<<endl;
}

void WeightManager::initializeFFWeights(const char *file)
{
  if (!file) return;
  if (strlen(file) == 0) return;
  if(verbose) cout << "initializing FF weights with file " << file << endl;
  doFFWeights = kTRUE;
  RooRealVar rrvr0("HQET_R0","HQET_R0",1.14);       
  RooRealVar rrvr1("HQET_R1","HQET_R1",1.401);       
  RooRealVar rrvr2("HQET_R2","HQET_R2",0.864);       
  RooRealVar rrvrho2("HQET_rho2","HQET_rho2",1.214); 
  RooRealVar rrvgsr("gSR","gSR",0); 

  RooArgSet ffweights(rrvr0,rrvr1,rrvr2,rrvrho2,rrvgsr);
  if (ffweights.readFromFile(file)) {
    cout << "Error reading FF weights from file " << file << endl;
    cout << "weights not updated" << endl;
    return;
  }
  HQET_R0 = rrvr0.getVal();
  HQET_R1 = rrvr1.getVal();
  HQET_R2 = rrvr2.getVal();
  HQET_rho2 = rrvrho2.getVal();
  gSR = rrvgsr.getVal();
}

void WeightManager::initializeDlnFFWeights(const char *file){
  if (!file) return;
  if (strlen(file) == 0) return;
  if(verbose) cout << "initializing Dlnu FF weights with file " << file << endl;
  doDlnFFWeights = kTRUE;
  RooRealVar rrvrhod2("rhoD2","rhoD2",1.18);    
  RooRealVar rrvdelta("Delta","Delta",1.);    
  RooRealVar rrvgsr("gSR","gSR",0); 

  RooArgSet dlnffweights(rrvrhod2,rrvdelta,rrvgsr);
  if (dlnffweights.readFromFile(file)) {
    cout << "Error reading Dlnu FF weights from file " << file << endl;
    cout << "weights not updated" << endl;
    return;
  }
  rhoD2 = rrvrhod2.getVal();
  Delta = rrvdelta.getVal();
  gSR = rrvgsr.getVal();
}

void WeightManager::initializeVariWeights(const char *folder){
  if (strlen(folder)==0) return;
  doVariWeights = kTRUE;
  nVari = 2;
  Variables[0] = "candEExtra"; Variables[1] = "candEExtranoMVA"; 
  //Variables[1] = "candMES"; Variables[2] = "candDeltaE";  Variables[3] = "candTagChargedMult";
  for(int vari=0; vari<nVari; vari++){
    TString wFileName(folder); wFileName+=Variables[vari]; wFileName+=".txt";
    fstream wFile;
    wFile.open(wFileName,fstream::in);
    wFile>>nbins[vari];
    for(int bin=0; bin<nbins[vari]; bin++){
      wFile>>ranges[vari][bin];
      for(int chan=0; chan<4; chan++) wFile>>variWeights[vari][chan][bin];
    }
    wFile>>ranges[vari][nbins[vari]];
    if(!wFileName.Contains("noMVA")) wFile>>BDTBias;
  }
}

void WeightManager::initializeBkgWeights(const char *file){
  if (!file) return;
  if (strlen(file)==0) return;
  if(verbose) cout << "initializing Bkg weights with file "<<file<<endl;
  doBkgWeights = kTRUE;

  fstream BkgFile;
  Float_t dummy;
  BkgFile.open(file,fstream::in);
  for(int cand=0; cand<4; cand++){
    BkgFile>>dummy>>wBkg[cand];
  }
  fstream mESBkgFile;
  TString NamemES = file; NamemES += "mES";
  mESBkgFile.open(NamemES,fstream::in);
  for(int cand=0; cand<4; cand++){
    mESBkgFile>>dummy>>wBkgmES[cand];
  }
  
}

void WeightManager::initializePlWeights(const char *file){
  if (!file) return;
  if (strlen(file)==0) return;
  if(verbose) cout << "initializing p*l weights with file "<<file<<endl;
  doPlWeights = kTRUE;

  fstream PlFileD;
  Float_t dummy;
  PlFileD.open(file,fstream::in);
  for(int cand=0; cand<4; cand++){
    for(int i=0; i<1000; i++){
      PlFileD>>dummy>>wPlD[cand][i];
    }
  }
}


void WeightManager::initializePi0Weights(const char *file)
{
  if (!file) return;
  if (strlen(file) == 0) return;
  if(verbose) cout << "initializing Pi0 weights with file " << file << endl;
  doPi0Weights = kTRUE;
  RooRealVar rrvint("intercept","intercept",1.0);
  RooRealVar rrvslope("slope","slope",0.0);
  RooArgSet pi0weights(rrvint,rrvslope);
  if (pi0weights.readFromFile(file)) {
    cout << "Error reading Pi0 weights from file " << file << endl;
    cout << "weights not updated" << endl;
    return;
  }
  pi0int = rrvint.getVal();
  pi0slope = rrvslope.getVal();
}

void WeightManager::initializeSigMC(const char *msg)
{
  if (strlen(msg) == 0) return;
  if (sscanf(msg,"%d",&sigMode) == 1)
    doSigMC = kTRUE;
  if(verbose) cout << "doing SigMCWeights " << sigMode << endl;
}

WeightManager::~WeightManager(){}


Float_t WeightManager::getUnsimWeight(Int_t mode)
{
  if (sigMode == 0) return 0;
  //normalized to D(*)lnu yields in generic MC
  //ex: D0lnu reconstructed in D0 channel in 970 fb: 1790 events
  //    1762 are in modes common to signal MC
  //    D0lnu->D0 in signal MC (reweighted with SigMC mode 9): sum = 10.8679
  //    resulting weight = 10.8679/1762 = 0.00617
  Float_t f = 0.0;
  if (mode == 1||mode==5) f = 0.00617;
  if (mode == 2||mode==6) f = 0.00605;
  if (mode == 3||mode==7) f = 0.02020;
  if (mode == 4||mode==8) f = 0.01492;
  if (doBFWeights) {
    f *= wDssGlobal;
  }
  return f;
}

Float_t WeightManager::Pi0Weight(Int_t MCType, Int_t candType, Int_t candDstarType, Float_t ppi0, Float_t dssppi0){
  Float_t res = 1.0;
  Float_t tw;
  if (ppi0 > 0) {
    tw = pi0int + pi0slope*ppi0;
    if (candType==2||candType==4||candType==6||candType==8)
      res *= tw;
    else res /= tw;
  }
  if (dssppi0 > 0) {
    tw = pi0int + pi0slope*dssppi0;
    if (candType > 4)
      res *= tw;
    else res /= tw;
    if (candType < 5 && dssppi0 < 0.3) res *= corrSoftPi;
  }
  return res;
}


void WeightManager::initializeBFWeights(const char *file){
  wDssGlobal = 1.0;
  if (!file) return;
  if (strlen(file) == 0) return;
  if(verbose) cout << "initializing BF weights with file " << file << endl;
  doBFWeights = kTRUE;

  RooRealVar rrvD0("weightD0","weightD0",1.0);
  RooRealVar rrvDs0("weightDs0","weightDs0",1.0);
  RooRealVar rrvDp("weightDp","weightDp",1.0);
  RooRealVar rrvDsp("weightDsp","weightDsp",1.0);

  RooRealVar rrvD10("weightD10","weightD10",1.0);
  RooRealVar rrvD20("weightD20","weightD20",1.0);
  RooRealVar rrvD00("weightD00","weightD00",1.0);
  RooRealVar rrvD1prime0("weightD1prime0","weightD1prime0",1.0);
  RooRealVar rrvDstarpi("weightDstarpi","weightDstarpi",1.0);
  RooRealVar rrvDstar0pi0("weightDstar0pi0","weightDstar0pi0",1.0);
  RooRealVar rrvDpi("weightDpi","weightDpi",1.0);
  RooRealVar rrvD0pi0("weightD0pi0","weightD0pi0",1.0);
  RooRealVar rrvD1p("weightD1p","weightD1p",1.0);
  RooRealVar rrvD2p("weightD2p","weightD2p",1.0);
  RooRealVar rrvD0p("weightD0p","weightD0p",1.0);
  RooRealVar rrvD1primep("weightD1primep","weightD1primep",1.0);
  RooRealVar rrvDstarpi0("weightDstarpi0","weightDstarpi0",1.0);
  RooRealVar rrvDstar0pi("weightDstar0pi","weightDstar0pi",1.0);
  RooRealVar rrvDpi0("weightDpi0","weightDpi0",1.0);
  RooRealVar rrvD0pi("weightD0pi","weightD0pi",1.0);

  RooRealVar rrvwDoubleSL("wDoubleSL","wDoubleSL",1.0);
  RooRealVar rrvwFakeLep("wFakeLep","wFakeLep",1.0);
  RooRealVar rrvwUnknown("wUnknown","wUnknown",1.0);

  RooRealVar rrvwDsTaunu("wDsTaunu","wDsTaunu",1.0);
  RooRealVar rrvwDsEtalnu("wDsEtalnu","wDsEtalnu",1.0);
  RooRealVar rrvwDsEtaprimelnu("wDsEtaprimelnu","wDsEtaprimelnu",1.0);
  RooRealVar rrvwDsphilnu("wDsphilnu","wDsphilnu",1.0);
  RooRealVar rrvwDsmunu("wDsmunu","wDsmunu",1.0);

  RooRealVar rrvwDzbDs("wDzbDs","wDzbDs",1.0);
  RooRealVar rrvwDzbDsstar("wDzbDsstar","wDzbDsstar",1.0);
  RooRealVar rrvwDstarzbDs("wDstarzbDs","wDstarzbDs",1.0);
  RooRealVar rrvwDstarzbDsstar("wDstarzbDsstar","wDstarzbDsstar",1.0);
  RooRealVar rrvwDstarstarDs("wDstarstarDs","wDstarstarDs",1.0);
  RooRealVar rrvwDpDzbKz("wDpDzbKz","wDpDzbKz",1.0);
  RooRealVar rrvwDzbDstarpKz("wDzbDstarpKz","wDzbDstarpKz",1.0);
  RooRealVar rrvwDstarzbDpKz("wDstarzbDpKz","wDstarzbDpKz",1.0);
  RooRealVar rrvwDstarzbDstarpKz("wDstarzbDstarpKz","wDstarzbDstarpKz",1.0);
  RooRealVar rrvwDzDzbKp("wDzDzbKp","wDzDzbKp",1.0);
  RooRealVar rrvwDstarzbDstarzKp("wDstarzbDstarzKp","wDstarzbDstarzKp",1.0);
  RooRealVar rrvwDstarzDzbKp("wDstarzDzbKp","wDstarzDzbKp",1.0);
  RooRealVar rrvwDstarmDpKp("wDstarmDpKp","wDstarmDpKp",1.0);
  RooRealVar rrvwDzba1p("wDzba1p","wDzba1p",1.0);
  RooRealVar rrvwDstarzba1p("wDstarzba1p","wDstarzba1p",1.0);
  RooRealVar rrvwDzbRhop("wDzbRhop","wDzbRhop",1.0);
  RooRealVar rrvwDstarzbRhop("wDstarzbRhop","wDstarzbRhop",1.0);

  RooRealVar rrvwDmDs("wDmDs","wDmDs",1.0);
  RooRealVar rrvwDstarmDs("wDstarmDs","wDstarmDs",1.0);
  RooRealVar rrvwDmDsstar("wDmDsstar","wDmDsstar",1.0);
  RooRealVar rrvwDstarmDsstar("wDstarmDsstar","wDstarmDsstar",1.0);
  RooRealVar rrvwDmDzKp("wDmDzKp","wDmDzKp",1.0);
  RooRealVar rrvwDmDstarzKp("wDmDstarzKp","wDmDstarzKp",1.0);
  RooRealVar rrvwDstarmDstarzKp("wDstarmDstarzKp","wDstarmDstarzKp",1.0);
  RooRealVar rrvwDpDmKz("wDpDmKz","wDpDmKz",1.0);
  RooRealVar rrvwDstarmDpKz("wDstarmDpKz","wDstarmDpKz",1.0);
  RooRealVar rrvwDstarmDstarpKz("wDstarmDstarpKz","wDstarmDstarpKz",1.0);
  RooRealVar rrvwDstarzDzbKz("wDstarzmDzbKz","wDstarzDzbKz",1.0);
  RooRealVar rrvwDstarpDstarm("wDstarpDstarm","wDstarpDstarm",1.0);
  RooRealVar rrvwDpDstarm("wDpDstarm","wDpDstarm",1.0);
  RooRealVar rrvwDpDm("wDpDm","wDpDm",1.0);
  RooRealVar rrvwDma1p("wDma1p","wDma1p",1.0);
  RooRealVar rrvwDstarma1p("wDstarma1p","wDstarma1p",1.0);
  RooRealVar rrvwDzbRhoz("wDzbRhoz","wDzbRhoz",1.0);
  RooRealVar rrvwDmRhop("wDmRhop","wDmRhop",1.0);
  RooRealVar rrvwDstarzbRhoz("wDstarzbRhoz","wDstarzbRhoz",1.0);
  RooRealVar rrvwDstarmRhop("wDstarmRhop","wDstarmRhop",1.0);

  RooArgSet bfweights(rrvD10,rrvD20,rrvD00,rrvD1prime0,rrvDstarpi,rrvDstar0pi0,rrvDpi,rrvD0pi0);
  bfweights.add(rrvD1p);
  bfweights.add(rrvD2p);
  bfweights.add(rrvD0p);
  bfweights.add(rrvD1primep);
  bfweights.add(rrvDstarpi0);
  bfweights.add(rrvDstar0pi);
  bfweights.add(rrvDpi0);
  bfweights.add(rrvD0pi);

  bfweights.add(rrvD0);
  bfweights.add(rrvDs0);
  bfweights.add(rrvDp);
  bfweights.add(rrvDsp);

  bfweights.add(rrvwDoubleSL);
  bfweights.add(rrvwFakeLep);
  bfweights.add(rrvwUnknown);
  bfweights.add(rrvwDsTaunu);
  bfweights.add(rrvwDsEtalnu);
  bfweights.add(rrvwDsEtaprimelnu);
  bfweights.add(rrvwDsphilnu);
  bfweights.add(rrvwDsmunu);
  bfweights.add(rrvwDzbDs);
  bfweights.add(rrvwDzbDsstar);
  bfweights.add(rrvwDstarzbDs);
  bfweights.add(rrvwDstarzbDsstar);
  bfweights.add(rrvwDstarstarDs);
  bfweights.add(rrvwDpDzbKz);
  bfweights.add(rrvwDzbDstarpKz);
  bfweights.add(rrvwDstarzbDpKz);
  bfweights.add(rrvwDstarzbDstarpKz);
  bfweights.add(rrvwDzDzbKp);
  bfweights.add(rrvwDstarzbDstarzKp);
  bfweights.add(rrvwDstarzDzbKp);
  bfweights.add(rrvwDstarmDpKp);
  bfweights.add(rrvwDzba1p);
  bfweights.add(rrvwDstarzba1p);
  bfweights.add(rrvwDzbRhop);
  bfweights.add(rrvwDstarzbRhop);

  bfweights.add(rrvwDmDs);
  bfweights.add(rrvwDstarmDs);
  bfweights.add(rrvwDmDsstar);
  bfweights.add(rrvwDstarmDsstar);
  bfweights.add(rrvwDmDzKp);
  bfweights.add(rrvwDmDstarzKp);
  bfweights.add(rrvwDstarmDstarzKp);
  bfweights.add(rrvwDpDmKz);
  bfweights.add(rrvwDstarmDpKz);
  bfweights.add(rrvwDstarmDstarpKz);
  bfweights.add(rrvwDstarzDzbKz);
  bfweights.add(rrvwDstarpDstarm);
  bfweights.add(rrvwDpDstarm);
  bfweights.add(rrvwDpDm);
  bfweights.add(rrvwDma1p);
  bfweights.add(rrvwDstarma1p);
  bfweights.add(rrvwDzbRhoz);
  bfweights.add(rrvwDmRhop);
  bfweights.add(rrvwDstarzbRhoz);
  bfweights.add(rrvwDstarmRhop);

  if (bfweights.readFromFile(file)) {
    cout << "Error reading BF weights from file " << file << endl;
    cout << "weights not updated" << endl;
    return;
  }

  wD0 = rrvD0.getVal();
  wDs0 = rrvDs0.getVal();
  wDp = rrvDp.getVal();
  wDsp = rrvDsp.getVal();

  wD10 = rrvD10.getVal();
  wD20 = rrvD20.getVal();
  wD00 = rrvD00.getVal();
  wD1prime0 = rrvD1prime0.getVal();
  wDstarpi = rrvDstarpi.getVal();
  wDstar0pi0 = rrvDstar0pi0.getVal();
  wDpi = rrvDpi.getVal();
  wD0pi0 = rrvD0pi0.getVal();
  wD1p = rrvD1p.getVal();
  wD2p = rrvD2p.getVal();
  wD0p = rrvD0p.getVal();
  wD1primep = rrvD1primep.getVal();
  wDstarpi0 = rrvDstarpi0.getVal();
  wDstar0pi = rrvDstar0pi.getVal();
  wDpi0 = rrvDpi0.getVal();
  wD0pi = rrvD0pi.getVal(); 

  wDoubleSL = rrvwDoubleSL.getVal();
  wFakeLep = rrvwFakeLep.getVal();
  wUnknown = rrvwUnknown.getVal();
  wDsTaunu = rrvwDsTaunu.getVal();
  wDsEtalnu = rrvwDsEtalnu.getVal();
  wDsEtaprimelnu = rrvwDsEtaprimelnu.getVal();
  wDsphilnu = rrvwDsphilnu.getVal();
  wDsmunu = rrvwDsmunu.getVal();
  wDzbDs = rrvwDzbDs.getVal();
  wDzbDsstar = rrvwDzbDsstar.getVal();
  wDstarzbDs = rrvwDstarzbDs.getVal();
  wDstarzbDsstar = rrvwDstarzbDsstar.getVal();
  wDstarstarDs = rrvwDstarstarDs.getVal();
  wDpDzbKz = rrvwDpDzbKz.getVal();
  wDzbDstarpKz = rrvwDzbDstarpKz.getVal();
  wDstarzbDpKz = rrvwDstarzbDpKz.getVal();
  wDstarzbDstarpKz = rrvwDstarzbDstarpKz.getVal();
  wDzDzbKp = rrvwDzDzbKp.getVal();
  wDstarzbDstarzKp = rrvwDstarzbDstarzKp.getVal();
  wDstarzDzbKp = rrvwDstarzDzbKp.getVal();
  wDstarmDpKp = rrvwDstarmDpKp.getVal();
  wDzba1p = rrvwDzba1p.getVal();
  wDstarzba1p = rrvwDstarzba1p.getVal();
  wDzbRhop = rrvwDzbRhop.getVal();
  wDstarzbRhop = rrvwDstarzbRhop.getVal();

  wDmDs = rrvwDmDs.getVal();
  wDstarmDs = rrvwDstarmDs.getVal();
  wDmDsstar = rrvwDmDsstar.getVal();
  wDstarmDsstar = rrvwDstarmDsstar.getVal();
  wDmDzKp = rrvwDmDzKp.getVal();
  wDmDstarzKp = rrvwDmDstarzKp.getVal();
  wDstarmDstarzKp = rrvwDstarmDstarzKp.getVal();
  wDpDmKz = rrvwDpDmKz.getVal();
  wDstarmDpKz = rrvwDstarmDpKz.getVal();
  wDstarmDstarpKz = rrvwDstarmDstarpKz.getVal();
  wDstarzDzbKz = rrvwDstarzDzbKz.getVal();
  wDstarpDstarm = rrvwDstarpDstarm.getVal();
  wDpDstarm = rrvwDpDstarm.getVal();
  wDpDm = rrvwDpDm.getVal();
  wDstarma1p = rrvwDstarma1p.getVal();
  wDma1p = rrvwDma1p.getVal();
  wDzbRhoz = rrvwDzbRhoz.getVal();
  wDmRhop = rrvwDmRhop.getVal();
  wDstarzbRhoz = rrvwDstarzbRhoz.getVal();
  wDstarmRhop = rrvwDstarmRhop.getVal();

  wDssGlobal = (wD10 + wD20 + wD00 + wD1prime0 + wDstarpi + wDstar0pi0 + wDpi + wD0pi0 +
		wD1p + wD2p + wD0p + wD1primep + wDstarpi0 + wDstar0pi + wDpi0 + wD0pi) / 16.0;

  fstream BFIter; BFIter.open("../DonutUtils/BFIteration.txt",fstream::in);
  TString dummy; 
  for(int i=0; i<8; i++) BFIter>>dummy>>BFIterRatio[i];
}

