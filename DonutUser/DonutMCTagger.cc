#include "BaBar/BaBar.hh"

#include "DonutUser/DonutMCTagger.hh"

#include <assert.h>
#include <iostream>
#include <math.h>

#include "AbsEnv/AbsEnv.hh"
#include "AbsEvent/AbsEvent.hh"
#include "AbsEventTag/AbsEventTag.hh"
#include "Beta/EventInfo.hh"
#include "BetaCoreTools/BtaPrintTree.hh"
#include "BetaCoreTools/BtaMcAssoc.hh"
#include "CLHEP/Alist/AIterator.h"
#include "ErrLogger/ErrLog.hh"
#include "GenEnv/GenEnv.hh"
#include "HepTuple/TupleManager.h"
#include "HepTuple/Histogram.h"
#include "ProxyDict/IfdStrKey.hh"
#include "ProxyDict/Ifd.hh"
#include "XslFFReweighting/XSLKin.hh"

using std::cout;
using std::endl;

DonutMCTagger::DonutMCTagger( const char* const theName, 
			const char* const theDescription )
  : AppModule( theName, theDescription ),
    _truthList ("MCTruthList"      , this, "MCTruth"), 
    _Verbose("Verbose",this,false)
{
  commands()->append( &_truthList );
  commands()->append( &_Verbose );
}

DonutMCTagger::~DonutMCTagger( )
{
}

AppResult
DonutMCTagger::endJob( AbsEvent* anEvent )
{
  return AppResult::OK;
}

AppResult
DonutMCTagger::beginJob( AbsEvent* anEvent )
{
  HepTupleManager* manager = gblEnv->getGen()->ntupleManager();
  assert (manager!=0);
  mcTypeHisto = manager -> histogram("mcType",15,-0.5,14.5);

  return AppResult::OK;
}

AppResult
DonutMCTagger::event( AbsEvent* anEvent )
{
  if(_Verbose.value())    cout<<"MCTagger manuelf: We start event"<<endl;
  tagMC(anEvent);
  AbsEventTag* tag = Ifd<AbsEventTag>::get( anEvent );
  tag -> putInt  (MCType      , "MCType"    );
  tag -> putInt  (MCSubmode   , "MCSubmode" );
  tag -> putInt  (MCCombmode  , "MCCombmode");
  tag -> putInt  (MCPions     , "MCPions"   );// Number of pions from SigB and Dss decay. +10 if pi0
  tag -> putInt  (MCTaumode   , "MCTaumode" );
  tag -> putInt  (MCBrem      , "MCBrem"    );
  tag -> putInt  (MCScatter   , "MCScatter" );
  tag -> putInt  (MCD         , "MCD"       );
  tag -> putInt  (MC2Body     , "MC2Body"   );
  tag -> putInt  (MCDoubleSL  , "MCDoubleSL");
  tag -> putInt  (nTrueP      , "nTrueP"    );
  tag -> putInt  (nTrueKL     , "nTrueKL"   );
  tag -> putInt  (nTrueNu     , "nTrueNu"   );
  tag -> putInt  (uidBTag     , "uidBTag"   );
  tag -> putInt  (uidD        , "uidD"      );
  tag -> putInt  (uidLep      , "uidLep"    );
  tag -> putInt  (uidNu       , "uidNu"     );
  tag -> putInt  (uidPi0      , "uidPi0"    );
  tag -> putInt  (uidDssPi0   , "uidDssPi0" );
  tag -> putFloat(truePLep    , "truePLep"  );
  tag -> putFloat(trueTauE    , "trueTauE"  );
  tag -> putFloat(trueLepE    , "trueLepE"  );
  tag -> putFloat(truePD      , "truePD"    );
  tag -> putFloat(trueDmass   , "trueDmass" );
  tag -> putFloat(trueCTL     , "trueCTL"   );
  tag -> putFloat(trueCTV     , "trueCTV"   );
  tag -> putFloat(trueChi     , "trueChi"   );
  tag -> putFloat(trueQ2      , "trueQ2"    );
  tag -> putFloat(trueCTL2    , "trueCTL2"  );
  tag -> putInt  (trueLepCharge,"trueLepCharge");
  tag -> putFloat(trueTauFlight,"trueTauFlight");
  tag -> putFloat(trueMM2     , "trueMM2"   );
  tag -> putInt  (isVub       , "isVub"     ); // 1 if there is a lepton and not a D
  tag -> putInt  (isBzero     , "isBzero"   );
  tag -> putInt  (MCDssmode   , "MCDssmode" );
  tag -> putInt  (MCUnsim     , "MCUnsim"   );
  // manuelf variables
  tag -> putInt  (nMufromPi     ,  "nMufromPi");
  tag -> putInt  (ntrueLep      ,  "ntrueLep");
  tag -> putInt  (trueMotherLep1,  "trueMotherLep1");
  tag -> putInt  (trueMotherLep2,  "trueMotherLep2");
  tag -> putInt  (trueMotherLep3,  "trueMotherLep3");
  tag -> putInt  (trueLundLep1  ,  "trueLundLep1");
  tag -> putInt  (trueLundLep2  ,  "trueLundLep2");
  tag -> putInt  (trueLundLep3  ,  "trueLundLep3");
  tag -> putFloat(truePLep1     ,  "truePLep1");
  tag -> putFloat(truePLep2     ,  "truePLep2");
  tag -> putFloat(truePLep3     ,  "truePLep3");

  tag -> putInt  (trueTagPi      ,  "trueTagPi");
  tag -> putInt  (trueTagK       ,  "trueTagK"); 
  tag -> putInt  (trueTagKs      ,  "trueTagKs");
  tag -> putInt  (trueTagPi0     ,  "trueTagPi0");
  tag -> putInt  (trueTagE       ,  "trueTagE"); 
  tag -> putInt  (trueTagMu      ,  "trueTagMu");
  tag -> putInt  (trueTagDs      ,  "trueTagDs");
  tag -> putInt  (trueTagD       ,  "trueTagD"); 
  tag -> putInt  (trueTagYPi     ,  "trueTagYPi");
  tag -> putInt  (trueTagYK      ,  "trueTagYK"); 
  tag -> putInt  (trueTagYKs     ,  "trueTagYKs");
  tag -> putInt  (trueTagYPi0    ,  "trueTagYPi0");
  tag -> putInt  (trueTagElse    ,  "trueTagElse");

  mcTypeHisto -> accumulate(MCType);
  if(_Verbose.value())    cout<<"MCTagger manuelf: We finish event"<<endl;

  return AppResult::OK;
}

void DonutMCTagger::tagMC(AbsEvent* anEvent)
{
  truthList = Ifd<HepAList< BtaCandidate > >::get(anEvent, _truthList.value());

  int ret = 0;
  int DType = 0;
  int lepType = 0;
  int pion = 0;
  int nE = 0;
  int nMu = 0;
  int nTau = 0;
  int tauFromDs = 0;
  int lepFromD = 0;
  
  int eLund = (int)PdtLund::e_minus;
  int muLund = (int)PdtLund::mu_minus;
  int piLund = (int)PdtLund::pi_plus;
  int pi0Lund = (int)PdtLund::pi0;
  int KLund = (int)PdtLund::K_plus;

  MCType = MCSubmode = MCCombmode = MCPions = MCTaumode = MCBrem =
    MCScatter = MCD = MCDoubleSL = nTrueP = nTrueKL = nTrueNu = isBzero = MCDssmode = MCUnsim = -1;
  isVub = 0;
  trueLepCharge = 0;
  truePLep = truePD = trueDmass = trueCTL = trueCTV = trueChi = trueQ2 = trueTauFlight = 
    trueTauE = trueLepE = trueCTL2 = trueMM2 = -999.;
  ntrueLep = trueMotherLep1 = trueMotherLep2 = trueMotherLep2 = trueLundLep1 =
    trueMotherLep2 = trueMotherLep3 = nMufromPi = 0;
  truePLep1 = truePLep2 = truePLep3 = 0.;
  trueTagPi= trueTagK= trueTagKs= trueTagPi0= trueTagE= trueTagMu= -9;
  trueTagDs= trueTagD=0;
  trueTagYPi= trueTagYK= trueTagYKs= trueTagYPi0= trueTagElse-9;


  if (!truthList) return;
  if (truthList -> length() < 4) return;

  MCType =  MCSubmode = MCCombmode = MCPions = MCBrem = MCScatter = MCD = 
    MC2Body = MCDoubleSL = nTrueP = nTrueKL = nTrueNu = 0;
  uidBTag = uidD = uidLep = uidNu = uidPi0 = uidDssPi0 = -1;

  printTree = new BtaPrintTree();
  printTree -> addTerminal(PdtLund::pi0);
  printTree -> addTerminal(PdtLund::pi_plus);
  printTree -> addTerminal(PdtLund::pi_minus);
  printTree -> addTerminal(PdtLund::K_plus);
  printTree -> addTerminal(PdtLund::K_minus);
  printTree -> addTerminal(PdtLund::mu_plus);
  printTree -> addTerminal(PdtLund::mu_minus);
  printTree -> setOutputFunc(BtaPrintTree::printName);

  BtaCandidate *aTru, *theB(0), *tagB(0);
  int nSL = 0;
  HepAListIterator<BtaCandidate> iterTruth(*truthList);
  if(_Verbose.value()) cout<<"MCTagger manuelf: Loop over all true particles in tagMC"<<endl;
  while (aTru = iterTruth()) {
    int ParticleLund = (int)aTru -> pdtEntry() -> lundId();
    // No need to continue if the particle comes from a final particle
    if (aTru->theMother()) {
      int MotherLund = abs(aTru->theMother()->pdtEntry()->lundId());
      if(MotherLund==piLund && abs(ParticleLund) == (int)PdtLund::mu_minus) nMufromPi++;
      if(MotherLund==eLund||MotherLund==muLund||MotherLund==piLund||MotherLund==pi0Lund||MotherLund==KLund)
	continue;
    }
      
    //If we have a B, set isBzero, check if it has a D and a lepton
    if (ParticleLund == (int)PdtLund::B_minus ||
	ParticleLund == (int)PdtLund::B_plus ||
	ParticleLund == (int)PdtLund::B0 ||
	ParticleLund == (int)PdtLund::anti_B0) {
      if (_Verbose.value()){
	cout<<printTree->print(*aTru);
	cout<<endl;
      }
      if (ParticleLund == (int)PdtLund::B_minus ||
	  ParticleLund == (int)PdtLund::B_plus)
	isBzero = 0;
      else isBzero = 1;
      //Check if it has a D and a lepton
      if (isSemiLep(aTru)) {
	theB = aTru;
	nSL ++;
      } else {
	uidBTag = (int)aTru->uid();
	tagB = aTru;
      }
      if (aTru -> nDaughters() == 2) {
	if (MC2Body > 0) MC2Body += 1000 * TwoBodyMode(aTru);
	else MC2Body = TwoBodyMode(aTru);
      }
    }
    //Set number of protons, klongs, neutrinos close to the vertex
    if (abs(ParticleLund) == (int)PdtLund::p_plus) {
      if (aTru->theMother()->decayVtx()->point().mag() < 10) nTrueP ++;
    } //proton
    if (abs(ParticleLund) == (int)PdtLund::K_L0) {
      if (aTru->theMother()->decayVtx()->point().mag() < 10) nTrueKL ++;
    } //klong
    if (abs(ParticleLund) == (int)PdtLund::nu_e  ||
	abs(ParticleLund) == (int)PdtLund::nu_mu ||
	abs(ParticleLund) == (int)PdtLund::nu_tau) {
      if (aTru->theMother()->decayVtx()->point().mag() < 10) nTrueNu ++;
    } //neutrino
    if (abs(ParticleLund) == eLund) { // elec not coming from tau
      if (aTru->theMother()) {
	if (abs((int)aTru->theMother()->pdtEntry()->lundId()) != (int)PdtLund::tau_minus)
	  if (aTru->theMother()->decayVtx()->point().mag() < 10) {
	    nE ++;
	    int mLund = abs((int)aTru->theMother()->pdtEntry()->lundId());
	    if (mLund==(int)PdtLund::D_plus||mLund==(int)PdtLund::D0||mLund==(int)PdtLund::D_s_plus) lepFromD = 1;
	  }
      }
    } 
    if (abs(ParticleLund) == (int)PdtLund::mu_minus) { //muon not coming from tau
      if (aTru->theMother()) {
	if (abs((int)aTru->theMother()->pdtEntry()->lundId()) != (int)PdtLund::tau_minus)
	  if (aTru->theMother()->decayVtx()->point().mag() < 10) {
	    nMu ++;
	    int mLund = abs((int)aTru->theMother()->pdtEntry()->lundId());
	    if (mLund==(int)PdtLund::D_plus||mLund==(int)PdtLund::D0||mLund==(int)PdtLund::D_s_plus) lepFromD = 1;
	  }
      }
    }
    if (abs(ParticleLund) == (int)PdtLund::tau_minus) { //tau
      if (aTru->theMother()) {
	if (aTru->theMother()->decayVtx()->point().mag() < 10) nTau ++;
	if (abs((int)aTru->theMother()->pdtEntry()->lundId()) == (int)PdtLund::D_s_plus) tauFromDs = 1;
      }
    }
    if(aTru->theMother() && ((abs(ParticleLund)==(int)PdtLund::mu_minus) || (abs(ParticleLund)==eLund))){
      int MotherLund = aTru->theMother()->pdtEntry()->lundId();
      if (aTru->theMother()->decayVtx()->point().mag() < 10) {
	ntrueLep++;
	switch(ntrueLep){
	case 1:
	  truePLep1 = aTru->p();
	  trueLundLep1 = ParticleLund;
	  trueMotherLep1 = MotherLund;
	  break;
	case 2:
	  truePLep2 = aTru->p();
	  trueLundLep2 = ParticleLund;
	  trueMotherLep2 = MotherLund;
	  break;
	case 3:
	  truePLep3 = aTru->p();
	  trueLundLep3 = ParticleLund;
	  trueMotherLep3 = MotherLund;
	  break;
	}
      }
    }
  }//end while
  iterTruth.rewind();

  // Setting of MCCombmode, types of combinatoric background
  if (tauFromDs) MCCombmode = 1;
  else if (nTau) MCCombmode = 2;
  else if (nE + nMu == 2) {
    if (nE == 2 && nMu == 0) MCCombmode = 3;
    else if (nE == 0 && nMu == 2) MCCombmode = 4;
    else if (nE == 1 && nMu == 1) MCCombmode = 5;
  }
  else if (nE + nMu == 1) {
    MCCombmode = 8;
    if (nMu == 1) MCCombmode += 1;
    if (lepFromD) MCCombmode -= 2;
  }
  else if (nE + nMu > 0) {
    MCCombmode = 10;
  }
  else if (nTrueKL) MCCombmode = 11;
  else MCCombmode = 12;

  if (nSL == 2) MCDoubleSL = 1;
  
  // Looping over the tag B
  if(tagB==0) return;
  int kcount = 0, kscount = 0, picount = 0, pi0count = 0, ecount = 0, mucount = 0;
  bool result = count(tagB,kcount,kscount,picount,pi0count,ecount,mucount);
  trueTagElse = (int) result;
  trueTagYPi = picount; trueTagYK = kcount; trueTagYKs = kscount;  trueTagYPi0 = pi0count;
  HepAListIterator<BtaCandidate> iterDau = tagB -> daughterIterator();
  BtaCandidate* aDaughter;
  if(_Verbose.value())  cout<<"MCTagger manuelf: Loop over the daughters of the tag B. ret "<<ret<<endl;
  while (aDaughter = iterDau()) {
    int dauLund = (int)aDaughter->pdtEntry()->lundId();
    if (abs(dauLund) == (int)PdtLund::D0 || 
	abs(dauLund) == (int)PdtLund::D_plus || 
	abs(dauLund) == (int)PdtLund::J_psi || 
	abs(dauLund) == (int)PdtLund::D_s_plus) {
      if(trueTagDs){
	trueTagElse = false;
	continue;
      }
      trueTagDs = trueTagD = dauLund;
      kcount = 0; kscount = 0; picount = 0; pi0count = 0; ecount = 0; mucount = 0;
      result = count(aDaughter,kcount,kscount,picount,pi0count,ecount,mucount);
      trueTagPi = picount; trueTagK = kcount; trueTagKs = kscount; trueTagPi0 = pi0count;
      trueTagE = ecount; trueTagMu = mucount;
      if(result==false) trueTagElse = 0;
    }
    if (abs(dauLund) == (int)PdtLund::D_star0 || 
	abs(dauLund) == (int)PdtLund::D_star_plus || 
	abs(dauLund) == (int)PdtLund::D_s_star_plus) {
      trueTagDs = dauLund;
      HepAListIterator<BtaCandidate> iterDDau = aDaughter -> daughterIterator();
      BtaCandidate* DDau;
      while (DDau = iterDDau()) {
	int DdauLund = (int)DDau->pdtEntry()->lundId();
	if (abs(DdauLund) == (int)PdtLund::D0 || 
	    abs(DdauLund) == (int)PdtLund::D_plus || 
	    abs(DdauLund) == (int)PdtLund::D_s_plus) {
	  if(trueTagD){
	    trueTagElse = false;
	    continue;
	  }
	  trueTagD = (int)DDau->pdtEntry()->lundId();
	  kcount = 0; kscount = 0; picount = 0; pi0count = 0; ecount = 0; mucount = 0;
	  result = count(DDau,kcount,kscount,picount,pi0count,ecount,mucount);
	  if(result==false) trueTagElse = 0;
	  trueTagPi = picount; trueTagK = kcount; trueTagKs = kscount; trueTagPi0 = pi0count;
	  trueTagE = ecount; trueTagMu = mucount;
	}
	if (abs(DdauLund) == (int)PdtLund::pi_plus)
	  trueTagYPi -= 1;
	if (abs(DdauLund) == (int)PdtLund::pi0)
	  trueTagYPi0 -= 1;
      }
      iterDDau.rewind();
    }
  }
  iterDau.rewind();
  if(trueTagPi>=0){
    trueTagYPi -= trueTagPi; trueTagYK -= trueTagK; trueTagYKs -= trueTagKs; trueTagYPi0 -= trueTagPi0;}
  if(_Verbose.value()) { 
    cout<<"YPi: "<<trueTagYPi<<", YK: "<<trueTagYK<<", YKs: "<<trueTagYKs<<", Else: "<<trueTagElse<<endl;
    cout<<"Ds: "<<trueTagDs<<", D: "<<trueTagD<<", Pi: "<<trueTagPi<<", K: "<<trueTagK<<endl<<endl;
  }
  if (nSL != 1) { // Only continues if it is a semileptonic decay
    MCType = MCSubmode = MCPions = MCBrem = MCScatter = MCD = 0;
    truePLep = truePD = 0.0;
    uidBTag = uidD = uidLep = -1;
    return;
  }
  iterDau = theB -> daughterIterator();
  if(_Verbose.value())  cout<<"MCTagger manuelf: Loop over the daughters of the signal B"<<endl;
  while (aDaughter = iterDau()) {
    // If it is one the D(*)'s, set theD, DType, MCSubmode and uidD
    if ((int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::D0 || 
	(int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::anti_D0) {
      if (!DType) {
	theD = aDaughter;
	DType = 1;
	MCSubmode = D0submode(aDaughter);
	uidD = (int)aDaughter->uid();
      }
      else {
	DType = -1;
	MCSubmode = 0;
	uidD = -1;
      }
    }
    if ((int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::D_star0 || 
	(int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::anti_D_star0) {
      if (!DType) {
	theD = aDaughter;
	DType = 2;
	MCSubmode = Ds0submode(aDaughter);
	uidD = (int)aDaughter->uid();
      }
      else {
	DType = -1;
	MCSubmode = 0;
	uidD = -1;
      }
    }
    if ((int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::D_plus || 
	(int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::D_minus) {
      if (!DType) {
	theD = aDaughter;
	DType = 3;
	MCSubmode = Dsubmode(aDaughter);
	uidD = (int)aDaughter->uid();
      }
      else {
	DType = -1;
	MCSubmode = 0;
	uidD = -1;
      }
    }
    if ((int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::D_star_plus || 
	(int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::D_star_minus) {
      if (!DType) {
	theD = aDaughter;
	DType = 4;
	MCSubmode = Dssubmode(aDaughter);
	uidD = (int)aDaughter->uid();
      }
      else {
	DType = -1;
	MCSubmode = 0;
	uidD = -1;
      }
    }
    // The same with the Dss. MCSubmode is the D(*) decay from Dss->D(*)pi
    int abspdt = (int)aDaughter->pdtEntry()->lundId();
    if (abspdt < 0) abspdt *= -1;
    if (abspdt == (int)PdtLund::D_0_star_plus ||
	abspdt == (int)PdtLund::D_0_star0 ||
	abspdt == (int)PdtLund::D_1_plus ||
	abspdt == (int)PdtLund::D_10 ||
	abspdt == (int)PdtLund::D_2_star_plus ||
	abspdt == (int)PdtLund::D_2_star0 ||
	abspdt == (int)PdtLund::D_prime_1_plus ||
	abspdt == (int)PdtLund::D_prime_10 ||
	abspdt == (int)PdtLund::D_2S_plus ||
	abspdt == (int)PdtLund::D_2S0 ||
	abspdt == (int)PdtLund::D_star_2S_plus ||
	abspdt == (int)PdtLund::D_star_2S0) {
      if (!DType) {
	theD = aDaughter;
	DType = 5;
	MCSubmode = Dstarstarsubmode(aDaughter);
      }
      else {
	DType = -1;
	MCSubmode = 0;
      }
    }

    if ((int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::e_minus || 
	(int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::e_plus) {
      if (!lepType) {
	theLepton = aDaughter;
	lepType = 1;
	trueLepCharge = (aDaughter->charge() > 0) ? 1 : -1;
      }
      else {
	lepType = -1;
      }
    }
    if ((int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::mu_minus || 
	(int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::mu_plus) {
      if (!lepType) {
	theLepton = aDaughter;
	lepType = 2;
	trueLepCharge = (aDaughter->charge() > 0) ? 1 : -1;
      }
      else {
	lepType = -1;
      }
    }
    if ((int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::tau_minus ||
	(int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::tau_plus) {
      if (!lepType) {
	theLepton = aDaughter;
	lepType = 3;
	trueLepCharge = (aDaughter->charge() > 0) ? 1 : -1;
      }
      else lepType = -1;
    }
    if (abs((int)aDaughter->pdtEntry()->lundId()) == (int)PdtLund::nu_e ||
	abs((int)aDaughter->pdtEntry()->lundId()) == (int)PdtLund::nu_mu) {
      uidNu = (int)aDaughter->uid();
    }
    if ((int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::pi0 ||
	(int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::pi_minus ||
	(int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::pi_plus) {
      pion = 1;
      if ((int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::pi0) {
	MCPions += 10;
	uidDssPi0 = (int)aDaughter->uid();
      }
      else MCPions += 1;
    }
    if (abs((int)aDaughter->pdtEntry()->lundId()) == (int)PdtLund::eta) {
      MCPions += 100;
    }
    if ((int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::gamma)
      MCBrem = 1;
  }// End loop over daughter of the signal B

  // If I found multiple D's or leptons
  if (DType < 0 || lepType < 0) {
    MCType = MCSubmode = MCBrem = MCScatter = MCD = 0;
    return;
  }

  MCD = (int)theD->pdtEntry()->lundId();

  // Setting the pointer theFinalLepton to the lepton that we'd detect. It used to crash with tau->pi
  theFinalLepton = 0;
  if (abs((int)theLepton->pdtEntry()->lundId()) == (int)PdtLund::tau_minus) {
    MCTaumode = 0;
    HepAListIterator<BtaCandidate> iterDauTau = theLepton->daughterIterator();
    BtaCandidate *aLepDau;
    if(_Verbose.value()) cout<<"MCTagger manuelf: While over Tau daughters"<<endl;
    while (aLepDau = iterDauTau()) {
      if(_Verbose.value())  cout<<"MCTagger manuelf: The lepton daughter is "<<aLepDau->pdtEntry()->lundId()<<endl;
      if (abs((int)aLepDau->pdtEntry()->lundId()) == (int)PdtLund::e_minus) {
	theFinalLepton = aLepDau;
	MCTaumode = 1;
      }
      if (abs((int)aLepDau->pdtEntry()->lundId()) == (int)PdtLund::mu_minus) {
	theFinalLepton = aLepDau;
	MCTaumode = 2;
      }
      if ((int)aLepDau->pdtEntry()->lundId() == (int)PdtLund::gamma)
	MCBrem = 1;
    }
    trueTauFlight = theLepton->decayVtx()->point().distanceTo(theLepton->theMother()->decayVtx()->point());
    HepLorentzVector pcmtau = theLepton->p4();
    pcmtau.boost(-(theLepton->theMother()->p4().boostVector()));
    trueTauE = pcmtau.e();
    if(theFinalLepton){ 
      HepLorentzVector pmiss = theB->p4() - theD->p4() - theFinalLepton->p4();
      trueMM2 = pmiss.mag2();
    }
  } //tau?
  else theFinalLepton = theLepton;

  if(_Verbose.value())  cout<<"MCTagger manuelf: Storing the true properties of the final lepton"<<endl;
  if (theFinalLepton) {
    if (theFinalLepton->nDaughters() > 0) MCScatter = 1;
    truePLep = theFinalLepton->p();
    uidLep = (int)theFinalLepton->uid();
    HepLorentzVector pcmlep = theFinalLepton->p4();
    pcmlep.boost(-(theB->p4().boostVector()));
    trueLepE = pcmlep.e();
  }
 if(_Verbose.value())  cout<<"MCTagger manuelf: Storing D kinematics (one angle and q^2)"<<endl;
  truePD = theD->p();
  trueDmass = theD->mass();
  if ((DType==1 || DType==3) && (lepType==1 || lepType==2 || lepType==3) && !pion) {
    calcKinVars(theB->p4(),theLepton->p4(),theD->p4());
    // If it is a tau we also store the helicity angle for the final lepton
    if (lepType==3 && theFinalLepton) {
      XSLKin k(theB->p4(),theFinalLepton->p4(),theD->p4());
      trueCTL2 = k.theta_l();
    }
  } //Dlnu FF calc
 if(_Verbose.value())  cout<<"MCTagger manuelf: Storing D* kinematics (3 angles and q^2)"<<endl;
  if ((DType==2 || DType==4 || DType==5) && (lepType==1 || lepType==2 || lepType==3) && !pion) {
    HepAListIterator<BtaCandidate> iterDauDstar = theD->daughterIterator();
    BtaCandidate* aDstarDau;
    BtaCandidate *theDstarDau(0);
    while (aDstarDau = iterDauDstar()) {
      int dauLund = abs((int)aDstarDau->pdtEntry()->lundId());
      if (dauLund == (int)PdtLund::D0 || dauLund == (int)PdtLund::D_plus ||
	  dauLund == (int)PdtLund::D_star_plus || dauLund == (int)PdtLund::D_star0) theDstarDau = aDstarDau;
    }
    iterDauDstar.rewind();
    calcKinVars(theB->p4(),theLepton->p4(),theD->p4(),theDstarDau->p4());
    if (lepType==3 && theFinalLepton) {
      XSLKin k(theB->p4(),theFinalLepton->p4(),theD->p4(),theDstarDau->p4());
      trueCTL2 = k.ctl();
    }
  } //D*lnu FF calc
  if(_Verbose.value()) cout<<"MCTagger manuelf: Kinematics calculated"<<endl;

  //if i found a pion
  if (pion) {
    MCType = 13;
    MCDssmode = DType;
    if (MCPions == 1) MCDssmode += 4;
    return;
  }

  if (DType > 0 && lepType > 0) {
    if (DType == 1 && lepType == 1) MCType = 1;
    if (DType == 2 && lepType == 1) MCType = 2;
    if (DType == 1 && lepType == 2) MCType = 3;
    if (DType == 2 && lepType == 2) MCType = 4;
    if (DType == 1 && lepType == 3) MCType = 5;
    if (DType == 2 && lepType == 3) MCType = 6;
    if (DType == 3 && lepType == 1) MCType = 7;
    if (DType == 4 && lepType == 1) MCType = 8;
    if (DType == 3 && lepType == 2) MCType = 9;
    if (DType == 4 && lepType == 2) MCType = 10;
    if (DType == 3 && lepType == 3) MCType = 11;
    if (DType == 4 && lepType == 3) MCType = 12;
    if (DType == 5)                 MCType = 14;
    if (MCPions >= 100)             MCType = 15;
    }

  if (!MCType) cout << "UNKNOWN MCType! DType = " << DType <<
		 " and MCD = " << MCD <<
		 " and lepType = " << lepType << endl;

}

// 1 if the B has a D and a lepton, 0 otherwise. Sets isVub to 1 if there is a lepton
// and there is not a D 
bool DonutMCTagger::isSemiLep(BtaCandidate* aB)
{
  HepAListIterator<BtaCandidate> iterDau = aB -> daughterIterator();
  BtaCandidate* aDau;
  int lid;
  int charm = 0, lep = 0;
  while (aDau = iterDau()) {
    lid = (int)aDau -> pdtEntry() -> lundId();
    if (lid < 0) lid *= -1;
    if (lid == (int)PdtLund::e_minus || lid == (int)PdtLund::mu_minus || lid == (int)PdtLund::tau_minus) {
      lep = 1;
    }
    if (lid == (int)PdtLund::D_plus ||
	lid == (int)PdtLund::D0 ||
	lid == (int)PdtLund::D_star_plus ||
	lid == (int)PdtLund::D_star0 ||
	lid == (int)PdtLund::D_0_star_plus ||
	lid == (int)PdtLund::D_0_star0 ||
	lid == (int)PdtLund::D_1_plus ||
	lid == (int)PdtLund::D_10 ||
	lid == (int)PdtLund::D_2_star_plus ||
	lid == (int)PdtLund::D_2_star0 ||
	lid == (int)PdtLund::D_prime_1_plus ||
	lid == (int)PdtLund::D_prime_10 ||
	lid == (int)PdtLund::D_2S_plus ||
	lid == (int)PdtLund::D_2S0 ||
	lid == (int)PdtLund::D_star_2S_plus ||
	lid == (int)PdtLund::D_star_2S0) {
      charm = 1;
    }
  }
  iterDau.rewind();
  if (charm && lep) return true;
  if (lep) isVub = 1;
  return false;
}

int DonutMCTagger::Dsubmode(BtaCandidate* aD)
{
  int kcount = 0, kscount = 0, picount = 0, pi0count = 0, ecount = 0, mucount = 0;

  bool ret = count(aD,kcount,kscount,picount,pi0count,ecount,mucount);

  MCUnsim = 1;
  if (!ret) return 0;
  //something other than K/pi/l in tree

  MCUnsim = 0;
  if (kcount == 1 && picount == 2 && kscount == 0 && pi0count == 0) return 1;
  if (kcount == 1 && picount == 2 && kscount == 0 && pi0count == 1) return 2;
  if (kcount == 0 && picount == 1 && kscount == 1 && pi0count == 0) return 3;
  if (kcount == 0 && picount == 1 && kscount == 1 && pi0count == 1) return 4;
  if (kcount == 0 && picount == 3 && kscount == 1 && pi0count == 0) return 5;
  MCUnsim = 1;
  if (kcount == 2 && picount == 1 && kscount == 0 && pi0count == 0) return 6;
  if (kcount == 1 && picount == 0 && kscount == 1 && pi0count == 0) return 7;
  return 0;
}
// Mode of D0 decay. MCUnnsim gets set to 1 if it is unknown
int DonutMCTagger::D0submode(BtaCandidate* aD)
{
  int kcount = 0, kscount = 0, picount = 0, pi0count = 0, ecount = 0, mucount = 0;

  bool ret = count(aD,kcount,kscount,picount,pi0count,ecount,mucount);

  MCUnsim = 1;
  if (!ret) return 0;
  //something other than K/pi in tree

  MCUnsim = 0;
  if (kcount == 1 && picount == 1 && kscount == 0 && pi0count == 0) return 1;
  if (kcount == 1 && picount == 1 && kscount == 0 && pi0count == 1) return 2;
  if (kcount == 1 && picount == 3 && kscount == 0 && pi0count == 0) return 3;
  if (kcount == 0 && picount == 2 && kscount == 1 && pi0count == 0) return 4;
  if (kcount == 0 && picount == 2 && kscount == 1 && pi0count == 1) return 5;
  MCUnsim = 1;
  if (kcount == 0 && picount == 0 && kscount == 2 && pi0count == 0) return 6;
  if (kcount == 0 && picount == 0 && kscount == 1 && pi0count == 1) return 7;
  if (kcount == 0 && picount == 2 && kscount == 0 && pi0count == 0) return 8;
  if (kcount == 2 && picount == 0 && kscount == 0 && pi0count == 0) return 9;
  return 0;
}
// Counts number of K, Ks, pi and pi0 among the daughters of aCand
// ret is false if there is something that is not K/pi/l
bool DonutMCTagger::count(BtaCandidate* aCand, int &kcount, int &kscount,
			  int &picount, int &pi0count, int &ecount, int &mucount){
  int eLund = (int)PdtLund::e_minus;
  int muLund = (int)PdtLund::mu_minus;
  int piLund = (int)PdtLund::pi_plus;
  int pi0Lund = (int)PdtLund::pi0;
  int KLund = (int)PdtLund::K_plus;
  int K0Lund = (int)PdtLund::K_S0;
  int aLund = abs((int)aCand->pdtEntry()->lundId());

  if (aLund == KLund) {
    kcount ++;
    return true;
  }
  if (aLund == K0Lund) {
    kscount ++;
    return true;
  }
  if (aLund == piLund) {
    picount ++;
    return true;
  }
  if (aLund == pi0Lund) {
    pi0count ++;
    return true;
  }
  if (aLund == eLund) {
    ecount ++;
    return true;
  }
  if (aLund == muLund) {
    mucount ++;
    return true;
  }

  if (aLund == (int)PdtLund::gamma) return true;
  //don't worry about FSR for mode ID

  if (aCand -> nDaughters() == 0) return false;
  //stable final-state particle, not K/pi/l

  bool ret = true;
  BtaCandidate* aDau;
  HepAListIterator<BtaCandidate> iterDau = aCand -> daughterIterator();
  while (aDau = iterDau()) {
    if (!count(aDau,kcount,kscount,picount,pi0count,ecount,mucount)) ret = false;
  }
  iterDau.rewind();
  return ret;
}


int DonutMCTagger::Dssubmode(BtaCandidate* aD)
{
  int dsmode=-1,dmode=-1,pi0;
  pi0 = 0;

  BtaCandidate* aDaughter;
  HepAListIterator<BtaCandidate> iterDau = aD -> daughterIterator();
  while (aDaughter = iterDau()) {
    if ((int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::D0 ||
	(int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::anti_D0) {
      dsmode = 1;
      dmode = D0submode(aDaughter);
    }
    if ((int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::D_plus ||
	(int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::D_minus) {
      dsmode = 2;
      dmode = Dsubmode(aDaughter);
    }
    if ((int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::pi0) {
      pi0 = 1;
      uidPi0 = (int)aDaughter->uid();
    }
  }
  iterDau.rewind();

  if (dsmode == 2 && pi0 == 0) {dsmode = 3; MCUnsim = 1;}
  return (100*dsmode) + dmode;
}

int DonutMCTagger::Ds0submode(BtaCandidate* aD)
{
  int dsmode=-1,dmode=-1,pi0;
  pi0 = 0;

  BtaCandidate* aDaughter;
  HepAListIterator<BtaCandidate> iterDau = aD -> daughterIterator();
  while (aDaughter = iterDau()) {
    if ((int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::D0 ||
	(int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::anti_D0) {
      dmode = D0submode(aDaughter);
    }
    if ((int)aDaughter->pdtEntry()->lundId() == (int)PdtLund::pi0) {
      pi0 = 1;
      uidPi0 = (int)aDaughter->uid();
    }
  }
  iterDau.rewind();

  dsmode = 1;
  if (pi0 == 0) dsmode = 2;
  return (100*dsmode) + dmode;
}

int DonutMCTagger::Dstarstarsubmode(BtaCandidate* aD)
{
  int submode=0;
  int pdt;

  BtaCandidate *aDaughter;
  HepAListIterator<BtaCandidate> iterDau = aD -> daughterIterator();
  while (aDaughter = iterDau()) {
    pdt = (int)aDaughter -> pdtEntry() -> lundId();
    if (pdt < 0) pdt *= -1;
    if (pdt == (int)PdtLund::D0)          {MCDssmode = 1; submode = D0submode (aDaughter);}
    if (pdt == (int)PdtLund::D_star0)     {MCDssmode = 2; submode = Ds0submode(aDaughter);}
    if (pdt == (int)PdtLund::D_plus)      {MCDssmode = 3; submode = Dsubmode  (aDaughter);}
    if (pdt == (int)PdtLund::D_star_plus) {MCDssmode = 4; submode = Dssubmode (aDaughter);}
    if (pdt == (int)PdtLund::pi_plus)     MCPions += 1;
    if (pdt == (int)PdtLund::pi0)         MCPions += 10;
    if (pdt == (int)PdtLund::eta)         MCPions += 100;
    if (pdt == (int)PdtLund::pi0)         uidDssPi0 = (int)aDaughter->uid();
  }
  if (MCPions==1) MCDssmode += 4;
  return submode;
}

// Returns MCTwoBodyMode, a number whose last 2 digits tell which kind of particle
// Last digit is: 1 pi, 2 rho, 3 K+, 4 K*+ or the same format as second to last: 1 D, 2 D*, 3 D**, 4 Ds, 5 Ds*, 
// 6 Ds** if bigger than 100. If >1000, 2 two body B's 
int DonutMCTagger::TwoBodyMode(BtaCandidate* aB)
{
  int dmode = 0, dmode2 = 0, hmode = 0, lid;
  BtaCandidate *aDau;
  HepAListIterator<BtaCandidate> iterDau = aB -> daughterIterator();
  while (aDau = iterDau()) {
    lid = (int)aDau -> pdtEntry() -> lundId();
    if (lid < 0) lid *= -1;
    if (lid == (int)PdtLund::D0 || 
	lid == (int)PdtLund::D_plus)
      if (dmode) dmode2 = 1;
      else dmode = 1;
    if (lid == (int)PdtLund::D_star0 ||
	lid == (int)PdtLund::D_star_plus)
      if (dmode) dmode2 = 2;
      else dmode = 2;
    if (lid == (int)PdtLund::D_0_star_plus ||
	lid == (int)PdtLund::D_0_star0 ||
	lid == (int)PdtLund::D_1_plus ||
	lid == (int)PdtLund::D_10 ||
	lid == (int)PdtLund::D_2_star_plus ||
	lid == (int)PdtLund::D_2_star0 ||
	lid == (int)PdtLund::D_prime_1_plus ||
	lid == (int)PdtLund::D_prime_10 ||
	lid == (int)PdtLund::D_2S_plus ||
	lid == (int)PdtLund::D_2S0 ||
	lid == (int)PdtLund::D_star_2S_plus ||
	lid == (int)PdtLund::D_star_2S0) 
      if (dmode) dmode2 = 3;
      else dmode = 3;
    if (lid == (int)PdtLund::D_s_plus)
      if (dmode) dmode2 = 4;
      else dmode = 4;
    if (lid == (int)PdtLund::D_s_star_plus)
      if (dmode) dmode2 = 5;
      else dmode = 5;
    if (lid == (int)PdtLund::D_s0_star_plus ||
	lid == (int)PdtLund::D_s1_plus ||
	lid == (int)PdtLund::D_s2_star_plus)
      if (dmode) dmode2 = 6;
      else dmode = 6;

    if (lid == (int)PdtLund::pi_plus) hmode = 1;
    if (lid == (int)PdtLund::rho_plus) hmode = 2;
    if (lid == (int)PdtLund::K_plus) hmode = 3;
    if (lid == (int)PdtLund::K_star_plus) hmode = 4;
  }
  iterDau.rewind();
  if (dmode > 0 && hmode > 0) return (10 * dmode) + hmode;
  if (dmode2 > 0) return 100 + 10*dmode + dmode2;
  return 0;
}

// Calculates the true angles and q^2 for a D*
void DonutMCTagger::calcKinVars(HepLorentzVector pB, HepLorentzVector pLep, HepLorentzVector pDstar, HepLorentzVector pD)
{
  XSLKin k(pB,pLep,pDstar,pD);
  trueCTL = k.ctl();
  trueCTV = k.ctv();
  trueChi = k.chi();
  trueQ2 = k.q2();
}

// Calculates the true angle and q^2 for a D
void DonutMCTagger::calcKinVars(HepLorentzVector pB, HepLorentzVector pLep, HepLorentzVector pD)
{
  XSLKin k(pB,pLep,pD);
  trueCTL = k.theta_l();
  trueQ2 = k.q2();
  //storing theta_l in a varible normally used for cos(theta_l)
  //there will never be conflicts (distinct MCTypes); just know that interpretation changes
}
