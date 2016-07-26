#ifndef DONUTMCTAGGER_HH
#define DONUTMCTAGGER_HH

#include "Framework/AppModule.hh"

#include "AbsParm/AbsParmIfdStrKey.hh"
#include "Beta/BtaCandidate.hh"
#include "CLHEP/Alist/AList.h"
#include "ErrLogger/ErrLog.hh"
//#include "FrameUtil/FwkString.hh"
#include "ProxyDict/IfdKey.hh"

class HepHistogram;
class BtaPrintTree;

class DonutMCTagger : public AppModule {

//--------------------
// Instance Members --
//--------------------

public:

  // Constructors
  DonutMCTagger( const char* const theName, 
		  const char* const theDescription );
  
  // Destructor
  virtual ~DonutMCTagger( );
  
  // Operations
  
  virtual AppResult           beginJob( AbsEvent* anEvent );
  virtual AppResult           event   ( AbsEvent* anEvent );
  virtual AppResult           endJob  ( AbsEvent* anEvent );

  // input lists
  AbsParmIfdStrKey _truthList;
  AbsParmBool      _Verbose;

  // the lists themselves
  HepAList<BtaCandidate>* truthList;

  //global stuff
  HepHistogram* mcTypeHisto;

  //results
  BtaCandidate *theD, *theLepton, *theFinalLepton;
  int MCType,MCSubmode,MCCombmode,MCPions,MCTaumode,MCBrem,MCScatter,
    MCD,MC2Body,nTrueP,nTrueKL,nTrueNu,isVub,isBzero,MCDoubleSL,MCDssmode,MCUnsim;
  int uidBTag,uidD,uidLep,uidNu,uidPi0,uidDssPi0,trueLepCharge;
  float truePLep,truePD,trueDmass,trueCTL,trueCTV,trueChi,trueQ2,trueTauFlight,trueTauE,
    trueLepE,trueCTL2,trueMM2;
  //results manuelf
  int ntrueLep, trueMotherLep1, trueMotherLep2, trueMotherLep3, trueLundLep1,
    trueLundLep2, trueLundLep3, nMufromPi;
  float truePLep1, truePLep2, truePLep3;
  // manuelf: tag variables
  int trueTagPi, trueTagK, trueTagKs, trueTagPi0, trueTagE, trueTagMu, trueTagDs, trueTagD;
  int trueTagYPi, trueTagYK, trueTagYKs, trueTagYPi0;
  int trueTagElse;

  //functions
  void tagMC(AbsEvent* anEvent);
  int D0submode (BtaCandidate* aD);
  int Dsubmode  (BtaCandidate* aD);
  int Ds0submode(BtaCandidate* aD);
  int Dssubmode (BtaCandidate* aD);
  int Dstarstarsubmode(BtaCandidate* aD);
  bool count(BtaCandidate* aCand, int &kcount, int &kscount, int &picount, int &pi0count, 
	     int &ecount, int &mucount);
  bool isSemiLep(BtaCandidate* aB);
  int TwoBodyMode(BtaCandidate *aB);
  void calcKinVars(HepLorentzVector pB, HepLorentzVector pLep, HepLorentzVector pDstar, HepLorentzVector pD);
  void calcKinVars(HepLorentzVector pB, HepLorentzVector pLep, HepLorentzVector pD);

private:
  BtaPrintTree *printTree;
};


#endif
