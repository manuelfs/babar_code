#ifndef DONUTEVETO_HH
#define DONUTEVETO_HH

#include "Framework/AppModule.hh"

#include "AbsParm/AbsParmIfdStrKey.hh"
#include "Beta/BtaCandidate.hh"
#include "CLHEP/Alist/AList.h"
#include "ErrLogger/ErrLog.hh"
//#include "FrameUtil/FwkString.hh"
#include "ProxyDict/IfdKey.hh"

class HepHistogram;
class BtaPrintTree;

class DonutEVeto : public AppModule {

//--------------------
// Instance Members --
//--------------------

public:

  // Constructors
  DonutEVeto( const char* const theName, 
		  const char* const theDescription );
  
  // Destructor
  virtual ~DonutEVeto( );
  
  // Operations
  
  virtual AppResult           beginJob( AbsEvent* anEvent );
  virtual AppResult           event   ( AbsEvent* anEvent );
  virtual AppResult           endJob  ( AbsEvent* anEvent );

  // input lists
  AbsParmIfdStrKey _rawList;
  AbsParmIfdStrKey _goodList;
  AbsParmIfdStrKey _rawTrackList;
  AbsParmIfdStrKey _goodTrackList;

  // the lists themselves
  HepAList<BtaCandidate>* rawList;
  HepAList<BtaCandidate>* goodList;
  HepAList<BtaCandidate>* rawTrackList;
  HepAList<BtaCandidate>* goodTrackList;

  //global stuff
  HepHistogram* eVetoHist;
  HepHistogram* nTrkSel;
  HepHistogram* nTrkVeto;
  HepHistogram* nGamSel;
  HepHistogram* nGamVeto;
  HepHistogram* qTot;
  HepHistogram* eTot;
  HepHistogram* qTotVeto;
};

#endif
