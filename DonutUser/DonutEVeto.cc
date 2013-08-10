#include "BaBar/BaBar.hh"

#include "DonutUser/DonutEVeto.hh"

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

DonutEVeto::DonutEVeto( const char* const theName, 
			const char* const theDescription )
  : AppModule( theName, theDescription )
  , _rawList ("rawList"      , this, "CalorNeutral")
  , _goodList("goodList"     , this, "RecoilGoodPhotons")
  , _rawTrackList ("rawTrackList"      , this, "ChargedTracks")
  , _goodTrackList("goodTrackList"     , this, "RecoilGoodTracks")
{
  commands()->append( &_rawList );
  commands()->append( &_goodList );
  commands()->append( &_rawTrackList );
  commands()->append( &_goodTrackList );
}

DonutEVeto::~DonutEVeto( )
{
}

AppResult
DonutEVeto::endJob( AbsEvent* anEvent )
{
  return AppResult::OK;
}

AppResult
DonutEVeto::beginJob( AbsEvent* anEvent )
{
  HepTupleManager* manager = gblEnv->getGen()->ntupleManager();
  assert (manager!=0);
  eVetoHist = manager -> histogram("eVeto",100,0,1);
  nTrkSel = manager -> histogram("nTrkSel",30,0,30);
  nTrkVeto = manager -> histogram("nTrkVeto",30,0,30);
  nGamSel = manager -> histogram("nGamSel",50,0,50);
  nGamVeto = manager -> histogram("nGamVeto",50,0,50);
  qTot = manager -> histogram("qTot",20,-10,10);
  eTot = manager -> histogram("eTot",200,0,10);
  qTotVeto = manager -> histogram("qTotVeto",20,-10,10);

  return AppResult::OK;
}

AppResult
DonutEVeto::event( AbsEvent* anEvent )
{
  rawList  = Ifd< HepAList<BtaCandidate> >::get(anEvent, _rawList.value());
  goodList = Ifd< HepAList<BtaCandidate> >::get(anEvent, _goodList.value());
  rawTrackList  = Ifd< HepAList<BtaCandidate> >::get(anEvent, _rawTrackList.value());
  goodTrackList = Ifd< HepAList<BtaCandidate> >::get(anEvent, _goodTrackList.value());

  HepAListIterator<BtaCandidate> iterRaw (*rawList);
  HepAListIterator<BtaCandidate> iterGood(*goodList);
  HepAListIterator<BtaCandidate> iterRawTrack (*rawTrackList);
  HepAListIterator<BtaCandidate> iterGoodTrack(*goodTrackList);

  float etot = 0.0, eveto = 0.0;
  int qtot = 0, qveto = 0;
  BtaCandidate *aGam;

  //cout<<"EVeto manuelf: Before whiles"<<endl;
   while (aGam = iterRaw()) {
    eveto += aGam->energy();
  }
  iterRaw.rewind();
  while (aGam = iterGood()) {
    eveto -= aGam->energy();
    etot += aGam->energy();
  }
  iterGood.rewind();

  while (aGam = iterRawTrack()) {
    qveto += (int)aGam->charge();
  }
  iterRawTrack.rewind();
  while (aGam = iterGoodTrack()) {
    qveto -= (int)aGam->charge();
    qtot += (int)aGam->charge();
  }
  iterGoodTrack.rewind();

  eVetoHist -> accumulate(eveto);
  nTrkSel -> accumulate(goodTrackList -> length());
  nTrkVeto -> accumulate(rawTrackList -> length() - goodTrackList -> length());
  nGamSel -> accumulate(goodList -> length());
  nGamVeto -> accumulate(rawList -> length() - goodList -> length());
  qTot -> accumulate(qtot);
  eTot -> accumulate(etot);
  qTotVeto -> accumulate(qveto);
  //cout<<"EVeto manuelf: End event"<<endl;

  return AppResult::OK;
}

