#ifndef DONUTANALYSIS_HH
#define DONUTANALYSIS_HH

#include "AbsParm/AbsParmIfdStrKey.hh"
#include "BdbTime/BdbTime.hh"
#include "BetaCoreTools/BtaMcAssoc.hh"
#include "DonutUtils/weightManager.hh"
#include "Framework/AppModule.hh"
#include "Framework/AbsParmBool.hh"
#include "Framework/AbsParmDouble.hh"
#include "Framework/AbsParmGeneral.hh"
#include "UsrData/UsrData.hh"
#include "PDT/PdtLund.hh"
#include "TMVA/Reader.h"

template <class T> class HepAList;
template <class T> class HepAListIterator;
class BtaCandidate;
class BtaPrintTree;
class HepTuple;
class HepHistogram;
class EventInfo;
class AbsEventID;
class BdbTime;

class DonutAnalysis : public AppModule 
{
public:
  DonutAnalysis( const char* const, const char* const );
  virtual ~DonutAnalysis( );
  virtual AppResult beginJob( AbsEvent* );
  virtual AppResult event( AbsEvent* );
  virtual AppResult endJob  ( AbsEvent* );    
  // manuelf
  int _GroupEmpty[90][59];
  double _TruePurity[90][59];
private:
  AbsParmIfdStrKey keyUpsilonList;
  AbsParmIfdStrKey keyTrackList;
  AbsParmIfdStrKey keyPhotonList;
  AbsParmIfdStrKey keyelecList;
  AbsParmIfdStrKey keymuonList;
  AbsParmIfdStrKey keykaonList;
  AbsParmIfdStrKey keypionList;
  AbsParmIfdStrKey keyGTLList;
  AbsParmIfdStrKey keyGTVLList;
  AbsParmIfdStrKey keyPi0List;
  AbsParmIfdStrKey keySoftPi0List;
  AbsParmIfdStrKey keyKsList;
  AbsParmIfdStrKey keyMCTruthList;
  AbsParmIfdStrKey keyeventInfoList;
  AbsParmIfdStrKey keyeventIDList;
  AbsParmIfdStrKey keycandData;
  AbsParmBool      _doVertex;
  AbsParmBool      _doMass;
  AbsParmBool      _doSideband;
  AbsParmBool      _deltaESideband;
  AbsParmBool      _skipBPlus;
  AbsParmBool      _Verbose;
  AbsParmGeneral<string> _EmptyFile;
  AbsParmBool      _doBestBEEx;
  AbsParmDouble    _deltaEcut;
  AbsParmDouble    _mEScut;
  AbsParmGeneral<string> _StemFolder;
  AbsParmGeneral<int>   _maxExtraTracks;

  HepAList<BtaCandidate> *UpsilonList;
  HepAList<BtaCandidate> *TrackList;
  HepAList<BtaCandidate> *PhotonList;
  HepAList<BtaCandidate> *AllTrackList;
  HepAList<BtaCandidate> *AllPhotonList;
  HepAList<BtaCandidate> *elecList;
  HepAList<BtaCandidate> *muonList;
  HepAList<BtaCandidate> *kaonList;
  HepAList<BtaCandidate> *pionList;
  HepAList<BtaCandidate> *GTLList;
  HepAList<BtaCandidate> *GTVLList;
  HepAList<BtaCandidate> *Pi0List;
  HepAList<BtaCandidate> *SoftPi0List;
  HepAList<BtaCandidate> *KsList;
  HepAList<BtaCandidate> *MCTruthList;

  HepAListIterator<BtaCandidate> *iterElec;
  HepAListIterator<BtaCandidate> *iterMuon;
  HepAListIterator<BtaCandidate> *iterKaon;
  HepAListIterator<BtaCandidate> *iterPion;
  HepAListIterator<BtaCandidate> *iterTracks;
  HepAListIterator<BtaCandidate> *iterAllTracks;
  HepAListIterator<BtaCandidate> *iterAllPhotons;
  HepAListIterator<BtaCandidate> *iterPhotons;
  HepAListIterator<BtaCandidate> *iterGTL;
  HepAListIterator<BtaCandidate> *iterGTVL;
  HepAListIterator<BtaCandidate> *iterPi0;
  HepAListIterator<BtaCandidate> *iterSoftPi0;
  HepAListIterator<BtaCandidate> *iterSoftPi02;
  HepAListIterator<BtaCandidate> *iterKs;
  HepAListIterator<BtaCandidate> *iterMCTruth;

  EventInfo *eventInfo;
  AbsEventID *_eventID;
  BtaPrintTree *printTree;
  BdbTime *run4Start;

  UsrVariable<int> _decayMode;
  UsrVariable<float> _purity,_intpurity,_mES,_deltaE;

  TString method; 
  TMVA::Reader * readBestB;
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
  Float_t impi0, ie1pi0, idmpi0, ieextrapi0, ippi0, imm2pi0, ipmisspi0;
  int tempRejectedTracks, tempRejectedPhotons, tempTagChargedMult;
  float tempBTagDmass, tempBTagDeltam;
  float tempEExtra, tempMES, tempDeltaE;
  float tempDmass, tempDeltam, tempCosT;
  bool mEScorrection, isCocktail;
  WeightManager *myWM, *myWMFF;

  HepTuple *cands;
  HepTuple *dmass;
  HepHistogram *effhist[4];
  HepHistogram *meshist[4],*dehist[4];
  HepHistogram *pidhists[4];

  float d0cutslow[10],d0cutshigh[10],dccutslow[10],dccutshigh[10];
  float ds0d0cutslow[10],ds0d0cutshigh[10],dscd0cutslow[10],dscd0cutshigh[10],dscdccutslow[10],dscdccutshigh[10];
  float ds0cutslow[2],ds0cutshigh[2],dsccutslow[2],dsccutshigh[2];

  float emcfwd,emcbwd,trkfwd,trkbwd;

  bool isElectron(BtaCandidate* aTrk);
  bool isMuon(BtaCandidate* aTrk);
  bool isKaon(BtaCandidate* aTrk);
  bool isPion(BtaCandidate* aTrk);
  bool isGTL(BtaCandidate* aTrk);
  bool isGTVL(BtaCandidate* aTrk);
  int DType(BtaCandidate* theD);
  int D0Type(BtaCandidate* theD);
  bool checkD0(float mass, int mode);
  bool checkD(float mass, int mode);
  bool checkDs0(float dsmass,float deltam,int dmode,int dsmode,int kill);
  bool checkDs(float dsmass,float deltam,int dmode,int dsmode,int kill);
  bool sigcheckD0(float mass, int mode);
  bool sigcheckD(float mass, int mode);
  bool sigcheckDs0(float dsmass,float deltam,int dmode,int dsmode,int kill);
  bool sigcheckDs(float dsmass,float deltam,int dmode,int dsmode,int kill);
  void countMult(BtaCandidate* aCand,int &nchg, int &nneut);
  int combMode(BtaCandidate *aB, int &bmode, int &dsmode, int &MClongCombD);
  void nTruthMatched(BtaCandidate *X, int &ChargedNon, int &NeutralNon, int &nCharged,
				  int &nNeutral,BtaMcAssoc* truthMap, int momB);
  int motherB(BtaCandidate *X);
  int expo(int x, int n);
};

#endif
