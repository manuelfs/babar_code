//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DonutAnalysis.cc,v 1.6 2009/09/03 01:16:20 manuelf Exp $
//
// Description:
//      DonutAnalysis - Core code to build the ntuples for the 
//      B->D(*)TauNu analysis
//
// Author List:
//      Michael Mazur (original)                  UC Santa Barbara
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      09/08/19 manuelf -- Added MCComblong with the numbers of all particles in combinatoric
//      09/08/19 manuelf -- Added MCComblong with the numbers of all particles in combinatoric
//                          Added the mES correction and the weights
//      09/06/10 manuelf -- Added the 8 MVA for the D**
//      09/05/25 manuelf -- Took candRejectedTracks/Photons out of the MVA
//                          Cleaned up some comments and unused smear
//                          Added candM2Tru and uncommented some other variables
//      09/05/14 manuelf -- Added the 8 MVA that I intend to use
//                          candMES has now the mES correction
//                          Cleaned up the warnings
//      09/05/05 manuelf -- Initial committed code with some TMVA
//------------------------------------------------------------------------


#include "BaBar/BaBar.hh"
#include "DonutUser/DonutAnalysis.hh"

#include "AbsEnv/AbsEnv.hh"
#include "AbsEvent/AbsEvent.hh"
#include "AbsEvent/AbsEventID.hh"
#include "AbsEventTag/AbsEventTag.hh"
#include "Beta/BtaCandidate.hh"
#include "Beta/BtaConstraint.hh"
#include "Beta/EventInfo.hh"
#include "BetaCoreTools/BtaBooster.hh"
#include "BetaCoreTools/BtaBVariables.hh"
#include "BetaCoreTools/BtaOpMakeTree.hh"
#include "BetaCoreTools/BtaPrintTree.hh"
#include "BetaCoreTools/BtaThrust.hh"
#include "BetaCoreTools/BtaTreeNavigator.hh"
#include "BetaMC/BtaMcAssocGHit.hh"
#include "BetaPid/PidWeight.hh"
#include "BtaEnv/BtaEnv.hh"
#include "CLHEP/Alist/AIterator.h"
#include "CLHEP/Alist/AList.h"
#include "CLHEP/Random/RandEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/Random.h"
#include "DonutUtils/PidCorrectMesMean.hh"
#include "EffTable/EffTable.hh"
#include "EidData/EidCondKeyTriplet.hh"
#include "EidData/EidEventTriplet.hh"
#include "ErrLogger/ErrLog.hh"
#include "GenEnv/GenEnv.hh"
#include "HepTuple/HTValOrderedVector.h"
#include "HepTuple/Histogram.h"
#include "HepTuple/Tuple.h"
#include "HepTuple/TupleManager.h"
#include "PDT/Pdt.hh"
#include "PDT/PdtEntry.hh"
#include "ProbTools/probab.hh"
#include "ProxyDict/IfdKey.hh"
#include "ProxyDict/IfdStrKey.hh"
#include "RooFitCore/RooGlobalFunc.hh"
#include "VtxBase/VtxAbsAlgorithm.hh"
#include "VtxFitter/VtxFitterOper.hh"
#include <map>

using std::cout;
using std::endl;

DonutAnalysis::DonutAnalysis( const char* const name, 
		  const char* const description )
  : AppModule( name, description ),
    keyUpsilonList("UpsilonList",this,"SmpUpsilonDefault"),
    keyTrackList("TrackList",this,"ChargedTracks"),
    keyPhotonList("PhotonList",this,"CalorNeutral"),
    keyelecList("elecList",this,"PidLHElectrons"),
    keymuonList("muonList",this,"muMicroTight"),
    keykaonList("kaonList",this,"KLHTight"),
    keypionList("pionList",this,"piLHLoose"),
    keyGTLList("GTLList",this,"GoodTracksLoose"),
    keyGTVLList("GTVLList",this,"GoodTracksVeryLoose"),
    keyPi0List("Pi0List",this,"pi0AllDefault"),
    keySoftPi0List("SoftPi0List",this,"pi0VeryLoose"),
    keyKsList("KsList",this,"KsDefault"),
    keyMCTruthList("MCTruthList",this,"MCTruth"),
    keyeventInfoList("eventInfoList",this,"Default"),
    keyeventIDList("eventIDList",this,"AbsEventID"),
    keycandData("candData",this,"VcbCandData"),
    _doVertex("doVertex",this,false),
    _doMass("doMass",this,false),
    _doSideband("doSideband",this,false),
    _deltaESideband("deltaESideband",this,false),
    _skipBPlus("skipBPlus",this,false),
    _Verbose("Verbose",this,false), 
    _EmptyFile("EmptyFile", this, "" ),
    _doBestBEEx("doBestBEEx",this,true),
    _deltaEcut("deltaEcut",this,0.072),
    _mEScut("mEScut",this,5.27),
    _StemFolder("StemFolder",this,""),
    _maxExtraTracks("maxExtraTracks",this,0),
    _decayMode("decayMode"),
    _purity("purity"),
    _intpurity("intpurity"),
    _mES("mES"),
    _deltaE("deltaE")
{
  commands()->append( &keyUpsilonList ); 
  commands()->append( &keyTrackList ); 
  commands()->append( &keyPhotonList ); 
  commands()->append( &keyelecList ); 
  commands()->append( &keymuonList ); 
  commands()->append( &keykaonList ); 
  commands()->append( &keypionList ); 
  commands()->append( &keyGTLList ); 
  commands()->append( &keyGTVLList ); 
  commands()->append( &keyPi0List ); 
  commands()->append( &keySoftPi0List ); 
  commands()->append( &keyKsList ); 
  commands()->append( &keyMCTruthList ); 
  commands()->append( &keyeventInfoList ); 
  commands()->append( &keyeventIDList ); 
  commands()->append( &keycandData ); 
  commands()->append( &_doVertex ); 
  commands()->append( &_doMass ); 
  commands()->append( &_doSideband );
  commands()->append( &_deltaESideband );
  commands()->append( &_skipBPlus );
  commands()->append( &_Verbose );
  commands()->append( &_EmptyFile );
  commands()->append( &_doBestBEEx );
  commands()->append( &_deltaEcut );
  commands()->append( &_mEScut );
  commands()->append( &_StemFolder );
  commands()->append( &_maxExtraTracks );
}

DonutAnalysis::~DonutAnalysis( ){
}

AppResult DonutAnalysis::beginJob( AbsEvent* event )
{
  run4Start = new BdbTime(0);
  ErrMsg(routine)<< "begin job"<< endmsg;
  //1978 + cdb 25 year offset = 2003
  if (!(BdbTime::parseTime("01Sep1978","12:00:00",BdbTime::UTC,*run4Start))) {
    cout << "MAM: problem constructing a BdbTime!" << endl;
  }

  // manuelf: Initializing mES correction and weights
  TString folder = (TString)_StemFolder.value();
  mEScorrection = false;
  //if(folder.Contains("BSemiExclAdd-Run5-OnPeak")||folder.Contains("BSemiExclAdd-Run6-OnPeak")) mEScorrection = true;
  if(mEScorrection){
    cout<<"Doing the mES correction"<<endl;
    PidCorrectMesMean::initialize("../DonutUtils/mES_means_Runs123456_final.txt");
  } else cout<<"Not doing mES correction because folder is "<<folder<<endl;
  myWM  = new WeightManager("babar_code/Reweight/wFFBF.txt");
  int IsCocktail = 0;
  if(folder.Contains("SP-35") || folder.Contains("SP-36") || folder.Contains("SP-63")
      || folder.Contains("SP-55") || folder.Contains("SP-96")) IsCocktail = 1;
  myWMFF =  new WeightManager("babar_code/Reweight/wAllFF.txt",IsCocktail);

  // manuelf: Booking of the TMVA Reader
  cout<<"Booking the MVA readers"<<endl;
  method = "BDT";
  if (!_doBestBEEx.value()){ 
    string Mname = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaBestB_BDT.weights.txt";
    readBestB = new TMVA::Reader("!V");
    readBestB->AddVariable("candEExtra", &tempEExtra);
    readBestB->AddVariable("candMES", &tempMES);
    readBestB->AddVariable("candDmass", &tempDmass);
    readBestB->AddVariable("candDeltam", &tempDeltam);
    readBestB->AddVariable("candBTagDeltam", &tempBTagDeltam);
    readBestB->AddVariable("candBTagDmass", &tempBTagDmass);
    readBestB->AddVariable("candDeltaE", &tempDeltaE);
    readBestB->AddVariable("candRejectedTracks", &tempRejectedTracks);
    readBestB->BookMVA(method,Mname);
  }

  readCombD0 = new TMVA::Reader("V:Silent");
  readCombDs0 = new TMVA::Reader("!V:Silent");
  readCombDp = new TMVA::Reader("!V:Silent");
  readCombDsp = new TMVA::Reader("!V:Silent");
  readDlD0 = new TMVA::Reader("V:Silent");
  readDlDs0 = new TMVA::Reader("V:Silent");
  readDlDp = new TMVA::Reader("V:Silent");
  readDlDsp = new TMVA::Reader("V:Silent");

  readCombD0->AddVariable("candEExtra", &tempEExtra);
  readCombD0->AddVariable("candMES", &tempMES);
  readCombD0->AddVariable("candDmass", &tempDmass);
  readCombD0->AddVariable("candTagChargedMult", &tempTagChargedMult);
  readCombD0->AddVariable("candBTagDeltam", &tempBTagDeltam);
  readCombD0->AddVariable("candBTagDmass", &tempBTagDmass);
  readCombD0->AddVariable("candDeltaE", &tempDeltaE);
  readCombD0->AddVariable("candCosT", &tempCosT);

  TString file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaCombD0_BDT.weights.txt";
  readCombD0->BookMVA(method,file);

  readDlD0->AddVariable("candEExtra", &tempEExtra);
  readDlD0->AddVariable("candMES", &tempMES);
  readDlD0->AddVariable("candDmass", &tempDmass);
  readDlD0->AddVariable("candTagChargedMult", &tempTagChargedMult);
  readDlD0->AddVariable("candBTagDeltam", &tempBTagDeltam);
  readDlD0->AddVariable("candBTagDmass", &tempBTagDmass);
  readDlD0->AddVariable("candDeltaE", &tempDeltaE);
  readDlD0->AddVariable("candCosT", &tempCosT);

  file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDlD0_BDT.weights.txt";
  readDlD0->BookMVA(method,file);

  readCombDs0->AddVariable("candEExtra", &tempEExtra);
  readCombDs0->AddVariable("candMES", &tempMES);
  readCombDs0->AddVariable("candDmass", &tempDmass);
  readCombDs0->AddVariable("candDeltam", &tempDeltam);
  readCombDs0->AddVariable("candTagChargedMult", &tempTagChargedMult);
  readCombDs0->AddVariable("candBTagDeltam", &tempBTagDeltam);
  readCombDs0->AddVariable("candBTagDmass", &tempBTagDmass);
  readCombDs0->AddVariable("candDeltaE", &tempDeltaE);
  readCombDs0->AddVariable("candCosT", &tempCosT);

  file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaCombDs0_BDT.weights.txt";
  readCombDs0->BookMVA(method,file);

  readDlDs0->AddVariable("candEExtra", &tempEExtra);
  readDlDs0->AddVariable("candMES", &tempMES);
  readDlDs0->AddVariable("candDmass", &tempDmass);
  readDlDs0->AddVariable("candDeltam", &tempDeltam);
  readDlDs0->AddVariable("candTagChargedMult", &tempTagChargedMult);
  readDlDs0->AddVariable("candBTagDeltam", &tempBTagDeltam);
  readDlDs0->AddVariable("candBTagDmass", &tempBTagDmass);
  readDlDs0->AddVariable("candDeltaE", &tempDeltaE);
  readDlDs0->AddVariable("candCosT", &tempCosT);

  file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDlDs0_BDT.weights.txt";
  readDlDs0->BookMVA(method,file);

  readCombDp->AddVariable("candEExtra", &tempEExtra);
  readCombDp->AddVariable("candMES", &tempMES);
  readCombDp->AddVariable("candDmass", &tempDmass);
  readCombDp->AddVariable("candTagChargedMult", &tempTagChargedMult);
  readCombDp->AddVariable("candBTagDeltam", &tempBTagDeltam);
  readCombDp->AddVariable("candBTagDmass", &tempBTagDmass);
  readCombDp->AddVariable("candDeltaE", &tempDeltaE);
  readCombDp->AddVariable("candCosT", &tempCosT);

  file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaCombDp_BDT.weights.txt";
  readCombDp->BookMVA(method,file);

  readDlDp->AddVariable("candEExtra", &tempEExtra);
  readDlDp->AddVariable("candMES", &tempMES);
  readDlDp->AddVariable("candDmass", &tempDmass);
  readDlDp->AddVariable("candTagChargedMult", &tempTagChargedMult);
  readDlDp->AddVariable("candBTagDeltam", &tempBTagDeltam);
  readDlDp->AddVariable("candBTagDmass", &tempBTagDmass);
  readDlDp->AddVariable("candDeltaE", &tempDeltaE);
  readDlDp->AddVariable("candCosT", &tempCosT);

  file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDlDp_BDT.weights.txt";
  readDlDp->BookMVA(method,file);

  readCombDsp->AddVariable("candEExtra", &tempEExtra);
  readCombDsp->AddVariable("candMES", &tempMES);
  readCombDsp->AddVariable("candDmass", &tempDmass);
  readCombDsp->AddVariable("candDeltam", &tempDeltam);
  readCombDsp->AddVariable("candTagChargedMult", &tempTagChargedMult);
  readCombDsp->AddVariable("candBTagDeltam", &tempBTagDeltam);
  readCombDsp->AddVariable("candBTagDmass", &tempBTagDmass);
  readCombDsp->AddVariable("candDeltaE", &tempDeltaE);
  readCombDsp->AddVariable("candCosT", &tempCosT);

  file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaCombDsp_BDT.weights.txt";
  readCombDsp->BookMVA(method,file);

  readDlDsp->AddVariable("candEExtra", &tempEExtra);
  readDlDsp->AddVariable("candMES", &tempMES);
  readDlDsp->AddVariable("candDmass", &tempDmass);
  readDlDsp->AddVariable("candDeltam", &tempDeltam);
  readDlDsp->AddVariable("candTagChargedMult", &tempTagChargedMult);
  readDlDsp->AddVariable("candBTagDeltam", &tempBTagDeltam);
  readDlDsp->AddVariable("candBTagDmass", &tempBTagDmass);
  readDlDsp->AddVariable("candDeltaE", &tempDeltaE);
  readDlDsp->AddVariable("candCosT", &tempCosT);

  file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDlDsp_BDT.weights.txt";
  readDlDsp->BookMVA(method,file);

  //===================== D** MVAs ===============================================
   readDssCombD0 = new TMVA::Reader("V:Silent");
   readDssCombD0->AddVariable("mpi0", &impi0);
   readDssCombD0->AddVariable("candDmass", &tempDmass);
   readDssCombD0->AddVariable("dmpi0", &idmpi0);
   readDssCombD0->AddVariable("eextrapi0", &ieextrapi0);
   readDssCombD0->AddVariable("ppi0", &ippi0);
   readDssCombD0->AddVariable("e1pi0", &ie1pi0);
   readDssCombD0->AddVariable("candCosT", &tempCosT);
   readDssCombD0->AddVariable("candMES", &tempMES);
   readDssCombD0->AddVariable("candDeltaE", &tempDeltaE);
   file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDssCombD0_BDT.weights.txt";
   readDssCombD0->BookMVA(method,file);

   readDssDlD0 = new TMVA::Reader("V:Silent");
   readDssDlD0->AddVariable("mpi0", &impi0);
   readDssDlD0->AddVariable("candDmass", &tempDmass);
   readDssDlD0->AddVariable("dmpi0", &idmpi0);
   readDssDlD0->AddVariable("eextrapi0", &ieextrapi0);
   readDssDlD0->AddVariable("ppi0", &ippi0);
   readDssDlD0->AddVariable("e1pi0", &ie1pi0);
   readDssDlD0->AddVariable("candCosT", &tempCosT);
   readDssDlD0->AddVariable("candMES", &tempMES);
   readDssDlD0->AddVariable("candDeltaE", &tempDeltaE);
   file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDssDlD0_BDT.weights.txt";
   readDssDlD0->BookMVA(method,file);

   readDssCombDs0 = new TMVA::Reader("V:Silent");
   readDssCombDs0->AddVariable("mpi0", &impi0);
   readDssCombDs0->AddVariable("candDmass", &tempDmass);
   readDssCombDs0->AddVariable("dmpi0", &idmpi0);
   readDssCombDs0->AddVariable("eextrapi0", &ieextrapi0);
   readDssCombDs0->AddVariable("ppi0", &ippi0);
   readDssCombDs0->AddVariable("e1pi0", &ie1pi0);
   readDssCombDs0->AddVariable("candCosT", &tempCosT);
   readDssCombDs0->AddVariable("candDeltam", &tempDmass);
   readDssCombDs0->AddVariable("candMES", &tempMES);
   readDssCombDs0->AddVariable("candDeltaE", &tempDeltaE);
   file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDssCombDs0_BDT.weights.txt";
   readDssCombDs0->BookMVA(method,file);

   readDssDlDs0 = new TMVA::Reader("V:Silent");
   readDssDlDs0->AddVariable("mpi0", &impi0);
   readDssDlDs0->AddVariable("candDmass", &tempDmass);
   readDssDlDs0->AddVariable("dmpi0", &idmpi0);
   readDssDlDs0->AddVariable("eextrapi0", &ieextrapi0);
   readDssDlDs0->AddVariable("ppi0", &ippi0);
   readDssDlDs0->AddVariable("e1pi0", &ie1pi0);
   readDssDlDs0->AddVariable("candCosT", &tempCosT);
   readDssDlDs0->AddVariable("candDeltam", &tempDmass);
   readDssDlDs0->AddVariable("candMES", &tempMES);
   readDssDlDs0->AddVariable("candDeltaE", &tempDeltaE);
   file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDssDlDs0_BDT.weights.txt";
   readDssDlDs0->BookMVA(method,file);

   readDssCombDp = new TMVA::Reader("V:Silent");
   readDssCombDp->AddVariable("mpi0", &impi0);
   readDssCombDp->AddVariable("candDmass", &tempDmass);
   readDssCombDp->AddVariable("dmpi0", &idmpi0);
   readDssCombDp->AddVariable("eextrapi0", &ieextrapi0);
   readDssCombDp->AddVariable("ppi0", &ippi0);
   readDssCombDp->AddVariable("e1pi0", &ie1pi0);
   readDssCombDp->AddVariable("candCosT", &tempCosT);
   readDssCombDp->AddVariable("candMES", &tempMES);
   readDssCombDp->AddVariable("candDeltaE", &tempDeltaE);
   file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDssCombDp_BDT.weights.txt";
   readDssCombDp->BookMVA(method,file);

   readDssDlDp = new TMVA::Reader("V:Silent");
   readDssDlDp->AddVariable("mpi0", &impi0);
   readDssDlDp->AddVariable("candDmass", &tempDmass);
   readDssDlDp->AddVariable("dmpi0", &idmpi0);
   readDssDlDp->AddVariable("eextrapi0", &ieextrapi0);
   readDssDlDp->AddVariable("ppi0", &ippi0);
   readDssDlDp->AddVariable("e1pi0", &ie1pi0);
   readDssDlDp->AddVariable("candCosT", &tempCosT);
   readDssDlDp->AddVariable("candMES", &tempMES);
   readDssDlDp->AddVariable("candDeltaE", &tempDeltaE);
   file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDssDlDp_BDT.weights.txt";
   readDssDlDp->BookMVA(method,file);

   readDssCombDsp = new TMVA::Reader("V:Silent");
   readDssCombDsp->AddVariable("mpi0", &impi0);
   readDssCombDsp->AddVariable("candDmass", &tempDmass);
   readDssCombDsp->AddVariable("dmpi0", &idmpi0);
   readDssCombDsp->AddVariable("eextrapi0", &ieextrapi0);
   readDssCombDsp->AddVariable("ppi0", &ippi0);
   readDssCombDsp->AddVariable("e1pi0", &ie1pi0);
   readDssCombDsp->AddVariable("candCosT", &tempCosT);
   readDssCombDsp->AddVariable("candDeltam", &tempDmass);
   readDssCombDsp->AddVariable("candMES", &tempMES);
   readDssCombDsp->AddVariable("candDeltaE", &tempDeltaE);
   file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDssCombDsp_BDT.weights.txt";
   readDssCombDsp->BookMVA(method,file);

   readDssDlDsp = new TMVA::Reader("V:Silent");
   readDssDlDsp->AddVariable("mpi0", &impi0);
   readDssDlDsp->AddVariable("candDmass", &tempDmass);
   readDssDlDsp->AddVariable("dmpi0", &idmpi0);
   readDssDlDsp->AddVariable("eextrapi0", &ieextrapi0);
   readDssDlDsp->AddVariable("ppi0", &ippi0);
   readDssDlDsp->AddVariable("e1pi0", &ie1pi0);
   readDssDlDsp->AddVariable("candCosT", &tempCosT);
   readDssDlDsp->AddVariable("candDeltam", &tempDmass);
   readDssDlDsp->AddVariable("candMES", &tempMES);
   readDssDlDsp->AddVariable("candDeltaE", &tempDeltaE);
   file = "/u/br/manuelf/releases/tmva/TMVA/macros/weights/mvaDssDlDsp_BDT.weights.txt";
   readDssDlDsp->BookMVA(method,file);

  // manuelf: Selection of good modes
  string Ename = "../DonutUser/"; Ename+=_EmptyFile.value(); Ename += ".txt";
  ifstream Empty(Ename.c_str());
  int  mode;
  int SeedMode, YMode;
  for(int i=0; i<90; i++)
    for(int j=0; j<59; j++)
      _GroupEmpty[i][j] = 0;
  if(Empty){
    cout<<name()<<": Killing modes found in "<<Ename<<endl;
    while(!Empty.eof() && Empty){
      Empty >> mode;
      SeedMode = (mode/100)%100; YMode = mode%100;
      _GroupEmpty[SeedMode][YMode] = 1;
    }
  }else{
    cout<<name()<<": File "<<Ename<<" not found. Using all modes"<<endl;
  }
  fstream textfile;
  textfile.open("../DonutUser/BSemiExclAddPurities.txt",fstream::in);
  double Signal, Bkg, Pur,PurGroup;
  while (textfile){
    textfile >> YMode >> SeedMode >> Signal >> Bkg >> Pur >> PurGroup;
    _TruePurity[SeedMode][YMode] = Pur;
  }


  HepTupleManager * manager = gblEnv->getGen()->ntupleManager();
  cands = manager->ntuple("cands",1);
  dmass = manager->ntuple("dmass",2);
  effhist[0] = manager->histogram("effhist0",15,0,15,8,0,8);
  effhist[1] = manager->histogram("effhist1",15,0,15,8,0,8);
  effhist[2] = manager->histogram("effhist2",15,0,15,8,0,8);
  effhist[3] = manager->histogram("effhist3",15,0,15,8,0,8);

  meshist[0] = manager->histogram("meshist0",100,5.2,5.3);
  meshist[1] = manager->histogram("meshist1",100,5.2,5.3);
  meshist[2] = manager->histogram("meshist2",100,5.2,5.3);
  meshist[3] = manager->histogram("meshist3",100,5.2,5.3);
  dehist[0] = manager->histogram("dehist0",100,-.12,.12);
  dehist[1] = manager->histogram("dehist1",100,-.12,.12);
  dehist[2] = manager->histogram("dehist2",100,-.12,.12);
  dehist[3] = manager->histogram("dehist3",100,-.12,.12);

  pidhists[0] = manager->histogram("elechist",7,0,7);
  pidhists[1] = manager->histogram("muonhist",7,0,7);
  pidhists[2] = manager->histogram("kaonhist",7,0,7);
  pidhists[3] = manager->histogram("pionhist",7,0,7);

  printTree = new BtaPrintTree();
  printTree -> addTerminal(PdtLund::pi0);
//   printTree -> addTerminal(PdtLund::pi_plus);
//   printTree -> addTerminal(PdtLund::pi_minus);
  printTree -> addTerminal(PdtLund::K_plus);
  printTree -> addTerminal(PdtLund::K_minus);
  printTree -> addTerminal(PdtLund::mu_plus);
  printTree -> addTerminal(PdtLund::mu_minus);
  printTree -> setOutputFunc(BtaPrintTree::printName);

  float sigmad = 4.0;
  float sigmads = 4.0;
  d0cutslow[0]  = 1.865 - sigmad * 0.0065;
  d0cutshigh[0] = 1.865 + sigmad * 0.0065;
  d0cutslow[1]  = 1.863 - sigmad * 0.0123;
  d0cutshigh[1] = 1.863 + sigmad * 0.0123;
  d0cutslow[2]  = 1.865 - sigmad * 0.004;
  d0cutshigh[2] = 1.865 + sigmad * 0.004;
  d0cutslow[3]  = 1.865 - sigmad * 0.006;
  d0cutshigh[3] = 1.865 + sigmad * 0.006;
  d0cutslow[4]  = 1.864 - sigmad * 0.0096;
  d0cutshigh[4] = 1.864 + sigmad * 0.0096;

  d0cutslow[5]  = 1.864 - sigmad * 0.0096;
  d0cutshigh[5] = 1.864 + sigmad * 0.0096;
  d0cutslow[6]  = 1.863 - sigmad * 0.0123;
  d0cutshigh[6] = 1.863 + sigmad * 0.0123;
  d0cutslow[7]  = 1.865 - sigmad * 0.0065;
  d0cutshigh[7] = 1.865 + sigmad * 0.0065;
  d0cutslow[8]  = 1.865 - sigmad * 0.0065;
  d0cutshigh[8] = 1.865 + sigmad * 0.0065;

  dccutslow[0]  = 1.869 - sigmad * 0.0053;
  dccutshigh[0] = 1.869 + sigmad * 0.0053;
  dccutslow[1]  = 1.867 - sigmad * 0.01;
  dccutshigh[1] = 1.867 + sigmad * 0.01;
  dccutslow[2]  = 1.870 - sigmad * 0.005;
  dccutshigh[2] = 1.870 + sigmad * 0.005;
  dccutslow[3]  = 1.869 - sigmad * 0.0094;
  dccutshigh[3] = 1.869 + sigmad * 0.0094;
  dccutslow[4]  = 1.869 - sigmad * 0.0058;
  dccutshigh[4] = 1.869 + sigmad * 0.0058;
  dccutslow[5]  = 1.869 - sigmad * 0.0053;
  dccutshigh[5] = 1.869 + sigmad * 0.0053;

  dccutslow[6]  = 1.870 - sigmad * 0.005;
  dccutshigh[6] = 1.870 + sigmad * 0.005;

  ds0d0cutslow[0]  = 1.865 - sigmads * 0.0065;
  ds0d0cutshigh[0] = 1.865 + sigmads * 0.0065;
  ds0d0cutslow[1]  = 1.863 - sigmads * 0.0123;
  ds0d0cutshigh[1] = 1.863 + sigmads * 0.0123;
  ds0d0cutslow[2]  = 1.865 - sigmads * 0.004;
  ds0d0cutshigh[2] = 1.865 + sigmads * 0.004;
  ds0d0cutslow[3]  = 1.865 - sigmads * 0.006;
  ds0d0cutshigh[3] = 1.865 + sigmads * 0.006;
  ds0d0cutslow[4]  = 1.864 - sigmads * 0.0096;
  ds0d0cutshigh[4] = 1.864 + sigmads * 0.0096;

  ds0d0cutslow[5]  = 1.864 - sigmad * 0.0096;
  ds0d0cutshigh[5] = 1.864 + sigmad * 0.0096;
  ds0d0cutslow[6]  = 1.863 - sigmad * 0.0123;
  ds0d0cutshigh[6] = 1.863 + sigmad * 0.0123;
  ds0d0cutslow[7]  = 1.865 - sigmad * 0.0065;
  ds0d0cutshigh[7] = 1.865 + sigmad * 0.0065;
  ds0d0cutslow[8]  = 1.865 - sigmad * 0.0065;
  ds0d0cutshigh[8] = 1.865 + sigmad * 0.0065;

  dscd0cutslow[0]  = 1.865 - sigmads * 0.0065;
  dscd0cutshigh[0] = 1.865 + sigmads * 0.0065;
  dscd0cutslow[1]  = 1.863 - sigmads * 0.0123;
  dscd0cutshigh[1] = 1.863 + sigmads * 0.0123;
  dscd0cutslow[2]  = 1.865 - sigmads * 0.004;
  dscd0cutshigh[2] = 1.865 + sigmads * 0.004;
  dscd0cutslow[3]  = 1.865 - sigmads * 0.006;
  dscd0cutshigh[3] = 1.865 + sigmads * 0.006;
  dscd0cutslow[4]  = 1.864 - sigmads * 0.0096;
  dscd0cutshigh[4] = 1.864 + sigmads * 0.0096;

  dscd0cutslow[5]  = 1.864 - sigmad * 0.0096;
  dscd0cutshigh[5] = 1.864 + sigmad * 0.0096;
  dscd0cutslow[6]  = 1.863 - sigmad * 0.0123;
  dscd0cutshigh[6] = 1.863 + sigmad * 0.0123;
  dscd0cutslow[7]  = 1.865 - sigmad * 0.0065;
  dscd0cutshigh[7] = 1.865 + sigmad * 0.0065;
  dscd0cutslow[8]  = 1.865 - sigmad * 0.0065;
  dscd0cutshigh[8] = 1.865 + sigmad * 0.0065;

  dscdccutslow[0]  = 1.869 - sigmads * 0.0053;
  dscdccutshigh[0] = 1.869 + sigmads * 0.0053;
  dscdccutslow[1]  = 1.867 - sigmads * 0.01;
  dscdccutshigh[1] = 1.867 + sigmads * 0.01;
  dscdccutslow[2]  = 1.870 - sigmads * 0.005;
  dscdccutshigh[2] = 1.870 + sigmads * 0.005;
  dscdccutslow[3]  = 1.869 - sigmads * 0.0094;
  dscdccutshigh[3] = 1.869 + sigmads * 0.0094;
  dscdccutslow[4]  = 1.869 - sigmads * 0.0058;
  dscdccutshigh[4] = 1.869 + sigmads * 0.0058;
  dscdccutslow[5]  = 1.869 - sigmads * 0.0053;
  dscdccutshigh[5] = 1.869 + sigmads * 0.0053;

  dscdccutslow[6]  = 1.870 - sigmad * 0.005;
  dscdccutshigh[6] = 1.870 + sigmad * 0.005;

  ds0cutslow[0]  = 0.134;
  ds0cutshigh[0] = 0.142 + sigmads * 0.0016;
  ds0cutslow[1]  = 0.125;
  ds0cutshigh[1] = 0.142 + sigmads * 0.0026;

  dsccutslow[0]  = 0.1455 - sigmads * 0.001;
  dsccutshigh[0] = 0.1455 + sigmads * 0.001;
  dsccutslow[1]  = 0.1407 - sigmads * 0.001;
  dsccutshigh[1] = 0.1407 + sigmads * 0.001;

  return AppResult::OK;
}

AppResult DonutAnalysis::endJob( AbsEvent* event ) 
{
  delete readBestB;
  delete printTree;
  delete readCombD0;
  delete readCombDs0;
  delete readCombDp;
  delete readCombDsp;
  delete readDlD0;
  delete readDlDs0;
  delete readDlDp;
  delete readDlDsp;
  ErrMsg(routine) << "end job" << endmsg;
  return AppResult::OK;
}

AppResult DonutAnalysis::event( AbsEvent* event )
{
  if(_Verbose.value()) cout<<"manuelf: Starting DonutAnalysis event"<<endl;
  UpsilonList = Ifd< HepAList<BtaCandidate> >::get(event, keyUpsilonList.value() );
  TrackList = Ifd< HepAList<BtaCandidate> >::get(event, keyTrackList.value() );
  PhotonList= Ifd< HepAList<BtaCandidate> >::get(event, keyPhotonList.value() );
  IfdStrKey CTkey("ChargedTracks"), CNkey("CalorNeutral");
  AllTrackList = Ifd< HepAList<BtaCandidate> >::get(event, CTkey );
  AllPhotonList = Ifd< HepAList<BtaCandidate> >::get(event, CNkey );
  elecList = Ifd< HepAList<BtaCandidate> >::get(event, keyelecList.value() );
  muonList = Ifd< HepAList<BtaCandidate> >::get(event, keymuonList.value() );
  kaonList = Ifd< HepAList<BtaCandidate> >::get(event, keykaonList.value() );
  pionList = Ifd< HepAList<BtaCandidate> >::get(event, keypionList.value() );
  GTLList = Ifd< HepAList<BtaCandidate> >::get(event, keyGTLList.value() );
  GTVLList = Ifd< HepAList<BtaCandidate> >::get(event, keyGTVLList.value() );
  Pi0List = Ifd< HepAList<BtaCandidate> >::get(event, keyPi0List.value() );
  SoftPi0List = Ifd< HepAList<BtaCandidate> >::get(event, keySoftPi0List.value() );
  KsList = Ifd< HepAList<BtaCandidate> >::get(event, keyKsList.value() );
  MCTruthList = Ifd< HepAList<BtaCandidate> >::get(event, keyMCTruthList.value() );
  HepAList< EventInfo >* infoList = 
    Ifd<HepAList< EventInfo > >::get(event, keyeventInfoList.value() );
  eventInfo = infoList->first();
  _eventID = Ifd<AbsEventID>::get(event, keyeventIDList.value());
  AbsEventTag* tag = Ifd<AbsEventTag>::get(event);
  static const IfdStrKey keyGHit("GHit");
  BtaMcAssoc* truthMap = Ifd< BtaMcAssoc >::get(event, keyGHit);
  static const IfdStrKey keyGTL("/physicstools/trk/svt_tables/all");
  const EffTable* tableGTL = Ifd<EffTable>::get(gblPEnv, keyGTL);
  if (!tableGTL) {
    cout << "##### could not get GTL table!" << endl;
  }
  std::map<int,PidWeight> *elecMap = Ifd< map< int,PidWeight> >::get(event,keyelecList.value());
  std::map<int,PidWeight> *muonMap = Ifd< map< int,PidWeight> >::get(event,keymuonList.value());
  std::map<int,PidWeight> *kaonMap = Ifd< map< int,PidWeight> >::get(event,keykaonList.value());
  std::map<int,PidWeight> *pionMap = Ifd< map< int,PidWeight> >::get(event,keypionList.value());

  if(_Verbose.value()) cout<<"manuelf: After list definitions in donutanalysis"<<endl;
  //from DonutMCTagger
  int MCType,MCSubmode,MCCombmode,MCPions,MCTaumode,MCBrem,MCScatter,MCD,MC2Body,MCDoubleSL,
    nTrueP,nTrueKL,nTrueNu,isVub,isBzero,uidBTag,uidD,uidLep,uidNu,uidPi0,uidDssPi0,
    MCCombB,MCCombDs,MCDssmode,MCUnsim;
  int trueLepMother=0,trueLepCharge;
  int ntrueLep, trueMotherLep1, trueMotherLep2, trueMotherLep3, trueLundLep1,
    trueLundLep2, trueLundLep3, nMufromPi;
  float truePLep1, truePLep2, truePLep3;
  int MCComblong = -1,MCComblongD = -1;
  // manuelf: tag variables
  int trueTagPi, trueTagK, trueTagKs, trueTagPi0, trueTagE, trueTagMu, trueTagDs, trueTagD;
  int trueTagYPi, trueTagYK, trueTagYKs, trueTagYPi0;
  int trueTagElse;
  float truePLep,truePD,trueDmass,trueCTL,trueCTV,trueChi,trueQ2,trueTauFlight;
  // Need MCType for eff calculations and need isBzero for skipBPlus mode
  // Leave the rest down below to speed up module
  tag -> getInt(MCType     ,"MCType");
  tag -> getInt(isBzero    ,"isBzero");
  if (_skipBPlus.value()) {
    if (isBzero == 0 ) return AppResult::OK;
  }
  iterElec = new HepAListIterator<BtaCandidate>(*elecList);
  iterMuon = new HepAListIterator<BtaCandidate>(*muonList);
  iterKaon = new HepAListIterator<BtaCandidate>(*kaonList);
  iterPion = new HepAListIterator<BtaCandidate>(*pionList);
  iterTracks = new HepAListIterator<BtaCandidate>(*TrackList);
  iterAllTracks = new HepAListIterator<BtaCandidate>(*AllTrackList);
  iterAllPhotons = new HepAListIterator<BtaCandidate>(*AllPhotonList);
  iterPhotons = new HepAListIterator<BtaCandidate>(*PhotonList);
  iterGTL = new HepAListIterator<BtaCandidate>(*GTLList);
  iterGTVL = new HepAListIterator<BtaCandidate>(*GTVLList);
  iterPi0 = new HepAListIterator<BtaCandidate>(*Pi0List);
  iterSoftPi0 = new HepAListIterator<BtaCandidate>(*SoftPi0List);
  iterSoftPi02 = new HepAListIterator<BtaCandidate>(*SoftPi0List);
  iterKs = new HepAListIterator<BtaCandidate>(*KsList);
  iterMCTruth = new HepAListIterator<BtaCandidate>(*MCTruthList);

  HepAListIterator<BtaCandidate> iterUpsilon(*UpsilonList);
  HepLorentzVector Upsilon = eventInfo -> cmFrame();
  Hep3Vector CMboost = Upsilon.boostVector();

  if(_Verbose.value()) cout<<"manuelf: Getting usrBlock"<<endl;
  UsrCandBlock *usrBlock = UsrIfd<UsrCandBlock>::get( event, keycandData.value());

  if (usrBlock == 0) 
    ErrMsg(fatal) << "could not find the usr data block " << keycandData.value() << endmsg;

  BtaOpMakeTree comb;
  VtxFitterOper *fitter(0);

  BtaCandidate *anUpsilon, *theBestUpsilon(0), *theBestBReco(0),*theBestLepton(0),*theBestD(0),
    *theBestBareD(0);
  HepLorentzVector theBestDP4, theBestMissP4;
  BtaCandidate *aBReco=0, *aLepton=0, *aDNoVtx=0, *aD=0;
  BtaCandidate *bareDNoVtx=0,*bareD=0,*thePisoft=0;
  BtaCandidate *aDauD=0, *aTrack=0, *aPhoton=0;

  //calculated locally
  int candType=-1,candDstarType=-1,candDType=-1,candIsMixed=-1,candExtraTracks=_maxExtraTracks.value(),
    candRejectedTracks=0,candRejectedPhotons=0;
  int candBntCha=-99,candBntNeu=-99,candDntCha=-99,candDntNeu=-99;
  int candBnCharged=-99,candBnNeutral=-99,candDnCharged=-99,candDnNeutral=-99;
  int candIsMu=-1,candLepCharge=-99,candDLund=0,candLepTru=-1,candDTru=-1,candBTru=-1,candIsBrem=-1, candLepIsPi=-1;
  int candFitStat=-1,candFitNdof=-1,candBMode=-1,candTagLeptons=-1,candTagChargedMult=-1,candTagNeutMult=-1;
  int SeedMode, YMode;
  float candIntPur=-1.0,candMES=-99.,candMES_orig=-99.,candDeltaE=-99.;
  float candBTagDeltaP,candBTagDeltaPx,candBTagDeltaPy,candBTagDeltaPz,candBTagDeltaE,candBTagYMass=-99.;
  float candBTagDmass=-9, candBTagDeltam = -9;
  float candDDeltaP,candDDeltaPx,candDDeltaPy,candDDeltaPz,candDDeltaE;
  float candLepDeltaP,candLepDeltaPx,candLepDeltaPy,candLepDeltaPz,candLepDeltaE;
  float candNuDeltaP,candNuDeltaPx,candNuDeltaPy,candNuDeltaPz,candNuDeltaE;
  float candBDPx,candBDPy,candBDPz,candBDE;
  float candMissPx,candMissPy,candMissPz,candMissE;
  float candDmass=-99.,candDeltam=-99.,candEExtra=99.9,candUsedEnergy=-99.,candM2=-99.,eventM2=-99.;
  float candKsMass=-99.,candPi0Mass=-99.,candPi0P=-99.,candMExtra=-99.,candMTot=-99.;
  float candPMiss=-99.,candEMiss=-99.,candThetaMiss=-99.,candPhiMiss=-99.,candThetaDL=-99.,candHelicity=-99.,candQ2=-99.;
  float candPLep=-99.,candPstarLep=-99.,candThetaLep=-99.,candBremE=-99.,candBremAngle=-99.;
  float candFitChi2=-99.,candFitProb=-99.,candPD=-99.,candPstarD=-99.,candEstarD=-99.;
  float candM2Tru=-99.,candM2TruBTag=-99.,candM2TruD=-99.,candM2TruLep=-99.;
  float candPisoftP=-99.,candPisoftTheta=-99.,candPisoftPhi=-99.,candPisoftE1=-99.,candPisoftE2=-99.,candPisoftM2=-99.;
  float candCosT = -99.;
  int nGTLHard,nGTLSoft,nOtherTracks,nCandK,nCandPi,candExtraPhotons=0,candUsedPhotons=-1;
  float trkWeight,candKSdxy,candKSTheta,candKSPT,candPPi0;
  float kaonWeight,pionWeight,elecWeight,muonWeight;
  HTValOrderedVector<float> candPKaons,candPPions;
  float wBF, wComb, wFF;

  int tempBntCha=-99,tempBntNeu=-99,tempDntCha=-99,tempDntNeu=-99;
  int tempBnCharged=-99,tempBnNeutral=-99,tempDnCharged=-99,tempDnNeutral=-99;
  int tempType,tempDstarType,tempDType,tempIsMixed,tempExtraTracks,tempIsBrem;
  int tempFitStat,tempFitNdof,tempExtraPhotons=0,tempUsedPhotons;
  float tempBremE,tempBremAngle,tempHelicity,tempUsedEnergy,
    tempKsMass,tempPi0Mass,tempPi0P,tempMExtra,tempMTot;
  float tempPisoftP=-99.,tempPisoftTheta=-99.,tempPisoftPhi=-99.,tempPisoftE1=-99.,tempPisoftE2=-99.,tempM2=-99.;
  float tempFitChi2;
  bool pid,wrongsign,pass=false,status;
  int ncands = 0,ncands1=0,ncands2=0,ncands3=0,ncands4=0,nups1=0,nups2=0,nups3=0,nups4=0;

  int npi0;
  HTValOrderedVector<float> cosEpi0, mpi0,ppi0,dmpi0,mm2pi0,pmisspi0,eextrapi0,e1pi0,e2pi0,q2pi0;
  HTValOrderedVector<int> bestepi0,bestmasspi0;
//   int nks;
//   HTValOrderedVector<float> mks,pks,vxks,vyks,vzks,mm2ks;
//   HTValOrderedVector<int> extratracksks;
  int ndpi0;
  HTValOrderedVector<float> mm2dpi0,eextradpi0,dmdpi0,pmissdpi0,p1dpi0,p2dpi0,m1dpi0,m2dpi0;
  HTValOrderedVector<int> bestdpi0;

  HTValOrderedVector<int> candBModeL, candBntChaL, candBntNeuL, candDntChaL, candDntNeuL, candTypeL,candLepTruL;  
  HTValOrderedVector<float> candEExtraL, candDeltaEL, candMvaBestBL;
  float candMvaBestB=-99.,candMvaComb=-99.,candMvaDl=-99.,candMvaDssComb=-99.,candMvaDssDl=-99., tempMvaBestB=-99.;
  float tempBTagYmass=99.;
  int tempTagNeutMult;

  bool hit[4][8];
  for (int a = 0 ; a < 4 ; a ++)
    for (int b = 0 ; b < 8 ; b ++)
      hit[a][b] = false;

  if(_Verbose.value()) cout<<"manuelf: Starting the While anUpsilon =============";
  while (anUpsilon = iterUpsilon()) {
    fitter = 0;
    HepAListIterator<BtaCandidate> iterDauUps = anUpsilon->daughterIterator();
    // manuelf: The Upsilon daugthers are always B, D, l, in this order
    aBReco  = iterDauUps();
    aDNoVtx = iterDauUps();
    aLepton = iterDauUps();
     
    status = true;
    status &= usrBlock->get(*anUpsilon,_decayMode);
    // manuelf: The empty modes are cut
    int mode = (int)_decayMode;
    SeedMode = (mode/100)%100; YMode = mode%100;
    if(_GroupEmpty[SeedMode][YMode]) continue;

    status &= usrBlock->get(*anUpsilon,_purity);
    status &= usrBlock->get(*anUpsilon,_intpurity);
    status &= usrBlock->get(*anUpsilon,_mES);
    status &= usrBlock->get(*anUpsilon,_deltaE);
    BtaBVariables BVars(aBReco->p4());
    tempDeltaE = _deltaE;
    TString folder = (TString)_StemFolder.value();
    if(mEScorrection){
      tempMES = PidCorrectMesMean::correctMes(_mES,(int)_eventID->run());	
    } else tempMES = _mES;
    if(folder.Contains("OffPeak")) tempMES += 0.02;
//     if (BVars.status() == BtaBVariables::GOOD){    
//       tempDeltaE = BVars.deltaE();
//       tempMES = BVars.m_ES();
//     }
    //break up upsilon, record topology
    // Checking if there is Bremsstrahlung electrons first
    if (aLepton->nDaughters() > 0) {
      HepAListIterator<BtaCandidate> iterDauE = aLepton->daughterIterator();
      BtaCandidate* atrk = iterDauE();
      BtaCandidate* aPhot = iterDauE();
      tempIsBrem = 1;
      tempBremE = aPhot->energy();
      tempBremAngle = atrk->p4().angle(aPhot->p3());
    }
    else {
      tempIsBrem = 0;
      tempBremE = tempBremAngle = -1;
    }
    tempType = 0;
    if (!aDNoVtx->pdtEntry()) {
      cout << "Error!!!!! Can't figure out what flavor of D I have!!!!!" << endl;
      continue;
    }
    if (abs((int)aDNoVtx->pdtEntry()->lundId()) == (int)PdtLund::D0)          tempType=1;
    if (abs((int)aDNoVtx->pdtEntry()->lundId()) == (int)PdtLund::D_star0)     tempType=2;
    if (abs((int)aDNoVtx->pdtEntry()->lundId()) == (int)PdtLund::D_plus)      tempType=3;
    if (abs((int)aDNoVtx->pdtEntry()->lundId()) == (int)PdtLund::D_star_plus) tempType=4;
    // Checking for the case of charged D(*) if the charge is right
    tempIsMixed = 0;
    if (tempType>2) {
      if ((int)aBReco->pdtEntry()->lundId() * (int)aDNoVtx->pdtEntry()->lundId() < 0) tempIsMixed = 1;
    }

    if (tempType==1) nups1 ++;
    if (tempType==2) nups2 ++;
    if (tempType==3) nups3 ++;
    if (tempType==4) nups4 ++;

    if (tempMES > 5.27) dehist[tempType-1] -> accumulate(tempDeltaE);
    if (fabs(tempDeltaE) < 0.05) meshist[tempType-1] -> accumulate(tempMES);
    // The D daughter of D* is stored, as well as the soft pion
    tempDstarType = 0;
    if (abs((int)aDNoVtx->pdtEntry()->lundId()) == (int)PdtLund::D_star0 || 
	abs((int)aDNoVtx->pdtEntry()->lundId()) == (int)PdtLund::D_star_plus) {
      HepAListIterator<BtaCandidate> iterDauD = aDNoVtx->daughterIterator();
      bareDNoVtx = iterDauD();
      thePisoft = iterDauD();
      if ((int)thePisoft->pdtEntry()->lundId() == (int)PdtLund::pi0) {
	if (aDNoVtx->charge() == 0) tempDstarType = 1;
	else tempDstarType = 2;
      }
      else {
	if (aDNoVtx->charge() == 0) tempDstarType = 2;
	else tempDstarType = 1;
      }
      tempPisoftP = thePisoft->p();
      tempPisoftTheta = thePisoft->p4().theta();
      tempPisoftPhi = thePisoft->p4().phi();
      if ((int)thePisoft->pdtEntry()->lundId() == (int)PdtLund::pi0) {
	tempPisoftE1 = tempPisoftE2 = -1;
	HepAListIterator<BtaCandidate> iterDauPi0 = thePisoft->daughterIterator();
	// The energy of the 2 most energetic photons is stored in order E1, E2
	while (BtaCandidate *aDauPi0 = iterDauPi0()) {
	  if (tempPisoftE1 < 0) {
	    tempPisoftE1 = tempPisoftE2 = aDauPi0->energy();
	  }
	  else if (aDauPi0->energy() > tempPisoftE1) tempPisoftE1 = aDauPi0->energy();
	  else tempPisoftE2 = aDauPi0->energy();
	}
	iterDauPi0.rewind();
      }
      else tempPisoftE1 = tempPisoftE2 = -1;
    }
    else bareDNoVtx = aDNoVtx;

    if (bareDNoVtx -> charge()) tempDType = DType(bareDNoVtx);
    else tempDType = D0Type(bareDNoVtx);

    if (!hit[tempType-1][0]) {
      hit[tempType-1][0] = true;
      effhist[tempType-1] -> accumulate(MCType,0);
    }

//     if(_Verbose.value()) {cout<<endl<<"New Upsilon. mES: "<<tempMES<<" dE: "<<tempDeltaE<<" Angle: "<<
// 			    aLepton->p4().theta()<<" Dtype: "<<tempType<<endl;
//     }
    //
    // Wrong sign for the neutral D. D0 goes with l- and antiD0 with l+
    // I do not understand the tempDType < 4 so it's commented out
    wrongsign = false;
    if (aDNoVtx->charge() == 0 ){//&& tempDType < 4) {
      if ((int)bareDNoVtx->pdtEntry()->lundId() == (int)PdtLund::D0 && aLepton->charge() > 0)
	wrongsign = true;
      if ((int)bareDNoVtx->pdtEntry()->lundId() == (int)PdtLund::anti_D0 && aLepton->charge() < 0)
	wrongsign = true;
    }
    if (wrongsign && fitter) {delete fitter; delete bareD; if (tempDstarType) delete aD;}
    if (wrongsign) continue;

    //
    //deltaE, mES and lepton angle cuts
    //
    bool passDE=false;
    if (tempType == 1 || tempType == 3) passDE = fabs(tempDeltaE) < _deltaEcut.value();
    if (tempType == 2 || tempType == 4) passDE = fabs(tempDeltaE) < _deltaEcut.value();
    if (_deltaESideband.value()) passDE = !passDE;

    if (!passDE && fitter) {delete fitter; delete bareD; if (tempDstarType) delete aD;}
    if (!passDE) continue;
    if (!hit[tempType-1][1]) {
      hit[tempType-1][1] = true;
      effhist[tempType-1] -> accumulate(MCType,1);
    }

    if (tempMES < _mEScut.value()) continue;
    if (aLepton->p4().theta() < 0.4 || aLepton->p4().theta() > 2.6) continue;

    if(_Verbose.value()) {cout<<endl<<"Printing the reconstructed upsilon. MCType "<<MCType<<endl;
    cout<<printTree -> print(*anUpsilon);
     cout<<endl;
    }
    //
    //lepton PID
    //
    pid = true;
    if (abs((int)aLepton->pdtEntry()->lundId()) == (int)PdtLund::mu_minus) {
      if (!isMuon(aLepton)) pid = false;
      //if (isElectron(aLepton)) pid = false;
      //if (isKaon(aLepton)) pid = false;
      if(_Verbose.value()) cout<<"manuelf: muon list: "<<pid;
      if (!isGTL(aLepton)) pid = false;
      if(_Verbose.value()) cout<<", GTL list: "<<pid<<endl;
    } //muon
    else {
      if (!isElectron(aLepton)) pid = false;
      //if (isMuon(aLepton)) pid = false;
      //if (isKaon(aLepton)) pid = false;
      if (!isGTL(aLepton)) pid = false;
      if (aLepton->p() < .3) pid = false;
    } //electron

    if ((!pid) && fitter) {delete fitter; delete bareD; if (tempDstarType) delete aD;}
    if (!pid) continue;
    if (!hit[tempType-1][2]) {
      hit[tempType-1][2] = true;
      effhist[tempType-1] -> accumulate(MCType,2);
    }

    // manuelf: Checking the PID's of the D daughters and that its momentum is < 4
    tempKsMass = tempPi0Mass = tempPi0P = 0;
    HepAListIterator<BtaCandidate> iterDauD = bareDNoVtx->daughterIterator();
    if(_Verbose.value()) cout<<"manuelf: Checking PID of D daughters. ";
    while (aDauD = iterDauD()) {
      if (abs((int)aDauD->pdtEntry()->lundId()) == (int)PdtLund::K_plus) {
	if (aDauD->p() > 4) pid = false;
	if (!isKaon(aDauD))    {pid = false;}
	//if (isElectron(aDauD)) {pid = false;}
	//if (isMuon(aDauD))     {pid = false;}
	if(_Verbose.value()) cout<<"K "<<pid<<" ";
      } //k
      if (abs((int)aDauD->pdtEntry()->lundId()) == (int)PdtLund::pi_plus) {
	if (aDauD->p() > 4) pid = false;
	if (!isPion(aDauD))    {pid = false;}
	//if (isElectron(aDauD)) {pid = false;}
	//if (isMuon(aDauD))     {pid = false;}
	//if (isKaon(aDauD))     {pid = false;}
	if(_Verbose.value()) cout<<"pi "<<pid<<" ";
      } //pi
      if ((int)aDauD->pdtEntry()->lundId() == (int)PdtLund::pi0) {
	HepLorentzVector rawP4;
	HepAListIterator<BtaCandidate> iterDauPi0 = aDauD->daughterIterator();
	while (BtaCandidate *aDauPi0 = iterDauPi0()) {
	  rawP4 += aDauPi0->p4();
	}
	iterDauPi0.rewind();
	if (rawP4.mag() < .125 || rawP4.mag() > .145) pid = false;
	if (aDauD->p() > 4) pid = false;
	tempPi0Mass = rawP4.mag();
	tempPi0P = aDauD->p();
	if(_Verbose.value()) cout<<"pi0 "<<pid<<" ";
      } //pi0
      if ((int)aDauD->pdtEntry()->lundId() == (int)PdtLund::K_S0) {
	tempKsMass = aDauD->mass();
	if (tempKsMass < 0.491 || tempKsMass > 0.506) pid = false;
	HepAListIterator<BtaCandidate> dauKS = aDauD->daughterIterator();
	while (BtaCandidate *aKSPion = dauKS()) {
	  if (aKSPion->p() > 3) pid = false;
	}
	if(_Verbose.value()) cout<<"KS "<<pid<<" ";
      } //Ks
    } //dauD
    // The soft pi coming from D*->Dpi has to have p<400MeV
    iterDauD.rewind();
    if (tempDstarType)
      if (thePisoft->p() > .4) pid = false;

    if ((!pid) && fitter) {delete fitter; delete bareD; if (tempDstarType) delete aD;}
    if (!pid) continue;
    if (!hit[tempType-1][3]) {
      hit[tempType-1][3] = true;
      effhist[tempType-1] -> accumulate(MCType,3);
    }


    // The kkk thingy kills a fraction randomly selected of the D* candidates
    tempDmass = aDNoVtx->mass();
    tempDeltam = aDNoVtx->mass() - bareDNoVtx->mass();
    int kkk = (MCType==2||MCType==4||MCType==6||MCType==8||MCType==10||MCType==12);
    if (MCType == 0 || MCType > 12)
      if (truthMap)
	kkk = (int)(truthMap->mcFromReco(aDNoVtx,0,true) != 0);
    kkk = 0;
    if(_Verbose.value()) cout<<"manuelf: D masses checks"<<endl;
    // Some D modes are rejected and the D masses and the D*-D mass difference are checked
    if (tempType == 1) pass = checkD0(tempDmass,tempDType);
    if (tempType == 2) pass = checkDs0(tempDmass,tempDeltam,tempDType,tempDstarType,kkk);
    if (tempType == 3) pass = checkD(tempDmass,tempDType);
    if (tempType == 4) pass = checkDs(tempDmass,tempDeltam,tempDType,tempDstarType,kkk);
    if ((!pass) && fitter) {delete fitter; delete bareD; if (tempDstarType) delete aD;}
    if (!pass) continue;
    if (!hit[tempType-1][4]) {
      hit[tempType-1][4] = true;
      effhist[tempType-1] -> accumulate(MCType,4);
    }
    if(_Verbose.value()) cout<<"manuelf: Passed D masses checks"<<endl;

    // Vertex the cands
    // The masses of the BReco, D and D* are constrained. In the latter case, there also a
    // beam constraint and D* candidate is created again from the D and pisoft
    BtaCandidate *copyBReco = new BtaCandidate(*aBReco);
    copyBReco -> invalidateFit();
    setMassConstraint(*aBReco);
    BtaCandidate *copyD;
    BtaCandidate *copyBareD(0);
    if (tempDstarType) {
      copyBareD = new BtaCandidate(*bareDNoVtx);
      copyBareD -> invalidateFit();
      setMassConstraint(*copyBareD);
      copyD = comb.create(*copyBareD,*thePisoft);
      copyD->setType(aDNoVtx->pdtEntry());
      setMassConstraint(*copyD);
      setBeamConstraint(*copyD,eventInfo,true);
    } else {
      copyD = new BtaCandidate(*aDNoVtx);
      copyD -> invalidateFit();
      setMassConstraint(*copyD);
    }
    // Creation of the signal B in candB with the lepton, the D and the nu calculated as the missing stuff
    BtaCandidate *copyLepton = new BtaCandidate(*aLepton);
    HepLorentzVector pNu = Upsilon - aBReco->p4() - aDNoVtx->p4() - aLepton->p4();
    BtaCandidate tempNu(pNu);
    tempNu.setType(Pdt::lookup(PdtLund::nu_e));
    BtaCandidate *candB = comb.create(*copyD,*copyLepton,tempNu);
    candB->setType(copyBReco->pdtEntry()->conjugate());
    //setMassConstraint(*candB);
    // Creation of the Upsilon, setting of the beam constraint and fitting
    BtaCandidate *candUps = comb.create(*aBReco,*candB);
    candUps->setType(Pdt::lookup(PdtLund::Upsilon_4S));
    //setMassConstraint(*candUps);
    setBeamConstraint(*candUps,eventInfo);
    setEnergyConstraint(*candUps,eventInfo);
    VtxFitterOper *upsfitter = new VtxFitterOper(*candUps,VtxAbsAlgorithm::TreeFitter);
    upsfitter->fit();
    BtaCandidate fittedUps = upsfitter -> getFittedTree();
    tempFitStat = fittedUps.decayVtx()->status();
    tempFitNdof = fittedUps.decayVtx()->nDof();
    tempFitChi2 = fittedUps.decayVtx()->chiSquared();
    // Cutting on failed fits
    //if (tempFitStat != 0) continue;
    if (tempFitStat == BtaAbsVertex::Failed) continue;
    if (tempFitStat == BtaAbsVertex::NonConverged) continue;
    if(_Verbose.value()) cout<<"manuelf: Passed Vertex fitting convergence"<<endl;
    // The fitted D and D* are stored
    if (tempDstarType) {
      BtaCandidate fittedD = upsfitter -> getFitted(*copyD);
      aD = new BtaCandidate(fittedD);
      BtaCandidate fittedBareD = upsfitter -> getFitted(*copyBareD);
      bareD = new BtaCandidate(fittedBareD);
    } else {
      BtaCandidate fittedD = upsfitter -> getFitted(*copyD);
      aD = new BtaCandidate(fittedD);
      bareD = aD;
    }
    // The missing mass and momentum are obtained from the fitted neutrino
    BtaCandidate pmissfitups = upsfitter -> getFitted(tempNu);
    tempM2 = pmissfitups.p4().mag2();
    HepLorentzVector pMiss = pmissfitups.p4();

    delete upsfitter;
    delete copyBReco; delete copyD; delete copyLepton; delete candB; delete candUps;
    if (copyBareD) delete copyBareD;

    if (!hit[tempType-1][5]) {
      hit[tempType-1][5] = true;
      effhist[tempType-1] -> accumulate(MCType,5);
    }

    //
    // extraTracks
    // Events with extra tracks are cut
    HepLorentzVector pextratracks;
    tempExtraTracks = 0;
    while (aTrack = (*iterTracks)()) {
      if (aTrack->overlaps(*aBReco)) continue;
      if (aTrack->overlaps(*aD)) continue;
      if (aTrack->overlaps(*aLepton)) continue;
      tempExtraTracks ++;
      pextratracks += aTrack->p4();
    }
    iterTracks->rewind();
    if (tempExtraTracks > candExtraTracks && fitter) {delete fitter; delete bareD; if (tempDstarType) delete aD;}
    if (tempExtraTracks > candExtraTracks) continue;
    if(_Verbose.value()) cout<<"manuelf: Passed tempExtraTracks > candExtraTracks"<<endl;
    if (tempExtraTracks == 0)
      if (!hit[tempType-1][6]) {
	hit[tempType-1][6] = true;
	effhist[tempType-1] -> accumulate(MCType,6);
      }

    // Calculation of CM lepton and D and q^2 cut
    HepLorentzVector thisB = Upsilon - aBReco->p4(); // Momentum of SigB
    HepLorentzVector qvec = thisB - aD->p4();
//     if(qvec.mag2()<=4.) continue;

    //================================ Finding truthmatched particles ===================================
    if (truthMap) {
      BtaCandidate *temptrueLepton;
      temptrueLepton = truthMap->mcFromReco(aLepton,0,true);
      if(!temptrueLepton && (abs((int)(aLepton->pdtEntry()->lundId()))==(int)PdtLund::e_minus)){
	HepAListIterator<BtaCandidate> dauglepIter = aLepton->daughterIterator();
	while (BtaCandidate* son = dauglepIter()) {
	  int pdt = son->pdtEntry()->lundId();
	  if(abs(pdt)==(int)PdtLund::e_minus) temptrueLepton = truthMap->mcFromReco(son,0,true);
	}
      }
      PdtLund::LundType sigBLund=PdtLund::null;
      if(aD->charge()==0){
	if(aLepton->charge()>0) sigBLund = PdtLund::B_plus;
	else sigBLund = PdtLund::B_minus;
      } else {
	if(aLepton->charge()>0) sigBLund = PdtLund::B0;
	else sigBLund = PdtLund::anti_B0;
      }
      int uidSig = 0, uidTag = 0, uid1 = 0, uid2 = 0;
      if (MCTruthList){
	BtaCandidate* aTru = 0;
	while (aTru = (*iterMCTruth)()) {
	  if ((int)aTru->pdtEntry()->lundId() == PdtLund::Upsilon_4S) {
	    HepAListIterator<BtaCandidate> dauXIter = aTru->daughterIterator();
	    BtaCandidate* dauX = 0;
	    PdtLund::LundType dauLund;
	    while (0 != (dauX = dauXIter())) {
	      dauLund = dauX->pdtEntry()->lundId();
	      if(abs(dauLund)==PdtLund::B0 || abs(dauLund)==PdtLund::B_plus){
		if(temptrueLepton){
		  int BLepId = motherB(temptrueLepton);
		  if(BLepId != dauX->uid()) uidTag = dauX->uid();
		  else uidSig = dauX->uid();
		  if(uid1==0) uid1 = dauX->uid();
		  else uid2 = dauX->uid();
		} else {
		  if(dauLund == sigBLund) uidSig = dauX->uid();
		  if(dauLund == aBReco->pdtEntry()->lundId()) uidTag = dauX->uid();
		  if(uid1==0) uid1 = dauX->uid();
		  else uid2 = dauX->uid();
		}
	      }
	      if(uid2!=0) break;
	    }
	  }
	}// While iterMCTruth
	if (uidSig ==0 || uidTag ==0 || uidSig == uidTag){
	  uidSig = uid1; uidTag = uid2;
	}
      }
      iterMCTruth->rewind();
      tempBntCha=tempBntNeu=tempDntCha=tempDntNeu=0;   
      tempBnCharged=tempBnNeutral=tempDnCharged=tempDnNeutral=0;   
      nTruthMatched(aD, tempDntCha, tempDntNeu, tempDnCharged, tempDnNeutral, truthMap, uidSig);
      nTruthMatched(aBReco, tempBntCha, tempBntNeu, tempBnCharged, tempBnNeutral, truthMap, uidTag);
    }
    BtaCandidate *trueD,*trueLepton,*trueB;
    if (truthMap) {
      trueD = truthMap->mcFromReco(bareD,0,true);
      trueB = truthMap->mcFromReco(aBReco,0,true);
      trueLepton = truthMap->mcFromReco(aLepton,0,true);
      if(!trueLepton && (abs((int)(aLepton->pdtEntry()->lundId()))==(int)PdtLund::e_minus)){
	HepAListIterator<BtaCandidate> dauglepIter = aLepton->daughterIterator();
	while (BtaCandidate* son = dauglepIter()) {
	  int pdt = son->pdtEntry()->lundId();
	  if(abs(pdt)==(int)PdtLund::e_minus) trueLepton = truthMap->mcFromReco(son,0,true);
	}
      }
	
    }
    else trueD = trueLepton = trueB = 0;

    //==================================== EExtra calculation  ===========================================
    // EExtra, EUsed, Extra mass and total p of the D plus pExtra are calculated
    tempEExtra = 0.0;
    tempExtraPhotons = 0;
    tempUsedEnergy = 0.0;
    tempUsedPhotons = 0;
    HepLorentzVector pextra;
    while (aPhoton = (*iterPhotons)()) {
      if (aPhoton->overlaps(*anUpsilon)) {
	tempUsedPhotons ++;
	tempUsedEnergy += aPhoton->energy();
      } else {
	tempExtraPhotons ++;
	tempEExtra += aPhoton->energy();
	pextra += aPhoton->p4();
      }
    }
    iterPhotons->rewind();
    tempMExtra = pextra.mag();
    pextra += aD->p4();
    tempMTot = pextra.mag();

    //=============================== Rejected tracks and photons  =======================================
    // Extra tracks and photons counted (wrt dirty lists)
    tempRejectedTracks = 0;
    while (aTrack = (*iterAllTracks)()) {
      if (aTrack->overlaps(*anUpsilon)) continue;
      tempRejectedTracks ++;
    }
    iterAllTracks->rewind();
    tempRejectedPhotons = 0;
    while (aTrack = (*iterAllPhotons)()) {
      if (aTrack->overlaps(*anUpsilon)) continue;
      tempRejectedPhotons ++;
    }
    iterAllPhotons->rewind();
    //==================================== Mass of the Y system  =========================================
    HepLorentzVector p4tag = aBReco->p4();
    HepAListIterator<BtaCandidate> iterDauTag = aBReco->daughterIterator();
    BtaCandidate *aDauTag = iterDauTag();
    HepLorentzVector p4tagseed = aDauTag->p4();
    int dauLund = (int)aDauTag->pdtEntry()->lundId();
    if (abs(dauLund) == (int)PdtLund::D0 || 
	abs(dauLund) == (int)PdtLund::D_plus || 
	abs(dauLund) == (int)PdtLund::J_psi || 
	abs(dauLund) == (int)PdtLund::D_s_plus) {
      tempBTagDmass = p4tagseed.m();
      tempBTagDeltam = 0;
    }
    if (abs(dauLund) == (int)PdtLund::D_star0 || 
	abs(dauLund) == (int)PdtLund::D_star_plus || 
	abs(dauLund) == (int)PdtLund::D_s_star_plus) {
      HepAListIterator<BtaCandidate> iterDDau = aDauTag -> daughterIterator();
      BtaCandidate* DDau;
      while (DDau = iterDDau()) {
	int DdauLund = (int)DDau->pdtEntry()->lundId();
	if (abs(DdauLund) == (int)PdtLund::D0 || 
	    abs(DdauLund) == (int)PdtLund::D_plus || 
	    abs(DdauLund) == (int)PdtLund::D_s_plus) {
	  tempBTagDmass = DDau->p4().m();
	  tempBTagDeltam = p4tagseed.m()-tempBTagDmass;
	}
      } 
      iterDDau.rewind();
    }   
    iterDauTag.rewind();
    p4tag -= p4tagseed;
    tempBTagYmass = p4tag.m();
    tempTagChargedMult = 0; tempTagNeutMult = 0;
    countMult(aBReco,tempTagChargedMult,tempTagNeutMult);

    //====================================================================================================
    //===================================== Best B selection =============================================
    //====================================================================================================

    bool keep = false;
    if (_doBestBEEx.value()){ 
      if((tempType<3&&tempEExtra<0.2)||(tempType==3&&tempEExtra<0.15)||(tempType==4&&tempEExtra<0.3)) 
	keep = true;
    } else { 
      tempMvaBestB = readBestB->EvaluateMVA("BestB");
      if(tempMvaBestB > 0) keep = true;
    }
    if(keep){
      ncands ++;
      if (tempType==1) ncands1 ++;
      if (tempType==2) ncands2 ++;
      if (tempType==3) ncands3 ++;
      if (tempType==4) ncands4 ++;
      candBModeL.append((int)_decayMode);
      candEExtraL.append(tempEExtra);
      candBntChaL.append(tempBntCha);
      candBntNeuL.append(tempBntNeu);
      candDntChaL.append(tempDntCha);
      candDntNeuL.append(tempDntNeu);
      candTypeL.append(tempType);
      candLepTruL.append((int)(trueLepton != 0));
      candDeltaEL.append(tempDeltaE);
      if (_doBestBEEx.value() == false) candMvaBestBL.append(tempMvaBestB);
    }
    // For the same mva (or EExtra), the one with the smallest DeltaE is kept
    if (_doBestBEEx.value()){ 
      if ((tempEExtra < candEExtra || tempEExtra == candEExtra && fabs(tempDeltaE)<fabs(candDeltaE)))
	keep = true;
      else keep = false;
     } else { 
      if (tempMvaBestB>candMvaBestB || tempMvaBestB==candMvaBestB && fabs(tempDeltaE)<fabs(candDeltaE))
	keep = true;
      else keep = false;
     }

    if ((!keep) && fitter) {delete fitter; delete bareD; if (tempDstarType) delete aD;}
    if (!keep) continue;
    if(_Verbose.value()) cout<<"manuelf: Passed all the EExtra stuff"<<endl;

//     if I ever want to use purity again, here it is
//     if(){   
//       double purity = _TruePurity[SeedMode][YMode];
//       if (keep && (purity>tempPuri || purity==tempPuri && tempEExtra < candEExtra)){
// 	keep = true;
// 	tempPuri = purity;
//       } else keep = false;
//     } 
    //====================================================================================================
    //====================================================================================================

    // Storing the true particles associated to the reconstructed ones

    Hep3Vector thisBBoost = thisB.boostVector();
    HepLorentzVector pLepCM(aLepton->p4());
    pLepCM.boost(-thisBBoost);
    HepLorentzVector pDCM(aD->p4());
    pDCM.boost(-thisBBoost);
    // Momentum tranfer (q) and addition of pMiss to the p of the pion from the D*
    HepLorentzVector pMissPisoft = pMiss;
    if (tempDstarType > 0) {
      pMissPisoft += thePisoft->p4();
    }
    HepLorentzVector pMissExtraTracks = pMiss - pextratracks;
    // Helicity angle in the W rest frame
    HepLorentzVector wstarCMB = qvec;
    wstarCMB.boost(-thisBBoost);
    HepLorentzVector pLepCMHelicity(aLepton->p4());
    pLepCMHelicity.boost(-thisBBoost);
    Hep3Vector WFromBBoost = wstarCMB.boostVector();
    pLepCMHelicity.boost(-WFromBBoost);
    tempHelicity = pLepCMHelicity.angle(wstarCMB.vect());

    //==================================== cos(Theta Thrust)  =========================================

    BtaBooster* upsBooster = new BtaBooster(Upsilon);
    // Thrust of the reconstructed (tag) B
    BtaTreeNavigator TreeTagB( *aBReco );
    std::list<BtaCandidate*> stableDaughtersSE = TreeTagB.getFinalCands();
    std::list<BtaCandidate*>::const_iterator stableDaughterSEIter;
    HepAList<BtaCandidate> daughterSE;
    for( stableDaughterSEIter = stableDaughtersSE.begin(); stableDaughterSEIter != stableDaughtersSE.end(); 
	 stableDaughterSEIter++) {
      daughterSE.append( **stableDaughterSEIter );
    }
    BtaThrustFitter SEthr;
    SEthr.boostAndCompute(daughterSE, upsBooster); 

    // Thrust of the rest of the event: the ChargedTracks and CalorNeutral not used in the Breco
    HepAList<BtaCandidate> restOfEvent;
    BtaCandidate *track(0);
    while (track = (*iterAllTracks)())
      if (!track->overlaps(*aBReco)) restOfEvent.append(track);
    iterAllTracks->rewind();    
    BtaCandidate *neutral(0);
    while (neutral = (*iterAllPhotons)())
      if (!neutral->overlaps(*aBReco)) restOfEvent.append(neutral);
    iterAllPhotons->rewind();
    BtaThrustFitter ROEthr;
    ROEthr.boostAndCompute(restOfEvent, upsBooster); 
    Hep3Vector thrustAxisOfROE = ROEthr.thrustAxis();

    candCosT = cos((SEthr.thrustAxis()).angle(thrustAxisOfROE));
    tempCosT = candCosT;
    delete upsBooster;


    //========================  Storing the values of the best candidate so far  =======================
    theBestUpsilon = anUpsilon;
    theBestBReco = aBReco;
    theBestLepton = aLepton;
    theBestD = aD;
    theBestBareD = bareD;
    theBestDP4 = aD->p4();
    theBestMissP4 = pMiss;

    candLepIsPi = 0;
    if(isPion(aLepton)) candLepIsPi += 1;
    if(isKaon(aLepton)) candLepIsPi += 10;
    candType = tempType;
    candDstarType= tempDstarType;
    candDType = tempDType;
    candIsMixed = tempIsMixed;
    candExtraTracks = tempExtraTracks;
    eventM2 = pMissExtraTracks.mag2();
    candIsMu = (bool)(abs((int)(aLepton->pdtEntry()->lundId()))==(int)PdtLund::mu_minus);
    candLepCharge = (int)(aLepton->charge());
    candDLund = (int)(aD->pdtEntry()->lundId());
    candLepTru = (int)(trueLepton != 0);
    if (trueLepton)
      if (trueLepton->theMother())
	trueLepMother = (int)trueLepton->theMother()->pdtEntry()->lundId();
    candDTru = (int)(trueD != 0);
    candBTru = (int)(trueB != 0);
    candDmass = tempDmass;
    candDeltam = tempDeltam;
    candFitStat = tempFitStat;
    candFitNdof = tempFitNdof;
    candFitChi2 = tempFitChi2;
    candFitProb = probab(tempFitNdof,tempFitChi2);
    if (candDstarType > 0) {
      candPisoftP = tempPisoftP;
      candPisoftTheta = tempPisoftTheta;
      candPisoftPhi = tempPisoftPhi;
      candPisoftE1 = tempPisoftE1;
      candPisoftE2 = tempPisoftE2;
      candPisoftM2 = pMissPisoft.mag2();
    } else {
      candPisoftP = candPisoftTheta = candPisoftPhi = candPisoftE1 = candPisoftE2 = candPisoftM2 = -1.;
    }
    candKsMass = tempKsMass;
    candPi0Mass = tempPi0Mass;
    candPi0P = tempPi0P;
    candPLep = aLepton->p();
    candThetaLep = aLepton->p4().theta();
    candPstarLep = pLepCM.rho();
    candEExtra = tempEExtra;
    candExtraPhotons = tempExtraPhotons;
    candUsedEnergy = tempUsedEnergy;
    candUsedPhotons = tempUsedPhotons;
    candMExtra = tempMExtra;
    candMTot = tempMTot;
    candBMode = _decayMode;
    candMES = tempMES;
    candMES_orig = _mES;
    candDeltaE = tempDeltaE;
    candIntPur = _TruePurity[SeedMode][YMode];
    candM2 = tempM2;
    candPMiss = pMiss.vect().mag();
    candEMiss = pMiss.t();
    candThetaMiss = pMiss.theta();
    candPhiMiss = pMiss.phi();
    candThetaDL = aD->p4().angle(aLepton->p3());
    candPD = aD->p();
    candPstarD = pDCM.rho();
    candEstarD = pDCM.e();
    candQ2 = qvec.mag2();
    candIsBrem = tempIsBrem;
    candBremE = tempBremE;
    candBremAngle = tempBremAngle;
    candHelicity = tempHelicity;
    candBntCha = tempBntCha;
    candBntNeu = tempBntNeu;
    candDntCha = tempDntCha;
    candDntNeu = tempDntNeu;
    candBnCharged = tempBnCharged;
    candBnNeutral = tempBnNeutral;
    candDnCharged = tempDnCharged;
    candDnNeutral = tempDnNeutral;
    candRejectedTracks = tempRejectedTracks;
    candRejectedPhotons = tempRejectedPhotons;
    candMvaBestB = tempMvaBestB;
    candBTagDmass = tempBTagDmass;
    candBTagYMass = tempBTagYmass;
    candBTagDeltam = tempBTagDeltam;
    candTagChargedMult = tempTagChargedMult;
    candTagNeutMult = tempTagNeutMult;

    tempRejectedPhotons -= tempExtraPhotons;
    if(candType==1){
      candMvaComb = readCombD0->EvaluateMVA(method);
      candMvaDl = readDlD0->EvaluateMVA(method);
    }
    if(candType==2){
      candMvaComb = readCombDs0->EvaluateMVA(method);
      candMvaDl = readDlDs0->EvaluateMVA(method);
    }
    if(candType==3){
      candMvaComb = readCombDp->EvaluateMVA(method);
      candMvaDl = readDlDp->EvaluateMVA(method);
    }
    if(candType==4){
      candMvaComb = readCombDsp->EvaluateMVA(method);
      candMvaDl = readDlDsp->EvaluateMVA(method);
    }
    if(_Verbose.value()) cout<<"manuelf: Storing Type "<<tempType<<" and isMu "<<candIsMu<<endl;

    if (fitter) {delete fitter; delete bareD; if (tempDstarType) delete aD;}
  } //iterUpsilon
  if(_Verbose.value()) cout<<"manuelf: End of iterUpsilon"<<endl;

  if (candType>-1 && candExtraTracks==0)
    effhist[candType-1] -> accumulate(MCType,7);

  // Getting the variables from MCTagger
  if (candType > -1) {
    tag -> getInt(MCSubmode  ,"MCSubmode");
    tag -> getInt(MCCombmode ,"MCCombmode");
    tag -> getInt(MCPions    ,"MCPions");
    tag -> getInt(MCTaumode  ,"MCTaumode");
    tag -> getInt(MCBrem     ,"MCBrem");
    tag -> getInt(MCScatter  ,"MCScatter");
    tag -> getInt(MCD        ,"MCD");
    tag -> getInt(MC2Body    ,"MC2Body");
    tag -> getInt(MCDoubleSL ,"MCDoubleSL");
    tag -> getInt(nTrueP     ,"nTrueP");
    tag -> getInt(nTrueKL    ,"nTrueKL");
    tag -> getInt(nTrueNu    ,"nTrueNu");
    tag -> getInt(isVub      ,"isVub");
    tag -> getInt(uidBTag    ,"uidBTag");
    tag -> getInt(uidD       ,"uidD");
    tag -> getInt(uidLep     ,"uidLep");
    tag -> getInt(uidNu      ,"uidNu");
    tag -> getInt(uidPi0     ,"uidPi0");
    tag -> getInt(uidDssPi0  ,"uidDssPi0");
    tag -> getFloat(truePLep ,"truePLep");
    tag -> getFloat(truePD   ,"truePD");
    tag -> getFloat(trueDmass,"trueDmass");
    tag -> getFloat(trueCTL  ,"trueCTL");
    tag -> getFloat(trueCTV  ,"trueCTV");
    tag -> getFloat(trueChi  ,"trueChi");
    tag -> getFloat(trueQ2   ,"trueQ2");
    tag -> getInt(trueLepCharge,"trueLepCharge");
    tag -> getFloat(trueTauFlight,"trueTauFlight");
    tag -> getInt(MCDssmode  ,"MCDssmode");
    tag -> getInt(MCUnsim    ,"MCUnsim");
    // manuelf variables
    tag -> getInt  (nMufromPi      ,  "nMufromPi");
    tag -> getInt  (ntrueLep       ,  "ntrueLep");
    tag -> getInt  (trueMotherLep1 ,  "trueMotherLep1");
    tag -> getInt  (trueMotherLep2 ,  "trueMotherLep2");
    tag -> getInt  (trueMotherLep3 ,  "trueMotherLep3");
    tag -> getInt  (trueLundLep1   ,  "trueLundLep1");
    tag -> getInt  (trueLundLep2   ,  "trueLundLep2");
    tag -> getInt  (trueLundLep3   ,  "trueLundLep3");
    tag -> getFloat(truePLep1      ,  "truePLep1");
    tag -> getFloat(truePLep2      ,  "truePLep2");
    tag -> getFloat(truePLep3      ,  "truePLep3");
    tag -> getInt  (trueTagPi      ,  "trueTagPi");
    tag -> getInt  (trueTagK       ,  "trueTagK");
    tag -> getInt  (trueTagKs      ,  "trueTagKs");
    tag -> getInt  (trueTagPi0     ,  "trueTagPi0");
    tag -> getInt  (trueTagE       ,  "trueTagE");
    tag -> getInt  (trueTagMu      ,  "trueTagMu");
    tag -> getInt  (trueTagDs      ,  "trueTagDs");
    tag -> getInt  (trueTagD       ,  "trueTagD");
    tag -> getInt  (trueTagYPi     ,  "trueTagYPi");
    tag -> getInt  (trueTagYK      ,  "trueTagYK");
    tag -> getInt  (trueTagYKs     ,  "trueTagYKs");
    tag -> getInt  (trueTagYPi0    ,  "trueTagYPi0");
    tag -> getInt  (trueTagElse    ,  "trueTagElse");

    // Checking if any lepton was used in the Breco
    candTagLeptons = 0;
    while (BtaCandidate* aLepton = (*iterElec)()) {
      if (aLepton->overlaps(*theBestBReco)) candTagLeptons ++;
    }
    iterElec->rewind();
    if(_Verbose.value()) cout<<"manuelf: Checking if muons were used in the BReco"<<endl;
    while (BtaCandidate* aLepton = (*iterMuon)()) {
      if (aLepton->overlaps(*theBestBReco)) candTagLeptons ++;
    }
    iterMuon->rewind();

    if(_Verbose.value()) cout<<"manuelf: Calculating Y system mass and counting of tracks and neutrals"<<endl;
    // tag-side
    // emc and trk variables have the extremes values of theta for the tracks and 
    // neutral clusters. Also, # of tracks and neutral particles is counted
    emcfwd = trkfwd = 9.0;
    emcbwd = trkbwd = -1.0;

    //
    // extra pi0 block
    // Lists are created with D**->Dpi0 and the new missing mass and energy without the pi0
    npi0 = 0;
    int theBestEpi0=-1,theBestMasspi0=-1;
    float bestEsofar=-999.9,bestMsofar=999.9;
    HepLorentzVector iP4pi0;
    while (BtaCandidate* aPi0 = (*iterSoftPi0)()) {
      if (aPi0->overlaps(*theBestUpsilon)) continue;

      // The momentum of the photons is stored in Ve1 and Ve2
      // and the total pi0 momentum in rawP4
      HepLorentzVector rawP4;
      HepAListIterator<BtaCandidate> iterDauPi0 = aPi0->daughterIterator();
      HepLorentzVector Ve1(0), Ve2(0);
      while (BtaCandidate *aDauPi0 = iterDauPi0()) {
	rawP4 += aDauPi0->p4();
	if(Ve1==0)
	  Ve1 = aDauPi0->p4();
	else
	  Ve2 = aDauPi0->p4();
      }
      iterDauPi0.rewind();
      Hep3Vector pi0CMboost = aPi0->p4().boostVector();
      HepLorentzVector maxEp(0);
      if(Ve2.e() < Ve1.e()){
	e1pi0.append(Ve1.e());
	e2pi0.append(Ve2.e());
	maxEp = Ve1;
      } else {
	e1pi0.append(Ve2.e());
	e2pi0.append(Ve1.e());
	maxEp = Ve2;
      }
      maxEp.boost(-pi0CMboost);
      cosEpi0.append(cos(maxEp.angle(aPi0->p4().vect())));
      mpi0.append(rawP4.mag());
      ppi0.append(aPi0->p());
      // Creation of a D** with the D and the current pi0
      HepLorentzVector dss = theBestDP4;
      dss += aPi0->p4();
      HepLorentzVector BsigP4 = Upsilon - theBestBReco->p4();
      HepLorentzVector qpi0vec = BsigP4 - dss;
      q2pi0.append(qpi0vec.mag2());
      dmpi0.append(dss.mag() - theBestDP4.mag());
      HepLorentzVector pi0miss = theBestMissP4;
      pi0miss -= aPi0->p4();
      mm2pi0.append(pi0miss.mag2());
      pmisspi0.append(pi0miss.vect().mag());
      eextrapi0.append(candEExtra - aPi0->energy());
      // The pi0 with the highest energy and mass closest to the PDG value are selected
      if (aPi0->energy() > bestEsofar) {
	theBestEpi0 = npi0;
	bestEsofar = aPi0->energy();
	iP4pi0 = aPi0->p4();
      }
      if (fabs(rawP4.mag() - 0.135) < bestMsofar) {
	theBestMasspi0 = npi0;
	bestMsofar = fabs(rawP4.mag() - 0.135);
      }
      npi0 ++;
    }// While aPi0
    iterSoftPi0->rewind();
    for (int q = 0 ; q < npi0 ; q ++) {
      bestepi0   .append((int)(q == theBestEpi0));
      bestmasspi0.append((int)(q == theBestMasspi0));
    }
    if(npi0>0){
	impi0 = mpi0[theBestEpi0];
	ie1pi0 = e1pi0[theBestEpi0];
	idmpi0 = dmpi0[theBestEpi0];
	ieextrapi0 = eextrapi0[theBestEpi0];
	ippi0 = ppi0[theBestEpi0];
	imm2pi0 = mm2pi0[theBestEpi0];
	ipmisspi0 = pmisspi0[theBestEpi0];
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
    //
    // extra pi0pi0 block
    // The same but with D** -> D pi0 pi0
    ndpi0 = 0;
    int theBestdpi0=-1;
    bestEsofar=-999.9;
    //mm2dpi0,eextradpi0,dmdpi0,pmissdpi0,bestdpi0;
    while (BtaCandidate* aPi0 = (*iterSoftPi0)()) {
      if (aPi0->overlaps(*theBestUpsilon)) continue;
      iterSoftPi02->rewind();
      while (BtaCandidate *aPi02 = (*iterSoftPi02)()) {
	if (aPi02->overlaps(*aPi0)) continue;
	if (aPi02->overlaps(*theBestUpsilon)) continue;
	HepLorentzVector rawP41,rawP42;
	HepLorentzVector dss = theBestDP4;
	dss += aPi0->p4()+aPi02->p4();
	dmdpi0.append(dss.mag() - theBestDP4.mag());
	HepLorentzVector pi0miss = theBestMissP4;
	pi0miss -= (aPi0->p4() + aPi02->p4());
	mm2dpi0.append(pi0miss.mag2());
	pmissdpi0.append(pi0miss.vect().mag());
	eextradpi0.append(candEExtra - aPi0->energy());
	p1dpi0.append(aPi0->p());
	p2dpi0.append(aPi02->p());
	HepAListIterator<BtaCandidate> iterDauPi01 = aPi0->daughterIterator();
	while (BtaCandidate *aDauPi0 = iterDauPi01())
	  rawP41 += aDauPi0->p4();
	iterDauPi01.rewind();
	m1dpi0.append(rawP41.mag());
	HepAListIterator<BtaCandidate> iterDauPi02 = aPi02->daughterIterator();
	while (BtaCandidate *aDauPi0 = iterDauPi02())
	  rawP42 += aDauPi0->p4();
	iterDauPi02.rewind();
	m2dpi0.append(rawP42.mag());
	if (aPi0->energy() + aPi02->energy() > bestEsofar) {
	  theBestdpi0 = ndpi0;
	  bestEsofar = aPi0->energy() + aPi02->energy();
	}
	ndpi0 ++;
      }
    }
    iterSoftPi0->rewind();
    for (int q = 0 ; q < ndpi0 ; q ++)
      bestdpi0.append((int)(q == theBestdpi0));


    //
    // extra Ks block
    //
//     nks = 0;
//     while (BtaCandidate* aKs = (*iterKs)()) {
//       if (aKs->overlaps(*theBestUpsilon)) continue;
//       mks.append(aKs->mass());
//       pks.append(aKs->p());
//       vxks.append(aKs->decayVtx()->point().x());
//       vyks.append(aKs->decayVtx()->point().y());
//       vzks.append(aKs->decayVtx()->point().z());
//       HepLorentzVector ksmiss = theBestMissP4;
//       ksmiss -= aKs->p4();
//       mm2ks.append(ksmiss.mag2());
//       int kset = 0;
//       while (BtaCandidate* ttrk = (*iterTracks)()) {
// 	if (ttrk->overlaps(*theBestUpsilon)) continue;
// 	if (ttrk->overlaps(*aKs)) continue;
// 	kset ++;
//       }
//       iterTracks->rewind();
//       extratracksks.append(kset);
//       nks ++;
//     }
//     iterKs->rewind();

    //
    // Recalculate mm2
    //

    candM2Tru = candM2TruBTag = candM2TruD = candM2TruLep = -99.;
    candBTagDeltaP = candBTagDeltaPx = candBTagDeltaPy = candBTagDeltaPz = candBTagDeltaE = -999;
    candDDeltaP = candDDeltaPx = candDDeltaPy = candDDeltaPz = candDDeltaE = -999;
    candLepDeltaP = candLepDeltaPx = candLepDeltaPy = candLepDeltaPz = candLepDeltaE = -999;
    candNuDeltaP = candNuDeltaPx = candNuDeltaPy = candNuDeltaPz = candNuDeltaE = -999;
    candBDPx = candBDPy = candBDPz = candBDE = -999;
    candMissPx = candMissPy = candMissPz = candMissE = -999;
    float truPPi0 = -999;
    float truDssPPi0 = -999;

    // Storing of the true particles
    BtaCandidate *truBTag(0),*truD(0),*truLep(0),*truNu(0),*aTru,*truUpsilon(0);
    if (MCTruthList)
      while (aTru = (*iterMCTruth)()) {
 	if ((int)aTru->pdtEntry()->lundId() == PdtLund::Upsilon_4S) {
	  truUpsilon = aTru;
// 	  if(MCType==0 && npi0>0 && abs(imm2pi0)<.5 && impi0>.125 && impi0<.145 && ieextrapi0<.5 &&
// 	     ippi0>.4 && ipmisspi0>.2 && candIsMu==1) {
// 	    printTree -> removeTerminal(PdtLund::K_plus);
// 	    printTree -> removeTerminal(PdtLund::K_minus);
// 	    printTree -> removeTerminal(PdtLund::pi_plus);
// 	    printTree -> removeTerminal(PdtLund::pi_minus);
// 	    printTree -> setOutputFunc(BtaPrintTree::printNameP3);
// 	    cout<<endl<<"Peaking combinatoric: candType "<<candType<<", mm2pi0 "<<imm2pi0<<", candM2 ";
// 	    cout<<candM2<<", pmisspi0 "<<ipmisspi0<<", ppi0 "<<ippi0<<", mass D** "<<idmpi0+candDmass<<endl;
// 	    cout<<"Reconstruted lepton: ("<<aLepton->p4().x()<<", "<<aLepton->p4().y()<<", "<<aLepton->p4().z()<<", ";
// 	    cout<<aLepton->p4().t()<<").   candLepTru: "<<candLepTru<<endl;
// 	    cout<<"Added pi0: ("<<iP4pi0.x()<<", "<<iP4pi0.y()<<", "<<iP4pi0.z()<<", ";
// 	    cout<<iP4pi0.t()<<") "<<endl;
// 	    cout<<"MC tree"<<endl;
// 	    cout<<printTree->print(*aTru)<<endl;
// 	    cout<<"Reconstructed tree"<<endl;
// 	    cout<<printTree->print(*theBestUpsilon)<<endl;
// 	    printTree -> addTerminal(PdtLund::K_plus);
// 	    printTree -> addTerminal(PdtLund::K_minus);
// 	    printTree -> addTerminal(PdtLund::pi_plus);
// 	    printTree -> addTerminal(PdtLund::pi_minus);
// 	    printTree -> setOutputFunc(BtaPrintTree::printName);
// 	    cout<<endl<<"MC tree without momenta"<<endl;
// 	    cout<<printTree->print(*aTru)<<endl;
// 	    cout<<"Reconstructed tree"<<endl;
// 	    cout<<printTree->print(*theBestUpsilon)<<endl<<endl;
// 	  }
	}
	if ((int)aTru->uid() == uidPi0) truPPi0 = aTru->p();
	if ((int)aTru->uid() == uidDssPi0) truDssPPi0 = aTru->p();
	if ((int)aTru->uid() == uidBTag) truBTag = aTru;
	if ((int)aTru->uid() == uidD   ) truD    = aTru;
	if ((int)aTru->uid() == uidLep ) truLep  = aTru;
	if ((int)aTru->uid() == uidNu  ) truNu   = aTru;
      }
    iterMCTruth->rewind();

    if (truUpsilon && truBTag && truD && truLep) {

      HepLorentzVector pmisstru = truUpsilon->p4()-truBTag->p4()-truD->p4()-truLep->p4();
      candM2Tru = pmisstru.mag2();

      HepLorentzVector pmisstru1 = theBestMissP4;
      pmisstru1 += theBestBReco->p4() - truBTag->p4();
      candM2TruBTag = pmisstru1.mag2();

      HepLorentzVector brecoCM = theBestBReco->p4();
      brecoCM.boost(-CMboost);
      HepLorentzVector trueBCM = truBTag->p4();
      trueBCM.boost(-CMboost);
      candBTagDeltaP = (brecoCM.vect() - trueBCM.vect()).mag();
      candBTagDeltaPx = brecoCM.px() - trueBCM.px();
      candBTagDeltaPy = brecoCM.py() - trueBCM.py();
      candBTagDeltaPz = brecoCM.pz() - trueBCM.pz();
      candBTagDeltaE  = brecoCM.e()  - trueBCM.e();
      HepLorentzVector pmisstru2 = theBestMissP4;
      pmisstru2 += theBestDP4 - truD->p4();
      candM2TruD = pmisstru2.mag2();

      HepLorentzVector theDCM = theBestDP4;
      theDCM.boost(-CMboost);
      HepLorentzVector truDCM = truD->p4();
      truDCM.boost(-CMboost);
      candDDeltaP = (theDCM.vect() - truDCM.vect()).mag();
      candDDeltaPx = theDCM.px() - truDCM.px();
      candDDeltaPy = theDCM.py() - truDCM.py();
      candDDeltaPz = theDCM.pz() - truDCM.pz();
      candDDeltaE  = theDCM.e()  - truDCM.e();
      HepLorentzVector pmisstru3 = theBestMissP4;
      pmisstru3 += theBestLepton->p4() - truLep->p4();
      candM2TruLep = pmisstru3.mag2();

      HepLorentzVector theLepCM = theBestLepton->p4();
      theLepCM.boost(-CMboost);
      HepLorentzVector truLepCM = truLep->p4();
      truLepCM.boost(-CMboost);
      candLepDeltaP = (theLepCM.vect() - truLepCM.vect()).mag();
      candLepDeltaPx = theLepCM.px() - truLepCM.px();
      candLepDeltaPy = theLepCM.py() - truLepCM.py();
      candLepDeltaPz = theLepCM.pz() - truLepCM.pz();
      candLepDeltaE  = theLepCM.e()  - truLepCM.e();

      HepLorentzVector theBDCM = brecoCM + theDCM;
      candBDPx = theBDCM.px();
      candBDPy = theBDCM.py();
      candBDPz = theBDCM.pz();
      candBDE  = theBDCM.e();

      HepLorentzVector pnuCM = theBestMissP4;
      pnuCM.boost(-CMboost);
      candMissPx = pnuCM.px();
      candMissPy = pnuCM.py();
      candMissPz = pnuCM.pz();
      candMissE  = pnuCM.e();
    }
//     if (truNu) {
//       HepLorentzVector pnuCM = theBestMissP4;
//       pnuCM.boost(-CMboost);
//       HepLorentzVector pnutruCM = truNu->p4();
//       pnutruCM.boost(-CMboost);
//       candNuDeltaP = (pnuCM.vect() - pnutruCM.vect()).mag();
//       candNuDeltaPx = pnuCM.px() - pnutruCM.px();
//       candNuDeltaPy = pnuCM.py() - pnutruCM.py();
//       candNuDeltaPz = pnuCM.pz() - pnutruCM.pz();
//       candNuDeltaE  = pnuCM.e()  - pnutruCM.e();
//     }

//     if (truBTag && truD && truLep && truNu) {
//       HepLorentzVector pCM3Tru = truBTag->p4() + truD->p4() + truLep->p4();
//       HepLorentzVector pCM4Tru = pCM3Tru + truNu->p4();
//       pCM3Tru.boost(-CMboost);
//       pCM4Tru.boost(-CMboost);

//       HepLorentzVector pCM3 = theBestBReco->p4() + theBestDP4 + theBestLepton->p4();
//       HepLorentzVector pCM4 = pCM3 + theBestMissP4;
//       pCM3.boost(-CMboost);
//       pCM4.boost(-CMboost);
//     }

    // Calculation of the B and Ds modes for combinatoric event in Monte Carlo
    MCCombB = MCCombDs = 0;
    if (MCType == 0 && truUpsilon) {
      BtaCandidate *B1(0), *B2(0);
      HepAListIterator<BtaCandidate> iterTrueB = truUpsilon->daughterIterator();
      B1 = iterTrueB();
      B2 = iterTrueB();
      iterTrueB.rewind();
      if (B1->p4().angle(theBestBReco->p3()) > B2->p4().angle(theBestBReco->p3()))
	MCComblong = combMode(B1,MCCombB,MCCombDs,MCComblongD);
      else MCComblong = combMode(B2,MCCombB,MCCombDs,MCComblongD);
    }

    //
    //systematic error related calculations
    //
    nGTLHard=nGTLSoft=nOtherTracks = 0;
    trkWeight=kaonWeight=pionWeight=elecWeight=muonWeight = 1.0;
    candKSdxy=candKSTheta=candKSPT=candPPi0 = -999.9;
    nCandK = nCandPi = 0;

    if(_Verbose.value()) cout<<"manuelf: systematic error related calculations"<<endl;     
    if (MCType > -1) {
      std::vector<double> trkVect(4,0);
      double trkCorr(0),trkErr(0);
      bool ok=false;
      trkVect[3] = GTVLList->length();

      //lepton
      BtaCandidate *tempLep = theBestLepton;
      if (candIsBrem) {
	HepAListIterator<BtaCandidate> iterDauLep = theBestLepton->daughterIterator();
	while (BtaCandidate *aDauL = iterDauLep()) {
	  if (abs((int)aDauL->pdtEntry()->lundId()) == (int)PdtLund::e_minus) tempLep = aDauL;
	}
	iterDauLep.rewind();
      }
      trkVect[0] = tempLep->pt();
      trkVect[1] = tempLep->p4().theta();
      trkVect[2] = tempLep->p4().phi();
      if (tableGTL) ok = tableGTL->get_efficiency_and_error(trkCorr,trkErr,trkVect);
      if (!ok) {
	cout << "######### problem with lepton trkEff: " << tempLep->p4() << endl;
      }
      trkWeight *= trkCorr;
      nGTLHard ++;
      if(_Verbose.value()) cout<<"manuelf: Here some e("<<elecMap
			       <<") and mu("<<muonMap<<") are accessed at "<<tempLep->uid()<<endl;
      if (candIsMu) {
	if(muonMap){
	  PidWeight pw = (*muonMap)[tempLep->uid()];
	  pidhists[1]->accumulate(pw._status);
	  if (pw._status==PidWeight::ok) muonWeight = pw._value;
	}
      } else {
	if(_Verbose.value()) cout<<"manuelf: elecMap accessed "<<elecMap<<endl;
	if(elecMap){
	  PidWeight pw = (*elecMap)[tempLep->uid()];
	  pidhists[0]->accumulate(pw._status);
	  if (pw._status==PidWeight::ok) elecWeight = pw._value;
	}
      }

      if(_Verbose.value()) cout<<"manuelf: Now onto the D trkEff"<<endl;     
      //daughters of the D
      HepAListIterator<BtaCandidate> iterBareD = theBestBareD->daughterIterator();
      while (aDauD = iterBareD()) {
	if (aDauD->pdtEntry()->lundId() == (int)PdtLund::K_S0) {
	  if(_Verbose.value()) cout<<"manuelf: Looking for Kshorts"<<endl;     
	  const BbrLorentzVectorErr & ip = gblEnv->getBta()->pepBeams()->interactionPoint();
	  HepPoint beamSpotPoint(ip.x(), ip.y(), ip.z());
	  candKSdxy = aDauD->docaXY(beamSpotPoint);
	  candKSTheta = aDauD->p4().theta();
	  candKSPT = aDauD->pt();
	} //KS
	else if (aDauD->pdtEntry()->lundId() == (int)PdtLund::pi0) {
	  candPPi0 = aDauD->p();
	} //pi0
	else {
	  if(_Verbose.value()) cout<<"manuelf: Looking into the tracks"<<endl;     
	  if (isGTL(aDauD)) {
	    trkVect[0] = aDauD->pt();
	    trkVect[1] = aDauD->p4().theta();
	    trkVect[2] = aDauD->p4().phi();
	    if (tableGTL) ok = tableGTL->get_efficiency_and_error(trkCorr,trkErr,trkVect);
	    if (!ok) {
	      cout << "######### problem with D trkEff: " << aDauD->p4() << endl;
	    }
	    trkWeight *= trkCorr;
	    if (aDauD->pt() > 0.2) nGTLHard ++;
	    else nGTLSoft ++;
	  } else {
	    if (isGTVL(aDauD)) trkWeight *= 0.995;
	    else trkWeight *= 0.9975;
	    nOtherTracks ++;
	  }
	  if(_Verbose.value()) cout<<"manuelf: Here some kaon("<<kaonMap<<") and pion("
				   <<pionMap<<") are accessed at " <<aDauD->uid()<<endl;
	  if (abs((int)aDauD->pdtEntry()->lundId()) == (int)PdtLund::K_plus) {
	    candPKaons.append(aDauD->p());
	    nCandK ++;
	    if(kaonMap){
	      PidWeight pw = (*kaonMap)[aDauD->uid()];
	      if (pw._status==PidWeight::ok) kaonWeight *= pw._value;
	      pidhists[2]->accumulate(pw._status);
	    }
	  } else  {
	    candPPions.append(aDauD->p());
	    nCandPi ++;
	    if(pionMap){
	      PidWeight pw = (*pionMap)[aDauD->uid()];
	      pidhists[3]->accumulate(pw._status);
	      if (pw._status==PidWeight::ok) pionWeight *= pw._value;
	    }
	  }
	} //k/pi
      }// While aDauD
      iterBareD.rewind();

      //pisoft
      if (candDstarType) {
	HepAListIterator<BtaCandidate> iterDauDstar = theBestD->daughterIterator();
	while (aDauD = iterDauDstar()) {
	  if (abs((int)aDauD->pdtEntry()->lundId()) == (int)PdtLund::pi_plus) {
	    nOtherTracks ++;
	  }
	}
	iterDauDstar.rewind();
      }
    } //if MC

    wBF = myWM->getEventWeight(candType,candDstarType,MCType,MCSubmode,MCDssmode,MCD,MCPions,
			       isBzero,0,MCUnsim,truPPi0,truDssPPi0,trueCTL,trueCTV,trueChi,trueQ2,
			       trueLepCharge,candM2);
    wFF = myWMFF->getEventWeight(candType,candDstarType,MCType,MCSubmode,MCDssmode,MCD,MCPions,
			       isBzero,0,MCUnsim,truPPi0,truDssPPi0,trueCTL,trueCTV,trueChi,trueQ2,
			       trueLepCharge,candM2);
    wComb = myWM->getCombWeight(MCCombB,MCCombDs,MCDoubleSL,candLepTru);

    if(_Verbose.value()) cout<<"manuelf: Start of the columning"<<endl;
    cands->column("runnum"   , (int)_eventID->run());
    cands->column("platform" , (int)_eventID->eventIdTriplet().platform());
    cands->column("partition", (int)_eventID->eventIdTriplet().partitionMask());
    cands->column("upperID"  , (int)_eventID->eventIdTriplet().timeStamp().binary().upper);
    cands->column("lowerID"  , (int)_eventID->eventIdTriplet().timeStamp().binary().lower);
    cands->column("isSP6"    , (int)(_eventID->condKeyTriplet().key() > *run4Start));
    cands->column("MCType",MCType);             // Type of SigB decay
    cands->column("MCSubmode",MCSubmode);       // True D mode from the SigB
    cands->column("MCCombmode",MCCombmode);     // Type of background
    cands->column("MCPions",MCPions);           // Number of pions from SigB. +10 if pi0
    cands->column("MCTaumode",MCTaumode);       // -1 no tau,0: tau no lepton,1 tau e, 2 tau mu
    cands->column("MCBrem",MCBrem);             // 1 if gamma from SigB or tau
    cands->column("MCScatter",MCScatter);       // 1 if final lepton has daughters
    cands->column("MCD",MCD);                   // Lund of the D from SigB
    cands->column("MC2Body",MC2Body);  // Last digit is: 1 pi, 2 rho, 3 K+, 4 K*+ or the same
                                       // format as second to last: 1 D, 2 D*, 3 D**, 4 Ds, 5 Ds*, 
                                       // 6 Ds** if bigger than 100. If >1000, 2 two body B's 
    cands->column("MCDoubleSL",MCDoubleSL);     // 1 if both B's are semileptonic
    cands->column("MCComblong",MCComblong);    // Indicates the signal B decay for combinatoric events: lepton/charmless
    cands->column("MCComblongD",MCComblongD);  // Indicates the signal B decay for combinatoric events: Charmed mesons
    cands->column("MCCombB",MCCombB);  // First digit: 1 pi+, 2 rho, 3 K+, 4 K0, 5 K*+, 9 more than
                                       // one charmless. 2nd and 3rd: 1 D0, 2 D+, 3 D*0, 4 D*+,
                                       // 5 D**0, 6 D**+, 7 Ds, 8 Ds*, 9 Ds** for the true SigB and
                                       // only for MCType == 0
    cands->column("MCCombDs",MCCombDs);// Ds into: 1 tau, 2 mu, 3 e, 4 phi l, 5 eta l, 6 eta' l
    cands->column("nTrueP",nTrueP);             // Number of protons
    cands->column("nTrueKL",nTrueKL);           // Number of Klongs
    cands->column("nTrueNu",nTrueNu);           // Number of neutrinos
    //     cands->column("isVub",isVub);               // 1 if there is a lepton
    cands->column("isBzero",isBzero);           // 1 if Y(4S) -> B0 B0bar
    cands->column("MCDssmode",MCDssmode);       // Dss daughter: 1 D0; 2 D*0; 3 D+; 4 D*+; +4 if MCPions==1
    cands->column("MCUnsim",MCUnsim);           // 1 if it is an unsimulated event
    cands->column("truePLep",truePLep);         // Momentum of the true lepton when it comes from SL B
    cands->column("trueLepMother",trueLepMother); // Lund ID of the true lepton mother and it comes from SL B
    cands->column("truePD",truePD);             // Momentum of the true D when it comes from SL B
    cands->column("trueDmass",trueDmass);       // Mass of the true D when it comes from SL B
    cands->column("trueCTL",trueCTL);           // True angle between l and the B meson, theta_l
    cands->column("trueCTV",trueCTV);           // True angle between D and the B meson, theta_V
    cands->column("trueChi",trueChi);           // True angle between D and the W planes, chi
    cands->column("trueQ2",trueQ2);             // True momentum transfer
    cands->column("trueLepCharge",trueLepCharge);// True charge of the lepton
    cands->column("trueTauFlight",trueTauFlight);// True tau flight
    cands->column("truePPi0",truPPi0);          // True pi0 momentum coming from a D*
    cands->column("trueDssPPi0",truDssPPi0);    // True pi0 momentum coming from a D**

    cands->column("nMufromPi"     ,nMufromPi);      // Number of muons whose mother is a pion
    cands->column("ntrueLep"      ,ntrueLep);      // Total number of true Leptons
    cands->column("trueMotherLep1",trueMotherLep1);// Lund of mother of lepton 1
    cands->column("trueMotherLep2",trueMotherLep2);// Lund of mother of lepton 3
    cands->column("trueMotherLep3",trueMotherLep3);// Lund of mother of lepton 3
    cands->column("trueLundLep1"  ,trueLundLep1);  // Lund of lepton 1
    cands->column("trueLundLep2"  ,trueLundLep2);  // Lund of lepton 2
    cands->column("trueLundLep3"  ,trueLundLep3);  // Lund of lepton 3
    cands->column("truePLep1"     ,truePLep1);     // Momentum of lepton 1
    cands->column("truePLep2"     ,truePLep2);     // Momentum of lepton 2
    cands->column("truePLep3"     ,truePLep3);     // Momentum of lepton 3


    cands->column("trueTagPi"     ,trueTagPi);     // # pi coming from the D in the tag B
    cands->column("trueTagK"      ,trueTagK);      // # K coming from the D in the tag B
    cands->column("trueTagKs"     ,trueTagKs);     // # Ks coming from the D in the tag B
    cands->column("trueTagPi0"    ,trueTagPi0);    // # pi0 coming from the D in the tag B
    cands->column("trueTagE"      ,trueTagE);      // # e coming from the D in the tag B
    cands->column("trueTagMu"     ,trueTagMu);     // # mu coming from the D in the tag B
    cands->column("trueTagDs"     ,trueTagDs);     // Lund of the D* or D in the tag B
    cands->column("trueTagD"      ,trueTagD);      // Lund of the D in the tag B
    cands->column("trueTagYPi"    ,trueTagYPi);    // # pi not coming from the D in the tag B
    cands->column("trueTagYK"     ,trueTagYK);     // # K not coming from the D in the tag B
    cands->column("trueTagYKs"    ,trueTagYKs);    // # Ks not coming from the D in the tag B
    cands->column("trueTagYPi0"   ,trueTagYPi0);   // # pi0 not coming from the D in the tag B
    cands->column("trueTagElse"   ,trueTagElse);   // 1 in the tag B decays only to pi/K/l
  
    cands->column("nUpsilon",UpsilonList->length()); // Number of upsilon candidates
    cands->column("nups1",nups1);               // Number of upsilon candidates with a D0
    cands->column("nups2",nups2);               // Number of upsilon candidates with a D*0
    cands->column("nups3",nups3);               // Number of upsilon candidates with a D+
    cands->column("nups4",nups4);               // Number of upsilon candidates with a D*+
    cands->column("ncands",ncands);             // Number of upsilon candidates after all skim cuts
    cands->column("ncands1",ncands1);           // Number of upsilon candidates after all skim cuts: D0
    cands->column("ncands2",ncands2);           // Number of upsilon candidates after all skim cuts: D*0
    cands->column("ncands3",ncands3);           // Number of upsilon candidates after all skim cuts: D+
    cands->column("ncands4",ncands4);           // Number of upsilon candidates after all skim cuts: D*+
    cands->column("candType",candType);         // 1 D0; 2 D*0; 3 D+; 4 D*+
    cands->column("candDstarType",candDstarType);// 0 D0/D+, 1 D*0->D0pi0/D*+->D0pi, 2 D*0->D0gamma/D*+->Dpi0
    cands->column("candDType",candDType);       // See 
    cands->column("candIsMixed",candIsMixed);   // 1 if the Lunds of BReco and sigD have different sign
    cands->column("candExtraTracks",candExtraTracks); // Number of not used tracks from RecoilGoodTracks
    cands->column("candRejectedTracks",candRejectedTracks); // Number of not used tracks from ChargedTracks
    cands->column("candRejectedPhotons",candRejectedPhotons-candExtraPhotons); // Difference between the not used
                                                           // photons of CalorNeutral list and RecoilGoodPhotons
    cands->column("wBF",wBF);                    // Weight due to the semileptonic BF
    cands->column("wcomb",wComb);                // Weight due to the combinatoric BF
    cands->column("wFF",wFF);                    // Weight due to the Form Factors
    cands->column("eventM2",eventM2);            // Mass squared of the momentum of the extra tracks
    cands->column("candIsMu",candIsMu);
    cands->column("candLepIsPi",candLepIsPi);
    cands->column("candLepCharge",candLepCharge);
    cands->column("candDLund",candDLund);
    cands->column("candLepTru",candLepTru);
    cands->column("candDTru",candDTru);
    cands->column("candBTru",candBTru);
    cands->column("candDmass",candDmass);
    cands->column("candDeltam",candDeltam);
    cands->column("candFitStat",candFitStat);
    cands->column("candFitNdof",candFitNdof);
    cands->column("candFitChi2",candFitChi2);
    cands->column("candFitProb",candFitProb);
    cands->column("candPisoftP",candPisoftP);
    cands->column("candPisoftTheta",candPisoftTheta);
    cands->column("candPisoftPhi",candPisoftPhi);
    cands->column("candPisoftE1",candPisoftE1);
    cands->column("candPisoftE2",candPisoftE2);
    cands->column("candPisoftM2",candPisoftM2);
    cands->column("candKsMass",candKsMass);
    cands->column("candPi0Mass",candPi0Mass);
    cands->column("candPi0P",candPi0P);
    cands->column("candPLep",candPLep);
    cands->column("candThetaLep",candThetaLep);
    cands->column("candPstarLep",candPstarLep);
    cands->column("candPD",candPD);
    cands->column("candPstarD",candPstarD);
    cands->column("candEstarD",candEstarD);
    cands->column("candIsBrem",candIsBrem);
    cands->column("candBremE",candBremE);
    cands->column("candBremAngle",candBremAngle);
    cands->column("candEExtra",candEExtra);
    cands->column("candExtraPhotons",candExtraPhotons);
    cands->column("candUsedEnergy",candUsedEnergy);
    cands->column("candUsedPhotons",candUsedPhotons);
    cands->column("candMExtra",candMExtra);
    cands->column("candMTot",candMTot);
    cands->column("candBMode",candBMode);
    cands->column("candMES",candMES);
    cands->column("candMES_orig",candMES_orig);
    cands->column("candDeltaE",candDeltaE);
    cands->column("candBTagDeltaP",candBTagDeltaP);
    cands->column("candBTagDeltaPx",candBTagDeltaPx);
    cands->column("candBTagDeltaPy",candBTagDeltaPy);
    cands->column("candBTagDeltaPz",candBTagDeltaPz);
    cands->column("candBTagDeltaE",candBTagDeltaE);
    cands->column("candBTagYMass",candBTagYMass);
    cands->column("candBTagDmass",candBTagDmass);
    cands->column("candBTagDeltam",candBTagDeltam);
    cands->column("candBntCha",candBntCha);
    cands->column("candBntNeu",candBntNeu);
    cands->column("candDntCha",candDntCha);
    cands->column("candDntNeu",candDntNeu);
    cands->column("candBnCharged",candBnCharged);
    cands->column("candBnNeutral",candBnNeutral);
    cands->column("candDnCharged",candDnCharged);
    cands->column("candDnNeutral",candDnNeutral);
    cands->column("candDDeltaP",candDDeltaP);
    cands->column("candDDeltaPx",candDDeltaPx);
    cands->column("candDDeltaPy",candDDeltaPy);
    cands->column("candDDeltaPz",candDDeltaPz);
    cands->column("candDDeltaE",candDDeltaE);
    cands->column("candLepDeltaP",candLepDeltaP);
    cands->column("candLepDeltaPx",candLepDeltaPx);
    cands->column("candLepDeltaPy",candLepDeltaPy);
    cands->column("candLepDeltaPz",candLepDeltaPz);
    cands->column("candLepDeltaE",candLepDeltaE);
    cands->column("candNuDeltaP",candNuDeltaP);
    cands->column("candNuDeltaPx",candNuDeltaPx);
    cands->column("candNuDeltaPy",candNuDeltaPy);
    cands->column("candNuDeltaPz",candNuDeltaPz);
    cands->column("candNuDeltaE",candNuDeltaE);
    cands->column("candBDPx",candBDPx);
    cands->column("candBDPy",candBDPy);
    cands->column("candBDPz",candBDPz);
    cands->column("candBDE",candBDE);
    cands->column("candMissPx",candMissPx);
    cands->column("candMissPy",candMissPy);
    cands->column("candMissPz",candMissPz);
    cands->column("candMissE",candMissE);
    cands->column("candIntPur",candIntPur);
    cands->column("candTagLeptons",candTagLeptons);
    cands->column("candTagChargedMult",candTagChargedMult);
    cands->column("candTagNeutMult",candTagNeutMult);
    cands->column("emcfwd",emcfwd);
    cands->column("emcbwd",emcbwd);
    cands->column("trkfwd",trkfwd);
    cands->column("trkbwd",trkbwd);
    cands->column("candM2",candM2);
    cands->column("candQ2",candQ2);
    cands->column("candPMiss",candPMiss);
    cands->column("candEMiss",candEMiss);
    cands->column("candThetaMiss",candThetaMiss);
    cands->column("candPhiMiss",candPhiMiss);
    cands->column("candThetaDL",candThetaDL);
    cands->column("candHelicity",candHelicity);
    cands->column("candM2Tru",candM2Tru);
    cands->column("candM2TruBTag",candM2TruBTag);
    cands->column("candM2TruD",candM2TruD);
    cands->column("candM2TruLep",candM2TruLep);
    cands->column("candCosT",candCosT);
    cands->column("candMvaBestB",candMvaBestB);
    cands->column("candMvaComb",candMvaComb);
    cands->column("candMvaDl",candMvaDl);
    cands->column("candMvaDssComb",candMvaDssComb);
    cands->column("candMvaDssDl",candMvaDssDl);
    cands->column("candMvaBestBL",candMvaBestBL,"ncands");
    cands->column("candBModeL",candBModeL,"ncands");
    cands->column("candEExtraL",candEExtraL,"ncands");
    cands->column("candBntChaL",candBntChaL,"ncands");
    cands->column("candBntNeuL",candBntNeuL,"ncands");
    cands->column("candDntChaL",candDntChaL,"ncands");
    cands->column("candDntNeuL",candDntNeuL,"ncands");
    cands->column("candTypeL",candTypeL,"ncands");
    cands->column("candLepTruL",candLepTruL,"ncands");
    cands->column("candDeltaEL",candDeltaEL,"ncands");
    cands->column("npi0",npi0);
    cands->column("cosEpi0",cosEpi0,"npi0");
    cands->column("q2pi0",q2pi0,"npi0");
    cands->column("mpi0",mpi0,"npi0");
    cands->column("ppi0",ppi0,"npi0");
    cands->column("dmpi0",dmpi0,"npi0");
    cands->column("mm2pi0",mm2pi0,"npi0");
    cands->column("pmisspi0",pmisspi0,"npi0");
    cands->column("eextrapi0",eextrapi0,"npi0");
    cands->column("e1pi0",e1pi0,"npi0");
    cands->column("e2pi0",e2pi0,"npi0");
    cands->column("bestepi0",bestepi0,"npi0");
    cands->column("bestmasspi0",bestmasspi0,"npi0");
    //     cands->column("nks",nks);
    //     cands->column("mks",mks,"nks");
    //     cands->column("pks",pks,"nks");
    //     cands->column("vxks",vxks,"nks");
    //     cands->column("vyks",vyks,"nks");
    //     cands->column("vzks",vzks,"nks");
    //     cands->column("mm2ks",mm2ks,"nks");
    //     cands->column("extratracksks",extratracksks,"nks");
    cands->column("ndpi0",ndpi0);
    cands->column("dmdpi0",dmdpi0,"ndpi0");
    cands->column("mm2dpi0",mm2dpi0,"ndpi0");
    cands->column("pmissdpi0",pmissdpi0,"ndpi0");
    cands->column("eextradpi0",eextradpi0,"ndpi0");
    cands->column("p1dpi0",p1dpi0,"ndpi0");
    cands->column("p2dpi0",p2dpi0,"ndpi0");
    cands->column("m1dpi0",m1dpi0,"ndpi0");
    cands->column("m2dpi0",m2dpi0,"ndpi0");
    cands->column("bestdpi0",bestdpi0,"ndpi0");
    cands->column("nGTLHard",nGTLHard);
    cands->column("nGTLSoft",nGTLSoft);
    cands->column("nOtherTracks",nOtherTracks);
    cands->column("trkWeight",trkWeight);
    cands->column("candKSdxy",candKSdxy);
    cands->column("candKSTheta",candKSTheta);
    cands->column("candKSPT",candKSPT);
    cands->column("nCandK",nCandK);
    cands->column("nCandPi",nCandPi);
    cands->column("candPKaons",candPKaons,"nCandK");
    cands->column("candPPions",candPPions,"nCandPi");
    cands->column("candPPi0",candPPi0);
    cands->column("kaonWeight",kaonWeight);
    cands->column("pionWeight",pionWeight);
    cands->column("elecWeight",elecWeight);
    cands->column("muonWeight",muonWeight);
    cands->dumpData();
  }

  delete iterElec; delete iterMuon; delete iterKaon; delete iterTracks;
  delete iterPhotons; delete iterGTL; delete iterGTVL; delete iterPi0; delete iterSoftPi0;
  delete iterKs; delete iterMCTruth; delete iterPion;

  return AppResult::OK;
}



bool DonutAnalysis::isElectron(BtaCandidate* aTrk)
{
  bool ans = false;
  BtaCandidate* aVeto;
  while (aVeto = (*iterElec)()) {
    if (aVeto -> overlaps(*aTrk)) ans = true;
  }
  iterElec->rewind();
  return ans;
}

bool DonutAnalysis::isMuon(BtaCandidate* aTrk)
{
  bool ans = false;
  BtaCandidate* aVeto;
  while (aVeto = (*iterMuon)()) {
    if (aVeto -> overlaps(*aTrk)) ans = true;
  }
  iterMuon->rewind();
  return ans;
}

bool DonutAnalysis::isKaon(BtaCandidate* aTrk)
{
  bool ans = false;
  BtaCandidate* aVeto;
  while (aVeto = (*iterKaon)()) {
    if (aVeto -> overlaps(*aTrk)) ans = true;
  }
  iterKaon->rewind();
  return ans;
}

bool DonutAnalysis::isPion(BtaCandidate* aTrk)
{
  bool ans = false;
  BtaCandidate* aVeto;
  while (aVeto = (*iterPion)()) {
    if (aVeto -> overlaps(*aTrk)) ans = true;
  }
  iterPion->rewind();
  return ans;
}

bool DonutAnalysis::isGTL(BtaCandidate* aTrk)
{
  bool ans = false;
  BtaCandidate* aVeto;
  while (aVeto = (*iterGTL)()) {
    if (aVeto -> overlaps(*aTrk)) ans = true;
  }
  iterGTL->rewind();
  return ans;
}

bool DonutAnalysis::isGTVL(BtaCandidate* aTrk)
{
  bool ans = false;
  BtaCandidate* aVeto;
  while (aVeto = (*iterGTVL)()) {
    if (aVeto -> overlaps(*aTrk)) ans = true;
  }
  iterGTVL->rewind();
  return ans;
}

int DonutAnalysis::DType(BtaCandidate* theD)
{
  HepAListIterator<BtaCandidate> iterDauD = theD -> daughterIterator();
  BtaCandidate* firstD = iterDauD();
  BtaCandidate* secondD = iterDauD();
  iterDauD.rewind();

  if (theD -> nDaughters() == 2) {
    if (abs((int)secondD -> pdtEntry() -> lundId()) == (int)PdtLund::K_plus) return 7; // Ks K
    else return 3; // Ks pi
  }
  if ((int)firstD -> pdtEntry() -> lundId() == (int)PdtLund::K_S0) {
    if (theD -> nDaughters() == 3) return 4; // Ks pi pi0
    if (theD -> nDaughters() == 4) return 5; // Ks 3pi
  }
  else {
    if (abs((int)secondD -> pdtEntry() -> lundId()) == (int)PdtLund::K_plus) return 6; // K K pi
    if (theD -> nDaughters() == 3) return 1; // K pi pi
    if (theD -> nDaughters() == 4) return 2; // K pi pi pi0
  }
  return -1;
}

int DonutAnalysis::D0Type(BtaCandidate* theD)
{
  HepAListIterator<BtaCandidate> iterDauD = theD -> daughterIterator();
  BtaCandidate* firstD = iterDauD();
  BtaCandidate* secondD = iterDauD();
  iterDauD.rewind();

  if (theD -> nDaughters() == 2) {
    if ((int)firstD -> pdtEntry() -> lundId() == (int)PdtLund::K_S0) {
      if (abs((int)secondD -> pdtEntry() -> lundId()) == (int)PdtLund::K_S0) return 6; // Ks Ks
      else return 7; // Ks pi0 
    }  else if (abs((int)firstD -> pdtEntry() -> lundId()) == (int)PdtLund::K_plus) {
      if (abs((int)secondD -> pdtEntry() -> lundId()) == (int)PdtLund::K_plus) return 9; // K K
      else return 1; // K pi
    } else return 8; // pi pi
  }
  if ((int)firstD -> pdtEntry() -> lundId() == (int)PdtLund::K_S0) {
    if (theD -> nDaughters() == 3) return 4; // Ks pi pi
    if (theD -> nDaughters() == 4) return 5; // Ks pi pi pi0
  }
  else {
    if (theD -> nDaughters() == 3) return 2; // K pi pi0
    if (theD -> nDaughters() == 4) return 3; // K 3pi
  }
  return -1;
}
// The Kspipipi0 mode is rejected and the D mass has to be within 4sigma of its value
bool DonutAnalysis::checkD0(float mass, int mode)
{
  if (mode == 6) return false;
  if (mode == 8) return false;
  if (_doSideband.value()) return !sigcheckD0(mass,mode);
  else return sigcheckD0(mass,mode);
}

bool DonutAnalysis::checkD(float mass, int mode)
{
//   if (mode == 4) return false; // It used to be 2, not sure why
//   if (mode == 5) return false;
  if (mode == 6) return false;
  if (_doSideband.value()) return !sigcheckD(mass,mode);
  else return sigcheckD(mass,mode);
}

bool DonutAnalysis::checkDs0(float dsmass,float deltam,int dmode,int dsmode,int kill)
{
  if (dmode == 6) return false;
  if (dmode == 8) return false;
  if (_doSideband.value()) return !sigcheckDs0(dsmass,deltam,dmode,dsmode,kill);
  else return sigcheckDs0(dsmass,deltam,dmode,dsmode,kill);
}

bool DonutAnalysis::checkDs(float dsmass,float deltam,int dmode,int dsmode,int kill)
{
  if (dsmode == 1) {
    if (dmode == 6) return false;
    if (dmode == 8) return false;
  }
  if (dsmode == 2) {
//     if (dmode == 4) return false;
//     if (dmode == 5) return false;
    if (dmode == 6) return false;
  }
  if (_doSideband.value()) return !sigcheckDs(dsmass,deltam,dmode,dsmode,kill);
  else return sigcheckDs(dsmass,deltam,dmode,dsmode,kill);
}

bool DonutAnalysis::sigcheckD0(float mass, int mode)
{
  if (mass < d0cutslow[mode-1] || mass > d0cutshigh[mode-1]) return false;
  return true;
}

bool DonutAnalysis::sigcheckD(float mass, int mode)
{
  if (mass < dccutslow[mode-1] || mass > dccutshigh[mode-1]) return false;
  return true;
}

bool DonutAnalysis::sigcheckDs0(float dsmass,float deltam,int dmode,int dsmode,int kill)
{
  float mass = dsmass - deltam;

  if (mass < ds0d0cutslow[dmode-1] || mass > ds0d0cutshigh[dmode-1]) return false;
  if (deltam < ds0cutslow[dsmode-1] || deltam > ds0cutshigh[dsmode-1]) return false;

  if (kill) {
    long seed = _eventID->eventIdTriplet().timeStamp().binary().lower;
    if (seed < 0) seed *= -1l;
    RandEngine engine(seed);
    RandFlat flat(engine);
    flat.getTheEngine()->setSeed(seed,0);
    float pcut = 0;
    if (dsmode==1) pcut = 372./461.413;
    if (dsmode==2) pcut = 335./389.783;
    if (flat.shoot() < pcut) return true;
    return false;
  }

  return true;
}

bool DonutAnalysis::sigcheckDs(float dsmass,float deltam,int dmode,int dsmode,int kill)
{
  float mass = dsmass - deltam;

  if (dsmode == 1) {
    if (mass < dscd0cutslow[dmode-1] || mass > dscd0cutshigh[dmode-1]) return false;
    if (deltam < dsccutslow[0] || deltam > dsccutshigh[0]) return false;

    if (kill) {
      long seed = _eventID->eventIdTriplet().timeStamp().binary().lower;
      if (seed < 0) seed *= -1l;
      RandEngine engine(seed);
      RandFlat flat(engine);
      flat.getTheEngine()->setSeed(seed,0);
      float pcut = 290./405.959;
      if (flat.shoot() < pcut) return true;
      return false;
    }
    return true;
  } //d0pi

  else {
    if (mass < dscdccutslow[dmode-1] || mass > dscdccutshigh[dmode-1]) return false;
    if (deltam < dsccutslow[1] || deltam > dsccutshigh[1]) return false;

    if (kill) {
      long seed = _eventID->eventIdTriplet().timeStamp().binary().lower;
      if (seed < 0) seed *= -1l;
      RandEngine engine(seed);
      RandFlat flat(engine);
      flat.getTheEngine()->setSeed(seed,0);
      float pcut = 76./103.176;
      if (flat.shoot() < pcut) return true;
      return false;
    }
    return true;
  } //dpi0
}

void DonutAnalysis::countMult(BtaCandidate* aCand,int &nchg, int &nneut)
{
  if (aCand -> nDaughters() == 0) {
    if (aCand -> charge() == 0) {
      nneut ++;
      if (aCand -> p4().theta() > emcbwd) emcbwd = aCand->p4().theta();
      if (aCand -> p4().theta() < emcfwd) emcfwd = aCand->p4().theta();
    }
    else {
      nchg ++;
      if (aCand -> p4().theta() > trkbwd) trkbwd = aCand->p4().theta();
      if (aCand -> p4().theta() < trkfwd) trkfwd = aCand->p4().theta();
    }
    return;
  } else {
    HepAListIterator<BtaCandidate> iterDau = aCand->daughterIterator();
    while (BtaCandidate *adau = iterDau()) {
      countMult(adau,nchg,nneut);
    }
    iterDau.rewind();
  }
}

int DonutAnalysis::expo(int x, int n){
  long base = 1;
  for(long i=0; i<n; i++)
    base *= x;
  return base;
}

int DonutAnalysis::combMode(BtaCandidate *aB, int &bmode, int &dsmode, int &MCComblongD)
{
  int dmode = 0, dmode2 = 0, hmode = 0, unknown = 0, lid;
  int npi = 0, nrho = 0, nk = 0, nkstar = 0, npi0 = 0, nrho0 = 0, nk0 = 0, nkstar0 = 0, ngamma = 0;
  int nel=0, nmu=0, ntau=0, nnu=0, nd0=0, ndp=0, nds0=0, ndsp=0, ndss0=0, ndssp=0, nds=0, ndss=0, ndsss=0;
  int na1=0, na10=0, nomega=0;
  BtaCandidate *aDau,*theDs(0),*theDsstar(0),*theDsstarstar(0);
  int sum = 0;
  MCComblongD = 0;

  HepAListIterator<BtaCandidate> iterDau = aB -> daughterIterator();
  while (aDau = iterDau()) {
    lid = (int)aDau -> pdtEntry() -> lundId();
    if (lid < 0) lid *= -1;

    if (lid == (int)PdtLund::D0){
      if (dmode) dmode2 = 1;
      else dmode = 1;
      nd0++;
    }else if (lid == (int)PdtLund::D_plus){
      if (dmode) dmode2 = 2;
      else dmode = 2;
      ndp++;
    }else if (lid == (int)PdtLund::D_star0){
      if (dmode) dmode2 = 3;
      else dmode = 3;
      nds0++;
    }else if (lid == (int)PdtLund::D_star_plus){
      if (dmode) dmode2 = 4;
      else dmode = 4;
      ndsp++;
    }else if (lid == (int)PdtLund::D_0_star0 ||
	      lid == (int)PdtLund::D_10 ||
	      lid == (int)PdtLund::D_2_star0 ||
	      lid == (int)PdtLund::D_prime_10 ||
	      lid == (int)PdtLund::D_2S0 ||
	      lid == (int)PdtLund::D_star_2S0) {
      if (dmode) dmode2 = 5;
      else dmode = 5;
      ndss0++;
    }else if (lid == (int)PdtLund::D_0_star_plus ||
	      lid == (int)PdtLund::D_1_plus ||
	      lid == (int)PdtLund::D_2_star_plus ||
	      lid == (int)PdtLund::D_prime_1_plus ||
	      lid == (int)PdtLund::D_2S_plus ||
	      lid == (int)PdtLund::D_star_2S_plus){
      if (dmode) dmode2 = 6;
      else dmode = 6;
      ndssp++;
    }else if (lid == (int)PdtLund::D_s_plus) {
      theDs = aDau;
      if (dmode) dmode2 = 7;
      else dmode = 7;
      nds++;
    } else if (lid == (int)PdtLund::D_s_star_plus) {
      theDsstar = aDau;
      if (dmode) dmode2 = 8;
      else dmode = 8;
      ndss++;
    } else if (lid == (int)PdtLund::D_s0_star_plus ||
	lid == (int)PdtLund::D_s1_plus ||
	lid == (int)PdtLund::D_s2_star_plus) {
      theDsstarstar = aDau;
      if (dmode) dmode2 = 9;
      else dmode = 9;
      ndsss++;
    } else if (lid == 22) ngamma++;
    else if (lid == (int)PdtLund::pi_plus) {
      npi ++;
      hmode = 1;
    } else if (lid == (int)PdtLund::pi0) {
      npi0 ++;
    } else if (lid == (int)PdtLund::rho_plus) {
      nrho ++;
      hmode = 2;
    } else if (lid == (int)PdtLund::rho0) {
      nrho0 ++;
    } else if (lid == (int)PdtLund::K_plus) {
      nk ++;
      hmode = 3;
    } else if (lid == (int)PdtLund::K0 || lid == (int)PdtLund::K_S0 || lid == (int)PdtLund::K_L0) {
      nk0 ++;
      hmode = 4;
    } else if (lid == (int)PdtLund::K_star_plus) {
      nkstar ++;
      hmode = 5;
    } else if (lid == (int)PdtLund::K_star0) {
      nkstar0 ++;
    } else if (lid == (int)PdtLund::a_1_plus) {
      na1 ++;
      hmode = 6;
    } else if (lid == (int)PdtLund::a_10) {
      na10 ++;
    } else if (lid == (int)PdtLund::omega) {
      nomega ++;
      hmode = 7;
    } else if (lid == (int)PdtLund::e_minus) {
      nel ++;
    } else if (lid == (int)PdtLund::mu_minus) {
      nmu ++;
    } else if (lid == (int)PdtLund::tau_minus) {
      ntau ++;
    } else if (lid==(int)PdtLund::nu_e||lid==(int)PdtLund::nu_mu||lid==(int)PdtLund::nu_tau) {
      nnu ++;
    } else unknown = 1;
  }
  iterDau.rewind();

  if (unknown) return -100;

  if (npi+nk+nrho+nkstar+npi0+nk0+nrho0+nkstar0+na1+na10+nomega > 1) {
    hmode = 9;
  }
  if((ngamma+nnu)==0){
    if (dmode2 > dmode) bmode = 100*hmode + 10*dmode2 + dmode;
    else bmode = 100*hmode+10*dmode+dmode2;
  }
  if (theDsstarstar) {
    HepAListIterator<BtaCandidate> id = theDsstarstar->daughterIterator();
    while (aDau = id()) {
      lid = aDau->pdtEntry()->lundId();
      if (lid < 0) lid *= -1;
      if (lid == (int)PdtLund::D_s_star_plus) theDsstar = aDau;
      if (lid == (int)PdtLund::D_s_plus) theDs = aDau;
    }
    id.rewind();
  }
  if (theDsstar) {
    HepAListIterator<BtaCandidate> id = theDsstar->daughterIterator();
    while (aDau = id()) {
      lid = aDau->pdtEntry()->lundId();
      if (lid < 0) lid *= -1;
      if (lid == (int)PdtLund::D_s_plus) theDs = aDau;
    }
    id.rewind();
  }

  int foundLep=0,foundMeson=0;
  if (theDs) {
    HepAListIterator<BtaCandidate> id = theDs->daughterIterator();
    while (aDau = id()) {
      lid = aDau->pdtEntry()->lundId();
      if (lid < 0) lid *= -1;
      if (lid == (int)PdtLund::e_minus) foundLep = 1;
      if (lid == (int)PdtLund::mu_minus) foundLep = 2;
      if (lid == (int)PdtLund::tau_minus) foundLep = 3;
      if (lid == (int)PdtLund::phi) foundMeson = 1;
      if (lid == (int)PdtLund::eta) foundMeson = 2;
      if (lid == (int)PdtLund::eta_prime) foundMeson = 3;
    }
    id.rewind();
  }

  if (foundLep==3 && foundMeson==0) dsmode = 1;
  if (foundLep==2 && foundMeson==0) dsmode = 2;
  if (foundLep==1 && foundMeson==0) dsmode = 3;
  if (foundLep>0&&foundLep<3&&foundMeson==1) dsmode = 4;
  if (foundLep>0&&foundLep<3&&foundMeson==2) dsmode = 5;
  if (foundLep>0&&foundLep<3&&foundMeson==3) dsmode = 6;
  if(npi>3||npi0>3||nk>3||nk0>3||nrho>3||nrho0>3||nkstar>3||nkstar0>3||nel>3||nmu>3||ntau>3||
     nd0>3||ndp>3||nds0>3||ndsp>3||ndss0>3||ndssp>3||nds>3||ndss>3||ndsss>3||ngamma>3||
     na1>3||na10>3||nomega>3) return -10;

  sum = nel+nmu*expo(4,1)+ntau*expo(4,2)+npi*expo(4,3)+npi0*expo(4,4)+nk*expo(4,5)+nk0*expo(4,6)
    +nrho*expo(4,7)+nrho0*expo(4,8)+nkstar*expo(4,9)+nkstar0*expo(4,10)+na1*expo(4,11)
    +na10*expo(4,12)+nomega*expo(4,13)+ngamma*expo(4,14);

  MCComblongD = nd0+ndp*expo(4,1)+nds0*expo(4,2)+ndsp*expo(4,3)+ndss0*expo(4,4)+ndssp*expo(4,5)
    +nds*expo(4,6)+ndss*expo(4,7)+ndsss*expo(4,8);

  return sum;
}

int DonutAnalysis::motherB(BtaCandidate *X){
  static const PdtLund::LundType B0Lund = PdtLund::B0;
  static const PdtLund::LundType BcLund = PdtLund::B_plus;

  PdtLund::LundType XLund = X->pdtEntry()->lundId();
  if(abs(XLund) == B0Lund || abs(XLund) == BcLund)
    return X->uid();
  else {
    BtaCandidate *Mom = X->theMother();
    if(Mom == 0) return 0;
    return motherB(Mom);
  }
}

// Returns the number of charged and neutral final particles (e,mu,pi,pi0,K,gamma) non-truthmatched
// or coming from the wrong B
void DonutAnalysis::nTruthMatched(BtaCandidate *X, int &ChargedNon, int &NeutralNon, int &nCharged,
				  int &nNeutral, BtaMcAssoc* truthMap, int momB){
  static const PdtLund::LundType eLund = PdtLund::e_minus;
  static const PdtLund::LundType muLund = PdtLund::mu_minus;
  static const PdtLund::LundType piLund = PdtLund::pi_plus;
  static const PdtLund::LundType pi0Lund = PdtLund::pi0;
  static const PdtLund::LundType KLund = PdtLund::K_plus;
  static const PdtLund::LundType gLund = PdtLund::gamma;

  if(X==0) return;

  PdtLund::LundType XLund = (PdtLund::LundType)abs(X->pdtEntry()->lundId());

  // Have we reached a stable particle yet?  
  if(XLund==eLund||XLund==muLund||XLund==piLund||XLund==pi0Lund||XLund==KLund||XLund==gLund){
    if(abs(X->charge())==0) nNeutral++;
    else nCharged++;
    BtaCandidate *trueX = truthMap->mcFromReco(X);
    if(trueX == 0 && XLund==eLund){
      HepAListIterator<BtaCandidate> dauglepIter = X->daughterIterator();
      while (BtaCandidate* son = dauglepIter()) {
	int pdt = son->pdtEntry()->lundId();
	if(abs(pdt)==(int)PdtLund::e_minus) trueX = truthMap->mcFromReco(son);
      }
    }
    if(trueX == 0 || momB != motherB(trueX)){
      if(abs(X->charge())==0) NeutralNon++;
      else ChargedNon++;
    }
  } else {		// Evaluate daughters of this intermediate particle
    HepAListIterator<BtaCandidate> dauXIter = X->daughterIterator();
    BtaCandidate* dauX = 0;
    while (0 != (dauX = dauXIter())) {
      nTruthMatched(dauX, ChargedNon, NeutralNon, nCharged, nNeutral, truthMap, momB);
    }
  }
}

