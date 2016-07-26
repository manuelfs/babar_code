//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtGen/Evt2HDM.hh
//
// Description:Implementation of the 2HDM model
// Class to handle semileptonic decays using the 2HDM
// model, as described in PRD 52 5 (1995) by 
// Isgur and Scora.  Electron, muon, and tau models
// are available.  Form factors, q2 and lepton energy
// spectra checked against code from Scora.
//
// Modification history:
//
//    DJL/RYD     September 25, 1996         Module created
//
//------------------------------------------------------------------------

#ifndef EVT2HDM_HH
#define EVT2HDM_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"
#include "EvtGenBase/EvtSemiLeptonicAmp.hh"

class EvtParticle;

class Evt2HDM:public  EvtDecayAmp  {

public:

  Evt2HDM();
  virtual ~Evt2HDM();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p);
  void initProbMax();
  void init();

private:
  EvtSemiLeptonicFF *thdmffmodel;
  EvtSemiLeptonicAmp *calcamp;
};

#endif

