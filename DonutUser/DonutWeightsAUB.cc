//------------------------------------------------------------------------
// File and Version Information:
//      $Id: DonutWeightsAUB.cc,v 1.1 2009/04/29 05:26:21 kelsey Exp $
//
// Description:
//	GeneratorsQA-style executable for validating the D(*)tau nu Monte
//	Carlo generator and form-factor reweighting code.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Michael Mazur <mazur@slac.stanford.edu>
//
// Copyright Information:
//	N.A.
//
// Revision History:
//	20090428  M. Kelsey -- Copied from uncommitted GeneratorsQA code
//		  into DonutUser pacakge.
//------------------------------------------------------------------------
//---------------
// BaBAr Headers:
//---------------
#include "BaBar/BaBar.hh"

//-------------
// C++ Headers:
//-------------
#include <signal.h>

//-----------------------
// This Class's Header --
//-----------------------
#include "Framework/AppFramework.hh"
#include "Framework/AppUserBuild.hh"

#include "CdbEvtModules/BdbSetTime.hh"
#include "CdbFwkModules/CdbApiInitModule.hh"
#include "DonutUser/DonutMCTagger.hh"
#include "DonutUser/DonutNtuple.hh"
#include "ErrLogger/ErrLog.hh"
#include "FrameLogger/AppActionErrLogger.hh"
#include "GenEnv/GenBuildEnv.hh"
#include "GenFwkInt/GfiSequence.hh"
#include "HbkTupleEnv/HbkTupleEnv.hh"
#include "RandControl/RacRandomControl.hh"
#include "RandControl/RacTestInput.hh"
#include "RecoUtils/AppActionHistoDir.hh"
#include "RecoUtils/AppActionSumTime.hh"
#include "RecoUtils/PdtInit.hh"
#include "RooTupleEnv/RooTupleEnv.hh"
#include "TruthTools/BtaLoadMcCandidates.hh"


AppUserBuild::AppUserBuild( AppFramework* theFramework )
    : AppBuild( theFramework ) {
  add(new GenBuildEnv( "GenBuildEnv", "Build General Environment" ));
  add(new HbkTupleEnv( "HbkTupleEnv", "Hbook" ));
  add(new RooTupleEnv( "RooTupleEnv", "Root" ));
  add(new PdtInit("PdtInit", "PdtInit")); // initialize PDT table
  add(new CdbApiInitModule("CdbApiInitModule", "Initialize CDB" ));
  add(new BdbSetTime("BdbSetTime", "Set the time before PepBuildEnv"));

  // Input module provided by RandControl:
  add(new RacTestInput);

  // Random control module:
  add(new RacRandomControl);

  // create generator modules and other needed stuff.
  GfiSequence(this);

  // Load BtaCandidates:
  add (new BtaLoadMcCandidates("BtaLoadMcCandidates",
			       "Load BtaCandidates from StdHepMC"));

  // Analysis modules:  
  add (new DonutMCTagger("DonutMCTagger", "DonutMCTagger"));
  add (new DonutNtuple("DonutNtuple","DonutNtuple"));

  // Framework actions:
  theFramework->actions()->append(new AppActionSumTime); 
  theFramework->actions()->append(new AppActionHistoDir); 
  theFramework->actions()->append(new AppActionErrLogger); 
}

AppUserBuild::~AppUserBuild() {}


