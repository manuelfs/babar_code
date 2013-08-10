//----------------------------------------------------------------------
// $Id: DonutUserAUB.cc,v 1.1 2009/04/27 17:46:20 manuelf Exp $
//
//----------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "Framework/AppUserBuild.hh"

// Standard analysis sequences everyone needs
#include "BetaMiniEnvSequences/BetaMiniEnvSequence.hh"
#include "BetaMiniSequences/BetaMiniSequence.hh"
#include "BetaMiniSequences/BetaMiniActions.hh"
#include "BetaMiniSequences/BetaMiniPhysicsSequence.hh"

// Local analysis modules
#include "BetaMicro/BtaSelectCand.hh"
#include "DonutUser/DonutMCTagger.hh"
#include "DonutUser/DonutEVeto.hh"
#include "DonutUser/DonutAnalysis.hh"


AppUserBuild::AppUserBuild( AppFramework* theFramework )
  : AppBuild(theFramework) {
  BetaMiniEnvSequence(this);		// Standard analysis sequences
  BetaMiniSequence(this);
  BetaMiniActions(theFramework);
  BetaMiniPhysicsSequence(this);

  add(new BtaSelectCand("GoodPhotonsVisibleESelection" , 
			       "Selection of good photons from BAD #633",
			       "BtaGoodNeutralSelector" ) );
  add(new BtaSelectCand("GoodTracksVisibleESelection", 
			       "Selection of good tracks from BAD #633",
			       "BtaGoodTrkSelector" ) );

  add(new DonutMCTagger("DonutMCTagger", "DonutMCTagger"));
  add(new DonutEVeto   ("DonutEVeto"   , "DonutEVeto"));
  add(new DonutAnalysis("DonutAnalysis", "DonutAnalysis"));
}

//--------------
// Destructor --
//--------------

AppUserBuild::~AppUserBuild( ) {}

