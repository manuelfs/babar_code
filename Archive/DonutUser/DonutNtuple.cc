// $Id: DonutNtuple.cc,v 1.1 2009/04/29 05:26:20 kelsey Exp $

#include "BaBar/BaBar.hh"
#include "DonutUser/DonutNtuple.hh"

#include "DonutUtils/weightManager.hh"
#include "AbsEvent/AbsEvent.hh"
#include "AbsEventTag/AbsEventTag.hh"
#include "AbsEnv/AbsEnv.hh"
#include "GenEnv/GenEnv.hh"
#include "HepTuple/Tuple.h"
#include "HepTuple/TupleManager.h"
#include "ProxyDict/Ifd.hh"


DonutNtuple::DonutNtuple(const char* const theName,
                       const char* const theDescription ) :
  AppModule( theName, theDescription ) {}

DonutNtuple::~DonutNtuple( ) {}

AppResult DonutNtuple::beginJob( AbsEvent* anEvent )
{
  HepTupleManager* manager = gblEnv->getGen()->ntupleManager();
  assert(manager != 0);
  _tuple = manager->ntuple("DonutNtuple Tuple");
  myWM = new WeightManager("dataFits/wFF");
  return AppResult::OK;
}

AppResult DonutNtuple::endJob( AbsEvent* anEvent )
{
  delete myWM;
  return AppResult::OK;
}

AppResult DonutNtuple::event( AbsEvent* anEvent )
{
  AbsEventTag* tag = Ifd<AbsEventTag>::get(anEvent);
  int MCType,MCSubmode,MCPions,MCBrem,MCScatter,MCD,isBzero,
    MCDssmode,MCUnsim,trueLepCharge;
  float truePLep,trueTauE,trueCTL,trueCTV,trueChi,trueQ2,trueLepE,trueCTL2,trueMM2;
  float weight;

  tag -> getInt(MCType     ,"MCType");
  tag -> getInt(MCSubmode  ,"MCSubmode");
  tag -> getFloat(trueCTL  ,"trueCTL");
  tag -> getFloat(trueCTV  ,"trueCTV");
  tag -> getFloat(trueChi  ,"trueChi");
  tag -> getFloat(trueQ2   ,"trueQ2");
  tag -> getInt(trueLepCharge,"trueLepCharge");
  tag -> getInt(MCDssmode  ,"MCDssmode");
  tag -> getInt(MCD        ,"MCD");
  tag -> getInt(MCPions    ,"MCPions");
  tag -> getInt(isBzero    ,"isBzero");
  tag -> getInt(MCUnsim    ,"MCUnsim");
  tag -> getFloat(truePLep ,"truePLep");
  tag -> getFloat(trueTauE ,"trueTauE");
  tag -> getFloat(trueLepE ,"trueLepE");
  tag -> getFloat(trueCTL2 ,"trueCTL2");
  tag -> getFloat(trueMM2  ,"trueMM2");
  tag -> getInt(MCBrem     ,"MCBrem");
  tag -> getInt(MCScatter  ,"MCScatter");

  weight = myWM->getEventWeight(0,0,MCType,MCSubmode,MCDssmode,MCD,MCPions,isBzero,0,
				MCUnsim,0.0,0.0,trueCTL,trueCTV,trueChi,trueQ2,trueLepCharge,0.0);

  _tuple->column("MCType",MCType);
  _tuple->column("MCSubmode",MCSubmode);
  _tuple->column("trueCTL",trueCTL);
  _tuple->column("trueCTV",trueCTV);
  _tuple->column("trueChi",trueChi);
  _tuple->column("trueQ2",trueQ2);
  _tuple->column("trueTauE",trueTauE);
  _tuple->column("trueLepE",trueLepE);
  _tuple->column("trueCTL2",trueCTL2);
  _tuple->column("trueMM2",trueMM2);
  _tuple->column("trueLepCharge",trueLepCharge);
  _tuple->column("MCDssmode",MCDssmode);
  _tuple->column("MCD",MCD);
  _tuple->column("MCPions",MCPions);
  _tuple->column("isBzero",isBzero);
  _tuple->column("MCUnsim",MCUnsim);
  _tuple->column("truePLep",truePLep);
  _tuple->column("MCBrem",MCBrem);
  _tuple->column("MCScatter",MCScatter);
  _tuple->column("weight",weight);
  _tuple->dumpData();

  return AppResult::OK;
}
