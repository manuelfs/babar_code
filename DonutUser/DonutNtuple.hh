#ifndef DTNINFO
#define DTNINFO
//$Id: DonutNtuple.hh,v 1.1 2009/04/29 05:26:21 kelsey Exp $

#include "Framework/AppModule.hh"

class AbsEvent;
class HepTuple;
class WeightManager;


class DonutNtuple : public AppModule {
public:
  DonutNtuple( const char* const theName, const char* const theDescription );
  virtual ~DonutNtuple();
  virtual AppResult beginJob( AbsEvent* anEvent );
  virtual AppResult event   ( AbsEvent* anEvent );
  virtual AppResult endJob  ( AbsEvent* anEvent );
  HepTuple * _tuple;
  WeightManager *myWM;
};

#endif
