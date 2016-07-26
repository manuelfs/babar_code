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
// Module: EvtGen/Evt2HDMFF.hh
//
// Description:Form factor routines specific to Evt2HDM
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVT2HDMFF_HH
#define EVT2HDMFF_HH

#include "EvtGenBase/EvtSemiLeptonicFF.hh"

class EvtId;

class Evt2HDMFF : public EvtSemiLeptonicFF {

public:
  Evt2HDMFF(double tBmH=0.);

  void getscalarff( EvtId parent, EvtId daught,
                       double t, double mass, double *fpf,
                       double *f0f );
  void getvectorff( EvtId parent, EvtId daught,
                       double t, double mass, double *a1f,
                       double *a2f, double *vf, double *a0f );
  void gettensorff( EvtId parent, EvtId daught,
                       double t, double mass, double *hf,
                       double *kf, double *bpf, double *bmf );


private:

    // getscalarff, getvectorff, and gettensorff call the
    // correct thdm form factor routine which computes 
    // form factors according to the 2HDM paper.


  void Evt2HDMFF3S1( EvtId parent, EvtId daught, 
                       double t, double mass, double *ff, double *gf,
                       double *apf, double *amf);
  void Evt2HDMFF23S1( EvtId parent, EvtId daught,
                        double t, double mass,double *fpf, double *gpf,
                        double *app, double *apm);
  void Evt2HDMFF3P1( EvtId parent, EvtId daught, 
                       double t, double mass,double *lf, double *qf,
                       double *cpf, double *cmf);
  void Evt2HDMFF3P0( EvtId parent, EvtId daught, 
                       double t, double mass, double *upf, double *umf);
  void Evt2HDMFF1S0( EvtId parent, EvtId daught, 
		       double t, double mass,double *fpf, double *fmf);
  void Evt2HDMFF21S0( EvtId parent, EvtId daught,
                        double t, double mass, double *fppf, double *fpmf);
  void Evt2HDMFF3P2( EvtId parent, EvtId daught, 
                       double t, double mass, double *h, double *k,
                       double *bp, double *bm);
  void Evt2HDMFF1P1( EvtId parent, EvtId daught, 
                       double t, double mass, double *rf, double *vf,
                       double *spf, double *smf);

  double EvtGetas( double mass );
  double EvtGetas( double mass,double mass1  );
  double EvtGetGammaji( double z );

  double _tBmH; // tanBeta/mH in GeV^-1

};

#endif


