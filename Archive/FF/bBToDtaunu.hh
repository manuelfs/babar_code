//------------------------------------------------------------------------
// File and Version Information:
//      $Id: RateCalc.hh,v 1.1 2012/03/04 00:33:02 manuelf Exp $
//
// Description:
//    Calculation of the B->DlNu rates for a lepton of mass ml
//    based on CLN FF from Tanaka, Watanabe (2010).
//    The angular dependence is from Korner, Shuler (1990).
//    Normalization is calculated for each set of FF values.
//    The New Physics dependence is parameterized in terms of
//    gSR = -mb (tanBeta/mH)^2 from hep-ph/1203.2654
//
//    FromSP8ToThisModel(q2, thetaL, ml) 
//       Re-weights SP8 MC to the CLN parameterization.
//    SetMasses(isBm) 
//       Sets the masses of the B and D depending on the charge.
//    Rate(isCLN, ml) 
//       Branching fraction.
//    Gamma_q2(q2, fplus, fminus, ml) 
//       q2 spectrum.
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//      Michael Mazur                             INFN Pisa
//
// Revision History:
//      12/05/10 manuelf -- Normalization validated with EvtGen, including Higgs
//                          Added the theta spectrum integrated over q2
//      12/05/09 manuelf -- Corrected a missing sin(thetaL)
//      12/05/03 manuelf -- Added the NP dependence in terms of gSR
//      12/03/30 manuelf -- Created off XSLBToDtaunu_CLN.cc by M. Mazur
//------------------------------------------------------------------------


#ifndef BTODTAUNU
#define BTODTAUNU

#include "TMath.h"
#include "TString.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

class BToDtaunu { 
public:
  BToDtaunu(double rho2=1.186, double Delta=1., double gSR=0);  
  ~BToDtaunu();

  // Primary constants from PDG and HFAG 2010
  static const double mTau = 1.7768;
  static const double mMu  = 0.10566;
  static const double mE   = 0.000511;
  static const double Vcb  = 0.03912;
  static const double V11  = 1.0816;
  static const double GF   = 0.000011664;
  static const double PI   = 3.14159265;
  static const double hbar = 6.582119e-25;
  static const double mb_quark = 4.20;  // [GeV] in the MSbar scheme, evaluated at mb. 
  static const double mc_quark = 0.901; // From Xing, Zhang, Zhou (2008)

  // Secondary constants
  double _mB, _mBISGW2, _BLifeTime;
  double _mD, _mDISGW2, _Dmaxq2; 
  double _rho2, _Delta, _gSR;
  double _isBm;


  void SetMasses(int isBm);
  void ComputeCLN(double q2, double &fplus, double &fminus);
  void ComputeISGW2(double q2, double &fplus, double &fminus);
  double EvtGetas (double massq, double massx);
  double EvtGetGammaji(double z);
  double Compute(double q2, double thetaL, int isCLN, double ml);
  double Normalization(double ml);
  double Gamma_q2tL(double q2, double thetaL, double fplus, double fminus, double ml);
  double FromSP8ToThisModel(double q2, double thetaL, double ml);

  double Compute(double q2, int isCLN, double ml);
  double ComputetL(double thetaL, int isCLN, double ml);


  // Functions only used to calculate the Normalization once
  double IntRate(double minX, double maxX, int isCLN, double thetaL, double ml, int nPoints=10000);
  double Gamma_q2(double q2, double fplus, double fminus, double ml);
  double Rate(int isCLN, double ml);
  double SubPoly(double C[3][3]);
  void Polynomial(double ml);
  double SumPoly(double Coef[3][3]);

  // Auxiliary functions
  void SetTanBetamH(double tBmH);
  void SetFF(double rho2=1.186, double Delta=1., double gSR=0);
  void SetNorm(double ml);
  double _Normalization;

};

#endif

