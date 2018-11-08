//------------------------------------------------------------------------
// File and Version Information:
//      $Id: RateCalc.hh,v 1.1 2012/03/04 00:33:02 manuelf Exp $
//
// Description:
//      RateCalc - Calculation of the B->D(*)TauNu rates based on
//                 hep-ph 1203.2654
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      12/03/16 manuelf -- Created 
//------------------------------------------------------------------------

#ifndef RATECALC
#define RATECALC

#include "TMath.h"
#include "TString.h"
#include "TMatrixT.h"
#include <fstream>
#include <iostream>
#define nVar 11

using namespace std;
using std::cout;
using std::endl;

class RateCalc {
public:
  RateCalc(TString NameFile="babar_code/NewPhysics/FFinputs_HFAG11");  
  ~RateCalc();

  // Primary constants
  static constexpr double mTau = 1.7768;
  static constexpr double mMu  = 0.10566;
  static constexpr double mE   = 0.000511;
  static constexpr double GF   = 0.000011664;
  static constexpr double pi   = 3.14159265;
  static constexpr double hbar = 6.582119e-25;
  double mB, mDs, mD, BLifeTime; 

  // Secondary constants
  double RDs, mB2, mDs2, Dsmaxq2, mD2, Dmaxq2;
  double gsR, gsL;
  double mb, mc, mb_Err, mc_Err; 
  int _isBm;
  enum FormFactors{irhoD2, iDelta, irhoDs2, iR1, iR2,
		   iCorrelation0, iCorrelation1, iCorrelation2, iR0, imb, imc};

  // ============================  R(D*) ========================================
  // Kinematic variables
  double wDs(double q2); double pDs(double q2);

  // D*TauNu FF parameters. 
  double Vcb, F1, rhoDs2, R[4], rhoDs2_Err, R_Err[4], Correlation[3]; 

  // D*TauNu Form Factors
  double hA1(double w); double R0(double w);  double R1(double w);  double R2(double w);
  double V(double q2);  double A0(double q2); double A1(double q2); double A2(double q2);

  // Helicity amplitudes
  double Hpp(double q2); double Hmm(double q2); double H00(double q2); double H0t(double q2); 
  
  // Rate
  double GammaDs_q2(double q2, double ml=mTau);

  // ============================  R(D) =========================================
  // DTauNu FF parameters. 
  double VcbG1, rhoD2, Delta, rhoD2_Err, Delta_Err; 

  // DTauNu Form Factor
  double G(double w);
  
  // Rate
  double GammaD_q2(double q2, double ml=mTau);
  double GammaD_q2_Tanaka(double q2, double ml=mTau);

  // =========================  Auxiliary functions  =============================
  double IntRate(double minX, double maxX, double ml, int isDs = 1, int nPoints=10000);
  double Rate(int isDs, double ml); 
  double RRate(int isDs, double ml=mE); 
  void Print();
  void SetMasses(int isBm=1);
  void Polynomial(double Coef[3], int isDs, double ml=mE);
  void AverPoly(double Coef[3], int isDs);
  TString RoundNumber(double n, int e, double d=1);
  double Gamma_q2(int isDs, double q2, double ml=mTau);

  void ReadFF(TString NameFile);
  void Errors(double Coef[3][3], int isDs, int nRep=10000);
  TMatrixT<double> Choleski(TMatrixT<double> A, int nRows);



};

#endif

