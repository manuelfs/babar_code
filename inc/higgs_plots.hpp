//------------------------------------------------------------------------
// File and Version Information:
//      $Id: HiggsPlot.hh,v 1.1 2012/03/04 00:33:02 manuelf Exp $
//
// Description:
//    Calculation of the B->TauNu branching fraction. 
//
//    PlotBF(double tBmH_max)
//       Plots the BF for different values of tanBeta/mH
//    SetMasses(isBm) 
//       Sets the masses of the B and D* depending on the charge.
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//      Michael Mazur                             INFN Pisa
//
// Revision History:
//      12/04/05 manuelf -- Created off maketaunuplots.hh by M. Mazur
//------------------------------------------------------------------------

#ifndef HIGGSPLOT
#define HIGGSPLOT

#include "TMath.h"
#include "TString.h"
#include "TSpline.h"
#include <fstream>
#include <iostream>
#define nVar 4

using namespace std;
using std::cout;
using std::endl;

class HiggsPlot {
public:
  HiggsPlot(TString NameFile="babar_code/NewPhysics/Measurements_BaBar");  
  ~HiggsPlot();

  // Primary constants from PDG and PRD 82 112007 (2010)
  static constexpr double mTau = 1.7768;
  static constexpr double mMu  = 0.10566;
  static constexpr double mE   = 0.000511;
  static constexpr double GF   = 0.000011664;
  static constexpr double PI   = 3.14159265;
  static constexpr double hbar = 6.582119e-25; // [GeV*s]

  // Secondary constants
  double Vub[2], fB[2], Measurement[nVar][2];
  double mB, BLifeTime;
  int _isBm, _isgS, _varyRD;

  TString DecayName[4];
  double RDCoef[2][3][3];
  TSpline3 *MeasuredRD[4];

  void SetMasses(int isBm);
  void PlotgSLPRL(int isgSR = 0, double tBmH_max=2.1);
  void PlotgSLTauNu(double b=0.12, double tBmH_max=2.5);
  void PlotPRL(int isPsfrag = 1);
  void PlotBF(int iDecay, double tBmH_max=1., double BF_max=-1);
  void PlotExclusion(int iDecay, double tBmH_max=1.2, TString Option="col");
  void PlotSExclusion(double rR_max=3.5, TString Option="col");
  void PlotChi2();
  double ProbChi2(int iDecay, double tBmH, double rL=0);
  double ProbChi2(int iDecay, double tBmH, double RD[2], double RDs[2]);
  double Chi2(double M1[2], double M2[2]);
  void BFTauNu(double tBmH, double BF[2], double ml=mTau);
  void RDCalc(double tBmH, double RD[2], int isDs);
  void Compute(double tBmH, double RD[2], int iDecay);
  double IntG(double mean, double sigma, double minX, double maxX);
  TString RoundNumber(double n, int e, double d=1);
  void ReadRD(TString NameFile);


};

#endif

