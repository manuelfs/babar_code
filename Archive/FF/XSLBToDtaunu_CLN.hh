//*****************************************************************
//   
//   Creation: Michael Mazur, INFN Pisa, July 2007
//   Based on code by: Martin Simard & David Cote, Universite de Montreal, September 2006 
//   From: EvtGenModels/EvtHQET2FF.cc
//   Description: form factors for B->Dlnu according to HQET with dispersive FF
//   Use dispersion relation parametrization from 
//   I.Caprini, L.Lellouch, M.Neubert, Nucl. Phys. B 530,153(1998)
//   Using decay distribution including lepton mass effects from Korner&Schuler, Z Phys C 46, 93 (1990).
//
//   Please also see general comments in the XSLEvtFFWeight.hh file or BAD#809

#ifndef XSLBTODTAUNU_CLN
#define XSLBTODTAUNU_CLN

#include "XslFFReweighting/XSLPseudoScalarFF.hh"
#include "CLHEP/Vector/LorentzVector.h"

class XSLKin;

class XSLBToDtaunu_CLN : public XSLPseudoScalarFF {
public:
  XSLBToDtaunu_CLN(double mB, double mD, double q2, double theta_l,double rho2=1.17,double v1_1=1.0);
  XSLBToDtaunu_CLN(HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector DLab,double rho2=1.17,double v1_1=1.0);
  XSLBToDtaunu_CLN( XSLKin* DecayKin, double rho2=1.17,double v1_1=1.0);
  virtual ~XSLBToDtaunu_CLN();
  
  virtual double FromISGW2ToThisModel();
  virtual double FromPHSPToThisModel();
  virtual double FromFLATQ2ToThisModel();
  virtual double FromSP4ToThisModel(){ return FromISGW2ToThisModel(); }
  virtual double FromSP5ToThisModel(){ return FromSP4ToThisModel(); }
  virtual double FromSP6ToThisModel(){ return FromSP4ToThisModel(); }
  virtual double FromSP7ToThisModel(){ return FromSP4ToThisModel(); }
  virtual double FromSP8ToThisModel(){ return FromSP4ToThisModel(); }

protected:
  virtual void Compute();
  void ComputeISGW2();
  void ComputeCLN();
  double EvtGetas (double massq, double massx);
  double EvtGetas (double mass);
  double EvtGetGammaji(double z);
  virtual double GetFplus(double q2);
  virtual void SetNormalizations(std::string mode);

  double _rho2;
  double _v1_1;

  double kin_p,kin_mu2;
  double isgw2_fplus,isgw2_fminus;
  double cln_fplus,cln_fminus;
};

#endif



