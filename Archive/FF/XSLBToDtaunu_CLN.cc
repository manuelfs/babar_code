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
#include "BaBar/BaBar.hh"

#include "XslFFReweighting/XSLBToDtaunu_CLN.hh"

#include "XslFFReweighting/XSLKin.hh"
#include "CLHEP/Vector/LorentzVector.h"
#include <assert.h>

using std::cout;
using std::endl;
using std::string;

XSLBToDtaunu_CLN::XSLBToDtaunu_CLN(double mB, double mD, double q2, double theta_l, double rho2,double v1_1) : 
  XSLPseudoScalarFF(mB,mD,q2,theta_l) 
{ 
  _rho2=rho2;
  _v1_1=v1_1;
  Compute();  
}


XSLBToDtaunu_CLN::XSLBToDtaunu_CLN(HepLorentzVector BLab, HepLorentzVector LepLab, HepLorentzVector DLab,double rho2,double v1_1) : 
  XSLPseudoScalarFF( BLab, LepLab, DLab) 
{ 
  _rho2=rho2;
  _v1_1=v1_1;
  Compute();  
}

XSLBToDtaunu_CLN::XSLBToDtaunu_CLN( XSLKin* DecayKin,double rho2,double v1_1 ): XSLPseudoScalarFF(DecayKin)
{ 
  _rho2=rho2;
  _v1_1=v1_1;
  Compute();  
}

XSLBToDtaunu_CLN::~XSLBToDtaunu_CLN()
{}

void
XSLBToDtaunu_CLN::Compute()
{
  kin_p = _mB*_mB*_mB*_mB + _mXu2*_mXu2 + _q2*_q2 -
    2.*_mB*_mB*_mXu2 - 2.*_mB*_mB*_q2 - 2.*_mXu2*_q2;
  if (kin_p < 0) kin_p = 0;
  else kin_p = sqrt(kin_p)/(2.*_mB);
  kin_mu2 = 1.777*1.777;

  ComputeISGW2();
  ComputeCLN();
}

void
XSLBToDtaunu_CLN::ComputeISGW2()
{
  //This code liberally stolen from EvtGenModels/EvtISGW2FF.cc
  double msb=5.2;
  double msd=0.33;
  double bb2=0.431*0.431;
  double mbb=5.31;
  double nf = 4.0;
  double msq=1.82;
  double bx2=0.45*0.45;
  double mbx=0.75*2.01+0.25*1.87;
  double nfp = 3.0;
  double mtb = msb + msd;
  double mtx = msq + msd;
  double mb = _mB;
  double mx = _mXu;
  double mup=1.0/(1.0/msq+1.0/msb);
  double mum=1.0/(1.0/msq-1.0/msb);
  double bbx2=0.5*(bb2+bx2);
  double tm=(mb-mx)*(mb-mx);
  double t = _q2;
  if (t > tm) t=0.99*tm;
  double wt=1.0+(tm-t)/(2.0*mbb*mbx);
  double mqm = 0.1;
  double r2=3.0/(4.0*msb*msq)+3*msd*msd/(2*mbb*mbx*bbx2) + 
    (16.0/(mbb*mbx*(33.0-2.0*nfp)))*
    log(EvtGetas(mqm,mqm)/EvtGetas(msq,msq));
  double f3 = sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,1.5) /
    ((1.0+r2*(tm-t)/12.0)*(1.0+r2*(tm-t)/12.0));
  double w = 1.0 + (( tm - t ) / ( 2.0* mb * mx ));
  double rcji = ( 1/sqrt(w*w -1 ))*log( w + sqrt( w*w -1 ));
  double al = (8.0 / ( 33.0 - 2.0*nfp ))*(w*rcji -1.0 );
  double ai = -1.0* ( 6.0/( 33.0 - 2.0*nf));  
  double cji = pow(( EvtGetas( msb,msb ) / EvtGetas( msq,msq ) ),ai);
  double zji = msq / msb;
  double gammaji = EvtGetGammaji( zji );
  double chiji = -1.0 - ( gammaji / ( 1- zji ));
  double betaji_fppfm = gammaji - (2.0/3.0)*chiji;
  double betaji_fpmfm = gammaji + (2.0/3.0)*chiji;
  double rfppfm = cji *(1.0 + betaji_fppfm*EvtGetas( msq,sqrt(msb*msq) )/PI);
  double rfpmfm = cji *(1.0 + betaji_fpmfm*EvtGetas( msq,sqrt(msb*msq) )/PI);
  double f3fppfm = f3*pow(( mbb / mtb ),-0.5)*pow((mbx/mtx),0.5);
  double f3fpmfm = f3*pow(( mbb / mtb ),0.5)*pow((mbx/mtx),-0.5);
  double fppfm = f3fppfm* rfppfm * ( 2.0 - ( ( mtx/msq)*(1- ( (msd*msq*bb2)
							      /(2.0*mup*mtx*bbx2)))));
  double fpmfm = f3fpmfm* rfpmfm * ( mtb/msq) * ( 1 - ( ( msd*msq*bb2)/
							( 2.0*mup*mtx*bbx2)));
  isgw2_fplus = (fppfm + fpmfm)/2.0;
  isgw2_fminus = (fppfm - fpmfm)/2.0;
}

double
XSLBToDtaunu_CLN::EvtGetas(double massq, double massx)
{
  double lqcd2 = 0.04;
  double nflav = 4;
  double temp = 0.6;
  if ( massx > 0.6 ) {
    if ( massq < 1.85 ) {
      nflav = 3.0;}
    temp = 12.0*PI / ( 33.0 - 2.0*nflav) /
      log( massx*massx/lqcd2);
  }
  return temp;
}

double
XSLBToDtaunu_CLN::EvtGetas(double mass)
{
  double lqcd2 = 0.04;
  double nflav = 4;
  double temp = 0.6;
  if ( mass > 0.6 ) {
    if ( mass < 1.85 ) {
      nflav = 3.0;}
    temp = 12.0*PI / ( 33.0 - 2.0*nflav) /
      log( mass*mass/lqcd2);
  }
  return temp;
}

double
XSLBToDtaunu_CLN::EvtGetGammaji(double z)
{
  double temp;
  temp = 2+((2.0*z)/(1-z))*log(z);
  temp = -1.0*temp;
  return temp;
}

void
XSLBToDtaunu_CLN::ComputeCLN()
{
  double w = (_mB*_mB+_mXu*_mXu-_q2)/(2.*_mB*_mXu);
  double z = (sqrt(w+1)-sqrt(2.))/(sqrt(w+1)+sqrt(2.));
  double v1 = _v1_1*(1.- 8.*_rho2*z + (51.*_rho2-10.)*z*z - (252.*_rho2-84.)*z*z*z);

  double wm1 = w - 1.0;
  double s1 = v1 * (1.0036 - 0.0068*wm1 + 0.0017*wm1*wm1 - 0.0013*wm1*wm1*wm1);

  cln_fplus = v1;
  cln_fminus = (s1 - v1) * (_mB*_mB-_mXu2) / _q2;
}

double 
XSLBToDtaunu_CLN::FromISGW2ToThisModel()
{
  //Rather than use an approximate normalization for ISGW2, i take advantage of the fact
  //that most of the constants (G_F, pi, Vcb) cancel when we take the ratio of
  //dGamma/d(kinvars) between CLN and ISGW2. The only part that doesn't cancel
  //completely is the product of currents L_munu*H^munu (see Eq 15 in K&S).
  //I simply calculate L*H for the two models, following Eq 18, with the appropriate
  //helicity projections. To do this, I need f+(q2) and f-(q2) for both models.

  double isgw2,cln;
  double isgw2_hl,isgw2_hs,isgw2_hsl,cln_hl,cln_hs,cln_hsl;
  double isgw2_h0,isgw2_ht,cln_h0,cln_ht;

  isgw2_h0  = 2.*_mB*kin_p*isgw2_fplus / sqrt(_q2);
  isgw2_ht  = ((_mB*_mB-_mXu2)*isgw2_fplus + _q2*isgw2_fminus) / sqrt(_q2);
  isgw2_hl  = isgw2_h0*isgw2_h0;
  isgw2_hs  = 3.*isgw2_ht*isgw2_ht;
  isgw2_hsl = isgw2_ht*isgw2_h0;
  isgw2 = 0.75*sin(_thL)*sin(_thL)*isgw2_hl +
    (kin_mu2/(2.*_q2)) * (1.5*cos(_thL)*cos(_thL)*isgw2_hl +
			  0.5*isgw2_hs + 3.0*cos(_thL)*isgw2_hsl);

  cln_h0  = 2.*_mB*kin_p*cln_fplus / sqrt(_q2);
  cln_ht  = ((_mB*_mB-_mXu2)*cln_fplus + _q2*cln_fminus) / sqrt(_q2);
  cln_hl  = cln_h0*cln_h0;
  cln_hs  = 3.*cln_ht*cln_ht;
  cln_hsl = cln_ht*cln_h0;
  cln = 0.75*sin(_thL)*sin(_thL)*cln_hl +
    (kin_mu2/(2.*_q2)) * (1.5*cos(_thL)*cos(_thL)*cln_hl +
			  0.5*cln_hs + 3.0*cos(_thL)*cln_hsl);
  return cln/isgw2/0.785;
  //factor 0.785 external normalization taken from EvtGen run
}

double 
XSLBToDtaunu_CLN::GetFplus(double q2)
{
  return 1.0;
}

void
XSLBToDtaunu_CLN::SetNormalizations(string mode)
{
  //i use a different norm strategy for the tau modes 
  _FLATQ2Normalization=1.0; 
  _PHSPNormalization=1.0; 
  _ISGW2Normalization=1.0;
}

double XSLBToDtaunu_CLN::FromPHSPToThisModel()
{ cout<<"XSLBToDtaunu_CLN::FromPHSPToThisModel not implemented!   returning 1.0"<<endl; return 1.0; }

double XSLBToDtaunu_CLN::FromFLATQ2ToThisModel()
{ cout<<"XSLBToDtaunu_CLN::FromFLATQ2ToThisModel not implemented!   returning 1.0"<<endl; return 1.0; }
