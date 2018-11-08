//------------------------------------------------------------------------
// File and Version Information:
//      $Id: RateCalc.hh,v 1.1 2012/03/04 00:33:02 manuelf Exp $
//
// Description:
//    Calculation of the B->DlNu rates based on CLN FF from hep-ph/1203.2654 
//    for a lepton of mass ml.
//    The angular dependence is from Korner, Shuler (1990).
//    Normalization is calculated for each set of FF values.
//    The New Physics dependence is parameterized in terms of
//    gSR = -mb (tanBeta/mH)^2 from hep-ph/1203.2654
//
//    FromSP8ToThisModel(q2, ctl,  ctv, chi, isDgamma, lplus, ml) 
//       Re-weights SP8 MC to the CLN parameterization.
//    SetMasses(isBm) 
//       Sets the masses of the B and D* depending on the charge.
//    Rate(isCLN, ml) 
//       Branching fraction.
//    Gamma_q2(q2, A1, V, A2, A0, ml) 
//       q2 spectrum.
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//      Michael Mazur                             INFN Pisa
//
// Revision History:
//      12/05/10 manuelf -- Normalization validated with EvtGen, including Higgs
//                          Added the theta spectrum integrated over q2 *NOT VALIDATED*
//      12/05/09 manuelf -- Corrected a missing sin(thetaL)sin(thetaV) and added D*->Dgamma
//      12/05/03 manuelf -- Added the NP dependence in terms of gSR
//      12/04/02 manuelf -- Created off XSLBToDstrtaunu_CLN.cc by M. Mazur
//------------------------------------------------------------------------

#include "ff_dstaunu.hpp"

BToDstaunu::BToDstaunu(double rho2, double R1, double R2, double R0, double gSR) {

  _rho2 = rho2;
  _R1 = R1; _R2 = R2; _R0 = R0; 
  _gSR   = gSR;
  SetMasses(1); // Charged B
}

void BToDstaunu::SetMasses(int isBm){
  if(isBm) {  // Charged B
    _mB  = 5.2792; _mBSP8 = 5.2791; _mDs = 2.00696; _mDsSP8 = 2.006620; _BLifeTime = 1.638e-12; 
  } else {   // Neutral B
    _mB  = 5.2795; _mBSP8 = 5.2794; _mDs = 2.01025; _mDsSP8 = 2.009936; _BLifeTime = 1.525e-12; 
  }
  _Dsmaxq2 = pow(_mB-_mDs,2);
  _isBm = isBm;
}

double BToDstaunu::FromSP8ToThisModel(double q2, double ctl, double ctv, double chi, int isDgamma, bool lplus, double ml) {

  double fix_gSR = _gSR; _gSR = 0;
  double SP8 = Compute(q2, ctl, ctv, chi, isDgamma, lplus, 0, ml);
  _gSR = fix_gSR;

  double cln = Compute(q2, ctl, ctv, chi, isDgamma, lplus, 1, ml);

  return cln/SP8/Normalization(ml);
}

double BToDstaunu::Compute(double q2, double ctl, double ctv, double chi, int isDgamma, bool lplus, int isCLN, double ml) {
  if(q2<= ml*ml || q2>=pow(_mB-_mDs,2)) return 0;
  double A1, V, A2, A0;
  if(isCLN) ComputeCLN(q2, A1, V, A2, A0);
  else{
    if(ml<mTau) ComputeLinearQ2(q2, A1, V, A2, A0);
    else        ComputeISGW2(q2, A1, V, A2, A0);
  }
  return Gamma_q2Angular(q2, ctl, ctv, chi, isDgamma, lplus, A1, V, A2, A0, ml);
}


double BToDstaunu::Gamma_q2Angular(double q2, double ctl, double ctv, double chi, int isDgamma, bool lplus,
				   double A1, double V, double A2, double A0, double ml){
  double pDs = _mB*_mB*_mB*_mB + _mDs*_mDs*_mDs*_mDs + q2*q2 -
    2.*_mB*_mB*_mDs*_mDs - 2.*_mB*_mB*q2 - 2.*_mDs*_mDs*q2;
  if (pDs < 0) pDs = 0;
  else pDs = sqrt(pDs)/(2.*_mB);

  double H0, Ht, Hplus, Hminus; 
  HadronicAmp(q2, A1, V, A2, A0,H0, Ht, Hplus, Hminus);

  //flip represents the spin flip term in Eq 24
  double flip  = ml*ml / (2.0 * q2);
  //these are just time-saving macros for sin(2*theta)
  double sin2tl = 2.*ctl*sqrt(1.-ctl*ctl);
  double sin2tv = 2.*ctv*sqrt(1.-ctv*ctv);
  //the 13 angular terms here correspond to the 13 lines of Eq 22, in the same order
  double ang1  = (3./8.)*(1.+ctl*ctl)*(3./4.)*(1.-ctv*ctv);                   // HV = H+2 + H-2
  double ang2  = (3./4.)*(1.-ctl*ctl)*(3./2.)*(ctv*ctv);		      // HL = H02
  double ang3  = (-3./4.)*(1.-ctl*ctl)*cos(2.*chi)*(3./4.)*(1.-ctv*ctv);      // HT = H+H-
  double ang4  = (-9./16.)*sin2tl*cos(chi)*sin2tv;			      // HI = (H+H0 + H-H0)/2
  double ang5  = (3./4.)*ctl*(3./4.)*(1.-ctv*ctv);			      // HP = H+2 - H-2
  double ang6  = (-9./8.)*sqrt(1.-ctl*ctl)*cos(chi)*sin2tv;		      // HA = (H+H0 - H-H0)/2
  double ang7  = (3./4.)*(1.-ctl*ctl)*(3./4.)*(1.-ctv*ctv);		      // HV = H+2 + H-2
  double ang8  = (3./2.)*ctl*ctl*(3./2.)*ctv*ctv;			      // HL = H02
  double ang9  = (3./4.)*(1.-ctl*ctl)*cos(2.*chi)*(3./4.)*(1.-ctv*ctv)*2.0;   // HT = H+H-
  double ang10 = (9./8.)*sin2tl*cos(chi)*sin2tv;			      // HI = (H+H0 + H-H0)/2
  double ang11 = (3./2.)*ctv*ctv*(1./2.);				      // HS = 3Ht2
  double ang12 = (3.)*ctl*(3./2.)*ctv*ctv;				      // HSL = HtH0
  double ang13 = (9./4.)*sqrt(1.-ctl*ctl)*cos(chi)*sin2tv;                    // HST = (H+Ht + H-Ht)/2
  if(isDgamma){
    ang1  = (3./8.)*(1.+ctl*ctl)*(3./4.)*(1.+ctv*ctv);                        // HV = H+2 + H-2	      
    ang2  = (3./4.)*(1.-ctl*ctl)*(3./2.)*(1.-ctv*ctv);			      // HL = H02	      
    ang3  = (3./4.)*(1.-ctl*ctl)*cos(2.*chi)*(3./4.)*(1.-ctv*ctv);	      // HT = H+H-	      
    ang4  = (9./16.)*sin2tl*cos(chi)*sin2tv;				      // HI = (H+H0 + H-H0)/2 
    ang5  = (3./4.)*ctl*(3./4.)*(1.+ctv*ctv);				      // HP = H+2 - H-2	      
    ang6  = (9./8.)*sqrt(1.-ctl*ctl)*cos(chi)*sin2tv;			      // HA = (H+H0 - H-H0)/2 
    ang7  = (3./4.)*(1.-ctl*ctl)*(3./4.)*(1.+ctv*ctv);			      // HV = H+2 + H-2	      
    ang8  = (3./2.)*ctl*ctl*(3./2.)*(1.-ctv*ctv);			      // HL = H02	      
    ang9  = (-3./4.)*(1.-ctl*ctl)*cos(2.*chi)*(3./4.)*(1.-ctv*ctv)*2.0;	      // HT = H+H-	      
    ang10 = (-9./8.)*sin2tl*cos(chi)*sin2tv;				      // HI = (H+H0 + H-H0)/2 
    ang11 = (3./2.)*(1.-ctv*ctv)*(1./2.);				      // HS = 3Ht2	      
    ang12 = (3.)*ctl*(3./2.)*(1.-ctv*ctv);				      // HSL = HtH0	      
    ang13 = (-9./4.)*sqrt(1.-ctl*ctl)*cos(chi)*sin2tv;			      // HST = (H+Ht + H-Ht)/2
  }
  if (lplus) {
    ang5 *= -1.0;
    ang6 *= -1.0;
  } //tau+ parity flip

  double Hu  = Hplus*Hplus + Hminus*Hminus;
  double Hl  = H0*H0;
  double Hp  = Hplus*Hplus - Hminus*Hminus;
  double Hs  = 3.0 * Ht*Ht;
  double Hsl = Ht*H0;
  double Hti = Hplus*Hminus;
  double Hi  = 0.5 * (Hplus*H0 + Hminus*H0);
  double Ha  = 0.5 * (Hplus*H0 - Hminus*H0);
  double Hst = 0.5 * (Hplus*Ht + Hminus*Ht);
  double gamma = ang1*Hu +
    ang2*Hl +
    ang3*Hti +
    ang4*Hi +
    ang5*Hp +
    ang6*Ha +
    ang7*flip*Hu +
    ang8*flip*Hl +
    ang9*flip*Hti +
    ang10*flip*Hi +
    ang11*flip*Hs +
    ang12*flip*Hsl +
    ang13*flip*Hst ;

  return GF*GF/pow(2*PI,3)*Vcb*Vcb/12./pow(_mB,2)*gamma*pow(q2-ml*ml,2)*pDs/q2;
}

void BToDstaunu::ComputeCLN(double q2, double &A1, double &V, double &A2, double &A0) {
  double w = (_mB*_mB+_mDs*_mDs-q2)/(2.*_mB*_mDs);
  double z = (sqrt(w+1.)-sqrt(2.))/(sqrt(w+1.)+sqrt(2.));
  double hA1 = F1*(1.- 8.*_rho2*z + (53.*_rho2-15.)*z*z - (231.*_rho2-91.)*z*z*z);
  double wm1 = w - 1.0;
  double RDs = 2*sqrt(_mB*_mDs)/(_mB+_mDs);

  A1 = hA1* RDs*(w+1)/2.;
  V  = hA1* (_R1 - 0.12*wm1 + 0.05*wm1*wm1) / RDs;
  A2 = hA1* (_R2 + 0.11*wm1 - 0.06*wm1*wm1) / RDs;
  A0 = hA1* (_R0 - 0.11*wm1 + 0.01*wm1*wm1) / RDs;
}

// The total rate for LinearQ2 is off because I do not include F1 here (to keep the code from EvtGen)
void BToDstaunu::ComputeLinearQ2(double q2, double &A1, double &V, double &A2, double &A0) {
  double w = (_mBSP8*_mBSP8+_mDsSP8*_mDsSP8-q2)/(2.*_mBSP8*_mDsSP8);
  double RDs = 2*sqrt(_mBSP8*_mDsSP8)/(_mBSP8+_mDsSP8);
  double rho2_SP8 = 0.77, R1_SP8 = 1.33, R2_SP8 = 0.92;
  double hA1 = 1-rho2_SP8*(w-1);

  A1 = hA1* RDs*(w+1)/2.;
  V  = hA1* R1_SP8/RDs;
  A2 = hA1* R2_SP8/RDs;
  A0 = 0.;
}

void BToDstaunu::ComputeISGW2(double q2, double &A1, double &V, double &A2, double &A0) {
  //This code liberally stolen from EvtGenModels/EvtISGW2FF.cc
  double msb=5.2;
  double msd=0.33;
  double bb2=0.431*0.431;
  double mbb=5.31;
  double nf = 4.0;
  double cf=0.989;
  double msq=1.82;
  double bx2=0.38*0.38;
  double mbx=0.75*2.01+0.25*1.87;
  double nfp = 3.0;
  double mtb=msb+msd;
  double mtx=msq+msd;
  double mup=1.0/(1.0/msq+1.0/msb);
  double mum=1.0/(1.0/msq-1.0/msb);
  double bbx2=0.5*(bb2+bx2);
  double mb=_mBSP8;
  double mx=_mDsSP8; 
  double tm=(mb-mx)*(mb-mx);
  double t = q2;
  if ( t > tm ) t = 0.99*tm;
  double wt=1.0+(tm-t)/(2.0*mbb*mbx);
  double mqm = 0.1;
  double r2=3.0/(4.0*msb*msq)+3*msd*msd/(2*mbb*mbx*bbx2) + 
    (16.0/(mbb*mbx*(33.0-2.0*nfp)))*
    log(EvtGetas(mqm,mqm)/EvtGetas(msq,msq));
  double ai = -1.0* ( 6.0/( 33.0 - 2.0*nf));  
  double cji = pow(( EvtGetas( msb,msb ) / EvtGetas( msq,msq ) ),ai);
  double zji = msq / msb;
  double gammaji = EvtGetGammaji( zji );
  double chiji = -1.0 - ( gammaji / ( 1- zji ));
  double betaji_g = (2.0/3.0)+gammaji;
  double betaji_f = (-2.0/3.0)+gammaji;
  double betaji_appam = -1.0-chiji+(4.0/(3.0*(1.0-zji)))+
    (2.0*(1+zji)*gammaji/(3.0*(1.0-zji)*(1.0-zji)));
  double betaji_apmam = (1.0/3.0)-chiji-(4.0/(3.0*(1.0-zji)))-
    (2.0*(1+zji)*gammaji/(3.0*(1.0-zji)*(1.0-zji)))+
    gammaji;
  double r_g = cji*(1+(betaji_g*EvtGetas( msq,sqrt(mb*msq) )/(PI)));
  double r_f = cji*(1+(betaji_f*EvtGetas( msq,sqrt(mb*msq) )/(PI)));
  double r_apmam = cji*(1+(betaji_apmam*EvtGetas( msq,sqrt(mb*msq) )/(PI)));
  double f3=sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,1.5)/
    ((1.0+r2*(tm-t)/12.0)*(1.0+r2*(tm-t)/12.0));
  double f3f=sqrt(mbx*mbb/(mtx*mtb))*f3;
  double f3g=sqrt(mtx*mtb/(mbx*mbb))*f3;
  double f3appam=sqrt(mtb*mtb*mtb*mbx/(mbb*mbb*mbb*mtx))*f3;
  double f3apmam=sqrt(mtx*mtb/(mbx*mbb))*f3;
  double ff=cf*mtb*(1+wt+msd*(wt-1)/(2*mup))*f3f*r_f;
  double gf=0.5*(1/msq-msd*bb2/(2*mum*mtx*bbx2))*f3g*r_g;
  double appam=cji*(msd*bx2*(1-msd*bx2/(2*mtb*bbx2))/ 
		    ((1+wt)*msq*msb*bbx2)-
		    betaji_appam*EvtGetas( msq,sqrt(msq*mb) )/
		    (mtb*PI))*f3appam;
  double apmam=-1.0*(mtb/msb-msd*bx2/(2*mup*bbx2)+wt*msd*mtb*bx2*
		     (1-msd*bx2/(2*mtb*bbx2))/((wt+1)*msq*msb*bbx2))*
    f3apmam*r_apmam/mtx;
  double apf=0.5*(appam+apmam);
  double amf=0.5*(appam-apmam);
  double mass = _mDs;
  V = (gf)*(mb+mass);
  A1 = (ff)/(mb+mass);
  A2 = -1.0*(apf)*(mb+mass);
  double a3f = ((mb+mass)/(2.0*mass))*(A1) - ((mb-mass)/(2.0*mass))*(A2);
  A0 = a3f + ( (t*amf)/(2.0*mass));
}

double BToDstaunu::EvtGetas(double massq, double massx) {
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

double BToDstaunu::EvtGetGammaji(double z) {
  double temp;
  temp = 2+((2.0*z)/(1-z))*log(z);
  temp = -1.0*temp;
  return temp;
}


// Calculated with Polynomials(ml)
// The difference between B- and B0 is 0.2% for Tau and less than 0.005 for mu/e, so it is a bit of an overkill
double BToDstaunu::Normalization(double ml){
  double CoefEBm[3][7] = {{0.133038, -2.74008e-09, 2.68753e-08, -0.000194059, 0.00326472, -0.0716297, 0.0135873}, 
			  {-0.0621645, 1.78822e-09, -1.74019e-08, 9.18807e-05, -0.00140314, 0.0386214, -0.00761311}, 
			  {0.00808777, -2.92815e-10, 2.8352e-09, -1.14938e-05, 0.000164962, -0.00543084, 0.00109709}};
  double CoefMuBm[3][7] = {{0.130545, -4.79136e-05, 0.000494882, -0.000192699, 0.00324618, -0.06957, 0.0131328}, 
			   {-0.0606082, 3.02863e-05, -0.000307856, 9.11299e-05, -0.00139347, 0.0373018, -0.00732002}, 
			   {0.00784021, -4.82226e-06, 4.85109e-05, -1.13868e-05, 0.000163633, -0.00521798, 0.00104961}};
  double CoefTauBm[3][7] = {{0.0180312, -6.26354e-05, 0.00136444, -2.1904e-05, 0.000545498, -0.00422738, 0.000522854}, 
			    {-0.00484354, 2.33609e-05, -0.00046664, 7.39039e-06, -0.000163964, 0.00141692, -0.000188988}, 
			    {0.000388608, -2.28211e-06, 4.30657e-05, -6.66758e-07, 1.37327e-05, -0.000128625, 1.80067e-05}};
  double CoefEB0[3][7] = {{0.123304, -2.53786e-09, 2.49757e-08, -0.00017901, 0.00302065, -0.0662345, 0.0125426}, 
			  {-0.0574465, 1.65215e-09, -1.61318e-08, 8.45318e-05, -0.00129468, 0.0356172, -0.0070095}, 
			  {0.00745349, -2.69863e-10, 2.62176e-09, -1.05471e-05, 0.000151808, -0.0049954, 0.00100752}};
  double CoefMuB0[3][7] = {{0.120993, -4.43569e-05, 0.000459713, -0.000177753, 0.00300347, -0.0643277, 0.0121224}, 
			   {-0.056007, 2.79674e-05, -0.000285253, 8.38397e-05, -0.00128574, 0.0343987, -0.00673927}, 
			   {0.0072251, -4.44187e-06, 4.48358e-05, -1.04487e-05, 0.000150583, -0.00479934, 0.000963856}};
  double CoefTauB0[3][7] = {{0.0166738, -5.76492e-05, 0.00126155, -2.01117e-05, 0.000502887, -0.00389332, 0.000480201}, 
			    {-0.00446049, 2.14227e-05, -0.000429823, 6.76147e-06, -0.000150599, 0.0013, -0.000172938}, 
			    {0.000356511, -2.08534e-06, 3.95239e-05, -6.07896e-07, 1.25685e-05, -0.000117581, 1.6419e-05}};
  double RateSP8=0, (*Coef)[7]=0;
  if(_isBm){
    if(ml<0.01) {     // Electron
      Coef = CoefEBm;
      RateSP8 = 0.0641614;  
    } else if(ml<1){  // Muon
      Coef = CoefMuBm;
      RateSP8 = 0.0635797;  
    } else {          // Tau
      Coef = CoefTauBm;
      RateSP8 = 0.0135739;  
    }
  } else {
    if(ml<0.01) {     // Electron
      Coef = CoefEB0;
      RateSP8 = 0.0596401;  
    } else if(ml<1){  // Muon
      Coef = CoefMuB0;
      RateSP8 = 0.0590979;  
    } else {          // Tau
      Coef = CoefTauB0;
      RateSP8 = 0.0125446;  
    }
  }
  double RateCLN;
  if(fabs(_gSR)<1e-6) RateCLN = SumPoly(Coef);  // The normalization is pre-calculated for the SM
  else RateCLN = Rate(1,ml);                    // With NP, the rate has to be integrated
  return RateCLN/RateSP8;
}

double BToDstaunu::Compute(double q2, int isCLN, double ml) {
  if(q2<= ml*ml || q2>=pow(_mB-_mDs,2)) return 0;
  double A1, V, A2, A0;
  if(isCLN) ComputeCLN(q2, A1, V, A2, A0);
  else{
    if(ml<mTau) ComputeLinearQ2(q2, A1, V, A2, A0);
    else        ComputeISGW2(q2, A1, V, A2, A0);
  }

  return Gamma_q2(q2, A1, V, A2, A0, ml);
}

void BToDstaunu::HadronicAmp(double q2, double A1, double V, double A2, double A0,
			     double &H0, double &Ht, double &Hplus, double &Hminus){
  double pDs = _mB*_mB*_mB*_mB + _mDs*_mDs*_mDs*_mDs + q2*q2 -
    2.*_mB*_mB*_mDs*_mDs - 2.*_mB*_mB*q2 - 2.*_mDs*_mDs*q2;
  if (pDs < 0) pDs = 0;
  else pDs = sqrt(pDs)/(2.*_mB);

  H0 = ((_mB*_mB-_mDs*_mDs-q2)*(_mB+_mDs)*A1 - 4*_mB*_mB*pDs*pDs*A2/(_mB+_mDs))/(2*_mDs*sqrt(q2));
  Ht = 2*_mB*pDs*A0/sqrt(q2)*(1+_gSR*q2/(mb_quark+mc_quark));
  Hplus  = (_mB+_mDs)*A1 + 2*_mB/(_mB+_mDs)*pDs*V; // The sign in the middle is reversed from 1203.2654
  Hminus = (_mB+_mDs)*A1 - 2*_mB/(_mB+_mDs)*pDs*V; // to match Korner-Shuler (acording to Mazur)
}


// The q2 spectrum (no thetaL) is only used to integrate the rates
double BToDstaunu::Gamma_q2(double q2, double A1, double V, double A2, double A0, double ml){
  double pDs = _mB*_mB*_mB*_mB + _mDs*_mDs*_mDs*_mDs + q2*q2 -
    2.*_mB*_mB*_mDs*_mDs - 2.*_mB*_mB*q2 - 2.*_mDs*_mDs*q2;
  if (pDs < 0) pDs = 0;
  else pDs = sqrt(pDs)/(2.*_mB);

  double H0, Ht, Hplus, Hminus; 
  HadronicAmp(q2, A1, V, A2, A0,H0, Ht, Hplus, Hminus);

  double Term1 = (Hplus*Hplus+Hminus*Hminus+H0*H0)*(1+ml*ml/(2*q2));
  double Term2 = 3*ml*ml/(2*q2)*Ht*Ht;
  double Factor = GF*GF*Vcb*Vcb*pDs*q2/(96*pow(PI,3)*_mB*_mB)*pow(1-ml*ml/q2,2);
  return Factor*(Term1+Term2);
}

double BToDstaunu::ComputetL(double ctl, int isCLN, double ml) {
  return IntRate(ml*ml,_Dsmaxq2,isCLN,ctl,ml,1000);
}

double BToDstaunu::Compute(double q2, double ctl, int isCLN, double ml) {
  if(q2<= ml*ml || q2>=pow(_mB-_mDs,2)) return 0;
  double A1, V, A2, A0;
  if(isCLN) ComputeCLN(q2, A1, V, A2, A0);
  else{
    if(ml<mTau) ComputeLinearQ2(q2, A1, V, A2, A0);
    else        ComputeISGW2(q2, A1, V, A2, A0);
  }

  return Gamma_q2tL(q2, ctl, A1, V, A2, A0, ml);
}

// This spectrum uses the formula from hep-ph/1203.2654, fixing the sign in Ht+H0*ctl
double BToDstaunu::Gamma_q2tL(double q2, double ctl, double A1, double V, double A2, double A0, double ml){
  double pDs = _mB*_mB*_mB*_mB + _mDs*_mDs*_mDs*_mDs + q2*q2 -
    2.*_mB*_mB*_mDs*_mDs - 2.*_mB*_mB*q2 - 2.*_mDs*_mDs*q2;
  if (pDs < 0) pDs = 0;
  else pDs = sqrt(pDs)/(2.*_mB);

  double H0, Ht, Hplus, Hminus, stl = sqrt(1-ctl*ctl); 
  HadronicAmp(q2, A1, V, A2, A0,H0, Ht, Hplus, Hminus);

  double Term1 = pow(1-ctl,2)*Hplus*Hplus + pow(1+ctl,2)*Hminus*Hminus + 2*stl*stl*H0*H0;
  double Term2 = ml*ml/q2*(stl*stl*(Hplus*Hplus+Hminus*Hminus) + 2*pow(Ht+H0*ctl,2));
  double Factor = GF*GF*Vcb*Vcb*pDs*q2/(256*pow(PI,3)*_mB*_mB)*pow(1-ml*ml/q2,2);
  return Factor*(Term1+Term2);
}

// Simpson integration
double BToDstaunu::IntRate(double minX, double maxX, int isCLN, double ctl, double ml, int nPoints){
  double intF = 0, x = minX, dx = (maxX-minX)/static_cast<double>(nPoints);
  double Fmin, Fval;
  if(ctl<-50) Fmin = Compute(x, isCLN, ml);
  else        Fmin = Compute(x, ctl, isCLN, ml);
  for(int step=0; step<nPoints; step++){
    if(ctl<-50) Fval = Compute(x+dx/2, isCLN, ml);
    else        Fval = Compute(x+dx/2, ctl, isCLN, ml);
    intF += dx/6.*(Fmin+4*Fval);
    x += dx;
    if(ctl<-50) Fmin = Compute(x, isCLN, ml);
    else        Fmin = Compute(x, ctl, isCLN, ml);
    intF +=dx/6.*Fmin;
  }
  
  return intF;
}

// Decay rate of B->D*lnu with respect to the total rate
double BToDstaunu::Rate(int isCLN, double ml) {
  double totalRate = hbar/_BLifeTime;
  return IntRate(ml*ml,_Dsmaxq2,isCLN,-99,ml)/totalRate;
}


// The integral of the rate is a polynomial of order 2 in _rho2 and _R0,_R1,_R2
// The coefficients are found once, and stored in Normalization
void BToDstaunu::Polynomial(double ml) {
  double F[3], Coef[3][7], fixRho2 = _rho2, fixR0 = _R0, fixR1 = _R1, fixR2 = _R2;

  _R0 = 0; _R1 = 0; _R2 = 0; 
  _rho2 = -1; F[0] = Rate(1,ml);
  _rho2 =  0; F[1] = Rate(1,ml);
  _rho2 =  1; F[2] = Rate(1,ml);
  Coef[0][0] = F[1];
  Coef[1][0] = (F[2]-F[0])/2;
  Coef[2][0] = -F[1]+(F[2]+F[0])/2;

  // Terms in _R0
  _rho2 = 0; _R1 = 0; _R2 = 0; 
  _R0 = -1; F[0] = Rate(1,ml);
  _R0 =  1; F[2] = Rate(1,ml);
  Coef[0][1] = (F[2]-F[0])/2;
  Coef[0][2] = -F[1]+(F[2]+F[0])/2;

  _rho2 =  1; _R0 =  1; double F11  = Rate(1,ml)-SubPoly(Coef, 0);
  _rho2 = -1; _R0 =  1; double Fm11 = Rate(1,ml)-SubPoly(Coef, 0);
  _rho2 =  1; _R0 = -1; double F1m1 = Rate(1,ml)-SubPoly(Coef, 0);
  _rho2 =  2; _R0 = -1; double F2m1 = Rate(1,ml)-SubPoly(Coef, 0);

  Coef[1][1] = (F11-4*F1m1+F2m1-Fm11)/4;
  Coef[2][1] = (F11+2*F1m1-F2m1+Fm11)/4;
  Coef[1][2] = (F11+4*F1m1-F2m1-Fm11)/4;
  Coef[2][2] = (F11-2*F1m1+F2m1+Fm11)/4;

  // Terms in _R1
  _rho2 = 0; _R0 = 0; _R2 = 0; 
  _R1 = -1; F[0] = Rate(1,ml);
  _R1 =  1; F[2] = Rate(1,ml);
  Coef[0][3] = (F[2]-F[0])/2;
  Coef[0][4] = -F[1]+(F[2]+F[0])/2;

  _rho2 =  1; _R1 =  1; F11  = Rate(1,ml)-SubPoly(Coef, 1);
  _rho2 = -1; _R1 =  1; Fm11 = Rate(1,ml)-SubPoly(Coef, 1);
  _rho2 =  1; _R1 = -1; F1m1 = Rate(1,ml)-SubPoly(Coef, 1);
  _rho2 =  2; _R1 = -1; F2m1 = Rate(1,ml)-SubPoly(Coef, 1);

  Coef[1][3] = (F11-4*F1m1+F2m1-Fm11)/4;
  Coef[2][3] = (F11+2*F1m1-F2m1+Fm11)/4;
  Coef[1][4] = (F11+4*F1m1-F2m1-Fm11)/4;
  Coef[2][4] = (F11-2*F1m1+F2m1+Fm11)/4;

  // Terms in _R2
  _rho2 = 0; _R0 = 0; _R1 = 0; 
  _R2 = -1; F[0] = Rate(1,ml);
  _R2 =  1; F[2] = Rate(1,ml);
  Coef[0][5] = (F[2]-F[0])/2;
  Coef[0][6] = -F[1]+(F[2]+F[0])/2;

  _rho2 =  1; _R2 =  1; F11  = Rate(1,ml)-SubPoly(Coef, 2);
  _rho2 = -1; _R2 =  1; Fm11 = Rate(1,ml)-SubPoly(Coef, 2);
  _rho2 =  1; _R2 = -1; F1m1 = Rate(1,ml)-SubPoly(Coef, 2);
  _rho2 =  2; _R2 = -1; F2m1 = Rate(1,ml)-SubPoly(Coef, 2);

  Coef[1][5] = (F11-4*F1m1+F2m1-Fm11)/4;
  Coef[2][5] = (F11+2*F1m1-F2m1+Fm11)/4;
  Coef[1][6] = (F11+4*F1m1-F2m1-Fm11)/4;
  Coef[2][6] = (F11-2*F1m1+F2m1+Fm11)/4;


  _rho2=1.214; _R1=1.401; _R2=0.864; _R0=1.1387; // Check that the calculation is correct
  cout<<"Rate is "<<Rate(1,ml)<<"\t Polynomial yields "<<SumPoly(Coef)<<endl;

  cout<<endl<<"double RateSP8 = "<<Rate(0,ml)<<";"<<endl;
  cout<<endl<<"double Coef[3][7] = {";
  for(int nRho=0; nRho<3; nRho++){
    cout<<"{";
    for(int nR=0; nR<7; nR++) {
      cout<<Coef[nRho][nR];
      if(nR<6)cout<<", ";
    }
    cout<<"}";
    if(nRho<2)cout<<", "<<endl;
  }
  cout<<"};"<<endl<<endl;


  _rho2 = fixRho2; _R0 = fixR0; _R1 = fixR1; _R2 = fixR2;
}

double BToDstaunu::SumPoly(double Coef[3][7]){
  double poly = 0;

  for(int nRho=0; nRho<3; nRho++){
    for(int nR0=0; nR0<3; nR0++)
      poly += Coef[nRho][nR0]*pow(_rho2,nRho)*pow(_R0,nR0);
    for(int nR1=1; nR1<3; nR1++)
      poly += Coef[nRho][2+nR1]*pow(_rho2,nRho)*pow(_R1,nR1);
    for(int nR2=1; nR2<3; nR2++)
      poly += Coef[nRho][4+nR2]*pow(_rho2,nRho)*pow(_R2,nR2);
  }

  return poly;
}

double BToDstaunu::SubPoly(double C[3][7], int indexR){
  if(indexR == 0)
    return C[0][0]+C[1][0]*_rho2+C[2][0]*_rho2*_rho2+C[0][1]*_R0+C[0][2]*_R0*_R0;
  else if(indexR == 1)
    return C[0][0]+C[1][0]*_rho2+C[2][0]*_rho2*_rho2+C[0][3]*_R1+C[0][4]*_R1*_R1;
  else if(indexR == 2)
    return C[0][0]+C[1][0]*_rho2+C[2][0]*_rho2*_rho2+C[0][5]*_R2+C[0][6]*_R2*_R2;

  cout<<"Error in SubPoly"<<endl;
  return 0;
}

BToDstaunu::~BToDstaunu() {
}
