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

#include "babar_code/FF/BToDtaunu.hh"

BToDtaunu::BToDtaunu(double rho2, double Delta, double gSR) {
  SetFF(rho2, Delta, gSR);
  SetMasses(1); // Charged B
}

void BToDtaunu::SetMasses(int isBm){
  if(isBm) { // Charged B
    _mB  = 5.2792; _mBISGW2 = 5.2791; _mD = 1.8648; _mDISGW2 = 1.8645; _BLifeTime = 1.638e-12;
  } else { // Neutral B
    _mB  = 5.2795; _mBISGW2 = 5.2794; _mD = 1.8696; _mDISGW2 = 1.8693; _BLifeTime = 1.525e-12;
  }
  _Dmaxq2 = pow(_mB-_mD,2);
  _isBm = isBm;
}

double BToDtaunu::FromSP8ToThisModel(double q2, double thetaL, double ml) {

  double fix_gSR = _gSR; _gSR = 0;
  double isgw2 = Compute(q2, thetaL,0, ml);
  _gSR = fix_gSR;

  double cln = Compute(q2, thetaL,1, ml);

  return cln/isgw2/Normalization(ml);
}

double BToDtaunu::Compute(double q2, double thetaL, int isCLN, double ml) {
  if(q2<= ml*ml || q2>=pow(_mB-_mD,2)) return 0;

  double fplus, fminus;
  if(isCLN) ComputeCLN(q2, fplus, fminus);
  else    ComputeISGW2(q2, fplus, fminus);

  return Gamma_q2tL(q2, thetaL, fplus, fminus, ml);
}

// New physics from Tanaka 2010
double BToDtaunu::Gamma_q2tL(double q2, double thetaL, double fplus, double fminus, double ml){
  double hl,hs,hsl,h0,ht;
  double ml2 = ml*ml;
  double pD = _mB*_mB*_mB*_mB + _mD*_mD*_mD*_mD + q2*q2 -
    2.*_mB*_mB*_mD*_mD - 2.*_mB*_mB*q2 - 2.*_mD*_mD*q2;
  if (pD < 0) pD = 0;
  else pD = sqrt(pD)/(2.*_mB);

  h0  = 2.*_mB*pD*fplus / sqrt(q2);
  ht  = ((_mB*_mB-_mD*_mD)*fplus + q2*fminus) / sqrt(q2)*(1+_gSR*q2/(mb_quark-mc_quark));
  hl  = h0*h0;
  hs  = 3.*ht*ht;
  hsl = ht*h0;
  double LH = 2./3*(q2-ml2)*(0.75*sin(thetaL)*sin(thetaL)*hl + (ml2/(2.*q2)) * 
 				 (1.5*cos(thetaL)*cos(thetaL)*hl +  0.5*hs + 3.0*cos(thetaL)*hsl));

  //cout<<"h0 "<<h0<<", ht "<<ht<<", NP factor "<<(1+_gSR*q2/(mb_quark-mc_quark))<<", LH "<<LH<<endl;
  return GF*GF/pow(2*PI,3)*Vcb*Vcb*(q2-ml2)*pD/(8*_mB*_mB*q2)*LH;
}

// S1 taken from Tanaka, Watanabe 2010
void BToDtaunu::ComputeCLN(double q2, double &fplus, double &fminus) {
  double w = (_mB*_mB+_mD*_mD-q2)/(2.*_mB*_mD);
  double z = (sqrt(w+1)-sqrt(2.))/(sqrt(w+1)+sqrt(2.));
  double V1 = V11*(1.- 8.*_rho2*z + (51.*_rho2-10.)*z*z - (252.*_rho2-84.)*z*z*z);

  double wm1 = w - 1.0;
  double S1 = V1 * (1 +_Delta*(-0.019 + 0.041*wm1 - 0.015*wm1*wm1));

  double RD = (2*sqrt(_mB*_mD))/(_mB+_mD);
  double F1 = V1/RD;
  double F0 = (1-q2/pow(_mB+_mD,2))/RD*S1;

  fplus = F1;
  fminus = (F0 - F1) * (_mB*_mB-_mD*_mD) / q2;
}

void BToDtaunu::ComputeISGW2(double q2, double &fplus, double &fminus) {
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
  double mb = _mBISGW2;      // B from pdg.table
  double mx = _mDISGW2;      // D from pdg.table
  double mup=1.0/(1.0/msq+1.0/msb);
  double bbx2=0.5*(bb2+bx2);
  double tm=(mb-mx)*(mb-mx);
  double t = q2;
  if (t > tm) t=0.99*tm;
  double mqm = 0.1;
  double r2=3.0/(4.0*msb*msq)+3*msd*msd/(2*mbb*mbx*bbx2) + 
    (16.0/(mbb*mbx*(33.0-2.0*nfp)))*
    log(EvtGetas(mqm,mqm)/EvtGetas(msq,msq));
  double f3 = sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,1.5) /
    ((1.0+r2*(tm-t)/12.0)*(1.0+r2*(tm-t)/12.0));
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
  fplus = (fppfm + fpmfm)/2.0;
  fminus = (fppfm - fpmfm)/2.0;
}

double BToDtaunu::EvtGetas(double massq, double massx) {
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

double BToDtaunu::EvtGetGammaji(double z) {
  double temp;
  temp = 2+((2.0*z)/(1-z))*log(z);
  temp = -1.0*temp;
  return temp;
}

// Calculated with Polynomials(ml)
// The difference between B- and B0 is 0.2% for Tau and less than 0.005 for mu/e, so it is a bit of an overkill
double BToDtaunu::Normalization(double ml){
  double CoefEBm[3][3] = {{0.0534102, -8.07824e-11, 2.61002e-13}, {-0.0318507, 4.15612e-11, -1.00026e-13}, 
			  {0.00499707, -6.00037e-12, 1.17653e-14}};
  double CoefMuBm[3][3] = {{0.0530911, -3.42366e-06, 1.11261e-08}, {-0.0315869, 1.758e-06, -4.26344e-09}, 
			   {0.00494503, -2.53443e-07, 4.95416e-10}};
  double CoefTauBm[3][3] = {{0.0113728, -0.000140859, 8.32051e-07}, {-0.00416981, 4.36796e-05, -2.25079e-07}, 
			    {0.000428858, -4.05144e-06, 1.89768e-08}};
  double CoefEB0[3][3] = {{0.0494652, -8.13116e-11, 2.47308e-13}, {-0.0293821, 4.30308e-11, -9.51478e-14}, 
			  {0.00459203, -6.37326e-12, 1.06442e-14}};
  double CoefMuB0[3][3] = {{0.0491692, -3.30057e-06, 1.05176e-08}, {-0.0291383, 1.71053e-06, -4.04567e-09}, 
			   {0.00454409, -2.48521e-07, 4.71801e-10}};
  double CoefTauB0[3][3] = {{0.0104972, -0.000130837, 7.75248e-07}, {-0.0038271, 4.0385e-05, -2.08983e-07}, 
			    {0.000391533, -3.72895e-06, 1.75561e-08}};
  double RateSP8=0, (*Coef)[3]=0;
  if(_isBm){
    if(ml<0.01) {     // Electron
      Coef = CoefEBm;
      RateSP8 = 0.0298156;  
    } else if(ml<1){  // Muon
      Coef = CoefMuBm;
      RateSP8 = 0.0296754;  
    } else {          // Tau
      Coef = CoefTauBm;
      RateSP8 = 0.00769757;  
    }
  } else {
    if(ml<0.01) {     // Electron
      Coef = CoefEB0;
      RateSP8 = 0.0276774;  
    } else if(ml<1){  // Muon
      Coef = CoefMuB0;
      RateSP8 = 0.027547;  
    } else {          // Tau
      Coef = CoefTauB0;
      RateSP8 = 0.00711795;  
    }
  }
  double RateCLN;
  if(fabs(_gSR)<1e-6) RateCLN = SumPoly(Coef);  // The normalization is pre-calculated for the SM
  else RateCLN = Rate(1,ml);                    // With NP, the rate has to be integrated

  return RateCLN/RateSP8;
}

double BToDtaunu::ComputetL(double thetaL, int isCLN, double ml) {
  return IntRate(ml*ml,_Dmaxq2,isCLN,thetaL,ml,1000);
}

double BToDtaunu::Compute(double q2, int isCLN, double ml) {
  if(q2<= ml*ml || q2>=pow(_mB-_mD,2)) return 0;
  double fplus, fminus;
  if(isCLN) ComputeCLN(q2, fplus, fminus);
  else ComputeISGW2(q2, fplus, fminus);

  return Gamma_q2(q2, fplus, fminus, ml);
}

// Simpson integration
double BToDtaunu::IntRate(double minX, double maxX, int isCLN, double thetaL, double ml, int nPoints){
  double intF = 0, x = minX, dx = (maxX-minX)/(double)nPoints;
  double Fmin, Fval;
  if(thetaL<-50) Fmin = Compute(x, isCLN, ml);
  else           Fmin = Compute(x, thetaL, isCLN, ml);
  for(int step=0; step<nPoints; step++){
    if(thetaL<-50) Fval = Compute(x+dx/2, isCLN, ml);
    else           Fval = Compute(x+dx/2, thetaL, isCLN, ml);
    intF += dx/6.*(Fmin+4*Fval);
    x += dx;
    if(thetaL<-50) Fmin = Compute(x, isCLN, ml);
    else           Fmin = Compute(x, thetaL, isCLN, ml);
    intF +=dx/6.*Fmin;
  }
  
  return intF;
}

// Decay rate of B->Dlnu with respect to the total rate
double BToDtaunu::Rate(int isCLN, double ml) {
  double totalRate = hbar/_BLifeTime;
  return IntRate(ml*ml,_Dmaxq2,isCLN,-99,ml)/totalRate;
}

// The q2 spectrum (no thetaL) is only used to integrate the rates
// New physics from Tanaka 2010
double BToDtaunu::Gamma_q2(double q2, double fplus, double fminus, double ml){
  double hl,hs,hsl,h0,ht;
  double ml2 = ml*ml;
  double pD = _mB*_mB*_mB*_mB + _mD*_mD*_mD*_mD + q2*q2 -
    2.*_mB*_mB*_mD*_mD - 2.*_mB*_mB*q2 - 2.*_mD*_mD*q2;
  if (pD < 0) pD = 0;
  else pD = sqrt(pD)/(2.*_mB);

  h0  = 2.*_mB*pD*fplus / sqrt(q2);
  ht  = ((_mB*_mB-_mD*_mD)*fplus + q2*fminus) / sqrt(q2)*(1+_gSR*q2/(mb_quark-mc_quark));
  hl  = h0*h0;
  hs  = 3.*ht*ht;
  hsl = ht*h0;
  double LH = 2./3*(q2-ml2) * (hl+(ml2/(2.*q2))*(hl+hs));
  return GF*GF/pow(2*PI,3)*Vcb*Vcb*(q2-ml2)*pD/(8*_mB*_mB*q2)*LH;
}

// The integral of the rate is a polynomial of order 2 in _rho2 and _Delta
// The coefficients are found once, and stored in Normaliztion
void BToDtaunu::Polynomial(double ml) {
  double F[3], Coef[3][3], fixRho2 = _rho2, fixDelta = _Delta;

  _Delta = 0;
  _rho2 = -1; F[0] = Rate(1,ml);
  _rho2 =  0; F[1] = Rate(1,ml);
  _rho2 =  1; F[2] = Rate(1,ml);
  Coef[0][0] = F[1];
  Coef[1][0] = (F[2]-F[0])/2;
  Coef[2][0] = -F[1]+(F[2]+F[0])/2;

  _rho2 = 0;
  _Delta = -1; F[0] = Rate(1,ml);   
  _Delta =  0; F[1] = Rate(1,ml);   // Redundant, I know. For symmetry!
  _Delta =  1; F[2] = Rate(1,ml);
  Coef[0][0] = F[1];
  Coef[0][1] = (F[2]-F[0])/2;
  Coef[0][2] = -F[1]+(F[2]+F[0])/2;

  _rho2 =  1; _Delta =  1; double F11  = Rate(1,ml)-SubPoly(Coef);
  _rho2 = -1; _Delta =  1; double Fm11 = Rate(1,ml)-SubPoly(Coef);
  _rho2 =  1; _Delta = -1; double F1m1 = Rate(1,ml)-SubPoly(Coef);
  _rho2 =  2; _Delta = -1; double F2m1 = Rate(1,ml)-SubPoly(Coef);

  Coef[1][1] = (F11-4*F1m1+F2m1-Fm11)/4;
  Coef[2][1] = (F11+2*F1m1-F2m1+Fm11)/4;
  Coef[1][2] = (F11+4*F1m1-F2m1-Fm11)/4;
  Coef[2][2] = (F11-2*F1m1+F2m1+Fm11)/4;

  _rho2=1.18; _Delta=1; // Check that the calculation is correct
  cout<<"Rate is "<<Rate(1,ml)<<"\t Polynomial yields "<<SumPoly(Coef)<<endl;

  cout<<endl<<"double RateSP8 = "<<Rate(0,ml)<<";"<<endl;
  cout<<endl<<"double Coef[3][3] = {";

  for(int nRho=0; nRho<3; nRho++){
    cout<<"{";
    for(int nDelta=0; nDelta<3; nDelta++) {
      cout<<Coef[nRho][nDelta];
      if(nDelta<2)cout<<", ";
    }
    cout<<"}";
    if(nRho<2)cout<<", ";
  }
  cout<<"};"<<endl<<endl;

  _rho2 = fixRho2; _Delta = fixDelta;
}

double BToDtaunu::SumPoly(double Coef[3][3]){
  double poly = 0;
  
  for(int nRho=0; nRho<3; nRho++)
    for(int nDelta=0; nDelta<3; nDelta++)
      poly += Coef[nRho][nDelta]*pow(_rho2,nRho)*pow(_Delta,nDelta);

  return poly;
}

double BToDtaunu::SubPoly(double C[3][3]){
  return C[0][0]+C[1][0]*_rho2+C[2][0]*_rho2*_rho2+C[0][1]*_Delta+C[0][2]*_Delta*_Delta;
}

void BToDtaunu::SettBmH(double tBmH){
  _gSR = -mb_quark*pow(tBmH,2);
} 

void BToDtaunu::SetFF(double rho2, double Delta, double gSR){
  _rho2  = rho2;
  _Delta = Delta;
  _gSR   = gSR;
}

BToDtaunu::~BToDtaunu() {
}

