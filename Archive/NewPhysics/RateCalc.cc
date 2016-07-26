#include "babar_code/NewPhysics/RateCalc.hh"
#include "TRandom3.h"


RateCalc::RateCalc(TString NameFile) { 

  SetMasses(1); // Charged B
  double r = mDs/mB;
  R[3]   = 0.97;
  R[0]   = (2*r+pow(1-r,2)*R[3]-(1-r)*R[2])/r/(1+r);

  ReadFF(NameFile);

  Vcb    = 0.0375;
  F1     = 0.921;
  VcbG1  = 0.03912*1.0816;
  rhoD2  = 1.186;

  mc   = 0.901;
  mb   = 4.20;
  gsL  = 0;
  gsR  = 0;
}

void RateCalc::SetMasses(int isBm){
  if(isBm) {
    mB  = 5.2792; mDs = 2.00696; mD  = 1.8648; BLifeTime = 1.638e-12; // Charged B
  } else {
    mB  = 5.2795; mDs = 2.01025; mD  = 1.8696; BLifeTime = 1.525e-12; // Neutral B
  }
  RDs = 2.*sqrt(mB*mDs)/(mB+mDs);
  mB2 = mB*mB; mDs2 = mDs*mDs; mD2 = mD*mD;
  Dsmaxq2 = pow(mB-mDs,2); Dmaxq2 = pow(mB-mD,2);
  _isBm = isBm;
}

// ============================  R(D*) ========================================

  // Kinematic variables
double RateCalc::wDs(double q2) {
  return (mB2+mDs2-q2)/(2*mB*mDs);
}

double RateCalc::pDs(double q2) {
  return sqrt(mB2*mB2+mDs2*mDs2+q2*q2-2*(mB2*mDs2+mDs2*q2+q2*mB2))/(2*mB);
}
 
  // D*TauNu Form Factors
double RateCalc::hA1(double w) {
  double z = (sqrt(w+1)-sqrt(2))/(sqrt(w+1)+sqrt(2));
  return F1*(1 - 8.*rhoDs2*z + (53.*rhoDs2-15.)*z*z - (231.*rhoDs2-91.)*z*z*z);
}

double RateCalc::R1(double w) {
  return R[1] - 0.12*(w-1) + 0.05*(w-1)*(w-1);
}

double RateCalc::R2(double w) {
  return R[2] + 0.11*(w-1) - 0.06*(w-1)*(w-1);
}

double RateCalc::R0(double w) {
  return R[0] - 0.11*(w-1) + 0.01*(w-1)*(w-1);
}

double RateCalc::A1(double q2) {
  double w = wDs(q2);
  return (w+1)*RDs*hA1(w)/2.;
}

double RateCalc::A0(double q2) {
  double w = wDs(q2);
  return R0(w)/RDs*hA1(w);
}

double RateCalc::V(double q2) {
  double w = wDs(q2);
  return R1(w)/RDs*hA1(w);
}

double RateCalc::A2(double q2) {
  double w = wDs(q2);
  return R2(w)/RDs*hA1(w);
}

double RateCalc::Hpp(double q2) {
  return (mB+mDs)*A1(q2) - 2*mB/(mB+mDs)*pDs(q2)*V(q2);
}

double RateCalc::Hmm(double q2) {
  return (mB+mDs)*A1(q2) + 2*mB/(mB+mDs)*pDs(q2)*V(q2);
}

double RateCalc::H00(double q2) {
  double Term1 = (mB2-mDs2-q2)*(mB+mDs)*A1(q2);
  double Term2 = 4*mB2*pow(pDs(q2),2)/(mB+mDs)*A2(q2);
  return (Term1-Term2)/(2*mDs*sqrt(q2));
}

double RateCalc::H0t(double q2) {
  return 2*mB*pDs(q2)/sqrt(q2)*A0(q2)*(1-(gsL-gsR)*q2/(mb+mc));
}

  // Rate
double RateCalc::GammaDs_q2(double q2, double ml) {
  double ml2 = ml*ml;
  if(q2<ml2 || q2>=Dsmaxq2){
    //cout<<q2<<" outside q2 range ("<<ml2<<", "<<Dsmaxq2<<")"<<endl; 
    return 0;
  }

  double Factor = GF*GF*Vcb*Vcb*pDs(q2)*q2*pow(1-ml2/q2,2)/(96*pi*pi*pi*mB2);
  double Term1 = (pow(Hpp(q2),2)+pow(Hmm(q2),2)+pow(H00(q2),2))*(1+ml2/2/q2);
  double Term2 = 3*ml2*pow(H0t(q2),2)/(2*q2);
  return Factor*(Term1+Term2);
}

// ============================  R(D) =========================================
  // DTauNu Form Factor
double RateCalc::G(double w){
  double z = (sqrt(w+1)-sqrt(2))/(sqrt(w+1)+sqrt(2));
  return 1 - 8*rhoD2*z + (51*rhoD2-10)*z*z - (252*rhoD2-84)*z*z*z;
}

  // Rate Tanaka, Watanabe 2010
double RateCalc::GammaD_q2_Tanaka(double q2, double ml) {
  double ml2 = ml*ml;
  if(q2<ml2 || q2>=Dmaxq2) return 0;
  
  double w = (mB*mB+mD*mD-q2)/(2.*mB*mD);
  double wm1 = w - 1.0;
  double S1 = G(w) * (1 +Delta*(-0.019 + 0.041*wm1 - 0.015*wm1*wm1));

  double RD = (2*sqrt(mB*mD))/(mB+mD);
  double F1 = G(w)/RD;
  double F0 = (1-q2/pow(mB+mD,2))/RD*S1;

  double fplus = F1;
  double fminus = (F0 - F1) * (mB*mB-mD*mD) / q2;

  double pD = mB*mB*mB*mB + mD*mD*mD*mD + q2*q2 -
    2.*mB*mB*mD*mD - 2.*mB*mB*q2 - 2.*mD*mD*q2;
  if (pD < 0) pD = 0;
  else pD = sqrt(pD)/(2.*mB);

  double h0  = 2.*mB*pD*fplus / sqrt(q2);
  double ht  = ((mB*mB-mD*mD)*fplus + q2*fminus) / sqrt(q2);
  double hl  = h0*h0;
  double hs  = 3.*ht*ht;
  double LH = 2./3*(q2-ml2) * (hl+(ml2/(2.*q2))*(hl+hs));
  return GF*GF/pow(2*pi,3)*VcbG1*VcbG1*(q2-ml2)*pD/(8*mB*mB*q2)*LH;

}

double RateCalc::GammaD_q2(double q2, double ml) {
  double ml2 = ml*ml;
  if(q2<ml2 || q2>=Dmaxq2) return 0;

  double dw_dq2 = 1/(2*mB*mD);
  double w = (mB2+mD2-q2)/(2*mB*mD);
  double t = mB2+mD2-2*w*mD*mB;
  double rhoV = 4*pow(1+mD/mB,2)*pow(mD/mB,3)*pow(w*w-1,3/2.)*pow(1-ml2/t,2)*
    (1+ml2/(2*t))*pow(G(w),2);
  double rhoS = 3*mB2/(2*t)/(1+ml2/(2*t))*(1+w)/(1-w)*pow(Delta,2);
  double rhoS_NP = pow(1+t*(gsL+gsR)/(mb-mc),2)*rhoS;
  double Factor = dw_dq2*GF*GF*VcbG1*VcbG1*pow(mB,5)/(192*pi*pi*pi);
  return Factor*rhoV*(1-ml2/mB2*rhoS_NP);
}


// =========================  Auxiliary functions  =============================

double RateCalc::Gamma_q2(int isDs, double q2, double ml) {
  if(0){
    //Delta = 1; Delta_Err = 1;
    return (isDs ? GammaDs_q2(q2,ml) : GammaD_q2_Tanaka(q2,ml));
  }
  return (isDs ? GammaDs_q2(q2,ml) : GammaD_q2(q2,ml));
}

// Simpson integration
double RateCalc::IntRate(double minX, double maxX, double ml, int isDs, int nPoints){
  double intF = 0, x = minX, dx = (maxX-minX)/(double)nPoints;
  double Fmin = Gamma_q2(isDs,x,ml);
  for(int step=0; step<nPoints; step++){
    intF += dx/6.*(Fmin+4*Gamma_q2(isDs,x+dx/2,ml));
    x += dx;
    Fmin = Gamma_q2(isDs,x,ml);
    intF +=dx/6.*Fmin;
  }
  
  return intF;
}

double RateCalc::Rate(int isDs, double ml){
  double maxq2 = (isDs ? Dsmaxq2 : Dmaxq2);
  double totalRate = hbar/BLifeTime;
  return IntRate(ml*ml, maxq2, ml, isDs)/totalRate;
}

double RateCalc::RRate(int isDs, double ml){
  return Rate(isDs, mTau)/Rate(isDs, ml);
}

void RateCalc::Print(){
  double Coef[3];

  Polynomial(Coef,1);
  cout<<"R(D*) = "<<RoundNumber(Coef[0],4)<<"*(1+"<<RoundNumber(Coef[1],4)<<"*mTau(gSR-gSL)+"
      <<RoundNumber(Coef[2],4)<<"*(mTau(gSR-gSL))^2)"<<endl;

  Polynomial(Coef,0);
  cout<<"R(D)  = "<<RoundNumber(Coef[0],4)<<"*(1+"<<RoundNumber(Coef[1],4)<<"*mTau(gSR+gSL)+"
      <<RoundNumber(Coef[2],4)<<"*(mTau(gSR+gSL))^2)"<<endl;
}

void RateCalc::Polynomial(double Coef[3], int isDs, double ml){
  double F[3], gsR_Fix = gsR, gsL_Fix = gsL;

  gsL = 0; 
  gsR = -1; F[0] = RRate(isDs, ml);
  gsR = 0;  F[1] = RRate(isDs, ml);
  gsR = 1;  F[2] = RRate(isDs, ml);

  gsR = gsR_Fix; gsL = gsL_Fix;

//   // 1203.2654 style
//   Coef[0] = F[1];
//   Coef[1] = (F[2]-F[0])/(2*F[1]*mTau);
//   Coef[2] = ((F[2]+F[0])/(2*F[1])-1.)/(mTau*mTau);

  Coef[0] = F[1];
  Coef[1] = (F[2]-F[0])/2;
  Coef[2] = -F[1]+(F[2]+F[0])/2;
}

void RateCalc::AverPoly(double Coef[3], int isDs){
  double ACoef[3];
  int fix_isBm = _isBm;
  for(int iC=0; iC<3; iC++) Coef[iC] = 0;

  SetMasses(1);
  Polynomial(ACoef, isDs, mE);
  for(int iC=0; iC<3; iC++) Coef[iC] += ACoef[iC];
  Polynomial(ACoef, isDs, mMu);
  for(int iC=0; iC<3; iC++) Coef[iC] += ACoef[iC];

  SetMasses(0);
  Polynomial(ACoef, isDs, mE);
  for(int iC=0; iC<3; iC++) Coef[iC] += ACoef[iC];
  Polynomial(ACoef, isDs, mMu);
  for(int iC=0; iC<3; iC++) Coef[iC] += ACoef[iC];

  for(int iC=0; iC<3; iC++) Coef[iC] /= 4.;
  SetMasses(fix_isBm);
}

TString RateCalc::RoundNumber(double n, int e, double d){
  if(d==0) return " - ";
  double neg = 1; if(n*d<0) neg = -1;
  double b = (int)(neg*n/d*pow(10.,(double)e)+0.5);
  b /= pow(10.,(double)e)*neg;
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(!result.Contains(".") && e != 0) result += ".";
  
  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<e-afterdot.Length(); i++)
    result += "0";
  return result;
}

void RateCalc::ReadFF(TString NameFile){
  fstream textFile; textFile.open(NameFile,fstream::in);
  TString dummy, Name;
  double dum, Value, Error;
  TString VarNames[]  = {"rhoD2", "Delta", "rhoDs2", "R1", "R2", 
			 "Corrrho2R1", "CorrR1R2", "CorrR2rho2", "R0", "mb", "mc"};
  double *VarValues[] ={&rhoD2, &Delta, &rhoDs2, &R[1], &R[2],
			 &Correlation[0], &Correlation[1], &Correlation[2],&R[0], &mb, &mc};
  double *VarErrors[] ={&rhoD2_Err, &Delta_Err, &rhoDs2_Err, &R_Err[1], &R_Err[2],
			&dum, &dum, &dum, &R_Err[0], &mb_Err, &mc_Err};

  while(textFile){
    textFile>>Name>>dummy>>Value>>dummy>>Error;
    for(int iVar=0; iVar<nVar; iVar++){
      if(Name == VarNames[iVar]) {
	*VarValues[iVar] = Value;
	*VarErrors[iVar] = Error;
      }
    }
  }
}

void RateCalc::Errors(double Coef[3][3], int isDs, int nRep){
  double fixFF[nVar], fZ[3], dum, C[3];//, mean[nVar], rms[nVar], corr=0;
//   TString VarNames[]  = {"rhoD2", "Delta", "rhoDs2", "R1", "R2", 
// 			 "Corrrho2R1", "CorrR1R2", "CorrR2rho2", "R0", "mb", "mc"};
  double *VarValues[] ={&rhoD2, &Delta, &rhoDs2, &R[1], &R[2],
			 &Correlation[0], &Correlation[1], &Correlation[2],&R[0], &mb, &mc};
  double *VarErrors[] ={&rhoD2_Err, &Delta_Err, &rhoDs2_Err, &R_Err[1], &R_Err[2],
			&dum, &dum, &dum, &R_Err[0], &mb_Err, &mc_Err};
  double valFF[3][3] = {{rhoDs2, R[1], R[2]},{rhoDs2_Err, R_Err[1], R_Err[2]},
			{Correlation[0],Correlation[1],Correlation[2]}};
  for(int iC=0; iC<3; iC++)
    for(int iT=0; iT<3; iT++)
      Coef[iT][iC] = 0; 
  for(int iVar=0; iVar<nVar; iVar++) {
    fixFF[iVar] = *VarValues[iVar];
    //mean[iVar] = 0; rms[iVar] = 0;
  }

  TRandom3 rand(0);
  TMatrixT<double> CovFF(3,3); //CovFF is the covariance matrix
  for(int iFF=0; iFF<3; iFF++) {
    CovFF(iFF,iFF) = pow(valFF[1][iFF],2);
    CovFF(iFF,(iFF+1)%3) = valFF[2][iFF]*valFF[1][iFF]*valFF[1][(iFF+1)%3];
    CovFF((iFF+1)%3,iFF) = CovFF(iFF,(iFF+1)%3);
  }
  CovFF = Choleski(CovFF,3);

  for(int rep=0; rep<nRep; rep++){
    for(int chan=0; chan<3; chan++) fZ[chan] = rand.Gaus(0, 1);
    for(int chan=0; chan<3; chan++){
      *VarValues[irhoDs2+chan] = fixFF[irhoDs2+chan];
      for(int col=0; col<3; col++) *VarValues[irhoDs2+chan] += fZ[col]*CovFF(chan,col);
    }
    for(int iVar=0; iVar<nVar; iVar++){
      if(iVar<irhoDs2 || iVar>iCorrelation2) 
	*VarValues[iVar] = rand.Gaus(fixFF[iVar], *VarErrors[iVar]);
//       mean[iVar] += *VarValues[iVar];
//       rms[iVar] += (*VarValues[iVar])*(*VarValues[iVar]);
//       if(iVar==iR2) corr += (*VarValues[iVar])*(*VarValues[iVar-1]);
    }
    AverPoly(C,isDs);
    C[1] *= -mb; C[2] *= mb*mb; // Translation from g_SR to (tanBeta/mH)^2
    for(int iC=0; iC<3; iC++){
      Coef[0][iC] += C[iC];
      Coef[1][iC] += C[iC]*C[iC];
      Coef[2][iC] += C[iC]*C[(iC+1)%3];
    }
  }
  for(int iVar=0; iVar<nVar; iVar++) *VarValues[iVar] = fixFF[iVar];
  AverPoly(C,isDs);
  C[1] *= -mb; C[2] *= mb*mb; // Translation from g_SR to (tanBeta/mH)^2
  double N = nRep;
  for(int iC=0; iC<3; iC++) 
    Coef[2][iC] = (N*Coef[2][iC] - Coef[0][iC]*Coef[0][(iC+1)%3])/
      sqrt(N*Coef[1][iC]        - pow(Coef[0][iC],2))/
      sqrt(N*Coef[1][(iC+1)%3] - pow(Coef[0][(iC+1)%3],2));
  for(int iC=0; iC<3; iC++){
    Coef[0][iC] /= N; 
    Coef[1][iC] =  sqrt((Coef[1][iC]-Coef[0][iC]*Coef[0][iC]*N)/(N-1)); 
    Coef[0][iC] = C[iC]; 
  }
  cout<<endl<<"  double Coef"<<isDs<<"[3][3] = {";
  for(int nRho=0; nRho<3; nRho++){
    cout<<"{";
    for(int nDelta=0; nDelta<3; nDelta++) {
      cout<<RoundNumber(Coef[nRho][nDelta],4);
      if(nDelta<2)cout<<", ";
    }
    cout<<"}";
    if(nRho<2)cout<<", ";
  }
  cout<<"};"<<endl<<endl;
//   cout<<"Correlation = "<<(corr*N-mean[iR1]*mean[iR2])/
//     sqrt(N*rms[iR1]-pow(mean[iR1],2))/sqrt(N*rms[iR2]-pow(mean[iR2],2))<<endl;
//   for(int iVar=0; iVar<nVar; iVar++) {
//     mean[iVar] /= N; 
//     rms[iVar] = sqrt((rms[iVar]-mean[iVar]*mean[iVar]*N)/(N-1));
//     cout<<VarNames[iVar]<<"\t = "<<mean[iVar]<<" +- "<<rms[iVar]<<endl;
//   }
}


// Cholesky-Banachiewicz algorithm
TMatrixT<double> RateCalc::Choleski(TMatrixT<double> A, int nRows){
  TMatrixT<double> L(nRows,nRows);
  
  for(int row=0; row<nRows; row++)
    for(int col=0; col<nRows; col++)
      L(row,col) = 0;
  
  for(int row=0; row<nRows; row++){
    for(int col=0; col<(row+1); col++){
      if(col==row) {
	for(int irow=0; irow<col; irow++) L(row,col) += L(col,irow)*L(col,irow);
	L(row,col) = sqrt(A(row,col) - L(row,col));
      } else {
	for(int irow=0; irow<col; irow++) L(row,col) += L(row,irow)*L(col,irow);
	L(row,col) = (A(row,col) - L(row,col))/L(col,col);
      }
    }
  }
  return L;
}

RateCalc::~RateCalc() {
}

