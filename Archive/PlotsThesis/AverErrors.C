//------------------------------------------------------------------------
// Description:
//    Average of the R(D) and R(D*) by Belle and BaBar
//
//    CalcRatio(double BTau[3], double Bl[2], double RD[3])
//       Ratio of D(*)taunu and D(*)lnu BF
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      12/04/11 manuelf -- Created 
//------------------------------------------------------------------------

#include "TString.h"
#include "TMatrixT.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

TString RoundNumber(double n, int e, double d=1);
void CalcRatio(double BTau[3], double Bl[2], double RD[3]);
void averageN(double RD[][3], int nRD, double rho, double averRD[]);
void PrintRD(double RD[3], int isDs, int digits=3);
double totError(double E1, double E2);

void AverErrors(double  BelleRho = 0.75, double BaBeRho = 0.25){
  TString RDNames[] = {"R(D) ", "R(D*)"};
  double B0Dslnu[] = {4.63, 0.27}, BmDslnu[2], BmDlnu[] = {2.23, 0.11}; // Dslnu from Belle, Dlnu from PDG2010
  double RD[2][4][3];
  TString MeasNames[] = {"Belle 2007", "Belle 2009", "Belle 2010", "BaBar 2012"};
  TString SepLine = "==========================================================================";
  for(int iE=0; iE<2; iE++) BmDslnu[iE] = B0Dslnu[iE]*1.638e-12/1.525e-12;

  // Belle 2007
  double BF_B0DsmTauNu_07[] = {2.02, 0.39, 0.37}; // BF(B0->D*-TauNu)
  CalcRatio(BF_B0DsmTauNu_07, B0Dslnu, RD[1][0]);

  // Belle 2009
  double RD_09[2][3]  = {{0.70, 0.19, 0.10}, {0.48, 0.21, 0.06}};  // R(D0)  and R(D+)
  double RDs_09[2][3] = {{0.47, 0.11, 0.07}, {0.48, 0.13, 0.05}};  // R(D*0) and R(D*+)
  averageN(RD_09,  2, BelleRho, RD[0][1]);
  averageN(RDs_09, 2, BelleRho, RD[1][1]);

  // Belle 2010
  double BF_BmD0TauNu_10[]  = {0.77, 0.22, 0.12}; // BF(B-->D0TauNu)
  double BF_BmDs0TauNu_10[] = {2.12, 0.28, 0.29}; // BF(B-->D*0TauNu)
  CalcRatio(BF_BmD0TauNu_10,  BmDlnu,  RD[0][2]);
  CalcRatio(BF_BmDs0TauNu_10, BmDslnu, RD[1][2]);

  // BaBar 2012
  RD[0][3][0] = 0.4402; RD[0][3][1] = 0.0577; RD[0][3][2] = 0.0423; 
  RD[1][3][0] = 0.3316; RD[1][3][1] = 0.0236; RD[1][3][2] = 0.0183; 
  // BaBar 2008
//   RD[0][3][0] = 0.416;  RD[0][3][1] = 0.117;  RD[0][3][2] = 0.052; 
//   RD[1][3][0] = 0.297;  RD[1][3][1] = 0.056;  RD[1][3][2] = 0.018; 

  for(int iMeas=0; iMeas<4; iMeas++){
    cout<<endl<<MeasNames[iMeas]<<endl<<SepLine<<endl;
    for(int isDs=0; isDs<2; isDs++){
      if(iMeas==0 && isDs==0) continue;
      PrintRD(RD[isDs][iMeas], isDs);
    }
  }

  // Belle averages
  double AverBelle[2][2][3];
  averageN(RD[0]+1, 2, BelleRho, AverBelle[0][0]);
  averageN(RD[1], 3, BelleRho, AverBelle[1][0]);
  cout<<endl<<"Belle averages"<<endl<<SepLine<<endl;
  for(int isDs=0; isDs<2; isDs++) PrintRD(AverBelle[isDs][0], isDs);
  for(int isDs=0; isDs<2; isDs++) {
    double BeSig = sqrt(pow(AverBelle[isDs][0][1],2)+pow(AverBelle[isDs][0][2],2));
    double BaSig = sqrt(pow(RD[isDs][3][1],2)+pow(RD[isDs][3][2],2));
    double RDsig = sqrt(pow(BaSig,2)+pow(BeSig,2)-BaBeRho*BaSig*BeSig);
    cout<<"Belle's "<< RDNames[isDs]<<" exceeds BaBar's by "<<RoundNumber(AverBelle[isDs][0][0]-RD[isDs][3][0],2,RDsig)<<" sigma"<<endl;
  }

  // BaBe averages
  double AverBaBe[2][3];
  cout<<endl<<"BaBe averages"<<endl<<SepLine<<endl;
  for(int isDs=0; isDs<2; isDs++) {
    for(int iRD=0; iRD<3; iRD++) AverBelle[isDs][1][iRD] = RD[isDs][3][iRD];
    averageN(AverBelle[isDs], 2, BaBeRho, AverBaBe[isDs]);
    PrintRD(AverBaBe[isDs], isDs, 4);
  }

  // The correlation average is probably wrong
  double sigBa[2] = {totError(RD[0][3][1], RD[0][3][2]), totError(RD[1][3][1], RD[1][3][2])};
  double sigBe[2] = {totError(AverBelle[0][0][1], AverBelle[0][0][2]), totError(AverBelle[1][0][1], AverBelle[1][0][2])};
  double sigTo[2] = {totError(sigBa[0], sigBe[0]), totError(sigBa[1], sigBe[1])};
  double rhoTot = -0.11*sigBa[0]*sigBa[1]/sigTo[0]/sigTo[1];
  cout<<"Correlation is "<<RoundNumber(rhoTot,2)<<endl<<endl;
}

double totError(double E1, double E2){
  return sqrt(pow(E1,2)+pow(E2,2));
}

void PrintRD(double RD[3], int isDs, int digits){
  TString RDNames[] = {"R(D)  = ", "R(D*) = "};
  double totalE = sqrt(pow(RD[1],2)+pow(RD[2],2));
  cout<<RDNames[isDs]<<RoundNumber(RD[0],digits)<<" +- "<<RoundNumber(RD[1],digits)<<" +- "<<
    RoundNumber(RD[2],digits)<<" => Total error is "<<RoundNumber(totalE,digits)<<", a "<<
    RoundNumber(totalE*100,1,RD[0])<<"% error"<<endl;
}

// Average of nRD measurements, all with a correlation of rho. 
// Uses Bob Kowalewski "Methodology and code used in HFAG semileptonic averages" (2005)
void averageN(double RD[][3], int nRD, double rho, double averRD[]){
  TMatrixT<double> StatE(nRD,nRD), SystE(nRD,nRD);
  
  for(int iR=0; iR<nRD; iR++){
    for(int iRD=0; iRD<3; iRD++) averRD[iRD] = 0;
    for(int jR=0; jR<nRD; jR++){
      if(iR==jR){
	StatE(iR, jR) = pow(RD[iR][1],2);
	SystE(iR, jR) = pow(RD[iR][2],2);
      } else {
	StatE(iR, jR) = 0;
	SystE(iR, jR) = rho*RD[iR][2]*RD[(iR+1)%nRD][2];
      }
    }
  }

  TMatrixT<double> TotE(StatE); TotE += SystE; TotE.Invert();
  TMatrixT<double> WStatW(nRD,nRD), WSystW(nRD,nRD);
  WStatW.Mult(TotE, StatE); WStatW *= TotE;
  WSystW.Mult(TotE, SystE); WSystW *= TotE;
  double sumW = 0;
  for(int iR=0; iR<nRD; iR++){
    for(int jR=0; jR<nRD; jR++){
      averRD[0] += TotE(iR,jR)*RD[iR][0];
      averRD[1] += WStatW(iR,jR);
      averRD[2] += WSystW(iR,jR);
      sumW += TotE(iR,jR);
    }
  }
  for(int iRD=0; iRD<3; iRD++){
    averRD[iRD] /= sumW;
    if(iRD>0) {
      averRD[iRD] /= sumW;
      averRD[iRD] = sqrt(averRD[iRD]);
    }
  }
}

void CalcRatio(double BTau[3], double Bl[2], double RD[3]){
  RD[0] = BTau[0]/Bl[0];
  RD[1] = BTau[1]/Bl[0];
  RD[2] = RD[0]*sqrt(pow(BTau[2]/BTau[0],2)+pow(Bl[1]/Bl[0],2));
}

TString RoundNumber(double n, int e, double d){
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
