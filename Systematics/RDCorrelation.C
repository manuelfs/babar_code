#include "babar_code/PlotsThesis/PlotUtils.cc"
#include "TString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMatrixT.h"
#include "TRandom3.h"
#include "TMath.h"
#include <fstream>
#include <iostream>

#define nRows 5
using namespace std;
using std::cout;
using std::endl;

void RDCorrelation(int nRep = 10000, int isDcIso=2){

  TString textName = "FitAll/fits/TextFinalIsoDataNe4x100.txt";
  TString correlName = "FitAll/fits/CorrelationIsoDataNe4x100.txt";
  if(isDcIso<2){
    textName.ReplaceAll("Iso","");
    correlName.ReplaceAll("Iso","");
  }
  double Yield[2][70], Error[70], Correl[22][11];
  ReadFitFile(textName, Yield, Error);
  ReadCorrelFile(correlName, Correl);

  int NInd[]   = {1, 9, 2,10,13}; // {Ntau, Nl, N*tau, N*l, N*lx}
  int RhoInd[] = {0, 2, 3, 4, 1}, dRho = 0;
  if(isDcIso==1) {
    dRho = 11;
    for(int row=0; row<5; row++) NInd[row] += 2;
  }

  TMatrixT<double> CovM(nRows, nRows); //CovM is the covariance matrix
  for(int row=0; row<nRows; row++)
    for(int col=0; col<nRows; col++)
      CovM(row,col) = Correl[dRho+RhoInd[row]][RhoInd[col]]*
	Error[NInd[row]]*Error[NInd[col]];
  CovM = Choleski(CovM, nRows);
  
  double RD[3][2], NYield[nRows], NZ[nRows], dRD[2];
  for(int row=0; row<3; row++)
    for(int col=0; col<2; col++) RD[row][col] = 0;

  TRandom3 rand(0); 
  for(int rep=0; rep<nRep; rep++){
    for(int row=0; row<nRows; row++) NZ[row] = rand.Gaus(0, 1);
    for(int row=0; row<nRows; row++){
      NYield[row] = Yield[0][NInd[row]];
      for(int col=0; col<nRows; col++) NYield[row] += NZ[col]*CovM(row,col);
    }
    dRD[0] = NYield[0]/NYield[1];
    dRD[1] = NYield[2]/(NYield[3]+NYield[4]);
    for(int row=0; row<2; row++){
      RD[0][row] += dRD[row];
      RD[1][row] += dRD[row]*dRD[row];
    }
    RD[2][0] += dRD[0]*dRD[1];
  }
  double N = (double)nRep;
  RD[2][0] = (N*RD[2][0] - RD[0][0]*RD[0][1])/
    sqrt(N*RD[1][0] - pow(RD[0][0],2))/
    sqrt(N*RD[1][1] - pow(RD[0][1],2));
  for(int row=0; row<2; row++){
    RD[0][row] /= N; 
    RD[1][row] = sqrt((RD[1][row]-RD[0][row]*RD[0][row]*N)/(N-1)); 
  }
  cout<<endl<<"RD is "<<RD[0][0]<<" +- "<<RD[1][0]<<"  and RD* is "<<
    RD[0][1]<<" +- "<<RD[1][1]<<", with a correlation of "<< RD[2][0]<<endl<<endl;  
}


