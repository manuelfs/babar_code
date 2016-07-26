#include "TString.h"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void Results_Table(TString textName){

  int isIso = 0;
  if(textName.Contains("Iso")) isIso = 1;
  double Yields[16][9], Errors[8][9];
  TString buffer, sError;
  fstream textFile; textFile.open(textName,fstream::in);
  int IndText[8] = {0,4,1,5,2,6,3,7};
  for(int i=0; i<8; i++){
    int index = IndText[i];
    for(int j=0; j<9; j++){
      textFile>>buffer>>Yields[index+8][j]>>Yields[index][j]>>buffer>>sError>>buffer;
      Errors[index][j] = sError.Atof();
      //cout<<"True: "<<Yields[index+8][j]<<"\tFit: "<<Yields[index][j]<<" +- "<< Errors[index][j]
      //<<"\tChannel "<<index<<endl;
      if(index>3 && j==6){
	textFile>>buffer;
	break;
      }
    }
    //cout<<endl;
  }

  double totY[2] = {0,0}, totE[2] = {0,0};
  for(int cand=0; cand<2; cand++){
//     double fMC = Yields[8+2*cand][3]/Yields[8+1+2*cand][1];
//     double fFit = Yields[2*cand][3]/Yields[1+2*cand][1];
//     double eMC = fMC*sqrt(1/Yields[8+2*cand][3]+1/Yields[8+1+2*cand][1])/3.;
//     double eFit = fFit*sqrt(pow(Errors[2*cand][3]/Yields[2*cand][3],2)
// 			    +pow(Errors[1+2*cand][1]/Yields[1+2*cand][1],2));

    double fMC = Yields[8+2*cand][2]/Yields[8+1+2*cand][0];
    double fFit = Yields[2*cand][2]/Yields[1+2*cand][0];
    double eMC = fMC*sqrt(1/Yields[8+2*cand][2]+1/Yields[8+1+2*cand][0])/3.;
    double eFit = fFit*sqrt(pow(Errors[2*cand][2]/Yields[2*cand][0],2)
			    +pow(Errors[1+2*cand][0]/Yields[1+2*cand][0],2));


    double eRat = fFit/fMC*sqrt(pow(eMC/fMC,2)+pow(eFit/fFit,2));
    cout<<"MC: "<<RoundNumber(fMC,3)<<"+- "<<RoundNumber(eMC,3)
	<<"\t Fit: "<<RoundNumber(fFit,3)<<"+- "<<RoundNumber(eFit,3)
	<<"\t Fit/MC: "<<RoundNumber(fFit/fMC,3)<<"+- "<<RoundNumber(eRat,3)<<endl;
  }
}
