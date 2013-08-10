#include "TString.h"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void CountYields(TString textName, int comp, int isDss=0){

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

  double truY[2] = {0,0}, totY[2] = {0,0}, totE[2] = {0,0};
  int isOther = -1;
  if(isDss==0) {
    if(comp==4) isOther = 0;
    if(comp==8) isOther = 6;
  }
  for(int cand=0; cand<4; cand++){
    int iY = 0; if(isIso) iY = cand%2;
    totY[iY] += Yields[cand+4*isDss][comp];
    truY[iY] += Yields[8+cand+4*isDss][comp];
    double err = Errors[cand+4*isDss][comp];
    if(isOther>=0) {
      double YSig = Yields[cand][comp], YDpi = Yields[cand+4][isOther];
      if(isIso) YSig += Yields[cand+2][comp]; 
      err = Errors[cand+4][isOther]*YSig/YDpi;
    }
    totE[iY] += pow(err,2);
  }
  cout<<endl<<"Yields:   True = "<<RoundNumber(truY[isIso],0)
      <<" \t Fitted = "<<RoundNumber(totY[isIso],0)<<" +- "<<RoundNumber(sqrt(totE[isIso]),0)
      <<" \t "<<RoundNumber(totY[isIso]-truY[isIso],1,sqrt(totE[isIso]))<<" sigma away"<<endl;
  cout<<endl;
}
