// Calculates yields after fit for m2>1 GeV^2
#include "TH2F.h"
#include "TString.h"
#include "TFile.h"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <iomanip>
#include <fstream>
#include <iostream>

using std::cout;
using std::endl;

void ConstraintErrors(){

  double Yields[70], TruYields[70], Errors[70];
  TString buffer, textName = "FitAll/fits/TextFinalDataNo2x100.txt", sError;
  fstream textFile; textFile.open(textName,fstream::in);
  for(int i=0; i<8; i++){
    for(int j=0; j<9; j++){
      textFile>>buffer;
      buffer.ReplaceAll("Yield[",""); buffer.ReplaceAll("]",""); 
      int index = buffer.Atoi();
      textFile>>TruYields[index]>>Yields[index]>>buffer>>sError>>buffer;
      Errors[index] = sError.Atof();
      //cout<<"True: "<<TruYields[index]<<"\tFit: "<<Yields[index]<<"\tIndex "<<index<<endl;
      if(i%2==1 && j==6){
	textFile>>buffer;
	break;
      }
    }
    //cout<<endl;
  }
  TString outName = "babar_code/Systematics/Text/ConstraintsError.txt";
  fstream outFile; outFile.open(outName,fstream::out);
  double errorD[2], errorDs[2];
  for(int cand=0; cand<2; cand++){
    errorD[cand] = (Yields[13+2*cand]/Yields[10+2*cand])/(TruYields[13+2*cand]/TruYields[10+2*cand])-1;
    errorDs[cand] = sqrt(pow(Errors[13+2*cand]/Yields[13+2*cand],2)+
			 pow(Errors[10+2*cand]/Yields[10+2*cand],2));
    outFile<<RoundNumber(fabs(errorD[cand]),4)<<endl<<RoundNumber(errorDs[cand],4)<<endl;
  }
  cout<<"Written "<<outName<<endl;
}



