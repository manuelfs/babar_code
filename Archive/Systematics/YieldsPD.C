// Calculates yields after fit for m2>1 GeV^2
#include "TH2F.h"
#include "TString.h"
#include "TFile.h"
#include <iomanip>
#include <fstream>
#include <iostream>

using std::cout;
using std::endl;

void YieldsPD(int IsoConst = 0){

  int IndPdf[2][9] = {{1,9,5,13,17,41,51,61,37},{21,25,29,45,55,65,33,-1,-1}};
  double Yields[70], TruYields[70], TailYields[70];
  TString buffer, textName = "FitAll/fits/TextFinalDataIter2.txt";
  fstream textFile; textFile.open(textName,fstream::in);
  for(int i=0; i<8; i++){
    for(int j=0; j<9; j++){
      textFile>>buffer;
      buffer.ReplaceAll("Yield[",""); buffer.ReplaceAll("]",""); 
      int index = buffer.Atoi();
      textFile>>TruYields[index]>>Yields[index]>>buffer>>buffer>>buffer;
      //cout<<"True: "<<TruYields[index]<<"\tFit: "<<Yields[index]<<"\tIndex "<<index<<endl;
      if(i%2==1 && j==6){
	textFile>>buffer;
	break;
      }
    }
    //cout<<endl;
  }
  double ChanYield[] = {0,0,0,0};
  double m2cut[2];
  TString m2CutName = "babar_code/Systematics/MmissCut.txt";
  fstream m2File; m2File.open(m2CutName,fstream::in);
  m2File>>m2cut[0]>>m2cut[1];
  TString outName = "FitAll/PD/Initial_RAll.txt";
  fstream outFile; outFile.open(outName,fstream::out);
  TString constName = "FitAll/PD/Constraints_RAll.txt";
  fstream constFile; constFile.open(constName,fstream::out);
  double m2min = -4, m2max = 12; 
  TFile hfile("keys/root/Fit/pdfKeys_1_Fit.root"); 
  TH2F*  h = (TH2F *)(hfile.Get("h2"))->Clone("ForBins");
  h->SetDirectory(0);
  int nM2bin = h->GetNbinsX(), nPlbin = h->GetNbinsY();
  int binCut0 = (int)((m2cut[0]-m2min)/(m2max-m2min)*(double)nM2bin)+1;
  int binCut1 = (int)((m2cut[1]-m2min)/(m2max-m2min)*(double)nM2bin);
  if(binCut1>nM2bin) binCut1 = nM2bin;
  //cout<<binCut0<<"\t"<<binCut1<<endl;
  for(int i=0; i<8; i++){
    for(int chan=0; chan<4; chan++){
      if(IndPdf[chan>3][i]<0) continue;
      int index = IndPdf[chan>3][i]+chan-4*(chan>3);
      TString hname = "keys/root/Fit/pdfKeys_"; hname += index; hname += "_Fit.root";
      TFile hfile(hname); 
      TString h2Name = "h2"; 
      TString pdfName = "pdf"; pdfName += index;
      TH2F*  h2 = (TH2F *)(hfile.Get(h2Name))->Clone(pdfName);
      h2->SetDirectory(0);
      double totInt = h2->Integral();
      double cutInt = h2->Integral(binCut0,binCut1,1,nPlbin);
      TailYields[index] = Yields[index]*cutInt/totInt;
      if(index>=5&&index<=8) 
	constFile<<index<<"\t"<<TailYields[index]/TailYields[index-5+2*(index%2)]<<"  +- 0.01"<<endl;
      //else if(index>=61) 
      //constFile<<index<<"\t"<<TailYields[index]/TailYields[index-10]<<"  +- 0.01"<<endl;
      else if(IsoConst && index>=3&&index<=4 ) 
	constFile<<index<<"\t"<<TailYields[index]/TailYields[index-2]<<"  +- 0.01"<<endl;
      else outFile<<index<<"\t"<<TailYields[index]<<endl;
      ChanYield[chan] += TailYields[index];
      h2->Delete();
    }
  }
  for(int chan=0; chan<4; chan++)cout<<chan<<" has "<<ChanYield[chan]<<endl;
  cout<<"Written "<<outName<<" and "<<constName<<endl;
}



