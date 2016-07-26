#include "babar_code/PlotsThesis/PlotUtils.cc"
#include "TString.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void IntegralPDF(TString hname, double &lowM, double &higM, TString histo="h2");
void WriteHisto(int IndPDF, int IndTimes);

void IntegralKeys(){
  
  double Yields[16][9], lowM, higM, iYields[70][2];
  TString buffer, textName = "FitAll/fits/TextFinalDataSi2x100.txt", hname, Yname, histo;
  fstream textFile; textFile.open(textName,fstream::in);
  int IndText[8] = {0,4,1,5,2,6,3,7};
  for(int i=0; i<8; i++){
    int index = IndText[i];
    for(int j=0; j<9; j++){
      textFile>>Yname>>Yields[index+8][j]>>Yields[index][j]>>buffer>>buffer>>buffer;
      //cout<<"True: "<<Yields[index+8][j]<<"\tFit: "<<Yields[index][j]<<"\tChannel "<<index<<endl;
      Yname.ReplaceAll("Yield[",""); Yname.ReplaceAll("]","");
      iYields[Yname.Atoi()][0] = Yields[index+8][j];
      iYields[Yname.Atoi()][1] = Yields[index][j];
      if(index>3 && j==6){
	textFile>>buffer;
	break;
      }
    }
    //cout<<endl;
  }

  double totDiff=0;
//   double Error[2][3] = {{0.1, 0.045, 0.13}, {0.1, 0.045, 0.1}};
//   int IndPDF[] = {17, 51, 1};
//   for(int cand=0; cand<2; cand++){
//     for(int i=0; i<3; i++){
//       int index = IndPDF[i]+2*cand;
//       hname = "keys/root/fitSigx100/pdfKeys_"; hname += index;
//       hname += "_Fit.root";
//       IntegralPDF(hname, lowM, higM);
//       double diff = iYields[index][1]*(lowM-higM)*Error[cand][i];
//       cout<<index<<" - Low M2: "<<RoundNumber(lowM,4)<<" \t High M2: "<<RoundNumber(higM,4)
// 	  <<" \t Yield diff: "<< RoundNumber(diff,1)<<endl;
//       totDiff += fabs(diff);
//     }
//     cout<<endl;
//   }
//   cout<<"Max total: "<<totDiff<<endl;

  TH1F hist("h","",100,-35,35);
  double meanY[2][3] = {{0,0,0},{0,0,0}}, rmsY[2][3] = {{0,0,0},{0,0,0}};
  int IndPDF[] = {1, 5, 13}, nRep = 140, maxInd[2][3] = {{0,0,0},{0,0,0}};
  //int IndPDF[] = {17, 51, 61}, nRep = 140, maxInd[2][3] = {{0,0,0},{0,0,0}};
  for(int cand=0; cand<2; cand++){
    for(int i=0; i<3; i++){
      double maxDiff=0;
      int index = IndPDF[i]+2*cand;
      hname = "keys/root/Times/hTimes_"; hname += index;
      hname += "_Fit.root";
      for(int rep=0; rep<nRep; rep++){
	histo = "hTimes"; histo += rep;
	IntegralPDF(hname, lowM, higM, histo);
	double diff = iYields[index][1]*(lowM-higM);
	meanY[cand][i] += diff;
	rmsY[cand][i] += diff*diff;
// 	if(fabs(diff)>maxDiff) {
// 	  maxInd[cand][i] = rep; maxDiff = fabs(diff);}
	if(diff>maxDiff) {
	  maxInd[cand][i] = rep; maxDiff = diff;}
	if(index==61) hist.Fill(diff);
      }
      meanY[cand][i] /= (double)nRep;
      rmsY[cand][i] = sqrt((rmsY[cand][i] - meanY[cand][i]*meanY[cand][i]*(double)nRep)
			   /((double)nRep-1));
      maxDiff -= fabs(meanY[cand][i]);
      cout<<index<<" -  Mean: "<<RoundNumber(meanY[cand][i],1)<<" \t RMS: "
	  <<RoundNumber(rmsY[cand][i],1)
	  <<" \t Max diff: "<< maxInd[cand][i]<<endl;//RoundNumber(maxDiff,1)<<endl;
      totDiff += pow(rmsY[cand][i],2);
    }
    cout<<endl;
  }
  cout<<"Total: "<<sqrt(totDiff)<<endl;
//   TCanvas can; can.cd();
//   hist.Draw();
//   can.SaveAs("test.eps");

}

void IntegralPDF(TString hname, double &lowM, double &higM, TString histo){
  TFile hfile(hname); hfile.cd();
  TH2F *h2 = (TH2F *)gDirectory->Get(histo);
  int nPlbin = h2->GetNbinsY(); 

  double totI = h2->Integral(); if(totI<=0) totI = 1;
  lowM = h2->Integral(345,500,1,nPlbin)/totI;
  higM = h2->Integral(564,875,1,nPlbin)/totI;

}

void MovePDFs(){
  TString IndPDF = "63", IndTimes = "6";
  int Indices[12][2] = {{17, 71},{51, 115},{61, 20},{1, 64},{5, 48},{13, 62},
			{19, 79},{53, 0},{63, 6},{3, 117},{7, 55},{15, 76}};
  for(int i=0; i<12; i++) WriteHisto(Indices[i][0], Indices[i][1]);
}


void WriteHisto(int IndPDF, int IndTimes){
  TString TimesName = "keys/root/Times/hTimes_"; TimesName += IndPDF;
  TimesName += "_Fit.root";
  TString hName = "hTimes"; hName += IndTimes;
  TString fName = "keys/root/fitSigx100Test/pdfKeys_"; fName += IndPDF;
  fName += "_Fit.root";

  TFile fTimes(TimesName); fTimes.cd();
  TH2F *h2 = (TH2F *)(fTimes.Get(hName))->Clone("h2");  
  h2->SetDirectory(0);
  fTimes.Close();
  TFile f(fName,"RECREATE"); 
  f.cd();
  h2->Write();
  f.Close();
  cout<<"Written "<<fName<<endl;
}
