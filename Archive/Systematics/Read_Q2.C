#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace TMath;
using namespace std;
using std::cout;
using std::endl;

void Read_Q2(TString hisTag = "BB", TString hName = "public_html/Test_Histograms.root", 
	     int row=0, int col=0){
  TFile saveHistos(hName); saveHistos.cd();
  TH1F *hTemp;

  int nRows = 2, nBins = 17, dig=1;
  int pad = nRows*col+row;
  cout<<"Reading "<<hisTag<<" from "<<hName<<endl<<"============================="<<endl;
  for(int bin=1; bin<=nBins; bin++){
    hName = hisTag; hName += "q2_"; hName += pad; hName += "_"; hName += bin;
    hTemp = (TH1F *)saveHistos.Get(hName);
    if(hTemp==0) {cout<<"No "<<hName<<endl; return;}
    int nBinsH = hTemp->GetNbinsX();
    double expY = hTemp->GetBinLowEdge(1)/0.3;
    if(hisTag=="Data") expY = hTemp->GetBinLowEdge(1)+50;
    double meanY = hTemp->GetMean(), errY = hTemp->GetRMS();
    cout<<bin<<": \t"<<RoundNumber(expY,dig)<<"\t "<<RoundNumber(meanY,dig)
	<<" +- "<<RoundNumber(errY,dig)<<"\t a "<<RoundNumber(meanY-expY,dig)
	<<" diff \t "<<RoundNumber(errY*100,1,expY)
	<<"% error \t Underflow "<<hTemp->GetBinContent(0)<<" \t Overflow "
	<<hTemp->GetBinContent(nBinsH+1)<<endl;
  }


}
