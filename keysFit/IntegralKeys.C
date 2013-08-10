#include "babar_code/PlotsThesis/PlotUtils.cc"
#include "TString.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TKey.h"
#include "TCanvas.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void IntegralKeys(double minM2=2.5, double maxM2=12, int maxRep=100){
  
  int delta = 9;
  double Yield[2][70],  Error[70];
  ReadFitFile("FitAll/fits/TextFinalDataNe2x100.txt", Yield, Error);
  TString chaName[2][4] = {{"\\Dz", "\\Dstarz", "\\Dp", "\\Dstarp"},{"\\Dstarz", "\\Dz", "\\Dstarp", "\\Dp"}};
  TString hName, rArrow = "\\Rightarrow";
  //TCanvas can; can.Divide(2,2);
  //TH1F *histo[4];
  //double MeanRMS[2][4] = {{0.0132, 0.0104, 0.0223, 0.0069},{0.0021, 0.0016, 0.0035, 0.0011}};
  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4, plCut = 1.4;
  TKey *key;
  cout<<endl;
  for(int chan=0; chan<4; chan++){
    for(int isXf=0; isXf<2; isXf++){
      if(isXf==0) delta = 9;
      else delta = 13;
      hName = "hname"; hName += chan+delta; 
      //double minH = MeanRMS[0][chan]-3*MeanRMS[1][chan], maxH = MeanRMS[0][chan]+3*MeanRMS[1][chan];
      //histo[chan] = new TH1F(hName,"", 50, minH*100, maxH*100);
      hName = "keys/Archive/12-05-01/Times/hTimes_"; hName += chan+delta; hName += "_Fit.root";
      TFile File(hName);
      TIter Fiter(File.GetListOfKeys());
      double mean = 0, rms = 0, meanTot=0; int rep=0;
      while ( (key = (TKey*)Fiter()) && rep<maxRep) {
	TH2F *h2 = (TH2F*)key->ReadObj();
	int nM2bin = h2->GetNbinsX(), nPlbin = h2->GetNbinsY(); 
	int iniBin = (int)((minM2-m2min)*(double)nM2bin/(m2max-m2min))+1;
	int finBin = (int)((maxM2-m2min)*(double)nM2bin/(m2max-m2min));
	int plBin  = (int)((plCut-plmin)*(double)nPlbin/(plmax-plmin));
	double fraction = h2->Integral(iniBin,finBin,1,plBin)/h2->Integral(1,nM2bin,1,plBin);
	mean += fraction; rms += fraction*fraction;
	meanTot += h2->Integral(1,nM2bin,1,plBin)/h2->Integral();
	//histo[chan] ->Fill(fraction*100);
	rep++;
      }
      double N = (double)(rep);
      meanTot /= N; mean /= N; rms = sqrt((rms - mean*mean*N)/(N-1));
      //     cout<<"PDF "<<chan+delta<<":\tFraction is ("<<RoundNumber(mean*100,2)<<" +- "<<RoundNumber(rms*100,2)<<")%   =>   "
      // 	<<RoundNumber(mean*Yield[0][chan+delta],1)<<" +- "<<RoundNumber(Error[chan+delta]*mean,1)
      // 	<<" +- "<<RoundNumber(rms*Yield[0][chan+delta],1)
      // 	<<" events out of "<<RoundNumber(Yield[0][chan+delta],1)<<" +- "<<RoundNumber(Error[chan+delta],1)<<endl;
      int iSig = chan+1+isXf*4;
      double Ytail = mean*Yield[0][chan+delta]*meanTot, Ystat = Error[chan+delta]*mean*meanTot;
      double Ysyst = rms*Yield[0][chan+delta]*meanTot;
      double sigY = sqrt(pow(Ysyst,2)+pow(Ystat,2))/Ytail*100;
      cout<<"$"<<chaName[isXf][chan]<<rArrow<<chaName[0][chan];
      if(!(isXf==0&&chan%2==1))cout<<"\t";
      cout<<"$\t& "<<RoundNumber(Yield[0][iSig],0)
	  <<" & "<<RoundNumber(Error[iSig],0)
	  <<" \t& "<<RoundNumber(Yield[0][chan+delta],0)<<" & "<<RoundNumber(Error[chan+delta],0)
	  <<" \t& "<<RoundNumber(mean*100,2)<<" & "<<RoundNumber(rms*100,2)
	  <<" \t& "<<RoundNumber(Ytail,1)<<" & "<<RoundNumber(Ystat,1)
	  <<" & "<<RoundNumber(Ysyst,1)<<" \t& "<<RoundNumber(sigY,1)
	  <<" \t& "<<RoundNumber(mean*Yield[0][chan+delta]*100*meanTot,1,Yield[0][chan+1]+Yield[0][chan+5])<<" \\\\"<<endl;
      //can.cd(chan+1);
      //histo[chan]->Draw();
    }
    cout<<"\\hline"<<endl;
  }
  cout<<endl;
  //can.SaveAs("Integral.eps");
  //for(int chan=0; chan<4; chan++)histo[chan] ->Delete();

}
