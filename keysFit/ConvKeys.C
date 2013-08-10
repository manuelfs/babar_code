#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TH2F.h"
#include "TSystem.h"
#include "TROOT.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;
double IntG(double mean, double sigma, double minX, double maxX);
void ConvKeys(TH2F *h2, double width); 

void ConvKeys(){
  TString inFolder = "keys/root/fitFinal/";
  TString outFolder = "keys/root/Fit/";

  double width[4];
  fstream ResolFile; ResolFile.open("../DonutUtils/ResolutionWidths.txt",fstream::in);
  for(int chan=0; chan<4; chan++) ResolFile>>width[chan];
  void *dir = gSystem->OpenDirectory(inFolder);
  int nfiles = 0;
  TString file;
  if (dir) {
    while ((file = gSystem->GetDirEntry(dir))){// && nfiles<4) {
      if (file=="") break;
      //if (!file.Contains("_10_Fit.root")) continue;
      if (!file.Contains("_Fit.root")) continue;
      file = inFolder + file;
      cout<<file<<endl;

      TString outFile(file);
      outFile.Remove(0,outFile.Last('/')+1);      // Remove what comes before the /
      TString sample = outFile;
      sample.Remove(sample.Last('_'),sample.Length());      
      sample.Remove(0,sample.Last('_')+1);  
      int sam = sample.Atoi();
      outFile = (outFolder + outFile);
      TFile hIn(file);
      TH2F *h2 = (TH2F *)(hIn.Get("h2"));
      TFile hOut(outFile,"RECREATE"); hOut.cd();
      if(sam>=9&&sam<=20 || sam>=41&&sam<=44) ConvKeys(h2, width[(sample.Atoi()-1)%4]);
      h2->Write();
      hIn.Close();
      hOut.Close();
      nfiles++;
    }
  }
}



double IntG(double mean, double sigma, double minX, double maxX){
  return (TMath::Erf((maxX-mean)/sigma/sqrt(2.))-TMath::Erf((minX-mean)/sigma/sqrt(2.)))/2.;
}


void ConvKeys(TH2F *h2, double width){
  int nM2bin = h2->GetNbinsX(), nPlbin = h2->GetNbinsY(); 
  double minX = h2->GetXaxis()->GetBinLowEdge(h2->GetXaxis()->GetFirst());
  double maxX = h2->GetXaxis()->GetBinLowEdge(h2->GetXaxis()->GetLast()+1);

  int PointsBin = 1;
  double nSigma = 3., val[2000];
  double wBin = (maxX-minX)/(double)nM2bin;
  double dx = wBin/((double)(PointsBin+1)), totAreaG = IntG(0,1,-nSigma,nSigma);


  for(int bin=1; bin<=nPlbin; bin++){
    for(int mbin=1; mbin<=nM2bin; mbin++) val[mbin-1] = h2->GetBinContent(mbin,bin);
    for(int mbin=1; mbin<=nM2bin; mbin++) {
      double valC = 0;
      for(int point=0; point<PointsBin; point++){
	double x = minX+(double)(mbin-1)*wBin+(double)(point+1)*dx;
	double minG = x-(double)nSigma*width, maxG = x+(double)nSigma*width;
	double AreaG = totAreaG, valH = 0;
	bool lessRange = false;
	if(minG<minX) {minG=minX; lessRange = true;} if(maxG>maxX) {maxG=maxX; lessRange = true;} 
	if(lessRange) AreaG = IntG(x,width,minG,maxG);
	int iniBin = (int)((minG-minX)/wBin)+1, finBin = (int)((maxG-minX)/wBin)+1;
	for(int binH=iniBin; binH<=finBin; binH++){
	  double minH = (double)(binH-1)*wBin+minX, maxH = (double)(binH)*wBin+minX;
	  if(minH<minG) minH=minG; if(maxH>maxG) maxH=maxG; 
	  valH += IntG(x, width, minH, maxH)*val[binH-1];
	}
	valC += valH/(AreaG*(double)PointsBin);
      }
      h2->SetBinContent(mbin,bin,valC);
    }
  }
}








