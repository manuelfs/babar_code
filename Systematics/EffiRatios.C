#include "TString.h"
#include "TChain.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

TString round(double n, int e, double d=1.);

void EffiRatios(int nFiles, int doW = 1){

  TString folder = "AWG82/systematics/BF/FitRoot/Fit"; 
  TString channels[] = {"$D^0$","$D^{*0}$","$D^+$","$D^{*+}$","$D$","$D^{*}$"};
  double BFratio[] = {0.0070/0.0224, 0.0160/0.0617, 0.0070/0.0207, 0.0160/0.0570, 
		      (0.0070+0.0070)/(0.0224+0.0207), (0.0160+0.0160)/(0.0617+0.0570)};  // SP8 values
  if(doW){
    BFratio[0] = 0.3;  BFratio[2] = 0.3;            // Re-weighted values
    BFratio[1] = 0.25; BFratio[3] = 0.25;           // Re-weighted values
    BFratio[4] = 0.3 ; BFratio[5] = 0.25;           // Re-weighted values
  }
  float weight;
  int candType, MCType;
  double mean[6], rms[6], meanBkg[8][2], rmsBkg[8][2], Bkg[8][2];
  for(int i=0; i<6; i++){mean[i] = 0; rms[i] = 0;}
  for(int cand=0; cand<8; cand++){
    for(int i=0; i<2; i++){
      meanBkg[cand][i] = 0; rmsBkg[cand][i] = 0;
    }
  }
  for(int file=0; file<nFiles; file++){
    TString genName = folder; genName += file; genName += "RAll_RunAll.root";
    TChain *genChain = new TChain("ntp1");
    genChain->Add(genName);
    double totW[2][6], totW2[2][6];
    for(int i=0; i<4; i++){
      for(int j=0; j<2; j++){
	totW[j][i] = 0; totW2[j][i] = 0; 
      }
    }
    for(int cand=0; cand<8; cand++)for(int i=0; i<2; i++)Bkg[cand][i] = 0;
    genChain->SetBranchAddress("weight",&weight);
    genChain->SetBranchAddress("candType",&candType);
    genChain->SetBranchAddress("MCType",&MCType);
    for(int entry=0; entry<genChain->GetEntries(); entry++){
      genChain->GetEvent(entry);
      if(doW==0) weight = 1;

      //if(MCType==5 || MCType==6 || MCType==11 || MCType==12) weight = 0.1; //To check the effect of the D(*)TauNu FF 

      if(candType<=2 && MCType==5) {totW[0][0] += weight; totW2[0][0] += weight*weight;}
      if(candType<=2 && MCType==6) {totW[0][1] += weight; totW2[0][1] += weight*weight;}
      if(candType>=3 && candType<=4 && MCType==11){totW[0][2] += weight; totW2[0][2] += weight*weight;}
      if(candType>=3 && candType<=4 && MCType==12){totW[0][3] += weight; totW2[0][3] += weight*weight;}

      if(candType<=2 && (MCType==1||MCType==3)) {totW[1][0] += weight; totW2[1][0] += weight*weight;}
      if(candType<=2 && (MCType==2||MCType==4)) {totW[1][1] += weight; totW2[1][1] += weight*weight;}
      if(candType>=3 && candType<=4 && (MCType==7||MCType==9)) {totW[1][2] += weight; totW2[1][2] += weight*weight;}
      if(candType>=3 && candType<=4 && (MCType==8||MCType==10)){totW[1][3] += weight; totW2[1][3] += weight*weight;}
      for(int cand=0; cand<8; cand++){
	if(MCType==0 && candType==cand+1) Bkg[cand][0] += weight;
	if(MCType<0 && candType==cand+1) Bkg[cand][1] += weight;
      }
    }
    for(int chan = 0; chan<2; chan++){
      for(int j=0; j<2; j++){
	totW[j][chan+4] = totW[j][chan]+totW[j][chan+2];
	totW2[j][chan+4] = totW2[j][chan]+totW2[j][chan+2];
      }
    }
    for(int i=0; i<6; i++){
      double n = totW[0][i];
      double N = totW[1][i];
      double ratio = -1;
      if(N!=0) {
	ratio = n/N/BFratio[i]*2;
      }
      mean[i] += ratio;
      rms[i] += ratio*ratio;
    }
    for(int cand=0; cand<8; cand++){
      for(int i=0; i<2; i++){
	meanBkg[cand][i] += Bkg[cand][i]; rmsBkg[cand][i] += Bkg[cand][i]*Bkg[cand][i];
      }
    }
    genChain->Delete();
  }
  for(int i=0; i<6; i++){
    mean[i] /= (double)nFiles;
    rms[i] = sqrt((rms[i] - mean[i]*mean[i]*nFiles)/(nFiles-1));
    cout<<"Ratio "<<channels[i]<<" is     \t"<<round(mean[i],4)<<" +- "<<round(rms[i],4)<<",\t a "<<
      round(rms[i]/mean[i]*100,2)<<" % error"<<endl;
  }
  cout<<endl;
  for(int i=0; i<2; i++){
    for(int cand=0; cand<8; cand++){
      meanBkg[cand][i] /= (double)nFiles;
      rmsBkg[cand][i] = sqrt((rmsBkg[cand][i] - meanBkg[cand][i]*meanBkg[cand][i]*nFiles)/(nFiles-1));
      cout<<"Bkg "<<channels[cand-4*(cand>3)]<<" is     \t"<<round(meanBkg[cand][i],1)<<" +- "<<round(rmsBkg[cand][i],1)<<",\t a "<<
	round(rmsBkg[cand][i]/meanBkg[cand][i]*100,2)<<" % error"<<endl;
    }
    cout<<endl;
  }
  cout<<endl;
}

TString round(double n, int e, double d){
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
