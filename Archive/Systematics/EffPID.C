#include "TString.h"
#include "TChain.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

#define nPbinsMax 50

int indexComp(int candType, int MCType);
TString round(double n, int e, double d=1.);

void EffPID(int TypeTable, TString base = "Nom2"){
  TString genName = "AWG82/ntuples/small/FitRAll"; genName+=base; genName+="_RunAll.root";
  TChain genChain("ntp1");
  genChain.Add(genName);
  double totW[2][12];
  for(int i=0; i<2; i++)
    for(int j=0; j<12; j++) totW[i][j] = 0; 

  TString PIDBase = "babar_code/Systematics/Text/RatioPID", buffer;
  fstream PIDFile[2]; 
  double Pweight[2][nPbinsMax], Prange[2][nPbinsMax]; int nPbins[2] = {0,0};
  for(int emu=0; emu<2; emu++){
    TString PIDName = PIDBase; PIDName += (2*TypeTable+emu); PIDName += ".txt";
    //cout<<PIDName<<endl;
    PIDFile[emu].open(PIDName,fstream::in);
    while(PIDFile[emu] && nPbins[emu]<50){
      PIDFile[emu]>>Prange[emu][nPbins[emu]]>>buffer>>buffer>>Pweight[emu][nPbins[emu]];
      //cout<<Prange[emu][nPbins[emu]]<<"\t"<<Pweight[emu][nPbins[emu]]<<endl;
      //Pweight[emu][nPbins[emu]] = 1;
      nPbins[emu]++;
    }
    nPbins[emu]--;
  }
  Prange[0][nPbins[0]] = 10.;
  Prange[1][nPbins[1]] = 99.;
  float weight, candPLep;
  int candType, MCType, candIsMu;
  genChain.SetBranchAddress("weight",&weight);
  genChain.SetBranchAddress("candPLep",&candPLep);
  genChain.SetBranchAddress("candType",&candType);
  genChain.SetBranchAddress("MCType",&MCType);
  genChain.SetBranchAddress("candIsMu",&candIsMu);
  for(int entry=0; entry<genChain.GetEntries(); entry++){
    genChain.GetEvent(entry);
    int index = indexComp(candType, MCType);
    if(index>=0){
      int Weighted = 0;
      for(int p=0; p<nPbins[candIsMu]; p++){
	if(candPLep>=Prange[candIsMu][p] && candPLep<Prange[candIsMu][p+1]) {
	  if(Weighted==1) cout<<candPLep<<"\t"<<candIsMu<<"\t"<<index<<endl;
	  totW[0][index] += weight*Pweight[candIsMu][p];
	  totW[1][index] += weight*weight*Pweight[candIsMu][p]*Pweight[candIsMu][p];
	  Weighted = 1;
	}
      }
      if(Weighted==0) cout<<candPLep<<"\t"<<candIsMu<<"\t"<<index<<endl;
    }
  }


  for(int i=0; i<2; i++)
    for(int j=4; j<6; j++) {
      totW[i][j] = totW[i][j-4]+totW[i][j-2]; 
      totW[i][j+6] = totW[i][j+2]+totW[i][j+4]; 
    }

  double BFratio[] = {0.0070/0.0224, 0.0160/0.0617, 0.0070/0.0207, 0.0160/0.0570, 
		      (0.0070+0.0070)/(0.0224+0.0207), (0.0160+0.0160)/(0.0617+0.0570)};  // SP8 values
  BFratio[0] = 0.3;  BFratio[2] = 0.3;            // Re-weighted values
  BFratio[1] = 0.25; BFratio[3] = 0.25;           // Re-weighted values
  BFratio[4] = 0.3 ; BFratio[5] = 0.25;           // Re-weighted values
  TString texName = "babar_code/Tables/texEffiRatios2.txt";
  fstream tex;
  tex.open(texName,fstream::out);

  TString texType = "\\bf{Signal }";  
  tex<<"\\begin{tabular}{lc}"<<endl<<"\\hline\\hline"<<endl;
  tex<<texType<<" & $\\epsilon_{\\text{sig}}/\\epsilon_{\\text{norm}}$ \\\\ \\hline"<<endl;
  TString channels[] = {"$D^0$","$D^{*0}$","$D^+$","$D^{*+}$","$D$","$D^{*}$"};
  double NomEff[] = {0.747491, 0.475598, 0.768119, 0.439973, 0.754137, 0.465025 };
  for(int i=0; i<6; i++){
    double n = 0, N = 0, n2 = 0, N2 = 0;
    n = totW[0][i];
    N = totW[0][i+6];
    n2 = totW[1][i];
    N2 = totW[1][i+6];
    double ratio = -1, err = -1;
    if(N!=0) {
      ratio = n/N/BFratio[i]*2;///(0.1778+0.1731);
      err = ratio*sqrt(n2/n/n+N2/N/N);
    }
    cout<<"Ratio "<<channels[i]<<" is     \t"<<round(ratio,4)<<" +- "<<round(err,4)<<",\t a "<<
      round((ratio-NomEff[i])/NomEff[i]*100,2)<<" % error"<<endl;
    tex<<channels[i]<<" & "<<round(ratio,3)<<" $\\pm$ "<<round(err,3)<<" \\\\"<<endl;
    if(i==3)cout<<endl;
  }
  tex<<"\\hline\\hline \\end{tabular}\\,\\,"<<endl;
  //cout<<endl<<texName<<" done"<<endl<<endl;
  
}

int indexComp(int candType, int MCType){
  if(candType<=2 && MCType==5) return 0;
  if(candType<=2 && MCType==6) return 1;
  if(candType>=3 && candType<=4 && MCType==11) return 2;
  if(candType>=3 && candType<=4 && MCType==12) return 3;

  if(candType<=2 && (MCType==1||MCType==3)) return 6;
  if(candType<=2 && (MCType==2||MCType==4)) return 7;
  if(candType>=3 && candType<=4 && (MCType==7||MCType==9)) return 8;
  if(candType>=3 && candType<=4 && (MCType==8||MCType==10)) return 9;

  return -1;
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
