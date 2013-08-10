#include "TString.h"
#include "TChain.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

TString round(double n, int e, double d=1.);

void ratiosEff(int doW = 0){
  TString genName = "AWG82/ntuples/small/FitRAllW1_RunAll.root";
  TChain genChain("ntp1");
  genChain.Add(genName);
  double totW[2][4], totW2[2][4];
  for(int i=0; i<4; i++){
    for(int j=0; j<2; j++){
      totW[j][i] = 0; totW2[j][i] = 0; 
    }
  }
  float weight;
  int candType, MCType;
  genChain.SetBranchAddress("weight",&weight);
  genChain.SetBranchAddress("candType",&candType);
  genChain.SetBranchAddress("MCType",&MCType);
  for(int entry=0; entry<genChain.GetEntries(); entry++){
    genChain.GetEvent(entry);
    if(doW==0) weight = 1;

    if(candType<=2 && MCType==5) {totW[0][0] += weight; totW2[0][0] += weight*weight;}
    if(candType<=2 && MCType==6) {totW[0][1] += weight; totW2[0][1] += weight*weight;}
    if(candType>=3 && candType<=4 && MCType==11){totW[0][2] += weight; totW2[0][2] += weight*weight;}
    if(candType>=3 && candType<=4 && MCType==12){totW[0][3] += weight; totW2[0][3] += weight*weight;}

    if(candType<=2 && (MCType==1||MCType==3)) {totW[1][0] += weight; totW2[1][0] += weight*weight;}
    if(candType<=2 && (MCType==2||MCType==4)) {totW[1][1] += weight; totW2[1][1] += weight*weight;}
    if(candType>=3 && candType<=4 && (MCType==7||MCType==9)) {totW[1][2] += weight; totW2[1][2] += weight*weight;}
    if(candType>=3 && candType<=4 && (MCType==8||MCType==10)){totW[1][3] += weight; totW2[1][3] += weight*weight;}
  }

  double BFratio[] = {0.313, 0.259, 0.338, 0.281};  // SP8 values
  if(doW){
    BFratio[0] = 0.3;  BFratio[2] = 0.3;            // Re-weighted values
    BFratio[1] = 0.25; BFratio[3] = 0.25;           // Re-weighted values
  }
  TString texName = "babar_code/Tables/texRatio";texName+="DEff"; texName += ".txt";
  fstream tex;
  tex.open(texName,fstream::out);

  TString texType = "\\bf{Signal }";  
  tex<<"\\begin{tabular}{lc}"<<endl<<"\\hline\\hline"<<endl;
  tex<<texType<<" & $\\epsilon_{\\text{sig}}/\\epsilon_{\\text{norm}}$ \\\\ \\hline"<<endl;
  TString channels[] = {"$D^0$","$D^{*0}$","$D^+$","$D^{*+}$"};
  for(int i=0; i<4; i++){
    double n = totW[0][i];
    double N = totW[1][i];
    double ratio = -1, err = -1;
    if(N!=0) {
      ratio = n/N/BFratio[i]*2;
      err = ratio*sqrt(totW2[0][i]/n/n+totW2[1][i]/N/N);
    }
    cout<<"Ratio "<<channels[i]<<" is     \t"<<round(ratio,3)<<" +- "<<round(err,3)<<endl;
    tex<<channels[i]<<" & "<<round(ratio,3)<<" $\\pm$ "<<round(err,3)<<" \\\\"<<endl;
  }
  tex<<"\\hline\\hline \\end{tabular}\\,\\,"<<endl;
  cout<<texName<<" done"<<endl;
  
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
