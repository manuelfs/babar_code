#include "TString.h"
#include "TChain.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

TString round(double n, int e, double d=1.);

void EffiRatios(TString base = "Newx100", TString Runs = "All", int doW = 0){
  TString genName = "AWG82/ntuples/small/FitRAll"; genName+=base; genName+="_Run"; genName+=Runs; genName+=".root";
  TChain genChain("ntp1");
  genChain.Add(genName);
  double totW[2][6], totW2[2][6], totN[2][2][6];
  for(int i=0; i<4; i++){
    for(int j=0; j<2; j++){
      totW[j][i] = 0; totW2[j][i] = 0; 
      for(int k=0; k<2; k++) totN[k][j][i] = 0;
    }
  }
  float weight, candPisoftP;
  int candType, candDstarType, MCType;
  genChain.SetBranchAddress("weight",&weight);
  genChain.SetBranchAddress("candType",&candType);
  genChain.SetBranchAddress("candDstarType",&candDstarType);
  genChain.SetBranchAddress("MCType",&MCType);
  genChain.SetBranchAddress("candPisoftP",&candPisoftP);
  double Bdata = 21.0713, Bmc = 17.5988, p0data = 0.0269, p0mc = 0.026;
  for(int entry=0; entry<genChain.GetEntries(); entry++){
    genChain.GetEvent(entry);
    if(doW) {
      if(candType==4 && candDstarType==1){
	if(candPisoftP<p0mc) weight = 0;
	else weight *= (1-1/(Bdata*(candPisoftP-p0data)+1))/(1-1/(Bmc*(candPisoftP-p0mc)+1));      
      }
      if((candType==4&& candDstarType==2)||(candType==2 && candDstarType==2))
	weight *= 0.9395 + 0.0118*candPisoftP;
    }

    if(candType<=2 && MCType==5) {
      totW[0][0] += weight; totW2[0][0] += weight*weight;
      if(candType==1) totN[0][0][0] += weight;
      else totN[1][0][0] += weight;
    }
    if(candType<=2 && MCType==6) {
      totW[0][1] += weight; totW2[0][1] += weight*weight;
      if(candType==2) totN[0][0][1] += weight;
      else totN[1][0][1] += weight;
    }
    if(candType>=3 && candType<=4 && MCType==11){
      totW[0][2] += weight; totW2[0][2] += weight*weight;
      if(candType==3) totN[0][0][2] += weight;
      else totN[1][0][2] += weight;
    }
    if(candType>=3 && candType<=4 && MCType==12){
      totW[0][3] += weight; totW2[0][3] += weight*weight;
      if(candType==4) totN[0][0][3] += weight;
      else totN[1][0][3] += weight;
    }

    if(candType<=2 && (MCType==1||MCType==3)) {
      totW[1][0] += weight; totW2[1][0] += weight*weight;
      if(candType==1) totN[0][1][0] += weight;
      else totN[1][1][0] += weight;
    }
    if(candType<=2 && (MCType==2||MCType==4)) {
      totW[1][1] += weight; totW2[1][1] += weight*weight;
      if(candType==2) totN[0][1][1] += weight;
      else totN[1][1][1] += weight;
    }
    if(candType>=3 && candType<=4 && (MCType==7||MCType==9)) {
      totW[1][2] += weight; totW2[1][2] += weight*weight;
      if(candType==3) totN[0][1][2] += weight;
      else totN[1][1][2] += weight;
    }
    if(candType>=3 && candType<=4 && (MCType==8||MCType==10)){
      totW[1][3] += weight; totW2[1][3] += weight*weight;
      if(candType==4) totN[0][1][3] += weight;
      else totN[1][1][3] += weight;
    }
  }
  for(int chan = 0; chan<2; chan++){
    for(int j=0; j<2; j++){
      totW[j][chan+4] = totW[j][chan]+totW[j][chan+2];
      totW2[j][chan+4] = totW2[j][chan]+totW2[j][chan+2];
      //for(int k=0; k<2; k++) totN[k][j][i] += totN[k][j][chan]+totN[k][j][chan+2];
    }
  }

  double BFratio[] = {0.0070/0.0224, 0.0160/0.0617, 0.0070/0.0207, 0.0160/0.0570, 
		      (0.0070+0.0070)/(0.0224+0.0207), (0.0160+0.0160)/(0.0617+0.0570)};  // SP8 values
  //if(doW){
    BFratio[0] = 0.3;  BFratio[2] = 0.3;            // Re-weighted values
    BFratio[1] = 0.25; BFratio[3] = 0.25;           // Re-weighted values
    BFratio[4] = 0.3 ; BFratio[5] = 0.25;           // Re-weighted values
    //}
  TString texName = "babar_code/Tables/texEffiRatios.txt";
  fstream tex;
  tex.open(texName,fstream::out);

  TString texType = "\\bf{Signal }";  
  tex<<"\\begin{tabular}{lc}"<<endl<<"\\hline\\hline"<<endl;
  tex<<texType<<" & $\\epsilon_{\\text{sig}}/\\epsilon_{\\text{norm}}$ \\\\ \\hline"<<endl;
  TString channels[] = {"$D^0$","$D^{*0}$","$D^+$","$D^{*+}$","$D$","$D^{*}$"};
  for(int i=0; i<6; i++){
    double n = totW[0][i];
    double N = totW[1][i];
    double ratio = -1, err = -1;
    if(N!=0) {
      ratio = n/(N/2.)/BFratio[i]/(0.1778+0.1731); 
      ratio = n/N/BFratio[i];
      err = ratio*sqrt(totW2[0][i]/n/n+totW2[1][i]/N/N);
      
    }
    cout<<"Ratio "<<channels[i]<<" is     \t"<<round(ratio,6)<<" +- "<<round(err,4)<<",\t a "<<
      round(err/ratio*100,2)<<" % error"<<endl;
    tex<<channels[i]<<" & "<<round(ratio,3)<<" $\\pm$ "<<round(err,3)<<" \\\\"<<endl;
    if(i==3)cout<<endl;
  }
  tex<<"\\hline\\hline \\end{tabular}\\,\\,"<<endl;
  //cout<<endl<<texName<<" done"<<endl<<endl;
  
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
