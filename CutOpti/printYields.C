#include "TCut.h"
#include "TString.h"
#include "TChain.h"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace std;
using std::cout;
using std::endl;

TString RoundNumber(double n, int e, double d=1);

void printYields(int candType, TString tag, TString sample, double mvaComb, double mvaDl, TString m2Cut = "yes"){

  TString folder = "small/";
  if(tag.Contains("old")) {tag.ReplaceAll("old",""); folder = "Oldsmall/";}
  TString genName = "AWG82/ntuples/"; genName += folder; genName += tag; genName += "_RunAll.root";
  TString udsName = "AWG82/ntuples/"; udsName += folder; udsName += "uds_RunAll.root";
  TString ccbarName = "AWG82/ntuples/"; ccbarName += folder; ccbarName += "ccbar_RunAll.root";
  TChain gen("ntp1"), uds("ntp1"), ccbar("ntp1");
  gen.Add(genName);
  uds.Add(udsName);
  ccbar.Add(ccbarName);
  TString candCut = "candType==";candCut += candType; 
  if(sample=="sig"){
    candCut+="&&candMvaDl>"; candCut+=mvaDl; candCut+="&&candMvaComb>";candCut+=mvaComb;
    candCut+="&&candPMiss>.2&&candQ2>4";
    if(m2Cut=="yes") candCut += "&&candM2>1.5";
  } else  {
    candCut+="&&candMvaDssDl>"; candCut+=mvaDl; candCut+="&&candMvaDssComb>";candCut+=mvaComb;
    candCut += "&&abs(candCosT)<.8"; 
    if(m2Cut=="yes") candCut += "&&mm2pi0>-0.5&&mm2pi0<1";
  }
  double NB0Run[2][6] = {{34878000, 101690000, 56035000, 166784000, 215168000, 130336000},
			 {44172000, 130056000, 67039000, 216196000, 291011000, 160817000}};
  double NBpRun[2][6] = {{34941000, 104188000, 57888000, 169801000, 215953000, 135224000},
			 {43741000, 128709000, 67786000, 213996000, 245874000, 124453000}};
  double nccbar[] = {55254000, 164146000, 88321000, 267308000, 343667000, 208664000};
  double nuds[] = {160514000, 451636000, 275869000, 421599000, 553604000, 327032000};
  double nMCB = 0, totuds = 0, totccbar = 0, fracRest = 1/3.;
  int index = 0;
  if(tag.Contains("R26")) {
    index = 1;
    fracRest = 1/4.;
  }
  if(tag.Contains("R24")){
    if(tag.Contains("Rests")) fracRest = 2/3.;
  }
  if(tag=="Add") fracRest = 2.;
  for(int i=1; i<7; i++){
    totuds += nuds[i-1];
    totccbar += nccbar[i-1];
    nMCB += NBpRun[index][i-1]+NB0Run[index][i-1]; 
  }
  double wuds = nMCB/totuds*2.09/1.05*fracRest*3/2; 
  double wccbar = nMCB/totccbar*1.3/1.05*fracRest*3/2; 
//   cout<<"totuds "<<RoundNumber(totuds,0,1e6)<<", totccbar "<<RoundNumber(totccbar,0,1e6)<<", nMCB "
//       <<RoundNumber(nMCB,0,1e6)<<", wuds "<<wuds<<", wccbar "<<wccbar<<endl;

  TCut sigCut[] = {"MCType==5","MCType==6", "MCType==11","MCType==12"};
  TCut normCut[] = {"(MCType>0&&MCType<5)","(MCType>0&&MCType<5)", "(MCType>6&&MCType<11)","(MCType>6&&MCType<11)"};
  TCut crossCut[] = {"MCType>6&&MCType<13","MCType>6&&MCType<13", "MCType>0&&MCType<7","MCType>0&&MCType<7"};
  TCut normDssCut[] = {"(MCType>0&&MCType<7)","(MCType>0&&MCType<7)", "(MCType>6&&MCType<13)","(MCType>6&&MCType<13)"};
  for(int i=0; i<4; i++){
    sigCut[i] += candCut;
    normCut[i] += candCut;
    normDssCut[i] += candCut;
    crossCut[i] += candCut;
  }
  TCut DssCut = "MCType>12"; DssCut += candCut;
  TCut combCut = "MCType==0"; combCut += candCut;
  
  double yield[10], total = gen.GetEntries(candCut);
  yield[0] = gen.GetEntries(sigCut[candType-1]);
  yield[1] = gen.GetEntries(normCut[candType-1]);
  yield[2] = gen.GetEntries(DssCut);
  yield[3] = gen.GetEntries(crossCut[candType-1]);
  yield[4] = gen.GetEntries(combCut);
  yield[5] = gen.GetEntries(normDssCut[candType-1]);
  yield[6] = uds.GetEntries(candCut)*wuds;
  yield[7] = ccbar.GetEntries(candCut)*wccbar;
  total += yield[6] + yield[7];
  double fp = 1/6.70384e+08, f0 = 1/6.71486e+08;
  if(tag.Contains("R24")) {fp *= 7.17995e+08; f0 *= 7.04891e+08;}
  if(tag.Contains("R26")) {fp *= (8.24559e+08)*3/4; f0 *= (9.09291e+08)*3/4;}
  if(tag.Contains("R2")){
    for(int i=0; i<8; i++){
      if(candType<3) yield[i] /= fp;
      else yield[i] /= f0;
    }
    if(candType<3) total /= fp;
    else total /= f0;
  }
  if(sample=="sig"){
    cout<<"S/(S+B): "<<RoundNumber(yield[0],2,total)<<", S/sqrt(S+B): "<<RoundNumber(yield[0],2,sqrt(total))<<"    ";
    cout<<"Signal: "<<RoundNumber(yield[0],0)<<", Norm: "<<RoundNumber(yield[1],0) <<", D**: "<<
      RoundNumber(yield[2],0)<<", cross: "<<RoundNumber(yield[3],0)
	<<", Comb: "<<RoundNumber(yield[4],0)<<", uds: "<<RoundNumber(yield[6],0)<<", ccbar: "<<RoundNumber(yield[7],0)<<endl;
  } else {
    cout<<"S/(S+B): "<<RoundNumber(yield[2],2,total)<<", S/sqrt(S+B): "<<RoundNumber(yield[2],2,sqrt(total))<<"    ";
    cout<<"D**: "<<RoundNumber(yield[2],0)<<", Norm: "<<RoundNumber(yield[5],0) <<", cross: "<<RoundNumber(yield[3],0)
	<<", Comb: "<<RoundNumber(yield[4],0)<<", uds: "<<RoundNumber(yield[6],0)<<", ccbar: "<<RoundNumber(yield[7],0)<<endl;
  }
}

TString RoundNumber(double n, int e, double d){
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
