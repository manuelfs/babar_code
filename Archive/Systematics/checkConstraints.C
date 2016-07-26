#include "TString.h"
#include "TChain.h"
#include "TCut.h"
#include "TMath.h"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include "../DonutUtils/cuts.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void checkConstraints(TString extraCut = ""){
  TChain c("ntp1");
  c.Add("AWG82/ntuples/small/RAll_RunAll.root");

  TString MCcuts[2][2] = {{"(MCType==2||MCType==4)","(MCType==8||MCType==10)"},
			  {"MCType==6", "MCType==12"}};
  TString Candcuts[2][2] = {{"candType==1","candType==2"},
			    {"candType==3","candType==4"}};
  TString MVAcuts[2] = {"candMvaDl>0.48","candMvaDl>0.41"};
  TCut cuts[2][2];
  double f[2][2], ef[2][2], n[2][2][2];

  for(int lep=0; lep<2; lep++){
    for(int mc=0; mc<2; mc++){
      for(int cand=0; cand<2; cand++){
	cuts[mc][lep] = basic;//PMiss+M2P;
	cuts[mc][lep] += MCcuts[mc][lep];
	if(cand==0) cuts[mc][lep] += extraCut;
	else cuts[mc][lep] += MVAcuts[cand];
	cuts[mc][lep] += Candcuts[lep][cand];
	n[mc][lep][cand] = c.GetEntries(cuts[mc][lep]);
      }
      f[mc][lep] = n[mc][lep][0]/n[mc][lep][1];
      ef[mc][lep] = f[mc][lep]*sqrt(1/n[mc][lep][0]+1/n[mc][lep][1]);
      cout<<RoundNumber(f[mc][lep],4)<<" +- "<<RoundNumber(ef[mc][lep],4)<<"\t ";
    }
    double rat = f[1][lep]/f[0][lep];
    double erat = rat*sqrt(pow(ef[1][lep]/f[1][lep],2)+pow(ef[0][lep]/f[0][lep],2));
    cout<<"Ratio: "<<RoundNumber(rat,4)<<" +- "<<RoundNumber(erat,4)<<endl;    
  }

}

