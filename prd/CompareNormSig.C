#include "DonutUtils/cuts.cc"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include "TH1F.h"
#include "TChain.h"
#include "TString.h"

using namespace std;
using std::cout;
using std::endl;

void CompareNormSig(TString variable, int doHadr = 1, int mvacuts = 1, int digits = 2){
  TString hName = "histo";
  TH1F *histo = new TH1F(hName,"",100,-50,50);
  TChain c("ntp1");
  c.Add("AWG82/ntuples/small/RAll_RunAll.root");
  TString Names[2][4] = {{"D0 l ", "D*0 l", "D+ l ", "D*+ l"},
			 {"D0 tau ", "D*0 tau", "D+ tau ", "D*+ tau"}};
  TCut NormSigCuts[2][4]  = {{"((MCType==1||MCType==3))",
			      "((MCType==2||MCType==4))",
			      "((MCType==7||MCType==9))",
			      "((MCType==8||MCType==10))"},
			     {"((MCType==5))",
			      "((MCType==6))",
			      "((MCType==11))",
			      "((MCType==12))"}};

  TString HNames[2][4] = {{"D0 tau(->lnunu) ", "D*0 tau(->lnunu)", 
			   "D+ tau(->lnunu) ", "D*+ tau(->lnunu)"},
			  {"D0 tau(->Hadr) ", "D*0 tau(->Hadr)", 
			   "D+ tau(->Hadr) ", "D*+ tau(->Hadr)"}};
  TCut HadrSigCuts[2][4]  = {{"(MCType==5&&MCTaumode>0)",
			      "(MCType==6&&MCTaumode>0)",
			      "(MCType==11&&MCTaumode>0)",
			      "(MCType==12&&MCTaumode>0)"},
			     {"(MCType==5&&MCTaumode==0)",
			      "(MCType==6&&MCTaumode==0)",
			      "(MCType==11&&MCTaumode==0)",
			      "(MCType==12&&MCTaumode==0)"}};


  cout<<endl;
  TCut allCuts[2], chaCuts[2][2];
  double mean[2][7], rms[2][7];
  for(int cut=0; cut<4; cut++){
    TString candCut = "candType=="; candCut += cut+1;
    for(int sig=0; sig<2; sig++){
      TCut thisCut = ""; thisCut += candCut;
      if(doHadr) {thisCut += HadrSigCuts[sig][cut]; Names[sig][cut] = HNames[sig][cut];}
      else thisCut += NormSigCuts[sig][cut]; 
      if(mvacuts) thisCut += MvaAll;
      c.Project(hName, variable, thisCut);
      mean[sig][cut] = histo->GetMean();
      rms[sig][cut] = histo->GetRMS();
      double N = histo->GetEntries();
      cout<<mean[sig][cut]<<", "<<rms[sig][cut]<<", "<<N<<endl;
      if(N) rms[sig][cut] /= sqrt(N);
      else rms[sig][cut] = 999.;
      double signi = (mean[0][cut]-mean[1][cut])/sqrt(pow(rms[0][cut],2)+pow(rms[1][cut],2));
      cout<<Names[sig][cut]<<":  "<<RoundNumber(mean[sig][cut], digits)<<" +- " 
	  <<RoundNumber(rms[sig][cut],digits)<<" \t ";

      if(sig==1) cout<<RoundNumber(signi,2)<<" sigma away";
      if(cut==0) allCuts[sig] = thisCut;
      else allCuts[sig] = allCuts[sig] || thisCut;
      if(cut%2==0) chaCuts[sig][cut>1] = thisCut;
      else chaCuts[sig][cut>1] = chaCuts[sig][cut>1] || thisCut;
    }
    cout<<endl;
  }
  cout<<endl;

  TString ChaNames[2][2][2] = {{{"D(*)0 l", "D(*)0 tau"},{"D(*)0 tau(->lnunu)", "D(*)0 tau(->Hadr)"}},
			    {{"D(*)+ l", "D(*)+ tau"},{"D(*)+ tau(->lnunu)", "D(*)+ tau(->Hadr)"}}};
  for(int cut=0; cut<2; cut++){
    for(int sig=0; sig<2; sig++){
      int index = 6+cut;
      //chaCuts[sig][cut].Print();
      c.Project(hName, variable, chaCuts[sig][cut]);
      mean[sig][index] = histo->GetMean(); 
      rms[sig][index] = histo->GetRMS();
      double N = histo->GetEntries();
      cout<<mean[sig][index]<<", "<<rms[sig][index]<<", "<<N<<endl;
      if(N) rms[sig][index] /= sqrt(N);
      else rms[sig][index] = 999.;
      double signi = (mean[0][index]-mean[1][index])/sqrt(pow(rms[0][index],2)+pow(rms[1][index],2));
      cout<<ChaNames[cut][doHadr][sig]<<":  "<<RoundNumber(mean[sig][index], digits)<<" +- " 
	  <<RoundNumber(rms[sig][index],digits)<<" \t ";
      if(sig==1) cout<<RoundNumber(signi,2)<<" sigma away";
    }
    cout<<endl;
  }
  cout<<endl;

  TString ANames[2][2] = {{"D(*) l", "D(*) tau"},{"D(*) tau(->lnunu)", "D(*) tau(->Hadr)"}};
  for(int sig=0; sig<2; sig++){
    c.Project(hName, variable, allCuts[sig]);
    mean[sig][5] = histo->GetMean(); 
    rms[sig][5] = histo->GetRMS();
    double N = histo->GetEntries();
    if(N) rms[sig][5] /= sqrt(N);
    else rms[sig][5] = 999.;
    double signi = (mean[0][5]-mean[1][5])/sqrt(pow(rms[0][5],2)+pow(rms[1][5],2));
    cout<<ANames[doHadr][sig]<<":  "<<RoundNumber(mean[sig][5], digits)<<" +- " 
	<<RoundNumber(rms[sig][5],digits)<<" \t ";
    if(sig==1) cout<<RoundNumber(signi,2)<<" sigma away";
  }
  cout<<endl<<endl;
  histo->Delete();
 
}
