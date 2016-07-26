#include "TMath.h"
#include "TString.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLegend.h"
#include "DonutUtils/cuts.cc"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;


void PlotMult(){
  TH1F *histo[3][2][2];
  TString hName, Bname[] = {"B0","Bp"}, Varname[] = {"Mult_", "Char_", "Neu_"};
  TString Typename[] = {"Nor_", "Sig_"};
  TString variable[] = {"candBnCharged+candBnNeutral","candBnCharged", "candBnNeutral"};
  TString options[] = {"hist","same"};
  int color[]={4,2};
  TString chaName[] = {"B^{+}_{tag} ", "B^{0}_{tag} "};
  TString varName[] = {"total multiplicity", "charged multiplicity", "neutral multiplicity"};
  TChain tree("ntp1");
  tree.Add("AWG82/ntuples/small/RAll_RunAll.root");
  
  gStyle->SetOptStat(0);
  double chi2[3][2];
  int ndof[3][2];
  TCut NormSigCuts[2][2]  = {{"(candType==1&&(MCType==1||MCType==3)||candType==2&&(MCType==2||MCType==4))",
			      "(candType==3&&(MCType==7||MCType==9)||candType==4&&(MCType==8||MCType==10))"},
			     {"(candType==1&&MCType==5||candType==2&&MCType==6)",
			      "(candType==3&&MCType==11||candType==4&&MCType==12)"}};

  TCanvas can("can","Multiplicities",800,1200); can.Divide(2,3);
  for(int var=0; var<3; var++)
    for(int typ=0; typ<2; typ++)
      for(int his=0; his<2; his++) {
  	hName = Varname[var]; hName += Bname[his]; hName += Typename[typ];
  	histo[var][his][typ] = new TH1F(hName, "", 12, 0, 12);
  	TCut thisCut = MvaAll; thisCut += NormSigCuts[typ][his];
  	tree.Project(hName, variable[var], thisCut);
  	histo[var][his][typ]->SetLineColor(color[typ]);
  	can.cd(his+2*var+1);
  	histo[var][his][typ]->Sumw2();
  	histo[var][his][typ]->Scale(10./histo[var][his][typ]->Integral());
  	histo[var][his][typ]->Draw(options[typ]);
	chi2[var][his] = 0; ndof[var][his] = -1;
	if(typ==1){
	  for(int bin=1; bin<=histo[var][his][typ]->GetNbinsX(); bin++){
	    double val1 = histo[var][his][1]->GetBinContent(bin);
	    double val2 = histo[var][his][0]->GetBinContent(bin);
	    double err1 = histo[var][his][1]->GetBinError(bin);
	    double err2 = histo[var][his][0]->GetBinError(bin);
	    double err = pow(err1,2)+pow(err2,2);
	    if(err){
	      chi2[var][his] += pow(val1-val2,2)/err;
	      ndof[var][his]++;
	    }
	  }
	  cout<<"chi2 "<<chi2[var][his]<<", ndof "<<ndof[var][his]<<", p "
	      <<RoundNumber(TMath::Prob(chi2[var][his], ndof[var][his])*100,1)<<" %"<<endl;
	  TString title = "#chi^{2}: "; title += RoundNumber(chi2[var][his],1);
	  title += "/"; title += ndof[var][his]; title += " #rightarrow p = ";
	  title += RoundNumber(TMath::Prob(chi2[var][his], ndof[var][his])*100,1);
	  title += "%";
	  histo[var][his][0]->SetTitle(title);
	} else {
	  TString xTitle = chaName[his]; xTitle += varName[var];
	  histo[var][his][typ]->SetTitleSize(0.05,"XY");
	  histo[var][his][typ]->SetXTitle(xTitle);
	  histo[var][his][typ]->SetYTitle("10 #times Density");
	}
      }
  double legW = 0.27, legH = 0.15;
  double legX = 0.8, legY = 0.86;
  TLegend leg(legX-legW, legY-legH, legX, legY);
  leg.SetTextSize(0.05); leg.SetFillColor(0); 
  leg.SetTextFont(132); leg.SetBorderSize(0);
  leg.AddEntry(histo[0][0][0],"Normalization");
  leg.AddEntry(histo[0][0][1],"Signal");
  can.cd(5);
  leg.Draw();
  can.SaveAs("multi.gif");

  for(int var=0; var<3; var++)
    for(int his=0; his<2; his++)  
      for(int typ=0; typ<2; typ++) histo[var][his][typ]->Delete();
}




