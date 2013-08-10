#include "TH1F.h"
#include "TCut.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TString.h"
#include "mycode/cuts.hh"
#include <fstream>
#include <iostream>
using std::cout;
using std::endl;

void PlotR22R24(){
  TChain r22("ntp1");
  TChain r24("ntp1");
  TChain r26("ntp1");
  r22.Add("AWG82/ntuples/small/Gen_AddRests_RunAll.root");
  r24.Add("AWG82/ntuples/small/R22Add_R24Rests_RunAll.root");
  r26.Add("AWG82/ntuples/small/Add_R26_Run1.root");
  r24.SetLineColor(2);
  r26.SetLineColor(4);
    double NB0Run[2][6] = {{36968000, 103498000, 50556000, 167994000, 244322000,  68148000},
  			 {34878000, 101690000, 56035000, 166784000, 215168000, 130336000}};
   double NBpRun[2][6] = {{37200000, 103124000, 49766000, 167466000, 244812000,  68016000},
 			 {34941000, 104188000, 57888000, 169801000, 215953000, 135224000}};
  double Nevents[3] = {0,0,0};
  for(int i=0; i<2; i++) for(int j=0; j<6; j++) Nevents[i] += NBpRun[i][j];
   TCut cuts[6] = {"MCType==5","MCType==6","(MCType==1||MCType==3)",
 		  "MCType==2||MCType==4","MCType>6", "MCType==0"};
//  TCut cuts[6] = {"MCType==11","MCType==12","(MCType==7||MCType==9)",
//		  "MCType==8||MCType==10","((MCType>0&&MCType<7)||MCType>12)", "MCType==0"};
  TString limits[6] = {"(15,-1,10)","(15,-1,10)","(50,-0.8,1.5)","(50,-0.8,2.5)","(25,-2,10)","(18,-2,11)"};
  //TString limits[6] = {"(13,-1.2,7.5)","(15,-1,10)","(15,-2,1.5)","(50,-0.8,1.5)","(25,-2,10)","(18,-3,10)"};
  TString titles[6] = {"D^{0} #tau #nu","D*^{0} #tau #nu","D^{0} l #nu","D*^{0} l #nu", "D** & Crossfeed","Combinatoric"};
  //TString titles[6] = {"D^{+} #tau #nu","D*^{+} #tau #nu","D^{+} l #nu","D*^{+} l #nu", "D** & Crossfeed","Combinatoric"};
  TH1F *h[3][6]; TLatex label; TLegend *leg[6];
  label.SetTextSize(0.07);
  TCanvas c("R22_R24","Comparison R22/R24",1700,1800);
  c.Divide(2,3,0.001,0.001);
  gStyle->SetOptStat(0);
  for(int i=0; i<6; i++){
    c.cd(i+1);
    TString var = "candM2>>histo"; 
    TString var24 = "candM2>>histo24"; 
    cuts[i] += "candType==1";
    TCut cutR24 = cuts[i];
    if(i==-3) {
      var = "mm2pi0>>histo"; 
      var24 = "mm2pi0>>histo24"; 
      cuts[i] += dssMvaAll22;
      cutR24 += dssMvaAll;
    } else {
      cuts[i] += MvaAll22;
      cutR24 += MvaAll22;
    }
    var += limits[i];
    var24 += limits[i];
    double n22 = r22.Draw(var,cuts[i]);
    TH1F *hTemp = (TH1F*)gDirectory->Get("histo");
    h[0][i] = (TH1F*)hTemp->Clone("r22_"+i);
    double n24 = r24.Draw(var24,cutR24);
    hTemp = (TH1F*)gDirectory->Get("histo24");
    h[1][i] = (TH1F*)hTemp->Clone("r24_"+i);
//     double n26 = r26.Draw(var,cuts[i]);
//     hTemp = (TH1F*)gDirectory->Get("histo");
//     h[2][i] = (TH1F*)hTemp->Clone("r26_"+i);
    double int0 = h[0][i]->Integral();
    if(int0) h[1][i]->Scale(int0/h[1][i]->Integral());
    //    h[2][i]->Scale(Nevents[0]/Nevents[2]);
    h[0][i]->SetTitle("");
    double maxi = h[0][i]->GetMaximum();
    if(maxi<h[1][i]->GetMaximum()) maxi = h[1][i]->GetMaximum();
    //    if(maxi<h[2][i]->GetMaximum()) maxi = h[2][i]->GetMaximum();
    h[0][i]->SetMaximum(maxi*1.1);
    h[0][i]->SetXTitle("m^{2}_{miss} [Gev^{2}]");
    h[0][i]->Draw();
    h[1][i]->Draw("same");
    //h[2][i]->Draw("same");
    if(i==-4)leg[i] = new TLegend(0.1,.75,0.38,0.9);
    else leg[i] = new TLegend(0.62,.75,0.9,0.9);
    leg[i]->SetFillColor(0);
    TString tag = "R22 ("; tag += (int)n22; tag += ")";
    leg[i]->AddEntry(h[0][i],tag);
    tag = "R24 ("; tag += (int)(n24*Nevents[0]/Nevents[1]); tag += ")";
    leg[i]->AddEntry(h[1][i],tag);
    //    tag = "R26 ("; tag += (int)(n26*Nevents[0]/Nevents[2]); tag += ")";
    //leg[i]->AddEntry(h[2][i],tag);
    leg[i]->Draw();
    label.SetNDC(kTRUE);
    label.DrawLatex(0.13,0.92,titles[i]);
  }
  c.SaveAs("babar_code/R24/R22_R24.eps");
}
