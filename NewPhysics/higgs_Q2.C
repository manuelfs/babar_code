#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "babar_code/Styles/Styles.cc"
#include <fstream>
#include <iostream>

#define nHis 4
using namespace std;
using std::cout;
using std::endl;

TString RoundNumber(double n, int e, double d=1);

void higgs_Tl(TString tBmH1 = "000", TString tBmH2 = "050", TString tBmH3 = "100"){

  Styles style; style.setPadsStyle(2); style.applyStyle();

  double legW = 0.2, legH = 0.27;
  double legX = style.PadLeftMargin+0.27, legY = 1-style.PadTopMargin-0.52;
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(style.LabelSize); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont); leg.SetBorderSize(0);

  double PI = 3.1415927;
  TChain *tree[nHis];
  TString higName[] = {"000", tBmH1, tBmH2, tBmH3}, fileName, legName[nHis], xTitle = "#theta_{l}";
  TString yTitle[] = {"Exp. D#tau#nu yield in data", "Exp. D*#tau#nu yield in data"};
  TCanvas can("can","tL comparison"); can.Divide(2,1);
  TH1F *Histo[2][nHis];
  double limX[] = {0,PI}, maxH[] = {0,0};
  int nBins = 25, colors[] = {2,4,28,1};
  TString Cuts[] = {"weight*(MCType==5||MCType==11)", "weight*(MCType==6||MCType==12)"};

  for(int his=0; his<nHis; his++){
    fileName = "AWG82/ntuples/small/FitRAllHigx"; fileName += higName[his];
    fileName += "_RunAll.root";
    tree[his] = new TChain("ntp1");
    tree[his]->Add(fileName);
    for(int isDs=0; isDs<=1; isDs++){
      fileName = "Histo"; fileName += isDs; fileName += his;
      Histo[isDs][his] = new TH1F(fileName,"",nBins,limX[0],limX[1]);
      Histo[isDs][his]->SetLineWidth(2);
      Histo[isDs][his]->SetLineColor(colors[his]);
      TString vari = "trueCTL";
      if(isDs) vari = "acos(trueCTL)";
      TString finalCuts = Cuts[isDs];
      if(his==1) {
	if(isDs) finalCuts = "weight*(MCType==4||MCType==10)";
	else finalCuts = "weight*(MCType==3||MCType==9)";
      }
      tree[his]->Project(fileName,vari,finalCuts);

      legName[his] = "t#beta/m_{H} = "; legName[his] += RoundNumber(higName[his].Atof(),1,100.);
      if(fabs(higName[his].Atof())<1e-6) legName[his] = "SM"; 
      if(his==1) {
	legName[his] += ", #mu"; 
	Histo[isDs][his]->Scale(Histo[isDs][his-1]->Integral()/Histo[isDs][his]->Integral());
      }
      if(isDs==1) leg.AddEntry(Histo[isDs][his], legName[his]);
      if(maxH[isDs]<Histo[isDs][his]->GetMaximum()) maxH[isDs]=Histo[isDs][his]->GetMaximum();
    }
  }

  for(int isDs=0; isDs<=1; isDs++){
    can.cd(isDs+1);
    for(int histo=0; histo<nHis; histo++){
      if(histo==0){
	Histo[isDs][histo]->SetMinimum(0); Histo[isDs][histo]->SetMaximum(maxH[isDs]*1.05);
	Histo[isDs][histo]->Draw();
	style.setTitles(Histo[isDs][histo],xTitle,yTitle[isDs]);
      } else Histo[isDs][histo]->Draw("same");
    }
    if(isDs==1) leg.Draw();
  }
      

  TString pName = "public_html/ExpData_Tl.eps"; 
  can.SaveAs(pName);
  for(int his=0; his<nHis; his++){
    tree[his]->Delete();
    for(int isDs=0; isDs<=1; isDs++) Histo[isDs][his]->Delete();
  }

}

void higgs_Q2(TString tBmH1 = "030", TString tBmH2 = "050", TString tBmH3 = "100"){

  Styles style; style.setPadsStyle(2); style.applyStyle();

  double legW = 0.2, legH = 0.27;
  double legX = style.PadLeftMargin+0.27, legY = 1-style.PadTopMargin-0.52;
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(style.LabelSize); leg.SetFillColor(0); 
  leg.SetTextFont(style.nFont); leg.SetBorderSize(0);

  TChain *tree[nHis];
  TString higName[] = {"000", tBmH1, tBmH2, tBmH3}, fileName, legName[nHis], xTitle = "q^{2} (GeV^{2})";
  TString yTitle[] = {"Exp. D#tau#nu yield in data", "Exp. D*#tau#nu yield in data"};
  TCanvas can("can","q2 comparison"); can.Divide(2,1);
  TH1F *Histo[2][nHis];
  double limX[] = {4,12}, maxH[] = {0,0};
  int nBins = 32, colors[] = {2,4,28,1};
  TString Cuts[] = {"weight*(MCType==5||MCType==11)", "weight*(MCType==6||MCType==12)"};

  for(int his=0; his<nHis; his++){
    fileName = "AWG82/ntuples/small/FitRAllHigx"; fileName += higName[his];
    fileName += "_RunAll.root";
    tree[his] = new TChain("ntp1");
    tree[his]->Add(fileName);
    for(int isDs=0; isDs<=1; isDs++){
      fileName = "Histo"; fileName += isDs; fileName += his;
      Histo[isDs][his] = new TH1F(fileName,"",nBins,limX[0],limX[1]);
      Histo[isDs][his]->SetLineWidth(2);
      Histo[isDs][his]->SetLineColor(colors[his]);
      tree[his]->Project(fileName,"trueQ2",Cuts[isDs]);

      legName[his] = "t#beta/m_{H} = "; legName[his] += RoundNumber(higName[his].Atof(),1,100.);
      if(fabs(higName[his].Atof())<1e-6) legName[his] = "SM"; 
      if(isDs==1) leg.AddEntry(Histo[isDs][his], legName[his]);
      if(maxH[isDs]<Histo[isDs][his]->GetMaximum()) maxH[isDs]=Histo[isDs][his]->GetMaximum();
    }
  }

  for(int isDs=0; isDs<=1; isDs++){
    can.cd(isDs+1);
    for(int histo=0; histo<nHis; histo++){
      if(histo==0){
	Histo[isDs][histo]->SetMinimum(0); Histo[isDs][histo]->SetMaximum(maxH[isDs]*1.05);
	Histo[isDs][histo]->Draw();
	style.setTitles(Histo[isDs][histo],xTitle,yTitle[isDs]);
      } else Histo[isDs][histo]->Draw("same");
    }
    if(isDs==1) leg.Draw();
  }
      

  TString pName = "public_html/ExpData_Q2.eps"; 
  can.SaveAs(pName);
  for(int his=0; his<nHis; his++){
    tree[his]->Delete();
    for(int isDs=0; isDs<=1; isDs++) Histo[isDs][his]->Delete();
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
