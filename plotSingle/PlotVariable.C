//------------------------------------------------------------------------
// Description:
//      PlotVariable - Plots un variable for the 4 channels in generics
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      15/10/10 manuelf -- Created from mycode/CSampleBF.C
//------------------------------------------------------------------------

#include "TString.h"
#include "TPad.h"
#include "TCut.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TChain.h"
#include "../DonutUtils/cuts.cc"
#include "babar_code/Styles/Styles.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

void PlotVariable(TString Variable="candM2", int nbins=20, float minX=-0.2, float maxX=0.3,
		  TString TagFile = "", TString extraCut="", TString sample = "sig"){

  TChain gen("ntp1"), gen2("ntp1"), gen3("ntp1");
  gen2.Add("AWG82/ntuples/small/RAll_RunAll.root");
  gen.Add("AWG82/ntuples/small/cocktail_RunAll.root");
  gen3.Add("AWG82/ntuples/small/Dpipi_RunAll.root");
//    gen.Add("AWG82/ntuples/small/uds_RunAll.root");
//    gen.Add("AWG82/ntuples/small/ccbar_RunAll.root");

  Styles style4; style4.setPadsStyle(4); style4.applyStyle();
  TString titles[] = {"D^{0}","D*^{0}","D^{+}","D*^{+}"};
  //TString cuts[] = {"(MCType==1||MCType==3)","(MCType==2||MCType==4)","(MCType==7||MCType==9)","(MCType==8||MCType==10)"};
  //TString cuts[] = {"(MCType==5)", "(MCType==6)","(MCType==11)", "(MCType==12)"};
  //TString cuts[] = {"(MCType>0&&MCType<7)","(MCType>0&&MCType<7)","(MCType>6&&MCType<13)","(MCType>6&&MCType<13)"};
  //TString cuts[] = {"(MCType>6&&MCType<13)","(MCType>6&&MCType<13)","(MCType>0&&MCType<7)","(MCType>0&&MCType<7)"};
  TString cuts[] = {"(MCType>13)","(MCType>13)","(MCType>13)","(MCType>13)"};
  TString legName[] = {"D** gen (", "D** #rightarrow D#pi (", "D** #rightarrow D#pi#pi ("};
  TString xtitle = Variable;
  if(Variable=="candM2" || Variable == "mm2pi0") xtitle = "m^{2}_{miss} [GeV^{2}]"; 
  if(Variable=="candPstarLep" ) xtitle = "p^{*}_{l} [GeV]"; 
  if(Variable=="candEExtra" ) xtitle = "E_{Extra} [GeV]"; 
  if(Variable=="candM2-candM2Tru") xtitle = "m^{2}_{miss} - m^{2}_{miss,true} [GeV^{2}]"; 
  TH1F *hVari[4],*hVari2[4],*hVari3[4];
  gStyle->SetOptStat(0); gStyle->SetStatFont(style4.nFont); 
  gStyle->SetStatY(1-style4.PadTopMargin); gStyle->SetStatX(1-style4.PadRightMargin);
  gStyle->SetStatW(0.3); gStyle->SetStatH(0.3); 
  TLegend *leg[4];
  double legW = 0.45, legH = 0.25;
  double legX = 1-style4.PadRightMargin-0.01, legY = 1-style4.PadTopMargin-0.01;
  int entries[3];
  TCanvas c("Canvas",Variable);
  c.Divide(2,2);
  for(int i=0; i<4; i++){
    c.cd(i+1);
    TString candCut = "candType=="; candCut += i+1; candCut += "&&"; candCut += cuts[i];
    if(extraCut != ""){candCut += "&&"; candCut += extraCut;}
    TCut totCut = MvaAll; if(sample.Contains("dss")) totCut = dssMvaAll;
    totCut += candCut;
    TString hname = Variable; hname += i; hname.ReplaceAll("-","_");
    TString vari = Variable; vari += ">>"; vari += hname;
    hVari[i] = new TH1F(hname,"",nbins,minX,maxX);
    entries[1] = gen.Draw(vari,totCut);
    hname = Variable; hname += (i+20); hname.ReplaceAll("-","_");
    vari = Variable; vari += ">>"; vari += hname;
    hVari2[i] = new TH1F(hname,"",nbins,minX,maxX);
    entries[0] = gen2.Draw(vari,totCut);
    hname = Variable; hname += (i+30); hname.ReplaceAll("-","_");
    vari = Variable; vari += ">>"; vari += hname;
    hVari3[i] = new TH1F(hname,"",nbins,minX,maxX);
    entries[2] = gen3.Draw(vari,totCut);
    hVari[i]->SetTitle("");
    hVari[i]->SetMinimum(0);
    hVari[i]->SetLineWidth(2);
    hVari[i]->SetLineColor(4);
    hVari[i]->Draw();
    hVari2[i]->Scale(hVari[i]->Integral()/hVari2[i]->Integral());
    hVari2[i]->SetLineWidth(2);
    hVari2[i]->SetLineColor(1);
    hVari2[i]->Draw("same");
    hVari3[i]->Scale(hVari[i]->Integral()/hVari3[i]->Integral());
    hVari3[i]->SetLineWidth(2);
    hVari3[i]->SetLineColor(2);
    hVari3[i]->Draw("same");
    //titles[i] += "  "; titles[i] += TagTitle;
    style4.setTitles(hVari[i],xtitle,"",titles[i]);
    leg[i] = new TLegend(legX-legW, legY-legH, legX, legY);
    TH1F *hAll[] = {hVari2[i], hVari[i], hVari3[i]};
    for(int j=0; j<3; j++){
      TString Name = legName[j]; Name += entries[j]; Name += ")";
      leg[i]->AddEntry(hAll[j],Name);
    }
    leg[i]->SetTextSize(style4.LabelSize*0.9); leg[i]->SetFillColor(0); leg[i]->SetTextFont(style4.nFont);
    leg[i]->SetBorderSize(0);
    if(!sample.Contains("noLeg")) leg[i]->Draw();
  }
  TString fileName = "babar_code/eps/";  fileName += Variable; fileName += TagFile; fileName += ".eps";
  c.SaveAs(fileName);
  for(int i=0; i<4; i++){
    hVari[i]->Delete();
    hVari2[i]->Delete();
    hVari3[i]->Delete();
  }
}

 
