#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TCut.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

TString RoundNumber(double n, int e, double d=1);
void getNumberB(TString genName, TString runs, double &totMCB, double &totdata, 
		double &totuds, double &totccbar, double &totOffdata); 

void EvtSel_TableCuts(int doDss = 0){

  TString NameTrees[3] = {"AWG82/ntuples/small/RAll_RunAll.root", 
			  "AWG82/ntuples/small/uds_RunAll.root", "AWG82/ntuples/small/ccbar_RunAll.root"};
  TChain gen("ntp1"), cont("ntp1");
  gen.Add(NameTrees[0]);
  for(int t=1; t<3; t++) cont.Add(NameTrees[t]);
  double totMCB = 0, totuds = 0, totccbar = 0, totdata = 0, totOffdata = 0;
  getNumberB(NameTrees[0], "All", totMCB, totdata, totuds, totccbar, totOffdata);
  double wuds = totMCB/totuds*2.09/1.05;     
  double wMC = totdata/totMCB; 

   TCut MvaBp = "(candMvaDl>0.458&&candType==1)||(candMvaDl>0.48&&candType==2)";      // onDx200
   TCut MvaB0 = "(candMvaDl>0.360&&candType==3)||(candMvaDl>0.41&&candType==4)";      // onDx200
//   TCut MvaBp = "(candMvaDl>0.63&&candType==1)||(candMvaDl>0.48&&candType==2)";
//   TCut MvaB0 = "(candMvaDl>0.51&&candType==3)||(candMvaDl>0.41&&candType==4)";
  TCut Mva = MvaBp||MvaB0;
  TCut other = "candMES>5.27&&candPMiss>0.2&&candQ2>4&&candM2>1.5";
  TCut ee = "(candType<3&&candEExtra<.2||candType==3&&candEExtra<.15||candType==4&&candEExtra<.3)";

  TCut dssacc = "candMES>5.27&&mm2pi0>-4&&mm2pi0<12&&candPstarLep>0&&candPstarLep<2.4&&candMES>5.2&&candMES<5.3";
  TCut mpi0 = "mpi0>.125&&mpi0<.145";
  TCut Q2Pmiss = "candQ2>4&&pmisspi0>.2";
  TCut Mpi0 = "mpi0>.12&&mpi0<.15";
  TCut mm2pi0 = "mm2pi0>-0.5&&mm2pi0<1.5";
  TCut dssee = "eextrapi0<.5&&ppi0>.4";
  TCut cosT = "abs(candCosT)<.8";
  TCut dssMvaDl = "candMvaDssDl>-0.45";
  TCut dssMvaComb="candMvaDssComb>-0.35&&candType<3||candMvaDssComb>-0.3&&candType==3||candMvaDssComb>-0.35&&candType==4";
  TCut dssMva = dssMvaDl+dssMvaComb;
  TCut dsseeAll = dssacc+mpi0+dssee+mm2pi0;

  //TCut rowCuts[2][7] = {{"candM2>1.5", "candPMiss>0.2&&candQ2>4", Mva, "candM2>1.5",other+ee,"",""},
  TCut rowCuts[2][7] = {{"candMES>5.27&&candM2>1.5", "candPMiss>0.2&&candQ2>4", Mva, "candM2>1.5",other+ee,"",""},
			{dssacc, Q2Pmiss, Mpi0, cosT, dssMva, mm2pi0, dsseeAll}};
  TCut colCuts[2][5] = {{"(candType<3&&MCType>4&&MCType<7||candType>2&&MCType>10&&MCType<13)",
			 "(candType<3&&MCType>0&&MCType<5||candType>2&&MCType>6&&MCType<11)",
			 "MCType>12", "(candType>2&&MCType>0&&MCType<7||candType<3&&MCType>6&&MCType<13)",
			 "MCType==0"},
			{"MCType>12",
			 "(candType<3&&MCType>0&&MCType<7||candType>2&&MCType>6&&MCType<13)",
			 "(candType>2&&MCType>0&&MCType<7||candType<3&&MCType>6&&MCType<13)",
			 "MCType==0","MCType<0"}};
  TString labels[2][7] = {{"$\\Btag\\ds\\ell$","Kin.","BDT",
			  "$\\mmiss>1.5\\text{ GeV}^2$", "Previous","",""},
			  {"$\\Btag\\ds\\pi^0\\ell$","Kin.", 
			   "$120<m_{\\pi^0}<150\\text{ MeV}$","$|\\text{cos}\\Delta\\theta_T|<0.8$",
			   "BDT", "$-0.5<\\mmiss<1.5\\text{ GeV}^2$", "Previous"}};

  fstream TextFile;
  TextFile.open("public_html/EvtSel_TableCuts.txt",fstream::out);
  TextFile<<"Signal sample"<<endl<<endl;
  TH1F *hCount = new TH1F("hCount","",100,-4,12);
  int nRows = 5, nCols = 5; 
  if(doDss) {nRows = 7; nCols = 4;}
  for(int row=0; row<nRows; row++) {
    TCut rowCut = rowCuts[doDss][0]; 
    if(row < (nRows-1)) for(int cut=1; cut<=row; cut++) rowCut += rowCuts[doDss][cut]; 
    else rowCut = rowCuts[doDss][nRows-1];
    //rowCut += "candType==4";
    TCut WrowCut = rowCut; WrowCut *= "weight";
    gen.Draw("candM2>>hCount",WrowCut);
    double nTotal = hCount->Integral(), nSignal = 0;
    cont.Draw("candM2>>hCount",WrowCut);
    double nCont = hCount->Integral()*wuds;
    nTotal += nCont;
    TextFile<<labels[doDss][row]<<"\t & ";
    for(int col=0; col < nCols; col++) {
      TCut totCut = rowCut; totCut += colCuts[doDss][col]; totCut *= "weight";
      gen.Draw("candM2>>hCount",totCut);
      double n = hCount->Integral();
      if(col==0) nSignal = n;
      TextFile<<RoundNumber(n*100,1,nTotal)<<" & ";
    }
    TextFile<<RoundNumber(nCont*100,1,nTotal)<<" & "<<RoundNumber(nSignal*wMC,0)<<" & ";
    TextFile<<RoundNumber((nTotal-nSignal)*wMC,0)<<" \\\\"<<endl;
    if(row==nRows-2) TextFile<<"\\hline"<<endl;
  }
  doDss = 1; nRows = 7; nCols = 4;
  TextFile<<endl<<endl<<"Dpi0 sample"<<endl<<endl;
  for(int row=0; row<nRows; row++) {
    TCut rowCut = rowCuts[doDss][0]; 
    if(row < (nRows-1)) for(int cut=1; cut<=row; cut++) rowCut += rowCuts[doDss][cut]; 
    else rowCut = rowCuts[doDss][nRows-1];
    //rowCut += "candType==4";
    TCut WrowCut = rowCut; WrowCut *= "weight";
    gen.Draw("candM2>>hCount",WrowCut);
    double nTotal = hCount->Integral(), nSignal = 0;
    cont.Draw("candM2>>hCount",WrowCut);
    double nCont = hCount->Integral()*wuds;
    nTotal += nCont;
    TextFile<<labels[doDss][row]<<"\t & ";
    for(int col=0; col < nCols; col++) {
      TCut totCut = rowCut; totCut += colCuts[doDss][col]; totCut *= "weight";
      gen.Draw("candM2>>hCount",totCut);
      double n = hCount->Integral();
      if(col==0) nSignal = n;
      TextFile<<RoundNumber(n*100,1,nTotal)<<" & ";
    }
    TextFile<<RoundNumber(nCont*100,1,nTotal)<<" & "<<RoundNumber(nSignal*wMC,0)<<" & ";
    TextFile<<RoundNumber((nTotal-nSignal)*wMC,0)<<" \\\\"<<endl;
    if(row==nRows-2) TextFile<<"\\hline"<<endl;
  }
  hCount->Delete();
}

void getNumberB(TString genName, TString runs, double &totMCB, double &totdata, 
		double &totuds, double &totccbar, double &totOffdata){
  double NBR24[2][6] = {{  34878000, 101690000,  56035000, 166784000, 215168000, 130336000},
			{  34941000, 104188000,  57888000, 169801000, 215953000, 135224000}};
  double NBR26[2][6] = {{  78999000, 235345000, 120771000, 389670000, 509712000, 290169000},
			{  78560000, 246656000, 122374000, 383657000, 526675000, 291862000}};
  double nccbar[] =     {  55254000, 164146000,  88321000, 267308000, 343667000, 208664000};
  double nuds[] =       { 160514000, 451636000, 275869000, 421599000, 553604000, 327032000};
  double ndata[] =      { 22556256.9, 68438426.0, 35763257.9, 111429669.4, 147620363.4, 85194672.2};
  double lumOffdata[] = {2621.575, 7030.748, 2495.880, 10228.272, 14546.087, 7887.3};
  double fracRest24 = 2/3., fracRest26 = 1.;
  totMCB = 0; totuds = 0; totccbar = 0; totdata = 0; totOffdata = 0;
  if(genName.Contains("R24Rest")) fracRest24 = 1/3.;
  if(genName.Contains("R26Rest")) fracRest26 = 1/7.;
  for(int i=0; i<6; i++){
    TString Run = ""; Run += i+1;
    if(runs=="All" || runs.Contains(Run)){
      totdata += ndata[i];			                     //  471.003M events
      if(genName.Contains("RAll") || genName.Contains("R24")) 
	totMCB += (NBR24[0][i]+NBR24[1][i])*fracRest24;              //  948.591M events
      if(genName.Contains("RAll") || genName.Contains("R26")) 
	totMCB += (NBR26[0][i]+NBR26[1][i])*fracRest26;              // 3274.450M events
    }
    totOffdata += lumOffdata[i];                                     //   44.810k pb
    totuds += nuds[i];				                     // 2190.254M events
    totccbar += nccbar[i];			                     // 1127.360M events
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
