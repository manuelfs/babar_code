#include "../DonutUtils/cuts.cc"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TCut.h"
#include "TBranch.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;
using std::cout;
using std::endl;

struct Pair {
  int a;
  TString b;
};
long expo(long x, int n);
int compare (const void *const first, const void *const second);
TString round(double n, int e, double d=1.);

void CombAbundance(){
  TString files = "AWG82/ntuples/small/RAll_RunAll.root";
  TChain c("ntp1");
  c.Add(files);
  TCut cut = "candQ2>4&&candEExtra>1.2&&candEExtra<2.4";
  cut = MvaAll; cut += "candM2>5&&(candType%2==0)";
  TTree *t = c.CopyTree(cut);
  //TTree *t = c.CopyTree(cut);
  cout<<"Using "<<files<<endl;

  int MCType, MCComblong, MCComblongD, isBzero;

  TBranch *b_MCType;
  TBranch *b_MCComblong;
  TBranch *b_MCComblongD;
  TBranch *b_isBzero;

  t->SetBranchStatus("*",0);
  t->SetBranchStatus("MCType",1);
  t->SetBranchStatus("MCComblong",1);
  t->SetBranchStatus("MCComblongD",1);
  t->SetBranchStatus("isBzero",1);
  t->SetBranchAddress("MCType",&MCType,&b_MCType);
  t->SetBranchAddress("MCComblong",&MCComblong,&b_MCComblong);
  t->SetBranchAddress("MCComblongD",&MCComblongD,&b_MCComblongD);
  t->SetBranchAddress("isBzero",&isBzero,&b_isBzero);

  TString taFile = "babar_code/Reweight/Tables/NewCombinatoric_Distribution.txt"; 
  fstream Table;
  Table.open(taFile,fstream::out);
  long nentries = t->GetEntries();
  //nentries=100;
  int n[24][2], decay2[24][24][2], decay3[24][24][24][2], decay4[24][24][24][24][2];
  int nBComb[]={0,0}, B0, nbody[20];
  for(int i=0; i<24; i++)
    for(int j=0; j<24; j++)
      for(int k=0; k<2; k++){
	decay2[i][j][k] = 0;
	for(int m=0; m<24; m++){
	  decay3[i][j][m][k] = 0;
	  for(int p=0; p<24; p++)
	    decay4[i][j][m][p][k] = 0;
	}
      }
  for(int i=0; i<20; i++) nbody[i]=0;

  TString names[] = {"$D^0$","$D^{\\pm}$","$D^{*0}$","$D^{*\\pm}$","$D^{**0}$","$D^{**\\pm}$","$D_s^{\\pm}$",
		     "$D_s^{*\\pm}$","$D_s^{**\\pm}$","$e$","$\\mu$","$\\tau$","$\\pi^{\\pm}$","$\\pi^0$",
		     "$K^{\\pm}$","$K^0$","$\\rho^{\\pm}$","$\\rho^0$",
		     "$K^{*\\pm}$","$K^{*0}$","$a_1^{\\pm}$","$a_1^0$","$\\omega$","$\\gamma$"};
  for (int j = 0 ; j < nentries ; j ++) {
    t->LoadTree(j);
    t->GetEntry(j);
    if (j%40000 == 0) cout<<"Entry "<<j<<" of "<<nentries<<endl;

    if(MCType!=0) continue;
    if(isBzero) B0=1;
    else B0 = 0;
    nBComb[B0]++;
    if(MCComblong>0 || MCComblongD>0) {
      int total=0;
      for(int i=0; i<9; i++){
	long mask = expo(2,i*2)+expo(2,i*2+1);
	long num = expo(2,i*2);
	n[i][B0] = (MCComblongD & mask)/num;
	total += n[i][B0];
      }
      
      for(int i=0; i<15; i++){
	long mask = expo(2,i*2)+expo(2,i*2+1);
	long num = expo(2,i*2);
	n[i+9][B0] = (MCComblong & mask)/num;
	total += n[i+9][B0];
      }
      nbody[total]++;
      if(total==2){
	int first = -1;
	for(int i=0; i<24; i++){
	  if(n[i][B0]==1){
	    if(first==-1) first = i;
	    else {
	      decay2[first][i][B0]++;
	      break;
	    }
	  }else if(n[i][B0]==2){
	    decay2[i][i][B0]++;
	    break;
	  }
	}
      }
      if(total==3){
	int first = -1, second = -1;
	for(int i=0; i<24; i++){
	  if(n[i][B0]==1){
	    if(first==-1) first = i;
	    else if (second==-1) second = i;
	    else {
	      decay3[first][second][i][B0]++;
	      break;
	    }
	  }else if(n[i][B0]==2){
	    if(first==-1) {
	      first = i;
	      second = i;
	    }else{
	      decay3[first][i][i][B0]++;
	      break;
	    }
	  } else if(n[i][B0]==3){
	      decay3[i][i][i][B0]++;
	      break;
	  }
	}
      }
      if(total==4){
	int first = -1, second = -1, third = -1;
	for(int i=0; i<24; i++){
	  if(n[i][B0]==1){
	    if(first==-1) first = i;
	    else if (second==-1) second = i;
	    else if (third==-1) third = i;
	    else {
	      decay4[first][second][third][i][B0]++;
	      break;
	    }
	  }else if(n[i][B0]==2){
	    if(first==-1) {
	      first = i;
	      second = i;
	    }else if (second==-1){
	      second = i;
	      third = i;
	    }else {
	      decay4[first][second][i][i][B0]++;
	      break;
	    }
	  } else if(n[i][B0]==3){
	    if(first==-1) {
	      first = i;
	      second = i;
	      third = i;
	    }else {
	      decay4[first][i][i][i][B0]++;
	      break;
	    } 
	  } else if(n[i][B0]==4){
	      decay4[i][i][i][i][B0]++;
	      break;
	  }
	}
      }
    } // if it is a tagged combinatoric event
  } // end of entries for

  Table<<"\\begin{landscape} \\begin{tabular}{|l||cccc|cc|ccc|ccc|cccc|cccc|cccc|}"<<endl<<"\\hline"<<endl;
  for(int i=0; i<24; i++)
    Table<<" & "<<names[i];
  Table<<" \\\\ \\hline \\hline"<<endl;
  for(int i=0; i<24; i++){
    Table<<names[i]<<" ";
    for(int j=0; j<24; j++){
      Table<<" & ";
      if(j>=i) Table<<decay2[i][j][0]+decay2[i][j][1];
    }
    Table<<" \\\\"<<endl;
    if(i==3||i==5||i==8||i==11||i==15||i==19) Table<<"\\hline "<<endl;
  }
  Table<<"\\hline \\end{tabular} \\end{landscape}"<<endl<<endl;

  Pair list2Bp[24*24],list2B0[24*24],list3Bp[24*24*24],list3B0[24*24*24];
  Pair listBp[24*24*24+24*24],listB0[24*24*24+24*24];
  int ntot = 0, ntotB0 = 0, nmax = 24*24*24+24*24;
  for(int i=0; i<24; i++){
    for(int j=0; j<24; j++){
      list2Bp[i+j*24].a = -decay2[i][j][0];
      list2Bp[i+j*24].b = names[i]; list2Bp[i+j*24].b += " "; list2Bp[i+j*24].b += names[j];
      list2B0[i+j*24].a = -decay2[i][j][1];
      list2B0[i+j*24].b = names[i]; list2B0[i+j*24].b += " "; list2B0[i+j*24].b += names[j];

      listBp[ntot].a = -decay2[i][j][0];
      listBp[ntot].b = names[i]; listBp[ntot].b += " "; listBp[ntot].b += names[j];
      ntot++;
      listB0[ntotB0].a = -decay2[i][j][1];
      listB0[ntotB0].b = names[i]; listB0[ntotB0].b += " "; listB0[ntotB0].b += names[j];
      ntotB0++;
      for(int k=0; k<24; k++){
	list3Bp[i+j*24+k*24*24].a = -decay3[i][j][k][0];
	list3Bp[i+j*24+k*24*24].b = names[i]; list3Bp[i+j*24+k*24*24].b += " "; list3Bp[i+j*24+k*24*24].b += names[j];
	list3Bp[i+j*24+k*24*24].b += " "; list3Bp[i+j*24+k*24*24].b += names[k];
	list3B0[i+j*24+k*24*24].a = -decay3[i][j][k][1];
	list3B0[i+j*24+k*24*24].b = names[i]; list3B0[i+j*24+k*24*24].b += " "; list3B0[i+j*24+k*24*24].b += names[j];
	list3B0[i+j*24+k*24*24].b += " "; list3B0[i+j*24+k*24*24].b += names[k];

	if(decay3[i][j][k][0]>0){
	  listBp[ntot].a = -decay3[i][j][k][0];
	  listBp[ntot].b = names[i]; listBp[ntot].b += " "; listBp[ntot].b += names[j];
	  listBp[ntot].b += " "; listBp[ntot].b += names[k];
	  ntot++;
	}
	if(decay3[i][j][k][1]>0){
	  listB0[ntotB0].a = -decay3[i][j][k][1];
	  listB0[ntotB0].b = names[i]; listB0[ntotB0].b += " "; listB0[ntotB0].b += names[j];
	  listB0[ntotB0].b += " "; listB0[ntotB0].b += names[k];
	  ntotB0++;
	}
	for(int m=0; m<24; m++){
	  if(decay4[i][j][k][m][0]>0){
	    listBp[ntot].a = -decay4[i][j][k][m][0];
	    listBp[ntot].b = names[i]; listBp[ntot].b += " "; listBp[ntot].b += names[j];
	    listBp[ntot].b += " "; listBp[ntot].b += names[k]; listBp[ntot].b += " "; listBp[ntot].b += names[m];
	    ntot++;
	  }
	  if(decay4[i][j][k][m][1]>0){
	    listB0[ntotB0].a = -decay4[i][j][k][m][1];
	    listB0[ntotB0].b = names[i]; listB0[ntotB0].b += " "; listB0[ntotB0].b += names[j];
	    listB0[ntotB0].b += " "; listB0[ntotB0].b += names[k]; listB0[ntotB0].b += " "; listB0[ntotB0].b += names[m];
	    ntotB0++;
	  }
	  if(ntot>=nmax || ntotB0>=nmax) {cout<<"Too many non-zero entries"<<endl; break;}
	}
	if(ntot>=nmax || ntotB0>=nmax) {cout<<"Too many non-zero entries"<<endl; break;}
      }
      if(ntot>=nmax || ntotB0>=nmax) {cout<<"Too many non-zero entries"<<endl; break;}
    }
  }
  cout<<"Filled "<<ntot<<" out of "<<nmax<<" elements"<<endl;

  qsort(list2Bp, sizeof list2Bp / sizeof *list2Bp, sizeof *list2Bp, compare);
  qsort(list2B0, sizeof list2B0 / sizeof *list2B0, sizeof *list2B0, compare);
  qsort(list3Bp, sizeof list3Bp / sizeof *list3Bp, sizeof *list3Bp, compare);
  qsort(list3B0, sizeof list3B0 / sizeof *list3B0, sizeof *list3B0, compare);
  qsort(listBp, ntot, sizeof *listBp, compare);
  qsort(listB0, ntotB0, sizeof *listB0, compare);

  Table<<"\\begin{tabular}{rlr|lr}"<<endl<<"\\multicolumn{3}{c|}{$\\boldsymbol{B^+("<<nBComb[0]<<")}$} & ";
  Table<<"\\multicolumn{2}{c}{$\\boldsymbol{B^0("<<nBComb[1]<<")}$}"<<" \\\\"<<endl;
  Table<<"& Decay & \\% & Decay & \\% \\\\ \\hline \\hline"<<endl;
  double B0done = 0., Bpdone = 0;
  for(int i=0; i<60; i++){
    Table<<i+1<<" & "<<listBp[i].b<<" & "<<round((double)(-listBp[i].a)/nBComb[0]*100.,2)<<" & ";
    Table<<listB0[i].b<<" & "<<round((double)(-listB0[i].a)/nBComb[1]*100.,2)<<" \\\\"<<endl;
    Bpdone -= listBp[i].a; B0done -= listB0[i].a;
  }
  Table<<"\\hline & Total done & "<<round(Bpdone/nBComb[0]*100,2)<<" & Total done & ";
  Table<<round(B0done/nBComb[1]*100,2)<<"\\\\"<<endl;
  Table<<"\\end{tabular}"<<endl<<endl;

  for(int i=0; i<20; i++) cout<<nbody[i]<<", ";
  cout<<endl<<"Combinatoric distribution stored in "<<taFile<<endl;
 
  
}
/*  sum = nel+nmu*expo(4,1)+ntau*expo(4,2)+npi*expo(4,3)+npi0*expo(4,4)+nk*expo(4,5)+nk0*expo(4,6)
    +nrho*expo(4,7)+nrho0*expo(4,8)+nkstar*expo(4,9)+nkstar0*expo(4,10)+na1*expo(4,11)
    +na10*expo(4,12)+nomega*expo(4,13)+ngamma*expo(4,14);

  MCComblongD = nd0+ndp*expo(4,1)+nds0*expo(4,2)+ndsp*expo(4,3)+ndss0*expo(4,4)+ndssp*expo(4,5)
    +nds*expo(4,6)+ndss*expo(4,7)+ndsss*expo(4,8);*/

long expo(long x, int n){
  long base = 1;
  for(long i=0; i<n; i++)
    base *= x;
  return base;
}
int compare (const void *const first, const void *const second){
  if (((const struct Pair *)first)->a > ((const struct Pair *)second)->a)
    return 1;
  else if (((const struct Pair *)first)->a < ((const struct Pair *)second)->a)
    return -1;
  else
    return 0;
}
TString round(double n, int e, double d){
  if(d==0) return " - ";
  double b = (int)(n/d*expo(10,e)+0.5);
  b /= expo(10,e);
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(!result.Contains(".")) result += ".";
  
  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<e-afterdot.Length(); i++)
    result += "0";
  return result;
}


