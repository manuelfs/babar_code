//---------------------------------------------------------------------------------
// Description:
//      Second stage of the mode optimization
//      Macro that calculates "babar_code/txt/yield"+Run+".txt" with the yields
//      of signal and background for each BSemiExclAdd mode
//      I typically cut on EExtra/MVA when making the table
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      09/05/15 manuelf -- Adaptation and clean up from Ntuple.C
//---------------------------------------------------------------------------------

#define yieldDssList_cxx
#include "yieldDssList.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

TString round(double n, double d);

void yieldDssList::Loop(bool doCut) {
   if (fChain == 0) return;

   int type[90][59][6];
   int SeedMode, YMode;
   fstream table;
   TString filename = "babar_code/txt/yieldDss"; filename += Run;
   if(doCut)
     filename += "_cuts.txt";
   else
     filename += ".txt";
   table.open(filename,fstream::out);

   int mode;
   for(int i=0; i<87; i++){
     if(i==2) i = 10; if(i==19) i = 20; if(i==27) i = 30;
     if(i==49) i = 50; if(i==59) i = 60; if(i==62) i = 71;
     if(i==73) i = 80; if(i==85) i = 86;
     for(int j=1; j<54; j++){
       for(int k=0; k<5; k++)
	 type[i][j][k] = 0;
     }
   }

   Long64_t nbytes = 0, nb = 0;
   Long64_t nentries = fChain->GetEntries();
   //nentries = 10000;
   cout <<"Entries are "<<nentries<<endl;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (jentry%10000 == 0) cout<<"Entry "<<jentry<<" of "<<nentries<<endl;

      if(doCut){
	int isBestPi0 = -1;
	for(int i=0; i<npi0; i++) {
	  if(bestepi0[i] == 1) isBestPi0 = i;
	}
	if(isBestPi0==-1) continue;
	if(mpi0[isBestPi0]<=0.125 || mpi0[isBestPi0]>=0.145) continue;
	if(ppi0[isBestPi0]<=0.4) continue;
	if(eextrapi0[isBestPi0]>=.5) continue;
	if(pmisspi0[isBestPi0]<=0.2) continue;
	if(mm2pi0[isBestPi0]<=-4 || mm2pi0[isBestPi0]>=12) continue;
      }
      SeedMode = (candBMode/100)%100;
      YMode = candBMode%100;
      if(MCType==0 || candLepTru==0){
	type[SeedMode][YMode][4]++;
      }else{
	if(MCType>12)
	  type[SeedMode][YMode][0]++;
	if(candType<3){
	  if(MCType>6&&MCType<13)
	    type[SeedMode][YMode][3]++;
	  if(MCType>0&&MCType<5)
	    type[SeedMode][YMode][1]++;
	  if(MCType>4&&MCType<7)
	    type[SeedMode][YMode][1]++;
	}else{
	  if(MCType>0&&MCType<7)
	    type[SeedMode][YMode][3]++;
	  if(MCType>6&&MCType<11)
	    type[SeedMode][YMode][1]++;
	  if(MCType>10&&MCType<13)
	    type[SeedMode][YMode][1]++;
	}
      }
   }
   for(int i=0; i<87; i++){
     if(i==2) i = 10; if(i==19) i = 20; if(i==27) i = 30;
     if(i==49) i = 50; if(i==59) i = 60;if(i==62) i = 71;
     if(i==73) i = 80; if(i==85) i = 86;
     for(int j=1; j<54; j++){
       int B = type[i][j][2]+type[i][j][3]+type[i][j][4];
       int S = type[i][j][0];
       if(!(S+type[i][j][1]+B)) continue;
       table<<i<<"   "<<j<<"\t|\t"<<S<<"\t"<<type[i][j][1]<<"\t"<<type[i][j][2];
       table<<"\t"<<type[i][j][3]<<"\t"<<type[i][j][4]<<"\t"<<round((double)S,(double)(S+B));
       table<<endl;
     }
     table<<endl;
   }

   table.close();
   cout<<"Yields for each mode stored in "<<filename<<endl;
}


TString round(double n, double d){
  if(d==0) return " - ";
  double b = ((int)(n/d *100+.5));
  b/=100.;
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(b<10){
    if(result.Length() == 1)
      result += ".00";
    if(result.Length() == 3)
      result += "0";
  }
  return result;
}
