#define Ntuple_cxx
#include "Ntuple.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

TString round(double n, double d);

void Ntuple::Loop(bool doCut) {
   if (fChain == 0) return;

   TString MMiss = "candPMiss>.2";
   TString Q2 = "candQ2>4";
   TString ee1 = "candType<3&&candEExtra<.2";
   TString ee2 = "candType==3&&candEExtra<.15";
   TString ee3 = "candType==4&&candEExtra<.3";
   TString ee = (ee1||ee2||ee3);
   TString basic = MMiss+Q2+ee;
   int type[90][59][6];
   int SeedMode, YMode;
   fstream table;
   if(doCut)
     table.open("modes/txt/TableRun"+Run+"_truthmatched_cuts.txt",fstream::out);
   else
     table.open("modes/txt/TableRun"+Run+"_truthmatched.txt",fstream::out);
   int mode;
   for(int i=0; i<87; i++){
     if(i==2) i = 10; if(i==19) i = 20; if(i==27) i = 30;
     if(i==49) i = 50; if(i==59) i = 60; if(i==62) i = 71;
     if(i==73) i = 80; if(i==85) i = 86;
     for(int j=1; j<54; j++){
       for(int k=0; k<6; k++)
	 type[i][j][k] = 0;
     }
   }

   Long64_t nbytes = 0, nb = 0;
   Long64_t nentries = fChain->GetEntriesFast();
   //nentries = 10000;
   cout <<"Entries are "<<nentries<<endl;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (jentry%10000 == 0) cout<<"Entry "<<jentry<<" of "<<nentries<<endl;

      if(doCut){
	if(candPMiss<=.2) continue;
	if(candQ2<=4) continue;
	if(candType<3&&candEExtra>=.2 || candType==3&&candEExtra>=.15 || candType==4&&candEExtra>=.3) continue;
      }
      if(!candLepTru) continue;
      SeedMode = (candBMode/100)%100;
      YMode = candBMode%100;
      //0 is track truthmatched, 1 is track average, 2 is neutral truthmatched, 3 is track average, 4 is total
      if(candBntCha==0&&candDntCha==0) type[SeedMode][YMode][0]++;
      if(candBntNeu==0) type[SeedMode][YMode][2]++;
      type[SeedMode][YMode][1] += candBntCha + candDntCha;
      type[SeedMode][YMode][3] += candBntNeu;
      type[SeedMode][YMode][4]++;
   }
//    table<<"D   Y          Trk     Tpuri    Tave    Neu    Npuri    Nave"<<endl;
//    table<<"============================================================"<<endl;
   for(int i=0; i<87; i++){
     if(i==2) i = 10; if(i==19) i = 20; if(i==27) i = 30;
     if(i==49) i = 50; if(i==59) i = 60; if(i==62) i = 71;
     if(i==73) i = 80; if(i==85) i = 86;
     for(int j=1; j<54; j++){
       int goodT = type[i][j][0];
       int goodN = type[i][j][2];
       int nevents = type[i][j][4];
       if(!goodT) continue;
       table<<i<<"   "<<j<<"\t|\t"<<goodT<<"\t"<<round((double)goodT,(double)nevents)<<"\t";
       table<<round((double)type[i][j][1],(double)nevents);
       table<<"\t"<<goodN<<"\t"<<round((double)goodN,(double)nevents)<<"\t"<<round((double)type[i][j][3],(double)nevents);
       table<<endl;
     }
       table<<endl;
   }
   table.close();
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
