#include "TString.h"
#include "TSystem.h"
#include <fstream>
#include <iostream>

using namespace std; 
using std::cout;
using std::endl;

void initialize(double smoothings[][7], int counts[][7]);
bool areEqual(double a, double b);
TString RoundNumber(double n, int e, double d=1);

void parseCrossV(const char* folder, int Fit_Cro_Lim=3){

  void *dir = gSystem->OpenDirectory(gSystem->ExpandPathName(folder));
  TString logname, buffer;
  double smoothings[70][7]; // 0 smoo_m, 1 ini_m, 2 fin_m, 3 smoo_p, 4 ini_p, 5 fin_p
  int counts[70][7];        // 0 mTimes, 1 pTimes, 2 nFinished, 3 nMaximum
  initialize(smoothings,counts);
  if (dir) {
    int ndir = 0;
    while ((logname = gSystem->GetDirEntry(dir)) && ndir < 222) {
      if (logname=="") break;
      if (logname=="." || logname==".." || !logname.Contains("log")) continue;
      TString filename = folder; filename += logname;
      TString sam = logname; sam.Remove(0,sam.First('_')+1);
      sam.Remove(sam.First('.'),sam.Length());
      int Sam = sam.Atoi();
      fstream file;
      file.open(filename,fstream::in);
      double minCV = 1e10, smoo_m = -99, smoo_p = -99, TotminCV = 1e10;
      int n_m = 1, n_p = 1, nFinished = 0;
      while(file){
	file>>buffer;
	if(buffer == "Minimum") counts[Sam][3]++;
	if(buffer == "File") {
	  counts[Sam][4]++;
	  if(counts[Sam][4]>1)  {
	    for(int j=0; j<7; j++) smoothings[Sam][j] = -1;
	    nFinished = 0; minCV = 999;
	  }
	}
	if(buffer == "Finished") {
	  nFinished++;
	  if(nFinished == 1){
	    for(int i = 0; i < 7; i++) file>>buffer;
	    smoothings[Sam][2] = buffer.Atof();
	    for(int i = 0; i < 3; i++) file>>buffer;
	    smoothings[Sam][5] = buffer.Atof();
	  }else{
	    for(int i = 0; i < 10; i++) file>>buffer;
	    smoothings[Sam][5] = buffer.Atof();
	  }
	}
	if(buffer.Contains("/") && buffer.Contains(",")) {
	  int posToCV = 3;
	  TString pos_m = buffer;
	  pos_m.Remove(pos_m.First('/'),pos_m.Length());
	  if(pos_m=="1" && nFinished==0){
	    buffer.Remove(0,buffer.First('/')+1);
	    buffer.Remove(buffer.Length()-1,buffer.Length());
	    n_m = buffer.Atoi();
	    file>>buffer; posToCV--;
	    buffer.Remove(0,buffer.First('/')+1);
	    buffer.Remove(buffer.Length()-1,buffer.Length());
	    n_p = buffer.Atoi();
	  }
	  file>>buffer;
	  if(!(buffer=="CV") && !(buffer.Contains('/') && buffer.Contains(':'))) continue;
	  posToCV--;
	  for(int i = 0; i < posToCV; i++) file>>buffer;
	  double CV = buffer.Atof();
	  if(CV <= minCV){
	    minCV = CV;
	    for(int i = 0; i < 4; i++) file>>buffer;
	    smoo_m = buffer.Atof();
	    for(int i = 0; i < 3; i++) file>>buffer;
	    smoo_p = buffer.Atof();
	    if(pos_m=="1" && nFinished==0){
	      smoothings[Sam][1] = smoo_m; smoothings[Sam][4] = smoo_p; 
	    }
	  }
	  if(CV <= TotminCV)  TotminCV = CV;
	} // buffer contains '/' or '.'
      }
      counts[Sam][0] = n_m; counts[Sam][1] = n_p; counts[Sam][2] = nFinished; 
      smoothings[Sam][0] = smoo_m; smoothings[Sam][3] = smoo_p; 
      if(TotminCV<minCV) smoothings[Sam][6] = 1;
      ndir++;
    }
  } else cout<<"The directory does not exist"<<endl;
  for(int i = 0; i<69; i++){
    if(smoothings[i][0]>0){
      int rep_m = counts[i][0]+1, rep_p = counts[i][1]+1;
      double ini_m = smoothings[i][0], fin_m = smoothings[i][0]; 
      double ini_p = smoothings[i][3], fin_p = smoothings[i][3];
      if(counts[i][2]==0){  // It didn't even finish the first smoo_p round
	ini_m = smoothings[i][0]*0.7; 
	fin_m = smoothings[i][0]*1.4;
	ini_p = smoothings[i][3]*0.9;
	fin_p = smoothings[i][3]*2;
	rep_p = (int) sqrt( ((double)rep_m) / 2); 
	rep_m = 2*rep_p;
      } else { 
	if(areEqual(smoothings[i][0], smoothings[i][1])) ini_m = smoothings[i][0]*0.6;
	if(areEqual(smoothings[i][0], smoothings[i][2])) {
	  ini_m = smoothings[i][0]*0.9;
	  fin_m = smoothings[i][0]*1.6;
	}
	if(areEqual(smoothings[i][3], smoothings[i][4])) ini_p = smoothings[i][3]*0.6;
	if(areEqual(smoothings[i][3], smoothings[i][5])) {
	  ini_p = smoothings[i][3]*0.9;
	  fin_p = smoothings[i][3]*1.8;
	}
	if(counts[i][3]==0){
	  rep_p = (int) sqrt( ((double)(counts[i][2]*rep_m)) / 2);
	  rep_m = 2*rep_p;
	}
      }
      if(Fit_Cro_Lim==1){
	double smoo_m = smoothings[i][0], smoo_p = smoothings[i][3];
	cout<<RoundNumber(smoo_m*100,0)<<" \t"<<RoundNumber(smoo_m*smoo_p*100,0)<<endl;
// 	cout<<"bsub -q $q -o $logdir/logFit_"<<i<<".log FitPdfSampleKEYS "<<i<<
// 	  "\t"<<RoundNumber(smoo_m*100,0)<<" "<<
// 	  RoundNumber(smoo_p*100,0)<<"\t700 200 "<<RoundNumber(smoo_m*60,0)<<" "<<
// 	  RoundNumber(smoo_p*60,0)<<"\t"<<RoundNumber(stitch1,1)<<"  "<<RoundNumber(stitch2,1);
// 	if(i>=12&&i<16 || i==20 || i==22) cout<<" noRot";
// 	cout<<"\t # Sample "<<i;
// 	cout<<endl; 
      }
      if(Fit_Cro_Lim==2){//smoothings[i][4]<=smoothings[i][5]){
	cout<<"bsub -q $"; if(i>=9&&i<=13 || i==15 || i==21) cout<<"x"; else cout<<"q";
	cout<<" -o $logdir/logCrossV_"<<i<<".log CrossValidation "<<i<<
	  " \t"<<rep_m<<" "<<RoundNumber(ini_m*90,0)<<" "<<
	  //	  RoundNumber(fin_m*110,0)<<" \t"<<rep_p<<" "<<RoundNumber(ini_m*ini_p*80,0)<<" "<<
	  //	  RoundNumber(fin_m*fin_p*120,0);
	  RoundNumber(fin_m*110,0)<<" \t"<<rep_p<<" "<<RoundNumber(ini_p*80,0)<<" "<<
	  RoundNumber(fin_p*120,0);
	cout<<"\tnoRot RAllDssonDx200   \t# Sample "<<i<<"c";
	cout<<endl; 
      }
      if(Fit_Cro_Lim==3){//smoothings[i][4]<=smoothings[i][5]){
	TString var[] = {"Mmiss: ", "|", "|", "\tPl: ", "|", "|"};
	cout<<"Sample "<<i<<" -\t";
	int Indices[] = {1,0,2,4,3,5};
	for(int j=0; j<6; j++) cout<<var[j]<<" "<<RoundNumber(smoothings[i][Indices[j]],2)<<" ";
	cout<<"\t Counts - "; 
	for(int j=0; j<5; j++) cout<<counts[i][j]<<", ";
	if(smoothings[i][6]>0) cout<<"\t- There was another minimum";
	if(smoothings[i][0]>smoothings[i][1] && smoothings[i][0]<smoothings[i][2] &&
	   smoothings[i][3]>smoothings[i][4] && smoothings[i][3]<smoothings[i][5])cout<<"\t Done";
	cout<<endl; 
      }
      if(i==32 || i==-20) cout<<endl; 
    } 
  }
  
} // 0 smoo_m, 1 ini_m, 2 fin_m, 3 smoo_p, 4 ini_p, 5 fin_p


TString RoundNumber(double n, int e, double d){
  if(d==0) return " - ";
  double b = (int)(n/d*pow(10.,(double)e)+0.5);
  b /= pow(10.,(double)e);
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(!result.Contains(".") && e != 0) result += ".";
  
  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<e-afterdot.Length(); i++)
    result += "0";
  return result;
}

bool areEqual(double a, double b){
  bool equal = true;
  if(fabs(a-b)>0.00001) equal = false;

  return equal;
}

void initialize(double smoothings[][7], int counts[][7]){
  for(int i = 0; i<70; i++)
    for(int j = 0; j<7; j++){
      smoothings[i][j] = -1; 
      counts[i][j] = 0;
    }
  return;
}
