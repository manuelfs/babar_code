//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: VaryConstraints.cc,v 1.5 2012/08/23 02:22:18 manuelf Exp $
//
// Description:
//      VaryConstraints - Varies the cross-feed constraints by Poisson fluctuating
//      the yields
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      10/04/28 manuelf -- Created
//------------------------------------------------------------------------

#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TIterator.h"
#include "TRandom3.h"
#include "DonutUtils/KeysUtils.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

Bool_t isFloated(int channel);
int floatDinDss=1, IsoConst=0;

int main(int argc, char *argv[]){
  if (argc < 2 || argc > 5 ) {
    cout << "USAGE: VaryConstraints RootFile [IsoConst=0] [Type=vary] [nTimes=10,000]" << endl
	 << "TYPES: Table, writes values - vary, 4 Stat - Crossfeed, 10 Crossfeed - Iso, 11 Isospin "<<endl;;
    return 0;
  }
  // Setting the input parameters
  TString Rootname = argv[1];
  TString ISO = "0";
  if (argc>2) {ISO = argv[2]; IsoConst = ISO.Atoi();}
  TString Type = "vary";
  if (argc>3) Type = argv[3];
  int ntimes = 10000;
  if (argc>4) {TString temp = argv[4]; ntimes = temp.Atoi();}

  if(Type.Contains("Iso")) IsoConst = 1;
  TChain *treeData = new TChain("ntp1");
  treeData->Add(Rootname);
  int IndPdf[2][9] = {{1,9,5,13,17,41,51,61,37},{21,25,29,45,55,65,33,-1,-1}};

  Int_t MCType,candType;
  float candM2,candPstarLep, weight=-1;
  treeData->SetBranchAddress("MCType",&MCType);
  treeData->SetBranchAddress("candType",&candType);
  treeData->SetBranchAddress("candM2",&candM2);
  treeData->SetBranchAddress("candPstarLep",&candPstarLep);
  treeData->SetBranchAddress("weight",&weight);

  long entries = treeData->GetEntries();
  double truYield[70], truYieldFix[70], Entries[70];
  int IndfYield[70];
  for(int i=0; i<70; i++) {truYieldFix[i] = 0; Entries[i] = 0;}
  TRandom3 rand(0); 
  for (int evt = 0 ; evt < entries; evt ++) {
    treeData->GetEvent(evt);
    if(candM2<-4||candM2>12 || candPstarLep<0||candPstarLep>2.4 || weight<0||weight>100) continue;
    for(int i=0; i<70; i++) 
      if(isSample(i,MCType,candType)) {
	truYieldFix[i] += weight;
	Entries[i]++;
      }
  }
  TString outName = "babar_code/Systematics/Text/ConstraintsError.txt", buffer;
  fstream outFile; outFile.open(outName,fstream::in);
  double CrossfeedError[4], IsoRandD = 1., IsoRandDs = 1.;
  for(int cand=0; cand<4; cand++){
    outFile>>buffer;
    CrossfeedError[cand] = buffer.Atof();
  }
  TString textName = "keys/text/"; textName+=Type; textName+="_"; 
  if(IsoConst) textName +="Iso"; 
  textName+="Const.txt";
  fstream os;
  os.open(textName,fstream::out);
  if(Type == "Table") ntimes = 1;
  for(int rep = 0; rep < ntimes; rep++){
    for(int i=0; i<70; i++) {
      if(Type == "Table") truYield[i] = truYieldFix[i];
      else if(Type == "Crossfeed"){
	if(i>=5&&i<=8 || i==14 || i==16){
	  truYield[i] = truYieldFix[i]*rand.Gaus(1,CrossfeedError[(i-1)%4-(i>10)]);
	  //if(rep==0) cout<<i<<": yield "<<truYieldFix[i]<<" with RMS "<<CrossfeedError[(i-1)%4-(i>10)]<<endl;
	  if(truYield[i]<0) truYield[i] = 0;
	} else truYield[i] = truYieldFix[i];
      } else if(Type == "Iso"){
	if(i==3 || i==11 || i==15){
	  if(i==3) IsoRandD = rand.Gaus(1,0.034);
	  truYield[i] = truYieldFix[i]*IsoRandD;
	  if(truYield[i]<0) truYield[i] = 0;
	} else if(i==4 || i==12){
	  if(i==4) IsoRandDs = rand.Gaus(1,0.036);
	  truYield[i] = truYieldFix[i]*IsoRandDs;
	  if(truYield[i]<0) truYield[i] = 0;
	} else truYield[i] = truYieldFix[i];
      } else truYield[i] = rand.PoissonD(Entries[i])*truYieldFix[i]/Entries[i];
    }
    for(int i=0; i<9; i++){
      for(int chan=0; chan<8; chan++){
	if(IndPdf[chan>3][i]<0) continue;
	int index = IndPdf[chan>3][i]+chan-4*(chan>3);
	if(!isFloated(index)) {
	  if(rep==0){
	    if(IsoConst){
	      if(index>=5 && index<=6 || index==14) IndfYield[index] = index-5+2*(index%2);
	      else if(index>=7 && index<=8 || index==16) IndfYield[index] = index-5+2*(index%2)-2;
	      else if(index== 3 || index== 4 || index==11 || index==12 || index==15 || 
		      index==23 || index==24 || floatDinDss&&index>=27&&index<=28) IndfYield[index] = index-2;
	      else if(index>=17 && index<=20) IndfYield[index] = index+4;    // D**lnu all floated 
	      //else if(index>=17 && index<=18) IndfYield[index] = index+4;        // D**lnu not all floated 
	      //else if(index>=19 && index<=20) IndfYield[index] = index+2;        // D**lnu not all floated 
	      else if(index>=25 && index<=26) IndfYield[index] = index-16+(index%2)*4;
	      else if(index>=27 && index<=28) IndfYield[index] = index-18+(index%2)*4;
	      else if(index>=29 && index<=30){ 
		if(floatDinDss) IndfYield[index] = index-4;
		else IndfYield[index] = index-20-(index%2);
	      } else if(index>=31 && index<=32){ 
		if(floatDinDss) IndfYield[index] = index-6;
		else IndfYield[index] = index-22-(index%2);
	      }
	      else if(index>=37 && index<=40) IndfYield[index] = index-4;
	      //else if(index>=35 && index<=36) IndfYield[index] = index-2;
	      //else if(index>=37 && index<=38) IndfYield[index] = index-4;
	      //else if(index>=39 && index<=40) IndfYield[index] = index-6;
	    } else {
	      if(index>=5 && index<=8 || index==14 || index==16) IndfYield[index] = index-5+2*(index%2);
	      else if(index>=17 && index<=20) IndfYield[index] = index+4;
	      else if(index>=25 && index<=28) IndfYield[index] = index-16+(index%2)*4;
	      else if(index>=29 && index<=32){ 
		if(floatDinDss) IndfYield[index] = index-4;
		else IndfYield[index] = index-20-(index%2==0);
	      }
	      else if(index>=37 && index<=40) IndfYield[index] = index-4;
	      else if(index>=35 && index<=36) IndfYield[index] = index-2;
	    } //IsoConst
	  }
	  double n = truYield[index], d = truYield[IndfYield[index]];
	  double error = sqrt(n/d/d+n*n/d/d/d);
	  if(d<=0) d=1;
	  os<<index<<"\t"<<RoundNumber(n,6,d);
	  if(Type=="Table") os<<"  +-  "<<RoundNumber(error,6);
	  os<<endl;
	}
      }
    }
    os<<endl;
  }
  cout<<"Written "<<textName<<endl;
}

Bool_t isFloated(int i){
  if(IsoConst){
    if(IsoConst==2){
      if(i>=1&&i<=2 || i>=9&&i<=13 || i==15 || i>=21&&i<=24 || floatDinDss==1&&i>=25&&i<=28 || i>=33&&i<=36 || i>=41&&i<=68) return true;
      else return false;
    } else {
      if(i>=1&&i<=2 || i>=9&&i<=10 || i==13 || i>=21&&i<=24 || floatDinDss==1&&i>=25&&i<=26 || i>=33&&i<=36 || i>=41&&i<=68) return true;
      else return false;
    }
  }
  if(i>=1&&i<=4 || i>=9&&i<=13 || i==15 || i>=21&&i<=24 || floatDinDss==1&&i>=25&&i<=28 || i>=33&&i<=36 || i>=41&&i<=68) return true;
  return false;
}


