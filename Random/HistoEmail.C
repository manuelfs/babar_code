#include "TCanvas.h"
#include "TH1F.h"
#include "TPad.h"
#include "TString.h"
#include "babar_code/Styles/Styles.cc"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;
#define nmaxEmail 1000

void HistoEmail(int isOut = 0){

  Styles style; style.setPadsStyle(2); style.applyStyle();
  TString buf, Sender[nmaxEmail];
  TString months[] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", 
		      "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
  int Date[nmaxEmail][3], MonEmail[2][50];
  for(int mon=1; mon<=48; mon++) for(int em=0; em<2; em++) MonEmail[mon][em] = 0;
  TString sSender = "ltr", baseName = "babar_code/Random/Emails/Gmail";
  int nFiles = 7;
  if(isOut){
    sSender = "To:";
    baseName = "babar_code/Random/Emails/Gout";
    nFiles = 5;
  }
  for(int file=1; file<nFiles; file++){
    TString fileName = baseName; fileName += file; 
    fstream textFile; textFile.open(fileName,fstream::in);
    int nWords = 0, nEmail=0;
    while(textFile && nEmail <= 100){
      textFile>>buf;
      if(buf.Contains("checkbox")) {
	// Sender
  	while(!buf.Contains(sSender)) textFile>>buf;
 	Sender[nEmail] = buf; Sender[nEmail].Remove(0,Sender[nEmail].First('>')+1);
	if(isOut) Sender[nEmail] = "";
  	while(!buf.Contains("<")) {
	  textFile>>buf;
	  Sender[nEmail] += " "; Sender[nEmail] += buf;
	}
	Sender[nEmail].Remove(Sender[nEmail].First('<'), Sender[nEmail].Sizeof());
	// Date
  	while(!buf.Contains("nowrap")) textFile>>buf;
	while(!buf.Contains("nbsp") && !buf.Contains("/"))textFile>>buf;
	if(!buf.Contains("/")){
	  Date[nEmail][0] = 11;
	  for(int mon=1; mon<=12; mon++) if(buf.Contains(months[mon-1])){Date[nEmail][1] = mon; break;}
	  buf.Remove(0,9); Date[nEmail][2] = buf.Atoi();
	} else {
	  TString date = buf; date.Remove(date.First("/"),date.Sizeof()); Date[nEmail][1] = date.Atoi();
	  date = buf; date.Remove(0,date.First("/")+1); 
	  date.Remove(date.First("/"),date.Sizeof()); Date[nEmail][2] = date.Atoi();
	  date = buf; date.Remove(0,date.Last('/')+1); Date[nEmail][0] = date.Atoi();
	}
	if(Sender[nEmail].Contains("achel")) MonEmail[0][Date[nEmail][1]+(Date[nEmail][0]-8)*12]++;
//  	cout<<Date[nEmail][0]<<"/"<<Date[nEmail][1]<<"/"<<Date[nEmail][2]
// 	    <<"\t"<<Sender[nEmail]<<endl;
	nEmail++;
      }
      //cout<<buf<<endl;
      nWords++;
    }
  }
  gStyle->SetPadLeftMargin(0.1); gStyle->SetPadBottomMargin(0.12);
  TCanvas can("can","RD results",800,600);
  can.Divide(1,2); TPad *cPad; can.cd(1);
  TH1F histo("histo","",39,0,39);
  for(int mon=1; mon<=39; mon++) histo.SetBinContent(mon,MonEmail[0][mon+5]);
  histo.Draw();
  can.SaveAs("babar_code/Random/RM.eps");

}



