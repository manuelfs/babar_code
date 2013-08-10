#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TLine.h"
#include "TLatex.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;
using std::cout;
using std::endl;

TString RoundNumber(double n, int e, double d=1);
void readPulls(TString textName, double Pulls[8][7], double Errors[8][7]);

void PlotYiels(int half = 1, double maxPull = 57){
  //TString folders[4] = {"All","TruthMatch_Rests_0_2","Rests_0_2","R22_R24"};
  //TString folders[4] = {"Sep_Rests_0_2","Rests","Rests_0_2","R22_R24"};
  TString folders[4] = {"Rests6","Rests8","Rests2","Rests4"};
  //TString files[6] = {"Rests_0_2","Rests_1_2","Rests2_0_2","Rests2_1_2","MVAs_0_2","MVAs_1_2"};
  //TString files[6] = {"Rests_1_4","Rests_2_4","Rests_3_4","Rests_4_4","Rest","Rest2"};
  //TString files[6] = {"Rests_2_6","Rests_3_6","Rests_5_6","Rests_6_6","Rest","Rest2"};
  TString files[4][8] = {{"Rests_2_6","Rests_3_6","Rests_5_6","Rests_6_6",
			  "Rests_5_8","Rests_6_8","Rests_7_8","Rests_8_8"}
			 ,{"Rests_1_8","Rests_2_8","Rests_3_8","Rests_4_8",
			   "Rests_5_8","Rests_6_8","Rests_7_8","Rests_8_8"}
			 ,{"Rest","Rest2","Rests_3_8","Rests_4_8",
			   "Rests_5_8","Rests_6_8","Rests_7_8","Rests_8_8"}
			 ,{"Rests_1_4","Rests_2_4","Rests_3_4","Rests_4_4",
			   "Rests_5_8","Rests_6_8","Rests_7_8","Rests_8_8"}};
  int maxFolder[] = {4,8,2,4};
  vector<Double_t> vPulls[8][7];
  double Pulls[5][12][8][7], Errors[5][12][8][7], tPulls[8][7], tErrors[8][7];
  int folder1 = 0, folder2 = 2;
  if(!half){folder1 = 2; folder2 = 4;}
  for(int i=0; i<8; i++) {
    for(int rest=folder1; rest<folder2; rest++){
      if(i>=maxFolder[rest]) continue;
      TString textName = "FitAll/"; textName+=folders[rest]; textName+="/TextFinal"; 
      textName += files[rest][i]; 
      textName+=".txt";
      cout<<"Reading "<<textName<<endl;
      readPulls(textName, tPulls, tErrors);
      for(int j=0; j<7; j++)for(int m=0; m<8; m++) {
	Pulls[rest][i][m][j] = tPulls[m][j];
	Errors[rest][i][m][j] = tErrors[m][j];
	//if(rest != 1 && rest != 3) 
	  vPulls[m][j].push_back(tPulls[m][j]);
      }
    }
//     TString textName = "FitAll/Rests/TextFinal"; 
//     textName+=folders[i]; textName+=".txt";
//     readPulls(textName, tPulls);
//     for(int j=0; j<7; j++)for(int m=0; m<8; m++) Pulls[4][i][m][j] = tPulls[m][j];
  }
  for(int j=0; j<7; j++)for(int m=0; m<8; m++) sort(vPulls[m][j].begin(), vPulls[m][j].end());
  TString titles[4] = {"Signal:  #sigma_{stat} ~ ","Normalization:  #sigma_{stat} ~ ",
		       "Feed-down:  #sigma_{stat} ~ ","D**:  #sigma_{stat} ~ "};
  TString tags[4] = {"D*^{0}","D^{0}","D*^{+}","D^{+}"};
  int indices[2][4] = {{0,1,3,0},{0,0,0,4}};
  int colors[4] = {4,9,4,9};
  TCanvas c("Yields","Yield error",1000,700);
  c.Divide(2,2);
  gStyle->SetOptStat(0);
  TMarker mark; TLine lin; TLine linh; TLatex label;label.SetTextSize(0.08);
  mark.SetMarkerStyle(8); mark.SetMarkerSize(0.7);
  TH1F* h[4];
  for(int i=0; i<4; i++) {
    c.cd(i+1);
    double Emean = 0; int nErrors = 0;
    TString hname = "histo"; hname += i;
    h[i] = new TH1F(hname,"",20,-maxPull,maxPull);
    h[i]->SetMaximum(4);
    h[i]->SetXTitle("Bias");
    h[i]->GetXaxis()->SetLabelSize(0.07);
    h[i]->GetXaxis()->SetTitleOffset(0.6);
    h[i]->GetXaxis()->SetTitleSize(0.057);
    h[i]->GetYaxis()->SetNdivisions(0);
    h[i]->Draw();
    for(double x=-(int)(maxPull/20)*20; x<maxPull; x+=20) {
      if(x==0){
	lin.SetLineStyle(1);
	lin.SetLineWidth(2);
      } else {
	lin.SetLineStyle(2);
	lin.SetLineWidth(1);
      }
      lin.DrawLine(x,0,x,4);
    }
    for(int chan=0; chan<4; chan++){
      label.DrawLatex(-1.2*maxPull,3.4-(double)chan,tags[chan]);
      int IndCha = indices[1][i]+chan, IndRow = indices[0][i];
      if(IndRow==3 && IndCha%2==0) continue;
      double mean = 0, rms = 0, mheight = 3.81-(double)chan-0.44;
      int vSize = vPulls[IndCha][IndRow].size();
      for(int pulli=0; pulli<vSize; pulli++) {
	mean += vPulls[IndCha][IndRow][pulli];
	rms += vPulls[IndCha][IndRow][pulli]*vPulls[IndCha][IndRow][pulli];
      }
      mean /= vSize; rms = sqrt((rms/vSize - mean*mean)/vSize);
      mark.SetMarkerSize(0.6);mark.SetMarkerColor(2);
      mark.DrawMarker(mean, mheight); mark.SetMarkerColor(28);
      linh.SetLineColor(2); linh.DrawLine(mean-rms,mheight,mean+rms, mheight);

      double labLeft = -51; if(mean<0) labLeft -= 1.5; if(fabs(mean)>=10) labLeft -= 2; 
      label.SetTextSize(0.05);
      TString labMean = RoundNumber(mean,1); 
      label.DrawLatex(labLeft,mheight-0.12,labMean);
      labMean = " #pm  "; labMean += RoundNumber(rms,1);
      label.DrawLatex(-45,mheight-0.12,labMean);
      label.SetTextSize(0.08);
      for(int fold=folder1; fold<folder2; fold++){
	double minp=99., maxp=-99., height = 3.81-(double)chan-(double)(fold-folder1)*0.22;
	for(int rest=0; rest<8; rest++){
	  if(rest>=maxFolder[fold]) continue;
	  double pull = Pulls[fold][rest][IndCha][IndRow];
	  Emean += Errors[fold][rest][IndCha][IndRow]; nErrors++;
	  if(pull>maxp) maxp = pull;
	  if(pull<minp) minp = pull;
	  if(fold==4)mark.SetMarkerColor(4);
	  else mark.SetMarkerColor(colors[fold]);
	  mark.DrawMarker(pull,height);
	}
	if(fold==4)linh.SetLineColor(4);
	else linh.SetLineColor(colors[fold]);
	if(!(IndRow==3 && IndCha%2==0)) linh.DrawLine(minp,height,maxp, height);
      }
    }
    label.SetTextSize(0.07);
    titles[i] += (int)(Emean/nErrors+1); titles[i] += " events";
    label.DrawLatex(-0.9*maxPull,4.2,titles[i]);
    label.SetTextSize(0.08);    
  }
  TString pName = "FitAll/PlotYields"; pName += half; pName += ".eps";
  c.SaveAs(pName);
}

void readPulls(TString textName, double Pulls[8][7], double Errors[8][7]){
  double Yields[16][7];
  for(int i=0; i<8; i++)
    for(int j=0; j<7; j++){
      Pulls[i][j] = 0;
      Errors[i][j] = 0;
    }
  TString buffer, buffer2;
  fstream textFile; textFile.open(textName,fstream::in);
  int begChan = 0, endChan = 4;
  for(int i=2*begChan; i<2*endChan; i++){
    int index = i/2+(i%2)*4;
    for(int j=0; j<7; j++){
      textFile>>buffer>>Yields[index+8][j]>>Yields[index][j]>>buffer>>buffer2>>buffer;
      if(buffer != "-" || buffer != "fixed") Pulls[index][j] = Yields[index][j]-Yields[index+8][j];
      if(buffer2 != "-") Errors[index][j] = buffer2.Atof();
      //cout<<"True: "<<Yields[index+8][j]<<"\tFit: "<<Yields[index][j]<<"\tPull: "<<Pulls[index][j]<<"\tChannel "<<index<<endl;
      if(index>3 && j==4){
	textFile>>buffer;
	break;
      }
    }
    //cout<<endl;
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
