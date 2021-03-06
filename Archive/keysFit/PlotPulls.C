#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"
#include "babar_code/Styles/Styles.cc"
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;
using namespace TMath; 
using std::cout;
using std::endl;

TString RoundNumber(double n, int e, double d=1);
void readPulls(TString textName, double Pulls[8][9], double Errors[8][9], double TotalYield[8][2], int plotYields);
double truYields[8][9];
int rmsPull_1 = 0;
TString fileTrue = "FitAll/fits/TextFinalRAll.txt";

void PlotPulls(int isHalf = 1,double maxPull = 2.8, int plotYields=0, int doR24=0, int doRMS1=1){

  rmsPull_1 = doRMS1;
  if(rmsPull_1 && isHalf) rmsPull_1++;
  Styles style6; style6.setPadsStyle(6); 
  vector<Double_t> vPulls[8][9], vYield[8][2];
  double Pulls[20][20][8][9], Errors[20][20][8][9], tPulls[8][9], tErrors[8][9], TotalYield[8][2];
  readPulls(fileTrue, tPulls, tErrors, TotalYield, plotYields);

  int nFiles[] = {2,7,7}, iniFolder = 1, nFolders = 3;
  if(doR24==1) {iniFolder = 0; nFolders = 1;}
  if(doR24==2) {iniFolder = 0; nFolders = 3;}
   TString folder = "FitAll/fits/";
   //if(isHalf) folder = "FitAll/TextFilesHalf/";
   //  TString folder = "FitAll/Archive/Out-of-the-box/TextFiles/";
   //if(isHalf) folder = "FitAll/Archive/Out-of-the-box/TextFilesHalf/";
  TString RelName[] = {"TextFinalR24","TextFinalR26","TextFinalR26"};
  for(int i=iniFolder; i<nFolders; i++) {
    int maxFiles = nFiles[i]*(isHalf+1);
    int iniRest = (i==2)*3*(isHalf+1), finRest = maxFiles-(i==1)*4*(isHalf+1);
    for(int rest=iniRest; rest<finRest; rest++){
      TString textName = folder; textName+=RelName[i]; textName+="Rests_"; 
      textName += rest+1; textName+="_"; textName += maxFiles; textName+=".txt";
      //cout<<"Reading "<<textName<<endl;
      readPulls(textName, tPulls, tErrors, TotalYield, plotYields);
      for(int j=0; j<8; j++){
	vYield[j][0].push_back(TotalYield[j][0]);
	vYield[j][1].push_back(TotalYield[j][1]);
	for(int m=0; m<8; m++) {
	  Pulls[i][rest][m][j] = tPulls[m][j];
	  Errors[i][rest][m][j] = tErrors[m][j];
	  vPulls[m][j].push_back(tPulls[m][j]);
	}
      }
    }
  }
  readPulls(fileTrue, tPulls, tErrors, TotalYield, plotYields);
  for(int j=0; j<8; j++){
    double mean = 0, rms = 0;
    int vSize = vYield[j][0].size();
    for(int pulli=0; pulli<vSize; pulli++) {
      double diff = vYield[j][1][pulli]-vYield[j][0][pulli];
      mean += diff;
      rms += diff*diff;
    }
    mean /= vSize; rms = sqrt((rms - mean*mean*vSize)/(vSize-1)/vSize);
    cout<<j<<": "<<RoundNumber(mean,2)<<" +- "<<RoundNumber(rms,2)<<endl;
  }
  for(int j=0; j<8; j++)for(int m=0; m<8; m++) sort(vPulls[m][j].begin(), vPulls[m][j].end());
  TString titles[6] = {"Signal","Normalization", "Feed-down","D** #Rightarrow D#pi^{0}",
		       "D* #Rightarrow D*#pi^{0}","Comb. #Rightarrow D#pi^{0}"};
//   TString titles[6] = {"Signal:  <#sigma_{stat}> = ","Normalization:  <#sigma_{stat}> = ",
// 		       "Feed-down:  <#sigma_{stat}> = ","D** #Rightarrow D#pi^{0}:  <#sigma_{stat}> = ",
// 		       "D* #Rightarrow D*#pi^{0}:  <#sigma_{stat}> = ","Comb. #Rightarrow D#pi^{0}:  <#sigma_{stat}> = "};
  TString tags[4] = {"D^{0}","D*^{0}","D^{+}","D*^{+}"};
  int indices[2][6] = {{0,1,3,0,1,4},{0,0,0,4,4,4}};
  //int indices[2][6] = {{0,1,3,4,5,6},{0,0,0,0,0,0}};
  int colors[4] = {4,9,9,9};
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.25); if(plotYields) gStyle->SetPadLeftMargin(0.22);
  gStyle->SetPadRightMargin(0.01); 
  gStyle->SetTitleFont(style6.nFont,"xy");          // Set the all 2 axes title font
  gStyle->SetLabelFont(style6.nFont,"xy");          // Set the all 2 axes label font
  gStyle->SetTitleSize(style6.TitleSize,"xy");      // Set the 2 axes title size
  gStyle->SetTitleOffset(style6.xTitleOffset*0.93,"x");     
  gStyle->SetPadBottomMargin(style6.PadBottomMargin*0.8); gStyle->SetPadTopMargin(0.11); 
  gStyle->SetNdivisions(306, "x");
  TCanvas c("Pulls","Pulls",style6.CanvasW,style6.CanvasH);
  c.Divide(2,3);
  TMarker mark; TLine lin; TLine linh; 
  TLatex label; label.SetTextSize(0.08); label.SetTextFont(style6.nFont);
  mark.SetMarkerStyle(8); mark.SetMarkerSize(0.75);
  TH1F* h[6];
  TH1F hPull("hPull","Pulls",100,-3.5,3.5);
  double chi2=0, rmsAll=0, meanAll=0; int ndof=0, vSizeAll=0;
  for(int i=0; i<6; i++) {
    c.cd(i+1);
    double Emean = 0; int nErrors = 0;
    TString hname = "histo"; hname += i;hname += isHalf;hname += plotYields;
    h[i] = new TH1F(hname,"",20,-maxPull,maxPull);
    h[i]->SetMaximum(4);
    h[i]->SetXTitle("Pull");
    if(plotYields) h[i]->SetXTitle("Bias (Events)");
    h[i]->GetXaxis()->SetLabelSize(0.07);
    h[i]->GetYaxis()->SetNdivisions(0);
    h[i]->Draw();
    double minX = (int)(0.1-maxPull), dX = 1;
    if(plotYields){
      dX = 20;
      //if(maxPull<45) dX = 15;
      minX = -(int)(maxPull/dX)*dX;
    }
    for(double x=minX; x<maxPull; x += dX) {
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
      int IndCha = indices[1][i]+chan, IndRow = indices[0][i];
      if(IndRow==3 && IndCha%2==1) continue;
      label.DrawLatex(-1.35*maxPull,3.45-(double)chan,tags[chan]);
      double height = 3.75-(double)chan;
      for(int fold=iniFolder; fold<nFolders; fold++){
	double minp=99., maxp=-99.;
	int maxFiles = nFiles[fold]*(isHalf+1);
	int iniRest = (fold==2)*3*(isHalf+1), finRest = maxFiles-(fold==1)*4*(isHalf+1);
	for(int rest=iniRest; rest<finRest; rest++){
	  double pull = Pulls[fold][rest][IndCha][IndRow];
	  Emean += Errors[fold][rest][IndCha][IndRow]; nErrors++;
	  if(pull>maxp) maxp = pull;
	  if(pull<minp) minp = pull;
	  mark.SetMarkerColor(colors[fold]); mark.SetMarkerSize(0.6);
	  if(fabs(pull)<maxPull) mark.DrawMarker(pull,height);
	  //cout<<fold<<", "<<rest<<", "<<IndCha<<", "<<IndRow<<": Plotting pull "<<pull<<endl;
	}
	if(maxp>maxPull) maxp = maxPull; if(minp<-maxPull) minp = -maxPull; 
	linh.SetLineColor(colors[fold]);
	if(!(IndRow==3 && IndCha%2==1)) linh.DrawLine(minp,height,maxp, height);
	height -= 0.17;
      }
      double mean = 0, rms = 0;
      int vSize = vPulls[IndCha][IndRow].size();
      for(int pulli=0; pulli<vSize; pulli++) {
	hPull.Fill(vPulls[IndCha][IndRow][pulli]);
	mean += vPulls[IndCha][IndRow][pulli];
	rms += vPulls[IndCha][IndRow][pulli]*vPulls[IndCha][IndRow][pulli];
      }
      meanAll += mean; rmsAll += rms; vSizeAll += vSize;
      mean /= vSize; rms = sqrt((rms - mean*mean*vSize)/(vSize-1)/vSize); // Uncertainty on the mean: sigma/sqrt(n)
      //if(rmsPull_1) rms = 1/sqrt(vSize);
      if(rmsPull_1) rms = 0.906/sqrt(vSize);
      chi2 += pow(mean/rms,2); ndof++;
      mark.SetMarkerSize(0.7);mark.SetMarkerColor(2);
      mark.DrawMarker(mean, height); mark.SetMarkerColor(28);
      mark.SetMarkerSize(0.8);mark.SetMarkerColor(28);
      //mark.DrawMarker(tPulls[IndCha][IndRow], height-0.18); mark.SetMarkerColor(28);
      linh.SetLineColor(2); linh.SetLineWidth(2);
      linh.DrawLine(mean-rms,height,mean+rms, height); linh.SetLineWidth(1);
      TString labMean = RoundNumber(mean,2-plotYields); labMean += " #pm "; labMean += RoundNumber(rms,2-plotYields);
      label.SetTextSize(0.065);label.SetTextAlign(31);
      //if(fabs(mean) > 2*rms) {label.SetTextColor(2);label.SetTextFont(22);}
      label.DrawLatex(-1.07*maxPull,3.75-(double)chan-0.17*3-0.1,labMean);label.SetTextSize(0.08);
      label.SetTextAlign(11);label.SetTextColor(1);label.SetTextFont(style6.nFont);
    }
    label.SetTextSize(0.07);
    //titles[i] += (int)(Emean/nErrors+1); 
    //if(plotYields) titles[i] += " events"; else titles[i] += "%";
    label.DrawLatex(-0.8*maxPull,4.2,titles[i]);
    label.SetTextSize(0.08);    
  }
  c.cd(0);label.SetTextSize(0.022);label.SetTextAlign(23);
//   TString chiText = "#splitline{#chi^{2}: "; chiText += RoundNumber(chi2,1); chiText += "/";
//   chiText += ndof; chiText += " = "; chiText += RoundNumber(chi2,2,ndof);
//   chiText += "}{Prob. = "; chiText += RoundNumber(Prob(chi2,ndof)*100,2); chiText += "%}";
  TString chiText = "#chi^{2} = "; chiText += RoundNumber(chi2,1); chiText += "/";
  chiText += ndof; chiText += "; p = "; chiText += RoundNumber(Prob(chi2,ndof)*100,2); chiText += "%";
  label.DrawLatex(0.54,0.991,chiText);
  TString pName = "FitAll/PlotPulls"; pName += isHalf; pName += plotYields; pName += doR24;
  pName += ".eps";
  c.SaveAs(pName);
  for(int i=0; i<6; i++) h[i]->Delete();
  for(int i=0; i<8; i++) {
    for(int j=0; j<8; j++) vPulls[i][j].clear(); 
    for(int j=0; j<2; j++) vYield[i][j].clear();
  }
  TF1 fGaus("Gaussian","gaus",-4,4);
  hPull.Fit(&fGaus,"L Q");
  double Mean = fGaus.GetParameter(1), eMean = fGaus.GetParError(1);
  double Sigma = fGaus.GetParameter(2), eSigma = fGaus.GetParError(2);
  cout<<"hPull has "<<hPull.GetEntries()<<" entries: \t Mean is "<<RoundNumber(Mean,3)<<" +- "<<
    RoundNumber(eMean,3)<<" and Sigma "<<RoundNumber(Sigma,3)<<" +- "<<RoundNumber(eSigma,3)<<
    "\t\t Prob: "<< RoundNumber(Prob(chi2,ndof)*100,2)<<"%"<<endl;
}

void readPulls(TString textName, double Pulls[8][9], double Errors[8][9], double TotalYield[8][2], int plotYields){
  double Yields[16][9];
  for(int i=0; i<8; i++){
    TotalYield[i][0] = 0;
    TotalYield[i][1] = 0;
    for(int j=0; j<9; j++){
      Pulls[i][j] = 0;
      Errors[i][j] = 0;
    }
  }
  TString buffer, buffer2;
  fstream textFile; textFile.open(textName,fstream::in);
  int begChan = 0, endChan = 4, IndText[8] = {0,4,1,5,2,6,3,7};
  for(int i=2*begChan; i<2*endChan; i++){
    int index = IndText[i];
    for(int j=0; j<9; j++){
      textFile>>buffer>>Yields[index+8][j]>>Yields[index][j]>>buffer>>buffer2>>buffer;
      if(rmsPull_1){
	if(!textName.Contains(fileTrue)) Yields[index+8][j] = truYields[index][j]/(double)rmsPull_1;
	else truYields[index][j] = Yields[index+8][j];
      }
      TotalYield[index][0] += Yields[index+8][j]; 
      TotalYield[index][1] += Yields[index][j]; 
      if(buffer != "-" || buffer != "fixed") {
	if(plotYields) Pulls[index][j] = Yields[index][j]-Yields[index+8][j];
	else {
	  if(rmsPull_1==0) Pulls[index][j] = buffer.Atof();
	  else Pulls[index][j] = (Yields[index][j]-Yields[index+8][j])/buffer2.Atof();
	}
      }
      if(buffer2 != "-") {
	Errors[index][j] = buffer2.Atof();
	if(plotYields==0) Errors[index][j] /= (Yields[index+8][j]/100.);
      }
      //cout<<"True: "<<Yields[index+8][j]<<"\tFit: "<<Yields[index][j]<<"\tPull: "<<Pulls[index][j]<<"\tChannel "<<index<<endl;
      if(index>3 && j==6){
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
