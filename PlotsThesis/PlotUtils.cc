#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMatrixT.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

// Reads a text file with the results of a fit
void ReadFitFile(TString textName, double Yield[2][70], double Error[70]){
  TString buffer, sError;
  fstream textFile; textFile.open(textName,fstream::in);
  for(int i=0; i<8; i++){
    for(int j=0; j<9; j++){
      textFile>>buffer;
      buffer.ReplaceAll("Yield[",""); buffer.ReplaceAll("]",""); 
      int index = buffer.Atoi();
      textFile>>Yield[1][index]>>Yield[0][index]>>buffer>>sError>>buffer;
      if(sError.Atof()>0) Error[index] = sError.Atof();
      else Error[index] = 0;
      //cout<<"True: "<<Yield[1][index]<<"\tFit: "<<Yield[0][index]<<" +- "<<Error[index]<<"\tIndex "<<index<<endl;
      if(i%2==1 && j==6){
	textFile>>buffer;
	break;
      }
    }
    //cout<<endl;
  }

}

void ReadCorrelFile(TString textName, double Correl[22][11]){
  int nRows = 22; if(textName.Contains("Iso")) nRows = 11;
  TString buffer;
  fstream textFile; textFile.open(textName,fstream::in);
  for(int row=0; row<nRows; row++){
    textFile>>buffer;    
    for(int col=0; col<11; col++){
      textFile>>Correl[row][col];
      //cout<< Correl[row][col]<<"\t";
    }
    //cout<<endl;
  }

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

TString RoundNumber(double n, int e, double d=1){
  if(d==0) return " - ";
  double neg = 1; if(n*d<0) neg = -1;
  n /= neg*d; n += 0.5*pow(10.,-e);
  int n_int = (int)n;
  int n_dec = (int)((1+n-n_int)*pow(10.,e));
  TString s_dec = ""; s_dec += n_dec; s_dec.Remove(0,1);
  TString result=""; 
  if(neg<0) result+="-";
  result+= n_int;
  if(e>0) {
    result+="."; result+=s_dec;
  }
  
  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<e-afterdot.Length(); i++)
    result += "0";
  return result;
}

// Cholesky-Banachiewicz algorithm
TMatrixT<double> Choleski(TMatrixT<double> A, int nRows){
  TMatrixT<double> L(nRows,nRows);
  
  for(int row=0; row<nRows; row++)
    for(int col=0; col<nRows; col++)
      L(row,col) = 0;
  
  for(int row=0; row<nRows; row++){
    for(int col=0; col<(row+1); col++){
      if(col==row) {
	for(int irow=0; irow<col; irow++) L(row,col) += L(col,irow)*L(col,irow);
	L(row,col) = sqrt(A(row,col) - L(row,col));
      } else {
	for(int irow=0; irow<col; irow++) L(row,col) += L(row,irow)*L(col,irow);
	L(row,col) = (A(row,col) - L(row,col))/L(col,col);
      }
    }
  }
  return L;
}


TH1F binHisto(TH2F *h2, int nbins, double xlow, double xhigh, double miny, double maxy, 
	      TString hname, TString var){
  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4;
  int nm2bin = h2->GetNbinsX(); 
  int nplbin = h2->GetNbinsY(); 
  double h2xmin = m2min, h2xmax = m2max, h2ymin = plmin, h2ymax = plmax;
  int nh2xbin = nm2bin, nh2ybin = nplbin;
  if(var=="candPstarLep"){
    h2xmin = plmin, h2xmax = plmax; 
    h2ymin = m2min, h2ymax = m2max; 
    nh2xbin = nplbin; nh2ybin = nm2bin;}
  double h2width = (h2xmax-h2xmin)/(double)nh2xbin;
  double h1width = (xhigh-xlow)/(double)nbins;
  double h2yini = (miny-h2ymin)*(double)nh2ybin/(h2ymax-h2ymin)+1;
  double h2yfin = (maxy-h2ymin)*(double)nh2ybin/(h2ymax-h2ymin)+1;
  int h2yini_i = (int)h2yini, h2yfin_i = (int)h2yfin; 
  double iniyfrac = h2yini_i+1 - h2yini, finyfrac = h2yfin - h2yfin_i;
  if(h2yini_i<1) h2yini_i = 1;	        
  if(h2yfin_i>nh2ybin) h2yfin_i = nh2ybin;

  TH1F h1(hname,"",nbins, xlow, xhigh);
  for(int bin = 1; bin <= nbins; bin++){
    double binvar = 0, tempvar = 0;
    double binmin = xlow+(double)(bin-1)*h1width;
    double h2xini = (binmin-h2xmin)/h2width+1, h2xfin = (binmin+h1width-h2xmin)/h2width+1;
    int h2xini_i = (int)h2xini, h2xfin_i = (int)h2xfin; 
    double inifrac = h2xini_i+1 - h2xini, finfrac = h2xfin - h2xfin_i;
    if(h2xini_i<1) h2xini_i = 1;	        
    if(h2xfin_i>nh2xbin) h2xfin_i = nh2xbin;
    if(h2xini_i==h2xfin_i){inifrac = h1width/h2width; finfrac = 1.;}
    for(int h2xbin = h2xini_i; h2xbin<h2xfin_i+1; h2xbin++){
      for(int h2ybin = h2yini_i; h2ybin < h2yfin_i+1; h2ybin++){
	if(var=="candPstarLep") tempvar = h2->GetBinContent(h2ybin, h2xbin);
	else tempvar = h2->GetBinContent(h2xbin, h2ybin);
	if(h2xbin == h2xini_i) tempvar *= inifrac;
	if(h2xbin == h2xfin_i) tempvar *= finfrac;
	if(h2ybin == h2yini_i) tempvar *= iniyfrac;
	if(h2ybin == h2yfin_i) tempvar *= finyfrac;
	binvar += tempvar;
      }
    }
    h1.SetBinContent(bin,binvar);
    h1.SetBinError(bin,0);
  }

  return h1;
}


