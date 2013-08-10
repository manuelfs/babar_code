// Generates random numbers distributed as an n-dimensional Gaussian, including correlations
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TMatrixTSym.h"
#include "TF3.h"

#define numberFF 3
TMatrixTSym<double> CovFF(numberFF); //CovFF will be both the covariance matrix, and its inverse
double valFF[3][3] = {{1.16,1.37,0.845}, {0.05, 0.06, 0.038},{0.71, -0.84, -0.83}};
Double_t Gaus3(Double_t *x, Double_t *par);

void GenCorrGaus(int nPoints = 10000, double corr=0.71){
  valFF[2][0] = corr;
  for(int iFF=0; iFF<numberFF; iFF++) {
    CovFF(iFF,iFF) = pow(valFF[1][iFF],2);
    CovFF(iFF,(iFF+1)%3) = valFF[2][iFF]*valFF[1][iFF]*valFF[1][(iFF+1)%3];
    CovFF((iFF+1)%3,iFF) = CovFF(iFF,(iFF+1)%3);
  }
  CovFF.Invert(); 
  TF3 FFGaus("FFGaus", Gaus3, 0,2, 0,2, 0,2, 1);
  FFGaus.SetParameter(0,1);
  FFGaus.SetNpx(200); FFGaus.SetNpy(200); FFGaus.SetNpz(200); // Setting more integration points
  TCanvas can("can","",600,300);can.Divide(2,1);
  int nbins = 200;
  double minX = 0.9, maxX = 1.6;
  TH2F h2("h2","", nbins, minX, maxX, nbins, minX, maxX);
  TH1F h1("h1","", nbins, minX, maxX);
  double x=0,y=0,z=0; 
  for(int point = 1; point<=nPoints; point++) {
    FFGaus.GetRandom3(x,y,z);
    h2.Fill(x,y);
    h1.Fill(x);
  }
  can.cd(1); h2.Draw("contz");
  can.cd(2); h1.Draw();
  can.SaveAs("FFVariations.eps");
    
}

// 3-dimensional Gaussian with inverse covariance matrix CovFF
Double_t Gaus3(Double_t *x, Double_t *par){
  double y[3], total = 0, vCov[] = {0,0,0};

  for(int iFF=0; iFF<numberFF; iFF++) {
    y[iFF] = x[iFF]-valFF[0][iFF];
    for(int jFF=0; jFF<numberFF; jFF++) 
      vCov[iFF] += y[jFF]*CovFF(iFF,jFF);
  }
  for(int iFF=0; iFF<numberFF; iFF++) total += vCov[iFF]*y[iFF];
  return par[0]*exp(-total/2);
}


