#include "TString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include <fstream>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;

Double_t alphas = 0.26;
Double_t realz = 0.29, epsc = 0.35845, epsb = 0.10395;  // CLN paper
//Double_t realz = 0.2145, epsc = 0.5549, epsb = 0.119;  // MSbar, mb scale
Double_t lam_bar = 0.48;
static const double pi = 3.14159265;
static const double mB  = 5.2792, mDs = 2.00696, mD  = 1.8648; // Charged B
double RDs = 2.*sqrt(mB*mDs)/(mB+mDs);

Double_t r(Double_t w){
  Double_t ans = log(w + sqrt(w*w-1));
  return ans / sqrt(w*w-1);
}

Double_t H(Double_t w,Double_t z){
  Double_t ans;
  ans = 2.*(w-1.)*z*(1.+z)*log(z);
  ans -= (w+1. - 2.*w*(2.*w+1.)*z + (5.*w+2.*w*w-1.)*z*z - 2.*z*z*z) * r(w);
  ans *= z / (1.-2.*w*z+z*z) / (1.-2.*w*z+z*z);
  ans += (z*(1.-log(z)-z))/(1.-2.*w*z+z*z);
  return ans;
}

Double_t H5(Double_t w,Double_t z){
  Double_t ans;
  ans = 2.*(w+1.)*z*(1-z)*log(z);
  ans -= (w-1. - 2.*w*(2.*w-1.)*z + (5.*w-2.*w*w+1.)*z*z - 2.*z*z*z) * r(w);
  ans *= z / (1.-2.*w*z+z*z) / (1.-2.*w*z+z*z);
  ans += (z*(1.-log(z)+z))/(1.-2.*w*z+z*z);
  return ans;
}

Double_t C1(){
  return 1.4;
}

Double_t C15(Double_t w){
  return C1() * (1. - 4.*alphas*r(w)/(3*pi));
}

Double_t C25(Double_t w){
  return C1() * -2.*alphas*H5(w,1./realz)/(3*pi);
}

Double_t C35(Double_t w){
  return C1() * 2.*alphas*H5(w,realz)/(3*pi);
}

Double_t L1(Double_t w){
  return 0.72*(w-1)*lam_bar;
}

Double_t L2(Double_t w){
  return -0.16*(w-1)*lam_bar;
}

Double_t L3(Double_t w){
  return -0.24*lam_bar;
}

Double_t L4(Double_t w){
  return 0.24*lam_bar;
}

Double_t L5(Double_t w){
  return -1.*lam_bar;
}

Double_t L6(Double_t w){
  return -3.24*lam_bar/(w+1);
}

Double_t hv(Double_t w){
  return C1() + epsc*(L2(w)-L5(w)) + epsb*(L1(w)-L4(w));
}

Double_t ha1(Double_t w){
  return C15(w) + epsc*(L2(w) - (w-1)*L5(w)/(w+1)) + epsb*(L1(w) - (w-1)*L4(w)/(w+1));
}

Double_t ha2(Double_t w){
  return C25(w) + epsc*(L3(w)+L6(w));
}

Double_t ha3(Double_t w){
  return C15(w) + C35(w) + epsc*(L2(w)-L3(w)-L5(w)+L6(w)) + epsb*(L1(w)-L4(w));
}

Double_t A3(Double_t w){
  double Term1 = mB/(mB+mDs)*(w+1)*ha1(w);
  double Term2 = (mB-mDs)/(2*mDs)*(ha3(w)+mDs/mB*ha2(w));
  return (mB+mDs)/(2*sqrt(mB*mDs))*(Term1-Term2);
}

Double_t A0(Double_t w){
  double q2 = mB*mB+mDs*mDs-w*2*mB*mDs;
  return A3(w)+q2/(4*mB*mDs)*sqrt(mB/mDs)*(ha3(w)-mDs/mB*ha2(w));
}

Double_t R0(Double_t w){
  return A0(w)*RDs/ha1(w);
}

Double_t R1(Double_t w){
  return hv(w)/ha1(w);
}

Double_t R2(Double_t w){
  return (ha3(w) + mDs/mB*ha2(w))/ha1(w);
}

Double_t R4(Double_t w){
  return (ha3(w) - (2.007/5.28)*ha2(w))/ha1(w);
}

Double_t Rtest(Double_t w){
  return ha2(w) / ha1(w);
}

void FF(double lam=0.48, double al=0.26){

  alphas = al;
  lam_bar = lam;

  TCanvas can;

  int nBins = 1000;
  TH1F *h = new TH1F("h","",nBins,0,0.5);
  for (int i=1 ; i<=nBins ; i++)
    h->SetBinContent(i,R0(1+h->GetBinCenter(i)));
  h->Fit("pol2");
  h->SetLineColor(kRed);
  h->SetLineWidth(4);
  h->GetXaxis()->SetTitle("w-1");
  h->GetYaxis()->SetTitle("R_{0}");
  h->GetXaxis()->SetNdivisions(205);
  h->GetYaxis()->SetNdivisions(205);
  h->Draw();
  can.SaveAs("babar_code/FF/eps/R2.eps");
}
