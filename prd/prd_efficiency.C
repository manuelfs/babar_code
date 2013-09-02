#include "babar_code/PlotsThesis/PlotUtils.cc"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TChain.h"
#include "TStyle.h"
#include "TString.h"
#include "TRandom3.h"
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;
using std::cout;
using std::endl;

int findMode(double x, double cfmode[], int Nmodes) {
  int minMode=0, maxMode=Nmodes;
  int midMode = (maxMode-minMode)/2+minMode, rep=0;
  while(maxMode>minMode+1 && rep <20){
    // cout<<x<<" is between cfmode["<<minMode<<"] = "<<cfmode[minMode]
    // 	<<" and cfmode["<<maxMode<<"] = "<<cfmode[maxMode]
    // 	<<", with cfmode["<<midMode<<"] = "<<cfmode[midMode]<<endl;
    if(x<cfmode[midMode]) maxMode = midMode;
    else minMode = midMode;
    midMode = (maxMode-minMode)/2+minMode;
    rep++;
  }
  return minMode;
}

void prd_efficiency(long Nevents=1200, long Nrep=10, double effmin=0.0001, double effmax=0.8, 
		    int seed=1){
  TString filename = "babar_code/prd/fraction_eff.txt";
  ifstream infile(filename);
  TRandom3 rand(seed); // 1: first seed constant - 0: time dependent seed.
  double fmode[1000], cfmode[1000], eff[1000], f_eff, norm=0, AverEff=0;
  int Nmodes=0, imode;
  while(infile>> f_eff){
    eff[Nmodes] = rand.Uniform(effmin, effmax);
    //eff[Nmodes] = rand.Exp(effmax);
    //eff[Nmodes] = rand.Gaus(effmin,effmax);
    if(eff[Nmodes]<=0) eff[Nmodes] = 1e-6;
    if(eff[Nmodes]>1) eff[Nmodes] = 1.;
    //fmode[Nmodes] = f_eff/eff[Nmodes];
    fmode[Nmodes] = f_eff;
    norm += fmode[Nmodes];
    Nmodes++;
  }
  cfmode[0] = 0;
  for(int mode=0; mode<Nmodes; mode++) {
    fmode[mode] /= norm;
    cfmode[mode+1] = fmode[mode] + cfmode[mode];
    AverEff += fmode[mode]*eff[mode];
  }

  double Ntotal = Nevents/AverEff, nSig = 4.3;
  double sig_paper = sqrt(Nevents);
  int nBins = 60, nperBin = 2*nSig*sig_paper/nBins;
  int minX = Nevents-nSig*sig_paper+(2*nSig*sig_paper-nperBin*nBins)/2, maxX = minX+nperBin*nBins;
  if(nperBin<1){
    minX = Nevents-nBins/2; maxX = Nevents+nBins/2; nperBin = 1;
  }
  TString title = "Entries/"; title += nperBin;
  gStyle->SetOptStat(0);
  TCanvas can("can","Efficiencies",1400,800); 
  TH1F hN("hN","",nBins, minX, maxX);
  hN.SetXTitle("Simulated N_{sel}");
  hN.SetYTitle(title);

  double mean = 0, rms = 0;
  for(int rep=0; rep<Nrep; rep++){
    //Ntotal = rand.PoissonD(Nevents)/AverEff;
    double Npass=0;
    for(int eve=0; eve < Ntotal; eve++) {
      imode = findMode(rand.Uniform(1), cfmode, Nmodes);
      if(rand.Uniform(1) < eff[imode]) Npass++;
    }
    mean += Npass; rms += Npass*Npass;
    hN.Fill(Npass);
  }
  double DNrep = (double)Nrep;
  mean /= DNrep; rms = sqrt((rms - mean*mean*DNrep)/(DNrep-1));

  int digits = 2;
  double sig_sim = rms;
  cout<<"Paper: "<<RoundNumber(Nevents,0)<<" +- "<<RoundNumber(sig_paper,digits)
      <<"  -  Simulation: "<<RoundNumber(mean, digits)<<" +- "<<RoundNumber(sig_sim, digits)
      <<"  \t  Ratio: "<< RoundNumber(sig_paper,4,sig_sim)<<endl;

  hN.SetLineColor(4); hN.SetLineWidth(2);
  hN.SetMaximum(1.15*hN.GetMaximum());
  hN.Draw("");

  TLatex label; label.SetNDC(kTRUE); label.SetTextFont(132); label.SetTextSize(0.036);
  title = "#bar{N_{sel}} = ";title += Nevents; title += ", #epsilon = ";
  title += RoundNumber(AverEff*100,3); title += "%  #Rightarrow  N_{gen} = "; 
  title += RoundNumber(Ntotal,0);
  label.DrawLatex(0.129, 0.84, title); 

  title = "#sqrt{#bar{N_{sel}}} = ";  title += RoundNumber(sig_paper,1); 
  label.DrawLatex(0.129, 0.76, title); 

  title = "#sigma(N_{sel}) = ";  title += RoundNumber(sig_sim,1); 
  label.DrawLatex(0.128, 0.70, title); 

  title = "#Rightarrow  #frac{#sqrt{#bar{N_{sel}}}}{#sigma(N_{sel})} = "; 
  title += RoundNumber(sig_paper,3,sig_sim);
  label.DrawLatex(0.235, 0.74, title); 

  can.SaveAs("babar_code/prd/plot_efficiency.gif"); 

}

void calc_fractions(){
  TChain chain("ntp1");
  chain.Add("AWG82/ntuples/small/RAll_RunAll.root");

  int candBMode;
  chain.SetBranchAddress("candBMode",&candBMode);

  double events[10000];
  for(int eve=0; eve < 10000; eve++) events[eve] = 0;
  for(int eve=0; eve < chain.GetEntries(); eve++){
    chain.GetEvent(eve);
    events[candBMode-10000]++;
  }

  vector<Double_t> effModes;
  TString filename = "babar_code/prd/fraction_eff.txt";
  ofstream outfile(filename);
  double total = 0;
  for(int eve=0; eve < 10000; eve++) total += events[eve];
  for(int eve=0; eve < 10000; eve++) {
    if(events[eve] > 0) {
      double effi = events[eve]*100/total;     
      effModes.push_back(effi);
    }
  }
  sort(effModes.rbegin(),effModes.rend()); // Sort in descending order
  gStyle->SetOptStat(0);
  TCanvas can("can","Efficiencies",1400,800); 
  TH1F hN("hN","", effModes.size(), 0, effModes.size());
  hN.SetXTitle("B_{tag} decay mode");
  hN.SetYTitle("Abundance (%)");

  for(unsigned int mode=0; mode < effModes.size(); mode++) {
    outfile << RoundNumber(effModes[mode], 5) << " ";
    hN.SetBinContent(mode+1,effModes[mode]);
    if((mode+1)%10 == 0) outfile<<endl;
  }
  cout<<"Printed "<<effModes.size()<<" in "<<filename<<endl;
  outfile.close();

  hN.SetLineColor(2); hN.SetLineWidth(2);
  hN.Draw("");
  can.SaveAs("babar_code/prd/plot_abundance.gif"); 

}

