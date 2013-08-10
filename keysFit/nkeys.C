// #include "TString.h"
// #include "TFile.h"
// #include "TH1F.h"
// #include "TH2F.h"
// #include "TCanvas.h"
// #include "TChain.h"
// #include "TStyle.h"
// #include "TSystem.h"
// #include "TRandom3.h"
// #include "RooRealVar.h"
// #include "RooAddPdf.h"
// #include "RooAbsReal.h"
// #include "RooDataSet.h"
// //#include "RooNDKeysPdf.h"
// #include <fstream>
// #include <iostream>

// using namespace std;
// using namespace RooFit;
// using std::cout;
// using std::endl;

void nkeys(TString sam="0", TString smooth = "100", TString adaptive = "av",
	  int nM2bin = 1000, int nPlbin = 300, int Times = 1) {
  gSystem->Load("libHtml");
  gSystem->Load("libMinuit");
  gSystem->Load("libRooFitCore.so");
  gSystem->Load("libRooFitModels.so");
  using namespace RooFit;
  time_t start,end;
  time (&start);
  double dif;

  int Sam = sam.Atoi();
  TString Base = "ND"; Base += adaptive; Base += "_"; Base += sam; Base += "_"; 
  Base += nM2bin*1000+nPlbin; Base += "_"; Base += smooth; 

  TString inputfile = "fitSamples/pdfSample"; inputfile += sam; inputfile += ".root";
  cout << "File = " << inputfile << endl;	
  TChain c("ntp1");
  c.Add(inputfile);

  Int_t MCType,MCSubmode,MCDssmode,MCD,MCPions,MCCombB,MCCombDs,MCDoubleSL,
    candTruLep,candDstarType,isBzero,isSP6,MCTaumode,trueLepCharge;
  Float_t candM2,candPstarLep;
  Float_t truePPi0,trueDssPPi0,CTL,CTV,Chi,Q2,trueDmass;
//   c.SetBranchAddress("MCType",&MCType);
//   c.SetBranchAddress("MCSubmode",&MCSubmode);
//   c.SetBranchAddress("MCDssmode",&MCDssmode);
//   c.SetBranchAddress("MCD",&MCD);
//   c.SetBranchAddress("MCPions",&MCPions);
//   c.SetBranchAddress("MCCombB",&MCCombB);
//   c.SetBranchAddress("MCCombDs",&MCCombDs);
//   c.SetBranchAddress("MCDoubleSL",&MCDoubleSL);
//   c.SetBranchAddress("candLepTru",&candTruLep);
//   c.SetBranchAddress("candDstarType",&candDstarType);
//   c.SetBranchAddress("isBzero",&isBzero);
//   c.SetBranchAddress("isSP6",&isSP6);
//   c.SetBranchAddress("MCTaumode",&MCTaumode);
//   c.SetBranchAddress("truePPi0",&truePPi0);
//   c.SetBranchAddress("trueDssPPi0",&trueDssPPi0);
//   c.SetBranchAddress("trueCTL",&CTL);
//   c.SetBranchAddress("trueCTV",&CTV);
//   c.SetBranchAddress("trueChi",&Chi);
//   c.SetBranchAddress("trueQ2",&Q2);
//   c.SetBranchAddress("trueLepCharge",&trueLepCharge);
//   c.SetBranchAddress("trueDmass",&trueDmass);
  c.SetBranchAddress("candM2",&candM2);
  c.SetBranchAddress("candPstarLep",&candPstarLep);
  TRandom3 rand; 
  TCanvas mm("mm","KEYS fits to mmiss-pstarl",1200,800);
  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4;
  double xlow = m2min,xhigh = m2max, ylow = plmin,yhigh = plmax;
  Int_t nbinx = 80,nbiny = 80;
  if (Sam==0 || Sam==2 || Sam==10 || Sam==12 || Sam == 20 || Sam == 23 || Sam == 26 || Sam == 29) {
    xlow = -2; xhigh = 4;
    if (Sam > 12) {nbinx = 40; nbiny = 40;}
  } else if (Sam==1 || Sam==11) { xlow = -2; xhigh = 6;}
  if (Sam==6 || Sam==7 || Sam==16 || Sam==17) {nbinx = 40; nbiny = 40;}
  if (Sam==8 || Sam==18) {xhigh = 4; nbinx = 40; nbiny = 40;}
  if (Sam==9 || Sam==19) {nbinx = 40; nbiny = 40;}
  if (Sam==21 || Sam==22 || Sam==24 || Sam==25 || Sam==27 || Sam==28 || Sam==30 || Sam==31) {
    xhigh = 8; nbinx = 40; nbiny = 20;}
  if (Sam > 31) {nbinx = 40; nbiny = 20;}
  double entries = c.GetEntries(), rejected = 0;Double_t hIntegral=-99.;
  TH2F *h2;
  h2 = new TH2F("h2","KEYS",nM2bin,m2min,m2max,nPlbin,plmin,plmax);
  TString hname = "AWG82/results/keys/root/";
  if(Times==1) hname += "single";
  else hname += "Times";
  hname += "/hKeys"; hname += Base; hname += ".root";
  TFile* hfile = new TFile(hname,"RECREATE"); 
  hfile->Close(); hfile->Delete();
  //m2min=xlow; m2max = xhigh;
  for(int rep = 0; rep < Times; rep++){
    RooRealVar mmiss2("candM2","candM2",m2min,m2max);
    RooRealVar pstarl("candPstarLep","candPstarLep",plmin,plmax);
    RooArgSet myVars(mmiss2,pstarl);
    RooDataSet  data("data","data",myVars);
    if(rep%10==0) cout<<"Doing repetition "<<rep<<" of "<< Times <<endl; 
    for (int evt = 0 ; evt < entries; evt ++) {
      if(Times>1) {
	c.GetEvent((int)rand.Uniform(entries));
      } else c.GetEvent(evt);
      if(candM2<m2min || candM2>m2max || candPstarLep<plmin || candPstarLep>plmax) {
	rejected++; continue;}
      mmiss2.setVal(candM2);
      pstarl.setVal(candPstarLep);
      data.add(RooArgSet(mmiss2,pstarl));
    }
    double smoo = smooth.Atof()/100.;
    RooNDKeysPdf DPpdf("DPpdf","DPpdf",RooArgList(mmiss2,pstarl),data,adaptive,smoo,4);
    time (&end);dif = difftime (end,start);
    if(rep%50==0) cout<<dif<<" seconds after finding the KEYS function"<<endl;
    time (&start);
    TString hNameTimes = "hTimes"; hNameTimes += rep;
    if(Times>1) h2 = new TH2F(hNameTimes,"KEYS",nM2bin,m2min,m2max,nPlbin,plmin,plmax);
    DPpdf.fillHistogram(h2,RooArgList(mmiss2,pstarl));
    hfile = new TFile(hname,"UPDATE"); 
    h2->Write();
    hfile->Close(); hfile->Delete();
    time (&end);dif = difftime (end,start);
    if(rep%50==0) cout<<dif<<" seconds after making histogram"<<endl;
    time (&start);
    if(Times>1) h2->Delete();
  }
  cout<<hname<<" saved"<<endl; 
  if(Times>1) {return;}
  entries -= rejected;

  TString M2titles[] = {"0 < p*_{l} < 1 GeV","1 < p*_{l} < 1.4 GeV","1.4 < p*_{l} < 1.8 GeV",
		      "1.8 < p*_{l} < 2.4 GeV","0 < p*_{l} < 2.4 GeV"};
  TString Pltitles[] = {"-4 < m^{2}_{miss} < 1 GeV^{2}","1 < m^{2}_{miss} < 12 GeV^{2}",
			"-4 < m^{2}_{miss} < 12 GeV^{2}"};
  TString M2cuts[] = {"candPstarLep<1","candPstarLep>1&&candPstarLep<1.4",
		      "candPstarLep>1.4&&candPstarLep<1.8","candPstarLep>1.8&&candPstarLep<2.4", ""};
  TString Plcuts[] = {"candM2<1","candM2>=1",""};
  double PlLimits[] = {0, 1, 1.4, 1.8, 2.4};
  double M2Limits[] = {0., 5., 16.};
  int binlim[5]; int plbinlim[3];
  for(int i=0;i<5;i++) {
    PlLimits[i] = PlLimits[i]/2.4*(double)nPlbin;
    binlim[i] = (int)PlLimits[i];
    if(i<3){
      M2Limits[i] = M2Limits[i]/16.*(double)nM2bin;
      plbinlim[i] = (int)M2Limits[i];
    }
  }

  TString psname = "AWG82/results/keys/eps/epsKeys"; psname += Base; psname += ".ps";
  mm.Print(psname+"[");
  TH1F *hm2[5], *m2[5], *hpl[3], *pl[3];
  TString M2names[5], Plnames[3];
  for(int i=0;i<5;i++){
    M2names[i] = "hm2_"; M2names[i] += i;
    hm2[i] = new TH1F(M2names[i],M2titles[i],nM2bin,m2min,m2max); 
    if(i<3) {
      Plnames[i] = "hpl_"; Plnames[i] += i;
      hpl[i] = new TH1F(Plnames[i],Pltitles[i],nPlbin,plmin,plmax); 
    }
  }
  
  hIntegral = h2->Integral();
  for(int i=0;i<5;i++){
    TString hname = "m2"; hname += i;
    TString vari = "candM2>>"; vari+=hname; vari+="("; vari+= nbinx; vari+=",";vari+= xlow; 
    vari+=",";vari+= xhigh; vari+=")";
    c.Draw(vari,M2cuts[i]);
    m2[i] = (TH1F*)gDirectory->Get(hname);
    m2[i]->SetXTitle("m^{2}_{miss} [GeV^{2}]");
    m2[i]->SetTitle(M2titles[i]);
    m2[i]->Sumw2();
    m2[i]->SetMarkerStyle(20);
    m2[i]->SetMarkerSize(.6);
    gStyle->SetOptStat(0);
    if(i<4){
      for(int j=1; j<nM2bin+1; j++){
	double binVal = 0;
	for(int binp = binlim[i]+1; binp < binlim[i+1]+1; binp++){
	  binVal+=h2->GetBinContent(j,binp)*entries*nM2bin*(xhigh-xlow)/nbinx/(m2max-m2min)/hIntegral;
	}
	hm2[i]->SetBinContent(j,binVal);
      }
    } 
    hm2[i]->SetLineColor(4);
    hm2[i]->SetLineWidth(1);
    if(i<4) hm2[4]->Add(hm2[i]);
  }
  for(int i=0;i<3;i++){
    TString hname = "pl"; hname += i;
    TString vari = "candPstarLep>>"; vari+=hname; vari+="("; vari+= nbiny; vari+=",";vari+= ylow; 
    vari+=",";vari+= yhigh; vari+=")";
    c.Draw(vari,Plcuts[i]);
    pl[i] = (TH1F*)gDirectory->Get(hname);
    pl[i]->SetXTitle("p*_{l} [GeV]");
    pl[i]->SetTitle(Pltitles[i]);
    pl[i]->Sumw2();
    pl[i]->SetMarkerStyle(20);
    pl[i]->SetMarkerSize(.6);
    gStyle->SetOptStat(0);
    if(i<2){
      for(int j=1; j<nPlbin+1; j++){
	double binVal = 0;
	for(int binp = plbinlim[i]+1; binp < plbinlim[i+1]+1; binp++){
	  binVal += h2->GetBinContent(binp,j)*entries*nPlbin/nbiny/hIntegral;	
	}
	hpl[i]->SetBinContent(j,binVal);
      }
    }
    hpl[i]->SetLineColor(4);
    hpl[i]->SetLineWidth(1);
    if(i<2) hpl[2]->Add(hpl[i]);
  }
  m2[4]->Draw("e0"); m2[4]->Draw("e1 same"); hm2[4]->Draw("c same"); mm.Print(psname);
  pl[2]->Draw("e0"); pl[2]->Draw("e1 same"); hpl[2]->Draw("c same"); mm.Print(psname);
  for(int i=0;i<4;i++){
    m2[i]->Draw("e0"); m2[i]->Draw("e1 same"); hm2[i]->Draw("c same"); mm.Print(psname);
  }
  for(int i=0;i<2;i++){
    pl[i]->Draw("e0"); pl[i]->Draw("e1 same"); hpl[i]->Draw("c same"); mm.Print(psname);
  }
  mm.Print(psname+"]");
  TCanvas all6("all6","KEYS fits to mmiss-pstarl",2000,2400);
  all6.cd();
  all6.Divide(2,3);
  all6.cd(1);
  m2[4]->Draw("e0"); m2[4]->Draw("e1 same"); hm2[4]->Draw("c same");
  all6.cd(2);
  pl[2]->Draw("e0"); pl[2]->Draw("e1 same"); hpl[2]->Draw("c same");
  for(int i=0;i<4;i++){
    all6.cd(i+3);
    m2[i]->Draw("e0"); m2[i]->Draw("e1 same"); hm2[i]->Draw("c same"); 
  }
  TString epsname = "AWG82/results/keys/eps/EpsKeys";
  epsname += Base; epsname += ".eps";
  all6.SaveAs(epsname);
  time (&end);dif = difftime (end,start);
  h2->Delete();
  cout<<dif<<" seconds after plotting data. Written "<<psname<<endl;
  return; 

} 

