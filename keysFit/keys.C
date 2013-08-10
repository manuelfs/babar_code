//#include "Roo2DKeysPdf.h"
void rotate(double &mrot, double &prot, double r11, double r12, double Xmean, double Ymean, double m2, double pl);

void keys(TString sam="0", TString smooth = "100", int Times = 1, TString doRot = "NoR", TString adaptive = "av",
	  int nM2bin = 200, int nPlbin = 200) {
  gSystem->Load("libHtml");
  gSystem->Load("libMinuit");
  gSystem->Load("libRooFitCore.so");
  gSystem->Load("libRooFitModels.so");
  using namespace RooFit;
  time_t start,end;
  time (&start);
  double dif;

  int Sam = sam.Atoi();
  //if(Sam==0 || Sam==2 || Sam==10 || Sam==11 || Sam==12) doRot = "Rot";
  TString Base = adaptive; Base += doRot; Base += sam; Base += "_"; Base += nM2bin*1000+nPlbin; 
  Base += "_"; Base += smooth; 

  TString inputfile = "fitSamples/pdfSample"; inputfile += sam; inputfile += ".root";
  cout << "File = " << inputfile << endl;	
  TChain c("ntp1");
  c.Add(inputfile);

  Int_t MCType,MCSubmode,MCDssmode,MCD,MCPions,MCCombB,MCCombDs,MCDoubleSL,
    candTruLep,candDstarType,isBzero,isSP6,MCTaumode,trueLepCharge;
  Float_t candM2,candPstarLep;
  Float_t truePPi0,trueDssPPi0,CTL,CTV,Chi,Q2,trueDmass;
  c.SetBranchAddress("MCType",&MCType);
  c.SetBranchAddress("MCSubmode",&MCSubmode);
  c.SetBranchAddress("MCDssmode",&MCDssmode);
  c.SetBranchAddress("MCD",&MCD);
  c.SetBranchAddress("MCPions",&MCPions);
  c.SetBranchAddress("MCCombB",&MCCombB);
  c.SetBranchAddress("MCCombDs",&MCCombDs);
  c.SetBranchAddress("MCDoubleSL",&MCDoubleSL);
  c.SetBranchAddress("candLepTru",&candTruLep);
  c.SetBranchAddress("candDstarType",&candDstarType);
  c.SetBranchAddress("isBzero",&isBzero);
  c.SetBranchAddress("isSP6",&isSP6);
  c.SetBranchAddress("MCTaumode",&MCTaumode);
  c.SetBranchAddress("truePPi0",&truePPi0);
  c.SetBranchAddress("trueDssPPi0",&trueDssPPi0);
  c.SetBranchAddress("trueCTL",&CTL);
  c.SetBranchAddress("trueCTV",&CTV);
  c.SetBranchAddress("trueChi",&Chi);
  c.SetBranchAddress("trueQ2",&Q2);
  c.SetBranchAddress("trueLepCharge",&trueLepCharge);
  c.SetBranchAddress("trueDmass",&trueDmass);
  c.SetBranchAddress("candM2",&candM2);
  c.SetBranchAddress("candPstarLep",&candPstarLep);
  TRandom3 rand; int ran;
  TCanvas mm("mm","KEYS fits to mmiss-pstarl",1200,800);
  gStyle->SetPalette(1);
  double All = 5;
  int nbins = 120;
  TH2F rotated("rotated","Rotated m^{2}_{miss}-p*_{l}",nbins,-All,All,nbins,-All,All);
  TH2F ori("ori","Original m^{2}_{miss}-p*_{l}",nbins,-All,All,nbins,-All,All);
  double r11, r12, Xmean, Ymean;
  double mmp, plp, mmax=0,mmin=0,pmax=0,pmin=0;
  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4;
  if(doRot == "Rot"){
    c.Draw("candPstarLep:candM2>>cov(800,-4,12,800,0,2.4)","","contz");
    TH2F *totcov = (TH2F*)gDirectory->Get("cov");
    double xx = totcov->GetRMS(1); xx = xx*xx;
    double yy = totcov->GetRMS(2); yy = yy*yy;
    double xy = totcov->GetCovariance();
    Xmean = totcov->GetMean(1);
    Ymean = totcov->GetMean(2);
    double lambda = (-sqrt(xx*xx-2*xx*yy+4*xy*xy+yy*yy)+xx+yy)/2;
    double lambda2 = (sqrt(xx*xx-2*xx*yy+4*xy*xy+yy*yy)+xx+yy)/2;
    if(lambda2>lambda) lambda = lambda2;
    r11 = (lambda-yy)/xy;
    r12 = -1/sqrt(r11*r11+1);
    r11 = -r11/sqrt(r11*r11+1);
    if(r12*r11>0&&r12<0 || r12*r11<0&&r11<0){
      r12 = -r12;
      r11 = -r11;
    }
    cout<<"RMSx "<<xx<<", RMSy "<<yy<<", lambda "<<lambda<<" and covariance "<<xy<<endl;
    rotate(mmp,plp,r11,r12,Xmean,Ymean,m2min,0);
    if(mmp<mmin)mmin=mmp;  if(plp<pmin)pmin=plp;
    if(mmp>mmax)mmax=mmp;  if(plp>pmax)pmax=plp;
    rotate(mmp,plp,r11,r12,Xmean,Ymean,m2min,2.4);
    if(mmp<mmin)mmin=mmp;  if(plp<pmin)pmin=plp;
    if(mmp>mmax)mmax=mmp;  if(plp>pmax)pmax=plp;
    rotate(mmp,plp,r11,r12,Xmean,Ymean,m2max,0);
    if(mmp<mmin)mmin=mmp;  if(plp<pmin)pmin=plp;
    if(mmp>mmax)mmax=mmp;  if(plp>pmax)pmax=plp;
    rotate(mmp,plp,r11,r12,Xmean,Ymean,m2max,2.4);
    if(mmp<mmin)mmin=mmp;  if(plp<pmin)pmin=plp;
    if(mmp>mmax)mmax=mmp;  if(plp>pmax)pmax=plp;
  }else{
    mmax=m2max; mmin=m2min; pmax=plmax; pmin=plmin;
  }
  double entries = c.GetEntries();
  //entries = 4;
  //if(entries<10000) Times = 500;
  //else Times = 100;
  TH2F *h2;
  for(int rep = 0; rep < Times; rep++){
    RooRealVar mmiss2("candM2","candM2",mmin,mmax);
    RooRealVar pstarl("candPstarLep","candPstarLep",pmin,pmax);
    RooArgSet myVars(mmiss2,pstarl);
    RooDataSet  data("data","data",myVars);
    cout<<"Doing repetition "<<rep<<" of "<< Times <<endl; 
    for (int evt = 0 ; evt < entries; evt ++) {
      ran = rand.Uniform(entries);
      if(Times>1) c.GetEvent(ran);
      else c.GetEvent(evt);
      if(doRot == "Rot"){
	double Mx = candM2-Xmean, Py = candPstarLep-Ymean;
	mmp = r11*(Mx)+r12*(Py);
	plp = -r12*(Mx)+r11*(Py);
	ori.Fill(Mx,Py);
	if(rep==0) rotated.Fill(mmp,plp);
      }else{
	mmp = candM2;
	plp = candPstarLep;
      }
      mmiss2.setVal(mmp);
      pstarl.setVal(plp);
      //     if (MCType == 0)
      //       totWeight.setVal(myWM->getCombWeight(MCCombB,MCCombDs,MCDoubleSL,candTruLep));
      //     else
      //       totWeight.setVal(myWM->getEventWeight(candType,candDstarType,MCType,MCSubmode,MCDssmode,MCD,MCPions,
      // 					    isBzero,isSP6,MCTaumode,truePPi0,trueDmass,CTL,CTV,Chi,Q2,
      // 					    trueLepCharge,candM2));
      data.add(RooArgSet(mmiss2,pstarl));
    }
    //data.setWeightVar(totWeight);
    cout<<"Data added"<<endl;
    if(doRot == "Rot" && rep==0){
      TString rotName = "AWG82/results/keys/root/rotation/Rotation_"; rotName += Base; rotName += ".root";
      TFile* hRot = new TFile(rotName,"RECREATE"); 
      hRot->cd();
      ori.Write(); rotated.Write();
      hRot->Close(); hRot->Delete();
      cout<<"("<<r11<<", "<<r12<<") and covariance "<<rotated.GetCovariance()<<endl;
    }
    double smoo = smooth.Atof()/100.;
    Roo2DKeysPdf DPpdf("DPpdf","DPpdf",mmiss2,pstarl,data,adaptive,smoo);
    time (&end);dif = difftime (end,start);
    cout<<dif<<" seconds after finding the KEYS function"<<endl;
    time (&start);
    TH2F *h2r;
    h2 = new TH2F("h2","KEYS",nM2bin,m2min,m2max,nPlbin,plmin,plmax);
    if(doRot == "Rot"){
      h2r = new TH2F("h2r","KEYS rotated",nM2bin,mmin,mmax,nPlbin,pmin,pmax);
      DPpdf.fillHistogram(h2r,RooArgList(mmiss2,pstarl));
      for(int i=1; i<nM2bin+1; i++){
	for(int j=1; j<nPlbin+1; j++){
	  double M2 = (m2max-m2min)/nM2bin*(i-0.5)+m2min, Pl = (plmax-plmin)/nPlbin*(j-0.5)+plmin;
	  rotate(mmp,plp,r11,r12,Xmean,Ymean,M2,Pl);
	  int ibin = (mmp-mmin)*nM2bin/(mmax-mmin)+0.5;
	  int jbin = (plp-pmin)*nPlbin/(pmax-pmin)+0.5;
	  h2->SetBinContent(i,j,h2r->GetBinContent(ibin,jbin));
	}
      }
      h2r->Delete();
    }else{
      DPpdf.fillHistogram(h2,RooArgList(mmiss2,pstarl));
    } 
    TString hname = "AWG82/results/keys/root/";
    if(Times==1) hname += "single";
    else hname += "new";
    hname += "/hKeys"; hname += Base; 
    if(Times>1){
      hname += "_"; hname += rep;
    }
    hname += ".root";
    TFile* hfile = new TFile(hname,"RECREATE"); 
    h2->Write();
    hfile->Close();
    cout<<"KEYS histogram saved in "<<hname<<endl;
    //   RooDataHist* Rdh2 = new RooDataHist("Rdh2","KEYS",RooArgList(mmiss2,pstarl),h2);
    //   RooHistPdf* Rh2 = new RooHistPdf("Rh2","KEYS",RooArgList(mmiss2,pstarl),*Rdh2,2);
    time (&end);dif = difftime (end,start);
    cout<<dif<<" seconds after making histogram"<<endl;
    time (&start);
    if(Times>1) h2->Delete();;
  }
  if(Times>1) return;
  Float_t xlow,xhigh;
  Int_t nbinx,nbiny;
  xlow = m2min;
  xhigh = m2max;
  nbinx = 80;
  nbiny = 80;
  if (Sam==0 || Sam==2 || Sam==10 || Sam==12 || Sam == 20 || Sam == 23 || Sam == 26 || Sam == 29) {
    xlow = -2; xhigh = 4;
    if (Sam > 12) {nbinx = 40; nbiny = 40;}
  }
  else if (Sam==1 || Sam==11) {
    xlow = -2; xhigh = 6;
  }
  if (Sam==6 || Sam==7 || Sam==16 || Sam==17) {
    nbinx = 40; nbiny = 40;
  }
  if (Sam==8 || Sam==18) {
    xhigh = 4; nbinx = 40; nbiny = 40;
  }
  if (Sam==9 || Sam==19) {
    nbinx = 40; nbiny = 40;
  }
  if (Sam==21 || Sam==22 || Sam==24 || Sam==25 || Sam==27 || Sam==28 || Sam==30 || Sam==31) {
    xhigh = 8; nbinx = 40; nbiny = 20;
  }
  if (Sam > 31) {
    nbinx = 40; nbiny = 20;
  }

  TString M2titles[] = {"0 < p*_{l} < 1 GeV","1 < p*_{l} < 1.4 GeV","1.4 < p*_{l} < 1.8 GeV",
		      "1.8 < p*_{l} < 2.4 GeV","0 < p*_{l} < 2.4 GeV"};
  TString Pltitles[] = {"-4 < m^{2}_{miss} < 1 GeV^{2}","1 < m^{2}_{miss} < 12 GeV^{2}",
			"-4 < m^{2}_{miss} < 12 GeV^{2}"};
  TString M2cuts[] = {"candPstarLep<1","candPstarLep>1&&candPstarLep<1.4",
		      "candPstarLep>1.4&&candPstarLep<1.8","candPstarLep>1.8&&candPstarLep<2.4", ""};
  TString Plcuts[] = {"candM2<1","candM2>=1",""};
  double limits[] = {0, 1, 1.4, 1.8, 2.4};
  int binlim[5];
  for(int i=0;i<5;i++) binlim[i] = limits[i]/(plmax-plmin)*nPlbin;

  TString psname = "AWG82/results/keys/eps/epsKeys";
  psname += Base; psname += ".ps";
  double tot = h2->Integral();
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
    m2[i]->SetMarkerSize(1);
    gStyle->SetOptStat(0);
    if(i<4){
      for(int j=1; j<nM2bin+1; j++){
	double binVal = 0;
	for(int binp = binlim[i]+1; binp < binlim[i+1]+1; binp++){
	  binVal += h2->GetBinContent(j,binp)*entries*nM2bin*(xhigh-xlow)/nbinx/16/tot;	
	}
	hm2[i]->SetBinContent(j,binVal);
      }
    } 
    hm2[i]->SetLineColor(4);
    hm2[i]->SetLineWidth(2);
    if(i<4) hm2[4]->Add(hm2[i]);
  }
  int plbinlim[3] = {0,nM2bin*5/16,nM2bin};
  for(int i=0;i<3;i++){
    TString hname = "pl"; hname += i;
    TString vari = "candPstarLep>>"; vari+=hname; vari+="("; vari+= nbiny; vari+=",0,2.4)";
    c.Draw(vari,Plcuts[i]);
    pl[i] = (TH1F*)gDirectory->Get(hname);
    pl[i]->SetXTitle("p*_{l} [GeV]");
    pl[i]->SetTitle(Pltitles[i]);
    pl[i]->Sumw2();
    pl[i]->SetMarkerStyle(20);
    pl[i]->SetMarkerSize(1);
    gStyle->SetOptStat(0);
    if(i<2){
      for(int j=1; j<nPlbin+1; j++){
	double binVal = 0;
	for(int binp = plbinlim[i]+1; binp < plbinlim[i+1]+1; binp++){
	  binVal += h2->GetBinContent(binp,j)*entries*nPlbin/nbiny/tot;	
	}
	hpl[i]->SetBinContent(j,binVal);
      }
    }
    hpl[i]->SetLineColor(4);
    hpl[i]->SetLineWidth(2);
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
  time (&end);dif = difftime (end,start);
  h2->Delete();
  cout<<dif<<" seconds after plotting data. Written "<<psname<<endl;
  return; 

} 


void rotate(double &mrot, double &prot, double r11, double r12, double Xmean, double Ymean, double m2, double pl){
    double Mx = m2-Xmean, Py = pl-Ymean;
    mrot = r11*(Mx)+r12*(Py);
    prot = -r12*(Mx)+r11*(Py);

}
