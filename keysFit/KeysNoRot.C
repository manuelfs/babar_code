//#include "Roo2DKeysPdf.h"

void keys(TString sam="0", TString smooth = "100") {
  gSystem->Load("libHtml");
  gSystem->Load("libMinuit");
  gSystem->Load("libRooFitCore.so");
  gSystem->Load("libRooFitModels.so");
  using namespace RooFit;
  time_t start,end;
  time (&start);
  double dif;
  RooRealVar mmiss2("candM2","candM2",-4,12);
  RooRealVar pstarl("candPstarLep","candPstarLep",0.,2.4);

  RooArgSet myVars(mmiss2,pstarl);

  TString inputfile = "fitSamples/pdfSample"; inputfile += sam; inputfile += ".root";
  cout << "File = " << inputfile << endl;	
  TChain c("ntp1");
  c.Add(inputfile);
  RooDataSet  data("data","data",myVars);

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
  int transform = 1;
  TCanvas mm("mm","KEYS fits to mmiss-pstarl",1200,800);
  gStyle->SetPalette(1);
  double All = 8;
  TH2F rotated("rotated","Rotated m^{2}_{miss}-p*_{l}",200,-All,All,200,-All,All);
  TH2F ori("ori","Original m^{2}_{miss}-p*_{l}",200,-All,All,200,-All,All);
  //TH2F totcov("ori2","Original m^{2}_{miss}-p*_{l}",200,-All,All,200,-All,All);
  double r11, r12, Xmean, Ymean;
  double x[] = {-2,-1,1,2};
  double y[] = {-4,-2,2,4};
  if(transform ==1){
    c.Draw("candPstarLep:candM2>>cov(200,-4,12,200,0,2.4)","","contz");
    TH2F *totcov = (TH2F*)gDirectory->Get("cov");
    //for(int i=0;i<4;i++)totcov.Fill(x[i],y[i]);
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
  }
  double mmp, plp;
  double entries = c.GetEntries();
  //entries = 4;
  for (int evt = 0 ; evt < entries; evt ++) {
    ran = rand.Uniform(entries);
    //c.GetEvent(ran);
    c.GetEvent(evt);
    double Mx = candM2-Xmean, Py = candPstarLep-Ymean;
    mmp = r11*(Mx)+r12*(Py);
    plp = -r12*(Mx)+r11*(Py);
    ori.Fill(Mx,Py);
    rotated.Fill(mmp,plp);
    mmiss2.setVal(candM2);
    pstarl.setVal(candPstarLep);
//     if (MCType == 0)
//       totWeight.setVal(myWM->getCombWeight(MCCombB,MCCombDs,MCDoubleSL,candTruLep));
//     else
//       totWeight.setVal(myWM->getEventWeight(candType,candDstarType,MCType,MCSubmode,MCDssmode,MCD,MCPions,
// 					    isBzero,isSP6,MCTaumode,truePPi0,trueDmass,CTL,CTV,Chi,Q2,
// 					    trueLepCharge,candM2));
    data.add(RooArgSet(mmiss2,pstarl));
  }
  //data.setWeightVar(totWeight);
  ori.Draw("contz");
  mm.SaveAs("original.eps");
  rotated.Draw("contz");
  mm.SaveAs("rotated.eps");
  cout<<"("<<r11<<", "<<r12<<") and covariance "<<rotated.GetCovariance()<<endl;
  return;

  double smoo = smooth.Atof()/100.;
  Roo2DKeysPdf DPpdf("DPpdf","DPpdf",mmiss2,pstarl,data,"av",smoo);
  time (&end);dif = difftime (end,start);
  cout<<dif<<" seconds after finding the KEYS function"<<endl;
  time (&start);
  int ntotbin = 800;
  TH2F *h2 = new TH2F("h2","KEYS",ntotbin,-4,12,ntotbin,0,2.4);
  DPpdf.fillHistogram(h2,RooArgList(mmiss2,pstarl));
  TString hname = "AWG82/results/keys/root/hKeys"; hname += sam; hname += "_"; hname += smooth; hname += ".root";
  TFile* hfile = new TFile(hname,"RECREATE"); 
  h2->Write();
  hfile->Close();
  cout<<"KEYS histogram saved in "<<hname<<endl;
  RooDataHist* Rdh2 = new RooDataHist("Rdh2","KEYS",RooArgList(mmiss2,pstarl),h2);
  RooHistPdf* Rh2 = new RooHistPdf("Rh2","KEYS",RooArgList(mmiss2,pstarl),*Rdh2,2);
  time (&end);dif = difftime (end,start);
  cout<<dif<<" seconds after making histogram"<<endl;
  time (&start);

  Float_t xlow,xhigh;
  Int_t nbinx,nbiny,Sam = sam.Atoi();
  xlow = -4;
  xhigh = 12;
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
  TString Pltitles[] = {"-4 < m^{2}_{miss} < 1.5 GeV^{2}","1.5 < m^{2}_{miss} < 12 GeV^{2}",
			"-4 < m^{2}_{miss} < 12 GeV^{2}"};
  TString M2cuts[] = {"candPstarLep<1","candPstarLep>1&&candPstarLep<1.4",
		      "candPstarLep>1.4&&candPstarLep<1.8","candPstarLep>1.8&&candPstarLep<2.4", ""};
  TString Plcuts[] = {"candM2<1.5","candM2>=1.5",""};
  double limits[] = {0, 1, 1.4, 1.8, 2.4};
  int binlim[5];
  for(int i=0;i<5;i++) binlim[i] = limits[i]/2.4*ntotbin;

  TString psname = "AWG82/results/keys/eps/eps2Keys"; psname+=sam; psname+="_";
  psname += smooth; psname += ".ps";
  double tot = 0;
  mm.Print(psname+"[");
  TH1F *hm2[5], *m2[5], *hpl[3], *pl[3];
  TString M2names[5], Plnames[3];
  for(int i=0;i<5;i++){
    M2names[i] = "hm2_"; M2names[i] += i;
    hm2[i] = new TH1F(M2names[i],M2titles[i],ntotbin,-4,12); 
    if(i<3) {
      Plnames[i] = "hpl_"; Plnames[i] += i;
      hpl[i] = new TH1F(Plnames[i],Pltitles[i],ntotbin,0,2.4); 
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
      for(int j=1; j<ntotbin+1; j++){
	double binVal = 0;
	for(int binp = binlim[i]+1; binp < binlim[i+1]+1; binp++){
	  binVal += h2->GetBinContent(j,binp)*entries*ntotbin*(xhigh-xlow)/nbinx/16;	
	}
	hm2[i]->SetBinContent(j,binVal);
      }
    } 
    hm2[i]->SetLineColor(4);
    hm2[i]->SetLineWidth(2);
    if(i<4) hm2[4]->Add(hm2[i]);
  }
  int plbinlim[3] = {0,ntotbin*5.5/16,ntotbin};
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
      for(int j=1; j<ntotbin+1; j++){
	double binVal = 0;
	for(int binp = plbinlim[i]+1; binp < plbinlim[i+1]+1; binp++){
	  binVal += h2->GetBinContent(binp,j)*entries*ntotbin/nbiny;	
	}
	hpl[i]->SetBinContent(j,binVal);
      }
    }
    hpl[i]->SetLineColor(4);
    hpl[i]->SetLineWidth(2);
    if(i<2) hpl[2]->Add(hpl[i]);
  }
  m2[4]->Draw("e1"); hm2[4]->Draw("c same"); mm.Print(psname);
  pl[2]->Draw("e1"); hpl[2]->Draw("c same"); mm.Print(psname);
  for(int i=0;i<4;i++){
    m2[i]->Draw("e1"); hm2[i]->Draw("c same"); mm.Print(psname);
  }
  for(int i=0;i<2;i++){
    pl[i]->Draw("e1"); hpl[i]->Draw("c same"); mm.Print(psname);
  }
  mm.Print(psname+"]");
  time (&end);dif = difftime (end,start);
  cout<<dif<<" seconds after plotting data. Written "<<psname<<endl;
  return; 

} 
