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
  TRandom3 rand; int ran, entries = c.GetEntries();
  for (int evt = 0 ; evt < entries; evt ++) {
    ran = rand.Uniform(entries);
    c.GetEvent(ran);
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


  double smoo = smooth.Atof()/100.;
  Roo2DKeysPdf DPpdf("DPpdf","DPpdf",mmiss2,pstarl,data,"av",smoo);
  time (&end);dif = difftime (end,start);
  cout<<dif<<" seconds after finding the KEYS function"<<endl;
  return;
  time (&start);
  TH2F *h2 = new TH2F("h2","KEYS",100,-4,12,100,0,2.4);
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
  RooDataSet *proj=(RooDataSet*)Rh2->generate(RooArgSet(mmiss2,pstarl),2000);
  time (&end);dif = difftime (end,start);
  cout<<dif<<" seconds after generating toy"<<endl;
  time (&start);

  Float_t xlow,xhigh;
  Int_t nbinx,nbiny,Sam = sam.Atoi();
  xlow = -4;
  xhigh = 12;
  nbinx = 80;
  
  if (Sam==0 || Sam==2 || Sam==10 || Sam==12 || Sam == 20 || Sam == 23 || Sam == 26 || Sam == 29) {
    xlow = -2; xhigh = 4;
    if (Sam > 12) nbinx = 40; 
  } else if (Sam==1 || Sam==11) {
    xlow = -2; xhigh = 6;
  }
  if (Sam==6 || Sam==7 || Sam==16 || Sam==17)  nbinx = 40; 
  if (Sam==8 || Sam==18) {
    xhigh = 4; nbinx = 40; 
  }
  if (Sam==9 || Sam==19)  nbinx = 40; 
  if (Sam==21 || Sam==22 || Sam==24 || Sam==25 || Sam==27 || Sam==28 || Sam==30 || Sam==31) {
    xhigh = 8; nbinx = 40; 
  }
  if (Sam > 31)  nbinx = 40; 
  TString hnames[] = {"0 < p*_{l} < 1","1 < p*_{l} < 1.4","1.4 < p*_{l} < 1.8","1.8 < p*_{l} < 2.4"};
  TString cuts[] = {"candPstarLep<1","candPstarLep>1&&candPstarLep<1.4",
		    "candPstarLep>1.4&&candPstarLep<1.8",
		    "candPstarLep>1.8&&candPstarLep<2.4"};
  RooDataSet *toys[4], *datas[4];
  RooPlot *mmslice[4];
  for(int i=0;i<4;i++){
    toys[i] = (RooDataSet*)proj->reduce(cuts[i]);
    datas[i] = (RooDataSet*)data.reduce(cuts[i]);
    mmslice[i] = mmiss2.frame(xlow,xhigh,nbinx);
    mmslice[i]->SetTitle(hnames[i]);
  }
  time (&end);dif = difftime (end,start);
  cout<<dif<<" seconds after reducing"<<endl;
  time (&start);

  TCanvas mm("mm","Missing mass");
  //mm.Divide(2,2);
  for(int i=0;i<4;i++){
    //mm.cd(i+1);
    datas[i]->plotOn(mmslice[i],DataError(RooAbsData::SumW2));
    Rh2->plotOn(mmslice[i],ProjWData(pstarl,*toys[i]));
    time (&end);dif = difftime (end,start);
    mmslice[i]->Draw();
    mmslice[i]->SetXTitle("m^{2}_{miss} [GeV^{2}]");
    TString epsname = "AWG82/results/keys/eps/epsKeys"; epsname+=sam; epsname+="_";
    epsname += smooth; epsname += "_"; epsname += i+1; epsname += ".eps";
    mm.SaveAs(epsname);
    cout<<dif<<" seconds after plotting data "<<i<<endl;
    time (&start);
  }
  return; 

} 
