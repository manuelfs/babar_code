void NDKeys(TString sam="18",int nM2bin = 100, int nPlbin = 100){
  gSystem->Load("libHtml");
  gSystem->Load("libMinuit");
  gSystem->Load("libRooFitCore.so");
  gSystem->Load("libRooFitModels.so");
  using namespace RooFit;

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

  double m2min = -4, m2max = 4, plmin = 0, plmax = 2.4;
  RooRealVar mmiss2("candM2","candM2",m2min,m2max);
  RooRealVar pstarl("candPstarLep","candPstarLep",plmin,plmax);
  RooArgSet myVars(mmiss2,pstarl);
  RooDataSet  data("data","data",myVars);
  double entries = c.GetEntries();
  //entries = 10;
  for (int evt = 0 ; evt < entries; evt ++) {
    c.GetEvent(evt);
    mmiss2.setVal(candM2);
    pstarl.setVal(candPstarLep);
      data.add(RooArgSet(mmiss2,pstarl));
  }


  RooNDKeysPdf keys("DPpdf","DPpdf",RooArgList(mmiss2,pstarl),data,"a",1,3);
  TH2F *h2= new TH2F("h2","KEYS",nM2bin,m2min,m2max,nPlbin,plmin,plmax);
  keys.fillHistogram(h2,RooArgList(mmiss2,pstarl));
  TCanvas mm("mm","KEYS fits to mmiss-pstarl",1200,800);
  h2->Draw("lego");
  mm.SaveAs("testKeys.eps");
}


