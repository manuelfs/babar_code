void test(){

//   TFile f1("~/KKKs/KKKs-RooDalitz-a51/workdir/mKKSplotspp.root");
//   TFile f2("~/KKKs/KKKs-RooDalitz-a51/workdir/mKmKsSplotspp.root");
//   TFile f3("~/KKKs/KKKs-RooDalitz-a51/workdir/mKpKsSplotspp.root");
  //  TFile f4("~/KKKs/KKKs-RooDalitz-a51/workdir/mKKSplotspp.root");

 //  TFile f1("~/KKK/KKK-RooDalitz-a51-test/workdir/mKKlowSplots.root");
//   TFile f2("~/KKK/KKK-RooDalitz-a51-test/workdir/mKKhighSplots.root");

//   TFile f3("~/KKsKs/KKsKs-RooDalitz-a51/workdir/mvarSplots.root");
//   TFile f4("~/KKsKs/KKsKs-RooDalitz-a51/workdir/thvarSplots.root");

//   TCanvas* canv = new TCanvas("canv","canv");
//   canv->Print("c1.ps[");



//   RooPlot* p1 = (RooPlot*)f1.Get("totalBkgmKK");
//   p1->Draw();
//   canv->Print("c1.ps");

//   RooPlot* p2 = (RooPlot*)f2.Get("totalBkgmKKhigh");
//   p2->Draw();
//   canv->Print("c1.ps");



//   RooPlot* p3 = (RooPlot*)f3.Get("totalBkgmKK");
//   p3->Draw();
//   canv->Print("c1.ps");


//   RooPlot* p4 = (RooPlot*)f4.Get("totalBkgthvar");
//   p4->Draw();
//   canv->Print("c1.ps");








// //   RooPlot* p3 = (RooPlot*)f3.Get("totalBkgmKK");
// //   p3->Draw();
// //   canv->Print("c1.ps");

//   canv->Print("c1.ps]");
 
  


//   TFile f1("~/KKKs/KKKs-RooDalitz-a51/workdir/mesSplotspp.root");
//   TFile f2("~/KKKs/KKKs-RooDalitz-a51/workdir/dESplotspp.root");
//   TFile f3("~/KKKs/KKKs-RooDalitz-a51/workdir/mKKSplotspp.root");
//   TFile f4("~/KKKs/KKKs-RooDalitz-a51/workdir/mKmKsSplotspp.root");
//   TFile f5("~/KKKs/KKKs-RooDalitz-a51/workdir/mKpKsSplotspp.root");


//   TCanvas* canv = new TCanvas("canv","canv");
//   canv->Print("c1.ps[");

//   RooPlot* p1 = (RooPlot*)f1.Get("signalppMes");
//   p1->Draw();
//   canv->Print("c1.ps");

//   RooPlot* p2 = (RooPlot*)f2.Get("signalppdE");
//   p2->Draw();
//   canv->Print("c1.ps");


//   RooPlot* p3 = (RooPlot*)f3.Get("signalppmKK");
//   p3->Draw();
//   canv->Print("c1.ps");


//   RooPlot* p4 = (RooPlot*)f4.Get("signalppmKK");
//   p4->Draw();
//   canv->Print("c1.ps");


//   RooPlot* p5 = (RooPlot*)f5.Get("signalppmKK");
//   p5->Draw();
//   canv->Print("c1.ps");




// //   RooPlot* p3 = (RooPlot*)f3.Get("totalBkgmKK");
// //   p3->Draw();
// //   canv->Print("c1.ps");

//   canv->Print("c1.ps]");


  RooRealVar* mes = new RooRealVar("mES","mES",5.28,5.26,5.29);
  RooRealVar* de = new RooRealVar("dE","dE",0.0,-0.2,0.2);

  RooRealVar* mes0 = new RooRealVar("mes0","mes0",5.282,5.26,5.29);
  RooRealVar* de0 = new RooRealVar("de0","de0",0.03,-0.2,0.2);

  // RooRealVar* rho = new RooRealVar("rho","rho",5,0,10);
  

  RooArgList xvec(*mes,*de);

  RooArgList muvec(*mes0,*de0);

  TMatrixDSym* corr = new TMatrixDSym(2);

  double rho = 0.9;
  
  (*corr)(0,0) = 1.0;
  (*corr)(0,1) = rho;
  (*corr)(1,0) = rho;  
  (*corr)(1,1) = 1.0;
 
  




  RooMultiVarGaussian* a = new RooMultiVarGaussian("pdf","pdf",xvec,muvec,*corr);
  
  //  a->Print("V");

  cout<<"generating 100 events..."<<endl;

  RooDataSet* dset = a->generate(xvec,1000);

  dset->tree()->Draw("mES:dE");


}
