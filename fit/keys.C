#include "Roo2DKeysPdf.h"

void keys() {

  RooRealVar mmiss2("candM2","candM2",-4,12);
  RooRealVar pstarl("candPstarLep","candPstarLep",0.,2.4);

  RooArgSet myVars(mmiss2,pstarl);

  TString inputfile = "fitSamples/pdfSample3.root";
  cout << "File = " << inputfile << endl;	
  TChain c("ntp1");
  c.Add(inputfile);
  RooDataSet  data("data","data",&c,myVars);

  Roo2DKeysPdf DPpdf("DPpdf","DPpdf",mmiss2,pstarl,data,"mav",0.4);

  RooPlot* mppframe=mpp.frame(50);
  data->plotOn(mppframe);
  DPpdf.plotOn(mppframe);
  mppframe->Draw();
  return; 

} 
