void mESCorr(){
  TH1F h("h","",4,25,225);
  h.SetBinContent(1,0.88727); h.SetBinError(1,0.06417); 
  h.SetBinContent(2,0.92976); h.SetBinError(2,0.05066); 
  h.SetBinContent(3,0.93569); h.SetBinError(3,0.03959); 
  h.SetBinContent(4,0.93735); h.SetBinError(4,0.03262); 
  TF1* line = new TF1("line","[0]+[1]*x",-25,225);
  line->SetParameter(0,0.85); line->SetParameter(1,0.0002); 
  h.Fit("line");


  TH1F h("h","",6,-75,225);
  //h.SetBinContent(1,1.34106); h.SetBinError(1,0.44771); 
  h.SetBinContent(2,1.00911); h.SetBinError(2,0.13262); 
  h.SetBinContent(3,0.88727); h.SetBinError(3,0.06417); 
  h.SetBinContent(4,0.92976); h.SetBinError(4,0.05066); 
  h.SetBinContent(5,0.93569); h.SetBinError(5,0.03959); 
  h.SetBinContent(6,0.93735); h.SetBinError(6,0.03262); 
  TF1* line = new TF1("line","[0]+[1]*x",-100,225);
  line->SetParameter(0,0.85); line->SetParameter(1,0.0002); 
  h.Fit("line");



}
