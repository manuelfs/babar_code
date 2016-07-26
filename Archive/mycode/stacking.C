#include "TLegend.h"

void stacking(TString Svari="candQ2",float hmin=-4,float hmax=14,TString title="",TString xtitle="",
	      TString fname="Q2",float cut=4,TString legwhere="left",bool doLog=false){
  TString isAdd = "Add";
  // Loading files
  TChain ch("ntp1");
  if(isAdd==""){
    ch.Add("AWG82/ntuples/SP-1235-BSemiExcl-Run2-R22d-v06/*");
    ch.Add("AWG82/ntuples/SP-1235-BSemiExcl-Run4-R22d-v06/*");
    ch.Add("AWG82/ntuples/SP-1237-BSemiExcl-Run2-R22d-v06/*");
    ch.Add("AWG82/ntuples/SP-1237-BSemiExcl-Run4-R22d-v06/*");
  } else {
    ch.Add("AWG82/ntuples/SP-1235-BSemiExclAdd-Run2-R22f-v01/*");
    ch.Add("AWG82/ntuples/SP-1235-BSemiExclAdd-Run4-R22f-v01/*");
    ch.Add("AWG82/ntuples/SP-1237-BSemiExclAdd-Run2-R22f-v01/*");
    ch.Add("AWG82/ntuples/SP-1237-BSemiExclAdd-Run4-R22f-v01/*");
  }

  ch.SetMakeClass(1);
  Float_t vari;
  Int_t MCType;
  Int_t candType;
  Float_t         candPMiss;
  Float_t         candQ2;
  Float_t         candEExtra;
  TBranch *b_vari;   
  TBranch *b_MCType;  
  TBranch *b_candType;  
  TBranch        *b_candQ2;   //!
  TBranch        *b_candPMiss;   //!
  TBranch        *b_candEExtra;   //!
  ch.SetBranchAddress(Svari, &vari, &b_vari);
  ch.SetBranchAddress("MCType", &MCType, &b_MCType);
  ch.SetBranchAddress("candType", &candType, &b_candType);
  ch.SetBranchAddress("candEExtra", &candEExtra, &b_candEExtra);
  ch.SetBranchAddress("candQ2", &candQ2, &b_candQ2);
  ch.SetBranchAddress("candPMiss", &candPMiss, &b_candPMiss);

  TCut MMiss = "candPMiss>.2";
  TCut Q2 = "candQ2>4";
  TCut ee1 = "candType<3&&candEExtra<.2";
  TCut ee2 = "candType==3&&candEExtra<.15";
  TCut ee3 = "candType==4&&candEExtra<.3";
  TCut ee = (ee1||ee2||ee3);
  TCut basic= MMiss+Q2+ee;

  TH1F *h[10];
  for(int i=0; i<10; i++){
    TString hname = "h"; hname += i;
    h[i] = new TH1F(hname,"",60,hmin,hmax);
  }
  THStack hs("hs",title);

  int hnum = 9;
  Long64_t nentries = ch.GetEntries();
  //   nentries = 1000;
  Long64_t ntotal = 0;
  for(Long64_t entry=0; entry<nentries; entry++){
    if(!(entry%5000)) 
      cout<<"Done "<<entry<<" of "<<nentries<<" entries"<<endl;
    ch.LoadTree(entry);
    ch.GetEntry(entry);
    if(!(candPMiss>.2&&candQ2>4&&(candType<3&&candEExtra<.2||candType==3&&candEExtra<.15||candType==4&&candEExtra<.3))) continue;
    if(!(candType==1||candType==3)) continue;
    switch(MCType){
    case 0:
      hnum = 9;
      break;
    case 1:
    case 3:
      hnum = 4;
      break;
    case 2:
    case 4:
      hnum = 5;
      break;
    case 5:
      hnum = 0;
      break;
    case 6:
      hnum = 1;
      break;
    case 7:
    case 9:
      hnum = 6;
      break;
    case 8:
    case 10:
      hnum = 7;
      break;
    case 11:
      hnum = 2;
      break;
    case 12:
      hnum = 3;
      break;
    case 13:
    case 14:
      hnum = 8;
      break;
    }
    h[hnum]->Fill(vari); ntotal++;
  }

  TLegend *leg;
  if(legwhere=="left")
    leg= new TLegend(0.1,0.45,0.22,0.9);
  else
    leg= new TLegend(0.78,0.45,0.9,0.9);

  leg->SetFillColor(0);
  int colors[] = {2,3,4,5,6,7,8,9,10,11}; //12,13,17,15,16
  TString lnames[] = {"D^{0} #tau #nu","D*^{0} #tau #nu","D^{+} #tau #nu",
		      "D*^{+} #tau #nu","D^{0} l #nu","D*^{0} l #nu",
		      "D^{+} l #nu","D*^{+} l #nu","D** l #nu","Bkg"};
  for(int i=0; i<10; i++){
    h[i]->SetFillColor(colors[i]);
    hs.Add(h[i]);
    leg->AddEntry(h[i],lnames[i]);
  }
  TCanvas c("Yields","Yields",600,400);
  hs.Draw();
  leg->Draw();
  hs.GetXaxis()->SetTitle(xtitle);
  float hMaxi = 1.05*hs.GetMaximum();
  if(doLog){
    gPad->SetLogy(1);
    fname += "Log";
    hMaxi = pow(hs.GetMaximum(),1.05);
  }
  TLine *line= new TLine(0,0,0,0);
  line->SetLineColor(4);
  line->SetLineStyle(2);
  line->SetLineWidth(3.5);
  //  line->DrawLine(cut,0,cut,hMaxi);
  leg->Draw();
  TString stotal = "Entries: "; stotal += ntotal;
  TLatex *label = new TLatex();
  label->SetNDC(kTRUE);
  label->SetTextSize(0.04);
  label->DrawLatex(.72,0.93,stotal);
  
  if(isAdd==""){
    c.SaveAs("mycode/eps/NNVeryLoose/"+fname+".eps");
  } else {
    c.SaveAs("mycode/eps/Add/"+fname+".eps");
  }

  for(int i=0; i<10; i++)
    h[i]->Delete();
}


