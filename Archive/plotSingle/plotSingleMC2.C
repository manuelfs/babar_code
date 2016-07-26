//---------------------------------------------------------------------------------
// Description:
//      Plots the a single variable for the 4 channels
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      09/06/30 manuelf -- Upgrade
//      09/05/31 manuelf -- Creation
//---------------------------------------------------------------------------------

#include "mycode/cuts.hh"
TString round(double n, int e, double d=1.);

void plotSingleMC2(TString Tag="Dl"){

  //TString mfile = "AWG82/ntuples/small/Add_Truth30Rest_Run123456.root";
  TString mfile = "AWG82/ntuples/small/Add_Aug20Rest_Run123456.root";
  TChain gen("ntp1");
  gen.Add(mfile);
  TChain sig("ntp1");
  sig.Add("AWG82/ntuples/small/*Full*All.root");
  sig.SetLineColor(2);
  TChain old("ntp1");
  old.Add("AWG82/ntuples/small/BDT_KM_RunAll.root");
  old.SetLineColor(4);


  TCut modeCut[2];
  TString decay[2];
  TString limits[2];
  if(Tag=="Dl"){
    modeCut[0] = "(MCType==1||MCType==3)";
    modeCut[1] = "(MCType==7||MCType==9)";
    decay[0] = "D^{0} (e/#mu) #nu";
    decay[1] = "D^{+} (e/#mu) #nu";
    limits[0] = "(50,-0.3,0.5)";
    limits[1] = "(50,-1.6,0.7)";
  }else if(Tag=="Dsl"){
    modeCut[0] = "(MCType==2||MCType==4)";
    modeCut[1] = "(MCType==8||MCType==10)";
    decay[0] = "D*^{0} (e/#mu) #nu";
    decay[1] = "D*^{+} (e/#mu) #nu";
    limits[0] = "(50,-0.4,2)";
    limits[1] = "(50,-0.4,0.7)";
  }else if(Tag=="Dtau"){
    modeCut[0] = "(MCType==5)";
    modeCut[1] = "(MCType==11)";
    decay[0] = "D^{0} #tau #nu";
    decay[1] = "D^{+} #tau #nu";
    limits[0] = "(50,-0.5,9)";
    limits[1] = "(50,-1,8)";
  }else if(Tag=="Dstau"){
    modeCut[0] = "(MCType==6)";
    modeCut[1] = "(MCType==12)";
    decay[0] = "D*^{0} #tau #nu";
    decay[1] = "D*^{+} #tau #nu";
    limits[0] = "(50,-0.5,9)";
    limits[1] = "(50,-0.5,9)";
  }else if(Tag=="Dssl"){
    modeCut[0] = "(MCType>12&&MCTaumode==-1)";
    modeCut[1] = "(MCType>12&&MCTaumode==-1)";
    decay[0] = "D** (e/#mu) #nu";
    decay[1] = "D** (e/#mu) #nu";
    limits[0] = "(50,-0.5,9)";
    limits[1] = "(50,-0.5,9)";
  }else if(Tag=="Dsstau"){
    modeCut[0] = "(MCType>12&&MCTaumode>-1)";
    modeCut[1] = "(MCType>12&&MCTaumode>-1)";
    decay[0] = "D** #tau #nu";
    decay[1] = "D** #tau #nu";
    limits[0] = "(50,-0.5,9)";
    limits[1] = "(50,-0.5,9)";
  }else if(Tag=="Dsstaul"){
    modeCut[0] = "(MCType>12)";
    modeCut[1] = "(MCType>12)";
    decay[0] = "D** (e/#mu/#tau) #nu";
    decay[1] = "D** (e/#mu/#tau) #nu";
    limits[0] = "(50,-0.5,9)";
    limits[1] = "(50,-0.5,9)";
  }else if(Tag=="pi0Dssl"){
    modeCut[0] = "(MCType>12&&MCTaumode==-1)";
    modeCut[1] = "(MCType>12&&MCTaumode==-1)";
    decay[0] = "D** (e/#mu) #nu";
    decay[1] = "D** (e/#mu) #nu";
    limits[0] = "(50,-0.5,2)";
    limits[1] = "(50,-0.7,2)";
  }else if(Tag=="pi0Dsstau"){
    modeCut[0] = "(MCType>12&&MCTaumode>-1)";
    modeCut[1] = "(MCType>12&&MCTaumode>-1)";
    decay[0] = "D** #tau #nu";
    decay[1] = "D** #tau #nu";
    limits[0] = "(50,-0.5,7)";
    limits[1] = "(50,-0.5,7)";
  }else if(Tag=="pi0Dsstaul"){
    modeCut[0] = "(MCType>12)";
    modeCut[1] = "(MCType>12)";
    decay[0] = "D** (e/#mu/#tau) #nu";
    decay[1] = "D** (e/#mu/#tau) #nu";
    limits[0] = "(50,-0.5,7)";
    limits[1] = "(50,-0.7,7)";
  }else if(Tag=="pi0D0"){
    modeCut[0] = "(MCType==1||MCType==3)";
    decay[0] = "D^{0} (e/#mu) #nu";
    limits[0] = "(50,-2.5,1)";
    limits[1] = "(50,-2.7,1)";
  }else if(Tag=="pi0Ds0"){
    modeCut[0] = "(MCType==2||MCType==4)";
    decay[0] = "D*^{0} (e/#mu) #nu";
    limits[0] = "(50,-2.5,1)";
    limits[1] = "(50,-2.7,1)";
  }else {
    cout<<"Select one: D0l, Ds0l, Dpl, Dspl, D0tau, Ds0tau, Dptau, Dsptau, Dss, pi0Dss, pi0D0, pi0Ds0"<<endl;
    return;
  }
  TString Dnames[] = {"D^{0}","D*^{0}","D^{+}","D*^{+}"};
  TH1F *mmiss[4][3];
  TLatex *label = new TLatex();
  label->SetNDC(kTRUE);
  label->SetTextSize(0.08);
  TCanvas c("cMvaDl","Plot for the different channels",1000,600);
  c.Divide(2,2);
  gStyle->SetOptStat(0);
  int a=0;
  TCut Btag = "(candBMode>11000&&candBMode<11400||candBMode>12000&&candBMode<12500)";
  for(int i=0; i<4; i++){
    a = i%2; 
    c.cd(i+1);
    TString cuts = "candType=="; cuts += i+1; 
    TCut tot = modeCut[i>1]+MvaAll; tot += cuts;  
    TString vari = "candM2>>h"; vari+=limits[a]; 
    if(Tag.Contains("pi0")){
      tot = modeCut[i>1]+dssMvaAll; tot += cuts;
      vari = "mm2pi0>>h"; vari+=limits[a]; 
    }
    //tot += Btag; 
    TCut tot2 = tot; //tot2+="MCTaumode==-1"; tot2*="wBF";
    tot*="wBF";
    gen.Draw(vari,tot);
    TH1F *h = (TH1F*)gDirectory->Get("h");
    mmiss[i][0] = (TH1F*)h->Clone("cand"+i);
    sig.Draw(vari,tot);
    TH1F *h = (TH1F*)gDirectory->Get("h");
    mmiss[i][1] = (TH1F*)h->Clone("candMva"+i);
    //tot *= "wFF"; tot *= "wBF";
    old.Draw(vari,tot);
    TH1F *h = (TH1F*)gDirectory->Get("h");
    mmiss[i][2] = (TH1F*)h->Clone("candMva"+i);
    double ngen = mmiss[i][0]->Integral();
    double nsig = mmiss[i][1]->Integral();
    double nold = mmiss[i][2]->Integral();
    if(nsig) {
      if(ngen) mmiss[i][1]->Scale(ngen/nsig);
      else if(nold) mmiss[i][2]->Scale(nsig/nold);
    }      
    if(nold && ngen) mmiss[i][2]->Scale(ngen/nold);

    double maxi = mmiss[i][0]->GetMaximum();
    if(maxi<mmiss[i][1]->GetMaximum()) maxi = mmiss[i][1]->GetMaximum();
    if(maxi<mmiss[i][2]->GetMaximum()) maxi = mmiss[i][2]->GetMaximum();

    mmiss[i][0]->SetMaximum(1.15*maxi);
    mmiss[i][0]->GetXaxis()->SetTitleOffset(.7);
    mmiss[i][0]->GetXaxis()->SetTitleSize(0.059);
    mmiss[i][0]->GetYaxis()->SetLabelSize(0.062);
    mmiss[i][0]->GetXaxis()->SetLabelSize(0.056);
    mmiss[i][0]->GetYaxis()->SetNdivisions(5+100*2);
    mmiss[i][0]->GetXaxis()->SetNdivisions(6+100*2);
    mmiss[i][0]->SetLineWidth(2);
    mmiss[i][0]->SetTitle("");
    mmiss[i][0]->SetXTitle("m^{2}_{miss} [GeV^{2}]");
    mmiss[i][0]->Draw();
    mmiss[i][1]->SetLineWidth(2);
    mmiss[i][1]->SetLineColor(2);
    mmiss[i][2]->Draw("same");
    mmiss[i][0]->Draw("same");
    mmiss[i][1]->Draw("same");
    TString tag = Dnames[i]; if(Tag.Contains("pi0")) tag += "#pi^{0}";
    tag+= " channel: "; tag+= decay[i>1]; 
    label->DrawLatex(.11,0.92,tag);
    TLegend *leg;
    if(Tag=="pi0Ds0"||Tag=="pi0D0")
      leg = new TLegend(0.1,.66,0.43,0.9);
    else
      leg = new TLegend(0.57,.66,0.9,0.9);
    leg->SetTextSize(0.058);
    leg->SetFillColor(0);
    TString legtag = "BrecoAdd ("; legtag+=round(ngen,0); legtag+=")";
    leg->AddEntry(mmiss[i][0],legtag);
    legtag = "Breco ("; legtag+=round(nold,0); legtag+=")";
    leg->AddEntry(mmiss[i][2],legtag);
    legtag = "Cocktail ("; legtag+=round(nsig,0); legtag+=")";
    leg->AddEntry(mmiss[i][1],legtag);
    leg->Draw();
  }
  c.SaveAs("babar_code/eps/"+Tag+"MC2compMmissBF.eps");
}


TString round(double n, int e, double d){
  if(d==0) return " - ";
  double b = (int)(n/d*pow(10,e)+0.5);
  b /= pow(10,e);
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(b<10 && e==2){
    if(result.Length() == 1)
      result += ".00";
    if(result.Length() == 3)
      result += "0";
  }
  if(e==1 && !result.Contains(".")) result += ".0";
  return result;
}
