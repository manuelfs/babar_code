#include "TString.h"
#include "TPad.h"
#include "TCut.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TH1F.h"
#include "THStack.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TChain.h"
#include <fstream>
#include <iostream>

using namespace TMath;
using namespace std;
using std::cout;
using std::endl;

void formatHisto(TH1F *h);
TString RoundNumber(double n, int e, double d=1);

void CSampleChiara(TString Variable="MissMass", int nbins=25, float minX=-0.4, float maxX=0.7, 
		   double q2Cut = 4, double EexCut = 1, bool doScale=true, TString sample = "Chi",
		   TString Channel = "D0"){
  int chan = 0;
  if(Channel == "D0") chan = 1;
  if(Channel == "Dp") chan = 2;

  TString treeName = "MyTuple"; treeName += Channel; treeName += "Excl";
  TString folder = "AWG82/ntuples/dlopes/Eex_125_"; folder += Channel;
  TChain BB(treeName);
  TChain Data(treeName);
  TString BBName = folder; BBName += "/*BB*";
  TString DataName = folder; DataName += "/*Data*";
  BB.Add(BBName);
  Data.Add(DataName);

  TChain BBMan("ntp1");
  TChain DataMan("ntp1");
  TChain uds("ntp1");
  TChain ccbar("ntp1");
  BBMan.Add("AWG82/ntuples/small/Add_R24Rests_RunAll.root");
  DataMan.Add("AWG82/ntuples/small/Adddata_RunAll.root");
  uds.Add("AWG82/ntuples/small/uds_RunAll.root");
  ccbar.Add("AWG82/ntuples/small/ccbar_RunAll.root");
  double NBRun[] = {36968000+37200000., 103498000+103124000., 50556000+49766000., 
		   167994000+167466000., 244322000+244812000., 68148000+68016000.};
  double ndata[] = {22389980.4, 67394307.5, 35569248.8, 110449802.7, 
		    147190396.5, 84767412.6};
  double NBRunMan[2][6] = {{34878000, 101690000, 56035000, 166784000, 215168000, 130336000},
			{34941000, 104188000, 57888000, 169801000, 215953000, 135224000}};
  double nccbar[] = {55254000, 164146000, 88321000, 267308000, 343667000, 208664000};
  double nuds[] = {160514000, 451636000, 275869000, 421599000, 553604000, 327032000};
  double ndataMan[] = {22556256.9, 68438426.0, 35763257.9, 111429669.4, 147620363.4, 85194672.2};
  double totMCB = 0, totdata = 0, totuds = 0, totccbar = 0, Bentries = 0, Dentries = 0, Byield[6];
  for(int i=0; i<6; i++){
    if(sample=="Chi"){
      totdata += ndata[i];	                      //  467.761M events
      totMCB += NBRun[i];                             // 1341.870M events
    } else {
      totdata += ndataMan[i];			      //  471.003M events
      totMCB += (NBRunMan[0][i]+NBRunMan[1][i])*2/3.; //  948.591M events
    }
    totuds += nuds[i];				      // 2190.254M events
    totccbar += nccbar[i];			      // 1127.360M events
  }
  double wuds = totMCB/totuds*2.09/1.05;              // 0.862068
  double wccbar = totMCB/totccbar*1.3/1.05;           // 1.041766
  double wMC = totdata/totMCB;                        // 0.348589 Chi; 0.496529 Man

  TString Titles[] = {"D^{0}: MC/data = ","D^{+}: MC/data = "};
  TString VariableM = "candM2NF";
  if(Variable == "Eextra") VariableM = "candEExtra";
  if(Variable == "plepcms") VariableM = "candPstarLep";
  double mB = 5.279, mDp = 1.869;
  double wCut = (mB*mB+mDp*mDp-q2Cut)/(2*mB*mDp);
  TString wCut_s = "w2>"; wCut_s += wCut; wCut_s += "&&Eextra<"; wCut_s += EexCut;
  TCut baseCut = "mESB>5.27"; baseCut += wCut_s;
  TString samCut[2][4] = {{"TrueDs0==1","TrueD0==1","TrueDss==1",
			  "(TrueDsp==1||TrueDp==1||Cascade==1||fake==1)"},
			 {"TrueDsp==1","TrueDp==1","TrueDss==1",
			  "(TrueDs0==1||TrueD0==1||Cascade==1||fake==1)"}};
  TString wCutM_s = "candQ2<"; wCutM_s += q2Cut; wCutM_s += "&&candEExtra<"; wCutM_s += EexCut;
  wCutM_s += "&&candType=="; wCutM_s += (chan + (chan==2));
  TCut baseCutM = "candPMiss>0.2"; baseCutM += wCutM_s;
  TString samCutM[2][4] = {{"(MCType==2||MCType==4||MCType==6)","(MCType==1||MCType==3||MCType==5)",
			    "MCType>12","(MCType>6&&MCType<13||MCType==0)"},
			   {"(MCType==8||MCType==10||MCType==12)","(MCType==7||MCType==9||MCType==11)",
			    "MCType>12","MCType<7"}};

  TString xtitle = Variable;
  if(Variable=="MissMass") xtitle = "m^{2}_{miss} [GeV^{2}]"; 
  if(Variable=="plepcms") xtitle = "p*_{l} [GeV]"; 
  if(Variable=="Eextra") xtitle = "E_{Extra} [GeV]"; 
  TString PosLeg = "r";
  bool AddConti = true;
  TString legTag[] = {"D*l#nu (", "Dl#nu (", "D**l#nu (", "Bkg (", "Cont ("};
  TLatex *label = new TLatex(); label->SetNDC(kTRUE);
  int colors[] = {9,2,28,3,5};
  TH1F *hdata, *hdata2, *hStack[6], *hMCsum=0;
  THStack *hs; hs = new THStack("Stack","");
  gROOT->SetStyle("Plain");gStyle->SetOptStat(0);
  TCanvas c("dataMC","data Vs MC",450,350);
  TPad *p1 = new TPad("p1","",0,.25,1,1);
  p1->Draw(); p1->cd();
  hdata = new TH1F("data","",nbins,minX,maxX);
  hdata->Sumw2();
  if(sample=="Chi") Data.Project("data",Variable,baseCut);
  else DataMan.Project("data",VariableM,baseCutM);
  Dentries = hdata->Integral();
  for(int sam=0; sam<4; sam++){
    TString hname = "BB"; hname += sam;
    hStack[sam] = new TH1F(hname,"",nbins,minX,maxX);
    hStack[sam]->Sumw2();
    if(sample=="Chi") {
      TCut BBCut = baseCut; BBCut += samCut[chan-1][sam]; //BBCut *= "weightB0";
      BB.Project(hname,Variable,BBCut);
    } else {
      TCut BBCutM = baseCutM; BBCutM += samCutM[chan-1][sam]; 
      BBMan.Project(hname,VariableM,BBCutM);
    }
    hStack[sam]->SetFillColor(colors[sam]);
    hStack[sam]->SetLineColor(colors[sam]);
    Byield[sam] = hStack[sam]->Integral();
    hStack[sam]->Scale(wMC);
    Bentries += Byield[sam];
    Byield[sam] = Byield[sam]*wMC;
  }
  if(VariableM=="candM2NF") VariableM="candM2"; 
  hStack[4] = new TH1F("BB4","",nbins,minX,maxX);
  hStack[4]->Sumw2();
  uds.Project("BB4",VariableM,baseCutM);
  hStack[5] = new TH1F("BB5","",nbins,minX,maxX);
  hStack[5]->Sumw2();
  ccbar.Project("BB5",VariableM,baseCutM);
  hStack[5]->Scale(wccbar/wuds);
  hStack[4]->Add(hStack[5]);
  hStack[4]->SetFillColor(colors[4]);
  hStack[4]->SetLineColor(colors[4]);
  Byield[4] = hStack[4]->Integral()*wuds;
  hStack[4]->Scale(wMC);
  if(AddConti) Bentries += Byield[4];
  Byield[4] = Byield[4]*wMC;

  int samIni = 3;
  if(AddConti) samIni = 4;
  double MCdata = Bentries/Dentries*wMC;
  for(int sam=samIni; sam>=0; sam--){
    if(doScale) hStack[sam]->Scale(1/MCdata);
    if(sam==samIni) hMCsum = (TH1F*)hStack[sam]->Clone("MCSum");
    else hMCsum->Add(hStack[sam]);
    hs->Add(hStack[sam],"hist");
  }

  formatHisto(hdata);
  float maxi = hdata->GetMaximum();
  if(hs->GetMaximum()>maxi) maxi = hs->GetMaximum();
  hdata->SetMaximum(1.15*maxi);
  hs->Draw("");
  hs->GetYaxis()->SetLabelSize(0.062);
  hs->GetXaxis()->SetLabelSize(0.075);
  hs->GetYaxis()->SetNdivisions(6+100*2);
  hs->GetXaxis()->SetNdivisions(7+100*2);
  hdata->SetMarkerSize(1);
  hdata->Draw("a");
  hs->Draw("same");
  hdata->Draw("same");

  double errRatio = sqrt(Bentries/Dentries/Dentries+Bentries*Bentries/Dentries/Dentries/Dentries/Dentries*Dentries);
  TString RatioTitle = Titles[chan-1]; RatioTitle += RoundNumber(Bentries*wMC,2,Dentries);
  RatioTitle+=" #pm "; RatioTitle += RoundNumber(errRatio*wMC,2);
  label->SetTextSize(0.075);  label->DrawLatex(0.12,0.915,RatioTitle);
  if(PosLeg!="no"){
    TLegend *leg;
    if(PosLeg=="l") leg = new TLegend(0.1,.5,0.34,0.9);
    else leg = new TLegend(0.68,.52,0.9,0.9);
    leg->SetTextSize(0.058);
    leg->SetFillColor(0);
    TString dleg = "data ("; 
    dleg += RoundNumber(Dentries,0); dleg += ")";
    leg->AddEntry(hdata,dleg);
    for(int sam=0; sam<5; sam++){
      dleg = legTag[sam]; dleg += RoundNumber(Byield[sam],0); dleg += ")";
      leg->AddEntry(hStack[sam], dleg);
    }
    leg->Draw();
  }

  c.cd();
  TPad *p2 = new TPad("p2","",0,0.0,1,.2);
  p2->Draw(); p2->cd();
  hdata2 = (TH1F *)hdata->Clone("data2");
  hdata2->Divide(hMCsum);
  hdata2->SetMaximum(1.4);
  hdata2->SetMinimum(0.6);
  hdata2->GetYaxis()->SetNdivisions(2+100*2);
  hdata2->SetTitle("");
  hdata2->SetXTitle("");
  hdata2->GetYaxis()->SetLabelSize(0.22);
  hdata2->GetXaxis()->SetLabelSize(0.01);
  hdata2->SetMarkerSize(1);
  hdata2->Draw("E0 P");
  int dof=0; double chi2=0;
  for(int bin = 1; bin<hdata2->GetNbinsX()+1; bin++){
    double vBin = hdata2->GetBinContent(bin), eBin = hdata2->GetBinError(bin);
    if(vBin != 0 && eBin != 0){
      chi2 += pow((vBin-1.)/eBin,2);
      dof++;
    }
  }
  if(doScale) dof--;
  double ChiProb = Prob(chi2,dof);
  TLine l;
  l.SetLineColor(2);
  l.DrawLine(minX, 1.0, maxX, 1.0);
  c.cd();
  TString ChiProb_s = RoundNumber(100*ChiProb,1); ChiProb_s += "%";
  label->SetTextSize(0.048);
  label->DrawLatex(0.91,0.164,"#tilde#chi^{2} =");
  label->DrawLatex(0.908,0.115,RoundNumber(chi2,2,(double)dof));
  label->DrawLatex(0.91,0.069,"p =");
  label->DrawLatex(0.908,0.02,ChiProb_s);
  label->SetTextSize(0.06);
  label->DrawLatex(0.68,0.21,xtitle);
  TString Scaled_s = "Not scaled"; 
  if(doScale) Scaled_s = "Scaled";
  label->SetTextSize(0.04);
  label->DrawLatex(0.1,0.21,Scaled_s);
  TString Sample_s = "Chiara"; 
  if(sample=="Man") Sample_s = "Manuel";
  label->DrawLatex(0.85,0.95,Sample_s);

  TString epsName = "babar_code/CSample/"; epsName += Variable; epsName += "_"; 
  epsName += sample; epsName += "_"; epsName += Channel; 
  if(AddConti)epsName += "_Cont"; 
  epsName += ".eps";
  c.SaveAs(epsName);
  hdata->Delete(); hdata2->Delete(); hMCsum->Delete(); hs->Delete();
  p1->Delete(); p2->Delete(); 
  for(int sam=0; sam<6; sam++)hStack[sam]->Delete();
 
}

void formatHisto(TH1F *h){
  h->SetMinimum(0);
  h->SetTitle("");
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.1);
  h->GetYaxis()->SetNdivisions(3+100*2);
  h->GetXaxis()->SetNdivisions(5+100*2);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.25);
  h->GetXaxis()->SetLabelSize(0.2);
  h->GetXaxis()->SetTitleOffset(1.);
  return;
}


TString RoundNumber(double n, int e, double d){
  if(d==0) return " - ";
  double neg = 1; if(n*d<0) neg = -1;
  double b = (int)(neg*n/d*pow(10.,(double)e)+0.5);
  b /= pow(10.,(double)e)*neg;
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(!result.Contains(".") && e != 0) result += ".";
  
  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<e-afterdot.Length(); i++)
    result += "0";
  return result;
}


