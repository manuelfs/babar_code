#include "babar_code/NewPhysics/HiggsPlot.hh"
#include "babar_code/Styles/Styles.cc"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TBox.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMatrixT.h"

#define Nsigma 6

HiggsPlot::HiggsPlot(TString NameFile) {
  SetMasses(1);  // Charged B

  // From hep-ph/0908.3470
  Vub[0] = 0.00352; Vub[1] = 0.00011;
  fB[0]  = 0.196;   fB[1]  = 0.011; // [GeV]

  // Coefficients found with RateCalc::Errors with 10,000 variations
  double Coef0[3][3] = {{0.2971, -3.2475, 16.9180}, {0.0165, 0.3152, 1.9629}, {-0.9279, -0.9567, 0.7891}};
  double Coef1[3][3] = {{0.2520, -0.2298,  0.6427}, {0.0028, 0.0289, 0.0846}, {-0.9457, -0.9854, 0.9036}};
  for(int isDs=0; isDs<2; isDs++)
    for(int isError=0; isError<3; isError++)
      for(int iC=0; iC<3; iC++)
	if(isDs) RDCoef[isDs][isError][iC] = Coef1[isError][iC];
	else     RDCoef[isDs][isError][iC] = Coef0[isError][iC];

  _isgS = 1; _varyRD = 1;

  ReadRD(NameFile);

  DecayName[0] = "TauNu";   DecayName[1] = "DTauNu"; 
  DecayName[2] = "DsTauNu"; DecayName[3] = "DDsTauNu"; 

}
void HiggsPlot::SetMasses(int isBm){
  if(isBm) {  // Charged B
    mB  = 5.2792; BLifeTime = 1.638e-12; 
  } else {   // Neutral B
    mB  = 5.2795; BLifeTime = 1.525e-12; 
  }
  BLifeTime /= hbar;
  _isBm = isBm;
}

void HiggsPlot::PlotChi2(){
  Styles style; style.setPadsStyle(1); 
  style.PadLeftMargin = 0.14; style.PadBottomMargin = 0.17; 
  style.yTitleOffset = 0.7;   style.xTitleOffset = 0.7; 
  style.TitleSize = 0.1; style.LabelSize = 0.1;
  style.applyStyle();
  TCanvas can("can","Likelihood scan");
  TPad *cPad = (TPad *)can.cd(0);
  TString xTitle[] = {"R(D)", "R(D*)"};

  double nSig = 4.2;
  int nBins = 60;
  _varyRD = 0;

  double RD[2], RDs[2], frac=0.018;
  Compute(0,RD,1); Compute(0,RDs,2);
  double Limit[2][2] = {{Measurement[1][0]-nSig*Measurement[1][1], Measurement[1][0]+nSig*Measurement[1][1]},
			{Measurement[2][0]-nSig*Measurement[2][1], Measurement[2][0]+nSig*Measurement[2][1]}};
  TH2F hChi2("hChi2","",nBins, Limit[0][0],Limit[0][1], nBins, Limit[1][0],Limit[1][1]);

  for(int iRD=1; iRD<=nBins; iRD++){
    for(int iRDs=1; iRDs<=nBins; iRDs++){
      RD[0] = hChi2.GetXaxis()->GetBinCenter(iRD);
      RDs[0] = hChi2.GetYaxis()->GetBinCenter(iRDs);
      hChi2.SetBinContent(iRD, iRDs, ProbChi2(3, 0, RD, RDs));
      //cout<<1-ProbChi2(3, 0, RD, RDs)<<endl;
    }
    //cout<<endl;
  }

  hChi2.SetXTitle(xTitle[0]); hChi2.SetYTitle(xTitle[1]); 
  //int colors[] = {kBlue+2, kBlue+1, kBlue-4, kBlue-7, kBlue-9, kBlue-10, 0};
  int colors[] = {kBlue-3, kBlue-4, kBlue-7, kBlue-9, kBlue-10, 0};
  gStyle->SetPalette(Nsigma, colors);
  double zCont[Nsigma];
  for(int ns=1; ns<=Nsigma; ns++) zCont[ns] = IntG(0,1,-ns,ns);//cout<<endl<<zCont[ns]<<endl;}
  hChi2.SetContour(Nsigma,zCont);
  hChi2.GetXaxis()->CenterTitle(true); hChi2.GetYaxis()->CenterTitle(true);
  hChi2.SetLabelOffset(0.009,"xy");
  hChi2.SetLabelSize(0.09,"xy");
  hChi2.SetTitleOffset(1,"x");
  hChi2.SetTitleOffset(0.6,"y");
  hChi2.Draw("cont4");

  double padL[2][2] = {{cPad->GetLeftMargin(),   1-cPad->GetRightMargin()}, // Margins are reversed for x-y!
		       {cPad->GetBottomMargin(), 1-cPad->GetTopMargin()}};
  double SMyield[2];
  Compute(0,RD,1); Compute(0,RDs,2); RD[1] = RDs[0];
  TLine line; line.SetLineWidth(2); line.SetLineColor(1);
  for(int chan=0; chan<2; chan++){
    double range = Limit[chan][1]-Limit[chan][0];
    SMyield[chan] = padL[chan][0] + (padL[chan][1]-padL[chan][0])*(RD[chan]-Limit[chan][0])/range;
  }
  line.DrawLineNDC(SMyield[0]-frac, SMyield[1], SMyield[0]+frac, SMyield[1]);
  line.DrawLineNDC(SMyield[0], SMyield[1]-2*frac, SMyield[0], SMyield[1]+2*frac);

  TLatex latex; latex.SetNDC(kTRUE); latex.SetTextAlign(33); latex.SetTextSize(style.TitleSize);
  latex.DrawLatex(SMyield[0]-frac,SMyield[1]-frac,"SM");

  line.SetLineColor(0);
  for(int chan=0; chan<2; chan++){
    double range = Limit[chan][1]-Limit[chan][0];
    SMyield[chan] = padL[chan][0] + (padL[chan][1]-padL[chan][0])*(Measurement[1+chan][0]-Limit[chan][0])/range;
  }
  line.DrawLineNDC(SMyield[0]-frac, SMyield[1], SMyield[0]+frac, SMyield[1]);
  line.DrawLineNDC(SMyield[0], SMyield[1]-2*frac, SMyield[0], SMyield[1]+2*frac);
  
  TH1F *histo[Nsigma];
  double legW = 0.08, legH = 0.07*Nsigma;
  double legX = 1-style.PadRightMargin-0.02, legY = 1-style.PadTopMargin-0.02;
  TLegend *leg = new TLegend(legX-legW, legY-legH, legX, legY);
  leg->SetTextSize(0.08); leg->SetFillColor(0); leg->SetTextFont(style.nFont);
  for(int ileg=0; ileg<Nsigma-1; ileg++) {
    TString label = "histo"; label += ileg+1; 
    histo[ileg] = new TH1F(label,"histo",10,0,10);
    histo[ileg]->SetLineColor(colors[ileg]);histo[ileg]->SetFillColor(colors[ileg]);
    label = " "; label += ileg+1; label += "#sigma";
    leg->AddEntry(histo[ileg],label);
  }
  leg->Draw(); 

  
  TString pName = "public_html/Higgs_Chi2.eps"; 
  can.SaveAs(pName);

  for(int ileg=0; ileg<Nsigma-1; ileg++) histo[ileg]->Delete();
}

void HiggsPlot::PlotPRL(int isPsfrag){
  double bottom = 0.12;
  Styles style; style.setPadsStyle(-2); 
  style.PadBottomMargin = 2*bottom/(1+bottom);
  style.applyStyle();

  int nBins = 1000, iDecay, colors[2][2] = {{2,1},{kBlue-10,kBlue+1}};
  int iniHig = 0, iHig = iniHig, finHig = 100, dHig = 5, bin, nFiles = (finHig-iniHig)/dHig+1;
  double tBmH_max = 1, PadLimit[2][2] = {{(1-bottom)/2.,1}, {0, (1+bottom)/2.}};
  double maxRD[] = {0.82, 0.47}, minRD[] = {0.08, 0.17}, TopMargin[] = {0.01, 0};
  TString hName, epsName = "public_html/PRL_Higgs.eps", label, padName;
  TString TagDecay[] = {"R(D)", "R(D*)"};
  
  TBox box; box.SetLineColor(4);box.SetFillColor(4);
  TLine line; line.SetLineStyle(1); line.SetLineColor(4); line.SetLineWidth(3);
  TCanvas can("can","RD Vs Higgs");
  TPad *Pads[2];
  TH1F *hBF[2][2], *hMeas[2][2];
  for(int his=0; his<1; his++) {
    for(int isDs=0; isDs<2; isDs++){
      hName = "hBF"; hName += his; hName += isDs;
      hBF[isDs][his] = new TH1F(hName,"",nBins,0,tBmH_max);
      hBF[isDs][his]->SetFillColor(colors[0][his]); hBF[isDs][his]->SetLineColor(colors[0][his]);      
      hBF[isDs][his]->SetLineWidth(2);
      hName = "hMeas"; hName += his; hName += isDs;
      hMeas[isDs][his] = new TH1F(hName,"",nBins,0,tBmH_max);
      hMeas[isDs][his]->SetFillColor(colors[1][his]); 
      hMeas[isDs][his]->SetLineColor(colors[1][his]); hMeas[isDs][his]->SetLineWidth(3);
    }
  }
  double RD[2][30], errorRD[2][30], List_tBmH[30], SysUncert[] = {0.095, 0.053};
  TString folder = "FitAll/fits/TextFinalIsoDataHi2x", text = ""; 
  for(int file=0; file<nFiles; file++){
    TString fileName = folder; 
    if(iHig<100) fileName += "0";
    if(iHig<10) fileName += "0";
    fileName += iHig; fileName += ".txt";
    fstream textFile; textFile.open(fileName,fstream::in);
    bin = 0;
    while(!text.Contains("R(D")) {
      textFile>>text; 
      bin++;
      if(bin>1000) {cout<<"R(D) not found in "<<fileName<<endl; return;}
    }
    for(int cand=0; cand<2; cand++){
      if(cand>0) textFile>>text;
      textFile>>text; RD[cand][file] = text.Atof();
      textFile>>text; textFile>>text; text.ReplaceAll(",",""); errorRD[cand][file] = text.Atof();
      for(int i=0; i<12; i++) textFile>>text; 
      //cout<<RD[file][cand]<<" #pm "<<errorRD[cand][file]<<endl;
      errorRD[cand][file] = sqrt(pow(errorRD[cand][file],2)+pow(RD[cand][file]*SysUncert[cand],2)); // Syst. error added
      List_tBmH[file] = (double)iHig/100.;
    }
    //cout<<endl;
    iHig += dHig;
  }

  double tBmH, BF[2];
  for(int isDs=0; isDs<2; isDs++){
    can.cd(0);
    padName = "Pad"; padName += isDs;
    Pads[isDs] = new TPad(padName,"",0,PadLimit[isDs][0], 1, PadLimit[isDs][1]); 
    Pads[isDs]->SetTopMargin(TopMargin[isDs]);
    Pads[isDs]->Draw(); Pads[isDs]->cd();
    iDecay = isDs+1;
    for(int bin=1; bin<=nBins; bin++){
      tBmH = hBF[isDs][0]->GetBinCenter(bin);
      Compute(tBmH,BF,iDecay);
      hBF[isDs][0]->SetBinContent(bin, BF[0]);
      hBF[isDs][0]->SetBinError(bin, BF[1]);
      if(_varyRD){
	hMeas[isDs][0]->SetBinContent(bin, MeasuredRD[isDs]->Eval(hMeas[isDs][0]->GetBinCenter(bin)));
	hMeas[isDs][0]->SetBinError(bin, MeasuredRD[isDs+2]->Eval(hMeas[isDs][0]->GetBinCenter(bin)));
      } else {
	hMeas[isDs][0]->SetBinContent(bin, Measurement[isDs+1][0]);
	hMeas[isDs][0]->SetBinError(bin, Measurement[isDs+1][1]);
      }
    }
    hName += isDs;
    hBF[isDs][1] = (TH1F*)hBF[isDs][0]->Clone(hName);
    hName += isDs;
    hMeas[isDs][1] = (TH1F*)hMeas[isDs][0]->Clone(hName);
    for(int bin=1; bin<=nBins; bin++){hBF[isDs][1]->SetBinError(bin,0);hMeas[isDs][1]->SetBinError(bin,0);}
    //hMeas[isDs][0]->SetFillStyle(3002);
    hBF[isDs][1]->SetFillColor(0); hBF[isDs][1]->SetLineColor(colors[0][1]); 
    hBF[isDs][1]->SetLineWidth(2);
    hMeas[isDs][1]->SetFillColor(0); hMeas[isDs][1]->SetLineColor(colors[1][1]);      
    hMeas[isDs][1]->SetLineWidth(3);
    hBF[isDs][0]->GetXaxis()->CenterTitle(true); hBF[isDs][0]->GetYaxis()->CenterTitle(true); 
    hBF[isDs][0]->SetMinimum(minRD[isDs]); hBF[isDs][0]->SetMaximum(maxRD[isDs]);
    hBF[isDs][0]->Draw("e3");

    hMeas[isDs][0]->Draw("e3 same"); hMeas[isDs][1]->Draw("c same");
    hBF[isDs][0]->Draw("e3 same");   hBF[isDs][1]->Draw("c same");
    TString xTitle = "tan#beta/m_{H^{+}} (GeV^{-1})";
    if(isPsfrag) xTitle = "t";
    if(isDs==0) xTitle = "";
    style.setTitles(hBF[isDs][0],xTitle,TagDecay[isDs]);
    hBF[isDs][0]->Draw("axis same");
  }
  if(isPsfrag) epsName.ReplaceAll(".eps","psfrag.eps");
  can.SaveAs(epsName);
  for(int isDs=0; isDs<2; isDs++){
    //Pads[isDs]->Delete();
    for(int his=0; his<2; his++) {
      hBF[isDs][his]->Delete();
      hMeas[isDs][his]->Delete();
    }
  }
}

void HiggsPlot::PlotgSLPRL(int isgSR, double tBmH_max){
  double b = 0.12;
  Styles style; style.setPadsStyle(-2); style.TitleSize = 0.1;
  style.PadBottomMargin = 2*b/(1+b);
  style.applyStyle();
  int nBins = 1000, iDecay, colors[2][2] = {{2,1},{kBlue-10,kBlue+1}};
  double PadLimit[2][2] = {{(1-b)/2.,1}, {0, (1+b)/2.}};
  double maxRD[] = {0.82, 0.47}, minRD[] = {0.1, 0.21}, TopMargin[] = {0.01, 0};
  TString hName, epsName = "public_html/Conclusions_gSR.eps", label, padName;
  //TString TagDecay[] = {"R(D)", "R(D*)"}, xTitle = "-g_{SR} (GeV^{-1})";
  TString TagDecay[] = {"R(D)", "R(D*)"}, xTitle = "g";
  TBox box; box.SetLineColor(4);box.SetFillColor(4);
  TLine line; line.SetLineStyle(1); line.SetLineColor(4); line.SetLineWidth(3);
  TCanvas can("can","RD Vs Higgs");
  TPad *Pads[2];
  TH1F *hBF[3][2], *hMeas[3][2];
  if(isgSR==0){
    colors[0][0] = kGreen+1;
    RDCoef[1][0][1] = -RDCoef[1][0][1];
    epsName.ReplaceAll("gSR","gSL");
    //xTitle = "-g_{SL} (GeV^{-1})";
  }
  for(int his=0; his<1; his++) {
    for(int isDs=0; isDs<3; isDs++){
      hName = "hBF"; hName += his; hName += isDs;
      hBF[isDs][his] = new TH1F(hName,"",nBins,0,tBmH_max);
      hBF[isDs][his]->SetFillColor(colors[0][his]); hBF[isDs][his]->SetLineColor(colors[0][his]);      
      hBF[isDs][his]->SetLineWidth(2);
      hName = "hMeas"; hName += his; hName += isDs;
      hMeas[isDs][his] = new TH1F(hName,"",nBins,0,tBmH_max);
      hMeas[isDs][his]->SetFillColor(colors[1][his]); hMeas[isDs][his]->SetLineColor(colors[1][his]);      
      hMeas[isDs][his]->SetLineWidth(3);
    }
  }

  double tBmH, BF[2], gSL, mTau = 1.7768;
  for(int isDs=0; isDs<2; isDs++){
    can.cd(0);
    padName = "Pad"; padName += isDs;
    Pads[isDs] = new TPad(padName,"",0,PadLimit[isDs][0], 1, PadLimit[isDs][1]); 
    Pads[isDs]->SetTopMargin(TopMargin[isDs]);
    Pads[isDs]->Draw(); Pads[isDs]->cd();
    iDecay = isDs+1;
    for(int bin=1; bin<=nBins; bin++){
      gSL = hBF[isDs][0]->GetBinCenter(bin)/mTau;
      tBmH = sqrt(gSL/4.2);
      Compute(tBmH,BF,iDecay);
      hBF[isDs][0]->SetBinContent(bin, BF[0]);
      hBF[isDs][0]->SetBinError(bin, BF[1]);
      if(_varyRD){
	hMeas[isDs][0]->SetBinContent(bin, MeasuredRD[isDs]->Eval(tBmH));
	hMeas[isDs][0]->SetBinError(bin, MeasuredRD[isDs+2]->Eval(tBmH));
      } else {
	hMeas[isDs][0]->SetBinContent(bin, Measurement[isDs+1][0]);
	hMeas[isDs][0]->SetBinError(bin, Measurement[isDs+1][1]);
      }
    }

    hName += isDs;
    hBF[isDs][1] = (TH1F*)hBF[isDs][0]->Clone(hName);
    hName += isDs;
    hMeas[isDs][1] = (TH1F*)hMeas[isDs][0]->Clone(hName);
    for(int bin=1; bin<=nBins; bin++){hBF[isDs][1]->SetBinError(bin,0);hMeas[isDs][1]->SetBinError(bin,0);}
    //hMeas[isDs][0]->SetFillStyle(3002);
    hBF[isDs][1]->SetFillColor(0); hBF[isDs][1]->SetLineColor(colors[0][1]); 
    hBF[isDs][1]->SetLineWidth(2);
    hMeas[isDs][1]->SetFillColor(0); hMeas[isDs][1]->SetLineColor(colors[1][1]);      
    hMeas[isDs][1]->SetLineWidth(3);


    hBF[isDs][0]->GetXaxis()->CenterTitle(true); hBF[isDs][0]->GetYaxis()->CenterTitle(true); 
    hBF[isDs][0]->SetMinimum(minRD[isDs]);
    hBF[isDs][0]->SetMaximum(maxRD[isDs]);
    hBF[isDs][0]->Draw("e3");

    hMeas[isDs][0]->Draw("e3 same"); hMeas[isDs][1]->Draw("c same");
    hBF[isDs][0]->Draw("e3 same");
    hBF[isDs][1]->Draw("c same");
    if(0){
      RDCoef[isDs][0][1] = -RDCoef[isDs][0][1];
      for(int bin=1; bin<=nBins; bin++){
	gSL = hBF[isDs+1][0]->GetBinCenter(bin);
	tBmH = sqrt(gSL/4.2);
	Compute(tBmH,BF,iDecay);
	hBF[isDs+1][0]->SetBinContent(bin, BF[0]);
	hBF[isDs+1][0]->SetBinError(bin, BF[1]);
      }
      hBF[isDs+1][0]->SetFillColor(8); hBF[isDs][0]->SetLineColor(8);
      hName += isDs+1;
      hBF[isDs+1][1] = (TH1F*)hBF[isDs+1][0]->Clone(hName);
      for(int bin=1; bin<=nBins; bin++)hBF[isDs+1][1]->SetBinError(bin,0);
      hBF[isDs+1][1]->SetFillColor(0);hBF[isDs+1][1]->SetLineColor(1); hBF[isDs+1][1]->SetLineWidth(2);

      hBF[isDs][0]->SetMinimum(minRD[isDs]);
      hBF[isDs][0]->SetMaximum(maxRD[isDs]);
      hBF[isDs+1][0]->Draw("e3");
      hMeas[isDs][0]->Draw("e3 same"); hMeas[isDs][1]->Draw("c same");
      hBF[isDs+1][0]->Draw("e3 same");
      hBF[isDs+1][0]->Draw("e3 same");
      hBF[isDs+1][1]->Draw("c same");
      RDCoef[isDs][0][1] = -RDCoef[isDs][0][1];
    }
    TString Title = xTitle;
    if(isDs==0) Title = "";
    style.setTitles(hBF[isDs][0],Title,TagDecay[isDs]);
    hBF[isDs][0]->Draw("axis same");
  }
  can.SaveAs(epsName);
  for(int isDs=0; isDs<3; isDs++){
    //Pads[isDs]->Delete();
    for(int his=0; his<2; his++) {
      if(isDs==2 && his == 1) continue;
      hBF[isDs][his]->Delete();
    }
    hMeas[isDs][0]->Delete();
  }
}

void HiggsPlot::PlotgSLTauNu(double b, double tBmH_max){
  Styles style; style.setPadsStyle(-2); 
  style.PadBottomMargin = 2*b/(1+b);
  style.applyStyle();
  int nBins = 1000;
  double PadLimit[2][2] = {{(1-b)/2.,1}, {0, (1+b)/2.}};
  double maxRD[] = {45, 45}, minRD[] = {0., 0.}, TopMargin[] = {0.01, 0};
  TString hName, epsName = "public_html/Higgs_TauNu_gSL.eps", label, padName;
  TString TagDecay[] = {"BF(B#rightarrow#tau#nu) (10^{-5})", "BF(B#rightarrow#tau#nu) (10^{-5})"};
  TBox box; box.SetLineColor(4);box.SetFillColor(4);
  TLine line; line.SetLineStyle(1); line.SetLineColor(4); line.SetLineWidth(3);
  TCanvas can("can","RD Vs Higgs");
  TPad *Pads[2];
  TH1F *hBF[4][2];
  for(int his=0; his<1; his++) {
    for(int isDs=0; isDs<4; isDs++){
      hName = "hBF"; hName += his; hName += isDs;
      hBF[isDs][his] = new TH1F(hName,"",nBins,0,tBmH_max);
    }
  }
  Vub[1] = 0.0003;
  double tBmH, BF[2], gSL;
  for(int isDs=0; isDs<2; isDs++){
    can.cd(0);
    if(isDs==0) Vub[0] = 0.00313; else Vub[0] = 0.00431;
    padName = "Pad"; padName += isDs;
    Pads[isDs] = new TPad(padName,"",0,PadLimit[isDs][0], 1, PadLimit[isDs][1]); 
    Pads[isDs]->SetTopMargin(TopMargin[isDs]);
    Pads[isDs]->Draw(); Pads[isDs]->cd();
    for(int bin=1; bin<=nBins; bin++){
      gSL = hBF[isDs][0]->GetBinCenter(bin);
      tBmH = sqrt(gSL/4.2);
      Compute(tBmH,BF,0);
      hBF[isDs][0]->SetBinContent(bin, BF[0]);
      hBF[isDs][0]->SetBinError(bin, BF[1]);
    }
    hBF[isDs][0]->SetFillColor(2); hBF[isDs][0]->SetLineColor(2);
    hName += isDs;
    hBF[isDs][1] = (TH1F*)hBF[isDs][0]->Clone(hName);
    for(int bin=1; bin<=nBins; bin++)hBF[isDs][1]->SetBinError(bin,0);
    hBF[isDs][1]->SetFillColor(0);hBF[isDs][1]->SetLineColor(1); hBF[isDs][1]->SetLineWidth(2);

    hBF[isDs][0]->SetMinimum(minRD[isDs]);
    hBF[isDs][0]->SetMaximum(maxRD[isDs]);
    hBF[isDs][0]->Draw("e3");

    box.SetFillStyle(3002); 
    box.DrawBox(0,Measurement[0][0]-Measurement[0][1],
		tBmH_max,Measurement[0][0]+Measurement[0][1]);
    line.DrawLine(0,Measurement[0][0],tBmH_max,Measurement[0][0]);
    hBF[isDs][0]->Draw("e3 same");
    hBF[isDs][1]->Draw("c same");
      _isgS = -1;
      for(int bin=1; bin<=nBins; bin++){
	gSL = hBF[isDs+2][0]->GetBinCenter(bin);
	tBmH = sqrt(gSL/4.2);
	Compute(tBmH,BF,0);
	hBF[isDs+2][0]->SetBinContent(bin, BF[0]);
	hBF[isDs+2][0]->SetBinError(bin, BF[1]);
      }
      hBF[isDs+2][0]->SetFillColor(8); hBF[isDs][0]->SetLineColor(8);
      hName += isDs+2;
      hBF[isDs+2][1] = (TH1F*)hBF[isDs+2][0]->Clone(hName);
      for(int bin=1; bin<=nBins; bin++)hBF[isDs+2][1]->SetBinError(bin,0);
      hBF[isDs+2][1]->SetFillColor(0);hBF[isDs+2][1]->SetLineColor(1); hBF[isDs+2][1]->SetLineWidth(2);

      hBF[isDs][0]->SetMinimum(minRD[isDs]);
      hBF[isDs][0]->SetMaximum(maxRD[isDs]);
      hBF[isDs+2][0]->Draw("e3 same");
      hBF[isDs+2][0]->Draw("e3 same");
      hBF[isDs+2][1]->Draw("c same");
      _isgS = 1;

    style.setTitles(hBF[isDs][0],"-g (GeV^{-1})",TagDecay[isDs]);
  }
  can.SaveAs(epsName);
  for(int isDs=0; isDs<4; isDs++){
    //Pads[isDs]->Delete();
    for(int his=0; his<2; his++) hBF[isDs][his]->Delete();
  }
}

void HiggsPlot::PlotBF(int iDecay, double tBmH_max, double BF_max){
  if(iDecay<0 || iDecay>2) {cout<<"iDecay must be 0, 1 or 2"<<endl; return;}
  Styles style; style.setPadsStyle(-1); style.applyStyle();
  int nBins = 1000;
  TString hName, epsName = "public_html/Higgs_BF_TEMP_BaBar.eps", label, Llabel, Rlabel;
  TString yTitle[] = {"BF(B#rightarrow#tau#nu) (10^{-5})", "R(D)", "R(D*)"};
  TString TagDecay[] = {"BF", "R(D)", "R(D*)"};
  TCanvas can;
  TH1F *hBF[2];
  for(int his=0; his<1; his++) {
    hName = "hBF"; hName += his;
    hBF[his] = new TH1F(hName,"",nBins,0,tBmH_max);
  }
  double tBmH, BF[2];
  for(int bin=1; bin<=nBins; bin++){
    tBmH = hBF[0]->GetBinCenter(bin);
    Compute(tBmH,BF,iDecay);
    hBF[0]->SetBinContent(bin, BF[0]);
    hBF[0]->SetBinError(bin, BF[1]);
  }
  hBF[0]->SetFillColor(2); hBF[0]->SetLineColor(2);
  hName += "1";
  hBF[1] = (TH1F*)hBF[0]->Clone(hName);
  for(int bin=1; bin<=nBins; bin++)hBF[1]->SetBinError(bin,0);
  hBF[1]->SetFillColor(0);hBF[1]->SetLineColor(1); hBF[1]->SetLineWidth(2);

  TBox box; box.SetLineColor(4);box.SetFillColor(4);
  TLine line; line.SetLineStyle(1); line.SetLineColor(4); line.SetLineWidth(3);
  if(BF_max>0) hBF[0]->SetMaximum(BF_max);
  hBF[0]->Draw("e3");
  box.SetFillStyle(3002); 
  box.DrawBox(0,Measurement[iDecay][0]-Measurement[iDecay][1],
	      tBmH_max,Measurement[iDecay][0]+Measurement[iDecay][1]);
  line.DrawLine(0,Measurement[iDecay][0],tBmH_max,Measurement[iDecay][0]);
  hBF[0]->Draw("e3 same");
  hBF[1]->Draw("c same");

  Compute(0,BF,iDecay);
  label = "#splitline{"; label += TagDecay[iDecay]; label += "_{SM} = ";
  if(iDecay==0){
    label+="(";label+=RoundNumber(BF[0],1); label+=" #pm ";
    label+=RoundNumber(BF[1],1); label+=")#times10^{-5}}{BF_{exp} = (";
    label+=RoundNumber(Measurement[iDecay][0],1); label+=" #pm ";
    label+=RoundNumber(Measurement[iDecay][1],1); label+=")#times10^{-5}}";
    Llabel = ""; Rlabel = label;
  } else {
    label+=RoundNumber(BF[0],3); label+=" #pm ";
    label+=RoundNumber(BF[1],3); label+="}{"; label += TagDecay[iDecay]; label += "_{exp} = ";
    label+=RoundNumber(Measurement[iDecay][0],3); label+=" #pm ";
    label+=RoundNumber(Measurement[iDecay][1],3); label+="}";
    Rlabel = ""; Llabel = label;
  }
  style.setTitles(hBF[0],"tan#beta/m_{H^{+}} (GeV^{-1})",yTitle[iDecay],Llabel,Rlabel);
  epsName.ReplaceAll("TEMP",DecayName[iDecay]);
  can.SaveAs(epsName);
  for(int his=0; his<2; his++) hBF[his]->Delete();
}

void HiggsPlot::PlotExclusion(int iDecay, double tBmH_max, TString Option){
  //_varyRD = 0;
  if(iDecay<0 || iDecay>3) {cout<<"iDecay must be 0, 1, 2 or 3"<<endl; return;}
  Styles style; style.setPadsStyle(-1); style.applyStyle();
  int nBins = 1000, nBins2D = 200;
  double LikelytBmH=-1, Pvalue, Pmin=1;
  TString hName, epsName = "public_html/Higgs_Exclusion1D_TEMP_BaBar.eps",label;
  TCanvas can("can","Exclusion 2HDM");
  TH1F h1D("h1D","",nBins,0,tBmH_max);
  h1D.GetXaxis()->SetTitle("tan#beta/m_{H^{+}} (GeV^{-1})");
  h1D.GetYaxis()->SetTitle("Exclusion probability  (%)");
  double tBmH;
  for(int bin=1; bin<=nBins; bin++){
    tBmH = h1D.GetBinCenter(bin);
    Pvalue = ProbChi2(iDecay, tBmH);
    if(Pvalue<Pmin) {LikelytBmH = tBmH; Pmin=Pvalue;}
    h1D.SetBinContent(bin, 100*Pvalue);
  }
  h1D.Draw("c");
  Pvalue = ProbChi2(iDecay, 0);
  label = "#splitline{SM excluded at "; label += RoundNumber(sqrt(2)*TMath::ErfInverse(Pvalue),1); 
  label +=" #sigma}{"; label += "Most likely tan#beta/m_{H^{+}} = "; 
  label += RoundNumber(LikelytBmH,2); label += "}";
  TLatex latex; latex.SetNDC(kTRUE); latex.SetTextAlign(31); latex.SetTextSize(style.TextSize);
  double xLabel = 1-style.PadRightMargin-0.02; 
  if(iDecay==2) {
    xLabel = style.PadLeftMargin+0.03;
    latex.SetTextAlign(11);
  }
  latex.DrawLatex(xLabel,style.PadBottomMargin+0.1,label);

  //////////////////    Plot 2D  ///////////////////////
  int nSigmas[] = {3,4,5, 1}, Nsig=3; double dl[5];
  for(int ns=0; ns<Nsig+1; ns++) dl[ns] = IntG(0,1,-nSigmas[ns],nSigmas[ns]);
  TLine line; line.SetLineStyle(2); line.SetLineColor(1); line.SetLineWidth(1);
  line.DrawLine(0,dl[3]*100,tBmH_max,dl[3]*100);
  line.DrawLine(0,dl[0]*100,tBmH_max,dl[0]*100);

  epsName.ReplaceAll("TEMP",DecayName[iDecay]);
  can.SaveAs(epsName);

  int colors[] = {kBlue-10, kBlue-9, kBlue-7, kBlue-4, 0};
  gStyle->SetPalette(Nsig, colors);
  TH2F h2D("h2D","",nBins2D,0,1000,nBins2D,0,110);
  h2D.GetYaxis()->CenterTitle(true); h2D.GetXaxis()->CenterTitle(true);
  h2D.GetYaxis()->SetTitle("tan#beta"); h2D.GetXaxis()->SetTitle("m_{H^{+}} (GeV)");
  h2D.GetYaxis()->SetTitle("B"); h2D.GetXaxis()->SetTitle("M");
  for(int binH=1; binH<=nBins2D; binH++){
    for(int binB=1; binB<=nBins2D; binB++){
      tBmH = h2D.GetYaxis()->GetBinCenter(binB)/h2D.GetXaxis()->GetBinCenter(binH);
      h2D.SetCellContent(binH, binB, ProbChi2(iDecay, tBmH));
    }
  }
  h2D.SetContour(Nsig,dl);
  h2D.Draw(Option);
  TH1F *histo[4];
  double legW = 0.2, legH = 0.28;
  double legX = 1-style.PadRightMargin-0.02, legY = style.PadBottomMargin+legH+0.04;
  TLegend *leg = new TLegend(legX-legW, legY-legH, legX, legY);
  leg->SetTextSize(0.075); leg->SetFillColor(0); leg->SetTextFont(style.nFont);
  for(int ileg=0; ileg<Nsig; ileg++) {
    TString label = "histo"; label += ileg+1; 
    histo[ileg] = new TH1F(label,"histo",10,0,10);
    histo[ileg]->SetLineColor(colors[ileg]);histo[ileg]->SetFillColor(colors[ileg]);
    label = "Excl. at "; label += nSigmas[ileg]; label += "#sigma";
    leg->AddEntry(histo[ileg],label);
  }
  leg->Draw();
  epsName.ReplaceAll("1D","2D");
  can.SaveAs(epsName);
  for(int his=0; his<Nsig; his++) histo[his]->Delete();
}

double HiggsPlot::ProbChi2(int iDecay, double tBmH){
  double RD[2], RDs[2], chi2;
  double valMeas[3][2] = {{Measurement[0][0], Measurement[1][0]},
			  {MeasuredRD[0]->Eval(tBmH), MeasuredRD[2]->Eval(tBmH)},
			  {MeasuredRD[1]->Eval(tBmH), MeasuredRD[3]->Eval(tBmH)}};
  int ndof=1;
  if(iDecay<3) {
    Compute(tBmH,RD,iDecay);
    chi2 = Chi2(RD,valMeas[iDecay]);
  } else {
    Compute(tBmH,RD,1); Compute(tBmH,RDs,2);
    //if(tBmH==0) {RD[0] = 0.316; RD[1] = 0.014;}  // MILC 2012
    //if(tBmH==0) {RD[0] = 0.302; RD[1] = 0.016;}  // Tanaka 2010
    //if(tBmH==0) {RD[0] = 0.310; RD[1] = 0.020;}  // Nierste 2008
    TMatrixT<double> CovM(2,2), DiffRD(1,2), MDiff(1,2); //CovM is the covariance matrix
    if(_varyRD){
      CovM(0,0) = pow(RD[1],2) +pow(valMeas[1][1],2);
      CovM(1,1) = pow(RDs[1],2)+pow(valMeas[2][1],2);
      CovM(1,0) = Measurement[3][0]*valMeas[1][1]*valMeas[2][1];
      DiffRD(0,0) = valMeas[1][0]-RD[0]; DiffRD(0,1) = valMeas[2][0]-RDs[0];
    } else {
      CovM(0,0) = pow(RD[1],2) +pow(Measurement[1][1],2);
      CovM(1,1) = pow(RDs[1],2)+pow(Measurement[2][1],2);
      CovM(1,0) = Measurement[3][0]*Measurement[1][1]*Measurement[2][1];
      DiffRD(0,0) = Measurement[1][0]-RD[0]; DiffRD(0,1) = Measurement[2][0]-RDs[0];
    }
    CovM(0,1) = CovM(1,0);
    CovM.Invert();
    MDiff.Mult(DiffRD,CovM);
    chi2 = MDiff(0,0)*DiffRD(0,0) + MDiff(0,1)*DiffRD(0,1);
    ndof++;
    if(tBmH==0 && _varyRD) 
      cout<<endl<<"SM chi2 is "<<RoundNumber(chi2,2)<<" with a p-value of "<<TMath::Prob(chi2,ndof)<<". This is "
	  <<RoundNumber(sqrt(2)*TMath::ErfInverse(1 - TMath::Prob(chi2,ndof)),2)<<" sigma away"<<endl<<endl;
  }

  return 1 - TMath::Prob(chi2,ndof);
}

double HiggsPlot::ProbChi2(int iDecay, double tBmH, double RD[2], double RDs[2]){
  double chi2;
  double valMeas[3][2] = {{Measurement[0][0], Measurement[1][0]},
			  {MeasuredRD[0]->Eval(tBmH), MeasuredRD[2]->Eval(tBmH)},
			  {MeasuredRD[1]->Eval(tBmH), MeasuredRD[3]->Eval(tBmH)}};
  int ndof=1;
  if(iDecay<3) {
    chi2 = Chi2(RD,valMeas[iDecay]);
  } else {
    TMatrixT<double> CovM(2,2), DiffRD(1,2), MDiff(1,2); //CovM is the covariance matrix
    if(_varyRD){
      CovM(0,0) = pow(RD[1],2) +pow(valMeas[1][1],2);
      CovM(1,1) = pow(RDs[1],2)+pow(valMeas[2][1],2);
      CovM(1,0) = Measurement[3][0]*valMeas[1][1]*valMeas[2][1];
      DiffRD(0,0) = valMeas[1][0]-RD[0]; DiffRD(0,1) = valMeas[2][0]-RDs[0];
    } else {
      CovM(0,0) = pow(RD[1],2) +pow(Measurement[1][1],2);
      CovM(1,1) = pow(RDs[1],2)+pow(Measurement[2][1],2);
      CovM(1,0) = Measurement[3][0]*Measurement[1][1]*Measurement[2][1];
      DiffRD(0,0) = Measurement[1][0]-RD[0]; DiffRD(0,1) = Measurement[2][0]-RDs[0];
    }
    CovM(0,1) = CovM(1,0);
    CovM.Invert();
    MDiff.Mult(DiffRD,CovM);
    chi2 = MDiff(0,0)*DiffRD(0,0) + MDiff(0,1)*DiffRD(0,1);
    ndof++;
    if(tBmH==0 && _varyRD) 
      cout<<endl<<"SM chi2 is "<<RoundNumber(chi2,2)<<" with a p-value of "<<TMath::Prob(chi2,ndof)<<". This is "
	  <<RoundNumber(sqrt(2)*TMath::ErfInverse(1 - TMath::Prob(chi2,ndof)),2)<<" sigma away"<<endl<<endl;
  }

  return 1 - TMath::Prob(chi2,ndof);
}

double HiggsPlot::Chi2(double M1[2], double M2[2]){
  double totSigma = pow(M1[1],2)+pow(M2[1],2);

  return pow(M1[0]-M2[0],2)/totSigma;
}

void HiggsPlot::Compute(double tBmH, double RD[2], int iDecay){
    if(iDecay==0) BFTauNu(tBmH,RD);
    else RDCalc(tBmH,RD,iDecay-1);
}

void HiggsPlot::BFTauNu(double tBmH, double BF[2], double ml){
  if(_isBm!=1) {
    cout<<"B->TauNu not possible for neutral B"<<endl;
    BF[0] = 0; BF[1] = 0;
    return;
  }

  double gS = pow(mB*tBmH,2);
  double BF_SM = GF*GF*mB*ml*ml/(8*PI)*pow(1-ml*ml/mB/mB,2)*pow(fB[0],2)*pow(Vub[0],2)*BLifeTime;
  BF[0] = BF_SM*pow(1-_isgS*gS,2)*1e5;
  BF[1] = BF[0]*2*sqrt(pow(fB[1]/fB[0],2)+pow(Vub[1]/Vub[0],2));
}

void HiggsPlot::RDCalc(double tBmH, double RD[2], int isDs){
  RD[0] = RDCoef[isDs][0][0] + RDCoef[isDs][0][1]*pow(tBmH,2) + RDCoef[isDs][0][2]*pow(tBmH,4);
  if(_isgS!=1&&isDs) RD[0] = RDCoef[isDs][0][0] - RDCoef[isDs][0][1]*pow(tBmH,2) + RDCoef[isDs][0][2]*pow(tBmH,4);

  double CovM[3][3], dRD[3] = {1, pow(tBmH,2), pow(tBmH,4)};
  for(int iM=0; iM<3; iM++){
    CovM[iM][iM] = pow(RDCoef[isDs][1][iM],2);
    CovM[iM][(iM+1)%3] = RDCoef[isDs][2][iM]*RDCoef[isDs][1][iM]*RDCoef[isDs][1][(iM+1)%3];
    CovM[(iM+1)%3][iM] = CovM[iM][(iM+1)%3];
  }
  RD[1] = 0;
  for(int iM=0; iM<3; iM++)
    for(int jM=0; jM<3; jM++)
      RD[1] += dRD[iM]*CovM[iM][jM]*dRD[jM];
  RD[1] = sqrt(RD[1]);
}

double HiggsPlot::IntG(double mean, double sigma, double minX, double maxX){
  return (TMath::Erf((maxX-mean)/sigma/sqrt(2.))-TMath::Erf((minX-mean)/sigma/sqrt(2.)))/2.;
}

TString HiggsPlot::RoundNumber(double n, int e, double d){
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

void HiggsPlot::ReadRD(TString NameFile){
  fstream textFile; textFile.open(NameFile,fstream::in);
  TString dummy, Name, VarNames[]  = {"BFTauNu", "RD", "RDs", "Correl."};
  double Value, Error;

  while(textFile){
    textFile>>Name>>dummy>>Value>>dummy>>Error;
    for(int iVar=0; iVar<nVar; iVar++){
      if(Name == VarNames[iVar]) {
	Measurement[iVar][0] = Value;
	Measurement[iVar][1] = Error;
      }
    }
  }
  // Obtained with PlotEffi(2,5) from babar_code/NewPhysics/PlotRate.C
//   double dEffError[2][21] = {{0.00, 0.00, 0.00, 0.00, 0.23, 0.43, 0.68, 1.22, 1.79, 2.08, 2.16, 
// 			      2.15, 2.11, 2.06, 2.02, 1.98, 1.95, 1.92, 1.90, 1.88, 1.86}, 
// 			     {0.00, 0.00, 0.00, 0.03, 0.06, 0.10, 0.15, 0.21, 0.28, 0.37, 0.47, 
// 			      0.58, 0.70, 0.81, 0.92, 1.02, 1.11, 1.18, 1.24, 1.28, 1.31}};

  double RelError[2][21] = {{2.10, 2.10, 2.10, 2.10, 2.12, 2.17, 2.27, 2.57, 3.04, 3.31, 3.39, 
			     3.37, 3.34, 3.29, 3.25, 3.22, 3.19, 3.16, 3.14, 3.13, 3.11}, 
			    {1.36, 1.36, 1.36, 1.36, 1.36, 1.37, 1.37, 1.38, 1.40, 1.43, 1.47, 
			     1.52, 1.58, 1.65, 1.73, 1.80, 1.87, 1.93, 1.98, 2.02, 2.05}};


  int iniHig = 0, iHig = iniHig, finHig = 100, dHig = 5, bin, nFiles = (finHig-iniHig)/dHig+1;
  double higRD[4][30], List_tBmH[30], SysUncert[3][2] = {{0.096, 0.055}, {0.018, 0.012}, {0.044, 0.020}};
  for(int isDs=0; isDs<2; isDs++)
    for(int sys=0; sys<3; sys++) SysUncert[sys][isDs] *= SysUncert[sys][isDs];
  TString folder = "FitAll/fits/TextFinalIsoDataHi2x", text = ""; 
  for(int file=0; file<nFiles; file++){
    TString fileName = folder; 
    if(iHig<100) fileName += "0";
    if(iHig<10) fileName += "0";
    fileName += iHig; fileName += ".txt";
    fstream textFile; textFile.open(fileName,fstream::in);
    bin = 0;
    while(!text.Contains("R(D")) {
      textFile>>text; 
      bin++;
      if(bin>1000) {cout<<"R(D) not found in "<<fileName<<endl; return;}
    }
    for(int cand=0; cand<2; cand++){
      if(cand>0) textFile>>text;
      textFile>>text; higRD[cand][file] = text.Atof();
      textFile>>text; textFile>>text; text.ReplaceAll(",",""); higRD[cand+2][file] = text.Atof();
      for(int i=0; i<12; i++) textFile>>text; 
      // Total uncertainty = Stat + Syst* + Eff + PDFs
      double ErrSystStar2 = pow(higRD[cand][file],2)*(SysUncert[0][cand]-SysUncert[1][cand]-SysUncert[2][cand]);
      double ErrEff2 = pow(higRD[cand][file]*RelError[cand][file]/RelError[cand][0],2)*SysUncert[1][cand];
      double ErrPDFs2 = pow(higRD[cand][file]*RelError[cand][file]/RelError[cand][0],2)*SysUncert[2][cand];

//       cout<<higRD[cand][file]<<" +- "<<RoundNumber(higRD[cand+2][file],3)<<" +- "
// 	  <<RoundNumber(sqrt(ErrSystStar2),3)
// 	  <<" +- "<< RoundNumber(sqrt(ErrEff2),3)<<" +- "<< RoundNumber(sqrt(ErrPDFs2),3)
// 	  <<",\t  a "<<RoundNumber(sqrt(ErrEff2)*100,2,higRD[cand][file])<<" +- "
// 	  << RoundNumber(sqrt(ErrPDFs2)*100,2,higRD[cand][file])<<" % error"<<endl;

      higRD[cand+2][file] = sqrt(pow(higRD[cand+2][file],2) + ErrSystStar2 + ErrEff2 + ErrPDFs2); // Syst. errors added

      List_tBmH[file] = (double)iHig/100.;
    }
    //cout<<endl;
    iHig += dHig;
  }
  for(int isDs=0; isDs<4; isDs++){
    text = "SplineRD"; text += isDs;
    MeasuredRD[isDs] = new TSpline3(text,List_tBmH,higRD[isDs],nFiles,"ble1",0); 
  }
  
  TString fileName = "babar_code/NewPhysics/BaBar_RDx_2HDM.txt";
  fstream thdmFile; thdmFile.open(fileName,fstream::out);
  thdmFile<<endl<<"tanB/mH\t\t     R(D)\t\t     R(D*)"<<endl<<endl;

  int nPoints = 1000;
  for(int point=0; point<=nPoints; point++){
    double tBmH = point/(double)nPoints;
    thdmFile<<RoundNumber(tBmH,3)<<"\t\t"<<RoundNumber(MeasuredRD[0]->Eval(tBmH),3)<<" +- "
	    <<RoundNumber(MeasuredRD[2]->Eval(tBmH),3)<<" \t\t"<<RoundNumber(MeasuredRD[1]->Eval(tBmH),3)<<" +- "
	    <<RoundNumber(MeasuredRD[3]->Eval(tBmH),3)<<endl;

  }
  cout<<"Written "<<fileName<<endl;


  fileName = "babar_code/NewPhysics/BaBar_RDx_2HDM_Theory.txt";
  fstream thdmFile2; thdmFile2.open(fileName,fstream::out);
  thdmFile2<<endl<<"tanB/mH\t\t     R(D)\t\t     R(D*)"<<endl<<endl;

  double BF_D[2], BF_Ds[2];
  for(int point=0; point<=nPoints; point++){
    double tBmH = point/(double)nPoints;
    Compute(tBmH,BF_D,1); Compute(tBmH,BF_Ds,2);
    thdmFile2<<RoundNumber(tBmH,3)<<"\t\t"<<RoundNumber(BF_D[0],3)<<" +- "
	     <<RoundNumber(BF_D[1],3)<<" \t\t"<<RoundNumber(BF_Ds[0],3)<<" +- "
	     <<RoundNumber(BF_Ds[1],3)<<endl;
  }
  cout<<"Written "<<fileName<<endl;

}

HiggsPlot::~HiggsPlot() {
  for(int isDs=0; isDs<4; isDs++) if(MeasuredRD[isDs]) MeasuredRD[isDs]->Delete();
}
