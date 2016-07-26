#include "TMath.h"
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>

#define nRows 2
#define nPads 6
#define nMaxBins 20
using namespace TMath;
using namespace std;
using std::cout;
using std::endl;

void SumHistos(TH1F *hSum, TH1F *h1, TH1F *h2=0, TH1F *h3=0, TH1F *h4=0, TH1F *h5=0, TH1F *h6=0, TH1F *h7=0);
double CalcChi2(TH1F *h1, TH1F *h2);
void CopyHisto(TH1F *h1, TH1F *h2, int onlyError=0);
void FitBB(double Err_BB[nMaxBins][2][2], int row, TH1F *hData[2], TH1F *hComp[2][2], TH1F *hBB[2][2], 
	   TH1F *hDss[2][2], TH1F *hNorm[2], TH1F *hDtau[2][2], TH1F *hDstau[2][2]);
void FitDss(TH1F *hDssPDF[3][2], int row, TH1F *hData[2], TH1F *hComp[2][2], TH1F *hBB[2][2], 
	    TH1F *hDss[2][2], TH1F *hNorm[2], TH1F *hDtau[2][2], TH1F *hDstau[2][2]);
void FitSignal(TH1F *hData[2], TH1F *hTot[2], TH1F *hComp[2][2], TH1F *hDtau[2][2], TH1F *hDstau[2][2]);

void Chi2_Q2(int nReps=2, TString doFit="BkgDssSig", int firstCol=0){
  if(nReps>10000) {cout<<"Too many reps"<<endl; return;}
  TString Folder = "keys/eps/CSample/root/";
  TString fName[3] = {"Q2", "Q2_030", "Q2_045"}, hName;
  int nPDFs = 9, nCols = 3, nBins = 17;
  TH1F *hComp[3][2][2], *hData[3][2], *hYields[nPads], *hMC[nPads][10], *hTemp, *hDssPDF[3][2];
  TH1F *hBB[3][2][2], *hDss[3][2][2], *hNorm[3][2], *hDtau[3][2][2], *hDstau[3][2][2];
  TH1F *hTot[3][2];
  double Err_BB[nMaxBins][2][2], limQ2[] = {4, 12.5};
  TString histosName = "public_html/Test_Histograms"; histosName+=doFit; histosName+=".root";
  TFile saveHistos(histosName, "RECREATE");
  cout<<"Writing "<<histosName<<endl;
  
  // Reading BB uncertainties
  TString textName = "public_html/Uncertainties_Bkg.txt";
  fstream textFile; textFile.open(textName,fstream::in);
  for(int row=0; row<nRows; row++){
    for(int bin=0; bin<nBins; bin++)
      textFile>>hName>>Err_BB[bin][row][0]>>Err_BB[bin][row][1];
  }
  textFile.close();

  gStyle->SetStatW(0.38); gStyle->SetStatH(0.14); gStyle->SetOptStat(111110);
  TCanvas can("can","",1000,700); can.Divide(3,2);
  //TPad *cPad;
  // Reading D**(l/tau)nu histograms
  textName = "public_html/Uncertainties_BkgDss.txt";
  textFile.open(textName,fstream::in);
  for(int row=0; row<nRows; row++){
    for(int dss=0; dss<3; dss++){
      hName = "hDssPDF"; hName += dss; hName += row;
      hDssPDF[dss][row] = new TH1F(hName, "", nBins, limQ2[0], limQ2[1]);
    }
    for(int bin=1; bin<nBins; bin++){
      textFile>>hName;
      double val[2];
      for(int dss=0; dss<3; dss++){
	textFile>>val[0]>>val[1];
	hDssPDF[dss][row]->SetBinContent(bin, val[0]);
	//hDssPDF[dss][row]->SetBinError(bin, val[1]);
	hDssPDF[dss][row]->SetBinError(bin, 0);
      }
    }
    for(int dss=0; dss<3; dss++)
      hDssPDF[dss][row]->Scale(1/hDssPDF[dss][row]->Integral());
  }

  // Chi2 for each pad
  for(int col=firstCol; col<nCols; col++){ 
    histosName = Folder; histosName += fName[col]; histosName += ".root";
    TFile fHistos(histosName);
    for(int row=0; row<nRows; row++){
      int pad = nRows*col+row;
      int candi = 2+row;
      fHistos.cd(); 
      hName = "Data"; hName += candi;
      hTemp = (TH1F *)fHistos.Get(hName);
      hName = "Data_"; hName += pad;
      if(hTemp) hData[col][row]  = (TH1F*)hTemp->Clone(hName);
      else continue;
      hData[col][row]->SetDirectory(0);
      hName = "hYields";  hName += pad;
      hTemp = (TH1F *)fHistos.Get("hYields");
      if(hTemp) hYields[pad]  = (TH1F*)hTemp->Clone(hName);
      else continue; 
      for(int pdf=0; pdf<nPDFs; pdf++){
	int ipdf = pdf;
	if(pdf==8) ipdf = 7;
	if(pdf==7) ipdf = 8;
	hName = "MC"; hName += candi+1; hName += ipdf;
 	double yield = 0;
	yield += hYields[pad]->GetBinContent(ipdf+10*row);
	yield += hYields[pad]->GetBinContent(ipdf+10*row+20);
	if(yield<=0){
	  cout<<pdf<<": ipdf "<<ipdf<<", yield "<<yield<<", row "<<row<<", col "<<col<<endl;
	  hMC[pad][pdf]  = (TH1F*)hMC[pad][pdf-1]->Clone(hName);
	} else {
	  hTemp = (TH1F *)fHistos.Get(hName);
	  hName = "MC_";  hName += pad; hName += pdf;
	  if(hTemp) hMC[pad][pdf]  = (TH1F*)hTemp->Clone(hName);
	  else {cout<<row<<", "<<col<<": PDF "<<pdf<<endl; return;}
	  if(hMC[pad][pdf]->Integral()<=0) {
	    cout<<pdf<<": ipdf "<<ipdf<<", yield "<<yield<<", row "<<row<<", col "<<col<<", candi "<<candi<<endl;
	    hName = "MC"; hName += candi+1-2; hName += ipdf;
	    hTemp = (TH1F *)fHistos.Get(hName);
	    hMC[pad][pdf]->Add(hTemp);
	  }
	}
	hMC[pad][pdf]->SetDirectory(0);
	hMC[pad][pdf]->Scale(yield/hMC[pad][pdf]->Integral());
	if(pdf==0) hTot[col][row] = (TH1F*)hMC[pad][pdf]->Clone();
	else hTot[col][row]->Add(hMC[pad][pdf]);
      }
      // BB background - Histograms with last index 0 are to be kept constant
      hBB[col][row][0] = (TH1F*) hMC[pad][0]->Clone();
      for(int pdf=1; pdf<3; pdf++) hBB[col][row][0]->Add(hMC[pad][pdf]);
      hBB[col][row][1] = (TH1F*) hBB[col][row][0]->Clone();

      // D**lnu
      hDss[col][row][0] = (TH1F*) hMC[pad][4]->Clone();
      hDss[col][row][1] = (TH1F*) hMC[pad][4]->Clone();

      // Normalization
      hNorm[col][row] = (TH1F*) hMC[pad][3]->Clone();
      hNorm[col][row]->Add(hMC[pad][5]); hNorm[col][row]->Add(hMC[pad][6]);

      // Signal
      hDstau[col][row][0] = (TH1F*) hMC[pad][7]->Clone();
      hDstau[col][row][1] = (TH1F*) hMC[pad][7]->Clone();
      hDtau[col][row][0] = (TH1F*) hMC[pad][8]->Clone();
      hDtau[col][row][1] = (TH1F*) hMC[pad][8]->Clone();

      // Chi2
      hComp[col][row][0] = (TH1F*) hBB[col][row][0]->Clone();
      hComp[col][row][0]->Add(hDss[col][row][1]);
      hComp[col][row][0]->Add(hNorm[col][row]);
      hComp[col][row][1] = (TH1F*) hComp[col][row][0]->Clone();
      hComp[col][row][0]->SetDirectory(0);
      hComp[col][row][1]->SetDirectory(0);

      hComp[col][row][1]->Add(hDtau[col][row][1]);
      hComp[col][row][1]->Add(hDstau[col][row][1]);

    }

    //////////// Chi2 optimization /////////////////////////
    for(int rep=0; rep<nReps; rep++){
      cout<<"Doing repetition "<<rep+1<<" of "<<nReps<<endl;
      if(doFit.Contains("Sig")) FitSignal(hData[col], hTot[col], hComp[col], hDtau[col], hDstau[col]);
      for(int row=0; row<nRows; row++){
	if(doFit.Contains("Dss")) 
	  FitDss(hDssPDF, row, hData[col], hComp[col], hBB[col], hDss[col], hNorm[col], hDtau[col], hDstau[col]);
	if(doFit.Contains("Bkg")) 
	  FitBB(Err_BB, row, hData[col], hComp[col], hBB[col], hDss[col], hNorm[col], hDtau[col], hDstau[col]);

	SumHistos(hComp[col][row][1], hBB[col][row][1], hDss[col][row][1], hNorm[col][row], 
		  hDtau[col][row][1], hDstau[col][row][1]);
	int ndof = 14; if(row) ndof = 12;
	double chi2 = CalcChi2(hData[col][row], hComp[col][row][1]);
	double Pchi2 = TMath::Prob(chi2, ndof);      
	int digits = 1; if(col==2&&row==0) digits=4;
	cout<<col<<","<<row<<": Chi2 "<<RoundNumber(chi2,1)<<",\t p-value "<<RoundNumber(Pchi2*100,digits)<<endl;
	can.cd(col+nCols*row+1);

	hComp[col][row][1]->SetFillStyle(0); hComp[col][row][1]->SetLineColor(2);
	//hComp[col][row][0]->SetFillStyle(0); hComp[col][row][0]->SetLineColor(4);
	hData[col][row]->SetMarkerStyle(20); hData[col][row]->SetMarkerSize(0.6); 
	hData[col][row]->Draw(); 
	//hComp[col][row][0]->Draw("hist same"); 
	hComp[col][row][1]->Draw("hist same");
	//hData[col][row]->Draw("e0 same");
	saveHistos.cd();
	hName = "BBOri"; hName += col; hName += row;
	hBB[col][row][0]->Write(hName); 
	hName = "BBFit"; hName += col; hName += row;
	hBB[col][row][1]->Write(hName); 
	hName = "DssOri"; hName += col; hName += row;
	hDss[col][row][0]->Write(hName); 
	hName = "DssFit"; hName += col; hName += row;
	hDss[col][row][1]->Write(hName); 
	hName = "DtauOri"; hName += col; hName += row;
	hDtau[col][row][0]->Write(hName); 
	hName = "DtauFit"; hName += col; hName += row;
	hDtau[col][row][1]->Write(hName); 
	hName = "DstauOri"; hName += col; hName += row;
	hDstau[col][row][0]->Write(hName); 
	hName = "DstauFit"; hName += col; hName += row;
	hDstau[col][row][1]->Write(hName); 
	hName = "Data"; hName += col; hName += row;
	hData[col][row]->Write(hName); 
	hName = "MCFit"; hName += col; hName += row;
	hComp[col][row][1]->Write(hName); 
      }
    }
  } // Columns

  TString pName = "public_html/Test_Chi2_"; pName+= doFit; pName += ".eps"; 
  can.SaveAs(pName);
  saveHistos.Close(); 
}


////////// BB fitting ////////////////////////////////////////////////////////////////
void FitBB(double Err_BB[nMaxBins][2][2], int row, TH1F *hData[2], TH1F *hComp[2][2], TH1F *hBB[2][2], 
	   TH1F *hDss[2][2], TH1F *hNorm[2], TH1F *hDtau[2][2], TH1F *hDstau[2][2]) {
  SumHistos(hComp[row][1], hDss[row][1], hNorm[row], hDtau[row][1], hDstau[row][1]);
  hComp[row][1]->Add(hData[row],-1);
	
  for(int bin=1; bin<=hBB[row][1]->GetNbinsX(); bin++) {
    double errBin = Err_BB[bin-1][row][1]/sqrt(12)*2;
    double loError = Err_BB[bin-1][row][0]-errBin;
    double hiError = Err_BB[bin-1][row][0]+errBin;
    if(loError > 1-errBin) loError = 1-errBin;
    if(hiError < 1+errBin) hiError = 1+errBin;
    double val = hBB[row][0]->GetBinContent(bin);
    double fitVal = -hComp[row][1]->GetBinContent(bin);
    if(fitVal > hiError*val) fitVal = hiError*val;
    if(fitVal < loError*val) fitVal = loError*val;
    hBB[row][1]->SetBinContent(bin, fitVal);
    //cout<<bin<<": val "<<val<<", fitVal "<<fitVal<<", hiError "<<hiError<<", loError "<<loError<<endl;
  }
}


////////// D** fitting ////////////////////////////////////////////////////////////////
void FitDss(TH1F *hDssPDF[3][2], int row, TH1F *hData[2], TH1F *hComp[2][2], TH1F *hBB[2][2], 
	    TH1F *hDss[2][2], TH1F *hNorm[2], TH1F *hDtau[2][2], TH1F *hDstau[2][2]) {
  TH1F *hBest = (TH1F*) hDss[row][0]->Clone();
  int nYield = 21, nF = 11;
  double dY = 1/(double)(nYield-1), dF = 1/(double)(nF-1);
  double YDss[2][2] = {{274.31, 87.967}, {34.994, 15.619}}, minChi2 = 1e9;
  double nDss = hDss[row][0]->Integral(), bestF[3]={-1,-1,-1};
  double sigDsstau[] = {-0.2, 0.5}, sigDpipi = 0.3, sigY = 0.15;
  for(double fDssY=0; fDssY<=1.01; fDssY+=dY){
    for(double fDssTau=0; fDssTau<=1.01; fDssTau+=dF){
      for(double fDpipi=0; fDpipi<=1.01; fDpipi+=dF){
	CopyHisto(hDss[row][1], hDss[row][0]);
	double YDssTau = nDss*(sigDsstau[0]+fDssTau*(sigDsstau[1]-sigDsstau[0]))/(YDss[0][row]/YDss[1][row]+1);
	hDss[row][1]->Add(hDssPDF[1][row], YDssTau);        // Adding D**taunu 
	//double YDpipi  = nDss*(sigDss[1]*(-1+fDpipi*2));
	double YDpipi  = nDss*sigDpipi*fDpipi;
	hDss[row][1]->Add(hDssPDF[2][row], YDpipi);         // Adding D**(->Dpipi)lnu 
	double Ytot = nDss*(1+(sigY*(-1+fDssY*2)));
	hDss[row][1]->Scale(Ytot/hDss[row][1]->Integral()); // Scaling the yield

	SumHistos(hComp[row][1], hBB[row][1], hDss[row][1], hNorm[row], hDtau[row][1], hDstau[row][1]);
	double chi2 = CalcChi2(hData[row], hComp[row][1]);
	if(chi2<minChi2){
	  CopyHisto(hBest, hDss[row][1]);
	  minChi2 = chi2;
	  bestF[0] = fDssY; bestF[1] = fDssTau;  bestF[2] = fDpipi; 
	  //cout<<"Best D**lnu: "<<minChi2<<", "<<bestF[0]<<", "<<bestF[1]<<", "<<bestF[2]<<endl;
	}
      }
    }
  }
  cout<<"Best D**lnu: "<<bestF[0]<<", "<<bestF[1]<<", "<<bestF[2]<<endl;
  CopyHisto(hDss[row][1], hBest);
  CopyHisto(hDss[row][1], hDss[row][0],1);
}

////////// Signal fitting ////////////////////////////////////////////////////////////////
void FitSignal(TH1F *hData[2], TH1F *hTot[2], TH1F *hComp[2][2], TH1F *hDtau[2][2], TH1F *hDstau[2][2]){
  double nDstau[2], nDtau[2], fYSig[2], bestFSig[2] = {1, 1};
  int row = 1;
  nDstau[row] = hTot[row]->Integral()-hComp[row][0]->Integral()-hDtau[row][0]->Integral();
  fYSig[1] = nDstau[row]/hDstau[row][0]->Integral();
  if(fYSig[1]<0.3 || fYSig[1]>1.6)
    cout<<" Tot "<<hTot[row]->Integral()<<", Comp "<<hComp[row][0]->Integral()
	<<", Dtau fu "<<hDtau[row][0]->Integral()<<" -> nDstau "<<nDstau[row]<<", f "<<fYSig[1]<<endl;
  row = 0;
  nDtau[row] = hTot[row]->Integral()-hComp[row][0]->Integral()-hDstau[row][0]->Integral();
  fYSig[0] = nDtau[row]/hDtau[row][0]->Integral();
  if(fYSig[0]<0.3 || fYSig[0]>1.6) 
    cout<<" Tot "<<hTot[row]->Integral()<<", Comp "<<hComp[row][0]->Integral()
	<<", Dstau fd "<<hDstau[row][0]->Integral()<<" -> nDtau "<<nDtau[row]<<", f "<<fYSig[0]<<endl;

  double df = 0.01, inif = 0.4, finf = 1.6, minChi2 = 1e10;
  for(double fTau = inif*fYSig[0]; fTau<finf*fYSig[0]; fTau += df*fYSig[0]){
    for(double fsTau = inif*fYSig[1]; fsTau<finf*fYSig[1]; fsTau += df*fYSig[1]){
      double chi2=0;
      for(int row=0; row<nRows; row++){
	CopyHisto(hDtau[row][1],  hDtau[row][0]);	
	CopyHisto(hDstau[row][1], hDstau[row][0]);	
	hDtau[row][1]->Scale(fTau); 
	hDstau[row][1]->Scale(fsTau);

	CopyHisto(hComp[row][1], hComp[row][0]);
	hComp[row][1]->Add(hDtau[row][1]);
	hComp[row][1]->Add(hDstau[row][1]);
	chi2 += CalcChi2(hData[row], hComp[row][1]);      
      }
      if(minChi2>chi2){
	//cout<<minChi2<<" > "<<chi2<<" - fTau "<<fTau<<", fsTau "<<fsTau<<endl;
	minChi2 = chi2;
	bestFSig[0] = fTau; bestFSig[1] = fsTau;
      }
    }
  }

  for(int row=0; row<nRows; row++){
    CopyHisto(hDtau[row][1],  hDtau[row][0]);	
    CopyHisto(hDstau[row][1], hDstau[row][0]);	
    hDtau[row][1]->Scale(bestFSig[0]); 
    hDstau[row][1]->Scale(bestFSig[1]);
  }
}


void SumHistos(TH1F *hSum, TH1F *h1, TH1F *h2, TH1F *h3, TH1F *h4, TH1F *h5, TH1F *h6, TH1F *h7){
  CopyHisto(hSum, h1);	
  if(h2) hSum->Add(h2); if(h3) hSum->Add(h3); if(h4) hSum->Add(h4);
  if(h5) hSum->Add(h5); if(h6) hSum->Add(h6); if(h7) hSum->Add(h7);
}

double CalcChi2(TH1F *hData, TH1F *hMC){
  double chi2 = 0; int ndof = -2;
  for(int bin=1; bin<=hData->GetNbinsX(); bin++) {
    double vBin = hData->GetBinContent(bin)-hMC->GetBinContent(bin);
    double eBin = sqrt(pow(hData->GetBinError(bin),2)+pow(hMC->GetBinError(bin),2));
    if(vBin != 0 && hData->GetBinError(bin) >= sqrt(8)){
      chi2 += pow(vBin/eBin,2);
      ndof++;
    }
  }
  //return TMath::Prob(chi2,ndof);
  return chi2;
}

void CopyHisto(TH1F *h1, TH1F *h2, int onlyError){
  for(int bin=1; bin<=h2->GetNbinsX(); bin++) {
    if(onlyError==0) h1->SetBinContent(bin, h2->GetBinContent(bin));
    h1->SetBinError(bin, h2->GetBinError(bin));
  }
}
