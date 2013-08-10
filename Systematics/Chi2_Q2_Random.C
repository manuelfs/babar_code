#include "TMath.h"
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TRandom3.h"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#define nPads 6
#define nMaxBins 20
using namespace TMath;
using namespace std;
using std::cout;
using std::endl;

double CalcChi2(TH1F *h1, TH1F *h2);
void CopyHisto(TH1F *h1, TH1F *h2, int onlyError=0);

void Chi2_Q2(int nReps=1, TString isRandom="BkgDssScaFit"){
  if(nReps>10000) {cout<<"Too many reps"<<endl; return;}
  TString Folder = "keys/eps/CSample/root/";
  TString fName[3] = {"Q2", "Q2_030", "Q2_045"}, hName;
  int nPDFs = 9, nRows = 2, nCols = 3, nBins = 17;
  TH1F *hComp[3][2][2], *hData[3][2], *hYields[nPads], *hMC[nPads][10], *hTemp, *hDssPDF[3][2];
  TH1F *hChi2[3][2], *hBB[3][2][2], *hDss[3][2][2], *hNorm[3][2], *hDtau[3][2][2], *hDstau[3][2][2];
  TH1F *hTot[3][2], *hQ2Yields[7][3][2][nMaxBins];
  double Err_BB[nMaxBins][2][2], limQ2[] = {4, 12.5}, NomPchi2[3][2], nDss[]={0,0};
  double YDss[2][2] = {{274.31, 87.967}, {34.994, 15.619}};
  vector<Double_t> Vchi2[3][2];
  TString histosName = "public_html/Test_Histograms"; histosName+=isRandom; histosName+=".root";
  TFile saveHistos(histosName, "RECREATE");
  
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
  TPad *cPad;
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
  for(int col=2; col<nCols; col++){ 
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
      nDss[row] = hDss[col][row][0]->Integral();

      // Normalization
      hNorm[col][row] = (TH1F*) hMC[pad][3]->Clone();
      hNorm[col][row]->Add(hMC[pad][5]); hNorm[col][row]->Add(hMC[pad][6]);

      // Signal
      hDstau[col][row][0] = (TH1F*) hMC[pad][7]->Clone();
      hDstau[col][row][1] = (TH1F*) hMC[pad][7]->Clone();
      hDtau[col][row][0] = (TH1F*) hMC[pad][8]->Clone();
      hDtau[col][row][1] = (TH1F*) hMC[pad][8]->Clone();

      // Chi2
      hComp[col][row][0] = (TH1F*) hMC[pad][8]->Clone();
      hComp[col][row][1] = (TH1F*) hMC[pad][8]->Clone();
      CopyHisto(hComp[col][row][1], hBB[col][row][0]);	
      hComp[col][row][1]->Add(hDss[col][row][1]);
      hComp[col][row][1]->Add(hNorm[col][row]);
      double fScale = 0.7;
      for(int bin=1; bin<=nBins; bin++) {
	int nBinsYields = 100;
	hName = "BBq2_"; hName += pad; hName += "_"; hName += bin;
	double val = hBB[col][row][0]->GetBinContent(bin);
	hQ2Yields[0][col][row][bin] = new TH1F(hName,"",nBinsYields, val*(1-fScale), val*(1+fScale));
	hName = "Dssq2_"; hName += pad; hName += "_"; hName += bin;
	val = hDss[col][row][0]->GetBinContent(bin);
	hQ2Yields[1][col][row][bin] = new TH1F(hName,"",nBinsYields, val*(1-fScale), val*(1+fScale));
	hName = "Dtauq2_"; hName += pad; hName += "_"; hName += bin;
	val = hDtau[col][row][0]->GetBinContent(bin);
	hQ2Yields[2][col][row][bin] = new TH1F(hName,"",nBinsYields, val*(1-fScale), val*(1+fScale));
	hName = "Dstauq2_"; hName += pad; hName += "_"; hName += bin;
	val = hDstau[col][row][0]->GetBinContent(bin);
	hQ2Yields[3][col][row][bin] = new TH1F(hName,"",nBinsYields, val*(1-fScale), val*(1+fScale));
	hName = "Dataq2_"; hName += pad; hName += "_"; hName += bin;
	val = hData[col][row]->GetBinContent(bin)-hComp[col][row][1]->GetBinContent(bin);
	if(row==0) val -= hDstau[col][row][0]->GetBinContent(bin);
	else  val -= hDtau[col][row][0]->GetBinContent(bin);
	hQ2Yields[4][col][row][bin] = new TH1F(hName,"",nBinsYields, val-50, val+50);
      }

      hComp[col][row][1]->Add(hDtau[col][row][1]);
      hComp[col][row][1]->Add(hDstau[col][row][1]);
      int ndof = 14; if(row) ndof = 12;
      NomPchi2[col][row] = TMath::Prob(CalcChi2(hData[col][row], hComp[col][row][1]), ndof);      
      hName = "hChi2"; hName += pad;
      //NomPchi2[col][row] = 0;
      //if(col==2 && row==0) hChi2[col][row] = new TH1F(hName, "",100, NomPchi2[col][row]*100, 1);
      //else hChi2[col][row] = new TH1F(hName, "",100, NomPchi2[col][row]*100, 100);
      if(col==2 && row==0) hChi2[col][row] = new TH1F(hName, "",100, 0, 1);
      else hChi2[col][row] = new TH1F(hName, "",100, 0, 100);
      hChi2[col][row]->SetDirectory(0);

    }

    // Repetitions
    int doUni = 1;
    double vGaus=0;
    TRandom3 rand(1); // 1 makes the first seed constant. 0 is a time dependent seed.
    //for(int rep=0; rep<nReps; rep++){

    int doSyst = 1, rep=0;
    double dx=0.25;
    for(double fBB=0; fBB<=1; fBB+=dx){
    for(double fDssTau=0; fDssTau<=1; fDssTau+=dx){
    for(double fDpipi=0; fDpipi<=1; fDpipi+=dx){
    for(double fDssY=0; fDssY<=1; fDssY+=dx){

      for(int row=0; row<nRows; row++){
	//////////////// BB variation ////////////////////////////////////////////////////////////////
	for(int bin=1; bin<=nBins; bin++) {
	  double errBin = Err_BB[bin-1][row][1];
	  double loError = Err_BB[bin-1][row][0]-errBin;
	  double hiError = Err_BB[bin-1][row][0]+errBin;
	  if(loError > 1-errBin) loError = 1-errBin;
	  if(hiError < 1+errBin) hiError = 1+errBin;
	  double mean = (hiError+loError)/2.;
	  double val = hBB[col][row][0]->GetBinContent(bin);
	  if(doUni==0) vGaus = rand.Gaus(mean, hiError-mean);
	  else         vGaus = rand.Uniform(loError, hiError);
	  if(vGaus<0) vGaus=0;
	  hBB[col][row][1]->SetBinContent(bin, val*vGaus);
	  if(doSyst) hBB[col][row][1]->SetBinContent(bin, val*(loError+fBB*(hiError-loError)));
	}
	if(isRandom.Contains("Bkg")) CopyHisto(hComp[col][row][0], hBB[col][row][1]);	
	else CopyHisto(hComp[col][row][0], hBB[col][row][0]);	

	////////// D**(l/tau)nu variation ////////////////////////////////////////////////////////////////
	double sigDss[] = {0.5, 0.3, 0.15};
	CopyHisto(hDss[col][row][1], hDss[col][row][0]);
	if(doUni==0) vGaus = rand.Gaus(0, sigDss[0]);
	else         vGaus = rand.Uniform(-sigDss[0], sigDss[0]);
	if(vGaus<-1) vGaus=-1;
	// Adding D**taunu scaled by G(0,sigma)
	if(doSyst) 
	  hDss[col][row][1]->Add(hDssPDF[1][row], 
				 (sigDss[0]*(-1+fDssTau*2))*nDss[row]/(YDss[0][row]/YDss[1][row]+1)); 
	else hDss[col][row][1]->Add(hDssPDF[1][row], vGaus*nDss[row]/(YDss[0][row]/YDss[1][row]+1)); 

	if(doUni==0) vGaus = rand.Gaus(0, sigDss[1]);
	else         vGaus = rand.Uniform(-sigDss[1], sigDss[1]);
	if(vGaus<0) vGaus=0;
	if(doSyst) 
	hDss[col][row][1]->Add(hDssPDF[2][row], (sigDss[1]*(-1+fDssTau*2))*nDss[row]); // Adding D**(->Dpipi)lnu 
	else hDss[col][row][1]->Add(hDssPDF[2][row], vGaus*nDss[row]); // Adding D**(->Dpipi)lnu 

	if(doUni==0) vGaus = rand.Gaus(1, sigDss[2]);
	else         vGaus = rand.Uniform(1-sigDss[2], 1+sigDss[2]);
	if(vGaus<0) vGaus=0;
	if(doSyst) 
	  hDss[col][row][1]->Scale((1+(sigDss[2]*(-1+fDssTau*2)))*nDss[row]/hDss[col][row][1]->Integral());
	else hDss[col][row][1]->Scale(vGaus*nDss[row]/hDss[col][row][1]->Integral());
	//cout<<vGaus<<", "<<hDss[col][row][0]->Integral()<<", "<<hDss[col][row][1]->Integral()<<endl;
	CopyHisto(hDss[col][row][1], hDss[col][row][0], 1); // Copying the uncertainties
	if(isRandom.Contains("Dss")) hComp[col][row][0]->Add(hDss[col][row][1]);
	else hComp[col][row][0]->Add(hDss[col][row][0]);
	hComp[col][row][0]->Add(hNorm[col][row]);
      }
      ////////// Signal fitting ////////////////////////////////////////////////////////////////
      double nDstau[2], nDtau[2], fDstau, fDtau;
      int row = 1;
      nDstau[row] = hTot[col][row]->Integral()-hComp[col][row][0]->Integral()-hDtau[col][row][0]->Integral();
      fDstau = nDstau[row]/hDstau[col][row][0]->Integral();
      if(fDstau<0.3 || fDstau>1.6)
	cout<<rep<<", "<<col<<" Tot "<<hTot[col][row]->Integral()<<", Comp "<<hComp[col][row][0]->Integral()
	    <<", Dtau fu "<<hDtau[col][row][0]->Integral()<<" -> nDstau "<<nDstau[row]<<", f "<<fDstau<<endl;
      row = 0;
      nDtau[row] = hTot[col][row]->Integral()-hComp[col][row][0]->Integral()-hDstau[col][row][0]->Integral();
      fDtau = nDtau[row]/hDtau[col][row][0]->Integral();
      if(fDtau<0.3 || fDtau>1.6) 
	cout<<rep<<", "<<col<<" Tot "<<hTot[col][row]->Integral()<<", Comp "<<hComp[col][row][0]->Integral()
	    <<", Dstau fd "<<hDstau[col][row][0]->Integral()<<" -> nDtau "<<nDtau[row]<<", f "<<fDtau<<endl;

      double df = 0.05, inif = 0.4, finf = 1.6, minChi2 = 1e10;
      double bestfTau = 0, bestfsTau = 0;
      //cout<<endl<<rep+1<<" of "<<nReps<<": col "<<col<<endl;
      for(double fTau = inif*fDtau; fTau<finf*fDtau; fTau += df*fDtau){
	for(double fsTau = inif*fDstau; fsTau<finf*fDstau; fsTau += df*fDstau){
	  double chi2=0;
	  for(int row=0; row<nRows; row++){
	    CopyHisto(hDtau[col][row][1],  hDtau[col][row][0]);	
	    CopyHisto(hDstau[col][row][1], hDstau[col][row][0]);	
	    hDtau[col][row][1]->Scale(fTau); 
	    hDstau[col][row][1]->Scale(fsTau);

	    CopyHisto(hComp[col][row][1], hComp[col][row][0]);
	    hComp[col][row][1]->Add(hDtau[col][row][1]);
	    hComp[col][row][1]->Add(hDstau[col][row][1]);
	    chi2 += CalcChi2(hData[col][row], hComp[col][row][1]);      
	  }
	  if(minChi2>chi2){
	    //cout<<minChi2<<" > "<<chi2<<" - fTau "<<fTau<<", fsTau "<<fsTau<<endl;
	    minChi2 = chi2;
	    bestfTau = fTau; bestfsTau = fsTau;
	  }
	}
      }
      if(!isRandom.Contains("Fit")) {bestfTau = fDtau; bestfsTau = fDstau;}
      if(bestfTau<=inif*fDtau || bestfTau>=finf*fDtau || bestfsTau<=inif*fDstau || bestfsTau>=finf*fDstau)
	cout<<"Rep "<<rep<<", col "<<col<<": bestfTau "<<bestfTau<<", bestfsTau "<<bestfsTau<<endl;
      for(int row=0; row<nRows; row++){
	CopyHisto(hDtau[col][row][1],  hDtau[col][row][0]);	
	CopyHisto(hDstau[col][row][1], hDstau[col][row][0]);	
	hDtau[col][row][1]->Scale(bestfTau); 
	hDstau[col][row][1]->Scale(bestfsTau);

	CopyHisto(hComp[col][row][1], hComp[col][row][0]);
	int doSca = 0;
	if(isRandom.Contains("Sca")) doSca = 1;
	for(int bin=1; bin<=nBins; bin++) {
	  hQ2Yields[0][col][row][bin]->Fill(hBB[col][row][1]->GetBinContent(bin));
	  hQ2Yields[1][col][row][bin]->Fill(hDss[col][row][1]->GetBinContent(bin));
	  hQ2Yields[2][col][row][bin]->Fill(hDtau[col][row][doSca]->GetBinContent(bin));
	  hQ2Yields[3][col][row][bin]->Fill(hDstau[col][row][doSca]->GetBinContent(bin));
	  double val = hData[col][row]->GetBinContent(bin)-hComp[col][row][1]->GetBinContent(bin);
	  if(row==0) val -= hDstau[col][row][doSca]->GetBinContent(bin);
	  else  val -= hDtau[col][row][doSca]->GetBinContent(bin);
	  hQ2Yields[4][col][row][bin]->Fill(val);
	}
	//cout<<"doSca "<<doSca<<", bestfTau "<<bestfTau<<", bestsfTau "<<bestfsTau<<endl;
	hComp[col][row][1]->Add(hDtau[col][row][doSca]);
	hComp[col][row][1]->Add(hDstau[col][row][doSca]);
	int ndof = 14; if(row) ndof = 12;
	double Pchi2 = TMath::Prob(CalcChi2(hData[col][row], hComp[col][row][1]), ndof);      
	hChi2[col][row]->Fill(Pchi2*100);
	Vchi2[col][row].push_back(Pchi2);

      }
      rep++;
    } // Repetitions
    }
    }
    }
    nReps=rep;
    for(int row=0; row<nRows; row++){
      int digits = 1; if(col==2&&row==0) digits=4;
      sort(Vchi2[col][row].begin(), Vchi2[col][row].end());
      for(int rep=0; rep<nReps; rep++){
	if(Vchi2[col][row][rep]>NomPchi2[col][row]){
	  int oneSigma = (int)(0.6827*(nReps-rep))+rep;
	  cout<<RoundNumber(Vchi2[col][row][oneSigma]*100,digits)<<": "<<rep<<", "<<oneSigma
	      <<" \t Best "<<RoundNumber(Vchi2[col][row][nReps-1]*100,digits)<<endl;
	  break;
	}
	if(rep==nReps-1) cout<<RoundNumber(Vchi2[col][row][rep]*100,digits)<<": "<<rep<<" - Final one "<<endl;
      }
      cPad = (TPad*)can.cd(col+nCols*row+1);
      //cPad->SetLogx(1);
      hChi2[col][row]->Draw(); 

      saveHistos.cd();
      for(int bin=1; bin<=nBins; bin++) 
	for(int his=0; his<5; his++)
	  hQ2Yields[his][col][row][bin]->Write();

    }
  } // Columns

  TString pName = "public_html/Test_Chi2_"; pName+= isRandom; pName += ".eps"; 
  can.SaveAs(pName);
  for(int row=0; row<nRows; row++){
    for(int dss=0; dss<3; dss++){
      //hDssPDF[dss][row]->Delete();
    }
    for(int col=0; col<nCols; col++){
      //hChi2[col][row]->Delete();
      //int pad = nRows*col+row;
      //hBB[col][row][0]->Delete();

      
    }
  }
  saveHistos.Close(); 
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
