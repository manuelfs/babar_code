#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "babar_code/PlotsThesis/PlotUtils.cc"
#include "babar_code/prdPlots/TSVDUnfold.h"
#include "babar_code/FF/BToDtaunu.cc"
#include "babar_code/FF/BToDstaunu.cc"
#include <fstream>
#include <iostream>

#define nPads 6
using namespace TMath;
using namespace std;
using std::cout;
using std::endl;

void UnfoldQ2(int kreg = 0, int whichPlot = 0){
  TString histosName = "keys/eps/CSample/root/Q2_Folded.root";
  TFile *fHistos2=0; fHistos2 = new TFile(histosName);
  histosName = "keys/eps/CSample/root/Q2_Unfolded.root";
  TFile fHistos(histosName, "RECREATE");
  TString Folder = "AWG82/ntuples/small/FitRAll", treeName[] = {"Newx100", "Higx030", "Higx045"};
  //TString Folder = "AWG82/ntuples/small/FitRAll", treeName[] = {"Newx100", "Newx100", "Newx100"};
  TChain *MC[3];

  TCanvas can("can","",1400,800); can.Divide(3,2);TPad *cPad;
  TString hName;
  TString mcCuts[] = {"weight*((MCType==5&&candType==1||MCType==11&&candType==3)&&candM2>1.5)",
		      "weight*((MCType==6&&candType==2||MCType==12&&candType==4)&&candM2>1.5)"};
  TH1D *hbini[nPads], *hxini[nPads], *hxdat[nPads], *hbdat[nPads], *hTemp, *hQ2[nPads], *hEff[nPads];
  TH1D *hD[nPads], *hD2[nPads];
  TH2D *hAdet[nPads], *hBcov[nPads], *hXtau[nPads];

  BToDtaunu Dtaunu; BToDstaunu Dstaunu;
  int nRows = 2, nCols = 3, nBins = 18, nValues = 50;
  double limQ2[] = {3.5, 12.5}, tBmH[] = {0, 0.3, 0.45};
  //double limQ2[] = {3.5, 12.5}, tBmH[] = {0, 0,0};
  if(kreg==0){nBins=17; limQ2[0] = 4;}
  for(int col=0; col<nCols; col++){
    hName = Folder; hName += treeName[col]; hName += "_RunAll.root";
    MC[col] = new TChain("ntp1");
    MC[col]->Add(hName);
    for(int row=0; row<nRows; row++){
      int pad = nRows*col+row;
      cPad = (TPad *)can.cd(col+nCols*row+1);
      //cPad->SetLogy(1);

      hName = "xini"; hName += pad; 
      hxini[pad] = new TH1D(hName,"", nBins, limQ2[0], limQ2[1]);
      hxini[pad]->Sumw2(); 
      MC[col]->Project(hName,"trueQ2",mcCuts[row]);
      hName = "bini"; hName += pad; 
      hbini[pad] = new TH1D(hName,"", nBins, limQ2[0], limQ2[1]);
      hbini[pad]->Sumw2(); 
      MC[col]->Project(hName,"candQ2",mcCuts[row]);
      hName = "Adet"; hName += pad;
      hAdet[pad] = new TH2D(hName,"", nBins, limQ2[0], limQ2[1], nBins, limQ2[0], limQ2[1]);
      MC[col]->Project(hName,"trueQ2:candQ2",mcCuts[row]);
      hName = "Data_"; hName += pad; 
      hTemp = (TH1D *)fHistos2->Get(hName);
      if(hTemp) {
	if(kreg){
	  hName = "Data"; hName += pad;
	  hbdat[pad] = new TH1D(hName,"", nBins, limQ2[0], limQ2[1]);
	  hName = "Bcov"; hName += pad;
	  hBcov[pad] = new TH2D(hName,"", nBins, limQ2[0], limQ2[1], nBins, limQ2[0], limQ2[1]);
	  for(int binx=1; binx<=nBins; binx++){
	    double val = 0, err = 0;
	    if(binx>1){
	      val = hTemp->GetBinContent(binx-1);
	      err = hTemp->GetBinError(binx-1);
	    }
	    hbdat[pad]->SetBinContent(binx, val);
	    hbdat[pad]->SetBinError(binx, err);
	    for(int biny=1; biny<=nBins; biny++){
	      if(biny==binx) hBcov[pad]->SetBinContent(binx, biny, pow(err,2));
	      else hBcov[pad]->SetBinContent(binx, biny, 0);
	    }
	  }
	} else hbdat[pad] = (TH1D *)hTemp->Clone();
      } else continue;
      hName = "MC_"; hName += pad; hName += "0"; 
      hTemp = (TH1D *)fHistos2->Get(hName);
      double MCevents; 
      if(hTemp) MCevents = hTemp->Integral();
      else continue;
      hbini[pad]->Scale(MCevents/hbini[pad]->Integral());

      hName = "hQ2"; hName += pad; 
      hQ2[pad] = new TH1D(hName,"", nBins, limQ2[0], limQ2[1]);
      Dtaunu._gSR  = -Dtaunu.mb_quark*pow(tBmH[col],2);
      Dstaunu._gSR = -Dtaunu.mb_quark*pow(tBmH[col],2);      
      for(int bin=1; bin<=nBins; bin++){
	double q2 = hQ2[pad]->GetBinLowEdge(bin), val=0;
	double dq2 = (hQ2[pad]->GetBinLowEdge(bin)-q2)/(double)nValues;
	q2 += dq2/2.;
	for(int ival=0; ival<nValues; ival++) {
	  if(row==0) val += 1e16*Dtaunu.Compute(q2,1, Dtaunu.mTau);
	  else       val += 1e16*Dstaunu.Compute(q2,1, Dtaunu.mTau);
	  q2 += dq2;
	}
	hQ2[pad]->SetBinContent(bin,val); hQ2[pad]->SetBinError(bin,0);
      }
      hQ2[pad]->Scale(MCevents/hQ2[pad]->Integral());
      hEff[pad] = (TH1D*)hQ2[pad]->Clone();
      if(kreg) hEff[pad]->Divide(hxini[pad]);
      else hEff[pad]->Divide(hbini[pad]);

      if(kreg){
	TSVDUnfold unfold(hbdat[pad], hBcov[pad], hbini[pad], hxini[pad], hAdet[pad]);
	hxdat[pad] = unfold.Unfold(kreg-row);
	hXtau[pad] = unfold.GetXtau();
	for(int binx=1; binx<=nBins; binx++)
	  hxdat[pad]->SetBinError(binx, sqrt(hXtau[pad]->GetBinContent(binx, binx)));
	hD2[pad] = (TH1D*)unfold.GetD();
	if(col==1&&row==0){
	  fstream textFile; textFile.open("public_html/CovMatrix.txt",fstream::out);
	  for(int binx=1; binx<nBins; binx++){
	    for(int biny=1; biny<nBins; biny++)
	      textFile<<RoundNumber(hXtau[pad]->GetBinContent(binx, biny),0)<<"\t& ";
	    textFile<<" \\\\ "<<endl;
	  }
	  textFile<<endl<<endl;
	  for(int binx=1; binx<nBins; binx++){
	    for(int biny=1; biny<nBins; biny++)
	      textFile<<RoundNumber(hBcov[pad]->GetBinContent(binx, biny),0)<<"\t& ";
	    textFile<<" \\\\ "<<endl;
	  }
	}
      } else {
	hxdat[pad] = (TH1D *)hbdat[pad]->Clone();
      }
      double Devents = hbdat[pad]->Integral();
      double scaleData = Devents/hxdat[pad]->Integral();
      hxdat[pad]->Multiply(hEff[pad]); hxdat[pad]->Scale(scaleData); 
      hQ2[pad]->Scale(scaleData); 
      //hbdat[pad]->Multiply(hEff[pad]); hbdat[pad]->Scale(Devents/hbdat[pad]->Integral()); 
      double maxH = hxdat[pad]->GetMaximum();

      switch(whichPlot){
      case 0:
	hxdat[pad]->SetMaximum(maxH*1.3); hxdat[pad]->SetMinimum(0);
	hxdat[pad]->Draw("hist"); hbdat[pad]->Draw("e0 same");
	break;
      case 1:
	hQ2[pad]->Draw("hist"); 
	if(kreg) hxini[pad]->Draw("e0 same");
	else hbini[pad]->Draw("e0 same");
	break;
      case 2:
	hD[pad] = (TH1D*)hD2[pad]->Clone();
	hD[pad]->SetDirectory(0);
	hD[pad]->Draw();
	break;
      }

      fHistos.cd();
      hName = "UnfoldedData_"; hName += pad;
      hxdat[pad]->Write(hName);
      hName = "UnfoldedMC_"; hName += pad;
      hQ2[pad]->Write(hName);
      hName = "FoldedMC_"; hName += pad;
      hbini[pad]->Write(hName);
    }
  }
  hName = "public_html/UnfoldQ2_"; hName += kreg; hName += "_"; hName += whichPlot; hName += ".eps";
  can.SaveAs(hName);

  fHistos2->Close(); fHistos2->Delete();
  fHistos.Close(); 
}


