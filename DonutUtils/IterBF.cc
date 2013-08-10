//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: IterBF.cc,v 1.2 2012/08/23 02:22:17 manuelf Exp $
//
// Description:
//      IterBF - Reads the output of a fit, updates BFIterations.txt with
//      the fitted values of RD, RDs and the D^{*}lnu BF, re-generates
//      FitRAll_RunAll.root, and with that the fitsamples 25-32
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      11/04/22 manuelf -- Created
//------------------------------------------------------------------------

#include "TCut.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "DonutUtils/cuts.cc"
#include "DonutUtils/KeysUtils.cc"
#include <fstream>
#include <iostream>

using std::cout;
using std::endl;

TH1F *_hGlobal;
Double_t HistoConv(Double_t *x, Double_t *par);
void ResolutionFitRAll(TString txtFile, double minX, double maxX, TString FileTag);
void EventMixRAll(TString FileTag);
void GenPdfRAll();
void fillTree(TTree *inputTree, TTree *outputTree, Int_t thecut, Int_t dss, TString sepTru_Print, double xData);

int main(int argc, char *argv[]){
  if (argc < 2 || argc > 5) {
    cout << "USAGE: IterBF txtFile [FileTag=RAllNewx100_RunAll] [minX=-0.3] [maxX=0.4]" << endl;
    return 1;
  }
  TString txtFile = argv[1];
  TString FileTag = "RAllNewx100_RunAll";
  if(argc>2) {FileTag = argv[2];} 
  double minX = -0.3;
  if(argc>3) {TString temp_s = argv[3]; minX = temp_s.Atof();} 
  double maxX = 0.4;
  if(argc>4) {TString temp_s = argv[4]; maxX = temp_s.Atof();} 

//    cout<<endl<<"Starting iteration: Fitting the unsimulated Mmiss resolution\n"<<
//      "========================================================================"<<endl;
//    ResolutionFitRAll(txtFile, minX, maxX, FileTag);
//   return 1;

  double Yields[16][9];
  TString buffer;
  fstream textFile; textFile.open(txtFile,fstream::in);
  int begChan = 0, endChan = 4, IndText[8] = {0,4,1,5,2,6,3,7};
  for(int i=2*begChan; i<2*endChan; i++){
    int index = IndText[i];
    for(int j=0; j<9; j++){
      textFile>>buffer>>Yields[index+8][j]>>Yields[index][j]>>buffer>>buffer>>buffer;
      if(index>3 && j==6){
	textFile>>buffer;
	break;
      }
    }
  }
  TString names[8] = {"D0", "Ds0", "Dp", "Dsp", "RD", "RDs","Ds0=>D0","Dsp=>Dp"}; 
  fstream iterFile; iterFile.open("../DonutUtils/BFIteration.txt",fstream::out);
  for(int i=0; i<4; i++) iterFile<<names[i]<<"\t"<<RoundNumber(Yields[i][1],5,Yields[i+8][1])<<endl;
  for(int i=0; i<2; i++) iterFile<<names[i+4]<<"\t"<<
			   RoundNumber(Yields[i][0]*Yields[i+8][1],5,Yields[i+8][0]*Yields[i][1])<<endl;
  for(int i=0; i<2; i++) iterFile<<names[i+6]<<"\t"<<
			   RoundNumber(Yields[2*i][3]*Yields[2*i+1+8][1],5,Yields[2*i+1][1]*Yields[2*i+8][3])<<endl;
			   
  iterFile.close();


  cout<<endl<<"Re-generating ntuple with the new D(*)lnu and D(*)taunu BF weights\n"<<
    "==========================================================================="<<endl;
  EventMixRAll(FileTag);
//   GenPdfRAll();
  return 1;
}

void EventMixRAll(TString FileTag){
  float candM2,candPstarLep, weight=-1, BFIterRatio[8], MCRD, MCRDs;
  int MCType, candIsMu, candType, MCTaumode, MCD;

  //TString RootFile = "AWG82/ntuples/small/Fit"; RootFile += FileTag; RootFile += "_RunAll.root";
  TString RootFile = "AWG82/ntuples/small/Fit"; RootFile += FileTag; RootFile += ".root";
  TChain inputTree("ntp1"); inputTree.Add(RootFile);
  inputTree.SetBranchAddress("weight",&weight);
  inputTree.SetBranchAddress("candM2",&candM2);
  inputTree.SetBranchAddress("candPstarLep",&candPstarLep);
  inputTree.SetBranchAddress("candType",&candType);
  inputTree.SetBranchAddress("MCType",&MCType);
  inputTree.SetBranchAddress("MCTaumode",&MCTaumode);
  inputTree.SetBranchAddress("MCD",&MCD);
  inputTree.SetBranchAddress("candIsMu",&candIsMu);
  inputTree.SetBranchAddress("MCRD",&MCRD);
  inputTree.SetBranchAddress("MCRDs",&MCRDs);

  TTree *outputTree = inputTree.CloneTree(0);
  outputTree -> Branch("weight",&weight,"weight/F");
  outputTree -> Branch("candType",&candType,"candType/I");
  outputTree -> Branch("MCType",&MCType,"MCType/I");
  outputTree -> Branch("MCD",&MCD,"MCD/I");
  outputTree -> Branch("MCTaumode",&MCTaumode,"MCTaumode/I");
  outputTree -> Branch("candIsMu",&candIsMu,"candIsMu/I");
  outputTree -> Branch("candM2",&candM2,"candM2/F");
  outputTree -> Branch("candPstarLep",&candPstarLep,"candPstarLep/F");
  outputTree -> Branch("MCRD",&MCRD,"MCRD/F");
  outputTree -> Branch("MCRDs",&MCRDs,"MCRDs/F");


  fstream BFIter; BFIter.open("../DonutUtils/BFIteration.txt",fstream::in);
  TString dummy; 
  for(int i=0; i<8; i++) BFIter>>dummy>>BFIterRatio[i];
  long iEntries = inputTree.GetEntries();
  for(int entry=0; entry<iEntries; entry++){
    inputTree.GetEvent(entry);
    if (MCType == 1 || MCType == 3 || MCType == 5) weight *= BFIterRatio[0];
    if (MCType == 2 || MCType == 4 || MCType == 6) weight *= BFIterRatio[1];
    if (MCType == 7 || MCType == 9 || MCType == 11) weight *= BFIterRatio[2];
    if (MCType == 8 || MCType == 10 || MCType == 12) weight *= BFIterRatio[3];
    if (MCType == 5)   weight *= BFIterRatio[4];
    if (MCType == 6)   weight *= BFIterRatio[5];
    if (MCType == 11)  weight *= BFIterRatio[4];
    if (MCType == 12)  weight *= BFIterRatio[5];
    if ((MCType == 2 || MCType == 4 || MCType == 6)&&candType==1) weight *= BFIterRatio[6];
    if ((MCType == 8 || MCType == 10 || MCType == 12)&&candType==3) weight *= BFIterRatio[7];
    MCRD = 0.3 * BFIterRatio[4]; MCRDs = 0.25 * BFIterRatio[5];
    outputTree -> Fill();
  }
  TString nameFolder = RootFile;
  nameFolder.ReplaceAll("Fit","FitIter");
  TFile f(nameFolder,"RECREATE");  f.cd();
  outputTree->Write();
  f.Write();  f.Close();
  cout<<"Written "<<nameFolder<<endl;
}

// void EventMixRAll(TString FileTag){
//   TString RootFile = "AWG82/ntuples/small/RAll_RunAll.root";
//   TString scale2Data = "yes_YES";
//   TString weightName = "babar_code/Reweight/wTotal.txt";

//   float candM2,candPstarLep, mm2pi0, weight=-1;
//   int MCType, candIsMu, candType, MCTaumode, MCD;
//   TTree *outputTree = new TTree("ntp1","cands");
//   outputTree -> Branch("weight",&weight,"weight/F");
//   outputTree -> Branch("candM2",&candM2,"candM2/F");
//   outputTree -> Branch("candPstarLep",&candPstarLep,"candPstarLep/F");
//   outputTree -> Branch("candType",&candType,"candType/I");
//   outputTree -> Branch("MCType",&MCType,"MCType/I");
//   outputTree -> Branch("MCD",&MCD,"MCD/I");
//   outputTree -> Branch("MCTaumode",&MCTaumode,"MCTaumode/I");
//   outputTree -> Branch("candIsMu",&candIsMu,"candIsMu/I");

//   double totMCB = 0, totuds = 0, totccbar = 0, totdata = 0, totOffdata = 0;
//   double totDpipi = 6874000+6799000, totDss = 7049000+7036000;
//   getNumberB(RootFile, "All", totMCB, totdata, totuds, totccbar, totOffdata);
//   double wuds = totMCB/totuds*2.09/1.05;     
//   double wMC = totdata/totMCB; if(!scale2Data.Contains("yes")) wMC = 1;
//   if(RootFile.Contains("_4_Run") || RootFile.Contains("_14_Run")) wuds /= 2.;
//   TString NameTrees[4] = {RootFile, "AWG82/ntuples/small/uds_RunAll.root", "AWG82/ntuples/small/ccbar_RunAll.root",
// 			  "AWG82/ntuples/small/GenDss_RunAll.root"};
//   TCut Cuts[2] = {MvaAll, dssMvaAll};
//   for(int cut=0; cut<2; cut++){
//     Cuts[cut] += MEScut;
//     for(int tree=0; tree<3; tree++){
//       if((RootFile.Contains("Data") ||scale2Data.Contains("NO")) && tree>0) continue;
//       double dEntries;
//       int iEntries, isCocktail=0;
//       if(tree>0) isCocktail=-1;
//       TTree *inputTree = WeightedTree(NameTrees[tree], dEntries, weightName,isCocktail,Cuts[cut]);
//       inputTree -> SetBranchAddress("weight",&weight);
//       inputTree -> SetBranchAddress("mm2pi0",&mm2pi0);
//       inputTree -> SetBranchAddress("candM2",&candM2);
//       inputTree -> SetBranchAddress("candPstarLep",&candPstarLep);
//       inputTree -> SetBranchAddress("candType",&candType);
//       inputTree -> SetBranchAddress("MCType",&MCType);
//       inputTree -> SetBranchAddress("MCTaumode",&MCTaumode);
//       inputTree -> SetBranchAddress("MCD",&MCD);
//       inputTree -> SetBranchAddress("candIsMu",&candIsMu);
//       iEntries = inputTree->GetEntries();
//       for(int entry=0; entry<iEntries; entry++){
// 	inputTree -> GetEvent(entry);
// 	if(tree<3) weight *= wMC;
// 	else {
// 	  if(MCType!=14) continue;
// 	  if(scale2Data.Contains("yes")) {
// 	    if(MCType==14) weight *= 0.01/(totDss/totdata);
// 	    else if(MCType==15) weight *= 0.01/(totDpipi*0.25/totdata);
// 	    else if(MCType==16) weight *= 0.01/(totDpipi*0.75/totdata);
// 	    else cout<<"Error. In GenDss we have MCType "<<MCType<<endl;
// 	  }
// 	  MCType = 16;
// 	}
// 	if(cut==1) {
// 	  candM2 = mm2pi0;
// 	  candType += 4;
// 	}
// 	if(tree==1 || tree==2) weight *= wuds;
// 	if(RootFile.Contains("Data")) weight = 1.;
// 	outputTree -> Fill();
//       }
//       inputTree->Delete();
//     }
//   }

//   TString nameFolder = RootFile;
//   RootFile.Remove(0,RootFile.Last('/')+1);
//   nameFolder.Remove(nameFolder.Last('/')+1,nameFolder.Length());
//   nameFolder += "FitIter"; nameFolder += RootFile;
//   nameFolder.ReplaceAll("RAll",FileTag);
//   TFile f(nameFolder,"RECREATE");
//   f.cd();
//   outputTree->Write();
//   f.Write();
//   f.Close();
//   cout<<"Written "<<nameFolder<<endl;
// }

void ResolutionFitRAll(TString txtFile, double minX, double maxX, TString FileTag){
  Styles style; style.setPadsStyle(4); style.applyStyle();
  TCanvas can; can.Divide(2,2);
  int isDss = 0, nbins = 35;//, nbinsConv = 200;  // Enough bins to plot it smoothly (option "c")
  if(isDss){ minX = -1.5; maxX = 2.5;}
  TString cuts = "(candM2>"; cuts+=minX; cuts+="&&candM2<"; cuts+=maxX; cuts+="&&candType==";
  TString fName = "HistoConvo", gName = "hGlobal", dName = "hData";
  TString DataName = "AWG82/ntuples/small/Fit"; DataName += FileTag; DataName += ".root";
  DataName.ReplaceAll("RAll","Data");
  TChain data("ntp1");
  data.Add(DataName);
  //data.Add("AWG82/ntuples/small/FitR24Rests_1_2_RunAll.root");
  double Yields[16][9];
  TString buffer; 
  fstream textFile; textFile.open(txtFile,fstream::in);
  int begChan = 0, endChan = 4, IndText[8] = {0,4,1,5,2,6,3,7};
  for(int i=2*begChan; i<2*endChan; i++){
    int index = IndText[i];
    for(int j=0; j<9; j++){
      textFile>>buffer>>Yields[index+8][j]>>Yields[index][j]>>buffer>>buffer>>buffer;
      if(index>3 && j==6){
	textFile>>buffer;
	break;
      }
    }
  }
  TH2F *pdf[70];
  for (int i = 1 ; i <= 68 ; i ++) {
    if(i==49)i=51; if(i==59)i=61; 
    TString hname = "keys/root/fitNoConv/pdfKeys_"; hname += i; hname += "_Fit.root";
    TFile hfile(hname); 
    TString pdfName = "pdf"; pdfName += i;
    pdf[i] = (TH2F *)(hfile.Get("h2"))->Clone(pdfName);    
    pdf[i]->SetDirectory(0);
  }
  int nm2bin = pdf[1]->GetNbinsX(); 
  double minx=-4, maxx=12, miny=0, maxy=2.4, valResW[4][2];
  int Indpdf[8][9] = {{61,51,41,17,37,13,9,5,1},   {62,52,42,18,38,14,10,6,2},
		      {63,53,43,19,39,15,11,7,3},  {64,54,44,20,40,16,12,8,4},
		      {65,55,45,29,25,21,33,-1,-1},{66,56,46,30,26,22,34,-1,-1},
		      {67,57,47,31,27,23,35,-1,-1},{68,58,48,32,28,24,36,-1,-1}};

  int IndYield[3][9] = {{7,6,5,4,8,3,1,2,0},{7,6,5,4,8,3,1,2,0},{5,4,3,2,1,0,6,-1,-1}};
  TF1 *hConvolution = new TF1(fName,HistoConv,minX,maxX,1);;
  hConvolution->SetParameter(0,0.025); 
  hConvolution->SetParLimits(0,0,0.1);
  TH1F *hData[4], *hKeys[4], pdfProj[4][9];
  for(int pad=0; pad<4; pad++){
    can.cd(pad+1);
    int channel = pad+4*isDss;
    TString totcuts = cuts; totcuts += channel+1; totcuts += ")*weight";
    int type = 0;
    if(!isDss) {if(pad%2==1) type=1;}
    else type=2;
    for(int i=0; i<9; i++) {
      if(Indpdf[channel][i]<0) break;
      TString pdfName = "pdfProj"; pdfName += pad; pdfName += i;
      double scale = Yields[channel][IndYield[type][i]]/pdf[Indpdf[channel][i]]->Integral();
      pdfProj[pad][i] = binHisto(pdf[Indpdf[channel][i]],nm2bin,minx,maxx,miny,maxy,pdfName,"candM2");
      pdfProj[pad][i].Scale(scale);
      if(i>0) pdfProj[pad][i].Add(&pdfProj[pad][i-1]);
    }
    _hGlobal = &pdfProj[pad][8-2*isDss];
    TString dName = "hData"; dName += pad;
    hData[pad] = new TH1F(dName,"",nbins,minX,maxX);
    hData[pad]->Sumw2();
    TString vari = "candM2>>"; vari += dName;
    data.Draw(vari,totcuts);
    data.UnbinnedFit(fName,"candM2",totcuts,"M E Q");
    valResW[pad][0] = hConvolution->GetParameter(0); valResW[pad][1] = hConvolution->GetParError(0);
    dName += "Keys";
    hKeys[pad] = (TH1F*)((TH1F*)hConvolution->GetHistogram())->Clone(dName);
    double scale = hData[pad]->Integral()*(double)hKeys[pad]->GetNbinsX()/hKeys[pad]->Integral()/(double)nbins;
    hKeys[pad]->Scale(scale);
    hKeys[pad]->Draw("same c");
  }

  TString names[] = {"D0 ", "Ds0", "Dp ", "Dsp"}; 
  double ve=0, e=0, ee=0;
  for(int chan=0; chan<4; chan++){
    cout<<names[chan]<<": "<<RoundNumber(valResW[chan][0]*1000,1)<<" +- "<<RoundNumber(valResW[chan][1]*1000,1)<<endl;
    if(valResW[chan][1]){
      ve += valResW[chan][0]/valResW[chan][1];
      e  += 1/valResW[chan][1];
      ee += 1/valResW[chan][1]/valResW[chan][1];
    } else cout<<"Error 0"<<endl;
  }
  cout<<"Average Unsimulated Mmiss resolution:\t"<<RoundNumber(1000*ve,1,e)<<" +- "
      <<RoundNumber(1000*1/sqrt(ee),1)<<endl<<endl;
  //fstream ResolFile; ResolFile.open("../DonutUtils/ResolutionWidths.txt",fstream::out);
  //for(int chan=0; chan<4; chan++) ResolFile << RoundNumber(ve,4,e)<<endl;
  TString epsName = "keys/eps/Resolution/IterConv.eps";
  can.SaveAs(epsName);
}

//par[1] the normalization
Double_t HistoConv(Double_t *x, Double_t *par){
  int nbins = _hGlobal->GetNbinsX();
  double minX = _hGlobal->GetXaxis()->GetBinLowEdge(_hGlobal->GetXaxis()->GetFirst());
  double maxX = _hGlobal->GetXaxis()->GetBinLowEdge(_hGlobal->GetXaxis()->GetLast()+1);
  double nSigma = 4., wBin = (maxX-minX)/(double)nbins;
  double width = par[0];

  double minG = x[0]-(double)nSigma*width, maxG = x[0]+(double)nSigma*width;
  if(minG<minX) {minG=minX;} if(maxG>maxX) {maxG=maxX;} 
  double AreaG = IntG(x[0],width,minG,maxG), valH = 0;
  int iniBin = (int)((minG-minX)/wBin)+1, finBin = (int)((maxG-minX)/wBin)+1;
  for(int binH=iniBin; binH<=finBin; binH++){
    double minH = (double)(binH-1)*wBin+minX, maxH = (double)(binH)*wBin+minX;
    if(minH<minG) minH=minG; if(maxH>maxG) maxH=maxG; 
    valH += IntG(x[0], width, minH, maxH)*_hGlobal->GetBinContent(binH);
  }

  return valH/AreaG;
}

void GenPdfRAll(){
  TString ntupleTag = "IterRAll";
  TString Print_SCALE = "yes_YES";

  TString ntupleName = "AWG82/ntuples/small/Fit"; ntupleName += ntupleTag; ntupleName += "_RunAll.root";
  TChain *chain = new TChain("ntp1");
  chain->Add(ntupleName);
  cout<<endl<<"Storing the new PDF samples from "<<ntupleName<<" with "<<chain->GetEntries()
      <<" entries\n====================================================================================================================="<<endl<<endl;
  double totMCB = 0, totuds = 0, totccbar = 0, totdata = 0, totOffdata = 0;
  getNumberB(ntupleName, "All", totMCB, totdata, totuds, totccbar, totOffdata);
  double xData=totMCB/totdata;

  TFile *f;
  TString pdfFolder = "fitSamples/"; pdfFolder += ntupleTag; 
  gSystem->mkdir(pdfFolder);
  for (int i = 25 ; i <= 32 ; i ++) {
    int isDss = 1;
    TString pdfName = pdfFolder; pdfName += "/pdfSample";pdfName += i; pdfName += ".root";
    f = new TFile(pdfName,"RECREATE");
    TTree *t = new TTree("ntp1","cands");
    fillTree(chain,t,i,isDss, Print_SCALE, xData);
    t->Write();
    t->ResetBranchAddresses();
    f->Write();
    f->Close();
  } 
}

void fillTree(TTree *inputTree, TTree *outputTree, Int_t thecut, Int_t isDss, TString Print_SCALE, double xData){
  bool doPrint_SCALE = true;
  if(Print_SCALE.Contains("no")) doPrint_SCALE = false;
  Float_t myCandM2;
  Int_t MCType,candType;
  Float_t candM2,candPstarLep,weight;

  inputTree -> SetBranchAddress("candPstarLep",  &candPstarLep);            
  inputTree -> SetBranchAddress("candM2",        &candM2);                  
  inputTree -> SetBranchAddress("candType",      &candType);                
  inputTree -> SetBranchAddress("MCType",        &MCType);                  
  inputTree -> SetBranchAddress("weight",        &weight);                       

  outputTree -> Branch("MCType",&MCType,"MCType/I");
  outputTree -> Branch("candPstarLep",&candPstarLep,"candPstarLep/F");
  outputTree -> Branch("candM2",&myCandM2,"candM2/F");
  outputTree -> Branch("candType",&candType,"candType/I");
  outputTree -> Branch("weight",&weight,"weight/F");

  bool Sep12 = false;
  Int_t n1 = 0, n2 = 0, nevt = (Int_t)inputTree -> GetEntries();
  for (int j = 0 ; j < nevt ; j ++) {
    inputTree -> GetEvent(j);

    if(!isSample(thecut, MCType, candType)) continue;
    if(Print_SCALE.Contains("YES")) weight *= xData;
    myCandM2 = candM2;
    outputTree -> Fill();
    if(!Sep12) n1++;
    else {
      if(thecut>=17&&thecut<25){
	if(MCType==13) n1++;
	else n2++;
      }
      if(thecut>=24&&thecut<29){
	if(MCType%2==1) n1++;
	else n2++;
      }
    }
  }
}

