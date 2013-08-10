#include "TChain.h"
#include "TString.h"
#include "TCut.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TPluginManager.h"
#include "TMVAGui.C"
#include <iostream>

#ifndef __CINT__
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#endif

void mvaDonut(TString Type = "Dl", int iChannel = 1, TString Sample = "Sig") {

  TString Channels[] = {"D0","Ds0","Dp","Dsp"};
  TString fname = "mva"; if(Sample=="Dss") fname += Sample; 
  fname += Type; fname += Channels[iChannel-1];
  TString outfileName = fname; outfileName += ".root";
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  int isDss = 0;
  if(Sample=="Dss") isDss=1;
  TMVA::Factory *factory = new TMVA::Factory( fname, outputFile, 
					      Form("!V:!Silent:%sColor", gROOT->IsBatch()?"!":"") );

  TChain c("ntp1");
  c.Add("~/releases/ntuplePID50/workdir/AWG82/ntuples/small/Add_R24MVA_RunAll.root");

  TString sigCuts[] = {"(MCType==1||MCType==3||MCType==5)", "(MCType==2||MCType==4||MCType==6)",
		       "(MCType==7||MCType==9||MCType==11)", "(MCType==8||MCType==10||MCType==12)",
		       "MCType>12"};
  TString bkgCuts[2][2] = {{"MCType>6", "(MCType>0&&MCType<7||MCType>12)"},
			   {"MCType>0&&MCType<13","MCType>0&&MCType<13"}};
  TString sigStr = "candLepTru==1&&"; if(isDss) sigStr += "pmisspi0"; else sigStr += "candPMiss";
  sigStr += ">0.2&&candType=="; sigStr += iChannel; sigStr += "&&";
  if(isDss) sigStr += sigCuts[4]; else sigStr += sigCuts[iChannel-1];
  TString bkgStr = "candType=="; bkgStr += iChannel; bkgStr += "&&";
  if(isDss) bkgStr += "pmisspi0"; else bkgStr += "candPMiss";
  bkgStr += ">0.2&&";
  if(Type=="Dl") bkgStr += bkgCuts[isDss][(iChannel-1)/2];
  else bkgStr += "MCType==0";
  if(isDss) {
    sigStr += "&&eextrapi0<1";
    bkgStr += "&&eextrapi0<1";
  } else {
    sigStr += "&&candEExtra<0.8";
    bkgStr += "&&candEExtra<0.8";
  }
  TCut sigCut = "1", bkgCut = "1", mycuts = "", mycutb = "";
  sigCut += sigStr; bkgCut += bkgStr;

//   --- Base ---
//   int nSig = 9, nDpi0 = 10;
//   TString sigVari[] = {"candEExtra","candMES","candDmass","candDeltam","candTagChargedMult","candBTagDeltam",
// 		       "candBTagDmass","candDeltaE","candCosT"};
//   TString Dpi0Vari[] = {"mpi0","candDmass","dmpi0","eextrapi0","ppi0","e1pi0","candCosT","candDeltam",
// 			"candMES","candDeltaE"};

//   --- NoDmNoMp0 ---
//   int nSig = 8, nDpi0 = 9;
//   TString sigVari[] = {"candEExtra","candMES","candDmass","candDeltam","candTagChargedMult",
// 		       "candBTagDmass","candDeltaE","candCosT"};
//   TString Dpi0Vari[] = {"candDmass","dmpi0","eextrapi0","ppi0","e1pi0","candCosT","candDeltam",
// 			"candMES","candDeltaE"};
//   sigCuts[4] = "MCType>12&&mpi0>.125&&mpi0<.145";

//  ---  NoMes ---
//   int nSig = 8, nDpi0 = 9;
//   TString sigVari[] = {"candEExtra","candDmass","candDeltam","candTagChargedMult","candBTagDeltam",
// 		       "candBTagDmass","candDeltaE","candCosT"};
//   TString Dpi0Vari[] = {"mpi0","candDmass","dmpi0","eextrapi0","ppi0","e1pi0","candCosT","candDeltam",
// 			"candDeltaE"};

//   --- NoMulYesDm ---
//   int nSig = 8, nDpi0 = 11;
//   TString sigVari[] = {"candEExtra","candMES","candDmass","candDeltam","candBTagDeltam",
// 		       "candBTagDmass","candDeltaE","candCosT"};
//   TString Dpi0Vari[] = {"mpi0","candDmass","dmpi0","eextrapi0","ppi0","e1pi0","candCosT","candDeltam",
// 			"candMES","candDeltaE","candBTagDeltam"};

//   --- Final ---
//   int nSig = 7, nDpi0 = 8;
//   TString sigVari[] = {"candEExtra","candDmass","candDeltam","candTagChargedMult",
// 		       "candBTagDmass","candDeltaE","candCosT"};
//   TString Dpi0Vari[] = {"candDmass","dmpi0","eextrapi0","ppi0","e1pi0","candCosT","candDeltam",
// 			"candDeltaE"};
//   sigCuts[4] = "MCType>12&&mpi0>.12&&mpi0<.15";

//   --- Fin ---
//   int nSig = 7, nDpi0 = 8;
//   TString sigVari[] = {"candEExtra","candMES","candDmass","candDeltam",
// 		       "candBTagDmass","candDeltaE","candCosT"};
//   TString Dpi0Vari[] = {"candDmass","dmpi0","eextrapi0","ppi0","e1pi0","candCosT","candDeltam",
// 			"candDeltaE"};
//   sigCuts[4] = "MCType>12&&mpi0>.12&&mpi0<.15";

//   --- NoMesMul ---
//   int nSig = 7, nDpi0 = 8;
//   TString sigVari[] = {"candEExtra","candDmass","candDeltam","candBTagDeltam",
// 		       "candBTagDmass","candDeltaE","candCosT"};
//   TString Dpi0Vari[] = {"candDmass","dmpi0","eextrapi0","ppi0","e1pi0","candCosT","candDeltam",
// 			"candDeltaE"};
//   sigCuts[4] = "MCType>12&&mpi0>.12&&mpi0<.15";

//   --- No3 ---
//   int nSig = 6, nDpi0 = 8;
//   TString sigVari[] = {"candEExtra","candDmass","candDeltam",
// 		       "candBTagDmass","candDeltaE","candCosT"};
//   TString Dpi0Vari[] = {"candDmass","dmpi0","eextrapi0","ppi0","e1pi0","candCosT","candDeltam",
// 			"candDeltaE"};
//   sigCuts[4] = "MCType>12&&mpi0>.12&&mpi0<.15";

//   --- End ---
  int nSig = 8, nDpi0 = 8;
  TString sigVari[] = {"candEExtra","candDmass","candDeltam","candTagChargedMult","candBTagDeltam",
		       "candBTagDmass","candDeltaE","candCosT"};
  TString Dpi0Vari[] = {"candDmass","dmpi0","eextrapi0","ppi0","e1pi0","candCosT","candDeltam",
			"candDeltaE"};
  sigCuts[4] = "MCType>12&&mpi0>.12&&mpi0<.15";

//   --- EndNoDm ---
//   int nSig = 7, nDpi0 = 8;
//   TString sigVari[] = {"candEExtra","candDmass","candTagChargedMult","candBTagDeltam",
// 		       "candBTagDmass","candDeltaE","candCosT"};
//   TString Dpi0Vari[] = {"candDmass","dmpi0","eextrapi0","ppi0","e1pi0","candCosT","candDeltam",
// 			"candDeltaE"};
//   sigCuts[4] = "MCType>12&&mpi0>.12&&mpi0<.15";

  factory->SetInputTrees(&c, sigCut, bkgCut);
  if(isDss==0){
    for(int vari = 0; vari < nSig; vari++){
      if(sigVari[vari]=="candDeltam" && iChannel%2==1) continue;
      char variChar = 'F';
      if(sigVari[vari]=="candTagChargedMult") variChar = 'I';
      factory->AddVariable(sigVari[vari], variChar);
    }
  } else {
    for(int vari = 0; vari < nDpi0; vari++){
      if(Dpi0Vari[vari]=="candDeltam" && iChannel%2==1) continue;
      factory->AddVariable(Dpi0Vari[vari], 'F');
    }
  }
  factory->SetWeightExpression("wFF");
  factory->PrepareTrainingAndTestTree( mycuts, mycutb,
				       "NSigTest=100:NBkgTest=100:SplitMode=Random:NormMode=NumEvents:!V" );

  factory->BookMethod( TMVA::Types::kBDT, "BDT", 
		       "!H:!V:NTrees=500:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrength=2.5" );

  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();    
   
  // Save the output
  outputFile->Close();
  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  delete factory;

  // Launch the GUI for the root macros
  //if (!gROOT->IsBatch()) TMVAGui( outfileName );
  gROOT->ProcessLine(".q");
}
