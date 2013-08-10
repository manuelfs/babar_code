#ifndef CUTS_HH
#define CUTS_HH

#include "TCut.h"
#include "TString.h"

//float BDTCuts[4] = {0.63, 0.48, 0.51, 0.41};
float BDTCuts[4] = {0.458, 0.48, 0.36, 0.41};               // onDx200
//float BDTCuts[4] = {0.7207, 565, 6, 525};               //0.5x
//float BDTCuts[4] = {0.6853, 0.5358, 0.57305, 0.487};    //0.67x
//float BDTCuts[4] = {0.6616, 0.514, 0.5468, 0.451};      //0.8x
//float BDTCuts[4] = {0.5938, 0.4546, 0.474, 0.3659};     //1.2x
//float BDTCuts[4] = {0.5413, 0.4133, 0.432, 0.3171};     //1.5x
//float BDTCuts[4] = {0.458, 0.344, 0.36, 0.265};         //2.0x
//float BDTCuts[4] = {0.381, 0.2804, 0.308, 0.224};       //2.5x
//float BDTCuts[4] = {0.3065, 0.235, 0.253, 0.194};             //3.0x

//float BDTCuts[4] = {0.3, 0.3, 0.3, 0.3};
//float BDTCuts[4] = {0.43, 0.28, 0.31, 0.21};
//float BDTCuts[4] = {0.48, 0.33, 0.36, 0.26};
//float BDTCuts[4] = {0.53, 0.38, 0.41, 0.31};
//float BDTCuts[4] = {0.51, 0.28, 0.36, 0.21};
float BDTCutspi0[2][4] = {{-0.45, -0.45, -0.45, -0.45},{-0.35, -0.35, -0.3, -0.35}};
//float BDTCutspi0[2][4] = {{-0.5, -0.5, -0.5, -0.5},{-0.4, -0.4, -0.35, -0.4}};
//float BDTCutspi0[2][4] = {{-0.55, -0.55, -0.55, -0.55},{-0.45, -0.45, -0.10, -0.45}};
// TString MvaCut = "";
// for(int chan=1; chan<=4; chan++){
//   MvaCut += "(candType=="; MvaCut += chan;
//   MvaCut += "&&candMvaDl>"; MvaCut += BDTCuts[chan-1];
//   if(chan<4) MvaCut += ")||";
//   else MvaCut += ")";
// }

TCut M2Signal = "candM2>1.5";
TCut MEScut = "candMES>5.27";
TCut M2P = "candPstarLep>0&&candPstarLep<2.4&&candM2>-4&&candM2<12";
TCut PMiss = "candPMiss>.2";
TCut Q2 = "candQ2>4";
TCut lowQ2 = "candQ2<=4";
TCut ee1 = "candType<3&&candEExtra<.2";
TCut ee2 = "candType==3&&candEExtra<.15";
TCut ee3 = "candType==4&&candEExtra<.3";
TCut ee = (ee1||ee2||ee3);
//    TCut MvaBp = "(candMvaDl>0.63&&candType==1)||(candMvaDl>0.48&&candType==2)";
//    TCut MvaB0 = "(candMvaDl>0.51&&candType==3)||(candMvaDl>0.41&&candType==4)";
   TCut MvaBp = "(candMvaDl>0.458&&candType==1)||(candMvaDl>0.48&&candType==2)";      // onDx200
   TCut MvaB0 = "(candMvaDl>0.360&&candType==3)||(candMvaDl>0.41&&candType==4)";      // onDx200
//TCut MvaBp = "(candMvaDl>0.7207&&candType==1)||(candMvaDl>0.565&&candType==2)";      //0.5x
//TCut MvaB0 = "(candMvaDl>0.6&&candType==3)||(candMvaDl>0.525&&candType==4)";         //0.5x
//TCut MvaBp = "(candMvaDl>0.6853&&candType==1)||(candMvaDl>0.5358&&candType==2)";     //0.67x
//TCut MvaB0 = "(candMvaDl>0.57305&&candType==3)||(candMvaDl>0.487&&candType==4)";     //0.67x
//TCut MvaBp = "(candMvaDl>0.6616&&candType==1)||(candMvaDl>0.514&&candType==2)";      //0.8x
//TCut MvaB0 = "(candMvaDl>0.5468&&candType==3)||(candMvaDl>0.451&&candType==4)";      //0.8x
//TCut MvaBp = "(candMvaDl>0.5938&&candType==1)||(candMvaDl>0.4546&&candType==2)";     //1.2x
//TCut MvaB0 = "(candMvaDl>0.474&&candType==3)||(candMvaDl>0.3659&&candType==4)";      //1.2x
//TCut MvaBp = "(candMvaDl>0.5413&&candType==1)||(candMvaDl>0.4133&&candType==2)";     //1.5x
//TCut MvaB0 = "(candMvaDl>0.4320&&candType==3)||(candMvaDl>0.3171&&candType==4)";     //1.5x
// TCut MvaBp = "(candMvaDl>0.458&&candType==1)||(candMvaDl>0.344&&candType==2)";       //2.0x
// TCut MvaB0 = "(candMvaDl>0.360&&candType==3)||(candMvaDl>0.265&&candType==4)";       //2.0x
//TCut MvaBp = "(candMvaDl>0.381&&candType==1)||(candMvaDl>0.2804&&candType==2)";      //2.5x
//TCut MvaB0 = "(candMvaDl>0.308&&candType==3)||(candMvaDl>0.2240&&candType==4)";      //2.5x
//    TCut MvaBp = "(candMvaDl>0.3065&&candType==1)||(candMvaDl>0.235&&candType==2)";      //3.0x
//    TCut MvaB0 = "(candMvaDl>0.253&&candType==3)||(candMvaDl>0.194&&candType==4)";       //3.0x
//TCut MvaBp = "(candMvaDl>0.3&&candType==1)||(candMvaDl>0.3&&candType==2)";
//TCut MvaB0 = "(candMvaDl>0.3&&candType==3)||(candMvaDl>0.3&&candType==4)";
//TCut MvaBp = "(candMvaDl>0.43&&candType==1)||(candMvaDl>0.28&&candType==2)";
//TCut MvaB0 = "(candMvaDl>0.31&&candType==3)||(candMvaDl>0.21&&candType==4)";
//TCut MvaBp = "(candMvaDl>0.48&&candType==1)||(candMvaDl>0.33&&candType==2)";
//TCut MvaB0 = "(candMvaDl>0.36&&candType==3)||(candMvaDl>0.26&&candType==4)";
//TCut MvaBp = "(candMvaDl>0.53&&candType==1)||(candMvaDl>0.38&&candType==2)";
//TCut MvaB0 = "(candMvaDl>0.41&&candType==3)||(candMvaDl>0.31&&candType==4)";
//TCut MvaBp = "(candMvaDl>0.51&&candType==1)||(candMvaDl>0.28&&candType==2)";
//TCut MvaB0 = "(candMvaDl>0.36&&candType==3)||(candMvaDl>0.21&&candType==4)";
TCut Mva = MvaBp||MvaB0;
TCut MvaBp51 = "(candMvaDl>0.48&&candType==1)||(candMvaDl>0.59&&candType==2)";
TCut MvaB051 = "(candMvaDl>0.38&&candType==3)||(candMvaDl>0.40&&candType==4)";
TCut Mva51 = MvaBp51||MvaB051;
TCut MvaBp50 = "(candMvaDl>0.44&&candType==1)||(candMvaDl>0.48&&candMvaComb>0.40&&candType==2)";
TCut MvaB050 = "(candMvaDl>0.38&&candType==3)||(candMvaDl>0.37&&candMvaComb>0.38&&candType==4)";
TCut Mva50 = MvaBp50||MvaB050;
TCut Mva22 = "(candMvaDl>0.15&&candType==1)||(candMvaDl>0.2 &&candType==2)||(candMvaDl>0.42&&candType==3)||(candMvaDl>0.17&&candType==4)";
TCut cosT = "abs(candCosT)<.8";
TCut basic = PMiss+Q2;
TCut MvaAll = basic+Mva+M2P;
TCut MvaAll51 = basic+Mva51+M2P;

TCut CSq2 = PMiss+Mva+M2P+"candQ2<5";
TCut CSpl = PMiss+Mva+M2P+Q2+"candPstarLep>1.5";

TCut mpi0 = "mpi0>.125&&mpi0<.145";
TCut Mpi0 = "mpi0>.12&&mpi0<.15";
TCut dssee = "eextrapi0<.5&&ppi0>.4";

TCut dssacc = "mm2pi0>-4&&mm2pi0<12&&candPstarLep>0&&candPstarLep<2.4&&candMES>5.2&&candMES<5.3&&pmisspi0>.2";
TCut bestpi0 = "bestepi0==1";
TCut dssMvaDl = "candMvaDssDl>-0.45";
TCut dssMvaComb="candMvaDssComb>-0.35&&candType<3||candMvaDssComb>-0.3&&candType==3||candMvaDssComb>-0.35&&candType==4";
//TCut dssMvaDl = "candMvaDssDl>-0.5";
//TCut dssMvaComb="candMvaDssComb>-0.4&&candType<3||candMvaDssComb>-0.35&&candType==3||candMvaDssComb>-0.4&&candType==4";
//TCut dssMvaDl = "candMvaDssDl>-0.55";
//TCut dssMvaComb="candMvaDssComb>-0.45&&candType<3||candMvaDssComb>-0.4&&candType==3||candMvaDssComb>-0.45&&candType==4";
TCut dssMva = dssMvaDl+dssMvaComb;
TCut dssMvaDl51 = "candMvaDssDl>-0.15&&candType==1||candMvaDssDl>-0.16&&candType==2||candMvaDssDl>-0.05&&candType==3||candMvaDssDl>-0.15&&candType==4";
TCut dssMvaComb51="candMvaDssComb>-0.10&&candType==1||candMvaDssComb>-0.12&&candType==2||candMvaDssComb>-0.15&&candType==3||candMvaDssComb>-0.10&&candType==4";
TCut dssMva51 = dssMvaDl51+dssMvaComb51;
TCut dssMvaBp22 = "(candMvaDssComb>0.13&&candMvaDssDl>0.15&&candType==1)||(candMvaDssComb>0.13&&candMvaDssDl>0.15&&candType==2)";
TCut dssMvaB022 = "(candMvaDssComb>0.15&&candMvaDssDl>0.11&&candType==3)||(candMvaDssComb>0.10&&candMvaDssDl>0.17&&candType==4)";
TCut dssMva22 = dssMvaBp22||dssMvaB022;
TCut dss = dssacc+bestpi0;
TCut dssMvaAll = dss+cosT+Q2+Mpi0+dssMva+!MvaAll;
TCut dssMvaAll51 = dss+cosT+dssMva51+!MvaAll51;
TCut dsseeAll = dss+mpi0+dssee;

TCut dsssig = "MCType>=13&&MCType<=14";
TCut dsscomb = "MCType==0";
TCut dssDl = "MCType>0&&MCType<13";

#endif	// CUTS_HH

