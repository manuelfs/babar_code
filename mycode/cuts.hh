#include "TCut.h"

TCut M2P = "candPstarLep>0&&candPstarLep<2.4&&candM2>-4&&candM2<12";
TCut PMiss = "candPMiss>.2";
TCut Q2 = "candQ2>4";
TCut ee1 = "candType<3&&candEExtra<.2";
TCut ee2 = "candType==3&&candEExtra<.15";
TCut ee3 = "candType==4&&candEExtra<.3";
TCut ee = (ee1||ee2||ee3);
TCut MvaBp = "(candMvaDl>0.48&&candType==1)||(candMvaDl>0.59&&candType==2)";
TCut MvaB0 = "(candMvaDl>0.38&&candType==3)||(candMvaDl>0.40&&candType==4)";
TCut Mva = MvaBp||MvaB0;
TCut MvaBp50 = "(candMvaDl>0.44&&candType==1)||(candMvaDl>0.48&&candMvaComb>0.40&&candType==2)";
TCut MvaB050 = "(candMvaDl>0.38&&candType==3)||(candMvaDl>0.37&&candMvaComb>0.38&&candType==4)";
TCut Mva50 = MvaBp||MvaB0;
TCut Mva22 = "(candMvaDl>0.15&&candType==1)||(candMvaDl>0.2 &&candType==2)||(candMvaDl>0.42&&candType==3)||(candMvaDl>0.17&&candType==4)";
TCut cosT = "abs(candCosT)<.8";
TCut basic = PMiss+Q2;
TCut MvaAll = basic+Mva+M2P;
TCut MvaAll22 = basic+Mva22+M2P;

TCut CSq2 = PMiss+Mva+M2P+"candQ2<5";
TCut CSpl = PMiss+Mva+M2P+Q2+"candPstarLep>1.5";

TCut mpi0 = "mpi0>.125&&mpi0<.145";
TCut dssee = "eextrapi0<.5&&ppi0>.4";

TCut dssacc = "mm2pi0>-4&&mm2pi0<12&&candPstarLep>0&&candPstarLep<2.4&&candMES>5.2&&candMES<5.3&&pmisspi0>.2";
TCut bestpi0 = "bestepi0==1";
TCut dssMvaDl = "candMvaDssDl>-0.15&&candType==1||candMvaDssDl>-0.16&&candType==2||candMvaDssDl>-0.05&&candType==3||candMvaDssDl>-0.15&&candType==4";
TCut dssMvaComb="candMvaDssComb>-0.10&&candType==1||candMvaDssComb>-0.12&&candType==2||candMvaDssComb>-0.15&&candType==3||candMvaDssComb>-0.10&&candType==4";
TCut dssMva = dssMvaDl+dssMvaComb;
TCut dssMvaBp22 = "(candMvaDssComb>0.13&&candMvaDssDl>0.15&&candType==1)||(candMvaDssComb>0.13&&candMvaDssDl>0.15&&candType==2)";
TCut dssMvaB022 = "(candMvaDssComb>0.15&&candMvaDssDl>0.11&&candType==3)||(candMvaDssComb>0.10&&candMvaDssDl>0.17&&candType==4)";
TCut dssMva22 = dssMvaBp22||dssMvaB022;
TCut dss = dssacc+bestpi0;
TCut dssMvaAll = dss+cosT+Q2+dssMva+!MvaAll;
TCut dssMvaAll22 = dss+cosT+dssMva22+!MvaAll22;
TCut dsseeAll = dss+mpi0+dssee;

TCut dsssig = "MCType>12";
TCut dsscomb = "MCType==0";
TCut dssDl = "MCType>0&&MCType<13";


