//everything ends in a cut called "basic", or "basic2" without the mES cut

TCut acc = "candM2>-4&&candM2<12&&candPstarLep>0&&candPstarLep<2.4&&candMES>5.2&&candMES<5.3&&candPMiss>.2";
TCut kin1 = "(candType==1||candType==3)&&candPstarLep<2.31";
TCut kin2 = "(candType==2||candType==4)&&candPstarLep<2.26";
TCut kin = (kin1||kin2);
TCut fit = "candFitProb>-.001";
TCut mes = "candMES>5.27";
TCut bbb1 = "candExtraTracks==0&&candThetaLep>.4&&candThetaLep<2.6&&candQ2>4";
TCut ee1 = "candType<3&&candEExtra<.2";
TCut ee2 = "candType==3&&candEExtra<.15";
TCut ee3 = "candType==4&&candEExtra<.3";
// TCut ee1 = "candType<3&&candEExtra<.35";
// TCut ee2 = "candType==3&&candEExtra<.25";
// TCut ee3 = "candType==4&&candEExtra<.45";
TCut ee = (ee1||ee2||ee3);
TCut basic = acc+bbb1+fit+mes+ee+kin;
TCut basic2 = acc+bbb1+fit+ee+kin;
TCut basicnoee = acc+bbb1+fit+mes+kin;
TCut kscut = "extratracksks<3&&sqrt(vxks*vxks+vyks*vyks)<40&&abs(vzks)<40&&pks>.3&&mks>.491&&mks<.506";
TCut bbbks = "candThetaLep>.4&&candThetaLep<2.6&&candQ2>4";

TCut mpi0 = "mpi0>.125&&mpi0<.145";
TCut dssee1 = "candType==1&&bestepi0&&eextrapi0<.5&&ppi0>.4";
TCut dssee2 = "candType==2&&bestepi0&&eextrapi0<.5&&ppi0>.4";
TCut dssee3 = "candType==3&&bestepi0&&eextrapi0<.5&&ppi0>.4";
TCut dssee4 = "candType==4&&bestepi0&&eextrapi0<.5&&ppi0>.4";

TCut dssacc = "mm2pi0>-4&&mm2pi0<12&&candPstarLep>0&&candPstarLep<2.4&&candMES>5.2&&candMES<5.3&&pmisspi0>.2";
TCut dssee = (dssee1||dssee2||dssee3||dssee4);
TCut dss = dssacc+bbb1+fit+mes+mpi0+dssee+kin;
TCut dss2 = dssacc+bbb1+fit+mpi0+dssee+kin;

TCut dpi0 = "bestdpi0&&m1dpi0>0.125&&m1dpi0<0.145&&m2dpi0>0.125&&m2dpi0<0.145";
